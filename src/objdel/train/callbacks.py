"""Level 2 controller for Phases 2-3 (plan D9, §9): best-checkpoint tracking, the hack audit on held-out evaluation
with rollback, and guard switches that change TRL's reward weights and KL coefficient while training runs.

Held-out evaluation is TRL's own: GRPOConfig(eval_strategy="steps", eval_steps=N, eval_on_start=True) with held-out
prompts (objdel.train.dataset.grpo_rows, validation split). TRL samples G = 8 completions per prompt at the training
temperature and calls the same TRLReward functions, so on_evaluate receives the typed held-out curves
(eval_rewards/<term>/mean, eval_error/<code>, eval_defect/<class>); eval_on_start evaluates the SFT weights at step 0.

Controller(trainer, metric, alarm, baseline, audit):
- hack audit (plan §9, break-it B4-4 / issue I21): TRLReward records the completions it scores; the records are
  cleared at every step end, so at on_evaluate they are exactly the held-out ones (transformers runs on_step_end,
  then logging and evaluation). They are gathered over processes (all take the same decision) and audited by
  objdel.eval.alarms_for against `baseline` (an audit of the SFT checkpoint; default: the step-0 evaluation's).
  Statistics are eval_audit/<structure or dataset>/<statistic>; eval_audit/fired lists the checks that fire now and
  not on the baseline (rolling back cannot fix what the SFT checkpoint already fails). Trainer.evaluate logs its
  metrics before on_evaluate, so the controller logs the numeric statistics itself as one more log entry of that
  step (log_history and report_to, through on_log), and the audit event in controller.jsonl holds them all
  (REVIEW.md #3).
- alarm(metrics) -> reason | None, default: eval_audit/fired. Checked BEFORE `best` is updated (break-it B2-1): an
  alarmed evaluation never becomes the rollback target. Hedges score a higher IoU_tc than honest output (M8 demo), so
  "metric improved and alarm fired" is the expected case.
- best: the trainable parameters go to <output_dir>/best/trainable-<rank>.pt when `metric` (default
  eval_rewards/iou_tc/mean, legality-gated, objdel.train.rewards) improves on an evaluation without alarm; the
  initial weights are saved at construction, so an early rollback restores the SFT weights (break-it B2-2). A
  rollback restores the parameters and clears the optimizer state (AdamW re-initialises lazily); the LR is constant
  (D8).
- switch(reward_weights={name: w}, beta=b): applied at the next optimizer step's start, before its generation
  (GRPOTrainer reads self.reward_weights and self.beta there; tests/test_train_trl.py). A weight w on a weight-0
  logger (e.g. g_area) adds w (term - group mean) to every advantage: the guard.
- by hand: <output_dir>/control.json, e.g. {"reward_weights": {"g_area": 0.5}, "beta": 0.01, "rollback": "reason"},
  read by the main process when its mtime changes and broadcast to every process. Write it atomically (temporary
  file, then os.replace). An unparsable file (editor mid-save) is re-read at the next step; invalid content is
  logged and ignored (break-it B2-3). Neither stops training.
- events: every best / switch / rollback / audit / control_error, also appended to <output_dir>/controller.jsonl
  on the main process.

Refused: num_generations or num_generations_eval != REWARD_V4.group (break-it B2-5); raising beta from 0 without a
reference model (TRL builds none when beta = 0 at init without PEFT, lab t_stack T3).

References:
- transformers TrainerCallback / TrainerControl: https://huggingface.co/docs/transformers/main_classes/callback
- transformers 5.20 Trainer loop (on_step_end, then _maybe_log_save_evaluate; evaluate logs, then on_evaluate;
  Trainer.log appends to log_history and calls on_log):
  https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/trainer.py
- TRL (version pinned in pyproject.toml) GRPOTrainer (reward_weights, beta, ref_model, num_generations_eval;
  log writes the completions table):
  https://github.com/huggingface/trl/blob/v1.17.1/trl/trainer/grpo_trainer.py
- accelerate gather_object / broadcast_object_list (no-ops in one process):
  https://huggingface.co/docs/accelerate/main/en/package_reference/utilities
- os.replace is atomic on POSIX: https://docs.python.org/3/library/os.html#os.replace
- J. Yao et al. 2026, Multimodal Reward Hacking in Reinforcement Learning (failures measured against the SFT model,
  monitored during training): https://arxiv.org/abs/2607.09492
- Lab: experiments/phase1/lab/t_stack/t3_grpo.py (controller callback), lab/t_dynamics/REPORT.md (typed alarm +
  rollback; the alarm lags, hence rollback to the best held-out checkpoint).
"""
import json
from pathlib import Path

import torch
from accelerate.utils import broadcast_object_list, gather_object
from peft import PeftModel
from transformers import Trainer, TrainerCallback

from objdel.eval import alarms_for
from objdel.reward import REWARD_V4

FIRED = "eval_audit/fired"


def audit_alarm(metrics: dict) -> str | None:
    """The default alarm: the audit checks that fire now and did not fire on the baseline."""
    return ", ".join(metrics.get(FIRED, [])) or None


def _flat(audit: dict) -> dict:
    """eval_audit/<structure or dataset>/<statistic> for every scalar statistic of the audit."""
    rows = {**audit["structures"], **audit["datasets"]}
    return {f"eval_audit/{s}/{k}": v for s, st in rows.items() if st for k, v in st.items()
            if v is None or isinstance(v, int | float)}


class Controller(TrainerCallback):
    def __init__(self, trainer, metric: str = "eval_rewards/iou_tc/mean", alarm=audit_alarm,
                 baseline: dict | None = None, audit: bool = True):
        self.trainer, self.metric, self.alarm, self.baseline = trainer, metric, alarm, baseline
        self.out = Path(trainer.args.output_dir)
        self.events, self._pending = [], {}
        self._control_mtime = self._control_bad = None
        self.reward = next((f.trl_reward for f in getattr(trainer, "reward_funcs", ()) if hasattr(f, "trl_reward")),
                           None) if audit else None
        if self.reward is not None:
            sizes = (trainer.num_generations, trainer.num_generations_eval)
            if any(n != REWARD_V4.group for n in sizes):
                raise ValueError(f"num_generations / num_generations_eval = {sizes}: reward v4 scores groups of "
                                 f"exactly {REWARD_V4.group}")
            self.reward.recorded = []
        self.best = {"step": trainer.state.global_step, "value": float("-inf"), "path": self._save()}  # initial

    def switch(self, reward_weights: dict[str, float] | None = None, beta: float | None = None) -> None:
        t = self.trainer
        reward_weights = None if reward_weights is None else {k: float(v) for k, v in dict(reward_weights).items()}
        beta = None if beta is None else float(beta)
        unknown = set(reward_weights or ()) - set(t.reward_func_names)
        if unknown:
            raise ValueError(f"unknown reward functions {sorted(unknown)}; known: {t.reward_func_names}")
        if beta and t.beta == 0.0 and t.ref_model is None and not isinstance(t.model, PeftModel):
            raise ValueError("beta was 0 at init without PEFT, so there is no reference model to raise it against")
        self._pending |= {k: v for k, v in (("reward_weights", reward_weights), ("beta", beta)) if v is not None}

    def rollback(self, reason: str) -> None:
        """Restore the best (or, before any, the initial) trainable parameters and clear the optimizer state."""
        saved = torch.load(self.best["path"], map_location="cpu")
        with torch.no_grad():
            for n, p in self.trainer.model.named_parameters():
                if n in saved:
                    p.copy_(saved[n])
        self.trainer.optimizer.state.clear()
        self._event("rollback", to_step=self.best["step"], reason=reason)

    def on_step_begin(self, args, state, control, **kw):
        c = broadcast_object_list([self._read_control() if args.process_index == 0 else None])[0]
        if c is not None:
            try:
                self.switch(c.get("reward_weights"), c.get("beta"))
            except (TypeError, ValueError) as e:
                self._event("control_error", error=str(e))
            else:
                if c.get("rollback"):
                    self.rollback(str(c["rollback"]))
        if self._pending:
            t, p = self.trainer, self._pending
            if "reward_weights" in p:
                w = t.reward_weights.clone()
                for name, v in p["reward_weights"].items():
                    w[t.reward_func_names.index(name)] = v
                t.reward_weights = w
            if "beta" in p:
                t.beta = p["beta"]
            self._event("switch", **p)
            self._pending = {}

    def on_step_end(self, args, state, control, **kw):
        if self.reward is not None:
            self.reward.recorded.clear()  # training completions; the held-out ones arrive after this

    def on_evaluate(self, args, state, control, metrics=None, **kw):
        metrics = {} if metrics is None else metrics
        if self.reward is not None:
            records = gather_object(list(self.reward.recorded))  # a copy: in one process gather returns its input
            self.reward.recorded.clear()
            if records:
                metrics.update(self._audit(records, state.global_step))
        reason = self.alarm(metrics) if self.alarm else None
        if reason:  # before `best` (B2-1): an alarmed evaluation never becomes the rollback target
            self.rollback(reason)
            return
        value = metrics.get(self.metric)
        if value is not None and value > self.best["value"]:
            self.best = {"step": state.global_step, "value": float(value), "path": self._save()}
            self._event("best", value=float(value))

    def _audit(self, records, step: int) -> dict:
        audit = alarms_for(records, self.reward.caches, self.baseline)
        if self.baseline is None and step == 0:  # eval_on_start: the initial (SFT) weights
            self.baseline = audit
        new = sorted(set(audit["fired"]) - set(self.baseline["fired"] if self.baseline else ()))
        stats = _flat(audit)
        self._event("audit", fired=audit["fired"], new=new, completions=len(records), stats=stats)
        # Trainer.evaluate logged its metrics before on_evaluate, so the statistics get one more log entry here.
        # Trainer.log, not GRPOTrainer.log, which would rewrite and reprint the completions table.
        Trainer.log(self.trainer, {k: v for k, v in stats.items() if v is not None})
        return stats | {FIRED: new}

    def _read_control(self) -> dict | None:
        """control.json if it changed since the last read and parses as a JSON object, else None."""
        path = self.out / "control.json"
        try:
            mtime = path.stat().st_mtime_ns
        except OSError:
            return None
        if mtime == self._control_mtime:
            return None
        try:
            c = json.loads(path.read_text())
        except (OSError, ValueError) as e:  # half-written (editor mid-save): retry at the next step, log once
            if mtime != self._control_bad:
                self._control_bad = mtime
                self._event("control_error", error=f"control.json unreadable, retrying: {e}")
            return None
        self._control_mtime = mtime
        if not isinstance(c, dict):
            self._event("control_error", error="control.json is not a JSON object")
            return None
        return c

    def _save(self) -> str:
        path = self.out / "best" / f"trainable-{self.trainer.args.process_index}.pt"
        path.parent.mkdir(parents=True, exist_ok=True)
        torch.save({n: p.detach().cpu() for n, p in self.trainer.model.named_parameters() if p.requires_grad}, path)
        return str(path)

    def _event(self, kind: str, **detail) -> None:
        e = {"step": self.trainer.state.global_step, "kind": kind, **detail}
        self.events.append(e)
        if self.trainer.args.process_index == 0:
            self.out.mkdir(parents=True, exist_ok=True)
            with open(self.out / "controller.jsonl", "a") as f:
                f.write(json.dumps(e) + "\n")


class TraceStep(TrainerCallback):
    """torch.profiler on one optimizer step after `wait` steps and one profiler warm-up step (CPU ops, and CUDA
    kernels on a GPU): <out>/trace.json.gz and ops.txt (top ops by device and CPU time). train_sft `profile:`,
    g_probe --trace."""

    def __init__(self, out: str, wait: int):
        self.out, self.wait = Path(out), wait

    def on_train_begin(self, args, state, control, **kw):
        from torch.profiler import ProfilerActivity, profile, schedule

        acts = [ProfilerActivity.CPU] + [ProfilerActivity.CUDA] * torch.cuda.is_available()
        self.prof = profile(activities=acts, schedule=schedule(wait=self.wait, warmup=1, active=1),
                            on_trace_ready=self.write)  # fmt: skip
        self.prof.start()

    def on_step_end(self, args, state, control, **kw):
        self.prof.step()

    def on_train_end(self, args, state, control, **kw):
        self.prof.stop()

    def write(self, prof):
        self.out.mkdir(parents=True, exist_ok=True)
        prof.export_chrome_trace(str(self.out / "trace.json.gz"))
        ops = prof.key_averages()
        (self.out / "ops.txt").write_text("\n\n".join(ops.table(sort_by=k, row_limit=40) for k in
                                                       ("self_device_time_total", "self_cpu_time_total")))
