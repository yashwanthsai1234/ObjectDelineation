"""Level 2 controller for Phases 2-3 (plan D9, §9): best-checkpoint tracking, the hack audit on held-out evaluation
with rollback, and guard switches that change TRL's reward weights and KL coefficient while training runs.

Held-out evaluation is TRL's own: GRPOConfig(eval_strategy="steps", eval_steps=N, eval_on_start=True) with an
eval_dataset of held-out prompts (objdel.train.dataset.grpo_rows on the validation split). TRL samples G = 8
completions per prompt at the training temperature and calls the same TRLReward functions, so on_evaluate receives
the typed held-out curves (eval_rewards/<term>/mean, eval_error/<code>, eval_defect/<class>; GRPOTrainer.log
updates the metrics dict in place before on_evaluate runs). eval_on_start evaluates the initial (SFT) weights at
step 0.

Controller(trainer, metric, alarm, baseline, audit):
- hack audit (plan §9, break-it B4-4 / issue I21): TRLReward records the completions it scores (key, text,
  truncated, W, H, names); the records are cleared at every step end, so at on_evaluate they are exactly the held-out
  evaluation's completions (transformers runs on_step_end, then logging and evaluation). They are gathered over
  processes (every process then takes the same decision) and audited by objdel.eval.alarms_for (eval.metrics per
  item on the prepared caches, then hack_audit) against `baseline` (an audit of the SFT checkpoint; by default
  the audit of the step-0 evaluation). The audit statistics are eval_audit/<structure or dataset>/<statistic>;
  eval_audit/fired lists the checks that fire now and did not fire on the baseline (a check the SFT checkpoint
  already fails cannot be fixed by rolling back to it). transformers logs the evaluation metrics before it calls
  on_evaluate (Trainer.evaluate: self.log(output.metrics), then callback_handler.on_evaluate), so adding them to
  the metrics dict reaches only the alarm and evaluate()'s return value. The controller therefore records them
  itself: the audit event's `stats` (controller.jsonl) holds every statistic, and the numeric ones are logged as
  one more log entry of the evaluation's step (state.log_history and the report_to backends, through on_log)
  (code review REVIEW.md #3).
- alarm(metrics) -> reason | None; the default fires on eval_audit/fired. The alarm is checked BEFORE `best`
  is updated (break-it B2-1): an evaluation that fires it never becomes the best checkpoint, and the rollback
  restores the last evaluation that did not fire. The hedges the audit exists for score a higher IoU_tc than
  honest output (M8 demo), so "metric improved and alarm fired" is the expected case.
- best: the trainable parameters are saved to <output_dir>/best/trainable-<rank>.pt when `metric` (default
  eval_rewards/iou_tc/mean, the legality-gated IoU_tc logger of objdel.train.rewards) improves
  (greater is better) on an evaluation without alarm; the initial weights are saved there at construction, so a
  rollback before any evaluation (alarm or control file) restores the initial (SFT) weights instead of failing
  (break-it B2-2). A rollback restores the parameters and clears the optimizer state (AdamW re-initialises its
  moments lazily); the learning-rate schedule is constant (D8).
- switch(reward_weights={name: w}, beta=b): applied at the next optimizer step's start, before its generation.
  GRPOTrainer reads self.reward_weights in _generate_and_score_completions and self.beta when computing reference
  log-probabilities and the loss, so the change takes effect at that step (tests/test_train_trl.py). A weight w on
  a weight-0 logger (e.g. g_area) adds w (term - group mean) to every advantage: the guard.
- by hand: <output_dir>/control.json, e.g. {"reward_weights": {"g_area": 0.5}, "beta": 0.01, "rollback": "reason"},
  read by the main process whenever its modification time changes and broadcast to every process (all processes
  apply the same switch). Write it atomically (write a temporary file in the same directory, then os.replace).
  A file that does not parse (an editor mid-save) is skipped and re-read at the next step; a parsed file with
  invalid content is logged and ignored (break-it B2-3). Neither stops training.
- events: every best / switch / rollback / audit / control_error, also appended to <output_dir>/controller.jsonl
  on the main process.

The controller refuses a trainer whose num_generations or num_generations_eval is not REWARD_V4.group (break-it
B2-5). Raising beta from 0 needs a reference model: TRL builds none when beta = 0 at init without PEFT (lab
t_stack T3); switch refuses it.

References:
- transformers TrainerCallback / TrainerControl (on_step_end, on_evaluate with the metrics dict; callbacks are
  read-only apart from TrainerControl): https://huggingface.co/docs/transformers/main_classes/callback
- transformers 5.20 Trainer loop (on_step_end, then _maybe_log_save_evaluate; eval_on_start before the first
  epoch; evaluate logs its metrics, then calls on_evaluate; Trainer.log appends to log_history and calls on_log):
  https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/trainer.py
- TRL 1.16 GRPOTrainer (reward_weights, beta, ref_model, num_generations_eval; log writes the completions table):
  https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_trainer.py
- accelerate gather_object / broadcast_object_list (no-ops in one process):
  https://huggingface.co/docs/accelerate/main/en/package_reference/utilities
- os.replace is atomic on POSIX (a reader sees the old or the new file, never a mix):
  https://docs.python.org/3/library/os.html#os.replace
- J. Yao et al. 2026, Multimodal Reward Hacking in Reinforcement Learning (newly rewarded failures measured against
  the SFT model; monitor them during training): https://arxiv.org/abs/2607.09492
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
        self.events: list[dict] = []
        self._pending: dict = {}
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
        # transformers logged the evaluation metrics before on_evaluate (Trainer.evaluate: self.log(output.metrics),
        # then callback_handler.on_evaluate), so the statistics are logged here as one more log_history / report_to
        # entry. Trainer.log, not GRPOTrainer.log, which would rewrite and reprint the completions table.
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
