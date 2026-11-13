"""The body of a full GRPO run (train_rl): `run_rl(stage)` trains an objdel.core.spec.RlStage for --steps, saves every
--save-steps and keeps the tuned model; `merge_last(stage)` finishes a run that stopped early.

Every refusal comes before any model or engine load, the session budget first (session_budget: the flags, else the
environment run.sh exports; BudgetFlagError when a budget has no session start or rate).

The trainer is the G-probe's (objdel.rl.build.build_trainer): FastGRPO with OpenOnce and ImageMemo, vLLM colocate
drawing, scoring on the GPU next to the CUDA checker, G = 8, 32 LIDC train prompts per step, the probe LR unless --lr.

Rows. A rows root (build_rl_rows: manifest.json has `schedule`) is trained in its file order, which is the schedule.
Before any model or engine loads, train_rl refuses (ValueError) a manifest whose SHA-256 is not
--rows-manifest-sha256, rows that are not steps x per_step keys with the schedule's ids_sha256, a sampler that does
not draw per_step prompts per optimizer step in order (TRL's RepeatSampler: generation_batch_size // num_generations,
steps_per_generation == gradient_accumulation_steps, num_iterations 1), and --steps beyond the schedule or not whole
blocks. shuffle_dataset is False, so an HF resume skips the consumed batches and keeps the order. A root without
`schedule` (data/prepared) trains shuffled with disagreement items x4.

Saves. Every --save-steps it writes <out>/trainer/checkpoint-N (LoRA adapter, optimizer, scheduler, RNG, trainer
state), so --resume continues from the last one after a lost box (a run without --resume refuses an --out that holds
checkpoints; --resume refuses one past --steps). At the end the unmerged adapter goes to <out>/final_adapter, then
the LoRA is merged into the base and saved with its processor to <out>/final, a plain model folder (select_ckpt,
final_eval and vLLM read it unchanged). --merge-last does the same from the last checkpoint of a run that did not
reach its end (rl.json: partial true), and nothing when <out>/final exists (save_final writes final.tmp and renames
it, so final is complete; a missing rl.json is written with steps null).

Budget. --budget-usd stops and saves (callbacks.BudgetStop) when the session spend --fixed-usd + hours since
--session-start x --rate-usd-per-h reaches it; a flag left out is read from OBJDEL_SESSION_START, OBJDEL_RATE_USD_PER_H
and OBJDEL_FIXED_USD (run.sh exports them from the active run) by session_budget, not by the stage.

Held-out curve (off by default): --heldout VIEW evaluates --heldout-items val-shape items of that eval view every
--heldout-every steps, 8 answers each; rl.json "eval" holds eval/iou_tc and eval/allmiss_share.

<out>/rl.json: steps, LR, numerics, opsm_start_delta, beta, reward (name, fingerprint, c_topo, loc_width, area_gate),
rows provenance, timing, `report` (run_report: the clipped_ratio slope, the truncated_grad_share median and max) and
TRL's log. --reward v5 trains on reward v5 and adds callbacks.LocGuard (stop and save; loc_guard_stopped).
"""

import functools
import gc
import json
import math
import numbers
import os
import shutil
import sys
import time
from pathlib import Path

from objdel.core.hashing import bytes_sha256
from objdel.core.io import plain, read_json, write_json
from objdel.policy.errors import Refusal, exits_on_refusal
from objdel.rl import build

usage = functools.partial(Refusal.usage, "train_rl")


@exits_on_refusal
def checked(a):
    """The RlStage `a` with f_check resolved (--f-check, --no-f-check, else on under the f_check.PRESETS), after the
    rules that need objdel.rl (numerics preset names, the --grpo keys train_rl owns); exit 2 on a broken one."""
    from objdel.rl import config as configs
    from objdel.rl import f_check

    numerics = a.numerics or configs.DEFAULT_NUMERICS
    if numerics not in configs.NUMERICS:
        raise usage(f"--numerics {a.numerics}: one of {sorted(configs.NUMERICS)}")
    if a.f_check and numerics not in f_check.PRESETS:
        raise usage(f"--f-check: only under --numerics {' or '.join(f_check.PRESETS)}, not {numerics}")
    a = a.model_copy(update={"f_check": a.f_check or (numerics in f_check.PRESETS and not a.no_f_check)})
    if a.f_check and numerics == "lora-seqmask" and not a.merge_last and a.b1_adapter is None:
        raise usage("--f-check under --numerics lora-seqmask needs --b1-adapter (the comparator's B1 adapter)")
    if a.opsm_start_delta is not None and (numerics != "bf16-tis-opsm" or a.opsm_start_delta <= 0):
        raise usage("--opsm-start-delta is a float > 0, valid only with --numerics bf16-tis-opsm")
    if "off_policy_mask_threshold" in a.grpo:
        raise usage("--grpo cannot set off_policy_mask_threshold: use --opsm-start-delta")
    if a.grpo.get("logging_steps", 1) != 1:  # seq_gaps and cut_seq_gaps need one block per step
        raise usage("--grpo logging_steps: must stay 1 (seq_gaps.f32 holds one block per optimizer step)")
    if (why := stage_refusal(a.model, a.dataset, a.no_ckpt_layers)) is not None:
        raise usage(why)
    return a


def stage_refusal(model: str, dataset: str, no_ckpt_layers: int | None) -> str | None:
    """Why an RL stage (train_rl, the G-probe) cannot run: --no-ckpt-layers above the model's text decoder layers,
    then rl_refusal(dataset); else None."""
    from objdel.policy import models as policy_models

    if no_ckpt_layers and no_ckpt_layers > (n := policy_models.get(model).text_layers):
        return f"--no-ckpt-layers {no_ckpt_layers}: {model} has {n} text decoder layers"
    return rl_refusal(dataset)


def rl_refusal(dataset: str) -> str | None:
    """Why RL cannot train on the dataset (unknown, or fewer raters per item than the reward needs), else None."""
    from objdel.core import datasets
    from objdel.verify.reward.cache import MIN_RATERS

    if dataset not in datasets.DATASETS:
        return f"--dataset {dataset}: not one of {sorted(datasets.DATASETS)} (objdel.core.datasets)"
    if (k := datasets.get(dataset).min_raters) < MIN_RATERS:
        return (f"--dataset {dataset}: {k} rater(s) per item; the RL reward needs {MIN_RATERS} (its band and "
                "topology vote), and there is no single-rater mode")
    return None


def check_flags(a) -> None:
    """Refusals that need no model: the completion limit (set from the outline count, never above the model's
    budget) and the runaway cap (C11: a capped answer is cut mid-text, so masking it would bias the update toward
    more capping)."""
    from objdel.policy import models as policy_models

    over = a.grpo
    budget = policy_models.get(a.model).max_completion
    if over.get("max_completion_length", budget) > budget:
        raise ValueError(f"max_completion_length {over['max_completion_length']} is above the model's {budget}")
    if a.runaway_cap and over.get("mask_truncated_completions", True):
        raise ValueError("--runaway-cap needs --grpo '{\"mask_truncated_completions\": false}'")


def numerics_overrides(a) -> dict:
    """The grpo_config arguments of --numerics and --opsm-start-delta (absent: the preset's own threshold)."""
    from objdel.rl import config as rl_config

    return {"numerics": a.numerics or rl_config.DEFAULT_NUMERICS} | ({} if a.opsm_start_delta is None
                                       else {"off_policy_mask_threshold": a.opsm_start_delta})


class BudgetFlagError(ValueError):
    """--budget-usd without a session start or a rate (neither a flag nor its environment variable)."""


def session_budget(a, environ=os.environ):
    """callbacks.BudgetStop of the budget flags, or None without --budget-usd. A session flag left out is read from
    its environment variable (run.sh exports them; empty counts as unset); the flag wins. Each value is finite and
    >= 0, the rate > 0."""
    if a.budget_usd is None:
        return None
    got = {}
    for field, var in (("session_start", "OBJDEL_SESSION_START"), ("rate_usd_per_h", "OBJDEL_RATE_USD_PER_H"),
                       ("fixed_usd", "OBJDEL_FIXED_USD")):
        v = getattr(a, field)
        if v is None and environ.get(var):
            try:
                v = float(environ[var])
            except ValueError:
                v = math.nan
            if not math.isfinite(v) or v < 0 or (v == 0 and field == "rate_usd_per_h"):
                raise BudgetFlagError(f"{var}={environ[var]!r}: needs a finite number >= 0 (the rate > 0)")
        got[field] = v
    for field, var in (("session_start", "OBJDEL_SESSION_START"), ("rate_usd_per_h", "OBJDEL_RATE_USD_PER_H")):
        if got[field] is None:
            raise BudgetFlagError(f"--budget-usd needs --{field.replace('_', '-')} or {var} (run.sh exports it)")
    from objdel.policy.callbacks import BudgetStop

    return BudgetStop(got["rate_usd_per_h"], a.budget_usd, got["session_start"], fixed_usd=got["fixed_usd"] or 0.0)


def rows_schedule(root: Path, manifest_sha256: str | None) -> dict | None:
    """The manifest's `schedule` (None for a root without one), after checking the manifest's SHA-256 against the
    pinned value; a schedule root needs it."""
    raw = (Path(root) / "manifest.json").read_bytes()
    schedule = json.loads(raw).get("schedule")
    if schedule is not None and not (isinstance(schedule, dict) and all(
            k in schedule for k in ("steps", "per_step", "block_steps", "ids_sha256"))):
        raise ValueError(f"{root}/manifest.json: `schedule` needs steps, per_step, block_steps and ids_sha256")
    if schedule is not None and manifest_sha256 is None:
        raise ValueError(f"{root}: a rows root with a schedule needs --rows-manifest-sha256")
    if manifest_sha256 is not None and (got := bytes_sha256(raw)) != manifest_sha256:
        raise ValueError(f"{root}/manifest.json SHA-256 {got} is not --rows-manifest-sha256 {manifest_sha256}")
    return schedule


def check_rows(keys, schedule: dict) -> None:
    """The rows are the schedule: steps x per_step keys whose newline-joined SHA-256 is ids_sha256 (a dropped,
    swapped or reordered row would shift every later step)."""
    keys, n = list(keys), schedule["steps"] * schedule["per_step"]
    if len(keys) != n:
        raise ValueError(f"rows: {len(keys)} rows, the schedule holds steps x per_step = {n}")
    if bytes_sha256("\n".join(keys).encode()) != schedule["ids_sha256"]:
        raise ValueError("rows: ids_sha256 differs from the schedule's (a row dropped, swapped or reordered)")


def check_sampler(cfg, schedule: dict, steps: int) -> None:
    """TRL's RepeatSampler on unshuffled rows draws the schedule's per_step prompts per optimizer step in file order
    only with this shape; --steps must be whole blocks within the schedule."""
    if cfg.generation_batch_size // cfg.num_generations != schedule["per_step"]:
        raise ValueError(f"generation_batch_size {cfg.generation_batch_size} // num_generations "
                         f"{cfg.num_generations} is not the schedule's per_step {schedule['per_step']}")
    if cfg.steps_per_generation != cfg.gradient_accumulation_steps:
        raise ValueError(f"steps_per_generation {cfg.steps_per_generation} != gradient_accumulation_steps "
                         f"{cfg.gradient_accumulation_steps}: one generation batch must be one optimizer step")
    if cfg.num_iterations != 1:
        raise ValueError(f"num_iterations {cfg.num_iterations}: the schedule gives each step its own prompts")
    if steps < 1:
        raise ValueError(f"--steps {steps}: needs >= 1")
    if steps > schedule["steps"]:
        raise ValueError(f"--steps {steps} > the schedule's {schedule['steps']} steps")
    if steps % schedule["block_steps"]:
        raise ValueError(f"--steps {steps} is not a whole number of {schedule['block_steps']}-step blocks")


def heldout_rows(view: Path, n: int, dataset: str):
    """n val-shape items of an eval view (a seeded draw in record order), as GRPO rows: val-select items never enter
    training's held-out curve."""

    from objdel.data.splits import val_half
    from objdel.policy import rows as policy_rows

    policy_rows.check_root(view, dataset, "val")
    shape = {r["item_id"] for r in policy_rows.records(view, dataset, "val") if val_half(r["group_id"]) == "shape"}
    rows = [r for r in policy_rows.prompts(view, dataset, "val") if r["key"] in shape]
    if not rows:
        raise ValueError(f"--heldout {view}: no val-shape RL item")
    if n < len(rows):
        rows = policy_rows.draw_subset(rows, n, 0)
    return policy_rows.as_dataset(rows)


def save_final(model, proc, out: Path) -> None:
    """The unmerged LoRA adapter to out/final_adapter, then the merged model and its processor to out/final.tmp,
    renamed to out/final at the end: a folder named final is always complete."""
    from objdel.policy import models as policy_models

    policy_models.save_merged(model, proc, out / "final", adapter_to=out / "final_adapter", tmp=out / "final.tmp")


def last_checkpoint(out: Path) -> tuple[str | None, int]:
    """The newest checkpoint-N under out/trainer and N, (None, 0) when there is none. transformers'
    get_last_checkpoint: it does not check that the checkpoint is complete."""
    from transformers.trainer_utils import get_last_checkpoint

    last = get_last_checkpoint(str(out / "trainer")) if (out / "trainer").is_dir() else None
    return last, int(Path(last).name.split("-")[1]) if last else 0


def cut_seq_gaps(trainer_dir: Path, step: int) -> None:
    """seq_gaps.f32 back to the first `step` steps (rows a crashed run logged after its last checkpoint go; a fresh
    run, step 0, deletes it), so step i of the file stays optimizer step i after --resume."""
    f, meta = trainer_dir / "seq_gaps.f32", trainer_dir / "seq_gaps.json"
    if not step or not meta.is_file() or not f.is_file():
        f.unlink(missing_ok=True)
        meta.unlink(missing_ok=True)
        return
    m = read_json(meta)
    with open(f, "r+b") as fh:
        fh.truncate(min(f.stat().st_size, step * m["rows_per_step"] * 3 * 4))
    m["steps"] = f.stat().st_size // (m["rows_per_step"] * 3 * 4)
    meta.write_text(json.dumps(m) + "\n")


def run_report(log: list[dict]) -> dict:
    """Run-report line: the least-squares slope of completions/clipped_ratio over the steps (None below two
    steps) and the median and max of completions/truncated_grad_share (None when the run logged none)."""
    import numpy as np

    def series(key):
        return [(h["step"], h[key]) for h in log if key in h]

    clipped, share = series("completions/clipped_ratio"), [v for _, v in series("completions/truncated_grad_share")]
    return {"clipped_ratio_slope": float(np.polyfit(*zip(*clipped, strict=True), 1)[0]) if len(clipped) > 1 else None,
            "truncated_grad_share": {"median": float(np.median(share)), "max": float(max(share))} if share else None}


def merge_last(a) -> int:
    """--merge-last: <out>/final from the last checkpoint when the run ended before its own save."""
    out = Path(a.out)
    final = out / "final"
    if (final / "config.json").is_file() and any(final.glob("*.safetensors")):  # save_final renames it when complete
        print(f"{final} exists: nothing to merge")
        if (out / "rl.json").is_file():
            return 0
        return save_json(out / "rl.json", {  # the box died between save_final and save_json: the step is unknown
            "model": a.model, "steps": None, "partial": None, "ckpt": str(a.ckpt), "final": str(final),
            "recovered": "final complete, rl.json missing"}, keys=("steps", "partial", "recovered"))
    shutil.rmtree(out / "final.tmp", ignore_errors=True)  # the box died inside save_final
    last, step = last_checkpoint(out)
    if last is None:
        print(f"--merge-last: no checkpoint in {out / 'trainer'}", file=sys.stderr)
        return 1
    from objdel.policy import models as policy_models

    save_final(policy_models.adapted(a.ckpt, last), policy_models.processor(a.model, a.ckpt), out)
    return save_json(out / "rl.json", {
        "model": a.model, "steps": step, "partial": True, "merged_from": last,
        "ckpt": str(a.ckpt),
        "final": str(final),
    }, keys=("steps", "partial", "merged_from"))


def with_caps(train, a):
    """C11: the rows plus a `cap` column (FastGRPO sends it as each prompt's max_tokens); T_i from the row's own
    record."""
    from objdel.policy import models as policy_models
    from objdel.policy import rows as policy_rows
    from objdel.rl import runaway

    tok = policy_models.processor(a.model, a.ckpt).tokenizer
    caps = runaway.row_caps(train["key"], policy_rows.records(a.train, a.dataset, "train"),
                            lambda text: policy_rows.target_tokens(tok, text),
                            a.grpo.get("max_completion_length", policy_models.get(a.model).max_completion),
                            a.runaway_cap, runaway.MARGIN)
    return train.add_column("cap", caps)


def f_check_callback(a, trainer, train):
    """f_check.FCheck of the run's preset. fp16-seqmask: NaN/inf and the engine dtype. lora-seqmask: the comparator,
    restore base hash and timing measured now, before step 1 (f_check.run_comparator on the first generation
    batch's prompts, rows 0-31 in file order, scored by way (b); its output-level arms only when --ckpt is
    rl-full-3ep, f_check.comparator_skip); then the LoRA loads counted on add_lora and the
    packed-layout equality after sync 1 and after the first sleep/wake cycle. An expansion error writes the record
    and raises here."""
    import tempfile

    import torch

    from objdel.rl import f_check

    vg = trainer.vllm_generation
    if a.numerics == "fp16-seqmask":
        engine = vg.llm.llm_engine
        return f_check.FCheck(Path(a.out), a.numerics,
                              {"engine_dtype_ok": lambda: engine.model_config.dtype == torch.float16})
    model, G = vg.model, trainer.num_generations
    prompts = trainer.comparator_prompts([train[i] for i in range(trainer.args.generation_batch_size // G)], G)
    tmp = Path(tempfile.mkdtemp(prefix="objdel-lora-b1-", dir=Path(trainer.lora_dir).parent))
    fixed = f_check.run_comparator(vg, model, prompts, model.config.image_token_id, Path(a.b1_adapter), tmp,
                                   lambda outs, adapter: trainer.way_b_gap(prompts, outs, adapter),
                                   skip=f_check.comparator_skip(a.ckpt))
    cb = f_check.FCheck(Path(a.out), a.numerics, {}, fixed=fixed)
    if not fixed.get("prompt_expansion", True):
        cb.finish()  # writes f_check.json and raises
    engine = vg.llm.llm_engine
    cb.loads = engine.add_lora = f_check.count_calls(engine.add_lora)  # wrapped after the comparator's own load
    probe = f_check.AfterSyncs(vg, lambda: f_check.lora_equal(vg, model), 2 if vg.enable_sleep_mode else 3)
    cb.probes = {"lora_equal": probe.ok}
    return cb


def finite(x):
    """NaN / inf -> None, recursively. A caller applies it to the values that may be NaN by design (a Spearman CI of
    constant gains, a trainer log), so save_json, which refuses NaN, gets null there."""
    if isinstance(x, dict):
        return {k: finite(v) for k, v in x.items()}
    if isinstance(x, list | tuple):
        return [finite(v) for v in x]
    return None if isinstance(x, numbers.Real) and not math.isfinite(x) else x


def save_json(path: Path, doc: dict, keys: tuple[str, ...]) -> int:
    """Strict JSON (numpy values as plain numbers and lists) to path, through write_atomic; prints doc's `keys` on
    one line. A NaN or inf anywhere raises ValueError before anything is written (U25)."""
    write_json(path, doc, sort_keys=False, strict=True, default=plain)
    print(json.dumps({k: doc[k] for k in keys if k in doc}, ensure_ascii=False, default=plain))
    return 0


class Step1Mismatch(RuntimeError):
    """Step 1 of this run differs from --expect-step1's arm (ExpectStep1); nothing is saved."""


def run_rl(a) -> dict:
    """One full run of the RlStage `a` (the body of train_rl after parsing): resume from the last checkpoint, train,
    merge the LoRA, save <out>/final, write <out>/rl.json; returns the rl.json document."""
    from objdel.policy import models as policy_models
    from objdel.policy import rows as policy_rows
    from objdel.policy.callbacks import ExpectStep1, LocGuard
    from objdel.rl import record
    from objdel.rl import rows as rl_rows
    from objdel.verify import store as serve
    from objdel.verify.reward.spec import REWARDS, reward_fingerprint

    stop = session_budget(a)  # every refusal before any model or engine load, the budget first
    check_flags(a)
    out = Path(a.out)
    step1 = ExpectStep1.from_files(a.expect_step1, a.encode_check, out) if a.expect_step1 else None
    if a.record_batches:
        record.open_dir(a.record_batches)
    schedule = rows_schedule(a.train, a.rows_manifest_sha256)
    policy_rows.check_root(a.train, a.dataset, "train")
    train = rl_rows.grpo_rows(a.train, a.dataset, "train",
                              oversample=1 if schedule is not None else policy_rows.DISAGREEMENT_OVERSAMPLE)
    if schedule is not None:
        check_rows(train["key"], schedule)
    if a.runaway_cap:
        train = with_caps(train, a)
    heldout = heldout_rows(a.heldout, a.heldout_items, a.dataset) if a.heldout else None
    last, start = last_checkpoint(out)
    if not a.resume and last:
        raise ValueError(f"{last} exists: a fresh run needs an --out without checkpoints (or pass --resume)")
    if start > a.steps:
        raise ValueError(f"--resume: {last} is past --steps {a.steps} (an older run in this --out?)")
    if a.resume and last is None:
        print(f"--resume: no checkpoint in {out / 'trainer'}; starting at step 0", file=sys.stderr)
    cut_seq_gaps(out / "trainer", start)
    roots = ([Path(a.heldout)] if a.heldout else []) + [Path(a.train)]
    for root in roots[::-1]:
        serve.preload(root, [a.dataset])
    gc.freeze()  # preloaded store: leave it out of the GC scans (it ran inside the retired worker pool)
    over = {"save_strategy": "steps", "save_steps": a.save_steps} | numerics_overrides(a)
    over |= {"shuffle_dataset": False} if schedule is not None else {}
    over |= {"eval_strategy": "steps", "eval_steps": a.heldout_every, "num_generations_eval": 8} if a.heldout else {}
    check = (lambda cfg: check_sampler(cfg, schedule, a.steps)) if schedule is not None else None
    spec = REWARDS[a.reward]
    guard = LocGuard() if a.reward == "v5" else None
    trainer, proc, lr = build.build_trainer(a.model, Path(a.ckpt), train, out, roots, a.steps, a, eval_dataset=heldout,
                                      check=check, calls=[c for c in (stop, guard, step1) if c], reward=spec,
                                      record_dir=a.record_batches, **over)
    if a.f_check:
        trainer.add_callback(f_check_callback(a, trainer, train))
    t0 = time.time()
    trainer.train(resume_from_checkpoint=last)
    t1 = time.time()
    if step1 and step1.mismatch:
        raise Step1Mismatch(f"step 1 differs from {a.expect_step1}: {step1.mismatch}")
    save_final(trainer.accelerator.unwrap_model(trainer.model), proc, out)
    done = trainer.state.global_step
    log, cfg = trainer.state.log_history, trainer.args
    rows = schedule is not None and {"root": str(a.train), "revision": a.rows_revision,
                                     "manifest_sha256": a.rows_manifest_sha256, "ids_sha256": schedule["ids_sha256"]}
    doc = {
        "model": a.model, "steps": done, "lr": lr, "ckpt": str(a.ckpt), "resumed_from": last,
        "final": str(out / "final"),
        "numerics": numerics_overrides(a)["numerics"],
        "opsm_start_delta": a.opsm_start_delta if a.opsm_start_delta is not None else cfg.off_policy_mask_threshold,
        "beta": cfg.beta, "mask_truncated_completions": cfg.mask_truncated_completions,
        "max_completion_length": cfg.max_completion_length,
        "loss_scale": policy_models.get(a.model).max_completion / cfg.max_completion_length,  # dr_grpo divides by L
        "reward": {"name": spec.name, "fingerprint": reward_fingerprint(spec), "c_topo": spec.c_topo,
                   "loc_width": getattr(spec, "loc_width", None), "area_gate": getattr(spec, "area_gate", None)},
        "rows": rows, "train_rows": len(train), "budget_stopped": bool(stop and stop.tripped),
        "loc_guard_stopped": bool(guard and guard.tripped),
        "eval": finite([{"step": h["step"], "eval/iou_tc": h["eval_rewards/iou_tc/mean"],
                         "eval/allmiss_share": h["eval_allmiss_share"]}
                        for h in log if "eval_rewards/iou_tc/mean" in h]),
        "timing_s": {"train": t1 - t0, "per_step": (t1 - t0) / max(1, done - start)},
        "report": finite(run_report(log)), "log": finite(log),  # a trainer log may hold NaN (a logged 0/0): null
    }
    save_json(out / "rl.json", doc, keys=("steps", "lr", "resumed_from", "final", "budget_stopped", "report"))
    return doc
