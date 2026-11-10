"""Full GRPO run: train a checkpoint for --steps, save every --save-steps, keep the tuned model.

    python -I -m experiments.phase2.train_rl --model qwen3-vl-2b --ckpt runs/r1/branches/e2 --train data/prepared
        --out runs/rl --steps 500 --save-steps 50 [--resume]
    python -I -m experiments.phase2.train_rl --model qwen3-vl-2b --ckpt {rd}/ckpt --train {data}/rows/lidc-train-b1
        --rows-revision <hex40> --rows-manifest-sha256 <hex64> --out {rd}/b1-rows --steps 260 --budget-usd 2.05
    python -I -m experiments.phase2.train_rl --merge-last --model qwen3-vl-2b --ckpt {rd}/ckpt --out {rd}/b1-rows
    python -I -m experiments.phase2.train_rl --spec specs/b1_rows.yaml --stage train [--rd RD] [--data DATA]

The flags are an objdel.core.spec.RlStage (cli_stage: one flag per field, --grpo as JSON, --profile-window S,E);
train_rl adds the rules that need objdel.train (the numerics preset names, the --grpo keys it owns).

The trainer is the G-probe's (g_probe.make_trainer): FastGRPO with OpenOnce and ImageMemo, vLLM colocate drawing,
scoring on the GPU next to the CUDA checker, G = 8, 32 LIDC train prompts per step, the probe LR unless --lr.

Rows. A rows root (build_rl_rows: manifest.json has `schedule`) is trained in its file order, which is the schedule:
before any model or engine loads, train_rl refuses (ValueError) a manifest whose SHA-256 is not
--rows-manifest-sha256, rows that are not steps x per_step keys with the schedule's ids_sha256, a sampler that does
not draw per_step prompts per optimizer step in order (TRL's RepeatSampler: generation_batch_size // num_generations,
steps_per_generation == gradient_accumulation_steps, num_iterations 1), and --steps beyond the schedule or not whole
blocks. shuffle_dataset is False, so an HF resume skips the consumed batches and keeps the order. A root without
`schedule` (data/prepared) trains shuffled with disagreement items x4 until B1.8.

Saves. Every --save-steps it writes <out>/trainer/checkpoint-N: the LoRA adapter, optimizer, scheduler, RNG and
trainer state, so --resume continues from the last one after a lost box. At the end the unmerged adapter goes to
<out>/final_adapter, then the LoRA is merged into the base and saved with its processor to <out>/final, a plain
model folder (select_ckpt, final_eval and vLLM read it unchanged). --merge-last does the same from the last
checkpoint of a run that did not reach its end (rl.json: partial true), and nothing when <out>/final exists.

Budget. --budget-usd stops and saves (callbacks.BudgetStop) when the session spend --fixed-usd + hours since
--session-start x --rate-usd-per-h reaches it; a flag left out is read from OBJDEL_SESSION_START, OBJDEL_RATE_USD_PER_H
and OBJDEL_FIXED_USD (run.sh exports them from the active run) at run time, by session_budget, not by the stage.

Held-out curve (off by default): --heldout VIEW evaluates --heldout-items val-shape items (splits.val_half) of that
eval view every --heldout-every steps, 8 answers each; rl.json "eval" holds eval/iou_tc and eval/allmiss_share.

<out>/rl.json: steps, LR, numerics, opsm_start_delta (None outside bf16-tis-opsm), beta, reward (name, fingerprint,
c_topo, loc_width, area_gate), rows provenance, timing and TRL's log. --reward v5 trains on reward v5 and adds
callbacks.LocGuard (stop and save; loc_guard_stopped).
"""
import gc
import hashlib
import json
import math
import os
import shutil
import sys
import time
from pathlib import Path

# lora-seqmask needs vLLM LoRA serving (option F plumbing, B2.2), not built yet
NEEDS_FASTGRPO = ("lora-seqmask",)


def parse(argv=None):
    """The RlStage of the flags (or of --spec/--stage), and the rules that need objdel.train; exit 2 on a broken one."""
    from objdel.core.spec import RlStage, cli_stage
    from objdel.train import configs

    a = cli_stage(RlStage, sys.argv[1:] if argv is None else argv)

    def error(msg: str):
        print(f"train_rl: error: {msg}", file=sys.stderr)
        sys.exit(2)

    if a.numerics is not None and a.numerics not in configs.NUMERICS:
        error(f"--numerics {a.numerics}: one of {sorted(configs.NUMERICS)}")
    if a.numerics in NEEDS_FASTGRPO:
        error(f"--numerics {a.numerics} needs FastGRPO plumbing (vLLM LoRA serving, B2.2), not built yet")
    if a.opsm_start_delta is not None and (a.numerics != "bf16-tis-opsm" or a.opsm_start_delta <= 0):
        error("--opsm-start-delta is a float > 0, valid only with --numerics bf16-tis-opsm")
    if "off_policy_mask_threshold" in a.grpo:
        error("--grpo cannot set off_policy_mask_threshold: use --opsm-start-delta")
    return a


def check_flags(a) -> None:
    """Refusals that need no model: the completion limit (C7: set from the outline count, never above the model's
    budget) and the runaway cap (C11: a capped answer is cut mid-text, so masking it would bias the update toward
    more capping)."""
    from objdel.train import configs

    over = a.grpo
    budget = configs.MAX_COMPLETION[a.model]
    if over.get("max_completion_length", budget) > budget:
        raise ValueError(f"max_completion_length {over['max_completion_length']} is above the model's {budget}")
    if a.runaway_cap and over.get("mask_truncated_completions", True):
        raise ValueError("--runaway-cap needs --grpo '{\"mask_truncated_completions\": false}'")
    if a.runaway_cap and a.group_requests:
        raise ValueError("--runaway-cap and --group-requests exclude each other: one n = G request has one max_tokens")


def numerics_overrides(a) -> dict:
    """The grpo_config arguments of --numerics and --opsm-start-delta (absent: the preset's own threshold)."""
    from objdel.train import configs

    return {"numerics": a.numerics or configs.DEFAULT_NUMERICS} | ({} if a.opsm_start_delta is None
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
    from objdel.train.callbacks import BudgetStop

    return BudgetStop(got["rate_usd_per_h"], a.budget_usd, got["session_start"], fixed_usd=got["fixed_usd"] or 0.0)


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def rows_schedule(root: Path, manifest_sha256: str | None) -> dict | None:
    """The manifest's `schedule` (None for a root without one), after checking the manifest's SHA-256 against the
    pinned value; a schedule root needs it."""
    raw = (Path(root) / "manifest.json").read_bytes()
    schedule = json.loads(raw).get("schedule")
    if schedule is not None and manifest_sha256 is None:
        raise ValueError(f"{root}: a rows root with a schedule needs --rows-manifest-sha256")
    if manifest_sha256 is not None and _sha(raw) != manifest_sha256:
        raise ValueError(f"{root}/manifest.json SHA-256 {_sha(raw)} is not --rows-manifest-sha256 {manifest_sha256}")
    return schedule


def check_rows(keys, schedule: dict) -> None:
    """The rows are the schedule: steps x per_step keys whose newline-joined SHA-256 is ids_sha256 (a dropped,
    swapped or reordered row would shift every later step)."""
    keys, n = list(keys), schedule["steps"] * schedule["per_step"]
    if len(keys) != n:
        raise ValueError(f"rows: {len(keys)} rows, the schedule holds steps x per_step = {n}")
    if _sha("\n".join(keys).encode()) != schedule["ids_sha256"]:
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


def heldout_rows(view: Path, n: int):
    """n val-shape items of an eval view (a seeded draw in record order), as GRPO rows: val-select items never enter
    training's held-out curve."""
    import numpy as np
    from datasets import Dataset, Image

    from objdel.data.splits import val_half
    from objdel.train import dataset

    shape = {r["item_id"] for r in dataset.records(view, "lidc", "val") if val_half(r["group_id"]) == "shape"}
    rows = [r for r in dataset.prompts(view, "lidc", "val") if r["key"] in shape]
    if not rows:
        raise ValueError(f"--heldout {view}: no val-shape RL item")
    if n < len(rows):
        rows = [rows[i] for i in sorted(np.random.default_rng(0).choice(len(rows), n, replace=False))]
    return Dataset.from_list(rows).cast_column("image", Image())


def save_final(model, proc, out: Path) -> None:
    """The unmerged LoRA adapter to out/final_adapter, then the merged model and its processor to out/final."""
    if hasattr(model, "merge_and_unload"):
        model.save_pretrained(out / "final_adapter")
        model = model.merge_and_unload()  # LoRA -> plain weights
    model.save_pretrained(out / "final")
    proc.save_pretrained(out / "final")


def cut_seq_gaps(trainer_dir: Path, step: int) -> None:
    """seq_gaps.f32 back to the first `step` steps (rows a crashed run logged after its last checkpoint go; a fresh
    run, step 0, deletes it), so step i of the file stays optimizer step i after --resume."""
    f, meta = trainer_dir / "seq_gaps.f32", trainer_dir / "seq_gaps.json"
    if not step or not meta.is_file() or not f.is_file():
        f.unlink(missing_ok=True)
        meta.unlink(missing_ok=True)
        return
    m = json.loads(meta.read_text())
    with open(f, "r+b") as fh:
        fh.truncate(min(f.stat().st_size, step * m["rows_per_step"] * 3 * 4))
    m["steps"] = f.stat().st_size // (m["rows_per_step"] * 3 * 4)
    meta.write_text(json.dumps(m) + "\n")


def merge_last(a) -> int:
    """--merge-last: <out>/final from the last checkpoint when the run ended before its own save."""
    from transformers.trainer_utils import get_last_checkpoint

    from experiments.phase2.g_probe import save_json

    out = Path(a.out)
    final = out / "final"
    if (final / "config.json").is_file() and any(final.glob("*.safetensors")) and (out / "rl.json").is_file():
        print(f"{final} exists: nothing to merge")
        return 0
    shutil.rmtree(final, ignore_errors=True)  # a half-written final (the box died inside save_final)
    last = get_last_checkpoint(str(out / "trainer")) if (out / "trainer").is_dir() else None
    if last is None:
        print(f"--merge-last: no checkpoint in {out / 'trainer'}", file=sys.stderr)
        return 1
    from peft import PeftModel
    from transformers import AutoModelForImageTextToText, AutoProcessor

    from objdel.train import configs

    model = PeftModel.from_pretrained(AutoModelForImageTextToText.from_pretrained(a.ckpt), last)
    save_final(model, AutoProcessor.from_pretrained(a.ckpt, **configs.processor_kwargs(a.model, "lidc")), out)
    return save_json(out / "rl.json", {
        "model": a.model, "steps": int(Path(last).name.split("-")[1]), "partial": True, "merged_from": last,
        "ckpt": str(a.ckpt),
        "final": str(final),
    }, keys=("steps", "partial", "merged_from"))


def with_caps(train, a):
    """C11: the rows plus a `cap` column (FastGRPO sends it as each prompt's max_tokens); T_i from the row's own record."""
    from transformers import AutoProcessor

    from objdel.train import configs, dataset, runaway

    tok = AutoProcessor.from_pretrained(a.ckpt, **configs.processor_kwargs(a.model, "lidc")).tokenizer
    caps = runaway.row_caps(train["key"], dataset.records(a.train, "lidc", "train"),
                            lambda text: dataset.target_tokens(tok, text), a.grpo.get("max_completion_length", configs.MAX_COMPLETION[a.model]),
                            a.runaway_cap, runaway.MARGIN)
    return train.add_column("cap", caps)


def main(argv=None) -> int:
    a = parse(argv)
    if a.merge_last:
        return merge_last(a)
    from transformers.trainer_utils import get_last_checkpoint

    from experiments.phase2.g_probe import make_trainer, save_json
    from objdel.train import configs, dataset, serve

    check_flags(a)
    schedule = rows_schedule(a.train, a.rows_manifest_sha256)  # every refusal before any model or engine load
    try:
        stop = session_budget(a)
    except BudgetFlagError as e:
        print(f"train_rl: error: {e}", file=sys.stderr)
        return 2
    out = Path(a.out)
    train = dataset.grpo_rows(a.train, "lidc", "train", oversample=1 if schedule else dataset.DISAGREEMENT_OVERSAMPLE)
    if schedule:
        check_rows(train["key"], schedule)
    if a.runaway_cap:
        train = with_caps(train, a)
    heldout = heldout_rows(a.heldout, a.heldout_items) if a.heldout else None
    last = get_last_checkpoint(str(out / "trainer")) if a.resume and (out / "trainer").is_dir() else None
    if a.resume and last is None:
        print(f"--resume: no checkpoint in {out / 'trainer'}; starting at step 0", file=sys.stderr)
    start = int(Path(last).name.split("-")[1]) if last else 0
    cut_seq_gaps(out / "trainer", start)
    for root in [a.train] + ([a.heldout] if a.heldout else []):
        serve.preload(root, ["lidc"])
    gc.freeze()  # preloaded store: leave it out of the GC scans (was inside serve.fork_pool)
    over = {"save_strategy": "steps", "save_steps": a.save_steps} | numerics_overrides(a)
    over |= {"shuffle_dataset": False} if schedule else {}
    over |= {"eval_strategy": "steps", "eval_steps": a.heldout_every, "num_generations_eval": 8} if a.heldout else {}
    check = (lambda cfg: check_sampler(cfg, schedule, a.steps)) if schedule else None
    from objdel.reward import REWARDS, reward_fingerprint
    from objdel.train.callbacks import LocGuard

    spec = REWARDS[a.reward]
    guard = LocGuard() if a.reward == "v5" else None
    trainer, proc, lr = make_trainer(a, Path(a.ckpt), train, a.steps, eval_dataset=heldout, check=check,
                                     calls=[c for c in (stop, guard) if c], reward=spec, **over)
    t0 = time.time()
    trainer.train(resume_from_checkpoint=last)
    t1 = time.time()
    save_final(trainer.accelerator.unwrap_model(trainer.model), proc, out)
    done = trainer.state.global_step
    log, cfg = trainer.state.log_history, trainer.args
    rows = schedule and {"root": str(a.train), "revision": a.rows_revision, "manifest_sha256": a.rows_manifest_sha256,
                         "ids_sha256": schedule["ids_sha256"]}
    return save_json(out / "rl.json", {
        "model": a.model, "steps": done, "lr": lr, "ckpt": str(a.ckpt), "resumed_from": last,
        "final": str(out / "final"),
        "numerics": numerics_overrides(a)["numerics"],
        "opsm_start_delta": a.opsm_start_delta if a.opsm_start_delta is not None else cfg.off_policy_mask_threshold,
        "beta": cfg.beta, "mask_truncated_completions": cfg.mask_truncated_completions,
        "max_completion_length": cfg.max_completion_length,
        "loss_scale": configs.MAX_COMPLETION[a.model] / cfg.max_completion_length,  # dr_grpo divides by L (C7)
        "reward": {"name": spec.name, "fingerprint": reward_fingerprint(spec), "c_topo": spec.c_topo,
                   "loc_width": getattr(spec, "loc_width", None), "area_gate": getattr(spec, "area_gate", None)},
        "rows": rows, "train_rows": len(train), "budget_stopped": bool(stop and stop.tripped),
        "loc_guard_stopped": bool(guard and guard.tripped),
        "eval": [{"step": h["step"], "eval/iou_tc": h["eval_rewards/iou_tc/mean"],
                  "eval/allmiss_share": h["eval_allmiss_share"]} for h in log if "eval_rewards/iou_tc/mean" in h],
        "timing_s": {"train": t1 - t0, "per_step": (t1 - t0) / max(1, done - start)},
        "log": log,
    }, keys=("steps", "lr", "resumed_from", "final", "budget_stopped"))


if __name__ == "__main__":
    sys.exit(main())
