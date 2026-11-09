"""Full GRPO run: train a checkpoint for --steps, save every --save-steps, keep the tuned model.

    python -I -m experiments.phase2.train_rl --model qwen3-vl-2b --ckpt runs/r1/branches/e2 --train data/prepared
        --out runs/rl --steps 500 --save-steps 50 [--resume]
    python -I -m experiments.phase2.train_rl --model qwen3-vl-2b --ckpt {rd}/ckpt --train {data}/rows/lidc-train-b1
        --rows-revision <hex40> --rows-manifest-sha256 <hex64> --out {rd}/b1-rows --steps 260 --budget-usd 2.05
    python -I -m experiments.phase2.train_rl --merge-last --model qwen3-vl-2b --ckpt {rd}/ckpt --out {rd}/b1-rows

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
--session-start x --rate-usd-per-h reaches it; run.sh exports OBJDEL_SESSION_START, OBJDEL_RATE_USD_PER_H and
OBJDEL_FIXED_USD from the active run, the flags' defaults.

Held-out curve (off by default): --heldout VIEW evaluates --heldout-items val-shape items (splits.val_half) of that
eval view every --heldout-every steps, 8 answers each; rl.json "eval" holds eval/iou_tc and eval/allmiss_share.

<out>/rl.json: steps, LR, numerics, beta, rows provenance, timing and TRL's log.
"""
import argparse
import gc
import hashlib
import json
import os
import re
import sys
import time
from pathlib import Path

NUMERICS = "bf16-seqmask"  # today's values (bf16, sequence_mask, cap 3) until B2.2 adds --numerics


def _hex(n: int):
    def check(s: str) -> str:
        if not re.fullmatch(f"[0-9a-f]{{{n}}}", s):
            raise argparse.ArgumentTypeError(f"{s!r}: needs exactly {n} lowercase hex digits")
        return s
    return check


def parse(argv=None) -> argparse.Namespace:
    env = os.environ.get
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0], allow_abbrev=False)
    ap.add_argument("--model", required=True, help="configs.MODELS key")
    ap.add_argument("--ckpt", type=Path, required=True, help="the checkpoint to tune (a plain model folder)")
    ap.add_argument("--train", type=Path, help="prepared or rows root with lidc/train.jsonl and caches")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--steps", type=int, help="optimizer steps in total (a resumed run included)")
    ap.add_argument("--save-steps", type=int, default=50, help="checkpoint every N steps (default 50)")
    ap.add_argument("--shared-prefix", action="store_true", help="pack each group's prompt once per chunk (FastGRPO)")
    ap.add_argument("--lr", type=float, help="learning rate (default: the G-probe's, LoRA or full FT)")
    ap.add_argument("--resume", action="store_true", help="continue from the last <out>/trainer/checkpoint-N")
    ap.add_argument("--merge-last", action="store_true",
                    help="no training: merge the last checkpoint into --ckpt as <out>/final unless it exists")
    ap.add_argument("--rows-revision", type=_hex(40), help="HF commit of the rows root (provenance in rl.json)")
    ap.add_argument("--rows-manifest-sha256", type=_hex(64), help="SHA-256 of the rows root's manifest.json")
    ap.add_argument("--budget-usd", type=float, help="stop and save when the session spend reaches this")
    ap.add_argument("--session-start", type=float, default=env("OBJDEL_SESSION_START"),
                    help="unix time the rental started (default: $OBJDEL_SESSION_START)")
    ap.add_argument("--rate-usd-per-h", type=float, default=env("OBJDEL_RATE_USD_PER_H"),
                    help="USD per hour, bandwidth included (default: $OBJDEL_RATE_USD_PER_H)")
    ap.add_argument("--fixed-usd", type=float, default=env("OBJDEL_FIXED_USD", "0"),
                    help="one-off USD of the session (default: $OBJDEL_FIXED_USD, else 0)")
    ap.add_argument("--heldout", type=Path, help="eval view for the held-out curve (val-shape items only)")
    ap.add_argument("--heldout-items", type=int, default=64)
    ap.add_argument("--heldout-every", type=int, default=50)
    ap.add_argument("--grpo", default="{}", help="JSON overrides of configs.grpo_config")
    ap.add_argument("--profile-window", type=int, nargs=2, metavar=("S", "E"),
                    help="NVTX step ranges and a cudaProfilerApi capture of steps S..E (nsys_capture.sh; S >= 2)")
    ap.add_argument("--trace", help="torch.profiler trace of one step to this dir: step S with --profile-window, "
                                    "else step 3 (the nsys fallback)")
    a = ap.parse_args(argv)
    if not a.merge_last and (a.train is None or a.steps is None):
        ap.error("--train and --steps are required (except with --merge-last)")
    if a.budget_usd is not None and a.session_start is None:
        ap.error("--budget-usd needs a session start: --session-start or OBJDEL_SESSION_START (run.sh exports it)")
    if a.budget_usd is not None and a.rate_usd_per_h is None:
        ap.error("--budget-usd needs a rate: --rate-usd-per-h or OBJDEL_RATE_USD_PER_H (run.sh exports it)")
    if a.profile_window and not 2 <= a.profile_window[0] <= a.profile_window[1]:
        ap.error(f"--profile-window {a.profile_window}: needs 2 <= S <= E")
    return a


def budget_stop(a):
    """callbacks.BudgetStop of the budget flags, or None without --budget-usd."""
    if a.budget_usd is None:
        return None
    from objdel.train.callbacks import BudgetStop

    return BudgetStop(a.rate_usd_per_h, a.budget_usd, a.session_start, fixed_usd=a.fixed_usd)


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


def merge_last(a) -> int:
    """--merge-last: <out>/final from the last checkpoint when the run ended before its own save."""
    from transformers.trainer_utils import get_last_checkpoint

    from experiments.phase2.g_probe import save_json

    if (a.out / "final").is_dir():
        print(f"{a.out / 'final'} exists: nothing to merge")
        return 0
    last = get_last_checkpoint(str(a.out / "trainer")) if (a.out / "trainer").is_dir() else None
    if last is None:
        print(f"--merge-last: no checkpoint in {a.out / 'trainer'}", file=sys.stderr)
        return 1
    from peft import PeftModel
    from transformers import AutoModelForImageTextToText, AutoProcessor

    from objdel.train import configs

    model = PeftModel.from_pretrained(AutoModelForImageTextToText.from_pretrained(a.ckpt), last)
    save_final(model, AutoProcessor.from_pretrained(a.ckpt, **configs.processor_kwargs(a.model, "lidc")), a.out)
    return save_json(a.out / "rl.json", {
        "model": a.model, "steps": int(Path(last).name.split("-")[1]), "partial": True, "merged_from": last,
        "ckpt": str(a.ckpt),
        "final": str(a.out / "final"),
    }, keys=("steps", "partial", "merged_from"))


def main(argv=None) -> int:
    a = parse(argv)
    if a.merge_last:
        return merge_last(a)
    from transformers.trainer_utils import get_last_checkpoint

    from experiments.phase2.g_probe import make_trainer, save_json
    from objdel.train import dataset, serve

    schedule = rows_schedule(a.train, a.rows_manifest_sha256)  # every refusal before any model or engine load
    stop = budget_stop(a)
    train = dataset.grpo_rows(a.train, "lidc", "train", oversample=1 if schedule else dataset.DISAGREEMENT_OVERSAMPLE)
    if schedule:
        check_rows(train["key"], schedule)
    heldout = heldout_rows(a.heldout, a.heldout_items) if a.heldout else None
    last = get_last_checkpoint(str(a.out / "trainer")) if a.resume and (a.out / "trainer").is_dir() else None
    if a.resume and last is None:
        print(f"--resume: no checkpoint in {a.out / 'trainer'}; starting at step 0", file=sys.stderr)
    for root in [a.train] + ([a.heldout] if a.heldout else []):
        serve.preload(root, ["lidc"])
    gc.freeze()  # preloaded store: leave it out of the GC scans (was inside serve.fork_pool)
    over = {"save_strategy": "steps", "save_steps": a.save_steps}
    over |= {"shuffle_dataset": False} if schedule else {}
    over |= {"eval_strategy": "steps", "eval_steps": a.heldout_every, "num_generations_eval": 8} if a.heldout else {}
    check = (lambda cfg: check_sampler(cfg, schedule, a.steps)) if schedule else None
    trainer, proc, lr = make_trainer(a, a.ckpt, train, a.steps, eval_dataset=heldout, check=check,
                                     calls=[stop] if stop else [], **over)
    t0 = time.time()
    trainer.train(resume_from_checkpoint=last)
    t1 = time.time()
    save_final(trainer.accelerator.unwrap_model(trainer.model), proc, a.out)
    done, start = trainer.state.global_step, int(Path(last).name.split("-")[1]) if last else 0
    log, cfg = trainer.state.log_history, trainer.args
    rows = schedule and {"root": str(a.train), "revision": a.rows_revision, "manifest_sha256": a.rows_manifest_sha256,
                         "ids_sha256": schedule["ids_sha256"]}
    return save_json(a.out / "rl.json", {
        "model": a.model, "steps": done, "lr": lr, "ckpt": str(a.ckpt), "resumed_from": last,
        "final": str(a.out / "final"),
        "numerics": NUMERICS, "beta": cfg.beta, "mask_truncated_completions": cfg.mask_truncated_completions,
        "rows": rows, "train_rows": len(train), "budget_stopped": bool(stop and stop.tripped),
        "eval": [{"step": h["step"], "eval/iou_tc": h["eval_rewards/iou_tc/mean"],
                  "eval/allmiss_share": h["eval_allmiss_share"]} for h in log if "eval_rewards/iou_tc/mean" in h],
        "timing_s": {"train": t1 - t0, "per_step": (t1 - t0) / max(1, done - start)},
        "log": log,
    }, keys=("steps", "lr", "resumed_from", "final", "budget_stopped"))


if __name__ == "__main__":
    sys.exit(main())
