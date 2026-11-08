"""Full GRPO run: train a checkpoint for --steps, save every --save-steps, keep the tuned model.

    python -I -m experiments.phase2.train_rl --model qwen3-vl-2b --ckpt runs/r1/branches/e2 --train data/prepared
        --out runs/rl --steps 500 --save-steps 50 [--resume]

The trainer is the G-probe's (g_probe.make_trainer): FastGRPO with OpenOnce and ImageMemo, vLLM colocate drawing,
scoring on the GPU next to the CUDA checker, G = 8, 32 LIDC train prompts per step (disagreement items x4), the
probe LR unless --lr. Every --save-steps it writes <out>/trainer/checkpoint-N: the LoRA adapter, optimizer,
scheduler, RNG and trainer state, so --resume continues from the last one after a lost box. At the end the LoRA
is merged into the base and saved with its processor to <out>/final, a plain model folder (train_sft's branch
layout: select_ckpt and vLLM read it unchanged). <out>/rl.json: steps, LR, timing and TRL's log.
"""
import argparse
import sys
import time
from pathlib import Path


def main(argv=None) -> int:
    from objdel.train import serve

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--model", required=True, help="configs.MODELS key")
    ap.add_argument("--ckpt", type=Path, required=True, help="the SFT checkpoint to tune (a plain model folder)")
    ap.add_argument("--train", type=Path, required=True, help="prepared root with lidc/train.jsonl and caches")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--steps", type=int, required=True, help="optimizer steps in total (a resumed run included)")
    ap.add_argument("--save-steps", type=int, default=50, help="checkpoint every N steps (default 50)")
    ap.add_argument("--lr", type=float, help="learning rate (default: the G-probe's, LoRA or full FT)")
    ap.add_argument("--resume", action="store_true", help="continue from the last <out>/trainer/checkpoint-N")
    ap.add_argument("--procs", type=int, default=serve.pool_size(), help="scoring workers (default: usable CPUs - 1)")
    ap.add_argument("--grpo", default="{}", help="JSON overrides of configs.grpo_config")
    a = ap.parse_args(argv)

    from transformers.trainer_utils import get_last_checkpoint

    from experiments.phase2.g_probe import make_trainer, save_json
    from objdel.train import dataset

    last = get_last_checkpoint(str(a.out / "trainer")) if a.resume and (a.out / "trainer").is_dir() else None
    if a.resume and last is None:
        print(f"--resume: no checkpoint in {a.out / 'trainer'}; starting at step 0", file=sys.stderr)
    serve.preload(a.train, ["lidc"])
    pool = serve.fork_pool(a.procs)  # before the model and CUDA (D14)
    train = dataset.grpo_rows(a.train, "lidc", "train", oversample=dataset.DISAGREEMENT_OVERSAMPLE)
    trainer, proc, lr = make_trainer(a, a.ckpt, train, a.steps, save_strategy="steps", save_steps=a.save_steps)
    t0 = time.time()
    trainer.train(resume_from_checkpoint=last)
    t1 = time.time()
    model = trainer.accelerator.unwrap_model(trainer.model)
    model = model.merge_and_unload() if hasattr(model, "merge_and_unload") else model  # LoRA -> plain weights
    model.save_pretrained(a.out / "final")
    proc.save_pretrained(a.out / "final")
    pool.shutdown()
    done, start = trainer.state.global_step, int(Path(last).name.split("-")[1]) if last else 0
    return save_json(a.out / "rl.json", {
        "steps": done, "lr": lr, "ckpt": str(a.ckpt), "resumed_from": last, "final": str(a.out / "final"),
        "timing_s": {"train": t1 - t0, "per_step": (t1 - t0) / max(1, done - start)},
        "log": trainer.state.log_history,
    }, keys=("steps", "lr", "resumed_from", "final"))


if __name__ == "__main__":
    sys.exit(main())
