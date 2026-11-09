"""Write eval_run JSONL ({"item_id", "samples", "truncated"}) from a checkpoint or from a GRPO run's completions.

    # checkpoint -> K samples per validation item (TRL's own generation path, objdel.train.sample.sample)
    .venv/bin/python -I -m experiments.phase1.sample_eval --checkpoint data/runs/x/sft --model smolvlm-256m \
        --dataset lidc --out data/runs/x/val_k64.jsonl [--k 64] [--split val] [--limit 64] [--prepared data/prepared]
    # GRPO completions tables (log_completions) -> JSONL, keeping the items of one split (objdel.train.sample.export)
    .venv/bin/python -I -m experiments.phase1.sample_eval --completions data/runs/x/grpo/completions \
        --out data/runs/x/heldout.jsonl [--split val]

Then: .venv/bin/python -I -m experiments.phase1.eval_run --completions <out> --out <report dir>.
The checkpoint must be a full model directory (merge LoRA adapters first). Sampling uses grpo_config's settings
(temperature 1, max_completion_length, vLLM for LFM / Qwen); the Phase 2 gates need K a multiple of 8.
Reads untrusted data: run with python -I.

References:
- objdel.train.sample (TRL GRPOTrainer._generate, completions tables): TRL (version pinned in pyproject.toml)
  https://github.com/huggingface/trl/blob/v1.17.1/trl/trainer/grpo_trainer.py ; log_completions:
  https://huggingface.co/docs/trl/main/en/grpo_trainer
- transformers 5.20 utils.is_torch_bf16_gpu_available (on CUDA: torch.cuda.is_bf16_supported) and the
  TrainingArguments bf16 check:
  https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/utils/import_utils.py ;
  https://docs.pytorch.org/docs/2.14/generated/torch.cuda.is_bf16_supported.html
"""
import argparse
import tempfile
from pathlib import Path


def bf16() -> bool:
    """bf16 only on a CUDA card that supports it (the predicate TrainingArguments checks there, so a pre-Ampere card
    gets fp32 instead of an error); fp32 elsewhere, the Mac included."""
    import torch
    return torch.cuda.is_available() and torch.cuda.is_bf16_supported()


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--checkpoint", help="model directory (or hub id) to sample from")
    src.add_argument("--completions", type=Path, help="a GRPO run's completions/ directory")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--prepared", type=Path, default=Path("data/prepared"))
    ap.add_argument("--split", default="val")
    ap.add_argument("--model", help="objdel.train.configs.MODELS key (with --checkpoint)")
    ap.add_argument("--dataset", choices=("riga", "lidc"), help="with --checkpoint")
    ap.add_argument("--k", type=int, default=64)
    ap.add_argument("--limit", type=int, help="first N RL items of the split")
    a = ap.parse_args()

    from objdel.geometry import SCHEMAS
    from objdel.train import dataset as D
    from objdel.train import sample

    if a.completions:
        keys = {r["item_id"] for ds in SCHEMAS for r in D.records(a.prepared, ds, a.split)}
        rows = sample.export(sorted(a.completions.glob("completions_*.parquet")), keys)
    else:
        if not (a.model and a.dataset):
            ap.error("--checkpoint needs --model and --dataset")
        from transformers import AutoModelForImageTextToText, AutoProcessor
        from trl import GRPOTrainer

        from objdel.train import configs
        proc = AutoProcessor.from_pretrained(configs.MODELS[a.model], **configs.processor_kwargs(a.model, a.dataset))
        model = AutoModelForImageTextToText.from_pretrained(a.checkpoint)
        rows_ds = D.grpo_rows(a.prepared, a.dataset, a.split)
        if a.limit:
            rows_ds = rows_ds.select(range(min(a.limit, len(rows_ds))))
        with tempfile.TemporaryDirectory() as tmp:
            args = configs.grpo_config(a.model, tmp, reward_weights=None, log_completions=False, bf16=bf16())
            trainer = GRPOTrainer(model=model, reward_funcs=lambda completions, **_: [0.0] * len(completions),
                                  args=args, train_dataset=rows_ds, processing_class=proc)
            rows = sample.sample(trainer, rows_ds, a.k)
    sample.write(rows, a.out)
    print(f"{len(rows)} items, {sum(len(r['samples']) for r in rows)} samples -> {a.out}")


if __name__ == "__main__":
    main()
