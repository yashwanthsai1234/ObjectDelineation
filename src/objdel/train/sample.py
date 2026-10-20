"""Completions to the evaluation harness: experiments/phase1/eval_run.py reads JSONL rows
{"item_id", "samples": [text, ...], "truncated": [bool, ...]} (plan §3.1 arrows "GRPO completions -> eval" and
"checkpoint -> K samples per validation item -> eval", break-it B4-3).

- sample(trainer, rows, k): k completions per prompt row through the GRPOTrainer's own generation
  (GRPOTrainer._generate in eval mode: the same chat template, image injection (prepare_multimodal_messages),
  sampling settings and backend, HF generate or vLLM, as training), in calls of whole groups of G = 8 identical
  prompts (TRL's vLLM path generates n = num_generations_eval completions per G-th prompt). Truncation is the
  reward's: objdel.train.rewards.truncated on rewards.stop_ids(model, processing class), TRL's EOS + pad rule.
- export(paths, keys): TRL's completions tables (log_completions=True writes
  <output_dir>/completions/completions_<step>.parquet with every log_extra column; TRLReward logs key, truncated
  and group) as rows, one per run of consecutive groups of one key in a table. TRL keeps the last
  generation_batch_size completions and writes the table at every log; an evaluation's log writes the same
  step's file again, holding held-out completions (possibly after the step's training ones), so pass `keys`
  (e.g. the validation ids) to keep one split.
- write(rows, path): the JSONL file.

References:
- TRL 1.16 GRPOTrainer: _generate, _generate_and_score_completions (prepare_multimodal_messages, is_truncated),
  log (completions table, log_extra columns, deque of generation_batch_size):
  https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_trainer.py
- TRL GRPO docs (log_completions, log_extra): https://huggingface.co/docs/trl/main/en/grpo_trainer
- TRL vLLM generation (one request with n = num_generations per G-th prompt):
  https://github.com/huggingface/trl/blob/v1.16.0/trl/generation/vllm_generation.py
"""
import json
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
from trl.data_utils import prepare_multimodal_messages

from objdel.reward import REWARD_V4

from . import rewards

G = REWARD_V4.group


def sample(trainer, rows: Iterable[dict], k: int) -> list[dict]:
    """k completions (a multiple of G) per row of objdel.train.dataset.grpo_rows (prompt, image, key)."""
    if k <= 0 or k % G:
        raise ValueError(f"k = {k}: sample whole groups of {G}")
    stop = rewards.stop_ids(trainer.model, trainer.processing_class)  # the reward's rule
    n = max(G, trainer.args.per_device_eval_batch_size // G * G)  # completions per generation call
    was_training = trainer.model.training
    trainer.model.eval()  # TRL then generates as in held-out evaluation (num_generations_eval)
    out = []
    try:
        for r in rows:
            prompt = r["prompt"]
            if r.get("image") is not None:
                prompt = prepare_multimodal_messages(prompt, images=[r["image"]])
            texts, truncated = [], []
            for start in range(0, k, n):
                _, ids, _, completions, *_ = trainer._generate([prompt] * min(n, k - start))
                texts += [rewards._text(c) for c in completions]
                truncated += [rewards.truncated(x, stop) for x in ids]
            out.append({"item_id": r["key"], "samples": texts, "truncated": truncated})
    finally:
        trainer.model.train(was_training)
    return out


def export(paths: Iterable, keys: Iterable[str] | None = None) -> list[dict]:
    """Rows of TRL completions tables, in file and table order; ValueError unless the table is whole groups."""
    keys = None if keys is None else set(keys)
    rows = []
    for path in paths:
        df = pd.read_parquet(path, columns=["step", "completion", "key", "truncated", "group"])
        last = None
        for i in range(0, len(df), G):
            b = df.iloc[i:i + G]
            if len(b) != G or b["key"].nunique() != 1 or b["group"].nunique() != 1:
                raise ValueError(f"{path}: rows {i}..{i + len(b) - 1} are not one group of {G} completions")
            key = str(b["key"].iloc[0])
            if keys is not None and key not in keys:
                last = None
                continue
            texts = [rewards._text(c if isinstance(c, str) else list(c)) for c in b["completion"]]
            if last is None or last["item_id"] != key:
                last = {"item_id": key, "step": int(b["step"].iloc[0]), "samples": [], "truncated": []}
                rows.append(last)
            last["samples"] += texts
            last["truncated"] += [bool(t) for t in b["truncated"]]
    return rows


def write(rows: Iterable[dict], path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for r in rows:
            f.write(json.dumps(r) + "\n")
