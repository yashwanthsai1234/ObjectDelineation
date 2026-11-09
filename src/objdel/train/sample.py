"""Completions to the evaluation harness: experiments/phase1/eval_run.py reads JSONL rows
{"item_id", "samples": [text, ...], "truncated": [bool, ...]} (plan §3.1 arrows "GRPO completions -> eval" and
"checkpoint -> K samples per validation item -> eval", break-it B4-3).

- sample(trainer, rows, k, items): k completions per prompt row through the GRPOTrainer's own generation
  (GRPOTrainer._generate in eval mode: the same chat template, image injection (prepare_multimodal_messages),
  sampling settings and backend, HF generate or vLLM, as training), `items` rows x k prompts per call, each row's k
  as whole groups of G = 8 identical prompts (TRL's vLLM path generates n = num_generations_eval completions per
  G-th prompt; HF generate pads the batch on the left, as in training). Truncation is the reward's:
  objdel.train.rewards.truncated on rewards.stop_ids(model, processing class), TRL's EOS + pad rule.
- export(paths, keys): TRL's completions tables (log_completions=True writes
  <output_dir>/completions/completions_<step>.parquet with every log_extra column; TRLReward logs key, truncated
  and group) as rows, one per run of consecutive groups of one key in a table. TRL keeps the last
  generation_batch_size completions and writes the table at every log; an evaluation's log writes the same
  step's file again, holding held-out completions (possibly after the step's training ones), so pass `keys`
  (e.g. the validation ids) to keep one split.
- write(rows, path): the JSONL file.
- export_prompts / import_completions: the JSONL bridge to experiments/phase2/sample_vllm.py, which samples with
  vLLM 0.32 in its own venv (plan D13, C-BRIDGE; promoted from archive/phase2_lab/B5
  export_prompts.py and vllm_parity.py). The
  prompt is TRL GRPO's generation prompt as text (prepare_multimodal_messages, then the processor's
  apply_chat_template with add_generation_prompt=True, tokenize=False; B5: vLLM's prompt ids equal the training
  env's in 64 / 64 LFM items and through the renderer for Qwen3.5 and SmolVLM), the image path, the processor
  kwargs minus padding_side (vLLM's mm_processor_kwargs) and every sampling field (configs.sampling, seed,
  stop_token_ids). Completions come back as token ids; truncation is rewards.truncated, the reward's rule.
  sample() (HF generate through TRL) stays the reference and the Qwen3.5 fallback.

References:
- TRL (version pinned in pyproject.toml) GRPOTrainer: _generate, _generate_and_score_completions
  (prepare_multimodal_messages, is_truncated),
  log (completions table, log_extra columns, deque of generation_batch_size):
  https://github.com/huggingface/trl/blob/v1.17.1/trl/trainer/grpo_trainer.py
- TRL GRPO docs (log_completions, log_extra): https://huggingface.co/docs/trl/main/en/grpo_trainer
- TRL vLLM generation (one request with n = num_generations per G-th prompt):
  https://github.com/huggingface/trl/blob/v1.17.1/trl/generation/vllm_generation.py
"""
import json
from collections.abc import Iterable
from pathlib import Path

import pandas as pd
from trl.data_utils import prepare_multimodal_messages

from objdel.reward import REWARD_V4

from . import configs, rewards

G = REWARD_V4.group


def sample(trainer, rows: Iterable[dict], k: int, items: int = 16) -> list[dict]:
    """k completions (a multiple of G) per row of objdel.train.dataset.grpo_rows (prompt, image, key), `items` rows
    per generation call (a 450M model on a 4090 is overhead-bound at one row per call)."""
    if k <= 0 or k % G:
        raise ValueError(f"k = {k}: sample whole groups of {G}")
    stop = rewards.stop_ids(trainer.model, trainer.processing_class)  # the reward's rule
    rows = list(rows)
    was_training = trainer.model.training
    trainer.model.eval()  # TRL then generates as in held-out evaluation (num_generations_eval)
    out = []
    try:
        for start in range(0, len(rows), items):
            chunk = rows[start:start + items]
            prompts = [r["prompt"] if r.get("image") is None else
                       prepare_multimodal_messages(r["prompt"], images=[r["image"]]) for r in chunk]
            _, ids, _, completions, *_ = trainer._generate([p for p in prompts for _ in range(k)])
            for j, r in enumerate(chunk):
                part = slice(j * k, (j + 1) * k)
                out.append({"item_id": r["key"], "samples": [rewards._text(c) for c in completions[part]],
                            "truncated": [rewards.truncated(x, stop) for x in ids[part]]})
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


def _path(image) -> str:
    return image if isinstance(image, str) else image["path"]


def export_prompts(rows, path, model: str, dataset: str, n: int, processor=None, stop_ids=None,
                   seed: int = 0) -> None:
    """rows: prompt, image (a path; grpo_rows' decoded Image column is read as paths), key (unique) -> the prompts
    JSONL. processor: AutoProcessor of the sampled checkpoint (default: MODELS[model] with processor_kwargs)."""
    kw = configs.processor_kwargs(model, dataset)
    if processor is None:
        from transformers import AutoProcessor
        processor = AutoProcessor.from_pretrained(configs.MODELS[model], **kw)
    if hasattr(rows, "cast_column"):  # a datasets.Dataset: paths, not decoded images
        from datasets import Image
        rows = rows.cast_column("image", Image(decode=False))
    sampling = configs.sampling(model) | {"seed": seed}
    if stop_ids is not None:
        sampling["stop_token_ids"] = sorted(stop_ids)
    mm = {k: v for k, v in kw.items() if k != "padding_side"}
    seen, out = set(), []
    for r in rows:
        if r["key"] in seen:
            raise ValueError(f"{r['key']}: exported twice (ids are item keys)")
        seen.add(r["key"])
        image = _path(r["image"])
        text = processor.apply_chat_template([prepare_multimodal_messages(r["prompt"], images=[image])],
                                             add_generation_prompt=True, tokenize=False)
        out.append({"id": r["key"], "prompt": text[0] if isinstance(text, list) else text, "image": image,
                    "processor_kwargs": mm, "sampling": sampling, "n": n})
    write(out, path)


def import_completions(path, stop_ids, n: int = G, ids: Iterable[str] | None = None) -> list[dict]:
    """sample_vllm's completions -> eval_run rows {"item_id", "samples", "truncated", "tokens"} (tokens: completion
    length in token ids, for eval.metrics' token stats; samples in k order, items in order of first appearance).
    ValueError unless every id has exactly the exported n samples k = 0 .. n-1, each once, and (with ids, the
    exported keys) the ids are exactly those (BK2-3, BK3-7: a cut or doubled file would change K silently)."""
    by: dict[str, dict[int, dict]] = {}
    with open(path) as f:
        for c in map(json.loads, f):
            if by.setdefault(c["id"], {}).setdefault(c["k"], c) is not c:
                raise ValueError(f"{c['id']}: sample k = {c['k']} twice")
    if ids is not None and set(by) != set(ids):
        raise ValueError(f"{path}: missing: {sorted(set(ids) - set(by))[:5]}, not exported: "
                         f"{sorted(set(by) - set(ids))[:5]}")
    rows = []
    for key, cs in by.items():
        if sorted(cs) != list(range(n)):
            raise ValueError(f"{key}: {len(cs)} samples k = {sorted(cs)[:3]}.., want {n} samples k = 0 .. {n - 1}")
        cs = [cs[k] for k in range(n)]
        rows.append({"item_id": key, "samples": [c["text"] for c in cs],
                     "truncated": [rewards.truncated(c["token_ids"], stop_ids) for c in cs],
                     "tokens": [len(c["token_ids"]) for c in cs]})
    return rows
