"""train_rl --record-batches (plan C0): each generation batch of a check step as one file that
`torch.load(weights_only=True)` reads, so G-EQ2 and G-EQ3 can replay it on the Mac. A record holds only tensors and
dicts of tensors: prompt and answer ids and masks, the image grids, RowMeta and the rope positions of HF
(`get_rope_index` on the full rows). Pixels are left out (1.6 GB per batch)."""

import dataclasses
from pathlib import Path

import torch

TENSORS = ("prompt_ids", "prompt_mask", "completion_ids", "completion_mask", "image_grid_thw")


def open_dir(path) -> Path:
    """The folder to record into, created. A folder that holds anything raises FileExistsError before a write."""
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f"--record-batches {path}: not empty (nothing is written)")
    path.mkdir(parents=True, exist_ok=True)
    return path


def positions(base, batch: dict) -> torch.Tensor:
    """(3, rows, width) rope positions of HF's get_rope_index (`base`: the Qwen-VL module) on prompt + answer."""
    ids = torch.cat([batch["prompt_ids"], batch["completion_ids"]], 1).cpu()
    mask = torch.cat([batch["prompt_mask"], batch["completion_mask"]], 1).cpu()
    pos, _ = base.get_rope_index(ids, mm_token_type_ids=batch["mm_token_type_ids"].cpu(),
                                 image_grid_thw=batch["image_grid_thw"].cpu(), attention_mask=mask)
    return pos


def batch_record(batch: dict, meta, base) -> dict:
    """The record of a generation batch: TENSORS, num_images, `meta` (RowMeta fields as tensors, when given) and
    `positions` (when `base` is given and the batch has mm_token_type_ids)."""
    rec = {k: batch[k].cpu() for k in TENSORS if k in batch}
    rec["num_images"] = torch.as_tensor(batch["num_images"])
    if meta is not None:
        rec["meta"] = {f.name: torch.from_numpy(getattr(meta, f.name).copy()) for f in dataclasses.fields(meta)}
    if base is not None and "mm_token_type_ids" in batch:
        rec["positions"] = positions(base, batch)
    return rec
