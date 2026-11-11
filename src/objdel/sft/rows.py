"""SFT rows (sft_rows) and the frozen view rows (view_rows) from the prepared set (objdel.policy.rows holds the
records and the image check)."""

import hashlib
import json
from pathlib import Path

import cv2
import numpy as np
from datasets import Dataset, Image
from PIL import Image as PILImage

from objdel.core.io import read_json
from objdel.geom.shape import canonicalize
from objdel.geom.text import SCHEMAS, parse, serialize
from objdel.policy.prompt import _user
from objdel.policy.rows import _image, _verified, records, target_tokens
from objdel.verify.store import voting


def sft_rows(root, dataset: str, split: str, tokenizer=None, max_tokens: int | None = None) -> Dataset:
    """One row per voting rater's target; with tokenizer and max_tokens, answers over the budget are dropped."""
    rows = []
    for r in _verified(root, records(root, dataset, split)):
        voters = voting(root, r, dataset)
        for i, (rater, text) in enumerate(zip(r["raters"], r["targets"], strict=True)):
            if text is None or i not in voters or (max_tokens and target_tokens(tokenizer, text) + 1 > max_tokens):
                continue  # max_tokens + 1: the end-of-turn token
            answer = {"role": "assistant", "content": [{"type": "text", "text": text}]}
            rows.append({"messages": [_user(dataset), answer], "image": _image(root, r["item_id"]),
                         "key": r["item_id"], "rater": rater})
    return Dataset.from_list(rows).cast_column("image", Image())


class _Render:
    """The view rows' transform: each row's spec -> its model image, augment.render on the row's context (C-VIEW).
    A class, so DataLoader workers started with spawn can pickle it."""

    def __init__(self, view_dir):
        self.view_dir = str(view_dir)

    def __call__(self, batch: dict) -> dict:
        from objdel.data.targets import augment  # agent V's module (M1)
        specs = [json.loads(s) for s in batch.pop("spec")]
        batch["image"] = [PILImage.fromarray(augment.render(s, augment.load_context(self.view_dir, s["item_id"])))
                          for s in specs]
        return batch


class _ViewRows(Dataset):
    """A column read by name (`ds["length"]`) skips the transform: transformers' group_by_length sampler reads the
    length column every epoch, which with a plain with_transform renders every image (impl note M2)."""

    def __getitem__(self, key):
        if isinstance(key, str):
            return Dataset.__getitem__(self.with_format(None), key)
        return super().__getitem__(key)


def _inside(ring: np.ndarray, other: np.ndarray) -> bool:
    """Ring lies inside other (parts of a legal target never cross, so any vertex strictly inside decides)."""
    poly = other.astype(np.float32).reshape(-1, 1, 2)
    return any(cv2.pointPolygonTest(poly, (float(x), float(y)), False) > 0 for x, y in ring)


def _no_holes(text: str, names) -> str:
    """A0: every holes list emptied, canonical text (objdel.geom.text round trip). A part inside another part's hole
    (an island) is dropped with the hole: its outer ring would otherwise nest in the filled part's (illegal)."""
    objs = parse(text, names).objs
    keep = {n: [(o, []) for o, _ in parts if not any(_inside(o, q) for q, _ in parts if q is not o)]
            for n, parts in objs.items()}
    return serialize(canonicalize(keep), names)


def _contexts(view_dir: Path, manifest: dict, items: set[str]) -> None:
    """Every row's context file must hash to the manifest's contexts[item_id] (SHA-256 of the npz bytes,
    written by build_views). Checked once at load, before training starts: 0.3 s (warm) to 1.2 s (cold) for the
    7,341 LIDC train contexts (223 MB, impl note M2), instead of failing inside a DataLoader worker hours in."""
    want = manifest.get("contexts")
    if want is None:
        raise ValueError(f"{view_dir}: manifest.json has no context hashes (rebuild the view)")
    for item in sorted(items):
        path = view_dir / "context" / f"{item}.npz"
        if item not in want or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != want[item]:
            raise ValueError(f"{path}: context missing or not the one the manifest records")


def view_rows(view_dir, tokenizer=None, max_tokens: int | None = None, a0: bool = False) -> Dataset:
    """SFT rows of a training view directory (C-ROWS): messages (as sft_rows), image (rendered on access in the
    DataLoader workers), key, rater, row_id, epoch (the view epoch), length (target tokens with a tokenizer, else
    characters; the group_by_length column). The rows are the view's (voting raters only, non-simple targets
    repeated); rows.jsonl must hash to the manifest's view_hash and every context to its contexts entry.
    max_tokens drops answers over the budget (as sft_rows; it needs the tokenizer); a0 empties every holes
    list and keeps the outermost parts."""
    if max_tokens and tokenizer is None:
        raise ValueError("max_tokens is a token budget: pass the tokenizer")
    view_dir = Path(view_dir)
    text = (view_dir / "rows.jsonl").read_bytes()
    manifest = read_json(view_dir / "manifest.json")
    if hashlib.sha256(text).hexdigest() != manifest["view_hash"]:
        raise ValueError(f"{view_dir}: rows.jsonl does not hash to the manifest's view_hash")
    rows = []
    for spec in map(json.loads, text.splitlines()):
        target = spec.pop("target")
        if a0:
            target = _no_holes(target, SCHEMAS[spec["dataset"]])
        n = target_tokens(tokenizer, target) if tokenizer else len(target)
        if max_tokens and n + 1 > max_tokens:
            continue  # + 1: the end-of-turn token
        answer = {"role": "assistant", "content": [{"type": "text", "text": target}]}
        rows.append({"messages": [_user(spec["dataset"]), answer], "spec": json.dumps(spec), "key": spec["item_id"],
                     "rater": spec["rater"], "row_id": spec["row_id"], "epoch": spec["epoch"], "length": n})
    _contexts(view_dir, manifest, {r["key"] for r in rows})
    ds = Dataset.from_list(rows)
    return _ViewRows(ds.data, ds.info, ds.split).with_transform(_Render(view_dir))
