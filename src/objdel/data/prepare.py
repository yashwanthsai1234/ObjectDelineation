"""Prepare frames, SFT targets and expert caches once, with hashes (the prepared set the row builders read).

prepare(store, out, prompts=...) reads the item store (objdel.data.store) and writes, per item that has pixels (every
registered dataset the store holds, objdel.core.datasets):
- images/<item_id>.png: the model input of objdel.data.targets.frame (uint8 RGB; PNG is lossless);
- caches/<item_id>/<object>.npz: objdel.verify.reward.cache.ExpertCache.save of the frame's rater stack (items with >= 2
  raters, rl_eligible); loaded lazily by objdel.verify.store.cache_loader, never serialized into a dataset;
- <dataset>/<split>.jsonl, one record per item: item_id, group_id, W, H, box, scale, raters, targets (one canonical
  text per rater, objdel.data.targets.build.build_target with the reward's size rule, None if uncertifiable), rl,
  disagreement (the raters' joint topology answers differ), the SHA-256 of the image pixels and each cache digest;
- manifest.json: SHA-256 of every jsonl (which pins every image and cache through their hashes), per split the
  target-set hash (objdel.data.targets.build.target_set_hash), counts (items, targets, RL items and targets,
  disagreement items, items skipped for lack of pixels, target errors), the store manifest's SHA-256, the reward
  fingerprint, the prompt table the caller passes (the policy's prompts; data names no policy module) and library
  versions. Re-running gives the same manifest (rows sorted by item_id, content hashes).
Items without a model image (Item.image None: LIDC series without CT) are skipped and counted. window_record
builds the same record for one window of an item (RL rows, eval views).

References:
- Python multiprocessing Pool: https://docs.python.org/3/library/multiprocessing.html
"""

import hashlib
import json
import platform
from collections.abc import Iterable, Mapping
from importlib.metadata import version
from itertools import batched
from multiprocessing import get_context
from pathlib import Path

import cv2
import numpy as np

from objdel.core import datasets
from objdel.core.hashing import pixels_sha256
from objdel.core.io import read_json, write_json, write_jsonl
from objdel.data import store
from objdel.data.targets.build import build_target, min_areas, target_set_hash
from objdel.data.targets.frame import frame
from objdel.verify.reward.cache import ExpertCache, rl_eligible
from objdel.verify.reward.score import rater_answers
from objdel.verify.reward.spec import fingerprint

PINS = ("numpy", "numba", "opencv-python-headless", "objdel")  # the library versions every data manifest records
_BATCH = 12  # items framed per pool.map: bounds main-process memory (RIGA: 1.8 GB at 6, 2.4 GB at 32)


def versions(pins=PINS) -> dict[str, str]:
    """The manifests' `versions`: the Python version and the installed version of each package in `pins`."""
    return {"python": platform.python_version(), **{p: version(p) for p in pins}}


def map_batched(pool, fn, jobs, size: int = _BATCH) -> list:
    """pool.map(fn) over `jobs` in chunks of `size`, results in job order; a chunk at a time bounds the parent's
    memory (jobs is consumed lazily)."""
    out = []
    for chunk in batched(jobs, size):
        out += pool.map(fn, chunk)
    return out


def target_triples(recs) -> list[tuple[str, str, str]]:
    """(item_id, rater, target text) of every rater with a target, in record and rater order."""
    return [(r["item_id"], rater, t) for r in recs for rater, t in zip(r["raters"], r["targets"], strict=True) if t]


def set_summary(recs, skipped: int = 0, uniq=None) -> dict:
    """The manifest record of one set: items, skipped_no_image, the target counts, target_set_sha256 and errors,
    rl_items, rl_targets and disagreement_items. `uniq` (rows roots, where a window can repeat in `recs`): the
    distinct records, which give targets, target_set_sha256, target_errors and rl_targets; items, rl_items and
    disagreement_items count `recs`."""
    uniq = recs if uniq is None else uniq
    triples = target_triples(uniq)
    return {"items": len(recs), "skipped_no_image": skipped, "targets": len(triples),
            "target_set_sha256": target_set_hash(triples),
            "target_errors": [e for r in uniq for e in r["errors"]],
            "rl_items": sum(r["rl"] for r in recs),
            "rl_targets": sum(t is not None for r in uniq if r["rl"] for t in r["targets"]),
            "disagreement_items": sum(r["disagreement"] for r in recs)}


def _check_id(out, item_id: str) -> None:
    """ValueError when item_id would write outside out: a NUL, an absolute path or a ".." part (ids of
    untrusted data name files: images/<id>.png, caches/<id>/)."""
    images = Path(out, "images").resolve()
    if "\0" in item_id or Path(item_id).is_absolute() or ".." in Path(item_id).parts or not (
            images / item_id).resolve().is_relative_to(images):
        raise ValueError(f"unsafe id {item_id!r}: it would write outside {out}")


def read_image(path) -> np.ndarray | None:
    """The PNG at `path` as uint8 RGB (H, W, 3) (cv2 decode, BGR -> RGB); None when cv2 cannot read it."""
    im = cv2.imread(str(path), cv2.IMREAD_COLOR)
    return None if im is None else cv2.cvtColor(im, cv2.COLOR_BGR2RGB)


def image_matches(path, sha: str) -> bool:
    """True when the PNG at `path` reads and its pixels hash (pixels_sha256) to `sha`, the record's image_sha256."""
    im = read_image(path)
    return im is not None and pixels_sha256(im) == sha


def write_png(path, rgb: np.ndarray) -> None:
    """Write uint8 RGB (H, W, 3) as a PNG at `path` (cv2, RGB -> BGR), making its folder; OSError when cv2 fails."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(path), cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR)):  # cv2 reports failure by value
        raise OSError(f"{path}: cv2.imwrite failed")


def _prepare(job) -> dict:
    """One framed item -> its record; writes its image and caches under out. The job carries the structure names, so a
    worker never reads the dataset registry."""
    out, item_id, names, group_id, raters, f = job
    _check_id(out, item_id)
    ma = min_areas(f.masks)  # the reward's size rule
    targets, errors = [], []
    for r, rater in enumerate(raters):
        try:
            targets.append(build_target({n: f.masks[n][r] for n in names}, names, None, ma).text)
        except ValueError as e:
            targets.append(None)
            errors.append([item_id, rater, str(e)])
    write_png(Path(out, "images", item_id + ".png"), f.image)
    rec = {"item_id": item_id, "group_id": group_id, "W": f.W, "H": f.H, "box": list(map(int, f.box)),
           "scale": f.scale, "raters": list(raters), "targets": targets, "rl": rl_eligible(raters),
           "disagreement": False, "image_sha256": pixels_sha256(f.image), "caches": {}, "errors": errors}
    if rec["rl"]:
        caches = {n: ExpertCache.build(f.masks[n]) for n in names}
        cdir = Path(out, "caches", item_id)
        cdir.mkdir(parents=True, exist_ok=True)
        rec["caches"] = {n: c.save(cdir / f"{n}.npz") for n, c in caches.items()}
        rec["disagreement"] = len(rater_answers(caches, names)) > 1  # voting raters only
    return rec


def window_record(out, item, window=None, *, f=None) -> dict:
    """The `_prepare` record of one window of `item`, writing its image and caches under `out`.

    window (an objdel.data.windows.Window, RL rows): framed at (dx, dy) and side, keyed by `window.window_id` (so
    two windows of one item never share an image or a cache folder), plus `source_item_id`, `cell` and
    `view = {kind, draw, offset, side, fallback: false, kept_side: false, op: "id"}`.
    window None (eval views): the given frame `f`, keyed by the item's own id; the caller adds its `view`."""
    names = datasets.get(item.dataset).names
    if window is None:
        return _prepare((str(out), item.item_id, names, item.group_id, item.raters, f))
    f = frame(item, offset=(window.dx, window.dy), side=window.side)
    rec = _prepare((str(out), window.window_id, names, item.group_id, item.raters, f))
    return rec | {"source_item_id": item.item_id, "cell": f"{window.cell[0]}{window.cell[1]}",
                  "view": {"kind": window.kind, "draw": window.draw, "offset": [window.dx, window.dy],
                           "side": window.side, "fallback": False, "kept_side": False, "op": "id"}}


def _jobs(store_root, out, dataset, split, skipped):
    for it in store.load(dataset, split, store_root):
        if it.image is None:
            skipped.append(it.item_id)
            continue
        yield str(out), it.item_id, datasets.get(dataset).names, it.group_id, it.raters, frame(it)


def prepare(store_root, out_root, splits: Iterable[str] = ("train", "val"), procs: int = 3, *,
            prompts: Mapping[str, str]) -> dict:
    """Build images, targets and caches of every item with pixels; returns (and writes) the manifest, which stores
    `prompts` (the caller's prompt table: data names no prompt module)."""
    out = Path(out_root)
    out.mkdir(parents=True, exist_ok=True)
    store_manifest = Path(store_root, "manifest.json").read_bytes()
    held = {rel.split("/")[0] for rel in json.loads(store_manifest)["files"]}
    manifest = {"store_manifest_sha256": hashlib.sha256(store_manifest).hexdigest(),
                "reward_fingerprint": fingerprint(), "prompts": dict(prompts),
                "versions": versions(), "files": {}, "sets": {}}
    with get_context("spawn").Pool(procs) as pool:  # a fork of a parent that used CUDA fails in the engine
        for dataset in (d for d in datasets.DATASETS if d in held):
            for split in splits:
                skipped = []
                recs = map_batched(pool, _prepare, _jobs(store_root, out, dataset, split, skipped))
                rel = f"{dataset}/{split}.jsonl"
                (out / dataset).mkdir(exist_ok=True)
                sha = write_jsonl(out / rel, recs)
                manifest["files"][rel] = {"sha256": sha}
                manifest["sets"][f"{dataset}/{split}"] = set_summary(recs, len(skipped))
    write_json(out / "manifest.json", manifest)
    return read_json(out / "manifest.json")  # the manifest as written (JSON types)
