"""Training data: prepare frames, SFT targets and expert caches once, then build TRL rows from them.

prepare(store, out) reads the item store (objdel.data.store) and writes, per item that has pixels:
- images/<item_id>.png: the model input of objdel.targets.frame (uint8 RGB; PNG is lossless);
- caches/<item_id>/<object>.npz: objdel.reward.ExpertCache.save of the frame's rater stack (items with >= 2 raters,
  rl_eligible); loaded lazily by cache_loader, never serialized into a dataset;
- <dataset>/<split>.jsonl, one record per item: item_id, group_id, W, H, box, scale, raters, targets (one canonical
  text per rater, objdel.targets.build_target with the reward's size rule, None if uncertifiable), rl,
  disagreement (the raters' joint topology answers differ), the SHA-256 of the image pixels and each cache digest;
- manifest.json: SHA-256 of every jsonl (which pins every image and cache through their hashes), per split the
  target-set hash (objdel.targets.target_set_hash), counts (items, targets, RL items and targets, disagreement
  items, items skipped for lack of pixels, target errors), the store manifest's SHA-256, the reward fingerprint and
  library versions. Re-running gives the same manifest (rows sorted by item_id, content hashes).
Items without a model image (Item.image None: RIGA BinRushed1 reconstructions, break-it B3-1; LIDC series without CT)
are skipped and counted.

Rows (Hugging Face datasets, conversational, one {"type": "image"} placeholder that TRL fills from the "image"
column, decoded lazily from the PNG path):
- sft_rows: messages = [user: image + PROMPTS[dataset], assistant: target text], one row per voting rater with a
  target. A rater whose whole mask is below the reward's size filter abstains from the topology vote (break-it B2-4,
  objdel.reward.rater_answers): its own answer can never be topo-any, so SFT does not teach it (10 / 23,469 LIDC
  train targets; RL keeps the items). All voting raters are rows (not one random rater per epoch): an epoch then
  shows each item's whole rater distribution, which is the topology histogram the RL reward matches; it needs no
  custom sampler. Issue I5: with tokenizer and
  max_tokens (the model's max_completion_length), targets longer than the budget (+1 end-of-turn token) are
  dropped, so SFT never teaches an answer that GRPO generation would truncate.
- grpo_rows: prompt (the same user turn), image, key (item_id, the cache key), W, H, names; only RL-eligible items.
  oversample > 1 repeats disagreement items (plan §10). TRL's RepeatSampler permutes rows uniformly and takes no
  weights, so a repeated row is the simplest sampler-free weighting; the factor is decided in
  experiments/phase1/impl/M7_train.md (DISAGREEMENT_OVERSAMPLE).
Both row builders decode every image they reference and compare its pixels with the record's image_sha256
(break-it B3-5): rows are the only way images reach a trainer, so a swapped or edited PNG cannot be trained on
(ValueError). Cost: about 1 ms per image with 3 threads.
cache_loader checks every cache against the record's digest when loading it (ExpertCache.load(expected_digest));
a cache built under another reward fingerprint raises StaleCacheError ("re-run prepare_rl"), distinct from a
corrupted or swapped file (break-it B2-8).

Flip augmentation (issue I7) is not implemented: Phase 1 runs no training. If Phase 2 adds it, flip the image and
the rater masks and rebuild targets (build_target) and caches (ExpertCache.build) from the flipped masks; never
mirror polygon vertices, because the top-left pixel-centre tie rule is not mirror-symmetric.

References:
- TRL dataset formats (conversational, vision: "image" column + {"type": "image"}): https://huggingface.co/docs/trl/main/en/dataset_formats
- TRL GRPO (image column, extra columns reach reward functions): https://huggingface.co/docs/trl/main/en/grpo_trainer
- datasets Image feature (paths decoded on access): https://huggingface.co/docs/datasets/main/en/image_load
- Python multiprocessing Pool: https://docs.python.org/3/library/multiprocessing.html
- OpenCV Python bindings release the GIL around every call (PyAllowThreads, modules/python/src2/cv2_util.hpp):
  https://github.com/opencv/opencv/blob/4.x/modules/python/src2/cv2_util.hpp
- Plan §4.3 (caches stored with a hash), §4.4 (targets), §10 (oversample disagreement items).
"""
import hashlib
import json
import platform
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from importlib.metadata import version
from itertools import batched
from multiprocessing import Pool
from pathlib import Path

import cv2
import numpy as np
from datasets import Dataset, Image

from objdel.data import store
from objdel.geometry import SCHEMAS
from objdel.reward import ExpertCache, fingerprint, rater_answers, rl_eligible
from objdel.targets import build_target, frame, min_areas, target_set_hash

PROMPTS = {
    "riga": 'Outline the optic disc and the optic cup in this fundus crop. Answer with JSON only: '
            '{"disc":[{"outer":[[x,y],...],"holes":[]}],"cup":[{"outer":[[x,y],...],"holes":[]}]}, '
            'integers 0-1000 relative to the crop, outer rings clockwise.',
    "lidc": 'Outline the lung nodule in this CT crop. Answer with JSON only: '
            '{"nodule":[{"outer":[[x,y],...],"holes":[[[x,y],...]]}]}, '
            'integers 0-1000 relative to the crop, outer rings clockwise, holes counter-clockwise.',
}
DISAGREEMENT_OVERSAMPLE = 4  # decided in experiments/phase1/impl/M7_train.md
_PINS = ("numpy", "numba", "opencv-python-headless", "objdel")
_BATCH = 12  # items framed per pool.map: bounds main-process memory (RIGA: 1.8 GB at 6, 2.4 GB at 32)


def _pixels_sha(image: np.ndarray) -> str:
    return hashlib.sha256(repr(image.shape).encode() + image.tobytes()).hexdigest()


def _prepare(job) -> dict:
    """One framed item -> its record; writes its image and caches under out."""
    out, item_id, dataset, group_id, raters, f = job
    names = SCHEMAS[dataset]
    ma = min_areas(f.masks)  # the reward's size rule
    targets, errors = [], []
    for r, rater in enumerate(raters):
        try:
            targets.append(build_target({n: f.masks[n][r] for n in names}, names, None, ma).text)
        except ValueError as e:
            targets.append(None)
            errors.append([item_id, rater, str(e)])
    path = Path(out, "images", item_id + ".png")
    path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(path), cv2.cvtColor(f.image, cv2.COLOR_RGB2BGR))
    rec = {"item_id": item_id, "group_id": group_id, "W": f.W, "H": f.H, "box": list(map(int, f.box)),
           "scale": f.scale, "raters": list(raters), "targets": targets, "rl": rl_eligible(raters),
           "disagreement": False, "image_sha256": _pixels_sha(f.image), "caches": {}, "errors": errors}
    if rec["rl"]:
        caches = {n: ExpertCache.build(f.masks[n]) for n in names}
        cdir = Path(out, "caches", item_id)
        cdir.mkdir(parents=True, exist_ok=True)
        rec["caches"] = {n: c.save(cdir / f"{n}.npz") for n, c in caches.items()}
        rec["disagreement"] = len(rater_answers(caches, names)) > 1  # voting raters only (B2-4)
    return rec


def _jobs(store_root, out, dataset, split, skipped):
    for it in store.load(dataset, split, store_root):
        if it.image is None:
            skipped.append(it.item_id)
            continue
        yield str(out), it.item_id, dataset, it.group_id, it.raters, frame(it)


def prepare(store_root, out_root, splits: Iterable[str] = ("train", "val"), procs: int = 3) -> dict:
    """Build images, targets and caches of every item with pixels; returns (and writes) the manifest."""
    out = Path(out_root)
    out.mkdir(parents=True, exist_ok=True)
    manifest = {"store_manifest_sha256": hashlib.sha256(Path(store_root, "manifest.json").read_bytes()).hexdigest(),
                "reward_fingerprint": fingerprint(), "prompts": PROMPTS,
                "versions": {"python": platform.python_version(), **{p: version(p) for p in _PINS}},
                "files": {}, "sets": {}}
    with Pool(procs) as pool:
        for dataset in SCHEMAS:
            for split in splits:
                skipped, recs = [], []
                for jobs in batched(_jobs(store_root, out, dataset, split, skipped), _BATCH):
                    recs += pool.map(_prepare, jobs)
                rel = f"{dataset}/{split}.jsonl"
                (out / dataset).mkdir(exist_ok=True)
                text = "".join(json.dumps(r, sort_keys=True) + "\n" for r in recs)
                (out / rel).write_text(text)
                triples = [(r["item_id"], rater, t) for r in recs
                           for rater, t in zip(r["raters"], r["targets"], strict=True) if t]
                manifest["files"][rel] = {"sha256": hashlib.sha256(text.encode()).hexdigest()}
                manifest["sets"][f"{dataset}/{split}"] = {
                    "items": len(recs), "skipped_no_image": len(skipped), "targets": len(triples),
                    "target_set_sha256": target_set_hash(triples),
                    "target_errors": [e for r in recs for e in r["errors"]],
                    "rl_items": sum(r["rl"] for r in recs),
                    "rl_targets": sum(t is not None for r in recs if r["rl"] for t in r["targets"]),
                    "disagreement_items": sum(r["disagreement"] for r in recs)}
    doc = json.dumps(manifest, indent=1, sort_keys=True) + "\n"
    (out / "manifest.json").write_text(doc)
    return json.loads(doc)  # the manifest as written (JSON types)


def records(root, dataset: str, split: str) -> list[dict]:
    with open(Path(root, dataset, f"{split}.jsonl")) as f:
        return [json.loads(line) for line in f]


def _user(dataset: str) -> dict:
    return {"role": "user", "content": [{"type": "image"}, {"type": "text", "text": PROMPTS[dataset]}]}


def _image(root, item_id: str) -> str:
    return str(Path(root, "images", item_id + ".png").resolve())


def _check(root, rec: dict) -> None:
    im = cv2.imread(_image(root, rec["item_id"]), cv2.IMREAD_COLOR)
    if im is None or _pixels_sha(cv2.cvtColor(im, cv2.COLOR_BGR2RGB)) != rec["image_sha256"]:
        raise ValueError(f"{rec['item_id']}: prepared image missing or not the one recorded (image_sha256)")


def _verified(root, recs: list[dict]) -> list[dict]:
    """recs, after checking every image's pixels against its record (cv2 decodes without the GIL: 3 threads)."""
    with ThreadPoolExecutor(3) as pool:
        list(pool.map(lambda r: _check(root, r), recs))
    return recs


def _tokens(tokenizer, text: str) -> int:
    return len(tokenizer(text, add_special_tokens=False)["input_ids"])


def _caches(root, key: str, digests: dict[str, str]) -> dict[str, ExpertCache]:
    return {n: ExpertCache.load(Path(root, "caches", key, f"{n}.npz"), d) for n, d in digests.items()}


def _voting(root, rec: dict, dataset: str) -> set[int]:
    """Indices of the raters with a topology vote (objdel.reward.rater_answers; all raters of an item without
    caches, which has one rater)."""
    if not rec["caches"]:
        return set(range(len(rec["raters"])))
    return {i for rs in rater_answers(_caches(root, rec["item_id"], rec["caches"]), SCHEMAS[dataset]).values()
            for i in rs}


def sft_rows(root, dataset: str, split: str, tokenizer=None, max_tokens: int | None = None) -> Dataset:
    """One row per voting rater's target (B2-4); with tokenizer and max_tokens, answers over the budget are dropped
    (I5)."""
    rows = []
    for r in _verified(root, records(root, dataset, split)):
        voting = _voting(root, r, dataset)
        for i, (rater, text) in enumerate(zip(r["raters"], r["targets"], strict=True)):
            if text is None or i not in voting or (max_tokens and _tokens(tokenizer, text) + 1 > max_tokens):
                continue  # max_tokens + 1: the end-of-turn token
            answer = {"role": "assistant", "content": [{"type": "text", "text": text}]}
            rows.append({"messages": [_user(dataset), answer], "image": _image(root, r["item_id"]),
                         "key": r["item_id"], "rater": rater})
    return Dataset.from_list(rows).cast_column("image", Image())


def grpo_rows(root, dataset: str, split: str, oversample: int = 1) -> Dataset:
    """RL prompts: RL-eligible items, disagreement items repeated `oversample` times."""
    rows = [{"prompt": [_user(dataset)], "image": _image(root, r["item_id"]), "key": r["item_id"], "W": r["W"],
             "H": r["H"], "names": list(SCHEMAS[dataset])}
            for r in _verified(root, [r for r in records(root, dataset, split) if r["rl"]])
            for _ in range(oversample if r["disagreement"] else 1)]
    return Dataset.from_list(rows).cast_column("image", Image())


def cache_loader(root, maxsize: int = 64) -> Callable[[str], dict[str, ExpertCache]]:
    """key -> {object: ExpertCache}, loaded on first use (LRU; a RIGA item's caches hold ~9 MB) and checked against
    the digests in the prepared records (ValueError; StaleCacheError after a reward fingerprint change)."""
    root = Path(root)
    digests = {r["item_id"]: r["caches"] for ds in SCHEMAS for p in sorted(root.glob(f"{ds}/*.jsonl"))
               for r in records(root, ds, p.stem) if r["rl"]}

    @lru_cache(maxsize)
    def load(key: str) -> dict[str, ExpertCache]:
        return _caches(root, key, digests[key])

    return load
