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
- prompt_rows / grpo_rows: prompt (the same user turn), image, key (item_id, the cache key), W, H, names; only
  RL-eligible items (prompt_rows(rl_only=False): every item, for probe views).
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
from PIL import Image as PILImage

from objdel.data import store
from objdel.geometry import SCHEMAS, canonicalize, parse, serialize
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
    if not cv2.imwrite(str(path), cv2.cvtColor(f.image, cv2.COLOR_RGB2BGR)):  # I45: cv2 reports failure by value
        raise OSError(f"{path}: cv2.imwrite failed")
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


def prompts(root, dataset: str, split: str, rl_only: bool = True, oversample: int = 1, subset=None) -> list[dict]:
    """Prompt rows (prompt, image path, key, W, H, names) of the records (images checked against them): RL-eligible
    items only, or with rl_only=False every item (probe views, whose records are rl: false); disagreement items
    repeated `oversample` times; subset = (n, seed): n records drawn once with that seed, in record order."""
    recs = [r for r in records(root, dataset, split) if r["rl"] or not rl_only]
    if subset and subset[0] < len(recs):
        n, seed = subset
        recs = [recs[i] for i in sorted(np.random.default_rng(seed).choice(len(recs), n, replace=False))]
    return [{"prompt": [_user(dataset)], "image": _image(root, r["item_id"]), "key": r["item_id"], "W": r["W"],
             "H": r["H"], "names": list(SCHEMAS[dataset])}
            for r in _verified(root, recs) for _ in range(oversample if r["disagreement"] else 1)]


def prompt_rows(root, dataset: str, split: str, rl_only: bool = True, oversample: int = 1) -> Dataset:
    """prompts() as a datasets.Dataset with decoded images (TRL GRPO)."""
    return Dataset.from_list(prompts(root, dataset, split, rl_only, oversample)).cast_column("image", Image())


def grpo_rows(root, dataset: str, split: str, oversample: int = 1) -> Dataset:
    """RL prompts: RL-eligible items, disagreement items repeated `oversample` times."""
    return prompt_rows(root, dataset, split, True, oversample)


class OpenOnce:
    """GRPO rows whose images are opened once per generation batch (perf: the 16 x 256 opens per step).

    TRL's RepeatSampler hands the DataLoader the same 256 indices (32 images x 8 answers) once per micro-batch, 16 times
    a step, and transformers fetches all 16 before the step (get_batch_samples); TRL generates from the first and reuses
    its answers for the other 15. A datasets.Dataset opens a PNG for every row of every fetch: 4,096 opens for 32 images
    (P1: 12.96 s a step with the GPU idle). OpenOnce reads the rows without their images, opens each distinct image once
    with the column's own decoder (the same pixels), shares it among the rows that name it, and hands the same list back
    while the indices repeat. TRL copies the messages it edits (prepare_multimodal_messages), so a shared row is never
    changed in place. Main process only (dataloader_num_workers=0): each worker would keep its own last stack."""

    def __init__(self, rows: Dataset):
        self.feature = rows.features["image"]
        self.rows = rows.cast_column("image", Image(decode=False))
        self.last: tuple[tuple, list] = ((), [])

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, i: int) -> dict:
        return self.__getitems__([i])[0]

    def __getitems__(self, indices: list[int]) -> list[dict]:
        key = tuple(indices)
        if key and key == self.last[0]:
            return self.last[1]
        opened, out = {}, []
        for row in self.rows.__getitems__(list(indices)):
            v = row["image"]
            k = v["path"] or v["bytes"]
            if k not in opened:
                opened[k] = self.feature.decode_example(v)
            out.append(row | {"image": opened[k]})
        self.last = (key, out)
        return out


def cache_loader(root, maxsize: int | None = 64) -> Callable[[str], dict[str, ExpertCache]]:
    """key -> {object: ExpertCache}, checked against the digests in the prepared records (ValueError;
    StaleCacheError after a reward fingerprint change). The records are the files the manifest lists, each checked
    against its SHA-256 (I45: no glob, so a stray or edited jsonl is never read). maxsize: LRU size (a RIGA item's
    caches hold ~9 MB), loaded on first use; None loads every cache now (preload: the LIDC train set is ~0.33 GB)."""
    root = Path(root)
    digests = {}
    for rel, entry in json.loads((root / "manifest.json").read_text())["files"].items():
        text = (root / rel).read_bytes()
        if hashlib.sha256(text).hexdigest() != entry["sha256"]:
            raise ValueError(f"{root / rel}: not the file the manifest records")
        digests |= {r["item_id"]: r["caches"] for r in map(json.loads, text.splitlines()) if r["rl"]}
    if maxsize is None:
        return {k: _caches(root, k, d) for k, d in digests.items()}.__getitem__

    @lru_cache(maxsize)
    def load(key: str) -> dict[str, ExpertCache]:
        return _caches(root, key, digests[key])

    return load


class _Render:
    """The view rows' transform: each row's spec -> its model image, augment.render on the row's context (C-VIEW).
    A class, so DataLoader workers started with spawn can pickle it."""

    def __init__(self, view_dir):
        self.view_dir = str(view_dir)

    def __call__(self, batch: dict) -> dict:
        from objdel.targets import augment  # agent V's module (M1)
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
    """A0: every holes list emptied, canonical text (objdel.geometry round trip). A part inside another part's hole
    (an island) is dropped with the hole: its outer ring would otherwise nest in the filled part's (illegal, BK1-4)."""
    objs = parse(text, names).objs
    keep = {n: [(o, []) for o, _ in parts if not any(_inside(o, q) for q, _ in parts if q is not o)]
            for n, parts in objs.items()}
    return serialize(canonicalize(keep), names)


def _contexts(view_dir: Path, manifest: dict, items: set[str]) -> None:
    """BK1-3: every row's context file must hash to the manifest's contexts[item_id] (SHA-256 of the npz bytes,
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
    max_tokens drops answers over the budget (I5, as sft_rows; it needs the tokenizer); a0 empties every holes
    list and keeps the outermost parts (plan D11)."""
    if max_tokens and tokenizer is None:
        raise ValueError("max_tokens is a token budget: pass the tokenizer")
    view_dir = Path(view_dir)
    text = (view_dir / "rows.jsonl").read_bytes()
    manifest = json.loads((view_dir / "manifest.json").read_text())
    if hashlib.sha256(text).hexdigest() != manifest["view_hash"]:
        raise ValueError(f"{view_dir}: rows.jsonl does not hash to the manifest's view_hash")
    rows = []
    for spec in map(json.loads, text.splitlines()):
        target = spec.pop("target")
        if a0:
            target = _no_holes(target, SCHEMAS[spec["dataset"]])
        n = _tokens(tokenizer, target) if tokenizer else len(target)
        if max_tokens and n + 1 > max_tokens:
            continue  # + 1: the end-of-turn token
        answer = {"role": "assistant", "content": [{"type": "text", "text": target}]}
        rows.append({"messages": [_user(spec["dataset"]), answer], "spec": json.dumps(spec), "key": spec["item_id"],
                     "rater": spec["rater"], "row_id": spec["row_id"], "epoch": spec["epoch"], "length": n})
    _contexts(view_dir, manifest, {r["key"] for r in rows})
    ds = Dataset.from_list(rows)
    return _ViewRows(ds.data, ds.info, ds.split).with_transform(_Render(view_dir))
