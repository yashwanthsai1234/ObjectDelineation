"""The doctor-store loader: the expert caches of a prepared root for the processes that score (plan D14, §3.2):
select_ckpt, g_probe, final_eval, k3_harness, GRPO and the verifier service; never train_sft, which scores nothing.
Scoring runs on the CUDA checker (objdel.verify.reward.score.score_batch) in the one process that holds the device
(D7); no worker pool exists.

1. `preload(root, datasets)` loads every RL cache of those datasets' splits under a prepared root (data/prepared
   or an eval view, same layout) into STORE, before the model loads. Checks: each records file's SHA-256 against
   manifest.json, each cache against its record's digest (ExpertCache.load). LIDC-only reads no RIGA file.
2. `TRLReward(caches=store.loader(root), stop_ids)` scores on preloaded caches. A `CacheRef(key, root)` names one
   item's caches and pickles as its key and root.
3. `cache_loader(root)` loads one item's caches on first use (LRU), checked against the manifest's records;
   `voting(root, rec, dataset)` names the raters with a topology vote.

    store.preload("data/prepared", ("lidc",))
    reward = TRLReward(store.loader("data/prepared"), stop_ids)

Several roots (prepared val plus eval views) share item ids, so a CacheRef always names its root.
"""

import hashlib
import json
from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache, partial
from pathlib import Path

from objdel.core.io import read_json
from objdel.geom.text import SCHEMAS
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import rater_answers

STORE: dict[str, dict[str, dict[str, ExpertCache]]] = {}  # {resolved root: {item key: {object: ExpertCache}}}


def preload(root, datasets: Iterable[str], threads: int = 4) -> int:
    """Loads the RL caches of every `<dataset>/*.jsonl` under root into STORE; returns the number of items.
    ValueError when a dataset has no records, a records file is not the manifest's or a cache not its record's."""
    root = Path(root).resolve()
    files = read_json(root / "manifest.json")["files"]
    digests = {}
    for ds in datasets:
        paths = sorted(root.glob(f"{ds}/*.jsonl"))
        if not paths:
            raise ValueError(f"{root}: no {ds} records")
        for p in paths:
            data = p.read_bytes()
            if files.get(p.relative_to(root).as_posix(), {}).get("sha256") != hashlib.sha256(data).hexdigest():
                raise ValueError(f"{p}: SHA-256 differs from manifest.json")
            digests |= {r["item_id"]: r["caches"] for r in map(json.loads, data.splitlines()) if r["rl"]}

    def load(item):
        key, ds = item
        return key, {n: ExpertCache.load(root / "caches" / key / f"{n}.npz", d) for n, d in ds.items()}

    with ThreadPoolExecutor(threads) as ex:
        STORE.setdefault(str(root), {}).update(ex.map(load, digests.items()))
    return len(digests)


class CacheRef(Mapping):
    """{object: ExpertCache} of one preloaded item, resolved through STORE; pickles as (key, root)."""

    __slots__ = ("key", "root")

    def __init__(self, key: str, root: str):
        if key not in STORE.get(root, {}):
            raise KeyError(f"{key}: not preloaded from {root}")
        self.key, self.root = key, root

    def __getitem__(self, name):
        return STORE[self.root][self.key][name]

    def __iter__(self):
        return iter(STORE[self.root][self.key])

    def __len__(self):
        return len(STORE[self.root][self.key])

    def __reduce__(self):
        return CacheRef, (self.key, self.root)


def loader(*roots) -> Callable[[str], CacheRef]:
    """`caches` bound to preloaded roots: a key resolves in the first root that holds it, so a key no root holds
    raises KeyError. train_rl --heldout puts the val view first: the prepared root holds the same val item ids
    (centred frames) beside its train items, and no train id is in a val view."""
    return partial(_first, [str(Path(r).resolve()) for r in roots])


def _first(paths: list[str], key: str) -> CacheRef:
    return CacheRef(key, next((p for p in paths if key in STORE.get(p, {})), paths[0]))


def _caches(root, key: str, digests: dict[str, str]) -> dict[str, ExpertCache]:
    return {n: ExpertCache.load(Path(root, "caches", key, f"{n}.npz"), d) for n, d in digests.items()}


def voting(root, rec: dict, dataset: str) -> set[int]:
    """Indices of the raters with a topology vote (objdel.verify.reward.score.rater_answers; all raters of an item
    without caches, which has one rater)."""
    if not rec["caches"]:
        return set(range(len(rec["raters"])))
    return {i for rs in rater_answers(_caches(root, rec["item_id"], rec["caches"]), SCHEMAS[dataset]).values()
            for i in rs}


def cache_loader(root, maxsize: int | None = 64) -> Callable[[str], dict[str, ExpertCache]]:
    """key -> {object: ExpertCache}, checked against the digests in the prepared records (ValueError;
    StaleCacheError after a reward fingerprint change). The records are the files the manifest lists, each checked
    against its SHA-256 (I45: no glob, so a stray or edited jsonl is never read). maxsize: LRU size (a RIGA item's
    caches hold ~9 MB), loaded on first use; None loads every cache now (preload: the LIDC train set is ~0.33 GB)."""
    root = Path(root)
    digests = {}
    for rel, entry in read_json(root / "manifest.json")["files"].items():
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
