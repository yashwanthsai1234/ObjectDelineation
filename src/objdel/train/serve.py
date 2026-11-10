"""Reward serving for the processes that score (plan D14, §3.2): select_ckpt, g_probe, final_eval, k3_harness, GRPO;
never train_sft, which scores nothing. Scoring runs on the CUDA checker (objdel.reward.score_batch) in the one process
that holds the device (D7); no worker pool exists.

1. `preload(root, datasets)` loads every RL cache of those datasets' splits under a prepared root (data/prepared
   or an eval view, same layout) into STORE, before the model loads. Checks: each records file's SHA-256 against
   manifest.json, each cache against its record's digest (ExpertCache.load). LIDC-only reads no RIGA file.
2. `TRLReward(caches=serve.loader(root), stop_ids)` scores on preloaded caches. A `CacheRef(key, root)` names one
   item's caches and pickles as its key and root.

    serve.preload("data/prepared", ("lidc",))
    reward = TRLReward(serve.loader("data/prepared"), stop_ids)

Several roots (prepared val plus eval views) share item ids, so a CacheRef always names its root.
"""

import hashlib
import json
from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path

from objdel.reward import ExpertCache

STORE: dict[str, dict[str, dict[str, ExpertCache]]] = {}  # {resolved root: {item key: {object: ExpertCache}}}


def preload(root, datasets: Iterable[str], threads: int = 4) -> int:
    """Loads the RL caches of every `<dataset>/*.jsonl` under root into STORE; returns the number of items.
    ValueError when a dataset has no records, a records file is not the manifest's or a cache not its record's."""
    root = Path(root).resolve()
    files = json.loads((root / "manifest.json").read_text())["files"]
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
