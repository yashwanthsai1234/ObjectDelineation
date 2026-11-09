"""Reward serving for the processes that score (plan D14, §3.2): select_ckpt, g_probe, k3_harness, Phase 3 GRPO;
never train_sft, which scores nothing. Scoring itself runs on the CUDA checker (objdel.reward.score_batch); the pool
below only spreads the exam's item-level metrics (objdel.eval.metrics) over CPU workers.

1. `preload(root, datasets)` loads every RL cache of those datasets' splits under a prepared root (data/prepared
   or an eval view, same layout) into STORE, before the model loads. Checks: each records file's SHA-256 against
   manifest.json, each cache against its record's digest (ExpertCache.load). LIDC-only reads no RIGA file.
2. `fork_pool(n)`: gc.freeze(), then a ProcessPoolExecutor on the fork context whose n workers are forked at once
   by n warm-up tasks (with fork, the executor starts every worker at the first submit and never respawns), so they
   share STORE copy-on-write. Refused once CUDA is initialized. Workers run OpenCV single-threaded. Workers see
   STORE as it was at the fork: preload every root first; `preload` raises RuntimeError while a pool is alive (until
   its shutdown). Every module a worker runs must be imported before the fork: on macOS a worker that imports a
   module after the fork can crash the pool (S1-11). fork_pool imports SCORERS and `modules` (e.g. the script that
   defines the worker function, when it is not already imported) itself first.
3. `TRLReward(caches=serve.loader(root), stop_ids)` scores on preloaded caches; pool jobs (select_ckpt._score) carry
   `CacheRef(key, root)`, which pickles as its key and root (bytes, not the 96 MB of masks a RIGA batch's
   ExpertCaches pickle to, re-running both EDTs in the worker). A worker uploads the caches it scores to its own
   device context (ExpertCache.device).

A dead worker (OOM kill) makes every pending and later call raise BrokenProcessPool within a second: the run fails
loudly and resumes from a checkpoint. multiprocessing.Pool would hang instead; forkserver or spawn re-import the
training script's torch in every worker (lab B4 §2).

    serve.preload("data/prepared", ("lidc",))
    pool = serve.fork_pool(4)                 # before the model / CUDA; exam metrics only
    reward = TRLReward(serve.loader("data/prepared"), stop_ids)

Several roots (prepared val plus eval views) share item ids, so a CacheRef always names its root.

References: https://docs.python.org/3/library/concurrent.futures.html#processpoolexecutor ,
https://docs.python.org/3/library/gc.html#gc.freeze , archive/phase2_lab/B4/REPORT.md.
"""

import gc
import hashlib
import importlib
import json
import os
import sys
import weakref
from collections.abc import Callable, Iterable, Mapping
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from functools import partial
from multiprocessing import get_context
from pathlib import Path

import cv2

from objdel.reward import ExpertCache

STORE: dict[str, dict[str, dict[str, ExpertCache]]] = {}  # {resolved root: {item key: {object: ExpertCache}}}
SCORERS = ("objdel.reward", "objdel.eval.metrics", "objdel.eval.probes")  # what pool workers run; imported pre-fork
_POOLS = weakref.WeakSet()  # pools forked from STORE and not shut down (BK2-11)


def preload(root, datasets: Iterable[str], threads: int = 4) -> int:
    """Loads the RL caches of every `<dataset>/*.jsonl` under root into STORE; returns the number of items.
    ValueError when a dataset has no records, a records file is not the manifest's or a cache not its record's;
    RuntimeError while a pool from fork_pool is alive (its workers would keep the old STORE, BK2-11)."""
    if _POOLS:
        raise RuntimeError("preload after fork_pool: its workers keep STORE as it was at the fork; preload every root "
                           "first, or shut the pool down")
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


class _Pool(ProcessPoolExecutor):
    def shutdown(self, wait=True, *, cancel_futures=False):
        _POOLS.discard(self)
        super().shutdown(wait, cancel_futures=cancel_futures)


CGROUP = Path("/sys/fs/cgroup")
POOL_CAP = 16  # each worker holds a forked copy of the preloaded expert caches


def quota(cgroup: Path) -> float | None:
    """The container's CPU quota in CPUs: cgroup v2 cpu.max, else v1 cpu.cfs_quota_us / cpu.cfs_period_us (the P1
    Vast box); None when unlimited or unreadable."""
    try:
        q, period = (cgroup / "cpu.max").read_text().split()[:2]
        return None if q == "max" else int(q) / int(period)
    except (OSError, ValueError):
        pass
    for d in ("cpu,cpuacct", "cpu"):
        try:
            q = int((cgroup / d / "cpu.cfs_quota_us").read_text())
            return None if q < 0 else q / int((cgroup / d / "cpu.cfs_period_us").read_text())
        except (OSError, ValueError):
            continue
    return None


def cpus(cgroup: Path = CGROUP) -> int:
    """CPUs this process may use: the affinity mask, bounded by the container's quota (at least 1)."""
    n = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 1
    q = quota(cgroup)
    return n if q is None else min(n, max(1, int(q)))


def pool_size(cgroup: Path = CGROUP) -> int:
    """Scoring workers (EV-5): the usable CPUs minus one for the main process, capped at POOL_CAP."""
    return max(1, min(POOL_CAP, cpus(cgroup) - 1))


def fork_pool(n: int, modules: Iterable[str] = ()) -> ProcessPoolExecutor:
    """n workers forked now, after the preload and before CUDA (RuntimeError once torch has initialized CUDA), with
    SCORERS and `modules` imported first (S1-11)."""
    torch = sys.modules.get("torch")
    if torch is not None and torch.cuda.is_initialized():
        raise RuntimeError("fork_pool after CUDA initialization: create the pool before the model loads (D14)")
    for name in (*SCORERS, *modules):
        importlib.import_module(name)
    gc.freeze()  # the children's GC then leaves the shared objects' headers (and pages) alone
    ex = _Pool(n, mp_context=get_context("fork"), initializer=cv2.setNumThreads, initargs=(0,))
    _POOLS.add(ex)
    for f in [ex.submit(os.getpid) for _ in range(n)]:
        f.result()
    return ex
