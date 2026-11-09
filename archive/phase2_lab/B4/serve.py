"""B4 prototype of reward serving: preloaded caches (H-A6-4) and a process pool forked after the preload (H-A6-2).

Design (lab only; nothing here changes src/):
- STORE: module-level {item key: {object: ExpertCache}}, filled once by `preload` in the process that will fork.
- CacheRef(key): a read-only Mapping {object: ExpertCache} that resolves through STORE and pickles as its key
  alone. TRLReward(caches=ref_loader, pool=...) then builds jobs whose pickles are a few hundred bytes (text, flags,
  names, the key, W, H): ExpertCache.__reduce__ would otherwise ship the rater masks and re-run both EDTs in the
  worker (that is why M6's process pool lost to serial).
- Pools:
  * `fork_pool(n)`: ProcessPoolExecutor on the "fork" context, warmed with one no-op task per worker right away.
    With "fork", ProcessPoolExecutor launches every worker at the first submit and never respawns, so the warm-up
    pins the fork to this moment (after the preload, before torch/CUDA initialization), and a dead worker makes
    the next call raise BrokenProcessPool instead of hanging (multiprocessing.Pool.map waits forever when a worker
    is SIGKILLed, e.g. by the OOM killer).
  * `forkserver_pool(n, keys)`: the "forkserver" context with this module and `serve_preload` preloaded in the
    single-threaded fork server, which fills STORE itself (B4_ROOT / B4_KEYS env); safe to create after CUDA init.
  * gc.freeze() before forking, so the cyclic GC in the children does not write to the parent's objects' headers
    (Python docs, gc.freeze).
- Workers set cv2.setNumThreads(0) (one OpenCV thread per worker; the results do not depend on it).

    from experiments.phase2.lab.B4 import serve
    serve.preload(root, keys, threads=2)        # parent, before torch work
    pool = serve.fork_pool(4)
    reward = TRLReward(serve.ref_loader, stop_ids, pool=pool)

References: https://docs.python.org/3/library/multiprocessing.html#contexts-and-start-methods ,
https://docs.python.org/3/library/concurrent.futures.html#processpoolexecutor ,
https://docs.python.org/3/library/gc.html#gc.freeze
"""
import gc
import json
import os
import pickle
import time
from collections.abc import Iterable, Mapping
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing import get_context
from pathlib import Path

import cv2

from objdel.reward import ExpertCache

STORE: dict[str, dict[str, ExpertCache]] = {}
LOADS = {"n": 0}  # cache loads performed in this process (preload or lazy)


class CacheRef(Mapping):
    """{object: ExpertCache} of one item, resolved through STORE; pickles as the key."""

    __slots__ = ("key",)

    def __init__(self, key: str):
        self.key = key

    def __getitem__(self, name):
        return STORE[self.key][name]

    def __iter__(self):
        return iter(STORE[self.key])

    def __len__(self):
        return len(STORE[self.key])

    def __reduce__(self):
        return CacheRef, (self.key,)


def ref_loader(key: str) -> CacheRef:
    if key not in STORE:
        raise KeyError(f"{key}: not preloaded")
    return CacheRef(key)


def rl_digests(root, datasets=("riga", "lidc"), splits=("train",)) -> dict[str, dict[str, str]]:
    """{item key: {object: digest}} of the RL items, as objdel.train.dataset.cache_loader reads them (read here
    without importing objdel.train.dataset, whose `datasets` import the fork server does not need)."""
    out = {}
    for ds in datasets:
        for sp in splits:
            with open(Path(root, ds, f"{sp}.jsonl")) as f:
                out |= {r["item_id"]: r["caches"] for r in map(json.loads, f) if r["rl"]}
    return out


def _load(root, key, digests):
    return key, {n: ExpertCache.load(Path(root, "caches", key, f"{n}.npz"), d) for n, d in digests.items()}


def preload(root, digests: Mapping[str, Mapping[str, str]], threads: int = 2) -> float:
    """Fill STORE with every listed item's caches (digest-checked, as cache_loader); returns the wall time in s."""
    t0 = time.perf_counter()
    if threads <= 1:  # no thread at all (the fork server must stay single-threaded)
        results = (_load(root, *kv) for kv in digests.items())
        for key, cs in results:
            STORE[key] = cs
            LOADS["n"] += len(cs)
        return time.perf_counter() - t0
    with ThreadPoolExecutor(threads) as ex:
        for key, cs in ex.map(lambda kv: _load(root, *kv), list(digests.items())):
            STORE[key] = cs
            LOADS["n"] += len(cs)
    return time.perf_counter() - t0


def _init_worker():
    cv2.setNumThreads(0)
    gc.enable()


def _noop(_=None):
    return os.getpid()


def _noop_chunk(jobs):
    """Same pickling as _score_chunk both ways (jobs in, one small dict per job out), no scoring."""
    return [{"format_ok": 0.0, "error": "", "defects": {}} for _ in jobs]


def fork_pool(n: int) -> ProcessPoolExecutor:
    gc.freeze()
    ex = ProcessPoolExecutor(n, mp_context=get_context("fork"), initializer=_init_worker)
    list(ex.map(_noop, range(n)))  # launches all n workers now (fork context: no dynamic spawning later)
    return ex


def forkserver_pool(n: int, root, keys: Iterable[str]) -> ProcessPoolExecutor:
    """Workers forked from a fork server that preloads the caches of `keys` itself (serve_preload)."""
    ctx = get_context("forkserver")
    kfile = Path(os.environ.get("B4_SCRATCH", "/tmp"), f"b4_keys_{os.getpid()}.json")
    kfile.write_text(json.dumps(sorted(keys)))
    os.environ["B4_ROOT"], os.environ["B4_KEYS"] = str(root), str(kfile)
    ctx.set_forkserver_preload(["experiments.phase2.lab.B4.serve_preload"])
    ex = ProcessPoolExecutor(n, mp_context=ctx, initializer=_init_worker)
    list(ex.map(_noop, range(n)))
    return ex


def _self_mem(delay: float) -> tuple[int, dict]:
    time.sleep(delay)  # keeps this worker busy so the other tasks land on the other workers
    return os.getpid(), mem_mb(os.getpid())


def workers_mem(ex: ProcessPoolExecutor, n: int) -> dict[int, dict]:
    """{pid: mem_mb} reported by each worker itself (macOS denies USS of another process without root)."""
    out = {}
    for _ in range(5):
        out |= dict(ex.map(_self_mem, [0.2] * (2 * n)))
        if len(out) >= n:
            break
    return out


def worker_pids(ex: ProcessPoolExecutor) -> list[int]:
    return sorted(ex._processes)  # private, lab use only


def mem_mb(pid: int) -> dict:
    """RSS and USS (and PSS on Linux) in MB; Linux reads /proc/<pid>/smaps_rollup."""
    import psutil
    p = psutil.Process(pid)
    out = {"rss": p.memory_info().rss / 2**20}
    try:
        full = p.memory_full_info()
        out["uss"] = full.uss / 2**20
        if hasattr(full, "pss"):
            out["pss"] = full.pss / 2**20
    except (psutil.AccessDenied, AttributeError):
        pass
    return out


def cpu_s(pids) -> float:
    """User + system CPU seconds of the given processes (psutil)."""
    import psutil
    return sum(sum(psutil.Process(p).cpu_times()[:2]) for p in pids)


def pickled_bytes(jobs) -> int:
    return len(pickle.dumps(jobs, protocol=pickle.HIGHEST_PROTOCOL))
