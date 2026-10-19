"""Score a GRPO batch: per-prompt caches, rollouts scored serially, on threads, or on processes.

Results are returned in input order and every per-rollout computation is independent and
order-free, so the output does not depend on the number or kind of workers.
"""
from __future__ import annotations

import os
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import cv2
import numpy as np

from .reward import DEFAULT_CONFIG, ExpertCache, RewardConfig, reward, score

TERM_FIELDS = ("reward", "iou", "band", "topo_any", "b0", "b1", "n_spill", "n_miss", "depth",
               "iou_tc", "topo_freq", "sig_id")


def _row(mask: np.ndarray, cache: ExpertCache, cfg: RewardConfig) -> tuple:
    t = score(mask, cache)
    return (reward(t, cfg), t.iou, t.band[cfg.band_formula], t.topo_any, t.b0, t.b1, t.n_spill, t.n_miss, t.depth,
            t.iou_tc, t.topo_freq, float(hash(t.sig) % 1_000_003))  # tuple-of-int hashes are not salted


def _timed_row(mask, cache, cfg):
    """(terms, wall ns, CPU ns of this thread)."""
    t0, c0 = time.perf_counter_ns(), time.thread_time_ns()
    r = _row(mask, cache, cfg)
    return r, time.perf_counter_ns() - t0, time.thread_time_ns() - c0


def _noop(_):
    return os.getpid()


# ---- process workers: caches are shipped once per worker, masks as packed bits ----
_W: dict = {}


def _init_worker(caches, cfg):
    cv2.setNumThreads(1)
    _W["caches"], _W["cfg"] = caches, cfg


def _proc_task(args):
    i, bits, shape = args
    mask = np.unpackbits(bits, count=shape[0] * shape[1]).reshape(shape).astype(bool)
    return _timed_row(mask, _W["caches"][i], _W["cfg"])


class BatchScorer:
    """Scores (prompt index, mask) pairs against fixed per-prompt caches."""

    def __init__(self, caches: list[ExpertCache], cfg: RewardConfig = DEFAULT_CONFIG,
                 workers: int = 1, backend: str = "thread"):
        cv2.setNumThreads(1)  # parallelism comes from the pool; avoid oversubscription
        self.caches, self.cfg, self.workers, self.backend = caches, cfg, workers, backend
        self._pool = None
        if backend not in ("thread", "process"):
            raise ValueError(backend)
        if workers > 1 and backend == "thread":
            self._pool = ThreadPoolExecutor(workers)
        elif backend == "process":  # also with one worker (cross-process determinism tests)
            self._pool = ProcessPoolExecutor(workers, initializer=_init_worker, initargs=(caches, cfg))
            list(self._pool.map(_noop, range(workers * 4)))  # start and warm up the workers

    def __call__(self, prompt_idx, masks) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Returns (terms (N, len(TERM_FIELDS)) float64, per-rollout wall ns (N,), per-rollout CPU ns (N,))."""
        n = len(masks)
        if self._pool is None:
            out = [_timed_row(m, self.caches[i], self.cfg) for i, m in zip(prompt_idx, masks)]
        elif self.backend == "thread":
            out = list(self._pool.map(lambda a: _timed_row(a[1], self.caches[a[0]], self.cfg),
                                      zip(prompt_idx, masks), chunksize=max(1, n // (4 * self.workers))))
        else:
            args = [(i, np.packbits(m), m.shape) for i, m in zip(prompt_idx, masks)]
            out = list(self._pool.map(_proc_task, args, chunksize=max(1, n // (4 * self.workers))))
        terms = np.array([r for r, _, _ in out], dtype=np.float64)
        wall = np.array([t for _, t, _ in out], dtype=np.int64)
        cpu = np.array([c for _, _, c in out], dtype=np.int64)
        return terms, wall, cpu

    def close(self):
        if self._pool is not None:
            self._pool.shutdown()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def advantages_v2(terms: np.ndarray, group: int = 8, c_topo: float = 2.0, lam: float = 0.5,
                  legal: np.ndarray | None = None) -> np.ndarray:
    """Reward-v2 advantages for a batch scored by BatchScorer (rows ordered prompt by prompt, `group` each)."""
    from .advantages import dm_loo_topology, reward_v2

    f = {n: i for i, n in enumerate(TERM_FIELDS)}
    T = terms.reshape(-1, group, terms.shape[1])
    lg = np.ones(T.shape[:2]) if legal is None else np.asarray(legal, float).reshape(T.shape[:2])
    valid = (T[..., f["topo_any"]] > 0) & (T[..., f["iou_tc"]] >= 0.5)
    r_topo = dm_loo_topology(T[..., f["sig_id"]], T[..., f["topo_freq"]], valid, lam=lam)
    R = np.stack([T[..., f["iou_tc"]], T[..., f["band"]], lg], -1)
    a, _ = reward_v2(R, np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25]), ("iou", "band", "legal"),
                     r_topo, c_topo=c_topo)
    return a.reshape(-1)


def digest(terms: np.ndarray) -> str:
    """SHA-256 of the float64 bytes (bit-for-bit comparison of a batch of scores)."""
    import hashlib

    return hashlib.sha256(np.ascontiguousarray(terms, dtype=np.float64).tobytes()).hexdigest()
