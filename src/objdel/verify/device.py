"""The pixel work of objdel.verify.metrics on the CUDA checker (Implementation.md D7): one objdel.geom.engine.eval_sets
call per item rasterizes its K answers and compares them, per object, with the raters and the vote of the object's
device tables (ExpertCache.device). Only integer counts come back (ItemCounts); metrics() makes every score from them
with the formulas it used on host masks, so the values do not change (B14 of D.spec.md).

Host cuts (D.spec.md §3.2):
- soft Dice at threshold t over n masks: a pixel is in when (masks with it) / n >= t in float64, one rounding, as
  numpy's mean of a bool stack; soft_cut(n, t) = the smallest integer c in [0, n] with float64(c) / n >= t, so
  "count >= soft_cut" is the same test;
- NSD at tolerance tau: nsd_cut(tau) = the largest integer q >= 0 with sqrt(float64(q)) <= tau, so a squared corner
  distance d2 is beyond tau exactly when d2 > nsd_cut(tau);
- duplicates: samples i and j are duplicates when their format and legal flags are equal and, for every object,
  inter_raw[i, j] == area_raw[i] == area_raw[j] (their raw rasters are equal).

The engine splits a call into chunks of at most MAX_SCRATCH_BYTES of scratch; the results do not depend on the split.
"""
import math
import time
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np

from objdel.geom import engine as _gpucheck
from objdel.geom.shape import pack_objects
from objdel.geom.text import GRID, Parsed
from objdel.verify.reward.cache import ExpertCache, device_many

MAX_SCRATCH_BYTES = 256 << 20  # eval_sets scratch cap per chunk [estimate; D1 arena measurement]
eval_sets_s = 0.0  # wall time of this process's eval_sets calls (final_eval resets it per view: D-G14 record)


@dataclass(frozen=True)
class ItemCounts:
    """One object of one item: M = K + k + 1 masks (K drawn samples, k raters, the vote)."""
    inter: np.ndarray       # (M, M) popcount(a & b); the diagonal holds the areas
    inter_raw: np.ndarray   # (K, K) the same between the raw sample rasters (before the legal gate)
    area_raw: np.ndarray    # (K,)
    band_inter: np.ndarray  # (M, M) the same between the boundary-IoU bands (metrics.BIOU_D)
    soft: np.ndarray        # (T, 3) |P_t|, |R_t|, |P_t & R_t| per soft-Dice threshold
    surf: np.ndarray        # (K + k, 3) surfels of the drawn samples and raters by length (sqrt(2)/2, 1, sqrt(2))
    nsd_out: np.ndarray     # (K, k, 2, 3) per (sample, rater), direction, length: surfels beyond the NSD cut
    vote_iou_tc: float      # IoU_tc of the pixel vote of the legal samples (objdel.verify.select.vote_mask)


def soft_cut(n: int, t: float) -> int:
    """The smallest c in [0, n] with float64(c) / n >= t (n + 1 never happens for t <= 1)."""
    c = min(n, max(0, math.ceil(t * n)))
    while c > 0 and (c - 1) / n >= t:
        c -= 1
    while c < n and c / n < t:
        c += 1
    return c


def nsd_cut(tau: float) -> int:
    """The largest integer q >= 0 with sqrt(float64(q)) <= tau. ValueError for a negative, NaN or infinite tau, or one
    above 2^26 px (q beyond 2^52, where float64(q) stops being exact)."""
    if not (0.0 <= tau <= 2.0 ** 26):
        raise ValueError(f"nsd_cut: tau must be in [0, 2^26], got {tau}")
    q = math.floor(tau * tau)
    while q > 0 and math.sqrt(q) > tau:
        q -= 1
    while math.sqrt(q + 1) <= tau:
        q += 1
    return q


def duplicate_keys(format_ok: Sequence[float], legal: Sequence[float], counts: Mapping[str, ItemCounts]) -> list[int]:
    """Per sample, the first sample it duplicates (itself when none): equal flags and equal raw rasters of every
    object (module docstring). Mask equality is transitive, so the classes are those of a hash of the rasters."""
    keys = []
    for i in range(len(format_ok)):
        keys.append(next(j for j in range(i + 1) if j == i or (
            format_ok[j] == format_ok[i] and legal[j] == legal[i]
            and all(c.inter_raw[i, j] == c.area_raw[i] == c.area_raw[j] for c in counts.values()))))
    return keys


def item_counts(parsed: Sequence[Parsed], names: Sequence[str], caches: Mapping[str, ExpertCache], W: int, H: int,
                legal: Sequence[bool], *, thresholds: Sequence[float], nsd_tau: Mapping[str, float],
                band_d: Mapping[str, int], max_scratch_bytes: int = MAX_SCRATCH_BYTES) -> dict[str, ItemCounts]:
    """ItemCounts per object of one item: K parsed answers (an unparsed one draws nothing), their legality (the
    reward's gate), the object caches (one rater count), the soft-Dice thresholds, the NSD tolerance and the band
    width per object. One engine call."""
    global eval_sets_s
    names = tuple(names)
    K, k = len(parsed), caches[names[0]].k
    device_many([caches[n] for n in names])
    objs = [(p.objs[n] if p.ok else [], W, H) for n in names for p in parsed]
    cuts = np.array([[[soft_cut(K, t), soft_cut(k, t)] for t in thresholds]] * len(names), np.int32)
    t0 = time.perf_counter()
    r = _gpucheck.eval_sets(*pack_objects(objs), GRID, np.array([caches[n].device() for n in names], np.int64), K,
                            np.tile(np.asarray(legal, np.uint8), len(names)), cuts,
                            np.array([nsd_cut(nsd_tau[n]) for n in names], np.int64),
                            np.array([band_d[n] for n in names], np.int32), max_scratch_bytes=max_scratch_bytes)
    eval_sets_s += time.perf_counter() - t0
    return {n: ItemCounts(r["inter"][s], r["inter_raw"][s], r["area_raw"][s], r["band_inter"][s], r["soft"][s],
                          r["surf"][s], r["nsd_out"][s], float(r["vote_iou_tc"][s])) for s, n in enumerate(names)}
