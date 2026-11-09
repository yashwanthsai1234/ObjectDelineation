"""Per-prompt expert cache: everything reward v4 needs that depends only on the rater masks.

Built once per (item, object) offline and reused by every rollout of every step:
- vote (strict majority), union U, core I (intersection), contour length P of U;
- the two exact Euclidean distance transforms of the egg-yolk band: distance to U (spill) and to the
  complement of I (core miss);
- tau = max(1 px, P95 of the pooled symmetric rater-pair contour distances), as Nikolov et al. set the
  surface-Dice tolerance from inter-rater variation;
- T5, the self-calibrated F5 tolerance: each rater j is scored against the band of the other m - 1 raters,
  q_j = 2 D_j / P_j (mean squared excursion beyond the band), and with S_m the second moment of a new
  expert's excursion beyond the [min, max] band of m i.i.d. N(0, 1) raters,
      T5^2 = max(mean_j(q_j) * S_m / S_(m-1) / x*, (tau_floor * g(m))^2),   x* = 1 / s* - 1,
  so a new expert's expected F5 = 1 / (1 + q / T5^2) is s* = 0.95 at any object size and rater count;
  S_m is a property of N(0, 1), frozen as constants for m <= 8 (below);
- per-rater size-filtered topology signatures (the rater topology histogram f), rater areas (A_min, A_max)
  and the size rule min_area(|U|).

Abstaining raters (B2-4): a rater whose every component is below the size filter has the empty signature (),
i.e. no topology at the reward's resolution. LIDC readers outline every nodule >= 3 mm on each slice where it
appears (Armato et al. 2011), so such a mask is the cross-section of a real nodule on an end slice, not an answer
"nothing here"; as a topology vote it would let a <= 4 px speck match it. Such a rater abstains: it is left out of
the topology histogram and of IoU_tc's rater set (an annotator who gives no label is dropped from the label
likelihood, as in Dawid-Skene with missing labels), and stays in the band, the vote and the area range. When every
rater abstains, all of them vote (the tiny object's unanimous answer is "below the filter").

Items with one rater are excluded from RL (`rl_eligible`): their band has no width and their topology set is
a single answer. `build` refuses them.

Bit-for-bit reproducibility across CPUs (M9 F1, F2; impl/M9_x86.md, impl/S2.md). The EDT is scipy's exact Maurer
transform: integer squared distances, a float64 sqrt, rounded to float32. That is the correctly rounded float32
sqrt (double rounding is innocuous for sqrt with 53 >= 2 * 24 + 2 bits), equal to OpenCV's own trueDistTrans on
arm64 (103,311 fixture masks, impl/S2.md). OpenCV's x86 wheel routes DIST_MASK_PRECISE to Intel IPP, which is not
correctly rounded and varies between runs (opencv#24082), so it is not used. S_m from quad over scipy.stats differs
by 1 ULP between arm64 and x86 (m = 3, 7), so S_1..S_8 are the reference values (equal to quad to 1e-12, tested); a
larger rater count falls back to quad.

Storage: `save` writes the rater masks, tau and T5 (the only costly values) to a compressed .npz together with
the reward fingerprint (spec.fingerprint: CACHE_SPEC, size rule, numpy / scipy / OpenCV versions) and a SHA-256 of
the content. `load` refuses a corrupted file (ValueError), a file whose digest differs from `expected_digest`
(ValueError), and a file built under another fingerprint (StaleCacheError, a ValueError) unless `rebuild=True`, which
rebuilds tau and T5 from the verified rater masks.

References:
- S. Nikolov et al., "Clinically applicable segmentation of head and neck anatomy for radiotherapy" (surface
  Dice; tolerance from inter-rater variation), JMIR 2021. https://arxiv.org/abs/1809.04430 ;
  https://github.com/google-deepmind/surface-distance
- A. G. Cohn, N. M. Gotts, "The 'egg-yolk' representation of regions with indeterminate boundaries", in Geographic
  Objects with Indeterminate Boundaries (GISDATA 2), Taylor & Francis, 1996. https://www.researchgate.net/publication/2526626
- S. K. Warfield, K. H. Zou, W. M. Wells, "Simultaneous truth and performance level estimation (STAPLE)",
  IEEE TMI 23(7), 2004. https://doi.org/10.1109/TMI.2004.828354
- scipy.ndimage.distance_transform_edt (exact EDT; C. R. Maurer, R. Qi, V. Raghavan, "A linear time algorithm for
  computing exact Euclidean distance transforms of binary images in arbitrary dimensions", IEEE TPAMI 25(2), 2003;
  squared distances, then np.sqrt):
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html ;
  https://github.com/scipy/scipy/blob/main/scipy/ndimage/src/ni_morphology.c ;
  https://github.com/scipy/scipy/blob/main/scipy/ndimage/_morphology.py
- Double rounding is innocuous for sqrt when the intermediate precision is about twice the final one (binary64 ->
  binary32): S. A. Figueroa, "When is double rounding innocuous?", ACM SIGNUM Newsletter 30(3), 1995,
  https://doi.org/10.1145/221332.221334 ; P. Roux, "Innocuous double rounding of basic arithmetic operations",
  J. Formalized Reasoning 7(1), 2014, https://jfr.unibo.it/article/view/4359
- OpenCV's x86 IPP distance transform (used for DIST_MASK_PRECISE when single-threaded or < 2^14 px) and its
  run-to-run rounding: https://github.com/opencv/opencv/blob/5.x/hal/ipp/src/distancetransform_ipp.cpp ;
  https://github.com/opencv/opencv/issues/24082 ; partial fix https://github.com/opencv/opencv/pull/27432
- P. F. Felzenszwalb, D. P. Huttenlocher, "Distance transforms of sampled functions", Theory of Computing 8,
  2012 (OpenCV's own exact DIST_L2 + DIST_MASK_PRECISE path, equal to scipy's output):
  https://theoryofcomputing.org/articles/v008a019/ ;
  https://github.com/opencv/opencv/blob/4.x/modules/imgproc/include/opencv2/imgproc.hpp
- S. G. Armato III et al., "The Lung Image Database Consortium (LIDC) and Image Database Resource Initiative
  (IDRI)", Med. Phys. 38(2), 2011 (nodules >= 3 mm outlined on every slice, < 3 mm marked by a centroid only):
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3041807/
- Missing labels in annotator models (only observed labels enter the likelihood): Stan User's Guide, "Data coding
  models" (Dawid-Skene), https://mc-stan.org/docs/2_32/stan-users-guide/data-coding-models.html ; R. Passonneau,
  B. Carpenter, "The benefits of a model of annotation", TACL 2014, https://aclanthology.org/Q14-1025/
- Lab derivation and evidence: archive/phase1_lab/t_reward/rewardcore/band.py, results/ITERATION3.md (R5).
"""
import hashlib
import math
import os
import weakref
from collections import defaultdict
from functools import cache, cached_property
from itertools import combinations
from pathlib import Path

import cv2
import numpy as np
from scipy import integrate, ndimage, stats

from objdel import _gpucheck
from objdel.checker.topology import min_area

from .spec import CACHE_SPEC, fingerprint

_CROSS = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], np.uint8)
_BIG = np.float32(1e6)  # distance when no target pixel exists in the crop
_P95_PAIR_OVER_SIGMA = 1.959963984540054 * np.sqrt(2.0)  # P95 of |X - Y|, X, Y ~ N(0, 1)
_X_STAR = 1.0 / CACHE_SPEC.expert_score - 1.0
# S_m, m = 1..8 (module docstring): _excursion_quad's values on the reference arm64 Mac (the stored caches' T5)
S_M = {m: float.fromhex(h) for m, h in enumerate((
    "0x1.0000000000001p+1", "0x1.9062e2bc50596p-1", "0x1.cb7070234b5dap-2", "0x1.35e135a79b622p-2",
    "0x1.c92c7ee05f565p-3", "0x1.64f6947907d49p-3", "0x1.21ded7fbee6b0p-3", "0x1.e4791addcadb8p-4"), start=1)}


def rl_eligible(stack) -> bool:
    """Single-rater items are excluded from RL (lab ITERATION3, R5)."""
    return len(stack) >= 2


def voters(sigs) -> tuple[int, ...]:
    """Raters with a topology vote: a non-empty size-filtered signature; every rater when all abstain (module
    docstring)."""
    return tuple(i for i, s in enumerate(sigs) if s) or tuple(range(len(sigs)))


def _drop_tables(pid: int, handle: int) -> None:
    if os.getpid() == pid:  # a forked child never frees its parent's tables
        _gpucheck.drop_tables(handle)


def device_many(caches) -> None:
    """Builds the device tables of every cache that has none in this process, in one engine call
    (objdel._gpucheck.add_tables_from_raters: one upload, one host wait), and stores each one's handle and returned
    scalars. The reward scores a step's answers after one such call for the step's cold caches (score_many)."""
    pid = os.getpid()
    cold = list({id(c): c for c in caches if c._device is None or c._device[0] != pid}.values())
    if not cold:
        return
    r = _gpucheck.add_tables_from_raters([c.raters.view(np.uint8) for c in cold],
                                         np.array([c.t5 for c in cold], np.float64),
                                         np.array([c.min_area for c in cold], np.int64))
    ss, sig, a = r["sig_start"], r["sig"].tolist(), r["rater_area_start"]
    for i, c in enumerate(cold):
        h = int(r["handles"][i])
        weakref.finalize(c, _drop_tables, pid, h).atexit = False  # the process's exit frees the device
        c._device = (pid, h)
        c._built = (int(r["perimeter"][i]), tuple(tuple(sig[ss[j]:ss[j + 1]]) for j in range(a[i], a[i + 1])))


class StaleCacheError(ValueError):
    """The cache file was built under another reward fingerprint (CACHE_SPEC, size rule or a numpy / scipy / OpenCV
    version): re-run prepare_rl, or load it with rebuild=True."""


def boundary(m: np.ndarray) -> np.ndarray:
    """Foreground pixels with a 4-adjacent background pixel (outside the array counts as background)."""
    m8 = np.ascontiguousarray(m, dtype=np.uint8)
    return m8 > cv2.erode(m8, _CROSS, borderType=cv2.BORDER_CONSTANT, borderValue=0)


def edt(nonzero: np.ndarray) -> np.ndarray:
    """Exact Euclidean distance (float32) from each nonzero pixel to the nearest zero pixel; zero pixels get 0."""
    src = np.ascontiguousarray(nonzero, dtype=np.uint8)
    if src.all():  # no zero pixel: scipy's distances would be meaningless (OpenCV returned its 65536 sentinel)
        return np.full(src.shape, _BIG, np.float32)
    return ndimage.distance_transform_edt(src).astype(np.float32)  # correctly rounded on every CPU (docstring)


def band_depth(mask, union, core, dist_out_union, dist_in_core) -> float:
    """D = sum over error pixels of (d - 1/2): spill (mask outside U, d = distance to U) plus core miss
    (I outside mask, d = distance to the complement of I). float64 sums in a fixed order (deterministic)."""
    ds = dist_out_union[mask & ~union]
    dm = dist_in_core[core & ~mask]
    return (float(np.sum(ds, dtype=np.float64) - 0.5 * ds.size)
            + float(np.sum(dm, dtype=np.float64) - 0.5 * dm.size))


def excursion_moment(m: int) -> float:
    """S_m = E[((X - max Y)+)^2] + E[((min Y - X)+)^2] for X, Y_1..Y_m i.i.d. N(0, 1): S_M for m <= 8, else quad."""
    return S_M[m] if m in S_M else _excursion_quad(m)


@cache
def _excursion_quad(m: int) -> float:
    pdf_max = lambda y: m * stats.norm.pdf(y) * stats.norm.cdf(y) ** (m - 1)
    tail2 = lambda y: (1 + y * y) * stats.norm.sf(y) - y * stats.norm.pdf(y)  # E[((X - y)+)^2]
    return 2.0 * integrate.quad(lambda y: pdf_max(y) * tail2(y), -12, 12, epsabs=1e-12)[0]


def _crop(stack: np.ndarray) -> np.ndarray:
    """The stack on the union's bounding box plus one background pixel on every side. Contour distances and band
    distances are unchanged: every relevant target pixel lies inside the box."""
    rows, cols = np.flatnonzero(stack.any((0, 2))), np.flatnonzero(stack.any((0, 1)))
    if not rows.size:
        return stack
    return np.pad(stack[:, rows[0]:rows[-1] + 1, cols[0]:cols[-1] + 1], ((0, 0), (1, 1), (1, 1)))


def _contour_distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Distance from each contour pixel of a to the nearest contour pixel of b."""
    ba, bb = boundary(a), boundary(b)
    if not bb.any():
        return np.full(np.count_nonzero(ba), _BIG, np.float32)
    return edt(~bb)[ba]


def rater_tau(stack: np.ndarray) -> float:
    """max(tau_floor, P95 of the pooled symmetric contour distances over all rater pairs) in px."""
    st = _crop(stack)
    d = np.concatenate([np.concatenate([_contour_distances(a, b), _contour_distances(b, a)])
                        for a, b in combinations(st, 2)])
    raw = float(np.percentile(d, CACHE_SPEC.tau_quantile)) if d.size else 0.0
    return max(CACHE_SPEC.tau_floor, raw)


def _band(stack: np.ndarray):
    s = stack.sum(0, dtype=np.int16)
    union, core = s > 0, s == len(stack)
    return union, core, max(1, np.count_nonzero(boundary(union))), edt(~union), edt(core)


def f5_tolerance(stack: np.ndarray) -> float:
    """T5 in px (module docstring); needs >= 2 raters."""
    st, m = _crop(stack), len(stack)
    q = []
    for j in range(m):
        union, core, perimeter, d_out, d_in = _band(np.delete(st, j, axis=0))
        q.append(2.0 * band_depth(st[j], union, core, d_out, d_in) / perimeter)
    q_m = float(np.mean(q)) * excursion_moment(m) / excursion_moment(m - 1)
    floor = CACHE_SPEC.tau_floor * np.sqrt(excursion_moment(m) / _X_STAR) / _P95_PAIR_OVER_SIGMA
    return max(float(np.sqrt(q_m / _X_STAR)), float(floor))


def _rater_stack(stack) -> np.ndarray:
    stack = np.ascontiguousarray(stack, dtype=bool)
    if stack.ndim != 3 or not rl_eligible(stack):
        raise ValueError(f"need a (k >= 2, H, W) rater stack, got shape {stack.shape}")
    if len(stack) > 32:  # the engine's tables hold at most 32 raters
        raise ValueError(f"at most 32 raters, got {len(stack)}")
    return stack


def _digest(raters: np.ndarray, tau: float, t5: float, fp: str) -> str:
    h = hashlib.sha256(repr(raters.shape).encode())
    h.update(np.packbits(raters).tobytes())
    h.update(np.array([tau, t5], np.float64).tobytes())
    h.update(fp.encode())
    return h.hexdigest()


class ExpertCache:
    """Rater-only part of reward v4 for one object of one item (module docstring). Rater order must be the
    same in every object's cache of an item (score_completion's joint topology uses it).

    Construction and `load` are host-only (numpy counts, no engine call: a parent must not open a CUDA context before
    it forks). The device tables, and the scalars that come back with them (perimeter, rater signatures, voters), are
    built by the first `device()` (or `device_many`) in the process; the distance maps on the host only when read."""

    def __init__(self, raters: np.ndarray, tau: float, t5: float):
        if not (math.isfinite(t5) and t5 > 0):
            raise ValueError(f"t5 must be finite and > 0, got {t5}")
        raters = _rater_stack(raters)
        self.raters, self.tau, self.t5, self.k = raters, float(tau), float(t5), len(raters)
        s = raters.sum(0, dtype=np.int16)
        self.union, self.core, self.vote = s > 0, s == self.k, 2 * s > self.k
        self.n_vote = int(np.count_nonzero(self.vote))
        self.rater_areas = tuple(int(np.count_nonzero(m)) for m in raters)
        self.min_area = min_area(int(np.count_nonzero(self.union)))
        x, y, w, h = cv2.boundingRect(self.union.view(np.uint8))
        self.box = (y, y + h, x, x + w)  # union bounding box (empty if the union is empty)
        self.digest = _digest(raters, self.tau, self.t5, fingerprint())
        self._device = None  # (pid, handle) of the tables on the checker's device
        self._built = None   # (perimeter, rater_sigs) returned with the tables

    def device(self) -> int:
        """Handle of this cache's tables on the checker's device, which the reward scores against: built there from
        the rater masks at the first call in this process (device_many; a forked child builds its own), resident until
        the cache is garbage-collected."""
        device_many([self])
        return self._device[1]

    @property
    def perimeter(self) -> int:
        """Contour length of the union: max(1, its 4-adjacent boundary pixels) (from the device tables)."""
        self.device()
        return self._built[0]

    @property
    def rater_sigs(self) -> tuple[tuple[int, ...], ...]:
        """Per rater, the size-filtered topology signature at min_area (from the device tables)."""
        self.device()
        return self._built[1]

    @property
    def voters(self) -> tuple[int, ...]:
        return voters(self.rater_sigs)

    @property
    def same_sig(self) -> dict[tuple[int, ...], list[int]]:
        """Signature -> voting raters who drew it (the rater topology histogram)."""
        same = defaultdict(list)
        for i in self.voters:
            same[self.rater_sigs[i]].append(i)
        return dict(same)

    @cached_property
    def dist_out_union(self) -> np.ndarray:
        return edt(~self.union)

    @cached_property
    def dist_in_core(self) -> np.ndarray:
        return edt(self.core)

    def table_args(self) -> tuple:
        """The 13 arguments of objdel._gpucheck.add_tables for this cache (the checker benchmark's replay dumps them,
        experiments/phase2/perf/bench_checker.py)."""
        sigs = self.rater_sigs
        return (self.raters.view(np.uint8), self.union.view(np.uint8), self.core.view(np.uint8),
                self.vote.view(np.uint8), np.where(self.union, self.dist_in_core, self.dist_out_union),
                np.array(self.rater_areas, np.int64), np.isin(np.arange(self.k), self.voters).astype(np.uint8),
                np.cumsum([0, *map(len, sigs)], dtype=np.int64), np.array([v for s in sigs for v in s], np.int32),
                self.n_vote, self.perimeter, self.t5, self.min_area)

    def __reduce__(self):  # pickle only what build computed; the rest is rederived (exactly) on unpickling
        return ExpertCache, (self.raters, self.tau, self.t5)

    @classmethod
    def build(cls, stack) -> "ExpertCache":
        stack = _rater_stack(stack)
        return cls(stack, rater_tau(stack), f5_tolerance(stack))

    def save(self, path) -> str:
        """Compressed .npz with the content's SHA-256; returns the SHA-256."""
        fp = fingerprint()
        with open(path, "wb") as f:
            np.savez_compressed(f, raters=self.raters, tau=self.tau, t5=self.t5, fingerprint=fp, sha256=self.digest)
        return self.digest

    @classmethod
    def load(cls, path, expected_digest: str | None = None, rebuild: bool = False) -> "ExpertCache":
        """The cache saved at `path`. ValueError if the file is corrupted or its digest is not `expected_digest`
        (e.g. the prepared record's); StaleCacheError if it was built under another fingerprint, unless rebuild=True:
        then tau and T5 are rebuilt from the verified masks (the result's digest is the current fingerprint's)."""
        with np.load(Path(path), allow_pickle=False) as z:
            raters, tau, t5 = z["raters"].astype(bool), float(z["tau"]), float(z["t5"])
            fp, digest = str(z["fingerprint"]), str(z["sha256"])
        if _digest(raters, tau, t5, fp) != digest:
            raise ValueError(f"{path}: SHA-256 mismatch (corrupted or edited cache)")
        if expected_digest is not None and digest != expected_digest:
            raise ValueError(f"{path}: cache and expected digests differ ({digest[:12]} vs {expected_digest[:12]})")
        if fp != fingerprint():
            if not rebuild:
                raise StaleCacheError(f"{path}: built under another reward fingerprint ({fp[:12]}, now "
                                      f"{fingerprint()[:12]}: CACHE_SPEC, size rule or numpy/scipy/OpenCV changed); "
                                      "re-run prepare_rl or load with rebuild=True")
            return cls.build(raters)
        return cls(raters, tau, t5)
