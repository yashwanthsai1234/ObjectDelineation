"""Egg-yolk band scores: how far a mask leaves the rater band [core I, union U].

Error pixels are spill (pred outside U) and core miss (I outside pred). Each error pixel x
has a depth d(x) >= 1: the Euclidean distance (pixel centres) to the nearest pixel of U
(spill) or of the complement of I (miss); d - 1/2 is its distance beyond the band edge.
Notation: A = number of error pixels, D = sum over error pixels of (d - 1/2),
P = contour length of U in px, tau = per-image tolerance in px.

For a uniform normal excursion of e px along the whole contour, A = e P and D = e^2 P / 2,
so every formula below except F0 is a function of e / tau only (scale-free):

    F0 (a, current) = 1 - A / |U|                     (scale-biased: area in the denominator)
    F1 (b)          = 1 - A / (tau P)                 = 1 - e / tau
    F2 (c)          = 1 - 2 D / (tau^2 P)             = 1 - (e / tau)^2
    F3              = 1 / (1 + 2 D / (tau^2 P))       = 1 / (1 + (e / tau)^2)
    F4 (d, derived) = 1 / (1 + 2 D / (T4^2 P)),       T4 = tau * g(m)
    F5 (d', derived)= 1 / (1 + 2 D / (T5^2 P)),       T5 self-calibrated from the raters

F0..F2 are clipped to [0, 1]. F3 is F2 made smooth and non-saturating (keeps ranking signal
among bad rollouts). F4 adds a rater-count-aware tolerance. F5 is the recommended formula
(see REPORT.md, K8). Model: at each contour point the
raters' boundary positions are i.i.d. N(0, sigma^2) along the normal. The band of m raters
spans [min, max] of m draws, so a new expert X leaves it by (X - max)+ or (min - X)+, with
second moment S_m * sigma^2 (S_1 = 2, S_2 = 0.78, S_3 = 0.45, S_5 = 0.22, S_6 = 0.17): the
band of few raters is narrow and a new expert legitimately leaves it more. tau = P95 of
pairwise contour distances = 2.772 sigma in the model. Choosing
    tau_m^2 = S_m sigma^2 / x*,   x* = 1 / s* - 1,   s* = 0.95
makes a new expert's expected F4 equal s* for any m and any size; errors beyond expert
variability then decay as 1 / (1 + (excursion / tau_m)^2). Nothing is fitted to data.

F5 drops the Gaussian P95 -> sigma conversion (real rater disagreement is not Gaussian and
differs between datasets). Per image, each rater j is scored against the band of the other
m - 1 raters: q_j = 2 D_j / P_j is that rater's mean squared excursion. The expected value
for a new expert against all m raters is q_m = mean_j(q_j) * S_m / S_(m-1), and
T5^2 = max(q_m / x*, (tau_floor g(m))^2). Needs m >= 2; with one rater T5 = T4.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy import integrate, stats

FORMULAS = ("F0", "F1", "F2", "F3", "F4", "F5")
EXPERT_SCORE = 0.95  # s*: expected F4 of a new expert's mask
P95_PAIR_OVER_SIGMA = 1.959963984540054 * np.sqrt(2.0)  # P95 of |X - Y|, X, Y ~ N(0, 1)


@lru_cache(maxsize=None)
def excursion_second_moment(m: int) -> float:
    """S_m = E[((X - max Y)+)^2] + E[((min Y - X)+)^2], X, Y_1..Y_m i.i.d. N(0, 1)."""
    if m < 1:
        raise ValueError("need at least one rater")
    pdf_max = lambda y: m * stats.norm.pdf(y) * stats.norm.cdf(y) ** (m - 1)  # noqa: E731
    tail2 = lambda y: (1 + y * y) * stats.norm.sf(y) - y * stats.norm.pdf(y)  # E[((X - y)+)^2]  # noqa: E731
    return 2.0 * integrate.quad(lambda y: pdf_max(y) * tail2(y), -12, 12, epsabs=1e-12)[0]


def tolerance_factor(m: int, expert_score: float = EXPERT_SCORE) -> float:
    """g(m) with tau_m = tau * g(m) (see module docstring)."""
    x_star = 1.0 / expert_score - 1.0
    return float(np.sqrt(excursion_second_moment(m) / x_star) / P95_PAIR_OVER_SIGMA)


@dataclass(frozen=True)
class BandErrors:
    n_spill: int
    n_miss: int
    depth_spill: float  # sum of (d - 1/2) over spill pixels
    depth_miss: float  # sum of (d - 1/2) over missed core pixels

    @property
    def area(self) -> int:
        return self.n_spill + self.n_miss

    @property
    def depth(self) -> float:
        return self.depth_spill + self.depth_miss


def band_errors(mask: np.ndarray, union: np.ndarray, core: np.ndarray,
                dist_out_union: np.ndarray, dist_in_core: np.ndarray) -> BandErrors:
    spill = mask & ~union
    miss = core & ~mask
    ds = dist_out_union[spill]
    dm = dist_in_core[miss]
    # float64 accumulation in a fixed order: deterministic
    return BandErrors(int(ds.size), int(dm.size),
                      float(np.sum(ds, dtype=np.float64) - 0.5 * ds.size),
                      float(np.sum(dm, dtype=np.float64) - 0.5 * dm.size))


def band_scores(err: BandErrors, union_area: int, perimeter: int, tau: float,
                tol4: float, tol5: float) -> dict[str, float]:
    """All candidate formulas (tol4 = T4, tol5 = T5 in px, see module docstring)."""
    if union_area == 0:
        ok = float(err.area == 0)
        return {f: ok for f in FORMULAS}
    sq = 2.0 * err.depth / perimeter  # mean squared excursion (px^2)
    x = sq / (tau * tau)
    return {
        "F0": max(0.0, 1.0 - err.area / union_area),
        "F1": max(0.0, 1.0 - err.area / (tau * perimeter)),
        "F2": max(0.0, 1.0 - x),
        "F3": 1.0 / (1.0 + x),
        "F4": 1.0 / (1.0 + sq / (tol4 * tol4)),
        "F5": 1.0 / (1.0 + sq / (tol5 * tol5)),
    }


@dataclass(frozen=True)
class Band:
    """Band geometry of a rater set (the per-prompt part of the band score)."""
    union: np.ndarray
    core: np.ndarray
    n_union: int
    perimeter: int
    dist_out_union: np.ndarray
    dist_in_core: np.ndarray

    @classmethod
    def of(cls, stack: np.ndarray) -> "Band":
        from .distances import boundary, edt

        s = stack.sum(0, dtype=np.int16)
        union, core = s > 0, s == len(stack)
        return cls(union, core, int(np.count_nonzero(union)), max(1, int(np.count_nonzero(boundary(union)))),
                   edt(~union), edt(core))

    def errors(self, mask: np.ndarray) -> BandErrors:
        return band_errors(mask, self.union, self.core, self.dist_out_union, self.dist_in_core)


def self_calibrated_sq(stack: np.ndarray) -> float:
    """Mean over raters j of 2 D_j / P_j, rater j scored against the band of the others (needs >= 2 raters)."""
    from .distances import crop_to_union

    stack = crop_to_union(np.asarray(stack, bool))
    out = []
    for j in range(len(stack)):
        b = Band.of(np.delete(stack, j, axis=0))
        out.append(2.0 * b.errors(stack[j]).depth / b.perimeter)
    return float(np.mean(out))


def expert_sigma(stack: np.ndarray) -> float:
    """Model-equivalent per-point rater spread (px): sqrt(mean_j q_j / S_(m-1)). Needs >= 2 raters."""
    return float(np.sqrt(self_calibrated_sq(stack) / excursion_second_moment(len(stack) - 1)))


def tolerances(stack: np.ndarray, tau: float, tau_floor: float,
               sigma_single: float | None = None) -> tuple[float, float]:
    """(T4, T5) in px for a rater stack of m masks with tolerance tau.

    With one rater, T5 needs an external spread estimate sigma_single (e.g. a size prior);
    without it T5 = T4.
    """
    m = len(stack)
    t4 = tau * tolerance_factor(m)
    x_star = 1.0 / EXPERT_SCORE - 1.0
    if m < 2:
        if sigma_single is None:
            return t4, t4
        t5 = max(sigma_single * float(np.sqrt(excursion_second_moment(1) / x_star)), tau_floor * tolerance_factor(1))
        return t4, t5
    q_m = self_calibrated_sq(stack) * excursion_second_moment(m) / excursion_second_moment(m - 1)
    t5 = max(float(np.sqrt(q_m / x_star)), tau_floor * tolerance_factor(m))
    return t4, t5
