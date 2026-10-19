"""Reward core: per-prompt expert cache + per-rollout terms + scalar reward.

Everything that depends only on the rater masks is computed once per prompt (ExpertCache)
and reused by all rollouts of the GRPO group (and across steps): vote/union/core masks,
contour length of the union, the two distance transforms used by the band, tau and the
band tolerances, and the set of rater topologies. It never changes for a prompt, so it can be
built offline for the whole dataset. Per rollout the work is a few boolean passes, two gathers
and two CCL passes, all restricted to the bounding box of (prediction | rater union).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from .band import band_errors, band_scores, tolerances
from .distances import boundary, edt, rater_tau
from .topology import bbox, instance_topology_full


@dataclass
class ExpertCache:
    k: int
    vote: np.ndarray  # strict majority
    union: np.ndarray
    core: np.ndarray
    n_vote: int
    n_union: int
    perimeter: int  # contour length of the union (px)
    dist_out_union: np.ndarray  # float32, distance to the union (0 inside)
    dist_in_core: np.ndarray  # float32, distance to the complement of the core (0 outside)
    tau: float
    tau_raw: float
    tol4: float  # F4 tolerance (px)
    tol5: float  # F5 tolerance (px)
    signatures: Counter = field(default_factory=Counter)  # per-instance topology -> #raters
    bettis: tuple = ()
    union_box: tuple | None = None  # (y0, y1, x0, x1) of the union
    raters: np.ndarray | None = None  # (k, H, W) bool rater masks (for the topology-conditional IoU)
    rater_sigs: tuple = ()  # per-instance topology signature of each rater
    rater_areas: tuple = ()

    @classmethod
    def build(cls, stack: np.ndarray, tau_floor: float = 1.0,
              tau_prior: Callable[[float], float] | None = None,
              sigma_prior: Callable[[float], float] | None = None) -> "ExpertCache":
        """tau = max(tau_floor, P95 of pooled symmetric rater-pair contour distances).

        With a single rater there are no pairs: tau = max(tau_floor, tau_prior(r_eq)) and the F5
        spread = sigma_prior(r_eq) if priors are given (r_eq = sqrt(rater area / pi)), else floors.
        """
        stack = np.asarray(stack, bool)
        k = len(stack)
        s = stack.sum(0, dtype=np.int16)
        union, core, vote = s > 0, s == k, 2 * s > k
        tau, raw = rater_tau(stack, tau_floor)
        sigma_single = None
        if k < 2:
            r_eq = float(np.sqrt(stack.sum() / k / np.pi))
            if tau_prior is not None:
                tau = max(tau_floor, float(tau_prior(r_eq)))
            if sigma_prior is not None:
                sigma_single = float(sigma_prior(r_eq))
        tol4, tol5 = tolerances(stack, tau, tau_floor, sigma_single)
        topos = [instance_topology_full(m) for m in stack]
        return cls(
            k=k, vote=vote, union=union, core=core,
            n_vote=int(np.count_nonzero(vote)), n_union=int(np.count_nonzero(union)),
            perimeter=max(1, int(np.count_nonzero(boundary(union)))),
            dist_out_union=edt(~union), dist_in_core=edt(core),
            tau=float(tau), tau_raw=raw, tol4=tol4, tol5=tol5,
            signatures=Counter(t.signature() for t in topos),
            bettis=tuple(t.betti for t in topos),
            union_box=bbox(union),
            raters=stack,
            rater_sigs=tuple(t.signature() for t in topos),
            rater_areas=tuple(int(np.count_nonzero(m)) for m in stack),
        )


@dataclass(frozen=True)
class Terms:
    iou: float  # vs strict-majority mask
    band: dict  # formula name -> score
    topo_any: float  # 1 if the per-instance topology equals some rater's
    topo_freq: float  # share of raters with this per-instance topology
    topo_near: float  # max over raters of 1 / (1 + |d b0| + |d b1|)
    b0: int
    b1: int
    n_spill: int
    n_miss: int
    depth: float  # sum of (d - 1/2) over error pixels
    area: int
    iou_tc: float = 0.0  # topology-conditional IoU: max IoU over raters with the same topology (else all raters)
    sig: tuple = ()  # per-instance topology signature of the mask


def _roi(mask: np.ndarray, cache: ExpertCache) -> tuple[slice, slice]:
    """Bounding box of pred | union: every term depends only on pixels inside it."""
    boxes = [b for b in (bbox(mask), cache.union_box) if b is not None]
    if not boxes:
        return slice(0, 0), slice(0, 0)
    return (slice(min(b[0] for b in boxes), max(b[1] for b in boxes)),
            slice(min(b[2] for b in boxes), max(b[3] for b in boxes)))


def score(mask: np.ndarray, cache: ExpertCache) -> Terms:
    m = np.asarray(mask, bool)
    if m.shape != cache.union.shape:
        raise ValueError(f"mask shape {m.shape} != rater shape {cache.union.shape}")
    roi = _roi(m, cache)
    m = m[roi]
    inter = int(np.count_nonzero(m & cache.vote[roi]))
    area = int(np.count_nonzero(m))
    uni = area + cache.n_vote - inter
    iou = inter / uni if uni else 1.0
    err = band_errors(m, cache.union[roi], cache.core[roi], cache.dist_out_union[roi], cache.dist_in_core[roi])
    band = band_scores(err, cache.n_union, cache.perimeter, cache.tau, cache.tol4, cache.tol5)
    topo = instance_topology_full(m)
    sig = topo.signature()
    b0, b1 = topo.betti
    hits = cache.signatures.get(sig, 0)
    near = max(1.0 / (1 + abs(b0 - t0) + abs(b1 - t1)) for t0, t1 in cache.bettis)
    return Terms(iou, band, float(hits > 0), hits / cache.k, near, b0, b1, err.n_spill, err.n_miss, err.depth, area,
                 topo_conditional_iou(m, roi, area, sig, cache), sig)


def topo_conditional_iou(m: np.ndarray, roi, area: int, sig: tuple, cache: ExpertCache) -> float:
    """Max IoU over the raters whose per-instance topology equals the mask's (all raters if none).

    Topology-neutral: a sample of a minority topology is compared with the raters who drew that
    topology, not with a vote mask that lacks it (iteration 2, results/iter2/topo_gap.txt).
    """
    same = [i for i, s in enumerate(cache.rater_sigs) if s == sig] or range(cache.k)
    best = 0.0
    for i in same:
        inter = int(np.count_nonzero(m & cache.raters[i][roi]))
        uni = area + cache.rater_areas[i] - inter
        best = max(best, inter / uni if uni else 1.0)
    return best


@dataclass(frozen=True)
class RewardConfig:
    w_iou: float = 0.3
    w_band: float = 0.3
    w_topo: float = 0.4
    band_formula: str = "F4"
    topo_term: str = "any"  # 'any' | 'freq' | 'near'
    topo_gate_iou: float = 0.5  # topology credit only when IoU >= gate (no credit for a far-off blob)


DEFAULT_CONFIG = RewardConfig()


def reward(t: Terms, cfg: RewardConfig = DEFAULT_CONFIG) -> float:
    topo = getattr(t, "topo_" + cfg.topo_term) if t.iou >= cfg.topo_gate_iou else 0.0
    return cfg.w_iou * t.iou + cfg.w_band * t.band[cfg.band_formula] + cfg.w_topo * topo


def term_vector(t: Terms, band_formula: str = "F4") -> dict[str, float]:
    return {"iou": t.iou, "band": t.band[band_formula], "topo": t.topo_any}
