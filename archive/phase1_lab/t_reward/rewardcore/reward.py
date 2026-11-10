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

from .band import (EXPERT_SCORE, Band, band_errors, band_scores, excursion_second_moment, expert_sigma,
                   tolerance_factor, tolerances)
from .distances import boundary, edt, rater_tau
from .topology import bbox, filtered_signature, instance_topology_full

MIN_AREA_PX, MIN_AREA_FRAC = 5, 0.005


def signature_min_area(union_area: int) -> int:
    """THE size filter (census, targets, reward, eval): components and holes below
    max(5 px, ceil(0.5 % of the rater union)) do not count in the topology signature.
    = 5 px on every LIDC item, about 100 px on a RIGA disc (the Phase 0 values)."""
    return max(MIN_AREA_PX, int(np.ceil(MIN_AREA_FRAC * union_area)))


def area_guard(area: int, cache: "ExpertCache", tolerance: str = "zero") -> float:
    """Smooth inflation/deflation guard from ALL raters' areas (in RL every rater is known, so any
    rater-like answer is inside [A_min, A_max] by construction):

        g = exp(-max(0, log(A / A_hi)) - max(0, log(A_lo / A)))  =  min(1, A_hi / A) * min(1, A / A_lo)

    tolerance 'zero': A_hi = A_max, A_lo = A_min (no constants).
    tolerance 'tauP': A_hi = A_max + tau * P_U, A_lo = A_min - tau * P_U (boundary tolerance tau along
    the union contour; no deflation penalty if A_lo <= 0).
    """
    a_min, a_max = min(cache.rater_areas), max(cache.rater_areas)
    slack = cache.tau * cache.perimeter if tolerance == "tauP" else 0.0
    a_hi, a_lo = a_max + slack, a_min - slack
    if area <= 0:
        return 0.0
    g = min(1.0, a_hi / area)
    if a_lo > 0:
        g *= min(1.0, area / a_lo)
    return float(g)


def rl_eligible(stack: np.ndarray) -> bool:
    """Single-rater policy (iteration 3, R5): items with fewer than 2 raters are EXCLUDED from RL.

    With one rater the band has zero width (its tolerance needs a dataset prior; LOO k = 2 fails
    K8 with a pooled prior), and the topology set is a single answer with f = 1, so the set-valued
    reward degenerates. Such items stay usable for SFT and eval (F4/F5 with a per-dataset prior).
    """
    return len(stack) >= 2


def signature_id(sig: tuple) -> int:
    """Deterministic integer id of a signature (CRC32; independent of PYTHONHASHSEED)."""
    import zlib

    return zlib.crc32(repr(tuple(int(v) for v in sig)).encode())


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
    min_area: int = 1  # signature size filter (signature_min_area of the union)
    class_bands: dict = field(default_factory=dict)  # signature -> (Band, T) for minority topology classes

    @classmethod
    def build(cls, stack: np.ndarray, tau_floor: float = 1.0,
              tau_prior: Callable[[float], float] | None = None,
              sigma_prior: Callable[[float], float] | None = None,
              topology_conditional_band: bool = False) -> "ExpertCache":
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
        min_area = signature_min_area(int(np.count_nonzero(union)))
        topos = [instance_topology_full(m) for m in stack]
        sigs = [filtered_signature(m, min_area, t) for m, t in zip(stack, topos)]
        return cls(
            k=k, vote=vote, union=union, core=core,
            n_vote=int(np.count_nonzero(vote)), n_union=int(np.count_nonzero(union)),
            perimeter=max(1, int(np.count_nonzero(boundary(union)))),
            dist_out_union=edt(~union), dist_in_core=edt(core),
            tau=float(tau), tau_raw=raw, tol4=tol4, tol5=tol5,
            signatures=Counter(sigs),
            bettis=tuple(t.betti for t in topos),
            union_box=bbox(union),
            raters=stack,
            rater_sigs=tuple(sigs),
            rater_areas=tuple(int(np.count_nonzero(m)) for m in stack),
            min_area=min_area,
            class_bands=_class_bands(stack, sigs, tau_floor) if topology_conditional_band else {},
        )


def _class_bands(stack: np.ndarray, sigs: list, tau_floor: float) -> dict:
    """Topology-conditional band (iteration 3): for every topology class drawn by a strict subset
    of the raters, the band of that class's raters and its tolerance T = sigma * sqrt(S_mc / x*),
    sigma = the item's model-equivalent rater spread (all raters, F5 self-calibration), mc = class
    size. A rollout with that topology is band-scored against its own class (as IoU_tc is), so
    consensus geometry + a minority topology ("Frankenstein" output) gets no free band credit.
    For a class equal to all raters this is exactly F5, so unanimous items are unchanged.
    NOT used in the frozen reward (iteration 3): it barely helps the hole adversaries (44 -> 54 %)
    and costs topology-class neutrality on 251 items (+0.038 vs +0.009); kept as an option."""
    if len(set(sigs)) < 2:
        return {}
    sigma = expert_sigma(stack)
    x_star = 1.0 / EXPERT_SCORE - 1.0
    out = {}
    for sg in set(sigs):
        members = stack[[i for i, s in enumerate(sigs) if s == sg]]
        mc = len(members)
        t = max(sigma * float(np.sqrt(excursion_second_moment(mc) / x_star)), tau_floor * tolerance_factor(mc))
        out[sg] = (Band.of(members), t)
    return out


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
    iou_tc: float = 0.0  # topology-conditional IoU: max IoU over raters with the same topology (else IoU vs vote)
    sig: tuple = ()  # size-filtered per-instance topology signature (the one that is scored)
    sig_raw: tuple = ()  # unfiltered signature (logging: sub-filter specks / pinholes)
    g_area: float = 1.0  # area guard, tolerance 'zero' (the frozen v4 form)
    g_area_tauP: float = 1.0  # area guard, tolerance tau * P_U (reported alternative)
    iou_tc_mean: float = 0.0  # mean IoU over the raters with the same topology (else IoU vs vote)


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
    topo = instance_topology_full(m)
    sig = filtered_signature(m, cache.min_area, topo)
    if sig in cache.class_bands:  # minority topology class: band of that class's raters
        cb, t = cache.class_bands[sig]
        full = np.zeros(cache.union.shape, bool)
        full[roi] = m
        err = cb.errors(full)
        band = band_scores(err, cb.n_union, cb.perimeter, cache.tau, t, t)
    else:
        err = band_errors(m, cache.union[roi], cache.core[roi], cache.dist_out_union[roi], cache.dist_in_core[roi])
        band = band_scores(err, cache.n_union, cache.perimeter, cache.tau, cache.tol4, cache.tol5)
    b0, b1 = topo.betti
    hits = cache.signatures.get(sig, 0)
    near = max(1.0 / (1 + abs(b0 - t0) + abs(b1 - t1)) for t0, t1 in cache.bettis)
    return Terms(iou, band, float(hits > 0), hits / cache.k, near, b0, b1, err.n_spill, err.n_miss, err.depth, area,
                 topo_conditional_iou(m, roi, area, sig, cache, iou), sig, topo.signature(),
                 area_guard(area, cache, "zero"), area_guard(area, cache, "tauP"),
                 topo_conditional_iou(m, roi, area, sig, cache, iou, reduce="mean"))


def topo_conditional_iou(m: np.ndarray, roi, area: int, sig: tuple, cache: ExpertCache, iou_vote: float,
                         reduce: str = "max") -> float:
    """Max IoU over the raters whose per-instance topology equals the mask's; IoU vs the vote if no
    rater drew that topology (iteration 3: a wrong topology must not keep a near-1 IoU by matching
    the rater it was derived from).

    Topology-neutral: a sample of a minority topology is compared with the raters who drew that
    topology, not with a vote mask that lacks it (iteration 2, results/iter2/topo_gap.txt).
    """
    same = [i for i, s in enumerate(cache.rater_sigs) if s == sig]
    if not same:
        return iou_vote
    vals = []
    for i in same:
        inter = int(np.count_nonzero(m & cache.raters[i][roi]))
        uni = area + cache.rater_areas[i] - inter
        vals.append(inter / uni if uni else 1.0)
    # reduce='max': closest rater of that topology (iteration 2/3); 'mean': expected IoU vs a random
    # rater of that topology (iteration 4: hedges far from the opposite-extreme raters score lower)
    return float(max(vals) if reduce == "max" else np.mean(vals))


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
