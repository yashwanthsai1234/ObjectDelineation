"""Superseded and candidate reward variants for the A/B tests (EXECUTION.md §8), on packed rollouts.

Not production code: objdel.reward implements reward v4 only. Variants:
- phase0 (AB2): the Phase 0 reward (deleted `objdel.reward.annotators`): IoU vs the strict-majority vote,
  band 1 - (missed core + spill) / |U|, nearest-Betti topology gated on IoU >= 0.5, all on size-filtered
  masks; per-group standardize each term, mix 0.3 / 0.3 / 0.4, standardize again (`group_advantages`).
- v2 (AB2): lab ITERATION2 (t_stack/prod.py): IoU_tc with a max-over-all-raters fallback, F5, legality gate,
  per-term group z with floors, batch-level z of the geometric sum, un-debiased and unclipped DM-loo topology
  term mean-centred with c_topo = 2.
- v4_tauP (AB4): g_area with the tau * P tolerance of lab ITERATION4.
- v4_hole_tc / v4_hole_all (AB3): v4 plus a hole-matching term (weight 0.3, sigma floor 0.05): predicted holes
  matched to consensus holes by Hungarian assignment on IoU (scipy linear_sum_assignment), score = sum of
  matched IoUs / max(#predicted, #consensus), 1 if both are empty. Consensus holes = pixels inside a hole of a
  strict majority of the raters, components >= min_area. 'tc': raters with the prediction's signature
  (topology-conditional, like IoU_tc; all raters if none); 'all': all raters.

References:
- scipy.optimize.linear_sum_assignment (modified Jonker-Volgenant, Crouse 2016):
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html
- Hole matching as instances matched by IoU: docs/research/notes/Phase 1 stack and algorithms/mask_algorithms.md
  Q4; Metrics Reloaded https://arxiv.org/abs/2206.01653
"""
from functools import partial

import numpy as np
from scipy import ndimage
from scipy.optimize import linear_sum_assignment

from objdel.checker.topology import betti, size_filter
from objdel.reward import ExpertCache, score_mask

from .reward_harness import FIELDS, G, advantages, packed

AB_FIELDS = (*FIELDS, "iou_max_all", "p0_iou", "p0_band", "p0_near", "g_tauP", "hole_tc", "hole_all")
_I = {f: i for i, f in enumerate(AB_FIELDS)}
_STATE: dict[str, dict] = {}  # derived state of the current cache, keyed by its digest


def _holes(m: np.ndarray, area: int) -> list[np.ndarray]:
    h = ndimage.binary_fill_holes(m) & ~m
    lab, n = ndimage.label(h)
    sizes = np.bincount(lab.ravel(), minlength=n + 1)
    return [lab == i for i in range(1, n + 1) if sizes[i] >= area]


def _consensus(hole_stack: np.ndarray, area: int) -> list[np.ndarray]:
    if not len(hole_stack):
        return []
    maj = 2 * hole_stack.sum(0) > len(hole_stack)
    lab, n = ndimage.label(maj)
    sizes = np.bincount(lab.ravel(), minlength=n + 1)
    return [lab == i for i in range(1, n + 1) if sizes[i] >= area]


def _state(cache: ExpertCache) -> dict:
    s = _STATE.get(cache.digest)
    if s is None:
        _STATE.clear()  # one item at a time: keep memory bounded
        f = np.stack([size_filter(m, cache.min_area) for m in cache.raters])
        cnt = f.sum(0)
        rh = np.stack([np.any(_holes(m, cache.min_area) or [np.zeros_like(m)], axis=0) for m in cache.raters])
        s = dict(p0_vote=2 * cnt > len(f), p0_core=cnt == len(f), p0_union=cnt > 0, p0_betti=[betti(m) for m in f],
                 rater_holes=rh, cons_all=_consensus(rh, cache.min_area), cons_by_sig={})
        _STATE[cache.digest] = s
    return s


def hole_match(pred_holes: list[np.ndarray], cons: list[np.ndarray]) -> float:
    if not pred_holes and not cons:
        return 1.0
    if not pred_holes or not cons:
        return 0.0
    iou = np.array([[(p & c).sum() / (p | c).sum() for c in cons] for p in pred_holes])
    r, c = linear_sum_assignment(iou, maximize=True)
    return float(iou[r, c].sum() / max(len(pred_holes), len(cons)))


def pack_ab(mask, cache: ExpertCache) -> tuple:
    t = score_mask(mask, cache)
    s = _state(cache)
    m = np.asarray(mask, bool)
    ious = [((m & r).sum() / (m | r).sum()) if (m | r).any() else 1.0 for r in cache.raters]
    iou_max_all = t.iou_tc if t.topo_freq > 0 else float(max(ious))
    mf = size_filter(m, cache.min_area)
    u = (mf | s["p0_vote"]).sum()
    p0_iou = float((mf & s["p0_vote"]).sum() / u) if u else 1.0
    n_u = s["p0_union"].sum()
    err = (s["p0_core"] & ~mf).sum() + (mf & ~s["p0_union"]).sum()
    p0_band = max(0.0, 1 - err / n_u) if n_u else float(err == 0)
    b0, b1 = betti(mf)
    p0_near = max(1 / (1 + abs(b0 - t0) + abs(b1 - t1)) for t0, t1 in s["p0_betti"])
    a = t.area
    slack = cache.tau * cache.perimeter
    lo, hi = min(cache.rater_areas), max(cache.rater_areas)
    g_tauP = 0.0 if a == 0 else min(1.0, (hi + slack) / a) * (min(1.0, a / (lo - slack)) if lo > slack else 1.0)
    ph = _holes(m, cache.min_area)
    same = cache.same_sig.get(t.sig)
    if same:
        key = tuple(same)
        if key not in s["cons_by_sig"]:
            s["cons_by_sig"][key] = _consensus(s["rater_holes"][list(same)], cache.min_area)
        cons_tc = s["cons_by_sig"][key]
    else:
        cons_tc = s["cons_all"]
    return (*packed(t), iou_max_all, p0_iou, p0_band, p0_near, g_tauP, hole_match(ph, cons_tc),
            hole_match(ph, s["cons_all"]))


def _std(x: np.ndarray) -> np.ndarray:
    sd = x.std(1, keepdims=True)
    return (x - x.mean(1, keepdims=True)) / np.where(sd > 0, sd, 1.0)


def adv_phase0(T: np.ndarray) -> np.ndarray:
    iou = T[..., _I["p0_iou"]]
    topo = np.where(iou >= 0.5, T[..., _I["p0_near"]], 0.0)
    return _std(0.3 * _std(iou) + 0.3 * _std(T[..., _I["p0_band"]]) + 0.4 * _std(topo))


def adv_v2(T: np.ndarray) -> np.ndarray:
    iou = T[..., _I["iou_max_all"]]
    X = np.stack([iou, T[..., _I["f5"]], np.ones(T.shape[:2])], -1)
    z = (X - X.mean(1, keepdims=True)) / np.maximum(X.std(1, keepdims=True), (0.05, 0.05, 0.25))
    geo = (z * (0.3, 0.3, 0.1)).sum(-1)
    sd = geo.std()
    geo = (geo - geo.mean()) / (sd if sd > 1e-8 else 1.0)  # batch-level z
    sid = T[..., _I["sig_id"]]
    pi = np.maximum((sid[:, :, None] == sid[:, None, :]).sum(-1) - 1, 0.5) / (G - 1)
    valid = (T[..., _I["topo_any"]] > 0) & (iou >= 0.5)
    r = np.where(valid, 1 + 0.5 * (np.log(np.where(valid, T[..., _I["topo_freq"]], 1.0)) - np.log(pi)), 0.0)
    return geo + 2.0 * (r - r.mean(1, keepdims=True))


def adv_v4(T: np.ndarray) -> np.ndarray:
    return advantages(T[..., :len(FIELDS)])


def adv_v4_tauP(T: np.ndarray) -> np.ndarray:
    X = T[..., :len(FIELDS)].copy()
    X[..., 2] = T[..., _I["g_tauP"]]
    return advantages(X)


def adv_v4_hole(T: np.ndarray, field: str) -> np.ndarray:
    h = T[..., _I[field]]
    z = (h - h.mean(1, keepdims=True)) / np.maximum(h.std(1, keepdims=True), 0.05)
    return advantages(T[..., :len(FIELDS)]) + 0.3 * z


SCHEMES = {"phase0": adv_phase0, "v2": adv_v2, "v4": adv_v4, "v4_tauP": adv_v4_tauP,
           "v4_hole_tc": partial(adv_v4_hole, field="hole_tc"), "v4_hole_all": partial(adv_v4_hole, field="hole_all")}
