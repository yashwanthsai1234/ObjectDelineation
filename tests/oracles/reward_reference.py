"""Slow, independent reference of reward v4 (tests only): cache values, per-object terms, completion terms and
group advantages, written from the definition (lab ITERATION4.md, plan §3.2) without sharing code with
objdel.reward.

Independence: scipy.ndimage distance transforms and morphology instead of OpenCV; the topology oracle
(scipy labeling of each component on its own) instead of Spaghetti; full-frame computations instead of the
production bounding-box crops; S_m by a fine trapezoid rule with math.erf instead of scipy.integrate.quad;
b(f, G) from scipy.stats.binom instead of an exact sum; advantages by explicit loops over groups.
"""
import math
from functools import cache
from itertools import combinations

import numpy as np
from scipy import ndimage, stats

from objdel.geometry import parse, to_mask

from .topology_oracle import signature

S_STAR, TAU_FLOOR = 0.95, 1.0
W = {"iou_tc": 0.3, "f5": 0.3, "legal": 0.1, "g_area": 0.3}
FLOOR = {"iou_tc": 0.05, "f5": 0.05, "legal": 0.25, "g_area": 0.05}
C_TOPO, LAM, EPS, DROP = 2.0, 0.5, 0.1, 1e-2
CROSS = ndimage.generate_binary_structure(2, 1)


def boundary(m):
    m = np.asarray(m, bool)
    return m & ~ndimage.binary_erosion(np.pad(m, 1), CROSS)[1:-1, 1:-1]


def edt(nonzero):
    """Distance to the nearest zero pixel (1e6 if none)."""
    nonzero = np.asarray(nonzero, bool)
    return np.full(nonzero.shape, 1e6) if nonzero.all() else ndimage.distance_transform_edt(nonzero)


def pad_crop(stack, pad=3):
    """Crop to the union box with `pad` zero pixels (distances inside the box do not change)."""
    ys, xs = np.nonzero(stack.any(0))
    return np.pad(stack[:, ys.min():ys.max() + 1, xs.min():xs.max() + 1], ((0, 0), (pad, pad), (pad, pad)))


@cache
def s_moment(m: int) -> float:
    y = np.linspace(-12, 12, 400_001)
    pdf = np.exp(-y * y / 2) / math.sqrt(2 * math.pi)
    cdf = 0.5 * (1 + np.vectorize(math.erf)(y / math.sqrt(2)))
    tail2 = (1 + y * y) * (1 - cdf) - y * pdf
    f = m * pdf * cdf ** (m - 1) * tail2
    return 2.0 * float(np.sum((f[1:] + f[:-1]) / 2) * (y[1] - y[0]))


def depth(mask, union, core):
    d_out, d_in = edt(~union), edt(core)
    spill, miss = mask & ~union, core & ~mask
    return float((d_out[spill] - 0.5).sum() + (d_in[miss] - 0.5).sum())


def cache(stack) -> dict:
    stack = np.asarray(stack, bool)
    k = len(stack)
    st = pad_crop(stack)
    d = []
    for a, b in combinations(st, 2):
        for p, q in ((a, b), (b, a)):
            d += list(edt(~boundary(q))[boundary(p)])
    tau = max(TAU_FLOOR, float(np.percentile(d, 95)) if d else 0.0)
    q = []
    for j in range(k):
        rest = np.delete(st, j, 0)
        u, c = rest.any(0), rest.all(0)
        q.append(2 * depth(st[j], u, c) / max(1, boundary(u).sum()))
    x_star = 1 / S_STAR - 1
    t5 = max(math.sqrt(np.mean(q) * s_moment(k) / s_moment(k - 1) / x_star),
             TAU_FLOOR * math.sqrt(s_moment(k) / x_star) / (1.959963984540054 * math.sqrt(2)))
    union = stack.any(0)
    a = max(5, math.ceil(0.005 * int(union.sum())))
    return dict(raters=stack, k=k, tau=tau, t5=t5, union=union, core=stack.all(0),
                vote=2 * stack.sum(0) > k, perimeter=max(1, int(boundary(union).sum())), min_area=a,
                sigs=[signature(m, a) for m in stack], areas=[int(m.sum()) for m in stack])


def host_tables(raters, t5: float, min_area: int) -> tuple:
    """The 13 objdel._gpucheck.add_tables arguments of one cache, from the definition on the host (the oracle of
    add_tables_from_raters, D7b): union, core, strict-majority vote, the distance map (inside the union the distance
    to the complement of the core, outside it the distance to the union; 1e6 when no such pixel exists), the
    perimeter max(1, |boundary(U)|), rater areas, voters (non-empty signature, all if none) and signatures."""
    st = np.asarray(raters, bool)
    k = len(st)
    union, core = st.any(0), st.all(0)
    vote = 2 * st.sum(0) > k
    dist = np.where(union, edt(core), edt(~union)).astype(np.float32)
    sigs = [signature(m, min_area) for m in st]
    voter = np.array([len(s) > 0 for s in sigs]) if any(sigs) else np.ones(k, bool)
    return (st.view(np.uint8), union.view(np.uint8), core.view(np.uint8), vote.view(np.uint8), dist,
            st.reshape(k, -1).sum(1).astype(np.int64), voter.view(np.uint8),
            np.cumsum([0, *map(len, sigs)], dtype=np.int64), np.array([v for s in sigs for v in s], np.int32),
            int(vote.sum()), max(1, int(boundary(union).sum())), float(t5), int(min_area))


def votes(c) -> list[int]:
    """Raters with a topology vote: non-empty filtered signature (all-sub-filter raters abstain), all if none."""
    v = [i for i in range(c["k"]) if len(c["sigs"][i]) > 0]
    return v if v else list(range(c["k"]))


def iou(a, b):
    u = (a | b).sum()
    return float((a & b).sum() / u) if u else 1.0


def object_terms(mask, c) -> dict:
    m = np.asarray(mask, bool)
    sig = signature(m, c["min_area"])
    same = [i for i in votes(c) if c["sigs"][i] == sig]
    iou_tc = max(iou(m, c["raters"][i]) for i in same) if same else iou(m, c["vote"])
    f5 = 1 / (1 + (2 * depth(m, c["union"], c["core"]) / c["perimeter"]) / c["t5"] ** 2)
    a, lo, hi = int(m.sum()), min(c["areas"]), max(c["areas"])
    g = 0.0 if a == 0 else min(1, hi / a) * (min(1, a / lo) if lo else 1.0)
    return dict(iou_tc=iou_tc, f5=f5, g_area=g, sig=sig, topo_freq=len(same) / len(votes(c)))


def completion_terms(text, truncated, names, caches, W_, H_) -> dict:
    p = parse(text, names, truncated)
    if not p.ok:
        return dict(format_ok=0.0, legal=0.0, iou_tc=0.0, iou_min=0.0, f5=0.0, g_area=0.0, topo_any=0.0, topo_freq=0.0)
    geo = [to_mask(p.objs[n], W_, H_) for n in names]
    ts = [object_terms(g.mask, caches[n]) for g, n in zip(geo, names, strict=True)]
    k = caches[names[0]]["k"]
    voting = [r for r in range(k) if all(r in votes(caches[n]) for n in names)] or list(range(k))
    joint = np.mean([all(caches[n]["sigs"][r] == t["sig"] for n, t in zip(names, ts, strict=True)) for r in voting])
    means = {key: float(np.mean([t[key] for t in ts])) for key in ("iou_tc", "f5", "g_area")}
    return dict(format_ok=1.0, legal=float(all(g.legal for g in geo)), topo_any=float(joint > 0),
                topo_freq=float(joint), iou_min=float(min(t["iou_tc"] for t in ts)), **means)


def loo_bias(f, G):
    n = np.arange(G)
    return math.log(f) - float(np.sum(stats.binom.pmf(n, G - 1, f) * np.log(np.maximum(n, 0.5) / (G - 1))))


def advantages(terms: dict, G: int = 8) -> np.ndarray:
    """Explicit per-group loops over the definition."""
    n = len(terms["legal"])
    out = np.zeros(n)
    for g0 in range(0, n, G):
        idx = range(g0, g0 + G)
        legal = [terms["legal"][i] >= 0.5 for i in idx]
        ans = [terms["sig_id"][i] if legal[j] else "FAIL" for j, i in enumerate(idx)]  # illegal: one shared answer
        a = [0.0] * G
        for key in W:
            x = [terms[key][i] if legal[j] else 0.0 for j, i in enumerate(idx)]
            mu = sum(x) / G
            sd = max(math.sqrt(sum((v - mu) ** 2 for v in x) / G), FLOOR[key])
            a = [a[j] + W[key] * (x[j] - mu) / sd for j in range(G)]
        r = []
        for j, i in enumerate(idx):
            valid = legal[j] and terms["topo_any"][i] > 0 and terms.get("iou_min", terms["iou_tc"])[i] >= 0.5
            if not valid:
                r.append(0.0)
                continue
            n_same = ans.count(ans[j])
            pi = max(n_same - 1, 0.5) / (G - 1)
            f = terms["topo_freq"][i]
            r.append(max(EPS, 1 + LAM * (math.log(f) - loo_bias(f, G) - math.log(pi))))
        mean_r = sum(r) / G
        a = [a[j] + C_TOPO * (r[j] - mean_r) for j in range(G)]
        spread = max(max(terms[key][i] for i in idx) - min(terms[key][i] for i in idx) for key in W)
        if spread < DROP and len(set(ans)) == 1:
            a = [0.0] * G
        out[g0:g0 + G] = a
    return out
