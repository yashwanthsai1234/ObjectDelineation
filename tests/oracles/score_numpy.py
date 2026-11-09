"""The host scorer that the CUDA checker's scorer replaced (tests only): objdel.reward.score.score_mask and
band_depth as they were before scoring moved onto the device (git ed74aca), unchanged, so the tests can demand
bit-identical results from the device on any mask. Masks and the per-completion combination come from the caller.

Reward v5 (B3.2) adds the centroid sums, from np.nonzero of the full-frame mask, and the location terms, with the
same float64 arithmetic as objdel.reward.score (location, gated_loc).

Numpy and OpenCV: bounding-box crop, Spaghetti-labeled size-filtered signature (objdel.checker.topology), float64
sums of the float32 band distances.
"""
import cv2
import numpy as np

from objdel.checker.topology import signature


def iou(a, area_a, b, area_b) -> float:
    inter = int(np.count_nonzero(a & b))
    union = area_a + area_b - inter
    return inter / union if union else 1.0


def area_guard(area: int, cache) -> float:
    if area <= 0:
        return 0.0
    a_min, a_max = min(cache.rater_areas), max(cache.rater_areas)
    return min(1.0, a_max / area) * (min(1.0, area / a_min) if a_min > 0 else 1.0)


def band_depth(mask, union, core, dist_out_union, dist_in_core) -> float:
    ds = dist_out_union[mask & ~union]
    dm = dist_in_core[core & ~mask]
    return (float(np.sum(ds, dtype=np.float64) - 0.5 * ds.size)
            + float(np.sum(dm, dtype=np.float64) - 0.5 * dm.size))


def loc_raw(sum_x, sum_y, area, cache, width=1.0) -> float:
    ref = cache.vote if cache.n_vote > 0 else cache.union
    ys, xs = np.nonzero(ref)
    if area == 0 or xs.size == 0:
        return 0.0
    d2 = (sum_x / area - int(xs.sum()) / xs.size) ** 2 + (sum_y / area - int(ys.sum()) / xs.size) ** 2
    return float(np.exp(-d2 / (2 * width**2 * (int(np.count_nonzero(cache.union)) / np.pi))))


def gated_loc(t: dict, miss_iou=0.1) -> float:
    if t["iou_tc"] >= miss_iou or t["topo_freq"] == 0:
        return 0.0
    return (1 - t["iou_tc"] / miss_iou) * t["g_area"] * t["loc_raw"]


def score_mask(mask, cache) -> dict:
    """{iou_tc, f5, g_area, sig, topo_freq, area, sum_x, sum_y, loc_raw} of one object's mask against its
    ExpertCache."""
    m = np.ascontiguousarray(mask, dtype=bool)
    ys, xs = np.nonzero(m)
    sx, sy = int(xs.sum()), int(ys.sum())
    x, y, w, h = cv2.boundingRect(m.view(np.uint8))
    boxes = [b for b in (cache.box, (y, y + h, x, x + w)) if b[1] > b[0]]
    roi = (slice(min(b[0] for b in boxes), max(b[1] for b in boxes)),
           slice(min(b[2] for b in boxes), max(b[3] for b in boxes))) if boxes else (slice(0, 0), slice(0, 0))
    m = m[roi]
    area = int(np.count_nonzero(m))
    sig = signature(m, cache.min_area)
    same = cache.same_sig.get(sig, ())
    if same:
        iou_tc = max(iou(m, area, cache.raters[i][roi], cache.rater_areas[i]) for i in same)
    else:
        iou_tc = iou(m, area, cache.vote[roi], cache.n_vote)
    d = band_depth(m, cache.union[roi], cache.core[roi], cache.dist_out_union[roi], cache.dist_in_core[roi])
    f5 = 1.0 / (1.0 + (2.0 * d / cache.perimeter) / (cache.t5 * cache.t5))
    return {"iou_tc": iou_tc, "f5": f5, "g_area": area_guard(area, cache), "sig": sig,
            "topo_freq": len(same) / len(cache.voters), "area": area, "sum_x": sx, "sum_y": sy,
            "loc_raw": loc_raw(sx, sy, area, cache)}


def score_completion(text, truncated, names, caches, W, H) -> dict:
    """score_completion's terms, as the host scorer combined them: masks from objdel.geometry.to_masks (the checker's
    raster), objects scored by score_mask above."""
    from objdel.geometry import parse, to_masks
    from objdel.reward import NUMERIC, rater_answers, sig_id

    p = parse(text, names, truncated)
    if not p.ok:
        return dict.fromkeys(NUMERIC, 0.0) | {"sig_id": -1.0, "error": p.error, "defects": {}}
    geoms = to_masks([(p.objs[n], W, H) for n in names])
    a = rater_answers(caches, names)
    legal, defects, ts = True, {}, []
    for n, g in zip(names, geoms, strict=True):
        legal &= g.legal
        for d, c in g.defects.items():
            defects[d] = defects.get(d, 0) + c
        ts.append(score_mask(g.mask, caches[n]))
    sigs = tuple(t["sig"] for t in ts)
    joint = len(a.get(sigs, ())) / sum(map(len, a.values()))
    return {"format_ok": 1.0, "legal": float(legal), "iou_tc": float(np.mean([t["iou_tc"] for t in ts])),
            "iou_min": float(min(t["iou_tc"] for t in ts)), "f5": float(np.mean([t["f5"] for t in ts])),
            "g_area": float(np.mean([t["g_area"] for t in ts])), "topo_any": float(joint > 0), "topo_freq": joint,
            "sig_id": sig_id(sigs), "loc": float(np.mean([gated_loc(t) for t in ts])), "error": "",
            "defects": defects}
