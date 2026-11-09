"""A4 sandbox: the human ceiling on our val frames, for the Phase 2 gates and the precision probe.

Per item (frame resolution, the reward's frame): every pair of voting raters (non-empty masks): IoU, ASSD, HD95,
NSD at 1 and 2 px (pixel-boundary approximation, as MONAI / DeepMind surface-distance on a unit grid); per rater the
leave-one-out IoU_tc analogue (best IoU to another rater with the same size-filtered signature, else IoU to the
others' strict-majority vote) and leave-one-out D_max (best IoU to any other rater). Also object sizes (equivalent
diameter) and, per size, the IoU a uniform 1 / 2 px boundary error costs a disc of that size.
Run: .venv/bin/python -I rater_ceiling.py --out rater_ceiling.json
"""
import argparse
import json
from itertools import combinations

import numpy as np
from scipy import ndimage

from objdel.checker import signature
from objdel.data import store
from objdel.targets import frame, min_areas


def boundary(m):
    return m & ~ndimage.binary_erosion(m, structure=ndimage.generate_binary_structure(2, 1), border_value=0)


def surf(a, b):
    """Directed boundary distances a->b and b->a (px), on a padded canvas."""
    pad = 8
    a, b = np.pad(a, pad), np.pad(b, pad)
    ba, bb = boundary(a), boundary(b)
    da = ndimage.distance_transform_edt(~bb)[ba]
    db = ndimage.distance_transform_edt(~ba)[bb]
    return da, db


def iou(a, b):
    u = np.count_nonzero(a | b)
    return 1.0 if u == 0 else np.count_nonzero(a & b) / u


def run(dataset, every):
    pairs, loo, sizes = [], [], []
    for i, it in enumerate(store.load(dataset, "val", "data/store")):
        if i % every:
            continue
        fr = frame(it)
        ma = min_areas(fr.masks)
        for n, st in fr.masks.items():
            k = st.shape[0]
            voters = [r for r in range(k) if signature(st[r], ma[n]) != ()]  # whole mask below the filter abstains
            for r in voters:
                sizes.append((n, 2 * np.sqrt(st[r].sum() / np.pi)))
            for r, s in combinations(voters, 2):
                da, db = surf(st[r], st[s])
                d = np.concatenate([da, db])
                pairs.append({"obj": n, "iou": iou(st[r], st[s]), "assd": float(d.mean()),
                              "hd95": float(max(np.percentile(da, 95), np.percentile(db, 95))),
                              "nsd1": float((d <= 1).mean()), "nsd2": float((d <= 2).mean()),
                              "diam": float(2 * np.sqrt((st[r].sum() + st[s].sum()) / 2 / np.pi))})
            if len(voters) >= 2:
                for r in voters:
                    others = [s for s in voters if s != r]
                    sr = signature(st[r], ma[n])
                    same = [s for s in others if signature(st[s], ma[n]) == sr]
                    vote = st[others].sum(0) * 2 > len(others)
                    loo.append({"obj": n, "iou_tc": max(iou(st[r], st[s]) for s in same) if same else iou(st[r], vote),
                                "d_max": max(iou(st[r], st[s]) for s in others), "iou_vote": iou(st[r], vote),
                                "topo_match": bool(same), "k": len(voters)})
    return pairs, loo, sizes


def q(x, ps=(10, 25, 50, 75, 90)):
    return {f"p{p}": round(float(np.percentile(x, p)), 3) for p in ps}


def disc_iou(dia, e):
    r = dia / 2
    return max(r - e, 0) ** 2 / r ** 2  # a disc shrunk by e px (inflated: r^2 / (r + e)^2, similar)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="rater_ceiling.json")
    a = ap.parse_args()
    out = {}
    for ds, every in (("riga", 1), ("lidc", 1)):
        pairs, loo, sizes = run(ds, every)
        res = {}
        for n in sorted({p["obj"] for p in pairs}):
            P = [p for p in pairs if p["obj"] == n]
            L = [x for x in loo if x["obj"] == n]
            D = np.array([s for o, s in sizes if o == n])
            res[n] = {"pairs": len(P), **{m: q([p[m] for p in P]) for m in ("iou", "assd", "hd95", "nsd1", "nsd2")},
                      "loo_raters": len(L), "loo_iou_tc": q([x["iou_tc"] for x in L]),
                      "loo_d_max": q([x["d_max"] for x in L]), "loo_iou_vote": q([x["iou_vote"] for x in L]),
                      "loo_topo_match": float(np.mean([x["topo_match"] for x in L])),
                      "loo_iou_tc_ge_0.5": float(np.mean([x["iou_tc"] >= 0.5 for x in L])),
                      "diam_px": q(D),
                      "disc_iou_at_err": {f"{e}px": {f"d{int(d)}": round(disc_iou(d, e), 3)
                                                     for d in np.percentile(D, (10, 50, 90))} for e in (1, 2)}}
            # by size bins (LIDC)
            if ds == "lidc":
                bins = [(0, 6), (6, 10), (10, 15), (15, 25), (25, 1e9)]
                res[n]["by_diam"] = {f"{lo}-{hi}": {"pairs": len(sub), "iou_med": round(float(np.median([p["iou"] for p in sub])), 3),
                                                    "assd_med": round(float(np.median([p["assd"] for p in sub])), 3),
                                                    "hd95_med": round(float(np.median([p["hd95"] for p in sub])), 3),
                                                    "nsd2_med": round(float(np.median([p["nsd2"] for p in sub])), 3)}
                                     for lo, hi in bins if (sub := [p for p in P if lo <= p["diam"] < hi])}
        out[ds] = res
        print(ds, json.dumps(res, indent=1), flush=True)
    json.dump(out, open(a.out, "w"), indent=1)
