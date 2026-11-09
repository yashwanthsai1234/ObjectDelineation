"""A4 sandbox: what per-vertex coordinate error, translation and scale error cost on real LIDC / RIGA targets.

For val targets (certified polygons, the SFT answers): (a) every vertex + N(0, s px), (b) whole polygon shifted by
d px in a random direction, (c) scaled about its centroid by (1 + e). Reports the legal share and, among legal
outputs, IoU and ASSD to the rater's mask. Calibrates the precision probe (which error budget a model may have).
Run: .venv/bin/python -I noise_model.py --out noise_model.json
"""
import argparse
import json

import numpy as np
from scipy import ndimage

from objdel.data import store
from objdel.geometry import GRID, canonicalize, to_mask
from objdel.targets import build_target, frame, min_areas

rng = np.random.default_rng(1)


def boundary(m):
    return m & ~ndimage.binary_erosion(m, structure=ndimage.generate_binary_structure(2, 1), border_value=0)


def assd(a, b):
    a, b = np.pad(a, 8), np.pad(b, 8)
    ba, bb = boundary(a), boundary(b)
    if not ba.any() or not bb.any():
        return None
    return float(np.concatenate([ndimage.distance_transform_edt(~bb)[ba], ndimage.distance_transform_edt(~ba)[bb]]).mean())


def iou(a, b):
    u = np.count_nonzero(a | b)
    return 1.0 if u == 0 else np.count_nonzero(a & b) / u


def apply(parts, f):
    return canonicalize({"x": [(f(o), [f(h) for h in hs]) for o, hs in parts]})["x"]


def run(ds, every, obj, max_items):
    res, nv = {}, []
    seen = 0
    for i, it in enumerate(store.load(ds, "val", "data/store")):
        if i % every:
            continue
        seen += 1
        if seen > max_items:
            break
        fr = frame(it)
        names = list(fr.masks)
        ma = min_areas(fr.masks)
        W = fr.W
        k = GRID / W
        for r in range(fr.masks[names[0]].shape[0]):
            masks = {n: fr.masks[n][r] for n in names}
            if not masks[obj].any():
                continue
            try:
                t = build_target(masks, names, min_area=ma)
            except ValueError:
                continue
            parts, m = t.objs[obj], masks[obj]
            nv.append(t.info[obj]["vertices"])
            allv = np.concatenate([np.concatenate([o, *hs]) for o, hs in parts]).astype(float)
            c = allv.mean(0)
            moves = {f"vertex_{s}": (lambda s: lambda v: np.rint(v + rng.normal(0, s * k, v.shape)).astype(np.int64))(s)
                     for s in (0.25, 0.5, 1.0, 2.0)}
            for d in (0.5, 1.0, 2.0):
                ang = rng.uniform(0, 2 * np.pi)
                moves[f"shift_{d}"] = (lambda dx, dy: lambda v: np.rint(v + [dx, dy]).astype(np.int64))(
                    d * k * np.cos(ang), d * k * np.sin(ang))
            for e in (-0.2, -0.1, 0.1, 0.2):
                moves[f"scale_{e}"] = (lambda e: lambda v: np.rint(c + (v - c) * (1 + e)).astype(np.int64))(e)
            for name, f in moves.items():
                g = to_mask(apply(parts, f), W, fr.H)
                row = res.setdefault(name, {"legal": [], "iou": [], "assd": []})
                row["legal"].append(g.legal)
                if g.legal:
                    row["iou"].append(iou(g.mask, m))
                    a = assd(g.mask, m)
                    if a is not None:
                        row["assd"].append(a)
    out = {name: {"legal": round(float(np.mean(v["legal"])), 3),
                  "iou_p50": round(float(np.median(v["iou"])), 3) if v["iou"] else None,
                  "iou_p10": round(float(np.percentile(v["iou"], 10)), 3) if v["iou"] else None,
                  "assd_p50": round(float(np.median(v["assd"])), 3) if v["assd"] else None, "n": len(v["legal"])}
           for name, v in res.items()}
    return out, {"vertices_p50": float(np.median(nv)), "vertices_p90": float(np.percentile(nv, 90))}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="noise_model.json")
    a = ap.parse_args()
    out = {}
    for ds, every, obj, mx in (("lidc", 3, "nodule", 400), ("riga", 2, "cup", 32)):
        out[ds] = run(ds, every, obj, mx)
        print(ds, json.dumps(out[ds], indent=1), flush=True)
    json.dump(out, open(a.out, "w"), indent=1)
