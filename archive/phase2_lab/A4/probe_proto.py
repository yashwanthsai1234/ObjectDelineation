"""A4 sandbox: precision-probe prototype. Calibrates the probe metric (ASSD / HD95 / NSD / IoU) with shape-blind
baselines, on (1) real LIDC val reader masks and (2) synthetic LIDC-sized shapes rendered with our exact rasterizer.

Baselines per object: 'circle' = equal-area disc at the centroid (knows location + size, no shape); 'bbox_ellipse'
= ellipse inscribed in the bounding box (what a box-grounding VLM gives, Phase 0 Q6); 'hull' = convex hull;
'jitter_s' = the true polygon with every vertex moved by N(0, s px) (s = 0.5, 1, 2): the coordinate-noise model.
Run: .venv/bin/python -I probe_proto.py --out probe_proto.json
"""
import argparse
import json

import cv2
import numpy as np
from scipy import ndimage

from objdel.data import store
from objdel.geometry import GRID, canonicalize, to_mask
from objdel.targets import frame

rng = np.random.default_rng(0)


def boundary(m):
    return m & ~ndimage.binary_erosion(m, structure=ndimage.generate_binary_structure(2, 1), border_value=0)


def surf(a, b):
    a, b = np.pad(a, 8), np.pad(b, 8)
    ba, bb = boundary(a), boundary(b)
    if not ba.any() or not bb.any():
        return None
    da, db = ndimage.distance_transform_edt(~bb)[ba], ndimage.distance_transform_edt(~ba)[bb]
    d = np.concatenate([da, db])
    return {"assd": float(d.mean()), "hd95": float(max(np.percentile(da, 95), np.percentile(db, 95))),
            "nsd1": float((d <= 1).mean()), "nsd2": float((d <= 2).mean())}


def iou(a, b):
    u = np.count_nonzero(a | b)
    return 1.0 if u == 0 else np.count_nonzero(a & b) / u


def baselines(m):
    H, W = m.shape
    ys, xs = np.nonzero(m)
    out = {}
    c = np.zeros_like(m, np.uint8)
    r = np.sqrt(len(xs) / np.pi)
    yy, xx = np.mgrid[:H, :W]
    out["circle"] = (xx + .5 - (xs.mean() + .5)) ** 2 + (yy + .5 - (ys.mean() + .5)) ** 2 <= r * r
    x0, x1, y0, y1 = xs.min(), xs.max() + 1, ys.min(), ys.max() + 1
    out["bbox_ellipse"] = (((xx + .5 - (x0 + x1) / 2) / ((x1 - x0) / 2)) ** 2
                           + ((yy + .5 - (y0 + y1) / 2) / ((y1 - y0) / 2)) ** 2) <= 1
    pts = np.c_[xs, ys].astype(np.int32)
    corners = np.concatenate([pts, pts + [1, 0], pts + [0, 1], pts + [1, 1]])
    hull = cv2.convexHull(corners.astype(np.int32))
    cv2.fillPoly(c, [hull], 1)
    out["hull"] = c.astype(bool)
    return out


def jittered(parts, W, H, s):
    """Vertices moved by N(0, s px) on the 0-1000 grid, re-rasterized; illegal results count as empty."""
    k = GRID / W
    new = [(np.rint(o + rng.normal(0, s * k, o.shape)).astype(np.int64),
            [np.rint(h + rng.normal(0, s * k, h.shape)).astype(np.int64) for h in hs]) for o, hs in parts]
    g = to_mask(canonicalize({"x": new})["x"], W, H)
    return g.mask if g.legal else None


def synth_shape(diam, S):
    """Smooth star-shaped polygon (radial Fourier, 3-5 harmonics) of equivalent diameter ~diam px, centred with
    jitter, on an S x S frame; 10 % get a hole (as LIDC cavitary nodules), returned as grid-unit parts."""
    n = 64
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    rr = 1 + sum(rng.uniform(0, 0.25 / k) * np.cos(k * th + rng.uniform(0, 2 * np.pi)) for k in range(2, 6))
    rr *= diam / 2
    cx, cy = S / 2 + rng.uniform(-S / 8, S / 8, 2)
    k = GRID / S
    outer = np.rint(np.c_[cx + rr * np.cos(th), cy + rr * np.sin(th)] * k).astype(np.int64)
    holes = []
    if rng.random() < 0.1 and diam >= 10:
        hr = diam * rng.uniform(0.12, 0.22)
        holes = [np.rint(np.c_[cx + hr * np.cos(th[::4]), cy + hr * np.sin(th[::4])] * k).astype(np.int64)]
    return canonicalize({"x": [(outer, holes)]})["x"]


def score(rows, name, gt, pred):
    s = surf(gt, pred)
    if s is not None:
        rows.setdefault(name, []).append({"iou": iou(gt, pred), **s})


def summ(rows):
    return {k: {m: {"p50": round(float(np.median([r[m] for r in v])), 3),
                    "p90": round(float(np.percentile([r[m] for r in v], 90)), 3)} for m in v[0]} | {"n": len(v)}
            for k, v in rows.items()}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="probe_proto.json")
    ap.add_argument("--every", type=int, default=3)
    a = ap.parse_args()
    real, diams = {}, []
    for i, it in enumerate(store.load("lidc", "val", "data/store")):
        if i % a.every:
            continue
        fr = frame(it)
        st = fr.masks["nodule"]
        for r in range(st.shape[0]):
            m = st[r]
            if m.sum() < 5:
                continue
            diams.append(2 * np.sqrt(m.sum() / np.pi))
            for k, b in baselines(m).items():
                score(real, k, m, b)
            # one other rater as the 'prediction' (human), for reference
            for s2 in range(st.shape[0]):
                if s2 != r and st[s2].sum() >= 5:
                    score(real, "other_reader", m, st[s2])
                    break
    synth, illegal = {}, {}
    for _ in range(600):
        d = float(rng.choice(diams))
        S = 50 if d < 20 else 100
        parts = synth_shape(max(d, 4.0), S)
        g = to_mask(parts, S, S)
        if not g.legal or g.mask.sum() < 5:
            continue
        for k, b in baselines(g.mask).items():
            score(synth, k, g.mask, b)
        for s in (0.5, 1.0, 2.0):
            jm = jittered(parts, S, S, s)
            illegal[s] = illegal.get(s, 0) + (jm is None)
            if jm is not None:
                score(synth, f"jitter_{s}", g.mask, jm)
    out = {"real_lidc_val": summ(real), "synthetic": summ(synth),
           "jitter_illegal_of_600": illegal, "diam_px": {p: round(float(np.percentile(diams, p)), 2) for p in (10, 50, 90)}}
    print(json.dumps(out, indent=1))
    json.dump(out, open(a.out, "w"), indent=1)
