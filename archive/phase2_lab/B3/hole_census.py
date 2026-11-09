"""Empirical hole and part geometry of the real non-simple LIDC targets (P5.1; judge MUST-FIX 1).

    .venv/bin/python -I -m experiments.phase2.lab.B3.hole_census --census <scratch>/census.jsonl \
        --out <scratch>/hole_records.jsonl --summary <lab>/hole_census.json

For every non-simple item (census.py) and every voting rater: the rater mask on the model frame after the reward's
size rule (checker.size_filter at min_area(|union|)), then per hole (4-connected background component not touching
the frame border) and per part (8-connected component):
- hole: area px; fraction = hole area / filled area of its part; offset = |hole centroid - filled part centroid| /
  equivalent radius of the filled part; wall = min Euclidean distance (px centres) from a hole pixel to a pixel
  outside the filled part, minus 1 (0 = a one-pixel-thin wall); the part's filled equivalent diameter;
- parts: count, area ratio (2nd largest / largest), gap = min centre distance between the two largest parts minus 1
  (background pixels between them), the object's equivalent diameter.
Reads untrusted data: run with python -I.
"""
import argparse
import json

import numpy as np
from scipy import ndimage

from objdel.checker import size_filter
from objdel.data import store
from objdel.targets import frame

C4 = ndimage.generate_binary_structure(2, 1)
C8 = ndimage.generate_binary_structure(2, 2)


def holes_of(m):
    """Hole label image (4-connected background components not touching the border) and count."""
    bg, n = ndimage.label(~m, structure=C4)
    border = set(np.unique(np.r_[bg[0], bg[-1], bg[:, 0], bg[:, -1]]).tolist())
    keep = [i for i in range(1, n + 1) if i not in border]
    out = np.zeros_like(bg)
    for k, i in enumerate(keep, 1):
        out[bg == i] = k
    return out, len(keep)


def geometry(m):
    parts, n_parts = ndimage.label(m, structure=C8)
    hl, n_holes = holes_of(m)
    rec = {"area": int(m.sum()), "diam": float(2 * np.sqrt(m.sum() / np.pi)), "n_parts": n_parts,
           "n_holes": n_holes, "holes": []}
    for h in range(1, n_holes + 1):
        hm = hl == h
        # the part enclosing the hole: the part label most frequent on the hole's 4-neighbour ring
        ring = ndimage.binary_dilation(hm, structure=C4) & ~hm
        pid = int(np.bincount(parts[ring][parts[ring] > 0]).argmax())
        filled = ndimage.binary_fill_holes(parts == pid)
        fy, fx = np.nonzero(filled)
        hy, hx = np.nonzero(hm)
        r_eq = np.sqrt(filled.sum() / np.pi)
        dist_out = ndimage.distance_transform_edt(np.pad(filled, 1))[1:-1, 1:-1]
        rec["holes"].append({
            "area": int(hm.sum()), "fraction": float(hm.sum() / filled.sum()),
            "offset": float(np.hypot(hx.mean() - fx.mean(), hy.mean() - fy.mean()) / r_eq),
            "wall": float(dist_out[hm].min() - 1), "part_diam": float(2 * r_eq)})
    if n_parts >= 2:
        sizes = np.bincount(parts.ravel())[1:]
        a, b = np.argsort(sizes)[::-1][:2] + 1
        gap = ndimage.distance_transform_edt(parts != a)[parts == b].min() - 1
        rec.update(part_ratio=float(sizes[b - 1] / sizes[a - 1]), gap=float(gap),
                   part_diams=[float(2 * np.sqrt(sizes[a - 1] / np.pi)), float(2 * np.sqrt(sizes[b - 1] / np.pi))])
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default="data/store")
    ap.add_argument("--census", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()
    with open(a.census) as fh:
        cen = {r["item_id"]: r for r in map(json.loads, fh) if r["nonsimple"]}
    recs = []
    for split in ("train", "val", "test"):
        for it in store.load("lidc", split, a.store):
            c = cen.get(it.item_id)
            if c is None:
                continue
            st = frame(it).masks["nodule"]
            for i in c["voters"]:
                m = size_filter(st[i], c["min_area"])
                g = geometry(m)
                g.update(item_id=it.item_id, split=split, nodule=c["nodule"], rater=i, min_area=c["min_area"],
                         sig=c["sigs"][i], k_vote=len(c["voters"]),
                         n_raters_hole=sum(sum(c["sigs"][j]) > 0 for j in c["voters"]),
                         n_raters_multi=sum(len(c["sigs"][j]) > 1 for j in c["voters"]))
                recs.append(g)
        print(split, len(recs), flush=True)
    with open(a.out, "w") as fh:
        for r in recs:
            fh.write(json.dumps(r) + "\n")
    summ = {}
    q = (10, 25, 50, 75, 90)
    for split in ("train", "val", "test", "all"):
        rs = [r for r in recs if split in ("all", r["split"])]
        hs = [h for r in rs for h in r["holes"]]
        mp = [r for r in rs if r["n_parts"] >= 2]
        pct = lambda v: {p: round(float(np.percentile(v, p)), 3) for p in q} if len(v) else None  # noqa: E731
        summ[split] = {
            "voting_targets_of_nonsimple_items": len(rs),
            "targets_with_holes": sum(r["n_holes"] > 0 for r in rs),
            "targets_multipart": len(mp),
            "hole_count_per_holed_target": {str(k): sum(r["n_holes"] == k for r in rs)
                                            for k in sorted({r["n_holes"] for r in rs if r["n_holes"]})},
            "holes": len(hs),
            "hole_area_px": pct([h["area"] for h in hs]),
            "hole_area_bins_le10_11to30_gt30": [sum(h["area"] <= 10 for h in hs), sum(10 < h["area"] <= 30 for h in hs),
                                                sum(h["area"] > 30 for h in hs)],
            "hole_area_tertiles": [round(float(x), 1) for x in np.percentile([h["area"] for h in hs], [100 / 3, 200 / 3])]
            if hs else None,
            "hole_fraction": pct([h["fraction"] for h in hs]),
            "hole_offset": pct([h["offset"] for h in hs]),
            "hole_wall_px": pct([h["wall"] for h in hs]),
            "holed_part_diam_px": pct([h["part_diam"] for h in hs]),
            "parts_count": {str(k): sum(r["n_parts"] == k for r in mp) for k in sorted({r["n_parts"] for r in mp})},
            "part_ratio": pct([r["part_ratio"] for r in mp]),
            "part_gap_px": pct([r["gap"] for r in mp]),
            "multipart_object_diam_px": pct([r["diam"] for r in mp]),
            "hole_rater_agreement": {f"{k}": sum(r["n_raters_hole"] == k for r in rs if r["n_holes"])
                                     for k in range(1, 5)},
        }
    print(json.dumps(summ, indent=1))
    with open(a.summary, "w") as fh:
        json.dump(summ, fh, indent=1)


if __name__ == "__main__":
    main()
