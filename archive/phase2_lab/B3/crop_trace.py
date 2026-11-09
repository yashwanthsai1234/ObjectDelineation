"""H-A6-9 sandbox: trace rings on the mask cropped to its bounding box (crop before polygonize / connectedComponents)
inside objdel.targets.build_target, by monkeypatching objdel.targets.build._rings in this process only (src/ is not
edited). Pass criteria: prereg.py (P4).

    .venv/bin/python -I -m experiments.phase2.lab.B3.crop_trace --out <lab>/crop_trace.json

Every train/val target of RIGA and LIDC (as experiments/phase1/build_targets: frame, min_areas, BuildConfig()) is
built twice in one process, baseline and cropped, in alternating order per target (ABBA) so machine load hits both
arms alike. Reports per-file and overall target_set_hash of both arms and the paired per-target time ratio.

Why the texts cannot change: polygonize already takes an offset; boundary_edges pads by one pixel itself, so the
crop needs no margin; the ring start key (y * (W + 1) + x) * 4 + d orders vertices lexicographically by (y, x, d)
for any width, so a translation keeps every ring's start and order; cv2.connectedComponents labels in raster order,
which a crop preserves (labels only group holes with their outer ring anyway).
"""
import argparse
import json
import time

import cv2
import numpy as np

import objdel.targets.build as B
from objdel.checker import polygonize
from objdel.data import store
from objdel.geometry import SCHEMAS
from objdel.targets import BuildConfig, frame, min_areas, target_set_hash
from objdel.targets.simplify import lattice_path

from . import prereg

BASE = B._rings


def rings_cropped(mask):
    m = np.ascontiguousarray(mask, dtype=bool)
    x0, y0, w, h = cv2.boundingRect(m.view(np.uint8))
    if w == 0:
        return []
    sub = np.ascontiguousarray(m[y0:y0 + h, x0:x0 + w])
    _, label = cv2.connectedComponents(sub.view(np.uint8), connectivity=8)
    out = []
    for r in polygonize(sub, offset=(x0, y0)):
        path, dirs = lattice_path(r.pts)
        x, y = path[0] + B._RIGHT[dirs[0]]
        out.append((path, dirs, int(label[y - y0, x - x0]), r.hole))
    return out


def build(src, names, ma, cropped):
    B._rings = rings_cropped if cropped else BASE
    t0 = time.perf_counter()
    t = B.build_target(src, names, BuildConfig(), ma)
    return t.text, time.perf_counter() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default="data/store")
    ap.add_argument("--out", required=True)
    ap.add_argument("--every", type=int, default=1, help="every n-th item (1 = all)")
    a = ap.parse_args()
    res, rows_all = {"files": {}}, {"base": [], "crop": []}
    t_all = time.time()
    for ds in ("riga", "lidc"):
        names = SCHEMAS[ds]
        for split in ("train", "val"):
            rows = {"base": [], "crop": []}
            times, k = [], 0
            for i, it in enumerate(store.load(ds, split, a.store)):
                if i % a.every:
                    continue
                f = frame(it)
                ma = min_areas(f.masks)
                for r, rater in enumerate(it.raters):
                    src = {n: f.masks[n][r] for n in names}
                    order = (False, True) if k % 2 == 0 else (True, False)
                    k += 1
                    got = {}
                    for cropped in order:
                        try:
                            got[cropped] = build(src, names, ma, cropped)
                        except ValueError as e:
                            got[cropped] = (f"error: {e}", 0.0)
                    rows["base"].append((it.item_id, rater, got[False][0]))
                    rows["crop"].append((it.item_id, rater, got[True][0]))
                    times.append((got[False][1], got[True][1]))
            B._rings = BASE
            rel = f"{ds}/{split}"
            t = np.array(times)
            ratio = t[:, 1] / t[:, 0]
            e = {"targets": len(times), "hash_base": target_set_hash(rows["base"]),
                 "hash_crop": target_set_hash(rows["crop"]),
                 "texts_differ": sum(x != y for x, y in zip(rows["base"], rows["crop"], strict=True)),
                 "base_ms_median": round(float(np.median(t[:, 0])) * 1e3, 3),
                 "crop_ms_median": round(float(np.median(t[:, 1])) * 1e3, 3),
                 "paired_ratio_median": round(float(np.median(ratio)), 4),
                 "paired_ratio_p10_p90": [round(float(np.percentile(ratio, q)), 4) for q in (10, 90)],
                 "total_s_base": round(float(t[:, 0].sum()), 1), "total_s_crop": round(float(t[:, 1].sum()), 1)}
            if a.every == 1:
                e["base_matches_data_targets"] = e["hash_base"] == prereg.P4_FILE_HASH[rel]
            res["files"][rel] = e
            rows_all["base"] += rows["base"]
            rows_all["crop"] += rows["crop"]
            print(rel, json.dumps(e), f"{time.time() - t_all:.0f}s", flush=True)
    res["overall_base"] = target_set_hash(rows_all["base"])
    res["overall_crop"] = target_set_hash(rows_all["crop"])
    res["verdict"] = {
        "hash_unchanged": res["overall_crop"] == res["overall_base"]
        and (a.every != 1 or res["overall_crop"] == prereg.P4_HASH),
        "riga_speedup": 1 - res["files"]["riga/train"]["paired_ratio_median"] >= prereg.P4_RIGA_SPEEDUP}
    res["wall_s"] = round(time.time() - t_all, 1)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res["verdict"]), res["overall_base"], res["overall_crop"])


if __name__ == "__main__":
    main()
