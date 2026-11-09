"""Compare two mask directories (0/127/254 PNGs): bit-exact equality and disc/cup IoU per common file.

    python -I compare.py <dir_a> <dir_b> <out.json> [name_filter_substring]
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from PIL import Image


def iou(a, b):
    u = np.logical_or(a, b).sum()
    return float(np.logical_and(a, b).sum() / u) if u else 1.0


def one(args):
    a, b = args
    x, y = np.asarray(Image.open(a)), np.asarray(Image.open(b))
    return {"file": a.name, "exact": bool(np.array_equal(x, y)),
            "disc_iou": iou(x >= 100, y >= 100), "cup_iou": iou(x >= 200, y >= 200)}


def main(da, db, out, filt=""):
    A = {p.name: p for p in Path(da).glob("*.png") if filt in p.name}
    B = {p.name: p for p in Path(db).glob("*.png") if filt in p.name}
    common = sorted(set(A) & set(B))
    with Pool(7) as p:
        rows = p.map(one, [(A[n], B[n]) for n in common], chunksize=16)
    d = np.array([r["disc_iou"] for r in rows]) if rows else np.array([1.0])
    c = np.array([r["cup_iou"] for r in rows]) if rows else np.array([1.0])
    s = {"only_a": len(set(A) - set(B)), "only_b": len(set(B) - set(A)), "common": len(common),
         "exact": sum(r["exact"] for r in rows),
         "disc_iou": {"min": float(d.min()), "p1": float(np.percentile(d, 1)), "median": float(np.median(d))},
         "cup_iou": {"min": float(c.min()), "p1": float(np.percentile(c, 1)), "median": float(np.median(c))},
         "worst": sorted(rows, key=lambda r: min(r["disc_iou"], r["cup_iou"]))[:8],
         "only_b_files": sorted(set(B) - set(A))[:40], "only_a_files_sample": sorted(set(A) - set(B))[:10]}
    Path(out).write_text(json.dumps(s, indent=1))
    print(json.dumps({k: v for k, v in s.items() if k not in ("worst",)}, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:])
