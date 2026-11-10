"""Sensitivity of encodings to one wrong digit / one +-1 unit error. Usage: -I perturb.py data/targets out.json"""

import glob
import json
import os
import sys

import numpy as np

from objdel.geometry import SCHEMAS, parse, to_mask

rng = np.random.default_rng(0)
root, outp = sys.argv[1], sys.argv[2]
recs = []
for f in sorted(glob.glob(os.path.join(root, "*", "val.jsonl"))):
    ds = f.split(os.sep)[-2]
    recs += [dict(json.loads(l), ds=ds) for l in open(f)]
idx = rng.choice(len(recs), 2000, replace=False)


def iou(a, b):
    u = (a | b).sum()
    return 1.0 if u == 0 else (a & b).sum() / u


res = {}
for i in idx:
    r = recs[i]
    names = SCHEMAS[r["ds"]]
    p = parse(r["text"], names)
    W = r["W"]
    k = 1000 // W
    for n in names:
        parts = p.objs[n]
        if not parts:
            continue
        ref = to_mask(parts, W, W).mask
        # pick one ring (outer of a random part) and one coordinate
        pi = rng.integers(len(parts))
        ring = parts[pi][0]
        vi = rng.integers(len(ring))
        ax = rng.integers(2)
        for enc in ("abs1k", "abs_px", "delta1k", "delta_px"):
            unit = k if enc.endswith("px") else 1
            a = ring // unit if enc.endswith("px") else ring.copy()
            e = a if enc.startswith("abs") else np.vstack([a[:1], np.diff(a, axis=0)])
            for err in ("pm1", "digit"):
                e2 = e.copy()
                v = int(e2[vi, ax])
                if err == "pm1":
                    v2 = v + int(rng.choice([-1, 1]))
                else:
                    s = str(abs(v))
                    j = rng.integers(len(s))
                    d = rng.choice([c for c in "0123456789" if c != s[j]])
                    s2 = (s[:j] + d + s[j + 1 :]).lstrip("0") or "0"
                    v2 = int(s2) * (-1 if v < 0 else 1)
                e2[vi, ax] = v2
                a2 = e2 if enc.startswith("abs") else np.cumsum(e2, axis=0)
                ring2 = np.clip(a2 * unit, 0, 1000)
                parts2 = list(parts)
                parts2[pi] = (ring2, parts[pi][1])
                g = to_mask(parts2, W, W)
                key = f"{r['ds']}/{enc}/{err}"
                res.setdefault(key, []).append((iou(ref, g.mask), g.legal))
out = {}
for key, v in sorted(res.items()):
    a = np.array(v, float)
    out[key] = dict(n=len(a), iou_mean=a[:, 0].mean(), iou_p5=np.percentile(a[:, 0], 5), illegal=1 - a[:, 1].mean())
    print(
        f"{key:24s} n={len(a):4d} IoU mean {a[:, 0].mean():.4f} p5 {np.percentile(a[:, 0], 5):.4f} illegal {1 - a[:, 1].mean():.3f}"
    )
json.dump(out, open(outp, "w"), indent=1)
