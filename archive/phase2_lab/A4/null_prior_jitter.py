"""A4 sandbox: the centred-disc prior (null_prior.py constants) on frames whose window is shifted at random so the
raters' union stays inside (translation jitter of the crop, the proposed augmentation and eval protocol).
Run: .venv/bin/python -I null_prior_jitter.py --prior null_prior.json --out null_prior_jitter.json
"""
import argparse
import json
from collections import defaultdict

import numpy as np

from objdel.data import store
from objdel.reward import ExpertCache, score_mask
from objdel.targets import frame
from objdel.targets.frame import _crop

rng = np.random.default_rng(0)


def disc(S, d):
    yy, xx = np.mgrid[:S, :S] + 0.5
    return (xx - S / 2) ** 2 + (yy - S / 2) ** 2 <= (d / 2) ** 2


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--prior", required=True)
    ap.add_argument("--out", default="null_prior_jitter.json")
    a = ap.parse_args()
    priors = json.load(open(a.prior))
    out = {}
    for ds in ("lidc", "riga"):
        prior = {int(k): v for k, v in priors[ds]["prior_diam_by_S"].items()}
        res = defaultdict(list)
        for it in store.load(ds, "val", "data/store"):
            fr = frame(it)
            if fr.scale != 1:
                continue
            S, (x0, y0, _, _) = fr.W, fr.box
            u = np.logical_or.reduce([s.any(0) for s in it.masks.values()])
            ys, xs = np.nonzero(u)
            dx = int(rng.integers(-(x0 + S - xs.max() - 1), xs.min() - x0 + 1))  # keep the union inside
            dy = int(rng.integers(-(y0 + S - ys.max() - 1), ys.min() - y0 + 1))
            pS = prior.get(S) or prior[min(prior, key=lambda s: abs(s - S))]
            for n, st in it.masks.items():
                m = _crop(st, x0 + dx, y0 + dy, S)
                assert m.any(0).sum() == st.any(0).sum()
                res[n].append(float(score_mask(disc(S, pS[n]), ExpertCache.build(m)).iou_tc))
        out[ds] = {n: {"items": len(v), "iou_tc_median": round(float(np.median(v)), 3),
                       "iou_tc_ge_0.5": round(float(np.mean(np.array(v) >= 0.5)), 3)} for n, v in res.items()}
        print(ds, out[ds], flush=True)
    json.dump(out, open(a.out, "w"), indent=1)
