"""A4 sandbox: image-blind null models against the Phase 2 geometry gates.

The model frame is centred on the raters' union box with a margin (targets.frame), so the crop itself tells a model
where the object is and, through the side S, roughly how big it is. 'prior_disc': a disc (RIGA: disc and cup
discs) at the frame centre whose diameter is the TRAIN median diameter for that frame side S (and, for the cup, the
train median cup/disc ratio). Never looks at the image. Scored with objdel.reward.score_mask (IoU_tc) on val, the
same path as the gates (median IoU_tc; share with IoU_tc >= 0.5).
Run: .venv/bin/python -I null_prior.py --out null_prior.json
"""
import argparse
import json
from collections import defaultdict

import numpy as np

from objdel.data import store
from objdel.reward import ExpertCache, score_mask
from objdel.targets import frame


def diam(m):
    return 2 * np.sqrt(m.sum() / np.pi)


def disc(S, d):
    yy, xx = np.mgrid[:S, :S] + 0.5
    return (xx - S / 2) ** 2 + (yy - S / 2) ** 2 <= (d / 2) ** 2


def fit(ds, every):
    by = defaultdict(lambda: defaultdict(list))
    for i, it in enumerate(store.load(ds, "train", "data/store")):
        if i % every:
            continue
        fr = frame(it)
        for n, st in fr.masks.items():
            for m in st:
                if m.sum() >= 5:
                    by[fr.W][n].append(diam(m))
    return {S: {n: float(np.median(v)) for n, v in d.items()} for S, d in by.items()}


def evaluate(ds, prior):
    res = defaultdict(list)
    for it in store.load(ds, "val", "data/store"):
        fr = frame(it)
        S = fr.W
        pS = prior.get(S) or prior[min(prior, key=lambda s: abs(s - S))]
        for n, st in fr.masks.items():
            c = ExpertCache.build(st)
            res[n].append(float(score_mask(disc(S, pS[n]), c).iou_tc))
    return {n: {"items": len(v), "iou_tc_median": round(float(np.median(v)), 3),
                "iou_tc_ge_0.5": round(float(np.mean(np.array(v) >= 0.5)), 3)} for n, v in res.items()}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="null_prior.json")
    a = ap.parse_args()
    out = {}
    for ds, every in (("lidc", 10), ("riga", 4)):
        prior = fit(ds, every)
        out[ds] = {"prior_diam_by_S": {str(k): v for k, v in sorted(prior.items())}, "val": evaluate(ds, prior)}
        print(ds, json.dumps(out[ds], indent=1), flush=True)
    json.dump(out, open(a.out, "w"), indent=1)
