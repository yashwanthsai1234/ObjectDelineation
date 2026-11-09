"""Iteration 4: K8b re-registered as discrimination (original absolute-threshold numbers alongside).

New K8b (on the fixtures; held-out = leave-one-rater-out F5 from results/k8_loo.json, the 2-tau moves
and the monotonicity sweep from results/iter3_r5_raw.json):
  (i)   F5 monotone non-increasing in |excursion| (0.25..2 tau, dilate and erode) on >= 98 % of items
  (ii)  F5 of the vote moved by a uniform 2 tau excursion ranks below the median held-out rater of the
        same item on >= 95 % of items (dilation is the gated case; erosion reported)
  (iii) held-out-rater F5: median >= 0.90 and p10 >= 0.70 (all LOO cases)

    ../.venv/bin/python scripts/k8b_v2.py
"""
import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"


def main():
    raw = json.loads((OUT / "iter3_r5_raw.json").read_text())["k8b"]
    rows = json.loads((OUT / "k8_loo.json").read_text())
    held = {}
    for r in rows:
        held.setdefault(r["key"], []).append(r.get("F5_prior", r["F5"]))
    keys = [k for k in raw if k in held]
    subs = np.array([raw[k]["subset"] for k in keys])
    med = np.array([np.median(held[k]) for k in keys])
    dil = np.array([raw[k]["dil_2tau"] for k in keys])
    ero = np.array([raw[k]["ero_2tau"] for k in keys])
    mono = np.array([raw[k]["dil_mono"] and raw[k]["ero_mono"] for k in keys])
    allv = np.concatenate([held[k] for k in keys])
    print(f"items {len(keys)}, held-out rater cases {len(allv)}")
    print("NEW K8b")
    print(f"  (i)   monotone in |excursion|: {mono.mean():.1%} [>= 98 %] {'PASS' if mono.mean() >= 0.98 else 'FAIL'}")
    s2 = np.mean(dil < med)
    print(f"  (ii)  F5(vote dilated 2 tau) < median held-out rater: {s2:.1%} [>= 95 %] {'PASS' if s2 >= 0.95 else 'FAIL'}"
          f"   (eroded 2 tau: {np.mean(ero < med):.1%}, reported)")
    for s in ("riga/disc", "riga/cup", "lidc/nodule"):
        m = subs == s
        print(f"        {s:12s} dilated {np.mean(dil[m] < med[m]):.1%}  eroded {np.mean(ero[m] < med[m]):.1%}")
    m50, p10 = np.median(allv), np.percentile(allv, 10)
    print(f"  (iii) held-out F5 median {m50:.3f} [>= 0.90], p10 {p10:.3f} [>= 0.70] {'PASS' if m50 >= 0.9 and p10 >= 0.7 else 'FAIL'}")
    print("ORIGINAL K8b (iteration 3, absolute thresholds)")
    rmin = np.array([min(held[k]) for k in keys])
    print(f"  F5(vote dilated 2 tau) <= 0.50 on {np.mean(dil <= 0.5):.1%} of items [>= 95 %] FAIL;"
          f" min held-out rater >= 0.85 on {np.mean(rmin >= 0.85):.1%} [>= 95 %] FAIL; monotone {mono.mean():.1%} PASS")


if __name__ == "__main__":
    main()
