"""Task 2a, hypothesis H5: does inter-rater boundary disagreement scale with object size?

For every fixture item and every rater pair: P95 and mean of the pooled symmetric contour
distances. Regress log(P95) on log(r_eq), r_eq = sqrt(mean rater area / pi).
slope ~ 0: disagreement is flat in pixels; slope ~ 1: proportional to size.

    ../.venv/bin/python scripts/h5_disagreement.py
"""
import csv
import sys
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from rewardcore import fixtures  # noqa: E402
from rewardcore.distances import symmetric_surface_distances  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"


def item_rows():
    rows = []
    for ds, key, stack in fixtures.load_all():
        area = float(stack.reshape(len(stack), -1).sum(1).mean())
        p95s, means = [], []
        for a, b in combinations(stack, 2):
            d = symmetric_surface_distances(a, b)
            p95s.append(float(np.percentile(d, 95)))
            means.append(float(d.mean()))
        rows.append(dict(subset=fixtures.subset(ds, key), key=key, k=len(stack), area=area,
                         r_eq=float(np.sqrt(area / np.pi)), p95=float(np.mean(p95s)), assd=float(np.mean(means))))
    return rows


def slope_ci(x, y, rng, n_boot=2000):
    lx, ly = np.log(x), np.log(y)
    slope = np.polyfit(lx, ly, 1)[0]
    idx = rng.integers(0, len(x), (n_boot, len(x)))
    boots = [np.polyfit(lx[i], ly[i], 1)[0] for i in idx]
    return slope, np.percentile(boots, [2.5, 97.5])


def main():
    rows = item_rows()
    OUT.mkdir(exist_ok=True)
    with open(OUT / "h5_items.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    rng = np.random.default_rng(0)
    groups = {s: [r for r in rows if r["subset"] == s] for s in ("riga/disc", "riga/cup", "lidc/nodule")}
    groups["all"] = rows
    print("Mean over rater pairs of the per-pair P95 / mean symmetric contour distance (px)")
    print(f"{'subset':12s} {'n':>4} {'r_eq med [p10,p90]':>22} {'P95 med [p10,p90]':>20} {'P95/r_eq med':>12} "
          f"{'rho(P95,area)':>13} {'slope logP95~logr [95% CI]':>28} {'slope ASSD':>10}")
    for name, g in groups.items():
        r = np.array([x["r_eq"] for x in g])
        p = np.maximum(np.array([x["p95"] for x in g]), 0.5)  # 0 when raters coincide; floor for the log
        a = np.array([x["assd"] for x in g]) + 0.05
        sl, ci = slope_ci(r, p, rng)
        sla, _ = slope_ci(r, a, rng, 200)
        rho = stats.spearmanr(p, r).statistic
        print(f"{name:12s} {len(g):4d} {np.median(r):7.1f} [{np.percentile(r, 10):5.1f},{np.percentile(r, 90):6.1f}]"
              f" {np.median(p):6.2f} [{np.percentile(p, 10):5.2f},{np.percentile(p, 90):6.2f}] {np.median(p / r):12.3f}"
              f" {rho:13.2f} {sl:10.2f} [{ci[0]:.2f}, {ci[1]:.2f}] {sla:10.2f}")
    print("\nBinned by r_eq (quintiles within subset): median P95 px | median P95/r_eq")
    for name in ("riga/disc", "riga/cup", "lidc/nodule"):
        g = groups[name]
        r = np.array([x["r_eq"] for x in g])
        p = np.array([x["p95"] for x in g])
        edges = np.percentile(r, [0, 20, 40, 60, 80, 100])
        cells = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            sel = (r >= lo) & (r <= hi)
            cells.append(f"r {np.median(r[sel]):5.1f}: {np.median(p[sel]):5.2f} | {np.median(p[sel] / r[sel]):.3f}")
        print(f"  {name:12s} " + "   ".join(cells))
    lidc = groups["lidc/nodule"]
    for k in (2, 3, 4):
        g = [x for x in lidc if x["k"] == k]
        print(f"  lidc k={k}: n={len(g)}, median P95 {np.median([x['p95'] for x in g]):.2f} px, median r_eq {np.median([x['r_eq'] for x in g]):.1f}")


if __name__ == "__main__":
    main()
