"""Iteration 2, task 3: K8 with rater-count adjustment (from results/k8_loo.json).

(1) partial Spearman: Pearson correlation of rank(score) and rank(area) after linearly removing
    rank(k) from both (k = raters per item);
(2) stratified: Spearman rho within each k stratum (k = 2, 3, 4 LIDC; 6 RIGA), combined with
    n-weighted Fisher z.
95 % CIs: cluster bootstrap over items (1000 resamples). Single-reference LOO cases (k = 2) use the
cross-fitted size prior, as in iteration 1.

    ../.venv/bin/python scripts/k8_adjusted.py
"""
import json
from pathlib import Path

import numpy as np
from scipy import stats

OUT = Path(__file__).resolve().parents[1] / "results"
FORMULAS = ("F0", "F1", "F2", "F3", "F4", "F5")


def partial_spearman(v, a, k):
    rv, ra, rk = (stats.rankdata(x) for x in (v, a, k))
    X = np.column_stack([np.ones_like(rk), rk])
    res = lambda y: y - X @ np.linalg.lstsq(X, y, rcond=None)[0]  # noqa: E731
    return float(np.corrcoef(res(rv), res(ra))[0, 1])


def stratified(v, a, k):
    zs, ns, per = [], [], {}
    for kk in np.unique(k):
        m = k == kk
        if m.sum() > 3 and np.std(v[m]) > 0:
            r = stats.spearmanr(v[m], a[m]).statistic
            per[int(kk)] = r
            zs.append(np.arctanh(np.clip(r, -0.999, 0.999)))
            ns.append(m.sum())
    return float(np.tanh(np.average(zs, weights=ns))), per


def main():
    rows = json.loads((OUT / "k8_loo.json").read_text())
    keys = sorted({r["key"] for r in rows})
    idx_of = {kk: [] for kk in keys}
    for j, r in enumerate(rows):
        idx_of[r["key"]].append(j)
    a = np.array([r["area"] for r in rows])
    k = np.array([r["k"] for r in rows])
    rng = np.random.default_rng(0)
    boots_idx = [np.concatenate([idx_of[keys[i]] for i in rng.integers(0, len(keys), len(keys))]) for _ in range(1000)]
    print(f"LOO cases {len(rows)}, items {len(keys)}, strata k: {dict(zip(*np.unique(k, return_counts=True)))}")
    print(f"{'formula':7s} {'pooled rho':>10s} {'partial rho | k [95% CI]':>28s} {'stratified rho [95% CI]':>27s}  per-stratum rho (k=2,3,4,6)  median  verdict")
    for f in FORMULAS:
        v = np.array([r.get(f + "_prior", r[f]) for r in rows])
        pooled = stats.spearmanr(v, a).statistic
        pr = partial_spearman(v, a, k)
        st, per = stratified(v, a, k)
        pb = [partial_spearman(v[i], a[i], k[i]) for i in boots_idx]
        sb = [stratified(v[i], a[i], k[i])[0] for i in boots_idx]
        (p_lo, p_hi), (s_lo, s_hi) = np.percentile(pb, [2.5, 97.5]), np.percentile(sb, [2.5, 97.5])
        ok = "PASS" if max(abs(pr), abs(st)) <= 0.1 and np.median(v) >= 0.9 else "FAIL"
        print(f"{f:7s} {pooled:+10.3f} {pr:+10.3f} [{p_lo:+.3f}, {p_hi:+.3f}] {st:+10.3f} [{s_lo:+.3f}, {s_hi:+.3f}]  "
              + " ".join(f"{per.get(kk, float('nan')):+.2f}" for kk in (2, 3, 4, 6)) + f"   {np.median(v):.3f}  {ok}")


if __name__ == "__main__":
    main()
