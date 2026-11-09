"""Iteration 3 (R5): K8 re-registered per stratum, K8b band discrimination, RIGA rater-subsample test.

K8  (per stratum): |Spearman rho(F5, area)| <= 0.10 in riga/disc, riga/cup, lidc k=2, k=3, k=4 (LOO,
    from results/k8_loo.json; k=2 = single remaining rater). Pooled rho reported, not gated.
    The k=2 stratum is also re-scored with a LIDC-only sigma prior (2-fold cross-fit by item).
K8b (discrimination), on >= 95 % of items: F5 of the vote moved by a uniform 2 tau excursion
    (dilation; erosion reported too) <= 0.50; F5 of every held-out rater (LOO) >= 0.85; and F5 is
    monotone non-increasing in |excursion| (0.25..2 tau, both directions) on >= 98 % of items.
S_m test: RIGA (6 raters), m in {2, 3, 4, 5} reference raters drawn at random + 1 held-out expert
    (3 draws per image); median F5 of the held-out expert in [0.92, 0.98] for every m.

    ../.venv/bin/python scripts/k8_r5.py
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from rewardcore import fixtures, perturb  # noqa: E402
from rewardcore.band import excursion_second_moment, expert_sigma, tolerance_factor  # noqa: E402
from rewardcore.reward import ExpertCache, score  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
LEVELS = (0.25, 0.5, 1.0, 1.5, 2.0)
X_STAR = 1 / 0.95 - 1


def strata(rows, f="F5"):
    out = {}
    for name, sel in (("riga/disc", lambda r: r["subset"] == "riga/disc"), ("riga/cup", lambda r: r["subset"] == "riga/cup"),
                      ("lidc k=2", lambda r: r["subset"] == "lidc/nodule" and r["k"] == 2),
                      ("lidc k=3", lambda r: r["subset"] == "lidc/nodule" and r["k"] == 3),
                      ("lidc k=4", lambda r: r["subset"] == "lidc/nodule" and r["k"] == 4)):
        g = [r for r in rows if sel(r)]
        v = np.array([r.get(f + "_prior", r[f]) for r in g])
        a = np.array([r["area"] for r in g])
        out[name] = (len(g), stats.spearmanr(v, a).statistic, np.median(v))
    return out


def lidc_sigma_prior_k2(rows):
    """Re-score the k=2 LOO cases with a sigma prior fitted on LIDC items only (other fold)."""
    lidc = fixtures.load("lidc")
    keys = sorted(lidc)
    sig = {k: expert_sigma(lidc[k]) for k in keys}
    r_eq = {k: float(np.sqrt(lidc[k].reshape(len(lidc[k]), -1).sum(1).mean() / np.pi)) for k in keys}
    fold = {k: i % 2 for i, k in enumerate(keys)}
    fits = {}
    for p in (0, 1):
        ks = [k for k in keys if fold[k] != p]
        alpha, logc = np.polyfit(np.log([r_eq[k] for k in ks]), np.log([max(sig[k], 0.3) for k in ks]), 1)
        fits[p] = (np.exp(logc), alpha)
    out = []
    for r in rows:
        if r["subset"] != "lidc/nodule" or r["k"] != 2:
            continue
        c, a = fits[fold[r["key"]]]
        t5 = max(c * r["r_eq_rest"] ** a * np.sqrt(excursion_second_moment(1) / X_STAR), tolerance_factor(1))
        out.append((r["area"], 1 / (1 + (2 * r["depth"] / r["perimeter"]) / t5 ** 2)))
    a, v = np.array(out).T
    return fits, len(v), stats.spearmanr(v, a).statistic, np.median(v)


def k8b_item(args):
    ds, key, st = args
    cache = ExpertCache.build(st)
    out = {"subset": fixtures.subset(ds, key)}
    for sign, name in ((1, "dil"), (-1, "ero")):
        seq = [score(perturb.grow(cache.vote, sign * c * cache.tau), cache).band["F5"] for c in LEVELS]
        out[name + "_2tau"] = seq[-1]
        out[name + "_mono"] = all(b <= a + 1e-12 for a, b in zip([1.0] + seq, seq))
    return key, out


def subsample_item(args):
    key, st, seed = args
    rng = np.random.default_rng(seed)
    res = {}
    for m in (2, 3, 4, 5):
        vals = []
        for _ in range(3):
            pick = rng.permutation(len(st))[:m + 1]
            cache = ExpertCache.build(st[pick[1:]])
            vals.append(score(st[pick[0]], cache).band["F5"])
        res[m] = vals
    return key, res


def main():
    rows = json.loads((OUT / "k8_loo.json").read_text())
    print("K8 per stratum (LOO; |rho| <= 0.10 each). pooled rho reported, not gated:",
          f"{stats.spearmanr([r.get('F5_prior', r['F5']) for r in rows], [r['area'] for r in rows]).statistic:+.3f}")
    for f in ("F5", "F4", "F0"):
        st = strata(rows, f)
        print(f"  {f}: " + " | ".join(f"{k} n={n} rho {rho:+.3f} med {med:.3f} {'PASS' if abs(rho) <= 0.10 else 'FAIL'}"
                                      for k, (n, rho, med) in st.items()))
    fits, n, rho, med = lidc_sigma_prior_k2(rows)
    print(f"  F5 lidc k=2 with a LIDC-only cross-fitted sigma prior {fits}: n={n} rho {rho:+.3f} med {med:.3f}"
          f" {'PASS' if abs(rho) <= 0.10 else 'FAIL'}")

    items = fixtures.load_all()
    with Pool(1) as pool:
        k8b = dict(pool.map(k8b_item, items, chunksize=16))
    loo_min = {}
    for r in rows:
        v = r.get("F5_prior", r["F5"])
        loo_min[r["key"]] = min(loo_min.get(r["key"], 1.0), v)
    keys = list(k8b)
    dil = np.array([k8b[k]["dil_2tau"] for k in keys])
    ero = np.array([k8b[k]["ero_2tau"] for k in keys])
    mono = np.array([k8b[k]["dil_mono"] and k8b[k]["ero_mono"] for k in keys])
    rmin = np.array([loo_min[k] for k in keys])
    print("\nK8b discrimination (share of items)")
    print(f"  F5(vote dilated 2 tau) <= 0.50: {np.mean(dil <= 0.5):.1%} (median {np.median(dil):.3f})  [target >= 95 %]"
          f" {'PASS' if np.mean(dil <= 0.5) >= 0.95 else 'FAIL'}")
    print(f"  F5(vote eroded 2 tau)  <= 0.50: {np.mean(ero <= 0.5):.1%} (median {np.median(ero):.3f})  (reported)")
    print(f"  min over held-out raters F5 >= 0.85: {np.mean(rmin >= 0.85):.1%} (median of item-min {np.median(rmin):.3f})"
          f"  [target >= 95 %] {'PASS' if np.mean(rmin >= 0.85) >= 0.95 else 'FAIL'}")
    print(f"  monotone in |excursion| (dilate and erode): {np.mean(mono):.1%}  [target >= 98 %] {'PASS' if np.mean(mono) >= 0.98 else 'FAIL'}")
    subs = np.array([k8b[k]["subset"] for k in keys])
    for s in ("riga/disc", "riga/cup", "lidc/nodule"):
        m = subs == s
        print(f"    {s:12s} dil<=0.5 {np.mean(dil[m] <= 0.5):.1%}  rater-min>=0.85 {np.mean(rmin[m] >= 0.85):.1%}  mono {np.mean(mono[m]):.1%}")
    # per held-out rater (not per item) for context
    allv = np.array([r.get("F5_prior", r["F5"]) for r in rows])
    print(f"  per held-out rater (all LOO cases): share >= 0.85 {np.mean(allv >= 0.85):.1%}")

    riga = fixtures.load("riga")
    with Pool(1) as pool:
        sub = dict(pool.map(subsample_item, [(k, st, j) for j, (k, st) in enumerate(riga.items())], chunksize=8))
    print("\nS_m correction, RIGA rater subsamples: median F5 of a held-out expert vs m reference raters [0.92, 0.98]")
    for m in (2, 3, 4, 5):
        v = np.concatenate([sub[k][m] for k in sub])
        vd = np.concatenate([sub[k][m] for k in sub if k.endswith("disc")])
        vc = np.concatenate([sub[k][m] for k in sub if k.endswith("cup")])
        med = np.median(v)
        print(f"  m={m}: median {med:.3f} (disc {np.median(vd):.3f}, cup {np.median(vc):.3f}; p10 {np.percentile(v, 10):.3f})"
              f" {'PASS' if 0.92 <= med <= 0.98 else 'FAIL'}")
    (OUT / "iter3_r5_raw.json").write_text(json.dumps({"k8b": k8b, "sub": sub}))


if __name__ == "__main__":
    main()
