"""C2 power and primary-metric choice on the stratified topology probe (P5.3, P5.5; judge MUST-FIX 1 + 3).

    .venv/bin/python -I -m experiments.phase2.lab.B3.power_c2 --out <lab>/power_c2.json

Simulation (no data): 3 hole-area bins x 120 one-hole items, pooled hole recall = train-share-weighted mean of the
bins (shares from hole_census: 0.4386 / 0.3614 / 0.2000), stratified bootstrap (B = 2,000), R replicates.
Per item i, seed s (2 seeds), the model's probability of emitting the hole at T = 1 is
    q_is = logistic(mu_b + u_i + v_s + w_is),  u_i ~ N(0, su^2) (item difficulty), v_s ~ N(0, sv^2) (seed),
    w_is ~ N(0, 0.5^2) (item x seed),
and the decodes are: greedy = 1[q_is > 0.5] (greedy picks the mode of a binary choice); mean-of-8 = Bin(8, q)/8;
any-of-8 = 1 - Bin(8, q) == 0. Per-item statistic = mean over the 2 seeds. A0 emits no holes, so the paired
(A0' - A0) difference on hole items equals A0''s statistic; its CI is reported as such.
Reports: median CI half-width per metric (true recall 0.5, su in {0.5, 1.5, 3}); per-bin half-widths; the per-bin
floor (point estimate >= 0.30) false-fail rate at true 0.5 and false-pass rate at true 0.15; coverage of the
pooled CI for the arm's seed-population mean when sv in {0, 0.3, 0.6} (2 seeds cannot estimate seed variance).
"""
import argparse
import json

import numpy as np

from . import prereg

W = np.array([0.4386, 0.3614, 0.2000])
N_BIN = prereg.P5_PER_BIN


def sim_stats(mu, su, sv, rng, n_bin=N_BIN, seeds=2):
    """Per-bin arrays of per-item statistics (mean over seeds) for each decode."""
    out = {"greedy": [], "mean8": [], "any8": []}
    v = rng.normal(0, sv, seeds)
    for b in range(3):
        u = rng.normal(0, su, n_bin)
        q = 1 / (1 + np.exp(-(mu[b] + u[:, None] + v[None, :] + rng.normal(0, 0.5, (n_bin, seeds)))))
        k = rng.binomial(8, q)
        out["greedy"].append((q > 0.5).mean(1))
        out["mean8"].append((k / 8).mean(1))
        out["any8"].append((k > 0).mean(1))
    return out


def pooled(per_bin, rng, B=2000):
    est = sum(w * x[rng.integers(0, len(x), (B, len(x)))].mean(1) for w, x in zip(W, per_bin, strict=True))
    point = sum(w * x.mean() for w, x in zip(W, per_bin, strict=True))
    lo, hi = np.percentile(est, [2.5, 97.5])
    return point, lo, hi


def bin_ci(x, rng, B=2000):
    e = x[rng.integers(0, len(x), (B, len(x)))].mean(1)
    lo, hi = np.percentile(e, [2.5, 97.5])
    return (hi - lo) / 2


def mu_for(target, su, metric, rng, n=200_000):
    """mu such that the metric's population mean (no seed effect) is target, by bisection."""
    u = rng.normal(0, su, n) + rng.normal(0, 0.5, n)
    lo, hi = -10.0, 10.0
    for _ in range(50):
        m = (lo + hi) / 2
        q = 1 / (1 + np.exp(-(m + u)))
        val = {"greedy": (q > 0.5).mean(), "mean8": q.mean(), "any8": (1 - (1 - q) ** 8).mean()}[metric]
        lo, hi = (m, hi) if val < target else (lo, m)
    return (lo + hi) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--R", type=int, default=300)
    a = ap.parse_args()
    rng = np.random.default_rng(0)
    res = {"weights": W.tolist(), "n_per_bin": N_BIN, "half_width": {}, "per_bin_half_width": {}, "floor": {},
           "coverage": {}}
    for su in (0.5, 1.5, 3.0):
        hw = {m: [] for m in ("greedy", "mean8", "any8")}
        bw = {m: [] for m in ("greedy", "mean8", "any8")}
        mus = {m: mu_for(0.5, su, m, rng) for m in hw}
        for _ in range(a.R):
            for m in hw:
                st = sim_stats([mus[m]] * 3, su, 0.3, rng)[m]
                _, lo, hi = pooled(st, rng)
                hw[m].append((hi - lo) / 2)
                bw[m].append(np.mean([bin_ci(x, rng) for x in st]))
        res["half_width"][f"su={su}"] = {m: round(float(np.median(v)), 4) for m, v in hw.items()}
        res["per_bin_half_width"][f"su={su}"] = {m: round(float(np.median(v)), 4) for m, v in bw.items()}
        print(su, res["half_width"][f"su={su}"], res["per_bin_half_width"][f"su={su}"], flush=True)
    # floor: false fail at true 0.5, false pass at true 0.15 (su = 1.5, sv = 0.3)
    for m in ("greedy", "mean8"):
        f = {}
        for true in (0.5, 0.4, 0.15, 0.2):
            mu = mu_for(true, 1.5, m, rng)
            est = [sim_stats([mu] * 3, 1.5, 0.3, rng)[m][0].mean() for _ in range(2000)]
            f[f"true={true}"] = {"P(point<0.30)": round(float(np.mean(np.array(est) < prereg.P5_FLOOR)), 4)}
        res["floor"][m] = f
    print("floor", res["floor"], flush=True)
    # coverage of the pooled CI for the seed-population mean (mean8, su = 1.5)
    mu = mu_for(0.5, 1.5, "mean8", rng)
    for sv in (0.0, 0.3, 0.6):
        # the seed-population mean: average over many seeds
        u = rng.normal(0, 1.5, 200_000) + rng.normal(0, 0.5, 200_000) + rng.normal(0, sv, 200_000)
        truth = float((1 / (1 + np.exp(-(mu + u)))).mean())
        cov = []
        for _ in range(a.R):
            _, lo, hi = pooled(sim_stats([mu] * 3, 1.5, sv, rng)["mean8"], rng)
            cov.append(lo <= truth <= hi)
        res["coverage"][f"sv={sv}"] = round(float(np.mean(cov)), 4)
    print("coverage", res["coverage"], flush=True)
    h = res["half_width"]["su=1.5"]
    res["P5_3_pass"] = bool(h["mean8"] <= prereg.P5_HALF_WIDTH)
    g, m8 = h["greedy"], h["mean8"]
    res["P5_5_primary"] = "mean8" if (m8 <= g or abs(m8 - g) <= 0.1 * min(m8, g)) else "greedy"
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps({k: res[k] for k in ("P5_3_pass", "P5_5_primary")}))


if __name__ == "__main__":
    main()
