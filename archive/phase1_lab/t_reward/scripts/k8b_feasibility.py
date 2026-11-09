"""Iteration 3 (R5): is K8b feasible for ANY band of the F5 family?

F5 = h(x) with x = (2D/P)/T5^2. For each item: x_d = x of the vote dilated by 2 tau (all raters),
x_r = max over held-out raters of their LOO x. With a global calibration constant c (equivalent to
changing s*) and shape h, an item satisfies K8b iff h(c x_r) >= 0.85 and h(c x_d) <= 0.50.
Reported: best share of items over c, for h = 1/(1+x^p) (p = 1, 2, 3) and exp(-x); and the share of
items where the held-out rater's x even exceeds the dilated vote's (infeasible for every monotone h).

    ../.venv/bin/python scripts/k8b_feasibility.py
"""
import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results"


def main():
    raw = json.loads((OUT / "iter3_r5_raw.json").read_text())["k8b"]
    rows = json.loads((OUT / "k8_loo.json").read_text())
    xr = {}
    for r in rows:
        f = r.get("F5_prior", r["F5"])
        xr[r["key"]] = max(xr.get(r["key"], 0.0), 1 / max(f, 1e-12) - 1)
    keys = [k for k in raw if k in xr]
    xd = np.array([1 / max(raw[k]["dil_2tau"], 1e-12) - 1 for k in keys])
    xh = np.array([xr[k] for k in keys])
    subs = np.array([raw[k]["subset"] for k in keys])
    print(f"items {len(keys)}; held-out rater worse than the 2-tau-dilated vote (x_r >= x_d): {np.mean(xh >= xd):.1%}")
    shapes = {"1/(1+x)": (1.0, lambda p: p), "1/(1+x^2)": (2.0, None), "1/(1+x^3)": (3.0, None), "exp(-x)": (None, None)}
    cs = np.exp(np.linspace(np.log(1e-3), np.log(1e3), 2001))
    for name, (p, _) in shapes.items():
        if p is None:
            lo_ok, hi_ok = -np.log(0.85), np.log(2.0)  # c x_r <= lo_ok, c x_d >= hi_ok
        else:
            lo_ok, hi_ok = (1 / 0.85 - 1) ** (1 / p), 1.0
        share = [np.mean((c * xh <= lo_ok) & (c * xd >= hi_ok)) for c in cs]
        j = int(np.argmax(share))
        c = cs[j]
        per = {s: np.mean(((c * xh <= lo_ok) & (c * xd >= hi_ok))[subs == s]) for s in sorted(set(subs))}
        print(f"  h = {name:10s} best share {share[j]:.1%} at c = {c:.3g} (c = 1 is F5 as calibrated)  per subset "
              + ", ".join(f"{s} {v:.1%}" for s, v in per.items()))
    print(f"  F5 as calibrated (c = 1, h = 1/(1+x)): {np.mean((xh <= 1 / 0.85 - 1) & (xd >= 1)):.1%}")


if __name__ == "__main__":
    main()


def per_item_bound():
    """Upper bound with a per-item optimal tolerance: an item can satisfy K8b under h iff
    x_d / x_r >= h^-1(0.5) / h^-1(0.85)."""
    raw = json.loads((OUT / "iter3_r5_raw.json").read_text())
    rows = json.loads((OUT / "k8_loo.json").read_text())
    xr = {}
    for r in rows:
        f = r.get("F5_prior", r["F5"])
        xr[r["key"]] = max(xr.get(r["key"], 0.0), 1 / max(f, 1e-12) - 1)
    k8b = raw["k8b"]
    ratio = np.array([(1 / max(k8b[k]["dil_2tau"], 1e-12) - 1) / max(xr[k], 1e-12) for k in k8b])
    print("\nper-item optimal tolerance (upper bound on K8b share): need x_d / x_r >= threshold")
    for name, thr in (("1/(1+x)", 1 / (1 / 0.85 - 1)), ("exp(-x)", np.log(2) / -np.log(0.85)),
                      ("1/(1+x^2)", 1 / np.sqrt(1 / 0.85 - 1)), ("1/(1+x^3)", 1 / (1 / 0.85 - 1) ** (1 / 3))):
        print(f"  {name:10s} threshold {thr:5.2f}: {np.mean(ratio >= thr):.1%} of items")
    print(f"  ratio quantiles p5/p25/p50: {np.percentile(ratio, [5, 25, 50]).round(2).tolist()}")
    sub = raw["sub"]
    print("S_m m=5 exact median:", float(np.median(np.concatenate([sub[k]["5"] for k in sub]))))


if __name__ == "__main__":
    per_item_bound()
