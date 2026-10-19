"""Task 2c: K8 summary from results/k8_loo.json and k8_synth.json (run k8_band.py first).

- pooled Spearman rho(score, area) with a 95 % cluster-bootstrap CI (resampling items);
- synthetic errors: rho(score, area) WITHIN each subset (same protocol / rater count), and
  monotonicity restricted to items where every level could be generated.

    ../.venv/bin/python scripts/k8_summary.py
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

OUT = Path(__file__).resolve().parents[1] / "results"
FORMULAS = ("F0", "F1", "F2", "F3", "F4", "F5")
SUBSETS = ("riga/disc", "riga/cup", "lidc/nodule")


def value(r, f):
    return r.get(f + "_prior", r[f])


def loo_table(rows, rng, n_boot=1000):
    keys = sorted({r["key"] for r in rows})
    by_key = {k: [] for k in keys}
    for j, r in enumerate(rows):
        by_key[r["key"]].append(j)
    area = np.array([r["area"] for r in rows])
    print("K8 LOO (single-reference items use the cross-fitted size prior). rho with 95% cluster-bootstrap CI")
    for f in FORMULAS:
        v = np.array([value(r, f) for r in rows])
        rho = stats.spearmanr(v, area).statistic
        boots = []
        for _ in range(n_boot):
            idx = np.concatenate([by_key[keys[i]] for i in rng.integers(0, len(keys), len(keys))])
            boots.append(stats.spearmanr(v[idx], area[idx]).statistic)
        lo, hi = np.percentile(boots, [2.5, 97.5])
        ok = "PASS" if abs(rho) <= 0.1 and np.median(v) >= 0.9 else "FAIL"
        print(f"  {f}: rho {rho:+.3f} [{lo:+.3f}, {hi:+.3f}]  median {np.median(v):.3f}  mean {v.mean():.3f}"
              f"  p10 {np.percentile(v, 10):.3f}  -> {ok}")
    print("  per k (F5):", {k: f"med {np.median([value(r, 'F5') for r in rows if r['k'] == k]):.3f}"
                           for k in sorted({r['k'] for r in rows})})


def synth_table(rows):
    print("\nSynthetic errors vs all raters: median score and rho(score, area) within subset")
    fams = []
    for r in rows:
        if r["fam"] not in fams:
            fams.append(r["fam"])
    for f in ("F0", "F1", "F5"):
        print(f"  {f}")
        for fam in fams:
            levels = sorted({r["level"] for r in rows if r["fam"] == fam})
            cells = []
            for lv in levels:
                sub_cells = []
                for s in SUBSETS:
                    sel = [r for r in rows if r["fam"] == fam and r["level"] == lv and r["subset"] == s]
                    v = np.array([r[f] for r in sel])
                    a = np.array([r["area"] for r in sel])
                    rho = stats.spearmanr(v, a).statistic if len(sel) > 2 and v.std() > 0 else float("nan")
                    sub_cells.append(f"{np.median(v):.2f}/{rho:+.2f}")
                cells.append(f"{lv}: " + " ".join(sub_cells))
            print(f"    {fam:9s} " + " | ".join(cells))
    print("  (cells: level: disc med/rho, cup med/rho, nodule med/rho)")
    print("\nMonotonicity (items with every level generated; per sign of the level)")
    for f in FORMULAS:
        cells = []
        for fam in fams:
            levels = sorted({r["level"] for r in rows if r["fam"] == fam})
            by_item = {}
            for r in rows:
                if r["fam"] == fam:
                    by_item.setdefault(r["key"], {})[r["level"]] = r[f]
            n = ok = 0
            for v in by_item.values():
                if len(v) != len(levels):
                    continue
                for sgn in (-1, 1):
                    seq = [v[lv] for lv in sorted((lv for lv in levels if np.sign(lv) == sgn), key=abs)]
                    if len(seq) >= 2:
                        n += 1
                        ok += all(b <= a + 1e-12 for a, b in zip(seq, seq[1:]))
            cells.append(f"{fam} {ok}/{n}")
        print(f"  {f}: " + ", ".join(cells))


def main():
    rng = np.random.default_rng(0)
    loo_table(json.loads((OUT / "k8_loo.json").read_text()), rng)
    synth_table(json.loads((OUT / "k8_synth.json").read_text()))


if __name__ == "__main__":
    sys.exit(main())
