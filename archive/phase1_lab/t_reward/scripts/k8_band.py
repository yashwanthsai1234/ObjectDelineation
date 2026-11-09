"""Task 2b, KPI K8: scale-aware band, leave-one-rater-out on the fixtures + synthetic errors.

LOO: each rater's mask is scored against the band of the remaining raters, tau from the
remaining raters' pairwise contour distances (no leakage). When one rater remains there are
no pairs; tau then comes from a size prior tau = c * r_eq^alpha fitted (log-log) on the
OTHER half of the items (2-fold cross-fit by item), or from the floor only.

KPI: Spearman |rho|(band score, object area) <= 0.1 and median score >= 0.9.

    ../.venv/bin/python scripts/k8_band.py
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402
from scipy import stats  # noqa: E402

from rewardcore import fixtures, perturb  # noqa: E402
from rewardcore.band import FORMULAS, expert_sigma  # noqa: E402
from rewardcore.reward import ExpertCache, score  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
FLOOR = 1.0


def fit_prior(rows, field):
    """log-log fit of a per-item spread (rater-pair P95 or F5 sigma) on r_eq (items with >= 2 raters)."""
    r = np.array([x["r_eq"] for x in rows])
    p = np.maximum(np.array([x[field] for x in rows]), 0.3)
    alpha, logc = np.polyfit(np.log(r), np.log(p), 1)
    return float(np.exp(logc)), float(alpha)


def power_law(c, alpha):
    return lambda r: c * r ** alpha


def item_sigma(args):
    return expert_sigma(args[2])


def loo_item(args):
    ds, key, stack, prior, sprior = args
    out = []
    area = float(stack.reshape(len(stack), -1).sum(1).mean())
    for i in range(len(stack)):
        rest = np.delete(stack, i, axis=0)
        cache = ExpertCache.build(rest, FLOOR)
        t = score(stack[i], cache)
        row = dict(subset=fixtures.subset(ds, key), key=key, k=len(stack), rater=i, area=area,
                   tau=cache.tau, tau_raw=cache.tau_raw, iou=t.iou, err_area=t.n_spill + t.n_miss, depth=t.depth,
                   perimeter=cache.perimeter, n_union=cache.n_union, tol4=cache.tol4, tol5=cache.tol5,
                   **{f: t.band[f] for f in FORMULAS})
        if len(rest) == 1:  # size-prior tau instead of the floor
            r_eq = float(np.sqrt(rest.sum() / np.pi))
            cache_p = ExpertCache.build(rest, FLOOR, tau_prior=power_law(*prior), sigma_prior=power_law(*sprior))
            tp = score(stack[i], cache_p)
            row.update({f + "_prior": tp.band[f] for f in FORMULAS}, tau_prior=cache_p.tau, r_eq_rest=r_eq)
        out.append(row)
    return out


PERTS = {
    "grow": [-3, -2, -1, 1, 2, 3],  # px, signed boundary move
    "shift": [1, 2, 3, 4],  # px, diagonal-ish shift (d, d//2)
    "speck_gap": [3, 4, 8, 16],  # 4-px speck centred this far outside the rater union
    "hole_area": [1, 4, 9, 16],  # hole of this many px inside the object
    "grow_tau": [-1.0, -0.5, 0.5, 1.0, 2.0],  # boundary move by c * tau px
    "shift_tau": [0.25, 0.5, 1.0, 2.0],  # shift by c * tau px (diagonal-ish)
}


def synth_item(args):
    ds, key, stack = args
    cache = ExpertCache.build(stack, FLOOR)
    base = cache.vote
    rng = np.random.default_rng(int.from_bytes(key.encode()[:8].ljust(8, b"\0"), "little"))
    rows = []
    area = float(stack.reshape(len(stack), -1).sum(1).mean())
    for fam, levels in PERTS.items():
        for lv in levels:
            if fam == "grow":
                m = perturb.grow(base, lv)
            elif fam == "shift":
                m = perturb.shift(base, lv, lv // 2)
            elif fam == "speck_gap":
                m = perturb.add_speck(base, 4, lv, rng, ref=cache.union)
            elif fam == "hole_area":
                m = perturb.punch_hole(base, lv, rng, centre=True)
            elif fam == "grow_tau":
                m = perturb.grow(base, lv * cache.tau)
            else:
                d = int(round(lv * cache.tau))
                m = perturb.shift(base, d, d // 2) if d else None
            if m is None:
                continue
            t = score(m, cache)
            rows.append(dict(subset=fixtures.subset(ds, key), key=key, k=len(stack), area=area, tau=cache.tau, fam=fam,
                             level=lv, iou=t.iou, err_area=t.n_spill + t.n_miss, depth=t.depth,
                             perimeter=cache.perimeter, **{f: t.band[f] for f in FORMULAS}))
    return rows


def summarize_loo(rows):
    lines = []
    for variant in ("floor", "prior", "k>=3 only"):
        lines.append(f"\n[LOO tau variant: {variant}]  Spearman rho(score, area) | median score | p10 score")
        for f in FORMULAS:
            vals, areas, subs = [], [], []
            for r in rows:
                if variant == "k>=3 only" and r["k"] < 3:
                    continue
                v = r[f + "_prior"] if (variant == "prior" and f + "_prior" in r) else r[f]
                vals.append(v)
                areas.append(r["area"])
                subs.append(r["subset"])
            vals, areas, subs = np.array(vals), np.array(areas), np.array(subs)
            rho = stats.spearmanr(vals, areas).statistic
            parts = [f"{f}: ALL rho {rho:+.3f} med {np.median(vals):.3f} p10 {np.percentile(vals, 10):.3f} (n={len(vals)})"]
            for s in ("riga/disc", "riga/cup", "lidc/nodule"):
                sel = subs == s
                rs = stats.spearmanr(vals[sel], areas[sel]).statistic
                parts.append(f"{s} rho {rs:+.2f} med {np.median(vals[sel]):.3f}")
            lines.append("  " + " | ".join(parts))
    return lines


def summarize_synth(rows):
    lines = ["\n[synthetic errors vs the full rater set] median score by level; rho(score, area) at each level;"
             " share of items where the score is monotone non-increasing in |level|"]
    for f in FORMULAS:
        lines.append(f"  {f}:")
        for fam, levels in PERTS.items():
            cells = []
            for lv in levels:
                sel = [r for r in rows if r["fam"] == fam and r["level"] == lv]
                v = np.array([r[f] for r in sel])
                a = np.array([r["area"] for r in sel])
                cells.append(f"{lv}: {np.median(v):.2f} (rho {stats.spearmanr(v, a).statistic:+.2f})")
            # monotonicity per item
            by_item = {}
            for r in rows:
                if r["fam"] == fam:
                    by_item.setdefault(r["key"], []).append((abs(r["level"]), np.sign(r["level"]), r[f]))
            mono = n = 0
            for v in by_item.values():
                for sgn in (-1, 1):
                    seq = [s for lv, sg, s in sorted(v) if sg == sgn]
                    if len(seq) >= 2:
                        n += 1
                        mono += all(b <= a + 1e-12 for a, b in zip(seq, seq[1:]))
            lines.append(f"    {fam:10s} " + " | ".join(cells) + f" | monotone {mono}/{n}")
    return lines


def main():
    items = fixtures.load_all()
    h5 = {(r["key"]): r for r in __import__("csv").DictReader(open(OUT / "h5_items.csv"))}
    for r in h5.values():
        r["r_eq"], r["p95"] = float(r["r_eq"]), float(r["p95"])
    with Pool(8) as pool:
        for (_, k, _), sig in zip(items, pool.map(item_sigma, items, chunksize=8)):
            h5[k]["sigma"] = sig
        # 2-fold cross-fit of the size priors: fold by item index parity; fold p uses a fit on the other fold
        folds = [[h5[k] for j, (_, k, _) in enumerate(items) if j % 2 == p] for p in (0, 1)]
        priors = [fit_prior(folds[1], "p95"), fit_prior(folds[0], "p95")]
        spriors = [fit_prior(folds[1], "sigma"), fit_prior(folds[0], "sigma")]
        print("size priors (c, alpha) per fold: P95", priors, "sigma", spriors)
        args = [(ds, k, st, priors[j % 2], spriors[j % 2]) for j, (ds, k, st) in enumerate(items)]
        loo = [r for rows in pool.map(loo_item, args, chunksize=8) for r in rows]
        syn = [r for rows in pool.map(synth_item, [(ds, k, st) for ds, k, st in items], chunksize=8) for r in rows]
    (OUT / "k8_loo.json").write_text(json.dumps(loo))
    (OUT / "k8_synth.json").write_text(json.dumps(syn))
    taus = np.array([r["tau"] for r in loo])
    print(f"LOO cases: {len(loo)}; tau median {np.median(taus):.2f} px, share at floor {np.mean(taus == FLOOR):.2f}")
    for line in summarize_loo(loo) + summarize_synth(syn):
        print(line)


if __name__ == "__main__":
    main()
