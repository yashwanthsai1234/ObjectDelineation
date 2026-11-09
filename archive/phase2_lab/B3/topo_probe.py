"""Topology probe design and power (H-A1-7): candidate probes, their train cost, and cluster-bootstrap CIs of the probe
metrics on synthetic predictions at plausible rates and within-nodule correlations.

    .venv/bin/python -I -m experiments.phase2.lab.B3.topo_probe --census <scratch>/census.jsonl --out <lab>/topo_probe.json

Inputs: census.jsonl (census.py). No images are read. Candidate probes:
- val: the val split's non-simple items;
- val+train_all: plus every train non-simple patient (A1's proposal, 366 raw-rule items);
- val+train_half: plus half of the train non-simple patients, dealt by objdel.data.splits.deal_weighted on their
  non-simple item counts (seeded), so about half of train's non-simple items stay in SFT;
- synthetic: N independent generated items with exact topology (precision_probe.py generator), one per background.
Metrics per item (greedy output vs the voting raters' size-filtered signatures):
- nonsimple_hit (primary; A1's "rater-topology hit"): the prediction's signature is one of the non-simple rater
  signatures (synthetic: equals the ground truth);
- topo_any: the prediction's signature is some rater's (the entry-gate metric; a simple blob passes whenever a rater
  drew a simple shape);
- hole recall (items where a voting rater drew a hole: the prediction has a hole), hole precision (over probe and
  simple-control items: predictions with a hole on items where some rater drew one), part recall / part-count accuracy.
CIs: percentile cluster bootstrap over nodules (Field & Welsh 2007: the cluster bootstrap is consistent under both
the random-effect and the transformation model), B = 2,000, for simulated predictions with a nodule random effect on
the logit scale (latent ICC rho = s^2 / (s^2 + pi^2 / 3)); half-width = (hi - lo) / 2, median over R simulated
datasets. The empirical ICC of the raters' non-simple share within nodules is reported as a plausible rho.
"""
import argparse
import json
from collections import Counter, defaultdict

import numpy as np

from objdel.data.splits import deal_weighted

from . import prereg

SIMPLE = [0]
RHOS = (0.0, 0.15, 0.3, 0.6, 0.9)
RATES = (0.3, 0.5, 0.7)


def load(path):
    with open(path) as fh:
        return [json.loads(line) for line in fh]


def anova_icc(values, groups):
    """One-way ANOVA ICC(1) with unequal group sizes (Donner 1986 n0)."""
    by = defaultdict(list)
    for v, g in zip(values, groups, strict=True):
        by[g].append(v)
    gs = [np.array(v, float) for v in by.values()]
    k, N = len(gs), sum(len(g) for g in gs)
    grand = np.concatenate(gs).mean()
    msb = sum(len(g) * (g.mean() - grand) ** 2 for g in gs) / (k - 1)
    msw = sum(((g - g.mean()) ** 2).sum() for g in gs) / (N - k)
    n0 = (N - sum(len(g) ** 2 for g in gs) / N) / (k - 1)
    return float((msb - msw) / (msb + (n0 - 1) * msw))


def simulate(cluster_sizes, rate, rho, rng):
    """0/1 predictions with a per-cluster logit random effect giving latent ICC rho."""
    s = np.sqrt(rho / (1 - rho) * np.pi ** 2 / 3) if rho > 0 else 0.0
    mu = np.log(rate / (1 - rate))
    out = []
    for m in cluster_sizes:
        p = 1 / (1 + np.exp(-(mu + s * rng.standard_normal())))
        out.append(rng.random(m) < p)
    # the logit-normal marginal equals `rate` only at rate = 0.5 (the pre-registered point); elsewhere it is pulled
    # toward 0.5, which slightly overstates the half-width
    return out


def boot_halfwidth(clusters, rng, B=2000):
    sums = np.array([c.sum() for c in clusters], float)
    ns = np.array([len(c) for c in clusters], float)
    idx = rng.integers(0, len(clusters), (B, len(clusters)))
    est = sums[idx].sum(1) / ns[idx].sum(1)
    lo, hi = np.percentile(est, [2.5, 97.5])
    return (hi - lo) / 2, sums.sum() / ns.sum()


def power(cluster_sizes, rng, R=200):
    out = {}
    for rho in RHOS:
        for rate in RATES:
            hw = [boot_halfwidth(simulate(cluster_sizes, rate, rho, rng), rng)[0] for _ in range(R)]
            out[f"rho={rho}|p={rate}"] = round(float(np.median(hw)), 4)
    return out


def paired_power(probe, rng, R=200, p1=0.5, p0_parts=0.5, rho=0.6):
    """A0 vs A0': A0 never emits holes (hit 0 on hole-only items), hits multi-part items at p0_parts; A0' hits every
    non-simple item at p1. Shared nodule effect. Half-width of the paired difference (A0' - A0) by cluster bootstrap."""
    by = defaultdict(list)
    for r in probe:
        hole = any(sum(r["sigs"][i]) > 0 for i in r["voters"])
        part = any(len(r["sigs"][i]) > 1 for i in r["voters"])
        by[r["nodule"]].append((hole, part))
    s = np.sqrt(rho / (1 - rho) * np.pi ** 2 / 3)
    lg = lambda p: np.log(p / (1 - p))  # noqa: E731
    hws, diffs = [], []
    for _ in range(R):
        d_sum, n = [], []
        for items in by.values():
            u = s * rng.standard_normal()
            a1 = rng.random(len(items)) < 1 / (1 + np.exp(-(lg(p1) + u)))
            a0 = np.array([part and not hole for hole, part in items]) & (
                rng.random(len(items)) < 1 / (1 + np.exp(-(lg(p0_parts) + u))))
            d_sum.append(float(a1.sum() - a0.sum()))
            n.append(len(items))
        d_sum, n = np.array(d_sum), np.array(n, float)
        idx = rng.integers(0, len(n), (2000, len(n)))
        est = d_sum[idx].sum(1) / n[idx].sum(1)
        lo, hi = np.percentile(est, [2.5, 97.5])
        hws.append((hi - lo) / 2)
        diffs.append(d_sum.sum() / n.sum())
    return {"halfwidth_median": round(float(np.median(hws)), 4), "diff_mean": round(float(np.mean(diffs)), 4),
            "rho": rho, "p1": p1, "p0_parts": p0_parts}


def describe(probe, train_all, held_patients):
    tr_ns = [r for r in train_all if r["nonsimple"]]
    kept_tr = [r for r in train_all if r["patient"] not in held_patients]
    ns_t = lambda rs: sum(tuple(r["sigs"][i]) != (0,) for r in rs for i in r["voters"])  # noqa: E731
    return {
        "items": len(probe), "nodules": len({r["nodule"] for r in probe}),
        "patients": len({r["patient"] for r in probe}),
        "unanimous": sum(r["unanimous_nonsimple"] for r in probe),
        "majority_nonsimple": sum(r["majority_nonsimple"] for r in probe),
        "hole_items": sum(any(sum(r["sigs"][i]) > 0 for i in r["voters"]) for r in probe),
        "multipart_items": sum(any(len(r["sigs"][i]) > 1 for i in r["voters"]) for r in probe),
        "cluster_size_max": max(Counter(r["nodule"] for r in probe).values()),
        "train_patients_removed": len(held_patients),
        "train_items_removed": len(train_all) - len(kept_tr),
        "train_targets_removed": sum(len(r["voters"]) for r in train_all) - sum(len(r["voters"]) for r in kept_tr),
        "train_nonsimple_items_kept": f"{sum(r['nonsimple'] for r in kept_tr)} / {len(tr_ns)}",
        "train_nonsimple_targets_kept_share": round(ns_t(kept_tr) / ns_t(train_all), 4),
        "train_unanimous_kept": sum(r["unanimous_nonsimple"] for r in kept_tr),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--probe-list", default=None, help="write the chosen real probe's item ids + patients here")
    ap.add_argument("--R", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    rows = load(a.census)
    tr = [r for r in rows if r["split"] == "train"]
    va = [r for r in rows if r["split"] == "val"]
    ns = [r for r in rows if r["nonsimple"]]

    # plausible within-nodule correlation: the raters' non-simple share across the slices of non-simple nodules
    pool_nod = {r["nodule"] for r in ns}
    pool = [r for r in rows if r["nodule"] in pool_nod]
    share = [np.mean([r["sigs"][i] != SIMPLE for i in r["voters"]]) for r in pool]
    icc = {"rater_nonsimple_share_all_slices_of_nonsimple_nodules": anova_icc(share, [r["nodule"] for r in pool]),
           "rater_nonsimple_share_nonsimple_items": anova_icc(
               [np.mean([r["sigs"][i] != SIMPLE for i in r["voters"]]) for r in ns], [r["nodule"] for r in ns]),
           "majority_nonsimple_nonsimple_items": anova_icc([float(r["majority_nonsimple"]) for r in ns],
                                                           [r["nodule"] for r in ns])}

    tr_ns_pat = Counter(r["patient"] for r in tr if r["nonsimple"])
    half = deal_weighted(sorted(tr_ns_pat), dict(tr_ns_pat), (("probe", 0.5), ("sft", 0.5)), a.seed)
    held_half = {p for p, s in half.items() if s == "probe"}
    designs = {
        "val": ([r for r in va if r["nonsimple"]], set()),
        "val+train_half": ([r for r in va if r["nonsimple"]]
                           + [r for r in tr if r["nonsimple"] and r["patient"] in held_half], held_half),
        "val+train_all": ([r for r in va if r["nonsimple"]] + [r for r in tr if r["nonsimple"]], set(tr_ns_pat)),
    }
    out = {"icc": icc, "designs": {}, "prereg": {"half_width": prereg.P1_HALF_WIDTH,
                                                 "min_train_nonsimple_kept": prereg.P1_MIN_TRAIN_NONSIMPLE_KEPT}}
    for name, (probe, held) in designs.items():
        sizes = list(Counter(r["nodule"] for r in probe).values())
        d = describe(probe, tr, held)
        d["binomial_halfwidth_p0.5"] = round(1.96 * np.sqrt(0.25 / len(probe)), 4)
        d["cluster_bootstrap_halfwidth"] = power(sizes, rng, a.R)
        d["paired_A0_vs_A0prime"] = {f"rho={rho}": paired_power(probe, rng, a.R // 2, rho=rho) for rho in (0.3, 0.6)}
        hw = d["cluster_bootstrap_halfwidth"]
        d["verdict"] = {
            "ci_rho0.15": hw["rho=0.15|p=0.5"] <= prereg.P1_HALF_WIDTH,
            "ci_rho0.3": hw["rho=0.3|p=0.5"] <= prereg.P1_HALF_WIDTH,
            "ci_rho0.6": hw["rho=0.6|p=0.5"] <= prereg.P1_HALF_WIDTH,
            "leak": "0 by construction (held-out patients removed from SFT train)",
            "train_kept": d["train_nonsimple_targets_kept_share"] >= prereg.P1_MIN_TRAIN_NONSIMPLE_KEPT}
        out["designs"][name] = d
        print(name, json.dumps(d), flush=True)
    # synthetic: independent items, one per background patch and shape (clusters of size 1)
    for n in (160, 240, 320, 480):
        out["designs"][f"synthetic_nonsimple_{n}"] = {
            "items": n, "cluster_bootstrap_halfwidth": power([1] * n, rng, a.R),
            "binomial_halfwidth_p0.5": round(1.96 * np.sqrt(0.25 / n), 4)}
        print(n, out["designs"][f"synthetic_nonsimple_{n}"], flush=True)
    if a.probe_list:
        probe = designs["val+train_half"][0]
        with open(a.probe_list, "w") as fh:
            json.dump({"design": "val+train_half", "seed": a.seed,
                       "held_out_train_patients": sorted(held_half),
                       "items": sorted(r["item_id"] for r in probe)}, fh, indent=0)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
