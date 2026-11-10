"""C2 analysis (PLAN §5.2, D20, M8): A0' pooled hole recall on the topology probe, stratified bootstrap, per-seed
checks and the pre-registered outcome (experiments.phase2.prereg).

    .venv/bin/python -I -m experiments.phase2.c2_analyze --seeds s1/report.json s2/report.json s3/report.json \
        [s4/report.json] --out c2.json

Per seed (A0' arm) the report.json of `select_ckpt --only topology`, read at probes.topology (probes.score_topology
plus select_ckpt's "ids"): {"bins": [[per-item hole hit, mean over K = 8 at T = 1] x 120] x 3 hole-area bins
(<= 10 / 11-30 / > 30 px), in probe record order; "false_nonsimple": share of simple control samples drawn
non-simple; "legal": legal share; "ids": the sampled probe items}. Every seed must have sampled the same items in
the same order (so the per-item means pair item by item; exit otherwise). Per-item statistic = mean over seeds; pooled recall =
train-share-weighted mean of the bins; CI = within-bin item bootstrap (B3 power_c2.pooled, B = 2,000, percentile).
Clauses (4) and (5) use the mean over seeds. Reported, not gated: per-seed CIs, the across-seed t-interval, the
seed-discordant flag (per-seed CIs on both sides of 0.50, §5.0) and seeds_both_sides (point estimates, §5.2).
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats

from objdel.eval.probes import stratified_ci
from objdel.eval.signal import DIGITS

from . import prereg as P

WEIGHTS = P.C2["bin_weights"]


def pooled_ci(bins, seed: int = 0, boot: int = 2000) -> tuple[float, float, float]:
    """(point, lo, hi) of the weighted pooled mean, bins resampled within themselves (probes.stratified_ci = B3
    power_c2.pooled), rounded to 12 decimals so an exact 0.30 / 0.50 meets its line (BK2-4)."""
    return stratified_ci(bins, WEIGHTS, np.random.default_rng(seed), boot)


def analyze(seeds: list[dict], seed: int = 0, boot: int = 2000) -> dict:
    n = len(seeds)
    per_seed = [pooled_ci(s["bins"], seed, boot) for s in seeds]
    bins = [np.mean([np.asarray(s["bins"][b], float) for s in seeds], 0) for b in range(len(WEIGHTS))]
    point, lo, hi = pooled_ci(bins, seed, boot)
    pts = [p[0] for p in per_seed]
    fns, legal = (float(np.mean([s[k] for s in seeds])) for k in ("false_nonsimple", "legal"))
    bin_points = [round(float(b.mean()), DIGITS) for b in bins]
    clauses = P.c2_clauses(point, lo, pts, bin_points, fns, legal)
    outcome, conditional = P.c2_outcome(point, clauses)
    line = P.C2["pass_line"]
    t = stats.t.ppf(0.975, n - 1) * np.std(pts, ddof=1) / np.sqrt(n)
    return {
        "seeds": n,
        "pooled": [point, lo, hi],
        "bins": bin_points,
        "per_seed": [list(p) for p in per_seed],
        "seed_t_interval": [float(np.mean(pts) - t), float(np.mean(pts) + t)],
        "false_nonsimple": fns,
        "legal": legal,
        "clauses": clauses,
        "seed_discordant": bool(any(p[1] > line for p in per_seed) and any(p[2] < line for p in per_seed)),
        "seeds_both_sides": bool(min(pts) < line <= max(pts)),
        "fourth_seed_eligible": P.fourth_seed_eligible(pts, point, {k: v for k, v in clauses.items() if k != "2"}),
        "outcome": outcome,
        "conditional": conditional,
    }


def load(paths) -> list[dict]:
    """Each seed's probes.topology from select_ckpt --only topology's report.json, items checked equal."""
    seeds = [json.loads(Path(p).read_text())["probes"]["topology"] for p in paths]
    if any(s["ids"] != seeds[0]["ids"] for s in seeds):
        raise SystemExit(f"c2_analyze: the seeds sampled different topology probe items: {list(paths)}")
    return seeds


def main():
    ap = argparse.ArgumentParser(description="C2: A0' pooled hole recall, per-seed checks and outcome.")
    ap.add_argument("--seeds", nargs="+", required=True, help="per-seed select_ckpt --only topology report.json "
                    "(3, or 4 with the fourth seed)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    res = analyze(load(a.seeds), a.seed, a.boot)
    Path(a.out).write_text(json.dumps(res, indent=1))
    print(json.dumps({k: res[k] for k in ("pooled", "clauses", "outcome", "conditional")}))


if __name__ == "__main__":
    main()
