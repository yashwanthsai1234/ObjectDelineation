"""Iteration 3 (R3 guard design): area excess over the largest rater, in band-tolerance units.

x = max(0, area - max rater area) / (T5 * P_U). Union, hull, honest jittered rollouts and a
held-out largest rater (a new expert that happens to draw the largest mask). 300 random items.

    ../.venv/bin/python scripts/guard_probe.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from adversaries import hull  # noqa: E402
from rewardcore import fixtures, perturb  # noqa: E402
from rewardcore.reward import ExpertCache  # noqa: E402


def main():
    rng = np.random.default_rng(0)
    res = {k: [] for k in ("union", "hull", "honest 1px (largest rater)", "honest 0.5tau (largest rater)",
                           "honest 0.5tau (any rater)", "held-out largest rater (k>=3)")}
    items = fixtures.load_all()
    for j in rng.choice(len(items), 300, replace=False):
        ds, key, st = items[j]
        c = ExpertCache.build(st)
        amax = max(c.rater_areas)
        s = c.tol5 * c.perimeter

        def x(a, amax=amax, s=s):
            return max(0, a - amax) / s
        res["union"].append(x(c.n_union))
        res["hull"].append(x(hull(c.union).sum()))
        jm = int(np.argmax(c.rater_areas))
        a = max(1.0, 0.5 * c.tau)
        res["honest 1px (largest rater)"].append(np.mean([x(perturb.jitter(st[jm], 1.0, rng, 3.0).sum()) for _ in range(4)]))
        res["honest 0.5tau (largest rater)"].append(np.mean([x(perturb.jitter(st[jm], a, rng, max(3, c.tau)).sum()) for _ in range(4)]))
        res["honest 0.5tau (any rater)"].append(
            np.mean([x(perturb.jitter(st[rng.integers(len(st))], a, rng, max(3, c.tau)).sum()) for _ in range(4)]))
        if len(st) >= 3:
            rest = ExpertCache.build(np.delete(st, jm, axis=0))
            res["held-out largest rater (k>=3)"].append(max(0, c.rater_areas[jm] - max(rest.rater_areas)) / (rest.tol5 * rest.perimeter))
    for k, v in res.items():
        v = np.array(v)
        print(f"{k:32s} n={len(v):3d} median {np.median(v):.3f} p90 {np.percentile(v, 90):.3f} p99 {np.percentile(v, 99):.3f}")


if __name__ == "__main__":
    main()
