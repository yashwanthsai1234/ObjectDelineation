"""Iteration 3 (R3): hole adversaries with the frozen candidate (global F5 band, debiased DM-loo,
c_topo = 2, lambda = 0.5, reward_v3), on every item that has them (fixtures + full-LIDC train).

minority_hole = vote with the minority rater's holes (topology matches that rater).
holes_stamped = vote minus every rater's holes.
For minority_hole two conditions are reported separately:
  (a) below the minority rater's OWN jittered mask (same topology, so the DM term is identical:
      pure geometry) - the Frankenstein test;
  (b) below the median honest rollout - conflicts with distribution matching when the group
      under-samples that (legitimate) rater topology, because DM-loo then boosts every rollout
      that carries it, the minority rater's own mask included.

    ../.venv/bin/python scripts/hole_adversaries.py
"""
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from adversaries import G, N_GROUPS, advantages, ci, item_data, load_items  # noqa: E402


def main():
    rng = np.random.default_rng(0)
    items = [(ds, k, st) for ds, k, st in load_items() if ds == "lidc"]
    for band_kw, label in (({"topology_conditional_band": False}, "global band (frozen)"),
                           ({"topology_conditional_band": True}, "topology-conditional band")):
        for amp in ("1px", "0.5tau"):
            with Pool(3) as pool:
                data = pool.map(item_data, [(ds, k, st, amp, j, band_kw) for j, (ds, k, st) in enumerate(items)], chunksize=8)
            print(f"\n[{label}; honest = rater masks + {amp} jitter]")
            for name in ("minority_hole", "holes_stamped"):
                T, it_of = [], []
                for d in data:
                    if name not in d["adv"] or (name == "minority_hole" and d["minority_rater"] is None):
                        continue
                    for _ in range(N_GROUPS):
                        draws = rng.integers(0, d["k"], G - 1)
                        if name == "minority_hole":
                            draws[0] = d["minority_rater"]
                        T.append([d["adv"][name]] + [d["pools"][r][rng.integers(len(d["pools"][r]))] for r in draws])
                        it_of.append(d["key"])
                T = np.array(T, float)
                a = np.concatenate([advantages(T[s:s + 256], "v3", 2.0) for s in range(0, len(T), 256)])
                below_med = a[:, 0] < np.median(a[:, 1:], axis=1)
                out = [("< median honest", below_med)]
                if name == "minority_hole":
                    out = [("(a) < own rater's jittered mask", a[:, 0] < a[:, 1]), ("(b) < median honest", below_med),
                           ("(a) and (b)", (a[:, 0] < a[:, 1]) & below_med)]
                for lab, ok in out:
                    lo, hi = ci(ok, it_of, rng)
                    print(f"  {name:14s} items {len(set(it_of)):3d}  {lab:34s} {ok.mean():6.1%} [{lo:6.1%}, {hi:6.1%}]")


if __name__ == "__main__":
    main()
