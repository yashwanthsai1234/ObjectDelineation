"""Iteration 3 (R3): can any IoU_tc / band weighting make the minority-hole adversary lose to the
minority rater's own jittered mask (criterion (a), pure geometry: same topology, same DM term)?
Raw comparison of the two geometric terms per group, global band, 1-px honest jitter.

    ../.venv/bin/python scripts/hole_weights.py
"""
import sys
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

from adversaries import item_data, load_items  # noqa: E402


def main():
    items = [(ds, k, st) for ds, k, st in load_items() if ds == "lidc"]
    with Pool(3) as pool:
        data = pool.map(item_data, [(ds, k, st, "1px", j, {"topology_conditional_band": False}) for j, (ds, k, st) in enumerate(items)],
                        chunksize=8)
    rows = []
    for d in data:
        if d["minority_rater"] is None or "minority_hole" not in d["adv"]:
            continue
        adv = d["adv"]["minority_hole"]
        for t in d["pools"][d["minority_rater"]]:
            rows.append((adv[1] - t[1], adv[2] - t[2]))  # adversary - own jittered mask: (IoU_tc, F5)
    r = np.array(rows)
    print(f"pairs (adversary vs a jittered mask of the minority rater): {len(r)}")
    print(f"  IoU_tc: adversary lower in {np.mean(r[:, 0] < 0):.1%}; median diff {np.median(r[:, 0]):+.3f}")
    print(f"  F5:     adversary lower in {np.mean(r[:, 1] < 0):.1%}; median diff {np.median(r[:, 1]):+.3f}")
    for w_iou, w_band in ((0.3, 0.3), (0.45, 0.15), (0.6, 0.0)):
        print(f"  w_iou {w_iou}, w_band {w_band}: raw weighted sum lower in {np.mean(w_iou * r[:, 0] + w_band * r[:, 1] < 0):.1%}")


if __name__ == "__main__":
    main()
