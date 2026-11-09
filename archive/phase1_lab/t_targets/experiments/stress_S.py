"""Builder success vs crop side S on pinch-heavy random patches embedded in an S x S crop
(pixel = 1000/S grid units). Reports fail rate (no legal topology-preserving target found)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from targets.builder import LIDC_SCHEMA, BuildConfig, build_target  # noqa: E402

rng = np.random.default_rng(7)
patches = []
for _ in range(150):
    h, w = rng.integers(4, 16, 2)
    m = rng.random((h, w)) < rng.uniform(0.3, 0.7)
    patches.append(m)
cfgs = {"hyb3": BuildConfig("minseg_hybrid", {"tau": 3}), "lossless": BuildConfig("minseg_exact", {}, iou_min=1.0),
        "exact_rings": BuildConfig("exact_rings", {}, iou_min=1.0)}
print("| S | grid units / px | config | fail % | fallback % | bevels / mask |")
print("|---|---|---|---|---|---|")
for S in [250, 499, 500, 501, 640, 999, 1000, 1024]:
    for name, cfg in cfgs.items():
        fails = fb = nb = 0
        for m in patches:
            big = np.zeros((S, S), bool)
            big[S // 3:S // 3 + m.shape[0], S // 2:S // 2 + m.shape[1]] = m
            _, objs, info = build_target({"nodule": big}, LIDC_SCHEMA, cfg)
            fails += objs["nodule"].status == "fail"
            fb += objs["nodule"].status == "fallback"
            nb += info["bevels"]
        print(f"| {S} | {1000 / S:.2f} | {name} | {100 * fails / len(patches):.1f} | {100 * fb / len(patches):.1f} | "
              f"{nb / len(patches):.2f} |", flush=True)
