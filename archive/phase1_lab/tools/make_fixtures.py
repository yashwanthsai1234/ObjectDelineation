"""Build small fixture sets of real expert masks for Phase 1 module tests.

    python -I make_fixtures.py <repo_root> <out_dir>
Writes riga_experts.npz and lidc_experts.npz: keys '<item>__<structure>' -> uint8 (k, h, w) stacks.
"""
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from PIL import Image

repo, out = Path(sys.argv[1]), Path(sys.argv[2])
sys.path.insert(0, str(repo))  # -I drops cwd; add the repo explicitly
from experiments.phase0.lidc_census import nodule_clusters  # noqa: E402

rng = random.Random(0)
TRUE = re.compile(r"^(.+)_image(\d+)-(\d)\.png$")
groups = defaultdict(dict)
for f in (repo / "data/riga_true").glob("*.png"):
    if m := TRUE.match(f.name):
        groups[(m.group(1), int(m.group(2)))][int(m.group(3))] = f
keys = sorted(k for k in groups if len(groups[k]) >= 2)
rng.shuffle(keys)
riga = {}
for sub, idx in keys[:200]:
    ex = [np.asarray(Image.open(f)) for _, f in sorted(groups[(sub, idx)].items())]
    disc = np.stack([e >= 100 for e in ex]); cup = np.stack([e >= 200 for e in ex])
    u = disc.any(0); ys, xs = np.nonzero(u)
    h, w = np.ptp(ys) + 1, np.ptp(xs) + 1
    y0, y1 = max(0, ys.min() - h // 2), min(u.shape[0], ys.max() + h // 2 + 1)
    x0, x1 = max(0, xs.min() - w // 2), min(u.shape[1], xs.max() + w // 2 + 1)
    name = f"{sub}_{idx}"
    riga[f"{name}__disc"] = disc[:, y0:y1, x0:x1].astype(np.uint8)
    riga[f"{name}__cup"] = cup[:, y0:y1, x0:x1].astype(np.uint8)
np.savez_compressed(out / "riga_experts.npz", **riga)
print("riga items", len(riga) // 2)

xmls = sorted((repo / "data/raw/lidc_xml").rglob("*.xml"))
rng.shuffle(xmls)
lidc = {}
for p in xmls:
    if len(lidc) >= 600:
        break
    try:
        got = nodule_clusters(p)
    except Exception:
        continue
    if got is None:
        continue
    series, sessions, zi, clusters, _ = got
    for ci, cl in enumerate(clusters):
        for k in sorted({k for sl in cl["per_reader"].values() for k in sl}):
            masks = [sl[k] for sl in cl["per_reader"].values() if k in sl]
            if len(masks) < 2:
                continue
            raw = np.any(masks, axis=0)
            if not raw.any():
                continue
            ys, xs = np.nonzero(raw)
            sl = (slice(max(0, ys.min() - 8), ys.max() + 9), slice(max(0, xs.min() - 8), xs.max() + 9))
            lidc[f"{str(series)[-10:]}_{ci}_{k}__nodule"] = np.stack([m[sl] for m in masks]).astype(np.uint8)
np.savez_compressed(out / "lidc_experts.npz", **lidc)
print("lidc slices", len(lidc))
