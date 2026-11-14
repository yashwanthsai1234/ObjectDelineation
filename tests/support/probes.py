"""Synthetic inputs of the probe builders (objdel.data.probes): lung patches, diameters, hole / two-part records."""
import numpy as np


def lung_patches(n=6, side=125, seed=0):
    rng = np.random.default_rng(seed)
    groups = [f"LIDC-IDRI-{i:04d}" for i in range(n)]
    return [(g, rng.normal(-850, 40, (side, side)).astype(np.int16)) for g in groups]


def diameters(n=4000, seed=0):  # log-normal fitted to the LIDC train p10 / p50 (6.9 / 12.2 px)
    return np.exp(np.random.default_rng(seed).normal(np.log(12.2), 0.444, n))


# hole / two-part records with the real p50s (hole_census.json): part diam 24 px, wall 1.24 px, offset 0.46
HOLES = [
    {"part_diam": d, "area": a, "wall": w, "offset": o}
    for d, a, w, o in [
        (24, 6, 1.2, 0.45),
        (22, 8, 0.41, 0.5),
        (26, 5, 1.0, 0.4),
        (25, 15, 1.5, 0.4),
        (28, 22, 2.0, 0.3),
        (24, 12, 1.24, 0.46),
        (34, 40, 2.5, 0.3),
        (40, 55, 3.0, 0.25),
        (30, 35, 1.5, 0.35),
    ]
]
MULTIS = [{"diam": d, "part_ratio": r, "gap": g} for d, r, g in [(18, 0.57, 6.8), (22, 0.4, 3.0), (16, 0.8, 12)]]
