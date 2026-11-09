"""Realistic GRPO rollouts for K3/K4/K10: 1-3 rings, 16-200 vertices, square crops 64-512 px.

Shapes are smooth blobs (low-frequency radial noise) like model outputs for organs /
nodules; ring 2/3 is a hole (ccw) or a second part. A rollout comes with a ground-truth
mask: the exact raster of a slightly perturbed reference shape.
"""
from __future__ import annotations

import numpy as np

import polycore as pc


def blob(rng, n, cx, cy, r, amp=0.25, ccw=False):
    k = np.arange(1, 5)
    a = rng.normal(0, amp / k)
    ph = rng.uniform(0, 2 * np.pi, k.size)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False) + rng.uniform(0, 2 * np.pi / n)
    rad = r * (1 + (a[:, None] * np.cos(k[:, None] * t + ph[:, None])).sum(0))
    rad = np.clip(rad, 0.3 * r, None)
    pts = np.c_[cx + rad * np.cos(t), cy + rad * np.sin(t)]
    pts = np.clip(np.rint(pts), 0, 1000).astype(np.int64)
    keep = np.r_[True, np.any(np.diff(pts, axis=0) != 0, axis=1)]
    pts = pts[keep]
    return pts[::-1] if ccw else pts


def make_rollouts(n=256, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        S = int(rng.integers(64, 513))
        nr = int(rng.integers(1, 4))
        cx, cy = rng.integers(350, 651, 2)
        R = float(rng.uniform(150, 330))
        nv = lambda: int(rng.integers(16, 201))  # noqa: E731
        rings = [blob(rng, nv(), cx, cy, R)]
        roles = [-1]
        for _ in range(nr - 1):
            if rng.random() < 0.5:   # hole
                rings.append(blob(rng, nv(), cx + rng.integers(-20, 21), cy + rng.integers(-20, 21), R * rng.uniform(.15, .3), ccw=True))
                roles.append(0)
            else:                    # second part elsewhere (may overlap: realistic model error)
                rings.append(blob(rng, nv(), rng.integers(100, 901), rng.integers(100, 901), R * rng.uniform(.1, .3)))
                roles.append(-1)
        xy, off = pc.pack(rings)
        ref = [blob(rng, 64, cx + rng.integers(-15, 16), cy + rng.integers(-15, 16), R)]
        gt = pc.rasterize(ref, S, S)["evenodd"].astype(np.uint8)
        out.append(dict(rings=[r.tolist() for r in rings], xy=xy, off=off, roles=np.array(roles, np.int32), S=S, gt=gt))
    return out


def batch_args(rolls):
    return ([r["xy"] for r in rolls], [r["off"] for r in rolls], [r["S"] for r in rolls], [r["S"] for r in rolls],
            [r["roles"] for r in rolls], [r["gt"] for r in rolls])


def stats_hash(d):
    import hashlib
    h = hashlib.sha256()
    for k in sorted(d):
        v = d[k]
        h.update(k.encode())
        h.update(np.ascontiguousarray(v).tobytes() if isinstance(v, np.ndarray) else repr(v).encode())
    return h.hexdigest()
