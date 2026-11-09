"""Slow, exact raster oracle for the checker's masks (objdel._gpucheck), and the K1 shape families (tests only).

Derived from the Phase 1 lab (archive/phase1_lab/common/oracle.py, archive/phase1_lab/t_polycore/tests/helpers.py).
The spec, literally: pixel (c, r) is inside iff its centre is inside; a ray from the centre towards +x counts
edge (x0,y0)->(x1,y1) iff min(y0,y1) <= py < max(y0,y1) and the crossing x is strictly greater than px; the
winding adds +1 for edges going down the screen, -1 going up. Even-odd = odd crossings, nonzero = winding != 0.
Everything is scaled by (2W, 2H) so that vertices and centres ((2c+1) G, (2r+1) G) are integers.

- `winding_and_parity`: pure Python, pixel by pixel (the spec itself);
- `winding_np`: a literal numpy vectorisation of the same formula (exact below 2^30, asserted), proven
  equal to the pure oracle by a test, used for the 10k-shape K1 run.
"""
from __future__ import annotations

import random

import numpy as np
from oracles.validity_oracle import (
    SIZES,
    annulus,
    half_pixel_walk,
    maybe_rev,
    mutate,
    random_grid_polygon,
    random_polygon_soup,
    random_simple_polygon,
    rect_staircase,
    star,
    tiny,
)

G = 1000


def winding_and_parity(rings, W, H, G=G):
    """(winding, crossings) int64 (H, W) at the pixel centres, pure Python integers."""
    scaled = [[(x * 2 * W, y * 2 * H) for x, y in ring] for ring in rings]
    wind = np.zeros((H, W), np.int64)
    cross = np.zeros((H, W), np.int64)
    for r in range(H):
        py = (2 * r + 1) * G
        for ring in scaled:
            n = len(ring)
            for i in range(n):
                (x0, y0), (x1, y1) = ring[i], ring[(i + 1) % n]
                if y0 == y1 or not (min(y0, y1) <= py < max(y0, y1)):
                    continue
                for c in range(W):
                    px = (2 * c + 1) * G
                    num = (x0 - px) * (y1 - y0) + (py - y0) * (x1 - x0)  # (crossing x - px) * (y1 - y0)
                    if (num > 0) if y1 > y0 else (num < 0):
                        wind[r, c] += 1 if y1 > y0 else -1
                        cross[r, c] += 1
    return wind, cross


def winding_np(rings, W, H, G=G):
    """Literal numpy vectorisation of winding_and_parity."""
    wind = np.zeros((H, W), np.int64)
    cross = np.zeros((H, W), np.int64)
    px = (2 * np.arange(W, dtype=np.int64) + 1) * G
    pys = (2 * np.arange(H, dtype=np.int64) + 1) * G
    for ring in rings:
        n = len(ring)
        s = [(x * 2 * W, y * 2 * H) for x, y in ring]
        for i in range(n):
            (x0, y0), (x1, y1) = s[i], s[(i + 1) % n]
            assert max(abs(x0), abs(y0), abs(x1), abs(y1), 2 * W * G, 2 * H * G) < 2**30
            if y0 == y1:
                continue
            rows = np.nonzero((min(y0, y1) <= pys) & (pys < max(y0, y1)))[0]
            if rows.size == 0:
                continue
            num = (x0 - px[None, :]) * (y1 - y0) + (pys[rows][:, None] - y0) * (x1 - x0)
            hit = (num > 0) if y1 > y0 else (num < 0)
            wind[rows] += (1 if y1 > y0 else -1) * hit
            cross[rows] += hit
    return wind, cross


def oracle_winding(rings, W, H, G=G, pure=False):
    """The vectorised oracle when exact (all scaled values < 2^30), else the pure one."""
    m = max((max(abs(x), abs(y)) for r in rings for x, y in r), default=0)
    if pure or not (m * 2 * max(W, H) < 2**30 and 2 * max(W, H) * G < 2**30):
        return winding_and_parity(rings, W, H, G)
    return winding_np(rings, W, H, G)


def ring_ranges(rings, W, H, G=G, pure=False):
    """Per ring (min(0, min winding), max(0, max winding)) of the ring alone."""
    out = []
    for ring in rings:
        w, _ = oracle_winding([ring], W, H, G, pure)
        out.append((min(0, int(w.min())), max(0, int(w.max()))))
    return out


def on_boundary_centres(rings, W, H, G=G):
    """Bool (H, W): the pixel centre lies exactly on some edge (exact integer test)."""
    B = np.zeros((H, W), bool)
    X, Y = np.meshgrid((2 * np.arange(W, dtype=np.int64) + 1) * G, (2 * np.arange(H, dtype=np.int64) + 1) * G)
    for ring in rings:
        n = len(ring)
        s = [(x * 2 * W, y * 2 * H) for x, y in ring]
        for i in range(n):
            (x0, y0), (x1, y1) = s[i], s[(i + 1) % n]
            cr = (x1 - x0) * (Y - y0) - (y1 - y0) * (X - x0)
            B |= (cr == 0) & (min(x0, x1) <= X) & (X <= max(x0, x1)) & (min(y0, y1) <= Y) & (Y <= max(y0, y1))
    return B


K1_FAMILIES = ("star", "soup", "grid", "halfpixel", "holes", "multipart", "tiny", "degenerate", "rect")


def k1_case(i, seed=12345):
    """Deterministic K1 case i -> (family, rings, W, H); family = K1_FAMILIES[i % 9]."""
    rng = random.Random(seed * 1_000_003 + i)
    fam = K1_FAMILIES[i % len(K1_FAMILIES)]
    W, H = rng.choice(SIZES), rng.choice(SIZES)
    if fam == "star":
        rings = [maybe_rev(rng, random_simple_polygon(rng, rng.randint(3, 40)))]
    elif fam == "soup":
        rings = [random_polygon_soup(rng, rng.randint(3, 25))]
    elif fam == "grid":
        rings = [random_grid_polygon(rng, rng.randint(3, 16), step=rng.choice([50, 100, 125, 200, 250]))]
    elif fam == "halfpixel":
        rings = [half_pixel_walk(rng, W, H, rng.randint(3, 14)) for _ in range(rng.randint(1, 2))]
    elif fam == "holes":
        cx, cy, R = rng.randint(300, 700), rng.randint(300, 700), rng.randint(60, 290)
        if rng.random() < 0.3:
            rings = annulus(cx, cy, R, int(R * rng.uniform(.2, .8)), rng.randint(4, 40), hole_ccw=rng.random() < .5)
        else:
            rings = [star(rng, rng.randint(5, 40), cx, cy, R)] + [
                maybe_rev(rng, star(rng, rng.randint(3, 20), cx + rng.randint(-5, 5), cy + rng.randint(-5, 5),
                                    R * rng.uniform(.05, .29))) for _ in range(rng.randint(1, 2))]
    elif fam == "multipart":
        rings = [maybe_rev(rng, star(rng, rng.randint(3, 30), rng.randint(150, 850), rng.randint(150, 850),
                                     rng.randint(20, 300))) for _ in range(rng.randint(2, 4))]
    elif fam == "tiny":
        rings = [tiny(rng, W, H) for _ in range(rng.randint(1, 3))]
    elif fam == "degenerate":
        base = rng.choice([lambda: random_simple_polygon(rng, rng.randint(3, 20)),
                           lambda: random_grid_polygon(rng, rng.randint(3, 10), 100),
                           lambda: random_polygon_soup(rng, rng.randint(3, 10))])()
        rings = [mutate(rng, base)]
        if rng.random() < 0.3:
            rings.append([(rng.randint(0, G), rng.randint(0, G))] * rng.randint(1, 2))  # 1-2 vertex ring
        if rng.random() < 0.3:
            dx, dy = rng.randint(-400, 400), rng.randint(-400, 400)
            rings = [[(x + dx, y + dy) for x, y in r] for r in rings]  # out of range
    else:
        rings = [rect_staircase(rng)]
    return fam, rings, W, H
