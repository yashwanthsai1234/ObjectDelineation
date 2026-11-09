"""Exact reference implementations ("oracles") for Phase 1 module tests.

Slow, simple and exact (Python integers only). Candidate implementations are
tested against these, never the other way round.

Conventions (the spec every candidate must match)
- A polygon is given on an integer grid G x G (G = 1000 by default) covering a crop
  of W x H pixels. Grid point (x, y) maps to pixel coordinates (x * W / G, y * H / G).
- Pixel (col c, row r) is the unit square [c, c+1) x [r, r+1); it is "inside" iff its
  centre (c + 1/2, r + 1/2) is inside.
- Inside-ness uses the half-open crossing rule: a ray from the centre towards +x;
  an edge (x0,y0)->(x1,y1) counts iff min(y0,y1) <= py < max(y0,y1) and the crossing
  x is strictly greater than px. Even-odd: odd number of crossings. Nonzero: the sum
  of +1 (edge going down, y increasing) / -1 (going up) is nonzero.
  Centres on an edge are thus resolved deterministically (top-left-like rule).
- A shape is a list of rings; each ring is a list of (x, y) integer tuples, implicitly
  closed. Fill rules apply to all rings of a shape together.
"""
from __future__ import annotations

import random

import numpy as np


def _scaled(rings, W, H, G):
    # Scale so pixel centres and vertices are integers: x' = x * 2W, centre' = (2c+1) * G.
    return [[(x * 2 * W, y * 2 * H) for x, y in ring] for ring in rings]


def winding_and_parity(rings, W, H, G=1000):
    """Return (winding, crossings) int arrays of shape (H, W) for the pixel centres."""
    srings = _scaled(rings, W, H, G)
    wind = np.zeros((H, W), np.int64)
    cross = np.zeros((H, W), np.int64)
    for r in range(H):
        py = (2 * r + 1) * G
        for ring in srings:
            n = len(ring)
            for i in range(n):
                (x0, y0), (x1, y1) = ring[i], ring[(i + 1) % n]
                if y0 == y1 or not (min(y0, y1) <= py < max(y0, y1)):
                    continue
                d = 1 if y1 > y0 else -1
                for c in range(W):
                    px = (2 * c + 1) * G
                    # crossing x = x0 + (py - y0) * (x1 - x0) / (y1 - y0) > px ?
                    num = (x0 - px) * (y1 - y0) + (py - y0) * (x1 - x0)
                    if (num > 0) if (y1 - y0) > 0 else (num < 0):
                        wind[r, c] += d
                        cross[r, c] += 1
    return wind, cross


def rasterize(rings, W, H, G=1000, rule="evenodd"):
    wind, cross = winding_and_parity(rings, W, H, G)
    return (cross % 2 == 1) if rule == "evenodd" else (wind != 0)


# ---------- exact segment predicates (integers) ----------

def orient(a, b, c):
    v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (v > 0) - (v < 0)


def on_segment(a, b, p):
    return orient(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def segments_relation(a, b, c, d):
    """'none' | 'cross' (proper interior crossing) | 'touch' (share a point, incl. endpoints/collinear overlap)."""
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    if o1 * o2 < 0 and o3 * o4 < 0:
        return "cross"
    if on_segment(a, b, c) or on_segment(a, b, d) or on_segment(c, d, a) or on_segment(c, d, b):
        return "touch"
    return "none"


def twice_area(ring):
    n = len(ring)
    return sum(ring[i][0] * ring[(i + 1) % n][1] - ring[(i + 1) % n][0] * ring[i][1] for i in range(n))


def ring_problems(ring):
    """Exact legality labels of a single ring (sorted list; [] = simple polygon with nonzero area).

    Consecutive duplicate vertices are removed first; a closing duplicate is allowed.
    Labels: fewer_than_3_points, zero_area, self_crossing, self_touching (incl. collinear overlap,
    a vertex on a non-adjacent edge, spikes that fold back on the previous edge).
    """
    pts = [p for i, p in enumerate(ring) if i == 0 or p != ring[i - 1]]
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts = pts[:-1]
    if len(set(pts)) < 3:
        return ["fewer_than_3_points"]
    probs = set()
    if twice_area(pts) == 0:
        probs.add("zero_area")
    n = len(pts)
    for i in range(n):
        a, b = pts[i], pts[(i + 1) % n]
        for j in range(i + 1, n):
            c, d = pts[j], pts[(j + 1) % n]
            adjacent = j == i + 1 or (i == 0 and j == n - 1)
            if adjacent:
                # adjacent edges share exactly one endpoint; anything more is a fold-back overlap
                shared = b if j == i + 1 else a
                other_ab = a if shared == b else b
                other_cd = d if shared == c else c
                if orient(a, b, other_cd) == 0 and (on_segment(a, b, other_cd) or on_segment(c, d, other_ab)):
                    probs.add("self_touching")
                continue
            rel = segments_relation(a, b, c, d)
            if rel == "cross":
                probs.add("self_crossing")
            elif rel == "touch":
                probs.add("self_touching")
    return sorted(probs)


# ---------- generators ----------

def random_simple_polygon(rng: random.Random, n=12, G=1000, cx=None, cy=None, rmax=None):
    """Star-shaped simple polygon: sorted angles, random radii (rejects degenerate)."""
    cx = rng.randint(200, 800) if cx is None else cx
    cy = rng.randint(200, 800) if cy is None else cy
    rmax = rng.randint(30, 190) if rmax is None else rmax
    while True:
        angles = sorted(rng.uniform(0, 2 * np.pi) for _ in range(n))
        pts = []
        for a in angles:
            r = rng.uniform(0.3 * rmax, rmax)
            p = (min(G, max(0, round(cx + r * np.cos(a)))), min(G, max(0, round(cy + r * np.sin(a)))))
            if not pts or p != pts[-1]:
                pts.append(p)
        if len(set(pts)) >= 3 and ring_problems(pts) == []:
            return pts


def random_polygon_soup(rng: random.Random, n=10, G=1000):
    """Arbitrary vertex order: usually self-intersecting."""
    return [(rng.randint(0, G), rng.randint(0, G)) for _ in range(n)]


def random_grid_polygon(rng: random.Random, n=10, step=50, G=1000):
    """Vertices on a coarse grid: many collinear / touching / degenerate cases."""
    return [(rng.randrange(0, G + 1, step), rng.randrange(0, G + 1, step)) for _ in range(n)]


def annulus(cx=500, cy=500, r_out=300, r_in=120, n=24, hole_ccw=True):
    out = [(round(cx + r_out * np.cos(2 * np.pi * k / n)), round(cy + r_out * np.sin(2 * np.pi * k / n))) for k in range(n)]
    inn = [(round(cx + r_in * np.cos(2 * np.pi * k / n)), round(cy + r_in * np.sin(2 * np.pi * k / n))) for k in range(n)]
    return [out, inn[::-1] if hole_ccw else inn]
