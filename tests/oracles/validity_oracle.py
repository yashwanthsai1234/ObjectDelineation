"""Slow, exact validity oracle for objdel.checker.polycore.validate, and the shape generators of its tests.

Pure Python integers, all-pairs loops, no shared code with polycore (tests only). Derived from the Phase 1 lab
(experiments/phase1/lab/common/oracle.py and lab/t_polycore/tests/helpers.py).

- `ref_validate(rings, roles)`: every polycore defect class, counted pair by pair;
- `ring_problems(ring)`: the lab's label spec of one ring (K2 compares against it), the labels of `ref_validate`;
- generators: the lab's ring families (star, soup, grid, lattice walks, degenerate mutations, ...).

Rings are lists of (x, y) int tuples on the G x G grid, implicitly closed. Outer rings are clockwise on
screen (positive `twice_area`), holes counter-clockwise.
"""
from __future__ import annotations

import math
import random

G = 1000
DEFECTS = ("fewer_than_3_points", "zero_area", "self_crossing", "vertex_pinch", "collinear_overlap", "spike",
           "ring_ring_crossing", "ring_ring_touching", "hole_outside_outer", "nested_holes", "nested_outers",
           "wrong_winding", "out_of_range")
SELF_TOUCHING = ("vertex_pinch", "collinear_overlap", "spike")  # union = ring_problems' "self_touching"


# ------------------------------------------------------------------ exact predicates
def orient(a, b, c) -> int:
    v = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return (v > 0) - (v < 0)


def on_segment(a, b, p) -> bool:
    return orient(a, b, p) == 0 and min(a[0], b[0]) <= p[0] <= max(a[0], b[0]) and \
        min(a[1], b[1]) <= p[1] <= max(a[1], b[1])


def segments_relation(a, b, c, d) -> str:
    """'none' | 'cross' (proper interior crossing) | 'touch' (any shared point, incl. collinear overlap)."""
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    if o1 * o2 < 0 and o3 * o4 < 0:
        return "cross"
    if on_segment(a, b, c) or on_segment(a, b, d) or on_segment(c, d, a) or on_segment(c, d, b):
        return "touch"
    return "none"


def twice_area(ring) -> int:
    n = len(ring)
    return sum(ring[i][0] * ring[(i + 1) % n][1] - ring[(i + 1) % n][0] * ring[i][1] for i in range(n))


def clean(ring) -> list:
    """Drop consecutive repeats and a closing duplicate."""
    pts = [tuple(p) for i, p in enumerate(ring) if i == 0 or tuple(p) != tuple(ring[i - 1])]
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts = pts[:-1]
    return pts


def _adjacent_fold(a, b, c, d, j_follows_i) -> bool:
    """Adjacent edges ab, cd share one vertex; True iff they fold back over each other (a spike)."""
    shared = b if j_follows_i else a
    other_ab = a if shared == b else b
    other_cd = d if shared == c else c
    return orient(a, b, other_cd) == 0 and (on_segment(a, b, other_cd) or on_segment(c, d, other_ab))


def ring_problems(ring) -> list[str]:
    """Sorted labels of one ring ([] = simple with nonzero area): fewer_than_3_points, zero_area,
    self_crossing, self_touching (shared point of non-adjacent edges, or a fold-back of adjacent ones)."""
    return labels(list(ref_validate([ring])[0][0].values()))


def labels(per_ring_row) -> list[str]:
    """Map one row of polycore's per-ring counts to ring_problems-style labels."""
    c = dict(zip(DEFECTS, per_ring_row, strict=True))
    if c["fewer_than_3_points"]:
        return ["fewer_than_3_points"]
    out = [name for name, hit in (("self_crossing", c["self_crossing"]),
                                  ("self_touching", any(c[k] for k in SELF_TOUCHING)),
                                  ("zero_area", c["zero_area"])) if hit]
    return sorted(out)


# ------------------------------------------------------------------ full reference of polycore.validate
def _locate2(ring, q) -> int:
    """q in doubled coordinates: -1 outside, 0 on the boundary, 1 inside (even-odd)."""
    n, inside = len(ring), False
    for i in range(n):
        a = (2 * ring[i][0], 2 * ring[i][1])
        b = (2 * ring[(i + 1) % n][0], 2 * ring[(i + 1) % n][1])
        if on_segment(a, b, q):
            return 0
        if (a[1] <= q[1] < b[1]) or (b[1] <= q[1] < a[1]):
            num = (a[0] - q[0]) * (b[1] - a[1]) + (q[1] - a[1]) * (b[0] - a[0])
            if (num > 0) if b[1] > a[1] else (num < 0):
                inside = not inside
    return 1 if inside else -1


def _witnesses(ring):
    """Every vertex and edge midpoint, in doubled coordinates."""
    n = len(ring)
    for i in range(n):
        a, b = ring[i], ring[(i + 1) % n]
        yield (2 * a[0], 2 * a[1])
        yield (a[0] + b[0], a[1] + b[1])


def ref_validate(rings, roles=None, G=G):
    """All-pairs reference of polycore.validate -> (per-ring dict of counts, pair_cross, pair_touch)."""
    R = len(rings)
    cnt = [dict.fromkeys(DEFECTS, 0) for _ in range(R)]
    cl = [clean(r) for r in rings]
    fewer = [len(set(c)) < 3 for c in cl]
    pcross = [[0] * R for _ in range(R)]
    ptouch = [[0] * R for _ in range(R)]
    for r, ring in enumerate(rings):
        cnt[r]["out_of_range"] = sum(1 for x, y in ring if not (0 <= x <= G and 0 <= y <= G))
        if fewer[r]:
            cnt[r]["fewer_than_3_points"] = 1
            continue
        pts, n = cl[r], len(cl[r])
        if twice_area(pts) == 0:
            cnt[r]["zero_area"] = 1
        for i in range(n):
            a, b = pts[i], pts[(i + 1) % n]
            for j in range(i + 1, n):
                c, d = pts[j], pts[(j + 1) % n]
                if j == i + 1 or (i == 0 and j == n - 1):
                    cnt[r]["spike"] += _adjacent_fold(a, b, c, d, j == i + 1)
                    continue
                rel = segments_relation(a, b, c, d)
                if rel == "cross":
                    cnt[r]["self_crossing"] += 1
                elif rel == "touch":
                    overlap = orient(a, b, c) == 0 and orient(a, b, d) == 0 and _overlap_positive(a, b, c, d)
                    cnt[r]["collinear_overlap" if overlap else "vertex_pinch"] += 1
    for r in range(R):
        for q in range(r + 1, R):
            if fewer[r] or fewer[q]:
                continue
            A, B = cl[r], cl[q]
            for i in range(len(A)):
                for j in range(len(B)):
                    rel = segments_relation(A[i], A[(i + 1) % len(A)], B[j], B[(j + 1) % len(B)])
                    m = pcross if rel == "cross" else ptouch if rel == "touch" else None
                    if m is not None:
                        m[r][q] += 1
                        m[q][r] += 1
    hole = [roles is not None and roles[r] >= 0 for r in range(R)]
    for r in range(R):
        cnt[r]["ring_ring_crossing"] = sum(pcross[r])
        cnt[r]["ring_ring_touching"] = sum(ptouch[r])
        if fewer[r]:
            continue
        a2 = twice_area(cl[r])
        cnt[r]["wrong_winding"] = int((not hole[r] and a2 < 0) or (hole[r] and a2 > 0))
        if hole[r] and not fewer[roles[r]]:
            cnt[r]["hole_outside_outer"] = int(any(_locate2(cl[roles[r]], w) == -1 for w in _witnesses(cl[r])))

    def region(o, q):  # outer o with its holes, even-odd
        loc = _locate2(cl[o], q)
        if loc == 0:
            return 0
        inside = loc > 0
        for h in range(R):
            if hole[h] and roles[h] == o and not fewer[h]:
                lh = _locate2(cl[h], q)
                if lh == 0:
                    return 0
                inside ^= lh > 0
        return 1 if inside else -1

    for a in range(R):
        for b in range(a + 1, R):
            if fewer[a] or fewer[b]:
                continue
            if hole[a] and hole[b] and roles[a] == roles[b]:
                key, loc = "nested_holes", (lambda k, q: _locate2(cl[k], q))
            elif not hole[a] and not hole[b]:
                key, loc = "nested_outers", region
            else:
                continue
            if any(loc(a, w) == 1 for w in _witnesses(cl[b])) or any(loc(b, w) == 1 for w in _witnesses(cl[a])):
                cnt[a][key] += 1
                cnt[b][key] += 1
    return cnt, pcross, ptouch


def _overlap_positive(a, b, c, d) -> bool:
    k = 0 if a[0] != b[0] else 1
    return min(max(a[k], b[k]), max(c[k], d[k])) > max(min(a[k], b[k]), min(c[k], d[k]))


# ------------------------------------------------------------------ generators (lab families)
SIZES = (7, 37, 64, 125, 250, 333, 512, 1000)


def random_simple_polygon(rng: random.Random, n=12, cx=None, cy=None, rmax=None):
    """Star-shaped simple polygon: sorted angles, random radii (rejects degenerate draws)."""
    cx = rng.randint(200, 800) if cx is None else cx
    cy = rng.randint(200, 800) if cy is None else cy
    rmax = rng.randint(30, 190) if rmax is None else rmax
    while True:
        pts = []
        for a in sorted(rng.uniform(0, 2 * math.pi) for _ in range(n)):
            r = rng.uniform(0.3 * rmax, rmax)
            p = (min(G, max(0, round(cx + r * math.cos(a)))), min(G, max(0, round(cy + r * math.sin(a)))))
            if not pts or p != pts[-1]:
                pts.append(p)
        if len(set(pts)) >= 3 and ring_problems(pts) == []:
            return pts


def random_polygon_soup(rng, n=10):
    """Arbitrary vertex order: usually self-intersecting."""
    return [(rng.randint(0, G), rng.randint(0, G)) for _ in range(n)]


def random_grid_polygon(rng, n=10, step=50):
    """Vertices on a coarse grid: many collinear, touching and degenerate cases."""
    return [(rng.randrange(0, G + 1, step), rng.randrange(0, G + 1, step)) for _ in range(n)]


def annulus(cx, cy, r_out, r_in, n, hole_ccw=True):
    def circle(r):
        return [(round(cx + r * math.cos(2 * math.pi * k / n)), round(cy + r * math.sin(2 * math.pi * k / n)))
                for k in range(n)]
    inner = circle(r_in)
    return [circle(r_out), inner[::-1] if hole_ccw else inner]


def star(rng, n, cx, cy, rmax):
    pts = []
    for a in sorted(rng.uniform(0, 2 * math.pi) for _ in range(n)):
        r = rng.uniform(0.3 * rmax, rmax)
        p = (round(cx + r * math.cos(a)), round(cy + r * math.sin(a)))
        if not pts or p != pts[-1]:
            pts.append(p)
    return pts


def rect_staircase(rng):
    cx, cy, m, R = rng.randint(200, 800), rng.randint(200, 800), rng.randint(2, 12), rng.randint(20, 200)
    ang = sorted(rng.uniform(0, 2 * math.pi) for _ in range(m))
    Q = [(round(cx + R * rng.uniform(.5, 1) * math.cos(a)), round(cy + R * rng.uniform(.5, 1) * math.sin(a)))
         for a in ang]
    out = []
    for k in range(m):
        out += [Q[k], (Q[(k + 1) % m][0], Q[k][1])]
    return out


def half_pixel_walk(rng, W, H, n):
    """Vertices on the half-pixel lattice (pixel corners and centres, rounded to the grid)."""
    kx, ky = rng.randint(0, 2 * W), rng.randint(0, 2 * H)
    pts = []
    for _ in range(n):
        kx = min(2 * W, max(0, kx + rng.randint(-3, 3)))
        ky = min(2 * H, max(0, ky + rng.randint(-3, 3)))
        pts.append((round(kx * G / (2 * W)), round(ky * G / (2 * H))))
    return pts


def tiny(rng, W, H):
    s = rng.uniform(1, 3)
    bx, by = s * G / W, s * G / H
    x0, y0 = rng.uniform(0, G - bx), rng.uniform(0, G - by)
    return [(round(x0 + rng.uniform(0, bx)), round(y0 + rng.uniform(0, by))) for _ in range(rng.randint(3, 6))]


def maybe_rev(rng, ring):
    return ring[::-1] if rng.random() < 0.5 else ring


def mutate(rng, ring):
    """Inject degeneracies: repeated vertices, collinear midpoints, spikes, closing duplicates, pinches."""
    r = list(ring)
    for _ in range(rng.randint(1, 4)):
        if not r:
            break
        k, op = rng.randrange(len(r)), rng.randrange(6)
        if op == 0:
            r.insert(k, r[k])
        elif op == 1:
            a, b = r[k], r[(k + 1) % len(r)]
            if (a[0] + b[0]) % 2 == 0 and (a[1] + b[1]) % 2 == 0:
                r.insert(k + 1, ((a[0] + b[0]) // 2, (a[1] + b[1]) // 2))
        elif op == 2:
            a = r[k]
            r[k + 1:k + 1] = [(a[0] + rng.randint(-60, 60), a[1] + rng.randint(-60, 60)), a]
        elif op == 3:
            r.append(r[0])
        elif op == 4:
            r.insert(k, r[rng.randrange(len(r))])
        else:
            r[k] = (r[k][0] + rng.choice([-1, 0, 1]), r[k][1] + rng.choice([-1, 0, 1]))
    return r


K2_FAMILIES = ("star", "soup", "grid", "smallgrid", "degenerate", "rect", "halfpixel")


def k2_ring(i, seed=777):
    """Deterministic K2 case i -> (family, ring)."""
    rng = random.Random(seed * 1_000_003 + i)
    fam = K2_FAMILIES[i % len(K2_FAMILIES)]
    if fam == "star":
        r = maybe_rev(rng, random_simple_polygon(rng, rng.randint(3, 40)))
    elif fam == "soup":
        r = random_polygon_soup(rng, rng.randint(1, 20))
    elif fam == "grid":
        r = random_grid_polygon(rng, rng.randint(1, 14), step=rng.choice([50, 100, 125, 200, 250, 500]))
    elif fam == "smallgrid":
        k = rng.randint(1, 4)
        r = [(rng.randint(0, k), rng.randint(0, k)) for _ in range(rng.randint(1, 10))]
    elif fam == "degenerate":
        r = mutate(rng, rng.choice([lambda: random_simple_polygon(rng, rng.randint(3, 15)),
                                    lambda: random_grid_polygon(rng, rng.randint(3, 8), 250),
                                    lambda: rect_staircase(rng)])())
    elif fam == "rect":
        r = rect_staircase(rng)
    else:
        r = half_pixel_walk(rng, 125, 125, rng.randint(3, 14))
    return fam, r


def multi_ring(rng):
    """1-2 outers with 0-3 holes each, then optional reversal, mutation and coarse snapping -> (rings, roles)."""
    cx, cy, R = rng.randint(250, 750), rng.randint(250, 750), rng.randint(40, 250)
    rings, roles = [], []
    for _ in range(rng.randint(1, 2)):
        o = len(rings)
        ox, oy = cx + rng.randint(-200, 200), cy + rng.randint(-200, 200)
        rings.append(star(rng, rng.randint(3, 16), ox, oy, R))
        roles.append(-1)
        for _ in range(rng.randint(0, 3)):
            rings.append(star(rng, rng.randint(3, 10), ox + rng.randint(-R, R) // 2, oy + rng.randint(-R, R) // 2,
                              R * rng.uniform(.05, .6))[::-1])
            roles.append(o)
    if rng.random() < .3:
        k = rng.randrange(len(rings))
        rings[k] = rings[k][::-1]
    if rng.random() < .3:
        k = rng.randrange(len(rings))
        rings[k] = mutate(rng, rings[k])
    if rng.random() < .2:
        rings = [[(x // 50 * 50, y // 50 * 50) for x, y in r] for r in rings]
    return rings, roles
