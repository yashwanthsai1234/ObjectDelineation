"""Shared test helpers: shape generators, a vectorised exact oracle, a Python
reference for the fine-grained validity classes."""
from __future__ import annotations

import math
import random

import numpy as np

import oracle as O

SIZES = (7, 37, 64, 125, 250, 333, 512, 1000)
G = 1000


# ------------------------------------------------------------------ oracle (vectorised)
def winding_np(rings, W, H, G=G):
    """Literal numpy vectorisation of oracle.winding_and_parity (same integer formula).

    Exact while every scaled quantity stays < 2**30 (asserted); use the pure oracle otherwise.
    """
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
            py = pys[rows][:, None]
            num = (x0 - px[None, :]) * (y1 - y0) + (py - y0) * (x1 - x0)
            hit = (num > 0) if (y1 - y0) > 0 else (num < 0)
            wind[rows] += (1 if y1 > y0 else -1) * hit
            cross[rows] += hit
    return wind, cross


def fits_np(rings, W, H, G=G):
    m = max((max(abs(x), abs(y)) for r in rings for x, y in r), default=0)
    return m * 2 * max(W, H) < 2**30 and 2 * max(W, H) * G < 2**30


def oracle_winding(rings, W, H, G=G, pure=False):
    if pure or not fits_np(rings, W, H, G):
        return O.winding_and_parity(rings, W, H, G)
    return winding_np(rings, W, H, G)


def on_boundary_centres(rings, W, H, G=G):
    """Bool (H, W): pixel centre lies exactly on some edge (exact integer test)."""
    B = np.zeros((H, W), bool)
    px = (2 * np.arange(W, dtype=np.int64) + 1) * G
    pys = (2 * np.arange(H, dtype=np.int64) + 1) * G
    X, Y = np.meshgrid(px, pys)
    for ring in rings:
        n = len(ring)
        s = [(x * 2 * W, y * 2 * H) for x, y in ring]
        for i in range(n):
            (x0, y0), (x1, y1) = s[i], s[(i + 1) % n]
            cr = (x1 - x0) * (Y - y0) - (y1 - y0) * (X - x0)
            B |= (cr == 0) & (np.minimum(x0, x1) <= X) & (X <= np.maximum(x0, x1)) & \
                 (np.minimum(y0, y1) <= Y) & (Y <= np.maximum(y0, y1))
    return B


# ------------------------------------------------------------------ generators
def star(rng, n, cx, cy, rmax, G=G):
    ang = sorted(rng.uniform(0, 2 * math.pi) for _ in range(n))
    pts = []
    for a in ang:
        r = rng.uniform(0.3 * rmax, rmax)
        p = (round(cx + r * math.cos(a)), round(cy + r * math.sin(a)))
        if not pts or p != pts[-1]:
            pts.append(p)
    return pts


def rect_staircase(rng):
    cx, cy = rng.randint(200, 800), rng.randint(200, 800)
    m = rng.randint(2, 12)
    ang = sorted(rng.uniform(0, 2 * math.pi) for _ in range(m))
    R = rng.randint(20, 200)
    Q = [(round(cx + R * rng.uniform(.5, 1) * math.cos(a)), round(cy + R * rng.uniform(.5, 1) * math.sin(a))) for a in ang]
    out = []
    for k in range(m):
        out.append(Q[k])
        out.append((Q[(k + 1) % m][0], Q[k][1]))
    return out


def half_pixel_walk(rng, W, H, n, G=G):
    """Vertices on the half-pixel lattice (pixel corners and centres, rounded to the grid)."""
    kx, ky = rng.randint(0, 2 * W), rng.randint(0, 2 * H)
    pts = []
    for _ in range(n):
        kx = min(2 * W, max(0, kx + rng.randint(-3, 3)))
        ky = min(2 * H, max(0, ky + rng.randint(-3, 3)))
        pts.append((round(kx * G / (2 * W)), round(ky * G / (2 * H))))
    return pts


def tiny(rng, W, H, G=G):
    s = rng.uniform(1, 3)
    bx, by = s * G / W, s * G / H
    x0, y0 = rng.uniform(0, G - bx), rng.uniform(0, G - by)
    return [(round(x0 + rng.uniform(0, bx)), round(y0 + rng.uniform(0, by))) for _ in range(rng.randint(3, 6))]


def maybe_rev(rng, ring):
    return ring[::-1] if rng.random() < 0.5 else ring


def mutate(rng, ring):
    """Inject degeneracies: repeated vertices, collinear midpoints, spikes, closing duplicates."""
    r = list(ring)
    for _ in range(rng.randint(1, 4)):
        if not r:
            break
        k = rng.randrange(len(r))
        op = rng.randrange(6)
        if op == 0:
            r.insert(k, r[k])                                      # consecutive repeat
        elif op == 1:
            a, b = r[k], r[(k + 1) % len(r)]
            if (a[0] + b[0]) % 2 == 0 and (a[1] + b[1]) % 2 == 0:
                r.insert(k + 1, ((a[0] + b[0]) // 2, (a[1] + b[1]) // 2))  # redundant collinear vertex
        elif op == 2:
            a = r[k]
            r[k + 1:k + 1] = [(a[0] + rng.randint(-60, 60), a[1] + rng.randint(-60, 60)), a]  # spike out and back
        elif op == 3:
            r.append(r[0])                                         # closing duplicate
        elif op == 4:
            j = rng.randrange(len(r))
            r.insert(k, r[j])                                      # non-consecutive repeat (pinch)
        else:
            r[k] = (r[k][0] + rng.choice([-1, 0, 1]), r[k][1] + rng.choice([-1, 0, 1]))
    return r


K1_FAMILIES = ("star", "soup", "grid", "halfpixel", "holes", "multipart", "tiny", "degenerate", "rect")


def k1_case(i, seed=12345):
    """Deterministic K1 case i -> (family, rings, W, H)."""
    rng = random.Random(seed * 1_000_003 + i)
    fam = K1_FAMILIES[i % len(K1_FAMILIES)]
    W, H = rng.choice(SIZES), rng.choice(SIZES)
    if fam == "star":
        rings = [maybe_rev(rng, O.random_simple_polygon(rng, rng.randint(3, 40)))]
    elif fam == "soup":
        rings = [O.random_polygon_soup(rng, rng.randint(3, 25))]
    elif fam == "grid":
        rings = [O.random_grid_polygon(rng, rng.randint(3, 16), step=rng.choice([50, 100, 125, 200, 250]))]
    elif fam == "halfpixel":
        rings = [half_pixel_walk(rng, W, H, rng.randint(3, 14)) for _ in range(rng.randint(1, 2))]
    elif fam == "holes":
        cx, cy, R = rng.randint(300, 700), rng.randint(300, 700), rng.randint(60, 290)
        if rng.random() < 0.3:
            rings = O.annulus(cx, cy, R, int(R * rng.uniform(.2, .8)), rng.randint(4, 40), hole_ccw=rng.random() < .5)
        else:
            outer = star(rng, rng.randint(5, 40), cx, cy, R)
            rings = [outer] + [maybe_rev(rng, star(rng, rng.randint(3, 20), cx + rng.randint(-5, 5), cy + rng.randint(-5, 5), R * rng.uniform(.05, .29)))
                               for _ in range(rng.randint(1, 2))]
    elif fam == "multipart":
        rings = [maybe_rev(rng, star(rng, rng.randint(3, 30), rng.randint(150, 850), rng.randint(150, 850), rng.randint(20, 300)))
                 for _ in range(rng.randint(2, 4))]
    elif fam == "tiny":
        rings = [tiny(rng, W, H) for _ in range(rng.randint(1, 3))]
    elif fam == "degenerate":
        base = rng.choice([lambda: O.random_simple_polygon(rng, rng.randint(3, 20)),
                           lambda: O.random_grid_polygon(rng, rng.randint(3, 10), 100),
                           lambda: O.random_polygon_soup(rng, rng.randint(3, 10))])()
        rings = [mutate(rng, base)]
        if rng.random() < 0.3:
            rings.append([(rng.randint(0, G), rng.randint(0, G))] * rng.randint(1, 2))   # 1-2 vertex ring
        if rng.random() < 0.3:
            dx, dy = rng.randint(-400, 400), rng.randint(-400, 400)
            rings = [[(x + dx, y + dy) for x, y in r] for r in rings]                   # out of range
    else:
        rings = [rect_staircase(rng)]
    return fam, rings, W, H


K2_FAMILIES = ("star", "soup", "grid", "smallgrid", "degenerate", "rect", "halfpixel")


def k2_ring(i, seed=777):
    rng = random.Random(seed * 1_000_003 + i)
    fam = K2_FAMILIES[i % len(K2_FAMILIES)]
    if fam == "star":
        r = maybe_rev(rng, O.random_simple_polygon(rng, rng.randint(3, 40)))
    elif fam == "soup":
        r = O.random_polygon_soup(rng, rng.randint(1, 20))
    elif fam == "grid":
        r = O.random_grid_polygon(rng, rng.randint(1, 14), step=rng.choice([50, 100, 125, 200, 250, 500]))
    elif fam == "smallgrid":
        k = rng.randint(1, 4)
        r = [(rng.randint(0, k), rng.randint(0, k)) for _ in range(rng.randint(1, 10))]
    elif fam == "degenerate":
        r = mutate(rng, rng.choice([lambda: O.random_simple_polygon(rng, rng.randint(3, 15)),
                                    lambda: O.random_grid_polygon(rng, rng.randint(3, 8), 250),
                                    lambda: rect_staircase(rng)])())
    elif fam == "rect":
        r = rect_staircase(rng)
    else:
        r = half_pixel_walk(rng, 125, 125, rng.randint(3, 14))
    return fam, r


# ------------------------------------------------------------------ python reference for fine classes
def _clean(ring):
    pts = [p for i, p in enumerate(ring) if i == 0 or p != ring[i - 1]]
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts = pts[:-1]
    return pts


def _overlap_len_pos(a, b, c, d):
    k = 0 if a[0] != b[0] else 1
    lo = max(min(a[k], b[k]), min(c[k], d[k]))
    hi = min(max(a[k], b[k]), max(c[k], d[k]))
    return hi > lo


def _touch_kind(a, b, c, d):
    if O.orient(a, b, c) == 0 and O.orient(a, b, d) == 0 and _overlap_len_pos(a, b, c, d):
        return "collinear_overlap"
    return "vertex_pinch"


def _locate2(ring, q):
    """q in doubled coords; -1 out, 0 on boundary, 1 in (even-odd)."""
    n = len(ring)
    inside = False
    for i in range(n):
        a = (2 * ring[i][0], 2 * ring[i][1])
        b = (2 * ring[(i + 1) % n][0], 2 * ring[(i + 1) % n][1])
        if O.on_segment(a, b, q):
            return 0
        if (a[1] <= q[1] < b[1]) or (b[1] <= q[1] < a[1]):
            num = (a[0] - q[0]) * (b[1] - a[1]) + (q[1] - a[1]) * (b[0] - a[0])
            if (num > 0) if b[1] > a[1] else (num < 0):
                inside = not inside
    return 1 if inside else -1


def _witnesses(ring):
    n = len(ring)
    for i in range(n):
        a, b = ring[i], ring[(i + 1) % n]
        yield (2 * a[0], 2 * a[1])
        yield (a[0] + b[0], a[1] + b[1])


def ref_validate(rings, roles=None, G=G):
    """Slow independent Python reference of polycore.validate -> (ring counts dict list, pair_cross, pair_touch)."""
    import polycore as pc
    R = len(rings)
    names = pc.DEFECTS
    cnt = [dict.fromkeys(names, 0) for _ in range(R)]
    clean = [_clean(r) for r in rings]
    fewer = [len(set(c)) < 3 for c in clean]
    pc_ = np.zeros((R, R), int)
    pt_ = np.zeros((R, R), int)
    for r, ring in enumerate(rings):
        cnt[r]["out_of_range"] = sum(1 for x, y in ring if not (0 <= x <= G and 0 <= y <= G))
        if fewer[r]:
            cnt[r]["fewer_than_3_points"] = 1
            continue
        pts = clean[r]
        n = len(pts)
        if O.twice_area(pts) == 0:
            cnt[r]["zero_area"] = 1
        for i in range(n):
            a, b = pts[i], pts[(i + 1) % n]
            for j in range(i + 1, n):
                c, d = pts[j], pts[(j + 1) % n]
                if j == i + 1 or (i == 0 and j == n - 1):
                    shared = b if j == i + 1 else a
                    oab = a if shared == b else b
                    ocd = d if shared == c else c
                    if O.orient(a, b, ocd) == 0 and (O.on_segment(a, b, ocd) or O.on_segment(c, d, oab)):
                        cnt[r]["spike"] += 1
                    continue
                rel = O.segments_relation(a, b, c, d)
                if rel == "cross":
                    cnt[r]["self_crossing"] += 1
                elif rel == "touch":
                    cnt[r][_touch_kind(a, b, c, d)] += 1
    for r in range(R):
        for q in range(r + 1, R):
            if fewer[r] or fewer[q]:
                continue
            A, B = clean[r], clean[q]
            for i in range(len(A)):
                a, b = A[i], A[(i + 1) % len(A)]
                for j in range(len(B)):
                    c, d = B[j], B[(j + 1) % len(B)]
                    rel = O.segments_relation(a, b, c, d)
                    if rel == "cross":
                        pc_[r, q] += 1
                        pc_[q, r] += 1
                    elif rel == "touch":
                        pt_[r, q] += 1
                        pt_[q, r] += 1
    for r in range(R):
        cnt[r]["ring_ring_crossing"] = int(pc_[r].sum())
        cnt[r]["ring_ring_touching"] = int(pt_[r].sum())
    is_hole = [roles is not None and roles[r] >= 0 for r in range(R)]
    for r in range(R):
        if fewer[r]:
            continue
        a2 = O.twice_area(clean[r])
        if (not is_hole[r] and a2 < 0) or (is_hole[r] and a2 > 0):
            cnt[r]["wrong_winding"] = 1
        if is_hole[r] and not fewer[roles[r]]:
            if any(_locate2(clean[roles[r]], w) == -1 for w in _witnesses(clean[r])):
                cnt[r]["hole_outside_outer"] = 1

    def region(o, q):
        loc = _locate2(clean[o], q)
        if loc == 0:
            return 0
        inside = loc > 0
        for h in range(R):
            if is_hole[h] and roles[h] == o and not fewer[h]:
                lh = _locate2(clean[h], q)
                if lh == 0:
                    return 0
                inside ^= lh > 0
        return 1 if inside else -1

    for a in range(R):
        for b in range(a + 1, R):
            if fewer[a] or fewer[b]:
                continue
            if is_hole[a] and is_hole[b] and roles[a] == roles[b]:
                if any(_locate2(clean[a], w) == 1 for w in _witnesses(clean[b])) or \
                   any(_locate2(clean[b], w) == 1 for w in _witnesses(clean[a])):
                    cnt[a]["nested_holes"] += 1
                    cnt[b]["nested_holes"] += 1
            elif not is_hole[a] and not is_hole[b]:
                if any(region(a, w) == 1 for w in _witnesses(clean[b])) or \
                   any(region(b, w) == 1 for w in _witnesses(clean[a])):
                    cnt[a]["nested_outers"] += 1
                    cnt[b]["nested_outers"] += 1
    return cnt, pc_, pt_


def counts_as_dicts(res):
    import polycore as pc
    return [dict(zip(pc.DEFECTS, map(int, row))) for row in res["ring"]]
