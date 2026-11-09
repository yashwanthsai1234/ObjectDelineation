"""Fast exact geometry kernels (numba, int64 only) matching ../common/oracle.py.

- raster(rings, W, H, G): pixel-centre, half-open crossing rule, even-odd and winding.
- legality(rings): exact ring_problems per ring + ring-ring touch/cross, with offending edge pairs.

Rings are lists of (k, 2) int64 arrays of grid vertices (implicitly closed).
"""
from __future__ import annotations

import numpy as np
from numba import njit


def pack(rings):
    """List of (k,2) arrays -> (xy (N,2) int64, starts (R+1,) int64)."""
    rings = [np.asarray(r, dtype=np.int64).reshape(-1, 2) for r in rings]
    starts = np.zeros(len(rings) + 1, np.int64)
    for i, r in enumerate(rings):
        starts[i + 1] = starts[i] + len(r)
    xy = np.concatenate(rings) if rings else np.zeros((0, 2), np.int64)
    return np.ascontiguousarray(xy, dtype=np.int64), starts


# ---------------------------------------------------------------- rasterizer

@njit(cache=True)
def _raster_kernel(xy, starts, W, H, G):
    cross = np.zeros((H, W + 1), np.int32)
    wind = np.zeros((H, W + 1), np.int32)
    R = starts.shape[0] - 1
    for r in range(R):
        s, e = starts[r], starts[r + 1]
        n = e - s
        for i in range(n):
            x0 = xy[s + i, 0] * 2 * W
            y0 = xy[s + i, 1] * 2 * H
            x1 = xy[s + (i + 1) % n, 0] * 2 * W
            y1 = xy[s + (i + 1) % n, 1] * 2 * H
            if y0 == y1:
                continue
            d = 1 if y1 > y0 else -1
            ylo = min(y0, y1)
            yhi = max(y0, y1)
            # rows with ylo <= (2r+1)G < yhi
            # smallest r: (2r+1)G >= ylo  -> r >= (ylo - G) / (2G)
            rlo = -((-(ylo - G)) // (2 * G))
            if rlo < 0:
                rlo = 0
            # largest r: (2r+1)G < yhi -> (2r+1)G <= yhi-1
            rhi = (yhi - 1 - G) // (2 * G)
            if rhi > H - 1:
                rhi = H - 1
            # normalise so dy > 0
            if y1 > y0:
                ax, ay, bx, by = x0, y0, x1, y1
            else:
                ax, ay, bx, by = x1, y1, x0, y0
            dy = by - ay
            dx = bx - ax
            D = G * dy
            for row in range(rlo, rhi + 1):
                py = (2 * row + 1) * G
                N = ax * dy + (py - ay) * dx  # crossing X = N / dy ; pixel c counts iff (2c+1)G*dy < N
                # count c >= 0 with (2c+1)*D < N  <=>  2c+1 <= ceil(N/D)-1
                cn = -((-N) // D)  # ceil(N / D)
                k = (cn - 2) // 2 + 1
                if k <= 0:
                    continue
                if k > W:
                    k = W
                cross[row, 0] += 1
                cross[row, k] -= 1
                wind[row, 0] += d
                wind[row, k] -= d
    return cross, wind


@njit(cache=True)
def _prefix(acc, W, odd):
    H = acc.shape[0]
    out = np.zeros((H, W), np.bool_)
    for r in range(H):
        s = 0
        for c in range(W):
            s += acc[r, c]
            out[r, c] = (s & 1) == 1 if odd else s != 0
    return out


def raster(rings, W, H, G=1000, rule="evenodd"):
    xy, starts = pack(rings)
    cross, wind = _raster_kernel(xy, starts, int(W), int(H), int(G))
    if rule == "evenodd":
        return _prefix(cross, int(W), True)
    return _prefix(wind, int(W), False)


# ---------------------------------------------------------------- legality

@njit(cache=True)
def _orient(ax, ay, bx, by, cx, cy):
    v = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    return (v > 0) - (v < 0)


@njit(cache=True)
def _onseg(ax, ay, bx, by, px, py):
    return (_orient(ax, ay, bx, by, px, py) == 0 and min(ax, bx) <= px <= max(ax, bx)
            and min(ay, by) <= py <= max(ay, by))


@njit(cache=True)
def _relation(ax, ay, bx, by, cx, cy, dx, dy):
    """0 none, 1 touch, 2 cross (oracle.segments_relation)."""
    o1 = _orient(ax, ay, bx, by, cx, cy)
    o2 = _orient(ax, ay, bx, by, dx, dy)
    o3 = _orient(cx, cy, dx, dy, ax, ay)
    o4 = _orient(cx, cy, dx, dy, bx, by)
    if o1 * o2 < 0 and o3 * o4 < 0:
        return 2
    if (_onseg(ax, ay, bx, by, cx, cy) or _onseg(ax, ay, bx, by, dx, dy)
            or _onseg(cx, cy, dx, dy, ax, ay) or _onseg(cx, cy, dx, dy, bx, by)):
        return 1
    return 0


@njit(cache=True)
def _legality_kernel(xy, starts, max_pairs):
    """Rings must already be free of consecutive duplicates and closing repeats.

    Returns per-ring flags (R, 4): [fewer_than_3, zero_area, self_crossing, self_touching],
    inter-ring flag counts [touch, cross], and up to max_pairs offending (edge_a, edge_b, kind) rows
    (edge index = global vertex index of the edge's start).
    """
    R = starts.shape[0] - 1
    N = xy.shape[0]
    flags = np.zeros((R, 4), np.int64)
    inter = np.zeros(2, np.int64)
    pairs = np.zeros((max_pairs, 3), np.int64)
    npairs = 0
    ring_of = np.zeros(N, np.int64)
    nxt = np.zeros(N, np.int64)
    for r in range(R):
        s, e = starts[r], starts[r + 1]
        n = e - s
        for i in range(s, e):
            ring_of[i] = r
            nxt[i] = i + 1 if i + 1 < e else s
        # distinct points
        if n < 3:
            flags[r, 0] = 1
        else:
            a2 = 0
            for i in range(s, e):
                j = nxt[i]
                a2 += xy[i, 0] * xy[j, 1] - xy[j, 0] * xy[i, 1]
            if a2 == 0:
                flags[r, 1] = 1
    # sweep over edges sorted by min x
    xmin = np.empty(N, np.int64)
    xmax = np.empty(N, np.int64)
    ymin = np.empty(N, np.int64)
    ymax = np.empty(N, np.int64)
    for i in range(N):
        j = nxt[i]
        xmin[i] = min(xy[i, 0], xy[j, 0])
        xmax[i] = max(xy[i, 0], xy[j, 0])
        ymin[i] = min(xy[i, 1], xy[j, 1])
        ymax[i] = max(xy[i, 1], xy[j, 1])
    order = np.argsort(xmin, kind="mergesort")
    for oi in range(N):
        i = order[oi]
        ax, ay = xy[i, 0], xy[i, 1]
        bx, by = xy[nxt[i], 0], xy[nxt[i], 1]
        for oj in range(oi + 1, N):
            j = order[oj]
            if xmin[j] > xmax[i]:
                break
            if ymin[j] > ymax[i] or ymax[j] < ymin[i]:
                continue
            cx, cy = xy[j, 0], xy[j, 1]
            dx, dy = xy[nxt[j], 0], xy[nxt[j], 1]
            ri, rj = ring_of[i], ring_of[j]
            kind = 0
            if ri == rj:
                n = starts[ri + 1] - starts[ri]
                adjacent = nxt[i] == j or nxt[j] == i
                if adjacent:
                    if n < 3:
                        continue
                    # oracle: shared endpoint; fold-back overlap if collinear and the other end on the other edge
                    if nxt[i] == j:  # a->b, b->d  shared b
                        if _orient(ax, ay, bx, by, dx, dy) == 0 and (
                                _onseg(ax, ay, bx, by, dx, dy) or _onseg(cx, cy, dx, dy, ax, ay)):
                            kind = 3
                    else:  # c->d(=a), a->b shared a
                        if _orient(ax, ay, bx, by, cx, cy) == 0 and (
                                _onseg(ax, ay, bx, by, cx, cy) or _onseg(cx, cy, dx, dy, bx, by)):
                            kind = 3
                    if kind == 3:
                        flags[ri, 3] = 1
                else:
                    rel = _relation(ax, ay, bx, by, cx, cy, dx, dy)
                    if rel == 2:
                        flags[ri, 2] = 1
                        kind = 2
                    elif rel == 1:
                        flags[ri, 3] = 1
                        kind = 1
            else:
                rel = _relation(ax, ay, bx, by, cx, cy, dx, dy)
                if rel == 2:
                    inter[1] += 1
                    kind = 2
                elif rel == 1:
                    inter[0] += 1
                    kind = 1
            if kind and npairs < max_pairs:
                pairs[npairs, 0] = i
                pairs[npairs, 1] = j
                pairs[npairs, 2] = kind
                npairs += 1
    return flags, inter, pairs[:npairs]


def clean_ring(ring):
    """Drop consecutive duplicates and a closing repeat (oracle convention)."""
    r = np.asarray(ring, dtype=np.int64).reshape(-1, 2)
    if len(r) == 0:
        return r
    keep = np.r_[True, np.any(r[1:] != r[:-1], axis=1)]
    r = r[keep]
    if len(r) > 1 and np.all(r[0] == r[-1]):
        r = r[:-1]
    return r


def ring_problems(ring):
    """Exact per-ring labels identical to oracle.ring_problems."""
    r = clean_ring(ring)
    if len({tuple(p) for p in r.tolist()}) < 3:
        return ["fewer_than_3_points"]
    xy, starts = pack([r])
    flags, _, _ = _legality_kernel(xy, starts, 1)
    out = []
    if flags[0, 1]:
        out.append("zero_area")
    if flags[0, 2]:
        out.append("self_crossing")
    if flags[0, 3]:
        out.append("self_touching")
    return sorted(out)


def legality(rings, max_pairs=4096):
    """Exact legality of a set of rings (one object). Returns dict with per-ring problems,
    inter-ring touch/cross counts and offending edge pairs as (ring, edge, ring, edge, kind)."""
    rings = [clean_ring(r) for r in rings]
    probs = []
    for r in rings:
        if len({tuple(p) for p in r.tolist()}) < 3:
            probs.append(["fewer_than_3_points"])
        else:
            probs.append(None)
    xy, starts = pack(rings)
    flags, inter, pairs = _legality_kernel(xy, starts, max_pairs)
    for k in range(len(rings)):
        if probs[k] is None:
            p = []
            if flags[k, 1]:
                p.append("zero_area")
            if flags[k, 2]:
                p.append("self_crossing")
            if flags[k, 3]:
                p.append("self_touching")
            probs[k] = sorted(p)
    ring_of = np.repeat(np.arange(len(rings)), np.diff(starts))
    loc = [(int(ring_of[a]), int(a - starts[ring_of[a]]), int(ring_of[b]), int(b - starts[ring_of[b]]), int(kd))
           for a, b, kd in pairs]
    ok = all(not p for p in probs) and inter.sum() == 0
    return {"ok": bool(ok), "ring_problems": probs, "ring_touch": int(inter[0]), "ring_cross": int(inter[1]),
            "pairs": loc}


@njit(cache=True)
def _point_in_ring(px2, py2, xy, s, e):
    """Even-odd test of point (px2/2, py2/2) given doubled coords against ring vertices (doubled inside)."""
    n = e - s
    inside = False
    for i in range(n):
        x0, y0 = 2 * xy[s + i, 0], 2 * xy[s + i, 1]
        x1, y1 = 2 * xy[s + (i + 1) % n, 0], 2 * xy[s + (i + 1) % n, 1]
        if (y0 <= py2 < y1) or (y1 <= py2 < y0):
            num = (x0 - px2) * (y1 - y0) + (py2 - y0) * (x1 - x0)
            if (num > 0) if (y1 - y0) > 0 else (num < 0):
                inside = not inside
    return inside


def point_in_ring(pt, ring):
    xy, starts = pack([ring])
    return bool(_point_in_ring(2 * int(pt[0]), 2 * int(pt[1]), xy, 0, len(xy)))


def twice_area(ring):
    r = np.asarray(ring, dtype=np.int64)
    return int(np.sum(r[:, 0] * np.roll(r[:, 1], -1) - np.roll(r[:, 0], -1) * r[:, 1]))
