"""Ring simplifiers. Kept vertices are always lattice points of the exact crack-edge ring (pixel
corners), so after mapping to the 0..G grid every vertex is the rounded image of a pixel corner.

Input of every simplifier: a closed unit-step lattice path (n, 2) int64 (all lattice points of an exact
ring from objdel.checker.polygonize; foreground on the right of every step, y down), its step
directions (n,) in {0 E, 1 S, 2 W, 3 N}, and the crop frame (W, H, G). Output: sorted kept indices.

All admissibility tests run in the *target frame*: X = gx * 2W, Y = gy * 2H, where (gx, gy) is the
rounded grid image of a lattice point. Pixel centres are then the integers ((2c+1) G, (2r+1) G),
exactly as in common/oracle.py, so every test is exact int64 arithmetic and already includes the
grid rounding.

Methods
- minseg (potrace polygon stage, corner-anchored, no vertex adjustment): optimal cycle by
  lexicographic (number of segments, sum of squared deviations) over admissible segments i->j:
  the sub-path i..j does not use all four directions, and
    linf  : the ray from i through j meets the L-inf box of half-size t around every skipped
            lattice point (potrace's straightness test; t_px pixels, potrace = 1/2, or tau grid units),
    exact : the segment's line strictly separates the inside and outside pixel centres adjacent to
            every unit step of the sub-path (lossless: no pixel centre changes side),
    hybrid: exact OR linf(tau).
- douglas_peucker(eps grid units), visvalingam(area grid units^2) on grid coordinates.
- corners_only: the exact ring (fallback, lossless for crops < G px).
"""
from __future__ import annotations

import heapq
import math

import numpy as np
from numba import njit

DIRV = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]], np.int64)  # E S W N (y down)


def lattice_path(corners):
    """Exact ring corners (k, 2) -> all unit-step lattice points (n, 2) and step directions (n,)."""
    c = np.asarray(corners, dtype=np.int64)
    nxt = np.roll(c, -1, axis=0)
    d = nxt - c
    lens = np.abs(d).sum(axis=1)
    assert np.all((d[:, 0] == 0) | (d[:, 1] == 0)) and np.all(lens > 0), "ring not rectilinear"
    unit = np.sign(d)
    idx = np.repeat(np.arange(len(c)), lens)
    off = np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens)
    pts = c[idx] + unit[idx] * off[:, None]
    u = unit[idx]
    dirs = np.where(u[:, 0] == 1, 0, np.where(u[:, 1] == 1, 1, np.where(u[:, 0] == -1, 2, 3)))
    return pts.astype(np.int64), dirs.astype(np.int64)


def to_grid(px, W, H, G=1000):
    """Pixel-corner coordinates -> grid integers, round half up (exact integer arithmetic)."""
    px = np.asarray(px, np.int64).reshape(-1, 2)
    return np.stack([(2 * px[:, 0] * G + W) // (2 * W), (2 * px[:, 1] * G + H) // (2 * H)], axis=1)


def frame(path, dirs, W, H, G=1000):
    """Target-frame vertices P and adjacent inside/outside pixel centres CI, CO (all int64)."""
    g = to_grid(path, W, H, G)
    P = g * np.array([2 * W, 2 * H], np.int64)
    d = DIRV[dirs]
    rgt = np.stack([-d[:, 1], d[:, 0]], axis=1)  # screen-right normal (foreground side)
    CI = G * (2 * path + d + rgt)
    CO = G * (2 * path + d - rgt)
    return P, CI.astype(np.int64), CO.astype(np.int64)


# ------------------------------------------------------------------ admissibility (numba)

@njit(cache=True)
def _wrap_to(a, ref):
    two = 2.0 * math.pi
    return a - two * math.floor((a - ref + math.pi) / two)


@njit(cache=True)
def _narrow(lo, hi, have, ref, amin, amax):
    mid = 0.5 * (amin + amax)
    if not have:
        ref = mid
        have = True
    sh = _wrap_to(mid, ref) - mid
    return max(lo, amin + sh), min(hi, amax + sh), have, ref


@njit(cache=True)
def _cross(ax, ay, bx, by):
    return ax * by - ay * bx


@njit(cache=True)
def _in_arc(vx, vy, lx, ly, hx, hy):
    """v inside the closed arc lo..hi (arc narrower than pi), exact."""
    return _cross(lx, ly, vx, vy) >= 0 and _cross(vx, vy, hx, hy) >= 0


@njit(cache=True)
def _box_arc(dx, dy, tx, ty):
    """Extreme corner directions (lo, hi) of the box centred at (dx, dy) with half sides (tx, ty),
    seen from the origin (origin outside the box)."""
    cx = np.array([dx - tx, dx + tx, dx - tx, dx + tx])
    cy = np.array([dy - ty, dy - ty, dy + ty, dy + ty])
    lo = 0
    hi = 0
    for a in range(4):
        okl = True
        okh = True
        for c in range(4):
            if _cross(cx[a], cy[a], cx[c], cy[c]) < 0:
                okl = False
            if _cross(cx[c], cy[c], cx[a], cy[a]) < 0:
                okh = False
        if okl:
            lo = a
        if okh:
            hi = a
    return cx[lo], cy[lo], cx[hi], cy[hi]


@njit(cache=True)
def _admissible(P, dirs, CI, CO, tx, ty, use_linf, use_exact, maxspan):
    """adm[i, l]: segment P[i] -> P[(i+l)%n] admissible (l = 1..L); pen[i, l]: its penalty.

    linf uses potrace's ray semantics: the direction P[j]-P[i] must lie in the (exact, integer) cone of
    directions whose ray from P[i] meets every box around P[i+1..j-1]. exact checks every adjacent
    centre with integer cross products. Penalty = sum of squared distances (pixels) of the skipped
    lattice points to the segment's line, via prefix sums (O(1) per segment)."""
    n = P.shape[0]
    # prefix sums over the doubled path in pixel-ish float units
    Q = P.astype(np.float64) / 2000.0
    sx_ = np.zeros(2 * n + 1)
    sy_ = np.zeros(2 * n + 1)
    sxx = np.zeros(2 * n + 1)
    syy = np.zeros(2 * n + 1)
    sxy = np.zeros(2 * n + 1)
    for k in range(2 * n):
        x, y = Q[k % n, 0], Q[k % n, 1]
        sx_[k + 1] = sx_[k] + x
        sy_[k + 1] = sy_[k] + y
        sxx[k + 1] = sxx[k] + x * x
        syy[k + 1] = syy[k] + y * y
        sxy[k + 1] = sxy[k] + x * y
    reach = np.ones(n, np.int64)
    rows = []
    for i in range(n):
        cnt = np.zeros(4, np.int64)
        full = True
        lx = ly = hx = hy = 0
        alive_l = use_linf
        lo_e, hi_e, have_e, ref_e = -1e18, 1e18, False, 0.0
        alive_e = use_exact
        px, py = P[i, 0], P[i, 1]
        oks = np.zeros(min(n, maxspan) + 1, np.bool_)
        oks[1] = True
        l = 1
        while l < n - 1 and l < maxspan:
            j = (i + l) % n
            cnt[dirs[(i + l - 1) % n]] += 1
            if cnt[0] > 0 and cnt[1] > 0 and cnt[2] > 0 and cnt[3] > 0:
                break
            if not (alive_l or alive_e):
                break
            reach[i] = l
            vx, vy = P[j, 0] - px, P[j, 1] - py
            if l >= 2 and not (vx == 0 and vy == 0):
                ok = False
                if alive_l and (full or _in_arc(vx, vy, lx, ly, hx, hy)):
                    ok = True
                if not ok and alive_e:
                    ok = True
                    for m in range(0, l):
                        e = (i + m) % n
                        if vx * (CI[e, 1] - py) - vy * (CI[e, 0] - px) <= 0:
                            ok = False
                            break
                        if vx * (CO[e, 1] - py) - vy * (CO[e, 0] - px) >= 0:
                            ok = False
                            break
                oks[l] = ok
            # constraints that apply to candidates > l
            if alive_l and (abs(vx) > tx or abs(vy) > ty):
                ax_, ay_, bx_, by_ = _box_arc(vx, vy, tx, ty)
                if full:
                    lx, ly, hx, hy = ax_, ay_, bx_, by_
                    full = False
                else:
                    if _in_arc(ax_, ay_, lx, ly, hx, hy):
                        nlx, nly = ax_, ay_
                    elif _in_arc(lx, ly, ax_, ay_, bx_, by_):
                        nlx, nly = lx, ly
                    else:
                        alive_l = False
                        nlx, nly = 0, 0
                    if alive_l:
                        if _in_arc(bx_, by_, lx, ly, hx, hy):
                            nhx, nhy = bx_, by_
                        elif _in_arc(hx, hy, ax_, ay_, bx_, by_):
                            nhx, nhy = hx, hy
                        else:
                            alive_l = False
                            nhx, nhy = 0, 0
                    if alive_l:
                        lx, ly, hx, hy = nlx, nly, nhx, nhy
                        if _cross(lx, ly, hx, hy) < 0:
                            alive_l = False
            if alive_e:
                th = math.atan2(CI[j, 1] - py, CI[j, 0] - px)
                lo_e, hi_e, have_e, ref_e = _narrow(lo_e, hi_e, have_e, ref_e, th - math.pi, th)
                th = math.atan2(CO[j, 1] - py, CO[j, 0] - px)
                lo_e, hi_e, have_e, ref_e = _narrow(lo_e, hi_e, have_e, ref_e, th, th + math.pi)
                if lo_e > hi_e + 1e-9:
                    alive_e = False
            l += 1
        rows.append(oks[: reach[i] + 1].copy())
    L = 1
    for i in range(n):
        L = max(L, reach[i])
    adm = np.zeros((n, L + 1), np.bool_)
    pen = np.zeros((n, L + 1), np.float64)
    for i in range(n):
        r = rows[i]
        ax, ay = Q[i, 0], Q[i, 1]
        for l in range(1, r.shape[0]):
            if not r[l]:
                continue
            adm[i, l] = True
            if l < 2:
                continue
            j = (i + l) % n
            dx, dy = Q[j, 0] - ax, Q[j, 1] - ay
            a0, a1 = i + 1, i + l  # skipped points i+1 .. i+l-1 in the doubled prefix arrays
            m = float(l - 1)
            Sxx = sxx[a1] - sxx[a0] - 2 * ax * (sx_[a1] - sx_[a0]) + m * ax * ax
            Syy = syy[a1] - syy[a0] - 2 * ay * (sy_[a1] - sy_[a0]) + m * ay * ay
            Sxy = (sxy[a1] - sxy[a0] - ax * (sy_[a1] - sy_[a0]) - ay * (sx_[a1] - sx_[a0]) + m * ax * ay)
            # sum (dx*y' - dy*x')^2 / |d|^2 with x' = x - ax, y' = y - ay
            num = dx * dx * Syy - 2 * dx * dy * Sxy + dy * dy * Sxx
            pen[i, l] = max(num, 0.0) / (dx * dx + dy * dy)
    return adm, pen


@njit(cache=True)
def _best_cycle(adm, pen):
    """Optimal cycle (min #segments, then min penalty) over the admissible graph. Returns kept indices."""
    n, L1 = adm.shape
    L = L1 - 1
    maxl = np.ones(n, np.int64)
    for i in range(n):
        for l in range(L, 0, -1):
            if adm[i, l]:
                maxl[i] = l
                break
    # choose the anchor point i* with the fewest straddling segment starts
    best_cnt = 1 << 60
    best_i = 0
    for istar in range(n):
        c = 0
        for back in range(0, L + 1):
            a = (istar - back) % n
            if maxl[a] > back:
                c += 1
        if c < best_cnt:
            best_cnt = c
            best_i = istar
    INF = 1 << 60
    best_count = INF
    best_pen = 1e300
    best_path = np.zeros(1, np.int64)
    cnt = np.empty(n + 1, np.int64)
    pv = np.empty(n + 1, np.float64)
    par = np.empty(n + 1, np.int64)
    for back in range(0, L + 1):
        a = (best_i - back) % n
        if maxl[a] <= back:
            continue
        for q in range(n + 1):
            cnt[q] = INF
            pv[q] = 1e300
            par[q] = -1
        cnt[0] = 0
        pv[0] = 0.0
        for q in range(n):
            if cnt[q] == INF:
                continue
            src = (a + q) % n
            for l in range(1, min(L, n - q) + 1):
                if not adm[src, l]:
                    continue
                r = q + l
                c2 = cnt[q] + 1
                p2 = pv[q] + pen[src, l]
                if c2 < cnt[r] or (c2 == cnt[r] and p2 < pv[r]):
                    cnt[r] = c2
                    pv[r] = p2
                    par[r] = q
        if cnt[n] < best_count or (cnt[n] == best_count and pv[n] < best_pen):
            best_count = cnt[n]
            best_pen = pv[n]
            m = cnt[n]
            out = np.empty(m, np.int64)
            r = n
            k = m - 1
            while r > 0:
                q = par[r]
                out[k] = (a + q) % n
                k -= 1
                r = q
            best_path = out
    return np.sort(best_path)


def apply_pins(adm, pins):
    """Forbid segments that skip over a pinned lattice index (the polygon must pass through pins)."""
    if pins is None or len(pins) == 0:
        return adm
    n, L1 = adm.shape
    p = np.unique(np.asarray(pins, np.int64) % n)
    i = np.arange(n)
    k = np.searchsorted(p, i, side="right")
    nxt = np.where(k < len(p), p[np.minimum(k, len(p) - 1)], p[0] + n)
    d = nxt - i  # steps to the next pin strictly after i
    adm = adm & (np.arange(L1)[None, :] <= d[:, None])
    adm[:, 1] = True
    return adm


def minseg(path, dirs, W, H, G=1000, mode="linf", t_px=None, tau=None, maxspan=100000, pins=None):
    """mode linf: tolerance t_px pixels (potrace: 0.5) or tau grid units; exact; hybrid (exact or tau).
    pins: lattice indices the polygon must keep (used by the local repair)."""
    n = len(path)
    if n <= 4:
        return np.arange(n)
    P, CI, CO = frame(path, dirs, W, H, G)
    if t_px is not None:
        tx = ty = int(round(t_px * 2 * G))
    elif tau is not None:
        tx, ty = int(round(tau * 2 * W)), int(round(tau * 2 * H))
    else:
        tx = ty = 0
    use_linf = mode in ("linf", "hybrid")
    use_exact = mode in ("exact", "hybrid")
    adm, pen = _admissible(P, dirs, CI, CO, tx, ty, use_linf, use_exact, maxspan)
    return _best_cycle(apply_pins(adm, pins), pen)


# ------------------------------------------------------------------ Douglas-Peucker / Visvalingam

@njit(cache=True)
def _dp_chain(path, i0, i1, n, eps2, keep):
    # iterative stack over (start, end) offsets on the cyclic chain
    stack = [(i0, i1)]
    while len(stack):
        a, b = stack.pop()
        if b - a < 2:
            continue
        ax, ay = path[a % n, 0], path[a % n, 1]
        bx, by = path[b % n, 0], path[b % n, 1]
        dx, dy = bx - ax, by - ay
        L2 = dx * dx + dy * dy
        best = -1.0
        bi = -1
        for k in range(a + 1, b):
            px, py = path[k % n, 0], path[k % n, 1]
            if L2 == 0:
                d2 = float((px - ax) ** 2 + (py - ay) ** 2)
            else:
                tt = ((px - ax) * dx + (py - ay) * dy) / L2
                if tt < 0:
                    d2 = float((px - ax) ** 2 + (py - ay) ** 2)
                elif tt > 1:
                    d2 = float((px - bx) ** 2 + (py - by) ** 2)
                else:
                    c = (dx * (py - ay) - dy * (px - ax))
                    d2 = c * c / L2
            if d2 > best:
                best = d2
                bi = k
        if best > eps2:
            keep[bi % n] = True
            stack.append((a, bi))
            stack.append((bi, b))


def douglas_peucker(path, dirs, W, H, G=1000, eps=5.0):
    """Cyclic DP on grid coordinates, eps in grid units."""
    path = to_grid(path, W, H, G)
    n = len(path)
    if n <= 4:
        return np.arange(n)
    d2 = ((path - path[0]) ** 2).sum(axis=1)
    far = int(np.argmax(d2))
    keep = np.zeros(n, bool)
    keep[0] = keep[far] = True
    _dp_chain(path, 0, far, n, float(eps) ** 2, keep)
    _dp_chain(path, far, n, n, float(eps) ** 2, keep)
    return np.nonzero(keep)[0]


def visvalingam(path, dirs, W, H, G=1000, area=25.0):
    """Visvalingam-Whyatt on grid coordinates, area in grid units^2."""
    path = to_grid(path, W, H, G)
    n = len(path)
    if n <= 4:
        return np.arange(n)
    P = path.astype(np.float64)
    prev = np.roll(np.arange(n), 1)
    nxt = np.roll(np.arange(n), -1)
    alive = np.ones(n, bool)

    def tri(i):
        a, b, c = P[prev[i]], P[i], P[nxt[i]]
        return abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2.0

    cur = np.array([tri(i) for i in range(n)])
    heap = [(cur[i], i) for i in range(n)]
    heapq.heapify(heap)
    alive_n = n
    while heap and alive_n > 3:
        a, i = heapq.heappop(heap)
        if not alive[i] or a != cur[i]:
            continue
        if a >= area:
            break
        alive[i] = False
        alive_n -= 1
        p, q = prev[i], nxt[i]
        nxt[p], prev[q] = q, p
        for j in (p, q):
            cur[j] = tri(j)
            heapq.heappush(heap, (cur[j], j))
    return np.nonzero(alive)[0]


def corners_only(path, dirs, *a, **k):
    """No simplification: the exact ring's corners (fallback)."""
    return np.nonzero(dirs != np.roll(dirs, 1))[0]


METHODS = {
    "exact_rings": corners_only,
    "dp": douglas_peucker,
    "vw": visvalingam,
    "minseg_linf": lambda p, d, W, H, G=1000, **k: minseg(p, d, W, H, G, "linf", **k),
    "minseg_exact": lambda p, d, W, H, G=1000, **k: minseg(p, d, W, H, G, "exact", **k),
    "minseg_hybrid": lambda p, d, W, H, G=1000, **k: minseg(p, d, W, H, G, "hybrid", **k),
}
