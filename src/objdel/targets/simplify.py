"""Minimum-vertex simplification of one exact crack-edge ring (potrace polygon stage OR lossless straight segments).

Input: the ring's closed unit-step lattice path p_0..p_{n-1} (pixel corners, y down, foreground on the right of
every step). Output: the indices of the kept corners; every kept vertex is a pixel corner, so with a crop side
dividing 1000 it is an exact grid point (no vertex adjustment, unlike potrace).

A shortcut p_i -> p_j (2 <= j - i, cyclically) is admissible when the steps i..j-1 do not use all four directions
and either
- *lossless*: the shortcut's line has every pixel centre on the inner side of those steps strictly on its right and
  every outer centre strictly on its left, so no pixel centre changes side (a digital straight segment test), or
- *L-inf tau*: the ray from p_i through p_j meets the L-inf box of half-side tau grid units around every skipped
  point (potrace's straightness test, with tau instead of 1/2 px).
Both tests narrow an exact integer cone of feasible directions as j grows (potrace §2.2.1: "a constraint consists of
at most two inequalities and can be updated in constant time"), so all shortcuts cost O(n * longest shortcut).
They run in the target frame X = 2W·gx, Y = 2H·gy of the rounded grid points, where pixel centres are the odd
multiples of G: every predicate is an int64 cross product and already includes the grid rounding. Admissibility is a
heuristic; exactness comes from the builder's certification on the checker.

The polygon is the optimal cycle in the shortcut graph, lexicographic in (number of segments, penalty), as in
potrace §2.2.3, where the penalty is the sum of squared pixel distances of the skipped corners to the segment's line
(exact integer numerator from prefix sums, one IEEE division, so platform-independent). Potrace fixes vertex 0
("non-cyclic version" in its bestpolygon); here the cycle is optimal: every cycle has a vertex at the anchor a* or a
segment over it, so one shortest-path run per such start suffices.

References:
    P. Selinger, "Potrace: a polygon-based tracing algorithm", 2003, §2.2: https://potrace.sourceforge.net/potrace.pdf
    I. Debled-Rennesson, J.-P. Reveillès, "A linear algorithm for segmentation of digital curves", IJPRAI 9(4), 1995,
        https://doi.org/10.1142/S0218001495000249 (arithmetical DSS; DGtal: https://dgtal.org/doc/2.1/moduleArithDSSReco.html)
    J.-O. Lachaud, A. Vialard, F. de Vieilleville, "Fast, accurate and convergent tangent estimation on digital
        contours" / maximal DSS: https://arxiv.org/abs/0906.2716
    H. Imai, M. Iri, min-# approximation as a shortest path; survey and O(n^2) min-#: https://arxiv.org/abs/cs/0209007
    archive/phase1_lab/t_targets/targets/simplify.py (the lab kernel this replaces; same admissibility)
"""
import numpy as np
from numba import njit

STEPS = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]], np.int64)  # E S W N on screen (y down)


def lattice_path(corners) -> tuple[np.ndarray, np.ndarray]:
    """Rectilinear ring corners (k, 2) -> every unit-step lattice point (n, 2) and step directions (n,) in 0..3."""
    c = np.asarray(corners, np.int64)
    d = np.roll(c, -1, 0) - c
    length = np.abs(d).sum(1)
    idx = np.repeat(np.arange(len(c)), length)
    u = np.sign(d)[idx]
    off = np.arange(length.sum()) - np.repeat(np.cumsum(length) - length, length)
    dirs = np.select([u[:, 0] == 1, u[:, 1] == 1, u[:, 0] == -1], [0, 1, 2], 3)
    return c[idx] + u * off[:, None], dirs


def to_grid(p, W: int, H: int, G: int = 1000) -> np.ndarray:
    """Pixel corners -> 0..G grid, round half up (exact; the identity scaling when W and H divide G)."""
    p = np.asarray(p, np.int64).reshape(-1, 2)
    return np.stack([(2 * p[:, 0] * G + W) // (2 * W), (2 * p[:, 1] * G + H) // (2 * H)], 1)


@njit(cache=True)
def _cross(ax, ay, bx, by):
    return ax * by - ay * bx


@njit(cache=True)
def _within(vx, vy, lx, ly, hx, hy):
    """v in the closed arc lo..hi (narrower than pi)."""
    return _cross(lx, ly, vx, vy) >= 0 and _cross(vx, vy, hx, hy) >= 0


@njit(cache=True)
def _meet(c, ax, ay, bx, by):
    """Intersect the cone c = [alive, lx, ly, hx, hy, bounded] in place with the closed arc a..b (< pi)."""
    if c[5] == 0:
        c[1], c[2], c[3], c[4], c[5] = ax, ay, bx, by, 1
        return
    lx, ly, hx, hy = c[1], c[2], c[3], c[4]
    if _within(ax, ay, lx, ly, hx, hy):
        lx, ly = ax, ay
    elif not _within(lx, ly, ax, ay, bx, by):
        c[0] = 0
        return
    if _within(bx, by, c[1], c[2], hx, hy):
        hx, hy = bx, by
    elif not _within(hx, hy, ax, ay, bx, by):
        c[0] = 0
        return
    c[1], c[2], c[3], c[4] = lx, ly, hx, hy
    if _cross(lx, ly, hx, hy) < 0:
        c[0] = 0


@njit(cache=True)
def _meet_box(c, dx, dy, tx, ty):
    """Intersect c with the directions whose ray from the origin meets the box (dx +- tx, dy +- ty)."""
    xs = (dx - tx, dx + tx, dx - tx, dx + tx)
    ys = (dy - ty, dy - ty, dy + ty, dy + ty)
    lo = hi = 0
    for a in range(4):
        is_lo = is_hi = True
        for b in range(4):
            is_lo &= _cross(xs[a], ys[a], xs[b], ys[b]) >= 0
            is_hi &= _cross(xs[b], ys[b], xs[a], ys[a]) >= 0
        lo = a if is_lo else lo
        hi = a if is_hi else hi
    _meet(c, xs[lo], ys[lo], xs[hi], ys[hi])


@njit(cache=True)
def _graph(P, CI, CO, p, dirs, tx, ty):
    """adm[i, l]: shortcut i -> i + l admissible; pen[i, l]: its penalty (squared pixel distances)."""
    n = P.shape[0]
    rows = []
    L = 1
    for i in range(n):
        seen = np.zeros(4, np.int64)
        lin = np.array([tx > 0 or ty > 0, 0, 0, 0, 0, 0], np.int64)  # potrace cone
        dss = np.array([1, 0, 0, 0, 0, 0], np.int64)  # lossless cone
        ok = np.zeros(n, np.bool_)
        ok[1] = True
        last = 1
        for l in range(1, n - 1):
            j, e = (i + l) % n, (i + l - 1) % n
            seen[dirs[e]] = 1
            if seen.sum() == 4:
                break
            _meet(dss, CO[e, 0] - P[i, 0], CO[e, 1] - P[i, 1], CI[e, 0] - P[i, 0], CI[e, 1] - P[i, 1])
            if dss[0] == 0 and lin[0] == 0:
                break
            last = l
            vx, vy = P[j, 0] - P[i, 0], P[j, 1] - P[i, 1]
            if l >= 2 and (vx != 0 or vy != 0):
                ok[l] = (lin[0] == 1 and (lin[5] == 0 or _within(vx, vy, lin[1], lin[2], lin[3], lin[4]))) or (
                    dss[0] == 1 and _cross(dss[1], dss[2], vx, vy) > 0 and _cross(vx, vy, dss[3], dss[4]) > 0)
            if lin[0] == 1 and (abs(vx) > tx or abs(vy) > ty):  # a box around the origin constrains nothing
                _meet_box(lin, vx, vy, tx, ty)
        rows.append(ok[:last + 1].copy())
        L = max(L, last)
    S = np.zeros((2 * n + 1, 5), np.int64)  # prefix sums of x, y, x^2, y^2, xy over two laps
    for k in range(2 * n):
        x, y = p[k % n, 0], p[k % n, 1]
        S[k + 1] = S[k] + np.array([x, y, x * x, y * y, x * y])
    adm = np.zeros((n, L + 1), np.bool_)
    pen = np.zeros((n, L + 1), np.float64)
    for i in range(n):
        xi, yi = p[i, 0], p[i, 1]
        for l in range(1, rows[i].shape[0]):
            adm[i, l] = rows[i][l]
            if adm[i, l] and l >= 2:
                j, a, b, m = (i + l) % n, i + 1, i + l, l - 1
                vx, vy = p[j, 0] - xi, p[j, 1] - yi
                sx, sy = S[b, 0] - S[a, 0], S[b, 1] - S[a, 1]
                sxx = S[b, 2] - S[a, 2] - 2 * xi * sx + m * xi * xi
                syy = S[b, 3] - S[a, 3] - 2 * yi * sy + m * yi * yi
                sxy = S[b, 4] - S[a, 4] - xi * sy - yi * sx + m * xi * yi
                pen[i, l] = (vx * vx * syy - 2 * vx * vy * sxy + vy * vy * sxx) / (vx * vx + vy * vy)
    return adm, pen


@njit(cache=True)
def _best_cycle(adm, pen):
    n, L = adm.shape[0], adm.shape[1] - 1
    reach = np.ones(n, np.int64)
    for i in range(n):
        for l in range(L, 0, -1):
            if adm[i, l]:
                reach[i] = l
                break
    anchor, fewest = 0, n + 1  # the point with the fewest shortcuts passing over it
    for s in range(n):
        c = 0
        for back in range(L + 1):
            c += reach[(s - back) % n] > back
        if c < fewest:
            anchor, fewest = s, c
    INF = 1 << 60
    best_n, best_pen, best = INF, np.inf, np.arange(n)
    cnt, pv, par = np.empty(n + 1, np.int64), np.empty(n + 1), np.empty(n + 1, np.int64)
    for back in range(L + 1):
        a = (anchor - back) % n
        if reach[a] <= back:
            continue
        cnt[:], pv[:] = INF, np.inf
        cnt[0], pv[0] = 0, 0.0
        for q in range(n):
            if cnt[q] < INF:
                s = (a + q) % n
                for l in range(1, min(L, n - q) + 1):
                    if adm[s, l]:
                        r, c2, p2 = q + l, cnt[q] + 1, pv[q] + pen[s, l]
                        if c2 < cnt[r] or (c2 == cnt[r] and p2 < pv[r]):
                            cnt[r], pv[r], par[r] = c2, p2, q
        if cnt[n] < best_n or (cnt[n] == best_n and pv[n] < best_pen):
            best_n, best_pen = cnt[n], pv[n]
            best = np.empty(best_n, np.int64)
            r = n
            for k in range(best_n - 1, -1, -1):
                r = par[r]
                best[k] = (a + r) % n
    return np.sort(best)


def simplify(path, dirs, W: int, H: int, tau: float, G: int = 1000) -> np.ndarray:
    """Sorted indices of the kept corners; tau > 0 admits L-inf tau shortcuts (grid units), tau = 0 lossless only."""
    n = len(path)
    if n <= 4:
        return np.arange(n)
    d = STEPS[dirs]
    right = np.stack([-d[:, 1], d[:, 0]], 1)  # the foreground side of each step
    P = to_grid(path, W, H, G) * np.array([2 * W, 2 * H])
    CI, CO = G * (2 * path + d + right), G * (2 * path + d - right)  # pixel centres inside / outside each step
    return _best_cycle(*_graph(P, CI, CO, path, dirs, round(tau * 2 * W), round(tau * 2 * H)))
