"""Vendored copy of the numba prototypes from the repo notes
(docs/research/notes/Phase 1 stack and algorithms/prototypes/poly_fillcheck.py `scan`,
poly_pool.py `scan_bbox`/`reward`, poly_timing.py `seg_quadratic`), unchanged except:
cache disabled (so nothing is written next to the repo files) and nogil on the kernels.
Square crops only (S x S, G = 1000), exactly as in the prototypes."""
import numba as nb
import numpy as np


@nb.njit(nogil=True)
def scan(V, S):  # V: (n,2) int64 in x2000 frame; returns winding per pixel
    n = V.shape[0]; W = np.zeros((S, S + 1), np.int32)
    for k in range(n):
        ax, ay = V[k]; bx, by = V[(k + 1) % n]
        if ay == by: continue
        s = 1 if by > ay else -1
        y0, y1 = min(ay, by), max(ay, by)
        jlo = max(0, -((-(y0 - 1000)) // 2000)); jhi = min(S - 1, (y1 - 1 - 1000) // 2000)
        if ay < by: x0, yy0, dx, dy = ax, ay, bx - ax, by - ay
        else: x0, yy0, dx, dy = bx, by, ax - bx, ay - by
        for j in range(jlo, jhi + 1):
            Y = 1000 * (2 * j + 1)
            num = (Y - yy0) * dx + x0 * dy
            kmax = (num - 1) // (1000 * dy)
            imax = (kmax - 1) // 2
            if imax < 0: continue
            if imax > S - 1: imax = S - 1
            W[j, 0] += s; W[j, imax + 1] -= s
    out = np.empty((S, S), np.int32)
    for j in range(S):
        acc = 0
        for i in range(S):
            acc += W[j, i]; out[j, i] = acc
    return out


@nb.njit(nogil=True)
def scan_bbox(V, S, W):
    n = V.shape[0]
    for k in range(n):
        ax, ay = V[k]; bx, by = V[(k + 1) % n]
        if ay == by: continue
        s = 1 if by > ay else -1
        y0, y1 = min(ay, by), max(ay, by)
        jlo = max(0, -((-(y0 - 1000)) // 2000)); jhi = min(S - 1, (y1 - 1 - 1000) // 2000)
        if ay < by: x0, yy0, dx, dy = ax, ay, bx - ax, by - ay
        else: x0, yy0, dx, dy = bx, by, ax - bx, ay - by
        for j in range(jlo, jhi + 1):
            num = (1000 * (2 * j + 1) - yy0) * dx + x0 * dy
            imax = ((num - 1) // (1000 * dy) - 1) // 2
            if imax < 0: continue
            if imax > S - 1: imax = S - 1
            W[j, 0] += s; W[j, imax + 1] -= s


@nb.njit(nogil=True)
def seg_quadratic(P):  # 0 ok, 1 touch, 2 cross
    n = P.shape[0]; worst = 0
    for i in range(n):
        ax, ay = P[i]; bx, by = P[(i + 1) % n]
        for j in range(i + 1, n):
            cx, cy = P[j]; dx, dy = P[(j + 1) % n]
            d1 = (dx-cx)*(ay-cy)-(dy-cy)*(ax-cx); d2 = (dx-cx)*(by-cy)-(dy-cy)*(bx-cx)
            d3 = (bx-ax)*(cy-ay)-(by-ay)*(cx-ax); d4 = (bx-ax)*(dy-ay)-(by-ay)*(dx-ax)
            adj = (j == i + 1) or (i == 0 and j == n - 1)
            if ((d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0)) and ((d3 > 0 and d4 < 0) or (d3 < 0 and d4 > 0)):
                return 2
            if adj:
                if d3 == 0 and d4 == 0 and ((bx-ax)*(dx-cx)+(by-ay)*(dy-cy)) < 0: worst = 1
                continue
            for d, px, py, qx, qy, rx, ry in ((d1, cx, cy, dx, dy, ax, ay), (d2, cx, cy, dx, dy, bx, by), (d3, ax, ay, bx, by, cx, cy), (d4, ax, ay, bx, by, dx, dy)):
                if d == 0 and min(px, qx) <= rx <= max(px, qx) and min(py, qy) <= ry <= max(py, qy): worst = 1
    return worst


@nb.njit(nogil=True)
def reward(rings, offs, S, gt):
    """poly_pool.reward extended with the O(n^2) validity of every ring (realistic full reward)."""
    W = np.zeros((S, S + 1), np.int32)
    worst = 0
    for r in range(offs.shape[0] - 1):
        P = rings[offs[r]:offs[r + 1]]
        scan_bbox(P * 2 * S, S, W)
        worst = max(worst, seg_quadratic(P))
    inter_eo = 0; uni_eo = 0; overlap = 0
    for j in range(S):
        acc = 0
        for i in range(S):
            acc += W[j, i]; eo = (acc & 1) == 1; nz = acc != 0; g = gt[j, i]
            inter_eo += eo and g; uni_eo += eo or g; overlap += eo != nz
    return inter_eo / max(uni_eo, 1), overlap, worst
