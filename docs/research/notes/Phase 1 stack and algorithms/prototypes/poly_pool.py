import time, numpy as np, numba as nb
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from poly_fillcheck import gen

@nb.njit(nogil=True, cache=True)
def scan_bbox(V, S, W):  # accumulate signed winding difference into W (S, S+1)
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

@nb.njit(nogil=True, cache=True)
def reward(rings, offs, S, gt):  # rings: concatenated int64 (m,2); offs ring offsets
    W = np.zeros((S, S + 1), np.int32)
    for r in range(offs.shape[0] - 1):
        scan_bbox(rings[offs[r]:offs[r + 1]] * 2 * S, S, W)
    inter_eo = 0; uni_eo = 0; overlap = 0
    for j in range(S):
        acc = 0
        for i in range(S):
            acc += W[j, i]; eo = (acc & 1) == 1; nz = acc != 0; g = gt[j, i]
            inter_eo += eo and g; uni_eo += eo or g; overlap += eo != nz
    return inter_eo / max(uni_eo, 1), overlap

rng = np.random.default_rng(0); S = 512
jobs = []
for _ in range(256):
    rs = [gen("star", rng) for _ in range(3)]
    offs = np.cumsum([0] + [len(r) for r in rs]).astype(np.int64)
    jobs.append((np.concatenate(rs), offs))
gt = np.zeros((S, S), np.bool_); gt[100:400, 100:400] = True
reward(jobs[0][0], jobs[0][1], S, gt)
for nt in (1, 2, 4, 8):
    t = time.perf_counter()
    with ThreadPoolExecutor(nt) as ex: list(ex.map(lambda j: reward(j[0], j[1], S, gt), jobs))
    print(f"threads={nt}: 256 rollouts (3 rings each, S={S}) in {1e3*(time.perf_counter()-t):.1f} ms")
