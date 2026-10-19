import sys, time, numpy as np, numba as nb, shapely
sys.path.insert(0, ".")
from poly_fillcheck import brute, scan, gen, libs

@nb.njit(cache=True)
def onedge(V, S):
    n = V.shape[0]; B = np.zeros((S, S), np.bool_)
    for k in range(n):
        ax, ay = V[k]; bx, by = V[(k + 1) % n]
        i0 = max(0, (min(ax, bx) // 1000 - 1) // 2); i1 = min(S - 1, (max(ax, bx) // 1000) // 2 + 1)
        j0 = max(0, (min(ay, by) // 1000 - 1) // 2); j1 = min(S - 1, (max(ay, by) // 1000) // 2 + 1)
        for j in range(j0, j1 + 1):
            Y = 1000 * (2 * j + 1)
            if Y < min(ay, by) or Y > max(ay, by): continue
            for i in range(i0, i1 + 1):
                X = 1000 * (2 * i + 1)
                if X < min(ax, bx) or X > max(ax, bx): continue
                if (bx - ax) * (Y - ay) - (by - ay) * (X - ax) == 0: B[j, i] = True
    return B

S = int(sys.argv[1]); N = int(sys.argv[2]); rng = np.random.default_rng(1)
st = {}; bm = 0; bpx = 0
for t in range(N):
    kind = ("star", "random", "rect")[t % 3]; P = gen(kind, rng)
    if len(np.unique(P, axis=0)) < 3: continue
    V = P * 2 * S; W = scan(V, S); B = onedge(V, S); bpx += int(B.sum())
    if t < 200 and S <= 256:
        Wb, Bb = brute(V, S); bm += int((Wb != W).sum()) + int((Bb != B).sum())
    eo = (W & 1).astype(bool); nz = W != 0; selfx = bool((eo != nz).any())
    for name, m in libs(P, S).items():
        for rule, ref in (("eo", eo), ("nz", nz)):
            d = m != ref; off = int((d & ~B).sum()); on = int((d & B).sum())
            s = st.setdefault((name, rule, "selfoverlap" if selfx else "no-overlap"), [0, 0, 0, 0])
            s[0] += 1; s[1] += off > 0; s[2] += off; s[3] += on
print(f"S={S} N={N}: brute-vs-scan mismatches={bm}; total on-edge pixel centres={bpx}")
print(f"{'lib':26s} rule class        polys  polys_offB_diff offB_px onB_px")
for k in sorted(st):
    s = st[k]; print(f"{k[0]:26s} {k[1]:4s} {k[2]:12s} {s[0]:5d} {s[1]:6d} {s[2]:9d} {s[3]:7d}")
