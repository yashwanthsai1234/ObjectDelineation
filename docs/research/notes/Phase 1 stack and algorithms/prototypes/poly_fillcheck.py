"""Compare polygon rasterizers against an exact integer pixel-centre reference.
Vertices: ints on 0..1000 grid; pixel coord x = v*S/1000; pixel (i,j) centre (i+.5, j+.5).
Exact frame (x2000): V = 2*S*v, C = 1000*(2i+1)."""
import sys, time, numpy as np, numba as nb
from PIL import Image, ImageDraw
import skimage.draw, cv2, shapely
from matplotlib.path import Path
from rasterio.features import rasterize
from rasterio.transform import Affine

@nb.njit(cache=True)
def brute(V, S):  # V: (n,2) int64 in x2000 frame; returns winding per pixel
    n = V.shape[0]; W = np.zeros((S, S), np.int32); B = np.zeros((S, S), np.bool_)
    for j in range(S):
        Y = 1000 * (2 * j + 1)
        for i in range(S):
            X = 1000 * (2 * i + 1); w = 0
            for k in range(n):
                ax, ay = V[k]; bx, by = V[(k + 1) % n]
                cr = (bx - ax) * (Y - ay) - (by - ay) * (X - ax)
                if cr == 0 and min(ax, bx) <= X <= max(ax, bx) and min(ay, by) <= Y <= max(ay, by):
                    B[j, i] = True
                if ay <= Y < by and cr > 0: w += 1
                elif by <= Y < ay and cr < 0: w -= 1
            W[j, i] = w
    return W, B

@nb.njit(cache=True)
def scan(V, S):  # exact scanline + difference array (signed accumulation), O(rows*edges + pixels)
    n = V.shape[0]; W = np.zeros((S, S + 1), np.int32)
    for k in range(n):
        ax, ay = V[k]; bx, by = V[(k + 1) % n]
        if ay == by: continue
        s = 1 if by > ay else -1
        y0, y1 = min(ay, by), max(ay, by)
        # rows j with y0 <= 1000(2j+1) < y1
        jlo = max(0, -((-(y0 - 1000)) // 2000)); jhi = min(S - 1, (y1 - 1 - 1000) // 2000)
        if ay < by: x0, yy0, dx, dy = ax, ay, bx - ax, by - ay
        else: x0, yy0, dx, dy = bx, by, ax - bx, ay - by
        for j in range(jlo, jhi + 1):
            Y = 1000 * (2 * j + 1)
            num = (Y - yy0) * dx + x0 * dy          # x_cross*dy ; pixel i counts if 1000(2i+1)*dy < num
            kmax = (num - 1) // (1000 * dy)         # largest odd-or-even k with k*1000*dy < num
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

def gen(kind, rng):
    n = int(rng.integers(3, 120))
    if kind == "star":
        c = rng.integers(200, 800, 2); r = rng.integers(20, 200)
        t = np.sort(rng.uniform(0, 2 * np.pi, n)); rad = r * rng.uniform(0.3, 1, n)
        P = np.c_[c[0] + rad * np.cos(t), c[1] + rad * np.sin(t)]
    elif kind == "random":
        c = rng.integers(200, 800, 2); P = c + rng.integers(-200, 200, (n, 2))
    else:  # rectilinear staircase-ish with long axis-aligned runs
        c = rng.integers(200, 800, 2); m = max(2, n // 2)
        t = np.sort(rng.uniform(0, 2 * np.pi, m)); rad = rng.integers(20, 200) * rng.uniform(0.5, 1, m)
        Q = np.c_[c[0] + rad * np.cos(t), c[1] + rad * np.sin(t)]
        P = np.empty((2 * m, 2)); P[0::2] = Q; P[1::2, 0] = np.roll(Q[:, 0], -1); P[1::2, 1] = Q[:, 1]
    return np.clip(np.rint(P), 0, 1000).astype(np.int64)

def libs(P, S):
    x = P[:, 0] * S / 1000.0; y = P[:, 1] * S / 1000.0
    out = {}
    for off in (0.0, -0.5):
        im = Image.new("1", (S, S), 0); ImageDraw.Draw(im).polygon(list(zip(x + off, y + off)), fill=1)
        out[f"PIL fill off={off}"] = np.asarray(im, bool)
        im = Image.new("1", (S, S), 0); ImageDraw.Draw(im).polygon(list(zip(x + off, y + off)), fill=1, outline=1)
        out[f"PIL fill+outline off={off}"] = np.asarray(im, bool)
    m = np.zeros((S, S), bool); rr, cc = skimage.draw.polygon(y - 0.5, x - 0.5, (S, S)); m[rr, cc] = True
    out["skimage centre@int"] = m
    m = np.zeros((S, S), np.uint8)
    cv2.fillPoly(m, [np.rint((np.c_[x, y] - 0.5) * 256).astype(np.int32)], 1, lineType=cv2.LINE_8, shift=8)
    out["cv2.fillPoly shift=8"] = m.astype(bool)
    jj, ii = np.mgrid[0:S, 0:S]; C = np.c_[ii.ravel() + 0.5, jj.ravel() + 0.5]
    out["mpl contains_points"] = Path(np.c_[x, y]).contains_points(C).reshape(S, S)
    g = shapely.Polygon(np.c_[x, y])
    out["shapely contains_xy"] = shapely.contains_xy(g, C[:, 0], C[:, 1]).reshape(S, S)
    out["rasterio rasterize"] = rasterize([(g, 1)], out_shape=(S, S), transform=Affine.identity(), fill=0, dtype="uint8").astype(bool)
    return out

if __name__ == "__main__":
    S = int(sys.argv[1]); N = int(sys.argv[2]); rng = np.random.default_rng(0)
    stats = {}; nbrute_mismatch = 0; nb_bound = 0
    for t in range(N):
        kind = ("star", "random", "rect")[t % 3]; P = gen(kind, rng)
        if len(np.unique(P, axis=0)) < 3: continue
        V = P * 2 * S
        Wr = scan(V, S)
        if t < 300 and S <= 256:
            Wb, B = brute(V, S); nbrute_mismatch += int((Wb != Wr).sum()); nb_bound += int(B.sum())
        eo = (Wr & 1).astype(bool); nz = Wr != 0
        valid = shapely.Polygon(P).is_valid
        for name, m in libs(P, S).items():
            for rule, ref in (("eo", eo), ("nz", nz)):
                k = (name, rule, kind if kind != "random" else ("random-valid" if valid else "random-invalid"))
                d = int((m != ref).sum()); s = stats.setdefault(k, [0, 0, 0]); s[0] += 1; s[1] += d > 0; s[2] += d
    print(f"S={S} N={N} brute-vs-scan mismatched pixels={nbrute_mismatch} (boundary centres seen={nb_bound})")
    for k in sorted(stats):
        s = stats[k]; print(f"{k[0]:28s} {k[1]} {k[2]:15s} polys={s[0]:5d} polys_with_diff={s[1]:5d} diff_px={s[2]}")
