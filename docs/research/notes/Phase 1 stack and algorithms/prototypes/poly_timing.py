import sys, time, numpy as np, shapely, cv2, skimage.draw
sys.path.insert(0, "."); sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[5] / "src"))
from poly_fillcheck import scan, gen, libs
from PIL import Image, ImageDraw
from matplotlib.path import Path
from rasterio.features import rasterize
from rasterio.transform import Affine
from objdel.checker.polygon_io import ring_problems
import numba as nb

@nb.njit(cache=True)
def seg_quadratic(P):  # exact int64 O(n^2) all-pairs incl. adjacent collinear overlap; returns code 0 ok 1 touch 2 cross
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
            if adj:  # adjacent edges share a vertex: only collinear back-tracking (spike) is a problem
                if d3 == 0 and d4 == 0 and ((bx-ax)*(dx-cx)+(by-ay)*(dy-cy)) < 0: worst = 1
                continue
            for d, px, py, qx, qy, rx, ry in ((d1, cx, cy, dx, dy, ax, ay), (d2, cx, cy, dx, dy, bx, by), (d3, ax, ay, bx, by, cx, cy), (d4, ax, ay, bx, by, dx, dy)):
                if d == 0 and min(px, qx) <= rx <= max(px, qx) and min(py, qy) <= ry <= max(py, qy): worst = 1
    return worst

def t(f, reps=50):
    f(); ts = []
    for _ in range(reps):
        a = time.perf_counter(); f(); ts.append(time.perf_counter() - a)
    return np.median(ts) * 1e6

rng = np.random.default_rng(3); S = 512
for nv in (16, 64, 200):
    c = np.array([500, 500]); tt = np.sort(rng.uniform(0, 2*np.pi, nv)); r = 300 * rng.uniform(0.6, 1, nv)
    P = np.rint(np.c_[c[0] + r*np.cos(tt), c[1] + r*np.sin(tt)]).astype(np.int64)
    x = P[:, 0]*S/1000.; y = P[:, 1]*S/1000.; V = P*2*S
    jj, ii = np.mgrid[0:S, 0:S]; cx = ii.ravel()+.5; cy = jj.ravel()+.5; C = np.c_[cx, cy]
    g = shapely.Polygon(np.c_[x, y]); gi = shapely.Polygon(P)
    def pil():
        im = Image.new("1", (S, S), 0); ImageDraw.Draw(im).polygon(list(zip(x, y)), fill=1); return np.asarray(im)
    def sk():
        m = np.zeros((S, S), bool); rr, cc2 = skimage.draw.polygon(y-.5, x-.5, (S, S)); m[rr, cc2] = 1
    def cvf():
        m = np.zeros((S, S), np.uint8); cv2.fillPoly(m, [np.rint((np.c_[x, y]-.5)*256).astype(np.int32)], 1, shift=8)
    res = {
        "exact scan (numba, eo+nz from one winding pass)": t(lambda: scan(V, S)),
        "PIL ImageDraw.polygon": t(pil),
        "cv2.fillPoly": t(cvf),
        "skimage.draw.polygon": t(sk, 5),
        "matplotlib contains_points": t(lambda: Path(np.c_[x, y]).contains_points(C), 5),
        "shapely.contains_xy (S^2 pts)": t(lambda: shapely.contains_xy(g, cx, cy), 5),
        "rasterio.rasterize": t(lambda: rasterize([(g, 1)], out_shape=(S, S), transform=Affine.identity())),
        "shapely.is_valid": t(lambda: shapely.is_valid(gi)),
        "shapely.is_valid_reason": t(lambda: shapely.is_valid_reason(gi)),
        "repo ring_problems (pure-Python O(n^2))": t(lambda: ring_problems(P), 3),
        "exact O(n^2) seg test (numba int64)": t(lambda: seg_quadratic(P)),
    }
    m = scan(V, S)
    print(f"n={nv} S={S} area_px={(m!=0).sum()}")
    for k, v in res.items(): print(f"   {k:50s} {v:10.1f} us")
# over/under fill split for PIL and cv2
rng = np.random.default_rng(5); over = {"PIL": [0, 0], "cv2": [0, 0]}
for k in range(200):
    P = gen("star", rng); S2 = 250; W = scan(P*2*S2, S2) != 0; L = libs(P, S2)
    for name, key in (("PIL", "PIL fill off=0.0"), ("cv2", "cv2.fillPoly shift=8")):
        over[name][0] += int((L[key] & ~W).sum()); over[name][1] += int((~L[key] & W).sum())
print("star polys S=250: [lib-only px, ref-only px]", over)
