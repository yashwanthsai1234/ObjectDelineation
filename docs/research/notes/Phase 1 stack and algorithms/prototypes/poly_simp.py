import sys, numpy as np, shapely, potrace, time
from scipy import ndimage as ndi
sys.path.insert(0, "."); sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[5] / "src"))
from poly_fillcheck import scan
from objdel.checker.polygonize import polygonize
from objdel.checker.topology import betti

def raster_rings(rings, S):  # even-odd over all rings (pixel coords, float -> exact via x2000 frame on 0..1000 grid)
    W = np.zeros((S, S), np.int32)
    for r in rings:
        if len(r) >= 3: W += scan(np.asarray(r, np.int64) * 2 * S, S)
    return (W & 1).astype(bool)

def to_grid(pts, S): return np.rint(np.asarray(pts, float) * 1000.0 / S).astype(np.int64)

def blob(S, rng, sigma):
    f = ndi.gaussian_filter(rng.standard_normal((S, S)), sigma)
    m = f > np.quantile(f, rng.uniform(0.55, 0.8))
    lab, n = ndi.label(m, structure=np.ones((3, 3)))
    if n == 0: return m
    sizes = ndi.sum(m, lab, range(1, n + 1)); return lab == (1 + int(np.argmax(sizes)))

def stats(name, rings_px, mask, S, acc, t):
    g = [to_grid(r, S) for r in rings_px]
    # dedupe consecutive duplicate vertices created by quantisation
    g2 = []
    for r in g:
        keep = np.r_[True, np.any(np.diff(r, axis=0) != 0, axis=1)]; r = r[keep]
        if len(r) > 1 and (r[0] == r[-1]).all(): r = r[:-1]
        g2.append(r)
    m2 = raster_rings(g2, S)
    iou = (m2 & mask).sum() / max(1, (m2 | mask).sum())
    topo = betti(m2) == betti(mask)
    lr = [shapely.LinearRing(r) for r in g2 if len(r) >= 3]
    simple = all(shapely.is_simple(l) for l in lr)
    noncross = True
    for i in range(len(lr)):
        for j in range(i + 1, len(lr)):
            if shapely.intersects(lr[i], lr[j]): noncross = False
    a = acc.setdefault(name, dict(v=[], iou=[], topo=[], ok=[], t=[]))
    a["v"].append(sum(len(r) for r in g2)); a["iou"].append(iou); a["topo"].append(topo)
    a["ok"].append(simple and noncross); a["t"].append(t)

rng = np.random.default_rng(7); acc = {}
for k in range(int(sys.argv[2])):
    S = int(sys.argv[1]); mask = blob(S, rng, rng.uniform(4, 12))
    t0 = time.perf_counter(); rings = polygonize(mask); t1 = time.perf_counter()
    stats("exact pixel-edge rings", [r.pts for r in rings], mask, S, acc, t1 - t0)
    for tol in (0.5, 1.0):
        for pt in (True, False):
            t0 = time.perf_counter()
            out = [np.asarray(shapely.simplify(shapely.LinearRing(r.pts), tol, preserve_topology=pt).coords)[:-1] for r in rings]
            stats(f"DP per-ring tol={tol} preserve_topology={pt}", out, mask, S, acc, time.perf_counter() - t0)
    t0 = time.perf_counter()
    curves = potrace.Bitmap(~mask).trace(turdsize=0, turnpolicy=potrace.POTRACE_TURNPOLICY_BLACK, alphamax=0.0, opticurve=False)
    out = [np.array([(s.c.x, s.c.y) for s in c.segments]) for c in curves]
    stats("potrace polygon stage (alphamax=0)", out, mask, S, acc, time.perf_counter() - t0)
print(f"S={sys.argv[1]} masks={sys.argv[2]}  (vertices after rounding to 0..1000 grid; IoU/topology from exact even-odd re-raster)")
for n, a in acc.items():
    print(f"{n:45s} verts_mean={np.mean(a['v']):7.1f} IoU_min={np.min(a['iou']):.4f} IoU_mean={np.mean(a['iou']):.4f} topo_kept={np.mean(a['topo']):.2f} simple&noncross={np.mean(a['ok']):.2f} t_ms={1e3*np.median(a['t']):.2f}")
