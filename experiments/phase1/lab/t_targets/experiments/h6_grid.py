"""H6: does the 0..1000 grid mapping lose information, and do crop sizes dividing 1000 avoid it?

(a) rounding-only effect: simplify in the *pixel* frame (vertices on pixel corners, as a pixel-space
    simplifier would), rasterize exactly before rounding (G = S) and after rounding to the 0..1000 grid;
    count flipped pixels, topology changes and rings that become illegal, per crop side S.
(b) identity round trip: exact rings (no simplification) -> grid -> re-raster == mask ?
Masks: fixture masks padded to squares of side S in a sweep (native side, divisors of 1000, and
non-divisors), plus synthetic smooth blobs with holes and parts.
    python experiments/h6_grid.py results/h6.json
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from objdel.checker.polygonize import polygonize  # noqa: E402
from targets import fastgeom as F  # noqa: E402
from targets.builder import betti  # noqa: E402
from targets.simplify import douglas_peucker, lattice_path, minseg, to_grid  # noqa: E402

FIX = HERE.parent.parent / "fixtures"


def blobs(rng, S, n):
    out = []
    while len(out) < n:
        sig = rng.uniform(0.04, 0.12) * S
        f = ndimage.gaussian_filter(rng.standard_normal((S, S)), sig)
        m = f > np.quantile(f, rng.uniform(0.6, 0.9))
        m[:2] = m[-2:] = False
        m[:, :2] = m[:, -2:] = False
        if m.any():
            out.append(m)
    return out


def pad_to(m, S):
    H, W = m.shape
    if max(H, W) > S:
        return None
    o = np.zeros((S, S), bool)
    o[(S - H) // 2:(S - H) // 2 + H, (S - W) // 2:(S - W) // 2 + W] = m
    return o


def evaluate(m, S, simp):
    rings_px = []
    for r in polygonize(m):
        p, d = lattice_path(r.pts)
        if simp == "identity":
            K = np.nonzero(d != np.roll(d, 1))[0]
        elif simp == "dp1px":
            K = douglas_peucker(p, d, S, S, G=S, eps=1.0)  # G=S: grid == pixel lattice
        else:  # potrace-like, 1 px, pixel frame
            K = minseg(p, d, S, S, G=S, mode="linf", t_px=1.0)
        rings_px.append(p[K])
    pre = F.raster(rings_px, S, S, G=S)
    rings_g = [to_grid(r, S, S, 1000) for r in rings_px]
    post = F.raster(rings_g, S, S, G=1000)
    legal_pre = F.legality(rings_px)["ok"]
    legal_post = F.legality(rings_g)["ok"]
    flips = int((pre != post).sum())
    u = np.logical_or(pre, post).sum()
    return {"flips": flips, "iou_round": float(np.logical_and(pre, post).sum() / u) if u else 1.0,
            "topo_change": betti(pre) != betti(post), "illegal_new": bool(legal_pre and not legal_post),
            "src_iou_post": float(np.logical_and(post, m).sum() / np.logical_or(post, m).sum()),
            "identity_exact": bool(np.array_equal(post, m)) if simp == "identity" else None}


def main(out):
    rng = np.random.default_rng(0)
    res = []
    Ss = [100, 125, 128, 200, 250, 256, 333, 384, 400, 500, 512, 640, 1000, 1024, 1200]
    for S in Ss:
        for i, m in enumerate(blobs(rng, S, 12)):
            for simp in ("identity", "dp1px", "minseg1px"):
                r = evaluate(m, S, simp)
                res.append({"src": "blob", "S": S, "div": 1000 % S == 0, "simp": simp, **r})
    z = np.load(FIX / "riga_experts.npz")
    keys = sorted(z.files)[:40]
    for k in keys:
        m0 = z[k][0].astype(bool)
        n = max(m0.shape)
        for S in sorted({n, 500 if n <= 500 else 1000, n + 7, 1000, 1250}):
            m = pad_to(m0, S)
            if m is None:
                continue
            for simp in ("identity", "dp1px", "minseg1px"):
                r = evaluate(m, S, simp)
                res.append({"src": "riga", "S": S, "div": 1000 % S == 0, "simp": simp, **r})
    z = np.load(FIX / "lidc_experts.npz")
    for k in sorted(z.files)[::6]:
        m0 = z[k][0].astype(bool)
        n = max(m0.shape)
        for S in sorted({n, 50 if n <= 50 else 100, n + 3}):
            m = pad_to(m0, S)
            for simp in ("identity", "minseg1px"):
                r = evaluate(m, S, simp)
                res.append({"src": "lidc", "S": S, "div": 1000 % S == 0, "simp": simp, **r})
    Path(out).write_text(json.dumps(res))


if __name__ == "__main__":
    main(sys.argv[1])
