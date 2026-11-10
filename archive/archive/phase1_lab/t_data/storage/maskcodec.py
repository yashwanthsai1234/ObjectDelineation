"""Bit-exact mask codecs for the Parquet store.

poly: exact rectilinear rings from objdel.checker.polygonize (pixel-corner
      vertices, outer rings and holes). Decoded by an even-odd scanline over
      the vertical edges, so decode(encode(m)) == m bit for bit (8-connected
      foreground, as the checker defines it).
      Stored as {"h","w","rings": [[x0,y0,x1,y1,...], ...]} with int32 vertices.
rle:  row-major run lengths starting with a background run (COCO-like, but
      row-major), as a list of uint32.
"""
import hashlib
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(REPO / "src")]
from objdel.checker import polygonize  # noqa: E402
from objdel.checker.certify import crop  # noqa: E402


def poly_encode(mask):
    m = np.asarray(mask, bool)
    c, off = crop(m)
    rings = polygonize(c, off) if c.size else []
    return {"h": int(m.shape[0]), "w": int(m.shape[1]),
            "rings": [r.pts.astype(np.int32).ravel().tolist() for r in rings]}


def poly_decode(enc):
    """Even-odd fill from the vertical edges, vectorised, inside the rings' bbox."""
    h, w = enc["h"], enc["w"]
    out = np.zeros((h, w), bool)
    if not enc["rings"]:
        return out
    P = [np.asarray(f, np.int64).reshape(-1, 2) for f in enc["rings"]]
    allp = np.concatenate(P)
    x0, y0 = allp.min(0)
    x1, y1 = allp.max(0)
    t = np.zeros((y1 - y0 + 1, x1 - x0 + 1), np.int32)
    for p in P:
        q = np.roll(p, -1, axis=0)
        v = p[:, 0] == q[:, 0]
        x = p[v, 0] - x0
        a, b = np.minimum(p[v, 1], q[v, 1]) - y0, np.maximum(p[v, 1], q[v, 1]) - y0
        np.add.at(t, (a, x), 1)
        np.add.at(t, (b, x), -1)
    cover = np.cumsum(t, axis=0)[:-1, :-1]  # vertical edges crossing each pixel row, per column x
    out[y0:y1, x0:x1] = (np.cumsum(cover, axis=1) & 1).astype(bool)
    return out


def rle_encode(mask):
    m = np.asarray(mask, bool)
    flat = m.ravel()
    change = np.flatnonzero(flat[1:] != flat[:-1]) + 1
    bounds = np.concatenate([[0], change, [flat.size]])
    runs = np.diff(bounds)
    if flat.size and flat[0]:
        runs = np.concatenate([[0], runs])
    return {"h": int(m.shape[0]), "w": int(m.shape[1]), "counts": runs.astype(np.uint32).tolist()}


def rle_decode(enc):
    counts = np.asarray(enc["counts"], np.int64)
    vals = np.arange(counts.size) % 2 == 1
    return np.repeat(vals, counts).reshape(enc["h"], enc["w"])


def mask_sha(mask):
    m = np.asarray(mask, bool)
    return hashlib.sha256(np.packbits(m).tobytes() + f"{m.shape}".encode()).hexdigest()
