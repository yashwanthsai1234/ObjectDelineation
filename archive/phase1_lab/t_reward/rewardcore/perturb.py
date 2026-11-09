"""Controlled mask perturbations: boundary moves, shifts, specks, pinholes, bridges, jitter."""
from __future__ import annotations

import numpy as np

from .distances import edt
from .synth import smooth_noise


def signed_distance(m: np.ndarray) -> np.ndarray:
    """> 0 inside (distance to background), <= 0 outside (minus distance to foreground)."""
    m = np.asarray(m, bool)
    return np.where(m, edt(m), -edt(~m)).astype(np.float32)


def dilate(m: np.ndarray, d: float) -> np.ndarray:
    """Euclidean dilation by d px (pixels within distance d of the mask)."""
    return np.asarray(m, bool) | (edt(~np.asarray(m, bool)) <= d) if d > 0 else np.asarray(m, bool).copy()


def erode(m: np.ndarray, d: float) -> np.ndarray:
    """Euclidean erosion by d px (pixels farther than d from the background)."""
    m = np.asarray(m, bool)
    return edt(m) > d if d > 0 else m.copy()


def grow(m: np.ndarray, d: float) -> np.ndarray:
    """Signed boundary move: d > 0 dilates, d < 0 erodes."""
    return dilate(m, d) if d >= 0 else erode(m, -d)


def shift(m: np.ndarray, dy: int, dx: int) -> np.ndarray:
    m = np.asarray(m, bool)
    out = np.zeros_like(m)
    h, w = m.shape
    ys, yd = (slice(0, h - dy), slice(dy, h)) if dy >= 0 else (slice(-dy, h), slice(0, h + dy))
    xs, xd = (slice(0, w - dx), slice(dx, w)) if dx >= 0 else (slice(-dx, w), slice(0, w + dx))
    out[yd, xd] = m[ys, xs]
    return out


def _disk_offsets(area: int) -> np.ndarray:
    """Pixel offsets of the most compact digital disk with exactly `area` pixels (4-connected)."""
    r = int(np.ceil(np.sqrt(area / np.pi))) + 2
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    order = np.lexsort((xx.ravel(), yy.ravel(), np.hypot(yy, xx).ravel()))
    return np.stack([yy.ravel()[order][:area], xx.ravel()[order][:area]], 1)


def add_speck(m: np.ndarray, area: int, gap: float, rng: np.random.Generator,
              ref: np.ndarray | None = None) -> np.ndarray | None:
    """Add a detached compact speck of `area` px whose centre is about `gap` px outside `ref`
    (default: the mask itself; pass the rater union to place it beyond the band).

    The speck never touches (8-adjacency) the mask or ref. Returns None if no room in the crop.
    """
    m = np.asarray(m, bool)
    off = _disk_offsets(area)
    reach = np.abs(off).max()
    dist = edt(~(m | (np.asarray(ref, bool) if ref is not None else m)))
    h, w = m.shape
    cand = np.argwhere((np.abs(dist - gap) <= 0.5) & (dist > reach + 1.5))
    cand = cand[(cand[:, 0] >= reach) & (cand[:, 0] < h - reach) & (cand[:, 1] >= reach) & (cand[:, 1] < w - reach)]
    if not len(cand):
        return None
    c = cand[rng.integers(len(cand))]
    out = m.copy()
    out[c[0] + off[:, 0], c[1] + off[:, 1]] = True
    return out


def punch_hole(m: np.ndarray, area: int, rng: np.random.Generator, centre: bool = False) -> np.ndarray | None:
    """Remove a compact `area`-px region strictly inside the mask (a new hole). None if no room."""
    m = np.asarray(m, bool)
    off = _disk_offsets(area)
    reach = np.abs(off).max()
    depth = edt(m)
    if centre:
        cand = np.argwhere(depth == depth.max())
    else:
        cand = np.argwhere(depth > reach + 1.5)
    if not len(cand) or depth.max() <= reach + 1.5:
        return None
    c = cand[rng.integers(len(cand))]
    out = m.copy()
    out[c[0] + off[:, 0], c[1] + off[:, 1]] = False
    return out


def bridge(m: np.ndarray, a: np.ndarray, b: np.ndarray, width: int = 1) -> np.ndarray:
    """Join regions a and b of mask m with a straight strip between their closest pixels."""
    pa, pb = np.argwhere(a), np.argwhere(b)
    # closest pair (brute force over contour-sized sets is fine for fixtures)
    ia = pa[np.linspace(0, len(pa) - 1, min(len(pa), 400)).astype(int)]
    ib = pb[np.linspace(0, len(pb) - 1, min(len(pb), 400)).astype(int)]
    d = ((ia[:, None, :] - ib[None, :, :]) ** 2).sum(-1)
    i, j = np.unravel_index(np.argmin(d), d.shape)
    p, q = ia[i].astype(float), ib[j].astype(float)
    n = int(np.ceil(np.abs(q - p).max())) * 2 + 1
    out = np.asarray(m, bool).copy()
    for t in np.linspace(0, 1, n):
        y, x = np.round(p + t * (q - p)).astype(int)
        out[max(0, y - width // 2): y + (width + 1) // 2, max(0, x - width // 2): x + (width + 1) // 2] = True
    return out


def jitter(m: np.ndarray, amp: float, rng: np.random.Generator, corr: float = 6.0) -> np.ndarray:
    """Smooth random normal displacement of the boundary by up to ~amp px (may change topology; check)."""
    m = np.asarray(m, bool)
    noise = smooth_noise(rng, *m.shape, sigma=corr)
    return signed_distance(m) + amp * noise > 0.5
