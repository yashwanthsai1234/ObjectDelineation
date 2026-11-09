"""Random and adversarial binary masks for topology tests and benchmarks."""
from __future__ import annotations

import numpy as np
from scipy import ndimage


def smooth_noise(rng: np.random.Generator, h: int, w: int, sigma: float) -> np.ndarray:
    f = ndimage.gaussian_filter(rng.standard_normal((h, w)), sigma, mode="wrap")
    return f / (f.std() or 1.0)


def blob(rng: np.random.Generator, h: int, w: int) -> np.ndarray:
    """One smooth object near the centre, sometimes with holes (polygon-like)."""
    y, x = np.mgrid[:h, :w]
    r = np.hypot((y - h / 2) / (h / 2), (x - w / 2) / (w / 2))
    f = 0.6 - r + 0.25 * smooth_noise(rng, h, w, max(1.0, min(h, w) / 8))
    m = f > 0
    if rng.random() < 0.5:  # punch a few holes
        for _ in range(rng.integers(1, 4)):
            cy, cx, rr = rng.integers(0, h), rng.integers(0, w), rng.uniform(0.5, max(1.0, min(h, w) / 10))
            m &= np.hypot(y - cy, x - cx) > rr
    return m


def noisy(rng: np.random.Generator, h: int, w: int) -> np.ndarray:
    """Thresholded smooth noise: many components, holes and islands-in-holes."""
    return smooth_noise(rng, h, w, rng.uniform(0.7, 4.0)) > rng.normal(0, 0.3)


def bernoulli(rng: np.random.Generator, h: int, w: int) -> np.ndarray:
    """i.i.d. pixels: diagonal contacts and one-pixel holes everywhere."""
    return rng.random((h, w)) < rng.uniform(0.2, 0.8)


def checkerboard(h: int, w: int, phase: int = 0) -> np.ndarray:
    y, x = np.mgrid[:h, :w]
    return (y + x + phase) % 2 == 0


def nested_rings(h: int, w: int, step: int = 2) -> np.ndarray:
    """Concentric square rings: island in hole in island ... (Chebyshev distance bands)."""
    y, x = np.mgrid[:h, :w]
    d = np.minimum(np.minimum(y, x), np.minimum(h - 1 - y, w - 1 - x))
    return (d // step) % 2 == 0


def diamond_rings(h: int, w: int) -> np.ndarray:
    """8-connected (diagonal-only) rings enclosing 4-connected holes; L1 distance bands of width 1."""
    y, x = np.mgrid[:h, :w]
    d = np.abs(y - h // 2) + np.abs(x - w // 2)
    return (d % 3 == 0) & (d > 0)


def adversarial(rng: np.random.Generator, h: int, w: int) -> np.ndarray:
    kind = rng.integers(0, 6)
    if kind == 0:
        return checkerboard(h, w, int(rng.integers(0, 2)))
    if kind == 1:
        return nested_rings(h, w, int(rng.integers(1, 4)))
    if kind == 2:
        return diamond_rings(h, w)
    if kind == 3:  # holes and objects touching the crop edge
        m = noisy(rng, h, w)
        m[0, :] = rng.random(w) < 0.5
        m[:, -1] = rng.random(h) < 0.5
        return m
    if kind == 4:  # sparse dots on a frame: many tiny components inside one hole
        m = rng.random((h, w)) < 0.1
        m[[0, -1], :] = True
        m[:, [0, -1]] = True
        return m
    return ~noisy(rng, h, w)


GENERATORS = {"blob": blob, "noisy": noisy, "bernoulli": bernoulli, "adversarial": adversarial}


def random_mask(rng: np.random.Generator, min_size: int = 8, max_size: int = 256) -> tuple[str, np.ndarray]:
    name = list(GENERATORS)[int(rng.integers(0, len(GENERATORS)))]
    h, w = (int(v) for v in rng.integers(min_size, max_size + 1, size=2))
    return name, GENERATORS[name](rng, h, w)
