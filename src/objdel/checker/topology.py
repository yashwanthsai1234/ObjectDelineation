"""Betti numbers by connected-component labelling, and the size filter.

Convention: 8-connected foreground, 4-connected background. A hole is a
background component that does not touch the image border.
"""
import numpy as np
from scipy import ndimage

FG8 = np.ones((3, 3), dtype=bool)
BG4 = ndimage.generate_binary_structure(2, 1)


def betti(mask):
    """(beta0, beta1) of a 2D binary mask."""
    m = np.asarray(mask, dtype=bool)
    _, b0 = ndimage.label(m, structure=FG8)
    # Pad with background so every border-touching background region joins one component.
    _, nbg = ndimage.label(~np.pad(m, 1), structure=BG4)
    return b0, nbg - 1


def size_filter(mask, min_size):
    """Drop foreground components and fill holes smaller than min_size pixels."""
    m = np.asarray(mask, dtype=bool).copy()
    if min_size <= 1:
        return m
    lab, n = ndimage.label(m, structure=FG8)
    if n:
        sizes = np.bincount(lab.ravel())
        small = sizes < min_size
        small[0] = False
        m[small[lab]] = False
    bg, nbg = ndimage.label(~np.pad(m, 1), structure=BG4)
    if nbg > 1:
        sizes = np.bincount(bg.ravel())
        outer = bg[0, 0]
        small = sizes < min_size
        small[0] = False
        small[outer] = False
        m[small[bg[1:-1, 1:-1]]] = True
    return m


def betti3d(vol):
    """(beta0, beta2) of a 3D binary volume: 26-connected objects, 6-connected cavities."""
    v = np.asarray(vol, dtype=bool)
    _, b0 = ndimage.label(v, structure=np.ones((3, 3, 3), dtype=bool))
    _, nbg = ndimage.label(~np.pad(v, 1), structure=ndimage.generate_binary_structure(3, 1))
    return b0, nbg - 1
