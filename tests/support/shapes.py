"""Boolean masks drawn on a grid."""
import numpy as np


def disc(h, w, cy, cx, r):
    """bool (h, w): pixels whose squared distance to (cy, cx) is at most r * r."""
    y, x = np.ogrid[:h, :w]
    return (y - cy) ** 2 + (x - cx) ** 2 <= r * r


def annulus(h, w, cy, cx, r, hole=0):
    """bool (h, w): the disc minus a hole of radius `hole`. hole=0 also drops the exact centre pixel (d2 > 0): the
    synthetic items and the stored fixtures were built with it, so it stays."""
    y, x = np.ogrid[:h, :w]
    d2 = (y - cy) ** 2 + (x - cx) ** 2
    return (d2 <= r * r) & (d2 > hole * hole)
