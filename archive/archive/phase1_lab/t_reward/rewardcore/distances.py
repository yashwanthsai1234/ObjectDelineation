"""Boundaries, Euclidean distance transforms and inter-rater surface distances (pixel units)."""
from __future__ import annotations

from itertools import combinations

import cv2
import numpy as np

CROSS = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], np.uint8)
BIG = np.float32(1e6)


def boundary(m: np.ndarray) -> np.ndarray:
    """Contour pixels: foreground pixels with a 4-adjacent background pixel (outside the crop = background).

    Under 8/4 connectivity these pixels form the 8-connected contour; their count approximates
    the Euclidean perimeter in px.
    """
    m8 = np.ascontiguousarray(m, dtype=np.uint8)
    er = cv2.erode(m8, CROSS, borderType=cv2.BORDER_CONSTANT, borderValue=0)
    return (m8 > er)


def edt(nonzero: np.ndarray) -> np.ndarray:
    """Exact Euclidean distance from every nonzero pixel to the nearest zero pixel (float32).

    Zero pixels get 0. If there is no zero pixel at all, every pixel gets BIG.
    """
    src = np.ascontiguousarray(nonzero, dtype=np.uint8)
    if src.all():
        return np.full(src.shape, BIG, np.float32)
    return cv2.distanceTransform(src, cv2.DIST_L2, cv2.DIST_MASK_PRECISE)


def directed_surface_distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Distance from each contour pixel of a to the nearest contour pixel of b."""
    ba, bb = boundary(a), boundary(b)
    if not ba.any():
        return np.zeros(0, np.float32)
    if not bb.any():
        return np.full(int(ba.sum()), BIG, np.float32)
    return edt(~bb)[ba]


def symmetric_surface_distances(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.concatenate([directed_surface_distances(a, b), directed_surface_distances(b, a)])


def pair_p95(a: np.ndarray, b: np.ndarray, q: float = 95.0) -> float:
    d = symmetric_surface_distances(a, b)
    return float(np.percentile(d, q)) if d.size else 0.0


def crop_to_union(stack: np.ndarray, pad: int = 1) -> np.ndarray:
    """Crop a rater stack to the union's bounding box plus `pad` background pixels.

    Contour distances, and distances to the union or to the complement of the core, are
    unchanged: every relevant zero pixel lies inside the crop (clamping a point to the box
    never increases its distance to a point inside the box).
    """
    u = stack.any(0)
    rows, cols = np.flatnonzero(u.any(1)), np.flatnonzero(u.any(0))
    if not rows.size:
        return stack
    y0, y1 = max(0, rows[0] - pad), min(u.shape[0], rows[-1] + 1 + pad)
    x0, x1 = max(0, cols[0] - pad), min(u.shape[1], cols[-1] + 1 + pad)
    out = np.zeros((len(stack), y1 - y0 + 2 * pad, x1 - x0 + 2 * pad), bool)  # pad even at crop borders
    out[:, pad:pad + y1 - y0, pad:pad + x1 - x0] = stack[:, y0:y1, x0:x1]
    return out


def rater_tau(stack: np.ndarray, floor: float = 1.0, q: float = 95.0) -> tuple[float, float]:
    """(tau, raw) with raw = q-th percentile of the pooled symmetric contour distances over all rater pairs.

    tau = max(floor, raw). With fewer than two raters raw is NaN and tau = floor.
    """
    if len(stack) < 2:
        return float(floor), float("nan")
    stack = crop_to_union(np.asarray(stack, bool))
    d = np.concatenate([symmetric_surface_distances(a, b) for a, b in combinations(stack, 2)])
    raw = float(np.percentile(d, q)) if d.size else 0.0
    return max(float(floor), raw), raw
