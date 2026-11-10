"""Host reference of objdel.eval.metrics' item-level comparisons: the code metrics() ran on host masks before D7
moved its pixel work to the engine (objdel.eval.device), kept as the oracle of gc.eval_sets and as the subject of A5's
boundary-IoU and NSD reference tests (the official boundary-iou-api and google-deepmind/surface-distance formulas).
"""
import hashlib

import cv2
import numpy as np
from scipy import ndimage

from objdel.eval.select import vote_mask
from objdel.reward import score_mask

_D = 0.5 * np.sqrt(2.0)  # surfel length per 2 x 2 neighbourhood code (surface-distance lookup table, 1 px spacing)
_SURFEL = np.array([0, _D, _D, 1, _D, 1, 2 * _D, _D, _D, 2 * _D, 1, _D, 1, _D, _D, 0])
_CLASSES = ((1, 2, 4, 7, 8, 11, 13, 14), (3, 5, 10, 12), (6, 9))  # codes of length sqrt(2)/2, 1, sqrt(2)


def joint_box(masks: np.ndarray) -> tuple[slice, slice] | None:
    """Bounding box (row slice, column slice) of the foreground of all masks (n, H, W); None if there is none."""
    rows, cols = np.flatnonzero(masks.any((0, 2))), np.flatnonzero(masks.any((0, 1)))
    return (slice(rows[0], rows[-1] + 1), slice(cols[0], cols[-1] + 1)) if rows.size else None


def pairwise_iou(masks: np.ndarray) -> np.ndarray:
    """IoU of every pair of masks (n, H, W); two empty masks: 1."""
    x = masks.reshape(len(masks), -1).astype(np.float64)
    inter = x @ x.T
    area = np.diag(inter)
    union = area[:, None] + area[None, :] - inter
    return np.where(union > 0, inter / np.where(union > 0, union, 1.0), 1.0)


def band(m: np.ndarray, d: int) -> np.ndarray:
    """The inner d px band of a mask (H, W): the official mask_to_boundary (1 px zero pad, d erosions); bool."""
    m = np.asarray(m, bool)
    e = cv2.erode(np.pad(m, 1).view(np.uint8), np.ones((3, 3), np.uint8), iterations=d)[1:-1, 1:-1]
    return m & ~e.view(bool)


def band_iou(masks: np.ndarray, d: int) -> np.ndarray:
    """Boundary IoU of every pair of masks (n, H, W): the IoU of their bands, on the joint bounding box."""
    box = joint_box(masks)
    if box is None:
        return np.ones((len(masks), len(masks)))
    return pairwise_iou(np.stack([band(m[box], d) for m in masks]))


def boundary_iou(a: np.ndarray, b: np.ndarray, d: int) -> float:
    """Boundary IoU of two masks at band width d px."""
    return float(band_iou(np.stack([np.asarray(a, bool), np.asarray(b, bool)]), d)[0, 1])


def surface(m: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(2 x 2 code on the (H + 1, W + 1) corner grid, surfel length (0 off the surface), distance of every corner to
    the surface (inf for an empty mask))."""
    p = np.pad(np.asarray(m, bool), 1).astype(np.uint8)
    code = 8 * p[:-1, :-1] + 4 * p[:-1, 1:] + 2 * p[1:, :-1] + p[1:, 1:]
    border = (code != 0) & (code != 15)
    dist = ndimage.distance_transform_edt(~border) if border.any() else np.full(border.shape, np.inf)
    return code, _SURFEL[code], dist


def nsd(a: np.ndarray, b: np.ndarray, tau: float) -> float:
    """Normalized surface Dice of two masks at tolerance tau px (two empty: 1)."""
    (_, la, da), (_, lb, db) = surface(a), surface(b)
    total = la.sum() + lb.sum()
    return min(1.0, max(0.0, 1.0 - (la[db > tau].sum() + lb[da > tau].sum()) / total)) if total else 1.0


def dup_hash(format_ok, legal, raw: dict) -> str:
    """The pre-D7 duplicate key: SHA-256 of the flags and the packed raw masks of every object."""
    h = hashlib.sha256(bytes([int(format_ok), int(legal)]))
    for m in raw.values():
        h.update(np.packbits(m).tobytes())
    return h.hexdigest()


def counts(raw, drawn, raters, vote, thresholds, tau, d) -> dict:
    """gc.eval_sets' outputs of one set from host rasters: raw / drawn (K, H, W), raters (k, H, W), vote (H, W)."""
    K, k = len(raw), len(raters)
    allm = np.concatenate([drawn, raters, vote[None]]).astype(np.int64)
    flat = lambda x: np.asarray(x, np.int64).reshape(len(x), -1)
    rf, bf = flat(raw), flat([band(m, d) for m in allm])
    inter_raw = rf @ rf.T
    p, r = drawn.mean(0), raters.mean(0)
    soft = [[(p >= t).sum(), (r >= t).sum(), ((p >= t) & (r >= t)).sum()] for t in thresholds]
    surf = [surface(m) for m in np.concatenate([drawn, raters])]
    by_class = lambda code, sel: [int((np.isin(code, c) & sel).sum()) for c in _CLASSES]
    nsd_out = np.zeros((K, k, 2, 3), np.int64)
    for i in range(K):
        for j in range(k):
            (ci, _, di), (cj, _, dj) = surf[i], surf[K + j]
            nsd_out[i, j] = [by_class(ci, dj > tau), by_class(cj, di > tau)]
    return {"inter": flat(allm) @ flat(allm).T, "inter_raw": inter_raw, "area_raw": np.diag(inter_raw),
            "band_inter": bf @ bf.T, "soft": np.array(soft, np.int64),
            "surf": np.array([by_class(c, np.ones(c.shape, bool)) for c, _, _ in surf], np.int64), "nsd_out": nsd_out}


def vote_iou_tc(raw, legal, cache) -> float:
    """IoU_tc of the pixel vote of the legal raw masks (metrics' pixel_vote_iou_tc before D7)."""
    return float(score_mask(vote_mask(raw, legal), cache).iou_tc)
