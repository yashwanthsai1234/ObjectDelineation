"""Exact integer polygon core: rasterize (even-odd and nonzero), graded validity, signed area.

Python face of the C++17 header ``cpp/polycore/polycore.hpp`` (nanobind module ``objdel._polycore``).
All kernels are exact integer code and run with the GIL released, so Python threads scale.

A *shape* is a list of rings, each a sequence of integer (x, y) vertices on the G x G grid (G = 1000)
over a W x H crop, or an already packed ``(xy, off)`` pair from :func:`pack`. Pixel (c, r) is inside
iff its centre is inside under the half-open crossing rule (top-left ties); see the header for the
full spec. Outer rings are clockwise on screen (positive :func:`signed_area2`), holes counter-clockwise.

Caps: at most ``MAX_RINGS`` rings per shape, ``MAX_RING_VERTICES`` vertices per ring and
1 <= W, H <= ``MAX_DIM``; anything beyond raises ``ValueError``.

References:
    https://nanobind.readthedocs.io/en/latest/functions.html#gil-management
    https://scikit-build-core.readthedocs.io/en/latest/configuration/editable.html
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from objdel._polycore import DEFECTS, MAX_DIM, MAX_RING_VERTICES, MAX_RINGS
from objdel._polycore import rasterize as _rasterize
from objdel._polycore import score_batch as _score_batch
from objdel._polycore import signed_area2 as _signed_area2
from objdel._polycore import validate as _validate

__all__ = [
    "DEFECTS", "MAX_DIM", "MAX_RINGS", "MAX_RING_VERTICES",
    "pack", "rasterize", "validate", "signed_area2", "score_batch",
]

_I32 = np.iinfo(np.int32)


def _as_i32(a, ndim: int) -> np.ndarray:
    """Contiguous int32 array of shape (k, 2) (ndim 2: vertices) or (k,) (ndim 1: offsets, roles). Refuses other
    shapes (ValueError; nothing is reshaped), non-integers (TypeError) and values outside int32 (OverflowError:
    astype would wrap). An empty sequence is an empty array of the right shape."""
    a = np.asarray(a)
    if a.shape == (0,) or (ndim == 2 and a.shape == (0, 2)):  # [] is an empty ring or list
        return np.zeros((0, 2) if ndim == 2 else (0,), np.int32)
    if a.ndim != ndim or (ndim == 2 and a.shape[1] != 2):
        raise ValueError(f"expected {'a (k, 2) ring' if ndim == 2 else 'a 1-D array'}, got shape {a.shape}")
    if a.dtype.kind not in "iu":
        raise TypeError(f"expected integers, got {a.dtype}")
    if a.dtype != np.int32 and a.size and (a.min() < _I32.min or a.max() > _I32.max):
        raise OverflowError("values must fit in int32")
    return np.ascontiguousarray(a, np.int32)


def pack(rings: Sequence) -> tuple[np.ndarray, np.ndarray]:
    """List of rings of (x, y) ints -> (xy int32 (N, 2), off int32 (R+1,)); ring r = xy[off[r]:off[r+1]]."""
    parts = [_as_i32(r, 2) for r in rings]
    off = np.zeros(len(parts) + 1, np.int32)
    off[1:] = np.cumsum([len(p) for p in parts])
    return (np.concatenate(parts) if parts else np.zeros((0, 2), np.int32)), off


def _packed(shape) -> tuple[np.ndarray, np.ndarray]:
    if isinstance(shape, tuple) and len(shape) == 2 and isinstance(shape[0], np.ndarray):
        return _as_i32(shape[0], 2), _as_i32(shape[1], 1)
    return pack(shape)


def _roles(roles) -> np.ndarray | None:
    return None if roles is None else _as_i32(roles, 1)


def _mask(gt) -> np.ndarray | None:
    if gt is None:
        return None
    gt = np.ascontiguousarray(gt)
    if gt.dtype == np.bool_:
        gt = gt.view(np.uint8)  # zero-copy for bool and uint8 masks
    return gt if gt.dtype == np.uint8 else (gt != 0).view(np.uint8)


def rasterize(shape, W: int, H: int, G: int = 1000) -> dict:
    """Exact masks of all rings together.

    Returns ``evenodd`` and ``nonzero`` bool (H, W); ``ring_wmin``/``ring_wmax`` int32 (R,), the range of
    each ring's own winding (both signs = bow-tie or figure-8); ``area_evenodd``, ``area_nonzero``,
    ``px_multi`` (pixels with |winding| >= 2) and ``rings_mixed_sign``.
    """
    xy, off = _packed(shape)
    return _rasterize(xy, off, W, H, G)


def validate(shape, roles=None, G: int = 1000) -> dict:
    """Graded legality; never repairs.

    roles: None (every ring is an outer ring) or one int per ring, -1 for an outer ring and p for a hole
    of outer ring p. Returns ``total`` int64 (13,) and ``per_ring`` int32 (R, 13) defect counts in
    :data:`DEFECTS` order (crossing/touching/overlap classes count edge pairs), plus ``area2``,
    ``n_clean``, ``pair_cross`` and ``pair_touch`` (R, R).
    """
    xy, off = _packed(shape)
    return _validate(xy, off, _roles(roles), G)


def signed_area2(ring) -> int:
    """Twice the signed shoelace area of one ring (exact Python int); > 0 means clockwise on screen."""
    hi, lo = _signed_area2(_as_i32(ring, 2))
    return (hi << 64) + lo


def score_batch(shapes: Sequence, Ws: Sequence[int], Hs: Sequence[int], roles=None, gts=None,
                G: int = 1000, threads: int = 0) -> dict:
    """Rasterize statistics (+ intersections with an optional gt mask) and validity of many shapes.

    Runs on an internal C++ pool of ``threads`` threads (0 or 1 = inline). Results are written by index,
    so they do not depend on the thread count. Returns per-shape arrays ``area_evenodd``,
    ``area_nonzero``, ``px_multi``, ``rings_mixed_sign``, ``gt_area``, ``inter_evenodd``,
    ``inter_nonzero`` and ``defects`` (n, 13) (= ``validate(...)["total"]``).
    """
    packed = [_packed(s) for s in shapes]
    return _score_batch([p[0] for p in packed], [p[1] for p in packed], list(Ws), list(Hs),
                        None if roles is None else [_roles(r) for r in roles],
                        None if gts is None else [_mask(g) for g in gts], G, threads)
