"""Exact integer polygon core (C++17 header-only, nanobind bindings).

Shapes are passed as two int32 arrays: ``xy`` (N, 2) vertices on the G x G grid
and ``off`` (R+1,) ring offsets, ring r = xy[off[r]:off[r+1]]. ``pack`` builds
them from a list of rings. All kernels release the GIL.

Conventions: see cpp/polycore/polycore.hpp (identical to phase1_lab/common/oracle.py).
"""
from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

from ._polycore import (  # noqa: F401
    DEFECTS,
    MAX_RINGS_VALIDATE,
    hardware_threads,
    rasterize as _rasterize,
    score as _score,
    score_batch as _score_batch,
    signed_area2 as _signed_area2,
    validate as _validate,
)

__all__ = [
    "DEFECTS", "SELF_TOUCHING", "pack", "rasterize", "validate", "signed_area2",
    "ring_labels", "score", "score_batch", "hardware_threads",
]

#: Sub-classes whose union is oracle.ring_problems' "self_touching".
SELF_TOUCHING = ("vertex_pinch", "collinear_overlap", "spike")
_IDX = {name: i for i, name in enumerate(DEFECTS)}


def pack(rings: Iterable[Sequence[Sequence[int]]]) -> tuple[np.ndarray, np.ndarray]:
    """List of rings of (x, y) ints -> (xy int32 (N, 2), off int32 (R+1,)).

    Raises OverflowError / ValueError for values that do not fit int32 or malformed rings.
    """
    rings = [np.asarray(r, dtype=np.int64).reshape(-1, 2) if len(r) else np.zeros((0, 2), np.int64)
             for r in rings]
    off = np.zeros(len(rings) + 1, np.int64)
    if rings:
        off[1:] = np.cumsum([len(r) for r in rings])
    xy = np.concatenate(rings) if rings else np.zeros((0, 2), np.int64)
    info = np.iinfo(np.int32)
    if xy.size and (xy.min() < info.min or xy.max() > info.max):
        raise OverflowError("coordinates must fit in int32")
    if off[-1] > info.max:
        raise OverflowError("too many vertices")
    return np.ascontiguousarray(xy, np.int32), off.astype(np.int32)


def as_i32(a) -> np.ndarray:
    """Contiguous int32 copy/view of an integer array; refuses floats and out-of-range values
    (a plain astype would silently truncate or wrap)."""
    a = np.asarray(a)
    if a.size == 0:
        return np.ascontiguousarray(a, np.int32)
    if a.dtype.kind not in "iu":
        raise TypeError(f"expected an integer array, got {a.dtype}")
    if a.dtype != np.int32:
        info = np.iinfo(np.int32)
        if a.min() < info.min or a.max() > info.max:
            raise OverflowError("values must fit in int32")
    return np.ascontiguousarray(a, np.int32)


def _as_arrays(shape):
    if isinstance(shape, tuple) and len(shape) == 2 and isinstance(shape[0], np.ndarray):
        return as_i32(shape[0]), as_i32(shape[1])
    return pack(shape)


def rasterize(shape, W: int, H: int, G: int = 1000, winding: bool = False) -> dict:
    """shape: list of rings or (xy, off). Returns dict with 'evenodd', 'nonzero' bool (H, W),
    'ring_wmin', 'ring_wmax' (per ring), 'rings_mixed_sign', 'px_multi', areas, optional 'winding'."""
    xy, off = _as_arrays(shape)
    return _rasterize(xy, off, W, H, G, winding)


def validate(shape, roles=None, G: int = 1000, method: str = "sweep") -> dict:
    """Graded legality; see DEFECTS. roles: None (all outer) or per-ring -1 / outer index."""
    xy, off = _as_arrays(shape)
    if roles is not None:
        roles = as_i32(roles)
    return _validate(xy, off, roles, G, method)


def signed_area2(ring) -> int:
    """Twice the signed area of one ring (exact); > 0 = clockwise on screen (y down)."""
    a = np.asarray(ring, dtype=np.int64).reshape(-1, 2)
    if a.size and (a.min() < np.iinfo(np.int32).min or a.max() > np.iinfo(np.int32).max):
        raise OverflowError("coordinates must fit in int32")
    return _signed_area2(np.ascontiguousarray(a, np.int32))


def ring_labels(counts_row) -> list[str]:
    """Map one row of validate()['ring'] to oracle.ring_problems-style labels (self part only)."""
    c = np.asarray(counts_row)
    if c[_IDX["fewer_than_3_points"]]:
        return ["fewer_than_3_points"]
    out = []
    if c[_IDX["self_crossing"]]:
        out.append("self_crossing")
    if any(c[_IDX[k]] for k in SELF_TOUCHING):
        out.append("self_touching")
    if c[_IDX["zero_area"]]:
        out.append("zero_area")
    return sorted(out)


def score(xy, off, W, H, roles=None, gt=None, G: int = 1000) -> dict:
    """One rollout: rasterize (+ IoU counts vs gt) and validate, GIL released.
    xy/off/roles must already be int32 (no conversion on this hot path); gt uint8 (H, W)."""
    return _score(xy, off, W, H, roles, gt, G)


def score_batch(xys, offs, Ws, Hs, roles=None, gts=None, G: int = 1000, threads: int = 0) -> dict:
    """Many rollouts on an internal C++ thread pool (threads <= 1 runs inline)."""
    return _score_batch(list(xys), list(offs), list(Ws), list(Hs),
                        None if roles is None else list(roles),
                        None if gts is None else list(gts), G, threads)
