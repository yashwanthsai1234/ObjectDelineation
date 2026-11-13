"""Canonical form of polygon objects, and their masks and legality from the CUDA checker.

Canonical form (screen coordinates, x right, y down; `signed_area2` > 0 means clockwise on screen):
- an outer ring runs clockwise (signed area > 0), a hole counter-clockwise (< 0); a zero-area ring has no
  orientation, so both directions are candidates;
- a ring starts at its top-most, then left-most vertex (minimum (y, x)); if that vertex occurs more than once,
  every occurrence is a candidate start;
- among the candidates the ring with the smallest (y, x) vertex sequence wins;
- holes of a part, and parts of an object, are sorted by the (y, x) vertex sequence of their canonical ring
  (parts: outer ring, then the list of holes). For legal shapes, whose starts differ, this is the row-major order
  of the start vertices.
This is the convention of the repo's polygonize. Clockwise-from-the-top-left order follows
PolyFormer and VectorLLM (which pick the vertex nearest the image corner; we use row-major order because it is
exact on the integer grid and has no distance ties). Numerically it is the GeoJSON right-hand rule (exterior
positive shoelace area, holes negative); on a y-down screen that exterior appears clockwise.

References:
    https://arxiv.org/abs/2302.07387 (PolyFormer §3: clockwise, start nearest the top-left, polygons ordered by start)
    https://arxiv.org/abs/2507.04664 (VectorLLM: clockwise from the top-left corner; a fixed start helps SFT)
    https://datatracker.ietf.org/doc/html/rfc7946#section-3.1.6 (right-hand rule: exterior CCW, holes CW, y up)
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

import numpy as np

from objdel.geom import engine as _gpucheck
from objdel.geom.text import GRID, Objects, Part

_I32 = np.iinfo(np.int32)


class Geometry(NamedTuple):
    """One rasterized object: even-odd bool mask, legality and the non-zero defect counts."""
    mask: np.ndarray  # bool (H, W), even-odd fill sampled at pixel centres
    legal: bool
    defects: dict[str, int]  # non-zero _gpucheck.DEFECTS counts


def signed_area2(ring) -> int:
    """Twice the signed shoelace area of one ring of integer (x, y) vertices (exact Python int); > 0 means clockwise
    on screen (y down)."""
    x, y = np.asarray(ring).reshape(-1, 2).T.tolist()
    return sum(x[i - 1] * y[i] - x[i] * y[i - 1] for i in range(len(x)))


def _key(ring: np.ndarray) -> list:
    return ring[:, ::-1].tolist()  # (y, x) per vertex: row-major order


def _canonical_ring(ring, hole: bool) -> np.ndarray:
    r = np.asarray(ring, np.int64).reshape(-1, 2)
    a2 = signed_area2(r)
    if a2 == 0:
        directions = (r, r[::-1])
    else:
        directions = (r[::-1],) if (a2 > 0) == hole else (r,)
    candidates = []
    for d in directions:
        top = d[:, 1] == d[:, 1].min()
        starts = np.flatnonzero(top & (d[:, 0] == d[top, 0].min()))
        candidates += [np.roll(d, -int(i), axis=0) for i in starts]
    return min(candidates, key=_key)


def canonicalize(objs: Objects) -> Objects:
    """Canonical copy of objects (inputs are not modified); see the module docstring for the order."""
    out = {}
    for name, parts in objs.items():
        canon = [(_canonical_ring(o, False), sorted((_canonical_ring(h, True) for h in hs), key=_key))
                 for o, hs in parts]
        out[name] = sorted(canon, key=lambda p: (_key(p[0]), [_key(h) for h in p[1]]))
    return out


def to_mask(parts: list[Part], W: int, H: int) -> Geometry:
    """Parts of one object on a W x H crop -> even-odd mask and legality (each hole checked against its outer), on
    the CUDA checker (to_masks of one object).

    No parts is a legal empty mask: "no object" is a well-formed answer, judged by the reward's topology terms."""
    return to_masks([(parts, W, H)])[0]


def unpack_rows(words: np.ndarray, H: int, W: int) -> np.ndarray:
    """The checker's packed mask of one object -> bool (H, W). Each row is ceil(W / 32) uint32 words; pixel x is bit
    x % 32 (LSB first) of word x // 32, and the bits past W in a row's last word are ignored."""
    rows = np.ascontiguousarray(words, "<u4").view(np.uint8).reshape(H, -1)
    return np.ascontiguousarray(np.unpackbits(rows, axis=1, bitorder="little")[:, :W]).view(bool)


def _ring(r) -> np.ndarray:
    """One ring as a contiguous int32 (k, 2) array. Refuses other shapes (ValueError; nothing is reshaped),
    non-integers (TypeError) and values outside int32 (OverflowError: astype would wrap); [] is an empty ring."""
    a = np.asarray(r)
    if a.shape in ((0,), (0, 2)):
        return np.zeros((0, 2), np.int32)
    if a.ndim != 2 or a.shape[1] != 2:
        raise ValueError(f"expected a (k, 2) ring, got shape {a.shape}")
    if a.dtype.kind not in "iu":
        raise TypeError(f"expected integers, got {a.dtype}")
    if a.dtype != np.int32 and (a.min() < _I32.min or a.max() > _I32.max):
        raise OverflowError("values must fit in int32")
    return np.ascontiguousarray(a, np.int32)


def pack_rings(rings: Sequence) -> tuple[np.ndarray, np.ndarray]:
    """Rings of integer (x, y) vertices -> (xy int32 (N, 2), ring_start int64 (R + 1,)); ring r is
    xy[ring_start[r]:ring_start[r + 1]]."""
    parts = [_ring(r) for r in rings]
    start = np.zeros(len(parts) + 1, np.int64)
    start[1:] = np.cumsum([len(p) for p in parts])
    return (np.concatenate(parts) if parts else np.zeros((0, 2), np.int32)), start


def pack_objects(objs: Sequence[tuple[list[Part], int, int]]) -> tuple[np.ndarray, ...]:
    """Objects, each (parts, W, H), as the CUDA checker's packed arrays (xy (P, 2) int32, ring_start, shape_start
    int64, roles, W, H int32): one shape per object, its outer rings and holes in order; a hole's role is its outer
    ring's index in the object."""
    rings, roles, shape_start = [], [], [0]
    for parts, _, _ in objs:
        first = len(rings)
        for outer, holes in parts:
            roles += [-1] + [len(rings) - first] * len(holes)
            rings += [outer, *holes]
        shape_start.append(len(rings))
    xy, ring_start = pack_rings(rings)
    dims = np.array([(W, H) for _, W, H in objs], np.int32).reshape(-1, 2)
    return (xy, ring_start, np.array(shape_start, np.int64), np.array(roles, np.int32),
            np.ascontiguousarray(dims[:, 0]), np.ascontiguousarray(dims[:, 1]))


def defects_of(row) -> dict[str, int]:
    """{defect: count} of the nonzero counts of one engine result row (`r["defects"][i]`), in DEFECTS order."""
    return {d: int(c) for d, c in zip(_gpucheck.DEFECTS, row.tolist(), strict=True) if c}


def to_masks(objs: Sequence[tuple[list[Part], int, int]]) -> list[Geometry]:
    """to_mask of many objects, each (parts, W, H), on the CUDA checker: one objdel.geom.engine.check_polygons call for
    all of them, which returns the even-odd masks bit-packed (unpack_rows). On a CUDA box it runs on the GPU (`make
    cuda`; policy.chain builds it before any scoring stage); elsewhere the same source is compiled for the CPU. Masks
    and defects equal those of the archived CPU checker, polycore (recorded in tests/goldens/polycore.npz:
    tests/verify/test_gpu_scoring.py, tests/geom/test_checker_host.py). The reward never calls this: it scores on the
    device (objdel.verify.reward.score.score_many); targets, probes and the exam's item-level metrics read these
    masks."""
    r = _gpucheck.check_polygons(*pack_objects(objs), GRID)
    out = []
    for i, (_, W, H) in enumerate(objs):
        defects = defects_of(r["defects"][i])
        mask = unpack_rows(r["evenodd"][r["word_start"][i]:r["word_start"][i + 1]], H, W)
        out.append(Geometry(mask, not defects, defects))
    return out
