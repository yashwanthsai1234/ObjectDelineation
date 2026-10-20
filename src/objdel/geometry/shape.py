"""Canonical form of polygon objects and their mask, both on top of objdel.checker.polycore (plan D1, D2, D10).

Canonical form (screen coordinates, x right, y down; `signed_area2` > 0 means clockwise on screen):
- an outer ring runs clockwise (signed area > 0), a hole counter-clockwise (< 0); a zero-area ring has no
  orientation, so both directions are candidates;
- a ring starts at its top-most, then left-most vertex (minimum (y, x)); if that vertex occurs more than once,
  every occurrence is a candidate start;
- among the candidates the ring with the smallest (y, x) vertex sequence wins;
- holes of a part, and parts of an object, are sorted by the (y, x) vertex sequence of their canonical ring
  (parts: outer ring, then the list of holes). For legal shapes, whose starts differ, this is the row-major order
  of the start vertices, as in the lab target builder.
This is the convention of the repo's polygonize and of the lab targets. Clockwise-from-the-top-left order follows
PolyFormer and VectorLLM (which pick the vertex nearest the image corner; we use row-major order because it is
exact on the integer grid and has no distance ties). Numerically it is the GeoJSON right-hand rule (exterior
positive shoelace area, holes negative); on a y-down screen that exterior appears clockwise.

References:
    https://arxiv.org/abs/2302.07387 (PolyFormer §3: clockwise, start nearest the top-left, polygons ordered by start)
    https://arxiv.org/abs/2507.04664 (VectorLLM: clockwise from the top-left corner; a fixed start helps SFT)
    https://datatracker.ietf.org/doc/html/rfc7946#section-3.1.6 (right-hand rule: exterior CCW, holes CW, y up)
    experiments/phase1/lab/t_targets/targets/builder.py (canon_ring, canonical_parts)
    experiments/phase1/EXECUTION.md §3 (polycore contract: pack, rasterize, validate with roles, DEFECTS)
"""
from __future__ import annotations

from typing import NamedTuple

import numpy as np

from objdel.checker import polycore

from .text import GRID, Objects, Part


class Geometry(NamedTuple):
    mask: np.ndarray  # bool (H, W), even-odd fill sampled at pixel centres
    legal: bool
    defects: dict[str, int]  # non-zero polycore.DEFECTS counts


def _key(ring: np.ndarray) -> list:
    return ring[:, ::-1].tolist()  # (y, x) per vertex: row-major order


def _canonical_ring(ring, hole: bool) -> np.ndarray:
    r = np.asarray(ring, np.int64).reshape(-1, 2)
    a2 = polycore.signed_area2(r)
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
    """Parts of one object on a W x H crop -> even-odd mask and legality (each hole checked against its outer).

    No parts is a legal empty mask: "no object" is a well-formed answer, judged by the reward's topology terms."""
    rings, roles = [], []
    for outer, holes in parts:
        roles += [-1] + [len(rings)] * len(holes)
        rings += [outer, *holes]
    shape = polycore.pack(rings)
    counts = polycore.validate(shape, roles=np.array(roles, np.int32), G=GRID)["total"]
    defects = {d: int(c) for d, c in zip(polycore.DEFECTS, counts, strict=True) if c}
    mask = np.asarray(polycore.rasterize(shape, W, H, G=GRID)["evenodd"], bool)
    return Geometry(mask, not defects, defects)
