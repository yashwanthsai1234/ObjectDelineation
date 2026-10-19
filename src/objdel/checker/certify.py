"""One call that measures a mask three independent ways and checks they agree."""
from dataclasses import dataclass, field

import numpy as np

from .bitquad import bitquads
from .polygonize import mask_hash, polygonize
from .topology import betti, size_filter


@dataclass
class Report:
    b0: int
    b1: int
    area: int
    perimeter: int
    corners: int
    chi8: int
    hash: str
    certified: bool
    failures: list = field(default_factory=list)
    rings: list = field(default_factory=list, repr=False)


def crop(mask):
    """Crop to the foreground bounding box. Topology is unchanged; functions pad internally."""
    m = np.asarray(mask, dtype=bool)
    rows, cols = np.any(m, axis=1), np.any(m, axis=0)
    if not rows.any():
        return m[:0, :0], (0, 0)
    r0, r1 = np.argmax(rows), len(rows) - np.argmax(rows[::-1])
    c0, c1 = np.argmax(cols), len(cols) - np.argmax(cols[::-1])
    return m[r0:r1, c0:c1], (int(c0), int(r0))


def certify(mask, min_size=0, keep_rings=False):
    m = np.asarray(mask, dtype=bool)
    if min_size > 1:
        m = size_filter(m, min_size)
    m, offset = crop(m)
    q = bitquads(m)
    b0, b1 = betti(m)
    rings = polygonize(m, offset)
    outer = sum(1 for r in rings if not r.hole)
    holes = len(rings) - outer
    fails = []
    if b0 - b1 != q.chi8:
        fails.append(f"euler: b0-b1={b0 - b1} chi8={q.chi8}")
    if (outer, holes) != (b0, b1):
        fails.append(f"rings: outer/holes={outer}/{holes} betti={b0}/{b1}")
    if sum(r.area for r in rings) != q.area or q.area != int(m.sum()):
        fails.append(f"area: rings={sum(r.area for r in rings)} quads={q.area} pixels={int(m.sum())}")
    if sum(len(r.pts) for r in rings) != q.corners:
        fails.append(f"corners: rings={sum(len(r.pts) for r in rings)} quads={q.corners}")
    return Report(b0, b1, q.area, q.perimeter, q.corners, q.chi8, mask_hash(rings),
                  not fails, fails, rings if keep_rings else [])
