"""One call that measures a mask three independent ways and checks they agree.

The three measurements: Gray's bit-quad counts (checker.bitquad), the labeled topology (checker.topology, OpenCV
Spaghetti labeling) and the rectilinear rings (checker.polygonize). They must agree on the Euler number
chi8 = b0 - b1 (Euler-Poincare formula for 8-connected foreground, 4-connected background), on the number of outer
rings and holes, on the area and on the ring vertices.

References:
- S. B. Gray, "Local properties of binary images in two dimensions", IEEE Trans. Computers C-20(5):551-561, 1971.
  https://doi.org/10.1109/T-C.1971.223289
- T. Y. Kong, A. Rosenfeld, "Digital topology: introduction and survey", CVGIP 48(3):357-393, 1989 (the (8, 4)
  pairing and the Euler number). https://doi.org/10.1016/0734-189X(89)90147-3
- Sources of the three methods: the References lists of checker.bitquad, checker.topology and checker.polygonize.
"""
from dataclasses import dataclass

from .bitquad import bitquads
from .polygonize import polygonize
from .topology import as_mask, bbox, betti, size_filter


@dataclass
class Report:
    b0: int
    b1: int
    area: int
    perimeter: int
    corners: int
    chi8: int
    certified: bool
    failures: list


def crop(mask):
    """Crop to the foreground bounding box; returns the crop and its offset (x, y). Topology is unchanged; functions
    pad internally. Also used by the lab target builders (archive/phase1_lab)."""
    m = as_mask(mask)
    x, y, w, h = bbox(m)
    return m[y:y + h, x:x + w], (x, y)


def certify(mask, min_size=0):
    m, offset = crop(size_filter(mask, min_size) if min_size > 1 else mask)
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
    return Report(b0, b1, q.area, q.perimeter, q.corners, q.chi8, not fails, fails)
