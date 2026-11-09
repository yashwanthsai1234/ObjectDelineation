"""One call that measures a mask three independent ways and checks they agree.

The three measurements: Gray's bit-quad counts and the rectilinear rings, both from one CUDA engine call
(checker.rings.trace: the bit-quads count 2x2 windows of the packed words, the rings come from run labelling and
crack tracing), and the labeled topology (checker.topology: the engine's pixel union-find label chain, a second
labeller, so the ring count is never checked against the rings' own run labeller). They must agree on the Euler
number chi8 = b0 - b1 (Euler-Poincare formula for 8-connected foreground, 4-connected background), on the number of
outer rings and holes, on the area and on the ring vertices. From the bit-quads (n1, n2, nD, n3, n4) [Gray 1971; Pratt]:

    chi8      = (n1 - n3 - 2 nD) / 4
    area      = (n1 + 2 n2 + 3 n3 + 4 n4 + 2 nD) / 4
    perimeter = n1 + n2 + n3 + 2 nD      (unit boundary edges)
    corners   = n1 + n3 + 2 nD           (rectilinear ring vertices)

References:
- S. B. Gray, "Local properties of binary images in two dimensions", IEEE Trans. Computers C-20(5):551-561, 1971.
  https://doi.org/10.1109/T-C.1971.223289
- W. K. Pratt, "Digital Image Processing", 3rd ed., Wiley, 2001, pp. 593-595 (bit quads).
- T. Y. Kong, A. Rosenfeld, "Digital topology: introduction and survey", CVGIP 48(3):357-393, 1989 (the (8, 4)
  pairing and the Euler number). https://doi.org/10.1016/0734-189X(89)90147-3
- Sources of the three methods: the References lists of checker.topology and of the numpy references
  tests/oracles/bitquad_oracle.py and tests/oracles/polygonize_oracle.py.
"""
from dataclasses import dataclass

import numpy as np

from .rings import trace
from .topology import as_mask, bbox, filtered


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
    pad internally. Also used by the lab target builders (experiments/phase1/lab)."""
    m = as_mask(mask)
    x, y, w, h = bbox(m)
    return m[y:y + h, x:x + w], (x, y)


def certify_many(masks, min_size=0) -> list[Report]:
    """certify of each mask: the size filter and Betti numbers of all in one topology call, the rings and bit-quads of
    all in one trace call. min_size: one value, or one per mask (<= 1: no filter)."""
    sizes = [min_size] * len(masks) if np.ndim(min_size) == 0 else list(min_size)
    if len(sizes) != len(masks):
        raise ValueError(f"{len(sizes)} min_size values for {len(masks)} masks")
    fb = filtered(masks, [max(int(s), 0) for s in sizes])
    ms = [crop(f)[0] for f, _, _ in fb]
    return [_report(m, t, b0, b1) for m, t, (_, b0, b1) in zip(ms, trace(ms), fb, strict=True)]


def certify(mask, min_size=0) -> Report:
    return certify_many([mask], min_size)[0]


def _report(m, t, b0: int, b1: int) -> Report:
    n1, n2, nd, n3, n4 = t.quads
    chi4, area4 = n1 - n3 - 2 * nd, n1 + 2 * n2 + 3 * n3 + 4 * n4 + 2 * nd  # 4 chi8 and 4 area
    corners = n1 + n3 + 2 * nd
    outer = sum(1 for r in t.rings if not r.hole)
    holes = len(t.rings) - outer
    ring_area, ring_corners, pixels = sum(r.area for r in t.rings), sum(len(r.pts) for r in t.rings), int(m.sum())
    fails = []
    if chi4 != 4 * (b0 - b1):
        fails.append(f"euler: b0-b1={b0 - b1} 4*chi8={chi4}")
    if (outer, holes) != (b0, b1):
        fails.append(f"rings: outer/holes={outer}/{holes} betti={b0}/{b1}")
    if 4 * ring_area != area4 or area4 != 4 * pixels:
        fails.append(f"area: rings={ring_area} 4*quads={area4} pixels={pixels}")
    if ring_corners != corners:
        fails.append(f"corners: rings={ring_corners} quads={corners}")
    return Report(b0, b1, area4 // 4, n1 + n2 + n3 + 2 * nd, corners, chi4 // 4, not fails, fails)
