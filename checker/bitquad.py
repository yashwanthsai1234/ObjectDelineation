"""Bit-quad (2x2 pixel pattern) counts and Gray's 1971 invariants.

Each 2x2 window over the zero-padded mask gets a code TL=1, TR=2, BL=4, BR=8.
Q1 = one pixel set, Q2 = two adjacent, QD = two diagonal (codes 6, 9),
Q3 = three, Q4 = four. With 8-connected foreground and 4-connected background:

    chi8      = (n1 - n3 - 2 nD) / 4
    area      = (n1 + 2 n2 + 3 n3 + 4 n4 + 2 nD) / 4
    perimeter = n1 + n2 + n3 + 2 nD      (unit boundary edges)
    corners   = n1 + n3 + 2 nD           (rectilinear ring vertices)
"""
from dataclasses import dataclass

import numpy as np

Q1 = (1, 2, 4, 8)
Q2 = (3, 5, 10, 12)
QD = (6, 9)
Q3 = (7, 11, 13, 14)
Q4 = (15,)


@dataclass(frozen=True)
class BitQuads:
    hist: tuple
    n1: int
    n2: int
    nD: int
    n3: int
    n4: int

    @property
    def chi8(self):
        v = self.n1 - self.n3 - 2 * self.nD
        assert v % 4 == 0
        return v // 4

    @property
    def chi4(self):
        v = self.n1 - self.n3 + 2 * self.nD
        assert v % 4 == 0
        return v // 4

    @property
    def area(self):
        v = self.n1 + 2 * self.n2 + 3 * self.n3 + 4 * self.n4 + 2 * self.nD
        assert v % 4 == 0
        return v // 4

    @property
    def perimeter(self):
        return self.n1 + self.n2 + self.n3 + 2 * self.nD

    @property
    def corners(self):
        return self.n1 + self.n3 + 2 * self.nD


def quad_codes(mask):
    """16-valued code for every 2x2 window of the padded mask, shape (H+1, W+1)."""
    m = np.pad(np.asarray(mask, dtype=bool), 1).astype(np.uint8)
    return m[:-1, :-1] | (m[:-1, 1:] << 1) | (m[1:, :-1] << 2) | (m[1:, 1:] << 3)


def bitquads(mask):
    hist = np.bincount(quad_codes(mask).ravel(), minlength=16)
    pick = lambda codes: int(hist[list(codes)].sum())
    return BitQuads(tuple(int(h) for h in hist), pick(Q1), pick(Q2), pick(QD), pick(Q3), pick(Q4))
