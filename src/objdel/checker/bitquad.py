"""Bit-quad (2x2 pixel pattern) counts and Gray's 1971 invariants.

Each 2x2 window over the zero-padded mask gets a code TL=1, TR=2, BL=4, BR=8.
Q1 = one pixel set, Q2 = two adjacent, QD = two diagonal (codes 6, 9),
Q3 = three, Q4 = four. With 8-connected foreground and 4-connected background:

    chi8      = (n1 - n3 - 2 nD) / 4
    area      = (n1 + 2 n2 + 3 n3 + 4 n4 + 2 nD) / 4
    perimeter = n1 + n2 + n3 + 2 nD      (unit boundary edges)
    corners   = n1 + n3 + 2 nD           (rectilinear ring vertices)

chi8, area and perimeter are Gray's bit-quad formulas [Gray 1971; Pratt]. corners counts the windows where the
boundary turns: one for Q1 and Q3, two for a diagonal QD, none for Q2 and Q4.

References:
- S. B. Gray, "Local properties of binary images in two dimensions", IEEE Trans. Computers C-20(5):551-561, 1971.
  https://doi.org/10.1109/T-C.1971.223289
- W. K. Pratt, "Digital Image Processing", 3rd ed., Wiley, 2001, pp. 593-595 (bit quads; the reference of Octave's
  bweuler): https://octave.sourceforge.io/image/function/bweuler.html
- B. Yao, L. He, S. Kang, Y. Chao, X. Zhao, "A novel bit-quad-based Euler number computing algorithm",
  SpringerPlus 4:735, 2015 (restates Gray's algorithm). https://doi.org/10.1186/s40064-015-1511-8
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


def bitquads(mask):
    m = np.pad(np.asarray(mask, dtype=bool), 1).astype(np.uint8)
    codes = m[:-1, :-1] | (m[:-1, 1:] << 1) | (m[1:, :-1] << 2) | (m[1:, 1:] << 3)  # (H+1, W+1) windows
    hist = np.bincount(codes.ravel(), minlength=16)
    pick = lambda codes: int(hist[list(codes)].sum())
    return BitQuads(pick(Q1), pick(Q2), pick(QD), pick(Q3), pick(Q4))
