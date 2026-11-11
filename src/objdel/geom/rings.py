"""Mask -> rectilinear rings on the CUDA engine (objdel.geom.engine.check_masks: nvcc build on the box, CPU build of the
same source on the Mac). polygonize is trace of one mask.

The rings are the cracks between pixels, clockwise on screen (y down) around foreground, so outer rings have positive
shoelace area and holes negative; at a saddle vertex the walk turns left (8-connected foreground, 4-connected
background). Rings come in the order of their first edge, each from its minimum edge key (y * (W + 1) + x) * 4 + d.
A ring's part is the raster rank (by first pixel) of the 8-connected foreground component it bounds: its own for an
outer ring, the surrounding one for a hole. One engine call also gives each mask's bit-quad counts and Betti numbers
from the same run labeller (so certify takes b0, b1 from an independent one).
The numpy references are tests/oracles/polygonize_oracle.py and bitquad_oracle.py; the engine equals them on all
37,093 LIDC crops and the fuzz families (tests/experiments/test_rings.py).
"""
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from objdel.geom import engine as _gpucheck
from objdel.geom.topology import as_mask, pack_masks


@dataclass(frozen=True)
class Ring:
    """One rectilinear ring of a mask: corner vertices in pixel coordinates, signed area and part rank."""
    pts: np.ndarray  # (k, 2) int64 corner vertices (x, y), canonical start first, offset applied
    area: int  # signed: > 0 outer, < 0 hole
    part: int  # raster rank of the bounded 8-component within its mask

    @property
    def hole(self):
        """True for a hole ring (negative signed area)."""
        return self.area < 0


@dataclass(frozen=True)
class MaskRings:
    """Rings, bit-quad counts and Betti numbers of one mask, as returned by `trace`."""
    rings: list[Ring]
    quads: tuple[int, int, int, int, int]  # bit-quad counts n1, n2, nD, n3, n4
    b0: int  # 8-connected foreground components
    b1: int  # bounded 4-connected background components


def trace(masks: Sequence, offsets: Sequence[tuple[int, int]] | None = None) -> list[MaskRings]:
    """Rings, bit-quads and Betti numbers of 2-D masks (any nonzero value is foreground, sides 0..4096; larger raise
    ValueError) in one engine call. offsets[i] (x, y) is added to mask i's corners, e.g. after cropping."""
    ms = [as_mask(m) for m in masks]
    if not ms:
        return []
    offs = np.zeros((len(ms), 2), np.int64) if offsets is None else np.asarray(offsets, np.int64)
    if offs.shape != (len(ms), 2):
        raise ValueError(f"offsets shape {offs.shape}, need ({len(ms)}, 2)")
    out = _gpucheck.check_masks(*pack_masks(ms))
    rs, ps, pts = out["ring_start"], out["pts_start"], out["pts"].astype(np.int64)
    area, part = out["ring_area"].tolist(), out["ring_part"].tolist()
    quads, betti = out["quads"].tolist(), out["betti"].tolist()
    return [MaskRings([Ring(pts[ps[r]:ps[r + 1]] + offs[i], area[r], part[r]) for r in range(rs[i], rs[i + 1])],
                      tuple(quads[i]), *betti[i]) for i in range(len(ms))]


def polygonize(mask, offset=(0, 0)) -> list[Ring]:
    """Rings of one 2-D mask; offset (x, y) is added to the corners."""
    return trace([mask], [offset])[0].rings
