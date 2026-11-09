"""Mask -> rectilinear rings on the CUDA engine (objdel._gpucheck.check_masks: nvcc build on the box, CPU build of the
same source on the Mac). The polygonizer seam of the plan (B1.4a polygonize, D3 trace).

The rings are the cracks between pixels, clockwise on screen (y down) around foreground, so outer rings have positive
shoelace area and holes negative; at a saddle vertex the walk turns left (8-connected foreground, 4-connected
background). Rings come in the order of their first edge, each from its minimum edge key (y * (W + 1) + x) * 4 + d.
A ring's part is the raster rank (by first pixel) of the 8-connected foreground component it bounds: its own for an
outer ring, the surrounding one for a hole. One engine call also gives each mask's bit-quad counts and Betti numbers
from the same run labeller (so certify takes b0, b1 from an independent one).
The numpy references are tests/oracles/polygonize_oracle.py and bitquad_oracle.py; the engine equals them on all
37,093 LIDC crops and the fuzz families (tests/test_rings.py; docs/implementation2.0/kernels/K1/results).
"""
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from objdel import _gpucheck

from .topology import as_mask


@dataclass(frozen=True)
class Ring:
    pts: np.ndarray  # (k, 2) int64 corner vertices (x, y), canonical start first, offset applied
    area: int  # signed: > 0 outer, < 0 hole
    part: int  # raster rank of the bounded 8-component within its mask

    @property
    def hole(self):
        return self.area < 0


@dataclass(frozen=True)
class MaskRings:
    rings: list[Ring]
    quads: tuple[int, int, int, int, int]  # bit-quad counts n1, n2, nD, n3, n4
    b0: int  # 8-connected foreground components
    b1: int  # bounded 4-connected background components


def trace(masks: Sequence, offsets: Sequence[tuple[int, int]] | None = None) -> list[MaskRings]:
    """Rings, bit-quads and Betti numbers of 2-D masks (any nonzero value is foreground, sides 0..4096; larger raise
    ValueError) in one engine call. offsets[i] (x, y) is added to mask i's corners, e.g. after cropping."""
    ms = [as_mask(m) for m in masks]
    raw = np.zeros((len(ms), 2), np.int64) if offsets is None else np.asarray(offsets)
    if raw.size == 0 and not ms:
        raw = np.zeros((0, 2), np.int64)
    if raw.dtype.kind not in "iuf":
        raise ValueError(f"offsets dtype {raw.dtype}, need real numbers")
    if raw.shape != (len(ms), 2):
        raise ValueError(f"offsets shape {raw.shape}, need ({len(ms)}, 2)")
    if raw.size and not (np.issubdtype(raw.dtype, np.integer) or np.array_equal(raw, np.round(raw))):
        raise ValueError("offsets must be whole numbers")
    if raw.size and np.abs(raw).max() > 2**62:  # corners (<= 4097) plus offset must fit in int64
        raise ValueError("offsets out of range")
    offs = raw.astype(np.int64)
    start = np.zeros(len(ms) + 1, np.int64)
    start[1:] = np.cumsum([m.size for m in ms])
    bits = np.concatenate([m.view(np.uint8).ravel() for m in ms]) if ms else np.zeros(0, np.uint8)
    out = _gpucheck.check_masks(bits, np.array([m.shape for m in ms], np.int32).reshape(-1, 2), start)
    rs, ps, pts = out["ring_start"], out["pts_start"], out["pts"].astype(np.int64)
    area, part = out["ring_area"].tolist(), out["ring_part"].tolist()
    quads, betti = out["quads"].tolist(), out["betti"].tolist()
    return [MaskRings([Ring(pts[ps[r]:ps[r + 1]] + offs[i], area[r], part[r]) for r in range(rs[i], rs[i + 1])],
                      tuple(quads[i]), *betti[i]) for i in range(len(ms))]


def polygonize(mask, offset=(0, 0)) -> list[Ring]:
    """Rings of one 2-D mask; offset (x, y) is added to the corners."""
    return trace([mask], [offset])[0].rings
