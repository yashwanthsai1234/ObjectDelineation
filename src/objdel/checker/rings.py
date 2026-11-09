"""Mask -> rectilinear rings on the CUDA engine (objdel._gpucheck.check_masks: nvcc build on the box, CPU build of the
same source on the Mac). The polygonizer seam of the plan (B1.4a; D3 adds trace and MaskRings behind it).

The rings are the cracks between pixels, clockwise on screen (y down) around foreground, so outer rings have positive
shoelace area and holes negative; at a saddle vertex the walk turns left (8-connected foreground, 4-connected
background). Rings come in the order of their first edge, each from its minimum edge key (y * (W + 1) + x) * 4 + d.
The numpy reference with the same rings and its sources is tests/oracles/polygonize_oracle.py; the engine equals it
on all 37,093 LIDC crops and the fuzz families (tests/test_rings.py; docs/implementation2.0/kernels/K1/results).
"""
from dataclasses import dataclass

import numpy as np

from objdel import _gpucheck

from .topology import as_mask


@dataclass
class Ring:
    pts: np.ndarray  # (k, 2) int64 corner vertices (x, y), canonical start first
    area: int  # signed: > 0 outer, < 0 hole

    @property
    def hole(self):
        return self.area < 0


def polygonize(mask, offset=(0, 0)) -> list[Ring]:
    """Rings of a 2-D mask (any nonzero value is foreground). offset (x, y) is added to vertices, e.g. after
    cropping. Sides above 4096 raise ValueError (the engine's cap)."""
    m = as_mask(mask)
    out = _gpucheck.check_masks(m.view(np.uint8).ravel(), np.array([m.shape], np.int32),
                                np.array([0, m.size], np.int64))
    pts, ps, area = out["pts"].astype(np.int64) + np.asarray(offset, np.int64), out["pts_start"], out["ring_area"]
    return [Ring(pts[ps[r]:ps[r + 1]], int(area[r])) for r in range(len(area))]
