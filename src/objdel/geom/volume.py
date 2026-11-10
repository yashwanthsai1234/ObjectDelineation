"""3-D topology of voxel volumes: slice masks stacked along z, and the volume's Betti numbers.

A CT nodule is delineated slice by slice; its outlines, rasterized (`objdel.geom.shape.to_masks`) and stacked,
form a voxel mask. Each voxel is a closed unit cube, so parts are 26-connected and the background is 6-connected
(the dual pair that keeps the counts consistent, as 8/4 per slice):
- b0, parts: 26-connected foreground components;
- b2, cavities: 6-connected background components that do not reach the border;
- b1, tunnels: from the Euler characteristic chi = V - E + F - C of the cube complex, chi = b0 - b1 + b2.
`topology_any` is used by tests; phase G (3D) uses it for topology against any doctor.
"""

from collections.abc import Iterable, Mapping

import numpy as np
from scipy import ndimage

PARTS = np.ones((3, 3, 3), bool)  # 26-connectivity
BACKGROUND = ndimage.generate_binary_structure(3, 1)  # 6-connectivity


def stack(slices: Mapping[int, np.ndarray]) -> tuple[int, np.ndarray]:
    """(z0, volume): 2-D masks keyed by slice index placed at volume[z - z0]; missing slices in between are empty."""
    shapes = {np.shape(m) for m in slices.values()}
    if len(shapes) != 1:
        raise ValueError(f"slices must have the same shape, got {sorted(shapes)}")
    z0, z1 = min(slices), max(slices)
    vol = np.zeros((z1 - z0 + 1, *shapes.pop()), bool)
    for z, m in slices.items():
        vol[z - z0] = np.asarray(m, bool)
    return z0, vol


def euler(vol: np.ndarray) -> int:
    """Euler characteristic of the union of the closed unit cubes at the True voxels: V - E + F - C, each cell counted
    once when any cube containing it is present."""
    v = np.pad(np.asarray(vol, bool), 1)

    def cells(axes) -> int:  # cells spanning `axes`: present when any of the 2^(3 - len(axes)) cubes around is
        c = v
        for ax in {0, 1, 2} - set(axes):
            c = c[(slice(None),) * ax + (slice(1, None),)] | c[(slice(None),) * ax + (slice(None, -1),)]
        return int(c.sum())

    vertices = cells(())
    edges = sum(cells((a,)) for a in range(3))
    faces = sum(cells(f) for f in ((0, 1), (0, 2), (1, 2)))
    return vertices - edges + faces - int(v.sum())


def betti3d(vol) -> tuple[int, int, int]:
    """(b0, b1, b2) of a 3-D binary volume: parts, tunnels, cavities."""
    v = np.asarray(vol, bool)
    if not v.any():
        return 0, 0, 0
    _, b0 = ndimage.label(v, structure=PARTS)
    _, background = ndimage.label(~np.pad(v, 1), structure=BACKGROUND)
    b2 = background - 1  # the outer background is one component
    return b0, b0 + b2 - euler(v), b2


def topology_any(pred, readers: Iterable) -> bool:
    """True when some reader's volume has the predicted volume's (b0, b1, b2)."""
    return betti3d(pred) in {betti3d(r) for r in readers}
