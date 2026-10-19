"""Exact per-instance 2D topology: components (8-connected) and the holes of each one.

Convention: 8-connected foreground, 4-connected background. A hole is a bounded
4-connected background component. "Holes of component k" are the holes that
component k would have if it were the only foreground (an island inside a hole
of k does not change k's count; the island's own holes are counted for the island).

Algorithm (two OpenCV Spaghetti labelings, O(pixels), no per-component loop):
0. Crop to the foreground bounding box (holes cannot leave it) and pad by one pixel.
1. Label the padded crop with 8-connectivity -> components.
2. Label the padded background with 4-connectivity; the label of pixel (0, 0)
   is the outside, every other background label is a hole.
3. The first raster pixel p of a hole has a foreground pixel directly above it
   (a background pixel there would be 4-adjacent and earlier in raster order).
   That pixel belongs to the component enclosing the hole: islands inside the
   hole lie strictly below the hole's top row. bincount then gives holes per
   component. Only hole pixels whose upper neighbour is foreground are
   candidates, so the per-hole minimum runs over a boundary-sized set.

`instance_topology_euler` is an independent second method (per-component
bit-quad Euler number; every 2x2 window holds pixels of at most one 8-connected
component), used as a certificate in tests.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import cv2
import numpy as np


def label(img_u8: np.ndarray, connectivity: int) -> tuple[int, np.ndarray]:
    """OpenCV Spaghetti labeling (Spaghetti4C for connectivity 4). Returns (n_labels incl. 0, labels)."""
    return cv2.connectedComponentsWithAlgorithm(img_u8, connectivity, cv2.CV_32S, cv2.CCL_SPAGHETTI)


@dataclass(frozen=True)
class InstanceTopology:
    labels: np.ndarray  # int32 labels of the bounding-box crop, 0 = background, 1..n = components
    origin: tuple[int, int]  # (row, col) of labels[0, 0] in the input mask
    areas: np.ndarray  # (n,) int64, pixels per component
    holes: np.ndarray  # (n,) int64, holes per component
    hole_areas: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))  # background pixels per hole

    @property
    def n(self) -> int:
        return int(self.areas.size)

    @property
    def betti(self) -> tuple[int, int]:
        return self.n, int(self.holes.sum())

    def full_labels(self, shape: tuple[int, int]) -> np.ndarray:
        out = np.zeros(shape, np.int32)
        y, x = self.origin
        out[y:y + self.labels.shape[0], x:x + self.labels.shape[1]] = self.labels
        return out

    def as_list(self) -> list[tuple[int, int]]:
        """Canonical (area, holes) list: sorted by area then holes, descending."""
        return sorted(zip(self.areas.tolist(), self.holes.tolist()), reverse=True)

    def signature(self, min_area: int = 0) -> tuple[int, ...]:
        """Sorted per-component hole counts of components with area >= min_area.

        This is the set-valued topology key: (0,) = one simply connected part,
        (0, 0) = two parts, (1,) = one part with one hole.
        """
        keep = self.areas >= min_area
        return tuple(sorted(self.holes[keep].tolist(), reverse=True))


_EMPTY = np.zeros(0, np.int64)


def bbox(m: np.ndarray) -> tuple[int, int, int, int] | None:
    """(y0, y1, x0, x1) of the foreground, or None if empty."""
    rows = np.flatnonzero(m.any(1))
    if not rows.size:
        return None
    cols = np.flatnonzero(m.any(0))
    return int(rows[0]), int(rows[-1]) + 1, int(cols[0]), int(cols[-1]) + 1


def instance_topology_full(mask: np.ndarray) -> InstanceTopology:
    m = np.asarray(mask, dtype=bool)
    box = bbox(m)
    if box is None:
        return InstanceTopology(np.zeros((0, 0), np.int32), (0, 0), _EMPTY, _EMPTY)
    y0, y1, x0, x1 = box
    h, w = y1 - y0, x1 - x0
    mp = np.zeros((h + 2, w + 2), np.uint8)  # crop to the bounding box, pad with background
    mp[1:-1, 1:-1] = m[y0:y1, x0:x1]
    n1, lab = label(mp, 8)
    if n1 == 2:
        areas = np.array([np.count_nonzero(mp)], np.int64)
    else:
        areas = np.bincount(lab.ravel(), minlength=n1)[1:].astype(np.int64)
    holes = np.zeros(n1, np.int64)
    nb, bl = label(1 - mp, 4)  # label 0 = foreground pixels, outside = bl[0, 0]
    if nb > 2:  # more than one background component -> at least one hole
        w2 = w + 2
        outside = bl[0, 0]
        # hole pixels whose upper neighbour is foreground (each hole's first raster pixel is one of them)
        cand = (bl[1:] != outside) & (mp[1:] == 0) & (mp[:-1] == 1)
        idx = np.flatnonzero(cand) + w2  # flat index in the padded crop, ascending
        hl = bl.ravel()[idx]
        hole_ids, first = np.unique(hl, return_index=True)  # first (= raster-first) candidate of each hole
        enclosing = lab.ravel()[idx[first] - w2]
        holes += np.bincount(enclosing, minlength=n1)
        hole_areas = np.bincount(bl.ravel(), minlength=nb)[hole_ids].astype(np.int64)
        return InstanceTopology(lab[1:-1, 1:-1], (y0, x0), areas, holes[1:], hole_areas)
    return InstanceTopology(lab[1:-1, 1:-1], (y0, x0), areas, holes[1:])


def size_filter(mask: np.ndarray, min_area: int) -> np.ndarray:
    """Remove 8-connected components and fill 4-connected holes with fewer than min_area pixels."""
    m = np.asarray(mask, dtype=bool)
    if min_area <= 1:
        return m.copy()
    box = bbox(m)
    if box is None:
        return m.copy()
    y0, y1, x0, x1 = box
    mp = np.zeros((y1 - y0 + 2, x1 - x0 + 2), np.uint8)
    mp[1:-1, 1:-1] = m[y0:y1, x0:x1]
    n1, lab = label(mp, 8)
    small = np.bincount(lab.ravel(), minlength=n1) < min_area
    small[0] = False
    mp[small[lab]] = 0
    nb, bl = label(1 - mp, 4)
    small = np.bincount(bl.ravel(), minlength=nb) < min_area
    small[0] = False
    small[bl[0, 0]] = False
    mp[small[bl]] = 1
    out = m.copy()
    out[y0:y1, x0:x1] = mp[1:-1, 1:-1].astype(bool)
    return out


def filtered_signature(mask: np.ndarray, min_area: int, topo: InstanceTopology | None = None) -> tuple[int, ...]:
    """Per-instance topology signature ignoring components and holes smaller than min_area pixels.

    Exact: when anything is below the filter the mask is size-filtered and its topology recomputed
    (filling a small hole can merge an island into its enclosing component).
    """
    topo = instance_topology_full(mask) if topo is None else topo
    if not np.any(topo.areas < min_area) and not np.any(topo.hole_areas < min_area):
        return topo.signature()
    return instance_topology_full(size_filter(mask, min_area)).signature()


def instance_topology(mask: np.ndarray) -> list[tuple[int, int]]:
    """List of (component area, holes in that component), canonical order (area, holes descending)."""
    return instance_topology_full(mask).as_list()


# --- independent method: per-component bit-quad Euler number -------------------------------

# chi8 contribution x4 of each 2x2 code (TL=1, TR=2, BL=4, BR=8): Q1 +1, Q3 -1, diagonal -2.
_CHI8_X4 = np.zeros(16, np.int64)
_CHI8_X4[[1, 2, 4, 8]] = 1
_CHI8_X4[[7, 11, 13, 14]] = -1
_CHI8_X4[[6, 9]] = -2


def instance_topology_euler(mask: np.ndarray) -> InstanceTopology:
    m = np.ascontiguousarray(mask, dtype=bool)
    mp = np.pad(m, 1).astype(np.uint8)
    n1, lab = label(mp, 8)
    code = mp[:-1, :-1] | (mp[:-1, 1:] << 1) | (mp[1:, :-1] << 2) | (mp[1:, 1:] << 3)
    win = np.maximum(np.maximum(lab[:-1, :-1], lab[:-1, 1:]), np.maximum(lab[1:, :-1], lab[1:, 1:]))
    hist = np.bincount((win.astype(np.int64) * 16 + code).ravel(), minlength=16 * n1).reshape(n1, 16)
    chi4 = hist @ _CHI8_X4
    assert np.all(chi4 % 4 == 0)
    holes = 1 - chi4[1:] // 4
    areas = np.bincount(lab.ravel(), minlength=n1)[1:].astype(np.int64)
    return InstanceTopology(lab[1:-1, 1:-1], (0, 0), areas, holes.astype(np.int64))
