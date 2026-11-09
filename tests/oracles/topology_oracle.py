"""Slow, independent oracles for objdel.checker.topology (tests only).

Convention: 8-connected foreground, 4-connected background. A hole is a bounded 4-connected background
component; the holes of a component are the holes it has when it is the only foreground.

- scipy oracles label each component on its own (no hole-to-component assignment trick);
- `instance_topology_bfs` is pure Python (no labeling library at all), for small masks.
"""
from collections import deque

import numpy as np
from scipy import ndimage

FG8 = np.ones((3, 3), dtype=bool)
BG4 = ndimage.generate_binary_structure(2, 1)


def instance_topology(mask) -> list[tuple[int, int]]:
    """Sorted (area, holes) per 8-connected component: each component alone, padded, its bounded 4-regions."""
    m = np.asarray(mask, dtype=bool)
    lab, _ = ndimage.label(m, structure=FG8)
    out = []
    for k, sl in enumerate(ndimage.find_objects(lab), start=1):
        comp = lab[sl] == k  # holes cannot leave the component's bounding box
        _, nbg = ndimage.label(~np.pad(comp, 1), structure=BG4)
        out.append((int(comp.sum()), nbg - 1))
    return sorted(out, reverse=True)


def hole_areas(mask) -> list[int]:
    """Pixel counts of the bounded 4-connected background components, descending."""
    bg, n = ndimage.label(~np.pad(np.asarray(mask, dtype=bool), 1), structure=BG4)
    sizes = np.bincount(bg.ravel(), minlength=n + 1)
    return sorted((int(s) for k, s in enumerate(sizes) if k not in (0, bg[0, 0])), reverse=True)


def size_filter(mask, min_area: int) -> np.ndarray:
    """Drop 8-connected components with < min_area pixels, then fill holes with < min_area pixels."""
    m = np.asarray(mask, dtype=bool).copy()
    lab, _ = ndimage.label(m, structure=FG8)
    sizes = np.bincount(lab.ravel())
    small = sizes < min_area
    small[0] = False
    m[small[lab]] = False
    bg, _ = ndimage.label(~np.pad(m, 1), structure=BG4)
    sizes = np.bincount(bg.ravel())
    small = sizes < min_area
    small[[0, bg[0, 0]]] = False
    m[small[bg[1:-1, 1:-1]]] = True
    return m


def signature(mask, min_area: int = 0) -> tuple[int, ...]:
    return tuple(sorted((h for _, h in instance_topology(size_filter(mask, min_area))), reverse=True))


N8 = [(dy, dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dy or dx]
N4 = [(-1, 0), (1, 0), (0, -1), (0, 1)]


def _flood(grid: np.ndarray, nbrs) -> tuple[np.ndarray, int]:
    h, w = grid.shape
    lab = np.zeros((h, w), np.int64)
    n = 0
    for r in range(h):
        for c in range(w):
            if grid[r, c] and not lab[r, c]:
                n += 1
                lab[r, c] = n
                q = deque([(r, c)])
                while q:
                    y, x = q.popleft()
                    for dy, dx in nbrs:
                        yy, xx = y + dy, x + dx
                        if 0 <= yy < h and 0 <= xx < w and grid[yy, xx] and not lab[yy, xx]:
                            lab[yy, xx] = n
                            q.append((yy, xx))
    return lab, n


def instance_topology_bfs(mask) -> list[tuple[int, int]]:
    """Pure-Python BFS version of `instance_topology`, for small masks."""
    m = np.asarray(mask, dtype=bool)
    lab, n = _flood(m, N8)
    out = []
    for k in range(1, n + 1):
        comp = lab == k
        _, nbg = _flood(~np.pad(comp, 1), N4)
        out.append((int(comp.sum()), nbg - 1))
    return sorted(out, reverse=True)


# ---- the OpenCV versions of objdel.checker.topology before D4b (same values and types as the engine versions) ----
import cv2  # noqa: E402

from objdel.checker.topology import InstanceTopology, as_mask, bbox  # noqa: E402

_EMPTY = np.zeros(0, np.int64)
_EMPTY.flags.writeable = False


def _label(img: np.ndarray, connectivity: int) -> tuple[int, np.ndarray]:
    # opencv-python 5.0.0 rejects bool input and the `ccltype=` keyword, hence uint8 and positional arguments.
    return cv2.connectedComponentsWithAlgorithm(img, connectivity, cv2.CV_32S, cv2.CCL_SPAGHETTI)


def _crop_pad(m: np.ndarray):
    x, y, w, h = box = bbox(m)
    if w == 0:
        return None, box
    return cv2.copyMakeBorder(m[y:y + h, x:x + w].view(np.uint8), 1, 1, 1, 1, cv2.BORDER_CONSTANT, value=0), box


def cv2_instance_topology(mask) -> InstanceTopology:
    mp, _ = _crop_pad(as_mask(mask))
    if mp is None:
        return InstanceTopology(_EMPTY, _EMPTY, _EMPTY)
    n1, lab = _label(mp, 8)
    areas = np.array([np.count_nonzero(mp)]) if n1 == 2 else np.bincount(lab.ravel(), minlength=n1)[1:]
    holes = np.zeros(n1, np.int64)
    hole_areas = _EMPTY
    nb, bl = _label(1 - mp, 4)
    if nb > 2:
        w2 = mp.shape[1]
        outside = bl[0, 0]
        idx = np.flatnonzero((bl[1:] != outside) & (mp[1:] == 0) & (mp[:-1] == 1)) + w2
        hole_ids, first = np.unique(bl.ravel()[idx], return_index=True)
        holes += np.bincount(lab.ravel()[idx[first] - w2], minlength=n1)
        hole_areas = np.sort(np.bincount(bl.ravel(), minlength=nb)[hole_ids])[::-1]
    areas, holes = areas.astype(np.int64), holes[1:]
    order = np.lexsort((-holes, -areas))
    return InstanceTopology(areas[order], holes[order], hole_areas.astype(np.int64))


def cv2_betti(mask) -> tuple[int, int]:
    t = cv2_instance_topology(mask)
    return int(t.areas.size), int(t.holes.sum())


def cv2_size_filter(mask, min_area: int) -> np.ndarray:
    out = as_mask(mask)
    mp, (x, y, w, h) = _crop_pad(out)
    if mp is None:
        return out
    n1, lab = _label(mp, 8)
    small = np.bincount(lab.ravel(), minlength=n1) < min_area
    small[0] = False
    mp[small[lab]] = 0
    nb, bl = _label(1 - mp, 4)
    small = np.bincount(bl.ravel(), minlength=nb) < min_area
    small[[0, bl[0, 0]]] = False
    mp[small[bl]] = 1
    out[y:y + h, x:x + w] = mp[1:-1, 1:-1]
    return out


def cv2_signature(mask, min_area: int = 0) -> tuple[int, ...]:
    t = cv2_instance_topology(mask)
    if (t.areas < min_area).any() or (t.hole_areas < min_area).any():
        t = cv2_instance_topology(cv2_size_filter(mask, min_area))
    return tuple(sorted(t.holes.tolist(), reverse=True))


def cv2_components(mask, min_area: int = 0):
    """The engine's raw topology fields of the size-filtered mask: per 8-component in raster order of its first pixel
    (area, holes, [hole areas in raster order of their first pixels]), from OpenCV labels ranked by first pixel."""
    mp, _ = _crop_pad(cv2_size_filter(mask, min_area))
    if mp is None:
        return []
    n1, lab = _label(mp, 8)
    nb, bl = _label(1 - mp, 4)
    flat, bflat, w2 = lab.ravel(), np.where(mp.ravel() == 0, bl.ravel(), -1), mp.shape[1]
    ids, first = np.unique(flat, return_index=True)
    comps = ids[ids > 0][np.argsort(first[ids > 0])]
    out = {int(k): [int(n), 0, []] for k, n in zip(comps, np.bincount(flat, minlength=n1)[comps], strict=True)}
    hids, hfirst = np.unique(bflat, return_index=True)
    keep = (hids >= 0) & (hids != bl[0, 0])
    sizes = np.bincount(bflat[bflat >= 0], minlength=nb)
    for h, f in sorted(zip(hids[keep].tolist(), hfirst[keep].tolist(), strict=True), key=lambda t: t[1]):
        c = out[int(flat[f - w2])]
        c[1] += 1
        c[2].append(int(sizes[h]))
    return [(a, b, hs) for a, b, hs in out.values()]
