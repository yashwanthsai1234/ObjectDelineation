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
