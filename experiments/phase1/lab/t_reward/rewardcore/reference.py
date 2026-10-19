"""Slow, obviously-correct references used only by tests."""
from __future__ import annotations

from collections import deque

import numpy as np
from scipy import ndimage

FG8 = np.ones((3, 3), dtype=bool)
BG4 = ndimage.generate_binary_structure(2, 1)


def instance_topology_bruteforce(mask: np.ndarray) -> list[tuple[int, int]]:
    """Per component: the component alone, padded, its bounded 4-connected background regions."""
    m = np.asarray(mask, dtype=bool)
    lab, n = ndimage.label(m, structure=FG8)
    out = []
    for k, sl in enumerate(ndimage.find_objects(lab), start=1):
        comp = lab[sl] == k  # the component alone, cropped to its bounding box (holes cannot leave it)
        _, nbg = ndimage.label(~np.pad(comp, 1), structure=BG4)
        out.append((int(comp.sum()), nbg - 1))
    return sorted(out, reverse=True)


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


N8 = [(dy, dx) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if dy or dx]
N4 = [(-1, 0), (1, 0), (0, -1), (0, 1)]


def instance_topology_bfs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Pure-Python BFS version (no scipy), for small masks: independent of any labeling library."""
    m = np.asarray(mask, dtype=bool)
    lab, n = _flood(m, N8)
    out = []
    for k in range(1, n + 1):
        comp = lab == k
        _, nbg = _flood(~np.pad(comp, 1), N4)
        out.append((int(comp.sum()), nbg - 1))
    return sorted(out, reverse=True)
