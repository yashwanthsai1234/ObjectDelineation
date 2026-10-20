import numpy as np
import pytest

from objdel.checker import betti, bitquads, certify, polygonize

N4 = [(1, 0), (-1, 0), (0, 1), (0, -1)]
N8 = N4 + [(1, 1), (1, -1), (-1, 1), (-1, -1)]


def flood_count(g, val, nb):
    H, W = g.shape
    seen = np.zeros_like(g, dtype=bool)
    n = 0
    for i in range(H):
        for j in range(W):
            if g[i, j] == val and not seen[i, j]:
                n += 1
                st = [(i, j)]
                seen[i, j] = True
                while st:
                    y, x = st.pop()
                    for dy, dx in nb:
                        yy, xx = y + dy, x + dx
                        if 0 <= yy < H and 0 <= xx < W and g[yy, xx] == val and not seen[yy, xx]:
                            seen[yy, xx] = True
                            st.append((yy, xx))
    return n


def reference(m):
    p = np.pad(m.astype(np.int64), 1)
    b0 = flood_count(p, 1, N8)
    b1 = flood_count(p, 0, N4) - 1
    edges = int(np.abs(np.diff(p, axis=0)).sum() + np.abs(np.diff(p, axis=1)).sum())
    return b0, b1, edges


def test_random_grids_brute_force():
    rng = np.random.default_rng(0)
    for t in range(3000):
        h, w = rng.integers(1, 14, size=2)
        m = rng.random((h, w)) < rng.uniform(0.15, 0.85)
        b0, b1, edges = reference(m)
        assert betti(m) == (b0, b1)
        q = bitquads(m)
        assert q.chi8 == b0 - b1 and q.area == m.sum() and q.perimeter == edges
        r = certify(m)
        assert r.certified, (t, r.failures, m.astype(int))
        assert (r.b0, r.b1) == (b0, b1)


def grid(rows):
    return np.array([[c == "#" for c in row] for row in rows])


@pytest.mark.parametrize("rows,expect", [
    ([""], (0, 0)),
    (["#"], (1, 0)),
    (["###", "#.#", "###"], (1, 1)),
    (["#####", "#...#", "#.#.#", "#...#", "#####"], (2, 1)),  # island in a hole
    (["#.", ".#"], (1, 0)),  # diagonal pixels: one 8-connected object
    ([".#.", "#.#", ".#."], (1, 1)),  # diagonal ring: 4-connected hole survives
    (["##.##", "#...#", "##.##"], (2, 0)),  # gap to the border: no hole
    (["#.#", "...", "#.#"], (4, 0)),
    (["####", "#..#", "#..#", "####"], (1, 1)),
])
def test_known_shapes(rows, expect):
    m = grid(rows) if rows != [""] else np.zeros((3, 3), bool)
    r = certify(m)
    assert r.certified, r.failures
    assert (r.b0, r.b1) == expect


def test_ring_orientation_and_canonical_start():
    m = grid(["###", "#.#", "###"])
    rings = polygonize(m)
    outer = [r for r in rings if not r.hole]
    hole = [r for r in rings if r.hole]
    assert outer[0].area == 9 and hole[0].area == -1
    assert tuple(outer[0].pts[0]) == (0, 0) and tuple(hole[0].pts[0]) == (1, 1)
    assert len(outer[0].pts) == 4 and len(hole[0].pts) == 4


def spiral(n):
    m = np.zeros((n, n), bool)
    x = y = 0
    dx, dy = 1, 0
    lo_x, lo_y, hi_x, hi_y = 0, 2, n - 1, n - 1
    for _ in range(n * n):
        m[y, x] = True
        nx, ny = x + dx, y + dy
        if dx == 1 and nx > hi_x:
            dx, dy, hi_x = 0, 1, hi_x - 2
        elif dy == 1 and ny > hi_y:
            dx, dy, hi_y = -1, 0, hi_y - 2
        elif dx == -1 and nx < lo_x:
            dx, dy, lo_x = 0, -1, lo_x + 2
        elif dy == -1 and ny < lo_y:
            dx, dy, lo_y = 1, 0, lo_y + 2
        if lo_x > hi_x + 2 or lo_y > hi_y + 2:
            break
        x, y = x + dx, y + dy
    return m


@pytest.mark.parametrize("seed,p,n", [(1, 0.5, 300), (2, 0.3, 300), (3, 0.7, 300), (4, 0.59, 400)])
def test_large_random(seed, p, n):
    m = np.random.default_rng(seed).random((n, n)) < p
    b0, b1, edges = reference(m)
    r = certify(m)
    assert r.certified, r.failures
    assert (r.b0, r.b1, r.perimeter) == (b0, b1, edges)


def test_long_ring_spiral():
    m = spiral(201)
    r = certify(m)
    assert r.certified, r.failures
    assert (r.b0, r.b1) == reference(m)[:2]
    assert max(len(x.pts) for x in polygonize(m)) > 100


def test_non_canonical_bool_bytes():
    """I1: bitquads, polygonize and certify read a bool array whose True bytes are 2 or 255 as the canonical mask
    (their numpy operations compare or cast; topology normalises at its boundary, tested in test_topology.py)."""
    m = np.zeros((9, 9), bool)
    m[1:8, 1:8] = True
    m[4, 4] = False
    for byte in (2, 255):
        nb = (m.view(np.uint8) * np.uint8(byte)).view(bool)
        assert bitquads(nb) == bitquads(m)
        assert [r.pts.tolist() for r in polygonize(nb)] == [r.pts.tolist() for r in polygonize(m)]
        assert certify(nb) == certify(m) and certify(nb, min_size=2) == certify(m, min_size=2)
