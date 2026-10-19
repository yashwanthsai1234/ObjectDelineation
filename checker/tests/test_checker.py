import numpy as np
import pytest

from checker import AnnotatorSet, betti, bitquads, certify, group_advantages, polygonize, size_filter

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


def test_hash_depends_on_position_not_crop():
    a = np.zeros((10, 10), bool)
    a[2:5, 2:5] = True
    b = np.zeros((20, 20), bool)
    b[2:5, 2:5] = True
    c = np.zeros((10, 10), bool)
    c[3:6, 2:5] = True
    assert certify(a).hash == certify(b).hash != certify(c).hash


def test_size_filter():
    m = np.zeros((12, 12), bool)
    m[1:10, 1:10] = True
    m[4:6, 4:6] = False  # hole of 4 px
    m[11, 11] = True  # speck of 1 px
    assert betti(m) == (2, 1)
    assert betti(size_filter(m, 5)) == (1, 0)
    assert betti(size_filter(m, 2)) == (1, 1)


def make(W, H, ops):
    g = np.zeros((H, W), bool)
    for v, x0, y0, x1, y1 in ops:
        g[max(y0, 0):y1, max(x0, 0):x1] = bool(v)
    return g


def test_reproduces_artifact_grpo_demo():
    W, H = 22, 14
    annot = [make(W, H, [[1, 3, 2, 19, 12], [0, 8, 5, 13, 9]]),
             make(W, H, [[1, 3, 2, 18, 12], [0, 8, 5, 14, 10]]),
             make(W, H, [[1, 4, 3, 19, 12], [0, 9, 5, 13, 9]])]
    samples = {
        "faithful": [[1, 4, 2, 19, 12], [0, 8, 5, 13, 9]],
        "hole_filled": [[1, 3, 2, 19, 12]],
        "extra_speck": [[1, 4, 2, 19, 12], [0, 8, 5, 13, 9], [1, 0, 0, 2, 1]],
        "over_dilated": [[1, 2, 1, 20, 13], [0, 9, 6, 12, 8]],
        "broken_ring": [[1, 4, 2, 19, 12], [0, 8, 5, 13, 9], [0, 15, 2, 16, 12]],
        "slightly_thin": [[1, 4, 3, 18, 11], [0, 8, 5, 13, 9]],
    }
    aset = AnnotatorSet(annot)
    scores = [aset.score(make(W, H, ops)) for ops in samples.values()]
    adv_iou = group_advantages(scores, 1, 0, 0)
    adv = group_advantages(scores)
    expect_iou = [1.01, 0.52, 0.89, -1.38, 0.35, -1.38]
    expect = [1.94, -0.45, 0.16, -1.36, -0.36, 0.08]
    assert np.allclose(adv_iou, expect_iou, atol=0.01), adv_iou
    assert np.allclose(adv, expect, atol=0.01), adv
    assert scores[1].errors["missing_holes"] == 1
    assert scores[2].errors["extra_objects"] == 1


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
    r = certify(m, keep_rings=True)
    assert r.certified, r.failures
    assert (r.b0, r.b1) == reference(m)[:2]
    assert max(len(x.pts) for x in r.rings) > 100


def test_polygon_io():
    from checker.polygon_io import as_points, extract_json, rasterize, ring_problems
    reply = 'Sure!\n```json\n{"disc": [[10, 10], [90, 10], [90, 90], [10, 90]], "cup": [0, 0, 5]}\n```'
    d = extract_json(reply)
    assert ring_problems(as_points(d["disc"])) == []
    assert as_points(d["cup"]) is None and ring_problems(None) == ["unparseable"]
    assert "self_crossing" in ring_problems(as_points([[0, 0], [10, 10], [10, 0], [0, 10]]))
    assert ring_problems(as_points([[0, 0], [10, 0], [10, 10], [5, -5], [0, 10]])) == ["self_crossing"]
    assert ring_problems(as_points([[0, 0], [1, 1]])) == ["fewer_than_3_points"]
    assert ring_problems(as_points([[0, 0], [10, 0], [10, 10], [0, 10], [0, 0]])) == []
    m = rasterize(as_points(d["disc"]), 100, 100)
    assert certify(m).certified and betti(m) == (1, 0) and m[50, 50] and not m[5, 5]
    assert extract_json("no json here") is None
    stair = as_points([350, 350, 350, 450, 450, 450, 450, 550, 350, 550, 350, 650])
    assert ring_problems(stair) == ["self_touching"]
    assert ring_problems(as_points([[0, 0], [1200, 0], [0, 10]]), bounds=(0, 1000)) == ["out_of_range"]
