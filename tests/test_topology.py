"""objdel.checker.topology: exact per-instance topology and the signature size rule (plan D3, §4.2, gate R6)."""
import cv2
import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays
from oracles import topology_oracle as oracle
from scipy import ndimage

from objdel.checker import as_mask, betti, instance_topology, min_area, signature, size_filter
from objdel.checker.topology import bbox
from tests.oracles.bitquad_oracle import bitquads

N_RANDOM = 10_000


def pairs(t):
    return list(zip(t.areas.tolist(), t.holes.tolist(), strict=True))


def check(m):
    """Instance topology equals the oracle, in canonical order; the bit-quad Euler number certifies it."""
    t = instance_topology(m)
    ref = oracle.instance_topology(m)
    assert pairs(t) == ref
    assert t.hole_areas.tolist() == oracle.hole_areas(m)
    assert t.holes.sum() == t.hole_areas.size and t.areas.sum() == np.count_nonzero(m)
    assert bitquads(m).chi8 == sum(1 - h for _, h in ref)  # chi = b0 - b1, an independent method
    assert betti(m) == (len(ref), sum(h for _, h in ref))
    return t


def grid(rows):
    return np.array([[c == "#" for c in row] for row in rows], dtype=bool)


def disk(n, r, cy, cx):
    y, x = np.mgrid[:n, :n]
    return np.hypot(y - cy, x - cx) <= r


# ---------------------------------------------------------------- generators (random and adversarial)

def smooth(rng, h, w, sigma):
    f = ndimage.gaussian_filter(rng.standard_normal((h, w)), sigma, mode="wrap")
    return f / (f.std() or 1.0)


def random_mask(rng, h, w):
    """One of ten mask families: densities 5-95 %, 1-px features, diagonal touches, holes at the crop edge,
    islands in holes in islands, many tiny parts inside one hole."""
    y, x = np.mgrid[:h, :w]
    kind = int(rng.integers(0, 10))
    if kind == 0:  # i.i.d. pixels: diagonal contacts and 1-px holes everywhere
        return rng.random((h, w)) < rng.choice([0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 0.95])
    if kind == 1:  # thresholded smooth noise: many parts, holes, islands in holes
        return smooth(rng, h, w, rng.uniform(0.7, 4.0)) > rng.normal(0, 0.5)
    if kind == 2:  # complement of the above: holes and frames touching the border
        return ~(smooth(rng, h, w, rng.uniform(0.7, 4.0)) > rng.normal(0, 0.5))
    if kind == 3:  # one smooth object with punched holes, some of them 1 px
        r = np.hypot((y - h / 2) / (h / 2), (x - w / 2) / (w / 2))
        m = 0.6 - r + 0.25 * smooth(rng, h, w, max(1.0, min(h, w) / 8)) > 0
        for _ in range(rng.integers(0, 4)):
            cy, cx, rr = rng.integers(0, h), rng.integers(0, w), rng.uniform(0.0, max(1.0, min(h, w) / 8))
            m &= np.hypot(y - cy, x - cx) > rr
        return m
    if kind == 4:  # checkerboard: every pixel touches its neighbours only diagonally
        return (y + x + rng.integers(0, 2)) % 2 == 0
    if kind == 5:  # concentric square rings: island in hole in island ...
        d = np.minimum(np.minimum(y, x), np.minimum(h - 1 - y, w - 1 - x))
        return (d // rng.integers(1, 4)) % 2 == 0
    if kind == 6:  # diagonal-only (8-connected) rings around 4-connected holes
        d = np.abs(y - h // 2) + np.abs(x - w // 2)
        return (d % 3 == 0) & (d > 0)
    if kind == 7:  # sparse dots on a frame: many tiny parts inside one hole
        m = rng.random((h, w)) < 0.1
        m[[0, -1], :] = m[:, [0, -1]] = True
        return m
    if kind == 8:  # noise with random rows/columns of the crop edge set
        m = smooth(rng, h, w, rng.uniform(0.7, 3.0)) > 0
        m[0, :] = rng.random(w) < 0.5
        m[:, -1] = rng.random(h) < 0.5
        return m
    m = rng.random((h, w)) < 0.5  # 1-px structures: lines, single pixels and diagonals
    return m & ((y % 2 == 0) | (x % 3 == 0))


def random_masks(n, seed):
    rng = np.random.default_rng(seed)
    for i in range(n):
        hi = 200 if i % 10 == 0 else 48
        h, w = (int(v) for v in rng.integers(1, hi + 1, size=2))
        yield random_mask(rng, h, w)


# ---------------------------------------------------------------- exactness

@pytest.mark.parametrize("rows, expected", [
    (["#"], [(1, 0)]),
    (["#.", ".#"], [(2, 0)]),  # diagonal pixels: one 8-connected component
    ([".#.", "#.#", ".#."], [(4, 1)]),  # diagonal ring: its 4-connected hole survives
    (["###", "#.#", "###"], [(8, 1)]),
    (["#####", "#...#", "#.#.#", "#...#", "#####"], [(16, 1), (1, 0)]),  # island in a hole
    (["##.##", "#...#", "##.##"], [(5, 0), (5, 0)]),  # gap to the border: no hole
    (["#.#", "...", "#.#"], [(1, 0)] * 4),
    (["#.#.", ".#.#", "#.#.", ".#.#"], [(8, 2)]),  # checkerboard: two interior 4-isolated pixels
])
def test_known_shapes(rows, expected):
    assert pairs(check(grid(rows))) == expected


def test_island_in_a_hole_does_not_count_for_the_enclosing_component():
    m = np.zeros((9, 9), bool)
    m[1:8, 1:8] = True
    m[2:7, 2:7] = False  # one hole of the outer ring ...
    m[3, 3:6] = m[5, 3:6] = m[3:6, 3] = m[3:6, 5] = True  # ... holding an island ring with its own hole
    t = check(m)
    assert pairs(t) == [(24, 1), (8, 1)] and t.hole_areas.tolist() == [16, 1]


def test_random_masks_vs_oracle():
    shapes = set()
    for m in random_masks(N_RANDOM, seed=2026):
        check(m)
        shapes.add(m.shape)
    assert len(shapes) > 1000


def test_small_masks_vs_pure_python_oracle():
    rng = np.random.default_rng(5)
    for _ in range(300):
        m = random_mask(rng, *(int(v) for v in rng.integers(1, 13, size=2)))
        assert pairs(instance_topology(m)) == oracle.instance_topology_bfs(m)


@pytest.mark.data
@pytest.mark.parametrize("name", ["riga_experts", "lidc_experts", pytest.param("lidc_full", marks=pytest.mark.slow)])
def test_fixture_masks_vs_oracle(fixtures_dir, name):
    n = 0
    with np.load(fixtures_dir / f"{name}.npz") as z:
        for key in z.files:
            stack = z[key].astype(bool)
            union = stack.any(0)
            a = min_area(int(union.sum()))
            y, x = np.nonzero(union)
            for m in stack[:, max(y.min() - 2, 0):y.max() + 3, max(x.min() - 2, 0):x.max() + 3]:  # oracle speed
                check(m)
                check(~m)  # complement: a frame with holes, touching the border
                assert signature(m, a) == oracle.signature(m, a)
                n += 1
    assert n > 1000


# ---------------------------------------------------------------- edge cases and determinism

@pytest.mark.parametrize("shape", [(0, 0), (0, 4), (1, 1), (1, 7), (7, 1), (5, 6)])
def test_empty_mask(shape):
    m = np.zeros(shape, bool)
    t = instance_topology(m)
    assert t.areas.size == t.holes.size == t.hole_areas.size == 0
    assert signature(m) == signature(m, 50) == () and betti(m) == (0, 0)
    assert np.array_equal(size_filter(m, 5), m)


@pytest.mark.parametrize("shape", [(1, 1), (1, 9), (9, 1), (6, 7)])
def test_full_mask(shape):
    m = np.ones(shape, bool)
    assert pairs(check(m)) == [(m.size, 0)]
    assert signature(m) == (0,)


def test_bbox_is_the_foreground_bounding_box():
    """bbox (shared by the topology functions and certify) equals the numpy box; (0, 0, 0, 0) when empty."""
    masks = [np.zeros(s, bool) for s in [(0, 0), (0, 4), (3, 0), (1, 1), (5, 6)]] + [np.ones((6, 7), bool)]
    for m in [*masks, *random_masks(300, seed=17)]:
        ys, xs = np.nonzero(m)
        ref = (0, 0, 0, 0) if ys.size == 0 else (xs.min(), ys.min(), xs.max() - xs.min() + 1, ys.max() - ys.min() + 1)
        assert bbox(as_mask(m)) == tuple(int(v) for v in ref)


def test_mask_touching_every_border():
    m = np.ones((6, 8), bool)
    m[2:4, 3:5] = False  # one hole
    m[0, 0] = m[-1, -1] = False  # dents in the corners are not holes
    assert pairs(check(m)) == [(42, 1)] and check(m).hole_areas.tolist() == [4]
    assert pairs(check(~m)) == [(4, 0), (1, 0), (1, 0)]


def test_one_row_and_one_column():
    row = np.array([[1, 1, 0, 1, 0, 0, 1, 1, 1]], bool)
    assert pairs(check(row)) == [(3, 0), (2, 0), (1, 0)]
    assert pairs(check(row.T)) == [(3, 0), (2, 0), (1, 0)]


def test_deterministic_across_calls_and_opencv_threads():
    rng = np.random.default_rng(9)
    masks = [random_mask(rng, 512, 384) for _ in range(30)]  # tall enough for OpenCV's parallel labeling
    before = cv2.getNumThreads()
    try:
        out = {}
        for threads in (0, 4, 0):  # 0 = sequential on every backend (1 is a no-op on macOS's GCD backend, I14)
            cv2.setNumThreads(threads)
            out.setdefault(threads, []).append([
                (t.areas.tobytes(), t.holes.tobytes(), t.hole_areas.tobytes(), signature(m, 7))
                for m in masks for t in [instance_topology(m)]])
        assert out[0][0] == out[0][1] == out[4][0]
    finally:
        cv2.setNumThreads(before)


def test_opencv_set_num_threads_zero_is_sequential():
    """I14: setNumThreads(0) disables OpenCV's threading on every parallel backend; setNumThreads(1) does not on the
    GCD backend (macOS wheels), where getNumThreads() stays at the core count."""
    before = cv2.getNumThreads()
    try:
        cv2.setNumThreads(0)
        assert cv2.getNumThreads() == 1
    finally:
        cv2.setNumThreads(before)


def _hole_square():
    m = np.zeros((12, 12), bool)
    m[2:10, 2:10] = True
    m[5, 5] = False  # a 1-px hole
    m[0, 11] = True  # a 1-px speck
    return m


@pytest.mark.parametrize("byte", [2, 128, 255])
def test_bool_masks_whose_true_bytes_are_not_one(byte):
    """I1 / B1-5: numpy keeps the raw bytes of a bool array (np.asarray(PIL "1" image, bool) stores 255); any
    nonzero byte is foreground, so every function gives the canonical mask's answer."""
    m = _hole_square()
    nb = (m.view(np.uint8) * np.uint8(byte)).view(bool)
    assert nb.view(np.uint8).max() == byte and np.array_equal(nb, m)
    assert pairs(instance_topology(nb)) == pairs(instance_topology(m)) == [(64 - 1, 1), (1, 0)]
    assert instance_topology(nb).hole_areas.tolist() == [1]
    assert signature(nb) == signature(m) == (1, 0) and signature(nb, 2) == signature(m, 2) == (0,)
    assert betti(nb) == betti(m) == (2, 1)
    for a in (0, 2):
        f = size_filter(nb, a)
        assert f.dtype == bool and f.view(np.uint8).max() <= 1 and np.array_equal(f, size_filter(m, a))
    c = as_mask(nb)
    assert c.dtype == bool and c.flags.c_contiguous and c.view(np.uint8).max() == 1 and np.array_equal(c, m)


def test_mask_values_and_shapes():
    """Any nonzero value is foreground (NaN included), for every dtype and layout (Fortran order, strided views);
    anything but 2-D is refused."""
    m = _hole_square()
    ref = pairs(instance_topology(m))
    for v in (m * np.uint8(255), m * -3, m * 0.5, np.where(m, np.nan, 0.0), m.astype(object), m.tolist(),
              np.asfortranarray(m * np.int16(7)), (m.view(np.uint8) * 3).view(bool)[::-1][::-1],
              np.pad(m, 3)[3:-3, 3:-3], np.pad(m, ((0, 0), (0, 5)))[:, :12]):
        assert pairs(instance_topology(v)) == ref and np.array_equal(as_mask(v), m)
    for bad in (np.zeros((2, 2, 2), bool), np.zeros(4, bool), np.bool_(True), np.zeros((1, 3, 3), np.uint8)):
        for f in (instance_topology, signature, betti, as_mask, lambda x: size_filter(x, 5)):
            with pytest.raises(ValueError, match="2-D"):
                f(bad)
    assert as_mask(np.zeros((2, 3, 4), bool), ndim=3).shape == (2, 3, 4)


# ---------------------------------------------------------------- size rule (R6)

@pytest.mark.parametrize("union_area, expected", [
    (0, 5), (1, 5), (1000, 5), (1001, 6), (1200, 6), (1201, 7), (20_000, 100), (20_001, 101), (10**9, 5 * 10**6),
])
def test_min_area_rule(union_area, expected):
    assert min_area(union_area) == expected  # max(5, ceil(0.5 % of the union)), exact integer arithmetic


def stamp(m, area, cy, cx, value):
    """Set `area` pixels nearest to (cy, cx) (a compact 4-connected blob) to `value`."""
    y, x = np.mgrid[-10:11, -10:11]
    order = np.lexsort((x.ravel(), y.ravel(), np.hypot(y, x).ravel()))[:area]
    out = m.copy()
    out[cy + y.ravel()[order], cx + x.ravel()[order]] = value
    return out


@pytest.fixture
def raters():
    return np.stack([disk(60, 18, 30, 30), disk(60, 19, 30, 30), disk(60, 20, 30.5, 29.5)])


def signature_set(stack):
    a = min_area(int(stack.any(0).sum()))
    return {signature(m, a) for m in stack}, a


def test_sub_filter_specks_and_pinholes_do_not_change_the_signature_set(raters):
    sigs, a = signature_set(raters)
    assert sigs == {(0,)} and a == 7  # |U| = 1,257 px -> ceil(6.3) = 7
    noisy = raters.copy()
    noisy[0] = stamp(noisy[0], a - 1, 4, 4, True)  # detached speck below the filter
    noisy[1] = stamp(noisy[1], a - 1, 30, 30, False)  # pinhole below the filter
    noisy[2] = stamp(stamp(noisy[2], 1, 55, 4, True), 1, 20, 40, False)
    assert signature_set(noisy) == (sigs, a)
    assert {signature(m) for m in noisy} == {(0, 0), (1,), (1, 0)}  # the raw signatures do see them


def test_above_filter_components_and_holes_change_the_signature_set(raters):
    _, a = signature_set(raters)
    noisy = raters.copy()
    noisy[0] = stamp(noisy[0], a, 4, 4, True)  # detached part at the filter size
    noisy[1] = stamp(noisy[1], a, 30, 30, False)  # hole at the filter size
    assert signature_set(noisy) == ({(0,), (0, 0), (1,)}, a)


def test_a_hole_is_measured_after_specks_inside_it_are_removed():
    m = np.zeros((12, 12), bool)
    m[1:11, 1:11] = True
    m[4:7, 4:7] = False  # 9-px hole ...
    m[5, 5] = True  # ... holding a 1-px speck: 8 background pixels
    t = instance_topology(m)
    assert pairs(t) == [(91, 1), (1, 0)] and t.hole_areas.tolist() == [8]
    assert signature(m, 9) == (1,) == oracle.signature(m, 9)  # speck removed first, then the hole has 9 px
    assert signature(m, 10) == (0,)


def test_a_hole_is_measured_by_its_own_pixels_not_by_what_it_encloses():
    """A thin moat (32 px) around a large island (49 px): with min_area 40 the moat is filled and the island
    merges into the enclosing ring. Measured by enclosed area (81 px), the moat would stay a hole."""
    m = np.zeros((13, 13), bool)
    m[1:12, 1:12] = True
    m[2:11, 2:11] = False
    m[3:10, 3:10] = True
    assert pairs(instance_topology(m)) == [(49, 0), (40, 1)] and instance_topology(m).hole_areas.tolist() == [32]
    assert signature(m, 40) == (0,) == oracle.signature(m, 40)
    assert signature(m, 32) == (1, 0) == oracle.signature(m, 32)


@settings(max_examples=300, deadline=None)
@given(arrays(bool, st.tuples(st.integers(1, 24), st.integers(1, 24))), st.integers(0, 12))
def test_size_filter_vs_oracle_and_idempotent(m, a):
    f = size_filter(m, a)
    assert np.array_equal(f, oracle.size_filter(m, a))
    assert np.array_equal(size_filter(f, a), f)
    t = instance_topology(f)
    assert t.areas.size == 0 or (t.areas.min() >= a and (t.hole_areas.size == 0 or t.hole_areas.min() >= a))
    assert signature(m, a) == oracle.signature(m, a)
    assert signature(m, a) == tuple(sorted(t.holes.tolist(), reverse=True))


def test_random_masks_size_filter_vs_oracle():
    rng = np.random.default_rng(17)
    for m in random_masks(2000, seed=7):
        a = int(rng.integers(1, 30))
        assert np.array_equal(size_filter(m, a), oracle.size_filter(m, a))
        assert signature(m, a) == oracle.signature(m, a)


# ---------------------------------------------------------------- symmetry (Hypothesis)

small_masks = arrays(bool, st.tuples(st.integers(1, 20), st.integers(1, 20)))


@settings(max_examples=300, deadline=None)
@given(small_masks)
def test_small_masks_vs_bfs(m):
    assert pairs(check(m)) == oracle.instance_topology_bfs(m)


@settings(max_examples=300, deadline=None)
@given(small_masks, st.integers(0, 7), st.tuples(*[st.integers(0, 5)] * 4), st.integers(0, 6))
def test_invariant_under_dihedral_maps_translation_and_padding(m, k, pad, a):
    w = np.rot90(m, k % 4)
    w = w.T if k >= 4 else w  # the 8 symmetries of the square grid preserve 8- and 4-adjacency
    v = np.pad(w, ((pad[0], pad[1]), (pad[2], pad[3])))
    ref, t = instance_topology(m), instance_topology(v)
    for f in ("areas", "holes", "hole_areas"):
        assert np.array_equal(getattr(t, f), getattr(ref, f))
    assert signature(v, a) == signature(m, a)
    inner = size_filter(v, a)[pad[0]:v.shape[0] - pad[1], pad[2]:v.shape[1] - pad[3]]
    assert np.array_equal(inner, size_filter(w, a)) and not size_filter(v, a).sum() - inner.sum()
