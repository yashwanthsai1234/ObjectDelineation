"""Hypothesis properties of objdel.checker.polycore on adversarial lattices: oracle agreement, translation,
mirror (top-left ties), ring reversal, start vertex, validity invariances, determinism across threads."""
import math
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from oracles import raster_oracle as RO
from oracles import validity_oracle as VO

from objdel.checker import polycore as pc

G = 1000
D = {n: i for i, n in enumerate(pc.DEFECTS)}
SETTINGS = settings(max_examples=200, deadline=None, suppress_health_check=[HealthCheck.too_slow])

# Coordinates biased towards lattices that put pixel centres exactly on edges and vertices.
coord = st.one_of(st.integers(0, G), st.integers(0, 40).map(lambda k: k * 25),
                  st.integers(0, 16).map(lambda k: k * 62 + 31), st.integers(0, 8).map(lambda k: k * 125))
ring = st.lists(st.tuples(coord, coord), min_size=0, max_size=12)
shape = st.lists(ring, min_size=1, max_size=3)
dims = st.sampled_from([7, 8, 10, 16, 25, 37, 40])


def masks(rings, W, H):
    r = pc.rasterize(rings, W, H)
    return r["evenodd"], r["nonzero"], r


@SETTINGS
@given(shape, dims, dims)
def test_matches_pure_oracle(rings, W, H):
    eo, nz, r = masks(rings, W, H)
    w, c = RO.winding_and_parity(rings, W, H)
    assert np.array_equal(eo, c % 2 == 1) and np.array_equal(nz, w != 0)
    ranges = list(zip(r["ring_wmin"].tolist(), r["ring_wmax"].tolist(), strict=True))
    assert ranges == RO.ring_ranges(rings, W, H, pure=True)


@SETTINGS
@given(shape, st.sampled_from([8, 10, 20, 25, 40, 50, 125, 250]), st.sampled_from([8, 10, 25, 40, 125]),
       st.integers(-3, 3), st.integers(-3, 3))
def test_translation_by_whole_pixels(rings, W, H, kx, ky):
    gx, gy = G // math.gcd(G, W), G // math.gcd(G, H)  # smallest grid step that is a whole number of pixels
    sx, sy = kx * gx * W // G, ky * gy * H // G
    moved = [[(x + kx * gx, y + ky * gy) for x, y in r] for r in rings]
    for a, b in zip(masks(rings, W, H)[:2], masks(moved, W, H)[:2], strict=True):
        r0, r1, c0, c1 = max(0, -sy), min(H, H - sy), max(0, -sx), min(W, W - sx)
        if r0 < r1 and c0 < c1:
            assert np.array_equal(b[r0 + sy:r1 + sy, c0 + sx:c1 + sx], a[r0:r1, c0:c1])


@SETTINGS
@given(shape, dims, dims)
def test_mirror_differs_only_at_centres_on_edges(rings, W, H):
    """x -> G - x maps centres onto centres; the top-left tie rule is not mirror symmetric, so the masks may
    differ exactly at centres lying on an edge, and nowhere else."""
    on_edge = RO.on_boundary_centres(rings, W, H)
    for mirrored, flip in (([[(G - x, y) for x, y in r] for r in rings], lambda m: m[:, ::-1]),
                           ([[(x, G - y) for x, y in r] for r in rings], lambda m: m[::-1, :])):
        for a, b in zip(masks(rings, W, H)[:2], masks(mirrored, W, H)[:2], strict=True):
            assert not ((flip(b) != a) & ~on_edge).any()


@SETTINGS
@given(shape, dims, dims)
def test_reversal_keeps_masks_and_negates_ring_winding(rings, W, H):
    eo, nz, a = masks(rings, W, H)
    eo2, nz2, b = masks([r[::-1] for r in rings], W, H)
    assert np.array_equal(eo, eo2) and np.array_equal(nz, nz2)
    assert np.array_equal(b["ring_wmin"], -a["ring_wmax"]) and np.array_equal(b["ring_wmax"], -a["ring_wmin"])
    assert b["px_multi"] == a["px_multi"]


@SETTINGS
@given(shape, dims, dims, st.integers(0, 11))
def test_start_vertex_invariance(rings, W, H, k):
    rot = [r[k % len(r):] + r[:k % len(r)] if r else r for r in rings]
    for a, b in zip(masks(rings, W, H)[:2], masks(rot, W, H)[:2], strict=True):
        assert np.array_equal(a, b)
    va, vb = pc.validate(rings), pc.validate(rot)
    assert np.array_equal(va["per_ring"], vb["per_ring"]) and np.array_equal(va["pair_touch"], vb["pair_touch"])


def _strip(res, *drop):
    m = res["per_ring"].copy()
    for d in drop:
        m[:, D[d]] = 0
    return m


@SETTINGS
@given(shape, st.integers(-500, 500), st.integers(-500, 500), st.integers(1, 7))
def test_validate_affine_invariance(rings, tx, ty, s):
    """Translation, integer scaling and 90-degree rotation keep every class but out_of_range; mirroring and
    reversal keep every class but wrong_winding, which flips for every ring with nonzero area."""
    base = pc.validate(rings)
    scaled = [[(s * x + tx, s * y + ty) for x, y in r] for r in rings]
    for moved in (scaled, [[(G - y, x) for x, y in r] for r in rings]):
        assert np.array_equal(_strip(base, "out_of_range"), _strip(pc.validate(moved), "out_of_range"))
    sel = (base["area2"] != 0) & (base["per_ring"][:, D["fewer_than_3_points"]] == 0)
    for other in (pc.validate([[(G - x, y) for x, y in r] for r in rings]), pc.validate([r[::-1] for r in rings])):
        keep = ("wrong_winding", "out_of_range")
        assert np.array_equal(_strip(base, *keep), _strip(other, *keep))
        ww = D["wrong_winding"]
        assert np.array_equal(other["per_ring"][sel, ww], 1 - base["per_ring"][sel, ww])


@SETTINGS
@given(ring)
def test_labels_match_oracle_on_lattices(r):
    assert VO.labels(pc.validate([r])["per_ring"][0]) == VO.ring_problems(r)


@settings(max_examples=40, deadline=None)
@given(st.lists(st.tuples(shape, dims), min_size=1, max_size=24))
def test_score_batch_deterministic_across_threads(batch):
    shapes, sizes = [s for s, _ in batch], [d for _, d in batch]
    ref = pc.score_batch(shapes, sizes, sizes, threads=1)
    for t in (0, 3, 8):
        out = pc.score_batch(shapes, sizes, sizes, threads=t)
        assert all(np.array_equal(out[k], ref[k]) for k in ref)
    with ThreadPoolExecutor(3) as ex:
        single = list(ex.map(lambda sd: pc.validate(sd[0])["total"], batch))
    assert np.array_equal(np.stack(single), ref["defects"])
