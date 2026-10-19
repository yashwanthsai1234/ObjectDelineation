"""Property tests (hypothesis): oracle agreement on adversarial lattices, translation,
mirror (tie-rule caveat), ring reversal, start-vertex rotation, validate invariances."""
import math

from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import oracle as O
import polycore as pc
from helpers import on_boundary_centres

G = 1000
D = {n: i for i, n in enumerate(pc.DEFECTS)}
SETTINGS = settings(max_examples=300, deadline=None, suppress_health_check=[HealthCheck.too_slow])

# Coordinates biased towards lattices that put pixel centres exactly on edges/vertices.
coord = st.one_of(
    st.integers(0, G),
    st.integers(0, 40).map(lambda k: k * 25),
    st.integers(0, 16).map(lambda k: k * 62 + 31),
    st.integers(0, 8).map(lambda k: k * 125),
)
point = st.tuples(coord, coord)
ring = st.lists(point, min_size=0, max_size=12)
shape = st.lists(ring, min_size=1, max_size=3)
dims = st.sampled_from([7, 8, 10, 16, 25, 37, 40])


def W_of(rings, W, H):
    return pc.rasterize(rings, W, H, winding=True)["winding"]


@SETTINGS
@given(shape, dims, dims)
def test_matches_pure_oracle(rings, W, H):
    r = pc.rasterize(rings, W, H, winding=True)
    w, c = O.winding_and_parity(rings, W, H)
    assert (r["winding"] == w).all()
    assert (r["evenodd"] == (c % 2 == 1)).all() and (r["nonzero"] == (w != 0)).all()


@SETTINGS
@given(shape, st.sampled_from([8, 10, 20, 25, 40, 50, 125, 250]), st.sampled_from([8, 10, 25, 40, 125]),
       st.integers(-3, 3), st.integers(-3, 3))
def test_translation_by_whole_pixels(rings, W, H, kx, ky):
    gx, gy = G // math.gcd(G, W), G // math.gcd(G, H)        # grid step = whole number of pixels
    sx, sy = kx * gx * W // G, ky * gy * H // G             # pixel shift
    moved = [[(x + kx * gx, y + ky * gy) for x, y in r] for r in rings]
    a, b = W_of(rings, W, H), W_of(moved, W, H)
    # b[r + sy, c + sx] == a[r, c] wherever both indices are inside the crop
    r0, r1 = max(0, -sy), min(H, H - sy)
    c0, c1 = max(0, -sx), min(W, W - sx)
    if r0 < r1 and c0 < c1:
        assert (b[r0 + sy:r1 + sy, c0 + sx:c1 + sx] == a[r0:r1, c0:c1]).all()


@SETTINGS
@given(shape, dims, dims)
def test_mirror_only_differs_on_edge_centres(rings, W, H):
    """x -> G - x maps pixel centres onto pixel centres; the top-left tie rule is not mirror
    symmetric, so masks may differ exactly at centres lying on an edge, nowhere else."""
    for mir, flip in (([[(G - x, y) for x, y in r] for r in rings], lambda m: m[:, ::-1]),
                      ([[(x, G - y) for x, y in r] for r in rings], lambda m: m[::-1, :])):
        a = pc.rasterize(rings, W, H)
        b = pc.rasterize(mir, W, H)
        B = on_boundary_centres(rings, W, H)
        for k in ("evenodd", "nonzero"):
            diff = flip(b[k]) != a[k]
            assert not (diff & ~B).any()


@SETTINGS
@given(shape, dims, dims)
def test_reversal_negates_winding(rings, W, H):
    a = pc.rasterize(rings, W, H, winding=True)
    b = pc.rasterize([r[::-1] for r in rings], W, H, winding=True)
    assert (b["winding"] == -a["winding"]).all()
    assert (a["evenodd"] == b["evenodd"]).all() and (a["nonzero"] == b["nonzero"]).all()
    assert (b["ring_wmin"] == -a["ring_wmax"]).all()


def _rot(r, k):
    return r[k % len(r):] + r[:k % len(r)] if r else r


@SETTINGS
@given(shape, dims, dims, st.integers(0, 11))
def test_start_vertex_invariance(rings, W, H, k):
    rot = [_rot(r, k) for r in rings]
    assert (W_of(rings, W, H) == W_of(rot, W, H)).all()
    a, b = pc.validate(rings), pc.validate(rot)
    assert (a["ring"] == b["ring"]).all() and (a["pair_touch"] == b["pair_touch"]).all()


def _strip(res, *drop):
    m = res["ring"].copy()
    for d in drop:
        m[:, D[d]] = 0
    return m


@SETTINGS
@given(shape, st.integers(-500, 500), st.integers(-500, 500), st.integers(1, 7))
def test_validate_affine_invariance(rings, tx, ty, s):
    """Translation, integer scaling, 90-degree rotation keep every class (except out_of_range);
    mirroring and reversal keep every class except wrong_winding, which flips."""
    base = pc.validate(rings)
    tr = pc.validate([[(s * x + tx, s * y + ty) for x, y in r] for r in rings])
    assert (_strip(base, "out_of_range") == _strip(tr, "out_of_range")).all()
    rot = pc.validate([[(G - y, x) for x, y in r] for r in rings])
    assert (_strip(base, "out_of_range") == _strip(rot, "out_of_range")).all()
    for other in (pc.validate([[(G - x, y) for x, y in r] for r in rings]), pc.validate([r[::-1] for r in rings])):
        assert (_strip(base, "wrong_winding", "out_of_range") == _strip(other, "wrong_winding", "out_of_range")).all()
        sel = (base["area2"] != 0) & (base["ring"][:, D["fewer_than_3_points"]] == 0)
        ww = D["wrong_winding"]
        assert (other["ring"][sel, ww] == 1 - base["ring"][sel, ww]).all()


@SETTINGS
@given(ring)
def test_labels_match_oracle_on_lattices(r):
    assert pc.ring_labels(pc.validate([r])["ring"][0]) == O.ring_problems(r)
