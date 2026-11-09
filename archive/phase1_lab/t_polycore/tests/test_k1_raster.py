"""K1: rasterize == oracle (both fill rules), bit for bit, incl. on-edge centres."""
import random

import numpy as np
import pytest

import oracle as O
import polycore as pc
from conftest import env_int
from helpers import half_pixel_walk, k1_case, on_boundary_centres, winding_np
from k1_runner import run

N_K1 = env_int("POLYCORE_K1_N", 10000)
N_PURE = env_int("POLYCORE_K1_PURE_N", 400)


def test_vectorised_oracle_is_the_oracle():
    """winding_np is a literal vectorisation of oracle.winding_and_parity; prove it on 300 cases."""
    for i in range(300):
        _, rings, W, H = k1_case(i, seed=99)
        W, H = min(W, 64), min(H, 37)
        w0, c0 = O.winding_and_parity(rings, W, H)
        w1, c1 = winding_np(rings, W, H)
        assert (w0 == w1).all() and (c0 == c1).all(), i


def test_k1_10k_shapes_vectorised_oracle():
    st, fails = run(N_K1)
    print(st)
    assert st["cases"] == N_K1
    assert fails == [], fails[:10]
    assert st["bad_winding_px"] == st["bad_evenodd_px"] == st["bad_nonzero_px"] == st["bad_ring_flags"] == 0
    assert st["sampled_on_edge_centres"] > 0 and st["eo_ne_nz_cases"] > 0 and st["mixed_sign_rings"] > 0


def test_k1_pure_oracle_subset():
    st, fails = run(N_PURE, pure=True, seed=4242)
    assert fails == [], fails[:10]


def test_top_left_rule_square_on_centres():
    # W=H=10, G=1000: centres at 50, 150, ...; square edges pass exactly through centres.
    sq = [(50, 50), (150, 50), (150, 150), (50, 150)]
    r = pc.rasterize([sq], 10, 10)
    exp = np.zeros((10, 10), bool)
    exp[0, 0] = True  # left + top edges in, right + bottom out
    assert (r["evenodd"] == exp).all() and (r["nonzero"] == exp).all()
    assert (r["evenodd"] == O.rasterize([sq], 10, 10)).all()


def test_shared_edge_partition():
    """A rectangle cut by a polyline into A and B (shared path traversed in opposite
    directions): w(A) + w(B) == w(rectangle) at every pixel, i.e. centres on the shared
    path go to exactly one side (no gaps, no double counting)."""
    rng = random.Random(3)
    for _ in range(300):
        W, H = rng.choice([7, 37, 64, 125, 250]), rng.choice([7, 37, 64, 125])
        x0, x1 = sorted(rng.sample(range(0, 1001, 8), 2))
        y0, y1 = sorted(rng.sample(range(0, 1001, 8), 2))
        mid = [(rng.randrange(x0, x1 + 1, 4) if rng.random() < .5 else rng.randint(x0, x1),
                rng.randrange(y0, y1 + 1, 4) if rng.random() < .5 else rng.randint(y0, y1))
               for _ in range(rng.randint(0, 5))]
        path = [(x0, y0)] + mid + [(x1, y1)]
        a = path + [(x1, y0)]
        b = (path + [(x0, y1)])[::-1]
        rect = [(x0, y0), (x0, y1), (x1, y1), (x1, y0)]
        wa = pc.rasterize([a], W, H, winding=True)["winding"]
        wb = pc.rasterize([b], W, H, winding=True)["winding"]
        wr = pc.rasterize([rect], W, H, winding=True)["winding"]
        assert (wa + wb == wr).all()


def test_mixed_sign_flag_bowtie_vs_xor():
    """Bow-tie: eo == nz everywhere (XOR flag misses it) but the ring has +1 and -1 lobes."""
    bow = [(100, 100), (900, 900), (900, 100), (100, 900)]
    r = pc.rasterize([bow], 64, 64)
    assert (r["evenodd"] == r["nonzero"]).all()
    assert r["ring_wmin"][0] == -1 and r["ring_wmax"][0] == 1 and r["rings_mixed_sign"] == 1
    loop2 = [(100, 100), (900, 100), (900, 900), (100, 900)] * 2  # doubly wound
    r = pc.rasterize([loop2], 64, 64)
    assert r["px_multi"] > 0 and r["area_evenodd"] == 0 and r["ring_wmax"][0] == 2


@pytest.mark.parametrize("W,H", [(8, 8), (40, 25), (125, 125), (250, 100), (500, 20)])
def test_on_edge_centres_are_exercised(W, H):
    """Half-pixel lattice walks for this crop (exact when W, H divide 500): many centres lie
    exactly on edges / vertices; the result must still equal the oracle."""
    rng = random.Random(W * H)
    seen = 0
    for _ in range(60):
        rings = [half_pixel_walk(rng, W, H, rng.randint(3, 14)) for _ in range(rng.randint(1, 2))]
        seen += int(on_boundary_centres(rings, W, H).sum())
        r = pc.rasterize(rings, W, H, winding=True)
        w, c = winding_np(rings, W, H)
        assert (r["winding"] == w).all() and (r["evenodd"] == (c % 2 == 1)).all()
    assert seen >= 30  # 37..267 exact on-edge centres per crop size
