"""Garbage in must never crash: errors are Python exceptions, otherwise results stay exact
(huge coordinates take the __int128 path and are checked against the pure oracle)."""
import random

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

import oracle as O
import polycore as pc
from polycore import _polycore as C

I32 = np.iinfo(np.int32)
OK_ERRORS = (ValueError, TypeError, OverflowError)


def _rand_xy(rng, n):
    mode = int(rng.integers(4))
    if mode == 0:
        return rng.integers(I32.min, I32.max, (n, 2), dtype=np.int64, endpoint=True).astype(np.int32)
    if mode == 1:
        return rng.choice(np.array([I32.min, I32.min + 1, -1, 0, 1, 1000, I32.max - 1, I32.max], np.int32), (n, 2))
    if mode == 2:
        return rng.integers(-3, 4, (n, 2)).astype(np.int32)
    return rng.integers(-5000, 6000, (n, 2)).astype(np.int32)


def test_random_garbage_never_crashes():
    rng = np.random.default_rng(0)
    n_ok = n_err = 0
    for it in range(4000):
        n = int(rng.integers(0, 30))
        xy = _rand_xy(rng, n)
        R = int(rng.integers(0, 6))
        if rng.random() < 0.8:
            off = np.sort(rng.integers(0, n + 1, R + 1)).astype(np.int32)
        else:  # malformed offsets
            off = rng.integers(-3, n + 4, R + 1).astype(np.int32)
        W, H = int(rng.integers(-2, 70)), int(rng.integers(-2, 70))
        G = int(rng.choice([0, 1, 7, 1000, 2**30, 2**31 - 1]))
        roles = None
        if R and rng.random() < 0.3:
            roles = rng.integers(-2, R + 1, R).astype(np.int32)
        for fn in (lambda: C.rasterize(xy, off, W, H, G, True),
                   lambda: C.validate(xy, off, roles, G, "sweep"),
                   lambda: C.score(xy, off, W, H, roles, None, G)):
            try:
                fn()
                n_ok += 1
            except OK_ERRORS:
                n_err += 1
    assert n_ok > 1000 and n_err > 100


def test_malformed_inputs_raise():
    xy = np.zeros((4, 2), np.int32)
    good = np.array([0, 4], np.int32)
    bad = [
        lambda: C.rasterize(xy, np.array([0, 5], np.int32), 8, 8, 1000, False),     # off past end
        lambda: C.rasterize(xy, np.array([2, 1], np.int32), 8, 8, 1000, False),     # decreasing
        lambda: C.rasterize(xy, np.array([-1, 2], np.int32), 8, 8, 1000, False),    # negative
        lambda: C.rasterize(xy, np.array([], np.int32), 8, 8, 1000, False),         # no offsets
        lambda: C.rasterize(np.zeros(3, np.int32), good, 8, 8, 1000, False),        # odd flat length
        lambda: C.rasterize(np.zeros((2, 3), np.int32), good, 8, 8, 1000, False),   # bad shape
        lambda: C.rasterize(xy, good, 0, 8, 1000, False),
        lambda: C.rasterize(xy, good, 8, 5000, 1000, False),
        lambda: C.rasterize(xy, good, 8, 8, 0, False),
        lambda: C.rasterize(xy, good, 8, 8, 2**31, False),
        lambda: C.rasterize(xy.astype(np.float64), good, 8, 8, 1000, False),        # wrong dtype
        lambda: C.validate(xy, good, np.array([0], np.int32), 1000, "sweep"),        # self-parent
        lambda: C.validate(xy, good, None, 1000, "nope"),
        lambda: C.validate(xy, np.zeros(2000, np.int32), None, 1000, "sweep"),       # too many rings
        lambda: C.score(xy, good, 8, 8, None, np.zeros((7, 8), np.uint8), 1000),     # gt shape
        lambda: pc.pack([[(2**40, 0)]]),
        lambda: C.score_batch([xy], [good, good], [8], [8], None, None, 1000, 2),
    ]
    for f in bad:
        with pytest.raises(OK_ERRORS):
            f()


def test_empty_and_tiny_rings():
    for rings in ([], [[]], [[(5, 5)]], [[(1, 1), (9, 9)]], [[], [(3, 3)], []]):
        r = pc.rasterize(rings, 9, 9)
        assert not r["evenodd"].any() and not r["nonzero"].any()
        v = pc.validate(rings)
        assert v["ring"].shape == (len(rings), len(pc.DEFECTS))
        assert (v["ring"][:, 0] == 1).all()


def test_huge_coordinates_exact_against_pure_oracle():
    """int32-extreme coordinates go through the __int128 path; compare to oracle (Python ints)."""
    rng = random.Random(1)
    for _ in range(300):
        pick = [I32.min, I32.min + 1, -10**9, -1, 0, 1, 37, 500, 999, 1000, 10**9, I32.max]
        rings = [[(rng.choice(pick) if rng.random() < .5 else rng.randint(-2000, 3000),
                   rng.choice(pick) if rng.random() < .5 else rng.randint(-2000, 3000))
                  for _ in range(rng.randint(1, 7))] for _ in range(rng.randint(1, 2))]
        W, H = rng.randint(1, 9), rng.randint(1, 9)
        G = rng.choice([1000, 7, 2**30])
        r = pc.rasterize(rings, W, H, G=G, winding=True)
        w, c = O.winding_and_parity(rings, W, H, G)
        assert (r["winding"] == w).all(), (rings, W, H, G)
        assert (r["evenodd"] == (c % 2 == 1)).all()
        for ring in rings:
            assert pc.ring_labels(pc.validate([ring])["ring"][0]) == O.ring_problems(ring), ring
            assert pc.signed_area2(ring) == O.twice_area(ring)


@settings(max_examples=400, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(st.lists(st.lists(st.tuples(st.integers(I32.min, I32.max), st.integers(I32.min, I32.max)), max_size=8),
                max_size=4),
       st.integers(1, 12), st.integers(1, 12))
def test_hypothesis_int32_rings(rings, W, H):
    r = pc.rasterize(rings, W, H, winding=True)
    w, _ = O.winding_and_parity(rings, W, H)
    assert (r["winding"] == w).all()
    v = pc.validate(rings)
    for k, ring in enumerate(rings):
        assert pc.ring_labels(v["ring"][k]) == O.ring_problems(ring)
