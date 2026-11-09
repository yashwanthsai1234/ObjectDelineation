"""K1/K2-style differential tests: fast kernels == exact oracle."""
import random

import numpy as np
import oracle as O
import pytest

from targets import fastgeom as F


def _shapes(rng):
    kind = rng.randrange(5)
    if kind == 0:
        return [O.random_simple_polygon(rng, rng.randint(3, 14))]
    if kind == 1:
        return [O.random_polygon_soup(rng, rng.randint(3, 10))]
    if kind == 2:
        return [O.random_grid_polygon(rng, rng.randint(3, 10), step=rng.choice([20, 50, 100, 125, 200]))]
    if kind == 3:
        return list(O.annulus(rng.randint(300, 700), rng.randint(300, 700), rng.randint(150, 290), rng.randint(20, 140),
                              rng.randint(4, 16), hole_ccw=rng.random() < 0.5))
    return [O.random_grid_polygon(rng, rng.randint(3, 8), step=rng.choice([25, 50])) for _ in range(rng.randint(2, 3))]


@pytest.mark.parametrize("seed", range(6))
def test_raster_matches_oracle(seed):
    rng = random.Random(seed)
    for _ in range(250):
        rings = _shapes(rng)
        W = rng.choice([1, 3, 8, 10, 20, 25, 37, 40, 50, 64])
        H = rng.choice([1, 5, 10, 23, 40, 50])
        for rule in ("evenodd", "nonzero"):
            ref = O.rasterize(rings, W, H, rule=rule)
            got = F.raster(rings, W, H, rule=rule)
            assert np.array_equal(ref, got), (rings, W, H, rule)


@pytest.mark.parametrize("seed", range(4))
def test_ring_problems_match_oracle(seed):
    rng = random.Random(100 + seed)
    for _ in range(2500):
        k = rng.randrange(4)
        if k == 0:
            ring = O.random_simple_polygon(rng, rng.randint(3, 12))
        elif k == 1:
            ring = O.random_polygon_soup(rng, rng.randint(2, 9))
        elif k == 2:
            ring = O.random_grid_polygon(rng, rng.randint(2, 9), step=rng.choice([100, 200, 250, 500]))
        else:  # staircase-ish lattice rings with repeats and spikes
            ring = [(rng.randint(0, 4), rng.randint(0, 4)) for _ in range(rng.randint(3, 9))]
        assert F.ring_problems(ring) == O.ring_problems(ring), ring
