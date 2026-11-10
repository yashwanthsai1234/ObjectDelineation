"""Break-it round 1 on the CUDA engine (CPU build): edge cases that the engine or its Python wrappers get wrong.

Each test states the wanted behaviour and was xfail(strict=True) until the code was fixed. The engine itself matched
the oracles on every probed case (empty batch, 0 x N, 1 x N, one pixel, all foreground, duplicate / collinear /
frame-edge / outside / int32-extreme vertices, bow-ties, holes touching or equal to the outer ring, G = 1 and 2^30,
zero tables, sum_x / sum_y of empty masks, arena reuse after trim); the findings are at its seams.
"""
import numpy as np
import pytest

from objdel.geom import engine as gc
from objdel.geom.rings import trace
from objdel.reward.cache import ExpertCache


def _cache(fill) -> ExpertCache:
    r = np.zeros((2, 8, 8), bool)
    r[:, 2:6, 2:6] = fill
    return ExpertCache(r, 1.0, 1.0)


def test_stale_handle_is_refused_after_its_slot_is_reused():
    a, b = _cache(True), _cache(False)
    h = gc.add_tables(*a.table_args())
    gc.drop_tables(h)
    h2 = gc.add_tables(*b.table_args())
    try:
        m = np.zeros((8, 8), np.uint8)
        with pytest.raises(ValueError):
            gc.score_masks(m.ravel(), np.array([[8, 8]], np.int32), np.array([0, 64], np.int64),
                           np.array([h], np.int64))
    finally:
        gc.drop_tables(h2)


def test_cache_refuses_more_raters_than_the_engine_holds():
    with pytest.raises(ValueError):
        ExpertCache(np.ones((33, 4, 4), bool), 1.0, 1.0)


def test_trace_offset_overflow_raises():
    with pytest.raises((ValueError, OverflowError)):
        trace([np.ones((1, 1), bool)], [(2**63 - 1, 0)])


def test_trace_refuses_misshaped_offsets():
    with pytest.raises(ValueError):
        trace([np.ones((1, 1), bool), np.ones((1, 1), bool)], [(1, 2, 3, 4)])


def test_trace_refuses_fractional_offsets():
    with pytest.raises(ValueError):
        trace([np.ones((1, 1), bool)], [(0.7, 0)])
