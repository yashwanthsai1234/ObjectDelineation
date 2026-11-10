"""Break-it round 2 on the CUDA engine's CPU build (checker.rings.trace, reward.cache): each test shows one defect,
now fixed."""
import numpy as np
import pytest

from objdel.geom.rings import polygonize, trace
from objdel.verify.reward.cache import ExpertCache


def test_complex_offset_refused():
    with pytest.raises((ValueError, TypeError)):
        polygonize(np.ones((2, 2), bool), (1 + 1j, 0))


@pytest.mark.parametrize("t5", [float("nan"), 0.0, float("inf"), -1.0])
def test_bad_t5_refused_at_construction(t5):
    st = np.zeros((2, 6, 6), bool)
    st[:, 1:5, 1:5] = True
    with pytest.raises(ValueError):
        ExpertCache(st, 1.0, t5)


def test_empty_batch_with_empty_offsets():
    assert trace([], []) == []
