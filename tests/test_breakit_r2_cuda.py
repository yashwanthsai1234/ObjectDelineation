"""Break-it round 2 on the CUDA engine's CPU build (checker.rings.trace, reward.cache): each test shows one defect and
is xfail(strict) until it is fixed."""
import numpy as np
import pytest

from objdel.checker import polygonize, trace
from objdel.reward.cache import ExpertCache


@pytest.mark.xfail(strict=True, reason="BR2-cuda-1: trace accepts complex offsets and drops the imaginary part "
                                       "silently (the BR1-cuda-5 whole-number check passes 1+1j)")
def test_complex_offset_refused():
    with pytest.raises((ValueError, TypeError)):
        polygonize(np.ones((2, 2), bool), (1 + 1j, 0))


@pytest.mark.xfail(strict=True, reason="BR2-cuda-2: ExpertCache accepts t5 = NaN, 0 or inf; the engine refuses it "
                                       "only at the first device() call, inside score_many mid-step")
@pytest.mark.parametrize("t5", [float("nan"), 0.0, float("inf"), -1.0])
def test_bad_t5_refused_at_construction(t5):
    st = np.zeros((2, 6, 6), bool)
    st[:, 1:5, 1:5] = True
    with pytest.raises(ValueError):
        ExpertCache(st, 1.0, t5)


@pytest.mark.xfail(strict=True, reason="BR2-cuda-3: trace([], []) raises ValueError: an empty batch with an empty "
                                       "offsets list is refused (shape (0,) != (0, 2))")
def test_empty_batch_with_empty_offsets():
    assert trace([], []) == []
