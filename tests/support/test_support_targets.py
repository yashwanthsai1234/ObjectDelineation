"""The helpers of tests/support/targets.py agree with the spelled-out library calls they replace."""
import numpy as np

from objdel.data.targets.build import build_target
from objdel.verify.reward.cache import ExpertCache
from tests.support.shapes import disc
from tests.support.targets import mask_text, rater_texts


def test_expert_cache_min_area_is_min_areas_for_one_name():
    stack = np.stack([disc(50, 50, 25, 25, 10), disc(50, 50, 25, 26, 11), disc(50, 50, 24, 25, 9)])
    ma = {"nodule": ExpertCache.build(stack).min_area}
    expected = [build_target({"nodule": m}, ("nodule",), min_area=ma).text for m in stack]
    assert rater_texts({"nodule": stack}) == expected


def test_mask_text_is_a_one_rater_rater_text():
    m = disc(50, 50, 25, 25, 10)
    assert [mask_text(m)] == rater_texts({"nodule": m[None]})
