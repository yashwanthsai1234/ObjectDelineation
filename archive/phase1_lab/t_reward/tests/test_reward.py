"""Task 3: set-valued topology term and the scalar reward."""
import numpy as np
import pytest

from rewardcore import perturb
from rewardcore.reward import ExpertCache, RewardConfig, reward, score


def blob(h, w, boxes):
    m = np.zeros((h, w), bool)
    for y0, y1, x0, x1 in boxes:
        m[y0:y1, x0:x1] = True
    return m


@pytest.fixture
def mixed_raters():
    """Two raters draw one part, one rater draws two parts."""
    one = blob(40, 60, [(10, 30, 10, 50)])
    two = blob(40, 60, [(10, 30, 10, 28), (10, 30, 31, 50)])
    return ExpertCache.build(np.stack([one, one, two])), one, two


def test_set_valued_topology(mixed_raters):
    cache, one, two = mixed_raters
    t1, t2 = score(one, cache), score(two, cache)
    assert (t1.topo_any, t1.topo_freq) == (1.0, pytest.approx(2 / 3))
    assert (t2.topo_any, t2.topo_freq) == (1.0, pytest.approx(1 / 3))
    holed = one.copy()
    holed[20, 30] = False  # 1-px pinhole: below the signature size filter (R6), not a topology change
    t3 = score(holed, cache)
    assert t3.topo_any == 1.0 and (t3.b0, t3.b1) == (1, 1) and t3.sig_raw == (1,)
    holed[18:21, 29:32] = False  # 9 px >= min_area (5): a real hole no rater drew
    t4 = score(holed, cache)
    assert cache.min_area == 5 and t4.topo_any == 0.0 and t4.sig == (1,)


def test_topology_gate_and_weights(mixed_raters):
    cache, one, _ = mixed_raters
    far = perturb.shift(one, 0, 25)  # IoU < 0.5: no topology credit
    t = score(far, cache)
    assert t.iou < 0.5 and t.topo_any == 1.0
    cfg = RewardConfig(band_formula="F5")
    assert reward(t, cfg) == pytest.approx(cfg.w_iou * t.iou + cfg.w_band * t.band["F5"])
    assert reward(score(one, cache), cfg) == pytest.approx(1.0)


def test_empty_prediction_and_shape_check(mixed_raters):
    cache, one, _ = mixed_raters
    t = score(np.zeros_like(one), cache)
    assert t.iou == 0.0 and t.topo_any == 0.0 and t.area == 0
    with pytest.raises(ValueError):
        score(np.zeros((10, 10), bool), cache)


def test_single_rater_cache_uses_priors():
    m = blob(50, 50, [(10, 40, 10, 40)])
    c0 = ExpertCache.build(m[None])
    c1 = ExpertCache.build(m[None], tau_prior=lambda r: 3.0, sigma_prior=lambda r: 1.0)
    assert np.isnan(c0.tau_raw) and c0.tau == 1.0 and c0.tol5 == c0.tol4
    assert c1.tau == 3.0 and c1.tol5 != c1.tol4
