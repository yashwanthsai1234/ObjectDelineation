"""B3.4 (reward v5 advantages): RLOO, linearity, identical answers, drop rule; T7-T9, T12, T14, T15
(docs/implementation2.0/specs/B.spec.md s5.2). Its `group` and `R_of` helpers serve test_reward_v5_shaping too."""
import dataclasses

import numpy as np

from objdel.verify.reward.advantages import _groups, _topology, advantages_v5, rewards_v5
from objdel.verify.reward.spec import REWARD_V5


def group(**kw):
    """A group of 8 synthetic answers: legal single misses, one topology a doctor drew; kw overrides per term."""
    base = dict(format_ok=1.0, legal=1.0, iou_tc=0.0, iou_min=0.0, f5=0.5, g_area=1.0, topo_any=1.0, topo_freq=1.0,
                sig_id=7.0, loc=0.0)
    T = {k: np.full(8, v, np.float64) for k, v in base.items()}
    for k, v in kw.items():
        T[k] = np.asarray(v, np.float64)
    if "iou_min" not in kw:
        T["iou_min"] = T["iou_tc"]
    return T


def R_of(T, spec=REWARD_V5):
    G = _groups(T, len(T["iou_tc"]), ("loc",))
    return rewards_v5(G, _topology(G, len(T["iou_tc"]), spec), spec).reshape(-1)


def test_t7_t8_rloo_sums_to_zero_and_is_linear():
    rng = np.random.default_rng(0)
    T = group(iou_tc=rng.uniform(0, 0.09, 8), loc=rng.uniform(0, 1, 8), f5=rng.uniform(0, 1, 8),
              g_area=rng.uniform(0, 1, 8), topo_any=np.zeros(8))
    a = advantages_v5(T)
    assert abs(a.sum()) < 1e-9 and np.ptp(a) > 0
    R = R_of(T)
    assert np.allclose(a, 8 / 7 * (R - R.mean()), atol=1e-12)
    spec2 = dataclasses.replace(REWARD_V5, alpha_loc=3 * REWARD_V5.alpha_loc, miss_iou=1.0, w_f5=3 * REWARD_V5.w_f5,
                                w_area=3 * REWARD_V5.w_area)
    T3 = dict(T, iou_tc=3 * T["iou_tc"])  # every term of R times 3 (r_topo = 0)
    assert np.allclose(advantages_v5(T3, spec=spec2), 3 * a, atol=1e-12)


def test_t9_identical_answers_get_zero_without_nan():
    a = advantages_v5(group(iou_tc=np.full(8, 0.6), topo_any=np.ones(8)))
    assert np.all(a == 0) and np.all(np.isfinite(a))


def test_t12_closest_miss_gets_the_largest_positive_advantage():
    r2, ramp = 7.0**2, 1 - 0.02 / REWARD_V5.miss_iou
    d = np.array([5, 10, 20, 20, 20, 20, 20, 20], float)
    loc = ramp * np.exp(-d**2 / (2 * r2))
    a = advantages_v5(group(iou_tc=np.full(8, 0.02), loc=loc, topo_any=np.zeros(8)))
    assert a[0] > 0 and a[0] == a.max() and a[0] > a[1] > a[2]


def test_t14_far_miss_group_is_not_dropped():
    a = advantages_v5(group(loc=[0.05, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01], topo_any=np.zeros(8)))
    assert a[0] > 0 and np.any(a != 0)


def test_t15_drop_rule_tests_legal():
    a = advantages_v5(group(iou_tc=np.full(8, 0.3), legal=[1, 1, 1, 1, 1, 1, 0, 0], topo_any=np.zeros(8)))
    assert np.all(a[6:] < 0) and np.all(a[:6] > 0)
