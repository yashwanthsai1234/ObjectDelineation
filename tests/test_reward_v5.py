"""B3.2-B3.4, B3.6: reward v5 (location term, RLOO advantages) and its anti-gaming tests T1-T12, T14, T15
(docs/implementation2.0/specs/B.spec.md s5.2; T10 is tests/test_scorer_sums.py, T13 left with adv_scale)."""
import dataclasses
import json

import numpy as np
import pytest
from oracles import score_numpy as host

from objdel.reward import REWARD_V5, ExpertCache, advantages_v5, score_completion, score_mask, sig_id
from objdel.reward.advantages import _groups, _topology, rewards_v5
from objdel.reward.score import gated_loc, location
from objdel.reward.spec import RewardV5

N = 96
C = (48, 30)  # nodule centre (row, column)


def disk(r, cy, cx):
    yy, xx = np.mgrid[:N, :N]
    return (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r


@pytest.fixture(scope="module")
def cache():
    """Three raters, discs of radius 6-8 around C (one component each)."""
    return ExpertCache.build(np.stack([disk(r, *C) for r in (6, 7, 8)]))


def row(t, legal=1.0):
    """score_completion's numeric terms of a single-object answer from its ObjectTerms."""
    return dict(format_ok=1.0, legal=legal, iou_tc=t.iou_tc, iou_min=t.iou_tc, f5=t.f5, g_area=t.g_area,
                topo_any=float(t.topo_freq > 0), topo_freq=t.topo_freq, sig_id=sig_id((t.sig,)),
                loc=gated_loc(t.iou_tc, t.g_area, t.topo_freq, t.loc_raw))


def stack(rows):
    return {k: np.array([r[k] for r in rows], np.float64) for k in rows[0]}


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


# ------------------------------------------------------------------ spec (B3.3)

def test_alpha_loc_below_miss_iou_is_required():
    with pytest.raises(ValueError):
        RewardV5(c_topo=1.0, alpha_loc=0.1, miss_iou=0.1)  # a full-location miss would tie a hit (math-stats s11)
    assert REWARD_V5.alpha_loc < REWARD_V5.miss_iou and REWARD_V5.name == "reward_v5"


# ------------------------------------------------------------------ location term (B3.2)

def test_location_matches_the_oracle_and_the_sums(cache):
    m = disk(7, C[0], C[1] + 20)
    t = score_mask(m, cache)
    assert dataclasses.asdict(t) == host.score_mask(m, cache)
    ys, xs = np.nonzero(m)
    assert (t.sum_x, t.sum_y) == (int(xs.sum()), int(ys.sum()))
    d = 20.0  # centroid offset; vote centroid is C (symmetric discs)
    assert t.loc_raw == pytest.approx(np.exp(-d * d / (2 * cache.radius2)), rel=1e-12)


def test_location_is_zero_for_an_empty_mask():
    assert location(0, 0, 0, (10, 10, 1), 3.0) == 0.0


def test_t5_loc_strictly_decreases_with_distance(cache):
    ts = [score_mask(disk(7, C[0], C[1] + d), cache) for d in (16, 20, 26, 34, 44)]
    assert all(t.iou_tc == 0 and t.g_area > 0 and t.topo_freq > 0 for t in ts)
    loc = [gated_loc(t.iou_tc, t.g_area, t.topo_freq, t.loc_raw) for t in ts]
    assert all(a > b > 0 for a, b in zip(loc[:-1], loc[1:], strict=True))


def test_t2_huge_centred_outline_gets_almost_no_location(cache):
    """Under the smooth g_area a huge outline keeps g_area = A_max / A > 0, so T2's precondition g_area = 0 cannot
    be built; the gates still bound it: IoU ~ 1/k and g_area ~ 1/k give loc <= miss_iou / 4 (at k = 2 / miss_iou)."""
    for r in (20, 24, 29):  # inside the frame
        t = score_mask(disk(r, *C), cache)
        assert t.loc_raw > 0.99
        assert gated_loc(t.iou_tc, t.g_area, t.topo_freq, t.loc_raw) <= REWARD_V5.miss_iou / 4


def test_t3_speck_at_the_centroid_gets_zero(cache):
    m = np.zeros((N, N), bool)
    m[C[0]:C[0] + 2, C[1]:C[1] + 2] = True  # below the size filter: empty signature, no doctor drew it
    t = score_mask(m, cache)
    assert t.loc_raw > 0.9 and gated_loc(t.iou_tc, t.g_area, t.topo_freq, t.loc_raw) == 0.0


def test_t4_two_blobs_around_the_centroid_get_zero(cache):
    m = disk(6, C[0], C[1] - 20) | disk(6, C[0], C[1] + 20)
    t = score_mask(m, cache)
    assert t.iou_tc == 0 and t.loc_raw > 0.99 and t.topo_freq == 0
    assert gated_loc(t.iou_tc, t.g_area, t.topo_freq, t.loc_raw) == 0.0


def test_t6_illegal_or_unparseable_answers_get_zero(cache):
    caches = {"nodule": cache}
    bad = score_completion("not json", False, ("nodule",), caches, N, N)
    assert bad["loc"] == 0.0 and bad["legal"] == 0.0
    t = row(score_mask(disk(7, C[0], C[1] + 20), cache), legal=0.0)
    T = stack([t] * 8)
    assert np.all(R_of(T) == 0.0)


def test_completion_loc_matches_the_oracle(cache):
    caches = {"nodule": cache}
    x0, x1, y0, y1 = [round(v * 1000 / N) for v in (C[1] + 15, C[1] + 27, C[0] - 6, C[0] + 6)]
    text = json.dumps({"nodule": [{"outer": [[x0, y0], [x1, y0], [x1, y1], [x0, y1]], "holes": []}]},
                      separators=(",", ":"))
    got = score_completion(text, False, ("nodule",), caches, N, N)
    assert got == host.score_completion(text, False, ("nodule",), caches, N, N) and got["loc"] > 0


# ------------------------------------------------------------------ reward and advantages (B3.4)

@pytest.mark.parametrize("loc_raw", [0.0, 0.3, 1.0])
@pytest.mark.parametrize("g_area", [0.2, 1.0])
def test_t1_r_continuous_at_miss_iou_and_non_decreasing(loc_raw, g_area):
    m = REWARD_V5.miss_iou
    iou = np.r_[np.linspace(0, 1, 401), m - 1e-12, m]
    iou.sort()
    loc = [gated_loc(x, g_area, 1.0, loc_raw) for x in iou]
    R = R_of(group(iou_tc=np.resize(iou, 8 * 51)[:408], loc=np.resize(loc, 408)[:408], g_area=np.full(408, g_area),
                   f5=np.full(408, 0.5), legal=np.ones(408), topo_any=np.zeros(408), topo_freq=np.full(408, 1.0),
                   sig_id=np.full(408, 7.0), format_ok=np.ones(408)))
    R = R[:iou.size]
    assert np.all(np.diff(R) >= -1e-12)
    j = np.searchsorted(iou, m)
    assert abs(R[j] - R[j - 1]) < 1e-9


def test_t11_hit_beats_any_miss():
    m = REWARD_V5.miss_iou
    hits = np.linspace(m, 1, 8)
    misses = np.linspace(0, m - 1e-9, 8)
    Rh = R_of(group(iou_tc=hits, loc=np.zeros(8), topo_any=np.zeros(8)))
    Rm = R_of(group(iou_tc=misses, loc=[gated_loc(x, 1.0, 1.0, 1.0) for x in misses], topo_any=np.zeros(8)))
    assert Rh.min() > Rm.max()


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
