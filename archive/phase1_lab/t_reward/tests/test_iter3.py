"""Iteration 3: single-rater policy (R5), signature size filter (R6), IoU_tc vote fallback (R3)."""
import pickle
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from rewardcore import fixtures, perturb
from rewardcore.reward import ExpertCache, rl_eligible, score, signature_min_area
from rewardcore.topology import filtered_signature

ROOT = Path(__file__).resolve().parents[1]


def disk(n, r, cy=None, cx=None):
    y, x = np.mgrid[:n, :n]
    return np.hypot(y - (n / 2 if cy is None else cy), x - (n / 2 if cx is None else cx)) <= r


# ---------------- R5: single-rater policy and deterministic k = 1 cache ----------------

def test_single_rater_items_are_excluded_from_rl():
    st = fixtures.load("lidc")
    k2 = next(v for v in st.values() if len(v) == 2)
    assert rl_eligible(k2) and not rl_eligible(k2[:1])


def _cache_fingerprint(c):
    return (c.tau, c.tol4, c.tol5, c.perimeter, c.min_area, tuple(c.rater_sigs),
            c.dist_out_union.tobytes(), c.dist_in_core.tobytes(), c.vote.tobytes())


def test_single_rater_cache_is_deterministic():
    m = next(iter(fixtures.load("lidc").values()))[:1]
    prior = dict(tau_prior=lambda r: 0.75 * r ** 0.55, sigma_prior=lambda r: 0.48 * r ** 0.40)
    a, b = ExpertCache.build(m, **prior), ExpertCache.build(m, **prior)
    assert _cache_fingerprint(a) == _cache_fingerprint(b)
    assert _cache_fingerprint(pickle.loads(pickle.dumps(a))) == _cache_fingerprint(a)
    code = ("import sys; sys.path.insert(0, %r)\n"
            "from rewardcore import fixtures\nfrom rewardcore.reward import ExpertCache\n"
            "m = next(iter(fixtures.load('lidc').values()))[:1]\n"
            "c = ExpertCache.build(m, tau_prior=lambda r: 0.75 * r ** 0.55, sigma_prior=lambda r: 0.48 * r ** 0.40)\n"
            "print(repr((c.tau, c.tol4, c.tol5, c.perimeter, c.min_area)))\n" % str(ROOT))
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == repr((a.tau, a.tol4, a.tol5, a.perimeter, a.min_area))


# ---------------- R6: one signature size filter everywhere ----------------

def test_filter_rule():
    assert signature_min_area(100) == 5 and signature_min_area(20000) == 100
    assert signature_min_area(0) == 5


@pytest.fixture
def raters():
    return np.stack([disk(60, 18), disk(60, 19), disk(60, 20)])


def _stamp(m, area, cy, cx, value):
    off = perturb._disk_offsets(area)
    out = m.copy()
    out[cy + off[:, 0], cx + off[:, 1]] = value
    return out


def test_sub_filter_speck_on_a_rater_does_not_change_the_topology_set(raters):
    base = ExpertCache.build(raters)
    A = base.min_area  # max(5, ceil(0.5 % of |U| = 1257)) = 7
    assert A == 7
    noisy = raters.copy()
    noisy[0] = _stamp(noisy[0], A - 1, 4, 4, True)  # detached speck below the filter
    noisy[1] = _stamp(noisy[1], A - 1, 30, 30, False)  # pinhole below the filter
    assert ExpertCache.build(noisy).min_area == A
    assert ExpertCache.build(noisy).signatures == base.signatures == {(0,): 3}


def test_above_filter_component_changes_the_topology_set(raters):
    A = ExpertCache.build(raters).min_area
    noisy = raters.copy()
    noisy[0] = _stamp(noisy[0], A, 4, 4, True)  # detached part at the filter size
    noisy[1] = _stamp(noisy[1], A, 30, 30, False)  # hole at the filter size
    assert ExpertCache.build(noisy).signatures == {(0,): 1, (0, 0): 1, (1,): 1}


def test_topo_any_invariant_to_sub_filter_specks(raters):
    cache = ExpertCache.build(raters)
    rng = np.random.default_rng(0)
    pred = raters[1].copy()
    t0 = score(pred, cache)
    for _ in range(50):
        p = pred.copy()
        a = int(rng.integers(1, cache.min_area))
        sp = perturb.add_speck(p, a, float(rng.uniform(4, 8)), rng)
        hole = perturb.punch_hole(p, int(rng.integers(1, cache.min_area)), rng)
        for q in (sp, hole):
            if q is None:
                continue
            t = score(q, cache)
            assert (t.topo_any, t.sig) == (t0.topo_any, t0.sig)
            assert t.sig_raw != t0.sig_raw  # still logged
    big = perturb.add_speck(pred, cache.min_area, 6.0, rng)
    assert score(big, cache).topo_any == 0.0


def test_filtered_signature_matches_brute_force_filter():
    from objdel.checker.topology import size_filter as repo_filter

    from rewardcore import synth
    from rewardcore.topology import instance_topology_full

    rng = np.random.default_rng(1)
    for _ in range(500):
        _, m = synth.random_mask(rng, 8, 64)
        a = int(rng.integers(1, 20))
        assert filtered_signature(m, a) == instance_topology_full(repo_filter(m, a)).signature()


# ---------------- R3: IoU_tc fallback ----------------

def test_iou_tc_falls_back_to_vote_for_a_topology_no_rater_drew(raters):
    cache = ExpertCache.build(raters)
    wrong = perturb.punch_hole(raters[0], 3 * cache.min_area, np.random.default_rng(2), centre=True)
    t = score(wrong, cache)
    assert t.topo_any == 0.0 and t.iou_tc == t.iou  # not max over raters (which would be ~1 vs rater 0)


# ---------------- iteration 4: area guard ----------------

def test_area_guard_is_one_inside_the_rater_range_and_penalizes_hedges(raters):
    from rewardcore.reward import area_guard

    cache = ExpertCache.build(raters)
    a_min, a_max = min(cache.rater_areas), max(cache.rater_areas)
    for m in raters:
        assert score(m, cache).g_area == 1.0
    assert score(cache.vote, cache).g_area == 1.0
    t_union = score(cache.union, cache)
    assert t_union.g_area == pytest.approx(a_max / cache.n_union)
    t_core = score(cache.core, cache)
    assert t_core.g_area == pytest.approx(np.count_nonzero(cache.core) / a_min)
    assert area_guard(2 * a_max, cache) == pytest.approx(0.5) and area_guard(a_min // 2, cache) < 0.51
    assert area_guard(0, cache) == 0.0
    assert t_union.g_area_tauP == 1.0  # the tau*P tolerance never touches the union


def test_reward_v4_term_guard_ranks_inflation_lower():
    from rewardcore.advantages import dm_loo_topology, reward_v4

    rng = np.random.default_rng(0)
    R = np.stack([rng.uniform(0.85, 0.9, (1, 8)), rng.uniform(0.9, 0.95, (1, 8)), np.ones((1, 8)), np.ones((1, 8))], -1)
    R[0, 0, 3] = 0.8  # rollout 0 is 25 % larger than the largest rater
    r = dm_loo_topology(np.zeros((1, 8)), np.ones((1, 8)), np.ones((1, 8), bool))
    names = ("iou", "band", "legal", "g_area")
    W, F = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25])
    a_term, c_term = reward_v4(R, W, F, names, r, guard="term")
    a_gate, _ = reward_v4(R, W, F, names, r, guard="gate")
    R1 = R.copy()
    R1[0, 0, 3] = 1.0
    a_none, _ = reward_v4(R1, W, F, names, r, guard="term")
    assert a_term[0, 0] < a_none[0, 0] and a_gate[0, 0] < a_none[0, 0]
    assert c_term.shape[-1] == 5  # iou, band, legal, g_area, topology
