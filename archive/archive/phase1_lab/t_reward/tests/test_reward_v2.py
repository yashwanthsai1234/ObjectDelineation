"""Iteration 2: topology-conditional IoU, DM-loo topology term, reward-v2 advantages."""
import numpy as np
import pytest

from rewardcore import fixtures
from rewardcore.advantages import dm_loo_bias, dm_loo_topology, reward_v2
from rewardcore.reward import ExpertCache, score, signature_id
from rewardcore.topology import instance_topology_full

W2, F2, NAMES2 = np.array([0.3, 0.3, 0.1]), np.array([0.05, 0.05, 0.25]), ("iou", "band", "legal")


def ring_and_disk():
    y, x = np.mgrid[:40, :40]
    d = np.hypot(y - 19.5, x - 19.5)
    disk = d <= 12
    ring = disk & (d > 3)
    return disk, ring


def test_topology_conditional_iou_compares_with_same_topology_raters():
    disk, ring = ring_and_disk()
    cache = ExpertCache.build(np.stack([disk, disk, disk, ring]))  # 3:1, the vote has no hole
    t_ring, t_disk = score(ring, cache), score(disk, cache)
    assert t_ring.iou < 0.95 and t_ring.iou_tc == 1.0  # vote-IoU penalizes the minority hole, IoU_tc does not
    assert t_disk.iou_tc == 1.0
    speck = disk.copy()
    speck[2, 2] = True  # topology no rater has: falls back to all raters
    assert score(speck, cache).iou_tc == pytest.approx(disk.sum() / speck.sum())


def test_dm_loo_values():
    sig = np.array([[0, 0, 0, 0, 0, 0, 1, 2]])
    freq = np.array([[0.75] * 6 + [0.25, 0.0]])
    valid = np.array([[1] * 7 + [0]], bool)
    r = dm_loo_topology(sig, freq, valid, lam=0.5)
    assert r[0, 7] == 0.0
    assert r[0, 0] == pytest.approx(1 + 0.5 * (np.log(0.75) - dm_loo_bias(0.75, 8) - np.log(5 / 7)))
    assert r[0, 6] == pytest.approx(1 + 0.5 * (np.log(0.25) - dm_loo_bias(0.25, 8) - np.log(0.5 / 7)))
    raw = dm_loo_topology(sig, freq, valid, lam=0.5, debias=False)
    assert raw[0, 6] == pytest.approx(1 + 0.5 * (np.log(0.25) - np.log(0.5 / 7)))


@pytest.mark.parametrize("freqs, lam", [((0.75, 0.25), 0.5), ((0.5, 0.5), 0.5), ((0.5, 0.25, 0.25), 0.5),
                                        ((2 / 3, 1 / 3), 0.5), ((5 / 6, 1 / 6), 0.25)])
def test_debiased_topology_reward_has_equal_expectation_per_class_at_target(freqs, lam):
    """Exact enumeration over all groups of G = 8 drawn from the rater distribution (pi = f): the
    expected topology reward of a rollout is the same for every rater answer (no class favoured).
    Holds where the eps clip does not bind (f = 1/6 at lam = 0.5 hits it: bracket 0.068 < 0.1)."""
    from itertools import product
    from math import prod

    G, K = 8, len(freqs)
    exp_r = np.zeros(K)
    for combo in product(range(K), repeat=G):
        p = prod(freqs[c] for c in combo)
        sig = np.array([combo])
        f = np.array([[freqs[c] for c in combo]])
        r = dm_loo_topology(sig, f, np.ones((1, G), bool), lam=lam)[0]
        for c in range(K):
            m = np.array(combo) == c
            exp_r[c] += p * r[m].sum()
    per_rollout = exp_r / (G * np.array(freqs))  # E[r | answer c]
    assert np.ptp(per_rollout) < 1e-9


def test_topology_term_is_mean_centred_not_z_scored():
    rng = np.random.default_rng(0)
    R = np.stack([rng.uniform(0.8, 0.9, (4, 8)), rng.uniform(0.8, 0.9, (4, 8)), np.ones((4, 8))], -1)
    r_topo = np.ones((4, 8))
    r_noisy = r_topo + rng.uniform(0, 1e-4, r_topo.shape)
    a0, _ = reward_v2(R, W2, F2, NAMES2, r_topo, c_topo=2.0)
    a1, _ = reward_v2(R, W2, F2, NAMES2, r_noisy, c_topo=2.0)
    assert np.abs(a1 - a0).max() < 1e-3  # no dark-room amplification
    R[:, 0, 2] = 0.0  # illegal -> gated out of the topology term too
    a2, c2 = reward_v2(R, W2, F2, NAMES2, r_topo, c_topo=2.0)
    assert np.all(c2[:, 0, -1] < 0) and np.allclose(a2.mean(1), a2.mean())


def test_no_topology_class_favoured_for_rater_groups():
    """Groups drawn i.i.d. from the raters' own masks: v2 advantage per class ~ 0 (v1 favours the vote's class)."""
    items = []
    for _, _, st in fixtures.load_all():
        sigs = [instance_topology_full(m).signature() for m in st]
        if len(set(sigs)) > 1:
            cache = ExpertCache.build(st)
            mode = max(set(sigs), key=sigs.count)
            items.append(([score(m, cache) for m in st], [s != mode for s in sigs]))
    assert len(items) >= 10
    rng = np.random.default_rng(1)
    diff = []
    for terms, minority in items:
        a_min, a_mode = [], []
        for _ in range(200):
            d = rng.integers(0, len(terms), (1, 8))
            T = [[terms[j] for j in d[0]]]
            sig = np.array([[signature_id(t.sig) for t in T[0]]])
            freq = np.array([[t.topo_freq for t in T[0]]])
            valid = np.array([[t.topo_any > 0 and t.iou_tc >= 0.5 for t in T[0]]])
            R = np.array([[[t.iou_tc, t.band["F5"], 1.0] for t in T[0]]])
            a, _ = reward_v2(R, W2, F2, NAMES2, dm_loo_topology(sig, freq, valid), c_topo=2.0)
            for g, j in enumerate(d[0]):
                (a_min if minority[j] else a_mode).append(a[0, g])
        diff.append(np.mean(a_min) - np.mean(a_mode))
    assert abs(np.mean(diff)) < 0.15  # measured +0.08 (slight pro-rare lean of the LOO estimator)


@pytest.mark.parametrize("lam", [0.25, 0.5, 1.0, 2.0])
def test_valid_always_beats_invalid_for_all_reachable_states(lam):
    """Every reachable (f, n): f in rater shares j/k (k = 2..6), n(a) in 1..G (G = 8). The topology
    reward of a valid rollout is strictly above that of an invalid one (0), and so is its
    mean-centred advantage contribution in any group mixing the two."""
    G = 8
    shares = sorted({j / k for k in range(2, 7) for j in range(1, k + 1)})
    for f in shares:
        for n in range(1, G):  # n rollouts with the valid answer a, G - n invalid
            sig = np.array([[0] * n + list(range(1, G - n + 1))])
            freq = np.array([[f] * n + [0.0] * (G - n)])
            valid = np.array([[True] * n + [False] * (G - n)])
            r = dm_loo_topology(sig, freq, valid, lam=lam)
            assert r[0, :n].min() > r[0, n:].max() == 0.0
            centred = r - r.mean(1, keepdims=True)
            assert centred[0, :n].min() > centred[0, n:].max()
        r = dm_loo_topology(np.zeros((1, G), int), np.full((1, G), f), np.ones((1, G), bool), lam=lam)
        assert r.min() >= 0.1  # n = G: the clip binds only when lam (log f - log pi) < -0.9
