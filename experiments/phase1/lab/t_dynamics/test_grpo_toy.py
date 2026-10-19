import sys

import numpy as np
import pytest
from hypothesis import given, settings, strategies as st
from hypothesis.extra.numpy import arrays

import grpo_toy as gt

REPO = "/Users/sreeyashwanthsaivenkatesh/Workspace/projects/CUDA/ObjectDelineation/src"


@given(arrays(float, (3, 4, 8), elements=st.floats(-5, 5)), st.sampled_from(["std", "mean"]))
@settings(max_examples=200, deadline=None)
def test_advantages_sum_to_zero_per_group(r, mode):
    A = gt.standardize(r, mode=mode)
    assert np.allclose(A.sum(-1), 0, atol=1e-9)
    if mode == "std":
        sd = r.std(-1)
        nz = sd > 1e-6
        assert np.allclose(A.std(-1)[nz], 1, atol=1e-6)


def test_zero_variance_group_gives_zero_advantage():
    A = gt.standardize(np.ones((1, 1, 8)))
    assert np.all(A == 0)


def test_tiny_variance_is_amplified_dark_room():
    r = np.ones(8)
    r[0] += 1e-6
    A = gt.standardize(r)
    assert abs(A[0]) > 2  # 1e-6 reward difference becomes an order-1 advantage


def test_gdpo_matches_repo_group_advantages():
    sys.dont_write_bytecode = True
    sys.path.insert(0, REPO)
    from objdel.reward.annotators import group_advantages

    class S:
        def __init__(s, iou, band, topo):
            s.iou, s.band, s.topo_near = iou, band, topo

    rng = np.random.default_rng(0)
    for _ in range(50):
        iou = rng.uniform(0.3, 1, 8)
        band = rng.uniform(0, 1, 8)
        topo = rng.integers(0, 2, 8).astype(float)
        ref = group_advantages([S(*x) for x in zip(iou, band, topo)])
        ours = gt.gdpo_advantages(iou[None], band[None], topo[None])[0]
        assert np.allclose(ref, ours)


def _surrogate(theta, a, A, logp_old, ref, cfg):
    pi = gt.softmax(theta)
    lp = np.take_along_axis(np.log(pi), a, -1)
    rho = np.exp(lp - logp_old)
    obj = np.minimum(rho * A, np.clip(rho, 1 - cfg.eps_low, 1 + cfg.eps_high) * A)
    r = np.exp(np.take_along_axis(np.broadcast_to(np.log(ref), theta.shape), a, -1) - lp)
    obj = obj - cfg.beta * (r - np.log(r) - 1)
    return obj.mean(-1).sum()


@pytest.mark.parametrize("beta", [0.0, 0.1])
@pytest.mark.parametrize("shift", [0.0, 0.3, 1.5])
def test_gradient_matches_finite_differences(beta, shift):
    rng = np.random.default_rng(1)
    cfg = gt.Cfg(beta=beta)
    theta_old = rng.standard_normal((2, 3, gt.K))
    pi_old = gt.softmax(theta_old)
    a = gt.sample(pi_old, 8, rng)
    logp_old = np.take_along_axis(np.log(pi_old), a, -1)
    A = gt.standardize(rng.standard_normal(a.shape))
    ref = gt.softmax(rng.standard_normal((3, gt.K)))
    theta = theta_old + shift * rng.standard_normal(theta_old.shape)
    g = gt.grad_surrogate(theta, a, A, logp_old, ref, cfg)
    h = 1e-6
    num = np.zeros_like(theta)
    for idx in np.ndindex(theta.shape):
        tp, tm = theta.copy(), theta.copy()
        tp[idx] += h
        tm[idx] -= h
        num[idx] = (_surrogate(tp, a, A, logp_old, ref, cfg) - _surrogate(tm, a, A, logp_old, ref, cfg)) / (2 * h)
    assert np.allclose(g, num, atol=1e-5)


def test_exact_kl_gradient_matches_fd():
    rng = np.random.default_rng(2)
    theta = rng.standard_normal((1, 1, gt.K))
    ref = gt.softmax(rng.standard_normal((1, gt.K)))
    cfg = gt.Cfg(beta=1.0, kl="exact")
    a = np.zeros((1, 1, 8), int)
    g = gt.grad_surrogate(theta, a, np.zeros((1, 1, 8)), np.log(gt.softmax(theta))[..., [0] * 8], ref, cfg)

    def negkl(t):
        p = gt.softmax(t)
        return -(p * (np.log(p) - np.log(ref))).sum()

    h = 1e-6
    num = np.array([(negkl(theta + h * np.eye(gt.K)[k]) - negkl(theta - h * np.eye(gt.K)[k])) / (2 * h)
                    for k in range(gt.K)])
    assert np.allclose(g[0, 0], num, atol=1e-6)


def test_update_raises_prob_of_positive_advantage_answers():
    rng = np.random.default_rng(3)
    cfg = gt.Cfg(reward="set", lr=0.3, mu=1)
    freq = gt.hist({(1, 0): 1})[None]
    theta = np.zeros((1, 1, gt.K))
    bias = np.zeros((1, gt.K))
    p0 = gt.softmax(theta)[0, 0, gt.IDX[(1, 0)]]
    a, A, _ = gt.step(theta, bias, freq, gt.softmax(np.zeros((1, gt.K))), cfg, rng)
    p1 = gt.softmax(theta)[0, 0, gt.IDX[(1, 0)]]
    if (A > 0).any():
        assert p1 > p0


def test_clipping_blocks_gradient():
    cfg = gt.Cfg()
    theta = np.zeros((1, 1, gt.K))
    theta[0, 0, 0] = 2.0  # pi(0) now much larger than pi_old(0)
    a = np.zeros((1, 1, 8), int)
    logp_old = np.full((1, 1, 8), np.log(1 / gt.K))
    g = gt.grad_surrogate(theta, a, np.ones((1, 1, 8)), logp_old, None, cfg)
    assert np.allclose(g, 0)


def test_kl_pulls_toward_reference_when_no_advantage():
    rng = np.random.default_rng(4)
    ref = gt.softmax(rng.standard_normal((1, gt.K)))
    theta = np.zeros((4, 1, gt.K))
    cfg = gt.Cfg(beta=1.0, lr=0.5, mu=1)
    for _ in range(3000):
        pi = gt.softmax(theta)
        a = gt.sample(pi, 64, rng)
        lp = np.take_along_axis(np.log(pi), a, -1)
        theta += cfg.lr * gt.grad_surrogate(theta, a, np.zeros(a.shape), lp, ref, cfg)
    assert 0.5 * np.abs(gt.softmax(theta).mean(0)[0] - ref[0]).sum() < 0.05


def test_dm_reward_is_flat_at_target():
    f = gt.hist({(1, 0): 3, (1, 1): 1})
    a = np.array([[[gt.IDX[(1, 0)]] * 6 + [gt.IDX[(1, 1)]] * 2]])
    lp = np.log(f[a])
    r = gt.topo_reward("dm", a, f[None], lp)
    assert np.allclose(gt.standardize(r), 0)


def test_rich_get_richer_under_set_reward():
    """Expected GRPO update with two equally rewarded valid answers widens their gap."""
    rng = np.random.default_rng(5)
    S = 4000
    f = gt.hist({(1, 0): 2, (1, 1): 2})[None]
    p = gt.sft_policy(f, 0.2)
    theta = np.broadcast_to(np.log(p), (S, 1, gt.K)).copy()
    cfg = gt.Cfg(reward="set", lr=0.1, mu=1)
    theta[:, :, gt.IDX[(1, 0)]] += 0.3  # (1,0) slightly more likely
    before = theta[:, 0, gt.IDX[(1, 0)]] - theta[:, 0, gt.IDX[(1, 1)]]
    gt.step(theta, np.zeros((S, gt.K)), f, None, cfg, rng)
    after = theta[:, 0, gt.IDX[(1, 0)]] - theta[:, 0, gt.IDX[(1, 1)]]
    d = after - before
    assert d.mean() > 3 * d.std() / np.sqrt(S)


def test_metrics_zero_at_target():
    f = gt.hist({(1, 0): 2, (1, 1): 1, (2, 0): 1})
    m = gt.metrics(f[None], f[None])
    assert np.isclose(m["tv"], 0) and np.isclose(m["js"], 0, atol=1e-12) and np.isclose(m["minret"], 1)


def test_near_matches_repo_formula():
    rng = np.random.default_rng(6)
    for name, h in gt.SUITE.items():
        f = gt.hist(h)
        a = rng.integers(0, gt.K, (1, 1, 8))
        r = gt.topo_reward("near", a, f[None], None)
        for g in range(8):
            b0, b1 = gt.ANSWERS[a[0, 0, g]]
            want = max(1 / (1 + abs(b0 - t0) + abs(b1 - t1)) for t0, t1 in h)
            assert np.isclose(r[0, 0, g], want)


def test_coverage_reward_sum_is_lambda_times_distinct_valid():
    f = gt.hist({(1, 0): 2, (1, 1): 2})
    a = np.array([[[gt.IDX[(1, 0)]] * 5 + [gt.IDX[(1, 1)]] * 2 + [gt.IDX[(0, 0)]]]])
    r = gt.topo_reward("cov", a, f[None], None, lam=0.5)
    assert np.isclose(r.sum(), 7 + 0.5 * 2)


def test_invalid_answers_get_lowest_reward_in_all_variants():
    f = gt.hist({(1, 0): 3, (1, 1): 1})
    rng = np.random.default_rng(7)
    for name in ["set", "near", "freq", "cov", "uniq", "share", "dm", "set+dm", "set+dm_emp", "set+dm_loo"]:
        for _ in range(50):
            a = rng.choice([gt.IDX[(1, 0)], gt.IDX[(1, 1)], gt.IDX[(0, 0)]], (1, 1, 8))
            lp = np.log(np.full(a.shape, 0.3))
            r = gt.topo_reward(name, a, f[None], lp, lam=0.5)
            inv = a == gt.IDX[(0, 0)]
            if inv.any() and (~inv).any():
                assert r[inv].max() < r[~inv].min(), name
