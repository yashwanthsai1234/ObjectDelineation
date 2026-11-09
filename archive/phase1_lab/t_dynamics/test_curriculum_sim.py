from dataclasses import replace

import numpy as np

import curriculum_sim as cs


def test_group_signal_rewards_fixing_the_error_type():
    rng = np.random.default_rng(0)
    p = np.tile([[0.0, 0.0, 0.0, 0.0, 0.5]], (64, 1))
    e = cs.rollouts(p, 8, rng)
    c, _ = cs.group_signal(e, np.full(5, 0.2))
    assert c[:, 4].mean() > 0.3 and np.allclose(c[:, :4], 0)


def test_group_signal_zero_for_constant_groups():
    e = np.ones((4, 8, 5))
    c, _ = cs.group_signal(e, np.full(5, 0.2))
    assert np.allclose(c, 0)


def test_advantages_inside_signal_are_zero_mean():
    rng = np.random.default_rng(1)
    e = cs.rollouts(rng.random((16, 5)), 8, rng)
    r = 1 - (e * 0.2).sum(-1)
    sd = r.std(-1, keepdims=True)
    A = (r - r.mean(-1, keepdims=True)) / np.where(sd > 0, sd, 1)
    assert np.allclose(A.sum(-1), 0)


def test_run_is_deterministic_and_learns():
    w = replace(cs.World(), T=200)
    a = cs.run(w, "uniform", 3)
    b = cs.run(w, "uniform", 3)
    assert np.array_equal(a["test_fail"], b["test_fail"])
    assert a["test_fail"][-1] < 0.7 * a["test_fail"][0]


def test_same_seed_same_world_across_samplers():
    w = replace(cs.World(), T=0)
    a = cs.run(w, "uniform", 5)
    b = cs.run(w, "typed-balance", 5)
    assert np.isclose(a["test_fail"][0], b["test_fail"][0])


def test_memorisation_does_not_transfer():
    w = replace(cs.World(), T=300, alpha_m=0.5, alpha=(0,) * 5)
    r = cs.run(w, "uniform", 0)
    assert r["train_fail"][-1] < r["train_fail"][0] - 0.1
    assert np.isclose(r["test_fail"][-1], r["test_fail"][0])


def test_samplers_return_valid_batches():
    rng = np.random.default_rng(2)
    w = cs.World()
    init = rng.random((w.n_train, cs.K)) * 0.5
    for name, mk in cs.SAMPLERS.items():
        s = mk(w.n_train, init, rng, w)
        idx = s.draw(w.B)
        assert len(set(idx.tolist())) == w.B and idx.max() < w.n_train


def test_alarm_fires_on_rising_val_error_and_not_on_flat():
    a = cs.Alarm()
    fired = False
    for t in range(30):
        fired |= a.update(t, 0.5 + 0.01 * t, np.array([0.1, 0.1, 0.1 + 0.01 * t, 0.1, 0.1]), 100)
    assert fired
    rng = np.random.default_rng(3)
    fires = 0
    for _ in range(200):  # flat val errors (pure noise), 60 evals each
        b = cs.Alarm()
        fires += any(b.update(t, 0.5 + 0.01 * t, 0.1 + 0.03 * rng.standard_normal(5), 100) for t in range(60))
    assert fires / 200 < 0.05


def test_hack_raises_train_reward_but_hurts_heldout():
    w = replace(cs.world_for("S5 + reward hack"), T=600)
    r = cs.run(w, "uniform", 0)
    assert r["H"][-1] > 1
    assert r["test_fail"][-1] > r["test_fail"].min() + 0.05
