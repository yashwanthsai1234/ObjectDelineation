"""Task 4: group-advantage schemes."""
import numpy as np
import pytest

from rewardcore.advantages import SCHEMES

NAMES = ("iou", "band", "topo", "legal")
W = np.array([0.3, 0.3, 0.3, 0.1])
FLOORS = np.array([0.05, 0.05, 0.25, 0.25])
LOW_VAR = np.array([False, False, True, True])
KW = dict(floors=FLOORS, low_var=LOW_VAR, names=NAMES)


def batch(rng, b=16):
    """Groups of 8 with terms that agree on the order: output j is better than output j+1."""
    q = np.linspace(0.95, 0.4, 8)
    R = np.empty((b, 8, 4))
    R[..., 0] = q + rng.normal(0, 0.005, (b, 8))
    R[..., 1] = np.clip(q + 0.05 + rng.normal(0, 0.005, (b, 8)), 0, 1)
    R[..., 2] = (np.arange(8) < 6).astype(float)
    R[..., 3] = 1.0
    return R


@pytest.mark.parametrize("name", list(SCHEMES))
def test_agreeing_terms_keep_the_order(name):
    R = batch(np.random.default_rng(0))
    a, _ = SCHEMES[name](R, W, **KW)
    assert np.all(np.diff(a, axis=1) < 0)


@pytest.mark.parametrize("name", [n for n in SCHEMES if n != "gated"])
def test_contributions_sum_to_advantage_and_constant_term_contributes_nothing(name):
    R = batch(np.random.default_rng(1))
    a, c = SCHEMES[name](R, W, **KW)
    assert np.allclose(c.sum(-1), a)
    if name not in ("gdpo", "gdpo_batch", "gdpo_floor", "gdpo_mean", "gdpo_floor_gated"):
        assert np.allclose(c[..., 3], 0)  # legality exactly constant
    else:
        # final rescaling spreads the centring evenly; the term's own share is zero before it
        assert np.ptp(c[..., 3], axis=1).max() < 1e-12


def test_dark_room_gdpo_amplifies_noise_floors_do_not():
    rng = np.random.default_rng(2)
    R = batch(rng)
    Rn = R.copy()
    Rn[..., 3] = 1 - rng.uniform(0, 1e-4, Rn.shape[:2])  # numeric noise on a constant term
    change = {}
    for name in SCHEMES:
        a0, _ = SCHEMES[name](R, W, **KW)
        a1, _ = SCHEMES[name](Rn, W, **KW)
        change[name] = np.abs(a1 - a0).max()
    assert change["gdpo"] > 0.1  # 1e-4 of noise becomes a sizeable advantage
    for name in ("grpo", "drgrpo", "gdpo_floor", "gdpo_mean", "gated", "gdpo_floor_gated"):
        assert change[name] < 1e-2, (name, change[name])


def test_gates():
    R = batch(np.random.default_rng(3))
    R[:, 0, 3] = 0.0  # best output is illegal
    a, _ = SCHEMES["gdpo_floor_gated"](R, W, **KW)
    assert np.all(a[:, 0] == a.min(1))
    R = batch(np.random.default_rng(4))
    R[:, :, 0] = 0.3  # IoU below the gate: topology must not matter
    R2 = R.copy()
    R2[..., 2] = 1 - R2[..., 2]
    a1, _ = SCHEMES["gdpo_floor_gated"](R, W, **KW)
    a2, _ = SCHEMES["gdpo_floor_gated"](R2, W, **KW)
    assert np.allclose(a1, a2)
