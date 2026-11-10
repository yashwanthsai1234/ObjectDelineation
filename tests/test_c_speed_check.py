"""c_speed_check --sync: the pass rules of G-SY1..G-SY3 (the box run itself needs the real engine)."""

import pytest

torch = pytest.importorskip("torch")
from experiments.phase2 import c_speed_check as C  # noqa: E402

MB = 1 << 20


def test_sync_passes_a_fast_flat_series():
    v = C.sync_verdict([0.05] * 20, [100 * MB] * 50, equal=[True] * 3)
    assert v["passed"] and v["median_s"] == 0.05


@pytest.mark.parametrize(
    "kw",
    [
        {"seconds": [0.3] * 20},  # G-SY2: median over 0.2 s
        {"mem": [100 * MB] * 49 + [102 * MB]},  # G-SY3: a leak of 2 MB
        {"equal": [True, False, True]},  # G-SY1: one push differs from TRL's
        {"seconds": [0.05] * 19 + [0.5]},  # G-SY3: seconds not stable
    ],
)
def test_sync_fails_each_rule(kw):
    args = {"seconds": [0.05] * 20, "mem": [100 * MB] * 50, "equal": [True] * 3} | kw
    assert not C.sync_verdict(**args)["passed"]


def test_hash_params_sees_one_changed_value():
    a = {"w": torch.ones(4, dtype=torch.bfloat16)}
    b = {"w": torch.ones(4, dtype=torch.bfloat16)}
    assert C.hash_params(a.items()) == C.hash_params(b.items())
    b["w"][2] = 2
    assert C.hash_params(a.items()) != C.hash_params(b.items())


def test_parse_takes_sync_and_the_train_rl_flags():
    sync, a = C.parse(["--sync", "--model", "qwen3vl2b", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", "1"])
    assert sync and a.fast_sync is False
