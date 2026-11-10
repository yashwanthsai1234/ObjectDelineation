"""B3.5b: LocGuard stops only on a rising location term with a falling IoU_tc, with a false-stop rate < 1 % on
flat, autocorrelated per-step means (Implementation.md B3.5b)."""
import json
import pathlib

import numpy as np
import pytest
from transformers import TrainerControl

from objdel.policy.callbacks import LocGuard, slope_nw

ROOT = pathlib.Path(__file__).resolve().parents[2]
SD = 0.038  # per-step sd of rewards/iou_tc around its trend, full RL run (runs/rl/train_rl/rl.json)


def ar1(rng, n, phi, sd):
    e = rng.normal(0, sd * np.sqrt(1 - phi * phi), n)
    y = np.empty(n)
    y[0] = rng.normal(0, sd)
    for t in range(1, n):
        y[t] = phi * y[t - 1] + e[t]
    return y


def stops(loc, iou) -> int | None:
    g = LocGuard()
    for t, (a, b) in enumerate(zip(loc, iou, strict=True)):
        if g.add(a, b):
            return t + 1
    return None


def false_stop_rate(phi, runs=400, n=300, seed=0):
    rng = np.random.default_rng(seed)
    return np.mean([stops(0.2 + ar1(rng, n, phi, SD), 0.4 + ar1(rng, n, phi, SD)) is not None for _ in range(runs)])


def test_slope_nw_recovers_a_line():
    b, se = slope_nw(3.0 + 0.5 * np.arange(50) + 1e-3 * np.sin(np.arange(50)))
    assert b == pytest.approx(0.5, abs=1e-4) and se < 1e-3


def test_flat_noise_rarely_stops():
    assert false_stop_rate(0.0) < 0.01


def test_rising_location_with_falling_iou_stops_within_100_steps():
    rng = np.random.default_rng(1)
    t = np.arange(300)
    at = stops(0.1 + 0.002 * t + rng.normal(0, SD, 300), 0.5 - 0.001 * t + rng.normal(0, SD, 300))
    assert at is not None and at <= 100


def test_rising_both_does_not_stop():
    rng = np.random.default_rng(2)
    t = np.arange(300)
    assert stops(0.1 + 0.002 * t + rng.normal(0, SD, 300), 0.4 + 0.001 * t + rng.normal(0, SD, 300)) is None


def test_fewer_than_100_steps_never_stop():
    t = np.arange(99)
    assert stops(0.1 + 0.01 * t, 0.5 - 0.01 * t) is None


def test_autocorrelated_flat_runs_rarely_stop():
    """AR(1) coefficient of the full run's per-step rewards/iou_tc (detrended), when the record is present."""
    p = ROOT / "runs/rl/train_rl/rl.json"
    if not p.is_file():
        pytest.skip("local-only record runs/rl/train_rl/rl.json is absent")
    v = np.array([h["rewards/iou_tc/mean"] for h in json.loads(p.read_text())["log"] if "rewards/iou_tc/mean" in h])
    x = np.arange(v.size)
    e = v - np.polyval(np.polyfit(x, v, 1), x)
    phi = float(np.corrcoef(e[:-1], e[1:])[0, 1])
    assert 0 < phi < 1
    assert false_stop_rate(phi) < 0.01


def test_on_log_reads_the_off_centre_kinds_and_stops_training():
    g, c = LocGuard(window=2, lag=0), TrainerControl()
    g.block = lambda loc, iou: True  # every block passes: trips at the 4th logged step
    for _ in range(3):
        g.on_log(None, None, c, logs={"kind/shifted/loc": 0.1, "kind/shifted/iou_tc": 0.2, "kind/rescaled/loc": 0.3,
                                      "kind/rescaled/iou_tc": 0.4})
        g.on_log(None, None, c, logs={"eval_kind/shifted/loc": 0.1})  # eval logs are ignored
    assert not c.should_training_stop and g.loc == [0.2] * 3 and g.iou == [0.30000000000000004] * 3
    g.on_log(None, None, c, logs={"kind/shifted/loc": 0.1, "kind/shifted/iou_tc": 0.2})
    assert c.should_training_stop and c.should_save
