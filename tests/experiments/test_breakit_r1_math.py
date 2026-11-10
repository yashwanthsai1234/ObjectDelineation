"""BREAK-IT round 1, mathematical correctness: edge cases of eval.signal, eval.ceiling, gate_power and
reward.advantages against Implementation.md (A1, A2, B1.9) and notes/research/math-stats.md (sections 1 and 11).

Each test pins one proven finding, now fixed. Fixed seeds, no data folder."""

import numpy as np
import pytest

from experiments.phase2 import gate_power as GP
from objdel.verify import ceiling as C
from objdel.verify import stats as signal
from objdel.verify.reward.advantages import advantages_v4
from objdel.verify.reward.cache import ExpertCache


def _group_terms(iou_last, f5_last):
    """One GRPO group of 8: seven legal rollouts, the last one parsed but illegal with the given raw terms."""
    return {"format_ok": [1.0] * 8, "legal": [1.0] * 7 + [0.0],
            "iou_tc": [0.4, 0.5, 0.6, 0.7, 0.5, 0.6, 0.4, iou_last],
            "f5": [0.3, 0.5, 0.7, 0.6, 0.4, 0.5, 0.6, f5_last],
            "g_area": [0.6, 0.7, 0.5, 0.4, 0.3, 0.6, 0.7, 0.0],
            "topo_any": [1.0] * 8, "topo_freq": [1.0] * 8, "sig_id": [7.0] * 8}


def test_illegal_rollout_nan_term_does_not_poison_the_group():
    a = advantages_v4(_group_terms(np.nan, np.inf), 8)
    assert np.isfinite(a).all(), a


def test_ci_with_a_nan_item_is_not_a_finite_interval_around_nan():
    x = np.array([0.5] * 9 + [np.nan])
    m, lo, hi = signal.ci(x, signal.draws([f"i{j}" for j in range(10)], boot=200, seed=0))
    assert not (np.isnan(m) and np.isfinite(lo) and np.isfinite(hi)), (m, lo, hi)


def test_b1_rules_report_the_empty_subset_draw_count():
    far, near = np.array([True, False, False, False]), np.array([False, True, False, False])
    # 5 patient draws over items of patients (0, 1, 1, 1): only the last draw holds the far item
    W = np.array([[0, 1, 1, 1]] * 4 + [[1, 1, 1, 1]], float)
    b1 = np.array([[0.5, 0.6, 0.6, 0.6]])
    base = np.array([[0.05, 0.6, 0.6, 0.6]])
    out = GP.b1_rules(b1, base, far, near, W)
    assert any("nan" in k for k in out), sorted(out)


def test_b3_rules_count_a_one_patient_both_found_run_as_a_fail():
    groups = np.array([0, 0, 1, 1])
    x = np.array([[0.5, 0.6, 0.0, 0.0]])  # both arms find only patient 0's items
    out = GP.b3_rules(x, x.copy(), groups)
    assert not out["both_found_iou"][0]


def test_cluster_t_empty_input_gives_its_own_error():
    with pytest.raises(ValueError, match="cluster_t"):
        signal.cluster_t(np.array([]), np.array([]))


def _disc(r, cy=16, cx=16, n=32):
    y, x = np.mgrid[:n, :n]
    return (y - cy) ** 2 + (x - cx) ** 2 <= r * r


def test_ceiling_item_with_answers_and_no_truncated_flags():
    cache = ExpertCache.build(np.stack([_disc(6), _disc(7), _disc(6, 17, 17)]))
    row = C.item(cache, ["not an outline", "still not"], None, 32, 32)
    assert row["model"] == 0.0
