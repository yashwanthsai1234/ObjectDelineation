"""BREAK-IT round 2, mathematical correctness: edge cases of eval.ceiling, eval.signal, eval.metrics and
reward.advantages against Implementation.md and notes/research/math-stats.md (sections 1 and 11).

Each test pins one proven finding (xfail strict until fixed). Fixed seeds, no data folder, independent reference
code where a value is compared."""

import importlib

import numpy as np
import pytest
from scipy import stats

from objdel.verify import ceiling as C
from objdel.verify import stats as signal
from objdel.verify.reward.advantages import advantages_v4
from objdel.verify.reward.cache import ExpertCache

M = importlib.import_module("objdel.verify.metrics")  # objdel.eval re-exports the function under the same name


def _disc(r, cy=16, cx=16, n=32):
    y, x = np.mgrid[:n, :n]
    return (y - cy) ** 2 + (x - cx) ** 2 <= r * r


def test_ceiling_model_ignores_raw_terms_of_illegal_answers(monkeypatch):
    cache = ExpertCache.build(np.stack([_disc(6), _disc(7), _disc(6, 17, 17)]))

    def fake_score_many(jobs):  # every answer parsed but illegal; raw IoU_tc kept as NaN, as score.py may log it
        return [({"iou_tc": np.nan, "legal": 0.0}, {}) for _ in jobs]

    monkeypatch.setattr(C.score, "score_many", fake_score_many)
    with np.errstate(invalid="ignore"):
        row = C.item(cache, ["a", "b"], [False, False], 32, 32)
    assert row["model"] == 0.0 and row["model_all_k"] == 0.0, row


def _topo_group(sig_last):
    """Seven legal valid rollouts (answers 7 x4, 8 x3) and one parsed illegal rollout with signature sig_last."""
    return {"format_ok": [1.0] * 8, "legal": [1.0] * 7 + [0.0],
            "iou_tc": [0.6, 0.7, 0.8, 0.6, 0.7, 0.8, 0.9, 0.5], "f5": [0.5] * 8, "g_area": [0.6] * 8,
            "topo_any": [1.0] * 8, "topo_freq": [0.5] * 8, "sig_id": [7.0] * 4 + [8.0] * 3 + [sig_last]}


def test_illegal_rollout_signature_does_not_move_legal_advantages():
    a, b = advantages_v4(_topo_group(7.0), 8), advantages_v4(_topo_group(9.0), 8)
    np.testing.assert_allclose(a[:7], b[:7])


def test_spearman_ci_counts_or_refuses_nan_draws():
    x, y = np.array([0.0] * 9 + [1.0]), np.arange(10.0)
    D = signal.draws([str(i) for i in range(10)], 2000, 0)
    with np.errstate(invalid="ignore"), pytest.warns(stats.ConstantInputWarning):
        nan_draws = int(np.isnan([stats.spearmanr(x[d], y[d]).statistic for d in D]).sum())
        out = signal.spearman_ci(x, y)
    assert nan_draws > 0  # the fixture does hit the case (673 of 2,000 draws)
    rho, lo, hi = out[:3]
    # math-stats section 1: the NaN-draw count goes beside every CI; a percentile CI must contain its point
    assert len(out) > 3 or not (np.isfinite(lo) and np.isfinite(hi)), out
    assert not np.isfinite(lo) or lo <= rho <= hi, out  # a refused interval (NaN bounds) carries the count


def test_pass_at_k_zero_passes_is_zero():
    assert M.pass_at_k(4, 0, 8) == 0.0
