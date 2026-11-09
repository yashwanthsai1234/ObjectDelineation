"""BREAK-IT round 2, mathematical correctness: edge cases of eval.ceiling, eval.signal, eval.metrics and
reward.advantages against Implementation.md and notes/research/math-stats.md (sections 1 and 11).

Each test pins one proven finding (xfail strict until fixed). Fixed seeds, no data folder, independent reference
code where a value is compared."""

import importlib

import numpy as np
import pytest
from scipy import stats

from objdel.eval import ceiling as C
from objdel.eval import signal
from objdel.reward.advantages import advantages_v4
from objdel.reward.cache import ExpertCache

M = importlib.import_module("objdel.eval.metrics")  # objdel.eval re-exports the function under the same name


def _disc(r, cy=16, cx=16, n=32):
    y, x = np.mgrid[:n, :n]
    return (y - cy) ** 2 + (x - cx) ** 2 <= r * r


@pytest.mark.xfail(strict=True, reason="BR2-math-1: ceiling.item gates illegal answers by multiplication "
                   "(iou_tc * legal); a parsed illegal answer with a NaN iou_tc (kept for logging) makes model NaN")
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


@pytest.mark.xfail(strict=True, reason="BR2-math-2: r_topo's pi_loo counts the signature of an illegal rollout, so "
                   "the legal rollouts' advantages depend on what an illegal rollout drew (eval counts it as FAIL)")
def test_illegal_rollout_signature_does_not_move_legal_advantages():
    a, b = advantages_v4(_topo_group(7.0), 8), advantages_v4(_topo_group(9.0), 8)
    np.testing.assert_allclose(a[:7], b[:7])


@pytest.mark.xfail(strict=True, reason="BR2-math-3: signal._ci drops NaN bootstrap replicates (spearman on a "
                   "constant draw) with nanquantile when the point is finite: no NaN-draw count, and the interval "
                   "can exclude its own point")
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
    assert lo <= rho <= hi, out


@pytest.mark.xfail(strict=True, reason="BR2-math-4: pass_at_k returns 1.0 when k > n even with 0 passes "
                   "(the n - c < k shortcut); with no pass pass@k is 0")
def test_pass_at_k_zero_passes_is_zero():
    assert M.pass_at_k(4, 0, 8) == 0.0
