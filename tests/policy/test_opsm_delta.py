"""U68: the live OPSM delta callback (B2.4 rule 4)."""
import json
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("torch")  # the [dev]-only CI job has no ML stack

import torch

from objdel.policy.callbacks import OpsmDelta

ROWS = 256


def run(tmp_path, gaps, advs=None, start=0.7, state0=0):
    """Feeds one (kl, adv) block per step through the callback; returns (trainer, callback, deltas per step)."""
    t = SimpleNamespace(off_policy_mask_threshold=start)
    cb = OpsmDelta(t, str(tmp_path))
    cb.on_train_begin(None, SimpleNamespace(global_step=state0), None)
    seen = []
    for i, g in enumerate(gaps, state0 + 1):
        a = torch.full((ROWS,), -1.0) if advs is None else advs[i - state0 - 1]
        cb.add(torch.as_tensor(g, dtype=torch.float32), a)
        cb.on_step_end(None, SimpleNamespace(global_step=i), None)
        seen.append(t.off_policy_mask_threshold)
    return t, cb, seen


def test_delta_is_the_planted_start_until_step_20_then_q99_5_of_all_20_x_256_values_and_frozen(tmp_path):
    g = np.random.default_rng(0).exponential(0.01, (25, ROWS))
    t, _, seen = run(tmp_path, g)
    want = float(np.quantile(g[:20].ravel(), 0.995))
    assert seen[:19] == [0.7] * 19
    assert np.isclose(seen[19], want, rtol=1e-5) and seen[19:] == [seen[19]] * 6
    assert np.isclose(json.loads((tmp_path / "opsm_delta.json").read_text())["delta"], want, rtol=1e-5)


def test_resume_after_step_20_reads_the_file_and_changes_nothing(tmp_path):
    (tmp_path / "opsm_delta.json").write_text(json.dumps({"delta": 0.123, "guarded": False}))
    _, _, seen = run(tmp_path, np.full((30, ROWS), 5.0), state0=40)
    assert seen == [0.123] * 30


def test_guard_recomputes_once_when_the_masked_share_stays_above_2_percent_for_20_steps(tmp_path):
    g = np.full((60, ROWS), 0.01)
    g[20:] = 0.05  # gap grows after step 20: delta (0.01) masks every negative-advantage sequence
    t, cb, seen = run(tmp_path, g)
    assert seen[19] == np.float32(0.01) and seen[38] == seen[19]
    assert np.isclose(seen[39], 0.05) and seen[39:] == [seen[39]] * 21  # step 40 = 20th masked step, then frozen
    assert json.loads((tmp_path / "opsm_delta.json").read_text())["guarded"] is True


def test_no_guard_when_the_masked_share_dips_below_2_percent(tmp_path):
    g = np.full((60, ROWS), 0.01)
    g[20:] = 0.05
    neg, pos = torch.full((ROWS,), -1.0), torch.ones(ROWS)
    adv = [neg] * 30 + [pos] + [neg] * 19 + [pos] + [neg] * 9  # steps 31 and 51 mask nothing: no 20-step streak
    t, _, seen = run(tmp_path, g, adv)
    assert seen[-1] == seen[19]


def test_set_is_ratio_feeds_the_watch_in_train_mode_only():
    from tests.experiments.test_nvtx import fake_trainer

    t = fake_trainer(vllm_importance_sampling_mode="token_truncate", vllm_importance_sampling_clip_min=None,
                     vllm_importance_sampling_clip_max=2.0)  # fmt: skip
    t.opsm_watch = SimpleNamespace(cur=[])
    t.opsm_watch.add = lambda kl, adv: t.opsm_watch.cur.append((kl, adv))
    inputs = {"sampling_per_token_logps": torch.zeros(2, 3), "completion_mask": torch.ones(2, 3),
              "advantages": torch.tensor([1.0, -1.0])}
    t.set_is_ratio(torch.zeros(2, 3), inputs)
    assert len(t.opsm_watch.cur) == 1 and t.opsm_watch.cur[0][1].tolist() == [1.0, -1.0]
    t.model.training = False
    t.set_is_ratio(torch.zeros(2, 3), inputs)
    assert len(t.opsm_watch.cur) == 1
