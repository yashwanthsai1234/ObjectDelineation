"""B3.5: TRLReward with reward v5 (spec=REWARD_V5): advantage_v5 aggregator, gated loc logger, kind/<k>/loc and
completions/truncated_grad_share logs; grpo_config follows the spec's group. The default stays reward v4."""
import numpy as np
from test_train_rewards import CACHES, EOS, NAMES, PAD, W, batch, call_all

from objdel.rl import reward as R
from objdel.verify.reward import advantages as reward_advantages
from objdel.verify.reward.advantages import advantages_v5, legal_gate, topology_reward
from objdel.verify.reward.score import score_batch
from objdel.verify.reward.spec import REWARD_V4, REWARD_V5


def reward_v5():
    return R.TRLReward(CACHES.__getitem__, stop_ids={EOS, PAD}, spec=REWARD_V5)


def test_default_is_v4_and_v5_names_its_aggregator():
    assert R.TRLReward(CACHES.__getitem__, {EOS}).aggregator == R.AGGREGATOR == "advantage_v4"
    funcs, weights = reward_v5().funcs()
    assert [f.__name__ for f in funcs] == [*reward_advantages.LOGGED, "advantage_v5"] and weights[-1] == 1.0
    assert "loc" in reward_advantages.LOGGED and "loc" in reward_advantages.GATED and REWARD_V4.group == REWARD_V5.group


def test_v5_values_equal_score_batch_and_advantages_v5():
    kw, texts, keys = batch(truncate=(1,))
    got = call_all(reward_v5(), kw)
    ref = score_batch([(t, i == 1, NAMES, CACHES[k], W, W) for i, (t, k) in enumerate(zip(texts, keys, strict=True))])
    assert got["advantage_v5"] == advantages_v5(ref, 8, REWARD_V5).tolist()
    assert got["loc"] == (ref["loc"] * legal_gate(ref)).tolist()
    assert got["r_topo"] == topology_reward(ref, 8, REWARD_V5).tolist()
    a = np.array(got["advantage_v5"]).reshape(-1, 8)
    assert np.allclose(a.sum(1), 0, atol=1e-12) and np.ptp(a) > 0


def test_v5_logs_kind_loc_and_truncated_grad_share():
    kw, _, _ = batch(truncate=(1,))
    kw["completion_ids"][1] = [5] * 30  # a long truncated row
    logs = {}
    call_all(reward_v5(), kw, log_metric=logs.__setitem__, kind=["shifted"] * 8 + ["centred"] * 8)
    assert {"kind/shifted/loc", "kind/centred/loc", "kind/shifted/iou_tc"} <= set(logs)
    assert 0 <= logs["completions/truncated_grad_share"] <= 1


def test_grpo_config_follows_the_spec_group(tmp_path):
    from objdel.rl import config as configs

    c = configs.grpo_config("qwen3-vl-2b", str(tmp_path), reward=REWARD_V5, use_cpu=True, use_vllm=False)
    assert c.num_generations == REWARD_V5.group and len(c.reward_weights) == len(reward_advantages.LOGGED) + 1


def test_i26_tilt_stays_removed(tmp_path):
    """I26 (B3.8): schedule takes no tilt, TRLReward logs no cell/<ij>/group_sd, num_generations is the one group."""
    import inspect

    import pytest

    from objdel.data import windows
    from objdel.rl import config as configs

    assert not any("tilt" in p for p in inspect.signature(windows.schedule).parameters)
    kw, _, _ = batch()
    logs = {}
    call_all(reward_v5(), kw, log_metric=logs.__setitem__, cell=["11"] * 8 + ["23"] * 8)
    assert "cell/11/iou_tc" in logs and not any("group_sd" in k for k in logs)
    with pytest.raises(ValueError, match="exactly 8"):
        configs.grpo_config("qwen3-vl-2b", str(tmp_path), reward=REWARD_V5, num_generations=16, use_cpu=True,
                            use_vllm=False)
