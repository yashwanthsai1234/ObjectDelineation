"""E1.1: the pre-registered GRPO learning rates and the code's one learning-rate function cannot drift."""
from experiments.phase2.prereg import G_PROBE
from objdel.policy import models as policy_models
from objdel.rl import config as rl_config


def test_rl_lr_equals_the_pre_registered_rate_for_every_model():
    for m in policy_models.MODELS:
        assert rl_config.rl_lr(m) == G_PROBE["lr_lora" if m in policy_models._LORA else "lr_full"], m


def test_full_fine_tuning_rate_stays_a_frozen_record():
    """CT-62: CLEAN8 retires full fine-tuning; the pre-registered full-FT rate is a record and does not change."""
    assert G_PROBE["lr_full"] == 5e-6 and policy_models._LORA == frozenset(policy_models.MODELS)


def test_grpo_config_defaults_to_rl_lr():
    for m in policy_models.MODELS:
        assert rl_config.grpo_config(m, "out").learning_rate == rl_config.rl_lr(m), m


def test_rl_lr_refuses_an_unknown_model_and_an_override_wins():
    """U13: no silent 5e-6 for a typo; grpo_config's explicit learning_rate wins over rl_lr."""
    import pytest

    with pytest.raises(ValueError, match="qwen3-vl-2"):
        rl_config.rl_lr("qwen3-vl-2")
    assert rl_config.grpo_config("qwen3-vl-2b", "out", learning_rate=1e-3).learning_rate == 1e-3
