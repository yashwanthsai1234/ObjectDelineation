"""E1.1: the pre-registered GRPO learning rates and the code's one learning-rate function cannot drift."""
from experiments.phase2.prereg import G_PROBE
from objdel.train import configs


def test_rl_lr_equals_the_pre_registered_rate_for_every_model():
    for m in configs.MODELS:
        assert configs.rl_lr(m) == G_PROBE["lr_lora" if m in configs._LORA else "lr_full"], m


def test_full_fine_tuning_rate_stays_a_frozen_record():
    """CT-62: CLEAN8 retires full fine-tuning; the pre-registered full-FT rate is a record and does not change."""
    assert G_PROBE["lr_full"] == 5e-6 and configs._LORA == frozenset(configs.MODELS)


def test_grpo_config_defaults_to_rl_lr():
    for m in configs.MODELS:
        assert configs.grpo_config(m, "out").learning_rate == configs.rl_lr(m), m
