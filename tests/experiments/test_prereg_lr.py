"""E1.1: the pre-registered GRPO learning rates and the code's one learning-rate function cannot drift."""
from experiments.phase2.prereg import G_PROBE
from objdel.train import configs


def test_rl_lr_equals_the_pre_registered_rate_for_every_model():
    for m in configs.MODELS:
        assert configs.rl_lr(m) == G_PROBE["lr_lora" if m in configs._LORA else "lr_full"], m


def test_grpo_config_defaults_to_rl_lr():
    for m in configs.MODELS:
        assert configs.grpo_config(m, "out").learning_rate == configs.rl_lr(m), m
