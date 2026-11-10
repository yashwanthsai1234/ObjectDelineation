"""perf fixes F1-F3 (experiments/phase2/perf/fix_test.md): G-probe generates with vLLM colocate, the GRPO step keeps
256 answers in larger micro-batches, and the per-micro-batch host syncs are off."""

from pathlib import Path

import pytest

pytest.importorskip("trl")
from objdel.policy import models as policy_models  # noqa: E402
from objdel.rl import config as rl_config  # noqa: E402
from objdel.sft import config as sft_config  # noqa: E402
from objdel.verify.reward.spec import REWARD_V4  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KV_PER_TOKEN = 112 * 1024  # Qwen3-VL-2B bf16: 28 layers x 8 KV heads x 128 x (K, V) x 2 bytes (PLAN D13)
PROMPT = 330  # tokens: chat text + 256 image tokens (R1 collator)


def test_g_probe_generates_with_vllm_colocate():
    src = (ROOT / "src/objdel/rl/probe.py").read_text()
    assert '"use_vllm": False' not in src  # HF generate of 256 sequences overflows the card (issues RL-1)
    c = rl_config.grpo_config("qwen3-vl-2b", "o", max_steps=20, save_strategy="no")
    assert c.use_vllm and c.vllm_mode == "colocate" and c.vllm_enable_sleep_mode
    assert 0 < c.vllm_gpu_memory_utilization <= 0.4


@pytest.mark.parametrize("card_gib", [24, 32])
def test_vllm_kv_budget_holds_one_runaway_group(card_gib):
    """vLLM pages the KV cache and schedules what fits; its share of the card, after the weights, must hold at least
    one whole group of G answers that all run to MAX_COMPLETION."""
    c = rl_config.grpo_config("qwen3-vl-2b", "o")
    weights = 2.13e9 * 2  # bf16
    kv = c.vllm_gpu_memory_utilization * card_gib * 2**30 - weights
    group = REWARD_V4.group * (PROMPT + policy_models.MAX_COMPLETION["qwen3-vl-2b"]) * KV_PER_TOKEN
    assert kv >= group


@pytest.mark.parametrize("model", list(policy_models.MODELS))
def test_grpo_step_keeps_256_answers_in_larger_micro_batches(model):
    c = rl_config.grpo_config(model, "o")
    g = REWARD_V4.group
    assert c.per_device_train_batch_size * c.gradient_accumulation_steps == 256 == c.generation_batch_size
    assert c.per_device_train_batch_size % g == 0 and c.per_device_train_batch_size > g


@pytest.mark.parametrize("model", list(policy_models.MODELS))
def test_configs_turn_off_the_per_micro_batch_syncs(model, tmp_path):
    for c in (sft_config.sft_config(model, str(tmp_path)), rl_config.grpo_config(model, str(tmp_path))):
        assert c.logging_nan_inf_filter is False and c.accelerator_config.non_blocking is True


def test_trainers_skip_flop_counting():
    """SFT: sft.train's Trainer subclass; GRPO: rl.build.build_trainer builds trl_adapter.FastGRPO
    (tests/rl/test_trl_adapter_equivalence.py), for the G-probe and for the full run (rl.run)."""
    src = (ROOT / "src/objdel/sft/train.py").read_text()
    assert "def floating_point_ops(self, inputs):" in src and "        return 0" in src
    build = (ROOT / "src/objdel/rl/build.py").read_text()
    assert "cls = GRPOTrainer if flags.plain else rl_trl_adapter.FastGRPO" in build
    assert "build.build_trainer(a.model, Path(a.ckpt), train" in (ROOT / "src/objdel/rl/run.py").read_text()


def test_vllm_context_is_bounded():
    """TRL's default vllm_max_model_length is the model's (262,144 for Qwen3-VL): vLLM sizes its first KV block and
    graph capture for it and fails at start. Our longest sequence: ~400 prompt tokens (a 512 px image is 256) + 1536."""
    c = rl_config.grpo_config("qwen3-vl-2b", "o")
    assert policy_models.MAX_COMPLETION["qwen3-vl-2b"] + 512 <= c.vllm_max_model_length <= 4096


def test_g_probe_takes_the_profile_arm_flags(capsys):
    from experiments.phase2 import g_probe

    with pytest.raises(SystemExit) as e:
        g_probe.main(["--help"])
    out = capsys.readouterr().out
    assert e.value.code == 0 and all(f in out for f in ("--plain", "--trace", "--compile"))


def test_grpo_config_keeps_the_p1_measured_choice():
    """P1 on an RTX 5090 (2026-11-06, experiments/phase2/perf/p1_results.md): micro-batch 16 x 16 at 44.3 s/step beat
    32 x 8 (59.1); no checkpointing ran out of 31.45 GiB; compiled layers (47.9) and CUDA Graphs (crash) stay off."""
    from objdel.rl import config as rl_config

    c = rl_config.grpo_config("qwen3-vl-2b", "/tmp/o")
    assert (c.per_device_train_batch_size, c.gradient_accumulation_steps) == (16, 16)
    assert c.gradient_checkpointing is True and c.torch_compile is False
