"""perf fixes F1-F3 (experiments/phase2/perf/fix_test.md): G-probe generates with vLLM colocate, the GRPO step keeps
256 answers in larger micro-batches, and the per-micro-batch host syncs are off."""

from pathlib import Path

import pytest

pytest.importorskip("trl")
from objdel.reward.spec import REWARD_V4  # noqa: E402
from objdel.train import configs  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KV_PER_TOKEN = 112 * 1024  # Qwen3-VL-2B bf16: 28 layers x 8 KV heads x 128 x (K, V) x 2 bytes (PLAN D13)
PROMPT = 330  # tokens: chat text + 256 image tokens (R1 collator)


def test_g_probe_generates_with_vllm_colocate():
    src = (ROOT / "experiments/phase2/g_probe.py").read_text()
    assert '"use_vllm": False' not in src  # HF generate of 256 sequences overflows the card (issues RL-1)
    c = configs.grpo_config("qwen3-vl-2b", "o", max_steps=20, save_strategy="no")
    assert c.use_vllm and c.vllm_mode == "colocate" and c.vllm_enable_sleep_mode
    assert 0 < c.vllm_gpu_memory_utilization <= 0.4


@pytest.mark.parametrize("card_gib", [24, 32])
def test_vllm_kv_budget_holds_one_runaway_group(card_gib):
    """vLLM pages the KV cache and schedules what fits; its share of the card, after the weights, must hold at least
    one whole group of G answers that all run to MAX_COMPLETION."""
    c = configs.grpo_config("qwen3-vl-2b", "o")
    weights = 2.13e9 * 2  # bf16
    kv = c.vllm_gpu_memory_utilization * card_gib * 2**30 - weights
    group = REWARD_V4.group * (PROMPT + configs.MAX_COMPLETION["qwen3-vl-2b"]) * KV_PER_TOKEN
    assert kv >= group


@pytest.mark.parametrize("model", list(configs.MODELS))
def test_grpo_step_keeps_256_answers_in_larger_micro_batches(model):
    c = configs.grpo_config(model, "o")
    g = REWARD_V4.group
    assert c.per_device_train_batch_size * c.gradient_accumulation_steps == 256 == c.generation_batch_size
    assert c.per_device_train_batch_size % g == 0 and c.per_device_train_batch_size > g


@pytest.mark.parametrize("model", list(configs.MODELS))
def test_configs_turn_off_the_per_micro_batch_syncs(model, tmp_path):
    for c in (configs.sft_config(model, str(tmp_path)), configs.grpo_config(model, str(tmp_path))):
        assert c.logging_nan_inf_filter is False and c.accelerator_config.non_blocking is True


def test_trainers_skip_flop_counting():
    for f in ("experiments/phase2/train_sft.py", "experiments/phase2/g_probe.py"):
        src = (ROOT / f).read_text()
        assert "floating_point_ops = lambda inputs: 0" in src or (
            "def floating_point_ops(self, inputs):" in src and "            return 0" in src), f
