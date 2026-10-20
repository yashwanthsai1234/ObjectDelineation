"""objdel.train.configs: plan D8 (GRPO), SFT settings, D12 image sizes, the I5 completion budget."""
import pytest

pytest.importorskip("trl")

from objdel.reward import REWARD_V4  # noqa: E402
from objdel.train import configs  # noqa: E402
from objdel.train.rewards import LOGGED  # noqa: E402


def test_models_are_the_screened_three():
    assert configs.MODELS == {"smolvlm-256m": "HuggingFaceTB/SmolVLM-256M-Instruct",
                              "lfm2.5-vl-450m": "LiquidAI/LFM2.5-VL-450M",
                              "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct"}


@pytest.mark.parametrize("model", list(configs.MODELS))
def test_grpo_config_is_plan_d8(model, tmp_path):
    c = configs.grpo_config(model, str(tmp_path))
    assert (c.loss_type, c.epsilon, c.epsilon_high, c.beta) == ("dr_grpo", 0.2, 0.28, 1e-3)
    assert c.mask_truncated_completions and c.lr_scheduler_type == "constant" and c.warmup_steps == 0
    assert c.num_generations == REWARD_V4.group and c.num_generations_eval is None and c.temperature == 1.0
    assert (c.scale_rewards, c.multi_objective_aggregation) == ("none", "sum_then_normalize")
    assert c.reward_weights == [0.0] * len(LOGGED) + [1.0]
    assert c.max_completion_length == configs.MAX_COMPLETION[model] == {"lfm2.5-vl-450m": 1024}.get(model, 1536)
    assert c.use_vllm == (model != "smolvlm-256m")  # TRL #6294: vLLM corrupts SmolVLM prompts
    assert not c.remove_unused_columns and c.log_completions
    # whole groups on every process, at any world size (tests/test_train_trl.py simulates 1, 2, 4 processes)
    assert c.per_device_train_batch_size % c.num_generations == 0
    assert c.per_device_eval_batch_size % c.num_generations == 0
    assert c.generation_batch_size == 256 and c.steps_per_generation == c.gradient_accumulation_steps


@pytest.mark.parametrize("model", list(configs.MODELS))
def test_sft_config(model, tmp_path):
    c = configs.sft_config(model, str(tmp_path))
    assert c.assistant_only_loss and c.max_length is None and not c.packing
    assert c.learning_rate == (5e-4 if configs.lora(model) else 5e-5)


def test_lora_only_where_the_plan_needs_it():
    assert configs.lora("smolvlm-256m") is None and configs.lora("lfm2.5-vl-450m") is None
    q = configs.lora("qwen3-vl-2b")
    assert (q.r, q.lora_alpha, q.target_modules) == (16, 32, "all-linear")


def test_image_sizes_are_plan_d12():
    assert configs.processor_kwargs("smolvlm-256m", "riga") == {"size": {"longest_edge": 1024}, "padding_side": "left"}
    assert configs.processor_kwargs("smolvlm-256m", "lidc") == {"size": {"longest_edge": 512}, "padding_side": "left"}
    for m in ("lfm2.5-vl-450m", "qwen3-vl-2b"):
        assert configs.processor_kwargs(m, "riga") == {"padding_side": "left"}


@pytest.mark.parametrize("supported", [True, False])
def test_bf16_only_where_the_hardware_supports_it(supported, tmp_path, monkeypatch):
    """CPU-only Linux (CI train job): no bf16, so TrainingArguments does not raise; CUDA Ampere+/MPS: bf16."""
    monkeypatch.setattr(configs, "is_torch_bf16_gpu_available", lambda: supported)
    # use_cpu: TrainingArguments re-checks the real hardware, and accepts bf16 (CPU autocast) only with use_cpu
    for c in (configs.grpo_config("smolvlm-256m", str(tmp_path), use_cpu=True),
              configs.sft_config("smolvlm-256m", str(tmp_path), use_cpu=True)):
        assert c.bf16 is supported and not c.fp16
    c = configs.grpo_config("smolvlm-256m", str(tmp_path), bf16=False, max_steps=3)  # overrides win
    assert (c.bf16, c.max_steps) == (False, 3)
    assert not configs.sft_config("smolvlm-256m", str(tmp_path), fp16=True, use_cpu=True).bf16


def test_configs_build_on_this_machine(tmp_path):
    """Without overrides (CPU-only Linux CI included): the real hardware check decides bf16."""
    for m in configs.MODELS:
        assert configs.grpo_config(m, str(tmp_path)).bf16 == configs.is_torch_bf16_gpu_available()
        assert configs.sft_config(m, str(tmp_path)).bf16 == configs.is_torch_bf16_gpu_available()
