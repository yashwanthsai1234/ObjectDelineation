"""objdel.train.configs: plan D8 (GRPO), SFT settings, D12 image sizes, the I5 completion budget; CLEAN7 (the 256M
model retired: CT-05, CT-06, CT-07)."""
import re
import subprocess
from dataclasses import asdict
from pathlib import Path

import pytest

pytest.importorskip("trl")

from objdel.reward import REWARD_V4  # noqa: E402
from objdel.train import configs  # noqa: E402
from objdel.train.rewards import LOGGED  # noqa: E402

OLD = "smol" + "vlm-256m"  # the retired key, built so this file never holds the word


def test_models_and_budgets():
    """CT-05: no exact set (a 4B row may join); every model has a budget, the retired key is gone."""
    assert {"qwen3-vl-2b", "lfm2.5-vl-450m"} <= set(configs.MODELS)
    assert configs.MODELS["qwen3-vl-2b"] == "Qwen/Qwen3-VL-2B-Instruct"
    assert OLD not in configs.MODELS and OLD not in configs.MAX_COMPLETION
    assert set(configs.MODELS) == set(configs.MAX_COMPLETION)


@pytest.mark.parametrize("model", list(configs.MODELS))
def test_grpo_config_is_plan_d8(model, tmp_path):
    c = configs.grpo_config(model, str(tmp_path))
    assert (c.loss_type, c.epsilon, c.epsilon_high, c.beta) == ("dr_grpo", 0.2, 0.28, 1e-3)
    assert c.mask_truncated_completions and c.lr_scheduler_type == "constant" and c.warmup_steps == 0
    assert c.num_generations == REWARD_V4.group and c.num_generations_eval is None and c.temperature == 1.0
    assert (c.scale_rewards, c.multi_objective_aggregation) == ("none", "sum_then_normalize")
    assert c.reward_weights == [0.0] * len(LOGGED) + [1.0]
    budget = {"lfm2.5-vl-450m": 1024}.get(model, 1536)
    assert c.max_completion_length == configs.MAX_COMPLETION[model] == budget
    assert c.use_vllm is True
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
    assert configs.lora("lfm2.5-vl-450m") is None
    q = configs.lora("qwen3-vl-2b")
    # Phase 2 (plan D6 / D7, M2): language model only, vision tower and projector excluded
    assert (q.r, q.lora_alpha, q.target_modules, q.exclude_modules) == (16, 32, "all-linear", configs.VISION)


def test_processors_keep_their_default_image_sizes():
    """CT-05 (D12 after CLEAN7): left padding only, for every model and any dataset string (no size key, no
    lookup on the dataset)."""
    for m in configs.MODELS:
        for d in ("lidc", "riga", "anything"):
            assert configs.processor_kwargs(m, d) == {"padding_side": "left"}
    with pytest.raises(KeyError):
        configs.sampling(OLD)  # the retired key has no budget


GRPO_KEYS = ("loss_type", "epsilon", "epsilon_high", "beta", "num_generations", "per_device_train_batch_size",
             "max_completion_length", "vllm_gpu_memory_utilization", "vllm_enable_sleep_mode", "vllm_max_model_length",
             "learning_rate", "reward_weights")
GRPO_PRE_CLEAN = {  # captured at 08e6455 (asdict of grpo_config(m, out, use_cpu=True, use_vllm=False))
    m: {"loss_type": "dr_grpo", "epsilon": 0.2, "epsilon_high": 0.28, "beta": 0.001, "num_generations": 8,
        "per_device_train_batch_size": 16, "max_completion_length": n, "vllm_gpu_memory_utilization": 0.35,
        "vllm_enable_sleep_mode": True, "vllm_max_model_length": 4096, "learning_rate": lr,
        "reward_weights": [0.0] * 7 + [1.0]}
    for m, n, lr in (("qwen3-vl-2b", 1536, 1e-05), ("lfm2.5-vl-450m", 1024, 1e-06))}


@pytest.mark.parametrize("model", sorted(GRPO_PRE_CLEAN))
def test_grpo_config_values_and_the_use_vllm_override(model, tmp_path):
    """CT-06: vLLM for every model unless the caller overrides it; the plan values are the pre-CLEAN ones."""
    c = configs.grpo_config(model, str(tmp_path), use_cpu=True, use_vllm=False)
    assert c.use_vllm is False
    assert {k: v for k, v in asdict(c).items() if k in GRPO_KEYS} == GRPO_PRE_CLEAN[model]
    assert configs.grpo_config(model, str(tmp_path), use_vllm=True).use_vllm is True
    assert configs.grpo_config(model, str(tmp_path), use_cpu=True, use_vllm=None).use_vllm is None  # pre-CLEAN
    assert configs.grpo_config(model, str(tmp_path), use_cpu=True, use_vllm=0).use_vllm == 0  # pre-CLEAN
    with pytest.raises(ValueError):
        configs.grpo_config(model, str(tmp_path), num_generations=16)


SMOL = re.compile("smol" + "vlm|hugging" + "facetb", re.I)
LFM = re.compile(r"(?<![a-z])" + "l" + "fm", re.I)
TEXT = (".py", ".yaml", ".md", ".toml")


def test_forbidden_word_needles():
    """CT-07 checker part (the lfm needle is for CLEAN8's Wave 3 scan)."""
    for t in ("selfmade", "shelfmark", "aGVsZm1Gx"):
        assert not LFM.search(t)
    for t in ("LFM2.5-VL", "lfm2.5-vl-450m", "# uses LFM"):
        assert LFM.search(t)
    for t in ("Smol" + "VLM", "smol" + "vlm_256m", "SMOL" + "VLM", "Hugging" + "FaceTB/x"):
        assert SMOL.search(t)


def scanned() -> list[str]:
    """Text sources of src/, experiments/phase2/*.py and tests/ (git ls-files -co, so no __pycache__)."""
    root = Path(__file__).resolve().parents[1]
    files = subprocess.run(["git", "-C", str(root), "ls-files", "-z", "-co", "--exclude-standard", "src", "tests",
                            "experiments/phase2"], capture_output=True, text=True, check=True).stdout.split("\0")
    keep = [f for f in files if f.endswith(TEXT) and (root / f).is_file()
            and (f.startswith(("src/", "tests/")) or re.fullmatch(r"experiments/phase2/[^/]+\.py", f))]
    return [f for f in keep if SMOL.search((root / f).read_text(errors="replace"))]


def test_no_retired_model_text_in_code():
    """CT-07 repo part, Wave 1 (CLEAN-G10): no file in the scope names the retired 256M model."""
    assert scanned() == []


@pytest.mark.parametrize("supported", [True, False])
def test_bf16_only_where_the_hardware_supports_it(supported, tmp_path, monkeypatch):
    """CPU-only Linux (CI train job): no bf16, so TrainingArguments does not raise; CUDA Ampere+/MPS: bf16."""
    monkeypatch.setattr(configs, "is_torch_bf16_gpu_available", lambda: supported)
    # use_cpu: TrainingArguments re-checks the real hardware, and accepts bf16 (CPU autocast) only with use_cpu
    for c in (configs.grpo_config("qwen3-vl-2b", str(tmp_path), use_cpu=True),
              configs.sft_config("qwen3-vl-2b", str(tmp_path), use_cpu=True)):
        assert c.bf16 is supported and not c.fp16
    c = configs.grpo_config("qwen3-vl-2b", str(tmp_path), bf16=False, max_steps=3)  # overrides win
    assert (c.bf16, c.max_steps) == (False, 3)
    assert not configs.sft_config("qwen3-vl-2b", str(tmp_path), fp16=True, use_cpu=True).bf16


def test_configs_build_on_this_machine(tmp_path):
    """Without overrides (CPU-only Linux CI included): the real hardware check decides bf16."""
    for m in configs.MODELS:
        assert configs.grpo_config(m, str(tmp_path)).bf16 == configs.is_torch_bf16_gpu_available()
        assert configs.sft_config(m, str(tmp_path)).bf16 == configs.is_torch_bf16_gpu_available()
