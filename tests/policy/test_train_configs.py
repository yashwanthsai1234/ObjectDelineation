"""objdel.policy.models, rl.config and sft.config (the old train/configs.py): plan D8 (GRPO), SFT settings, D12 image
sizes, the I5 completion budget; CLEAN7 (the 256M model retired: CT-05, CT-06, CT-07); CLEAN8 (the 450M model and full
fine-tuning retired: CT-62, CT-88)."""
import ast
import dataclasses
import re
import subprocess
from dataclasses import asdict
from pathlib import Path

import pytest

pytest.importorskip("trl")

from objdel.policy import models as policy_models  # noqa: E402
from objdel.rl import config as rl_config  # noqa: E402
from objdel.sft import config as sft_config  # noqa: E402
from objdel.verify.reward.advantages import LOGGED  # noqa: E402
from objdel.verify.reward.spec import REWARD_V4  # noqa: E402

OLD = "smol" + "vlm-256m"  # the retired key, built so this file never holds the word
OLD8 = "l" + "fm2.5-vl-450m"  # CLEAN8's retired key, built the same way
ROOT = Path(__file__).resolve().parents[2]
KEPT_MODELS = {"qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct"}  # CT-62: MODELS of the parent commit minus OLD8
KEPT_MAX = {"qwen3-vl-2b": 1536}  # CT-62: MAX_COMPLETION of the parent commit minus OLD8
NEEDS_ARCHIVE = pytest.mark.skipif(not (ROOT / "archive").exists(), reason="archive/ is local-only (git-ignored)")


def test_models_and_budgets():
    """CT-05: no exact set (a 4B row may join); every model has a budget, the retired keys are gone."""
    assert policy_models.MODELS["qwen3-vl-2b"] == "Qwen/Qwen3-VL-2B-Instruct"
    for old in (OLD, OLD8):
        assert old not in policy_models.MODELS and old not in policy_models.MAX_COMPLETION
    assert set(policy_models.MODELS) == set(policy_models.MAX_COMPLETION)


@pytest.mark.parametrize("model", list(policy_models.MODELS))
def test_grpo_config_is_plan_d8(model, tmp_path):
    c = rl_config.grpo_config(model, str(tmp_path))
    assert (c.loss_type, c.epsilon, c.epsilon_high, c.beta) == ("dr_grpo", 0.2, 0.28, 1e-3)
    assert c.mask_truncated_completions and c.lr_scheduler_type == "constant" and c.warmup_steps == 0
    assert c.num_generations == REWARD_V4.group and c.num_generations_eval is None and c.temperature == 1.0
    assert (c.scale_rewards, c.multi_objective_aggregation) == ("none", "sum_then_normalize")
    assert c.reward_weights == [0.0] * len(LOGGED) + [1.0]
    assert c.max_completion_length == policy_models.MAX_COMPLETION[model]
    assert c.use_vllm is True
    assert not c.remove_unused_columns and c.log_completions
    # whole groups on every process, at any world size (tests/rl/test_train_trl.py simulates 1, 2, 4 processes)
    assert c.per_device_train_batch_size % c.num_generations == 0
    assert c.per_device_eval_batch_size % c.num_generations == 0
    assert c.generation_batch_size == 256 and c.steps_per_generation == c.gradient_accumulation_steps


@pytest.mark.parametrize("model", list(policy_models.MODELS))
def test_sft_config(model, tmp_path):
    c = sft_config.sft_config(model, str(tmp_path))
    assert c.assistant_only_loss and c.max_length is None and not c.packing
    assert c.learning_rate == 5e-4  # every model is LoRA


def test_lora_for_every_model():
    q = policy_models.lora("qwen3-vl-2b")
    # Phase 2 (plan D6 / D7, M2): language model only, vision tower and projector excluded
    assert (q.r, q.lora_alpha, q.target_modules, q.exclude_modules) == (16, 32, "all-linear", policy_models.VISION)


def test_processors_keep_their_default_image_sizes():
    """CT-05 (D12 after CLEAN7): left padding only, for every model and any dataset string (no size key, no
    lookup on the dataset)."""
    for m in policy_models.MODELS:
        for d in ("lidc", "riga", "anything"):
            assert policy_models.processor_kwargs(m, d) == {"padding_side": "left"}
    with pytest.raises(KeyError):
        policy_models.sampling(OLD)  # the retired key has no budget


GRPO_KEYS = ("loss_type", "epsilon", "epsilon_high", "beta", "num_generations", "per_device_train_batch_size",
             "max_completion_length", "vllm_gpu_memory_utilization", "vllm_enable_sleep_mode", "vllm_max_model_length",
             "learning_rate", "reward_weights")
GRPO_PRE_CLEAN = {  # captured at 08e6455 (asdict of grpo_config(m, out, use_cpu=True, use_vllm=False))
    m: {"loss_type": "dr_grpo", "epsilon": 0.2, "epsilon_high": 0.28, "beta": 0.001, "num_generations": 8,
        "per_device_train_batch_size": 16, "max_completion_length": n, "vllm_gpu_memory_utilization": 0.35,
        "vllm_enable_sleep_mode": True, "vllm_max_model_length": 4096, "learning_rate": lr,
        "reward_weights": [0.0] * 8 + [1.0]}  # B3.5 adds the weight-0 loc logger
    for m, n, lr in (("qwen3-vl-2b", 1536, 5e-05),)}  # E1: rl_lr, pre-registered; CLEAN8 drops the 450M row


@pytest.mark.parametrize("model", sorted(GRPO_PRE_CLEAN))
def test_grpo_config_values_and_the_use_vllm_override(model, tmp_path):
    """CT-06: vLLM for every model unless the caller overrides it; the plan values are the pre-CLEAN ones."""
    c = rl_config.grpo_config(model, str(tmp_path), use_cpu=True, use_vllm=False)
    assert c.use_vllm is False
    assert {k: v for k, v in asdict(c).items() if k in GRPO_KEYS} == GRPO_PRE_CLEAN[model]
    assert rl_config.grpo_config(model, str(tmp_path), use_vllm=True).use_vllm is True
    assert rl_config.grpo_config(model, str(tmp_path), use_cpu=True, use_vllm=None).use_vllm is None  # pre-CLEAN
    assert rl_config.grpo_config(model, str(tmp_path), use_cpu=True, use_vllm=0).use_vllm == 0  # pre-CLEAN
    with pytest.raises(ValueError):
        rl_config.grpo_config(model, str(tmp_path), num_generations=16)


SMOL = re.compile("smol" + "vlm|hugging" + "facetb", re.I)
R450 = re.compile(r"(?<![a-z])" + "l" + "fm|liquid" + "ai", re.I)
TEXT = (".py", ".yaml", ".md", ".toml")


def test_forbidden_word_needles():
    """CT-07 checker part: the needles match the retired names and no ordinary word."""
    for t in ("selfmade", "shelfmark", "aGVsZm1Gx"):
        assert not R450.search(t)
    for t in ("L" + "FM2.5-VL", OLD8, "# uses " + "L" + "FM", "Liquid" + "AI/x"):
        assert R450.search(t)
    for t in ("Smol" + "VLM", "smol" + "vlm_256m", "SMOL" + "VLM", "Hugging" + "FaceTB/x"):
        assert SMOL.search(t)


def scanned() -> list[str]:
    """Text sources of src/, experiments/phase2/*.py and tests/ (git ls-files -co, so no __pycache__)."""
    root = Path(__file__).resolve().parents[2]
    files = subprocess.run(["git", "-C", str(root), "ls-files", "-z", "-co", "--exclude-standard", "src", "tests",
                            "experiments/phase2"], capture_output=True, text=True, check=True).stdout.split("\0")
    keep = [f for f in files if f.endswith(TEXT) and (root / f).is_file()
            and (f.startswith(("src/", "tests/")) or re.fullmatch(r"experiments/phase2/[^/]+\.py", f))]
    checker = "tests/experiments/test_final_sweep.py"  # tests the final sweep's own retired-model rule
    return [f for f in keep if f != checker
            and any(r.search((root / f).read_text(errors="replace")) for r in (SMOL, R450))]


def test_no_retired_model_text_in_code():
    """CT-07 repo part, Waves 1 and 3 (CLEAN-G10): no file in the scope names a retired model."""
    assert scanned() == []


def _defs(path: Path) -> set[str]:
    return {n.name for n in ast.walk(ast.parse(path.read_text())) if isinstance(n, ast.FunctionDef)}


def test_retired_450m_model_and_full_fine_tuning_gone():
    """CT-62: the model tables are the parent commit's minus the retired key (every other row survives); every model
    is LoRA; the full fine-tuning branches are gone (the pre-registered record: tests/experiments/test_prereg_lr.py)."""
    assert policy_models.MODELS == KEPT_MODELS and policy_models.MAX_COMPLETION == KEPT_MAX
    assert policy_models._LORA == frozenset(policy_models.MODELS) and policy_models._GRAD_CKPT
    assert rl_config.rl_lr("qwen3-vl-2b") == 5e-5
    for f in ("policy/models.py", "rl/config.py", "sft/config.py"):
        assert not {"param_groups", "dtype"} & _defs(ROOT / "src/objdel" / f)
    sft = ROOT / "src/objdel/sft/train.py"
    assert "_load_from_checkpoint" not in _defs(sft) and "param_groups" not in sft.read_text()
    assert not (ROOT / "experiments/phase2/configs" / "r1_l" "fm.yaml").exists()  # the folder left git in E6.6


ARCHIVED = ("r0_parity.py", "hf_rate.py", "c2_analyze.py")  # CLEAN8 step 4 (project_r1.py went in E6.6)


def test_phase2_run_scripts_left_the_live_tree():
    """CT-88, the part a clone can check: the scripts and their plans are gone and nothing live names them."""
    plans = ("r0.yaml", "r1_l" "fm.yaml")
    assert not [n for n in ARCHIVED if (ROOT / "experiments/phase2" / n).exists()]
    assert not [n for n in plans if (ROOT / "experiments/phase2/configs" / n).exists()]
    files = subprocess.run(["git", "ls-files", "src", "experiments", "tests", "specs"], cwd=ROOT,
                           capture_output=True, text=True, check=True).stdout.split()
    mods = [f"experiments.phase2.{n[:-3]}" for n in ARCHIVED]
    assert [f for f in files if (ROOT / f).is_file() and any(m in (ROOT / f).read_text(errors="ignore")
                                                              for m in mods)] == []


@NEEDS_ARCHIVE
def test_phase2_run_scripts_are_archived():
    """CT-88, the owner's checkout: the scripts and plans sit in archive/phase2_runs/, the plans with a README row."""
    assert all((ROOT / "archive/phase2_runs" / n).is_file() for n in ARCHIVED)
    readme = (ROOT / "archive/phase2_runs/configs/README.md").read_text()
    assert all((ROOT / "archive/phase2_runs/configs" / n).is_file() and f"| {n} |" in readme
               for n in ("r0.yaml", "r1_l" "fm.yaml"))


@pytest.mark.parametrize("supported", [True, False])
def test_bf16_only_where_the_hardware_supports_it(supported, tmp_path, monkeypatch):
    """CPU-only Linux (CI train job): no bf16, so TrainingArguments does not raise; CUDA Ampere+/MPS: bf16."""
    monkeypatch.setattr(policy_models, "is_torch_bf16_gpu_available", lambda: supported)
    # use_cpu: TrainingArguments re-checks the real hardware, and accepts bf16 (CPU autocast) only with use_cpu
    for c in (rl_config.grpo_config("qwen3-vl-2b", str(tmp_path), use_cpu=True),
              sft_config.sft_config("qwen3-vl-2b", str(tmp_path), use_cpu=True)):
        assert c.bf16 is supported and not c.fp16
    c = rl_config.grpo_config("qwen3-vl-2b", str(tmp_path), bf16=False, max_steps=3)  # overrides win
    assert (c.bf16, c.max_steps) == (False, 3)
    assert not sft_config.sft_config("qwen3-vl-2b", str(tmp_path), fp16=True, use_cpu=True).bf16


def test_configs_build_on_this_machine(tmp_path):
    """Without overrides (CPU-only Linux CI included): the real hardware check decides bf16."""
    for m in policy_models.MODELS:
        assert rl_config.grpo_config(m, str(tmp_path)).bf16 == policy_models.is_torch_bf16_gpu_available()
        assert sft_config.sft_config(m, str(tmp_path)).bf16 == policy_models.is_torch_bf16_gpu_available()


# B2.2 numerics presets
def test_default_preset_is_the_b1_control(tmp_path):
    """No numerics argument: today's TRL values, written out (sequence mask, cap 3, no OPSM, no vllm_lora)."""
    assert rl_config.DEFAULT_NUMERICS == "bf16-seqmask"
    c = rl_config.grpo_config("qwen3-vl-2b", str(tmp_path))
    assert (c.vllm_importance_sampling_mode, c.vllm_importance_sampling_clip_max) == ("sequence_mask", 3.0)
    assert c.vllm_importance_sampling_clip_min is None and c.off_policy_mask_threshold is None
    assert not c.fp16 and not getattr(c, "vllm_lora", False)


def test_numerics_names():
    assert set(rl_config.NUMERICS) == {"bf16-seqmask", "bf16-tis", "bf16-tis-opsm", "fp16-seqmask", "lora-seqmask"}


@pytest.mark.parametrize("name", list(rl_config.NUMERICS))
def test_every_preset_builds_a_grpo_config(name, tmp_path):
    c = rl_config.grpo_config("qwen3-vl-2b", str(tmp_path), numerics=name)
    assert not (c.bf16 and c.fp16)


def test_tis_presets_use_token_truncation_at_2(tmp_path):
    for name in ("bf16-tis", "bf16-tis-opsm"):
        c = rl_config.grpo_config("qwen3-vl-2b", str(tmp_path), numerics=name)
        assert (c.vllm_importance_sampling_mode, c.vllm_importance_sampling_clip_max) == ("token_truncate", 2.0)
    assert rl_config.grpo_config("qwen3-vl-2b", str(tmp_path), numerics="bf16-tis").off_policy_mask_threshold is None
    assert rl_config.grpo_config("qwen3-vl-2b", str(tmp_path),
        numerics="bf16-tis-opsm").off_policy_mask_threshold == 0.5


def test_opsm_start_delta_override(tmp_path):
    c = rl_config.grpo_config("qwen3-vl-2b", str(tmp_path), numerics="bf16-tis-opsm", off_policy_mask_threshold=0.9)
    assert c.off_policy_mask_threshold == 0.9


def test_fp16_preset_is_fp16_never_with_bf16_and_keeps_the_control_is_mode(tmp_path):
    c = rl_config.grpo_config("qwen3-vl-2b", str(tmp_path), numerics="fp16-seqmask")
    assert c.fp16 and not c.bf16
    assert (c.vllm_importance_sampling_mode, c.vllm_importance_sampling_clip_max) == ("sequence_mask", 3.0)


def test_lora_preset_marks_the_config_and_never_passes_the_key_to_trl(tmp_path):
    c = rl_config.grpo_config("qwen3-vl-2b", str(tmp_path), numerics="lora-seqmask")
    assert c.vllm_lora is True
    assert "vllm_lora" not in {f.name for f in dataclasses.fields(c)}  # an attribute of ours, not a TRL field
    assert (c.vllm_importance_sampling_mode, c.vllm_importance_sampling_clip_max) == ("sequence_mask", 3.0)


def test_unknown_preset_is_refused(tmp_path):
    with pytest.raises(ValueError, match="bf16-fp8"):
        rl_config.grpo_config("qwen3-vl-2b", str(tmp_path), numerics="bf16-fp8")


@pytest.mark.gpu
def test_g05_f_check_record():
    """G05 box part: the F check record that `train_rl --f-check` wrote (`<run>/*/f_check.json`, $OBJDEL_RUN_DIR):
    `passed` and every check flag hold. The run itself already raised on a failure at step 15 or earlier."""
    import json
    import os
    from pathlib import Path

    from objdel.rl import f_check

    files = sorted(Path(os.environ.get("OBJDEL_RUN_DIR", "runs")).glob("*/f_check.json"))
    if not files:
        pytest.skip("no f_check.json: no run used lora-seqmask or fp16-seqmask yet")
    for f in files:
        f_check.check_record(json.loads(f.read_text()))
