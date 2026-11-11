"""GRPO configuration: grpo_config, rl_lr, the numerics presets (B2.2, ADR-0013) and the soft overlong check (C7).
The reasons for each setting are in objdel.policy.models."""

from typing import TYPE_CHECKING

from objdel.policy import models as policy_models
from objdel.policy.models import NO_SYNC, _bf16
from objdel.verify.reward.advantages import LOGGED
from objdel.verify.reward.spec import REWARD_V4

if TYPE_CHECKING:
    from trl import GRPOConfig


LONGEST_OUTLINE = 1182  # LIDC train doctor outline, end token counted (A1 F9)


def check_soft_overlong(limit: int, cache: int) -> None:
    """DAPO's soft overlong penalty starts at limit - cache tokens: it must start above the longest doctor outline
    + 32, or a correct doctor-length answer is punished (C7 option O)."""
    if limit - cache < LONGEST_OUTLINE + 32:
        raise ValueError(f"soft overlong penalty starts at {limit - cache} tokens; "
                         f"it needs limit - cache >= {LONGEST_OUTLINE + 32}")


def rl_lr(model: str) -> float:
    """The GRPO learning rate of every real run: 5e-5 (every model of policy.models.MODELS is LoRA, ADR-0008). An
    unknown model raises: the full fine-tuning rate 5e-6 is a frozen record (prereg.G_PROBE), never a default."""
    if model not in policy_models.MODELS:
        raise ValueError(f"rl_lr: unknown model {model!r}; the models are {sorted(policy_models.MODELS)}")
    return policy_models.MODELS[model].rl_lr


# Numerics presets (B2.2, ADR-0013): the importance-sampling mode and dtype of a run. The shared parts are named first
# (a literal cannot refer to itself). Only the keys a preset changes are listed; vllm_lora is ours (see grpo_config).
_SEQMASK = {"vllm_importance_sampling_mode": "sequence_mask", "vllm_importance_sampling_clip_max": 3.0,
            "vllm_importance_sampling_clip_min": None}  # TRL's defaults, written out

_TIS = {"vllm_importance_sampling_mode": "token_truncate", "vllm_importance_sampling_clip_max": 2.0,
        "vllm_importance_sampling_clip_min": None}  # C = 2, SkyRL's token default

NUMERICS = {
    "bf16-seqmask": _SEQMASK,  # ADR option C: the control (B1 and B3)
    "bf16-tis": _TIS,  # option B
    "bf16-tis-opsm": {**_TIS, "off_policy_mask_threshold": 0.5},  # option E; train_rl --opsm-start-delta overrides 0.5
    "fp16-seqmask": {**_SEQMASK, "bf16": False, "fp16": True},  # option A numerics, the control's IS mode
    "lora-seqmask": {**_SEQMASK, "vllm_lora": True},  # option F, the control's IS mode
}

DEFAULT_NUMERICS = "bf16-seqmask"


def grpo_config(model: str, output_dir: str, reward=REWARD_V4, numerics: str = DEFAULT_NUMERICS,
                **overrides) -> "GRPOConfig":
    """reward: the reward spec (its group sets num_generations). numerics: a key of NUMERICS; explicit overrides win
    over the preset. vllm_lora is not a TRL field: it is set as an attribute of the returned config, for FastGRPO to
    read (False elsewhere)."""
    if numerics not in NUMERICS:
        raise ValueError(f"numerics {numerics!r}: not one of {sorted(NUMERICS)}")
    preset = dict(NUMERICS[numerics])
    lora = preset.pop("vllm_lora", False)
    overrides = preset | overrides
    g = reward.group
    for k in ("num_generations", "num_generations_eval"):
        if overrides.get(k) not in (None, g):
            raise ValueError(f"{k}={overrides[k]}: {reward.name} scores groups of exactly {g} (its group)")
    from trl import GRPOConfig  # here, as in sft_config

    m = policy_models.get(model)
    cfg = GRPOConfig(**{
        "output_dir": output_dir, "bf16": _bf16(overrides), "num_generations": g,
        "per_device_train_batch_size": m.groups_per_device * g, "gradient_accumulation_steps": m.grad_accum,
        "per_device_eval_batch_size": m.groups_per_device * g, "gradient_checkpointing": m.grad_ckpt,
        "max_completion_length": m.max_completion, "temperature": 1.0,
        "loss_type": "dr_grpo", "epsilon": 0.2, "epsilon_high": 0.28, "beta": 1e-3, "mask_truncated_completions": True,
        "learning_rate": rl_lr(model), "lr_scheduler_type": "constant", "warmup_steps": 0,
        "multi_objective_aggregation": "sum_then_normalize", "scale_rewards": "none",
        "reward_weights": [0.0] * len(LOGGED) + [1.0],  # TRLReward.funcs() order: loggers, then AGGREGATOR
        "use_vllm": True, "vllm_mode": "colocate", "vllm_gpu_memory_utilization": m.vllm_train[0],
        "vllm_enable_sleep_mode": True, "vllm_max_model_length": m.vllm_train[1], **NO_SYNC,
        "remove_unused_columns": False,  # key, W, H, names reach the reward functions
        "log_completions": True, "report_to": "none", "seed": 0} | overrides)
    train, spg = cfg.per_device_train_batch_size, cfg.steps_per_generation
    b = {"per_device_train_batch_size": train, "per_device_eval_batch_size": cfg.per_device_eval_batch_size,
         "per_device_train_batch_size * steps_per_generation": train * spg}
    for k, v in b.items():  # TRLReward scores whole groups: refuse here, not at the first held-out eval
        if v % g:
            raise ValueError(f"{k}={v}: must be a multiple of the group size {g} ({reward.name}'s group)")
    cfg.vllm_lora = lora
    return cfg
