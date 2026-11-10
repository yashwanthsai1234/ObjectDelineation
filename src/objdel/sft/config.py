"""SFT configuration: sft_config and the WSD schedule fractions WARMUP and DECAY (D5). The reasons for each setting
are in objdel.policy.models."""

from typing import TYPE_CHECKING

from objdel.policy.models import _GRAD_CKPT, _LORA, NO_SYNC, _bf16

if TYPE_CHECKING:
    from trl import SFTConfig


WARMUP, DECAY = 0.03, 0.15  # D5: of the stable run's steps; of the branch point's steps


def sft_config(model: str, output_dir: str, steps: int | None = None, branch_at: int | None = None,
               **overrides) -> "SFTConfig":
    """The bundle's SFT settings (§6.0). Schedule (D5, WSD): linear warmup over 3 % of the stable run's `steps`,
    then constant; `branch_at=s` is a decay branch from the stable run's step-s checkpoint, decaying linearly to 0
    over 15 % of s (max_steps = s + that), identical to the stable run up to s. steps=None (smoke tests): 3 %
    warmup of max_steps, no decay. Loader (D4): group_by_length on the `length` column, 6 persistent workers,
    prefetch 2 (with dataloader_num_workers=0 neither, as TrainingArguments requires; BK3-1); gradient checkpointing
    only for Qwen3-VL-2B. lr: 5e-4 with LoRA (the models in _LORA), else 5e-5; an explicit learning_rate wins."""
    if branch_at is not None and steps is None:
        raise ValueError("a decay branch needs the stable run's steps")
    warmup = WARMUP if steps is None else round(WARMUP * steps)
    stable_end = branch_at if branch_at is not None else steps
    decay = 0 if branch_at is None else round(DECAY * branch_at)
    sched = {"num_decay_steps": decay, "decay_type": "linear", "min_lr_ratio": 0.0}  # stable = the rest
    from trl import SFTConfig  # here: the sampler and evaluation import this module without trl.trainer (E4)

    workers = overrides.get("dataloader_num_workers", 6)
    return SFTConfig(**{
        "output_dir": output_dir, "bf16": _bf16(overrides), "assistant_only_loss": True, "max_length": None,
        "packing": False, "learning_rate": 5e-4 if model in _LORA else 5e-5,
        "lr_scheduler_type": "warmup_stable_decay", "warmup_steps": warmup, "lr_scheduler_kwargs": sched,
        **({"max_steps": stable_end + decay} if steps is not None else {}),
        "train_sampling_strategy": "group_by_length", "length_column_name": "length",
        "dataloader_num_workers": workers, "dataloader_persistent_workers": workers > 0,
        "dataloader_prefetch_factor": 2 if workers else None,
        **NO_SYNC,
        "gradient_checkpointing": model in _GRAD_CKPT,
        "remove_unused_columns": False, "report_to": "none", "seed": 0, "data_seed": 0} | overrides)
