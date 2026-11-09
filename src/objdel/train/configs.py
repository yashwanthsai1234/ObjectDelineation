"""TRL (version pinned in pyproject.toml) SFT and GRPO configurations for the screened models
(plan D8, D10, D12; lab t_stack).

GRPO (D8): loss_type="dr_grpo" (constant length normaliser, no std scaling), clip epsilon 0.2 / epsilon_high 0.28
(DAPO clip-higher), beta = 1e-3 from the start (with full fine-tuning a beta of 0 at init loads no reference
model, so raising it later crashes; lab t_stack T3), mask_truncated_completions, constant learning rate,
G = num_generations = 8 (REWARD_V4.group); any other num_generations or num_generations_eval raises ValueError,
because TRLReward would score a TRL group of 16 as two v4 groups of 8 without noticing (break-it B2-5). Reward:
TRLReward's weight-0 loggers plus the aggregator with multi_objective_aggregation="sum_then_normalize" and
scale_rewards="none" (objdel.train.rewards).
One step = 32 prompts x 8 completions (generation_batch_size 256, plan §3.2): per_device_train_batch_size 16 (two
groups) with 16 accumulation steps on one GPU (perf F2, 2026-11-06: half the CPU launch work per answer of 8 x 32;
the profiling run measures 8 / 16 / 32); per_device sizes are multiples of G so every process's reward call holds
whole groups at any world size. vLLM colocate (perf F1): vllm_gpu_memory_utilization 0.35, sleep mode on (vLLM
frees its memory during the update). num_generations_eval = None (= 8): held-out evaluation samples groups like
training, which is what the Phase 2 GRPO-signal gate (duplicates, sigma above the floor at T = 1) measures.

max_completion_length: 1,536 tokens for Qwen3-VL, 1,024 for LFM2.5-VL (5.2 vs 8.2 tokens per vertex,
D10). Issue I5: on the final prepared set 1 of 3,282 RIGA train targets is over budget for each of the three
models and 0 of 23,469 LIDC train targets (prepare_rl's token_budget, re-counted with each tokenizer; 10 / 3,546
before BinRushed1 lost its model image, B3-1). The budget stays; such answers are dropped from SFT
(objdel.train.dataset.sft_rows with the model's tokenizer and MAX_COMPLETION), see experiments/phase1/impl/M7_train.md.

Generation: vLLM colocate for every model (each image placeholder expands to one contiguous run, lab t_stack T5);
a CPU caller passes use_vllm=False for HF generate.

Both configs (perf F3): logging_nan_inf_filter off (an isnan/isinf host sync on every micro-batch; a NaN loss still
shows at the log steps) and non-blocking host-to-GPU copies (pinned memory is on).

SFT: assistant-only loss (vision datasets need transformers >= 5.18; TRL swaps in its own training template for
the processor's chat template, qwen3_vl_training for Qwen3-VL and lfm2_training for LFM2.5-VL), max_length=None
(truncation can cut image tokens), no packing. Full fine-tuning lr 5e-5, LoRA r16 on the language model's linear
layers lr 5e-4 (lab T4). LoRA only for Qwen3-VL-2B (full fine-tuning with AdamW needs ~34 GB for weights, gradients
and optimizer state, more than a 24 GB card has; LoRA with a frozen bf16 base peaks far below); LFM2.5-VL-450M
is trained in full. GRPO learning rate: `rl_lr` (5e-5 LoRA, 5e-6 full), the default of `grpo_config` and the value
`g_probe.make_trainer` uses; `tests/experiments/test_prereg_lr.py` pins it to the pre-registered `prereg.G_PROBE`.
Precision: bf16 exactly where transformers accepts it (transformers.utils.is_torch_bf16_gpu_available: CUDA with
bf16 support, i.e. Ampere or newer such as the Vast cards, MPS on macOS >= 14, XPU / NPU / TPU), else fp32. TRL's
own default (bf16=True unless fp16) makes TrainingArguments raise "Your setup doesn't support bf16/gpu" on a
CPU-only machine such as the Linux CI train job. fp16=True or an explicit bf16 override wins.

Images (D12): every model keeps its processor's default image size (256 tokens at 512 px). Left padding for
generation.

References:
- TRL GRPO trainer and config: https://huggingface.co/docs/trl/main/en/grpo_trainer ;
  https://github.com/huggingface/trl/blob/v1.17.1/trl/trainer/grpo_config.py
- TRL SFT trainer (VLMs: max_length=None; assistant_only_loss with vision datasets): https://huggingface.co/docs/trl/main/en/sft_trainer
- Dr. GRPO, "Understanding R1-Zero-Like Training": https://arxiv.org/abs/2503.20783
- DAPO (clip-higher 0.28, dynamic sampling): https://arxiv.org/abs/2503.14476
- transformers 5.20 TrainingArguments.__post_init__ (bf16 check) and utils.is_torch_bf16_gpu_available:
  https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/training_args.py
- Lab: archive/phase1_lab/t_stack/REPORT.md, results/INTEGRATION.md (K12, R7).
"""
import re

from peft import LoraConfig
from transformers.utils import is_torch_bf16_gpu_available
from trl import GRPOConfig, SFTConfig

from objdel.reward import REWARD_V4

from .rewards import LOGGED

MODELS = {
    "lfm2.5-vl-450m": "LiquidAI/LFM2.5-VL-450M",
    "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct",
}
MAX_COMPLETION = {"lfm2.5-vl-450m": 1024, "qwen3-vl-2b": 1536}
_LORA = {"qwen3-vl-2b"}
_GRAD_CKPT = {"qwen3-vl-2b"}
# vision tower + projector of every model (recorded module names in tests/test_p2_train.py;
# archive/phase2_lab/B1/common.VISION)
VISION = r".*(visual|vision_model|vision_tower|connector|multi_modal_projector).*"
WARMUP, DECAY = 0.03, 0.15  # D5: of the stable run's steps; of the branch point's steps


def dtype(lora: bool) -> str:
    """from_pretrained(dtype=...). Full fine-tuning loads fp32 weights, AdamW's master copy, and trains under
    bf16 autocast (sft_config's bf16=True): updates below bf16 resolution are not lost (S1-13). Only a LoRA base,
    which is frozen, loads in bf16."""
    return "bfloat16" if lora else "float32"


def lora(model: str) -> LoraConfig | None:
    """LoRA r16 on the language model's linear layers (vision tower and projector excluded, D6 / D7) for the models
    trained with adapters, else None (full fine-tuning)."""
    if model not in _LORA:
        return None
    return LoraConfig(r=16, lora_alpha=32, target_modules="all-linear", exclude_modules=VISION, task_type="CAUSAL_LM")


def param_groups(model, base_lr: float) -> list[dict]:
    """Optimizer groups for full fine-tuning (D7): the language model at base_lr, vision tower + projector at 0.1x;
    frozen parameters are left out. torch.optim.AdamW(param_groups(model, lr)) -> Trainer(optimizers=(opt, None))."""
    lm, vision = [], []
    for name, p in model.named_parameters():
        if p.requires_grad:
            (vision if re.fullmatch(VISION, name) else lm).append(p)
    return [{"params": lm, "lr": base_lr}, {"params": vision, "lr": 0.1 * base_lr}]


def sampling(model: str, api: str = "vllm") -> dict:
    """Every sampling field explicit (D13): T 1 and no filter or penalty, overriding the cards' defaults (Qwen3-VL
    T 0.7 / top_p 0.8 / top_k 20; LFM T 0.1 / min_p 0.15 / repetition 1.05). One source of values,
    spelled for api "vllm" (SamplingParams), "grpo" (GRPOConfig) or "hf" (GenerationConfig), see as_names."""
    s = {"temperature": 1.0, "top_p": 1.0, "top_k": 0, "min_p": 0.0, "repetition_penalty": 1.0,
         "presence_penalty": 0.0, "frequency_penalty": 0.0, "max_tokens": MAX_COMPLETION[model]}
    return as_names(s, api)


_RENAME = {"grpo": {"max_tokens": "max_completion_length"},
           "hf": {"max_tokens": "max_new_tokens", "stop_token_ids": "eos_token_id"}}


def as_names(s: dict, api: str) -> dict:
    """vLLM SamplingParams fields (e.g. a prompts JSONL's "sampling") -> api's names. GRPOConfig and GenerationConfig
    have no presence / frequency penalty, so both must be 0 (ValueError otherwise); the request seed is dropped
    (HF generate seeds torch), and HF samples with do_sample=True."""
    if api == "vllm":
        return dict(s)
    if s.get("presence_penalty") or s.get("frequency_penalty"):
        raise ValueError(f"{api} has no presence_penalty / frequency_penalty: {s}")
    out = {_RENAME[api].get(k, k): v for k, v in s.items()
           if k not in ("presence_penalty", "frequency_penalty", "seed")}
    return out | {"do_sample": True} if api == "hf" else out


def processor_kwargs(model: str, dataset: str) -> dict:
    """AutoProcessor.from_pretrained(MODELS[model], **processor_kwargs(model, dataset))."""
    return {"padding_side": "left"}


def _bf16(overrides: dict) -> bool:
    """bf16 where TrainingArguments' own check (__post_init__) accepts it, never together with fp16."""
    return is_torch_bf16_gpu_available() and not overrides.get("fp16", False)


NO_SYNC = {"logging_nan_inf_filter": False, "accelerator_config": {"non_blocking": True}}  # perf F3


def sft_config(model: str, output_dir: str, steps: int | None = None, branch_at: int | None = None,
               **overrides) -> SFTConfig:
    """The bundle's SFT settings (§6.0). Schedule (D5, WSD): linear warmup over 3 % of the stable run's `steps`,
    then constant; `branch_at=s` is a decay branch from the stable run's step-s checkpoint, decaying linearly to 0
    over 15 % of s (max_steps = s + that), identical to the stable run up to s. steps=None (smoke tests): 3 %
    warmup of max_steps, no decay. Loader (D4): group_by_length on the `length` column, 6 persistent workers,
    prefetch 2 (with dataloader_num_workers=0 neither, as TrainingArguments requires; BK3-1); gradient checkpointing
    only for Qwen3-VL-2B. bf16 is autocast over the fp32 weights of dtype() for full fine-tuning. lr: 5e-4 with
    LoRA (the models in _LORA), else 5e-5; an explicit learning_rate wins."""
    if branch_at is not None and steps is None:
        raise ValueError("a decay branch needs the stable run's steps")
    warmup = WARMUP if steps is None else round(WARMUP * steps)
    stable_end = branch_at if branch_at is not None else steps
    decay = 0 if branch_at is None else round(DECAY * branch_at)
    sched = {"num_decay_steps": decay, "decay_type": "linear", "min_lr_ratio": 0.0}  # stable = the rest
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


def rl_lr(model: str) -> float:
    """The GRPO learning rate of every real run: 5e-5 for LoRA models, 5e-6 for full fine-tuning."""
    return 5e-5 if model in _LORA else 5e-6


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


def grpo_config(model: str, output_dir: str, numerics: str = DEFAULT_NUMERICS, **overrides) -> GRPOConfig:
    """numerics: a key of NUMERICS; explicit overrides win over the preset. vllm_lora is not a TRL field: it is set
    as an attribute of the returned config, for FastGRPO to read (False elsewhere)."""
    if numerics not in NUMERICS:
        raise ValueError(f"numerics {numerics!r}: not one of {sorted(NUMERICS)}")
    preset = dict(NUMERICS[numerics])
    lora = preset.pop("vllm_lora", False)
    overrides = preset | overrides
    g = REWARD_V4.group
    for k in ("num_generations", "num_generations_eval"):
        if overrides.get(k) not in (None, g):
            raise ValueError(f"{k}={overrides[k]}: reward v4 scores groups of exactly {g} (REWARD_V4.group)")
    cfg = GRPOConfig(**{
        "output_dir": output_dir, "bf16": _bf16(overrides), "num_generations": g,
        # P1, RTX 5090 (2026-11-06, perf/p1_results.md): 16 x 16 44.3 s/step; 32 x 8 59.1; no checkpointing: OOM
        "per_device_train_batch_size": 2 * g, "gradient_accumulation_steps": 16, "per_device_eval_batch_size": 2 * g,
        "gradient_checkpointing": True,
        "max_completion_length": MAX_COMPLETION[model], "temperature": 1.0,
        "loss_type": "dr_grpo", "epsilon": 0.2, "epsilon_high": 0.28, "beta": 1e-3, "mask_truncated_completions": True,
        "learning_rate": rl_lr(model), "lr_scheduler_type": "constant", "warmup_steps": 0,
        "multi_objective_aggregation": "sum_then_normalize", "scale_rewards": "none",
        "reward_weights": [0.0] * len(LOGGED) + [1.0],  # TRLReward.funcs() order: loggers, then AGGREGATOR
        "use_vllm": True, "vllm_mode": "colocate", "vllm_gpu_memory_utilization": 0.35,
        "vllm_enable_sleep_mode": True, "vllm_max_model_length": 4096, **NO_SYNC,
        "remove_unused_columns": False,  # key, W, H, names reach the reward functions
        "log_completions": True, "report_to": "none", "seed": 0} | overrides)
    cfg.vllm_lora = lora
    return cfg

