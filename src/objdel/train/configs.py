"""TRL 1.16 SFT and GRPO configurations for the three screened models (plan D8, D10, D12; lab t_stack).

GRPO (D8): loss_type="dr_grpo" (constant length normaliser, no std scaling), clip epsilon 0.2 / epsilon_high 0.28
(DAPO clip-higher), beta = 1e-3 from the start (with full fine-tuning a beta of 0 at init loads no reference
model, so raising it later crashes; lab t_stack T3), mask_truncated_completions, constant learning rate,
G = num_generations = 8 (REWARD_V4.group); any other num_generations or num_generations_eval raises ValueError,
because TRLReward would score a TRL group of 16 as two v4 groups of 8 without noticing (break-it B2-5). Reward:
TRLReward's weight-0 loggers plus the aggregator with multi_objective_aggregation="sum_then_normalize" and
scale_rewards="none" (objdel.train.rewards).
One step = 32 prompts x 8 completions (generation_batch_size 256, plan §3.2): per_device_train_batch_size 8 with 32
accumulation steps on one GPU; per_device sizes are multiples of G so every process's reward call holds whole
groups at any world size. num_generations_eval = None (= 8): held-out evaluation samples groups like training,
which is what the Phase 2 GRPO-signal gate (duplicates, sigma above the floor at T = 1) measures.

max_completion_length: 1,536 tokens for Qwen3-VL and SmolVLM, 1,024 for LFM2.5-VL (5.2 vs 8.2 tokens per vertex,
D10). Issue I5: on the final prepared set 1 of 3,282 RIGA train targets is over budget for each of the three
models and 0 of 23,469 LIDC train targets (prepare_rl's token_budget, re-counted with each tokenizer; 10 / 3,546
before BinRushed1 lost its model image, B3-1). The budget stays; such answers are dropped from SFT
(objdel.train.dataset.sft_rows with the model's tokenizer and MAX_COMPLETION), see experiments/phase1/impl/M7_train.md.

Generation: HF generate for SmolVLM, because TRL + vLLM re-expands SmolVLM/Idefics3's pre-expanded image tokens
(issue #6294 and fix PR #6718 still open on 2026-10-19); vLLM colocate for LFM2.5-VL and Qwen3-VL, whose image
placeholder expands to one contiguous run (lab t_stack T5).

SFT: assistant-only loss (vision datasets need transformers >= 5.18), max_length=None (truncation can cut image
tokens), no packing. Full fine-tuning lr 5e-5, LoRA r16 on all linear layers lr 5e-4 (lab T4). LoRA only for
Qwen3-VL-2B (a 24 GB card does not hold full fine-tuning of 2B parameters with AdamW plus generation); SmolVLM-256M
and LFM2.5-VL-450M are trained in full. GRPO learning rates: 1e-6 full, 1e-5 LoRA (lab K12 used 1e-5 with LoRA).
Precision: bf16 exactly where transformers accepts it (transformers.utils.is_torch_bf16_gpu_available: CUDA with
bf16 support, i.e. Ampere or newer such as the Vast cards, MPS on macOS >= 14, XPU / NPU / TPU), else fp32. TRL's
own default (bf16=True unless fp16) makes TrainingArguments raise "Your setup doesn't support bf16/gpu" on a
CPU-only machine such as the Linux CI train job. fp16=True or an explicit bf16 override wins.

Images (D12): SmolVLM resizes every image so its longest edge equals size["longest_edge"] (64 tokens per 512 px
tile plus a global tile): 1,024 for RIGA crops (320 tokens; 500 px crops would otherwise be compressed into 64
tokens), 512 for LIDC (the frame image is already upscaled to 512, 64 tokens). LFM2.5-VL and Qwen3-VL keep their
processors' defaults (256 tokens at 512 px). Left padding for generation.

References:
- TRL GRPO trainer and config: https://huggingface.co/docs/trl/main/en/grpo_trainer ;
  https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_config.py
- TRL SFT trainer (VLMs: max_length=None; assistant_only_loss with vision datasets): https://huggingface.co/docs/trl/main/en/sft_trainer
- Dr. GRPO, "Understanding R1-Zero-Like Training": https://arxiv.org/abs/2503.20783
- DAPO (clip-higher 0.28, dynamic sampling): https://arxiv.org/abs/2503.14476
- TRL issue #6294 / PR #6718 (SmolVLM + vLLM): https://github.com/huggingface/trl/issues/6294 ,
  https://github.com/huggingface/trl/pull/6718
- SmolVLM processor (longest_edge, 512 px tiles): https://huggingface.co/HuggingFaceTB/SmolVLM-256M-Instruct ;
  sandbox measurement in experiments/phase1/impl/M7_train.md
- transformers 5.20 TrainingArguments.__post_init__ (bf16 check) and utils.is_torch_bf16_gpu_available:
  https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/training_args.py
- Lab: experiments/phase1/lab/t_stack/REPORT.md, results/INTEGRATION.md (K12, R7).
"""
from peft import LoraConfig
from transformers.utils import is_torch_bf16_gpu_available
from trl import GRPOConfig, SFTConfig

from objdel.reward import REWARD_V4

from .rewards import LOGGED

MODELS = {
    "smolvlm-256m": "HuggingFaceTB/SmolVLM-256M-Instruct",
    "lfm2.5-vl-450m": "LiquidAI/LFM2.5-VL-450M",
    "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct",
}
MAX_COMPLETION = {"smolvlm-256m": 1536, "lfm2.5-vl-450m": 1024, "qwen3-vl-2b": 1536}
_LORA = {"qwen3-vl-2b"}


def lora(model: str) -> LoraConfig | None:
    """LoRA r16 on all linear layers for the models trained with adapters, else None (full fine-tuning)."""
    if model not in _LORA:
        return None
    return LoraConfig(r=16, lora_alpha=32, target_modules="all-linear", task_type="CAUSAL_LM")


def processor_kwargs(model: str, dataset: str) -> dict:
    """AutoProcessor.from_pretrained(MODELS[model], **processor_kwargs(model, dataset))."""
    size = {"size": {"longest_edge": {"riga": 1024, "lidc": 512}[dataset]}} if model == "smolvlm-256m" else {}
    return size | {"padding_side": "left"}


def _bf16(overrides: dict) -> bool:
    """bf16 where TrainingArguments' own check (__post_init__) accepts it, never together with fp16."""
    return is_torch_bf16_gpu_available() and not overrides.get("fp16", False)


def sft_config(model: str, output_dir: str, **overrides) -> SFTConfig:
    return SFTConfig(**{
        "output_dir": output_dir, "bf16": _bf16(overrides), "assistant_only_loss": True, "max_length": None,
        "packing": False,
        "learning_rate": 5e-4 if model in _LORA else 5e-5, "lr_scheduler_type": "constant", "warmup_steps": 0,
        "remove_unused_columns": False, "report_to": "none", "seed": 0} | overrides)


def grpo_config(model: str, output_dir: str, **overrides) -> GRPOConfig:
    g = REWARD_V4.group
    for k in ("num_generations", "num_generations_eval"):
        if overrides.get(k) not in (None, g):
            raise ValueError(f"{k}={overrides[k]}: reward v4 scores groups of exactly {g} (REWARD_V4.group)")
    return GRPOConfig(**{
        "output_dir": output_dir, "bf16": _bf16(overrides), "num_generations": g, "per_device_train_batch_size": g,
        "gradient_accumulation_steps": 32, "per_device_eval_batch_size": 2 * g,
        "max_completion_length": MAX_COMPLETION[model], "temperature": 1.0,
        "loss_type": "dr_grpo", "epsilon": 0.2, "epsilon_high": 0.28, "beta": 1e-3, "mask_truncated_completions": True,
        "learning_rate": 1e-5 if model in _LORA else 1e-6, "lr_scheduler_type": "constant", "warmup_steps": 0,
        "multi_objective_aggregation": "sum_then_normalize", "scale_rewards": "none",
        "reward_weights": [0.0] * len(LOGGED) + [1.0],  # TRLReward.funcs() order: loggers, then AGGREGATOR
        "use_vllm": model != "smolvlm-256m", "vllm_mode": "colocate",
        "remove_unused_columns": False,  # key, W, H, names reach the reward functions
        "log_completions": True, "report_to": "none", "seed": 0} | overrides)

