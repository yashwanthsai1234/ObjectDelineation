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
import re

from peft import LoraConfig
from transformers.utils import is_torch_bf16_gpu_available
from trl import GRPOConfig, SFTConfig
from trl import chat_template_utils as ct

from objdel.reward import REWARD_V4

from .rewards import LOGGED

MODELS = {
    "smolvlm-256m": "HuggingFaceTB/SmolVLM-256M-Instruct",
    "lfm2.5-vl-450m": "LiquidAI/LFM2.5-VL-450M",
    "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct",
    "qwen3.5-0.8b": "Qwen/Qwen3.5-0.8B",
}
MAX_COMPLETION = {"smolvlm-256m": 1536, "lfm2.5-vl-450m": 1024, "qwen3-vl-2b": 1536, "qwen3.5-0.8b": 1024}
TEMPLATES = {"qwen3.5-0.8b": "qwen3_5_nothink_training"}  # TRL's assistant-only template the plan names
_LORA = {"qwen3-vl-2b"}
_BF16_BASE = {"qwen3.5-0.8b"}  # a frozen LoRA base may load in bf16 (B1); trained weights never do (S1-13)
_GRAD_CKPT = {"qwen3-vl-2b"}
# vision tower + projector of every model (recorded module names in tests/test_p2_train.py; lab/B1/common.VISION)
VISION = r".*(visual|vision_model|vision_tower|connector|multi_modal_projector).*"
WARMUP, DECAY = 0.03, 0.15  # D5: of the stable run's steps; of the branch point's steps


def dtype(model: str, lora: bool) -> str:
    """from_pretrained(dtype=...). Full fine-tuning loads fp32 weights, AdamW's master copy, and trains under
    bf16 autocast (sft_config's bf16=True): updates below bf16 resolution are not lost (S1-13). Only a LoRA base,
    which is frozen, loads in bf16, and only where B1 ran it so (Qwen3.5)."""
    return "bfloat16" if lora and model in _BF16_BASE else "float32"


def lora(model: str, force: bool = False) -> LoraConfig | None:
    """LoRA r16 on the language model's linear layers (vision tower and projector excluded, D6 / D7) for the models
    trained with adapters (or any model with force, e.g. Qwen3.5 if R0's peak memory says so), else None."""
    if model not in _LORA and not force:
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
    """Every sampling field explicit (D13): T 1 and no filter or penalty, overriding the cards' defaults (Qwen3.5
    T 0.7 / top_p 0.8 / top_k 20 / presence 1.5; LFM T 0.1 / min_p 0.15 / repetition 1.05). One source of values,
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


def training_template(model: str, processor) -> str | None:
    """The chat template SFTTrainer swaps in for assistant-only loss (TRL get_training_chat_template), checked
    against TEMPLATES (Qwen3.5: non-thinking, prefix-consistent with the generation prompt, B1 §4.5, B5 Q2);
    None for models without an entry (TRL's own choice)."""
    if model not in TEMPLATES:
        return None
    want = getattr(ct, TEMPLATES[model] + "_chat_template")
    try:
        got = ct.get_training_chat_template(processor)
    except ValueError:
        got = None
    if got != want:
        raise ValueError(f"{model}: the processor's chat template does not map to TRL's {TEMPLATES[model]}")
    return got


def processor_kwargs(model: str, dataset: str) -> dict:
    """AutoProcessor.from_pretrained(MODELS[model], **processor_kwargs(model, dataset))."""
    size = {"size": {"longest_edge": {"riga": 1024, "lidc": 512}[dataset]}} if model == "smolvlm-256m" else {}
    return size | {"padding_side": "left"}


def _bf16(overrides: dict) -> bool:
    """bf16 where TrainingArguments' own check (__post_init__) accepts it, never together with fp16."""
    return is_torch_bf16_gpu_available() and not overrides.get("fp16", False)


def sft_config(model: str, output_dir: str, steps: int | None = None, branch_at: int | None = None,
               lora: bool | None = None, **overrides) -> SFTConfig:
    """The bundle's SFT settings (§6.0). Schedule (D5, WSD): linear warmup over 3 % of the stable run's `steps`,
    then constant; `branch_at=s` is a decay branch from the stable run's step-s checkpoint, decaying linearly to 0
    over 15 % of s (max_steps = s + that), identical to the stable run up to s. steps=None (smoke tests): 3 %
    warmup of max_steps, no decay. Loader (D4): group_by_length on the `length` column, 6 persistent workers,
    prefetch 2 (with dataloader_num_workers=0 neither, as TrainingArguments requires; BK3-1); gradient checkpointing
    only for Qwen3-VL-2B. bf16 is autocast over the fp32 weights of dtype() for full fine-tuning. lr: 5e-4 with
    LoRA (lora=None: the model's default mode, as lora(model)), else 5e-5; an explicit learning_rate wins."""
    if branch_at is not None and steps is None:
        raise ValueError("a decay branch needs the stable run's steps")
    warmup = WARMUP if steps is None else round(WARMUP * steps)
    stable_end = branch_at if branch_at is not None else steps
    decay = 0 if branch_at is None else round(DECAY * branch_at)
    sched = {"num_decay_steps": decay, "decay_type": "linear", "min_lr_ratio": 0.0}  # stable = the rest
    workers = overrides.get("dataloader_num_workers", 6)
    return SFTConfig(**{
        "output_dir": output_dir, "bf16": _bf16(overrides), "assistant_only_loss": True, "max_length": None,
        "packing": False, "learning_rate": 5e-4 if (model in _LORA if lora is None else lora) else 5e-5,
        "lr_scheduler_type": "warmup_stable_decay", "warmup_steps": warmup, "lr_scheduler_kwargs": sched,
        **({"max_steps": stable_end + decay} if steps is not None else {}),
        "train_sampling_strategy": "group_by_length", "length_column_name": "length",
        "dataloader_num_workers": workers, "dataloader_persistent_workers": workers > 0,
        "dataloader_prefetch_factor": 2 if workers else None,
        "gradient_checkpointing": model in _GRAD_CKPT,
        "remove_unused_columns": False, "report_to": "none", "seed": 0, "data_seed": 0} | overrides)


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

