"""TRL (version pinned in pyproject.toml) SFT and GRPO configurations for the screened models
(plan D8, D10, D12; lab t_stack).

GRPO (D8): loss_type="dr_grpo" (constant length normaliser, no std scaling), clip epsilon 0.2 / epsilon_high 0.28
(DAPO clip-higher), beta = 1e-3 from the start (lab t_stack T3), mask_truncated_completions, constant learning rate,
G = num_generations = 8 (REWARD_V4.group); any other num_generations or num_generations_eval raises ValueError,
because TRLReward would score a TRL group of 16 as two v4 groups of 8 without noticing (break-it B2-5). Reward:
TRLReward's weight-0 loggers plus the aggregator with multi_objective_aggregation="sum_then_normalize" and
scale_rewards="none" (objdel.rl.reward).
One step = 32 prompts x 8 completions (generation_batch_size 256, plan §3.2): per_device_train_batch_size 16 (two
groups) with 16 accumulation steps on one GPU (perf F2, 2026-11-06: half the CPU launch work per answer of 8 x 32;
the profiling run measures 8 / 16 / 32); per_device sizes are multiples of G so every process's reward call holds
whole groups at any world size. vLLM colocate (perf F1): vllm_gpu_memory_utilization 0.35, sleep mode on (vLLM
frees its memory during the update). num_generations_eval = None (= 8): held-out evaluation samples groups like
training, which is what the Phase 2 GRPO-signal gate (duplicates, sigma above the floor at T = 1) measures.

Owner decisions (C7, C8): beta stays 1e-3 in the control (B1) and is 0.0 for a B run that sets it through
`--grpo '{"beta": 0.0}'` (with beta 0 TRL loads no reference adapter and skips the reference pass; with LoRA the
reference is the adapter-off model, so beta can be raised again later). A completion limit L below the model's
budget goes in through `--grpo '{"max_completion_length": L}'`: train_rl refuses L above the budget, records the
dr_grpo loss scale budget / L in rl.json, and the soft overlong penalty (DAPO) needs `check_soft_overlong`. Truncated
rows stay masked by default; the runaway cap (objdel.rl.runaway) needs them unmasked.

max_completion_length: 1,536 tokens for Qwen3-VL (5.2 tokens per vertex, D10). Issue I5: on the final prepared
set 1 of 3,282 RIGA train targets is over budget and 0 of 23,469 LIDC train targets (prepare_rl's token_budget;
10 / 3,546 before BinRushed1 lost its model image, B3-1). The budget stays; such answers are dropped from SFT
(objdel.sft.rows.sft_rows with the model's tokenizer and MAX_COMPLETION), see experiments/phase1/impl/M7_train.md.

Generation: vLLM colocate for every model (each image placeholder expands to one contiguous run, lab t_stack T5);
a CPU caller passes use_vllm=False for HF generate.

Both configs (perf F3): logging_nan_inf_filter off (an isnan/isinf host sync on every micro-batch; a NaN loss still
shows at the log steps) and non-blocking host-to-GPU copies (pinned memory is on).

SFT: assistant-only loss (vision datasets need transformers >= 5.18; TRL swaps in its own training template for
the processor's chat template, qwen3_vl_training for Qwen3-VL), max_length=None (truncation can cut image tokens),
no packing. Every model is LoRA (ADR-0008): r16 on the language model's linear layers, lr 5e-4 (lab T4), over a
frozen bf16 base (full fine-tuning of Qwen3-VL-2B with AdamW needs ~34 GB for weights, gradients and optimizer state,
more than a 24 GB card has). GRPO learning rate: `rl_lr`, the default of `grpo_config` and the value
`rl.build.build_trainer` uses; `tests/experiments/test_prereg_lr.py` pins it to the pre-registered `prereg.G_PROBE`.
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

from peft import LoraConfig
from transformers.utils import is_torch_bf16_gpu_available

MODELS = {
    "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct",
}

MAX_COMPLETION = {"qwen3-vl-2b": 1536}

_LORA = frozenset(MODELS)  # every model is LoRA (ADR-0008)

_GRAD_CKPT = {"qwen3-vl-2b"}

# vision tower + projector of every model (recorded module names in tests/test_p2_train.py;
# archive/phase2_lab/B1/common.VISION)
VISION = r".*(visual|vision_model|vision_tower|connector|multi_modal_projector).*"


def lora(model: str) -> LoraConfig | None:
    """LoRA r16 on the language model's linear layers (vision tower and projector excluded, D6 / D7) for the models
    trained with adapters (every model of MODELS), else None."""
    if model not in _LORA:
        return None
    return LoraConfig(r=16, lora_alpha=32, target_modules="all-linear", exclude_modules=VISION, task_type="CAUSAL_LM")


def sampling(model: str, api: str = "vllm") -> dict:
    """Every sampling field explicit (D13): T 1 and no filter or penalty, overriding the cards' defaults (Qwen3-VL
    T 0.7 / top_p 0.8 / top_k 20). One source of values,
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


def compile_text_layers(model, mode: str) -> int:
    """Stage F arm: torch.compile each text decoder layer in place (regional compilation: one code object, BUCKETS
    shapes; parameter names unchanged). mode: "default" or "reduce-overhead" (CUDA Graphs). Returns the layer count."""
    import torch._dynamo.config as dc

    # one code object serves every layer; its variants are BUCKETS widths x grad on/off (update vs old/ref) x adapter
    # on/off (the ref pass): past recompile_limit (default 8) Dynamo runs eager with only a warning
    dc.recompile_limit = max(dc.recompile_limit, 64)
    dc.accumulated_recompile_limit = max(dc.accumulated_recompile_limit, 1024)
    layers = next(m for n, m in model.named_modules() if n.endswith("language_model.layers"))
    for layer in layers:
        layer.compile(mode=mode, dynamic=False)
    return len(layers)


def stop_ids(model, processing_class) -> set[int]:
    """The ids that end a completion for TRL: generation_config.eos_token_id, the tokenizer's EOS and pad."""
    tok = getattr(processing_class, "tokenizer", processing_class)
    eos = model.generation_config.eos_token_id
    return {*(eos if isinstance(eos, list) else [eos]), tok.eos_token_id, tok.pad_token_id} - {None}


def truncated(ids, stop) -> bool:
    """TRL's rule: a completion is truncated when its last id is not a stop id (an empty one has no end either)."""
    return not ids or ids[-1] not in stop


def _text(completion) -> str:
    """Conversational completions are [{"role": "assistant", "content": str | [parts]}]."""
    if isinstance(completion, list):
        content = completion[0]["content"]
        return content if isinstance(content, str) else "".join(p.get("text", "") for p in content)
    return completion
