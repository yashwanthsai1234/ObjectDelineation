"""Per-model constants and helpers shared by SFT, RL and sampling: MODELS (Hugging Face ids), MAX_COMPLETION, the
LoRA config, the sampling fields, the stop ids and the truncation rule.

- Sampling: every field is explicit (T = 1, no filter, no penalty), overriding the model card's defaults, spelled for
  vLLM, GRPOConfig or HF generate by `sampling` and `as_names`.
- Completion budget: 1,536 tokens for Qwen3-VL (about 5.2 tokens per vertex). SFT drops targets that would not fit
  (objdel.sft.rows.sft_rows).
- LoRA: r16 on the language model's linear layers over a frozen bf16 base; the vision tower and projector are excluded
  (VISION).
- Truncation: TRL's rule, a completion is truncated when its last id is not a stop id (`stop_ids`, `truncated`).
- Images: every model keeps its processor's default image size (256 tokens at 512 px); left padding for generation.
- NO_SYNC: TrainingArguments that avoid a host sync per micro-batch.
"""

from peft import LoraConfig
from transformers.utils import is_torch_bf16_gpu_available

MODELS = {
    "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct",
}

MAX_COMPLETION = {"qwen3-vl-2b": 1536}

_LORA = frozenset(MODELS)  # every model is LoRA (ADR-0008)

_GRAD_CKPT = {"qwen3-vl-2b"}

# vision tower + projector of every model (recorded module names in tests/sft/test_p2_train.py)
VISION = r".*(visual|vision_model|vision_tower|connector|multi_modal_projector).*"


def lora(model: str) -> LoraConfig | None:
    """LoRA r16 on the language model's linear layers (vision tower and projector excluded) for the models
    trained with adapters (every model of MODELS), else None."""
    if model not in _LORA:
        return None
    return LoraConfig(r=16, lora_alpha=32, target_modules="all-linear", exclude_modules=VISION, task_type="CAUSAL_LM")


def tokenizer(model: str):
    """(tokenizer, max completion tokens) of a model, from the local HF cache: DG12's arguments of
    objdel.data.quality.audit."""
    from transformers import AutoTokenizer

    return AutoTokenizer.from_pretrained(MODELS[model]), MAX_COMPLETION[model]


def sampling(model: str, api: str = "vllm") -> dict:
    """Every sampling field explicit: T 1 and no filter or penalty, overriding the cards' defaults (Qwen3-VL
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
