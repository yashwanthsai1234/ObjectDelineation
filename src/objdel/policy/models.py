"""The model registry and the helpers shared by SFT, RL and sampling: MODELS (one `Model` entry per model key: the
Hugging Face id, the completion budget, the LoRA and learning-rate settings, the processor keywords), the LoRA config,
the sampling fields, the stop ids and the truncation rule. Adding a model = one entry in MODELS; code reads it at call
time through `get`.

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

from dataclasses import dataclass

from peft import LoraConfig
from transformers.utils import is_torch_bf16_gpu_available

# vision tower + projector of every model (recorded module names in tests/sft/test_p2_train.py)
VISION = r".*(visual|vision_model|vision_tower|connector|multi_modal_projector).*"


@dataclass(frozen=True)
class Model:
    """A base model and the settings tuned for it. Every model is LoRA (ADR-0008)."""
    hf_id: str  # Hugging Face repo id of the base
    max_completion: int  # completion budget, tokens (module docstring)
    text_layers: int  # text decoder layers: the bound of --no-ckpt-layers (objdel.rl.run.checked)
    lora_r: int = 16
    lora_alpha: int = 32
    lora_exclude: str = VISION  # regex of the modules LoRA leaves out
    grad_ckpt: bool = True  # gradient checkpointing, SFT and GRPO (GRPO without it: OOM on the RTX 5090)
    sft_lr: float = 5e-4
    rl_lr: float = 5e-5  # the GRPO rate (prereg.G_PROBE["lr_lora"])
    processor: tuple[tuple[str, str], ...] = (("padding_side", "left"),)  # AutoProcessor keywords
    # GRPO batch on the reference card (RTX 5090, 2026-11-06, perf/p1_results.md: 16 x 16 44.3 s/step; 32 x 8 59.1):
    # groups per device (x the reward's group = the per-device batch) and gradient accumulation steps
    groups_per_device: int = 2
    grad_accum: int = 16
    vllm_train: tuple[float, int] = (0.35, 4096)  # colocated GRPO engine: GPU memory share, max model length
    vllm_eval: tuple[float, int] = (0.85, 2048)  # the eval sampler's engine (objdel.policy.sampler.VllmSampler)


MODELS: dict[str, Model] = {
    "qwen3-vl-2b": Model("Qwen/Qwen3-VL-2B-Instruct", max_completion=1536, text_layers=28),
}

# views of the library's entries, derived at import
MAX_COMPLETION = {k: m.max_completion for k, m in MODELS.items()}
_LORA = frozenset(MODELS)
_GRAD_CKPT = {k for k, m in MODELS.items() if m.grad_ckpt}


def get(model: str) -> Model:
    """The entry of a model key; KeyError naming the registered models."""
    if model not in MODELS:
        raise KeyError(f"unknown model {model!r}; registered: {sorted(MODELS)} (objdel.policy.models)")
    return MODELS[model]


def lora(model: str) -> LoraConfig:
    """LoRA on the language model's linear layers (vision tower and projector excluded), rank and alpha of the
    model's entry."""
    m = get(model)
    return LoraConfig(r=m.lora_r, lora_alpha=m.lora_alpha, target_modules="all-linear", exclude_modules=m.lora_exclude,
                      task_type="CAUSAL_LM")


def tokenizer(model: str):
    """(tokenizer, max completion tokens) of a model, from the local HF cache: DG12's arguments of
    objdel.data.quality.audit."""
    from transformers import AutoTokenizer

    m = get(model)
    return AutoTokenizer.from_pretrained(m.hf_id), m.max_completion


def sampling(model: str, api: str = "vllm") -> dict:
    """Every sampling field explicit: T 1 and no filter or penalty, overriding the cards' defaults (Qwen3-VL
    T 0.7 / top_p 0.8 / top_k 20). One source of values,
    spelled for api "vllm" (SamplingParams), "grpo" (GRPOConfig) or "hf" (GenerationConfig), see as_names."""
    s = {"temperature": 1.0, "top_p": 1.0, "top_k": 0, "min_p": 0.0, "repetition_penalty": 1.0,
         "presence_penalty": 0.0, "frequency_penalty": 0.0, "max_tokens": get(model).max_completion}
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


def processor_kwargs(model: str) -> dict:
    """AutoProcessor.from_pretrained(path, **processor_kwargs(model)): the entry's processor keywords."""
    return dict(get(model).processor)


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
