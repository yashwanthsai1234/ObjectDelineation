"""Phase C box tests (plan C.9.4). G-RC0 and G-GR1 load the start checkpoint `$OBJDEL_RUN_DIR/ckpt` (HF revision
`rl-full-3ep`, fetched by the plan's item 0) on the pinned vLLM and run by node id in their own plan steps; U-RC6 checks
on the Mac that the same model's config, tokenizer and processor load (marker `data`: network and an HF token).
`gpu` tests skip on the CPU build; a missing checkpoint fails them, never skips."""

import os
from pathlib import Path

import pytest

from experiments.phase2.hf_push import REPO
from objdel.geom import engine as _gpucheck

CUDA = pytest.mark.skipif(_gpucheck.build() != "cuda", reason="CPU build: the CUDA build runs on the box")
REVISION = "rl-full-3ep"
FILES = ("config.json", "generation_config.json", "tokenizer.json", "tokenizer_config.json", "processor_config.json",
         "chat_template.jinja")  # what U-RC6 downloads: no weights
PROMPT = "Describe a nodule in one sentence."


def check_model_files(path):
    """U-RC6: the folder's config, tokenizer and processor load, and the chat template renders one text prompt."""
    from transformers import AutoConfig, AutoProcessor, AutoTokenizer

    cfg = AutoConfig.from_pretrained(path)
    assert cfg.architectures == ["Qwen3VLForConditionalGeneration"]
    assert cfg.text_config.max_position_embeddings >= 256  # the engine's max_model_len
    tokenizer = AutoTokenizer.from_pretrained(path)
    assert tokenizer("Outline of object A:")["input_ids"] and tokenizer.pad_token_id is not None
    processor = AutoProcessor.from_pretrained(path)
    assert getattr(processor, "image_processor", None) is not None and processor.image_processor.patch_size == 16
    chat = [{"role": "user", "content": [{"type": "text", "text": PROMPT}]}]
    assert PROMPT in processor.apply_chat_template(chat, tokenize=False, add_generation_prompt=True)


@pytest.mark.data
def test_vllm_check_model_loads(tmp_path):
    """U-RC6: the model of G-RC0 and G-GR1 (HF revision rl-full-3ep) loads on the Mac from its small files only."""
    from huggingface_hub import get_token, snapshot_download

    token = get_token()
    if token is None:
        pytest.skip("no HF token set: this check needs the owner's token for the private model repo")
    snapshot_download(REPO, revision=REVISION, allow_patterns=list(FILES), local_dir=tmp_path, token=token)
    check_model_files(tmp_path)


def engine():
    """The pinned vLLM on the start checkpoint, text-only prompts (plan C.9.4 model)."""
    from vllm import LLM

    ckpt = Path(os.environ["OBJDEL_RUN_DIR"]) / "ckpt"
    assert (ckpt / "config.json").is_file(), f"{ckpt}: no start checkpoint (fetch_hf.sh <run> {REVISION} puts it there)"
    return LLM(model=str(ckpt), max_model_len=256, gpu_memory_utilization=0.3, enforce_eager=True)


@pytest.mark.gpu
@CUDA
def test_g_rc0_per_row_caps():
    """G-RC0: a SamplingParams list gives each answer its own max_tokens (ignore_eos: every answer reaches its cap,
    so a shared cap would show)."""
    from vllm import SamplingParams

    caps = [4, 8, 16]
    params = [SamplingParams(max_tokens=c, ignore_eos=True, temperature=1.0) for c in caps]
    out = engine().generate([PROMPT] * len(caps), sampling_params=params)
    assert [len(o.outputs[0].token_ids) for o in out] == caps


@pytest.mark.gpu
@CUDA
def test_g_gr1_requests_independent():
    """G-GR1: one request with n = 8 per prompt gives independent samples: at temperature 1.0 a prompt's 8 answers
    average at least 2 distinct texts over 32 prompts."""
    from vllm import SamplingParams

    prompts = [f"{PROMPT} Case {i}." for i in range(32)]
    out = engine().generate(prompts, sampling_params=SamplingParams(n=8, temperature=1.0, max_tokens=32))
    assert all(len(o.outputs) == 8 for o in out)
    assert sum(len({c.text for c in o.outputs}) for o in out) / len(out) >= 2
