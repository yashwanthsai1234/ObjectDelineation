"""Phase C box tests on the pinned vLLM (Implementation.md C.9.4): G-GR1 and G-RC0 load the start checkpoint
`$OBJDEL_RUN_DIR/ckpt` (HF revision rl-full-3ep, fetched by fetch_hf.sh) and run by node id in their own plans
(c_speed_c6.yaml step 0; the last, soft step of b3_v5.yaml). A missing checkpoint fails, never skips."""

import os
from pathlib import Path

import pytest

from objdel import _gpucheck

CUDA = pytest.mark.skipif(_gpucheck.build() != "cuda", reason="CPU build: the CUDA build runs on the box")


@pytest.fixture(scope="module")
def engine():
    """vLLM on the start checkpoint as the C.9.4 model note loads it: text-only, eager, 30 % of the GPU."""
    from vllm import LLM

    ckpt = Path(os.environ["OBJDEL_RUN_DIR"]) / "ckpt"
    assert (ckpt / "config.json").is_file(), f"{ckpt}: the start checkpoint is missing (fetch_hf.sh <run> rl-full-3ep)"
    return LLM(model=str(ckpt), max_model_len=256, gpu_memory_utilization=0.3, enforce_eager=True)


@pytest.mark.gpu
@CUDA
def test_g_gr1_requests_independent(engine):
    """G-GR1 (C6): one request with n = 8 at temperature 1.0 gives 8 independent answers: over 32 prompts the
    average number of distinct answers per prompt is at least 2."""
    from vllm import SamplingParams

    prompts = [f"Item {i}: describe the outline of object number {i * 7 + 3} in one sentence." for i in range(32)]
    outs = engine.generate(prompts, SamplingParams(n=8, temperature=1.0, max_tokens=32, seed=None))
    distinct = [len({o.text for o in out.outputs}) for out in outs]
    assert all(len(out.outputs) == 8 for out in outs)
    assert sum(distinct) / len(distinct) >= 2, distinct


@pytest.mark.gpu
@CUDA
def test_g_rc0_per_row_caps(engine):
    """G-RC0 (C11 step 3): a SamplingParams list as long as the prompts gives each answer its own max_tokens; with
    ignore_eos every answer has exactly its cap, so a shared cap would show."""
    from vllm import SamplingParams

    caps = [4, 8, 16]
    params = [SamplingParams(temperature=1.0, max_tokens=c, ignore_eos=True) for c in caps]
    outs = engine.generate(["Outline of object A:", "Outline of object B:", "Outline of object C:"], params)
    assert [len(out.outputs[0].token_ids) for out in outs] == caps
