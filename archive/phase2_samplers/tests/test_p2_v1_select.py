"""Archived with experiments/phase2/sample_vllm.py (E4 step 3): tests of the script itself."""


def test_sample_vllm_turns_off_the_flashinfer_sampler_before_importing_vllm():
    """R1 2026-11-05: FlashInfer's sampler JIT-compiles with nvcc (absent from /opt/vllm); vLLM reads env at import."""
    import inspect

    src = inspect.getsource(V.main)
    assert V.ENV == {"VLLM_USE_FLASHINFER_SAMPLER": "0"}
    assert -1 < src.index("os.environ.update(ENV)") < src.index("import vllm")


