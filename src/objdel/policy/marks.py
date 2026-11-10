"""NVTX ranges for an Nsight Systems capture (plan A6), off unless OBJDEL_NVTX=1 and CUDA is present.

Fork safe: the import reads only the environment. It imports no torch and asks no CUDA question, because
`torch.cuda.is_available()` runs cuInit (unless PYTORCH_NVML_BASED_CUDA_CHECK=1), which poisons a later fork. The
CUDA query runs at the first `enabled()` call (the first range) and is cached for the process.

OBJDEL_NVTX_SYNC (read once, here): "1" (default) keeps FastGRPO.timed's end-of-phase device sync, as production
runs; "0" times each phase with a CUDA event pair instead, so a trace shows the real overlap between phases.

Range names are fixed (tests/experiments/test_nvtx.py scans for them): step, loader, micro, batch, generate, sync, wake,
wake_weights, wake_kv, sleep, rewards, checker, logps, forward, optimizer, log, vision_cache, rope_cache,
rope_template, image_memo, encoder_once.
Plus `generate/tail` (C11, see `tail`), opened inside `generate`.
"""

import contextlib
import os

_ON = os.environ.get("OBJDEL_NVTX") == "1"
SYNC = os.environ.get("OBJDEL_NVTX_SYNC", "1") != "0"
_cuda: bool | None = None  # torch.cuda.is_available(), asked once


def enabled() -> bool:
    """True when NVTX is on (`OBJDEL_NVTX=1`) and CUDA is available; the CUDA query is asked once."""
    global _cuda
    if not _ON:
        return False
    if _cuda is None:
        import torch

        _cuda = torch.cuda.is_available()
    return _cuda


def push(name: str) -> None:
    """Push an NVTX range when enabled."""
    if enabled():
        import torch

        torch.cuda.nvtx.range_push(name)


def pop() -> None:
    """Pop the current NVTX range when enabled."""
    if enabled():
        import torch

        torch.cuda.nvtx.range_pop()


@contextlib.contextmanager
def _range(name):
    push(name)
    try:
        yield
    finally:
        pop()


def range(name: str):  # noqa: A001 (the plan's name; torch.cuda.nvtx has `range` too)
    """A context manager: range_push(name) on entry and range_pop() on exit (also on an exception) when enabled,
    else a no-op."""
    return _range(name) if enabled() else contextlib.nullcontext()


@contextlib.contextmanager
def tail(engine, below: int):
    """The decode tail (C11): a `generate/tail` range from the first `engine.step()` that leaves fewer than `below`
    unfinished requests to the end of the block, so a trace shows the seconds spent on the few longest answers (C0
    step 6b checks its logged-length estimate against it). The patched `step` lives on the instance only inside the
    block. A no-op when NVTX is off. Needs an engine with `step()` and `get_num_unfinished_requests()`, which
    `LLM._run_engine` loops on [unverified for vLLM 0.32.0: the box]."""
    if not enabled():
        yield
        return
    step, thin = engine.step, False

    def counted(*args, **kwargs):
        nonlocal thin
        out = step(*args, **kwargs)
        if not thin and engine.get_num_unfinished_requests() < below:
            push("generate/tail")
            thin = True
        return out

    engine.step = counted
    try:
        yield
    finally:
        del engine.step  # the instance patch only; the class's step is back
        if thin:
            pop()


class SyncCount:
    """Result holder of count_syncs: the number of synchronizing CUDA operations seen in the block."""
    syncs = 0  # synchronizing CUDA operations torch reported inside the block


@contextlib.contextmanager
def count_syncs():
    """Counts torch's synchronizing CUDA operations inside the block (C0 step 1): torch.cuda.set_sync_debug_mode("warn")
    and every warning recorded (simplefilter "always": Python would record one code line's warning once), counting
    the messages that say "synchroniz"; the previous mode comes back on exit, also on an exception. A nested block
    keeps its own count. Without CUDA a no-op (0). torch's mode is experimental: not every sync warns."""
    import warnings

    import torch

    out = SyncCount()
    if not torch.cuda.is_available():
        yield out
        return
    prev = torch.cuda.get_sync_debug_mode()
    with warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter("always")
        torch.cuda.set_sync_debug_mode("warn")
        try:
            yield out
        finally:
            torch.cuda.set_sync_debug_mode(prev)
            out.syncs = sum("synchroniz" in str(w.message) for w in seen)
