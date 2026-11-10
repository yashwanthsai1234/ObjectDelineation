"""NVTX ranges for an Nsight Systems capture (plan A6), off unless OBJDEL_NVTX=1 and CUDA is present.

Fork safe: the import reads only the environment. It imports no torch and asks no CUDA question, because
`torch.cuda.is_available()` runs cuInit (unless PYTORCH_NVML_BASED_CUDA_CHECK=1), which poisons a later fork, and
serve.fork_pool's guard sees only `torch.cuda.is_initialized()`. The CUDA query runs at the first `enabled()` call
(the first range, after the fork) and is cached for the process.

OBJDEL_NVTX_SYNC (read once, here): "1" (default) keeps FastGRPO.timed's end-of-phase device sync, as production
runs; "0" times each phase with a CUDA event pair instead, so a trace shows the real overlap between phases.

Range names are fixed (tests/test_nvtx.py scans for them): step, loader, micro, batch, generate, sync, wake, sleep,
rewards, checker, logps, forward, optimizer, log, vision_cache, rope_cache, rope_template (C2), image_memo.
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
