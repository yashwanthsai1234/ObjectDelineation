"""TrainerCallbacks of the live runs: TraceStep (one torch.profiler step), NvtxSteps (NVTX step ranges and the
nsys capture window) and BudgetStop (stop and save at a session spend).

The Phase 2 hack-audit callback (audit of held-out evaluation, best-checkpoint rollback, live guard switches) is
archived in archive/controller/ (CLEAN6): no live run attached it. The held-out curve is train_rl --heldout (TRL's own
evaluation, logged to rl.json).

References:
- transformers TrainerCallback / TrainerControl: https://huggingface.co/docs/transformers/main_classes/callback
"""
import math
import time
from pathlib import Path

import torch
from transformers import TrainerCallback

from objdel.train import nvtx


class TraceStep(TrainerCallback):
    """torch.profiler on one optimizer step after `wait` steps and one profiler warm-up step (CPU ops, and CUDA
    kernels on a GPU): <out>/trace.json.gz and ops.txt (top ops by device and CPU time). train_sft `profile:`,
    g_probe --trace."""

    def __init__(self, out: str, wait: int):
        self.out, self.wait = Path(out), wait

    def on_train_begin(self, args, state, control, **kw):
        from torch.profiler import ProfilerActivity, profile, schedule

        acts = [ProfilerActivity.CPU] + [ProfilerActivity.CUDA] * torch.cuda.is_available()
        self.prof = profile(activities=acts, schedule=schedule(wait=self.wait, warmup=1, active=1),
                            on_trace_ready=self.write)  # fmt: skip
        self.prof.start()

    def on_step_end(self, args, state, control, **kw):
        self.prof.step()

    def on_train_end(self, args, state, control, **kw):
        self.prof.stop()

    def write(self, prof):
        self.out.mkdir(parents=True, exist_ok=True)
        prof.export_chrome_trace(str(self.out / "trace.json.gz"))
        ops = prof.key_averages()
        (self.out / "ops.txt").write_text("\n\n".join(ops.table(sort_by=k, row_limit=40) for k in
                                                       ("self_device_time_total", "self_cpu_time_total")))


class NvtxSteps(TrainerCallback):
    """NVTX range `step` over each step period (one on_step_end to the next: the previous step's log, this step's
    loader, micro-batches, grad clip and optimizer) and `optimizer` from on_pre_optimizer_step to on_optimizer_step
    (objdel.train.nvtx: no-ops unless OBJDEL_NVTX=1 and CUDA). With window = (s, e) and CUDA, the cudaProfilerApi
    capture (nsys --capture-range=cudaProfilerApi, nsys_capture.sh) holds exactly the step ranges of steps s..e: it
    starts at the end of step s - 1, between the pop and the next push, and stops at the end of step e after the pop
    (or at on_train_end, when the run ends inside the window). train_rl --profile-window S E."""

    def __init__(self, window: tuple[int, int] | None):
        if window is not None and not 2 <= window[0] <= window[1]:
            raise ValueError(f"profile window {window}: needs 2 <= start <= end (the capture starts at the end of "
                             "step start - 1, and TraceStep's wait is start - 2)")
        self.window, self.open, self.opt_open, self.capturing, self.done = window, False, False, False, False

    def _in_window(self, step: int) -> bool:  # a range, not == s - 1: a run resumed inside the window still captures
        return bool(self.window) and not self.done and self.window[0] - 1 <= step < self.window[1]

    def on_train_begin(self, args, state, control, **kw):
        if self._in_window(state.global_step) and torch.cuda.is_available():
            torch.cuda.profiler.start()
            self.capturing = True
        nvtx.push("step")
        self.open = True

    def on_pre_optimizer_step(self, args, state, control, **kw):
        nvtx.push("optimizer")
        self.opt_open = True

    def on_optimizer_step(self, args, state, control, **kw):
        if self.opt_open:
            nvtx.pop()
            self.opt_open = False

    def on_step_end(self, args, state, control, **kw):
        if self.open:
            nvtx.pop()
        if self.window and not self.done and torch.cuda.is_available():
            if self._in_window(state.global_step) and not self.capturing:
                torch.cuda.profiler.start()
                self.capturing = True
            elif state.global_step == self.window[1] and self.capturing:
                self._stop()
        nvtx.push("step")
        self.open = True

    def on_train_end(self, args, state, control, **kw):
        if self.open:
            nvtx.pop()
            self.open = False
        if self.capturing:
            self._stop()

    def _stop(self):
        torch.cuda.profiler.stop()
        self.capturing, self.done = False, True


class BudgetStop(TrainerCallback):
    """Stop and save once the session spend fixed_usd + (now - session_start) h x rate_usd_per_h reaches cap_usd
    (train_rl --budget-usd). With run.sh's exports (rate = dph + 6 x inet_up_cost, fixed = 15 x inet_down_cost,
    session_start = the launch time) this is the amount vast.policy.accrued charges, so the guard trips before the
    rental's own stop rules and the evaluation steps after train_rl still run. Once tripped it stays tripped."""

    def __init__(self, rate_usd_per_h: float, cap_usd: float, session_start: float, fixed_usd: float = 0.0,
                 clock=time.time):
        if not all(map(math.isfinite, (rate_usd_per_h, cap_usd, session_start, fixed_usd))):
            raise ValueError(f"budget stop: rate, cap, start and fixed must be finite numbers, got {rate_usd_per_h}, "
                             f"{cap_usd}, {session_start}, {fixed_usd}")
        if rate_usd_per_h <= 0 or cap_usd < 0:
            raise ValueError(f"budget stop: rate {rate_usd_per_h} $/h must be > 0 and cap ${cap_usd} >= 0")
        if session_start > clock() + 60:
            raise ValueError(f"budget stop: session start {session_start} is in the future (ms? wrong clock?)")
        self.rate, self.cap, self.start, self.fixed = rate_usd_per_h, cap_usd, session_start, fixed_usd
        self.clock, self.tripped = clock, False

    def spend(self) -> float:
        return self.fixed + max(0.0, self.clock() - self.start) / 3600 * self.rate

    def on_step_end(self, args, state, control, **kw):
        self.tripped = self.tripped or not self.spend() < self.cap
        if self.tripped:
            control.should_training_stop = control.should_save = True
