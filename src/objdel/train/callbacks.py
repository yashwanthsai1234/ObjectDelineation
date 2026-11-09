"""TrainerCallbacks of the live runs: TraceStep (one torch.profiler step), NvtxSteps (NVTX step ranges and the
nsys capture window), BudgetStop (stop and save at a session spend), OpsmDelta (the live OPSM delta),
NoCkptLayers (C9: no checkpointing on K layers) and LocGuard (stop and save when reward v5's location term
rises while IoU_tc falls on the off-centre rows; Implementation.md B3.5b).

The Phase 2 hack-audit callback (audit of held-out evaluation, best-checkpoint rollback, live guard switches) is
archived in archive/controller/ (CLEAN6): no live run attached it. The held-out curve is train_rl --heldout (TRL's own
evaluation, logged to rl.json).

References:
- transformers TrainerCallback / TrainerControl: https://huggingface.co/docs/transformers/main_classes/callback
- Newey-West HAC standard error (Bartlett weights): Newey & West 1987, https://doi.org/10.2307/1913610
"""
import json
import math
import time
from collections import deque
from pathlib import Path

import numpy as np
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
        """Start the profiler."""
        from torch.profiler import ProfilerActivity, profile, schedule

        acts = [ProfilerActivity.CPU] + [ProfilerActivity.CUDA] * torch.cuda.is_available()
        self.prof = profile(activities=acts, schedule=schedule(wait=self.wait, warmup=1, active=1),
                            on_trace_ready=self.write)  # fmt: skip
        self.prof.start()

    def on_step_end(self, args, state, control, **kw):
        """Advance the profiler one step."""
        self.prof.step()

    def on_train_end(self, args, state, control, **kw):
        """Stop the profiler."""
        self.prof.stop()

    def write(self, prof):
        """Write the Chrome trace and the top-op tables to `out` when the profiler trace is ready."""
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
        """Open the `step` range."""
        if self._in_window(state.global_step) and torch.cuda.is_available():
            torch.cuda.profiler.start()
            self.capturing = True
        nvtx.push("step")
        self.open = True

    def on_pre_optimizer_step(self, args, state, control, **kw):
        """Open the `optimizer` range."""
        nvtx.push("optimizer")
        self.opt_open = True

    def on_optimizer_step(self, args, state, control, **kw):
        """Close the `optimizer` range."""
        if self.opt_open:
            nvtx.pop()
            self.opt_open = False

    def on_step_end(self, args, state, control, **kw):
        """Close the `step` range, start or stop the CUDA capture at the window ends, and open the next `step` range."""
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
        """Close any open range and stop a running capture."""
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
        """Session spend so far in USD."""
        return self.fixed + max(0.0, self.clock() - self.start) / 3600 * self.rate

    def on_step_end(self, args, state, control, **kw):
        """After each step, trip the guard once spend reaches the cap; once tripped, request a stop with a save."""
        self.tripped = self.tripped or not self.spend() < self.cap
        if self.tripped:
            control.should_training_stop = control.should_save = True


class OpsmDelta(TrainerCallback):
    """B2.4 rule 4: the OPSM delta (trainer.off_policy_mask_threshold, read by TRL at each loss call) is the preset's
    start value for steps 1-20, then q99.5 of the 20 x 256 per-sequence seq_kl values of those steps, written to
    <out>/opsm_delta.json. Guard, once: if the share of sequences OPSM masks (advantage < 0 and seq_kl > delta)
    stays above 2 % for 20 consecutive steps, delta becomes q99.5 of the last 20 steps. Then it never changes. A
    resume reads the file. FastGRPO.set_is_ratio calls `add` in train mode (trainer.opsm_watch)."""

    WINDOW, SHARE = 20, 0.02

    def __init__(self, trainer, out_dir: str):
        self.trainer, self.path = trainer, Path(out_dir) / "opsm_delta.json"
        self.cur, self.ring, self.streak, self.final = [], deque(maxlen=self.WINDOW), 0, False

    def add(self, kl: torch.Tensor, adv: torch.Tensor):
        """One micro-batch's per-sequence seq_kl and advantage, on the device."""
        self.cur.append(torch.stack([kl.float(), adv.float()]))

    def on_train_begin(self, args, state, control, **kw):
        """On a resume past the window, take the delta saved in opsm_delta.json (final: no guard)."""
        if state.global_step >= self.WINDOW and self.path.is_file():
            saved = json.loads(self.path.read_text())
            self.trainer.off_policy_mask_threshold, self.final = saved["delta"], True  # the ring is lost: no guard

    def on_step_end(self, args, state, control, **kw):
        """Close the step's rows; set the delta at step 20, then watch the masked share (the guard)."""
        if not self.cur:
            return
        step, self.cur = torch.cat(self.cur, 1), []
        if self.final:
            return
        self.ring.append(step)
        if state.global_step == self.WINDOW:
            self.set_delta(guarded=False)
        elif state.global_step > self.WINDOW:
            masked = (step[1] < 0) & (step[0] > self.trainer.off_policy_mask_threshold)
            self.streak = self.streak + 1 if masked.float().mean().item() > self.SHARE else 0
            if self.streak >= self.WINDOW:
                self.set_delta(guarded=True)

    def set_delta(self, guarded: bool):
        """Delta = q99.5 of the ring's seq_kl; write opsm_delta.json."""
        delta = torch.quantile(torch.cat([s[0] for s in self.ring]), 0.995).item()
        self.trainer.off_policy_mask_threshold, self.final = delta, guarded
        self.path.write_text(json.dumps({"delta": delta, "guarded": guarded}) + "\n")
        print(f"OpsmDelta: delta = {delta:.6g} (guard: {guarded})")


class NoCkptLayers(TrainerCallback):
    """C9 arm (a), train_rl --no-ckpt-layers K: gradient checkpointing off on the last K text decoder layers (their
    activations kept: less recompute, more memory). At train begin, because Trainer.train re-enables checkpointing
    on every layer just before it."""

    def __init__(self, k: int):
        self.k = k

    def on_train_begin(self, args, state, control, model=None, **kw):
        """Turn checkpointing off on the last K decoder layers; ValueError if K or a layer is unfit."""
        from transformers.modeling_layers import GradientCheckpointingLayer

        layers = next(m for n, m in model.named_modules() if n.endswith("language_model.layers"))
        if not 0 < self.k <= len(layers) or not all(isinstance(x, GradientCheckpointingLayer) for x in layers):
            raise ValueError(f"--no-ckpt-layers {self.k}: needs 1..{len(layers)} checkpointable decoder layers")
        for layer in list(layers)[-self.k:]:
            layer.gradient_checkpointing = False


def slope_nw(y, lag: int = 5) -> tuple[float, float]:
    """(least-squares slope of y over its index, its Newey-West standard error with Bartlett weights up to lag)."""
    y = np.asarray(y, np.float64)
    x = np.arange(y.size) - (y.size - 1) / 2
    sxx = float(x @ x)
    b = float(x @ (y - y.mean())) / sxx
    u = x * (y - y.mean() - b * x)
    s = float(u @ u) + 2 * sum((1 - k / (lag + 1)) * float(u[k:] @ u[:-k]) for k in range(1, lag + 1))
    return b, float(np.sqrt(max(s, 0.0))) / sxx


class LocGuard(TrainerCallback):
    """Reward v5 run guard (train_rl --reward v5): the model learns to move outlines toward the nodule without
    hitting it. Per logged training step it takes the mean over the off-centre kinds (shifted, rescaled) of
    kind/<k>/loc and kind/<k>/iou_tc. At the end of every block of `window` such steps it tests the block: (a) the
    location slope has t > 3, (b) the upper one-sided 95 % bound of the IoU_tc slope is < 0, both with Newey-West
    standard errors (lag 5; per-step means are autocorrelated). It stops and saves when two consecutive blocks
    pass, so never before 2 x window steps. Once tripped it stays tripped."""

    KINDS = ("shifted", "rescaled")

    def __init__(self, window: int = 50, lag: int = 5):
        self.window, self.lag = window, lag
        self.loc: list[float] = []
        self.iou: list[float] = []
        self.passed: list[bool] = []
        self.tripped = False

    def block(self, loc, iou) -> bool:
        """True when one block shows outlines moving off-centre (loc slope t > 3) while IoU falls (95 % one-sided)."""
        b, se = slope_nw(loc, self.lag)
        c, se_c = slope_nw(iou, self.lag)
        return bool(se > 0 and b / se > 3 and c + 1.645 * se_c < 0)

    def add(self, loc: float, iou: float) -> bool:
        """Record one step's off-centre means; True once the guard has tripped."""
        self.loc.append(loc)
        self.iou.append(iou)
        if len(self.loc) % self.window == 0:
            w = slice(-self.window, None)
            self.passed.append(self.block(self.loc[w], self.iou[w]))
            self.tripped = self.tripped or self.passed[-2:] == [True, True]
        return self.tripped

    def on_log(self, args, state, control, logs=None, **kw):
        """Feed the logged off-centre means to the guard; request a stop with a save once it trips."""
        keys = [k for k in self.KINDS if f"kind/{k}/loc" in (logs or {})]
        if keys and self.add(float(np.mean([logs[f"kind/{k}/loc"] for k in keys])),
                             float(np.mean([logs[f"kind/{k}/iou_tc"] for k in keys]))):
            control.should_training_stop = control.should_save = True
