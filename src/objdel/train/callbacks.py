"""TrainerCallbacks of the live runs: TraceStep (one torch.profiler step), NvtxSteps (NVTX step ranges and the
nsys capture window), BudgetStop (stop and save at a session spend) and OpsmDelta (the live OPSM delta).

The Phase 2 hack-audit callback (audit of held-out evaluation, best-checkpoint rollback, live guard switches) is
archived in archive/controller/ (CLEAN6): no live run attached it. The held-out curve is train_rl --heldout (TRL's own
evaluation, logged to rl.json).

References:
- transformers TrainerCallback / TrainerControl: https://huggingface.co/docs/transformers/main_classes/callback
"""
import json
import time
from collections import deque
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

    def on_train_begin(self, args, state, control, **kw):
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
            if state.global_step == self.window[0] - 1 and not self.capturing:
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
        if rate_usd_per_h <= 0:
            raise ValueError(f"budget stop: rate {rate_usd_per_h} $/h must be > 0")
        self.rate, self.cap, self.start, self.fixed = rate_usd_per_h, cap_usd, session_start, fixed_usd
        self.clock, self.tripped = clock, False

    def spend(self) -> float:
        return self.fixed + max(0.0, self.clock() - self.start) / 3600 * self.rate

    def on_step_end(self, args, state, control, **kw):
        self.tripped = self.tripped or self.spend() >= self.cap
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
        if state.global_step >= self.WINDOW and self.path.is_file():
            saved = json.loads(self.path.read_text())
            self.trainer.off_policy_mask_threshold, self.final = saved["delta"], True  # the ring is lost: no guard

    def on_step_end(self, args, state, control, **kw):
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
        delta = torch.quantile(torch.cat([s[0] for s in self.ring]), 0.995).item()
        self.trainer.off_policy_mask_threshold, self.final = delta, guarded
        self.path.write_text(json.dumps({"delta": delta, "guarded": guarded}) + "\n")
        print(f"OpsmDelta: delta = {delta:.6g} (guard: {guarded})")
