"""B2.2 F step 5: the F check, run by `train_rl --f-check` inside the first run that uses `lora-seqmask` or
`fp16-seqmask` (one process, the run's own vLLM engine). The callback counts NaN/inf and the vLLM LoRA loads per RL
step over steps 1-15, runs the preset's probes, and writes `<out>/f_check.json`: the preset, each check's value, and
`passed`, the AND of every pass of that preset (`restart_preset` reads only `preset` and `passed`). On a failure it
raises after writing the file, so the plan step fails and the run restarts under B2.4 rule 0's fallback.

Fields under `fp16-seqmask`: nan, engine_dtype_ok. Under `lora-seqmask`: nan, comparator, base_hash_ok, lora_equal,
load_count, lora_over_base. The probes are callables (the engine reads); the comparator and the timing ratio are
measured before step 1 and passed in as values.
"""
import json
import math
from pathlib import Path

import numpy as np
from transformers import TrainerCallback

STEPS = 15  # the check covers optimizer steps 1-15
PRESETS = ("fp16-seqmask", "lora-seqmask")
HALF = 0.5  # comparator: the LoRA gap above the floor is at most this share of the merged gap
FLOOR = 0.002  # below this merged gap above the floor the comparator cannot tell F from no F
LORA_OVER_BASE_MAX = 1.15  # timing: median of the LoRA / base generate-time ratios


def comparator(f_merged: float, f_lora: float, f0: float) -> dict:
    """Output-level check: mean |logp gap| of the merged push, of the adapter push, and of the base at d = 0 (the
    floor). ratio = (f_lora - f0) / (f_merged - f0) passes at <= 0.5; a merged gap above the floor under 0.002 is not
    informative (ok, no ratio)."""
    gap = f_merged - f0
    out = {"f_merged": f_merged, "f_lora": f_lora, "f0": f0}
    if gap < FLOOR:
        return out | {"ratio": None, "informative": False, "ok": True}
    ratio = (f_lora - f0) / gap
    return out | {"ratio": ratio, "informative": True, "ok": ratio <= HALF}


def lora_over_base(t_lora: list[float], t_base: list[float]) -> float:
    """The timing value: the median of the t_A / t_B ratios (A = slot 2, B = `lora_request=None`), one per A B pair."""
    return float(np.median([a / b for a, b in zip(t_lora, t_base, strict=True)]))


def passes(rec: dict) -> bool:
    """The AND of every pass of the record's preset."""
    if rec["preset"] == "fp16-seqmask":
        return rec["nan"] == 0 and rec["engine_dtype_ok"]
    return (rec["nan"] == 0 and rec["comparator"]["ok"] and rec["base_hash_ok"] and rec["lora_equal"]
            and rec["load_count"] == 1 and rec["lora_over_base"] <= LORA_OVER_BASE_MAX)


def check_record(rec: dict) -> None:
    """Assert a written record: `passed` is the AND of the flags (never trusted alone) and it is True."""
    assert rec["passed"] == passes(rec), "passed differs from the AND of the checks"
    assert rec["passed"], {k: v for k, v in rec.items() if k != "passed"}


def count_calls(fn):
    """fn with a `calls` counter (the LoRA load count per RL step wraps `llm_engine.add_lora`)."""
    def wrapped(*a, **k):
        wrapped.calls += 1
        return fn(*a, **k)
    wrapped.calls = 0
    return wrapped


class FCheck(TrainerCallback):
    """Counts NaN/inf in the logs and the loads per step over steps 1-15, then writes `<out>/f_check.json`.
    probes: {name: callable} of the preset's checks (engine_dtype_ok; lora_equal, load_count's counter is `loads`);
    fixed: values measured before step 1 (comparator, base_hash_ok, lora_over_base)."""

    def __init__(self, out, preset: str, probes: dict, loads=None, fixed: dict | None = None):
        if preset not in PRESETS:
            raise ValueError(f"--f-check: no check for preset {preset!r}")
        self.out, self.preset, self.probes, self.loads, self.fixed = Path(out), preset, probes, loads, fixed or {}
        self.nan, self.counts, self.seen, self.done = 0, [], 0, False

    def on_log(self, args, state, control, logs=None, **kw):
        """NaN/inf count over the logged values of steps 1-15."""
        if state.global_step <= STEPS:
            self.nan += sum(isinstance(v, float) and not math.isfinite(v) for v in (logs or {}).values())

    def on_step_end(self, args, state, control, **kw):
        """The loads of this step; the record at step 15 or at the first failure."""
        if self.loads is not None:
            self.counts.append(self.loads.calls - self.seen)
            self.seen = self.loads.calls
        if not self.done and (state.global_step >= STEPS or self.nan or (self.counts and self.counts[-1] != 1)):
            self.finish()

    def finish(self) -> None:
        """Write f_check.json; raise when a check failed."""
        self.done = True
        rec = {"preset": self.preset, "nan": self.nan} | {k: p() for k, p in self.probes.items()} | self.fixed
        if self.loads is not None:
            rec["load_count"] = next((c for c in self.counts if c != 1), 1)
        rec["passed"] = bool(passes(rec))
        self.out.mkdir(parents=True, exist_ok=True)
        (self.out / "f_check.json").write_text(json.dumps(rec, indent=1, default=float) + "\n")
        if not rec["passed"]:
            raise RuntimeError(f"F check failed ({self.out / 'f_check.json'}): {rec}")
