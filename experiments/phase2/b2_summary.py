"""B2.4 decision (ADR-0013 rules 0-4, $0) and the B2.3 fp16-arm summary.

    python -I -m experiments.phase2.b2_summary decide --b1 RD/b1-rows --merge-rounding RD/merge_rounding.json
        [--fp16 RD/b2-fp16] --out RD/decision.json
    python -I -m experiments.phase2.b2_summary summary --fp16 RD/b2-fp16 --b1 RD/b1-rows

The functions are pure over the B1 run's records: `<run>/rl.json` (the `cf/*` stats of B2.1 per step),
`<run>/trainer/seq_gaps.f32` (steps x rows x 3 float32: seq_kl, advantage, valid_tokens) and B2.0's table.
Preset names are configs.NUMERICS keys. `start_preset` is the choice before the run, `restart_preset` the choice
after the F check (`f_check.json`) failed; both return a Decision whose start_delta is the number for
`train_rl --opsm-start-delta` (None for every preset but bf16-tis-opsm).
"""
import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

REL_ERR_LORA = 0.2  # rule 0: B2.0 bf16 relative error at step 100, at or above it -> lora-seqmask
FP16_HALF = 0.5  # rule 1: fp16 logp difference at most this share of B1's
OPSM_VETO = 0.02  # rule 3: share of sequences masked at delta 0.5 above which the start delta is raised
CLAMP_NOTE = 0.01  # rule 3: share of tokens clamped at C = 2 above which the share is recorded
ESS_VETO = 0.5  # rule 3: median ESS at delta 0.1 below which bf16-tis-opsm falls back to bf16-seqmask
START_DELTA = 0.5  # TRL's docstring example value of off_policy_mask_threshold
CF_STEPS = (50, 260)  # B1 steps whose counterfactual stats rules 2-3 read (the adapter has drifted)
CF_KEYS = ("cf/opsm_share@0.5", "cf/tok_clamp2_share", "cf/ess_E@0.1")
FP16_STEPS = (41, 60)


@dataclass(frozen=True)
class Decision:
    preset: str
    start_delta: float | None
    notes: tuple[str, ...] = ()


def fp16_passes(fp16: dict | None, b1: dict) -> bool:
    """Rule 1: the fp16 arm passed Gate B2's fp16 rules and its logp difference over steps 41-60 is at most half of
    B1's. Speed is a report line only (two rentals, one dtype per vLLM engine)."""
    return bool(fp16 and fp16["passed_gate"] and fp16["logp_diff_41_60"] <= FP16_HALF * b1["logp_diff_41_60"])


def rules_2_3(cf: dict, b1: dict) -> Decision:
    """bf16-tis-opsm with its vetoes over B1 steps 50-260; bf16-seqmask when the ESS veto fires."""
    if cf["cf/ess_E@0.1"] < ESS_VETO:
        return Decision("bf16-seqmask", None, (f"cf/ess_E@0.1 {cf['cf/ess_E@0.1']} < {ESS_VETO}: bf16-seqmask",))
    notes = []
    delta = START_DELTA
    if cf["cf/opsm_share@0.5"] > OPSM_VETO:
        delta = float(np.quantile(b1["seq_kl_50_260"], 0.999))
        notes.append(f"cf/opsm_share@0.5 {cf['cf/opsm_share@0.5']} > {OPSM_VETO}: start delta = q99.9 seq_kl = {delta}")
    if cf["cf/tok_clamp2_share"] > CLAMP_NOTE:
        notes.append(f"cf/tok_clamp2_share {cf['cf/tok_clamp2_share']} > {CLAMP_NOTE}: recorded, preset kept")
    return Decision("bf16-tis-opsm", delta, tuple(notes))


def start_preset(cf: dict, rel_err_100: float, fp16: dict | None, b1: dict) -> Decision:
    """The preset the first run after Gate B3 starts under. Reads the B2.0 error, never an F result."""
    if rel_err_100 >= REL_ERR_LORA:
        return Decision("lora-seqmask", None, (f"B2.0 relative error {rel_err_100} >= {REL_ERR_LORA}",))
    if fp16_passes(fp16, b1):
        return Decision("fp16-seqmask", None)
    return rules_2_3(cf, b1)


def restart_preset(cf: dict, f_check: dict, fp16: dict | None, b1: dict) -> Decision:
    """The preset a failed F check restarts under; reads only f_check["preset"] and f_check["passed"]. A restart
    never returns the preset that just failed."""
    if f_check["passed"]:
        return Decision(f_check["preset"], None)
    if f_check["preset"] == "lora-seqmask" and fp16_passes(fp16, b1):
        return Decision("fp16-seqmask", None)
    return rules_2_3(cf, b1)


def opsm_counterfactual(gaps: np.ndarray) -> tuple[float, float]:
    """Rule 4 priced on B1: delta = q99.5 of seq_kl over the first 20 steps; the share of steps 50-260 sequences
    with advantage < 0 and seq_kl > delta (the ones OPSM would mask)."""
    delta = float(np.quantile(gaps[:20, :, 0], 0.995))
    window = gaps[CF_STEPS[0] - 1:CF_STEPS[1]]
    return delta, float(((window[:, :, 1] < 0) & (window[:, :, 0] > delta)).mean())


def load_seq_gaps(trainer_dir: Path) -> np.ndarray:
    """(steps, rows, 3) float32 from seq_gaps.f32 and the row count in seq_gaps.json."""
    d = Path(trainer_dir)
    rows = json.loads((d / "seq_gaps.json").read_text())["rows_per_step"]
    flat = np.fromfile(d / "seq_gaps.f32", dtype="<f4")
    if flat.size % (rows * 3):
        raise ValueError(f"{d}/seq_gaps.f32: {flat.size} values is not whole steps of {rows} x 3")
    return flat.reshape(-1, rows, 3)


def _window(log: list[dict], lo: int, hi: int) -> list[dict]:
    got = {h["step"]: h for h in log if "step" in h}
    missing = [s for s in (lo, hi) if s not in got]
    if missing:
        raise ValueError(f"the log has no entry for step {missing[0]} (steps {lo}-{hi} are needed)")
    return [got[s] for s in range(lo, hi + 1) if s in got]


def cf_medians(log: list[dict]) -> dict:
    """Median over B1 steps 50-260 of each cf/* stat rules 2-3 read; ValueError names a missing stat or step."""
    rows = _window(log, *CF_STEPS)
    for k in CF_KEYS:
        if any(k not in h for h in rows):
            raise ValueError(f"the B1 log lacks {k} in steps {CF_STEPS[0]}-{CF_STEPS[1]}")
    return {k: float(np.median([h[k] for h in rows])) for k in CF_KEYS}


def b1_inputs(rl: dict, gaps: np.ndarray) -> dict:
    """The `b1` argument of the decision functions."""
    rows = _window(rl["log"], *FP16_STEPS)
    return {"logp_diff_41_60": _mean(rows, "sampling/sampling_logp_difference/mean"),
            "s_per_step": rl["timing_s"]["per_step"],
            "seq_kl_50_260": gaps[CF_STEPS[0] - 1:CF_STEPS[1], :, 0].ravel()}  # fmt: skip


def rel_err_100(path: Path) -> float:
    """B2.0 table: the step-100 bf16 relative error."""
    for row in json.loads(Path(path).read_text())["rows"]:
        if row["step"] == 100:
            return row["bf16"]["rel_err"]
    raise ValueError(f"{path}: no step 100 row")


def _mean(rows: list[dict], key: str) -> float:
    return float(np.mean([h[key] for h in rows]))


def decide(b1_run: Path, merge_rounding: Path, fp16: dict | None = None) -> dict:
    """Both presets (ADR-0013), computed before the run: the start preset, and the one a failed F check restarts
    under (the F check can only fail under lora-seqmask or fp16-seqmask)."""
    rl = json.loads((Path(b1_run) / "rl.json").read_text())
    gaps = load_seq_gaps(Path(b1_run) / "trainer")
    cf, b1, err = cf_medians(rl["log"]), b1_inputs(rl, gaps), rel_err_100(merge_rounding)
    start = start_preset(cf, err, fp16, b1)
    restart = restart_preset(cf, {"preset": start.preset, "passed": False}, fp16, b1)
    delta, share = opsm_counterfactual(gaps)
    return {"rel_err_100": err, "cf": cf, "logp_diff_41_60_b1": b1["logp_diff_41_60"], "fp16": fp16,
            "start": asdict(start), "restart": asdict(restart),
            "rule4": {"delta_q995_steps_1_20": delta, "masked_share_steps_50_260": share}}  # fmt: skip


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("decide", help="ADR-0013 rules 0-4 on the B1 run's records")
    d.add_argument("--b1", type=Path, required=True, help="the B1 run folder (rl.json, trainer/seq_gaps.f32)")
    d.add_argument("--merge-rounding", type=Path, required=True)
    d.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    out = decide(a.b1, a.merge_rounding)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
