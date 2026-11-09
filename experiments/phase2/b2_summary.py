"""B2.4 decision (ADR-0013 rules 0-4, $0) and the B2.3 fp16-arm summary.

    python -I -m experiments.phase2.b2_summary decide --b1 RD/b1-rows --merge-rounding RD/merge_rounding.json
        [--fp16 RD/b2-fp16] --out RD/decision.json
    python -I -m experiments.phase2.b2_summary summary --fp16 RD/b2-fp16 --b1 RD/b1-rows

The fp16 arm's bootstrap needs the per-prompt `rewards/iou_tc/by_prompt` (32 values per step) in each rl.json log
entry; TRL's own log holds only the step mean (see impl-log/B2.md).

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
FP16_BOUND = -0.02  # Gate B2: the 95 % lower bound of (fp16 - B1) rewards/iou_tc must exceed it
BY_PROMPT = "rewards/iou_tc/by_prompt"
WINDOW_KEYS = ("sampling/importance_sampling_ratio/mean", "sampling/sampling_logp_difference/mean",
               "rewards/iou_tc/mean")  # plus every cf/ and opsm_ key B2.1 logs


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


def window_means(rl: dict) -> dict:
    """Mean over steps 41-60 of the B2.1 metrics (cf/*, opsm_*) and of WINDOW_KEYS present in the log."""
    rows = _window(rl["log"], *FP16_STEPS)
    keys = [k for k in rows[0] if k in WINDOW_KEYS or k.startswith(("cf/", "sampling/opsm_"))]
    return {k: _mean(rows, k) for k in keys}


def nan_inf_count(rl: dict) -> int:
    """Non-finite numbers in the whole log (logging_nan_inf_filter is off, so a NaN shows up there)."""
    return sum(1 for h in rl["log"] for v in h.values() if isinstance(v, float) and not np.isfinite(v))


def paired_bootstrap_lower(a: np.ndarray, b: np.ndarray, seed: int = 0, n: int = 10_000) -> float:
    """95 % lower bound (2.5th percentile) of mean(a - b), resampling the paired prompts."""
    if a.shape != b.shape or a.ndim != 1:
        raise ValueError(f"the arms need the same prompts in the same order (640 over steps 41-60): {a.shape} vs {b.shape}")
    d = a - b
    idx = np.random.default_rng(seed).integers(0, d.size, (n, d.size))
    return float(np.quantile(d[idx].mean(axis=1), 0.025))


def t_lower(d: np.ndarray) -> float:
    """One-sided 95 % t lower bound of the mean of the per-step differences (a second line: autocorrelated)."""
    from scipy import stats

    return float(d.mean() - stats.t.ppf(0.95, d.size - 1) * d.std(ddof=1) / np.sqrt(d.size))


def _read(run: Path) -> dict:
    path = Path(run) / "rl.json"
    if not path.is_file():
        raise FileNotFoundError(f"{path}: no rl.json (is the run finished?)")
    return json.loads(path.read_text())


def _by_prompt(rl: dict, name: str) -> np.ndarray:
    rows = _window(rl["log"], *FP16_STEPS)
    if any(BY_PROMPT not in h for h in rows):
        raise ValueError(f"{name}: steps {FP16_STEPS[0]}-{FP16_STEPS[1]} lack {BY_PROMPT}")
    return np.array([h[BY_PROMPT] for h in rows], dtype=float).ravel()


def summary(fp16_run: Path, b1_run: Path) -> dict:
    """B2.3: the fp16 arm against the B1 run over steps 41-60; prints it and returns it."""
    fp, ref = _read(fp16_run), _read(b1_run)
    a, b = _by_prompt(fp, "fp16 run"), _by_prompt(ref, "B1 run")
    steps = FP16_STEPS[1] - FP16_STEPS[0] + 1
    per_step = (a - b).reshape(steps, -1).mean(axis=1)
    out = {"means": {"fp16": window_means(fp), "b1": window_means(ref)},
           "nan_inf": nan_inf_count(fp), "s_per_step": {"fp16": fp["timing_s"]["per_step"],
                                                        "b1": ref["timing_s"]["per_step"]},
           "prompts": int(a.size), "bootstrap_lower": paired_bootstrap_lower(a, b),
           "t_lower": t_lower(per_step)}  # fmt: skip
    print(json.dumps(out, indent=1))
    print(f"paired prompt bootstrap 95 % lower bound of (fp16 - B1) rewards/iou_tc: {out['bootstrap_lower']:.4f}"
          f" (Gate B2 needs > {FP16_BOUND}); one-sided t bound over steps: {out['t_lower']:.4f} (second line only)")
    return out


def fp16_record(fp16_run: Path, b1_run: Path) -> dict:
    """The `fp16` argument of the decision functions: Gate B2's fp16 rules (the run's f_check.json `passed`, which
    holds NaN/inf = 0 and the engine dtype, and the bootstrap bound) and the logp difference over steps 41-60."""
    out = summary(fp16_run, b1_run)
    f_check = json.loads((Path(fp16_run) / "f_check.json").read_text())
    return {"passed_gate": bool(f_check["passed"] and out["bootstrap_lower"] > FP16_BOUND),
            "logp_diff_41_60": out["means"]["fp16"]["sampling/sampling_logp_difference/mean"],
            "s_per_step": out["s_per_step"]["fp16"]}


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
    d.add_argument("--fp16", type=Path, help="the fp16 arm's run folder (rl.json, f_check.json); absent: rule 1 skipped")
    d.add_argument("--out", type=Path, required=True)
    s = sub.add_parser("summary", help="B2.3: the fp16 arm against the B1 run over steps 41-60")
    s.add_argument("--fp16", type=Path, required=True)
    s.add_argument("--b1", type=Path, required=True)
    a = ap.parse_args(argv)
    if a.cmd == "summary":
        summary(a.fp16, a.b1)
        return 0
    out = decide(a.b1, a.merge_rounding, fp16_record(a.fp16, a.b1) if a.fp16 else None)
    a.out.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
