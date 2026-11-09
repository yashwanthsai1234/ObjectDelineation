"""B2.4 decision rules 0-4 (U51) and the B2.3 fp16-arm summary (I23): pure functions over the B1 run's records."""
import json

import numpy as np
import pytest

from experiments.phase2 import b2_summary as b2
from experiments.phase2.b2_summary import restart_preset, start_preset

CF = lambda opsm05=0.01, clamp=0.005, ess_e=0.9: {  # noqa: E731  medians over B1 steps 50-260
    "cf/opsm_share@0.5": opsm05, "cf/tok_clamp2_share": clamp, "cf/ess_E@0.1": ess_e}
FP = lambda ok, diff, sps: dict(passed_gate=ok, logp_diff_41_60=diff, s_per_step=sps)  # noqa: E731


def F(preset="lora-seqmask", passed=True, **change):
    """A full f_check.json record (B2.2 F step 5); `passed` is set by train_rl, never recomputed here."""
    good = dict(preset=preset, nan=0, lora_equal=True, load_count=1, base_hash_ok=True, lora_over_base=1.05,
                comparator=dict(f_merged=0.0176, f_lora=0.0096, f0=0.0096, ratio=0.0, informative=True, ok=True))
    return good | change | dict(passed=passed)


CMP_051 = dict(f_merged=0.0176, f_lora=0.0137, f0=0.0096, ratio=0.51, informative=True, ok=False)
SEQ_KL = np.random.default_rng(0).gamma(2.0, 0.05, 211 * 256).astype(np.float32)
B1 = dict(logp_diff_41_60=0.020, s_per_step=20.0, seq_kl_50_260=SEQ_KL)


@pytest.mark.parametrize("rel_err,fp,want", [
    (0.531, None, "lora-seqmask"),  # no F result exists yet: rule 0 fires on the B2.0 error alone
    (0.2, None, "lora-seqmask"),  # boundary: >= 0.2
    (0.531, FP(True, 0.010, 21.0), "lora-seqmask"),  # rule 0 comes before rule 1
    (0.1999, None, "bf16-tis-opsm"),  # rule 0 skipped, no fp16 arm
    (0.1999, FP(True, 0.010, 40.0), "fp16-seqmask"),  # rule 1, exactly half; 2x the B1 s/step is a report line
    (0.1999, FP(True, 0.011, 20.0), "bf16-tis-opsm"),  # fp16 above half
    (0.1999, FP(False, 0.001, 20.0), "bf16-tis-opsm"),  # fp16 failed its gate (NaN, bootstrap bound)
])
def test_start_preset_reads_only_the_b20_error(rel_err, fp, want):
    assert start_preset(CF(), rel_err_100=rel_err, fp16=fp, b1=B1).preset == want


@pytest.mark.parametrize("f,fp,want", [
    (F(), None, "lora-seqmask"),  # passed: the start preset is kept
    (F(passed=False, comparator=CMP_051), None, "bf16-tis-opsm"),  # a comparator-only failure never keeps F
    (F(passed=False, comparator=CMP_051), FP(True, 0.010, 21.0), "fp16-seqmask"),
    (F(passed=False, lora_over_base=1.151), FP(True, 0.010, 40.0), "fp16-seqmask"),
    (F(passed=False, base_hash_ok=False), FP(True, 0.011, 20.0), "bf16-tis-opsm"),
    (F(passed=False), None, "bf16-tis-opsm"),  # every field flag true, `passed` false: the flag decides
    (F("fp16-seqmask", passed=False), FP(True, 0.010, 21.0), "bf16-tis-opsm"),  # never fp16 again
    (F("fp16-seqmask"), FP(True, 0.010, 21.0), "fp16-seqmask"),
])
def test_restart_preset_after_the_f_check(f, fp, want):
    assert restart_preset(CF(), f_check=f, fp16=fp, b1=B1).preset == want
    assert want == f["preset"] or not f["passed"]  # a restart never returns the failed preset
    assert want != f["preset"] or f["passed"]


def test_restart_preset_veto():
    assert restart_preset(CF(ess_e=0.499), f_check=F(passed=False), fp16=None, b1=B1).preset == "bf16-seqmask"


def test_vetoes_on_e_never_maximise_ess():
    d = start_preset(CF(opsm05=0.021), rel_err_100=0.1, fp16=None, b1=B1)
    assert d.preset == "bf16-tis-opsm" and isinstance(d.start_delta, float)
    assert d.start_delta == pytest.approx(float(np.quantile(SEQ_KL, 0.999)))
    assert start_preset(CF(opsm05=0.020), rel_err_100=0.1, fp16=None, b1=B1).start_delta == 0.5
    assert restart_preset(CF(opsm05=0.021), f_check=F(passed=False), fp16=None, b1=B1).start_delta == d.start_delta
    assert start_preset(CF(), rel_err_100=0.531, fp16=None, b1=B1).start_delta is None  # lora-seqmask: no OPSM
    assert start_preset(CF(ess_e=0.499), rel_err_100=0.1, fp16=None, b1=B1).start_delta is None  # seqmask: no OPSM
    clamp = start_preset(CF(clamp=0.011), rel_err_100=0.1, fp16=None, b1=B1)
    assert clamp.preset == "bf16-tis-opsm" and any("clamp" in n for n in clamp.notes)  # recorded, preset kept
    assert start_preset(CF(clamp=0.010), rel_err_100=0.1, fp16=None, b1=B1).notes == ()
    assert start_preset(CF(ess_e=0.499), rel_err_100=0.1, fp16=None, b1=B1).preset == "bf16-seqmask"
    assert start_preset(CF(ess_e=0.5), rel_err_100=0.1, fp16=None, b1=B1).preset == "bf16-tis-opsm"


# the records the decision reads
def write_gaps(d, steps=300, rows=4, seed=0):
    g = np.random.default_rng(seed).normal(size=(steps, rows, 3)).astype("<f4")
    (d / "seq_gaps.f32").write_bytes(g.tobytes())
    (d / "seq_gaps.json").write_text(json.dumps({"rows_per_step": rows, "steps": steps}))
    return g


def test_load_seq_gaps_reads_the_f32_file_by_the_json_row_count(tmp_path):
    g = write_gaps(tmp_path)
    assert np.array_equal(b2.load_seq_gaps(tmp_path), g) and b2.load_seq_gaps(tmp_path).shape == (300, 4, 3)


def test_load_seq_gaps_refuses_a_torn_file(tmp_path):
    write_gaps(tmp_path)
    with open(tmp_path / "seq_gaps.f32", "ab") as f:
        f.write(b"\0\0\0\0")
    with pytest.raises(ValueError, match="whole steps"):
        b2.load_seq_gaps(tmp_path)


def log_entry(step, **kv):
    return {"step": step, "sampling/sampling_logp_difference/mean": 0.02, "rewards/iou_tc/mean": 0.6, **kv}


def test_cf_medians_use_steps_50_to_260_only():
    log = [log_entry(s, **{"cf/opsm_share@0.5": 0.9 if s in (49, 261) else 0.01 + s * 1e-6, "cf/tok_clamp2_share": 0.0,
                           "cf/ess_E@0.1": 0.8}) for s in range(1, 301)]
    cf = b2.cf_medians(log)
    assert cf["cf/opsm_share@0.5"] == pytest.approx(0.01 + 155 * 1e-6) and cf["cf/ess_E@0.1"] == 0.8


def test_cf_medians_name_a_missing_stat_and_a_missing_step():
    with pytest.raises(ValueError, match=r"cf/ess_E@0.1"):
        b2.cf_medians([log_entry(s, **{"cf/opsm_share@0.5": 0.0, "cf/tok_clamp2_share": 0.0}) for s in range(1, 300)])
    with pytest.raises(ValueError, match="step 260"):
        b2.cf_medians([log_entry(s, **{"cf/opsm_share@0.5": 0.0, "cf/tok_clamp2_share": 0.0, "cf/ess_E@0.1": 1.0})
                       for s in range(1, 260)])


def test_b1_inputs(tmp_path):
    g = write_gaps(tmp_path)
    log = [log_entry(s, **{"sampling/sampling_logp_difference/mean": 0.01 if 41 <= s <= 60 else 9.0})
           for s in range(1, 301)]
    b1 = b2.b1_inputs({"log": log, "timing_s": {"per_step": 21.5}}, g)
    assert b1["logp_diff_41_60"] == pytest.approx(0.01) and b1["s_per_step"] == 21.5
    assert np.array_equal(b1["seq_kl_50_260"], g[49:260, :, 0].ravel())


def test_rule4_counterfactual_applies_the_first_20_steps_q995_to_steps_50_260():
    g = np.zeros((300, 10, 3), dtype=np.float32)
    g[:20, :, 0] = np.arange(200).reshape(20, 10) / 199.0  # q99.5 of 0..1 on 200 values
    g[49:260, :, 0] = 2.0  # every sequence above delta in the window
    g[49:260, :5, 1] = -1.0  # half of them with advantage < 0
    delta, share = b2.opsm_counterfactual(g)
    assert delta == pytest.approx(float(np.quantile(g[:20, :, 0], 0.995))) and share == pytest.approx(0.5)


def test_rule0_error_is_the_step_100_bf16_rel_err(tmp_path):
    p = tmp_path / "m.json"
    rows = [{"step": 50, "bf16": {"rel_err": 0.6}}, {"step": 100, "bf16": {"rel_err": 0.53}}]
    p.write_text(json.dumps({"rows": rows}))
    assert b2.rel_err_100(p) == 0.53
    p.write_text(json.dumps({"rows": [{"step": 50, "bf16": {"rel_err": 0.6}}]}))
    with pytest.raises(ValueError, match="step 100"):
        b2.rel_err_100(p)
