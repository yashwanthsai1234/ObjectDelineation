"""B2.4 decision rules 0-4 (U51): `start_preset` and `restart_preset`, pure functions over the B1 run's records
(the records and the fp16-arm summary stay in tests/test_b2_summary.py)."""
import numpy as np
import pytest

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
