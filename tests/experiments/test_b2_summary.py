"""B2.3 fp16-arm summary (I23) and the B2.4 decision's records: the B1 run's rl.json and seq_gaps (the rules
themselves are tests/experiments/test_phase_b_decisions.py)."""
import json

import numpy as np
import pytest

from tools.studies import b2_summary as b2


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


# B2.3 fp16-arm summary (I23)
def run_record(steps=70, iou=0.6, per_prompt=None, nan_at=None, per_step=20.0, seed=0):
    rng = np.random.default_rng(seed)
    log = []
    for s in range(1, steps + 1):
        h = log_entry(s, **{"rewards/iou_tc/mean": iou if not (41 <= s <= 60) else 0.5, "cf/ess_C": 0.5,
                            "sampling/importance_sampling_ratio/mean": 1.0, "loss": 0.1})
        if per_prompt is not None:
            h["rewards/iou_tc/by_prompt"] = per_prompt(rng, s)
        if s == nan_at:
            h["loss"] = float("nan")
        log.append(h)
    return {"log": log, "timing_s": {"per_step": per_step}}


def test_window_means_use_steps_41_to_60_only():
    m = b2.window_means(run_record())
    assert m["rewards/iou_tc/mean"] == pytest.approx(0.5) and m["cf/ess_C"] == 0.5
    assert m["sampling/sampling_logp_difference/mean"] == pytest.approx(0.02)


def test_nan_inf_count_covers_the_whole_run():
    assert b2.nan_inf_count(run_record()) == 0
    assert b2.nan_inf_count(run_record(nan_at=3)) == 1


def test_paired_bootstrap_lower_bound_matches_numpy_with_the_same_seed():
    rng = np.random.default_rng(1)
    a, b = rng.normal(0.55, 0.2, 640), rng.normal(0.6, 0.2, 640)
    d = a - b
    idx = np.random.default_rng(0).integers(0, 640, (2000, 640))
    want = float(np.quantile(d[idx].mean(axis=1), 0.025))
    assert b2.paired_bootstrap_lower(a, b, seed=0, n=2000) == pytest.approx(want)


def test_paired_bootstrap_refuses_unpaired_arms():
    with pytest.raises(ValueError, match="640"):
        b2.paired_bootstrap_lower(np.zeros(640), np.zeros(639))


def test_t_bound_is_the_one_sided_95_percent_bound_of_the_per_step_differences():
    from scipy import stats

    d = np.linspace(-0.02, 0.01, 20)
    want = d.mean() - stats.t.ppf(0.95, 19) * d.std(ddof=1) / np.sqrt(20)
    assert b2.t_lower(d) == pytest.approx(want)


def by_prompt(shift):
    return lambda rng, s: [float(x) for x in 0.6 + shift + 0.01 * np.sin(np.arange(32) + s)]


def test_summary_prints_means_bounds_and_sps(tmp_path, capsys):
    for name, rec in (("fp16", run_record(per_prompt=by_prompt(0.0), per_step=25.0)),
                      ("b1", run_record(per_prompt=by_prompt(0.0), per_step=21.0))):
        (tmp_path / name).mkdir()
        (tmp_path / name / "rl.json").write_text(json.dumps(rec))
    out = b2.summary(tmp_path / "fp16", tmp_path / "b1")
    assert out["bootstrap_lower"] == pytest.approx(0.0, abs=1e-9) and out["nan_inf"] == 0
    assert out["s_per_step"] == {"fp16": 25.0, "b1": 21.0} and out["prompts"] == 640
    assert "bootstrap" in capsys.readouterr().out


def test_summary_ignores_steps_outside_the_window(tmp_path):
    fp, ref = run_record(per_prompt=by_prompt(0.0)), run_record(per_prompt=by_prompt(0.0))
    for h in fp["log"]:
        if not 41 <= h["step"] <= 60:
            h["rewards/iou_tc/by_prompt"] = [0.0] * 32
            h["rewards/iou_tc/mean"] = 9.0
    for name, rec in (("fp16", fp), ("b1", ref)):
        (tmp_path / name).mkdir()
        (tmp_path / name / "rl.json").write_text(json.dumps(rec))
    out = b2.summary(tmp_path / "fp16", tmp_path / "b1")
    assert out["bootstrap_lower"] == pytest.approx(0.0, abs=1e-9) and out["means"]["fp16"]["rewards/iou_tc/mean"] == 0.5


def test_summary_without_per_prompt_rewards_says_so(tmp_path):
    for name in ("fp16", "b1"):
        (tmp_path / name).mkdir()
        (tmp_path / name / "rl.json").write_text(json.dumps(run_record()))
    with pytest.raises(ValueError, match="by_prompt"):
        b2.summary(tmp_path / "fp16", tmp_path / "b1")


def test_summary_missing_file_names_it(tmp_path):
    with pytest.raises(FileNotFoundError, match="rl.json"):
        b2.summary(tmp_path / "nope", tmp_path / "b1")


def test_fp16_record_is_the_gate_b2_input_of_rule_1(tmp_path):
    fp, ref = run_record(per_prompt=by_prompt(-0.019)), run_record(per_prompt=by_prompt(0.0))
    for h in fp["log"]:
        h["sampling/sampling_logp_difference/mean"] = 0.009
    for name, rec in (("fp16", fp), ("b1", ref)):
        (tmp_path / name).mkdir()
        (tmp_path / name / "rl.json").write_text(json.dumps(rec))
    (tmp_path / "fp16" / "f_check.json").write_text(json.dumps({"passed": True}))
    r = b2.fp16_record(tmp_path / "fp16", tmp_path / "b1")
    assert r["passed_gate"] is True and r["logp_diff_41_60"] == pytest.approx(0.009) and r["s_per_step"] == 20.0
    (tmp_path / "fp16" / "f_check.json").write_text(json.dumps({"passed": False}))
    assert b2.fp16_record(tmp_path / "fp16", tmp_path / "b1")["passed_gate"] is False  # NaN or dtype failure
    for h in fp["log"]:
        h["rewards/iou_tc/by_prompt"] = [x - 0.05 for x in h["rewards/iou_tc/by_prompt"]]
    (tmp_path / "fp16" / "rl.json").write_text(json.dumps(fp))
    (tmp_path / "fp16" / "f_check.json").write_text(json.dumps({"passed": True}))
    assert b2.fp16_record(tmp_path / "fp16", tmp_path / "b1")["passed_gate"] is False  # bound below -0.02


def make_b1(tmp_path, err=0.531, opsm05=0.01):
    run = tmp_path / "b1"
    (run / "trainer").mkdir(parents=True)
    write_gaps(run / "trainer", steps=300, rows=8)
    log = [log_entry(s, **{"cf/opsm_share@0.5": opsm05, "cf/tok_clamp2_share": 0.0, "cf/ess_E@0.1": 0.9})
           for s in range(1, 301)]
    (run / "rl.json").write_text(json.dumps({"log": log, "timing_s": {"per_step": 20.0}}))
    mr = tmp_path / "merge_rounding.json"
    mr.write_text(json.dumps({"rows": [{"step": 100, "bf16": {"rel_err": err}}]}))
    return run, mr


def test_decide_writes_both_presets_and_the_rule_4_counterfactual(tmp_path):
    run, mr = make_b1(tmp_path)
    out = b2.decide(run, mr)
    assert out["start"]["preset"] == "lora-seqmask" and out["start"]["start_delta"] is None
    assert out["restart"]["preset"] == "bf16-tis-opsm" and out["restart"]["start_delta"] == 0.5  # never the failed one
    assert 0.0 <= out["rule4"]["masked_share_steps_50_260"] <= 1.0 and out["rel_err_100"] == 0.531


def test_decide_cli_writes_the_file(tmp_path):
    run, mr = make_b1(tmp_path, err=0.1)
    assert b2.main(["decide", "--b1", str(run), "--merge-rounding", str(mr), "--out", str(tmp_path / "d.json")]) == 0
    assert json.loads((tmp_path / "d.json").read_text())["start"]["preset"] == "bf16-tis-opsm"
