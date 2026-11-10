"""B2.2 F step 5, B2.4 rule 0 (U76, Mac twin of G05): the F check record, its pass rules and the `train_rl --f-check`
flag, with fake probes (no engine). The comparator arms, the packed-layout equality and the base restore need vLLM
and are not built yet."""
import json
from types import SimpleNamespace

import pytest

from experiments.phase2 import train_rl
from objdel.train import f_check as fc

BASE = ["--model", "qwen3-vl-2b", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "20"]


def lora_fixed(**change):
    return {"comparator": fc.comparator(0.0176, 0.0096, 0.0096), "base_hash_ok": True, "lora_over_base": 1.05,
            "lora_equal": True} | change


def run(tmp_path, preset="lora-seqmask", steps=15, nan_at=None, loads_per_step=1, **kw):
    """Drive the callback like TRL: one on_log and one on_step_end per optimizer step."""
    loads = fc.count_calls(lambda: None)
    probes = ({"engine_dtype_ok": lambda: kw.pop("dtype_ok", True)} if preset == "fp16-seqmask" else {})
    cb = fc.FCheck(tmp_path, preset, probes, loads=loads if preset == "lora-seqmask" else None,
                   fixed=kw.pop("fixed", lora_fixed() if preset == "lora-seqmask" else {}))
    for step in range(1, steps + 1):
        for _ in range(loads_per_step):
            loads()
        cb.on_log(None, SimpleNamespace(global_step=step), None, {"loss": float("nan") if step == nan_at else 0.1})
        cb.on_step_end(None, SimpleNamespace(global_step=step), None)
    return json.loads((tmp_path / "f_check.json").read_text())


def test_all_checks_pass_writes_the_record_after_step_15_and_the_run_goes_on(tmp_path):
    rec = run(tmp_path, steps=20)
    assert rec["passed"] is True and rec["preset"] == "lora-seqmask" and rec["nan"] == 0 and rec["load_count"] == 1
    fc.check_record(rec)


@pytest.mark.parametrize("fixed,flag", [
    ({"lora_equal": False}, "lora_equal"), ({"base_hash_ok": False}, "base_hash_ok"),
    ({"lora_over_base": 1.151}, "lora_over_base"),
    ({"comparator": fc.comparator(0.0176, 0.0140, 0.0096)}, "comparator")])
def test_a_failed_check_writes_passed_false_then_raises(tmp_path, fixed, flag):
    with pytest.raises(RuntimeError, match="F check failed"):
        run(tmp_path, fixed=lora_fixed(**fixed))
    rec = json.loads((tmp_path / "f_check.json").read_text())
    assert rec["passed"] is False and not (rec[flag]["ok"] if flag == "comparator" else rec[flag] is True)


def test_nan_at_step_9_raises_there_and_two_loads_in_one_step_fail(tmp_path):
    with pytest.raises(RuntimeError):
        run(tmp_path, nan_at=9)
    assert json.loads((tmp_path / "f_check.json").read_text())["nan"] == 1
    with pytest.raises(RuntimeError):
        run(tmp_path / "two", loads_per_step=2)
    assert json.loads((tmp_path / "two/f_check.json").read_text())["load_count"] == 2


def test_fp16_record_is_nan_and_the_engine_dtype(tmp_path):
    rec = run(tmp_path / "ok", "fp16-seqmask")
    assert rec == {"preset": "fp16-seqmask", "nan": 0, "engine_dtype_ok": True, "passed": True}
    with pytest.raises(RuntimeError):
        run(tmp_path / "bad", "fp16-seqmask", dtype_ok=False)
    assert json.loads((tmp_path / "bad/f_check.json").read_text())["passed"] is False


@pytest.mark.parametrize("f0,merged,lora,ok,informative", [
    (0.0096, 0.0176, 0.0096, True, True),
    (0.0096, 0.0176, 0.0140, False, True),  # ratio 0.55
    (0.0096, 0.0176, 0.0130, True, True),  # ratio 0.425
    (0.0096, 0.0115, 0.0115, True, False)])  # merged - f0 = 0.0019: cannot tell F from no F
def test_comparator_subtracts_the_floor(f0, merged, lora, ok, informative):
    c = fc.comparator(merged, lora, f0)
    assert (c["ok"], c["informative"]) == (ok, informative) and (c["ratio"] is None) == (not informative)


def test_timing_ratio_is_the_median_and_1_15_passes():
    assert fc.lora_over_base([1.15] * 6, [1.0] * 6) == 1.15
    assert fc.passes({"preset": "lora-seqmask", "nan": 0, "load_count": 1, "base_hash_ok": True, "lora_equal": True,
                      "comparator": {"ok": True}, "lora_over_base": 1.15})
    assert fc.lora_over_base([1.0, 1.0, 9.0], [1.0, 1.0, 1.0]) == 1.0


def test_check_record_never_trusts_passed_alone():
    good = {"preset": "lora-seqmask", "nan": 0, "load_count": 1, "base_hash_ok": True, "lora_equal": True,
            "comparator": {"ok": True}, "lora_over_base": 1.0, "passed": True}
    fc.check_record(good)
    with pytest.raises(AssertionError):
        fc.check_record(good | {"lora_equal": False})  # one flag false, `passed` true
    with pytest.raises(AssertionError):
        fc.check_record(good | {"lora_equal": False, "passed": False})  # consistent, but a failure


def test_f_check_flag_is_on_by_default_under_fp16_and_an_argparse_error_under_bf16(capsys):
    assert train_rl.parse([*BASE, "--numerics", "fp16-seqmask"]).f_check is True
    assert train_rl.parse([*BASE, "--numerics", "fp16-seqmask", "--no-f-check"]).f_check is False
    assert train_rl.parse(BASE).f_check is False
    with pytest.raises(SystemExit):
        train_rl.parse([*BASE, "--f-check"])
    assert "--f-check" in capsys.readouterr().err
