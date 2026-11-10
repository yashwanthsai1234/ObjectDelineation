"""U-CF5 (plan C0, E-EQ6 / G-EQ6): train_rl --expect-step1 FILE --encode-check FILE stops the fast arm after step 1
when step 1 differs from another arm's, before step 2 starts."""

import json
from types import SimpleNamespace

import pytest

torch = pytest.importorskip("torch")
from objdel.policy import callbacks  # noqa: E402

R5 = {"rewards/iou_tc/mean": 0.5467, "rewards/legal/mean": 0.8984, "completions/mean_length": 117.89,
      "grad_norm": 0.003264, "sampling/sampling_logp_difference/mean": 0.009515}  # runs/r5/shared_rl.json step 1


def test_exact_keys_must_be_equal_and_the_two_float_keys_follow_the_encode_bits():
    """E-EQ6: grad_norm 0.003266 against 0.003264 (R5's measured pair) is 6.1e-4 apart: it passes at the 1e-3
    rule (encode bits differ) and fails the 4 significant digits rule (5e-4, encode bits equal)."""
    near = R5 | {"grad_norm": 0.003266}
    assert callbacks.step1_mismatches(R5, R5, True) == {}
    assert callbacks.step1_mismatches(R5, near, False) == {}
    assert callbacks.step1_mismatches(R5, near, True) == {"grad_norm": [0.003264, 0.003266]}
    far = R5 | {"rewards/iou_tc/mean": 0.5468, "completions/mean_length": 117.9}
    assert set(callbacks.step1_mismatches(R5, far, False)) == {"rewards/iou_tc/mean", "completions/mean_length"}
    off = R5 | {"sampling/sampling_logp_difference/mean": 0.01}  # 5 % apart: no rule passes it
    assert set(callbacks.step1_mismatches(R5, off, False)) == {"sampling/sampling_logp_difference/mean"}


def test_four_significant_digits_across_a_decade_edge_use_the_looser_rule():
    """0.1000 and 0.09993 round to 4 digits 7e-4 apart, in two decades: tolerated at 1e-3; two values in one decade
    that far apart are not."""
    a = R5 | {"grad_norm": 0.1}
    assert callbacks.step1_mismatches(a, a | {"grad_norm": 0.09993}, True) == {}
    b = R5 | {"grad_norm": 0.3000}
    assert set(callbacks.step1_mismatches(b, b | {"grad_norm": 0.29979}, True)) == {"grad_norm"}  # 7e-4, one decade


def write(path, doc):
    path.write_text(json.dumps(doc))
    return path


def rl_json(path, step1=R5, steps=(1, 2)):
    return write(path, {"log": [{"step": s, **step1, "loss": 0.0} for s in steps]})


def test_from_files_reads_step_1_of_the_other_arm_and_the_encode_bits(tmp_path):
    cb = callbacks.ExpectStep1.from_files(rl_json(tmp_path / "rl.json", steps=(2, 1)),
                                          write(tmp_path / "enc.json", {"bits_equal": True}), tmp_path / "out")
    assert cb.base == R5 and cb.bits_equal is True


@pytest.mark.parametrize(("enc", "what"), [
    (None, "missing"), ("not json", "unreadable"), ({}, "bits_equal"), ({"bits_equal": "yes"}, "bits_equal")])
def test_a_bad_encode_check_file_is_a_start_up_error(tmp_path, enc, what):
    f = tmp_path / "enc.json"
    if isinstance(enc, dict):
        write(f, enc)
    elif enc:
        f.write_text(enc)
    with pytest.raises((FileNotFoundError, ValueError)):
        callbacks.ExpectStep1.from_files(rl_json(tmp_path / "rl.json"), f, tmp_path / "out")


def test_a_bad_expect_step1_file_is_a_start_up_error(tmp_path):
    enc = write(tmp_path / "enc.json", {"bits_equal": True})
    with pytest.raises(FileNotFoundError):
        callbacks.ExpectStep1.from_files(tmp_path / "none.json", enc, tmp_path / "out")
    for doc in ({"log": []}, {"log": [{"step": 1, "grad_norm": 1.0}]}, {"nolog": 1}):
        with pytest.raises(ValueError, match="step 1"):
            callbacks.ExpectStep1.from_files(write(tmp_path / "x.json", doc), enc, tmp_path / "out")


def run(cb, step1, step=1):
    """One on_log call as the Trainer makes it: the log entry is the history's last."""
    control = SimpleNamespace(should_training_stop=False)
    state = SimpleNamespace(global_step=step, log_history=[{"step": step, **step1}])
    cb.on_log(None, state, control, logs=step1)
    return control


def test_equal_step_1_goes_on_and_a_mismatch_writes_the_file_and_stops(tmp_path):
    cb = callbacks.ExpectStep1(R5, True, tmp_path / "out")
    assert run(cb, R5).should_training_stop is False and cb.mismatch == {}
    assert not (tmp_path / "out/step1_mismatch.json").exists()
    cb = callbacks.ExpectStep1(R5, True, tmp_path / "out")
    assert run(cb, R5 | {"rewards/iou_tc/mean": 0.4}).should_training_stop is True
    doc = json.loads((tmp_path / "out/step1_mismatch.json").read_text())
    assert doc["mismatch"] == {"rewards/iou_tc/mean": [0.5467, 0.4]} and doc["bits_equal"] is True
    assert doc["base"] == R5 and doc["got"]["rewards/iou_tc/mean"] == 0.4


def test_only_the_step_1_train_log_is_compared_and_only_once(tmp_path):
    cb = callbacks.ExpectStep1(R5, True, tmp_path / "out")
    assert run(cb, {"eval_loss": 1.0}, step=1).should_training_stop is False and cb.mismatch is None  # an eval log
    assert run(cb, R5 | {"grad_norm": 9.0}, step=2).should_training_stop is False and cb.mismatch is None
    run(cb, R5, step=1)
    assert cb.mismatch == {}
    assert run(cb, R5 | {"grad_norm": 9.0}, step=1).should_training_stop is False  # compared once


def test_the_same_values_with_bits_equal_false_go_on(tmp_path):
    near = R5 | {"grad_norm": 0.003266}
    assert run(callbacks.ExpectStep1(R5, False, tmp_path / "a"), near).should_training_stop is False
    assert run(callbacks.ExpectStep1(R5, True, tmp_path / "b"), near).should_training_stop is True
