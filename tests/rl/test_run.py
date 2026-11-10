"""E2.3 (U54): run_rl checks the session budget first, before any rows, model or engine load; the library takes typed
arguments (build_trainer reads TrainerFlags, not a namespace)."""
import inspect
import json

import numpy as np
import pytest

from objdel.core.spec import RlStage, TrainerFlags
from objdel.rl import build, run


def test_run_rl_refuses_a_budget_without_a_session_start_before_anything_loads(tmp_path, monkeypatch):
    for var in ("OBJDEL_SESSION_START", "OBJDEL_RATE_USD_PER_H", "OBJDEL_FIXED_USD"):
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setattr(run, "rows_schedule", lambda *a: pytest.fail("rows read before the budget check"))
    monkeypatch.setattr(build, "build_trainer", lambda *a, **k: pytest.fail("trainer built"))
    stage = RlStage(model="qwen3-vl-2b", ckpt="c", train=str(tmp_path / "missing"), out=str(tmp_path / "o"),
                    steps=20, budget_usd=2.05, rate_usd_per_h=0.5)
    with pytest.raises(run.BudgetFlagError, match="OBJDEL_SESSION_START"):
        run.run_rl(stage)


def test_build_trainer_takes_typed_arguments():
    params = inspect.signature(build.build_trainer).parameters
    assert list(params)[:7] == ["model", "ckpt", "rows", "out", "roots", "steps", "flags"]
    assert params["flags"].annotation is TrainerFlags


def _no_constants(c):
    raise AssertionError(f"non-standard JSON constant {c}")


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf"), np.float64("nan"), np.float32("inf"),
                                 {"a": [1.0, {"b": float("nan")}]}, [np.array([0.5, np.nan])]])
def test_save_json_rejects_nan_and_inf_and_writes_nothing(tmp_path, bad):
    """U25: a NaN or inf anywhere raises ValueError; an existing file stays as it was and no temp file is left."""
    path = tmp_path / "out" / "rl.json"
    with pytest.raises(ValueError):
        run.save_json(path, {"x": bad}, ())
    assert not path.parent.exists()
    path.parent.mkdir()
    path.write_text("old\n")
    with pytest.raises(ValueError):
        run.save_json(path, {"x": bad}, ())
    assert path.read_text() == "old\n" and [p.name for p in path.parent.iterdir()] == ["rl.json"]


def test_save_json_writes_numpy_values_as_plain_json(tmp_path, capsys):
    """U25: numpy scalars and arrays (nested, 0-d, empty) become plain JSON numbers, booleans and lists."""
    doc = {"f": np.float32(0.5), "i": np.int64(3), "b": np.bool_(True), "a": np.arange(3),
           "nested": [{"z": np.array(2.0)}, np.zeros((0,))], "m": np.ones((2, 1), np.int32)}
    run.save_json(tmp_path / "x.json", doc, ("i", "b"))
    got = json.loads((tmp_path / "x.json").read_text(), parse_constant=_no_constants)
    assert got == {"f": 0.5, "i": 3, "b": True, "a": [0, 1, 2], "nested": [{"z": 2.0}, []], "m": [[1], [1]]}
    assert type(got["b"]) is bool and type(got["i"]) is int
    assert json.loads(capsys.readouterr().out) == {"i": 3, "b": True}


def test_finite_turns_the_nan_a_caller_expects_into_null(tmp_path):
    """U25: the callers whose values are NaN by design (g_probe's CIs, a trainer log) pass them through finite."""
    doc = run.finite({"ci": [0.1, np.float32("nan"), (float("inf"), 2)], "n": 3})
    assert doc == {"ci": [0.1, None, [None, 2]], "n": 3}
    run.save_json(tmp_path / "g.json", doc, ())
    assert json.loads((tmp_path / "g.json").read_text())["ci"][1] is None
