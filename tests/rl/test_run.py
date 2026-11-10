"""E2.3 (U54): run_rl checks the session budget first, before any rows, model or engine load; the library takes typed
arguments (build_trainer reads TrainerFlags, not a namespace)."""
import inspect

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
