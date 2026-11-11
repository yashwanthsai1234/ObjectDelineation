"""Break-it round 2, the RL loop off the GPU: edge cases of the round-1 fixes and the parts round 1 did not reach
(TraceStep, merge_last, stale checkpoints, the --grpo overrides, the env defaults). Each test proves one finding
and ran xfail(strict) until it was fixed."""

import sys
from types import SimpleNamespace

import pytest

from objdel.rl import run as rl_run


def _trace_steps(tmp_path, wait: int, start: int, steps: int) -> list[int]:
    """Run TraceStep over global steps start+1..steps; return the global steps at which it wrote a trace."""
    from objdel.policy.callbacks import TraceStep

    cb, state, wrote = TraceStep(str(tmp_path / "trace"), wait=wait), SimpleNamespace(global_step=start), []
    cb.write = lambda prof: wrote.append(state.global_step)
    cb.on_train_begin(None, state, None)
    for s in range(start + 1, steps + 1):
        state.global_step = s
        cb.on_step_end(None, state, None)
    cb.on_train_end(None, state, None)
    return wrote


def test_trace_step_profiles_one_step_only(tmp_path):
    assert len(_trace_steps(tmp_path, wait=1, start=0, steps=12)) == 1


def test_trace_step_matches_nvtx_window_after_resume(tmp_path):
    s = 5  # build_trainer: TraceStep(wait=S - 2); NvtxSteps opens the capture at the end of step S - 1 = 4
    assert _trace_steps(tmp_path, wait=s - 2, start=3, steps=12)[:1] == [s]


def test_merge_last_keeps_a_complete_final_without_rl_json(tmp_path, monkeypatch):
    import transformers  # noqa: F401
    import transformers.trainer_utils  # noqa: F401  (imported before the fake peft below)
    from transformers import AutoModelForImageTextToText, AutoProcessor  # noqa: F401

    import objdel.rl.run  # noqa: F401
    pass  # noqa: F401  (imports the real peft)
    from objdel.cli.train import rl as train_rl

    out = tmp_path / "out"
    (out / "trainer" / "checkpoint-250").mkdir(parents=True)
    (out / "final_adapter").mkdir()
    (out / "final").mkdir()
    (out / "final" / "config.json").write_text("{}")
    (out / "final" / "model.safetensors").write_bytes(b"weights of step 260")

    class Peft:
        @staticmethod
        def from_pretrained(*a, **k):
            raise RuntimeError("fake: no merge in the test")

    monkeypatch.setitem(sys.modules, "peft", SimpleNamespace(PeftModel=Peft))
    a = train_rl.parse(["--merge-last", "--model", "qwen3-vl-2b", "--ckpt", str(tmp_path / "c"), "--out", str(out)])
    try:
        rl_run.merge_last(a)
    except (RuntimeError, OSError):  # the fake merge, or no base model at --ckpt
        pass
    assert (out / "final" / "model.safetensors").read_bytes() == b"weights of step 260"


@pytest.fixture
def stub(tmp_path, monkeypatch):
    """train_rl.main with no model, engine or preload (as round 1's stub); body(start, steps) runs in train()."""
    import gc

    from objdel.cli.train import rl as train_rl
    from objdel.rl import build
    from objdel.rl import rows as dataset
    from objdel.verify import store as serve

    monkeypatch.setattr(serve, "preload", lambda root, ds: None)
    monkeypatch.setattr(dataset, "grpo_rows", lambda *a, **k: ["k1", "k2"])
    monkeypatch.setattr(rl_run, "save_final", lambda model, proc, out: None)
    monkeypatch.setattr(gc, "freeze", lambda: None)
    state = SimpleNamespace(global_step=0, log_history=[])
    body = {"fn": lambda start, steps: None}

    def train(resume_from_checkpoint=None):
        start = int(resume_from_checkpoint.rsplit("-", 1)[1]) if resume_from_checkpoint else 0
        body["fn"](start, state.max_steps)
        state.global_step = state.max_steps

    trainer = SimpleNamespace(train=train, model=None, accelerator=SimpleNamespace(unwrap_model=lambda m: m),
                              state=state, args=SimpleNamespace(beta=1e-3, mask_truncated_completions=True))

    def make(model, ckpt, rows, out, roots, steps, flags, **k):
        state.max_steps = steps
        return trainer, None, 5e-5

    monkeypatch.setattr(build, "build_trainer", make)
    (tmp_path / "root").mkdir()
    (tmp_path / "root/manifest.json").write_text('{"sets": {"lidc/train": {}}}')  # the train set check_root looks for
    out = tmp_path / "out"
    base = ["--model", "qwen3-vl-2b", "--ckpt", str(tmp_path / "c"), "--out", str(out), "--train",
            str(tmp_path / "root")]
    return SimpleNamespace(main=train_rl.main, base=base, out=out, body=body)


class _Crash(Exception):
    pass


def test_fresh_run_then_resume_ignores_stale_checkpoints(stub):
    trainer_dir = stub.out / "trainer"
    (trainer_dir / "checkpoint-250").mkdir(parents=True)  # an earlier run in the same --out

    def crash(start, steps):
        (trainer_dir / "checkpoint-5").mkdir()
        raise _Crash

    stub.body["fn"] = crash
    try:
        stub.main([*stub.base, "--steps", "10", "--save-steps", "5"])
    except _Crash:
        pass
    except Exception:  # a fixed version may refuse the stale --out up front: also fine
        return
    seen = []
    stub.body["fn"] = lambda start, steps: seen.append(start)
    stub.main([*stub.base, "--steps", "10", "--save-steps", "5", "--resume"])
    assert seen == [5]


def test_seq_gaps_stay_one_block_per_step_with_logging_steps_override():
    """BR2-rl-5 fix: logging_steps != 1 would write several optimizer steps per seq_gaps block (file step i != step
    i after cut_seq_gaps), so train_rl refuses the override before any load."""
    from objdel.cli.train import rl as train_rl

    with pytest.raises(SystemExit):
        train_rl.parse(["--model", "qwen3-vl-2b", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", "10",
                        "--grpo", '{"logging_steps": 5}'])


@pytest.mark.parametrize("var", ["OBJDEL_FIXED_USD", "OBJDEL_SESSION_START", "OBJDEL_RATE_USD_PER_H"])
def test_empty_budget_env_does_not_break_a_run_without_budget(monkeypatch, var):
    """BR2-rl-6, fixed at E6.4: the stage reads no environment; session_budget reads it only under --budget-usd."""
    from objdel.cli.train import rl as train_rl

    monkeypatch.setenv(var, "")
    a = train_rl.parse(["--model", "qwen3-vl-2b", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", "10"])
    assert a.budget_usd is None
