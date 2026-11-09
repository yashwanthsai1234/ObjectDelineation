"""Break-it round 1, the RL loop off the GPU: edge cases of train_rl, the callbacks, the reward functions and the
TRL contract. Each test proves one finding and is a strict xfail until the fix lands."""

import json
import math
from types import SimpleNamespace

import numpy as np
import pytest


@pytest.fixture
def stub(tmp_path, monkeypatch):
    """train_rl.main with no model, engine or preload: make_trainer returns a fake trainer whose train() runs
    `body(start)` (start = the resumed checkpoint's step, 0 on a fresh run)."""
    import experiments.phase2.g_probe as g_probe
    from experiments.phase2 import train_rl
    from objdel.train import dataset, serve

    monkeypatch.setattr(serve, "preload", lambda root, ds: None)
    monkeypatch.setattr(dataset, "grpo_rows", lambda *a, **k: ["k1", "k2"])
    monkeypatch.setattr(train_rl, "save_final", lambda model, proc, out: None)
    import gc

    monkeypatch.setattr(gc, "freeze", lambda: None)
    state = SimpleNamespace(global_step=0, log_history=[])
    body = {"fn": lambda start, steps: None}

    def train(resume_from_checkpoint=None):
        start = int(resume_from_checkpoint.rsplit("-", 1)[1]) if resume_from_checkpoint else 0
        body["fn"](start, state.max_steps)
        state.global_step = state.max_steps

    trainer = SimpleNamespace(train=train, model=None, accelerator=SimpleNamespace(unwrap_model=lambda m: m),
                              state=state, args=SimpleNamespace(beta=1e-3, mask_truncated_completions=True,
                                                                  off_policy_mask_threshold=None,
                                                                  max_completion_length=1536))

    def make(a, ckpt, train_rows, steps, **k):
        state.max_steps = steps
        return trainer, None, 5e-5

    monkeypatch.setattr(g_probe, "make_trainer", make)
    (tmp_path / "root").mkdir()
    (tmp_path / "root/manifest.json").write_text("{}")
    out = tmp_path / "out"
    base = ["--model", "qwen3-vl-2b", "--ckpt", str(tmp_path / "c"), "--out", str(out), "--train",
            str(tmp_path / "root")]
    return SimpleNamespace(main=train_rl.main, base=base, out=out, body=body)


@pytest.mark.xfail(strict=True, reason="BR1-rl-1: --steps 0 or negative is accepted; TRL max_steps <= 0 means "
                                       "train by num_train_epochs (3 full epochs of GPU spend), not 0 steps")
@pytest.mark.parametrize("steps", ["0", "-5"])
def test_steps_must_be_positive(steps):
    from experiments.phase2 import train_rl

    with pytest.raises(SystemExit):
        train_rl.parse(["--model", "qwen3-vl-2b", "--ckpt", "c", "--out", "o", "--train", "t", "--steps", steps])


@pytest.mark.xfail(strict=True, reason="BR1-rl-1: check_sampler accepts --steps 0 on a schedule root")
def test_check_sampler_refuses_zero_steps():
    from experiments.phase2 import train_rl

    cfg = SimpleNamespace(generation_batch_size=256, num_generations=8, steps_per_generation=16,
                          gradient_accumulation_steps=16, num_iterations=1)
    with pytest.raises(ValueError):
        train_rl.check_sampler(cfg, {"per_step": 32, "steps": 260, "block_steps": 20}, 0)


@pytest.mark.xfail(strict=True, reason="BR1-rl-2: --resume from checkpoint-N keeps the seq_gaps.f32 rows of the "
                                       "steps after N that the crashed run logged; the resumed run appends them again, "
                                       "so step i of the file is no longer optimizer step i")
def test_resume_mid_run_keeps_one_seq_gaps_row_block_per_step(stub):
    from objdel.train.fast import write_seq_gaps

    rows = 4
    trainer_dir = stub.out / "trainer"

    def run(start, steps):  # FastGRPO.log appends one block per logged step
        for s in range(start + 1, steps + 1):
            write_seq_gaps(str(trainer_dir), [float(s)] * rows * 3, rows)

    # first run: logs steps 1..7, saves checkpoint-5, then the box dies before checkpoint-10
    run(0, 7)
    (trainer_dir / "checkpoint-5").mkdir(parents=True)
    stub.body["fn"] = run
    assert stub.main([*stub.base, "--steps", "10", "--save-steps", "5", "--resume"]) == 0
    meta = json.loads((trainer_dir / "seq_gaps.json").read_text())
    data = np.fromfile(trainer_dir / "seq_gaps.f32", "<f4").reshape(-1, rows, 3)[:, 0, 0]
    assert meta["steps"] == 10
    assert data.tolist() == [float(s) for s in range(1, 11)]


@pytest.mark.xfail(strict=True, reason="BR1-rl-3: --merge-last treats any <out>/final directory as done, also one a "
                                       "crash left half-written (no weights): exit 0, nothing merged, no rl.json")
def test_merge_last_does_not_accept_a_half_written_final(tmp_path, capsys):
    from experiments.phase2 import train_rl

    out = tmp_path / "out"
    (out / "trainer" / "checkpoint-40").mkdir(parents=True)
    (out / "final_adapter").mkdir(parents=True)  # save_final got this far
    (out / "final").mkdir()  # save_pretrained created the folder, then the box died
    a = train_rl.parse(["--merge-last", "--model", "qwen3-vl-2b", "--ckpt", str(tmp_path / "c"), "--out", str(out)])
    try:
        rc = train_rl.merge_last(a)
    except Exception:  # a fixed version may try the merge on the fake checkpoint: any loud outcome is fine
        return
    assert rc != 0 or (out / "rl.json").exists()


def _control():
    return SimpleNamespace(should_training_stop=False, should_save=False)


@pytest.mark.xfail(strict=True, reason="BR1-rl-4: BudgetStop never trips (and never refuses) with a session start in "
                                       "the future, e.g. a millisecond timestamp, or a NaN rate or cap: spend() clamps "
                                       "the hours to 0 / compares NaN, so the budget guard is silently off")
@pytest.mark.parametrize("case", ["start_in_ms", "nan_rate", "nan_cap"])
def test_budget_stop_cannot_be_silently_off(case):
    from objdel.train.callbacks import BudgetStop

    now = 1_800_000_000.0
    kw = {"rate_usd_per_h": 1.0, "cap_usd": 2.0, "session_start": now - 10 * 3600, "clock": lambda: now}
    if case == "start_in_ms":
        kw["session_start"] = (now - 10 * 3600) * 1000
    elif case == "nan_rate":
        kw["rate_usd_per_h"] = math.nan
    else:
        kw["cap_usd"] = math.nan
    try:
        cb = BudgetStop(**kw)
    except ValueError:
        return
    c = _control()
    cb.on_step_end(None, SimpleNamespace(global_step=1), c)  # 10 h at $1/h is far past a $2 cap
    assert c.should_training_stop and c.should_save


@pytest.mark.xfail(strict=True, reason="BR1-rl-5: NvtxSteps never opens the nsys capture when the run resumes at "
                                       "step S - 1 or inside the window (the start test is only global_step == S - 1 "
                                       "at a step end), so a resumed profiling run captures nothing")
@pytest.mark.parametrize("resume_at", [11, 12])
def test_nvtx_capture_window_survives_resume(monkeypatch, resume_at):
    import torch

    from objdel.train.callbacks import NvtxSteps

    calls = []
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch.cuda.profiler, "start", lambda: calls.append("start"))
    monkeypatch.setattr(torch.cuda.profiler, "stop", lambda: calls.append("stop"))
    cb, state = NvtxSteps((12, 14)), SimpleNamespace(global_step=resume_at)
    cb.on_train_begin(None, state, _control())
    for s in range(resume_at + 1, 16):
        state.global_step = s
        cb.on_step_end(None, state, _control())
    cb.on_train_end(None, state, _control())
    assert calls == ["start", "stop"]


def _fake_terms(n, iou):
    T = {k: np.ones(n) for k in ("format_ok", "legal", "f5", "g_area", "topo_any", "topo_freq", "iou_min")}
    T |= {"iou_tc": np.asarray(iou, float), "sig_id": np.zeros((n, 1)), "error": [""] * n, "defects": [()] * n}
    return T


@pytest.mark.xfail(strict=True, reason="BR1-rl-6: one NaN term of a legal answer (e.g. a 0/0 IoU from the checker) "
                                       "turns its whole group's advantages into NaN and TRLReward hands them to TRL: "
                                       "the loss and the LoRA weights go NaN with no error")
def test_nan_term_does_not_reach_trl_as_nan_advantages(monkeypatch):
    from objdel.train import rewards

    g = 8
    iou = np.linspace(0.2, 0.9, g)
    iou[3] = np.nan
    monkeypatch.setattr(rewards, "_score", lambda jobs: _fake_terms(len(jobs), iou))
    r = rewards.TRLReward(lambda k: {}, stop_ids=[2])
    funcs, _ = r.funcs()
    agg = funcs[-1]
    try:
        adv = agg(completions=["x"] * g, completion_ids=[[1, 2]] * g, key=["a"] * g, W=[64] * g, H=[64] * g,
                  names=[["o"]] * g)
    except ValueError:
        return
    assert np.isfinite(adv).all()


@pytest.mark.xfail(strict=True, reason="BR1-rl-7: grpo_config accepts a per_device_eval_batch_size that is not a "
                                       "multiple of G (the rewards docstring says configs guarantees it); the run then "
                                       "dies at its first held-out eval, after --heldout-every paid steps")
def test_grpo_config_refuses_eval_batch_not_multiple_of_group(tmp_path):
    from objdel.train import configs

    with pytest.raises(ValueError):
        configs.grpo_config("qwen3-vl-2b", str(tmp_path), per_device_eval_batch_size=12, use_vllm=False)


@pytest.mark.xfail(strict=True, reason="BR1-rl-8: trl_adapter.check passes when TRL stops passing log_metric / "
                                       "log_extra to reward functions; TRLReward then silently drops error/defect/"
                                       "allmiss metrics and the completions-table columns sample.export needs")
def test_contract_catches_dropped_reward_kwargs(monkeypatch):
    import trl.trainer.grpo_trainer as grpo

    from objdel.train import trl_adapter

    def _calculate_rewards(self, inputs, prompts, completions, completion_ids_list):  # same signature, no kwargs
        return None

    monkeypatch.setattr(grpo.GRPOTrainer, "_calculate_rewards", _calculate_rewards)
    with pytest.raises(trl_adapter.TRLContractError):
        trl_adapter.check()


def test_stub_runs_a_fresh_short_run(stub):  # the harness itself works (not an xfail)
    assert stub.main([*stub.base, "--steps", "2"]) == 0
    assert json.loads((stub.out / "rl.json").read_text())["steps"] == 2
