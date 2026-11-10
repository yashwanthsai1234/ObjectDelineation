"""Controller tests moved from tests/test_breakit_reward.py (CLEAN6, from commit f531f86). Archived: never collected, linted or
run. To run them, restore the Controller from the tag pre-clean-2026-11-08 (README.md)."""
from types import SimpleNamespace

import pytest
from callbacks import Controller  # archive/controller/callbacks.py


def fake_trainer(tmp_path):
    torch = pytest.importorskip("torch")
    model = torch.nn.Linear(2, 2)
    opt = torch.optim.AdamW(model.parameters())
    model(torch.ones(1, 2)).sum().backward()
    opt.step()
    return SimpleNamespace(args=SimpleNamespace(output_dir=str(tmp_path), process_index=0),
                           state=SimpleNamespace(global_step=0), model=model, optimizer=opt,
                           reward_func_names=["g_area", "advantage_v4"], reward_weights=torch.tensor([0.0, 1.0]),
                           beta=1e-3, ref_model=None)


def snapshot(model):
    return {n: p.detach().clone() for n, p in model.named_parameters()}


def test_b2_1_alarm_rolls_back_to_the_checkpoint_before_the_hacked_evaluation(tmp_path):
    torch = pytest.importorskip("torch")
    t = fake_trainer(tmp_path)
    metric = "eval_rewards/iou_tc/mean"
    ctl = Controller(t, metric=metric, alarm=lambda m: "area_median" if m.get("hedging") else None)
    ctl.on_evaluate(t.args, t.state, None, metrics={metric: 0.80})  # honest checkpoint, best
    honest = snapshot(t.model)
    with torch.no_grad():
        for p in t.model.parameters():
            p.add_(1.0)  # the policy learns to hedge (union / hull / core: higher IoU_tc, M8 demo)
    t.state.global_step = 50
    ctl.on_evaluate(t.args, t.state, None, metrics={metric: 0.90, "hedging": True})
    assert all(torch.equal(p, honest[n]) for n, p in t.model.named_parameters())  # rollback to the honest one
    assert ctl.best["value"] == 0.80


def test_b2_2_alarm_without_a_best_checkpoint_does_not_crash_training(tmp_path):
    pytest.importorskip("torch")  # the Controller saves checkpoints with torch ([train] extra)
    t = fake_trainer(tmp_path)
    ctl = Controller(t, metric="eval_rewards/iou_tc/mean", alarm=lambda m: "duplicates")
    ctl.on_evaluate(t.args, t.state, None, metrics={"eval_rewards/advantage_v4/mean": 0.0})  # metric missing
    assert ctl.events and ctl.events[-1]["kind"] != "best"


def test_b2_3_malformed_control_file_does_not_crash_training(tmp_path):
    pytest.importorskip("torch")  # the Controller saves checkpoints with torch ([train] extra)
    t = fake_trainer(tmp_path)
    ctl = Controller(t)
    (tmp_path / "control.json").write_text('{"reward_weights": {"g_area": 0.5')  # editor mid-save
    ctl.on_step_begin(t.args, t.state, None)
