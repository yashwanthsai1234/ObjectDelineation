"""Controller tests moved from tests/test_train_trl.py (CLEAN6, from commit f531f86). Archived: never collected, linted or
run. To run them, restore the Controller from the tag pre-clean-2026-11-08 (README.md)."""
import json
import os

import numpy as np
import pytest
import torch
from test_train_trl import (G, HEDGE, HONEST, _squares, advantages_v4, legal_gate, make_trainer,
                            scripted_trainer, terms_in_order, tiny)  # live helpers (tests/test_train_trl.py)
from transformers import TrainerCallback

from callbacks import Controller, audit_alarm  # archive/controller/callbacks.py


class At(TrainerCallback):
    """Runs `fn` at the end of optimizer step `step` (so it lands before the next step's generation)."""

    def __init__(self, step, fn):
        self.step, self.fn = step, fn

    def on_step_end(self, args, state, control, **kw):
        if state.global_step == self.step:
            self.fn()


@pytest.mark.parametrize("how", ["api", "control_file"])
def test_guard_switch_changes_the_next_rewards_and_beta(how, tiny, tmp_path, monkeypatch):
    trainer, captured = make_trainer(tiny, tmp_path, monkeypatch, steps=2)
    ctl = Controller(trainer)
    betas = []
    orig = trainer._compute_loss
    trainer._compute_loss = lambda *a, **k: betas.append(trainer.beta) or orig(*a, **k)
    if how == "api":
        switch = lambda: ctl.switch(reward_weights={"g_area": 2.0}, beta=0.5)
    else:
        switch = lambda: (tmp_path / "control.json").write_text(
            json.dumps({"reward_weights": {"g_area": 2.0}, "beta": 0.5}))
    trainer.add_callback(At(1, switch))
    trainer.add_callback(ctl)
    trainer.train()
    (k0, a0), (k1, a1) = captured
    assert np.abs(a0 - advantages_v4(terms_in_order(k0))).max() <= 1e-6
    T = terms_in_order(k1)
    ga = (T["g_area"] * legal_gate(T)).reshape(-1, G)  # the logger is legality-gated, as advantages_v4's terms
    expect = advantages_v4(T) + 2.0 * (ga - ga.mean(1, keepdims=True)).reshape(-1)
    assert np.abs(a1 - expect).max() <= 1e-5
    assert betas == [1e-3, 0.5]
    assert [e["kind"] for e in ctl.events] == ["switch"] and ctl.events[0]["step"] == 1


def test_raising_beta_without_a_reference_model_is_refused(tiny, tmp_path, monkeypatch):
    trainer, _ = make_trainer(tiny, tmp_path, monkeypatch, beta=0.0)
    with pytest.raises(ValueError, match="reference model"):
        Controller(trainer).switch(beta=0.1)


def params(model):
    return {n: p.detach().clone() for n, p in model.named_parameters() if p.requires_grad}


def test_rollback_restores_the_best_checkpoint(tiny, tmp_path, monkeypatch):
    trainer, _ = make_trainer(tiny, tmp_path, monkeypatch, steps=1)
    trainer.train()  # creates the optimizer state
    metric = "eval_rewards/iou_tc/mean"
    ctl = Controller(trainer, metric=metric, alarm=lambda m: "iou dropped" if m[metric] < 0.2 else None)
    args, state, control = trainer.args, trainer.state, trainer.control
    ctl.on_evaluate(args, state, control, metrics={metric: 0.5})
    best = params(trainer.model)
    with torch.no_grad():
        for p in trainer.model.parameters():
            p.add_(1.0)
    ctl.on_evaluate(args, state, control, metrics={metric: 0.3})  # worse, no alarm: kept, not saved
    assert any(not torch.equal(best[n], p) for n, p in params(trainer.model).items())
    assert ctl.best["value"] == 0.5
    assert len(trainer.optimizer.state) > 0
    ctl.on_evaluate(args, state, control, metrics={metric: 0.1})  # alarm: roll back
    assert all(torch.equal(best[n], p) for n, p in params(trainer.model).items())
    assert len(trainer.optimizer.state) == 0
    assert [e["kind"] for e in ctl.events] == ["best", "rollback"] and ctl.events[1]["reason"] == "iou dropped"


def test_held_out_evaluation_feeds_the_controller(tiny, tmp_path, monkeypatch):
    trainer, _ = make_trainer(tiny, tmp_path, monkeypatch, steps=1, eval_prompts=(2, 4))
    ctl = Controller(trainer, audit=False)  # hand-set terms, no caches to audit
    trainer.add_callback(ctl)
    trainer.train()
    assert ctl.best is not None and ctl.best["step"] == 1
    ev = [h for h in trainer.state.log_history if "eval_rewards/iou_tc/mean" in h]
    assert ev and ev[0]["eval_rewards/iou_tc/mean"] == pytest.approx(ctl.best["value"])
    assert os.path.isfile(ctl.best["path"])


def test_audit_alarm_rolls_back_a_hedging_policy_to_the_last_good_checkpoint(char_model, tmp_path):
    policy = [lambda k: [(t, True) for t in _squares(320, 400)]]  # SFT: honest, further off (lower IoU_tc)
    trainer = scripted_trainer(char_model, tmp_path, policy, eval_on_start=True)
    seen, snap = [], {}
    ctl = Controller(trainer, alarm=lambda m: seen.append(dict(m)) or audit_alarm(m))

    def at_step_end(step):
        if step == 1:
            snap.update(params(trainer.model))
            policy[0] = lambda k: HONEST  # the step-1 evaluation improves: best
        if step == 2:
            policy[0] = lambda k: HEDGE if k.startswith("lidc/v") else HONEST  # the step-2 evaluation hedges

    trainer.add_callback(At(1, lambda: at_step_end(1)))
    trainer.add_callback(At(2, lambda: at_step_end(2)))
    trainer.add_callback(ctl)
    trainer.train()
    assert len(seen) == 3  # evaluations at steps 0 (eval_on_start: the baseline), 1 and 2
    assert all(m["eval_audit/lidc/nodule/a_amax_median"] < 1.0 for m in seen[:2])
    assert seen[2]["eval_audit/lidc/nodule/a_amax_median"] == pytest.approx(2.5, abs=0.1)
    assert seen[0]["eval_audit/fired"] == seen[1]["eval_audit/fired"] == []
    assert {"area_median", "area_tail", "inflation"} <= set(seen[2]["eval_audit/fired"])
    assert ctl.baseline is not None and ctl.baseline["fired"] == []
    kinds = [(e["step"], e["kind"]) for e in ctl.events]
    assert (2, "rollback") in kinds and (2, "best") not in kinds
    assert ctl.best["step"] == 1 and ctl.events[-1]["to_step"] == 1
    assert all(torch.equal(snap[n], p) for n, p in params(trainer.model).items())  # the step-1 weights
    assert len(trainer.optimizer.state) == 0


def test_audit_alarms_that_fire_on_the_baseline_do_not_roll_back(char_model, tmp_path):
    """A check the SFT checkpoint already fails cannot be fixed by rolling back to it (e.g. issue I17)."""
    trainer = scripted_trainer(char_model, tmp_path, [lambda k: HEDGE], steps=1, eval_on_start=True)
    ctl = Controller(trainer)
    trainer.add_callback(ctl)
    trainer.train()
    audits = [e for e in ctl.events if e["kind"] == "audit"]
    assert len(audits) == 2 and "area_median" in audits[1]["fired"] and audits[1]["new"] == []
    assert not any(e["kind"] == "rollback" for e in ctl.events)


class Logs(TrainerCallback):
    """What report_to backends receive: every on_log call."""

    def __init__(self):
        self.seen = []

    def on_log(self, args, state, control, logs=None, **kw):
        self.seen.append(dict(logs))


def test_audit_statistics_reach_the_logs_and_controller_jsonl(char_model, tmp_path):
    """REVIEW #3: transformers logs the evaluation metrics before on_evaluate (trainer.py 5.20.1, evaluate:
    self.log(output.metrics), then callback_handler.on_evaluate), so the audit statistics must be logged by the
    controller itself and stored in its audit event."""
    trainer = scripted_trainer(char_model, tmp_path, [lambda k: HEDGE], steps=1, eval_on_start=True)
    ctl, logs = Controller(trainer), Logs()
    trainer.add_callback(ctl)
    trainer.add_callback(logs)
    trainer.train()
    stat = "eval_audit/lidc/nodule/a_amax_median"
    hist = [h for h in trainer.state.log_history if stat in h]
    assert [h["step"] for h in hist] == [0, 1]
    assert all(h[stat] == pytest.approx(2.5, abs=0.1) for h in hist)
    assert all(isinstance(v, int | float) for h in hist for k, v in h.items() if k.startswith("eval_audit/"))
    assert sum(stat in x for x in logs.seen) == 2  # on_log: what report_to backends get
    events = [json.loads(e) for e in (tmp_path / "controller.jsonl").read_text().splitlines()]
    audits = [e for e in events if e["kind"] == "audit"]
    assert [e["step"] for e in audits] == [0, 1]
    assert all(e["stats"][stat] == pytest.approx(2.5, abs=0.1) for e in audits)
    logged = {k: v for k, v in hist[1].items() if k.startswith("eval_audit/")}
    assert logged == {k: v for k, v in audits[1]["stats"].items() if v is not None}  # the same numbers


def _spiked(x0, side=400):
    """The square of _squares(x0, side, 1) with a collinear spike on its top edge: IoU_tc as the square, illegal."""
    m, y1 = x0 + side // 2, x0 + side
    return json.dumps({"nodule": [{"outer": [[x0, x0], [m, x0], [m, x0 - 100], [m, x0], [x0 + side, x0],
                                             [x0 + side, y1], [x0, y1]], "holes": []}]}, separators=(",", ":"))


def test_an_illegal_high_iou_evaluation_never_becomes_best(char_model, tmp_path):
    """R2 follow-up: the default best-checkpoint metric is legality-gated, so a policy whose held-out answers are
    illegal (spiked squares with IoU_tc 1) does not displace a legal, less accurate one."""
    policy = [lambda k: [(t, True) for t in _squares(320, 400)]]  # legal, off by 8 px
    trainer = scripted_trainer(char_model, tmp_path, policy, steps=1, eval_on_start=True)
    ctl = Controller(trainer, alarm=None)  # selection only
    trainer.add_callback(At(1, lambda: policy.__setitem__(0, lambda k: [(_spiked(200), True)] * G)))
    trainer.add_callback(ctl)
    trainer.train()
    evals = [h for h in trainer.state.log_history if "eval_rewards/iou_tc/mean" in h]
    assert [h["step"] for h in evals] == [0, 1]
    assert evals[1]["eval_rewards/legal/mean"] == 0.0 and evals[1]["eval_rewards/iou_tc/mean"] == 0.0
    assert ctl.best["step"] == 0 and not any(e["kind"] == "best" and e["step"] == 1 for e in ctl.events)


def test_controller_refuses_groups_other_than_eight(char_model, tmp_path):
    trainer = scripted_trainer(char_model, tmp_path, [lambda k: HONEST], steps=1)
    trainer.num_generations_eval = 16  # e.g. a GRPOConfig built by hand with num_generations_eval=16
    with pytest.raises(ValueError, match="groups of"):
        Controller(trainer)


def test_half_written_control_file_is_retried_and_bad_content_is_ignored(tiny, tmp_path, monkeypatch):
    trainer, _ = make_trainer(tiny, tmp_path, monkeypatch)
    ctl = Controller(trainer, audit=False)
    f = tmp_path / "control.json"
    args, state, control = trainer.args, trainer.state, trainer.control
    f.write_text('{"beta": 0.5')  # editor mid-save
    ctl.on_step_begin(args, state, control)
    ctl.on_step_begin(args, state, control)  # logged once
    assert [e["kind"] for e in ctl.events] == ["control_error"] and trainer.beta == 1e-3
    f.write_text('{"beta": 0.5}')
    os.utime(f, ns=(1, 1))  # same mtime as nothing seen before: the completed file is read
    ctl.on_step_begin(args, state, control)
    assert trainer.beta == 0.5 and ctl.events[-1]["kind"] == "switch"
    for bad in ('{"reward_weights": {"nope": 1}}', '{"beta": "high"}', "[1, 2]"):
        f.write_text(bad)
        os.utime(f, ns=(len(bad), len(bad)))
        ctl.on_step_begin(args, state, control)
        assert ctl.events[-1]["kind"] == "control_error" and trainer.beta == 0.5
    events = (tmp_path / "controller.jsonl").read_text().splitlines()
    assert [json.loads(e)["kind"] for e in events] == [e["kind"] for e in ctl.events]


def test_control_file_rollback_before_any_evaluation_restores_the_initial_weights(tiny, tmp_path, monkeypatch):
    trainer, _ = make_trainer(tiny, tmp_path, monkeypatch, steps=1)
    initial = params(trainer.model)
    ctl = Controller(trainer, audit=False)
    trainer.train()
    assert any(not torch.equal(initial[n], p) for n, p in params(trainer.model).items())
    (tmp_path / "control.json").write_text(json.dumps({"rollback": "by hand"}))
    ctl.on_step_begin(trainer.args, trainer.state, trainer.control)
    assert all(torch.equal(initial[n], p) for n, p in params(trainer.model).items())
    assert ctl.events[-1] | {"step": None} == {"step": None, "kind": "rollback", "to_step": 0, "reason": "by hand"}
