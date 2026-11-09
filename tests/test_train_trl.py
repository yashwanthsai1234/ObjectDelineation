"""TRL (version pinned in pyproject.toml) integration on a tiny offline model: R2 (TRL advantages = advantages_v4),
group contiguity per process
(world size 2 simulated with TRL's RepeatSampler and accelerate's sharding), callbacks (live guard switches,
best-checkpoint rollback, held-out evaluation, the hack audit on held-out completions with rollback, control-file
robustness), the completions table -> eval JSONL exporter and the checkpoint sampler (objdel.train.sample)."""
import json
import os

import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")

import pandas as pd  # noqa: E402
from accelerate.data_loader import prepare_data_loader  # noqa: E402
from datasets import Dataset  # noqa: E402
from tokenizers import Tokenizer, models, pre_tokenizers  # noqa: E402
from torch.utils.data import DataLoader  # noqa: E402
from transformers import LlamaConfig, LlamaForCausalLM, PreTrainedTokenizerFast, TrainerCallback  # noqa: E402
from trl import GRPOConfig, GRPOTrainer  # noqa: E402
from trl.trainer.utils import RepeatSampler  # noqa: E402

from objdel.reward import ExpertCache, advantages_v4, legal_gate, sig_id  # noqa: E402
from objdel.train import configs, sample  # noqa: E402
from objdel.train import rewards as R  # noqa: E402
from objdel.train.callbacks import Controller, audit_alarm  # noqa: E402

G = 8
A, B = sig_id(((0,),)), sig_id(((1,),))


def g(iou, f5, legal=(1,) * G, sig=(A,) * G, freq=(1.0,) * G, g_area=None):
    sig = np.asarray(sig, float)
    return {"format_ok": np.ones(G), "legal": np.asarray(legal, float), "iou_tc": np.asarray(iou, float),
            "f5": np.asarray(f5, float), "g_area": np.asarray(g_area if g_area is not None else f5, float) * 0.98,
            "topo_any": (np.asarray(freq) > 0).astype(float), "topo_freq": np.asarray(freq, float), "sig_id": sig}


HAND = {  # hand-set reward-v4 terms per prompt (one group each)
    0: g([.91, .88, .93, .90, .86, .92, .89, .95], [.97, .95, .98, .96, .93, .97, .94, .99],
         legal=[1, 1, 1, 0, 1, 1, 1, 1]),  # one illegal rollout
    1: g([.9] * G, [.95] * G),  # identical: dropped
    2: g([.85, .87, .83, .88, .86, .84, .80, .89], [.93, .94, .92, .95, .93, .92, .90, .96],
         sig=[A] * 6 + [B, A], freq=[2 / 3] * 6 + [1 / 3, 2 / 3]),  # a minority topology
    3: g([.900, .901, .903, .902, .904, .900, .905, .902],
         [.950, .951, .952, .950, .953, .951, .950, .952]),  # spread < 1e-2, one topology: dropped
    4: g([.70, .55, .82, .40, .65, .77, .30, .60], [.80, .70, .90, .50, .75, .85, .45, .72],
         sig=[A, A, A, B, A, A, B, A], freq=[1, 1, 1, 0, 1, 1, 0, 1]),  # two answers no rater drew
}


def fake_score_batch(jobs):
    """Stands in for objdel.reward.score_batch: the hand terms of each job's prompt (caches = {"key": prompt})."""
    seen, rows = {}, []
    for job in jobs:
        k = job[3]["key"]
        i = seen[k] = seen.get(k, -1) + 1
        rows.append({n: v[i % G] for n, v in HAND[k].items()})
    out = {n: np.array([r[n] for r in rows]) for n in rows[0]}
    return out | {"error": [""] * len(rows), "defects": [{}] * len(rows)}


@pytest.fixture(scope="module")
def tiny(tmp_path_factory):
    vocab = {w: i for i, w in enumerate(["<pad>", "<eos>", "prompt", *map(str, range(40))])}
    tk = Tokenizer(models.WordLevel(vocab, unk_token="<pad>"))
    tk.pre_tokenizer = pre_tokenizers.WhitespaceSplit()
    tok = PreTrainedTokenizerFast(tokenizer_object=tk, pad_token="<pad>", eos_token="<eos>", padding_side="left")
    torch.manual_seed(0)
    cfg = LlamaConfig(vocab_size=len(vocab), hidden_size=16, intermediate_size=32, num_hidden_layers=1,
                      num_attention_heads=2, num_key_value_heads=2, pad_token_id=0, eos_token_id=1, bos_token_id=1)
    path = tmp_path_factory.mktemp("tiny") / "model"
    LlamaForCausalLM(cfg).save_pretrained(path)  # TRL loads the KL reference model from the policy's path
    return path, tok


def make_trainer(tiny, out, monkeypatch, prompts=(0, 1, 2, 3, 4), steps=1, eval_prompts=None, **kw):
    monkeypatch.setattr(R, "score_batch", fake_score_batch)
    path, tok = tiny
    rw = R.TRLReward(lambda key: {"key": key}, stop_ids={tok.eos_token_id, tok.pad_token_id})
    funcs, weights = rw.funcs()
    row = lambda k: {"prompt": f"prompt {k}", "key": k, "W": 50, "H": 50, "names": ["nodule"]}
    args = dict(output_dir=str(out), max_steps=steps, per_device_train_batch_size=G * len(prompts),
                num_generations=G, max_completion_length=4, multi_objective_aggregation="sum_then_normalize",
                scale_rewards="none", loss_type="dr_grpo", reward_weights=weights, beta=1e-3, use_cpu=True,
                bf16=False, report_to="none", logging_steps=1, save_strategy="no", seed=0,
                remove_unused_columns=False, lr_scheduler_type="constant", learning_rate=1e-2) | kw
    if eval_prompts:
        args |= dict(eval_strategy="steps", eval_steps=1, per_device_eval_batch_size=G * len(eval_prompts))
    trainer = GRPOTrainer(model=LlamaForCausalLM.from_pretrained(path), reward_funcs=funcs, args=GRPOConfig(**args),
                          train_dataset=Dataset.from_list([row(k) for k in prompts]), processing_class=tok,
                          eval_dataset=Dataset.from_list([row(k) for k in eval_prompts]) if eval_prompts else None)
    captured = []
    orig = trainer._generate_and_score_completions

    def hook(inputs):
        out = orig(inputs)
        if trainer.model.training:
            captured.append(([x["key"] for x in inputs], out["advantages"].double().numpy().copy()))
        return out

    trainer._generate_and_score_completions = hook
    return trainer, captured


def terms_in_order(keys):
    T = fake_score_batch([(None, None, None, {"key": k}) for k in keys])
    return {n: v for n, v in T.items() if n not in ("error", "defects")}


def test_r2_trl_advantages_equal_advantages_v4(tiny, tmp_path, monkeypatch):
    trainer, captured = make_trainer(tiny, tmp_path, monkeypatch)
    trainer.train()
    ((keys, adv),) = captured
    blocks = [keys[i:i + G] for i in range(0, len(keys), G)]
    assert all(len(set(b)) == 1 for b in blocks) and sorted(b[0] for b in blocks) == [0, 1, 2, 3, 4]
    T = terms_in_order(keys)
    ref = advantages_v4(T)
    assert np.abs(adv - ref).max() <= 1e-6  # R2
    order = [b[0] for b in blocks]
    by_group = adv.reshape(-1, G)
    assert np.all(by_group[order.index(1)] == 0) and np.all(by_group[order.index(3)] == 0)  # dropped: exactly 0
    assert np.abs(by_group[order.index(4)]).max() > 0.5
    for name in ("format_ok", "legal", "iou_tc", "f5", "g_area"):  # per-term logs: the terms, legality-gated
        logged = np.array(list(trainer._logs["rewards"][name]))
        assert np.abs(logged - T[name] * (legal_gate(T) if name in R.GATED else 1)).max() <= 1e-6
    last = [h for h in trainer.state.log_history if "loss" in h][-1]
    assert {f"rewards/{n}/mean" for n in R.LOGGED} <= set(last)
    assert not any(isinstance(v, float) and np.isnan(v) for v in last.values())


def first_local_batch(cfg, world, rank, n_prompts=200):
    """What one process's reward call receives first: TRL's RepeatSampler stream, sharded by accelerate."""
    gbs = cfg.per_device_train_batch_size * world * cfg.steps_per_generation  # GRPOConfig at this world size
    sampler = RepeatSampler(range(n_prompts), mini_repeat_count=cfg.num_generations, batch_size=gbs // G,
                            repeat_count=cfg.num_iterations * cfg.steps_per_generation, shuffle=True, seed=cfg.seed)
    dl = DataLoader(range(n_prompts), batch_size=cfg.per_device_train_batch_size * cfg.steps_per_generation,
                    sampler=sampler, collate_fn=list)
    return next(iter(prepare_data_loader(dl, num_processes=world, process_index=rank, put_on_device=False)))


def call_reward_on(keys):
    rw = R.TRLReward(lambda key: {"key": key}, stop_ids={1})
    f = rw.funcs()[0][-1]
    n = len(keys)
    return f(prompts=[""] * n, completions=["x"] * n, completion_ids=[[1]] * n, key=keys, W=[50] * n, H=[50] * n,
             names=[["nodule"]] * n)


@pytest.mark.parametrize("model", list(configs.MODELS))
def test_every_process_receives_whole_groups(model, tmp_path, monkeypatch):
    monkeypatch.setattr(R, "score_batch", lambda jobs: fake_score_batch(
        [(None, None, None, {"key": k % 5}) for k in (j[3]["key"] for j in jobs)]))
    cfg = configs.grpo_config(model, str(tmp_path))
    for world in (1, 2, 4):
        local = [first_local_batch(cfg, world, rank) for rank in range(world)]
        prompts = [k for b in local for k in b[::G]]
        n_prompts = cfg.per_device_train_batch_size * world * cfg.steps_per_generation // G
        assert len(prompts) == len(set(prompts)) == n_prompts
        for keys in local:
            assert all(len(set(keys[i:i + G])) == 1 for i in range(0, len(keys), G))
            call_reward_on(keys)  # accepted


def test_a_batch_that_splits_groups_across_processes_is_refused(tmp_path):
    cfg = GRPOConfig(output_dir=str(tmp_path), per_device_train_batch_size=4, gradient_accumulation_steps=1,
                     num_generations=8, generation_batch_size=8,  # TRL accepts it: 8 completions over 2 processes
                     bf16=False)  # TRL's bf16 default raises on CPU-only machines (CI train job)
    cfg.steps_per_generation = 1
    keys = first_local_batch(cfg, 2, 0)
    assert len(keys) == 4
    with pytest.raises(ValueError, match="whole groups"):
        call_reward_on(keys)


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


# ------------------------------------------------------------------ hack audit, exporter, sampler (real caches)

def _squares(x0, side, n=G):
    """n masks of a square of `side` grid units at (x0 + 20i, x0); 20 grid units = 1 px on 50 px."""
    return [json.dumps({"nodule": [{"outer": [[a, x0], [a + side, x0], [a + side, x0 + side], [a, x0 + side]],
                                    "holes": []}]}, separators=(",", ":")) for a in range(x0, x0 + 20 * n, 20)]


def _caches():
    m = np.zeros((3, 50, 50), bool)
    for i, (a, b) in enumerate([(10, 30), (10, 30), (10, 32)]):
        m[i, a:b, a:b] = True  # rater areas 400, 400, 484 px
    c = {"nodule": ExpertCache.build(m)}
    return {k: c for k in ("lidc/t1", "lidc/t2", "lidc/v1", "lidc/v2")}


def scripted_trainer(char_model, out, policy, steps=2, **kw):
    """GRPOTrainer from the production grpo_config and TRLReward on real caches; generation is replaced by
    policy[0](key) -> G (text, ends with EOS) pairs, so the held-out completions are known exactly."""
    path, tok = char_model
    caches = _caches()
    row = lambda k: {"prompt": "x", "key": k, "W": 50, "H": 50, "names": ["nodule"]}
    train = Dataset.from_list([row(k) for k in ("lidc/t1", "lidc/t2")])
    val = Dataset.from_list([row(k) for k in ("lidc/v1", "lidc/v2")])
    model = LlamaForCausalLM.from_pretrained(path)
    funcs, _ = R.TRLReward(caches.__getitem__, R.stop_ids(model, tok)).funcs()
    args = configs.grpo_config("smolvlm-256m", str(out), **{
        "max_steps": steps, "per_device_train_batch_size": 2 * G, "gradient_accumulation_steps": 1,
        "per_device_eval_batch_size": 2 * G, "use_cpu": True, "bf16": False, "logging_steps": 1,
        "save_strategy": "no", "eval_strategy": "steps", "eval_steps": 1, "learning_rate": 1e-2} | kw)
    trainer = GRPOTrainer(model=model, reward_funcs=funcs, args=args, train_dataset=train, eval_dataset=val,
                          processing_class=tok)
    pending, orig = [], trainer._generate_and_score_completions

    def hook(inputs):
        pending[:] = [x["key"] for x in inputs]
        return orig(inputs)

    def generate(prompt_ids, images, multimodal_fields, num_generations, has_tool_images=False):
        ids, seen = [], {}
        for k in pending:
            i = seen[k] = seen.get(k, -1) + 1
            text, eos = policy[0](k)[i % G]
            ids.append(tok(text, add_special_tokens=False)["input_ids"] + ([tok.eos_token_id] if eos else []))
        return ids, None, None  # (completion_ids, logprobs, sampling_mask): TRL 1.17

    trainer._generate_and_score_completions = hook
    trainer._generate_single_turn = generate
    return trainer


HONEST = [(t, True) for t in _squares(160, 400)]   # the raters' 20 px square, shifted by -2 to +5 px
HEDGE = [(t, True) for t in _squares(60, 700)]     # 35 px: A / A_max = 2.5 (dilated union)


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


def test_completions_table_exports_to_eval_jsonl(char_model, tmp_path):
    texts = {"lidc/t1": HONEST, "lidc/t2": [(t, i % 3 != 0) for i, (t, _) in enumerate(HEDGE)],
             "lidc/v1": HONEST[::-1], "lidc/v2": HEDGE}
    trainer = scripted_trainer(char_model, tmp_path, [texts.__getitem__], steps=1)
    trainer.train()
    path = tmp_path / "completions" / "completions_00001.parquet"
    df = pd.read_parquet(path)
    assert {"key", "truncated", "group"} <= set(df.columns)  # B4-3
    rows = sample.export([path], keys=["lidc/v1", "lidc/v2"])  # this step's evaluation rewrote the table
    assert sorted(r["item_id"] for r in rows) == ["lidc/v1", "lidc/v2"]
    for r in rows:
        assert r["samples"] == [t for t, _ in texts[r["item_id"]]] and r["truncated"] == [False] * G
        assert r["step"] == 1
    sample.write(rows, tmp_path / "x.jsonl")
    assert [json.loads(line) for line in (tmp_path / "x.jsonl").read_text().splitlines()] == rows
    train_rows = sample.export([path], keys=["lidc/t2"])
    if train_rows:  # the step's training completions, when still in TRL's table
        assert train_rows[0]["truncated"] == [i % 3 == 0 for i in range(G)]


def test_export_refuses_a_table_that_is_not_whole_groups(tmp_path):
    df = pd.DataFrame({"step": [1] * G, "completion": ["a"] * G, "key": ["k"] * 4 + ["j"] * 4,
                       "truncated": [False] * G, "group": [0] * G})
    df.to_parquet(tmp_path / "c.parquet")
    with pytest.raises(ValueError, match="not one group"):
        sample.export([tmp_path / "c.parquet"])


def test_sampler_returns_k_completions_per_row_with_trls_truncation_rule(char_model, tmp_path, monkeypatch):
    script = {"lidc/v1": [(t, i != 2) for i, (t, _) in enumerate(HONEST)], "lidc/v2": HEDGE}
    trainer = scripted_trainer(char_model, tmp_path, [script.__getitem__], steps=1)
    tok, key = trainer.processing_class, [None]

    def generate(prompt_ids, *a, **kw):  # one call of 16 prompts (per_device_eval_batch_size)
        # (completion_ids, logprobs, sampling_mask): TRL 1.17
        return [tok(t, add_special_tokens=False)["input_ids"] + ([tok.eos_token_id] if e else [])
                for t, e in (script[key[0]] * 2)[:len(prompt_ids)]], None, None

    trainer._generate_single_turn = generate
    trainer._generate_and_score_completions = None  # the sampler neither scores nor trains
    mode, rows = trainer.model.training, []
    for key[0] in ("lidc/v1", "lidc/v2"):
        rows += sample.sample(trainer, [{"prompt": "x", "key": key[0]}], k=16)
    assert trainer.model.training == mode
    assert [r["item_id"] for r in rows] == ["lidc/v1", "lidc/v2"]
    assert rows[0]["samples"] == [t for t, _ in script["lidc/v1"]] * 2
    assert rows[0]["truncated"] == [i % G == 2 for i in range(16)] and rows[1]["truncated"] == [False] * 16
    with pytest.raises(ValueError, match="whole groups"):
        sample.sample(trainer, [{"prompt": "x", "key": "lidc/v1"}], k=12)
    stops = []  # the reward's stop ids (REVIEW #7): one rule for the sampler and the reward
    monkeypatch.setattr(R, "stop_ids", lambda model, proc: stops.append((model, proc)) or set())
    (row,) = sample.sample(trainer, [{"prompt": "x", "key": "lidc/v1"}], k=8)
    assert stops == [(trainer.model, trainer.processing_class)] and row["truncated"] == [True] * G


def test_sampler_generates_with_the_real_model(tiny, tmp_path):
    """Unscripted: TRL's HF generation of the tiny model (max_completion_length 4, so most samples truncate)."""
    path, tok = tiny
    args = GRPOConfig(output_dir=str(tmp_path), num_generations=G, per_device_train_batch_size=G,
                      per_device_eval_batch_size=G, max_completion_length=4, use_cpu=True, bf16=False,
                      report_to="none", remove_unused_columns=False)
    trainer = GRPOTrainer(model=LlamaForCausalLM.from_pretrained(path), reward_funcs=lambda completions, **_: [0.0] *
                          len(completions), args=args, processing_class=tok,
                          train_dataset=Dataset.from_list([{"prompt": "prompt 1"}]))
    ids, real = [], trainer._generate
    trainer._generate = lambda prompts: (out := real(prompts), ids.extend(out[1]))[0]
    (row,) = sample.sample(trainer, [{"prompt": "prompt 3", "key": "lidc/x"}], k=16)
    assert len(row["samples"]) == len(ids) == 16 and all(isinstance(t, str) for t in row["samples"])
    assert all(len(x) <= 4 for x in ids) and any(row["truncated"])
    assert row["truncated"] == [x[-1] not in R.stop_ids(trainer.model, tok) for x in ids]


def test_controller_refuses_groups_other_than_eight(char_model, tmp_path):
    trainer = scripted_trainer(char_model, tmp_path, [lambda k: HONEST], steps=1)
    trainer.num_generations_eval = 16  # e.g. a GRPOConfig built by hand with num_generations_eval=16
    with pytest.raises(ValueError, match="groups of"):
        Controller(trainer)


# ------------------------------------------------------------------ control file and alarms (no generation)

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
