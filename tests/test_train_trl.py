"""TRL (version pinned in pyproject.toml) integration on a tiny offline model: R2 (TRL advantages = advantages_v4),
group contiguity per process
(world size 2 simulated with TRL's RepeatSampler and accelerate's sharding), the completions table -> eval JSONL
exporter and the checkpoint sampler (objdel.train.sample). The Controller tests are in archive/controller/tests/."""
import json

import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")

import pandas as pd  # noqa: E402
from accelerate.data_loader import prepare_data_loader  # noqa: E402
from datasets import Dataset  # noqa: E402
from tokenizers import Tokenizer, models, pre_tokenizers  # noqa: E402
from torch.utils.data import DataLoader  # noqa: E402
from transformers import LlamaConfig, LlamaForCausalLM, PreTrainedTokenizerFast  # noqa: E402
from trl import GRPOConfig, GRPOTrainer  # noqa: E402
from trl.trainer.utils import RepeatSampler  # noqa: E402

from objdel.reward import ExpertCache, advantages_v4, legal_gate, sig_id  # noqa: E402
from objdel.train import configs, sample  # noqa: E402
from objdel.train import rewards as R  # noqa: E402

G = 8
A, B = sig_id(((0,),)), sig_id(((1,),))


def g(iou, f5, legal=(1,) * G, sig=(A,) * G, freq=(1.0,) * G, g_area=None):
    sig = np.asarray(sig, float)
    return {"format_ok": np.ones(G), "legal": np.asarray(legal, float), "iou_tc": np.asarray(iou, float),
            "f5": np.asarray(f5, float), "g_area": np.asarray(g_area if g_area is not None else f5, float) * 0.98,
            "topo_any": (np.asarray(freq) > 0).astype(float), "topo_freq": np.asarray(freq, float), "sig_id": sig,
            "loc": np.zeros(G)}


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


# ------------------------------------------------------------------ exporter, sampler (real caches)

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
    args = configs.grpo_config("qwen3-vl-2b", str(out), use_vllm=False, **{  # HF generate on the tiny model
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


# ------------------------------------------------------------------ CLEAN6: the Controller is archived (CT-08)

def test_controller_is_gone_and_the_live_callbacks_stay():
    import re
    import subprocess
    from pathlib import Path

    import objdel.eval
    from objdel.train import callbacks

    assert callbacks.TraceStep and callbacks.NvtxSteps and callbacks.BudgetStop
    for name in ("Contr" + "oller", "audit_" + "alarm", "FI" + "RED"):
        with pytest.raises(AttributeError):
            getattr(callbacks, name)
    assert all(hasattr(objdel.eval, n) for n in ("hack_audit", "gates", "metrics", "summarize", "vision"))
    assert not hasattr(objdel.eval, "alarms" + "_for")
    root = Path(__file__).resolve().parents[1]
    files = subprocess.run(["git", "-C", str(root), "ls-files", "-z", "-co", "--exclude-standard", "src",
                            "experiments/phase2"], capture_output=True, text=True, check=True).stdout.split("\0")
    words = re.compile(r"\bContr" + r"oller\b|\balarms" + r"_for\b|\.recor" + r"ded\b")
    hits = [f for f in files if (f.startswith("src/") or f.count("/") == 2) and f.endswith(".py")
            and (root / f).is_file() and words.search((root / f).read_text())]
    assert hits == []
