"""Break-it B4: end-to-end wiring of PLAN §3.1/§3.2 and the build / CI / docs around it.

The unmarked end-to-end test runs the whole chain without data in ~15 s: a synthetic item store (objdel.data.store)
-> objdel.data.prepare.prepare (targets.frame, targets.build_target, reward.ExpertCache) -> grpo_rows (RIGA and LIDC
rows mixed) -> the production TRLReward on the prepared caches inside a real TRL GRPOTrainer (configs.grpo_config,
tiny offline Llama, char-level tokenizer, a scripted policy that emits rater targets, a truncated target, junk and
an empty answer) -> TRL's held-out evaluation -> the same completions through experiments/phase1/eval_run.py
(store frames, eval.metrics) -> eval.gates.phase2_entry.
Tests `b4_<n>` reproduce the confirmed breakdowns, now fixed (experiments/phase1/impl/breakit/B4_wiring.md).
"""
import ast
import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PREPARED = ROOT / "data" / "prepared"
G = 8


# ------------------------------------------------------------------ end to end, no data

def _script(rec):
    """Eight (text, ends with EOS) completions per prompt; z004 gets eight junk answers (an all-failure group)."""
    t = [x for x in rec["targets"] if x]
    if rec["item_id"].endswith("z004"):
        return [("{", True)] * G
    empty = '{"nodule":[]}' if rec["item_id"].startswith("lidc") else '{"disc":[],"cup":[]}'
    return [(t[0], True), (t[-1], True), (t[0], False), ("{", True), (t[0][:-3], True), (empty, True),
            (t[0], True), (t[-1], True)]


@pytest.fixture(scope="module")
def chain(tmp_path_factory, synthetic_items, char_model):
    pytest.importorskip("trl")
    pytest.importorskip("pyarrow")
    from datasets import concatenate_datasets
    from transformers import LlamaForCausalLM
    from trl import GRPOTrainer

    from objdel.data import prepare as data_prepare
    from objdel.data import store
    from objdel.policy import models as policy_models
    from objdel.policy import prompt as policy_prompt
    from objdel.policy import rows as policy_rows
    from objdel.rl import config as configs
    from objdel.rl import reward as rl_reward
    from objdel.rl import rows as rl_rows
    from objdel.verify import store as verify_store

    work = tmp_path_factory.mktemp("b4")
    store.write(synthetic_items(("train", "val"), unanimous=True), work / "store")
    data_prepare.prepare(work / "store", work / "prep", splits=("train", "val"), procs=1, prompts=policy_prompt.PROMPTS)
    recs = {r["item_id"]: r for d in ("riga", "lidc") for s in ("train",
        "val") for r in policy_rows.records(work / "prep", d, s)}
    path, tok = char_model
    model = LlamaForCausalLM.from_pretrained(path)

    def rows(split):  # the text-only tiny model takes no image: drop the column, keep every reward column
        ds = concatenate_datasets([rl_rows.grpo_rows(work / "prep", d, split) for d in ("riga", "lidc")])
        return ds.remove_columns("image").map(lambda r: {"prompt": "x"})

    train, val = rows("train"), rows("val")
    funcs, _ = rl_reward.TRLReward(verify_store.cache_loader(work / "prep"), policy_models.stop_ids(model, tok)).funcs()
    args = configs.grpo_config("qwen3-vl-2b", str(work / "run"), use_vllm=False, learning_rate=1e-6, max_steps=2,
                               per_device_train_batch_size=G * len(train), gradient_accumulation_steps=1,
                               per_device_eval_batch_size=G * len(val), use_cpu=True, bf16=False,
                               logging_steps=1, save_strategy="no", eval_strategy="steps", eval_steps=1)
    trainer = GRPOTrainer(model=model, reward_funcs=funcs, args=args, train_dataset=train, eval_dataset=val,
                          processing_class=tok)
    script = {k: _script(r) for k, r in recs.items() if r["rl"]}
    pending, calls = [], []
    orig = trainer._generate_and_score_completions

    def hook(inputs):
        pending[:] = [x["key"] for x in inputs]
        out = orig(inputs)
        calls.append({"train": trainer.model.training, "keys": list(pending),
                      "adv": out["advantages"].double().cpu().numpy().copy()})
        return out

    def scripted(prompt_ids, images, multimodal_fields, num_generations, has_tool_images=False):
        ids, seen = [], {}
        for k in pending:
            i = seen[k] = seen.get(k, -1) + 1
            text, eos = script[k][i % G]
            ids.append(tok(text, add_special_tokens=False)["input_ids"] + ([tok.eos_token_id] if eos else []))
        return ids, None, None  # (completion_ids, logprobs, sampling_mask): TRL 1.17

    trainer._generate_and_score_completions = hook
    trainer._generate_single_turn = scripted
    trainer.train()
    return {"work": work, "recs": recs, "script": script, "calls": calls, "trainer": trainer,
            "val_keys": list(val["key"])}


def _reference_advantages(chain, keys):
    """advantages_v4 of the scripted texts, scored independently of TRL on the prepared caches."""
    from objdel.geom.text import SCHEMAS
    from objdel.verify import store as verify_store
    from objdel.verify.reward.advantages import advantages_v4
    from objdel.verify.reward.score import score_batch
    load, rec, seen, jobs = verify_store.cache_loader(chain["work"] / "prep"), chain["recs"], {}, []
    for k in keys:
        i = seen[k] = seen.get(k, -1) + 1
        text, eos = chain["script"][k][i % G]
        jobs.append((text, not eos, SCHEMAS[k.split("/")[0]], load(k), rec[k]["W"], rec[k]["H"]))
    return advantages_v4(score_batch(jobs)), score_batch(jobs)


def test_e2e_store_to_trl_to_eval_and_gates(chain):
    from experiments.phase1 import eval_run
    train_calls = [c for c in chain["calls"] if c["train"]]
    eval_calls = [c for c in chain["calls"] if not c["train"]]
    assert len(train_calls) == 2 and len(eval_calls) == 2
    for c in chain["calls"]:  # R2 on the production path: TRL advantages = advantages_v4 of the real reward
        ref, T = _reference_advantages(chain, c["keys"])
        assert len(c["adv"]) == len(c["keys"])  # one reward per completion
        assert np.abs(c["adv"] - ref).max() <= 1e-5
        by_key = {k: c["adv"][i:i + G] for i, k in zip(range(0, len(c["keys"]), G), c["keys"][::G], strict=True)}
        assert {k.split("/")[0] for k in by_key} == {"riga", "lidc"}  # mixed datasets in one reward call
        junk = [k for k in by_key if k.endswith("z004")]
        assert junk and all(np.all(by_key[k] == 0.0) for k in junk)  # an all-failure group gives exactly 0
        for i, k in enumerate(c["keys"]):
            text, eos = chain["script"][k][i % G]
            if eos and text in chain["recs"][k]["targets"]:  # a rater's own target: legal, perfect, a rater mode
                assert T["legal"][i] == 1 and T["iou_tc"][i] >= 0.99 and T["topo_any"][i] == 1
            if not eos:
                assert T["error"][i] == "truncated" and T["format_ok"][i] == 0
    logs = chain["trainer"].state.log_history
    train_log = [h for h in logs if "loss" in h][-1]
    assert {"error/truncated", "error/incomplete", "rewards/advantage_v4/mean", "rewards/iou_tc/mean"} <= set(train_log)
    assert not any(isinstance(v, float) and np.isnan(v) for h in logs for v in h.values())
    evals = [h for h in logs if "eval_rewards/iou_tc/mean" in h]
    assert [h["step"] for h in evals] == [1, 2]  # TRL's held-out evaluation at every step
    # the same held-out completions through the eval runner: store frames, eval.metrics, report, gates
    rows = [{"item_id": k, "samples": [t for t, _ in chain["script"][k]],
             "truncated": [not e for _, e in chain["script"][k]]} for k in chain["val_keys"]]
    frames = eval_run.load_frames(chain["work"] / "store", "val", [r["item_id"] for r in rows])
    results = eval_run.evaluate(rows, frames)
    trl_eval = {k: eval_calls[0]["adv"][i:i + G] for i, k in zip(range(0, len(eval_calls[0]["keys"]), G),
                                                                  eval_calls[0]["keys"][::G], strict=True)}
    for r in results:  # eval harness advantages = the advantages TRL trained on (float32)
        assert np.abs(np.array(r["samples"]["advantage"]) - trl_eval[r["item_id"]]).max() <= 1e-5
    report = eval_run.build_report(results)
    json.dumps(report, allow_nan=False)
    assert set(report["gates"]["gates"]) == {"G-format", "G-geometry", "G-prior-b", "G-precision", "G-topology-a",
                                             "G-signal", "G-vision"}  # prereg.C1_GATES (plan §5.1)
    assert report["gates"]["gates"]["G-format"]["pass"] is False
    assert report["gates"]["pass"] is False  # junk and truncated samples: format 62.5 %
    assert set(report["summary"]["datasets"]) == {"riga", "lidc"}
    assert eval_run.markdown(report, "b4").count("| riga") >= 2


def test_b4_3_logged_completions_carry_item_key_and_truncation(chain):
    import pandas as pd
    files = sorted((chain["work"] / "run" / "completions").glob("completions_*.parquet"))
    assert files
    cols = set(pd.read_parquet(files[-1]).columns)
    assert "key" in cols and "truncated" in cols




def test_eval_runner_edges(chain):
    """Empty JSONL -> a failing report, not a crash; unknown ids -> KeyError naming them; an item without pixels
    is still scorable (frames need masks only)."""
    from experiments.phase1 import eval_run
    store_root = chain["work"] / "store"
    report = eval_run.build_report(eval_run.evaluate([], eval_run.load_frames(store_root, "val", [])))
    assert report["gates"]["pass"] is False and report["summary"]["datasets"] == {}
    assert eval_run.markdown(report, "empty")
    with pytest.raises(KeyError, match="not in the val store"):
        eval_run.load_frames(store_root, "val", ["lidc/9.9.9/n000/z001"])
    k = "lidc/1.2.2/n000/z002"
    res = eval_run.evaluate([{"item_id": k, "samples": ['{"nodule":[]}', "{"]}],
                            eval_run.load_frames(store_root, "val", [k]))
    assert res[0]["samples"]["format_ok"] == [1.0, 0.0]


# ------------------------------------------------------------------ gates, docs, scripts, CI (static)

def _summary(ds):
    """The part of objdel.verify.metrics.summarize's per-dataset output that the gates read, all perfect."""
    obj = {"iou_tc_median": 1.0}
    objects = {"riga": {"disc": obj, "cup": obj}, "lidc": {"nodule": obj}}[ds]
    return {"format": 1.0, "legal": 1.0, "iou_tc_pass": 1.0, "group": G, "dup_rate": 0.0, "signal_share": 1.0,
            "objects": objects}


def _report(datasets):
    """A perfect plan §5.1 report (schema: objdel.verify.gates): shifted + rescaled views, prior, probes, G-vision."""
    summary = {"datasets": {d: _summary(d) for d in datasets}}
    names = {"riga": ("disc", "cup"), "lidc": ("nodule",)}
    prior = {d: {n: {"diff": 0.1, "ci": [0.05, 0.15]} for n in names[d]} for d in datasets}
    return {"views": {v: {"summary": summary} for v in ("shifted", "rescaled")},
            "prior": {"shifted": prior, "rescaled": prior},
            "probes": {"precision": {"pass": True}, "topology": {"pass": True}},
            "vision": {"ratio": 1.0, "tv": 0.01, "pass": True}, "g_probe": {"verdict": "PASS"}}


def test_gates_pass_a_perfect_two_dataset_report():
    from objdel.verify import gates
    g = gates.phase2_entry(_report(("riga", "lidc")))
    assert g["pass"] and g["datasets"] == {"lidc": True, "riga": True} and g["warnings"] == []


def test_b4_5_gates_accept_a_model_dropped_from_lidc():
    from objdel.verify import gates
    r = _report(("riga",))
    r["probes"] = {"precision": {"pass": False, "capacity": {"pass": False}}, "topology": {"pass": False}}
    assert gates.phase2_entry(r)["pass"]  # LIDC (and its probes) not judged


def test_b4_6_readme_quick_start_imports_exist():
    import importlib
    code = "\n".join(re.findall(r"```python\n(.*?)```", (ROOT / "README.md").read_text(), re.S))
    for node in ast.walk(ast.parse(code)):
        if isinstance(node, ast.ImportFrom) and node.module.startswith("objdel"):
            mod = importlib.import_module(node.module)
            for a in node.names:
                assert hasattr(mod, a.name), f"{node.module}.{a.name}"


@pytest.mark.parametrize("script", ["fetch_lidc_pixels"])
def test_b4_7_scripts_without_argparse(script):
    """Every script: tests/test_experiment_scripts.py."""
    r = subprocess.run([sys.executable, "-I", "-m", f"experiments.phase1.{script}", "--help"], cwd=ROOT,
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0 and "usage" in r.stdout


WORKFLOW = (ROOT / ".github" / "workflows" / "tests.yml").read_text()
PYPROJECT = tomllib.loads((ROOT / "pyproject.toml").read_text())


def test_b4_1_setup_uv_is_pinned_to_an_existing_ref():
    refs = re.findall(r"uses:\s*astral-sh/setup-uv@(\S+)", WORKFLOW)
    assert refs and all(re.fullmatch(r"v\d+\.\d+\.\d+|[0-9a-f]{40}", r) for r in refs)


def _ci_requirements() -> tuple[set[str], set[str]]:
    """(distribution names installed by the CI jobs' `uv pip install -e ".[...]"`, the extras)."""
    extras = set()
    for spec in re.findall(r'uv pip install -e "\.\[([^\]]*)\]"', WORKFLOW):
        extras |= {e.strip() for e in spec.split(",")}
    optional = PYPROJECT["project"]["optional-dependencies"]
    reqs = PYPROJECT["project"]["dependencies"] + [r for e in extras for r in optional[e]]
    return {re.split(r"[<>=!~\[ ;]", r, maxsplit=1)[0].lower() for r in reqs}, extras


def test_objdel_data_imports_pyarrow_at_module_level():
    tree = ast.parse((ROOT / "src" / "objdel" / "data" / "store.py").read_text())
    top = {a.name.split(".")[0] for n in tree.body if isinstance(n, ast.Import) for a in n.names}
    assert "pyarrow" in top  # and objdel.data.targets.frame and objdel.data.audit (via eval_run) import objdel.data


def test_b4_2_ci_installs_what_objdel_data_imports():
    reqs, _ = _ci_requirements()
    assert "pyarrow" in reqs


def test_b4_8_ci_runs_the_trl_integration_tests():
    reqs, _ = _ci_requirements()
    assert {"torch", "trl", "datasets", "peft"} <= reqs


# ------------------------------------------------------------------ real data: a small sample through the chain

@pytest.mark.data
def test_real_lidc_sample_prepared_equals_store_and_eval_frames(store_dir):
    """store.load -> frame -> build_target / ExpertCache == data/prepared records and caches == eval_run frames, and
    each rater target scores as a perfect, legal answer on both cache sets."""
    pytest.importorskip("datasets")
    if not (PREPARED / "manifest.json").is_file():
        pytest.skip("data/prepared absent")
    from experiments.phase1 import eval_run
    from objdel.data import store
    from objdel.data.targets.build import build_target
    from objdel.data.targets.frame import frame
    from objdel.geom.text import SCHEMAS
    from objdel.geom.topology import min_area
    from objdel.policy import rows as policy_rows
    from objdel.verify import store as verify_store
    from objdel.verify.reward.cache import ExpertCache
    from objdel.verify.reward.score import score_completion
    recs = {r["item_id"]: r for r in policy_rows.records(PREPARED, "lidc", "val")}
    load, names, seen = verify_store.cache_loader(PREPARED, maxsize=4), SCHEMAS["lidc"], []
    for it in store.load("lidc", "val", store_dir):  # the first 4 items and the first disagreement item
        rec = recs[it.item_id]
        if len(seen) >= 4 and not rec["disagreement"]:
            continue
        f = frame(it)
        assert (f.W, f.H, list(f.box), f.scale) == (rec["W"], rec["H"], rec["box"], rec["scale"])
        fresh = {n: ExpertCache.build(f.masks[n]) for n in names}
        assert {n: c.digest for n, c in fresh.items()} == rec["caches"]
        ma = {n: min_area(int(f.masks[n].any(0).sum())) for n in names}
        for r, text in enumerate(rec["targets"]):
            assert text == build_target({n: f.masks[n][r] for n in names}, names, None, ma).text
            for caches in (load(it.item_id), fresh):
                s = score_completion(text, False, names, caches, f.W, f.H)
                assert s["legal"] == 1 and s["topo_any"] == 1 and s["iou_tc"] >= 0.99
        seen.append(it.item_id)
        if len(seen) == 5:
            break
    assert len(seen) == 5 and recs[seen[-1]]["disagreement"]
    frames = eval_run.load_frames(store_dir, "val", seen)
    for k in seen:
        assert all(np.array_equal(frames[k]["masks"][n], load(k)[n].raters) for n in names)
