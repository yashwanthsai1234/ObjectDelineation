"""experiments.phase2.train_rl: the full GRPO run. It trains for --steps, saves a resumable checkpoint every
--save-steps, resumes from the last one (a Vast box can die mid-run), and leaves the tuned model as a plain merged
model folder (the layout train_sft's branches use, so select_ckpt and vLLM read it unchanged)."""

import json
import os
import shutil
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from test_p2_scripts import GRPO_SMOKE, ROOT, _module, _tiny_vlm


@pytest.mark.slow
def test_full_rl_saves_resumes_and_leaves_a_merged_tuned_model(tmp_path, synthetic_items):
    pytest.importorskip("trl")
    from objdel.data import store
    from objdel.train import dataset

    model = _tiny_vlm(tmp_path / "tiny", "qwen3-vl-2b")  # LoRA (configs.lora), like the real run
    store.write(synthetic_items(("train", "val"), unanimous=True), tmp_path / "store")
    dataset.prepare(tmp_path / "store", tmp_path / "prepared", procs=1)
    out, env = tmp_path / "rl", os.environ | {"HF_HUB_OFFLINE": "1", "CUDA_VISIBLE_DEVICES": ""}

    def run(steps, *extra):
        r = _module("experiments.phase2.train_rl", "--model", "qwen3-vl-2b", "--ckpt", str(model), "--train",
                    str(tmp_path / "prepared"), "--out", str(out), "--steps", str(steps), "--save-steps", "1",
                    "--grpo", json.dumps(GRPO_SMOKE), *extra, env=env)
        assert r.returncode == 0, r.stderr[-4000:]
        return json.loads((out / "rl.json").read_text())

    doc = run(2)
    for k in (1, 2):  # resumable: adapter, optimizer, scheduler, RNG and trainer state
        files = {f.name for f in (out / f"trainer/checkpoint-{k}").iterdir()}
        assert {"adapter_model.safetensors", "optimizer.pt", "scheduler.pt", "trainer_state.json"} <= files, files
    assert doc["steps"] == 2 and doc["resumed_from"] is None
    assert [h["step"] for h in doc["log"] if "reward" in h] == [1, 2]

    doc = run(3, "--resume")  # continues at step 3, not from step 1
    assert doc["resumed_from"].endswith("checkpoint-2") and doc["steps"] == 3
    assert [h["step"] for h in doc["log"] if "reward" in h] == [1, 2, 3]
    state = json.loads((out / "trainer/checkpoint-3/trainer_state.json").read_text())
    assert state["global_step"] == 3

    final = {f.name for f in (out / "final").iterdir()}  # merged: a plain model folder, no adapter files
    assert {"model.safetensors", "config.json", "processor_config.json", "tokenizer.json"} <= final, final
    assert not {f for f in final if f.startswith("adapter_")}, final
    from transformers import AutoModelForImageTextToText

    assert AutoModelForImageTextToText.from_pretrained(out / "final") is not None


@pytest.mark.slow
def test_full_rl_runs_with_the_shared_prefix_forward(tmp_path, synthetic_items):
    """--shared-prefix: TRL's real loop (generation, old/ref passes, update micro-batches of whole groups) on the
    packed forward, to a merged tuned model."""
    pytest.importorskip("trl")
    from objdel.data import store
    from objdel.train import dataset

    model = _tiny_vlm(tmp_path / "tiny", "qwen3-vl-2b")
    store.write(synthetic_items(("train", "val"), unanimous=True), tmp_path / "store")
    dataset.prepare(tmp_path / "store", tmp_path / "prepared", procs=1)
    out, env = tmp_path / "rl", os.environ | {"HF_HUB_OFFLINE": "1", "CUDA_VISIBLE_DEVICES": ""}
    r = _module("experiments.phase2.train_rl", "--model", "qwen3-vl-2b", "--ckpt", str(model), "--train",
                str(tmp_path / "prepared"), "--out", str(out), "--steps", "2", "--save-steps", "100",
                "--grpo", json.dumps(GRPO_SMOKE), "--shared-prefix", env=env)  # fmt: skip
    assert r.returncode == 0, r.stderr[-4000:]
    doc = json.loads((out / "rl.json").read_text())
    steps = [h for h in doc["log"] if "reward" in h]
    assert [h["step"] for h in steps] == [1, 2]
    assert all(0 < h["shared/token_ratio"] < 1 for h in steps)  # the packed path ran, on fewer tokens
    assert (out / "final/model.safetensors").exists()


SPY = """
import json, sys
from objdel.data import splits
from objdel.train import rewards
splits.val_half = lambda g: "shape"  # the held-out curve's val-shape rule on synthetic patients
orig = rewards.TRLReward._terms
def spy(self, completions, completion_ids, key, *a, **k):
    if self._last[0] is not completions:
        with open(sys.argv[1], "a") as f:
            f.write(json.dumps(sorted(set(key))) + "\\n")
    return orig(self, completions, completion_ids, key, *a, **k)
rewards.TRLReward._terms = spy
from experiments.phase2 import train_rl
sys.exit(train_rl.main(sys.argv[2:]))
"""


@pytest.fixture
def rows_run(tmp_path, synthetic_items):
    """run(out, *flags) -> (CompletedProcess, keys per reward call): train_rl on a 4-step rows root (one prompt per
    step, blocks of 2; test_train_rl_guards.schedule_root) from a tiny LoRA Qwen3-VL, the reward calls' keys
    recorded."""
    pytest.importorskip("trl")
    from test_train_rl_guards import SMOKE, schedule_root

    from objdel.data import store
    from objdel.train import dataset

    model = _tiny_vlm(tmp_path / "tiny", "qwen3-vl-2b")
    store.write(synthetic_items(("train", "val"), unanimous=True), tmp_path / "store")
    dataset.prepare(tmp_path / "store", tmp_path / "prepared", procs=1)
    root, keys, sha = schedule_root(tmp_path / "prepared", tmp_path / "rows")
    env = os.environ | {"HF_HUB_OFFLINE": "1", "CUDA_VISIBLE_DEVICES": ""}

    n = iter(range(100))

    def run(out, *flags, train=True):
        log = out.parent / f"{out.name}.{next(n)}.keys"  # one file per call
        args = ["--model", "qwen3-vl-2b", "--ckpt", str(model), "--out", str(out),
                "--grpo", json.dumps(SMOKE), *flags]
        if train:
            args += ["--train", str(root), "--rows-revision", "a" * 40, "--rows-manifest-sha256", sha]
        r = subprocess.run([sys.executable, "-I", "-c", SPY, str(log), *args], cwd=ROOT, env=env,
                           capture_output=True, text=True, timeout=900)
        calls = [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []
        return r, calls

    return SimpleNamespace(run=run, keys=keys, sha=sha, root=root, model=model, prepared=tmp_path / "prepared")


@pytest.mark.slow
def test_rows_run_trains_in_schedule_order_and_records_provenance(tmp_path, rows_run):
    """I15, I39 (B1 side), I41 (adapter before merge), I49: a rows run trains the schedule's keys in file order,
    rl.json holds the rows provenance, and --merge-last rebuilds final/ from the last checkpoint."""
    out = tmp_path / "b1"
    r, calls = rows_run.run(out, "--steps", "2", "--save-steps", "1")
    assert r.returncode == 0, r.stderr[-4000:]
    assert calls == [[k] for k in rows_run.keys[:2]]  # step s trains schedule row s
    doc = json.loads((out / "rl.json").read_text())
    assert doc["rows"] == {"root": str(rows_run.root), "revision": "a" * 40, "manifest_sha256": rows_run.sha,
                           "ids_sha256": json.loads((rows_run.root / "manifest.json").read_text())["schedule"]
                           ["ids_sha256"]}
    assert doc["numerics"] == "bf16-seqmask" and doc["beta"] == 0.0 and doc["lr"] == 5e-5
    assert doc["mask_truncated_completions"] is True and doc["steps"] == 2
    assert (out / "final_adapter/adapter_config.json").exists()  # the unmerged LoRA, before the merge
    assert not [f for f in os.listdir(out / "final") if f.startswith("adapter_")]

    (out / "final").rename(out / "final_ref")  # a run that died after checkpoint-2, before its own save
    (out / "rl.json").unlink()
    shutil.rmtree(out / "final_adapter")
    r, _ = rows_run.run(out, "--merge-last", train=False)
    assert r.returncode == 0, r.stderr[-4000:]
    doc = json.loads((out / "rl.json").read_text())
    assert doc["steps"] == 2 and doc["partial"] is True and doc["merged_from"].endswith("checkpoint-2")
    assert (out / "final_adapter/adapter_config.json").exists()
    from safetensors.numpy import load_file

    got, ref = (load_file(out / d / "model.safetensors") for d in ("final", "final_ref"))
    assert got.keys() == ref.keys() and all(np.array_equal(got[k], ref[k]) for k in ref)
    stamp = (out / "final/model.safetensors").stat().st_mtime_ns
    r, _ = rows_run.run(out, "--merge-last", train=False)  # final/ exists: nothing changes
    assert r.returncode == 0 and (out / "final/model.safetensors").stat().st_mtime_ns == stamp
    r, _ = rows_run.run(tmp_path / "empty", "--merge-last", train=False)
    assert r.returncode != 0 and "no checkpoint" in r.stderr


@pytest.mark.slow
def test_resume_keeps_the_schedule_order(tmp_path, rows_run):
    """I33: 4 steps stopped after step 2 and resumed train the same keys per step as an unbroken run, and the
    resumed run checks ids_sha256 again."""
    r, whole = rows_run.run(tmp_path / "whole", "--steps", "4", "--save-steps", "2")
    assert r.returncode == 0, r.stderr[-4000:]
    assert whole == [[k] for k in rows_run.keys]
    out = tmp_path / "broken"
    r, first = rows_run.run(out, "--steps", "2", "--save-steps", "2")
    assert r.returncode == 0, r.stderr[-4000:]
    r, rest = rows_run.run(out, "--steps", "4", "--save-steps", "2", "--resume")
    assert r.returncode == 0, r.stderr[-4000:]
    assert first + rest == whole
    assert json.loads((out / "rl.json").read_text())["resumed_from"].endswith("checkpoint-2")
    m = json.loads((rows_run.root / "manifest.json").read_text())
    m["schedule"]["ids_sha256"] = "0" * 64
    (rows_run.root / "manifest.json").write_text(json.dumps(m))
    r, _ = rows_run.run(out, "--steps", "4", "--resume")  # pinned SHA no longer matches: refused before training
    assert r.returncode != 0 and "SHA-256" in r.stderr


@pytest.mark.slow
def test_heldout_curve_logs_eval_iou_and_allmiss(tmp_path, rows_run):
    """I44: with --heldout, FastGRPO runs TRL's evaluation on val-shape items of the view every --heldout-every
    steps and rl.json holds eval/iou_tc and eval/allmiss_share per evaluation."""
    out = tmp_path / "h"
    r, calls = rows_run.run(out, "--steps", "2", "--save-steps", "100", "--heldout", str(rows_run.prepared),
                            "--heldout-items", "64", "--heldout-every", "2")
    assert r.returncode == 0, r.stderr[-4000:]
    doc = json.loads((out / "rl.json").read_text())
    assert [e["step"] for e in doc["eval"]] == [2]
    assert 0 <= doc["eval"][0]["eval/iou_tc"] <= 1 and 0 <= doc["eval"][0]["eval/allmiss_share"] <= 1
    val = ["lidc/1.2.2/n000/z001", "lidc/1.2.2/n000/z004"]
    assert sorted({k for c in calls for k in c if k in val}) == val  # both val RL items were scored


# --- CLEAN2: no scoring pool in train_rl; gc.freeze once after the last preload (CT-33 to CT-36) -------------------

OPTIONS = ["-h", "--help", "--model", "--ckpt", "--train", "--out", "--steps", "--save-steps", "--shared-prefix",
           "--lr", "--resume", "--merge-last", "--rows-revision", "--rows-manifest-sha256", "--budget-usd",
           "--session-start", "--rate-usd-per-h", "--fixed-usd", "--heldout", "--heldout-items", "--heldout-every",
           "--procs", "--grpo", "--profile-window", "--trace"]  # pre-CLEAN literal (ap._actions at 08e6455)
LANE_OPTIONS = {"--shared-prefix": ["--fast-sync", "--group-requests", "--no-ckpt-layers"],  # C5, C6, C9
                "--grpo": ["--runaway-cap", "--reward", "--numerics", "--opsm-start-delta"],  # C11, B3, B2.2
                "--trace": ["--count-syncs", "--host-meta", "--check-meta-every",  # C0, C4
                            "--gpu-pixels", "--allow-cpu-pixels"]}  # C1
RL_KEYS = ["beta", "budget_stopped", "ckpt", "eval", "final", "log", "lr", "mask_truncated_completions", "model",
           "numerics", "resumed_from", "rows", "steps", "timing_s", "train_rows"]  # pre-CLEAN literal, same stub
LANE_RL_KEYS = ["loc_guard_stopped", "loss_scale", "max_completion_length", "opsm_start_delta", "reward"]  # B2.2, B3


@pytest.fixture
def stub(tmp_path, monkeypatch):
    """train_rl.main in process with the preload, the trainer and every way to start a process pool stubbed:
    run(*extra) -> exit code; calls records preload / freeze / train in order."""
    import concurrent.futures
    import gc
    import multiprocessing
    import multiprocessing.pool

    import experiments.phase2.g_probe as g_probe
    from experiments.phase2 import train_rl
    from objdel.train import dataset, serve

    calls = []

    def boom(*a, **k):
        raise AssertionError("train_rl started a process pool")

    for owner, name in [(multiprocessing, "Pool"), (multiprocessing.pool, "Pool"), (multiprocessing, "get_context"),
                        (concurrent.futures, "ProcessPoolExecutor"), (os, "fork"), (serve, "fork_pool")]:
        monkeypatch.setattr(owner, name, boom)
    monkeypatch.setattr(serve, "preload", lambda root, ds: calls.append("preload"))
    monkeypatch.setattr(gc, "freeze", lambda: calls.append("freeze"))
    monkeypatch.setattr(dataset, "grpo_rows", lambda *a, **k: ["k1", "k2"])
    monkeypatch.setattr(train_rl, "heldout_rows", lambda view, n: ["h1"])  # val_half lands with L2 (B1.3)
    saver = SimpleNamespace(save_pretrained=lambda p: None)
    trainer = SimpleNamespace(train=lambda resume_from_checkpoint=None: calls.append("train"), model=saver,
                              accelerator=SimpleNamespace(unwrap_model=lambda m: m),
                              state=SimpleNamespace(global_step=2, log_history=[{"step": 1, "reward": 0.5}]),
                              args=SimpleNamespace(beta=0.0, mask_truncated_completions=True,
                                                off_policy_mask_threshold=None, max_completion_length=1536))
    monkeypatch.setattr(g_probe, "make_trainer", lambda *a, **k: (trainer, saver, 5e-5))
    (tmp_path / "root").mkdir()
    (tmp_path / "root/manifest.json").write_text("{}")
    base = ["--model", "qwen3-vl-2b", "--ckpt", str(tmp_path / "c"), "--out", str(tmp_path / "out")]
    train = ["--train", str(tmp_path / "root"), "--steps", "2"]
    return SimpleNamespace(run=lambda *extra: train_rl.main([*base, *extra]), calls=calls, base=base, train=train,
                           out=tmp_path / "out", main=train_rl.main)


@pytest.mark.parametrize("extra", [["--procs", "2"], ["--procs=2"], ["--proc", "2"], ["--pro", "2"], ["--procs"],
                                   ["-procs", "2"]])
def test_procs_is_an_argparse_error_before_any_preload(stub, capsys, extra):
    """CT-33: every spelling of the old pool flag (an abbreviation included) exits 2 before a preload or a freeze."""
    with pytest.raises(SystemExit) as e:
        stub.run(*stub.train, *extra)
    assert e.value.code == 2 and "unrecognized arguments" in capsys.readouterr().err
    assert stub.calls == []


def test_help_has_no_procs_and_keeps_the_other_options(stub, capsys, monkeypatch):
    """CT-33: --help names no pool; the option list is the pre-CLEAN one minus --procs (plus the lane flags)."""
    import argparse

    seen, orig = [], argparse.ArgumentParser.parse_args
    monkeypatch.setattr(argparse.ArgumentParser, "parse_args",
                        lambda self, *a, **k: seen.append(self) or orig(self, *a, **k))
    with pytest.raises(SystemExit) as e:
        stub.main(["--help"])
    out = capsys.readouterr().out
    assert e.value.code == 0 and "procs" not in out and "scoring workers" not in out
    want = [x for o in OPTIONS if o != "--procs" for x in [o, *LANE_OPTIONS.get(o, [])]]
    assert [s for act in seen[0]._actions for s in act.option_strings] == want


def test_no_pool_and_one_freeze_after_the_last_preload(stub):
    """CT-34, CT-35, CT-36: no process pool; gc.freeze runs once, after the last preload and before training; rl.json
    keeps its pre-CLEAN keys."""
    import gc

    frozen = gc.get_freeze_count()
    assert stub.run(*stub.train) == 0
    assert stub.calls == ["preload", "freeze", "train"]
    with_heldout = stub.calls.clear() or stub.run(*stub.train, "--heldout", str(stub.out.parent / "root"))
    assert with_heldout == 0 and stub.calls == ["preload", "preload", "freeze", "train"]
    assert sorted(json.loads((stub.out / "rl.json").read_text())) == sorted(RL_KEYS + LANE_RL_KEYS)
    stub.calls.clear()
    assert stub.run(*stub.train) == 0 and stub.calls.count("freeze") == 1  # each main call freezes once
    assert gc.get_freeze_count() == frozen


def test_no_freeze_on_a_failed_start(stub, monkeypatch):
    """CT-35 failure paths: a parse error or a failed preload freezes nothing."""
    from objdel.train import serve

    with pytest.raises(SystemExit):
        stub.run("--train", "t")  # no --steps
    def missing(root, ds):
        raise FileNotFoundError(root)
    monkeypatch.setattr(serve, "preload", missing)
    with pytest.raises(FileNotFoundError):
        stub.run(*stub.train)
    assert "freeze" not in stub.calls


def test_train_rl_source_names_no_pool():
    """CT-34 static part and the CLEAN-G17 source gate."""
    import ast

    src = (ROOT / "experiments/phase2/train_rl.py").read_text()
    names = {n.id for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Name)}
    names |= {n.attr for n in ast.walk(ast.parse(src)) if isinstance(n, ast.Attribute)}
    assert not names & {"fork_pool", "pool_size", "Pool", "pool"}
    assert '"--procs"' not in src and src.count("gc.freeze()") == 1


def test_no_stale_procs_for_train_rl():
    """CT-84: no RL plan passes --procs to train_rl; the runbook names the three pool users."""
    for f in sorted((ROOT / "experiments/phase2/configs").glob("*.yaml")):
        for line in f.read_text().splitlines():
            assert not ("train_rl" in line and "--procs" in line), (f.name, line)
    text = (ROOT / "docs/runbooks/troubleshooting.md").read_text()
    row = next(x for x in text.splitlines() if "BrokenProcessPool" in x)
    assert all(m in row for m in ("select_ckpt", "g_probe", "final_eval")) and "train_rl" not in row
