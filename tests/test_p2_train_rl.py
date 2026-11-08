"""experiments.phase2.train_rl: the full GRPO run. It trains for --steps, saves a resumable checkpoint every
--save-steps, resumes from the last one (a Vast box can die mid-run), and leaves the tuned model as a plain merged
model folder (the layout train_sft's branches use, so select_ckpt and vLLM read it unchanged)."""

import json
import os

import pytest
from test_p2_scripts import GRPO_SMOKE, _module, _tiny_vlm


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
                    "--procs", "1", "--grpo", json.dumps(GRPO_SMOKE), *extra, env=env)
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
                str(tmp_path / "prepared"), "--out", str(out), "--steps", "2", "--save-steps", "100", "--procs", "1",
                "--grpo", json.dumps(GRPO_SMOKE), "--shared-prefix", env=env)  # fmt: skip
    assert r.returncode == 0, r.stderr[-4000:]
    doc = json.loads((out / "rl.json").read_text())
    steps = [h for h in doc["log"] if "reward" in h]
    assert [h["step"] for h in steps] == [1, 2]
    assert all(0 < h["shared/token_ratio"] < 1 for h in steps)  # the packed path ran, on fewer tokens
    assert (out / "final/model.safetensors").exists()


@pytest.mark.slow
def test_full_rl_trains_on_shifted_windows_beside_the_centred_frames(tmp_path, synthetic_items):
    """--train-views: RL rows from the policy frames and from shifted / rescaled windows of the same train items
    (build_eval_views --split train), each scored on its own window's rater masks (rl_plateau_diagnosis_2026-11-08)."""
    pytest.importorskip("trl")
    from experiments.phase2 import build_eval_views as BE
    from objdel.data import store
    from objdel.train import dataset

    model = _tiny_vlm(tmp_path / "tiny", "qwen3-vl-2b")
    store.write(synthetic_items(("train", "val"), unanimous=True), tmp_path / "store")
    dataset.prepare(tmp_path / "store", tmp_path / "prepared", procs=1)
    BE.build(tmp_path / "store", tmp_path / "views", ("lidc",), "train", seed=0, procs=1)
    views = [str(tmp_path / "views" / f"lidc-train-{v}") for v in ("shifted", "rescaled")]
    out, env = tmp_path / "rl", os.environ | {"HF_HUB_OFFLINE": "1", "CUDA_VISIBLE_DEVICES": ""}
    r = _module("experiments.phase2.train_rl", "--model", "qwen3-vl-2b", "--ckpt", str(model), "--train",
                str(tmp_path / "prepared"), "--train-views", *views, "--out", str(out), "--steps", "2", "--save-steps",
                "100", "--procs", "1", "--grpo", json.dumps(GRPO_SMOKE), env=env)  # fmt: skip
    assert r.returncode == 0, r.stderr[-4000:]
    doc = json.loads((out / "rl.json").read_text())
    assert [h["step"] for h in doc["log"] if "reward" in h] == [1, 2]
    assert doc["train_rows"] == {"prepared": doc["train_rows"]["prepared"], "views": doc["train_rows"]["views"]}
    assert doc["train_rows"]["views"] > 0 and doc["train_rows"]["prepared"] > 0
