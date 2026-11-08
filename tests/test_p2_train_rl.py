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
