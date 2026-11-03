"""V1 readiness of train_sft: deterministic mode without full_determinism (PLAN D24), the run manifest's timing
block, C2's depth read from R2's g_probe.json (branches_from) and runner-only config keys (project_r1, the R0 stop
rule: test_p2_hf_eval)."""

import json
import os
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "experiments/phase2/configs"
C2 = sorted(CONFIGS.glob("c2_*.yaml")) + [CONFIGS / "rx8.yaml"]


def _yaml(path):
    return pytest.importorskip("yaml").safe_load(path.read_text())


@pytest.fixture
def TS():
    pytest.importorskip("trl")
    from experiments.phase2 import train_sft

    return train_sft


def _view(tmp_path, epochs):
    (tmp_path / "v").mkdir()
    (tmp_path / "v/manifest.json").write_text(json.dumps({"config": {"epochs": epochs}}))
    return str(tmp_path / "v")


# ---------------------------------------------------------------- 1. deterministic mode (D24)


def test_no_config_uses_full_determinism():
    for p in CONFIGS.glob("*.yaml"):
        assert "full_determinism" not in _yaml(p).get("sft", {}), p.name  # it sets CUDA_LAUNCH_BLOCKING=1


def test_deterministic_is_applied_before_anything_loads(TS, tmp_path, monkeypatch):
    import torch

    monkeypatch.delenv("CUBLAS_WORKSPACE_CONFIG", raising=False)
    a = TS.RunArgs("lfm2.5-vl-450m", _view(tmp_path, 0), str(tmp_path / "o"), vision=False, deterministic=True)
    try:
        with pytest.raises(ValueError, match="view epoch"):
            TS.train(a)
        assert os.environ["CUBLAS_WORKSPACE_CONFIG"] == ":4096:8"
        assert torch.are_deterministic_algorithms_enabled() and torch.is_deterministic_algorithms_warn_only_enabled()
        monkeypatch.setenv("CUBLAS_WORKSPACE_CONFIG", ":16:8")  # an explicit setting wins
        with pytest.raises(ValueError, match="view epoch"):
            TS.train(TS.RunArgs("m", a.view, str(tmp_path / "o"), vision=False, deterministic=True, run_dir=""))
        assert os.environ["CUBLAS_WORKSPACE_CONFIG"] == ":16:8"
    finally:
        torch.use_deterministic_algorithms(False)
    monkeypatch.delenv("CUBLAS_WORKSPACE_CONFIG")
    with pytest.raises(ValueError, match="view epoch"):
        TS.train(TS.RunArgs("m", a.view, str(tmp_path / "o"), vision=False, run_dir=""))
    assert "CUBLAS_WORKSPACE_CONFIG" not in os.environ and not torch.are_deterministic_algorithms_enabled()


# ---------------------------------------------------------------- 2. timing block


def test_timing_excludes_the_warmup_steps(TS, monkeypatch):
    import torch

    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    t = TS.timing([10.0, 12.0, 14.0, 15.0, 16.0, 17.0], 32)  # step-end clocks: steps 4-6 take 1 s each
    assert t == {"steps": 6, "warmup_excluded": 3, "s_per_step": 1.0, "samples_per_s": 32.0, "batch": 32,
                 "max_memory_gib": None}  # fmt: skip
    assert TS.timing([1.0, 2.0, 3.0], 32)["s_per_step"] is None  # nothing after the warm-up
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch.cuda, "max_memory_allocated", lambda: 3 * 2**30)
    assert TS.timing([], 32) == {
        "steps": 0, "warmup_excluded": 3, "s_per_step": None, "samples_per_s": None, "batch": 32, "max_memory_gib": 3.0
    }  # fmt: skip


def test_peak_memory_is_none_on_cpu(TS, monkeypatch):
    import torch

    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    assert TS.peak_gib() is None
    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setattr(torch.cuda, "max_memory_allocated", lambda: 2**29)
    assert TS.peak_gib() == 0.5


def test_alloc_conf_is_set_before_anything_loads(TS, tmp_path, monkeypatch):
    """PLAN D4: expandable segments (fragmentation), read by the CUDA caching allocator at its first allocation."""
    monkeypatch.setenv("PYTORCH_ALLOC_CONF", "")  # restored after the test
    monkeypatch.delenv("PYTORCH_ALLOC_CONF")
    a = TS.RunArgs("lfm2.5-vl-450m", _view(tmp_path, 0), str(tmp_path / "o"), vision=False)
    with pytest.raises(ValueError, match="view epoch"):
        TS.train(a)
    assert os.environ["PYTORCH_ALLOC_CONF"] == "expandable_segments:True"
    monkeypatch.setenv("PYTORCH_ALLOC_CONF", "max_split_size_mb:64")  # an explicit setting wins
    with pytest.raises(ValueError, match="view epoch"):
        TS.train(TS.RunArgs("m", a.view, str(tmp_path / "o"), vision=False, run_dir=""))
    assert os.environ["PYTORCH_ALLOC_CONF"] == "max_split_size_mb:64"


def _tiny_args(TS, tmp_path, synthetic_items, model_path=None, **kw):
    """RunArgs of a CPU run on a 2-epoch synthetic view (branches 0.5 / 1 at 10 steps/epoch unless overridden)."""
    from test_p2_scripts import CPU_SFT, _tiny_vlm

    from experiments.phase2 import build_views as BV
    from objdel.data import store

    store.write(synthetic_items(("train",), unanimous=True), tmp_path / "store")
    BV.build(tmp_path / "store", tmp_path / "view", "lidc", "train", epochs=2, seed=3, procs=1)
    o = {"branches": [0.5, 1.0], "steps_per_epoch": 10, "vision": False, "rolling_minutes": 1e3, "sft": CPU_SFT,
         "model_path": model_path or str(_tiny_vlm(tmp_path / "m"))}
    return TS.RunArgs("lfm2.5-vl-450m", str(tmp_path / "view"), str(tmp_path / "o"), **(o | kw))


def _manifest_of(tmp_path):
    return json.loads((tmp_path / "o/manifest.json").read_text())


@pytest.mark.slow
def test_manifest_records_timing(TS, tmp_path, synthetic_items):
    assert TS.train(_tiny_args(TS, tmp_path, synthetic_items, branches=[], steps_per_epoch=6)) == 0
    t = _manifest_of(tmp_path)["timing"]
    assert t["steps"] == 6 and t["s_per_step"] > 0 and t["samples_per_s"] == pytest.approx(1 / t["s_per_step"])
    assert t["batch"] == 1 and t["rows_per_epoch"] > 6 and t["max_memory_gib"] is None


@pytest.mark.slow
def test_resume_keeps_the_first_timing(TS, tmp_path, synthetic_items):
    """A resumed stable segment is partial (and warm-up again): the timing measured first stands."""
    a = _tiny_args(TS, tmp_path, synthetic_items, branches=[], steps_per_epoch=6)
    assert TS.train(a) == 0
    man = _manifest_of(tmp_path) | {"status": "running", "timing": {"s_per_step": 1.25}}
    (tmp_path / "o/manifest.json").write_text(json.dumps(man))
    assert TS.train(_tiny_args(TS, tmp_path, synthetic_items, a.model_path, branches=[], steps_per_epoch=6,
                               resume=True)) == 0  # fmt: skip
    assert _manifest_of(tmp_path)["timing"] == {"s_per_step": 1.25}


@pytest.mark.slow
def test_each_trainer_is_freed_before_the_next_one_loads(TS, tmp_path, synthetic_items, monkeypatch):
    """V1 R1 (2026-11-02, RTX 4090): the branch loop kept the stable and e0.5 trainers alive (model, grads, fp32 AdamW
    state, ~7 GB each), so branch e1 hit CUDA OOM at its first step. R0's 30-step timing arms build one trainer.
    Also: rolling/ (full checkpoints with AdamW state) is gone once every stable/e<b> exists, and each phase's peak
    CUDA memory is in the manifest (R0's 30-step arm peaked at 16.4 GiB, R1's stable phase at 20.9 GiB)."""
    import gc
    import weakref

    import trl

    models, rolling = [], []

    class Spy(trl.SFTTrainer):
        def __init__(self, *a, **k):
            gc.collect()
            assert all(r() is None for r in models), "an earlier trainer's model is still referenced"
            rolling.append((tmp_path / "o/rolling").exists())
            super().__init__(*a, **k)
            models.append(weakref.ref(self.model))

    monkeypatch.setattr(trl, "SFTTrainer", Spy)
    assert TS.train(_tiny_args(TS, tmp_path, synthetic_items)) == 0 and len(models) == 3  # stable + two branches
    assert rolling == [False, False, False]
    assert _manifest_of(tmp_path)["peak_memory_gib"] == {"stable": None, "e0.5": None, "e1": None}


@pytest.mark.slow
def test_resume_refuses_a_lost_stable_copy(TS, tmp_path, synthetic_items):
    """A crash between the rolling save at a branch point and its stable/e<b> copy: resuming would train past the
    point and never write stable/e<b>, so the branch would fail hours later. Refused up front."""
    a = _tiny_args(TS, tmp_path, synthetic_items, resume=True)
    (tmp_path / "o/rolling/checkpoint-5").mkdir(parents=True)
    (tmp_path / "o/rolling/checkpoint-5/trainer_state.json").write_text("{}")
    view = json.loads((tmp_path / "view/manifest.json").read_text())
    old = {"config": {k: getattr(a, k) for k in TS.SAME_RUN}, "status": "running",
           "view": {"view_hash": view["view_hash"]}}
    (tmp_path / "o/manifest.json").write_text(json.dumps(old))
    with pytest.raises(ValueError, match=r"checkpoint-5 is past branch e0.5 .*stable/e0.5 is missing"):
        TS.train(a)


@pytest.mark.slow
def test_a_vision_abort_in_a_decay_saves_no_branch(TS, tmp_path, synthetic_items, monkeypatch):
    """G-vision aborts inside branch e0.5's decay: the branch is unfinished, so branches/e0.5 must not exist."""
    out = tmp_path / "o"
    monkeypatch.setattr(TS, "preload", lambda cfg, model: {"items": [], "set_sha256": "s", "reference": 1.0})
    monkeypatch.setattr(TS, "vision_check", lambda model, proc, vis: {
        "action": "abort" if (out / "decay").exists() else "ok", "ratio": 0.1, "tv": 0.001})
    a = _tiny_args(TS, tmp_path, synthetic_items, vision={"prepared": "p", "reference": 1.0}, rolling_minutes=0.0)
    assert TS.train(a) == TS.ABORTED
    man = _manifest_of(tmp_path)
    assert man["status"] == "aborted" and man["abort"]["step"] == 6  # e0.5 = step 5, decay 1 step
    assert not (out / "branches").exists() and (out / "stable/e0.5").is_dir()


class _Stop(Exception):
    pass


def _resume_branch(TS, tmp_path, synthetic_items, monkeypatch, base, ckpt):
    """ckpt (a saved model) as stable/e1, a finished stable phase; train() up to the moment the e1 branch trainer has
    loaded it: returns that trainer's model state_dict (the decay itself is not run)."""
    import trl
    from transformers.trainer_callback import TrainerState

    loaded = {}

    class Spy(trl.SFTTrainer):
        def _inner_training_loop(self, *a, **k):  # called after _load_from_checkpoint and the trainer state
            loaded.update(self.model.state_dict())
            raise _Stop

    monkeypatch.setattr(trl, "SFTTrainer", Spy)
    a = _tiny_args(TS, tmp_path, synthetic_items, model_path=str(base), branches=[1.0])
    shutil.copytree(ckpt, tmp_path / "o/stable/e1")
    TrainerState(global_step=10, max_steps=10).save_to_json(str(tmp_path / "o/stable/e1/trainer_state.json"))
    with pytest.raises(_Stop):
        TS.train(a)
    return loaded


@pytest.mark.slow
def test_branch_resume_restores_the_trained_vision_tower(TS, tmp_path, synthetic_items, monkeypatch):
    """V1 R1: LFM2.5-VL saves its tower as model.vision_tower.vision_model.*, the live module is model.vision_tower.*;
    Trainer's load_state_dict(strict=False) skipped all 198 tower tensors with a warning, so every branch (and every
    crash resume) restarted the tower from base. Real weights: the tiny LFM config does not rename them."""
    torch = pytest.importorskip("torch")
    from huggingface_hub import snapshot_download
    from transformers import AutoModelForImageTextToText

    monkeypatch.setenv("HF_HUB_OFFLINE", "1")
    try:
        src = snapshot_download("LiquidAI/LFM2.5-VL-450M", local_files_only=True)
    except Exception:
        pytest.skip("LFM2.5-VL-450M not in the local HF cache")
    model = AutoModelForImageTextToText.from_pretrained(src, dtype=torch.float32)
    with torch.no_grad():
        model.get_parameter("model.vision_tower.embeddings.patch_embedding.weight").add_(0.5)  # "trained"
    model.save_pretrained(tmp_path / "ckpt")
    want = {k: v.detach().clone() for k, v in model.state_dict().items() if ".vision_tower." in k}
    del model
    loaded = _resume_branch(TS, tmp_path, synthetic_items, monkeypatch, src, tmp_path / "ckpt")
    assert len(want) > 190 and all(torch.equal(loaded[k], v) for k, v in want.items())
    assert torch.equal(loaded["lm_head.weight"], loaded["model.language_model.embed_tokens.weight"])  # still tied


@pytest.mark.slow
def test_branch_resume_refuses_a_checkpoint_missing_weights(TS, tmp_path, synthetic_items, monkeypatch):
    import safetensors.torch
    from test_p2_scripts import _tiny_vlm

    base = _tiny_vlm(tmp_path / "m")
    shutil.copytree(base, tmp_path / "ckpt")
    sd = safetensors.torch.load_file(tmp_path / "ckpt/model.safetensors")
    sd.pop(next(k for k in sd if "vision" in k))
    safetensors.torch.save_file(sd, tmp_path / "ckpt/model.safetensors", metadata={"format": "pt"})
    with pytest.raises(ValueError, match="missing"):
        _resume_branch(TS, tmp_path, synthetic_items, monkeypatch, base, tmp_path / "ckpt")


@pytest.mark.parametrize("name", ["r1_qwen3vl2b", "r0_qwen3vl_arm", "c2_a0p_seed2"])
def test_qwen3_vl_configs_train_as_one_recipe(name):
    """The 30-step arm and the C2 runs train exactly as R1 does: LM-only LoRA r16 (a 24 GB card cannot hold full
    fine-tuning of 2B parameters), micro-batch 4 x 8 (152k-token vocabulary logits), checkpointing on (the model's
    default, so no `gradient_checkpointing` key may switch it off)."""
    c = _yaml(CONFIGS / f"{name}.yaml")
    assert c["model"] == "qwen3-vl-2b" and "lora" not in c and "gradient_checkpointing" not in c["sft"]
    assert (c["sft"]["per_device_train_batch_size"], c["sft"]["gradient_accumulation_steps"]) == (4, 8)


# ---------------------------------------------------------------- 4. C2 depth from R2 (branches_from)


def test_c2_configs_take_their_depth_from_g_probe():
    for p in C2:
        c = _yaml(p)
        assert c["branches_from"] == "r1-qwen3vl2b/g_probe/g_probe.json" and "branches" not in c, p.name
        assert c["stages"] == ["train", "topology"], p.name


def _probe(tmp_path, **doc):
    (tmp_path / "r1-lfm/g_probe").mkdir(parents=True)
    (tmp_path / "r1-lfm/g_probe/g_probe.json").write_text(json.dumps(doc))
    return "r1-lfm/g_probe/g_probe.json"


def test_branches_from_reads_the_final_verdict_checkpoint(TS, tmp_path):
    ckpt = "/workspace/runs/v1/r1-lfm/branches/e0.5"  # D12(c): the repair branch e0.5 got the final verdict
    rel = _probe(tmp_path, verdict="INCONCLUSIVE", ckpt=ckpt, probed=["e1", "e0.5"])
    a = TS.RunArgs("lfm2.5-vl-450m", _view(tmp_path, 0), "c2", run_dir=str(tmp_path), branches_from=rel, vision=False)
    with pytest.raises(ValueError, match="decays to 0.575 epochs"):  # resolved to [0.5] before the view check
        TS.train(a)
    assert a.branches == [0.5]


@pytest.mark.parametrize(
    "doc", [None, {"verdict": None, "skipped": "gates"}, {"ckpt": "b/e1"}, {"verdict": "PASS", "ckpt": "c"}]
)
def test_branches_from_refuses_without_a_verdict(TS, tmp_path, doc):
    rel = _probe(tmp_path, **doc) if doc else "r1-lfm/g_probe/g_probe.json"
    a = TS.RunArgs("lfm2.5-vl-450m", _view(tmp_path, 3), "c2", run_dir=str(tmp_path), branches_from=rel, vision=False)
    with pytest.raises(SystemExit, match="C2 must not start"):
        TS.train(a)
    assert not (tmp_path / "c2").exists()


@pytest.mark.parametrize("path", C2 + [CONFIGS / "r1_lfm.yaml"], ids=lambda p: p.stem)
def test_runner_only_keys_parse(TS, path, tmp_path):
    from trl import TrlParser

    extra = {"sampler": "hf", "stages": ["train"]}
    (tmp_path / "c.yaml").write_text(json.dumps(_yaml(path) | extra))  # JSON is YAML
    (a,) = TrlParser([TS.RunArgs]).parse_args_and_config(["--config", str(tmp_path / "c.yaml")])
    assert a.sampler == "hf" and a.stages == ["train"]

