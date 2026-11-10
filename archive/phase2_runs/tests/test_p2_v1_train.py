"""Archived by CLEAN8 from tests/test_p2_v1_train.py (retired scripts and plans; nothing runs these)."""


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
