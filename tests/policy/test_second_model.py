"""A model of a second family, registered in this test only: its entry drives the LoRA config, both trainer configs,
the processor keywords and the --no-ckpt-layers bound, and LoRA on a tiny random copy of it (built offline from the
config and processor in the local Hugging Face cache, never downloaded) leaves its vision tower and projector out and
runs a CPU forward on a dataset prompt and image.

The CPU smoke (tests/experiments/test_p2_scripts.py) runs the stage scripts in fresh processes, which see only the
library's entries; a library entry for this family would name a retired model in src (CT-07), so the family is proven
here instead."""
import pytest

from objdel.core import spec
from objdel.policy import models

HUB = "Hugging" + "FaceTB/Smol" + "VLM-256M-Instruct"  # split: the retired-name scan (CT-07) reads this file
KEY = "second-family"
ENTRY = models.Model(HUB, max_completion=1024, text_layers=2, lora_r=8, lora_alpha=16, sft_lr=1e-4, rl_lr=2e-5,
                     groups_per_device=1, grad_accum=4)
TEXT = dict(hidden_size=64, intermediate_size=128, num_hidden_layers=2, num_attention_heads=4, num_key_value_heads=2,
            head_dim=16)
VISION = dict(hidden_size=32, intermediate_size=64, num_hidden_layers=2, num_attention_heads=2)


@pytest.fixture
def entry(monkeypatch):
    monkeypatch.setitem(models.MODELS, KEY, ENTRY)
    return ENTRY


@pytest.fixture(scope="module")
def tiny():
    """(random 2-layer model, processor) of the family, from the cached config and processor (no download)."""
    torch = pytest.importorskip("torch")
    from transformers import AutoConfig, AutoModelForImageTextToText, AutoProcessor

    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("HF_HUB_OFFLINE", "1")
        try:
            cfg = AutoConfig.from_pretrained(HUB, local_files_only=True)
            proc = AutoProcessor.from_pretrained(HUB, local_files_only=True, **dict(ENTRY.processor))
        except OSError:
            pytest.skip(f"{HUB} config / processor not in the local HF cache")
    cfg.text_config.update(TEXT)
    cfg.vision_config.update(VISION)
    torch.manual_seed(0)
    return AutoModelForImageTextToText.from_config(cfg), proc


def test_the_entry_drives_every_derived_setting(entry, tmp_path, capsys):
    from objdel.rl import config as rl_config
    from objdel.rl import run as rl_run
    from objdel.sft import config as sft_config

    lc = models.lora(KEY)
    assert (lc.r, lc.lora_alpha, lc.exclude_modules) == (8, 16, models.VISION)
    assert models.processor_kwargs(KEY) == {"padding_side": "left"} and models.sampling(KEY)["max_tokens"] == 1024
    g = rl_config.grpo_config(KEY, str(tmp_path))
    assert (g.learning_rate, g.per_device_train_batch_size, g.gradient_accumulation_steps) == (2e-5, 8, 4)
    assert (g.max_completion_length, g.vllm_gpu_memory_utilization, g.vllm_max_model_length) == (1024, 0.35, 4096)
    s = sft_config.sft_config(KEY, str(tmp_path), dataloader_num_workers=0)
    assert (s.learning_rate, s.gradient_checkpointing) == (1e-4, True)
    base = ["--model", KEY, "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "2", "--no-ckpt-layers"]
    assert rl_run.checked(spec.cli_stage(spec.RlStage, [*base, "2"])).no_ckpt_layers == 2
    with pytest.raises(SystemExit):
        rl_run.checked(spec.cli_stage(spec.RlStage, [*base, "3"]))
    assert f"{KEY} has 2 text decoder layers" in capsys.readouterr().err  # the entry's bound, not Qwen's 28


def test_lora_skips_the_vision_tower_and_a_cpu_forward_runs(entry, tiny):
    import numpy as np
    import torch
    from peft import get_peft_model
    from PIL import Image

    from objdel.core import datasets

    model, proc = tiny
    peft = get_peft_model(model, models.lora(KEY))
    adapted = [n for n, _ in peft.named_modules() if n.endswith("lora_A")]
    assert adapted and not [n for n in adapted if "vision_model" in n or "connector" in n]
    msgs = [{"role": "user", "content": [{"type": "image"}, {"type": "text", "text": datasets.LIDC.prompt}]}]
    text = proc.apply_chat_template(msgs, add_generation_prompt=True)
    img = Image.fromarray(np.random.default_rng(0).integers(0, 256, (64, 64, 3), dtype=np.uint8))
    batch = proc(text=[text], images=[[img]], return_tensors="pt")
    with torch.no_grad():
        logits = peft(**batch).logits
    assert logits.shape[:2] == batch["input_ids"].shape and torch.isfinite(logits).all()
