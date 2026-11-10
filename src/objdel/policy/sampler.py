"""Sampling helpers shared by RL and SFT. Today: `stop_ids` of a checkpoint directory (E5 moves the JSONL bridge of
objdel.train.sample here)."""

from pathlib import Path


def stop_ids(path: Path, proc) -> set[int]:
    """rewards.stop_ids of a checkpoint directory. Base weights may have no generation_config.json (eos only in
    text_config): then the model config's, as from_pretrained falls back for a loaded model."""
    from types import SimpleNamespace

    from transformers import AutoConfig, GenerationConfig

    from objdel.train import rewards

    try:
        gen = GenerationConfig.from_pretrained(path)
    except OSError:
        gen = GenerationConfig.from_model_config(AutoConfig.from_pretrained(path))  # reads text_config's eos too
    return rewards.stop_ids(SimpleNamespace(generation_config=gen), proc)
