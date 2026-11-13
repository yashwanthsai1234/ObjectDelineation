"""objdel.policy.models.processor and load_pair (MT-01, MT-10): one place loads a model entry's processor and an
inference (model, processor) pair, with each caller's own keywords; no model is downloaded (from_pretrained is
patched at class level, as the callers' tests do)."""
from types import SimpleNamespace

import pytest

from objdel.policy import models

transformers = pytest.importorskip("transformers")


def test_processor_uses_the_entry_keywords_and_hub_id(monkeypatch):
    seen = []
    monkeypatch.setattr(transformers.AutoProcessor, "from_pretrained", lambda path, **kw: seen.append((path, kw)))
    models.processor("qwen3-vl-2b")
    models.processor("qwen3-vl-2b", "/ckpt")
    want = {"padding_side": "left"}
    assert seen == [("Qwen/Qwen3-VL-2B-Instruct", want), ("/ckpt", want)]
    with pytest.raises(KeyError, match="unknown model"):
        models.processor("nope", "/ckpt")


class _Model:
    def __init__(self, log):
        self.log = log

    def to(self, device):
        self.log.append(("to", device))
        return self

    def eval(self):
        self.log.append(("eval",))
        return self


def test_load_pair_passes_each_callers_dtype_device_and_keywords(monkeypatch):
    log = []
    monkeypatch.setattr(transformers.AutoModelForImageTextToText, "from_pretrained",
                        lambda path, **kw: log.append(("model", path, kw)) or _Model(log))
    monkeypatch.setattr(transformers.AutoProcessor, "from_pretrained",
                        lambda path, **kw: log.append(("proc", path, kw)) or SimpleNamespace())
    models.load_pair("/s", dtype="bf16", device="cuda", proc_kw={"padding_side": "left"})  # HF sampler
    assert log == [("model", "/s", {"dtype": "bf16"}), ("to", "cuda"), ("eval",), ("proc", "/s", {"padding_side": "left"})]
    log.clear()
    models.load_pair("/b", dtype="float32", device_map="cpu")  # G-vision in select: the saved processor as is
    assert log == [("model", "/b", {"dtype": "float32", "device_map": "cpu"}), ("eval",), ("proc", "/b", {})]
