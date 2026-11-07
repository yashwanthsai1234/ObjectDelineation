"""fast.ImageMemo: the processor's image step once per distinct image per generation batch.

TRL runs Qwen3-VL's image processor over all 256 rows twice a step (prompt ids before vLLM, pixel values for the
update), though a step has 32 distinct images; dataset.OpenOnce hands both passes the same image object for the 8 rows
of an image, so a memo keyed by the object serves both passes from one preparation per image."""

import numpy as np
import PIL.Image
import pytest
import torch
from transformers import Qwen2VLImageProcessor

from objdel.train import fast


def image(seed, h=64, w=64):
    return PIL.Image.fromarray(np.random.default_rng(seed).integers(0, 256, (h, w, 3), dtype=np.uint8))


@pytest.fixture
def counted():
    ip = Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)  # Qwen3-VL's patch and merge
    seen, real = [], ip.__class__.__call__

    class Counting(ip.__class__):
        def __call__(self, images, *a, **k):
            seen.append(len(images))
            return real(self, images, *a, **k)

    ip.__class__ = Counting
    return ip, seen


def test_memo_gives_the_processors_own_output_for_rows_that_repeat_images(counted):
    ip, _ = counted
    a, b, c = image(0), image(1), image(2, 96, 64)  # c: another size, so another grid
    rows = [a, a, b, c, a, b, c, a]
    want = Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)(rows, return_tensors="pt")
    for got in (fast.ImageMemo(ip)(rows, return_tensors="pt"), fast.ImageMemo(ip)([[x] for x in rows])):
        assert set(got.keys()) == set(want.keys())
        for k in want:
            assert torch.equal(got[k], want[k]), k


def test_each_distinct_image_is_prepared_once_until_the_memo_is_cleared(counted):
    ip, seen = counted
    memo = fast.ImageMemo(ip)
    a, b = image(0), image(1)
    memo([a] * 8 + [b] * 8)  # pass 1 (prompt ids)
    memo([[a]] * 8 + [[b]] * 8, return_tensors="pt")  # pass 2 (pixel values): other kwargs
    memo([a] * 8 + [b] * 8)  # pass 1's kwargs again: nothing new
    assert seen == [2, 2]  # 2 distinct images per kwargs, not 16
    memo.clear()
    memo([a] * 8)
    assert seen == [2, 2, 1]


def test_memo_passes_every_other_attribute_through(counted):
    ip, _ = counted
    memo = fast.ImageMemo(ip)
    assert (memo.merge_size, memo.patch_size) == (2, 16)


def test_fast_grpo_installs_one_memo_and_clears_it_before_each_generation(monkeypatch):
    proc = type("P", (), {"image_processor": Qwen2VLImageProcessor()})()
    memo = fast.ImageMemo.install(proc)
    assert proc.image_processor is memo and fast.ImageMemo.install(proc) is memo  # idempotent
    assert fast.ImageMemo.install(type("T", (), {})()) is None  # a text-only processing class: nothing to memoise
    memo.memo[("x",)] = "stale"
    t = object.__new__(fast.FastGRPO)
    t.images, t.use_vllm = memo, False
    monkeypatch.setattr(fast.GRPOTrainer, "_generate", lambda self, prompts: (None, [], None, []))
    monkeypatch.setattr(fast.FastGRPO, "timed", lambda self, phase: __import__("contextlib").nullcontext())
    monkeypatch.setattr(fast.FastGRPO, "prefetch_rewards", lambda self, *a: None)
    fast.FastGRPO._generate(t, ["p"])
    assert memo.memo == {}
