"""fast.ImageMemo: the processor's image step once per distinct image per generation batch.

TRL runs Qwen3-VL's image processor over all 256 rows twice a step (prompt ids before vLLM, pixel values for the
update), though a step has 32 distinct images; dataset.OpenOnce hands both passes the same image object for the 8 rows
of an image, so a memo keyed by the object serves both passes from one preparation per image."""

import numpy as np
import PIL.Image
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from transformers import Qwen2VLImageProcessor  # noqa: E402

from objdel.train import fast  # noqa: E402


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


def test_a_saved_processor_holds_the_plain_image_processor_not_the_memo(tmp_path, monkeypatch):
    """TRL saves the processor with every checkpoint; ProcessorMixin.to_dict deep-copies its attributes. A copy of
    the memo is the plain image processor (no cached pixels), so the saved files equal those saved without it."""
    import copy
    import os

    from transformers import AutoProcessor

    monkeypatch.setenv("HF_HUB_OFFLINE", "1")
    try:
        proc = AutoProcessor.from_pretrained("Qwen/Qwen3-VL-2B-Instruct")
    except OSError:
        pytest.skip("Qwen/Qwen3-VL-2B-Instruct processor not in the local HF cache")
    proc.save_pretrained(tmp_path / "plain")
    memo = fast.ImageMemo.install(proc)
    memo([image(0)])  # something cached
    plain = copy.deepcopy(memo)
    assert type(plain) is type(memo.inner) and plain.to_dict() == memo.inner.to_dict()
    proc.save_pretrained(tmp_path / "memo")
    for f in sorted(os.listdir(tmp_path / "plain")):
        assert (tmp_path / "memo" / f).read_bytes() == (tmp_path / "plain" / f).read_bytes(), f
    assert proc.image_processor is memo  # saving left the memo installed


def test_gpu_memo_equals_the_processor_and_counts_one_miss_per_distinct_image(counted):  # C1
    ip, seen = counted
    imgs = [image(i) for i in range(32)]
    rows = [x for x in imgs for _ in range(8)]  # 256 rows, 8 per image
    want = Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)(rows, return_tensors="pt")
    memo = fast.ImageMemo(ip, gpu=("torch", "cpu"))
    got = memo(rows, return_tensors="pt")
    assert all(torch.equal(got[k], want[k]) for k in want) and seen == []  # never the CPU processor
    assert (memo.misses, memo.hits) == (32, 224)
    memo.clear()
    memo.prefill(imgs, return_tensors="pt")
    memo(rows, return_tensors="pt")
    assert (memo.misses, memo.hits) == (32, 256)


def test_per_image_memo_returns_each_distinct_images_patches_once_with_the_row_grid(counted):  # C1 --encoder-once
    ip, _ = counted
    imgs = [image(i) for i in range(4)]
    rows = [x for x in imgs for _ in range(8)]
    plain = Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)
    got = fast.ImageMemo(ip, gpu=("torch", "cpu"), per_image=True)(rows, return_tensors="pt")
    assert torch.equal(got["pixel_values"], plain(imgs, return_tensors="pt")["pixel_values"])
    assert torch.equal(got["image_grid_thw"], plain(rows, return_tensors="pt")["image_grid_thw"])


def test_an_unsupported_image_falls_back_to_the_processor_and_is_counted(counted):
    ip, seen = counted
    rows = [image(0, 50, 50)] * 2  # 50 x 50: the processor resizes it
    want = Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)(rows, return_tensors="pt")
    memo = fast.ImageMemo(ip, gpu=("torch", "cpu"))
    got = memo(rows, return_tensors="pt")
    assert all(torch.equal(got[k], want[k]) for k in want) and memo.fallback_calls == 1 and seen == [1]
