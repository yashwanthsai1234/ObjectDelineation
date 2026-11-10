"""gpu_pixels.patchify (plan C1): the processor's pixel values and grid, bit for bit, from uint8 images, on any device.

The Mac runs the CPU device only; the CUDA run is the box check (c_speed_check --pixels)."""

import numpy as np
import PIL.Image
import pytest

torch = pytest.importorskip("torch")
from transformers import Qwen2VLImageProcessor  # noqa: E402

from objdel.policy import gpu_pixels  # noqa: E402


@pytest.fixture(scope="module")
def ip():
    return Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)  # Qwen3-VL's patch and merge


@pytest.fixture(autouse=True)
def picked_torch(monkeypatch):
    """backend="auto" is what pick_backend picked; these tests call patchify without one."""
    monkeypatch.setattr(gpu_pixels, "picked", "torch")


def rand(seed, h=64, w=64):
    return np.random.default_rng(seed).integers(0, 256, (h, w, 3), dtype=np.uint8)


def digits(k, h=64, w=64):
    """Pixel = byte k of the element's flat index: a wrong index map shows as a wrong value (K3 DESIGN 6)."""
    return ((np.arange(h * w * 3, dtype=np.int64) >> (8 * k)) & 255).astype(np.uint8).reshape(h, w, 3)


def same(ip, arrays, backend="torch", **kw):
    want = ip([PIL.Image.fromarray(a).convert("RGB") for a in arrays], return_tensors="pt")
    pv, grid = gpu_pixels.patchify(arrays, ip, device="cpu", backend=backend, **kw)
    return torch.equal(pv, want["pixel_values"]) and torch.equal(grid, want["image_grid_thw"])


@pytest.mark.parametrize("backend", gpu_pixels.BACKENDS)
def test_patchify_equals_the_processor_on_random_flat_ramp_and_digit_images(ip, backend):
    ramp = np.broadcast_to(np.arange(256, dtype=np.uint8).reshape(16, 16)[..., None], (16, 16, 3))
    ramp = np.tile(ramp, (6, 4, 1))  # 96 x 64: every value on every channel, a second size
    cases = [rand(i) for i in range(32)] + [np.zeros((64, 64, 3), np.uint8), np.full((64, 64, 3), 255, np.uint8),
                                            np.ascontiguousarray(ramp)] + [digits(k, 96, 64) for k in range(3)]
    assert same(ip, cases, backend)


def test_patchify_on_a_512_image_and_a_uint8_tensor_batch(ip):
    a = rand(7, 512, 512)
    assert same(ip, [a])
    pv, grid = gpu_pixels.patchify(torch.from_numpy(np.stack([a, a])), ip, device="cpu")
    want = ip([PIL.Image.fromarray(a)] * 2, return_tensors="pt")
    assert torch.equal(pv, want["pixel_values"]) and grid.tolist() == [[1, 32, 32]] * 2


def test_a_grey_image_equals_the_processor_on_its_rgb_conversion(ip):
    g = np.random.default_rng(3).integers(0, 256, (64, 96), dtype=np.uint8)
    want = ip([PIL.Image.fromarray(g, mode="L").convert("RGB")], return_tensors="pt")
    for x in (PIL.Image.fromarray(g, mode="L"), g):
        pv, grid = gpu_pixels.patchify([x], ip, device="cpu")
        assert torch.equal(pv, want["pixel_values"]) and torch.equal(grid, want["image_grid_thw"])


def test_sizes_the_processor_would_resize_and_other_modes_are_refused(ip):
    with pytest.raises(gpu_pixels.UnsupportedImage):
        gpu_pixels.patchify([rand(0, 500, 500)], ip)
    with pytest.raises(gpu_pixels.UnsupportedImage):
        gpu_pixels.patchify([PIL.Image.new("RGBA", (64, 64))], ip)
    with pytest.raises(gpu_pixels.UnsupportedImage):
        gpu_pixels.patchify([rand(0)], Qwen2VLImageProcessor(patch_size=14, temporal_patch_size=2, merge_size=2))


def test_the_grid_stays_on_the_cpu_for_another_device(ip):
    _, grid = gpu_pixels.patchify([rand(0)], ip, device="meta")
    assert grid.device.type == "cpu"


def test_the_divisor_is_a_device_tensor_never_a_scalar(ip):  # U-PX16
    mean, std = gpu_pixels.fused(ip, "meta")
    assert std.device.type == "meta" and std.dim() > 0 and mean.dim() > 0


def test_a_scalar_reciprocal_differs_from_the_true_divide_it_must_not_replace():  # U-PX17
    x = torch.arange(256, dtype=torch.float32) - 127.5
    assert not torch.equal(x / 127.5, x * (1 / 127.5))


def test_pick_backend_returns_a_backend_that_matches_on_cpu_and_none_for_an_unsupported_image(ip):
    assert gpu_pixels.pick_backend([rand(0), digits(0)], ip, "cpu") in gpu_pixels.BACKENDS
    assert gpu_pixels.pick_backend([rand(0, 500, 500)], ip, "cpu") is None


def test_auto_is_the_default_and_means_the_backend_that_passed_the_self_check(ip, monkeypatch):
    """C1 interface: backend="auto" = the first backend that passed pick_backend; none picked yet is an error."""
    monkeypatch.setattr(gpu_pixels, "picked", None)
    with pytest.raises(RuntimeError, match="pick_backend"):
        gpu_pixels.patchify([rand(0)], ip)
    with pytest.raises(ValueError, match="bogus"):
        gpu_pixels.patchify([rand(0)], ip, backend="bogus")
    chosen = gpu_pixels.pick_backend([rand(0), digits(0)], ip, "cpu")
    assert gpu_pixels.picked == chosen
    assert same(ip, [rand(1), digits(1)], backend="auto") and same(ip, [rand(2)], backend=chosen)
    calls = []
    real = gpu_pixels._patchify
    monkeypatch.setattr(gpu_pixels, "_patchify", lambda *a: calls.append(1) or real(*a))
    monkeypatch.setattr(gpu_pixels, "picked", "torch")
    gpu_pixels.patchify([rand(3)], ip, device="cpu")  # the default backend is "auto"
    assert calls == [1]


def test_a_failed_self_check_clears_the_picked_backend(ip, monkeypatch):
    monkeypatch.setattr(gpu_pixels, "picked", "torch")
    assert gpu_pixels.pick_backend([rand(0, 500, 500)], ip, "cpu") is None
    assert gpu_pixels.picked is None


def test_hostile_values_equal_the_processor_single_and_batched(ip):  # U-PX2
    check = np.indices((64, 64)).sum(0) % 2 * 255
    levels = np.arange(64 * 64 * 3) % 256
    corner_lo, corner_hi = np.zeros((64, 64, 3), np.uint8), np.zeros((64, 64, 3), np.uint8)
    corner_lo[0, 0], corner_hi[-1, -1] = 255, 255
    cases = [np.zeros((64, 64, 3), np.uint8), np.full((64, 64, 3), 255, np.uint8),
             np.repeat(check[..., None], 3, 2).astype(np.uint8),
             np.random.default_rng(0).choice([127, 128], (64, 64, 3)).astype(np.uint8),
             levels.astype(np.uint8).reshape(64, 64, 3), corner_lo, corner_hi]  # fmt: skip
    per_level = [np.stack([np.full((64, 64), v), np.full((64, 64), (v + 85) % 256), np.full((64, 64), (v + 170) % 256)],
                          -1).astype(np.uint8) for v in range(256)]  # fmt: skip
    for a in cases:
        assert same(ip, [a])
    assert same(ip, cases) and same(ip, per_level)


def test_constants_and_settings_come_from_the_processor(ip):  # U-PX6
    custom = Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2, image_mean=[0.1, 0.2, 0.3],
                                   image_std=[0.7, 0.5, 0.3], rescale_factor=1 / 250)  # fmt: skip
    a = [rand(0), rand(1, 96, 64)]
    assert same(custom, a) and not torch.equal(gpu_pixels.patchify(a, custom)[0], gpu_pixels.patchify(a, ip)[0])
    for kw in ({"merge_size": 1}, {"merge_size": 4}, {"temporal_patch_size": 1}, {"do_rescale": False},
               {"do_normalize": False}):  # fmt: skip
        other = Qwen2VLImageProcessor(**{"patch_size": 16, "temporal_patch_size": 2, "merge_size": 2} | kw)
        with pytest.raises(gpu_pixels.UnsupportedImage):
            gpu_pixels.patchify(a, other)
    for kw in ({"do_convert_rgb": False}, {"resample": PIL.Image.Resampling.NEAREST}):  # equal, never wrong
        assert same(Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2, **kw), a)
    import inspect

    src = inspect.getsource(gpu_pixels)
    assert not any(c in src for c in ("0.48145466", "0.26862954", "0.5, 0.5, 0.5"))


@pytest.mark.parametrize("dtype", [torch.bfloat16, torch.float16])
def test_out_dtype_is_one_cast_of_the_float32_result_and_the_grid_stays_int64_on_the_cpu(ip, dtype):  # U-PX7
    a = [rand(0), rand(1, 96, 64)]
    ref, _ = gpu_pixels.patchify(a, ip)
    pv, grid = gpu_pixels.patchify(a, ip, out_dtype=dtype)
    assert pv.dtype == dtype and torch.equal(pv, ref.to(dtype))
    assert grid.dtype == torch.int64 and grid.device.type == "cpu"
    _, grid = gpu_pixels.patchify(a, ip, device="meta", out_dtype=dtype)
    assert grid.dtype == torch.int64 and grid.device.type == "cpu"


def test_no_image_raises_value_error_and_an_empty_prefill_never_calls_patchify(ip, monkeypatch):  # U-PX9
    for empty in ([], torch.zeros((0, 512, 512, 3), dtype=torch.uint8)):
        with pytest.raises(ValueError, match="no image"):
            gpu_pixels.patchify(empty, ip)
    from objdel.policy import caches

    monkeypatch.setattr(gpu_pixels, "patchify", lambda *a, **k: pytest.fail("patchify called"))
    memo = caches.ImageMemo(ip, gpu=("torch", "cpu"))
    memo.prefill([], return_tensors="pt")
    assert (memo.misses, memo.memo) == (0, {})


def test_inputs_are_untouched_never_aliased_and_layout_free(ip):  # U-PX10
    import warnings

    a = rand(0, 64, 96)
    before = a.copy()
    pv, _ = gpu_pixels.patchify([a], ip)
    assert np.array_equal(a, before) and not np.shares_memory(pv.numpy(), a)
    pv[:] = 0  # the result is its own memory
    assert np.array_equal(a, before)
    want = gpu_pixels.patchify([a], ip)[0]
    wide = np.ascontiguousarray(np.repeat(a, 2, 1))[:, ::2]  # a non-contiguous view of the same values
    ro = a.copy()
    ro.setflags(write=False)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # no "non-writable" warning
        for x in (wide, np.asfortranarray(a), ro):
            assert not x.flags.c_contiguous or not x.flags.writeable
            assert torch.equal(gpu_pixels.patchify([x], ip)[0], want)
    assert np.array_equal(ro, before)


def test_duplicate_images_keep_their_order(ip):  # U-PX11
    a, b, c = rand(0), rand(1), rand(2, 96, 64)
    assert same(ip, [a, b, a, c])


def test_both_temporal_slices_of_every_patch_row_are_identical(ip):  # U-PX12
    pv, _ = gpu_pixels.patchify([rand(0), rand(1, 96, 64)], ip)
    rows = pv.view(pv.size(0), 3, 2, 16 * 16)  # C, t, p * p
    assert torch.equal(rows[:, :, 0], rows[:, :, 1])


@pytest.mark.data
def test_real_lidc_views_and_a_riga_grey_view_equal_the_processor_or_raise(ip):  # U-PX13
    from pathlib import Path

    root = Path(__file__).resolve().parents[2] / "data"
    kinds = [root / "prepared/images/lidc", root / "eval_views/lidc-val-shifted/images",
             root / "eval_views/lidc-val-rescaled/images"]  # fmt: skip
    files = [f for k in kinds for f in sorted(k.rglob("*.png"))[:8]]
    riga = sorted((root / "prepared/images/riga").rglob("*.png"))[:1]
    if len(files) < 24 or not riga:
        pytest.skip("no LIDC centred/shifted/rescaled views or RIGA images under data/")
    assert same(ip, [np.asarray(PIL.Image.open(f).convert("RGB")) for f in files])
    grey = PIL.Image.open(riga[0]).convert("L")
    with pytest.raises(gpu_pixels.UnsupportedImage):  # 500 x 500: the processor resizes it, so patchify refuses
        gpu_pixels.patchify([grey], ip)
    crop = grey.crop((0, 0, 480, 480))  # a RIGA grey view the processor keeps as it is
    want = ip([crop.convert("RGB")], return_tensors="pt")
    pv, grid = gpu_pixels.patchify([crop], ip)
    assert torch.equal(pv, want["pixel_values"]) and torch.equal(grid, want["image_grid_thw"])


def test_an_untrusted_png_header_with_a_huge_size_raises_before_its_pixels_are_read(ip, tmp_path):  # U-PX14
    import struct
    import warnings
    import zlib

    path = tmp_path / "huge.png"
    PIL.Image.new("RGB", (8, 8)).save(path)
    data = bytearray(path.read_bytes())
    ihdr = data[12:29]  # b"IHDR" + 13 bytes
    ihdr[4:12] = struct.pack(">II", 12000, 12000)  # 144 M pixels, 432 MB of RGB if read
    data[12:29] = ihdr
    data[29:33] = struct.pack(">I", zlib.crc32(bytes(ihdr)))
    path.write_bytes(bytes(data))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PIL.Image.DecompressionBombWarning)
        img = PIL.Image.open(path)
    assert img.size == (12000, 12000)
    img.load = lambda: pytest.fail("pixels read")
    with pytest.raises(gpu_pixels.UnsupportedImage, match="resize"):
        gpu_pixels.patchify([img], ip)


def test_staging_is_a_fresh_pinned_buffer_per_call_on_cuda_and_plain_on_the_cpu(ip, monkeypatch):  # U-PX15
    calls, real = [], torch.empty

    class Staged(Exception):
        pass

    def spy(*a, pin_memory=False, **k):
        calls.append(pin_memory)
        if pin_memory:
            raise Staged  # the Mac has no CUDA to copy to: the staging allocation is what is checked
        return real(*a, **k)

    monkeypatch.setattr(torch, "empty", spy)
    gpu_pixels.patchify([rand(0)], ip, device="cpu")
    assert calls == [False]
    for _ in range(2):
        with pytest.raises(Staged):
            gpu_pixels.patchify([rand(0)], ip, device="cuda", backend="hf")
    assert calls == [False, True, True]
    assert not [k for k, v in vars(gpu_pixels).items() if isinstance(v, torch.Tensor)]  # no module-level buffer
