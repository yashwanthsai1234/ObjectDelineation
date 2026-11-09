"""gpu_pixels.patchify (plan C1): the processor's pixel values and grid, bit for bit, from uint8 images, on any device.

The Mac runs the CPU device only; the CUDA run is the box check (c_speed_check --pixels)."""

import numpy as np
import PIL.Image
import pytest

torch = pytest.importorskip("torch")
from transformers import Qwen2VLImageProcessor  # noqa: E402

from objdel.train import gpu_pixels  # noqa: E402


@pytest.fixture(scope="module")
def ip():
    return Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)  # Qwen3-VL's patch and merge


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
