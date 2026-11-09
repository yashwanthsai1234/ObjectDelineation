"""Qwen-VL pixel values built on the GPU from uint8 images (plan C1): one H2D copy of the uint8 pixels, then the
processor's fused rescale-normalize and patchify on the device, bit for bit the processor's own output.

    uint8 (N, H, W, 3|1) --pinned H2D--> float32 - mean' / std' --> patchify (HF's reshape/permute) --> (rows, 1536)
    image_grid_thw: (N, 3) int64 on the CPU, from (H, W) only (the processor reads it per image on the host)

mean' and std' come from the processor's own `_fuse_mean_std_and_rescale_factor` as (3, 1, 1) tensors on the device.
They are never Python numbers: CUDA divides by a CPU scalar as a multiply by its reciprocal, which differs from a
true divide in 126 of the 256 uint8 values when the divisor is 127.5 (Qwen3-VL-2B: mean = std = 0.5).

backend "hf" is the processor's own torch path on a device tensor; "torch" is the same ops written here.
`pick_backend` returns the first backend whose output equals the processor's CPU output on given images.
"""

import numpy as np
import torch

BACKENDS = ("hf", "torch")
CHECKED = {"patch_size": 16, "merge_size": 2, "temporal_patch_size": 2}  # Qwen3-VL's; others are not checked here


class UnsupportedImage(ValueError):
    """An image or a processor setting outside what patchify reproduces exactly; the caller uses the processor."""


def _array(x) -> np.ndarray:
    """(H, W, 3) or (H, W) uint8 of one image (PIL RGB / L, array or tensor)."""
    if isinstance(x, torch.Tensor):
        x = x.numpy()
    elif not isinstance(x, np.ndarray) and x.mode not in ("RGB", "L"):  # a PIL image
        raise UnsupportedImage(f"image mode {x.mode}: only RGB and L")
    a = np.asarray(x)
    if a.dtype != np.uint8 or a.ndim not in (2, 3) or (a.ndim == 3 and a.shape[2] != 3):
        raise UnsupportedImage(f"image {a.dtype} {a.shape}: needs uint8 (H, W, 3) or (H, W)")
    return a


def _check(ip, h: int, w: int) -> None:
    from transformers.models.qwen2_vl.image_processing_qwen2_vl import smart_resize

    for k, v in CHECKED.items():
        if getattr(ip, k, None) != v:
            raise UnsupportedImage(f"processor {k} {getattr(ip, k, None)}, checked only {v}")
    if not (ip.do_rescale and ip.do_normalize):
        raise UnsupportedImage("processor without rescale + normalize: not the fused path")
    if ip.do_resize and smart_resize(h, w, CHECKED["patch_size"] * CHECKED["merge_size"],
                                     ip.size["shortest_edge"], ip.size["longest_edge"]) != (h, w):  # fmt: skip
        raise UnsupportedImage(f"{h} x {w}: the processor would resize it")


def fused(ip, device) -> tuple[torch.Tensor, torch.Tensor]:
    """(mean', std') float32 (3, 1, 1) on `device`, from the processor's own fusion of rescale and normalize."""
    mean, std, _ = ip._fuse_mean_std_and_rescale_factor(
        do_normalize=True, image_mean=tuple(ip.image_mean), image_std=tuple(ip.image_std), do_rescale=True,
        rescale_factor=ip.rescale_factor, device=torch.device(device))  # fmt: skip
    return mean.float().view(-1, 1, 1), std.float().view(-1, 1, 1)


def _patchify(x: torch.Tensor, p: int, m: int, t: int) -> torch.Tensor:
    """HF Qwen2-VL patchify of (N, C, H, W) float: (N * grid_h * grid_w, C * t * p * p)."""
    n, c, h, w = x.shape
    gh, gw = h // p, w // p
    x = x.reshape(n, c, gh // m, m, p, gw // m, m, p).permute(0, 2, 5, 3, 6, 1, 4, 7)
    return x.unsqueeze(6).expand(-1, -1, -1, -1, -1, -1, t, -1, -1).reshape(n * gh * gw, c * t * p * p)


def patchify(images, image_processor, *, device=None, out_dtype=torch.float32, backend="torch"):
    """(pixel_values (sum t*h*w, 1536) on `device`, image_grid_thw (N, 3) int64 on the CPU), equal to
    image_processor(list(images), return_tensors="pt"). images: a uint8 tensor (N, H, W, 3) / (N, H, W), or a
    sequence of PIL RGB / L images or uint8 (H, W, 3) / (H, W) arrays; a grey image becomes R = G = B, as PIL's
    convert("RGB") does. Raises UnsupportedImage for another mode, a size the processor would resize, or unchecked
    processor settings."""
    if backend not in BACKENDS:
        raise ValueError(f"backend {backend!r}: one of {BACKENDS}")
    ip, device = image_processor, torch.device(device or "cpu")
    arrays = [_array(x) for x in images]
    for a in arrays:
        _check(ip, *a.shape[:2])
    p, m, t = (CHECKED[k] for k in ("patch_size", "merge_size", "temporal_patch_size"))
    grid = torch.tensor([[1, a.shape[0] // p, a.shape[1] // p] for a in arrays], dtype=torch.int64)
    groups: dict[tuple, list[int]] = {}
    for i, a in enumerate(arrays):
        groups.setdefault(a.shape, []).append(i)
    if backend == "torch":
        mean, std = fused(ip, device)
    out = [None] * len(arrays)
    for shape, idx in groups.items():
        host = torch.empty((len(idx), *shape), dtype=torch.uint8, pin_memory=device.type == "cuda")
        for j, i in enumerate(idx):
            host[j].numpy()[...] = arrays[i]
        x = host.to(device, non_blocking=True)
        if x.ndim == 3:  # grey: R = G = B on the device
            x = x[..., None].expand(-1, -1, -1, 3)
        x = x.permute(0, 3, 1, 2)
        if backend == "hf":
            pv = ip(list(x.contiguous()), return_tensors="pt", device=device)["pixel_values"]
        else:
            pv = _patchify(x.float().sub_(mean).div_(std), p, m, t)
        rows = pv.split([int(grid[i].prod()) for i in idx])
        for j, i in enumerate(idx):
            out[i] = rows[j]
    return torch.cat(out).to(out_dtype), grid


def pick_backend(images, image_processor, device) -> str | None:
    """The first of BACKENDS whose patchify equals the processor's CPU output (both tensors, bit for bit) on
    `images`, or None (then the caller keeps the CPU processor)."""
    import PIL.Image

    pil = [x if isinstance(x, PIL.Image.Image) else PIL.Image.fromarray(_array(x)) for x in images]
    want = image_processor([x.convert("RGB") for x in pil], return_tensors="pt")
    for b in BACKENDS:
        try:
            pv, grid = patchify(images, image_processor, device=device, backend=b)
        except UnsupportedImage:
            return None
        except Exception:  # a backend that does not run here (e.g. a processor without a device argument)
            continue
        if torch.equal(pv.cpu(), want["pixel_values"]) and torch.equal(grid, want["image_grid_thw"]):
            return b
    return None
