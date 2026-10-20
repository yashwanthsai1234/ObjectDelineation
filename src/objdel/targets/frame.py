"""The model's view of an item: a square crop around every rater's union, at the resolution the reward scores.

Policy (plan §4.4, lab t_targets H6):
- side s = the larger side of the union box (all raters, all structures), plus a margin of
  max(margin_px, ceil(margin_pct % of s)) on each side;
- integer factor f = ceil(need / max_side); the frame side S is the smallest divisor of 1000 that is >= the reduced
  need and >= min_side, so S <= 500 and one pixel is exactly 1000 / S grid units: every pixel corner is a grid point
  and the target round trip mask -> polygon -> mask has no rounding error (lab H6: crops whose side does not divide
  1000 flipped pixels and dropped LIDC IoU to 0.874). The context grows to reach S; the object is never resampled
  for that;
- the crop window of f*S item pixels is centred on the union box and shifted inside the image when it fits, else
  padded (masks False, image black);
- f > 1: masks are reduced by f x f block majority (ties, possible for even f, go to the foreground) and this reduced
  mask is the ground truth for SFT and reward alike; the image by the f x f block mean. Block voting is monotone, so
  cup inside disc stays inside. It can change topology (lab: 7.8 % of rims when forced); a topology-preserving
  downsampler exists (Chen & Peng 2024, integer programming) but is not needed because the reduced mask is the
  ground truth everywhere. None of the RIGA or LIDC items needs f > 1 (largest union 309 px);
- the model image is uint8 RGB. LIDC Hounsfield units go through the lung window [-1000, 400] (Phase 0) and the
  uint8 crop is upscaled to image_side (512) with Pillow's bicubic filter; masks, W and H stay at the frame
  resolution. Pillow resamples 8-bit images in fixed point (22-bit integer coefficients, integer accumulation), so
  the image is bit-identical across CPUs, while cv2.resize(INTER_CUBIC) differed by +-1 between x86 (with or without
  IPP) and arm64 (M9 finding 3, fix F3). Both sample at pixel centres (output centre (x + 1/2) * in / out), so the
  image lines up with the mask scaled by image_side / S; the kernels differ (Keys a = -0.5 in Pillow, -0.75 in
  OpenCV), which changes pixel values by a few levels but not the geometry. Accepted image layouts:
  gray (H, W) or (H, W, 1), repeated to RGB; RGB (H, W, 3); RGBA (H, W, 4), whose alpha channel is dropped (the
  stores hold opaque photos and CT slices); anything else raises ValueError (break-it B1-6);
- masks: any nonzero value or byte is foreground (objdel.checker.as_mask), so the frame's masks, which targets and
  reward caches read, are canonical bools (issue I1).

Frame semantics: masks are (k, H, W) with W = H = S; the reward rasterizes at W x H. A grid vertex g maps to frame
pixel g * W / 1000 and to item pixel box[0] + scale * g * W / 1000 (x; y likewise with box[1]). box = (x0, y0, x1,
y1) is half-open in item pixels, x1 - x0 = scale * W, and may extend past the image (padding).

References:
    experiments/phase1/lab/t_targets/REPORT.md (crop policy H6) and targets/crops.py (pad_div, block reduce)
    M. Chen, Y. Peng, "Topology-preserving downsampling of binary images", ECCV 2024: https://arxiv.org/abs/2407.17786
    Pillow resampling (bicubic a = -0.5, centre = (xx + 0.5) * scale, PRECISION_BITS = 22 fixed point for 8-bit):
        https://github.com/python-pillow/Pillow/blob/main/src/libImaging/Resample.c ,
        https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.Image.resize
    OpenCV INTER_CUBIC (A = -0.75, same half-pixel centres), replaced:
        https://github.com/opencv/opencv/blob/4.x/modules/imgproc/src/resize.cpp
    experiments/phase1/impl/M9_x86.md (finding 3: cv2 vs Pillow pixels on x86 and arm64)
    experiments/phase0/baseline_medsam.py (LUNG_WINDOW = (-1000, 400))
"""
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
from PIL import Image

from objdel.checker import as_mask
from objdel.data import Item

DIVISORS = tuple(d for d in range(1, 1001) if 1000 % d == 0)
LUNG_WINDOW = (-1000, 400)


@dataclass(frozen=True)
class FrameConfig:
    max_side: int = 500  # a divisor of 1000
    min_side: int = 50
    margin_px: int = 8
    margin_pct: int = 25
    image_side: int | None = None  # upscale the model image to this side (LIDC), else the frame side


FRAMES = {"riga": FrameConfig(), "lidc": FrameConfig(image_side=512)}


class Frame(NamedTuple):
    image: np.ndarray | None  # uint8 (h, w, 3) model input, or None when the item has no pixels
    masks: dict[str, np.ndarray]  # bool (k, H, W) per structure: the ground truth of targets and reward
    W: int
    H: int
    box: tuple[int, int, int, int]  # x0, y0, x1, y1 in item pixels, half-open
    scale: int  # item pixels per frame pixel


def _window(lo, hi, side, size):
    """Start of a window of `side` px around [lo, hi): centred, shifted inside [0, size) when it fits."""
    start = (lo + hi - side) // 2
    return min(max(start, 0), size - side) if side <= size else start


def _crop(a, x0, y0, side):
    """a[..., y0:y0+side, x0:x0+side] over the last two axes, zero-padded outside a."""
    H, W = a.shape[-2:]
    out = np.zeros(a.shape[:-2] + (side, side), a.dtype)
    ys, xs = slice(max(y0, 0), min(y0 + side, H)), slice(max(x0, 0), min(x0 + side, W))
    if ys.start < ys.stop and xs.start < xs.stop:
        out[..., ys.start - y0:ys.stop - y0, xs.start - x0:xs.stop - x0] = a[..., ys, xs]
    return out


def _blocks(a, f):
    """Sums over f x f blocks of the last two axes."""
    *lead, h, w = a.shape
    return a.reshape(*lead, h // f, f, w // f, f).sum((-3, -1), dtype=np.int64)


def frame(item: Item, cfg: FrameConfig | None = None) -> Frame:
    cfg = cfg or FRAMES[item.dataset]
    stacks = {n: as_mask(s, 3) for n, s in item.masks.items()}
    union = np.logical_or.reduce([s.any(0) for s in stacks.values()])
    ys, xs = np.nonzero(union)
    if not ys.size:
        raise ValueError(f"{item.item_id}: no rater drew anything")
    side = int(max(ys[-1] - ys[0], xs.max() - xs.min())) + 1
    need = side + 2 * max(cfg.margin_px, -(-side * cfg.margin_pct // 100))
    f = -(-need // cfg.max_side)
    S = min(d for d in DIVISORS if d >= max(-(-need // f), cfg.min_side))
    Hi, Wi = union.shape
    x0, y0 = int(_window(xs.min(), xs.max() + 1, f * S, Wi)), int(_window(ys[0], ys[-1] + 1, f * S, Hi))
    masks = {n: _crop(s, x0, y0, f * S) for n, s in stacks.items()}
    if f > 1:
        masks = {n: 2 * _blocks(m, f) >= f * f for n, m in masks.items()}
    image = item.image
    if image is not None:
        if image.ndim not in (2, 3) or image.ndim == 3 and image.shape[2] not in (1, 3, 4):
            raise ValueError(f"{item.item_id}: image must be gray, RGB or RGBA (H, W[, 1|3|4]), got {image.shape}")
        if image.dtype != np.uint8:  # Hounsfield units
            lo, hi = LUNG_WINDOW
            image = ((np.clip(image.astype(np.int64), lo, hi) - lo) * 255 + (hi - lo) // 2) // (hi - lo)
        image = _crop(np.atleast_3d(image)[..., :3].transpose(2, 0, 1), x0, y0, f * S)  # channels first, no alpha
        if f > 1:
            image = (_blocks(image, f) + f * f // 2) // (f * f)  # block mean, rounded half up
        image = np.ascontiguousarray(np.repeat(image, 3 // len(image), 0).transpose(1, 2, 0), np.uint8)  # gray: x3
        if cfg.image_side and cfg.image_side != S:  # Pillow, not cv2.INTER_CUBIC: bit-identical on every CPU (F3)
            image = np.array(Image.fromarray(image).resize((cfg.image_side,) * 2, Image.Resampling.BICUBIC))  # writable
    return Frame(image, masks, S, S, (x0, y0, x0 + f * S, y0 + f * S), f)
