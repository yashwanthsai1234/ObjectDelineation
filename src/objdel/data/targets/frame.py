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
- masks: any nonzero value or byte is foreground (objdel.geom.topology.as_mask), so the frame's masks, which targets and
  reward caches read, are canonical bools (issue I1).

Frame semantics: masks are (k, H, W) with W = H = S; the reward rasterizes at W x H. A grid vertex g maps to frame
pixel g * W / 1000 and to item pixel box[0] + scale * g * W / 1000 (x; y likewise with box[1]). box = (x0, y0, x1,
y1) is half-open in item pixels, x1 - x0 = scale * W, and may extend past the image (padding).

Views (Phase 2, D10; plan §5.0): frame(item, offset=(dx, dy), side=S, window=(lo, hi)) moves the scale-1 window of
side S (default: the policy's S) by (dx, dy) item pixels from where the policy would centre a window of that side,
and windows Hounsfield units with (lo, hi) instead of LUNG_WINDOW. A moved or resized window must hold the rater
union with >= CLEARANCE background pixels on every side and meet no other outlined nodule (meta["other_nodules"],
B3-2), else ValueError; the defaults reproduce the policy frame byte for byte and skip those checks.

References:
    archive/phase1_lab/t_targets/REPORT.md (crop policy H6) and targets/crops.py (pad_div, block reduce)
    M. Chen, Y. Peng, "Topology-preserving downsampling of binary images", ECCV 2024: https://arxiv.org/abs/2407.17786
    Pillow resampling (bicubic a = -0.5, centre = (xx + 0.5) * scale, PRECISION_BITS = 22 fixed point for 8-bit):
        https://github.com/python-pillow/Pillow/blob/main/src/libImaging/Resample.c ,
        https://pillow.readthedocs.io/en/stable/reference/Image.html#PIL.Image.Image.resize
    OpenCV INTER_CUBIC (A = -0.75, same half-pixel centres), replaced:
        https://github.com/opencv/opencv/blob/4.x/modules/imgproc/src/resize.cpp
    experiments/phase1/impl/M9_x86.md (finding 3: cv2 vs Pillow pixels on x86 and arm64)
    archive/phase0/baseline_medsam.py (LUNG_WINDOW = (-1000, 400))
"""
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
from PIL import Image

from objdel.data.store import Item
from objdel.geom.box import boxes_meet
from objdel.geom.topology import as_mask

DIVISORS = tuple(d for d in range(1, 1001) if 1000 % d == 0)
LUNG_WINDOW = (-1000, 400)
CLEARANCE = 2  # px between the rater union and the edge of a moved or resized window (D10)


@dataclass(frozen=True)
class FrameConfig:
    """Frame policy: largest and smallest side and margins (pixels, percent) and the model image side."""
    max_side: int = 500  # a divisor of 1000
    min_side: int = 50
    margin_px: int = 8
    margin_pct: int = 25
    image_side: int | None = None  # upscale the model image to this side (LIDC), else the frame side


FRAMES = {"riga": FrameConfig(), "lidc": FrameConfig(image_side=512)}


class Frame(NamedTuple):
    """A framed item: model image, bool masks per structure, side W = H in frame pixels, item-pixel box and scale."""
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


def union_box(item: Item) -> tuple[int, int, int, int]:
    """Half-open (x0, y0, x1, y1) box of every rater's union over all structures; ValueError when empty."""
    union = np.logical_or.reduce([as_mask(s, 3).any(0) for s in item.masks.values()])
    ys, xs = np.nonzero(union.any(1)), np.nonzero(union.any(0))
    if not ys[0].size:
        raise ValueError(f"{item.item_id}: no rater drew anything")
    return int(xs[0][0]), int(ys[0][0]), int(xs[0][-1]) + 1, int(ys[0][-1]) + 1


def centred(item: Item, side: int, box=None) -> tuple[int, int]:
    """Top-left corner of the policy's window of `side` px around the union (shifted inside the image if it fits)."""
    x0, y0, x1, y1 = box or union_box(item)
    Hi, Wi = next(iter(item.masks.values())).shape[-2:]
    return int(_window(x0, x1, side, Wi)), int(_window(y0, y1, side, Hi))


def crowded(item: Item, box) -> bool:
    """Does the half-open window box meet another outlined nodule (objdel.data.lidc.crowded, B3-2)?"""
    return boxes_meet(box, item.meta.get("other_nodules", ()))


def pixels(image: np.ndarray, x0: int, y0: int, side: int, window=None, noise=None) -> np.ndarray:
    """(c, side, side) crop at (x0, y0) of an (H, W[, c]) image, zero (black) outside it, channels first, alpha
    dropped. Hounsfield units (any dtype but uint8) get `noise` (side, side) added, then the HU window (default
    LUNG_WINDOW), inside the image only. Shared by frame and augment.render, so the two cannot drift."""
    chw = _crop(np.atleast_3d(image)[..., :3].transpose(2, 0, 1), x0, y0, side)
    if image.dtype == np.uint8:
        return chw
    H, W = image.shape[:2]
    y, x = np.arange(y0, y0 + side), np.arange(x0, x0 + side)
    inside = ((y >= 0) & (y < H))[:, None] & ((x >= 0) & (x < W))
    chw = chw.astype(np.int64) + (0 if noise is None else noise)
    return np.where(inside, hu_window(chw, LUNG_WINDOW if window is None else window), 0)


def rgb(chw: np.ndarray, side: int | None) -> np.ndarray:
    """(1 | 3, S, S) 0..255 -> uint8 (S, S, 3), upscaled to `side`."""
    return upscale(np.ascontiguousarray(np.repeat(chw, 3 // len(chw), 0).transpose(1, 2, 0), np.uint8), side)


def hu_window(a: np.ndarray, window=LUNG_WINDOW) -> np.ndarray:
    """Hounsfield units -> 0..255 (int64), rounded half up: the lung window of Phase 0 for the default; ValueError
    for an empty or inverted window (break-it BK1-6)."""
    lo, hi = window
    if hi <= lo:
        raise ValueError(f"HU window {window}: hi must exceed lo")
    return ((np.clip(a.astype(np.int64), lo, hi) - lo) * 255 + (hi - lo) // 2) // (hi - lo)


def upscale(image: np.ndarray, side: int | None) -> np.ndarray:
    """uint8 (S, S, 3) -> (side, side, 3) with Pillow bicubic (bit-identical on every CPU, F3); no-op for None / S."""
    if not side or side == image.shape[0]:
        return image
    return np.array(Image.fromarray(image).resize((side, side), Image.Resampling.BICUBIC))  # writable


def policy_side(box, cfg: FrameConfig) -> tuple[int, int]:
    """(f, S) of the policy frame of a union box: reduction factor and frame side (module docstring)."""
    extent = max(box[2] - box[0], box[3] - box[1])
    need = extent + 2 * max(cfg.margin_px, -(-extent * cfg.margin_pct // 100))
    f = -(-need // cfg.max_side)
    return f, min(d for d in DIVISORS if d >= max(-(-need // f), cfg.min_side))


def frame(item: Item, cfg: FrameConfig | None = None, *, offset=(0, 0), side: int | None = None,
          window=None) -> Frame:
    """Frame an item: the window, model image and masks in frame pixels (module docstring).

    Raises ValueError for a `side` that is not a divisor of 1000 in range, or when a moved window leaves the union
    without clearance or meets another outlined nodule.
    """
    cfg = cfg or FRAMES[item.dataset]
    stacks = {n: as_mask(s, 3) for n, s in item.masks.items()}
    ux0, uy0, ux1, uy1 = ub = union_box(item)
    f, S = policy_side(ub, cfg)
    dx, dy = (int(v) for v in offset)
    moved = (dx, dy) != (0, 0) or side is not None and (side, f) != (S, 1)
    if side is not None:
        if side not in DIVISORS or not cfg.min_side <= side <= cfg.max_side:
            raise ValueError(f"{item.item_id}: side {side} is not a divisor of 1000 in "
                             f"[{cfg.min_side}, {cfg.max_side}]")
        f, S = 1, int(side)
    x0, y0 = centred(item, f * S, ub)
    x0, y0 = x0 + dx, y0 + dy
    if moved:
        c = CLEARANCE
        if not (ux0 - x0 >= c and uy0 - y0 >= c and x0 + S - ux1 >= c and y0 + S - uy1 >= c):
            raise ValueError(f"{item.item_id}: the union {ub} is not inside the window {(x0, y0, x0 + S, y0 + S)} "
                             f"with {c} px clearance")
        if crowded(item, (x0, y0, x0 + S, y0 + S)):
            raise ValueError(f"{item.item_id}: the window {(x0, y0, x0 + S, y0 + S)} meets another outlined nodule")
    masks = {n: _crop(s, x0, y0, f * S) for n, s in stacks.items()}
    if f > 1:
        masks = {n: 2 * _blocks(m, f) >= f * f for n, m in masks.items()}
    image = item.image
    if image is not None:
        if image.ndim not in (2, 3) or image.ndim == 3 and image.shape[2] not in (1, 3, 4):
            raise ValueError(f"{item.item_id}: image must be gray, RGB or RGBA (H, W[, 1|3|4]), got {image.shape}")
        image = pixels(image, x0, y0, f * S, window)
        if f > 1:
            image = (_blocks(image, f) + f * f // 2) // (f * f)  # block mean, rounded half up
        image = rgb(image, cfg.image_side)  # Pillow, not cv2.INTER_CUBIC: bit-identical on every CPU (F3)
    return Frame(image, masks, S, S, (x0, y0, x0 + f * S, y0 + f * S), f)
