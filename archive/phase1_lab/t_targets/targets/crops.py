"""Crop policies for target building (all exact: padding only, no resampling).

native : the fixture crop as given (W x H arbitrary, possibly > 1000 px).
pad_sq : pad to a square of side max(W, H) (arbitrary side).
pad_div: pad to the smallest square side S >= max(W, H) with 1000 % S == 0 (H6 divisor crop);
         crops larger than 500 px are first reduced by an integer box factor f (block majority, ties
         -> foreground) so that S <= 500. (Only used when an object does not fit 500 px.)
resize500 / resize250: nearest-neighbour resize of the padded square to 500 / 250 (lossy reference).
"""
import numpy as np
from PIL import Image

DIVISORS = [d for d in range(1, 1001) if 1000 % d == 0]


def smallest_divisor_at_least(n):
    return min(d for d in DIVISORS if d >= n)


def _pad(m, S):
    H, W = m.shape
    out = np.zeros((S, S), bool)
    y0, x0 = (S - H) // 2, (S - W) // 2
    out[y0:y0 + H, x0:x0 + W] = m
    return out


def _block_reduce(m, f):
    H, W = m.shape
    Hp, Wp = -(-H // f) * f, -(-W // f) * f
    p = np.zeros((Hp, Wp), np.int32)
    p[:H, :W] = m
    s = p.reshape(Hp // f, f, Wp // f, f).sum(axis=(1, 3))
    return s * 2 >= f * f


def apply_crop(masks, policy):
    H, W = next(iter(masks.values())).shape
    if policy == "native":
        return masks, (H, W)
    n = max(H, W)
    if policy == "pad_sq":
        return {k: _pad(v, n) for k, v in masks.items()}, (n, n)
    if policy == "pad_div":
        f = 1
        while -(-n // f) > 500:
            f += 1
        if f > 1:
            masks = {k: _block_reduce(v, f) for k, v in masks.items()}
            H, W = next(iter(masks.values())).shape
            n = max(H, W)
        S = smallest_divisor_at_least(n)
        return {k: _pad(v, S) for k, v in masks.items()}, (S, S)
    if policy.startswith("resize"):
        S = int(policy[6:])
        out = {}
        for k, v in masks.items():
            sq = _pad(v, n)
            out[k] = np.asarray(Image.fromarray(sq.astype(np.uint8) * 255).resize((S, S), Image.NEAREST)) > 127
        return out, (S, S)
    raise ValueError(policy)
