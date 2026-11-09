"""The pad_div crop of the phase 1 lab (archive/phase1_lab, t_targets/targets/crops.py), ported so tests need no lab.

Masks larger than 500 px are first reduced by the smallest integer factor f that makes them fit (block majority, a
tie counts as foreground). Then each mask is centred in a square whose side is the smallest divisor of 1000 that
is at least the larger mask side."""

import hashlib

import numpy as np

DIVISORS = [d for d in range(1, 1001) if 1000 % d == 0]


def _block_reduce(m, f):
    H, W = m.shape
    Hp, Wp = -(-H // f) * f, -(-W // f) * f
    p = np.zeros((Hp, Wp), np.int32)
    p[:H, :W] = m
    return p.reshape(Hp // f, f, Wp // f, f).sum(axis=(1, 3)) * 2 >= f * f


def _pad(m, S):
    H, W = m.shape
    out = np.zeros((S, S), bool)
    y0, x0 = (S - H) // 2, (S - W) // 2
    out[y0 : y0 + H, x0 : x0 + W] = m
    return out


def pad_div(masks):
    n = max(next(iter(masks.values())).shape)
    f = 1
    while -(-n // f) > 500:
        f += 1
    if f > 1:
        masks = {k: _block_reduce(v, f) for k, v in masks.items()}
        n = max(next(iter(masks.values())).shape)
    S = min(d for d in DIVISORS if d >= n)
    return {k: _pad(v, S) for k, v in masks.items()}


def mask_sha(masks):
    """SHA-256 hex of the shapes and bytes of the masks, in the order of the dict."""
    h = hashlib.sha256()
    for m in masks.values():
        h.update(str(m.shape).encode() + np.ascontiguousarray(m, bool).tobytes())
    return h.hexdigest()
