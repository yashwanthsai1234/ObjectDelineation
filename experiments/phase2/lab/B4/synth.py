"""Synthetic RIGA-like items for the K3/K4 timing harness on hosts that must not see RIGA data (H-K3SYN).

An item is a 500 x 500 frame with 6 raters' disc and cup masks: smooth closed outlines r(t) = R (1 + sum of low
harmonics) with a per-item shape, per-rater radius and shape perturbations, the cup inside the disc and offset from
its centre, filled with cv2.fillPoly at 1/16 px. The parameter ranges are the 5-95 % quantiles of aggregate
statistics of the 64-item RIGA K3 prompt set (riga_stats.py; no RIGA pixel, mask or id is used here):

    disc mean radius 94-124 px, height/width 0.99-1.20, centred within +-1 % of the frame; cup/disc radius ratio
    0.31-0.68, cup offset -0.13..0.11 R (x) and -0.07..0.07 R (y), cup height/width 0.77-1.18; rater radius s.d.
    1.8 % (disc) and 13 % (cup) of the mean; band width (|U| - |I|) / P: disc 7.6-15 px, cup 12-27 px; one
    component, no hole; completion length (make_jobs) 624-1,444 characters, median 814.

Only Python floats and math.cos/sin feed integer vertex coordinates (rounded to 1/16 px), so the corpus is meant to
be identical on arm64 and x86-64: `corpus_digest` is compared across both (criterion PS-5).

    from experiments.phase2.lab.B4.synth import corpus
    items = corpus(64, seed=0)        # [{"disc": (6, 500, 500) bool, "cup": ...}, ...]
"""
import hashlib
import math

import cv2
import numpy as np

H = W = 500
K = 6
N_VERT = 720  # outline samples before rasterization (smooth at 1/16 px)
DISC_HF, CUP_HF = 0.0055, 0.013  # high-harmonic roughness: sets the traced vertex count (completion length)


def _outline(cy, cx, r, asp, harm, n=N_VERT) -> np.ndarray:
    """(n, 2) int32 (x, y) at 1/16 px of r(t) = r (1 + sum_k a_k cos(k t + p_k)), y stretched by asp."""
    pts = []
    for i in range(n):
        t = 2.0 * math.pi * i / n
        rr = r * (1.0 + sum(a * math.cos(k * t + p) for k, a, p in harm))
        pts.append((round(16.0 * (cx + rr * math.cos(t))), round(16.0 * (cy + asp * rr * math.sin(t)))))
    return np.array(pts, np.int32)


def _fill(poly: np.ndarray) -> np.ndarray:
    m = np.zeros((H, W), np.uint8)
    cv2.fillPoly(m, [poly], 1, lineType=cv2.LINE_8, shift=4)
    return m.astype(bool)


def _harm(rng, ks, amp) -> list[tuple[int, float, float]]:
    return [(k, float(rng.normal(0.0, amp)), float(rng.uniform(0.0, 2.0 * math.pi))) for k in ks]


def item(rng) -> dict[str, np.ndarray]:
    """One synthetic item: {"disc": (K, H, W) bool, "cup": (K, H, W) bool}."""
    R = float(rng.uniform(87.0, 116.0))
    asp = float(rng.uniform(0.99, 1.20))
    cy, cx = H / 2 + float(rng.uniform(-4.0, 4.0)), W / 2 + float(rng.uniform(-4.0, 4.0))
    disc_shape = _harm(rng, (2, 3, 4), 0.02)
    rho = float(rng.uniform(0.29, 0.64))
    cdy, cdx = float(rng.uniform(-0.07, 0.07)) * R, float(rng.uniform(-0.13, 0.11)) * R
    casp = float(rng.uniform(0.77, 1.18))
    cup_shape = _harm(rng, (2, 3), 0.04)
    discs, cups = [], []
    for _ in range(K):
        rj = R * (1.0 + float(rng.normal(0.0, 0.018)))
        oy, ox = (float(v) for v in rng.normal(0.0, 1.5, 2))
        dh = [(k, a + float(rng.normal(0.0, 0.025)), p) for k, a, p in disc_shape] + _harm(rng, range(5, 21), DISC_HF)
        d = _fill(_outline(cy + oy, cx + ox, rj, asp, dh))
        rc = rho * R * (1.0 + float(rng.normal(0.0, 0.13)))
        ch = [(k, a + float(rng.normal(0.0, 0.05)), p) for k, a, p in cup_shape] + _harm(rng, range(4, 15), CUP_HF)
        c = _fill(_outline(cy + cdy + 2 * oy, cx + cdx + 2 * ox, rc, casp, ch)) & d
        discs.append(d)
        cups.append(c)
    return {"disc": np.stack(discs), "cup": np.stack(cups)}


def corpus(n: int, seed: int = 0) -> list[dict[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    return [item(rng) for _ in range(n)]


def corpus_digest(items) -> str:
    h = hashlib.sha256()
    for it in items:
        for name in ("disc", "cup"):
            h.update(np.packbits(it[name]).tobytes())
    return h.hexdigest()
