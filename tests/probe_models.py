"""Synthetic probe answerers and helpers used only by the tests (moved out of objdel.eval.probes): exact, noisy,
shifted, centred-prior and hole-free (A0) answers to a probe record, and the image / mask hashing helpers."""

import cv2
import numpy as np

from objdel.eval import probes as P
from objdel.geom.text import GRID, parse
from objdel.geom.topology import min_area
from objdel.targets import build_target


def image_sha(path):
    return P.pixels_sha256(cv2.cvtColor(cv2.imread(str(path), cv2.IMREAD_COLOR), cv2.COLOR_BGR2RGB))


def mask_text(m):
    return build_target({"nodule": m}, P.NAMES, min_area={"nodule": min_area(int(m.sum()))}).text


def parts(rec):
    return parse(rec["targets"][0], P.NAMES).objs["nodule"]


def model_exact(rec, rng):
    return rec["targets"][0]


def model_noise(sigma):
    def f(rec, rng):
        k = GRID / rec["W"]
        jit = lambda a: np.clip(np.rint(a + rng.normal(0, sigma * k, a.shape)).astype(np.int64), 0, GRID)  # noqa: E731
        return P._text([(jit(o), [jit(h) for h in hs]) for o, hs in parts(rec)])

    return f


def model_shift(rec, rng):
    k = GRID // rec["W"]
    dx, dy = [(k, 0), (-k, 0), (0, k), (0, -k)][rng.integers(4)]
    return P._text([(o + [dx, dy], [h + [dx, dy] for h in hs]) for o, hs in parts(rec)])


def model_centred_prior(prior_diam):
    def f(rec, rng):
        S = rec["W"]
        return P._disc_text(S / 2, S / 2, prior_diam[min(prior_diam, key=lambda s: abs(s - S))] / 2, S)

    return f


def model_a0(rec, rng):
    return P._text([(o, []) for o, _ in parts(rec)])

