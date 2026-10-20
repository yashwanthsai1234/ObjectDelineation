"""SFT targets and the model frame (plan D11, §4.4).

frame:    item -> Frame: the square crop (side <= 500, a divisor of 1000) the model sees and the reward scores.
build:    one rater's masks on a frame -> Target: certified (polycore), canonical polygon text; the reward's size rule
          per structure (min_areas); target-set hash.
simplify: min-# ring simplification (potrace polygon stage OR lossless digital straight segments), numba kernels.
"""
from .build import BuildConfig, Target, build_target, min_areas, target_set_hash
from .frame import FRAMES, Frame, FrameConfig, frame

__all__ = [
    "FRAMES", "BuildConfig", "Frame", "FrameConfig", "Target", "build_target", "frame", "min_areas", "target_set_hash",
]
