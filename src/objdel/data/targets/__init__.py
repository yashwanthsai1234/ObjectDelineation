"""SFT targets and the model frame (plan D11, §4.4).

frame:    item -> Frame: the square crop (side <= 500, a divisor of 1000) the model sees and the reward scores.
build:    one rater's masks on a frame -> Target: certified (checker), canonical polygon text; the reward's size rule
          per structure (min_areas); target-set hash.
augment:  Phase 2 views: D4 ops, window draws (shift, scale), CT window jitter, context crops, render (C-VIEW).
simplify: min-# ring simplification (potrace polygon stage OR lossless digital straight segments), numba kernels.
"""
