"""SFT targets and the model frame.

frame:    item -> Frame: the square crop (side <= 500, a divisor of 1000) the model sees and the reward scores.
build:    one rater's masks on a frame -> Target: certified (checker), canonical polygon text; the reward's size rule
          per structure (min_areas); target-set hash.
augment:  views: D4 ops, window draws (shift, scale), CT window jitter, context crops, render.
simplify: min-# ring simplification (potrace polygon stage OR lossless digital straight segments), numba kernels.
"""
