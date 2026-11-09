"""Independent largest-remainder (Hamilton) rounding for tests (plan B.9.0): plain Python, ties to the lowest index.

lrr(x, n): floor every entry, then give the n - sum(floors) leftover units to the entries with the largest
fractional parts; equal fractional parts go to the lower index first. Returns an int numpy array (so tests can
compare with `.tolist()` or scale it)."""
import math

import numpy as np


def lrr(x, n: int) -> np.ndarray:
    x = [float(v) for v in x]
    out = [math.floor(v) for v in x]
    order = sorted(range(len(x)), key=lambda i: (-(x[i] - out[i]), i))
    for i in order[:n - sum(out)]:
        out[i] += 1
    return np.array(out, dtype=np.int64)
