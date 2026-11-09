"""Independent pixel check of a train window for tests (plan B.9.0 `brute_admissible`; DG4, DG5, DG18).

A half-open window box (x0, y0, x1, y1) in slice pixels is admissible when it lies inside the slice (no padding),
every doctor pixel is inside it with >= 2 px of window on each side, and no pixel of another outlined nodule's box
(item.meta["other_nodules"]) is inside it. Pixel masks only: no box arithmetic shared with the code under test."""
import numpy as np


def brute_admissible(item, box) -> bool:
    x0, y0, x1, y1 = (int(v) for v in box)
    stack = next(iter(item.masks.values()))
    H, W = stack.shape[-2:]
    if x0 < 0 or y0 < 0 or x1 > W or y1 > H:
        return False
    inner = np.zeros((H, W), bool)
    inner[y0 + 2:y1 - 2, x0 + 2:x1 - 2] = True
    doctors = np.logical_or.reduce([np.asarray(s, bool).reshape(-1, H, W).any(0) for s in item.masks.values()])
    if (doctors & ~inner).any():
        return False
    window = np.zeros((H, W), bool)
    window[y0:y1, x0:x1] = True
    for a, b, c, d in item.meta.get("other_nodules", ()):
        other = np.zeros((H, W), bool)
        other[max(b, 0):max(d, 0), max(a, 0):max(c, 0)] = True
        if (other & window).any():
            return False
    return True
