"""One box-overlap rule: used by the LIDC loader and by view selection."""
from collections.abc import Iterable


def boxes_meet(box, others: Iterable) -> bool:
    """Does the half-open window box = (x0, y0, x1, y1) overlap any (a, b, c, d) of others?"""
    x0, y0, x1, y1 = box
    return any(a < x1 and x0 < c and b < y1 and y0 < d for a, b, c, d in others)
