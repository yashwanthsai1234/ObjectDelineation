"""Polygon text from a model -> validated rings -> raster mask."""
import json
import re

import numpy as np
from PIL import Image, ImageDraw


def extract_json(text):
    """First JSON object in a model's reply, or None. Tolerates ```json fences."""
    text = re.sub(r"```(?:json)?", "", text)
    start = text.find("{")
    while start != -1:
        depth = 0
        for i in range(start, len(text)):
            depth += {"{": 1, "}": -1}.get(text[i], 0)
            if depth == 0:
                try:
                    return json.loads(text[start:i + 1])
                except json.JSONDecodeError:
                    break
        start = text.find("{", start + 1)
    return None


def as_points(value):
    """[[x, y], ...] or [x, y, x, y, ...] -> (n, 2) float array, or None."""
    try:
        a = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        return None
    if a.ndim == 1 and a.size % 2 == 0:
        a = a.reshape(-1, 2)
    if a.ndim != 2 or a.shape[1] != 2:
        return None
    return a


def _cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _on_segment(p, q, r):
    return min(p[0], q[0]) <= r[0] <= max(p[0], q[0]) and min(p[1], q[1]) <= r[1] <= max(p[1], q[1])


def _segments_meet(p1, p2, p3, p4):
    """'cross' for a proper crossing, 'touch' for touching or collinear overlap, None if apart."""
    d1, d2 = _cross(p3, p4, p1), _cross(p3, p4, p2)
    d3, d4 = _cross(p1, p2, p3), _cross(p1, p2, p4)
    if d1 * d2 < 0 and d3 * d4 < 0:
        return "cross"
    for d, a, b, c in ((d1, p3, p4, p1), (d2, p3, p4, p2), (d3, p1, p2, p3), (d4, p1, p2, p4)):
        if d == 0 and _on_segment(a, b, c):
            return "touch"
    return None


def ring_problems(pts, bounds=None):
    """Reasons a ring is invalid; empty list if it is a simple polygon.

    bounds = (lo, hi) flags coordinates outside the allowed range.
    """
    if pts is None:
        return ["unparseable"]
    probs = []
    if bounds is not None and (pts.min() < bounds[0] or pts.max() > bounds[1]):
        probs.append("out_of_range")
    p = pts[np.r_[True, np.any(np.diff(pts, axis=0) != 0, axis=1)]]
    if len(p) > 1 and np.all(p[0] == p[-1]):
        p = p[:-1]
    if len(p) < 3:
        return probs + ["fewer_than_3_points"]
    area2 = np.sum(p[:, 0] * np.roll(p[:, 1], -1) - np.roll(p[:, 0], -1) * p[:, 1])
    if not area2:
        probs.append("zero_area")
    n = len(p)
    meet = set()
    for i in range(n):
        for j in range(i + 2, n):
            if i == 0 and j == n - 1:
                continue
            m = _segments_meet(p[i], p[(i + 1) % n], p[j], p[(j + 1) % n])
            if m:
                meet.add(m)
    if "cross" in meet:
        probs.append("self_crossing")
    elif "touch" in meet:
        probs.append("self_touching")
    return probs


def rasterize(pts, width, height):
    """Fill a ring given in pixel coordinates (outline included)."""
    img = Image.new("1", (width, height), 0)
    if pts is not None and len(pts) >= 3:
        ImageDraw.Draw(img).polygon([tuple(map(float, q)) for q in pts], fill=1, outline=1)
    return np.asarray(img, dtype=bool)
