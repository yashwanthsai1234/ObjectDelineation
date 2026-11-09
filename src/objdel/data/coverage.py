"""Train-vs-exam coverage: does the training view contain what an evaluation view asks (R1's rescaled gap)?

A row is a dict of categorical properties of one target outline (props). distribution() gives each property's share per
category over many rows; gaps() lists every evaluation category with a share >= EVAL_MIN that training has below
TRAIN_MAX. cell_coverage() compares scheduled windows with a (size, offset) cell law (B1.3, gates DG7-DG9, DG16).
Pure functions, no I/O.
"""

import json
from collections import Counter

import numpy as np
from scipy.stats import ks_2samp

from .windows import cell_shares

EVAL_MIN = 0.05  # an eval category this common ...
TRAIN_MAX = 0.01  # ... must have at least this much training share
KEYS = ("W", "ext", "pos", "topo", "verts", "grid10", "op")
GRID = 1000  # outline coordinates are integers 0..1000 of the crop


def _bucket(x, edges: tuple, names: tuple) -> str:
    return names[sum(x >= e for e in edges)]


def props(target: str, W: int, op: str) -> dict | None:
    """Properties of one target text ({"nodule": [{"outer", "holes"}]}), or None when it holds no outline.
    ext: longer side of the outline's box as a share of the crop; pos: how far its centre is from the crop centre."""
    parts = next(iter(json.loads(target).values()), [])
    pts = [p for part in parts for p in part["outer"]]
    if not pts:
        return None
    xs, ys = [p[0] for p in pts], [p[1] for p in pts]
    ext = max(max(xs) - min(xs), max(ys) - min(ys))
    off = max(abs((max(xs) + min(xs)) / 2 - GRID / 2), abs((max(ys) + min(ys)) / 2 - GRID / 2))
    coords = [c for part in parts for ring in [part["outer"], *part["holes"]] for p in ring for c in p]
    return {
        "W": W,
        "ext": _bucket(ext, (100, 150, 300, 600), ("<10%", "10-15%", "15-30%", "30-60%", ">=60%")),
        "pos": _bucket(off, (150, 300), ("centre", "middle", "edge")),
        "topo": "parts>1" if len(parts) > 1 else "hole" if any(p["holes"] for p in parts) else "simple",
        "verts": _bucket(len(pts), (8, 16, 31), ("<8", "8-15", "16-30", ">30")),
        "grid10": all(c % 10 == 0 for c in coords),
        "op": op,
    }


def distribution(rows: list[dict]) -> dict[str, dict]:
    """{property: {category: share of rows}}."""
    out = {}
    for k in KEYS:
        c = Counter(r[k] for r in rows)
        out[k] = {cat: n / len(rows) for cat, n in c.items()}
    return out


def gaps(train: dict, exam: dict) -> list[dict]:
    """Every (property, category) the exam has at >= EVAL_MIN whose training share is < TRAIN_MAX."""
    return [{"key": k, "category": cat, "train": train[k].get(cat, 0.0), "exam": share}
            for k in KEYS for cat, share in sorted(exam[k].items(), key=lambda kv: str(kv[0]))
            if share >= EVAL_MIN and train[k].get(cat, 0.0) < TRAIN_MAX]  # fmt: skip


def mix_coverage(rows, target, target_offsets) -> dict:
    """Worst cell ratio (cells with target share >= 1 %), total variation and offset KS (scipy.stats.ks_2samp) of the
    windows `rows` against the (4, 4) cell law `target` and its offsets; no rows: 0, 1, 1 (a failure, never NaN)."""
    if not rows:
        return {"rows": 0, "worst_ratio": 0.0, "tv": 1.0, "ks": 1.0}
    share, t = cell_shares(rows), np.asarray(target, np.float64)
    keep = t >= 0.01
    return {"rows": len(rows), "worst_ratio": float((share[keep] / t[keep]).min()),
            "tv": float(np.abs(share - t).sum() / 2),
            "ks": float(ks_2samp([w.offset for w in rows], target_offsets).statistic)}


def cell_coverage(scheduled, target, target_offsets) -> dict:
    """mix_coverage of the off-centre (shifted and rescaled) windows of `scheduled`; centred windows are ignored."""
    return mix_coverage([w for w in scheduled if w.kind != "centred"], target, target_offsets)
