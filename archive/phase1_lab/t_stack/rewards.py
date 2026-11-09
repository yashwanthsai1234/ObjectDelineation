"""objdel-backed reward functions in TRL's batched signature: f(completions, **columns) -> list[float | None].

Ground truth is looked up by the dataset's `key` column in a module-level registry (masks are too big for a
datasets column). Coordinates are 0-1000 relative to the image, as in the SFT target format.
"""
import numpy as np

from objdel.checker.polygon_io import as_points, extract_json, rasterize, ring_problems
from objdel.reward import AnnotatorSet

GT = {}  # key -> AnnotatorSet of the side x side disc masks
SIDE = {}  # key -> image side in px


def register(key, masks):
    GT[key] = AnnotatorSet(masks)
    SIDE[key] = masks[0].shape[0]


def _text(c):
    # conversational completions are [{"role": "assistant", "content": "..."}]
    if isinstance(c, list):
        c = c[0]["content"] if c else ""
        if isinstance(c, list):  # content parts
            c = "".join(p.get("text", "") for p in c if isinstance(p, dict))
    return c


def parse(text):
    """-> (outer points (n,2) in 0-1000 or None, list of problems)."""
    obj = extract_json(text)
    if not isinstance(obj, dict) or not isinstance(obj.get("disc"), list) or not obj["disc"]:
        return None, ["unparseable"]
    part = obj["disc"][0]
    if not isinstance(part, dict):
        return None, ["bad_schema"]
    pts = as_points(part.get("outer"))
    return pts, ring_problems(pts, bounds=(0, 1000))


def _score(text, key):
    pts, probs = parse(text)
    if pts is None or len(pts) < 3:
        return None, probs
    side = SIDE[key]
    m = rasterize(pts * (side - 1) / 1000, side, side)
    return GT[key].score(m), probs


def format_reward(completions, key, **kw):
    """1 = parses into the schema as a simple polygon, 0.5 = parses but invalid ring, 0 = unparseable."""
    out = []
    for c in completions:
        pts, probs = parse(_text(c))
        out.append(0.0 if pts is None else (1.0 if not probs else 0.5))
    return out


def iou_reward(completions, key, **kw):
    """IoU against the experts' strict-majority mask (0 if unparseable)."""
    scores = [_score(_text(c), k)[0] for c, k in zip(completions, key)]
    return [s.iou if s is not None else 0.0 for s in scores]


def topo_reward(completions, key, gate=0.5, **kw):
    """topo_near, gated on IoU >= gate (as objdel.reward.group_advantages does)."""
    out = []
    for c, k in zip(completions, key):
        s, _ = _score(_text(c), k)
        out.append(float(s.topo_near) if s is not None and s.iou >= gate else 0.0)
    return out


def legality_const(completions, **kw):
    """Constant 1.0: stands in for a legality term that is ~always satisfied (NaN / zero-variance check)."""
    return [1.0] * len(completions)


def gdpo_like(rewards_per_func, weights, G, eps=1e-4):
    """Re-implementation of TRL 1.16 normalize_then_sum (grpo_trainer.py ~L2790) for offline checks."""
    r = np.asarray(rewards_per_func, float).reshape(-1, G, len(weights))
    z = (r - r.mean(1, keepdims=True)) / (r.std(1, ddof=1, keepdims=True) + eps)
    s = (z.reshape(-1, len(weights)) * np.asarray(weights)).sum(1)
    return (s - s.mean()) / (s.std(ddof=1) + eps)
