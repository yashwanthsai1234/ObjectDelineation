"""Shared helpers for the Phase 1 stack smoke tests (sandbox only, not repo code).

Run every script with `.venv/bin/python -I <script>` from this directory with HF_HOME=./hf.
RIGA files are untrusted third-party data; they are only decoded with PIL here.
"""
import json
import os
import random
import re
import resource
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))  # -I drops the script dir; re-add it explicitly (only our own scripts live here)
os.environ.setdefault("HF_HOME", str(HERE / "hf"))
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

REPO = Path(__file__).resolve().parents[3]
RAW, TRUE_DIR = REPO / "data/raw/riga", REPO / "data/riga_true"
QWEN3VL_LOCAL = Path.home() / ".cache/huggingface/hub/models--Qwen--Qwen3-VL-2B-Instruct/snapshots"

MODELS = {
    "smolvlm-256m": "HuggingFaceTB/SmolVLM-256M-Instruct",
    "smolvlm-500m": "HuggingFaceTB/SmolVLM-500M-Instruct",
    "lfm2.5-vl-450m": "LiquidAI/LFM2.5-VL-450M",
    "qwen3.5-0.8b": "Qwen/Qwen3.5-0.8B",
    "qwen3-vl-2b": "Qwen/Qwen3-VL-2B-Instruct",
}


def model_path(name):
    """Hub id, or the already-cached local snapshot for Qwen3-VL-2B (read-only, avoids a 4 GB re-download)."""
    if name == "qwen3-vl-2b":
        snaps = sorted(QWEN3VL_LOCAL.glob("*"))
        if snaps:
            return str(snaps[0])
    return MODELS[name]


def peak_rss_gb():
    # macOS reports ru_maxrss in bytes
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9


def mps_mem_gb():
    try:
        import torch
        if torch.backends.mps.is_available():
            return torch.mps.driver_allocated_memory() / 1e9
    except Exception:
        pass
    return float("nan")


class Timer:
    def __enter__(self):
        self.t = time.time()
        return self

    def __exit__(self, *a):
        self.s = time.time() - self.t


# ---------------------------------------------------------------- data
TRUE_RE = re.compile(r"^(.+)_image(\d+)-(\d)\.png$")


def _riga_index():
    from experiments.phase0.riga_outlines import jobs as riga_jobs
    primes = {(sub.replace("/", "_"), idx): p for sub, idx, p, _, _ in riga_jobs(str(RAW), None)}
    experts = defaultdict(dict)
    for f in TRUE_DIR.glob("*.png"):
        if m := TRUE_RE.match(f.name):
            experts[(m.group(1), int(m.group(2)))][int(m.group(3))] = f
    keys = sorted(k for k in experts if k in primes and len(experts[k]) >= 2)
    return primes, experts, keys


def riga_samples(n, seed=0, side=448, jitter=True):
    """n fundus crops of side x side px around the optic disc.

    The crop is a square of ~3x the disc size; with jitter the disc centre is placed at a random
    position in the middle 60% of the crop (so 'answer = centre box' does not score well).
    Returns dicts with image (PIL), disc_masks (list of bool arrays side x side, one per expert), key.
    """
    primes, experts, keys = _riga_index()
    rng = random.Random(seed)
    keys = rng.sample(keys, n)
    out = []
    for key in keys:
        ex = [np.asarray(Image.open(f)) for _, f in sorted(experts[key].items())]
        ex = [e[..., 0] if e.ndim == 3 else e for e in ex]
        discs = [e >= 100 for e in ex]
        anyd = np.any(discs, 0)
        ys, xs = np.nonzero(anyd)
        cy, cx = (ys.min() + ys.max()) / 2, (xs.min() + xs.max()) / 2
        size = max(ys.max() - ys.min(), xs.max() - xs.min())
        S = int(3 * size)
        fx, fy = (rng.uniform(0.2, 0.8), rng.uniform(0.2, 0.8)) if jitter else (0.5, 0.5)
        x0, y0 = int(cx - fx * S), int(cy - fy * S)
        img = Image.open(primes[key]).convert("RGB")
        crop = img.crop((x0, y0, x0 + S, y0 + S)).resize((side, side), Image.BICUBIC)
        masks = []
        for d in discs:
            m = Image.fromarray(d.astype(np.uint8) * 255).crop((x0, y0, x0 + S, y0 + S)).resize((side, side), Image.NEAREST)
            masks.append(np.asarray(m) > 127)
        out.append({"image": crop, "disc_masks": masks, "key": f"{key[0]}_image{key[1]}"})
    return out


def mask_bbox(m):
    ys, xs = np.nonzero(m)
    return [int(xs.min()), int(ys.min()), int(xs.max()), int(ys.max())]


def box_iou(a, b):
    ix = max(0, min(a[2], b[2]) - max(a[0], b[0]))
    iy = max(0, min(a[3], b[3]) - max(a[1], b[1]))
    inter = ix * iy
    ua = (a[2] - a[0]) * (a[3] - a[1]) + (b[2] - b[0]) * (b[3] - b[1]) - inter
    return inter / ua if ua > 0 else 0.0


def mask_to_polygon(mask, n=24, scale=1000):
    """Majority mask -> n-vertex polygon (radial sampling from centroid), ints in 0..scale."""
    h, w = mask.shape
    ys, xs = np.nonzero(mask)
    cy, cx = ys.mean(), xs.mean()
    pts = []
    for t in np.linspace(0, 2 * np.pi, n, endpoint=False):
        r = 0.0
        while True:
            y, x = int(round(cy + (r + 1) * np.sin(t))), int(round(cx + (r + 1) * np.cos(t)))
            if not (0 <= y < h and 0 <= x < w) or not mask[y, x]:
                break
            r += 1
        pts.append([int(round((cx + r * np.cos(t)) / (w - 1) * scale)), int(round((cy + r * np.sin(t)) / (h - 1) * scale))])
    return pts


def polygon_json(outer, holes=()):
    return json.dumps({"disc": [{"outer": outer, "holes": list(holes)}]}, separators=(",", ":"))


PROMPT_POLY = ("Outline the optic disc in this retinal fundus image. Reply with JSON only: "
               '{"disc":[{"outer":[[x,y],...],"holes":[]}]} with integer coordinates 0-1000 '
               "relative to the image width and height.")
