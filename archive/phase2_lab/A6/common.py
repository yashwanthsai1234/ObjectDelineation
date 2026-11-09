"""Shared helpers of the A6 performance lab (experiments/phase2/lab/A6): prepared-set loading, timers, profiles.

Run every script as a module with the isolated interpreter (data/ is untrusted):
    .venv/bin/python -I -m experiments.phase2.lab.A6.<script> [--out <scratch dir>]
"""
import cProfile
import json
import os
import platform
import pstats
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from objdel.reward import ExpertCache

PREPARED = Path("data/prepared")


def records(dataset: str, split: str, root: Path = PREPARED) -> list[dict]:
    with open(root / dataset / f"{split}.jsonl") as f:
        return [json.loads(line) for line in f]


def caches(rec: dict, root: Path = PREPARED) -> dict[str, ExpertCache]:
    return {n: ExpertCache.load(root / "caches" / rec["item_id"] / f"{n}.npz", d) for n, d in rec["caches"].items()}


def stats_ms(xs) -> dict:
    a = np.asarray(xs, np.float64)
    return {"n": int(a.size), "median": float(np.median(a)), "p95": float(np.percentile(a, 95)),
            "p99": float(np.percentile(a, 99)), "max": float(a.max()), "mean": float(a.mean()), "sum": float(a.sum())}


def fmt(s: dict) -> str:
    return f"{s['median']:.2f} / {s['p95']:.2f} / {s['max']:.2f}"


class Timer:
    """Accumulates wall time per label (ms) with perf_counter."""

    def __init__(self):
        self.t = defaultdict(float)
        self.n = defaultdict(int)

    def add(self, label: str, dt_s: float):
        self.t[label] += 1e3 * dt_s
        self.n[label] += 1


def profile(fn, *args, top: int = 25, **kw):
    """(result, rows) with rows = top functions by tottime: (ms tottime, ms cumtime, calls, 'file:func')."""
    pr = cProfile.Profile()
    pr.enable()
    out = fn(*args, **kw)
    pr.disable()
    st = pstats.Stats(pr)
    rows = []
    for (file, line, func), (_, nc, tt, ct, _) in st.stats.items():
        short = file.split("site-packages/")[-1].split("src/")[-1]
        rows.append((1e3 * tt, 1e3 * ct, nc, f"{short}:{line}:{func}"))
    rows.sort(reverse=True)
    return out, rows[:top], 1e3 * st.total_tt


def by_module(rows_all) -> dict:
    agg = defaultdict(float)
    for tt, _, _, name in rows_all:
        agg[name.split(":")[0]] += tt
    return dict(sorted(agg.items(), key=lambda kv: -kv[1]))


def env() -> dict:
    return {"python": platform.python_version(), "machine": platform.machine(), "cpu": platform.processor(),
            "loadavg_1m": os.getloadavg()[0], "time": time.strftime("%Y-%m-%d %H:%M:%S")}


def dump(obj, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=1, default=float) + "\n")
