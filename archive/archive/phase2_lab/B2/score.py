"""B2 scoring: re-score stored K = 8 samples with the frozen Phase 1 reward (objdel.eval.metrics: score_objects,
legality gate, advantages_v4, REWARD_V4 floors) on each item's own view, and keep per-item arrays for the GRPO-signal
metrics of prereg.py.

    .venv/bin/python -I -m experiments.phase2.lab.B2.score [--runs <name> ...]

Inputs: <scratch>/phase2/B1/runs/<run>/samples.json (B1) and <scratch>/phase2/B2/runs/<run>/samples.json (interp.py).
Output: <scratch>/phase2/B2/scored/<run>.json, per view: item_id, R (K floats), legal, iou_tc, advantage, sigma, dup,
distinct (distinct texts / K). Reads untrusted data (the val store): run with python -I.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from objdel.eval import metrics

from ..B1 import common as C
from .prereg import DEPTH_RUNS, SEED_RUNS, WEIGHTS

B2 = C.SCRATCH.parent / "B2"


def run_dir(run: str) -> Path:
    for root in (C.SCRATCH / "runs", B2 / "runs"):
        if (root / run / "samples.json").exists():
            return root / run
    raise FileNotFoundError(run)


def score_view(gen: list[dict], items: list[dict]) -> dict:
    out = {k: [] for k in ("item_id", "R", "legal", "iou_tc", "advantage", "sigma", "dup", "distinct", "format")}
    for g, x in zip(gen, items, strict=True):
        assert g["item_id"] == x["item_id"]
        m = metrics(
            g["samples"],
            {n: c.raters for n, c in x["caches"].items()},
            x["W"],
            x["H"],
            x["names"],
            g["truncated"],
            x["caches"],
        )
        s = m["samples"]
        R = sum(w * np.asarray(s[t]) for w, t in zip(WEIGHTS, ("iou_tc", "f5", "legal", "g_area"), strict=True))
        out["item_id"].append(x["item_id"])
        out["R"].append([float(v) for v in R])
        out["legal"].append(s["legal"])
        out["format"].append(s["format_ok"])
        out["iou_tc"].append(s["iou_tc"])
        out["advantage"].append(s["advantage"])
        out["sigma"].append(float(m["groups"]["signal"][0]))
        out["dup"].append(float(m["groups"]["dup"][0]))
        out["distinct"].append(len(set(g["samples"])) / len(g["samples"]))
    return out


def main():
    default = sorted({r for p in (*DEPTH_RUNS.values(), *SEED_RUNS.values()) for r in p})
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="*", default=default)
    a = ap.parse_args()
    (B2 / "scored").mkdir(parents=True, exist_ok=True)
    V = None
    for run in a.runs:
        samples = json.loads((run_dir(run) / "samples.json").read_text())
        n = len(samples["centred"])
        if V is None or len(V["centred"]) != n:
            V = C.views("lidc", C.eval_records("lidc", n))  # B1's items and shifts (seed 0)
        res = {view: score_view(samples[view], V[view]) for view in ("centred", "shifted")}
        (B2 / "scored" / f"{run}.json").write_text(json.dumps(res))
        print(run, {v: round(float(np.mean(res[v]["legal"])), 3) for v in res}, flush=True)


if __name__ == "__main__":
    main()
