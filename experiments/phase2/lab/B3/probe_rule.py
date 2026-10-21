"""Exploratory (NOT pre-registered): does a paired criterion against the centroid disc make the precision probe's
geometry part discriminate a location-and-size-only blob? Same probe set as precision_probe (seed 2026, calibrated
shape, N items). Per item in the 6-16 px bins: IoU of each synthetic model minus the centroid disc's IoU (illegal
answer = IoU 0). Reports the median paired gain, the win share, and A4's rule with tighter geometry thresholds.

    .venv/bin/python -I -m experiments.phase2.lab.B3.probe_rule --census <scratch>/census.jsonl \
        --calib <lab>/precision_probe.json --out <lab>/probe_rule.json
"""
import argparse
import json

import numpy as np

from . import precision_probe as PP
from . import probe_lib as P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True)
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=1200)
    a = ap.parse_args()
    with open(a.calib) as fh:
        shape = tuple(json.load(fh)["calibration"]["shape_star"])
    diams, prior = PP.train_diams(a.census)
    items = PP.draw_items(a.n, shape, diams, np.random.default_rng(2026))
    sel = [it for it in items if 6 <= it["diam"] < 16]
    disc = np.array([P.score(it, P.model_centroid_disc(it, None)).get("iou", 0.0) for it in sel])
    models = {"exact": P.model_exact, "sigma_0.25": P.model_noise(0.25), "sigma_0.5": P.model_noise(0.5),
              "sigma_0.75": P.model_noise(0.75), "sigma_1": P.model_noise(1.0), "shift_1px": P.model_shift,
              "centroid_disc": P.model_centroid_disc}
    out = {"items_6_16": len(sel), "disc_iou_median": float(np.median(disc)), "models": {}}
    for name, m in models.items():
        rng = np.random.default_rng(7)
        rows = [P.score(it, m(it, rng)) for it in sel]
        iou = np.array([r.get("iou", 0.0) if r["legal"] else 0.0 for r in rows])
        lg = [r for r in rows if r["legal"]]
        gain = iou - disc
        out["models"][name] = {
            "legal": float(np.mean([r["legal"] for r in rows])),
            "iou_median": float(np.median([r["iou"] for r in lg])),
            "assd_median": float(np.median([r["assd"] for r in lg])),
            "paired_gain_median": float(np.median(gain)), "win_share": float(np.mean(gain > 0)),
            "tight_rule_iou085_assd045": bool(np.median([r["iou"] for r in lg]) >= 0.85
                                              and np.median([r["assd"] for r in lg]) <= 0.45)}
        print(name, out["models"][name], flush=True)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
