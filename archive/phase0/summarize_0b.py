"""Phase 0b numbers: how often outputs have good IoU but a topology no annotator gave.

    .venv/bin/python -I -m experiments.phase0.summarize_0b experiments/phase0/results
"""
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ERRORS = ("extra_objects", "missing_objects", "extra_holes", "missing_holes")


def load(path):
    return list(csv.DictReader(open(path))) if Path(path).exists() else []


def block(rows, iou_cuts=(0.5, 0.8)):
    if not rows:
        return None
    iou = np.array([float(r["iou"]) for r in rows])
    wrong = np.array([float(r["topo_freq"]) == 0 for r in rows])
    out = {"n": len(rows), "iou_median": round(float(np.median(iou)), 3),
           "topology_wrong": {"n": int(wrong.sum()), "pct": round(100 * wrong.mean(), 2)}}
    for c in iou_cuts:
        sel = iou >= c
        out[f"iou_ge_{c}"] = {"n": int(sel.sum()), "topology_wrong": int((sel & wrong).sum()),
                              "pct_of_all": round(100 * (sel & wrong).mean(), 2),
                              "pct_of_these": round(100 * (sel & wrong).sum() / max(1, sel.sum()), 2)}
    if rows and "reader_iou" in rows[0]:
        ref = np.array([float(r["reader_iou"]) for r in rows])
        sel = iou >= ref
        out["reader_iou_median"] = round(float(np.median(ref)), 3)
        out["iou_ge_reader_iou"] = {"n": int(sel.sum()), "topology_wrong": int((sel & wrong).sum()),
                                    "pct_of_all": round(100 * (sel & wrong).mean(), 2),
                                    "pct_of_these": round(100 * (sel & wrong).sum() / max(1, sel.sum()), 2)}
    if wrong.any():
        out["iou_of_topology_wrong_median"] = round(float(np.median(iou[wrong])), 3)
        out["error_types"] = dict(Counter(e for r, w in zip(rows, wrong) if w for e in ERRORS if int(r[e]) > 0))
    band = np.array([float(r["band"]) for r in rows])
    out["band_median"] = round(float(np.median(band)), 3)
    out["band_below_0.9_pct"] = round(100 * float((band < 0.9).mean()), 2)
    return out


def main(root):
    root = Path(root)
    s = {}
    riga = load(root / "medsam_riga.csv")
    for st in ("disc", "cup", "rim"):
        s[f"medsam_riga_{st}"] = block([r for r in riga if r["structure"] == st])
    lidc = [r for p in sorted(root.glob("medsam_lidc*.csv")) for r in load(p)]
    s["medsam_lidc_all_slices"] = block(lidc)
    s["medsam_lidc_readers_agree"] = block([r for r in lidc if r["experts_agree"] == "True"])
    s["medsam_lidc_readers_disagree"] = block([r for r in lidc if r["experts_agree"] == "False"])
    for mode in ("box", "polygon"):
        v = load(root / f"vlm_{mode}" / "scores.csv")
        if not v:
            continue
        images = {(r["subset"], r["image"]) for r in v}
        s[f"vlm_{mode}"] = {
            "images": len(images),
            "parsed_pct": round(100 * np.mean([r["parsed"] == "True" for r in v if r["structure"] == "disc"]), 1),
            "valid_disc": dict(Counter(r["valid"] for r in v if r["structure"] == "disc")),
            "valid_cup": dict(Counter(r["valid"] for r in v if r["structure"] == "cup")),
            **{st: block([r for r in v if r["structure"] == st]) for st in ("disc", "cup", "rim")},
        }
    (root / "summary_0b.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
