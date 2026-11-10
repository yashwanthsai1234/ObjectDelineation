"""POST-HOC (not pre-registered) re-scoring of zeroshot.py's generic-prompt answers: Qwen3.5 answered the generic
{"bbox": ...} 0-1000 prompt with its native "bbox_2d" key, and LFM2.5-VL with [0, 1] floats. Lenient parse: either
key, and a box whose 4 values are all <= 1 is read as [0, 1]. Tells whether the generic-prompt failures are
convention-following failures or localisation failures. Reads untrusted data: run with python -I.

    .venv/bin/python -I -m experiments.phase2.lab.B1.zs_posthoc
"""
import json
import re

import numpy as np

from . import common as C
from .zeroshot import NUM, QUERY, parse_box, prompts


def lenient(text: str):
    m = re.search(rf'"bbox(?:_2d)?"\s*:\s*\[\s*({NUM})\s*,\s*({NUM})\s*,\s*({NUM})\s*,\s*({NUM})', text)
    if not m:
        return None
    b = [float(x) for x in m.groups()]
    b = b if max(b) <= 1 else [v / 1000 for v in b]
    return b if all(0 <= v <= 1 for v in b) and b[0] < b[2] and b[1] < b[3] else None


def main():
    out = {}
    for ds, n in (("riga", 64), ("lidc", 64)):
        V = C.views(ds, C.eval_records(ds, n))
        struct = QUERY[ds][2]
        for arm in C.ARMS:
            raw = json.loads((C.SCRATCH / "zeroshot" / f"{arm}_raw.json").read_text())
            for view, items in V.items():
                gt = [C.box_of(x["caches"][struct].union, x["W"]) for x in items]
                boxes = [lenient(s) for s in raw[f"{ds}/{view}/generic"]]
                iou = np.array([C.box_iou(b, g) if b else 0.0 for b, g in zip(boxes, gt, strict=True)])
                out[f"{arm}/{ds}/{view}/generic_lenient"] = {
                    "parse": round(float(np.mean([b is not None for b in boxes])), 3),
                    "box_iou_median": C.ci(lambda i, v=iou: np.median(v[i]), len(items))}
        for view, items in V.items():  # paired Qwen3.5 - LFM2.5-VL, native prompts (strict parse, as zeroshot.py)
            gt = [C.box_of(x["caches"][struct].union, x["W"]) for x in items]
            ious = {}
            for arm in ("qwen3.5-0.8b", "lfm2.5-vl-450m"):
                raw = json.loads((C.SCRATCH / "zeroshot" / f"{arm}_raw.json").read_text())
                _, key, scale = prompts(arm, ds)["native"]
                bx = [parse_box(t, key, scale) for t in raw[f"{ds}/{view}/native"]]
                ious[arm] = np.array([C.box_iou(b, g) if b else 0.0 for b, g in zip(bx, gt, strict=True)])
            d = ious["qwen3.5-0.8b"] - ious["lfm2.5-vl-450m"]
            out[f"{ds}/{view}/native: qwen - lfm"] = {"mean": C.ci(lambda i, d=d: d[i].mean(), len(d)),
                                                      "median_diff": C.ci(lambda i: np.median(ious["qwen3.5-0.8b"][i])
                                                                          - np.median(ious["lfm2.5-vl-450m"][i]), len(d))}
    print(json.dumps(out, indent=1))
    (C.SCRATCH / "zeroshot" / "posthoc_generic.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
