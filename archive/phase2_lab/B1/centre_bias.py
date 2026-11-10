"""POST-HOC diagnostic: do trained models localise the nodule or draw at the frame centre? For every parsed sample,
the vertex centroid (frame px) is compared with the frame centre and with the centroid of the raters' union; on the
shifted view the two differ by the shift. Reads untrusted data: run with python -I.

    .venv/bin/python -I -m experiments.phase2.lab.B1.centre_bias [--run-suffix _n150]
"""
import argparse
import json

import numpy as np

from objdel.geometry import parse

from . import common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-suffix", default="")
    a = ap.parse_args()
    V = C.views("lidc", C.eval_records("lidc", 64))
    out = {}
    for arm in C.ARMS:
        d = C.SCRATCH / "runs" / f"{arm}_s0{a.run_suffix}"
        if not (d / "samples.json").exists():
            continue
        S = json.loads((d / "samples.json").read_text())
        out[arm] = {}
        for view, rows in S.items():
            to_centre, to_object, closer = [], [], []
            for x, r in zip(V[view], rows, strict=True):
                ys, xs = np.nonzero(x["caches"]["nodule"].union)
                obj = np.array([xs.mean() + 0.5, ys.mean() + 0.5])
                centre = np.array([x["W"] / 2, x["H"] / 2])
                for s in r["samples"]:
                    p = parse(s, ("nodule",))
                    if not p.ok or not p.objs["nodule"]:
                        continue
                    v = np.concatenate([o for o, _ in p.objs["nodule"]]) * x["W"] / 1000
                    c = v.mean(0)
                    to_centre.append(np.linalg.norm(c - centre))
                    to_object.append(np.linalg.norm(c - obj))
                    closer.append(to_centre[-1] < to_object[-1])
            out[arm][view] = {"parsed": len(closer), "median_px_to_frame_centre": round(float(np.median(to_centre)), 2),
                              "median_px_to_object_centroid": round(float(np.median(to_object)), 2),
                              "share_closer_to_frame_centre": round(float(np.mean(closer)), 3)}
    shift = [np.hypot(*x["shift"]) for x in V["shifted"]]
    out["shift_px_median"] = round(float(np.median(shift)), 2)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
