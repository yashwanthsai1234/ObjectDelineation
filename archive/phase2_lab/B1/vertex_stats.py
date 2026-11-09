"""Vertices per parsed sample (outer + hole rings) and tokens per vertex, per arm and view, against the raters' targets
of the same eval items (objdel.geometry.parse; prepared records). Reads untrusted data: run with python -I.

    .venv/bin/python -I -m experiments.phase2.lab.B1.vertex_stats [--run-suffix _n150]
"""
import argparse
import json

import numpy as np

from objdel.geometry import parse
from objdel.train import dataset as D

from . import common as C


def verts(text: str) -> int | None:
    p = parse(text, ("nodule",))
    return None if not p.ok else sum(len(o) + sum(len(h) for h in hs) for o, hs in p.objs["nodule"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-suffix", default="")
    a = ap.parse_args()
    recs = {r["item_id"]: r for r in D.records(C.PREPARED, "lidc", "val")}
    out = {}
    for arm in C.ARMS:
        d = C.SCRATCH / "runs" / f"{arm}_s0{a.run_suffix}"
        if not (d / "samples.json").exists():
            continue
        S = json.loads((d / "samples.json").read_text())
        out[arm] = {}
        for view, rows in S.items():
            v, tpv = [], []
            for r in rows:
                for s, n in zip(r["samples"], r["tokens"], strict=True):
                    k = verts(s)
                    if k:
                        v.append(k)
                        tpv.append((n - 1) / k)  # minus the end-of-turn token
            tv = [verts(t) for r in rows for t in recs[r["item_id"]]["targets"] if t]
            out[arm][view] = {"parsed": len(v), "vertices_median": float(np.median(v)),
                              "target_vertices_median": float(np.median(tv)),
                              "tokens_per_vertex_median": round(float(np.median(tpv)), 2)}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
