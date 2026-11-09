"""Idempotence of lossy targets: T1 = build(M); T2 = build(raster(T1)); ... until T_{k+1} == T_k.
Reports: fraction with T2 == T1 (plain idempotence), fixpoint reached within 6 steps, and the IoU /
topology of the fixpoint vs the ORIGINAL mask M.
    python experiments/idempotence.py out.json cfg [cfg...]"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from run_methods import FIX, make_cfg  # noqa: E402
from targets.builder import (LIDC_SCHEMA, RIGA_DERIVED, RIGA_SCHEMA, betti, build_target,  # noqa: E402
                             parse, rasterize_parts)


def one(job):
    cname, ds, key, r = job
    z = np.load(FIX / f"{ds}_experts.npz")
    cfg = make_cfg(cname)
    if ds == "riga":
        M = {"disc": z[f"{key}__disc"][r].astype(bool), "cup": z[f"{key}__cup"][r].astype(bool)}
        schema, derived = RIGA_SCHEMA, RIGA_DERIVED
    else:
        M = {"nodule": z[f"{key}__nodule"][r].astype(bool)}
        schema, derived = LIDC_SCHEMA, ()
    H, W = next(iter(M.values())).shape
    texts = [build_target(M, schema, cfg, derived)[0]]
    for _ in range(6):
        parts = parse(texts[-1], schema)
        ras = {k: rasterize_parts(parts[k], W, H) for k in schema}
        texts.append(build_target(ras, schema, cfg, derived)[0])
        if texts[-1] == texts[-2]:
            break
    conv = texts[-1] == texts[-2]
    parts = parse(texts[-1], schema)
    ras = {k: rasterize_parts(parts[k], W, H) for k in schema}
    ious = {k: float((ras[k] & M[k]).sum() / (ras[k] | M[k]).sum()) for k in schema}
    topo = all(betti(ras[k]) == betti(M[k]) for k in schema) and all(
        betti(fn(ras)) == betti(fn(M)) for _, fn in derived)
    nv0 = sum(len(o) + sum(len(h) for h in hs) for k in schema for o, hs in parse(texts[0], schema)[k])
    nv1 = sum(len(o) + sum(len(h) for h in hs) for k in schema for o, hs in parts[k] for _ in [0])
    return {"cfg": cname, "ds": ds, "idem": texts[1] == texts[0], "converged": conv, "steps": len(texts) - 1,
            "fix_iou_min": min(ious.values()), "fix_topo": topo, "nv0": nv0, "nv_fix": nv1}


if __name__ == "__main__":
    out, cfgs = sys.argv[1], sys.argv[2:]
    jobs = []
    zr = np.load(FIX / "riga_experts.npz")
    items = sorted({k.rsplit("__", 1)[0] for k in zr.files})
    zl = np.load(FIX / "lidc_experts.npz")
    for c in cfgs:
        jobs += [(c, "riga", it, 0) for it in items[:100]]
        jobs += [(c, "lidc", k.rsplit("__", 1)[0], 0) for k in sorted(zl.files)[::2]]
    with Pool(4) as p:
        res = p.map(one, jobs, chunksize=4)
    Path(out).write_text(json.dumps(res))
    for c in cfgs:
        for ds in ("riga", "lidc"):
            rs = [r for r in res if r["cfg"] == c and r["ds"] == ds]
            print(c, ds, len(rs), "idem %.3f conv %.3f fix_topo %.3f fix_iou_min %.4f" % (
                np.mean([r["idem"] for r in rs]), np.mean([r["converged"] for r in rs]),
                np.mean([r["fix_topo"] for r in rs]), min(r["fix_iou_min"] for r in rs)))
