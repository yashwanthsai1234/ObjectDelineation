"""Build targets for every fixture mask (each rater separately) with every method config.

    python experiments/run_methods.py <out.jsonl> [crop=native|pad_div|pad_sq] [configs...]
Writes one JSON line per (config, item, rater) with per-object metrics and the target text.
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from targets.builder import (LIDC_SCHEMA, RIGA_DERIVED, RIGA_SCHEMA, BuildConfig,  # noqa: E402
                             build_target)
from targets.crops import apply_crop  # noqa: E402

FIX = HERE.parent.parent / "fixtures"

CONFIGS = {
    "exact_rings": ("exact_rings", {}),
    "dp2": ("dp", {"eps": 2}), "dp3": ("dp", {"eps": 3}), "dp5": ("dp", {"eps": 5}), "dp8": ("dp", {"eps": 8}),
    "vw10": ("vw", {"area": 10}), "vw25": ("vw", {"area": 25}), "vw50": ("vw", {"area": 50}),
    "vw100": ("vw", {"area": 100}),
    "potrace": ("minseg_linf", {"t_px": 0.5}),
    "linf2": ("minseg_linf", {"tau": 2}), "linf3": ("minseg_linf", {"tau": 3}), "linf4": ("minseg_linf", {"tau": 4}),
    "exact": ("minseg_exact", {}),
    "hyb2": ("minseg_hybrid", {"tau": 2}), "hyb3": ("minseg_hybrid", {"tau": 3}),
    "hyb4": ("minseg_hybrid", {"tau": 4}), "hyb5": ("minseg_hybrid", {"tau": 5}),
}
# variants: <base>@<iou_min>@<repair>, e.g. hyb3@0.99@auto



def make_cfg(cname):
    base, *rest = cname.split("@")
    method, params = CONFIGS[base]
    iou_min = float(rest[0]) if rest else (1.0 if base == "exact" else 0.98)
    repair = rest[1] if len(rest) > 1 else "best"
    return BuildConfig(method, params, iou_min=iou_min, repair=repair)


def jobs(cfgnames, crop):
    out = []
    for ds in __import__("os").environ.get("DS", "riga,lidc").split(","):
        z = np.load(FIX / f"{ds}_experts.npz")
        if ds == "riga":
            items = sorted({k.rsplit("__", 1)[0] for k in z.files})
            for it in items:
                for r in range(z[f"{it}__disc"].shape[0]):
                    for c in cfgnames:
                        out.append((c, ds, it, r, crop))
        else:
            for k in sorted(z.files):
                for r in range(z[k].shape[0]):
                    for c in cfgnames:
                        out.append((c, ds, k.rsplit("__", 1)[0], r, crop))
    return out


_Z = {}


def run(job):
    cname, ds, it, r, crop = job
    if ds not in _Z:
        _Z[ds] = np.load(FIX / f"{ds}_experts.npz")
    z = _Z[ds]
    cfg = make_cfg(cname)
    if ds == "riga":
        masks = {"disc": z[f"{it}__disc"][r].astype(bool), "cup": z[f"{it}__cup"][r].astype(bool)}
        schema, derived = RIGA_SCHEMA, RIGA_DERIVED
    else:
        masks = {"nodule": z[f"{it}__nodule"][r].astype(bool)}
        schema, derived = LIDC_SCHEMA, ()
    masks, (H, W) = apply_crop(masks, crop)
    text, objs, info = build_target(masks, schema, cfg, derived)
    per = {}
    for k, o in objs.items():
        per[k] = {"status": o.status, "iou": o.iou, "b_src": o.betti_src, "b_out": o.betti_out,
                  "topo": o.betti_src == o.betti_out, "legal": o.legal, "vertices": o.vertices,
                  "parts": len(o.parts), "holes": sum(len(h) for _, h in o.parts),
                  "text": json.dumps({k: json.loads(text)[k]}, separators=(",", ":"))}
    return {"cfg": cname, "ds": ds, "item": it, "rater": r, "crop": crop, "W": W, "H": H, "text": text,
            "objs": per, "derived_ok": info["derived_ok"], "time": info["time"], "rounds": info["rounds"],
            "raw": info["raw"], "bevels": info["bevels"], "schema": info["schema"]}


if __name__ == "__main__":
    out = Path(sys.argv[1])
    crop = sys.argv[2] if len(sys.argv) > 2 else "native"
    cfgs = sys.argv[3:] or list(CONFIGS)
    js = jobs(cfgs, crop)
    with Pool(int(__import__("os").environ.get("NPROC", 6))) as pool, out.open("w") as f:
        for i, rec in enumerate(pool.imap(run, js, chunksize=8)):
            f.write(json.dumps(rec) + "\n")
            if i % 2000 == 0:
                print(i, "/", len(js), flush=True)
