"""B1 geometry audit. Run from REPO: .venv/bin/python -I b1_geometry.py OUTDIR
Reads data/store (masks+meta only, no pixels). Simulates the exact windows build_eval_views draws (seed 0) for every
split, and the train-row mix of rl_shifted.yaml. Writes OUTDIR/geom.json and geom_rows.npz (per-item arrays)."""
import json, sys, hashlib
from pathlib import Path
import numpy as np, pyarrow.parquet as pq
from dataclasses import replace
from objdel.data.store import Item, decode_masks
from objdel.targets import augment, frame
from objdel.targets.frame import union_box, policy_side, FRAMES, centred

def load(split):
    pf = pq.ParquetFile(f"data/store/lidc/{split}.parquet")
    for b in pf.iter_batches(batch_size=64, columns=["item_id","dataset","group_id","split","raters","width","masks","meta"]):
        for r in b.to_pylist():
            yield Item(r["item_id"], r["dataset"], r["group_id"], r["split"], None,
                       {m["name"]: decode_masks(m["bits"], r["width"]) for m in r["masks"]}, r["raters"], json.loads(r["meta"]))

def windows(it, seed=0):
    """-> dict kind -> (dx,dy,S,fallback,kept)"""
    S0 = frame(it).W
    out = {"centred": (0, 0, S0, False, False)}
    rng = augment.row_rng(seed, "shifted", it.item_id)
    got = augment.draw(it, rng, side=S0)
    (dx, dy), S = got or ((0, 0), S0)
    out["shifted"] = (dx, dy, S, got is None, False)
    rng = augment.row_rng(seed, "rescaled", it.item_id)
    got = augment.neighbour_side(it, S0, rng)
    kept = got is None
    got = got or augment.draw(it, rng, side=S0)
    (dx, dy), S = got or ((0, 0), S0)
    out["rescaled"] = (dx, dy, S, got is None, kept)
    return out

def main(out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    res = {}
    for split in ("train", "val", "test"):
        rows = []
        for it in load(split):
            ub = union_box(it)
            k = len(it.raters)
            masks = it.masks["nodule"]
            areas = masks.reshape(k, -1).sum(1)
            m = it.meta
            ser = m["series"]; nod = it.item_id.split("/")[2]
            for kind, (dx, dy, S, fb, kept) in windows(it).items():
                f = frame(it, offset=(dx, dy), side=S) if kind != "centred" else frame(it)
                inside = all(int(f.masks["nodule"][r].sum()) == int(areas[r]) for r in range(k))
                x0, y0 = f.box[:2]
                cx, cy = (ub[0]+ub[2])/2, (ub[1]+ub[3])/2
                wx, wy = x0 + f.W*f.scale/2, y0 + f.H*f.scale/2
                ext = max(ub[2]-ub[0], ub[3]-ub[1])
                pad = max(0, -f.box[0]) + max(0, -f.box[1]) + max(0, f.box[2]-512) + max(0, f.box[3]-512)
                other = [(a,b,c,d) for a,b,c,d in m.get("other_nodules", ())]
                crowd = any(a < f.box[2] and f.box[0] < c and b < f.box[3] and f.box[1] < d for a,b,c,d in other)
                rows.append((split, kind, it.item_id, it.group_id, ser, nod, k, f.W, f.scale,
                             float(np.hypot(cx-wx, cy-wy))/f.W, float(max(abs(cx-wx), abs(cy-wy)))/f.W, ext/f.W,
                             ext, int(areas.min()), int(areas.max()), int(inside), int(pad > 0), int(crowd), int(fb), int(kept),
                             int(f.box[0]), int(f.box[1]), int(dx), int(dy), len(other)))
        res[split] = rows
        print(split, len(rows), flush=True)
    cols = ["split","kind","item_id","group_id","series","nodule","k","S","scale","off_euc","off_cheb","ext_frac","ext_px",
            "area_min","area_max","masks_inside","padded","crowded","fallback","kept","x0","y0","dx","dy","n_other"]
    flat = [r for s in res.values() for r in s]
    json.dump({"cols": cols, "rows": flat}, open(out/"geom.json", "w"))
if __name__ == "__main__":
    main(sys.argv[1])
