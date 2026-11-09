"""B1 pixel/duplicate audit. Run from REPO: .venv/bin/python -I b1_pixels.py OUTDIR
Reads data/store with pixels. Per split and window kind (same draws as build_eval_views seed 0): intensity statistics of
the model image, nodule contrast, HU inside the nodule; sha check of built views against a fresh render; exact and
near duplicate check of CT slices / crops across splits."""
import json, sys, hashlib
from pathlib import Path
import numpy as np, pyarrow.parquet as pq
from dataclasses import replace
from scipy import ndimage as ndi
from objdel.data.store import Item, decode_masks, _from_npy
from objdel.targets import augment, frame, FRAMES
from objdel.targets.frame import union_box
sys.path.insert(0, str(Path(__file__).parent))
from b1_geometry import windows
out = Path(sys.argv[1])
cfg0 = replace(FRAMES["lidc"], image_side=None)

def load(split):
    pf = pq.ParquetFile(f"data/store/lidc/{split}.parquet")
    for b in pf.iter_batches(batch_size=16):
        for r in b.to_pylist():
            yield Item(r["item_id"], r["dataset"], r["group_id"], r["split"], _from_npy(r["image"]),
                       {m["name"]: decode_masks(m["bits"], r["width"]) for m in r["masks"]}, r["raters"], json.loads(r["meta"]))

recs = {}
for root, split, kinds in (("data/prepared/lidc", "train", ["centred"]), ("data/prepared/lidc", "val", ["centred"]),
                           ("data/views/lidc-train-shifted/lidc", "train", ["shifted"]), ("data/views/lidc-train-rescaled/lidc", "train", ["rescaled"]),
                           ("data/eval_views/lidc-val-shifted/lidc", "val", ["shifted"]), ("data/eval_views/lidc-val-rescaled/lidc", "val", ["rescaled"])):
    for l in open(f"{root}/{split}.jsonl"):
        r = json.loads(l); recs[(split, kinds[0], r["item_id"])] = r
rows, dup = [], {}
shaok = {}
for split in ("train", "val", "test"):
    for it in load(split):
        img = it.image
        hs = hashlib.sha256(img.tobytes()).hexdigest()
        ub = union_box(it)
        cx, cy = (ub[0]+ub[2])//2, (ub[1]+ub[3])//2
        crop = img[max(cy-16,0):cy+16, max(cx-16,0):cx+16]
        dup[it.item_id] = (split, it.group_id, hs, hashlib.sha256(crop.astype(np.int16).tobytes()).hexdigest(),
                           np.pad(crop, ((0, 32-crop.shape[0]), (0, 32-crop.shape[1])), constant_values=-1000).astype(np.float32))
        union = np.logical_or.reduce(it.masks["nodule"], axis=0)
        hu_in = img[union]
        for kind, (dx, dy, S, fb, kept) in windows(replace(it, image=None)).items():
            f = frame(it, cfg0, offset=(dx, dy), side=S) if kind != "centred" else frame(it, cfg0)
            g = f.image[..., 0].astype(np.float32)
            m = np.logical_or.reduce(f.masks["nodule"], axis=0)
            ring = ndi.binary_dilation(m, iterations=4) & ~ndi.binary_dilation(m, iterations=1)
            rows.append((split, kind, S, float(g.mean()), float(g.std()), float((g == 0).mean()), float((g == 255).mean()),
                         float(g[m].mean() - g[ring].mean()) if ring.any() else np.nan, float(np.median(hu_in)), float(np.percentile(hu_in, 10))))
            key = (split, kind, it.item_id)
            if key in recs and (len(shaok.setdefault((split, kind), [])) < 150):
                f2 = frame(it, offset=(dx, dy), side=S) if kind != "centred" else frame(it)
                from objdel.train.dataset import _pixels_sha
                shaok[(split, kind)].append((_pixels_sha(f2.image) == recs[key]["image_sha256"], list(f2.box) == recs[key]["box"]))
    print(split, flush=True)
json.dump({"cols": ["split","kind","S","mean","std","black","white","contrast","hu_med","hu_p10"], "rows": rows}, open(out/"pix.json", "w"))
json.dump({f"{a}/{b}": [sum(x for x, _ in v), sum(y for _, y in v), len(v)] for (a, b), v in shaok.items()}, open(out/"shacheck.json", "w"))
# duplicates
ids = list(dup)
sp = np.array([dup[i][0] for i in ids]); gp = [dup[i][1] for i in ids]
from collections import defaultdict
for idx, name in ((2, "slice_sha"), (3, "crop_sha")):
    d = defaultdict(set)
    for i in ids: d[dup[i][idx]].add((dup[i][0], dup[i][1]))
    print(name, "hashes", len(d), "shared across >1 patient", sum(len({p for _, p in v}) > 1 for v in d.values()), "shared across >1 split", sum(len({s for s, _ in v}) > 1 for v in d.values()))
X = np.stack([dup[i][4].reshape(-1) for i in ids]); X = (X - X.mean(1, keepdims=True)); X /= (np.linalg.norm(X, axis=1, keepdims=True) + 1e-6)
res = {}
for a, b in (("test", "train"), ("val", "train"), ("test", "val")):
    ia, ib = np.where(sp == a)[0], np.where(sp == b)[0]
    C = X[ia] @ X[ib].T
    best = C.max(1)
    res[f"{a}->{b}"] = {"max_cos_median": float(np.median(best)), "p99": float(np.percentile(best, 99)), "n>=0.95": int((best >= 0.95).sum()), "n>=0.98": int((best >= 0.98).sum()), "n>=0.995": int((best >= 0.995).sum())}
# within-train baseline (other patient only)
ia = np.where(sp == "train")[0]; C = X[ia] @ X[ia].T
g = np.array([gp[i] for i in ia]); C[g[:, None] == g[None, :]] = -1
res["train->train_other_patient"] = {"max_cos_median": float(np.median(C.max(1))), "n>=0.95": int((C.max(1) >= 0.95).sum()), "n>=0.98": int((C.max(1) >= 0.98).sum())}
json.dump(res, open(out/"neardup.json", "w"), indent=1); print(json.dumps(res, indent=1))
