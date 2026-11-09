"""Final review r1 (B): R5 val-select percentile CI width by item, nodule and patient (B = 4,000, seed 0)."""
import hashlib, json, sys
import numpy as np
G = json.load(open(sys.argv[1])); c = {k: i for i, k in enumerate(G["cols"])}
iou = {json.loads(l)["item_id"]: json.loads(l)["iou_tc"] for l in open(sys.argv[2])}
def half(g): return "shape" if hashlib.sha256(g.encode()).digest()[0] % 2 == 0 else "select"
rs = [r for r in G["rows"] if r[0] == "val" and r[1] == "shifted" and half(r[c["group_id"]]) == "select"]
v = np.array([iou[r[2]] for r in rs]); rng = np.random.default_rng(0)
out = {"items": len(rs), "mean": round(float(v.mean()), 4)}
for name, key in (("item", lambda r: r[2]), ("nodule", lambda r: r[c["series"]] + r[c["nodule"]]), ("patient", lambda r: r[c["group_id"]])):
    ks = sorted({key(r) for r in rs}); idx = {k: i for i, k in enumerate(ks)}; g = np.array([idx[key(r)] for r in rs])
    s = np.bincount(g, v, len(ks)); n = np.bincount(g, None, len(ks)); m = []
    for _ in range(4000):
        w = np.bincount(rng.integers(0, len(ks), len(ks)), None, len(ks)); m.append((w * s).sum() / (w * n).sum())
    lo, hi = np.percentile(m, [2.5, 97.5]); out[name] = {"clusters": len(ks), "ci": [round(lo, 3), round(hi, 3)], "width": round(hi - lo, 3)}
print(json.dumps(out))
