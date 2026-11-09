"""Phase B review: val-shape / val-select patient split (SHA-256(group_id) first byte even -> shape) and baselines."""
import hashlib, json, sys
import numpy as np
from scipy.stats import ks_2samp
G = json.load(open(sys.argv[1])); ITEMS = sys.argv[2]
c = {k: i for i, k in enumerate(G["cols"])}
rows = G["rows"]
iou = {}
for line in open(ITEMS):
    r = json.loads(line); iou[r["item_id"]] = (r["iou_tc"], r["legal"])
def half(gid): return "shape" if hashlib.sha256(gid.encode()).digest()[0] % 2 == 0 else "select"
SE, OE = (10, 20, 40), (0.15, 0.2, 0.3)
def cell(r): return (int(np.searchsorted(SE, r[c["ext_px"]], side="right")), int(np.searchsorted(OE, r[c["off_euc"]], side="right")))
def shares(rs):
    m = np.zeros((4, 4))
    for r in rs: m[cell(r)] += 1
    return m / m.sum()
val = [r for r in rows if r[0] == "val" and r[1] == "shifted"]
out = {}
for h in ("shape", "select", "all"):
    rs = [r for r in val if h == "all" or half(r[c["group_id"]]) == h]
    v = np.array([iou[r[2]][0] for r in rs]); off = np.array([r[c["off_euc"]] for r in rs]); leg = np.array([iou[r[2]][1] for r in rs])
    far, near = off > 0.3, off < 0.15
    out[h] = {"patients": len({r[c["group_id"]] for r in rs}), "items": len(rs), "nodules": len({r[c["series"]] + r[c["nodule"]] for r in rs}),
              "iou": round(v.mean(), 3), "miss": round((v < 0.1).mean(), 3), "found_iou": round(v[v >= 0.1].mean(), 3),
              "near_n": int(near.sum()), "near_iou": round(v[near].mean(), 3), "far_n": int(far.sum()),
              "far_iou": round(v[far].mean(), 3), "far_miss": round((v[far] < 0.1).mean(), 3), "legal": round(leg.mean(), 3)}
tgt = shares([r for r in val if half(r[c["group_id"]]) == "shape"])
toff = [r[c["off_euc"]] for r in val if half(r[c["group_id"]]) == "shape"]
def cov(rs):
    s = shares(rs); m = tgt >= 0.01
    return {"worst_ratio": round(float((s[m] / tgt[m]).min()), 3), "tv": round(float(0.5 * np.abs(s - tgt).sum()), 3),
            "ks": round(float(ks_2samp([r[c["off_euc"]] for r in rs], toff).statistic), 3)}
tr = [r for r in rows if r[0] == "train"]
out["coverage_vs_val_shape"] = {"train_shifted": cov([r for r in tr if r[1] == "shifted"]),
                                "train_shifted_plus_rescaled_unweighted": cov([r for r in tr if r[1] != "centred"])}
# item 27: decomposition of val IoU gain from the far-miss gate (val-select)
print(json.dumps(out, indent=1))
