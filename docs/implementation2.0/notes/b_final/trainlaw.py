"""Final review r1 (B): coverage against the train-split shifted law and val-select; offset medians; 3-view mix vs val-select."""
import hashlib, json, sys
import numpy as np
from scipy.stats import ks_2samp
G = json.load(open(sys.argv[1])); c = {k: i for i, k in enumerate(G["cols"])}; rows = G["rows"]
def half(gid): return "shape" if hashlib.sha256(gid.encode()).digest()[0] % 2 == 0 else "select"
SE, OE = (10, 20, 40), (0.15, 0.2, 0.3)
def cell(r): return (int(np.searchsorted(SE, r[c["ext_px"]], side="right")), int(np.searchsorted(OE, r[c["off_euc"]], side="right")))
def shares(rs):
    m = np.zeros((4, 4))
    for r in rs: m[cell(r)] += 1
    return m / m.sum()
def cov(rs, ref):
    t = shares(ref); s = shares(rs); m = t >= 0.01
    return {"n": len(rs), "worst_ratio": round(float((s[m] / t[m]).min()), 3), "tv": round(float(0.5 * np.abs(s - t).sum()), 3),
            "ks": round(float(ks_2samp([r[c["off_euc"]] for r in rs], [r[c["off_euc"]] for r in ref]).statistic), 3)}
val = [r for r in rows if r[0] == "val" and r[1] == "shifted"]
vsh = [r for r in val if half(r[c["group_id"]]) == "shape"]; vse = [r for r in val if half(r[c["group_id"]]) == "select"]
tr = [r for r in rows if r[0] == "train"]
trs = [r for r in tr if r[1] == "shifted"]; troff = [r for r in tr if r[1] != "centred"]; trall = tr
out = {"cols": G["cols"], "kinds": sorted({r[1] for r in tr}),
 "vs_train_law": {"today_offcentre": cov(troff, trs), "today_all_kinds": cov(trall, trs), "val_all": cov(val, trs)},
 "vs_val_select": {"train_law": cov(trs, vse), "val_shape": cov(vsh, vse), "today_offcentre": cov(troff, vse), "today_all_kinds": cov(trall, vse)},
 "vs_val_all": {"train_law": cov(trs, val), "today_all_kinds": cov(trall, val)},
 "offset_median": {"train_law": round(float(np.median([r[c["off_euc"]] for r in trs])), 3),
                   "today_all_kinds": round(float(np.median([r[c["off_euc"]] for r in trall])), 3),
                   "val_all": round(float(np.median([r[c["off_euc"]] for r in val])), 3)},
 "share_off_gt_0.2": {"today_all_kinds": round(float(np.mean([r[c["off_euc"]] > 0.2 for r in trall])), 3),
                      "val_all": round(float(np.mean([r[c["off_euc"]] > 0.2 for r in val])), 3)}}
print(json.dumps(out, indent=1))
