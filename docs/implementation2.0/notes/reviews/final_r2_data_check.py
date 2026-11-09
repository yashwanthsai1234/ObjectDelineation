import json, sys, collections
import numpy as np
from scipy.stats import ks_2samp
G = json.load(open(sys.argv[1])); c = {k: i for i, k in enumerate(G["cols"])}; rows = G["rows"]
SE, OE = (10, 20, 40), (0.15, 0.2, 0.3)
S0 = {r[c["item_id"]]: r[c["S"]] for r in rows if r[1] == "centred" and r[0] == "train"}
def cell(r, apparent):
    e = r[c["ext_px"]]
    if apparent and r[1] == "rescaled": e = e * S0[r[c["item_id"]]] / r[c["S"]]
    return (int(np.searchsorted(SE, e, side="right")), int(np.searchsorted(OE, r[c["off_euc"]], side="right")))
def shares(rs, ap):
    m = np.zeros((4, 4))
    for r in rs: m[cell(r, ap)] += 1
    return m / m.sum()
def cov(rs, ref, ap=True):
    t = shares(ref, False); s = shares(rs, ap); m = t >= 0.01
    return round(float((s[m] / t[m]).min()), 3), round(float(0.5 * np.abs(s - t).sum()), 3), round(float(ks_2samp([r[c["off_euc"]] for r in rs], [r[c["off_euc"]] for r in ref]).statistic), 3)
tr = [r for r in rows if r[0] == "train"]; trs = [r for r in tr if r[1] == "shifted"]; troff = [r for r in tr if r[1] != "centred"]
print("today offcentre raw", cov(troff, trs, False), "apparent", cov(troff, trs, True))
print("today all raw", cov(tr, trs, False), "apparent", cov(tr, trs, True))
resc = [r for r in tr if r[1] == "rescaled"]
print("rescaled cells raw vs apparent changed share", np.mean([cell(r, False) != cell(r, True) for r in resc]))
# null: random patient subsamples of train shifted vs full train law
bypat = collections.defaultdict(list)
for r in trs: bypat[r[c["group_id"]]].append(r)
pats = sorted(bypat); rng = np.random.default_rng(0)
for npat, obs in ((78, "val_all 0.376/0.077/0.061"), (38, "val_select vs train law ?")):
    W, T, K = [], [], []
    for _ in range(500):
        sub = [r for p in rng.choice(pats, npat, replace=False) for r in bypat[p]]
        w, t, k = cov(sub, trs, False); W.append(w); T.append(t); K.append(k)
    print(npat, "patients: worst p5/p50", np.percentile(W, [5, 50]).round(3), "TV p50/p95", np.percentile(T, [50, 95]).round(3), "KS p50/p95", np.percentile(K, [50, 95]).round(3), "obs", obs)
val = [r for r in rows if r[0] == "val" and r[1] == "shifted"]
print("val all vs train law", cov(val, trs, False))
print("---- size-bin feasibility for the rescaled kind (apparent size)")
t = shares(trs, False)
print("target size-bin shares", t.sum(1).round(3), "offset-bin shares", t.sum(0).round(3))
ap = shares(resc, True); print("rescaled draw0 apparent size-bin shares", ap.sum(1).round(3))
# items that can EVER give apparent size bin k under either neighbour side: S0/S in {0.5 (S0=50->100), 0.8 (100->125), 2 (100->50), 1.25(125->100)}
from collections import Counter
nb = {50: [100], 100: [50, 125], 125: [100, 200]}
cnt = Counter()
for r in tr:
    if r[1] != "centred": continue
    s0 = r[c["S"]]; e = r[c["ext_px"]]
    bins = {int(np.searchsorted(SE, e * s0 / s, side="right")) for s in nb.get(s0, [])}
    for b in bins: cnt[b] += 1
print("items able to reach each apparent size bin (any neighbour):", dict(sorted(cnt.items())), "of", sum(1 for r in tr if r[1]=='centred'))
