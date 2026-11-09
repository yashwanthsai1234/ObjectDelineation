"""Legacy LOO d_max ceiling on the store (slice masks), both pooling orders, all raters vs voters."""
import sys
import numpy as np
from objdel.data import store
from objdel.reward.cache import voters, signature, min_area

root = sys.argv[1]
for split in ("val", "test"):
    per_item_all, pooled_all, per_item_v, pooled_v, n = [], [], [], [], 0
    for it in store.load("lidc", split, root):
        for name, st in it.masks.items():
            st = np.asarray(st, bool)
            k = len(st)
            if k < 2:
                continue
            n += 1
            ious = np.zeros((k, k))
            for a in range(k):
                for b in range(k):
                    if a != b:
                        u = np.count_nonzero(st[a] | st[b])
                        ious[a, b] = np.count_nonzero(st[a] & st[b]) / u if u else 1.0
            dmax = [max(ious[r, j] for j in range(k) if j != r) for r in range(k)]
            ma = min_area(int(np.count_nonzero(st.any(0))))
            vs = voters(tuple(signature(m, ma) for m in st))
            per_item_all.append(np.mean(dmax)); pooled_all += dmax
            dv = [dmax[r] for r in vs]
            per_item_v.append(np.mean(dv)); pooled_v += dv
    print(split, n, "all: item %.4f pooled %.4f | voters: item %.4f pooled %.4f" % (
        np.mean(per_item_all), np.mean(pooled_all), np.mean(per_item_v), np.mean(pooled_v)))
