"""Legacy LOO d_max ceiling per item on val/test, joined with R5 items: all vs found (K-mean IoU_tc >= 0.1)."""
import sys, json
import numpy as np
from objdel.data import store

root, runs = sys.argv[1], sys.argv[2]
for split in ("val", "test"):
    leg = {}
    for it in store.load("lidc", split, root):
        for name, st in it.masks.items():
            st = np.asarray(st, bool); k = len(st)
            d = []
            for r in range(k):
                best = 0.0
                for j in range(k):
                    if j != r:
                        u = np.count_nonzero(st[r] | st[j])
                        best = max(best, np.count_nonzero(st[r] & st[j]) / u if u else 1.0)
                d.append(best)
            leg[it.item_id] = (np.mean(d), d)
    rows = [json.loads(l) for l in open(f"{runs}/eval_{split}/items.jsonl")]
    key = "iou_tc"
    found = [r for r in rows if r[key] >= 0.1]
    def item(rs): return np.mean([leg[r["item_id"]][0] for r in rs])
    def pooled(rs): return np.mean([x for r in rs for x in leg[r["item_id"]][1]])
    print(split, len(rows), "all item %.4f pooled %.4f | found(%d) item %.4f pooled %.4f | model all %.4f found %.4f" % (
        item(rows), pooled(rows), len(found), item(found), pooled(found),
        np.mean([r[key] for r in rows]), np.mean([r[key] for r in found])))
