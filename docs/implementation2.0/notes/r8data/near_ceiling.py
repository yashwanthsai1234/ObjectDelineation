"""Legacy LOO d_max ceiling and R5 IoU per offset subset on val (offset = union-box centre vs window centre / side)."""
import sys, json
import numpy as np
from objdel.data import store
from objdel.targets.frame import union_box

repo = sys.argv[1]
view = {json.loads(l)["item_id"]: json.loads(l) for l in open(f"{repo}/data/eval_views/lidc-val-shifted/lidc/val.jsonl")}
r5 = {json.loads(l)["item_id"]: json.loads(l) for l in open(f"{repo}/runs/r5/eval/eval_val/items.jsonl")}
rows = []
for it in store.load("lidc", "val", f"{repo}/data/store"):
    st = np.asarray(it.masks["nodule"], bool); k = len(st)
    d = [max(np.count_nonzero(st[r] & st[j]) / max(np.count_nonzero(st[r] | st[j]), 1) for j in range(k) if j != r) for r in range(k)]
    x0, y0, x1, y1 = union_box(it)
    b = view[it.item_id]["box"]
    off = np.hypot((x0 + x1) / 2 - (b[0] + b[2]) / 2, (y0 + y1) / 2 - (b[1] + b[3]) / 2) / (b[2] - b[0])
    rows.append((off, np.mean(d), r5[it.item_id]["iou_tc"]))
a = np.array(rows)
for name, m in (("near<0.15", a[:, 0] < 0.15), ("far>0.3", a[:, 0] > 0.3), ("all", a[:, 0] >= 0)):
    print(name, int(m.sum()), "ceiling %.4f  R5 %.4f  misses %.4f" % (a[m, 1].mean(), a[m, 2].mean(), (a[m, 2] < 0.1).mean()))
