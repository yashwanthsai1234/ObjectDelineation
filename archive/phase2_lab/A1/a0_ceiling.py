"""Ceiling of A0 variants vs A0' on the real targets. -I a0_ceiling.py data/targets"""

import json
import os
import sys

import numpy as np

from objdel.checker.polycore import signed_area2
from objdel.geometry import parse, to_mask

root = sys.argv[1]
for split in ("train", "val"):
    f = os.path.join(root, "lidc", f"{split}.jsonl")
    items = {}
    stats = dict(n=0, multi=0, holes=0, both=0)
    ious = {"A0single": [], "A0multi": []}
    for l in open(f):
        r = json.loads(l)
        p = parse(r["text"], ("nodule",))
        parts = p.objs["nodule"]
        W = r["W"]
        nh = sum(len(h) for _, h in parts)
        mp = len(parts) > 1
        stats["n"] += 1
        stats["multi"] += mp
        stats["holes"] += nh > 0
        stats["both"] += mp and nh > 0
        it = items.setdefault(r["item_id"], [0, 0])
        it[0] += 1
        it[1] += mp or nh > 0
        if mp or nh:
            full = to_mask(parts, W, W).mask
            big = max(parts, key=lambda q: signed_area2(q[0]))
            for name, pp in (("A0single", [(big[0], [])]), ("A0multi", [(o, []) for o, _ in parts])):
                m = to_mask(pp, W, W).mask
                ious[name].append((full & m).sum() / (full | m).sum())
    any_t = sum(1 for v in items.values() if v[1] > 0)
    all_t = sum(1 for v in items.values() if v[1] == v[0])
    print(
        split,
        stats,
        "items",
        len(items),
        "items with >=1 rater non-simple",
        any_t,
        "items all raters non-simple",
        all_t,
    )
    for k, v in ious.items():
        v = np.array(v)
        print(
            "  ",
            k,
            f"on non-simple targets: IoU mean {v.mean():.3f} median {np.median(v):.3f} p5 {np.percentile(v, 5):.3f}",
        )
