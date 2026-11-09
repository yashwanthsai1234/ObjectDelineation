"""IoU ceiling / legality of re-quantizing targets to a coarser grid G (0..G). -I coarse_grid.py data/targets"""

import json
import os
import sys

import numpy as np

from objdel.geometry import SCHEMAS, parse, to_mask

root = sys.argv[1]
rng = np.random.default_rng(1)
for ds in ("riga", "lidc"):
    recs = [json.loads(l) for l in open(os.path.join(root, ds, "val.jsonl"))]
    recs = [recs[i] for i in rng.choice(len(recs), min(1500, len(recs)), replace=False)]
    for G in (100, 250, 500):
        ious, legal, dup = [], 0, 0
        n = 0
        for r in recs:
            p = parse(r["text"], SCHEMAS[ds])
            W = r["W"]
            for name in SCHEMAS[ds]:
                parts = p.objs[name]
                if not parts:
                    continue
                q = lambda a: (np.rint(a * G / 1000) * 1000 / G).round().astype(np.int64)
                parts2 = [(q(o), [q(h) for h in hs]) for o, hs in parts]
                ref = to_mask(parts, W, W).mask
                g = to_mask(parts2, W, W)
                ious.append((ref & g.mask).sum() / max(1, (ref | g.mask).sum()))
                legal += g.legal
                n += 1
        ious = np.array(ious)
        print(
            f"{ds} G={G}: n={n} IoU mean {ious.mean():.4f} p5 {np.percentile(ious, 5):.4f} min {ious.min():.4f} legal {legal / n:.4f}"
        )
