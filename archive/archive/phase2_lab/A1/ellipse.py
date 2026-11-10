import json
import os
import sys

import cv2
import numpy as np

from objdel.geometry import parse, to_mask

root = sys.argv[1]
res = {"disc": [], "cup": []}
for l in open(os.path.join(root, "riga", "val.jsonl")):
    r = json.loads(l)
    p = parse(r["text"], ("disc", "cup"))
    W = r["W"]
    for n in ("disc", "cup"):
        m = to_mask(p.objs[n], W, W).mask
        cs, _ = cv2.findContours(m.astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        c = max(cs, key=cv2.contourArea)
        if len(c) < 5:
            continue
        e = cv2.fitEllipse(c)
        em = np.zeros_like(m, np.uint8)
        cv2.ellipse(em, e, 1, -1)
        em = em.astype(bool)
        res[n].append((m & em).sum() / (m | em).sum())
for n, v in res.items():
    v = np.array(v)
    print(n, len(v), f"ellipse IoU mean {v.mean():.4f} median {np.median(v):.4f} p5 {np.percentile(v, 5):.4f}")
