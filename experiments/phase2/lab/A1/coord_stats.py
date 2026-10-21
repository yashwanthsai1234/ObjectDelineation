import glob
import json
import os
import sys

import numpy as np

from objdel.geometry.text import SCHEMAS, parse

root = sys.argv[1]
st = {}
for f in sorted(glob.glob(os.path.join(root, "*", "*.jsonl"))):
    ds = f.split(os.sep)[-2]
    for line in open(f):
        r = json.loads(line)
        p = parse(r["text"], SCHEMAS[ds])
        W = r["W"]
        d = st.setdefault(ds, dict(vals=[], dmax=[], nparts=[], nholes=[], W=[], nv=[], steps=[]))
        for n in SCHEMAS[ds]:
            d["nparts"].append(len(p.objs[n]))
            d["nholes"].append(sum(len(h) for _, h in p.objs[n]))
            for o, hs in p.objs[n]:
                for a in [o, *hs]:
                    d["vals"].append(a.ravel())
                    d["W"].append(W)
                    d["nv"].append(len(a))
                    e = np.abs(np.diff(np.vstack([a, a[:1]]), axis=0)) // (1000 // W)
                    d["steps"].append(e.max(axis=1))
for ds, d in st.items():
    v = np.concatenate(d["vals"])
    steps = np.concatenate(d["steps"])
    np_ = np.array(d["nparts"])
    nh = np.array(d["nholes"])
    print(
        ds,
        "coords",
        v.size,
        "range",
        v.min(),
        v.max(),
        "==0:",
        (v == 0).sum(),
        "==1000:",
        (v == 1000).sum(),
        "| parts/obj dist",
        dict(zip(*np.unique(np_, return_counts=True))),
        "| holes>0:",
        (nh > 0).sum(),
        "of",
        nh.size,
        "| edge length px (chebyshev) p50/p95/max",
        np.percentile(steps, [50, 95]).tolist(),
        steps.max(),
        "| frac edges <=9px",
        round(float((steps <= 9).mean()), 4),
        "<=99px",
        round(float((steps <= 99).mean()), 4),
    )
