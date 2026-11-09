"""Topology census of every LIDC item (train / val / test) on the model frame, for the topology probe (H-A1-7).

    .venv/bin/python -I -m experiments.phase2.lab.B3.census --out <scratch>/census.jsonl --summary <lab>/census.json

Per item: patient (group_id, the split unit of objdel.data.splits), nodule (series + cluster: the slices of one
nodule are correlated and form the bootstrap cluster), each rater's size-filtered signature at the reward's rule
min_area(|union|) and its raw signature, the voting raters (reward.cache.voters), the union's equivalent diameter.
Reads untrusted data: run with python -I.
"""
import argparse
import json
from collections import Counter

import numpy as np

from objdel.checker import signature, topology
from objdel.data import store
from objdel.reward.cache import voters
from objdel.targets import frame

SIMPLE = (0,)


def item_row(it):
    f = frame(it)
    st = f.masks["nodule"]
    union = st.any(0)
    ma = topology.min_area(int(union.sum()))
    sigs = [signature(m, ma) for m in st]
    raw = [signature(m) for m in st]
    vote = voters(sigs)
    vs = [sigs[i] for i in vote]
    return {"item_id": it.item_id, "split": it.split, "patient": it.group_id,
            "nodule": it.item_id.rsplit("/", 1)[0], "k": len(sigs), "voters": list(vote),
            "sigs": [list(s) for s in sigs], "raw": [list(s) for s in raw], "min_area": ma, "W": f.W,
            "diam": float(2 * np.sqrt(union.sum() / np.pi)),
            "rater_diam": [float(2 * np.sqrt(m.sum() / np.pi)) for m in st],
            "nonsimple": any(s != SIMPLE for s in vs),
            "unanimous_nonsimple": len(set(vs)) == 1 and vs[0] != SIMPLE,
            "majority_nonsimple": sum(s != SIMPLE for s in vs) * 2 > len(vs),
            "raw_nonsimple": any(tuple(s) != SIMPLE for s in raw)}


def summarize(rows):
    out = {}
    for split in ("train", "val", "test"):
        rs = [r for r in rows if r["split"] == split]
        ns = [r for r in rs if r["nonsimple"]]
        tg = [(r, i) for r in rs for i in r["voters"]]
        out[split] = {
            "items": len(rs), "patients": len({r["patient"] for r in rs}), "nodules": len({r["nodule"] for r in rs}),
            "targets_voting": len(tg),
            "nonsimple_targets": sum(tuple(r["sigs"][i]) != SIMPLE for r, i in tg),
            "raw_nonsimple_targets": sum(tuple(s) != SIMPLE for r in rs for s in r["raw"]),
            "raw_nonsimple_items": sum(r["raw_nonsimple"] for r in rs),
            "nonsimple_items": len(ns), "nonsimple_nodules": len({r["nodule"] for r in ns}),
            "nonsimple_patients": len({r["patient"] for r in ns}),
            "unanimous_nonsimple_items": sum(r["unanimous_nonsimple"] for r in rs),
            "majority_nonsimple_items": sum(r["majority_nonsimple"] for r in rs),
            "hole_items": sum(any(sum(r["sigs"][i]) > 0 for i in r["voters"]) for r in ns),
            "multipart_items": sum(any(len(r["sigs"][i]) > 1 for i in r["voters"]) for r in ns),
            "nonsimple_items_per_patient_top10": sorted(Counter(r["patient"] for r in ns).values())[::-1][:10],
            "unanimous_sigs": dict(Counter(str(r["sigs"][r["voters"][0]]) for r in rs if r["unanimous_nonsimple"])),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--store", default="data/store")
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()
    rows = []
    with open(a.out, "w") as fh:
        for split in ("train", "val", "test"):
            for it in store.load("lidc", split, a.store):
                r = item_row(it)
                rows.append(r)
                fh.write(json.dumps(r) + "\n")
            print(split, len(rows), flush=True)
    s = summarize(rows)
    print(json.dumps(s, indent=1))
    with open(a.summary, "w") as fh:
        json.dump(s, fh, indent=1)


if __name__ == "__main__":
    main()
