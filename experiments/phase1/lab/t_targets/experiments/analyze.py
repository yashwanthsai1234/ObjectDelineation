"""Summarise a results jsonl (optionally with tokens) into markdown tables.
    python experiments/analyze.py results/native.tok.jsonl [--by-object]"""
import gzip
import json
import sys
from collections import defaultdict

import numpy as np

recs = [json.loads(l) for l in (gzip.open if sys.argv[1].endswith(".gz") else open)(sys.argv[1], "rt")]
by_obj = "--by-object" in sys.argv
tok = "tokens" in recs[0]
order = []
for r in recs:
    if r["cfg"] not in order:
        order.append(r["cfg"])

rows = defaultdict(list)
items = defaultdict(list)
for r in recs:
    items[(r["cfg"], r["ds"])].append(r)
    for k, o in r["objs"].items():
        rows[(r["cfg"], k if by_obj else r["ds"])].append((r, k, o))


def pct(x):
    return f"{100 * np.mean(x):.2f}"


groups = sorted({g for _, g in rows}, key=lambda g: ["riga", "disc", "cup", "lidc", "nodule"].index(g))
print("| config | set | n | K5 topo % (raw) | K5 topo % | legal % (raw) | legal % | IoU min | IoU p5 | IoU med "
      "| verts med / p95 / max | Qwen tok/obj med / p95 / max | LFM2.5 med / p95 | Smol med | ok / rep / fb / fail "
      "| CPU s/item med / p95 |")
print("|" + "---|" * 16)
for c in order:
    for g in groups:
        L = rows.get((c, g))
        if not L:
            continue
        ds = L[0][0]["ds"]
        topo_raw = [r["raw"][k]["topo"] for r, k, o in L]
        topo = [o["topo"] for r, k, o in L]
        if ds == "riga" and (not by_obj or g == "cup"):
            topo_raw += [r["raw"]["rim"]["topo"] for r in items[(c, ds)]]
            topo += [r["derived_ok"]["rim"] for r in items[(c, ds)]]
        leg_raw = [r["raw"][k]["legal"] for r, k, o in L]
        leg = [o["legal"] for r, k, o in L]
        iou = np.array([o["iou"] for r, k, o in L])
        v = np.array([o["vertices"] for r, k, o in L])
        st = defaultdict(int)
        for r, k, o in L:
            st[o["status"]] += 1
        t = np.array([r["time"] for r in items[(c, ds)]])
        tq = tl = ts = "-"
        if tok:
            q = np.array([o["tokens"]["qwen3vl"] for r, k, o in L])
            lf = np.array([o["tokens"]["lfm25"] for r, k, o in L])
            sm = np.array([o["tokens"]["smolvlm"] for r, k, o in L])
            tq = f"{np.median(q):.0f} / {np.percentile(q, 95):.0f} / {q.max()}"
            tl = f"{np.median(lf):.0f} / {np.percentile(lf, 95):.0f}"
            ts = f"{np.median(sm):.0f}"
        print(f"| {c} | {g}{' (+rim)' if ds == 'riga' and (not by_obj or g == 'cup') else ''} | {len(L)} | "
              f"{pct(topo_raw)} | {pct(topo)} | {pct(leg_raw)} | {pct(leg)} | {iou.min():.4f} | "
              f"{np.percentile(iou, 5):.4f} | {np.median(iou):.4f} | {np.median(v):.0f} / {np.percentile(v, 95):.0f} / "
              f"{v.max()} | {tq} | {tl} | {ts} | {st['ok']} / {st['repaired']} / {st['fallback']} / {st['fail']} | "
              f"{np.median(t):.3f} / {np.percentile(t, 95):.3f} |")
if tok:
    print("\nItem-level target length (whole JSON answer), Qwen3-VL tokens: p50 / p95 / p99 / max")
    for c in order:
        for ds in ("riga", "lidc"):
            q = np.array([r["tokens"]["qwen3vl"] for r in items[(c, ds)]])
            if len(q):
                print(f"- {c} {ds}: {np.median(q):.0f} / {np.percentile(q, 95):.0f} / {np.percentile(q, 99):.0f} / {q.max()}")
