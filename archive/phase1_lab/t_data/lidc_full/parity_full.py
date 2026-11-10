"""Parity of the fixed parser (lidc_fix) vs pylidc 0.2.3 on N scans, with the aligned outline rule.

pylidc-aligned mask of an annotation on a slice (rule P2 expressed with pylidc parts):
    (pylidc inclusion interiors  U  inclusion contour points)  minus  (exclusion interiors U exclusion contour points)
where interiors are pylidc's own matplotlib Path.contains_points test on pixel
centres (what Annotation.boolean_mask uses). Ours is PIL polygon fill + outline.

    python -I parity_full.py <n_scans|all> <out_json> [workers=4] [seed=0]
"""
import json
import random
import sys
import warnings
from collections import Counter, defaultdict
from multiprocessing import get_context
from pathlib import Path

HERE = Path(__file__).resolve().parent
TD = HERE.parent
REPO = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(REPO / "src"), str(HERE)]

import numpy as np  # noqa: E402

import lidc_fix  # noqa: E402

XML_ROOT = REPO / "data/raw/lidc_xml"


def shim():
    for a, t in (("int", int), ("bool", bool), ("float", float)):
        if not hasattr(np, a):
            setattr(np, a, t)
    warnings.filterwarnings("ignore")


def key_of(rois):
    return frozenset((round(z, 2), bool(inc), tuple(map(tuple, pts))) for z, inc, pts in rois if len(pts) > 1)


def iou(a, b):
    u = np.logical_or(a, b).sum()
    return 1.0 if u == 0 else float(np.logical_and(a, b).sum() / u)


def pylidc_aligned(pa, zvals):
    """{slice index: 512x512 mask} for one pylidc annotation under rule P2."""
    from matplotlib.path import Path as MPath
    ii, jj = np.mgrid[0:512, 0:512]
    pts = np.c_[ii.ravel(), jj.ravel()]
    inc, exc = defaultdict(lambda: np.zeros((512, 512), bool)), defaultdict(list)
    for c in pa.contours:
        C = c.to_matrix(include_k=False)  # (i=row, j=col)
        k = int(np.abs(zvals - c.image_z_position).argmin())
        if c.inclusion:
            inc[k] |= _poly(C, pts, MPath)
        else:
            exc[k].append(C)
    for k, Cs in exc.items():
        for C in Cs:
            if k in inc:
                inc[k] &= ~_poly(C, pts, MPath)
    return {k: m for k, m in inc.items() if m.any()}


def _poly(C, pts, MPath):
    if (C[0] != C[-1]).any():
        Cc = np.vstack([C, C[:1]])
    else:
        Cc = C
    i0, j0 = C.min(0)
    i1, j1 = C.max(0)
    m = np.zeros((512, 512), bool)
    sub = (pts[:, 0] >= i0) & (pts[:, 0] <= i1) & (pts[:, 1] >= j0) & (pts[:, 1] <= j1)
    inside = MPath(Cc, closed=True).contains_points(pts[sub])
    m.ravel()[np.flatnonzero(sub)[inside]] = True
    m[C[:, 0], C[:, 1]] = True  # contour points belong to the polygon (rule P2)
    return m


def analyse(args):
    series, xml = args
    shim()
    import contextlib
    import io
    import pylidc as pl
    scan = pl.query(pl.Scan).filter(pl.Scan.series_instance_uid == series).first()
    zvals = np.asarray(scan.slice_zvals, float)
    ours = lidc_fix.clusters(XML_ROOT / xml, zvals)
    anns = ours["anns"]
    by_key = defaultdict(list)
    for i, a in enumerate(anns):
        by_key[key_of(a["rois"])].append(i)
    panns = list(scan.annotations)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        pcl = scan.cluster_annotations(verbose=True)
    p2o, used, unmatched = {}, Counter(), 0
    for pa in panns:
        rois = [(c.image_z_position, bool(c.inclusion), [tuple(int(v) for v in p.split(",")) for p in c.coords.split("\n")])
                for c in pa.contours]
        k = key_of(rois)
        if by_key.get(k):
            p2o[pa.id] = by_key[k][min(used[k], len(by_key[k]) - 1)]
            used[k] += 1
        else:
            unmatched += 1
    our_part = {frozenset(c["member_ids"]) for c in ours["clusters"]}
    pyl_part = {frozenset(p2o[a.id] for a in c if a.id in p2o) for c in pcl} - {frozenset()}
    rs = lambda part: sorted(tuple(sorted({anns[i]["reader"] for i in c})) for c in part)  # noqa: E731
    out = {"series": series, "our_anns": len(anns), "pylidc_anns": len(panns), "pylidc_unmatched": unmatched,
           "our_unmatched": len(set(range(len(anns))) - set(p2o.values())),
           "our_nodules": len(ours["clusters"]), "pylidc_nodules": len(pcl),
           "partition_equal": our_part == pyl_part, "reader_sets_equal": rs(our_part) == rs(pyl_part),
           "min_tol_hit": "Failed to reduce" in buf.getvalue(), "ious": [], "exact": 0, "slices": 0,
           "exclusion_slices": 0, "exclusion_exact": 0}
    pmask = {p2o[pa.id]: pylidc_aligned(pa, zvals) for pa in panns if pa.id in p2o}
    for c in ours["clusters"]:
        if frozenset(c["member_ids"]) not in pyl_part:
            continue
        for r in c["readers"]:
            ids = [i for i in c["member_ids"] if anns[i]["reader"] == r]
            ks = set().union(*(set(anns[i]["slices"]) | set(pmask[i]) for i in ids))
            for k in ks:
                o = np.zeros((512, 512), bool)
                p = np.zeros((512, 512), bool)
                for i in ids:
                    o |= anns[i]["slices"].get(k, False)
                    p |= pmask[i].get(k, False)
                v = iou(o, p)
                has_exc = any(not inc and abs(z - zvals[k]) < 0.01 for i in ids for z, inc, _ in anns[i]["rois"])
                out["ious"].append(round(v, 5))
                out["slices"] += 1
                out["exact"] += v == 1.0
                out["exclusion_slices"] += has_exc
                out["exclusion_exact"] += has_exc and v == 1.0
    return out


def main(n, out_json, workers=4, seed=0):
    sp = lidc_fix.load_series_patient(TD / "series_patient.csv")
    series = sorted(sp)
    if n != "all":
        rng = random.Random(int(seed))
        series = sorted(rng.sample(series, int(n)))
    with get_context("spawn").Pool(int(workers)) as pool:
        res = pool.map(analyse, [(s, sp[s]["xml_file"]) for s in series], chunksize=4)
    ious = np.array([v for r in res for v in r["ious"]])
    s = {"scans": len(res),
         "annotations": {"ours": sum(r["our_anns"] for r in res), "pylidc": sum(r["pylidc_anns"] for r in res),
                         "pylidc_unmatched": sum(r["pylidc_unmatched"] for r in res),
                         "ours_unmatched": sum(r["our_unmatched"] for r in res)},
         "nodule_count_equal_scans": sum(r["our_nodules"] == r["pylidc_nodules"] for r in res),
         "partition_equal_scans": sum(r["partition_equal"] for r in res),
         "reader_sets_equal_scans": sum(r["reader_sets_equal"] for r in res),
         "min_tol_hit_scans": sum(r["min_tol_hit"] for r in res),
         "mismatch_scans": [r["series"][-12:] for r in res if not (r["partition_equal"] and r["reader_sets_equal"]
                                                                  and r["our_nodules"] == r["pylidc_nodules"])],
         "reader_slices": int(ious.size), "exact_pct": round(100 * float((ious == 1).mean()), 3),
         "iou": {"min": float(ious.min()), "p1": float(np.percentile(ious, 1)), "p5": float(np.percentile(ious, 5)),
                 "median": float(np.median(ious)), "mean": float(ious.mean())},
         "iou_ge_0.95_pct": round(100 * float((ious >= 0.95).mean()), 3),
         "exclusion_slices": sum(r["exclusion_slices"] for r in res),
         "exclusion_exact": sum(r["exclusion_exact"] for r in res)}
    Path(out_json).write_text(json.dumps(s, indent=1))
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:])
