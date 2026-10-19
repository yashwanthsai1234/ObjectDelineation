"""Full LIDC multi-reader fixture (all 1,018 series, XML only) with the fixed parser.

For every nodule cluster with >= 2 readers and every slice with >= 2 reader masks:
uint8 stack (k, h, w) of the reader masks, cropped to the union bbox + 8 px.
Metadata per item: series, patient, split, cluster, z index (true slice index),
z (mm), k, readers in cluster, per-reader per-instance topology signature
(components < 5 px dropped and holes < 5 px filled, objdel size_filter; 8-conn
foreground / 4-conn background), and a disagreement flag.

    python -I build_fixture.py <out_npz> <out_meta_csv> <summary_json> [workers=4]
"""
import csv
import json
import sys
import warnings
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve().parent
TD = HERE.parent
REPO = Path("/Users/sreeyashwanthsaivenkatesh/Workspace/projects/CUDA/ObjectDelineation")
sys.path[:0] = [str(REPO / "src"), str(HERE)]

import numpy as np  # noqa: E402
from scipy import ndimage  # noqa: E402

import lidc_fix  # noqa: E402
from objdel.checker import betti, size_filter  # noqa: E402
from objdel.checker.topology import FG8  # noqa: E402

XML_ROOT = REPO / "data/raw/lidc_xml"
MIN_SIZE = 5
MARGIN = 8
ZCACHE = HERE / "slice_zvals.npz"


def slice_zvals_all():
    """series -> sorted slice z positions, from pylidc's bundled DB (cached)."""
    if ZCACHE.exists():
        z = np.load(ZCACHE, allow_pickle=False)
        return {k: z[k] for k in z.files}
    for a, t in (("int", int), ("bool", bool), ("float", float)):
        if not hasattr(np, a):
            setattr(np, a, t)
    warnings.filterwarnings("ignore")
    import pylidc as pl
    out = {s.series_instance_uid: np.asarray(s.slice_zvals, float) for s in pl.query(pl.Scan)}
    np.savez_compressed(ZCACHE, **out)
    return out


def signature(m):
    """Per-instance topology: 'n:h1,h2,..' (components, each one's hole count, sorted)."""
    f = size_filter(m, MIN_SIZE)
    lab, n = ndimage.label(f, structure=FG8)
    holes = sorted(betti(lab == i)[1] for i in range(1, n + 1))
    return f"{n}:" + ",".join(map(str, holes))


def one(args):
    series, sp, zvals, split = args
    got = lidc_fix.clusters(XML_ROOT / sp["xml_file"], zvals)
    zs = np.asarray(sorted(zvals))
    items, rows = {}, []
    summ = {"series": series, "sessions": got["sessions"], "clusters": len(got["clusters"]),
            "multi_reader_clusters": 0, **got["stats"]}
    for ci, cl in enumerate(got["clusters"]):
        if len(cl["readers"]) < 2:
            continue
        summ["multi_reader_clusters"] += 1
        for k in sorted({k for pr in cl["per_reader"].values() for k in pr}):
            readers = [r for r in cl["readers"] if k in cl["per_reader"][r]]
            if len(readers) < 2:
                continue
            masks = [cl["per_reader"][r][k] for r in readers]
            ys, xs = np.nonzero(np.any(masks, axis=0))
            y0, y1 = max(0, ys.min() - MARGIN), min(lidc_fix.GRID, ys.max() + MARGIN + 1)
            x0, x1 = max(0, xs.min() - MARGIN), min(lidc_fix.GRID, xs.max() + MARGIN + 1)
            stack = np.stack([m[y0:y1, x0:x1] for m in masks]).astype(np.uint8)
            sigs = [signature(m) for m in stack.astype(bool)]
            key = f"{series[-12:]}_c{ci}_k{k}__nodule"  # same "<item>__<structure>" layout as lidc_experts.npz
            items[key] = stack
            rows.append({"key": key, "series": series, "patient": sp["patient_id"], "split": split,
                         "cluster": ci, "z_index": k, "z_mm": round(float(zs[k]), 3), "k": len(readers),
                         "readers_in_cluster": len(cl["readers"]), "sessions": got["sessions"],
                         "reader_ids": ";".join(str(r + 1) for r in readers), "signatures": "|".join(sigs),
                         "disagree": len(set(sigs)) > 1, "areas": ";".join(str(int(s.sum())) for s in stack),
                         "crop_y0": int(y0), "crop_x0": int(x0), "h": int(y1 - y0), "w": int(x1 - x0),
                         "dup_reader_in_cluster": cl["dup_reader"]})
    return items, rows, summ


def main(out_npz, out_csv, out_json, workers=4):
    sp = lidc_fix.load_series_patient(TD / "series_patient.csv")
    split = {r["patient_id"]: r["split"] for r in csv.DictReader(open(TD / "splits/out/lidc_patients.csv"))}
    zv = slice_zvals_all()
    jobs = [(s, sp[s], zv[s], split[sp[s]["patient_id"]]) for s in sorted(sp)]
    with Pool(int(workers)) as pool:
        res = pool.map(one, jobs, chunksize=4)
    items, rows, summs = {}, [], []
    for it, rs, sm in res:
        assert not set(it) & set(items)
        items.update(it)
        rows.extend(rs)
        summs.append(sm)
    tmp = Path(str(out_npz) + ".part.npz")
    np.savez_compressed(tmp, **items)
    tmp.replace(out_npz)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    report = {"series": len(summs), "items": len(rows),
              "single_point_rois_dropped": sum(s["single_point_rois_dropped"] for s in summs),
              "exclusion_contours_applied": sum(s["exclusions"] for s in summs),
              "annotations_empty_after_raster": sum(s["empty_after_raster"] for s in summs),
              "max_dz_mm": max(s["max_dz"] for s in summs),
              "series_dz_gt_0.01": sum(s["max_dz"] > 0.01 for s in summs),
              "multi_reader_clusters": sum(s["multi_reader_clusters"] for s in summs),
              "per_split": {}}
    nod_dis = Counter()
    for r in rows:
        if r["disagree"]:
            nod_dis[(r["split"], r["series"], r["cluster"])] = 1
    nods = {(r["split"], r["series"], r["cluster"]) for r in rows}
    for name in ("train", "val", "test"):
        rs = [r for r in rows if r["split"] == name]
        report["per_split"][name] = {
            "patients": len({r["patient"] for r in rs}), "nodules": sum(1 for n in nods if n[0] == name),
            "slices": len(rs), "disagreement_slices": sum(r["disagree"] for r in rs),
            "disagreement_nodules": sum(1 for n in nod_dis if n[0] == name),
            "k_hist": dict(sorted(Counter(r["k"] for r in rs).items())),
            "k_hist_disagree": dict(sorted(Counter(r["k"] for r in rs if r["disagree"]).items()))}
    report["signature_hist"] = dict(Counter(s for r in rows for s in r["signatures"].split("|")).most_common(12))
    Path(out_json).write_text(json.dumps(report, indent=1, default=str))
    print(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main(*sys.argv[1:])
