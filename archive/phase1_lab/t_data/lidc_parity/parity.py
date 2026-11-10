"""LIDC parity: our XML parser + clustering (experiments.phase0.lidc_census) vs pylidc 0.2.3.

pylidc's bundled SQLite DB (annotations, contours, slice z values, patient ids)
is enough for clustering and Annotation.boolean_mask; no DICOM is read.

    ../../.venv/bin/python -I parity.py [n_slices_series=30] [n_random=6]
"""
import contextlib
import csv
import io
import json
import random
import sys
import warnings
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from multiprocessing import get_context
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path[:0] = [str(REPO), str(REPO / "src")]
OUT = Path(__file__).resolve().parent
XML_ROOT = REPO / "data/raw/lidc_xml"
SLICES = REPO / "data/lidc_slices"
NS = {"n": "http://www.nih.gov"}

# numpy-2 shim for pylidc 0.2.3 (uses removed aliases); set before pylidc code runs.
for _a, _t in (("int", int), ("bool", bool), ("float", float)):
    if not hasattr(np, _a):
        setattr(np, _a, _t)
warnings.filterwarnings("ignore")


QUIRKS = [
    "pylidc 0.2.3 cluster_annotations default tol = scan.slice_thickness (mm) applied to distances between contour points in voxel-INDEX units (i, j, k), not mm; docs/notes saying tol=pixel_spacing are wrong for 0.2.3",
    "pylidc returns a bare list (ignores return_distance_matrix) when a scan has <= 1 annotation",
    "pylidc needs numpy<2 aliases (np.int/np.bool/np.float) and pkg_resources (setuptools<81); DB-only use (clustering, boolean_mask, slice_zvals, patient_id) works without DICOM",
    "pylidc boolean_mask(include_contour_points=True) does NOT add outline pixels, it only skips removing them; ours = that OR outline pixels (exact on 99.57% of slices)",
    "exclusions: ours removes exclusion interior and keeps its outline pixels; pylidc removes interior and the exclusion contour points",
    "single-point ROIs inside >=3 mm nodules (15 ROIs / 8 annotations in sample, incl. 1 exclusion): pylidc DB drops them, our parser rasterizes them as isolated 1-px blobs (spurious b0 component; census MIN_SIZE_2D=5 hides it)",
    "6 series appear in 2-6 XML files each (188/082-086, 087-090, 091-093, 094-096, 189/028-033, 173-174): the reading sessions are byte-identical duplicates, so census first-file dedupe loses nothing; 161 has a resubmitted correction (session 1 differs) which census and pylidc both use",
    "188/011.xml (LIDC-IDRI-0566) has 8 readingSessions; sessions 5-8 hold only <3 mm/non-characterised points at z=-94.4,-71.0 which lie on no slice of this series (max 3.5 mm off, -71.0 beyond the last slice) -> misfiled reads; census counts sessions=8 and excludes it from the 4-session stats",
    "189/110.xml (LIDC-IDRI-0700) has 3 readingSessions",
    "our z_index (step = min diff of annotated z) disagrees with true slice-index offsets in 6/67 scans (overlapping recon: slice_thickness 1.25 vs spacing 0.625, or the misfiled 011 z values); 2D masks unaffected, 3D adjacency/betti3d gaps are",
    "044.xml / 191.xml exist in 4 folders each (8 files); all parse and match pylidc exactly; 186/044 and 188/044 contain only single-point (<3 mm) nodules without characteristics -> 0 contoured nodules on both sides",
    "servicingRadiologistID is anonymised ('anon'/''/'anonymous'/numeric); no duplicated named IDs within a file; 255 reading sessions contain no nodule at all (readers who marked nothing)",
    "XML contains no PatientID: series->patient mapping must come from pylidc DB (or DICOM); 1018 series / 1010 patients; 8 patients have two series",
]

def series_files():
    """series uid -> list of XML paths (sorted, as lidc_census.main sees them) and census-chosen file."""
    by = defaultdict(list)
    for p in sorted(XML_ROOT.rglob("*.xml")):
        r = ET.parse(p).getroot()
        if r.tag.endswith("LidcReadMessage"):
            by[r.find("n:ResponseHeader/n:SeriesInstanceUid", NS).text.strip()].append(p)
    chosen = {}
    for s, ps in by.items():  # census rule: first file, unless a "correction" file exists
        chosen[s] = next((p for p in ps if "correction" in p.name), ps[0])
    return by, chosen


def ann_key(rois):
    """Exact contour identity: frozenset of (z, inclusion, point tuple).

    Single-point ROIs are ignored: pylidc's DB drops 1-point contours that sit inside
    a >=3 mm nodule, our parser keeps them (and rasterizes them as a 1-pixel blob)."""
    return frozenset((round(z, 2), inc, tuple(map(tuple, pts))) for z, inc, pts in rois if len(pts) > 1)


def iou(a, b):
    u = np.logical_or(a, b).sum()
    return 1.0 if u == 0 else float(np.logical_and(a, b).sum() / u)


def analyse(args):
    series, xml_path = args
    import pylidc as pl
    from experiments.phase0 import lidc_census as lc

    scan = pl.query(pl.Scan).filter(pl.Scan.series_instance_uid == series).first()
    out = {"series": series, "xml": str(Path(xml_path).relative_to(XML_ROOT)), "patient_id": scan.patient_id,
           "slice_thickness": scan.slice_thickness, "pixel_spacing": scan.pixel_spacing,
           "slice_spacing": float(scan.slice_spacing)}
    _, sessions = lc.parse(xml_path)
    got = lc.nodule_clusters(xml_path)
    _, _, zi, ours, _ = got
    out["sessions"] = len(sessions)

    # --- our annotations, in the order nodule_clusters builds them ---
    our_anns = []  # (reader, noduleID, key, rois, large, nonempty)
    for r, nods in enumerate(sessions):
        for nod in nods:
            k = ann_key(nod["rois"])
            nonempty = None
            if nod["large"]:
                sl, _ = lc.rasterize(nod, zi)
                nonempty = any(m.any() for m in sl.values())
            our_anns.append({"reader": r, "nid": nod["id"], "key": k, "rois": nod["rois"],
                             "large": nod["large"], "nonempty": nonempty,
                             "npts": max(len(p) for _, _, p in nod["rois"]),
                             "single_pt_rois": sum(len(p) == 1 for _, _, p in nod["rois"]) if nod["large"] else 0,
                             "has_chars": None})
    # characteristics presence per nodule (to explain the "large" rule)
    root = ET.parse(xml_path).getroot()
    flat = [nod for s in root.findall("n:readingSession", NS) for nod in s.findall("n:unblindedReadNodule", NS)
            if nod.findall("n:roi", NS) and any(roi.findall("n:edgeMap", NS) for roi in nod.findall("n:roi", NS))]
    for a, nod in zip(our_anns, flat):
        a["has_chars"] = nod.find("n:characteristics", NS) is not None
    # map our cluster members -> our ann index via (reader, voxel keys)
    vox_of = {}
    for i, a in enumerate(our_anns):
        if a["large"] and a["nonempty"]:
            sl, _ = lc.rasterize({"rois": a["rois"]}, zi)
            sl = {k: m for k, m in sl.items() if m.any()}
            vox_of.setdefault((a["reader"], frozenset(map(tuple, lc.voxels(sl).tolist()))), []).append(i)
    our_part = []
    for cl in ours:
        ids = []
        for m in cl["members"]:
            cand = vox_of[(m["reader"], frozenset(m["keys"]))]
            ids.append(cand[0] if len(cand) == 1 else tuple(cand))
        our_part.append(ids)

    # --- pylidc ---
    panns = list(scan.annotations)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        # quirk: for N <= 1 annotations pylidc ignores return_distance_matrix and returns a bare list
        pclusters = scan.cluster_annotations(verbose=True)
    out["pylidc_min_tol_hit"] = "Failed to reduce" in buf.getvalue()
    out["pylidc_cluster_gt4"] = sum(len(c) > 4 for c in pclusters)
    # match pylidc annotations to ours by exact contour identity
    by_key = defaultdict(list)
    for i, a in enumerate(our_anns):
        by_key[a["key"]].append(i)
    p2o, unmatched_p = {}, []
    used = Counter()
    for pa in panns:
        rois = [(c.image_z_position, bool(c.inclusion),
                 [tuple(int(v) for v in pt.split(",")) for pt in c.coords.split("\n")]) for c in pa.contours]
        k = ann_key(rois)
        cand = by_key.get(k, [])
        if cand:
            j = cand[min(used[k], len(cand) - 1)]
            used[k] += 1
            p2o[pa.id] = j
            if our_anns[j]["nid"] != pa._nodule_id:
                out.setdefault("nodule_id_mismatch", 0)
                out["nodule_id_mismatch"] += 1
        else:
            unmatched_p.append(pa.id)
    out["pylidc_anns"] = len(panns)
    out["pylidc_anns_unmatched"] = len(unmatched_p)
    our_used = [i for i, a in enumerate(our_anns) if a["large"] and a["nonempty"]]
    matched_o = set(p2o.values())
    out["our_anns"] = len(our_used)
    only_ours = [i for i in our_used if i not in matched_o]
    only_pyl_in_ours = [our_anns[p2o[pa.id]] for pa in panns if pa.id in p2o and p2o[pa.id] not in set(our_used)]
    out["our_only_anns"] = len(only_ours)
    out["single_point_rois_in_large"] = sum(our_anns[i]["single_pt_rois"] for i in our_used)
    out["anns_with_single_point_rois"] = sum(our_anns[i]["single_pt_rois"] > 0 for i in our_used)
    out["pylidc_only_anns"] = len(only_pyl_in_ours) + len(unmatched_p)
    out["pylidc_only_reasons"] = dict(Counter(
        ("no_characteristics" if not a["has_chars"] else "single_point" if a["npts"] == 1 else
         "empty_raster" if a["large"] and not a["nonempty"] else "other") for a in only_pyl_in_ours))
    out["our_only_detail"] = [(our_anns[i]["reader"], our_anns[i]["nid"], our_anns[i]["npts"]) for i in only_ours]

    # partitions on the common annotation set
    common = set(p2o.values()) & set(our_used)

    def norm(part):
        s = set()
        for c in part:
            flat_ids = []
            for x in c:
                flat_ids.extend(x if isinstance(x, tuple) else [x])
            f = frozenset(i for i in flat_ids if i in common)
            if f:
                s.add(f)
        return s

    pid_to_ann = {pa.id: pa for pa in panns}
    pyl_part = [[p2o[a.id] for a in c if a.id in p2o] for c in pclusters]
    op, pp = norm(our_part), norm(pyl_part)
    out["our_nodules"] = len(ours)
    out["pylidc_nodules"] = len(pclusters)
    out["nodules_on_common"] = (len(op), len(pp))
    out["partition_equal"] = op == pp
    out["partition_diff"] = {"ours_not_pylidc": [sorted(x) for x in op - pp],
                             "pylidc_not_ours": [sorted(x) for x in pp - op]}
    # reader sets per cluster (on equal clusters they coincide by construction); duplicated readers
    out["pylidc_clusters_dup_reader"] = sum(
        1 for c in pyl_part if len(c) != len({our_anns[i]["reader"] for i in c}))
    out["our_clusters_dup_reader"] = sum(
        1 for c in ours if len(c["members"]) != len({m["reader"] for m in c["members"]}))
    rs_o = sorted(tuple(sorted({our_anns[i]["reader"] for i in c})) for c in op)
    rs_p = sorted(tuple(sorted({our_anns[i]["reader"] for i in c})) for c in pp)
    out["reader_sets_equal"] = rs_o == rs_p

    # --- z matching ---
    zvals = scan.slice_zvals
    out["scan_dup_z"] = int(len(zvals) - len(np.unique(np.round(zvals, 3))))
    xz = sorted({round(z, 2) for nods in sessions for n in nods for z, _, _ in n["rois"]})
    kz = {z: int(np.abs(zvals - z).argmin()) for z in xz}
    out["z_max_mismatch_mm"] = float(max((abs(zvals[kz[z]] - z) for z in xz), default=0.0))
    out["z_not_on_slice"] = int(sum(abs(zvals[kz[z]] - z) > 0.01 for z in xz))
    # our z index (from XML z values) vs true slice index offsets
    if xz:
        z0 = xz[0]
        out["z_index_offset_mismatch"] = sum((zi[z] - zi[z0]) != (kz[z] - kz[z0]) for z in xz)
    # SOP UID <-> z consistency inside the XML
    sop_z = defaultdict(set)
    z_sop = defaultdict(set)
    for roi in root.iter("{http://www.nih.gov}roi"):
        z = round(float(roi.find("n:imageZposition", NS).text), 2)
        sop = roi.find("n:imageSOP_UID", NS).text.strip()
        sop_z[sop].add(z)
        z_sop[z].add(sop)
    out["sop_multi_z"] = sum(len(v) > 1 for v in sop_z.values())
    out["z_multi_sop"] = sum(len(v) > 1 for v in z_sop.values())
    # pylidc dicom_file_name per contour vs the scan's file list (consistency without DICOM)
    out["exclusion_contours_xml"] = sum(1 for nods in sessions for n in nods for _, inc, _ in n["rois"] if not inc)
    out["exclusion_contours_pylidc"] = sum(1 for pa in panns for c in pa.contours if not c.inclusion)

    # --- masks: per reader per slice within equal clusters ---
    slices_rows = []
    our_raster = {}
    for i in common:
        sl, _ = lc.rasterize({"rois": our_anns[i]["rois"]}, zi)
        our_raster[i] = sl  # {z index: mask}
    inv_zi = {v: z for z, v in zi.items()}
    pmask = {}
    for pa in panns:
        if pa.id not in p2o or p2o[pa.id] not in common:
            continue
        bb = pa.bbox_matrix()
        res = {}
        for flag in (True, False):
            m = pa.boolean_mask(include_contour_points=flag)
            per_z = {}
            for kk in range(m.shape[2]):
                full = np.zeros((512, 512), bool)
                full[bb[0, 0]:bb[0, 1] + 1, bb[1, 0]:bb[1, 1] + 1] = m[:, :, kk]
                per_z[round(float(zvals[bb[2, 0] + kk]), 2)] = full
            res[flag] = per_z
        pmask[p2o[pa.id]] = res
    # per-reader unions per cluster
    for c in op & pp:
        for r in sorted({our_anns[i]["reader"] for i in c}):
            ids = [i for i in c if our_anns[i]["reader"] == r]
            zs = set()
            for i in ids:
                zs |= {inv_zi[k] for k in our_raster[i]}
                zs |= set(pmask[i][True])
            for z in sorted(zs):
                o = np.zeros((512, 512), bool)
                for i in ids:
                    if zi.get(z) in our_raster[i]:
                        o |= our_raster[i][zi[z]]
                pt = np.zeros((512, 512), bool)
                pf = np.zeros((512, 512), bool)
                for i in ids:
                    pt |= pmask[i][True].get(z, np.zeros((512, 512), bool))
                    pf |= pmask[i][False].get(z, np.zeros((512, 512), bool))
                # contour (outline) points of this reader on this slice
                outline = np.zeros((512, 512), bool)
                has_exc = False
                for i in ids:
                    for zz, inc, pts in our_anns[i]["rois"]:
                        if round(zz, 2) == z:
                            outline |= lc.draw(pts, outline_only=True) if len(pts) > 1 else lc.draw(pts)
                            has_exc |= not inc
                slices_rows.append({
                    "series": series[-12:], "reader": r, "z": z, "multi_ann": len(ids) > 1, "exclusion": has_exc,
                    "area_ours": int(o.sum()), "area_pyl_incl": int(pt.sum()), "area_pyl_excl": int(pf.sum()),
                    "iou_incl": iou(o, pt), "iou_excl": iou(o, pf),
                    "iou_ours_minus_outline_vs_excl": iou(o & ~outline, pf),
                    "iou_ours_vs_incl_plus_outline": iou(o, pt | outline),
                    "diff_px_incl": int((o ^ pt).sum()), "diff_px_excl": int((o ^ pf).sum()),
                    "outline_px": int(outline.sum()),
                })
    out["n_slice_rows"] = len(slices_rows)
    return out, slices_rows


def select(by, seed=0, n_slices=30, n_random=6):
    import pylidc as pl
    rng = random.Random(seed)
    have = sorted(p.stem for p in SLICES.glob("*.npz"))
    pick = {s: "lidc_slices" for s in sorted(rng.sample(have, min(n_slices, len(have))))}
    for s, ps in by.items():
        if any(p.name in ("044.xml", "191.xml") for p in ps):
            pick.setdefault(s, "xml_044_191")
        if len(ps) > 1:
            pick.setdefault(s, "multi_file_series")
    pp = defaultdict(list)
    for sc in pl.query(pl.Scan).all():
        pp[sc.patient_id].append(sc.series_instance_uid)
    two = {k: sorted(v) for k, v in pp.items() if len(v) > 1}
    for v in two.values():
        for s in v:
            pick.setdefault(s, "two_scan_patient")
    for s in ("011.xml", "110.xml"):  # 8- and 3-session files
        for ser, ps in by.items():
            if any(p.parent.name in ("188", "189") and p.name == s for p in ps):
                pass
    rest = sorted(set(by) - set(pick))
    for s in rng.sample(rest, n_random):
        pick[s] = "random"
    return pick, two, pp


def main(n_slices=30, n_random=6):
    import pylidc as pl
    by, chosen = series_files()
    pick, two, pp = select(by, n_slices=int(n_slices), n_random=int(n_random))
    # odd session counts
    for ser, ps in by.items():
        if any((p.parent.name, p.name) in (("188", "011.xml"), ("189", "110.xml")) for p in ps):
            pick.setdefault(ser, "odd_session_count")
    # series -> patient for every series in our XML set
    db = {sc.series_instance_uid: sc.patient_id for sc in pl.query(pl.Scan).all()}
    with open(OUT.parent / "series_patient.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["series_uid", "patient_id", "xml_file", "all_xml_files"])
        for s in sorted(by):
            w.writerow([s, db.get(s, ""), str(chosen[s].relative_to(XML_ROOT)),
                        ";".join(str(p.relative_to(XML_ROOT)) for p in by[s])])
    jobs = [(s, str(chosen[s])) for s in sorted(pick)]
    with get_context("spawn").Pool(7) as pool:
        res = pool.map(analyse, jobs, chunksize=1)
    scans = []
    rows = []
    for (o, r) in res:
        o["reason"] = pick[o["series"]]
        scans.append(o)
        rows.extend(r)
    (OUT / "scans.json").write_text(json.dumps(scans, indent=1, default=str))
    with open(OUT / "slices.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    def dist(key, sel=None):
        v = np.array([r[key] for r in rows if sel is None or sel(r)])
        return {"n": int(v.size), "min": float(v.min()), "p5": float(np.percentile(v, 5)),
                "median": float(np.median(v)), "mean": float(v.mean()), "pct_exact_1": float(100 * np.mean(v == 1.0))}

    small = lambda r: r["area_ours"] < 50
    summary = {
        "scans": len(scans), "selection": dict(Counter(pick.values())),
        "pylidc_db": {"scans": len(db), "patients": len(set(db.values())),
                      "xml_series": len(by), "xml_not_in_db": len(set(by) - set(db)), "db_not_in_xml": len(set(db) - set(by))},
        "two_scan_patients": two,
        "multi_file_series": {s: [str(p.relative_to(XML_ROOT)) for p in ps] for s, ps in by.items() if len(ps) > 1},
        "nodule_count_equal_scans": sum(o["our_nodules"] == o["pylidc_nodules"] for o in scans),
        "partition_equal_scans": sum(o["partition_equal"] for o in scans),
        "reader_sets_equal_scans": sum(o["reader_sets_equal"] for o in scans),
        "annotations": {"ours": sum(o["our_anns"] for o in scans), "pylidc": sum(o["pylidc_anns"] for o in scans),
                        "pylidc_unmatched_contours": sum(o["pylidc_anns_unmatched"] for o in scans),
                        "our_only": sum(o["our_only_anns"] for o in scans),
                        "pylidc_only": sum(o["pylidc_only_anns"] for o in scans),
                        "pylidc_only_reasons": dict(sum((Counter(o["pylidc_only_reasons"]) for o in scans), Counter())),
                        "nodule_id_mismatch": sum(o.get("nodule_id_mismatch", 0) for o in scans),
                        "single_point_rois_in_large_ours": sum(o["single_point_rois_in_large"] for o in scans),
                        "anns_with_single_point_rois": sum(o["anns_with_single_point_rois"] for o in scans)},
        "pylidc_min_tol_hit_scans": [o["series"][-12:] for o in scans if o["pylidc_min_tol_hit"]],
        "pylidc_clusters_gt4": sum(o["pylidc_cluster_gt4"] for o in scans),
        "pylidc_clusters_dup_reader": sum(o["pylidc_clusters_dup_reader"] for o in scans),
        "our_clusters_dup_reader": sum(o["our_clusters_dup_reader"] for o in scans),
        "z": {"max_mismatch_mm": max(o["z_max_mismatch_mm"] for o in scans),
              "z_not_on_slice": sum(o["z_not_on_slice"] for o in scans),
              "scans_dup_slice_z": sum(o["scan_dup_z"] > 0 for o in scans),
              "scans_z_index_offset_mismatch": [o["series"][-12:] for o in scans if o.get("z_index_offset_mismatch")],
              "sop_multi_z": sum(o["sop_multi_z"] for o in scans), "z_multi_sop": sum(o["z_multi_sop"] for o in scans)},
        "exclusion_contours": {"xml": sum(o["exclusion_contours_xml"] for o in scans),
                               "pylidc": sum(o["exclusion_contours_pylidc"] for o in scans)},
        "slice_masks": {
            "iou_vs_pylidc_include_contour_points_True": dist("iou_incl"),
            "iou_vs_pylidc_include_contour_points_False": dist("iou_excl"),
            "iou_ours_minus_outline_vs_pylidc_False": dist("iou_ours_minus_outline_vs_excl"),
            "iou_ours_vs_pylidc_True_plus_outline": dist("iou_ours_vs_incl_plus_outline"),
            "iou_incl_small_lt50px": dist("iou_incl", small),
            "iou_excl_small_lt50px": dist("iou_excl", small),
            "median_area_ours": float(np.median([r["area_ours"] for r in rows])),
            "median_area_ratio_pylTrue_over_ours": float(np.median([r["area_pyl_incl"] / max(1, r["area_ours"]) for r in rows])),
            "median_area_ratio_pylFalse_over_ours": float(np.median([r["area_pyl_excl"] / max(1, r["area_ours"]) for r in rows])),
            "slices_with_exclusion": sum(r["exclusion"] for r in rows),
            "iou_incl_exclusion_slices": dist("iou_incl", lambda r: r["exclusion"]) if any(r["exclusion"] for r in rows) else None,
        },
        "quirks": QUIRKS,
        "scan_mismatches": [{k: o[k] for k in ("series", "xml", "reason", "our_nodules", "pylidc_nodules", "partition_diff",
                                                "pylidc_only_reasons", "our_only_detail")}
                            for o in scans if not o["partition_equal"] or o["our_nodules"] != o["pylidc_nodules"]
                            or o["pylidc_only_anns"] or o["our_only_anns"]],
    }
    (OUT / "SUMMARY.json").write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps({k: v for k, v in summary.items() if k not in ("scan_mismatches",)}, indent=1, default=str))


if __name__ == "__main__":
    main(*sys.argv[1:3])
