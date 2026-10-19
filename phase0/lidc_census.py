"""Phase 0 census on LIDC-IDRI radiologist annotations (XML only, no CT needed).

Each CT series has up to 4 reading sessions (final, unblinded reads). Nodules
>= 3 mm are outlined slice by slice; inclusion=FALSE contours mark regions to
exclude (holes); nodules < 3 mm are a single point. Contours are pixel lists on
the 512 x 512 slice grid.

Raster convention: an inclusion contour covers its filled polygon including the
outline pixels; an exclusion removes its polygon interior (outline pixels stay).

    .venv/bin/python -m phase0.lidc_census data/raw/lidc_xml phase0/results/lidc
"""
import csv
import json
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from checker import betti, betti3d, certify, size_filter

NS = {"n": "http://www.nih.gov"}
GRID = 512
MIN_SIZE_2D = 5


def text(el, tag):
    t = el.find(f"n:{tag}", NS)
    return t.text.strip() if t is not None and t.text else None


def draw(points, outline_only=False):
    img = Image.new("1", (GRID, GRID), 0)
    d = ImageDraw.Draw(img)
    if len(points) == 1:
        d.point(points, fill=1)
    elif outline_only:
        d.line(points + [points[0]], fill=1)
    else:
        d.polygon(points, fill=1, outline=1)
    return np.asarray(img, dtype=bool)


def parse(path):
    root = ET.parse(path).getroot()
    if not root.tag.endswith("LidcReadMessage"):
        return None
    series = text(root.find("n:ResponseHeader", NS), "SeriesInstanceUid")
    sessions = []
    for s in root.findall("n:readingSession", NS):
        nodules = []
        for nod in s.findall("n:unblindedReadNodule", NS):
            rois = []
            for roi in nod.findall("n:roi", NS):
                pts = [(int(e.find("n:xCoord", NS).text), int(e.find("n:yCoord", NS).text))
                       for e in roi.findall("n:edgeMap", NS)]
                if pts:
                    rois.append((float(text(roi, "imageZposition")), text(roi, "inclusion") == "TRUE", pts))
            if rois:
                large = nod.find("n:characteristics", NS) is not None and any(len(p) > 1 for _, _, p in rois)
                nodules.append({"id": text(nod, "noduleID"), "large": large, "rois": rois})
        sessions.append(nodules)
    return series, sessions


def z_index(sessions):
    zs = sorted({round(z, 2) for nods in sessions for n in nods for z, _, _ in n["rois"]})
    if len(zs) < 2:
        return {z: 0 for z in zs}
    diffs = np.diff(zs)
    step = float(np.min(diffs[diffs > 0.05])) if np.any(diffs > 0.05) else 1.0
    return {z: int(round((z - zs[0]) / step)) for z in zs}


def rasterize(nodule, zi):
    """{z index: 2D mask} for one reader's nodule."""
    inc, exc = defaultdict(lambda: np.zeros((GRID, GRID), bool)), defaultdict(list)
    for z, inclusion, pts in nodule["rois"]:
        k = zi[round(z, 2)]
        if inclusion:
            inc[k] |= draw(pts)
        else:
            exc[k].append(pts)
    n_exc = 0
    for k, polys in exc.items():
        for pts in polys:
            if len(pts) >= 3:
                interior = draw(pts) & ~draw(pts, outline_only=True)
                n_exc += int(interior.any())
                inc[k] &= ~interior
    return dict(inc), n_exc


def voxels(slices):
    out = []
    for k, m in slices.items():
        y, x = np.nonzero(m)
        out.append(np.stack([np.full(y.size, k), y, x], axis=1))
    return np.concatenate(out) if out else np.zeros((0, 3), int)


def to_volume(slices, box):
    z0, y0, x0, z1, y1, x1 = box
    v = np.zeros((z1 - z0 + 1, y1 - y0 + 1, x1 - x0 + 1), bool)
    for k, m in slices.items():
        v[k - z0] |= m[y0:y1 + 1, x0:x1 + 1]
    return v


def nodule_clusters(path):
    """Readers' nodule outlines grouped into clusters across readers.

    Returns None for non-CT files, else (series, sessions, z index map, clusters, stats).
    Each cluster: members (annotations), box (z0, y0, x0, z1, y1, x1), readers,
    small_readers, and per_reader {reader: {z index: 2D mask}}.
    """
    parsed = parse(path)
    if parsed is None:
        return None
    series, sessions = parsed
    zi = z_index(sessions)
    anns = []  # (reader, nodule, slices, voxel keys)
    smalls = []  # (reader, z index, x, y)
    stats = {"cert_fail": 0, "slices_certified": 0, "exclusions_drawn": 0}
    for r, nods in enumerate(sessions):
        for nod in nods:
            if nod["large"]:
                sl, n_exc = rasterize(nod, zi)
                sl = {k: m for k, m in sl.items() if m.any()}
                if not sl:
                    continue
                stats["exclusions_drawn"] += n_exc
                for m in sl.values():
                    c = certify(m)
                    stats["slices_certified"] += 1
                    stats["cert_fail"] += not c.certified
                vox = voxels(sl)
                keys = set(map(tuple, vox.tolist()))
                anns.append({"reader": r, "slices": sl, "keys": keys, "exclusions": n_exc})
            else:
                z, _, pts = nod["rois"][0]
                smalls.append((r, zi[round(z, 2)], pts[0][0], pts[0][1]))

    # Cluster large annotations across readers by shared voxels (union-find).
    parent = list(range(len(anns)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    owner = {}
    for i, a in enumerate(anns):
        for key in a["keys"]:
            j = owner.setdefault(key, i)
            if j != i:
                parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i in range(len(anns)):
        groups[find(i)].append(i)

    clusters = []
    for idx in groups.values():
        members = [anns[i] for i in idx]
        keys = set().union(*(a["keys"] for a in members))
        arr = np.array(sorted(keys))
        box = (*arr.min(axis=0), *arr.max(axis=0))
        readers = sorted({a["reader"] for a in members})
        small_readers = set()
        for r, k, x, y in smalls:
            if r not in readers and box[0] - 1 <= k <= box[3] + 1 and box[1] - 3 <= y <= box[4] + 3 and box[2] - 3 <= x <= box[5] + 3:
                small_readers.add(r)
        per_reader = {}
        for r in readers:
            sl = defaultdict(lambda: np.zeros((GRID, GRID), bool))
            for a in members:
                if a["reader"] == r:
                    for k, m in a["slices"].items():
                        sl[k] |= m
            per_reader[r] = dict(sl)
        clusters.append({"members": members, "box": box, "readers": readers,
                         "small_readers": small_readers, "per_reader": per_reader})
    return series, sessions, zi, clusters, stats


def analyse(path):
    got = nodule_clusters(path)
    if got is None:
        return None
    series, sessions, _, found, stats = got
    clusters = []
    for cl in found:
        members, box, readers, small_readers = cl["members"], cl["box"], cl["readers"], cl["small_readers"]
        per_reader = []
        for r in readers:
            mine = [a for a in members if a["reader"] == r]
            sl = cl["per_reader"][r]
            row = {"reader": r, "annotations": len(mine), "exclusions": sum(a["exclusions"] for a in mine),
                   "slices": len(sl), "voxels": int(sum(m.sum() for m in sl.values()))}
            for ms in (0, MIN_SIZE_2D):
                fs = {k: (size_filter(m, ms) if ms else m) for k, m in sl.items()}
                fs = {k: m for k, m in fs.items() if m.any()}
                topo = {k: betti(m) for k, m in fs.items()}
                row[f"b0_3d_{ms}"], row[f"b2_3d_{ms}"] = betti3d(to_volume(fs, box)) if fs else (0, 0)
                row[f"slices_multi_{ms}"] = sum(1 for b0, _ in topo.values() if b0 > 1)
                row[f"slices_hole_{ms}"] = sum(1 for _, b1 in topo.values() if b1 > 0)
                row[f"slice_topo_{ms}"] = topo
            per_reader.append(row)
        # Slice-level comparison where at least two readers drew the slice.
        slice_cmp = {}
        for ms in (0, MIN_SIZE_2D):
            both = disagree = 0
            for k in {k for p in per_reader for k in p[f"slice_topo_{ms}"]}:
                ts = [p[f"slice_topo_{ms}"][k] for p in per_reader if k in p[f"slice_topo_{ms}"]]
                if len(ts) >= 2:
                    both += 1
                    disagree += len(set(ts)) > 1
            slice_cmp[ms] = (both, disagree)
        c = {"series": series, "file": str(path.name), "sessions": len(sessions),
             "readers_outlined": len(readers), "readers_small_only": len(small_readers),
             "readers_none": len(sessions) - len(readers) - len(small_readers),
             "split_by_a_reader": any(p["annotations"] > 1 for p in per_reader),
             "max_voxels": max(p["voxels"] for p in per_reader)}
        for ms in (0, MIN_SIZE_2D):
            b0s = [p[f"b0_3d_{ms}"] for p in per_reader]
            holes = [p[f"slices_hole_{ms}"] > 0 for p in per_reader]
            multi = [p[f"slices_multi_{ms}"] > 0 for p in per_reader]
            c[f"b0_3d_set_{ms}"] = ";".join(map(str, b0s))
            c[f"b0_3d_disagree_{ms}"] = len(set(b0s)) > 1
            c[f"hole_any_{ms}"] = any(holes)
            c[f"hole_disagree_{ms}"] = len(set(holes)) > 1
            c[f"multi_slice_disagree_{ms}"] = len(set(multi)) > 1
            c[f"cavity_any_{ms}"] = any(p[f"b2_3d_{ms}"] > 0 for p in per_reader)
            c[f"slices_compared_{ms}"], c[f"slices_disagree_{ms}"] = slice_cmp[ms]
        c["exclusions"] = sum(p["exclusions"] for p in per_reader)
        clusters.append(c)
    return {"series": series, "file": str(path), "clusters": clusters, **stats}


def main(root, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    files = sorted(Path(root).rglob("*.xml"))
    with Pool() as pool:
        results = [r for r in pool.map(analyse, files, chunksize=4) if r]
    # A top-level "resubmitted-correction" file replaces the reading of the same series.
    by_series = {}
    for r in results:
        if r["series"] not in by_series or "correction" in Path(r["file"]).name:
            by_series[r["series"]] = r
    results = list(by_series.values())
    clusters = [c for r in results for c in r["clusters"]]
    keys = list(dict.fromkeys(k for c in clusters for k in c))
    with open(out / "clusters.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(clusters)

    def share(sel, key):
        return {"n": sum(1 for c in sel if c[key]), "of": len(sel),
                "pct": round(100 * sum(1 for c in sel if c[key]) / max(1, len(sel)), 2)}

    four = [c for c in clusters if c["sessions"] == 4]
    multi = [c for c in four if c["readers_outlined"] >= 2]
    s = {"series": len(results), "clusters": len(clusters),
         "slices_certified": sum(r["slices_certified"] for r in results),
         "certificate_failures": sum(r["cert_fail"] for r in results),
         "exclusion_interiors_drawn": sum(r["exclusions_drawn"] for r in results),
         "clusters_4_sessions": len(four),
         "readers_outlined_hist": dict(sorted(Counter(c["readers_outlined"] for c in four).items())),
         "presence_disagree": {"n": sum(1 for c in four if c["readers_outlined"] < 4), "of": len(four)},
         "outlined_by_all_4": sum(1 for c in four if c["readers_outlined"] == 4),
         "multi_reader_clusters": len(multi),
         "split_by_a_reader": share(multi, "split_by_a_reader")}
    for ms in (0, MIN_SIZE_2D):
        s[f"min_size_{ms}"] = {
            "b0_3d_disagree": share(multi, f"b0_3d_disagree_{ms}"),
            "hole_any": share(multi, f"hole_any_{ms}"),
            "hole_disagree": share(multi, f"hole_disagree_{ms}"),
            "multi_slice_disagree": share(multi, f"multi_slice_disagree_{ms}"),
            "cavity_any": share(multi, f"cavity_any_{ms}"),
            "any_topology_disagree": {
                "n": sum(1 for c in multi if c[f"b0_3d_disagree_{ms}"] or c[f"hole_disagree_{ms}"]
                         or c[f"multi_slice_disagree_{ms}"] or c["split_by_a_reader"]),
                "of": len(multi)},
            "slices_compared": sum(c[f"slices_compared_{ms}"] for c in multi),
            "slices_topology_disagree": sum(c[f"slices_disagree_{ms}"] for c in multi),
        }
    (out / "summary.json").write_text(json.dumps(s, indent=2))
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:3])
