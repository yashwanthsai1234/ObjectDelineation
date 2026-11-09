"""Fixed LIDC XML parser + clustering (sandbox copy; based on experiments.phase0.lidc_census).

Fixes over Phase 0:
  P1  Single-point ROIs inside >= 3 mm nodules are dropped (pylidc's DB drops
      them too); they rasterized as stray 1-px components.
  P2  ONE outline rule, applied to every contour: a contour's pixels belong to
      the closed polygon it draws (PIL polygon fill with outline). Inclusion
      contours add their polygon (interior + outline); exclusion contours
      subtract their polygon (interior + outline). All inclusions of a slice
      are applied before its exclusions (as pylidc does).
      Inclusion side = Phase 0 behaviour (bit-identical on all slices without
      exclusions), so Phase 0 numbers and fixtures stay comparable. Exclusion
      side changes from Phase 0 (which kept the exclusion outline): treating
      both contour types alike is the only symmetric choice, and it equals
      pylidc's exclusion handling. Known bias: LIDC readers drew outlines just
      outside the nodule, so masks are ~1 px "generous"; it matters most for
      nodules < 50 px.
  P3  True slice index k from the scan's slice z positions (pylidc DB,
      Scan.slice_zvals), not from the spacing of annotated z values.
  P4  Patient id from series_patient.csv (XMLs carry none).
Clustering: union-find over shared voxels (z index, y, x), as in Phase 0.
"""
import csv
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

NS = {"n": "http://www.nih.gov"}
GRID = 512


def text(el, tag):
    t = el.find(f"n:{tag}", NS)
    return t.text.strip() if t is not None and t.text else None


def draw(points):
    img = Image.new("1", (GRID, GRID), 0)
    d = ImageDraw.Draw(img)
    if len(points) == 1:
        d.point(points, fill=1)
    else:
        d.polygon(points, fill=1, outline=1)
    return np.asarray(img, dtype=bool)


def parse(path):
    """-> (series, sessions, stats). sessions: list of reader sessions, each a list of nodules
    {id, large, rois: [(z, inclusion, pts)]}. Single-point ROIs of large nodules are dropped (P1)."""
    root = ET.parse(path).getroot()
    if not root.tag.endswith("LidcReadMessage"):
        return None
    series = text(root.find("n:ResponseHeader", NS), "SeriesInstanceUid")
    sessions, dropped = [], 0
    for s in root.findall("n:readingSession", NS):
        nodules = []
        for nod in s.findall("n:unblindedReadNodule", NS):
            rois = []
            for roi in nod.findall("n:roi", NS):
                pts = [(int(e.find("n:xCoord", NS).text), int(e.find("n:yCoord", NS).text))
                       for e in roi.findall("n:edgeMap", NS)]
                if pts:
                    rois.append((float(text(roi, "imageZposition")), text(roi, "inclusion") == "TRUE", pts))
            if not rois:
                continue
            large = nod.find("n:characteristics", NS) is not None and any(len(p) > 1 for _, _, p in rois)
            if large:
                n0 = len(rois)
                rois = [r for r in rois if len(r[2]) > 1]
                dropped += n0 - len(rois)
            nodules.append({"id": text(nod, "noduleID"), "large": large, "rois": rois})
        sessions.append(nodules)
    return series, sessions, {"single_point_rois_dropped": dropped}


def slice_index(zvals):
    zv = np.asarray(sorted(zvals), float)

    def k_of(z):
        k = int(np.abs(zv - z).argmin())
        return k, float(abs(zv[k] - z))
    return k_of


def rasterize(nodule, k_of):
    """{k: mask} for one annotation under rule P2; also max |z - slice z| and exclusion count."""
    inc, exc = defaultdict(lambda: np.zeros((GRID, GRID), bool)), defaultdict(list)
    dz = 0.0
    for z, inclusion, pts in nodule["rois"]:
        k, d = k_of(z)
        dz = max(dz, d)
        if inclusion:
            inc[k] |= draw(pts)
        else:
            exc[k].append(pts)
    n_exc = 0
    for k, polys in exc.items():
        for pts in polys:
            if k in inc:
                inc[k] &= ~draw(pts)
                n_exc += 1
    return {k: m for k, m in inc.items() if m.any()}, dz, n_exc


def clusters(path, zvals):
    """Large-nodule annotations clustered across readers. Returns dict or None."""
    got = parse(path)
    if got is None:
        return None
    series, sessions, stats = got
    k_of = slice_index(zvals)
    anns = []
    stats.update({"max_dz": 0.0, "exclusions": 0, "empty_after_raster": 0})
    for r, nods in enumerate(sessions):
        for nod in nods:
            if not nod["large"]:
                continue
            sl, dz, ne = rasterize(nod, k_of)
            stats["max_dz"] = max(stats["max_dz"], dz)
            stats["exclusions"] += ne
            if not sl:
                stats["empty_after_raster"] += 1
                continue
            keys = set()
            for k, m in sl.items():
                y, x = np.nonzero(m)
                keys.update(zip([k] * y.size, y.tolist(), x.tolist()))
            anns.append({"reader": r, "nid": nod["id"], "rois": nod["rois"], "slices": sl, "keys": keys})
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
    out = []
    for idx in sorted(groups.values(), key=lambda g: min(min(k for k, _, _ in anns[i]["keys"]) for i in g)):
        members = [anns[i] for i in idx]
        per_reader = defaultdict(dict)
        for a in members:
            for k, m in a["slices"].items():
                per_reader[a["reader"]][k] = per_reader[a["reader"]].get(k, np.zeros((GRID, GRID), bool)) | m
        out.append({"member_ids": idx, "readers": sorted(per_reader), "per_reader": dict(per_reader),
                    "dup_reader": len(members) > len(per_reader)})
    return {"series": series, "sessions": len(sessions), "anns": anns, "clusters": out, "stats": stats}


def load_series_patient(path):
    return {r["series_uid"]: r for r in csv.DictReader(open(path))}


def chosen_xml(sp_row, xml_root):
    return Path(xml_root) / sp_row["xml_file"]
