"""Fetch LIDC CT slices that carry reader contours, one series at a time.

Selects the series with a reader topology disagreement (Phase 0a, filtered) plus
an equal number of random control series. For each: download the DICOM zip from
TCIA, keep the annotated slices in Hounsfield units, delete the DICOM.

    .venv/bin/python -m phase0.lidc_fetch data/raw/lidc_xml phase0/results/lidc/clusters.csv data/lidc_slices
"""
import csv
import json
import os
import random
import shutil
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pydicom

API = "https://services.cancerimagingarchive.net/nbia-api/services/v1/getImage?SeriesInstanceUID="
NS = {"n": "http://www.nih.gov"}


def select(clusters_csv, seed=0):
    rows = [r for r in csv.DictReader(open(clusters_csv)) if r["sessions"] == "4" and int(r["readers_outlined"]) >= 2]
    flag = lambda r: any(r[k] == "True" for k in ("b0_3d_disagree_5", "hole_disagree_5", "multi_slice_disagree_5", "split_by_a_reader"))
    dis = sorted({r["series"] for r in rows if flag(r)})
    rest = sorted({r["series"] for r in rows} - set(dis))
    ctrl = sorted(random.Random(seed).sample(rest, len(dis)))
    return dis, ctrl


def sop_to_z(xml_root):
    out = {}
    for p in Path(xml_root).rglob("*.xml"):
        root = ET.parse(p).getroot()
        if not root.tag.endswith("LidcReadMessage"):
            continue
        series = root.find("n:ResponseHeader/n:SeriesInstanceUid", NS).text.strip()
        m = out.setdefault(series, {})
        for roi in root.iter("{http://www.nih.gov}roi"):
            m[roi.find("n:imageSOP_UID", NS).text.strip()] = round(float(roi.find("n:imageZposition", NS).text), 2)
    return out


def fetch(series, sops, out_dir, tmp_root):
    tmp = tmp_root / series
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    zpath = tmp / "series.zip"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(API + series, timeout=300) as r, open(zpath, "wb") as fh:
                shutil.copyfileobj(r, fh)
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(10 * (attempt + 1))
    slices, meta = {}, {}
    with zipfile.ZipFile(zpath) as zf:
        for name in zf.namelist():
            if not name.lower().endswith(".dcm"):
                continue
            with zf.open(name) as fh:
                ds = pydicom.dcmread(fh)
            uid = str(ds.SOPInstanceUID)
            if uid not in sops:
                continue
            hu = ds.pixel_array.astype(np.float32) * float(ds.get("RescaleSlope", 1)) + float(ds.get("RescaleIntercept", 0))
            slices[f"z{sops[uid]:.2f}"] = hu.astype(np.int16)
            meta = {"pixel_spacing": [float(v) for v in ds.PixelSpacing], "slice_thickness": float(ds.get("SliceThickness", 0) or 0)}
    (out_dir / f"{series}.json").write_text(json.dumps({**meta, "slices": len(slices), "annotated_sops": len(sops)}))
    part = tmp / "slices.npz"
    with open(part, "wb") as fh:
        np.savez_compressed(fh, **slices)
    os.replace(part, out_dir / f"{series}.npz")  # readers never see a partial file
    shutil.rmtree(tmp, ignore_errors=True)
    return len(slices)


def main(xml_root, clusters_csv, out, workers=3):
    workers = int(workers)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    tmp_root = out.parent / "lidc_tmp"
    dis, ctrl = select(clusters_csv)
    (out / "selection.json").write_text(json.dumps({"disagree": dis, "control": ctrl}, indent=1))
    zmap = sop_to_z(xml_root)
    todo = [s for s in dis + ctrl if not (out / f"{s}.npz").exists()]
    (out / "DONE").unlink(missing_ok=True)

    def one(s):
        t = time.time()
        try:
            n = fetch(s, zmap[s], out, tmp_root)
            return f"{s[-12:]} {n} slices {time.time() - t:.0f}s"
        except Exception as e:
            return f"{s[-12:]} FAILED {type(e).__name__}: {e}"

    with ThreadPoolExecutor(workers) as pool:
        for i, msg in enumerate(pool.map(one, todo), 1):
            print(f"[{i}/{len(todo)}] {msg}", flush=True)
    (out / "DONE").write_text(str(len(dis + ctrl)))


if __name__ == "__main__":
    main(*sys.argv[1:5])
