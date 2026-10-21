"""Fetch CT pixels for every LIDC series that has reader contours (Phase 0 fetched only a selection).

Per series: download the DICOM zip from TCIA, keep the annotated slices in Hounsfield units (keys "z<mm, 2 decimals>",
as `objdel.data.lidc.ct_slices` reads them), delete the DICOM. Series already in the output directory are skipped,
so a run can be resumed; `--refetch-incomplete` re-fetches series whose sidecar JSON reports fewer slices than
annotated z positions. `build_data` picks the slices up without code changes.

Slice matching (B3-4): a DICOM is kept when its SOPInstanceUID is one the XML's contours cite. A contour z position
that no SOP-matched DICOM covers falls back to the DICOM whose ImagePositionPatient z is within 0.01 mm of it, the
rule the store uses to map contours to slices (`lidc.ct_slices`, P3); a z shared by two DICOMs is ambiguous and
skipped. Two kept DICOMs whose z positions share one archive key ("z<mm, 2 decimals>") would overwrite each
other, and HU values outside int16 or not integral would be wrapped or truncated by the int16 cast: both raise
(I45) instead, and the series is reported FAILED. LIDC-IDRI-0365's XML cites no SOP UID of its TCIA series, and readers 2 and 4 of LIDC-IDRI-0017 cite five
absent UIDs at z positions the other readers' UIDs cover.

    .venv/bin/python -I -m experiments.phase1.fetch_lidc_pixels data/raw/lidc_xml data/lidc_slices [--workers 3]
        [--series UID ...] [--refetch-incomplete]

References:
- TCIA NBIA REST API (getImage by SeriesInstanceUID): https://wiki.cancerimagingarchive.net/display/Public/NBIA+Search+REST+API+Guide
- DICOM ImagePositionPatient (0020,0032), the patient-frame position of the first voxel:
  https://dicom.innolitics.com/ciods/ct-image/image-plane/00200032
- pylidc maps contours to slices by z, not by SOP UID (Scan.slice_zvals, Contour.image_z_position):
  https://pylidc.github.io/contour.html
- DICOM PS3.5 §9.1 UID encoding (digits and dots, <= 64 characters): every series uid, from the XML
  (`objdel.data.lidc.sop_z`) or from --series (`lidc.uid`), is checked before it names a directory or enters the URL
  (review #1: a planted "../.." uid made the rmtree below delete an arbitrary directory):
  https://dicom.nema.org/medical/dicom/current/output/chtml/part05/chapter_9.html
- urllib.parse.quote (percent-encode the uid as one query value): https://docs.python.org/3/library/urllib.parse.html#urllib.parse.quote
"""
import argparse
import json
import os
import shutil
import time
import urllib.parse
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from objdel.data import lidc

API = "https://services.cancerimagingarchive.net/nbia-api/services/v1/getImage?SeriesInstanceUID="
Z_TOL = 0.01  # mm, as objdel.data.lidc.ct_slices


def choose(headers, sops) -> tuple[dict[str, float], int]:
    """({DICOM name: z key} to keep, how many of them matched by z only).

    headers: (name, SOPInstanceUID, ImagePositionPatient z) per DICOM; sops: the XML's {SOP UID: contour z}."""
    keep = {name: sops[uid] for name, uid, _ in headers if uid in sops}
    covered = set(keep.values())
    by_z = 0
    for z in sorted(set(sops.values()) - covered):
        near = [name for name, _, dz in headers if abs(dz - z) < Z_TOL]
        if len(near) == 1 and near[0] not in keep:
            keep[near[0]] = z
            by_z += 1
    return keep, by_z


def z_keys(keep: dict[str, float]) -> dict[str, str]:
    """{DICOM name: archive key}; ValueError when two DICOMs would share a key (one would overwrite the other)."""
    keys = {name: f"z{z:.2f}" for name, z in keep.items()}
    if len(set(keys.values())) < len(keys):
        dup = sorted(k for k in set(keys.values()) if list(keys.values()).count(k) > 1)
        raise ValueError(f"two DICOMs share the z key(s) {dup}")
    return keys


def hounsfield(pixels: np.ndarray, slope: float, intercept: float) -> np.ndarray:
    """Stored values -> int16 HU; ValueError when a value is not an integer or does not fit int16 (I45: the cast
    wrapped and truncated silently)."""
    hu = pixels.astype(np.float64) * slope + intercept
    if not np.array_equal(hu, np.rint(hu)):
        raise ValueError("HU values are not integers (RescaleSlope / RescaleIntercept)")
    if hu.size and (hu.min() < np.iinfo(np.int16).min or hu.max() > np.iinfo(np.int16).max):
        raise ValueError(f"HU values [{hu.min()}, {hu.max()}] do not fit int16")
    return hu.astype(np.int16)


def fetch(series, sops, out_dir: Path, tmp_root: Path) -> int:
    tmp = tmp_root / series
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    zpath, url = tmp / "series.zip", API + urllib.parse.quote(series, safe="")
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=300) as r, open(zpath, "wb") as fh:
                shutil.copyfileobj(r, fh)
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(10 * (attempt + 1))
    import pydicom  # [data] extra; imported here so `--help` works without it

    slices, meta = {}, {}
    with zipfile.ZipFile(zpath) as zf:
        names = [n for n in zf.namelist() if n.lower().endswith(".dcm")]
        headers = []
        for name in names:
            with zf.open(name) as fh:
                ds = pydicom.dcmread(fh, stop_before_pixels=True)
            headers.append((name, str(ds.SOPInstanceUID), float(ds.ImagePositionPatient[2])))
        keep, by_z = choose(headers, sops)
        for name, key in sorted(z_keys(keep).items()):
            with zf.open(name) as fh:
                ds = pydicom.dcmread(fh)
            slices[key] = hounsfield(ds.pixel_array, float(ds.get("RescaleSlope", 1)),
                                     float(ds.get("RescaleIntercept", 0)))
            meta = {"pixel_spacing": [float(v) for v in ds.PixelSpacing],
                    "slice_thickness": float(ds.get("SliceThickness", 0) or 0)}
    (out_dir / f"{series}.json").write_text(json.dumps({
        **meta, "slices": len(slices), "annotated_sops": len(sops), "annotated_z": len(set(sops.values())),
        "matched_by_z": by_z}))
    part = tmp / "slices.npz"
    with open(part, "wb") as fh:
        np.savez_compressed(fh, **slices)
    os.replace(part, out_dir / f"{series}.npz")  # readers never see a partial file
    shutil.rmtree(tmp, ignore_errors=True)
    return len(slices)


def incomplete(out: Path, zmap) -> list[str]:
    """Fetched series with fewer slices than annotated z positions."""
    todo = []
    for s, sops in zmap.items():
        side = out / f"{s}.json"
        if side.exists() and json.loads(side.read_text()).get("slices", 0) < len(set(sops.values())):
            todo.append(s)
    return sorted(todo)


def series_uid(s: str) -> str:
    try:
        return lidc.uid(s)
    except ValueError as e:  # argparse shows an ArgumentTypeError's message, a ValueError's only as "invalid value"
        raise argparse.ArgumentTypeError(str(e)) from None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("xml_root")
    ap.add_argument("out")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--series", nargs="*", default=[], type=series_uid, help="(re-)fetch exactly these series uids")
    ap.add_argument("--refetch-incomplete", action="store_true")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    zmap = lidc.sop_z(a.xml_root)  # ValueError on a series uid that is not a DICOM UID, before any disk or network
    if a.series:  # exactly these
        todo = sorted(a.series)
    else:
        todo = {s for s, sops in zmap.items() if sops and not (out / f"{s}.npz").exists()}
        todo = sorted(todo | set(incomplete(out, zmap) if a.refetch_incomplete else []))

    def one(s):
        t = time.time()
        try:
            return f"{s[-12:]} {fetch(s, zmap[s], out, out.parent / 'lidc_tmp')} slices {time.time() - t:.0f}s"
        except Exception as e:  # keep going; a re-run retries the failures
            return f"{s[-12:]} FAILED {type(e).__name__}: {e}"

    with ThreadPoolExecutor(a.workers) as pool:
        for i, msg in enumerate(pool.map(one, todo), 1):
            print(f"[{i}/{len(todo)}] {msg}", flush=True)


if __name__ == "__main__":
    main()
