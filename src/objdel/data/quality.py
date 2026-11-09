"""Data quality checks of the LIDC store (Implementation.md B1.7, gate DG17).

Contour-to-slice check of every series. `lidc.slice_index` puts each contour on the scan slice nearest to its z, so a
series with an ambiguous z can put a contour on the wrong slice. Per series of the store, from the XML and the fetch
sidecars only (no DICOM):
- (a) ``shared_z``: one z cited by two or more SOP UIDs (two images at one position);
- (b) ``far_rois``: a contour z farther than half the slice spacing (median gap of the scan's z values) from every
  slice;
- (d) ``sidecar``: the fetch sidecar has fewer or more kept slices than cited SOP UIDs, or matched slices by z
  (``lidc.ct_slices`` drops a slice key that is not within 0.01 mm of a scan slice, so the contour loses its image).
Every series with a hit goes to step 2 (the DICOM header check, `dicom_z`). Train hits must be dropped
(`train_exclusions.json`) or mapped (fixed, store rebuilt, re-run lists no hit); val and test hits keep their items
in the frozen views and are listed as known label issues.

References:
- LIDC-IDRI known issues ("Annotation and Markup Issues" item 6: LIDC-IDRI-0396, 139.xml, SOP Instance UID at
  position 1420, fixed 2018-06-28): https://www.cancerimagingarchive.net/collection/lidc-idri/
- Third-party list of duplicate-z series (0085, 0146, 0418, 0572, 0979):
  https://huggingface.co/datasets/YongchengYAO/LIDC-IDRI-Lite
- DICOM ImagePositionPatient (0020,0032): https://dicom.innolitics.com/ciods/ct-image/image-plane/00200032
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from . import lidc
from .store import SPLITS

KNOWN_0396 = ("LIDC-IDRI-0396", 1420.0)  # TCIA known issue 6: patient, contour z of the corrected SOP UID


def file_sha256(path) -> str:
    with open(path, "rb") as fh:
        return hashlib.file_digest(fh, "sha256").hexdigest()


def _sidecar_hit(path: Path) -> dict | None:
    if not path.exists():
        return None
    j = json.loads(path.read_text())
    s = {"slices": j["slices"], "annotated_sops": j["annotated_sops"], "matched_by_z": j.get("matched_by_z", 0)}
    return s if s["slices"] != s["annotated_sops"] or s["matched_by_z"] > 0 else None


def z_report(store_root, xml_root, sidecars, scans=None, previous=None) -> dict:
    """One entry per LIDC series of the store at ``store_root`` (split, patient, item count, hits (a), (b), (d)).

    ``step2``: the series with a hit; ``known_label_issues``: the item ids of val and test hit series; ``mapped``: the
    train hit series of ``previous`` (an older report) that list no hit now; ``check_0396``: the SOP UIDs that
    LIDC-IDRI-0396's XML cites at z 1420.0 and whether the fetch found every cited SOP UID (by UID, none by z)."""
    store_root, sidecars = Path(store_root), Path(sidecars)
    scans = lidc.scans() if scans is None else scans
    sop_z = lidc.sop_z(xml_root)
    series, item_ids = {}, {}
    for split in SPLITS:
        for row in pq.read_table(store_root / "lidc" / f"{split}.parquet", columns=["item_id", "group_id"]).to_pylist():
            s = row["item_id"].split("/")[1]
            series.setdefault(s, {"split": split, "patient": row["group_id"], "items": 0})["items"] += 1
            item_ids.setdefault(s, []).append(row["item_id"])
    for s, e in series.items():
        by_z = {}
        for sop, z in sop_z.get(s, {}).items():
            by_z.setdefault(z, []).append(sop)
        e["shared_z"] = {str(z): sorted(sops) for z, sops in sorted(by_z.items()) if len(sops) > 1}
        zvals = scans[s].zvals
        half = float(np.median(np.diff(zvals))) / 2
        e["far_rois"] = [{"sop": sop, "z": z, "distance_mm": round(d, 2), "spacing_mm": round(2 * half, 2)}
                         for sop, z in sorted(sop_z.get(s, {}).items(), key=lambda kv: (kv[1], kv[0]))
                         if (d := float(np.abs(zvals - z).min())) > half]
        e["sidecar"] = _sidecar_hit(sidecars / f"{s}.json")
    hits = sorted(s for s, e in series.items() if e["shared_z"] or e["far_rois"] or e["sidecar"])
    for s in hits:
        series[s]["item_ids"] = sorted(item_ids[s])
    old = [s for s in (previous or {}).get("step2", []) if previous["series"][s]["split"] == "train"]
    patient, z0396 = KNOWN_0396
    c = next((s for s, sc in scans.items() if sc.patient == patient), None)
    side = json.loads((sidecars / f"{c}.json").read_text()) if c else {}
    return {"train_sha256": file_sha256(store_root / "lidc" / "train.parquet"), "series": series, "step2": hits,
            "mapped": sorted(s for s in old if s not in hits),
            "known_label_issues": {sp: sorted(i for s in hits if series[s]["split"] == sp for i in item_ids[s])
                                   for sp in ("val", "test")},
            "check_0396": c and {"series": c,
                                 "sop_uids_at_1420": sorted(u for u, z in sop_z.get(c, {}).items() if z == z0396),
                                 "all_cited_sops_fetched": side.get("slices") == side.get("annotated_sops")
                                 and side.get("matched_by_z", 0) == 0}}


def train_hits_not_excluded(report, exclusions) -> dict[str, list[str]]:
    """Train hit series -> its item ids that are not keys of ``exclusions`` (series with none left out)."""
    out = {}
    for s in report["step2"]:
        if report["series"][s]["split"] == "train" and (m := [i for i in report["series"][s]["item_ids"]
                                                              if i not in exclusions]):
            out[s] = m
    return out


def dg17(report, exclusions, train_sha256) -> list[str]:
    """Gate DG17: failure lines, [] passes. A report on another train store is stale; a train hit is resolved only
    when every store item of its series is a key of ``exclusions`` (a mapped series lists no hit any more). Val and
    test hits never fail it (they are known label issues)."""
    if report["train_sha256"] != train_sha256:
        return [f"stale report: train store {train_sha256}, report {report['train_sha256']}"]
    return [f"{s}: {len(m)} of {report['series'][s]['items']} items not in train_exclusions.json"
            for s, m in train_hits_not_excluded(report, exclusions).items()]


def dicom_z(folder, zvals) -> dict:
    """Step 2 on the DICOM headers in ``folder`` (``*.dcm``, untrusted: run under ``python -I``): z positions
    (ImagePositionPatient[2], 0.01 mm) shared by two or more SOP Instance UIDs, and the scan z values that
    `lidc.slice_index` maps contours to but no DICOM has."""
    import pydicom  # [data] extra

    by_z = {}
    for p in sorted(Path(folder).glob("*.dcm")):
        ds = pydicom.dcmread(p, stop_before_pixels=True)
        by_z.setdefault(round(float(ds.ImagePositionPatient[2]), 2), []).append(str(ds.SOPInstanceUID))
    return {"duplicate_z": {str(z): sorted(u) for z, u in sorted(by_z.items()) if len(u) > 1},
            "slice_z_without_dicom": [float(z) for z in zvals if not any(abs(z - d) < 0.01 for d in by_z)]}
