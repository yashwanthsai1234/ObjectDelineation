"""B1.7: contour-to-slice check of every LIDC series of the store, and gate DG17 (Implementation.md B1.7).

Writes a JSON report with one entry per series (`objdel.data.quality.z_report`), prints every hit, and exits 1 when
DG17 fails (a train hit whose series is not fully in the exclusion list, or a report on another train store).
``--drop-train-hits`` adds every store item of each train hit series to the exclusion list (reason "DG17 <series
suffix>"); nothing else removes an item. ``--previous`` names the report before a rebuild: its train hits that list no
hit now are reported under ``mapped``, and a changed val or test store stops the run (the frozen views pair on those
items). ``--dicom DIR`` runs step 2 on ``DIR/<series>/*.dcm`` for every hit series that has a folder; without it no
DICOM is read. Run with ``python -I`` (the XML, sidecars and DICOM are downloaded data):

    .venv/bin/python -I -m experiments.phase2.check_series_z [--drop-train-hits] [--previous OLD] [--dicom DIR]
"""
import argparse
import datetime
import json
from pathlib import Path

from objdel.core.io import read_json
from objdel.data import lidc, quality


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--store", type=Path, default=Path("data/store"))
    ap.add_argument("--xml", type=Path, default=Path("data/raw/lidc_xml"))
    ap.add_argument("--slices", type=Path, default=Path("data/lidc_slices"), help="fetch sidecars <series>.json")
    ap.add_argument("--exclusions", type=Path, default=Path("reference/train_exclusions.json"))
    ap.add_argument("--out", type=Path, default=Path("runs/b1/series_z.json"))
    ap.add_argument("--previous", type=Path, help="the report before a rebuild of the store")
    ap.add_argument("--dicom", type=Path, help="step 2: folder of <series>/*.dcm headers")
    ap.add_argument("--drop-train-hits", action="store_true")
    a = ap.parse_args(argv)
    previous = read_json(a.previous) if a.previous else None
    scans = lidc.scans()
    r = quality.z_report(a.store, a.xml, a.slices, scans, previous)
    r["store_sha256"] = {sp: quality.file_sha256(a.store / "lidc" / f"{sp}.parquet") for sp in ("val", "test")}
    if previous and previous["store_sha256"] != r["store_sha256"]:
        raise SystemExit("the val or test store changed since --previous: the frozen views pair on those items")
    exclusions = read_json(a.exclusions) if a.exclusions.exists() else {}
    if a.drop_train_hits:
        today = datetime.date.today().isoformat()
        for s in quality.train_hits_not_excluded(r, exclusions):
            for i in r["series"][s]["item_ids"]:
                exclusions.setdefault(i, {"reason": f"DG17 {s[-6:]}", "date": today})
        a.exclusions.write_text(json.dumps(dict(sorted(exclusions.items())), indent=1) + "\n")
    open_hits = quality.train_hits_not_excluded(r, exclusions)
    for s in r["step2"]:
        e = r["series"][s]
        if e["split"] == "train":
            e["resolution"] = None if s in open_hits else "dropped"
        if a.dicom and (a.dicom / s).is_dir():
            e["dicom"] = quality.dicom_z(a.dicom / s, scans[s].zvals)
        print(f"{e['split']:5} {e['patient']} {s[-6:]} items {e['items']}: shared_z {len(e['shared_z'])}, "
              f"far_rois {len(e['far_rois'])}, sidecar {e['sidecar']}, resolution {e.get('resolution', 'keep')}")
    print(f"{len(r['series'])} series, {len(r['step2'])} with a hit; mapped {r['mapped']}; "
          f"known label issues: val {len(r['known_label_issues']['val'])}, test {len(r['known_label_issues']['test'])}"
          f" items; 0396: {r['check_0396']}")
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(r, indent=1) + "\n")
    fails = quality.dg17(r, exclusions, r["train_sha256"])
    print("DG17 pass" if not fails else "DG17 fail:\n  " + "\n  ".join(fails))
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
