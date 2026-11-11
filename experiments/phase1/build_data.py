"""Build the LIDC item store, splits and K11 checks into data/store/.

    .venv/bin/python -I -m experiments.phase1.build_data [--out data/store] [--workers 2] [--seed 0]

Reads data/raw/lidc_xml, data/lidc_slices (CT pixels where fetched) and pylidc's bundled database.
Writes <out>/lidc/{train,val,test}.parquet, manifest.json (SHA-256 per file) and build_report.json (counts, disagreement
counts, K11 results, timings). LIDC copied outlines count once, and items whose model frame
(objdel.data.targets.frame) overlaps another outlined nodule on the slice are not stored.

References:
- The per-dataset rules and their sources: the docstrings of objdel.data.lidc (LIDC-IDRI, Armato et al. 2011,
  https://doi.org/10.1118/1.3528204), objdel.data.splits and objdel.data.store.
- multiprocessing.Pool.map (chunksize; results in input order): https://docs.python.org/3/library/multiprocessing.html
"""
import argparse
import csv
import hashlib
import json
import shutil
import time
from collections import Counter
from multiprocessing import get_context
from pathlib import Path

from objdel.data import lidc, splits, store
from objdel.data.store import Item
from objdel.data.targets.frame import frame
from objdel.geom.topology import min_area, signature

REPO = Path(__file__).resolve().parents[2]
LAB_SPLIT = REPO / "data" / "fixtures" / "lidc_full_meta.csv"  # the lab's patient split, for the plan's 313 / 251


def crowded(s):
    """The item's model frame shows another outlined nodule of the slice."""
    it = Item("", "lidc", "", "train", None, {"nodule": s.full()}, [str(r) for r in s.readers], {})
    return lidc.crowded(s, frame(it).box)


def slice_flags(s):
    """Report counts of one slice: readers' per-instance topologies differ under the 5-px filter / the plan's size
    rule (with and without merged-reader clusters), copied layers, and a rater whose size-filtered signature
    is () and whether that rater alone makes the disagreement."""
    rule = {signature(m, min_area(int(s.stack.any(0).sum()))) for m in s.stack}
    return {"disagree_5px": len({signature(m, 5) for m in s.stack}) > 1, "disagree_rule": len(rule) > 1,
            "disagree_rule_no_dup_reader": len(rule) > 1 and not s.dup_reader,
            "with_copied_raters": bool(s.copied), "copied_layers": len(s.copied), "sub_filter_rater": () in rule,
            "disagree_only_from_sub_filter_rater": len(rule - {()}) == 1 < len(rule)}


def lidc_series(job):
    """(series, sessions, stored slices, their flags, slices left out because their frame is crowded)."""
    series, path, zvals = job
    _, sessions = lidc.parse(path)
    sl = lidc.slices(sessions, zvals)
    kept = [s for s in sl if not crowded(s)]
    return series, len(sessions), kept, [slice_flags(s) for s in kept], len(sl) - len(kept)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "data" / "store"))
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    out, raw = Path(a.out), REPO / "data" / "raw"
    t0 = time.time()
    rep = {"seed": a.seed, "timing_s": {}}

    scans = lidc.scans()
    files = lidc.series_files(raw / "lidc_xml")
    assert set(files) == set(scans) and len(scans) == 1018
    with get_context("spawn").Pool(a.workers) as pool:
        # Pass 1: what the splits need: LIDC slices and their disagreement.
        lidc_res = pool.map(lidc_series, [(s, files[s], scans[s].zvals) for s in sorted(files)], chunksize=8)
        per_patient = {}
        for series, _, sl, fl, _ in lidc_res:
            n, d = per_patient.get(scans[series].patient, (0, 0))
            per_patient[scans[series].patient] = (n + len(sl), d + sum(f["disagree_5px"] for f in fl))
        patients = {s.patient: splits.Patient(*per_patient.get(s.patient, (0, 0))) for s in scans.values()}
        lidc_split = splits.lidc_split(patients, a.seed)
        rep["timing_s"]["pass1"] = round(time.time() - t0, 1)

        # Pass 2: stream every item into the store.
        stats = {"image_sha": {}, "split_of": {}, "lidc": Counter()}

        def items():
            for series, n_sessions, sl, fl, _ in lidc_res:
                scan = scans[series]
                npz = REPO / "data" / "lidc_slices" / f"{series}.npz"
                ct = lidc.ct_slices(npz, scan.zvals) if npz.exists() else {}
                split = lidc_split[scan.patient]
                for s, f in zip(sl, fl, strict=True):
                    it = lidc.item(series, scan, s, split, ct.get(s.k), n_sessions,
                                   str(files[series].relative_to(raw / "lidc_xml")))
                    stats["lidc"].update({f"{split}/{k}": v for k, v in
                                          {"items": 1, "with_image": it.image is not None, **f}.items()})
                    stats["split_of"][it.item_id] = split
                    if it.image is not None:
                        stats["image_sha"][it.item_id] = hashlib.sha256(it.image.tobytes()).hexdigest()
                    yield it

        if out.exists():
            shutil.rmtree(out)
        manifest = store.write(items(), out)
    rep["timing_s"]["total"] = round(time.time() - t0, 1)

    # Counts and K11.
    lab_split = {}
    if LAB_SPLIT.exists():
        lab_split = {r["patient"]: r["split"] for r in csv.DictReader(open(LAB_SPLIT))}
    lab_train = [f for series, _, _, fl, _ in lidc_res for f in fl if lab_split.get(scans[series].patient) == "train"]
    crowded_per_split = Counter()
    for series, _, _, _, n in lidc_res:
        crowded_per_split[lidc_split[scans[series].patient]] += n
    patient_series = {s: scans[s].patient for s in scans}
    rep.update({
        "manifest": manifest,
        "lidc": {"series": len(lidc_res), "patients": len(patients), "items": sum(len(sl) for _, _, sl, _,
            _ in lidc_res),
                 "dropped_crowded_frame": dict(sorted(crowded_per_split.items())),
                 "patients_per_split": dict(sorted(Counter(lidc_split.values()).items())),
                 **{k: stats["lidc"][k] for k in sorted(stats["lidc"])},
                 **{f"lab_split_train_{k}": sum(f[k] for f in lab_train)
                    for k in ("disagree_5px", "disagree_rule", "disagree_rule_no_dup_reader")}},
        "k11": {"lidc_patient_leaks": splits.leaks({s: lidc_split[p] for s, p in patient_series.items()},
                                                   patient_series),
                "lidc_identical_slices_across_splits": splits.leaks(stats["split_of"], stats["image_sha"])},
    })
    (out / "build_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in rep.items() if k != "manifest"}, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
