"""M4: build the item store (RIGA + LIDC), splits and K11 checks into data/store/.

    .venv/bin/python -I -m experiments.phase1.build_data [--out data/store] [--workers 2] [--seed 0]

Reads data/raw/riga, data/raw/lidc_xml, data/lidc_slices (CT pixels where fetched) and pylidc's bundled database.
Writes <out>/{riga,lidc}/{train,val,test}.parquet, manifest.json (SHA-256 per file) and build_report.json
(counts, stage histograms, disagreement counts, K11 results, timings). The RIGA licence forbids redistribution:
the store stays local.

Wave-6 fixes (experiments/phase1/impl/fixes/F4.md): RIGA items with a reconstructed reference have no image (B3-1);
LIDC copied outlines count once (B3-3) and items whose model frame (objdel.targets.frame) overlaps another
outlined nodule on the slice are not stored (B3-2); RIGA near-duplicates also by ORB + RANSAC inliers (B3-8).
Review fixes (experiments/phase1/impl/fixes/R1.md): the report counts the images the store gets (review #4).

References:
- The per-dataset rules and their sources: the docstrings of objdel.data.riga (RIGA, Almazroa et al. 2018,
  https://doi.org/10.1117/12.2293584), objdel.data.lidc (LIDC-IDRI, Armato et al. 2011,
  https://doi.org/10.1118/1.3528204), objdel.data.splits and objdel.data.store.
- multiprocessing.Pool.map (chunksize; results in input order): https://docs.python.org/3/library/multiprocessing.html
- Design notes and K11 results: experiments/phase1/impl/M4_data.md, experiments/phase1/impl/fixes/F4.md.
"""
import argparse
import csv
import hashlib
import json
import shutil
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from objdel.data import Item, lidc, riga, splits, store
from objdel.geom.topology import min_area, signature
from objdel.targets import frame

REPO = Path(__file__).resolve().parents[2]
LAB_SPLIT = REPO / "data" / "fixtures" / "lidc_full_meta.csv"  # the lab's patient split, for the plan's 313 / 251


def riga_features(photo):
    ref = riga.reference(photo)[0]
    return riga.descriptor(ref), riga.keypoints(ref)


def orb_block(job):
    """Inliers of every pair (i < j) between two blocks of photographs (features sent once per block pair)."""
    (ia, fa), (ib, fb) = job
    feats = dict(zip(ia, fa)) | dict(zip(ib, fb))
    pairs = [(i, j) for i in ia for j in ib if i < j]
    return splits.keypoint_duplicates(feats, threshold=0, pairs=pairs)


def riga_item(stats, photo, extracted, split, group, duplicates) -> Item:
    """The store item of one extracted photograph, counted into the report from the item itself (review #4: the
    count used kind == "prime", but riga.item keeps the F7 reconstruction's image, so test said 94 for 95)."""
    ref, kind, masks = extracted
    stats["riga_stages"].update(m[2] for m in masks.values())
    stats["riga_reference"][kind] += 1
    it = riga.item(photo, ref, kind, masks, split, group, duplicates)
    stats["riga_with_image"][split] += it.image is not None
    return it


def crowded(s):
    """B3-2: the item's model frame shows another outlined nodule of the slice."""
    it = Item("", "lidc", "", "train", None, {"nodule": s.full()}, [str(r) for r in s.readers], {})
    return lidc.crowded(s, frame(it).box)


def slice_flags(s):
    """Report counts of one slice: readers' per-instance topologies differ under the 5-px filter / the plan's size
    rule (with and without merged-reader clusters), copied layers, and (B2-4) a rater whose size-filtered signature
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

    photos = riga.discover(raw / "riga")
    ids = [p.image_id for p in photos]
    scans = lidc.scans()
    files = lidc.series_files(raw / "lidc_xml")
    assert len(photos) == 750 and set(files) == set(scans) and len(scans) == 1018
    with Pool(a.workers) as pool:
        # Pass 1: what the splits need. RIGA near-duplicate groups; LIDC slices and their disagreement.
        feats = pool.map(riga_features, photos, chunksize=4)
        ncc = splits.near_duplicates(np.stack([d for d, _ in feats]), threshold=-np.inf)
        blocks = [list(range(b, min(b + 94, len(ids)))) for b in range(0, len(ids), 94)]
        jobs = [((a, [feats[i][1] for i in a]), (b, [feats[i][1] for i in b]))
                for x, a in enumerate(blocks) for b in blocks[x:]]
        orb = sorted(p for res in pool.map(orb_block, jobs) for p in res)
        del feats, jobs
        dups = [(ids[i], ids[j], c) for i, j, c in ncc if c >= splits.NCC_DUP]
        kdups = [(ids[i], ids[j], n) for i, j, n in orb if n >= splits.ORB_DUP]
        ncc_pairs = {(x, y) for x, y, _ in dups}
        groups = splits.group_ids(ids, [(x, y) for x, y, _ in dups + kdups])
        sub = {p.image_id: p.subset for p in photos}
        riga_split = splits.riga_split(sub, groups, a.seed)
        lidc_res = pool.map(lidc_series, [(s, files[s], scans[s].zvals) for s in sorted(files)], chunksize=8)
        per_patient = {}
        for series, _, sl, fl, _ in lidc_res:
            n, d = per_patient.get(scans[series].patient, (0, 0))
            per_patient[scans[series].patient] = (n + len(sl), d + sum(f["disagree_5px"] for f in fl))
        patients = {s.patient: splits.Patient(*per_patient.get(s.patient, (0, 0))) for s in scans.values()}
        lidc_split = splits.lidc_split(patients, a.seed)
        rep["timing_s"]["pass1"] = round(time.time() - t0, 1)

        # Pass 2: extract RIGA masks and stream every item into the store.
        stats = {"riga_stages": Counter(), "riga_reference": Counter(), "image_sha": {}, "split_of": {},
                 "lidc": Counter(), "riga_with_image": Counter()}
        members = {}
        for i in ids:
            members.setdefault(groups[i], []).append(i)

        def items():
            step = 2 * a.workers  # bounded batches: a RIGA result is ~56 MB, an unbounded imap backlog is not
            for c in range(0, len(photos), step):
                for photo, ex in zip(photos[c:c + step], pool.map(riga.extract, photos[c:c + step])):
                    g = groups[photo.image_id]
                    yield riga_item(stats, photo, ex, riga_split[photo.image_id], g,
                                    [i for i in members[g] if i != photo.image_id])
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
    riga_counts = Counter(riga_split.values())
    test_groups = {groups[i] for i in ids if riga_split[i] == "test"}
    cross = splits.crossing([(ids[i], ids[j], c) for i, j, c in ncc], riga_split)
    kcross = splits.crossing([(ids[i], ids[j], n) for i, j, n in orb], riga_split)
    val_by_subset = Counter(sub[i] for i in ids if riga_split[i] == "val")
    groups_by_subset = Counter(sub[g] for g in {groups[i] for i in ids if riga_split[i] != "test"})
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
        "riga": {"photos": len(photos), "masks": sum(stats["riga_stages"].values()),
            "splits": dict(sorted(riga_counts.items())),
                 "test_unique_photos": len(test_groups), "stages": dict(stats["riga_stages"].most_common()),
                 "reference": dict(stats["riga_reference"]),
                 "with_image": dict(sorted(stats["riga_with_image"].items())),
                 "near_duplicate_pairs": len({(x, y) for x, y, _ in dups + kdups}),
                 "near_duplicate_pairs_ncc": len(dups), "near_duplicate_pairs_keypoints": len(kdups),
                 "keypoint_only_pairs": sorted([x, y, n] for x, y, n in kdups if (x, y) not in ncc_pairs),
                 "ncc_only_pairs": sorted(ncc_pairs - {(x, y) for x, y, _ in kdups}),
                 "min_keypoint_inliers_of_ncc_pairs": min(n for i, j, n in orb if (ids[i], ids[j]) in ncc_pairs),
                 "same_photo_pairs": sum(c >= 0.95 for *_, c in dups),
                 "duplicate_groups": sum(len(v) > 1 for v in members.values()),
                 "val_by_subset": dict(sorted(val_by_subset.items())),
                 "groups_by_subset": dict(sorted(groups_by_subset.items())),
                 "split_note": "I10: the validation carve-out deals 10 % of the duplicate *groups* of each subset "
                               "(largest remainder over the seeded hash order of group ids; val_by_subset vs "
                               "groups_by_subset), and a group joins val with all its photographs, so train / val photo "
                               "counts follow the grouping: 591 / 64 with the NCC groups only, 592 / 63 with the "
                               "keypoint groups (B3-8). The plan's 590 / 65 is the lab's draw; every draw keeps each "
                               "subset's val share within one group of 10 %."},
        "lidc": {"series": len(lidc_res), "patients": len(patients), "items": sum(len(sl) for _, _, sl, _,
            _ in lidc_res),
                 "dropped_crowded_frame": dict(sorted(crowded_per_split.items())),
                 "patients_per_split": dict(sorted(Counter(lidc_split.values()).items())),
                 **{k: stats["lidc"][k] for k in sorted(stats["lidc"])},
                 **{f"lab_split_train_{k}": sum(f[k] for f in lab_train)
                    for k in ("disagree_5px", "disagree_rule", "disagree_rule_no_dup_reader")}},
        "k11": {"riga_group_leaks": splits.leaks(riga_split, groups),
                "riga_cross_split_near_duplicates": [c for c in cross if c[2] >= splits.NCC_DUP],
                "riga_max_cross_split_ncc": round(max(c for *_, c in cross), 4),
                "riga_max_ncc_below_threshold": round(max(c for *_, c in ncc if c < splits.NCC_DUP), 4),
                "riga_cross_split_keypoint_duplicates": [c for c in kcross if c[2] >= splits.ORB_DUP],
                "riga_max_cross_split_keypoint_inliers": max(n for *_, n in kcross),
                "riga_max_keypoint_inliers_below_threshold": max(n for *_, n in orb if n < splits.ORB_DUP),
                "lidc_patient_leaks": splits.leaks({s: lidc_split[p] for s, p in patient_series.items()},
                    patient_series),
                "lidc_identical_slices_across_splits": splits.leaks(stats["split_of"], stats["image_sha"])},
    })
    (out / "build_report.json").write_text(json.dumps(rep, indent=1, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in rep.items() if k != "manifest"}, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
