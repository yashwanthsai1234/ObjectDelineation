"""Deterministic, seeded, group-aware split builder (prototype).

LIDC: unit = patient (PatientID from the pylidc DB; the XMLs carry none), so the
8 two-scan patients stay together by construction. 70/10/20 by patient,
stratified: patients without disagreement by multi-reader slice-count bin
(count dealing); patients with reader-topology-disagreement slices by weighted
greedy on their disagreement-slice count, so those slices split ~70/10/20.

RIGA: literature split. test = Magrabia (95); train = BinRushed1-4 + MESSIDOR
(655) with a 10 % validation carve-out grouped by image (near-duplicate groups
from the dedup audit merged into one group) and stratified by subset.

Assignment: within each stratum, groups are ordered by SHA-256(seed, group id)
and dealt into splits by a largest-remainder quota schedule, so the result
depends only on the seed and the population (no RNG state, no library version).

    python -I build_splits.py --seed 0 --out <dir> --lidc-clusters <clusters.csv> --series-patient <csv>
        --riga-images <riga images.csv> [--riga-dups <dup_groups.json>]
"""
import argparse
import csv
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path

LIDC_FRACTIONS = (("train", 0.7), ("val", 0.1), ("test", 0.2))
RIGA_VAL_FRACTION = 0.1


def key(seed, gid):
    return hashlib.sha256(f"{seed}:{gid}".encode()).hexdigest()


def deal(groups, fractions, seed):
    """groups: list of ids (one stratum). Returns {id: split} hitting fractions per stratum.

    Largest-remainder schedule: the k-th group goes to the split whose running
    count is furthest below its target k*f, ties broken by the fraction order.
    """
    order = sorted(groups, key=lambda g: key(seed, g))
    counts = Counter()
    out = {}
    for k, g in enumerate(order, 1):
        name = max(fractions, key=lambda nf: (k * nf[1] - counts[nf[0]], nf[1]))[0]
        counts[name] += 1
        out[g] = name
    return out


def deal_weighted(groups, weight, fractions, seed):
    """Weighted greedy: heaviest group first (ties by seeded hash) to the split
    furthest below its target share of the total weight. Balances disagreement
    slices, which are concentrated in few patients."""
    order = sorted(groups, key=lambda g: (-weight[g], key(seed, g)))
    total = sum(weight[g] for g in groups) or 1
    have = Counter()
    out = {}
    for g in order:
        name = max(fractions, key=lambda nf: ((nf[1] * total - have[nf[0]]) / nf[1], nf[1]))[0]
        have[name] += weight[g]
        out[g] = name
    return out


def lidc_units(clusters_csv, series_patient_csv):
    patient = {r["series_uid"]: r["patient_id"] for r in csv.DictReader(open(series_patient_csv))}
    units = defaultdict(lambda: {"series": set(), "clusters": 0, "slices": 0, "disagree": 0, "multi_clusters": 0})
    for s, p in patient.items():
        units[p]["series"].add(s)
    for r in csv.DictReader(open(clusters_csv)):
        u = units[patient[r["series"]]]
        u["clusters"] += 1
        if int(r["readers_outlined"]) >= 2:
            u["multi_clusters"] += 1
            u["slices"] += int(r["slices_compared_5"])
            u["disagree"] += int(r["slices_disagree_5"])
    return units


def lidc_stratum(u):
    if u["disagree"] > 0:
        return "dis=1"  # weighted by disagreement slices, see deal_weighted
    s = u["slices"]
    size = "0" if s == 0 else "1-10" if s <= 10 else "11-40" if s <= 40 else "41+"
    return f"dis=0|slices={size}"


def build_lidc(units, seed):
    strata = defaultdict(list)
    for p, u in units.items():
        strata[lidc_stratum(u)].append(p)
    assign = {}
    for st in sorted(strata):
        if st == "dis=1":
            w = {p: units[p]["disagree"] for p in strata[st]}
            assign.update(deal_weighted(strata[st], w, LIDC_FRACTIONS, seed))
        else:
            assign.update(deal(strata[st], LIDC_FRACTIONS, seed))
    rows = []
    for p in sorted(units):
        u = units[p]
        rows.append({"patient_id": p, "split": assign[p], "stratum": lidc_stratum(u), "n_series": len(u["series"]),
                     "series_uids": ";".join(sorted(u["series"])), "clusters": u["clusters"],
                     "multi_reader_clusters": u["multi_clusters"], "multi_reader_slices": u["slices"],
                     "disagree_slices": u["disagree"]})
    return rows


def riga_groups(images_csv, dups_json=None):
    """Group id per image: the image itself, or its near-duplicate group (dedup ids use '_' for '/')."""
    imgs = list(csv.DictReader(open(images_csv)))
    gid = {f"{r['subset']}/image{r['image']}": f"riga:{r['subset']}/image{r['image']}" for r in imgs}
    back = {k.replace("/", "_", k.count("/") - 1): k for k in gid}
    if dups_json and Path(dups_json).exists():
        for rep, members in sorted(json.load(open(dups_json)).items()):
            for m in members:
                gid[back[m]] = "riga:dup:" + back[rep]
    return imgs, gid


def build_riga(images_csv, dups_json, seed):
    imgs, gid = riga_groups(images_csv, dups_json)
    rows = []
    for r in imgs:
        iid = f"{r['subset']}/image{r['image']}"
        test = r["subset"].startswith("Magrabia")
        rows.append({"image_id": iid, "group_id": gid[iid], "subset": r["subset"], "source": "Magrabia" if test else r["subset"].split("/")[0],
                     "split": "test" if test else None, "n_masks": int(r["masks"]), "ref": r["ref"]})
    # A duplicate group that touches the test set goes entirely to test (never seen in training).
    test_groups = {r["group_id"] for r in rows if r["split"] == "test"}
    for r in rows:
        if r["group_id"] in test_groups:
            r["split"] = "test"
    # Validation carve-out: groups stratified by the subset of their first member.
    by_group = defaultdict(list)
    for r in rows:
        if r["split"] is None:
            by_group[r["group_id"]].append(r)
    strata = defaultdict(list)
    for g, rs in by_group.items():
        strata[sorted(r["subset"] for r in rs)[0]].append(g)
    fr = (("train", 1 - RIGA_VAL_FRACTION), ("val", RIGA_VAL_FRACTION))
    assign = {}
    for st in sorted(strata):
        assign.update(deal(strata[st], fr, seed))
    for r in rows:
        if r["split"] is None:
            r["split"] = assign[r["group_id"]]
    rows.sort(key=lambda r: r["image_id"])
    return rows


def write_csv(path, rows):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    data = buf.getvalue().encode()
    Path(path).write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def stats_lidc(rows):
    s = {}
    for name, _ in LIDC_FRACTIONS:
        rs = [r for r in rows if r["split"] == name]
        s[name] = {"patients": len(rs), "series": sum(r["n_series"] for r in rs),
                   "two_scan_patients": sum(r["n_series"] > 1 for r in rs),
                   "patients_with_multi_reader_slices": sum(r["multi_reader_slices"] > 0 for r in rs),
                   "multi_reader_slices": sum(r["multi_reader_slices"] for r in rs),
                   "disagree_slices": sum(r["disagree_slices"] for r in rs),
                   "patients_with_disagree": sum(r["disagree_slices"] > 0 for r in rs)}
    tot = {k: sum(s[n][k] for n in s) for k in ("patients", "multi_reader_slices", "disagree_slices")}
    for n in s:
        for k in tot:
            s[n][f"{k}_share"] = round(s[n][k] / max(1, tot[k]), 4)
    return s


def stats_riga(rows):
    s = {}
    for name in ("train", "val", "test"):
        rs = [r for r in rows if r["split"] == name]
        s[name] = {"images": len(rs), "groups": len({r["group_id"] for r in rs}), "masks": sum(r["n_masks"] for r in rs),
                   "by_subset": dict(sorted(Counter(r["subset"] for r in rs).items()))}
    return s


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--lidc-clusters", required=True)
    ap.add_argument("--series-patient", required=True)
    ap.add_argument("--riga-images", required=True)
    ap.add_argument("--riga-dups", default=None)
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    lidc = build_lidc(lidc_units(a.lidc_clusters, a.series_patient), a.seed)
    riga = build_riga(a.riga_images, a.riga_dups, a.seed)
    hashes = {"lidc_patients.csv": write_csv(out / "lidc_patients.csv", lidc),
              "riga_images.csv": write_csv(out / "riga_images.csv", riga)}
    manifest = {"seed": a.seed, "files": hashes, "lidc": stats_lidc(lidc), "riga": stats_riga(riga),
                "lidc_strata": {st: dict(sorted(Counter(r["split"] for r in lidc if r["stratum"] == st).items()))
                                for st in sorted({r["stratum"] for r in lidc})}}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
