"""K11 split tests: zero group leakage, zero cross-split near-duplicates, determinism, stratification.

    ../../.venv/bin/python -m pytest -q test_splits.py
"""
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pytest

HERE = Path(__file__).parent
TD = HERE.parent
REPO = Path("/Users/sreeyashwanthsaivenkatesh/Workspace/projects/CUDA/ObjectDelineation")
sys.path.insert(0, str(HERE))
import build_splits  # noqa: E402

INPUTS = ["--lidc-clusters", str(REPO / "experiments/phase0/results/lidc/clusters.csv"),
          "--series-patient", str(TD / "series_patient.csv"),
          "--riga-images", str(TD / "riga/out/images.csv"),  # fixed extractor output
          "--riga-dups", str(TD / "riga/dedup/dup_groups.json")]
PAIRS = TD / "riga/dedup/pairs.json"


def build(out, seed=0):
    build_splits.main(["--seed", str(seed), "--out", str(out)] + INPUTS)
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(out).glob("*.csv"))}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("s0")
    build(out)
    lidc = list(csv.DictReader(open(out / "lidc_patients.csv")))
    riga = list(csv.DictReader(open(out / "riga_images.csv")))
    return out, lidc, riga


# ---- leakage ---------------------------------------------------------------

def test_lidc_patient_unique_and_complete(built):
    _, lidc, _ = built
    pids = [r["patient_id"] for r in lidc]
    assert len(pids) == len(set(pids))
    sp = list(csv.DictReader(open(TD / "series_patient.csv")))
    assert set(pids) == {r["patient_id"] for r in sp}


def test_lidc_series_in_exactly_one_split_and_two_scan_patients_together(built):
    _, lidc, _ = built
    split_of_series = defaultdict(set)
    for r in lidc:
        for s in r["series_uids"].split(";"):
            split_of_series[s].add(r["split"])
    assert all(len(v) == 1 for v in split_of_series.values())
    sp = list(csv.DictReader(open(TD / "series_patient.csv")))
    patient_splits = defaultdict(set)
    for r in sp:
        patient_splits[r["patient_id"]] |= split_of_series[r["series_uid"]]
    assert all(len(v) == 1 for v in patient_splits.values())
    two = [p for p, n in Counter(r["patient_id"] for r in sp).items() if n > 1]
    assert len(two) == 8


def test_riga_literature_split_and_group_leakage(built):
    _, _, riga = built
    assert len(riga) == 750
    assert all(r["split"] == "test" for r in riga if r["subset"].startswith("Magrabia"))
    assert sum(r["split"] == "test" for r in riga) == 95
    assert sum(r["split"] in ("train", "val") for r in riga) == 655
    groups = defaultdict(set)
    for r in riga:
        groups[r["group_id"]].add(r["split"])
    assert all(len(v) == 1 for v in groups.values())


def test_no_near_duplicates_across_splits(built):
    _, _, riga = built
    split = {r["image_id"].replace("/", "_", r["image_id"].count("/") - 1): r["split"] for r in riga}
    pairs = json.load(open(PAIRS))
    crossing = [p for p in pairs if split[p["a"]] != split[p["b"]]]
    assert crossing == []


# ---- determinism -----------------------------------------------------------

def test_same_seed_identical_hashes(tmp_path):
    assert build(tmp_path / "a", 0) == build(tmp_path / "b", 0)


def test_different_seed_differs(tmp_path):
    assert build(tmp_path / "a", 0) != build(tmp_path / "c", 1)


# ---- stratification --------------------------------------------------------

TARGET = dict(build_splits.LIDC_FRACTIONS)


def test_lidc_strata_within_one_patient_of_target(built):
    _, lidc, _ = built
    by = defaultdict(Counter)
    for r in lidc:
        by[r["stratum"]][r["split"]] += 1
    for st, c in by.items():
        if st == "dis=1":  # weighted by disagreement slices, tested below
            continue
        n = sum(c.values())
        for name, f in TARGET.items():
            assert abs(c[name] - f * n) <= 1.0 + 1e-9, (st, dict(c))


def test_lidc_disagreement_in_every_split_within_tolerance(built):
    _, lidc, _ = built
    tot = sum(int(r["disagree_slices"]) for r in lidc)
    totp = sum(int(r["disagree_slices"]) > 0 for r in lidc)
    for name, f in TARGET.items():
        rs = [r for r in lidc if r["split"] == name]
        d = sum(int(r["disagree_slices"]) for r in rs)
        assert d > 0
        assert abs(d / tot - f) <= 0.02, (name, d / tot)
        assert abs(sum(int(r["disagree_slices"]) > 0 for r in rs) / totp - f) <= 0.03
        assert abs(len(rs) / len(lidc) - f) <= 0.01


def test_riga_val_carveout_stratified(built):
    _, _, riga = built
    by = defaultdict(Counter)
    for r in riga:
        if r["split"] in ("train", "val"):
            by[r["subset"]][r["split"]] += 1
    for sub, c in by.items():
        n = sum(c.values())
        assert c["val"] > 0, sub
        assert abs(c["val"] - 0.1 * n) <= 1.0 + 1e-9, (sub, dict(c))


# ---- negative controls: the checks must catch injected leakage ----------------

def crossing_pairs(riga):
    split = {r["image_id"].replace("/", "_", r["image_id"].count("/") - 1): r["split"] for r in riga}
    return [p for p in json.load(open(PAIRS)) if split[p["a"]] != split[p["b"]]]


def test_negative_control_near_duplicate_detected(built):
    _, _, riga = built
    pairs = json.load(open(PAIRS))
    a = pairs[0]["a"]
    bad = [dict(r) for r in riga]
    for r in bad:
        if r["image_id"].replace("/", "_", r["image_id"].count("/") - 1) == a:
            r["split"] = "val" if r["split"] != "val" else "train"
    assert crossing_pairs(riga) == [] and len(crossing_pairs(bad)) >= 1


def patient_leaks(series_split):
    sp = list(csv.DictReader(open(TD / "series_patient.csv")))
    by = defaultdict(set)
    for r in sp:
        by[r["patient_id"]].add(series_split[r["series_uid"]])
    return [p for p, v in by.items() if len(v) > 1]


def test_patient_leak_check_and_negative_control(built):
    _, lidc, _ = built
    series_split = {s: r["split"] for r in lidc for s in r["series_uids"].split(";")}
    assert patient_leaks(series_split) == []
    two = next(r for r in lidc if int(r["n_series"]) > 1)
    s1, s2 = two["series_uids"].split(";")
    bad = dict(series_split)
    bad[s2] = "test" if bad[s1] != "test" else "train"  # what a per-series splitter could do
    assert patient_leaks(bad) == [two["patient_id"]]
