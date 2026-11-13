"""K11 split tests: group integrity, determinism, stratification and an injected-leak negative control. Synthetic
populations; the real store is checked by the data-backed tests at the end."""
import csv
import json
from collections import Counter

import numpy as np
import pytest

from objdel.data import splits, store

TARGET = dict(splits.LIDC_FRACTIONS)


def lidc_population(seed=0, n=1010, n_dis=120):
    """Patients like LIDC: many with no multi-reader slices, a heavy-tailed disagreement count on a minority."""
    rng = np.random.default_rng(seed)
    pats = {}
    for i in range(n):
        slices = int(rng.choice([0, rng.integers(1, 11), rng.integers(11, 41), rng.integers(41, 120)]))
        dis = int(rng.pareto(1.5) + 1) if i < n_dis and slices else 0
        pats[f"LIDC-IDRI-{i:04d}"] = splits.Patient(slices=slices, disagree=min(dis, slices))
    return pats


@pytest.fixture(scope="module")
def lidc():
    pats = lidc_population()
    return pats, splits.lidc_split(pats, seed=0)


# ---- leakage -------------------------------------------------------------------------------------------------

def test_lidc_every_patient_assigned_once(lidc):
    pats, assign = lidc
    assert set(assign) == set(pats) and set(assign.values()) == {"train", "val", "test"}


def test_lidc_series_follow_their_patient(lidc):
    _, assign = lidc
    series = {f"1.2.{i}": p for i, p in enumerate(sorted(assign))} | {"9.9.1": "LIDC-IDRI-0007"}  # a two-scan patient
    split_of = {s: assign[p] for s, p in series.items()}
    assert splits.leaks(split_of, series) == []


# ---- determinism ---------------------------------------------------------------------------------------------

def test_same_seed_identical(lidc):
    assert splits.lidc_split(lidc_population(), seed=0) == lidc[1]


def test_different_seed_differs(lidc):
    assert splits.lidc_split(lidc_population(), seed=1) != lidc[1]


# ---- stratification ------------------------------------------------------------------------------------------

def test_lidc_strata_within_one_patient_of_target(lidc):
    pats, assign = lidc
    by = {}
    for p, u in pats.items():
        by.setdefault(splits.lidc_stratum(u), Counter())[assign[p]] += 1
    for st, c in by.items():
        if st == "dis=1":
            continue  # weighted by disagreement slices, tested below
        n = sum(c.values())
        for name, f in TARGET.items():
            assert abs(c[name] - f * n) <= 1.0 + 1e-9, (st, dict(c))


def test_lidc_disagreement_in_every_split_within_tolerance(lidc):
    pats, assign = lidc
    tot = sum(u.disagree for u in pats.values())
    totp = sum(u.disagree > 0 for u in pats.values())
    for name, f in TARGET.items():
        us = [u for p, u in pats.items() if assign[p] == name]
        assert abs(sum(u.disagree for u in us) / tot - f) <= 0.02
        assert abs(sum(u.disagree > 0 for u in us) / totp - f) <= 0.03
        assert abs(len(us) / len(pats) - f) <= 0.01


# ---- negative control ------------------------------------------------------------------------------------------

def test_negative_control_patient_leak_detected(lidc):
    _, assign = lidc
    series = {"1.1": "LIDC-IDRI-0007", "1.2": "LIDC-IDRI-0007", "2.1": "LIDC-IDRI-0008"}
    split_of = {s: assign[p] for s, p in series.items()}
    bad = dict(split_of, **{"1.2": "test" if split_of["1.1"] != "test" else "train"})  # a per-series splitter
    assert splits.leaks(split_of, series) == [] and splits.leaks(bad, series) == ["LIDC-IDRI-0007"]


# ---- data-backed (K11 on the built store) -------------------------------------------------------------------------

@pytest.mark.data
def test_store_k11(store_dir):
    """0 two-scan patients split, 0 identical CT slices across splits, splits at 70/10/20 patients."""
    import hashlib

    from objdel.data import lidc
    split_of, group_of, img_hash = {}, {}, {}
    for split in ("train", "val", "test"):
        for it in store.load("lidc", split, store_dir):
            split_of[it.item_id], group_of[it.item_id] = it.split, it.group_id
            if it.image is not None:
                img_hash[it.item_id] = hashlib.sha256(it.image.tobytes()).hexdigest()
    assert splits.leaks(split_of, group_of) == []
    k11 = json.loads((store_dir / "build_report.json").read_text())["k11"]
    assert all(k11[k] == [] for k in ("lidc_patient_leaks", "lidc_identical_slices_across_splits"))
    assert splits.leaks(split_of, img_hash) == []
    scans = lidc.scans()
    series_split = {i.split("/")[1]: s for i, s in split_of.items() if i.startswith("lidc/")}
    assert splits.leaks(series_split, {s: scans[s].patient for s in series_split}) == []  # incl. two-scan patients
    shares = Counter(s for i, s in split_of.items() if i.startswith("lidc/"))
    assert shares["train"] > shares["test"] > shares["val"]


@pytest.mark.data
def test_lidc_train_disagreement_counts_on_the_lab_split(store_dir, fixtures_dir):
    """The plan's numbers on the lab's patient split were 313 train slices disagreeing under the 5-px filter and 251
    under the size rule max(5, ceil(0.5 % of the rater union)) once clusters with a merged reader are left out (as
    the lab's reward fixture did; 308 with them). After the wave-6 rules (copied layers count once, crowded frames
    left out; fixes/F4.md) the same counts are 296 and 242 (291 with merged-reader clusters)."""
    from objdel.geom.topology import min_area, signature
    lab_split = {r["patient"]: r["split"] for r in csv.DictReader(open(fixtures_dir / "lidc_full_meta.csv"))}
    five = rule = 0
    for split in ("train", "val", "test"):
        for it in store.load("lidc", split, store_dir):
            if lab_split[it.group_id] != "train":
                continue
            m = it.masks["nodule"]
            five += len({signature(x, 5) for x in m}) > 1
            rule += not it.meta["dup_reader_in_cluster"] and len(
                {signature(x, min_area(int(m.any(0).sum()))) for x in m}) > 1
    assert (five, rule) == (296, 242)
