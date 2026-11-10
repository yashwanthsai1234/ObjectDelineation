"""K11 split tests (ported from archive/phase1_lab/t_data/splits/test_splits.py): group integrity,
no cross-split near-duplicates, determinism, stratification, and
injected-leak negative controls. Synthetic populations; the real store is checked by the data-backed tests at
the end."""
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy import ndimage

from objdel.data import riga, splits, store

ROOT = Path(__file__).resolve().parents[2]
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


def riga_population():
    subsets = {}
    for sub, n in [("BinRushed/BinRushed1", 50), ("BinRushed/BinRushed2", 47), ("BinRushed/BinRushed3", 47),
                   ("BinRushed/BinRushed4", 51), ("MESSIDOR", 460), ("Magrabia/MagrabiFemale", 47),
                   ("Magrabia/MagrabiaMale", 48)]:
        for i in range(1, n + 1):
            subsets[f"{sub}/image{i:03d}"] = sub
    pairs = [("MESSIDOR/image001", "MESSIDOR/image002"), ("BinRushed/BinRushed2/image005", "MESSIDOR/image010"),
             ("Magrabia/MagrabiFemale/image003", "Magrabia/MagrabiaMale/image004"),
             ("MESSIDOR/image020", "Magrabia/MagrabiaMale/image007")]  # a train photo duplicated in test
    return subsets, splits.group_ids(sorted(subsets), pairs), pairs


@pytest.fixture(scope="module")
def lidc():
    pats = lidc_population()
    return pats, splits.lidc_split(pats, seed=0)


@pytest.fixture(scope="module")
def rig():
    subsets, groups, pairs = riga_population()
    return subsets, groups, pairs, splits.riga_split(subsets, groups, seed=0)


# ---- leakage -------------------------------------------------------------------------------------------------

def test_lidc_every_patient_assigned_once(lidc):
    pats, assign = lidc
    assert set(assign) == set(pats) and set(assign.values()) == {"train", "val", "test"}


def test_lidc_series_follow_their_patient(lidc):
    _, assign = lidc
    series = {f"1.2.{i}": p for i, p in enumerate(sorted(assign))} | {"9.9.1": "LIDC-IDRI-0007"}  # a two-scan patient
    split_of = {s: assign[p] for s, p in series.items()}
    assert splits.leaks(split_of, series) == []


def test_riga_literature_split_and_group_leakage(rig):
    subsets, groups, _, assign = rig
    assert all(assign[i] == "test" for i, s in subsets.items() if s.startswith("Magrabia"))
    assert assign["MESSIDOR/image020"] == "test"  # duplicated in Magrabia: the whole group goes to test
    assert sum(v == "test" for v in assign.values()) == 96
    assert splits.leaks(assign, groups) == []


def test_no_near_duplicates_across_splits(rig):
    _, _, pairs, assign = rig
    assert splits.crossing(pairs, assign) == []


# ---- determinism ---------------------------------------------------------------------------------------------

def test_same_seed_identical(lidc, rig):
    subsets, groups, _, assign = rig
    assert splits.lidc_split(lidc_population(), seed=0) == lidc[1]
    assert splits.riga_split(subsets, groups, seed=0) == assign


def test_different_seed_differs(lidc, rig):
    subsets, groups, _, assign = rig
    assert splits.lidc_split(lidc_population(), seed=1) != lidc[1]
    assert splits.riga_split(subsets, groups, seed=1) != assign


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


def test_riga_val_carveout_stratified(rig):
    subsets, _, _, assign = rig
    by = {}
    for i, s in subsets.items():
        if assign[i] != "test":
            by.setdefault(s, Counter())[assign[i]] += 1
    for sub, c in by.items():
        n = sum(c.values())
        assert c["val"] > 0 and abs(c["val"] - 0.1 * n) <= 1.0 + 1e-9, (sub, dict(c))


# ---- near-duplicate detection and negative controls --------------------------------------------------------------

def texture(seed, shift=(0, 0)):
    rng = np.random.default_rng(seed)
    img = ndimage.gaussian_filter(rng.random((420, 620)), 3)
    img = np.roll(img, shift, axis=(0, 1))[10:394, 20:596]
    g = (255 * (img - img.min()) / np.ptp(img)).astype(np.uint8)
    return np.stack([g, g // 2, g // 4], -1)


def test_near_duplicates_found_by_vessel_texture():
    imgs = [texture(0), texture(1), texture(0, shift=(6, -10)), texture(2)]  # 2 is image 0 moved by (6, -10) px
    pairs = splits.near_duplicates(np.stack([riga.descriptor(i) for i in imgs]))
    assert [(a, b) for a, b, _ in pairs] == [(0, 2)] and pairs[0][2] >= splits.NCC_DUP


def reframed(img, angle=2.5, scale=0.96, shift=(30, -5)):
    import cv2
    h, w = img.shape[:2]
    A = cv2.getRotationMatrix2D((w / 2, h / 2), angle, scale)
    A[:, 2] += shift
    return cv2.warpAffine(img, A, (w, h), borderMode=cv2.BORDER_REFLECT)


def test_reframed_same_eye_found_by_keypoints():  # B3-8: rotation / scale / shift defeat the NCC, not keypoints
    imgs = [texture(0), texture(1), reframed(texture(0)), texture(2)]
    ncc = splits.near_duplicates(np.stack([riga.descriptor(i) for i in imgs]))
    assert (0, 2) not in [(a, b) for a, b, _ in ncc]  # the shift-only NCC misses the reframed photograph
    feats = [riga.keypoints(i) for i in imgs]
    hits = splits.keypoint_duplicates(feats)
    assert [(a, b) for a, b, _ in hits] == [(0, 2)] and hits[0][2] >= splits.ORB_DUP
    every = splits.keypoint_duplicates(feats, threshold=0)
    assert [(a, b) for a, b, _ in every] == [(a, b) for a in range(4) for b in range(a + 1, 4)]
    assert max(n for a, b, n in every if (a, b) != (0, 2)) < splits.ORB_DUP // 2
    got = splits.keypoint_duplicates(feats, threshold=0, pairs=[(2, 0), (1, 3)])  # given pairs, in order
    assert [(a, b) for a, b, _ in got] == [(2, 0), (1, 3)] and got[0][2] >= splits.ORB_DUP > got[1][2]
    assert splits.keypoint_duplicates([riga.keypoints(i) for i in imgs]) == hits  # deterministic


def test_keypoints_of_a_blank_image_match_nothing():
    blank = np.zeros((384, 576, 3), np.uint8)
    feats = [riga.keypoints(blank), riga.keypoints(texture(0))]
    assert splits.keypoint_duplicates(feats, threshold=0) == [(0, 1, 0)]


def test_negative_control_near_duplicate_detected(rig):
    _, _, pairs, assign = rig
    bad = dict(assign)
    a = pairs[0][0]
    bad[a] = "val" if bad[a] != "val" else "train"
    assert splits.crossing(pairs, assign) == [] and splits.crossing(pairs, bad) == [pairs[0]]


def test_negative_control_patient_leak_detected(lidc):
    _, assign = lidc
    series = {"1.1": "LIDC-IDRI-0007", "1.2": "LIDC-IDRI-0007", "2.1": "LIDC-IDRI-0008"}
    split_of = {s: assign[p] for s, p in series.items()}
    bad = dict(split_of, **{"1.2": "test" if split_of["1.1"] != "test" else "train"})  # a per-series splitter
    assert splits.leaks(split_of, series) == [] and splits.leaks(bad, series) == ["LIDC-IDRI-0007"]


# ---- data-backed (K11 on the built store) -------------------------------------------------------------------------

@pytest.mark.data
def test_store_k11(store_dir):
    """0 groups across splits, 0 two-scan patients split, 0 cross-split near-duplicates (RIGA vessel texture, LIDC
    identical CT slices), splits at 70/10/20 patients."""
    import hashlib

    from objdel.data import lidc
    split_of, group_of, desc, img_hash = {}, {}, {}, {}
    for ds in ("riga", "lidc"):
        for split in ("train", "val", "test"):
            for it in store.load(ds, split, store_dir):
                split_of[it.item_id], group_of[it.item_id] = it.split, it.group_id
                if ds == "riga" and it.image is not None:  # B3-1: reconstructed references have no image
                    desc[it.item_id] = riga.descriptor(it.image)
                elif it.image is not None:
                    img_hash[it.item_id] = hashlib.sha256(it.image.tobytes()).hexdigest()
    assert splits.leaks(split_of, group_of) == []
    ids = sorted(desc)
    pairs = [(ids[a], ids[b]) for a, b, _ in splits.near_duplicates(np.stack([desc[i] for i in ids]))]
    # 37 NCC pairs in the build (from the raw references); 9 involve a BinRushed1 photograph, which has no image
    assert len(pairs) == 28 and splits.crossing(pairs, split_of) == []
    assert all(group_of[a] == group_of[b] for a, b in pairs)
    k11 = json.loads((store_dir / "build_report.json").read_text())["k11"]  # built from raw, incl. BinRushed1
    assert all(k11[k] == [] for k in ("riga_group_leaks", "riga_cross_split_near_duplicates",
                                      "riga_cross_split_keypoint_duplicates", "lidc_patient_leaks",
                                      "lidc_identical_slices_across_splits"))
    assert splits.leaks(split_of, img_hash) == []
    scans = lidc.scans()
    series_split = {i.split("/")[1]: s for i, s in split_of.items() if i.startswith("lidc/")}
    assert splits.leaks(series_split, {s: scans[s].patient for s in series_split}) == []  # incl. two-scan patients
    shares = Counter(s for i, s in split_of.items() if i.startswith("lidc/"))
    assert shares["train"] > shares["test"] > shares["val"]


@pytest.mark.data
@pytest.mark.slow
def test_store_k11_keypoints(store_dir):
    """B3-8: no cross-split pair of RIGA photographs (with an image) reaches the keypoint duplicate threshold, and the
    break-it pair BinRushed3/image004 (val group) + BinRushed4/image002 shares one group and split."""
    feats, split_of, group_of = [], [], []
    for split in ("train", "val", "test"):
        for it in store.load("riga", split, store_dir):
            if it.image is not None:
                feats.append(riga.keypoints(it.image))
                split_of.append(it.split)
                group_of.append((it.item_id, it.group_id))
    pairs = [(i, j) for i in range(len(feats)) for j in range(i + 1, len(feats)) if split_of[i] != split_of[j]]
    assert splits.keypoint_duplicates(feats, pairs=pairs) == []
    g = dict(group_of)
    assert g["riga/BinRushed/BinRushed3/image004"] == g["riga/BinRushed/BinRushed4/image002"]


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
