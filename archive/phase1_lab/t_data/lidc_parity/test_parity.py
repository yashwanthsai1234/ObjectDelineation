"""Parity gates: our LIDC XML parser/clustering vs pylidc 0.2.3 (run parity.py first).

    ../../.venv/bin/python -I -m pytest -q test_parity.py
"""
import csv
import json
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
S = json.loads((HERE / "SUMMARY.json").read_text())
ROWS = list(csv.DictReader(open(HERE / "slices.csv")))


def test_sample_size():
    assert S["scans"] >= 50


def test_annotation_sets_match():
    a = S["annotations"]
    assert a["ours"] == a["pylidc"]
    assert a["pylidc_unmatched_contours"] == 0 and a["our_only"] == 0 and a["pylidc_only"] == 0
    assert a["nodule_id_mismatch"] == 0


def test_clustering_parity():
    assert S["nodule_count_equal_scans"] == S["scans"]
    assert S["partition_equal_scans"] == S["scans"]
    assert S["reader_sets_equal_scans"] == S["scans"]


def test_mask_convention_explained():
    # ours == pylidc boolean_mask(include_contour_points=True) OR outline pixels, on >= 99.5 % of slices
    m = S["slice_masks"]["iou_ours_vs_pylidc_True_plus_outline"]
    assert m["pct_exact_1"] >= 99.5 and m["p5"] == 1.0
    # and the raw conventions really differ (outline-inclusion matters ~5 % of area at median)
    assert S["slice_masks"]["iou_vs_pylidc_include_contour_points_True"]["median"] < 0.97


def test_only_known_zero_iou_slices():
    # ours non-empty / pylidc empty only for 1-pixel single-point ROIs or all-outline tiny contours
    zero = [r for r in ROWS if float(r["iou_incl"]) == 0]
    assert all(int(r["area_ours"]) <= 4 and int(r["area_pyl_incl"]) == 0 for r in zero)


def test_z_matching():
    z = S["z"]
    assert z["sop_multi_z"] == 0 and z["z_multi_sop"] == 0 and z["scans_dup_slice_z"] == 0
    assert z["z_not_on_slice"] <= 2  # only 188/011.xml (sessions 5-8 belong to another series)


def test_series_patient_map():
    rows = list(csv.DictReader(open(HERE.parent / "series_patient.csv")))
    assert len(rows) == 1018 and all(r["patient_id"].startswith("LIDC-IDRI-") for r in rows)
    assert len({r["series_uid"] for r in rows}) == 1018
    per = {}
    for r in rows:
        per.setdefault(r["patient_id"], []).append(r["series_uid"])
    assert len(per) == 1010
    assert sorted(p for p, v in per.items() if len(v) > 1) == sorted(S["two_scan_patients"])
    assert len(S["two_scan_patients"]) == 8


@pytest.mark.slow
def test_live_rerun_one_scan():
    import parity
    by, ch = parity.series_files()
    s = "1.3.6.1.4.1.14519.5.2.1.6279.6001.100398138793540579077826395208"
    o, rows = parity.analyse((s, str(ch[s])))
    assert o["partition_equal"] and o["our_nodules"] == o["pylidc_nodules"] == 3
    assert all(r["iou_ours_vs_incl_plus_outline"] == 1.0 for r in rows)
