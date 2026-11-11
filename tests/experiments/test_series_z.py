"""B1.7 contour-to-slice check of every series and gate DG17 (Implementation.md B1.7; tests I17, I35, I37, I45) on a
five-series fixture: A (train) has two SOP UIDs at z 123.5, B (train) an ROI 2 mm from the nearest slice at 1.25 mm
spacing, C (train, LIDC-IDRI-0396) is clean, T (test) has a sidecar with fewer slices than cited SOPs, V (val) a
sidecar with slices matched by z."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from objdel.data import lidc, quality, store
from tools.data import check_series_z

A, B, C, T, V = (f"1.2.840.{i}00001" for i in range(1, 6))  # suffix (last 6 digits) "100001", ...
SPLIT = {A: "train", B: "train", C: "train", T: "test", V: "val"}
PATIENT = {A: "LIDC-IDRI-0013", B: "LIDC-IDRI-0566", C: "LIDC-IDRI-0396", T: "LIDC-IDRI-0659", V: "LIDC-IDRI-0085"}
NITEMS = {A: 3, B: 2, C: 1, T: 2, V: 2}
ZVALS = np.arange(120.0, 130.0, 1.25)  # 120.0, 121.25, ..., 128.75
ROIS = {A: [(123.5, "1.9.1"), (123.5, "1.9.2"), (124.75, "1.9.3")],
        B: [(121.25, "1.9.4"), (129.3, "1.9.10"), (129.45, "1.9.11"), (130.75, "1.9.5")],  # 0.55, 0.7, 2.0 mm off
        C: [(1420.0, "1.9.6")],
        T: [(122.5, "1.9.7")],
        V: [(126.0, "1.9.8")]}
SIDECAR = {A: (3, 3, 0), B: (4, 4, 0), C: (1, 1, 0), T: (1, 2, 0), V: (1, 1, 1)}  # slices, annotated_sops, by z


def message(series, rois) -> str:
    body = "".join(f"<roi><imageZposition>{z}</imageZposition><imageSOP_UID>{sop}</imageSOP_UID>"
                   "<inclusion>TRUE</inclusion><edgeMap><xCoord>1</xCoord><yCoord>1</yCoord></edgeMap></roi>"
                   for z, sop in rois)
    return ('<?xml version="1.0" encoding="UTF-8"?><LidcReadMessage xmlns="http://www.nih.gov"><ResponseHeader>'
            f"<SeriesInstanceUid>{series}</SeriesInstanceUid></ResponseHeader><readingSession>"
            f"<unblindedReadNodule>{body}</unblindedReadNodule></readingSession></LidcReadMessage>")


def ids(series) -> list[str]:
    return [f"lidc/{series}/n000/z{k:03d}" for k in range(NITEMS[series])]


def write_store(root, rois=ROIS):
    items = [store.Item(i, "lidc", PATIENT[s], SPLIT[s], None, {"nodule": np.ones((2, 4, 4), bool)}, ["1", "2"], {})
             for s in sorted(SPLIT) for i in ids(s)]
    store.write(sorted(items, key=lambda it: it.item_id), root / "store")
    for s in SPLIT:
        (root / "xml" / f"{s}.xml").write_text(message(s, rois[s]))
        slices, sops, by_z = SIDECAR[s]
        (root / "slices" / f"{s}.json").write_text(json.dumps({"slices": slices, "annotated_sops": sops,
                                                                 "matched_by_z": by_z}))


@pytest.fixture
def data(tmp_path):
    (tmp_path / "xml").mkdir()
    (tmp_path / "slices").mkdir()
    write_store(tmp_path)
    return tmp_path


SCANS = {s: lidc.Scan(PATIENT[s], ZVALS if s != C else np.array([1417.5, 1420.0, 1422.5]), 0.7, 1.25)
         for s in SPLIT}


def report(data, previous=None):
    return quality.z_report(data / "store", data / "xml", data / "slices", SCANS, previous)


def sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_one_entry_per_store_series(data):  # I35: every series of the store, not a list
    r = report(data)
    assert set(r["series"]) == set(SPLIT)
    assert {s: (e["split"], e["patient"], e["items"]) for s, e in r["series"].items()} == {
        s: (SPLIT[s], PATIENT[s], NITEMS[s]) for s in SPLIT}
    assert r["train_sha256"] == sha(data / "store/lidc/train.parquet")


def test_shared_z_lists_both_sop_uids(data):  # I17 (a)
    r = report(data)
    assert r["series"][A]["shared_z"] == {"123.5": ["1.9.1", "1.9.2"]}
    assert all(not r["series"][s]["shared_z"] for s in (B, C, T, V))


def test_roi_beyond_half_the_slice_spacing_is_listed(data):  # I17 (b)
    r = report(data)
    assert r["series"][B]["far_rois"] == [{"sop": "1.9.11", "z": 129.45, "distance_mm": 0.7, "spacing_mm": 1.25},
                                          {"sop": "1.9.5", "z": 130.75, "distance_mm": 2.0, "spacing_mm": 1.25}]
    assert all(not r["series"][s]["far_rois"] for s in (A, C, T, V))


def test_sidecar_hits_and_every_hit_go_to_step_2(data):  # I17 (d)
    r = report(data)
    assert r["series"][T]["sidecar"] == {"slices": 1, "annotated_sops": 2, "matched_by_z": 0}
    assert r["series"][V]["sidecar"] == {"slices": 1, "annotated_sops": 1, "matched_by_z": 1}
    assert all(r["series"][s]["sidecar"] is None for s in (A, B, C))
    assert r["step2"] == sorted([A, B, T, V])
    assert r["series"][C] == {"split": "train", "patient": PATIENT[C], "items": 1, "shared_z": {}, "far_rois": [],
                              "sidecar": None}


def test_0396_check_reports_the_sop_uid_at_1420(data):  # I17 (c), TCIA known issue 6
    assert report(data)["check_0396"] == {"series": C, "sop_uids_at_1420": ["1.9.6"], "all_cited_sops_fetched": True}
    (data / "slices" / f"{C}.json").write_text(json.dumps({"slices": 0, "annotated_sops": 1}))
    assert report(data)["check_0396"]["all_cited_sops_fetched"] is False
    # LIDC-IDRI-0396 has no store item today: its XML is still checked (the series comes from the scan list)
    s = "1.2.840.600001"
    (data / "xml" / "139.xml").write_text(message(s, [(1420.0, "1.9.12"), (1422.5, "1.9.13")]))
    (data / "slices" / f"{s}.json").write_text(json.dumps({"slices": 2, "annotated_sops": 2, "matched_by_z": 0}))
    scans = {**SCANS, C: lidc.Scan("LIDC-IDRI-0001", SCANS[C].zvals, 0.7, 1.25), s: SCANS[C]}
    r = quality.z_report(data / "store", data / "xml", data / "slices", scans)
    assert s not in r["series"]
    assert r["check_0396"] == {"series": s, "sop_uids_at_1420": ["1.9.12"], "all_cited_sops_fetched": True}


def test_val_and_test_hits_are_known_label_issues_not_gate_failures(data):  # I17, I37, I45
    r = report(data)
    assert r["known_label_issues"] == {"val": ids(V), "test": ids(T)}
    full = {i: {"reason": "x", "date": "2026-11-09"} for s in (A, B) for i in ids(s)}
    assert quality.dg17(r, full, r["train_sha256"]) == []


def test_dg17_counts_train_hits_until_every_item_of_the_series_is_excluded(data):  # I17
    r = report(data)
    assert [line.split(":")[0] for line in quality.dg17(r, {}, r["train_sha256"])] == [A, B]
    part = {i: {} for i in ids(A)[:-1] + ids(B)}
    assert quality.dg17(r, part, r["train_sha256"]) == [f"{A}: 1 of 3 items not in train_exclusions.json"]
    assert quality.dg17(r, part | {ids(A)[-1]: {}}, r["train_sha256"]) == []


def test_dg17_fails_a_report_on_another_train_store(data):  # I17: stale report
    r = report(data)
    full = {i: {} for s in (A, B) for i in ids(s)}
    assert quality.dg17(r, full, "0" * 64) == [f"stale report: train store {'0' * 64}, report {r['train_sha256']}"]


def test_mapped_series_is_resolved_only_when_the_rerun_lists_no_hit(data):  # I17 "mapped"
    old = report(data)
    assert old["mapped"] == []
    still = report(data, previous=old)
    assert A in still["step2"] and still["mapped"] == []
    assert any(line.startswith(A) for line in quality.dg17(still, {}, still["train_sha256"]))
    write_store(data, {**ROIS, A: [(123.5, "1.9.1"), (124.75, "1.9.3")]})  # rebuilt: A's contours on one SOP per z
    new = report(data, previous=old)
    assert A not in new["step2"] and new["mapped"] == [A]
    assert quality.dg17(new, {i: {} for i in ids(B)}, new["train_sha256"]) == []


def run(data, *extra):
    return check_series_z.main(["--store", str(data / "store"), "--xml", str(data / "xml"), "--slices",
                                str(data / "slices"), "--exclusions", str(data / "exclusions.json"), "--out",
                                str(data / "out/report.json"), *extra])


@pytest.fixture
def cli(data, monkeypatch):
    monkeypatch.setattr(lidc, "scans", lambda db=None: SCANS)
    import pydicom

    def no_dicom(*a, **k):
        raise AssertionError("a DICOM was opened without --dicom")

    monkeypatch.setattr(pydicom, "dcmread", no_dicom)
    (data / "exclusions.json").write_text(json.dumps({"lidc/x/n000/z000": {"reason": "DG3 a|b", "date": "2026-11-01"}}))
    return data


def test_cli_reports_without_changing_any_input_and_fails_dg17(cli, capsys):  # I17, I37, I45
    before = {p: sha(p) for p in cli.rglob("*") if p.is_file()}
    assert run(cli) == 1
    assert {p: sha(p) for p in before} == before  # store, XML, sidecars, exclusions untouched: the report drops nothing
    r = json.loads((cli / "out/report.json").read_text())
    assert r["step2"] == sorted([A, B, T, V])
    assert r["series"][A]["resolution"] is None and "resolution" not in r["series"][T]
    out = capsys.readouterr().out
    assert "DG17 fail" in out and A in out


def test_cli_drop_writes_every_train_hit_item_to_the_exclusion_list(cli, capsys):  # I17 "dropped"
    assert run(cli, "--drop-train-hits") == 0
    ex = json.loads((cli / "exclusions.json").read_text())
    assert ex["lidc/x/n000/z000"] == {"reason": "DG3 a|b", "date": "2026-11-01"}
    assert {i: e["reason"] for i, e in ex.items() if e["reason"].startswith("DG17")} == {
        **{i: "DG17 100001" for i in ids(A)}, **{i: "DG17 200001" for i in ids(B)}}
    r = json.loads((cli / "out/report.json").read_text())
    assert (r["series"][A]["resolution"], r["series"][B]["resolution"]) == ("dropped", "dropped")
    assert "DG17 pass" in capsys.readouterr().out
    assert run(cli) == 0  # a re-run keeps the drops


def test_dicom_step_lists_duplicate_positions_and_slice_z_without_a_dicom(tmp_path):  # step 2 (fixture DICOM)
    from pydicom.dataset import Dataset, FileMetaDataset
    from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian

    folder = tmp_path / A
    folder.mkdir()
    for k, (sop, z) in enumerate([("1.9.1", 123.5), ("1.9.2", 123.5), ("1.9.3", 124.75)]):
        ds = Dataset()
        ds.SOPInstanceUID, ds.ImagePositionPatient = sop, [0.0, 0.0, z]
        ds.file_meta = FileMetaDataset()
        ds.file_meta.TransferSyntaxUID, ds.file_meta.MediaStorageSOPInstanceUID = ExplicitVRLittleEndian, sop
        ds.file_meta.MediaStorageSOPClassUID = CTImageStorage
        ds.save_as(folder / f"{k}.dcm", enforce_file_format=True)
    assert quality.dicom_z(folder, np.array([123.5, 124.75, 126.0])) == {
        "duplicate_z": {"123.5": ["1.9.1", "1.9.2"]}, "slice_z_without_dicom": [126.0]}
