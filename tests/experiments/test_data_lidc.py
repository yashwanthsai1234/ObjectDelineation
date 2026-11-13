"""LIDC parser: rules P1-P6 on crafted XML, clustering, slice records, patient/z lookup; data-backed parity with the
lab fixture and with pylidc's database contours (sample), the item contract and the B3 regressions on the store."""
import csv
import importlib.util
import json
import random
import sqlite3
from collections import Counter
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import pytest
from hypothesis import given
from hypothesis import strategies as st

from objdel.data import lidc, store
from tests.support.paths import ROOT


def square(x0, y0, x1, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def roi_xml(z, points, inclusion=True):
    edges = "".join(f"<edgeMap><xCoord>{x}</xCoord><yCoord>{y}</yCoord></edgeMap>" for x, y in points)
    return (f"<roi><imageZposition>{z}</imageZposition><imageSOP_UID>s{z}</imageSOP_UID>"
            f"<inclusion>{'TRUE' if inclusion else 'FALSE'}</inclusion>{edges}</roi>")


def nodule_xml(rois, characteristics=True, nid="N"):
    chars = "<characteristics><subtlety>3</subtlety></characteristics>" if characteristics else ""
    rois = "".join(roi_xml(*r) for r in rois)
    return f"<unblindedReadNodule><noduleID>{nid}</noduleID>{chars}{rois}</unblindedReadNodule>"


def lidc_xml(sessions, series="1.2.3", root="LidcReadMessage"):
    body = "".join(f"<readingSession><servicingRadiologistID>anon</servicingRadiologistID>{''.join(s)}</readingSession>"
                   for s in sessions)
    return (f'<?xml version="1.0" encoding="UTF-8"?><{root} xmlns="http://www.nih.gov"><ResponseHeader>'
            f"<SeriesInstanceUid>{series}</SeriesInstanceUid></ResponseHeader>{body}</{root}>")


def write(tmp_path, sessions, name="001.xml", **kw):
    p = tmp_path / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(lidc_xml(sessions, **kw))
    return p


ZVALS = np.array([-12.5, -10.0, -7.5, -5.0, -2.5])


def test_parse_keeps_contoured_nodules_and_drops_single_points(tmp_path):  # P1
    big = nodule_xml([(-10.0, square(2, 2, 6, 6)), (-7.5, [(4, 4)])])  # a 1-point contour inside a >= 3 mm nodule
    small = nodule_xml([(-10.0, [(40, 40)])], characteristics=False)  # a < 3 mm nodule: one point, no contour
    series, sessions = lidc.parse(write(tmp_path, [[big, small], []]))
    assert series == "1.2.3" and len(sessions) == 2 and sessions[1] == []
    assert [[(r.z, r.inclusion, r.points) for r in nod] for nod in sessions[0]] == [
        [(-10.0, True, tuple(square(2, 2, 6, 6)))]]


def test_parse_ignores_other_messages(tmp_path):
    assert lidc.parse(write(tmp_path, [], root="IdriReadMessage")) is None


def test_outline_pixels_belong_to_the_mask():  # P2 (inclusion side)
    m = lidc.draw(square(2, 2, 6, 6))
    assert m.shape == (lidc.GRID, lidc.GRID) and m.sum() == 25 and m[2:7, 2:7].all()
    assert m.view(np.uint8).max() == 1  # canonical bools: True is byte 1, not PIL's 255


def test_hole_survives_into_the_topology_signature(tmp_path):
    from objdel.geom.topology import signature
    ring = nodule_xml([(-10.0, square(0, 0, 10, 10)), (-10.0, square(3, 3, 7, 7), False)])
    _, sessions = lidc.parse(write(tmp_path, [[ring]]))
    assert signature(lidc.rasterize(sessions[0][0], ZVALS)[1], 5) == (1,)


def test_exclusions_remove_their_outline_after_all_inclusions(tmp_path):  # P2 (exclusion side)
    excl = (-10.0, square(3, 3, 7, 7), False)
    nod = nodule_xml([excl, (-10.0, square(0, 0, 10, 10))])  # the exclusion comes first in the file
    _, sessions = lidc.parse(write(tmp_path, [[nod]]))
    masks = lidc.rasterize(sessions[0][0], ZVALS)
    assert list(masks) == [1]
    expect = np.zeros((lidc.GRID, lidc.GRID), bool)
    expect[0:11, 0:11] = True
    expect[3:8, 3:8] = False
    assert np.array_equal(masks[1], expect)


def test_slice_index_from_scan_z_positions():  # P3: annotated z values skip slices
    assert [lidc.slice_index(ZVALS, z) for z in (-10.0, -5.0, -12.5, -2.5)] == [1, 3, 0, 4]


def test_clusters_by_shared_voxels(tmp_path):
    r1 = nodule_xml([(-10.0, square(10, 10, 20, 20)), (-7.5, square(10, 10, 20, 20))])
    r2 = nodule_xml([(-7.5, square(15, 15, 25, 25))])
    r3 = nodule_xml([(-5.0, square(100, 100, 110, 110))])
    r3b = nodule_xml([(-12.5, square(200, 200, 210, 210))])
    _, sessions = lidc.parse(write(tmp_path, [[r1], [r2], [r3, r3b]]))
    cl = lidc.clusters(sessions, ZVALS)
    assert [(sorted(c), dup) for c, dup in cl] == [([2], False), ([0, 1], False), ([2], False)]  # by first slice
    assert sorted(cl[1][0][0]) == [1, 2] and sorted(cl[1][0][1]) == [2]


def test_two_annotations_of_one_reader_merge_and_are_flagged(tmp_path):
    a = nodule_xml([(-10.0, square(10, 10, 20, 20))])
    b = nodule_xml([(-10.0, square(18, 18, 30, 30))])  # same reader, overlapping: one nodule
    _, sessions = lidc.parse(write(tmp_path, [[a, b], [a]]))
    ((per_reader, dup),) = lidc.clusters(sessions, ZVALS)
    merged = lidc.draw(square(10, 10, 20, 20)) | lidc.draw(square(18, 18, 30, 30))
    assert dup and np.array_equal(per_reader[0][1], merged)
    assert lidc.slices(sessions, ZVALS)[0].dup_reader


def test_slices_are_multi_reader_crops(tmp_path):
    r1 = nodule_xml([(-10.0, square(10, 10, 20, 20)), (-7.5, square(10, 10, 20, 20))])
    r2 = nodule_xml([(-7.5, square(15, 15, 25, 25))])
    _, sessions = lidc.parse(write(tmp_path, [[r1], [r2], []]))
    (s,) = lidc.slices(sessions, ZVALS)  # slice k=1 has one reader only
    assert (s.cluster, s.k, s.readers, s.y0, s.x0) == (0, 2, [1, 2], 10 - lidc.MARGIN, 10 - lidc.MARGIN)
    assert s.stack.shape == (2, 16 + 2 * lidc.MARGIN, 16 + 2 * lidc.MARGIN)
    assert s.stack[0].sum() == 121 and s.stack[1].sum() == 121 and (s.stack[0] & s.stack[1]).sum() == 36
    full = s.full()
    assert full.shape == (2, lidc.GRID, lidc.GRID) and full[0, 10:21, 10:21].all() and full.sum() == 242


def test_series_files_prefers_correction(tmp_path):
    a = write(tmp_path, [], name="x/157/161.xml", series="9.9")
    c = write(tmp_path, [], name="161-resubmitted-correction.xml", series="9.9")
    b1 = write(tmp_path, [], name="y/188/082.xml", series="8.8")
    write(tmp_path, [], name="y/188/083.xml", series="8.8")
    write(tmp_path, [], name="z/cxr.xml", root="IdriReadMessage")
    assert lidc.series_files(tmp_path) == {"9.9": c, "8.8": b1} and a.exists()


def test_scans_from_pylidc_database(tmp_path):
    db = tmp_path / "pylidc.sqlite"
    con = sqlite3.connect(db)
    con.execute("create table scans (id integer, series_instance_uid varchar, patient_id varchar, "
                "pixel_spacing float, slice_thickness float)")
    con.execute("create table zvals (id integer, scan_id integer, val float)")
    con.execute("insert into scans values (1, '1.2.3', 'LIDC-IDRI-0001', 0.7, 2.5)")
    con.executemany("insert into zvals values (?, 1, ?)", [(1, -5.0), (2, -12.5), (3, -7.5)])
    con.commit()
    con.close()
    s = lidc.scans(db)["1.2.3"]
    assert s.patient == "LIDC-IDRI-0001" and s.zvals.tolist() == [-12.5, -7.5, -5.0]
    assert (s.pixel_spacing, s.slice_thickness) == (0.7, 2.5)


def test_ct_slices_matched_by_z(tmp_path):
    a, b = np.full((4, 4), -1000, np.int16), np.full((4, 4), 40, np.int16)
    np.savez(tmp_path / "1.2.3.npz", **{"z-10.00": a, "z-5.00": b, "z3.00": b})
    got = lidc.ct_slices(tmp_path / "1.2.3.npz", ZVALS)
    assert sorted(got) == [1, 3] and np.array_equal(got[1], a) and got[3].dtype == np.int16


def test_sop_z_maps_every_roi_of_every_file(tmp_path):
    write(tmp_path, [[nodule_xml([(-10.0, square(1, 1, 5, 5)), (-7.504, [(3, 3)])])]], name="a/1.xml")
    write(tmp_path, [[nodule_xml([(-5.0, square(1, 1, 5, 5))], characteristics=False)]], name="b/2.xml")
    write(tmp_path, [], name="c/cxr.xml", root="IdriReadMessage")
    assert lidc.sop_z(tmp_path) == {"1.2.3": {"s-10.0": -10.0, "s-7.504": -7.5, "s-5.0": -5.0}}


def test_fetch_matches_slices_by_sop_uid_then_by_z():  # B3-4 (LIDC-IDRI-0365: no cited SOP UID exists on TCIA)
    pytest.importorskip("pydicom")
    from objdel.cli.data.fetch_lidc import choose
    headers = [("a.dcm", "u1", -10.0), ("b.dcm", "u2", -7.5), ("c.dcm", "u3", -5.004), ("d.dcm", "u4", -2.5),
               ("e.dcm", "u5", -2.5), ("f.dcm", "u6", 0.0)]
    sops = {"u1": -10.0, "x1": -10.0, "x2": -5.0, "x3": -2.5, "x4": -1.0}  # x*: UIDs absent from the series
    keep, by_z = choose(headers, sops)
    # -10 by UID (x1 at the same z adds nothing); -5 by z within 0.01 mm; -2.5 ambiguous (two DICOMs); -1 absent
    assert keep == {"a.dcm": -10.0, "c.dcm": -5.0} and by_z == 1
    assert choose(headers, {"u2": -7.5}) == ({"b.dcm": -7.5}, 0)


REAL_UID = "1.3.6.1.4.1.14519.5.2.1.6279.6001.100225287222365663678666836860"  # 64 characters, the LIDC maximum


@pytest.mark.parametrize("bad", [
    "../../../victim_dir", "/Users/someone/work", "1.2.3/../../x", "1.2.3?x=1&y=2", "1.2.3%2F..", "1.2.3\n", " 1.2.3",
    "", "1", "1.", ".1.2", "1..2", "1.02.3", "01.2", "1.2.3\x00", "١.٢",  # Arabic-Indic digits match \d
    REAL_UID + "1", ".".join(["1"] * 33)])  # 65 characters
def test_uid_refuses_everything_outside_the_dicom_grammar(bad):  # review #1, DICOM PS3.5 §9.1
    with pytest.raises(ValueError, match="DICOM UID"):
        lidc.uid(bad)


@pytest.mark.parametrize("good", ["1.2.3", "0.0", "1.2.840.10008.1.2.1", REAL_UID])
def test_uid_accepts_dicom_uids(good):
    assert lidc.uid(good) == good


@given(st.from_regex(r"[0-9./%?&\\a-z\x00 -]{0,70}", fullmatch=True) | st.text(max_size=70))
def test_an_accepted_uid_is_one_plain_path_component_and_url_safe(s):
    from urllib.parse import quote
    try:
        lidc.uid(s)
    except ValueError:
        return
    assert Path("/root", s).parent == Path("/root") and s not in (".", "..") and quote(s, safe="") == s


@pytest.mark.parametrize("planted", ["../../victim", "VICTIM"])
def test_planted_series_uid_is_refused_by_every_xml_reader(tmp_path, planted):  # review #1
    series = planted if planted != "VICTIM" else str(tmp_path / "victim")  # relative traversal / absolute path
    p = write(tmp_path / "xml", [[nodule_xml([(-7.5, square(1, 1, 5, 5))])]], name="bad.xml", series=series)
    for read in (lambda: lidc.sop_z(tmp_path / "xml"), lambda: lidc.series_files(tmp_path / "xml"),
                 lambda: lidc.parse(p)):
        with pytest.raises(ValueError, match="bad.xml.*DICOM UID"):
            read()


@pytest.mark.parametrize("planted", ["../../victim", "VICTIM"])
def test_fetch_never_touches_a_path_named_by_a_planted_uid(tmp_path, monkeypatch, planted):  # review #1, no network
    pytest.importorskip("pydicom")
    from objdel.cli.data import fetch_lidc as F
    victim = tmp_path / "victim"
    victim.mkdir()
    (victim / "keep.txt").write_text("x")
    out = tmp_path / "data" / "lidc_slices"  # tmp root: data/lidc_tmp, so data/lidc_tmp/../../victim is the victim
    series = planted if planted != "VICTIM" else str(victim)
    assert (out.parent / "lidc_tmp" / series).resolve() == victim.resolve()
    write(tmp_path / "xml", [[nodule_xml([(-7.5, square(1, 1, 5, 5))])]], series=series)
    calls = []
    monkeypatch.setattr(F.urllib.request, "urlopen", lambda *a, **k: calls.append(a) or 1 / 0)
    with pytest.raises(ValueError, match="DICOM UID"):
        F.main([str(tmp_path / "xml"), str(out)])
    assert (victim / "keep.txt").read_text() == "x" and calls == [] and not (out.parent / "lidc_tmp").exists()


def test_fetch_cli_series_must_be_dicom_uids(tmp_path, monkeypatch, capsys):  # review #1: the other way in
    pytest.importorskip("pydicom")
    from objdel.cli.data import fetch_lidc as F
    monkeypatch.setattr(F.urllib.request, "urlopen", lambda *a, **k: 1 / 0)
    write(tmp_path / "xml", [[nodule_xml([(-7.5, square(1, 1, 5, 5))])]])
    with pytest.raises(SystemExit):
        F.main([str(tmp_path / "xml"), str(tmp_path / "out"), "--series", "../../x"])
    assert "DICOM UID" in capsys.readouterr().err


def test_fetch_url_quotes_the_uid(tmp_path, monkeypatch):  # review #1: the UID goes into the query as one value
    pytest.importorskip("pydicom")
    from urllib.parse import parse_qs, urlsplit

    from objdel.cli.data import fetch_lidc as F
    urls = []
    monkeypatch.setattr(F.urllib.request, "urlopen", lambda url, **k: urls.append(url) or 1 / 0)
    monkeypatch.setattr(F.time, "sleep", lambda s: None)
    with pytest.raises(ZeroDivisionError):
        F.fetch(REAL_UID, {}, tmp_path / "out", tmp_path / "tmp")
    assert len(urls) == 4 and parse_qs(urlsplit(urls[0]).query) == {"SeriesInstanceUID": [REAL_UID]}


def test_item_format(tmp_path):
    r1 = nodule_xml([(-7.5, square(10, 10, 20, 20))])
    r2 = nodule_xml([(-7.5, square(12, 12, 20, 20))])
    _, sessions = lidc.parse(write(tmp_path, [[r1], [r2]]))
    (s,) = lidc.slices(sessions, ZVALS)
    scan = lidc.Scan("LIDC-IDRI-0007", ZVALS, 0.7, 2.5)
    it = lidc.item("1.2.3", scan, s, "val", None, sessions=2, xml="157/161.xml")
    assert (it.item_id, it.dataset, it.group_id, it.split, it.raters) == (
        "lidc/1.2.3/n000/z002", "lidc", "LIDC-IDRI-0007", "val", ["1", "2"])
    assert it.masks["nodule"].shape == (2, 512, 512) and it.image is None
    assert it.meta == {"series": "1.2.3", "patient": "LIDC-IDRI-0007", "cluster": 0, "slice": 2, "z_mm": -7.5,
                       "readers_in_cluster": ["1", "2"], "dup_reader_in_cluster": False, "sessions": 2,
                       "xml": "157/161.xml", "pixel_spacing_mm": 0.7, "slice_thickness_mm": 2.5,
                       "copied_raters": {}, "other_nodules": []}


def test_identical_layers_are_one_rater(tmp_path):  # B3-3: a copied outline is not a second opinion
    a, b = square(10, 10, 20, 20), square(12, 12, 22, 22)
    r1 = nodule_xml([(-7.5, a), (-5.0, a)])
    r2 = nodule_xml([(-7.5, b), (-5.0, b)])
    r3 = nodule_xml([(-7.5, a), (-5.0, square(11, 10, 20, 20))])  # copies r1 on z -7.5 only
    _, sessions = lidc.parse(write(tmp_path, [[r1], [r2], [r3]]))
    s2, s3 = lidc.slices(sessions, ZVALS)
    assert (s2.k, s2.readers, s2.copied) == (2, [1, 2], {3: 1}) and len(s2.stack) == 2
    assert (s3.k, s3.readers, s3.copied) == (3, [1, 2, 3], {})
    assert s2.cluster_readers == [1, 2, 3]
    it = lidc.item("1.2.3", lidc.Scan("P", ZVALS, 0.7, 2.5), s2, "train", None, sessions=3, xml="x.xml")
    assert it.raters == ["1", "2"] and it.meta["copied_raters"] == {"3": "1"}


def test_a_slice_left_with_one_reader_after_copies_is_dropped(tmp_path):  # B3-3
    a = square(10, 10, 20, 20)
    _, sessions = lidc.parse(write(tmp_path, [[nodule_xml([(-7.5, a)])], [nodule_xml([(-7.5, a)])]]))
    assert lidc.slices(sessions, ZVALS) == []


def test_other_nodules_on_the_slice_are_listed(tmp_path):  # B3-2
    near = square(40, 10, 45, 14)  # another nodule, outlined by one reader only
    r1 = nodule_xml([(-7.5, square(10, 10, 20, 20))])
    r2 = nodule_xml([(-7.5, square(12, 12, 20, 20))])
    far, far2 = (nodule_xml([(-7.5, square(x, 300, 310, 305))], nid="far") for x in (300, 301))
    _, sessions = lidc.parse(write(tmp_path, [[r1, nodule_xml([(-7.5, near)]), far], [r2, far2],
                                              [nodule_xml([(-10.0, square(40, 10, 45, 14))])]]))
    (s, s_far) = lidc.slices(sessions, ZVALS)
    assert s.others == [(40, 10, 46, 15), (300, 300, 311, 306)]  # half-open boxes; the z -10 nodule is elsewhere
    assert s_far.others == [(10, 10, 21, 21), (40, 10, 46, 15)]
    assert lidc.crowded(s, (0, 0, 40, 40)) is False  # half-open: x 40 is outside
    assert lidc.crowded(s, (0, 0, 41, 40)) is True
    assert lidc.crowded(s, (290, 290, 300, 300)) is False and lidc.crowded(s, (290, 290, 301, 301)) is True


# ---- data-backed -----------------------------------------------------------------------------------------------

@pytest.mark.data
def test_store_lidc_equals_lab_fixture(store_dir, fixtures_dir):
    """Every stored slice, bit for bit, against the lab's full fixture (fixed parser, pylidc z values), minus the
    wave-6 rules: a copied layer (P5) is gone and equals its kept twin; fixture slices left out are those with fewer
    than two distinct layers (P5) or whose frame showed another nodule (P6, counted in build_report)."""
    fx = np.load(fixtures_dir / "lidc_full.npz")
    meta = {r["key"]: r for r in csv.DictReader(open(fixtures_dir / "lidc_full_meta.csv"))}
    seen = set()
    for split in ("train", "val", "test"):
        for it in store.load("lidc", split, store_dir):
            key = f"{it.meta['series'][-12:]}_c{it.meta['cluster']}_k{it.meta['slice']}__nodule"
            r = meta[key]
            y0, x0, h, w = (int(r[c]) for c in ("crop_y0", "crop_x0", "h", "w"))
            readers, layers = r["reader_ids"].split(";"), fx[key].astype(bool)
            copied = it.meta["copied_raters"]
            assert it.raters == [e for e in readers if e not in copied] and it.meta["patient"] == r["patient"]
            assert it.meta["dup_reader_in_cluster"] == (r["dup_reader_in_cluster"] == "True")
            keep = [readers.index(e) for e in it.raters]
            assert np.array_equal(it.masks["nodule"][:, y0:y0 + h, x0:x0 + w], layers[keep]), key
            assert it.masks["nodule"].sum() == layers[keep].sum(), key  # nothing outside the lab crop
            for e, twin in copied.items():
                assert np.array_equal(layers[readers.index(e)], layers[readers.index(twin)]), key
            seen.add(key)
    rest = set(meta) - seen
    distinct = {k: len({fx[k][i].tobytes() for i in range(len(fx[k]))}) for k in rest}
    crowded = sum(json.loads((store_dir / "build_report.json").read_text())["lidc"]["dropped_crowded_frame"].values())
    assert len(meta) == 10921 and sum(d >= 2 for d in distinct.values()) == crowded
    assert sum(d < 2 for d in distinct.values()) == 80  # P5: slices left with one distinct reader


@pytest.mark.data
def test_store_lidc_val_rows_keep_the_item_contract(store_dir):
    """Val split, every row: masks (k, 512, 512), k = len(raters) >= 2 unique raters from the cluster's readers, no
    empty layer and (B3-3) no two identical layers, int16 512 x 512 HU, group_id = patient = pylidc patient, z_mm on
    the slice; B3-2: the model frame overlaps no other outlined nodule of the slice."""
    from objdel.data.targets.frame import frame
    scans, n = lidc.scans(), 0
    for it in store.load("lidc", "val", store_dir):
        m, k, meta = it.masks["nodule"], len(it.raters), it.meta
        assert m.shape == (k, 512, 512) and k >= 2 and len(set(it.raters)) == k and m.any((1, 2)).all()
        assert len({x.tobytes() for x in m}) == k and set(it.raters) <= set(meta["readers_in_cluster"])
        sc = scans[meta["series"]]
        assert it.group_id == meta["patient"] == sc.patient and sc.zvals[meta["slice"]] == meta["z_mm"]
        assert it.split == "val" and it.image.dtype == np.int16 and it.image.shape == (512, 512)
        x0, y0, x1, y1 = frame(it).box
        assert not any(a < x1 and x0 < c and b < y1 and y0 < d for a, b, c, d in it.meta["other_nodules"])
        n += 1
    assert n > 1000


@pytest.mark.data
def test_store_lidc_copied_reader_of_0811_is_one_rater(store_dir):  # B3-3 regression on the break-it item
    iid = "lidc/1.3.6.1.4.1.14519.5.2.1.6279.6001.220596530836092324070084384692/n000/z243"
    (r,) = pq.read_table(store_dir / "lidc" / "train.parquet", filters=[("item_id", "==", iid)]).to_pylist()
    assert "3" not in r["raters"] and json.loads(r["meta"])["copied_raters"] == {"3": "1"}


@pytest.mark.data
def test_store_lidc_crowded_break_it_item_is_left_out(store_dir):  # B3-2 regression: LIDC-IDRI-0942 n000/z092
    series = "1.3.6.1.4.1.14519.5.2.1.6279.6001.137763212752154081977261297097"
    ids = set(pq.read_table(store_dir / "lidc" / "train.parquet", columns=["item_id"]).column(0).to_pylist())
    assert f"lidc/{series}/n000/z092" not in ids


@pytest.mark.data
def test_store_lidc_every_item_has_pixels(store_dir):  # B3-4: z fallback of fetch_lidc_pixels (LIDC-IDRI-0365)
    rep = json.loads((store_dir / "build_report.json").read_text())["lidc"]
    assert all(rep[f"{s}/with_image"] == rep[f"{s}/items"] for s in ("train", "val", "test"))


def pylidc_rule_mask(contours):
    """Oracle: pylidc's own rasterizer (matplotlib Path.contains_points on pixel centres, as in
    pylidc Annotation.boolean_mask) with rule P2: outline points added for inclusions, removed with exclusions."""
    from matplotlib.path import Path as MPath
    yy, xx = np.mgrid[:lidc.GRID, :lidc.GRID]
    pts = np.c_[xx.ravel(), yy.ravel()]
    m = np.zeros((lidc.GRID, lidc.GRID), bool)
    for inclusion in (True, False):
        for inc, xy in contours:
            if inc != inclusion:
                continue
            xy = np.asarray(xy)
            inside = MPath(np.vstack([xy, xy[:1]]), closed=True).contains_points(pts).reshape(m.shape)
            inside[xy[:, 1], xy[:, 0]] = True
            m = m | inside if inclusion else m & ~inside
    return m


@pytest.mark.data
def test_parity_with_pylidc_database_sample():
    """Lab: 99.97 % of reader-slice masks exact over all 1,018 scans. Here: 25 random scans, per annotation slice."""
    pytest.importorskip("matplotlib")
    xml_root = ROOT / "data" / "raw" / "lidc_xml"
    if not xml_root.is_dir():
        pytest.skip("LIDC XML absent")
    if importlib.util.find_spec("pylidc") is None:  # BK4-20: a [dev]-only env with data/; pylidc's import is broken
        pytest.skip("pylidc not installed")
    db = lidc.pylidc_db()
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    files = lidc.series_files(xml_root)
    sample = random.Random(0).sample(sorted(files), 25)
    exact = total = 0
    for series in sample:
        scan = lidc.scans(db)[series]
        _, sessions = lidc.parse(files[series])
        ours = {}
        for nods in sessions:
            for nod in nods:
                key = frozenset((round(r.z, 2), r.inclusion, r.points) for r in nod)
                ours.setdefault(key, []).append(nod)
        rows = con.execute(
            "select a.id, c.inclusion, c.image_z_position, c.coords from annotations a "
            "join scans s on a.scan_id = s.id join contours c on c.annotation_id = a.id "
            "where s.series_instance_uid = ?", (series,)).fetchall()
        theirs = {}
        for aid, inc, z, coords in rows:
            pts = tuple(tuple(int(v) for v in p.split(",")) for p in coords.split("\n") if p)
            theirs.setdefault(aid, []).append((round(z, 2), bool(inc), pts))
        assert Counter({k: len(v) for k, v in ours.items()}) == Counter(
            frozenset(c for c in cs if len(c[2]) > 1) for cs in theirs.values()), series
        for key, nods in ours.items():
            got = lidc.rasterize(nods[0], scan.zvals)
            for k in got:
                z = round(float(scan.zvals[k]), 2)
                want = pylidc_rule_mask([(inc, pts) for zz, inc, pts in key if zz == z])
                total += 1
                exact += bool(np.array_equal(got[k], want))
    assert total > 100 and exact / total >= 0.995, (exact, total)
