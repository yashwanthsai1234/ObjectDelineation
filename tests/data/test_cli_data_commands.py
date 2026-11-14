"""The data commands in process (objdel.cli.data.targets, prepare, probes, store, views, fetch_lidc) on synthetic
stores and crafted XML: what each writes and prints, and what it refuses. Process pools run inline
(tests.support.proc.inline_context), so the worker functions run here; model tokenizers are a character tokenizer
(one token per character), so no model files or network are needed."""
import json
import sys
from dataclasses import replace
from functools import partial

import numpy as np
import pytest

from objdel.core import datasets
from objdel.data import store
from tests.support.proc import inline_context


def char_tokenizer(text, add_special_tokens=False):
    """One token per character (the shape of a transformers tokenizer call)."""
    return {"input_ids": [0] * len(text)}


def _argv(monkeypatch, prog, *args):
    monkeypatch.setattr(sys, "argv", [prog, *map(str, args)])


def _jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


@pytest.fixture
def char_tokens(monkeypatch):
    """objdel.policy.models.tokenizer -> a character tokenizer for every model (budgets unchanged)."""
    from objdel.policy import models

    monkeypatch.setattr(models, "tokenizer", lambda m: (char_tokenizer, models.MODELS[m].max_completion))
    return models.MODELS


@pytest.fixture
def train_val_store(tmp_path, synthetic_items):
    store.write(synthetic_items(("train", "val")), tmp_path / "store")
    return tmp_path / "store"


# ------------------------------------------------------------------------------------------------------ targets
def test_targets_model_tokens_counts_every_model(monkeypatch, char_tokens):
    from objdel.cli.data import targets

    monkeypatch.setattr(targets, "_TOK", {})
    assert targets._model_tokens('{"nodule": []}') == {k: 14 for k in char_tokens}


def test_targets_main_writes_the_library_targets(tmp_path, monkeypatch, capsys, char_tokens, train_val_store):
    pytest.importorskip("transformers")  # a pinned version the manifest records
    from objdel.cli.data import targets
    from objdel.data.targets.build import BuildConfig, rater_targets, target_set_hash
    from objdel.data.targets.frame import frame
    from objdel.geom.text import parse

    monkeypatch.setattr(targets, "_TOK", {})
    monkeypatch.setattr(targets, "get_context", inline_context)
    out = tmp_path / "targets"
    _argv(monkeypatch, "targets", "--store", train_val_store, "--out", out, "--procs", 1)
    targets.main()

    man = json.loads((out / "manifest.json").read_text())
    summary = json.loads((out / "summary.json").read_text())
    every = []
    for ds in datasets.DATASETS:
        names = datasets.get(ds).names
        for split in ("train", "val"):
            want, frames = [], {}
            for it in store.load(ds, split, train_val_store):
                f = frame(it)
                frames[it.item_id] = f
                _, built = rater_targets(f.masks, names, len(it.raters), BuildConfig())
                want += [(it.item_id, r, t.text) for r, t in zip(it.raters, built, strict=True)
                         if not isinstance(t, ValueError)]
            rows = _jsonl(out / ds / f"{split}.jsonl")
            assert [(r["item_id"], r["rater"], r["text"]) for r in rows] == want and want
            for r in rows:
                f = frames[r["item_id"]]
                assert (r["W"], r["H"], r["box"], r["scale"]) == (f.W, f.H, [int(v) for v in f.box], f.scale)
                assert parse(r["text"], names).ok
            assert man["files"][f"{ds}/{split}.jsonl"] == {"targets": len(want), "sha256": target_set_hash(want)}
            s = summary[f"{ds}/{split}"]
            assert s["items"] == len(frames) and s["targets"] == len(want) and s["errors"] == [] and s["parse"] == 1.0
            # character tokenizer: the answer length is the text length
            assert all(s["answer_tokens"][k]["100"] == max(len(t) for *_, t in want) for k in char_tokens)
            assert set(s["objects"]) == set(names)
            assert all(o["n"] == len(want) and o["legal"] == 1.0 and o["iou_min"] > 0.9 for o in s["objects"].values())
            every += want
    assert man["target_set_sha256"] == summary["target_set_sha256"] == target_set_hash(every)
    assert man["store_manifest_sha256"] and set(man["versions"]) >= set(targets.PINS)
    assert '"target_set_sha256"' in capsys.readouterr().out


def test_targets_summary_keeps_errors_out_of_the_statistics():
    from objdel.cli.data import targets

    ok = {"item_id": "a", "rater": "1", "parse": True, "time": 0.5, "tokens": {k: 10 for k in targets.MODELS},
          "obj": {"nodule": {"k5": True, "raw": False}}}
    s = targets.summarize([ok, {"item_id": "a", "rater": "2", "error": "MAX_RINGS"}],
                          [{"item_id": "a", "W": 64, "scale": 1}])
    assert s["errors"] == [("a", "2", "MAX_RINGS")] and s["targets"] == 2 and s["frame_side"] == {64: 1}
    assert s["build_s"] == {"median": 0.5, "p95": 0.5, "max": 0.5, "total": 0.5}
    assert s["objects"] == {"nodule": {"n": 1, "K5": 1.0, "raw_topology": 0.0}}  # no "iou": the short form
    assert targets._pct([], 50) is None


# ------------------------------------------------------------------------------------------------------ prepare
def test_prepare_main_prints_sets_and_the_token_budget(tmp_path, monkeypatch, capsys, char_tokens, train_val_store):
    from objdel.cli.data import prepare
    from objdel.data import prepare as data_prepare
    from objdel.policy import rows as policy_rows

    monkeypatch.setattr(data_prepare, "get_context", inline_context)
    out = tmp_path / "prepared"
    _argv(monkeypatch, "prepare", "--store", train_val_store, "--out", out, "--procs", 1)
    prepare.main()
    man = json.loads((out / "manifest.json").read_text())
    assert set(man["sets"]) == {f"{d}/{s}" for d in datasets.DATASETS for s in ("train", "val")}
    printed = capsys.readouterr().out
    assert json.dumps(man["sets"], indent=1) in printed and json.dumps(man["files"], indent=1) in printed

    # the I5 budget (answer tokens + 1): a budget at the shortest LIDC target puts the longer ones over it
    from objdel.policy import models
    lengths = {ds: [len(t) + 1 for r in policy_rows.records(out, ds, "train") for t in r["targets"] if t]
               for ds in datasets.DATASETS}
    budget = min(lengths["lidc"])
    monkeypatch.setattr(models, "MODELS", {m: replace(e, max_completion=budget) for m, e in char_tokens.items()})
    got = prepare.token_budget(out)
    assert set(got) == {f"{m}/{d}" for m in char_tokens for d in datasets.DATASETS}
    for key, v in got.items():
        n = lengths[key.split("/")[1]]
        assert v == {"targets": len(n), "budget": budget, "over": sum(x > budget for x in n), "max": max(n),
                     "over_lengths": sorted(x for x in n if x > budget)}
    assert all(0 < v["over"] < v["targets"] for k, v in got.items() if k.endswith("/lidc"))


# ------------------------------------------------------------------------------------------------------- probes
def _holes():
    """The probe tests' hole records with the filled-part fraction that objdel.data.probes.geometry also records."""
    from tests.support.probes import HOLES

    return [h | {"fraction": h["area"] / (np.pi * h["part_diam"] ** 2 / 4)} for h in HOLES]


@pytest.fixture
def small_probes(monkeypatch):
    """objdel.data.probes with synthetic stats and patches and small topology / fidelity draws; returns the calls."""
    from objdel.data import probes as P
    from tests.support.probes import MULTIS, diameters, lung_patches

    calls = []
    monkeypatch.setattr(P, "lidc_stats", lambda root: calls.append(("stats", str(root))) or
                        {"diams": diameters(), "holes": _holes(), "multis": MULTIS, "prior": {}})
    monkeypatch.setattr(P, "patches", lambda root, split, n: calls.append(("patches", split, n)) or lung_patches())
    monkeypatch.setattr(P, "build_topology", partial(P.build_topology, per_bin=3, n_two=3, n_simple=3))
    monkeypatch.setattr(P, "fidelity", partial(P.fidelity, n_holes=40, n_two=20))
    return calls


def test_probes_build_writes_both_probes_and_the_fidelity_check(tmp_path, monkeypatch, capsys, small_probes):
    from objdel.cli.data import probes
    from objdel.data import probes as P
    from tests.support.probes import MULTIS

    _argv(monkeypatch, "probes", "build", "--store", tmp_path / "s", "--out", tmp_path / "ev", "--n-test", 12,
          "--n-train", 7)
    probes.main()
    assert small_probes == [("stats", str(tmp_path / "s")), ("patches", "val", 200), ("patches", "train", 400),
                            ("patches", "val", 200)]
    assert len(P.records(tmp_path / "ev" / "probe-precision")) == 12
    assert len(P.records(tmp_path / "ev" / "probe-precision", "train")) == 7
    assert len(P.records(tmp_path / "ev" / "probe-topology")) == 3 * 3 + 3 + 3
    man = json.loads((tmp_path / "ev" / "probe-topology" / "manifest.json").read_text())
    assert man["sets"]["lidc/val"]["fidelity_seed_2027"] == P.fidelity(_holes(), MULTIS, 2027)
    out = capsys.readouterr().out
    assert "stats: 4000 diameters, 9 holes, 3 multi-part" in out and "precision train: 7 items" in out


def test_probes_build_only_topology_and_score(tmp_path, monkeypatch, capsys, small_probes):
    from objdel.cli.data import probes
    from objdel.data import probes as P

    _argv(monkeypatch, "probes", "build", "--out", tmp_path / "ev", "--only", "topology")
    probes.main()
    assert not (tmp_path / "ev" / "probe-precision").exists() and [c[1] for c in small_probes[1:]] == ["val"]
    root = tmp_path / "ev" / "probe-topology"
    rng = np.random.default_rng(0)
    rows = [{"item_id": r["item_id"], "samples": [P.model_centroid_disc(r, rng)]} for r in P.records(root)]
    (tmp_path / "c.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
    capsys.readouterr()
    _argv(monkeypatch, "probes", "score", "--probe", "topology", "--completions", tmp_path / "c.jsonl", "--root", root,
          "--out", tmp_path / "score.json")
    probes.main()
    want = json.loads(json.dumps(P.score_topology(rows, root), default=float))
    assert json.loads((tmp_path / "score.json").read_text()) == want
    assert json.loads(capsys.readouterr().out) == want


def test_probes_cli_refuses_an_unknown_probe(monkeypatch):
    from objdel.cli.data import probes

    _argv(monkeypatch, "probes", "score", "--probe", "shape", "--completions", "c.jsonl")
    with pytest.raises(SystemExit) as e:
        probes.main()
    assert e.value.code == 2


# -------------------------------------------------------------------------------------------------------- store
ZVALS = np.array([-12.5, -10.0, -7.5, -5.0, -2.5])
N_SERIES = 1018  # the command refuses any other series count


@pytest.fixture
def lidc_repo(tmp_path, monkeypatch):
    """A repo root with two crafted series under data/raw/lidc_xml, CT pixels of the first, and a lab split; the
    pylidc database is 1018 scans (the two plus 1016 without nodules, all on one empty XML). Returns (repo, ct)."""
    from objdel.cli.data import store as cli
    from objdel.data import lidc
    from tests.support.lidc_xml import nodule_xml, square, write

    repo = tmp_path / "repo"
    xml = repo / "data" / "raw" / "lidc_xml"
    # 1.2.1: nodule A by readers 1 and 2 on z -10 and -7.5; reader 1's nodule B next to A on z -10 crowds that slice
    a1 = nodule_xml([(-10.0, square(100, 100, 120, 120)), (-7.5, square(100, 100, 120, 120))], nid="A")
    a2 = nodule_xml([(-10.0, square(102, 102, 121, 121)), (-7.5, square(101, 101, 119, 121))], nid="A")
    write(xml, [[a1, nodule_xml([(-10.0, square(126, 100, 136, 120))], nid="B")], [a2]], name="a/1.xml",
          series="1.2.1")
    # 1.2.2: a ring, a solid outline and a copy of the ring on z -5: topologies disagree, reader 3 is a copy
    ring = nodule_xml([(-5.0, square(200, 200, 230, 230)), (-5.0, square(210, 210, 220, 220), False)])
    write(xml, [[ring], [nodule_xml([(-5.0, square(201, 201, 230, 231))])], [ring]], name="b/2.xml", series="1.2.2")
    empty = write(tmp_path, [[]], name="empty.xml", series="1.3.0")
    pad = {f"1.3.{i}": empty for i in range(N_SERIES - 2)}
    scans = {s: lidc.Scan(f"LIDC-IDRI-{i:04d}", ZVALS, 0.7, 2.5) for i, s in enumerate(["1.2.1", "1.2.2", *pad])}
    real = lidc.series_files
    monkeypatch.setattr(lidc, "series_files", lambda root: real(root) | pad)
    monkeypatch.setattr(lidc, "scans", lambda: scans)
    ct = np.random.default_rng(0).integers(-1000, 400, (512, 512), dtype=np.int16)
    (repo / "data" / "lidc_slices").mkdir(parents=True)
    np.savez(repo / "data" / "lidc_slices" / "1.2.1.npz", **{"z-7.50": ct})
    (repo / "lab.csv").write_text("patient,split\nLIDC-IDRI-0000,val\nLIDC-IDRI-0001,train\n")
    monkeypatch.setattr(cli, "REPO", repo)
    monkeypatch.setattr(cli, "LAB_SPLIT", repo / "lab.csv")
    monkeypatch.setattr(cli, "get_context", inline_context)
    return repo, ct


def test_store_main_builds_items_splits_and_report(lidc_repo, capsys):
    from objdel.cli.data import store as cli

    repo, ct = lidc_repo
    out = repo / "data" / "store"
    out.mkdir(parents=True)
    (out / "stale.txt").write_text("from an old build")  # the command replaces the whole directory
    cli.main(["--out", str(out), "--workers", "1"])
    assert not (out / "stale.txt").exists()
    rep = json.loads((out / "build_report.json").read_text())
    items = {it.item_id: it for rel in rep["manifest"]["files"]
             for it in store.load("lidc", rel.split("/")[1].split(".")[0], out)}
    assert sorted(items) == ["lidc/1.2.1/n000/z002", "lidc/1.2.2/n000/z003"]
    a, b = items["lidc/1.2.1/n000/z002"], items["lidc/1.2.2/n000/z003"]
    assert np.array_equal(a.image, ct) and b.image is None
    assert a.group_id == "LIDC-IDRI-0000" and a.meta["xml"] == "a/1.xml" and a.meta["sessions"] == 2
    assert b.raters == ["1", "2"] and b.meta["copied_raters"] == {"3": "1"}
    lid = rep["lidc"]
    assert (lid["series"], lid["patients"], lid["items"]) == (N_SERIES, N_SERIES, 2)
    assert sum(lid["dropped_crowded_frame"].values()) == 1 and sum(lid["patients_per_split"].values()) == N_SERIES
    flags = lambda key: sum(v for k, v in lid.items() if k.endswith("/" + key))  # noqa: E731
    assert (flags("items"), flags("with_image"), flags("disagree_5px"), flags("copied_layers")) == (2, 1, 1, 1)
    # the lab split puts only 1.2.2's patient in train: its one slice disagrees
    assert lid["lab_split_train_disagree_5px"] == lid["lab_split_train_disagree_rule"] == 1
    assert rep["k11"] == {"lidc_patient_leaks": [], "lidc_identical_slices_across_splits": []}
    printed = json.loads(capsys.readouterr().out)
    assert "manifest" not in printed and printed["lidc"] == lid


def test_store_main_refuses_a_partial_lidc_database(lidc_repo, monkeypatch):
    from objdel.cli.data import store as cli
    from objdel.data import lidc

    repo, _ = lidc_repo
    some = dict(list(lidc.scans().items())[:-1])
    monkeypatch.setattr(lidc, "scans", lambda: some)
    with pytest.raises(AssertionError):
        cli.main(["--out", str(repo / "data" / "store")])
    assert not (repo / "data" / "store").exists()


def test_store_slice_flags_of_a_sub_filter_rater():
    from objdel.cli.data import store as cli
    from objdel.data import lidc

    speck = np.zeros((3, 30, 30), bool)
    speck[:2, 5:25, 5:25] = True
    speck[2, 14:16, 14:16] = True  # 4 px: below the size rule, its filtered signature is ()
    s = lidc.Slice(0, 1, [1, 2, 3], [1, 2, 3], False, 200, 200, speck, {}, [])
    assert cli.slice_flags(s) == {
        "disagree_5px": True, "disagree_rule": True, "disagree_rule_no_dup_reader": True, "with_copied_raters": False,
        "copied_layers": 0, "sub_filter_rater": True, "disagree_only_from_sub_filter_rater": True}
    assert not cli.crowded(s) and cli.crowded(lidc.Slice(0, 1, [1, 2], [1, 2], False, 200, 200, speck[:2], {},
                                                         [(222, 200, 240, 230)]))


# -------------------------------------------------------------------------------------------------------- views
def test_views_main_inline_rows_manifest_and_report(tmp_path, monkeypatch, capsys, synthetic_items):
    from collections import Counter

    from objdel.cli.data import views
    from objdel.core.hashing import file_sha256

    store.write(synthetic_items(), tmp_path / "store")
    (tmp_path / "prepared").mkdir()
    (tmp_path / "prepared" / "manifest.json").write_text("{}")
    monkeypatch.setattr(views, "get_context", inline_context)
    out = tmp_path / "views" / "lidc-train"
    views.main(["--store", str(tmp_path / "store"), "--prepared", str(tmp_path / "prepared"), "--out", str(out),
                "--epochs", "2", "--seed", "3", "--procs", "1"])
    m = json.loads((out / "manifest.json").read_text())
    rows = _jsonl(out / "rows.jsonl")
    # z001: a solid disc (one row) and a ring (non-simple: 4 rows); z003: one rater; z002 has no pixels (skipped)
    assert Counter((r["epoch"], r["item_id"][-4:], r["rater"]) for r in rows) == {
        (e, i, k): n for e in (0, 1) for i, k, n in (("z001", "1", 1), ("z001", "2", 4), ("z003", "1", 1))}
    assert m["counts"]["rows"] == len(rows) == 12 and m["counts"]["skipped_no_image"] == 1
    assert m["counts"]["nonsimple_rows"] == 8 and m["counts"]["errors"] == 0
    assert sum(v["rows"] for v in m["counts"]["per_op"].values()) == len(rows)
    assert all(v["exact"] == v["rows"] for v in m["counts"]["per_op"].values())  # LIDC targets are pixel-exact
    assert m["view_hash"] == file_sha256(out / "rows.jsonl")
    assert m["contexts"] == {i: file_sha256(out / "context" / f"{i}.npz") for i in {r["item_id"] for r in rows}}
    assert m["sources"] == {"store_manifest_sha256": file_sha256(tmp_path / "store" / "manifest.json"),
                            "prepared_manifest_sha256": file_sha256(tmp_path / "prepared" / "manifest.json")}
    printed = capsys.readouterr().out
    assert json.JSONDecoder().raw_decode(printed)[0] == m["counts"]
    assert json.loads(views.report_path(out).read_text())["context_files"] == 2


# --------------------------------------------------------------------------------------------------- fetch_lidc
def _dicom_zip(slices, slope=1.0, intercept=-1024.0):
    """Bytes of a zip of CT DICOMs, one per (SOPInstanceUID, z, stored uint16 pixels)."""
    import io
    import zipfile

    from pydicom.dataset import Dataset, FileMetaDataset
    from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for k, (sop, z, px) in enumerate(slices):
            ds = Dataset()
            ds.SOPInstanceUID, ds.SOPClassUID, ds.ImagePositionPatient = sop, CTImageStorage, [0.0, 0.0, z]
            ds.Rows, ds.Columns = px.shape
            ds.BitsAllocated, ds.BitsStored, ds.HighBit, ds.PixelRepresentation = 16, 16, 15, 0
            ds.SamplesPerPixel, ds.PhotometricInterpretation = 1, "MONOCHROME2"
            ds.RescaleSlope, ds.RescaleIntercept, ds.PixelSpacing, ds.SliceThickness = slope, intercept, [0.7, 0.7], 2.5
            ds.PixelData = px.astype("<u2").tobytes()
            ds.file_meta = FileMetaDataset()
            ds.file_meta.TransferSyntaxUID, ds.file_meta.MediaStorageSOPInstanceUID = ExplicitVRLittleEndian, sop
            ds.file_meta.MediaStorageSOPClassUID = CTImageStorage
            one = io.BytesIO()
            ds.save_as(one, enforce_file_format=True)
            zf.writestr(f"dicom/{k}.dcm", one.getvalue())
    return buf.getvalue()


@pytest.fixture
def tcia(tmp_path, monkeypatch):
    """One XML series (1.2.1, contours on z -10 and -7.5) and a fake TCIA: urlopen serves `tcia.body` for the
    series URL (failing `tcia.fail` times first); sleeps are skipped. Returns the fake (calls in .urls)."""
    import io

    from objdel.cli.data import fetch_lidc
    from tests.support.lidc_xml import nodule_xml, square, write

    pytest.importorskip("pydicom")
    nod = nodule_xml([(-10.0, square(10, 10, 20, 20)), (-7.5, square(10, 10, 20, 20))])
    write(tmp_path / "xml", [[nod]], name="1.xml", series="1.2.1")
    pixels = {z: np.random.default_rng(int(-z * 10)).integers(0, 2000, (8, 8)) for z in (-10.0, -7.5, -5.0)}
    # SOP s-10.0 is the one the XML cites (roi_xml); -7.5 is matched by z only; -5.0 has no contour
    body = _dicom_zip([("s-10.0", -10.0, pixels[-10.0]), ("1.9.75", -7.5, pixels[-7.5]),
                       ("1.9.5", -5.0, pixels[-5.0])])

    class Fake:
        urls, fail, px = [], 0, pixels

        def __call__(self, url, timeout):
            self.urls.append(url)
            if self.fail:
                self.fail -= 1
                raise OSError("connection reset")
            return io.BytesIO(body)

    fake = Fake()
    monkeypatch.setattr(fetch_lidc.urllib.request, "urlopen", fake)
    monkeypatch.setattr(fetch_lidc.time, "sleep", lambda s: None)
    return fake


@pytest.mark.filterwarnings("ignore:Invalid value for VR UI")  # the XML's SOP uids ("s-10.0", roi_xml)
def test_fetch_lidc_main_keeps_annotated_slices_in_hu(tmp_path, capsys, tcia):
    from objdel.cli.data import fetch_lidc

    tcia.fail = 2  # two dropped connections, then the download
    out = tmp_path / "slices"
    fetch_lidc.main([str(tmp_path / "xml"), str(out), "--workers", "1"])
    assert tcia.urls == [fetch_lidc.API + "1.2.1"] * 3
    with np.load(out / "1.2.1.npz") as z:
        assert sorted(z.files) == ["z-10.00", "z-7.50"]
        for key, zmm in (("z-10.00", -10.0), ("z-7.50", -7.5)):
            assert z[key].dtype == np.int16 and np.array_equal(z[key], tcia.px[zmm] - 1024)
    assert json.loads((out / "1.2.1.json").read_text()) == {
        "pixel_spacing": [0.7, 0.7], "slice_thickness": 2.5, "slices": 2, "annotated_sops": 2, "annotated_z": 2,
        "matched_by_z": 1}
    assert not (tmp_path / "lidc_tmp" / "1.2.1").exists()  # the DICOM is deleted
    assert "[1/1] 1.2.1 2 slices" in capsys.readouterr().out
    fetch_lidc.main([str(tmp_path / "xml"), str(out)])  # fetched series are skipped
    assert len(tcia.urls) == 3 and capsys.readouterr().out == ""


@pytest.mark.filterwarnings("ignore:Invalid value for VR UI")  # the XML's SOP uids ("s-10.0", roi_xml)
def test_fetch_lidc_refetches_incomplete_and_reports_failures(tmp_path, capsys, tcia):
    from objdel.cli.data import fetch_lidc

    out = tmp_path / "slices"
    out.mkdir()
    (out / "1.2.1.npz").write_bytes(b"")
    (out / "1.2.1.json").write_text(json.dumps({"slices": 1}))  # fewer than the 2 annotated z positions
    assert fetch_lidc.incomplete(out, {"1.2.1": {"a": -10.0, "b": -7.5}, "1.2.9": {"c": 1.0}}) == ["1.2.1"]
    tcia.fail = 4  # every attempt fails: the series is reported, the run goes on
    fetch_lidc.main([str(tmp_path / "xml"), str(out), "--refetch-incomplete"])
    assert len(tcia.urls) == 4 and "[1/1] 1.2.1 FAILED OSError: connection reset" in capsys.readouterr().out
    fetch_lidc.main([str(tmp_path / "xml"), str(out), "--series", "1.2.1"])
    assert json.loads((out / "1.2.1.json").read_text())["slices"] == 2


def test_fetch_lidc_refuses_a_planted_series_uid(tmp_path, capsys):
    from objdel.cli.data import fetch_lidc

    with pytest.raises(SystemExit) as e:
        fetch_lidc.main([str(tmp_path / "xml"), str(tmp_path / "out"), "--series", "../.."])
    assert e.value.code == 2 and "not a DICOM UID" in capsys.readouterr().err
