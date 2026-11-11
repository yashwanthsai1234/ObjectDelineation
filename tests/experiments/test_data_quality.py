"""objdel.data.quality gates DG1-DG18 (B1.3; plan B.9 U18-U24, U55, U69, U74, I03, I04, I18, I30, I31, E04).

The fixture: a tmp store of the 30 synthetic train items of test_data_windows (each with its own HU pixels) plus 6 val
and 6 test items, the val and test views of build_eval_views, and a rows root (pool and schedule seed 4) of one
13-step block (2 rows per step: 13 centred and the 13 draw-0 shifted windows, so the off-centre rows are exactly the
train shifted law)."""
import hashlib
import json
import shutil
import subprocess
import sys
import types
from dataclasses import replace

import numpy as np
import pytest
from test_data_windows import lidc_item, store30

from objdel.data import quality as Q
from objdel.data import store
from objdel.data import windows as Wn
from objdel.data.targets.frame import frame

SHARES = {"centred": 0.5, "shifted": 0.5, "rescaled": 0.0}
PLAN = {"steps": 13, "per_step": 2, "block_steps": 13, "seed": 4, "draws": 1, "shares": SHARES,
        "disagreement_weight": 1.0}


def hu(i, n=200):
    return np.random.default_rng(1000 + i).integers(-1000, 400, (n, n), dtype=np.int16)


def all_items():
    train = [replace(it, image=hu(i)) for i, it in enumerate(store30())]
    other = []
    for s, split in ((2, "val"), (3, "test")):
        for i in range(6):
            it = lidc_item(f"lidc/{s}.2.{i}/n000/z001", [(100, 100, 12), (100, 100, 13)], group=f"LIDC-IDRI-{s}00{i}")
            other.append(replace(it, split=split, image=hu(100 * s + i)))
    return train, other


def write_store(root, items):
    store.write(sorted(items, key=lambda it: (it.split, it.item_id)), root)
    return root


def write_rows(root, items, scheduled, plan=PLAN, excluded=()):
    """A rows root as B1.4 step 3 writes it: one prepare._prepare record per scheduled window, keyed by window_id."""
    from objdel.data import prepare as D
    by_id, recs = {it.item_id: it for it in items}, []
    for w in scheduled:
        it = by_id[w.item_id]
        f = frame(it, offset=(w.dx, w.dy), side=w.side)
        r = D._prepare((str(root), w.window_id, ("nodule",), it.group_id, it.raters, f))
        recs.append(r | {"source_item_id": w.item_id, "cell": f"{w.cell[0]}{w.cell[1]}",
                         "view": {"kind": w.kind, "draw": w.draw, "offset": [w.dx, w.dy], "side": w.side,
                                  "fallback": False, "kept_side": False, "op": "id"}})
    (root / "lidc").mkdir(parents=True, exist_ok=True)
    (root / "lidc" / "train.jsonl").write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in recs))
    ex = json.dumps({i: {"reason": "test"} for i in excluded}).encode()
    sched = plan | {"excluded_ids": sorted(excluded), "exclusions_sha256": hashlib.sha256(ex).hexdigest()}
    (root / "manifest.json").write_text(json.dumps({"schedule": sched}))
    return root


def edit_records(root, fn):
    path = root / "lidc" / "train.jsonl"
    recs = [json.loads(x) for x in path.read_text().splitlines()]
    fn(recs)
    path.write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in recs))


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """(tmp dir, store root, rows root, train items, pool, schedule, eval views root)."""
    from objdel.cli.data import eval_views as B
    tmp = tmp_path_factory.mktemp("q")
    train, other = all_items()
    st = write_store(tmp / "store", train + other)
    P = Wn.pool(train, draws=1, seed=PLAN["seed"])  # B1.4: one seed for the pool and the schedule
    law = [w for w in P.windows if w.kind == "shifted"]
    sched = Wn.schedule(P.windows, Wn.cell_shares(law), 13, 2, SHARES, PLAN["seed"], 1.0, 13)
    rows = write_rows(tmp / "rows", train, sched)
    views = tmp / "eval_views"
    B.build(st, views, split="val", procs=1)
    B.build(st, views, split="test", procs=1)
    B.build(st, views, split="val", procs=1, kinds=("edge",))
    return types.SimpleNamespace(tmp=tmp, store=st, rows=rows, train=train, pool=P, sched=sched, views=views)


@pytest.fixture
def rows(world, tmp_path):
    """A copy of the clean rows root to plant failures in."""
    return shutil.copytree(world.rows, tmp_path / "rows")


# --------------------------------------------------------------------------------------------- import (U69)
def test_quality_imports_no_torch_and_no_transformers():
    code = "import objdel.data.quality, sys; assert 'transformers' not in sys.modules and 'torch' not in sys.modules"
    subprocess.run([sys.executable, "-I", "-c", code], check=True)


# ----------------------------------------------------------------------------------- near-duplicate score (U74)
def seed_unit(crop):  # docs/implementation2.0/seeds/b1_scripts/b1_pixels.py:67
    x = crop.astype(np.float32).reshape(1, -1)
    x = x - x.mean(1, keepdims=True)
    return (x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-6))[0]


def test_near_score_is_a_pearson_correlation_of_the_padded_crop():
    rng = np.random.default_rng(0)
    a = rng.integers(-1000, 400, (32, 32)).astype(np.int16)
    assert float(Q.unit(a) @ Q.unit(a + 50)) == pytest.approx(1.0, abs=1e-5)
    assert float(Q.unit(a) @ Q.unit(-a)) == pytest.approx(-1.0, abs=1e-5)
    img = rng.integers(-1000, 400, (100, 100)).astype(np.int16)
    c = Q.crop32(img, (0, 0, 10, 10))  # centre (5, 5): rows and columns 0..20 of the slice, padded below and right
    assert c.shape == (32, 32) and c.dtype == np.int16
    assert np.array_equal(c[:21, :21], img[:21, :21]) and (c[21:] == -1000).all() and (c[:, 21:] == -1000).all()
    for _ in range(20):
        x, y = rng.integers(-1000, 400, (2, 32, 32)).astype(np.int16)
        assert float(Q.unit(x) @ Q.unit(y)) == pytest.approx(float(seed_unit(x) @ seed_unit(y)), abs=1e-6)


def unit2(score):
    """Two unit vectors in R^1024 whose dot product is `score`."""
    a, b = np.zeros(1024, np.float32), np.zeros(1024, np.float32)
    a[0], b[0], b[1] = 1, score, np.sqrt(1 - score * score)
    return a, b


def test_near_duplicates_list_cross_split_pairs_at_or_above_098():  # U24
    a, b = unit2(0.98)
    c, d = (np.roll(v, 4) for v in unit2(0.9799))  # each pair in its own two axes
    meta = [("t1", "train", "P1"), ("v1", "val", "P2"), ("t2", "train", "P3"), ("s2", "test", "P4"),
            ("t3", "train", "P5"), ("t4", "train", "P6")]
    far = np.zeros(1024, np.float32)
    far[5] = 1
    x, y = unit2(0.99)
    x, y = np.roll(x, 10), np.roll(y, 10)
    near = Q.near_duplicates(meta, np.stack([a, b, c, d, x, y]))
    assert [(p["a"], p["b"]) for p in near["pairs"]] == [("t1", "v1")]
    assert near["pairs"][0]["score"] == pytest.approx(0.98, abs=1e-6)
    assert near["within_train"] == 2  # t3 and t4: information only
    assert Q.near_duplicates(meta[:2], np.stack([a, far]))["pairs"] == []


def test_a_near_pair_is_reviewed_by_an_owner_ok_or_a_drop_listed_in_the_exclusions():
    p = {"a": "t1", "b": "v1", "splits": "train-val"}
    assert not Q.reviewed(p, {}, {})
    assert Q.reviewed(p, {"t1 v1": "owner OK, 2026-11-09"}, {})
    assert not Q.reviewed(p, {"t1 v1": "drop train item"}, {})
    assert Q.reviewed(p, {"t1 v1": "drop train item"}, {"t1": {"reason": "DG3 t1 v1"}})
    assert not Q.reviewed({"a": "v1", "b": "s1", "splits": "val-test"}, {"v1 s1": "drop train item"}, {"v1": {}})


# ----------------------------------------------------------------------------------- leakage (I03, I30)
def test_a_clean_store_has_no_leakage(world):
    got = Q.leakage(world.store, {}, {}, world.tmp / "audit")
    assert {k: got[k] for k in ("patients", "series", "items", "slice_sha", "crop_sha")} == dict.fromkeys(
        ("patients", "series", "items", "slice_sha", "crop_sha"), 0)
    assert got["near"] == []


def test_planted_leakage_fails_dg1_and_dg2(tmp_path):
    train, other = all_items()
    a = other[0]
    moved = replace(a, item_id=a.item_id.replace("/z001", "/z009"), split="train")  # its patient and series in train
    got = Q.leakage(write_store(tmp_path / "s1", train + other + [moved]), {}, {}, tmp_path / "o")
    assert got["patients"] == 1 and got["series"] == 1 and got["items"] == 0
    twin = replace(train[3], image=train[0].image)  # one slice, two patients (items 0-2 are patient 0)
    got = Q.leakage(write_store(tmp_path / "s2", [*train[:3], twin, *train[4:], *other]), {}, {}, tmp_path / "o")
    assert got["slice_sha"] == 1 and got["crop_sha"] == 0
    img, (x0, y0, x1, y1) = train[3].image.copy(), Q.union_box(train[3])
    cx, cy = (x0 + x1) // 2, (y0 + y1) // 2
    img[cy - 16:cy + 16, cx - 16:cx + 16] = Q.crop32(train[0].image, Q.union_box(train[0]))  # one crop, two patients
    got = Q.leakage(write_store(tmp_path / "s3", [*train[:3], replace(train[3], image=img), *train[4:], *other]), {},
                    {}, tmp_path / "o")
    assert got["slice_sha"] == 0 and got["crop_sha"] == 1
    same = replace(train[1], image=train[0].image)  # one patient: not a leak
    got = Q.leakage(write_store(tmp_path / "s4", [train[0], same, *train[2:], *other]), {}, {}, tmp_path / "o")
    assert got["slice_sha"] == 0 and got["crop_sha"] == 0
    with pytest.raises(KeyError, match="lidc/val.parquet"):  # no val file: an error, never a vacuous pass
        Q.leakage(write_store(tmp_path / "s5", train), {}, {}, tmp_path / "o")


def test_a_train_crop_copied_into_val_is_a_near_pair_until_reviewed(tmp_path):  # I30
    train, other = all_items()
    t0 = train[0]
    v = other[0]
    ub = (100 - 13, 100 - 13, 101 + 13, 101 + 13)
    img = v.image.copy()
    c = Q.crop32(t0.image, Q.union_box(t0))
    cx, cy = (ub[0] + ub[2]) // 2, (ub[1] + ub[3]) // 2
    img[cy - 16:cy + 16, cx - 16:cx + 16] = c + 7  # the same crop, 7 HU brighter: score 1, other bytes
    st = write_store(tmp_path / "s", [*train, replace(v, image=img), *other[1:]])
    got = Q.leakage(st, {}, {}, tmp_path / "audit")
    (p,) = got["near"]
    assert (p["a"], p["b"], p["splits"], p["reviewed"]) == (t0.item_id, v.item_id, "train-val", False)
    assert p["score"] >= 0.98 and (tmp_path / "audit" / p["png"]).is_file()
    assert got["crop_sha"] == 0
    key = f"{t0.item_id} {v.item_id}"
    assert not Q.leakage(st, {key: "drop train item"}, {}, tmp_path / "audit")["near"][0]["reviewed"]
    ex = {t0.item_id: {"reason": f"DG3 {key}"}}
    assert Q.leakage(st, {key: "drop train item"}, ex, tmp_path / "audit")["near"][0]["reviewed"]
    assert t0.item_id not in {w.item_id for w in Wn.pool(train, 1, exclude=frozenset(ex)).windows}


# ----------------------------------------------------------------------------------- rows_leakage (U55)
def test_rows_leakage_names_each_planted_row(world, rows):
    views = Q.view_image_hashes(world.views)
    clean = Q.rows_leakage(rows, world.store, views, frozenset())
    assert clean == {"not_train": [], "val_test_patient": [], "view_image": [], "excluded": []}
    val = json.loads((world.views / "lidc-val-shifted/lidc/val.jsonl").read_text().splitlines()[0])
    test = json.loads((world.views / "lidc-test-shifted/lidc/test.jsonl").read_text().splitlines()[0])

    def plant(recs):
        recs[0]["source_item_id"] = val["item_id"]
        recs[1]["group_id"] = test["group_id"]
        recs[2]["image_sha256"] = val["image_sha256"]
        recs[3]["image_sha256"] = test["image_sha256"]
    edit_records(rows, plant)
    recs = Q.records(rows)
    got = Q.rows_leakage(rows, world.store, views, frozenset({recs[4]["source_item_id"]}))
    assert got == {"not_train": [recs[0]["item_id"]], "val_test_patient": [recs[1]["item_id"]],
                   "view_image": [recs[2]["item_id"], recs[3]["item_id"]],
                   "excluded": [r["item_id"] for r in recs if r["source_item_id"] == recs[4]["source_item_id"]]}


# ----------------------------------------------------------------------------------- window_checks (U20, I04, I18)
def test_window_checks_of_the_clean_root(world):
    got = Q.window_checks(world.rows, world.store)
    assert got["windows"] == 26 and got["inside"] == 26 and got["clear"] == 26
    assert (got["fallback"], got["kept_side"], got["padded"], got["mismatch"]) == (0, 0, 0, [])
    assert set(got) >= {"low_contrast"}


def test_a_window_cutting_one_doctor_pixel_fails_dg4_and_dg13(world, rows):
    def cut(recs):
        r = next(r for r in recs if r["view"]["kind"] == "shifted")
        it = next(i for i in world.train if i.item_id == r["source_item_id"])
        x1 = Q.union_box(it)[2]  # the window's right edge at x1 - 1 leaves the union's last column out
        r["box"] = [x1 - 1 - r["view"]["side"], r["box"][1], x1 - 1, r["box"][3]]
    edit_records(rows, cut)
    got = Q.window_checks(rows, world.store)
    assert got["inside"] == 25 and len(got["mismatch"]) == 1
    g = {x.name: x for x in Q.gates(clean_report() | {"windows": got})}
    assert not g["DG4"].ok and g["DG4"].value == pytest.approx(25 / 26) and not g["DG13"].ok


def test_clearance_one_fails_dg5_and_two_passes(world, rows):
    def move(recs, gap):
        r = next(r for r in recs if r["view"]["kind"] == "shifted")
        it = next(i for i in world.train if i.item_id == r["source_item_id"])
        x0 = Q.union_box(it)[0] - gap
        r["box"] = [x0, r["box"][1], x0 + r["view"]["side"], r["box"][3]]
    edit_records(rows, lambda recs: move(recs, 1))
    assert Q.window_checks(rows, world.store)["clear"] == 25
    edit_records(rows, lambda recs: move(recs, 2))
    assert Q.window_checks(rows, world.store)["clear"] == 26


def test_a_centred_policy_frame_with_another_nodule_fails_dg5(world, tmp_path):
    train, other = all_items()
    it = train[0]
    ub = Q.union_box(it)
    crowd = replace(it, meta={"other_nodules": [[ub[2] + 1, ub[1], ub[2] + 3, ub[3]]]})  # beside the union
    st = write_store(tmp_path / "s", [crowd, *train[1:], *other])
    c = Wn.Window(it.item_id, it.group_id, "n", "centred", 0, 0, 0, frame(it).W, frame(it).W, 0.0, 10, False)
    root = write_rows(tmp_path / "rows", [crowd], [c])
    got = Q.window_checks(root, st)
    assert got["windows"] == 1 and got["clear"] == 0


def test_fallback_and_kept_side_windows_fail_dg6(world, rows):
    def flag(recs):
        recs[0]["view"]["fallback"] = True
        recs[1]["view"]["kept_side"] = True
    edit_records(rows, flag)
    got = Q.window_checks(rows, world.store)
    assert (got["fallback"], got["kept_side"]) == (1, 1)
    assert not {x.name: x for x in Q.gates(clean_report() | {"windows": got})}["DG6"].ok
    assert Q.window_checks(world.rows, world.store)["kept_side"] == 0  # centred rows are not kept-side windows


def test_a_padded_window_fails_dg18(world, rows):
    def pad(recs):
        r = recs[0]
        r["box"] = [-5, r["box"][1], -5 + r["view"]["side"], r["box"][3]]
    edit_records(rows, pad)
    got = Q.window_checks(rows, world.store)
    assert got["padded"] == 1 and not {x.name: x for x in Q.gates(clean_report() | {"windows": got})}["DG18"].ok


def test_one_tampered_png_byte_is_exactly_one_dg13_mismatch(world, rows):  # I04
    rec = Q.records(rows)[5]
    png = rows / "images" / f"{rec['item_id']}.png"
    b = bytearray(png.read_bytes())
    b[len(b) // 2] ^= 0xFF
    png.write_bytes(bytes(b))
    got = Q.window_checks(rows, world.store)
    assert got["mismatch"] == [rec["item_id"]]
    edit_records(rows, lambda recs: recs.reverse())  # every window, whatever the order
    assert Q.window_checks(rows, world.store)["mismatch"] == [rec["item_id"]]


# ----------------------------------------------------------------------------------- schedule_report (I18, U22, U59)
def val_windows(world):
    return Q.view_windows(world.views / "lidc-val-shifted", world.store)


def report_of(world, sched=None):
    return Q.schedule_report(world.sched if sched is None else sched, world.pool,
                             {"val": val_windows(world)}, **{k: PLAN[k] for k in Q.SCHEDULE_ARGS})


def test_schedule_report_keys_and_coverage_against_the_train_shifted_law(world):
    got = report_of(world)
    law = [w for w in world.pool.windows if w.kind == "shifted" and w.draw == 0]
    want = Q.cell_coverage(world.sched, Wn.cell_shares(law), [w.offset for w in law])
    assert got["prefixes"] == [want] and want["tv"] == 0.0 and want["ks"] == 0.0 and want["worst_ratio"] == 1.0
    assert got["mix"]["rows"] == 26
    assert got["s0_share"] == 1.0 and got["kinds"] == {"centred": 13, "shifted": 13, "rescaled": 0}
    assert got["count_mismatch"] == [] and got["unknown"] == 0
    assert got["dropped"] == dict(world.pool.dropped)
    for k in ("disagreement_share", "max_disagreement_reuse", "top10_nodule_share", "cells", "rescaled"):
        assert k in got
    v = got["val"]["val"]
    assert set(v) == {"patients", "coverage", "null"} and v["patients"] == 6
    assert set(v["null"]) == {"worst_ratio", "tv", "ks"}


def test_dg11_fails_when_one_row_moves_to_another_cell(world):
    s = list(world.sched)
    i = next(i for i, w in enumerate(s) if w.kind == "shifted")
    a, b = s[i].cell
    s[i] = replace(s[i], offset=(0.05, 0.17, 0.25, 0.35)[(b + 1) % 4])  # same kind, next offset bin
    got = report_of(world, s)
    assert len(got["count_mismatch"]) == 2 and got["kinds"] == report_of(world)["kinds"]  # one cell short, one over
    assert f"cell {(a, b)} has" in got["count_mismatch"][0] + got["count_mismatch"][1]


def test_a_centred_only_schedule_fails_dg7_and_dg8(world):
    s = [w for w in world.sched if w.kind == "centred"]
    g = {x.name: x for x in Q.gates(clean_report() | {"schedule": report_of(world, s)})}
    assert not g["DG7"].ok and not g["DG8"].ok and not g["DG9"].ok


def test_a_first_block_off_target_fails_dg7_on_its_prefix(world):
    good = report_of(world)
    bad = dict(good, prefixes=[{"rows": 13, "worst_ratio": 0.5, "tv": 0.3, "ks": 0.2}, *good["prefixes"]])
    g = {x.name: x for x in Q.gates(clean_report() | {"schedule": bad})}
    assert not g["DG7"].ok and g["DG7"].value == 0.5


def test_schedule_report_reads_the_apparent_size(world, monkeypatch):  # U59
    before = report_of(world)
    monkeypatch.setattr(Wn.Window, "apparent_px", property(lambda w: 5.0))
    assert report_of(world)["val"]["val"]["coverage"] != before["val"]["val"]["coverage"]


# ----------------------------------------------------------------------------------- tokens (U23), agreement (I31)
class Tok:
    def __init__(self, n):
        self.n = n

    def __call__(self, texts, add_special_tokens=False):
        assert not add_special_tokens
        return {"input_ids": [[0] * (self.n if i == 0 else 3) for i, _ in enumerate(texts)]}


@pytest.mark.parametrize("n,ok", [(1536 - 102, True), (1536 - 101, False), (1536 - 100, False)])
def test_dg12_fencepost(world, n, ok):
    got = Q.tokens(world.rows, Tok(n), 1536)
    assert got["max"] == n and got["errors"] == 0
    assert {x.name: x for x in Q.gates(clean_report() | {"tokens": got})}["DG12"].ok is ok


def test_dg12_fails_on_a_target_error(world, rows):
    edit_records(rows, lambda recs: recs[0]["errors"].append(["x", "1", "bad"]))
    got = Q.tokens(rows, Tok(10), 1536)
    assert got["errors"] == 1 and not {x.name: x for x in Q.gates(clean_report() | {"tokens": got})}["DG12"].ok


def stub_ceiling(monkeypatch, values):
    """objdel.verify.ceiling with a loo_caches whose left-out doctor r scores values[r % len(values)] (score_mask
    patched to return the stub cache's value)."""
    from objdel.verify.reward import score

    def loo_caches(cache):
        return [types.SimpleNamespace(index=r, cache=values[r % len(values)]) for r in range(cache.k)]
    monkeypatch.setitem(sys.modules, "objdel.verify.ceiling", types.SimpleNamespace(loo_caches=loo_caches))
    monkeypatch.setattr(score, "score_mask", lambda mask, cache: types.SimpleNamespace(iou_tc=cache))


def test_agreement_reports_the_doctor_loo_iou_per_cell(tmp_path, monkeypatch):
    it = lidc_item("lidc/1.2.0/n000/z001", [(100, 100, 12), (100, 100, 13), (100, 100, 14)])
    S = frame(it).W
    w = [Wn.Window(it.item_id, it.group_id, "n", "centred", 0, 0, 0, S, S, 0.0, 29, False),
         Wn.Window(it.item_id, it.group_id, "n", "shifted", 0, 5, 0, S, S, 0.1, 29, False)]
    root = write_rows(tmp_path / "rows", [it], w)
    stub_ceiling(monkeypatch, [0.9, 0.8, 0.1])  # three doctors, one outlier: mean 0.6
    got = Q.agreement(root, low_contrast={w[0].window_id})
    assert got["windows"] == 2 and got["mean"] == pytest.approx(0.6)
    assert got["share_lt_0.3"] == 0.0 and got["share_lt_0.3_low_contrast"] == 0.0
    stub_ceiling(monkeypatch, [0.2, 0.3, 0.25])
    got = Q.agreement(root, low_contrast={w[0].window_id})
    assert got["share_lt_0.3"] == 1.0 and got["share_lt_0.3_low_contrast"] == 0.5
    assert got["cells"] == {"20": {"n": 2, "mean": pytest.approx(0.25), "share_lt_0.3": 1.0}}


# ----------------------------------------------------------------------------------- frozen views (DG15)
def test_frozen_views_need_equal_content_and_target_set_hashes(world, tmp_path):
    entries = {v: {"content_sha256": Q.view_content_sha256(world.views / v),
                   "target_set_sha256": Q.target_set_sha256(world.views / v)} for v in Q.FROZEN}
    assert Q.frozen(world.views, entries) == {}
    assert set(Q.frozen(world.views, {})) == set(Q.FROZEN)
    views = shutil.copytree(world.views, tmp_path / "v")
    m = views / "lidc-val-shifted/manifest.json"
    m.write_text(json.dumps(json.loads(m.read_text()) | {"versions": {}, "reward_fingerprint": "x"}))
    assert Q.frozen(views, entries) == {}
    path = views / "lidc-val-shifted/lidc/val.jsonl"
    recs = [json.loads(x) for x in path.read_text().splitlines()]
    for edit in (lambda r: r[0]["view"]["offset"].__setitem__(0, r[0]["view"]["offset"][0] + 1),  # I05: +1 px
                 lambda r: r.__setitem__(slice(0, 2), r[1::-1]),  # two records reordered
                 lambda r: r[0].__setitem__("image_sha256", "0" * 64)):  # one image changed
        bad = [json.loads(json.dumps(r)) for r in recs]
        edit(bad)
        path.write_text("".join(json.dumps(r) + "\n" for r in bad))
        assert Q.frozen(views, entries) == {"lidc-val-shifted": "content_sha256 differs"}


# ----------------------------------------------------------------------------------- gates (U19, U21)
def clean_report():
    cov = {"rows": 100, "worst_ratio": 1.0, "tv": 0.0, "ks": 0.0}
    return {
        "leakage": {"patients": 0, "series": 0, "items": 0, "slice_sha": 0, "crop_sha": 0, "near": [],
                    "within_train": 0},
        "rows": {"not_train": [], "val_test_patient": [], "view_image": [], "excluded": [], "exclusions_stale": False},
        "windows": {"windows": 10, "inside": 10, "clear": 10, "fallback": 0, "kept_side": 0, "padded": 0,
                    "mismatch": [], "low_contrast": []},
        "schedule": {"prefixes": [cov], "mix": cov, "s0_share": 1.0, "count_mismatch": [], "unknown": 0,
                     "disagreement_share": 0.1, "val": {"val": {"coverage": cov}}},
        "tokens": {"max": 100, "limit": 1536, "errors": 0},
        "agreement": {"share_lt_0.3": 0.0},
        "frozen": {},
        "z": [],
    }


ORDER = ["DG1", "DG1b", "DG2", "DG3", "DG4", "DG5", "DG6", "DG7", "DG8", "DG9", "DG9b", "DG10", "DG11", "DG12",
         "DG13", "DG14", "DG14b", "DG15", "DG16", "DG17", "DG18"]


def test_gates_are_the_21_in_order_and_only_hard_failures_fail():
    g = Q.gates(clean_report())
    assert [x.name for x in g] == ORDER
    assert {x.name for x in g if not x.hard} == {"DG9b", "DG14", "DG14b", "DG16"}
    assert all(x.ok for x in g) and Q.failures(g) == []
    assert all(isinstance(x.limit, str) and isinstance(x.ok, bool) for x in g)
    bad = clean_report()
    bad["frozen"] = {"lidc-test-shifted": "missing"}
    bad["schedule"]["val"]["val"]["coverage"] = {"rows": 1, "worst_ratio": 0.0, "tv": 1.0, "ks": 1.0}
    assert [x.name for x in Q.failures(Q.gates(bad))] == ["DG15"]
    del bad["tokens"]
    with pytest.raises(KeyError):
        Q.gates(bad)


@pytest.mark.parametrize("name,key,ok_value,bad_value", [
    ("DG7", "worst_ratio", 0.8, 0.7999), ("DG8", "tv", 0.10, 0.1001), ("DG9", "ks", 0.06, 0.0601)])
def test_coverage_limits_are_exact(name, key, ok_value, bad_value):
    for v, ok in ((ok_value, True), (bad_value, False)):
        r = clean_report()
        r["schedule"]["prefixes"] = [{"rows": 100, "worst_ratio": 1.0, "tv": 0.0, "ks": 0.0, key: v}]
        assert {x.name: x for x in Q.gates(r)}[name].ok is ok


def test_policy_side_limit_is_exact(world):
    s = [w for w in world.sched if w.kind != "centred"]
    at_s0 = [replace(w, kind="shifted") for w in s[:1]] * 75
    off = [replace(w, kind="rescaled", side=w.s0 * 2) for w in s[:1]]
    for n, ok in ((75, True), (74, False)):
        rows = at_s0[:n] + off * (100 - n)
        r = clean_report()
        r["schedule"]["s0_share"] = report_of(world, rows)["s0_share"]
        assert {x.name: x for x in Q.gates(r)}["DG10"].ok is ok


# ----------------------------------------------------------------------------------- audit_data CLI (E04)
@pytest.fixture
def cli(world, tmp_path, monkeypatch):
    from objdel.cli.data import audit as A
    from objdel.data import splits
    from objdel.policy import models
    stub_ceiling(monkeypatch, [0.9, 0.8])
    monkeypatch.setattr(models, "tokenizer", lambda model: (Tok(50), 1536))
    monkeypatch.setattr(splits, "val_half", lambda g: "shape" if int(g[-1]) % 2 else "select")  # fixture patients
    ex = tmp_path / "train_exclusions.json"
    ex.write_text("{}")
    frozen = tmp_path / "frozen_views.json"
    frozen.write_text(json.dumps({v: {"content_sha256": Q.view_content_sha256(world.views / v),
                                      "target_set_sha256": Q.target_set_sha256(world.views / v)} for v in Q.FROZEN}))
    z = {"train_sha256": Q.file_sha256(world.store / "lidc/train.parquet"), "series": {}, "step2": []}
    (tmp_path / "z.json").write_text(json.dumps(z))

    def run(rows, *extra):
        out = tmp_path / "audit.json"
        code = A.main(["--rows", str(rows), "--store", str(world.store), "--eval-views", str(world.views),
                       "--exclusions", str(ex), "--review", str(tmp_path / "none.json"), "--frozen", str(frozen),
                       "--series-z", str(tmp_path / "z.json"), "--out", str(tmp_path / "png"), "--json", str(out),
                       *extra])
        return code, json.loads(out.read_text())
    return types.SimpleNamespace(run=run, ex=ex, z=tmp_path / "z.json")


def test_audit_data_passes_a_clean_root(cli, world, capsys):
    code, doc = cli.run(world.rows)
    assert code == 0 and [g["name"] for g in doc["gates"]] == ORDER
    out = capsys.readouterr().out
    assert all(name in out for name in ORDER)


def test_audit_data_fails_on_a_hard_gate_and_names_it(cli, world, rows):
    edit_records(rows, lambda recs: recs[0]["view"].update(fallback=True))
    code, doc = cli.run(rows)
    assert code == 1 and [g["name"] for g in doc["gates"] if not g["ok"]] == ["DG6"]


def test_audit_data_dg1b_and_dg17_failures(cli, world):
    rec = Q.records(world.rows)[0]
    cli.ex.write_text(json.dumps({rec["source_item_id"]: {"reason": "x"}}))
    code, doc = cli.run(world.rows)
    assert code == 1 and "DG1b" in [g["name"] for g in doc["gates"] if not g["ok"]]
    cli.ex.write_text("{}")
    z = json.loads(cli.z.read_text())
    s = rec["source_item_id"].split("/")[1]
    z.update(step2=[s], series={s: {"split": "train", "items": 1, "item_ids": [rec["source_item_id"]]}})
    cli.z.write_text(json.dumps(z))
    code, doc = cli.run(world.rows)
    assert code == 1 and [g["name"] for g in doc["gates"] if not g["ok"]] == ["DG17"]


def test_report_gates_never_change_the_exit_code(cli, world, monkeypatch):
    stub_ceiling(monkeypatch, [0.1, 0.1])  # every doctor disagrees: DG14b is a report
    code, doc = cli.run(world.rows)
    assert code == 0
