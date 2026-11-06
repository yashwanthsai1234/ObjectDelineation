"""Phase 2 M1: frames with offset / side / window, dihedral ops, training views (build_views), eval views
(build_eval_views), and the M9 fetch_lidc_pixels fixes (I45). CI tests use the synthetic items of conftest; tests
marked `data` read data/store and data/prepared (skipped when absent)."""
import hashlib
import json
from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from objdel.checker import signature
from objdel.data import store
from objdel.geometry import SCHEMAS, parse, to_mask
from objdel.targets import FRAMES, LUNG_WINDOW, FrameConfig, augment, build_target, frame, min_areas, target_set_hash
from objdel.targets.frame import CLEARANCE, DIVISORS, union_box, upscale

ROOT = Path(__file__).resolve().parents[1]
PREPARED, STORE = ROOT / "data" / "prepared", ROOT / "data" / "store"
BUNDLE, EVAL_VIEWS = ROOT / "data" / "views" / "lidc-bundle-v3", ROOT / "data" / "eval_views"  # the bundle trained on
FRAME_SHA = {  # targets.frame before the Phase 2 change (box, masks, image bytes)
    "lidc/1.2.1/n000/z001": "d70579514adab44459d219e6a20e28496a5d970cdb95606050afc020ee9f8d92",
    "lidc/1.2.1/n000/z002": "b6df0facfcecd3c313c4efa31de3cc7b6b16f9dde59a9c8596d586fa4fd80aba",
    "lidc/1.2.1/n000/z003": "2e557b0bf981b398038634f533c3b83b2d8718f490cc2833ea6db647a0da464d",
    "lidc/1.2.1/n000/z004": "263aebcd28d605109585774ceac640b8895cb6c996b7f6a4be51a8c17ae68cd8",
    "lidc/1.2.2/n000/z001": "abdc0ccc24865fb0617f7a29e6501215eb82618e7870ca9b279fecc15e833e21",
    "lidc/1.2.2/n000/z003": "7975b9867df7753a82f53b66012a01407043a31de1e7fc1a4293a5bd0bfc42a6",
    "lidc/1.2.2/n000/z004": "f1f85d2b162415cc68d5ca75efef6055c06c479e99bd1b845a710324af61d262",
    "riga/MESSIDOR/image001": "5d394aa8a2bfd5641c23af8cc75bf968a836a1ae91f42450e84b9a0ebd2e5205",
    "riga/MESSIDOR/image002": "1fbfd9551c7538f023b148f9117760e581f8439b255c9f7ec3abfca0c981e6cf",
}


def _sha(f):
    return hashlib.sha256(repr(f.box).encode() + b"".join(m.tobytes() for m in f.masks.values())
                          + (f.image.tobytes() if f.image is not None else b"")).hexdigest()


def _pixels_sha(image):  # objdel.train.dataset._pixels_sha (that module needs the train extra)
    return hashlib.sha256(repr(image.shape).encode() + image.tobytes()).hexdigest()


@pytest.fixture(scope="module")
def items(synthetic_items):
    return [it for it in synthetic_items(("train", "val"), unanimous=True) if it.image is not None]


def _reference(it, offset, side, window, op, sigma=0.0, seed=0):
    """Offline reference of render: HU noise on the whole image, then frame at frame resolution, op, upscale."""
    img = it.image
    if sigma:
        S, f = side, frame(replace(it, image=None), offset=offset, side=side)
        x0, y0 = f.box[:2]
        noise = np.rint(augment.row_rng("noise", seed).normal(0, sigma, (S, S))).astype(np.int64)
        img = img.astype(np.int64)
        H, W = img.shape
        ys, xs = slice(max(y0, 0), min(y0 + S, H)), slice(max(x0, 0), min(x0 + S, W))
        img[ys, xs] += noise[ys.start - y0:ys.stop - y0, xs.start - x0:xs.stop - x0]
    cfg = FrameConfig(**{**FRAMES[it.dataset].__dict__, "image_side": None})
    f = frame(replace(it, image=img), cfg, offset=offset, side=side, window=window)
    out = np.moveaxis(augment.apply(np.moveaxis(f.image, -1, 0), op), 0, -1)
    from objdel.targets.frame import upscale
    return upscale(np.ascontiguousarray(out), FRAMES[it.dataset].image_side)


# ------------------------------------------------------------------------------------------------- frame(offset, side)
def test_frame_defaults_are_byte_identical_to_phase1(items):
    for it in items:
        S = frame(it).W
        for f in (frame(it), frame(it, offset=(0, 0), side=None, window=None), frame(it, side=S),
                  frame(it, window=LUNG_WINDOW)):
            assert _sha(f) == FRAME_SHA[it.item_id]


def test_offset_moves_the_window_and_the_crop(items):
    it = next(i for i in items if i.item_id.endswith("z001"))
    base, f = frame(it), frame(it, offset=(5, -7))
    assert f.box == (base.box[0] + 5, base.box[1] - 7, base.box[2] + 5, base.box[3] - 7)
    x0, y0, x1, y1 = f.box
    assert np.array_equal(f.masks["nodule"], it.masks["nodule"][:, y0:y1, x0:x1])
    from objdel.targets.frame import hu_window, upscale
    crop = hu_window(it.image[y0:y1, x0:x1]).astype(np.uint8)
    assert np.array_equal(f.image, upscale(np.repeat(crop[..., None], 3, 2), 512))


def test_side_and_offset_refuse_inadmissible_windows(items):
    it = next(i for i in items if i.item_id.endswith("z001"))  # union 25 px wide, at 38..63
    x0, y0, x1, y1 = union_box(it)
    S = frame(it).W
    with pytest.raises(ValueError, match="divisor"):
        frame(it, side=60)
    with pytest.raises(ValueError, match="divisor"):
        frame(it, side=1000)
    with pytest.raises(ValueError, match="clearance"):
        frame(it, FrameConfig(min_side=20, image_side=512), side=25)  # a divisor, too small for union + 2 x 2 px
    cx = frame(it).box[0]
    left = x0 - CLEARANCE - cx  # the largest dx that keeps 2 px on the left
    frame(it, offset=(left, 0))
    with pytest.raises(ValueError, match="clearance"):
        frame(it, offset=(left + 1, 0))
    other = replace(it, meta={"other_nodules": [[x1 + 4, y0, x1 + 8, y1]]})
    frame(other)  # the policy frame is never refused (the store dropped crowded items)
    edge = x1 + 4 - S - cx  # the window's right edge at the other nodule's left edge
    assert frame(other, offset=(edge, 0)).box[2] == x1 + 4
    with pytest.raises(ValueError, match="other outlined nodule"):
        frame(other, offset=(edge + 1, 0))


def test_window_changes_lidc_pixels_only(items):
    lidc = next(i for i in items if i.dataset == "lidc")
    riga = next(i for i in items if i.dataset == "riga")
    assert not np.array_equal(frame(lidc, window=(-1100, 300)).image, frame(lidc).image)
    assert np.array_equal(frame(riga, window=(-1100, 300)).image, frame(riga).image)


# ------------------------------------------------------------------------------------------------- augment
def test_ops_are_the_dihedral_group():
    a = np.arange(2 * 3 * 4).reshape(2, 3, 4)
    out = {op: augment.apply(a, op) for op in augment.OPS}
    assert len({o.tobytes() + bytes(o.shape) for o in out.values()}) == 8
    for op, o in out.items():
        assert np.array_equal(augment.apply(o, augment.INVERSE[op]), a) and o.flags.c_contiguous
    assert np.array_equal(out["fh"], a[..., ::-1]) and np.array_equal(out["tr"], a.swapaxes(1, 2))
    assert augment.DATASET_OPS["riga"] == ("id", "fh")


@pytest.mark.parametrize("seed", range(40))
def test_draw_invariants(items, seed):
    for it in items:
        rng = np.random.default_rng(seed)
        got = augment.draw(it, rng)
        assert got is not None
        (dx, dy), S = got
        assert 1000 % S == 0 and S <= 500 and S in DIVISORS
        f = frame(it, offset=(dx, dy), side=S)  # raises when the union is not inside with 2 px or crowded
        x0, y0, x1, y1 = union_box(it)
        assert f.scale == 1 and x0 - f.box[0] >= 2 and f.box[2] - x1 >= 2 and y0 - f.box[1] >= 2 and f.box[3] - y1 >= 2


def test_scale_draw_covers_every_eval_side(items):
    """D10: S is uniform over the policy side and its LADDER neighbours, so the draw gives every side the shifted
    (policy S) and rescaled (augment.neighbour_side) eval views can ask for, each about equally often."""
    for it in items:
        S0, n = frame(it).W, 300
        drawn = Counter(augment.draw(it, np.random.default_rng(s))[1] for s in range(n))
        rescaled = {got[1] for s in range(100) if (got := augment.neighbour_side(it, S0, np.random.default_rng(s)))}
        assert {S0} | rescaled <= set(drawn) <= set(augment.sides(it))
        assert min(drawn.values()) > n / len(drawn) / 2, (it.item_id, drawn)
    it = next(i for i in items if i.item_id.endswith("z001"))  # 25 px union, policy S 50: the margin draw gave 50 / 100
    assert augment.sides(it) == [50, 100]
    rg = next(i for i in items if i.dataset == "riga")  # 65 px union, policy S 100: 50 cannot hold 65 + 2 x 2 px
    assert augment.sides(rg) == [50, 100, 125]
    assert {augment.draw(rg, np.random.default_rng(s))[1] for s in range(200)} == {100, 125}


def test_draw_never_shows_another_nodule(items):
    it = next(i for i in items if i.item_id.endswith("z001"))
    x0, y0, x1, y1 = union_box(it)
    other = replace(it, meta={"other_nodules": [[x1 + 4, 0, x1 + 9, 100]]})
    for s in range(100):
        (dx, dy), S = augment.draw(other, np.random.default_rng(s))
        assert frame(other, offset=(dx, dy), side=S).box[2] <= x1 + 4


def test_render_equals_the_offline_reference(items):
    for it in items:
        for s in range(6):
            rng = np.random.default_rng(s)
            (dx, dy), S = augment.draw(it, rng)
            op = augment.DATASET_OPS[it.dataset][s % len(augment.DATASET_OPS[it.dataset])]
            j = augment.jitter(rng) if it.dataset == "lidc" else {"window": None, "noise_sigma": 0, "noise_seed": 0}
            f = frame(replace(it, image=None), offset=(dx, dy), side=S)
            spec = {"dataset": it.dataset, "side": S, "box": list(f.box), "op": op, **j}
            ctx = _context(it)
            ref = _reference(it, (dx, dy), S, j["window"] and tuple(j["window"]), op, j["noise_sigma"],
                             j["noise_seed"])
            assert np.array_equal(augment.render(spec, ctx), ref), (it.item_id, s)


def _context(it, tmp=None):
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        augment.write_context(Path(d, "context", it.item_id + ".npz"), it, frame(it).box)
        return augment.load_context(d, it.item_id)


def test_render_at_zero_jitter_is_the_policy_frame(items):
    for it in items:
        f = frame(it)
        spec = {"dataset": it.dataset, "side": f.W, "box": list(f.box), "op": "id", "window": None,
                "noise_sigma": 0, "noise_seed": 0}
        assert np.array_equal(augment.render(spec, _context(it)), f.image)


def test_context_holds_the_store_pixels(items):
    for it in items:
        c = _context(it)
        x0, y0 = int(c["x0"]), int(c["y0"])
        h, w = c["pixels"].shape[:2]
        assert c["pixels"].dtype == it.image.dtype and np.array_equal(c["pixels"], it.image[y0:y0 + h, x0:x0 + w])
        assert tuple(c["shape"]) == it.image.shape[:2]


def test_render_refuses_a_window_outside_the_context(items):
    it = next(i for i in items if i.dataset == "lidc")
    c, f = _context(it), frame(it)
    cut = f.box[1] + 1 - int(c["y0"])  # the context now starts one row below the window
    c = {**c, "pixels": c["pixels"][cut:], "y0": c["y0"] + cut}
    with pytest.raises(ValueError, match="context"):
        augment.render({"dataset": "lidc", "side": f.W, "box": list(f.box), "op": "id"}, c)


def _certified(text, src, names, ma):
    p = parse(text, names)
    if not p.ok:
        return False, False
    ok, exact = True, True
    for n in names:
        g = to_mask(p.objs[n], src[n].shape[1], src[n].shape[0])
        u = np.count_nonzero(g.mask | src[n])
        ok &= g.legal and np.count_nonzero(g.mask & src[n]) / max(u, 1) >= 0.99 and \
            signature(g.mask, ma[n]) == signature(src[n], ma[n])
        exact &= bool(np.array_equal(g.mask, src[n]))
    return ok, exact


def test_lidc_variants_are_pixel_exact_and_certified(items):
    for it in (i for i in items if i.dataset == "lidc"):
        (dx, dy), S = augment.draw(it, np.random.default_rng(1))
        f = frame(it, offset=(dx, dy), side=S)
        ma = min_areas(f.masks)
        for r in range(len(it.raters)):
            ras = {}
            for op in augment.OPS:
                src = {"nodule": augment.apply(f.masks["nodule"][r], op)}
                t = build_target(src, ("nodule",), None, ma)
                assert _certified(t.text, src, ("nodule",), ma) == (True, True)
                g = to_mask(parse(t.text, ("nodule",)).objs["nodule"], S, S).mask
                ras[op] = augment.apply(g, augment.INVERSE[op])
            assert all(np.array_equal(v, ras["id"]) for v in ras.values())


# ------------------------------------------------------------------------------------------------- build_views
@pytest.fixture(scope="module")
def views(tmp_path_factory, synthetic_items):
    from experiments.phase2 import build_views as BV
    root = tmp_path_factory.mktemp("views")
    store.write(synthetic_items(), root / "store")
    out = [BV.build(root / "store", root / f"v{k}" / "lidc-train", "lidc", "train", epochs=2, seed=3,
                    nonsimple_repeat=4, procs=k) for k in (1, 2)]
    return root, out


def _rows(d):
    return [json.loads(line) for line in (d / "rows.jsonl").read_text().splitlines()]


def test_view_rows_schema_order_and_repeats(views, synthetic_items):
    root, (m, _) = views
    d = root / "v1" / "lidc-train"
    rows = _rows(d)
    keys = {"row_id", "epoch", "item_id", "dataset", "rater", "rater_index", "op", "offset", "side", "box", "window",
            "noise_sigma", "noise_seed", "simple", "rep", "target"}
    assert all(set(r) == keys for r in rows)
    assert rows == sorted(rows, key=lambda r: (r["epoch"], r["item_id"], r["rater_index"], r["rep"]))
    # z001: a solid disc (simple) and a ring (non-simple, x4); z003: one rater; z002 has no pixels (skipped)
    per = Counter((r["epoch"], r["item_id"][-4:], r["rater"]) for r in rows)
    assert per == {(e, i, k): n for e in (0, 1) for i, k, n in
                   (("z001", "1", 1), ("z001", "2", 4), ("z003", "1", 1))}
    ring = [r for r in rows if r["item_id"].endswith("z001") and r["rater"] == "2"]
    assert not any(r["simple"] for r in ring) and all(r["simple"] for r in rows if r not in ring)
    assert len({(tuple(r["offset"]), r["side"], r["op"], r["noise_seed"]) for r in ring}) == len(ring)
    assert m["view_hash"] == hashlib.sha256((d / "rows.jsonl").read_bytes()).hexdigest()
    assert m["counts"]["rows"] == len(rows) and m["view"] == "lidc-train"
    for r in rows:
        lo, hi = r["window"]
        assert -360 <= (lo + hi) / 2 <= -240 and 1320 <= hi - lo <= 1480 and 0 <= r["noise_sigma"] <= 10
        assert r["op"] in augment.OPS and r["side"] in (50, 100)


def test_view_targets_are_certified_per_op(views, synthetic_items):
    root, (m, _) = views
    its = {i.item_id: i for i in synthetic_items()}
    for r in _rows(root / "v1" / "lidc-train"):
        it = its[r["item_id"]]
        f = frame(it, offset=tuple(r["offset"]), side=r["side"])
        assert list(f.box) == r["box"]
        src = {"nodule": augment.apply(f.masks["nodule"][r["rater_index"]], r["op"])}
        assert _certified(r["target"], src, ("nodule",), min_areas(frame(it).masks)) == (True, True)
    assert all(v["exact"] == v["rows"] for v in m["counts"]["per_op"].values())  # build_target: exact for LIDC


def test_two_view_builds_are_byte_identical(views):
    root, (m1, m2) = views
    a, b = root / "v1" / "lidc-train", root / "v2" / "lidc-train"
    fa = sorted(p.relative_to(a) for p in a.rglob("*") if p.is_file())
    assert fa == sorted(p.relative_to(b) for p in b.rglob("*") if p.is_file())
    assert len(fa) == 4  # rows, manifest, 2 contexts: timing lives outside the view (BK1-7)
    assert all((a / p).read_bytes() == (b / p).read_bytes() for p in fa)
    assert (root / "v1" / "build_reports" / "lidc-train.json").is_file()


def test_view_manifest_pins_every_context(views):
    """BK1-3 (write side): manifest contexts = {item_id: SHA-256 of context/<item_id>.npz}, one per item with rows."""
    root, (m, _) = views
    d = root / "v1" / "lidc-train"
    assert set(m["contexts"]) == {r["item_id"] for r in _rows(d)}
    for item_id, sha in m["contexts"].items():
        assert hashlib.sha256((d / "context" / f"{item_id}.npz").read_bytes()).hexdigest() == sha
    assert json.loads((d / "manifest.json").read_text())["contexts"] == m["contexts"]


def test_build_views_defaults_to_three_epochs():
    """S1-8: R1's 2-epoch branch decays into epoch 3, so a default build must hold 3 epochs."""
    import inspect

    from experiments.phase2 import build_views as BV
    assert inspect.signature(BV.build).parameters["epochs"].default == 3


def test_view_rows_render_like_the_reference(views, synthetic_items):
    root, _ = views
    d = root / "v1" / "lidc-train"
    its = {i.item_id: i for i in synthetic_items()}
    for r in _rows(d)[::3]:
        ref = _reference(its[r["item_id"]], tuple(r["offset"]), r["side"], tuple(r["window"]), r["op"],
                         r["noise_sigma"], r["noise_seed"])
        assert np.array_equal(augment.render(r, augment.load_context(d, r["item_id"])), ref)


def test_build_views_help():
    from experiments.phase2 import build_views as BV
    with pytest.raises(SystemExit):
        BV.main(["--help"])


# ------------------------------------------------------------------------------------------------- build_eval_views
@pytest.fixture(scope="module")
def eval_views(tmp_path_factory, synthetic_items):
    pytest.importorskip("datasets")
    from experiments.phase2 import build_eval_views as BE
    root = tmp_path_factory.mktemp("eval_views")
    store.write(synthetic_items(("train", "val")), root / "store")
    return root, BE.build(root / "store", root / "out", ("lidc", "riga"), "val", seed=0, procs=1)


def test_eval_views_read_through_the_prepared_layout(eval_views):
    from objdel.train import dataset as D
    root, manifests = eval_views
    assert not (root / "out" / "riga-val-rotated").exists()  # fundus orientation is physiological (D10)
    for name in ("lidc-val-shifted", "lidc-val-rescaled", "lidc-val-rotated", "riga-val-shifted", "riga-val-rescaled"):
        d = root / "out" / name
        ds = name.split("-")[0]
        m = json.loads((d / "manifest.json").read_text())
        assert m["sets"][f"{ds}/val"]["target_errors"] == [] and m["view"] == name
        recs = D.records(d, ds, "val")
        assert recs and all("view" in r for r in recs)
        load = D.cache_loader(d)
        for r in recs:
            S = r["view"]["side"]
            assert r["W"] == S and r["scale"] == 1
            if r["rl"]:
                assert set(load(r["item_id"])) == set(SCHEMAS[ds])
        if ds == "lidc":
            assert len(D.grpo_rows(d, ds, "val")) == sum(r["rl"] for r in recs)  # images match image_sha256
            D.sft_rows(d, ds, "val")


def test_eval_views_shift_and_rescale(eval_views, synthetic_items):
    from objdel.train import dataset as D
    root, _ = eval_views
    its = {i.item_id: i for i in synthetic_items(("train", "val"))}
    for ds in ("lidc", "riga"):
        sh, rs = (D.records(root / "out" / f"{ds}-val-{k}", ds, "val") for k in ("shifted", "rescaled"))
        for a, b in zip(sh, rs, strict=True):
            it = its[a["item_id"]]
            S0 = frame(it).W
            assert a["view"]["side"] == S0 and not a["view"]["fallback"]
            assert b["view"]["side"] != S0 or b["view"]["kept_side"]
            f = frame(it, offset=tuple(b["view"]["offset"]), side=b["view"]["side"])
            assert list(f.box) == b["box"] and b["image_sha256"] == _pixels_sha(f.image)
            assert all(b["targets"]) and not b["view"]["fallback"]
        assert any(r["view"]["offset"] != [0, 0] for r in sh)


def test_two_eval_builds_are_identical(eval_views, tmp_path):
    from experiments.phase2 import build_eval_views as BE
    root, _ = eval_views
    BE.build(root / "store", tmp_path, ("lidc",), "val", seed=0, procs=2)
    for name in ("lidc-val-shifted", "lidc-val-rescaled", "lidc-val-rotated"):
        a, b = root / "out" / name, tmp_path / name
        assert (a / "manifest.json").read_bytes() == (b / "manifest.json").read_bytes()


def test_rotated_eval_view_is_the_shifted_view_turned_r90(eval_views, synthetic_items):
    """Same items and windows as the shifted view; image and masks rotated by r90 at frame resolution (then the
    upscale), so every target rasterizes to the rotated shifted mask."""
    from objdel.train import dataset as D
    root, _ = eval_views
    its = {i.item_id: i for i in synthetic_items(("train", "val"))}
    sh, rot = (D.records(root / "out" / f"lidc-val-{k}", "lidc", "val") for k in ("shifted", "rotated"))
    assert [r["item_id"] for r in sh] == [r["item_id"] for r in rot]
    for a, b in zip(sh, rot, strict=True):
        assert (a["view"]["offset"], a["view"]["side"], b["view"]["op"]) == (b["view"]["offset"], b["view"]["side"],
                                                                             "r90")
        it, S = its[a["item_id"]], b["view"]["side"]
        cfg = FrameConfig(**{**FRAMES["lidc"].__dict__, "image_side": None})
        f = frame(it, cfg, offset=tuple(b["view"]["offset"]), side=S)
        img = np.moveaxis(augment.apply(np.moveaxis(f.image, -1, 0), "r90"), 0, -1)
        assert b["image_sha256"] == _pixels_sha(upscale(np.ascontiguousarray(img), FRAMES["lidc"].image_side))
        for r, t in enumerate(b["targets"]):
            g = to_mask(parse(t, ("nodule",)).objs["nodule"], S, S).mask
            assert np.array_equal(g, augment.apply(f.masks["nodule"][r], "r90"))


def test_build_eval_views_help():
    from experiments.phase2 import build_eval_views as BE
    with pytest.raises(SystemExit):
        BE.main(["--help"])


# --------------------------------------------------------------------------------------------- M9: fetch_lidc_pixels
def test_fetch_refuses_hu_truncation():
    from experiments.phase1.fetch_lidc_pixels import hounsfield
    assert hounsfield(np.array([[0, 3000]], np.uint16), 1.0, -1024.0).tolist() == [[-1024, 1976]]
    assert hounsfield(np.array([[0, 2]], np.uint16), 0.5, -1024.0).dtype == np.int16
    with pytest.raises(ValueError, match="int16"):
        hounsfield(np.array([[40000]], np.uint16), 1.0, 0.0)
    with pytest.raises(ValueError, match="integer"):
        hounsfield(np.array([[1]], np.uint16), 0.5, 0.0)


def test_fetch_refuses_a_shared_z_key():
    from experiments.phase1.fetch_lidc_pixels import z_keys
    assert z_keys({"a.dcm": -10.0, "b.dcm": -7.5}) == {"a.dcm": "z-10.00", "b.dcm": "z-7.50"}
    with pytest.raises(ValueError, match="z-10.00"):
        z_keys({"a.dcm": -10.0, "b.dcm": -10.001})


# ------------------------------------------------------------------------------------------------- real data
def _need(*paths):
    for p in paths:
        if not p.exists():
            pytest.skip(f"{p} not found")


def _val_items(ds, every=1):
    return [it for i, it in enumerate(store.load(ds, "val", STORE)) if i % every == 0 and it.image is not None]


@pytest.mark.data
@pytest.mark.parametrize("ds", ["lidc", "riga"])
def test_crop_before_trace_keeps_the_prepared_target_hash(ds):
    _need(STORE, PREPARED / "manifest.json")
    m = json.loads((PREPARED / "manifest.json").read_text())
    rows = []
    for it in _val_items(ds):
        f = frame(it)
        ma = min_areas(f.masks)
        for r, rater in enumerate(it.raters):
            try:
                rows.append((it.item_id, rater, build_target({n: f.masks[n][r] for n in SCHEMAS[ds]}, SCHEMAS[ds],
                                                             None, ma).text))
            except ValueError:
                pass
    assert target_set_hash(rows) == m["sets"][f"{ds}/val"]["target_set_sha256"]


@pytest.mark.data
def test_hu_context_reproduces_the_prepared_lidc_png():
    _need(STORE, PREPARED / "lidc" / "val.jsonl")
    import cv2
    recs = {json.loads(line)["item_id"]: json.loads(line) for line in open(PREPARED / "lidc" / "val.jsonl")}
    for it in _val_items("lidc", every=10):
        rec, f = recs[it.item_id], frame(it)
        assert _pixels_sha(f.image) == rec["image_sha256"]  # offset (0, 0) = Phase 1 frame and PNG hash
        ctx = _context(it)
        assert ctx["pixels"].dtype == np.int16
        x0, y0 = int(ctx["x0"]), int(ctx["y0"])
        assert np.array_equal(ctx["pixels"], it.image[y0:y0 + ctx["pixels"].shape[0], x0:x0 + ctx["pixels"].shape[1]])
        img = augment.render({"dataset": "lidc", "side": f.W, "box": list(f.box), "op": "id",
                              "window": list(LUNG_WINDOW), "noise_sigma": 0, "noise_seed": 0}, ctx)
        png = cv2.cvtColor(cv2.imread(str(PREPARED / "images" / (it.item_id + ".png"))), cv2.COLOR_BGR2RGB)
        assert np.array_equal(img, png) and _pixels_sha(img) == rec["image_sha256"]


@pytest.mark.data
def test_real_lidc_variants_are_certified_and_pixel_exact():
    """Every variant certifies and is pixel-exact (the raster is the mask), at the policy side and on jittered windows
    (BK1-1: before the exactness rung, 11 of 24,701 view rows at S = 200, where tau = 3 is 0.6 px, were not)."""
    _need(STORE)
    n, exact = 0, Counter()
    for it in _val_items("lidc", every=40):
        policy = frame(it)
        ma = min_areas(policy.masks)
        (dx, dy), S = augment.draw(it, np.random.default_rng(n))
        for kind, f in (("policy", policy), ("jittered", frame(it, offset=(dx, dy), side=S))):
            for r in range(len(it.raters)):
                for op in augment.OPS:
                    src = {"nodule": augment.apply(f.masks["nodule"][r], op)}
                    ok, ex = _certified(build_target(src, ("nodule",), None, ma).text, src, ("nodule",), ma)
                    assert ok
                    exact[kind] += ex
                    exact[kind + "/n"] += 1
                    n += 1
    assert n > 500 and exact["policy"] == exact["policy/n"] and exact["jittered"] == exact["jittered/n"]


@pytest.mark.data
def test_real_view_is_pixel_exact_and_pins_its_contexts():
    d = BUNDLE
    _need(d / "manifest.json")
    m = json.loads((d / "manifest.json").read_text())
    assert all(v["exact"] == v["rows"] and not v["error"] for v in m["counts"]["per_op"].values())
    assert m["config"]["epochs"] >= 3 and len(m["contexts"]) == m["counts"]["items"]
    assert not (d / "build_report.json").exists()


@pytest.mark.data
def test_written_view_contexts_equal_the_store():
    d = BUNDLE
    _need(d / "rows.jsonl", STORE)
    rows = _rows(d)
    want = {r["item_id"]: [] for r in rows[::2000]}
    for r in rows:
        want.get(r["item_id"], []).append(r)
    for it in store.load("lidc", "train", STORE):
        if it.item_id in want:
            c = augment.load_context(d, it.item_id)
            x0, y0 = int(c["x0"]), int(c["y0"])
            assert np.array_equal(c["pixels"], it.image[y0:y0 + c["pixels"].shape[0], x0:x0 + c["pixels"].shape[1]])
            for r in want.pop(it.item_id):
                ref = _reference(it, tuple(r["offset"]), r["side"], tuple(r["window"]), r["op"], r["noise_sigma"],
                                 r["noise_seed"])
                assert np.array_equal(augment.render(r, c), ref)
            if not want:
                break
    assert not want


def _eval_recs(kind):
    return [json.loads(line) for line in open(EVAL_VIEWS / f"lidc-val-{kind}" / "lidc" / "val.jsonl")]


@pytest.mark.data
def test_eval_view_windows_are_unchanged():
    """The pre-registered eval views are fixed (plan §5.0): augment.draw with the policy side and neighbour_side still
    give every recorded window, so a change to the training draw cannot move them."""
    _need(STORE, *(EVAL_VIEWS / f"lidc-val-{k}" / "manifest.json" for k in ("shifted", "rescaled")))
    kinds = ("shifted", "rescaled")
    seed = {k: json.loads((EVAL_VIEWS / f"lidc-val-{k}" / "manifest.json").read_text())["seed"] for k in kinds}
    want = {(k, r["item_id"]): (r["view"]["offset"], r["view"]["side"]) for k in kinds for r in _eval_recs(k)}
    got = {}
    for it in store.load("lidc", "val", STORE):
        bare = replace(it, image=None)
        S0 = frame(bare).W
        for k in kinds:
            rng = augment.row_rng(seed[k], k, it.item_id)
            w = (augment.neighbour_side(bare, S0, rng) if k == "rescaled" else None) or augment.draw(bare, rng, side=S0)
            got[k, it.item_id] = (list(w[0]), w[1])
    assert got == want


@pytest.mark.data
def test_bundle_covers_every_eval_view():
    """D10 revision: no category (side, nodule extent as share of the frame, position) that an eval view has at >= 5 %
    is under 1 % of training rows (epoch 0), and each such share is at least a third of the eval view's (R1 bundle:
    nodules under 10 % of the frame 0.5 % of rows vs 27.1 % of the rescaled val)."""
    from objdel.eval import coverage
    kinds = ("shifted", "rescaled", "rotated")
    _need(BUNDLE / "rows.jsonl", *(EVAL_VIEWS / f"lidc-val-{k}" / "lidc" / "val.jsonl" for k in kinds))
    keys = ("W", "ext", "pos")
    rows = [coverage.props(r["target"], r["side"], r["op"]) for r in _rows(BUNDLE) if r["epoch"] == 0]
    train = coverage.distribution([r for r in rows if r])
    for kind in kinds:
        exam = coverage.distribution([p for r in _eval_recs(kind) for t in r["targets"] if t
                                      if (p := coverage.props(t, r["W"], r["view"]["op"]))])
        assert not [g for g in coverage.gaps(train, exam) if g["key"] in keys], kind
        for k in keys:
            for cat, share in exam[k].items():
                assert share < coverage.EVAL_MIN or train[k].get(cat, 0) >= share / 3, (kind, k, cat)
