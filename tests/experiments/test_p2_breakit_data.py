"""Break-it BK1 (Phase 2 Stage 1): the data path store -> objdel.data.targets.frame (offset / side / window) -> augment
(ops, render, contexts) -> build_views / build_eval_views -> dataset.view_rows (lazy render, A0, length, max_tokens,
hash check) -> TRL vision collator -> vast/bundle. A test marked xfail(strict) exposes a defect recorded in
docs/archive/experiments/phase2/impl/BK1.md;
every other test is a regression guard. CI tests use synthetic items; tests marked `data` read data/store and
data/views (skipped when absent)."""
import hashlib
import json
import shutil
import tarfile
from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from objdel.core import datasets
from objdel.core.hashing import file_sha256
from objdel.data import store
from objdel.data.store import Item
from objdel.data.targets import augment
from objdel.data.targets.build import BuildConfig, build_target, min_areas
from objdel.data.targets.frame import CLEARANCE, FrameConfig, centred, crowded, frame, union_box, upscale
from objdel.geom.shape import to_mask
from objdel.geom.text import parse
from objdel.geom.topology import signature
from tests.support.prepared import prepare_items

INVERSE = {"id": "id", "r90": "r270", "r180": "r180", "r270": "r90", "fh": "fh", "fv": "fv", "tr": "tr", "at": "at"}

ROOT = Path(__file__).resolve().parents[2]
NOD = ("nodule",)


def _ellipse(h, w, cy, cx, a, b=None, th=0.0):
    y, x = np.ogrid[:h, :w]
    X = (x - cx) * np.cos(th) + (y - cy) * np.sin(th)
    Y = -(x - cx) * np.sin(th) + (y - cy) * np.cos(th)
    return (X / a) ** 2 + (Y / (b or a)) ** 2 <= 1


def _shapes(h=100, w=100, cy=50, cx=50):
    """Multi-part with a hole (an L with a hole plus a blob), and a ring with an island in its hole."""
    lhole = np.zeros((h, w), bool)
    lhole[cy - 12:cy + 12, cx - 12:cx - 2] = True
    lhole[cy + 4:cy + 12, cx - 12:cx + 12] = True
    lhole[cy - 8:cy - 4, cx - 9:cx - 5] = False  # a hole in the L
    lhole |= _ellipse(h, w, cy - 8, cx + 8, 3)  # a second part
    island = _ellipse(h, w, cy, cx, 14) & ~_ellipse(h, w, cy, cx, 9) | _ellipse(h, w, cy, cx, 3)
    return lhole, island


def _lidc(i, masks, shape=(100, 100), meta=None):
    """LIDC item whose CT is bright (300 HU) exactly on rater 0's mask, -900 HU elsewhere."""
    stack = np.stack(masks)
    img = np.where(stack[0], 300, -900).astype(np.int16)
    return Item(f"lidc/1.2.{i:02d}/n000/z001", "lidc", f"LIDC-IDRI-{i:04d}", "train", img, {"nodule": stack},
                [str(r + 1) for r in range(len(masks))], meta or {})


def _items():
    lhole, island = _shapes()
    e = lambda cy, cx, a=4, h=100, w=100: _ellipse(h, w, cy, cx, a)  # noqa: E731
    its = [_lidc(1, [lhole, lhole]), _lidc(2, [island]),
           _lidc(3, [e(4, 4)]), _lidc(4, [e(95, 50)]), _lidc(5, [e(50, 96, 3)]), _lidc(6, [e(97, 97, 2)]),
           _lidc(7, [e(0, 50, 3)]),  # union on the image edge: every moved window pads
           _lidc(8, [_ellipse(260, 260, 130, 130, 45, 30, 0.4)] * 2, (260, 260)),  # S up to 200
           _lidc(9, [e(50, 50, 6)], meta={"other_nodules": [[64, 40, 70, 60]]}),  # a crowded neighbour
           _lidc(10, [e(30, 70, 5), e(30, 70, 6)]), _lidc(11, [e(70, 30, 7)]), _lidc(12, [e(20, 20, 5)])]
    return its


def _reference(it, row):
    """Offline reference of render: noise on the whole HU image, frame at frame resolution, op, Pillow upscale."""
    S, (dx, dy) = row["side"], row["offset"]
    img = it.image
    if row.get("noise_sigma"):
        x0, y0 = frame(replace(it, image=None), offset=(dx, dy), side=S).box[:2]
        noise = np.rint(augment.row_rng("noise", row["noise_seed"]).normal(0, row["noise_sigma"], (S, S)))
        img = img.astype(np.int64)
        H, W = img.shape
        ys, xs = slice(max(y0, 0), min(y0 + S, H)), slice(max(x0, 0), min(x0 + S, W))
        img[ys, xs] += noise[ys.start - y0:ys.stop - y0, xs.start - x0:xs.stop - x0].astype(np.int64)
    cfg = FrameConfig(**{**datasets.get(it.dataset).frame.__dict__, "image_side": None})
    f = frame(replace(it, image=img), cfg, offset=(dx, dy), side=S, window=row["window"] and tuple(row["window"]))
    out = np.moveaxis(augment.apply(np.moveaxis(f.image, -1, 0), row["op"]), 0, -1)
    return upscale(np.ascontiguousarray(out), datasets.get(it.dataset).frame.image_side)


def _context(it, tmp):
    augment.write_context(Path(tmp, "context", it.item_id + ".npz"), it, frame(it).box)
    return augment.load_context(tmp, it.item_id)


def _row(it, rng, **over):
    got = augment.draw(it, rng)
    assert got is not None
    (dx, dy), S = got
    ds = datasets.get(it.dataset)
    jit = augment.jitter(rng, ds) if ds.window_jitter else {"window": None, "noise_sigma": 0, "noise_seed": 0}
    ops = ds.ops
    row = {"dataset": it.dataset, "offset": [dx, dy], "side": S, "op": ops[int(rng.integers(len(ops)))], **jit}
    row |= over
    row["box"] = list(frame(replace(it, image=None), offset=tuple(row["offset"]), side=row["side"]).box)
    return row


# ------------------------------------------------------------------------------------------- frame / augment (numpy)
def test_border_windows_pad_and_render_like_the_reference(tmp_path):
    """Unions on or next to an image edge: moved windows extend past the image; render == offline reference."""
    for it in _items()[2:7]:
        ctx = _context(it, tmp_path)
        padded = 0
        for s in range(12):
            row = _row(it, np.random.default_rng(s), op=augment.OPS[s % 8])
            x0, y0, x1, y1 = row["box"]
            padded += x0 < 0 or y0 < 0 or x1 > 100 or y1 > 100
            assert np.array_equal(augment.render(row, ctx), _reference(it, row)), (it.item_id, s)
        assert padded, it.item_id


def test_offsets_at_the_clearance_limit():
    """Exactly CLEARANCE px between union and window edge is admissible on every side; one more px is not."""
    it = _items()[10]
    ux0, uy0, ux1, uy1 = union_box(it)
    for S in (50, 100, 125):
        cx, cy = centred(it, S)
        lo = (ux1 + CLEARANCE - S - cx, uy1 + CLEARANCE - S - cy)
        hi = (ux0 - CLEARANCE - cx, uy0 - CLEARANCE - cy)
        for dx, dy in [(lo[0], lo[1]), (hi[0], hi[1]), (lo[0], hi[1]), (hi[0], lo[1])]:
            m = frame(it, offset=(dx, dy), side=S).masks["nodule"].any(0)
            assert m.sum() == it.masks["nodule"].any(0).sum()  # the union is whole
            ys, xs = np.nonzero(m.any(1))[0], np.nonzero(m.any(0))[0]
            gaps = (xs[0], ys[0], S - 1 - xs[-1], S - 1 - ys[-1])
            assert min(gaps) == CLEARANCE and sorted(gaps)[1] == CLEARANCE  # two sides exactly at the limit
        for dx, dy in [(lo[0] - 1, 0), (hi[0] + 1, 0), (0, lo[1] - 1), (0, hi[1] + 1)]:
            with pytest.raises(ValueError, match="clearance"):
                frame(it, offset=(dx, dy), side=S)


def test_side_across_divisors_keeps_every_mask_pixel():
    it = _items()[7]  # 80+ px union: fits 125 up, not 50 or 100
    total = {n: s.sum((1, 2)) for n, s in it.masks.items()}
    for S in augment.LADDER:
        if S < max(union_box(it)[2] - union_box(it)[0], union_box(it)[3] - union_box(it)[1]) + 2 * CLEARANCE:
            with pytest.raises(ValueError):
                frame(it, side=S)
            continue
        f = frame(it, side=S)
        assert (f.W, f.scale, f.box[2] - f.box[0]) == (S, 1, S)
        assert all(np.array_equal(f.masks[n].sum((1, 2)), total[n]) for n in total)
    for bad in (40, 64, 600, 1000, 1):
        with pytest.raises(ValueError, match="divisor"):
            frame(it, side=bad)


@pytest.mark.parametrize("S", [50, 100, 125])
def test_ops_round_trip_targets_with_holes_and_parts(S):
    """Every D4 op x (multi-part + hole, ring + island): the target of the transformed mask rasterizes back, through
    the inverse op, to the original mask pixel for pixel (S <= 125), with the same signature."""
    lhole, island = _shapes(S, S, S // 2, S // 2)
    for m in (lhole, island):
        ma = min_areas({"nodule": m[None]})
        for op in augment.OPS:
            src = augment.apply(m, op)
            assert np.array_equal(augment.apply(src, INVERSE[op]), m)
            t = build_target({"nodule": src}, NOD, None, ma)
            g = to_mask(parse(t.text, NOD).objs["nodule"], S, S)
            assert g.legal and np.array_equal(augment.apply(g.mask, INVERSE[op]), m), op
            assert signature(g.mask, ma["nodule"]) == signature(m, ma["nodule"])
    a = np.arange(6 * 7).reshape(6, 7)  # non-square: ops and inverses still compose to the identity
    assert all(np.array_equal(augment.apply(augment.apply(a, op), INVERSE[op]), a) for op in augment.OPS)


def test_riga_nested_cup_under_hflip_stays_inside_the_disc():
    disc, cup = _ellipse(200, 200, 100, 90, 60, 50, 0.3), _ellipse(200, 200, 95, 80, 25, 20, 0.3)
    ma = min_areas({"disc": disc[None], "cup": cup[None]})
    for op in datasets.get("riga").ops:
        src = {"disc": augment.apply(disc, op), "cup": augment.apply(cup, op)}
        t = build_target(src, ("disc", "cup"), None, ma)
        p = parse(t.text, ("disc", "cup")).objs
        d, c = (to_mask(p[n], 200, 200).mask for n in ("disc", "cup"))
        assert not (c & ~d).any()  # cup inside disc


def test_hu_extremes_noise_and_window_bounds_render_like_the_reference(tmp_path):
    it = _items()[0]
    img = it.image.copy()
    img[40:45, 40:60], img[60:62, 30:70] = -32768, 32767  # int16 extremes inside the window
    it = replace(it, image=img)
    ctx = _context(it, tmp_path)
    (lc, lr), (wc, wr) = datasets.LIDC.window_jitter
    for k, (level, width) in enumerate([(lc - lr, wc - wr), (lc + lr, wc + wr), (lc - lr, wc + wr)]):
        lo = level - width // 2
        row = _row(it, np.random.default_rng(k), window=[lo, lo + width], noise_sigma=datasets.LIDC.noise_max,
                   noise_seed=2**31 - 1)
        out = augment.render(row, ctx)
        assert out.dtype == np.uint8 and np.array_equal(out, _reference(it, row))
    rng = np.random.default_rng(0)
    for _ in range(2000):
        j = augment.jitter(rng, datasets.LIDC)
        lo, hi = j["window"]
        assert abs(lo + (hi - lo) // 2 - lc) <= lr and abs(hi - lo - wc) <= wr and 0 <= j["noise_sigma"] <= 10


def test_degenerate_hu_window_is_refused():
    it = _items()[0]
    for w in [(0, 0), (400, -1000)]:
        with pytest.raises(ValueError):
            frame(it, window=w)


def test_contexts_cover_every_admissible_window(tmp_path):
    """Any window the scale + shift draw can give (including the extreme corners at the largest side) and the
    policy window lies inside the item's context crop (render never refuses a drawn row)."""
    for it in _items():
        ctx = _context(it, tmp_path)
        rng = np.random.default_rng(1)
        for _ in range(60):
            augment.render(_row(it, rng), ctx)
        ux0, uy0, ux1, uy1 = union_box(it)
        S = max(augment.sides(it))
        cx, cy = centred(it, S)
        for dx in (ux1 + CLEARANCE - S - cx, ux0 - CLEARANCE - cx):
            for dy in (uy1 + CLEARANCE - S - cy, uy0 - CLEARANCE - cy):
                if not crowded(it, (cx + dx, cy + dy, cx + dx + S, cy + dy + S)):
                    augment.render(_row(it, rng, offset=[dx, dy], side=S), ctx)


def test_small_unions_are_drawn_small_in_frame():
    """The rescaled eval view moves a small union's S from 50 to 100, where it fills <= 15 % of the frame. The margin
    draw ([8 px, 60 % of the extent]) pinned every union under 23 px at S = 50 (min_side), so training had 0.5 % of
    rows under 10 % of the frame against the rescaled val's 27 %; the side draw puts about half of them at S = 100."""
    for it in (_items()[i] for i in (2, 8, 10, 11)):  # 9-15 px unions (one crowded, one at the image corner)
        x0, y0, x1, y1 = union_box(it)
        assert frame(it).W == 50 and max(x1 - x0, y1 - y0) <= 15
        drawn = Counter(augment.draw(it, np.random.default_rng(s))[1] for s in range(400))
        assert set(drawn) == {50, 100} and 0.4 < drawn[100] / 400 < 0.6, (it.item_id, drawn)


def test_s200_targets_are_exact_with_tau_below_half_a_pixel():
    """S1-1 root cause, the fix's evidence: at S = 200 a pixel is 5 grid units, so tau = 3 is 0.6 px. Capping tau at
    0.45 px (2.25 grid units) makes the same rotated ellipses pixel-exact."""
    for seed in range(8):
        m = _s200_ellipse(seed)
        ma = min_areas({"nodule": m[None]})
        t = build_target({"nodule": m}, NOD, BuildConfig(tau=0.45 * 1000 / 200), ma)
        assert np.array_equal(to_mask(parse(t.text, NOD).objs["nodule"], 200, 200).mask, m), seed


def _s200_ellipse(seed):
    rng = np.random.default_rng(seed)
    a, b = rng.uniform(8, 25, 2)
    return _ellipse(200, 200, 100, 100, a, b, rng.uniform(0, np.pi))


def test_s200_lidc_targets_are_pixel_exact():
    m = _s200_ellipse(0)  # 6 pixels flip at tau = 3 (26 / 40 such ellipses are not exact; 0 / 40 at S = 125)
    t = build_target({"nodule": m}, NOD, None, min_areas({"nodule": m[None]}))
    assert np.array_equal(to_mask(parse(t.text, NOD).objs["nodule"], 200, 200).mask, m)


# ------------------------------------------------------------------------------------------------ build_views
@pytest.fixture(scope="module")
def bk_views(tmp_path_factory):
    from objdel.cli.data import views as BV
    root = tmp_path_factory.mktemp("bk1views")
    its = sorted(_items(), key=lambda i: i.item_id)
    store.write(its, root / "store")
    m3 = BV.build(root / "store", root / "e3" / "lidc-train", epochs=3, seed=5, procs=1)
    m1 = BV.build(root / "store", root / "e1" / "lidc-train", epochs=1, seed=5, procs=2)
    return root, {i.item_id: i for i in its}, m3, m1


def _rows(d):
    return [json.loads(line) for line in (Path(d) / "rows.jsonl").read_text().splitlines()]


def test_epoch_prefix_and_process_count_determinism(bk_views):
    """Epoch 0 of a 3-epoch, 1-process build = the 1-epoch, 2-process build (12 items > imap chunksize 8); contexts
    byte-identical; every epoch draws fresh windows."""
    root, its, m3, m1 = bk_views
    r3, r1 = _rows(root / "e3/lidc-train"), _rows(root / "e1/lidc-train")
    assert [r for r in r3 if r["epoch"] == 0] == r1
    for p in (root / "e1/lidc-train/context").rglob("*.npz"):
        assert p.read_bytes() == (root / "e3/lidc-train/context" / p.relative_to(root / "e1/lidc-train/context")
                                  ).read_bytes()
    per = {e: [(r["offset"], r["side"], r["op"]) for r in r3 if r["epoch"] == e] for e in range(3)}
    assert per[0] != per[1] != per[2] and len(per[0]) == len(per[1]) == len(per[2])
    assert m3["view_hash"] == file_sha256(root / "e3/lidc-train/rows.jsonl")


def test_row_counts_voting_raters_and_repeats(bk_views):
    root, its, m3, _ = bk_views
    want = 0
    for it in its.values():
        f = frame(replace(it, image=None))
        ma = min_areas(f.masks)
        sigs = [signature(m, ma["nodule"]) for m in f.masks["nodule"]]
        want += sum(1 if s == (0,) else 4 for s in sigs if s)
    assert m3["counts"]["rows"] == 3 * want and m3["counts"]["errors"] == 0


def test_every_row_renders_aligned_with_its_mask_and_target(bk_views):
    """End to end per row: box = frame(offset, side).box; render = offline reference; the rendered image (sampled at
    frame pixel centres) is bright exactly on rater 0's transformed mask; the target is certified against it."""
    root, its, _, _ = bk_views
    view = root / "e3/lidc-train"
    sides = set()
    for row in _rows(view):
        it = its[row["item_id"]]
        f = frame(replace(it, image=None), offset=tuple(row["offset"]), side=row["side"])
        assert list(f.box) == row["box"]
        img = augment.render(row, augment.load_context(view, row["item_id"]))
        assert np.array_equal(img, _reference(it, row)), row["row_id"]
        S = row["side"]
        sides.add(S)
        idx = ((np.arange(S) + 0.5) * img.shape[0] / S).astype(int)
        mask = augment.apply(f.masks["nodule"][row["rater_index"]], row["op"])
        if row["rater_index"] == 0:
            assert np.array_equal(img[idx][:, idx, 0] > 128, augment.apply(f.masks["nodule"][0], row["op"]))
        g = to_mask(parse(row["target"], NOD).objs["nodule"], S, S).mask
        assert np.count_nonzero(g & mask) >= 0.99 * np.count_nonzero(g | mask)
        if S <= 125:
            assert np.array_equal(g, mask), row["row_id"]
    assert {50, 100, 200} <= sides


# ------------------------------------------------------------------------------------------------ view_rows
@pytest.fixture
def D():
    pytest.importorskip("datasets")
    from objdel.sft import rows

    return rows


class Chars:
    """One token per character."""

    def __call__(self, text, add_special_tokens=False):
        return {"input_ids": list(range(len(text)))}


def test_max_tokens_boundary_counts_the_end_token(bk_views, D):
    view = bk_views[0] / "e1/lidc-train"
    full = D.view_rows(view, Chars())
    n = max(full["length"])
    assert sum(x == n for x in full["length"]) == len(full) - len(D.view_rows(view, Chars(), max_tokens=n))
    assert len(D.view_rows(view, Chars(), max_tokens=n + 1)) == len(full)


def test_max_tokens_without_a_tokenizer_is_refused(bk_views, D):
    with pytest.raises((ValueError, TypeError)):
        D.view_rows(bk_views[0] / "e1/lidc-train", None, max_tokens=1024)


def test_a0_on_multi_part_with_holes_keeps_parts_and_fills(bk_views, D):
    view = bk_views[0] / "e1/lidc-train"
    ds, a0 = D.view_rows(view), D.view_rows(view, a0=True)
    rows = {r["row_id"]: r for r in _rows(view)}
    seen = 0
    for x, y in zip(ds, a0, strict=True):
        if x["key"] != "lidc/1.2.01/n000/z001":
            continue
        S = rows[x["row_id"]]["side"]
        full = parse(x["messages"][1]["content"][0]["text"], NOD).objs["nodule"]
        filled = parse(y["messages"][1]["content"][0]["text"], NOD).objs["nodule"]
        assert len(full) == len(filled) == 2 and any(h for _, h in full) and not any(h for _, h in filled)
        g, h = to_mask(full, S, S).mask, to_mask(filled, S, S)
        assert h.legal and (g & ~h.mask).sum() == 0 and h.mask.sum() > g.sum()
        assert y["length"] < x["length"]
        seen += 1
    assert seen


def test_a0_island_inside_a_hole_stays_legal(bk_views, D):
    view = bk_views[0] / "e1/lidc-train"
    rows = {r["row_id"]: r for r in _rows(view)}
    for y in D.view_rows(view, a0=True):
        if y["key"] == "lidc/1.2.02/n000/z001":
            S = rows[y["row_id"]]["side"]
            assert to_mask(parse(y["messages"][1]["content"][0]["text"], NOD).objs["nodule"], S, S).legal


def test_view_hash_refuses_reordered_or_extended_rows(bk_views, D, tmp_path):
    for edit in (lambda t: t + "\n", lambda t: "".join(reversed(t.splitlines(True)))):
        v = tmp_path / str(hash(edit))
        shutil.copytree(bk_views[0] / "e1/lidc-train", v)
        (v / "rows.jsonl").write_text(edit((v / "rows.jsonl").read_text()))
        with pytest.raises(ValueError, match="view_hash"):
            D.view_rows(v)


def test_view_rows_refuse_an_edited_context(bk_views, D, tmp_path):
    v = tmp_path / "v"
    shutil.copytree(bk_views[0] / "e1/lidc-train", v)
    p = next((v / "context").rglob("*.npz"))
    with np.load(p) as z:
        arrays = {k: z[k] for k in z.files}
    np.savez(p.with_suffix(""), **{**arrays, "pixels": -arrays["pixels"]})  # same shape and origin, other pixels
    with pytest.raises(ValueError):
        ds = D.view_rows(v)
        [ds[i]["image"] for i in range(len(ds))]


# ------------------------------------------------------------------------------------------------ bundle
@pytest.fixture(scope="module")
def bundled(tmp_path_factory, synthetic_items):
    pytest.importorskip("datasets")
    from objdel.cli.data import views as BV
    from objdel.infra import bundle as B
    root = tmp_path_factory.mktemp("bk1bundle")
    data = root / "data"
    prepare_items(synthetic_items(("train", "val")), root, out="data/prepared")  # both datasets, as data/prepared does
    BV.build(root / "store", data / "views" / "lidc-train", procs=1)
    for n in B.LIDC_EVAL:
        (data / "eval_views" / n).mkdir(parents=True)
        (data / "eval_views" / n / "manifest.json").write_text("{}")
    (data / "checker_eval").mkdir()
    (data / "checker_eval" / "manifest.json").write_text("{}")
    tar = B.build("lidc", data, root / "out")
    box = root / "box"
    with tarfile.open(tar) as t:
        names = t.getnames()
        t.extractall(box, filter="data")
    return box, names


def test_bundled_prepared_root_preloads_on_the_box(bundled):
    from objdel.verify import store as serve
    assert serve.preload(bundled[0] / "prepared", ["lidc"]) > 0


def test_bundled_prepared_root_loads_with_cache_loader(bundled):
    from objdel.verify import store as verify_store
    verify_store.cache_loader(bundled[0] / "prepared")


def test_bundle_carries_no_wall_clock_file(bundled):
    assert not [n for n in bundled[1] if n.endswith("build_report.json")]


# ------------------------------------------------------------------------------------------------ real data
REAL = ROOT / "data" / "views" / "lidc-bundle-v3"  # the bundle trained on


def _need(*paths):
    for p in paths:
        if not p.exists():
            pytest.skip(f"{p} not found")


@pytest.mark.data
def test_real_view_rows_match_the_store_reference():
    """The first 25 store items: every row's box, render (vs the offline reference on the store's HU) and target
    certification; the manifest's hash and counts."""
    _need(REAL / "rows.jsonl", ROOT / "data" / "store" / "manifest.json")
    text = (REAL / "rows.jsonl").read_bytes()
    man = json.loads((REAL / "manifest.json").read_text())
    assert hashlib.sha256(text).hexdigest() == man["view_hash"]
    rows = [json.loads(line) for line in text.splitlines()]
    assert len(rows) == man["counts"]["rows"]
    by = {}
    for r in rows:
        by.setdefault(r["item_id"], []).append(r)
    n = 0
    for k, it in enumerate(store.load("lidc", "train", ROOT / "data" / "store")):
        if k >= 25:
            break
        bare = replace(it, image=None)
        ma = min_areas(frame(bare).masks)
        ctx = augment.load_context(REAL, it.item_id)
        for r in by.get(it.item_id, []):
            f = frame(bare, offset=tuple(r["offset"]), side=r["side"])
            assert list(f.box) == r["box"]
            assert np.array_equal(augment.render(r, ctx), _reference(it, r)), r["row_id"]
            src = augment.apply(f.masks["nodule"][r["rater_index"]], r["op"])
            g = to_mask(parse(r["target"], NOD).objs["nodule"], r["side"], r["side"])
            assert g.legal and signature(g.mask, ma["nodule"]) == signature(src, ma["nodule"])
            n += 1
    assert n > 50


@pytest.mark.data
def test_real_a0_targets_stay_legal():
    pytest.importorskip("datasets")
    from objdel.sft import rows as sft_rows
    _need(REAL / "rows.jsonl")
    for line in (REAL / "rows.jsonl").read_text().splitlines():
        r = json.loads(line)
        if not r["simple"]:
            t = sft_rows._no_holes(r["target"], datasets.get("lidc").names)
            assert to_mask(parse(t, NOD).objs["nodule"], r["side"], r["side"]).legal, r["row_id"]
