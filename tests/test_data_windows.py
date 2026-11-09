"""objdel.data.windows (B1.1 window pool and cells, B1.2 scheduler; plan B.9 U03-U17, U59, U60, P01-P06, I02, F02,
R03). Synthetic LIDC items are disks on 200 x 200 slices; tests marked `data` read data/store (skipped when absent)."""
import numpy as np
import pytest

from objdel.data import windows as Wn


def win(kind="shifted", draw=0, *, item_id="lidc/p0/n000/z001", side=100, s0=100, offset=0.0, extent_px=20,
        disagreement=False):
    return Wn.Window(item_id, "p0", item_id.rsplit("/", 1)[0], kind, draw, 0, 0, side=side, s0=s0, offset=offset,
                     extent_px=extent_px, disagreement=disagreement)


# ---------------------------------------------------------------------------------------------- cell, offset_of (U03)
@pytest.mark.parametrize("extent,off,want", [(9, .149, (0, 0)), (10, .15, (1, 1)), (19, .1999, (1, 1)),
                                             (20, .2, (2, 2)), (39, .2999, (2, 2)), (40, .3, (3, 3)),
                                             (0, 0.0, (0, 0)), (10**6, 10.0, (3, 3)), (10, 0.0, (1, 0)),
                                             (9, .3, (0, 3))])
def test_cell_edges_go_to_upper_bin(extent, off, want):
    assert Wn.cell(extent, off) == want


def test_cell_casts_float32_to_float64():
    # B.9.12 item 2: float64(float32(0.15)) = 0.15000000596... > 0.15, so it is in offset bin 1
    assert float(np.float32(0.15)) > 0.15 and Wn.cell(5, np.float32(0.15)) == (0, 1)
    assert Wn.cell(np.float32(9.999999), 0.0) == (0, 0)


@pytest.mark.parametrize("union,window,want", [((40, 40, 60, 60), (0, 0, 100, 100), 0.0),
                                               ((60, 40, 80, 60), (0, 0, 100, 100), 0.2),
                                               ((10, 10, 20, 20), (0, 0, 40, 40), 0.1767767),
                                               ((0, 0, 1, 1), (0, 0, 3, 3), 0.4714045),
                                               ((0, 0, 1, 1), (0, 0, 2, 2), 0.3535534)])
def test_offset_of_literals(union, window, want):  # U04: half-open boxes, S = the window's side
    assert Wn.offset_of(union, window) == pytest.approx(want, abs=1e-6)


def test_offset_of_properties():  # P01
    rng = np.random.default_rng(0)
    for _ in range(200):
        x0, y0 = rng.integers(0, 300, 2)
        u = (x0, y0, x0 + rng.integers(1, 40), y0 + rng.integers(1, 40))
        wx, wy, S = rng.integers(0, 300), rng.integers(0, 300), int(rng.integers(50, 200))
        w = (wx, wy, wx + S, wy + S)
        o = Wn.offset_of(u, w)
        assert o >= 0
        t = int(rng.integers(-50, 50))
        assert Wn.offset_of(tuple(v + t for v in u), tuple(v + t for v in w)) == pytest.approx(o)
        assert Wn.offset_of((u[1], u[0], u[3], u[2]), (w[1], w[0], w[3], w[2])) == pytest.approx(o)
        assert Wn.offset_of(tuple(3 * v for v in u), tuple(3 * v for v in w)) == pytest.approx(o)
    assert Wn.offset_of((10, 10, 30, 30), (0, 0, 40, 40)) == 0.0


# --------------------------------------------------------------------------------------------- Window (U05, U59)
def test_window_id_format():
    item = "lidc/LIDC-IDRI-0001/n003/z017"
    for kind, draw in (("centred", 0), ("shifted", 2), ("rescaled", 11)):
        w = win(kind, draw, item_id=item)
        assert w.window_id == item + "~" + kind + str(draw)
        assert w.window_id.rsplit("/", 1)[0] == item.rsplit("/", 1)[0]
        assert "/" not in w.window_id[len(item):]
    assert win(item_id="lidc/p~1/s").item_id == "lidc/p~1/s"


def test_apparent_size_axis():
    assert win(extent_px=30).cell[0] == 2  # shifted: side == S0, apparent = extent
    assert win("rescaled", side=200, s0=100, extent_px=30).apparent_px == 15.0
    assert win("rescaled", side=200, s0=100, extent_px=30).cell[0] == 1
    assert win("centred", extent_px=30).cell[0] == win(extent_px=30).cell[0]
    assert win(offset=0.15).cell == (2, 1)


# --------------------------------------------------------------------------------------------- cell_shares (U12)
def test_cell_shares():
    ws = [win(extent_px=5, offset=0.16)] * 3 + [win(extent_px=25, offset=0.5)]
    got = Wn.cell_shares(ws)
    assert got.shape == (4, 4) and got.dtype == np.float64
    assert got[0, 1] == 0.75 and got[2, 3] == 0.25 and abs(got.sum() - 1) < 1e-12
    assert np.array_equal(Wn.cell_shares(ws[::-1]), got)
    with pytest.raises(ValueError):
        Wn.cell_shares([])


# --------------------------------------------------------------------------------------------------- pool (B1.1)
def disk(cy, cx, r, hole=0, n=200):
    y, x = np.ogrid[:n, :n]
    d2 = (y - cy) ** 2 + (x - cx) ** 2
    return (d2 <= r * r) & (d2 > hole * hole)


def lidc_item(item_id, discs, *, other=(), group="LIDC-IDRI-0000", n=200, image=True):
    """A synthetic LIDC slice: one doctor mask per disc (cy, cx, r[, hole]), other outlined nodules as boxes."""
    from objdel.data import Item
    hu = np.random.default_rng(len(item_id)).integers(-1000, 400, (n, n), dtype=np.int16) if image else None
    return Item(item_id, "lidc", group, "train", hu, {"nodule": np.stack([disk(*d, n=n) for d in discs])},
                [str(k + 1) for k in range(len(discs))], {"other_nodules": [list(b) for b in other]})


def store30():
    """30 items of 10 patients: 1-4 doctors, radii 3-29, some single-doctor, some with a ring doctor (disagreement),
    some with a second outlined nodule beside the target, some at the slice border."""
    rng, out = np.random.default_rng(7), []
    for i in range(30):
        k = 1 if i % 7 == 3 else int(rng.integers(2, 5))
        r = int(rng.integers(3, 30))
        cy, cx = (int(v) for v in rng.integers(r + 3, 197 - r, 2))
        if i % 10 == 9:
            cx = r + 3  # union 3 px from the left border
        discs = [(cy, cx, r + j % 2) for j in range(k)]
        if i % 6 == 5 and k >= 2:
            discs[1] = (cy, cx, r + 1, max(2, r // 2))  # a ring: another topology
        other = [(cx + r + 6, cy - 5, cx + r + 14, cy + 5)] if i % 5 == 4 else ()
        out.append(lidc_item(f"lidc/1.2.{i // 3}/n00{i % 3}/z{i:03d}", discs, other=other,
                             group=f"LIDC-IDRI-{i // 3:04d}"))
    return out


@pytest.fixture(scope="module")
def items30():
    return store30()


def box_of(it, w):
    """The window's box as the renderer frames it (frame raises if the window breaks clearance or crowding)."""
    from dataclasses import replace

    from objdel.targets import frame
    bare = replace(it, image=None)
    return frame(bare).box if w.kind == "centred" else frame(bare, offset=(w.dx, w.dy), side=w.side).box


DROPS = {"excluded", "not_rl", "shifted_none", "rescaled_none", "padded", "centred_none"}


def test_pool_counts_and_keys():  # U09
    its = [lidc_item(f"lidc/1.2.0/n000/z{i:03d}", [(100, 100, 10)] * (1 if i < 4 else 2)) for i in range(10)]
    p = Wn.pool(its, draws=4)
    assert set(p.dropped) == DROPS and all(type(v) is int for v in p.dropped.values())
    assert p.dropped["not_rl"] == 4 and {w.item_id for w in p.windows} == {it.item_id for it in its[4:]}
    one = [w for w in p.windows if w.item_id == its[5].item_id]
    assert [w.draw for w in one if w.kind == "centred"] == [0]
    assert sum(w.kind == "shifted" for w in one) <= 4 and sum(w.kind == "rescaled" for w in one) <= 4
    assert len({w.window_id for w in p.windows}) == len(p.windows)
    assert Wn.pool([], draws=4) == Wn.WindowPool((), dict.fromkeys(DROPS, 0))
    p0 = Wn.pool(its, draws=0)
    assert [w.kind for w in p0.windows] == ["centred"] * 6 and sum(p0.dropped.values()) == 4
    gone = its[6].item_id
    px = Wn.pool(its, draws=4, exclude=frozenset({gone}))
    assert px.dropped["excluded"] == 1 and px.windows == tuple(w for w in p.windows if w.item_id != gone)


@pytest.mark.parametrize("draws", [0, 1, 5])
def test_pool_structure_and_geometry(items30, draws):  # P02, P03 (DG4, DG5, DG6, DG10, DG18)
    from oracles.windows_oracle import brute_admissible

    from objdel.targets import frame
    p = Wn.pool(items30, draws=draws)
    by_id = {it.item_id: it for it in items30}
    rl = [it for it in items30 if len(it.raters) >= 2]
    assert p.dropped["not_rl"] == len(items30) - len(rl)
    keys = [(w.item_id, w.kind, w.draw) for w in p.windows]
    assert len(set(keys)) == len(keys)
    for w in p.windows:
        it = by_id[w.item_id]
        S0 = frame(it).W
        assert w.s0 == S0 and w.draw < max(draws, 1)
        assert (w.side == S0) == (w.kind != "rescaled")
        assert w.nodule == w.item_id.rsplit("/", 1)[0] and w.group_id == it.group_id
        assert brute_admissible(it, box_of(it, w)), w
    for kind in ("shifted", "rescaled"):
        made = sum(w.kind == kind for w in p.windows)
        assert made + p.dropped[f"{kind}_none"] <= draws * len(rl)
    assert sum(w.kind == "centred" for w in p.windows) + p.dropped["centred_none"] <= len(rl)
    assert len(p.windows) + sum(p.dropped[k] for k in ("shifted_none", "rescaled_none", "padded", "centred_none")) \
        == len(rl) * (1 + 2 * draws)


def test_pool_window_fields_match_geometry(items30):
    from dataclasses import replace

    from objdel.targets.frame import union_box
    for w in Wn.pool(items30, draws=2).windows:
        it = next(i for i in items30 if i.item_id == w.item_id)
        u = union_box(replace(it, image=None))
        box = box_of(it, w)
        assert w.extent_px == max(u[2] - u[0], u[3] - u[1])
        assert w.offset == pytest.approx(np.hypot((u[0] + u[2] - box[0] - box[2]) / 2,
                                                  (u[1] + u[3] - box[1] - box[3]) / 2) / (box[2] - box[0]))
        assert box[2] - box[0] == w.side


def test_clearance_boundary_and_padding():  # U06, U60, U07: inside the slice, clearance exactly 2 or dropped
    two = lidc_item("lidc/1.2.0/n000/z001", [(100, 7, 5)] * 2)  # doctor pixels from column 2
    p = Wn.pool([two], draws=200)
    sh = [w for w in p.windows if w.kind == "shifted"]
    assert sh and all(box_of(two, w)[0] == 0 for w in sh)  # the only admissible column: clearance exactly 2
    assert p.dropped["padded"] > 0 and [w.kind for w in p.windows].count("centred") == 1
    one = lidc_item("lidc/1.2.0/n000/z002", [(100, 6, 5)] * 2)  # doctor pixels from column 1
    q = Wn.pool([one], draws=50)
    assert q.windows == () and q.dropped["centred_none"] == 1  # the policy frame, moved inside, has clearance 1
    assert q.dropped["padded"] + q.dropped["shifted_none"] + q.dropped["rescaled_none"] == 100
    edge = lidc_item("lidc/1.2.0/n000/z003", [(100, 5, 5)] * 2)  # touches the border (x0 = 0)
    assert Wn.pool([edge], draws=20).windows == ()


def test_large_nodule_does_not_hang():  # U07: the 400-px nodule in 512
    from oracles.windows_oracle import brute_admissible
    big = lidc_item("lidc/1.2.0/n000/z001", [(256, 256, 200)] * 2, n=512, image=False)
    p = Wn.pool([big], draws=8)
    assert all(brute_admissible(big, box_of(big, w)) for w in p.windows)
    assert len(p.windows) + sum(p.dropped.values()) == 17


def test_adjacent_nodule_is_never_inside():  # U08: a second outlined nodule 3 px from the target
    it = lidc_item("lidc/1.2.0/n000/z001", [(100, 100, 10)] * 2, other=[(113, 90, 120, 110)])
    p = Wn.pool([it], draws=30)
    assert p.dropped["centred_none"] == 1  # the policy frame holds it: checked too
    assert all(w.kind != "centred" and box_of(it, w)[2] <= 113 for w in p.windows)
    walled = lidc_item("lidc/1.2.0/n000/z002", [(100, 100, 10)] * 2, other=[(113, 0, 120, 200), (80, 0, 88, 200)])
    q = Wn.pool([walled], draws=5)
    assert q.dropped["shifted_none"] == 5 and q.windows == ()


def test_rescaled_never_keeps_the_side():  # U10
    # union 10 px: S0 = 50, its one neighbour is 100; walls 45 px out let a 50 window pass but never a 100 one
    walls = [(50, 0, 54, 200), (146, 0, 150, 200)]
    it = lidc_item("lidc/1.2.0/n000/z001", [(100, 100, 5)] * 2, other=walls)
    p = Wn.pool([it], draws=6)
    assert p.dropped["rescaled_none"] == 6 and not any(w.kind == "rescaled" for w in p.windows)
    assert sum(w.kind == "shifted" for w in p.windows) == 6
    # union 36 px: S0 = 100 (neighbours 50 and 125); walls leave one 50 column and no 100 window: shifted none,
    # rescaled at 50, and no shifted window comes in as a substitute
    tight = lidc_item("lidc/1.2.0/n000/z002", [(100, 100, 18)] * 2, other=[(71, 0, 75, 200), (125, 0, 129, 200)])
    q = Wn.pool([tight], draws=4)
    assert q.dropped["shifted_none"] == 4 and q.dropped["rescaled_none"] < 4
    assert {(w.kind, w.side) for w in q.windows} == {("rescaled", 50)}
    assert not any(hasattr(w, "fallback") or hasattr(w, "kept_side") for w in p.windows + q.windows)


def test_pool_order_seed_and_prefix(items30):  # P04
    p = Wn.pool(items30, draws=3)
    assert p.windows == tuple(sorted(p.windows, key=lambda w: w.window_id))
    rng = np.random.default_rng(0)
    for _ in range(5):
        assert Wn.pool([items30[i] for i in rng.permutation(30)], draws=3) == p
    p5 = Wn.pool(items30, draws=5)
    assert tuple(w for w in p5.windows if w.draw < 3) == p.windows
    other = {(w.item_id, w.dx, w.dy) for w in Wn.pool(items30, draws=3, seed=1).windows if w.kind == "shifted"}
    mine = {(w.item_id, w.dx, w.dy) for w in p.windows if w.kind == "shifted" and w.draw == 0}
    items = {i for i, _, _ in mine}
    assert sum((i, dx, dy) not in other for i, dx, dy in mine) >= 0.9 * len(items)
    extra = lidc_item("lidc/1.2.99/n000/z001", [(100, 100, 8)] * 3)
    assert tuple(w for w in Wn.pool([*items30, extra], draws=3).windows if w.item_id != extra.item_id) == p.windows


def test_disagreement_is_per_slice(items30, tmp_path):  # U11: the dataset._prepare rule on the policy frame
    from objdel.targets import frame
    from objdel.train import dataset as D
    agree = lidc_item("lidc/1.2.50/n000/z010", [(100, 100, 12)] * 3)
    split = lidc_item("lidc/1.2.50/n000/z011", [(100, 100, 12), (100, 100, 12, 5), (100, 100, 12)])
    p = Wn.pool([*items30, agree, split], draws=2)
    for it in [*items30, agree, split]:
        flags = {w.disagreement for w in p.windows if w.item_id == it.item_id}
        if flags:
            rec = D._prepare((str(tmp_path), it.item_id, "lidc", it.group_id, it.raters, frame(it)))
            assert flags == {rec["disagreement"]}, it.item_id
    assert {w.disagreement for w in p.windows if w.item_id == split.item_id} == {True}
    assert {w.disagreement for w in p.windows if w.item_id == agree.item_id} == {False}
    assert any(w.disagreement for w in p.windows if w.item_id not in (agree.item_id, split.item_id))


def test_draw0_matches_build_eval_views(items30, tmp_path):  # I02 on synthetic items
    from experiments.phase2 import build_eval_views as V
    p = Wn.pool(items30, draws=1)
    got = {(w.item_id, w.kind): (w.dx, w.dy, w.side) for w in p.windows if w.kind != "centred"}
    from oracles.windows_oracle import brute_admissible
    compared = 0
    for it in items30:
        if len(it.raters) < 2:
            continue
        for kind in ("shifted", "rescaled"):
            view = V._record((str(tmp_path / kind), it, kind, 0))["view"]
            if view["fallback"] or view["kept_side"]:
                assert (it.item_id, kind) not in got
                continue
            want = (*view["offset"], view["side"])
            box = box_of(it, Wn.Window(it.item_id, "", "", kind, 0, want[0], want[1], want[2], 0, 0.0, 0, False))
            if brute_admissible(it, box):
                assert got[it.item_id, kind] == want
                compared += 1
            else:  # padded in the eval view: the pool drops it
                assert (it.item_id, kind) not in got
    assert compared >= 20  # 24 of 52 on these 200-px slices; the rest are padded


# ----------------------------------------------------------------------------------------- real store (data marker)
@pytest.mark.data
def test_draw0_matches_real_eval_views(store_dir):  # I02: 200 val items against data/eval_views (seed 0)
    import json
    from itertools import islice

    from objdel.data import store
    views = store_dir.parent / "eval_views"
    if not views.is_dir():
        pytest.skip("no data/eval_views")
    its = list(islice(store.load("lidc", "val", store_dir), 200))
    got = {(w.item_id, w.kind): (w.dx, w.dy, w.side) for w in Wn.pool(its, draws=1).windows}
    ids, n = {it.item_id for it in its}, 0
    for kind in ("shifted", "rescaled"):
        for line in (views / f"lidc-val-{kind}" / "lidc" / "val.jsonl").open():
            r = json.loads(line)
            if r["item_id"] not in ids or not r["rl"]:
                continue
            v = r["view"]
            if v["fallback"] or v["kept_side"]:
                assert (r["item_id"], kind) not in got
            elif (r["item_id"], kind) in got:  # else padded: dropped by the pool
                assert got[r["item_id"], kind] == (*v["offset"], v["side"])
                n += 1
    assert n >= 300


@pytest.mark.data
@pytest.mark.slow
def test_real_train_disagreement_count(store_dir):  # U11: 289 flagged slices (not 1,394, the nodule-level count)
    from objdel.data import store
    p = Wn.pool(store.load("lidc", "train", store_dir), draws=0)
    assert sum(w.disagreement for w in p.windows) == 289 and len(p.windows) == 7341
