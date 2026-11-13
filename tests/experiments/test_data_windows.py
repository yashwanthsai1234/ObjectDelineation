"""objdel.data.windows (B1.1 window pool and cells, B1.2 scheduler; plan B.9 U03-U17, U59, U60, P01-P06, I02, F02,
R03). Synthetic LIDC items are disks on 200 x 200 slices; tests marked `data` read data/store (skipped when absent)."""
from collections import Counter

import numpy as np
import pytest
from hypothesis import event, given, settings
from hypothesis import strategies as st
from oracles.largest_remainder import lrr

from objdel.data import windows as Wn
from tests.support.shapes import annulus


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
    return annulus(n, n, cy, cx, r, hole)


def lidc_item(item_id, discs, *, other=(), group="LIDC-IDRI-0000", n=200, image=True):
    """A synthetic LIDC slice: one doctor mask per disc (cy, cx, r[, hole]), other outlined nodules as boxes."""
    from objdel.data.store import Item
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

    from objdel.data.targets.frame import frame
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

    from objdel.data.targets.frame import frame
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

    from objdel.data.targets.frame import union_box
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


def test_disagreement_is_per_slice(items30, tmp_path):  # U11: the prepare._prepare rule on the policy frame
    from objdel.data import prepare as D
    from objdel.data.targets.frame import frame
    agree = lidc_item("lidc/1.2.50/n000/z010", [(100, 100, 12)] * 3)
    split = lidc_item("lidc/1.2.50/n000/z011", [(100, 100, 12), (100, 100, 12, 5), (100, 100, 12)])
    p = Wn.pool([*items30, agree, split], draws=2)
    for it in [*items30, agree, split]:
        flags = {w.disagreement for w in p.windows if w.item_id == it.item_id}
        if flags:
            rec = D._prepare((str(tmp_path), it.item_id, ("nodule",), it.group_id, it.raters, frame(it)))
            assert flags == {rec["disagreement"]}, it.item_id
    assert {w.disagreement for w in p.windows if w.item_id == split.item_id} == {True}
    assert {w.disagreement for w in p.windows if w.item_id == agree.item_id} == {False}
    assert any(w.disagreement for w in p.windows if w.item_id not in (agree.item_id, split.item_id))


def test_draw0_matches_build_eval_views(items30, tmp_path):  # I02 on synthetic items
    from objdel.cli.data import eval_views as V
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


@pytest.fixture(scope="module")
def train_pool(store_dir):
    from objdel.data import store
    return Wn.pool(store.load("lidc", "train", store_dir), draws=1)


@pytest.mark.data
@pytest.mark.slow
def test_real_train_disagreement_count(train_pool):  # U11: 289 flagged slices (not 1,394, the nodule-level count)
    centred = [w for w in train_pool.windows if w.kind == "centred"]
    assert sum(w.disagreement for w in centred) == 289 and len(centred) == 7341


@pytest.mark.data
@pytest.mark.slow
def test_real_train_schedule(train_pool):  # R03: DG7-DG10 against the train shifted law, every block prefix
    from scipy.stats import ks_2samp
    law = [w for w in train_pool.windows if w.kind == "shifted" and w.draw == 0]
    target = Wn.cell_shares(law)
    rows = Wn.schedule(train_pool.windows, target, 100, per_step=32, seed=0)
    assert [sum(w.kind == k for w in rows) for k in Wn.KINDS] == [480, 2400, 320]
    big = target >= 0.01
    for end in range(640, 3201, 640):
        off = [w for w in rows[:end] if w.kind != "centred"]
        got = Wn.cell_shares(off)
        ratio, tv = (got[big] / target[big]).min(), 0.5 * abs(got - target).sum()
        ks = ks_2samp([w.offset for w in off], [w.offset for w in law]).statistic
        s0 = sum(w.side == w.s0 for w in off) / len(off)
        print(f"steps {end // 32}: DG7 {ratio:.3f} DG8 {tv:.4f} DG9 {ks:.4f} DG10 {s0:.4f}")
        assert ratio >= 0.8 and tv <= 0.10 and ks <= 0.06 and s0 >= 0.75


# --------------------------------------------------------------------------------------------- schedule (B1.2)
SH = {"centred": 0.15, "shifted": 0.75, "rescaled": 0.10}
EDGE_SIZE, EDGE_OFF = (5, 15, 30, 60), (0.07, 0.17, 0.25, 0.4)  # one value inside each bin


def synth_pool(rng, n_items=400, per_item=6):  # the plan's sketch: random extents and offsets
    out = []
    for i in range(n_items):
        dis = bool(rng.random() < 0.2)
        for kind in ("shifted", "rescaled"):
            for d in range(per_item // 2):
                out.append(Wn.Window(f"lidc/p{i//3}/s{i}", f"p{i//3}", f"n{i}", kind, d, 0, 0, side=100, s0=100,
                                     offset=rng.random() * 0.5, extent_px=int(rng.integers(3, 80)), disagreement=dis))
        out.append(Wn.Window(f"lidc/p{i//3}/s{i}", f"p{i//3}", f"n{i}", "centred", 0, 0, 0, side=100, s0=100,
                             offset=0.0, extent_px=20, disagreement=dis))
    return tuple(sorted(out, key=lambda w: w.window_id))


def cell_pool(n_items=2000, draws=4, dis_share=0.2, seed=0, cells=None):
    """Windows whose cells are drawn uniformly (or from `cells`): every cell well stocked for both off-centre kinds."""
    rng, out = np.random.default_rng(seed), []
    for i in range(n_items):
        item, dis = f"lidc/p{i // 3}/n{i}/z1", bool(rng.random() < dis_share)
        out.append(win("centred", item_id=item, extent_px=EDGE_SIZE[i % 4], disagreement=dis))
        for kind in ("shifted", "rescaled"):
            for d in range(draws):
                a, b = cells[i % len(cells)] if cells else rng.integers(0, 4, 2)
                out.append(win(kind, d, item_id=item, extent_px=EDGE_SIZE[a], offset=EDGE_OFF[b], disagreement=dis))
    return tuple(sorted(out, key=lambda w: w.window_id))


@pytest.fixture(scope="module")
def big_pool():
    return Wn.WindowPool(cell_pool(), {})


@pytest.fixture
def uniform_target():
    return np.full((4, 4), 1 / 16)


def slot_set(w):  # the centred slots of a block are one set
    return (w.kind, None if w.kind == "centred" else w.cell)


@settings(max_examples=40, deadline=None)
@given(blocks=st.integers(1, 3), block_steps=st.sampled_from([1, 2, 4]), per_step=st.integers(1, 32),
       seed=st.integers(0, 2**16), dw=st.sampled_from([1.0, 4.0, 16.0]))
def test_schedule_counts_and_rules(blocks, block_steps, per_step, seed, dw):  # P05, P06, P14, DG11
    rng = np.random.default_rng(seed)
    ws = synth_pool(rng)
    target = rng.dirichlet(np.ones(16)).reshape(4, 4)
    steps = blocks * block_steps
    try:
        out = Wn.schedule(ws, target, steps, per_step=per_step, shares=SH, seed=seed,
                          disagreement_weight=dw, block_steps=block_steps)
    except ValueError as e:  # unfillable cell: must name kind (and cell)
        event("raised")
        assert any(k in str(e) for k in ("shifted", "rescaled", "centred"))
        return
    n, nb = steps * per_step, block_steps * per_step
    assert len(out) == n
    use = Counter(w.window_id for w in out)  # B1.2 step 4: others once, disagreement <= 4
    assert all(c == 1 or (c <= 4 and w.disagreement) for w in out for c in [use[w.window_id]])
    pools = Counter(slot_set(w) for w in ws)
    pdis = {k: sum(slot_set(w) == k and w.disagreement for w in ws) / pools[k] for k in pools}
    for b in range(blocks):  # exact counts in every block (so in every prefix of whole blocks)
        blk = out[b * nb:(b + 1) * nb]
        want = dict(zip(Wn.KINDS, lrr(np.array([SH[k] for k in Wn.KINDS]) * nb, nb), strict=True))
        for k in Wn.KINDS:
            assert sum(w.kind == k for w in blk) == want[k]
        for k in ("shifted", "rescaled"):
            rows = [w for w in blk if w.kind == k]
            got = np.zeros(16, int)
            for w in rows:
                got[w.cell[0] * 4 + w.cell[1]] += 1
            assert got.tolist() == lrr(target.reshape(-1) * len(rows), len(rows)).tolist()
        for key, slots in Counter(slot_set(w) for w in blk).items():  # B1.2 step 3 quota, every block
            p_ = pdis[key]
            q = dw * p_ / (1 + (dw - 1) * p_)
            got_dis = sum(slot_set(w) == key and w.disagreement for w in blk)
            assert got_dis == lrr(np.array([q, 1 - q]) * slots, slots)[0]
    for s in range(steps):  # one item per step
        ids = [w.item_id for w in out[s * per_step:(s + 1) * per_step]]
        assert len(set(ids)) == per_step
    assert out == Wn.schedule(ws, target, steps, per_step=per_step, shares=SH, seed=seed,
                              disagreement_weight=dw, block_steps=block_steps)  # deterministic


def test_feasible_pool_returns_rows(big_pool, uniform_target):  # a schedule that always raises fails here
    rows = Wn.schedule(big_pool.windows, uniform_target, 20, per_step=32, seed=0)
    assert len(rows) == 640 and len({w.window_id for w in rows}) == 640


def test_steps_must_be_whole_blocks(big_pool, uniform_target):
    with pytest.raises(ValueError, match="block"):
        Wn.schedule(big_pool.windows, uniform_target, 30, per_step=32, seed=0, block_steps=20)


@pytest.mark.parametrize("steps,per_step", [(20, 32), (40, 32), (20, 5)])
def test_kind_counts_match_table(steps, per_step, big_pool, uniform_target):
    rows = Wn.schedule(big_pool.windows, uniform_target, steps, per_step=per_step, seed=0)
    got = [sum(w.kind == k for w in rows) for k in Wn.KINDS]
    nb = 20 * per_step
    assert got == list(steps // 20 * lrr(np.array([.15, .75, .10]) * nb, nb))


@pytest.mark.parametrize("steps,per_step,want", [(10, 32, (48, 240, 32)), (7, 32, (34, 168, 22)), (1, 32, (5, 24, 3)),
                                                 (1, 1, (0, 1, 0)), (3, 5, (2, 11, 2)), (31, 1, (5, 23, 3))])
def test_kind_rounding_table(steps, per_step, want, big_pool, uniform_target):  # U13 (one block of `steps`)
    rows = Wn.schedule(big_pool.windows, uniform_target, steps, per_step=per_step, block_steps=steps)
    assert tuple(sum(w.kind == k for w in rows) for k in Wn.KINDS) == want


def test_kind_rounding_ties_and_step_spread(big_pool, uniform_target):  # U13 ties; B1.2 step 1 per-step spread
    third = dict.fromkeys(Wn.KINDS, 1 / 3)
    for n, want in ((10, (4, 3, 3)), (32, (11, 11, 10))):  # ties go to the lower index
        runs = {tuple(sum(w.kind == k for w in Wn.schedule(big_pool.windows, uniform_target, n, per_step=1,
                                                             shares=third, block_steps=n)) for k in Wn.KINDS)
                for _ in range(10)}
        assert runs == {want}
    rows = Wn.schedule(big_pool.windows, uniform_target, 20, per_step=32)
    for s in range(20):
        step = Counter(w.kind for w in rows[32 * s:32 * (s + 1)])
        assert step["centred"] in (4, 5) and step["rescaled"] in (3, 4)


@pytest.mark.parametrize("n,want", [(240, [15] * 16), (241, [16] + [15] * 15), (40, [3] * 8 + [2] * 8)])
def test_cell_rounding_table(n, want, big_pool, uniform_target):  # U13 cell counts
    only = {"centred": 0.0, "shifted": 1.0, "rescaled": 0.0}
    rows = Wn.schedule(big_pool.windows, uniform_target, 1, per_step=n, shares=only, block_steps=1)
    got = Counter(w.cell for w in rows)
    assert [got[a, b] for a in range(4) for b in range(4)] == want


def test_cell_rounding_matches_oracle(big_pool):  # U13: peaked targets equal lrr
    rng = np.random.default_rng(1)
    only = {"centred": 0.0, "shifted": 1.0, "rescaled": 0.0}
    for _ in range(200):
        t = rng.dirichlet(np.full(16, 0.3))
        n = int(rng.integers(1, 200))
        rows = Wn.schedule(big_pool.windows, t.reshape(4, 4), 1, per_step=n, shares=only, block_steps=1,
                           seed=int(rng.integers(1000)))
        got = Counter(w.cell for w in rows)
        assert [got[a, b] for a in range(4) for b in range(4)] == lrr(t * n, n).tolist()


def bad_target(kind):
    t = np.full((4, 4), 1 / 16)
    if kind == "shape":
        return np.full((3, 4), 1 / 12)
    t[0, 0] = {"neg": -1e-3, "nan": np.nan, "inf": np.inf, "off7": 1 / 16 + 1e-7, "off2e9": 1 / 16 + 2e-9}[kind]
    return t


@pytest.mark.parametrize("kind", ["shape", "neg", "nan", "inf", "off7", "off2e9"])
def test_target_validation(kind, big_pool):  # U14
    with pytest.raises(ValueError):
        Wn.schedule(big_pool.windows, bad_target(kind), 1, per_step=8, block_steps=1)


@pytest.mark.parametrize("shares", [{"centred": .15, "shifted": .85}, {**SH, "edge": 0.0},
                                    {"centred": -.1, "shifted": 1.0, "rescaled": .1},
                                    {"centred": .14, "shifted": .75, "rescaled": .10},
                                    {"centred": .16, "shifted": .75, "rescaled": .10}])
def test_shares_validation(shares, big_pool, uniform_target):  # U14
    with pytest.raises(ValueError):
        Wn.schedule(big_pool.windows, uniform_target, 1, per_step=8, shares=shares, block_steps=1)


@pytest.mark.parametrize("args", [dict(steps=0), dict(per_step=0), dict(block_steps=0), dict(disagreement_weight=0.0),
                                  dict(disagreement_weight=float("nan"))])
def test_size_validation(args, big_pool, uniform_target):  # U14: no hang, no empty list
    kw = {"steps": 20, "per_step": 8, **args}
    with pytest.raises(ValueError):
        Wn.schedule(big_pool.windows, uniform_target, kw.pop("steps"), **kw)


def test_target_sum_tolerance(big_pool):  # U14: 5e-10 off is accepted
    t = np.full((4, 4), 1 / 16)
    t[0, 0] += 5e-10
    assert len(Wn.schedule(big_pool.windows, t, 1, per_step=8, block_steps=1)) == 8


def test_unfillable_cell_names_kind_and_cell_and_does_not_substitute():  # U15 (the sketch)
    rng = np.random.default_rng(0)
    ws = tuple(w for w in synth_pool(rng) if not (w.kind == "shifted" and w.cell == (1, 2)))
    with pytest.raises(ValueError, match=r"shifted.*\(?1, ?2\)?"):
        Wn.schedule(ws, np.full((4, 4), 1 / 16), 4, per_step=32, block_steps=4)


def test_unfillable_counts_and_step_rule():  # U15, P05: exactly enough succeeds, one short raises
    only = {"centred": 0.0, "shifted": 1.0, "rescaled": 0.0}
    t = np.zeros((4, 4))
    t[2, 3] = 1.0
    ws = cell_pool(32, draws=1, dis_share=0.0, cells=[(2, 3)])
    assert len(Wn.schedule(ws, t, 4, per_step=8, shares=only, block_steps=4)) == 32
    short = tuple(w for w in ws if w.item_id != "lidc/p0/n0/z1")
    with pytest.raises(ValueError, match=r"shifted.*\(2, 3\).*other"):
        Wn.schedule(short, t, 4, per_step=8, shares=only, block_steps=4)
    two = cell_pool(2, draws=8, dis_share=0.0, cells=[(2, 3)])  # 16 windows, but only 2 items for 3 rows a step
    with pytest.raises(ValueError, match="shifted"):
        Wn.schedule(two, t, 1, per_step=3, shares=only, block_steps=1)


def test_zero_mass_cells_and_centred_cells(big_pool):  # U16
    t = np.zeros((4, 4))
    t[1, 1], t[3, 0] = 0.25, 0.75
    rows = Wn.schedule(big_pool.windows, t, 20, per_step=32)
    off = Counter(w.cell for w in rows if w.kind != "centred")
    assert set(off) == {(1, 1), (3, 0)}
    t1 = np.zeros((4, 4))
    t1[0, 3] = 1.0
    assert {w.cell for w in Wn.schedule(big_pool.windows, t1, 20, per_step=32) if w.kind != "centred"} == {(0, 3)}
    rows = Wn.schedule(big_pool.windows, t1, 20, per_step=32)
    assert {w.cell for w in rows if w.kind == "centred"} == {(a, 0) for a in range(4)}  # not binned by the target


def test_schedule_does_not_mutate(big_pool, uniform_target):  # U17
    ws, t = list(big_pool.windows), uniform_target.copy()
    before = (list(ws), t.copy(), dict(Wn.SHARES))
    Wn.schedule(ws, t, 20, per_step=32)
    assert (ws, dict(Wn.SHARES)) == before[::2] and np.array_equal(t, before[1])


def test_seed_changes_the_rows_not_the_counts(big_pool, uniform_target):  # P05
    a = Wn.schedule(big_pool.windows, uniform_target, 20, per_step=32, seed=0)
    b = Wn.schedule(big_pool.windows, uniform_target, 20, per_step=32, seed=1)
    assert a != b and {w.window_id for w in a} != {w.window_id for w in b}
    assert Counter(slot_set(w) for w in a) == Counter(slot_set(w) for w in b)
    assert Wn.schedule(big_pool.windows[::-1], uniform_target, 20, per_step=32, seed=0) == a  # input order


def test_same_rows_across_processes(tmp_path):  # P04, P05: no hash() order, no global RNG
    import subprocess
    import sys
    from pathlib import Path
    code = ("import sys, hashlib; sys.path[:0] = [sys.argv[1], sys.argv[1] + '/..']; import numpy as np;"
            "import test_data_windows as T;"
            "from objdel.data import windows as Wn; p = Wn.pool(T.store30(), draws=2);"
            "rows = Wn.schedule(T.cell_pool(300), np.full((4, 4), 1 / 16), 2, per_step=16, block_steps=2);"
            "print(hashlib.sha256(repr((p, [w.window_id for w in rows])).encode()).hexdigest())")
    tests = str(Path(__file__).parent)
    out = {subprocess.run([sys.executable, "-c", code, tests], capture_output=True, text=True, check=True,
                          env={"PYTHONHASHSEED": h, "PATH": ""}).stdout for h in ("0", "1")}
    assert len(out) == 1 and len(out.pop()) == 65


def test_disagreement_weight(big_pool, uniform_target):  # P06: quota per slot set, monotone, no front-loading
    ws = cell_pool(1500, dis_share=0.4, seed=3)
    by_w = {}
    for dw in (1.0, 4.0, 16.0, 64.0):
        rows = Wn.schedule(ws, uniform_target, 40, per_step=32, disagreement_weight=dw)
        by_w[dw] = rows
        blocks = [Counter(slot_set(w) for w in rows[i:i + 640] if w.disagreement) for i in (0, 640)]
        assert all(abs(blocks[0][k] - blocks[1][k]) <= 0 for k in blocks[0] | blocks[1])  # same quota every block
    counts = [Counter(slot_set(w) for w in r) for r in by_w.values()]
    assert all(c == counts[0] for c in counts)  # the weight only moves the quota inside a slot set
    dis = [Counter(slot_set(w) for w in r if w.disagreement) for r in by_w.values()]
    for lo, hi in zip(dis, dis[1:], strict=False):
        assert all(lo[k] <= hi[k] for k in lo | hi)
    share = sum(w.disagreement for w in by_w[1.0]) / 1280
    pool_share = sum(w.disagreement for w in ws if w.kind != "centred") / sum(w.kind != "centred" for w in ws)
    assert abs(share - pool_share) < 0.05  # weight 1: the pool's share (q = p per slot set)


def test_disagreement_reuse_limit():  # B1.2 step 4: a used-up disagreement cell reuses windows, at most 4 times
    only = {"centred": 0.0, "shifted": 1.0, "rescaled": 0.0}
    t = np.zeros((4, 4))
    t[1, 1] = 1.0
    ws = cell_pool(40, draws=1, dis_share=0.25, seed=5, cells=[(1, 1)])
    n_dis = sum(w.disagreement for w in ws if w.kind == "shifted")
    rows = Wn.schedule(ws, t, 8, per_step=4, shares=only, block_steps=8, disagreement_weight=8.0)
    use = Counter(w.window_id for w in rows if w.disagreement)
    assert sum(use.values()) > n_dis and max(use.values()) <= 4
    assert all(c == 1 for w, c in Counter(w.window_id for w in rows if not w.disagreement).items())
    with pytest.raises(ValueError, match="disagreement"):
        Wn.schedule(ws, t, 40, per_step=4, shares=only, block_steps=40, disagreement_weight=64.0)


def test_no_nodule_weighting():  # I25: no nodule-weight argument
    import inspect
    assert list(inspect.signature(Wn.schedule).parameters) == ["windows", "target", "steps", "per_step", "shares",
                                                                "seed", "disagreement_weight", "block_steps"]


@pytest.mark.slow
def test_schedule_speed():  # F02: about 58k windows, 260 steps x 32; doubling the pool grows time < 3x
    import time
    t = np.full((4, 4), 1 / 16)
    secs = []
    for n in (3600, 7200):
        ws = cell_pool(n, draws=8)
        s = time.perf_counter()
        Wn.schedule(ws, t, 260, per_step=32)
        secs.append(time.perf_counter() - s)
    assert secs[0] < 5 and secs[1] < 3 * secs[0] + 0.5, secs


def test_zero_share_kind_needs_no_windows(big_pool, uniform_target):  # a kind with 0 slots asks the pool for nothing
    no_centred = tuple(w for w in big_pool.windows if w.kind != "centred")
    rows = Wn.schedule(no_centred, uniform_target, 1, per_step=8, block_steps=1,
                       shares={"centred": 0.0, "shifted": 0.9, "rescaled": 0.1})
    assert len(rows) == 8 and not any(w.kind == "centred" for w in rows)
