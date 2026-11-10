"""Expert tables built on the device (Phase D, D7b; tests I04, I06, F02, F07, S01 part and GAP-19 of Implementation.md
D.9, B18 of D.spec.md).

`add_tables_from_raters` builds each cache's tables from its rater planes in one batched engine call: union, core,
vote, the distance map (two exact EDTs), perimeter, rater areas, signatures, voters and n_vote. The oracle is the
host build of tests/oracles/reward_reference.py (scipy EDT and morphology, scipy-labelled signatures) fed to
`add_tables`. The stored tables are compared through the scorer: every candidate mask scored against both handles
gives bit-identical terms, and on small frames every single-pixel candidate is scored, so each pixel of the distance
map is seen. `ExpertCache` construction and `load` make no engine call; the first `device()` builds the tables.
"""
import os
from pathlib import Path

import numpy as np
import pytest

from objdel import _gpucheck as gc
from objdel.reward import ExpertCache
from objdel.reward import cache as rc
from objdel.reward.score import score_many
from tests.oracles.reward_reference import host_tables

DATA = Path(os.environ.get("OBJDEL_DATA_ROOT", Path(__file__).resolve().parents[1] / "data"))
KEYS = ("area", "iou_tc", "f5", "g_area", "topo_freq", "sig_start", "sig")


def _disc(h, w, cy, cx, r, hole=0.0):
    y, x = np.mgrid[:h, :w]
    d = np.hypot(y - cy, x - cx)
    return (d <= r) & (d >= hole)


def _families():
    """(raters uint8 (k, H, W), t5, min_area) of 28 synthetic caches: k in {2, 3, 4, 32}; frames 1x1, 17x33, 40x51,
    128x128, 512x512; empty, all-ones, identical and tie-voting raters; min_area in {0, 1, 3, 5, 50, 100}."""
    rng = np.random.default_rng(7)
    out = []

    def add(st, t5, ma):
        out.append((np.ascontiguousarray(st, np.uint8), t5, ma))

    for a, b in ((0, 0), (0, 1), (1, 1)):
        add(np.array([[[a]], [[b]]]), 0.5, 0)
    for k, ma in ((2, 0), (3, 1), (4, 3), (32, 5)):
        for d in (0.05, 0.5, 0.95):
            add(rng.random((k, 17, 33)) < d, 0.25 + d, ma)
    add(np.zeros((3, 40, 51)), 1.0, 5)                                      # all empty: n_vote 0, no union
    add(np.ones((2, 17, 33)), 0.75, 3)                                      # all ones: no pixel outside the core
    add(np.stack([np.zeros((17, 33)), np.ones((17, 33))]), 0.5, 0)          # one empty, one full
    add(np.stack([_disc(128, 128, 60, 64, 30), _disc(128, 128, 66, 60, 31)]), 1.0, 5)   # k = 2: ties
    add(np.stack([_disc(128, 128, 64, 64, 20)] * 3), 3.7, 50)               # identical raters
    add(np.stack([_disc(128, 128, 64, 64, 40 + i, hole=10 + 2 * i) | _disc(128, 128, 10, 10, 2 + i)
                  for i in range(4)]), 2.0, 100)                            # holes, specks below the filter
    add(rng.random((32, 128, 128)) < 0.3, 0.5, 1)
    add(np.stack([_disc(512, 512, 250 + 3 * i, 260, 120 + 5 * i, hole=30) for i in range(4)]), 1.5, 100)
    add(np.stack([rng.random((40, 51)) < 0.02, _disc(40, 51, 39, 0, 6)]), 0.6, 0)    # sources only at the edges
    return out


FAMILIES = _families()


def _candidates(raters, rng):
    """Masks scored against both tables: the raters, union, core, vote, empty, full, random and shifted masks
    (>= 50), and every single-pixel mask on frames of at most 600 pixels."""
    st = raters.astype(bool)
    k, H, W = st.shape
    c = [*st, st.any(0), st.all(0), 2 * st.sum(0) > k, np.zeros((H, W), bool), np.ones((H, W), bool)]
    c += [rng.random((H, W)) < d for d in np.linspace(0.02, 0.98, 30)]
    c += [np.roll(st[i % k], (dy, dx), (0, 1)) for i, (dy, dx) in enumerate(rng.integers(-5, 6, (15, 2)))]
    if H * W <= 600:
        c += [np.eye(1, H * W, p, dtype=bool).reshape(H, W) for p in range(H * W)]
    return c


def _score(masks, handle):
    H, W = masks[0].shape
    n = len(masks)
    bits = np.concatenate([m.astype(np.uint8).reshape(-1) for m in masks])
    return gc.score_masks(bits, np.tile(np.array([[H, W]], np.int32), (n, 1)),
                          np.arange(n + 1, dtype=np.int64) * H * W, np.full(n, handle, np.int64))


def _build(fams):
    return gc.add_tables_from_raters([r for r, _, _ in fams], np.array([t for _, t, _ in fams], np.float64),
                                     np.array([m for _, _, m in fams], np.int64))


def _check_equal(fams, dev, rng):
    """The returned scalars equal the host build's; the terms of every candidate are bit-identical on both."""
    for i, (raters, t5, ma) in enumerate(fams):
        host = host_tables(raters, t5, ma)
        a, b = dev["rater_area_start"][i], dev["rater_area_start"][i + 1]
        assert dev["perimeter"][i] == host[10] and dev["n_vote"][i] == host[9], i
        assert np.array_equal(dev["rater_area"][a:b], host[5]), i
        assert np.array_equal(dev["voter"][a:b], host[6]), i
        ss = dev["sig_start"][a:b + 1]
        assert np.array_equal(ss - ss[0], host[7]) and np.array_equal(dev["sig"][ss[0]:ss[-1]], host[8]), i
        h = gc.add_tables(*host)
        try:
            cands = _candidates(raters, rng)
            got, want = _score(cands, dev["handles"][i]), _score(cands, h)
            for key in KEYS:
                assert got[key].tobytes() == want[key].tobytes(), (i, key)
        finally:
            gc.drop_tables(h)


def test_device_tables_equal_host_tables_in_one_call_with_one_wait():
    """I04: 28 caches of mixed frames and rater counts in one call (one host wait) equal the host build."""
    gc.reset_stats()
    dev = _build(FAMILIES)
    assert gc.stats()["waits"] == 1
    m = len(FAMILIES)
    assert len(dev["handles"]) == m and dev["rater_area_start"][-1] == sum(len(r) for r, _, _ in FAMILIES)
    try:
        _check_equal(FAMILIES, dev, np.random.default_rng(1))
    finally:
        for h in dev["handles"]:
            gc.drop_tables(int(h))


def test_batch_results_do_not_depend_on_the_batch():
    """Each cache alone gives the batch's outputs."""
    whole = _build(FAMILIES)
    for i, f in enumerate(FAMILIES[::5]):
        one = _build([f])
        j = 5 * i
        a, b = whole["rater_area_start"][j], whole["rater_area_start"][j + 1]
        ss = whole["sig_start"][a:b + 1]
        assert np.array_equal(one["sig_start"], ss - ss[0]) and np.array_equal(one["sig"], whole["sig"][ss[0]:ss[-1]])
        assert (one["perimeter"][0], one["n_vote"][0]) == (whole["perimeter"][j], whole["n_vote"][j])
        assert np.array_equal(one["voter"], whole["voter"][a:b])
        gc.drop_tables(int(one["handles"][0]))
    for h in whole["handles"]:
        gc.drop_tables(int(h))


def test_table_cycles_leave_no_device_bytes():
    """I19 for the device build: after warm-up, 20 build / drop cycles leave live_device_bytes unchanged."""
    batch = FAMILIES[:12]
    for h in _build(batch)["handles"]:
        gc.drop_tables(int(h))
    live = gc.stats()["live_device_bytes"]
    for _ in range(20):
        for h in _build(batch)["handles"]:
            gc.drop_tables(int(h))
    assert gc.stats()["live_device_bytes"] == live


def test_empty_batch_makes_no_wait():
    gc.reset_stats()
    r = gc.add_tables_from_raters([], np.zeros(0), np.zeros(0, np.int64))
    assert r["handles"].dtype == np.int64 and r["handles"].size == 0 and r["sig_start"].tolist() == [0]
    assert gc.stats()["waits"] == 0


@pytest.mark.parametrize("raters, t5, ma", [
    (np.zeros((1, 4, 4), np.uint8), 1.0, 5),       # k = 1
    (np.zeros((33, 4, 4), np.uint8), 1.0, 5),      # k = 33
    (np.zeros((2, 0, 4), np.uint8), 1.0, 5),       # H = 0
    (np.zeros((2, 4, 0), np.uint8), 1.0, 5),       # W = 0
    (np.zeros((2, 4097, 1), np.uint8), 1.0, 5),    # H > 4096
    (np.zeros((2, 4, 4), np.uint8), 0.0, 5),       # t5 <= 0
    (np.zeros((2, 4, 4), np.uint8), np.nan, 5),
    (np.zeros((2, 4, 4), np.uint8), 1.0, -1),      # min_area < 0
    (np.zeros((2, 4, 4), np.uint8), 1.0, 2 ** 31),
])
def test_bad_input_raises_and_leaves_no_tables(raters, t5, ma):
    """F02, F07: a bad item anywhere in the batch raises ValueError before any upload; resident() is unchanged."""
    good = np.ones((2, 4, 4), np.uint8)
    n0 = gc.resident()
    with pytest.raises(ValueError):
        gc.add_tables_from_raters([good, raters], np.array([1.0, t5]), np.array([5, ma], np.int64))
    assert gc.resident() == n0
    with pytest.raises(ValueError):
        gc.add_tables_from_raters([good], np.array([1.0, 1.0]), np.array([5], np.int64))  # lengths differ
    assert gc.resident() == n0


def test_k2_and_k32_are_accepted():
    r = gc.add_tables_from_raters([np.ones((2, 3, 3), np.uint8), np.ones((32, 3, 3), np.uint8)], np.ones(2),
                                  np.zeros(2, np.int64))
    assert r["rater_area"].tolist() == [9] * 34
    for h in r["handles"]:
        gc.drop_tables(int(h))


# ---------------------------------------------------------------- ExpertCache on the device tables

class _Spy:
    """objdel._gpucheck with every call counted."""

    def __init__(self):
        self.calls = []

    def __getattr__(self, name):
        f = getattr(gc, name)
        if not callable(f):
            return f

        def g(*a, **k):
            self.calls.append(name)
            return f(*a, **k)
        return g


def _stack():
    return np.stack([_disc(40, 51, 20, 25, r, hole=h) for r, h in ((10, 3), (11, 0), (12, 4))])


def test_construction_and_load_make_no_engine_call(tmp_path, monkeypatch):
    """I06 (C11): constructing and loading make no engine call and leave gc.stats() unchanged; the first device()
    builds the tables (one add_tables_from_raters call), the second makes none; perimeter, rater_sigs, voters and
    same_sig call device() first; the distance maps come from cache.edt (one engine edt call each) on first host
    access, equal to scipy's EDT."""
    c0 = ExpertCache.build(_stack())
    path = tmp_path / "c.npz"
    c0.save(path)
    spy = _Spy()
    monkeypatch.setattr(rc, "_gpucheck", spy)
    gc.reset_stats()
    s0 = gc.stats()
    c = ExpertCache.load(path)
    c2 = ExpertCache(c0.raters, c0.tau, c0.t5)
    assert spy.calls == [] and gc.stats() == s0
    assert (c.k, c.n_vote, c.rater_areas, c.min_area, c.box) == (3, int(c.vote.sum()), c0.rater_areas, 5, c0.box)
    n0 = gc.resident()
    assert c.perimeter == host_tables(c.raters, c.t5, c.min_area)[10]
    assert spy.calls == ["add_tables_from_raters"] and gc.resident() == n0 + 1
    h = c.device()
    assert c.device() == h and spy.calls == ["add_tables_from_raters"]
    for name in ("rater_sigs", "voters", "same_sig", "perimeter"):
        spy.calls.clear()
        getattr(c2, name)
        assert spy.calls == ["add_tables_from_raters"], name
        c2 = ExpertCache(c0.raters, c0.tau, c0.t5)
    assert (c.rater_sigs, c.voters, c.same_sig) == (c0.rater_sigs, c0.voters, c0.same_sig)
    spy.calls.clear()
    maps = (c.dist_out_union, c.dist_in_core, c.dist_out_union, c.dist_in_core)
    assert spy.calls == ["edt", "edt"]
    from scipy import ndimage
    for got, nz in zip(maps, (~c.union, c.core, ~c.union, c.core), strict=True):
        assert got.tobytes() == ndimage.distance_transform_edt(nz).astype(np.float32).tobytes()


def test_device_many_builds_only_cold_caches_in_one_call(monkeypatch):
    """GAP-19: one score_many over 32 cold caches makes one add_tables_from_raters call and one table-build wait;
    a second call makes none; the terms equal those of caches built one by one."""
    rng = np.random.default_rng(5)
    stacks = [np.stack([_disc(32, 32, 16 + rng.integers(-3, 4), 16, 6 + r) for r in range(3)]) for _ in range(32)]
    text = '{"nodule":[{"outer":[[300,300],[700,300],[700,700],[300,700]],"holes":[]}]}'
    jobs = lambda cs: [(text if i % 4 else "not json", False, ("nodule",), {"nodule": c}, 32, 32)
                       for i, c in enumerate(cs)]
    one_by_one = [ExpertCache(s, 1.0, 1.0) for s in stacks]
    for c in one_by_one:
        c.device()
    want = score_many(jobs(one_by_one))
    spy = _Spy()
    monkeypatch.setattr(rc, "_gpucheck", spy)
    cold = [ExpertCache(s, 1.0, 1.0) for s in stacks]
    gc.reset_stats()
    got = score_many(jobs(cold))
    cold_waits = gc.stats()["waits"]
    assert spy.calls == ["add_tables_from_raters"]
    gc.reset_stats()
    assert score_many(jobs(cold)) == got == want and got[1][0]["format_ok"] == 1.0
    assert spy.calls == ["add_tables_from_raters"] and cold_waits - gc.stats()["waits"] == 1


def test_no_cpu_labelling_or_edt_in_build_or_after_load(tmp_path, monkeypatch):
    """S01 part (D7b, D7): with scipy's EDT and labelling and OpenCV's connected components patched to raise,
    ExpertCache.build gives the same tau, T5 and distance maps, and a loaded cache builds its tables and scores as
    before."""
    import cv2
    from scipy import ndimage
    c0 = ExpertCache.build(_stack())
    c0.save(tmp_path / "c.npz")
    mask = _stack()[1]
    want = rc_score(mask, c0)

    def boom(*a, **k):
        raise AssertionError("CPU geometry at runtime")
    for mod, name in ((ndimage, "distance_transform_edt"), (ndimage, "label"), (cv2, "connectedComponents"),
                      (cv2, "connectedComponentsWithStats"), (cv2, "connectedComponentsWithAlgorithm")):
        monkeypatch.setattr(mod, name, boom)
    c = ExpertCache.load(tmp_path / "c.npz")
    assert rc_score(mask, c) == want and c.rater_sigs == c0.rater_sigs and c.perimeter == c0.perimeter
    b = ExpertCache.build(_stack())
    assert (b.tau, b.t5, b.digest) == (c0.tau, c0.t5, c0.digest)
    assert np.array_equal(b.dist_out_union, c0.dist_out_union) and np.array_equal(b.dist_in_core, c0.dist_in_core)


def rc_score(mask, cache):
    from objdel.reward.score import score_mask
    return score_mask(mask, cache)


def _real_stacks(root, limit):
    files = sorted(root.rglob("*.npz"))[:limit]
    if not files:
        pytest.skip(f"no caches under {root}")
    out = []
    for f in files:
        c = ExpertCache.load(f)
        out.append((c.raters.view(np.uint8), c.t5, c.min_area))
    return out


@pytest.mark.data
@pytest.mark.slow  # about 10 min on the Mac CPU build
@pytest.mark.parametrize("ds, limit", [("lidc", 1000), ("riga", 10 ** 6)])
def test_device_tables_equal_host_tables_on_real_caches(ds, limit):
    """I04 data variant (D-G10, Mac half): the first 1,000 LIDC caches by sorted path and every RIGA cache."""
    fams = _real_stacks(DATA / "prepared" / "caches" / ds, limit)
    rng = np.random.default_rng(2)
    for i in range(0, len(fams), 64):
        batch = fams[i:i + 64]
        dev = _build(batch)
        try:
            _check_equal(batch, dev, rng)
        finally:
            for h in dev["handles"]:
                gc.drop_tables(int(h))
