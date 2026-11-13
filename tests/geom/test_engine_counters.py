"""Engine counters, stage times and memory of objdel._gpucheck (Phase D, D1 and D5; tests I05, I08, I09, I14, I15, I19
and I20 of Implementation.md D.9).

`stats()` counts host waits, the drain syncs on every exit of a call, the crop rows the scorer's second label pair
relabels, the allocation calls split by kind and the growths of the pinned staging buffers after their one sizing;
it also reports the device bytes held (arenas, tables) and the pools' reserved bytes (CUDA build). The CPU build
counts the same events as the CUDA build: device memory as pool allocations (D5c: the engine's stream-ordered pool)
and host staging as synchronizing ones. `timing=True` adds per-stage times without adding a host wait.
"""
import dataclasses
import gc as pygc
import subprocess
import sys

import numpy as np
import pytest

from objdel.geom import engine as gc
from objdel.verify.reward.cache import ExpertCache
from tests.support.paths import ROOT

BASELINE = (4, 2)  # (waits, drain_syncs) per score_polygons call before D5 (gpucheck.cpp score_polygons, score_words)
EVENTS = ("waits", "drain_syncs", "label_pairs", "sync_alloc", "sync_free", "pool_allocs", "pool_frees", "pinned_grows",
          "eval_scratch_peak")
STATE = ("arena_bytes", "table_bytes", "live_device_bytes", "pool_reserved", "default_pool_reserved")
OUT_BOUND = 512 * 64 * 4  # the out staging bound (gpucheck.cpp kStageObjects x kStageComps signature entries)


def _disc(n, r):
    y, x = np.mgrid[:n, :n]
    return np.hypot(y - n / 2, x - n / 2) <= r


def _square_args(n, side, handle, a=200, b=800):
    """n copies of one square on the 0-1000 grid, scored on a side x side frame against `handle`."""
    xy = np.tile(np.array([[a, a], [b, a], [b, b], [a, b]], np.int32), (n, 1))
    W = np.full(n, side, np.int32)
    return (xy, np.arange(0, 4 * n + 1, 4, dtype=np.int64), np.arange(n + 1, dtype=np.int64), None, W, W.copy(), 1000,
            np.full(n, handle, np.int64))


def _caches():
    """Two expert caches (40 and 96 px frames)."""
    return (ExpertCache.build(np.stack([_disc(40, r) for r in (8, 9, 10)])),
            ExpertCache.build(np.stack([_disc(96, r) for r in (30, 31, 33)])))


@pytest.fixture()
def cache():
    """The two expert caches, their tables on the device for the test."""
    small, large = _caches()
    yield small, large
    del small, large


def _tables(rng, k=3, H=40, W=40):
    """add_tables arguments of a synthetic cache, made without an engine call (so a fresh process's first engine
    call can be add_tables)."""
    raters = (rng.random((k, H, W)) < 0.3).astype(np.uint8)
    uni, core = raters.max(0), raters.min(0)
    vote = (2 * raters.sum(0) > k).astype(np.uint8)
    return (raters, uni, core, vote, np.zeros((H, W), np.float32), raters.reshape(k, -1).sum(1).astype(np.int64),
            np.ones(k, np.uint8), np.arange(k + 1, dtype=np.int64), np.ones(k, np.int32), int(vote.sum()), 10, 2.0, 5)


def _delta(fn):
    gc.reset_stats()
    out = fn()
    return out, gc.stats()


def test_stats_keys_and_reset_keep_state():
    gc.check_masks(np.ones(4, np.uint8), np.array([[2, 2]], np.int32), np.array([0, 4], np.int64))
    st = gc.stats()
    assert set(st) == {*EVENTS, *STATE} and all(isinstance(v, int) and v >= 0 for v in st.values())
    gc.reset_stats()
    after = gc.stats()
    assert all(after[k] == 0 for k in EVENTS) and all(after[k] == st[k] for k in STATE)


def test_score_polygons_waits_drains_and_no_allocation_after_warm_up(cache):
    small, _ = cache
    args = _square_args(64, 40, small.device())
    want = gc.score_polygons(*args)  # warm-up: grows the arenas and the staging buffers once
    for timing in (False, True, False):
        res, st = _delta(lambda t=timing: gc.score_polygons(*args, timing=t))
        assert (st["pool_allocs"], st["pool_frees"], st["sync_alloc"], st["sync_free"]) == (0, 0, 0, 0)
        assert ("timing" in res) is timing
        assert all(res[k].tobytes() == want[k].tobytes() for k in want)
        assert (st["waits"], st["drain_syncs"]) == BASELINE  # timing reads its events after the existing waits
    times = gc.score_polygons(*args, timing=True)["timing"]
    assert list(times) == ["raster", "validity", "label", "count", "filter", "score"]
    assert all(isinstance(v, float) and v >= 0 for v in times.values())


def _filter_cases():
    """40 x 40 masks against min_area 5 (the small cache), with the crop rows (bh + 2) the size filter makes the
    second label pair relabel: removes nothing (0), removes a component, fills a hole, removes everything."""
    sq = np.zeros((40, 40), bool)
    sq[5:15, 5:15] = True  # 10 px high
    speck = sq.copy()
    speck[30, 30:32] = True  # 2 px part: dropped; box rows 5..30 = 26
    hole = sq.copy()
    hole[8:10, 8:10] = False  # 4 px hole: filled
    alone = np.zeros((40, 40), bool)
    alone[30, 30:32] = True  # all dropped; 1 row
    return [(sq, 0), (speck, 26 + 2), (hole, 10 + 2), (alone, 1 + 2)]


def _score_masks(masks, handle):
    bits = np.concatenate([m.astype(np.uint8).ravel() for m in masks])
    n = len(masks)
    return gc.score_masks(bits, np.tile(np.array([[40, 40]], np.int32), (n, 1)),
                          np.arange(0, 1600 * n + 1, 1600, dtype=np.int64), np.full(n, handle, np.int64))


def test_label_pairs_count_only_the_rows_the_size_filter_changed(cache):
    """D5 step 1 (I14): the second label pair runs per crop row only for objects the filter changed; results equal
    the host scorer in all four cases, alone and in one call."""
    from oracles import score_numpy as host

    from objdel.verify.reward.score import score_mask
    small, _ = cache
    cases = _filter_cases()
    for m, rows in cases:
        terms, st = _delta(lambda m=m: dataclasses.asdict(score_mask(m, small)))
        assert st["label_pairs"] == rows and terms == host.score_mask(m, small)
    _, st = _delta(lambda: _score_masks([m for m, _ in cases], small.device()))
    assert st["label_pairs"] == sum(r for _, r in cases)
    _, st = _delta(lambda: gc.score_polygons(*_square_args(1, 40, small.device(), a=0, b=0)))  # no crop
    assert st["label_pairs"] == 0


def test_label_skip_flags_do_not_go_stale(cache):
    """Arenas are reused and never cleared: a call whose objects were all changed, then one whose objects were not,
    gives the second call's own count and results (its flags are zeroed per call)."""
    small, _ = cache
    (sq, _), *changed = _filter_cases()
    want, st0 = _delta(lambda: _score_masks([sq] * 3, small.device()))
    _score_masks([m for m, _ in changed], small.device())
    got, st = _delta(lambda: _score_masks([sq] * 3, small.device()))
    assert st["label_pairs"] == st0["label_pairs"] == 0
    assert all(got[k].tobytes() == want[k].tobytes() for k in want)


def test_larger_call_allocates_and_its_repeat_does_not():
    """Arenas are process-wide and grow-only (no trim before D5c), so an earlier test's larger call (check_masks on a
    4096 x 4096 mask) would leave nothing to grow: the check runs in a fresh child process."""
    code = """
from tests.geom.test_engine_counters import _caches, _delta, _square_args, gc
small, large = _caches()
gc.score_polygons(*_square_args(64, 40, small.device()))
big = _square_args(256, 96, large.device())
_, st = _delta(lambda: gc.score_polygons(*big))
assert st["pool_allocs"] >= 1 and st["pool_allocs"] == st["pool_frees"]  # grow-only arenas: replace = 1 + 1
assert (st["sync_alloc"], st["sync_free"], st["pinned_grows"]) == (0, 0, 0)  # inside the pinned bound
_, st = _delta(lambda: gc.score_polygons(*big))
assert (st["pool_allocs"], st["pool_frees"], st["sync_alloc"], st["sync_free"]) == (0, 0, 0, 0)
"""
    p = subprocess.run([sys.executable, "-I", "-c", f"import sys; sys.path.insert(0, {str(ROOT)!r})\n{code}"],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]


def test_tables_count_their_bytes_and_drop_makes_no_wait():
    pygc.collect()  # drop the tables of earlier tests' garbage first, so only this cache's drop is counted
    c = ExpertCache.build(np.stack([_disc(40, r) for r in (8, 9)]))
    before = gc.stats()
    h, st = _delta(c.device)
    assert st["waits"] == 1 and st["pool_allocs"] == 1 and gc.resident() >= 1
    grown = st["table_bytes"] - before["table_bytes"]
    assert grown > 0 and st["live_device_bytes"] - before["live_device_bytes"] == grown
    gc.reset_stats()
    del c, h  # the cache's finalizer drops its tables
    pygc.collect()
    st = gc.stats()
    assert (st["waits"], st["drain_syncs"], st["pool_frees"]) == (0, 0, 1)
    assert st["table_bytes"] == before["table_bytes"]


def test_live_device_bytes_are_arenas_plus_tables(cache):
    small, _ = cache
    gc.score_polygons(*_square_args(64, 40, small.device()))
    gc.check_masks(np.ones(16, np.uint8), np.array([[4, 4]], np.int32), np.array([0, 16], np.int64))
    st = gc.stats()
    assert st["arena_bytes"] > 0 and st["live_device_bytes"] == st["arena_bytes"] + st["table_bytes"]


def test_a_call_that_throws_still_drains_both_lanes(cache):
    """No stream sync escapes the counted helper: the drain on a throwing exit counts too."""
    small, _ = cache
    gc.score_polygons(*_square_args(1, 40, small.device()))
    wrong = _square_args(1, 41, small.device())  # frame 41 against 40 px tables
    gc.reset_stats()
    with pytest.raises(ValueError, match="frame"):
        gc.score_polygons(*wrong)
    st = gc.stats()
    assert (st["waits"], st["drain_syncs"]) == (0, 2)


def test_failing_calls_leave_counters_sane(cache):
    small, _ = cache
    args = _square_args(8, 40, small.device())
    want = gc.score_polygons(*args)
    arena = gc.stats()["arena_bytes"]
    bad_hw = np.array([[3, 3]], np.int32)
    for _ in range(100):
        with pytest.raises(ValueError):
            gc.score_polygons(*_square_args(1, 41, small.device()))
        with pytest.raises(ValueError):
            gc.check_masks(np.ones(8, np.uint8), bad_hw, np.array([0, 8], np.int64))  # 3 x 3 != 8 pixels
    st = gc.stats()
    assert all(v >= 0 for v in st.values()) and st["arena_bytes"] == arena
    got = gc.score_polygons(*args)
    assert all(got[k].tobytes() == want[k].tobytes() for k in want)


def test_check_masks_and_check_polygons_stage_times():
    bits, hw, start = np.ones(64, np.uint8), np.array([[8, 8]], np.int32), np.array([0, 64], np.int64)
    want = gc.check_masks(bits, hw, start)
    got = gc.check_masks(bits, hw, start, timing=True)
    assert list(got.pop("timing")) == ["ccl", "rings", "trace"]
    assert all(got[k].tobytes() == want[k].tobytes() for k in want) and "timing" not in want
    sq = np.array([[0, 0], [9, 0], [9, 9], [0, 9]], np.int32)
    poly = (sq, np.array([0, 4], np.int64), np.array([0, 1], np.int64), None, np.array([5], np.int32),
            np.array([5], np.int32), 10)
    assert "timing" not in gc.check_polygons(*poly)
    assert list(gc.check_polygons(*poly, timing=True)["timing"]) == ["raster", "validity"]


def test_score_masks_stage_times(cache):
    small, _ = cache
    m = _disc(40, 9).astype(np.uint8).ravel()
    args = (m, np.array([[40, 40]], np.int32), np.array([0, m.size], np.int64), np.array([small.device()], np.int64))
    want = gc.score_masks(*args)
    got = gc.score_masks(*args, timing=True)
    assert list(got.pop("timing")) == ["label", "count", "filter", "score"]
    assert all(got[k].tobytes() == want[k].tobytes() for k in want)


def test_trim_is_a_noop_on_the_cpu_build_and_rejects_negative_bytes():
    st = gc.stats()
    gc.trim(0)
    if gc.build() == "cpu":
        assert gc.stats() == st
    with pytest.raises(ValueError):
        gc.trim(-1)


def test_add_tables_many_is_one_wait_and_equals_add_tables(cache):
    """I05: 5 tuples give handles that score as 5 add_tables handles; one wait per batch; an empty list makes no
    engine call; drop_tables makes no wait."""
    small, _ = cache
    rng = np.random.default_rng(1)
    args = [_tables(rng) for _ in range(5)]
    one = [gc.add_tables(*a) for a in args]
    resident = gc.resident()
    many, st = _delta(lambda: gc.add_tables_many(args))
    assert many.dtype == np.int64 and many.shape == (5,) and st["waits"] == 1 and gc.resident() == resident + 5
    for h1, h2 in zip(one, many, strict=True):
        a, b = gc.score_polygons(*_square_args(4, 40, h1)), gc.score_polygons(*_square_args(4, 40, h2))
        assert all(a[k].tobytes() == b[k].tobytes() for k in a)
    empty, st = _delta(lambda: gc.add_tables_many([]))
    assert empty.dtype == np.int64 and empty.shape == (0,) and st["waits"] == 0
    gc.reset_stats()
    for h in [*one, *many]:
        gc.drop_tables(h)
    st = gc.stats()
    assert st["waits"] == 0 and st["pool_frees"] == 10 and gc.resident() == resident - 5
    bad = list(args[0])
    bad[5] = bad[5][:2]  # rater_area of 2 for 3 raters
    with pytest.raises(ValueError):
        gc.add_tables_many([args[1], tuple(bad)])
    assert gc.resident() == resident - 5  # nothing of a failed batch stays


def test_no_device_leak_over_calls_and_table_cycles(cache):
    """I19: after warm-up, 50 score_polygons calls and 20 add_tables_many / drop_tables cycles leave the device bytes
    as they were, and they are the arenas plus the tables."""
    small, _ = cache
    rng = np.random.default_rng(2)
    args = _square_args(64, 40, small.device())
    batch = [_tables(rng) for _ in range(4)]
    gc.score_polygons(*args)
    for h in gc.add_tables_many(batch):  # warm-up of the table path
        gc.drop_tables(h)
    before, resident = gc.stats(), gc.resident()
    for _ in range(50):
        gc.score_polygons(*args)
    for _ in range(20):
        for h in gc.add_tables_many(batch):
            gc.drop_tables(h)
    st = gc.stats()
    assert st["live_device_bytes"] == before["live_device_bytes"] == st["arena_bytes"] + st["table_bytes"]
    assert st["table_bytes"] == before["table_bytes"] and gc.resident() == resident


def _fresh(code: str) -> None:
    p = subprocess.run([sys.executable, "-I", "-c", f"import sys; sys.path.insert(0, {str(ROOT)!r})\n{code}"],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-2000:]


I20 = """
import numpy as np
from tests.geom.test_engine_counters import _square_args, _tables, gc
rng = np.random.default_rng(0)
def answers(handles):  # 256 answers over the handles
    hs = np.resize(np.asarray(handles, np.int64), 256)
    a = list(_square_args(256, 40, 0))
    a[-1] = hs
    return a
def cold(size):  # D7b: a step's cold caches are built from their rater planes in one call
    st = [(rng.random((3, 40, 40)) < 0.3).astype(np.uint8) for _ in range(size)]
    return list(gc.add_tables_from_raters(st, np.full(size, 2.0), np.full(size, 5, np.int64))["handles"])
def first_tables():
    return cold(32)
if FIRST == "tables":
    hs = first_tables()
    gc.score_polygons(*answers(hs))
else:
    h0 = gc.add_tables(*_tables(rng))  # needed to score; still no staging entry before score_polygons
    gc.score_polygons(*answers([h0]))
    hs = first_tables()
base = gc.stats()
for size in (0, 1, 7, 32, 7):
    hs += cold(size)
    gc.score_polygons(*answers(hs))
st = gc.stats()
assert st["pinned_grows"] == 0, st
assert (st["sync_alloc"] - base["sync_alloc"], st["sync_free"] - base["sync_free"]) == (0, 0), st
"""


@pytest.mark.parametrize("first", ["tables", "score_polygons"])
def test_pinned_bound_does_not_depend_on_the_first_entry(first):
    """I20: the staging buffers are sized to their bounds inside the first engine call, whatever entry it is, so
    neither order of table batches and 256-answer score_polygons calls grows them afterwards."""
    _fresh(f"FIRST = {first!r}\n{I20}")


def test_a_call_above_the_out_bound_grows_the_staging_once():
    """I20 control: a score_polygons call whose defects (n x 13 x 8 B) exceed the out bound adds exactly 1 to
    pinned_grows (the counter is live), and its repeat adds nothing."""
    _fresh(f"""
import numpy as np
from tests.geom.test_engine_counters import _square_args, _tables, gc
h = gc.add_tables(*_tables(np.random.default_rng(0)))
n = {OUT_BOUND} // (13 * 8) + 1
gc.reset_stats()
gc.score_polygons(*_square_args(n, 40, h))
assert gc.stats()["pinned_grows"] == 1, gc.stats()
gc.reset_stats()
gc.score_polygons(*_square_args(n, 40, h))
assert gc.stats()["pinned_grows"] == 0, gc.stats()
""")


def test_engine_calls_from_threads_with_trim_and_reset():
    """I15: 4 threads x 50 check_masks calls on different batches while another thread calls trim(0) and
    reset_stats(): results equal the single-thread ones, no deadlock, counters stay non-negative."""
    import threading
    rng = np.random.default_rng(3)
    batches = [[rng.random((int(rng.integers(5, 60)), int(rng.integers(5, 60)))) < d for d in (.1, .5, .9)]
               for _ in range(4)]
    packed = []
    for ms in batches:
        bits = np.concatenate([m.astype(np.uint8).ravel() for m in ms])
        hw = np.array([m.shape for m in ms], np.int32)
        packed.append((bits, hw, np.cumsum([0, *(m.size for m in ms)], dtype=np.int64)))
    want = [gc.check_masks(*p) for p in packed]
    bad, stop = [], threading.Event()

    def work(k):
        for _ in range(50):
            got = gc.check_masks(*packed[k])
            if any(got[x].tobytes() != want[k][x].tobytes() for x in want[k]):
                bad.append(k)

    def chores():
        while not stop.is_set():
            gc.trim(0)
            gc.reset_stats()
            if any(v < 0 for v in gc.stats().values()):
                bad.append("stats")

    threads = [threading.Thread(target=work, args=(k,)) for k in range(4)]
    side = threading.Thread(target=chores)
    for t in (*threads, side):
        t.start()
    for t in threads:
        t.join(60)
    stop.set()
    side.join(60)
    assert not any(t.is_alive() for t in (*threads, side)) and bad == []
