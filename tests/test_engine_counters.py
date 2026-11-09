"""Engine counters and stage times of objdel._gpucheck (Phase D, D1; tests I08 and I09 of Implementation.md D.9).

`stats()` counts host waits, the drain syncs on every exit of a call, the crop rows the scorer's second label pair
relabels, and the allocation calls split by kind; it also reports the device bytes held (arenas, tables). The CPU
build counts the same events as the CUDA build: there, device memory counts as pool allocations (the D5c model) and
host staging as synchronizing ones. `timing=True` adds per-stage times without adding a host wait.
"""
import gc as pygc
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from objdel import _gpucheck as gc
from objdel.reward import ExpertCache

ROOT = Path(__file__).resolve().parents[1]
BASELINE = (4, 2)  # (waits, drain_syncs) per score_polygons call before D5 (gpucheck.cpp score_polygons, score_words)
EVENTS = ("waits", "drain_syncs", "label_pairs", "sync_alloc", "sync_free", "pool_allocs", "pool_frees")
STATE = ("arena_bytes", "table_bytes", "live_device_bytes", "contexts")


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
    assert after["contexts"] == int(gc.build() == "cuda")


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


def test_label_pairs_count_the_rows_of_the_second_label_pair(cache):
    """No skip yet (D5): the second label pair relabels every crop row, bh + 2 rows per object with pixels."""
    small, _ = cache
    _, st = _delta(lambda: gc.score_polygons(*_square_args(3, 40, small.device())))
    assert st["label_pairs"] == 3 * (24 + 2)  # the square [200, 800) of a 40 px frame is 24 px high
    empty = _square_args(1, 40, small.device(), a=0, b=0)  # zero area: no crop
    _, st = _delta(lambda: gc.score_polygons(*empty))
    assert st["label_pairs"] == 0


def test_larger_call_allocates_and_its_repeat_does_not():
    """Arenas are process-wide and grow-only (no trim before D5c), so an earlier test's larger call (check_masks on a
    4096 x 4096 mask) would leave nothing to grow: the check runs in a fresh child process."""
    code = """
from tests.test_engine_counters import _caches, _delta, _square_args, gc
small, large = _caches()
gc.score_polygons(*_square_args(64, 40, small.device()))
big = _square_args(256, 96, large.device())
_, st = _delta(lambda: gc.score_polygons(*big))
assert st["pool_allocs"] >= 1 and st["pool_allocs"] == st["pool_frees"]  # grow-only arenas: replace = 1 + 1
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
