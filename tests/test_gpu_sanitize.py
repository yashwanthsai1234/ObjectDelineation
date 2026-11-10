"""The sanitizer corpus (Implementation.md D8): small fixed engine calls that `run_box.sh --gpu-check` runs under each
compute-sanitizer tool (memcheck, initcheck, synccheck, racecheck; `-m "gpu and sanitizer"`), each inside its 300 s cap.
Synthetic inputs from fixed seeds, no data/. tests/test_gpu_determinism.py runs the same calls 20 times."""
import numpy as np
import pytest

from experiments.cuda_checker import compare as C
from objdel import _gpucheck as gc
from tests.test_cuda_checker_gpu import _gpus

pytestmark = [pytest.mark.gpu, pytest.mark.sanitizer,
              pytest.mark.skipif(_gpus() == 0, reason="no CUDA device (nvidia-smi -L)")]


def _spiral(n: int = 33) -> np.ndarray:
    """A one-pixel square spiral with one-pixel gaps: nested rings, each opened and joined to the next."""
    m = np.zeros((n, n), bool)
    for i in range(0, n // 2, 2):
        m[i, i:n - i] = m[n - 1 - i, i:n - i] = m[i:n - i, i] = m[i:n - i, n - 1 - i] = True
        m[i + 1, i] = False  # open this ring ...
        m[i + 1, i + 1] = True  # ... and join it to the next
    return m


def check_masks() -> dict:  # 200 masks: the fixed edge cases (checkerboard among them), a spiral, random families
    return C.run_masks(C.pack_masks([m for _, m in C.fuzz_masks(191, seed=31)] + [_spiral()]))


def check_polygons() -> dict:  # at most 500 fuzz polygons on the 0-1000 grid
    cases = [(r, rl, W, H) for _, r, rl, W, H, G in C.fuzz_polygons(20, seed=32) if G == 1000][:500]
    return C.run_polygons(C.pack_polygons(cases))


def score_polygons() -> dict:  # 64 completions on one cache
    from tests.test_engine_counters import _caches, _square_args
    small, _ = _caches()  # the cache keeps its tables alive for the call
    return gc.score_polygons(*_square_args(64, 40, small.device()))


def certify_polygons() -> dict:  # two batches; the second has a rim pair
    from tests.test_certify_polygons import fuzz_objects, raster, run, square
    objs = [(parts, W, H) for parts, W, H, G in fuzz_objects(200, seed=33) if G == 1000][:40]
    refs = [raster(*o) for o in objs]
    a = run(objs, refs, np.arange(len(objs)), np.full(len(objs), 5))
    disc = [([(square(100, 100, 900, 900), [])], 64, 64), ([(square(300, 300, 700, 700), [])], 64, 64)]
    b = run(disc, [raster(*o) for o in disc], [0, 1], [0, 0], pairs=[(0, 1)], ref_pairs=[(0, 1)], pair_min_area=[0])
    return {**a, **{f"rim_{k}": v for k, v in b.items()}}


def topology() -> dict:  # 100 masks at mixed min areas
    masks = [m for _, m in C.fuzz_masks(92, seed=34)]
    p = C.pack_masks(masks)
    ma = np.random.default_rng(34).choice([0, 1, 3, 10], len(masks)).astype(np.int64)
    return gc.topology(p["bits"], p["hw"], p["start"], ma, True)


def tables() -> dict:  # 8 caches: add_tables_from_raters, then drop_tables
    from tests.test_tables_device import FAMILIES, _build
    dev = _build(FAMILIES[:8])
    for h in dev.pop("handles"):
        gc.drop_tables(int(h))
    return dev


def eval_sets() -> dict:  # 8 sets of 8 samples against caches of 4 raters
    from tests.test_eval_device import _call, _sets
    return _call(_sets(808, 8, 4, n=8), 8)


def edt() -> dict:  # 16 maps of mixed sides and densities
    rng = np.random.default_rng(35)
    maps = [(rng.random(tuple(rng.integers(1, 90, 2))) < rng.uniform(0.05, 0.95)).astype(np.uint8)
            for _ in range(16)]
    return {f"map{i}": gc.edt(m) for i, m in enumerate(maps)}


CORPUS = {f.__name__: f for f in (check_masks, check_polygons, score_polygons, certify_polygons, topology, tables,
                                  eval_sets, edt)}


@pytest.mark.parametrize("name", CORPUS)
def test_san(name):
    out = CORPUS[name]()
    assert out and all(isinstance(v, np.ndarray) for v in out.values())


def test_san_no_leak():
    """D5c leak check: after warm-up, a round of the corpus leaves live_device_bytes unchanged and equal to the arenas
    plus the tables (memcheck runs without --leak-check full: the engine is never destroyed, by design)."""
    for f in CORPUS.values():
        f()
    live = gc.stats()["live_device_bytes"]
    for f in CORPUS.values():
        f()
    st = gc.stats()
    assert st["live_device_bytes"] == live == st["arena_bytes"] + st["table_bytes"]


def test_san_pool():
    """D5c (GAP-17, D-G15 box half): after warm-up the engine's own pool holds memory and keeps it over 50 calls of
    one shape (release threshold UINT64_MAX, no allocation), the device's default pool does not grow, the device
    bytes do not leak, and trim(0) gives unused pool memory back without breaking the next call."""
    from tests.test_engine_counters import _caches, _square_args
    small, _ = _caches()
    args = _square_args(256, 40, small.device())
    want = gc.score_polygons(*args)
    before = gc.stats()
    assert before["pool_reserved"] > 0
    for _ in range(50):
        gc.score_polygons(*args)
    st = gc.stats()
    assert st["pool_reserved"] == before["pool_reserved"]
    assert st["default_pool_reserved"] <= before["default_pool_reserved"]
    assert st["live_device_bytes"] == before["live_device_bytes"] == st["arena_bytes"] + st["table_bytes"]
    assert st["sync_alloc"] == before["sync_alloc"] and st["pinned_grows"] == before["pinned_grows"]
    gc.trim(0)
    got = gc.score_polygons(*args)
    assert all(got[k].tobytes() == want[k].tobytes() for k in want)
