"""objdel.checker.polycore: API contract, decided semantics, input caps, garbage input, batch determinism, GIL."""
import hashlib
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest

from objdel import _polycore as C
from objdel.checker import polycore as pc

SQ = [(100, 100), (900, 100), (900, 900), (100, 900)]  # clockwise on screen (area > 0)
HOLE = [(400, 400), (400, 600), (600, 600), (600, 400)]  # counter-clockwise
D = {n: i for i, n in enumerate(pc.DEFECTS)}
I32 = np.iinfo(np.int32)


def V(rings, roles=None):
    """Nonzero per-ring defect counts as dicts, plus the raw result."""
    res = pc.validate(rings, roles)
    return [{k: int(v) for k, v in zip(pc.DEFECTS, row, strict=True) if v} for row in res["per_ring"]], res


# ------------------------------------------------------------------ API contract
def test_defects_and_caps_constants():
    assert pc.DEFECTS == (
        "fewer_than_3_points", "zero_area", "self_crossing", "vertex_pinch", "collinear_overlap", "spike",
        "ring_ring_crossing", "ring_ring_touching", "hole_outside_outer", "nested_holes", "nested_outers",
        "wrong_winding", "out_of_range")
    assert (pc.MAX_RINGS, pc.MAX_RING_VERTICES, pc.MAX_DIM) == (10, 1000, 4096)


def test_pack():
    xy, off = pc.pack([SQ, HOLE[:3], []])
    assert xy.dtype == np.int32 and xy.shape == (7, 2) and off.dtype == np.int32
    assert off.tolist() == [0, 4, 7, 7] and xy[4].tolist() == list(HOLE[0])
    xy, off = pc.pack([])
    assert xy.shape == (0, 2) and off.tolist() == [0]
    with pytest.raises(OverflowError):
        pc.pack([[(2**31, 0)]])
    with pytest.raises(TypeError):
        pc.pack([[(0.5, 0)]])


@pytest.mark.parametrize("ring", [np.zeros((2, 3), np.int64), np.zeros((6, 1), np.int64), np.zeros((2, 2, 2), np.int64),
                                  np.zeros(6, np.int64), np.zeros((0, 3), np.int64), [[0, 0, 1], [1, 1, 1]]])
def test_rings_must_be_k_by_2(ring):
    """B1-2: a ring is (k, 2); other shapes used to be reshaped to (-1, 2) silently."""
    with pytest.raises(ValueError, match=r"\(k, 2\)"):
        pc.pack([ring])
    with pytest.raises(ValueError):
        pc.rasterize([ring], 4, 4)
    with pytest.raises(ValueError):
        pc.signed_area2(ring)


def test_packed_offsets_and_roles_must_be_1d():
    xy, off = pc.pack([SQ, HOLE])
    with pytest.raises(ValueError):
        pc.rasterize((xy, off[:, None]), 8, 8)
    with pytest.raises(ValueError):  # two rings in a tuple are not a packed (xy, off) pair
        pc.rasterize((np.array(SQ), np.array(HOLE)), 8, 8)
    with pytest.raises(ValueError):
        pc.validate([SQ, HOLE], roles=np.array([[-1], [0]]))
    assert pc.validate([SQ, HOLE], roles=[-1, 0])["total"].sum() == 0


def test_result_keys_shapes_and_dtypes():
    r = pc.rasterize([SQ, HOLE], 20, 10)
    assert r["evenodd"].shape == r["nonzero"].shape == (10, 20) and r["evenodd"].dtype == bool
    assert r["ring_wmin"].tolist() == [0, -1] and r["ring_wmax"].tolist() == [1, 0]
    assert r["area_evenodd"] == r["evenodd"].sum() and r["area_nonzero"] == r["nonzero"].sum()
    assert r["px_multi"] == 0 and r["rings_mixed_sign"] == 0
    v = pc.validate([SQ, HOLE], [-1, 0])
    assert v["total"].shape == (13,) and v["per_ring"].shape == (2, 13) and v["total"].sum() == 0
    assert v["area2"].tolist() == [pc.signed_area2(SQ), pc.signed_area2(HOLE)] and v["n_clean"].tolist() == [4, 4]
    assert v["pair_cross"].shape == v["pair_touch"].shape == (2, 2)


def test_packed_and_list_inputs_agree():
    shape = pc.pack([SQ, HOLE])
    a, b = pc.rasterize(shape, 33, 17), pc.rasterize([SQ, HOLE], 33, 17)
    assert np.array_equal(a["evenodd"], b["evenodd"]) and np.array_equal(a["nonzero"], b["nonzero"])
    a, b = pc.validate(shape, np.array([-1, 0])), pc.validate([SQ, HOLE], [-1, 0])
    assert np.array_equal(a["per_ring"], b["per_ring"])


def test_empty_and_degenerate_shapes():
    for rings in ([], [[]], [[(5, 5)]], [[(1, 1), (9, 9)]], [[], [(3, 3)], []]):
        r = pc.rasterize(rings, 9, 9)
        assert not r["evenodd"].any() and not r["nonzero"].any()
        v = pc.validate(rings)
        assert v["per_ring"].shape == (len(rings), 13) and (v["per_ring"][:, D["fewer_than_3_points"]] == 1).all()


# ------------------------------------------------------------------ raster rules
def test_top_left_rule_square_on_centres():
    # W = H = 10: centres at 50, 150, ...; the square's edges pass exactly through centres.
    r = pc.rasterize([[(50, 50), (150, 50), (150, 150), (50, 150)]], 10, 10)
    expected = np.zeros((10, 10), bool)
    expected[0, 0] = True  # left and top edges in, right and bottom edges out
    assert np.array_equal(r["evenodd"], expected) and np.array_equal(r["nonzero"], expected)


def test_shared_edge_partition():
    """A rectangle cut by a monotone staircase path into A and B: every centre of the rectangle lands in
    exactly one of A and B, also centres lying on the shared path (no gaps, no double counting)."""
    rng = np.random.default_rng(3)
    for _ in range(200):
        W, H = int(rng.choice([7, 37, 64, 125, 250])), int(rng.choice([7, 37, 64, 125]))
        x0, x1 = sorted(rng.choice(np.arange(0, 1001, 8), 2, replace=False).tolist())
        y0, y1 = sorted(rng.choice(np.arange(0, 1001, 8), 2, replace=False).tolist())
        k = int(rng.integers(0, 6))
        xs, ys = np.sort(rng.integers(x0, x1 + 1, k)), np.sort(rng.integers(y0, y1 + 1, k))
        path = [(x0, y0)] + [(int(x), int(y)) for x, y in zip(xs, ys, strict=True)] + [(x1, y1)]
        a = pc.rasterize([path + [(x1, y0)]], W, H)["evenodd"]
        b = pc.rasterize([path + [(x0, y1)]], W, H)["evenodd"]
        rect = pc.rasterize([[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]], W, H)["evenodd"]
        assert not (a & b).any() and np.array_equal(a | b, rect)


def test_mixed_sign_flag_and_multi_winding():
    bow = [(100, 100), (900, 900), (900, 100), (100, 900)]
    r = pc.rasterize([bow], 64, 64)
    assert np.array_equal(r["evenodd"], r["nonzero"])  # even-odd vs nonzero cannot see a bow-tie
    assert (r["ring_wmin"][0], r["ring_wmax"][0], r["rings_mixed_sign"]) == (-1, 1, 1)
    r = pc.rasterize([[(100, 100), (900, 100), (900, 900), (100, 900)] * 2], 64, 64)  # wound twice
    assert r["px_multi"] > 0 and r["area_evenodd"] == 0 and r["ring_wmax"][0] == 2


# ------------------------------------------------------------------ validity: targeted cases
def test_valid_square_and_annulus():
    assert V([SQ])[0] == [{}]
    assert V([SQ, HOLE], [-1, 0])[0] == [{}, {}]


def test_wrong_winding_uses_net_area():
    """Decided (plan §4.1): a self-crossing ring is judged by its net shoelace area."""
    c, _ = V([[(0, 0), (1000, 1000), (1000, 0), (0, 1000)]])
    assert c[0] == {"self_crossing": 1, "zero_area": 1}  # net area 0: neither orientation is wrong
    c, _ = V([[(0, 0), (1000, 1000), (1000, 0), (0, 800)]])
    assert c[0] == {"self_crossing": 1, "wrong_winding": 1}  # larger lobe counter-clockwise
    c, _ = V([[(0, 800), (1000, 0), (1000, 1000), (0, 0)]])
    assert c[0] == {"self_crossing": 1}
    c, _ = V([SQ[::-1], HOLE[::-1]], [-1, 0])
    assert c == [{"wrong_winding": 1}, {"wrong_winding": 1}]


def test_hole_crossing_its_outer_gets_both_flags():
    """Decided (plan §4.1): a hole crossing its outer is both crossing and (partly) outside."""
    crossing = [(800, 400), (800, 600), (980, 600), (980, 400)]
    c, res = V([SQ, crossing], [-1, 0])
    assert c[1] == {"ring_ring_crossing": 2, "hole_outside_outer": 1}
    assert c[0] == {"ring_ring_crossing": 2} and res["pair_cross"][0, 1] == 2
    assert res["total"][D["ring_ring_crossing"]] == 2  # a pair is counted once in total
    c, _ = V([SQ, [(950, 950), (950, 990), (990, 990), (990, 950)]], [-1, 0])
    assert c[1] == {"hole_outside_outer": 1}


def test_pinch_spike_and_overlap():
    fig8 = [(0, 0), (10, 0), (5, 5), (10, 10), (0, 10), (5, 5)]
    assert set(V([fig8])[0][0]) == {"vertex_pinch"}
    c, _ = V([[(0, 0), (10, 0), (10, 10), (5, 0), (0, 10)]])  # vertex on a non-adjacent edge
    assert c[0].get("vertex_pinch", 0) > 0 and "self_crossing" not in c[0]
    assert V([[(100, 100), (500, 100), (700, 100), (500, 100), (900, 100), (900, 900), (100, 900)]])[0][0].get("spike")
    assert V([[(100, 100), (900, 100), (900, 900), (950, 950), (900, 900), (100, 900)]])[0][0].get("spike") == 1
    c, _ = V([[(0, 0), (6, 0), (6, 4), (3, 4), (3, 0), (9, 0), (9, 8), (0, 8)]])  # edges 0 and 4 overlap
    assert c[0].get("collinear_overlap", 0) >= 1


def test_hole_touching_outer_at_one_point():
    hole = [(100, 100), (300, 400), (400, 300)]
    hole = hole if pc.signed_area2(hole) < 0 else hole[::-1]
    c, res = V([SQ, hole], [-1, 0])
    assert c == [{"ring_ring_touching": 4}, {"ring_ring_touching": 4}]  # 2 x 2 edge pairs at the corner
    assert res["pair_touch"][0, 1] == 4 and res["total"][D["ring_ring_touching"]] == 4


def test_nested_holes_islands_and_outers():
    h1 = [(200, 200), (200, 800), (800, 800), (800, 200)]
    h2 = [(300, 300), (300, 700), (700, 700), (700, 300)]
    assert V([SQ, h1, h2], [-1, 0, 0])[0][1:] == [{"nested_holes": 1}, {"nested_holes": 1}]
    island = [(350, 350), (650, 350), (650, 650), (350, 650)]
    assert V([SQ, h1, island], [-1, 0, -1])[0] == [{}, {}, {}]  # an outer inside a hole is legal
    c, res = V([SQ, island])
    assert c[0] == {"nested_outers": 1} and res["total"][D["nested_outers"]] == 1


def test_out_of_range_and_degenerate_rings():
    assert V([[(-5, 0), (1000, 0), (1001, 1000), (0, 1000)]])[0][0] == {"out_of_range": 2}
    assert V([[], [(5, 5)], [(1, 1), (2, 2)], [(1, 1), (2, 2), (1, 1), (2, 2)]])[0] == [{"fewer_than_3_points": 1}] * 4
    assert V([[(0, 0), (5, 0), (10, 0)]])[0] == [{"zero_area": 1, "spike": 2}]
    assert V([[(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]])[0] == [{}]  # closing duplicate allowed


def test_signed_area2_is_exact():
    assert pc.signed_area2(SQ) == 2 * 800 * 800 and pc.signed_area2(HOLE) == -2 * 200 * 200
    big = [(I32.min, I32.min), (I32.max, I32.min), (I32.max, I32.max)]
    exact = sum(big[i][0] * big[(i + 1) % 3][1] - big[(i + 1) % 3][0] * big[i][1] for i in range(3))
    assert pc.signed_area2(big) == exact > 2**63 and pc.signed_area2(big[::-1]) == -exact
    assert pc.signed_area2([]) == 0


def test_roles_are_checked():
    bad = [([SQ, HOLE], [-1, 5]),           # parent out of range
           ([SQ, HOLE], [-2, 0]),           # negative but not -1
           ([SQ, HOLE], [-1, 1]),           # own parent
           ([SQ, HOLE], [1, 0]),            # parent is a hole
           ([SQ, HOLE, HOLE], [-1, 0, 1]),  # parent is a hole
           ([SQ, HOLE], [-1])]              # one role per ring
    for rings, roles in bad:
        with pytest.raises(ValueError):
            pc.validate(rings, roles)
    assert pc.validate([SQ, HOLE], np.array([-1, 0]))["total"].sum() == 0


# ------------------------------------------------------------------ caps and garbage input
def _ring(n, cx=500, cy=500, r=400):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.c_[cx + r * np.cos(t), cy + r * np.sin(t)].round().astype(np.int64)


def test_caps_are_accepted_at_the_limit():
    rings = [_ring(1000, 100 + 80 * k, 500, 30) for k in range(10)]
    r = pc.rasterize(rings, 4096, 4096)
    assert r["evenodd"].shape == (4096, 4096) and r["area_evenodd"] > 0
    assert pc.validate(rings)["per_ring"].shape == (10, 13)
    assert pc.score_batch([rings], [4096], [4096], threads=2)["defects"].shape == (1, 13)


@pytest.mark.parametrize("case", ["rings", "vertices", "W", "H", "W0", "G", "G0"])
def test_beyond_caps_raise_value_error(case):
    rings, W, H, G = [SQ], 64, 64, 1000
    if case == "rings":
        rings = [SQ] * 11
    elif case == "vertices":
        rings = [_ring(1001)]
    elif case in ("W", "H", "W0"):
        W, H = {"W": (4097, 64), "H": (64, 4097), "W0": (0, 64)}[case]
    else:
        G = 2**30 + 1 if case == "G" else 0
    with pytest.raises(ValueError, match="10 rings|1000 vertices|4096|2\\^30"):
        pc.rasterize(rings, W, H, G)
    with pytest.raises(ValueError):
        pc.score_batch([[SQ], rings], [64, W], [64, H], G=G)
    if case in ("rings", "vertices", "G", "G0"):
        with pytest.raises(ValueError):
            pc.validate(rings, G=G)


def test_malformed_raw_inputs_raise():
    xy, good = np.zeros((4, 2), np.int32), np.array([0, 4], np.int32)
    bad = [
        lambda: C.rasterize(xy, np.array([0, 5], np.int32), 8, 8, 1000),   # offsets past the end
        lambda: C.rasterize(xy, np.array([2, 1], np.int32), 8, 8, 1000),   # decreasing
        lambda: C.rasterize(xy, np.array([-1, 2], np.int32), 8, 8, 1000),  # negative
        lambda: C.rasterize(xy, np.array([], np.int32), 8, 8, 1000),       # no offsets
        lambda: C.rasterize(np.zeros((2, 3), np.int32), good, 8, 8, 1000),  # not (N, 2)
        lambda: C.rasterize(xy.astype(np.float64), good, 8, 8, 1000),      # wrong dtype
        lambda: C.rasterize(xy.astype(np.int64), good, 8, 8, 1000),        # no silent narrowing
        lambda: C.validate(xy, good, np.array([0], np.int32), 1000),       # own parent
        lambda: C.score_batch([xy], [good, good], [8], [8], None, None, 1000, 2),  # length mismatch
        lambda: pc.score_batch([[SQ]], [8], [8], gts=[np.zeros((7, 8), bool)]),    # gt shape
    ]
    for f in bad:
        with pytest.raises((ValueError, TypeError)):
            f()


def test_random_garbage_never_crashes():
    rng = np.random.default_rng(0)
    n_ok = n_err = 0
    for _ in range(3000):
        n = int(rng.integers(0, 30))
        mode = int(rng.integers(3))
        lo, hi = [(I32.min, I32.max), (-3, 4), (-5000, 6000)][mode]
        xy = rng.integers(lo, hi, (n, 2), dtype=np.int64, endpoint=True).astype(np.int32)
        R = int(rng.integers(0, 13))
        off = (np.sort(rng.integers(0, n + 1, R + 1)) if rng.random() < .8 else rng.integers(-3, n + 4, R + 1))
        off = off.astype(np.int32)
        W, H = int(rng.integers(-2, 70)), int(rng.integers(-2, 70))
        G = int(rng.choice([0, 1, 7, 1000, 2**30, 2**31 - 1]))
        roles = rng.integers(-2, R + 1, R).astype(np.int32) if R and rng.random() < .3 else None
        for f, args in ((C.rasterize, (xy, off, W, H, G)), (C.validate, (xy, off, roles, G)),
                        (C.score_batch, ([xy], [off], [W], [H], [roles], None, G, 2))):
            try:
                f(*args)
                n_ok += 1
            except (ValueError, OverflowError, TypeError):
                n_err += 1
    assert n_ok > 1000 and n_err > 1000


# ------------------------------------------------------------------ batch, determinism, GIL
def _blob(rng, n, cx, cy, r, ccw=False):
    k = np.arange(1, 5)
    a, ph = rng.normal(0, 0.25 / k), rng.uniform(0, 2 * np.pi, k.size)
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    rad = np.clip(r * (1 + (a[:, None] * np.cos(k[:, None] * t + ph[:, None])).sum(0)), 0.3 * r, None)
    pts = np.clip(np.rint(np.c_[cx + rad * np.cos(t), cy + rad * np.sin(t)]), 0, 1000).astype(np.int64)
    pts = pts[np.r_[True, np.any(np.diff(pts, axis=0) != 0, axis=1)]]
    return pts[::-1] if ccw else pts


def rollouts(n=96, seed=7):
    """GRPO-like rollouts: 1-3 smooth rings (holes or second parts), 16-200 vertices, crops 64-512 px."""
    rng = np.random.default_rng(seed)
    shapes, sizes, roles, gts = [], [], [], []
    for _ in range(n):
        S = int(rng.integers(64, 513))
        cx, cy = rng.integers(350, 651, 2)
        R = float(rng.uniform(150, 330))
        rings, rl = [_blob(rng, int(rng.integers(16, 201)), cx, cy, R)], [-1]
        for _ in range(int(rng.integers(0, 3))):
            if rng.random() < 0.5:
                rings.append(_blob(rng, int(rng.integers(16, 201)), cx, cy, R * rng.uniform(.15, .3), ccw=True))
                rl.append(0)
            else:
                rings.append(_blob(rng, int(rng.integers(16, 201)), *rng.integers(100, 901, 2), R * .2))
                rl.append(-1)
        shapes.append(pc.pack(rings))
        sizes.append(S)
        roles.append(np.array(rl, np.int32))
        gts.append(pc.rasterize([_blob(rng, 64, cx + 9, cy - 7, R)], S, S)["evenodd"])
    return shapes, sizes, roles, gts


ROLL = rollouts()


def _digest(d):
    h = hashlib.sha256()
    for k in sorted(d):
        h.update(k.encode() + np.ascontiguousarray(d[k]).tobytes())
    return h.hexdigest()


def test_score_batch_equals_single_calls():
    shapes, S, roles, gts = ROLL
    b = pc.score_batch(shapes, S, S, roles, gts, threads=3)
    for i, (shape, s, rl, gt) in enumerate(zip(shapes, S, roles, gts, strict=True)):
        r, v = pc.rasterize(shape, s, s), pc.validate(shape, rl)
        assert b["area_evenodd"][i] == r["area_evenodd"] and b["area_nonzero"][i] == r["area_nonzero"]
        assert b["px_multi"][i] == r["px_multi"] and b["rings_mixed_sign"][i] == r["rings_mixed_sign"]
        assert b["inter_evenodd"][i] == (r["evenodd"] & gt).sum() and b["inter_nonzero"][i] == (r["nonzero"] & gt).sum()
        assert b["gt_area"][i] == gt.sum() and np.array_equal(b["defects"][i], v["total"])
    assert pc.score_batch([], [], [])["defects"].shape == (0, 13)
    as_u8 = pc.score_batch(shapes, S, S, roles, [g.astype(np.uint8) * 255 for g in gts])  # any nonzero is inside
    as_i64 = pc.score_batch(shapes, S, S, roles, [g.astype(np.int64) for g in gts])
    assert all(np.array_equal(b[k], x[k]) for x in (as_u8, as_i64) for k in b)
    nogt = pc.score_batch(shapes, S, S)
    assert (nogt["gt_area"] == 0).all() and np.array_equal(nogt["area_evenodd"], b["area_evenodd"])


def _batch(threads):
    shapes, S, roles, gts = ROLL
    return _digest(pc.score_batch(shapes, S, S, roles, gts, threads=threads))


def test_score_batch_deterministic_across_thread_counts():
    ref = _batch(1)
    assert all(_batch(t) == ref for t in (0, 2, 3, 8, 64))
    with ThreadPoolExecutor(3) as ex:  # concurrent Python callers share no state
        assert set(ex.map(_batch, [1, 2, 3] * 3)) == {ref}


def test_score_batch_deterministic_across_processes():
    code = "import sys; sys.path.insert(0, sys.argv[1]); import test_polycore as t; print(t._batch(int(sys.argv[2])))"
    here = str(Path(__file__).parent)
    run = [subprocess.run([sys.executable, "-c", code, here, t], capture_output=True, text=True, check=True)
           for t in ("1", "3")]
    outs = {r.stdout.strip() for r in run}
    assert outs == {_batch(1)}


_OOM = """
import resource, sys
import numpy as np
from objdel.checker import polycore as pc
# 10 rings x 1,000 vertices zigzagging the full height: ~5 M crossing keys (41 MB) per worker's scratch
shape = pc.pack([np.c_[np.arange(1000) * (r + 1) % 1001, np.arange(1000) % 2 * 1000] for r in range(10)])
vm = int(open("/proc/self/statm").read().split()[0]) * resource.getpagesize()
resource.setrlimit(resource.RLIMIT_AS, (vm + (24 << 20), resource.RLIM_INFINITY))
try:
    pc.score_batch([shape] * 8, [4096] * 8, [4096] * 8, threads=int(sys.argv[1]))
except MemoryError:
    print("MemoryError")
finally:  # the interpreter is not robust at the limit (e.g. at shutdown)
    resource.setrlimit(resource.RLIMIT_AS, (resource.RLIM_INFINITY, resource.RLIM_INFINITY))
"""


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="needs an enforced RLIMIT_AS (Linux)")
@pytest.mark.parametrize("threads", [1, 4, 64])
def test_score_batch_out_of_memory_raises_in_python(threads):
    """I8: std::bad_alloc in a worker thread used to call std::terminate (process abort); it must reach Python as
    MemoryError, as it already did inline (threads = 1)."""
    # No OpenBLAS thread servers: OpenCV's bundled OpenBLAS 0.3.15 jumps to a null pointer when a server thread's
    # lazy buffer allocation fails (blas_memory_alloc, seen in a core dump), which is not what this test is about.
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "1"}
    r = subprocess.run([sys.executable, "-c", _OOM, str(threads)], capture_output=True, text=True, timeout=120, env=env)
    assert r.returncode == 0 and r.stdout.strip() == "MemoryError", (r.returncode, r.stdout, r.stderr[-2000:])


def test_kernels_release_the_gil():
    """While one thread is inside a long kernel, another Python thread keeps running at full speed."""
    rng = np.random.default_rng(1)
    shape = pc.pack([np.c_[rng.integers(0, 1001, 1000), np.arange(1000) % 2 * 1000] for _ in range(10)])
    ticks, done = [0], threading.Event()

    def spin():
        while not done.is_set():
            ticks[0] += 1

    th = threading.Thread(target=spin)
    th.start()
    try:
        c0 = ticks[0]
        time.sleep(0.1)  # the main thread sleeps without the GIL: the spinner's full rate
        rate = (ticks[0] - c0) / 0.1
        c1, t0 = ticks[0], time.perf_counter()
        pc.rasterize(shape, 2048, 2048)
        dt, during = time.perf_counter() - t0, ticks[0] - c1
    finally:
        done.set()
        th.join()
    assert dt > 0.05 and during > 0.3 * rate * dt, (dt, during, rate)
