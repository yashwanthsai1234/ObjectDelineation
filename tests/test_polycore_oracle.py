"""K1 and K2 (plan §7): polycore equals the slow exact oracles in tests/oracles/.

Fast subsets run by default; the full 10,000-shape K1 and 10,000-ring K2 runs are marked slow
(`pytest -m slow tests/test_polycore_oracle.py -s` prints their statistics).
"""
import collections
import random
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pytest
from oracles import raster_oracle as RO
from oracles import validity_oracle as VO

from objdel.checker import polycore as pc

WORKERS = 3  # repo rule: at most 3 worker processes


# ------------------------------------------------------------------ K1: raster
def k1_chunk(args):
    """Compare cases [lo, hi) with the oracle -> (statistics, failing case ids)."""
    lo, hi, pure, seed = args
    st, fails = collections.Counter(), []
    for i in range(lo, hi):
        fam, rings, W, H = RO.k1_case(i, seed)
        res = pc.rasterize(rings, W, H)
        w, c = RO.oracle_winding(rings, W, H, pure=pure)
        eo, nz = c % 2 == 1, w != 0
        ranges = ([(min(0, int(w.min())), max(0, int(w.max())))] if len(rings) == 1
                  else RO.ring_ranges(rings, W, H, pure=pure))
        bad = {
            "bad_evenodd_px": int((res["evenodd"] != eo).sum()),
            "bad_nonzero_px": int((res["nonzero"] != nz).sum()),
            "bad_ring_flags": sum(int(a != res["ring_wmin"][k]) + int(b != res["ring_wmax"][k])
                                  for k, (a, b) in enumerate(ranges)),
            "bad_stats": int(res["area_evenodd"] != eo.sum()) + int(res["area_nonzero"] != nz.sum())
                         + int(res["px_multi"] != (np.abs(w) >= 2).sum())
                         + int(res["rings_mixed_sign"] != sum(a < 0 < b for a, b in ranges)),
        }
        st.update(bad)
        st.update({"cases": 1, f"fam:{fam}": 1, "pixels": W * H, "eo_ne_nz_cases": int((eo != nz).any()),
                   "mixed_sign_rings": sum(a < 0 < b for a, b in ranges)})
        if i % 10 == 0:  # exact on-edge centres, on a 10 % sample (costly)
            st["sampled_on_edge_centres"] += int(RO.on_boundary_centres(rings, W, H).sum())
        if any(bad.values()):
            fails.append(i)
    return st, fails


def run_k1(n, pure=False, seed=12345, workers=1, chunk=50):
    jobs = [(a, min(a + chunk, n), pure, seed) for a in range(0, n, chunk)]
    total, fails = collections.Counter(), []
    if workers == 1:
        results = list(map(k1_chunk, jobs))
    else:
        with ProcessPoolExecutor(workers) as ex:
            results = list(ex.map(k1_chunk, jobs))
    for st, f in results:
        total.update(st)
        fails += f
    return total, fails


def check_k1(st, fails, n):
    print(dict(st))
    assert st["cases"] == n and fails == [], fails[:10]
    assert st["bad_evenodd_px"] == st["bad_nonzero_px"] == st["bad_ring_flags"] == st["bad_stats"] == 0
    assert all(st[f"fam:{f}"] >= n // 9 for f in RO.K1_FAMILIES)
    assert st["sampled_on_edge_centres"] > 0 and st["eo_ne_nz_cases"] > 0 and st["mixed_sign_rings"] > 0


def test_vectorised_oracle_is_the_pure_oracle():
    for i in range(150):
        _, rings, W, H = RO.k1_case(i, seed=99)
        W, H = min(W, 64), min(H, 37)
        w0, c0 = RO.winding_and_parity(rings, W, H)
        w1, c1 = RO.winding_np(rings, W, H)
        assert np.array_equal(w0, w1) and np.array_equal(c0, c1), i


def test_k1_subset():
    st, fails = run_k1(450)
    check_k1(st, fails, 450)


@pytest.mark.slow
def test_k1_full_10k_shapes():
    st, fails = run_k1(10_000, workers=WORKERS)
    check_k1(st, fails, 10_000)


@pytest.mark.slow
def test_k1_pure_oracle_400_shapes():
    st, fails = run_k1(400, pure=True, seed=4242, workers=WORKERS, chunk=10)
    check_k1(st, fails, 400)


@pytest.mark.parametrize("W,H", [(8, 8), (40, 25), (125, 125), (250, 100), (500, 20)])
def test_on_edge_centres_are_exercised(W, H):
    """Half-pixel lattice walks put many centres exactly on edges and vertices; still equal to the oracle."""
    rng = random.Random(W * H)
    seen = 0
    for _ in range(60):
        rings = [VO.half_pixel_walk(rng, W, H, rng.randint(3, 14)) for _ in range(rng.randint(1, 2))]
        seen += int(RO.on_boundary_centres(rings, W, H).sum())
        r = pc.rasterize(rings, W, H)
        w, c = RO.winding_np(rings, W, H)
        assert np.array_equal(r["evenodd"], c % 2 == 1) and np.array_equal(r["nonzero"], w != 0)
    assert seen >= 30


def test_huge_coordinates_exact_against_pure_oracle():
    """int32-extreme coordinates and grids take the __int128 path; compare with the pure oracle."""
    rng = random.Random(1)
    i32 = np.iinfo(np.int32)
    pick = [i32.min, i32.min + 1, -10**9, -1, 0, 1, 37, 500, 999, 1000, 10**9, i32.max]
    for _ in range(300):
        rings = [[(rng.choice(pick) if rng.random() < .5 else rng.randint(-2000, 3000),
                   rng.choice(pick) if rng.random() < .5 else rng.randint(-2000, 3000))
                  for _ in range(rng.randint(1, 7))] for _ in range(rng.randint(1, 2))]
        W, H, G = rng.randint(1, 9), rng.randint(1, 9), rng.choice([1000, 7, 2**30])
        r = pc.rasterize(rings, W, H, G=G)
        w, c = RO.winding_and_parity(rings, W, H, G)
        assert np.array_equal(r["evenodd"], c % 2 == 1) and np.array_equal(r["nonzero"], w != 0), (rings, W, H, G)
        for ring in rings:
            assert VO.labels(pc.validate([ring])["per_ring"][0]) == VO.ring_problems(ring), ring
            assert pc.signed_area2(ring) == VO.twice_area(ring)


# ------------------------------------------------------------------ K2: validity
def run_k2(n):
    st, bad = collections.Counter(), []
    for i in range(n):
        fam, ring = VO.k2_ring(i)
        res = pc.validate([ring])
        mine, ref = VO.labels(res["per_ring"][0]), VO.ring_problems(ring)
        st.update([fam] + [f"label:{lab}" for lab in ref or ["valid"]])
        if mine != ref:
            bad.append((i, ring, mine, ref))
        if "fewer_than_3_points" not in ref:
            assert res["area2"][0] == VO.twice_area(VO.clean(ring)), i
    print(dict(st))
    f = n / 10_000  # coverage on 10k: ~4.5k crossing, 3k touching, 0.3k zero area, 0.7k fewer, 3.6k valid
    assert st["label:self_crossing"] > 2000 * f and st["label:self_touching"] > 1500 * f
    assert st["label:zero_area"] > 100 * f and st["label:fewer_than_3_points"] > 300 * f
    assert st["label:valid"] > 2000 * f
    assert not bad, bad[:5]


def test_k2_subset():
    run_k2(1400)


@pytest.mark.slow
def test_k2_full_10k_rings():
    run_k2(10_000)


def run_fine_classes(n):
    """Exact per-class pair counts (crossing, pinch, overlap, spike, ...) vs the all-pairs reference."""
    for i in range(n):
        _, ring = VO.k2_ring(i, seed=31)
        ref, _, _ = VO.ref_validate([ring])
        mine = [dict(zip(pc.DEFECTS, map(int, row), strict=True)) for row in pc.validate([ring])["per_ring"]]
        assert mine == ref, (i, ring)


def run_multi_ring(n):
    rng = random.Random(5)
    seen = collections.Counter()
    for i in range(n):
        rings, roles = VO.multi_ring(rng)
        res = pc.validate(rings, roles)
        ref, pcx, ptx = VO.ref_validate(rings, roles)
        assert [dict(zip(pc.DEFECTS, map(int, row), strict=True)) for row in res["per_ring"]] == ref, (i, rings, roles)
        assert res["pair_cross"].tolist() == pcx and res["pair_touch"].tolist() == ptx
        seen.update(k for k, v in zip(pc.DEFECTS, res["total"], strict=True) if v)
    print(dict(seen))
    for k in ("ring_ring_crossing", "ring_ring_touching", "hole_outside_outer", "nested_holes", "nested_outers",
              "wrong_winding", "self_crossing"):
        assert seen[k] > n // 150, k


def test_fine_classes_subset():
    run_fine_classes(300)


def test_multi_ring_subset():
    run_multi_ring(300)


@pytest.mark.slow
def test_fine_classes_3000_rings():
    run_fine_classes(3000)


@pytest.mark.slow
def test_multi_ring_1500_shapes():
    run_multi_ring(1500)
