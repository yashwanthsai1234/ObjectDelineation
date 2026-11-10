"""K1 and K2 (plan §7) on the CUDA checker: objdel.geom.engine.check_polygons equals the slow exact oracles in
tests/oracles/ (pure Python, no code shared with any checker). These were the CPU checker's oracle tests
(archive/cpu_checker/tests/test_polycore_oracle.py); the checker returns the even-odd mask and defect totals, so
those are compared. Fast subsets run by default; the full 10,000-shape K1 and 10,000-ring K2 runs are marked slow.
"""
import collections
import random

import numpy as np
import pytest
from oracles import raster_oracle as RO
from oracles import validity_oracle as VO

from experiments.cuda_checker import compare as C
from objdel.geom import engine as _gpucheck


def check(cases, G=1000) -> list[dict]:
    """[(rings, roles or None, W, H)] -> per case {defects (13,), evenodd (H, W)}, one checker call."""
    p = C.pack_polygons(cases)
    out = C.run_polygons(p, G)
    return [C.gpu_polygon(out, p, s) for s in range(len(cases))]


def test_defect_names_are_the_oracles():
    assert _gpucheck.DEFECTS == VO.DEFECTS


# ------------------------------------------------------------------ K1: raster
def run_k1(n, pure=False, seed=12345):
    st, fails = collections.Counter(), []
    for lo in range(0, n, 50):
        cases = [RO.k1_case(i, seed) for i in range(lo, min(lo + 50, n))]
        for i, (fam, rings, W, H), g in zip(range(lo, n), cases, check([(r, None, W, H) for _, r, W, H in cases]),
                                            strict=False):
            _, c = RO.oracle_winding(rings, W, H, pure=pure)
            bad = int((g["evenodd"] != (c % 2 == 1)).sum())
            st.update({"cases": 1, f"fam:{fam}": 1, "pixels": W * H, "bad_evenodd_px": bad})
            if i % 10 == 0:  # exact on-edge centres, on a 10 % sample (costly)
                st["sampled_on_edge_centres"] += int(RO.on_boundary_centres(rings, W, H).sum())
            if bad:
                fails.append(i)
    print(dict(st))
    assert st["cases"] == n and fails == [], fails[:10]
    assert all(st[f"fam:{f}"] >= n // 9 for f in RO.K1_FAMILIES) and st["sampled_on_edge_centres"] > 0


def test_k1_subset():
    run_k1(450)


@pytest.mark.slow
def test_k1_full_10k_shapes():
    run_k1(10_000)


@pytest.mark.slow
def test_k1_pure_oracle_400_shapes():
    run_k1(400, pure=True, seed=4242)


@pytest.mark.parametrize("W,H", [(8, 8), (40, 25), (125, 125), (250, 100), (500, 20)])
def test_on_edge_centres_are_exercised(W, H):
    """Half-pixel lattice walks put many centres exactly on edges and vertices; still equal to the oracle."""
    rng = random.Random(W * H)
    shapes = [[VO.half_pixel_walk(rng, W, H, rng.randint(3, 14)) for _ in range(rng.randint(1, 2))] for _ in range(60)]
    assert sum(int(RO.on_boundary_centres(rings, W, H).sum()) for rings in shapes) >= 30
    for rings, g in zip(shapes, check([(r, None, W, H) for r in shapes]), strict=True):
        assert np.array_equal(g["evenodd"], RO.winding_np(rings, W, H)[1] % 2 == 1)


def test_huge_coordinates_exact_against_pure_oracle():
    """int32-extreme coordinates and grids take the __int128 paths; compare with the pure oracles."""
    rng = random.Random(1)
    i32 = np.iinfo(np.int32)
    pick = [i32.min, i32.min + 1, -10**9, -1, 0, 1, 37, 500, 999, 1000, 10**9, i32.max]
    for _ in range(300):
        rings = [[(rng.choice(pick) if rng.random() < .5 else rng.randint(-2000, 3000),
                   rng.choice(pick) if rng.random() < .5 else rng.randint(-2000, 3000))
                  for _ in range(rng.randint(1, 7))] for _ in range(rng.randint(1, 2))]
        W, H, G = rng.randint(1, 9), rng.randint(1, 9), rng.choice([1000, 7, 2**30])
        g, = check([(rings, None, W, H)], G)
        assert np.array_equal(g["evenodd"], RO.winding_and_parity(rings, W, H, G)[1] % 2 == 1), (rings, W, H, G)
        assert g["defects"].tolist() == VO.ref_totals(rings, G=G), (rings, G)


# ------------------------------------------------------------------ K2: validity
def run_k2(n):
    st, bad = collections.Counter(), []
    rings = [VO.k2_ring(i) for i in range(n)]
    for i, ((fam, ring), g) in enumerate(zip(rings, check([([r], None, 16, 16) for _, r in rings]), strict=True)):
        mine, ref = VO.labels(g["defects"]), VO.ring_problems(ring)
        st.update([fam] + [f"label:{lab}" for lab in ref or ["valid"]])
        if mine != ref:
            bad.append((i, ring, mine, ref))
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
    """Exact per-class pair counts (crossing, pinch, overlap, spike, ...) of single rings."""
    rings = [VO.k2_ring(i, seed=31)[1] for i in range(n)]
    for ring, g in zip(rings, check([([r], None, 16, 16) for r in rings]), strict=True):
        assert g["defects"].tolist() == VO.ref_totals([ring]), ring


def run_multi_ring(n):
    rng = random.Random(5)
    shapes = [VO.multi_ring(rng) for _ in range(n)]
    seen = collections.Counter()
    for (rings, roles), g in zip(shapes, check([(r, rl, 64, 64) for r, rl in shapes]), strict=True):
        assert g["defects"].tolist() == VO.ref_totals(rings, roles), (rings, roles)
        seen.update(k for k, v in zip(VO.DEFECTS, g["defects"], strict=True) if v)
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


def test_unrecorded_cases_fall_back_to_the_exact_oracles():
    """R3 (RTX 5090, 2026-11-07): the K1 "holes" family places points with math.cos/sin and round, which differ in
    the last bit between macOS ARM and Linux x86, so 1 of 4,755 and 3 of 2,853 box-generated cases were never recorded
    on the Mac. Such cases take their reference from the exact integer oracles; recorded cases equal those oracles."""
    rings, roles = [[(100, 100), (900, 130), (870, 910), (120, 880)], [(400, 400), (420, 610), (610, 590)]], [-1, 0]
    assert C.case_key(rings, roles, 37, 29) not in C.goldens()
    ref = C.polygon_reference(rings, roles, 37, 29)
    assert ref["defects"].tolist() == VO.ref_totals(rings, roles)
    assert ref["evenodd"] == C.mask_sha(RO.winding_and_parity(rings, 37, 29)[1] % 2 == 1)
    as_numpy = [np.array(r, np.int64) for r in rings]  # object cases arrive as int64 arrays (test_gpu_scoring)
    got = C.polygon_reference(as_numpy, np.array(roles), 37, 29)
    assert got["defects"].tolist() == ref["defects"].tolist() and got["evenodd"] == ref["evenodd"]
    recorded = [c for c in C.fuzz_polygons(20, 0) if c[0] == "holes"][:5]
    for _, rings, _, W, H, G in recorded:
        roles = [-1] * len(rings)
        assert C.case_key(rings, roles, W, H, G) in C.goldens()
        ref = C.polygon_reference(rings, roles, W, H, G)
        assert ref["defects"].tolist() == VO.ref_totals(rings, G=G)
        assert ref["evenodd"] == C.mask_sha(RO.winding_and_parity(rings, W, H, G)[1] % 2 == 1)
