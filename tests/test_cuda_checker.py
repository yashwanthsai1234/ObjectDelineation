"""objdel._gpucheck (cuda/gpucheck.cpp) equals the reference checker bit for bit: polycore's even-odd and nonzero
masks and validate totals, and polygonize's rings and bitquad counts. Runs on whichever build is installed: the CPU
build by default, the CUDA build after `make cuda` (tests/test_cuda_checker_gpu.py then also checks it is one)."""
import random
import re
from pathlib import Path

import numpy as np
import pytest

from experiments.cuda_checker import compare as C
from objdel import _gpucheck as gc
from objdel.checker import certify

SQ = [(100, 100), (900, 100), (900, 900), (100, 900)]
HOLE = [(400, 400), (400, 600), (600, 600), (600, 400)]


def polygon_mismatches(cases):
    """[(family, rings, roles, W, H, G)] -> [(family, index, first differing field)], one batch per grid G."""
    bad = []
    for G in sorted({c[5] for c in cases}):
        sub = [(i, c) for i, c in enumerate(cases) if c[5] == G]
        p = C.pack_polygons([(r, rl, W, H) for _, (_, r, rl, W, H, _) in sub])
        out = C.run_polygons(p, G)
        for s, (i, (fam, rings, roles, W, H, _)) in enumerate(sub):
            ref = C.polygon_reference(rings, [-1] * len(rings) if roles is None else roles, W, H, G)
            if d := C.polygon_diff(C.gpu_polygon(out, p, s), ref):
                bad.append((fam, i, d))
    return bad


def mask_mismatches(masks):
    p = C.pack_masks(masks)
    out = C.run_masks(p)
    return [(i, d) for i, m in enumerate(masks) if (d := C.mask_diff(out, i, C.mask_reference(np.asarray(m, bool))))]


def test_comparison_reports_the_first_differing_field():
    """The equivalence checks have teeth: a perturbed output is caught, under the right field name."""
    p = C.pack_polygons([([SQ, HOLE], [-1, 0], 16, 16)])
    out = C.run_polygons(p)
    ref = C.polygon_reference([SQ, HOLE], [-1, 0], 16, 16)
    assert C.polygon_diff(C.gpu_polygon(out, p, 0), ref) is None
    for field, edit, want in [("defects", lambda g: g["defects"].__setitem__(0, 1), "legal"),
                              ("evenodd", lambda g: g["evenodd"].__setitem__((0, 0), True), "evenodd"),
                              ("nonzero", lambda g: g["nonzero"].__setitem__((0, 0), True), "nonzero")]:
        g = {k: v.copy() for k, v in C.gpu_polygon(out, p, 0).items()}
        edit(g)
        assert C.polygon_diff(g, ref) == want, field
    ref2 = {**ref, "defects": ref["defects"] + np.eye(13, dtype=np.int64)[2]}
    g = C.gpu_polygon(out, p, 0)
    assert C.polygon_diff(g, ref2) == "legal" and C.polygon_diff({**g, "defects": ref2["defects"]}, ref) == "legal"
    m = np.zeros((6, 6), bool)
    m[1:5, 1:5], m[2:4, 2:4] = True, False
    mo = C.run_masks(C.pack_masks([m]))
    ref = C.mask_reference(m)
    assert C.mask_diff(mo, 0, ref) is None
    for key, edit, want in [("quads", lambda o: o["quads"].__setitem__((0, 0), 9), "quads"),
                            ("pts", lambda o: o["pts"].__setitem__((0, 0), 2), "ring_sha"),
                            ("ring_area", lambda o: o["ring_area"].__setitem__(1, 3), "b0")]:  # the hole turns outer
        o = {k: v.copy() for k, v in mo.items()}
        edit(o)
        assert C.mask_diff(o, 0, ref) == want, key


# ------------------------------------------------------------------ polygons
def test_fuzz_families_match_polycore():
    """K1 raster families, K2 rings, multi-ring shapes with roles, rollouts, int32-extreme coordinates, caps."""
    cases = C.fuzz_polygons(60, seed=1)
    assert {c[0] for c in cases} >= {"star", "soup", "holes", "degenerate", "roles", "rollout", "huge", "caps"}
    assert polygon_mismatches(cases) == []


def test_every_defect_class_is_exercised():
    cases = [c for c in C.fuzz_polygons(60, seed=1) if c[5] == 1000]
    p = C.pack_polygons([(r, rl, W, H) for _, r, rl, W, H, _ in cases])
    seen = (C.run_polygons(p)["defects"] > 0).sum(0)
    assert (seen > 0).all(), dict(zip(C.polycore.DEFECTS, seen.tolist(), strict=True))


@pytest.mark.parametrize("rings, roles, W, H", [
    ([], None, 5, 4),                       # no rings: empty mask, legal
    ([[]], None, 3, 3),                     # one empty ring
    ([[(5, 5)]], None, 1, 1),               # one vertex
    ([[(0, 0), (1000, 1000)]], None, 2, 2),  # two vertices
    ([SQ, HOLE], [-1, 0], 1, 1),            # 1 x 1 crop
    ([SQ, HOLE], [-1, 0], 300, 7),          # wider than one 256-column window
    ([SQ, HOLE[::-1]], [-1, 0], 20, 20),    # hole with the wrong winding
    ([SQ, SQ], None, 16, 16),               # identical outers: overlap, nested, winding 2
    ([[(0, 0), (1000, 0), (1000, 1000), (0, 1000)]], None, 4096, 3),  # whole crop, the widest row
])
def test_edge_cases_match_polycore(rings, roles, W, H):
    assert polygon_mismatches([("edge", rings, roles, W, H, 1000)]) == []


def test_batch_results_do_not_depend_on_batching():
    cases = C.fuzz_polygons(10, seed=2)
    cases = [c for c in cases if c[5] == 1000]
    p = C.pack_polygons([(r, rl, W, H) for _, r, rl, W, H, _ in cases])
    whole = C.run_polygons(p)
    for s in range(0, len(cases), 37):
        one = C.run_polygons(C.pack_polygons([C.shape_of(p, s)]))
        g = C.gpu_polygon(whole, p, s)
        assert np.array_equal(one["defects"][0], g["defects"]) and np.array_equal(one["evenodd"], g["evenodd"].ravel())


def test_bad_polygon_input_raises():
    p = C.pack_polygons([([SQ, HOLE], [-1, 0], 8, 8)])
    run = lambda **kw: gc.check_polygons(*({**p, **kw}[k] for k in ("xy", "ring_start", "shape_start", "roles", "W",
                                                                        "H")), kw.get("G", 1000))
    bad = [dict(ring_start=np.array([1, 4, 8], np.int64)), dict(ring_start=np.array([0, 4, 9], np.int64)),
           dict(ring_start=np.array([0, 5, 4], np.int64)), dict(shape_start=np.array([0, 3], np.int64)),
           dict(roles=np.array([0, 0], np.int32)), dict(roles=np.array([-1], np.int32)),
           dict(W=np.array([0], np.int32)),
           dict(H=np.array([4097], np.int32)), dict(G=0), dict(G=2**30 + 1),
           dict(xy=p["xy"].astype(np.int64)), dict(shape_start=np.array([0, 2, 2], np.int64))]
    for kw in bad:
        with pytest.raises((ValueError, TypeError)):
            run(**kw)
    eleven = C.pack_polygons([([SQ] * 11, None, 8, 8)])
    with pytest.raises(ValueError):
        C.run_polygons(eleven)
    assert run()["defects"].sum() == 0


# ------------------------------------------------------------------ masks
def test_fuzz_masks_match_polygonize_and_bitquads():
    masks = [m for _, m in C.fuzz_masks(600, seed=3)]
    assert mask_mismatches(masks) == []


def test_mask_batch_is_split_into_chunks_without_changing_results():
    """Masks larger than one chunk of scratch go alone; neighbours in the batch are unaffected."""
    rng = np.random.default_rng(4)
    big = np.zeros((1500, 1500), bool)  # 4 x 1501^2 slots: more than one chunk
    big[100:200, 300:400], big[150, 350], big[1499, 1499] = True, False, True
    masks = [rng.random((5, 6)) < .5, big, rng.random((7, 3)) < .5]
    assert mask_mismatches(masks) == []


def test_bad_mask_input_raises():
    p = C.pack_masks([np.ones((3, 4), bool)])
    for bits, hw, start in [(p["bits"], np.array([[4, 4]], np.int32), p["start"]),
                            (p["bits"], p["hw"], np.array([0, 11], np.int64)),
                            (p["bits"], np.array([[-1, -12]], np.int32), p["start"]),
                            (p["bits"], np.array([[3, 4], [0, 0]], np.int32), p["start"])]:
        with pytest.raises(ValueError):
            gc.check_masks(bits, hw, start)


@pytest.mark.data
@pytest.mark.parametrize("name", ["lidc_experts", pytest.param("lidc_full", marks=pytest.mark.slow)])
def test_phase0_lidc_masks(fixtures_dir, name):
    """The Phase 0 LIDC masks: same rings (ring hash), Betti numbers and bit-quad counts as certify/polygonize."""
    z = np.load(fixtures_dir / f"{name}.npz", allow_pickle=False)
    masks = [m.astype(bool) for k in sorted(z.files) for m in z[k]]
    out = C.run_masks(C.pack_masks(masks))
    rng = random.Random(0)
    for i in rng.sample(range(len(masks)), min(len(masks), 2000)) if name == "lidc_experts" else range(len(masks)):
        ref = C.mask_reference(masks[i])
        assert C.mask_diff(out, i, ref) is None, i
        r = certify(masks[i])
        assert (r.b0, r.b1) == (ref["b0"], ref["b1"])


def test_cmake_sets_the_cuda_archs_on_the_target():
    """R0 v3 (2026-11-03): the default `native` was set after add_library(gpucheck), and a target copies
    CMAKE_CUDA_ARCHITECTURES only when it is created, so `make cuda` failed on the box. The docker check passed the
    archs on the command line and never saw it; no GPU-less build can, so the file is read."""
    text = (Path(__file__).resolve().parents[1] / "CMakeLists.txt").read_text()
    cuda = text[text.index("if(OBJDEL_CUDA)"):text.index("endif()", text.index("enable_language(CUDA)"))]
    archs = r'CUDA_ARCHITECTURES "\$\{CMAKE_CUDA_ARCHITECTURES\}"'
    assert re.search(r"set_target_properties\(gpucheck PROPERTIES[^)]*" + archs, cuda)
