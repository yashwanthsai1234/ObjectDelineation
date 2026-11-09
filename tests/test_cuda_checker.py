"""objdel._gpucheck (cuda/gpucheck.cpp) equals the reference checker bit for bit: the even-odd masks (bit-packed
rows) and defect totals of the archived CPU checker, polycore (recorded per case in tests/goldens/polycore.npz), and
polygonize's rings and bitquad counts. Runs on whichever build is installed:
the CPU build by default, the CUDA build after `make cuda` (tests/test_cuda_checker_gpu.py then also checks it is
one)."""
import random
import re
from pathlib import Path

import cv2
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
                              ("evenodd", lambda g: g["evenodd"].__setitem__((0, 0), True), "evenodd")]:
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
                            ("ring_area", lambda o: o["ring_area"].__setitem__(1, 3), "b0"),  # the hole turns outer
                            ("betti", lambda o: o["betti"].__setitem__((0, 1), 2), "betti"),
                            ("ring_part", lambda o: o["ring_part"].__setitem__(1, 1), "parts")]:
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
    assert (seen > 0).all(), dict(zip(gc.DEFECTS, seen.tolist(), strict=True))


def _strips(k, holes, span=(100, 900)):
    """k strips side by side sharing their long edges: k outers, or k holes of one outer (touching each other, and
    the outer too when they span it)."""
    w = 800 // k
    rects = [[(100 + i * w, span[0]), (100 + (i + 1) * w, span[0]), (100 + (i + 1) * w, span[1]),
              (100 + i * w, span[1])] for i in range(k)]
    if not holes:
        return rects, None
    return [[(100, 50), (100 + k * w, 50), (100 + k * w, 950), (100, 950)]] + [r[::-1] for r in rects], [-1] + [0] * k


def test_validity_adversarial_families():
    """I17 (D6): touching strips (ring-ring touching, nested holes and outers), corner-touching squares, +-2^31
    coordinates, the 10 x 1,000-vertex caps, and a straggler batch (cap shapes among small rollouts): the 13 defect
    totals equal polycore's recorded results or the exact oracles, and each shape of the batch equals its run alone."""
    lo, hi = -2**31, 2**31 - 1
    big = [(lo, lo), (hi, lo), (hi, hi), (lo, hi)]
    cases = [("strips", *_strips(k, False), 64, 64, 1000) for k in (2, 5, 10)]
    cases += [("strip_holes", *_strips(k, True, span), 64, 64, 1000) for k in (2, 5, 9) for span in ((100, 900),
                                                                                                    (50, 950))]
    cases += [("corners", [[(x, y), (x + 100, y), (x + 100, y + 100), (x, y + 100)] for x, y in
                           ((100, 100), (200, 200), (300, 100), (400, 200))], None, 32, 32, 1000),
              ("int32", [big, [(0, 0), (0, 10), (10, 10), (10, 0)]], [-1, 0], 9, 9, 2**30),
              ("int32", [big, big], None, 3, 3, 1)]
    sq = [[(500 - r, 500 - r), (500 + r, 500 - r), (500 + r, 500 + r), (500 - r, 500 + r)] for r in range(450, 0, -45)]
    cases += [("nested_outers", sq, None, 64, 64, 1000),                          # 45 nested pairs
              ("nested_holes", [sq[0]] + [r[::-1] for r in sq[1:]], [-1] + [0] * 9, 64, 64, 1000),  # 36 pairs
              ("hole_outside", [sq[5], sq[0][::-1]], [-1, 0], 64, 64, 1000)]
    assert polygon_mismatches(cases) == []
    caps = [c for c in C.fuzz_polygons(150, seed=1) if c[0] == "caps"]
    assert {c[3] for c in caps} == {64, 512, 4096}
    small = [c for c in C.fuzz_polygons(60, seed=3) if c[0] == "rollout"]
    batch = small[:30] + caps + small[30:]
    p = C.pack_polygons([(r, rl, W, H) for _, r, rl, W, H, _ in batch])
    whole = C.run_polygons(p)
    for s in range(len(batch)):
        q = C.pack_polygons([C.shape_of(p, s)])
        one, g = C.gpu_polygon(C.run_polygons(q), q, 0), C.gpu_polygon(whole, p, s)
        assert np.array_equal(one["defects"], g["defects"]) and np.array_equal(one["evenodd"], g["evenodd"]), s


@pytest.mark.parametrize("rings, roles, W, H", [
    ([], None, 5, 4),                       # no rings: empty mask, legal
    ([[]], None, 3, 3),                     # one empty ring
    ([[(5, 5)]], None, 1, 1),               # one vertex
    ([[(0, 0), (1000, 1000)]], None, 2, 2),  # two vertices
    ([SQ, HOLE], [-1, 0], 1, 1),            # 1 x 1 crop
    ([SQ, HOLE], [-1, 0], 300, 7),          # wider than 256 columns
    ([SQ, HOLE], [-1, 0], 1500, 5),         # wider than one 1024-column window
    ([SQ, HOLE], [-1, 0], 2049, 3),         # three windows, the last one word of one pixel
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
        q = C.pack_polygons([C.shape_of(p, s)])
        one, g = C.gpu_polygon(C.run_polygons(q), q, 0), C.gpu_polygon(whole, p, s)
        assert np.array_equal(one["defects"], g["defects"]) and np.array_equal(one["evenodd"], g["evenodd"])


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
def k1_masks(seed: int) -> list:
    """K1-draft's generators (docs/implementation2.0/kernels/K1/masks_gen.h) as bool masks: noise 0.05-0.95, the
    structured families at sides around one word (31/32/33) and two, 1xN, Nx1, 0xW and Hx0."""
    rng = np.random.default_rng(seed)
    out = [rng.random((s, int(rng.integers(1, 70)))) < p for p in (.05, .1, .3, .5, .7, .9, .95) for s in (31, 32, 33)]
    for h, w in ((31, 33), (32, 32), (33, 31), (64, 65)):
        out += k1_families(rng, h, w)
    return out + [np.ones((1, 77), bool), rng.random((1, 77)) < .5, rng.random((77, 1)) < .5,
                  np.zeros((0, 7), bool), np.zeros((7, 0), bool), np.zeros((0, 0), bool)]


def k1_families(rng, h, w) -> list:
    """Blobs with holes and islands, checkerboard (all saddles), nested squares, comb, spiral, diagonals, vertical
    and horizontal stripes, a frame of 1-pixel holes."""
    y, x = np.indices((h, w))
    d = np.minimum(np.minimum(x, y), np.minimum(w - 1 - x, h - 1 - y))
    frame = np.ones((h, w), bool)
    frame[1:h - 1:2, 1:w - 1:2] = False
    return [blobs(rng, h, w, 3, min(h, w) / 4), (x + y) % 2 == 1, (d // 2) % 2 == 0, (y == h - 1) | (x % 2 == 0),
            spiral(h, w), ((x + y) % 5 == 0) | ((x - y + 1000) % 7 == 0), x % 2 == 0, y % 2 == 0, frame]


def blobs(rng, h, w, k, rmax):
    """k random discs, then k smaller discs cut out: blob-like parts with holes and islands (masks_gen.h blobs)."""
    m = np.zeros((h, w), bool)
    y, x = np.indices((h, w)) + .5
    for t in range(2 * k):
        cx, cy, r = rng.random() * w, rng.random() * h, (.3 + .7 * rng.random()) * rmax * (1 if t < k else .45)
        m[(x - cx) ** 2 + (y - cy) ** 2 <= r * r] = t < k
    return m


def spiral(h, w):
    """One long thin spiral: one ring as long as the mask is large (masks_gen.h spiral)."""
    m = np.zeros((h, w), bool)
    x0, y0, x1, y1 = 0, 0, w - 1, h - 1
    while x0 <= x1 and y0 <= y1:
        m[y0, x0:x1 + 1] = m[y0:y1 + 1, x1] = True
        if y1 > y0:
            m[y1, x0:x1 + 1] = True
        if x0 + 2 <= x1:
            m[y0 + 2:y1 + 1, x0] = True
        if x0 + 2 <= x1 - 2 and y0 + 2 <= y1:
            m[y0 + 2, x0 + 1] = True
        x0, y0, x1, y1 = x0 + 2, y0 + 2, x1 - 2, y1 - 2
    return m


def mixed_masks(seed: int, n: int = 40) -> list:
    """n masks of every kind for the batching tests: the K1 families (empty ones included) and two larger ones."""
    ms = k1_masks(seed)
    return ms[::max(1, len(ms) // (n - 2))][:n - 2] + [spiral(128, 128), np.random.default_rng(seed).random((96, 120))
                                                       < .5]


def masks_out(masks, **kw) -> dict:
    p = C.pack_masks(masks)
    return gc.check_masks(p["bits"], p["hw"], p["start"], **kw)


def one_mask(out, i) -> tuple:
    """Mask i's slices of a check_masks result, as bytes (ring_part and pts relative to the mask)."""
    a, b = out["ring_start"][i:i + 2]
    pa, pb = out["pts_start"][a], out["pts_start"][b]
    return (out["quads"][i].tobytes(), out["betti"][i].tobytes(), out["ring_area"][a:b].tobytes(),
            out["ring_part"][a:b].tobytes(), (out["pts_start"][a:b + 1] - pa).tobytes(), out["pts"][pa:pb].tobytes())


def test_fuzz_masks_match_polygonize_and_bitquads():
    masks = [m for _, m in C.fuzz_masks(600, seed=3)]
    assert mask_mismatches(masks) == []


def test_k1_masks_match_the_oracle_with_parts_and_betti():
    """Rings, quads, ring parts (cv2) and Betti numbers (cv2 8-connected foreground, scipy 4-connected background)."""
    import scipy.ndimage as ndi
    masks = k1_masks(seed=7)
    out = masks_out(masks)
    assert out["ring_part"].dtype == np.int32 and out["betti"].dtype == np.int64
    assert out["betti"].shape == (len(masks), 2)
    assert mask_mismatches(masks) == []
    for i, m in enumerate(masks):
        if m.size == 0:
            assert out["betti"][i].tolist() == [0, 0]
            continue
        b0 = cv2.connectedComponents(m.astype(np.uint8), connectivity=8)[0] - 1
        b1 = ndi.label(~np.pad(m, 1))[1] - 1
        n1, _, nd, n3, _ = out["quads"][i].tolist()
        assert out["betti"][i].tolist() == [b0, b1] and b0 - b1 == (n1 - n3 - 2 * nd) // 4, i
        a, b = out["ring_start"][i:i + 2]
        assert sorted(set(out["ring_part"][a:b].tolist())) == list(range(b0)), i


def test_masks_single_pixel_literals():
    """U01: one ring (x, y) -> (x+1, y) -> (x+1, y+1) -> (x, y+1), area 1, part 0, quads (4, 0, 0, 0, 0),
    betti (1, 0)."""
    masks = [np.ones((1, 1), bool)]
    for k in range(25):
        m = np.zeros((5, 5), bool)
        m[k // 5, k % 5] = True
        masks.append(m)
    out = masks_out(masks)
    for i, m in enumerate(masks):
        y, x = np.argwhere(m)[0]
        a, b = out["ring_start"][i:i + 2]
        assert b - a == 1 and out["ring_area"][a] == 1 and out["ring_part"][a] == 0
        assert out["pts"][out["pts_start"][a]:out["pts_start"][a + 1]].tolist() == [[x, y], [x + 1, y], [x + 1, y + 1],
                                                                                    [x, y + 1]]
        assert out["quads"][i].tolist() == [4, 0, 0, 0, 0] and out["betti"][i].tolist() == [1, 0]


def test_masks_empty_and_degenerate_shapes():
    """U02: empty masks in one batch with others, one between two non-empty ones; neighbours equal their solo call."""
    m = np.zeros((6, 6), bool)
    m[1:5, 1:5], m[2, 2] = True, False
    empty = [np.zeros(s, bool) for s in ((0, 0), (0, 7), (7, 0), (1, 0), (1, 1), (3, 3))]
    masks = [m, *empty[:3], m.T, *empty[3:], m[::-1]]
    out = masks_out(masks)
    assert out["ring_start"][-1] == len(out["ring_area"]) and out["pts_start"][-1] == len(out["pts"])
    assert (np.diff(out["ring_start"]) >= 0).all() and (np.diff(out["pts_start"]) > 0).all()
    for i, mk in enumerate(masks):
        if not mk.any():
            assert out["ring_start"][i] == out["ring_start"][i + 1]
            assert out["quads"][i].tolist() == [0] * 5 and out["betti"][i].tolist() == [0, 0]
        else:
            assert one_mask(out, i) == one_mask(masks_out([mk]), 0)
    assert mask_mismatches(masks) == []
    none = masks_out([])
    assert [none[k].shape for k in ("quads", "ring_start", "pts_start", "ring_area", "pts", "ring_part", "betti")] == \
        [(0, 5), (1,), (1,), (0,), (0, 2), (0,), (0, 2)]


@pytest.mark.parametrize("h, w", [(1, 1), (1, 7), (7, 1), (1, 33), (33, 1), (31, 31), (32, 32), (33, 33), (64, 64),
                                  (65, 3)])
def test_masks_full_and_word_boundaries(h, w):
    """U03: padding bits past W are not foreground."""
    out = masks_out([np.ones((h, w), bool)])
    assert out["pts"].tolist() == [[0, 0], [w, 0], [w, h], [0, h]] and out["ring_area"].tolist() == [h * w]
    assert out["ring_part"].tolist() == [0] and out["betti"].tolist() == [[1, 0]]


def test_masks_connectivity_pair():
    """U04: foreground 8-connected, background 4-connected, the border is background."""
    diag, anti = np.eye(2, dtype=bool), np.eye(2, dtype=bool)[::-1]
    out = masks_out([diag, anti])
    assert out["betti"].tolist() == [[1, 0], [1, 0]] and out["ring_area"].tolist() == [2, 2]
    assert out["pts"][:out["pts_start"][1]].tolist() == [[0, 0], [1, 0], [1, 1], [2, 1], [2, 2], [1, 2], [1, 1], [0, 1]]
    hole = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 0]], bool)
    full7 = np.ones((7, 7), bool)
    full7[3, 3] = False
    corner = np.ones((5, 5), bool)
    corner[0, 0] = corner[1, 1] = False  # touches the border background only through a diagonal: still a hole
    split = np.ones((5, 5), bool)
    split[1, 1] = split[2, 2] = split[3, 3] = False  # three holes touching only diagonally
    out = masks_out([hole, full7, corner, split])
    assert out["betti"].tolist() == [[1, 1], [1, 1], [1, 1], [1, 3]]
    assert out["ring_area"][:4].tolist() == [8, -1, 49, -1] and out["ring_part"][:4].tolist() == [0, 0, 0, 0]
    assert mask_mismatches([diag, anti, hole, full7, corner, split]) == []


def test_masks_nested_parts_and_islands():
    """U05: 9x9 frame, its hole, an island inside. Rings go by their first east edge: the hole's are on its bottom
    side (y = 8), after the island's top (y = 4), so the order is frame, island, hole with parts 0, 1, 0. 21x21 nested
    frames: outers first (parts by first-pixel rank), then the holes from the innermost out."""
    m = np.ones((9, 9), bool)
    m[1:8, 1:8], m[4, 4] = False, True
    out = masks_out([m])
    assert out["ring_area"].tolist() == [81, 1, -49] and out["ring_part"].tolist() == [0, 1, 0]
    assert out["betti"].tolist() == [[2, 1]]
    nested = np.zeros((21, 21), bool)
    for k in range(0, 10, 2):
        nested[k:21 - k, k:21 - k], nested[k + 1:20 - k, k + 1:20 - k] = True, False
    out = masks_out([nested])
    assert out["ring_part"].tolist() == [0, 1, 2, 3, 4, 4, 3, 2, 1, 0] and out["betti"].tolist() == [[5, 5]]
    assert mask_mismatches([m, nested]) == []


def test_masks_touching_borders():
    """U06: the image border reads as background on all four sides."""
    cases = []
    for sl in (np.s_[0, :], np.s_[-1, :], np.s_[:, 0], np.s_[:, -1]):
        m = np.zeros((6, 6), bool)
        m[sl] = True
        cases.append(m)
    c = np.zeros((6, 6), bool)
    c[0, 0] = c[0, -1] = c[-1, 0] = c[-1, -1] = True
    frame = np.ones((6, 6), bool)
    frame[1:5, 1:5] = False
    out = masks_out([*cases, c, frame])
    assert out["betti"].tolist() == [[1, 0]] * 4 + [[4, 0], [1, 1]]
    assert out["ring_area"][-2:].tolist() == [36, -16]
    assert mask_mismatches([*cases, c, frame]) == []


def test_masks_ring_order_literal():
    """U16: rings in raster order of their first east edge, corners from the minimum edge key."""
    m = np.zeros((4, 5), bool)
    m[0, 3] = True              # island in the top row, right
    m[1:4, 0], m[3, 0:3] = True, True  # L shape, first pixel (0, 1)
    out = masks_out([m])
    rings = [out["pts"][out["pts_start"][r]:out["pts_start"][r + 1]].tolist() for r in range(2)]
    assert rings == [[[3, 0], [4, 0], [4, 1], [3, 1]], [[0, 1], [1, 1], [1, 3], [3, 3], [3, 4], [0, 4]]]
    assert out["ring_part"].tolist() == [0, 1]


def test_masks_split_into_chunks_equal_one_call():
    """40 mixed small masks with 4096-pixel chunks and a 64 KiB scratch cap (chunks halve down to one mask) give the
    bytes of the one-chunk call."""
    rng = np.random.default_rng(5)
    masks = [m for h, w in ((9, 11), (12, 12), (14, 13)) for m in k1_families(rng, h, w)] + k1_masks(5)[-6:]
    masks += [rng.random((12, 17)) < p for p in (.05, .2, .4, .5, .6, .8, .95)]
    assert len(masks) == 40
    one, split = masks_out(masks), masks_out(masks, chunk_pixels=4096, max_scratch_bytes=65536)
    assert all(one[k].tobytes() == split[k].tobytes() for k in one)


@pytest.mark.parametrize("chunk", [1, 64, 4096, 65536, 1 << 25])
@pytest.mark.parametrize("scratch", [4096, 65536, 1 << 30])
def test_masks_independent_of_batching(chunk, scratch):
    """P03: each mask alone, in one batch, in any order, among 200 random neighbours, under any chunk and scratch cap
    gives the same bytes; a cap too small for one mask raises ValueError, never wrong data."""
    masks = mixed_masks(3)
    solo = [one_mask(masks_out([m]), 0) for m in masks]
    try:
        out = masks_out(masks, chunk_pixels=chunk, max_scratch_bytes=scratch)
    except ValueError as e:
        assert scratch < 1 << 20 and "max_scratch_bytes" in str(e)
        return
    assert [one_mask(out, i) for i in range(len(masks))] == solo
    rng = np.random.default_rng(chunk % 97)
    perm = rng.permutation(len(masks))
    out = masks_out([masks[j] for j in perm], chunk_pixels=chunk, max_scratch_bytes=scratch)
    assert [one_mask(out, int(np.argwhere(perm == i)[0, 0])) for i in range(len(masks))] == solo
    if chunk == 4096 and scratch == 1 << 30:
        noise = [rng.random(tuple(rng.integers(1, 40, 2))) < rng.random() for _ in range(200)]
        out = masks_out(noise[:100] + masks + noise[100:])
        assert [one_mask(out, 100 + i) for i in range(len(masks))] == solo


def test_masks_seed_rows_do_not_change_results():
    """P04: seed rows 0 (one trace per ring) and 1, 2, 3, 16 give identical output, long rings included."""
    masks = k1_masks(9) + [spiral(128, 128), np.indices((96, 96))[1] % 2 == 0]
    want = masks_out(masks, seed_rows=16)
    for k in (0, 1, 2, 3):
        got = masks_out(masks, seed_rows=k)
        assert all(got[f].tobytes() == want[f].tobytes() for f in want), k
    assert all(want[f].tobytes() == v.tobytes() for f, v in masks_out(masks).items())


def test_masks_adversarial_at_512():
    """P05: a 512 spiral (one long ring), a lattice of 16,384 single pixels (ring k at pixel k), a full mask with a
    lattice of 1-pixel holes."""
    lattice = np.zeros((512, 512), bool)
    lattice[::4, ::4] = True
    holes = np.ones((512, 512), bool)
    holes[1:511:2, 1:511:2] = False
    out = masks_out([spiral(512, 512), lattice, holes])
    assert out["betti"].tolist() == [[1, 0], [16384, 0], [1, 255 * 255]]
    a, b = out["ring_start"][1:3]
    first = out["pts"][out["pts_start"][a:b]]
    assert first.tolist() == [[x, y] for y in range(0, 512, 4) for x in range(0, 512, 4)]
    assert out["ring_part"][a:b].tolist() == list(range(16384))
    assert mask_mismatches([spiral(512, 512), lattice, holes]) == []


def test_masks_stale_scratch_does_not_leak():
    """I07: dense noise fills the arenas, then sparse and empty batches reuse them (nothing is zeroed for free)."""
    rng = np.random.default_rng(6)
    dense = [rng.random((300, 300)) < .5 for _ in range(4)]
    sparse = [rng.random((40, 50)) < .02 for _ in range(30)]
    empty = [np.zeros((30, 30), bool), np.zeros((0, 4), bool)]
    for batches in ([dense, sparse, empty], [empty, sparse, dense], [dense, sparse, empty]):
        for ms in batches:
            assert mask_mismatches(ms) == []
        gc.trim(0)  # the third round runs on memory the pool handed back (CUDA build; no-op on the CPU build)


def test_masks_scratch_too_small_for_one_mask_raises():
    """F04: no truncation, no crash; the next call is right; chunk_pixels=1 runs every mask alone."""
    m = np.random.default_rng(8).random((512, 512)) < .5
    with pytest.raises(ValueError, match="max_scratch_bytes"):
        masks_out([m], max_scratch_bytes=1024)
    assert mask_mismatches([m]) == []
    masks = k1_masks(10)[:20]
    assert all(v.tobytes() == masks_out(masks)[k].tobytes() for k, v in masks_out(masks, chunk_pixels=1).items())


def test_masks_bad_options_raise():
    """F01, F05: the test keywords are checked before any work (a chunk never holds more than 2^30 pixels)."""
    p = C.pack_masks([np.ones((3, 4), bool)])
    for kw in (dict(chunk_pixels=0), dict(chunk_pixels=2**31), dict(max_scratch_bytes=0), dict(seed_rows=-1)):
        with pytest.raises(ValueError):
            gc.check_masks(p["bits"], p["hw"], p["start"], **kw)
    raw = (np.arange(12) * 23 % 256).astype(np.uint8).reshape(3, 4)  # every nonzero byte is foreground
    got = gc.check_masks(raw.ravel(), np.array([raw.shape], np.int32), np.array([0, raw.size], np.int64))
    assert all(got[k].tobytes() == v.tobytes() for k, v in masks_out([raw != 0]).items())


def test_masks_at_the_side_cap():
    """F03 (fast part): 4096 x 1 and 1 x 4096 equal the oracle; a side of 4097 raises a ValueError naming 4096."""
    rng = np.random.default_rng(11)
    assert mask_mismatches([rng.random((4096, 1)) < .5, rng.random((1, 4096)) < .5, np.ones((4096, 2), bool)]) == []
    with pytest.raises(ValueError, match="4096"):
        masks_out([np.ones((4097, 1), bool)])


def test_masks_larger_than_a_chunk_go_alone():
    """Masks larger than one chunk go alone; neighbours in the batch are unaffected."""
    rng = np.random.default_rng(4)
    big = np.zeros((1500, 1500), bool)
    big[100:200, 300:400], big[150, 350], big[1499, 1499] = True, False, True
    masks = [rng.random((5, 6)) < .5, big, rng.random((7, 3)) < .5]
    assert mask_mismatches(masks) == []
    assert all(v.tobytes() == masks_out(masks)[k].tobytes() for k, v in masks_out(masks, chunk_pixels=50).items())


def test_bad_masks_input_raises():
    p = C.pack_masks([np.ones((3, 4), bool)])
    for bits, hw, start in [(p["bits"], np.array([[4, 4]], np.int32), p["start"]),
                            (p["bits"], np.array([[3, 4]], np.int32), np.array([0, 11], np.int64)),
                            (p["bits"], np.array([[-1, -12]], np.int32), p["start"]),
                            (p["bits"], np.array([[3, 4], [0, 0]], np.int32), p["start"]),
                            (p["bits"], np.array([[3, 4]], np.int32), np.array([1, 12], np.int64)),
                            (p["bits"], np.array([[4, 3]], np.int32), np.array([0, 12], np.int64))]:
        if hw.tolist() == [[4, 3]]:  # swapped columns: a valid 4 x 3 mask of the same bytes
            assert gc.check_masks(bits, hw, start)["ring_area"].tolist() == [12]
            continue
        with pytest.raises(ValueError):
            gc.check_masks(bits, hw, start)


def test_no_test_hooks_in_this_build():
    """The suite's build carries no test-only compile definition (tests/variants.py builds those apart)."""
    assert gc.test_hooks() == ()


def test_slot_polygonizer_is_archived():
    """S04 (D2 part): the dense slot polygonizer lives in archive/ only."""
    root = Path(__file__).resolve().parents[1]
    arch = root / "archive" / "cuda_slot_polygonizer"
    assert (arch / "gpucheck_slots.cpp").is_file() and "InitSlots" in (arch / "gpucheck_slots.cpp").read_text()
    readme = (arch / "README.md").read_text()
    assert re.search(r"\b[0-9a-f]{7,40}\b", readme) and "K1" in readme
    src = (root / "cuda" / "gpucheck.cpp").read_text()
    assert not [s for s in ("InitSlots", "LabelRound", "RankRound", "check_mask_chunk", "jump_rounds", "kChunkSlots")
                if s in src]


def test_shape_totals_is_archived():
    """S04 (D6 part): the one-thread-per-shape validity totals live in archive/ only; the engine runs the warp form."""
    root = Path(__file__).resolve().parents[1]
    arch = root / "archive" / "cuda_validity"
    assert "struct ShapeTotals {" in (arch / "shape_totals.cpp").read_text()
    readme = (arch / "README.md").read_text()
    assert re.search(r"\b[0-9a-f]{7,40}\b", readme) and "ShapeTotalsWarp" in readme
    src = (root / "cuda" / "gpucheck.cpp").read_text()
    assert "struct ShapeTotals {" not in src and "ShapeTotals{" not in src
    assert "for_each_warp(st, c.n, ShapeTotalsWarp{" in src


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


@pytest.mark.data
@pytest.mark.slow
def test_lidc_crop_masks_counts_and_oracle():
    """R08: data/checker_eval/masks.npz in one call: 37,093 crops, 38,002 rings, 1,476,508 corners, every crop equal to
    the numpy oracle with its ring parts (cv2) and Betti numbers (CPU build; the box half is the GPU test)."""
    path = Path(__file__).resolve().parents[1] / "data" / "checker_eval" / "masks.npz"
    if not path.is_file():
        pytest.skip(f"{path} not found")
    z = dict(np.load(path, allow_pickle=False))
    out = gc.check_masks(z["bits"], z["hw"], z["start"])
    assert (len(out["quads"]), len(out["ring_area"]), len(out["pts"])) == (37093, 38002, 1476508)
    for k in range(len(z["hw"])):
        m = z["bits"][z["start"][k]:z["start"][k + 1]].reshape(z["hw"][k]).astype(bool)
        assert C.mask_diff(out, k, C.mask_reference(m)) is None, k


def test_cmake_sets_the_cuda_archs_on_the_target():
    """R0 v3 (2026-11-03): the default `native` was set after add_library(gpucheck), and a target copies
    CMAKE_CUDA_ARCHITECTURES only when it is created, so `make cuda` failed on the box. The docker check passed the
    archs on the command line and never saw it; no GPU-less build can, so the file is read."""
    text = (Path(__file__).resolve().parents[1] / "CMakeLists.txt").read_text()
    cuda = text[text.index("if(OBJDEL_CUDA)"):text.index("endif()", text.index("enable_language(CUDA)"))]
    archs = r'CUDA_ARCHITECTURES "\$\{CMAKE_CUDA_ARCHITECTURES\}"'
    assert re.search(r"set_target_properties\(gpucheck PROPERTIES[^)]*" + archs, cuda)
