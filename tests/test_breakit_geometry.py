"""Break-it B1: geometry wiring (parser, checker, topology, targets) under adversarial input.

Tests named `b1_<n>` reproduce confirmed breakdowns, now fixed (experiments/phase1/impl/breakit/B1_geometry.md);
the others are attacks that held. All are kept as regression tests.
"""
import random
import re
import threading
from fractions import Fraction

import numpy as np
import pytest
from oracles import raster_oracle as RO
from test_geometry_oracle import lab_parse

from objdel import _gpucheck
from objdel.checker.bitquad import bitquads
from objdel.checker.topology import instance_topology, min_area, signature, size_filter
from objdel.data import Item
from objdel.geometry import parse, pattern, serialize, shape, to_mask, to_masks
from objdel.reward.cache import ExpertCache
from objdel.reward.score import score_completion
from objdel.targets import build_target, frame

DOC = '{"nodule":[{"outer":[[0,0],[1000,0],[1000,1000]],"holes":[]}]}'
NODULE = ("nodule",)


def _square_with_hole() -> np.ndarray:
    m = np.zeros((50, 50), bool)
    m[10:40, 10:40] = True
    m[20:30, 20:30] = False
    return m


# ------------------------------------------------------------------ confirmed breakdowns (fixed, now regression tests)

@pytest.mark.parametrize("ws", ["　", " ", "\xa0", "\x85", "\x1c", "\x1f", "\x0b"])
def test_b1_1_only_json_whitespace_is_stripped(ws):
    assert parse(DOC, NODULE).ok
    p = parse(ws + DOC + ws, NODULE)
    assert not p.ok  # RFC 8259 whitespace is space, tab, LF, CR only; the rest is not_ascii / syntax


@pytest.mark.parametrize("ring", [np.array([[0, 0, 1000], [0, 1000, 1000]]),
                                  np.array([[0], [0], [1000], [0], [1000], [1000]]),
                                  np.zeros((2, 2, 2), np.int64)])
def test_b1_2_pack_rejects_rings_that_are_not_k_by_2(ring):
    with pytest.raises((ValueError, TypeError)):
        shape.pack_rings([ring])


def test_b1_3_default_min_area_target_keeps_signature_under_a_larger_union_rule():
    S = 500
    yy, xx = np.mgrid[:S, :S]
    m = (yy - 250) ** 2 + (xx - 250) ** 2 < 150 ** 2
    m &= ~((yy - 250) ** 2 / 14 ** 2 + (xx - 250) ** 2 / 8.4 ** 2 < 1)  # one elliptic hole of 371 px
    own = min_area(int(m.sum()))
    t = build_target({"nodule": m}, NODULE)  # default: min_area of this mask (352)
    g = to_mask(t.objs["nodule"], S, S).mask
    union_rule = 371  # a union of raters ~6 % larger than this rater's mask
    assert union_rule >= own
    assert signature(g, own) == signature(m, own)  # certified at the default
    assert signature(g, union_rule) == signature(m, union_rule)  # "never larger, hence stricter"


def test_b1_4_sub_filter_pinholes_do_not_make_a_mask_uncertifiable():
    yy, xx = np.ogrid[:100, :100]
    d = (yy - 50) ** 2 + (xx - 50) ** 2 < 40 ** 2
    d[30:70:4, 50] = False  # ten 1-px pinholes: 11 raw rings, filtered signature (0,)
    assert signature(d, 5) == (0,)
    # a certifiable answer exists: the filtered disc is legal, IoU 0.998 and keeps the filtered signature
    g = to_mask(build_target({"nodule": size_filter(d, 5)}, NODULE).objs["nodule"], 100, 100)
    assert g.legal and np.count_nonzero(g.mask & d) / np.count_nonzero(g.mask | d) >= 0.99
    assert signature(g.mask, 5) == signature(d, 5)
    build_target({"nodule": d}, NODULE)  # raises "11 rings > MAX_RINGS"


@pytest.mark.parametrize("byte", [2, 255])
def test_b1_5_targets_normalise_bool_bytes(byte):
    m = _square_with_hole()
    nb = (m.view(np.uint8) * np.uint8(byte)).view(bool)
    assert np.array_equal(nb, m)
    assert build_target({"nodule": nb}, NODULE).text == build_target({"nodule": m}, NODULE).text


@pytest.mark.parametrize("byte", [2, 255])
def test_b1_5_expert_cache_normalises_bool_bytes(byte):
    m = _square_with_hole()
    nb = (m.view(np.uint8) * np.uint8(byte)).view(bool)
    assert ExpertCache.build(np.stack([nb, nb])).rater_sigs == ExpertCache.build(np.stack([m, m])).rater_sigs


@pytest.mark.parametrize("channels", [2, 4])
def test_b1_6_frame_image_is_rgb_or_refused(channels):
    m = np.zeros((1, 80, 80), bool)
    m[0, 30:50, 30:50] = True
    item = Item("x", "riga", "g", "train", np.full((80, 80, channels), 200, np.uint8), {"disc": m, "cup": m}, ["a"], {})
    try:
        f = frame(item)
    except ValueError:
        return
    assert f.image.shape == (f.H, f.W, 3)


# ------------------------------------------------------------------ attacks that held (regression tests)

def test_parse_accepts_exactly_the_grammar_on_mutants():
    rng = random.Random(1)
    alphabet = list('[]{},:"0123456789outerholesdiscupn') + [" ", "-", ".", "e", "١", "０", "\x00", "\\", "1000", "00"]
    for _ in range(3000):
        names = rng.choice([("disc", "cup"), NODULE])
        objs = {n: [(np.array([[rng.choice([0, 1, 999, 1000, rng.randint(0, 1000)]) for _ in range(2)]
                               for _ in range(rng.choice([2, 3, 4, 7]))]), []) for _ in range(rng.choice([0, 1, 2]))]
                for n in names}
        s = serialize(objs, names)
        for _ in range(rng.choice([0, 1, 2, 3])):
            pos = rng.randint(0, len(s))
            op = rng.random()
            s = (s[:pos] + rng.choice(alphabet) + s[pos:] if op < 0.4 else s[:pos] + s[pos + 1:] if op < 0.8
                 else s[:pos] + rng.choice(alphabet) + s[pos + 1:])
        if any(c.isspace() and not c.isascii() for c in s):
            continue  # B1-1
        p, r = parse(s, names), lab_parse(s, names)  # independent recursive-descent reference
        assert p.ok == (r is not None), s
        if p.ok:
            assert serialize(p.objs, names) == serialize(r, names)


@pytest.mark.parametrize("text", [
    "[" * 50_000 + "]" * 50_000,
    '{"nodule":' * 20_000,
    '{"a":' * 20_000 + "1" + "}" * 20_000,
    '{"nodule":[{"outer":[' + ",".join(["[1000,1000]"] * 20_000) + ',[10000,1]],"holes":[]}]}',
    '{"nodule":[{"outer":[' + ",".join(["[1.5,1]"] * 20_000) + '],"holes":[]}]}',
    '{"nodule":[{"outer":[' + ",".join(["[1,1]"] * 20_000) + ',[01,1]],"holes":[]}]}',
    DOC[:-1] + "\ud800}",
    "﻿" + DOC,
])
def test_parse_is_total_on_hostile_texts(text):
    p = parse(text, NODULE)
    assert not p.ok and p.error
    assert re.fullmatch(pattern(NODULE), text) is None


def test_rasterize_exact_on_boundary_biased_shapes():
    rng = random.Random(7)
    for _ in range(150):
        W, H = rng.choice([1, 2, 3, 7, 13]), rng.choice([1, 2, 3, 7, 13])

        def coord(D):
            k = rng.random()
            if k < 0.3:
                return rng.choice([0, 1000])
            if k < 0.6:
                v = Fraction(rng.randint(0, 2 * D) * 1000, 2 * D)
                if v.denominator == 1:
                    return int(v)
            return rng.randint(-2, 1002)

        rings = [[(coord(W), coord(H)) for _ in range(rng.randint(3, 6))] for _ in range(rng.randint(1, 3))]
        g, = to_masks([([(np.array(x), []) for x in rings], W, H)])  # every ring an outer
        _, c = RO.winding_and_parity(rings, W, H)
        assert np.array_equal(g.mask, c % 2 == 1), (W, H, rings)


def test_checker_concurrent_python_threads_are_deterministic():
    rng = np.random.default_rng(0)
    objs = [([(rng.integers(-50, 1051, (int(rng.integers(3, 80)), 2)), []) for _ in range(int(rng.integers(1, 11)))],
             int(rng.integers(1, 120)), int(rng.integers(1, 120))) for _ in range(60)]
    ref = to_masks(objs)
    bad = []

    def work(k):
        for i in np.random.default_rng(k).permutation(len(objs)):
            g, = to_masks([objs[i]])
            if not (np.array_equal(g.mask, ref[i].mask) and g.defects == ref[i].defects):
                bad.append(i)
        b = _gpucheck.check_polygons(*shape.pack_objects(objs), 1000)
        if [{d: int(c) for d, c in zip(_gpucheck.DEFECTS, row, strict=True) if c} for row in b["defects"]] != \
                [r.defects for r in ref]:
            bad.append(-1)

    threads = [threading.Thread(target=work, args=(k,)) for k in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert not bad


def _sq(x0, y0, x1, y1, hole=False):
    r = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
    return r[::-1] if hole else r


@pytest.mark.parametrize("parts, legal", [
    ([(_sq(0, 0, 1000, 1000), [_sq(100, 100, 900, 900, True)]), (_sq(200, 200, 800, 800), [])], True),  # island
    ([(_sq(0, 0, 1000, 1000), [_sq(100, 100, 900, 900, True)]), (_sq(100, 200, 800, 800), [])], False),
    ([(_sq(0, 0, 1000, 1000), [_sq(100, 100, 900, 900, True), _sq(200, 200, 800, 800, True)])], False),
    ([(_sq(0, 0, 500, 500), []), (_sq(500, 500, 1000, 1000), [])], False),  # parts sharing a corner
    ([(_sq(0, 0, 1000, 1000), [np.array([[500, 0], [400, 300], [600, 300]])])], False),  # hole vertex on outer
    ([(np.array([[0, 0], [0, 0], [1000, 0], [1000, 1000], [0, 1000], [0, 0]]), [])], True),  # repeats
    ([(np.array([[0, 0], [500, 0], [1000, 0], [1000, 1000], [0, 1000]]), [])], True),  # straight-through
    ([(np.array([[5, 5]] * 4), [])], False),
])
def test_validity_corner_cases(parts, legal):
    assert to_mask(parts, 10, 10).legal == legal


@pytest.mark.parametrize("dtype", [np.uint8, np.int8, np.uint16, np.int64, np.uint64, np.float32, np.float64, object])
def test_topology_value_dtypes(dtype):
    m = _square_with_hole()
    assert signature(m.astype(dtype)) == signature(m) == (1,)
    assert signature(np.asfortranarray(m.astype(dtype))[::-1, ::-1]) == (1,)


def test_topology_chi_certificate_on_checkerboards_and_diagonals():
    chk = np.indices((40, 40)).sum(0) % 2 == 0
    t = instance_topology(chk)
    assert bitquads(chk).chi8 == int((1 - t.holes).sum())
    d = np.zeros((7, 7), bool)
    d[[1, 2, 3, 4, 5, 4, 3, 2], [3, 2, 1, 2, 3, 4, 5, 4]] = True  # diagonal-only diamond
    assert signature(d) == (1,) and instance_topology(d).hole_areas.tolist() == [5]


@pytest.mark.parametrize("a", [5, 6, 17])
def test_size_rule_exact_boundaries(a):
    m = np.zeros((40, 40), bool)
    m[5:35, 5:35] = True
    m[10, 10:10 + a] = False  # a hole of exactly a px stays
    assert signature(m, a) == (1,)
    m[10, 10 + a - 1] = True  # a - 1 px is filled
    assert signature(m, a) == (0,)
    s = np.zeros((40, 40), bool)
    s[2:30, 2:30] = True
    s[38, :a] = True
    assert signature(s, a) == (0, 0)
    s[38, a - 1] = False
    assert signature(s, a) == (0,)


def test_targets_recertify_through_the_reward():
    """mask -> target -> parse -> to_mask -> reward: legal, own topology, IoU_tc >= 0.99 (union min_area)."""
    rng = np.random.default_rng(5)
    for S in (50, 125, 250, 500):
        yy, xx = np.mgrid[:S, :S]
        for _ in range(3):
            c, r = rng.uniform(0.3, 0.7, 2) * S, rng.uniform(0.15, 0.3) * S
            a = ((yy - c[0]) ** 2 + (xx - c[1]) ** 2 < r ** 2) & ~((yy - c[0]) ** 2 + (xx - c[1]) ** 2 < (r / 3) ** 2)
            b = a.copy()
            b[int(c[0]), :] = False  # second rater cuts the ring: different topology
            stack = np.stack([a, b])
            cache = ExpertCache.build(stack)
            ma = {"nodule": min_area(int(stack.any(0).sum()))}
            assert ma["nodule"] == cache.min_area
            for k in range(2):
                t = build_target({"nodule": stack[k]}, NODULE, None, ma)
                s = score_completion(t.text, False, NODULE, {"nodule": cache}, S, S)
                assert s["format_ok"] and s["legal"] and s["topo_any"] and s["iou_tc"] >= 0.99


@pytest.mark.parametrize("shape", [(1, 1), (1, 500), (500, 1), (7, 333), (499, 13)])
def test_targets_on_degenerate_and_rectangular_crops(shape):
    H, W = shape
    for m in (np.ones(shape, bool), np.eye(H, W, dtype=bool), np.zeros(shape, bool)):
        t = build_target({"nodule": m}, NODULE)
        g = to_mask(parse(t.text, NODULE).objs["nodule"], W, H)
        assert g.legal and signature(g.mask, min_area(int(m.sum()))) == signature(m, min_area(int(m.sum())))
