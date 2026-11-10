"""SFT targets (plan D11, gates K5-K7, R1 target part): every target parses, is canonical and legal, keeps the
rater mask's topology (raw and size-filtered signature), has IoU >= 0.99 with it, and is deterministic."""
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import cv2
import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from scipy import ndimage

from objdel.geom.shape import canonicalize, to_mask
from objdel.geom.text import MAX_RINGS, parse, serialize
from objdel.geom.topology import min_area, signature, size_filter
from objdel.targets import BuildConfig, build_target, min_areas, target_set_hash
from tests.oracles.crops_oracle import mask_sha, pad_div
from tests.oracles.kept_oracle import kept
from tests.oracles.polygonize_oracle import polygonize


def check(masks, names, ma=None, cfg=None):
    """Build, then certify the text independently of the builder's own bookkeeping: the signature is compared under
    the size rule ma[name] (default: min_areas of the masks themselves, a one-rater union). Masks with more rings than
    the grammar allows (MAX_RINGS) after the size rule must be refused; masks with more raw rings are traced after
    the size rule, which may cost the IoU."""
    cfg = cfg or BuildConfig()
    ma = ma or min_areas({n: np.asarray(masks[n], bool)[None] for n in names})
    H, W = next(iter(masks.values())).shape
    if any(len(polygonize(size_filter(masks[n], ma[n]))) > MAX_RINGS for n in names):
        with pytest.raises(ValueError, match="MAX_RINGS"):
            build_target(masks, names, cfg, ma)
        return None, None
    if any(len(polygonize(masks[n])) > MAX_RINGS for n in names):
        try:
            t = build_target(masks, names, cfg, ma)
        except ValueError as e:
            assert "no certified target" in str(e)
            return None, None
    else:
        t = build_target(masks, names, cfg, ma)
    p = parse(t.text, names)
    assert p.ok, p.error
    assert serialize(canonicalize(p.objs), names) == t.text  # canonical
    assert serialize(t.objs, names) == t.text
    out = {}
    for n in names:
        g = to_mask(p.objs[n], W, H)
        assert g.legal, (n, g.defects)
        src = np.asarray(masks[n], bool)
        assert kept(g.mask, src, ma[n], ma[n]), n
        u = np.count_nonzero(g.mask | src)
        iou = 1.0 if u == 0 else np.count_nonzero(g.mask & src) / u
        assert iou >= cfg.iou_min, (n, iou)
        assert t.info[n]["iou"] == iou
        out[n] = g.mask
    if set(names) == {"disc", "cup"}:
        rim, rim_src = out["disc"] & ~out["cup"], masks["disc"] & ~masks["cup"]
        assert kept(rim, rim_src, ma["rim"], ma["rim"])
    return t, out


def max_component(m):
    n, lab, stats, _ = cv2.connectedComponentsWithStats(m.view(np.uint8), connectivity=8)
    return lab == 1 + np.argmax(stats[1:, cv2.CC_STAT_AREA])


def place(m, S, at=(0, 0)):
    big = np.zeros((S, S), bool)
    big[at[0]:at[0] + m.shape[0], at[1]:at[1] + m.shape[1]] = m
    return big


def disc(S, cy, cx, r):
    y, x = np.mgrid[:S, :S]
    return (y - cy) ** 2 + (x - cx) ** 2 <= r * r


CASES = {
    "one_pixel": np.ones((1, 1), bool),
    "hline": np.ones((1, 9), bool),
    "vline": np.ones((9, 1), bool),
    "diag_line": np.eye(7, dtype=bool),
    "pinch": np.array([[1, 0], [0, 1]], bool),
    "anti_pinch": np.array([[0, 1], [1, 0]], bool),
    "checker": np.add.outer(np.arange(5), np.arange(5)) % 2 == 0,
    "holes_touching_diagonally": np.array([[1, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 1]], bool),
    "hole_touching_outside": np.array([[1, 1, 1, 0], [1, 0, 1, 1], [1, 1, 1, 1]], bool),
    "staircase": np.tril(np.ones((7, 7), bool)),
}
CASES["island_in_hole"] = np.ones((5, 5), bool)
CASES["island_in_hole"][1:4, 1:4] = False
CASES["island_in_hole"][2, 2] = True


@pytest.mark.parametrize("name", sorted(CASES))
@pytest.mark.parametrize("S", [0, 50, 500])  # 0 = the mask's own frame
def test_edge_cases(name, S):
    m = CASES[name]
    if S:
        m = place(m, S, (S // 3, S // 2))
    check({"nodule": m}, ("nodule",))


def touching_ellipses(S, a1, b1, th1, a2, b2, th2):
    """Two ellipses moved so that their extreme pixels along (1, 1) and (-1, -1) touch at one corner (a pinch)."""
    y, x = np.mgrid[:S, :S]
    c = S // 2

    def ellipse(cy, cx, a, b, th):
        u, v = (x - cx) * np.cos(th) + (y - cy) * np.sin(th), (y - cy) * np.cos(th) - (x - cx) * np.sin(th)
        return (u / a) ** 2 + (v / b) ** 2 <= 1

    def moved(m, sign, to):
        ys, xs = np.nonzero(m)
        k = np.argmax(sign * (ys + xs))
        return np.roll(m, (to - ys[k], to - xs[k]), (0, 1))

    return moved(ellipse(c - 31, c - 31, a1, b1, th1), 1, c - 1) | moved(ellipse(c + 30, c + 30, a2, b2, th2), -1, c)


@pytest.mark.parametrize("params", [(12, 55, 2.0, 25, 11, 2.9), (12, 55, 2.0, 26, 9, 1.35), (40, 8, 1.0, 55, 13, 2.2)])
def test_pinch_between_slanted_edges(params):
    """At a pinch the simplified edges are slanted. The bevel must step along the pinch's own crack edges: the lab
    stepped one grid unit along the simplified edges, which at S = 500 can land on a pixel centre and flip it, so
    every rung failed and the object fell back to its exact rings (7 RIGA cups and discs in the store)."""
    m = touching_ellipses(500, *params)
    assert signature(m) == (0,)
    check({"nodule": m}, ("nodule",))
    t, out = check({"nodule": m}, ("nodule",), cfg=BuildConfig(tau=0.0, step=1.0, iou_min=1.0))  # lossless only
    assert t.info["nodule"]["tau"] == 0.0 and np.array_equal(out["nodule"], m)


def test_size_rule_not_raw_topology():
    """A jagged boundary (pen-stroke jitter, like RIGA BinRushed masks) whose tau = 3 polygon rasterizes with two
    sub-filter pinholes: certified under the size rule (plan §4.2) at tau = 3. Demanding the raw signature pushed it
    to tau = 1.25 with 350 vertices (and real BinRushed masks to 550 vertices, 6,000-token answers). The size rule's
    latitude applies to objects outside BuildConfig.exact (RIGA); LIDC nodules must be pixel-exact (BK1-1)."""
    S, rng = 500, np.random.default_rng(16)
    amp, k, jitter = rng.uniform(0.5, 4), int(rng.integers(10, 60)), rng.uniform(1, 4)
    y, x = np.mgrid[:S, :S]
    noise = cv2.GaussianBlur(rng.normal(size=(S, S)), (0, 0), 1.0) / 0.28
    m = np.hypot(y - 250, x - 250) <= 150 + amp * np.sin(k * np.arctan2(y - 250, x - 250)) + jitter * noise
    m = ndimage.binary_fill_holes(max_component(m))
    assert signature(m) == (0,)
    t, out = check({"nodule": m}, ("nodule",), {"nodule": 400}, BuildConfig(exact=()))
    assert signature(out["nodule"]) != (0,)  # the raw raster has sub-filter pinholes
    assert t.info["nodule"]["tau"] == 3.0 and t.info["nodule"]["vertices"] < 150


def test_rings_beyond_the_vertex_cap_are_not_targets():
    """A ring that needs > MAX_VERTICES vertices (here even losslessly: 2 per tooth) is not a target."""
    m = np.zeros((500, 500), bool)
    m[100:400, 100:400] = True
    m[100:400:2, 99] = m[99, 100:400:2] = m[100:400:2, 400] = m[400, 100:400:2] = True  # 600 one-pixel teeth
    with pytest.raises(ValueError, match="no certified target"):
        build_target({"nodule": m}, ("nodule",), BuildConfig(tau=0.0, step=1.0, iou_min=1.0), {"nodule": 5})


def test_island_in_hole_is_two_parts():
    m = CASES["island_in_hole"]
    assert signature(m) == (1, 0)
    t, _ = check({"nodule": m}, ("nodule",))
    assert len(t.objs["nodule"]) == 2 and sum(len(h) for _, h in t.objs["nodule"]) == 1


@pytest.mark.parametrize("S", [1, 7, 50, 500])
def test_empty_mask_is_the_empty_object(S):
    t, out = check({"nodule": np.zeros((S, S), bool)}, ("nodule",))
    assert t.text == '{"nodule":[]}' and not out["nodule"].any()


@pytest.mark.parametrize("S", [40, 125, 500])
def test_object_touching_every_border(S):
    m = np.ones((S, S), bool)
    m[S // 2, S // 2] = False  # one pinhole
    t, out = check({"nodule": m}, ("nodule",))
    assert np.array_equal(out["nodule"], m)
    c0, c1 = S // 2 * 1000 // S, (S // 2 + 1) * 1000 // S  # the pinhole's corners on the grid
    hole = [[c0, c0], [c0, c1], [c1, c1], [c1, c0]]  # counter-clockwise on screen, top-left start
    want = {"nodule": [{"outer": [[0, 0], [1000, 0], [1000, 1000], [0, 1000]], "holes": [hole]}]}
    assert t.text == json.dumps(want, separators=(",", ":"))


def test_rectangle_vertices_are_exact_pixel_corners():
    m = np.zeros((50, 50), bool)
    m[10:31, 5:45] = True  # columns 5..44, rows 10..30
    t, _ = check({"nodule": m}, ("nodule",))
    assert t.text == '{"nodule":[{"outer":[[100,200],[900,200],[900,620],[100,620]],"holes":[]}]}'


def test_large_smooth_objects_riga_like():
    S = 500
    d, c = disc(S, 250, 250, 240), disc(S, 240, 255, 120)  # disc almost fills the crop, cup inside
    t, _ = check({"disc": d, "cup": c}, ("disc", "cup"), {"disc": 300, "cup": 300, "rim": 300})
    assert list(t.objs) == ["disc", "cup"]
    assert t.info["disc"]["vertices"] < 60 and t.info["cup"]["vertices"] < 60  # simplified, not exact rings


def test_rim_topology_kept_when_cup_touches_disc():
    S = 250
    d, c = disc(S, 125, 125, 100), disc(S, 125, 160, 65)  # cup touches the disc border on the right
    c &= d
    check({"disc": d, "cup": c}, ("disc", "cup"))


def test_wiggly_boundary_needs_a_finer_tolerance():
    rng = np.random.default_rng(0)
    S = 500
    y, x = np.mgrid[:S, :S]
    a = np.arctan2(y - 250, x - 250)
    m = np.hypot(y - 250, x - 250) <= 120 + 6 * np.sin(23 * a) + rng.normal(0, 0.3, (S, S))
    t, _ = check({"nodule": m}, ("nodule",))
    assert t.info["nodule"]["tau"] < 3


@st.composite
def blobs(draw, max_side=18):
    H, W = draw(st.integers(1, max_side)), draw(st.integers(1, max_side))
    m = np.zeros((H, W), bool)
    for _ in range(draw(st.integers(1, 6))):
        kind = draw(st.sampled_from(["rect", "hole", "pixel", "diag", "checker", "hline", "vline"]))
        r0, c0 = draw(st.integers(0, H - 1)), draw(st.integers(0, W - 1))
        r1, c1 = draw(st.integers(r0, H - 1)), draw(st.integers(c0, W - 1))
        if kind == "rect":
            m[r0:r1 + 1, c0:c1 + 1] = True
        elif kind == "hole":
            m[r0:r1 + 1, c0:c1 + 1] = False
        elif kind == "pixel":
            m[r0, c0] = ~m[r0, c0]
        elif kind == "diag":
            for k in range(min(r1 - r0, c1 - c0) + 1):
                m[r0 + k, c0 + k] = True
        elif kind == "checker":
            m[r0:r1 + 1, c0:c1 + 1] |= np.add.outer(np.arange(r1 - r0 + 1), np.arange(c1 - c0 + 1)) % 2 == 0
        elif kind == "hline":
            m[r0, c0:c1 + 1] = True
        else:
            m[r0:r1 + 1, c0] = True
    return m


@settings(max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(m=blobs(), S=st.sampled_from([0, 50, 125, 333, 500]), a=st.integers(0, 8))
def test_random_masks(m, S, a):
    if S:
        m = place(m, S, (S // 3, S // 2))
    check({"nodule": m}, ("nodule",), {"nodule": a})


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(seed=st.integers(0, 2**32 - 1), S=st.sampled_from([100, 250, 500]))
def test_random_smooth_blobs(seed, S):
    """Thresholded smooth noise: components, holes and thin necks at realistic scales."""
    rng = np.random.default_rng(seed)
    z = lambda sigma: (lambda a: a / a.std())(cv2.GaussianBlur(rng.normal(size=(S, S)), (0, 0), S * sigma))
    m = (z(0.12) > rng.uniform(0, 1)) & ~(z(0.04) > 2.2)  # blobs (median 3 rings), some holes
    check({"nodule": m}, ("nodule",), {"nodule": min_area(int(m.sum()))})


def _brute_kept(a, b, lo, hi):
    return all(signature(a, t) == signature(b, t) for t in range(lo, hi + 1))


@settings(max_examples=300, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(a=blobs(12), b=blobs(12), lo=st.integers(0, 12), span=st.integers(0, 30))
def test_kept_over_a_range_of_rules_equals_brute_force(a, b, lo, span):
    """The interval loop (tests/oracles/kept_oracle.kept, build.py's _kept before D4) checks the signature only at the
    rules where one can change; it must agree with checking every integer rule in [lo, hi] (pairs of masks of one
    shape, including near-copies)."""
    H, W = min(a.shape[0], b.shape[0]), min(a.shape[1], b.shape[1])
    a, b = a[:H, :W], b[:H, :W]
    for x, y in ((a, b), (a, a ^ (b & np.eye(H, W, dtype=bool))), (a, size_filter(a, max(lo, 1)))):
        assert kept(x, y, lo, lo + span) == _brute_kept(x, y, lo, lo + span)


def test_kept_sees_a_hole_freed_by_a_removed_island():
    """The signature need not be monotone in the rule: a 28-px moat around a 36-px island is filled at rule 30, but
    at rule 40 the island is removed first and the hole (now 64 px) survives. The rules where a signature can change
    include the raw components' sizes + 1, not only the filtered holes' sizes + 1."""
    m = np.ones((12, 12), bool)
    m[2:10, 2:10] = False
    m[3:9, 3:9] = True
    assert (signature(m, 28), signature(m, 30), signature(m, 40), signature(m, 70)) == ((1, 0), (0,), (1,), (0,))
    g = m.copy()
    g[3, 3] = False  # island 35 px, moat 29 px: differs from m at rules 29 and 36 only
    assert kept(m, g, 30, 35) and _brute_kept(m, g, 30, 35)
    assert not kept(m, g, 30, 36) and not _brute_kept(m, g, 30, 36)
    assert not kept(m, g, 29, 29)


def test_ring_cap_applies_after_the_size_rule():
    """B1-4: sub-filter pinholes and specks do not count toward MAX_RINGS; the target is certified against the raw
    mask (size-filtered signature, IoU)."""
    yy, xx = np.ogrid[:100, :100]
    d = (yy - 50) ** 2 + (xx - 50) ** 2 < 40 ** 2
    d[30:70:4, 50] = False  # ten 1-px pinholes
    d[2, 2:96:8] = True  # twelve 1-px specks
    assert len(polygonize(d)) == 23
    t, _ = check({"nodule": d}, ("nodule",), {"nodule": 5})
    assert len(t.objs["nodule"]) == 1 and t.objs["nodule"][0][1] == []
    t, _ = check({"nodule": d}, ("nodule",))
    m = np.zeros((30, 60), bool)
    m[1:29, 1:59] = True
    for i in range(11):
        m[5:7, 3 + 5 * i:6 + 5 * i] = False  # eleven 6-px holes survive the size rule at 5
    with pytest.raises(ValueError, match="MAX_RINGS"):
        build_target({"nodule": m}, ("nodule",), None, {"nodule": 5})


@pytest.mark.parametrize("byte", [2, 255])
def test_bool_bytes_other_than_one(byte):
    """B1-5: a bool mask whose True bytes are not 1 gives the canonical mask's target."""
    m = disc(60, 30, 30, 20) & ~disc(60, 30, 30, 6)
    nb = (m.view(np.uint8) * np.uint8(byte)).view(bool)
    for ma in ({"nodule": 0}, {"nodule": 5}):
        assert build_target({"nodule": nb}, ("nodule",), None, ma).text == build_target({"nodule": m}, ("nodule",),
                                                                                      None, ma).text


def test_deterministic_and_input_untouched():
    m = disc(250, 120, 130, 90) & ~disc(250, 120, 130, 20)
    keep = m.copy()
    a = build_target({"nodule": m}, ("nodule",), None, {"nodule": 5}).text
    assert np.array_equal(m, keep)
    assert build_target({"nodule": m.copy()}, ("nodule",), None, {"nodule": 5}).text == a


def test_deterministic_across_processes():
    code = ("import hashlib, numpy as np; from objdel.targets import build_target;"
            "y, x = np.mgrid[:200, :200]; a = np.arctan2(y - 100, x - 100);"
            "m = np.hypot(y - 100, x - 100) < 80 + 5 * np.sin(7 * a); m[90:99, 60:64] = m[120, 130:140] = False;"
            "t = build_target({'nodule': m}, ('nodule',), None, {'nodule': 5});"
            "print(hashlib.sha256(t.text.encode()).hexdigest())")
    outs = {subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                           env={"PYTHONHASHSEED": str(s)}).stdout for s in (1, 2)}
    assert len(outs) == 1 and len(next(iter(outs)).strip()) == 64


def test_invalid_inputs():
    with pytest.raises(ValueError):
        build_target({"nodule": np.ones((501, 10), bool)}, ("nodule",), None, {"nodule": 5})  # beyond the 500 px
    with pytest.raises(ValueError):
        build_target({"disc": np.ones((5, 5), bool), "cup": np.ones((5, 6), bool)}, ("disc", "cup"), None,
                     {"disc": 5, "cup": 5, "rim": 5})
    with pytest.raises(ValueError):
        build_target({"nodule": np.ones((5, 5, 2), bool)}, ("nodule",), None, {"nodule": 5})
    with pytest.raises(ValueError):  # more than MAX_RINGS rings: the text would not parse
        build_target({"nodule": np.kron(np.ones((4, 4), bool), [[1, 0, 0], [0, 0, 0], [0, 0, 0]])}, ("nodule",), None,
                     {"nodule": 1})


def test_min_areas_is_the_rule_of_each_union_and_of_the_rim():
    d = np.stack([disc(80, 40, 40, r) for r in (30, 33)])
    c = np.stack([disc(80, 40, 40, r) for r in (12, 20)])
    assert min_areas({"disc": d, "cup": c}) == {"disc": min_area(int(d[1].sum())), "cup": min_area(int(c[1].sum())),
                                                 "rim": min_area(int((d[1] & ~c[0]).sum()))}
    assert min_areas({"nodule": c}) == {"nodule": min_area(int(c[1].sum()))}


def test_target_set_hash():
    rows = [("a", "1", '{"nodule":[]}'), ("b", "2", '{"nodule":[]}')]
    h = target_set_hash(rows)
    assert len(h) == 64 and h == target_set_hash(rows[::-1])
    assert h != target_set_hash([rows[0], ("b", "2", '{"nodule":[ ]}')])
    assert h != target_set_hash([rows[0], ("b", "3", '{"nodule":[]}')])


GOLDEN = Path(__file__).parent / "goldens" / "targets_lab.npz"


def test_targets_lab_golden_holds_no_masks():
    """CLEAN4 licence rule: per-case LIDC records, RIGA only as three totals, nothing 2-D, no bool, no pickle."""
    z = np.load(GOLDEN, allow_pickle=False)
    lidc = [k for k in z.files if k.startswith("lidc/")]
    assert len(lidc) == 329 and all(re.fullmatch(r"lidc/[^/]+/[^/]+", k) for k in lidc)
    assert sorted(set(z.files) - set(lidc)) == ["riga_cases", "riga_sha", "riga_vertices"]
    assert all(z[k].ndim <= 1 and z[k].dtype != bool and z[k].size <= 4096 for k in z.files)
    assert GOLDEN.stat().st_size < 1_000_000


@pytest.mark.data
def test_fixtures_differential_vs_lab_golden(fixtures_dir):
    """R1 target part on fixture crops against what the lab builder recorded (tests/goldens/targets_lab.npz): the
    crop hashes match, LIDC vertex counts are equal and texts identical except where a pinch is bevelled along its
    crack edges (1 of 329 here), and RIGA objects are never longer in total."""
    z = np.load(GOLDEN, allow_pickle=False)
    seen, same, prod_v, riga_n, riga_h = set(), 0, 0, 0, hashlib.sha256()
    for ds, step in (("riga", 20), ("lidc", 6)):
        zf = np.load(fixtures_dir / f"{ds}_experts.npz")
        names = ("disc", "cup") if ds == "riga" else ("nodule",)
        for it in sorted({k.rsplit("__", 1)[0] for k in zf.files})[::step]:
            st_ = {n: zf[f"{it}__{n}"].astype(bool) for n in names}
            for r in range(st_[names[0]].shape[0]):
                m = pad_div({n: st_[n][r] for n in names})
                t, _ = check(m, names)
                if ds == "riga":
                    riga_n += 1
                    riga_h.update(mask_sha(m).encode())
                    prod_v += sum(t.info[n]["vertices"] for n in names)
                    continue
                key = f"lidc/{it}/{r}"
                seen.add(key)
                assert mask_sha(m) == z[key]["crop_sha"], f"{key}: the ported crop differs from the lab crop"
                assert t.info["nodule"]["vertices"] == z[key]["vertices"], key
                same += hashlib.sha256(t.text.encode()).hexdigest() == z[key]["text_sha"]
    assert seen == {k for k in z.files if k.startswith("lidc/")}
    assert riga_n == z["riga_cases"] and riga_h.hexdigest() == z["riga_sha"]
    print("same", same, "lidc", len(seen), "lab_v", int(z["riga_vertices"]), "prod_v", prod_v)
    assert 100 * same >= 99 * len(seen) and prod_v <= z["riga_vertices"]
