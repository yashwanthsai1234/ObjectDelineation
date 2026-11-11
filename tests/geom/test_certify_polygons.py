"""D4 (B10, B11; gates D-G3, D-G4): the engine's ladder certificate `certify_polygons` against set algebra on an
independent raster (tests/oracles/raster_oracle.py: even-odd crossings at pixel centres), the scipy signature oracle
and check_polygons' legality; its rim pair queries and argument checks; build_target's ladder through it (one call
per name and level) and min_area required. Expected values never come from the
engine call under test."""
import random

import numpy as np
import pytest
from oracles.raster_oracle import oracle_winding
from oracles.topology_oracle import signature as oracle_signature

from objdel.data.targets import build as B
from objdel.data.targets.build import BuildConfig, build_target, min_areas
from objdel.geom import engine as _gpucheck
from objdel.geom.shape import pack_objects
from objdel.geom.text import GRID, MAX_RINGS
from tests.oracles import cuda_parity as C
from tests.oracles.kept_oracle import kept as kept_oracle


def pack_masks(ms):
    ms = [np.ascontiguousarray(m, np.uint8) for m in ms]
    start = np.zeros(len(ms) + 1, np.int64)
    start[1:] = np.cumsum([m.size for m in ms])
    bits = np.concatenate([m.ravel() for m in ms]) if ms else np.zeros(0, np.uint8)
    return bits, np.array([m.shape for m in ms], np.int32).reshape(-1, 2), start


def run(objs, refs, ref_of, min_area, G=GRID, pairs=None, ref_pairs=None, pair_min_area=None):
    """certify_polygons of objects [(parts, W, H)] against reference masks."""
    i32 = lambda a: None if a is None else np.asarray(a, np.int32).reshape(-1, 2)
    return _gpucheck.certify_polygons(
        *pack_objects(objs), G, *pack_masks(refs), np.asarray(ref_of, np.int32), np.asarray(min_area, np.int64),
        pairs=i32(pairs), ref_pairs=i32(ref_pairs),
        pair_min_area=None if pair_min_area is None else np.asarray(pair_min_area, np.int64))


def raster(parts, W, H, G=GRID):
    """Even-odd raster of an object's rings from the independent oracle."""
    rings = [np.asarray(r).reshape(-1, 2).tolist() for o, hs in parts for r in (o, *hs)]
    return oracle_winding(rings, W, H, G)[1] % 2 == 1 if rings else np.zeros((H, W), bool)


def legal(objs, G=GRID):
    return [not d.any() for d in _gpucheck.check_polygons(*pack_objects(objs), G)["defects"]]


def sl(out, key, i):
    s = out[key + "_start"]
    return tuple(out[key][s[i]:s[i + 1]].tolist())


def square(x0, y0, x1, y1, hole=False):
    r = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
    return r[::-1] if hole else r


def disc_ring(cx, cy, r, n=48):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.c_[cx + r * np.cos(t), cy + r * np.sin(t)].round().astype(np.int64)


def check_against_oracles(objs, refs, ref_of, min_area, out, G=GRID):
    ok = legal(objs, G)
    for i, (parts, W, H) in enumerate(objs):
        c = raster(parts, W, H, G) if ok[i] else np.zeros((H, W), bool)  # an illegal shape draws nothing
        r = np.asarray(refs[ref_of[i]], bool)
        inter, union = int((c & r).sum()), int((c | r).sum())
        assert bool(out["legal"][i]) == ok[i], i
        assert (int(out["inter"][i]), int(out["union"][i])) == (inter, union), i
        assert union == int(c.sum()) + int(r.sum()) - inter
        assert bool(out["equal"][i]) == np.array_equal(c, r), i
        assert sl(out, "sig", i) == oracle_signature(c, min_area[i]), i
        assert sl(out, "ref_sig", i) == oracle_signature(r, min_area[i]), i


def fuzz_objects(n, seed):
    """(parts, W, H, G) of the fuzz polygon families with valid roles, sides <= 512."""
    out = []
    for _, rings, roles, W, H, G in C.fuzz_polygons(20, seed=seed):
        roles = [-1] * len(rings) if roles is None else list(roles)
        if W > 512 or H > 512 or len(rings) > 10 or any(len(r) > 1000 for r in rings):
            continue
        if any(p != -1 and not (0 <= p < len(rings) and roles[p] == -1) for p in roles):
            continue
        parts = [(np.asarray(rings[i]).reshape(-1, 2), [np.asarray(rings[j]).reshape(-1, 2)
                                                        for j, p in enumerate(roles) if p == i])
                 for i, p in enumerate(roles) if p == -1]
        out.append((parts, W, H, G))
    random.Random(seed).shuffle(out)
    return out[:n]


def test_set_algebra_on_fuzz_polygons():
    """I01: legal, inter, union, equal, sig and ref_sig on 300 fuzz polygons (every family, G in {7, 1000, 2^30}),
    references random, the shape's own raster, empty or nested in it, shared by several shapes and shuffled."""
    rng = np.random.default_rng(1)
    cases = fuzz_objects(300, seed=4)
    assert len(cases) == 300
    by_g = {}
    for parts, W, H, G in cases:
        by_g.setdefault(G, []).append((parts, W, H))
    assert {7, 1000, 2**30} <= set(by_g)
    seen = set()
    for G, objs in by_g.items():
        refs = []
        for k, (parts, W, H) in enumerate(objs):
            own = raster(parts, W, H, G)
            kind = k % 4
            if kind == 0:
                refs.append(rng.random((H, W)) < rng.uniform(0.05, 0.6))
            elif kind == 1:
                refs.append(own)
            elif kind == 2:
                refs.append(np.zeros((H, W), bool))
            else:
                nested = own.copy()
                nested[H // 3:H // 3 + 2, W // 3:W // 3 + 2] ^= True
                refs.append(nested)
        order = rng.permutation(len(refs))
        refs = [refs[j] for j in order]
        ref_of = np.argsort(order).astype(np.int32)
        for k in range(1, len(objs), 7):  # a few shapes share the previous shape's reference
            if objs[k][1:] == objs[k - 1][1:]:
                ref_of[k] = ref_of[k - 1]
        min_area = [(0, 1, 4, 20)[k % 4] for k in range(len(objs))]
        out = run(objs, refs, ref_of, min_area, G)
        check_against_oracles(objs, refs, ref_of, min_area, out, G)
        seen |= {(bool(a), bool(b)) for a, b in zip(out["legal"], out["equal"], strict=True)}
    assert seen == {(True, True), (True, False), (False, True), (False, False)}  # illegal vs empty reference: equal


def test_identical_disjoint_empty_and_illegal():
    """U14: equal only for equal masks, both empty is equal with inter = union = 0, superset, disjoint, and an
    illegal shape (a bow-tie) has legal False and draws nothing."""
    W = H = 10
    sq = [(square(200, 200, 600, 600), [])]  # pixels 2..5 x 2..5
    own = np.zeros((H, W), bool)
    own[2:6, 2:6] = True
    big = np.zeros((H, W), bool)
    big[1:8, 1:8] = True
    far = np.zeros((H, W), bool)
    far[8:, 8:] = True
    bow = [(np.array([[100, 100], [800, 800], [800, 100], [100, 800]]), [])]
    objs = [(sq, W, H), ([], W, H), (sq, W, H), ([], W, H), (sq, W, H), (sq, W, H), (bow, W, H)]
    refs = [own, np.zeros((H, W), bool), big, far]
    ref_of = [0, 1, 1, 0, 3, 2, 0]
    out = run(objs, refs, ref_of, [0] * 7)
    assert out["legal"].tolist() == [True] * 6 + [False]
    assert out["equal"].tolist() == [True, True, False, False, False, False, False]
    assert out["inter"].tolist() == [16, 0, 0, 0, 0, 16, 0]
    assert out["union"].tolist() == [16, 0, 16, 16, 20, 49, 16]
    assert sl(out, "sig", 6) == () and sl(out, "ref_sig", 6) == (0,)
    check_against_oracles(objs, refs, ref_of, [0] * 7, out)


def test_rim_pairs():
    """I02: pair_sig == signature(a & ~b, pm) of the rasters and pair_ref_sig of the references (and-not, not xor):
    cup partly outside the disc, cup == disc, empty cup, a == b; shapes without a pair; pm in {0, 5, 50}."""
    S = 64
    disc = [(disc_ring(500, 500, 400), [])]
    cups = {"inside": [(disc_ring(520, 480, 180), [])], "outside": [(disc_ring(800, 500, 300), [])],
            "same": disc, "empty": []}
    objs = [(disc, S, S)] + [(c, S, S) for c in cups.values()] + [(cups["inside"], S, S)]
    rasters = [raster(p, S, S) for p, _, _ in objs]
    rng = np.random.default_rng(2)
    refs = [rasters[0] ^ (rng.random((S, S)) < 0.01), rasters[1], rasters[2] | rasters[0], rasters[0]]
    pairs = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 1)]
    ref_pairs = [(0, 1), (0, 2), (0, 3), (0, 0), (2, 1)]
    for pm in (0, 5, 50):
        out = run(objs, refs, [0, 1, 2, 3, 0, 1], [3] * 6, pairs=pairs, ref_pairs=ref_pairs,
                  pair_min_area=[pm] * len(pairs))
        for p, ((a, b), (ra, rb)) in enumerate(zip(pairs, ref_pairs, strict=True)):
            assert sl(out, "pair_sig", p) == oracle_signature(rasters[a] & ~rasters[b], pm), (pm, p)
            assert sl(out, "pair_ref_sig", p) == oracle_signature(refs[ra] & ~refs[rb], pm), (pm, p)
        assert sl(out, "pair_sig", 1) == oracle_signature(rasters[0] & ~rasters[2], pm) != ()  # cup partly outside
    out = run(objs[:1], refs[:1], [0], [0])
    for k in ("pair_sig", "pair_ref_sig"):
        assert out[k + "_start"].tolist() == [0] and out[k].size == 0


def test_argument_checks():
    """F02 (certify_polygons part): bad indices, sizes, min_area and the polygon caps raise ValueError."""
    objs = [([(square(0, 0, 500, 500), [])], 8, 8)]
    ref = [np.ones((8, 8), bool)]
    run(objs, ref, [0], [0], G=1)
    run(objs, ref, [0], [0], G=2**30)
    with pytest.raises(ValueError):
        run(objs, ref, [1], [0])  # ref_of out of range
    with pytest.raises(ValueError):
        run(objs, ref, [-1], [0])
    with pytest.raises(ValueError):
        run(objs, ref, [0], [-1])  # min_area < 0
    with pytest.raises(ValueError):
        run(objs, [np.ones((8, 9), bool)], [0], [0])  # reference size differs from the shape's frame
    with pytest.raises(ValueError):
        run(objs, ref, [0, 0], [0])  # one ref_of per shape
    with pytest.raises(ValueError):
        run(objs, ref, [0], [0, 0])  # one min_area per shape
    for G in (0, 2**30 + 1):
        with pytest.raises(ValueError):
            run(objs, ref, [0], [0], G=G)
    for side in (0, 4097):
        with pytest.raises(ValueError):
            run([([(square(0, 0, 500, 500), [])], side, 8)], [np.ones((8, max(side, 1)), bool)], [0], [0])
    two = objs + [([(square(0, 0, 500, 500), [])], 8, 9)]
    refs2 = ref + [np.ones((9, 8), bool)]
    run(two, refs2, [0, 1], [0, 0], pairs=[(0, 0)], ref_pairs=[(0, 0)], pair_min_area=[0])
    for kw in ({"pairs": [(0, 2)], "ref_pairs": [(0, 0)], "pair_min_area": [0]},  # shape index out of range
               {"pairs": [(0, 0)], "ref_pairs": [(0, 2)], "pair_min_area": [0]},  # reference index out of range
               {"pairs": [(0, 1)], "ref_pairs": [(0, 0)], "pair_min_area": [0]},  # a pair of two frames
               {"pairs": [(0, 0)], "ref_pairs": [(0, 1)], "pair_min_area": [0]},
               {"pairs": [(0, 0)], "ref_pairs": [(0, 0)], "pair_min_area": [-1]},
               {"pairs": [(0, 0)], "ref_pairs": [(0, 0), (0, 0)], "pair_min_area": [0]},
               {"pairs": [(0, 0)], "ref_pairs": None, "pair_min_area": None}):
        with pytest.raises(ValueError):
            run(two, refs2, [0, 1], [0, 0], **kw)


def test_signature_equality_is_the_kept_oracle():
    """P07: the ladder's topology test (sig == ref_sig at the reward's rule) equals the archived interval loop
    (tests/oracles/kept_oracle.py) at that rule, on fuzz polygons against random and near-copy references."""
    rng = np.random.default_rng(3)
    objs = [(p, W, H) for p, W, H, G in fuzz_objects(120, seed=9) if G == GRID]
    ok = legal(objs)
    refs, want, ma = [], [], []
    for k, (parts, W, H) in enumerate(objs):
        c = raster(parts, W, H) if ok[k] else np.zeros((H, W), bool)
        r = c ^ (rng.random((H, W)) < 0.02) if k % 2 else rng.random((H, W)) < 0.3
        m = int(rng.integers(0, 12))
        refs.append(r), ma.append(m), want.append(kept_oracle(c, r, m, m))
    out = run(objs, refs, range(len(objs)), ma)
    got = [sl(out, "sig", i) == sl(out, "ref_sig", i) for i in range(len(objs))]
    assert got == want and any(want) and not all(want)


class Spy:
    """Counts build_target's certify_polygons calls: the number of shapes per call."""

    def __init__(self, monkeypatch):
        self.calls, real = [], _gpucheck.certify_polygons

        def spy(*a, **k):
            self.calls.append(len(a[4]))
            return real(*a, **k)
        monkeypatch.setattr(B._gpucheck, "certify_polygons", spy)


def _disc(S, cy, cx, r):
    y, x = np.ogrid[:S, :S]
    return (y - cy) ** 2 + (x - cx) ** 2 <= r * r


def test_one_call_per_name_and_level(monkeypatch):
    """A nodule that passes its first rung is one call with one shape; a two-object item is one call per name."""
    m = _disc(60, 30, 30, 20)
    spy = Spy(monkeypatch)
    build_target({"nodule": m}, ("nodule",), None, min_areas({"nodule": m[None]}))
    assert spy.calls == [1]
    d, c = _disc(80, 40, 40, 30), _disc(80, 40, 44, 12)
    ma = min_areas({"disc": d[None], "cup": c[None]})
    spy = Spy(monkeypatch)
    build_target({"disc": d, "cup": c}, ("disc", "cup"), None, ma)
    assert spy.calls == [1, 1]


def test_rungs_of_one_call_never_exceed_the_ladder(monkeypatch):
    """I12: a wiggly mask climbs the ladder; every call starts at a new rung (shape counts add up to at most the
    ladder's 14 rungs) and the chosen rung is the first certified one (the D-R02 golden pins the texts)."""
    rng, S = np.random.default_rng(0), 500
    y, x = np.mgrid[:S, :S]
    m = np.hypot(y - 250, x - 250) <= 120 + 6 * np.sin(23 * np.arctan2(y - 250, x - 250)) + rng.normal(0, 0.3, (S, S))
    spy = Spy(monkeypatch)
    t = build_target({"nodule": m}, ("nodule",), None, {"nodule": 5})  # exact: rung 9 (tau 0.75)
    assert t.info["nodule"]["tau"] == 0.75 and spy.calls == [1, 2, 4, 7]  # rungs 0, 1-2, 3-6, 7-13
    assert sum(spy.calls) <= len(BuildConfig().ladder())


def test_min_area_is_required():
    """F06 (B10): a missing name raises ValueError naming the key; negative raises;
    0 and extra keys are accepted."""
    m = _disc(40, 20, 20, 10)
    with pytest.raises(ValueError, match="'nodule'"):
        build_target({"nodule": m}, ("nodule",))
    with pytest.raises(ValueError, match="'nodule'"):
        build_target({"nodule": m}, ("nodule",), None, {})
    with pytest.raises(ValueError):
        build_target({"nodule": m}, ("nodule",), None, {"nodule": -1})
    a = build_target({"nodule": m}, ("nodule",), None, {"nodule": 0})
    assert build_target({"nodule": m}, ("nodule",), None, {"nodule": 0, "cup": 7}).text == a.text


@pytest.mark.data
def test_ladder_rungs_of_real_builds_equal_the_host_path(store_dir):
    """Step 1 data variant: every ladder rung of real LIDC builds certified in one call equals the host path (raster
    from check_polygons, numpy IoU, the signature oracle)."""
    from itertools import islice

    from objdel.data import store
    from objdel.data.targets.frame import frame
    from objdel.geom.shape import to_masks
    from objdel.geom.text import MAX_VERTICES

    ladder, rungs = BuildConfig().ladder(), 0
    for it in islice(store.load("lidc", "train", store_dir), 100):
        f = frame(it)
        names, ma = ("nodule",), min_areas(f.masks)
        r = 0
        src = {k: f.masks[k][r] for k in names}
        H, W = src[names[0]].shape
        parts = {}
        for k in names:
            rings = B._rings(src[k])
            if len(rings) > MAX_RINGS:
                continue
            parts[k] = [p for p in (B._parts(rings, tau, W, H) for tau in ladder)
                        if all(len(x) <= MAX_VERTICES for o, hs in p for x in (o, *hs))]
            objs = [(p, W, H) for p in parts[k]]
            out = run(objs, [src[k]], [0] * len(objs), [ma[k]] * len(objs))
            for i, g in enumerate(to_masks(objs)):
                c = g.mask if g.legal else np.zeros_like(g.mask)
                assert bool(out["legal"][i]) == g.legal
                assert int(out["inter"][i]) == int((c & src[k]).sum())
                assert int(out["union"][i]) == int((c | src[k]).sum())
                assert bool(out["equal"][i]) == np.array_equal(c, src[k])
                assert sl(out, "sig", i) == oracle_signature(c, ma[k])
                assert sl(out, "ref_sig", i) == oracle_signature(src[k], ma[k])
            rungs += len(objs)
    assert rungs >= 1000
