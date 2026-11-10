"""D4b (B17, B8; gates D-G3, D-G4, D-G11): objdel.geom.topology on the engine's pixel label chain (gc.topology)
equals the OpenCV versions it replaced (tests/oracles/topology_oracle.py, cv2_*) field by field, in value and type,
at every min_area; the batched forms are one engine call; certify takes b0, b1 from it; no OpenCV or scipy labelling
runs in the target builder, the cache or certify."""
from pathlib import Path

import numpy as np
import pytest
from oracles.topology_oracle import cv2_betti, cv2_components, cv2_instance_topology, cv2_signature, cv2_size_filter

from objdel.geom import engine as _gpucheck
from objdel.geom import topology as T
from objdel.geom.certify import certify_many
from objdel.geom.topology import betti, instance_topology, signature, signatures, size_filter, topologies
from tests.oracles import cuda_parity as C

MIN_AREAS = (0, 1, 2, 3, 7, 10, 100, 10**6, 10**9)


def corpus(n=300, seed=0):
    """Fuzz families (holes, diagonals, checkerboards, noise, edge cases incl. 0 x 0) plus nested frames."""
    out = [m for _, m in C.fuzz_masks(n, seed)]
    sq = np.zeros((21, 21), bool)
    for i in range(0, 10, 2):
        sq[i:21 - i, i:21 - i] = True
        sq[i + 1:20 - i, i + 1:20 - i] = False
    return out + [sq, np.eye(9, dtype=bool), np.ones((1, 40), bool), np.zeros((7, 0), bool)]


def engine(masks, min_area, want_masks=False):
    ms = [np.ascontiguousarray(m, np.uint8) for m in masks]
    start = np.zeros(len(ms) + 1, np.int64)
    start[1:] = np.cumsum([m.size for m in ms])
    bits = np.concatenate([m.ravel() for m in ms]) if ms else np.zeros(0, np.uint8)
    hw = np.array([m.shape for m in ms], np.int32).reshape(-1, 2)
    return _gpucheck.topology(bits, hw, start, np.asarray(min_area, np.int64), want_masks)


def fields(out, i):
    """Mask i of a gc.topology result as cv2_components gives it."""
    cs, hs = out["comp_start"], out["hole_start"]
    return [(int(out["comp_area"][c]), int(out["comp_holes"][c]), out["hole_area"][hs[c]:hs[c + 1]].tolist())
            for c in range(cs[i], cs[i + 1])]


def test_engine_fields_equal_the_cv2_oracle():
    """Per mask: components in raster order of their first pixel with area and hole count, each component's hole
    areas, the signature, and (want_masks) the filtered mask; min_area varies per mask in one call."""
    masks = corpus()
    rng = np.random.default_rng(0)
    ma = [int(rng.choice(MIN_AREAS)) for _ in masks]
    out = engine(masks, ma, want_masks=True)
    assert set(out) == {"comp_start", "comp_area", "comp_holes", "hole_start", "hole_area", "sig_start", "sig",
                        "masks"}
    starts = np.r_[0, np.cumsum([m.size for m in masks])]
    for i, m in enumerate(masks):
        assert fields(out, i) == cv2_components(m, ma[i]), i
        assert tuple(out["sig"][out["sig_start"][i]:out["sig_start"][i + 1]].tolist()) == cv2_signature(m, ma[i])
        got = out["masks"][starts[i]:starts[i + 1]].reshape(m.shape)
        assert np.array_equal(got, cv2_size_filter(m, ma[i]).view(np.uint8)), i
    assert "masks" not in engine(masks[:3], [0, 0, 0])
    assert engine([], [])["comp_start"].tolist() == [0]


def test_public_functions_equal_the_cv2_versions_in_value_and_type():
    """P06: instance_topology, betti, size_filter and signature on 500 masks x 9 rules (a sample per mask)."""
    masks = corpus(500, seed=1)
    rng = np.random.default_rng(1)
    for m in masks:
        a, b = instance_topology(m), cv2_instance_topology(m)
        for x, y in zip((a.areas, a.holes, a.hole_areas), (b.areas, b.holes, b.hole_areas), strict=True):
            assert x.dtype == y.dtype == np.int64 and x.tolist() == y.tolist()
        assert betti(m) == cv2_betti(m) and all(type(v) is int for v in betti(m))
        for k in rng.choice(MIN_AREAS, 3, replace=False):
            k = int(k)
            f, g = size_filter(m, k), cv2_size_filter(m, k)
            assert f.dtype == g.dtype == bool and f.shape == g.shape and np.array_equal(f, g)
            s = signature(m, k)
            assert s == cv2_signature(m, k) and type(s) is tuple and all(type(v) is int for v in s)


def test_size_filter_literals():
    """U11: component A (2 px) and B (a 3 x 3 frame, 8 px, one 1-px hole): rules 0..9; a 9 x 9 square with a hole of
    m - 1, m and m + 1 px at rule m; a 1-px border component; three holes; nothing."""
    m = np.zeros((5, 7), bool)
    m[0, 0:2] = True  # A
    m[1:4, 3:6] = True
    m[2, 4] = False  # B and its hole
    want = {0: (1, 0), 1: (1, 0), 2: (0, 0), **{k: (0,) for k in range(3, 9)}, 9: ()}
    for k, sig in want.items():
        assert signature(m, k) == sig == cv2_signature(m, k), k
    assert size_filter(m, 2).sum() == 11 and size_filter(m, 3).sum() == 9 and not size_filter(m, 9).any()
    for hole in (8, 9, 10):  # rule 9
        sq = np.ones((9, 9), bool)
        sq.flat[[10 + (j // 7) * 9 + j % 7 for j in range(hole)]] = False
        assert signature(sq, 9) == ((1,) if hole >= 9 else (0,)) == cv2_signature(sq, 9), hole
    dot = np.zeros((6, 6), bool)
    dot[0, 5] = True
    assert signature(dot, 1) == (0,) and signature(dot, 2) == () and instance_topology(dot).areas.tolist() == [1]
    three = np.ones((5, 9), bool)
    three[2, 1] = three[2, 4] = three[2, 7] = False
    assert signature(three, 0) == (3,)
    assert signature(np.c_[three, np.zeros((5, 1), bool), np.ones((5, 1), bool)]) == (3, 0)
    assert signature(np.zeros((4, 4), bool), 0) == () and signature(np.zeros((0, 0), bool)) == ()


class Spy:
    def __init__(self, monkeypatch):
        self.n, real = 0, _gpucheck.topology

        def spy(*a, **k):
            self.n += 1
            return real(*a, **k)
        monkeypatch.setattr(T._gpucheck, "topology", spy)


def test_batched_forms_are_one_engine_call(monkeypatch):
    """I11: signatures and topologies of 200 masks (empty and 0 x 0 included) are one gc.topology call each and equal
    the per-mask results."""
    masks = corpus(200, seed=2)[:200] + [np.zeros((3, 3), bool), np.zeros((0, 0), bool)]
    each_sig = [signature(m, 3) for m in masks]
    each_top = [instance_topology(m) for m in masks]
    spy = Spy(monkeypatch)
    assert signatures(masks, 3) == each_sig
    got = topologies(masks)
    assert spy.n == 2
    assert [(t.areas.tolist(), t.holes.tolist(), t.hole_areas.tolist()) for t in got] == [
        (t.areas.tolist(), t.holes.tolist(), t.hole_areas.tolist()) for t in each_top]
    assert signatures([], 5) == [] and topologies([]) == []


def test_certify_takes_betti_from_the_pixel_chain(monkeypatch):
    """B8 leg 2: certify_many makes one gc.topology call for the b0, b1 of all masks (no OpenCV labelling)."""
    masks = corpus(40, seed=3)[:40]
    want = [(r.b0, r.b1) for r in certify_many(masks)]
    assert want == [cv2_betti(m) if m.size else (0, 0) for m in masks]
    spy = Spy(monkeypatch)
    certify_many(masks)
    assert spy.n == 1


def test_no_cpu_labelling_in_targets_cache_and_certify(monkeypatch):
    """S01 (D4b part): with OpenCV's connected-component functions and scipy.ndimage.label raising, build_target on a
    LIDC-like and a RIGA-like fixture, ExpertCache.build and certify run and give the unpatched results."""
    import cv2
    from scipy import ndimage

    from objdel.data.targets.build import build_target, min_areas
    from objdel.geom.certify import certify
    from objdel.verify.reward.cache import ExpertCache

    y, x = np.ogrid[:80, :80]
    disc, cup = (y - 40) ** 2 + (x - 40) ** 2 <= 900, (y - 40) ** 2 + (x - 44) ** 2 <= 144
    nod = disc & ~((y - 40) ** 2 + (x - 40) ** 2 <= 30)
    stacks = {"disc": np.stack([disc, disc]), "cup": np.stack([cup, cup & ~nod])}

    def run():
        return (build_target({"nodule": nod}, ("nodule",), None, min_areas({"nodule": nod[None]})).text,
                build_target({"disc": disc, "cup": cup}, ("disc", "cup"), None, min_areas(stacks)).text,
                ExpertCache.build(stacks["cup"]).rater_sigs, certify(nod, 5))
    want = run()

    def boom(*a, **k):
        raise AssertionError("CPU labelling called")
    for name in ("connectedComponents", "connectedComponentsWithAlgorithm", "connectedComponentsWithStats",
                 "connectedComponentsWithStatsWithAlgorithm"):
        monkeypatch.setattr(cv2, name, boom)
    monkeypatch.setattr(ndimage, "label", boom)
    assert run() == want


@pytest.mark.data
@pytest.mark.slow
def test_lidc_crops_equal_the_cv2_oracle():
    """The 37,093 LIDC crops of data/checker_eval/masks.npz at rules 0, 5 and each crop's own rule: one call each."""
    path = Path(__file__).resolve().parents[1] / "data" / "checker_eval" / "masks.npz"
    if not path.is_file():
        pytest.skip(f"{path} not found")
    z = dict(np.load(path, allow_pickle=False))
    masks = [z["bits"][z["start"][k]:z["start"][k + 1]].reshape(z["hw"][k]).astype(bool) for k in range(len(z["hw"]))]
    assert len(masks) == 37093
    for rule in (0, 5, None):
        ma = [T.min_area(int(m.sum())) if rule is None else rule for m in masks]
        out = engine(masks, ma)
        assert [i for i, m in enumerate(masks) if fields(out, i) != cv2_components(m, ma[i])] == []


@pytest.mark.data
def test_rater_masks_equal_the_cv2_oracle(store_dir):
    """1,000 LIDC and RIGA rater masks of the store, framed as the builders frame them, at rules 0, 5 and the item's
    rule (Mac only for RIGA)."""
    from itertools import islice

    from objdel.data import store
    from objdel.data.targets.build import min_areas
    from objdel.data.targets.frame import frame

    masks, rules = [], []
    for ds, n in (("lidc", 250), ("riga", 20)):
        for it in islice(store.load(ds, "train", store_dir), n):
            f = frame(it)
            ma = min_areas(f.masks)
            for name, st in f.masks.items():
                masks += list(st)
                rules += [ma[name]] * len(st)
    assert len(masks) >= 1000
    for rule in (0, 5, None):
        ma = rules if rule is None else [rule] * len(masks)
        out = engine(masks, ma)
        assert [i for i, m in enumerate(masks) if fields(out, i) != cv2_components(m, ma[i])] == []
