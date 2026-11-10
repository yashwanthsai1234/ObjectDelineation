"""Task 2: band formulas, tolerances, distance transforms, reward terms."""
import numpy as np
import pytest
from scipy import ndimage

from rewardcore import fixtures, perturb
from rewardcore.band import (FORMULAS, P95_PAIR_OVER_SIGMA, excursion_second_moment, tolerance_factor)
from rewardcore.distances import boundary, edt, rater_tau, symmetric_surface_distances
from rewardcore.reward import ExpertCache, score


def disk(n, r, cy=None, cx=None):
    y, x = np.mgrid[:n, :n]
    cy = n / 2 - 0.5 if cy is None else cy
    cx = n / 2 - 0.5 if cx is None else cx
    return np.hypot(y - cy, x - cx) <= r


def test_edt_matches_scipy():
    rng = np.random.default_rng(0)
    for _ in range(50):
        m = rng.random((int(rng.integers(5, 80)), int(rng.integers(5, 80)))) < 0.7
        if m.all():
            continue
        assert np.allclose(edt(m), ndimage.distance_transform_edt(m), atol=1e-4)


def test_boundary_is_four_neighbour_contour():
    rng = np.random.default_rng(1)
    for _ in range(50):
        m = rng.random((30, 40)) < 0.6
        ref = m & ~ndimage.binary_erosion(np.pad(m, 1), ndimage.generate_binary_structure(2, 1))[1:-1, 1:-1]
        assert np.array_equal(boundary(m), ref)


def test_surface_distances_and_tau():
    a = disk(60, 15)
    assert symmetric_surface_distances(a, a).max() == 0
    b = disk(60, 18)
    d = symmetric_surface_distances(a, b)
    assert 2 <= np.percentile(d, 95) <= 4
    tau, raw = rater_tau(np.stack([a, a]), floor=1.0)
    assert raw == 0 and tau == 1.0
    tau, raw = rater_tau(np.stack([a]), floor=1.5)
    assert np.isnan(raw) and tau == 1.5


def test_excursion_moments_match_monte_carlo():
    rng = np.random.default_rng(2)
    n = 400_000
    assert excursion_second_moment(1) == pytest.approx(2.0, abs=1e-9)  # E[(X - Y)^2] = 2
    for m in (2, 3, 5):
        y = rng.standard_normal((n, m))
        x = rng.standard_normal(n)
        mc = np.mean(np.maximum(0, x - y.max(1)) ** 2 + np.maximum(0, y.min(1) - x) ** 2)
        assert excursion_second_moment(m) == pytest.approx(mc, rel=0.02)
    g = [tolerance_factor(m) for m in range(1, 9)]
    assert all(a > b for a, b in zip(g, g[1:]))
    assert P95_PAIR_OVER_SIGMA == pytest.approx(2.7718, abs=1e-4)


def test_perfect_and_inside_band_scores_one():
    stack = np.stack([disk(80, 18), disk(80, 20), disk(80, 22)])
    cache = ExpertCache.build(stack)
    for m in (stack[0], stack[1], stack[2], cache.vote, cache.union, cache.core):
        t = score(m, cache)
        assert all(t.band[f] == 1.0 for f in FORMULAS)
        assert t.topo_any == 1.0


@pytest.mark.parametrize("formula", ["F1", "F2", "F3", "F4", "F5"])
def test_scale_free_for_uniform_excursion(formula):
    """Same excursion in units of tau at radius 10..80 -> nearly the same score; F0 is not."""
    vals = []
    for r in (10, 20, 40, 80):
        n = 2 * r + 40
        stack = np.stack([disk(n, r - 1), disk(n, r + 1)])
        cache = ExpertCache.build(stack)
        e = 0.5 * cache.tau
        vals.append(score(perturb.dilate(cache.union, e), cache).band[formula])
    assert max(vals) - min(vals) < 0.12, vals


def test_current_band_is_scale_biased():
    vals = []
    for r in (10, 80):
        n = 2 * r + 40
        stack = np.stack([disk(n, r - 1), disk(n, r + 1)])
        cache = ExpertCache.build(stack)
        vals.append(score(perturb.dilate(cache.union, 2), cache).band["F0"])
    assert vals[1] - vals[0] > 0.1


def test_monotone_in_dilation_erosion_shift():
    stack = fixtures.load("lidc")
    stack = next(v for v in stack.values() if len(v) == 4 and v[0].sum() > 150)
    cache = ExpertCache.build(stack)
    for fam in ("grow", "shrink", "shift"):
        prev = {f: 1.0 + 1e-9 for f in FORMULAS}
        for d in (1, 2, 3, 4):
            m = {"grow": perturb.dilate(cache.vote, d), "shrink": perturb.erode(cache.vote, d),
                 "shift": perturb.shift(cache.vote, d, 0)}[fam]
            t = score(m, cache)
            for f in FORMULAS:
                assert t.band[f] <= prev[f] + 1e-12, (fam, d, f)
                prev[f] = t.band[f]


def test_roi_crop_is_exact():
    """score() crops to bbox(pred | union); compare with a full-frame computation."""
    rng = np.random.default_rng(3)
    for st in list(fixtures.load("lidc").values())[:40]:
        cache = ExpertCache.build(st)
        for m in (st[0], perturb.add_speck(st[0], 4, 9, rng), perturb.shift(st[0], 3, -2), np.zeros_like(st[0])):
            if m is None:
                continue
            t = score(m, cache)
            spill, miss = m & ~cache.union, cache.core & ~m
            depth = (cache.dist_out_union[spill].astype(np.float64).sum() - 0.5 * spill.sum()
                     + cache.dist_in_core[miss].astype(np.float64).sum() - 0.5 * miss.sum())
            inter = (m & cache.vote).sum()
            uni = (m | cache.vote).sum()
            assert t.iou == (inter / uni if uni else 1.0)
            assert (t.n_spill, t.n_miss) == (spill.sum(), miss.sum())
            assert t.depth == pytest.approx(depth, rel=1e-12, abs=1e-9)


def test_far_speck_costs_more_than_near_speck():
    stack = np.stack([disk(120, 20), disk(120, 22)])
    cache = ExpertCache.build(stack)
    near, far = cache.vote.copy(), cache.vote.copy()
    c = 60
    near[c, c + 22 + 3] = True
    far[c, c + 22 + 30] = True
    tn, tf = score(near, cache), score(far, cache)
    assert tf.band["F5"] < tn.band["F5"] < 1.0
    assert tn.band["F1"] == tf.band["F1"]  # F1 is not distance-weighted
    assert tn.topo_any == tf.topo_any == 1.0  # 1-px specks are below the signature size filter (R6)
    assert tn.sig_raw == tf.sig_raw == (0, 0)


def test_cache_crop_is_exact():
    """tau and the F5 calibration are computed on the union's bounding box: same values bit for bit."""
    from itertools import combinations

    from rewardcore.band import Band, self_calibrated_sq

    for ds in ("riga", "lidc"):
        for st in list(fixtures.load(ds).values())[:15]:
            d = np.concatenate([symmetric_surface_distances(a, b) for a, b in combinations(st, 2)])
            assert rater_tau(st, 0.0)[1] == float(np.percentile(d, 95))
            q = []
            for j in range(len(st)):
                b = Band.of(np.delete(st, j, axis=0))
                q.append(2.0 * b.errors(st[j]).depth / b.perimeter)
            assert self_calibrated_sq(st) == float(np.mean(q))
