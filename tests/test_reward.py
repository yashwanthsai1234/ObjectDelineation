"""objdel.reward: reward v4 cache, terms, completion scoring, advantages, determinism (plan §4.3; no data needed).

Data-backed gates (K8, K8b, K9, R3, R4, R5, local gradient) are in test_reward_gates.py."""
import json
import math
import multiprocessing
import pickle
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from itertools import product

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from oracles import reward_harness as H
from oracles import reward_reference as ref

import objdel.reward as R
from objdel.geometry import to_mask
from objdel.reward import (
    REWARD_V4,
    ExpertCache,
    advantages_v4,
    loo_bias,
    rl_eligible,
    score_batch,
    score_completion,
    score_mask,
    topology_reward,
)
from objdel.reward import cache as cache_mod


def disk(n, r, cy=None, cx=None):
    y, x = np.mgrid[:n, :n]
    return np.hypot(y - (n / 2 - 0.5 if cy is None else cy), x - (n / 2 - 0.5 if cx is None else cx)) <= r


def ring(n, r_out, r_in):
    return disk(n, r_out) & ~disk(n, r_in)


@pytest.fixture(scope="module")
def three():
    """Three concentric disks of radius 18, 19, 20 on a 60 x 60 crop."""
    st_ = np.stack([disk(60, 18), disk(60, 19), disk(60, 20)])
    return st_, ExpertCache.build(st_)


@pytest.fixture(scope="module")
def mixed():
    """3 raters draw a disk, 1 rater a ring (hole of 3 px radius): the vote has no hole."""
    d, r = disk(40, 12), ring(40, 12, 3)
    st_ = np.stack([d, d, d, r])
    return st_, ExpertCache.build(st_)


# ------------------------------------------------------------------ frozen constants and the cache

def test_frozen_constants_match_plan():
    v = REWARD_V4
    assert (v.weights, v.floors, v.c_topo, v.lam, v.eps, v.group, v.drop_spread) == (
        (0.3, 0.3, 0.1, 0.3), (0.05, 0.05, 0.25, 0.05), 2.0, 0.5, 0.1, 8, 1e-2)
    c = R.CACHE_SPEC
    assert (v.iou_gate, c.expert_score, c.tau_floor, c.tau_quantile) == (0.5, 0.95, 1.0, 95.0)
    assert R.TERMS == ("iou_tc", "f5", "legal", "g_area")


def test_single_rater_items_are_not_rl_eligible():
    st_ = np.stack([disk(30, 8)])
    assert not rl_eligible(st_) and rl_eligible(np.concatenate([st_, st_]))
    with pytest.raises(ValueError):
        ExpertCache.build(st_)
    with pytest.raises(ValueError):
        ExpertCache.build(st_[0])


@pytest.mark.parametrize("stack", [
    np.stack([disk(60, 18), disk(60, 19), disk(60, 20)]),
    np.stack([disk(50, 10), disk(50, 12, 26, 23), ring(50, 12, 4)]),
    np.stack([disk(40, 6, 10, 10) | disk(40, 6, 28, 28), disk(40, 7, 11, 10)]),
    np.stack([disk(30, 14), disk(30, 13), disk(30, 15), disk(30, 12, 14, 16)]),  # touches the crop border
])
def test_cache_matches_reference(stack):
    c, r = ExpertCache.build(stack), ref.cache(stack)
    assert c.tau == pytest.approx(r["tau"], rel=1e-5) and c.t5 == pytest.approx(r["t5"], rel=1e-5)
    assert (c.perimeter, c.min_area, list(c.rater_sigs), list(c.rater_areas)) == (
        r["perimeter"], r["min_area"], r["sigs"], r["areas"])
    assert all(np.array_equal(getattr(c, a), r[a]) for a in ("vote", "union", "core"))


def test_identical_raters_hit_the_tau_and_t5_floors():
    c = ExpertCache.build(np.stack([disk(40, 10)] * 3))
    floor = math.sqrt(cache_mod.excursion_moment(3) / (1 / 0.95 - 1)) / (1.959963984540054 * math.sqrt(2))
    assert c.tau == 1.0 and c.t5 == pytest.approx(floor)


def test_excursion_moment_matches_monte_carlo_and_reference():
    rng = np.random.default_rng(2)
    assert cache_mod.excursion_moment(1) == pytest.approx(2.0, abs=1e-9)  # E[(X - Y)^2] = 2
    for m in (2, 3, 5):
        y, x = rng.standard_normal((400_000, m)), rng.standard_normal(400_000)
        mc = np.mean(np.maximum(0, x - y.max(1)) ** 2 + np.maximum(0, y.min(1) - x) ** 2)
        assert cache_mod.excursion_moment(m) == pytest.approx(mc, rel=0.02)
        assert cache_mod.excursion_moment(m) == pytest.approx(ref.s_moment(m), rel=1e-7)


def test_excursion_moment_is_frozen_and_equals_quad(monkeypatch):
    """M9 F2: S_m for m <= 8 are constants (the quad of scipy.stats differs by 1 ULP across CPUs for m = 3, 7), equal
    to a fresh quad to 1e-12 (bit-equal on the Mac that froze them); beyond the table quad is used."""
    for m in range(1, 9):
        assert cache_mod.excursion_moment(m) == pytest.approx(cache_mod._excursion_quad(m), rel=0, abs=1e-12)
    assert cache_mod.excursion_moment(3) == float.fromhex("0x1.cb7070234b5dap-2")
    assert cache_mod.excursion_moment(9) == cache_mod._excursion_quad(9)
    monkeypatch.setattr(cache_mod.integrate, "quad", lambda *a, **k: 1 / 0)  # not computed at run time
    assert [cache_mod.excursion_moment(m) for m in range(1, 9)] == [cache_mod.S_M[m] for m in range(1, 9)]


def test_edt_is_correctly_rounded_without_opencv(monkeypatch):
    """Exact EDT = the correctly rounded float32 sqrt of the integer squared distance to the nearest zero pixel
    (brute force). M9 F1: OpenCV's x86 IPP path is not (opencv#24082), so edt must not call it."""
    monkeypatch.setattr(cache_mod.cv2, "distanceTransform", lambda *a, **k: 1 / 0)
    rng = np.random.default_rng(0)
    for _ in range(30):
        m = rng.random((int(rng.integers(2, 30)), int(rng.integers(2, 30)))) < rng.uniform(0.5, 0.99)
        if m.all():
            continue
        zy, zx = np.nonzero(~m)
        yy, xx = np.mgrid[:m.shape[0], :m.shape[1]]
        d2 = ((yy[..., None] - zy) ** 2 + (xx[..., None] - zx) ** 2).min(-1)
        assert np.array_equal(cache_mod.edt(m), np.sqrt(d2.astype(np.float64)).astype(np.float32))
    full = cache_mod.edt(np.ones((4, 4), bool))
    assert (full == 1e6).all() and full.dtype == np.float32


def test_in_band_masks_score_one(three):
    st_, c = three
    for m in (*st_, c.vote, c.union, c.core):
        t = score_mask(m, c)
        assert t.f5 == 1.0 and t.topo_freq == 1.0
    for m in st_:
        assert score_mask(m, c).g_area == 1.0 and score_mask(m, c).iou_tc == 1.0


def test_save_load_roundtrip_tamper_and_fingerprint(tmp_path, three, monkeypatch):
    _, c = three
    digest = c.save(tmp_path / "c.npz")
    d = ExpertCache.load(tmp_path / "c.npz")
    assert digest == c.digest == d.digest and (d.tau, d.t5) == (c.tau, c.t5)
    for a in ("raters", "vote", "union", "core", "dist_out_union", "dist_in_core"):
        assert np.array_equal(getattr(c, a), getattr(d, a))
    with np.load(tmp_path / "c.npz") as z:
        parts = dict(z)
    np.savez_compressed(tmp_path / "bad.npz", **(parts | {"t5": np.float64(c.t5 * 1.01)}))
    with pytest.raises(ValueError, match="SHA-256"):
        ExpertCache.load(tmp_path / "bad.npz")
    with pytest.raises(ValueError, match="digests differ"):  # the file is not the one a prepared record names
        ExpertCache.load(tmp_path / "c.npz", expected_digest="0" * 64)
    assert ExpertCache.load(tmp_path / "c.npz", expected_digest=digest).digest == digest
    stale = parts | {"t5": np.float64(123.0), "fingerprint": np.str_("old")}
    old_digest = cache_mod._digest(parts["raters"], c.tau, 123.0, "old")
    np.savez_compressed(tmp_path / "old.npz", **(stale | {"sha256": np.str_(old_digest)}))
    with pytest.raises(R.StaleCacheError, match="fingerprint"):  # valid file of another reward fingerprint (B2-8)
        ExpertCache.load(tmp_path / "old.npz", expected_digest=old_digest)
    rebuilt = ExpertCache.load(tmp_path / "old.npz", expected_digest=old_digest, rebuild=True)
    assert rebuilt.t5 == c.t5 and rebuilt.digest == c.digest  # rebuilt from the verified masks
    monkeypatch.setattr(cache_mod, "fingerprint", lambda: "changed constants")
    with pytest.raises(R.StaleCacheError):
        ExpertCache.load(tmp_path / "c.npz")
    assert ExpertCache.load(tmp_path / "c.npz", rebuild=True).digest != c.digest  # under the new fingerprint


def test_fingerprint_covers_every_library_of_the_stored_values(monkeypatch):
    """tau uses numpy and OpenCV, T5 also scipy (quad, norm): a version change of any of them changes the
    fingerprint, so stored caches are refused as stale (B2-8)."""
    import cv2
    import scipy
    base = R.fingerprint()
    for mod in (np, cv2, scipy):
        with monkeypatch.context() as m:
            m.setattr(mod, "__version__", "0.0.0-test")
            assert R.fingerprint() != base, mod.__name__
    assert R.fingerprint() == base


def test_pickle_rederives_the_cache(three):
    _, c = three
    d = pickle.loads(pickle.dumps(c))
    assert d.digest == c.digest and np.array_equal(d.dist_in_core, c.dist_in_core) and d.same_sig == c.same_sig
    assert len(pickle.dumps(c)) < c.raters.nbytes + 2000  # only masks, tau, T5 travel


def test_cache_build_is_deterministic_across_interpreters(three):
    st_, c = three
    code = ("import numpy as np\nfrom objdel.reward import ExpertCache\n"
            "y, x = np.mgrid[:60, :60]\nd = np.hypot(y - 29.5, x - 29.5)\n"
            "print(ExpertCache.build(np.stack([d <= 18, d <= 19, d <= 20])).digest)")
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == c.digest == ExpertCache.build(st_).digest


# ------------------------------------------------------------------ per-object terms

def _random_masks(stack, rng, n):
    for _ in range(n):
        m = stack[int(rng.integers(len(stack)))]
        kind = int(rng.integers(0, 5))
        if kind == 0:
            yield H.jitter(m, float(rng.uniform(0.5, 3)), rng, 4.0)
        elif kind == 1:
            yield H.shift(m, int(rng.integers(-4, 5)), int(rng.integers(-4, 5)))
        elif kind == 2:
            yield H.grow(m, float(rng.uniform(-3, 3)))
        elif kind == 3:
            yield rng.random(m.shape) < 0.3
        else:
            yield np.zeros_like(m)


@pytest.mark.parametrize("seed", range(4))
def test_object_terms_match_reference(seed, three, mixed):
    rng = np.random.default_rng(seed)
    for st_, c in (three, mixed):
        rc = ref.cache(st_)
        for m in _random_masks(st_, rng, 25):
            t, r = score_mask(m, c), ref.object_terms(m, rc)
            assert t.sig == r["sig"] and t.topo_freq == r["topo_freq"]
            assert t.iou_tc == pytest.approx(r["iou_tc"], abs=1e-12)
            assert t.g_area == pytest.approx(r["g_area"], abs=1e-12)
            assert t.f5 == pytest.approx(r["f5"], rel=1e-5)


def test_iou_tc_compares_with_same_topology_raters_else_the_vote(mixed):
    st_, c = mixed
    t_ring = score_mask(st_[3], c)
    assert t_ring.iou_tc == 1.0 and t_ring.topo_freq == 0.25  # the minority hole is not judged by the vote
    assert score_mask(st_[0], c).iou_tc == 1.0 and score_mask(st_[0], c).topo_freq == 0.75
    two = st_[0] | disk(40, 2.5, 3, 3)  # a second part no rater drew: IoU vs the vote, not vs rater 0
    t = score_mask(two, c)
    assert t.topo_freq == 0 and t.iou_tc == pytest.approx(c.vote.sum() / two.sum())


def test_sub_filter_specks_are_invisible_to_topology(three):
    st_, c = three
    a = c.min_area
    assert a == max(5, math.ceil(0.005 * c.union.sum()))
    base = score_mask(st_[1], c)
    small = st_[1].copy()
    small[2, 2:4] = True  # 2-px speck
    hole = st_[1].copy()
    hole[30, 30] = hole[30, 31] = False  # 2-px pinhole
    for m in (small, hole):
        assert score_mask(m, c).sig == base.sig and score_mask(m, c).topo_freq == 1.0
    big = H.add_speck(st_[1], a, 6.0, np.random.default_rng(0))
    assert score_mask(big, c).topo_freq == 0.0


def test_sub_filter_raters_abstain_from_topology():
    """B2-4: a rater whose whole mask is below the size filter casts no topology vote: no histogram entry, no IoU_tc
    partner, shares over the voting raters; it stays in the band, vote and area range. All abstain: all vote."""
    big, tiny = disk(40, 6), np.zeros((40, 40), bool)
    tiny[19:21, 19:21] = True  # 4 px inside the 113 px disk; min_area = 5
    c = ExpertCache.build(np.stack([big, big, tiny]))
    assert c.rater_sigs == ((0,), (0,), ()) and c.voters == (0, 1) and c.same_sig == {(0,): [0, 1]}
    assert score_mask(big, c).topo_freq == 1.0 and min(c.rater_areas) == 4  # the area range keeps the tiny rater
    speck = score_mask(tiny, c)
    assert speck.sig == () and speck.topo_freq == 0.0 and speck.iou_tc == pytest.approx(4 / big.sum())  # vs the vote
    r = ref.cache(c.raters)
    for m in (big, tiny, disk(40, 7), np.zeros_like(big)):
        t, rt = score_mask(m, c), ref.object_terms(m, r)
        assert (t.sig, t.topo_freq, t.iou_tc) == (rt["sig"], rt["topo_freq"], pytest.approx(rt["iou_tc"], abs=1e-12))
    both = ExpertCache.build(np.stack([tiny, tiny]))
    assert both.voters == (0, 1) and score_mask(tiny, both).topo_freq == 1.0  # all abstain: the unanimous answer
    assert R.voters([(), ()]) == (0, 1) and R.voters([(1,), (), (0, 0)]) == (0, 2)


def test_rater_answers_are_joint_over_voters_and_check_rater_counts(riga_like):
    _, _, caches, _ = riga_like
    assert R.rater_answers(caches, ("disc", "cup")) == {((0,), (0,)): [0, 1], ((0,), (1,)): [2]}
    tiny = np.zeros((50, 50), bool)
    tiny[20:22, 20:22] = True
    cup = ExpertCache.build(np.stack([caches["cup"].raters[0], caches["cup"].raters[1], tiny]))
    assert R.rater_answers({"disc": caches["disc"], "cup": cup}, ("disc", "cup")) == {((0,), (0,)): [0, 1]}
    two = ExpertCache.build(caches["cup"].raters[:2])
    for pair in ({"disc": caches["disc"], "cup": two}, {"disc": two, "cup": caches["cup"]}):  # B2-9
        with pytest.raises(ValueError, match="rater counts"):
            R.rater_answers(pair, ("disc", "cup"))


def test_area_guard(three):
    st_, c = three
    a_min, a_max = min(c.rater_areas), max(c.rater_areas)
    assert score_mask(c.union, c).g_area == pytest.approx(a_max / c.union.sum())
    assert score_mask(c.core, c).g_area == pytest.approx(c.core.sum() / a_min)
    first = lambda n: (np.arange(c.union.size) < n).reshape(c.union.shape)  # noqa: E731  (n pixels)
    assert score_mask(first(2 * a_max), c).g_area == 0.5
    assert score_mask(first(a_min // 2), c).g_area == pytest.approx((a_min // 2) / a_min)
    assert score_mask(np.zeros_like(c.union), c).g_area == 0.0


def test_band_is_monotone_in_excursion(three):
    _, c = three
    for sign in (1, -1):
        f = [score_mask(H.grow(c.vote, sign * d * c.tau), c).f5 for d in (0, 0.5, 1, 1.5, 2, 3)]
        assert all(b <= a for a, b in zip(f, f[1:], strict=False)) and f[-1] < f[0]


def test_shape_mismatch_raises(three):
    with pytest.raises(ValueError):
        score_mask(np.zeros((10, 10), bool), three[1])


# ------------------------------------------------------------------ completion text

G20 = 20  # grid units per pixel on a 50 x 50 crop


def rect(x0, y0, x1, y1, hole=False):
    r = [[x0 * G20, y0 * G20], [x1 * G20, y0 * G20], [x1 * G20, y1 * G20], [x0 * G20, y1 * G20]]
    return r[::-1] if hole else r


def text(objs):
    return json.dumps({n: [{"outer": o, "holes": hs} for o, hs in parts] for n, parts in objs.items()},
                      separators=(",", ":"))


def raster(parts):
    return to_mask([(np.array(o), [np.array(h) for h in hs]) for o, hs in parts], 50, 50).mask


@pytest.fixture(scope="module")
def riga_like():
    """Two objects, three raters; rater 2 draws a cup with a hole (a joint topology of its own)."""
    discs = [[(rect(10, 10, 40, 40), [])], [(rect(11, 10, 40, 41), [])], [(rect(10, 11, 41, 40), [])]]
    cups = [[(rect(18, 18, 32, 32), [])], [(rect(18, 19, 32, 32), [])],
            [(rect(17, 18, 32, 33), [rect(23, 23, 27, 27, True)])]]
    caches = {"disc": ExpertCache.build(np.stack([raster(p) for p in discs])),
              "cup": ExpertCache.build(np.stack([raster(p) for p in cups]))}
    rcaches = {n: ref.cache(c.raters) for n, c in caches.items()}
    return discs, cups, caches, rcaches


def test_contract_keys_and_failures(riga_like):
    discs, cups, caches, _ = riga_like
    names = ("disc", "cup")
    ok = score_completion(text({"disc": discs[0], "cup": cups[0]}), False, names, caches, 50, 50)
    assert set(ok) == set(R.NUMERIC) | {"error", "defects"}
    assert ok["format_ok"] == ok["legal"] == 1.0 and ok["error"] == "" and ok["defects"] == {}
    (row, objs), = R.score_many([(text({"disc": discs[0], "cup": cups[0]}), False, names, caches, 50, 50)])
    assert row == ok  # one pass: terms and ObjectTerms
    assert objs == {n: score_mask(raster(p[0]), caches[n]) for n, p in (("disc", discs), ("cup", cups))}
    assert R.score_many([("hello", False, names, caches, 50, 50)])[0][1] == {}
    for bad, trunc, code in ((text({"disc": discs[0], "cup": cups[0]}), True, "truncated"),
                             ('{"disc":[],"cup":[', False, "incomplete"), ("hello", False, "syntax"),
                             (text({"cup": cups[0], "disc": discs[0]}), False, "keys")):
        r = score_completion(bad, trunc, names, caches, 50, 50)
        assert r["error"] == code and r["format_ok"] == r["legal"] == 0.0 and r["sig_id"] == -1.0
        assert all(r[k] == 0.0 for k in R.NUMERIC if k != "sig_id")


def test_illegal_completion_keeps_terms_but_is_not_legal(riga_like):
    discs, cups, caches, _ = riga_like
    bowtie = [[200, 200], [800, 800], [800, 200], [200, 800]]
    r = score_completion(text({"disc": [(bowtie, [])], "cup": cups[0]}), False, ("disc", "cup"), caches, 50, 50)
    assert r["format_ok"] == 1.0 and r["legal"] == 0.0 and r["defects"] and r["error"] == ""
    a = advantages_v4({k: np.repeat(r[k], 8) if k != "legal" else np.r_[0.0, np.ones(7)]
                       for k in R.NUMERIC} | {"iou_tc": np.r_[0.9, np.full(7, 0.5)]})
    b = advantages_v4({k: np.repeat(r[k], 8) if k != "legal" else np.r_[0.0, np.ones(7)]
                       for k in R.NUMERIC} | {"iou_tc": np.r_[0.1, np.full(7, 0.5)]})
    assert np.array_equal(a, b) and a[0] < a[1:].min()  # an illegal rollout's terms do not count


def test_empty_object_list_is_a_legal_empty_mask(riga_like):
    _, cups, caches, _ = riga_like
    r = score_completion(text({"disc": [], "cup": cups[0]}), False, ("disc", "cup"), caches, 50, 50)
    assert r["legal"] == 1.0 and r["topo_any"] == 0.0  # no rater drew "no disc"
    assert r["g_area"] == pytest.approx(0.5) and r["iou_tc"] == pytest.approx(0.5)  # disc: 0, cup: 1


def test_iou_min_gates_topology_per_object(riga_like):
    """B2-7: a perfect disc and a cup outside it average IoU_tc 0.5 but the cup's IoU_tc is 0: not a valid answer."""
    discs, cups, caches, _ = riga_like
    names = ("disc", "cup")
    missed = score_completion(text({"disc": discs[0], "cup": [(rect(0, 0, 9, 9), [])]}), False, names, caches, 50, 50)
    assert missed["topo_any"] == 1.0 and missed["iou_min"] == 0.0 and missed["iou_tc"] >= 0.5
    good = score_completion(text({"disc": discs[0], "cup": cups[0]}), False, names, caches, 50, 50)
    assert good["iou_min"] == 1.0
    T = {k: np.array([good[k]] * 7 + [missed[k]]) for k in R.NUMERIC}
    assert topology_reward(T)[7] == 0.0 and topology_reward(T)[0] > 0
    single = {k: v for k, v in T.items() if k != "iou_min"}  # terms without iou_min: gated on iou_tc (one object)
    assert topology_reward(single)[7] > 0


def test_joint_topology_over_objects(riga_like):
    discs, cups, caches, _ = riga_like
    names = ("disc", "cup")
    r_major = score_completion(text({"disc": discs[0], "cup": cups[1]}), False, names, caches, 50, 50)
    r_minor = score_completion(text({"disc": discs[1], "cup": cups[2]}), False, names, caches, 50, 50)
    assert r_major["topo_freq"] == pytest.approx(2 / 3) and r_minor["topo_freq"] == pytest.approx(1 / 3)
    assert r_major["sig_id"] != r_minor["sig_id"] and r_major["topo_any"] == r_minor["topo_any"] == 1.0


@pytest.mark.parametrize("seed", range(3))
def test_completion_terms_match_reference(seed, riga_like):
    discs, cups, caches, rcaches = riga_like
    rng = np.random.default_rng(seed)
    names = ("disc", "cup")
    for _ in range(10):
        objs = {}
        for n, pool in (("disc", discs), ("cup", cups)):
            (o, hs), = pool[int(rng.integers(3))]
            d = rng.integers(-40, 41, 2)
            objs[n] = [((np.array(o) + d).clip(0, 1000).tolist(), [(np.array(h) + d).tolist() for h in hs])]
        t = text(objs)
        a, b = score_completion(t, False, names, caches, 50, 50), ref.completion_terms(t, False, names, rcaches, 50, 50)
        for k in ("format_ok", "legal", "topo_any", "topo_freq", "iou_tc", "iou_min", "g_area"):
            assert a[k] == pytest.approx(b[k], abs=1e-12), k
        assert a["f5"] == pytest.approx(b["f5"], rel=1e-5)


def test_score_batch_matches_jobs_in_order_and_concurrent_calls(riga_like):
    discs, cups, caches, _ = riga_like
    names = ("disc", "cup")
    jobs = [(text({"disc": discs[i % 3], "cup": cups[(i * 7) % 3]}), i % 5 == 0, names, caches, 50, 50)
            for i in range(24)]
    serial = score_batch(jobs)
    assert [serial["topo_freq"][i] for i in range(24)] == [score_completion(*j)["topo_freq"] for j in jobs]
    with ThreadPoolExecutor(4) as pool:  # concurrent calls share the checker
        for t in pool.map(score_batch, [jobs] * 4):
            assert all(np.array_equal(serial[k], t[k]) for k in R.NUMERIC) and serial["error"] == t["error"]


# ------------------------------------------------------------------ advantages

FREQS = sorted({j / k for k in range(2, 7) for j in range(1, k + 1)})


@st.composite
def term_batches(draw):
    B = draw(st.integers(1, 6))
    n = B * 8
    fl = st.floats(0, 1, allow_nan=False)
    f = draw(st.lists(st.sampled_from(FREQS), min_size=n, max_size=n))
    sid = draw(st.lists(st.integers(0, 3), min_size=n, max_size=n))
    return {"iou_tc": np.array(draw(st.lists(fl, min_size=n, max_size=n))),
            "f5": np.array(draw(st.lists(fl, min_size=n, max_size=n))),
            "g_area": np.array(draw(st.lists(fl, min_size=n, max_size=n))),
            "legal": np.array(draw(st.lists(st.sampled_from([0.0, 1.0, 1.0, 1.0]), min_size=n, max_size=n))),
            "topo_any": np.array(draw(st.lists(st.sampled_from([0.0, 1.0, 1.0]), min_size=n, max_size=n))),
            "topo_freq": np.array(f), "sig_id": np.array(sid, float)}


@settings(max_examples=300, deadline=None)
@given(term_batches())
def test_advantages_match_reference(terms):
    terms["topo_freq"] = np.where(terms["topo_any"] > 0, terms["topo_freq"], 0.0)
    assert np.allclose(advantages_v4(terms), ref.advantages(terms), atol=1e-9, rtol=0)
    terms["iou_min"] = np.minimum(terms["iou_tc"], terms["f5"])  # a second object's IoU_tc (any value <= the mean)
    assert np.allclose(advantages_v4(terms), ref.advantages(terms), atol=1e-9, rtol=0)


def test_low_spread_groups_are_dropped():
    base = {"iou_tc": np.full(8, 0.9), "f5": np.full(8, 0.95), "g_area": np.ones(8), "legal": np.ones(8),
            "topo_any": np.ones(8), "topo_freq": np.full(8, 0.5), "sig_id": np.zeros(8)}
    flat = base | {"iou_tc": 0.9 + np.linspace(0, 0.009, 8)}
    assert (advantages_v4(flat) == 0).all()
    assert (advantages_v4(base | {"iou_tc": 0.9 + np.linspace(0, 0.011, 8)}) != 0).any()
    assert (advantages_v4(flat | {"sig_id": np.r_[1.0, np.zeros(7)]}) != 0).any()  # two topology answers: kept


def test_groups_are_independent_and_centred():
    rng = np.random.default_rng(0)
    terms = {k: rng.random(64) for k in ("iou_tc", "f5", "g_area")} | {
        "legal": (rng.random(64) > 0.1).astype(float), "topo_any": np.ones(64),
        "topo_freq": rng.choice(FREQS, 64), "sig_id": rng.integers(0, 3, 64).astype(float)}
    a = advantages_v4(terms)
    first = advantages_v4({k: v[:8] for k, v in terms.items()})
    assert np.array_equal(a[:8], first)  # no batch-level z-score
    assert np.abs(a.reshape(-1, 8).sum(1)).max() < 1e-12
    with pytest.raises(ValueError):
        advantages_v4({k: v[:12] for k, v in terms.items()})


@pytest.mark.parametrize("freqs, lam", [((0.75, 0.25), 0.5), ((0.5, 0.5), 0.5), ((0.5, 0.25, 0.25), 0.5),
                                        ((2 / 3, 1 / 3), 0.5)])
def test_debiased_topology_reward_favours_no_rater_answer_at_the_target(freqs, lam):
    """Exact enumeration over all groups of 8 drawn at pi = f: E[r_topo | answer] is equal for every rater
    answer (where the clip does not bind)."""
    from math import prod
    G, K = 8, len(freqs)
    exp_r = np.zeros(K)
    for combo in product(range(K), repeat=G):
        p = prod(freqs[c] for c in combo)
        terms = {"sig_id": np.array(combo, float), "topo_freq": np.array([freqs[c] for c in combo]),
                 **{k: np.ones(G) for k in ("topo_any", "iou_tc", "legal", "f5", "g_area")}}
        r = topology_reward(terms)
        for c in range(K):
            exp_r[c] += p * r[np.array(combo) == c].sum()
    assert np.ptp(exp_r / (G * np.array(freqs))) < 1e-9


def test_loo_bias_values():
    assert loo_bias(1.0, 8) == 0.0
    for f in FREQS:
        assert loo_bias(f, 8) == pytest.approx(ref.loo_bias(f, 8), abs=1e-12)
    assert [round(loo_bias(f, 8), 3) for f in (0.25, 0.5, 0.75)] == [0.145, 0.088, 0.028]  # lab ITERATION3


def test_valid_beats_invalid_at_every_reachable_state():
    """Every rater share f = j/k (k = 2..6) and every count n of the valid answer in the group: the valid
    rollouts' topology reward and advantage are strictly above the invalid ones' (geometry held equal)."""
    for f in FREQS:
        for n in range(1, 9):
            terms = {"sig_id": np.r_[np.zeros(n), np.arange(1, 9 - n)],
                     "topo_freq": np.r_[np.full(n, f), np.zeros(8 - n)],
                     "topo_any": np.r_[np.ones(n), np.zeros(8 - n)], "iou_tc": np.full(8, 0.8), "f5": np.full(8, 0.9),
                     "g_area": np.ones(8), "legal": np.ones(8)}
            r, a = topology_reward(terms), advantages_v4(terms)
            assert r[:n].min() >= REWARD_V4.eps
            if n < 8:
                assert r[:n].min() > r[n:].max() == 0.0 and a[:n].min() > a[n:].max()


@pytest.mark.parametrize("term, bound", [("iou_tc", 0.006), ("f5", 0.006), ("g_area", 0.006), ("legal", 0.0004)])
def test_dark_room_noise_on_a_constant_term_moves_advantages_little(term, bound):
    rng = np.random.default_rng(1)
    for _ in range(200):
        terms = {"iou_tc": rng.uniform(0.6, 0.95, 8), "f5": rng.uniform(0.6, 1, 8), "g_area": rng.uniform(0.8, 1, 8),
                 "legal": np.ones(8), "topo_any": np.ones(8), "topo_freq": np.full(8, 0.5),
                 "sig_id": rng.integers(0, 2, 8).astype(float)}
        terms[term] = np.ones(8)
        a0 = advantages_v4(terms)
        a1 = advantages_v4(terms | {term: 1 - rng.uniform(0, 1e-3, 8)})
        assert np.abs(a1 - a0).max() <= bound * (1 + 1e-9)


# ------------------------------------------------------------------ determinism (K10)

K10_SNIPPET = """
import hashlib, json
import numpy as np
from objdel.reward import NUMERIC, ExpertCache, advantages_v4, score_batch
def digest(t):
    h = hashlib.sha256()
    for k in NUMERIC:
        h.update(np.ascontiguousarray(t[k], np.float64).tobytes())
    h.update(np.ascontiguousarray(advantages_v4(t), np.float64).tobytes())
    return h.hexdigest()
G = 20
def rect(x0, y0, x1, y1):
    return [[x0 * G, y0 * G], [x1 * G, y0 * G], [x1 * G, y1 * G], [x0 * G, y1 * G]]
def raster(x0, y0, x1, y1):
    m = np.zeros((50, 50), bool); m[y0:y1, x0:x1] = True; return m
stack = np.stack([raster(10, 10, 30, 30), raster(11, 10, 31, 31), raster(9, 12, 30, 29)])
caches = {"nodule": ExpertCache.build(stack)}
rng = np.random.default_rng(5)
jobs = []
for i in range(64):
    x0, y0 = map(int, rng.integers(5, 15, 2)); x1, y1 = map(int, rng.integers(25, 40, 2))
    jobs.append((json.dumps({"nodule": [{"outer": rect(x0, y0, x1, y1), "holes": []}]}, separators=(",", ":")),
                 bool(i % 13 == 0), ("nodule",), caches, 50, 50))
"""


def _k10_run(mode):
    ns = {}
    exec(K10_SNIPPET, ns)
    jobs = ns["jobs"]
    if mode == "serial":
        t = score_batch(jobs)
    elif mode == "threads":
        with ThreadPoolExecutor(8) as pool:
            t = list(pool.map(score_batch, [jobs] * 8))[-1]
    else:
        with ProcessPoolExecutor(2, mp_context=multiprocessing.get_context("spawn")) as pool:
            t = pool.submit(score_batch, jobs).result()
    return ns["digest"](t)


def test_k10_identical_rewards_across_runs_threads_processes_interpreters():
    ref_digest = _k10_run("serial")
    assert _k10_run("serial") == _k10_run("threads") == _k10_run("processes") == ref_digest
    code = K10_SNIPPET + "print(digest(score_batch(jobs)))"
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, check=True,
                         env={"PYTHONHASHSEED": "123"})
    assert out.stdout.strip() == ref_digest


def test_legal_gate_is_the_one_legality_rule():
    """REVIEW #2: the legality gate shared by advantages_v4, topology_valid and eval: an illegal rollout's terms do not
    move its advantage (they are zeroed), and it is never topology-valid."""
    rng = np.random.default_rng(0)
    T = {k: rng.uniform(0.6, 1.0, 8) for k in ("iou_tc", "iou_min", "f5", "g_area", "topo_freq")}
    T |= {"legal": np.array([1.0] * 6 + [0.0, 0.0]), "topo_any": np.ones(8), "sig_id": np.zeros(8)}
    assert R.legal_gate(T).tolist() == [True] * 6 + [False] * 2
    assert not R.topology_valid(T)[6:].any() and R.topology_valid(T)[:6].all()
    U = {k: v.copy() for k, v in T.items()}
    for k in ("iou_tc", "f5", "g_area"):
        U[k][6:] = rng.uniform(0, 1, 2)  # illegal rollouts' terms are not seen by the advantage
    assert np.array_equal(R.advantages_v4(T), R.advantages_v4(U))
