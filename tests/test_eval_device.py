"""Evaluation pixel work on the device (Phase D, D7; tests U12, U13, I03, P03 (eval_sets part), D-R05 of
Implementation.md D.9, B6 and B14 of D.spec.md).

`gc.eval_sets` returns integer counts per set of K sample polygons against one cache's tables; the oracle is numpy on
the rasters (tests/oracles/eval_oracle.py: the host code metrics() ran before D7). `metrics()` builds every field from
the counts and equals its pre-D7 values (tests/goldens/d_r05_metrics.json) exactly, `nsd` within 1e-12, without
`objdel.geometry.to_masks` or any CPU labelling or EDT.
"""
import hashlib
import importlib
import json
import math
import os
from pathlib import Path

import cv2
import numpy as np
import pytest
from oracles import eval_oracle as O
from oracles.eval_items import BOWTIE, items, val_items
from scipy import ndimage

from objdel import _gpucheck as gc
from objdel.eval import device as D
from objdel.geometry import GRID, parse
from objdel.geometry.shape import pack_objects, to_masks
from objdel.reward import ExpertCache

M_ = importlib.import_module("objdel.eval.metrics")
DATA = Path(os.environ.get("OBJDEL_DATA_ROOT", Path(__file__).resolve().parents[1] / "data"))
GOLDEN = Path(__file__).parent / "goldens" / "d_r05_metrics.json"
RUNS = Path(__file__).resolve().parents[1] / "runs"


# ---------------------------------------------------------------- U12 host cuts

def test_soft_cut_brute_force():
    for n in range(1, 41):
        for t in (0.0, 1e-12, 0.1, 0.3, 0.5, 1 / 3, 2 / 3, 0.7, 0.9, 1.0):
            assert D.soft_cut(n, t) == min(c for c in range(n + 1) if c / n >= t), (n, t)
    assert (D.soft_cut(3, 1 / 3), D.soft_cut(10, 0.7), D.soft_cut(8, 0.0), D.soft_cut(8, 1.0)) == (1, 7, 0, 8)


def test_nsd_cut_brute_force():
    for tau in (0.0, 0.5, 1.0, math.sqrt(2), np.nextafter(math.sqrt(2), 0), 1.5, 2.0, 2.9999999999999996, 3.0,
                3.0000000000000004, 100.7):
        assert D.nsd_cut(float(tau)) == max(q for q in range(10 ** 5) if math.sqrt(q) <= tau), tau
    assert (D.nsd_cut(math.sqrt(2)), D.nsd_cut(float(np.nextafter(math.sqrt(2), 0))), D.nsd_cut(3.0)) == (2, 1, 9)
    assert D.nsd_cut(1e6) == 10 ** 12
    for bad in (-1.0, -1e-300, math.nan, math.inf, 2.0 ** 27):
        with pytest.raises(ValueError):
            D.nsd_cut(bad)


# ---------------------------------------------------------------- U13 duplicates

def _counts(raw_by_object):
    """ItemCounts stand-ins with the raw intersections of {object: bool (K, H, W)}."""
    out = {}
    for n, raw in raw_by_object.items():
        r = raw.reshape(len(raw), -1).astype(np.int64)
        inter = r @ r.T
        out[n] = D.ItemCounts(None, inter, np.diag(inter).copy(), None, None, None, None, 0.0)
    return out


def test_duplicate_rule():
    """U13: flags and raw rasters of every object; a subset is not equal; empty equals empty; one differing object
    breaks it; equals the old raw-hash classes."""
    z = np.zeros((6, 6), bool)
    a, b = z.copy(), z.copy()
    a[1:4, 1:4] = True
    b[1:4, 1:5] = True  # a is a subset of b
    one = z.copy()
    one[1:4, 1:4], one[5, 5] = True, True  # a plus one pixel
    cases = [  # (format_ok, legal, {object: K raw masks}, expected keys)
        ([1, 1, 1], [1, 1, 1], {"o": np.stack([a, a, b])}, [0, 0, 2]),
        ([1, 1, 1], [1, 0, 1], {"o": np.stack([a, a, a])}, [0, 1, 0]),
        ([1, 0, 1], [0, 0, 0], {"o": np.stack([a, a, a])}, [0, 1, 0]),
        ([1, 1, 1], [1, 1, 1], {"o": np.stack([a, one, a])}, [0, 1, 0]),
        ([1, 1, 1, 1], [1, 1, 1, 1], {"o": np.stack([a, b, b, a])}, [0, 1, 1, 0]),
        ([0, 0, 1], [0, 0, 1], {"o": np.stack([z, z, z])}, [0, 0, 2]),
        ([1, 1, 1], [1, 1, 1], {"disc": np.stack([b, b, b]), "cup": np.stack([a, one, a])}, [0, 1, 0]),
    ]
    for fmt, legal, raw, want in cases:
        assert D.duplicate_keys(fmt, legal, _counts(raw)) == want
        hashes = [O.dup_hash(fmt[i], legal[i], {n: m[i] for n, m in raw.items()}) for i in range(len(fmt))]
        assert [hashes.index(h) for h in hashes] == want


# ---------------------------------------------------------------- I03 eval_sets against numpy

def _cache(rng, k, H, W):
    y, x = np.mgrid[:H, :W]
    cy, cx = H / 2 + rng.uniform(-2, 2), W / 2 + rng.uniform(-2, 2)
    st = np.stack([np.hypot(y - cy - rng.uniform(-2, 2), x - cx - rng.uniform(-2, 2)) <= rng.uniform(3, min(H, W) / 3)
                   for _ in range(k)])
    if rng.random() < 0.3:
        st[0] &= ~(np.hypot(y - cy, x - cx) <= 2)
    return ExpertCache(st, 1.0, 1.0)


def _polys(rng, K, W, H, kind):
    """K random objects: rectangles, two-part shapes, an empty list and a bow-tie."""
    out = []
    for i in range(K):
        if kind == "empty" or i % 7 == 3:
            out.append([])
        elif kind == "same" and out:
            out.append(out[0])
        elif i % 7 == 5:
            out.append(parse(BOWTIE["nodule"], ("nodule",)).objs["nodule"])
        else:
            parts = []
            for _ in range(1 + (i % 3 == 2)):
                x0, y0 = rng.integers(0, 600, 2)
                x1, y1 = x0 + rng.integers(50, 400), y0 + rng.integers(50, 400)
                parts.append((np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]]), []))
            out.append(parts)
    return [(p, W, H) for p in out]


def _call(sets, K, T=(0.1, 0.5, 0.9), tau=2.0, d=2, cap=D.MAX_SCRATCH_BYTES):
    objs = [o for s in sets for o in s["objs"]]
    k = sets[0]["cache"].k
    cuts = np.array([[[D.soft_cut(K, t), D.soft_cut(k, t)] for t in T]] * len(sets), np.int32)
    return gc.eval_sets(*pack_objects(objs), GRID, np.array([s["cache"].device() for s in sets], np.int64), K,
                        np.concatenate([s["legal"] for s in sets]).astype(np.uint8), cuts,
                        np.full(len(sets), D.nsd_cut(tau), np.int64), np.full(len(sets), d, np.int32),
                        max_scratch_bytes=cap)


def _sets(seed, K, k, n=6):
    rng = np.random.default_rng(seed)
    sets = []
    for i in range(n):
        H, W = int(rng.integers(16, 120)), int(rng.integers(16, 120))
        if i == 0:
            H = W = 4 * int(rng.integers(4, 9))  # a word-aligned frame (W % 32 == 0 for 32)
        kind = "empty" if i == 1 else "same" if i == 2 else "random"
        objs = _polys(rng, K, W, H, kind)
        sets.append({"objs": objs, "cache": _cache(rng, k, H, W), "legal": rng.random(K) < 0.8})
    return sets


def _expect(s, K, T=(0.1, 0.5, 0.9), tau=2.0, d=2):
    c = s["cache"]
    raw = np.stack([g.mask for g in to_masks(s["objs"])])
    drawn = raw & s["legal"][:, None, None]
    return O.counts(raw, drawn, c.raters, c.vote, T, tau, d) | {
        "vote_iou_tc": O.vote_iou_tc(raw, s["legal"], c)}


@pytest.mark.parametrize("K, k", [(1, 2), (8, 4), (16, 32), (4, 3)])
def test_eval_sets_equals_numpy(K, k):
    """I03: 6 sets (an all-empty set, an all-identical set, random rectangles with illegal and empty samples):
    every output equals the numpy oracle on the rasters; inter symmetric with the areas on its diagonal."""
    sets = _sets(K * 100 + k, K, k)
    r = _call(sets, K)
    for i, s in enumerate(sets):
        want = _expect(s, K)
        for key, v in want.items():
            assert np.array_equal(r[key][i], v), (i, key)
        assert np.array_equal(r["inter"][i], r["inter"][i].T)
        assert np.array_equal(r["area_raw"][i], np.diag(r["inter_raw"][i]))


def test_eval_sets_literal_square():
    """I03 literal: a 2 x 2 square at (1, 1) on a 4 x 4 frame has 4 corner surfels of sqrt(2)/2 and one of length 1 at
    the middle of each side; a rater equal to it gives no surfel beyond any cut; all samples illegal draw nothing."""
    st = np.zeros((2, 4, 4), bool)
    st[:, 1:3, 1:3] = True
    c = ExpertCache(st, 1.0, 1.0)
    text = '{"nodule":[{"outer":[[250,250],[750,250],[750,750],[250,750]],"holes":[]}]}'
    objs = [(parse(text, ("nodule",)).objs["nodule"], 4, 4)] * 2
    s = {"objs": objs, "cache": c, "legal": np.array([True, True])}
    r = _call([s], 2, tau=0.0)
    assert r["surf"][0].tolist() == [[4, 4, 0]] * 4 and not r["nsd_out"].any()
    assert r["inter"][0].tolist() == [[4] * 5] * 5
    s["legal"] = np.array([False, False])
    r = _call([s], 2, tau=0.0)
    assert r["surf"][0][:2].tolist() == [[0, 0, 0]] * 2 and r["soft"][0][:, 0].tolist() == [0, 0, 0]
    assert r["nsd_out"][0, :, :, 1].tolist() == [[[4, 4, 0]] * 2] * 2  # every rater surfel: the sample has no surface
    assert r["inter_raw"][0].tolist() == [[4, 4], [4, 4]] and r["vote_iou_tc"][0] == 0.0


def test_eval_sets_scratch_cap_and_order_invariance():
    """P03 (eval_sets part, B6): caps 1 << 12, 1 << 16, 1 << 28 and the default give byte-identical outputs, and
    eval_scratch_peak after reset_stats is > 0 and <= the cap; a cap below one set raises ValueError."""
    sets = _sets(5, 8, 3, n=8)
    want = _call(sets, 8)
    peak_small = None
    for cap in (1 << 16, 1 << 18, 1 << 28):
        gc.reset_stats()
        try:
            got = _call(sets, 8, cap=cap)
        except ValueError:
            assert cap < 1 << 20
            continue
        peak = gc.stats()["eval_scratch_peak"]
        assert 0 < peak <= cap
        peak_small = peak_small or peak
        assert all(got[k].tobytes() == want[k].tobytes() for k in want), cap
    assert peak_small is not None
    with pytest.raises(ValueError, match="scratch"):
        _call(sets, 8, cap=1 << 12)
    rev = _call(sets[::-1], 8)
    assert all(rev[k][::-1].tobytes() == want[k].tobytes() for k in want)


def test_eval_sets_bad_inputs():
    sets = _sets(9, 2, 2, n=1)
    c = sets[0]["cache"]
    ok = dict(K=2)
    with pytest.raises(ValueError):
        _call(sets, 3)  # shapes are not K per set
    other = ExpertCache(np.zeros((3, 5, 5), bool), 1.0, 1.0)
    with pytest.raises(ValueError):
        _call(sets + [{"objs": _polys(np.random.default_rng(0), 2, 5, 5, "random"), "cache": other,
                       "legal": np.array([True, True])}], **ok)  # two rater counts in one call
    with pytest.raises(ValueError):
        _call(sets, 2, d=0)
    h = int(gc.add_tables_from_raters([c.raters.view(np.uint8)], np.array([1.0]), np.array([c.min_area]))["handles"][0])
    gc.drop_tables(h)
    gone = ExpertCache(c.raters, 1.0, 1.0)
    gone._device = h  # a stale handle (no finalizer: device() never built it)
    with pytest.raises(ValueError):
        _call([{**sets[0], "cache": gone}], 2)  # a dropped handle


# ---------------------------------------------------------------- D-R05 metrics() equality and S01

def _no_host_geometry(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("host pixels in metrics()")
    monkeypatch.setattr(M_, "to_masks", boom, raising=False)
    from objdel import geometry
    from objdel.geometry import shape
    for mod, name in ((geometry, "to_masks"), (shape, "to_masks"), (ndimage, "distance_transform_edt"),
                      (ndimage, "label"), (cv2, "connectedComponents"), (cv2, "connectedComponentsWithStats"),
                      (cv2, "connectedComponentsWithAlgorithm"), (cv2, "erode")):
        monkeypatch.setattr(mod, name, boom)


def _equal(got, want, where):
    """Every field equal; nsd within 1e-12."""
    assert got.keys() == want.keys(), where
    for key, v in want.items():
        if isinstance(v, dict):
            _equal(got[key], v, f"{where}.{key}")
        elif key == "nsd":
            assert np.allclose(got[key], v, rtol=0, atol=1e-12), where
        else:
            assert got[key] == v, f"{where}.{key}"


def test_metrics_equal_the_golden_without_host_pixels(monkeypatch):
    """D-R05 / S01: the 50 synthetic items (duplicates, illegal, unparseable, truncated, all-identical) give every
    field of the pre-D7 golden, with to_masks, cv2.erode, scipy labelling and EDT and OpenCV labelling raising."""
    g = json.loads(GOLDEN.read_text())
    its = items()
    h = hashlib.sha256()
    for texts, tr, raters, W, H, names in its:
        h.update(repr((texts, tr, W, H, names)).encode())
        for n in names:
            h.update(np.packbits(raters[n]).tobytes())
    assert h.hexdigest() == g["inputs_sha256"]
    caches = [{n: ExpertCache.build(raters[n]) for n in names} for _, _, raters, _, _, names in its]
    _no_host_geometry(monkeypatch)
    for i, ((texts, tr, raters, W, H, names), want) in enumerate(zip(its, g["metrics"], strict=True)):
        _equal(M_.metrics(texts, raters, W, H, names, truncated=tr, caches=caches[i]), want, f"item {i}")
    assert g["metrics"][1]["groups"]["dup"] == [0.875] and g["metrics"][1]["samples"]["legal"] == [1.0] * 8


@pytest.mark.data
def test_metrics_equal_the_golden_on_real_val_items(monkeypatch):
    """D-R05 data half: 20 LIDC and 20 RIGA val items (answers cycled from their targets, one illegal, one
    unparseable) equal runs/d7/r05_val_golden.json, written by the code before the switch."""
    path = RUNS / "d7" / "r05_val_golden.json"
    if not path.exists():
        pytest.skip(f"{path} missing (written once by the pre-D7 code)")
    want = json.loads(path.read_text())
    got = {}
    for ds in ("lidc", "riga"):
        for item_id, texts, tr, W, H, names, caches in val_items(DATA / "prepared", ds):
            got[item_id] = M_.metrics(texts, {}, W, H, names, tr, caches=caches)
    assert got.keys() == want.keys() and len(got) == 40
    for key in want:
        _equal(got[key], want[key], key)


# ---------------------------------------------------------------- I18 one device process (D7 step 6)

def _child_device(raters) -> str:
    try:
        ExpertCache(raters, 1.0, 1.0).device()
        return "built"
    except RuntimeError as e:
        return str(e)


def test_device_in_a_forked_child_raises_and_a_fresh_process_builds():
    """I18: a forked child inherits the pid recorded at import, so device() raises there and leaves the parent's
    tables alone; a spawned (fresh) process records its own pid and builds; the parent still scores."""
    import multiprocessing as mp

    st = np.zeros((2, 8, 8), bool)
    st[:, 2:6, 2:6] = True
    parent = ExpertCache(st, 1.0, 1.0)
    h = parent.device()
    with mp.get_context("fork").Pool(1) as p:
        assert p.apply(_child_device, (st,)) == "device() in a forked child"
    with mp.get_context("spawn").Pool(1) as p:
        assert p.apply(_child_device, (st,)) == "built"
    assert parent.device() == h and parent.perimeter > 0


def test_eval_scripts_score_in_one_process():
    """I18: final_eval, select_ckpt, g_probe and the ceiling CLI fork no pool; select_ckpt and g_probe have no
    --procs (the preflight dry run passes none)."""
    root = Path(__file__).resolve().parents[1]
    for name in ("final_eval", "select_ckpt", "g_probe", "ceiling"):
        src = (root / "experiments" / "phase2" / f"{name}.py").read_text()
        assert "Pool" not in src and '"--procs"' not in src, name
    assert '"--procs"' not in (root / "tests" / "test_p2_scripts.py").read_text()
