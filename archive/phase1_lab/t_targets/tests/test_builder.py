"""Property tests of the SFT target builder (pytest + hypothesis).

Every built target must be: legal under rule L (checked with the exact pure-Python oracle),
canonical, round-trip identical through JSON, topology-preserving on re-rasterization with the
oracle's rule, IoU >= iou_min (or an explicit exact-ring fallback that is lossless), deterministic,
and - for the lossless configuration - idempotent: building from the re-rasterized mask returns
the same text.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import oracle as O
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from objdel.checker import betti as betti_ref
from targets.builder import (LIDC_SCHEMA, RIGA_DERIVED, RIGA_SCHEMA, BuildConfig, build_target, is_canonical,
                             parse, rasterize_parts)
from targets.crops import apply_crop

CONFIGS = {
    "hyb3": BuildConfig("minseg_hybrid", {"tau": 3}, iou_min=0.98),
    "lossless": BuildConfig("minseg_exact", {}, iou_min=1.0),
    "dp5": BuildConfig("dp", {"eps": 5}, iou_min=0.98),
    "potrace": BuildConfig("minseg_linf", {"t_px": 0.5}, iou_min=0.98),
}


@st.composite
def masks(draw, max_side=18):
    H = draw(st.integers(1, max_side))
    W = draw(st.integers(1, max_side))
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
        elif kind == "diag":  # 8-connected diagonal line: pinches
            for k in range(min(r1 - r0, c1 - c0) + 1):
                m[r0 + k, c0 + k] = True
        elif kind == "checker":
            sub = (np.add.outer(np.arange(r1 - r0 + 1), np.arange(c1 - c0 + 1)) % 2) == 0
            m[r0:r1 + 1, c0:c1 + 1] |= sub
        elif kind == "hline":
            m[r0, c0:c1 + 1] = True
        else:
            m[r0:r1 + 1, c0] = True
    return m


def oracle_legal(parts, G=1000):
    rings = [r for o, hs in parts for r in [o, *hs]]
    for r in rings:
        assert all(0 <= v <= G for v in r.ravel())
        assert O.ring_problems([tuple(p) for p in r.tolist()]) == [], r.tolist()
    for i in range(len(rings)):
        for j in range(i + 1, len(rings)):
            a, b = rings[i].tolist(), rings[j].tolist()
            for k in range(len(a)):
                for q in range(len(b)):
                    rel = O.segments_relation(a[k], a[(k + 1) % len(a)], b[q], b[(q + 1) % len(b)])
                    assert rel == "none", (i, j, rel)
    for o, hs in parts:
        assert O.twice_area([tuple(p) for p in o.tolist()]) > 0  # outer clockwise on screen (y down)
        for h in hs:
            assert O.twice_area([tuple(p) for p in h.tolist()]) < 0
            # hole inside its outer: a hole vertex is strictly inside (rings are disjoint)
            x, y = h[0]
            cnt = 0
            ol = o.tolist()
            for k in range(len(ol)):
                (x0, y0), (x1, y1) = ol[k], ol[(k + 1) % len(ol)]
                if (y0 <= y < y1) or (y1 <= y < y0):
                    num = (x0 - x) * (y1 - y0) + (y - y0) * (x1 - x0)
                    if (num > 0) if (y1 - y0) > 0 else (num < 0):
                        cnt += 1
            assert cnt % 2 == 1


def check_target(m, cfg, crop="native"):
    masks_, (H, W) = apply_crop({"nodule": m}, crop)
    src = masks_["nodule"]
    text, objs, info = build_target(masks_, LIDC_SCHEMA, cfg)
    o = objs["nodule"]
    assert o.status != "fail", (o.status, o.problems)
    parts = parse(text, LIDC_SCHEMA)["nodule"]
    # round trip: parse -> identical rings, canonical, deterministic serialization
    assert len(parts) == len(o.parts)
    for (a, ah), (b, bh) in zip(parts, o.parts):
        assert np.array_equal(a, b) and len(ah) == len(bh) and all(np.array_equal(x, y) for x, y in zip(ah, bh))
    assert is_canonical(parts)
    assert json.dumps(json.loads(text), separators=(",", ":")) == text
    if src.any():
        oracle_legal(parts)
    # exact oracle re-rasterization
    rings = [r for p, hs in parts for r in [p, *hs]]
    ras = O.rasterize([[tuple(v) for v in r.tolist()] for r in rings], W, H) if rings else np.zeros((H, W), bool)
    assert np.array_equal(ras, rasterize_parts(parts, W, H))
    assert tuple(betti_ref(ras)) == tuple(betti_ref(src))
    u = np.logical_or(ras, src).sum()
    iou = 1.0 if u == 0 else np.logical_and(ras, src).sum() / u
    assert iou >= cfg.iou_min or o.status == "fallback"
    if o.status == "fallback":
        assert iou == 1.0
    return text, ras, src


@settings(max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(m=masks(), cfgname=st.sampled_from(sorted(CONFIGS)), crop=st.sampled_from(["native", "pad_div", "pad_sq"]))
def test_builder_legal_canonical_topology(m, cfgname, crop):
    check_target(m, CONFIGS[cfgname], crop)


@settings(max_examples=120, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(m=masks(), crop=st.sampled_from(["native", "pad_div"]))
def test_lossless_idempotent(m, crop):
    cfg = CONFIGS["lossless"]
    text, ras, src = check_target(m, cfg, crop)
    assert np.array_equal(ras, src)
    masks_, _ = apply_crop({"nodule": m}, crop)
    text2, _, _ = build_target({"nodule": ras}, LIDC_SCHEMA, cfg)
    assert text2 == text


@settings(max_examples=80, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(m=masks(), cfgname=st.sampled_from(["hyb3", "dp5"]))
def test_lossy_rebuild_is_legal_and_topology_stable(m, cfgname):
    """Lossy targets are not idempotent in general, but rebuilding from their raster stays valid."""
    cfg = CONFIGS[cfgname]
    _, ras, src = check_target(m, cfg)
    check_target(ras, cfg)


@settings(max_examples=60, deadline=None)
@given(m=masks(), cfgname=st.sampled_from(sorted(CONFIGS)))
def test_deterministic(m, cfgname):
    a = build_target({"nodule": m}, LIDC_SCHEMA, CONFIGS[cfgname])[0]
    b = build_target({"nodule": m.copy()}, LIDC_SCHEMA, CONFIGS[cfgname])[0]
    assert a == b


PINCH_CASES = {
    "diag2": np.array([[1, 0], [0, 1]], bool),
    "anti2": np.array([[0, 1], [1, 0]], bool),
    "checker": (np.add.outer(np.arange(5), np.arange(5)) % 2 == 0),
    "ring_with_diag_hole": np.array([[1, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 1]], bool),
    "hole_touching_outside": np.array([[1, 1, 1, 0], [1, 0, 1, 1], [1, 1, 1, 1]], bool),
    "one_pixel": np.ones((1, 1), bool),
    "thin_line": np.ones((1, 9), bool),
    "staircase": np.tril(np.ones((7, 7), bool)),
}


@pytest.mark.parametrize("name", sorted(PINCH_CASES))
@pytest.mark.parametrize("cfgname", ["hyb3", "lossless", "dp5"])
@pytest.mark.parametrize("crop", ["native", "pad_div"])
def test_pinches_and_thin_features(name, cfgname, crop):
    check_target(PINCH_CASES[name], CONFIGS[cfgname], crop)


def test_exact_rings_bevel_pinch():
    m = PINCH_CASES["checker"]
    text, objs, info = build_target({"nodule": m}, LIDC_SCHEMA, BuildConfig("exact_rings", {}, iou_min=1.0))
    assert info["bevels"] > 0 and objs["nodule"].status in ("ok", "fallback")
    check_target(m, BuildConfig("exact_rings", {}, iou_min=1.0))


def test_riga_schema_and_rim():
    z = np.load(Path(__file__).resolve().parents[2] / "fixtures" / "riga_experts.npz")
    item = sorted({k.rsplit("__", 1)[0] for k in z.files})[0]
    masks_ = {"disc": z[f"{item}__disc"][0].astype(bool), "cup": z[f"{item}__cup"][0].astype(bool)}
    text, objs, info = build_target(masks_, RIGA_SCHEMA, CONFIGS["hyb3"], RIGA_DERIVED)
    d = json.loads(text)
    assert list(d) == ["disc", "cup"] and isinstance(d["disc"][0][0], int)
    assert info["derived_ok"]["rim"]
    assert all(o.status != "fail" for o in objs.values())


def test_determinism_across_processes():
    code = (
        "import sys, hashlib, numpy as np; sys.path[:0] = sys.argv[1:3];"
        "from targets.builder import *;"
        "rng = np.random.default_rng(3); m = rng.random((30, 30)) < 0.55;"
        "print(hashlib.sha1(build_target({'nodule': m}, LIDC_SCHEMA, BuildConfig('minseg_hybrid', {'tau': 3}))[0]"
        ".encode()).hexdigest())")
    here = Path(__file__).resolve().parents[1]
    outs = {subprocess.run([sys.executable, "-c", code, str(here), str(here.parent / "common")],
                           capture_output=True, text=True, env={"PYTHONHASHSEED": str(s)}).stdout for s in (1, 2)}
    assert len(outs) == 1 and len(next(iter(outs))) == 41


def fast_check(m, cfg, S):
    """Large crops (grid unit close to a pixel): fast kernels (validated against the oracle)."""
    from targets import fastgeom as F
    from targets.builder import betti
    big = np.zeros((S, S), bool)
    big[S // 3:S // 3 + m.shape[0], S // 2:S // 2 + m.shape[1]] = m
    text, objs, info = build_target({"nodule": big}, LIDC_SCHEMA, cfg)
    o = objs["nodule"]
    parts = parse(text, LIDC_SCHEMA)["nodule"]
    rings = [r for p, hs in parts for r in [p, *hs]]
    if big.any():
        assert F.legality(rings)["ok"]
    ras = rasterize_parts(parts, S, S)
    assert betti(ras) == betti(big)
    u = np.logical_or(ras, big).sum()
    iou = 1.0 if u == 0 else np.logical_and(ras, big).sum() / u
    assert o.status != "fail" and (iou >= cfg.iou_min)
    return text, ras, big


@settings(max_examples=120, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(m=masks(max_side=14), S=st.sampled_from([125, 333, 499, 500]), cfgname=st.sampled_from(sorted(CONFIGS)))
def test_large_crops_up_to_500(m, S, cfgname):
    """Crops <= 500 px: one grid unit <= 1/2 px; the builder must always succeed (pinches included)."""
    fast_check(m, CONFIGS[cfgname], S)


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(m=masks(max_side=14), S=st.sampled_from([333, 500]))
def test_large_crops_lossless_idempotent(m, S):
    text, ras, big = fast_check(m, CONFIGS["lossless"], S)
    assert np.array_equal(ras, big)
    assert build_target({"nodule": ras}, LIDC_SCHEMA, CONFIGS["lossless"])[0] == text
