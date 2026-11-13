"""Scoring on the CUDA checker (objdel._gpucheck): polygon -> mask -> defects -> IoU_tc, F5, g_area and the
size-filtered signature, all on the device against expert tables resident there; only per-object numbers come back.
The reference is the host scorer it replaced (oracles/score_numpy.py, numpy + OpenCV): every term must be
bit-identical, on hypothesis-generated rater stacks and masks and on polygon completions. On a CUDA box the chain runs
this file on the GPU build before training (policy.chain); elsewhere it runs on the same source compiled for the CPU.
"""
import dataclasses
import gc
import json

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from oracles import score_numpy as host

from objdel.core.provenance import repo_root
from objdel.geom import engine as _gpucheck
from objdel.geom.shape import to_masks
from objdel.geom.text import parse, serialize
from objdel.infra.vast import chain as P
from objdel.verify.metrics import metrics
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_batch, score_many, score_mask
from tests.oracles import cuda_parity as C

ROOT = repo_root()


def disk(n, r, cy, cx):
    y, x = np.mgrid[:n, :n]
    return np.hypot(y - cy, x - cx) <= r


def parts_of(rings, roles):
    """A fuzz case's rings as to_mask parts, or None when a hole's role does not name an outer ring."""
    roles = [-1] * len(rings) if roles is None else list(roles)
    if any(p != -1 and not (0 <= p < len(rings) and roles[p] == -1) for p in roles):
        return None
    return [(np.asarray(rings[i]).reshape(-1, 2), [np.asarray(rings[j]).reshape(-1, 2)
                                                    for j, p in enumerate(roles) if p == i])
            for i, p in enumerate(roles) if p == -1]  # fmt: skip


def fuzz_objects():
    out = []
    for _, rings, roles, W, H, G in C.fuzz_polygons(20, seed=3):
        parts = parts_of(rings, roles)
        if G == 1000 and parts is not None and sum(1 + len(h) for _, h in parts) <= 10:
            out.append((parts, W, H))
    return out


def test_to_masks_equals_polycore():
    """Masks and defects equal those the archived CPU checker recorded (tests/goldens/polycore.npz)."""
    objs = fuzz_objects() + [([], 5, 4)]  # + an object with no parts: a legal empty mask
    assert len(objs) > 300
    gpu = to_masks(objs)
    assert [i for i, ((p, W, H), g) in enumerate(zip(objs, gpu, strict=True))
            if g.mask.shape != (H, W) or C.geometry_diff(g, C.polygon_reference(*C.object_case(p), W, H))] == []
    assert any(not g.legal for g in gpu) and any(g.legal and g.mask.any() for g in gpu)


def test_to_masks_of_nothing():
    assert to_masks([]) == [] and _gpucheck.build() in ("cuda", "cpu")


# ------------------------------------------------------------------ masks: device scorer vs the host scorer

def blobs(rng, H, W, n, holes=True):
    """n random ellipses (some with an elliptic hole), some 1-4 px specks."""
    y, x = np.mgrid[:H, :W]
    m = np.zeros((H, W), bool)
    for _ in range(n):
        cy, cx = rng.uniform(-2, H + 2), rng.uniform(-2, W + 2)
        ry, rx = rng.uniform(0.5, max(1.0, H / 2)), rng.uniform(0.5, max(1.0, W / 2))
        e = ((y - cy) / ry) ** 2 + ((x - cx) / rx) ** 2
        m |= e <= 1
        if holes and rng.random() < 0.4:
            m &= ~(e <= rng.uniform(0.05, 0.5))
    for _ in range(int(rng.integers(0, 4))):
        sy, sx = int(rng.integers(0, H)), int(rng.integers(0, W))
        m[sy:sy + int(rng.integers(1, 3)), sx:sx + int(rng.integers(1, 3))] = True
    return m


def rater_stack(rng, H, W):
    """k = 2..5 raters around one or two shared blobs: grown, shrunk, shifted, holed or speckled copies."""
    base = blobs(rng, H, W, int(rng.integers(1, 3)), holes=False)
    out = []
    for _ in range(int(rng.integers(2, 6))):
        r = np.roll(base, (int(rng.integers(-1, 2)), int(rng.integers(-1, 2))), (0, 1))
        r = r | (blobs(rng, H, W, 1, holes=False) & (rng.random() < 0.3))
        if rng.random() < 0.3:
            r &= ~blobs(rng, H, W, 1, holes=False)
        out.append(r)
    return np.stack(out)


def predictions(rng, stack, cache):
    """Masks a model could draw: rater masks and their union, core and vote, perturbed copies, noise, nothing."""
    k, H, W = stack.shape
    noise = rng.random((H, W))
    out = [*stack, cache.union, cache.core, cache.vote, np.zeros((H, W), bool), np.ones((H, W), bool),
           noise < 0.5, noise < 0.08, blobs(rng, H, W, int(rng.integers(1, 5)))]
    for r in stack[:2]:
        out += [r ^ (noise < 0.05), np.roll(r, int(rng.integers(-3, 4)), 1), r & ~blobs(rng, H, W, 2)]
    return out


def device_terms(m, cache) -> dict:
    return dataclasses.asdict(score_mask(m, cache))


@settings(max_examples=60, deadline=None)
@given(st.integers(0, 2**32 - 1), st.integers(1, 70), st.integers(1, 70))
def test_score_mask_equals_the_host_scorer(seed, H, W):
    """Every term bit-identical: IoU_tc (same-signature raters, else the vote), F5 (band depth: float64 sums),
    g_area, the size-filtered signature and its rater share; frames of every width around the 32-pixel words."""
    rng = np.random.default_rng(seed)
    stack = rater_stack(rng, H, W)
    try:
        cache = ExpertCache.build(stack)
    except ValueError:  # e.g. raters without any contour pair: not a cache the reward can have
        return
    for m in predictions(rng, stack, cache):
        assert device_terms(m, cache) == host.score_mask(m, cache)


@settings(max_examples=6, deadline=None)
@given(st.integers(0, 2**32 - 1), st.integers(180, 320), st.integers(180, 320))
def test_score_mask_on_large_frames_equals_the_host_scorer(seed, H, W):
    """Frames of 180-320 px (review 2026-11-07: the other cases stop at 70 px): many rows per crop, many labels
    merging at once in the device union-find, long band-depth sums. The chain runs it on the GPU and under memcheck."""
    rng = np.random.default_rng(seed)
    stack = rater_stack(rng, H, W)
    try:
        cache = ExpertCache.build(stack)
    except ValueError:
        return
    for m in predictions(rng, stack, cache):
        assert device_terms(m, cache) == host.score_mask(m, cache)


def test_score_mask_edge_cases_equal_the_host_scorer():
    """A full-frame core (its band distances are the 1e6 sentinel), components and holes exactly at the size
    filter, a checkerboard (one 8-connected part, many holes), single pixels, a frame of one row."""
    full = np.ones((2, 9, 40), bool)
    full[1, 0, 0] = False
    sq = np.zeros((3, 40, 37), bool)
    sq[:, 5:30, 4:33] = True
    sq[1, 10:12, 10:12] = False  # a 4 px hole: below min_area 5
    sq[2, 10:13, 10:12] = False  # a 6 px hole
    caches = [ExpertCache.build(s) for s in (full, sq)]
    fy, fx = np.mgrid[:9, :40]
    cases = [np.ones((9, 40), bool), (fy + fx) % 2 == 0, fy % 3 == 0]
    yy, xx = np.mgrid[:40, :37]
    sq_cases = [sq[1], sq[2], (yy + xx) % 2 == 0, (yy % 2 == 0) & (xx % 2 == 0), np.eye(40, 37, dtype=bool)]
    one = np.zeros((40, 37), bool)
    one[0, 36] = True
    for m in sq_cases[:2]:
        for size in (4, 5, 6):  # a part of exactly min_area - 1, min_area, min_area + 1 pixels beside it
            extra = m.copy()
            extra[35, 2:2 + size] = True
            sq_cases.append(extra)
    for c, ms in zip(caches, (cases, [*sq_cases, one]), strict=True):
        for m in ms:
            assert device_terms(m, c) == host.score_mask(m, c)
    line = np.stack([np.arange(70) < 30, np.arange(70) < 34])[:, None, :]
    c = ExpertCache.build(line)
    for m in (line[0, 0], ~line[1, 0], np.arange(70) % 3 == 0):
        assert device_terms(m[None], c) == host.score_mask(m[None], c)


def test_score_mask_refuses_another_frame():
    c = ExpertCache.build(np.stack([disk(20, r, 10, 10) for r in (4, 5)]))
    with pytest.raises(ValueError, match="shape"):
        score_mask(np.zeros((20, 21), bool), c)


def test_tables_are_uploaded_once_and_freed_with_the_cache():
    before = _gpucheck.resident()
    c = ExpertCache.build(np.stack([disk(20, r, 10, 10) for r in (4, 5)]))
    h = c.device()
    assert c.device() == h and _gpucheck.resident() == before + 1
    score_mask(c.vote, c)
    assert _gpucheck.resident() == before + 1
    del c
    gc.collect()
    assert _gpucheck.resident() == before


def test_resident_tables_are_bounded_and_rebuild_alike(monkeypatch):
    """device_many keeps at most MAX_TABLES table sets: the least recently used caches lose theirs (never one of the
    running call's), rebuild the same tables on their next use and score the same; every handle is dropped once."""
    from collections import OrderedDict

    from objdel.verify.reward import cache as rc
    monkeypatch.setattr(rc, "MAX_TABLES", 2)
    monkeypatch.setattr(rc, "_TABLES", OrderedDict())  # other tests' caches keep their tables
    dropped, drop = [], _gpucheck.drop_tables
    monkeypatch.setattr(_gpucheck, "drop_tables", lambda h: (dropped.append(h), drop(h)))
    gc.collect()
    before = _gpucheck.resident()
    cs = [ExpertCache.build(np.stack([disk(20, r + i, 10, 10) for r in (4, 5)])) for i in range(4)]
    want = [(score_mask(c.vote, c), c.rater_sigs, c.perimeter) for c in cs]
    assert [c._device for c in cs[:2]] == [None, None] and _gpucheck.resident() == before + 2
    rc.device_many([cs[2]])  # least recent first: 2, 3; touched, 3, 2; so building 0 drops 3
    rc.device_many([cs[0]])
    assert [c._device is None for c in cs] == [False, True, False, True] and _gpucheck.resident() == before + 2
    rc.device_many(cs)  # a call keeps all its caches, even past the bound
    assert None not in [c._device for c in cs] and _gpucheck.resident() == before + 4
    assert [(score_mask(c.vote, c), c.rater_sigs, c.perimeter) for c in cs] == want
    del cs, want
    gc.collect()
    assert _gpucheck.resident() == before and len(dropped) == len(set(dropped)) == 7


# ------------------------------------------------------------------ polygons: completions end to end

@pytest.fixture(scope="module")
def jobs():
    """LIDC (one object) and RIGA-like (two objects) jobs: fuzz outlines on the 0-1000 grid, an empty answer, an
    illegal bow-tie, a hole, garbage, a truncated and a wrong-schema completion."""
    nod = ExpertCache.build(np.stack([disk(40, r, 20, 21) for r in (8, 9, 10)]))
    disc = ExpertCache.build(np.stack([disk(60, r, 30, 30) for r in (20, 21, 22)]))
    cup = ExpertCache.build(np.stack([disk(60, r, 31, 30) for r in (8, 9, 9)]))
    lidc, riga = {"nodule": nod}, {"disc": disc, "cup": cup}
    texts = []
    for parts, _, _ in fuzz_objects():
        if all(np.asarray(r).min() >= 0 and np.asarray(r).max() <= 1000 and len(r) >= 3
               for o, hs in parts for r in (o, *hs)):
            texts.append(serialize({"nodule": parts}, ("nodule",)))
    texts = [t for t in texts if parse(t, ("nodule",)).ok][:150]
    texts += ['{"nodule":[]}', '{"nodule":[{"outer":[[100,100],[900,900],[900,100],[100,900]],"holes":[]}]}',
              '{"nodule":[{"outer":[[200,200],[800,200],[800,800],[200,800]],"holes":[[[400,400],[400,600],'
              '[600,600],[600,400]]]}]}', "garbage", '{"disc":[]}']
    out = [(t, False, ("nodule",), lidc, 40, 40) for t in texts]
    out.append((texts[0], True, ("nodule",), lidc, 40, 40))

    def square(a, b):
        return [(np.array([[a, a], [b, a], [b, b], [a, b]]), [])]

    out += [(serialize({"disc": square(150, 850), "cup": square(a, b)}, ("disc", "cup")), False, ("disc", "cup"), riga,
             60, 60) for a, b in ((350, 650), (400, 600), (100, 900))]
    return out


def test_score_batch_equals_the_host_scorer(jobs):
    gpu = score_batch(jobs)
    ref = [host.score_completion(*j) for j in jobs]
    for k, v in gpu.items():
        want = [r[k] for r in ref]
        assert (v.tolist() if isinstance(v, np.ndarray) else v) == want, k
    assert 0 < gpu["legal"].sum() < len(jobs) and gpu["format_ok"].sum() < len(jobs)


def test_score_many_keeps_job_order(jobs):
    one = [score_many([j])[0] for j in jobs]
    many = score_many(jobs)
    for (t1, o1), (t2, o2) in zip(one, many, strict=True):
        assert json.dumps(t1, sort_keys=True) == json.dumps(t2, sort_keys=True) and o1 == o2


def test_score_polygons_returns_numbers_only(jobs):
    """The device scorer's result: per-object defects, area, terms and signature, never a mask."""
    c = jobs[0][3]["nodule"]
    r = _gpucheck.score_polygons(np.array([[0, 0], [1000, 0], [1000, 1000]], np.int32), np.array([0, 3], np.int64),
                                 np.array([0, 1], np.int64), None, np.array([40], np.int32), np.array([40], np.int32),
                                 1000, np.array([c.device()], np.int64))
    assert set(r) == {"defects", "area", "sum_x", "sum_y", "iou_tc", "f5", "g_area", "topo_freq", "sig_start", "sig"}
    assert all(v.size <= 13 for v in r.values())
    with pytest.raises(ValueError, match="frame"):
        _gpucheck.score_polygons(np.zeros((0, 2), np.int32), np.array([0], np.int64), np.array([0, 0], np.int64),
                                 None, np.array([41], np.int32), np.array([40], np.int32), 1000,
                                 np.array([c.device()], np.int64))


def test_exam_metrics_equal_the_host_scorer(jobs):
    lidc = [j for j in jobs if j[2] == ("nodule",)][:16]
    texts, flags, caches = [j[0] for j in lidc], [j[1] for j in lidc], lidc[0][3]
    m = metrics(texts, {}, 40, 40, ("nodule",), flags, caches)
    ref = [host.score_completion(*j) for j in lidc]
    legal = np.array([r["legal"] for r in ref])
    assert m["samples"]["iou_tc"] == (np.array([r["iou_tc"] for r in ref]) * legal).tolist()
    assert m["samples"]["topo_any"] == [r["topo_any"] for r in ref]


def test_every_scoring_chain_builds_and_checks_the_cuda_checker_first():
    gpu, train, select, probe = P.chain(ROOT / "specs/r1_qwen3vl2b.yaml", "r1")
    assert gpu[:2] == ["bash", "-c"] and "run_box.sh --build-only && " in gpu[2]
    assert gpu[2].endswith("tests/verify/test_gpu_scoring.py")  # this file, on the GPU build, before training
    assert [train[1], select[1], probe[1]] == ["objdel.cli.train.sft", "objdel.cli.train.select",
                                               "objdel.cli.train.probe"]
