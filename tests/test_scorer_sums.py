"""Centroid sums in the CUDA scorer (B3.1; tests P12, I29, G01, G07 of Implementation.md B.9).

`score_polygons` and `score_masks` return per object `sum_x` (sum of the column indices of the mask's pixels) and
`sum_y` (sum of the row indices), int64, exact: they equal the np.nonzero sums. Every field returned before B3.1 is
unchanged: tests/goldens/b_g01.json was written by experiments/phase2/perf/make_b_goldens.py on the Mac CPU build
before the sums were added. Every object is scored against its own `b_tables` (the G01 tables rule).
"""
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pytest

from experiments.phase2.perf.make_b_goldens import b_sets, b_tables, random_masks, score_masks_of
from objdel import _gpucheck as G
from objdel.geometry.shape import unpack_rows
from tests.test_cuda_checker_gpu import _gpus

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "tests" / "goldens" / "b_g01.json"
BOX_INPUTS = ("checker_eval/masks.npz",)
on_gpu = pytest.mark.skipif(_gpus() == 0, reason="no CUDA device (nvidia-smi -L)")


def _masks_npz() -> Path:  # D8's box_input replaces this; a missing input fails (gpu tests never skip on the box)
    p = Path(os.environ.get("OBJDEL_DATA_ROOT", ROOT / "data")) / BOX_INPUTS[0]
    assert p.exists(), f"box input missing: {p}"
    return p


def _sums_ok(r, ms):
    assert r["sum_x"].dtype == np.int64 and r["sum_y"].dtype == np.int64 and len(r["sum_x"]) == len(ms)
    for i, m in enumerate(ms):
        ys, xs = np.nonzero(m)
        assert (int(r["sum_x"][i]), int(r["sum_y"][i])) == (int(xs.sum()), int(ys.sum())), i


def _old_fields_ok(r, gold, name):
    assert sorted(r) == sorted([*gold[name]["fields"], "sum_x", "sum_y"])  # old keys + the two sums
    for k, want in gold[name]["fields"].items():
        a = np.ascontiguousarray(r[k])
        assert [str(a.dtype), list(a.shape), hashlib.sha256(a.tobytes()).hexdigest()] == want, (name, k)


def test_sums_match_numpy_and_old_fields_are_unchanged():
    """P12, I29: the golden's `random` set (203 masks, frames 8..512, all-zero, all-one, one pixel at (63, 63));
    sums equal np.nonzero's; every other field equals the pre-B3.1 golden."""
    gold = json.loads(GOLDEN.read_text())
    ms = list(random_masks(np.random.default_rng(1)))
    sha, r = score_masks_of(ms)
    assert sha == gold["random"]["input_sha256"]
    _sums_ok(r, ms)
    _old_fields_ok(r, gold, "random")


def test_full_frame_and_corner_pixels_are_exact():
    m = np.zeros((40, 33), bool)
    m[0, 0] = m[39, 32] = True
    ms = [np.ones((512, 512), bool), m, np.ones((1, 300), bool), np.ones((300, 1), bool)]
    r = score_masks_of(ms)[1]
    assert int(r["sum_x"][0]) == int(r["sum_y"][0]) == 512 * (511 * 512 // 2) == 66_977_792
    assert (r["sum_x"][1:].tolist(), r["sum_y"][1:].tolist()) == ([32, 299 * 300 // 2, 0], [39, 0, 299 * 300 // 2])


def test_sums_above_2_to_32_do_not_wrap():
    """P12: a 4096 x 4096 all-ones mask: sum_x = sum_y = 4096 * (4095 * 4096 / 2) = 34,351,349,760 > 2^32."""
    r = score_masks_of([np.ones((4096, 4096), bool)])[1]
    assert int(r["sum_x"][0]) == int(r["sum_y"][0]) == 34_351_349_760


def test_polygon_and_its_raster_agree_and_holes_are_excluded():
    """I29: a polygon scored by score_polygons and its even-odd raster scored by score_masks give equal sums (and
    equal old fields); two disjoint blobs in one object give one value; a hole's pixels are not counted."""
    xy = np.array([[100, 100], [900, 100], [900, 900], [100, 900],   # outer
                   [300, 300], [300, 700], [700, 700], [700, 300],   # its hole
                   [0, 0], [50, 0], [50, 50], [0, 50],               # a second object: two blobs
                   [950, 950], [1000, 950], [1000, 1000], [950, 1000]], np.int32)
    ring_start = np.array([0, 4, 8, 12, 16], np.int64)
    shape_start = np.array([0, 2, 4], np.int64)
    roles = np.array([-1, 0, -1, -1], np.int32)
    W = H = np.array([40, 40], np.int32)
    ch = G.check_polygons(xy, ring_start, shape_start, roles, W, H, 1000)
    ms = [unpack_rows(ch["evenodd"][ch["word_start"][i]:ch["word_start"][i + 1]], 40, 40) for i in range(2)]
    caches = [b_tables(m) for m in ms]
    handles = np.array([c.device() for c in caches], np.int64)
    rp = G.score_polygons(xy, ring_start, shape_start, roles, W, H, 1000, handles)
    rm = G.score_masks(np.concatenate([m.view(np.uint8).reshape(-1) for m in ms]), np.array([[40, 40]] * 2, np.int32),
                       np.array([0, 1600, 3200], np.int64), handles)
    _sums_ok(rp, ms)
    for k in ("sum_x", "sum_y", "area", "iou_tc", "f5", "g_area", "topo_freq", "sig_start", "sig"):
        assert rp[k].tobytes() == rm[k].tobytes(), k
    assert not ms[0][20, 20] and ms[0][5, 5]  # the hole is empty, the ring full


@pytest.mark.gpu
@pytest.mark.sanitizer
@on_gpu
def test_g07_sanitize_sums():
    """G07: the only scorer test under compute-sanitizer (R-D2's sanitize step): a G01 subset, sums only."""
    sets = b_sets(_masks_npz(), n_masks=32, n_polys=64, n_crops=32)
    assert [len(ms) for _, ms, _ in sets.values()] == [32, 64, 64, 64, 32]
    for _, ms, r in sets.values():
        _sums_ok(r, ms)


@pytest.mark.gpu  # never under a sanitizer (too long, B.9.5 G07)
@on_gpu
def test_g01_golden():
    """G01: the CUDA build equals the CPU build on every field (old fields through the Mac golden, sums through
    np.nonzero) on the five G01 sets."""
    gold = json.loads(GOLDEN.read_text())
    sets = b_sets(_masks_npz())
    assert sorted(sets) == sorted(gold)
    for name, (sha, ms, r) in sets.items():
        assert sha == gold[name]["input_sha256"], f"{name}: input differs"
        _sums_ok(r, ms)
        _old_fields_ok(r, gold, name)
