"""GPU group of the CUDA checker, for the preflight on a rented box (skipped where `nvidia-smi -L` lists no GPU).

On a GPU box objdel._gpucheck must be the nvcc build (`make cuda`): these tests fail, not skip, when it is not.
They then run large batches in one call each (many shapes per kernel launch) and compare with the CPU checker's
recorded results (tests/goldens/polycore.npz);
tests/test_cuda_checker.py runs the rest of the equivalence suite on the same build."""
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from experiments.cuda_checker import compare as C
from objdel import _gpucheck as gc


def _gpus() -> int:
    if not shutil.which("nvidia-smi"):
        return 0
    try:
        out = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return 0
    return sum(line.startswith("GPU ") for line in out.splitlines())


pytestmark = [pytest.mark.gpu, pytest.mark.skipif(_gpus() == 0, reason="no CUDA device (nvidia-smi -L)")]


def test_module_is_the_cuda_build():
    assert gc.build() == "cuda", "objdel._gpucheck is the CPU build: run `make cuda` on the box first"
    assert gc.device_count() == _gpus()


def test_large_polygon_batches_match_polycore():
    cases = C.fuzz_polygons(250, seed=11)  # ~6,000 shapes, every family, a few at the 4096 x 4096 cap
    for G in sorted({c[5] for c in cases}):
        sub = [c for c in cases if c[5] == G]
        p = C.pack_polygons([(r, rl, W, H) for _, r, rl, W, H, _ in sub])
        out = C.run_polygons(p, G)
        for s, (fam, rings, roles, W, H, _) in enumerate(sub):
            ref = C.polygon_reference(rings, [-1] * len(rings) if roles is None else roles, W, H, G)
            assert C.polygon_diff(C.gpu_polygon(out, p, s), ref) is None, (fam, s)


def test_large_mask_batch_matches_polygonize():
    masks = [m for _, m in C.fuzz_masks(3000, seed=12)]
    rng = np.random.default_rng(12)
    masks += [rng.random((512, 512)) < .5, rng.random((2048, 2048)) < .1]  # many rings, deep union-find chains
    out = C.run_masks(C.pack_masks(masks))
    for i, m in enumerate(masks):
        assert C.mask_diff(out, i, C.mask_reference(np.asarray(m, bool))) is None, i


def test_lidc_crop_masks_match_the_oracle():
    """D-G3: all 37,093 LIDC crops in one call equal the numpy oracle (rings, quads, Betti numbers, ring parts), so a
    CUDA-only ring defect fails here. The crops ship in the box payload (data/checker_eval, OBJDEL_DATA_ROOT)."""
    path = Path(os.environ.get("OBJDEL_DATA_ROOT", Path(__file__).resolve().parents[1] / "data"))
    z = dict(np.load(path / "checker_eval" / "masks.npz", allow_pickle=False))
    out = gc.check_masks(z["bits"], z["hw"], z["start"])
    assert len(out["quads"]) == 37093 and len(out["ring_area"]) == 38002 and len(out["pts"]) == 1476508
    for k in range(len(z["hw"])):
        m = z["bits"][z["start"][k]:z["start"][k + 1]].reshape(z["hw"][k]).astype(bool)
        assert C.mask_diff(out, k, C.mask_reference(m)) is None, k


def test_repeated_calls_are_identical():
    p = C.pack_polygons([(r, rl, W, H) for _, r, rl, W, H, G in C.fuzz_polygons(20, seed=13) if G == 1000])
    a, b = C.run_polygons(p), C.run_polygons(p)
    assert all(np.array_equal(a[k], b[k]) for k in a)


def _shape_pool(n: int, seed: int):
    """n fuzz shapes (grid 1000) on random crops, mostly small, some wider than one 1024-column raster window, with
    their polycore references: (cases [(rings, roles, W, H)], references)."""
    rng = np.random.default_rng(seed)
    cases = []
    for _, rings, roles, _, _, G in C.fuzz_polygons(n // 15 + 1, seed=seed):
        if G != 1000:
            continue
        W, H = (int(rng.integers(1, 1500)), int(rng.integers(1, 40))) if rng.random() < .1 else \
            (int(rng.integers(1, 300)), int(rng.integers(1, 300)))
        cases.append((rings, [-1] * len(rings) if roles is None else list(roles), W, H))
    cases = cases[:n]
    return cases, [C.polygon_reference(r, rl, W, H) for r, rl, W, H in cases]


def _mismatches(cases, refs, idx):
    p = C.pack_polygons([cases[i] for i in idx])
    out = C.run_polygons(p)
    return [i for s, i in enumerate(idx) if C.polygon_diff(C.gpu_polygon(out, p, s), refs[i]) is not None]


def arena_reuse(calls: int, seed: int = 21):
    """Calls of random batch sizes 1-1024 over a pool of shapes of random sizes, so the arena and the staging buffers
    grow, and later calls run on memory that earlier, larger ones left behind (nothing is zeroed but the counters):
    every shape must equal polycore in every call."""
    cases, refs = _shape_pool(1200, seed)
    rng = np.random.default_rng(seed)
    bad = []
    for k in range(calls):
        size = [1, 1024, 3][k] if k < 3 else int(rng.integers(1, 1025))
        bad += _mismatches(cases, refs, rng.integers(0, len(cases), size).tolist())
    return bad


def two_stream_chunks(n: int = 1000, seed: int = 22):
    """A batch of n >= 64 shapes runs as two halves on two streams; the same shapes in calls of 32 run on one. Words
    (padding bits included), offsets and defects must be identical, and equal to polycore."""
    cases, refs = _shape_pool(n, seed)
    assert len(cases) == n
    big = C.run_polygons(C.pack_polygons(cases))
    small = [C.run_polygons(C.pack_polygons(cases[a:a + 32])) for a in range(0, n, 32)]
    words = np.concatenate([s["evenodd"] for s in small])
    starts = np.concatenate([[0]] + [s["word_start"][1:] + off for s, off in
                                     zip(small, np.cumsum([0] + [s["evenodd"].size for s in small[:-1]]), strict=True)])
    return (np.array_equal(big["evenodd"], words) and np.array_equal(big["word_start"], starts)
            and np.array_equal(big["defects"], np.concatenate([s["defects"] for s in small]))
            and _mismatches(cases, refs, list(range(n))) == [])


def test_arena_reuse_is_exact():
    assert arena_reuse(200) == []


def test_two_stream_chunks_equal_one_call():
    assert two_stream_chunks()
