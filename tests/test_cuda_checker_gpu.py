"""GPU group of the CUDA checker, for the preflight on a rented box (skipped where `nvidia-smi -L` lists no GPU).

On a GPU box objdel._gpucheck must be the nvcc build (`make cuda`): these tests fail, not skip, when it is not.
They then run large batches in one call each (many shapes per kernel launch) and compare with the CPU reference;
tests/test_cuda_checker.py runs the rest of the equivalence suite on the same build."""
import shutil
import subprocess

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


pytestmark = pytest.mark.skipif(_gpus() == 0, reason="no CUDA device (nvidia-smi -L)")


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
    masks += [rng.random((512, 512)) < .5, rng.random((2048, 2048)) < .1]  # long cycles: many jump rounds
    out = C.run_masks(C.pack_masks(masks))
    for i, m in enumerate(masks):
        assert C.mask_diff(out, i, C.mask_reference(np.asarray(m, bool))) is None, i


def test_repeated_calls_are_identical():
    p = C.pack_polygons([(r, rl, W, H) for _, r, rl, W, H, G in C.fuzz_polygons(20, seed=13) if G == 1000])
    a, b = C.run_polygons(p), C.run_polygons(p)
    assert all(np.array_equal(a[k], b[k]) for k in a)
