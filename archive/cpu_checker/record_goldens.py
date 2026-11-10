"""Records the CPU checker's results (polycore: 13 defect counts and the even-odd mask's SHA-256) for every polygon
case that the tests and experiments/cuda_checker/eval.py compare the CUDA checker with, into tests/goldens/polycore.npz.
Run once, on 2026-11-07, just before polycore was archived:
    .venv/bin/python -I archive/cpu_checker/record_goldens.py
It needs objdel.checker.polycore and objdel._polycore, which are no longer built (archive/cpu_checker/README.md): kept
as the record of how the goldens were made. Cases are keyed by experiments.cuda_checker.compare.case_key.
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "tests")]

from experiments.cuda_checker import compare as C  # noqa: E402
from objdel.checker import polycore  # noqa: E402  (archive/cpu_checker/src/objdel/checker/polycore.py)

RECORDED = {}


def record(rings, roles, W, H, G=1000) -> dict:
    shape = polycore.pack(rings)
    ref = dict(defects=polycore.validate(shape, roles=np.array(roles, np.int32), G=G)["total"].astype(np.int64),
               evenodd=C.mask_sha(polycore.rasterize(shape, W, H, G)["evenodd"]))
    RECORDED[C.case_key(rings, roles, W, H, G)] = ref
    return ref


def main():
    C.polygon_reference = record  # the test helpers below build their references through it
    import test_checker_host
    import test_cuda_checker
    import test_cuda_checker_gpu
    import test_gpu_scoring

    for n, seed in [(60, 1), (250, 11), (20, 0), (150, 0)]:  # test_cuda_checker, _gpu, eval --quick, eval
        for _, rings, roles, W, H, G in C.fuzz_polygons(n, seed):
            record(rings, [-1] * len(rings) if roles is None else list(roles), W, H, G)
    for n, seed in [(1200, 21), (1000, 22)]:  # test_cuda_checker_gpu: arena reuse, two streams
        test_cuda_checker_gpu._shape_pool(n, seed)
    for rings, roles, W, H in test_cuda_checker.test_edge_cases_match_polycore.pytestmark[0].args[1]:
        record(rings, [-1] * len(rings) if roles is None else roles, W, H)
    record([test_cuda_checker.SQ, test_cuda_checker.HOLE], [-1, 0], 16, 16)
    for parts, W, H in test_gpu_scoring.fuzz_objects() + [([], 5, 4)] + test_checker_host._objects(1, 400):
        record(*C.object_case(parts), W, H)
    keys = sorted(RECORDED)
    C.GOLDEN.parent.mkdir(exist_ok=True)
    np.savez_compressed(C.GOLDEN, key=np.frombuffer(b"".join(keys), np.uint8).reshape(-1, 16),
                        defects=np.stack([RECORDED[k]["defects"] for k in keys]),
                        evenodd=np.stack([np.frombuffer(bytes.fromhex(RECORDED[k]["evenodd"].decode()), np.uint8)
                                          for k in keys]))
    print(f"{len(keys)} cases -> {C.GOLDEN} ({C.GOLDEN.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
