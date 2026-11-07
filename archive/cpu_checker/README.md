# The CPU checker (polycore), archived 2026-11-07

**Nothing here is built, imported or tested any more.** Checking and scoring run in the CUDA checker,
`objdel._gpucheck` (`cuda/`, one source for the CPU and nvcc builds). See `docs/research/notes/gpu_scoring.md`.

## What it was

polycore was the Phase 1 exact integer polygon core (plan P1-D1, D14). It was a header-only C++17 library with a
nanobind extension:

- it rasterized integer rings on the 1000 grid (even-odd and nonzero masks, pixel-centre sampling, top-left ties);
- it validated them into 13 defect classes (crossings, touches, winding, nesting, range);
- it computed signed areas;
- it ran a threaded batch scorer.

It was the reward's checker in Phase 1 and early Phase 2: polycore made the masks and numpy / OpenCV scored them on a
CPU process pool. The CUDA checker was then ported from it with the same integer arithmetic and verified bit for bit
against it (experiments/cuda_checker/RESULTS.md). After that the CUDA checker took over the scoring as well, so
polycore had no runtime use left.

| Path | Was |
|---|---|
| `cpp/polycore/polycore.hpp` | the C++17 core (was `cpp/polycore/`) |
| `cpp/polycore/tests/fuzz_main.cpp` | native ASan/UBSan fuzz driver (was `make fuzz`) |
| `src/objdel/_polycore.cpp` | nanobind binding (module `objdel._polycore`) |
| `src/objdel/checker/polycore.py` | Python face: `pack`, `rasterize`, `validate`, `signed_area2`, `score_batch` |
| `tests/test_polycore*.py` | unit, property and K1/K2 oracle tests |
| `experiments/prepare_data.py` | built `data/checker_eval/` (LIDC polygons and masks with polycore's references) |
| `record_goldens.py` | recorded polycore's results for the CUDA checker's tests (below) |

## What replaced each use

- **Caps, defect classes and input checks** now live in `cuda/gpucheck.h`.
- **Masks and legality** (`objdel.geometry.to_mask`, used for targets and probes) come from `to_masks` on the CUDA
  checker.
- **Packing** is `objdel.geometry.shape.pack_rings`.
- **Signed area** is `objdel.geometry.shape.signed_area2`, in exact Python integers.
- **Reference results:** before the archive, `record_goldens.py` ran polycore on every polygon case the tests and
  `experiments/cuda_checker/eval.py` use. It stored the defect counts and the SHA-256 of each even-odd mask, keyed per
  case, in `tests/goldens/polycore.npz`. The tests now compare the CUDA checker with those recorded results
  (`experiments/cuda_checker/compare.py: polygon_reference`).
- **K1/K2 oracle tests:** the slow, exact Python oracles (`tests/oracles/`) now check the CUDA checker directly, in
  `tests/test_checker_oracle.py`.
- **`data/checker_eval/`** still holds polycore's references, and `eval.py` still reads them. The set cannot be rebuilt
  without polycore, so keep a copy.

Running any of this again needs the old build lines: `nanobind_add_module(_polycore ...)` in `CMakeLists.txt` and
`-I cpp`. They are in git history before the archive.
