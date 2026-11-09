# ADR-0001: One source file builds the checker for the CPU and for nvcc

**Status:** Accepted
**Date:** 2026-11-02 (commit `734b1ec`, first `cuda/gpucheck.cpp`)
**Deciders:** Repository owner

## Context

The reward needs an exact geometry checker. Phase 1 built one for the CPU, called polycore.
It worked, but scoring 256 answers per step on the CPU kept the GPU waiting.
A GPU checker must give the same result as polycore on every input. A second, separate implementation is a second place for bugs.

Forces:

- The checker must be bit-exact. A different result on the CPU and the GPU would make tests meaningless.
- The Mac has no GPU. The development machine must still run every checker test.
- The checker is about 1,800 lines. Two copies would drift.

## Decision

Write every algorithm once, in `cuda/gpucheck.cpp`, as `__host__ __device__` functors.
A function `for_each` runs a functor over independent items.
The C++ compiler builds `for_each` as a plain loop. nvcc builds it as a kernel launch.
Only a small "platform" section differs: memory, copies, streams and four atomic helpers.
The raster has one extra GPU-only step, the warp ballot.

`CMakeLists.txt` builds `gpucheck` as a static library so that nanobind flags never reach nvcc.
`make cuda` sets `OBJDEL_CUDA=ON` and compiles the same file with nvcc.

## Options considered

### Option A: One source, two compilers (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | Medium. Every algorithm must be a functor over independent items |
| Cost | Low. One implementation to maintain |
| Test path | The CPU build is the test path of the GPU code: one loop iteration per GPU thread |
| Team familiarity | Medium |

**Pros:** no drift between CPU and GPU. Every checker test runs on a Mac.
**Cons:** the CPU build is slow by design (all-pairs edge tests, no sort). It is not a production CPU checker.

### Option B: Keep polycore on the CPU and write separate CUDA kernels

| Dimension | Assessment |
|---|---|
| Complexity | High. Two code bases with equal results |
| Cost | Double the maintenance |
| Test path | Needs a GPU to test the kernels |

**Pros:** each side uses its natural algorithms (sweep line on the CPU).
**Cons:** equality must be proved on every change. Mac development cannot run the kernels.

### Option C: Python with numba CUDA

The lab prototyped numba in `experiments/phase1/lab/t_polycore/tests/numba_proto.py`.
It was not promoted. The lab file is the only record. This record does not rank it further.

## Trade-off analysis

Option A trades CPU speed for certainty. The CPU speed does not matter because production runs on the GPU.
Certainty matters because the reward is a proof object.
Option B would give better CPU speed and a heavier proof burden.

## Consequences

- Easier: one change, one place. The Mac runs the full checker suite.
- Easier: bit-for-bit equality tests. polycore results are recorded once in `tests/goldens/polycore.npz`.
- Harder: every algorithm must avoid data-dependent global state. Functors take items by index.
- Harder: the build must forbid FMA contraction on both compilers (`-ffp-contract=off`, `--fmad=false`), or float64 terms differ by one rounding.
- Revisit: 1,000-vertex shapes at batch 1,024 run slower than the 8-thread CPU path (`experiments/phase2/perf/p1_results.md`). The cause is one thread per shape in the witness tests.

## Evidence

| Check | Result | Source |
|---|---|---|
| CUDA vs CPU checker on 4 case sets | 195,809 of 195,809 cases identical (RTX 4090) | `experiments/cuda_checker/RESULTS.md` |
| Second architecture | 0 mismatches on 8,559 fuzz polygons and 6,024 fuzz masks (sm_86) | `experiments/phase2/perf/k1_results.md` |
| Throughput | 7.74 M polygons/s, RTX 5090, 16 vertices, batch 1,024 | `experiments/phase2/perf/p1_results.md` |

## Action items

1. [x] Port polycore to functors (`cuda/gpucheck.cpp`).
2. [x] Archive polycore and record its results (`archive/cpu_checker`, `tests/goldens/polycore.npz`).
3. [ ] Revisit the 1,000-vertex batch path (`experiments/phase2/perf/p1_results.md`; related costs in issue CK-4 of `perf/issues.md`).
