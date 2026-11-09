# CUDA engine reference: `gpucheck`

The checker is one C++ source file, `cuda/gpucheck.cpp`, with a header `cuda/gpucheck.h` and a Python binding `cuda/gpucheck_py.cpp`.
The C++ compiler builds it for the CPU. nvcc builds the same file for the GPU. The Python module is `objdel._gpucheck`.
The reasons are in [ADR-0001](../architecture/decisions/0001-single-source-cuda-checker.md) and [ADR-0004](../architecture/decisions/0004-scoring-on-the-gpu.md).

## 1. What the engine does

```text
   packed polygons ----> raster (even-odd mask) ----> validate (13 defects)
                              |
                              +--> scorer: area, IoU vs raters and vote, band depth, signature
   packed masks -------> bit quads, boundary slots, ring labels, corner ranks (polygonize)
   expert tables (once) -> resident on the device until dropped
```

| Job | Entry point | Output |
|---|---|---|
| Draw masks of outlines and validate them | `check_polygons` | Bit-packed even-odd masks, 13 defect counts per shape |
| Trace rings of masks | `check_masks` | Bit-quad counts, rings in `polygonize` order, corners |
| Keep expert tables on the device | `add_tables`, `drop_tables`, `resident` | A handle |
| Score outlines against tables | `score_polygons` | Per-object terms and signature. No mask leaves the device |
| Score masks against tables | `score_masks` | Same, for masks that do not come from an answer |

## 2. Python API (`objdel._gpucheck`)

All arrays must be C-contiguous NumPy arrays of the exact dtype. The binding uses `noconvert`, so a wrong dtype raises `TypeError` and is never converted silently.
The binding checks every cap before it releases the GIL. A bad input raises `ValueError` and never reaches a kernel.
The GIL is released during the checker call. Calls from several threads are serialised by one mutex.

| Function | Arguments | Returns |
|---|---|---|
| `build()` | none | `"cuda"` or `"cpu"` |
| `device_count()` | none | CUDA devices. 0 on a CPU build or without a driver |
| `DEFECTS` | constant | Tuple of the 13 defect names, in index order |
| `check_polygons(xy, ring_start, shape_start, roles, W, H, G)` | see below | dict: `word_start` (n+1), `evenodd` (uint32 words), `defects` (n, 13) |
| `check_masks(bits, hw, start)` | `bits` uint8 (all masks, row-major), `hw` int32 (n, 2) as (H, W), `start` int64 (n+1) | dict: `quads` (n, 5), `ring_start`, `pts_start`, `ring_area`, `pts` |
| `add_tables(raters, uni, core, vote, dist, rater_area, voter, sig_start, sig, n_vote, perimeter, t5, min_area)` | One `ExpertCache` | int handle |
| `drop_tables(handle)` | handle | none. `ValueError` for an unknown handle |
| `resident()` | none | Number of tables on the device |
| `score_polygons(xy, ring_start, shape_start, roles, W, H, G, tables)` | polygons plus one handle per object (int64) | dict: `defects`, `area`, `iou_tc`, `f5`, `g_area`, `topo_freq`, `sig_start`, `sig` |
| `score_masks(bits, hw, start, tables)` | masks plus handles | the same, `defects` all 0 |

### Packed polygon input

| Array | dtype and shape | Meaning |
|---|---|---|
| `xy` | int32 (P, 2) | Vertices of all rings, `x0, y0, x1, y1, ...` |
| `ring_start` | int64 (R+1) | Ring r is `xy[ring_start[r] : ring_start[r+1]]` |
| `shape_start` | int64 (n+1) | Shape s is rings `shape_start[s] : shape_start[s+1]`. At most 10 rings |
| `roles` | int32 (R) or `None` | -1 for an outer ring. p for a hole of outer ring p. `None` means all outer |
| `W`, `H` | int32 (n) | Frame size of each shape, 1 to 4096 |
| `G` | int | Grid. 1000 in production. 1 to 2^30 |

`objdel.geometry.shape.pack_objects` builds these arrays from parts.
The mask of shape s is `evenodd[word_start[s] : word_start[s+1]]`: H rows of ceil(W / 32) words. Pixel x is bit `x % 32` (LSB first) of word `x // 32`. Bits past W are 0 (`unpack_rows`).

### Input caps

| Cap | Value | Where checked |
|---|---|---|
| Frame side W, H | 1 to 4,096 | `check_dims` |
| Rings per shape | at most 10 (`kMaxRings`) | `check_offsets` |
| Vertices per ring | at most 1,000 (`kMaxRingVertices`) | `check_offsets` |
| Grid G | 1 to 2^30 | `check_grid` |
| Raters per table | 2 to 32 | `add_tables` binding |
| Table `dist` values | exactly 0, or in [1, 2^24] | `add_tables` binding. The exactness argument needs it |
| Table `perimeter`, `t5` | `perimeter >= 1`, `t5 > 0` and finite | `add_tables` binding |

### The 13 defect classes

Counts are per shape. Crossing, touching and overlap classes count edge pairs. Ring-ring classes count once per pair.

| # | Name | Meaning |
|---|---|---|
| 0 | `fewer_than_3_points` | Fewer than 3 distinct vertices after dropping repeats |
| 1 | `zero_area` | Shoelace sum is 0 |
| 2 | `self_crossing` | Non-adjacent edges of one ring cross properly |
| 3 | `vertex_pinch` | Non-adjacent edges touch at a point |
| 4 | `collinear_overlap` | Non-adjacent edges overlap along a segment |
| 5 | `spike` | Adjacent edges fold back on each other |
| 6 | `ring_ring_crossing` | Edges of two rings cross properly |
| 7 | `ring_ring_touching` | Edges of two rings touch or overlap |
| 8 | `hole_outside_outer` | A hole has a witness point strictly outside its outer ring |
| 9 | `nested_holes` | Sibling holes overlap or nest |
| 10 | `nested_outers` | Outer parts overlap or nest |
| 11 | `wrong_winding` | Outer with net area < 0, or hole with net area > 0 |
| 12 | `out_of_range` | Vertex with x or y outside [0, G] |

An answer is legal when all 13 counts are 0 for every object.

## 3. Kernels

Every algorithm is a `__host__ __device__` functor run by `for_each` over independent items.
On the CPU build `for_each` is a loop. On the GPU build it is a kernel with 128 threads per block.

| Stage | Items | Functor or kernel | Notes |
|---|---|---|---|
| Precision choice | one per shape | `Precision` | Picks int64 or `__int128` for raster and for validate from fixed bounds |
| Raster | one pixel row per warp | `raster_kernel` | Window of 1,024 columns. Edges crossing the row centre add their sign to the winding. A suffix sum gives the winding. Lanes build each 32-pixel word with `__ballot_sync`. Lane j keeps word j, so a window leaves in one coalesced store |
| Clean rings | one per ring | `CleanRings` | Drop repeats, count area and range |
| Edge pairs | one per edge | `EdgePairs` | Classify the edge against every later edge of its shape |
| Shape totals | one per shape | `ShapeTotals` | Ring-ring counts, winding, containment by exact witness points |
| Mask slots | one per boundary-edge slot | `InitSlots`, `Pred`, `LabelRound`, `RankInit`, `RankRound` | Slot = (vertex, direction). Successor from the 2x2 bit quad, left turn at saddles. Ring label = minimum slot by pointer jumping. Corner rank by list ranking |
| Mask output | per ring | `RingStats`, `ScatterCorners`, `GatherRings` | Rings in `polygonize` order |
| Row boxes | one per frame row | `RowBox` | First and last foreground column, count |
| Crop labelling | one per crop row | `CropFill`, `LabelInit`, `LabelMerge`, `LabelFlatten` | Union-find on the object's bounding box plus one pixel. 8-connected foreground, 4-connected background |
| Size filter | per crop | `Filter`, `Holes`, `Gather` | Drop parts below `min_area`, fill holes below `min_area`, relabel, count each hole for the part above its first raster pixel |
| Terms | one per object | `ScoreRows`, `ObjTerms`, `SigCompact` | IoU with every rater and the vote, band depth, area guard, signature |

A scoring call has three host round trips: the crop layout (5 numbers per object), the terms and signature lengths, then the signatures.

## 4. Memory model

```text
  Engine (one per process, never destroyed)
  mutex
  Lane 0: stream, dev arena, score arena, in (pinned), out (pinned)
  Lane 1: stream, dev arena, score arena, in (pinned), out (pinned)

  Registry: resident expert tables (one cudaMalloc each, handle = index, reused after drop)
```

- **Arena.** A `Block` is a grow-only buffer. `get(bytes)` keeps it when it is large enough, else replaces it by one of the next power of two up to 64 MiB, and exact above 64 MiB. Nothing is allocated per call. Arrays of one call sit at 256-byte aligned offsets (`Layout`).
- **Pinned staging.** The `in` and `out` buffers per lane are page-locked host memory (`cudaMallocHost`). A polygon call packs its input into the pinned buffer, does one copy in, runs kernels on the lane's stream, does one copy out, then synchronises.
- **Streams.** Two non-blocking streams at the lowest priority, so the trainer's kernels, queued at the same time, go first. A polygon batch of at least 64 shapes (`kSplitMin`) is cut in two halves on the two lanes. The copy and the host packing of one half overlap the kernels of the other. Scoring and mask calls use lane 0.
- **Drain.** A guard waits for both lanes on every exit of a call, also on an exception. So no copy stays in flight from a staging buffer that the next call refills.
- **Error state.** Each call clears the CUDA error state first. A failed call, for example an out-of-memory error, does not poison the next one.
- **Tables.** `add_tables` packs each plane into 32-pixel words (LSB leftmost) and copies them with one `cudaMalloc`. A 50 x 50 LIDC object takes about 13 KB, a 100 x 100 one about 51 KB, and all 7,341 LIDC train items about 140 MB (`docs/research/notes/gpu_scoring.md`).
- **Exit.** The engine and registry are never destroyed. At process exit the driver frees all memory. This avoids tear-down order problems with the CUDA runtime.
- **Forks.** `ExpertCache.device()` records the process id. A forked child uploads its own tables and never frees its parent's.

## 5. Exactness guarantees

| Quantity | Guarantee | How |
|---|---|---|
| Raster and defects | Equal to the exact oracles and to polycore | Integer arithmetic. `int64` or `__int128` by fixed bounds. Half-open crossing rule |
| Area, intersections, counts | Integers | No floating point |
| IoU | `inter / union` in float64. `1.0` when the union is 0 | Same operation as the host formula |
| Band depth D | Equal to any float64 summation order while the sum is below 2^30 | Distances are float32, 0 or at least 1, so multiples of 2^-23. Each row sums as an int64 count of 2^-23, rows add in order in int128 |
| F5, g_area, shares | Same IEEE operations in the same order as the Python code | No FMA: `-ffp-contract=off` for C++, `--fmad=false` for nvcc |
| Signature | Equal to `checker.topology.signature` | Same steps: open, close, relabel, count holes by first raster pixel, sort |

The ties of "max IoU over raters with the same signature" take the first rater of equal values, as Python's `max` does.

Measured equivalence:

| Set | Cases | Mismatches | Hardware | Source |
|---|---|---|---|---|
| Fuzz polygons | 8,559 | 0 | RTX 4090, RTX 5090 (sm_120), RTX 3090 (sm_86), Mac CPU | `experiments/cuda_checker/RESULTS.md`, `perf/p1_results.md`, `perf/k1_results.md` |
| Fuzz masks | 6,024 | 0 | RTX 4090, RTX 3090 | same |
| LIDC polygons | 144,133 | 0 | RTX 4090, RTX 5090 | `RESULTS.md`, `docs/research/notes/issues_solved_2026-11-07.md` |
| LIDC masks | 37,093 | 0 | RTX 4090, RTX 5090 | same |

## 6. Throughput (RTX 5090, `experiments/phase2/perf/p1_results.md`)

Polygons per second, validate plus both masks:

| Vertices | Frame | Batch | GPU checker | polycore loop (CPU) |
|---|---|---|---|---|
| 16 | 64 | 1 | 32,787 | 96,665 |
| 16 | 64 | 64 | 1,307,336 | 96,917 |
| 16 | 64 | 1,024 | **7,738,843** | 88,431 |
| 128 | 256 | 1,024 | 1,110,802 | 34,528 |
| 1,000 | 512 | 64 | 69,377 | 8,752 |
| 1,000 | 512 | 1,024 | 45,943 | 8,801 |

The first CUDA build (R0 v3, managed memory per call) ran 54,710 polygons per second at 16 vertices and batch 1,024. The arena, pinned staging and two streams gave 141x.

Known limit: at 1,000 vertices and batch 1,024 the GPU checker is below the 8-thread CPU path (80,839 per second), because the witness tests use one thread per shape (`p1_results.md`). Issue CK-4 in `experiments/phase2/perf/issues.md` lists related second-order costs, among them the O(E^2) edge pairs. Reward answers have far fewer vertices.

## 7. Build

| Item | Where | Notes |
|---|---|---|
| Project | `CMakeLists.txt` | `find_package(nanobind CONFIG REQUIRED)`. `nanobind_add_module(_gpucheck NB_STATIC NOMINSIZE STABLE_ABI ...)` |
| Core library | `add_library(gpucheck STATIC cuda/gpucheck.cpp)` | C++17. Static, so nanobind flags never reach nvcc |
| FMA | `-ffp-contract=off` (C++) and `--fmad=false` (CUDA) | Required for exactness |
| GPU switch | `-DOBJDEL_CUDA=ON` | Sets the file's language to CUDA, links `CUDA::cudart_static`. Architectures default to `native`. Override with `CMAKE_CUDA_ARCHITECTURES` |
| `NOMINSIZE` | nanobind option | The default `-Os` halves kernel speed |
| Wheel | `pyproject.toml` | scikit-build-core 1.1.1, nanobind 3.1.0. `wheel.py-api = "cp312"`: one `cp312-abi3` wheel per platform |
| Build directory | `build/<wheel tag>` | Persistent, so rebuilds are incremental |

Make targets (`Makefile`):

| Command | Effect |
|---|---|
| `make build` | Rebuild the extension if `CMakeLists.txt` or `cuda/` changed. Uses the stamp `build/.gpucheck-stamp` |
| `make cuda` | Reinstall with `OBJDEL_CUDA=ON` and check `build() == "cuda"`. `CUDA_ARCHS="80;90"` cross-compiles |
| `make wheels` | `cibuildwheel` for manylinux x86-64 and aarch64 in Docker |
| `bash cuda/docker_check.sh` | nvcc build check on the Mac in an amd64 container for sm_80, 86, 89, 90 and 120. With no GPU, calls must fail with a clean CUDA error |
| `bash experiments/cuda_checker/run_box.sh [--build-only]` | On a CUDA box: install nvcc 13.0.3 and compute-sanitizer 13.4.92 from PyPI, `make cuda`, GPU tests, equivalence suite, then the eval |

`run_box.sh` refuses to run when any RIGA file is under `data/`.
The module is linked against the CUDA 13.0 runtime. It needs only the driver at run time (driver 580 or newer).

## 8. Sanitizer status

| Tool | Scope | Result | Source |
|---|---|---|---|
| compute-sanitizer memcheck | 40 GPU tests, RTX 3090 | 0 errors | `experiments/phase2/perf/k1_results.md` |
| compute-sanitizer initcheck | same | 0 errors | same |
| memcheck and initcheck | 52 tests in 5 files, RTX 5090, full RL run chain | 0 errors | `docs/research/notes/full_rl_2026-11-08.md` |

Notes:

- The sanitizer runs with `PYTORCH_NO_CUDA_MEMORY_CACHING=1`, so torch's pool does not hide out-of-bounds accesses.
- Version 13.0.85 of the sanitizer wheel lacks its launcher executables. Pin `nvidia-cuda-sanitizer-api==13.4.92`, as `run_box.sh` does.
- `racecheck` and `synccheck` are not recorded in the files read.
- The `rl_full.yaml` chain runs memcheck and then initcheck as a soft step: a missing binary or any error fails the step.

## 9. How to change the engine safely

1. Write the algorithm as a functor over independent items. No shared mutable state except through the platform helpers (`atomic_add`, `atomic_min`, `popc`).
2. Keep integer arithmetic. If you add a float term, keep the operation order of the host formula and build without FMA.
3. Add the case to the oracle or the golden test first: see [Testing](../testing.md).
4. Run `make test`. The CPU build runs the full checker suite.
5. Run `make cuda` and the GPU group on a box: `tests/test_cuda_checker_gpu.py`, `tests/test_gpu_scoring.py`.
6. Run memcheck and initcheck on the box when you change memory use.
