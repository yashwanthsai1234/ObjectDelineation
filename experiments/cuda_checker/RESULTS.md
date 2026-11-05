# CUDA checker eval

What is checked: `objdel._gpucheck` (cuda/gpucheck.cpp, one source for the CPU and the GPU) against the CPU reference
checker, bit for bit, on every output the reward and certify use: the legality verdict and the 13 defect counts of
polycore `validate`, the even-odd and nonzero masks of polycore `rasterize`, and for masks the bit-quad counts and the
rings of `polygonize` (Betti numbers b0/b1 and the ring hash, SHA-256 of the rings). Command and case sets:
`experiments/cuda_checker/eval.py` (docstring). Data: `data/checker_eval/` from `prepare_data.py`, LIDC only.

Each run rewrites its own section below (`<build>` = cpu or cuda); the JSON is in `results/<build>.json`, the
smallest reproducing case of every mismatch in `data/checker_eval/repro/<build>/`.

## CPU run (Mac, CPU build)

<!-- BEGIN eval:cpu -->
Run 2026-11-02 21:55, commit f8aa4ef, build **cpu**, Darwin arm64, 8 CPUs; wall time 213 s.

| Case set | Cases | Mismatches | Exact match | First differing field (count) |
|---|---|---|---|---|
| (a) fuzz polygons | 8,559 | 0 | 100.000 % | none |
| (a) fuzz masks | 6,024 | 0 | 100.000 % | none |
| (b) LIDC polygons | 144,133 | 0 | 100.000 % | none |
| (c) LIDC Phase 0 masks | 37,093 | 0 | 100.000 % | none |

LIDC polygons by reference DEFECTS class (a case counts once per class it has):

| Class | Cases | Mismatches |
|---|---|---|
| legal | 138,257 | 0 |
| fewer_than_3_points | 4 | 0 |
| zero_area | 41 | 0 |
| self_crossing | 3,843 | 0 |
| vertex_pinch | 2,584 | 0 |
| collinear_overlap | 667 | 0 |
| spike | 1,249 | 0 |
| ring_ring_crossing | 861 | 0 |
| ring_ring_touching | 584 | 0 |
| hole_outside_outer | 284 | 0 |
| nested_holes | 4 | 0 |
| nested_outers | 723 | 0 |
| wrong_winding | 952 | 0 |

Polygons per second (validate + both masks; polycore_loop = the reward's per-shape path):

| Vertices | Crop | Batch | gpucheck | polycore loop | score_batch 1 thread | score_batch 8 threads |
|---|---|---|---|---|---|---|
| 16 | 64 | 1 | 100,847 | 98,561 | 110,339 | 114,286 |
| 16 | 64 | 64 | 62,004 | 54,773 | 277,407 | 157,514 |
| 16 | 64 | 1024 | 66,849 | 30,793 | 143,762 | 195,642 |
| 128 | 256 | 1 | 5,230 | 30,284 | 52,288 | 37,617 |
| 128 | 256 | 64 | 3,862 | 15,966 | 30,377 | 45,033 |
| 128 | 256 | 1024 | 3,334 | 4,917 | 12,173 | 109,624 |
| 1000 | 512 | 1 | 351 | 16,461 | 17,241 | 17,518 |
| 1000 | 512 | 64 | 162 | 5,521 | 7,400 | 24,136 |
| 1000 | 512 | 1024 | 137 | 2,853 | 4,571 | 36,609 |

Masks per second (bit-quads + rings; python_certify = bitquad + betti + polygonize; OpenCV = findContours + connectedComponents, a reference point, not the same output):

| Masks | Shape | Batch | gpucheck | Python certify path | OpenCV |
|---|---|---|---|---|---|
| lidc | 22x21 | 256 | 3,108 | 1,711 | 41,463 |
| lidc | 22x21 | 8192 | 2,006 | 2,620 | 54,307 |
| noise512 | 512x512 | 4 | 5 | 1 | 9 |
<!-- END eval:cpu -->

## GPU run (CUDA box, `make cuda` build)

Ran 2026-11-03 in R0 `v3` (Vast box 58939372, RTX 4090, British Columbia). Score and diagnosis below the table.
It runs as the last step of R0 (`experiments/phase2/configs/r0.yaml`): the data travel in the LIDC bundle (`checker_eval/`), the results come back in `runs/<run>/cuda_checker/`. Host: any Vast card with driver >= 580 (the Phase 2 image's floor; the module is linked against the CUDA
13.0 runtime), US/CA only. Upload the repo snapshot and `data/checker_eval/` only (19 MB, LIDC; never RIGA: run_box.sh
and eval.py refuse it). Then, from the repo root in the Phase 2 image:

```
bash experiments/cuda_checker/run_box.sh
```

It adds nvcc (`cuda-toolkit[nvcc,cccl]==13.0.3`, under the image's lock) and pytest/hypothesis to /opt/venv, runs
`make cuda` (native arch), the GPU test group and the equivalence tests, then this eval. Expected wall time ~10 min
(downloads ~2, build ~2, tests ~2, eval ~4; the Mac CPU eval took 213 s). Copy back `results/cuda.json`, this file,
and `data/checker_eval/repro/cuda/` if any mismatch is reported. `cuda/docker_check.sh` already ran the build half
on the Mac (amd64 container, no GPU) for sm_80/86/89/90/120.

Reading the CPU numbers above: the CPU build of gpucheck is the test path of the GPU code (one loop iteration per GPU
thread, all-pairs edge tests, no sort), not a production CPU checker; polycore stays the CPU checker. The Mac was
shared with other jobs during the run, so CPU timings are indicative only.

<!-- BEGIN eval:cuda -->
Run 2026-11-03 09:09, code `999969b` + the CMake fix below, build **cuda** on NVIDIA GeForce RTX 4090, 615.71.09, Linux x86_64, 64 CPUs; wall time 73 s.

| Case set | Cases | Mismatches | Exact match | First differing field (count) |
|---|---|---|---|---|
| (a) fuzz polygons | 8,559 | 0 | 100.000 % | none |
| (a) fuzz masks | 6,024 | 0 | 100.000 % | none |
| (b) LIDC polygons | 144,133 | 0 | 100.000 % | none |
| (c) LIDC Phase 0 masks | 37,093 | 0 | 100.000 % | none |

LIDC polygons by reference DEFECTS class (a case counts once per class it has):

| Class | Cases | Mismatches |
|---|---|---|
| legal | 138,257 | 0 |
| fewer_than_3_points | 4 | 0 |
| zero_area | 41 | 0 |
| self_crossing | 3,843 | 0 |
| vertex_pinch | 2,584 | 0 |
| collinear_overlap | 667 | 0 |
| spike | 1,249 | 0 |
| ring_ring_crossing | 861 | 0 |
| ring_ring_touching | 584 | 0 |
| hole_outside_outer | 284 | 0 |
| nested_holes | 4 | 0 |
| nested_outers | 723 | 0 |
| wrong_winding | 952 | 0 |

Polygons per second (validate + both masks; polycore_loop = the reward's per-shape path):

| Vertices | Crop | Batch | gpucheck | polycore loop | score_batch 1 thread | score_batch 8 threads |
|---|---|---|---|---|---|---|
| 16 | 64 | 1 | 134 | 54,558 | 69,365 | 69,197 |
| 16 | 64 | 64 | 6,917 | 51,462 | 169,100 | 110,755 |
| 16 | 64 | 1024 | 54,710 | 51,743 | 174,781 | 290,142 |
| 128 | 256 | 1 | 133 | 29,536 | 38,718 | 38,944 |
| 128 | 256 | 64 | 4,252 | 22,886 | 42,248 | 75,440 |
| 128 | 256 | 1024 | 4,244 | 23,016 | 42,870 | 145,895 |
| 1000 | 512 | 1 | 105 | 7,244 | 9,085 | 9,149 |
| 1000 | 512 | 64 | 1,409 | 6,421 | 8,457 | 32,127 |
| 1000 | 512 | 1024 | 1,209 | 6,505 | 8,597 | 45,125 |

Masks per second (bit-quads + rings; python_certify = bitquad + betti + polygonize; OpenCV = findContours + connectedComponents, a reference point, not the same output):

| Masks | Shape | Batch | gpucheck | Python certify path | OpenCV |
|---|---|---|---|---|---|
| lidc | 22x21 | 256 | 4,043 | 2,580 | 57,348 |
| lidc | 22x21 | 8192 | 6,389 | 2,551 | 56,581 |
| noise512 | 512x512 | 4 | 19 | 3 | 22 |
<!-- END eval:cuda -->

### Score (RTX 4090, R0 v3)

**Correctness: pass.** 195,809 of 195,809 cases equal the CPU reference bit for bit: every fuzz category, every LIDC
source (train, eval and probe targets, 58,028 model outputs) and all 12 DEFECTS classes present in the data. No repro
case was written.

**Speed: fail.** The GPU build is slower than the CPU checker it should replace, at every polygon size and batch:

| Workload | gpucheck (GPU) | best CPU path on the same box | GPU / CPU |
|---|---|---|---|
| 16 vertices, batch 1 | 134 /s | 69,365 /s (score_batch, 1 thread) | 0.002 |
| 16 vertices, batch 1024 | 54,710 /s | 290,142 /s (score_batch, 8 threads) | 0.19 |
| 128 vertices, batch 1024 | 4,244 /s | 145,895 /s | 0.03 |
| 1000 vertices, batch 1024 | 1,209 /s | 45,125 /s | 0.03 |
| LIDC masks 22x21, batch 8192 | 6,389 /s | 2,551 /s (Python certify path) | 2.5 |
| noise 512x512 masks, batch 4 | 19 /s | 3 /s (Python certify path) | 6.3 |

The masks beat the Python certify path, but not OpenCV (56,581 /s; a reference point, not the same output).

### Where it goes wrong

Read from the timings and `cuda/gpucheck.cpp`; not yet profiled with Nsight Systems (the next box step).

1. **A fixed cost of about 7.5 to 9.5 ms per call.** One polygon takes 7.5 ms at 16 vertices and 9.5 ms at 1000
   vertices; a batch of 64 small polygons takes 9.3 ms. The cost hardly depends on the work. Every call creates each
   `Buf` with `cudaMallocManaged` + `cudaMemset` + `cudaDeviceSynchronize`, copies the input into managed memory on
   the host (the pages then migrate to the GPU on first touch), reads the results back through page faults, and frees
   everything. Fix: allocate device buffers once and reuse them, copy with `cudaMemcpyAsync` from pinned host memory
   on one stream, and synchronize once per call.
2. **Brute-force loops in two steps.** `RasterRows` runs one thread per pixel row, and each row thread reads every
   edge of its shape: H x E tests per shape (512 x 1000 = 5 x 10^5). `EdgePairs` runs one thread per edge, and each
   thread tests every later edge of its shape: about E^2 / 2 tests (5 x 10^5 at 1000 vertices). polycore sorts the
   edges and sweeps, so it skips most of these tests. A batch of 1024 shapes of 1000 vertices takes 847 ms: about
   10^9 tests, each reading its edge from global memory (no shared-memory reuse). The raster thread also keeps a
   256-entry array, which the compiler puts in local (off-chip) memory. Fix: one block per shape with its edges in
   shared memory; for the raster, each edge writes its row crossings and the rows sort them; for the pairs, sort the
   edges by x-min and test only the overlapping window.
3. **Masks: 16 slot-sized buffers per call and log-depth pointer jumping.** Each mask allocates 16 arrays of
   4 (W + 1) (H + 1) entries (about 1 M per 512 x 512 mask) through the same managed path, then runs `jump_rounds`
   label and rank rounds, one kernel launch each. Item 1's fix applies here as well.

Item 1 alone explains the small-batch numbers (a batch of 16 edges also fills 16 of the 4090's 16,384 cores); items 1 and 2 together explain the large-polygon numbers. Fix order:
1, then 2, then rerun this eval on the next box (correctness must stay at zero mismatches).

### Box defect found and fixed in this run

`make cuda` failed first on the box: "CUDA_ARCHITECTURES is empty for target gpucheck". `CMakeLists.txt` set the
default `CMAKE_CUDA_ARCHITECTURES native` after `add_library(gpucheck)`, and a target copies that variable only when
it is created. `cuda/docker_check.sh` passed because it gives the archs on the command line. Fixed by setting
`CUDA_ARCHITECTURES` on the target; the rerun on the same box built for the 4090 and passed 21 tests (2 skipped).
