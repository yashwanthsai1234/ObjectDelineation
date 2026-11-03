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

Not run yet. It runs as the last step of R0 (`experiments/phase2/configs/r0.yaml`): the data travel in the LIDC bundle (`checker_eval/`), the results come back in `runs/<run>/cuda_checker/`. Host: any Vast card with driver >= 580 (the Phase 2 image's floor; the module is linked against the CUDA
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
(not run yet)
<!-- END eval:cuda -->
