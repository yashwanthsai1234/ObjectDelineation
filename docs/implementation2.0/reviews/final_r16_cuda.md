# Final review round 16: CUDA and GPU performance (Implementation.md)

Reviewer role: CUDA and GPU performance. I checked the kernels, GPU placement, syncs, memory and the rental plan.
I read the doc by sections (front 5-6, C.3, C0, C1, C5, C.6, D.4-D.8, A3, A.9.6, B.9.5, B.9.8) and opened the code
and kernel files that it names.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Review only. I installed nothing, because an install needs the owner's approval. It sent me to NVIDIA's own pages for two checks: the CUDA 13.x driver floor (R580, minor version compatibility) and the compute-sanitizer release notes (no hard driver floor for 2026.x; see web checks). No item came from it. |
| engineering:code-review | I re-read `gpucheck.cpp` (`alloc`/`release` :78-84, `Block` :135-150, `Lane` :152-156, `Drain` :177-184, `check_polygons` :1663-1692, `score_words` :1549-1641, `score_polygons` :1784-1806, `add_tables` :1720-1763) against D.4, D5 and D5c. The D5c pinned-bound rule holds for RL: the reward calls only `score_polygons` [code: src/objdel/reward/score.py:115], and `check_polygons` (whose `out` holds the even-odd planes, :1679) runs only in `to_masks` (eval, targets). No item. |
| engineering:debug | I traced `ExpertCache.build` after D7: `rater_tau` and `f5_tolerance` call `edt` 6 + 6 times for 3 raters [code: src/objdel/reward/cache.py:173, 176-182, 188, 191-197], and D7 sends `edt` to the device. This gave item 1. I also traced `final_eval` with "pool size 1": the pool still forks one worker, which makes the engine calls. This gave item 3. |
| engineering:system-design | I walked the R-D2 box chain (`run_box.sh --gpu-check` steps 0-6) and its minutes in D.8 and front 6.3 V3 item 1. This gave items 1 and 4. |
| engineering:architecture | I checked the R-D1 path (`r_d1.yaml` running the K1 and K2 runbooks) against the box image. The image has no `nvcc` on `PATH`. This gave item 2. |
| engineering:documentation | Each item gives the source file, the heading, the doc quote, the code or measurement, and one exact fix. |
| huggingface-skills:huggingface-best | No model choice is in my area this round. Nothing changed. |
| huggingface-skills:huggingface-papers | No new paper claim is in my area. Nothing changed. |

Web checks (2026-11-09):
- docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html: CUDA 13.x needs driver 580 or later.
  The page does not mention developer tools.
- docs.nvidia.com/compute-sanitizer/ReleaseNotes: no general minimum driver for 2026.x. One feature needs >= 595
  (shared-memory coredumps) and one needs >= 590 (a Windows fix). So the 13.4.92 sanitizer on an R580 host is not a
  proven risk. K1's box had driver 615.71.09 [measured: runs/k1/box.txt]. No item.
- pypi.org `cuda-toolkit 13.0.3`: extra `nvtx` = `nvidia-nvtx==13.0.85.*`, extra `nvcc` = `nvidia-cuda-nvcc==13.0.88.*`,
  extra `cccl` = `nvidia-cuda-cccl==13.0.85.*`. These match D.4 "Common toolchain". The box lock already pins
  `nvidia-nvtx==13.0.85` [code: experiments/phase2/docker/locks/vllm_native.txt:479], so D1's extra `nvtx` installs
  nothing new and moves no torch pin.

Re-verified and found correct (no item):
- D5 step 1 line references: `label(1), label(0)` :1614, `zero(c.aux, ...)` :1615, `Holes` :1616, `bh + 2` crop rows
  :1583. The skip is safe: when no filter changes an object, the first-pass fg and bg labels are still exact.
- D5c `out` bound for `score_polygons`: defects `L.out.get(c.words - c.defects)` (:1799), boxes (:1572), `ObjOut` and
  `sig_start` (:1622-1623), signatures (:1638). All four are in the doc's list.
- D.4 toolchain and NVTX: FindCUDAToolkit (3.26-3.31) creates `CUDA::nvtx3` with no header check, so D1's explicit
  `find_file` is needed. The `cudart_static` link stays [code: CMakeLists.txt:34].
- D6 / K2 digests `b158d8132b18be83` and `f92b9a3939cec10f` [code: K2/run_cpu.sh:3; K2/GPU_RUN.md:20-21]. D-G1 counts
  `37093 masks, 38002 rings, 1476508 corners` [measured: K1/results/cross_numpy_m3.txt:7].
- C.3.3 sizes: 32 x 512 x 512 x 3 B = 25.2 MB; 32 x 1,024 x 1,536 x 4 B = 201 MB; 8 rows per image gives 1.61 GB.
- V1 money: items 53 + 119.2 + 10.4 + 6.9 + 3.5 = 193.0 min; $0.10 + 193.0 x $0.01133 = $2.29. Worst case
  $2.49 <= $19.90 - $16.90.
- R5 eval scoring time is small: `total - load_and_sample` = 0.69 s (val) and 3.11 s (test)
  [measured: runs/r5/eval/eval_*/final_eval.json]. So scoring with one worker in V1 does not break the 26-min eval
  tail.

Not re-raised (already in earlier logs or reviews): the pinned bound and its first-entry rule (r13 cuda-1), the
`arena_bytes` definitions (r15 cuda-1), sanitizer failure stopping steps 5-6 (r15 cuda-3, minor), the post-E4
`VllmSampler` memory values (r15 cuda-2), the `drop_tables` GIL (r15 cuda-4), the stale kernel evidence paths
(r12 item 3), and the pageable-copy sync claim (r14 item 3).

Counts: blocking 0, major 1, minor 3.

---

## 1. MAJOR. D8 puts all GPU cases of `tests/test_scorer_sums.py` into the sanitizer corpus. That includes G01, which builds 10,881 expert caches, and after D7 each build makes about 12 device `edt` calls. Four sanitizer tools, each with a 300 s cap, will likely time out, so D-G5 fails on the paid box with no defect

Sources: D (heading "D8", the `tests/test_gpu_sanitize.py` bullet and the `sanitize` bullet; heading "D.6", row
D-G5), B (heading "B.9.5", rows G01 and G07; heading "B.9.8", the `test_scorer_sums.py` sketch).

Doc:
- D8: "create `tests/test_gpu_sanitize.py` (markers `gpu` and `sanitizer`), a small fixed corpus so each tool
  finishes inside its cap: ... Once B3.1 lands, its `tests/test_scorer_sums.py` joins this corpus; its GPU cases
  carry `gpu` and `sanitizer`."
- D8 `sanitize`: "Each tool runs as `timeout 300 "$CUDA_HOME/bin/compute-sanitizer" --tool T ... "$PY" -I -m pytest
  -q $D_GPU_TESTS -m "gpu and sanitizer"`".
- D-G5: "memcheck, initcheck, synccheck, racecheck: `ERROR SUMMARY: 0 errors` on ... (`tests/test_gpu_sanitize.py`,
  plus `tests/test_scorer_sums.py` after B3.1), each tool inside `timeout 300`".
- B.9.5 G07: "`compute-sanitizer` (memcheck, racecheck) over the `score_masks` sums test (G01 inputs and G01's tables
  rule) ... runs in R-D2's `make gpu-check` sanitize step (D8 corpus, markers `gpu` and `sanitizer`)".
- B.9.8 sketch: the only GPU case is `@pytest.mark.gpu def test_g01_golden()`, which calls the whole `b_sets`. The
  sketch does not give it the `sanitizer` marker, but D8 tells the engineer to add it.

Code and measurement:
- `b_sets` scores 203 random masks, 9,678 polygons and 1,000 crops, and calls `b_tables` -> `ExpertCache.build` once
  per object. On the Mac it took 100.9 s, and the cache builds alone took 8.6 + 90.6 + 0.6 s [measured: B.9.5 G01,
  notes/r10qa1/run_sketch.py, probe2.py, probe_poly.py]. That was the scipy EDT, before D7.
- After D7, `cache.edt` calls `_gpucheck.edt` (D7 Files). With 3 raters, `ExpertCache.build` calls `edt` 6 times in
  `rater_tau` (3 pairs x 2 [code: src/objdel/reward/cache.py:173, 176-182]) and 6 times in `f5_tolerance` (3 x
  `_band`, 2 `edt` each [code: cache.py:188, 191-197]). Each call is a separate engine call with its own H2D copy,
  kernels, D2H copy and wait. The sketch also calls `c.device()` once per cache (`add_tables_from_raters`, one call
  per cache). So G01 makes about 10,881 x 13 = 141,000 engine calls [arithmetic] under each tool, plus at least
  100 s of host work [measured, Mac].
- compute-sanitizer instruments every launch and every API call. racecheck and initcheck are the slowest tools. At
  even 1 ms of tool cost per call, the device part alone is about 141 s [estimate]. Add the host part and pytest's
  start-up, and a 300 s cap is very likely to be hit [estimate].

Why major and not minor: no $0 check measures sanitizer time. The Mac has no compute-sanitizer, and D8's `--dry` mode
only runs `--collect-only`. A timeout exits 124, D-G5 fails, and with today's `set -euo pipefail` (r15 cuda-3) the
D-G7/D-G8 benches and `d_gpu_test.yaml` step 2 do not run either. That is paid R-D2 time with no gate result. The doc
also contradicts itself: "a small fixed corpus so each tool finishes inside its cap" versus "its GPU cases carry
`gpu` and `sanitizer`", where the one GPU case is G01.

Fix:
- D8 (`tests/test_gpu_sanitize.py` bullet): replace "Once B3.1 lands, its `tests/test_scorer_sums.py` joins this
  corpus; its GPU cases carry `gpu` and `sanitizer`." with "Once B3.1 lands, `tests/test_scorer_sums.py` adds one
  test, `test_g07_sanitize_sums` (markers `gpu` and `sanitizer`): the three fixed 64 x 64 masks plus the first 29
  masks of `random_masks(np.random.default_rng(1))`, the first 64 polygons of each G from
  `compare.fuzz_polygons(530, seed=0)`, and the first 32 crops of `checker_eval/masks.npz` (through `box_input`),
  each with its `b_tables` handle. It checks only `sum_x` and `sum_y` against `np.nonzero` (no golden). `test_g01_golden`
  keeps the marker `gpu` only and runs in step 3, never under a sanitizer."
- D-G5: "plus `tests/test_scorer_sums.py::test_g07_sanitize_sums` after B3.1".
- B.9.5 G07: "over `test_g07_sanitize_sums` (a G01 subset: 32 masks, 3 x 64 polygons, 32 crops; G01's tables rule)".
- B.9.8 sketch: add the `test_g07_sanitize_sums` function with both markers.

## 2. MINOR. `r_d1.yaml` runs the K1 runbook steps 1-5, but K1's step 1 leaves the toolkit install as comments and step 4 calls a bare `compute-sanitizer`. The box image has neither `nvcc` nor the sanitizer on `PATH`

Sources: D (heading "D1", step 0, bullets "Create `experiments/phase2/configs/r_d1.yaml`" and "Add a preflight
check"; heading "D8", "Why one shell entry"), 00_front (heading "6.3", V2 item 4).

Doc:
- D1 step 0: "Create `experiments/phase2/configs/r_d1.yaml` (like `checker_gpu.yaml`). It runs K1-draft and K2-draft
  runbook steps 1-5, then the bench if D1 has merged, through `run.sh` ... No manual `scp`, no toolchain install by
  hand."
- D8: "the box image has no nvcc. Only `run_box.sh` installs the toolkit and exports `CUDACXX`, and only inside its
  own shell. Each `steps:` entry of a run config is its own `bash -c`".

Code:
- K1 step 1 runs only `nvcc --version || true`. The install lines are comments: `#   pip install
  "cuda-toolkit[nvcc,cccl]==13.0.3"`, `#   export CUDA_HOME=...`, `export NVCC=...` [code:
  docs/implementation2.0/kernels/K1/GPU_RUN.md:43-47].
- K1 step 2 builds with `NVCC=${NVCC:-nvcc}` [code: K1/build_gpu.sh:12]. K1 step 4 calls `timeout 300
  compute-sanitizer --tool memcheck ...` with no path [code: K1/GPU_RUN.md:87-88].
- The image sets `ENV PATH=/opt/venv/bin:$PATH` [code: experiments/phase2/docker/Dockerfile.vllm_native:31]. The pip
  `nvcc` and `compute-sanitizer` sit under `.../site-packages/nvidia/cu13/bin` (checker_gpu.yaml finds the sanitizer
  there by its full path [code: experiments/phase2/configs/checker_gpu.yaml:16]).
- K2's runbook (GPU_RUN.md:41-49) and K3's (GPU_RUN.md step 1) install the toolkit and export the paths. K1's does not.

So a literal port of K1 steps 1-5 fails at `build_gpu.sh`, and the `&&` chain [code: policy.py:490] then skips K2 and
the bench. Why minor: R-D1 is optional, and front 6.3 says "on the time-line order it is dropped". A careful engineer
could also copy the pattern of `checker_gpu.yaml`.

Fix (D1 step 0, new bullet after "Create `experiments/phase2/configs/r_d1.yaml`"): "`r_d1.yaml` step 0 is
`bash experiments/cuda_checker/run_box.sh --build-only` (it installs `cuda-toolkit[nvcc,cccl]==13.0.3` and the
13.4.92 sanitizer into `/opt/venv`). Each later runbook step is one `bash -c` that starts with
`C=/opt/venv/lib/python3.12/site-packages/nvidia/cu13; export CUDA_HOME=$C NVCC=$C/bin/nvcc PATH=$C/bin:$PATH`.
In K1 GPU_RUN.md, step 1 loses its commented install and step 4 calls `"$CUDA_HOME/bin/compute-sanitizer"`." Add
`compute-sanitizer --tool` without `$CUDA_HOME` or `$SAN` to the preflight grep of D1 step 0.

## 3. MINOR. "Pool size 1" does not mean "no forked worker opens its own CUDA context". The one worker is a forked child, and it makes the engine calls

Sources: 00_front (heading "6.2", "Why step 0"), A (heading "A3", the `final_eval` bullet "on a CUDA build ... caps
the pool at 1 worker"; heading "A2", step 2).

Doc:
- Front 6.2: "On a CUDA build `final_eval` uses pool size 1, as A2 step 2 does, so no forked worker opens its own
  CUDA context (A3 adds the cap)."
- A3: "`final_eval` caps the pool at 1 worker ..., as A2 step 2 does, so no forked worker opens its own CUDA
  context."

Code:
- `fork_pool(n)` forks n workers with `get_context("fork")` [code: src/objdel/train/serve.py:164-177].
  `final_eval` then scores on that pool [code: experiments/phase2/final_eval.py:69, 77], through
  `pool.map(_score, ...)` [code: experiments/phase2/select_ckpt.py:230].
- In the worker, `metrics()` calls `to_masks` -> `_gpucheck.check_polygons` [code: src/objdel/geometry/shape.py:143]
  and `score_many` -> `cache.device()` -> `_gpucheck.add_tables` [code: src/objdel/reward/cache.py:250]. On a CUDA
  build this opens a CUDA context in the child.

This is not a correctness bug. The parent has no CUDA context at the fork (`fork_pool` refuses after torch has
initialized CUDA [code: serve.py:167-169], and cache construction is host-only today). So the child's context is
legal, and it costs one more context beside vLLM and the trainer. But the stated reason is false, and D-G17 (which
counts contexts) applies only after D7 step 6. An engineer who reads "no forked worker opens a context" may wrongly
expect a single context in V1's `nvidia-smi`.

Fix (front 6.2 and A3, same words): "On a CUDA build `final_eval` caps the pool at 1 worker, as A2 step 2 does. The
worker is forked before the model loads, so the parent has no CUDA context at the fork. The worker then opens one
CUDA context of its own at its first engine call (one extra context, not 16). D7 step 6 later moves this work into
the main process (D-G17)."

## 4. MINOR. D.8 prices `pytest -m gpu` on `D_GPU_TESTS` at "about 5" minutes, but the list holds G01's full `b_sets` (100.9 s of Mac host work before D7) beside GPU02's corpora

Sources: D (heading "D.8", row R-D2), 00_front (heading "6.3", V3 item 1).

Doc: D.8 R-D2: "`pytest -m gpu` on D's box test list (`D_GPU_TESTS`, D8) about 5". `D_GPU_TESTS` includes
`tests/test_scorer_sums.py` once B3.1 lands (D8).

Measurement: `b_sets` took 100.9 s on the Mac CPU build [measured: B.9.5 G01]. After D7 its 10,881 cache builds each
add about 12 device `edt` calls (item 1). GPU02 (in the same list) scores 37,093 crops, 144,133 LIDC polygons and
1,000 caches with EDT against scipy on the box (D.9.7 GPU02). So 5 minutes is low, and the V3 item 1 range
(27-46 min) may be short. Why minor: V3 item 1 is optional and owner-approved, and its stop rule is a skip.

Fix (D.8 R-D2 row and front 6.3 V3 item 1 basis): "`pytest -m gpu` on `D_GPU_TESTS` about 8 (G01's `b_sets` about 2
[measured: 100.9 s on the Mac]; GPU02 and the R3 corpus about 6 [estimate])". Then V3 item 1 becomes 30-49 min,
$0.34-0.56 at r_s (30 x $0.01133 and 49 x $0.01133; item 0 holds b), and the "Items 0, 1, 3 and 4" row and the V3 amount rule are recomputed with 49 min.
