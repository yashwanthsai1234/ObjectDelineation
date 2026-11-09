# Final review, round 3: CUDA and GPU performance (whole Implementation.md)

Reviewer role: CUDA and GPU performance (NVIDIA docs, Nsight, Blackwell sm_120, CCCL). Target:
`REPO/docs/implementation2.0/Implementation.md` (14,506 lines, "Revision: final round 2"), plus the code, kernel and
runbook files it names. Earlier files read: `reviews/final_r1_cuda.md`, `reviews/final_r2_cuda.md`. Round-2 items
1-17, 19 and 20 are fixed in the doc, and this file does not repeat them. Item 18 is fixed on the D side only. Item
11 below repeats its B side.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new catalog skill. `compute-sanitizer-diagnose-cuda-bugs` is already named. The catalog has no Nsight Systems skill, so I read the NVIDIA Nsight Systems user guide and the CUPTI docs directly (item 4). Nothing installed. |
| engineering:code-review | Checked each GPU claim against the code: `cuda/gpucheck.cpp` (lanes, `alloc`, `Block`, `Drain`, `add_tables`, `score_words`, `ScoreRows`, line refs 418-1643), `CMakeLists.txt`, `Makefile`, `run_box.sh`, `docker_check.sh`, `checker_gpu.yaml`, `reward/cache.py`, `eval/metrics.py`, `serve.fork_pool`, TRL's `LLM(...)` call, and the K1/K3 build scripts. This found items 1, 2, 7, 8, 9, 10 and 12. |
| engineering:debug | Reproduce-first lens: which gates fail with no bug (items 1, 3, 4, 6), and which commands fail on the paid box but not in the Mac dry run (item 2). |
| engineering:system-design | Rental item order and data dependencies inside V2 (item 5). Host variance against fixed-number gates (item 3). |
| engineering:architecture | Dead stream lane against the owner's dead-code rule (item 9). One tested build recipe for K3 (item 8). |
| engineering:documentation | Stale numbers across sections (item 6), stale runbook facts (item 12), wrong line ref (item 10). |
| huggingface-skills:huggingface-best | No model choice in this review. TRL builds vLLM with no `enforce_eager` [code: .venv trl/generation/vllm_generation.py:356-372], so vLLM decodes inside CUDA graphs (input to item 4). |
| huggingface-skills:huggingface-papers | No CUDA claim here rests on an arXiv paper that needs a new check (D.2 CCL/EDT sources are journal papers). No item. |

Web and repo checks made for this review (2026-11-08):
- Nsight Systems user guide, `--cuda-graph-trace`: "Use `graph` to trace whole CUDA graphs without node activities ...
  is the default when available" [web: docs.nvidia.com/nsight-systems/UserGuide]. CUPTI `CUpti_ActivityGraphTrace`:
  "represents execution for a graph without giving visibility about the execution of its nodes"
  [web: docs.nvidia.com/cupti/api/structCUpti__ActivityGraphTrace.html]. A third-party nsys schema reader states that in
  graph mode "the inner kernels don't appear in the kernel table" (`CUPTI_ACTIVITY_KIND_KERNEL`) [web:
  docs.rs/veloq-nsys-query, inspect/graph.rs]. The Analysis Guide describes `cuda_gpu_trace` only as "CUDA kernels and
  memory operations" [web: docs.nvidia.com/nsight-systems/AnalysisGuide]. So graph rows in that report are not
  confirmed (item 4).
- `nvidia-smi` in a container: the process panel is empty while processes run, because the driver returns host PIDs
  [web: github.com/NVIDIA/nvidia-docker/issues/1460; forums.developer.nvidia.com/t/343361] (item 7).
- `cuda/docker_check.sh` copies a `cpp` folder. That folder was deleted in commit 9c5dd2e (2026-11-07). The script
  was last changed on 2026-11-06 (f68b32d) [measured: `git log`, read-only; `ls cpp` fails] (item 1).
- `gpucheck.cpp` has no use of `lanes[1]` [measured: grep] (item 9).
- K3's tested torch-extension build: `TORCH_CUDA_ARCH_LIST=12.0`, plus
  `-U__CUDA_NO_BFLOAT16_CONVERSIONS__ -U__CUDA_NO_HALF_CONVERSIONS__ -U__CUDA_NO_HALF_OPERATORS__`, and host
  `-ffp-contract=off` [code: docs/implementation2.0/kernels/K3/test_torch_gpu.py:41-45; K3/GPU_RUN.md step 4] (item 8).
- Recomputed: V1 195.1 min = $2.084 and 175.85 min = $1.879; money left $0.74 / $0.87; LoRA 17.4 M; touched
  tensors 1.409 B; pinned bound 256 x 1000 x 10 x 8 B = 20.48 MB; the B3.1 popcount column-sum formula is correct for
  LSB = leftmost. No item.

Counts: **blocking 0, major 5, minor 8.**

---

## Major

### 1. [major] D (D1 step 3, D.9 S06, D.1 item 10): `cuda/docker_check.sh` fails today, before it builds anything
- Problem: the $0 gate before R-D1 says "`cuda/docker_check.sh` passes with the NVTX ranges and stage events
  compiled". The S06 build gate says the same. The script runs
  `cp -r --parents CMakeLists.txt Makefile pyproject.toml README.md cpp cuda src experiments/cuda_checker/run_box.sh /work/`
  under `bash -euo pipefail` [code: cuda/docker_check.sh:17]. The repo has no `cpp` folder: it was archived in
  9c5dd2e, one day after the script's last change [measured: git log; `ls cpp`]. So `cp` exits 1, and the gate fails
  with no D bug. Round 2 fixed the same stale path in D-G7's `git archive`, but not here.
- Fix: in D1's file list, add "change `cuda/docker_check.sh`: drop `cpp` from the `cp -r --parents` list". In S06,
  add "the copy list names only paths that exist (`git ls-files` has a match for each one)". This is one line, and
  it runs on the Mac at $0.

### 2. [major] D (D8 Makefile `gpu-check`, D-G7, D8 `sanitize`): `make gpu-check` on the box has no nvcc, no `CUDACXX` and the wrong `PY`
- Problem: D8 defines `gpu-check` = `make cuda` + `pytest -m gpu` + `sanitize` + benches. `make cuda` runs CMake with
  `OBJDEL_CUDA=ON`, and CMake finds nvcc only through `$CUDACXX` [code: Makefile `cuda` target comment]. The box image
  has no nvcc. Only `run_box.sh` installs `cuda-toolkit` and exports `CUDACXX`, and only inside its own shell
  [code: run_box.sh:24-29]. A later YAML step or a `make` call does not inherit it (each `checker_gpu.yaml` step is
  its own `bash -c`). The Makefile's `PY` defaults to `.venv/bin/python`, but the box venv is `/opt/venv`. D-G7 says
  "with the `CUDACXX` that `run_box.sh` exports", but nothing exports it to `gpu-check`. D8's `sanitize` calls bare
  `compute-sanitizer`, which is not on `PATH`: today's config uses the full path under `nvidia/cu13/bin`
  [code: checker_gpu.yaml:15-16]. `make gpu-check DRY=1` only prints the chain on the Mac, so none of this shows
  before R-D2 is paid for.
- Fix: in D8, write `gpu-check` as one shell entry: `bash experiments/cuda_checker/run_box.sh --gpu-check`. Add a new
  mode to `run_box.sh` that does all of this in the shell where `CUDA_HOME` and `CUDACXX` are set: (1) the toolkit
  install; (2) `make cuda PY=$PY`; (3) `pytest -m gpu`; (4) `sanitize`, which uses `"$CUDA_HOME/bin/compute-sanitizer"`;
  (5) the old-engine build `make -C <dir> cuda PY=/opt/old/bin/python` with the same `CUDACXX`; (6) the benches. In
  `d_gpu_test.yaml`, call that mode with `PY=/opt/venv/bin/python`. Add to S02: under `DRY=1`, the printed chain
  starts with `run_box.sh`, and every `make cuda` line carries `PY=` and runs after `CUDACXX=` is set.

### 3. [major] C (C.6 G-SP1, C.8.1) and front (6.1 "Step time with C code"): the adoption gate is a fixed number, while the baseline arm is measured on the same box
- Problem: G-SP1 passes at "mean step time <= 20.0 s (mode A) or <= 20.7 s (mode B)". These values come from R5's
  23.43 s minus the estimated savings, and "the values are fixed before the rental". The rental also measures the
  C0 baseline arm (today's code, `--shared-prefix`, same host, steps 3-15). Both arms also run under `nsys profile`
  with `OBJDEL_NVTX=1` (A6 step 8), which R5 did not. The step holds a 3.0-3.2 s CPU phase, and R4/R5 hosts already
  differed. So host speed and profiler injection move the absolute value. On a host whose baseline is 25 s, a fast arm
  that saves the full 3.4 s still fails. On a faster host, a fast arm that saves 1 s passes. Gate C-adopt depends on
  this gate alone among the speed gates.
- Fix: make G-SP1 relative to the same-rental baseline arm. Pass: `fast_mean <= base_mean - 3.4 s` (mode A) or
  `<= base_mean - 2.7 s` (mode B). These are the C.8.1 midpoints, fixed before the rental. Report the absolute means
  and the ratio beside them. If `base_mean` is outside 23.43 s ± 10 %, record "host not comparable to R5". The
  relative rule still decides. The re-profile line becomes `base_mean - 2.25 s` (A) or `- 1.55 s` (B). The money
  rows (front 6 "Step time with C code") keep 20.0 s only as a planning value. Same edit in C.8.1's last sentence and
  in the C.6 protocol paragraph.

### 4. [major] A (A6 step 9 `gpu_busy`) and C (C0 step 4 `gpu_active.py`, G-SP6): the GPU-busy share omits vLLM's decode, which runs in CUDA graphs that the capture traces at graph level
- Problem: TRL builds vLLM with no `enforce_eager` [code: trl vllm_generation.py:356-372], so decode runs as CUDA
  graph launches. A6 keeps nsys's default `--cuda-graph-trace=graph` and says the graphs "show up as single graph
  launches". In graph mode, CUPTI records each graph run as a `CUPTI_ACTIVITY_KIND_GRAPH_TRACE` record, and the inner
  kernels are not in the kernel table [web: CUPTI GraphTrace doc; veloq-nsys-query]. `gpu_busy` (A6) and
  `gpu_active.py` (C0) take the union of "kernel + memcpy + memset intervals from `cuda_gpu_trace`". NVIDIA's Analysis
  Guide describes `cuda_gpu_trace` as "kernels and memory operations", and I could not confirm that it lists graph
  runs [web: Analysis Guide]. If it does not, the decode part of `generate` (about 5.6 s + 4.3 ms per token of the
  longest answer, C.1) counts as idle. Then the A6 baseline idle share is too high, and G-SP6 (>= 70 %) fails with no
  bug. Also, A6's `host_syncs` counts only runtime API names (`cudaStreamSynchronize`, ...). vLLM and Triton call the
  driver API (`cuStreamSynchronize`, `cuCtxSynchronize`, `cuEventSynchronize`, `cuMemcpyDtoH*`), so the vLLM-range
  count misses those calls.
- Fix: in A6 step 9 and C0 step 4, build the busy intervals from the exported SQLite (A6 already writes
  `$NAME.sqlite`): the union of `CUPTI_ACTIVITY_KIND_KERNEL`, `..._MEMCPY`, `..._MEMSET` and
  `CUPTI_ACTIVITY_KIND_GRAPH_TRACE` (`start`, `end`), when the table exists. Add one graph row to the
  `tests/goldens/nsys_small/` fixture, plus a test that it counts as busy. On the box, once, compare the result with
  `nsys stats --report cuda_gpu_trace` and record whether graph rows appear there. Add the driver-API names to
  `host_syncs` after the suffix strip. Report trainer and vLLM counts apart, as C0 already does.

### 5. [major] C (C1 box test, C3 step 0, G-EQ1, G-EQ4) and front (6.3 items 2-3): `c_speed_check` needs "recorded" batches that V2 creates only in the next item
- Problem: G-EQ1 uses "all 32 images of 5 recorded steps", and G-EQ4 and C3 step 0 use "one recorded generation
  batch", all inside `c_speed_check` (V2 item 2). The only box source of recorded batches is
  `--record-batches` in the fast arm's check steps 18-32, which is V2 item 3. Item 3 runs after item 2, and it should
  run only if the G-EQ checks of item 2 pass ("C-G-EQ fails: C code is not used by B"). The Mac batches from A3's val
  answers (C2 tests) are the only other source, but no text says that `c_speed_check` reads them, what they contain,
  or how they reach the box. G-EQ4's gradient check also needs a loss with advantages, and the text gives no source
  for them.
- Fix: in C1/C3 and front 6.3 item 2, write: "`c_speed_check` reads 5 batches recorded on the Mac by the G-EQ2/3 Mac
  half (`runs/a3/batches/*.pt`, the `--record-batches` format: prompt ids, completion ids, masks, grids, `RowMeta`,
  image ids, and the group-centred advantages from A3's saved rewards). Their val images come from the
  `lidc-val-shifted` view in the bundle. Both ship in the V2 code and data sync." Add a Mac test: the files load with
  `torch.load(weights_only=True)`, the image hashes match the bundle, and the count is 5. The box batches from steps
  18-32 still feed G-EQ2/3's ">= 15 from the box".

---

## Minor

### 6. [minor] C (C.8.2): money-left figures and one sentence are stale against the front
- Problem: C.8.2 says "The money left after V1 is $0.59 (V1 at 27.5 s) or $0.79 (V1 at 23.65 s)" and "the front's own
  V2 table still omits the start-up rows (CUDA review round 2 item 3)". The front now gives $0.74 / $0.87 (6.3, with
  item 5 skipped at 27.5 s), and its V2 table has the start-up rows.
- Fix: replace the figures with "$0.74 (V1 at 27.5 s, item 5 skipped) or $0.87 (V1 at 23.65 s) (front 6.3)". Delete
  the "front's own V2 table still omits" clause. The conclusion ("does not fit") stays.

### 7. [minor] D (D-G17 (a), D.9 B06): the `nvidia-smi` row-count check can pass with no information
- Problem: in a container, the driver returns host PIDs, and `nvidia-smi` often shows an empty process list while
  processes run [web: nvidia-docker issue 1460; NVIDIA forum 343361]. If the query returns 0 rows before scoring and
  0 rows during it, (a) passes even when a new context appears.
- Fix: (a) first asserts at least 1 row before scoring. If the count is 0, (a) uses the container's own view: the
  number of processes in `/proc` that hold `/dev/nvidia[0-9]*` open
  (`for p in /proc/[0-9]*; do ls -l $p/fd 2>/dev/null | grep -q '/dev/nvidia[0-9]' && echo $p; done | wc -l`),
  sampled every 5 s, which must not rise during scoring. Record which method ran.

### 8. [minor] C (C1 step 7) and E (E.5): the K3 production build flags differ from the build that K3 tested
- Problem: C1 step 7 and E.5 prescribe `CUDAExtension` with nvcc flags `--fmad=false -arch=sm_120`. K3's tested build
  uses `TORCH_CUDA_ARCH_LIST=12.0` and no `-arch`. It also adds `-U__CUDA_NO_BFLOAT16_CONVERSIONS__`,
  `-U__CUDA_NO_HALF_CONVERSIONS__` and `-U__CUDA_NO_HALF_OPERATORS__` (torch defines these macros by default), and it
  builds host code with `-ffp-contract=off` [code: K3/test_torch_gpu.py:41-45]. Without the `-U` flags, the bf16 output
  path may not compile. `-arch` next to torch's own `-gencode` flags is a second arch source.
- Fix: in both places, copy the tested recipe: `TORCH_CUDA_ARCH_LIST=12.0`;
  `extra_compile_args={"cxx": ["-O3", "-ffp-contract=off"], "nvcc": ["-O3", "--fmad=false", "-lineinfo",
  "-U__CUDA_NO_BFLOAT16_CONVERSIONS__", "-U__CUDA_NO_HALF_CONVERSIONS__", "-U__CUDA_NO_HALF_OPERATORS__"]}`. Remove
  `-arch=sm_120`.

### 9. [minor] D (engine counters, D-G8, D9 ADR): `lanes[1]` is dead, yet `Drain` syncs it on every call
- Problem: `Engine` makes two lanes, and every entry uses `e.lanes[0]` only [measured: grep `lanes[1]`, 0 hits]. `Drain`
  syncs both streams [code: gpucheck.cpp:177-184]. D-G8 keeps "2 `Drain` syncs" as the baseline, and D5c/D7 add
  entries "on the scorer lane", which is also lane 0. The owner rule is "remove dead code".
- Fix: in D1 (with the counters), drop the second lane: `Lane lane;` and `Drain` syncs one stream. D-G8 then reads
  `drain_syncs` = 1 per call (measured in D1, as now), and the D9 ADR states "one lane". If a second lane is wanted
  later, add it together with its first user.

### 10. [minor] D (D5 step 1): wrong line ref for the second label pair
- Problem: "The final `label(1), label(0)` pair (:1615)". The pair is at gpucheck.cpp:1614. Line 1615 is the
  `zero(c.aux, ...)` before `Holes`. An implementer who gates line 1615 gates the wrong call.
- Fix: write ":1614". Add: "the `zero(c.aux)` at :1615 and `Holes` always run".

### 11. [minor] B (B3.1 tests, B.9.5 G07): the "GPU test job (Phase D GPU CI)" is still named (round-2 item 18, B side not done)
- Problem: B3.1 ends "Run `compute-sanitizer` in the GPU test job (Phase D GPU CI) once", and G07 says "run once in
  the GPU test job". D8 has no GPU CI. The sanitizer runs in R-D2's `make gpu-check`, and D8 already adds
  `tests/test_scorer_sums.py` to that corpus.
- Fix: B3.1 and G07: "runs in R-D2's `make gpu-check` sanitize step (D8 corpus, markers `gpu` and `sanitizer`)".

### 12. [minor] D (D1 step 0): the runbook fix leaves a wrong box torch version in K3's runbook
- Problem: D1 step 0 edits the three `GPU_RUN.md` files for paths and price. `K3/GPU_RUN.md:7` also says the box image
  has "torch 2.14.1 in `/opt/venv`". The box has torch 2.13.0+cu130 [measured: runs/r5/env_check.json; C header].
  K3's torch extension is built against that torch, so the runbook names the wrong ABI.
- Fix: add to D1 step 0: "set the K3 runbook's torch line to the `vllm_native` lock (torch 2.13.0+cu130), cited to
  `experiments/phase2/docker/locks/vllm_native.txt`". The preflight check that greps for `/[p]rivate/tmp` also
  checks for `torch 2.14.1` under `docs/implementation2.0/kernels/*/GPU_RUN.md`.

### 13. [minor] D (D5c pinned bound, D8 sanitize): the bound covers only the input staging, and the sanitizer run has no `-I`
- Problem: (a) D5c sizes each pinned block at its first call, but it gives a bound only for `score_polygons` input
  (20.5 MB per name). `Lane::out` also grows: defects (n x 13 x 8 B), `ObjOut`, signatures [code: gpucheck.cpp:1799,
  1572, 1623, 1638]. A larger later call grows it, which runs `cudaFreeHost` (a device sync) and adds a `pinned_grows` count
  that the gate then reads as a fail. (b) D8 runs `python -m pytest -m "gpu and sanitizer"` without `-I`. If any
  sanitizer case reads caches from `data/` (the corpus names "8 caches"), this breaks the repo's `-I` rule for
  untrusted data.
- Fix: (a) state the `out` bound: `n_max x (13 x 8 + sizeof(ObjOut) + 8) + n_max x 10 x 32 x 4 B` for signatures
  (`MAX_RINGS` components x a stated hole cap), with `n_max = 256 x names`. (b) The sanitize command uses
  `"$PY" -I -m pytest`, and the corpus states that its caches are synthetic fixtures under `tests/`, never `data/`.
