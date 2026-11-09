# Final review, round 2: CUDA and GPU performance (whole Implementation.md)

Reviewer role: CUDA and GPU performance (NVIDIA docs, Nsight, Blackwell sm_120, CCCL). Target:
`REPO/docs/implementation2.0/Implementation.md` (13,065 lines), plus the code, kernel and runbook files it names.
Earlier file read: `reviews/final_r1_cuda.md`. Items that round 2 fixed are not repeated (r1 items 1-6, 8-12 and
15-25 are fixed or answered). Where an r1 fix created a new defect, the item below says so.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No new catalog skill fits beyond `compute-sanitizer-diagnose-cuda-bugs` (already named by the doc). It pointed me to NVIDIA's own package indexes for the Nsight Systems and NVTX checks (items 1, 12). Nothing installed. |
| engineering:code-review | Each GPU claim checked against code: `CMakeLists.txt`, `run_box.sh`, `gpucheck.cpp` (alloc, Drain, wait sites, add/drop_tables, ScoreRows), `fast.py` (`set_is_ratio`, `log`), `reward/cache.py` (`device()`), `pyproject.toml`. Found items 1, 7, 9, 15, 16. |
| engineering:debug | Reproduce-first lens on gates: which gates fail with no bug (items 4, 5), which cannot be measured as written (item 8), which commands fail at once (items 1, 9). |
| engineering:system-design | Process, memory and stream topology on one 32 GB card: the stream-ordered pool vs the device default pool (item 7), host snapshot size in mode A (item 14), pinned growth syncs (item 16), table builds per step (item 19). |
| engineering:architecture | Ordering and budget: V1/V2 contents against Phase A's text (item 2), V2 minutes against the ledger (item 3), runbook fix order vs R-D1 (item 10). |
| engineering:documentation | Implementability: one nsys command (item 11), one planning rate (item 17), stale "GPU CI" references (item 18). |
| huggingface-skills:huggingface-best | No model choice in this review. Used the Qwen3-VL-2B `preprocessor_config.json` from the local HF cache (mean = std = 0.5 on all channels), which makes item 5 likely. |
| huggingface-skills:huggingface-papers | No CUDA claim in the doc rests on an arXiv paper that needs a new check; D.2's CCL/EDT sources are journal papers. No item. |

Web and package checks made for this review (2026-11-08):
- PyPI JSON for `nvidia-nvtx 13.0.85` and the x86_64 wheel's file list: the headers are at
  `nvidia/cu13/include/nvtx3/` (for example `nvidia/cu13/include/nvtx3/nvtx3.hpp`). There is no `nvidia/nvtx` package
  directory [measured: `unzip -l nvidia_nvtx-13.0.85-...manylinux1_x86_64.whl`].
- PyTorch `aten/src/ATen/native/cuda/BinaryDivTrueKernel.cu` (main): "if the second operand is a CPU scalar, compute
  a * reciprocal(b). Note that this may lose one bit of precision" [web: raw.githubusercontent.com/pytorch/pytorch].
  On the 256 uint8 values, `(x - 127.5) / 127.5` and `(x - 127.5) * (1/127.5)` differ in 126 of 256 float32 results
  [measured: numpy, this review].
- NVIDIA devtools repo `devtools/repos/ubuntu2404/amd64/Packages` lists `nsight-systems-cli-2026.6.1`
  (233 MB, Depends: libc6, libglib2.0-0) [measured: curl]. The A6 pin URL (CUDA debian13 repo,
  `nsight-systems-2026.5.2...deb`) returns HTTP 200, 458 MB [measured: curl -sI].
- Nsight Systems user guide: `--cuda-memory-usage` "may cause significant runtime overhead"
  [web: archive.docs.nvidia.com/nsight-systems/2024.6/UserGuide].
- nvidia-smi inside a container reports host PIDs, which do not exist in the container's PID namespace
  [web: github.com/NVIDIA/nvidia-docker/issues/1460; forums.developer.nvidia.com/t/343361].
- compute-sanitizer 13.4.92 already ran on the R580 driver (580.159.03) with 0 errors [measured:
  runs/rl/cuda_checker/memcheck.txt, runs/k1/cuda_checker/memcheck2.txt; runs/r5/env_check.json]. No item.

Counts: **blocking 1, major 9, minor 10.**

---

## Blocking

### 1. [blocking] D (D1 files, D.9 S06, D.R item 12): the NVTX include path names a package that the CUDA 13 wheel does not have
- Problem: D1 says CMake adds "the include directory of the installed `nvidia-nvtx` wheel, found by asking the build
  Python for the `nvidia.nvtx` package path, with a clear `FATAL_ERROR` when it is missing", and test S06 checks
  that "CMake fails clearly when `nvidia.nvtx` is missing". The `nvidia-nvtx 13.0.85` wheel (the one
  `cuda-toolkit[nvtx]==13.0.3` pulls) has no `nvidia/nvtx` directory. Its headers are in
  `nvidia/cu13/include/nvtx3/` [measured: wheel file list]. So the check fires on every box and on
  `docker_check.sh`: `make cuda` stops, and with it every GPU config that calls `run_box.sh` (the reward engine
  build included). `run_box.sh` already sets `CUDA_HOME=nvidia/cu13` [code: run_box.sh:25], so the headers are on
  nvcc's default include path once the extra is installed.
- Fix: in `CMakeLists.txt` under `if(OBJDEL_CUDA)`: after `find_package(CUDAToolkit REQUIRED)`, add
  `if(NOT TARGET CUDA::nvtx3) message(FATAL_ERROR "nvtx3 headers not found under ${CUDAToolkit_INCLUDE_DIRS}: install cuda-toolkit[nvtx]") endif()`
  and `target_link_libraries(gpucheck PRIVATE CUDA::nvtx3)` (FindCUDAToolkit has had `CUDA::nvtx3` since CMake 3.25;
  the file needs 3.26). Remove "ask the build Python for `nvidia.nvtx`". Change S06 to: "`docker_check.sh` passes;
  with the `nvtx` extra uninstalled, CMake stops with the message above". Same words in D.R item 12.

---

## Major

### 2. [major] 00_front (6.2 item 4, open question 8) vs A (A3, A6 step 8, A-G6, A.8, G06): the V1 nsys window does not exist in Phase A, and Phase A still runs in "G1a"
- Problem: the front says V1 item 4 is "A3 ... plus the A6 nsys window (6 steps, capture steps 4-5)" (28 min, 10 of
  them for nsys), and the front review log cites "A.md nsys command `--profile-window 4 5`". No such command exists
  [measured: grep]. Phase A says the opposite: "The plan has no nsys step. The A6 Nsight baseline comes from the C0
  baseline arm of the same rental (A6 step 8)"; A-G6 and G06 read "the C0 baseline arm capture of G1a"
  (steps 16-17). A also says "A3 is the first item of G1a after setup" (A.3.2, A3 step 3, A.8), while the front puts
  A3 last in V1 with a $1.90 skip rule, and the C0 baseline arm in V2. So V1 bills 10 minutes for a capture no config
  defines, and A-G6 cannot pass until V2.
- Fix: pick one. Recommended: drop the nsys part from V1 item 4 (A3 = val + test eval only, 17.5 min sampling +
  5-8 min load = about 23 min [estimate]); A-G6's box part stays on the C0 baseline arm in V2. Then in A replace
  every "G1a" with "V1" for A3 and "V2" for the C0 capture, replace "A3 is the first item after setup" with "A3 runs
  after the B1 run in V1 (front 6.2 item 4, skip rule $1.90)", and in the back matter change open question 8 to
  "C0 baseline arm in V2". Recompute V1's total (about 197 min = $2.10 at 27.5 s).

### 3. [major] 00_front (6.3 V2 table, decision 1) and C (C.8.2): the arm minutes omit process start-up, so "items 1-3 fit if V1 runs near 23.65 s" is false
- Problem: C.8.2 takes the arm minutes from step time only: "R5 base arm 7 min 45 s for 15 steps" = 465.8 s =
  15 x 31.05 s, and "R5 shared arm 5 min 54 s" = 354.8 s, both `timing_s.train` [measured: runs/r5/base/rl.json,
  runs/r5/shared_rl.json]. Each arm is a new `train_rl` process (model load, vLLM init, compile and graph capture),
  and `c_speed_check` also loads the model. The ledger shows this cost: r4 `run` 14:07:06 to `teardown` 14:38:24 =
  31.3 min for 15.0 min of steps; r5 `run` 16:13:57 to `teardown` 17:11:15 = 57.3 min for 13.7 min of steps +
  17.5 min of eval [measured: ledger.csv rows 125-126, 133-134]. That is about 5-8 min per process [estimate].
  V2 items 2-3 start 4 processes, so about 20-30 min are missing: items 1-3 are about 95-100 min = $1.01-1.07 at
  $0.641/h, above the $0.79 left. Stop rule 6.4 #2 (90 % of the approved amount) would also stop new items at
  $0.71, before the check arm.
- Fix: add a "start-up per process" row (6 min [estimate from ledger r4/r5]) to C.8.2 and front 6.3, and recompute.
  Cut one start-up: the check arm continues in the fast arm's process after step 17 (steps 18-32 with
  `--check-rope-every 1 --check-meta-every 1 --record-batches` switched on at step 18 by a step-indexed flag; times
  after step 17 are not used). State in 6.3 and decision 1 that V2 items 1-3 do not fit $3 at either step time, so
  C adoption needs the owner's budget decision.

### 4. [major] C (C-G-EQ4, C3 step 0, C.7 row "Batched encode"): the 2e-3 max-logprob bound is not calibrated and can fail with no bug
- Problem: when C3 step 0 shows different bits, G-EQ4 adds "max abs per-token logprob difference between the batched
  and the chunked encode <= 2e-3" as part of Gate C-adopt. A last-bit change in bf16 vision features goes through 28
  decoder layers. A max over about 256 x 131 tokens is a tail statistic, and the doc gives no measurement of its size
  for a same-code perturbation (R5's 0.018 is a median between two different engines). Today's `VisionCache` already
  encodes a varying subset per call [code: fast.py:269-291], so today's own bits depend on the chunk grouping. The
  noise floor can be measured, but the gate does not do it.
- Fix: in C3 step 0 also encode the 32 images in two different chunk groupings of today's path (A and B) and run the
  reference pass on both. Record `floor_max` and `floor_p99` of |logprob(A) - logprob(B)|. G-EQ4's second check
  becomes: batched vs chunked-A `p99 <= 2 x floor_p99` and `max <= 2 x floor_max`, both recorded. If the floor is 0
  (bits equal), keep `torch.equal`. Drop the fixed 2e-3.

### 5. [major] C (C1 interface "torch" backend, C1 tests): a Python scalar divisor makes CUDA multiply by a reciprocal, so the Mac test passes and G-EQ1 fails on the box
- Problem: the "torch" backend is written as `x.float().sub_(mean').div_(std')`. Qwen3-VL has mean = std = 0.5 on all
  three channels [measured: local HF cache `preprocessor_config.json`], so the fused values are one number (127.5),
  and an implementer may pass Python floats. PyTorch's CUDA true-divide then computes `a * (1/b)` [web:
  BinaryDivTrueKernel.cu], while the CPU float32 path divides. `(x - 127.5) / 127.5` and `(x - 127.5) * (1/127.5)`
  differ in 126 of the 256 uint8 values [measured]. The Mac tests run with `device="cpu"`, so they pass. The failure
  appears only in G-EQ1 on the paid box.
- Fix: in the C1 interface write: "`mean'` and `std'` are float32 tensors of shape `(3, 1, 1)` (or broadcast to the
  patch layout) on `device`, never Python numbers or 0-dim CPU tensors". Add a Mac unit test that asserts
  `isinstance(std_, torch.Tensor) and std_.device == x.device and std_.dim() > 0` inside the torch backend, and a
  CPU test that shows the risk: `x / 127.5 != x * (1 / 127.5)` for some uint8 x (fails if someone "simplifies" it).

### 6. [major] C (C2 interfaces and step 2): "bad" rows are read once in `log`, yet they must take the HF fallback in the same step
- Problem: `rope_positions` returns `bad (B,) bool` with "No host read inside", and step 2 says "accumulate
  `bad.any()` into a device flag read once in `log`. Rows with `bad` take HF `get_rope_index`". The host cannot route
  a row to the fallback without reading `bad`, so either a sync is hidden in every chunk, or bad rows train with
  wrong positions until `log` sees the flag.
- Fix: compute `bad` on the host, with no device read, from data the host already has: `prompt_ids` per row (C2 memo)
  and the completion id lists from vLLM (C4 `RowMeta`). A row is bad when its prompt ids differ from the template ids,
  or when its completion ids contain an image/video/vision token id. Pass the host list of bad rows to
  `rope_positions`, which uses HF `get_rope_index` for those rows only. Keep the device check (step 3,
  `--check-rope-every`) as the guard; on a mismatch it raises.

### 7. [major] D (D5 step 5, D5c): `cudaMallocAsync` does not take a pool; the "engine-owned `cudaMemPool_t`" needs `cudaMallocFromPoolAsync`
- Problem: D5c says device `alloc`/`release` "use `cudaMallocAsync` / `cudaFreeAsync` on the lane stream from one
  engine-owned `cudaMemPool_t` with `cudaMemPoolAttrReleaseThreshold = UINT64_MAX`". `cudaMallocAsync` allocates from
  the device's current pool, not from a pool the engine created. An implementer then either (a) sets the threshold on
  the engine pool but allocates from the default pool, so memory is released at every sync and the counters still
  pass, or (b) calls `cudaDeviceSetMemPool`, which changes process-wide state for torch, vLLM and any library that
  uses the stream-ordered allocator.
- Fix: write "create the pool with `cudaMemPoolCreate` (location = the current device, `cudaMemAllocationTypePinned`);
  allocate with `cudaMallocFromPoolAsync(&p, bytes, pool, stream)`; free with `cudaFreeAsync`; never call
  `cudaDeviceSetMemPool` or change the default pool's attributes". Add a box test in `test_gpu_sanitize.py`: after
  warm-up, `cudaMemPoolGetAttribute(pool, cudaMemPoolAttrReservedMemCurrent)` is greater than 0 and stays the same over
  50 calls, and the default pool's reserved memory does not grow (exposed by `gc.stats()`).

### 8. [major] D (D-G17, D7 step 6): `nvidia-smi --query-compute-apps` cannot name container PIDs
- Problem: D-G17 passes when "every listed pid is the main process". In a Vast container, the driver reports host
  PIDs, which do not exist in the container's PID namespace [web: nvidia-docker issue 1460; NVIDIA forum 343361]. So
  no listed pid matches the main process: the gate fails with no bug, or an engineer weakens it to nothing.
- Fix: gate on what the container can see. (a) The number of rows of
  `nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader` is constant during `final_eval` and equals
  the count before the scoring starts (no new context appears). (b) Each pool worker (if any is left before CLEAN16)
  writes `torch.cuda.is_initialized()` and `gc.stats()["contexts"]` (new counter: 1 after the first engine call) to
  its log, and both must be 0. (c) Peak `memory.used` from `nvidia-smi --query-gpu=memory.used` <= 30 GB.

### 9. [major] D (D-G7): the command that builds the pre-D2 engine fails as written
- Problem: D-G7 builds the old engine from
  `git archive <pre-D2 sha> cuda cpp src CMakeLists.txt pyproject.toml` and `make cuda` in `/opt/old`. There is no
  `cpp` path in the repo: `git archive HEAD cuda cpp ...` exits with "fatal: pathspec 'cpp' did not match any files"
  [measured: read-only git archive]. The list also omits `Makefile` (needed by `make cuda`) and `README.md`
  (`readme = "README.md"` in `pyproject.toml:5`, which the build reads). This runs inside R-D2, on the paid box.
- Fix: `git archive -o old_engine.tar <pre-D2 sha> cuda src CMakeLists.txt pyproject.toml Makefile README.md`. On the
  box: `uv venv /opt/old && uv pip install --python /opt/old/bin/python numpy && make -C <dir> cuda PY=/opt/old/bin/python`.
  Add a Mac dry check to `make gpu-check DRY=1`: unpack the archive in a temp folder and run the CPU build
  (`uv pip install --python <tmp venv> -e <dir>`), so a missing file fails at $0.

### 10. [major] 00_front (section 5, 6.1) vs D (D1 steps 4-5): the runbook path fix comes after the run that needs it
- Problem: front section 5 says "Before any rental, each `GPU_RUN.md` must use repo paths", and the three runbooks
  still set `S[P]=/[p]rivate/tmp/...` or rsync from the scratchpad [measured: grep, K1 GPU_RUN.md:18, K2:33, K3:15]. The
  only task that edits them is D1 step 5, which runs after step 4 (R-D1). The runbooks also move files with manual
  `scp` and install the toolkit by hand, outside `sync_up.sh`/`run.sh` (no heartbeat, no ledger row per step).
- Fix: move the path edit to a new D1 step 0 (Mac, $0, before V2): `K1=$REPO/docs/implementation2.0/kernels/K1` (and
  K2, K3), work files in `$REPO/data/k1_bundle/`, price line citing front section 6. Add an R-D1 config
  `experiments/phase2/configs/r_d1.yaml` (like `checker_gpu.yaml`) that runs the K1 and K2 runbook steps through
  `run.sh`, and a preflight check that no file under `docs/implementation2.0/kernels/` contains `/[p]rivate/tmp`.

---

## Minor

### 11. [minor] C (C0 step 3) vs A (A6 step 8): two nsys command lines; C0's adds `--cuda-memory-usage=true`
- Problem: A6 says `nsys_capture.sh` is "the one nsys mechanism" and C0 calls it, but C0 step 3 writes its own
  `nsys profile ... --cuda-memory-usage=true`. NVIDIA documents that flag as "may cause significant runtime overhead"
  [web: Nsight Systems user guide]. It slows steps 16-17, which feed G-SP6 (GPU-active share) and A-G6's overhead
  check.
- Fix: C0 step 3: "`nsys_capture.sh {rd} <arm> 16 17 -- <train_rl command>` (A6 step 8); no other nsys flags". Peak
  memory comes from `torch.cuda.max_memory_allocated` and `nvidia-smi` (G-MEM), not from nsys.

### 12. [minor] A (A.0 skills row, A.R item 10) and C (C0 step 0): "the devtools apt repo stops at 2021.4.1" is false
- Problem: A6 checked the legacy path `devtools/repo-deb/amd64/Packages`. The current repo
  `devtools/repos/ubuntu2404/amd64/` lists `nsight-systems-cli-2026.6.1` (233 MB, Depends: libc6, libglib2.0-0)
  [measured: curl]. The decision (a pinned .deb, SHA-256, `dpkg-deb -x`) still holds, but the reason is wrong and
  the pinned full package is twice the size (458 MB, 27 GUI dependencies).
- Fix: replace the sentence with the measured fact. Optional: pin the CLI-only
  `NsightSystems-linux-cli-public-2026.6.1.218-39266636.deb` from the ubuntu2404 devtools repo with its SHA-256 (from
  the same `Packages` file), keeping `dpkg-deb -x` and `nsys status -e`.

### 13. [minor] D (D1 files and step 1): `GC_RANGE(name)` declares `r`, so two ranges in one scope do not compile
- Problem: `GC_RANGE(name)` = `nvtx3::scoped_range r{name}`. D1 puts one at the top of each entry and one per stage
  inside the same function; a second use in the same block redeclares `r`.
- Fix: `#define GC_CAT2(a, b) a##b`, `#define GC_CAT(a, b) GC_CAT2(a, b)`,
  `#define GC_RANGE(name) nvtx3::scoped_range GC_CAT(gc_range_, __LINE__){name}` on CUDA; empty on the CPU build. Put
  each stage range in its own `{ }` block so the range ends with the stage.

### 14. [minor] C (C5 step 5(b), G-EQ5): a full copy of every vLLM parameter may not fit on the device in mode A
- Problem: the equality check "copy every vLLM parameter (int16 views)" is about 2.13 B x 2 B = 4.3 GB [estimate].
  Mode A plans about 28 of 32 GB in use (C5 step 2), so a device copy can OOM on the paid box.
- Fix: copy to host memory one tensor at a time (`p.view(torch.int16).cpu()`), or compare per tensor on the device
  without a snapshot: store a per-tensor SHA-256 of the host copy after TRL's push, then hash again after the one-call
  push. State the choice in C5 step 5.

### 15. [minor] B (B2.1): quantiles "over the step's sequences" cannot use "the same min/max merge rules"
- Problem: `set_is_ratio` runs once per update micro-batch (16 per step) and merges `is_stats` by sums, min and max
  [code: fast.py:581-602]. B2.1 adds q50, q99, q99.5 of the per-sequence mean gap "over the step's sequences" and
  says to concatenate them to the same vector with the same merge rules. A quantile cannot be merged that way.
- Fix: keep the per-sequence mean gaps (and advantages for the OPSM counts) as a device list per step
  (`self.seq_gaps.append(t)`; 256 floats in all), and compute the quantiles once in `log` with one
  `torch.cat` + `torch.quantile` before the existing single `.tolist()`. Clear the list in `log`. Test: 16
  micro-batches give the same quantiles as one batch of all rows.

### 16. [minor] D (D5c, D-G15): pinned staging growth after warm-up is a device-wide sync, and the gate window cannot see it
- Problem: D5c keeps `cudaMallocHost` for `Block<false>` "but grows only at warm-up". `Block::get` grows whenever a
  call needs more bytes [code: gpucheck.cpp:135-150]; a larger batch at step 150 (long answers, a bigger `eval_sets`
  call) frees and reallocates, and `cudaFreeHost` synchronizes the device. D-G15 counts only steps 2-6.
- Fix: size the pinned staging at the first call from a stated bound (for example 256 answers x the `MAX_VERTICES`
  polygon bytes for `score_polygons`, the scratch cap for `eval_sets`), and add a counter `pinned_grows`. D-G15
  checks `pinned_grows == 0` over the whole B run (read at the end), not only steps 2-6.

### 17. [minor] C (C.8.2, C.R3 item 15): "the front's section 6 planning rate, $0.60/h" does not match the front
- Problem: front section 6 now sets the planning rate at $0.641/h. C.8.2 still calls $0.60 the front rate and uses it
  for its main column and the payback line (r1 item 14 residue).
- Fix: C.8.2: "Rate: front section 6 planning rate, $0.641/h; measured range $0.544-0.628/h". Recompute the C rows
  and the payback line at $0.641/h; drop the $0.60 column.

### 18. [minor] B (B3.1 tests, B.9.5 G07): "the GPU test job (Phase D GPU CI)" does not exist
- Problem: D8 says "NO `.github/workflows/gpu.yml` in Phase D"; the GPU chain is `make gpu-check` in R-D2. B3.1 and
  G07 still send the sanitizer run to a "GPU test job / Phase D GPU CI".
- Fix: "run in R-D2's `make gpu-check` sanitize step (D8), which includes `tests/test_scorer_sums.py`". Add the test
  to the `gpu and sanitizer` corpus in D8.

### 19. [minor] D (D7b, D-G8): table builds happen per cache in the reward thread, and the replay gate excludes them
- Problem: after D7b, the first `device()` call of a cache builds its tables on the device. In the B1 run each of the
  9,600 rows is used once (300 steps x 32), so every step builds about 32 new table sets. `device()` is a per-cache
  lazy call, so this is 32 separate `add_tables_from_raters` calls, each with a lock, an H2D copy and a wait, not
  the batched call D7b designs. D-G8's replay scores on caches that are already uploaded, so this cost is not gated.
- Fix: in `TRLReward.prefetch` (before scoring), collect the step's caches whose `_device` is unset and build them in
  one `add_tables_from_raters` call. Add to D-G8: "replay with cold caches (tables dropped before each repeat):
  median time <= 1.05 x the old `add_tables` path, and 1 wait per step for table builds".

### 20. [minor] C (C1 step 7) and E (E.5): "a separate torch extension `objdel._k3` built only with the `train` extra" cannot be built that way
- Problem: a PEP 517 extra does not change what the build backend compiles. `pyproject.toml` builds with
  scikit-build-core and `requires = ["scikit-build-core==1.1.1", "nanobind==3.1.0"]` [code: pyproject.toml:35-37]. A
  torch extension needs torch at build time, and build isolation would pull PyPI torch (Mac pin 2.14.1), while the box
  runs torch 2.13.0+cu130: an ABI mismatch.
- Fix: if K3a is adopted, build it on the box against the image's torch: `cuda/k3/setup.py`
  (`torch.utils.cpp_extension.CUDAExtension`, `extra_compile_args={"nvcc": ["--fmad=false", "-arch=sm_120"]}`) and
  `pip install --no-build-isolation cuda/k3` in the image build, or `torch.utils.cpp_extension.load` at start-up with
  the same flags. Remove "with the `train` extra". `gpu_pixels.py` imports it lazily and the self-check decides.
