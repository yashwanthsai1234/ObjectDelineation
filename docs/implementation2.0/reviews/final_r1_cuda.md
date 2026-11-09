# Final review, round 1: CUDA and GPU performance (whole Implementation.md)

Reviewer role: CUDA / GPU performance (NVIDIA docs, Nsight, Blackwell sm_120, CCCL). Target:
`REPO/docs/implementation2.0/Implementation.md` (11,117 lines), plus the code and kernel files it names. No earlier
`final_r*_*.md` file existed, so nothing is skipped as "already fixed".

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | Catalog check reused from `(planning probe, not kept: nv_catalog_live.json)` (no new install). Only `compute-sanitizer-diagnose-cuda-bugs` fits; it shaped items 4 and 18 (leak-check and initcheck behaviour of the D8 chain). Nothing installed. |
| engineering:code-review | Each GPU claim checked against code: `cuda/gpucheck.cpp` (Block, Drain, wait sites, label chain, ScoreRows), `CMakeLists.txt`, `run_box.sh`, `fast.py` (VisionCache, ImageMemo), `serve.py` (fork_pool), transformers `_deepstack_process`. Found the C4 staging race (item 9), the Drain syncs not counted (item 19), the uninitialised flag (item 18). |
| engineering:debug | Reproduce-first lens on the gates: which gates can fail with no bug (item 8, G-EQ4 batch shape), which can never pass (items 2, 4, 10). |
| engineering:system-design | Process and memory topology on one 32 GB card: forked pool workers with CUDA contexts (item 3), view-bank placement vs fork-before-CUDA (item 5), vLLM's own CPU lane (item 12). |
| engineering:architecture | Ordering and budget as decisions: G1 contents vs the D sequence (item 1), Gate C value (item 2), spend before the data run (item 11). |
| engineering:documentation | Implementability: duplicate flags and undefined env names across A6/C0 (item 6), scratchpad paths that vanish (item 13), one planning rate (item 14). |
| huggingface-skills:huggingface-best | No model choice in this review; confirmed the 2B/4B section has no GPU-memory claim beyond the gated 15-step probe. No change. |
| huggingface-skills:huggingface-papers | GPU CCL/EDT references in D.2 are journal papers (no arXiv ids); no correction needed. No item. |

Web checks: compute-sanitizer leak-check counts leaks in `ERROR SUMMARY` and reports only at context destroy
[web: https://developer.nvidia.com/blog/debugging-cuda-more-efficiently-with-nvidia-compute-sanitizer ;
https://forums.developer.nvidia.com/t/how-to-find-leaks-cuda-gdb-runs-out-of-memory-but-compute-sanitizer-runs-without-erros/246420];
racecheck is shared-memory only [web: https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/ ;
https://developer.nvidia.com/compute-sanitizer].

Counts: **blocking 3, major 11, minor 11.**

---

## Blocking

### 1. [blocking] 00_front (section 6.2, back-matter decision 3) vs D: G1 runs `make gpu-check` (R-D2) before the D code can exist
- Problem: G1 item 5 runs D's `make gpu-check` (build, `pytest -m gpu`, 4 sanitizers, determinism, G7/G8). That chain
  is R-D2 and needs D2-D8 merged (D.8: "R-D2 (D2-D8)"; about 16 engineer-days). D says "Nothing below may claim a GPU
  speed until [D1/R-D1] gives numbers" (D1), and D2 step 2, D5 step 4, D6 step 1 take inputs from R-D1. The front
  flow puts D after Gate A, i.e. after the B runs, while G1 runs before G2 (B1). So G1 as listed is not executable
  before the data run, and the $1.62 / cap $1.74 figure is wrong. G1 also omits R-D1's own
  `bench_checker.py --engine current --stages` (D1 step 4, about 5 min).
- Fix: split G1. G1 (before G2) = setup, A3, C items 6-7 (only if C1-C5 passed on the Mac), K1/K2 runbooks plus
  `bench_checker.py --engine current --stages` (this IS R-D1). Move G1 item 5 out of G1: "R-D2 runs first in the first
  B/C rental after D2-D8 merge (D.8)". Recompute the G1 minutes and cap, and the G2 step count, from the new list.
  Change back-matter decision 3 to "R-D1 inside G1; R-D2 inside the first B/C rental after D8".

### 2. [blocking] 00_front (sections 1 and the flow diagram) vs C.6: "Gate C <= 18 s" decides whether B uses the C code
- Problem: the flow diagram says "Gate C: <= 18 s/step, same step-1 rewards ... no -> run on old code", and section 1
  says "Target 18 s per step or less (Gate C)". In C, the must-pass gate for B to use C code is G-EQ1..6 plus G-SP1
  (<= 20.0 s mode A, <= 20.7 s mode B); 18 s is G-C and needs C7 and C8, which are owner decisions that run INSIDE
  Phase B (C.6, C.8.1: C1-C5 estimate 18.7-21.2 s). As written, the front gate can never pass before B, so B always
  runs on old code. Section 6.3 then uses "20 s (Gate C passed)", a third reading.
- Fix: in section 1 and the diagram write "Gate C-adopt: G-EQ1..6 pass and G-SP1 (<= 20.0 s mode A / <= 20.7 s
  mode B); G-C (<= 18 s with C7+C8) is measured inside Phase B". Use "G-SP1" (not "Gate C") in section 6.3.

### 3. [blocking] D (D7, D7b, I06 test) vs E (E.3.3, E.5, E8 row): device eval work in forked pool workers; nobody deletes the pool
- Problem: exam metrics run in `serve.fork_pool` workers today (`experiments/phase2/final_eval.py:69-77` ->
  `select_ckpt.score_items` -> `pool.map(_score, ...)`, "one CUDA checker call per item"; `POOL_CAP = 16`
  [code: src/objdel/train/serve.py:132]). D7 moves `metrics()` pixel work onto the device (`eval_sets`, scratch cap
  256 MiB) and D7b builds tables per process at the first `device()`; D's test I06 even requires "after `os.fork()` a
  child can call `device()`". So each worker opens its own CUDA context and its own tables and scratch: up to
  16 x (context + 256 MiB + tables) beside vLLM and the trainer, and G14 measures only the audit callback. E states
  "Phase D moves those metrics to the GPU and deletes `pool.py`" (E.3.3, E.5, E8 file table), but no D task deletes it.
- Fix: add D7 step 6 "Process topology": `metrics()` device work runs in ONE process (the main process after vLLM
  generation, or the trainer process for the audit), batching `eval_sets` over items; `final_eval`/`select_ckpt` stop
  using `fork_pool` for metrics; delete `verify/pool.py` (or `serve.fork_pool` before E) in D9 with its tests. Change
  I06 to "a forked child must NOT call `device()`; it raises". Add a gate (G17): during `final_eval` on the box,
  `nvidia-smi --query-compute-apps` shows one process with the engine context, and peak memory <= 30 GB.

---

## Major

### 4. [major] D (D8, G5): `--leak-check full` against an engine that is never destroyed
- Problem: D8 runs memcheck with `--leak-check full`. The engine is "Never destroyed" by design and its arenas are
  never freed [code: cuda/gpucheck.cpp:169-174, Block 135-150]; D5c adds a pool with release threshold `UINT64_MAX`.
  Leaks count in `ERROR SUMMARY`, and leak reports appear only when the context is destroyed [web: NVIDIA blog,
  forum above]. So G5 ("0 errors") either fails by design or reports a meaningless 0. Today's memcheck step does not
  use `--leak-check` [code: experiments/phase2/configs/checker_gpu.yaml:16-21].
- Fix: add `gc.shutdown()` (sync both lanes, free every `Block`, every table, `cudaMemPoolDestroy`, destroy the two
  streams) and call it from a session-finish fixture in `tests/test_gpu_sanitize.py`; or drop `--leak-check full` and
  gate instead on `gc.stats()` live device bytes == 0 after `drop_tables` + `gc.trim(0)`. State which in D8.

### 5. [major] 00_front (2.1, 2.4, 4 "Data format") and back matter Q5/decision 11: the view bank has no task, and "load to pinned memory (or cuda)" breaks fork-before-CUDA and the memory budget
- Problem: (a) No phase task creates `data/viewbank.py` (it appears only in the layer list, section 4 and owner
  questions). (b) The STORE preload must happen before `fork_pool`, and `fork_pool` refuses once CUDA is initialised
  [code: serve.py:59-61, 167-169; g_probe.py:157]. Pinning 2.2+ GB at preload (`cudaHostAlloc`) initialises CUDA
  before the fork, so the run raises on the paid box. (c) "or cuda": 22,023 views x 512 x 512 B = 5.8 GB of uint8 grey
  planes uncompressed [estimate: arithmetic]; mode A already plans about 28 of 32 GB (C5 step 2), so a resident bank
  does not fit; the 2.2 GB is the compressed parquet, which needs CPU decompression per step.
- Fix: either drop the view bank from the front architecture until a task owns it, or add a B1/T1 task with: pageable
  in-memory planes preloaded before the fork; per step, copy the 32 needed planes into C1's persistent pinned buffer
  (25.2 MB); `bank_device: cuda` only if `torch.cuda.mem_get_info` after step 1 shows >= bank size + 2 GB free.
  Replace "load to pinned memory (or cuda)" in 2.4 with "pageable host memory (pre-fork); 32 planes per step to a
  pinned buffer".

### 6. [major] A (A6) vs C (C0): two Nsight installs, an undefined env name, two profiler-window flags
- Problem: A6 measured that NVIDIA's devtools apt repo "stops at 2021.4.1, too old for an RTX 5090" and fetches the
  pinned `nsight-systems-2026.5.2` .deb (SHA-256 given) at run time with `dpkg-deb -x`. C0 step 0 bakes
  `nsight-systems-cli` "from NVIDIA's devtools apt repo (or the standalone .run)" into the image. C0 step 2 uses
  `OBJDEL_NVTX_SYNC=0 (A6)`, but A6 defines only `OBJDEL_NVTX` and keeps the sync unchanged (A6 step 2). A6 adds
  `--profile-window S E` (NvtxSteps calls `torch.cuda.profiler.start/stop`); C0 adds `--profile-steps A-B` that does the
  same; C10 keeps both. Two callers of `cudaProfilerStart` in one run can open the capture twice.
- Fix: one nsys mechanism: A6's pinned .deb, unpacked with `dpkg-deb -x` in the last Dockerfile layer (C0 step 0
  cites A6, not the devtools repo). Define `OBJDEL_NVTX_SYNC` in A6 step 2 (default 1 = today's sync; 0 = CUDA events
  read in `log`). Delete `--profile-steps`; C0 uses `--profile-window 16 17`.

### 7. [major] 00_front (6.2) vs A3 / C.6 / C.8.3: A3 runs last in G1, but its saved answers feed the C gates and it is the item most at risk from the cap
- Problem: G-EQ2 and G-EQ3 need ">= 50 batches" built on the Mac from "Phase A's saved answers" (C.6, C.8.3), and A3
  (G1 item 8, last) is what produces them. So C's equality gates cannot be complete when G1 items 6-7 run, and the
  speed numbers "count" only after the rental. A3 is also the only item that unblocks Gate A and D7; with 161.5 min of
  items under a 174-min cap (7 % slack, 40-min setup estimated against measured 18-31 min plus toolchain), a slip cuts
  A3 off.
- Fix: run A3 first after setup (items 8 -> 2). State that the Mac half of G-EQ2/G-EQ3 runs after the rental on A3's
  answers, and that C code is adopted for B only after that Mac half passes.

### 8. [major] C (C3, G-EQ4, G-EQ6): bit-equality of the encoder-once path depends on the GEMM batch shape
- Problem: today `VisionCache` encodes the new images of each forward chunk (about 2 images, 2,048 patches)
  [code: fast.py:269-283]; C3 encodes all 32 images in one call (32,768 patches). cuBLAS picks kernels by shape
  (tile, split-K, stream-K), so the vision features can differ in the last bits with no bug, and G-EQ4 (`torch.equal` of
  logprobs and LoRA grads) fails; G-EQ6 asks `grad_norm` equal to 4 digits, and R5's base vs shared arms already
  differ at the 4th digit (C.R item 18). The debug order in C.7 does not separate this from a real ordering bug.
- Fix: add C3 step 0 (box, part of `c_speed_check --encoder`): encode the 32 images once batched and once in the
  chunked grouping; if bits differ, G-EQ4 compares the encoder-once path against the old path fed with the SAME
  32-image features (proves the plumbing exactly), and adds a numeric bound for the batched features (max |logprob
  diff| <= 2^-8 relative, record it). G-EQ6 then gates rewards, legal and mean length exactly (vLLM sampling does not
  depend on the trainer encode) and reports `grad_norm` with a stated tolerance.

### 9. [major] C (C1 step 1, C4 step 4, C.3.3): pinned staging reuse race and a wrong premise
- Problem: C4 step 4 builds indices "into one pinned int64 staging tensor and copy it with `non_blocking=True`" per
  chunk; if that tensor is reused while the previous async H2D is still queued (the host runs ahead of the GPU), the
  next chunk's numpy write corrupts the indices silently. C1 step 1 says "No per-step `pin_memory=True` allocation
  (each is a cudaHostAlloc)"; torch's pinned allocations go through the caching host allocator, which reuses blocks
  and records an event per non-blocking copy, so it is both cheap and race-free [code: torch CachingHostAllocator;
  estimate: not measured here].
- Fix: C4 step 4: "allocate each chunk's staging with `torch.empty(n, dtype=torch.int64, pin_memory=True)` (torch
  caching host allocator; safe reuse by recorded events), fill from numpy, `.to(device, non_blocking=True)`". If a
  persistent buffer is kept, require a CUDA event per use and `event.synchronize()` before the next host write, plus a
  test that runs two chunks back to back with a delayed stream. Remove the cudaHostAlloc sentence from C1 step 1 (keep
  the persistent buffer only if C0 shows the caching allocator costs time).

### 10. [major] C (G-SP4, C4 step 6, C4b): "< 200 trainer-side host syncs" is not reachable with the optional tasks off
- Problem: TRL's 7 metric `.item()` calls per micro-batch are 112 per step [doc C.3.3; code grpo_trainer.py:3381-3404];
  they go away only with C4 step 6 (adopted only if >= 0.1 s saved). Qwen3-VL's DeepStack merge uses boolean
  indexing (`hidden_states[visual_pos_masks, :]`, 3 DeepStack layers) in every forward chunk [code:
  .venv transformers/models/qwen3_vl/modeling_qwen3_vl.py:847-855]; C4b (optional, last) is the only task that removes
  it. With 16 reference + 16 update forwards, 112 + about 96 already exceeds 200 [estimate]. G-SP4 has no stated
  consequence when it fails.
- Fix: derive G-SP4 from the adopted set: "< 200 if C4 step 6 and C4b are adopted; else <= 25 % of C0's count
  excluding the TRL metric reads and DeepStack indexing, which are reported apart". State that G-SP4 is report-only
  for B adoption (only G-EQ1..6 and G-SP1 decide, C.8.3).

### 11. [major] 00_front (6.2 item 4) vs C (C1 step 7, C.8.2), and the order of spend vs the owner's data-first rule
- Problem: G1 runs the K3 runbook unconditionally (15 min, $0.15), but C says K3 runs only if both "hf" and "torch"
  backends fail the box self-check (C1 step 7; C.8.2 lists it as optional). K1/K2 first GPU runs (R-D1) and the D items
  serve Phase D, which comes after the B runs. Together they take about $0.6-0.7 of the $3 before the data-fix run and
  cut B1 from 300 planned steps to 160-190 (section 6.3). The owner rule says dataset quality is the top priority.
- Fix: make G1 item 4 conditional ("only if `c_speed_check --pixels` fails both backends"). Add an owner decision row:
  "R-D1 in G1 (before B1) or inside the first B rental (after B1)"; default inside the B rental, so the B1 run gets
  about 60 more minutes [estimate: 0.13 + 0.15 + 0.30 dollars at $0.60/h].

### 12. [major] 00_front (2.1, 2.4, 2.5): the end-state diagrams hide vLLM's CPU image lane and show unadopted kernels
- Problem: 2.1 says "Pixels never come back to the CPU after the view bank is loaded" and 2.4 shows only "25 MB H2D"
  and "vision encoder once". C.3.1/C.3.3 show that vLLM still gets 256 PIL references, runs its own CPU HF
  preprocessing, does about 201 MB fp32 H2D per step and a duplicate vision encode; only the C6b investigation
  ([unverified] vLLM DeepStack embeds) could remove it. 2.4 labels the patchify "K3a" and 2.5 shows "fused log-prob
  K3b" and "32 x n=8", but C1 uses torch/HF ops first, K3(b) is "Not in Phase C" (C.5), and C6 is default-off.
- Fix: 2.4: add the vLLM lane ("CPU: hash + HF preprocess of 32 images; ~201 MB fp32 H2D; duplicate encode; removed
  only if C6b passes") and list it as a CPU item with its reason in 2.1. Replace "K3a" with "C1 patchify (HF/torch;
  K3a fallback)", drop "K3b" from 2.5, and mark "32 x n=8" as "C6 arm, default off".

### 13. [major] 00_front (5), B.5, C, D, E.5, CLEAN.5, DOC.5, kernel runbooks: kernel paths point at a session scratchpad
- Problem: section 5 says the drafts live in `docs/implementation2.0/kernels/`, but the body cites `docs/implementation2.0/kernels/...`
  30 times, and B.5 gives absolute `/[p]rivate/tmp/claude-501/...` paths. The copied runbooks still set
  `S[P]=/[p]rivate/tmp/...` and build the K1 bundle from `$REPO/data/k1_bundle/lidc` [code:
  docs/implementation2.0/kernels/K1/GPU_RUN.md:18-30]; the kernels README admits it. `/[p]rivate/tmp` is session-only, so
  section 6.1 step 3 and G1 items 2-4 cannot be run from the repo.
- Fix: replace every `docs/implementation2.0/kernels/Kn` with `docs/implementation2.0/kernels/Kn`. In each GPU_RUN.md set
  `K1=$REPO/docs/implementation2.0/kernels/K1` and write work files to `$REPO/data/k1_bundle/` (git-ignored; LIDC only),
  regenerated from `data/checker_eval/masks.npz` by `cross_numpy.py`. Add a repo test that no file under
  `docs/implementation2.0/` contains `/[p]rivate/tmp/`.

### 14. [major] 00_front (6), B.8.1, C.8.2, D.8: four hourly rates, one of them wrong, and four step times
- Problem: section 6 says "Rate: $0.440/h (R3) to $0.641/h"; $0.440 is R3's billed TOTAL, not a rate
  [measured: experiments/phase2/vast/ledger.csv line 100 "billed=$0.440"; launch dph 0.5437]. D.1 step 5 already
  corrects this to $0.544-0.628/h [measured: ledger rows r3, r4, r5]. B.8.1 plans at $0.641/h, section 6 at $0.60,
  C.8.2 at $0.593. Step times used for money: 17 s (B.8.1 "with Phase C"), 18 s (front Gate C), 20 s (front 6.3),
  23.65 s, 27.5 s (R3 risk row "base case"). The owner's budget decision (back matter 1) depends on which pair is used.
- Fix: one planning table in section 6: rate $0.60/h planning, $0.544-0.628/h measured range; step time 23.65 s base
  (R5), 27.5 s pessimistic (shifted rows, R3), 20.0 s only after G-SP1 passes; delete 17 s from B.8.1. Every phase
  cost table cites that table.

---

## Minor

### 15. [minor] 00_front (6.2 item 1, 6.4 rule 3) vs D (R13): stop rule 50 min vs 35 min
- Problem: front tears down if no first sync by 50 min; D R13 says 35 min (r5's good launch took 31 min).
- Fix: one rule in 6.4 (35 min, from the ledger), cited by D and G1 item 1.

### 16. [minor] 00_front (6.2): fixed `--cap-hours 2.9` assumes $0.60/h
- Problem: at the measured $0.628/h plus internet cost the cap is about $1.82, not $1.74, and G2's arithmetic uses
  $1.74.
- Fix: compute `--cap-hours = floor_to_0.05(1.74 / dph_total)` from the picked offer and write the result in the
  launch line.

### 17. [minor] 00_front (6.2 toolchain line) vs D (common toolchain, D1): missing `nvtx` extra
- Problem: front lists `cuda-toolkit[nvcc,cccl]==13.0.3`; D1 makes CMake stop with `FATAL_ERROR` when the
  `nvidia-nvtx` wheel is missing, and installs `[nvcc,cccl,nvtx]`.
- Fix: front: `cuda-toolkit[nvcc,cccl,nvtx]==13.0.3` "as `run_box.sh` does after D1".

### 18. [minor] D (D5 step 1): the `changed[o]` flag must be zeroed per call
- Problem: "The flag lives in the arena"; arenas are reused and never cleared [code: gpucheck.cpp:135-150]. A stale 1
  only costs a relabel, but reading never-written arena memory is an initcheck error, so G5 fails.
- Fix: add "`zero(changed, sizeof(i32) * n, st)` before `CropFill`" and a CPU-build test that two calls with different
  objects give the same skip pattern as fresh calls.

### 19. [minor] D (G8, engine counters): `Drain` syncs are not counted
- Problem: G8 expects 3 waits per `score_polygons` [code: gpucheck.cpp:1574, 1626, 1639], but every entry also runs
  `Drain`, which calls `cudaStreamSynchronize` on both lanes directly, not through `wait()` [code: gpucheck.cpp:177-184].
  Real host syncs are 5 per call; the counter would hide a new one.
- Fix: route `Drain` through a counted `wait()` (or count it as `drain_syncs`) and set G8 to the measured total.

### 20. [minor] C (C0 step 1): nsys validation of the sync counter includes non-torch syncs
- Problem: `set_sync_debug_mode` sees only torch syncs, but nsys `cuda_api_sum` in the update range also counts the
  engine's `cudaStreamSynchronize` from the `prefetch_rewards` thread (fast.py:629) and its `Drain`.
- Fix: compare per thread (nsys `cuda_api_trace` has the TID) and subtract `gc.stats()` waits over the same steps.

### 21. [minor] C (C0 step 6a): CPU sampling may be unavailable in a Vast container
- Problem: step 6a checks vLLM input time "against the osrt CPU samples"; nsys CPU sampling needs
  `perf_event_open`, often blocked in containers.
- Fix: say "if `nsys status -e` reports CPU sampling unavailable, use the NVTX `generate` start to the first vision
  kernel only".

### 22. [minor] C (C.8.2 budget note): "$0.60 left under the $17.50 limit"
- Problem: contradicts the binding $3 GPU budget in section 6.
- Fix: delete the sentence; cite section 6.

### 23. [minor] E (E.3.3, E.5) vs D (D3, D.7 R10) and front order: the polygonizer seam order
- Problem: E says E3a exposes `geom.polygonize` first and D "edits one function and no caller"; D3 creates
  `checker/rings.py`, edits callers (`targets/build.py`, `data/riga.py`, `certify.py`) and lands before E (R10, front
  flow).
- Fix: E.3.3/E.5: "D3 lands first on today's paths; E3a then moves `checker.rings` to `geom` unchanged".

### 24. [minor] DOC / CLEAN: kernel `build/` folders hold macOS binaries
- Problem: `docs/implementation2.0/kernels/K1/build` (3.5 MB, Mach-O and a `.dSYM`) and `K2/build` (10 MB) would enter
  git when DOC1 puts docs under version control [measured: du -sh].
- Fix: DOC1 adds `docs/implementation2.0/kernels/*/build/` to `.gitignore` (keep `results/` text logs), or CLEAN deletes
  the binaries before DOC1.

### 25. [minor] C (C3, C.3.2 steps 6 and 8): no guard for 32-image `pixel_values` with a 256-row grid
- Problem [estimate: a risk, not an observed defect]: with `per_image` True, pass 2 returns 32 images of patches with a 256-row CPU grid (C3.2 step 6). Any TRL
  code that splits `pixel_values` by `image_grid_thw` before our pop (e.g. `split_pixel_values_by_grid` is skipped only
  because the key is popped at the END of `_generate_and_score_completions`, C3.2 step 8) or that runs the reference
  pass outside `FastForward` would read misaligned rows with no error; the only guard is "raise if
  `row_image_ids` is None".
- Fix: add a start-up assertion and a Mac test: under `--encoder-once`, any call that reaches the vision tower with
  `pixel_values.shape[0] != image_grid_thw.prod(-1).sum()` raises (`EncoderMismatch`), and the C.9 integration test
  runs TRL's full `_generate_and_score_completions` on the tiny model with that check on.
