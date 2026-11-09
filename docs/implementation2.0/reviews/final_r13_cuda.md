# Final review round 13: CUDA and GPU performance (Implementation.md)

Reviewer role: CUDA / GPU performance (kernels, GPU placement, syncs, memory, Blackwell sm_120, rental plan).
Scope read: 00_front sections 5 and 6 (6.1-6.4), C.md C1-C6 and C.6, D.md in full (D.0-D.9), the E/CLEAN/DOC
lines that name the engine, and the round 12 review, verdict and review-log rows. Code opened: `cuda/gpucheck.cpp`,
`docs/implementation2.0/kernels/K1/k1poly.cpp`, `CMakeLists.txt`, `Makefile`, `experiments/cuda_checker/run_box.sh`,
`experiments/phase2/configs/checker_gpu.yaml`, `experiments/phase2/final_eval.py`, `experiments/phase2/select_ckpt.py`,
`src/objdel/train/serve.py`, `src/objdel/train/rewards.py`, `src/objdel/reward/{score,cache}.py`, the image lock
`experiments/phase2/docker/locks/vllm_native.txt`, and the run records `runs/k1/*`, `runs/r5/box.txt`,
`runs/r5/eval/eval_*/final_eval.json`.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Review only. I installed nothing from the catalog, because an install needs the owner's approval. It sent me to NVIDIA's own pages for two facts. The Compute Sanitizer API guide says the API "will fail with SANITIZER_ERROR_NOT_INITIALIZED if the CUDA driver version is not compatible with the Compute Sanitizer version". The CUDA 13 release notes say R580 is the minimum driver for 13.x minor-version compatibility. Together these gave item 3. |
| engineering:code-review | Checked the round-12 D2 fix line by line: the three-arena mapping, `GC_POISON_SCAN`, and the D5c arena list against `Block`/`Lane` (gpucheck.cpp:131-158) and K1's `Engine`/`scan` (k1poly.cpp:119-160). No new defect in those lines. One stale line ref (item 5). |
| engineering:debug | Followed the RL step's engine calls in order: prefetch, then D7b `add_tables_from_raters`, then `score_polygons`. I compared that order with the D5c rule "sized at its first call from a stated bound" and found item 1. I also traced the V1 `final_eval` pool path (item 4). |
| engineering:system-design | Walked the device-memory budget for RL and eval: the D5c private pool, tables never dropped in RL (`weakref.finalize` only, cache.py:256), the eval scratch cap, and G-MEM / D-G17 (c). The box toolchain pins match the image lock: `nvidia-nvtx==13.0.85` is already in vllm_native.txt:479, so D1's `nvtx` extra causes no conflict. |
| engineering:architecture | Checked the D5c pool decisions (private pool, no `cudaDeviceSetMemPool`, release threshold) and the "pageable for D7b uploads" decision against what the code stages through `Lane::out` today. This gave item 1's fix. |
| engineering:documentation | Each item has: source file and heading, the doc quote, the code or measurement, and one exact fix. |
| huggingface-skills:huggingface-best | No model choice in my area this round. Nothing changed. |
| huggingface-skills:huggingface-papers | D.2 cites CCL, EDT and contour papers. The rest are NVIDIA docs. No claim in my area needed a new paper check. Nothing changed. |

Web checks (2026-11-09):
- docs.nvidia.com/compute-sanitizer/SanitizerApiGuide: "Compute Sanitizer API calls will fail with
  SANITIZER_ERROR_NOT_INITIALIZED if the CUDA driver version is not compatible with the Compute Sanitizer version."
- docs.nvidia.com/cuda/cuda-toolkit-release-notes (13.5): "Existing CUDA 13.x applications run on drivers >=580 under
  CUDA minor version compatibility". The driver branches are 13.0 = R580 and 13.4 = R615.
- docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html: under minor-version compatibility,
  applications run "with limited feature-set". The page names no rule for developer tools.

Re-verified and found correct (no item):
- K1 needs three device blocks, and the round-12 fix adds `scan` to `Lane`. `scan`'s CUDA path grows only `tmp`
  (k1poly.cpp:136-152).
- `GC_POISON_SCAN` covers every byte CUB can write, because CUB temp storage starts at the block base and is smaller
  than 4(n+1) bytes.
- D5c pool setup. `cudaMemAllocationTypePinned` and a device location are correct, and a private pool leaves torch
  and vLLM untouched.
- `score_polygons` makes 4 `wait()` calls (gpucheck.cpp:1574, 1626, 1639, 1802) plus 2 `Drain` syncs (177-184).
- `check_mask_chunk` makes 3 waits (1020, 1027, 1055).
- The 256 answers per reward call match the TRLReward docstring and `score_many`, which makes one
  `score_polygons` call (score.py:112-115).
- The scorer floats come from integer sums (`ObjTerms`, gpucheck.cpp:1488-1516), with no transcendental function. So
  GPU02's Mac-golden bit equality is sound.
- The box image already holds `nvidia-nvtx==13.0.85` (lock:479), so the CMake `find_file` for `nvtx3/nvtx3.hpp` under
  `nvidia/cu13/include` is reachable.
- V1 cap 3.35 h, V3 amounts $2.55 and $3.23, and the D share 79 x 0.641 / 60 = $0.84 all check out.

Not repeated: round 12 minor items 2-4 (D.1 item 9 wording, kernel log paths, runbook line numbers) are still open in
the text. Round 11 minor items are open too.

Counts: blocking 0, major 1, minor 5.

---

## 1. MAJOR. The pinned-staging rule depends on which entry calls first. In RL, D7b's `add_tables_from_raters` runs before `score_polygons`, so `pinned_grows` becomes 1 at step 1 and D-G15 fails with no defect

Sources: D (heading "D5. Scorer labels ...", step 5 "D5c"; heading "D7b. Expert tables built on the device";
heading "D.6", row D-G15; heading "D.4", "Engine counters", `pinned_grows`).

Doc:
- D.md D5c: "So the engine sizes each pinned block at its first call from a stated bound:", then one bound per
  entry: "`score_polygons`: ... 20.5 MB per name", "`score_polygons` output staging (`Lane::out`, one block reused
  by all copies of a call) ...", "`eval_sets`: its scratch cap (D7)." Then: "Each later growth adds 1 to
  `pinned_grows`".
- D.md D.4 counters: "`pinned_grows`: growths of a pinned staging block after its first sizing (D5c)."
- D.md D5c: "The D7b batch upload stages its host bytes in plain host memory (`operator new`, pageable) ... It never
  uses a pinned `Block`". This covers the upload only.
- D.md D7b: "It returns the handles plus the scalars the host needs (perimeter, areas, signatures, voter flags,
  n_vote). One H2D copy and one wait per batch". No sentence says how these outputs reach the host. Also: "change
  `src/objdel/train/rewards.py` `TRLReward.prefetch` ...: before scoring, collect the step's caches that have no
  tables for this pid. Build them in ONE `add_tables_from_raters` call".
- D.md D-G15: "`pinned_grows == 0` over the whole run, read at the end."

Code:
- Every device-to-host result in the engine today goes through the lane's pinned `out` block, which grows on
  demand. Boxes: `ptr<i64>(L.out.get(sizeof(i64) * size_t(5 * n)), 0)` (gpucheck.cpp:1572). Object results:
  `u8* ho = L.out.get(res_pad + ...)` (:1623). Signatures: `to_host(L.out.get(sizeof(i32) * size_t(total)), ...)`
  (:1638). Defects: `L.out.get(c.words - c.defects)` (:1799).
- `Block::get` keeps the buffer only `if (bytes <= cap)`. Otherwise it frees the old buffer and allocates a new one
  (gpucheck.cpp:139-149), at a minimum of 4,096 B (:140).
- A port of the new entry would follow these lines and return perimeter, areas, signatures and voters through
  `L.out`.

Failure: the first engine call in the RL process is D7b's batch at step 1, with 32 cold caches. Its outputs are a
few KB (32 perimeters, 32 x k areas and voter flags, the rater signatures). Read per entry, the D5c rule sizes
`L.out` from that call, to 4,096 B. The step's `score_polygons` call then needs at least the defects, 256 x 13 x 8 =
26,624 B, plus the signatures. So `L.out` grows once: a `cudaFreeHost`, which synchronizes the device, then a
`cudaMallocHost`, and `pinned_grows` becomes 1. This is correct code, yet D-G15 fails in V3, and S13 `test_b04`
reports it on the Mac.
- No $0 check sees it. The D5c test, I08 and the GAP-18 case in I15 call only `score_polygons`. No Mac test runs
  `add_tables_from_raters` and then `score_polygons` in a fresh process.
- The text pushes the reader toward the per-entry reading. Its "`eval_sets`: its scratch cap" bullet makes sense
  only per entry: a global maximum would pin 256 MiB of host memory in every RL process for an entry that RL never
  calls.

Fix (D.md D5c, D7b, D.9):
- D5c: replace the bound list's lead-in with this text. "Each lane's pinned blocks have ONE bound each, set at that
  lane's first engine call, whatever the entry: `in` = the `score_polygons` input bound, `out` = the
  `score_polygons` output bound (both below). No other entry sizes or grows them inside these bounds." Delete the
  "`eval_sets`: its scratch cap" bullet, and add: "`eval_sets` copies its integer results into pageable host
  memory (it runs only in `final_eval`, and D-G15 does not read that process)."
- D7b: add "`add_tables_from_raters` copies its returned scalars (perimeter, areas, signatures, voter flags,
  `n_vote`) into pageable host memory (`std::vector`), inside the batch's one wait, never through `Lane::out`. This
  is the same reason as the upload: the cold-cache count changes from step to step."
- D.9.3 I05 (or GAP-18): add a fresh-subprocess Mac case. Run `add_tables_from_raters` on 32 synthetic caches, then
  `score_polygons` on 256 answers, then 5 rounds of (a cold batch of 0-32 caches plus 256 answers). Pass:
  `pinned_grows == 0`, and the `sync_alloc` and `sync_free` deltas over rounds 2-6 are 0.

## 2. MINOR. D.1 item 5 says the eval CPU share "is not measured", but the R5 records hold it: 0.69 s of 415 s (val) and 3.1 s of 634 s (test)

Sources: D (heading "D.1. Goal and why", item 5).

Doc: "The 2,200-item test eval took 634 s wall, generation included [measured: REPO/runs/r5/eval_test.log]. The
CPU share is not measured."

Measurement: `runs/r5/eval/eval_val/final_eval.json` holds `seconds = {load_and_sample: 414.687, total: 415.378}`.
`runs/r5/eval/eval_test/final_eval.json` holds `{630.456, 633.567}`. In `final_eval.py:76-78`, `score_items` (the
host `metrics()` work on the pool) runs after `t_gen` is taken. So host scoring took 0.69 s (val, 1,084 items) and
3.11 s (test, 2,200 items).

Fix: "Host scoring (`score_items`, after generation) took 0.69 s of the 415 s val eval and 3.11 s of the 634 s test
eval [measured: runs/r5/eval/eval_*/final_eval.json `seconds`]. D7 is not a speed item: its value is the owner rule
and fewer CPU dependencies." D.1 item 9 already says this for RL.

## 3. MINOR. compute-sanitizer 13.4.92 has passed only on a driver-615 box, while the offer filter accepts R580. Add a 10-second sanitizer smoke before R-D2's paid steps

Sources: 00_front (heading "6.1 Before any rental", step 5; "6.4", last paragraph), D (heading "D8", the `sanitize`
function and `--gpu-check` steps 1-4).

Doc:
- 00_front 6.1 step 5: "US/CA host, driver R580 or newer (CUDA 13)".
- 00_front 6.4: "The box toolchain for V2 and V3 is `cuda-toolkit[nvcc,cccl,nvtx]==13.0.3` and
  `nvidia-cuda-sanitizer-api==13.4.92`."

Measurement:
- The only passing sanitizer run is `runs/k1/cuda_checker/memcheck2.txt` ("ERROR SUMMARY: 0 errors", memcheck and
  initcheck). It ran on driver 615.71.09, CUDA 13.4 (`runs/k1/box.txt`).
- The RTX 5090 boxes of R5 had driver 580.159.03, CUDA 13.0 (`runs/r5/box.txt`).
- The Compute Sanitizer API guide says the API fails when "the CUDA driver version is not compatible with the
  Compute Sanitizer version". No NVIDIA page I found says whether a 13.4 sanitizer runs on R580.
- If it does not run, `--gpu-check` reaches step 4 after the build and the pytest minutes, and D-G5 cannot pass.

No measurement proves a failure, so this item is minor.

Fix (D.md D8 `--gpu-check`, step 2): right after the `build() == 'cuda'` assert, add
`timeout 60 "$CUDA_HOME/bin/compute-sanitizer" --tool memcheck --error-exitcode 99 "$PY" -I -c "import numpy as np,
objdel._gpucheck as g; g.check_masks(...)"` on one 8x8 mask. On a non-zero exit, print `nvidia-smi
--query-gpu=driver_version` and stop before step 3. In 00_front 6.1 step 5, add: "when V3 item 1 is approved, prefer
an offer with driver >= 615 (the only driver on which sanitizer 13.4.92 is measured, runs/k1/box.txt)".

## 4. MINOR. Front 6.2 says pool size 1 means "no forked worker opens its own CUDA context", but `fork_pool(1)` forks one worker, and that worker opens a context

Sources: 00_front (heading "6.2 V1", "Why step 0"; also 6.3 V4 "(pool size 1, as A3)").

Doc: "On a CUDA build `final_eval` uses pool size 1, as A2 step 2 does, so no forked worker opens its own CUDA
context (A3 adds the cap)."

Code:
- `fork_pool(n)` builds `_Pool(n, mp_context=get_context("fork"), ...)` and waits for all n workers
  (serve.py:164-178). So pool size 1 is one forked child.
- `final_eval` maps `_score` -> `metrics()` on that pool (final_eval.py:77; select_ckpt.py:172-176, 228). `metrics()`
  calls `to_masks`, and `to_masks` calls `_gpucheck.check_polygons` (D.1 item 5 [code: shape.py:136-148]). On a CUDA
  build, the child therefore opens one CUDA context.
- This is legal, because `fork_pool` refuses to fork after CUDA starts (serve.py:168). It costs one extra context
  beside vLLM, and scoring is serial: about 11 s for 1,084 items [estimate: 16 workers x 0.69 s, item 2].

The sentence is wrong, but nothing breaks.

Fix: "On a CUDA build `final_eval` uses pool size 1 (A2 step 2): exactly one forked worker opens a CUDA context.
This is legal, because the pool forks before the parent starts CUDA (serve.py:166-168). D7 step 6 later removes the
pool."

## 5. MINOR. D2 step 1 cites K1's `Engine` blocks at k1poly.cpp:162; the line is 158

Sources: D (heading "D2", step 1).

Doc: "K1-draft's `Engine` holds three device blocks, `Block<true> devA, devB, tmp;` [code: K1/k1poly.cpp:162]".

Code: `grep -n "Block<true> devA" k1poly.cpp` prints `158:  Block<true> devA, devB, tmp;`.

Fix: change `K1/k1poly.cpp:162` to `K1/k1poly.cpp:158`.

## 6. MINOR. D7b puts the batched table build in `TRLReward.prefetch` without saying which thread runs it. Today table uploads run on the prefetch thread

Sources: D (heading "D7b", second Files bullet).

Doc: "change `src/objdel/train/rewards.py` `TRLReward.prefetch` and the non-prefetch path of `_terms`: before
scoring, collect the step's caches that have no tables for this pid. Build them in ONE `add_tables_from_raters`
call through a helper `reward.cache.device_many(caches)`."

Code:
- Today `prefetch` only builds `jobs` and submits `score_batch` to the one-thread executor (rewards.py:111-119).
- `cache.device()`, the table upload with its engine wait, runs inside `score_many` (score.py:114), on that thread.
- If `device_many` is called in `prefetch`'s own body, its H2D copy, device build and wait move onto the main thread
  before the old-logprob forward, which is the step's critical path.
- The cost is a few ms per step [estimate], so this is minor.

Fix: "Call `device_many` first inside the job that `prefetch` submits to `self._thread` (a small wrapper around
`score_batch`), never in `prefetch`'s own body. The non-prefetch path of `_terms` calls it in line, as today's
`score_many` does."
