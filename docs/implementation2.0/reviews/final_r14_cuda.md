# Final review round 14: CUDA and GPU performance (Implementation.md)

Reviewer role: CUDA / GPU performance. I checked kernels, GPU placement, syncs, memory, Blackwell sm_120 and the
rental plan.

Scope read:
- 00_front sections 1 (time line, end state, order rules), 2.4, 2.5, 5 and 6 (6.1-6.4).
- C.md C0-C6b and C.6.
- D.md D.4 counters, D1, D2, D4, D4b, D5 (D5c), D6, D7, D7b, D8, D.6, D.8, the D.9.3 rows I05/I08/I15/I19/I20, and
  the round-13 review log.
- A.md A2 step 2, A3 (the `final_eval` changes) and A6 steps 8-9 (`nsys_capture.sh`, `nsys_baseline.py`).
- E.md E4. The round-13 CUDA review, its verdict, and every round-13 log table.

Code and records opened:
- `cuda/gpucheck.cpp` (`alloc`/`release` :78-84, `Block` :131-150, `Lane` :152-156, `Engine`/`engine()` :158-172,
  `Drain` :177-184, `score_words` :1560-1643, `check_polygons` :1663-1693, `add_tables` :1720-1763,
  `score_polygons`/`score_masks` :1784-1832), `cuda/gpucheck_py.cpp` (GIL release), `CMakeLists.txt`, `Makefile`.
- `docs/implementation2.0/kernels/K1/k1poly.cpp`, `K2/k2.cpp`, `K2/k2score.inc`: no hidden `cudaMalloc`, thrust or
  device-wide sync.
- `experiments/cuda_checker/run_box.sh`, `experiments/phase2/{final_eval,select_ckpt,train_rl,g_probe}.py`,
  `src/objdel/train/{serve,rewards,configs}.py`, `src/objdel/reward/cache.py`, `src/objdel/geometry/{text,shape}.py`.
- `.venv trl/generation/vllm_generation.py` (colocate `external_launcher`).
- `runs/r5/{env_check_vllm.log, r5_shared.log, eval/*}`, `runs/k1/checker_gpu.log`.
- `data/prepared/lidc/train.jsonl` (first record, read with `python3 -I`) and one LIDC image header (`file`).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Review only. I installed nothing, because an install needs the owner's approval. No catalog skill covers a custom engine's pinned staging or stream-ordered pool. It sent me to the NVIDIA pages listed under "Web checks". |
| engineering:code-review | I checked the round-13 D5c/D7b fix line by line against `Block::get`, every `L.in.get` and `L.out.get` site, and `check_polygons`/`score_polygons`/`score_masks`. The bound list matches the code: defects :1797, boxes :1572, `ObjOut`+`sig_start` :1622-1623, signatures :1638, inputs `c.in_bytes`/`a_rbox`/`o_parts`. No defect there. Two smaller gaps became items 1 and 3. |
| engineering:debug | I traced the RL reward path (prefetch thread -> `score_batch` -> `score_many` -> `score_polygons`, GIL released [code: gpucheck_py.cpp:127-220]) against `FastGRPO.log`, which reads `gc.stats()` on the trainer thread. This gave item 2. I traced the V1 `final_eval` path (fork before CUDA, `_score` in the worker) for item 4. |
| engineering:system-design | I walked the device-memory budget with real crop sizes. LIDC caches are crop-sized (record `H = W = 50`), the image is 512 x 512, and tables are tens of KB, so RL tables that are never dropped are not a risk. I checked the eval topology before and after E4 (in-process colocate vLLM vs `VllmSampler`'s EngineCore child, item 5). |
| engineering:architecture | I checked the D5c pool decisions again (private pool, `cudaMemAllocationTypePinned`, release threshold, no `cudaDeviceSetMemPool`). Also the "pageable for D7b" rule against the CUDA runtime sync rules (item 3). |
| engineering:documentation | Each item has the source file and heading, the doc quote, the code or measurement, and one exact fix. |
| huggingface-skills:huggingface-best | No model choice in my area. Nothing changed. |
| huggingface-skills:huggingface-papers | No new paper claim in my area this round. Nothing changed. |

Web checks (2026-11-09):
- docs.nvidia.com/cuda/cuda-runtime-api, "API synchronization behavior". Two rules: "For transfers from pageable
  host memory to device memory, a stream sync is performed before the copy is initiated", and "For transfers from
  device memory to pageable host memory, the function will return only once the copy has completed". Used in item 3.
- docs.nvidia.com/cuda/cuda-toolkit-release-notes (13.5): "Existing CUDA 13.x applications run on drivers >=580
  under CUDA minor version compatibility". Branches: 13.0 = R580 and 13.4 = R615.
- docs.nvidia.com/compute-sanitizer/ReleaseNotes. One known issue covers only drivers 615.00-615.70 with apps built
  before 13.4. It says nothing for R580. So round-13 minor item 3 (sanitizer on R580) stays a minor item with no new
  proof, and I do not re-raise it.

Re-verified and found correct (no item):
- D5c bounds. The `in` bound (256 x names x 1000 x 10 x 8 B = 20.48 MB per name) covers `c.in_bytes`, plus the
  `Obj` arrays that `score_words` stages through `L.in` (:1566, :1588), which are about 45 KB for n = 256. The `out`
  bound list matches the four `L.out.get` sites. `Block::get` rounds 20.48 MB up to 32 MiB (:140-142).
- I20's control case is reachable. A call of about 700 answers needs 72,800 B of defects, more than the 65,536 B
  `out` block, so `out` grows exactly once.
- After D5c, the engine and the K1/K2 drafts make no allocation outside `alloc`/`release`. The only CUB call is
  `DeviceScan::ExclusiveSum` with caller scratch (k1poly.cpp:140-142).
- The cudart link is static (`CUDA::cudart_static`, CMakeLists.txt). So the D-G7 old-engine venv needs only the
  driver at run time.
- TRL colocate builds vLLM with `distributed_executor_backend="external_launcher"` (vllm_generation.py:364). So
  training and today's `final_eval` sample in-process. Nsys `cudaProfilerApi` capture, `count_syncs()` and the C5
  one-call `load_weights` all see the same process.
- The LoRA config excludes the vision tower (`exclude_modules=VISION`, configs.py:95). So C3's `no_grad` encode is
  sound.
- V1 and V3 money: V1 193.0 min = $2.29 and cap 3.35 h. V3 B-S2 195.6 min = $2.32, `BudgetStop` $2.20, amount
  $2.55. Items 0/1/3/4 give $2.88 / $3.23. The tails fit in $0.35 (30.9 min) at any offer inside the limit.
- vLLM colocate `gpu_memory_utilization` is 0.35 (configs.py:192). So V1's one forked eval worker (a second context
  of about 0.5 GB) fits.

Not repeated: round-13 minor items 2-6 and the open round-11/12 minor items.

Counts: blocking 0, major 0, minor 5.

---

## 1. MINOR. D5c sizes the pinned blocks from "names per answer" and C_max, but no line says how these values reach the C++ engine at its first call

Sources: D (heading "D5. Scorer labels ...", step 5 "D5c").

Doc (D.md D5c):
- "The engine sizes `in` and `out` of every lane to these bounds once, inside the first engine call of the process,
  whatever that entry is".
- "`score_polygons`: 256 answers x names per answer x `MAX_VERTICES` x `MAX_RINGS` x 8 bytes".
- "the bound uses C_max = 64 components per object [estimate], or twice the largest per-object count in the D1
  replay set (`sig_start` steps in `runs/d1/bench.json`), whichever is larger".

Code:
- The engine gets no dataset or batch size at a call. `add_tables_from_raters` (the first RL call) has no answer
  count. `score_polygons(const Polygons& in, const i64* tables)` knows only its own n (gpucheck.cpp:1784).
- `names` depends on the dataset (`SCHEMAS = {"riga": ("disc", "cup"), "lidc": ("nodule",)}`, text.py:58).

Why minor: in RL (LIDC only, 256 answers) every plausible choice (names 1 or 2, C_max 64) passes I20 and D-G15.
A wrong guess grows a block only in a RIGA or offline process, and D-G15 does not read those processes.

Fix (D.md D5c): add one sentence after "The bounds:". "The bounds are compile-time constants in `gpucheck.cpp`:
`kStageAnswers = 256`, `kStageNames = 2` (the largest `SCHEMAS` entry), `kStageCompPerObj = 64`. D1 step 2 checks
`kStageCompPerObj >= 2 x` the largest per-object count of the replay set, and a unit test asserts
`kStageNames == max(len(v) for v in SCHEMAS.values())`."

## 2. MINOR. `gc.stats()` is read on the trainer thread while the reward thread is inside an engine call with the GIL released, and the doc does not say the counters are atomic

Sources: D (heading "D.4. Tasks", "Engine counters"; heading "D7", step 5; D.9.3 row I15).

Doc:
- D.4: "a process-wide `Stats` struct in `gpucheck.cpp`, outside the `__CUDACC__` split, counts these events".
  Nothing says atomic or under `e.mu`.
- D7 step 5: "`FastGRPO.log` ...: at the end of steps 1-6 it logs `engine/<key>` for every `gc.stats()` key ...
  Both are host reads of counters, so they add no CUDA sync".
- I15: "4 Python threads x 50 `check_masks` calls ... while another thread calls `trim(0)` and `reset_stats()`:
  ... counters non-negative".

Code:
- Every entry releases the GIL (`nb::gil_scoped_release`, gpucheck_py.cpp:127, 147, 186, 204, 220).
- The reward runs on its own thread (`ThreadPoolExecutor(1, thread_name_prefix="reward-prefetch")`,
  rewards.py:118).
- So `stats()` and `reset_stats()` can run at the same time as the helpers that increment the counters.
- If `Stats` holds plain `i64` fields, this is a C++ data race, which is undefined behaviour. If an engineer adds
  `e.mu` to `stats()` to fix it, `log` blocks for a whole engine call on the trainer's path. That is a hidden wait
  that D-G8 does not count.

Why minor: on x86-64, aligned 64-bit loads do not tear. The D-G15 deltas would still be right.

Fix (D.md D.4 counters): add "Every `Stats` field is `std::atomic<i64>`. Updates and reads use
`memory_order_relaxed`. `stats()` and `reset_stats()` never take `e.mu`. `trim()` takes `e.mu`, because it must
not run during a `Block` growth. `stats()` reads `pool_reserved` and `default_pool_reserved` only after
`e.ready`, and reports 0 before, so a `stats()` call never creates a CUDA context."

## 3. MINOR. D1 says no host sync escapes the counters, but D7b's pageable copies (and today's `add_tables`) block the host outside `wait()`

Sources: D (heading "D.4. Tasks", "Engine counters", `drain_syncs` bullet; heading "D7b", Files, first bullet;
heading "D5", step 5).

Doc:
- D.4: "D1 routes it through a counted helper, so no host sync escapes the counters."
- D7b: "`add_tables_from_raters` copies the returned scalars ... into pageable host memory (`std::vector`) after the
  batch's one wait, never through `Lane::out` ... A device-to-pageable copy 'might be synchronous with respect to
  host' ... the copy adds only its own time".
- D5c: "The D7b batch upload stages its host bytes in plain host memory (`operator new`, pageable)".

Code and docs:
- `to_device`/`to_host` are `cudaMemcpyAsync` (gpucheck.cpp:87, 91). `add_tables` copies from a `std::vector`
  today (:1726, :1743).
- The CUDA runtime API says that a pageable H2D copy does "a stream sync ... before the copy is initiated". A D2H
  copy to pageable memory "will return only once the copy has completed".
- So each D7b batch blocks the host three times: H2D, `wait`, and the D2H of the scalars. `waits` counts one.

The design choice is right: the stream is idle, and no `cudaMallocHost` happens. Only the claim and the I05 count
are not exact.

Fix: in D.4, replace "so no host sync escapes the counters" with "so every `cudaStreamSynchronize` is counted.
Pageable copies (`add_tables`, `add_tables_from_raters`, `eval_sets`) block the host on an idle stream and are not
counted". In I05, add "per batch: `waits` +1, and no `sync_alloc` or `sync_free`".

## 4. MINOR. A3's reason for pool size 1 is wrong: the one forked worker still opens its own CUDA context

Sources: A (heading "A3", the `final_eval` bullets), 00_front (heading "6.2", "Why step 0" note).

Doc:
- A.md A3: "on a CUDA build ... `final_eval` caps the pool at 1 worker ..., as A2 step 2 does, so no forked worker
  opens its own CUDA context."
- 00_front 6.2 says the same: "`final_eval` uses pool size 1 ... so no forked worker opens its own CUDA context".

Code:
- `fork_pool(n)` forks n workers (serve.py:173-176).
- `score_items` maps `_score` over the pool, "one CUDA checker call per item" (select_ckpt.py:224-230).
- So with n = 1, the scoring runs in the one child, and that child opens a second context next to the parent's
  vLLM and trainer context.

The cap is still useful: it gives 1 extra context, not up to 16. It is safe, because the fork happens before the
parent initialises CUDA (serve.py:167-169). The vLLM utilisation of 0.35 (configs.py:192) leaves room.

Fix (A.md A3 and 00_front 6.2): "... caps the pool at 1 worker, so at most one forked worker opens a CUDA context
(one extra context, about 0.5 GB, not up to 16). This holds until D7 step 6 moves scoring into the main process."

## 5. MINOR. D-G17 (a) states two different pass rules for the same check ("stays equal" vs "must not rise")

Sources: D (heading "D.6", row D-G17).

Doc:
- D-G17 (a) nvidia-smi method: "If it is >= 1, the count, sampled every 5 s, stays equal to it: no new context
  appears."
- D-G17 (a) /proc method: the count "must be >= 1 before scoring and must not rise during scoring."

Measurement: a plain `vllm.LLM(...)` on the box runs its engine in a child process (`(MainProcess pid=2130)` vs
`(EngineCore pid=2697)`, runs/r5/env_check_vllm.log:1-7). E4's `VllmSampler` is "`vllm.LLM` offline, in process"
(E.md E4 step 1). So on the post-E path (a raise that merges `after-raise` after E0), `final_eval` has a second GPU
process. That process may exit during the last view's scoring, if the sampler is freed first. Under the
"stays equal" rule, a drop fails D-G17 on correct code. Under the /proc rule, a drop passes.

Why minor: E4 keeps one `VllmSampler` for the whole view loop, so a drop is not expected. No measurement shows one.

Fix (D.md D-G17 (a)): replace "stays equal to it" with "never rises above it (a drop, for example a vLLM engine
process that exits, is allowed)". Then both methods use one rule.
