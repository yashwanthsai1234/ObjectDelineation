# Final review round 6: CUDA and GPU performance (Implementation.md, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Routing only. No catalog skill covers rental timing, fork safety or pool allocators, and nothing was installed. I checked the CUDA facts against NVIDIA docs directly. The stream-ordered allocator guide confirms D5c: release threshold default 0, `UINT64_MAX` stops shrinking at sync, and `cudaMallocFromPoolAsync(ptr, size, pool, stream)` exists. So D5c raised no item. |
| engineering:code-review | I read the code behind every claim I tested: `gpucheck.cpp` (alloc/release :78-84, `Block` :135-150, `Drain` :177-184, waits :1574/1626/1639/1690/1744/1802, entries :1663-1835), `CMakeLists.txt:15`, `Makefile:41-45`, `run_box.sh`, `policy.py` (`code_files` :374-383, `sync_up` :386-418, `chain` :430-460, launch :315-329, `run` :470-490), `serve.py:125-177`, `final_eval.py:62-79`, `select_ckpt.py:172-231`, `reward/cache.py:241-259`, `reward/score.py:92-115`, `K1/cross_numpy.py:16`, `K1/build_cpu.sh`. |
| engineering:debug | I followed each paid path step by step from launch to teardown: V1 (`b1_rows` step 0 to item 5), V2 (item 1 to items 2-4), V3 (item 0 to item 4) and V4. I also followed the D7 step 6 fork guard into A2's ceiling CLI. That gave items 1, 2, 4 and 5. |
| engineering:system-design | I checked GPU placement after each phase: one CUDA context per process, fork before init, pool workers, and engine pool versus torch and vLLM memory. Item 6 is a wording error, not a placement error. The pinned bounds in D5c and the pool memory hold up. |
| engineering:architecture | I checked that rules agree across phases: the D7 step 6 guard and CLEAN16 against the A2 fork pool and I19 (item 1); the time line step 5 against the D1 step 4 R-D1 order (item 5); and front 6 "the one place for launch configs" against V2-V4 (item 3). |
| engineering:documentation | Every blocking and major item quotes the source file and heading, plus the code line or ledger row that proves it. |
| huggingface-skills:huggingface-best | No model choice is in scope (Qwen3-VL-2B is fixed). It changed nothing. |
| huggingface-skills:huggingface-papers | No paper claim is in my area this round. It changed nothing. |

Web checks (WebFetch and Exa loaded): CUDA Programming Guide §4.3.4.2 (release threshold, `UINT64_MAX`, `cudaMemPoolTrimTo`). CUDA Runtime API 13.4, Stream Ordered Memory Allocator (`cudaMallocFromPoolAsync` signature; `cudaMemPoolAttrReservedMemCurrent`). Both confirm D5c as written. The runtime "API synchronization behavior" page says nothing about allocation calls. The doc's statement that `cudaFreeHost` synchronizes rests on the code comment and common practice. No item depends on it.

Already fixed and not repeated: r5 cuda-1 (kernel drafts ship through the `code_files` keep rule), cuda-2 (the V2 bundle command names both train views), cuda-3 (D-G15 now gates `sync_alloc`/`sync_free` and records the pool counters), cuda-4 (V4 plan file `e_v4.yaml` with a CUDA build step 0), and r4 cuda-1 to cuda-4.

## Items

### 1. BLOCKING: D7 step 6 breaks A2's fork-pool ceiling CLI and test I19, so `make test` fails when D7 lands

- Sources: D (D7 step 6; D.9 / D-G17), A (A2 step 2; A.9.3 I19), CLEAN (CLEAN16 step 1).
- Doc, D.md "D7. Evaluation pixel work on the device", step 6: "`ExpertCache.device()` raises `RuntimeError("device() in a forked child")` when `os.getpid()` differs from the pid recorded when `objdel.reward.cache` was imported. A forked child inherits the parent's record, so it fails loudly". The same step converts only "`final_eval`, `select_ckpt` and `g_probe`".
- Doc, A.md "A2. Fair ceiling", step 2: "Parallelise over items with `serve.fork_pool(serve.pool_size())`, as `final_eval` does ... Each worker uploads its own cache tables (`ExpertCache.device` is per process)". A.md A.9.3 "I19 ceiling pool equals serial": "The same items through `serve.fork_pool` with pool size 1 and 2: `ceiling_items.jsonl` equal line for line (each worker uploads its own cache tables)".
- Doc, CLEAN.md "CLEAN16", step 1: it deletes "`serve._Pool`, `fork_pool`, `pool_size` ..." and edits only "`select_ckpt`, `g_probe`, `final_eval`". It does not mention the ceiling CLI.
- Code: `fork_pool` imports `SCORERS` (so `objdel.reward.cache` too) in the parent, then forks with `get_context("fork")` [src/objdel/train/serve.py:170-173]. Scoring calls `cache.device()` in each worker [src/objdel/reward/cache.py:241-257; reward/score.py:114]. So every forked ceiling worker has the parent's recorded pid, and the guard raises on the Mac CPU build too: it checks the pid, not the build.
- Effect: once D7 step 6 lands, I19 fails, `make test` is red, and the repo rule ("`make test` must pass") stops D7 from merging. After CLEAN16, `experiments.phase2.ceiling` cannot even import `fork_pool`. Two rules both claim to be binding: the A2 parallel path and the D7 guard.
- Fix: in D7 step 6, add `experiments/phase2/ceiling.py` (and `eval/ceiling.py`) to the list of callers that stop using the pool. The ceiling then scores items in-process through `eval.device.item_counts` / `score_many` in batches. In A2 step 2, add: "until D7 step 6; after it, in-process". Replace I19 with "ceiling in-process equals the pre-D7 pool output (golden `ceiling_items.jsonl` from the A2 commit)". In CLEAN16 step 1, add `ceiling` to the list of `--procs` and pool users.

### 2. MAJOR: V2, V3 and V4 leave out the HF downloads that their training and eval steps read

- Sources: 00_front (6.3 "V2, V3 and V4"), C (C0 step 3, C1 `c_speed_check`), B (B1.5 "Rows-root download", B3.9), E (E.8 V4).
- Doc, 00_front.md 6.3, V2 item 1: "Build the V2 bundle on the Mac ...; launch, sync, `env_check`, toolchain and CUDA build". There is no `fetch_hf.sh`. But the C0 baseline arm is `train_rl --ckpt {rd}/ckpt` (A.md A6 step 8), and C.md C1 says `c_speed_check` must "Load the same `--ckpt` as the C0 baseline arm".
- Doc, 00_front.md 6.3, V3 item 0: "Launch, sync, `fetch_hf.sh <run> rl-full-3ep` ... and `fetch_hf.sh <run> b1-rows --as b1`". There is no `fetch_hf.sh <run> --rows <rows_revision> lidc-train-b1`. But B.md B3.9 says `b3_v5.yaml` uses "the same B1 rows root and `rows_revision`". B1.5 says: "A plan step on the box cannot download: the token reaches the box only on the stdin of an ssh command from the Mac [code: policy.py:524, 536]".
- Doc, 00_front.md 6.3 V4 and E.md E.8 V4 step 3: "the E4 sampler arm: `final_eval` ... K 8" needs a `--ckpt`, and no fetch is listed.
- Code: `cmd_fetch_hf` is the only way a checkpoint or the rows root reaches the box [experiments/phase2/vast/policy.py:527-536]. The LIDC bundle holds neither [bundle.py:52-70].
- Effect: an engineer who follows the 6.3 tables gets a V2 whose `c_speed_check` and C0 arm fail at model load. V3 item 4, the B3 run ($1.57), refuses its missing rows root after the build and start-up minutes. V4's eval has no checkpoint.
- Fix: give V2, V3 and V4 launch-line blocks like the V1 block in 6.2. V2: `sync_up.sh <run> && fetch_hf.sh <run> rl-full-3ep`. V3: also add `fetch_hf.sh <run> --rows <rows_revision> lidc-train-b1`, with its manifest SHA-256 check. V4: `fetch_hf.sh <run> rl-full-3ep` (or the checkpoint E4 names). Count the downloads in each item-0 basis.

### 3. MAJOR: V2, V3 and V4 have no stamped launch config, which `launch.sh` requires

- Sources: 00_front (6 "One source"; 6.1 step 4; 6.3), B (B3.9), E (E.8).
- Doc, 00_front.md section 6, "One source": "This section is the only place that lists ... the launch configs that `make preflight` stamps". 6.1 step 4 names only `b1_launch.yaml` ("the only V1 launch config") and says that plan files "never follow `preflight CONFIG=`". 6.3 names plan files only (`c_speed.yaml`, `r_d1.yaml`, `d_gpu_test.yaml`, `a3_test.yaml`, `b3_v5.yaml`, `e_v4.yaml`).
- Code: `cmd_launch` refuses unless `preflight.refusal(Path(a.config))` passes [experiments/phase2/vast/policy.py:315-316]. `preflight.load` builds `RunArgs(**raw)`, so a `steps:` plan raises `TypeError` (front 6.2 "Why a new field and not a plan file").
- Effect: the engineer must guess the V2, V3 and V4 launch config. For V3 the guess changes the result. The B3 run trains on the B1 rows root, and only a `rows_plan` launch config puts the rows-root SHA-256 and the "ROWS root" quality check into the stamp (front 6.2 "How the rows root enters the stamp"). A guessed `r1_qwen3vl2b.yaml` launches B3 with no rows-root check.
- Fix: in 6.3, name one launch config per rental. V2: `c_launch.yaml` (RunArgs, `view: views/lidc-bundle-v3`, no `rows_plan`). V3: `b3_launch.yaml` (RunArgs, `rows_plan: experiments/phase2/configs/b3_v5.yaml`, so its `train_rl` stage's rows SHA enters the stamp). V4: `e_v4_launch.yaml`, or the E6 `RunSpec` path. Add each to the 6.1 step 4 check and to the B.9 I19-style `preflight.load` test.

### 4. MAJOR: V4's 30-minute overhead is below the measured launch-to-first-`run` time, so the cap can end V4 before its eval

- Sources: 00_front (6.3 "V4"; section 6 inputs), E (E.8).
- Doc, 00_front.md 6.3 V4: "V4 is its own rental: 30 min overhead + 3 + 9 = about 42 min = $0.58 at r_s [estimate]. Why 30 min: launch to first sync alone took 17-31 min". The front's own inputs table says: "Launch to first `run` | 17.7 min (rl), 24.1 min (r4), 38.2 min (r5). It holds sync and the downloads". V2 item 1 and V3 item 0 both use 38.2 min.
- Measurement: ledger.csv rows 131 and 133: r5 `launch` 2026-11-08T15:35:45Z, first `run` 16:13:57Z = 38.2 min.
- Code: the cap from $0.58 is `floor_to_0.05((0.58 - 0.10) / 0.680)` = 0.70 h = 42 min. `run.sh` ends every job at `min(cap, deadline - now - MARGIN_S)` = launch + cap [policy.py:480; PULL_S = MARGIN_S = 600, policy.py:73].
- Effect: on an r5-like host the plan starts at minute 38.2 and has 3.8 min for build (3), `env_check` (1), the box tests (about 2) and the 415 s sampler eval. The eval never runs and the $0.58 buys no E-G10 result. The 9 min of checks also leave out the untimed eval-process imports, which the front counts elsewhere as about 2.6 min per eval process (6.2 notes: "about 5.2 min for the untimed imports of the two eval processes").
- Fix: price V4 on the same basis as V2 and V3: 38.2 (launch to first `run`, downloads included) + 3 build + 1 `env_check` + 2 tests + 7 eval + 2.6 imports, about 54 min. That is about $0.71 at r_s ($0.10 + 54 x $0.01133). Update the one V4 price in 6.3 and the "Totals" line. E.8 and 99_back point to it.

### 5. MAJOR: R-D1 (V2 item 4) is written for a pre-D2 tree, but the time line lands D2-D8 before V2; on that tree its first step fails

- Sources: 00_front (section 1 time line steps 5-6; 6.3 V2 item 4), D (D1 steps 0 and 4; D3 files).
- Doc, 00_front.md section 1 "Time line": "5. C1-C5, B2/B3 and D2-D8 code land on the Mac, on pre-E paths. 6. V2 and V3 run only after the owner raises the budget". The flow diagram says the same: "Mac: C1-C5 code, B2/B3 code, D2-D8 code" comes before "V2 ... R-D1 (K1, K2, bench)".
- Doc, D.md "D1", step 0: "`r_d1.yaml` step 1 runs `build_cpu.sh /workspace/src/cuda` ..., then `/opt/venv/bin/python -I cross_numpy.py build/cpu/k1_dump_cpu {rd}/k1/lidc --masks {data}/checker_eval/masks.npz`". Step 4: "If D1 has not merged before that rental, run only the two runbooks". This assumes V2 may come before D1 merges, which contradicts time line step 5.
- Doc, D.md "D3", files: "Drop `bitquads` from the public exports" and "move `src/objdel/checker/bitquad.py` to `tests/oracles/bitquad_oracle.py`". The list of oracle callers to fix does not include `docs/implementation2.0/kernels/K1/cross_numpy.py`.
- Code: `docs/implementation2.0/kernels/K1/cross_numpy.py:16` reads `from objdel.checker import bitquads, polygonize`. `run.sh` joins plan steps with `&&` [policy.py:490], so one failed step stops the chain.
- Effect: on the time-line tree, `r_d1.yaml` step 1 raises `ImportError`, and the K1 runbook, the K2 runbook and the bench never run. Even with the import fixed, `polygonize` and `check_masks` are by then the K1 port (D2/D3), so `cross_numpy` and `test_k1`/`bench_k1` (K1 against `check_masks` [K1/build_cpu.sh:2-4]) compare K1 with itself. The $0.37 buys no D-G1 evidence.
- Fix: choose one rule and write it in both places. Either (a) R-D1 runs only on a tree before D2 (owner tag `v2-pre-d2`, with C code but no D2/D3), or (b) R-D1 is dropped whenever D2 has merged before V2. Under (b), D-G1 and D-G2 move to D-G3 in R-D2, as D1 already allows for a dropped R-D1. If (a) is chosen, also add `cross_numpy.py` to D3's list of oracle callers (`from tests.oracles...`) so that the Mac half of K1 GPU_RUN still runs after D3.

### 6. MINOR: "pool size 1, so no forked worker opens its own CUDA context" is false

- Sources: 00_front (6.2 "Why step 0"), A (A2 step 2; A3), E (E.8 V4 step 3).
- Doc: "On a CUDA build `final_eval` uses pool size 1, as A2 step 2 does, so no forked worker opens its own CUDA context".
- Code: `fork_pool(1)` still forks one `ProcessPoolExecutor` worker [serve.py:173-176]. `score_items` maps `_score` over it [select_ckpt.py:231], and `metrics()` calls `to_masks` (`_gpucheck.check_polygons`) and `cache.device()` in that worker [eval/metrics.py:255; geometry/shape.py:143]. So the worker opens a second CUDA context. This works only because the fork comes before the parent initializes CUDA [final_eval.py:69 before `make_trainer`]. The result is correct (vLLM is at 0.35 utilization [configs.py:192]), so this item is minor.
- Fix: write "pool size 1: one forked worker, forked before CUDA starts, so at most one extra CUDA context (about 0.5 GB)". Or, better, pass the pool as `None` on a CUDA build and score in the main process. D7 step 6 later makes that the rule anyway.

### 7. MINOR: the `BudgetStop` rule text does not match the V1 value it computes

- Source: 00_front (6.4 rule 5).
- Doc: "Its value is: approved amount - $0.35 (the 26-min eval tail plus the cap rounding) - the $ of every later item in the same rental. V1: $2.40 - $0.35 = $2.05". V1 items 3-5 ($0.24) come later in the same rental, but they are inside the $0.35 tail and are not subtracted again.
- Fix: write "approved amount - $0.35 (tail and rounding, which already holds the eval items of the training item's own plan) - the $ of every later item outside that tail".

Totals: 1 blocking, 4 major, 2 minor.
