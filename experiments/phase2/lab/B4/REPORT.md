# B4: reward serving speed (H-A6-4 preload, H-A6-2 forked pool, H-K3SYN, H-A6-3 decision)

Lab agent B4, wave B, 2026-10-19. CPU and Docker only; the Mac GPU was not used.
- Lab code: `experiments/phase2/lab/B4/` (ruff-clean).
- Raw outputs: `<scratchpad>/phase2/B4/{mac,docker/io,docker2/io}/*.json` and `*.log`.

## 0. Pre-registration

The criteria are in `prereg.py` (originally `criteria.py`, mtime 18:02; the first timed run started about 18:10). Every script checks its numbers against that file.
- **Calibration done beforehand, no timing:** `riga_stats.py` computes aggregate quantiles of the 64-item RIGA K3 prompt set.
- **Workloads:** K3 batches from `reward_gates.make_jobs`, 32 prompts × 8 completions = 256. Three mixes: RIGA-only, mixed (40 % RIGA) and LIDC-only, all on `data/prepared` train caches.
- **Wall time:** the warm second run of `advantages_v4(score_batch(...))`.
- **Hash:** `reward_gates.digest`.
- **Added after the coordinator's message:** CPU time (this process plus the pool workers) is logged alongside wall time for every pool and serial run.

| ID | Criterion |
|---|---|
| P4-1 | Preload of all RL train caches ≤ 60 s (Mac, 2 threads) |
| P4-2 | Memory after preload ≤ 6 GB |
| P4-3 | 0 cache loads per step with preload; LRU-64 path ≥ 28 misses per step |
| P4-4 | LRU step reward call ≥ 1.5 s slower than preloaded (median) |
| P4-5 | Preloaded `TRLReward` call / plain `score_batch` ≤ 1.10 |
| P4-6 | Rewards identical between the LRU and preloaded paths |
| P2-1 | Mac RIGA-only speed-up at 4 workers ≥ 3.0 |
| P2-2 | Mac RIGA-only best pool median ≤ 125 ms |
| P2-3 | Hashes identical to serial everywhere |
| P2-4 | Worker USS (private memory) growth < 200 MB |
| P2-5 | Pickling and pipe overhead ≤ 10 % of the 4-worker wall |
| P2-6 | TRL call pattern: 1 pool dispatch per generation batch, outputs `==` serial, 10 steps |
| P2-7 | Docker relative speed-up: 2 workers ≥ 1.6, 4 workers ≥ 2.2 |
| PS-1 | Synthetic/real RIGA-only median ratio in [0.90, 1.10] and its bootstrap 95 % CI within [0.85, 1.15] |
| PS-2 | Mixed ratio in [0.90, 1.10] |
| PS-3 | \|Δ p99 per rollout\| ≤ 1 ms |
| PS-4 | Docker RIGA-only ratio in [0.90, 1.10] |
| PS-5 | Synthetic corpus digest identical on the Mac and in x86 Docker |
| H-3 rule | The C++ scorer is needed if P2-1 or P2-2 fails on the Mac, or if the projected x86 RIGA-only wall at 4 workers exceeds 125 ms |

**Load caveat.** The 8-core M3 was shared with B1 (MPS training), B3 and B5 throughout. Load averages were 11–48; the best window was 7.5 at the start of a run, rising to 12 during it. Every number below states its load. Absolute Mac times varied up to 3× within minutes, so speed-ups are judged from serial and pool runs interleaved per batch. Docker ran linux/amd64 under Rosetta ("VirtualApple @ 2.50GHz", `--cpus 3 --memory 4g --network none`). Docker timings are relative only.

## 1. H-4: preload every prepared cache

Prototype: `serve.preload` and `preload.py`. The test replays 10 TRL-like steps (a seeded permutation of RIGA train items, 32 prompts per step) through the real `TRLReward`, calling its 8 functions in TRL's order with TRL's keyword arguments.

| | Mac (load 23→25, 19:29–19:37) | Docker (emulated; 200 RIGA items + all LIDC) |
|---|---|---|
| Preload wall time | **72.1 s** (RIGA 547 items: 42.0 s; LIDC 7,341 items: 30.1 s), 2 threads | 122 s |
| Memory | Array bytes measured: RIGA 6.27 MB per item = 3.35 GB, LIDC 0.046 MB per item = 0.33 GB. Mac RSS read 0.09–0.54 GB because macOS compresses memory, so it is meaningless here | max RSS 2.32 GB; extrapolated to all 547 RIGA items ≈ 4.2–4.5 GB on Linux |
| Misses per step | LRU-64: **32 / 32** every step; preloaded: **0 loads** | the same |
| RIGA step reward call, median | LRU **4.79 s** vs preloaded **0.65 s** (Δ **4.20 s**) | Δ 2.78 s |
| LIDC step, median | LRU 0.28 s vs preloaded 0.11 s | — |
| Preloaded `TRLReward` / plain | 1.075 (RIGA), 0.955 (LIDC) | 0.96 |
| Same rewards | yes, all 20 steps; aggregator == `advantages_v4` | yes |

**Verdict: accept.**
- P4-1 **fails** (72 s > 60 s at load 23 with 2 threads); P4-2 through P4-6 pass.
- P4-2 is judged on the Linux extrapolation, because macOS RSS is not measurable under compression.
- The cost of misses is real: the LRU of 64 cannot hold 547 RIGA items drawn in permuted order, and the measured cost is 2× A6's ~2 s estimate.
- Preload cost is about one minute once per run. On Vast, use 4 or more threads.
- Optional later step (Phase 3): store the band distance arrays in an uncompressed, digest-checked memmap. Preload then becomes page-cache reads, and any start method shares the pages.

## 2. H-2: process pool forked after the preload

Prototype: `serve.py`.
- `CacheRef(key)` is a read-only mapping that resolves through a module-level store and pickles as its key alone.
- `fork_pool(n)` is a `ProcessPoolExecutor` on the fork context, warmed immediately.
- `forkserver_pool` is the alternative, using `serve_preload` inside the fork server.
- `TRLReward` itself is unchanged: `TRLReward(serve.ref_loader, stop_ids, pool=pool)`.

**Why M6's pool lost.** A RIGA-only batch pickles to **96.2 MB** when jobs carry `ExpertCache` objects, and each worker re-runs the EDTs on unpickling. With `CacheRef` it pickles to **211 KB** (dump 0.1–0.3 ms when the machine is not loaded); results are 52 KB.

### Mac speed-ups

Speed-up is median serial / median pool; "paired" is the median of per-batch ratios. CPU ratio is (parent + workers CPU) / serial CPU.

| Window | Workers | RIGA-only | Mixed | LIDC-only | CPU ratio (RIGA) |
|---|---|---|---|---|---|
| Run 1, load 21–22 (overlapped my Docker run) | 2 | 1.61 (paired 1.62) | 1.75 | 1.55 | — |
| | 4 | 2.06 (paired 2.91) | 2.24 | 2.47 | — |
| | 8 | 2.83 (paired 2.70), pool median **50.9 ms** | 3.56 | 2.79 | — |
| Run 2, load 7.5→12, fork only, no Docker | 2 | 1.94 (paired 2.05) | 2.14 | 1.50 | 0.99 |
| | 4 | **2.46 (paired 2.47)**; serial 546 ms, pool 222 ms | 2.77 | 2.44 | 1.05 |
| | 8 | 2.87 (paired 2.81); serial 176 ms, pool **61 ms** | 2.57 | 2.29 | 1.70 |

### Docker (emulated, 3 CPUs)

| Run | Workers | RIGA | Mixed | LIDC | CPU ratio |
|---|---|---|---|---|---|
| Run 1 | 2 / 4 / 8 | 2.00 / 2.47 / 1.62 | 1.73 / 2.40 / 1.90 | 1.49 / 2.43 / 2.41 | — |
| Run 2 | 2 / 4 / 8 | 2.00 / 3.04 / 2.04 | 2.07 / 3.18 / 1.82 | 2.55 / 3.10 / 1.51 | 2 workers 1.0–1.15; 4 workers 1.04–1.25; 8 workers 1.5–2.0 (more workers than CPUs) |

Forkserver gave similar speed-ups but had higher CPU ratios in Docker (1.3–2.1 at 4 workers).

### Other results
- **Exactness (P2-3):** hashes identical to serial in every batch, workload, worker count and start method, on the Mac and in Docker. The same holds through `TRLReward` (outputs `==`).
- **Memory (P2-4):**
  - Worker private memory (USS) grew 2–12 MB on the Mac and 6–16 MB in Docker.
  - In Linux with the full preload test, each of 4 workers showed RSS 2.8 GB but USS 12 MB and PSS 569 MB: the pages are shared copy-on-write.
  - Transparent huge pages were set to "always" in the Docker VM; no COW growth was seen.
- **Overhead (P2-5):** measured with a no-op scorer over the same chunks. At 4 workers it was 6.0 / 222 ms = 2.7 % (Mac quiet run), 8.9 / 154 = 5.8 % (Mac run 1), and 17–23 / 329–333 ms = 5–7 % (Docker).
- **TRL call pattern (P2-6):** exactly 1 dispatch per generation batch over 10 steps of 256 plus one eval-sized batch of 64. Outputs equal the serial `TRLReward`'s; worker pids are stable. Holds for fork and forkserver, on the Mac and in Docker.
- **Pool start:** fork 0.02–0.9 s; forkserver 2–10 s, because the server preloads.

### Robustness findings (`trl_pattern.py`, `fs_main_probe.py`)
1. **Worker death.** A SIGKILLed worker makes `ProcessPoolExecutor` raise `BrokenProcessPool` within 5–430 ms. `multiprocessing.Pool.map_async` never returns: no result after 30 s, on both Mac and Linux. Use `ProcessPoolExecutor`. Python documents this BrokenProcessPool behaviour.
2. **When workers fork.** With the fork context, `ProcessPoolExecutor` forks *all* workers at the first `submit` (`_safe_to_dynamically_spawn_children = False`) and never respawns. So warm it with n no-op tasks immediately after creation, before the model loads and CUDA starts. Otherwise the fork happens at the first reward call, after CUDA initialisation.
3. **Threads at fork time.**
   - On Linux the process already had 22–23 OS threads after its imports and the preload (OpenBLAS thread pool, pyarrow). The first Docker run showed Python 3.12's DeprecationWarning about forking a multi-threaded process.
   - After torch's intra-op pool had started, the warning appeared on both platforms, and scoring still completed.
   - CPython counts threads through `/proc/self/stat` or macOS `task_threads`.
   - OpenBLAS resets its own thread pool across fork via `pthread_atfork`. libgomp is not fork-safe, but the children never run OpenMP code.
   - Risk accepted with this mitigation: create the pool before CUDA starts, and keep `cv2.setNumThreads(0)` in the worker initializer.
4. **Forkserver re-imports the training script.**
   - Workers re-import the parent's `__main__` as `__mp_main__`.
   - With top-level torch and datasets imports, the workers **crashed on macOS** in 3 of 3 runs (Objective-C `+[NSCharacterSet initialize]` fork-safety abort, then `BrokenProcessPool` even with `OBJC_DISABLE_INITIALIZE_FORK_SAFETY`).
   - On Linux the same setup worked but pool start took 49.7 s against 27.6 s with a light main (emulated).
   - Preloading inside the fork server must also be thread-free.
   - Conclusion: forkserver only suits an import-light entry point. A TRL script is not one.
5. **macOS fork** from the main process worked in every run. The Mac is a development box only; production is Linux.

**Verdict H-2:**
- P2-2 through P2-7 pass. **P2-1 fails** on the Mac: 2.06 / paired 2.91 at load 21; 2.46 / 2.47 at load 7.5→12.
- The CPU ratio of about 1.0–1.05 at 2–4 workers shows the pool adds no inherent CPU work. The missing speed-up is waiting for cores on a machine running at load 12 or more with 4 performance cores, not overhead.
- Pending the Vast measurement, accept it as the Phase 2 design.

## 3. H-K3SYN: K3/K4 harness without RIGA data

`synth.py` generates 500 × 500 frames with 6 raters, smooth radial-harmonic disc and cup outlines, the cup inside the disc. Parameters come from RIGA quantiles only. After calibration:
- disc radius median 108 px (real 106);
- cup/disc ratio 0.47 (real 0.43);
- band width: disc 11.1 px (real 10.1), cup 16.0 px (real 19.2);
- completion length 861 characters (real 856).

`k3syn.py` compares real and synthetic; `k3_harness.py` is the Vast tool.

| Run | RIGA-only wall ratio (95 % CI) | Mixed wall ratio (CI) | p99 real / synthetic (ms) | CPU-time ratio, secondary (CI) | p50 per rollout, real / synthetic (ms) |
|---|---|---|---|---|---|
| Mac run 1, load 28–48 | **1.03** (0.79–1.33) | 1.24 (0.91–1.58) | 5.47 / 6.54 | — | 0.965 / 0.979 |
| Mac replication, load 11–15 | **0.94** (0.83–0.99) | 1.03 (0.82–1.24) | 18.1 / 19.4 | RIGA **0.972 (0.935–1.016)**, mixed **0.990 (0.897–1.119)** | 2.70 / 2.63 |
| Docker, emulated | **1.08** (0.81–1.21) | 0.99 (0.73–1.22) | 50.5 / 33.4 | — | 3.80 / 3.63 |

- PS-1a passes in all three runs.
- PS-1b fails in both Mac runs. PS-2 fails in run 1 and passes in the replication. PS-3 fails in both Mac runs (1.06 ms and 1.35 ms).
- PS-4 passes (1.08).
- PS-5 passes: the digest is `e31d8f1e…` on arm64 and on x86.

Per-batch wall ratios ranged from 0.2 to 5.6 under load. The wall-time CIs measure the machine's load, not the workload.

**Verdict: inconclusive on the pre-registered wall-time criteria; accepted on the load-robust measures** (CPU-time ratio CI inside ±7 %, p50 within 1.5–4.5 %, character counts within 1 %).

`k3_harness.py` ran end-to-end in Docker: all hashes equal, K3/K4 reported (they fail only because of emulation). It now interleaves serial and pool per batch. Run it on the first Vast box (`python -I -m experiments.phase2.lab.B4.k3_harness --out …`). On a datacenter host, also run `--real-riga` once to recheck PS-1 on a quiet machine.

## 4. H-3: the C++ batch scorer

**Decision.**
- The pre-registered rule says "needed" because P2-1 failed on the Mac. P2-2 passed: 51–61 ms at 8 workers, against the 125 ms target and the 250 ms gate.
- Projection to x86:
  - Serial RIGA-only on Vast is probably between M6's quiet-Mac 104 ms and A6's 241 ms (load 9).
  - With the measured 2.5× at 4 workers, the 4-worker wall would be about 40–100 ms, under 125 ms.
- **Recommendation:** do not build H-3 before the first Vast `k3_harness` run.
- **Build it if:**
  - the RIGA-only 4-worker paired wall exceeds 125 ms or the speed-up is under 3× on a quiet box;
  - or Phase 3 needs larger groups or more rollouts per step;
  - or the box has fewer than 4 spare cores beside the dataloader workers.

**Design, if built (cpp/polycore, about 600–900 lines plus differential tests):**
- **API (nanobind).**
  - `CacheTable.add(...)` registers each item and object once, taking the arrays the Python `ExpertCache` has already computed: raters (u8, k×H×W), union, core and vote (u8), `dist_out_union` and `dist_in_core` (f32), rater areas, rater signatures, voters, `min_area`, perimeter, t5 and box. C++ never computes an EDT, so the scipy-exact distances stay authoritative.
  - `score_rewards(texts, truncated, handles[n, n_obj], W, H, threads)` returns `format_ok`, `legal`, `iou_tc`, `iou_min`, `f5`, `g_area` (float64), per-object signatures, an error-code index, defects and a `fallback` mask.
  - Python then computes `topo_freq`, `topo_any` and `sig_id` (CRC32 of the `repr` of the signatures) from the signatures, with `rater_answers` cached per item, and runs `advantages_v4` unchanged.
- **Exactness.**
  - The parser fast path accepts only the canonical `serialize` grammar. Any other text falls back to the Python parser, so error codes, `too_long` and truncation keep Python's semantics.
  - Areas and intersections are integers. IoU is `double(inter) / double(union)`.
  - The band depth is a float64 sum of float32 distances; A6 showed it is exact in any order.
  - F5 and `g_area` use the same operation order as `score.py`. Compile with `-ffp-contract=off` and without `-ffast-math`.
  - Signatures use 8-connected components plus 4-connected background holes on the ROI, including the size-filter slow path, and are differential-tested against `topology.signature`.
- **Threads.** The polycore `score_batch` pattern: an internal `std::thread` pool, GIL released, per-thread workspace, results written by index, and the first worker exception rethrown after join.
- **Acceptance.**
  - 0 mismatched NUMERIC terms against `score_completion` on 10 batches × 3 workloads plus the M6 adversaries and a fuzz set.
  - Identical K10 hashes on arm64 and x86.
  - Serial ≥ 1.5× and 8 threads ≥ 5× faster than Python serial on RIGA-only.

## 5. Recommended reward-serving design for Phase 2

1. At trainer start, **before model or CUDA initialisation**, preload every RL cache of the train and eval splits into a dict. Use 2–4 threads and check digests as `cache_loader` does. Cost: about 1 min and about 4.5 GB on Linux; size the Vast box's RAM for it.
2. Call `gc.freeze()`. Then create `ProcessPoolExecutor(4, mp_context=get_context("fork"), initializer=cv2.setNumThreads(0))` and warm it with 4 no-op tasks so every worker forks now.
3. Use `TRLReward(caches=ref_loader, pool=pool)`. Jobs carry `CacheRef(key)`, and there is one dispatch per generation batch in chunks of G. `TRLReward` needs no change; the library needs a small `objdel.train` serve helper (about 60 lines) and a preload or unbounded mode for `cache_loader`.
4. A dead worker raises `BrokenProcessPool`: the run fails loudly and resumes from a checkpoint. Never use `multiprocessing.Pool`, and never forkserver or spawn with a heavy training script as `__main__`.
5. Size the pool at 4 workers, leaving cores for the 4 or more dataloader workers. Measure on Vast with `k3_harness`, which gives the K3/K4 gate on LIDC plus synthetic RIGA with paired pool speed-ups.

## 6. Sources

- Python multiprocessing start methods: macOS default `spawn` since 3.8 because system libraries may start threads; POSIX default `forkserver` since 3.14; fork DeprecationWarning since 3.12; `set_forkserver_preload` ignores ImportError. https://docs.python.org/3/library/multiprocessing.html
- `ProcessPoolExecutor`: BrokenProcessPool on abrupt worker exit; fork no longer the default from 3.14. https://docs.python.org/3/library/concurrent.futures.html
- CPython 3.12 `process.py` (fork launches all workers at first submit, no dynamic respawn): https://github.com/python/cpython/blob/3.12/Lib/concurrent/futures/process.py ; `forkserver.py` preload loop: https://github.com/python/cpython/blob/3.12/Lib/multiprocessing/forkserver.py ; fork warning thread count (`/proc/self/stat`, `task_threads`): https://github.com/python/cpython/blob/3.12/Modules/posixmodule.c
- `gc.freeze` (call `gc.disable` early, `gc.freeze` before fork, `gc.enable` in the child): https://docs.python.org/3/library/gc.html ; Instagram's copy-on-write experience: https://instagram-engineering.com/dismissing-python-garbage-collection-at-instagram-4dca40b29172 , https://medium.com/instagram-engineering/copy-on-write-friendly-python-garbage-collection-ad6ed5233ddf
- PyTorch multiprocessing (CUDA needs spawn or forkserver in children; fork with lock-holding threads can deadlock; set threads per worker): https://docs.pytorch.org/docs/2.14/notes/multiprocessing.html ; "Cannot re-initialize CUDA in forked subprocess": https://discuss.pytorch.org/t/cannot-re-initialize-cuda-in-forked-subprocess-on-network-to-device-operation/138090
- OpenBLAS `pthread_atfork` thread-pool reset: https://testgitea.trustie.net/floraachy/OpenBLAS/commit/138a8413903776242cb3c184df9838ba54d04044 ; 2025 post-fork deadlock regression and fix: https://github.com/OpenMathLib/OpenBLAS/issues/5520
- numba threading layers and fork safety (omp/libgomp not fork-safe): https://numba.readthedocs.io/en/stable/user/threading-layer.html
- macOS Objective-C fork-safety abort: https://wefearchange.org/2018/11/forkmacos.rst.html , https://blog.phusion.nl/why-ruby-app-servers-break-on-macos-high-sierra-and-what-can-be-done-about-it/
- glibc `M_MMAP_THRESHOLD` (128 KiB, dynamic up to 32 MiB): https://man7.org/linux/man-pages/man3/mallopt.3.html ; numpy huge-page madvise: https://numpy.org/doc/2.0/reference/global_state.html
- fork latency about 4 ms per GB of mapped memory: https://lkml.iu.edu/hypermail/linux/kernel/2209.3/04474.html
- TRL GRPO (rewards per generation batch on each process; async only for I/O): https://huggingface.co/docs/trl/main/en/grpo_trainer
- verl RewardLoopManager (`reward.num_workers` chunks): https://verl.readthedocs.io/en/latest/advance/reward_loop.html
- Docker Rosetta emulation: https://docs.docker.com/desktop/settings-and-maintenance/settings/

(Report text returned by the B4 agent and saved by the coordinator; the harness blocked the agent's write.)
