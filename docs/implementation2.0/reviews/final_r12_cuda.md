# Final review round 12: CUDA and GPU performance (Implementation.md)

Reviewer role: CUDA / GPU performance (kernels, GPU placement, syncs, memory, Blackwell sm_120, rental plan).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Review only. No catalog install (installing needs the owner's approval). It sent me to NVIDIA's own docs for the sync facts: I read the CUDA runtime "API synchronization behavior" page for the pageable-copy rule that D7b and D5c rely on (no item: D7b's one wait per batch stays correct). |
| engineering:code-review | Line-by-line check of the D2 port against K1's `Engine`/`scan` (k1poly.cpp:119-165, 780-907) and gpucheck's `Lane`/`Block` (gpucheck.cpp:131-160). This found item 1. |
| engineering:debug | Traced what the item 1 aliasing does on each build: `Block::get` returns the same base when the capacity is enough (gpucheck.cpp:139), and the CPU `scan` never touches `tmp` (k1poly.cpp:144). So the bug shows only on the GPU. Also traced what `time/rewards_s` measures (item 2). |
| engineering:system-design | Walked the reward path (prefetch thread, engine lock, lane 0, Drain) and the D5c pool path (lane-ordered `cudaFreeAsync`, tables freed only under `e.mu` after Drain). No defect apart from item 1's missing third arena. |
| engineering:architecture | Checked D5c's private pool (no `cudaDeviceSetMemPool`, release threshold UINT64_MAX) and the "two labellers" decision against the code. No defect. |
| engineering:documentation | Each item quotes the source file and heading, then the code line or measurement, then one exact fix. |
| huggingface-skills:huggingface-best | No model choice in my area this round. Nothing changed. |
| huggingface-skills:huggingface-papers | D.2 cites kernel and CCL papers and NVIDIA docs only. No arXiv claim in my area needed a new check. Nothing changed. |

Web checks (2026-11-09): docs.nvidia.com/cuda/cuda-runtime-api/api-sync-behavior.html. For a synchronous pageable
H2D copy, "a stream sync is performed before the copy is initiated". The async form "might be synchronous with
respect to host". D7b's plan (pageable staging, one wait per batch) fits this. No item.

Re-verified and found correct (no item):
- D-G8 baseline: 4 `wait()` calls (gpucheck.cpp:1802, 1574, 1626, 1639) plus 2 `Drain` syncs (177-184).
- D5 step 1: crop rows per object = `bh + 2` (gpucheck.cpp:1583). The skip is exact, because the relabel runs whenever either filter changed something.
- B3.1 column-sum identity: 32w·popc + 1/2/4/8/16 masked popcounts, LSB = leftmost pixel (gpucheck.cpp:1144).
- `ObjOut` = 48 B (gpucheck.cpp:1199-1202). D5c pinned bound: 20.5 MB per name. `plan()`'s `in_bytes` is the xy array plus small per-ring and per-shape arrays (gpucheck.cpp:1079-1082), so the bound holds.
- K2 `warp_kernel`: warp-uniform exit before `__ballot_sync`/`__shfl_xor_sync` (k2score.inc:46-49).
- sm_120 full occupancy: at most 40 registers at 128 threads (64 K / 1,536).
- C1 byte counts: 25.2 MB, 201 MB and 1.61 GB.
- V1 cap 3.35 h and `BudgetStop` $2.05.
- `vllm_gpu_memory_utilization` 0.35 (configs.py:192), so D-G17 (c) ≤ 30 GB is reachable.
- Round 11 CUDA items 1-7 (all minor) are still open in the text, for example D-G13 still uses 0.641. They are not repeated here.

Counts: blocking 0, major 1, minor 4.

---

## 1. MAJOR. The D2 port maps K1's three device blocks onto a Lane that has two; an alias for CUB's scan scratch passes every Mac test and corrupts the GPU result

Sources: D (heading "D2. Work-optimal polygonizer inside the engine (port K1-draft)", step 1; heading "D1", engine
counters list, D5 step 5 "D5c").

Doc:
- D.md, D2 step 1: "Use gpucheck's `Lane` (lowest-priority non-blocking stream [code: gpucheck.cpp:66-71]) and its
  arenas in place of K1-draft's `Engine`."
- D.md, D1 counters: "the grow-only `Block<true>` arenas `dev` and `score` [code: gpucheck.cpp:131-154]". D5c lists
  the same two arenas. Neither the doc nor K1 DESIGN.md:230 ("Reuse gpucheck's `Lane` ... instead of K1's own
  `Engine`") says where K1's third block goes.

Code:
- K1 needs three device blocks: `Block<true> devA, devB, tmp;` (k1poly.cpp:162). `tmp` is CUB's scan scratch:
  `void* t = tmp.get(bytes > 0 ? bytes : 1); cub::DeviceScan::ExclusiveSum(t, bytes, ...)` (k1poly.cpp:140-142).
  `chunk` calls `scan(..., E.tmp)` six times while layouts A and B are live (k1poly.cpp:824-825, 871-872, 876, 884).
- A gpucheck `Lane` has two: `Block<true> dev, score; Block<false> in, out;` (gpucheck.cpp:155-159).
- `Block::get` returns the same base pointer when the capacity is enough: `if (bytes <= cap) return p;`
  (gpucheck.cpp:139).
- The CPU build of `scan` ignores the scratch block: `(void)st, (void)tmp;` (k1poly.cpp:144).

Failure: say an engineer passes `L.dev` (layout A) or `L.score` (layout B) as `tmp`. On the GPU, CUB then writes its
temporary storage over the start of the live layout: A's input bytes, starts and sizes, or B's run arrays. The
rings come out wrong. On the CPU build nothing changes, so every planned $0 check stays green: the D2 fuzz and split
tests, `make test-order`, the D3 hash rebuild (D-G4) and `test_rings.py`. The bug appears first in R-D2 (V3 item 1,
paid), as a D-G3 mismatch, after the build and pytest minutes are spent. If R-D2 is not approved, the bug stays in
the shipped engine. The caught-before-spend rule does not apply, because no $0 check can see it.

Fix (D.md D2 step 1, D1 counters, D5c):
- D2 step 1: "Add `Block<true> scan;` to `Lane` (K1's `tmp`, CUB scratch only). Map K1 to the lane like this:
  `devA` -> `L.dev`, `devB` -> `L.score`, `tmp` -> `L.scan`, `in` -> `L.in`, `out` -> `L.out`. Never pass a block
  that holds live data to `scan`."
- D2 tests: add a test-only CPU define `GC_POISON_SCAN` (beside `K1_ORDER`). Under it, the CPU `scan` calls
  `tmp.get(4 * (n + 1))` and fills that span with 0xFF before the loop. `make test-order` then also runs
  `tests/test_cuda_checker.py -k masks` under `GC_POISON_SCAN`. If `scan` is given A or B, the run fails on the Mac.
- D1 counters and D5c: "the grow-only `Block<true>` arenas `dev`, `score` and `scan`". `scan` uses the same pool
  path (`cudaMallocFromPoolAsync` / `cudaFreeAsync` on the lane stream).

## 2. MINOR. D.1 item 9 calls `time/rewards_s` the scoring time, but it is only the wait for the prefetched result, and 23.4 s is the other arm's step

Sources: D (heading "D.1. Goal and why", item 9).

Doc: "Reward scoring takes a mean 0.0046 s (median 0.0045, range 0.0039-0.0054 s, 15 steps) of a 23.4 s step
[measured: REPO/runs/r5/base/rl.json `time/rewards_s`]."

Code and measurement:
- `time/rewards_s` times `_calculate_rewards` (fast.py:638-640). Scoring has already started on a background thread
  in `prefetch_rewards` (fast.py:626; rewards.py:111-119 `self._thread.submit(score_batch, jobs)`). So the value is
  the time left to wait, not the scoring time. The true scoring time is not measured.
- `runs/r5/base/rl.json`: 15 values, mean 0.004603 s, median 0.004509 s (these match). But
  `timing_s.per_step` = 31.05 s. The 23.4 s step is the shared-prefix arm (`runs/r5/shared_rl.json`).

The conclusion still holds, because the scoring overlaps the forwards. D-G8 records `time/rewards_s` "for
information only", so no gate changes.

Fix: "The reward wait after the prefetch (`time/rewards_s`, the part not hidden behind the forwards) is a mean
0.0046 s of a 31.05 s base-arm step [measured: runs/r5/base/rl.json]. The scoring itself runs on the prefetch thread
(fast.py:626; rewards.py:119), and its time is not measured; D1's replay bench measures it."

## 3. MINOR. Kernel evidence paths cite files that have moved; D9 and DOC tell the engineer to move files that are already moved

Sources: D ("D.0" text, "D.1" items 3 and 7, D6 step 1, D9 "kernel drafts"), DOC (kernel folders bullet).

Doc:
- D.md D.1 item 3: "nvcc 13.0.88 builds it for sm_120 [measured: K1/results/nvcc_docker_sm120.log]".
- D.md D.1 item 7 and D6 step 1: "[measured: ... K2/build/run_cpu_mac.txt]".
- D.md D9: "copy the K2 text logs from `docs/implementation2.0/kernels/K2/build/` to `K2/results/` (D1 step 4)". D.md
  line 9: "sit in `K2/build/` today".
- DOC.md: "Move the text evidence ... `K2/build/run_cpu_mac.txt`, ... `K1/results/nvcc_docker_sm120.log` ...".

Measurement (`ls`, 2026-11-09):
- `K2/build/` does not exist. `K2/results/` holds `run_cpu_mac.txt`, `nvcc.log.txt`, `ptxas.txt`, `mut_run.sh` and
  `mut_run2.py`.
- `K3/results/` holds the three `.log.txt` files.
- The K1 file is `K1/results/nvcc_docker_sm120.log.txt`.
- `K1/DESIGN.md:180` still cites `results/nvcc_docker_sm120.log`, with no `.txt`.

Fix: in D.md, cite `K1/results/nvcc_docker_sm120.log.txt` and `K2/results/run_cpu_mac.txt`. Delete the D9
"copy the K2 text logs" bullet and the line-9 sentence. In DOC.md, replace the move list with "done (2026-11-09):
check only that `K*/DESIGN.md` cite `results/` paths that exist". Fix the one K1 DESIGN.md citation (`.log` ->
`.log.txt`).

## 4. MINOR. Front section 5 and D1 step 0 describe runbook defects that are gone, and give line numbers that have moved

Sources: 00_front (heading "5. Custom kernels index"), D (heading "D1", step 0).

Doc:
- 00_front 5: "Today the three runbooks still name the session scratchpad [measured: grep, 2026-11-08: K1
  GPU_RUN.md:15-18, K2 GPU_RUN.md:33, K3 GPU_RUN.md:15]."
- D.md D1 step 0: "Delete every scp and rsync line (today K1 GPU_RUN.md:51, 126; K2 GPU_RUN.md:33, 117; K3
  GPU_RUN.md:16, 93-94)".

Measurement (grep, 2026-11-09):
- No scratchpad path is left in `docs/implementation2.0/kernels/` (outside `build/`).
- `K1/GPU_RUN.md:18-19` and `K2:32` already use `REPO=$(git rev-parse --show-toplevel)`.
- The scp and rsync lines left are K1:50, 125; K2:34, 118; K3:17, 94, 95.
- `torch 2.14.1` is still in K3:7.

Fix:
- 00_front 5: "The scratchpad paths are gone. The runbooks still hold scp/rsync lines and K3's `torch 2.14.1`
  [measured: grep, 2026-11-09]."
- D1 step 0: cite K1:50, 125; K2:34, 118; K3:17, 94-95, or drop the line numbers and keep the S09 and preflight grep
  as the check.
