# Final review round 15: CUDA and GPU performance (Implementation.md)

Reviewer role: CUDA and GPU performance. I checked kernels, GPU placement, syncs, memory and the rental plan.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Review only. I installed nothing, because an install needs the owner's approval. No catalog skill covers a custom engine's stream-ordered pool or pinned staging. It sent me to the CUDA runtime memory-pool page (web checks below). |
| engineering:code-review | I re-checked the D.4 counter text, D5 step 1 and D5c against `gpucheck.cpp`: wait sites :1574, :1626, :1639, :1802; `Block::get` :135-150; `alloc`/`release` :78-84; `add_tables` :1720-1763; `drop_tables` :1765-1776; `plan` :1070-1090 (`c.words - c.defects` is the defects array, so D5c's `out` bound is right). I also checked the bindings' GIL handling in `gpucheck_py.cpp`. This gave item 4. |
| engineering:debug | I traced the eval process after D7 step 6: `metrics()` still calls `score_many`, so `score_polygons` and `eval_sets` share one engine in `final_eval`. This gave item 1. I traced `ExpertCache.load` -> `__init__` -> `signature` across D4b and D7b. That path is safe: no paid run sits between the D4b and D7b merges. |
| engineering:system-design | I walked the GPU memory of the post-E4 eval process (in-process vLLM plus the engine's context, tables and scratch). This gave item 2. I also re-checked V1 step 0 (`run_box.sh --build-only`). The `vllm_native` image has no `/opt/requirements.txt`, so the 13.4 -> 13.0 compiler swap has no constraint conflict. The K1 log shows the swap, and `rl_full.yaml` trained after it. |
| engineering:architecture | I checked the D5c pool decision again with the CUDA reuse defaults (all three reuse policies are on by default). Cross-stream reuse is therefore safe without extra events. No item. |
| engineering:documentation | Each item gives the source file, the heading, the doc quote, the code or measurement, and one exact fix. |
| huggingface-skills:huggingface-best | No model choice is in my area this round. Nothing changed. |
| huggingface-skills:huggingface-papers | No new paper claim is in my area. Nothing changed. |

Web checks (2026-11-09):
- docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY__POOLS.html: `cudaMemPoolReuseFollowEventDependencies`,
  `cudaMemPoolReuseAllowOpportunistic` and `cudaMemPoolReuseAllowInternalDependencies` are "(default enabled)".
  `cudaMemPoolAttrReleaseThreshold` is "(default 0)". So D5c's explicit `UINT64_MAX` is needed, and the reuse of
  freed blocks across the two lanes is safe.
- vLLM engine args (vllm.readthedocs.io "Engine Arguments"): `gpu_memory_utilization` defaults to 0.9. Used in item 2.

Re-verified and found correct (no item):
- D-G8 baseline: today `score_polygons` makes 4 waits (:1802, :1574, :1626, :1639) and 2 `Drain` syncs (:177-184).
- D5 step 1 line references (:1614 the second label pair, :1615 `zero`, :1616 `Holes`, :1583 bh + 2 crop rows).
- D5c `out` bound: the `score_polygons` request is `L.out.get(c.words - c.defects)` (:1799). In `plan` the defects
  array sits directly before `words` (:1086), so the request is n x 13 x 8 B rounded up to 256 B, as the doc says.
- K2 `for_each_warp`: the exit is warp-uniform (`i = tid >> 5`), so `__ballot_sync(0xFFFFFFFF, ...)` is legal
  (K2/k2score.inc:44-47, 66-68).
- C5 parameter counts: q/k/v/o/gate/up/down = 50.3 M per layer x 28 = 1.41 B. The adapter is 622,592 x 28 = 17,432,576.
- C6 step 2: `max_num_seqs = 16 x 1 x 16 = 256` [code: .venv trl/trainer/grpo_trainer.py:1186-1188], and
  `max_num_batched_tokens=4096` [code: .venv trl/generation/vllm_generation.py:368].
- V1 step 0 runs as one `bash -c` step (B.md B1.5 step 0). So `cmd_run` does not prefix `/opt/venv/bin/python -I`
  to the second command (policy.py:484-487). `tests/test_gpu_scoring.py` passes under `-I` (12 passed, Mac, 2026-11-09).
- The sanitizer chain works on the box with 13.4.92: memcheck and initcheck show "ERROR SUMMARY: 0 errors"
  [measured: runs/k1/cuda_checker/memcheck2.txt].

Not re-raised (already in earlier logs or reviews): stale kernel evidence paths (r12 item 3), the runbook
line numbers (r12 item 4), the pool-size-1 wording (r13 item 4, r14 item 4), the `Stats` atomics (r14 item 2), the
pageable-copy sync claim (r14 item 3), D-G17 (a)'s two rules (r14 item 5), and the "names per answer" path to the
engine (r14 item 1).

Counts: blocking 0, major 1, minor 3.

---

## 1. MAJOR. D-G14's eval rule "`arena_bytes` <= `eval_sets` scratch cap + table bytes" cannot be evaluated. `final_eval` never records table bytes, and its `arena_bytes` also holds the `score_polygons` arenas, which the rule does not count

Sources: D (heading "D.6. Acceptance gates", row D-G14; heading "D7", step 5; heading "D.9.6", rows S12 and S13;
heading "D.9.7", row B04; heading "D.4. Tasks", "Engine counters").

Doc:
- D.6 D-G14: "it records the `eval_sets` time, `arena_bytes` and peak memory in `final_eval.json["engine"]` (D7
  step 5). ... `arena_bytes` <= the `eval_sets` scratch cap + table bytes. No OOM."
- D7 step 5: `final_eval` "records, per view, the `eval_sets` wall time, `gc.stats()["arena_bytes"]` and
  `torch.cuda.max_memory_reserved()` into `final_eval.json["engine"]`."
- S12: "`final_eval.json["engine"]` holds `eval_sets_s`, `arena_bytes`, `max_memory_reserved`".
- S13: `test_b04` reads "`final_eval.json["engine"]` (the B3 val eval, B3.9; `arena_bytes` <= cap + table bytes)".
- D.4 defines only "the current arena bytes, and `live_device_bytes` (device bytes allocated and not yet released:
  arenas, tables, pool blocks)". It does not say which blocks "arena bytes" sums. D7 does not say which arena
  holds the `eval_sets` scratch ("It runs on the scorer lane").

Code:
- After D7, `metrics()` still gets its training-path terms from `score_many` (D7 keeps that call; only `to_masks`
  goes) [code: src/objdel/eval/metrics.py:250-258]. `score_many` calls `_gpucheck.score_polygons`
  [code: src/objdel/reward/score.py:115].
- `score_polygons` sizes the checker arena `L.dev.get(c.size)` [code: cuda/gpucheck.cpp:1799], and `score_words`
  sizes the scorer arena `L.score.get(...)` twice [code: cuda/gpucheck.cpp:1565, 1592]. These grow-only arenas are
  never freed [code: gpucheck.cpp:131-150].

So, in the one eval process, `arena_bytes` = the `score_polygons` arenas + the `eval_sets` scratch (wherever it
lives). The right side holds "table bytes", a number that no record holds. `live_device_bytes` is not written
either. The engineer who writes `test_b04` must guess both sides. The guess decides pass or fail:
- If "arena" counts `dev`/`score`/`scan` and the scratch reaches its chunk cap, the rule fails on correct code
  whenever the `score_polygons` arenas exceed the table bytes.
- If the engineer reads "table bytes" as 0, the rule fails as soon as any scorer arena exists.

No $0 test exposes this. S12 checks only that the three keys are written, with fakes. The rule is read only on
the Mac after the paid V3 `pull.sh`. A wrong verdict there is a wrong gate result. This is an untestable gate.

Fix (D.md D.4, D7 step 5, D-G14, S12, S13, B04):
- D.4 "Engine counters": add "`arena_bytes` = the sum of `cap` over every `Block<true>` (`dev`, `score`, `scan`) of
  both lanes. `table_bytes` = the sum of the sizes of the live table blocks in the registry. `eval_scratch_peak` =
  the largest scratch that one `eval_sets` chunk asked for since `reset_stats()`. The CPU build counts the same."
- D7 Files: "`eval_sets` takes its scratch from `L.score` of lane 0 and sets `eval_scratch_peak`."
- D7 step 5 and S12: `final_eval.json["engine"]` holds `eval_sets_s`, `arena_bytes`, `table_bytes`,
  `eval_scratch_peak`, `live_device_bytes` and `max_memory_reserved`.
- D-G14 eval half and S13/B04: replace "`arena_bytes` <= the `eval_sets` scratch cap + table bytes" with
  "`eval_scratch_peak` <= `max_scratch_bytes` (the D7 cap), and `live_device_bytes` <= `arena_bytes` +
  `table_bytes` (no hidden allocation). No OOM. `arena_bytes` and `table_bytes` are recorded, not gated."

## 2. MINOR. After E4, `final_eval` builds `VllmSampler` in process, but no value of `gpu_memory_utilization` or `max_model_len` is given. D7 step 6 then runs the engine in the same process after generation

Sources: E (heading "E4", step 1 class sketch and step 3), D (heading "D7", step 6), 00_front (heading "6.3", "V4").

Doc:
- E4 step 1: `class VllmSampler:  # vllm.LLM offline, in process` with
  `__init__(self, ckpt, model, gpu_memory_utilization: float, max_model_len: int)`. The step also carries over
  `check_memory`.
- E4 step 3: "`final_eval` calls `VllmSampler` + `evaluate` ... one `VllmSampler` build per process". No value is given.
- D7 step 6: "`metrics()` device work runs in ONE process: the main process after vLLM generation (`final_eval`, ...)".
- Front 6.3 V4: "the sampler-arm `final_eval` (pool size 1, as A3)".

Code and web:
- Today's sampler defaults to 0.85 and requires 90 % free memory before start [code:
  experiments/phase2/sample_vllm.py:31, 79, 87]. vLLM's own default is 0.9 [web: vLLM "Engine Arguments"].
- After D7 step 6 no pool is left (`final_eval.py:69` loses its `fork_pool` call). So "pool size 1" in V4 is stale.

At 0.9 on 32,607 MiB, vLLM keeps about 29.3 GB. The eval process still needs its own CUDA context, the tables and
up to 256 MiB of `eval_sets` scratch, and the room left is small. Why minor: V4 is the only post-E4 eval in the plan
(V3 runs before E, on the trainer path at 0.35). V4's E-G13 would also show an OOM at once.

Fix: in E4 step 3, add "`final_eval` passes `gpu_memory_utilization=0.80` and `max_model_len=2048` (as
`sample_vllm.py`'s `--max-model-len` default), so at least 6 GB stay free for the scoring that D7 step 6 runs in
the same process". In front 6.3 V4, replace "(pool size 1, as A3)" with "(scoring in process, D7 step 6)".

## 3. MINOR. D8 `--gpu-check` does not say whether a sanitizer failure (step 4) stops steps 5-6. Under `set -euo pipefail` it does, so a sanitizer timeout also loses the paid D-G7/D-G8 benches and A's G02b/G03

Sources: D (heading "D8", the `--gpu-check` steps 3-6 and the `d_gpu_test.yaml` bullet).

Doc:
- D8: "Each tool runs as `timeout 300 "$CUDA_HOME/bin/compute-sanitizer" --tool T --error-exitcode 99 ...`. It
  writes each tool's exit code and wall time to `"$OUT/sanitize.json"`".
- D8 steps 4-6: "4. `sanitize` (above); 5. the old-engine build for D-G7 ...; 6. the benches."
- `d_gpu_test.yaml` step 2 (A's G02b and G03) follows step 1 in one `&&` chain.

Code: `run_box.sh` runs under `set -euo pipefail` [code: experiments/cuda_checker/run_box.sh:15]. `cmd_run` joins
plan steps with `&&` [code: experiments/phase2/vast/policy.py:490].

So a racecheck run that only hits its 300 s cap ends the mode before the old-engine build and the benches. Step 2
of `d_gpu_test.yaml` then does not run either, although the box is paid for. Why minor: front 6.4 rule 4 ("A gate
fails: stop that item") makes stopping a defensible reading, and no measurement shows a timeout.

Fix (D8, `sanitize` bullet): "`sanitize` records each tool's exit code in `sanitize.json` and returns 0. The mode
runs steps 5-6, then exits 1 if any tool's code was non-zero or 124 (timeout), and names the tool. D-G5 reads
`sanitize.json`."

## 4. MINOR. `drop_tables` (and `resident`) keep the GIL while they wait for the engine lock, but the scoring entries release the GIL while they hold it. A table drop on the trainer thread then stalls the whole interpreter for one scoring call

Sources: D (heading "D5", step 5 D5c: "`drop_tables` frees with `cudaFreeAsync` on the lane stream (no device-wide
sync)"; heading "D7b", `ExpertCache` bullet).

Code:
- `m.def("drop_tables", &gpucheck::drop_tables, nb::arg("handle"));` has no GIL release
  [code: cuda/gpucheck_py.cpp:193]. The scoring and upload entries release it
  [code: cuda/gpucheck_py.cpp:127, 147, 186, 204, 220].
- `drop_tables` takes `e.mu` [code: cuda/gpucheck.cpp:1767-1768].
- It runs from `weakref.finalize` [code: src/objdel/reward/cache.py:256], so on whichever thread the Python
  garbage collector runs.

While the reward thread is inside `score_polygons` (GIL released, `e.mu` held), a collection on the trainer thread
blocks in `drop_tables` and keeps the GIL. Every Python thread then waits for the whole scoring call. No deadlock
follows, because the scoring thread drops `e.mu` before it takes the GIL back. Why minor: in RL the caches stay
alive in the dataset, so finalizers rarely run.

Fix (D5c, last sentence on `drop_tables`): "The `drop_tables`, `resident`, `stats` and `trim` bindings release the
GIL (`nb::call_guard<nb::gil_scoped_release>()`) before they take `e.mu`."
