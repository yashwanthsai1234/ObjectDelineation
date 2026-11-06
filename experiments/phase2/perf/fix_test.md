# Fixes and their tests (2026-11-06)

The fix for each issue in `issues.md`, the test that proves it, and how we check it on the GPU. Target: R2 (the RL
G-probe on R1's e2) on an **RTX 5090**.

Rules for every fix:
- the code stays simple;
- old code that a fix replaces is deleted;
- `make test` and `make lint` pass;
- nothing is rented before `make preflight` passes and the owner sees the hours and dollars.

```
Mac (no cost)                                             Box (owner approves the cost)
F0 decisions ─► F1 ─► F2 ─► F3 ─► F4 ─► F5 ─► F6 ─► F7 script ═► F7 profiling run ─► R2
                RL gen  RL upd  syncs  checker  scoring  5090        (~1 h, quote first)
```

---

## F0. Owner decisions before the code

| # | Question | Recommendation | Why |
|---|---|---|---|
| D-a | Which sampler generates the RL answers: vLLM colocate (F1-A) or HF `generate` in chunks (F1-B)? Plan D13 says HF for the G-probe. | **F1-A, vLLM colocate** | HF: ~1 s per answer, ~2.3 h for the G-probe, and out of memory at 256 sequences. vLLM: ~0.04–0.09 s per answer, paged KV cache. R0 v3 parity passed on the token median for Qwen3-VL-2B, and R1's G-hf-check passed. TRL corrects the small vLLM-to-trainer difference (`vllm_importance_sampling_correction`, default on). |
| D-b | Card for R2 | RTX 5090 (owner, 2026-11-06) | 32 GB is more room for RL. The speed gain shows only after F2–F5. |

---

## F1. RL generation (issues RL-1, RL-2, EV-6)

### F1-A (recommended): vLLM colocate

Change:
- `experiments/phase2/g_probe.py:82`: delete the `"use_vllm": False` override. `grpo_config` already sets
  `use_vllm: True, vllm_mode: "colocate"`.
- `src/objdel/train/configs.py` `grpo_config`: add explicit values:
  - `vllm_gpu_memory_utilization` (0.35 on 32 GB: ~11 GB for weights and KV);
  - `vllm_enable_sleep_mode: True` (vLLM frees its memory during the update);
  - `vllm_max_model_length` = prompt + `MAX_COMPLETION`.
- The probe sets (`g_probe.score` → `sample.sample` → `trainer._generate`) then use vLLM through the same trainer.
  Nothing else changes.
- PLAN D13: add a revision line for this owner decision.

```
before:  HF generate 256 seq ──► KV 20–51 GiB ──► out of memory / ~250 s per step
after :  vLLM colocate 256 seq ─► paged KV in ~11 GB ─► ~10–25 s per step (estimate)
```

Tests (`tests/test_p2_perf.py`, new):
- `test_g_probe_generates_with_vllm_colocate`: build the G-probe config as `probe()` builds it. Assert:
  - `use_vllm` is True and `vllm_mode` is "colocate";
  - sleep mode is on;
  - `0 < vllm_gpu_memory_utilization <= 0.4`.
- `test_vllm_kv_budget_fits`: KV for 256 × (prompt + 1,536) tokens at 112 KB is larger than the vLLM budget. vLLM
  pages it (blocks, preemption), so the test asserts that the budget holds at least
  `G × (prompt + MAX_COMPLETION)` tokens. That is one whole group running away.
- The existing G-probe tests stay green (`test_p2_*g_probe*`).

Box check (F7): one G-probe arm with 3 steps. It must report:
- no out-of-memory error;
- vLLM start time on sm_120;
- the generate time per step;
- `vllm_importance_sampling` metrics in the log.

### F1-B (fallback, if D-a = HF): chunk the generation

Change:
- `grpo_config` for Qwen3-VL-2B with `use_vllm False`: `steps_per_generation = 2` with micro-batch 32 (F2), so each
  `generate` call has 64 sequences. This is the PLAN D13 limit: 13 GiB of KV at 1,856 tokens.
- Optional: `use_transformers_paged: True`, after a smoke test (paged attention in HF, not verified for VLMs).

Test: `test_hf_generation_chunk_fits`: the sequences per generate call × (prompt + 1,536) × 112 KB must be ≤ 13 GiB.

---

## F2. RL update micro-batch (issue RL-3)

Change: in `grpo_config`, use `per_device_train_batch_size = 4 × G = 32` answers and `gradient_accumulation_steps = 8`.
One step still has 256 answers and 32 prompts, as `prereg.G_PROBE` requires. Gradient checkpointing is kept on until F7
measures the peak memory, and then it is set from F7's result.

```
before: 32 micro-batches × 8 answers     (32× the CPU launch work)
after :  8 micro-batches × 32 answers    ( 8× the CPU launch work)
```

Tests:
- `test_grpo_step_keeps_256_answers`:
  - `per_device_train_batch_size × gradient_accumulation_steps == 256`;
  - the batch is a multiple of G;
  - TRL's generation batch is 256 (vLLM) or 64 per call (F1-B).
- The existing `grpo_config` group checks stay (`num_generations == REWARD_V4.group`).

Box check (F7): peak memory and update time per step for micro-batch 8, 16 and 32, with checkpointing on and off.

---

## F3. Fewer syncs (issues RL-4, SFT-3)

Change:
- `sft_config` and `grpo_config`: set `logging_nan_inf_filter: False`. This removes the isnan/isinf sync on every
  micro-batch. A NaN loss is still logged at log steps, and the G-vision check and the loss curve still show
  divergence.
- `sft_config` and `grpo_config`: set `accelerator_config: {"non_blocking": True}` (pinned memory is already on).
- The Trainer subclasses (`train_sft.Trainer`; the G-probe trainer): `floating_point_ops` returns 0. This removes the
  walk over ~2.5k PEFT modules on every micro-batch.

Tests:
- `test_configs_turn_off_the_per_step_syncs`: both configs have `logging_nan_inf_filter is False` and
  `accelerator_config.non_blocking is True`.
- `test_trainer_skips_flop_counting`: `Trainer.floating_point_ops(inputs) == 0`.
- The CPU smoke runs (`test_p2_v1_train.py`, slow) stay green.

Box check (F7): the `syncs` counts in `profile.json` (`aten::_local_scalar_dense`, `cudaStreamSynchronize`) fall
against the base arm.

---

## F4. CUDA checker memory (issues CK-1, CK-2, CK-3)

Change (`cuda/gpucheck.cpp`, `cuda/gpucheck_py.cpp`, `src/objdel/geometry/shape.py`):

1. **One arena** of device buffers, made at the first call and kept. It grows only when a batch is larger. Use
   `cudaMalloc`, not `cudaMallocManaged`. Delete the 20 per-call `Buf`s and their `cudaFree`.
2. **Pinned staging buffers** (`cudaMallocHost`), kept with the arena. Then one `cudaMemcpyAsync` in, the kernels, one
   `cudaMemcpyAsync` out, and **one** `cudaStreamSynchronize` per call.
3. **Return the even-odd mask only, bit-packed** (8 pixels per byte, row-major). Delete the nonzero raster and the
   useless `cudaMemset` of the masks. `shape.to_masks` unpacks with `np.unpackbits(..., bitorder="little")` into the
   same boolean arrays as today. Delete the extra host copy (`assign`).

```
before: 20× (alloc + memset + sync) │ page faults │ 2 rasters × 8 bits │ 20× free+sync   ≈ 8 ms + 0.6 GB/s
after : copy in │ kernels │ copy out (1 raster × 1 bit) │ 1 sync                         ≈ 50 µs + PCIe speed
```

Tests (Mac):
- `test_to_masks_unpacks_bit_packed_rows`: a fake `_gpucheck` returns known packed bytes. Assert that `to_masks`
  gives exactly the boolean mask, for widths that are not a multiple of 8 (63, 64, 65 px).
- `test_to_masks_matches_polycore_on_fixtures`: the existing bit-exact fixtures stay green on the CPU reference.
- Build check: `bash cuda/docker_check.sh`. It compiles for archs 80, 86, 89, 90 and 120 (sm_120 = RTX 5090) in
  Docker, so no GPU is needed.

Box check (F7):
- `experiments/cuda_checker/run_box.sh` on the 5090: correctness eval with **0 mismatches** (as R0 v3), then the timing
  sweep against the polycore loop;
- repeated calls with growing and shrinking batches give identical results (the arena reuse);
- `nsys` with `--cuda-um-cpu-page-faults=true` (if nsys runs in the container): **no** unified-memory page faults.

Done when: the fixed cost per call is ≤ 0.2 ms and the GPU is faster than the polycore loop at batch ≥ 64.

---

## F5. Score one RL step in one checker call (issues RL-5, CK-5)

Change (`src/objdel/train/rewards.py`, `src/objdel/reward/score.py`):
- The reward function parses all 256 answers of the step.
- It calls the checker **once**, in the main process.
- Then it sends the masks (not the polygons) to the CPU pool for `score_mask`.
- The pool workers no longer touch CUDA: there are no extra GPU contexts.

Delete the per-group checker call on the pool.

```
before: 32 groups ─► 4 processes ─► 32 checker calls (4 GPU contexts)
after : 256 answers ─► 1 checker call (main process) ─► masks ─► CPU pool scores
```

Tests:
- `test_reward_calls_the_checker_once_per_step`: count the calls of a fake `check_polygons` for a 256-answer reward
  call (must be 1).
- `test_reward_values_unchanged`: the same rewards as the old path on the synthetic fixtures (exact equality).
- `test_pool_workers_do_not_init_cuda`: the worker function does not import or call `_gpucheck`.

---

## F6. Switch to the RTX 5090

Change:
- `experiments/phase2/preflight.py:35`: replace `DPH_CEILING = 1.00` with the chosen card's ceiling from
  `policy.CARDS[gpu]["ceiling"]` ($1.50/h for RTX_5090). The card is one argument, with no second copy of the number.
- Launch with `offers.sh RTX_5090 --pick ...` (already supported). US/CA only (fixed rule).
- G-probe and profile timings: mark the 4090 numbers in `project_r1` / preflight as 4090. F7 measures the 5090.

Tests:
- `test_preflight_ceiling_follows_the_card`: RTX_5090 → 1.50, RTX_4090 → 1.00.
- The existing policy tests (`test_vast.py`, renamed in implementation.md B0) stay green.

Box check (F7, first steps):
- `image_check.sh`: CUDA, sm_120, bf16 peak;
- one vLLM start on sm_120 (FlashInfer sampler off, as now);
- `make cuda` builds `native` (sm_120).

---

## F7. Profiling run on the 5090 (measure before R2)

Script: generalize `experiments/phase2/profile_sft.py` into `experiments/phase2/profile.py`, with **SFT** and **GRPO**
arms. Delete the SFT-only file. The GRPO arms run the G-probe trainer for 3 steps on e2, on the shared timing and trace
code.

| Arm | Generation | Micro-batch | Checkpointing |
|---|---|---|---|
| g_vllm_mb8 | vLLM colocate | 8 × 32 | on |
| g_vllm_mb32 | vLLM colocate | 32 × 8 | on |
| g_vllm_mb32_nockpt | vLLM colocate | 32 × 8 | off (may run out of memory: recorded) |
| g_hf_chunk64 | HF, 64 sequences per call (F1-B) | 32 × 8 | on |

Per arm:
- s/step, split into **generate / score / update**;
- peak GiB;
- GPU busy share;
- kernels and syncs in the traced step;
- nvidia-smi utilization and power.

The split comes from TRL's `profiling_context` timings when they are in `log_history`; otherwise we add our own timers
around the three phases.

The chain on the box:
1. `image_check`
2. `make cuda` + `run_box.sh` (F4 check)
3. the GRPO arms
4. the SFT arms of `profile_sft`, only if there is time (no SFT run is planned)

Tests: as `tests/test_p2_profile.py`. Arm configs keep 256 answers per step. The trace parse (`busy`), the split
parse, and the refusal outside timing runs are tested.

Cost: the 5090 price is not checked yet. Run `offers.sh RTX_5090` (read-only, no rental) for a quote. *Estimate*:
~1 h at ≤ $1.50/h, so ≤ $1.50. The quote and the hours go to the owner before `launch`.

Done when: `profile.json` shows the fastest arm that fits in memory. That arm's settings go into `grpo_config` for R2.

---

## F8. Later, only if F7 shows that they matter

| Fix | Issue | Test / check |
|---|---|---|
| torch.compile on the decoder layers only (`dynamic=True`) | RL-3, SFT-4 | F7-style arm; loss of step 1 equal to eager within bf16 noise |
| Flash-attention-2 for the vision tower (image rebuild) | RL-4 | image allowlist; loss parity |
| One vLLM engine start per checkpoint in `select_ckpt` | EV-1 | prompt-id prefix per view; outputs equal to the per-view path |
| Pipeline scoring with the next branch's sampling; pool size from `nproc`; sampler timeout ~3× measured | EV-4, EV-5 | timing log per phase |
| A/B of per-request seeds in vLLM | EV-3 | same prompts, `seed=None` against seeded: throughput |
| Coalesced row-tile stores, edges in shared memory, 32-bit math (bit-exact) | CK-4 | `run_box.sh` 0 mismatches; Nsight Compute sectors per request |
