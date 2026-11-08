# CUDA streams, syncs and graphs: review 2026-11-08

Read: `fast.py`, `rewards.py`, `gpucheck.cpp`, installed TRL 1.17.1 and transformers 5.20.1 (the box has 5.19).
**[V]** = verified in code, a log or the docs. **[I]** = estimate. No code changed.

## 0. R5 numbers (shared_rl.json, steps 3-15) [V]

step 23.43 s; generate 11.11 (sync 1.27, wake 0.11 and sleep 0.29 included); batch − generate 4.65; update 7.67.
Two-point fit with `packing_wins` costs (base 1.0, shared ≈0.40): 12.37 = a + b, 7.67 = 0.4a + b, so **b ≈ 3-5 s does
not scale with FLOPs** [I]. Cross-check: 2,744 tokens × 2.8 GFLOP × 4 passes ≈ 0.2-0.3 s of GPU time per
micro-batch, against 0.48 s measured. After packing, fixed cost per micro-batch (launches, syncs, small kernels)
dominates the update.

## 1. Syncs per step (16 ref chunks + 16 micro-batches)

| Site | What | /step | Avoidable? |
|---|---|---|---|
| fast.py:74, :125 | `same.tolist()`, `lengths.tolist()` | 64 | yes: one stacked read |
| fast.py:196-197 | `torch.tensor(..., device=)`: pageable blocking H2D | 64 | yes: build on the device |
| fast.py:111 | `grid.prod(-1).tolist()` | 17 | yes: host grid |
| fast.py:324 | RopeCache keys `.cpu()` ×2, `grid.tolist()` | 51 | partly |
| fast.py:267, 272, 276 | VisionCache blocking `w.to`, `h.tolist()`, grid ×2 | 128 | 96 of them |
| qwen3_vl:297-299 | `lengths.tolist()` in 24 vision blocks × 16 encodes of 2 images | 384 | yes: one 32-image encode → 24 |
| qwen3_vl:853-854 | DeepStack boolean index, 3 layers | 192 | only by patching transformers |
| qwen3_vl:1044, 1078 | `split_sizes.tolist()`, placeholder check | 64 | check: `TRANSFORMERS_DISABLE_TORCH_CHECK=1` |
| trl utils.py:1687-1697 | fused head `[mask]` | ~80 | keep (saves head FLOPs) |
| grpo_trainer.py:3384-3405 | 7 loss `.item()` before backward | 112 | yes, needs a TRL override |
| grpo_trainer.py:2387-2405, 1706, 1734, 2941-2972 | metrics, rewards | ~90 | no: GPU is idle there |
| fast.py:649 | `timed()` synchronize ×7 | 7 | keep (<10 ms) |
| vllm_generation.py:501 | `empty_cache` | 1 | R4: sleep off saved only 0.05 s |

Total ≈1,250 [V count; I cost 0.25-0.6 s]. `wake_s` nests inside `sync_s`: do not add them.

## 2. Copies and side streams

- **Correction (`grpo_utilization_research` §4):** the 1.61 GB fp32 `pixel_values` is not pinned or non-blocking.
  ImageMemo concatenates it on the CPU (fast.py:473), then `Trainer._prepare_input` calls a plain `.to(device)`
  (trainer.py:2356). `accelerator_config.non_blocking` reaches only the DataLoader (PIL rows here). Pinning would
  not help: the GPU is idle then, and PyTorch's tutorial says pin-then-copy "can be up to twice as slow". Upload 32
  images instead (prior review #7).
- **Checker: already optimal.** It runs on a prefetch thread with no torch calls, on lowest-priority non-blocking
  streams (gpucheck.cpp:66-71) with pinned staging. It waits 5 ms.
- **Weight sync: streams do not help.** It is about 700 Python `load_weights` calls (vllm_generation.py:522); batch
  them (prior review #6).
- **GPU-idle window:** about 1.5-2 s of the 4.65 s is CPU work (processor pass 2, the pixel cat and copy, RopeCache)
  [I]. The 32-image vision encode, queued from ImageMemo's pass-1 rows before pass 2, could overlap it (0.2-0.4 s) once
  the vision syncs are gone.
- No side streams in the update or ref pass: backward runs on its forward's stream [V CUDA semantics], and the passes
  depend on each other.

## 3. CUDA graphs

- vLLM decode is graphed (FULL_AND_PIECEWISE, `enforce_eager=False`, r5_shared.log) and keeps its graphs across sleep.
- **Correction (P1):** the reduce-overhead crash was cudagraph trees' "output ... overwritten by a subsequent run"
  (graphs/log.txt:177). It came in the no-grad old-logprob pass: 28 layers share one compiled function. The
  `cumem_allocator.cpp:249` error prints *after* the traceback. Line 249 is `cuMemRelease` in vLLM's free path, so
  it is teardown. **No sleep-allocator conflict is shown.**
- Real blockers: packed widths change per micro-batch (pad to 256-multiples, about 14 shapes); graphed backward keeps
  activations in static pools, so it cannot run with checkpointing; pools survive `empty_cache`, so less memory at
  vLLM wake; the ref pass (adapter off) needs its own graphs or stays eager; autocast `cache_enabled=False` and no
  hooks [V torch docs]. The syncs sit outside the decoder layers, so per-layer capture is allowed.
- Graphable: decoder layers of update micro-batches, 0.5-1.5 s after change #1 [I].

## 4. Verdict and top 3

The checker and the prefetch use streams well; only vLLM uses graphs. The trainer runs on one stream with about 1,250
syncs, and fixed cost per micro-batch now dominates. **Not optimal.**

| # | Change | s/step | Risk | Lines |
|---|---|---|---|---|
| 1 | `--grpo '{"per_device_train_batch_size":32,"gradient_accumulation_steps":8}'`: 4 groups per micro-batch. Same update under dr_grpo (sum/(256·L)); spg 8, 256 rows; C10 holds. Halves micro-batches and ref chunks. P1's mb32 loss was padding, now gone. | **1.3-2.5** | low-medium (memory: one arm) | 0 + P1 pin test |
| 2 | `make_graphed_callables` per layer and width bucket, update only, checkpointing off | 0.5-1.5 | high | ~200 |
| 3 | Sync cleanup: one host read per call, `seg`/`idx` on the device, cache `w`, one 32-image encode, torch-check env | 0.3-0.8 | low | ~40 |

Order: 1 → 3 → measure → 2. The savings overlap. First add a `time/ref_s` timer: the ref/CPU split of the 4.65 s is
not measured.
