# P1: RL readiness on an RTX 5090 (2026-11-06)

One Vast rental (offer 48229496, US, RTX 5090 32 GB, CUDA 13.4 driver 615.71, $0.694/h), image
`vllm_native-6008baf45bb0`, chain `experiments/phase2/configs/p1_rl.yaml`, model R1's selected branch e2
(Qwen3-VL-2B + LoRA, merged; SHA-256 `41471a27…b0b9f9d`). Billed **$1.253** (owner's allocation for this test: $1.50).
Raw results: `runs/p1/` (pulled; not in git).

## Card and stack (`env_check.json`)

| Check | Result |
|---|---|
| Capability | (12, 0), sm_120 in torch's arch list |
| bf16 matmul | 238 TFLOPS measured (dense peak in the table: 209.5) |
| Triton | 3.7.1 compiles and runs a kernel |
| torch.compile | 2.13.0+cu130, small MLP forward + backward |
| vLLM | 0.32.0 engine on e2 starts, CUDA graphs captured |
| nvcc | not in the image; `run_box.sh` installs `cuda-toolkit 13.0.3` from pip |

Found on the box and fixed: Triton's `@jit` reads its kernel's source, so `env_check.sh`'s checks now run from a
file, not stdin (commit "env_check: checks run from a file").

## CUDA checker (`cuda_checker/results/cuda.json`)

Built by nvcc for sm_120 (warp-ballot raster, `__int128` paths included). GPU tests: 40 passed (arena reuse over 200
calls, two-stream chunks equal to one call, rewards equal to the recorded values, one checker call per reward call),
2 skipped (they read local data fixtures, which are not shipped to the box). Exactness: **0 mismatches** against polycore on 8,559 fuzz polygons and all 13
defect classes (138,257 legal shapes, every illegal class flagged exactly). compute-sanitizer is not in the pip
toolkit (the step is marked `~`).

Polygons per second (polycore loop and score_batch with 8 threads on the box's 128-core CPU for comparison):

| Vertices | Crop | Batch | GPU checker | polycore loop | score_batch 8t |
|---|---|---|---|---|---|
| 16 | 64 | 1 | 32,787 | 96,665 | 117,786 |
| 16 | 64 | 64 | 1,307,336 | 96,917 | 194,220 |
| 16 | 64 | 1,024 | **7,738,843** | 88,431 | 536,617 |
| 128 | 256 | 64 | 417,839 | 34,594 | 151,407 |
| 128 | 256 | 1,024 | 1,110,802 | 34,528 | 284,011 |
| 1,000 | 512 | 64 | 69,377 | 8,752 | 63,012 |
| 1,000 | 512 | 1,024 | 45,943 | 8,801 | 80,839 |

The first CUDA build (R0 v3, managed memory per call) did 54,710 polygons/s at 16 vertices and batch 1,024; the
arena + pinned + two-stream host side does 7.74 M (141x), and about 30 µs per call at batch 1. Open: 1,000-vertex
shapes at batch 1,024 fall below the 8-thread CPU path (one thread per shape for the witness tests; CK-4, Stage H).
Masks (LIDC 22x21, batch 8,192): 424,732 masks/s against 80,160 for OpenCV contours + CCL.

## RL step: the profile arms (`profile/profile.json`)

32 prompts x 8 answers per optimizer step, 4 steps per arm, step 1 left out; time split from FastGRPO's timers.

| Arm | Exit | s/step | Generate | Rewards (wait) | Old + ref logprobs | Update | Peak MiB |
|---|---|---|---|---|---|---|---|
| **fast** (FastGRPO, 16 x 16) | 0 | **44.3** | 7.6 | 0.003 | 12.3 | 24.4 | 19,357 |
| trl (TRL 1.17 as shipped) | 0 | 80.3 | | | | | 19,245 |
| mb32 (32 x 8) | 0 | 59.1 | 9.5 | 0.003 | 17.7 | 31.9 | 29,105 |
| compile (regional, default) | 0 | 47.9 | 8.6 | 0.003 | 15.5 | 23.9 | 19,295 |
| graphs (reduce-overhead) | 1 | | | | | | |
| no_ckpt | 1 (OOM) | | | | | | 30,146 |
| trace (fast + profiler) | 0 | (not timed) | | | | | 19,409 |

- **FastGRPO is 1.81x TRL as shipped** on the same card, model and data.
- **Rewards are off the critical path:** the prefetch thread finishes scoring 256 answers (one checker call, CPU
  pool) during the logprob passes; the reward call waits 3 ms.
- **On-policy sampling.** vLLM vs trainer token logprob difference, per step:

  | Arm | mean | max | IS ratio mean |
  |---|---|---|---|
  | fast | 0.009-0.010 | 0.32-0.46 | 0.99-1.02 |
  | trl | 0.056-0.082 | 23-25 | 0.65-0.81 |

  TRL as shipped sends vLLM processor-expanded prompt ids with the image, and vLLM's Qwen3-VL processor then sees
  other image placeholders than the trainer (TRL #6294): the samples come from other image features, and the
  importance-sampling correction masks or down-weights a third of them. FastGRPO's one pad per image (checked
  against TRL's ids on every call) gives bf16-level agreement.
- **Micro-batch 32** is slower: each chunk's window is its longest answer, so bigger chunks carry more padding.
- **No checkpointing** runs out of the card's 31.45 GiB beside vLLM's 35 %.
- **Compiled layers** gain 2 % on the update and lose more on the logprob passes (recompiles across widths, grad
  modes and adapter states): below Stage F's 10 % bar, off.
- **CUDA Graphs (reduce-overhead, per layer)** crash: every decoder layer replays the same graph, so a layer's
  output is overwritten by the next layer's replay before backward reads it, and the cudagraph pool conflicts with
  vLLM's sleep-mode allocator (`cumem_allocator.cpp: invalid argument`). Not adoptable with vLLM colocated; manual
  capture (F4) is the only graph route left.

## Where the update step's time goes (`trace`, one optimizer step under torch.profiler)

61.6 s with the profiler on; GPU busy 70 %; 798,774 kernel launches; 31.5 s of CUDA runtime API time; host syncs:
10,151 `cudaStreamSynchronize`, 7,233 `aten::_local_scalar_dense` (`.item()`), 2,904 `aten::nonzero`. Top kernels:
cutlass bf16 GEMMs (4.0 s, 3.4 s, 2.8 s), memory-efficient attention (2.8 s), a float cast copy (2.6 s).

The update is launch- and sync-bound. Next, in order of measured payoff:
1. **C5 precomputed position ids:** Qwen3-VL's `get_rope_index` runs per forward with `.item()` and boolean
   indexing (most of the scalar reads and nonzero calls): compute once per generation batch, slice per chunk.
2. **C10 drop the old-logprob pass:** with one generation per optimizer step the old logprobs equal the update's
   own; about half of the 12.3 s logprob share.
3. **F4 manual graph capture** of the per-bucket decoder forward, after C5 removes the syncs that block capture.

## Decision

`grpo_config` keeps micro-batch 16 x accumulation 16, gradient checkpointing on, no compile, no graphs (test:
`test_grpo_config_keeps_the_p1_measured_choice`). RL runs on `FastGRPO`.
