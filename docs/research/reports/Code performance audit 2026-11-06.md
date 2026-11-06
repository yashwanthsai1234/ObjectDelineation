# Code performance audit (2026-11-06)

Scope: the training path, the evaluation and sampling path, and the CUDA checker. Owner question: "why is the GPU at
~60 % in R1, and what are we not optimizing?" Sources:

- three read-only code audits, against transformers 5.20.1, trl 1.16.0, peft 0.21.2 and torch 2.14.1 (the box image
  pins the same versions);
- CPU measurements on the Mac;
- the research report `NVIDIA CUDA performance for VLM training.md`.

No GPU measurement in this audit yet. The GPU numbers below are estimates until `profile_sft` runs on a box.

## 1. Measured on the Mac (CPU)

| What | Result | Meaning |
|---|---|---|
| Data path per micro-batch of 4 (render + TRL collator, Qwen3-VL processor) | 19 ms (9 render + 10 collate) | 0.15 s of CPU work per optimizer step, spread over 6 workers. The GPU step is 3.4 s. **The data loader is not the limit.** |
| Padding in a length-grouped batch | 10 % | small |
| Loss tokens with TRL's training chat template | answer tokens + `<|im_end|>` + newline | the loss mask is correct; the model learns to stop |
| `view_rows` (98,804 rows) | 12.7 s once per run | negligible |

R1's manifest gives peak CUDA memory 6.1 GiB for every phase. The "20.9 GiB" comment in `train_sft.py` was from the LFM
full fine-tune; it is now corrected.

## 2. Training path: ranked findings

Budget (audit estimate): about 195 TFLOP per optimizer step. At 3.4 s that is about 57 TFLOPS average, or about 95 TFLOPS
while the GPU is busy. The bf16 dense peak is 165 TFLOPS. The idle 40 % is CPU dispatch and host syncs.

| # | Finding | Evidence | Estimated effect | Fix |
|---|---|---|---|---|
| T1 | Micro-batch 4: the CPU launches the same kernels per micro-batch at any batch size, so 8 micro-batches cost 8× the CPU work. R1's reason ("152k-vocab logits ~12 GB") is out of date: TRL 1.16 scores only label positions in chunks (`trl/trainer/utils.py:1548-1660`). | `r1_qwen3vl2b.yaml:11-12` | −30 to −45 % step time | micro-batch 8 × 4 or 16 × 2 (effective batch stays 32) |
| T2 | Gradient checkpointing on with ~17 GB free: recomputes the LM forward (~26 % of the compute, ~⅓ of the launches) | `configs.py:72,164` | −15 to −25 % | off, or `every_n_layers: 2` at micro-batch 8 |
| T3 | 60–80 host syncs per micro-batch: isnan/isinf per micro-batch, TRL metric `.item()`s, Qwen3-VL `get_rope_index` `.tolist()`/`.item()`, vision `.tolist()` per layer | `trainer.py:1891-1894`; `sft_trainer.py:1640-1653`; `modeling_qwen3_vl.py:984-1021, 297-299` | −5 to −15 % | `logging_nan_inf_filter: false`; keep the metrics on the GPU until a log step; position ids in the collator |
| T4 | Vision attention on sdpa: one call per image and one `.tolist()` per layer (24 layers) | `modeling_qwen3_vl.py:296-315` | −3 to −6 % | flash-attention-2 for the vision tower only (needs an image change) |
| T5 | fp32 LoRA adapters add casts in 196 layers | peft `layer.py:1100,1113` | −3 to −8 % | let torch.compile fuse them (keep fp32 for AdamW) |
| T6 | No torch.compile | — | −15 to −30 % on top of T1–T2 | compile the decoder layers only (`dynamic=True`) |
| T7 | FLOP counting walks ~2.5k PEFT modules per micro-batch | `trainer.py:1902, 4206-4209` | ~−1 % | `floating_point_ops` → 0 in the Trainer subclass |
| T8 | Blocking host-to-GPU copies (~25 MB fp32 pixels per micro-batch) | `trainer_pt_utils.py:1518` | <1 % | `accelerator_config: {non_blocking: true}` |

Checked and fine:
- fused AdamW (the default);
- grad-norm `.item()` only on log steps;
- data workers, pin_memory and prefetch;
- padding;
- the LM head (chunked Triton, ~7 % of compute).

Packing and padding-free are refused by TRL for vision data. Liger is not installed, and TRL forbids its fused
cross-entropy.

Expected: T1 + T2 → ~1.8–2.2 s/step; + T3 + T6 → ~1.4–1.6 s/step (from 3.4).

## 3. Evaluation and sampling path: ranked findings

| # | Finding | Evidence | Estimated effect | Fix |
|---|---|---|---|---|
| E1 | One vLLM engine start (~95 s) per sampler call, 8–9 per select run | `select_ckpt.py:136-163, 195-207` | ~14 min per select run | one child process per checkpoint for all its views |
| E2 | HF check: 512 answers × ~1 s, launch-bound, host stop check every token, DynamicCache `torch.cat` | `select_ckpt.py:284-302`; `hf_sample.py:62-79` | 8–9 min | larger batch, static cache (test first) |
| E3 | Per-request seeds in vLLM may sample row by row (unverified) | `sample_vllm.py:44,52` | up to 1.5–3× decode | A/B with `seed=None` (loses batch-independent reruns) |
| E4 | CPU fp32 `vision_check` and scoring run between branches while the GPU waits | `select_ckpt.py:335-344, 379-383` | 2–8 min | overlap with the next branch's sampling |
| E5 | CUDA checker called once per item from 4 processes (4 GPU contexts) | `metrics.py:246-253`; `select_ckpt.py:217` | unmeasured | batch items in one process |
| E6 | Pool fixed at 4 processes; sampler timeout ~2.5 h on a hang | `select_ckpt.py:361, 58, 206` | failure cost | size from `nproc`; timeout ~3× the measured time |

## 4. CUDA checker: why it is slower than the CPU

The kernels are not the problem. Memory handling is:

| # | Finding | Evidence | Fix |
|---|---|---|---|
| K1 | ~20 managed buffers per call, each `cudaMallocManaged` + `cudaMemset` + `cudaDeviceSynchronize`, then `cudaFree`: ~40 syncs per call | `gpucheck.cpp:85-103, 809-815` | one device arena kept across calls; one sync per call |
| K2 | No `cudaMemcpy` at all: every byte moves by unified-memory page faults (~0.6 GB/s measured in RESULTS) | `gpucheck.cpp:89,95, 823-825` | `cudaMalloc` + pinned staging + `cudaMemcpyAsync` |
| K3 | Returns two full u8 rasters (even-odd and nonzero); `to_masks` uses only even-odd; host copies again | `gpucheck.cpp:813, 823-824`; `gpucheck_py.cpp:33-36` | even-odd only, bit-packed or bbox-cropped |
| K4 | One thread per row: uncoalesced stores, `buf[256]` in local memory | `gpucheck.cpp:221-236` | a block per row tile, consecutive 32-bit words |
| K5 | Calls of 8 shapes from forked workers | `score.py:159,173` | one call per reward step / per view |

The RESULTS baseline is not the same work: `score_batch` returns areas only. The fair baseline is the polycore loop, and
against it the GPU is 1.06× at 16 vertices / batch 1024 and 0.18–0.19× at larger shapes. After K1–K3 the expected
fixed cost falls from ~8 ms to ~50 µs, and the GPU should pass 8 CPU threads at batch ≥ 32–64. At batch 8 the gain
stays small, because `score_mask` still runs on the CPU per object.

## 5. GPU profiling plan (next rental, not started)

`experiments/phase2/profile_sft.py` runs R1's settings as 16-step timing runs, one process per arm. The arms are:
- base
- mb8
- mb16
- mb8_lean (no isnan sync, non-blocking copies)
- no_ckpt
- no_ckpt_mb8
- no_ckpt_compile

It runs torch.profiler on one step after warm-up (the new `train_sft` `profile:` option) with nvidia-smi alongside.
`profile.json` gives per arm:
- s/step
- peak GiB
- GPU busy share of the traced step
- kernel count
- CUDA API time
- host-sync calls
- top kernels

Command: `run.sh p0 experiments.phase2.profile_sft`. Estimate: ~25 min of arms + boot and sync ≈ 0.6 h; at R1's
~$0.5/h that is ≈ $0.30 (cap 0.75 h, ≤ $0.40). Preflight first.

Not in this plan:
- nsys/ncu: not in the image, and counters need host permissions on Vast.
- Liger: not installed.
- vision flash-attention: needs an image change.
- The CUDA checker profile: after K1–K3, with `nsys --cuda-um-cpu-page-faults=true --cuda-um-gpu-page-faults=true`.
