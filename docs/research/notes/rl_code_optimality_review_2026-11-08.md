# RL code optimality review: one RTX 5090 (2026-11-08)

A review from fresh eyes of the GRPO training path (`src/objdel/train/fast.py`, `configs.py`, `rewards.py`, `serve.py`,
`sample.py`, `experiments/phase2/g_probe.py`, `train_rl.py`), the installed TRL 1.17.1 code that it drives, and the
CUDA checker (`cuda/gpucheck.cpp`, `src/objdel/reward/score.py`, `cache.py`). The goal is heavy GPU use. I changed no
code and ran nothing on a GPU.

Tags: **VERIFIED** = read in this repo or in the installed library (`.venv`: trl 1.17.1, transformers 5.20.1; the
image pins transformers 5.19.0 and torch 2.13.0, so line numbers can differ a little there). **LOG** = computed from
`runs/rl/train_rl/rl.json` (the 770-step full run). **INFERRED** = an estimate that nobody has measured.

## 0. Where the 27.5 s step goes (LOG, means over 770 steps)

| Part | Mean s | Source |
|---|---|---|
| generation (`time/generate_s`) | 8.47 | fit: **6.06 s + 4.1 ms × completions/max_length**, r = 0.92 |
| rest of `_generate_and_score_completions` (processor pass, H2D, ref pass, rewards, metrics) | 5.79 | `batch_s − generate_s` |
| update: 16 micro-batches of forward + backward | 12.73 | `step_time − batch_s`; fit 9.85 s + 0.027 s × mean completion length |
| outside `step_time` (dataloader, optimizer, log) | ~0.5 | `per_step` 27.47 − `step_time` 27.00 |

Facts from the log that drive this review:

- **The answers are short.** `completions/mean_length` = 105 tokens (80–228). The prompt is about 325 tokens: 256
  image tokens plus text (`cpu_gpu_trace.md`, W = 325 + 1536). **About 75 % of the tokens in every forward are prompt
  tokens. Each group of G = 8 rows has the same prompt, so 7/8 of those prompt tokens are repeated work.** The update
  fit agrees: its intercept, 9.85 s at zero completion length, is mostly prompt compute (slope 0.027 s / 256 rows
  ≈ 105 µs per token; 83k prompt tokens × 105 µs ≈ 8.8 s).
- **Update efficiency.** The update does about 1.0 PFLOP of useful work per step: 110k real tokens × about 8.8 GFLOP
  per token (forward, checkpoint recompute and input-gradient backward over 1.41 B non-embedding text parameters;
  INFERRED from the `runs/r3/ckpt/config.json` shapes). In 12.7 s that is ≈ 80 TFLOPS, about 38 % of the 209.5 TFLOPS
  bf16 dense peak.
- **Generation is mostly fixed cost.** 6.06 s of the 8.47 s does not depend on length. This cost includes the sleep and
  wake calls, the weight sync, processor pass 1 and vLLM's input handling. In 182 of 770 steps (24 %) an answer reached
  1,536 tokens: those steps spend 12.40 s on generation, the others 7.25 s.
- Scoring takes 0.005 s (`time/rewards_s`). The checker is not a bottleneck.

## 1. Ranked findings

| # | Where | Issue | Expected gain per step | Confidence |
|---|---|---|---|---|
| 1 | `fast.py:53-96` (FastForward), TRL `grpo_trainer.py:3141-3170, 2826-2860` | **The shared prompt is computed 8 times** in the ref pass and in each update forward and recompute. A shared-prefix forward ([Prefix Grouper](https://arxiv.org/abs/2506.05433)) packs each group as `[prompt | c1 … c8]`: each completion attends to the prompt and causally to itself, and its positions continue from the prompt's end. That gives 32×325 + 256×105 ≈ 37k tokens instead of 110k (2.9× less). The authors report the same outputs and gradients (verified here: LoRA dropout is 0, `configs.py:95`; TRL disables dropout). It also cuts activation memory about 3×, which may let checkpointing go (see #5). | update 12.7 → ~5 s, ref ~3 → ~1.2 s: **~8–10 s (30–35 %)** | high on the token arithmetic (LOG + VERIFIED shapes); medium on the realized time. **High effort.** |
| 2 | `configs.py:192-193` (`vllm_enable_sleep_mode: True`); TRL `vllm_generation.py:500-503, 719, 772` | **Sleep and wake on every step.** Level 2 discards the vLLM weights and KV cache. Each step then runs `empty_cache()` (frees the whole torch cache, so the next ref pass and update call cudaMalloc again), wakes the weights, reloads **all** weights (vision tower and embeddings too, though LoRA never changes them), wakes the KV cache and sleeps again. These are the 2–4 s idle gaps you measured. With sleep off, vLLM's 0.35 × 32 GB stays resident. Training peak (17 GB) + vLLM (11.2 GB) ≈ 28 GB of 32 GB. That probably fits, but the margin is thin. vLLM's real KV need is small: 256 × ~430 tokens with prefix caching ≈ 37k tokens ≈ 4 GiB KV + 4 GiB weights. So `vllm_gpu_memory_utilization` 0.28–0.30 probably suffices. | **1–3 s** | medium. It depends on what the 17 GB peak measured (torch peak or nvidia-smi); test once. Low effort (one config key). |
| 3 | `configs.py:188` (`beta: 1e-3`); TRL `grpo_trainer.py:2826-2857` | **The ref pass exists only for a KL term of weight 0.001.** It is a full no-grad forward over 256 rows, ≈ 0.32 PFLOP, about ¼ of the step's training FLOPs. With LoRA, beta = 0 loads no extra model, and switching back to beta > 0 later is cheap (adapter disable), unlike the full-FT crash case named in the `configs.py` docstring. Dr. GRPO and DAPO train with beta = 0. **This is a method decision for the owner, not a pure speed fix.** | **~3–4 s** (INFERRED: the ref share of the 5.8 s) | medium |
| 4 | `fast.py:29` (BUCKETS), `fast.py:43` | **Bucket rounding pads every chunk, and nothing uses the buckets.** They exist for torch.compile's shape limit, but compile is off in `train_rl` (`g_probe.py:95` runs only with `--compile`, and `train_rl.py` has no such flag). The rows are 325 + ~105 ≈ 430 tokens, so a chunk whose longest row is 189–443 tokens pads to 768 instead of 514–768. A chunk at 60–187 pads to 512. INFERRED average waste: 10–20 % of ref + update FLOPs. Fix: round the width up to a multiple of 32 or 64 (`FastForward.buckets = tuple(range(256, 2049, 64))`; `window()` already accepts any tuple). Log `padded tokens / real tokens` to measure it. | **~1–2 s** | high on correctness; medium on size (no per-row lengths in the log). **Trivial effort.** |
| 5 | `configs.py:186` (`gradient_checkpointing: True`, every one of the 28 layers) | **Full recompute although the peak is 17 of 32 GB.** P1 hit OOM with no checkpointing at all (`p1_results.md:61`, 30.1 GB). Partial checkpointing fits in between. transformers 5 has a per-layer flag (`GradientCheckpointingLayer.gradient_checkpointing`, `modeling_layers.py:77`), so `layer.gradient_checkpointing = False` on k of the 28 layers removes k/28 of the recompute, which is ~⅓ of the update FLOPs. This conflicts with #2 for the same memory. After #1, both probably fit. | ~1–1.5 s (k ≈ 10) | medium; tune k by measurement |
| 6 | TRL `vllm_generation.py:520-522` | **About 700 separate `load_weights([(name, param)])` calls per step.** TRL does one call per tensor for the FSDP or ZeRO gather. With PEFT merged on one GPU, one call over the whole generator is the same operation. vLLM's model `load_weights` builds its loader and name mapper on each call (INFERRED; vLLM is not installed in the local venv). With sleep off (#2), you can also skip the tensors that never change (`visual.*`, embeddings, norms): about ⅓ of the bytes. **First put a timer on `sync_weights`.** It is the largest unmeasured part of the 6 s generation intercept. | 0.2–1+ s | low-medium (unverified per-call cost) |
| 7 | `fast.py:292-309` (ImageMemo) + TRL `grpo_trainer.py:1779-1806` (pass 1), `2645-2646` (pass 2) | **Two 1.61 GB fp32 CPU concatenations and one pageable 1.61 GB H2D copy per step, for 32 distinct images (0.2 GB).** ImageMemo dedups the processing but still builds all 256 rows on the CPU in both passes, and pass 1's pixels are then thrown away. Fix: keep the 32 per-image rows on the GPU and build the 256-row tensor there with `index_select` (≈1 ms D2D). TRL's `.to(device)` is then a no-op. Keep `image_grid_thw` on the CPU: the Qwen processor reads it in Python to expand the image pads. A smaller fix is bf16 rows plus pinned memory (bit-identical, because the tower casts to bf16, and VisionCache already hashes int16). | 0.3–0.8 s of CPU time with the GPU idle | medium (copy rates INFERRED) |
| 8 | TRL `vllm_generation.py:677` (colocate: `n: 1`, 256 requests) | **Colocate sends 256 requests with n = 1. TRL's server path sends 32 requests with n = 8** (`:618-631`). vLLM must then hash and look up 256 images, run the input processing 256 times, and may prefill or encode a duplicate prompt twice when it schedules both in the same step (chunked prefill, 4,096 tokens). FastGRPO already wraps `llm.generate` (`fast.py:469-489`). That wrapper can group runs of G identical rows, call with `n = G`, and expand back into 256 per-row outputs. | 0.2–1 s | low-medium (vLLM internals unverified) |
| 9 | ref pass order: TRL `grpo_trainer.py:2850-2857` (`batch_size` = 16, rows in sampler order) | **The ref pass is neither sorted by length nor batched wider.** `by_length` (`fast.py:225`) only reorders the update. In the ref pass, one long row widens its whole chunk of 16. The pass runs without grad, so it can use chunks of 32–64 rows, with the rows sorted by length and the order undone at the end. It also launches fewer kernels. | 0.3–0.6 s | medium |
| 10 | `fast.py:131-148` (VisionCache) | **The vision tower runs 16 times on 2 images, not once on 32.** The ref pass meets 2 new images per chunk. Each encode does about 26 host syncs (`lengths.tolist()` in 24 blocks, `modeling_qwen3_vl.py:299`) and gives the GPU a small 2k-token batch. Encode all misses once, when FastForward is entered (one `keys()` over the whole batch, one encode). | 0.1–0.3 s | medium |
| 11 | `configs.py:75` (`max_completion_length` 1,536) | **The runaway tail.** A truncated answer is masked from the loss (`mask_truncated_completions`), yet it holds the decode to 1,536 steps: +5.15 s in 24 % of steps, ≈ 1.2 s/step on average. A cap of 1,200 (the LIDC target maximum is 1,182) saves ~0.35 s/step. In 23 of 770 steps, a terminated answer was longer than 1,200 tokens and would then count as truncated. **Owner decision; small gain.** | 0.3–1 s | high (LOG) |
| 12 | TRL `_compute_loss` `grpo_trainer.py:3376-3405`; fused head `trl/trainer/utils.py` (`hidden_states[mask]`); `fast.py:39-40, 62, 129, 133, 181` | **About 29 host syncs per micro-batch** (`cpu_gpu_trace.md` A10). The 7 `.item()` calls in `_compute_loss` drain the GPU before backward is queued. **The real cost is low.** Each micro-batch is ~0.7 s of GPU-bound work (16 × ~512 tokens, about 6k launches at 5–10 µs ≈ 30–60 ms of CPU). After a sync, the GPU loses only about one launch latency. Fix only if cheap, for example by keeping TRL's loss metrics as device tensors. | < 0.3 s | medium-high |
| 13 | text tokenization / `apply_chat_template` ×256 in both passes | 32 unique prompts are rendered and tokenized 512 times per step. A memo keyed by prompt text would avoid this. | 0.1–0.4 s | low |
| 14 | `train_rl.py:44` | `fork_pool(a.procs)` forks up to 16 workers that `train_rl` never uses. It costs memory only (copy-on-write STORE), no step time. Cleanup, not speed. | 0 | high |

Rough total if #1–#7 all land: from ~27.5 s to ~12–15 s per step. The items overlap, so they do not simply add. This
is INFERRED and needs measuring one item at a time.

Suggested order: #4 (one line) → #2 (one config key, one memory test) → timers on `sync_weights`, wake/sleep and the
ref pass, to split the 6 s and 5.8 s parts → #3 (owner) → #7, #9, #10 → #1 (the large refactor) → #5 after #1.

## 2. Notes on #1 (shared-prefix forward), the largest lever

- Rows arrive as contiguous groups of 8 identical prompts. `RopeCache.keys` (`fast.py:178`) already finds identical
  rows by exact ids and mask.
- Positions: take the prompt's 3-D mrope positions from RopeCache once per group. Completion token t of every member
  gets position `prompt_end + t` on all three axes. `RopeCache.pending` already injects position ids.
- Attention: [prompt causal] + [each ci → prompt + causal within ci]. Use a 4-D boolean mask under SDPA, or
  `flex_attention` with a block mask. Prefix Grouper's split (prefix self-attention, then suffix → prefix + self) also
  fits FlashAttention varlen.
- Logits: the first token of each ci is predicted from the **prompt's last hidden state**, not from the previous
  packed token. So gather the predecessor hidden state per completion token before the fused log-prob kernel. Do not
  call the model's fused head on the packed sequence as it is.
- DeepStack injects features at image positions, which are only in the prompt: those are shared too.
- Micro-batches become whole groups (2 groups = 16 rows, as now), sorted by group length instead of by row length.
  Truncated rows have `completion_mask` = 0 and can be dropped from the packing entirely.
- Test: on the CPU with a tiny Qwen3-VL config, compare logps and LoRA grads with the current path (`torch.allclose`).

## 3. CUDA checker and scoring: no action needed (VERIFIED)

- One checker call per step for all 256 answers (`score.py:115`; `rewards.py:119`, prefetched on a thread right after
  generation, `fast.py:443-453`), so it overlaps the processor pass and the ref pass.
- Per call: polygon path 2 lanes × (1 H2D, 1 memset, 5 kernels, 1 D2H, 1 stream sync); scoring path (`score_words`,
  `gpucheck.cpp:1549-1640`) ≈ 20 launches and 3 stream syncs. That is ≈ 30 launches per step, independent of the
  answer count.
- Pinned staging buffers and device arenas grow only (powers of two), so a call allocates nothing in steady state
  (`gpucheck.cpp:120-140`). Streams are non-blocking and lowest priority (`:66-72`), so trainer kernels go first.
  Syncs are on the checker's own streams only.
- Tables upload once per cache per process (`cache.py:241-258`): a cudaMalloc + sync the first time each item is
  scored (epoch 1 only, on the prefetch thread).
- CUDA Graphs would not help: the sizes change on every call, the host needs counts between stages, and the whole
  call takes ~5 ms. `SigScan` runs on one thread (`for_each(st, 1, …)`) over ~256–512 objects, which is trivial.

## 4. Already good (VERIFIED, so do not redo)

OpenOnce (32 PNG decodes per step instead of 4,096); ImageMemo dedup; C10 skips the old-logprob pass; VisionCache and
RopeCache; `by_length` for update micro-batches; `logging_nan_inf_filter` off; `floating_point_ops` returns 0; LoRA
input cast removed under autocast; IS statistics kept on the device; vLLM CUDA graphs plus compile on
(`enforce_eager=False` in `runs/rl/rl_full.log:39`); prefix caching on; fused Triton log-prob/entropy LM head (TRL
`add_fused_lm_head`, only label positions scored); bf16 base (`dtype` "auto" reads bf16 from the checkpoint config).

## 5. Things that look like wins but are not

- **TF32** (the warning in `rl_full.log:1`): the matmuls run under bf16 autocast, so TF32 changes nothing.
- **`expandable_segments`**: not set for RL (`train_sft.py:264` sets it only for SFT). With vLLM sleep mode it is
  likely refused (vLLM's CuMem allocator does not support it; INFERRED). Consider it only after #2.
- **More dataloader workers / pin_memory**: OpenOnce needs `num_workers = 0` (`fast.py:337-339`). The 32 decodes per
  step are ~0.1–0.2 s.
- **Removing every sync**: see #12. The micro-batches are GPU-bound, so the GPU idle time comes from the phase-level
  gaps (generation fixed cost, CPU phases), not from per-op syncs. nvidia-smi "util" counts any running kernel as busy.
  Track MFU (≈38 % in the update now) and the phase times instead.
