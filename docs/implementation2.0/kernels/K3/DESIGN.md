# K3: GPU kernels for the GRPO loop. Image patches, rope positions, fused log-prob and IS ratio

Folder: `docs/implementation2.0/kernels/K3/` (SP = this session's scratchpad). REPO was not edited. No GPU was rented.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Checked the live catalog (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`, saved as `docs/implementation2.0/notes/doc/nv_catalog.json`). Matches, **recommended only, not installed**: `compute-sanitizer-diagnose-cuda-bugs` (GPU Development), which became step 6 of GPU_RUN.md (memcheck + racecheck); `dali-dynamic-mode` (Training AI; NVIDIA DALI GPU image pipelines), listed as an alternative for (a) in the decision record; `tilegym-monkey-patch-kernels-to-transformers`, a pattern for wiring kernels into HF models (integration, section 8). |
| engineering:code-review | Self-review of every kernel. Found and fixed: a 16-byte local array in the HWC preprocess path (ptxas showed a stack frame; now shifts, 0 bytes); a `-std=c++17` flag that breaks torch 2.14 extensions (its headers need C++20); two test bugs (ulp check at the overflow edge, inf == inf in a stats check). |
| engineering:debug | The mismatch-isolation plan in section 10 (reproduce on one image or token, then bisect by stage: value map, index map, split, lane). |
| engineering:system-design | The requirements, data-flow diagram and load numbers in sections 1-2; trade-offs written down per kernel. |
| engineering:architecture | The decision record in section 9: custom CUDA vs Triton vs torch ops vs a library, per kernel, with a rule for when to drop K3 code. |
| engineering:documentation | GPU_RUN.md as a runbook (prerequisites, steps with times, pass criteria, failure table). |
| huggingface-skills:huggingface-best | No change. The model is fixed (Qwen3-VL-2B, owner decision), so the leaderboard workflow was not run. |
| huggingface-skills:huggingface-papers | Paper lookups through the HF papers API; the IDs below were checked there (2410.10989, 2411.09009, 2510.26788). |

Web research: NVIDIA Blackwell Tuning Guide and CUDA Programming Guide (compute capability 12.0 limits), the TIS blog,
the SkyRL off-policy docs, TRL issue #4493, the CCE paper and code, and the FP16 mismatch paper (citations in section 11).

## 1. Problem

One R5 GRPO step takes 23.4 s, and the GPU is idle about 40 % of it. The CPU phase alone takes 3.0-3.2 s [measured:
BRIEF §2, from R5]. Three pieces of it are in scope here:

- **(a) Image patches.** TRL runs the Qwen3-VL image processor twice per step. Each pass builds 262,144 × 1,536 fp32
  `pixel_values` (1.61 GB) on the CPU. The second pass sends them to the GPU as a blocking, pageable copy
  [measured: REPO/docs/research/notes/cpu_gpu_trace.md rows 2.3, 4.5]. One processor pass measured 0.75-1.31 s on the
  Mac [same file, row 2.3, MEASURED-MAC].
- **(b) Log-prob + IS ratio.** Each update micro-batch needs per-token log-probs over a 151,936-word vocabulary, then
  the vLLM importance ratio. TRL's fused head runs cuBLAS on (≤4,096 × 32,768) tiles into a buffer, then Triton
  online-logsumexp [code: .venv/.../trl/kernels/chunked_logprob.py:474-545]. The ratio and its 7 logged stats are
  about two dozen small torch ops [code: REPO/src/objdel/train/fast.py:351-362, 581-602]. The ratio's minimum is 0 in 98.1 %
  of 770 steps [measured: BRIEF §2, runs/rl/train_rl/rl.json].
- **(c) Rope positions.** `RopeCache.positions` reads ids and masks to the host, runs HF `get_rope_index` as a Python
  loop per row, and copies 11.4 MB back. That costs 8 host syncs, plus 5 syncs per micro-batch for its keys
  [measured: cpu_gpu_trace.md rows 6.3, 10.4]. Yet every prompt is the same text with the same 512 × 512 image.

Requirements: (a) **bit-equal** to the HF fast image processor in fp32; (b) never builds full-vocab fp32 logits, does
token-level TIS (and TRL's other three modes); (c) positions computed once per prompt layout. All three run on
sm_120 (RTX 5090, 32 GB) and make no host syncs. Non-goals: the backward pass of (b) (section 9 says why), and
training changes.

## 2. Status (what is proven now, and where)

| Claim | Evidence |
|---|---|
| (a) K3 CPU reference == HF processor, fp32 and bf16, 41 images (32 random, 0, 255, all-values, 3 digit images, checkerboard, ramp, grayscale-as-RGB): **0 of 64.5 M values differ** | [measured: results/test_hf_cpu.log.txt] (macOS, transformers 5.20.1, torch 2.14.1) |
| (a) host emulation of the GPU index math == the HF-structured reference (512² and 96 × 64, HWC and CHW) | [measured: results/test_cpu.log.txt (run on the Mac and in the linux/amd64 nvcc container)] |
| (c) template + extrapolation == C++ port == HF `get_rope_index`, real prompt (Lp 325), 16 rows × 1,868, answers 0..1,536 tokens, a mask hole; 3 bad rows flagged | [measured: results/test_hf_cpu.log.txt]; plain-torch variant == HF on 12 rows (one-off check, this session) |
| k3_expf ≤ 0.908 ulp (2 M point sweep), k3_logf ≤ 1.10 ulp (2 M random floats) | [measured: results/test_cpu.log.txt, results/docker_nvcc_check2.log.txt] |
| (b) CPU reference (the order the GPU replays) within 3.6e-6 (logp) and 7.4e-6 (entropy) of fp64; temperature 0.7 checked | [measured: results/test_cpu.log.txt] |
| Every kernel compiles for sm_120 with nvcc 13.0.88 (the toolchain the box installs); 0 spill bytes; lse 128 regs (tensor cores) / 81 regs (exact), 33,792 B static smem | [measured: results/docker_nvcc_check2.log.txt] |
| GPU bit-equality, HF/TRL comparisons on the GPU, timings | **not run yet**: GPU_RUN.md, steps 3-5 |

## 3. Data flow

```
 decoded PNGs (PIL, uint8 HWC)          prompt ids (CPU, from the processor)
        |  32 distinct images / step           |  one row, once per prompt layout
        v                                      v
 [pinned uint8 stack 25.2 MB] --H2D--> [K3 preprocess] --> pixel_values (GPU, 256 rows, fp32 or bf16)
                                       row -> image idx          |
                                                                  v
 input_ids / mask / mm types (GPU) --> [K3 rope] --> position_ids (3,B,W) --> RopeCache pre-hook --> model forward
                                       template (3 x Lp), end                                         |
                                                                                    final hidden states (bf16)
                                                                                                      |
                     no-grad passes (reference, eval): [K3 lse] --> logp, entropy                     |
                     grad pass (update): TRL ChunkedLogProbFunction (cuBLAS + Triton) --> logp <------+
                                                                  |
                     vLLM sampling logps, mask --> [K3 is_ratio] --> ratio (B,T) or (B,1) + 7 stats (device)
                                                                  |
                                                     FastGRPO.log: one host read per step
```

## 4. Algorithms and complexity

### (a) Preprocessing: `k3gpu::preprocess` (k3_kernels.cu)

HF's path for a 512 × 512 image does not resize: `smart_resize(512, 512, factor 32)` returns (512, 512), and
torchvision returns the input unchanged [code: torchvision/transforms/v2/functional/_geometry.py:313-314]. HF fuses
rescale into normalise: mean' = fp32(mean) × fp32(1/rescale) = 127.5, and std' = 127.5 likewise
[code: transformers/image_processing_backends.py:300-340]. Then torchvision computes `sub_(mean').div_(std')` in fp32
[code: torchvision/.../_misc.py:62-67]. Patchify is a reshape → permute(0,2,5,3,6,1,4,7) → expand(temporal 2) →
reshape [code: transformers/models/qwen2_vl/image_processing_qwen2_vl.py, `patchify`].

The kernel fuses everything into one gather pass. Each thread writes 4 consecutive features (same channel, temporal
copy and patch row = 4 adjacent pixels). It reads those pixels with three 4-byte loads (HWC) or one (CHW), applies
`k3_norm`, and writes one float4 (fp32) or one 8-byte pack (bf16 RNE). An optional `idx[row]` picks the source image,
so the generation batch's 256 rows (8 per image) come out in row order with no concat step.

- Work: n_out × 1,024 × 1,536 outputs. 256 rows = 402.7 M values. Writes 1.61 GB (fp32) or 0.81 GB (bf16); reads
  25.2 MB. For 32 distinct images: 201 MB fp32.
- Memory bound. At the RTX 5090's 1,792 GB/s datasheet bandwidth [third-party spec page, section 11], 256 rows take
  ≥0.9 ms in fp32 and ≥0.45 ms in bf16 [estimate]. The H2D of 25.2 MB pinned uint8 replaces 1.61 GB pageable fp32.
- bf16 output: the vision tower casts `pixel_values` to its weight dtype first
  [code: transformers/models/qwen3_vl/modeling_qwen3_vl.py:96-101]. A bf16 output (RNE, like torch's `.to`) gives
  the conv exactly the bits it sees today, at half the memory.

### (c) Rope positions: `k3gpu::rope_positions`

For a row whose real tokens (mask ≠ 0, in column order) are the template prompt (Lp tokens) followed by text tokens,
HF's grouping loop [code: modeling_qwen3_vl.py:925-1031] gives real token r the position `tmpl_pos[a][r]` for r < Lp,
and `tmpl_end + (r − Lp)` after. The answer extends the prompt's last text group, whose positions are an arange on
all three axes. Masked columns get 0 (HF starts from zeros). One block runs per row, with a CUB `BlockScan` of the
mask for r. Each thread checks its token against the template (ids and types for r < Lp; text type for r ≥ Lp) and
row length ≥ Lp. Failures OR into an error word (bits 1, 2, 4), so a wrong prompt is never silently positioned.

- Work: O(B·W) = 256 × 1,868 = 478 k positions → 11.5 MB int64 written. Reads 3 × 3.8 MB. One launch, no host read.
- The template comes from HF `get_rope_index` on one bare prompt row, **once per prompt layout** (once per run for
  LIDC-only training), not per batch.

### (b) Fused log-prob: `k3gpu::lse_logprob`

Grid = (⌈N/64⌉ token tiles) × (S vocab splits). A CTA has 256 threads (8 warps). For each 128-column vocab tile in its
split:

1. GEMM tile. Tensor-core mode: WMMA bf16 m16n16k16 with fp32 accumulators; warps 2 × 4 of 32 × 32; K-steps of 32
   through a 2-stage `cp.async` pipeline (zero-fill out of range). Exact mode: each thread runs 32 fma chains, k in
   increasing order, from the same staged tiles.
2. The 64 × 128 fp32 tile goes to shared memory (aliasing the A/B stages).
3. **Fixed-order epilogue.** Lane q of 4 per row owns columns 4i + q. It takes the max, then
   E = Σ exp(z − m_new) and X = Σ e·z in i order. Lanes combine as (l0 + l1) + (l2 + l3) through `__shfl_xor`, and
   the result folds into the running (m, s, x). The target logit is picked when its column is in the tile.

Only (m, s, x, target) per (split, token) leave the SM: a workspace of 16·S·N bytes. A finish kernel combines splits
in order s = 0..S−1 into logp = (z_t − M) − log S and entropy = M + log S − X/S. No (N, V) or (N, chunk) buffer
exists.

- FLOPs: 2·N·V·D = 15.3 TFLOP for N = 24,576 (16 rows × 1,536 tokens, the worst update micro-batch).
- Traffic: per CTA tile, 768 KiB of A+B staging feeds 33.5 MFLOP (42.7 FLOP/B) from L2. CTAs of one split are
  launched x-fastest, so they stream the same W tiles together and W (622 MB) is read from HBM about once per wave
  [estimate].
- Split count: S = min(64, ⌈4·SMs / token tiles⌉), so at least 4 CTAs per SM are in flight (`default_splits`).

### (b) IS ratio: `k3gpu::is_ratio`

This mirrors `fast.is_ratio` (TRL 1.17's four modes) element for element. diff = nan_to_num((logp − sampling) × mask),
with ±inf mapped to ±FLT_MAX as torch does. Token modes compute exp(diff) per token; sequence modes compute exp(Σ diff)
per row. `*_truncate` is clamp (token_truncate = **token-level TIS**, min(ρ, C) [paper: arXiv 2510.26788, eq. 9]);
`*_mask` zeroes ratios outside [lo, hi] (MIS). One block runs per row, plus a one-block stats kernel. The 7 stats
(d_sum, d_n, d_max, r_sum, r_n, r_min, r_max) are the ones `FastGRPO.set_is_ratio` keeps. All reductions are
fixed-shape trees. O(B·T), 2 launches, replacing about two dozen torch launches.

## 5. Memory and occupancy plan (sm_120: 128 KB smem/SM, ≤99 KB/block, 64 K regs/SM, ≤48 warps/SM, 170 SMs)

| Kernel | Threads | Regs | Smem | CTAs/SM (limit) | Notes |
|---|---|---|---|---|---|
| preprocess | 256 | 40 | 0 | 6 (warps) | grid-stride, ≤16 CTAs/SM launched; float4 / 8-byte stores, coalesced |
| rope | 256 | 31 | 1,184 B | 6 (warps) | one CTA per row (256 rows ≈ 1.5 waves) |
| lse tensor-core | 256 | 128 | 33,792 B | **2 (registers)** = 16 warps | smem would allow 3; 128 regs from 32 z values + 4 accumulator fragments; 0 spills |
| lse exact | 256 | 81 | 33,792 B | 3 (smem) | test and reference mode only |
| is_rows / is_stats | 256 | 33-39 | 1 KB | 6 | one CTA per row |

[measured: results/docker_nvcc_check2.log.txt (ptxas -v); limits: NVIDIA Blackwell Tuning Guide §1.4.1.1]. Every
allocation is static and below 48 KB, so no opt-in is needed and nothing approaches the 99 KB sm_120 limit. (CUTLASS
templates sized for sm_100 break at that limit [third-party wiki, section 11].)

## 6. Exactness argument

The rules are in `k3_common.cuh`: IEEE binary32 round-to-nearest everywhere; no fast-math; no FTZ; no implicit FMA
(`--fmad=false`, `-ffp-contract=off`; FMAs are written as `fmaf`); exp and log built only from correctly rounded
operations; fixed reduction orders that the CPU replays.

- **(a) fp32 equal to HF.** Each output is q = RN(RN(x − m′) / s′), with x = float(u8) exact and m′, s′ the fp32
  values HF computes. The torch wrapper gets them from HF's own `_fuse_mean_std_and_rescale_factor`. IEEE `-` and `/`
  are correctly rounded, so CPU torch (true division: verified bit-equal) and the GPU (`__fsub_rn`, `__fdiv_rn`) give
  the same q. The outputs lie in [−1, 1], with no subnormals, so FTZ settings cannot matter. The proof has two parts:
  (i) the value map is checked on its whole domain, 256 values × 3 channels; (ii) the index map does not depend on
  values. The value map is injective, so equality on the three **digit images** (pixel = byte k of the source
  element's flat index) proves the index map at every one of the 1,572,864 positions. Both parts pass against HF
  [measured: results/test_hf_cpu.log.txt]. bf16 = RNE of the same fp32 value, also equal.
- **(c)** Integer arithmetic only. The extrapolation formula equals HF's grouping for rows of the stated form, and
  the kernel checks the form (error bits).
- **(b) exact mode** is bit-equal to `k3ref::lse_ref` by construction: same fma chain, same lane/split order, same
  k3_expf/k3_logf. test_gpu checks it at the real head shape. **Tensor-core mode is not bit-exact, by design**: MMA
  accumulation order and rounding inside a k-block are not specified by NVIDIA, so no CPU can replay them. Its error
  against fp64 is measured (gate 2e-3 nats). That gate is 9× below the measured trainer-vLLM median per-token gap of
  0.018 [BRIEF §2].
- **(b) is_ratio**: the portable-exp variant is bit-equal to the CPU reference (test_gpu). The libdevice variant uses
  CUDA `expf`, which torch.exp on CUDA also uses (torch is built without fast-math). So the token modes should be
  bit-equal to `fast.is_ratio` (test_torch_gpu checks it). Sequence modes sum in another order than `torch.sum`, so
  expect last-bit differences.

## 7. Expected gains (estimates until GPU_RUN.md runs)

| Piece | Today | With K3 | Gain per step |
|---|---|---|---|
| (a) | 2 CPU processor passes + 1.61 GB pageable blocking H2D [measured: cpu_gpu_trace.md] | 25 MB pinned H2D + ≈1 ms kernel | the image part of 1.5-2.6 s of CPU time (two passes at 0.75-1.31 s each, tokenization included) [estimate; the HF-CPU vs GPU numbers in torch_gpu.json settle it] |
| (c) | host loop 30-50 ms + 8 syncs + 5 syncs per micro-batch [measured-Mac: cpu_gpu_trace.md] | one ≈0.05 ms launch, no sync | 30-50 ms of host loop and ≈90+ syncs (8 + 5 per forward chunk × 16 update chunks, plus the reference pass's) [estimate] |
| (b) lse | TRL tile buffer; fast cuBLAS | no buffer, deterministic; WMMA slower than cuBLAS is likely | memory; time decided by bench (rule in section 9) |
| (b) ratio | ~two dozen launches per micro-batch | 2 launches | <1 ms [estimate] |

## 8. Integration points (REPO file:function; nothing edited)

1. **(a)** `src/objdel/train/fast.py:ImageMemo.__call__` (lines 418-476). It keeps its memo per distinct image. It
   stacks the 32 images' uint8 arrays (`np.asarray(PIL)`; `dataset.OpenOnce.__getitems__`, dataset.py:275, already
   hands one object per image) into a pinned buffer, copies them asynchronously, and returns CUDA `pixel_values` from
   `preprocess(src, idx=row→image, bf16=True)`. TRL's `_prepare_inputs` (grpo_trainer.py:2643) is then a no-op move.
   Guard: use K3 only when `smart_resize(H, W) == (H, W)` and the image is RGB uint8; otherwise fall back to HF. At
   start-up, self-check one image both ways (bytes equal). Follow-up: `VisionCache.keys` (fast.py:262) hashes all of
   `pixel_values` every forward; with GPU patches it can key by image index.
2. **(c)** `fast.py:RopeCache.positions` (336-348) and `RopeCache.keys` (321-327). Build the template once from one
   row with HF `get_rope_index`; call `rope(...)` for all rows every time (no keys, no host reads). Fold the error word
   into the single per-step read in `FastGRPO.log` (604-616), with the HF path as fallback.
3. **(b) ratio** `fast.py:FastGRPO.set_is_ratio` (581-602) and `is_ratio` (351-362): one `is_ratio(...)` call
   returns the ratio and the 7 stats on the device. `log` reads them once (as now).
4. **(b) lse** `trl/trainer/utils.py:1693` (`_ChunkedLogProbFunction.apply` inside `add_fused_lm_head`). Wrap it so
   that when `torch.is_grad_enabled()` is False, and there is no `support_ids`/`score_ids`, softcap or bias,
   `logit_scale == 1`, and outputs ⊆ {entropy}, it calls K3 lse (the β = 1e-3 reference pass, evaluation). The grad
   pass keeps TRL's. A FastGRPO patch hook, like the existing `time_vllm_phases`, is the place.

## 9. Decision record: custom CUDA vs Triton vs torch ops vs a library

Context: the owner prefers a good library, even with a little customisation, over our own code, and wants GPU
everywhere.

| Piece | Option A (library / plain torch) | Option B (K3 CUDA) | Decision |
|---|---|---|---|
| (a) | HF processor with `device="cuda"` (accepted kwarg, bit-equal on CPU: checked); or NVIDIA DALI | one fused pass, row duplication and bf16 in the same pass, no PIL→tensor per image | **Measure both** (test_torch_gpu prints HF-on-CUDA bit-equality and time). Take A if it is bit-equal and within a few ms; B otherwise. |
| (c) | 5 torch ops (cumsum, gather, where), == HF (checked on CPU) | one launch + template validation | **A by default** (simplest; same exactness); B if the validation flags are wanted in production. |
| (b) lse | TRL `ChunkedLogProbFunction` (already fused, has backward); Apple CCE (Triton, backward) [paper: arXiv 2411.09009]; Liger [paper: arXiv 2410.10989] | forward only; no buffer; deterministic; CPU-replayable exact mode | **B for no-grad passes only if** bench shows it no slower than TRL; otherwise keep TRL. No custom backward: CCE/TRL already do it, and a deterministic backward without atomics costs a second GEMM pass. |
| (b) ratio | `fast.is_ratio` (torch) | fused ratio + stats | B only for the launch count; bit-equal in token modes. Low priority. |

Consequence: K3 is partly a correctness harness (CPU references, digit images, flags) that the library options can
also be checked against. Any K3 kernel that loses its row of this table gets archived, per the owner's rule.

## 10. Debugging plan (if the GPU disagrees)

1. Reproduce on the smallest case: one image (test_gpu prints counts; dump the first differing index).
2. (a): decode the index → (row, patch, feature) → (c, y, x) with `k3_patch_src`. Value differs but the digit images
   agree: the value map (look for fast-math/FTZ in the build line). Digit images differ: the index map.
3. (b) exact: compare per-split partials (m, s, x) with `lse_ref` (expose the workspace); the first split that
   differs isolates GEMM order (z) vs epilogue (E, X) vs finish. Then check `--fmad=false` reached nvcc (torch
   extension builds pass `extra_cuda_cflags`).
4. (c): err bits tell which check fired; positions differ with err == 0: compare the template with HF's on the bare
   prompt.
5. Prevention: the start-up self-checks in section 8 and the CPU tests in CI (`make test-cpu` needs no GPU).

## 11. Risks

1. **HF/torchvision drift** (fused mean/std, resize identity, patch order). Mitigation: mean′/std′ come from HF's
   function; a start-up bytes check; the digit-image test in CI.
2. **Other image sizes or modes** (shifted/rescaled train views must still render 512 × 512; RGBA/grayscale must be
   converted to RGB first). The kernel rejects bad shapes (`cudaErrorInvalidValue`); the wrapper must check
   `smart_resize`.
3. **Mixed prompts** (LIDC + RIGA in one batch): one template per launch today. Group rows per template or add a
   template table; the error word catches it.
4. **Tensor-core lse slower than cuBLAS.** WMMA with 64 × 128 tiles and 2 stages is not CUTLASS. The decision rule
   above applies; the larger-tile, 3-stage variant is the follow-up. GeForce fp32-accumulate tensor throughput is
   unverified (a third-party page says a reduced rate).
5. **torch.exp ≠ libdevice expf** on some build: token-mode bit-equality would fail by ≤1-2 ulp. Fallback: tolerance.
6. **Sequence-mode sums** differ from torch in the last bits; a ratio exactly at the cap could flip mask/keep (rare;
   it does not matter statistically).
7. **No backward** in K3 lse: only the no-grad passes can use it.
8. **Build environment**: torch 2.14 needs C++20 for extensions; the box image has no nvcc (install from PyPI, GPU_RUN
   step 1). test_gpu needs ~4 GB of host RAM for the 256-row reference.

Citations:
- NVIDIA Blackwell Tuning Guide (cc 12.0: 48 warps/SM, 128 KB smem/SM, 99 KB/block): https://docs.nvidia.com/cuda/blackwell-tuning-guide/
- CUDA Programming Guide, compute capabilities: https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/compute-capabilities.html
- F. Yao et al., "Your Efficient RL Framework Secretly Brings You Off-Policy RL Training" (TIS): https://fengyao.notion.site/off-policy-rl
- P. Qi, Z. Liu et al., Defeating the Training-Inference Mismatch via FP16 (token-TIS / seq-MIS definitions) [paper: arXiv 2510.26788]
- E. Wijmans et al., Cut Your Losses in Large-Vocabulary Language Models (CCE) [paper: arXiv 2411.09009]
- P.-L. Hsu et al., Liger Kernel [paper: arXiv 2410.10989]
- SkyRL off-policy correction docs (token vs sequence TIS): https://docs.skyrl.ai/docs/algorithms/off_policy_correction
- TRL issue #4493 (sequence-level TIS/MIS): https://github.com/huggingface/trl/issues/4493
- Third-party (not NVIDIA) sm_120 notes: 1,792 GB/s, 96 MB L2, fp32-accumulate rate: https://www.skillsdirectory.com/skills/kekzl-sm120-cuda-expert ; SM100 vs SM120 smem cliff: https://0xsero.github.io/blackwell-gpu-wiki/blackwell/sm100-vs-sm120/

## 12. Files

| File | What |
|---|---|
| `k3_common.cuh` | shared host/device formulas: bf16, k3_expf/k3_logf, k3_norm, k3_patch_src, k3_rope_pos, LSE fold/finish, IS diff/bound, tree ops |
| `k3_kernels.cu`, `k3_kernels.h` | the CUDA kernels and their stream-async C++ API |
| `k3_ref.cpp`, `k3_ref.h` | CPU reference: HF-structured preprocessing, C++ port of `get_rope_index`, exact LSE replay, fp64 LSE, IS replay |
| `test_cpu.cpp` | CPU tests (pass on macOS clang 21 and linux gcc) |
| `k3_cpu_cli.cpp`, `test_hf_cpu.py` | CPU reference vs Hugging Face (passes on the Mac) |
| `test_gpu.cu` | GPU vs CPU reference, adversarial inputs (`--quick` for compute-sanitizer) |
| `bench_gpu.cu` | micro-benchmarks at the R5 step's shapes |
| `k3_torch.cpp`, `test_torch_gpu.py` | torch binding; GPU vs HF/TRL/objdel + library alternatives + timings → JSON |
| `docker_nvcc_check.sh` | nvcc 13.0 sm_120 build check on the Mac (linux/amd64 container) |
| `Makefile`, `GPU_RUN.md` | build; the rented-GPU runbook |
| `build/*.log` | logs of the runs above |
