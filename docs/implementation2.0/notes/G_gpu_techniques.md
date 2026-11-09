# G: GPU techniques for the VLM RL loop (outside-world survey, 2026-11-08)

Scope: what exists outside this repo, what is proven, with numbers, mapped onto THIS loop (one RTX 5090, sm_120,
32 GB; Qwen3-VL-2B LoRA; TRL 1.17.1 colocate; vLLM; CUDA checker). No repo edits, no spend.
Tags: **[V]** checked in a primary source or in installed code; **[I]** my estimate; **[U]** unverified, measure first.
Baseline used for gains (R5, from `cuda_streams_graphs_review_2026-11-08.md` §0, numbers only): step 23.4 s =
generate 11.1 (weight sync 1.27 inside) + batch-minus-generate 4.65 + update 7.67. I did not inherit that note's
conclusions; I only avoid repeating its items (sync list, mb32, per-layer graphs, IPC, sleep levels).

## Skills used

| Skill | What it contributed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Pulled the live catalog (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`). Only the finder is installed. Matches for this task (recommend only, not installed; install needs owner OK): `dali-dynamic-mode` (Training AI), `nemo-mbridge-perf-cuda-graphs`, `nemo-mbridge-perf-sequence-packing`, `launch-nemo-rl` / `nemo-rl-docs`, `compute-sanitizer-diagnose-cuda-bugs`, `tilegym-*` (cuTile kernels), `warp-*`. No catalog skill for cuCIM, nvImageCodec, Nsight, CUTLASS or TensorRT, so those were researched from NVIDIA docs directly. Install form: `npx skills add nvidia/skills --skill <name> --agent claude-code --global --yes`. |
| huggingface-skills:huggingface-papers | HF papers search API found the two most relevant 2026 systems papers: **Rollplex** (2608.14498, cross-phase GPU sharing for VLM RL) and **System-Integrated Speculative Decoding in RL** (2604.26779), plus QaRL (2604.07853, FP8 rollout mismatch) and EfficientRollout (2606.18967). Also XGrammar metadata (2411.15100). |
| huggingface-skills:huggingface-best | Applied to the model-choice side effect on kernels: read `Qwen/Qwen3.5-2B` config (the main thread's leaderboard pick candidate). It is hybrid linear attention (3 of 4 layers `linear_attention`), head_dim 256, vocab 248,320. That changes several kernel conclusions below (§8). Leaderboard scores themselves: see `MAIN_skill_lookups.md`. |
| engineering:system-design | Step budget first (Amdahl): only phases with seconds in them get ranked; requirements = exact/bit-reproducible checker, GPU-first, one card, no spend; explicit trade-off table. |
| engineering:architecture | ADR-style decision for the polygonizer design (§1.4): options, dimensions, decision, consequences. |
| engineering:debug | "Low utilization" treated as a symptom: hypotheses -> the measurement that separates them (§6) before any fix. |
| engineering:code-review | Correctness lens on proposals: bit-exactness, determinism (atomic order), distribution mismatch (grammar masks, quantized rollouts), memory. Applied to the polygonizer header and the TRL vLLM path (read-only). |
| engineering:documentation | Reader-first layout: ranked table at the top of the findings, details after, every number with its URL. |

---

## 0. Ranked adoptable techniques for THIS loop

Gains are seconds per 23.4 s step **[I]** unless marked. Ordered by expected gain / (work x risk).

| # | Technique | Phase hit | Expected gain | Risk | Evidence |
|---|---|---|---|---|---|
| 1 | **Measure with NVTX + Nsight Systems timeline** (CUPTI trace, `--cuda-graph-trace=node`, `--capture-range=cudaProfilerApi`); GPU-metrics sampling only if the Vast host allows counters | all | 0 s, but picks 2-5 | low | §6 |
| 2 | **Rollplex-style overlap**: run the trainer's 32-image vision encode and the reference pass's prompt prefix (adapter off, no grad, exact) on a side stream inside vLLM's decode window, keep the prefix KV, then score only completion tokens | batch-minus-generate (4.65) | 0.5-1.5 s | medium-high (thread + GIL beside vLLM's loop, memory at 0.35 vLLM share, KV plumbing for packed rows) | Rollplex: 1.23-1.30x over serial colocation (32 H800, Qwen2.5-VL-32B) |
| 3 | **`detokenize=False` for vLLM** via `GRPOConfig.generation_kwargs` (TRL decodes from `completion_ids` itself) | generate (11.1) | 0-0.5 s [U] | low (incompatible only with `stop` strings) | vLLM `SamplingParams.detokenize`; PR #3749 ~20 % (old engine); PR #46422 1.8-2.5x when detok was per-step overhead |
| 4 | **Speculative decoding in rollouts** (vLLM `ngram`/suffix; no draft model exists below Qwen3-VL-2B) | generate | 0-2 s [U]; lossless | medium (acceptance on coordinate digits unknown; gains mostly in the shrinking-batch tail) | NeMo-RL + vLLM: 1.8x rollout throughput at 8B sync RL; EfficientRollout: SD pays only when batch is memory-bound |
| 5 | **GPU-resident uint8 images + transformers torchvision backend `device="cuda"`** for the trainer's processor pass | batch-minus-generate | 0.1-0.3 s | medium: must stay bit-identical to the pixels vLLM saw (vLLM preprocesses on CPU) | HF issue #50100: 30 x 512^2 Qwen3-VL images 2.0-2.7 ms GPU vs 44-74 ms CPU |
| 6 | **vLLM sleep + `expandable_segments:True` together** (now allowed) to cut trainer fragmentation; may unlock the bigger micro-batch | update | 0-1 s (indirect) [U] | low-medium | vLLM PR #40812 / `cumem.py` toggles expandable segments off only inside the sleep pool |
| 7 | **FA4 CuTe-DSL sm_120 varlen** (fwd+bwd) as the packed-row attention, if the torch `varlen_attn` path disappoints | update + ref | 0.3-1 s [U] | medium-high (out-of-tree PR/bundle; no deterministic bwd) | Dao-AILab PR #2634; HF `SecondNatureComputing/flash-attn-4-sm120` |
| 8 | **Grammar-constrained decoding** (the repo regex already fits xgrammar) | quality, not speed | -7 to -17 % ITL | high for RL math: masked sampler vs full-vocab trainer logprobs; needs sampling-mask replay | §3 |
| 9 | Liger fused-linear GRPO loss | update memory | ~0 s (TRL already chunks the head) | low | Liger blog: up to 40 % memory, no speed claim |
| 10 | FP8 / MXFP8 / NVFP4 training on sm_120 | update | ~0-0.5 s | high | TE support matrix: MXFP8 **No** on 12.0, NVFP4 partial, FP8 per-tensor yes |
| 11 | Green contexts / MPS / more priority streams | checker | 0 s | low | checker already waits ~5 ms; MPS needs multi-process |
| 12 | GPU token->coordinate parser, GPU polygon simplification, new CCL (BKE) | reward | 0 s | medium (bit-exact parity work) | §1, §3 |

Bottom line: the CUDA geometry side (checker, scorer, polygonizer) is already at the design frontier and is not on
the critical path (~5 ms wait per step). The seconds are in generation (11.1 s) and in serial prefix work that
could hide under decode (Rollplex pattern). Measure (#1) before #2-#5.

---

## 1. GPU polygonization, contour tracing, connected components

### 1.1 Libraries [V]

| Tool | What it gives | Fit for an exact crack-edge polygonizer |
|---|---|---|
| **NPP** `nppiLabelMarkersUF` + `nppiCompressedMarkerLabelsUFContoursGenerateGeometryLists` (CUDA 11.x-13.3) | Union-find CCL (Chen et al. "Optimized UF ... using GPUs"); ordered contour pixel lists, CW or CCW | **No.** Pixel-centre chains, not crack edges; lists are returned in **host** memory; the doc itself says "a few contour pixels can escape insertion into the ordered geometry list" and the bounding box "occasionally" misses pixels. Not exact. https://docs.nvidia.com/cuda/archive/11.4.4/npp/group__image__filter__compresssed__marker__labels__info.html ; label: https://docs.nvidia.com/cuda/archive/11.4.3/npp/group__image__filter__label__markers.html ; structs in NPP 13.3: https://docs.nvidia.com/cuda/archive/13.3.0/npp/nppdefs.html |
| **cuCIM** 26.08 (`cucim.skimage.measure.label`, 64-bit indexing added in #1087) | GPU CCL (CuPy-derived kernels, grayscale mode) | Labels only. **No GPU `find_contours`** (marching squares not ported; hipCIM's mirror of the support list says advanced measure functions like `marching_cubes` lack GPU acceleration). Marching squares also gives half-pixel interpolated, unordered-then-stitched contours, not crack edges. https://github.com/rapidsai/cucim/releases ; https://rocm.docs.amd.com/projects/hipCIM/en/latest/reference/supported-functionality.html ; skimage semantics: https://scikit-image.org/docs/stable/api/skimage.measure.html |
| **OpenCV CUDA** | `cv::cuda::connectedComponents` exists; `findContours` is CPU only (Suzuki) | No. |
| **YACCLAB** (Allegretti, Bolelli, Grana, IEEE TPDS 2019) | BUF and BKE block-based CCL, 2D+3D | BKE-IC 1.1x (3DPeS) to 1.4x (XDOCS) faster than KE in 2D; 1.4-2.5x in 3D vs BE; BUF beats BKE below ~5-10 % density. Measured on a Maxwell K2200 with data already on device. https://federicobolelli.it/media/publications/pdfs/2019tpds.pdf ; https://github.com/prittt/YACCLAB |

### 1.2 Papers on GPU contour tracing [V]

- **G-Contour** (ICCAD 2025, CUHK): boundary extraction -> GPU union-find (atomicMin, path compression, shared-memory
  tiles) -> start-point (inflection) selection -> parallel per-contour search. 10x average, 98x peak vs OpenCV
  on VLSI layouts. http://www.cse.cuhk.edu.hk/~byu/papers/C295-ICCAD2025-GContour.pdf
- **Parallel Suzuki** (Garcia-Molla et al., J. Supercomputing 2022): tile the image, trace per tile, stitch;
  >10x vs fastest CPU. https://link.springer.com/article/10.1007/s11227-021-04260-y
- **Triads** (Alonso-Jorda et al., J. Supercomputing, 2026-04): 3 kernels (triad extraction, connection, contour id);
  goal is "CPU-level speed without transfers", i.e. no big speedup claimed.
  https://link.springer.com/article/10.1007/s11227-026-08488-4
- **Crack codes on GPU** (Butt et al. 2015): ~25,000 contours of a 2048x768 map in ~29 ms on a 480-ALU GPU, ~40x
  OpenCV CPU. https://ph02.tci-thaijo.org/index.php/stej/article/view/251345

### 1.3 List ranking / pointer jumping [V]

- Wyllie pointer jumping is O(n log n) work; random-splitter / Helman-JaJa is O(n). Rehman, Kothapalli, Narayanan
  (ICS 2009): RHJ ranks a random 32 M-element list in ~1 s, 8-9x CPU; Wyllie "does not scale well to large lists".
  https://doi.org/10.1145/1542275.1542311 ; Yogaratnam (arXiv 1002.4482): splitter beats Wyllie, GPU ~20x sequential.
  https://arxiv.org/html/1002.4482
- For this loop: one ring is a few hundred to ~2k unit edges (512x512 slices, per-object bounding boxes), so Wyllie
  needs ~9-11 rounds and fits in shared memory; the O(log n) extra work is irrelevant at that size [I]. Splitters
  only pay off at 10^6+ elements.

### 1.4 ADR: polygonizer design (engineering:architecture)

**Context.** Requirement: exact, bit-reproducible crack-edge rings, same ring order and start vertex as the CPU
reference, GPU-only, masks never leave the device. Current design (header of `cuda/gpucheck.cpp`, read-only):
every unit boundary edge is a slot; successor from the 2x2 bit-quad (left turn at saddles = 8-connected foreground,
4-connected background); ring label = min slot on the cycle by pointer jumping; corner index by weighted list
ranking; rings ordered by first east edge in raster order.

| Option | Exact | Deterministic | Work | Verdict |
|---|---|---|---|---|
| A. Current: local successor + min-label pointer jumping + weighted list ranking | yes (integer, local rule) | yes: min is commutative, no order-dependent atomics | O(n log n), n = edges per object | **keep** |
| B. NPP contour geometry lists | no ("pixels can escape") | n/a | - | reject |
| C. Marching squares (skimage semantics) then stitch | half-pixel coords, open contours at borders | stitching order needs care | O(n) + sort | reject |
| D. G-Contour style: CCL + start points + one thread walks each ring | yes if starts are canonical | yes | O(n) but serial depth = ring length | equal speed at these sizes; no reason to switch |
| E. Tile-and-stitch Suzuki / triads | pixel chains, not crack edges | stitching order | O(n) | reject |

**Decision.** Keep A. The only listed improvement worth noting: if objects ever reach 10^5+ edges (whole-slice
masks, 3-D surfaces), switch the ranking to a ruling-set (Helman-JaJa) pass; nothing else in the literature is
both exact and faster. Code-review notes for any change: keep every atomic either commutative (min/add on
integers) or absent; never use float atomics in the scorer (order-dependent rounding breaks bit parity).

### 1.5 Simplification and rasterization [V]

- GPU Ramer-Douglas-Peucker: iterative fine-grained form (no recursion, one thread per point, segmented reductions):
  3.8-6.8x multi-core CPU, ~1000 fps on embedded GPUs (Cebrian and Moure 2020).
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7302271/ ; trajectory DP with scan/segmented scan 10-15x serial:
  https://doi.org/10.48550/arxiv.2004.13653
- For crack-edge rings the exact simplification is collinear-run removal (keep turns only), which the current
  "corner index" already is. RDP would make the output depend on a float tolerance; if ever needed, use integer
  cross products and a fixed tie-break to stay bit-reproducible [I].
- Rasterization: the repo's warp-per-row winding + `__ballot_sync` is the standard scanline approach; no library
  (NPP, cuCIM, OpenCV CUDA) offers an exact even-odd polygon raster with integer edge tests. Nothing to adopt.

---

## 2. GPU image pipeline

Facts [V]:
- **PNG has no NVIDIA GPU decoder.** nvImageCodec lists `png: BackendKind.CPU_ONLY: opencv_png_decoder`
  (https://docs.nvidia.com/cuda/nvimagecodec/samples/query_and_select_plugins.html). Latest: nvImageCodec 0.8.0
  beta (pydicom GPU plugin, HTJ2K transcoding; https://docs.nvidia.com/cuda/nvimagecodec/release_notes_v0.8.0.html),
  0.9.0 now consumed by DALI 2.3.0 (2026-08-28, https://github.com/NVIDIA/DALI/releases/tag/v2.3.0).
  `torchvision.io.decode_jpeg(device="cuda")` is JPEG only.
- **DALI** 2.2/2.3: dynamic (eager) mode with "transparent pipelining"/capture mode, nvImageCodec default decoder,
  DLPack upgraded, a torchvision-style API (crop, RandomResizedCrop). https://github.com/NVIDIA/DALI/releases/tag/v2.2.0
- **cuCIM** 26.06/26.08 uses nvImageCodec 0.8 for TIFF/JPEG2000/JPEG; nothing for PNG.
- **Transformers image processors**: the torchvision backend is the default, and "you can set the `device`
  argument" (`processor(images, return_tensors="pt", device="cuda")`).
  https://github.com/huggingface/transformers/blob/main/docs/source/en/main_classes/image_processor.md .
  Qwen-VL fast processors: "up to 3x on CPU and 26x on GPU" (Gozlan, HF, 2025-07). Measured in
  https://github.com/huggingface/transformers/issues/50100 (Qwen3-VL processor, 30 x 512x512 per call):
  CPU 37.8-51.8 ms min (v5.19 regression from in-place normalize: +20 %, page faults), GPU 2.72 ms -> 2.03 ms
  with a proposed preallocated patch buffer. vLLM's own `apply_hf_processor` for Qwen3-VL, 30 images:
  84.9 ms (5.17) -> 104.4 ms (5.19).

What this means here [I]:
- The images are 512x512 PNG (146 KB, 2.6 ms CPU decode each, repo note). Best GPU-first design is not a GPU
  decoder but **no decode at all**: decode the train set once to uint8 and keep it resident on the GPU
  (e.g. 30k crops x 512x512x1 B = ~7.9 GB is too much beside vLLM; per-epoch shards or a 1-channel crop store
  of the views in use). DLPack (`torch.from_dlpack`) moves them between CuPy/DALI/torch without copies.
- The trainer-side processor can run on the GPU, but **vLLM processes the same images on the CPU** in colocate
  (TRL passes PIL images in `multi_modal_data`, `trl/generation/vllm_generation.py` [V installed code]). The repo's
  invariant "vLLM samples from exactly the image features the trainer scores" then needs bit-equal pixels from
  two code paths. With Qwen3-VL (patch 16, merge 2, factor 32) a 512x512 image needs no resize, so only
  rescale+normalize differs; GPU vs CPU float rounding (FMA) may still differ in the last ulp. Test equality
  before adopting. Feeding vLLM precomputed image embeddings instead (vLLM `image_embeds` inputs) would remove
  vLLM's CPU preprocessing and its vision encode, but Qwen3-VL's DeepStack features make that path [U].
- Gain is small (0.1-0.3 s of the ~1.5-2 s GPU-idle CPU window). Worth it mainly under the GPU-first rule.

---

## 3. Text <-> geometry on the GPU; constrained decoding

### 3.1 Parsing token ids on the GPU [V tokenizer facts]

- Qwen3-VL's pre-tokenizer regex has a single `\p{N}` alternative: **every digit is its own token**. In
  `Qwen/Qwen3-VL-2B-Instruct/tokenizer.json`: digits '0'..'9' = ids 15..24, and **0 multi-digit tokens** in the
  151,643-token vocab (checked locally). Punctuation is merged: '[[' 15505, '],[' 14955, '],' 1125, ']]' 5053,
  '{"' 4913, '":' 788.
- So a GPU parser is a lookup table (token id -> digit value, or a short punctuation code) plus a segmented scan:
  number value by Horner's rule as an associative (multiply, add) pair, numbers paired into points, bracket depth
  by prefix sums, validity by per-row reductions. O(tokens), one kernel or a few torch ops.
- Value for THIS loop: ~0 s. The repo's parser is a linear regex + C `json.loads`, bounded at 16,384 characters
  (`src/objdel/geometry/text.py` header); 256 completions parse in milliseconds on a background thread. A GPU
  parser must reproduce every rejection code of the regex grammar bit-for-bit (differential test on all logged
  completions). Only do it if the owner wants the end-to-end "ids never leave the GPU" story.
- Cheaper win on the same theme: **`detokenize=False`** in vLLM (item #3). TRL calls
  `processing_class.batch_decode(completion_ids, ...)` itself (`grpo_trainer.py:2338-2341`, installed 1.17.1) and
  forwards `generation_kwargs` into `SamplingParams` (`vllm_generation.py:675-697`), so vLLM's own incremental
  detokenization is pure overhead. https://docs.vllm.ai/en/latest/api/vllm/sampling_params/ (field `detokenize`;
  "only exposed at the engine level"), https://github.com/vllm-project/vllm/pull/3749 (~20 % on A100 then),
  https://github.com/vllm-project/vllm/pull/46422 (1.80x-2.53x where per-step detok dominated, RTX PRO 6000).

### 3.2 Grammar-constrained decoding [V]

- Backends: vLLM `auto` prefers **xgrammar**, falls back to **guidance (llguidance)**, outlines, lm-format-enforcer.
  XGrammar: <40 us/token for JSON schema and CFG, up to 100x faster than prior engines (MLSys 2025,
  https://arxiv.org/pdf/2411.15100v3.pdf). llguidance: ~50 us CPU per mask for a 128k vocab, <1 % of masks >1 ms
  (https://github.com/guidance-ai/llguidance/blob/main/README.md).
- Measured cost in vLLM 0.31/0.32 (Ada laptop GPU, Qwen3-0.6B family): ITL +7 % at concurrency 1, **+17 % at
  concurrency 16**; regex grammars compile in ~15 ms (https://ttlab.hopto.org/blog/vllm-features-04-structured-outputs).
  The forum notes xgrammar's CPU-GPU sync overhead was worse on an RTX 4090 than on datacenter cards
  (https://discuss.vllm.ai/t/general-questions-on-structured-output-backend/1444). At 256 concurrent rows the
  per-step mask work grows with batch [I].
- The repo's grammar is already xgrammar-ready (no lookarounds; `GRPOConfig.vllm_structured_outputs_regex`).
- **RL correctness risk (code-review lens):** a grammar mask makes the sampler draw from a renormalized subset,
  while the trainer computes full-vocabulary log-probs. This is the same bias TRL has for top-p/top-k
  (https://github.com/huggingface/trl/issues/6789): with `logprobs_mode="processed_logprobs"` (set in TRL's
  colocate `LLM(...)`), vLLM's sampled log-prob exceeds the trainer's by -log(kept mass), which contaminates the
  vLLM importance-sampling correction. Fix pattern = **Keep Sampling Mask / distribution replay**: vLLM
  `return_sampling_mask` returns the kept token set per step and the trainer masks its logits the same way
  (https://docs.vllm.ai/en/latest/training/sampling_mask/; DeepSeek-V3.2 §3.3; prime-rl #3235, slime v0.3.1).
  TRL 1.17.1 already passes `return_sampling_mask` when configured [V installed code]. Whether vLLM's sampling mask
  captures the grammar bitmask (applied before the sampler) is [U]; test on one batch. Also: grammar-constrained
  decoding distorts the model's distribution (Grammar-Aligned Decoding, https://arxiv.org/html/2405.21047v2).
- Recommendation: use the grammar only if the logged format-failure rate is material; it is a quality tool, it
  costs decode time, and it needs the mask replay to keep the RL math exact.

---

## 4. Training-side kernels

- **Attention on sm_120** [V]: sm_120 keeps Ampere-style `mma.sync` (no WGMMA, no tcgen05/TMEM), 99 KB SMEM.
  FA2/FA3 cubins do not cover it; vLLM's default sm_120 attention backend is FlashInfer
  (https://github.com/vllm-project/vllm/pull/40110). FA4 CuTe-DSL sm_120 work: fwd (#2329), bwd (#2330), varlen
  (#2333) are in `fa4-v4.0.0.beta5` (https://github.com/Dao-AILab/flash-attention/releases/tag/fa4-v4.0.0.beta5);
  PR #2634 adds public-API routing: fwd+bwd dense/causal/varlen/GQA, bwd head dims 64/96/128/256,
  `torch.compile`-opaque and bit-identical to eager, **`deterministic=True` backward unsupported**, causal MHA at
  seqlen >= 8192 ~0.95x FA2, GQA "at parity or faster" (https://github.com/Dao-AILab/flash-attention/pull/2634).
  Bundle: https://huggingface.co/SecondNatureComputing/flash-attn-4-sm120. Qwen3-VL-2B uses GQA, head_dim 128 -> in
  the supported set. Compare against torch `varlen_attn` (already researched in the repo) with a 30-s benchmark.
- **Chunked/fused log-prob** [V]: Liger `LigerFusedLinearGRPOLoss` via `use_liger_loss=True`; up to 40 % less
  peak memory on Qwen3-0.6B, no speed claim (https://huggingface.co/blog/liger-grpo); it chunks on batch and needs
  padded [B, L, H] (https://github.com/Scicom-AI-Enterprise-Organization/small-ablation/tree/main/multipacking-grpo),
  which conflicts with the repo's shared-prefix packed rows. Liger PR #1160 streams over the vocab. TRL already
  avoids full logits, so expected speed gain ~0.
- **FP8 / MXFP8 / NVFP4 on sm_120** [V]: Transformer Engine support matrix: 12.0 = BF16 yes, FP8 per-tensor yes,
  FP8 block scaling yes (CUDA >= 12.9), **MXFP8 no**, NVFP4 partial (forward only, or training with RHT and
  stochastic rounding off); default recipe on 12.0 is `Float8CurrentScaling`
  (https://nvidia.github.io/TransformerEngine/support_matrix.html). TE 2.14 enabled sm_120 fused attention with
  cuDNN >= 9.18.1; TE 2.20 added NVFP4 SR/RHT split-quantization on SM120/121. Third-party SM120 NVFP4 polyfills:
  https://github.com/Harry-Chen/fp4_sm120. For a 2B LoRA run whose update is dominated by fixed per-micro-batch
  cost (repo measurement), FP8 GEMMs buy little and add a train/rollout mismatch. Skip.
- **FP8 rollouts** [V]: QaRL (2604.07853) shows quantized rollouts amplify the training-inference gap and produce
  "error tokens" in long outputs; needs a rollout-aligned training forward. Not worth it at 2B.
- **CUDA graphs and memory** [V]: vLLM sleep mode used to assert against `expandable_segments:True`; vLLM PR #40812
  (merged; `vllm/device_allocator/cumem.py` on main) now disables it only inside the sleep pool and restores it
  after, so the trainer can use expandable segments (less fragmentation) while vLLM sleeps
  (https://github.com/vllm-project/vllm/pull/40812, root cause https://github.com/pytorch/pytorch/issues/147851).
  Check the box's vLLM version contains it. TE 2.20 also reduced CUDA-graph memory retention in
  `make_graphed_callables` (#3427) and added `capture_time_hooks` (#2831), useful if per-layer graphs are tried.

---

## 5. Async overlap on one GPU

- **Rollplex** (2608.14498, HF papers): VLM RL runtimes run rollout, reference scoring and training in strict
  serial phases; the prefix (vision + prompt) part of the reference/training forward does not depend on the
  generated response, so it runs during rollout decode, which leaves compute idle. 1.23-1.30x over serial
  colocation, 1.57-2.24x over disaggregation (32 H800, Qwen2.5-VL-32B), synchronous on-policy semantics kept.
  Mapping here [I]: the exact, cheap version is the **reference pass prefix** (adapter off, no grad): compute its
  KV during decode, then run only the completion branches. The update pass prefix needs gradients through the
  LoRA weights, so only its vision features (frozen tower) can move. The 32-image vision encode moves for both.
- **Speculative decoding in RL** (2604.26779): NeMo-RL + vLLM, lossless, 1.8x rollout throughput at 8B synchronous;
  projected 2.5x end-to-end at 235B with async. EfficientRollout (2606.18967): SD helps only when the active batch
  has shrunk into the memory-bound regime; uses a quantized self-drafter and a toggle. For Qwen3-VL-2B there is no
  MTP head or smaller draft; vLLM `ngram`/suffix speculation is the only zero-training option [U acceptance].
- **Green contexts** [V]: driver API since CUDA 12.4, runtime API (`cudaGreenCtxCreate`,
  `cudaExecutionCtxStreamCreate`) since CUDA 13.1 (https://docs.nvidia.com/cuda/cuda-programming-guide/04-special-topics/green-contexts.html;
  NVIDIA blog 2026-10-06). Minimum partition on 9.0+: 8 SMs, multiple of 8. PyTorch: `torch.cuda.green_contexts`
  (experimental since the 2025-10 commit, CUDA >= 12.8; 2.15 docs add `SMPartition.split`, driver 13.1+):
  https://docs.pytorch.org/docs/2.15/generated/torch.cuda.green_contexts.SMPartition.html . Use here: bound the
  interference of Rollplex-style side work with vLLM decode (e.g. 16 of 170 SMs). Not needed for the checker.
- **Priority streams**: already used by the checker (lowest priority, non-blocking). Nothing to add.
- **MPS** [V]: Linux only, per-process clients; all accounting goes to the MPS server. Users run it on GeForce
  (Folding@home: 3090 9M -> 14M PPD), but official GeForce support is unclear. The loop is one process (TRL
  colocate, `distributed_executor_backend="external_launcher"`), so MPS does not apply.
  https://docs.nvidia.com/deploy/mps/latest/when-to-use-mps.html
- **Weight sync via CUDA IPC** [V]: vLLM has an `ipc` weight-transfer backend for same-GPU trainer/server, but in
  colocate the trainer already calls `load_weights` in-process (no copy across processes), so IPC only matters if
  the loop moves to server mode. Existing repo note covers it.

---

## 6. Measurement (engineering:debug lens)

Symptom: "GPU busy ~70 %" (repo profiler). Hypotheses and the measurement that separates them:

| Hypothesis | Measurement |
|---|---|
| H1 GPU idle on CPU work (processor, Python, syncs) | nsys timeline: gaps on the GPU row aligned with NVTX ranges on the CPU row; `cudaStreamSynchronize`/`cudaMemcpy` stalls |
| H2 GPU "busy" but under-filled (small kernels, low occupancy) | nsys GPU metrics: SMs Active, SM Issue, Tensor Active, warp occupancy per sample; ncu on the top 5 kernels |
| H3 decode tail: few long completions keep the step open | NVTX per vLLM generate + per-step running-sequence count (vLLM scheduler stats); completion length histogram |
| H4 memory churn (allocator, wake/sleep) | `torch.cuda.memory._record_memory_history` snapshot; nsys CUDA API row for `cuMemMap`/`cudaMalloc` |

Facts [V]:
- **"Utilization" is two things.** NVML/`nvidia-smi` utilization = fraction of time at least one kernel was running
  (busy %). Nsight Systems "SM Active" = fraction of SMs with work per sample; warp occupancy and Tensor Active go
  deeper (https://developer.nvidia.com/blog/measuring-the-gpu-occupancy-of-multi-stream-workloads/). A 100 % busy
  GPU can have 10 % SMs active. "Max utilization" for this loop = minimum step time; report SM Active x time and
  Tensor Active, not busy %.
- Nsight Systems samples GPU metrics up to 200 kHz (`--gpu-metrics-devices`, `--gpu-metrics-set`, 10 Hz-200 kHz);
  GeForce has its own metric sets (https://docs.nvidia.com/nsight-systems/UserGuide/index.html;
  https://forums.developer.nvidia.com/t/time-series-data-of-metrics/260078/2).
- **Counter permission risk on Vast:** GPU metrics sampling and Nsight Compute need performance-counter access
  (ERR_NVGPUCTRPERM). Inside a container, access must be enabled on the host or the container started with
  `--cap-add=SYS_ADMIN` (R565+: `CAP_PERFMON`; R610+: nvidia-capabilities `profiler-device`)
  (https://developer.nvidia.com/nvidia-development-tools-solutions-ERR_NVGPUCTRPERM-permission-issue-performance-counters).
  A rented container usually cannot change this. CUPTI **activity tracing** (kernel/memcpy timeline, NVTX) does
  not need counters, so the nsys timeline still works. Plan the preflight to test `ncu --version` + one tiny
  `ncu` run and fall back to timeline-only.
- NVTX: `torch.cuda.nvtx.range_push/pop` (or `torch.autograd.profiler.emit_nvtx`), and in C++ NVTX3 header-only
  `nvtx3::scoped_range` around each checker call; `nsys profile --capture-range=cudaProfilerApi` with
  `torch.cuda.profiler.start/stop` around steps 3-5 keeps traces small; `--cuda-graph-trace=node` shows vLLM's
  graphed decode kernels individually.

---

## 7. What I checked in installed code (read-only)

- TRL 1.17.1 `trl/generation/vllm_generation.py`: colocate `LLM(...)` sets `logprobs_mode="processed_logprobs"`,
  `max_num_batched_tokens=4096`, optional `return_sampling_mask`; colocate sends `n=1` per row (the 8 answers of a
  prompt are 8 requests); `generation_kwargs` update `SamplingParams`; `structured_outputs_regex` supported.
- `trl/trainer/grpo_trainer.py:2338-2341`: completions text comes from `batch_decode(completion_ids)`, not from
  vLLM's text.
- Repo configs: temperature 1.0, top_p 1.0, top_k 0 (`src/objdel/train/configs.py:187`), so there is no top-p
  truncation bias today; a grammar would introduce one.

## 8. Model-switch warning (from huggingface-best)

`Qwen/Qwen3.5-2B/config.json` [V]: `layer_types` = 3 `linear_attention` : 1 `full_attention`
(`full_attention_interval` 4), 24 layers, hidden 2048, head_dim 256, 2 KV heads, vocab 248,320. If the owner moves
to Qwen3.5: (a) shared-prefix packing by attention mask does not carry over to the linear-attention layers
(recurrent state needs varlen/`cu_seqlens` kernels, e.g. flash-linear-attention, and branch-from-prefix state
copies); (b) head_dim 256 backward on sm_120 exists only in the FA4 PR path (D256 listed), not FA2; (c) logits are
1.64x wider (248k vs 151k vocab), so the chunked head matters more; (d) the GPU token parser needs a fresh check of
digit tokenization (Qwen3.5's regex still has a single `\p{N}` per the pre-tokenizer survey at
https://blog.qsliu.dev/post/pretokenizer/ [V secondary]).

## Sources (primary first)

NPP: https://docs.nvidia.com/cuda/archive/11.4.4/npp/group__image__filter__compresssed__marker__labels__info.html ·
https://docs.nvidia.com/cuda/archive/13.3.0/npp/nppdefs.html · cuCIM: https://github.com/rapidsai/cucim/releases ·
https://github.com/rapidsai/cucim/pull/1087 · YACCLAB/BKE: https://federicobolelli.it/media/publications/pdfs/2019tpds.pdf ·
G-Contour: http://www.cse.cuhk.edu.hk/~byu/papers/C295-ICCAD2025-GContour.pdf · Suzuki GPU:
https://link.springer.com/article/10.1007/s11227-021-04260-y · triads: https://link.springer.com/article/10.1007/s11227-026-08488-4 ·
list ranking: https://doi.org/10.1145/1542275.1542311 , https://arxiv.org/html/1002.4482 · RDP GPU:
https://pmc.ncbi.nlm.nih.gov/articles/PMC7302271/ · nvImageCodec: https://docs.nvidia.com/cuda/nvimagecodec/release_notes_v0.8.0.html ,
https://docs.nvidia.com/cuda/nvimagecodec/samples/query_and_select_plugins.html · DALI: https://github.com/NVIDIA/DALI/releases/tag/v2.3.0 ·
HF processors: https://github.com/huggingface/transformers/issues/50100 , https://github.com/huggingface/transformers/pull/39591 ·
Qwen tokenizer: https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct (tokenizer.json, checked locally) · vLLM sampling params:
https://docs.vllm.ai/en/latest/api/vllm/sampling_params/ · vLLM sampling mask: https://docs.vllm.ai/en/latest/training/sampling_mask/ ·
TRL top-p bias: https://github.com/huggingface/trl/issues/6789 · XGrammar: https://arxiv.org/pdf/2411.15100v3.pdf · llguidance:
https://github.com/guidance-ai/llguidance · structured-output cost: https://ttlab.hopto.org/blog/vllm-features-04-structured-outputs ·
FA4 sm_120: https://github.com/Dao-AILab/flash-attention/pull/2634 , https://github.com/vllm-project/vllm/pull/40110 ·
TE: https://nvidia.github.io/TransformerEngine/support_matrix.html , https://docs.nvidia.com/deeplearning/transformer-engine/release_notes.html ·
Liger: https://huggingface.co/blog/liger-grpo · vLLM expandable segments: https://github.com/vllm-project/vllm/pull/40812 ·
green contexts: https://docs.nvidia.com/cuda/cuda-programming-guide/04-special-topics/green-contexts.html ,
https://docs.pytorch.org/docs/2.15/generated/torch.cuda.green_contexts.SMPartition.html · MPS:
https://docs.nvidia.com/deploy/mps/latest/when-to-use-mps.html · Nsight: https://docs.nvidia.com/nsight-systems/UserGuide/index.html ,
https://developer.nvidia.com/blog/measuring-the-gpu-occupancy-of-multi-stream-workloads/ , ERR_NVGPUCTRPERM page above ·
papers: https://huggingface.co/papers/2608.14498 (Rollplex), https://huggingface.co/papers/2604.26779 (SD in RL),
https://huggingface.co/papers/2606.18967 (EfficientRollout), https://huggingface.co/papers/2604.07853 (QaRL).
