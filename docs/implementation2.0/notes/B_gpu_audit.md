# B. GPU audit of the GRPO hot path: transfers, syncs, idle GPU (2026-11-08)

Read-only audit. No repo file changed, no git change, no GPU rented. Python that touched `data/` ran with `-I`.
Code read: `src/objdel/train/{fast,rewards,serve,dataset,configs,callbacks,sample}.py`, `src/objdel/reward/*`,
`src/objdel/geometry/*`, `src/objdel/checker/*`, `src/objdel/targets/build.py`, `cuda/gpucheck.{h,cpp}`,
`cuda/gpucheck_py.cpp`, `experiments/phase2/{train_rl,g_probe}.py`, installed TRL 1.17.1 / transformers 5.20.1 /
peft 0.21.2 / torch 2.14.1 in `.venv` (vLLM is not installed locally). Measurements: `runs/r5/{base,shared}_rl.json`
(R5, RTX 5090, 15 steps each), `runs/rl/train_rl/rl.json` (770-step full run), `experiments/cuda_checker/RESULTS.md`.
The notes in `docs/research/notes/` were treated as claims. Where I re-checked one, I say so. `cpu_gpu_trace.md` is
stale: it describes the code before GPU scoring. Its headline "the scoring is on the CPU" is no longer true
(`score.py:115` now calls `_gpucheck.score_polygons`).

Tags: **[M]** = measured number from a repo log, recomputed by me. **[V]** = verified by reading code. **[E]** =
estimate (arithmetic or judgment, not measured). **[U]** = unverified (vLLM internals; vLLM is not installed locally).

---

## Skills used

| Skill | What it contributed |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | Live catalog check (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`). Relevant entries, none installed, recommended only: `nemo-mbridge-perf-cuda-graphs` (I read its SKILL.md: TE-scoped per-layer graphs, warmup steps, RNG tracker; it shaped §4.4), `nemo-mbridge-perf-sequence-packing`, `compute-sanitizer-diagnose-cuda-bugs` (I read its SKILL.md: memcheck / racecheck / initcheck / synccheck routing; used for the polygonizer test plan in §3.2), `tilegym-monkey-patch-kernels-to-transformers`, `launch-nemo-rl` / `nemo-rl-docs` (disaggregated RL reference). There is no Nsight Systems skill in the catalog, so §4.7 follows plain nsys practice. |
| `engineering:code-review` | Review lens (correctness, performance, maintainability). It structured §5, "Bugs and risks", with severity per item. |
| `engineering:debug` | Reproduce, isolate, diagnose. I rebuilt the step timeline from the R5 logs and isolated the decode-tail cost with a regression (§1). |
| `engineering:system-design` | Requirements, then high-level design, then trade-offs. Used for §4 (one GPU and many GPUs), with each trade-off stated. |
| `engineering:architecture` | ADR form. §4.1 is a short ADR: synchronous colocate vs asynchronous generation. |
| `engineering:documentation` | Most useful information first: the TL;DR, then sections with file:line citations. Links instead of copies. |
| `huggingface-skills:huggingface-best` | Used through the main thread's lookups (`MAIN_skill_lookups.md`): Qwen3.5-0.8B/2B/4B (Apache-2.0) are candidates. I note a portability risk in §4.8: TRL ships a `tree_gated_delta_rule.py` for packed forwards (`trl/experimental/async_grpo/`), so the shared-prefix 4-D mask may not carry over to a hybrid linear-attention model [U]. |
| `huggingface-skills:trl-training` | Pointed to `AsyncGRPOTrainer` (`trl/experimental/async_grpo/`). I read its weight transfer: NCCL broadcast to a vLLM server (`weight_transfer.py:44-160`), `max_staleness` default 4 (`async_grpo_config.py:151`). This is the basis of §4.1 and §4.6. |

---

## TL;DR

1. **Where the step goes today** (R5 shared-prefix, steps 3-15, 23.43 s) [M]: generation 11.11 s (weight sync 1.27,
   wake about 0.1, sleep 0.29 inside it); rest of the batch phase 4.65 s; update 7.67 s; outside `step_time` 0.22 s.
   **About 9-10 s of the 23.4 s have the GPU idle or nearly idle** [E]:
   - About 4.7-5.7 s are pure CPU or handover phases.
   - About 4-5 s are the decode tail, where 1-2 sequences decode alone.
2. **The largest idle source is not a sync. It is the decode tail** [M]. In R5, generation time = 5.61 s + 4.3 ms × the
   batch's longest answer (r = 0.92, 13 steps). The mean answer is 131 tokens, but the longest averages 1,271 tokens.
   38 % of steps hit the 1,536 cap. A capped answer is masked from the loss anyway (`mask_truncated_completions`,
   `configs.py:188`).
3. **The hot-path transfers and syncs, largest first** [V]:
   - 1.61 GB of fp32 `pixel_values` per step: built twice on the CPU (`fast.py:473`) and copied once, blocking and
     pageable (`transformers/trainer.py:2356`).
   - About 700 per-tensor weight pushes through `load_weights` (`trl/generation/vllm_generation.py:522`), after an
     `empty_cache` (`:501`).
   - About 1,250 small host syncs per step (counted in §2).
   - vLLM's one sync per decode step [U].
4. **The CUDA polygonizer already exists** [V]: `_gpucheck.check_masks` (`gpucheck.cpp:751-1057`, binding
   `gpucheck_py.cpp:140-159`). It is bit-exact on 195,809 cases (RTX 4090, `experiments/cuda_checker/RESULTS.md`).
   **No library caller uses it.** `targets/build.py:93`, `checker/certify.py:46` and `data/riga.py:226` call the numpy
   `polygonize`. It is also slow: 4,043 masks/s at batch 256 on the RTX 4090 vs OpenCV's 57,348. Its work scales with
   4(W+1)(H+1) slots instead of edges, and each call waits on the host 3 times. §3.2 gives a v2 design: edges only,
   rings ordered on the device, one host wait. It is bit-exact by construction (same keys).
   - The polygonizer is **offline only** (target building). It is not in the RL step, so this change satisfies the
     owner rule but saves no RL step time.
5. **Every LIDC prompt is token-identical** [V, data]. All 7,341 RL images are 512×512 (I checked 3,000 PNG headers).
   So tokenization, image-pad expansion and the rope index of the prompt are constants per run. TRL still recomputes
   them for 256 rows, twice per step, partly on the CPU.
6. **Top actions** (all [E], measure one at a time):
   - Fused, awake weight push of only the LoRA-touched weights: −1.5 s.
   - 32-image GPU pixel path plus a prompt memo: −0.6-1.5 s.
   - Larger micro-batches (32 × 8): −1.3-2.5 s (a claim of the streams review; same update under dr_grpo, verified at
     `grpo_trainer.py:3286-3289`).
   - Partial gradient checkpointing: −1-2.5 s.
   - Owner decisions on the decode tail: a cap or asynchronous generation, −0.7 to −4 s.

---

## 1. Timeline of one RL step today

Source [M]: `runs/r5/shared_rl.json`, means of steps 3-15. Code: `train_rl --shared-prefix`, FastGRPO, vLLM 0.32
colocate with sleep level 2, G = 8, 32 prompts × 8 answers, 16 micro-batches of 16 rows, beta = 1e-3.

| # | Phase (in order) | Seconds | GPU state | Evidence |
|---|---|---|---|---|
| 0 | Dataloader (OpenOnce: 32 PNG decodes), optimizer step, `log` | 0.22 | idle (CPU) | `per_step` 23.65 − `step_time` 23.43 [M]; `dataset.py:296-308`; fused AdamW over 17.4 M LoRA params is tiny [E] |
| 1 | **generate_s** (FastGRPO `_generate`, `fast.py:618-627`) | **11.11** | | [M] |
| 1a | TRL `copy.deepcopy(prompts)` (`grpo_trainer.py:2291`), chat template ×256, image processor pass 1 (ImageMemo: 32 images), 1.61 GB CPU cat (`fast.py:473`). Only `prompt_ids` is used (`grpo_trainer.py:1832-1846`) | 0.3-0.8 | idle | [E]; the pixel output is discarded on the vLLM path [V] |
| 1b | `sync_weights`: `empty_cache` (`vllm_generation.py:501`), `wake_up(weights)` (`:502`), about 700 `load_weights` calls (`:520-522`) through `lora_merged_params` (`fast.py:696-717`), `reset_prefix_cache` | 1.27 | mostly idle. 4.3 GB of D2D is about 5 ms of HBM time [E]; the rest is Python and allocator work | `time/sync_s` [M] |
| 1c | `wake_up(kv_cache)` (`vllm_generation.py:719`) | about 0.1 | idle | `time/wake_s` 0.106 (mean of the 2 wake calls; one is nested in 1b) [M] |
| 1d | `llm.generate`: 256 requests with n = 1 (`vllm_generation.py:723-743`); vLLM multimodal processing of 256 PIL images, prefill (prefix cache on), decode with CUDA graphs (`FULL_AND_PIECEWISE`, `runs/r5/r5_shared.log:39`) | about 9 (residual) | decode: busy, then a **long tail at 1-2 sequences** | fit: 5.61 s + 4.3 ms × max_len, r = 0.92 [M]; max_len mean 1,271 vs mean length 131 [M] |
| 1e | `sleep(level=2)` (`vllm_generation.py:772`) | 0.29 | idle | `time/sleep_s` [M] |
| 1f | `batch_decode`, length metrics: 3 `torch.tensor(device)` + 8 `.item()` (`grpo_trainer.py:2370-2405`) | < 0.05 | idle | [V], [E] |
| 1g | `prefetch_rewards` starts scoring on a thread (`fast.py:629-636`, `rewards.py:111-119`) | (0.005, overlapped) | checker streams | `time/rewards_s` 0.005 [M] |
| 2 | **batch − generate** (rest of `_generate_and_score_completions`) | **4.65** | | [M] |
| 2a | Pad ids and H2D (`grpo_trainer.py:2553-2600`); processor pass 2 over 256 rows (`:2645`): ImageMemo misses again (§5 B1), 1.61 GB CPU cat, blocking pageable 1.61 GB H2D (`trainer.py:2356`); `RopeCache.positions` on CPU (`fast.py:336-348`, `get_rope_index` with a Python loop per row, `modeling_qwen3_vl.py:925-985`); reward and advantage metrics (`grpo_trainer.py:2941-2972`); split / shuffle / unsplit of the 1.61 GB pixels (`grpo_trainer.py:1609-1615`, `trl/trainer/utils.py:1176-1236`) | **3.0-3.2** | **idle** | two-point fit [E] (see the note below the table) |
| 2b | Reference pass: 16 packed chunks, no grad, adapter off (`grpo_trainer.py:2845-2846`, `fast.py:118-136`) | **1.4-1.7** | busy, with sync bubbles | two-point fit [E] |
| 2c | Reward wait (prefetched) | 0.005 | – | [M] |
| 3 | **update** (16 micro-batches: forward, checkpoint recompute, backward) | **7.67** | busy. About 4.5-5.0 s of it does not scale with tokens | `step_time − batch_s` [M]; two-point fit [E] |
| | **Total `step_time`** | **23.43** | | [M] |

Two-point fit [E]. Base vs shared arm of R5, same host, same seed:
- batch − generate: 7.16 → 4.65 s.
- update: 12.37 → 7.67 s.
- Logged packed-token ratio (`shared/token_ratio`): 0.361 [M]. `packing_wins` cost ratio: about 0.40.

Assume the CPU part is the same in both arms and the reference pass scales with the ratio. Then:
- Reference pass: 3.9-4.2 s (base), 1.4-1.7 s (shared).
- CPU part: 3.0-3.2 s.
- Update: a FLOP-scaling part a = 7.4-7.8 s and a fixed part b = 4.5-5.0 s.

Caveat: packed attention is dense over W² (SDPA's boolean mask, `fast.py:205-214`), so part of "fixed" may be
attention. Only NVTX ranges can split this (§4.7).

**GPU-idle budget per step** [E]: 0.22 + 0.3-0.8 + 1.27 + 0.1 + 0.29 + 3.0-3.2 ≈ **5.2-5.9 s fully idle**. Add **about
4-5 s of tail decode** (4.3 ms × (1,271 − about 250) steps at 1-2 sequences). That is about 40 % of the step. The
full run's nvidia-smi 67-73 % (`full_rl_2026-11-08.md`) is consistent: nvidia-smi counts any running kernel as
"busy", so a batch-1 decode counts as 100 %.

Full run for comparison [M]: generation = 6.06 s + 4.1 ms × max_len. Mean max_len 586. 24 % of steps contain a
runaway answer to 1,536. **In 19 % of steps the runaway hit 1,536 while every finished answer was < 600 tokens**, so
those steps spent about 4.1 ms × (1,536 − about 600) ≈ 3.8 s decoding one answer that is then masked out. Averaged
over the run that is about 0.7 s/step [E].

---

## 2. Every host↔device transfer and forced sync in the hot path, ranked by estimated cost

Per step = 1 generation batch: 1 reference call with 16 chunks, plus 16 update calls with 1 chunk each, in the
shared-prefix path. Sizes use W ≈ 325 + 1,271 ≈ 1,600 columns (the R5 mean longest answer) [E].

| Rank | Site (file:line) | What moves / why it syncs | Per step | Bytes | Est. cost / step |
|---|---|---|---|---|---|
| 1 | vLLM decode loop (inside `llm.generate`, `vllm_generation.py:743`) | sampled ids and logprob D2H, scheduler sync per decode step | ≈ max_len ≈ 1,271 | ≈ 256 × 8 B per step | Hidden in decode; the cost is the tail (§1) [U] |
| 2 | `fast.py:473` (ImageMemo `torch.cat`) ×2 passes, then `transformers/trainer.py:2356` `.to(device)` | 256-row fp32 `pixel_values` built on the CPU twice; pass 2 copied H2D, pageable, blocking (`accelerator_config.non_blocking` reaches only the DataLoader, `configs.py:145`) | 2 cats + 1 H2D | **2 × 1.61 GB CPU memcpy + 1.61 GB H2D** (1,024 rows × 1,536 × 4 B × 256) | 0.3-0.9 s, GPU idle [E] |
| 3 | `vllm_generation.py:501-528` via `fast.py:503` | `empty_cache` (device-wide sync), wake / map, about 700 `load_weights([(name, t)])` calls, each through vLLM's name mapper; W + delta as a new bf16 tensor per LoRA layer (`fast.py:717`) | 1 | 4.26 GB D2D (2.13 B params × 2 B) | 1.27 s measured [M], mostly Python [E] |
| 4 | `vllm_generation.py:719, 772` | wake KV / sleep level 2 (cuMem map / unmap) | 2 | – | 0.1 + 0.29 s [M] |
| 5 | `fast.py:324, 342-345` (RopeCache) | `input_ids.cpu()`, `attention_mask.cpu()` per call; on the reference call (256 misses) also ids, mm types, mask, grid `.cpu()`, CPU `get_rope_index` (Python `groupby` per row, `modeling_qwen3_vl.py:985`), `pos.to(device)` pageable | 17 calls; 1 heavy | about 6.6 MB D2H per heavy call + 9.8 MB H2D (3 × 256 × 1,600 × int64) | 0.1-0.3 s [E] |
| 6 | `modeling_qwen3_vl.py:853-854` (DeepStack `hidden_states[visual_pos_masks, :]` read and write) | boolean index → `nonzero` → sync; 3 layers × 2 × (16 ref + 16 update forwards + 16 recomputes) | 192-288 | tiny | about 0.05-0.15 s [E] |
| 7 | `modeling_qwen3_vl.py:297-299` (vision attention without flash: `lengths.tolist()` per block) | 24 blocks per encode; about 16 encodes per step (2 new images per reference chunk, `fast.py:279-282`) | about 384 | tiny | 0.05-0.1 s [E] |
| 8 | `fast.py:267, 272, 276` (VisionCache `keys`, `forward`) | `torch.randint(CPU).to(device)` pageable H2D every call; `h.tolist()`, `grid.tolist()`; `grid.prod(-1).tolist()` | 3-4 per forward × 32 | 12 KB H2D per call | 0.03-0.08 s [E] |
| 9 | `trl/trainer/utils.py:1694, 1697` (fused LM head `hidden_states[mask]`, `labels[mask]`) | 2 `nonzero` syncs per head call | 64 | tiny | 0.02-0.05 s [E] |
| 10 | `grpo_trainer.py:3376-3405` (`_compute_loss` metrics: kl, entropy, 3 clip means via `global_masked_mean(...).item()`, `nanmin(...).item()`, `nanmax(...).item()`) | 7 `.item()` **before backward is queued** | 112 | tiny | 0.05-0.2 s: the GPU drains before each backward [E] |
| 11 | `fast.py:74, 125` (`same.tolist()`, `lengths.tolist()`); `fast.py:196-197` (`torch.tensor(list, device=)` for `seg`, `idx`: pageable, blocking) | per packed chunk | 4 × 32 = 128 | about 90 KB H2D per chunk | 0.03-0.08 s [E] |
| 12 | `fast.py:68` (`torch.stack(pix)` to compare pixels between rows) | no sync, but a 100 MB transient tensor and compare per chunk (16 × 1,024 × 1,536 × 4 B) | 32 | 3.2 GB D2D per step | about 2 ms each [E]; memory spikes |
| 13 | `trl/trainer/utils.py:1176` (`split_pixel_values_by_grid` grid `.tolist()`); `:1231` unsplit `torch.cat` | 1 sync, 1.61 GB D2D, plus the pixels stay resident all step | 1 | 1.61 GB D2D | about 2 ms + 1.6 GB of memory [E] |
| 14 | `fast.py:111` (`grid.prod(-1).tolist()`) | per FastForward call | 17 | tiny | < 0.01 s |
| 15 | `grpo_trainer.py:1706/1731` (8 reward functions: `torch.tensor(list, device)`); `:2941-2949` (19 `.item()`); `:2955-2956, 2971-2972` (`.tolist()`) | reward plumbing | about 30 | 8 × 1 KB | < 0.01 s, GPU idle anyway |
| 16 | `grpo_trainer.py:2370-2405` | lengths `torch.tensor(device)` ×3 + 8 `.item()` | 11 | about 6 KB | < 0.01 s |
| 17 | `grpo_trainer.py:2553-2600` | prompt / completion ids, masks, sampling logps: pad on CPU, `.to(device)` | about 6 | about 0.7 + 3.3 + 1.6 MB | < 0.01 s |
| 18 | `fast.py:381, 397` (`by_length` / `by_group_length` `order.tolist()`) | order read to the host | 1 | 2 KB | negligible |
| 19 | `modeling_qwen3_vl.py:1044` (`split_sizes.tolist()`); `:1209-1210` (`image_mask[visual_pos_masks]`); placeholder `torch_compilable_check` (`:1076-1088`) | per model forward | about 3 × 32 | tiny | 0.02-0.05 s [E] |
| 20 | `fast.py:648-649` (`timed()` `torch.cuda.synchronize`) | 7 phases + nested sync / wake / sleep | about 10 | – | < 10 ms; keep for timing |
| 21 | `fast.py:607` (`FastGRPO.log`), Trainer `tr_loss.item()`, `grad_norm.item()` | per log (`logging_steps=1`) | 3 | – | negligible |
| 22 | gpucheck on the prefetch thread: `gpucheck.cpp:1800-1802` (polygon H2D + wait), `:1569-1574`, `:1624-1626`, `:1638-1639` (3 scorer round trips), `Drain` (`:177-185`) | own non-blocking streams, pinned staging; does **not** block the trainer's stream | 1 call, 4 waits + 2 drains | about 50 KB H2D, about 40 KB D2H [E] | 5 ms, overlapped [M] |
| 23 | `gpucheck.cpp:138-149` (`Block::get` growth: `cudaFree` + `cudaMalloc`); `:1741` (`add_tables` `cudaMalloc` + sync per new cache) | **`cudaFree` is device-synchronizing**: it stalls the trainer's stream too | growth only (first steps); `add_tables` in epoch 1 | about 20 KB per table | rare; a stall at each new item in epoch 1 [E] |

Total small syncs ≈ 1,200-1,300 per step [V count, E cost 0.3-0.8 s]. That agrees with
`cuda_streams_graphs_review_2026-11-08.md` §1 (≈1,250), which I re-derived. I add the recompute pass to the DeepStack
row (that note lists 192; with checkpointing it is up to 288).

---

## 3. CPU-side stages that could move to the GPU

### 3.1 Table (saving per RL step, [E] unless marked)

| Stage | Today | GPU design | Exact? | Saving |
|---|---|---|---|---|
| **Image pixels for the trainer** | 2 CPU passes build 256 fp32 rows (1.61 GB), 1 pageable H2D | Decode the 32 distinct images once (OpenOnce already does). Upload 32 × uint8 512×512×3 (25 MB, pinned). Run the Qwen3-VL fast image processor's math on the GPU: rescale and normalize (the fused mean / std form), then reshape into 1,024 × 1,536 patches. Build the 256-row tensor with one `index_select` (about 1 ms D2D), or better keep 32 rows plus a row→image index (VisionCache needs only distinct images) | Elementwise IEEE ops on the same fp32 inputs; must be proven bit-equal to the CPU processor in a test (torch CUDA kernels are not built with fast-math) [U until tested] | 0.3-0.9 s idle CPU time |
| **Processor pass 1** (only for `prompt_ids`) | 32 image processings + a 1.61 GB cat, thrown away on the vLLM path | Skip it. Every prompt is identical (512×512 → 256 merged tokens), so memo `prompt_ids` per (prompt text, image size) and return without pixels | Exact (same ids; FastGRPO's `one_pad` check at `fast.py:680-686` still guards vLLM's expansion) | 0.2-0.5 s |
| **Tokenization / chat template** (512 renders per step) | CPU, Rust tokenizer | No GPU needed: memo (1 distinct prompt in LIDC). A GPU tokenizer would not pay | Exact | 0.1-0.3 s |
| **Rope index** (`RopeCache`, CPU `get_rope_index`) | CPU Python loop over 256 rows, plus D2H and H2D | Closed form on the GPU. The prompt's (3, P) positions are one constant tensor per run (identical prompts). Answer token t gets `max(prompt_pos) + 1 + t` on all 3 axes. Padding positions as HF sets them (1 for masked tokens). In general (varied images): a segment scan over `mm_token_type_ids` (cumsum of text tokens, grid offsets for image runs), all integer | Integer, bit-exact by test against `get_rope_index` | 0.1-0.3 s; removes 6 transfers |
| **VisionCache hashing** | GPU hash + `tolist()` per forward | Key images by the ImageMemo row→image index (known on the host at no cost); encode all 32 images once on entry, as one 32-image call, not 16 calls of 2 | Exact (no hash, so no collision risk, §5 B3) | 0.05-0.15 s; −380 syncs |
| **Text parsing of completions** (regex + `json.loads`, `geometry/text.py:102-123`) | CPU on the prefetch thread, about 2 ms / 256 answers (claimed in `cpu_gpu_trace.md`), overlapped | Possible: a token-id DFA on the GPU. Table [vocab × states] → next state (about 151,936 × 30 B ≈ 4.5 MB). One block per answer runs a parallel prefix of state-transition maps (composition is associative), then emits ints straight into gpucheck's packed `xy` / `ring_start` arrays: no host polygon packing, no H2D | Must reproduce `parse`'s acceptance exactly. The error codes (`_diagnose`) are only logged; they can stay on the CPU for rejects | **about 0 s** (already hidden). Low priority despite the GPU rule |
| **Advantages** (`advantages.py`, numpy over 256 values) | CPU, microseconds | Could be fused into the scorer as one more kernel. But TRL needs Python lists (`rewards.py:145-147`) | float64 order must match | 0 s |
| **Data loading** (32 PNG decodes per step) | CPU main thread, about 0.1 s [E] | Pre-decode the train set once into a uint8 store: 7,341 × 512² × 3 = 5.8 GB pinned host, or GPU-resident if single-channel (1.9 GB) [U: channels]. Per step: gather 32 frames (25 MB). nvImageCodec does not decode PNG on the GPU (it uses CPU codecs) [U] | Exact (same pixels) | 0.1 s |
| **vLLM image input** | vLLM re-processes 256 PIL images (its mm cache may dedupe) [U] | Pass the trainer's 32 GPU vision embeddings as `image_embeds` (vLLM supports embedding inputs for Qwen-VL [U for Qwen3-VL DeepStack]) | Different code path; must equal the vLLM encode | 0.2-1 s [U] |
| **Polygonizer** (offline: `targets/build.py:93`, `certify.py:46`, `data/riga.py:226`) | numpy `polygonize` (`checker/polygonize.py`) | **Already exists on the GPU as `_gpucheck.check_masks`; it needs wiring and a v2 (§3.2)** | Bit-exact (195,809 / 195,809 on an RTX 4090 [M]) | Offline only; 0 s per RL step |

### 3.2 The CUDA polygonizer

#### What exists [V]

`gpucheck.cpp:751-1057` already ports `polygonize.py` + `bitquad.py`:
- **Slot.** Slot (vertex, direction) key `(y(W+1)+x)·4+d`, equal to polygonize's edge key (`polygonize.py:130`).
- **Successor.** Taken from the 2×2 bit-quad at the edge's end, left turn first (`InitSlots`, `:806-842`), equal to
  `link` (`polygonize.py:72-89`).
- **Ring label.** Min-label pointer jumping (`LabelRound`, `:845-858`).
- **Corner rank.** Weighted Wyllie list ranking on `pred` (`RankInit` / `RankRound`, `:869-896`).
- **Ring order.** Ring key = first east edge in raster order (`RingStats`, `:898-915`, `atomic_min` at `:908`). This
  equals polygonize's `np.unique(lab)` order, because `boundary_edges` enumerates top-side (east) edges first, in
  raster order (`polygonize.py:60-65`).
- **Output.** Corners scattered by rank (`ScatterCorners`, `:917-929`).

Exactness argument (my check):
- Every quantity is an integer.
- Atomics are only integer add and min, so the order of operations does not matter.
- The ring key is unique within a mask: one east edge leaves a vertex.
- Extra jump rounds are idempotent: `jump_rounds(max_slots)` ≥ `_jumps(n_edges)`.
- The start edge is always a corner. It arrives from below or from the right (d ∈ {2,3}) and leaves east or south
  (d ∈ {0,1}).

Test coverage: `tests/test_cuda_checker.py:131-158` (fuzz + Phase 0 LIDC masks) and
`experiments/cuda_checker/eval.py`.

#### What is wrong with it

1. **Not used.** Library callers import the numpy version (`checker/__init__.py:6`). The owner rule is not met in
   practice.
2. **Slow.** 4,043 masks/s at 22×21 × 256, and 19 masks/s at 512² × 4, on an RTX 4090. OpenCV does 57,348 and 22 [M].
   - Work ∝ slots: every kernel walks all 4(W+1)(H+1) slots of a vertex row serially (`for_row_slots`, `:795-804`).
     A 22×21 mask has 2,024 slots but about 60 edges [E].
   - The jump rounds use `log2(max_slots)+1` (12 here) instead of `log2(max_edges)+1` (7).
   - Parallelism is one thread per vertex row: 23 threads per mask.
   - 3 blocking host round trips per chunk (`:1019-1020`, `:1025-1027`, `:1054-1055`).
   - The ring sort runs on the host (`:1032-1034`).
3. Scratch memory is about 60 B per slot (9 i32 arrays + 3 i64), so ≤ 0.5 GB per chunk (`kChunkSlots`, `:937`).

#### v2 design (an addition to gpucheck.cpp, same functor style; numpy stays the oracle)

```
masks (bit-packed rows, like the scorer's words)
  K1 CountEdges   1 thread / (mask, vertex row, 32-col word): codes from rows y-1, y; per direction d,
                  popc(has_edge_d bits) -> cnt[row]            (word-parallel: east = fg(y) & ~fg(y-1), ...)
  K2 Scan         exclusive scan of cnt (CUB DeviceScan)        -> edge offsets; E_total
  K3 EmitEdges    same walk, write key=(y(W+1)+x)*4+d in key order -> keys[] sorted within each mask
                  (+ bit-quad n1..n4, nD with warp reduce + one atomic per warp)
  K4 Succ         1 thread / edge: end vertex quad code, turn order {left, straight, right} (as InitSlots);
                  successor index = binary search of its key in the mask's sorted key range
  K5 Label        ceil(log2 E_mask_max)+1 rounds of lab=min(lab,lab[jmp]); jmp=jmp[jmp]   (1 thread / edge)
  K6 Pred+Rank    pred[succ[e]]=e; weighted Wyllie on pred, stop at lab==e (as RankInit/RankRound)
  K7 RingStats    root=lab[e]: atomicAdd area2 (i64), corner count; atomicMin first-east key; root==e -> ring record
  K8 Order        ring sort key (mask<<32 | first_east_key)  CUB DeviceRadixSort (or segmented, rings are few)
                  exclusive scan of corner counts -> pts offsets; per-mask ring counts -> ring_start (scan)
  K9 Scatter      corner e -> pts[base[lab[e]] + w[e]] = (x, y)
  copy-out        ONE D2H of (quads, ring_start, pts_start, ring_area, pts) sized from E_total
                  (corners <= edges), or zero copies: return device buffers through DLPack to torch
```

- **Exactness.** Keys, the successor rule, labels (min key over the cycle), ranks and the order key are the same
  integers as polygonize's, so the output is identical by construction. Keep the existing fuzz, LIDC and RIGA tests,
  and add a v1-vs-v2 equality test on the 195,809-case set.
- **Expected speed** [E]:
  - Work per round: ×30 less on LIDC masks.
  - Rounds: 12 → 7.
  - Host waits: 3 → 1 (or 0 with DLPack).
  - For 256 small masks the call becomes launch-bound: about 25 launches + 1 wait ≈ 0.2-0.5 ms, i.e. ≫ 100k
    masks/s.
- **Wiring.**
  - `objdel.checker.polygonize` dispatches to `_gpucheck.check_masks` when `_gpucheck.build() == "cuda"` and returns
    the same `Ring` list.
  - Add a batched `polygonize_many(masks)`, so `targets/build.py` (per rater × object) and `dataset.prepare` batch
    a whole item, or many items, per call.
  - `certify` uses the returned `quads` too: bit-quads come free.
  - The numpy file stays as the oracle in `tests/`.
- **Checks.**
  - Run compute-sanitizer memcheck, racecheck (shared memory, if K3 uses it), initcheck and synccheck on the new
    kernels (the NVIDIA `compute-sanitizer-diagnose-cuda-bugs` routing). The repo already runs memcheck and
    initcheck before rentals (`full_rl_2026-11-08.md`).
  - ncu on K4/K5 for memory throughput.

---

## 4. Structural ideas for maximum utilization

### 4.1 ADR-1 (Proposed): synchronous colocate vs asynchronous generation on one RTX 5090

**Context.** Generation and training alternate. Both CPU phases (§1 rows 1a, 2a) and the decode tail leave the GPU
idle about 40 % of the step [E]. Weights are pushed every step (1.27 s) [M].

| Option | Complexity | Expected step | Policy lag | Main risk |
|---|---|---|---|---|
| A. Keep sync colocate, remove the idle phases (§4.2) | low-medium | 23.4 → about 16-18 s [E] | 0 | none new |
| B. Async on one GPU: vLLM server process + trainer process sharing the GPU (MPS), CUDA-IPC weight push, generation of batch k+1 overlaps training on batch k (TRL `AsyncGRPOTrainer` pattern, `max_staleness` = 1) | high | about max(gen ≈ 9 s, train ≈ 12 s) + contention ≈ 13-14 s [E] | 1 step (IS correction already on, `fast.py:581-602`) | NCCL cannot pair two ranks on one GPU, so TRL's `weight_transfer.py` NCCL path does not apply. Memory: 11.2 GB vLLM + 17 GB trainer peak ≈ 28 of 32 GB [M peak from `full_rl_2026-11-08.md`] |
| C. Partial rollouts (cap decode per step; unfinished answers continue next step) | high (not in TRL) | removes the tail | mixed | off-policy tokens inside one answer |

**Decision (proposed).** A first: it is cheap and keeps the update exactly on-policy. B only after A's measurements
show the remaining idle time is generation-bound. Revisit B when moving to two or more GPUs, where it becomes the
default (§4.6).

### 4.2 One GPU, synchronous: concrete changes, ranked by saving per step [E]

1. **Fused weight push, vLLM awake** (−1.5 s).
   - Sleep off: R4 measured awake = asleep in step time, and memory fit (`r4_speed_2026-11-08.md` [M]).
   - At startup, map each LoRA layer to its vLLM parameter slice (qkv_proj and gate_up_proj are stacked).
   - Per step, for each of the 196 LoRA layers: `slice.copy_(W_base)`, then
     `slice.add_(((B@A)*s).to(bf16))`. This is bit-identical to `lora_merged_params` (`fast.py:717`: bf16 W + bf16
     delta).
   - Skip the vision tower, embeddings and norms: about 1.45 of 4.26 GB, never changed.
   - About 15 GB of HBM traffic ≈ 10 ms, plus about 800 launches: about 20-50 ms in total, vs 1.27 + 0.1 + 0.29 s.
   - vLLM's CUDA graphs keep valid addresses because the copy is in place. Keep `reset_prefix_cache`.
2. **Decode tail** (owner decision; −0.7 s full-run regime to −4 s R5 regime).
   - (a) `max_completion_length` near the longest LIDC target. The review claims 1,182 tokens; I did not verify it.
   - (b) Overlap: drive vLLM's engine step by step and start the reference pass on finished groups on a second
     stream. Complex.
   - (c) B in §4.1.
   - Early abort of a non-terminated answer is not reward-neutral: TRL masks truncated answers from the loss but
     trains on terminated invalid ones (`grpo_trainer.py:2604-2611`).
3. **Pixel path + prompt memo** (−0.6-1.5 s): §3.1 rows 1-3.
4. **Micro-batch 32 × GA 8** (−1.3-2.5 s, claimed by the streams review [E]).
   - Halves the fixed per-micro-batch cost (b ≈ 4.5-5 s, §1).
   - Same update: dr_grpo loss = Σ / (B_mb · L_max) / GA (`grpo_trainer.py:3286-3289`), so the normaliser is
     256 · L_max either way.
   - Needs one memory test on the box.
5. **Partial gradient checkpointing** (−1-2.5 s).
   - Packing cut activations about 2.8× (token ratio 0.361 [M]).
   - Turn checkpointing off on k of the 28 layers (`GradientCheckpointingLayer.gradient_checkpointing`, per layer in
     transformers 5).
   - Competes for memory with change 1 (vLLM awake).
6. **Sync cleanup** (−0.3-0.8 s; rows 5-11 of §2):
   - one host read per call;
   - build `seg` / `idx` on the device, or pinned with `non_blocking`;
   - precompute DeepStack and fused-head indices once per forward (index_select instead of boolean indexing);
   - keep loss metrics as device tensors (a `_compute_loss` override, read once in `log`, as `is_stats` already does
     at `fast.py:595-615`);
   - encode the 32 images in one call;
   - `TRANSFORMERS_DISABLE_TORCH_CHECK=1` for the placeholder check [U name: the review's claim].
7. **Reference pass in 64-row chunks, no grad** (−0.2-0.5 s): fewer launches, longer packed rows.
8. **Block-sparse attention for the packed rows** (≤ −0.5 s).
   - FlexAttention block mask: prompt causal + each answer → prompt + itself.
   - Attention work about 4× lower than dense W² (462k vs 1.88M query-key pairs per group at the R5 sizes [E]).
   - Attention is about 11 % of a row's cost (n² / 12,000 at n ≈ 1,373, `fast.py:81-93`), so the gain is small.

Expected after 1, 3, 4, 6, 7: about 16-18 s/step [E]. With the tail fixed or B: 12-14 s [E].

### 4.3 CUDA IPC weight sync (for B, or for any vLLM in another process)

- The trainer allocates one persistent bf16 "push buffer" in buckets (for example 4 × 512 MB). It exports CUDA IPC
  handles once (`torch.multiprocessing` reductions / `cudaIpcGetMemHandle`). The vLLM worker opens them once.
- Each step: the trainer writes merged W + delta into the buckets on a side stream and records an event (exported via
  `cudaIpcGetEventHandle`). An RPC tells vLLM to wait on the event and copy each bucket slice into its parameters.
- Cost about 2 × 2.8 GB / 1.79 TB/s ≈ 3 ms + RPC [E]. NCCL is not possible between two processes on one GPU.
- On many GPUs, use NCCL broadcast (TRL `weight_transfer.py`).

### 4.4 CUDA Graphs

- **Decode: already on.** vLLM `FULL_AND_PIECEWISE`, `enforce_eager=False` (`runs/r5/r5_shared.log:39` [V]).
- **Training: per-layer graphs** (`torch.cuda.make_graphed_callables`, or TE-scoped graphs as in the NVIDIA
  `nemo-mbridge-perf-cuda-graphs` skill) for the update micro-batches only.
  - Pad packed widths to multiples of 256 (about 14 shapes).
  - Turn checkpointing off for graphed layers.
  - The reference pass (adapter off) stays eager.
  - Blockers [V in code]: the syncs sit outside the decoder layers (good), but the DeepStack boolean index sits
    *inside* layers 5, 11 and 17 (`modeling_qwen3_vl.py:834-854`). Those 3 layers need the index fix first.
  - Gain 0.5-1.5 s [E], high risk. Do it after §4.2 change 4 (bigger micro-batches already cut launch overhead).

### 4.5 Fused log-prob / entropy and FP8

- **Liger-style fused head: already in.** TRL `add_fused_lm_head` + Triton `_ChunkedLogProbFunction`
  (`trl/trainer/utils.py:1569-1720`, `trl/kernels/chunked_logprob.py`). It scores only label positions and never
  builds (B, T, 151,936) logits [V].
  - Remaining waste: the 2 boolean-index syncs (`utils.py:1694, 1697`). The packed path already knows the label
    positions (`fast.py:197` `idx`), so pass precomputed indices.
- **FP8 on Blackwell (sm_120).**
  - (a) vLLM FP8 weights for generation: decode is bandwidth-bound, so about 1.3-1.6× faster decode [E]. But the
    sampler then differs from the bf16 trainer, the IS ratio grows, and the KL spikes already seen
    (`full_rl_2026-11-08.md`: step 156) may worsen.
  - (b) FP8 for the frozen base GEMMs in the update (torchao float8, LoRA path in bf16): GEMMs 1.3-1.5× faster [E],
    but the update is now dominated by fixed cost (§1).
  - Both change numerics. Measure after the structural fixes.

### 4.6 Many GPUs

- **N GPUs, colocate DDP** (TRL as is, `per_device` sizes are multiples of G, `configs.py:185`).
  - LoRA gradients all-reduce 70 MB/step (17.4 M fp32 params [V shapes from `runs/p1/ckpt/config.json`]).
  - Scales roughly linearly, but keeps every idle phase.
- **Disaggregated async** (TRL `AsyncGRPOTrainer`): k vLLM GPUs (data parallel) + m trainer GPUs.
  - Weights by NCCL broadcast (`weight_transfer.py:136-160`).
  - RTX 5090s have no NVLink, and GeForce P2P is usually off, so NCCL goes through host memory: 4.3 GB ≈ 0.2-0.5 s
    [E].
  - Push only the LoRA (about 70 MB fp32 / 35 MB bf16) and use vLLM's native LoRA serving. Trade-off: LoRA kernels
    slow decode by maybe 10-30 % [E].
  - Ratio: generation ≈ 9 s vs training ≈ 12 s per 256 answers today [E], so about 1:1 GPUs. After the tail fix,
    generation shrinks: about 1:2.
- **Sequence packing across the data-parallel ranks**: sort by group length globally before sharding.

### 4.7 NVTX + Nsight Systems measurement plan (before any change)

1. **NVTX ranges** (no behaviour change):
   - `FastGRPO.timed` → `torch.cuda.nvtx.range_push/pop` per phase. Add new ranges: `pass1`, `pass2+H2D`, `rope`,
     `ref_chunk[i]`, `update_mb[i].fwd/bwd`, `optimizer`, `log`.
   - The prefetch thread: `reward_prefetch`.
   - gpucheck C++: `nvtxRangePushA` around each API call and kernel stage (nvtx3 is header-only).
   - vLLM: `enable_layerwise_nvtx_tracing` exists in its observability config (`runs/r5/r5_shared.log:39`).
2. **Capture.** Steps 4-5 only: `torch.cuda.profiler.start()/stop()` in a callback, and
   `nsys profile --capture-range=cudaProfilerApi --capture-range-end=stop -t cuda,nvtx,osrt,cublas --cuda-memory-usage=true --cudabacktrace=sync --python-sampling=true -s cpu`.
   Add `--gpu-metrics-devices=all` if the container allows it: SM active % and DRAM bandwidth per phase, which shows
   the decode tail.
3. **Reports.**
   - `nsys stats -r nvtx_sum,cuda_gpu_kern_sum,cuda_api_sum,cuda_gpu_mem_size_sum,cuda_gpu_mem_time_sum`.
   - Per-step KPIs: GPU-active % per NVTX range; idle gaps > 1 ms; `cudaStreamSynchronize` count; H2D / D2H bytes;
     update MFU (useful FLOPs / time); decode SM-active % vs batch size.
   - ncu only for gpucheck kernels (raster, the polygonizer rounds).
4. **Success criteria.**
   - Fully idle time per step < 1.5 s.
   - Syncs per step < 200.
   - H2D per step < 50 MB.
   - Tail decode at < 10 % SM active < 1 s.
5. **Cost.** One 20-minute box session [E]; the owner rule needs approval first.

### 4.8 Model swap note (huggingface-best)

Qwen3.5-2B (2.27 B params) fits the same budget. TRL ships a tree gated delta-rule kernel for packed forwards
(`trl/experimental/async_grpo/tree_gated_delta_rule.py`). That suggests linear-attention layers, where the
shared-prefix 4-D mask (`fast.py:205-216`, SDPA or eager only) would not apply as written [U]. Check this before any
swap.

---

## 5. Bugs and correctness risks

| # | Severity | Where | Finding |
|---|---|---|---|
| B1 | perf (claim false) | `fast.py:418-427` vs `trl grpo_trainer.py:2291` | **ImageMemo does not share pixels across passes.** TRL deep-copies the prompts (with their PIL images) before pass 1, so pass 1 sees copies and pass 2 sees the originals. `id()` differs, and pass 2 processes the 32 images again. The docstring ("prepares each image once") is wrong; it is twice. Within a pass the dedupe works (deepcopy keeps aliasing). |
| B2 | perf | `fast.py:473` | ImageMemo builds the full 256-row fp32 tensor in pass 1 too, where TRL discards it (1.61 GB memcpy, idle GPU). |
| B3 | low-probability silent error | `fast.py:262-272` | VisionCache keys are a linear hash mod 2^64 of the pixel bit patterns. A collision would reuse another image's features without any error. Probability is about 2^-64 per pair [E], so the risk is negligible but silent. Keying by the ImageMemo / OpenOnce identity removes it and the syncs (§3.1). |
| B4 | estimate only | `fast.py:85-93` | `packing_wins` slices `lens` into equal groups (`size = len(lens) // groups`). `prompt_groups` can merge two adjacent groups of the same item: disagreement items are repeated ×4 (`dataset.py:81`, `train_rl.py:51`), so one run can hold 16 rows. The cost model is then wrong (correctness is fine: the rows share prompt and image). |
| B5 | memory spike | `fast.py:68` | Comparing pixels by stacking 16 × 6.3 MB per chunk adds 100 MB transient tensors 32 times per step. Use the image identity instead. |
| B6 | doc / priority | `gpucheck.cpp:66-71` | "Lowest priority, the trainer's kernels go first" is not achieved. `cudaDeviceGetStreamPriorityRange`'s *least* priority is 0, the same as the default stream's, so the checker and the trainer have equal priority. To favour the trainer it would need a higher-priority (negative) stream. Harmless today: the checker takes 5 ms. |
| B7 | stall | `gpucheck.cpp:143-145`, `:1741`, `:1772` | `cudaFree` in `Block::get` growth and in `drop_tables` is device-synchronizing and can stall the trainer mid-pass. `add_tables` calls `cudaMalloc` + a wait for each new cache in epoch 1 (more with `--train-views`). Fix: pre-size from the batch caps; batch-upload all preloaded caches at startup. |
| B8 | dead weight | `train_rl.py:50, 61` | `fork_pool` forks up to 16 workers that train_rl never uses (the pool is only shut down). Memory only (re-verified). |
| B9 | semantics trap | `configs.py:188` + TRL `grpo_trainer.py:2604-2611` | A runaway answer to 1,536 tokens costs about 3.8 s of decode [E from the fit] and contributes nothing (masked). Any early-stop rule changes which answers are "truncated" (masked) vs "terminated invalid" (trained with a negative advantage). This needs an owner decision, not a silent speed fix. |
| B10 | resolved-but-unverified-on-GPU | `fast.py:495-504, 696-717` | The bf16 base-drift fix (no merge / unmerge) is tested on the Mac only (`full_rl_2026-11-08.md`). The R4 and R5 logs show step 1 identical to the full run, but no run checks parameters bit-for-bit after N syncs on the GPU. Add that assert to the next preflight. |
| B11 | stale doc | `docs/research/notes/cpu_gpu_trace.md` | Its headline (scoring on the CPU, 16 workers, mask pickling) predates `score_polygons`. Anyone reading it gets the wrong picture. Mark it superseded. |
| B12 | stale doc | `gpucheck.h:1-6` | The header advertises "the mask polygonizer of objdel.checker", but the library never calls it (§3.2). |

---

## 6. What I could not verify

- vLLM 0.32 internals: per-request multimodal processing and caching, and the per-decode-step syncs. vLLM is not in
  the local venv.
- Bit-equality of a GPU image processor vs the CPU one (needs a test).
- Box memory headroom for "awake + partial checkpointing + micro-batch 32" together (needs one GPU session).
- The split of the 4.65 s and 7.67 s phases comes from a two-point fit, not from a trace. §4.7 resolves it.
