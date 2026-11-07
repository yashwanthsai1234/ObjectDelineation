# CPU vs GPU in one GRPO step of the R2 G-probe (code as on branch `rl-ready`, 2026-11-07)

Scope: one optimizer step of `python -m experiments.phase2.g_probe --model qwen3-vl-2b --ckpt ...` on one RTX 5090.
Stack as installed: TRL 1.17.1, transformers 5.20.1, peft 0.21.2, accelerate 1.15.0, torch 2.14.1; vLLM 0.32.0 runs
only in the box image and is **not installed locally**, so I could not read vLLM's internals. Every vLLM-internal
claim below is marked UNVERIFIED.

Labels used below:
- MEASURED-R2 means the number comes from `runs/r2/r2_qwen3vl2b.log` (20 logged steps).
- MEASURED-P1 means the P1 rental on 2026-11-06: `runs/p1/profile/{fast,trace}`. That box had an EPYC 9554 CPU (R2
  had an EPYC 7742), and P1 ran code from before C5 and C10.
- MEASURED-MAC means I timed it in this checkout on the Mac (Apple Silicon, `python -I`, real `data/prepared`).
- ESTIMATE means arithmetic from the code plus those measurements.
- UNVERIFIED means I could not check it.

## 0. Settings in effect (checked in code)

| Setting | Value | Where |
|---|---|---|
| per_device_train_batch_size / gradient_accumulation_steps | 16 / 16 | `src/objdel/train/configs.py:185` |
| steps_per_generation / generation_batch_size | 16 / 256 (TRL default: steps_per_generation = GA) | `trl/trainer/grpo_config.py:1115-1117` |
| num_generations, prompts per step | 8, 32 | `configs.py:178,183` |
| gradient_checkpointing | True (use_reentrant=False by transformers default) | `configs.py:186`; `transformers/modeling_utils.py:3135` |
| bf16 | autocast; base weights bf16 (ckpt `config.json` dtype bf16), LoRA r16 all-linear on the text model (≈17.4 M fp32 params, ESTIMATE from shapes) | `configs.py:95,183`; `experiments/phase2/g_probe.py:81,83` |
| LR / steps | 5e-5 LoRA, 20 steps, logging_steps 1, log_completions False | `g_probe.py:84-86`; `experiments/phase2/prereg.py:61-69` |
| vLLM | colocate, gpu_memory_utilization 0.35, sleep mode on (level 2), max_num_seqs 256, max_num_batched_tokens 4096, logprobs=0, n=1 per request (256 requests) | `configs.py:192-193`; `trl/trainer/grpo_trainer.py:1169-1202`; `trl/generation/vllm_generation.py:368,377,677` |
| IS correction | mode `sequence_mask` (TRL default), clip_max 3.0; C10 skips the old-logprob pass | `grpo_config.py:330,340`; `src/objdel/train/fast.py:276-294` |
| optimizer | `adamw_torch_fused` (transformers default for torch ≥ 2.8), max_grad_norm 1.0 | `transformers/training_args.py:800-807` |
| dataloader | num_workers 0 (TrainingArguments default), identity collator, batch = 256 rows | `grpo_trainer.py:1011,1321-1341` |
| scoring pool | `serve.pool_size()`: min(16, CPUs-1). The R2 box has 128 CPUs and box.txt shows no cgroup quota, so **16 forked workers** | `src/objdel/train/serve.py:145-147`; `runs/r2/box.txt:45` |
| `_gpucheck` | **CUDA build** on the R2 box (`_gpucheck: cuda build, 1 device(s)`) | `runs/r2/r2_qwen3vl2b.log:72-73` |
| image tokens | 512×512 LIDC image becomes grid (1,32,32): 1,024 pixel rows × 1,536 fp32 and 256 merged tokens. Prompt is 325 tokens wide (MEASURED-MAC, Qwen2VLImageProcessor from the base model) | processor of `Qwen/Qwen3-VL-2B-Instruct`; the run loads it from the ckpt (UNVERIFIED that the ckpt's processor config is identical) |
| reward masks | LIDC crops are 50×50 (81 % of the oversampled rows), 100×100 (19 %), 125×125 (0.3 %) (MEASURED-MAC on `data/prepared/lidc/train.jsonl`) | |

## Headline

1. **The scoring is on the CPU. Only the checker runs on the GPU.**
   - The CUDA checker rasterizes and validates all 256 answers in one call. The bit-packed masks then come back to the
     host (≈190 kB) and are unpacked there with numpy.
   - Everything after that runs on the CPU in 16 forked worker processes or in the main process: parsing, IoU,
     signature/topology, band depth/F5, area guard, joint topology, sig_id, and advantages_v4/topology_reward.
   - The masks are pickled to the workers (≈0.8 MB per step, MEASURED-MAC). The results go back as Python lists and
     are copied to the GPU as 8 small tensors.
   - It costs about 10-50 ms per step and is fully hidden. `time/rewards_s` = 0.0037 s mean (MEASURED-R2).
2. **The largest CPU cost in the step is not scoring. It is TRL's dataloader.**
   - transformers' `get_batch_samples` fetches 16 generation batches per optimizer step. TRL's `RepeatSampler` makes
     those 16 batches the *same* 256 rows, and only the first one is used.
   - That is 4,096 row reads per step, each decoding a 512×512 PNG in the main process while the GPU sits idle.
   - MEASURED-P1: 16 calls of `_SingleProcessDataLoaderIter._next_data` took 12.96 s CPU in one step, and nvidia-smi
     shows a contiguous idle gap of about 14.5 s per step.
   - In R2, wall time minus TRL's `step_time` = **22.6 s per step**. The dataloader is most of that (ESTIMATE
     18-22 s, because the EPYC 7742 is slower).
3. **Two full Qwen3-VL image-processor passes over 256 images per step, both on the CPU.**
   - Each pass produces 1.61 GB of fp32 `pixel_values`. The first pass's pixel tensor is thrown away; vLLM gets PIL
     images and processes them again (UNVERIFIED for vLLM).
   - The second pass's 1.61 GB goes to the GPU as a blocking, pageable copy.
   - MEASURED-MAC: 0.75-1.5 s per pass.
4. Generation time is set by the longest answer: `time/generate_s` ≈ 4.1 s + 5.6 ms × completions/max_length, with
   r = 0.90 over the 20 steps (MEASURED-R2 fit).
5. The trainer forces about 1,200 host-device syncs per step: `.item()`, `.tolist()`, `.cpu()`, boolean indexing
   and blocking `.to()` (ESTIMATE from code, itemized in §2). vLLM adds roughly one per decode step (UNVERIFIED).

---

## A. The step in execution order

Columns: **Dev** = where the work runs. **H↔D** = data crossing host and device per step (sizes at max completion
length L = 1536, so the full width W = 325 + 1536 = 1861). **Syncs** = host-device synchronizations it causes.
A blocking `.to(cuda)` from pageable memory counts as one sync: torch copies, then synchronizes the stream.

### A0. Trainer loop: fetch the generation batch (before anything else)

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 0.1 | `get_batch_samples` | `next(epoch_iterator)` ×16 (GA = 16) | CPU (main thread) | none | 0 | `transformers/trainer.py:1856, 2261-2276` |
| 0.2 | RepeatSampler | yields 32 indices × 8 copies, the whole chunk **repeated 16 times**, so all 16 batches of 256 rows are identical | CPU | none | 0 | `trl/trainer/utils.py:925-944`; `grpo_trainer.py:1343-1379` |
| 0.3 | Dataset `__getitems__` | arrow → Python for 256 rows (prompt, key, W, H, names) and `PIL.Image.open(path)` + `image.load()`: **a full 512×512 PNG decode per row, 4,096 per step** | CPU | none (PIL objects; `pin_memory` touches only tensors) | 0 | `datasets/features/image.py:175,193`; rows from `src/objdel/train/dataset.py:214-235` |
| 0.4 | identity collate, accelerate shard | list of dicts passed through | CPU | none | 0 | `grpo_trainer.py:1011` |
| 0.5 | `_get_num_items_in_batch` | `"labels" in batch_samples[0]` (False) | CPU | none | 0 | `trainer.py:2278-2300` |

Cost: MEASURED-P1 12.96 s for 16 calls (810 ms each). MEASURED-MAC 0.75-0.88 s per 256-row fetch. Only batch 1 of 16
is used (`grpo_trainer.py:1609-1615`); TRL ignores the other 15.

### A1. Micro-step 0: `training_step` → `_prepare_inputs` → `_generate_and_score_completions`

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 1.1 | `Trainer.training_step` | `model.train()` (module walk) | CPU | – | 0 | `trainer.py:2023-2050` |
| 1.2 | `FastGRPO._prepare_inputs` | swaps TRL's shuffle for `by_length` (dr_grpo, GA % spg == 0) | CPU | – | 0 | `fast.py:296-308` |
| 1.3 | TRL `_prepare_inputs` | `_step % 16 == 0` → generate | CPU | – | 0 | `grpo_trainer.py:1591-1620` |
| 1.4 | `FastGRPO._generate_and_score_completions` | clears VisionCache/RopeCache, keeps the dataset rows, C10: sets `vllm_importance_sampling_correction=False` for this call | CPU | – | 0 | `fast.py:271-294` |
| 1.5 | build prompts | `prepare_multimodal_messages` inserts each row's PIL image into its message (256×) | CPU | – | 0 | `grpo_trainer.py:2505-2530` |

### A2. Generation (`FastGRPO._generate` → TRL `_generate`): `time/generate_s` 11.7 s mean, 13.7 s at step 1

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 2.1 | `patch_vllm_prompts` | idempotent check (wraps `llm.generate` once) | CPU | – | 0 | `fast.py:361-363, 392-412` |
| 2.2 | `copy.deepcopy(prompts)` | deep copy of 256 message lists with PIL images (MEASURED-MAC 0.017 s) | CPU | – | 0 | `grpo_trainer.py:2291` |
| 2.3 | `_tokenize_prompts` → `processor.apply_chat_template(tokenize=True)` | chat template, tokenization, and the **full Qwen3-VL image preprocessing of 256 images**. Outputs `pixel_values` (262,144 × 1,536 fp32 = **1.61 GB**), `image_grid_thw`, `mm_token_type_ids`. **Only `prompt_ids` is used**: the vLLM branch ignores `multimodal_fields`. MEASURED-MAC 0.75-1.31 s | CPU (torch CPU ops, multithreaded) | none | 0 | `grpo_trainer.py:1753-1830` (call at 1806); unused at `1832-1851` |
| 2.4 | `sync_weights` → `empty_cache()` | torch caching allocator frees cached blocks (cudaFree: device-wide sync) | GPU mem API | – | 1 device sync | `trl/generation/vllm_generation.py:501` |
| 2.5 | `llm.wake_up(tags=["weights"])` | vLLM sleep level 2 discarded the weights; remaps physical memory | GPU mem API | – | UNVERIFIED (cuMem calls) | `vllm_generation.py:502` |
| 2.6 | `merge_adapter()` | per LoRA layer (196): delta = B@A·scale (fp32), `W.data += delta` in bf16 | GPU | – | 0 | `vllm_generation.py:459`; `peft/tuners/lora/layer.py:1007` |
| 2.7 | param loop → `load_weights([(name, param)])` | **one call per parameter** (≈700 tensors including the frozen vision tower), ≈4.3 GB bf16 copied into vLLM's buffers. Python overhead per call is CPU | CPU loop + GPU D2D | D2D ≈4.3 GB (ESTIMATE from config shapes), no host crossing | UNVERIFIED | `vllm_generation.py:446-490, 520-522` |
| 2.8 | `unmerge_adapter()` | `W.data -= delta.to(bf16)` | GPU | – | 0 | `vllm_generation.py:480`; `peft/.../layer.py:1036` |
| 2.9 | `reset_prefix_cache` | vLLM bookkeeping | CPU | – | UNVERIFIED | `vllm_generation.py:528` |
| 2.10 | `wake_up(tags=["kv_cache"])` | KV cache memory remapped | GPU mem API | – | UNVERIFIED | `vllm_generation.py:719` |
| 2.11 | build vLLM requests | 256 `{"prompt_token_ids", "multi_modal_data": PIL}`; `one_pad_per_image` cuts each run of 256 image pads to 1 (Python list comprehension over ids) | CPU | – | 0 | `vllm_generation.py:723-731`; `fast.py:239-247, 402-403` |
| 2.12 | `llm.generate` | vLLM V1, in-process with `external_launcher` (UNVERIFIED). CPU: multimodal input processing of the PIL images (HF processor inside vLLM; 32 distinct images × 8 copies, so its mm cache may dedupe: UNVERIFIED), scheduler, input prep, incremental detokenization, logprob objects. GPU: ViT encoder, prefill (prefix cache on: `enable_prefix_caching=True` in the env-check log), decode with CUDA graphs, sampling (T = 1, no filters) | CPU + GPU | H2D: image tensors and per-step inputs (UNVERIFIED sizes). D2H: sampled ids + logprobs every decode step | ≈1 per decode step, ≈L_max steps (≈1,350-1,536) (UNVERIFIED) | `vllm_generation.py:743`; generate ≈ 4.1 s + 5.6 ms × L_max (MEASURED-R2 fit) |
| 2.13 | `one_pad` post-check | compares vLLM's expanded prompt ids with the trainer's (256 lists of 325) | CPU | – | 0 | `fast.py:404-409` |
| 2.14 | `extract_logprobs` + top-1 | Python loops over ≈32 k token logprob dicts (`sorted`, NaN check) | CPU | – | 0 | `vllm_generation.py:70-100, 747`; `grpo_trainer.py:1851` |
| 2.15 | `llm.sleep(level=2)` | frees vLLM weights and KV (level 2: nothing offloaded to the host) | GPU mem API | – | UNVERIFIED | `vllm_generation.py:772` |
| 2.16 | `batch_decode` | 256 completions to text (Rust tokenizer) | CPU | – | 0 | `grpo_trainer.py:2338` |
| 2.17 | length/clip metrics | `torch.tensor(..., device)` ×3 and `.item()` ×8 (num_tokens, mean/min/max length, clipped_ratio, 3 terminated stats) | CPU→GPU→CPU | H2D 3 × ≈2 kB | ≈11 | `grpo_trainer.py:2370-2405` |
| 2.18 | `timed("generate")` end | `torch.cuda.synchronize()` | – | – | 1 device sync | `fast.py:382-390` |

### A3. Reward prefetch starts (C9). Background thread plus 16 worker processes, overlapping A4-A6

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 3.1 | `prefetch_rewards` | collects key/W/H/names of the 256 rows; builds jobs `(text, truncated, names, CacheRef(key, root), W, H)` | CPU main thread | – | 0 | `fast.py:366, 369-376`; `src/objdel/train/rewards.py:113-125`; `serve.py:83-108` |
| 3.2 | submit | `ThreadPoolExecutor(1).submit(score_batch, jobs, pool)` | CPU | – | 0 | `rewards.py:119-121` |
| 3.3 | parse ×256 (prefetch thread, holds the GIL) | `text.strip`, compiled-regex `fullmatch`, `json.loads`, caps, `_ring` → int64 numpy arrays (MEASURED-MAC ≈2 ms for 256 target texts) | CPU | – | 0 | `src/objdel/reward/score.py:140`; `src/objdel/geometry/text.py:102-123` |
| 3.4 | `to_masks` packing | ring/role/shape_start lists, `polycore.pack` (`np.concatenate`), dims array | CPU | – | 0 | `src/objdel/geometry/shape.py:101-109`; `src/objdel/checker/polycore.py:55-60` |
| 3.5 | `_gpucheck.check_polygons`: binding | input checks: `check_offsets`, `check_dims`, `check_roles` per shape, then **releases the GIL** | CPU | – | 0 | `cuda/gpucheck_py.cpp:52-78` |
| 3.6 | engine host side | `word_start` per shape, `plan()` ×2 (256 ≥ kSplitMin 64, so two halves on two lanes), `Block::get`: grow-only arena and **pinned** staging (a cudaFree/cudaMalloc pair only when a buffer grows, then a device-wide sync; rare after the first steps), `pack()` copies into pinned staging | CPU | – | 0 per call (growth: device sync, rare) | `cuda/gpucheck.cpp:135-154, 1031-1070, 1109-1130` |
| 3.7 | engine device side | per lane: 1 async H2D, memset, 5 kernels (Precision, `raster_kernel` (one warp per row, ballot-packed words), CleanRings, EdgePairs, ShapeTotals), 1 async D2H. Lowest-priority, non-blocking streams | **GPU** | H2D ≈ tens of kB (vertices int32 + offsets). D2H ≈ 27 kB defects + ≈160 kB packed masks (50×50 → 400 B, 100×100 → 1.6 kB each) | 2 `cudaStreamSynchronize` on the checker's own streams + 2 more in `Drain` (they do not block the trainer's stream) | `gpucheck.cpp:63-102, 1073-1090, 1131-1139, 181-189` |
| 3.8 | result copy-out | pinned → `std::vector` copies, GIL retaken, wrapped zero-copy as numpy | CPU | – | 0 | `gpucheck.cpp:1133-1138`; `gpucheck_py.cpp:79-84` |
| 3.9 | per-object unpack | `unpack_rows`: `np.unpackbits` to bool H×W, defects dict, `Geometry` ×256 | CPU | – | 0 | `shape.py:86-90, 112-117` |
| 3.10 | pool dispatch | `pool.map(_score_checked, 32 chunks of 8)`. The queue feeder thread pickles each chunk: CacheRef as (key, root), **masks as bool arrays**, ≈0.8 MB per step (MEASURED-MAC) | CPU (main-process threads, need the GIL) | – | 0 | `score.py:178-183`; `serve.py:102-103` |
| 3.11 | worker: `rater_answers` | joint rater topology answers from cached `rater_sigs`/`voters` (pure Python; precomputable) | **CPU (worker)** | – | 0 | `score.py:102-113, 149` |
| 3.12 | worker: `score_mask` per object | `cv2.boundingRect` ROI of pred ∪ union box; `np.count_nonzero` area; **`topology.signature`** (`_crop_pad`: `cv2.boundingRect` + `cv2.copyMakeBorder`; `cv2.connectedComponentsWithAlgorithm` 8-conn foreground and 4-conn background; bincount/unique/lexsort; `size_filter` relabels twice when a component or hole is below min_area); `same_sig` lookup; **`_iou`** (`count_nonzero(a & b)`) per same-signature rater or against the vote; **`band_depth`** (`mask & ~union`, `core & ~mask`, gathers on the cached EDT maps, float64 sums); F5 formula; **`area_guard`**; topo_freq. MEASURED-MAC ≈7.5 ms for 256 serially | **CPU (worker)** | – | 0 | `score.py:61-94`; `src/objdel/checker/topology.py:75-77, 96-125, 139-163`; `src/objdel/reward/cache.py:130-136` |
| 3.13 | worker: `_score_one` aggregation | legal AND, defect merge, **joint topology** share, mean/min over objects, **`sig_id` = CRC32 of `repr(sigs)`** | **CPU (worker)** | – | 0 | `score.py:97-99, 147-164` |
| 3.14 | results back | dict per answer pickled back, unpickled by the executor's manager thread | CPU | – | 0 | `score.py:183` |
| 3.15 | `score_batch` assembly | float64 numpy arrays for 9 numeric keys + error/defects lists | CPU | – | 0 | `score.py:184-185` |

Whole prefetch wall time: ESTIMATE 10-50 ms. MEASURED-MAC serial CPU work for 256 answers is ≈15 ms (parse 2 ms,
CPU-build checker 3.5 ms, scoring 10.8 ms) plus IPC. It finishes while the main thread is still in A4/A5. In practice
it overlaps the CPU preprocessing more than the reference forward, which the `fast.py:369-372` docstring names.

### A4. Back on the main thread: tensors and the second processor pass (inside `time/batch_s`)

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 4.1 | pad prompt/completion ids and masks, sampling logps | Python lists → `torch.tensor` ×256 each → `pad` → `.to(device)` | CPU → GPU | H2D prompt_ids 0.67 MB, prompt_mask 0.67 MB, completion_ids 3.1 MB, completion_mask 3.1 MB, sampling_per_token_logps 1.6 MB (NaN mapping in Python first) | 5 (blocking `.to`) | `grpo_trainer.py:2547-2582` |
| 4.2 | `mask_truncated_completions` | `torch.tensor(is_truncated, device)` | CPU → GPU | 256 B | 1 | `grpo_trainer.py:2604-2609` |
| 4.3 | `num_items_in_batch`, concatenations | GPU sums and cats | GPU | – | 0 | `grpo_trainer.py:2613-2621` |
| 4.4 | **processor pass #2** | `apply_chat_template` (text only) ×256, then `processing_class(images=256, text=..., padding=True)`: the **full image preprocessing again** (MEASURED-MAC 1.0-1.5 s) | **CPU** | – | 0 | `grpo_trainer.py:2633-2645` |
| 4.5 | `Trainer._prepare_inputs(prompt_inputs)` | blocking `.to(device)` of every key, **pageable** (`configs.NO_SYNC`'s `non_blocking`/pin settings reach only accelerate's DataLoader, not this path) | CPU → GPU | **H2D 1.611 GB pixel_values fp32** + input_ids 0.67 MB + attention_mask 0.67 MB (both moved, then dropped) + mm_token_type_ids 0.67 MB + grid 6 kB | 5 | `grpo_trainer.py:2646-2647`; `transformers/trainer.py:2341-2356` |
| 4.6 | extend `mm_token_type_ids` | zeros for the completion part | GPU | – | 0 | `grpo_trainer.py:2686-2698` |

### A5. Old-logprob pass: skipped (C10)

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 5.1 | `disable_gradient_checkpointing` | toggles checkpointing off and back on (module walks) | CPU | – | 0 | `grpo_trainer.py:2750`; `trl/models/utils.py:380-397` |
| 5.2 | old logps | **not run**: GA 16 % (spg 16 × 1) == 0 and correction was set False by `FastGRPO` | – | – | – | `grpo_trainer.py:2758-2775`; `fast.py:276-294` |
| 5.3 | TRL IS ratio block | **not run** (correction False during this call) | – | – | – | `grpo_trainer.py:2777-2824` |

### A6. Reference pass (beta = 1e-3, PEFT: adapters disabled), 256 rows in 16 chunks of 16, no grad

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 6.1 | `use_adapter(model, None)` | disables the LoRA adapters | CPU | – | 0 | `grpo_trainer.py:2844-2856`; `trl/trainer/utils.py:1432` |
| 6.2 | `FastForward`: `rows = grid.prod(-1).tolist()` | pixel rows per image | GPU → CPU | 2 kB | 1 | `fast.py:60` |
| 6.3 | **RopeCache.positions** (all 256 rows miss) | `keys()`: `input_ids.cpu()`, `attention_mask.cpu()`, `grid.tolist()`, `tobytes` keys. Misses: GPU gather, then `.cpu()` of ids, mm_token_type_ids, mask, grid. **`get_rope_index` on the host copies** (per-row Python groupby/arange, MEASURED-MAC 0.03-0.05 s for 256 × 1861). `pos.to(device)` | CPU (+ D2H/H2D) | **D2H ≈19 MB** (5 × 3.8 MB int64), **H2D 11.4 MB** (3×256×1861 int64) | 8 | `fast.py:176-196`; `transformers/models/qwen3_vl/modeling_qwen3_vl.py:925-1007` |
| 6.4 | per chunk: `window()` | first real prompt column and longest answer, `.tolist()` | GPU → CPU | 8 B | 1 ×16 | `fast.py:32-43, 69` |
| 6.5 | per chunk: slicing | ids/mask/mm-type columns, pixel rows, grid; `rope.pending` set | CPU (index math), GPU views | – | 0 | `fast.py:70-80` |
| 6.6 | per chunk: TRL `_get_per_token_logps...` | `cum_rows[...]` `.item()` ×2, labels via masked_fill | GPU → CPU | – | 2 ×16 | `grpo_trainer.py:1517-1523, 1562-1563` |
| 6.7 | per chunk: model forward under autocast | **RopeCache pre-hook** injects `position_ids` (TRL's per-forward `get_rope_index` is skipped) | GPU | – | 0 | `fast.py:170-173`; `modeling_qwen3_vl.py:1228-1237` |
| 6.8 | `get_image_features` | `pixel_values.type(bf16)` (16 images: 100 MB read, 50 MB write); `split_sizes .tolist()` | GPU | – | 1 ×16 | `modeling_qwen3_vl.py:1039, 1044` |
| 6.9 | **VisionCache** | `keys()`: CPU generator `randint(1536)` then `.to(device)`; integer hash of pixel rows on the GPU; `h.tolist()` + `grid_thw.tolist()`; `rows .tolist()`. **2 new images per chunk**, so the frozen vision tower encodes 16 times × 2 images per step (32 distinct images) | GPU (+ tiny H2D/D2H) | H2D 12 kB per chunk | 4 ×16, plus per encode: `grid_thw.tolist()` in vision helpers and `lengths.tolist()` in each of 24 SDPA vision blocks, so ≈26 ×16 (the attention path is UNVERIFIED) | `fast.py:100-146`; `transformers/vision_utils.py:112`; `modeling_qwen3_vl.py:299` |
| 6.10 | placeholder check | `torch_compilable_check(tensor)` → `torch._check_tensor_all_with` | GPU → CPU | – | 1 ×16 | `modeling_qwen3_vl.py:1078-1081`; `transformers/utils/import_utils.py:1902-1903` |
| 6.11 | `masked_scatter` of image embeds | – | GPU | – | UNVERIFIED (0-1) | `modeling_qwen3_vl.py:1190` |
| 6.12 | causal mask (SDPA) | `fast_all(padding_mask)` inside an `if` | GPU → CPU | – | 1 ×16 | `transformers/masking_utils.py:226-228, 271` |
| 6.13 | 28 decoder layers + DeepStack | layers on GPU; `_deepstack_process` on layers 0-2 indexes with a boolean mask (read + write, each a `nonzero`) | GPU | – | 6 ×16 | `modeling_qwen3_vl.py:821-855` |
| 6.14 | fused LM head (Triton) | `hidden_states[mask]`, `labels[mask]` (boolean indexing), chunked log-prob kernel, masked_scatter | GPU | – | 2 ×16 | `trl/trainer/utils.py:1569-1720` (1694) |
| 6.15 | `F.pad`, `torch.cat` | back to full width | GPU | – | 0 | `fast.py:88-94` |

Cost: not separately timed in R2. ESTIMATE 4-6 s, from MEASURED-P1 "old + ref logprobs 12.3 s" before C10 removed
half of it, and before C5 removed the per-forward `get_rope_index`.

### A7. Reward call and advantages (`time/rewards_s` = 0.0037 s, MEASURED-R2)

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 7.1 | `_calculate_rewards` setup | `torch.zeros((256, 8), device)`; reward kwargs lists built from 256 dicts | CPU/GPU | – | 0 | `grpo_trainer.py:1648-1664` |
| 7.2 | reward func 1 (`format_ok`) → `TRLReward._terms` | group check, rebuilds jobs, compares the signature with the prefetch, `future.result()` (no wait in practice) | CPU | – | 0 | `rewards.py:131-145` |
| 7.3 | **`legal_gate`, GATED products, `topology_reward`, `advantages_v4`** | numpy: reshape (32, 8); z-scores with floors; `_topology` (pairwise sig_id equality, `np.unique`, `loo_bias` in Python with `@cache`, log bracket); flat-group zeroing (`ptp`); `.tolist()` of 8 outputs (MEASURED-MAC 0.23 ms) | **CPU** | – | 0 | `rewards.py:146-149`; `src/objdel/reward/advantages.py:50-105` |
| 7.4 | `log_metric` / `log_extra` | error/defect shares (Python counts), completions-table columns | CPU | – | 0 | `rewards.py:150-160` |
| 7.5 | 8 reward funcs → device | `torch.tensor(list_of_256, device)` ×8; funcs 2-8 return the cached result | CPU → GPU | H2D 8 × 1 kB | 8 | `grpo_trainer.py:1693-1706` |
| 7.6 | all-NaN row check | `if torch.isnan(...).all(1).any()` | GPU → CPU | – | 1 | `grpo_trainer.py:1734` |
| 7.7 | `timed("rewards")` end | `cuda.synchronize` | – | – | 1 | `fast.py:378-380, 389` |
| 7.8 | sum_then_normalize | weighted sum, group mean/std, advantages, nan_to_num (A = advantages_v4 in fp32) | GPU | – | 0 (boolean `[mask] = nan` dispatches to masked_fill: UNVERIFIED) | `grpo_trainer.py:2880-2925` |

### A8. Metrics and logs inside `_generate_and_score_completions`

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 8.1 | per-reward-function mean/std | `.item()` ×2 ×8 | GPU → CPU | tiny | 16 | `grpo_trainer.py:2940-2944` |
| 8.2 | reward, reward_std, frac_reward_zero_std | `.item()` ×3 | GPU → CPU | tiny | 3 | `grpo_trainer.py:2945-2949` |
| 8.3 | completion logs | `rewards_per_func[:, i].tolist()` ×8, advantages `.tolist()` | GPU → CPU | 9 × 1 kB | 9 | `grpo_trainer.py:2951-2956` |
| 8.4 | pending `log_metric` values | per name (error/* and defect/*, about 16): `torch.tensor(device)` + gather + `.tolist()` | CPU ↔ GPU | tiny | ≈32 (H2D + read each) | `grpo_trainer.py:2968-2974` |
| 8.5 | `_logs["images"]` | keeps references to the 256 PIL images (log_multimodal default) | CPU | – | 0 | `grpo_trainer.py:2976-2977` |
| 8.6 | TRL IS metrics | **skipped** (C10; FastGRPO keeps device-side sums instead) | – | – | 0 | `grpo_trainer.py:2979-3010` |
| 8.7 | output dict | GPU tensors + `num_images` list | – | – | 0 | `grpo_trainer.py:3017-3058` |
| 8.8 | `timed("batch")` end | `cuda.synchronize` | – | – | 1 | `fast.py:280-281, 389` |

### A9. Buffering the generation batch (TRL `_prepare_inputs`, after `time/batch_s`)

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 9.1 | `split_pixel_values_by_grid` | `image_grid_thw.prod(-1).tolist()`, `torch.split` into 256 views | GPU → CPU | 2 kB | 1 | `grpo_trainer.py:1611`; `trl/trainer/utils.py:1158-1187` |
| 9.2 | **`by_length`** (replaces the shuffle) | real-token counts, stable sort on the GPU, `order.tolist()`, every tensor gathered on the GPU, lists reordered | GPU (+1 read) | 2 kB | 1 | `fast.py:223-236` |
| 9.3 | `split_tensor_dict` → 16 micro-batches | slicing | CPU/GPU views | – | 0 | `trl/trainer/utils.py:1052-1084` |
| 9.4 | `unsplit_pixel_values_by_grid` | `torch.cat` of 16 images' pixel rows per micro-batch: **D2D copy of all 1.61 GB** | GPU | D2D 1.61 GB | 0 | `grpo_trainer.py:1614`; `trl/trainer/utils.py:1224-1237` |

### A10. Per micro-batch (×16): loss, forward, backward

| # | Stage | What it does | Dev | H↔D | Syncs (per micro-batch) | Code |
|---|---|---|---|---|---|---|
| 10.1 | `training_step` → `_prepare_inputs` | returns `_buffered_inputs[_step % 16]` | CPU | – | 0 | `grpo_trainer.py:1615` |
| 10.2 | `FastGRPO._compute_loss` | marks `loss_inputs` for C10 (no `importance_sampling_ratio` in the inputs) | CPU | – | 0 | `fast.py:310-316` |
| 10.3 | TRL `_compute_loss` prep | cat ids/masks | GPU | – | 0 | `grpo_trainer.py:3141-3148` |
| 10.4 | `FastForward` | `rows .tolist()`; **RopeCache** `keys()`: `.cpu()` of ids and mask (16×1861 int64) + `grid.tolist()`, all hits, so a GPU `torch.stack`; `window()`; slicing | CPU + GPU | D2H ≈0.48 MB | 1 + 3 + 1 = 5 | `fast.py:57-80, 176-196` |
| 10.5 | TRL inner | `.item()` ×2, labels | – | – | 2 | `grpo_trainer.py:1523` |
| 10.6 | forward (autocast bf16, gradient checkpointing on the 28 decoder layers) | same as 6.8-6.14: pixel cast, VisionCache (all hits: `torch.cat` of cached embeds), placeholder check, `fast_all`, DeepStack ×3, fused LM head + entropy. LoRA layers run without PEFT's fp32 input cast (`no_lora_input_cast`, under autocast) | GPU | H2D 12 kB (VisionCache weights) | 4 (vision) + 1 (split) + 1 (check) + 1 (sdpa) + 6 (deepstack) + 2 (fused head) = 15 | `fast.py:213-220, 264-269`; `modeling_qwen3_vl.py`; `peft/tuners/tuners_utils.py:2220-2235` |
| 10.7 | **C10 `set_is_ratio`** | IS ratio (`is_ratio`: exp of the masked logprob gap, sequence_mask, clip 3.0) from the detached update logprobs; 7 IS statistics accumulated **on the device** | GPU | – | 0 | `fast.py:199-210, 318-345` |
| 10.8 | loss | ratio, clip (0.2/0.28), KL to ref, × IS ratio, dr_grpo sum / (16 × 1536) / GA | GPU | – | 0 | `grpo_trainer.py:3178-3311` |
| 10.9 | loss metrics | `kl`, `entropy`, clip low/high/region (`global_masked_mean(...).item()` ×5), `low_min`, `high_max` (`.item()` ×2). **All happen before backward is launched**, so the GPU drains before the backward kernels are queued | GPU → CPU | tiny | 7 | `grpo_trainer.py:3376-3405` |
| 10.10 | `accelerator.backward(loss)` | autograd on its device thread: fused-head backward, recompute of the 28 checkpointed layers, LoRA grads (fp32) | GPU (launched from the autograd CPU thread) | – | 0 known (DeepStack and mask creation sit outside the checkpointed layers, so they are not re-run) | `trainer.py:2098` |
| 10.11 | `_tr_loss += loss` | GPU add (`logging_nan_inf_filter` off, so no isnan sync) | GPU | – | 0 | `trainer.py:1890-1903`; `configs.py:145` |
| 10.12 | `floating_point_ops` | returns 0 (no module walk) | CPU | – | 0 | `fast.py:414-415` |

Per micro-batch: ≈29 syncs. Cost: `step_time − time/batch_s` = 14.4 s mean for A9 + 16 × A10 (MEASURED-R2).

### A11. Optimizer step and logging (outside `step_time`)

| # | Stage | What it does | Dev | H↔D | Syncs | Code |
|---|---|---|---|---|---|---|
| 11.1 | `_clip_grad_norm` | foreach norm over the LoRA grads, scale (returns a tensor) | GPU | – | 0 | `trainer.py:1911, 2738-2744` |
| 11.2 | `optimizer.step()` | fused AdamW over ≈392 LoRA tensors | GPU | – | 0 | `trainer.py:1915` |
| 11.3 | scheduler, `zero_grad`, callbacks (DefaultFlow, Progress) | constant LR, grads set to None | CPU | – | 0 | `trainer.py:1919-1929` |
| 11.4 | `_maybe_log_save_evaluate` | `tr_loss .item()` (waits for the last backward + AdamW), `grad_norm.item()` | GPU → CPU | tiny | 2 | `trainer.py:2229, 2236` |
| 11.5 | `FastGRPO.log` | `is_stats.tolist()` (one read for the 7 IS metrics of the step) | GPU → CPU | 28 B | 1 | `fast.py:347-359` |
| 11.6 | TRL `log` | averages `_metrics`, prints (log_completions False: no table or parquet) | CPU | – | 0 | `grpo_trainer.py:3424-3437` |

---

## 1. Everything on the CPU, ranked

Ranked by impact on step time; the last column says whether the item breaks the rule "checking and scoring stay on
the GPU".

| Rank | CPU work | Time per step | Rule? | Evidence |
|---|---|---|---|---|
| 1 | **Dataloader**: 16 × 256-row fetches = 4,096 PNG decodes + arrow→Python; 15 of the 16 batches are discarded | MEASURED-P1 12.96 s; R2 ESTIMATE 18-22 s of the MEASURED 22.6 s outside `training_step` | no | P1 `trace/ops.txt` (`enumerate(DataLoader)#_SingleProcessDataLoaderIter._next_data`, 16 calls, 810 ms avg); P1 `fast/smi.csv` ≈14.5 s idle runs per step; MEASURED-MAC 0.75-0.88 s per fetch |
| 2 | **Qwen3-VL image preprocessing, pass #1** (`_tokenize_prompts`; output pixel tensor unused) | MEASURED-MAC 0.75-1.31 s; R2 ESTIMATE 1-2.5 s | no | `grpo_trainer.py:1806` |
| 3 | **Image preprocessing, pass #2** (forward inputs) | MEASURED-MAC 1.0-1.5 s; R2 ESTIMATE 1.5-3 s | no | `grpo_trainer.py:2645` |
| 4 | vLLM host side: scheduler, input prep, detokenization, logprob objects, mm processing/hashing of 256 PIL images, per-decode-step sync. The decode tail (one long answer, ≈5.6 ms per step) is probably launch/host-bound at small batch | part of the 5.6 ms × L_max slope (MEASURED-R2) and the ≈4.1 s intercept | no | UNVERIFIED split |
| 5 | Host side of ≈1,200 trainer syncs (`.item()`, `.tolist()`, `.cpu()`, `nonzero`) and kernel launch overhead (P1: 798,774 launches/step before C5/C10) | not separable from GPU time; MEASURED-P1 "update is launch- and sync-bound", GPU busy 70 % | no | P1 trace |
| 6 | Weight sync loop: ≈700 `load_weights` Python calls + PEFT merge/unmerge dispatch | UNVERIFIED (ESTIMATE 0.1-1 s) | no | `vllm_generation.py:520-522` |
| 7 | `extract_logprobs` (≈32 k dicts), `one_pad` checks, `batch_decode`, padding lists into tensors, NaN mapping | ESTIMATE 0.05-0.3 s | no | §A2/A4 |
| 8 | **RopeCache `get_rope_index` on host copies** (once per step for 256 rows) + keys | MEASURED-MAC 0.03-0.05 s | no (position ids, not scoring) | `fast.py:189-193` |
| 9 | **Answer parsing** (regex + `json.loads` + caps) ×256 | MEASURED-MAC ≈2 ms | **yes: part of "checking"; CPU** | `text.py:102-123` |
| 10 | **Checker host side**: packing, binding checks, staging copies, `np.unpackbits` unpack | ≈1-3 ms (ESTIMATE) | **yes: host half of the checker** | `shape.py:101-117`; `gpucheck_py.cpp:57-71`; `gpucheck.cpp:1054-1070` |
| 11 | **Mask scoring in 16 workers**: IoU_tc, signature/topology (OpenCV CCL), band depth / F5, area guard, topo_freq, rater_answers, joint topology, sig_id | MEASURED-MAC 10.8 ms serial for 256 | **yes: the scoring itself, on the CPU** | `score.py:61-164`; `topology.py:96-163`; `cache.py:130-136` |
| 12 | **Pickling masks to and results from the workers** | ≈0.8 MB + ≈50-80 kB per step; ESTIMATE a few ms | **yes (only exists because scoring is on the CPU)** | `score.py:182-183` |
| 13 | **`advantages_v4`, `topology_reward`, `legal_gate`, GATED products, error/defect shares** | MEASURED-MAC 0.23 ms | **yes: advantage/topology scoring on the CPU** | `advantages.py:50-105`; `rewards.py:146-155` |
| 14 | Python glue: TRL loops, `deepcopy` (0.017 s MAC), message building, metric dicts, `model.train()` | small | no | – |

Items 9-13 together are all of "checking and scoring" except the checker kernels. They cost ≈15-50 ms per step,
fully overlapped (§3), but they run on the CPU.

## 2. Every host↔device transfer and sync per step

### 2a. Transfers (L = 1536, W = 1861)

| Direction | What | Size per step | Pinned? | Code |
|---|---|---|---|---|
| H2D | `pixel_values` fp32 (processor #2) | **1.611 GB** | no (blocking, pageable) | `grpo_trainer.py:2646` → `trainer.py:2356` |
| H2D | processor input_ids, attention_mask (moved, then dropped), mm_token_type_ids, image_grid_thw | ≈2 MB | no | same |
| H2D | prompt_ids, prompt_mask, completion_ids, completion_mask, sampling logps, 2 × is_truncated | ≈9.2 MB | no | `grpo_trainer.py:2547-2609` |
| H2D | RopeCache positions (ref pass) | 11.4 MB | no | `fast.py:193` |
| H2D | VisionCache hash weights (`w.to(device)`, 32 forwards) | 32 × 12 kB = 0.4 MB | no | `fast.py:122` |
| H2D | checker inputs (vertices, offsets, dims) | ≈tens of kB | **yes (staging)**, async | `gpucheck.cpp:1079` |
| H2D | rewards (8 lists of 256) and metric tensors (≈19 small) | ≈10 kB | no | `grpo_trainer.py:1706, 2971`; `_generate` metrics |
| H2D | vLLM: image tensors and per-step scheduler inputs | UNVERIFIED | UNVERIFIED | vLLM |
| D2H | checker outputs: defects + packed even-odd masks | ≈190 kB | **yes**, async | `gpucheck.cpp:1089` |
| D2H | RopeCache keys/misses (ref pass) | ≈19 MB | – | `fast.py:179, 190-192` |
| D2H | RopeCache keys (update, 16 ×) | ≈7.6 MB | – | `fast.py:179` |
| D2H | scalars and lists: window, `.item()`, `.tolist()`, metric reads | ≈30 kB | – | §A |
| D2H | vLLM sampled ids + logprobs per decode step | UNVERIFIED (small) | – | vLLM |
| D2D | vLLM weight sync (whole model, every step) | ≈4.3 GB (ESTIMATE) | – | `vllm_generation.py:522` |
| D2D | `unsplit_pixel_values_by_grid` cat | 1.61 GB | – | `trl/trainer/utils.py:1231` |
| D2D | pixel fp32 → bf16 cast per forward (32 ×) | 32 × (100 MB read + 50 MB write) | – | `modeling_qwen3_vl.py:1039` |
| host only | processor #1's 1.61 GB pixel tensor (never moved), processor #2's CPU copy, masks to the pool (0.8 MB pickled) | – | – | – |

Totals: H2D ≈ 1.63 GB per step, almost all of it `pixel_values`; it is pageable and blocking (at 10-25 GB/s, ESTIMATE
0.07-0.2 s). D2H ≈ 27 MB. The reward path itself moves only ≈200 kB down and ≈10 kB up.

### 2b. Syncs (ESTIMATE from code; current code not traced on a GPU)

| Phase | Count | From |
|---|---|---|
| Generation-phase metrics | ≈11 | `grpo_trainer.py:2370-2405` (`.item()` + `torch.tensor(device)`) |
| `empty_cache` before wake-up | 1 (device-wide) | `vllm_generation.py:501` |
| vLLM decode loop | ≈1 per decode step, ≈1,350-1,536 | UNVERIFIED |
| Timers (`generate`, `rewards`, `batch`) | 3 (`cuda.synchronize`) | `fast.py:389` |
| Post-generation tensors + processor #2 H2D | 11 | §A4 |
| Reference pass | 1 + 8 + 16 × (1 + 2 + 1 + 4 + 1 + 1 + 6 + 2) ≈ 297, plus ≈26 × 16 ≈ 416 in the 16 vision-encode calls | §A6 |
| Reward call | 10 | §A7 |
| Reward/log metrics | ≈60 | §A8 |
| Buffering | 2 | §A9 |
| Update | 16 × ≈29 ≈ 464 | §A10 |
| Optimizer and logging | 3 | §A11 |
| Checker | 4 `cudaStreamSynchronize` on its own low-priority streams (they do not block the trainer stream); `cudaFree`/`cudaMalloc` only when a buffer grows | `gpucheck.cpp:102, 186, 1135` |
| **Total, trainer side** | **≈1,280** (+ vLLM's ≈1,500) | |

For comparison, MEASURED-P1 (one traced step, older code): 10,151 `cudaStreamSynchronize`, 7,233
`_local_scalar_dense`, 2,904 `nonzero`. C5 removed most of the per-forward `get_rope_index` reads since then; there is
no trace of the current code.

## 3. Overlapped vs critical path

Threads and processes in a step:
- main thread (trainer; vLLM's engine in-process, UNVERIFIED);
- the `reward-prefetch` thread (`rewards.py:119-121`);
- the ProcessPoolExecutor manager thread and the multiprocessing queue feeder thread (pickling);
- 16 forked scoring workers;
- the autograd engine thread (backward);
- torch intra-op threads (used by the processors' CPU torch ops).

CUDA streams:
- torch's default stream (the trainer);
- 2 lowest-priority non-blocking streams of the checker (`gpucheck.cpp:70-76`);
- vLLM's own streams (UNVERIFIED).

| Work | Overlapped? | Evidence |
|---|---|---|
| Scoring (parse, checker, unpack, pool scoring, IPC) | **Yes**. It starts right after generation (`fast.py:366`) and finishes in tens of ms, while the main thread does A4 (processor #2, 1+ s) and A6 (ref pass). The checker kernels run on low-priority streams while the GPU is mostly idle (A4 is CPU work). | `time/rewards_s` 0.0037 s mean over 20 steps (MEASURED-R2), which is the residual wait + numpy + 8 copies + sync; MEASURED-P1 "reward call waits 3 ms" |
| `advantages_v4` etc. | No (main thread, in the reward call), but ≈0.2 ms | MEASURED-MAC |
| Dataloader (16 fetches) | **No.** Critical path, GPU idle; `num_workers=0`, no prefetch of the next generation batch | P1 smi: ≈14.5 s idle gaps per step (0.5 s samples); P1 trace 12.96 s |
| Processor pass #1 | No (inside generate, before weight sync) | `grpo_trainer.py:2315` |
| Processor pass #2 + 1.6 GB H2D | No (only the reward prefetch overlaps it) | `grpo_trainer.py:2645-2646` |
| vLLM generation vs training | No (colocated, sequential: sleep/wake around each generate) | `vllm_generation.py:719, 772` |
| Ref pass | No; the scoring has already finished | – |
| Update forward vs backward | No: the 7 `.item()` calls at the end of `_compute_loss` drain the GPU before backward is queued, and the 5+15 syncs per forward keep the CPU from running ahead | `grpo_trainer.py:3376-3405` |
| Last micro-batch's backward + AdamW | Overlaps the CPU's optimizer/callback code until the `tr_loss.item()` in logging | `trainer.py:2229` |

## 4. Where the 58 s goes

All figures are seconds per optimizer step. The average is over 20 steps; train_runtime 1132 s / 20 = 56.6 s.
Step 1's tqdm interval was 58.4 s.

| Bucket | Step 1 | Mean | Source |
|---|---|---|---|
| **Wall per step** | 58.4 | 56.6 | MEASURED-R2 (tqdm, train_runtime) |
| **Outside `training_step`** (wall − TRL `step_time`) | 21.7 | 22.6 | MEASURED-R2 difference |
| ├ 16 dataloader fetches (4,096 PNG decodes, CPU, GPU idle) | – | **≈18-22** | ESTIMATE: MEASURED-P1 12.96 s on an EPYC 9554; R2's EPYC 7742 is slower; MEASURED-MAC 12-14 s for 16 fetches |
| └ grad clip, AdamW, zero_grad, callbacks, logging syncs, drain of the last backward | – | ≈0.5-3 | ESTIMATE |
| **`step_time`** (16 × `training_step`) | 36.7 | 34.0 | MEASURED-R2 |
| ├ `time/batch_s` (`_generate_and_score_completions`) | 23.1 | 19.6 | MEASURED-R2 |
| │ ├ `time/generate_s` | 13.7 | 11.7 | MEASURED-R2; fit 4.1 s + 5.6 ms × completions/max_length (r = 0.90); L_max hit 1536 in 13 of 20 steps |
| │ │ ├ fixed part ≈4.1 s: processor #1 (≈1-2.5), deepcopy (≈0.02-0.05), empty_cache + wake + merge + ≈700 load_weights + unmerge + prefix reset, KV wake, vLLM prefill + ViT, sleep, extract_logprobs/decode/metrics | | | ESTIMATE / UNVERIFIED split |
| │ │ └ decode tail ≈5.6 ms × L_max (≈8.6 s at 1536) | | | MEASURED-R2 slope |
| │ └ batch − generate | 9.35 | 7.9 | MEASURED-R2 difference |
| │   ├ processor #2 (CPU) + 1.6 GB pageable H2D + tensor padding | | ≈1.5-3.5 | ESTIMATE (MEASURED-MAC 1.0-1.5 s for the processor) |
| │   ├ ref pass (16 no-grad chunks, RopeCache 0.03-0.05 s CPU, 16 vision encodes) | | ≈4-6 | ESTIMATE |
| │   ├ reward wait + advantages + 8 copies | 0.003 | 0.0037 | MEASURED-R2 |
| │   └ metrics syncs, timers | | ≈0.1 | ESTIMATE |
| └ buffering (split, by_length, 1.6 GB cat) + 16 × forward/backward (gradient checkpointing, ≈29 syncs each) | 13.6 | 14.4 | MEASURED-R2 difference |

In plain terms: about 40 % of the step is the CPU decoding the same 256 PNGs 16 times while the GPU idles. About 20 %
is vLLM generation, dominated by the single longest answer. About 25 % is the update and about 10 % the ref pass.
About 5 % is two CPU image-preprocessing passes. Scoring is about 0.05 %.

## 5. What the GPU-scoring move must cover

These are all the per-step functions that today run on the host after the checker kernels. Line numbers are as on
this checkout.

1. `objdel.geometry.shape.to_masks` lines 112-117 and `unpack_rows` 86-90. Today the packed masks are copied D2H
   (`gpucheck.cpp:1089, 1137`) and expanded with `np.unpackbits`. GPU scoring must consume the bit-packed words
   **on the device** (keep the lane arena or hand over a device pointer) instead of the `std::vector` copy-out.
2. `objdel.reward.score.score_batch` 167-185 and `_score_checked` 188-190: the `pool.map` over 32 chunks, and the
   pickling of masks and CacheRef (`serve.py:102-103`). This is replaced by one batched GPU scoring call. The numpy
   assembly of the 9 numeric keys also moves.
3. `objdel.reward.score._score_one` 147-164:
   - legal AND over objects and the defects merge (the defects come from the checker already);
   - joint topology share `joint = len(a.get(sigs)) / sum(...)`;
   - mean/min over objects (iou_tc, iou_min, f5, g_area);
   - `topo_any`;
   - `sig_id` 97-99. It is CRC32 of `repr(sigs)`, a host string hash: it needs a device-friendly exact id that
     `advantages._topology` can compare with equality, or it stays on the host on tiny signature tuples.
4. `objdel.reward.score.rater_answers` 102-113: per-item constant (rater signatures, voters). Precompute it once per
   item into device tables; there is no per-step compute.
5. `objdel.reward.score.score_mask` 75-94, every term:
   - ROI from `cv2.boundingRect` of the prediction ∪ `cache.box` (81-84);
   - area `count_nonzero` (85);
   - **`objdel.checker.topology.signature`** 157-163, with `instance_topology` 105-125 (8-connected foreground CCL,
     4-connected background CCL with the outside label, per-component areas, holes per component via the
     first-raster-pixel rule, hole areas, canonical `lexsort`), `size_filter` 139-154 (remove small components, fill
     small holes, then relabel), and `min_area` (cached);
   - `cache.same_sig` lookup (87);
   - **`_iou`** 61-64 per same-signature rater (max) or against the vote (88-91);
   - **`objdel.reward.cache.band_depth`** 130-136. Spill = sum over `mask & ~union` of (dist_out_union − ½); core
     miss = sum over `core & ~mask` of (dist_in_core − ½). These are float64 sums "in a fixed order (deterministic)",
     so a GPU reduction must keep a fixed order to stay bit-exact with the tests;
   - F5 formula (93);
   - **`area_guard`** 67-72;
   - `topo_freq` (94).
6. Device-resident expert caches. Each `ExpertCache` (`cache.py:210-231`) holds `raters`, `rater_areas`, `vote`,
   `n_vote`, `union`, `core`, `dist_out_union`, `dist_in_core`, `perimeter`, `t5`, `box`, `min_area`, `same_sig`,
   `voters` and `rater_sigs`. Today they live in host `serve.STORE`, forked into the workers (LIDC train ≈0.33 GB per
   the `dataset.cache_loader` docstring). They must be uploaded once and indexed by item key. `boundary`/`edt` are
   offline only, not per step.
7. `objdel.train.rewards.TRLReward._terms` 131-162: `legal_gate`, GATED products (147), `topology_reward` (148),
   `advantages_v4` (149), the error/defect share counts for `log_metric` (150-155), and the `.tolist()` to lists.
8. `objdel.reward.advantages`: `_groups` 50-56, `legal_gate` 67-70, `topology_valid` 73-76, `_topology` 79-87
   (pairwise sig equality, `np.unique` + `loo_bias` 59-64, a Python `@cache` exact binomial sum that can be a small
   device table for G = 8), `topology_reward` 90-92, `advantages_v4` 95-105 (per-group z-scores with floors, `ptp`
   flat-group zeroing).
9. TRL's reward hand-off: `_calculate_rewards` requires each reward function to return a **Python list**
   (`grpo_trainer.py:1693-1706`: `torch.tensor(list, device)`). Keeping rewards on the GPU end to end needs
   `FastGRPO._calculate_rewards` (`fast.py:378-380`) to fill `rewards_per_func` from device tensors directly.
   Otherwise there is a D2H + H2D round trip and 8 syncs.
10. Still on the CPU by nature, which the owner should decide on explicitly:
    - text parsing (`text.parse`, strict grammar + `json.loads`);
    - `polycore.pack`, the binding's input checks (`gpucheck_py.cpp:57-71`) and `pack()` into staging
      (`gpucheck.cpp:1054-1070`).
    These are the host half of "checking".
11. Outside the training step but in the same g_probe run: the probe scoring (`select_ckpt.score_items` →
    `_score` → `objdel.eval.metrics.metrics` → `score_many`, plus NSD/GED/soft Dice/HM-IoU etc.) runs **inside the
    pool workers** (`experiments/phase2/select_ckpt.py:172-176, 223-231`; `src/objdel/eval/metrics.py:247-260`).

## 6. Side findings (outside the strict CPU/GPU accounting; flagged for the owner)

1. **Probe-path workers call the CUDA checker.** `score_items` maps `_score` over the 16 forked workers, and
   `metrics` → `score_many` → `to_masks` → `_gpucheck.check_polygons` then runs **in each worker**.
   - The pool was forked before CUDA started, so each worker that runs a probe job creates **its own CUDA context**
     with its own arena and pinned buffers. The `Engine` is never freed (`gpucheck.cpp:173-176`).
   - That happens at the two step-0 probe draws, before training. Up to 16 extra contexts then stay resident on the
     5090 through all 20 training steps.
   - This contradicts the docstrings ("workers never touch CUDA": `score.py:169-171`, `rewards.py:91-93`), which hold
     only for the training reward path. UNVERIFIED on the box (no per-process `nvidia-smi` capture). Check with
     `nvidia-smi --query-compute-apps=pid,used_memory`.
2. `configs.NO_SYNC`'s comment says "non-blocking host-to-GPU copies (pinned memory is on)" (`configs.py:27-28,
   145`). In GRPO this reaches only accelerate's DataLoader, whose batches here hold PIL images and no tensors. Every
   real transfer, including the 1.6 GB `pixel_values`, goes through `Trainer._prepare_input`'s blocking
   `data.to(device)` (`transformers/trainer.py:2356`) or TRL's `.to(device=device)`.
3. **bf16 merge/unmerge drift.** Each weight sync does `W += delta` then `W -= delta.to(bf16)` on bf16 base weights
   (`peft/tuners/lora/layer.py:1007, 1036`). That is not exactly reversible, so the frozen base (and with it the
   reference model = base with adapters disabled) can drift by rounding over 20 syncs. UNVERIFIED magnitude.
4. VisionCache encodes 2 images per reference chunk (16 encode calls, each ≈26 syncs on the SDPA vision path),
   rather than encoding the 32 distinct images once.
5. vLLM gets each image as PIL and processes it again inside vLLM (UNVERIFIED), so per step the same 256 images pass
   through the Qwen3-VL image processor three times in total (2 × TRL + vLLM), plus 4,096 PNG decodes.

## UNVERIFIED list

- vLLM internals: in-process engine with `external_launcher`; mm processor cache dedupe of the 8 copies; encoder
  cache; per-decode-step D2H/sync; sleep/wake memory API syncs; `load_weights` per-call cost.
- R2's dataloader time: not logged; inferred from wall − step_time and P1's trace on a different CPU.
- Split of the ref pass vs processor #2 inside the 7.9 s, and of the fixed 4.1 s inside generate.
- Exact sync counts in `masked_scatter`, the SDPA vision path and TRL's boolean `[mask] = nan`. I counted standard
  PyTorch behavior (boolean indexing → `nonzero` → sync; `torch._check_tensor_all_with` → `.item()`), not a trace of
  the current code.
- That the ckpt's preprocessor config equals the base model's (patch 16, merge 2, temporal 2, which gives 1,024 rows
  per 512² image). The local `runs/r2/ckpt` copy is partial.
- Attention implementation in training is SDPA (no `attn_implementation` is set in the repo; P1's trace shows
  `_efficient_attention_*`).

Scripts used for the MEASURED-MAC numbers: `/tmp/cputrace/rew.py` and inline `python -I` snippets. No repo file was
edited.
