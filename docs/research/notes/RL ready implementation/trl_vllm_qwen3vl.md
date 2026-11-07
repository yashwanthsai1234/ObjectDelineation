# TRL 1.17.1 + vLLM 0.32.0 + Qwen3-VL-2B: GRPO data flow, colocate, sm_120

Date: 2026-11-06. Scope: G-probe on one RTX 5090. No code was edited. Nothing was run on a GPU.

## Source key

- `G` = installed TRL `trainer/grpo_trainer.py`. `V` = `generation/vllm_generation.py`. `U` = `trainer/utils.py`. `C` = `trainer/grpo_config.py`.
  All under `<repo>/.venv/lib/python3.12/site-packages/trl/` (TRL 1.17.1, `trl-1.17.1.dist-info`).
- `TF519`, `TF5201` = `src/transformers/models/qwen3_vl/modeling_qwen3_vl.py` at tags v5.19.0 and v5.20.1
  (raw.githubusercontent.com/huggingface/transformers/v5.19.0/... and /v5.20.1/...). The installed Mac copy is byte-identical to TF5201 (`diff -q`).
- `VLLM@0.32` = file at tag v0.32.0 (raw.githubusercontent.com/vllm-project/vllm/v0.32.0/...). Release notes: https://github.com/vllm-project/vllm/releases/tag/v0.32.0
- Line numbers for TF519 are about +4 from TF5201 after line 150. Both are given where they differ.
- UNVERIFIED = not confirmed in source or by a run.

## Top findings (read these first)

1. TRL 1.17.1 still sends processor-expanded `prompt_token_ids` plus the raw PIL image to vLLM in colocate mode (`V:725-727`). The fix PRs are open, not merged: https://github.com/huggingface/trl/pull/6718 and /pull/6719; issue https://github.com/huggingface/trl/issues/6294 is open. TRL `main` has the same code (`main` `vllm_generation.py:729`).
2. For Qwen3-VL, vLLM's image replacement target is the single `<|image_pad|>` token (`VLLM@0.32 vllm/model_executor/models/qwen3_vl.py:1655-1678`, target at ~1741). The planner applies one update per image (`vllm/multimodal/processing/processor.py:692-783`, `_apply_matches` 808-833). So a pre-expanded prompt of N pads very likely becomes 2N-1 pads inside vLLM. UNVERIFIED by a run. It is my reading of the source. Check it with a CPU/GPU smoke test (see "Preflight checks").
3. TRL issue #6401 (Qwen3-VL-2B, one image per prompt, colocate, sleep mode) reports a crash in `get_rope_index` and an importance ratio near 1e-6 on mixed-size images. It trained fine on uniform images. https://github.com/huggingface/trl/issues/6401 (closed as duplicate of #6294).
   My reading: the crash is trainer-side. `mm_token_type_ids` comes from a second processor call with `padding=True` (`G:2645`) and is never re-padded except under `pad_to_multiple_of` (`G:2686-2698`). If the processor pads right while `prompt_ids` pad left (`G:2551-2556`), the types misalign whenever prompt lengths differ. UNVERIFIED as the cause of #6401. The repo already sets `padding_side="left"` (`src/objdel/train/configs.py:137`), which avoids it. The tokenizer default is right (`transformers/tokenization_utils_base.py:980` in the installed 5.20.1).
4. The generation batch is padded once to the longest prompt and the longest completion of all 256 rows (`G:2550-2569`). Nothing is trimmed per micro-batch. Shuffle plus split happen at `G:1611-1614`. Micro-batches are plain row slices.
5. `position_ids` are never passed. The model inputs are only `input_ids`, `attention_mask`, optional `pixel_values`, `image_grid_thw`, `mm_token_type_ids` (+ `support_ids`, `labels`, `fused_lm_head`) (`G:1500-1567`). Qwen3-VL therefore runs `compute_3d_position_ids` -> `get_rope_index` (a Python per-sample loop) on every forward. Passing a `(3, B, L)` `position_ids` skips it (`TF5201:1228`).
6. TRL logs NO "share of masked sequences" metric. The only IS metrics are `sampling/importance_sampling_ratio/{min,mean,max}` and `sampling/sampling_logp_difference/{mean,max}` (`G:2987-3013`). In `sequence_mask` mode masked sequences have ratio 0.
7. Masked sequences still get the KL term (`G:3268-3272`: ratio multiplies the policy loss, then `beta * kl` is added after).
8. vLLM 0.32 handles `expandable_segments` itself. `CuMemAllocator.use_memory_pool` turns it off while allocating the sleep-mode pool and restores it after (`VLLM@0.32 vllm/device_allocator/cumem.py:381-428`). The repo note "sleep mode does not work with expandable_segments" is out of date for 0.32.0. Keeping it unset is still the safe choice. A runtime test is needed to call it safe (UNVERIFIED).
9. `vllm_max_model_length` must be set. Default is `None` -> model config value 262144 for Qwen3-VL-2B. vLLM raises at init if the KV cache cannot hold one sequence of `max_model_len` (`VLLM@0.32 vllm/v1/core/kv_cache_utils.py:1104-1140`). KV per token is 112 KiB, so 262144 tokens need 28 GiB.
10. vLLM 0.32.0 requires `torch==2.13.0` and `transformers>=5.16.1,<5.20.0` (`requirements/cuda.txt:7`, `requirements/common.txt:10` at v0.32.0). Box pins (torch 2.13.0, transformers 5.19.0) fit. The Mac venv (transformers 5.20.1) is outside vLLM's range, but vLLM is not installed there. TRL 1.17.1 accepts vLLM 0.22.0 to 0.32.x only (`trl/import_utils.py:107`).

## 1. TRL 1.17.1 GRPOTrainer data flow

### 1.1 Entry points

- `training_step` -> HF Trainer -> `_prepare_inputs` (`G:1591`) -> `compute_loss` (`G:3060`) -> `_compute_loss` (`G:3141`).
- `_prepare_inputs` in train mode (`G:1604-1616`):
  - `generate_every = steps_per_generation * num_iterations`. Generation runs when `_step % generate_every == 0` or the buffer is empty.
  - Order: `_generate_and_score_completions` (`G:1610`) -> `split_pixel_values_by_grid` (`G:1611`) -> `shuffle_sequence_dict` (`G:1612`) -> `split_tensor_dict(..., steps_per_generation)` (`G:1613`) -> `unsplit_pixel_values_by_grid` per chunk (`G:1614`).
  - The micro-batch returned is `_buffered_inputs[_step % steps_per_generation]` (`G:1616`).
- Your config: `steps_per_generation` defaults to `gradient_accumulation_steps` = 16 (`C:1115-1117`). `generation_batch_size = 16 * 1 * 16 = 256` (`C:1117`). One generation per optimizer step, 16 micro-batches of 16 rows.
- Dataloader batch is `per_device_train_batch_size * steps_per_generation` (`G:1321-1340`). The sampler repeats each prompt `num_generations` times (`RepeatSampler`, `G:1343-1380`). So the 256 rows already contain each of the 32 prompts 8 times, and images repeat with them. UNVERIFIED that the group is contiguous in order; the `RepeatSampler` `mini_repeat_count` semantics say yes.

### 1.2 `_generate_and_score_completions` (`G:2429-3058`)

Order of work, with lines:

1. Read `prompt`, `image` or `images` (`G:2518-2538`). Build multimodal chat messages (`G:2540-2553`).
2. `_generate` (`G:2286`): `_tokenize_prompts` (`G:1753`) uses `processing_class.apply_chat_template(..., tokenize=True, return_dict=True)` (`G:1819-1827`). That returns image-pad-expanded `input_ids` as Python lists. Then `_generate_single_turn` (`G:1832`) calls vLLM with these ids and the images (`G:1840-1853`). Logprobs are the top-1 only (`G:1855`).
3. Pad to tensors on the GPU (`G:2550-2569`):
   - `prompt_ids`, `prompt_mask`: left pad (`G:2551-2556`).
   - `completion_ids`, `completion_mask`: right pad (`G:2559-2568`).
   - `sampling_per_token_logps`: right pad with 0.0; NaN where vLLM returned `None` (`G:2569-2580`).
   - `sampling_mask` (candidate set per token): padded with -1, shape `(B, T, K)` (`G:2581-2590`). Only exists with `top_p<1`, `top_k>0` or `min_p>0` in colocate with vLLM >= 0.28 (`G:1155-1160`). Your `top_k` default is 0 and `top_p` is 1.0 (`C:538-547`), so it is `None` unless you set them.
   - `tool_mask`: right pad with 1 (`G:2591-2597`). Only with tools or an env mask.
   - `pad_to_multiple_of` is `None` on CUDA (`C:1056-1059`, "auto" -> 128 only on Neuron).
4. `mask_truncated_completions` (you set True): a completion whose last id is not an EOS or pad id gets `completion_mask` zeroed (`G:2599-2608`).
5. `num_items_in_batch = gather(loss_mask.sum()).sum()` (`G:2610-2611`). Used by `cispo`, `dapo`, `vespo` only (`G:3300-3305`). Your `dr_grpo` ignores it.
6. `prompt_completion_ids = cat(prompt_ids, completion_ids)` and `attention_mask` (`G:2614-2617`). `logits_to_keep = completion_ids.size(1)` (`G:2620`). `batch_size = per_device_train_batch_size` (`G:2621`).
7. Image tensors: `prompt_inputs = self.processing_class(images=images, text=prompts_text, padding=True, return_tensors="pt")` (`G:2645`). This is a second, full processor pass over all 256 images. Keys other than `input_ids`, `attention_mask` become `forward_kwargs` (`G:2647`): `pixel_values`, `image_grid_thw`, `mm_token_type_ids`.
8. `mm_token_type_ids` is extended with zeros for the completion width: final shape `(B, P_proc + C)` (`G:2685-2698`).
9. Old, IS and ref logprobs under `torch.no_grad()` (`G:2754-2860`, see section 2).
10. Rewards, then advantages (`G:2862-2935`).
11. Metrics (`G:2938-3016`), then the output dict (`G:3018-3057`).

### 1.3 Keys in the output dict (`G:3018-3057`)

Per-completion-token tensors, shape `(B, C)` unless noted:

| Key | Shape | When | Line |
|---|---|---|---|
| `completion_ids` | (B, C) | always | `G:3021` |
| `completion_mask` | (B, C) | always | `G:3022` |
| `old_per_token_logps` | (B, C) | always with vLLM + IS correction (`G:2756-2761`) | `G:3027` |
| `ref_per_token_logps` | (B, C) | `beta != 0` (yours 1e-3) | `G:3034-3035` |
| `sampling_per_token_logps` | (B, C) | with vLLM logprobs (TRL asks `logprobs=0`, `G:1199`) | `G:3030-3031` |
| `importance_sampling_ratio` | (B, 1) in `sequence_*` modes; (B, C) in `token_*` modes | vLLM + IS correction | `G:3028-3029`, shape note `G:2800-2802` |
| `sampling_mask` | (B, C, K) | only with top-k/top-p/min-p | `G:3032-3033` |
| `tool_mask` | (B, C) | tools or `env_mask` | `G:3055-3056` |
| `mm_token_type_ids` | (B, P+C) | VLM Qwen (spans prompt AND completion) | `G:3047-3048` |

Per-sequence and per-prompt keys:

| Key | Shape | Note | Line |
|---|---|---|---|
| `prompt_ids`, `prompt_mask` | (B, P) | left padded | `G:3019-3020` |
| `advantages` | (B,) | `_compute_loss` unsqueezes to (B, 1) (`G:3173-3176`) | `G:3023` |
| `num_items_in_batch` | 0-dim | copied to every chunk by `split_tensor_dict` (`U:1052-1084`) | `G:3024` |
| `pixel_values` | (sum of patch rows, feat_dim) packed | for Qwen, rows = prod(grid_thw) per image | `G:3037-3038` |
| `image_grid_thw` | (num_images, 3) | | `G:3039-3040` |
| `num_images` | Python list, len B | | `G:3051-3053` |
| `token_type_ids`, `image_position_ids`, `spatial_shapes`, `pixel_attention_mask`, `image_sizes`, `num_tiles` | other VLM families | not Qwen | `G:3041-3054` |

### 1.4 Padding, device move, shuffle, split

- Tensors go to the GPU inside `pad(...).to(device=device)` (`G:2550-2569`). The whole generation batch (256 rows) lives on the GPU.
- Pad length is the max over all 256 rows. No per-micro-batch trim exists in `_prepare_inputs`, `split_tensor_dict` or `_compute_loss`. `_compute_loss` rebuilds `input_ids = cat(prompt_ids, completion_ids)` and `logits_to_keep = completion_ids.size(1)` from the full padded width (`G:3143-3147`).
  Consequence: one 1536-token runaway completion in any of the 256 rows makes every micro-batch run at the longest completion width (cost estimate, not measured).
- `shuffle_sequence_dict` (`U:1087-1115`): one `torch.randperm(batch_size)`. Tensors with `ndim>=1` are indexed; Python lists are re-ordered; 0-dim tensors are passed through. So GRPO groups are scattered across micro-batches.
- `split_tensor_dict` (`U:1052-1084`): `chunk_size = first_tensor.shape[0] // num_chunks`, equal slices. Lists are sliced, 0-dim tensors copied.

### 1.5 Pixel values and grids for Qwen

- Before shuffle: `split_pixel_values_by_grid` (`U:1158-1221`). With `image_grid_thw` present (Qwen path, `U:1178-1190`): `lengths = image_grid_thw.prod(-1)`, raises if `sum(lengths) != pixel_values.size(0)`, then splits `pixel_values` and `image_grid_thw` into Python lists, one entry per sample, using `num_images`.
- After shuffle and split: `unsplit_pixel_values_by_grid` (`U:1224-1262`) concatenates the lists back in each chunk.
- Inside the forward (`_get_per_token_logps_and_entropies`, `G:1482-1577`): for each inner chunk it slices `pixel_values[row_start:row_end]` and `image_grid_thw[img_start:img_end]` using cumulative `rows_per_sample` (`G:1527-1535`). Zero-image rows are handled by `num_images` (`G:1523-1526`).
- The inner chunk size is `batch_size` = 16 (`G:1510`, `G:2621`), so each forward sees 16 sequences and 16 images.

### 1.6 `_compute_loss` (`G:3141-3414`)

- Re-runs the forward with grad over the 16-row micro-batch (`G:3151-3170`), entropy on (`compute_entropy=True`).
- `old_per_token_logps = inputs.get("old_per_token_logps")`; falls back to `per_token_logps.detach()` (`G:3182-3184`).
- Ratio, then per-token loss (`G:3199-3262`). `dr_grpo` normalises by `B * max_completion_length` (`G:3282-3286`).
- IS multiplier applied at `G:3268-3269`; KL added at `G:3271-3272`.
- Metrics logged: `kl`, `entropy`, `clip_ratio/low_mean|high_mean|region_mean|low_min|high_max`, `policy_loss` only if entropy bonus on (`G:3335`).

## 2. Where reference and old logprobs are computed

- All three scorings use `_get_per_token_logps_and_entropies` with `batch_size = per_device_train_batch_size` = 16 (`G:2621`).
- Old logprobs: `G:2762-2772`, over the WHOLE 256-row generation batch in 16 inner chunks, under `torch.no_grad()` and with gradient checkpointing disabled (`G:2750`).
  They are always computed with vLLM + IS correction on (`G:2756-2761`).
- Ref logprobs: `G:2825-2860`.
  - `beta != 0`: with a PEFT model carrying a pretrained adapter, TRL adds a frozen `ref` adapter copied from `default` (`G:526-548`). Ref logprobs run under `use_adapter(model, adapter_name="ref")` (`G:2843-2859`), same 256 rows in chunks of 16. Cost: a second full forward pass over 256 rows per step.
  - Requires the model be passed as a PeftModel for this. If the adapter is attached via `peft_config` instead, ref = adapters disabled = base model (`G:2845`, `use_adapter(None)`). That is NOT your SFT policy. UNVERIFIED which path the repo uses; check `GRPOTrainer(model=...)`.
- Each step therefore runs 3 forward-heavy passes per row: old (no grad), ref (no grad), train (grad, 16 micro-batches). Plus vLLM generation.
- `fused_lm_head`: TRL patches the base model's forward (`G:1115-1124`, `U:1569-1700`). `model(**inputs, labels=labels, fused_lm_head=True)` (`G:1567`) calls `self.base_model(*args, **kwargs)` (`U:1647`, `backbone_attr="base_model"` for transformers >= 5, `U:1625-1628`), takes `last_hidden_state[:, :-1]`, and computes log-probs and entropy in tiles with the Triton kernel (`trl/kernels/chunked_logprob.py`). The `(B, T, vocab)` logits are never built. The `lm_head` runs only on positions where `labels != -100`, i.e. the completion tokens (`G:1563-1566`).
  The Triton kernel on sm_120: UNVERIFIED. There is a torch fallback (`_use_triton`, `U:61`).
- `position_ids`: not in `model_inputs` (`G:1500-1567`). So `Qwen3VLModel.forward` computes them (section 5).
- `support_ids` (the `sampling_mask`) is padded to prompt+completion width with -1 (`G:1556-1561`). Only with top-k/p/min-p.
- Sampling temperature: the fused head uses `temperature=self.temperature` (`G:1116`), so logits are scaled like vLLM's.

## 3. vLLM colocate with Qwen3-VL in TRL 1.17.1

### 3.1 Config that TRL builds (`V:298-386`, `G:1169-1202`)

| Setting | Value | Source |
|---|---|---|
| `gpu_memory_utilization` | `vllm_gpu_memory_utilization` (yours 0.35; default 0.3) | `V:359`, `C:658-665` |
| `max_model_len` | `vllm_max_model_length` (default None -> 262144) | `V:360`, `C:666-672` |
| `max_num_seqs` | `per_device_train_batch_size * vllm_tensor_parallel_size * steps_per_generation` = 16*1*16 = **256** | `G:1186-1188` |
| `max_num_batched_tokens` | hard-coded **4096** (comment: "memory profiler is misled by 32768") | `V:367-368` |
| `enable_prefix_caching` | not passed -> vLLM default | `V:356-373`; default true for generate models: `VLLM@0.32 vllm/engine/arg_utils.py:2993-3031`, `vllm/config/cache.py:141` |
| `enable_sleep_mode` | `vllm_enable_sleep_mode` | `V:362` |
| `logprobs_mode` | `"processed_logprobs"` | `V:370` |
| `distributed_executor_backend` | `"external_launcher"`, `seed = process_index // tp` | `V:364-366` |
| sampling | `n=1`, `logprobs=0`, `max_tokens=max_completion_length` | `V:677-684`, `G:1199` |
| `mm_processor_kwargs`, `limit_mm_per_prompt`, `enforce_eager` | not passed | `V:356-373` |

- TRL sends all 256 rows (8 duplicates of each prompt) in ONE `llm.generate` call with `n=1` (`V:719-742`). The 8 duplicates share prompt KV blocks via prefix caching (default on). Whether the 8 duplicate images also hit the multimodal processor cache: UNVERIFIED.
- `max_num_batched_tokens=4096` is not configurable from `GRPOConfig`. For a multimodal model vLLM still sizes the encoder budget as `max(max_num_encoder_input_tokens, max_tokens_per_mm_item)` (`VLLM@0.32 vllm/v1/core/encoder_cache_manager.py:356-361`). So an image above 4096 tokens does not fail on that budget. UNVERIFIED for the chunked-MM path.
- Default Qwen3-VL image limits: patch 16, merge 2, `shortest_edge` 65536 px, `longest_edge` 16,777,216 px (https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct/raw/main/preprocessor_config.json). One token covers 32x32 px, so the processor allows up to 16,384 tokens per image. vLLM profiles the encoder at the model's max tokens per item. That one-time spike counts against the 0.35 budget. Size: UNVERIFIED. The repo sets no `size` for Qwen3-VL (`configs.py:136-137` only for smolvlm), so this applies.
- Observed side issue: vLLM V1 pre-allocates a 16384-token encoder cache for multimodal models, which broke a 16 GB card (https://github.com/vllm-project/vllm/issues/46625, closed stale; vLLM 0.23/0.24, RTX 5080). Not a 32 GB problem.

### 3.2 Image-token double expansion (the #6294 / #6401 family)

- TRL passes `{"prompt_token_ids": ids, "multi_modal_data": {"image": img}}` (`V:725-727`). `ids` come from `apply_chat_template(tokenize=True)` and already hold N `<|image_pad|>` tokens.
- Issue #6294 (SmolVLM): vLLM re-applies prompt updates to token-id prompts, so the first pre-expanded token is expanded again (855 -> 1672 tokens). https://github.com/huggingface/trl/issues/6294. Fix PRs #6718 (prompt images) and #6719 (tool images) are OPEN, stacked, not merged. #6300 was closed unmerged.
- The #6294 body says "Qwen-style VLMs expand into a single contiguous `<|image_pad|>` run" and so pass TRL's server-mode run check. It does not claim colocate is safe for Qwen.
- vLLM 0.32.0 for Qwen3-VL: `PromptReplacement(modality="image", target=[hf_processor.image_token_id], replacement=[image_token_id]*num_tokens)` (`VLLM@0.32 qwen3_vl.py:1672-1678, 1741-1745`). `_apply_prompt_updates` (`processor.py:1699-1727`) uses `_apply_token_matches_with_placeholders` (851-905). The queue consumes one match per image item (`_plan_prompt_updates_with`, 692-783). The code has no "already expanded" check. Reading: the first pad becomes N pads and the other N-1 original pads stay as text tokens. UNVERIFIED by a run.
  Effect if true: vLLM generates from a prompt that has N-1 extra `<|image_pad|>` text tokens after the image, while the trainer scores the clean prompt. This is a standing train/sampler mismatch. It would show in `sampling/sampling_logp_difference/mean` and in `importance_sampling_ratio/mean` drifting from 1.
- #6401 data point (TRL 1.8.0, vLLM 0.23, transformers 5.13, Qwen3-VL-2B): after switching to `llm.chat` with a custom `rollout_func`, the IS ratio mean was ~0.74 and 2000+ steps ran clean. On the stock path with mixed sizes it was ~1e-6 (about 14 nats/token apart). Uniform CLEVR 480x320 worked on the stock path.
- Workaround in that issue: `rollout_func` calling `trainer.vllm_generation.llm.chat(conversations, params)` so vLLM builds its own prompt ids. See https://github.com/huggingface/trl/issues/6401 (body).
- PR #6718 approach: keep expanded ids for the training forward and give vLLM un-expanded ids (tokenizer-only render, `add_special_tokens=False`). https://github.com/huggingface/trl/pull/6718.

### 3.3 Sleep mode and weight reload

- TRL sleeps at init: `llm.sleep(level=2)` (`V:376-379`) and after every generation (`V:771-773`). Level is fixed at 2. Level 2 discards weights and KV cache (`VLLM@0.32 vllm/device_allocator/sleep_mode_backend.py:127-137`: only `runtime` and `cudagraph` tags are offloaded to host at level 2).
- vLLM saves buffers to CPU before level-2 sleep and restores them on weight wake-up (`VLLM@0.32 vllm/v1/worker/gpu_worker.py:298-346`). Parameters are NOT restored (no `--sleep-preserve-parameter-names`), so TRL must push weights again.
- TRL does: `empty_cache()`, `llm.wake_up(tags=["weights"])` (`V:500-503`), loads every parameter (`V:520-522`), `reset_prefix_cache()` (`V:528`). Then at generate: `wake_up(tags=["kv_cache"])` (`V:719`).
  - Weight sync is skipped when `global_step == _last_loaded_step` (`G:1836-1841`). It also re-syncs inside `generate` if the engine is asleep (`V:602-604`, comment cites https://github.com/vllm-project/vllm/issues/29341).
  - TRL issue #5312 (weights not synced with sleep mode on) is closed (2026-07-28). https://github.com/huggingface/trl/issues/5312.
- vLLM 0.32 release notes on sleep and RL: sleep releases workspace scratch (up to 4.5 GiB/GPU, #59156) and KV-init state (~4.6 GiB/GPU, #59158); `sleep(level=2)` keeps FP8 KV scales (#57163); opt-in `sleep_mode_offload_cudagraph` (default False, `vllm/config/model.py`, #59160). https://github.com/vllm-project/vllm/releases/tag/v0.32.0
- Open vLLM bugs near sleep: #65325/#65324 (server only: requests queued during level-2 sleep run before `reload_weights`; waking in steps avoids it, TRL does) https://github.com/vllm-project/vllm/issues/65324 ; #65769 (`CuMemAllocator.discard` keeps the CPU backup) https://github.com/vllm-project/vllm/issues/65769. Neither is on the colocate path TRL uses. UNVERIFIED for host RAM growth; watch RSS.
- How much GPU memory vLLM keeps while asleep (cudagraph pool, NCCL buffers outside the pool): UNVERIFIED. Measure `torch.cuda.mem_get_info()` right after init and after the first `sleep`.

### 3.4 LoRA weight sync

- `_iter_named_params` (`V:446-490`): for a PEFT model it does `model.merge_adapter()` (`V:459`), yields each parameter with the PEFT prefix removed (`removeprefix("base_model.model.")`, `.replace(".base_layer", "")`, skipping names containing `model.prefix` ("lora_") and `original_module`) (`V:469-478`), then `model.unmerge_adapter()` (`V:480`). Each parameter goes to `model_runner.model.load_weights([(name, param)])` (`V:522`).
- Names match vLLM's loader: HF-5 names `model.visual.*`, `model.language_model.*` map via `Qwen3VLForConditionalGeneration.hf_to_vllm_mapper` (`VLLM@0.32 qwen3_vl.py:1969-1975`; `load_weights` at 3227-3229). Tied embeddings are handled at `qwen3_vl.py:1921`. The 2B model has `tie_word_embeddings: true` (https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct/raw/main/config.json).
- Merge and unmerge in bf16 can leave small rounding drift in the base weights over many steps. UNVERIFIED for this setup. A cheap check: compare base weights at step 0 and step N with adapters disabled.
- Which adapter gets merged: the active one (`default`). The frozen `ref` adapter is not merged. UNVERIFIED in PEFT source; only inferred from `use_adapter` use at `G:2843-2845`.
- Step 0: `_last_loaded_step = -1` (`G:1168`, `G:1203`), so weights are pushed before the first generation. The vLLM init loads the base checkpoint from `model.name_or_path` (`V:357`), so the SFT adapter reaches vLLM only through this first sync.

### 3.5 Memory sizing for a 2B VLM, 32 GB, shared with training

Facts:
- `gpu_memory_utilization` is a fraction of TOTAL device memory for this engine. It covers weights, non-torch memory, peak activation, and profiled CUDA graph memory; graph profiling has been default since v0.21.0 (`VLLM@0.32 gpu_worker.py:735-745, 790-805`). Startup fails if free memory is below `total * util` (`vllm/v1/worker/utils.py:716-770`).
- Weights: `model.safetensors` is 4,255,140,312 bytes, 2,127,532,032 BF16 parameters (https://huggingface.co/api/models/Qwen/Qwen3-VL-2B-Instruct?blobs=true).
- KV per token (full attention, bf16): 28 layers x 8 KV heads x 128 head_dim x 2 (K and V) x 2 bytes = 114,688 bytes (112 KiB), from config.json.

Estimates (not measured):
- 0.35 x 32 GB is about 11 GB. After weights (4.0 GiB), the encoder profile spike, the activation peak and graphs, KV is likely 4 to 6 GiB, about 37k to 55k tokens.
- One 600-token prompt plus 1536 completion tokens is about 2136 tokens, 240 MB. So roughly 15 to 25 such sequences fit at once, versus `max_num_seqs=256`. vLLM will queue the rest. Prefix caching lets the 8 samples of a prompt share the prompt blocks.
- Existing repo claim "0.35 on 32 GB: ~11 GB for weights and KV" matches this arithmetic (`experiments/phase2/perf/fix_test.md:37`).
- TRL's own docs and defaults: 0.3 default (`C:658`). I found no TRL guidance for a VLM on 32 GB. UNVERIFIED that 0.35 is optimal; raising it costs training headroom only while vLLM is awake. During the optimizer step vLLM sleeps and frees weights + KV.
- CUDA graph sizes: default `max_cudagraph_capture_size = min(max_num_seqs * 2, 512)` = 512 here (`VLLM@0.32 vllm/config/vllm.py:2297, 2354-2364`; pattern `[1,2,4] + range(8,256,8) + range(256,max+1,16)` per `vllm/config/compilation.py:695-712`). About 50 graphs. Fine for a 2B model.
- `vllm_max_model_length`: set it to max prompt tokens + `max_completion_length` (TRL help text, `C:666-672`). Default would be 262144 and fail the KV check (finding 9). Prompt tokens depend on image size (1 token per 32x32 px). Take the max over the dataset.
- `max_tokens` is no longer reserved from the context window in 0.32 (#57035, release notes). So a prompt near `max_model_len` with 1536 new tokens can be cut by the length limit. UNVERIFIED exact behaviour; keep `max_model_len` >= prompt + 1536.

### 3.6 `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`

- vLLM 0.32 reads `PYTORCH_ALLOC_CONF`, then `PYTORCH_CUDA_ALLOC_CONF`, then `PYTORCH_HIP_ALLOC_CONF` (`VLLM@0.32 vllm/device_allocator/alloc_conf.py:20-44`).
- Inside the sleep-mode memory pool it rewrites the whole config with `expandable_segments:False`, then restores it (`cumem.py:381-428`). Comment cites https://github.com/pytorch/pytorch/issues/147851 (expandable segments and memory pools are incompatible).
- Open vLLM PRs: #54351 (unified-name handling in CuMem and KV-transfer guards) https://github.com/vllm-project/vllm/pull/54351 ; #43923 / #60584 (custom all-reduce crashes with expandable segments; needs TP/DP > 1, not your single GPU).
- Not verified: whether toggling the allocator config at run time while the trainer already holds expandable-segment memory is safe in torch 2.13. UNVERIFIED. The repo choice (G-probe leaves it unset, `experiments/phase2/perf/implementation.md:145`) is the safe default. To enable it, run a 3-step smoke test on the box first.

## 4. vLLM 0.32 on RTX 5090 (sm_120), CUDA 13

- Wheel: PyPI default is CUDA 13.0; images `vllm/vllm-openai:v0.32.0` (CUDA 13.0) and `-ubuntu2404` (release notes, install table). Pins: `torch==2.13.0`, `flashinfer-python==0.7.0.post1` (`requirements/cuda.txt:7,13`). Python >= 3.11. Source builds need GCC 13 (#60259).
- Driver: all CUDA 13.x use minor-version compatibility with driver >= 580 (NVIDIA CUDA toolkit release notes, Table 3; table 2 maps 13.0 to R580). https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/index.html. A user ran an RTX 5090 with driver 580.178.04, torch 2.14.0+cu130 and vLLM main (https://github.com/vllm-project/vllm/issues/65769). I have no source for a minimum for the exact box driver. UNVERIFIED for your host; check `nvidia-smi`.
- Attention backend on sm_120 (compute capability major 12, dense model, no MLA): priority list is FLASH_ATTN, FLASHINFER, TRITON_ATTN, FLEX_ATTENTION, TURBOQUANT (`VLLM@0.32 vllm/platforms/cuda.py:83-175`, the `else` branch for non-SM100). So FLASH_ATTN is tried first. Which FA version it uses on sm_120: UNVERIFIED. Qwen3-VL's vision encoder picks its own backend: UNVERIFIED.
- FlashInfer sampler: `VLLM_USE_FLASHINFER_SAMPLER` defaults on (`vllm/envs.py:878-882`), but the `TopKTopPSampler` uses it only when `logprobs_mode` is not a processed mode (`vllm/v1/sample/ops/topk_topp_sampler.py:220-226`). TRL sets `processed_logprobs` (`V:370`), so the native (PyTorch/Triton) sampler is used. FlashInfer's JIT arch failure on sm_12x with old CUDA toolkits (#42393) therefore does not hit the sampler path. It can still hit FlashInfer attention if selected.
- Capture sizes: see 3.5 (512 cap, ~50 graphs).
- Known sm_120 issues I found (none on your Qwen3-VL-2B bf16 path):
  - #46625 Qwen3-VL-8B-FP8 on RTX 5080 hung at first generate; vLLM 0.23, FP8 only. https://github.com/vllm-project/vllm/issues/46625 (closed stale).
  - #55775 MTP + FlashInfer long-context IMA on RTX 5090 (NVFP4, speculative decoding). Not used here. https://github.com/vllm-project/vllm/issues/55775
  - #66193 (open, 2026-11-07) FlashInfer autotune peak cuts KV 2.4x under an explicit `--gpu-memory-utilization` on Blackwell Server (regression from #56492). Check the logged "Available KV cache memory" at startup. https://github.com/vllm-project/vllm/issues/66193
  - #65533 degenerate output on sm_120 for GLM-4.7-AWQ at 0.25.1. Quantised model, not yours.
  - #58624 MoE decode 15% slower on SM12x. Not applicable (dense).
- I found no open vLLM issue that names Qwen3-VL-2B bf16 on an RTX 5090 with 0.32.0. The search used `gh search issues` for "RTX 5090 Qwen3-VL", "sm_120", "SM120", "Blackwell sm120". Absence of reports is not proof. First run on the box is the test.
- Release 0.32 notes on Qwen3-VL: M-RoPE cache 4x smaller (#59853), device-side input normalisation for Qwen3-VL (#51289), per-video options with the multimodal cache (#55203).

## 5. transformers 5.19.0 vs 5.20.1: `modeling_qwen3_vl.py`

Full-file diff v5.19.0 vs v5.20.1 is 3 hunks only:

1. Vision rotary embedding `forward` (TF519:148-155 -> TF5201:147-151): the `maybe_autocast(enabled=False)` wrapper is dropped and `freqs = position_ids[..., None].float() * self.inv_freq.to(device=x.device, dtype=torch.float)`. Class at TF5201:106. Also the text rotary `Qwen3VLTextRotaryEmbedding.__init__(self, config, device=None)` (TF519:358) becomes `__init__(self, config)` (TF5201:354). Both are device/autocast handling only. Expected no effect on logprobs (elementwise float32 ops, no matmul). UNVERIFIED numerically.
2. Import `maybe_autocast` removed (TF519:49).
3. `_tp_plan = {"lm_head": "colwise_gather_output"}` added to `Qwen3VLForConditionalGeneration` (TF5201:1270). Only for tensor parallel.

Also changed in the 5.20.1 modular and config files (not modeling): `configuration_qwen3_vl.py` uses `sub_configs_defaults` with `SubConfigSpec` instead of `sub_configs` (modular diff lines 26, 193-217; 29-line diff of the config file). `processing_qwen3_vl.py` and `modeling_vision_utils`/`vision_utils.py`, `image_processing_qwen2_vl*.py`, `utils/generic.py` are identical across the two tags. `masking_utils.py`, `modeling_layers.py`, `processing_utils.py`, `sdpa_attention.py`, `modeling_flash_attention_utils.py` differ, but only in torch<2.6 gates, MTP, `block_sequence_ids` and tracing guards (read of the diffs; none touches `create_causal_mask` for plain causal text).
So the functions below are the same in both versions apart from the line offset (+4 in TF519).

### 5.1 Vision tower

- `Qwen3VLVisionModel.forward(self, hidden_states, grid_thw, **kwargs)` (TF5201:677-737, `@capture_outputs`). Input `hidden_states` is `pixel_values` of shape `(seq_len, feat)`; `grid_thw` is `(num_images, 3)`.
- Returns `BaseModelOutputWithDeepstackFeatures(last_hidden_state, pooler_output, deepstack_features)` (class at TF5201:64-70; return at ~733-737):
  - `last_hidden_state`: pre-merger hidden states.
  - `pooler_output`: merged embeddings `self.merger(hidden_states)`.
  - `deepstack_features`: list, one tensor per `deepstack_visual_indexes` layer. For Qwen3-VL-2B the indexes are `[5, 11, 17]` of 24 layers (config.json `vision_config`).
- Vision position ids and attention seqlens come from `transformers/vision_utils.py` (`get_vision_position_ids`, `get_vision_interpolation_indices_and_weights`, `get_vision_attention_seqlens`; imports at TF5201:46-50). Identical in both tags.

### 5.2 `get_image_features`

- `Qwen3VLModel.get_image_features(self, pixel_values, image_grid_thw=None, **kwargs)` (TF5201:1033-1049; TF519:1037). Casts `pixel_values.type(self.visual.dtype)`, calls `self.visual(pixel_values, grid_thw=image_grid_thw, return_dict=True, **kwargs)`.
- It returns the same output object, with `pooler_output` replaced by a TUPLE of per-image tensors (split by `image_grid_thw.prod(-1) // spatial_merge_size**2`) and `deepstack_features` as a list of tuples. It does not return one tensor.
- `Qwen3VLForConditionalGeneration.get_image_features` (TF5201:1295) delegates to `self.model.get_image_features` (same signature).
- In `Qwen3VLModel.forward` (TF5201:1143-1262): `mm_encoder_outputs["image"] = self.get_image_features(pixel_values, image_grid_thw, return_dict=True, **kwargs)`; embeds are `torch.cat(pooler_output)` and scattered into `inputs_embeds` with `masked_scatter`. Deepstack features go to `deepstack_visual_embeds` (cat per layer) with `visual_pos_masks`. The text model adds them after decoder layers 0..2 (`_deepstack_process`, `TextModel.forward` at TF5201:762-850).
- `get_placeholder_mask` raises if the number of `<|image_pad|>` tokens differs from the number of image feature rows (TF5201:1051). This is the guard that a token/feature mismatch hits.

### 5.3 Position ids

- Name and signature in v5.19.0 and v5.20.1 (same):
  - `Qwen3VLModel.get_rope_index(self, input_ids, mm_token_type_ids, image_grid_thw=None, video_grid_thw=None, attention_mask=None, **kwargs) -> (position_ids (3,B,L), mrope_position_deltas (B,1))` (TF5201:925; TF519:929).
  - `Qwen3VLModel.compute_3d_position_ids(self, input_ids, inputs_embeds, image_grid_thw=None, video_grid_thw=None, attention_mask=None, past_key_values=None, mm_token_type_ids=None)` (TF5201:1092; TF519:1096).
  - `get_vision_position_ids(start_position, grid_thw, temp_merge_size=1, spatial_merge_size=1, time_interval=1.0, device=None)` (TF5201:873).
- `get_rope_index` needs `mm_token_type_ids` (0 text, 1 image, 2 video). It loops over the batch in Python, takes the attention-masked tokens, groups them with `itertools.groupby` on `.tolist()` of the token types (a device sync per sample), and builds positions per modality group (TF5201:~960-1018). Left padding is handled through `attention_mask`.
- `compute_3d_position_ids` raises `ValueError` if `image_grid_thw` or `video_grid_thw` is passed without `mm_token_type_ids` and `input_ids` is not None (TF5201:1102-1108). This is why TRL must carry `mm_token_type_ids`.
- Does passing `position_ids` skip the rope computation? Yes. `if position_ids is None: position_ids = self.compute_3d_position_ids(...)` (TF5201:1228-1237). Same in `Qwen3VLTextModel.forward` (`if position_ids is None` at TF5201:794; 2D ids are expanded to 4 rows; 3D with `shape[0]==4` splits text row 0 off; a plain `(3, B, L)` is used as the mrope ids).
- The model keeps `self.rope_deltas` as module state; `compute_3d_position_ids` sets it (TF5201:1114-1124). That is stale state after generation. The code skips stale deltas when `input_ids` is given without cache (comment at TF5201:~1119-1121).
- Cost (ESTIMATE, not measured): with 3 scoring passes x 256 rows (+ ref) per step, `get_rope_index` runs about 1,000 sample iterations per step on CPU. A precomputed `position_ids` would remove that. TRL has no hook for it (`G:1500-1567`), so it needs a subclass or a patch; not done here.

## 6. Metric names in TRL 1.17.1

Completion lengths (`G:2388-2405`), logged in `_generate` and averaged over the logging window (`G:3424-3441`):
- `num_tokens` (`G:2388`, cumulative tokens seen)
- `completions/mean_length`, `completions/min_length`, `completions/max_length` (`G:2391-2393`). When `tool_mask` is set, only model tokens are counted (`G:2368-2370`).
- `completions/clipped_ratio` (`G:2399`): share whose last token is not an EOS/pad id. This is the truncation rate.
- `completions/mean_terminated_length`, `completions/min_terminated_length`, `completions/max_terminated_length` (`G:2403-2405`). They are 0 when no completion terminated (`G:2400-2402`).
- Tool metrics `tools/call_frequency`, `tools/failure_frequency` only with tools.

vLLM importance sampling (`G:2979-3016`), only when `use_vllm and vllm_importance_sampling_correction`:
- `sampling/sampling_logp_difference/mean` (`G:2987`) and `/max` (`G:2990`): mean and max of `|old_per_token_logps - sampling_per_token_logps|` over valid tokens (NaN vLLM logprobs excluded, `G:2979-2986`). This is the best direct drift signal.
- `sampling/importance_sampling_ratio/min`, `/mean`, `/max` (`G:3007-3013`). In `sequence_*` modes the stats are over one ratio per sequence (`G:2995-2998`); in `token_*` modes over valid tokens.
- No metric for the share of masked sequences exists in TRL 1.17.1 or in TRL main (grep of `main` `grpo_trainer.py` shows the same five names at lines 2987-3013). With `sequence_mask`, a masked sequence has ratio exactly 0 (`G:2811-2823`: `masked_fill(invalid_mis_mask, 0.0)`). So `importance_sampling_ratio/min == 0` means at least one masked sequence. If unmasked ratios average near 1, then masked share is about `1 - mean`. That is an approximation (UNVERIFIED); an exact count needs a subclass that logs `(ratio == 0).float().mean()` after `G:2823`.
- Defaults: `vllm_importance_sampling_mode="sequence_mask"`, `clip_max=3.0`, `clip_min=None` (`C:330-347`, `C:939-966`). With `clip_min=None` only ratios above 3.0 are masked (`G:2806-2817`). The ratio is `exp(sum(old - sampling))` over the sequence (`G:2783-2796`), so long completions have a wide spread.

Other loss-side metrics: `kl`, `entropy`, `clip_ratio/low_mean`, `clip_ratio/high_mean`, `clip_ratio/region_mean`, `clip_ratio/low_min`, `clip_ratio/high_max` (`G:3390-3405`); `reward`, `reward_std`, `frac_reward_zero_std`, `rewards/<func>/mean|std` (`G:2942-2949`).

## Preflight checks (cheap, before renting; none run here)

1. Assert that vLLM's returned `output.prompt_token_ids` length equals TRL's `prompt_ids` length for the same row. `V:745` already collects `all_prompt_ids` but TRL discards it (`G:1840`: `_, completion_ids, ...`). A length of 2N-1 pads confirms finding 2.
2. Assert `mm_token_type_ids == (prompt_completion_ids == image_pad_id)` on a batch with mixed prompt lengths (confirms the left-pad alignment, finding 3).
3. On step 0, log `sampling/sampling_logp_difference/mean` and `importance_sampling_ratio/mean`. Mixed-image-size data is the #6401 trigger.
4. Set `vllm_max_model_length` to the dataset max prompt tokens + 1536, and confirm vLLM logs "Available KV cache memory" at startup.
5. Record GPU memory right after vLLM init, after the first `sleep`, and during the first optimizer step.

## Not verified

- Whether vLLM 0.32.0 really builds 2N-1 pads for pre-expanded Qwen3-VL prompts (finding 2).
- The root cause of the #6401 crash (finding 3).
- Triton chunked-logprob kernel on sm_120.
- FlashAttention version and the ViT attention backend that vLLM picks on sm_120.
- GPU memory vLLM keeps while asleep; host RAM for `runtime` and `cudagraph` backups.
- Safety of toggling `expandable_segments` at runtime in the colocate process.
- Merge/unmerge bf16 drift across many steps.
- The box's driver version.
- KV and activation numbers in 3.5 (arithmetic from config, not measured).
