# TRL 1.17.1 GRPO + vLLM 0.32.0 colocate: internals that plan tasks touch

Scope: A6, A7, B2.0-B2.4, C0-C11, E5. Date 2026-11-09.

How the facts were checked:
- TRL: the installed `.venv` TRL 1.17.1 is byte-identical to tag `v1.17.1` for `vllm_generation.py`,
  `grpo_trainer.py` and `grpo_config.py` (I diffed them against raw.githubusercontent). So `.venv` line numbers equal
  the tag's line numbers.
- vLLM: vLLM is not on the Mac. I read the source at tag `v0.32.0` (released 2026-10-18; `v0.33.0rc3` exists, so
  do not read `main`).
- transformers: the installed 5.20.1 wheel.

Base URLs used below:
- `T` = https://github.com/huggingface/trl/blob/v1.17.1/trl/
- `V` = https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/
- `HF` = https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/ (line numbers are from the
  installed wheel)

## 1. One RL step in colocate mode (our config: sleep on, IS on, dr_grpo)

Our config, from `configs.py:184-194`: per_device 16, GA 16, so `steps_per_generation` = 16 and the generation
batch is 256. Also G = 8, T = 1.0, top_k = 0, top_p = 1, max_completion 1536, `vllm_max_model_length` 4096 and
`vllm_gpu_memory_utilization` 0.35.

1. `_prepare_inputs` runs `_generate_and_score_completions` once per 16 micro-steps
   (`_step % (steps_per_generation*num_iterations) == 0`). Then it runs `split_pixel_values_by_grid`, then
   `shuffle_sequence_dict` (module global, imported from `trl.trainer.utils`; `fast.py` swaps it), then
   `split_tensor_dict` into 16, then `unsplit_pixel_values_by_grid`. Source: `T`trainer/grpo_trainer.py#L1591-L1621.
2. `_generate` -> `_tokenize_prompts` -> `_generate_single_turn`. If `state.global_step != _last_loaded_step`, it
   runs `vllm_generation.sync_weights()` under `profiling_context(self,"sync_weights")`
   (`T`trainer/grpo_trainer.py#L1834-L1849). `_last_loaded_step = -1` at init (#L1203), so step 1 always syncs.
3. `VLLMGeneration.generate` (`T`generation/vllm_generation.py#L561-L775):
   - if sleep mode is on and `_llm_weights_sleeping` is True, it calls `sync_weights()` again (#L602-L603);
   - `SamplingParams(n=1, repetition_penalty, temperature, top_p, top_k, min_p, max_tokens=max_completion_length,
     logprobs=0, **generation_kwargs)` (#L675-L697). Colocate always uses **n = 1**: TRL sends 256 rows, each prompt
     8 times;
   - `llm.wake_up(tags=["kv_cache"])` (#L718-L719);
   - rows are `{"prompt_token_ids": ids, "multi_modal_data": {"image": img}}`. The ids are the trainer
     processor's expanded ids (#L722-L731);
   - `llm.generate(vllm_prompts, sampling_params=..., use_tqdm=False)` under `profiling_context(...,"vLLM.generate")`;
   - `extract_logprobs`: a NaN logprob becomes `None` (#L70-L99);
   - `llm.sleep(level=2)`, then `_llm_weights_sleeping = True` (#L771-L773).
4. Back in the trainer: `logprobs = [[lp[0] for lp in seq] ...]` (top-1 = sampled token). `None` becomes NaN when
   `sampling_per_token_logps` is built (#L2569-L2584). `mask_truncated_completions` zeroes rows whose last id is not
   EOS/pad (#L2603-L2609).
5. The old-logprob pass runs when `GA % (spg*num_iterations) != 0` **or** (`use_vllm and
   vllm_importance_sampling_correction`) (#L2757-L2774). This is the C10 hook: FastGRPO turns the flag off around
   `super()` and sets the ratio inside the update forward. `self.use_vllm and self.vllm_importance_sampling_correction`
   appears verbatim at #L2760 and #L2777 (A7 behaviour (a)).
6. Rewards: `_calculate_rewards` (#L1648). Then advantages (section 4).
7. `_compute_loss` (#L3141) per micro-batch (section 4).

Initialisation facts (`T`trainer/grpo_trainer.py#L1150-L1218, `T`generation/vllm_generation.py#L298-L386):
- `LLM(model=model.name_or_path, tensor_parallel_size, gpu_memory_utilization, max_model_len, max_num_seqs=per_device*
  tp*steps_per_generation (=256), enable_sleep_mode, model_impl, distributed_executor_backend="external_launcher",
  seed=process_index//tp, max_num_batched_tokens=4096, logprobs_mode="processed_logprobs", quantization,
  trust_remote_code[, return_sampling_mask=True])`.
- **No `dtype`, `enable_lora`, `enable_prefix_caching`, `mm_processor_kwargs` or `enforce_eager` is passed.** B2.2
  (`vllm_dtype` / `vllm_kwargs`) must wrap the module-level name `trl.generation.vllm_generation.LLM` (imported at
  #L43) with `functools.partial`. That name is a new private coupling for A7.
- `vllm_model_impl` defaults to `"vllm"` in GRPOConfig (`T`trainer/grpo_config.py#L606). The VLLMGeneration
  default is `"auto"`, but the trainer passes the config value.
- `return_sampling_mask` is True only if IS is on, (top_p<1 or top_k>0 or min_p>0), colocate, T>0 and vLLM>=0.28
  (#L1153-L1160). With our T=1 / top_k=0 / top_p=1 it is **False**, so `sampling_mask` is None and the trainer
  normalises over the full vocabulary. That is correct for unfiltered sampling.
- GRPOConfig `__post_init__` raises if IS is on, `top_k == 0` and (`top_p < 1` or `min_p > 0`)
  (`T`trainer/grpo_config.py#L1083-L1093). A future R5 setting "T 0.7 / top_p 0.8 / top_k 20" is legal (top_k>0).
- The trainer computes logprobs with a fused LM head, patched in by `add_fused_lm_head(model.get_base_model(),
  temperature=self.temperature, cast_lm_head_to_fp32=..., outputs=("log_probs","entropy"))`
  (#L1114-L1124; `T`trainer/utils.py#L1569). The forward is `model(**inputs, labels=labels, fused_lm_head=True)`.
  It reads `outputs.log_probs` / `outputs.entropy`, and temperature is applied inside (#L1563-L1567). A
  `FastForward` window override must keep `labels` (-100 outside the completion) and `fused_lm_head=True`.

## 2. Weight sync (A7, C5, B2.2 option F)

TRL path (`T`generation/vllm_generation.py#L446-L528):
- `sync_weights()`: if colocate and sleep mode is on, it runs `empty_cache()`, `llm.wake_up(tags=["weights"])`,
  `_llm_weights_sleeping=False`. Then **one `load_weights` call per parameter**:
  `self.llm.llm_engine.model_executor.driver_worker.model_runner.model.load_weights([(name, param)])` (#L520-L522).
  Then `self.llm.reset_prefix_cache()`.
- `_iter_named_params()` with PEFT: `with self._dist.gather_params(all params): model.merge_adapter()`. It yields
  `(name.removeprefix("base_model.model.").replace(".base_layer",""), param.data)`, skips names that contain
  `model.prefix` ("lora_") or `original_module`, and strips `modules_to_save.default.`. It then calls
  `model.unmerge_adapter()` **in place on the bf16 base**. `fast.lora_merged_params` mirrors this name rule.
- `llm_engine.model_executor` exists only when the engine runs in-process
  (`V`v1/engine/llm_engine.py#L136-L138). `distributed_executor_backend="external_launcher"` forces
  `VLLM_ENABLE_V1_MULTIPROCESSING=0` (`V`config/parallel.py#L1031-L1033). A server-mode or `mp` engine has no such
  attribute.
- Public vLLM alternatives with the same effect: `llm.collective_rpc(method_or_callable, timeout=None, args=(),
  kwargs=None)` and `llm.apply_model(func)` (`V`entrypoints/llm.py#L574-L616). `apply_model(lambda m:
  m.load_weights(it))` is the public way to do C5's "one call". The function is cloudpickled only in multiprocess
  mode; in-process it is called directly [unverified for the external_launcher path, so test on the box].
- `Qwen3VLForConditionalGeneration.load_weights(weights)` uses `AutoWeightsLoader(self).load_weights(weights,
  mapper=hf_to_vllm_mapper)`. The mapper is `model.visual.` -> `visual.`, `lm_head.` -> `language_model.lm_head.`,
  `model.language_model.` -> `language_model.model.` (`V`model_executor/models/qwen3_vl.py#L1969-L1975). So TRL's
  HF names (`model.language_model.layers.N.self_attn.q_proj.weight`) are correct input. Packed modules: `qkv_proj`
  = q, k, v; `gate_up_proj` = gate, up (#L1949-L1959).
- `load_weights` takes any iterable. A generator keeps the peak at one merged tensor. A `list(...)` of all merged
  weights costs about 4.3 GB extra (2.13 B x 2 B).
- vLLM's own `GPUModelRunner.reload_weights` also calls `reset_lora_state()`, `reset_encoder_cache()` and
  `reset_mm_cache()` (`V`v1/worker/gpu_model_runner.py#L5287-L5304). TRL's direct `load_weights` does none of these.
  This is harmless for us, because the vision tower is frozen and LoRA is off in vLLM. Under option F, re-add the
  adapter yourself (section 3).
- With `enable_lora=True`, linear layers are wrapped in `BaseLayerWithLoRA`. That class defines `load_weights`,
  which forwards to `base_layer` (`V`lora/layers/base.py#L17-L40). So a base push by name still works under F, by
  design. Check on the box that the push changes vLLM's base tensors.

Upstream status of the merge-drift bug (checked with `gh api`, 2026-11-09):
- #6688 is **closed** (2026-10-16) by a maintainer: "The bf16 drift is real ... tracked in #7423, with a fix in
  #7427". The maintainer also said: keep the IS defaults (`sequence_mask`); the vLLM+PEFT end-to-end test "is still
  skipped". https://github.com/huggingface/trl/issues/6688
- #7423 "Merged-LoRA export changes frozen base weights" is **open**. https://github.com/huggingface/trl/issues/7423
- PR #7427 is **open, not merged**, with 0 comments (last update 2026-10-03). Its fix keeps a CPU copy of the base and
  writes it back after unmerge, at the cost of one device-to-host copy per sync.
  https://github.com/huggingface/trl/pull/7427
- v1.17.1 (2026-10-30) still has the in-place merge/unmerge (#L459, #L480). So `lora_merged_params` (no write to W)
  stays needed. The bf16 rounding of `W + d` that vLLM receives is the same in both methods. Only the trainer-side
  drift goes away.

## 3. Sleep mode (C5 modes A/B, B2.2-F)

From `V`entrypoints/llm.py#L836-L885 and `V`v1/worker/gpu_worker.py#L294-L353:
- `LLM.sleep(level=1, mode="abort")`:
  - level 0 pauses scheduling only;
  - level 1 offloads the weights to CPU and discards the KV cache;
  - level 2 discards the weights and the KV cache.
  `wake_up(tags=None) -> bool` takes tags from `("weights","kv_cache","scheduling")` and returns "fully awake".
- Level 2 saves `model.named_buffers()` to CPU before it sleeps. It restores them on any wake that includes
  "weights" (#L294-L353). Rotary caches therefore survive, but **parameters do not**: they must be re-pushed (as TRL
  does).
- The cumem backend backs up the `runtime` and `cudagraph` pools to CPU at every level. It keeps `weights` only at
  level 1 (`V`device_allocator/sleep_mode_backend.py#L127-L137). On wake, non-deferrable tags are always restored;
  `weights`/`kv_cache` are restored only when named (`V`device_allocator/cumem.py#L331-L364). On CUDA, discarded
  memory comes back **unzeroed** (only ROCm memsets). So reading weights after `wake_up(["weights"])` and before a
  push returns garbage. Never hash or compare vLLM tensors while the flag is True (C5 check (b)).
- `LLMEngine.sleep(level>=1)` also clears the frontend multimodal processor cache (`renderer.clear_mm_cache()`,
  `V`v1/engine/llm_engine.py#L418-L425). `EngineCore.sleep` pauses the scheduler with `clear_cache=True`, which
  clears the prefix cache (`V`v1/engine/core.py#L959-L978). So under mode B nothing carries over between steps
  (prefix cache, mm cache, KV).
- **LoRA slots are lost at level 2.** `load_lora_model` (the wrapping and the GPU slot tensors) runs inside
  `GPUModelRunner.load_model` (#L5067), which `Worker.load_model` calls under the `"weights"` memory pool
  (`V`v1/worker/gpu_worker.py#L587-L593). So level 2 discards the stacked `lora_a/lora_b` tensors. The worker manager
  still lists the adapter as active, and `activate_adapter` returns early for an active id
  (`V`lora/model_manager.py#L346-L352). Only `add_lora(LoRARequest(..., load_inplace=True))` forces a reload. It runs
  `remove_adapter`, which deactivates, then reloads and re-activates (`V`lora/worker_manager.py#L289-L322). B2.2-F
  step order (wake -> base push -> `add_lora(load_inplace=True)` -> `reset_prefix_cache`) is therefore required
  every sync in mode B.
- [unverified] The Punica mapping metadata is also built inside the LoRA manager. `set_adapter_mapping` skips the
  rebuild when the mapping and the slot layout are unchanged (`V`lora/model_manager.py "set_adapter_mapping"). The
  per-token mapping differs between batches, so a stale buffer is unlikely. The output-level comparator (B2.2 4b) is
  the guard.
- `sleep_mode_offload_cudagraph` (default False) and `sleep_preserve_parameter_names` (globs kept across level 2)
  are new in `ModelConfig` (`V`config/model.py#L366-L387). The second could keep frozen vision weights across
  sleep, which would shrink the mode-B push, but only if the sender omits them [unverified for Qwen3-VL].
- C5 rule "do not switch TRL to level 1": confirmed. Level 1 restores all weights from host RAM on every wake.

## 4. IS correction, loss, advantages (B2.1-B2.4, A7, C8, C10)

Config defaults (`T`trainer/grpo_config.py#L927-L966, #L1040-L1169):

| field | default | notes |
|---|---|---|
| `vllm_importance_sampling_correction` | True | |
| `vllm_importance_sampling_mode` | `"sequence_mask"` | `token_truncate`, `token_mask`, `sequence_truncate`, `sequence_mask` |
| `vllm_importance_sampling_clip_max` | 3.0 | |
| `vllm_importance_sampling_clip_min` | None | must be < clip_max; `*_truncate` needs at least one bound set |
| `vllm_importance_sampling_cap` | deprecated | copied into clip_max with a FutureWarning |
| `loss_type` | `"dapo"` | we set `dr_grpo` |
| `scale_rewards` | `"group"` | True -> "group", False -> "none"; we set "none" |
| `beta` | 0.0 | 0.0 loads no reference model |
| `vllm_mode` | `"colocate"` | |
| `vllm_enable_sleep_mode` | False | we set True |

Ratio (`T`trainer/grpo_trainer.py#L2777-L2825):

```
diff = nan_to_num((old_logps - sampling_logps) * mask, nan=0)      # mask = completion_mask [* tool_mask]
seq modes:  r = exp(sum_t diff)  shape (B,1)       token modes: r = exp(diff)  shape (B,T)
*_truncate: r = clamp(r, clip_min, clip_max)
*_mask:     r = r.masked_fill((r < clip_min or -inf) | (r > clip_max or +inf), 0)
loss:       per_token_loss *= r   (all loss types except vespo, #L3268-L3269)
```

- The sequence ratio is the exponential of a **sum**, not a mean. With cap 3, a row of L = 1271 tokens is dropped
  when the mean trainer-minus-vLLM gap exceeds ln 3 / 1271 = 8.6e-4 nats per token. This is why the merge gap
  matters for B2.
- A NaN sampling logprob gives diff 0, so r = 1. A fully truncated row (mask 0) gives r = 1 and contributes 0 to the
  loss.
- Logged keys: `sampling/sampling_logp_difference/{mean,max}` and `sampling/importance_sampling_ratio/{min,mean,max}`
  (#L2979-L3014). `fast.FastGRPO.log` rebuilds the same keys.
- The off-policy mask (OPSM) uses `off_policy_mask_threshold` (`self.` attribute, #L882; read in `_compute_loss`,
  #L3189-L3203). It uses `sampling_per_token_logps` when present.
- `logprobs_mode="processed_logprobs"`: vLLM returns logprobs after temperature and top-k/top-p
  (`V`config/model.py#L100-L106, #L287-L295). The vLLM default is `raw_logprobs`, so a hand-built `LLM` (comparator,
  B2.2 4b) must pass it too, or its logprobs differ from TRL's at T != 1.

Loss (`T`trainer/grpo_trainer.py#L3205-L3290):

```
coef_1 = exp(logp - old_logp)  (token level; old_logp = logp.detach() if no old pass)
coef_2 = clamp(coef_1, 1-eps_low, 1+eps_high)              # eps 0.2 / eps_high 0.28
per_token_loss = -min(coef_1*A, coef_2*A) [* r_IS] [+ beta*KL]
dr_grpo: loss = sum(per_token_loss*mask) / (B_micro * max_completion_length) / current_gradient_accumulation_steps
```

- The normaliser is the micro-batch row count times **1536**, not the token count. So the update is invariant to the
  micro-batch split (as `fast._prepare_inputs` assumes).
- KL: `exp(ref-logp) - (ref-logp) - 1`, times `coef_1` if `use_bias_correction_kl` (default True). C8 sets beta = 0,
  so no reference pass runs.
- Advantages (`sum_then_normalize`, #L2880-L2908): `A = r - mean_group(r)`. With "group" or "batch", A is divided by
  `(std + 1e-4)`; with "none" it is not divided. Rows where every reward is None get NaN rewards and then A = 0.

## 5. vLLM 0.32.0 offline APIs (C5, C6, C11, B2.2-F)

- `LLM.generate(prompts, sampling_params=None, *, use_tqdm=True, lora_request=None, priority=None,
  tokenization_kwargs=None, mm_processor_kwargs=None) -> list[RequestOutput]` (`V`entrypoints/llm.py#L418-L428).
  `sampling_params` is one object or a list of the same length as the prompts (C11's per-row caps).
- `SamplingParams` fields (`V`sampling_params.py): `n: int = 1` (#L282), `temperature 1.0`, `top_p 1.0`,
  `top_k 0`, `min_p 0.0`, `seed: int|None` (#L322), `max_tokens: int|None = 16` (#L334), `logprobs: int|None =
  None` (#L339; 0 = the sampled token only, N = up to N+1 entries).
  - Pitfall (comparator "fixed seed"): one shared `seed` on 8 copies of a prompt sent as n = 1 rows gives **8
    identical answers**. Use `seed = base + row` or one request with `n=8`.
  - Bit equality across different batch compositions needs `VLLM_BATCH_INVARIANT=1` plus explicit seeds
    (`V`envs.py#L636; warning at `V`sampling_params.py#L602-L610). Otherwise keep B2.2's floor arm.
- `n > 1` in V1: the frontend fans out to child requests (`ParentRequest`, `V`v1/engine/parallel_sampling.py#L13).
  Input processing runs once per parent, so C6 saves host work only, as the plan expects.
- LoRA engine arguments: `enable_lora: bool = False`; `max_lora_rank` must be one of
  {1,8,16,32,64,128,256,320,512} (default 16); `max_loras = 1`; `max_cpu_loras = None`, which means `max_loras`;
  `lora_dtype = "auto"`, which means the model dtype (`V`config/lora.py#L25-L45, #L120-L139). Also
  `enable_tower_connector_lora = False` (we need it off: LoRA is language-model only); `target_modules` restricts
  wrapping. `VLLM_ALLOW_RUNTIME_LORA_UPDATING` only gates the server endpoints. Offline `add_lora` does not need it.
- `LoRARequest(lora_name, lora_int_id>=1, lora_path, base_model_name=None, load_inplace=False,
  is_3d_lora_weight=False)`, a msgspec struct. **Equality and hash use `lora_name` only**
  (`V`lora/request.py).
- `llm.llm_engine.add_lora(req) / remove_lora(id) / list_loras() / pin_lora(id)`
  (`V`v1/engine/llm_engine.py#L469-L483). There is no `LLM.add_lora`; `llm_engine` is the way in.
- On every engine step, `_set_active_loras` -> `set_active_adapters` -> `add_adapter(req)` runs. A request with
  `load_inplace=True` re-reads its folder at **every step** (`V`lora/worker_manager.py#L289-L322). This confirms the
  plan: pass `load_inplace=False` on `generate`.
- **Two-adapter pitfall (B2.2 4b comparator, "policy" id 1 + "b1" id 2).** With the defaults `max_loras=1` and
  `max_cpu_loras=1`, adding "b1" evicts "policy" from the GPU slot **and** from the CPU cache (LRU,
  `V`lora/model_manager.py#L1378-L1420). The next "policy" request then reloads from its folder. Either set
  `max_loras=2` (it costs slot memory and CUDA-graph variants), or re-run the sync after the comparator and keep the
  tmpfs folder alive.
- The adapter folder needs `adapter_model.safetensors` (or `.bin`) plus `adapter_config.json`
  (`V`lora/lora_model.py#L239-L240). Keys are parsed after the `base_model.model.` prefix is stripped and then
  mapped through the model's `hf_to_vllm_mapper` (`V`lora/utils.py#L183-L201). Modules not in the expected set
  raise `ValueError` (#L273-L277).
- `LoRALayerWeights.optimize()` folds `scaling = alpha/r` into `lora_b` (`V`lora/lora_weights.py#L37-L43;
  packed variant #L284). This confirms the plan's `lora_b == 2*B.to(bf16)` check.
- Prefix caching: `enable_prefix_caching=None`, which means on when the model supports it
  (`V`engine/arg_utils.py#L580, #L3021-L3022). It is on for Qwen3-VL [unverified: `is_prefix_caching_supported` was
  not traced]. Block hashes include the multimodal hash, so equal image + text prompts share blocks.
- Encoder cache: `EncoderCacheManager` keys entries by `mm_hash` and shares them across requests
  (`V`v1/core/encoder_cache_manager.py#L83-L134). The 8 copies of a prompt run the vision tower once while the entry
  is cached. `LLMEngine.reset_encoder_cache()` exists (#L410). Call it only if vision weights change (never for us).
- Multimodal processor cache: `mm_processor_cache_gb = 4` (`V`config/multimodal.py#L357). It is keyed by content
  hash, so the HF processor runs once per distinct image, but each row's image is still hashed. Sleep clears it
  (section 3). Per-row `multi_modal_uuids={"image": [uuid]}` replaces content hashing
  (https://docs.vllm.ai/en/v0.32.0/features/multimodal_inputs.html). Our images are rebuilt each step, so the gain is
  only the hash cost [unverified size, so measure in C0 step 6a].

## 6. Qwen3-VL inputs in vLLM and transformers (C1, C2, C3)

Processor (Qwen/Qwen3-VL-2B-Instruct `preprocessor_config.json`, fetched 2026-11-09):

```
{"size":{"longest_edge":16777216,"shortest_edge":65536},"patch_size":16,"temporal_patch_size":2,"merge_size":2,
 "image_mean":[0.5,0.5,0.5],"image_std":[0.5,0.5,0.5],"image_processor_type":"Qwen2VLImageProcessorFast"}
```

- Resize: `smart_resize(h, w, factor=patch_size*merge_size=32, min_pixels=size.shortest_edge,
  max_pixels=size.longest_edge)` (`HF`models/qwen2_vl/image_processing_qwen2_vl.py#L150-L173, #L63-L90).
- `patchify` (#L175-L207) reshapes to `[B, C, gh/m, m, P, gw/m, m, P]`, permutes to `(0,2,5,3,6,1,4,7)`, inserts a
  temporal axis and expands it to `T=2` (a copy, not a second frame), then reshapes to `(B, gh*gw, C*T*P*P = 3*2*16*16
  = 1536)`. `image_grid_thw = [1, gh, gw]`. Tokens per image = `gh*gw/4`. C1's GPU `patchify` must produce this exact
  row order: merge-window-major, then channel-major inside a row.
- vLLM uses the same HF processor with the model defaults, because TRL passes no `mm_processor_kwargs`. Our
  `processor_kwargs` for Qwen3 is only `padding_side` (`configs.py:135-138`), so both sides resize alike. If image
  kwargs are ever added, pass the same ones to `LLM(mm_processor_kwargs=...)` through the B2.2 `LLM` partial. The
  `one_pad` expansion check (`fast.py:670-690`) catches any mismatch.

transformers 5.20.1 `Qwen3VLModel` (`HF`models/qwen3_vl/modeling_qwen3_vl.py):
- `get_rope_index(self, input_ids, mm_token_type_ids, image_grid_thw=None, video_grid_thw=None,
  attention_mask=None, **kw) -> (position_ids (3,B,L), mrope_position_deltas (B,1))` (#L925-L1031).
  **`mm_token_type_ids` is required.** In the forward, `compute_3d_position_ids` raises `ValueError` if
  `image_grid_thw` is given without `mm_token_type_ids` (#L1092-L1110). C2's `rope_template` must pass the
  processor's `mm_token_type_ids`.
- The position rule: text runs get `arange + cur` on all 3 axes. An image run gets t = cur, h = cur + row,
  w = cur + col on the merged grid (`get_vision_position_ids(cur, grid, 1, spatial_merge_size)`, #L873-L923). After
  the image, `cur += max(gh, gw) // merge` (not the token count). Padded columns stay 0, and
  `delta = max + 1 - len`. Host syncs: `.tolist()` plus `itertools.groupby` per row, and `.item()` per grid. This is
  why C2 caches per prompt layout.
- `get_image_features(pixel_values, image_grid_thw)` returns `BaseModelOutputWithDeepstackFeatures`. Its
  `pooler_output` is split per image (rows `t*h*w/merge^2`), and `deepstack_features` is a list (one entry per
  DeepStack level) of per-image splits (#L1033-L1050). Qwen3-VL-2B has 3 levels
  (`config.vision_config.deepstack_visual_indexes`) [unverified count; read the config].
- `Qwen3VLForConditionalGeneration.forward(..., mm_token_type_ids=None, logits_to_keep=0,
  mm_encoder_outputs: dict[str, BaseModelOutputWithPooling] | None = None, **kw)` passes `mm_encoder_outputs` through
  to the model (#L1309-L1378). `Qwen3VLModel.forward` **raises** if both `pixel_values` and `mm_encoder_outputs` are
  given (#L1161-L1164). It concatenates `pooler_output` with `torch.cat` in row order, scatters it into the
  image-pad positions, and stacks DeepStack per level (#L1183-L1226). C3 must drop `pixel_values` but keep
  `image_grid_thw` and `mm_token_type_ids`, which the rope needs unless C2 supplies `position_ids`.
- vLLM's Qwen3-VL vision output is `cat([merger_out] + deepstack_list, dim=1)`, shape
  `[tokens, hidden*(1+levels)]` (`V`model_executor/models/qwen3_vl.py#L895-L905). To feed precomputed embeddings:
  `LLM(..., enable_mm_embeds=True)` (`V`config/multimodal.py#L326) and
  `multi_modal_data={"image": {"image_embeds": T[tokens, hidden*(1+levels)], "image_grid_thw": T[n,3]}}`. The docs
  say that for Qwen3-VL, image_embeds "should contain both the base image embedding and deepstack features", and
  that a wrong shape may crash the engine (https://docs.vllm.ai/en/v0.32.0/features/multimodal_inputs.html;
  parsing at `V`model_executor/models/qwen3_vl.py#L2400-L2474). Not in any task today. It would be the route to
  "encode once for both vLLM and trainer" [unverified end to end].

## 7. Profiling (A6, C0)

- TRL's `profiling_context(trainer, name)` times with `perf_counter` and logs `"profiling/Time taken: <Class>.<name>"`
  **only to wandb, mlflow or trackio**, never to `_metrics` (`T`extras/profiling.py#L95-L131). The names are
  `sync_weights`, `vLLM.generate`, and the decorated methods (`_prepare_inputs` etc.). It emits no NVTX, so A6's own
  ranges are needed.
- The vLLM engine runs in the trainer process (external_launcher). So the A6 `nsys --capture-range=cudaProfilerApi`
  window started by `torch.cuda.profiler.start()` also covers vLLM kernels and CUDA-graph replays.
- vLLM NVTX switches:
  - `VLLM_NVTX_SCOPES_FOR_PROFILING=1` wraps the model-runner phases (`gpu_model_runner: preprocess / forward /
    postprocess / sample / bookkeep`) in `nvtx.annotate`. It imports the **`nvtx` PyPI package**, not
    `torch.cuda.nvtx` (`V`v1/utils.py#L784-L800; `V`envs.py#L2027-L2030). Add `nvtx` to the box lock if used.
  - `enable_layerwise_nvtx_tracing` (EngineArgs / ObservabilityConfig) "doesn't work with CUDA graphs enabled"
    (`V`config/observability.py#L86-L89). TRL leaves graphs on, so do not use it in timing arms.
  - `ProfilerConfig(profiler="torch"|"cuda"|"proton", torch_profiler_dir=...)` with `llm.start_profile(prefix,
    delay_iterations=, max_iterations=)` / `llm.stop_profile()` (`V`config/profiler.py#L66-L80;
    `V`entrypoints/llm.py#L801-L827). "cuda" calls `cudaProfilerStart` itself (`V`profiler/wrapper.py#L666-L680). It
    would collide with A6's single capture range, so do not combine the two.
- Log lines from vLLM ("Sleep mode freed X GiB", `V`v1/worker/gpu_worker.py#L322-L329) give the per-step
  wake/sleep memory for free.

## 8. Trackio (deferred, FR14)

- transformers `TrainingArguments`: `report_to="none"` by default; `project="huggingface"`, `trackio_space_id`,
  `trackio_bucket_id`, `trackio_static_space_id`, `run_name`, `hub_private_repo`. `TrackioCallback.setup` calls
  `trackio.init(project=, name=run_name, space_id=, resume="allow", private=, bucket_id=)`, then `on_log` calls
  `trackio.log({...rewritten logs, "train/global_step": step})` (`HF`integrations/integration_utils.py#L935-L1096;
  `HF`training_args.py#L1037-L1080). `space_id=None` logs locally only.
- In TRL with "trackio" in `report_to`, `ProfilingContext` also calls `trackio.log(metrics, step=...)` once **per
  profiled block** (16+ calls per step), and the `log_completions` table goes to trackio too
  (`T`trainer/grpo_trainer.py#L3456-L3461). Trackio is not installed in `.venv`. The plan's Mac smoke (Q7 at plan
  line 17093) is still the only check that FastGRPO's `time/*` keys survive [unverified].

## 9. Private names and couplings to add to A7 / E5 (beyond the 29 items)

- `trl.generation.vllm_generation.LLM` (module global, B2.2 partial), `VLLMGeneration._llm_weights_sleeping` (C5
  sleep rule), `VLLMGeneration.enable_sleep_mode`, `VLLMGeneration.mode`, `trl.generation.vllm_generation.empty_cache`
  (C5 copies TRL's wake order).
- vLLM (box `env_check.sh`, not `ALL_ITEMS`): `llm.llm_engine.model_executor.driver_worker.model_runner.model`,
  `llm.llm_engine.add_lora`, `LoRARequest(..., load_inplace=)`, `LLM.wake_up(tags=)`, `LLM.sleep(level=)`,
  `LLM.reset_prefix_cache`, `LLM.generate(lora_request=)`.
- The E5 attribute list is verified at the tag: `_step` (#L995, #L1582), `_buffered_inputs` (#L998),
  `num_iterations` (#L991); `shuffle_sequence_dict`, `split_tensor_dict`, `split_pixel_values_by_grid` and
  `unsplit_pixel_values_by_grid` live in `trl/trainer/utils.py` (#L1052-L1224).

Web sources also consulted: https://docs.vllm.ai/en/v0.32.0/features/multimodal_inputs.html,
https://huggingface.co/docs/trl/main/en/grpo_trainer, https://docs.vllm.ai/en/stable/features/sleep_mode/,
https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct/raw/main/preprocessor_config.json.
