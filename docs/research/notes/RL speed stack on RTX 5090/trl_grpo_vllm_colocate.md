# TRL GRPOTrainer with vLLM colocate on one GPU (LoRA, Qwen3-VL-2B)

Research date 2026-11-06. Sources read: TRL docs (main branch, rendered), TRL source at tags v1.16.0 and v1.17.0 and main, TRL release notes, vLLM docs, HF blog posts. Source text was fetched through a summarizing tool, so code quotes are as relayed by it. Check them against the installed package before relying on exact wording.

## Version fit: TRL 1.16 with vLLM 0.32 (pins of our image)

### Takeaway
TRL 1.16.0 declares vLLM support only up to 0.31.0. vLLM 0.32 support arrived in TRL 1.17.0 (2026-10-29), so the pin set "TRL 1.16 + vLLM 0.32" is outside TRL's tested range. Moving to TRL 1.17.x also brings an importance-sampling change and a new config validation error.

### Cited Findings
- TRL v1.16.0 `pyproject.toml` declares the vLLM extra as `vllm>=0.22.0,<=0.31.0`, `accelerate>=1.10.1`, `transformers>=4.56.2`, `peft>=0.13.0` — [TRL v1.16.0 pyproject](https://raw.githubusercontent.com/huggingface/trl/v1.16.0/pyproject.toml)
- v1.17.0 (2026-10-29) "adds support for vLLM 0.32.0", raises minimum torch to >=2.6.0, accelerate >=1.11.0, peft >=0.17.0 — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- v1.16.0 (2026-10-15) drops vLLM 0.21.0, fixes FSDP2 vLLM colocate weight sync for tied embeddings, and stops VLM policies in GRPO/RLOO from generating the image placeholder token — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- v1.17.1 (2026-10-30) is the latest patch seen (LFM2.5 tool-call fix only) — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- Current docs (main) state TRL supports vLLM `0.23.0` to `0.33.x` — [TRL vLLM integration docs](https://huggingface.co/docs/trl/main/en/vllm_integration). This conflicts with the v1.16.0 pyproject range (<=0.31.0); main is newer than 1.16.
- v1.17.0 breaking change: with vLLM and IS correction on (the default), `GRPOConfig` raises when `top_p < 1.0` or `min_p > 0.0` is set without `top_k > 0`. v1.17.0 also renormalizes log-probs over the top-k/top-p candidate set in colocate mode (vLLM >= 0.28) — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- v1.15.0: a PEFT adapter on `lm_head` now raises an error; drops Python 3.10; `use_liger_kernel` deprecated in GRPO — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- v1.14.0 (2026-09-25): vLLM range >=0.20,<=0.30; GRPO streams log-probs through the chunked path — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- v1.14.2: fixed a bug where `mask_truncated_completions=True` dropped every finished completion on the vLLM path — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)

### Inferences
- The image pin set (TRL 1.16 / vLLM 0.32) needs either TRL >= 1.17.0 or a deliberate override of the pip constraint. I did not find evidence of a failure at 1.16 + 0.32, only that it is untested. Preflight should import-check and run a 2-step smoke test.
- If we move to 1.17.x, set `top_k` (docs error text suggests e.g. 512) whenever `top_p < 1` or `min_p > 0`, or leave top_p=1, min_p=0.
- torch 2.13 / transformers 5.19 are above any minimum I found; I did not find a compatibility matrix for those.

### Gaps
- No source found confirming TRL 1.16 behaviour on vLLM 0.32 (neither success nor failure).
- I did not verify the torch 2.13 / transformers 5.19 combination against TRL.

## Colocate architecture and how TRL builds the engine

### Takeaway
In colocate mode (the default for `use_vllm=True`) TRL builds a vLLM `LLM` inside the trainer process with `distributed_executor_backend="external_launcher"`. `max_num_seqs = per_device_train_batch_size * vllm_tensor_parallel_size * steps_per_generation`; `max_num_batched_tokens` is hard-coded to 4096; `gpu_memory_utilization` is the config value (default 0.3). I did not find `enforce_eager`, `limit_mm_per_prompt`, or a prefix-caching option set in TRL's engine construction, so these are vLLM defaults.

### Cited Findings
- Colocate: "vLLM runs inside the trainer process and shares GPU memory with the training model... may lead to memory contention. This is the default mode." — [TRL vLLM integration](https://huggingface.co/docs/trl/main/en/vllm_integration)
- Engine built in `_init_vllm` with `max_model_len=self.max_model_length`, `gpu_memory_utilization`, `max_num_seqs`, `enable_sleep_mode`, `distributed_executor_backend="external_launcher"`, `seed=accelerator.process_index // tensor_parallel_size`, `max_num_batched_tokens=4096`, `logprobs_mode="processed_logprobs"` — [TRL generation/vllm_generation.py (main)](https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py)
- Code comment on `max_num_batched_tokens=4096`: "Latest vLLM v1 memory profiler is misled by the high default value (i.e., 32768)" — same source
- GRPOTrainer passes `max_num_seqs=args.per_device_train_batch_size * args.vllm_tensor_parallel_size * args.steps_per_generation`, `gpu_memory_utilization=args.vllm_gpu_memory_utilization`, `max_model_length=args.vllm_max_model_length`, `enable_sleep_mode=args.vllm_enable_sleep_mode` — [TRL v1.16.0 grpo_trainer.py](https://raw.githubusercontent.com/huggingface/trl/v1.16.0/trl/trainer/grpo_trainer.py)
- `vllm_max_model_length` default `None` (inferred from model config); help text: set it to at least max prompt length plus `max_completion_length`. `vllm_gpu_memory_utilization` default 0.3. `vllm_tensor_parallel_size` default 1. `vllm_model_impl` default `"vllm"` (alternative `transformers`) — [TRL v1.17.0 grpo_config.py](https://raw.githubusercontent.com/huggingface/trl/v1.17.0/trl/trainer/grpo_config.py)
- The TRL source does not set `enforce_eager`, `limit_mm_per_prompt`, any mm-processor-cache option, or `enable_prefix_caching` in the engine call; `reset_prefix_cache()` is called after each weight sync — [TRL generation/vllm_generation.py (main)](https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py) (absence reported by the fetch tool on the read portion; not a full-file grep)
- `max_completion_length` default is 512; `generation_kwargs` are passed to `SamplingParams` for vLLM and can override other generation params — [TRL v1.17.0 grpo_config.py](https://raw.githubusercontent.com/huggingface/trl/v1.17.0/trl/trainer/grpo_config.py)
- Colocate originated in TRL v0.18.0 via PR #3394 (vLLM external launcher) — [HF blog: No GPU left behind](https://huggingface.co/blog/vllm-colocate)
- Blog notes: TP with external_launcher broke on vLLM >= 0.8.0 (vLLM issue #15895), fixed by setting an explicit random seed; the blog's example config used `vllm_gpu_memory_utilization=0.3`, `max_completion_length=1024`, `num_generations=4`, `per_device_train_batch_size=4` — [HF blog](https://huggingface.co/blog/vllm-colocate)
- With PEFT and TP > 1, a barrier before generation avoids an NCCL race on non-NVLink hardware — [TRL generation/vllm_generation.py](https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py); also a TRL release note says a GRPO+vLLM colocate+PEFT hang on non-NVLink hardware was fixed — [search summary of TRL release notes](https://releases.sh/sources/src_g_3xbu2oRX52tM8mwXJzq)

### Inferences
- For our config (per_device_train_batch_size and steps_per_generation set so one generation batch = 256 completions on one GPU), `max_num_seqs` equals the number of completions per generation batch: `per_device_train_batch_size * steps_per_generation` (TP=1). With 256 completions, vLLM will try to run up to 256 concurrent sequences, so KV cache must hold them (see memory section). If we want fewer concurrent sequences, we must lower `per_device_train_batch_size * steps_per_generation`, which changes the generation batch size; there is no separate `max_num_seqs` knob in GRPOConfig as far as I found.
- `max_num_batched_tokens=4096` means prefill of 330-token prompts is chunked across ~12 prompts per step; this caps activation memory in the profiler but can limit prefill throughput. Overriding needs a TRL patch or a vLLM-side env/engine args hook; I did not find a config field for it.
- Because `enforce_eager` is not set, vLLM's default (CUDA graphs and torch.compile enabled) should apply. I did not find explicit confirmation of the default compile/graph behaviour inside TRL colocate on vLLM 0.32.
- `generation_kwargs` is the only supported route for sampling params; engine args beyond those above are not exposed (could not confirm there is no `vllm_kwargs`-style escape hatch in 1.16).

### Gaps
- Exact `max_model_len` the trainer should use for VLM prompts (330 + 1536 = 1866 tokens, set `vllm_max_model_length` >= ~1,900) is my arithmetic, not a TRL statement.
- Whether any mm processor cache or prefix caching is on by default in vLLM 0.32 for Qwen3-VL: not verified.
- No `enforce_eager` evidence at all in TRL; whether CUDA graph capture memory is counted inside `gpu_memory_utilization` on 0.32 not verified.

## Weight sync with PEFT LoRA

### Takeaway
TRL does not use vLLM LoRA adapters in GRPO colocate. For PEFT models it merges the adapter into the base weights (`merge_adapter()`), loads each tensor into the vLLM model through `model.load_weights([(name, param)])`, then `unmerge_adapter()`, and resets the prefix cache. It runs once per optimizer (global) step, skipped during gradient accumulation.

### Cited Findings
- Colocate `sync_weights()` loops over named params and calls `self.llm.llm_engine.model_executor.driver_worker.model_runner.model.load_weights([(name, param)])`, followed by `self.llm.reset_prefix_cache()` — [TRL generation/vllm_generation.py](https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py)
- For PEFT: inside `gather_params(list(model.parameters()))` it runs `model.merge_adapter()`, reads merged weights, then `model.unmerge_adapter()`. Adapter-only tensors, `modules_to_save` duplicates and `original_module` entries are skipped, and names are rewritten for vLLM — same source
- Comment: "we must gather the full model at once before merging, as merging adapters in a sharded manner is not supported" (relevant to FSDP/ZeRO-3, not single GPU) — same source
- Timing: `if self.state.global_step != self._last_loaded_step: self.vllm_generation.sync_weights(); self._last_loaded_step = self.state.global_step` — i.e. once per global step, avoiding reloads during grad accumulation — [TRL v1.16.0 grpo_trainer.py](https://raw.githubusercontent.com/huggingface/trl/v1.16.0/trl/trainer/grpo_trainer.py)
- Generation under gradient checkpointing: the code wraps HF-generate in `disable_gradient_checkpointing(...)`; for vLLM the engine is separate, so checkpointing does not touch vLLM generation — [TRL v1.16.0 grpo_trainer.py](https://raw.githubusercontent.com/huggingface/trl/v1.16.0/trl/trainer/grpo_trainer.py) (HF path shown; the vLLM implication is my inference)
- AsyncGRPO (separate trainer) has adapter-only vLLM sync with PEFT since v1.14.0; this is not GRPOTrainer — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- Axolotl reports LoRA adapter sync (~40 MB) ~200 ms vs 350 ms to 5 s for NCCL merged-weight sync, from server mode, not TRL colocate — search summary of [Axolotl vLLM serving docs](https://docs.axolotl.ai/docs/vllm_serving.html); the exact figures were relayed by a search tool, not read on the page
- HF async-RL survey lists NCCL broadcast at ~100-500 ms and LoRA adapter transfer as sub-millisecond; gives no TRL number — [HF blog: async RL landscape](https://huggingface.co/blog/async-rl-training-landscape)
- Correctness issue seen outside TRL: ms-swift PR #9946 fixes an IndexError in `ColumnParallelLinearWithLoRA.set_lora()` when syncing LoRA to vLLM in GRPO colocate for Qwen2/Qwen3, reproduced on vLLM 0.26.0 (name-mapping collision) — [ms-swift PR 9946](https://github.com/modelscope/ms-swift/pull/9946). This is the vLLM-LoRA-adapter path, which TRL colocate does not use.
- Community comments (unvalidated) on the HF blog reported a ZeRO-3 + LoRA error and a 4-bit QLoRA error with colocate — [HF blog](https://huggingface.co/blog/vllm-colocate)

### Inferences
- A 2B model's merge + per-tensor load is a GPU-local copy of ~4.4 GB bf16 (Qwen3-VL-2B, my estimate from parameter count, not sourced); expect it to be sub-second to a few seconds. I found no measured number for TRL colocate sync. Needs one timed measurement in our smoke run.
- The merged weights in vLLM are bf16 copies of base + delta, so LoRA rank has little effect on sync cost; the cost is set by full model size.
- Merge/unmerge in bf16 can accumulate rounding error across steps (inference; no source found), a possible source of train-vs-rollout mismatch that IS correction partly absorbs.
- Does vLLM receive the ViT weights too? The loop covers all named params, minus adapter ones; not verified for the vision tower of Qwen3-VL.

### Gaps
- No measured per-step sync time for a 2B model in TRL colocate.
- No TRL issue found specifically on LoRA merge correctness with Qwen3-VL.

## Sleep mode (vllm_enable_sleep_mode), memory sizing

### Takeaway
TRL's sleep mode always uses vLLM level 2: weights and KV cache are discarded (not copied to CPU), and TRL re-pushes weights from the training model on wake. It adds an extra full weight sync each generation but frees all vLLM GPU memory for the optimizer/backward phase. Default is off; default `vllm_gpu_memory_utilization` is 0.3.

### Cited Findings
- vLLM: "Level 1 sleep will offload the model weights and discard the KV cache" (weights backed up to CPU RAM). "Level 2 sleep will discard both the model weights and the KV cache"; buffers such as rope scaling tensors stay in CPU; level 2 is for updating weights as in RLHF — [vLLM sleep mode docs](https://docs.vllm.ai/en/latest/features/sleep_mode/)
- vLLM `wake_up(tags=["weights"])` / `wake_up(tags=["kv_cache"])` restore selectively; docs give no wake-up timings; sleep can release "up to 90%+ of GPU memory"; sleep supported on CUDA and ROCm; requires the default `cumem` allocator backend for some paths — [vLLM sleep mode docs](https://docs.vllm.ai/en/latest/features/sleep_mode/)
- TRL: `llm.sleep(level=2)` after init and after each generation; `wake_up(tags=["kv_cache"])` before generation; `wake_up(tags=["weights"])` at start of `sync_weights()`; weights are re-pushed from the training model because "vLLM's reload_weights would reload the initial checkpoint from disk" — [TRL generation/vllm_generation.py](https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py)
- TRL docs: sleep mode "offload[s] vLLM parameters and cache to CPU RAM during the optimization step"; "waking the vLLM engine adds host-device transfer latency" — [TRL reducing memory usage](https://huggingface.co/docs/trl/main/en/reducing_memory_usage). This wording says offload, but the source code uses level 2 (discard); the docs wording is loose. Source is more reliable.
- GRPOConfig doc: `vllm_enable_sleep_mode` default False: "If you still get out-of-memory errors, set it to True." `vllm_gpu_memory_utilization` default 0.3, colocate only; HF Space `trl-lib/recommend-vllm-memory` gives a recommendation, add +0.05 or +0.1 buffer if it fails — [TRL GRPO trainer docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- With sleep mode, set `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True` only via that variable name (vLLM reads only `PYTORCH_CUDA_ALLOC_CONF` to work around it); setting only `PYTORCH_ALLOC_CONF` is not enough — [TRL reducing memory usage](https://huggingface.co/docs/trl/main/en/reducing_memory_usage). Without sleep mode, expandable_segments is recommended for long online-RL runs to reduce fragmentation — same page.
- Known sleep issues (vLLM 0.8-era, from the HF blog): level-2 buffers (e.g. BatchNorm stats) not restored, fixed in vLLM PR #16889 with TRL keeping a buffer copy; a segfault at shutdown (vLLM issue #16993) was open at blog time — [HF blog](https://huggingface.co/blog/vllm-colocate). Status on 0.32 not checked.
- 72B example in the blog used level 2 with ZeRO-3 + optimizer CPU offload — [HF blog](https://huggingface.co/blog/vllm-colocate)
- Axolotl docs: colocate is "significantly slower" than server mode because generation and training cannot overlap; practical up to ~3B on 24 GB — search summary of [Axolotl vLLM serving](https://docs.axolotl.ai/docs/vllm_serving.html) (not read in full)

### Inferences
- Memory arithmetic (my estimates; Qwen3-VL-2B ~2B params, KV 112 KB/token bf16 given by the user): 256 sequences x (330 + up to 1,536) tokens worst case = ~478k tokens x 112 KB = ~53 GB, far above any cache. Typical (330 + ~300) x 256 = ~161k tokens = ~18 GB. With prefix sharing across n=8 (330 prompt tokens shared) typical need drops. So the KV pool will be the binding limit and vLLM will queue the long tail rather than fail. Budget is by `gpu_memory_utilization`.
- vLLM `gpu_memory_utilization` is a fraction of total GPU memory for the engine (weights + cuda graphs + KV) in the colocate process; without sleep mode that slice is held during training. On 24 GB, 0.3 = 7.2 GB, leaving ~16.8 GB for training (model ~4.4 GB bf16 + LoRA + activations + logits for 1.9k-token seqs). On 32 GB, 0.3 = 9.6 GB.
- With sleep mode level 2, vLLM's GPU memory is released during the update, so utilization can be set higher (e.g. 0.5-0.6) for a bigger KV pool during generation. The trade-off is an extra weights push per generation and wake overhead; no measured numbers found. Whether TRL's memory profiling at engine init (which sees the training model resident) correctly accounts for this is not verified; the profiling can be confused by other-process-resident memory (the HF Space tool exists for that reason).
- Do not set `PYTORCH_ALLOC_CONF` alone if using sleep mode.

### Gaps
- No measured wake-up/sleep latency for a 2B model; vLLM docs give none.
- No source for the best `vllm_gpu_memory_utilization` on 24/32 GB cards for a 2B VLM; only the default 0.3 and the HF Space tool.
- Status of sleep-related bugs on vLLM 0.32 / TRL 1.16-1.17 not found.

## Importance-sampling correction

### Takeaway
On by default in TRL (`vllm_importance_sampling_correction=True`, mode `sequence_mask`, `clip_max=3.0`). It corrects for the numeric gap between vLLM's sampling logprobs and the trainer's logprobs. It forces TRL to compute `old_per_token_logps` with the trainer model, adding a forward pass.

### Cited Findings
- Defaults: correction True; mode `"sequence_mask"`; options `token_truncate`, `token_mask`, `sequence_truncate`, `sequence_mask`; `clip_max=3.0`; `clip_min=None`; `vllm_importance_sampling_cap` deprecated (use `clip_max`), removed in v2.0.0 — [TRL v1.17.0 grpo_config.py](https://raw.githubusercontent.com/huggingface/trl/v1.17.0/trl/trainer/grpo_config.py); [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- Semantics: truncate (TIS) clamps ratios to [clip_min, clip_max]; mask (MIS) sets ratios outside the range to 0; `sequence_*` is sequence-level, `token_*` token-level — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer); code gate `if self.use_vllm and self.vllm_importance_sampling_correction:` with `sequence_level_is = mode in ["sequence_mask","sequence_truncate"]` and `masked_fill(..., value=0.0)` — [TRL v1.16.0 grpo_trainer.py](https://raw.githubusercontent.com/huggingface/trl/v1.16.0/trl/trainer/grpo_trainer.py)
- Validation: truncate modes need at least one clip bound; `clip_min >= clip_max` errors — [TRL v1.17.0 grpo_config.py](https://raw.githubusercontent.com/huggingface/trl/v1.17.0/trl/trainer/grpo_config.py)
- Old logprobs are computed with the training model when `gradient_accumulation_steps % (steps_per_generation * num_iterations) != 0` or when vLLM IS correction is on; otherwise `old_per_token_logps = None` (use current logps detached) — [TRL v1.16.0 grpo_trainer.py](https://raw.githubusercontent.com/huggingface/trl/v1.16.0/trl/trainer/grpo_trainer.py) (condition as relayed; paraphrased from fragments)
- Logged metrics (when correction on with vLLM): `sampling/sampling_logp_difference/{mean,max}`, `sampling/importance_sampling_ratio/{min,mean,max}` — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- vLLM engine is built with `logprobs_mode="processed_logprobs"` "so temperature scaling/logit tweaking affects the TIS log probs" — [TRL generation/vllm_generation.py](https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py)
- v1.17.0: top-k/top-p renormalization of vLLM logprobs in colocate mode (vLLM >= 0.28), because top-p/min-p filtering previously biased the ratio toward 0; requires `top_k > 0` if `top_p<1`/`min_p>0` — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- vLLM `sampling_mask` feature returns the token IDs surviving top-k/top-p/min-p so training can normalize over the same set; truncated rollout sampling vs full-vocab training softmax "violates the principles of importance sampling and destabilizes training" — [vLLM 0.32 sampling mask docs](https://docs.vllm.ai/en/v0.32.0/training/sampling_mask/)
- Motivation paper: Yao, Liu et al., "Your Efficient RL Framework Secretly Brings You Off-Policy RL Training": different implementations for rollout (vLLM) and training (FSDP) make "on-policy" RL implicitly off-policy; they propose truncated importance sampling; vLLM v1 needed a patch to return the probabilities actually sampled from — [paper page](https://liyuanlucasliu.github.io/publication/feng-off-2025); [search summary incl. NeurIPS abstract](https://neurips.cc/virtual/2025/126728). The Notion blog failed to load, so no magnitude numbers or recommended C were verified.
- `loss_type="dr_grpo"`: divides by a constant (recommended = max completion length) instead of sequence length; TRL default loss_type is `"dapo"` — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- Loss application: `if self.use_vllm and self.vllm_importance_sampling_correction and self.loss_type != "vespo":` ratio multiplies the loss — [TRL v1.16.0 grpo_trainer.py](https://raw.githubusercontent.com/huggingface/trl/v1.16.0/trl/trainer/grpo_trainer.py)

### Inferences
- In `sequence_mask` mode a whole completion is dropped (zero weight) if its sequence-level ratio exceeds 3.0. For our 1,536-token tail completions, sequence-level ratios are the product/mean over many tokens, so long completions are the most likely to be masked. Watch `importance_sampling_ratio/max` and the fraction masked; I did not find a metric reporting the masked fraction.
- For a bf16 VLM with LoRA merged into vLLM weights, the mismatch might be larger than for a plain LLM (inference only; no source).
- Cost: extra forward over 256 completions at ~1.9k max tokens for old logps, on top of the ref forward.
- Keep `top_p=1, min_p=0` to avoid the 1.17 validation error and the filtering bias.

### Gaps
- Exact formulas for sequence-level ratio (mean vs sum of log-ratios) in the TRL code were not captured.
- No measured mismatch magnitudes for Qwen3-VL-2B.
- Paper's recommended truncation constant and results not read.

## Multimodal (VLM) GRPO with vLLM

### Takeaway
VLM GRPO with vLLM works in colocate mode for Qwen2/2.5-VL, Gemma3, LLaVA-NeXT, SmolVLM2 per the docs; Qwen3-VL is not on the documented list and I found no Qwen3-VL-specific TRL issue. Known bug: processor-expanded image tokens are sent to vLLM for Idefics3-style models (issue #6294).

### Cited Findings
- Tested VLMs listed in docs: Gemma3-4B-it, LLaVA-NeXT (Mistral-7B), Qwen2-VL-2B-Instruct, Qwen2.5-VL-3B-Instruct, SmolVLM2-2.2B; dataset needs `prompt` (chat-templated) and `image`/`images` (PIL); "most models work with vLLM (server and colocate)"; example uses `--use_vllm --vllm_mode colocate --use_peft` — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- Colocate image passing: each prompt dict gets `row["multi_modal_data"] = {"image": img_list if len(img_list) > 1 else img_list[0]}` in `generate()` — [TRL generation/vllm_generation.py](https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py)
- Server-mode multimodal needs `--enable-scale-out` on vLLM >= 0.30 — [TRL vLLM integration](https://huggingface.co/docs/trl/main/en/vllm_integration)
- Issue #6294: TRL sends processor-expanded prompt token IDs to vLLM for SmolVLM (Idefics3); prompt grew 38 to 855 tokens (729 image tokens) and vLLM re-expands the placeholders. Colocate: corrupted rollouts with many empty/incoherent completions. Server: hard error "Found 9 runs of image tokens in the prompt but 1 images were processed." Versions TRL 1.5.0.dev0, vLLM 0.23.0, transformers 5.13.0; reproduced on main with vLLM 0.24.0. Open at fetch; fix PR #6718 (colocate and server), follow-up #6719 for tool images; merge status not stated — [TRL issue 6294](https://github.com/huggingface/trl/issues/6294)
- Workaround in the issue: pass unexpanded prompt IDs to vLLM and keep expanded IDs for the training forward; `use_vllm=False` gave coherent output — same issue
- v1.16.0 and v1.17.0: VLM policies can no longer generate the image placeholder token under vLLM — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- v1.15.0: assistant-only loss works on vision datasets — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)
- Docs tip for VLMs: LoRA on vision-language projection layers — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- `limit_mm_per_prompt` and mm processor cache do not appear in TRL's engine call — [TRL generation/vllm_generation.py](https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py)

### Inferences
- #6294 is specific to Idefics3-style expansion. Qwen-VL processors also expand image placeholders (our 256 image tokens), so Qwen3-VL could be affected by the same pre-expansion path, but vLLM's Qwen-VL handling may accept already-expanded IDs (Qwen2.5-VL is documented as working). I did not verify Qwen3-VL. A smoke check should compare vLLM completions (non-empty rate, parse rate) with an HF-generate baseline on a few prompts, as in the issue.
- With 8 completions per prompt the same image is processed 8 times by the vLLM mm processor unless an mm cache dedups it; the default mm processor cache state for 0.32 was not verified.
- Default `limit_mm_per_prompt` in vLLM is permissive (not verified); one image per prompt is fine.

### Gaps
- No Qwen3-VL-specific GRPO+vLLM colocate report found (search returned none).
- Merge status of PR #6718 not known; check whether it is in 1.16/1.17.
- Qwen3-VL not in TRL's documented VLM list.

## Generation batch semantics, old/ref logprobs

### Takeaway
One "generation batch" = `per_device_train_batch_size * num_processes * steps_per_generation` completions (the default `steps_per_generation = gradient_accumulation_steps`). Generation happens every `steps_per_generation * num_iterations` micro-steps, and one global optimizer step consumes `gradient_accumulation_steps` micro-batches. For 256 completions per step on one GPU, set `per_device_train_batch_size * gradient_accumulation_steps = 256` with `steps_per_generation = gradient_accumulation_steps` (or set `generation_batch_size=256`), `num_generations=8`.

### Cited Findings
- `generation_batch_size` default = `per_device_train_batch_size * num_processes * steps_per_generation`; mutually exclusive with `steps_per_generation`; `steps_per_generation` default `gradient_accumulation_steps` — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- `num_iterations` default 1 (mu); with mu=1 the clipped objective simplifies to the original — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- `_prepare_inputs`: `generate_every = steps_per_generation * num_iterations`; when `self._step % generate_every == 0` generate and split into `steps_per_generation` buffered chunks via `split_tensor_dict`; each micro-step takes one slice — [TRL v1.16.0 grpo_trainer.py](https://raw.githubusercontent.com/huggingface/trl/v1.16.0/trl/trainer/grpo_trainer.py)
- `beta` default 0.0: "the reference model is not loaded, reducing memory usage and improving training speed"; DeepSeek-R1 uses 0.001 — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- Reference logprobs for PEFT: no separate ref model; under `if self.beta != 0.0:` the code uses `use_adapter(model, adapter_name="ref" if "ref" in model.peft_config else None)` ("disabling adapters yields the base model"); with beta 0 `ref_per_token_logps = None` — [TRL v1.16.0 grpo_trainer.py](https://raw.githubusercontent.com/huggingface/trl/v1.16.0/trl/trainer/grpo_trainer.py)
- `kl` metric logged only if beta nonzero; `step_time`, `num_tokens`, `completions/*_length`, `completions/clipped_ratio`, `frac_reward_zero_std` are logged — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- `mask_truncated_completions` default False — [TRL v1.17.0 grpo_config.py](https://raw.githubusercontent.com/huggingface/trl/v1.17.0/trl/trainer/grpo_config.py)
- GRPO streams log-probs via a chunked path (v1.14.0) — [TRL releases](https://api.github.com/repos/huggingface/trl/releases?per_page=8)

### Inferences
- Per optimizer step with beta 1e-3 and IS correction on, forward passes over 256 sequences: old logps (trainer model, no grad), ref logps (adapter disabled, no grad), and the training forward + backward = three forwards and one backward. Ref cost can't be avoided at beta>0; setting beta=0 would remove it (a decision for the owner, not a recommendation from the sources).
- Larger `steps_per_generation` with `num_iterations=1` is not off-policy by itself in the sense of more updates per rollout: but since the optimizer steps once per `gradient_accumulation_steps` micro-batches and generation covers `steps_per_generation` micro-steps, `steps_per_generation > gradient_accumulation_steps` means rollouts are used for multiple optimizer steps (stale policy). Keep `steps_per_generation == gradient_accumulation_steps` for strictly on-policy 1 rollout per step (inferred from the code shape; not stated in docs).
- Because vLLM `max_num_seqs = per_device_train_batch_size * steps_per_generation` (TP=1) and that product is 256 for our config, the engine admits 256 concurrent sequences regardless of how the batch is split into micro-batches; micro-batch size only affects training memory.

### Gaps
- Exact `old_per_token_logps` condition text was paraphrased from fragments; verify in the installed source.
- No source measuring the ref/old forward cost for a VLM.

## Step-time splits and tuning tips

### Takeaway
I found no published per-phase step-time breakdown for GRPO with vLLM colocate for a <=3B model on one GPU. Available data: HF blog speedup ratios against non-colocated setups, an AAPO no-vLLM baseline, and qualitative statements. Our own profiling (step_time metric plus phase timers) is needed.

### Cited Findings
- HF blog: colocate up to 1.43x faster than plain TRL for 1.5B at the largest batch size; 1.35x for 7B; for 1.5B, increasing TP reduces colocate performance (small models should avoid TP) — [HF blog](https://huggingface.co/blog/vllm-colocate)
- Offline single-H100 generation throughput: 7B ~6,300 output tok/s, 32B ~1,200 tok/s; 512 rollouts (G=8 x 64 prompts) at 2K tokens ~3 min for 7B — [HF blog: async RL landscape](https://huggingface.co/blog/async-rl-training-landscape)
- TRL docs: generation with the model's own `generate` is the bottleneck in online methods, vLLM serves completions far faster — [TRL vLLM integration](https://huggingface.co/docs/trl/main/en/vllm_integration)
- AAPO paper baseline without vLLM: ~567 s/step for 1.5B GRPO with ~3000-token responses on 52.2 GiB — [AAPO, arXiv 2505.14264](https://arxiv.org/pdf/2505.14264) (via search summary)
- Axolotl: colocate is slower than server mode because generation and training cannot overlap — [Axolotl vLLM serving docs](https://docs.axolotl.ai/docs/vllm_serving.html) (search summary)
- Tuning facts: TRL hard-codes `max_num_batched_tokens=4096`; `reset_prefix_cache()` after every sync, so prefix cache only helps within one generation batch (the n=8 samples share prompts) — [TRL generation/vllm_generation.py](https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py)
- vLLM engine seed `accelerator.process_index // tensor_parallel_size`; for a single GPU this is 0; a fixed seed makes rollouts repeatable across runs — same source (derived from the code line)

### Inferences
- Estimate (rough, mine): 256 completions x ~300 avg tokens = ~77k tokens. A 2B model on a 4090/5090 under vLLM decodes at several thousand to ~tens of thousands tok/s aggregate at batch 256, so the average case is likely seconds to tens of seconds, but the 1,536-token stragglers set the generation wall time (decode of one sequence is latency-bound at tens to hundreds tok/s). Cutting `max_completion_length` or stopping the tail gives the most speed. Not measured.
- The training half (3 forwards + backward on 256 x ~630-1,900 token sequences with gradient checkpointing and an image-token prefix) probably dominates for short completions; verify via profiling.
- Prefix caching: vLLM V1 enables automatic prefix caching by default for most models (my knowledge, not verified for Qwen3-VL on 0.32; multimodal prefix caching depends on mm hashing). I did not verify.

### Gaps
- No measured generation/reward/update split for <=3B colocate on 24/32 GB cards.
- No tuning guidance found for `max_num_batched_tokens` in TRL (it is hard-coded).
- No data on RTX 4090/5090 throughput for Qwen3-VL-2B in vLLM.
- No info on how to set the vLLM engine's `enforce_eager`/cudagraph via TRL config.
