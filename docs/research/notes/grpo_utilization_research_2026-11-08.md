# GPU utilization in the colocated GRPO loop (FastGRPO, Qwen3-VL-2B LoRA, RTX 5090): research 2026-11-08

Scope: how to cut the idle and low-utilization seconds in one GRPO step. Inputs are the owner's measured facts (not
re-measured): 32 images x 8 = 256 answers per step; generate ~7-12 s; update (16 x 16, checkpointing) ~14-19 s;
rewards 0.01 s; step 25-35 s; 1 Hz nvidia-smi mean 67-73 %, 20-30 % of seconds below 50 %, 2-4 s drops to 0 %
several times per step; peak 17 GB of 32. TRL 1.17.1, vLLM 0.32.0, torch 2.13 cu130, transformers 5.19, LoRA r=16
all-linear (vision excluded), beta=1e-3, sleep mode on, `vllm_gpu_memory_utilization` 0.35-0.4.

Labels: **[V]** = verified in the source or page named (URL or installed file + line). **[I]** = my inference or
estimate, not measured. No code in the repo was changed. Prior note with overlapping background:
`docs/research/notes/RL speed stack on RTX 5090/trl_grpo_vllm_colocate.md`.

Line numbers below are from the installed TRL 1.17.1 (`.venv/lib/python3.12/site-packages/trl/...`), which is the
latest release (v1.17.1, 2026-10-30 [V] https://api.github.com/repos/huggingface/trl/releases/latest). TRL `main`
has the same colocate code for sleep/merge/load_weights [V] (fetched
https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py: same calls at main lines
390, 473, 497, 519, 539, 736, 789).

---

## Top 3 recommendations (ranked by seconds saved per unit of work and risk)

| # | Change | Est. saving / step | Work | Main risk |
|---|---|---|---|---|
| 1 | `beta=0.0`: drops the reference forward over 256 answers | ~2-4 s [I] | 1 config line | Owner decision: no KL term (TRL default, DAPO, Dr. GRPO also drop it) |
| 2 | Rewrite the weight sync: out-of-place merge (fixes a confirmed bf16 drift bug), push only the 196 changed linear weights, then A/B sleep off vs. sleep on with a bigger KV | ~1.5-3 s [I] (sleep off) + correctness fix | ~60 lines in FastGRPO (override `VLLMGeneration._iter_named_params` / `_init_vllm`) | Memory: sleep off needs ~11 GB resident; a smaller KV may slow generation |
| 3 | Spend free memory on less recompute: turn checkpointing off on part of the 28 text layers (per-layer flag), or 8 x 32 with no checkpointing | ~1.5-2.5 s [I] | ~10 lines + 2 test arms | Competes with #2 for the same ~15 GB; OOM on long batches |

Before any of these, add timers (no rental cost beyond the next planned run): wrap `sync_weights`, each
`llm.wake_up`/`llm.sleep`, the reference pass and `llm.generate` with `FastGRPO.timed` (it synchronizes, so GPU
time is charged correctly). TRL's own `profiling_context` timings go only to wandb/mlflow/trackio and use
`perf_counter` without a device sync [V] `trl/extras/profiling.py` 56-88, so they are not in the JSON logs now.
Also read two vLLM startup lines: "GPU KV cache size: N tokens" and "Maximum concurrency"; and the preemption count.

Expected combined result if 1 + 2 (sleep off) work: 25-35 s -> ~20-30 s per step [I]. 1 + 3 (sleep on): similar.
2 (sleep off) and 3 cannot both use the same 15 GB; measure both arms.

---

## 1. What TRL 1.17.1 colocate does in one step (sleep mode on)

Order of operations, all on the default stream, in the trainer process [V]:

1. `GRPOTrainer._generate_single_turn`: if `global_step` changed, `vllm_generation.sync_weights()`
   [V] `trainer/grpo_trainer.py` 1836-1841.
2. `VLLMGeneration.sync_weights` [V] `generation/vllm_generation.py` 492-528:
   - `empty_cache()` -> `torch.cuda.empty_cache()` (line 501): the trainer's cached blocks go back to the driver.
   - `llm.wake_up(tags=["weights"])` (502): vLLM maps physical memory for the weights again.
   - `_iter_named_params()` (446-490): `model.merge_adapter()` (459) in place, then yields **every** base
     parameter (text, embeddings, norms, and the frozen vision tower), then `model.unmerge_adapter()` (480).
   - For each tensor, one call `model_runner.model.load_weights([(name, param)])` (520-522): one Python call per
     tensor.
   - `llm.reset_prefix_cache()` (528).
3. `VLLMGeneration.generate` [V] 561-775: `llm.wake_up(tags=["kv_cache"])` (718-719), `llm.generate` on 256
   prompts with `n=1` (677: "vLLM on each GPU generates only 1 in colocate mode"), then `llm.sleep(level=2)`
   (771-773). TRL always uses level 2 (also at init, 376-377).
4. Back in the trainer: no old-logprob pass (FastGRPO C10 skips it), the **reference pass** with the adapter
   disabled because beta != 0 [V] `grpo_trainer.py` 2826-2858, rewards, then 16 update micro-batches. The trainer's
   allocator re-grows its pool after the `empty_cache` [I].

Engine flags are fixed in `_init_vllm`'s `LLM(...)` call: `max_num_batched_tokens=4096`, `enable_sleep_mode`,
`gpu_memory_utilization`, `max_num_seqs` = 16 x 1 x 16 = 256 [V] `vllm_generation.py` 352-372,
`grpo_trainer.py` 1185-1188. There is **no pass-through** for other engine kwargs (only `return_sampling_mask`), so
`enable_lora`, `sleep_preserve_parameter_names` or `weight_transfer_config` need an override of `_init_vllm` [V].

### Which parts are CPU-bound or serial

- `wake_up`/`sleep`: driver memory map/unmap calls, no kernels; level 2 also copies small buffers (rope tables) to
  and from CPU [V] vLLM docs https://docs.vllm.ai/en/latest/features/sleep_mode.html. GPU shows 0 % [I].
- `empty_cache` + later re-allocation: `cudaFree`/`cudaMalloc` synchronize the device [I, standard CUDA behaviour].
- Weight push: ~600+ Python `load_weights` calls (count is my estimate from the Qwen3-VL-2B config: 28 text layers
  x 11 tensors + 24 vision blocks x 12 + mergers/embeddings) [I]. Each is a name mapping in Python plus one small
  copy. The bytes are small for the GPU (~4.2 GB bf16 device-to-device, a few ms of bandwidth) [I]; the cost is
  host-side. The merge itself is 196 small GEMMs (28 layers x 7 projections) [I].
- vLLM input processing for 256 multimodal prompts (hashing, processor) runs on the CPU before the first prefill [I].
- All of this is a strict dependency chain (sync -> generate -> ref -> update). **CUDA streams cannot fill these
  gaps**: at 0 % no kernel exists to overlap with. The fix is to remove the work, not overlap it [I].

### Published wake times (closest to our size)

vLLM blog, A100, vLLM 0.11.0, TP=1 [V] https://vllm.ai/blog/2025-10-26-sleep-mode:

| Model | Level 1 wake | Level 2 wake |
|---|---|---|
| Qwen3-0.6B | 0.26 s | 0.85 s |
| Phi-3-vision-128k-instruct (4.2B) | 0.82 s | 2.58 s |

Sleep (offload/discard) times are not reported. CUDA graphs and the CuMemAllocator are kept across level-2 sleep;
level 2 skips graph capture on wake [V] same page. For Qwen3-VL-2B (2.1B params) I expect ~1-2 s of wake per step,
plus sleep, plus the push: this matches the 2-4 s drops to 0 % [I].

---

## 2. Sleep mode: keep it, drop it, or change level

### Facts [V]

- Level 1: offloads weights to CPU, discards KV. Level 2: discards weights and KV; buffers kept on CPU; meant for
  RLHF weight updates. After level 2 "send every parameter not selected by `sleep_preserve_parameter_names`,
  including unchanged ones." https://docs.vllm.ai/en/latest/features/sleep_mode.html
- `sleep_preserve_parameter_names` exists in **v0.32.0** (our version): glob patterns on `model.named_parameters()`;
  level 2 backs selected parameters up to pageable CPU memory and restores them on wake; "The trainer must omit these
  parameters from updates" [V] https://raw.githubusercontent.com/vllm-project/vllm/v0.32.0/docs/features/sleep_mode.md
  lines 76-100.
- TRL docs: `vllm_enable_sleep_mode` default False; "If you still hit OOM, set `vllm_enable_sleep_mode=True`"; GRPO
  docs suggest tuning `vllm_gpu_memory_utilization` first [V] https://huggingface.co/docs/trl/main/en/grpo_trainer.
- Repo fact: training-side CUDA Graphs (reduce-overhead) also failed on a conflict with vLLM's sleep-mode allocator
  (`cumem_allocator.cpp: invalid argument`) [V] `experiments/phase2/perf/p1_results.md` lines 82-85.

### Memory arithmetic [I]

- Qwen3-VL-2B config [V] https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct/resolve/main/config.json: text hidden
  2048, 28 layers, 16 heads, 8 KV heads, head_dim 128, MLP 6144, vocab 151,936, tied embeddings; vision 24 blocks of
  1024.
- KV per token = 2 x 28 x 8 x 128 x 2 B = 112 KiB. Weights ~4.2 GB bf16.
- vLLM at 0.35 x 31.8 GiB ≈ 11.1 GB. Minus weights (4.2) and activations/graphs (~1-1.5) leaves ~5.4 GB of KV ≈
  47k tokens. 256 answers of a few hundred tokens plus image prompts (~300 tokens; shared only on prefix-cache hits)
  is more than 47k. So generation is **probably already KV-bound** (vLLM runs the 256 in waves or preempts). Check
  the startup line and preemptions before deciding.
- Sleep off: training peak (~17 GB) + vLLM (~11 GB) + context ≈ 28-29 GB of 31.8. It fits, but with little margin.
  P1's no-checkpointing arm reached 30.1 GB before OOM (`p1_results.md` table).

### Options

| Option | What changes | Saves [I] | Cost / risk |
|---|---|---|---|
| A. Sleep off, util 0.30-0.35 | No wake/sleep, no `empty_cache`, no full re-push; push only changed tensors | ~1.5-3 s | ~11 GB held through the update; smaller KV may make generation slower; tighter memory blocks #3 |
| B. Sleep on (level 2), util 0.5-0.55 | Bigger KV: fewer waves/preemptions | Unknown; could cut the 7-12 s if KV-bound | Keeps wake/sleep cost; needs util above the current <= 0.4 rule |
| C. Level 2 + `sleep_preserve_parameter_names` for frozen tensors (`visual.*`, norms, embeddings) | Push only the 196 merged linears after wake | ~0.3-0.8 s | CPU backup copy of preserved tensors each cycle; TRL override needed |
| D. Level 1 | Weights backed up to CPU and restored | None for merged sync (all linears are overwritten anyway) | Extra 4 GB D2H + H2D per step: worse than level 2 for us |

Recommendation: do option A vs. option B as a two-arm test after the sync rewrite (section 3). Sleep off also removes
the allocator conflict that blocked training-side graphs (P1), so it re-opens the F4 manual-graph route [I].
Without sleep mode, `expandable_segments` can be used for the trainer to reduce fragmentation (the repo notes it
breaks sleep mode, `experiments/phase2/perf/implementation.md` line 145) [I that it is safe with vLLM in-process;
check on the box].

---

## 3. Weight sync for LoRA

### 3a. A confirmed correctness bug in the current path [V]

TRL merges the adapter **in place** and unmerges it by subtracting a rounded delta. In bf16 the subtraction does not
cancel the addition, so the frozen base weights drift at every sync.

- PEFT 0.21.2 (installed): merge `base_layer.weight.data += delta_weight`; unmerge
  `weight.data -= delta_weight.to(orig_dtype)` [V] `peft/tuners/lora/layer.py` 1008, 1033-1035.
- TRL maintainer (Oct 16, 2026): "The bf16 drift is real: every weight sync merges and unmerges the adapter in
  place, and the rounding doesn't cancel." [V] https://github.com/huggingface/trl/issues/6688
- Tracking issue and fix PR, both **open**: https://github.com/huggingface/trl/issues/7423,
  https://github.com/huggingface/trl/pull/7427 (the PR keeps a CPU copy of the adapted base weights and writes it
  back after unmerge: "one device-to-host copy of the adapted base weights per sync").
- My CPU simulation (not our model; synthetic 2048x2048 bf16 W, rms 0.02, LoRA delta ~0.8 % of W, 300 merge/unmerge
  cycles): 49 % of base entries changed; rms drift 0.4 % of rms(W), about half the size of the LoRA delta itself
  [measured on CPU, synthetic magnitudes]. Real magnitude depends on our delta size; the effect is the same in kind.
- Effect: trainer and vLLM both see the drifted base, so the IS ratio does not detect it; the saved adapter is
  applied to a different base at evaluation time (issue 6688 comment by rmaguado).

Fix that is cheaper than PR 7427 [I]: for each LoRA linear, yield `W + scale * (B @ A)` as a **temporary** tensor
(largest temp 6144 x 2048 bf16 = 25 MB) and never write `W`. No CPU copy, no drift. Single GPU, so no gather.

### 3b. Push less

- Only the 196 LoRA-targeted linear weights change between steps (28 layers x q, k, v, o, gate, up, down) [I from
  the LoRA config `configs.py` line 95: all-linear, vision excluded]. Norms, embeddings and the 24 vision blocks never
  change. Text linears ≈ 1.41B params ≈ 2.8 GB; the rest (~0.7B) need not be pushed **if vLLM keeps them** (sleep
  off, or level 2 with `sleep_preserve_parameter_names`). With plain level 2 all tensors must be re-sent [V vLLM docs].
- Call `load_weights` once with the full list instead of 196-600 single-item calls [I: fewer Python round trips;
  vLLM's `load_weights` accepts an iterable].

### 3c. vLLM LoRA adapters instead of merged weights

Facts [V]:
- vLLM's `Qwen3VLForConditionalGeneration` implements `SupportsLoRA` with packed `qkv_proj` / `gate_up_proj`
  mappings https://docs.vllm.ai/en/latest/api/vllm/model_executor/models/qwen3_vl/.
- vLLM LoRA loads from a **path** (`LoRARequest(name, int_id, path)`), with `load_inplace=True` to replace an adapter
  under the same name https://raw.githubusercontent.com/vllm-project/vllm/v0.32.0/docs/features/lora.md lines
  218-230. TRL's own client comment: "vLLM has no in-memory adapter API" [V] installed
  `trl/experimental/async_grpo/vllm_client.py` 107.
- TRL GRPOTrainer has no adapter sync in 1.17.1. Open PR https://github.com/huggingface/trl/pull/6007 adds it for
  server **and colocate** (colocate: `add_lora(LoRARequest(..., load_inplace=True))`; `reset_prefix_cache()` after
  every swap because vLLM keys its prefix cache on `lora_name` only, vLLM #42125). The PR author: "vLLM LoRA
  inference is slower than serving merged weights ... not a blanket win". Tested only on Qwen2.5-0.5B text.
- TRL `AsyncGRPOTrainer` (experimental, 1.17.1) already has adapter-only sync (save adapter to a dir, then
  `/v1/load_lora_adapter`) and rejects DoRA, `modules_to_save`, bias, lm_head/embedding targets [V] installed
  `trl/experimental/async_grpo/async_grpo_trainer.py` 709-800.
- HF blog (text model, multi-node): adapter-only sync, IS ratio stayed 0.9993-1.0004 over 126 syncs; warns that too
  few `max_loras` slots silently evict an adapter still in use https://huggingface.co/blog/asyncgrpo-lora-hfjobs.
- Open vLLM issue: LoRA shrink/expand kernels nondeterministic in serving
  https://github.com/vllm-project/vllm/issues/50602.

Assessment [I]: adapter is ~17.4M params (r=16 x 7 projections x 28 layers = 622,592 per layer) ≈ 35 MB bf16, so
the push becomes ~0.1-0.3 s (save to `/dev/shm` + load). But: (1) it only pays if base weights stay resident (sleep
off or level 1), (2) every decode step then runs extra LoRA kernels on all 196 projections; for a 2B model at batch
256 I expect a 10-40 % slower decode (no published number found), which could cost more than the 7-12 s
generation saves, (3) vLLM's LoRA path for Qwen3-VL is not tested in TRL. **Not in the top 3.** Revisit if PR 6007
merges, and measure decode tokens/s with and without `enable_lora` first.

### 3d. vLLM IPC weight transfer

vLLM 0.32 has a weight-transfer engine with an `ipc` backend for trainer and vLLM on the **same GPU** (CUDA IPC
handles, no copy through the host; `packed=True` bounds memory to 1 GiB chunks) [V]
https://docs.vllm.ai/en/stable/training/weight_transfer/ipc/ (doc present at the v0.32.0 tag). TRL main docs say
"Use `ipc` instead when the trainer and the server share a GPU" [V]
https://huggingface.co/docs/trl/main/en/vllm_integration, but the 1.17.1 `VLLMClient` only sets up NCCL [V] grep of
installed `trl/generation/vllm_client.py`. Useful only for a two-process (server on the same GPU) design (section 5).
In colocate the trainer already calls `load_weights` directly in-process, so IPC brings nothing there [I].

---

## 4. CUDA streams and CUDA Graphs

- Rewards: already overlapped (prefetch thread during the logprob passes; reward wait 3 ms in P1) [V]
  `p1_results.md`, `fast.py` `prefetch_rewards`. Nothing left to gain.
- H2D copies: 32 images per step, pinned and non-blocking already (configs.py F3) [V]. Negligible.
- The 0 % drops are host/driver work with no kernels (section 1). Streams do not help there [I].
- Inside the update the GPU was 70 % busy with 798,774 launches and 31.5 s of CUDA runtime API time per step in
  P1's trace: launch- and sync-bound [V] `p1_results.md` lines 87-91. So the **mean** 67-73 % is set mostly by the
  update, not only by the transitions [I]. Streams do not cut launch cost; graphs do.
- vLLM decode already uses CUDA graphs and keeps them across sleep [V] vLLM blog.
- Training-side graphs: reduce-overhead crashed (layer graph reuse + cumem allocator conflict) [V] P1. With sleep
  off, the allocator conflict is gone [I]; F4 manual capture (`torch.cuda.make_graphed_callables` per bucket, plan
  in `implementation.md` line 467) is the remaining route. High effort; do it after #1-#3.
- Cheap launch-count experiment [I]: `model_init_kwargs={"use_kernels": True}` (TRL's replacement for Liger layer
  kernels [V] `grpo_trainer.py` 1098-1102). Gain on Qwen3-VL / sm_120 unknown.

---

## 5. Async / off-policy options on one GPU

Facts [V]:
- TRL `AsyncGRPOTrainer` (experimental): vLLM **server**, NCCL weight transfer (or adapter-only via files),
  `max_staleness` default 4, masked IS between policy and sampler [V] installed `async_grpo_config.py`. No image or
  pixel handling in `async_grpo_trainer.py` / `async_rollout_worker.py` (0 matches for "image|pixel") [V grep]: it is
  **text-only**, so it cannot run our VLM task as shipped.
- TRL docs: "the vLLM server and the trainer must run on separate CUDA devices" for NCCL server mode [V]
  https://huggingface.co/docs/trl/main/en/vllm_integration.
- verl one-step-off-policy recipe: +23 % (FSDP2) / +40 % (Megatron) wall-clock vs colocated sync, but with
  **separate** rollout and training GPU pools; weight sync over NCCL "under 300ms" [V]
  https://verl.readthedocs.io/en/latest/advance/one_step_off.html.
- Asynchronous RLHF (Noukhovitch et al., ICLR 2025, arXiv 2410.18252): one-step off-policy matches on-policy
  win-rate; the speedup comes mostly from separate GPUs [V via search summary of the paper's abstract].

Assessment [I]: on one RTX 5090 the only overlap design is two processes (vLLM server + trainer) sharing the GPU.
Without MPS their kernels time-slice rather than run together; with MPS (needs the daemon, often not available inside
a rented container) they can co-run. Upper bound: step ≈ max(generate, update) instead of the sum, i.e. 25-35 s ->
~15-20 s, minus contention. Costs: both resident at once (~11 + 17 GB), weight sync across processes (IPC or adapter
files), 1-step stale samples (TRL's `sequence_mask` IS handles small gaps; P1 IS ratio mean was 0.99-1.02 on-policy),
and a VLM-aware async rollout path that TRL does not have. Too much work and risk for now; **not recommended** until
#1-#3 are done.

`steps_per_generation` = 2 x accumulation (generate 512 answers per 2 optimizer steps) halves sync/wake per step but
makes half the updates off-policy and brings back the old-logprob pass that C10 removed [V `fast.py` 362-369]. Net
likely negative [I].

---

## 6. Shorter training phase: packing, padding-free, recompute

- TRL 1.17.1 `GRPOConfig` has **no** `padding_free` option [V grep of installed `trainer/grpo_config.py`]. TRL
  already computes logprobs with a chunked fused LM head (no full logits) [V `grpo_trainer.py` 1113-1125], so Liger's
  fused loss gives nothing extra.
- FastGRPO already trims each chunk to its own length in 8 width buckets and sorts the generation batch by length
  (C1, `by_length`) [V `fast.py`]. Remaining padding inside a sorted 16-row chunk is small; standard padding-free
  packing would add little [I].
- A bigger lever is **prompt sharing**: each image prompt (~300 tokens [I]) is forwarded 8 times per pass. TRL's
  async trainer has trie packing that forwards a shared prefix once ("a token shared by several samples is forwarded
  once") with a custom tree attention [V] installed `async_grpo/packing.py` 28-60, 237. Porting it to Qwen3-VL
  (M-RoPE positions, deepstack visual injection) is a research project, not a quick win. Possible gain: ~30-40 % fewer
  tokens per pass if answers are ~400 tokens [I].
- Recompute: full checkpointing costs about one extra forward per update (fwd + recompute + bwd ≈ 4 forward units vs
  3) [I], so ~25 % of the update at most. transformers 5.x decoder layers are `GradientCheckpointingLayer`s with a
  per-layer `gradient_checkpointing` flag [V installed transformers 5.20.1 `modeling_layers.py` 58-83]: after
  `gradient_checkpointing_enable()`, set it to False on N of the 28 text layers. PyTorch selective checkpointing
  (`create_selective_checkpoint_contexts`, `CheckpointPolicy.MUST_SAVE` for `aten.mm`) is the finer tool and needs
  `use_reentrant=False` [V] https://docs.pytorch.org/docs/2.15/checkpoint.html. No-checkpointing hit 30.1 GB in P1
  [V]; half the layers would need roughly +6 GB [I], which fits only with sleep on.
- Alternative arm: micro-batch 8 x 32 with no checkpointing (half the activations per micro-batch) [I]. Also re-test
  32 x 8: P1's loss for 32 x 8 came from padding before `by_length` existed (`implementation.md` line 509) [V]; it
  needs ~29 GB, so sleep on only.

---

## 7. Reference pass (recommendation #1, details)

- beta != 0 runs a no-grad forward over all 256 answers with the adapter disabled every step [V]
  `grpo_trainer.py` 2826-2858. TRL default is `beta=0.0`: "the reference model is not loaded, reducing memory usage
  and improving training speed"; KL "is not essential" (Open-Reasoner-Zero, Dr. GRPO, DAPO) [V]
  https://huggingface.co/docs/trl/main/en/grpo_trainer.
- Size: in P1 "old + ref logprobs" was 12.3 s and C10 removed about half (the old pass) [V `p1_results.md`]. With
  VisionCache and RopeCache since then, the ref pass is probably ~2-4 s now [I]. The owner's numbers leave ~4 s per
  step outside "generate" and "update" (7-12 + 14-19 vs 25-35), which fits ref pass + sync/wake [I].
- Repo reason for 1e-3: with **full** fine-tuning a beta of 0 loads no reference model, so raising it later crashes
  [V `configs.py` lines 3-5]. With LoRA the reference is the adapter-off model, so beta can be raised later without a
  reload [I from TRL code: `use_adapter(model, adapter_name=None)`]. Owner decision.

---

## 8. Smaller items

- vLLM shares encoder outputs of the same image across requests by content hash [V]
  https://docs.vllm.ai/en/latest/api/vllm/v1/core/encoder_cache_manager/. So the 8 copies of an image are not
  encoded 8 times on the vLLM side, if they land in the cache in time [I on the same-step case].
- Passing 32 prompts with `n=8` instead of 256 prompts with `n=1` (as TRL's server path does, line 618-631) cuts
  vLLM's front-end work 8x; TRL's `extract_logprobs` and completion flattening handle `n>1` output order [V
  `vllm_generation.py` 70-100, 746]. Saving unknown, likely < 1 s [I]. FastGRPO's `patch_vllm_prompts` already wraps
  `llm.generate`, so it is the place for it.
- 1 Hz nvidia-smi is too coarse for 2-4 s events; one `nsys` or torch profiler step shows exact gap causes [I].

---

## Sources

- TRL 1.17.1 installed source (`.venv/lib/python3.12/site-packages/trl`), same as tag v1.17.1:
  https://github.com/huggingface/trl/tree/v1.17.1
- TRL main `vllm_generation.py`: https://raw.githubusercontent.com/huggingface/trl/main/trl/generation/vllm_generation.py
- TRL vLLM integration docs: https://huggingface.co/docs/trl/main/en/vllm_integration
- TRL GRPO docs: https://huggingface.co/docs/trl/main/en/grpo_trainer
- TRL issues/PRs: https://github.com/huggingface/trl/issues/6688, https://github.com/huggingface/trl/issues/7423,
  https://github.com/huggingface/trl/pull/7427, https://github.com/huggingface/trl/pull/6007
- vLLM sleep mode: https://docs.vllm.ai/en/latest/features/sleep_mode.html and v0.32.0 tag
  https://raw.githubusercontent.com/vllm-project/vllm/v0.32.0/docs/features/sleep_mode.md
- vLLM sleep mode blog: https://vllm.ai/blog/2025-10-26-sleep-mode
- vLLM LoRA (v0.32.0): https://raw.githubusercontent.com/vllm-project/vllm/v0.32.0/docs/features/lora.md
- vLLM Qwen3-VL model API: https://docs.vllm.ai/en/latest/api/vllm/model_executor/models/qwen3_vl/
- vLLM IPC weight transfer: https://docs.vllm.ai/en/stable/training/weight_transfer/ipc/
- vLLM encoder cache: https://docs.vllm.ai/en/latest/api/vllm/v1/core/encoder_cache_manager/
- vLLM LoRA nondeterminism: https://github.com/vllm-project/vllm/issues/50602
- HF blog, async GRPO + LoRA: https://huggingface.co/blog/asyncgrpo-lora-hfjobs
- Unsloth memory-efficient RL (shared weights, "Standby" keeps weights and frees only KV):
  https://unsloth.ai/docs/get-started/reinforcement-learning-rl-guide/memory-efficient-rl
- verl one-step off-policy: https://verl.readthedocs.io/en/latest/advance/one_step_off.html
- PyTorch checkpoint docs: https://docs.pytorch.org/docs/2.15/checkpoint.html
- Qwen3-VL-2B config: https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct/resolve/main/config.json
- PEFT 0.21.2 installed `peft/tuners/lora/layer.py`; transformers 5.20.1 installed `modeling_layers.py`
- Repo: `src/objdel/train/fast.py`, `src/objdel/train/configs.py`, `experiments/phase2/perf/p1_results.md`,
  `experiments/phase2/perf/implementation.md`
