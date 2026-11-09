# D. RL post-training infrastructure survey (state as of 2026-11-08)

Scope: is the current stack (TRL 1.17.1 GRPO + "FastGRPO" patches, colocated vLLM with sleep mode, one RTX 5090,
Qwen3-VL-2B + LoRA, 32 prompts x 8 answers, 23.6 s/step, exact CUDA geometry reward) the best RL infra, and what a
top-tier, future-proof stack looks like for 1 GPU now, many GPUs later, VLM RL now and 3D/CAD/physical-AI later.

Labels: **[V]** = checked in a primary source (URL, GitHub API, or the installed package in `.venv`). **[I]** = my
inference or estimate. **[R]** = read in this repo. I spent no money and changed nothing in the repo. I did not read
session transcripts, the memory folder, or `.claude/`.

---

## Skills used

| Skill | What it contributed |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | Checked the live NVIDIA skills catalog (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`). It has a `launch-nemo-rl` skill and a family of `nemo-mbridge-perf-*` skills (sequence packing, CUDA graphs, memory tuning, MoE VLM training) plus `nemo-automodel-*`. None is installed, and I installed nothing. These are worth installing when a NeMo-RL spike starts (§6, stage 3). The main thread's lookups add CT segmentation, Cosmos and Isaac skills for later 3D and physical-AI work. |
| `engineering:system-design` | Gave the shape of §5 (reward service): requirements, then components, API, scale and failure, then trade-offs, plus "what to revisit as it grows". |
| `engineering:architecture` | Gave the ADR format for the keep-vs-migrate decision (§6.1): options, a dimension table, consequences, action items. |
| `engineering:code-review` | A read-only review of `src/objdel/train/fast.py` for upgrade fragility (§1.3). It found 9 overrides of private TRL methods or monkeypatches, which is the main maintenance risk of staying on TRL. |
| `engineering:debug` | Reproduce, isolate, diagnose applied to the step time. It isolated where the 23.6 s goes from the repo's R4 and R5 logs (§1.2), and that breakdown drives the recommendations: generation is about 45 % of the step and is now the largest item. |
| `engineering:documentation` | Report layout: the answer first, then evidence, then the migration runbook. |
| `huggingface-skills:huggingface-best` | Model choice in context. The main thread's leaderboard lookup (`MAIN_skill_lookups.md`) shows the Qwen3.5 small family (0.8B/2B/4B/9B, natively multimodal, Apache-2.0, Feb 2026). Every framework below that supports VLMs now lists Qwen3.5, so framework choice and the next model choice are linked (§4, §6). |
| `huggingface-skills:trl-training` | Baseline TRL GRPO and AsyncGRPO semantics, cross-checked against the installed TRL 1.17.1 source (§2.1). |
| `system-design-skills:api-design` | The reward-service contract (§5): REST/JSON by default with no exotic protocol (the numbers do not force one), a stable error envelope, capped batch size, versioning. Retries are safe by construction because scoring is a pure function, so idempotency keys are not needed. |
| `huggingface-papers` (via the HF papers API, `huggingface.co/api/papers/<id>`) | Checked the title and date of 26 papers (GRPO … score centering) and read 7 abstracts (§3, §4). |

---

## 0. Answer in one screen

1. **On one GPU, the current stack is close to the best available. No framework gives a big single-GPU win for this
   workload.** Every major framework runs the same single-GPU design: synchronous GRPO, vLLM colocated, sleep or
   offload between phases, CUDA-IPC or in-process weight refit (verl hybrid engine, NeMo-RL `colocated.enabled:
   true` + CUDA IPC, OpenRLHF hybrid engine, Unsloth Standby). FastGRPO already has the advanced trainer-side tricks
   the frameworks advertise: shared-prefix forward (verl ships the same thing as `use_prefix_grouper`; TRL
   AsyncGRPO ships it as tree attention), no-merge LoRA sync, vision-feature cache, rope cache, skipping the
   redundant old-logprob pass, IS correction. Moving to verl or NeMo-RL on one GPU would cost weeks and give back
   about the same 20-25 s step [I].
2. **What it is missing is overlap, and overlap needs a second GPU.** Generation (about 10.5 s) and training (about
   13 s) run one after the other [R]. Fully-async frameworks overlap them, but the main async paths put vLLM and the
   trainer on separate GPUs: TRL AsyncGRPO ("must run on separate GPUs") [V], NeMo-RL async GRPO with vLLM
   ("disable colocated inference") [V], verl `fully_async_policy` ("Rollouter and Trainer use separate computing
   resources") [V], prime-rl ("full RL stack … requires 2 GPUs") [V]. On 2 GPUs the expected step is about
   max(10.5, 13) ≈ 13-14 s instead of 23.6 s [I], so about 1.7x throughput for 2x the GPU cost: per dollar it is
   neutral, and the gain is wall-clock.
3. **TRL's own async trainer cannot train on images yet.** `AsyncGRPOTrainer` docs: "Images are not supported yet …
   Multimodal training also needs the padding-free packing path to build 3D M-RoPE positions … which is not
   implemented" [V]. Its reward functions also run in a child process with `CUDA_VISIBLE_DEVICES=""` [V], which
   blocks the in-process CUDA checker. So for async VLM RL the realistic choice today is **verl** or **NeMo-RL**.
4. **Recommendation:** keep TRL + FastGRPO for single-GPU research now. In parallel, (a) pull the reward out into a
   framework-neutral GPU batch service, (b) make datasets, reward and eval framework-agnostic, and (c) run a
   time-boxed parity spike on **verl** (first choice: the biggest community, Qwen3-VL GRPO examples, LoRA over
   vLLM, fully-async, rollout correction, prefix grouper, checkpoint engine, Apache-2.0, v0.9.1) when the project
   moves to 2 or more GPUs. Keep **NeMo-RL** (v0.8.0, 2026-11-06) as the NVIDIA-native alternative for Megatron
   scale and physical-AI (NeMo Gym, Cosmos/Isaac ecosystem). Details and a staged path are in §6.
5. **Algorithm:** for a dense, exact, continuous reward like IoU, keep a **Dr. GRPO / DAPO-style token loss with no
   per-group std scaling** (the repo already uses `dr_grpo`, `scale_rewards="none"`, clip-higher 0.28 [R]). Add
   **masked IS (MIS) against the vLLM logprobs** always (the repo already does), and when async arrives, use
   **CISPO or GSPO plus score centering or MIS**. Consider **GDPO** if the four reward terms (format, legality, IoU,
   topology) start to fight. §4.

---

## 1. Where the current stack stands

### 1.1 What it is [R]

- `pyproject.toml`: torch 2.14.1, transformers 5.20.1, trl 1.17.1, peft 0.21.2. vLLM in colocate mode with sleep mode
  on, `vllm_gpu_memory_utilization` 0.35.
- `configs.grpo_config`: G = 8, `per_device_train_batch_size` 16, `gradient_accumulation_steps` 16,
  gradient checkpointing, `loss_type="dr_grpo"`, `epsilon` 0.2 / `epsilon_high` 0.28, `beta` 1e-3,
  `mask_truncated_completions`, `scale_rewards="none"`, `multi_objective_aggregation="sum_then_normalize"`.
- `fast.py` (732 lines): `FastForward` (per-chunk trimmed widths, optional shared-prefix packed forward),
  `VisionCache`, `RopeCache`, `ImageMemo`, `is_ratio`, `lora_merged_params` (no-merge sync, fixing a bf16 drift TRL
  has, trl #6688), `skip_old_logps` (C10), reward prefetch on a background thread.
- Reward: `objdel.reward.score.score_batch(jobs)`, one CUDA checker call per batch. Expert tables stay on the
  device. 0.004-0.005 s per step (`time/rewards_s`) [R].

### 1.2 Where the step goes (debug method: isolate from logs) [R]

R4 (2026-11-08, RTX 5090, before shared prefix), mean of steps 3-15: step 29.94 s = generate 10.48 + rest of batch
(processor, H2D, ref pass, rewards) 7.0 + update 12.46. Sync / wake / sleep = 1.28 / 0.12 / 0.25 s. Keeping vLLM
awake saved 0.05 s. R5 (`--shared-prefix`): 31.04 → **23.65 s/step**, `shared/token_ratio` 0.361, step 1
bit-identical in reward and grad norm.

Diagnosis:
- Shared prefix shrinks only the training forwards, so generation (about 10.5 s) is now about 45 % of the step [I
  from R4+R5].
- The answers are short (about 105-130 tokens, mean). The earlier review fitted generation as **6.06 s fixed + 4.1 ms
  × completion length** [R], so most generation time is fixed overhead: the processor pass, request handling for
  256 n=1 requests, weight sync, and the tail of runaway 1,536-token answers in 24 % of steps.
- The handovers are not the bottleneck (1.6 s).

What follows for infra: a framework switch does not remove the fixed generation overhead. The levers are (a) cheaper
rollout requests (n = G per request, a length cap, mm-processor caching), and (b) overlap (async, needs 2 GPUs).

### 1.3 Code-review lens: fragility of FastGRPO (read-only)

`FastGRPO` overrides or patches these TRL internals: `_get_per_token_logps_and_entropies`,
`_generate_and_score_completions`, `_prepare_inputs` (swaps `trl_grpo.shuffle_sequence_dict` globally inside a
try/finally), `_compute_loss`, `_generate`, `_calculate_rewards`, `vllm_generation._iter_named_params`,
`vllm_generation.llm.generate`, `sync_weights`/`wake_up`/`sleep` wrappers, plus the vision tower's `forward` [R].
All are private names. TRL ships often: v1.17.1 came out on 2026-10-30 [V]. Each TRL upgrade is a re-validation job.
The repo already compensates with strong parity tests (log-prob atol 1e-5, gradient atol 1e-4) [R].

Specific risks found:
- `packed_chunk` builds a dense W×W boolean mask and refuses flash attention [R]. Memory is O(W²) per group row,
  and SDPA's masked kernel skips no tiles. That is fine at about 1.4k-token packed rows, but it does not scale to
  long answers (the code falls back via `packing_wins`). TRL AsyncGRPO's tree attention uses Hub FlashAttention
  kernels (FA4 on Blackwell) [V], so the upstream version is the better long-term carrier of this idea, once it
  supports images.
- Name collision: "FastGRPO" is also an ICLR 2026 paper on speculative decoding for GRPO (arXiv 2509.21792) [V].
  Renaming the class (for example `ObjdelGRPO`) would avoid confusion in reports. This is cosmetic.

Verdict: the code is correct and well-tested. Its cost is coupling to TRL's private API, not quality.

---

## 2. Frameworks (each checked against its README or docs; versions and dates from the GitHub API on 2026-11-08)

Legend for "1-GPU": **Yes** = documented colocated single-GPU path; **2+** = docs require separate GPUs for RL.

| Framework | Latest release | License / stars | VLM RL | LoRA | 1-GPU | Async / off-policy | Weight sync | Rollout engine | Scale / maturity |
|---|---|---|---|---|---|---|---|---|---|
| **TRL** `GRPOTrainer` | v1.17.1, 2026-10-30 | Apache-2.0 / 19.7k | Yes (Qwen2/2.5/3-VL, Gemma3, SmolVLM2, LLaVA-NeXT) | Yes (PEFT) | Yes (colocate) | Sync. `steps_per_generation`, `num_iterations` give mild off-policy reuse | In-process `load_weights` (colocate); NCCL to `trl vllm-serve` (server) | vLLM | Single node to modest multi-node (accelerate/FSDP/ZeRO). Mature |
| **TRL** `AsyncGRPOTrainer` (experimental) | in 1.17.1 | same | **Text only** ("Images are not supported yet") | Yes; adapter-only sync with `--enable-lora` | **2+** ("must run on separate GPUs") | Fully async, `max_staleness` 4, MIS 0.5-5, **score centering on by default** | vLLM native NCCL weight transfer | vLLM server | FSDP2 only; "intentionally kept minimal" |
| **verl** (verl-project, ex-volcengine) | v0.9.1, 2026-09-20 | Apache-2.0 / 24.1k | Yes (Qwen2.5-VL, Kimi-VL; Qwen3-VL 8B/30B/235B GRPO examples; verl-Omni, verl-vla) | Yes: FSDP+PEFT (merged or adapter sync to vLLM); Megatron-Bridge LoRA incl. `vlm_lora` | Yes (quickstart `n_gpus_per_node=1`) | `one_step_off_policy`, `fully_async_policy` (staleness, partial rollout, streaming); Rollout Correction (TIS/MIS, bypass/decoupled) | Checkpoint Engine: naive (colocated), NCCL, NIXL, Mooncake, kimi_ckpt_engine | vLLM, SGLang, HF; vexact (bitwise-aligned) | 671B MoE, hundreds of GPUs. Most widely used open RL lib |
| **NVIDIA NeMo-RL** | v0.8.0, 2026-11-06 | Apache-2.0 / 2.1k | Yes (VLM GRPO; recipes for Gemma4-E4B geo3k, Nemotron Omni CLEVR/video; Qwen3.5 VLM) | Yes (Automodel + Megatron, GRPO/DPO/SFT) | Yes (default `grpo_math_1B.yaml`: 1 GPU, colocated) | Async GRPO (replay buffer, `max_trajectory_age_steps`, in-flight weight updates); v0.8 Single Controller (streaming, samplers), Megatron backend only | Colocated: CUDA IPC (vLLM via ZMQ+IPC). Separate: NCCL, NCCL reshard, sparse delta (ZMQ/S3), NIXL | vLLM, SGLang, Megatron Inference, TRT-LLM env | Trained Nemotron 3 Nano/Super/Ultra. Container on NGC |
| **slime** (THUDM) | v0.4.0, 2026-10-09 | Apache-2.0 / 8.8k | Yes (`examples/geo3k_vlm`, multimodal docs) | via derivatives (Miles) | not its target [I] | Fully async example, Data Buffer, delta weight sync | Delta weight sync, external engines | SGLang (+router); vLLM variant "vime" | Behind GLM-4.5 … GLM-5.3. Megatron-centric |
| **Miles** (RadixArk, slime fork) | v0.1.3, 2026-10-23 | Apache-2.0 / 3.3k | Yes ("LLM and VLM post-training") | Yes, multi-LoRA into SGLang | not its target [I] | Fully async, TITO, R3 router replay, fault tolerance | P2P RDMA weight transfer | SGLang | Trillion-parameter. MXFP8/NVFP4 RL on Blackwell |
| **OpenRLHF** | v0.11.2, 2026-09-14 | Apache-2.0 / 10.1k | Yes (0.10: VLM RLHF, multi-turn VLM) | Yes | Hybrid engine shares GPUs | `--train.async_enable`, TIS, FlashREINFORCE (2026-09) | NCCL / CUDA IPC | vLLM | 70B+, Ray + DeepSpeed. README now recommends "Molt" (unverified, repo not fetchable) |
| **AReaL** | v2.1.1, 2026-11-04 | Apache-2.0 / 5.9k | Yes (Qwen-VL processors, CLEVR count) | Yes (LoRA math example) | `scheduler.type=local` single node | Fully async by design (2.77x vs sync, boba²); 2.0 = microservices; KPop/IcePop masking | Separate weight-update service | SGLang (default), vLLM | Large-scale agentic RL |
| **ROLL** (Alibaba) | v0.4.0, 2026-09-29 | Apache-2.0 / 3.4k | Yes (Qwen2.5-VL RLVR pipeline; geo3k data) | Yes (incl. Megatron LoRA) | single-node guide | Async rollout and training (ROLL Flash), RollPacker | — | vLLM, SGLang | Ray, Megatron/FSDP2 |
| **SkyRL** (Berkeley/Anyscale) | skyrl-v0.4.0, 2026-10-01 | Apache-2.0 / 2.5k | Yes (Megatron Qwen3.5-VL geometry3k) | Yes (incl. multi-LoRA) | — | Fully async with in-flight weight updates (blog) | vLLM native RL APIs | vLLM | Implements the **Tinker API** on your own GPUs; 397B/501B recipes |
| **prime-rl** (PRIME Intellect) | v0.10.0, 2026-10-13 | Apache-2.0 / 2.2k | Yes (Qwen3.5 / 3.8 VLM image RL) | Yes (LoRA example on 1 H100) | **2+** for RL ("requires 2 GPUs") | Fully async by design | FSDP2 → vLLM | vLLM | 1000+ GPUs, decentralized; tested on RTX 3090/4090/5090 |
| **Unsloth** | v0.1.916-beta, 2026-11-06 | Apache-2.0 / 79k | Yes (Qwen3-VL, Gemma 3 GRPO/GSPO) | Yes (its core) | Yes, single-GPU first; **Standby** shares weights between vLLM and training | Sync | In-process weight sharing | vLLM (`fast_inference`) | Single GPU. vLLM cannot LoRA the vision layers |
| **EasyR1** (hiyouga, verl fork) | v0.3.2, 2025-09-18 (active main) | Apache-2.0 / 5.2k | **Yes, the point of the project** (Qwen2/2.5/3-VL) | Yes (`qwen3_vl_4b_geo3k_grpo_lora.sh`) | Yes (1×12 GB for 1.5B LoRA, its table) | Sync (hybrid engine) | verl hybrid engine | vLLM SPMD | Seg-Zero, VisionReasoner, R1-Track built on it |
| **LlamaFactory** | v0.9.5, 2026-05-30 | Apache-2.0 / 76k | SFT/DPO/PPO-style; VLM RL lives in EasyR1 | Yes | Yes | — | — | — | SFT tool. Not an RL infra choice |
| **Oat** (Sea AI) | v0.2.4, 2025-12-23 | Apache-2.0 / 0.7k | No | Yes (2025-10) | 2+ ("as few as two GPUs") | Actor-learner-oracle | — | vLLM | Research. Low activity since 2026-01 |
| **torchforge** (PyTorch) | none | BSD-3 / 0.7k | — | — | 2+ | — | Monarch | vLLM | **"Development paused … being consolidated in torchtitan"**. Do not adopt. Monarch itself is active (v0.7.0, 2026-10-10) |

Sources: README / docs URLs in §8. Release tags and dates come from `gh api repos/<r>/releases/latest`.

### 2.1 TRL detail (installed 1.17.1, read in `.venv`) [V]

- `loss_type` ∈ {`grpo`, `dr_grpo`, `dapo` (default), `bnpo`, `cispo`, `sapo`, `luspo`, `vespo`}.
  `importance_sampling_level` ∈ {`token`, `sequence`} (sequence = GSPO). `top_entropy_quantile` (80/20 rule),
  `entropy_coef`, `mask_truncated_completions`, `off_policy_mask_threshold`. `trl.experimental` also has `gmpo`,
  `a2po`, `async_grpo`, `async_distillation`, `openreward`, `harbor`.
- vLLM IS correction: `vllm_importance_sampling_mode` (default `sequence_mask`; `token_*`/`sequence_*` ×
  `*_truncate`/`*_mask`), clip max 3.0.
- AsyncGRPO: tree-attention packing ("a token shared by several samples is forwarded once"; "about 2.8-3x at a
  ratio of 3.4 on Qwen3-0.6B/1.7B, H100"), token-budget micro-batching, NCCL weight transfer via vLLM's native
  engine, LoRA adapter sync ("vLLM generates roughly 1.33-1.39x slower per token when it serves a LoRA adapter"),
  default `dtype="float32"` because of the BF16 mismatch finding. Text-only for VLMs.

### 2.2 What "single GPU" really means across frameworks [I, from the docs above]

All colocated single-GPU paths are synchronous: generate → refit → train → refit, with memory handed over by sleep
or offload. On one card, the step is the sum of phases in every framework. Async needs separate devices, except
NeMo-RL's Megatron-Inference colocated-async recipe (`...megatron_generation-colocated-reshard-async-gym.yaml`) [V].
That is Megatron-only and aimed at 30B MoE on 2×8 GPUs, not a 2B VLM on one card. So "the best single-GPU RL
infra" means the leanest sync loop, and that is what FastGRPO has built.

---

## 3. Systems techniques (2025-2026), with where each stands for this project

| Technique | Primary source | Status in the project | Worth it here? |
|---|---|---|---|
| **Fully async rollout/train decoupling** | AReaL (arXiv 2505.24298, 2.77x vs sync) [V]; PipelineRL (2509.19128, "≈2x faster learning", 128 H100) [V]; verl fully_async (2.35-2.67x on 128 GPUs) [V]; NeMo-RL async + Single Controller [V] | No (sync colocate) | **Yes at 2+ GPUs.** About 1.7x wall-clock [I]. Needs an off-policy correction (below) |
| **Partial rollouts** (pause long generations across weight updates, resume) | Kimi k1.5 (2501.12599) [V]; verl `async_training.partial_rollout` [V] | No | Low value: answers are short (about 130 tokens). The 1,536-token tail is better fixed with a length cap |
| **In-flight weight updates** | PipelineRL [V]; vLLM `pause_generation(mode="keep")` + weight transfer, tokens tagged with `weight_versions` [V]; NeMo-RL `in_flight_weight_updates` [V]; SkyRL blog | No | Only with async |
| **Off-policy correction: TIS / MIS / rejection** | Yao et al. "Your Efficient RL Framework Secretly Brings You Off-Policy RL Training" → NeurIPS 2025 "On the Rollout-Training Mismatch in Modern RL Systems" (TIS) [V]; Liu, Li et al. "When Speed Kills Stability" + Trust Region Masking (arXiv 2512.23075) [V via verl doc]; IcePop (Ring-1T, 2510.18855) [V]; AReaL KPop [V] | **Yes**: TRL `sequence_mask` IS from the update's own forward (C10); logged vLLM-vs-trainer gap 0.011 nats [R] | Keep. Required for any async |
| **Score centering** | arXiv 2609.20807 (2026-09-17): "drift" is the cause; additive correction; "matches or outperforms IS under quantization"; composes with IS [V]; default on in TRL AsyncGRPO [V] | No | Cheap to try when going async or FP8. Needs top-k sampler logprobs |
| **FP precision mismatch** | "Defeating the Training-Inference Mismatch via FP16" (2510.26788): BF16 rounding breaks consistency; FP16 fixes it [V]; Oat Precision-RL; Unsloth doc; vexact bitwise-aligned rollout, "KL at exactly zero", "~2x higher final AIME" [V]; vLLM `VLLM_BATCH_INVARIANT=1` (beta) [V] | BF16 + IS; gap 0.011 nats is small [R] | An FP16 A/B is a cheap experiment. Bitwise alignment (vexact) only matters on long runs at scale |
| **CUDA-IPC / NCCL weight sync** | vLLM native weight-transfer API: `init_weight_transfer_engine` / `start_weight_update` / `update_weights` / `finish_weight_update`, backends NCCL and IPC, packed tensors from NeMo-RL (blog 2026-05-28) [V]; NeMo-RL refit transports table [V]; verl Checkpoint Engine [V] | In-process `load_weights` (colocate), 1.28 s/step [R] | Moving to vLLM's IPC engine could cut part of the 1.28 s [I]. Also push only the 196 changed linears (repo note) |
| **Sparse / delta refit** | NeMo-DCR preprint (delta-compressed, "bitwise-exact sparse weight synchronization") [V README]; slime delta weight sync [V] | — | Multi-node only |
| **LoRA adapter-only sync** | TRL AsyncGRPO (`--enable-lora`, about 160 MB at r=32 vs 15 GB merged for 7B) [V]; verl `lora.merge=False` [V] | Merged push (no base write) | Trade-off: the adapter path decodes 1.33-1.39x slower [V]. Merged is right while generation dominates |
| **Sequence packing / padding-free** | verl `use_remove_padding`; NeMo-RL sequence packing; TRL AsyncGRPO token budget [V] | Per-chunk trimming to 64-wide windows [R] | Already mostly covered |
| **KV / prefix sharing across group samples (trainer side)** | Prefix Grouper (arXiv 2506.05433) [V]; verl `use_prefix_grouper` (FSDP only; not with remove_padding or dynamic bsz) [V]; TRL tree attention [V] | **Yes**: `--shared-prefix`, -24 % [R] | Done. Upstream versions validate the approach |
| **Prefix caching on the rollout side** (n samples share one prefill) | vLLM automatic prefix caching; TRL server path sends n=G per request, colocate sends 256 × n=1 [R, repo note] | Colocate n=1 | **Yes, cheap**: group requests as n=8 inside the existing `llm.generate` wrapper [I] |
| **Dynamic batching / token budget** | verl `use_dynamic_bsz`, TRL `token_budget` [V] | Fixed 16×16 | Minor at these lengths |
| **FP8 rollout** | verl FP8 rollout-only and end-to-end (Megatron) [V]; FP8-RL (arXiv 2601.18150): "up to 44% rollout throughput", with token-level TIS [V]; NeMo-RL FP8 [V]; Miles MXFP8/NVFP4 [V] | No | Not now: a 2B model's decode is overhead-bound, not GEMM-bound [I]. verl FP8 rollout needs the Megatron backend [V] |
| **Speculative decoding in rollouts** | SPEC-RL (2509.23232, reuses previous-epoch rollouts as drafts, 2-3x rollout time) [V search]; FastGRPO paper (2509.21792, concurrency-aware + online draft, 2.35-2.72x end-to-end, ICLR 2026) [V]; verl-SpeCo (2026-06) [V]; NeMo-RL speculative decoding (v0.6) [V] | No | Low: at 256-way concurrency speculation gains little (FastGRPO notes plain EAGLE gets only 1.09-1.31x) [V search]. SPEC-RL style reuse fits multi-epoch runs |
| **CUDA graphs for decode** | vLLM default (graph capture in decode); vLLM batch-invariant mode uses breakable CUDA graphs [V] | vLLM default | Already on |
| **Reward servers** | Oat "Oracle … remote service" (Mosec, dynamic batching) [V]; NeMo-RL Environments as Ray actors [V]; slime/Miles/AReaL rollout+reward services; TRL AsyncGRPO rewards in a CPU-only child [V] | In-process CUDA checker, background thread [R] | **Yes**: make it a service (§5). It is the key to framework portability |

---

## 4. Algorithms (2025-2026) and fit for a dense, exact, continuous reward (IoU-like)

Titles and dates checked via the HF papers API [V].

| Method | Paper (arXiv, date) | Core idea | Fit for continuous verifiable reward |
|---|---|---|---|
| GRPO | 2402.03300 (2024-02) | Group mean/std advantage, PPO clip, KL | Baseline. Per-group std scaling over-amplifies tiny IoU spreads |
| **Dr. GRPO** | 2503.20783 (2025-03) | Removes length and std normalization biases | **Good**. In use (`dr_grpo`, `scale_rewards="none"`) [R] |
| **DAPO** | 2503.14476 (2025-03) | Clip-higher, dynamic sampling, token-level loss, overlong shaping | Clip-higher in use (0.28). Dynamic sampling (drop zero-variance groups) rarely triggers with continuous rewards [I] |
| **GSPO** | 2507.18071 (2025-07) | Sequence-level importance ratio and clipping | Good for stability, esp. under async or mismatch; TRL `importance_sampling_level="sequence"` |
| GMPO | 2507.20673 (2025-07) | Geometric mean of token ratios | Stability variant; `trl.experimental.gmpo` |
| **CISPO** | MiniMax-M1, 2506.13585 (2025-06) | Clip the IS weight (stop-grad), keep every token's gradient | **ScaleRL's choice**; TRL `loss_type="cispo"` |
| RLOO | 2402.14740 (2024-02) | Leave-one-out baseline | Unbiased baseline, no std scaling: good with continuous rewards |
| REINFORCE++ (-baseline) | 2501.03262 (2025-01) | Group-mean baseline + global batch normalization | Validated at scale by ScaleRL and ProRL v2 (OpenRLHF README) [V]. Batch-level std is safer than per-group std for dense rewards [I] |
| VAPO | 2504.05118 (2025-04) | Value-based PPO with long-CoT fixes | Needs a critic; overkill for 130-token answers |
| Entropy mechanism: Clip-Cov / KL-Cov | 2505.22617 (2025-05) | Entropy collapse driven by covariance of logprob and advantage; clip or KL-penalize high-cov tokens | Use if entropy collapses (verl has both) |
| 80/20 high-entropy tokens | 2506.01939 (2025-06) | Train only the top-20 % entropy tokens | TRL `top_entropy_quantile`. Polygon JSON has many forced tokens, so it may help [I] |
| ProRL | 2505.24864 (2025-05) | Prolonged RL with KL control and reference resets | Relevant to plateaus over long runs |
| ScaleRL ("The Art of Scaling RL Compute") | 2510.13786 (2025-10) | 400k GPU-h study. Recipe: PipelineRL-8 async, CISPO, FP32 LM head, prompt-level loss averaging, batch-level advantage normalization, zero-variance filtering, no-positive-resampling, forced length interruption [V] | **Best single reference recipe**. "Both GSPO and CISPO substantially outperform DAPO" [V] |
| SAPO | 2511.20347 (2025-11) | Soft temperature gate instead of hard clip | TRL `sapo` |
| GDPO | 2601.05242 (2026-01) | Normalize each reward separately, then combine; GRPO on summed rewards collapses distinct combinations | **Relevant**: the project sums format + legality + IoU + topology before normalizing (`sum_then_normalize`) [R]. NeMo-RL supports GDPO [V] |
| LUSPO | 2602.05261 (2026-02) | Length-unbiased GSPO | TRL `luspo` |
| VESPO | 2602.10693 (2026-02) | Closed-form sequence-level IS reshaping; stable "up to 64x" staleness | For heavily async runs. TRL `vespo` |
| Score centering | 2609.20807 (2026-09) | Additive drift correction | §3 |
| FlashREINFORCE | OpenRLHF (2026-09) | Critic-free, single-rollout async with binary-KL trust region | Very new. Watch |

Fit summary for IoU-like rewards [I, reasoned from the above]:
- Continuous rewards almost never give zero-variance groups, so DAPO's dynamic sampling and ScaleRL's zero-variance
  filtering matter less than they do in math RL. Per-group std normalization is the main hazard: a group with IoUs
  0.61-0.63 gets advantages as large as a group with 0.2-0.9. Keep `scale_rewards="none"` (Dr. GRPO), or use
  batch-level std (REINFORCE++-baseline / ScaleRL).
- The reward is exact and cannot be exploited by noise, so the KL term mostly acts as a regularizer. `beta=0` (Dr.
  GRPO / DAPO default) would also remove the reference pass. That is an owner decision; the repo notes estimate it
  saves about 2-4 s.
- For async, ScaleRL's evidence favors **CISPO** (or GSPO) plus a mismatch correction (MIS or score centering).
- Multi-term reward: try **GDPO** (per-term normalization) as an A/B against the current aggregator.

---

## 5. The reward as a standalone GPU service (system-design + api-design)

### 5.1 Requirements
- Functional: score a batch of completions for named items. Input: item key / view, completion text, truncated flag,
  W, H, object names. Output: the existing terms (`format_ok`, `legal`, `iou_tc`, `iou_min`, `f5`, `g_area`,
  `topo_any`, `topo_freq`, `sig_id`, `error`, `defects`) plus the scalar reward and its spec version. This is
  `score_batch` today [R].
- Non-functional: p99 latency ≤ 50 ms per 256-completion batch on the same host (today 4-5 ms) [R]. Bit-exact with
  the in-process path. Expert masks never leave the device: only per-object numbers return, as today [R]. Usable
  from TRL, verl, NeMo-RL, prime-rl, and TRL AsyncGRPO's CPU-only reward child [V].
- Scale estimate [I]: 256 completions × about 1 KB JSON every 20 s on one GPU, so about 13 KB/s. At 64 GPUs fully
  async with about 20 generations/s per GPU: about 1.3k completions/s, about 1.3 MB/s. Plain HTTP/JSON handles
  this; no gRPC or streaming is needed (api-design: escalate only when a number forces it).

### 5.2 Design
```
 trainer / rollout workers (any framework)                 reward host (same GPU now, own GPU later)
 ┌──────────────────────────────┐   POST /v1/score (batch)   ┌───────────────────────────────────────────┐
 │ thin client: reward_fn(...)  │ ─────────────────────────▶ │ HTTP server (uvicorn/FastAPI or aiohttp)  │
 │  - TRL reward_funcs callable │                            │  └ micro-batcher: coalesce ≤ N ms / ≤ B   │
 │  - verl custom reward fn     │ ◀───────────────────────── │  └ CUDA checker (score_batch)             │
 │  - NeMo-RL Environment actor │   terms + reward + version │  └ ExpertCache tables resident on device  │
 │  - AsyncGRPO CPU child (HTTP)│                            │  └ /v1/spec, /healthz, /metrics           │
 └──────────────────────────────┘                            └───────────────────────────────────────────┘
```
- One process owns the CUDA context and the preloaded expert caches (today's `serve.preload`) [R]. On one GPU it
  shares the card with vLLM and the trainer. The work is a few ms of kernels and small memory, and CUDA time-slices
  contexts. MPS is optional [I].
- Micro-batcher: requests from many rollout workers are coalesced into one `score_batch` call (keeps the
  one-call-per-batch property under async streaming).
- Transport: same host → Unix-domain socket or localhost HTTP. Cross host → HTTP; add msgpack if JSON parsing ever
  shows up in profiles.

### 5.3 Contract
```
POST /v1/score
{ "spec": "reward_v4",                      // must match the server's loaded spec, else 409
  "items": [ { "id": "c_000123",            // caller's id, echoed back
               "key": "lidc/train/0042::shift_+40_-12", "W": 512, "H": 512,
               "names": ["nodule"], "text": "{...polygon json...}", "truncated": false } ] }  // ≤ 4096 items
200 { "spec": "reward_v4", "checker": "<git sha>", "cache_digest": "<sha256 of manifest>",
      "results": [ { "id": "c_000123", "reward": 0.73, "terms": { "format_ok": 1, "legal": 1, "iou_tc": 0.71, ... },
                     "error": "", "defects": {} } ] }
4xx/5xx { "error": { "code": "unknown_key" | "spec_mismatch" | "batch_too_large" | "device_error",
                     "message": "...", "request_id": "...", "retryable": false } }
GET /v1/spec    -> reward spec, weights, group size, checker sha, cache digest
GET /healthz, GET /metrics (batch size, p50/p99 latency, GPU ms, errors by code)
```
- Scoring is a pure function of (spec, cache, text), so retries are safe without idempotency keys.
- Versioning: additive fields only. A new reward spec is a new `spec` name, so one server can host
  `reward_v4` and `reward_v5` for A/B tests. The response carries `checker` and `cache_digest`, so every logged
  reward can be traced to exact code and data.
- Failure: a device error returns 5xx with `retryable:false`, and the run fails loudly (same policy as today's
  `BrokenProcessPool`) [R]. Caller timeout about 5 s.
- 3D/CAD later: the same contract, with `task` added (`"polygon2d" | "mesh3d" | "cad_brep"`) and a payload field per
  task. Each task gets its own checker kernel behind the same micro-batcher. This is the shape NeMo Gym, OpenEnv,
  verl-vla and Isaac-style environments expect: an environment returns rewards, the trainer does not compute them.

### 5.4 What to revisit as it grows
A separate reward GPU when kernels exceed about 5 % of a step. Sharding by dataset when expert caches exceed one
card's memory. A gRPC stream only if per-request overhead shows in profiles.

---

## 6. Recommendation

### 6.1 ADR-D1: RL framework for ObjectDelineation

**Status:** Proposed. **Date:** 2026-11-08. **Deciders:** repo owner.

**Context.** One RTX 5090 now. 2B VLM + LoRA. Short structured answers. A ms-scale exact GPU reward. Plans for
multi-GPU, larger or newer VLMs (Qwen3.5 family), and later 3D/CAD/physical tasks. FastGRPO took the step from 56.6
to 23.6 s [R] using about 730 lines of patches over private TRL APIs.

| Dimension | A: Keep TRL + FastGRPO | B: verl | C: NeMo-RL | D: slime / Miles | E: prime-rl |
|---|---|---|---|---|---|
| 1-GPU speed for this workload | Best known (measured) | About the same [I] | About the same [I] | not targeted | needs 2 GPUs |
| VLM async (2+ GPUs) | No (AsyncGRPO is text-only) | Yes (fully_async + VLM) | Yes (async GRPO; Single Controller VLM GRPO on Megatron) | Yes | Yes |
| LoRA + VLM | Yes | Yes (FSDP PEFT; Megatron `vlm_lora`) | Yes | Miles yes | Yes |
| Scale ceiling | low-mid | very high | very high | very high | very high |
| Off-policy correction | TRL IS modes | Rollout Correction (TIS/MIS/RS presets) | IS correction (required for async) | yes | yes |
| Ecosystem fit for physical AI | low | verl-vla, verl-Omni | NeMo Gym, NVIDIA stack (Cosmos, Isaac skills) | agent connectors | verifiers / Env Hub |
| Migration cost | 0 | medium (Hydra configs, Ray; port reward, data, Qwen3-VL processor quirks) | medium-high (Ray, uv/containers, Megatron for best async) | high | medium |
| Patch maintenance | high (private TRL API) | low-medium (official features) | low-medium | medium | medium |

**Decision.** A now. B as the multi-GPU target, chosen by a parity spike. C as the fallback or co-target when
Megatron scale, NVIDIA physical-AI tooling or NGC containers matter. Not D (Megatron + SGLang at trillion scale is
far beyond need). Not E (no colocated single-GPU path). Not torchforge (paused). Not Oat (inactive). Unsloth is
interesting for memory (Standby weight sharing), but it is another patch layer over TRL and is sync-only.

**Consequences.** Easier: keep shipping research at 23.6 s/step with zero migration risk. Harder: FastGRPO has to be
re-validated on each TRL upgrade (pin TRL; upgrade on purpose). Revisit when: (1) a run regularly gets 2 or more
GPUs; (2) TRL ships image support in AsyncGRPO (then TRL async becomes the cheapest path); (3) a TRL upgrade breaks
more than one FastGRPO hook; (4) 3D/CAD tasks need multi-turn environments.

### 6.2 Staged migration path

**Stage 0 (now, single GPU, no framework change).** Make the project framework-neutral:
1. Reward service per §5, with a thin client that is bit-identical to the in-process `score_batch` (golden test on
   `tests/goldens`). TRL keeps calling it in-process by default. The service mode is a flag.
2. A framework-neutral dataset export: Parquet/JSONL rows of `{prompt messages with image ref, item key, W, H,
   names}`. verl, NeMo-RL, EasyR1 and prime-rl all take a chat-format dataset plus image column [V for verl geo3k,
   EasyR1, NeMo ResponseDataset].
3. A framework-neutral eval: final_eval already scores a checkpoint folder [R], so any framework's merged output
   plugs in.
4. Pin TRL 1.17.1 and record the FastGRPO hook list (§1.3) as an upgrade checklist.

**Stage 1 (single GPU, cheap wins inside the current stack) [I estimates; each needs a measured A/B].**
- Group colocate requests as n = 8 per prompt in the existing `llm.generate` wrapper (prefix cache and processor
  run once per image): 0.2-1 s.
- `max_completion_length` 1,536 → about 1,200 (the LIDC target max is 1,182 [R]): about 0.3-1 s.
- `beta=0` (owner decision): about 2-4 s.
- Use vLLM's native IPC weight-transfer engine, or push only changed linears: part of the 1.28 s sync.
- Optional A/Bs: FP16 vs BF16 (mismatch), CISPO vs dr_grpo, GDPO vs `sum_then_normalize`.

**Stage 2 (2 GPUs, first async).** Spike verl `fully_async_policy` with Qwen3-VL-2B LoRA, 1 rollout GPU + 1 train
GPU, `algorithm.rollout_correction` (MIS preset), staleness ≤ 1, reward via the §5 client as a custom reward
function. Acceptance: (a) same reward curve as TRL at equal samples (±1 SE over 100 steps); (b) wall-clock per 256
samples ≤ 60 % of the TRL single-GPU step; (c) the bit-exact reward golden passes through the service. Budget: a few
GPU-hours on 2×RTX 5090 (US/CA host) [I]. In parallel, read verl's `use_prefix_grouper` limits: it is incompatible
with `use_remove_padding` and `use_dynamic_bsz` [V], so measure which of the two wins for 325-token prompts.

**Stage 3 (8+ GPUs / multi-node).** Pick verl (FSDP2 or Megatron backend, checkpoint-engine NCCL) or NeMo-RL
(Megatron, Single Controller streaming, NCCL-reshard or sparse refit) based on the stage-2 spike and which model is
next. Install the NVIDIA catalog skills `launch-nemo-rl` and `nemo-mbridge-perf-*` at this point (ask first). Use
ScaleRL's recipe as the default config (CISPO, batch-level advantage normalization, FP32 LM head, PipelineRL-style
off-policy k ≤ 8).

**Stage 4 (3D / CAD / physical AI).** Keep rewards as environment services (§5, `task` field). Candidate hosts:
verl-vla (VLA RL across simulators and robots) [V], NeMo-RL + NeMo Gym (environment isolation, multi-env) [V], and
OpenEnv/verifiers connectors (Miles, prime-rl, TRL `openenv_harness`) [V]. Exact geometric checkers (mesh IoU,
B-rep validity, Betti numbers in 3D) slot in like the 2D checker.

---

## 7. Open questions and caveats
- No framework's single-GPU VLM step time was measured here. "About the same" (§0.1, ADR table) is an inference from
  their shared colocated-sync design. It is the claim the stage-2 spike should test, and is the only way to settle
  it.
- verl and NeMo-RL VLM recipes target 8B-235B models. A Qwen3-VL-2B LoRA run with 32×8 groups on 1-2 consumer GPUs
  is not a published recipe there. Expect config work, especially image-token handling: the repo already found
  that TRL's default vLLM image path mismatched Qwen3-VL's pads [R].
- OpenRLHF's recommended successor "Molt" could not be fetched. Unverified.
- Unsloth's speed and VRAM claims are the vendor's own.
- The Qwen3.5 benchmark values in `MAIN_skill_lookups.md` were parsed from model-card HTML. Verify before choosing.

## 8. Sources (all fetched 2026-11-08)
- verl: https://github.com/volcengine/verl ; fully async https://verl.readthedocs.io/en/latest/advance/fully_async.html ;
  LoRA https://verl.readthedocs.io/en/latest/advance/ppo_lora.html ; rollout correction
  https://verl.readthedocs.io/en/latest/algo/rollout_corr.html ; FP8 https://verl.readthedocs.io/en/latest/low_precision/fp8.html ;
  prefix grouper https://github.com/verl-project/verl/tree/main/examples/prefix_grouper ; checkpoint engine
  https://github.com/verl-project/verl/blob/main/verl/checkpoint_engine/README.md ; vexact https://github.com/verl-project/vexact
- NeMo-RL: https://github.com/NVIDIA-NeMo/RL ; v0.8.0 release notes https://github.com/NVIDIA-NeMo/RL/releases/latest ;
  async https://docs.nvidia.com/nemo/rl/latest/guides/async-grpo.html ; GRPO guide https://docs.nvidia.com/nemo/rl/latest/guides/grpo.html ;
  refit https://github.com/NVIDIA-NeMo/RL/blob/main/docs/guides/refit.md ; 1-GPU default
  https://github.com/NVIDIA-NeMo/RL/blob/main/examples/configs/grpo_math_1B.yaml
- slime https://github.com/THUDM/slime ; Miles https://github.com/radixark/miles ; OpenRLHF https://github.com/OpenRLHF/OpenRLHF ;
  AReaL https://github.com/inclusionAI/AReaL ; ROLL https://github.com/alibaba/ROLL ; SkyRL https://github.com/NovaSky-AI/SkyRL ;
  prime-rl https://github.com/PrimeIntellect-ai/prime-rl ; torchforge https://github.com/meta-pytorch/torchforge ;
  Oat https://github.com/sail-sg/oat ; EasyR1 https://github.com/hiyouga/EasyR1
- TRL: https://huggingface.co/docs/trl/main/en/grpo_trainer ; https://huggingface.co/docs/trl/main/en/async_grpo_trainer ;
  installed source `.venv/lib/python3.12/site-packages/trl/{trainer/grpo_config.py,experimental/async_grpo/*}`
- Unsloth: https://unsloth.ai/docs/get-started/reinforcement-learning-rl-guide/vision-reinforcement-learning-vlm-rl.md ;
  https://unsloth.ai/docs/get-started/reinforcement-learning-rl-guide/memory-efficient-rl.md
- vLLM: https://vllm.ai/blog/2026-05-28-native-rl-apis ; https://docs.vllm.ai/en/latest/training/weight_transfer/ ;
  https://docs.vllm.ai/en/stable/training/async_rl/ ; https://docs.vllm.ai/en/stable/features/batch_invariance/
- Papers (HF API checked): 2402.03300, 2503.20783, 2503.14476, 2507.18071, 2507.20673, 2506.13585, 2402.14740,
  2501.03262, 2504.05118, 2505.22617, 2505.24864, 2511.20347, 2602.05261, 2602.10693, 2510.13786, 2506.01939,
  2609.20807, 2510.26788, 2505.24298, 2506.05433, 2509.19128, 2501.12599, 2510.18855, 2601.05242; plus FP8-RL
  2601.18150, SPEC-RL 2509.23232, FastGRPO 2509.21792, Trust Region Masking 2512.23075,
  TIS https://fengyao.notion.site/off-policy-rl , rl-collapse https://richardli.xyz/rl-collapse
- Repo: README.md; src/objdel/train/{fast.py,configs.py,serve.py}; src/objdel/reward/score.py;
  docs/research/notes/{rl_code_optimality_review,r4_speed,shared_prefix,grpo_utilization_research}_2026-11-08.md
