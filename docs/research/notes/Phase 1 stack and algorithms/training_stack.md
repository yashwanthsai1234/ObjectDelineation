# Training stack for SFT + GRPO of the smallest VLMs with a CPU-computed verifiable geometry reward (as of 2026-10-18)

Scope: pick the framework(s) and model candidates for (1) SFT so that a small VLM emits polygon JSON (rings + holes), then (2) GRPO-family RL with a multi-term reward computed on CPU by our exact geometry checker, run as about 10 reward arms × seeds on one rented GPU (RTX 4090 24 GB, or one A100/H100 80 GB), with a "Level 2 controller" that changes curriculum and KL during training. Development is on an Apple M3 16 GB (MPS, no CUDA).

Repo snapshot (GitHub API, fetched 2026-10-18/19, via `gh api repos/...`): stars / latest release / last push.

| Repo | Stars | Latest release (date) | Last push |
|---|---|---|---|
| verl (now `verl-project/verl`) | 23.9k | v0.9.1 (2026-09-20) | 2026-10-19 |
| huggingface/trl | 19.5k | v1.16.0 (2026-10-15) | 2026-10-16 |
| modelscope/ms-swift | 15.9k | v4.6.0 (2026-10-13) | 2026-10-18 |
| OpenRLHF | 10.1k | v0.11.2 (2026-09-14) | 2026-10-05 |
| THUDM/slime | 8.7k | v0.4.0 (2026-10-09) | 2026-10-19 |
| VLM-R1 (om-ai-lab) | 6.0k | v0.2.1 (2025-04-15) | 2026-07-07 |
| AReaL | 5.8k | v2.1.0 (2026-08-25) | 2026-10-19 |
| hiyouga/EasyR1 | 5.2k | v0.3.2 (2025-09-18) | 2026-09-19 (last code commit 2026-07-30; later commits are docs only) |
| alibaba/ROLL | 3.4k | v0.4.0 (2026-09-29) | 2026-10-19 |
| NovaSky-AI/SkyRL | 2.4k | skyrl-v0.4.0 (2026-10-01) | 2026-10-17 |
| prime-rl | 2.2k | v0.10.0 (2026-10-13) | 2026-10-18 |
| NVIDIA-NeMo/RL | 2.1k | v0.7.0 (2026-07-29) | 2026-10-19 |
| LLaMA-Factory | 75.6k | v0.9.5 (2026-05-30) | 2026-10-10 |
| unsloth | 78.2k | v0.1.908-beta (2026-10-15) | 2026-10-19 |
| vLLM | 94.1k | v0.32.0 (2026-10-18) | 2026-10-19 |
| SGLang | 37.3k | v0.5.22 (2026-10-19) | 2026-10-19 |

Source for the table: GitHub REST API (`/repos/{owner}/{repo}` and `/releases/latest`), e.g. [verl](https://github.com/verl-project/verl), [TRL](https://github.com/huggingface/trl), [EasyR1](https://github.com/hiyouga/EasyR1), [ms-swift](https://github.com/modelscope/ms-swift), [OpenRLHF](https://github.com/OpenRLHF/OpenRLHF), [AReaL](https://github.com/inclusionAI/AReaL), [ROLL](https://github.com/alibaba/ROLL), [slime](https://github.com/THUDM/slime), [prime-rl](https://github.com/PrimeIntellect-ai/prime-rl), [NeMo-RL](https://github.com/NVIDIA-NeMo/RL), [SkyRL](https://github.com/NovaSky-AI/SkyRL), [LLaMA-Factory](https://github.com/hiyouga/LLaMA-Factory), [Unsloth](https://github.com/unslothai/unsloth), [vLLM](https://github.com/vllm-project/vllm).

---

## RL frameworks: VLM GRPO support, small-model fit, reward plug-in, algorithms, controller hooks

### Takeaway
For this project TRL `GRPOTrainer` (v1.16.0) fits best. It is the only framework that (a) explicitly tests Qwen3-VL-2B and SmolVLM2 for VLM GRPO, (b) runs on a single GPU and on the Mac CPU for smoke tests, (c) logs each reward function separately, and (d) includes per-reward group normalization (GDPO, `multi_objective_aggregation="normalize_then_sum"`) plus Dr.GRPO, DAPO, GSPO, CISPO and SAPO as config switches. It has one blocking pitfall: **GRPO + vLLM corrupts SmolVLM/Idefics3 prompts** (issue #6294 is open and the fix PR #6718 is unmerged as of 2026-10-18). EasyR1 and verl are the proven high-throughput choices for **Qwen-VL only** (Seg-Zero and VisionReasoner were built on EasyR1). Neither supports SmolVLM.

### Comparison table

| Framework (version) | VLM GRPO | SmolVLM | Qwen3-VL(-2B) | Single GPU | Rollout engine | LoRA | Custom reward plug-in | Per-component reward logging | Algorithms (GRPO variants) | Mid-training hooks for controller |
|---|---|---|---|---|---|---|---|---|---|---|
| **TRL GRPOTrainer** (v1.16.0) | Yes; tested Gemma3, LLaVA-NeXT, Qwen2-VL, Qwen2.5-VL, **Qwen3-VL-2B**, **SmolVLM2-2.2B** | Yes with HF generation; **broken with vLLM** (#6294, fix #6718 open) | Yes (tested) | Yes (colocate vLLM or HF generation; continuous batching) | vLLM colocate/server; transformers continuous batching | Yes (`peft_config`) | List of Python callables; sync **or `async def`** (run with `asyncio.gather`); batched (lists in, list of floats out); dataset columns passed as kwargs; `None` allowed per sample | Yes: `rewards/{func_name}/mean`, `/std`; `log_metric`, `log_extra` callbacks inside reward fn | `loss_type`: grpo, dr_grpo, dapo (default), bnpo, cispo, sapo, luspo, vespo; GSPO via `importance_sampling_level="sequence"`; `scale_rewards` group/batch/none; **GDPO** via `multi_objective_aggregation="normalize_then_sum"`; RLOO is a separate trainer | Standard HF `TrainerCallback`; adaptive entropy built in; `sync_ref_model` (TR-DPO). A dynamic reward-weight scheduler PR was closed unmerged (#5657/#5670) |
| **EasyR1** (v0.3.2, Sep 2025) | Yes (core purpose) | Not listed | Yes (Qwen2/2.5/3-VL; `qwen3_vl_4b_geo3k_grpo(_lora).sh`) | Yes: 1.5B full BF16 1×24 GB; 3B LoRA 1×24 GB; 1.5B LoRA 1×12 GB | vLLM SPMD hybrid engine (Docker vllm 0.11.0) | Yes | Python file with `compute_score(reward_input) -> dict`; `REWARD_TYPE = "sequential"` or `"batch"` | Yes: dict keys (e.g. `overall`, `format`, `accuracy`) | GRPO, DAPO, Reinforce++, ReMax, RLOO, GSPO, CISPO (+ SAPO example script) | YAML config (OmegaConf); no documented runtime hook; would need code fork |
| **verl** (v0.9.1) | Yes (Qwen2.5-VL, Qwen3-VL 8B/30B/235B examples, Kimi-VL, GLM-4V, Qwen3.5 patches) | No model patch for Idefics3/SmolVLM | Examples are 8B+, but the qwen3_vl patch is generic | Possible but built for multi-GPU Ray clusters | vLLM ≥0.8.2 and SGLang | Yes | `custom_reward_function.path/name`, `compute_score(data_source, solution_str, ground_truth, extra_info)`; reward managers naive / batch / dapo / prime | Yes: dict return uses `score` and other keys go to `reward_extra_info` | GRPO, DAPO, GSPO, DrGRPO, RLOO, ReMax, REINFORCE++, PRIME, KL_Cov/Clip_Cov, PPO | Python trainer is subclassable; recipes in a separate `verl-recipe` repo |
| **ms-swift** (v4.6.0) | Yes (400+ multimodal models; Qwen2.5-VL/Qwen3-VL ViT LoRA sync) | Not verified | Yes | Yes (colocate) | vLLM colocate or external `swift rollout` server | Yes | ORM plugin (details not verified) | Yes: `reward/{name}/mean`,`/std` | GRPO documented; other variants not verified in the page fetched | Not verified |
| **Unsloth** (beta) | Yes: Qwen2.5-VL-7B, Qwen3-VL, Gemma-3-4B (Llama 3.2 Vision without vLLM) | Not mentioned | Yes | Yes (its focus) | vLLM `fast_inference=True` (requires `finetune_vision_layers=False`) | LoRA-centred | Uses TRL GRPOTrainer reward functions | As TRL | As TRL (GRPO/GSPO highlighted) | As TRL |
| **OpenRLHF** (v0.11.2) | Yes since 0.10 (2026-04), e.g. Qwen3.5; multi-turn VLM | Not listed | Not explicitly listed | Not documented | Ray + vLLM hybrid engine | Yes (LoRA/QLoRA) | Not verified | Not verified | GRPO/REINFORCE++ family (not verified in detail) | Not verified |
| **AReaL** (v2.1.0) | Yes (Geometry3K, CLEVR-Count examples) | No | Yes (FSDP and Megatron); Qwen3.5 too | Not documented | SGLang / vLLM; async RL | Yes | Not verified | Not verified | GRPO, GSPO, PPO, DAPO, LitePPO, Dr.GRPO, REINFORCE++, RLOO, SAPO, IcePop, KPop, M2PO | Async off-policy control |
| **ROLL** (v0.4.0) | Yes (Qwen2.5-VL RLVR, agentic VLM) | No | Not listed | Megatron-oriented | vLLM/SGLang | Yes | Not verified | Not verified | GSPO etc. | Not verified |
| **NeMo-RL** (v0.7.0) | Yes (VLM SFT + GRPO; Qwen3.5 VLM) | No (NeMo AutoModel lists SmolVLM for SFT, not RL) | Qwen3-VL listed for larger sizes | Yes (1-GPU LLM GRPO examples) | vLLM | Yes (LoRA GRPO since 2026-02) | Environment classes | Not verified | GRPO/GSPO/DAPO, **GDPO** listed | Not verified |
| **prime-rl** (v0.10.0) | Yes (Qwen3.5 VLMs) | No | Qwen3.5 rather than Qwen3-VL | Yes (single H100 example; needs NVIDIA GPU) | vLLM | Yes | `verifiers` environments | Not verified | Not verified | Not verified |
| **SkyRL** (v0.4.0) | Yes (Qwen3-VL-8B GRPO) | No | Yes (8B shown) | Not documented | vLLM/SGLang | Yes | Not verified | Not verified | Not verified | Not verified |
| **slime** (v0.4.0) | Not in README (multimodal via forks such as Relax) | No | No | No (Megatron + SGLang, large scale) | SGLang | n/a | n/a | n/a | n/a | n/a |
| **LLaMA-Factory** (v0.9.5) | **No GRPO** (SFT, RM, PPO, DPO, KTO, ORPO); it points to EasyR1 for GRPO | SFT only | SFT yes (`qwen3_vl` template, 2B–235B) | Yes | n/a | Yes | n/a | n/a | n/a | n/a |
| **VLM-R1** (research code) | Yes (REC/OVD, TRL-based open-r1-multimodal) | No | Qwen2.5-VL, InternVL | Multi-GPU scripts | vLLM | Yes | Python reward fns | Partly | GRPO | Fork-level |

### Cited Findings
- TRL GRPO VLM support: tested Gemma3, LLaVA-NeXT, Qwen2-VL-2B, Qwen2.5-VL-3B, **Qwen3-VL-2B-Instruct**, **SmolVLM2-2.2B-Instruct**. The dataset needs `prompt` plus `image`/`images` columns. — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- TRL vLLM: `use_vllm=True`, `vllm_mode="colocate"` (default, in-process, shares GPU) or `"server"`; `vllm_gpu_memory_utilization` default 0.3. Train/inference mismatch is corrected by truncated IS (`vllm_importance_sampling_correction=True`, mode `sequence_mask`, `clip_max=3.0`). — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- TRL reward functions receive `prompts`, `completions`, `completion_ids`, `trainer_state`, `log_extra`, `log_metric` and any dataset column as kwargs, and return a list of floats (or `None`). They can be sync or `async def` (run concurrently via `asyncio.gather`). Several can be given with `reward_weights`, and each is logged as `rewards/{reward_func_name}/mean` and `/std`. — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- TRL `multi_objective_aggregation`: default `"sum_then_normalize"`. `"normalize_then_sum"` normalizes each reward function within the group, sums with weights, then normalizes at batch level. This is the approach suggested by GDPO (arXiv 2601.05242). — [TRL GRPO docs (config section)](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- GDPO (ICML 2026) shows that applying GRPO to a sum of rewards collapses distinct reward combinations into identical advantages. It normalizes each reward separately before aggregation and beats GRPO on tool calling, math and coding. — [arXiv 2601.05242](https://arxiv.org/pdf/2601.05242); [NVIDIA research page](https://research.nvidia.com/labs/twn/publication/icml_2026_gdpo/)
- TRL loss types: `grpo` (not recommended, length bias), `dr_grpo` (constant = `max_completion_length`), `dapo` (default), `bnpo`, `cispo`, `sapo`, `luspo` (requires `importance_sampling_level="sequence"`), `vespo`. `scale_rewards` is `"group"` (default), `"batch"` or `False`. The doc notes Dr.GRPO recommends no std scaling. `mask_truncated_completions`, `top_entropy_quantile`, adaptive entropy (`use_adaptive_entropy`, `entropy_target` default 0.2) and `off_policy_mask_threshold` are also available. — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- TRL KL: `beta` defaults to `0.0`, and then the reference model is **not loaded**. DeepSeek-R1 used 0.001. `sync_ref_model` / `ref_model_mixup_alpha` / `ref_model_sync_steps` implement TR-DPO-style reference updates. — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- TRL GRPOConfig defaults: `num_generations=8`, `max_completion_length=512`, `epsilon=0.2`, `learning_rate=1e-6`, `gradient_checkpointing=True`, `bf16=True`. The single-GPU guidance is to decouple generation from optimization with `generation_batch_size` / `steps_per_generation`. `use_transformers_continuous_batching=True` requires `transformers>=5.8.0` and is described as more efficient than vLLM for N<32. — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- **TRL pitfall (open):** "GRPO + vLLM corrupts SmolVLM multimodal prompts from pre-expanded image tokens". TRL sends processor-expanded `prompt_token_ids`, and vLLM expands the image placeholders a second time, giving empty or incoherent rollouts in colocate mode and `ValueError: Found 9 runs of image tokens…` at step 0 in server mode (Idefics3 interleaves `<image>` runs with `<row_i_col_j>` tokens). The issue was reproduced on main on 2026-08-04, and the fix PRs #6718/#6719 are still open/unmerged. The report says Qwen-style VLMs expand into one contiguous `<|image_pad|>` run and so are unaffected in server mode. — [TRL issue #6294](https://github.com/huggingface/trl/issues/6294); [PR #6718](https://github.com/huggingface/trl/pull/6718)
- **TRL pitfall (open, 2026-10-16):** with Gemma3 and `vllm_mode="server"` on transformers 5.x, vLLM silently generates as if no image were present. This affects GRPO, RLOO and the other `VLLMGeneration` trainers. Another open PR (#7698) stops VLM policies from generating the image placeholder token under vLLM. — [TRL issue #7772](https://github.com/huggingface/trl/issues/7772); [PR #7698](https://github.com/huggingface/trl/pull/7698)
- TRL: PRs "Add dynamic reward weight scheduler support for GRPO and RLOO" (#5657, #5670, April 2026) were closed without merge. — GitHub API on [huggingface/trl PR 5657](https://github.com/huggingface/trl/pull/5657)
- EasyR1 is a "clean fork of veRL to support vision language models". Models: Qwen2-VL/Qwen2.5-VL/Qwen3-VL (plus Llama3/Qwen text). Algorithms: GRPO, DAPO, Reinforce++, ReMax, RLOO, GSPO, CISPO. Also padding-free training, LoRA, and W&B / SwanLab / MLflow / TensorBoard tracking. Requirements: transformers≥4.54, flash-attn≥2.4.3, vllm≥0.8.3, Docker `hiyouga/verl:ngc-th2.8.0-cu12.9-vllm0.11.0`. — [EasyR1 README](https://github.com/hiyouga/EasyR1)
- EasyR1 estimated hardware for GRPO: full FT AMP 1.5B = 2×24 GB, 3B = 4×40 GB; full FT BF16 1.5B = **1×24 GB**, 3B = 1×40 GB; LoRA AMP 1.5B = **1×12 GB**, 3B = **1×24 GB**, 7B = 2×32 GB. — [EasyR1 README](https://github.com/hiyouga/EasyR1)
- EasyR1 reward format: a Python file with `REWARD_TYPE = "sequential"` (or batch) and `compute_score(reward_input: dict) -> dict[str, float]` returning `{"overall": ..., "format": ..., "accuracy": ...}`. The example files are `r1v.py`, `math.py`, `dapo.py`, `android_gui.py`. — [EasyR1 examples/reward_function/r1v.py](https://github.com/hiyouga/EasyR1/blob/main/examples/reward_function/r1v.py)
- EasyR1 example scripts include `qwen3_vl_4b_geo3k_grpo.sh`, `qwen3_vl_4b_geo3k_grpo_lora.sh`, `qwen2_5_vl_3b_*`, plus `_dapo`, `_gspo`, `_cispo`, `_sapo` and `_reinforce` variants. No SmolVLM scripts exist. — [EasyR1 examples dir](https://github.com/hiyouga/EasyR1/tree/main/examples)
- verl: VLM support (Qwen2.5-VL, Kimi-VL). Algorithms PPO, GRPO, GSPO, ReMax, REINFORCE++, RLOO, PRIME, DAPO, DrGRPO, KL_Cov & Clip_Cov. vLLM and SGLang rollouts, LoRA RL, and function-based "verifiable" rewards. 2026 news covers VeRL-Tinker (2026/08), VeRL-Omni v0.2.0 (2026/08) and recipes moved to `verl-recipe` (2026/01). — [verl README](https://github.com/verl-project/verl)
- verl `examples/grpo_trainer` has VLM scripts for Qwen2.5-VL-7B and Qwen3-VL 8B/30B-A3B/235B-A22B. `verl/models/transformers/` has patches for `qwen2_vl`, `qwen3_vl`, `qwen3_5`, `glm4v`, `kimi_vl`, but none for Idefics3/SmolVLM. — GitHub API listing of [verl repo](https://github.com/verl-project/verl/tree/main/examples/grpo_trainer)
- verl custom reward: `custom_reward_function.path` / `.name`, signature `(data_source, solution_str, ground_truth, extra_info)`. — [verl docs](https://verl.readthedocs.io/en/latest/preparation/reward_function.html). In `naive.py`, if the score is a dict it takes `score["score"]` and appends the other keys to `reward_extra_info`. Reward managers: `naive`, `batch`, `dapo`, `prime`. The naive manager's comment notes it "scores samples serially in the driver process". — [verl reward_manager/naive.py](https://github.com/verl-project/verl/blob/main/verl/workers/reward_manager/naive.py)
- ms-swift GRPO: colocate mode (internal vLLM) or async/external mode via `swift rollout`. LoRA weight sync (`--vllm_enable_lora true`) and Qwen2.5-VL/Qwen3-VL ViT LoRA sync. `reward_weights`; logs `reward/{reward_func_name}/mean` and `/std`. — [ms-swift GRPO docs](https://swift.readthedocs.io/en/latest/Instruction/GRPO/GetStarted/GRPO.html). README: 400+ multimodal models, multimodal packing (+100% speed claimed). — [ms-swift README](https://github.com/modelscope/ms-swift)
- Unsloth VLM RL: Qwen3-VL, Qwen2.5-VL-7B, Gemma-3-4B; claims 1.5–2× faster and 90% less VRAM than FA2 setups. With `fast_inference=True` (vLLM), "vLLM does not support LoRA for vision/encoder layers", so `finetune_vision_layers=False` is required. — [Unsloth VLM RL docs](https://unsloth.ai/docs/new/vision-reinforcement-learning-vlm-rl)
- OpenRLHF 0.10 (2026-04) added VLM RLHF (e.g. Qwen3.5) and multi-turn VLM RL. LMM-R1 is a multimodal fork. — [OpenRLHF README](https://github.com/OpenRLHF/OpenRLHF)
- AReaL: Qwen2.5-VL, Qwen3-VL, Qwen3.5 (FSDP+Megatron), Gemma 3 (FSDP). Algorithm list GRPO, GSPO, PPO, DAPO, LitePPO, Dr.GRPO, REINFORCE++, RLOO, SAPO, IcePop, KPop, M2PO; all have async and sync (`max_head_offpolicyness=0`) versions. VLM examples are Geometry3K and CLEVR Count. — [AReaL README](https://github.com/inclusionAI/AReaL)
- ROLL: Qwen2.5-VL RLVR pipeline (2025-06), LoRA, GSPO, FSDP2 + Megatron LoRA (2026-02). A listed downstream project, SocioReasoner, does RL-optimized segmentation. — [ROLL README](https://github.com/alibaba/ROLL)
- NeMo-RL: VLM SFT and GRPO; Qwen3.5 dense/MoE LLM+VLM GRPO (2026-04); LoRA GRPO on DTensor and Megatron (2026-02); GDPO listed under algorithms; single-GPU GRPO examples (for LLMs). — [NeMo-RL README](https://github.com/NVIDIA-NeMo/RL)
- prime-rl: multimodal support for Qwen3.5 VLMs. "Need at least one NVIDIA GPU". The single-H100 example is a text model. — [prime-rl README](https://github.com/PrimeIntellect-ai/prime-rl)
- SkyRL: VLM RL support with a Qwen3-VL-8B GRPO example. — [Anyscale blog](https://www.anyscale.com/blog/vision-language-model-reinforcement-learning-skyrl); [AWS blog](https://aws.amazon.com/blogs/machine-learning/accelerate-multimodal-rl-training-with-skyrl-on-amazon-sagemaker-hyperpod/)
- slime README mentions no VLM support. Multimodal is via derived frameworks (Relax, Miles). — [slime README](https://github.com/THUDM/slime)
- LLaMA-Factory covers SFT/RM/PPO/DPO/KTO/ORPO and Qwen3-VL (2B–235B, template `qwen3_vl`). It announced EasyR1 for GRPO (2025-02-24). — [LLaMA-Factory README](https://github.com/hiyouga/LLaMA-Factory)

### Inferences
- **Per-term logging and GDPO matter most for the experiment design.** A design with about 10 reward arms × seeds and a controller driven by per-error-type logs needs each reward term (IoU, disagreement band, topology, legality, hole matching) to be its own reward function. That gives free per-term mean/std logging in TRL, or dict keys in EasyR1/verl. TRL's `normalize_then_sum` lets one arm be "GDPO" against "sum-then-normalize GRPO" with a single flag. That is a clean ablation axis no other framework offers as a flag (NeMo-RL lists GDPO, unverified detail).
- **Controller hooks:** no framework has a documented API for changing reward weights or KL mid-run (the TRL PR for that was closed). In TRL the practical route is a `TrainerCallback` that holds a reference to the trainer and mutates `trainer.reward_weights` / `trainer.beta` / the sampler's weights at `on_step_begin`, or a reward function that reads controller state from a shared object (rewards receive `trainer_state`). This is **unverified** and must be smoke-tested (see hypotheses). Pitfall: with `beta=0.0` at init, no reference model is loaded, so a controller that may raise KL later must start with a small non-zero beta (or LoRA, where the reference is the adapter-disabled base, according to TRL behaviour from prior knowledge; not re-verified here).
- **CPU reward throughput:** verl's naive manager scores serially in the driver. TRL calls each reward function once per batch, so the geometry checker should do its own `multiprocessing`/`concurrent.futures` pool inside one batched reward function (or use `async def` with `run_in_executor`) and then split the components into separate functions through a per-batch cache. EasyR1's `batch` reward type allows the same.
- **SmolVLM in RL frameworks:** only TRL (and Unsloth through TRL) handles SmolVLM/Idefics3 at all, and with vLLM it is currently broken. For SmolVLM arms, run with `use_vllm=False` (HF generate or `use_transformers_continuous_batching`) or apply the PR #6718 patch locally. With 256M/500M models, HF generation may be acceptable throughput-wise; this is untested.
- EasyR1 is in a slowdown: no release since 2025-09, and only docs commits since 2026-07-30. verl is the actively maintained upstream.

### Gaps
- Could not verify ms-swift's full GRPO algorithm list (DAPO/GSPO/CISPO/SAPO flags) or its reward plug-in API from the page fetched. My expectation that it supports these is from memory, not verified.
- Could not verify OpenRLHF, ROLL, SkyRL, prime-rl and NeMo-RL reward-function APIs or per-component logging; README-level evidence only.
- No source found for anyone running GRPO on SmolVLM-256M/500M (searches returned only model announcements).
- Whether mutating `trainer.reward_weights` / `trainer.beta` in a callback takes effect in TRL v1.16 is not verified.

---

## SFT tooling for small VLMs

### Takeaway
TRL `SFTTrainer` is the simplest choice because it uses the same stack and data format as GRPO (conversational `messages` + `image(s)`). It works with any `AutoModelForImageTextToText` VLM, including SmolVLM, LFM2.5-VL and Qwen3-VL. The main caveats are `max_length=None` for VLMs (truncation can cut image tokens) and `assistant_only_loss` for vision needing transformers ≥5.18. ms-swift advertises multimodal packing (+100% speed). LLaMA-Factory has a mature `qwen3_vl` SFT path. Liquid AI ships TRL and Unsloth SFT notebooks for LFM2.5-VL.

### Cited Findings
- TRL SFTTrainer "fully supports training Vision-Language Models" via an `image` or `images` column. The default collator is `DataCollatorForVisionLanguageModeling`. For VLMs "truncating may remove image tokens… set `max_length=None`". When the dataset has images, `skip_prepare_dataset` is forced True (on-the-fly preprocessing). — [TRL SFT docs](https://huggingface.co/docs/trl/main/en/sft_trainer)
- TRL SFT options: `packing` (strategies `bfd`, `bfd_split`, `wrapped`); `padding_free` (FlashAttention 2/3 only); `completion_only_loss`; `assistant_only_loss` ("With a vision dataset, assistant-only loss requires `transformers>=5.18.0`"); `loss_type="dft"`; `use_liger_kernel`; `activation_offloading`; `peft_config` (typical LoRA LR ≈1e-4). Defaults: `max_length=1024`, `learning_rate=2e-5`, `gradient_checkpointing=True`. — [TRL SFT docs](https://huggingface.co/docs/trl/main/en/sft_trainer)
- ms-swift: "multimodal packing technology to improve training speed by 100%+" and independent control of ViT/aligner/LLM. — [ms-swift README](https://github.com/modelscope/ms-swift)
- LLaMA-Factory supports Qwen3-VL 2B/4B/8B/30B/32B/235B (template `qwen3_vl`). VLM-R1 used LLaMA-Factory for its SFT baseline (`llamafactory-cli train …qwen2_5_vl_full_sft.yaml`). — [LLaMA-Factory README](https://github.com/hiyouga/LLaMA-Factory); [VLM-R1 README](https://github.com/om-ai-lab/VLM-R1)
- LFM2.5-VL-450M ships "SFT with LoRA using Unsloth" and "SFT with LoRA using TRL" notebooks. — [LFM2.5-VL-450M card](https://huggingface.co/LiquidAI/LFM2.5-VL-450M)
- SmolVLM-256M "supports fine-tuning via transformers" with a notebook tutorial. — [SmolVLM-256M card](https://huggingface.co/HuggingFaceTB/SmolVLM-256M-Instruct)

### Inferences
- Polygon JSON is long and variable, and images cost a fixed number of tokens (64 per 512 px patch for SmolVLM; about 256 for LFM2.5-VL and InternVL tiles). Padding waste is therefore dominated by completion-length variance. Packing for VLMs in TRL is not explicitly documented as supported with images. Start without packing and use `group_by_length`-style sampling; test packing on the Mac.
- Use the same chat template and prompt format in SFT and GRPO, and keep `max_length=None` in SFT. Set GRPO `max_completion_length` from the SFT polygon-length distribution (e.g. 99th percentile) and turn on `mask_truncated_completions=True`.
- Use one stack for both stages (TRL SFT → TRL GRPO, or LLaMA-Factory SFT → EasyR1 GRPO for Qwen-only). That avoids checkpoint/template conversion bugs.

### Gaps
- No benchmark found comparing SFT throughput across TRL / ms-swift / LLaMA-Factory / Unsloth for sub-1B VLMs.
- Whether TRL `packing=True` works with image datasets (SmolVLM/Qwen3-VL) is not documented on the page fetched.

---

## Smallest VLMs that could learn to output polygons / coordinates

### Takeaway
Qwen3-VL-2B (Apache-2.0, 0–1000 relative coordinate JSON, tested in TRL, supported by EasyR1/verl/vLLM) is the safest primary model. LFM2.5-VL-450M is the strongest **sub-0.5B** candidate: native normalized-bbox JSON grounding (RefCOCO-M 81.28), vLLM with LoRA, and TRL/Unsloth notebooks, but under a custom LFM1.0 licence. SmolVLM-256M/500M are the cheapest (Apache-2.0) but have no grounding prior and hit the TRL+vLLM bug. Florence-2 (0.23B/0.77B, MIT) natively outputs polygons for referring segmentation, but it is encoder-decoder, so no surveyed GRPO framework supports it (baseline only). Qwen3.5-0.8B (Feb–Mar 2026, natively multimodal, RefCOCO 79.3) is the newest and most interesting sub-1B option, with kernel/framework-version risk.

### Model table

| Model | Params | Image handling | Native grounding / coord convention | Licence | vLLM | RL / framework support |
|---|---|---|---|---|---|---|
| SmolVLM-256M-Instruct | 256M (SigLIP 93M + SmolLM2-135M; Idefics3 arch) | 512×512 patches, 64 tokens/patch; up to `longest_edge=N*512` | None documented | Apache-2.0 | Yes (`Idefics3ForConditionalGeneration`, LoRA ✓) | TRL only (vLLM path bugged) |
| SmolVLM-500M-Instruct | ~500M (same 93M encoder) | as above | None | Apache-2.0 | Yes | TRL only |
| SmolVLM2-2.2B | 2.2B | — | None | Apache-2.0 (assumed, not verified) | Yes (`SmolVLMForConditionalGeneration`) | TRL GRPO-tested |
| LFM2.5-VL-450M | 450M (SigLIP2 86M + LFM2.5-350M) | Native 512×512, tiling + thumbnail, ≤256 image tokens by default | bbox JSON `[{"label","bbox":[x1,y1,x2,y2]}]`, normalized [0,1]; RefCOCO-M 81.28 | LFM1.0 (custom) | Yes (`Lfm2VlForConditionalGeneration`, LoRA ✓) | TRL / Unsloth SFT notebooks; GRPO via TRL (untested) |
| Qwen3.5-0.8B | 0.8B, natively multimodal, Gated DeltaNet hybrid | — | RefCOCO 79.3 reported (coordinate convention not on card) | Apache-2.0 | Needed vLLM main at release | AReaL, prime-rl, NeMo-RL, OpenRLHF, verl `qwen3_5` patch |
| Qwen3-VL-2B-Instruct | 2B | Dynamic resolution | JSON `bbox_2d` / `point_2d`, **relative 0–1000** | Apache-2.0 | Yes (`Qwen3VLForConditionalGeneration`, LoRA ✓) | TRL (tested), EasyR1, verl, Unsloth, ms-swift, LLaMA-Factory |
| InternVL3.5-1B | 1.1B (InternViT 0.3B + 0.8B LM) | 448 px tiles, 1024 → 256 tokens via pixel shuffle | Box format not on card | Apache-2.0 | Yes (`InternVLChatModel`) | Trained with Cascade RL (MPO + GSPO); VLM-R1 supports InternVL |
| Florence-2-base / -large | 0.23B / 0.77B, encoder-decoder | — | `<loc_*>` quantized tokens, 1000 bins; polygon output for referring-expression segmentation | MIT | Only via official **bart-plugin** (not native V1) | No GRPO framework support found |
| Moondream 3 (Preview) | Not verified | — | Native **SVG-path** segmentation + bbox; RefCOCO val 83.2 mIoU (2026-03 update) | Moondream licence (terms not verified) | Yes (`Moondream3ForCausalLM`, no LoRA listed) | None found |
| FastVLM-0.5B (Apple) | Card says 0.8B | — | No grounding documented | `apple-amlr` | Card shows `vllm serve` (not in vLLM's list I fetched) | None |
| Gemma 3n | Not verified | — | — | Gemma terms | Yes (no LoRA, needs timm≥1.0.17) | — |
| PaliGemma 2 | Smallest 3B (from memory) | — | Seg tokens (from memory) | Gemma | Yes (LoRA ✓) | — |

### Cited Findings
- SmolVLM-256M: 93M SigLIP + SmolLM2-135M on the Idefics3 architecture; 512×512 patches at 64 visual tokens each; up to 2048×512 via `size={"longest_edge": N*512}`; "inference on one image with under 1GB of GPU RAM"; Apache 2.0; released April 2025 (arXiv 2504.05299); vLLM and SGLang supported. — [SmolVLM-256M card](https://huggingface.co/HuggingFaceTB/SmolVLM-256M-Instruct)
- SmolVLM-256M/500M use the SigLIP base patch-16/512 (93M) encoder instead of the 400M SO encoder. — [HF blog "SmolVLM grows smaller"](https://huggingface.co/blog/smolervlm)
- vLLM supported multimodal list includes `SmolVLMForConditionalGeneration` (LoRA ✓), `Idefics3ForConditionalGeneration` (LoRA ✓), `PaliGemmaForConditionalGeneration`, `Lfm2VlForConditionalGeneration`, `InternVLChatModel`, `Qwen3VLForConditionalGeneration` (all LoRA ✓), `Gemma3nForConditionalGeneration` (needs `timm>=1.0.17`), `Moondream3ForCausalLM`, `MiniCPMV`, `KimiVLForConditionalGeneration`. Florence-2 and FastVLM are not in the main list. — [vLLM supported models](https://docs.vllm.ai/en/latest/models/supported_models/)
- vLLM V1: "Whisper is supported natively. Other encoder-decoder models are supported via the plugin system"; "`Florence2ForConditionalGeneration` is supported via the official bart-plugin". — [vLLM V1 guide](https://docs.vllm.ai/en/latest/usage/v1_guide/)
- LFM2.5-VL-450M: 450M params; card says released 2025-11-28 (third-party summary says it refreshes LFM2-VL-450M on the LFM2.5-350M backbone with an 86M SigLIP2 encoder); native 512×512 with tiling + thumbnail; ≤256 image tokens default; bbox JSON with normalized [0,1] coordinates; RefCOCO-M 81.28 ("new capability, LFM2-VL-450M lacked this"); LFM1.0 licence; vLLM compatible. — [LFM2.5-VL-450M card](https://huggingface.co/LiquidAI/LFM2.5-VL-450M); [llmrun.dev summary](https://llmrun.dev/model/liquidai-lfm2-5-vl-450m)
- Qwen3-VL grounding uses relative coordinates 0–1000 (changed from Qwen2.5-VL's absolute pixels) and JSON `[{"bbox_2d":[x1,y1,x2,y2],"label":...}]`, plus `point_2d` for points. — [Qwen3-VL grounding docs (mintlify mirror)](https://mintlify.wiki/QwenLM/Qwen3-VL/capabilities/grounding-2d); [QwenLM/Qwen3-VL issue #1616](https://github.com/QwenLM/Qwen3-VL/issues/1616)
- Qwen3-VL-2B-Instruct: Apache 2.0, native 256K context; the card recommended transformers from source at release. — [Qwen3-VL-2B card](https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct)
- Qwen3.5-0.8B: Feb 2026, Apache 2.0, "Causal Language Model with Vision Encoder", RefCOCO 79.3/77.8, Gated DeltaNet + gated attention hybrid, vLLM/SGLang main branch required at release. — [Qwen3.5-0.8B card](https://huggingface.co/Qwen/Qwen3.5-0.8B). Small series (0.8B/2B/4B/9B) released early March 2026. — [Gigazine](https://gigazine.net/gsc_news/en/20260303-qwen-3-5-small)
- InternVL3.5-1B: 0.3B vision + 0.8B LM, released 2025-08-25, 448 px tiles at 256 tokens, Apache 2.0, vLLM/LMDeploy, Cascade RL (MPO offline + GSPO online). — [InternVL3.5-1B card](https://huggingface.co/OpenGVLab/InternVL3_5-1B)
- Florence-2: base 0.23B, large 0.77B, MIT, `trust_remote_code=True` for the Microsoft repo, 1000-bin quantized location tokens. — [Florence-2-base card](https://huggingface.co/microsoft/Florence-2-base)
- Moondream 3 Preview introduced segmentation as native SVG paths (Sept 2025). The 2026-03-10 update reports RefCOCO val 83.2, RefCOCO+ 79.1, RefCOCOg 80.7, RefCOCO-M 88.2 mIoU, and says local weights would be released. — [Moondream blog 2026-03-10](https://moondream.ai/blog/segmenting-update-2026-03-10); [Moondream segment docs](https://docs.moondream.ai/skills/segment)
- FastVLM-0.5B: `apple-amlr` licence, `trust_remote_code=True`, CVPR 2025; the card lists 0.8B params. — [FastVLM-0.5B card](https://huggingface.co/apple/FastVLM-0.5B)
- Related RL-for-segmentation work (all ≥3B):
  - **SimpleSeg** (Moonshot/NJU, arXiv 2601.19228, Jan 2026): the MLLM predicts polygon point sequences as text, trained SFT→RL with an IoU reward. Models are Qwen2.5-VL-7B and Kimi-VL-16B-A3B, with normalized 0–1 polygon coords. SimpleSeg-Kimi-VL reports 74.8 avg RES and 86.5 REC. — [SimpleSeg GitHub](https://github.com/songtianhui/SimpleSeg); [HF card](https://huggingface.co/sthui/SimpleSeg-Kimi-VL/blob/main/README.md); [arXiv](https://arxiv.org/abs/2601.19228)
  - **GenSeg-R1** (Feb 2026): Qwen3-VL-4B/8B + GRPO producing box + 2 keypoints for frozen SAM 2; RefCOCOg 0.7127 cIoU (8B). — [arXiv 2602.09701](https://arxiv.org/abs/2602.09701)
  - **Seg-Zero / VisionReasoner**: built on EasyR1/veRL, Qwen2-VL/Qwen2.5-VL + SAM2, IoU and L1 rewards; 7B recommended on 4×80 GB. — [Seg-Zero README](https://github.com/dvlab-research/Seg-Zero)
  - **VLM-R1**: Qwen2.5-VL-3B REC/OVD; reports reward hacking in OD; LoRA GRPO; LLaMA-Factory SFT baseline. — [VLM-R1 README](https://github.com/om-ai-lab/VLM-R1); [arXiv 2504.07615](https://arxiv.org/abs/2504.07615)

### Inferences
- No published GRPO result exists for any VLM under 2B producing polygons. The project would be first-of-kind at the 256M–500M scale, which supports the planned SFT-first warm start. SimpleSeg is the closest recipe (text polygons, SFT→RL, IoU reward) and validates the approach at 7B+.
- Coordinate convention choice: Qwen3-VL's 0–1000 integers are token-efficient and match its pretraining. LFM2.5-VL uses [0,1] floats, which cost more tokens per vertex. SmolVLM has no prior, so our own integer grid (e.g. 0–999) can be chosen freely. Tokenization cost per vertex should be measured per tokenizer (see hypotheses) because rollout cost scales with it.
- Florence-2 is a strong **non-RL baseline** (zero-shot polygons, 0.23B). It cannot be the GRPO policy in TRL: TRL trainers state "Only causal language models are supported" for model objects in SFT. GRPO is assumed to have the same restriction (not verified).
- Qwen3.5-0.8B needs flash-linear-attention / causal-conv1d-style kernels for Gated DeltaNet. On Mac CPU/MPS it would use slow fallbacks, and framework support skews to large-scale frameworks. Treat it as a stretch arm.

### Gaps
- Qwen3-VL-2B per-image token count (patch size / merge) not confirmed from the card fetched. From memory it is 16 px patches with 2×2 merge (32×32 px per token), unverified.
- SmolVLM-500M LM size (SmolLM2-360M) is from memory, not re-verified. Moondream 3 param count and licence, Gemma 3n param counts, and PaliGemma 2 sizes/seg format are not verified.
- The LFM2.5-VL-450M release date (card: 2025-11-28) could not be cross-checked.
- No vLLM benchmark found for LFM2.5-VL / SmolVLM rollouts.

---

## Rollout throughput (vLLM) and GRPO memory on 4090 / A100

### Takeaway
No direct published number exists for 0.5–2B **VLM** GRPO rollouts on a 4090. The nearest proxies are vLLM offline text-only throughput on an RTX 4090 of about 8.3k output tok/s (1.5B) and about 6.2k (3B), and about 3.5k for Qwen2.5-VL-7B. Memory is not the constraint for ≤500M models. For Qwen3-VL-2B, full-parameter GRPO fits on 1×24 GB only in pure BF16 (EasyR1 estimate for 1.5B), and LoRA is comfortable. An 80 GB card allows full FT plus colocated vLLM with headroom.

### Cited Findings
- RTX 4090, vLLM, 300 concurrent requests, 100 input / 600 output tokens (page updated 2026-10-07, vLLM version not stated): DeepSeek-R1-Distill-Qwen-1.5B 8,311 output tok/s offline (6,981 online); Qwen2.5-3B-Instruct 6,184 offline (5,922 online); Qwen2.5-VL-7B-Instruct 3,476 output tok/s offline. — [DatabaseMart RTX4090 vLLM benchmark](https://databasemart.com/blog/vllm-gpu-benchmark-rtx4090)
- Qwen3-VL-2B-Instruct about 2,409 tok/s on an RTX PRO 6000 Blackwell with vLLM (search snippet only; page not fetched; setup unknown). — [benchmarklist.com vlmbench_perf](https://benchmarklist.com/benchmarks/vlmbench_perf/)
- EasyR1 GRPO memory estimates: 1.5B full BF16 1×24 GB; 1.5B LoRA 1×12 GB; 3B LoRA 1×24 GB; 3B full BF16 1×40 GB. — [EasyR1 README](https://github.com/hiyouga/EasyR1)
- TRL: `beta=0.0` by default, so no reference model is loaded ("faster, less memory"). `vllm_gpu_memory_utilization=0.3` by default in colocate mode. — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- SmolVLM-256M inference on one image fits in under 1 GB of GPU RAM. — [SmolVLM-256M card](https://huggingface.co/HuggingFaceTB/SmolVLM-256M-Instruct)

### Inferences
- Rough step-time estimate (assumption, not measured): 64 prompts × G=8 = 512 rollouts × ~400 polygon tokens ≈ 200k tokens per step. At ~4–8k tok/s colocated on a 4090, that is roughly 25–50 s of generation per step for a 0.5–2B model, plus image prefill. The CPU geometry reward on 512 samples must stay at a few seconds to avoid becoming the bottleneck. Use a process pool on the rented box's many vCPUs.
- Memory arithmetic (assumption): Qwen3-VL-2B full FT with AdamW needs about 2B × (2 bf16 weights + 2 grads + 8 fp32 Adam states) ≈ 24 GB before activations and vLLM. On a 4090 use LoRA (or the EasyR1 pure-bf16 optimizer), and on A100/H100 80 GB full FT plus vLLM at about 0.3 utilization fits. A KL reference (beta>0) adds about 4 GB bf16 for a 2B model, or about 0 with LoRA (adapter-disabled base).
- For 256M/500M, a 4090 can run several seeds concurrently (multiple processes per GPU), which helps the 10 arms × seeds design. Generation will likely be the bottleneck, not training.

### Gaps
- No measured vLLM throughput for Qwen3-VL-2B, LFM2.5-VL-450M or SmolVLM on a 4090/A100 with image inputs. A 15-minute benchmark on the rented box is needed.
- No published GRPO memory numbers for sub-1B VLMs.

---

## Experiment tracking / config / data tooling (2026)

### Takeaway
W&B remains universally supported. Trackio (HF, local-first, SQLite + optional HF Space) is the recommended tracker for TRL and is natively integrated (`report_to="trackio"`). EasyR1 supports W&B, SwanLab, MLflow and TensorBoard. Configs are YAML/OmegaConf (EasyR1, verl-family) or `HfArgumentParser` dataclasses (TRL). Data is HF `datasets` (Arrow/Parquet) with an `image(s)` column of PIL images.

### Cited Findings
- Trackio is "the recommended tracking solution for TRL", natively integrated with all trainers via `report_to="trackio"`. It stores locally in SQLite and syncs to an HF Space; latest 0.25.1 (2026-04-26). — [TRL Trackio integration](https://huggingface.co/docs/trl/main/en/trackio_integration); [CERN PyPI mirror listing](https://simple-repository.app.cern.ch/project/trackio)
- TRL `log_completions=True` logs (prompt, completion) samples to wandb and/or trackio, `log_multimodal` includes images, and `log_completions_hub_repo` pushes completions to the Hub. SFTConfig exposes `trackio_space_id` / `trackio_bucket_id`. — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer); [TRL SFT docs](https://huggingface.co/docs/trl/main/en/sft_trainer)
- EasyR1: "Wandb & SwanLab & Mlflow & Tensorboard tracking"; YAML config at `examples/config.yaml`. — [EasyR1 README](https://github.com/hiyouga/EasyR1)
- TRL GRPO/SFT datasets: conversational `prompt` + `image`/`images` columns (PIL). — [TRL GRPO docs](https://huggingface.co/docs/trl/main/en/grpo_trainer)
- verl added RL-Insight (2026/07) for online observability of RL training. — [verl README](https://github.com/verl-project/verl)

### Inferences
- For 10 arms × seeds, W&B groups/sweeps or Trackio projects work, and TRL's per-reward `rewards/{name}/mean` keys map directly onto the controller's per-error-type signals. Log the raw per-sample component vectors with `log_extra`, written to a local Parquet file each step, so the controller and offline analysis read the same data.
- Store the dataset as Parquet via HF `datasets` with images embedded (or as paths) and ground-truth geometry as WKT/GeoJSON strings. WebDataset is unnecessary at this scale.

### Gaps
- verl's config system: I believe it is Hydra, from memory; not verified this session.
- No 2026 survey found on tracker market share.

---

## Candidate hypotheses to test locally on the Mac (CPU/MPS, no CUDA)

Local environment observed: system `python3` is **3.14.4** with no torch/transformers/trl installed. Use `uv venv --python 3.12` (assumption: 3.12 has the broadest torch/transformers/trl wheel coverage). vLLM, flash-attn, EasyR1 and verl are CUDA-centric and should not be expected to run on the Mac.

1. **H1: TRL GRPO + SmolVLM-256M runs end to end on CPU/MPS without vLLM.** `pip install trl peft transformers accelerate datasets pillow`. Run 2 GRPO steps: `use_vllm=False`, `num_generations=2`, `per_device_train_batch_size=2`, `max_completion_length=32`, fp32 (MPS bf16 can be flaky), 4 synthetic images, 3 dummy reward functions. *Expected:* completes in minutes and logs `rewards/<fn>/mean` for each function. If MPS errors, `use_cpu=True` works.
2. **H2: `multi_objective_aggregation="normalize_then_sum"` is robust to a zero-variance term.** Make one reward constant within a group (e.g. legality always 1). *Expected:* no NaN (std epsilon), and that term contributes 0 advantage. Compare advantages printed under `sum_then_normalize`.
3. **H3: Controller hook works via a TrainerCallback.** Mutate `trainer.reward_weights` and `trainer.beta` (initialized at 0.001 so the ref model exists) at step 1. *Expected:* the change appears in the logged `reward` and `kl`. If `beta` starts at 0.0, raising it later fails or does nothing (no ref model). This decides whether we need to subclass `GRPOTrainer`.
4. **H4: Curriculum by sampler re-weighting.** Override the train sampler (or wrap the dataset with a mutable weight table read by a custom `RepeatSampler`) and change weights mid-run. *Expected:* feasible but requires subclassing `_get_train_sampler`. Check that grouping of `num_generations` per prompt is preserved.
5. **H5: Async/batched CPU reward does not block.** Implement the geometry checker as one batched function using a `ProcessPoolExecutor`, exposing components via cached per-term functions. Time 512 synthetic polygons. *Expected:* about 1 s or less on M3 for simple polygons; confirms the reward will not bottleneck a GPU step.
6. **H6: Tokens per polygon by tokenizer.** Serialize a 40-vertex polygon with one hole in candidate formats (0–999 ints, 0–1 floats with 3 decimals, Qwen-style `[[x,y],...]`, compact delta encoding). Count tokens with the SmolVLM, LFM2.5-VL, Qwen3-VL and Qwen3.5 tokenizers (tokenizers only, no weights needed). *Expected:* integer grids need about 2–4× fewer tokens than floats. Pick the format that minimizes rollout length.
7. **H7: Image tokens per image.** Run the processors for SmolVLM-256M (`longest_edge=512` vs 1024), LFM2.5-VL-450M and Qwen3-VL-2B (`min/max_pixels`) on a 512×512 image. *Expected:* SmolVLM 64 tokens per 512 px patch (plus a global patch if split), LFM about 256 or fewer, Qwen3-VL about 256 at 512² (if 32 px per token). This drives prefill cost and SFT sequence length.
8. **H8: TRL SFT on 20 synthetic polygon samples with SmolVLM-256M (LoRA and full).** `max_length=None`, `assistant_only_loss` off. Overfit until the model emits parseable JSON. *Expected:* loss falls and greedy outputs parse after about 100–300 steps on CPU/MPS. Then test `packing=True` with images; it may error or be unsupported.
9. **H9: Zero-shot priors.** Run LFM2.5-VL-450M and Qwen3-VL-2B (fp32/bf16 on MPS, ~5–9 GB) on 5 images with a bbox prompt, and Florence-2-base with `<REFERRING_EXPRESSION_SEGMENTATION>`. *Expected:* Qwen3-VL and LFM emit valid bbox JSON; Florence-2 emits polygons; SmolVLM-256M emits no usable coordinates. This confirms SFT is required for SmolVLM and the size of the "prior advantage" for the others.
10. **H10: Reproduce the TRL SmolVLM+vLLM bug logic offline.** Run the diagnostic from issue #6294 (processor-expanded vs tokenizer-only prompt lengths) for SmolVLM-256M and Qwen3-VL-2B. It needs no vLLM. *Expected:* SmolVLM expanded prompt has many `<image>` runs, while Qwen has one contiguous `<|image_pad|>` run. This confirms SmolVLM arms need `use_vllm=False` or the PR #6718 patch on the GPU box.
11. **H11: Reward-function portability.** Write the checker adapter once and wrap it for TRL (`fn(completions, **kw) -> list[float]`), EasyR1 (`compute_score(reward_input) -> dict`) and verl (`compute_score(data_source, solution_str, ground_truth, extra_info) -> dict` with `score` key). Unit-test all three wrappers locally. *Expected:* trivial, and it keeps the option to move Qwen arms to EasyR1/verl if TRL throughput is insufficient.

## Recommendation

**Framework (ranked):**
1. **TRL `GRPOTrainer` + `SFTTrainer` (pin v1.16.0; watch PRs #6718 and #7698).** Single stack for SFT and RL, and the only framework with explicit Qwen3-VL-2B and SmolVLM tests. Per-reward logging and `log_extra`/`log_metric` from inside reward functions; sync/async batched Python rewards. Every algorithm arm needed is a config flag: `loss_type` dr_grpo/dapo/cispo/sapo, GSPO via sequence-level IS, `scale_rewards`, GDPO via `normalize_then_sum`. It runs on the Mac for smoke tests and uses standard callbacks for the Level 2 controller. Use vLLM colocate for Qwen3-VL/LFM arms. For SmolVLM arms use HF generation / continuous batching, or locally apply the PR #6718 patch, until the fix merges.
2. **EasyR1** (or upstream **verl** v0.9.x) as the throughput fallback **for Qwen3-VL-2B only**. It has a proven segmentation-RL lineage (Seg-Zero, VisionReasoner), a dict-valued multi-component reward, and documented 1×24 GB GRPO for 1.5B. It does not support SmolVLM or LFM, controller hooks require forking, and EasyR1 has had no release since 2025-09. Choose verl over EasyR1 if going this route long term.
3. **ms-swift / Unsloth** as alternatives. ms-swift has the broadest model zoo and multimodal packing for SFT (GRPO details unverified). Unsloth offers memory savings for Qwen3-VL LoRA GRPO but freezes vision layers under vLLM and has no SmolVLM support. OpenRLHF, AReaL, ROLL, NeMo-RL, prime-rl, SkyRL and slime are aimed at multi-GPU/Megatron/async scale and Qwen-family VLMs, so they are poor fits for single-GPU sub-1B VLM experiments.

**Models (ranked):**
1. **Qwen3-VL-2B-Instruct.** Apache-2.0 with a native 0–1000 JSON grounding prior. Supported everywhere (TRL-tested, EasyR1, verl, vLLM+LoRA). LoRA on a 4090, full FT on an 80 GB card.
2. **LFM2.5-VL-450M.** Best sub-0.5B option: native normalized-bbox JSON grounding (RefCOCO-M 81.28), vLLM+LoRA, TRL/Unsloth SFT notebooks. Check the LFM1.0 licence terms before any commercial use.
3. **SmolVLM-500M / SmolVLM-256M.** Cheapest Apache-2.0 option for many seeds × arms. No grounding prior, so SFT does all the work, and the vLLM path in TRL is currently bugged. Good for a model-size scaling axis.
4. **Qwen3.5-0.8B** (stretch). Newest, natively multimodal, RefCOCO 79.3, but carries DeltaNet kernel and framework-version risk.
5. **InternVL3.5-1B** (alternative mid-size option).
6. **Florence-2-base/large** as a non-RL polygon baseline (native polygon output, MIT). It is not trainable with the surveyed GRPO frameworks.
