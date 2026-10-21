# A2: SFT recipe and software stack for Phase 2 (as of 2026-10-19)

Question: what SFT recipe and stack teaches SmolVLM-256M, LFM2.5-VL-450M and Qwen3-VL-2B to emit polygon JSON (D10)
in about 1 GPU-hour per run on one Vast GPU, so that the checkpoint is a good GRPO start for Phases 3–4 (TRL 1.16,
transformers 5.20, peft 0.21, torch 2.14)?

Method: primary sources (framework docs and source, release pages, papers 2025–2026), the installed TRL/transformers
source in `.venv`, and four cheap CPU checks on `data/prepared` (scripts and outputs in the A2 scratch directory:
`seq_budget.py`/`seq_budget.json`, `eos_mask.py`, `ntl_proto.py`, `collate_time.py`). No MPS training was run
(wave B). The session's web-search budget ran out near the end; the gaps that left are listed in §10.

## 0. Bottom line

1. **Stay on TRL `SFTTrainer`.** No other framework covers all three models *and* gives a VLM speed feature that we
   can use. TRL 1.16 refuses packing and padding-free for vision datasets. LLaMA-Factory, Axolotl, Unsloth and
   ms-swift each miss at least one of our models, or have no VLM packing either. The GRPO path is TRL, so the SFT
   checkpoint, chat template and loss mask stay identical across both stages.
2. **Padding is the only real waste, and a sampler removes it.** Use length-grouped sampling
   (`train_sampling_strategy="group_by_length"` with a precomputed `length` column). On our data it cuts padding
   from 7–32 % of tokens to ≤ 0.2 %. Packing would add almost nothing beyond that, even if TRL allowed it for VLMs.
3. **The CPU must keep up.** Collating a batch of 8 takes 110–350 ms on one core, about as long as a SmolVLM GPU
   step. Use `dataloader_num_workers ≥ 4`.
4. **The recipe:**
   - full fine-tuning for SmolVLM and LFM, with the vision tower unfrozen (to be tested);
   - LoRA r16 on all linear layers for Qwen3-VL-2B, merged after SFT;
   - assistant-only NLL; a warmup–stable–decay schedule; 1–3 epochs over all raters' targets.
5. **Pick the checkpoint by its GRPO signal, not by its SFT score.** Use decay branches at several depths. The
   evidence: past the point where SFT scores peak, RL gains fall, and pre-RL entropy and pass@k predict the RL
   outcome. Choose with the Phase 2 GRPO-signal gate (duplicates, σ above the floor), pass@8 and held-out NLL.
6. **One loss addition is worth testing:** a place-value-weighted number-token loss (NTL-WAS, λ = 0.3). Prototype
   result: computing it from hidden states at label positions reproduces TRL's fused NLL to 3e-5. Rejected
   alternatives: DFT (entropy collapse and drift), label smoothing (deprecated in TRL, and it does not repair
   over-trained checkpoints), consensus targets.
7. **Holes are rare.** Only 1.35 % of LIDC train targets have a hole and 0.85 % have several parts. A0′ needs
   targeted oversampling or synthetic hole shapes (hypothesis H5).

## 1. What the cheap checks on this repo showed

### 1.1 TRL 1.16 behaviour that constrains the stack (installed source, `trl/trainer/sft_trainer.py` v1.16.0)

| Fact | Where |
|---|---|
| `packing=True` with a vision dataset raises "Packing is not supported for vision datasets". `padding_free=True` raises "Padding-free training is yet not supported for vision datasets". | sft_trainer.py L779–787 ([v1.16.0 source](https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/sft_trainer.py)) |
| Images are processed on the fly in the collator (`DataCollatorForVisionLanguageModeling`). `skip_prepare_dataset` is forced for image datasets. | [TRL SFT docs](https://huggingface.co/docs/trl/main/en/sft_trainer) |
| The loss runs through TRL's own fused LM head (`add_fused_lm_head`): per-token log-probabilities in tiles at label positions only, with no full logits. It also logs `entropy` and `mean_token_accuracy`. | sft_trainer.py L1142–1147; trainer/utils.py `add_fused_lm_head` |
| `use_liger_kernel` is deprecated, and Liger's fused linear cross-entropy is refused. Only Liger's layer kernels are applied. Hub kernels are recommended instead (`model_init_kwargs={"use_kernels": True}`). | sft_trainer.py L1125–1140 |
| `label_smoothing_factor` is deprecated in SFT: it needs full logits, so each step runs an extra forward pass. `compute_loss_func` is deprecated; the advice is to subclass `compute_loss`. | sft_trainer.py L1049–1063 |
| `loss_type ∈ {"nll", "dft"}`; DFT multiplies each token's NLL by its detached probability. | sft_trainer.py L1613–1615; [DFT paper](https://arxiv.org/abs/2508.05629) |
| Assistant-only loss swaps in a bundled training template with `{% generation %}` markers when the model's own template lacks them. It also warns when the end-of-turn token is outside the mask. | sft_trainer.py L938–955; `chat_template_utils.is_chat_template_stop_token_trained` |
| `gradient_checkpointing=True` is the default. Selective activation checkpointing is available (`gradient_checkpointing_kwargs={"selective": True}`). | [TRL memory docs](https://huggingface.co/docs/trl/main/en/reducing_memory_usage) |
| transformers 5.20 provides `lr_scheduler_type="warmup_stable_decay"` (`get_wsd_schedule`) and `train_sampling_strategy ∈ {random, group_by_length, batch_rebalance}` with `length_column_name`. | [training_args.py v5.20.1](https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/training_args.py); `optimization.get_wsd_schedule` |
| `full_determinism=True` calls `enable_full_determinism`. That sets `CUDA_LAUNCH_BLOCKING=1` (serialised kernel launches), `CUBLAS_WORKSPACE_CONFIG=:16:8`, `FLASH_ATTENTION_DETERMINISTIC=1` and `torch.use_deterministic_algorithms(True)`. | [trainer_utils.py v5.20.1](https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/trainer_utils.py) |

### 1.2 Sequence budget per model and dataset (train split, every voting rater's target; `seq_budget.py`)

Image tokens come from the real prepared images with the D12 processor settings. The GPU-time column is an
**estimate**: 8·N·tokens FLOPs (6N for forward and backward plus 2N for checkpoint recompute) at an assumed 25 %
MFU, using datasheet bf16 dense peaks of 165 TFLOPS (RTX 4090) and 989 TFLOPS (H100). The first Vast session's
measurement replaces it (ROADMAP §Compute).

| Model / dataset | Targets | Completion tokens p50 / p99 | Prompt tokens (incl. image) | Padding, random batch of 8 | Padding, length-grouped | Tokens per epoch | Est. min/epoch, 4090 | Est. min/epoch, H100 |
|---|---|---|---|---|---|---|---|---|
| SmolVLM-256M / RIGA | 3,282 | 317 / 423 | 417 | 7.2 % | 0.15 % | 2.43 M | 2.0 | 0.3 |
| SmolVLM-256M / LIDC | 23,469 | 104 / 408 | 140 | **32.0 %** | 0.05 % | 6.37 M | 5.3 | 0.9 |
| LFM2.5-VL-450M / RIGA | 3,282 | 208 / 274 | 294 (150–342) | 6.6 % | 0.14 % | 1.66 M | 2.4 | 0.4 |
| LFM2.5-VL-450M / LIDC | 23,469 | 71 / 261 | 332 | 16.0 % | 0.02 % | 9.86 M | 14.3 | 2.4 |
| Qwen3-VL-2B / RIGA | 3,282 | 310 / 416 | 274 (140–332) | 8.8 % | 0.19 % | 1.94 M | 13.1 | 2.2 |
| Qwen3-VL-2B / LIDC | 23,469 | 101 / 405 | 325 | 22.0 % | 0.03 % | 10.64 M | 72.2 (full FT) | 12.0 |

What this implies:
- **SmolVLM and LFM:** one epoch over both datasets costs roughly 7–17 min on a 4090, even at a pessimistic MFU.
  The 1 GPU-hour budget therefore holds about 3–6 epochs, or a learning-rate sweep.
- **Qwen3-VL-2B:**
  - With full fine-tuning, one epoch takes about 85 min on a 4090, against about 14 min on an H100.
  - LoRA saves the weight-gradient matmuls, roughly a quarter of the FLOPs. One epoch on a 4090 is still more than
    an hour.
  - So in a 1-hour run, Qwen either goes to an H100 or A100, or trains on an LIDC subsample.
- **Prompt share:** 47–79 % of tokens are prompt (image plus text), which carries no loss. The completion is only
  21–53 % of the tokens.

### 1.3 The stop token is trained for all three models (`eos_mask.py`, TRL's own template swap)

With `assistant_only_loss=True` and TRL's template swap, the loss span ends on the end-of-turn token for every model:
- SmolVLM: `<end_of_utterance>`, bundled `idefics3_training` template;
- LFM2.5-VL: `<|im_end|>`, native generation markers;
- Qwen3-VL: `<|im_end|>`, bundled `qwen3_vl_training` template.

TRL's `is_chat_template_stop_token_trained` probe returns True for all three. Without the swap, the SmolVLM and Qwen
templates have no markers and give **0 loss tokens**. Any custom collator (for example a future NTL trainer or an
augmentation collator) must pass `chat_template=get_training_chat_template(proc)` the way SFTTrainer does.

### 1.4 Collation is a CPU bottleneck at this model size (`collate_time.py`; Mac under load average 31, so upper bounds)

| | RIGA | LIDC |
|---|---|---|
| SmolVLM-256M | 353 ms | 111 ms |
| LFM2.5-VL-450M | 213 ms | 231 ms |
| Qwen3-VL-2B | 297 ms | 233 ms |

The table gives the median time to collate one batch of 8 (PIL decode, processor, template, assistant mask) in a
single process. The estimated GPU time for one SmolVLM step on RIGA (batch 8 × 740 tokens) is about 0.3 s on a 4090,
so collation alone matches it. The default `dataloader_num_workers=0` would starve the GPU. Use 4–8 workers with
`dataloader_persistent_workers=True`, on a box with at least 8 cores (ROADMAP already asks for that for the reward).
Shared with A6.

### 1.5 Coordinate tokenisation (relevant to number-token losses; `ntl_proto.py`)

- SmolVLM and Qwen3-VL split every digit, using 10 digit tokens.
- LFM2.5 has a single token for each of 0–999. `1000` becomes `100`+`0`, but no train target contains 1000, so
  this does not matter.
- Digits make up **58.9 %** of SmolVLM loss tokens on LIDC.
- The NTL prototype:
  - computes CE from the backbone's hidden states at label positions only, and NTL-WAS from the 10 digit-token
    rows of `lm_head`, renormalised;
  - reproduces TRL's fused-head NLL: 1.57401 vs 1.57404, fp32;
  - gives a finite NTL of 1.40, gradients flow, and a one-hot correct distribution scores 0;
  - peaks at 3.2 GB RSS on CPU (batch of 4).

  An `SFTTrainer.compute_loss` override is therefore cheap to write and keeps TRL's no-full-logits design.

## 2. Frameworks

| Framework (latest release) | SmolVLM-256M | LFM2.5-VL-450M | Qwen3-VL-2B | VLM packing / padding-free | GRPO path | Verdict |
|---|---|---|---|---|---|---|
| **TRL** v1.16.0 (2026-10-15) | yes (bundled idefics3 training template) | yes | yes | **no** (refused for vision datasets) | same library (D8) | **keep** |
| Unsloth (v0.1.909-beta, 2026-10-19) | not listed | LFM2.5-VL-1.6B notebook only (LoRA, lr 2e-4, vision layers off) | yes | padding-free packing on by default, but only text benchmarks are published | wraps TRL's GRPOTrainer; vLLM path cannot put LoRA on vision layers | reject (coverage, beta, patching risk) |
| LLaMA-Factory v0.9.5 (2026-05-30) | not listed | `lfm2_vl` template (listed sizes 1.2B/1.6B) | `qwen3_vl` | `neat_packing` | no GRPO (points to EasyR1) | reject (a second stack, no SmolVLM) |
| ms-swift v4.6.0 (2026-10-13) | only Idefics3-8B listed | not listed | yes | `padding_free` and `packing` for SFT ("packing offers faster training and more stable memory usage") | own GRPO | reject (coverage) |
| Axolotl v0.21.0 (2026-10-16) | SmolVLM2 | LFM2-VL | Qwen3-VL | "sample_packing: false" required for multimodal | GRPO on images listed | reject (no gain over TRL, a second config system) |
| torchtune | — | — | — | — | — | "no longer actively maintained" |
| NeMo AutoModel | SmolVLM | not listed | Qwen3-VL | not verified | NeMo-RL (no SmolVLM RL) | reject (heavy; coverage) |
| oumi v0.9 (2026-09-16) | SmolVLM-Instruct (2B) | not listed | Qwen3-VL | not verified | — | reject |

Sources for the table:
- TRL: [releases](https://github.com/huggingface/trl/releases), [SFT docs](https://huggingface.co/docs/trl/main/en/sft_trainer)
- Unsloth: [vision fine-tuning](https://unsloth.ai/docs/basics/vision-fine-tuning), [packing blog](https://unsloth.ai/docs/blog/3x-faster-training-packing), [LFM2.5 guide](https://unsloth.ai/docs/models/tutorials/lfm2.5), [VLM RL](https://unsloth.ai/docs/new/vision-reinforcement-learning-vlm-rl)
- LLaMA-Factory: [README](https://github.com/hiyouga/LLaMA-Factory)
- ms-swift: [supported models](https://github.com/modelscope/ms-swift/blob/main/docs/source_en/Instruction/Supported-models-and-datasets.md), [CLI parameters](https://github.com/modelscope/ms-swift/blob/main/docs/source_en/Instruction/Command-line-parameters.md)
- Axolotl: [multimodal docs](https://docs.axolotl.ai/docs/multimodal.html)
- torchtune: [README and issue #2883](https://github.com/meta-pytorch/torchtune/issues/2883)
- NeMo AutoModel: [VLM coverage](https://docs.nvidia.com/nemo/automodel/latest/model-coverage/vlm/index.html)
- oumi: [supported models](https://oumi.ai/docs/en/latest/resources/models/supported_models.html)
- Release dates come from the GitHub releases API.

**Why staying is right.** The one feature that would justify a switch is VLM packing (ms-swift) or padding-free
(Unsloth). §1.2 shows that a length-grouped sampler recovers almost all of that padding inside TRL: waste falls to
≤ 0.2 %. The remaining gain of packing is larger effective micro-batches. For models of 256M–450M, that is better
bought by raising `per_device_train_batch_size`; memory is not the constraint (ROADMAP §GPUs: SmolVLM SFT about
6 GB).

A switch would also cost:
- a second chat-template and loss-mask implementation, which risks a mismatch between SFT and GRPO tokenisation;
- re-verifying the B2/B3 break-it guarantees (`sft_rows` image hashes, stop tokens);
- the loss of TRL's `entropy` and token-accuracy logging.

## 3. Kernels and performance knobs (single GPU)

| Knob | Recommendation | Evidence |
|---|---|---|
| Attention | SDPA (default) on all cards. FlashAttention 2 through Hub kernels (`attn_implementation="kernels-community/flash-attn2@v2"`, pinned) only if profiling shows attention matters, which is plausible for Qwen only. Hub kernels match source-built FA, "ready to use in about 2.5 seconds". Pin `@v2`: the `v3` CUDA 12.8 builds fail in the backward pass for GQA models. | [TRL kernels hub](https://huggingface.co/docs/trl/main/en/kernels_hub) |
| FA3 / FA4 | Not needed. FA3 is for Hopper. FA4 is written for Blackwell (B200, CuTe-DSL). Our sequences are ≤ 2.7k tokens, where attention is a small share of FLOPs. | [FA4 blog, 2026-03](https://tridao.me/blog/2026/flash4/) |
| Liger | Do not use. TRL already fuses the LM head. Liger's remaining layer kernels (RMSNorm, SwiGLU) support `smolvlm` and `qwen3_vl` but not `lfm2_vl`. Its headline "+20 % throughput, −60 % memory" comes from 8B-scale multi-GPU runs. | [Liger README](https://github.com/linkedin/Liger-Kernel), [monkey_patch.py](https://github.com/linkedin/Liger-Kernel/blob/main/src/liger_kernel/transformers/monkey_patch.py) |
| Gradient checkpointing | Turn it **off** for SmolVLM and LFM on a 24 GB card. Recompute costs about 2N of the 8N FLOPs per token, so this saves about 25 % if memory fits (H-A2-10 measures it). Keep it on, optionally selective, for Qwen full FT. | [TRL memory docs](https://huggingface.co/docs/trl/main/en/reducing_memory_usage) |
| Sampler | `train_sampling_strategy="group_by_length"` with a `length` column (completion tokens) added in `sft_rows`. The sampler makes megabatches of 50 × the effective batch, so batches stay random at the megabatch level. | transformers `LengthGroupedSampler`; §1.2 |
| Dataloader | `dataloader_num_workers=4–8`, `dataloader_persistent_workers=True`, `dataloader_prefetch_factor=2–4`. | §1.4 |
| Precision | bf16 on Ampere and newer (already in `configs._bf16`). | [TRL speed docs](https://huggingface.co/docs/trl/main/en/speeding_up_training) |
| Allocator | `PYTORCH_ALLOC_CONF=expandable_segments:True`. Also set `PYTORCH_CUDA_ALLOC_CONF` for vLLM, which only reads that name. | [TRL memory docs](https://huggingface.co/docs/trl/main/en/reducing_memory_usage) |
| torch.compile | No. Every batch has a different shape (variable completion lengths, and RIGA crops change LFM/Qwen image-token counts), so it would recompile. Unmeasured; low value at ≤ 1 h runs. | inference |
| FSDP2 / DeepSpeed | No. Single GPU; every model fits (ROADMAP §GPUs). | ROADMAP |
| Activation offloading | Only as a last resort for Qwen full FT on 24 GB. | [TRL memory docs](https://huggingface.co/docs/trl/main/en/reducing_memory_usage) |

## 4. Recipe

### 4.1 Full fine-tuning vs LoRA, DoRA and rsLoRA

- **"LoRA Without Regret"** ([Thinking Machines, 2025](https://thinkingmachines.ai/blog/lora/)):
  - LoRA matches full FT when it is applied to **all** layers, MLP included, and the dataset is within LoRA's
    capacity;
  - the optimal LoRA LR is about **10× full FT** (about 15× for runs of about 100 steps);
  - LoRA "is less tolerant of large batch sizes";
  - for policy-gradient RL, "LoRA fully matches … even with ranks as low as 1".
- **This matches Phase 1's settings:** full FT at 5e-5 and LoRA r16 at 5e-4 (`configs.py`, lab t_stack T4).
- **DoRA** adds a learnable magnitude; it "can improve … especially at low ranks" but has "a bigger overhead than
  pure LoRA". **rsLoRA** scales by α/√r ([PEFT LoRA guide](https://huggingface.co/docs/peft/main/en/developer_guides/lora);
  [DoRA](https://arxiv.org/abs/2402.09353); [rsLoRA](https://arxiv.org/abs/2312.03732)). Neither matters at r16
  for format learning.
- **Decision:** full FT for SmolVLM and LFM; they are small and the D10 output format is new.
- **Qwen3-VL-2B:** LoRA r16 on all linear layers at 2e-4–5e-4, merged after SFT. GRPO then starts a fresh adapter
  on the merged weights. That GRPO adapter must stay off the vision tower while vLLM colocate is used: Unsloth
  documents that "vLLM does not support LoRA for vision/encoder layers"
  ([Unsloth VLM RL](https://unsloth.ai/docs/new/vision-reinforcement-learning-vlm-rl)).

### 4.2 Freeze the vision tower, or not

The evidence favours unfreezing for localisation:
- VectorLLM fine-tunes "both the vision expert and the large language model" for point-by-point building polygons
  ([2507.04664](https://arxiv.org/abs/2507.04664)).
- Visual Position Prompt reports that unfreezing the encoder "significantly enhances" grounding
  ([2503.15426](https://arxiv.org/abs/2503.15426)).
- VLM4VLA finds that the visual module is "the primary performance bottleneck"
  ([2601.03309](https://arxiv.org/abs/2601.03309)).

Against unfreezing:
- Qwen's own fine-tuning README sets `tune_mm_vision=False` (for image + video mixtures) and vision LR 1e-6 against
  LM LR 2e-7–1e-6 ([qwen-vl-finetune](https://github.com/QwenLM/Qwen3-VL/blob/main/qwen-vl-finetune/README.md)).
- Axolotl recommends a frozen vision tower for LoRA ([Axolotl](https://docs.axolotl.ai/docs/multimodal.html)).

Our images are fundus and CT crops, far from SigLIP's pretraining, and boundary precision is a gate (≤ 2 px). Default:
**unfreeze** (SmolVLM's encoder is SigLIP-base, 93M ([SmolVLM-256M blog](https://huggingface.co/blog/smolervlm));
LFM's SigLIP2 is 86M ([card](https://huggingface.co/LiquidAI/LFM2.5-VL-450M))). The vision tower may need its own
learning rate (0.1×), set through `SFTTrainer(optimizers=…)` param groups. Tested in H-A2-6.

### 4.3 Learning rate, schedule, epochs

- **SimpleSeg** (SFT→RL polygons, the closest recipe):
  - one epoch over 800k samples, peak LR 5e-5, cosine to 2e-6, 3 % warmup, weight decay 0.1, global batch 256,
    Muon;
  - RL uses a constant 2e-6 ([SimpleSeg 2601.19228](https://arxiv.org/abs/2601.19228), Tables 4–5).
- **Schedule:** use **WSD** (`warmup_stable_decay`): 3 % warmup, constant LR, then a 10–20 % decay. The constant
  phase can be branched: decay from any saved stable checkpoint ([Wen et al., WSD river-valley, ICLR 2025](https://arxiv.org/abs/2410.05192)).
  One run then yields several SFT depths (for example 0.5, 1, 2 and 3 epochs) at little extra cost. That is exactly
  what checkpoint selection for GRPO needs (§4.10).
- **Peak LR:**
  - full FT: 5e-5 (Phase 1 and SimpleSeg), with 2e-5 and 1e-4 as the sweep;
  - LoRA: 10× those.
- **Batch and epochs:** effective batch 32–64 sequences. 1–3 epochs over all raters' targets. Prefer more unique
  data (augmentation) over repeated epochs (§4.10).

### 4.4 Image-resolution budget

The D12 settings stand:
- SmolVLM `longest_edge`: 1024 on RIGA (417 prompt tokens), 512 on LIDC (140);
- LFM and Qwen: processor defaults (≤ 256 image tokens).

The precision probe gate decides whether LIDC needs more pixels per object. For SmolVLM the next step is
`longest_edge` 1024 on LIDC as well: 5 tiles, about 320 image tokens, roughly +60 % tokens per LIDC row (§1.2). That
is A1's question; A2 only notes the cost.

### 4.5 Loss

| Option | Verdict | Evidence |
|---|---|---|
| Assistant-only NLL (current) | **keep** | §1.3: stop token trained for all three |
| **NTL-WAS** (Wasserstein-1 between the predicted digit distribution and the true digit), λ = 0.3, **place-value weighted** (our extension) | **test (H-A2-2)** | NTL: λ = 0.3 best of {0.15, 0.3, 0.8, 2.0}; "1 % overhead" on the loss; matches a regression head on rJokes (RMSE 1.81 vs 2.01 for CE); combines with Gaussian CE (0.48 vs 0.43). It treats "each digit-level token as equally important", which is the gap our place weighting fills ([Zausinger et al., ICML 2025](https://arxiv.org/abs/2411.02083); [PMLR](https://proceedings.mlr.press/v267/zausinger25a.html); [ntloss](https://pypi.org/project/ntloss/)). Feasibility: §1.5 |
| Label smoothing | reject | deprecated in TRL SFT, extra forward pass (§1.1); "standard regularization approaches (KL penalties, label smoothing) didn't resolve the collapsed checkpoint" ([2606.18487](https://arxiv.org/abs/2606.18487)) |
| DFT (`loss_type="dft"`) | reject; negative control only | reweights by token probability (mode-seeking), lowers entropy; "lacks distributional anchoring, leading to progressive drift" ([ASFT, ICLR 2026](https://arxiv.org/abs/2509.23753); [DFT](https://arxiv.org/abs/2508.05629)) |
| GEM (entropy-regularised reverse-KL SFT) | park | keeps diversity and helps best-of-n ([ICLR 2025](https://arxiv.org/abs/2408.16673)); not in TRL; our diversity source is the rater mixture itself |
| CHORD (SFT as a weighted auxiliary loss inside RL) | park (Phase 3 idea) | [2508.11408](https://arxiv.org/abs/2508.11408) |

### 4.6 Multi-rater targets

The options are:
- (a) all voting raters as rows (current `sft_rows`);
- (b) one random rater per item per epoch (ROADMAP's wording);
- (c) a consensus target.

(a) and (b) give the same expected gradient: both minimise cross-entropy to the empirical rater distribution. At
equal steps they differ only in variance and bookkeeping. (a) needs no custom sampler, and "an epoch then shows each
item's whole rater distribution" (dataset.py docstring). (c) removes the very distribution D5's reward matches. In
segmentation, label sampling gives better calibration than a single ground truth:
- Ji et al.: training with "unique ground-truths would be over-confident", while label sampling "was better
  calibrated" ([CVPR 2021](https://openaccess.thecvf.com/content/CVPR2021/papers/Ji_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-Rater_Agreement_Modeling_CVPR_2021_paper.pdf)).
- Lemay et al.: STAPLE with conventional training is overconfident, and the best fusion is task-dependent
  ([2202.07550](https://arxiv.org/abs/2202.07550)).

**Keep (a).** H-A2-4 checks that it does not lose greedy IoU against a medoid-rater target. It should also check
that T = 1 samples do not splice two raters into one illegal ring.

### 4.7 Curriculum and synthetic pretraining

- **Easy-first subsets hurt RL:** training "on shortest examples boosted SFT performance faster … yet final RL
  performance was significantly lower" ([Quagmires, ICLR 2026](https://arxiv.org/abs/2510.01624)).
- **Broad weak-label pretraining before SFT helps:** in SimpleSeg it adds +4.6 gIoU after SFT and +13.0 after RL,
  and RL adds +9.7–10.5 over SFT alone (Table 3, [2601.19228](https://arxiv.org/abs/2601.19228)).
- **Synthetic geometry data teaches localisation:** a procedural engine of 200k synthetic shapes lifts
  Florence-2/Qwen-VL from zero-shot to 49 % IoU ([2604.02893](https://arxiv.org/abs/2604.02893)).
- **Our data has a rare-class problem:** 317 / 23,469 LIDC train targets have holes and 200 have several parts
  (§1.5 count). RIGA has neither.

Recommendation:
- no difficulty curriculum;
- a **targeted** fix for A0′ syntax: oversample hole and part targets, or add a stage-0 of synthetic random
  polygons with holes. These are certified by polycore, rendered as CT-like or fundus-like crops, and kept in the
  same D10 JSON. Tested in H-A2-5.

### 4.8 Mixing RIGA and LIDC

The prompts differ per dataset, so a single joint model is task-conditioned. At natural proportions:
- RIGA is 12 % of rows but about 25–28 % of tokens (SmolVLM 2.4M of 8.8M), because its targets are about 3× longer;
- TRL normalises the loss per token across the accumulated batch (`num_items_in_batch`), so RIGA's gradient share
  follows the token share, not the row share.

Default: one joint run at natural proportions. Gates are checked per dataset (`phase2_entry`). If RIGA misses its
gates while LIDC passes, oversample RIGA 2× rather than train separate models. TRL's `DatasetConfig.fraction` mixes
by share ([script_utils](https://huggingface.co/docs/trl/main/en/script_utils)), but our `sft_rows` concatenation with
row repeats is simpler and keeps the image-hash checks.

### 4.9 Avoiding format collapse (truncation, loops, never stopping)

- **Already in place:**
  - canonical start vertex and orientation (D10). Fixed starting rules give "slight performance improvements and
    stability" (VectorLLM), and without clockwise order SimpleSeg's sequences "fail to form a valid polygon";
  - stop token in the loss (§1.3);
  - over-budget targets dropped (I5).
- **Moderate vertex density:** SimpleSeg's performance is unimodal in token length (best at about 221 tokens,
  worse at 78 and 859). Our medians of 71–317 are in range.
- **Watch these in SFT logs:**
  - **greedy parse, legal and truncation rate on val** every N steps;
  - **TRL `entropy`**, which falls as SFT sharpens.

### 4.10 What makes SFT a good GRPO start

- **SFT is needed for format:** "SFT stabilizes the model's output format, enabling subsequent RL to achieve its
  performance gains" ([Chu et al. 2025](https://arxiv.org/abs/2501.17161)). Zero-shot, every model fails D10
  (t_stack T2), and SimpleSeg's no-SFT pretrain-only model scores 25.3 gIoU.
- **Too much SFT hurts RL:**
  - Quagmires: 4 epochs gave +8.75 % SFT accuracy but −1.43 % post-RL gains. Pass@64 and held-out loss are better
    proxies than SFT accuracy. One epoch on unique data beats more epochs on half the data
    ([2510.01624](https://arxiv.org/abs/2510.01624)).
  - Rank inversion: pre-RL pass@1 rises with SFT depth while peak GRPO pass@10 falls from 0.806 to 0.481, and
    pre-RL entropy correlates with the GRPO outcome (ρ = +0.69) ([2606.18487](https://arxiv.org/abs/2606.18487)).
  - Plasticity loss after excessive SFT ([2606.09932](https://arxiv.org/abs/2606.09932)).
  - For VLMs, SFT can lock models "into imitative, rigid … modes" ([2504.11468](https://arxiv.org/abs/2504.11468)).
- **Repairs:**
  - WiSE-FT (interpolate an early and a late SFT checkpoint) "almost completely recovers Pass@k while also
    improving Pass@1", and is a better RL start ([2504.10478](https://arxiv.org/abs/2504.10478));
  - TailSFT filters already-fitted sequences: +17 points pass@16 and +4 points GRPO pass@1 ([2608.25756](https://arxiv.org/abs/2608.25756)).
- **Why entropy matters:** RL trades entropy for reward (R = −a·e^H + b), so the starting entropy bounds headroom
  ([Cui et al. 2025](https://arxiv.org/abs/2505.22617)). On-policy RL also stays KL-close to its start
  ([RL's Razor](https://arxiv.org/abs/2509.04259)), so a sharp SFT start constrains GRPO.

**Checkpoint selection rule (proposal):**
- Among WSD branch checkpoints that pass the format gates (parse ≥ 98 %, legal ≥ 90 %) at **T = 1**, take the
  one with the highest GRPO-signal score: the share of groups with σ above the floor, then duplicates ≤ 10 %, then
  pass@8 (IoU_tc ≥ 0.5) and the lowest held-out NLL.
- Do not pick by greedy legal rate or IoU.
- WiSE-FT interpolation is the cheap repair when the latest branch is too sharp.

## 5. Reproducibility, checkpoints, tracking, configuration

- **Seeds:** set `seed` and `data_seed` explicitly; the samplers are seeded through them.
  - Do not use `full_determinism=True`: it forces `CUDA_LAUNCH_BLOCKING=1` (§1.1), which serialises kernel
    launches.
  - Optional deterministic mode: export `CUBLAS_WORKSPACE_CONFIG=:4096:8` *before* the process starts, then call
    `torch.use_deterministic_algorithms(True, warn_only=True)`. Its cost is measured once on Vast (H-A2-9).
  - PyTorch: "Completely reproducible results are not guaranteed across PyTorch releases, individual commits, or
    different platforms", and "deterministic operations are often slower" ([PyTorch randomness notes](https://docs.pytorch.org/docs/2.14/notes/randomness.html)).
  - Practical target: identical data order and run-to-run loss curves within noise, with versions and kernel
    revisions pinned (`flash-attn2@v2`), not bitwise identity.
- **Evaluation sampling:** generation results depend on the batch composition. The main cause of inference
  nondeterminism is "lack of batch invariance"; batch-invariant kernels cost about 1.6–2.1× ([Thinking Machines](https://thinkingmachines.ai/blog/defeating-nondeterminism-in-llm-inference/)).
  Evaluate with fixed batch composition and fixed seeds, and report a confidence interval over K samples rather than
  expecting bitwise repeats.
- **Checkpoints:**
  - `save_strategy="steps"` about every 10 min of wall time, `save_total_limit=3` for rolling checkpoints, with the
    optimizer state so an interrupted instance can resume;
  - plus the named WSD branch checkpoints, kept;
  - model-only copies (`save_only_model=True` in a separate final save) are what reach GRPO;
  - storage and transfer belong to A3.
- **Tracking:** trackio, `report_to="trackio"`. It is free, local-first (SQLite), has a W&B-compatible API, can sync
  to an HF Space or Bucket, and TRL has `trackio_space_id` / `trackio_bucket_id` fields
  ([trackio docs](https://huggingface.co/docs/trackio/index); SFTConfig signature).
  - Fallback: W&B if the owner already uses it.
  - Always also keep TRL's JSON `trainer_state.json` and our gate reports in the run directory, so the evidence
    does not depend on a SaaS.
  - TensorBoard and MLflow add nothing here.
- **Configuration:** keep `objdel.train.configs` as the single source of defaults. Experiment overrides go in YAML
  read by `TrlParser.parse_args_and_config` (`--config run.yaml`; "Command line arguments can override values set
  by the config file"; an `env:` block sets environment variables) ([TRL script_utils](https://huggingface.co/docs/trl/main/en/script_utils)).
  Write the resolved config, git SHA, prepared-manifest hash and `pip freeze` into each run directory. No Hydra.

## 6. Patterns and trade-offs

1. **The waste is mostly prompt tokens and CPU, not padding.**
   - 47–79 % of tokens are prompt.
   - Prompt tokens per row are constant within a dataset for SmolVLM and LIDC, so the image encoder does the same
     work again on every rater row of an item. RIGA has 6 rows per image and LIDC about 3.2.
   - A collator that groups an item's rater rows into one sequence would share the image. That is impossible
     without packing or a multi-answer format, so the cost is accepted.
   - Augmentation (flips, 90° rotations, re-targeted per dataset.py) turns the repeated images into unique ones.
     That matches "diverse unique > repeated epochs" ([Quagmires](https://arxiv.org/abs/2510.01624)).
2. **The SFT objective and the GRPO objective pull in different directions.**
   - SFT on a rater mixture is already a distribution-matching objective (forward KL to the raters), which is what
     D5 rewards.
   - Anything that sharpens it (DFT, consensus targets, many epochs) helps greedy IoU and hurts GRPO signal.
   - So the SFT knobs that matter most are *depth* and *selection*, not loss tweaks. NTL is the one tweak that adds
     geometry knowledge rather than sharpness.
3. **Cutting-edge yet safe:**
   - WSD branching combined with GRPO-signal-based checkpoint selection, and WiSE-FT as a fallback. Both are pure
     orchestration over TRL, with no new kernels.
   - The place-weighted NTL is the only new code in the loss, guarded by an A/B test.
4. **Risk concentration:** the riskiest component is not the trainer but the CPU path (collation, augmentation
   re-targeting). Measure GPU utilisation in the first Vast session (`nvidia-smi dmon`); below 70 % means the
   CPU-side fix comes first.

## 7. Ranked recommendations

1. **Stack:** TRL 1.16 `SFTTrainer`, `objdel.train.configs.sft_config`. Unchanged: assistant-only loss,
   `max_length=None`, no packing, bf16. Add:
   - `train_sampling_strategy="group_by_length"` with a `length` column;
   - `dataloader_num_workers=6`, `dataloader_persistent_workers=True`;
   - `gradient_checkpointing=False` for SmolVLM and LFM.
2. **SmolVLM-256M and LFM2.5-VL-450M:**
   - full FT, vision tower unfrozen (H-A2-6), AdamW fused, weight decay 0.0–0.1;
   - peak LR 5e-5 (sweep 2e-5, 1e-4), WSD with 3 % warmup and a 15 % linear decay to 0;
   - effective batch 32, 3 epochs of stable LR, with decay branches after 0.5, 1, 2 and 3 epochs.
3. **Qwen3-VL-2B:**
   - LoRA r16 on all linear layers (vision included for SFT), LR 2e-4–5e-4, same WSD, merged after SFT;
   - on an H100 or A100 for a 1-hour run, or subsample LIDC to about 30 % on a 4090.
4. **Targets:** all voting raters as rows (current), joint RIGA + LIDC at natural proportions. Hole and part targets
   oversampled ×4–8 for A0′ (H-A2-5).
5. **Checkpoint selection for GRPO:** the §4.10 rule; WiSE-FT α = 0.5 between the 0.5-epoch and the selected branch
   as a fallback (H-A2-7).
6. **Optional loss:** CE + 0.3 · place-weighted NTL-WAS, if H-A2-2 passes.
7. **Operations:**
   - trackio plus run-dir JSON;
   - TrlParser YAML overrides;
   - seed and data_seed set;
   - no `full_determinism`;
   - checkpoints every about 10 min with optimizer state.

## 8. Rejected alternatives

| Alternative | Reason |
|---|---|
| Switching to ms-swift, LLaMA-Factory, Axolotl or Unsloth for SFT | Model coverage gaps (§2), a second template and mask path, and the packing gain is recovered by length grouping (§1.2) |
| torchtune | Not maintained ([README](https://github.com/meta-pytorch/torchtune)) |
| Packing or padding-free for VLMs in TRL | Refused by TRL 1.16 for vision datasets; ≤ 0.2 % padding remains after grouping |
| Liger kernels | Deprecated path in TRL SFT; no `lfm2_vl`; fused cross-entropy duplicates TRL's fused head |
| FA3 / FA4 | Wrong hardware or no benefit at our sequence lengths |
| `full_determinism=True` | `CUDA_LAUNCH_BLOCKING=1` serialises the GPU |
| Label smoothing | Deprecated in TRL SFT (extra forward pass); does not fix over-trained checkpoints |
| DFT | Mode-seeking; drift and entropy collapse (ASFT); opposite of what GRPO needs |
| Consensus or STAPLE targets | Overconfident; removes the rater distribution that D5 matches |
| Easy-first curriculum (short or simple shapes first) | Faster SFT gains, worse RL ([Quagmires](https://arxiv.org/abs/2510.01624)) |
| Selecting the SFT checkpoint by greedy legal or IoU | Rank inversion under RLVR ([2606.18487](https://arxiv.org/abs/2606.18487)) |
| Qwen recipe LR 2e-7–1e-6 with frozen vision | Tuned for general mixed image and video data on large models, not for teaching a new output format; LoRA at 10× full-FT LR is the evidence-backed choice |
| DoRA | Overhead with no gain expected at r16 for format learning |
| torch.compile | Every batch has a different shape (recompiles); untested; low value for ≤ 1 h runs |

## 9. Hypotheses to test (wave B, Mac)

Common setup (unless stated otherwise):
- SmolVLM-256M on LIDC only (short sequences: 271 tokens mean; RIGA needs `longest_edge` 1024 and more memory).
- Full FT with micro-batch 2 and gradient accumulation; LoRA r16 if 7 GB is exceeded.
- MPS memory capped with `torch.mps.set_per_process_memory_fraction` as in `smoke_grpo`.
- Validation: a fixed list of 100 LIDC val items, including all 44 disagreement items and every item with a hole
  rater.
- Generation with `objdel.train.sample` (HF generate). Scoring with `objdel.eval.metrics` / `summarize` for parse,
  legal, truncation, IoU_tc, topology and duplicates, plus the GRPO-signal gate definitions as in K12.
- 2 seeds per arm where the pass criterion is a difference.

H-A2-1 to H-A2-8 need **the MPS slot** (≤ 7 GB, one at a time). H-A2-7 trains nothing (it only interpolates and
evaluates H-A2-3's checkpoints), but it still uses MPS for sampling. H-A2-9 and H-A2-10 run on the first Vast box.
Suggested order, cheapest and most decisive first: H-A2-3 (with H-A2-8 and H-A2-7 riding on it), then H-A2-5,
H-A2-4, H-A2-2, H-A2-6, H-A2-1.

| ID | Hypothesis | Test design | Pass criterion |
|---|---|---|---|
| **H-A2-1** | Length-grouped sampling speeds up SFT without changing learning. | 2 arms (random vs `group_by_length` with a `length` column), 150 steps, batch 8 (micro 2 × 4), same seed. Log wall time per step and the train/val NLL curve. | Grouped is ≥ 15 % faster in wall time per epoch-equivalent of tokens, **and** final val NLL is within 2 % of random. |
| **H-A2-2** | Place-weighted NTL-WAS (λ = 0.3) improves coordinate accuracy without hurting format or diversity. | CE vs CE + 0.3·NTL (the `ntl_proto.py` loss as an `SFTTrainer.compute_loss` override; the template from `get_training_chat_template`), 400 steps, 2,000 LIDC train items with all raters, 2 seeds. Evaluate greedy plus K = 8 at T = 1 on the 100 val items. | Median IoU_tc +0.02 or more, **or** mean boundary distance −10 % or more (greedy); parse and legal not lower by > 1 pt; pass@8 and σ-above-floor share not lower by > 2 pt. |
| **H-A2-3** | Deeper SFT raises greedy quality but lowers the GRPO signal (rank inversion appears at our scale). | One WSD run (peak 5e-5, stable LR); branch decays (40 steps) at 100, 200, 400 and 800 steps. For each branch: greedy legal and IoU; T = 1, G = 8 duplicates; σ-above-floor share; pass@8; token entropy; val NLL. Then 20 GRPO steps (K12 harness, LoRA) from the best-greedy and the best-signal branch. | Prediction confirmed if the best-signal branch ≠ the best-greedy branch **and** its 20-step mean aggregate reward gain ≥ the other's. If both branches coincide, the §4.10 rule costs nothing and stays. |
| **H-A2-4** | All-rater targets keep the rater distribution at no greedy cost compared with a single medoid-rater target. | Arm A: current `sft_rows` (all voting raters). Arm B: one row per item, the rater with the highest mean IoU to the others. Equal steps (300). Evaluate topology TV/JS against the rater histogram on the 44 disagreement items (K = 8, T = 1), legal, greedy IoU_tc. | TV(A) ≤ TV(B) − 0.05; legal(A) ≥ legal(B) − 1 pt; greedy IoU_tc(A) ≥ IoU_tc(B) − 0.01; spliced-rater illegal rings (illegal at T = 1) ≤ 2 %. |
| **H-A2-5** | A0′ hole syntax is under-learned at 1.35 % prevalence; oversampling fixes it. | Arms: hole and part targets ×1 vs ×8 (row repeats in `sft_rows`), 400 steps. Evaluate on val items with a hole or part rater (47 hole targets) and on hole-free items, K = 8 at T = 1. | ×8 at least doubles the rate of samples whose topology matches a hole rater on hole items; hole hallucination on hole-free items rises ≤ 1 pt; overall legal not lower by > 1 pt. |
| **H-A2-6** | Unfreezing the vision tower improves boundary precision. | Frozen vision vs unfrozen (same LR) vs unfrozen at 0.1× LR (optimizer param groups), 400 steps. Evaluate greedy IoU_tc and boundary error on the smallest-third LIDC val nodules, plus the synthetic precision probe if A4 provides it. | The best unfrozen arm beats frozen by ≥ 0.02 median IoU_tc or ≥ 10 % boundary error, with legal not lower; otherwise freeze (cheaper, less drift). |
| **H-A2-7** | WiSE-FT restores diversity of an over-trained checkpoint. | No new training: interpolate the H-A2-3 checkpoints (100-step and 800-step branches) at α ∈ {0.25, 0.5, 0.75} and evaluate as in H-A2-3. | Some α has legal (T = 1) ≥ late − 1 pt **and** σ-above-floor share ≥ early's. |
| **H-A2-8** (MPS, an extra arm of H-A2-3) | DFT sharpens: lower entropy, fewer distinct samples. | `loss_type="dft"`, 200 steps, from the same init and data order as H-A2-3's 200-step branch; measure TRL `entropy`, duplicates at T = 1 and σ share against that branch. | DFT's duplicate rate ≥ NLL's + 5 pt, or its σ share ≤ NLL's − 10 pt. A pass confirms the rejection in §8. |
| **H-A2-9** (first Vast box) | Deterministic mode is affordable. | SmolVLM SFT, 100 steps, twice with `CUBLAS_WORKSPACE_CONFIG=:4096:8` and `use_deterministic_algorithms(True)`, twice without. | Repeat runs bitwise-equal in loss with determinism on; overhead ≤ 15 % → enable for the decisive runs, else seeds only. |
| **H-A2-10** (first Vast box) | Turning gradient checkpointing off and using 6 workers saves ≥ 25 % SFT wall time on SmolVLM. | 2 × 2 grid (checkpointing on/off × workers 0/6), 100 steps each; log step time, GPU utilisation, peak memory. | Off + 6 workers ≥ 25 % faster than the TRL default with peak memory ≤ 20 GB; GPU utilisation ≥ 70 %. |

## 10. Gaps and unverified items

- **The web-search budget was exhausted mid-task.** Not searched:
  - published uses of NTL or Gaussian coordinate losses in VLM grounding papers;
  - measured SFT throughput of sub-1B VLMs across frameworks;
  - Hägele et al. (WSD for fixed budgets, [2405.18392](https://arxiv.org/abs/2405.18392)) and "LoRA Learns Less and
    Forgets Less" ([2405.09673](https://arxiv.org/abs/2405.09673)) are cited from memory and only for context;
  - RTX 4090 and H100 peak TFLOPS are datasheet values not fetched here
    ([Ada whitepaper](https://images.nvidia.com/aem-dam/Solutions/geforce/ada/nvidia-ada-gpu-architecture.pdf),
    [H100](https://www.nvidia.com/en-us/data-center/h100/)).
- **Rex-Omni** (SFT then GRPO on 0–999 coordinate tokens, [2510.12798](https://arxiv.org/abs/2510.12798)) was seen
  only at abstract level.
- **ms-swift `padding_free` for Qwen3-VL** was not verified beyond the CLI documentation's general SFT statement.
- **Freeze-vs-unfreeze evidence** for VLM4VLA on small models comes partly from a search snippet; the abstract
  confirms only "the visual module … is the primary performance bottleneck".
- **MFU of 25 %** for 256M–2B models is an assumption, and small models are often below it. The first Vast session's
  per-step time replaces §1.2's estimates.
- **The collation timings** (§1.4) were taken under a load average of about 31 on the Mac, so they are upper bounds.

## Sources

Framework and library:
- TRL: [SFT](https://huggingface.co/docs/trl/main/en/sft_trainer), [memory](https://huggingface.co/docs/trl/main/en/reducing_memory_usage), [speed](https://huggingface.co/docs/trl/main/en/speeding_up_training), [kernels hub](https://huggingface.co/docs/trl/main/en/kernels_hub), [script utils](https://huggingface.co/docs/trl/main/en/script_utils), [v1.16.0 source](https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/sft_trainer.py), [chat_template_utils](https://github.com/huggingface/trl/blob/v1.16.0/trl/chat_template_utils.py), [releases](https://github.com/huggingface/trl/releases)
- transformers v5.20.1: [training_args](https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/training_args.py), [trainer_utils](https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/trainer_utils.py), [optimization](https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/optimization.py)
- Unsloth: [vision](https://unsloth.ai/docs/basics/vision-fine-tuning), [packing](https://unsloth.ai/docs/blog/3x-faster-training-packing), [LFM2.5](https://unsloth.ai/docs/models/tutorials/lfm2.5), [VLM RL](https://unsloth.ai/docs/new/vision-reinforcement-learning-vlm-rl)
- Other frameworks: [LLaMA-Factory](https://github.com/hiyouga/LLaMA-Factory), [ms-swift models](https://github.com/modelscope/ms-swift/blob/main/docs/source_en/Instruction/Supported-models-and-datasets.md), [ms-swift CLI](https://github.com/modelscope/ms-swift/blob/main/docs/source_en/Instruction/Command-line-parameters.md), [Axolotl multimodal](https://docs.axolotl.ai/docs/multimodal.html), [torchtune](https://github.com/meta-pytorch/torchtune/issues/2883), [NeMo AutoModel VLMs](https://docs.nvidia.com/nemo/automodel/latest/model-coverage/vlm/index.html), [oumi](https://oumi.ai/docs/en/latest/resources/models/supported_models.html)
- Kernels: [Liger](https://github.com/linkedin/Liger-Kernel), [Liger monkey_patch](https://github.com/linkedin/Liger-Kernel/blob/main/src/liger_kernel/transformers/monkey_patch.py), [FlashAttention-4](https://tridao.me/blog/2026/flash4/)
- PEFT and tracking: [PEFT LoRA guide](https://huggingface.co/docs/peft/main/en/developer_guides/lora), [trackio](https://huggingface.co/docs/trackio/index)
- Models: [Qwen3-VL fine-tune README](https://github.com/QwenLM/Qwen3-VL/blob/main/qwen-vl-finetune/README.md), [LFM2.5-VL-450M](https://huggingface.co/LiquidAI/LFM2.5-VL-450M), [SmolVLM-256M blog](https://huggingface.co/blog/smolervlm), [SmolVLM paper](https://arxiv.org/abs/2504.05299)
- Reproducibility: [PyTorch reproducibility](https://docs.pytorch.org/docs/2.14/notes/randomness.html), [nondeterminism in inference](https://thinkingmachines.ai/blog/defeating-nondeterminism-in-llm-inference/)

Recipe papers:
- Polygon and coordinate VLMs: [SimpleSeg](https://arxiv.org/abs/2601.19228), [VectorLLM](https://arxiv.org/abs/2507.04664), [Rex-Omni](https://arxiv.org/abs/2510.12798)
- Fine-tuning methods: [LoRA Without Regret](https://thinkingmachines.ai/blog/lora/), [DoRA](https://arxiv.org/abs/2402.09353), [rsLoRA](https://arxiv.org/abs/2312.03732)
- Losses: [NTL / Regress, Don't Guess](https://arxiv.org/abs/2411.02083), [DFT](https://arxiv.org/abs/2508.05629), [ASFT](https://arxiv.org/abs/2509.23753), [GEM](https://arxiv.org/abs/2408.16673), [CHORD](https://arxiv.org/abs/2508.11408)
- SFT→RL interplay: [SFT Memorizes, RL Generalizes](https://arxiv.org/abs/2501.17161), [SFT or RL for LVLMs](https://arxiv.org/abs/2504.11468), [Quagmires](https://arxiv.org/abs/2510.01624), [WiSE-FT for reasoning](https://arxiv.org/abs/2504.10478), [SFT overtraining rank inversion](https://arxiv.org/abs/2606.18487), [Rejuvenating plasticity](https://arxiv.org/abs/2606.09932), [TailSFT](https://arxiv.org/abs/2608.25756), [RL's Razor](https://arxiv.org/abs/2509.04259), [Entropy mechanism](https://arxiv.org/abs/2505.22617)
- Schedules: [WSD](https://arxiv.org/abs/2410.05192)
- Multi-rater training: [Ji et al. CVPR 2021](https://openaccess.thecvf.com/content/CVPR2021/papers/Ji_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-Rater_Agreement_Modeling_CVPR_2021_paper.pdf), [Lemay et al.](https://arxiv.org/abs/2202.07550)
- Vision tower: [VLM4VLA](https://arxiv.org/abs/2601.03309), [Visual Position Prompt](https://arxiv.org/abs/2503.15426)
- Synthetic data: [synthetic geometry grounding](https://arxiv.org/abs/2604.02893)
