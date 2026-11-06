# Fast transformer fine-tuning on an RTX 4090 (Ada Lovelace): hardware ceilings and training-loop levers

Research date 2026-11-06. Researched ~20 tool calls. Several sources I wanted (HF transformers training-perf page, Unsloth blog, torch.compile on 4090, FA2 on Qwen3-VL specifics) were not retrievable; see Gaps.

## 1. RTX 4090 (AD102) hardware specs and peak Tensor Core throughput

### Takeaway
The 4090 has 128 SMs, 512 4th-gen Tensor Cores, 24 GB GDDR6X at 1008 GB/s and 72 MB L2. Dense bf16 Tensor peak (FP32 accumulate) is 165.2 TFLOPS (330.4 with sparsity); fp16 with FP16 accumulate is 330.3 dense, but bf16 has no FP16-accumulate mode. That 165 TFLOPS is the correct denominator for MFU in bf16 training.

### Cited Findings
- 128 SMs, 16384 CUDA cores, 4 Tensor Cores/SM = 512 4th-gen Tensor Cores, boost clock 2520 MHz, TSMC 4N, 450 W TGP, 76.3 B transistors — [NVIDIA Ada GPU Architecture whitepaper v2.1, Appendix A](https://images.nvidia.com/aem-dam/Solutions/Data-Center/l4/nvidia-ada-gpu-architecture-whitepaper-v2.1.pdf)
- Peak BF16 Tensor (FP32 accumulate): 165.2 dense / 330.4 with sparsity TFLOPS. Peak FP16 Tensor with FP32 accumulate: 165.2 / 330.4. Peak FP16 Tensor with FP16 accumulate: 330.3 / 660.6. TF32 Tensor: 82.6 / 165.2. INT8: 660.6 / 1321.2 TOPS — same whitepaper, Appendix A (peak rates are at boost clock; the second number is "effective" with the Sparsity feature)
- FP8 Tensor: 660.6 / 1321.2 TFLOPS with FP16 accumulate; 330.3 / 660.6 with FP32 accumulate (the FP32-accumulate figure was corrected in whitepaper v2.0.2) — same whitepaper, Appendix A footnote 4
- Non-Tensor FP32/FP16/BF16 (CUDA core) peak is 82.6 TFLOPS — same whitepaper
- Memory: 24 GB GDDR6X, 384-bit, 21 Gbps, 1008 GB/s; L2 cache 73728 KB (72 MB), vs 6144 KB on RTX 3090 Ti; L1/shared 16384 KB total; register file 32768 KB — same whitepaper
- 4th-gen Tensor Core doubles throughput for FP16, BF16, TF32 and INT8 versus the 3090 Ti and adds FP8 (halves storage, doubles throughput vs FP16) — same whitepaper, Ada architecture section
- Whitepaper body text also mentions a 96 MB L2 for the full AD102 chip; the shipped 4090 is 72 MB per Appendix A — same whitepaper (full die vs shipped SKU; Appendix A is the 4090 figure)
- Transformer Engine lists "FP8 training on NVIDIA Hopper, Ada, Blackwell" GPUs; says nothing about a minimum compute capability on that page — [NVIDIA TransformerEngine GitHub](https://github.com/NVIDIA/TransformerEngine)
- Compare: RTX 3090 Ti dense BF16 80 TFLOPS, same bandwidth 1008 GB/s — whitepaper Appendix A (so Ada roughly doubled compute at equal bandwidth, raising the arithmetic-intensity needed to be compute-bound)
- PyTorch TF32: `torch.backends.cuda.matmul.allow_tf32` defaults to False since 1.12; `torch.backends.cudnn.allow_tf32` defaults to True; TF32 Tensor Cores exist since Ampere; from PyTorch 2.9 the docs recommend the newer `fp32_precision` settings and do not allow mixing old and new — [PyTorch CUDA semantics notes](https://docs.pytorch.org/docs/2.15/notes/cuda.html)
- `set_float32_matmul_precision` default is `'highest'`; `'high'`/`'medium'` trade precision for speed — [PyTorch tuning guide](https://docs.pytorch.org/tutorials/recipes/recipes/tuning_guide.html)

### Inferences
- Ridge point for bf16 matmul: 165.2e12 / 1.008e12 ≈ 164 FLOP/byte. Matmuls with small M (tokens per micro-batch) times the hidden size can fall below this; with 2240 tokens per micro-batch (4 x 560) and hidden ~2048 the big linears are compute-bound, but elementwise ops (norms, RoPE, SwiGLU, dtype casts, LoRA's thin r=16 matmuls) are bandwidth-bound.
- TF32 only matters for fp32 matmuls; a pure bf16 run does not use it. Peak TF32 is half of bf16 (82.6 vs 165.2).
- FP8 on 4090 would, per whitepaper, double bf16 tensor peak (330 FP32-accumulate); but for LoRA on a 2B model with small GEMMs, casting overhead likely eats much of that (not measured; see Gaps).
- The "330 TFLOPS" fp16 headline requires FP16 accumulation, which PyTorch bf16 autocast does not use; do not use it as an MFU denominator.

### Gaps
- No NVIDIA-sourced measured sustained cuBLAS bf16 GEMM TFLOPS on a 4090 (power/clock throttling at 450 W); a realistic achievable fraction of the 165.2 peak was not found.
- Whether Transformer Engine FP8 gives an end-to-end speedup on a 4090 for a ~2B model: no benchmark found.

## 2. MFU estimation for LoRA fine-tuning

### Takeaway
Use MFU = (model FLOPs per token x tokens/s) / 165.2e12. For LoRA the frozen base needs no weight-gradient matmuls, so FLOPs per token are ~4N (forward 2N + backward activation-gradient 2N), or ~6N if gradient checkpointing re-runs the forward, not the 6N/8N of full fine-tuning.

### Cited Findings
- FlashAttention-2 reports full GPT training at up to 225 TFLOPs/s per A100, equal to 72% MFU (MFU as fraction of 312 TFLOPS peak) — [FlashAttention repo](https://github.com/Dao-AILab/flash-attention)
- Megatron selective recomputation: 530B model on 2240 A100s reached 54.2% MFU, 29% faster than 42.1% with full recomputation; full recomputation costs 30-40% execution-time overhead in their runs — [Korthikanti et al. 2022, arXiv 2205.05198](https://arxiv.org/abs/2205.05198v1) (via search summary of the paper's abstract/text; the "33%" figure is a simple derived estimate, not measured in the paper)
- LoRA paper claim: vs full fine-tuning of GPT-3 175B with Adam, 10,000x fewer trainable params, 3x less GPU memory, "higher training throughput" — [PEFT LoRA guide](https://huggingface.co/docs/peft/main/en/developer_guides/lora)

### Inferences (worked example for the user's run; assumptions marked)
- Tokens per optimizer step = 4 x 8 x ~560 = ~17,920; at 3.4 s/step that is ~5,270 tokens/s.
- Assume LM params N ~ 1.7e9 (assumption: Qwen3-VL-2B LM; verify from the config) plus frozen vision tower forward on the image patches (not counted here).
- With checkpointing: ~6N FLOPs/token = 1.02e10 x 5,270 ≈ 54 TFLOPS ≈ 33% MFU of 165.2 (hardware FLOPs). "Useful" model FLOPs (4N, excluding recompute) ≈ 36 TFLOPS ≈ 22% MFU. Attention FLOPs and LoRA FLOPs are small at 560 tokens and are ignored.
- Report both "model-FLOPs MFU" (4N) and "hardware-FLOPs utilization" (6N with recompute) so comparisons stay honest.
- For a 2B model at ~560-token sequences, 22-33% is plausible; the LM is small enough that elementwise/launch overhead is a noticeable share.

### Gaps
- No published MFU numbers for LoRA fine-tuning of a ~2B VLM on a 4090 were found; the worked example above is my arithmetic, not a measurement.

## 3. Gradient checkpointing: cost, selective checkpointing, when to turn it off

### Takeaway
Full checkpointing re-runs the forward in backward: ~30-40% measured overhead in Megatron for full-FT (theoretical 4/3). For LoRA with frozen base the relative overhead is larger (~6N vs 4N, ~50% more matmul FLOPs, my derivation), so with only 6.9 GB of 24 GB used, turning it off is the first lever to test.

### Cited Findings
- Full recomputation: 30-40% execution time overhead; selective recomputation (attention core only) reduces activation memory ~5x while cutting recompute overhead by over 90% — [arXiv 2205.05198](https://arxiv.org/abs/2205.05198v1)
- NVIDIA Megatron Bridge docs make attention-only selective recompute (`recompute_granularity="selective"`, `recompute_modules=["core_attn"]`) the usual starting point, and warn that recomputing MLP plus attention is worse than MLP alone due to double recompute — [NVIDIA Megatron Bridge, Activation Recompute](https://docs.nvidia.com/nemo/megatron-bridge/nightly/skills/perf-activation-recompute/SKILL.html)
- TRL: gradient checkpointing is on by default in all trainers (`gradient_checkpointing=False` disables it); SFTTrainer offers `gradient_checkpointing_kwargs={"selective": True}` saving the attention output so it is not recomputed, "recovers most of the checkpointing slowdown at long context", costs one hidden-state-sized tensor per layer, forces non-reentrant checkpointing — [TRL, Reducing Memory Usage](https://huggingface.co/docs/trl/main/en/reducing_memory_usage)
- TRL also offers `activation_offloading=True` (CPU offload via saved_tensors_hooks, "slightly increased training time") — same TRL page

### Inferences
- Frozen-base LoRA FLOPs: forward 2N, activation-grad backward 2N (weight-grad matmuls skipped) => 4N; with checkpointing +2N => 6N. So removing checkpointing could cut ~33% of step time on the matmul-bound share (upper bound; step also has vision tower, loss, optimizer, data time). Activation memory per micro-batch of 2240 tokens without checkpointing for a 28-layer, hidden-2048 model is likely a few GB to ~10+ GB; with 17 GB headroom it probably fits, but this must be measured (test mb=4 no-GC; if OOM, mb=2 x accum 16 no-GC, or selective).
- "Turn it off when peak memory is far below capacity" is the rule that follows; the user's 6.9 GB peak is the strongest evidence in this note-set.

### Gaps
- The HF transformers training-performance page (which states checkpointing slowdown, ~20-25%) could not be fetched (URL served the inference overview); no HF-primary number for the slowdown. No measured activation-memory-vs-checkpointing figure for Qwen3-VL-2B.

## 4. Micro-batch, gradient accumulation, padding, packing

### Takeaway
Larger micro-batches raise kernel size/utilization until the GEMMs are compute-bound; accumulation only fixes optimizer batch size, not kernel efficiency. Padding waste is removed by padding-free batching or packing (FA2 varlen), with measured 1.4-2x throughput on highly variable data; with near-uniform ~560-token samples the gain should be small.

### Cited Findings
- TRL: without packing or padding-free, a large `max_length` makes training inefficient because many tokens are padding; truncate sensibly — [TRL, Reducing Memory Usage](https://huggingface.co/docs/trl/main/en/reducing_memory_usage)
- TRL packing is SFT-only, needs FlashAttention variants; strategies `bfd` (default, best-fit-decreasing, overflow discarded), `bfd_split`, `wrapped`; `padding_free=True` flattens a batch into one sequence keeping sequences intact and is set with `attn_implementation` e.g. `kernels-community/flash-attn2`; `pad_to_multiple_of` can help alignment — same TRL page
- HF padding-free + FA2 collator: up to 2x tokens/s on FLAN and ~1.4x on OrcaMath, peak memory -20% / -6%, on 8 A100-80GB, llama2-7B/mistral-7B/granite-8B, 20K-sample subsets; works by `flash_attn_varlen_func` with cu_seqlens from position_ids; optimizer-step count unchanged; needs model support for `position_ids` — [HF blog: Improving Hugging Face Training Efficiency Through Packing with Flash Attention](https://huggingface.co/blog/packing-with-FA2)
- PyTorch `expandable_segments:True` reduces allocator fragmentation (TRL notes it for long runs) — TRL page above

### Inferences
- Speedups in the HF blog depend on length variance (FLAN highly variable, 2x; OrcaMath, more uniform, 1.4x). A dataset of ~560 tokens where 256 are fixed image tokens likely has low length variance, so expect the low end or less; measure padding fraction first (1 - sum(len)/(batch x max_len)).
- Packing/padding-free for a VLM requires that the image-token position handling (mRoPE position ids in Qwen3-VL) is compatible with flattened batches; not verified (see Gaps).
- Micro-batch 4 x 560 = 2240 tokens per forward is already well past where large GEMMs saturate; raising micro-batch (with no GC, memory permitting) mainly helps small/elementwise kernels and launch overhead. Gains from mb 4 -> 8 are probably modest; test.

### Gaps
- No source on Qwen3-VL specifically supporting `padding_free`/packing with its multimodal RoPE; no source quantifying throughput vs micro-batch on a 4090.

## 5. Attention backends

### Takeaway
On a 4090, FlashAttention-2 works (Ampere/Ada/Hopper, fp16/bf16, head dim up to 256); FlashAttention-3 is Hopper-only and FA4 targets Hopper/Blackwell. PyTorch SDPA is the documented default for Qwen3-VL in HF and also dispatches to flash kernels on Ada.

### Cited Findings
- FA2 supports Ampere, Ada, Hopper "e.g. A100, RTX 3090, RTX 4090, H100", fp16/bf16, head dim up to 256; backward for head dim > 192 needs A100/H100 except head dim 256 backward works on consumer GPUs without dropout since flash-attn 2.5.5; FA3 beta is Hopper (needs CUDA >= 12.3); FA4 (CuTeDSL) is for Hopper and Blackwell — [FlashAttention repo](https://github.com/Dao-AILab/flash-attention)
- HF Qwen3-VL docs use `attn_implementation="sdpa"` in the usage example; the text model forward accepts `FlashAttentionKwargs` (so the flash path is wired in) — [HF Qwen3-VL docs](https://huggingface.co/docs/transformers/main/en/model_doc/qwen3_vl)
- HF transformers: `attn_implementation="flash_attention_2"` or the Hub-kernel form `"kernels-community/flash-attn2"` (no local install needed); `use_kernels=True` replaces supported layers with Hub kernels — [HF transformers optimization overview](https://huggingface.co/docs/transformers/en/perf_train_gpu_one) (served the inference overview) and the Qwen3-VL page
- HF warns against calling `torch.compile(model)` outside `generate()` because of per-step recompiles — same HF overview

### Inferences
- Qwen3-VL-2B head dim is likely 128 (assumption; check config), well within FA2 limits on Ada. SDPA vs FA2 difference at 560 tokens is probably small; attention is a small share of FLOPs at this length (attention FLOPs scale with T^2, ~4·T·d per token per layer vs ~24·d^2 for linears; at T=560 and d=2048 it is ~5%).
- FA2 matters more for memory/padding-free than raw speed here. xformers: no source gathered; treat as unsupported/unneeded (see Gaps).

### Gaps
- No benchmark of FA2 vs SDPA vs eager on a 4090; no xformers source; no source on whether the vision tower in Qwen3-VL needs a separate attn_implementation setting.

## 6. LoRA overhead, fused LoRA kernels, DoRA, merging, Unsloth/Liger

### Takeaway
LoRA adds small extra GEMMs (r=16) and elementwise ops per linear layer; at r16 all-linear that is a modest overhead but a visible one because the thin GEMMs are bandwidth/launch-bound. DoRA is slower than LoRA (no number found). Unsloth and Liger report 1.2-2x speedups, all measured on A100/L4/H100 vs HF baselines; not verified on a 4090 or Qwen3-VL.

### Cited Findings
- PEFT: "DoRA introduces a bigger overhead than pure LoRA, so it is recommended to merge weights for inference"; `merge_and_unload()` removes adapter inference latency; DoRA does not work with mixed-adapter batching; `lora_dropout` default 0.0; `target_modules="all-linear"` excludes the output layer — [PEFT LoRA guide](https://huggingface.co/docs/peft/main/en/developer_guides/lora)
- LoRA+ claims up to 2x finetuning speedup (learning-rate-ratio, i.e. convergence not kernel speed) — same PEFT page
- Raschka (Llama-2-7B): all-linear LoRA is ~5x the trainable params of default q/v (20.3 M vs 4.2 M) and +2.4 GB memory (16.62 vs 14.18 GB, QLoRA setting) — [Sebastian Raschka, Practical Tips for Finetuning LLMs](https://magazine.sebastianraschka.com/p/practical-tips-for-finetuning-llms)
- Liger-Kernel README: ">20% throughput", "60% memory reduction", headline setup LLaMA-3-8B, batch 8, bf16, AdamW, grad checkpointing, FSDP1, 8 A100; kernel benchmarks on H100/B200; Qwen2-VL, Qwen2.5-VL, Qwen3, Qwen3.5 listed, Qwen3-VL not listed; these are the project's own measurements — [Liger-Kernel GitHub](https://github.com/linkedin/Liger-Kernel)
- Liger paper: ~20% throughput and ~60% memory vs HF on average (abstract); per-kernel: RMSNorm ~7x faster/3x less memory, RoPE ~8x faster/3x less memory, SwiGLU ~parity speed/1.6x less memory at length 16384 (kernel microbenchmarks, A100 80GB) — [arXiv 2410.10989](https://arxiv.org/abs/2410.10989) and [HTML version](https://arxiv.org/html/2410.10989)
- Unsloth claims (own blog posts): Gemma 2.43x faster than HF+FA2 with 57.5% less VRAM; Llama 3.2 2x faster, 60% less memory vs HF+FA2 (L4 GPU for VRAM) — [Unsloth blog via search listing](https://unsloth.ai/blog/gemma)
- Independent LLaMA-Factory benchmark (A100 40GB, LoRA r8, Llama2-7B): HF 2,155 tok/s/GPU; HF+FA2 2,556 (119%); Unsloth+FA2 3,400 (158%) — [LLaMA-Factory performance comparison](https://github.com/hiyouga/LLaMA-Factory/wiki/Performance-comparison) (via search summary)

### Inferences
- Honest range for fused kernels on an HF baseline already using SDPA/FA2: 1.2-1.6x end-to-end; the 2x headlines compare against weaker baselines or longer sequences, and none are 4090/Qwen3-VL results. Treat them as upper bounds.
- Since Liger lists Qwen2.5-VL/Qwen3 but not Qwen3-VL, applying Liger to Qwen3-VL may need per-component patching (text RMSNorm/SwiGLU/FLCE) and is a correctness risk; check Liger's current model list or HF `use_liger_kernel` support before relying on it.
- LoRA merging is irrelevant for training throughput; relevant only for inference/eval speed after training.

### Gaps
- No quantitative per-layer LoRA overhead number found (e.g. % step time at r16); no DoRA slowdown figure; no Unsloth benchmark on a 4090 from a primary source (4,090 figures in secondary sites unverified); Unsloth Qwen3-VL support not checked.

## 7. Fused linear cross-entropy for 152k vocabulary

### Takeaway
Avoiding materialized logits saves a tensor of size tokens x vocab; for 152k vocab and 2240 tokens that is only ~1.4 GB fp32 per micro-batch, so FLCE helps memory/speed modestly here and its main benefit would be allowing larger micro-batches.

### Cited Findings
- Liger paper motivating example: Gemma 256k vocab, batch 8, seq 4096 gives a 16.8 GB bf16 logit tensor causing a memory spike; FLCE chunks the linear+CE; kernel benchmark of CE at vocab 163,840: ~3x faster, ~5x less memory (A100 80GB); paper does not report FLCE-specific numbers at 128k+ — [arXiv 2410.10989 HTML](https://arxiv.org/html/2410.10989)
- End-to-end full-suite results (not FLCE alone, 4xA100 80GB, Alpaca, bf16, seq len 512): LLaMA3-8B bs64 +42.8% throughput/-54.8% memory; Qwen2 bs48 +25.5%/-56.8%; Gemma-7B bs48 +11.9%/-51.8%; Mistral-7B bs128 +27%/-21%; Phi3 bs128 +17%/-13% — same paper
- Liger README says fused linear kernels reduce memory "up to 80%" for post-training losses — [Liger-Kernel GitHub](https://github.com/linkedin/Liger-Kernel)

### Inferences
- Logits for the user's micro-batch: 4 x 560 x 151,936 x 4 B (fp32 upcast) ≈ 1.36 GB, plus gradient of same size and softmax temporaries: ~3-4 GB transient; with only ~17 GB headroom this is not a constraint, so FLCE is mainly a speed (reduced memory traffic: the lm_head+CE is bandwidth-heavy) and larger-batch enabler. The lm_head GEMM itself (2048 x 152k) is ~0.31B params, roughly 18% of 1.7B LM FLOPs, so the loss path is a non-trivial share of compute.
- Another cheap approach: compute loss only on completion (label != -100) positions, lowering lm_head work if prompt tokens (256 image tokens plus prompt) dominate the sequence; not sourced, a design idea to test.

### Gaps
- No FLCE-only speed number at 152k vocab; no 4090 measurement.

## 8. bitsandbytes 4/8-bit (QLoRA) vs bf16 LoRA

### Takeaway
QLoRA is slower than bf16 LoRA (dequantization on the fly) and is only worth it when bf16 weights do not fit. For a 2B model using 6.9 GB of 24 GB, it would cost speed for no needed memory.

### Cited Findings
- Llama-2-7B on LoRA r8: bf16 LoRA 1.85 h / 21.33 GB; QLoRA NF4 2.79 h / 14.18 GB, i.e. -33% memory, +39% runtime (GPU not stated; later runs on A100) — [Raschka, Practical Tips](https://magazine.sebastianraschka.com/p/practical-tips-for-finetuning-llms)
- Other sources claim 20-60% slowdown with no reproducible methodology (conflicting, low quality) — search results: [gamut.so](https://www.gamut.so/blog/lora-vs-qlora) and others; not relied upon
- HF: bitsandbytes 4-bit with `bnb_4bit_compute_dtype=torch.bfloat16` — [HF transformers optimization overview](https://huggingface.co/docs/transformers/en/perf_train_gpu_one)

### Inferences
- A ~39% slowdown from one well-known run is the best single anchor; the effect on a 2B model with smaller GEMMs is likely proportionally similar or worse (dequant is per-forward and per-recompute; with gradient checkpointing it is done twice). Not measured.

### Gaps
- No primary benchmark on a 4090 or for small models; bitsandbytes release-note speed claims not retrieved.

## 9. torch.compile for training

### Takeaway
PyTorch reported 43% average training speedup (A100, 163 models, AMP 51%) but noted lower gains on desktop GPUs such as the 3090; dynamic shapes cause recompiles and HF warns about compiling the whole model outside generate. Treat as a possible 10-30% gain needing measurement, with minutes of compile time.

### Cited Findings
- 163 models (46 HF, 61 TIMM, 56 TorchBench), A100: 43% faster training on average; float32 21%; AMP 51%; works on 93% of models; "speedups are lower on desktop-class GPUs such as the NVIDIA 3090"; dynamic shape support was "limited and a rapid work in progress" at the 2.0 announcement (a 4-to-256 shape LM example was up to 40% faster than eager) — [PyTorch 2.0 release page](https://pytorch.org/get-started/pytorch-2.0/)
- PyTorch tuning guide: `torch.compile` fuses kernels; `reduce-overhead`/`max-autotune` modes use CUDA graphs and may raise memory and not all models compile — [PyTorch tuning guide](https://docs.pytorch.org/tutorials/recipes/recipes/tuning_guide.html)
- HF: avoid `torch.compile(model)` outside `generate()` to prevent recompiling every step — [HF optimization overview](https://huggingface.co/docs/transformers/en/perf_train_gpu_one)

### Inferences
- Variable image grids / sequence lengths in a VLM mean many distinct shapes; recompiles can erase gains unless inputs are bucketed or padded to fixed lengths (`pad_to_multiple_of`) or `dynamic=True`, plus gradient checkpointing and PEFT add compile friction. Compile only the LM decoder layers, not the full multimodal forward, if tried.
- Given the 2.0-era caveat on desktop GPUs, expect less than the 43% headline; experiment only after cheaper levers (no GC, fused optimizer, FA2/SDPA) are done.

### Gaps
- No torch.compile training numbers on a 4090 or for LoRA/VLM; compile times not found (2.0 page is dated; may not reflect current Inductor, TorchInductor support statement there only lists Volta/Ampere).

## 10. Runtime settings: TF32/cuBLAS, cudnn.benchmark, pinned memory, DataLoader

### Takeaway
These are small wins for a bf16 transformer: TF32 only affects fp32 matmuls, cudnn.benchmark only affects convolutions (the vision patch embed), and pinned memory plus workers matter if the GPU waits on data.

### Cited Findings
- `pin_memory=True` recommended with GPUs for faster async H2D copies, `num_workers` default 0 loads synchronously (tune >0 to overlap), `zero_grad(set_to_none=True)` saves memory ops; `cudnn.benchmark=True` helps convnets but can hurt with highly variable input sizes; "up to 3x" overall for mixed precision with Tensor Cores on Volta+ — [PyTorch tuning guide](https://docs.pytorch.org/tutorials/recipes/recipes/tuning_guide.html)
- TF32 flags and defaults: see section 1 — [PyTorch CUDA notes](https://docs.pytorch.org/docs/2.15/notes/cuda.html)
- TRL `pad_to_multiple_of` for alignment-friendly shapes — [TRL](https://huggingface.co/docs/trl/main/en/reducing_memory_usage)

### Inferences
- Set `torch.set_float32_matmul_precision("high")` anyway to cover any fp32 ops (loss upcast, fp32 LoRA weights if kept fp32); keep LoRA params consistent dtype to avoid autocast casts (not sourced).
- Check GPU idle time with a profiler or nvidia-smi sampling: if the 3.4 s step has gaps (image preprocessing in main process, Qwen processor on CPU), workers/pre-tokenized caching beat any kernel change.
- Fused AdamW (`optim="adamw_torch_fused"`) is a standard Trainer option; LoRA optimizer cost is tiny (r16 params), so no meaningful gain expected (not sourced).

### Gaps
- No quantified speedup for pin_memory/non_blocking or worker counts; no source on `non_blocking` specifics (tuning guide page did not cover it).

## Summary of levers, ranked by likely value for the stated run (inference, to be measured)
1. Disable gradient checkpointing (peak 6.9 GB/24 GB): up to ~25-33% less step time if matmul-bound (derived from 4N vs 6N, section 3); fall back to selective checkpointing if OOM.
2. Check data pipeline/host overhead (workers, pin_memory, image preprocessing).
3. Padding fraction: if material, padding-free/FA2 varlen (1.4-2x reported only on high-variance data).
4. Fused kernels (Liger FLCE/RMSNorm/SwiGLU): 1.1-1.4x claimed on other GPUs/models; verify Qwen3-VL support.
5. torch.compile on the LM only: uncertain, shape-sensitive.
6. Avoid QLoRA and DoRA (slower); do not expect FP8/TE to be a drop-in win at this model size.
