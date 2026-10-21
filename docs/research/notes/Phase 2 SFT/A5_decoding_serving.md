# A5: decoding and serving for Phase 2 eval sampling and Phase 3 GRPO (as of 2026-10-19)

Question: how should Phase 2 draw K samples per validation item (temperature 1, groups of G = 8) from SFT
checkpoints of SmolVLM-256M, LFM2.5-VL-450M and Qwen3-VL-2B, and what must be ready for Phase 3 GRPO rollouts?
Sub-questions: which inference engine supports which model; constrained decoding for the D10 grammar
(`objdel.geometry.text.pattern`); prefix caching, batching and speculative decoding.

Method:
- Primary sources: vLLM v0.32.0, TRL v1.16.0, transformers v5.20.1 and SGLang source and release notes, read through
  the GitHub API on 2026-10-19; engine docs; papers.
- Sandbox checks on the Mac, in `<scratchpad>/phase2/A5/`:
  - `scripts/grammar_check.py`: xgrammar 0.2.8 and llguidance 1.9.1 against all 30,843 targets × 3 tokenizers.
  - `scripts/forced_and_batch.py`: forced-token share and batched mask cost.
  - `scripts/vllm_check.py` and `scripts/hf_check.py`: vLLM 0.32.0 macOS CPU wheel, SmolVLM-256M, fp32, against HF
    teacher forcing.
  - Outputs are in `data/*.json` and `data/*.log` in the same directory.
- No GPU was used. The session's web-search budget ran out partway through. Later facts were taken from source
  code and release pages through `gh api` / WebFetch, and §9 lists what remains unverified.

Labels:
- **[src]**: read in source code or release notes at the tag named.
- **[doc]**: official documentation.
- **[paper]**: abstract or full text fetched this session.
- **[sandbox]**: measured here.
- **[inferred]**: reasoning, not measured.

## 0. Bottom line

1. **For Phase 2 eval, sample with vLLM directly, not through TRL, for all three models, SmolVLM included.**
   - The SmolVLM bug (#6294) is in TRL's hand-off to vLLM, not in vLLM.
   - Sandbox, SmolVLM-256M on the vLLM 0.32 CPU wheel:
     - chat text plus a PIL image yields exactly the HF processor's `input_ids`: 89 = 89 tokens at `longest_edge`
       512, and 355 = 355 at 1024 (our LIDC and RIGA settings);
     - unconstrained vLLM logprobs equal HF teacher-forced logprobs to max |Δ| = 6.9e-5 nats (512 px, 1,088 sampled
       tokens) and 1.4e-4 nats (1024 px, 465 tokens), in fp32;
     - TRL's current hand-off (processor-expanded ids plus the image) roughly doubles the prompt (89 → 155,
       355 → 687 tokens), and the output is garbage. So the bug is real even for single-tile SmolVLM, which answers
       Phase 1 lab T5's open question.
   - A ~150-line offline sampler (`LLM.generate`, `SamplingParams(n=8, temperature=1, top_p=1, top_k=0)`, unexpanded
     chat text plus image, `mm_processor_kwargs` = our processor kwargs) writes `eval_run` JSONL directly.
2. **vLLM cannot be installed into the pinned `.venv`.**
   - vLLM 0.32.0 pins `torch==2.13.0` and `transformers>=5.16.1,<5.20.0`; 0.31.0 pins `<5.18.0`. The repo pins
     torch 2.14.1 and transformers 5.20.1.
   - TRL 1.16 supports vLLM only from 0.22 to 0.31, and 0.31's transformers ceiling (< 5.18) excludes the ≥ 5.18 that
     vision assistant-only loss needs (configs.py).
   - Resolving `vllm==0.32.0` in a fresh venv gave torch 2.13.0, transformers 5.19.0, xgrammar 0.2.8 and llguidance
     1.7.6 [sandbox].
   - So:
     - **Phase 2:** eval sampling runs in its own environment (venv or the `vllm/vllm-openai:v0.32.0` image),
       reading merged checkpoints from disk.
     - **Phase 3:** colocate needs a re-pinned stack or server mode. **The Phase 1 assumption "LFM / Qwen use vLLM
       colocate" has never been run** (lab T5 only counted image-token runs; I22), and as pinned it cannot be
       installed.
3. **Gate sampling is unconstrained.**
   - The format gate (parse ≥ 98 %, legal ≥ 90 %) measures exactly what grammar-constrained decoding would hide.
   - On the untrained base model, constrained decoding produced `{"nodule":[]}` in 8/8 samples at 512 px, and at
     1024 px in 5/8, the rest being tiny 3-point rings such as `[[12,13],[14,15],[15,16]]`. It was grammatical and
     useless: the constraint turns format failures into degenerate geometry.
   - Run a constrained arm only as a diagnostic (§6, H-A5-5).
4. **The D10 regex is ready for xgrammar and llguidance** (H-A1-9 passed in sandbox).
   - All 30,843 targets (3,930 RIGA + 26,913 LIDC), each tokenized canonically plus the model's stop token, are
     accepted by both engines for all three tokenizers.
   - Eleven kinds of malformed texts × 50 targets are all rejected: leading zero, 1001, space, missing `holes`,
     wrong key, trailing text, float, negative, two-point ring, truncation, Arabic-Indic digit.
   - Compile time: 2–9 ms (xgrammar), 2–13 ms (llguidance).
   - Mask cost per token on one Mac core: xgrammar p50 1.4–1.9 µs, p99 5–12 µs; llguidance p50 5–12 µs, p99 43–120 µs.
   - A step of 32 groups × 8 at ~300 tokens is ~77k masks, under 0.5 s of one core with xgrammar. vLLM overlaps mask
     computation with the GPU forward pass.
5. **Do not turn on TRL's `vllm_structured_outputs_regex` for GRPO as it stands.**
   - vLLM applies the grammar bitmask to the logits *before* computing logprobs: in both raw and processed modes
     (gpu_model_runner → sampler, v0.32.0 [src]).
   - TRL's trainer recomputes logprobs without the mask, so the per-token ratio π_train/π_vLLM = Z_t ≤ 1, the
     grammar-allowed mass.
   - On the base SmolVLM the sequence log-ratio was −18 to −111 nats [sandbox]. That is a ratio below 1e-8, which
     silently zeroes those samples under TRL's default IS correction (sequence_mask, C_max 3, no C_min).
   - If constrained rollouts are ever wanted, the trainer must apply the same mask (grammar-consistent log-probs;
     the CPU cost is measured in §5) or score with the unconstrained model on purpose, as GrammarRL does.
6. **Throughput: decode dominates.**
   - Prompts are 89 tokens (SmolVLM, 512 px), ~400 (SmolVLM, 1024 px) and ~300 (LFM / Qwen at 256 image tokens).
     Answers are 71–317 tokens (target medians).
   - Prefix caching across the group works, even for copies submitted together, but it saves only prefill. In the
     CPU sandbox at 1024 px (355-token prompt, block size 128):
     - an `n=8` request reported 256 cached tokens;
     - of 8 duplicates submitted together, 7 got 256 cached tokens (only full blocks are shared);
     - at 512 px the 89-token prompt is below one CPU block, so nothing was cached.
   - The vision encoder output is shared by `mm_hash` (EncoderCacheManager [src]).
   - Speculative decoding is not worth it here (§7).
   - What matters is continuous batching (no static-batch tail) and the truncation tail: a runaway sample costs
     1,536 tokens, 5–15× a normal answer.

## 1. Our requirements, made precise

- **Volume per full validation pass at K = 64** [inferred from `data/targets` medians]:

  | Split | Items | Samples (K = 64) | Tokens per sample (SmolVLM / LFM / Qwen) | Total decode tokens |
  |---|---|---|---|---|
  | RIGA val | 63 | 4,032 | 317 / 208 / 310 | — |
  | LIDC val | 1,084 | 69,376 | 104 / 71 / 101 | — |
  | Both | 1,147 | 73,408 | — | ≈ 8.5 M (SmolVLM) |

  Add 1,536 tokens per truncated sample: 2 % truncation adds ≈ 2.3 M. Cheaper K for most items is A4's call; the
  sampler must handle any multiple of 8.
- **Sampling distribution must equal the GRPO rollout distribution**: temperature 1, top_p 1, top_k off, min_p 0,
  repetition penalty 1, `max_tokens` = `MAX_COMPLETION[model]` (1,536 / 1,024 / 1,536), and the stop set equal to
  `rewards.stop_ids` (generation_config EOS + tokenizer EOS + pad). Truncated means the last id is not a stop id.
  - Qwen3-VL-2B's `generation_config.json` ships `temperature 0.7, top_p 0.8, top_k 20` [sandbox, cached config].
  - vLLM's OpenAI server fills unset request fields from the model's generation config by default. So every field
    must be set explicitly, or the server launched with `--generation-config vllm`. vLLM's
    `generation_config` defaults to `"auto"` (loaded from the model path) [src: vllm/config/model.py v0.32.0]; H-A5-2
    verifies the effective parameters.
- **Same image preprocessing as SFT**:
  - SmolVLM `size={"longest_edge": 1024}` on RIGA (320 image tokens) and 512 on LIDC (64 tokens);
  - LFM / Qwen processor defaults (`configs.processor_kwargs`).
  - vLLM's Idefics3 processor passes `mm_processor_kwargs` through to the HF image processor
    (`get_merged_mm_kwargs`, idefics3.py v0.32.0 [src]), and the sandbox confirmed equal ids at 512 px.

## 2. Inference engines × our models (support matrix)

Versions: vLLM 0.32.0 (2026-10-18), SGLang 0.5.22 (2026-10-19), TensorRT-LLM 1.4.0rc0 (2026-10-14), llama.cpp
b11747 (2026-10-19), LMDeploy 0.18.0 (2026-09-28), MLC-LLM main, transformers 5.20.1, TRL 1.16.0 (2026-10-15).

| Engine | SmolVLM-256M (Idefics3) | LFM2.5-VL-450M (Lfm2Vl, hybrid conv + attention) | Qwen3-VL-2B | Grammar engines | TRL GRPO integration | Fit for us |
|---|---|---|---|---|---|---|
| **vLLM 0.32** | Native `Idefics3ForConditionalGeneration` since 2024-11 (#9767); LoRA ✓; sandbox ✓ exact HF parity | Native since 2026-01-08 (#31758); LoRA ✓, PP ✓; `IsHybrid` → prefix cache in mamba "align" mode | Native since v0.11.0 (#24727); LoRA ✓; 4× smaller M-RoPE cache in 0.32 (#59853) | xgrammar 0.2.8 (pinned), llguidance <1.8, outlines_core 0.2.14, lm-format-enforcer 0.11.3; `auto` default | colocate + server; TRL supports ≤ 0.31 only | **Use** (standalone for Phase 2; Phase 3 per §8) |
| **HF transformers 5.20.1 `generate`** | ✓ (current SmolVLM GRPO path) | ✓ | ✓ | `logits_processor` only (xgrammar / llguidance ship HF processors) | `use_vllm=False` path; static batch, `disable_compile=True` | Reference and fallback; slow tail |
| **transformers continuous batching** (`generate_batch`) | ✗: CB has no image path in 5.20.1; TRL raises "does not support multimodal models" | ✗ | ✗ | — | `use_transformers_continuous_batching` is text-only | Reject |
| **SGLang 0.5.22** | No native Idefics3/SmolVLM model file (only `idefics2.py`, the MiniCPM vision part); generic Transformers backend has a `MultiModalMixin`, untested here | Native `lfm2_vl.py` (since 2026-04) | Native `qwen3_vl.py` | xgrammar (default), llguidance, outlines | None in TRL; verl / AReaL / slime use it (no SmolVLM) | Alternative for LFM / Qwen only |
| **TensorRT-LLM 1.4.0rc0** | ✗ | ✗ | ✓ PyTorch backend, guided decoding ✓, KV reuse ✓, MM embedding cache | xgrammar, llguidance | None | Reject (one model, build cost) |
| **llama.cpp b11747** | ✓ GGUF `ggml-org/SmolVLM-256M-Instruct-GGUF` | ✓ (`LiquidAI/LFM2.5-VL-450M-GGUF`, conversion `lfm2.py`) | ✓ (`qwen3vl.cpp`) | GBNF; llguidance (`-DLLAMA_LLGUIDANCE=ON`) | None | Deployment / CPU demos only |
| **LMDeploy 0.18** | ✗ | ✗ | ✓ (2B–235B) | its own | None | Reject |
| **MLC-LLM** | ✗ | ✗ | ✗ (has `qwen2_5_vl`, `qwen35`) | xgrammar | None | Reject |

Notes:
- The vLLM supported-models page lists Idefics3, SmolVLM (for SmolVLM2), Lfm2Vl and Qwen3VL as T + I generative
  models with LoRA [doc].
- SmolVLM-256M-Instruct's config says `Idefics3ForConditionalGeneration` [sandbox, cached config].
- vLLM's own SmolVLM example uses `mm_processor_kwargs={"max_image_size": {"longest_edge": 384}}`, the precedent
  for passing processor kwargs [src].

## 3. TRL ↔ vLLM: what actually happens (TRL 1.16.0)

- **Colocate** builds `LLM(model=model.name_or_path, gpu_memory_utilization=…(0.3 default in GRPOConfig),
  max_num_batched_tokens=4096, logprobs_mode="processed_logprobs", distributed_executor_backend="external_launcher",
  enable_sleep_mode=…)`. It sends **`n=1` per prompt, with the G copies of each prompt as separate requests**
  ({"prompt_token_ids", "multi_modal_data"}) [src: trl/generation/vllm_generation.py].
  - No `mm_processor_kwargs` and no `hf_overrides` are forwarded (issue #6776 open).
  - vLLM therefore uses the processor config **in the model directory**. Save the processor with our kwargs next to
    the SFT checkpoint, or SmolVLM on RIGA gets the default 2,048-px tiling [inferred; H-A5-3].
- **Server mode** sends one request with `n = num_generations` per G-th prompt.
  - Images go through `/inference/v1/generate`, which needs `--enable-scale-out` (vLLM ≥ 0.30), plus
    `--logprobs-mode processed_logprobs`.
  - Weights are pushed over NCCL, or **`ipc` when trainer and server share a GPU** [doc: TRL vLLM integration].
  - Server mode on one GPU in a **separate venv** is the clean way around the torch / transformers pin conflict
    [inferred; H-A5-7].
- **#6294 (SmolVLM)**: open; fix PRs #6718 / #6719 open. Maintainer, 2026-10-13: "low priority: it only affects
  Idefics3 / SmolVLM… changing what we send to vLLM for every VLM needs care" [src: PR comments].
  - The PR's mechanism is to render the chat text untokenized and send the tokenizer-only ids, keeping the expanded
    ids for the loss. Our sandbox validates that mechanism: vLLM's own processing of text plus image gives HF's ids
    exactly.
  - Expect no merge in Phase 2. For Phase 3, either subclass `_generate_single_turn` / `_tokenize_prompts` locally
    (≈ the PR diff), or keep HF generation for SmolVLM.
- **Encoder cache after weight sync (TRL #7198, open)**: `sync_weights` resets the prefix cache but not vLLM's
  encoder cache. With a trainable vision tower (full fine-tuning of SmolVLM and LFM) and recurring images, rollouts
  use stale vision embeddings. Call `llm.reset_encoder_cache()` after each sync (callback), or freeze the vision
  tower in GRPO.
- **HF path**: `GenerationConfig(…, top_k=self.top_k (0), suppress_tokens=[image_token_id], disable_compile=True)`.
  - Completions are cut after the first EOS.
  - Static batches, so a batch runs until its longest member stops: a runaway sample holds the whole batch for 1,536
    steps.
  - `cache_implementation` is passable, but compile is disabled [src: grpo_trainer.py L1155–1173, L1840–1890].
- **IS correction defaults**:
  - `vllm_importance_sampling_correction=True`, `mode="sequence_mask"`, `clip_max=3.0`, `clip_min=None` [src:
    grpo_config.py].
  - TRL warns (does not refuse) when a server's `--logprobs-mode` is not `processed_logprobs` [src:
    grpo_trainer.py].

## 4. Constrained decoding: engines and our grammar

**The grammar** (D10, `pattern(names)`) is a Rust-compatible regex with no lookarounds. vLLM rejects NUL bytes in
regexes and compiles them with `xgr.Grammar.from_regex` [src: backend_xgrammar.py]. Sandbox, over all three tokenizers
and both schemas:

| Tokenizer | Vocab | Stop | Targets accepted (xgr / llg) | Compile xgr / llg (ms) | xgr mask µs p50 / p99 | llg mask µs p50 / p99 | Tokens with 1 legal choice |
|---|---|---|---|---|---|---|---|
| SmolVLM-256M | 49,280 | 49279 `<end_of_utterance>` | 30,843 / 30,843 | 5.0 / 3.2 | 1.4 / 9.2 (RIGA), 1.4 / 6.5 (LIDC) | 6.2 / 104, 6.3 / 112 | 12.3–12.5 % |
| LFM2.5-VL-450M | 64,400 | 7 `<\|im_end\|>` | 30,843 / 30,843 | 8.9 / 7.4 | 1.9 / 11.2, 1.8 / 11.7 | 12.2 / 104, 11.5 / 120 | 34–35 % |
| Qwen3-VL-2B | 151,669 | 151645, 151643 | 30,843 / 30,843 | 7.1 / 13.2 | 1.5 / 5.2, 1.7 / 5.4 | 5.3 / 43, 6.2 / 54 | 12.0–12.3 % |

Notes on the table:
- Mask timings come from the first 200 targets per schema (21.8k–82.5k masks).
- Rare first-call outliers reached 33–66 ms (xgrammar) and 18–111 ms (llguidance), once per matcher warm-up.
- Legal tokens per step: median 10 (digit-per-token vocabularies) or 3 (LFM). At number positions LFM allows up to
  1,000 legal tokens: its three-digit tokens.
- Malformed rejection: 11 kinds × 50 targets per schema, 0 accepted by either engine, all tokenizers.

**Engine notes:**
- **xgrammar**: vLLM's and SGLang's default; v0.32 pins 0.2.8, the version tested here.
  - XGrammar claims "up to 100x speedup" and near-zero end-to-end overhead [paper 2411.15100].
  - XGrammar-2 (2026-01) adds JIT compilation, cross-grammar caching, an Earley-based matcher and "repetition
    compression", with >6× over other engines [paper 2601.04426]. That last item is relevant if we ever bound
    repetitions (`{2,999}`).
- **llguidance**: "~50 µs of CPU time per token for a 128k tokenizer", negligible startup; in vLLM ≥ 0.8.2, SGLang ≥
  0.4.4, llama.cpp, TensorRT-LLM [doc: README]. It is ~4–10× slower per mask than xgrammar here, but still far below
  a GPU decode step.
- **Outlines / lm-format-enforcer**: in vLLM. lm-format-enforcer uses Python `re` syntax, the others Rust-style
  [doc: vLLM structured outputs]. No reason to prefer them.
- **What the grammar cannot do**:
  - It cannot make polygons legal: simplicity, ring orientation, holes inside outer and part disjointness are
    checked by polycore, not by the grammar.
  - It cannot prevent truncation. The grammar admits unbounded vertex lists, so a looping model still hits
    `max_tokens`.
  - It does not enforce the parser caps (≤ 1,000 vertices, ≤ 10 rings). These are unreachable within 1,536 tokens
    anyway (≈ 187 vertices at 8.2 tokens each).
- **Stop tokens**: vLLM 0.32 passes `sampling_params.all_stop_token_ids` to the matcher as `override_stop_tokens`
  [src: v1/structured_output/__init__.py]. So SmolVLM's `<end_of_utterance>` (generation config), which differs from
  the tokenizer's EOS, terminates correctly. In standalone xgrammar or llguidance use, pass the stop ids explicitly
  (as `grammar_check.py` does).

## 5. Constrained decoding and RL: literature and our numbers

- **Distortion**: grammar-constrained decoding samples from q(y) ∝ ∏ p(y_t | y<t) · 1[allowed], renormalized per
  step. That is not p conditioned on grammaticality: outputs are "grammatical but appear with likelihoods that are
  not proportional" to the LLM's (Grammar-Aligned Decoding, ASAp) [paper 2405.21047, NeurIPS 2024].
  - Our base-model demo is the extreme case: q put all 8 samples on `{"nodule":[]}`, the shortest grammatical answer
    [sandbox].
- **Format restriction cost**: "a significant decline in LLMs reasoning abilities under format restrictions…
  stricter format constraints generally lead to greater performance degradation" [paper 2408.02442]. This is less
  relevant for us: the answer *is* the format, and nothing is reasoned outside it.
- **RL with constrained rollouts**: GrammarRL (Tuccio et al., 2026-09) generates under the constraint but computes
  "every sequence probability used to score or train on these candidates… with the *unconstrained* model". RLOO, up
  to +22.8 BLEU over greedy constrained decoding [paper 2609.39869].
- **RL to learn the format instead**: ThinkJSON trains schema adherence with GRPO (1.5B; ~20 h on 8×H100)
  [paper 2502.14905]. SchemaRL: RL plus a fine-grained schema validator improves valid-JSON rates [paper 2502.18878].
  Our v4 reward already has `format_ok` / `legal` terms. That is the learn-it path, and D16 chose it.
- **The TRL mismatch, measured** (vLLM 0.32 CPU, base SmolVLM-256M, fp32, 8 samples each; `hf_check.py`):

  | Arm | px | Tokens | vLLM − HF logprob per token | Sequence log-ratio log(q/p) |
  |---|---|---|---|---|
  | Unconstrained | 512 | 1,088 | mean 2.6e-7, max \|Δ\| 6.9e-5 | 0.000 for all 8 |
  | Unconstrained | 1024 | 465 | mean −4.2e-7, max \|Δ\| 1.4e-4 | 0.000 for all 8 |
  | Constrained (D10 regex) | 512 | 158 | mean 1.50, median 0.10, max 13.5 nats | 18.0–103.4 nats |
  | Constrained (D10 regex) | 1024 | 169 | mean 2.33, median 0.32, max 17.0 nats | 18.3–110.5 nats |

  TRL's IS weight is p/q = e^(−18…−111). Sequence_mask with C_max = 3 never masks a ratio ≤ 1, so these samples get
  near-zero weight. With the correction off, the update uses p's gradient on q's samples: biased.
  - After SFT the gap should be small wherever the model already follows the grammar (Z_t ≈ 1). It is exactly large
    on the samples where constraint changed the outcome, which are the ones whose learning signal matters (H-A5-6
    measures it on the SFT checkpoint).
- **Grammar-consistent GRPO** (option, not a plan):
  - Apply the same xgrammar mask in the trainer's log-prob computation (log q = logit − logsumexp over allowed
    tokens), so π_train and π_vLLM both equal q.
  - Cost: masks for one 256-rollout step took 0.24–1.46 s on one core [sandbox], including writing dense bitmasks
    (0.7–3.4 GB padded, so materialize per micro-batch, or keep sparse allowed lists: ≤ 15 tokens per step except
    LFM number positions).
  - It needs a TRL subclass (`_get_per_token_logps_and_entropies`) and per-sample matchers. The policy then never
    learns the format, so deployment must also constrain.

## 6. Does constraining at eval hide what the gates should see?

Yes, by construction:
- Under the constraint, `parse` becomes 1 − truncated, and `error/<code>` collapses to `truncated`.
- The parse ≥ 98 % gate becomes untestable, and legal ≥ 90 % measures only geometry.

Therefore:
1. **Gate runs are unconstrained** (same as GRPO's planned rollouts).
2. **Diagnostic constrained arm** (cheap, same sampler, `StructuredOutputsParams(regex=pattern(names))`), on the same
   items and K. Report:
   - legal and IoU_tc under the constraint;
   - the share of unconstrained format failures whose item reaches IoU ≥ 0.5 under the constraint ("repairable");
   - E[−log Z] per token, from teacher forcing (H-A5-6).

   These numbers decide whether constrained or grammar-consistent rollouts are worth building for Phase 3. They are
   never gate inputs.

## 7. Sampling efficiency

- **Prefix caching across G = 8**:
  - vLLM V1 enables prefix caching by default (`enable_prefix_caching: bool = True`, config/cache.py v0.32 [src]).
    The image hash enters the block hash [doc: design/prefix_caching].
  - Sandbox (CPU backend, block size 128), as `num_cached_tokens`:

    | Prompt | `n=8` request | 8 duplicates submitted together | Repeat after completion |
    |---|---|---|---|
    | 1024 px, 355 tokens | 256 | 0 for the first, 256 for the other 7 | 256 |
    | 512 px, 89 tokens | 0 (less than one block) | 0 for all | — |

  - So same-step copies do share full blocks. On GPU the block size is 16, so nearly the whole prompt is shared
    [inferred; H-A5-4]. TRL colocate's G separate copies therefore pay prefill about once.
  - Prompt share: SmolVLM-512 has 89 prompt tokens versus ~104 answer tokens, so prefill is ≲ 15 % of FLOPs, and
    decode, being bandwidth-bound, dominates time. Expected gain ≤ 10 % [inferred; H-A5-4 measures it on GPU].
- **Vision encoder**: shared across requests with the same `mm_hash`; "Cache is enabled to share embeddings of same
  multimodal data item" [src: v1/core/encoder_cache_manager.py]. The `mm_processor_cache_gb` = 4 default also skips
  re-running the HF image processor. So G copies cost one vision-encoder pass.
- **Batching**:
  - vLLM continuously batches up to `max_num_seqs`; the default depends on hardware and usage context, so set it
    explicitly. TRL colocate caps `max_num_batched_tokens=4096`, which limits prefill chunking, not decode
    concurrency.
  - For the standalone sampler, submit all items at once and let the scheduler fill.
  - HF static batching wastes `1 − mean_len / max_len` of each batch. With one truncated sample in a batch of 16,
    that is up to 80–93 % [inferred from 104–317 versus 1,536 tokens].
- **KV memory: a constraint only for Qwen in colocate mode** [inferred from cached model configs]:
  - Qwen3-VL-2B: 28 layers × 8 KV heads × 128 dims × 2 (K, V) × 2 bytes = 112 KiB per token. So 256 sequences ×
    1,900 tokens ≈ 54 GB at full length, ≈ 17 GB at mean lengths (~600 tokens).
  - SmolVLM-256M: 30 × 3 × 64 × 2 × 2 = 23 KiB per token.
  - LFM2.5-VL-450M: 6 attention layers × 8 × 64 × 2 × 2 = 12 KiB per token, plus conv state.
  - On a 24 GB card in colocate mode at 0.3 utilization (7.2 GB, including vLLM's own 4.3 GB bf16 weight copy),
    Qwen's KV pool is ~3 GB ≈ 27k tokens: ~14 full-length or ~45 mean-length sequences out of 256. So Qwen needs
    sleep mode, a higher utilization or server mode. The small models do not.
- **Speculative decoding — reject for Phase 2/3**:
  - Rollouts and evals run hundreds of concurrent sequences. Speculative decoding spends extra compute to cut
    per-sequence latency, which pays at low batch, not for throughput at high batch [inferred].
  - No EAGLE / MTP drafter is known for these three models; the 0.32 release notes name other targets [src; no
    exhaustive check].
  - N-gram / prompt lookup has nothing to copy: coordinates are not in the prompt.
  - Grammar-forced tokens are 12 % (SmolVLM, Qwen) and 34 % (LFM) of steps [sandbox]. Jump-forward could skip them,
    but vLLM does not implement it. It also changes tokenization relative to the policy's own (RL-unsafe), and
    xgrammar's `find_jump_forward_string` averaged only 0.6–1.4 characters per step here.

## 8. Recommendations (ranked)

1. **Phase 2 eval sampler = standalone vLLM 0.32.0, offline `LLM.generate`, separate environment.**
   - Lab code in `experiments/phase2/`, not `src/`.
   - Inputs: merged checkpoint dir, the model's processor kwargs (`configs.processor_kwargs`), chat text rendered by
     our processor with `tokenize=False` (identical template to SFT/GRPO), PIL image.
   - `SamplingParams(n=8, temperature=1.0, top_p=1.0, top_k=0, min_p=0.0, repetition_penalty=1.0,
     max_tokens=MAX_COMPLETION[model], stop_token_ids=sorted(stop_ids), seed=…, logprobs=0)`, with K/8 requests per
     item.
   - Output: `eval_run` JSONL with `truncated` = TRL's rule on `token_ids`, cross-checked against
     `finish_reason == "length"`.
   - `dtype=bfloat16` on Ampere+ (as training), `enable_prefix_caching` on (free), `gpu_memory_utilization` 0.85
     (no trainer on the GPU).
   - Must pass H-A5-1 / H-A5-2 before its numbers count.
2. **Keep `objdel.train.sample` (TRL / HF generate) as the reference and fallback.**
   - Use it for the parity check (H-A5-2) and wherever the vLLM env is unavailable.
   - Raise `per_device_eval_batch_size` to the largest multiple of 8 that fits, to amortize the static tail.
   - It is also the only path that reproduces SmolVLM's Phase 3 rollouts exactly if HF generation is kept there.
3. **Gates unconstrained; constrained arm as a diagnostic** (§6), with the same sampler and the D10 regex via
   xgrammar (`backend="auto"` picks it; llguidance is the fallback, also verified).
4. **Phase 3 preparation, in order:**
   - (a) Decide the environment: either re-pin to TRL 1.16 + vLLM 0.32 + transformers 5.19 + torch 2.13 (one env,
     colocate; TRL warns on 0.32), or keep the training env and run vLLM **server mode on the same GPU in its own
     env** with `ipc` weight transfer. H-A5-7 decides.
   - (b) Save the processor config with our kwargs into every checkpoint dir that vLLM loads (colocate cannot
     forward `mm_processor_kwargs`).
   - (c) Reset vLLM's encoder cache after each weight sync when the vision tower trains (#7198).
   - (d) For SmolVLM, either HF generation (measure step time first: Phase 2 cost gate) or a local subclass
     implementing #6718's unexpanded-ids hand-off.
   - (e) Rollouts unconstrained; never `vllm_structured_outputs_regex` without grammar-consistent trainer logprobs
     (§5).
   - (f) Keep `vllm_importance_sampling_correction` on (bf16 vLLM versus bf16 trainer kernels differ even
     unconstrained).
5. **Grammar artifacts to keep**: `pattern(names)` as is, and the sandbox `grammar_check.py` as a regression check
   whenever the grammar or a tokenizer changes (now in `experiments/phase2/lab/A5/`).

## 9. Rejected alternatives

| Alternative | Why rejected |
|---|---|
| Install vLLM into the pinned `.venv` | Pins conflict: vLLM 0.32 torch==2.13.0, transformers <5.20; repo torch 2.14.1, transformers 5.20.1; TRL 1.16 ≤ vLLM 0.31 (transformers <5.18) |
| Sample Phase 2 evals through TRL + vLLM colocate | Breaks SmolVLM (#6294, reproduced); needs a trainer object and the conflicting env; colocate cannot pass `mm_processor_kwargs` |
| transformers continuous batching | Text-only in 5.20.1; TRL raises for processors |
| SGLang | No native SmolVLM; no TRL integration; adds a second engine for no gain over vLLM on LFM / Qwen |
| TensorRT-LLM, LMDeploy, MLC | Support at most Qwen3-VL; engine builds per checkpoint |
| llama.cpp for eval | Per-checkpoint GGUF conversion, limited batched sampling; keep for deployment demos |
| Constrained decoding at the gates | Hides the format failures the gate measures (§6) |
| `vllm_structured_outputs_regex` in GRPO as is | Logprob mismatch: TRL weights constrained samples by ∏ Z_t (measured 1e-8 to 1e-48 on base SmolVLM) or trains biased |
| ASAp / SMC grammar-aligned sampling | Fixes distortion at multiple-sample cost; irrelevant if we sample unconstrained |
| Speculative decoding, jump-forward | High-batch small-model decode; no drafters; forced tokens ≈ 12–34 % and RL-unsafe retokenization |
| Outlines, lm-format-enforcer | Slower or Python-regex dialect; no benefit over xgrammar / llguidance |

## 10. Hypotheses to test

Each hypothesis has a test design and a pass criterion. P2 means Phase 2 (first Vast box) and P3 means Phase 3
preparation.

**H-A5-1 (P2, GPU): the standalone vLLM sampler sees the same prompt as SFT for every model × dataset.**
- **Test:** for 32 val items per (model, dataset), with the SFT processor kwargs, compare `RequestOutput.prompt_token_ids` with the HF processor's `input_ids` (same chat text and image).
- **Pass:** identical for 100 % of items in all 6 combinations (sandbox: SmolVLM at 512 and 1024 ✓, one synthetic image each).

**H-A5-2 (P2, GPU): the vLLM samples follow the HF / TRL distribution.**
- **Test, part (a), log-probability parity:** teacher-force 64 vLLM samples per model (from an SFT checkpoint) through HF in the training dtype. Log the per-token |Δ logprob|.
- **Pass (a):** the median per-token |Δ| ≤ 0.02 nats in bf16 (≤ 1e-3 in fp32), and the 99th-percentile per-sequence |Δ| ≤ 0.5 nats.
- **Test, part (b), distribution parity:** sample 32 items × K = 64 with both samplers. Compare parse, legal, median IoU_tc, mean length and the duplicate share.
- **Pass (b):** every metric's difference lies within its bootstrap 95 % CI, and the effective sampling parameters echoed by the engine are T = 1, top_p = 1, top_k off. The Qwen generation_config defaults must not leak into vLLM.

**H-A5-3 (P3): TRL colocate picks up the processor kwargs from the checkpoint dir.**
- **Test:** save SmolVLM / LFM / Qwen SFT checkpoints with `processor.save_pretrained`, including the kwargs. Run 1 GRPO step in colocate. Assert that the vLLM prompt length equals the trainer's expanded prompt length after un-expansion.
- **Pass:** equal for all three models. Otherwise, add `mm_processor_kwargs` through a TRL subclass.

**H-A5-4 (P2, GPU): prefix and encoder caching across a G = 8 group.**
- **Test:** on a 4090, submit 256 requests (32 items × 8) in three ways: (i) all together as `n=8`, (ii) all together as duplicates, (iii) staggered (the first copy, then the other 7). Record `num_cached_tokens`, encoder-cache hits (vLLM metrics) and wall time.
- **Pass:** in (i) and (ii), ≥ 7 of 8 copies show `num_cached_tokens` ≥ prompt length − 16. The CPU sandbox showed sharing within a step at block size 128. Adopt (iii) only if it is ≥ 10 % faster end-to-end; otherwise submit everything at once.

**H-A5-5 (P2): constrained-arm diagnostic on the first SFT checkpoint.**
- **Test:** same items and K, with and without the regex. Report the parse, legal, IoU_tc and truncation share of both, and the "repairable" share (§6).
- **Decision rule, not a pass/fail:**
  - If the unconstrained parse is ≥ 98 % and the constrained legal ≤ the unconstrained legal + 2 points, then constrained rollouts are irrelevant for Phase 3.
  - If the unconstrained parse is < 98 % but the constrained legal is ≥ 90 % with IoU within 0.02, then build H-A5-6 and the grammar-consistent GRPO option.

**H-A5-6 (P3): the grammar mass gap on SFT samples is small.**
- **Test:** for 2,000 unconstrained samples of the SFT checkpoint, teacher-force and compute −log Z_t per token, where Z_t is the softmax mass on xgrammar-allowed tokens (masks in ~2 µs per token). Report the per-sequence Σ.
- **Pass:** median Σ ≤ 0.05 nats, and ≥ 95 % of sequences ≤ 1 nat. Then the IS weights under constraint would be ≥ e^−1, and constrained rollouts plus TRL's IS correction would be benign. Otherwise grammar-consistent log-probabilities are mandatory before any constrained rollouts.

**H-A5-7 (P3, GPU): a resolvable Phase 3 environment exists.**
- **Test:** two candidates.
  - (A) A single venv with trl 1.16.0, vllm 0.32.0, transformers 5.19.x and torch 2.13.0: run `make test`, the K12 smoke (4 GRPO steps) with LFM and with Qwen in colocate, and reward and target hash checks.
  - (B) The training venv as pinned, plus a vLLM 0.32 server in its own venv on the same GPU (`--weight-transfer-config '{"backend":"ipc"}'`, `--enable-scale-out`, `--logprobs-mode processed_logprobs`): 4 GRPO steps.
- **Pass:** a candidate passes when its steps complete, the weight sync changes the vLLM logprobs, and the hashes are unchanged. Prefer (A) if both pass and (A)'s step time is ≤ (B)'s.

**H-A5-8 (P2, GPU): HF generate is too slow for full-K evals; vLLM is not.**
- **Test:** for SmolVLM-256M SFT on the 4090, sample LIDC val (1,084 items × K = 16) with `objdel.train.sample` (eval batch 64) and with the vLLM sampler. Record wall time and decode tokens per second.
- **Pass (to make vLLM the default):** vLLM ≥ 5× faster. If both finish a K = 64 full pass in under 30 min, keep HF for simplicity.

**H-A5-9 (P3, GPU): vision-tower training plus vLLM needs an encoder-cache reset.**
- **Test:** full-FT SmolVLM or LFM GRPO, 3 steps. Compare vLLM logprobs for a recurring image with and without `reset_encoder_cache()` after the sync, against an HF teacher-forced reference.
- **Pass:** with the reset, |Δ| matches H-A5-2. Without it, a drift larger than H-A5-2's bound confirms #7198 for us.

## 11. Sandbox log (reproduce)

All files are in `<scratchpad>/phase2/A5/`, where `<scratchpad>` is this session's scratch directory.
The four scripts are also copied, ruff-clean, to `experiments/phase2/lab/A5/` (`grammar_check.py`,
`forced_and_batch.py`, `vllm_check.py`, `hf_check.py`). They need the scratch venvs described below, not the repo
`.venv`.

Setup and data export:
- Venvs: `venv` (xgrammar 0.2.8, llguidance 1.9.1, transformers 5.20.2, torch 2.14.1) and `vvenv` (vLLM
  0.32.0+cpu macOS wheel from the GitHub release; it resolved torch 2.13.0, transformers 5.19.0, llguidance 1.7.6).
- Targets were exported with the repo venv under `-I` to `data/targets_texts.json` (texts and patterns only). Parse
  check: 0 failures.

Commands:
1. `HF_HUB_OFFLINE=1 venv/bin/python -I scripts/grammar_check.py data/targets_texts.json`
   - Took 215 s; peak RSS ≈ 0.5 GB.
   - Writes `data/grammar_report.json`.
2. `HF_HUB_OFFLINE=1 venv/bin/python -I scripts/forced_and_batch.py data/targets_texts.json`
   - Measures forced-token shares and the 256-rollout batched mask cost (0.24–1.46 s single-thread).
3. `HF_HUB_OFFLINE=1 EDGE=512 vvenv/bin/python -I scripts/vllm_check.py data/targets_texts.json data/vllm_512.json`
   - Took 242 s; peak footprint 0.79 GB.
   - **Needs a local patch:** in the vLLM 0.32.0 macOS CPU wheel, the override `_compute_slot_mapping_kernel_impl`
     in `vllm/utils/cpu_triton_utils.py` receives the Triton kernel's arguments (starting with `num_reqs`)
     positionally. Its signature lacks `num_reqs`, so the first decode step dies with "Boolean value of Tensor…
     ambiguous". The sandbox patched that function to pick its tensors and block size by position. This is a CPU-only
     packaging bug, irrelevant on CUDA.
4. `vvenv/bin/python -I scripts/hf_check.py data/vllm_512.json` produces the parity and mismatch numbers in §5.

Results at 512 px:
- HF prompt 89 tokens; vLLM text + image prompt 89 tokens, ids identical.
- Expanded-ids hand-off: 155 tokens, gibberish output.
- `n=8`: 0 cached tokens. 8 duplicates: 0 cached tokens each. The prompt is under one 128-token CPU block.
- Unconstrained 8 × ≤ 160 tokens: 51 s on CPU. Constrained: 13.7 s, all `{"nodule":[]}`.

Results at 1024 px (`EDGE=1024`; took 332 s):
- Prompts: 355 = 355 tokens, ids identical. Expanded-ids hand-off: 687 tokens, gibberish output.
- Caching: `n=8` got 256 cached tokens; duplicates got [0, 256 × 7]; the repeat after completion got 256. The block
  size is 128.
- Constrained samples: 5/8 `{"nodule":[]}`; the others are tiny rings.
- Parity in `data/vllm_1024_hf.json`: unconstrained max |Δ| 1.4e-4 over 465 tokens; constrained sequence log-ratios
  of 18.3–110.5 nats.

## Sources

- TRL v1.16.0 source:
  - [grpo_config.py](https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_config.py) (vLLM fields, IS
    correction defaults, top_k = 0, continuous batching);
  - [vllm_generation.py](https://github.com/huggingface/trl/blob/v1.16.0/trl/generation/vllm_generation.py)
    (colocate `LLM(...)` arguments, n = 1 per copy, `StructuredOutputsParams(regex=…)`, `processed_logprobs`);
  - [grpo_trainer.py](https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_trainer.py) (CB refuses
    multimodal; HF GenerationConfig; EOS masking);
  - [pyproject.toml](https://github.com/huggingface/trl/blob/v1.16.0/pyproject.toml) (`vllm>=0.22.0,<=0.31.0`).
- TRL issues and PRs:
  - [#6294](https://github.com/huggingface/trl/issues/6294) and [PR #6718](https://github.com/huggingface/trl/pull/6718)
    (open; maintainer comment 2026-10-13); [PR #6719](https://github.com/huggingface/trl/pull/6719);
  - [#7198](https://github.com/huggingface/trl/issues/7198) (encoder cache not reset);
  - [#6776](https://github.com/huggingface/trl/issues/6776) (no hf_overrides in colocate).
- TRL docs: [vLLM integration](https://huggingface.co/docs/trl/main/en/vllm_integration) (modes, ipc weight transfer
  on a shared GPU, `--enable-scale-out`, `processed_logprobs`, version range).
- vLLM v0.32.0:
  - [release notes](https://github.com/vllm-project/vllm/releases/tag/v0.32.0) (transformers <5.20, macOS CPU wheel);
  - [requirements/common.txt](https://github.com/vllm-project/vllm/blob/v0.32.0/requirements/common.txt) and
    [cuda.txt](https://github.com/vllm-project/vllm/blob/v0.32.0/requirements/cuda.txt) (torch==2.13.0, xgrammar
    0.2.8, llguidance <1.8);
  - [v0.31.0 common.txt](https://github.com/vllm-project/vllm/blob/v0.31.0/requirements/common.txt)
    (transformers <5.18);
  - [gpu_model_runner.py](https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/v1/worker/gpu_model_runner.py)
    (bitmask before sampling) and [sampler.py](https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/v1/sample/sampler.py)
    (raw / processed logprobs);
  - [backend_xgrammar.py](https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/v1/structured_output/backend_xgrammar.py),
    [structured_output/__init__.py](https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/v1/structured_output/__init__.py)
    (stop tokens);
  - [encoder_cache_manager.py](https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/v1/core/encoder_cache_manager.py),
    [config/cache.py](https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/config/cache.py),
    [config/multimodal.py](https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/config/multimodal.py);
  - [idefics3.py](https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/model_executor/models/idefics3.py),
    [lfm2_vl.py](https://github.com/vllm-project/vllm/blob/v0.32.0/vllm/model_executor/models/lfm2_vl.py) (IsHybrid),
    [vision_language_offline.py](https://github.com/vllm-project/vllm/blob/v0.32.0/examples/generate/multimodal/vision_language_offline.py);
  - [PR #31758](https://github.com/vllm-project/vllm/pull/31758) (LFM2-VL);
    [v0.11.0 notes](https://github.com/vllm-project/vllm/releases/tag/v0.11.0) (Qwen3-VL).
- vLLM docs: [supported models](https://docs.vllm.ai/en/latest/models/supported_models/),
  [structured outputs](https://docs.vllm.ai/en/latest/features/structured_outputs/),
  [prefix caching design](https://docs.vllm.ai/en/latest/design/prefix_caching/).
- transformers v5.20.1:
  [continuous_batching/continuous_api.py](https://github.com/huggingface/transformers/blob/v5.20.1/src/transformers/generation/continuous_batching/continuous_api.py)
  (no image inputs).
- SGLang: [models dir](https://github.com/sgl-project/sglang/tree/main/python/sglang/srt/models),
  [structured outputs](https://docs.sglang.io/advanced_features/structured_outputs.html),
  [releases](https://github.com/sgl-project/sglang/releases).
- TensorRT-LLM: [supported models](https://github.com/NVIDIA/TensorRT-LLM/blob/main/docs/source/models/supported-models.md),
  [guided decoding](https://github.com/NVIDIA/TensorRT-LLM/blob/main/docs/source/features/guided-decoding.md).
- llama.cpp: [multimodal.md](https://github.com/ggml-org/llama.cpp/blob/master/docs/multimodal.md),
  [mtmd models](https://github.com/ggml-org/llama.cpp/tree/master/tools/mtmd/models);
  [LFM2.5-VL-450M-GGUF](https://huggingface.co/LiquidAI/LFM2.5-VL-450M-GGUF).
- LMDeploy: [supported models](https://github.com/InternLM/lmdeploy/blob/main/docs/en/supported_models/supported_models.md).
  MLC-LLM: [model dir](https://github.com/mlc-ai/mlc-llm/tree/main/python/mlc_llm/model).
- Grammar engines:
  - XGrammar [arXiv 2411.15100](https://arxiv.org/abs/2411.15100);
    XGrammar-2 [arXiv 2601.04426](https://arxiv.org/abs/2601.04426);
  - [llguidance README](https://github.com/guidance-ai/llguidance);
  - JSONSchemaBench [arXiv 2501.10868](https://arxiv.org/abs/2501.10868).
- Constrained decoding and RL:
  - Grammar-Aligned Decoding [arXiv 2405.21047](https://arxiv.org/abs/2405.21047);
  - Let Me Speak Freely? [arXiv 2408.02442](https://arxiv.org/abs/2408.02442);
  - GrammarRL [arXiv 2609.39869](https://arxiv.org/html/2609.39869);
  - ThinkJSON [arXiv 2502.14905](https://arxiv.org/abs/2502.14905);
  - SchemaRL [arXiv 2502.18878](https://arxiv.org/abs/2502.18878).
- Repo: `src/objdel/geometry/text.py` (`pattern`), `src/objdel/train/{configs,sample,rewards}.py`,
  `experiments/phase1/RESULTS.md` §7–8, `experiments/phase1/lab/t_stack/REPORT.md` (T5),
  `docs/research/notes/Phase 2 SFT/A1_models_tokens.md` (H-A1-8 / H-A1-9).
