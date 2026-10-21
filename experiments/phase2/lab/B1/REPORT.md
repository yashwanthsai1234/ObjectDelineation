# B1: model learnability screen (H-A1-1, H-A1-2, H-A1-8)

Wave B task B1 of `experiments/phase2/RESEARCH.md`. Mac M3, 16 GB, MPS. Run 2026-10-19 to 2026-10-20. Code: `experiments/phase2/lab/B1/`. Outputs: `<scratchpad>/phase2/B1/`.

Question: which small VLMs learn our polygon format (D10) fastest and best? The answer sets the Phase 2 screening ladder. A1 proposed LFM2.5-VL-450M → Qwen3.5-0.8B → Qwen3-VL-2B, with SmolVLM-256M demoted to an optional control.

## 0. Bottom line

1. **No arm uses the image after SFT on centred crops.**
   - On randomly shifted crops (coordinator request, A4 findings 1–3), every arm's median IoU_tc is 0.
   - Every arm is significantly below the image-blind prior on shifted crops (paired mean −0.026 to −0.045, every CI below 0).
   - 94–98 % of parsed samples land closer to the frame centre (median 2–3 px from it) than to the nodule (median 14.5 px away; the shift has median 15 px).
   - The models learned A4's centring shortcut. Crop-jitter training is mandatory (B6).
2. **The configured adapter destroys vision for Qwen3.5 and SmolVLM.**
   - `configs.lora` is LoRA r16 on all linear layers, vision tower and projector included, at lr 5e-4.
   - Cosine distance between images' projector outputs, base → after SFT:
     - Qwen3.5: 0.297 → 0.019 (36 steps) → 0.002 (150 steps);
     - SmolVLM: 0.171 → 0.003 → 0.0002;
     - LFM keeps its vision: 0.471 → 0.321 → 0.373.
   - The trained Qwen3.5 and SmolVLM give the same next-token distribution for the real crop, a black image and another item's crop (TV ≤ 0.0007, against 0.02–0.08 for the base models).
   - With the vision tower and projector frozen (LoRA on the language model only), Qwen3.5 keeps its features and its output depends on the image again (TV 0.008–0.015). Legality is unchanged (82 %) and centred IoU too (0.20 vs 0.19).
   - **Do not put LoRA on the vision tower at this learning rate.** Freeze it, or give it its own low learning rate (H-A2-6 becomes P1).
3. **Format learnability:**
   - At 36 steps × 8 rows, Qwen3.5 reaches 82–83 % legal samples at T = 1, against 55–59 % / 42–46 % for LFM (seeds 0 / 1) and 46–50 % for SmolVLM. All differences against Qwen3.5 are significant.
   - At 150 steps, LFM catches up on centred legality (83 % vs 85 %, not significant) and leads on centred median IoU_tc: 0.46 vs 0.40, paired difference +0.065 [0.026, 0.104].
   - SmolVLM is last at 150 steps (legal 67 %, IoU 0.28; significantly below both).
   - **No arm meets "learns the format"** (parse ≥ 95 %, legal ≥ 90 %) within 150 steps on the Mac.
4. **Zero-shot priors are weak for every model.**
   - Qwen3.5's native `bbox_2d` prior on shifted RIGA disc crops: median box IoU 0.41, which beats LFM by +0.13 [0.08, 0.18] but not the image-blind prior.
   - On LIDC no model localises (median ≤ 0.07).
   - Both Qwen3.5 and LFM answer in their own native convention even when the prompt asks for ours.
5. **Ladder recommendation:**
   - **LFM2.5-VL-450M → Qwen3.5-0.8B → Qwen3-VL-2B.** Screen both of the first two in the first Vast session; B1 cannot separate them at depth.
   - Drop SmolVLM-256M from the Vast screen.
   - Required changes: vision tower frozen or at a low learning rate, and crop-jitter training.

## 1. Pre-registration

Fixed before any measurement in `prereg.py`; `analyze.py` imports its constants. The amendments are listed in §1.5.

### 1.1 Common design

- **Arms:** SmolVLM-256M-Instruct, LFM2.5-VL-450M, Qwen3.5-0.8B (non-thinking, the model's default).
  - Qwen3.5 is not in `objdel.train.configs`. `common.ARMS` maps it to Qwen3-VL-2B's settings: processor defaults plus left padding, and the configured adapter.
  - TRL 1.16 swaps in its training templates (`smolvlm_training`, `qwen3_5_nothink_training`); LFM's own template has the markers.
  - Checked before the runs: for all three, the SFT text of a row starts with the generation prompt. For Qwen3.5 that prompt includes the empty `<think>\n\n</think>\n\n` block.
- **Eval items:**
  - 64 LIDC val items: a seeded random subset (seed 0) of the RL-eligible items (≥ 2 raters).
  - Test 1 also uses all 58 RL-eligible RIGA val items.
- **Views:**
  - centred: `targets.frame`.
  - shifted: the same S × S window moved by a uniform offset that keeps the raters' union inside (A4's `null_prior_jitter` procedure), one draw per item.
  - The image is cut exactly as `targets.frame` does. At offset (0, 0) the view is asserted equal to the frame and to the prepared PNG's pixel hash.
- **Prior reference row:** A4's image-blind centred disc, sized by the train median for the frame side S, scored with `reward.score_mask`.
- **Statistics:** 95 % percentile bootstrap over items, B = 2,000. Paired comparisons resample the same items for both arms; "significant" means the CI excludes 0.

### 1.2 Criteria

- **T1, H-A1-1 (zero-shot box prior):**
  - Procedure: greedy decoding, at most 96 new tokens. Native prompt per model card, plus a generic `{"bbox":[...]}` 0–1000 prompt. Strict parse.
  - ACCEPT on shifted RIGA disc crops with the native prompt: parse ≥ 90 %, median box IoU ≥ 0.6, and paired (model − prior) CI lower bound > 0.
  - REJECT: median < 0.4 on the shifted crops of both datasets.
  - Otherwise inconclusive.
- **T2, H-A1-2 (learnability), judged on the shifted view:**
  - Learns the format: parse ≥ 95 % and legal ≥ 90 %.
  - Uses the image: paired (item-mean IoU_tc − prior) CI lower bound > 0.
  - Ranking: by pooled median IoU_tc (illegal = 0).
  - H-A1-2 ACCEPT iff (i) LFM is not significantly below Qwen3.5, and (ii) Qwen3.5 is significantly above SmolVLM, or SmolVLM fails the format criterion while Qwen3.5 passes it. REJECT iff contradicted significantly. Otherwise inconclusive.
  - A1's rule: drop SmolVLM if it is more than 0.05 below the best arm.
- **T3, H-A1-8 (sampling legality):**
  - Pass: legal ≥ 80 % and duplicates ≤ 10 % on both views.
  - LFM zero-led tokens are a material failure mode iff ≥ 1 % of samples fail with `leading_zero`.

### 1.3 Training (identical across arms except as amended)

- TRL `SFTTrainer` via `configs.sft_config`, rows from `dataset.sft_rows` (every voting rater's target). That gives 3,242 rows from the same 1,000 LIDC train items (seed 0).
- Adapter: `configs.lora` (r16, α 32, all linear layers).
- lr 5e-4, constant after 10 warmup steps; 8 rows per step; fp32; MPS capped at 7.0 GB.

### 1.4 Evaluation

- K = 8 samples per item at T = 1 (top_p 1, top_k 0, no penalties), through TRL GRPOTrainer's HF generation path: `common.sample_k`, which is `objdel.train.sample.sample`'s loop plus the token ids. No reference model is loaded (β = 0).
- Scored with `objdel.eval.metrics` / `summarize`, which use the Phase 2 gate definitions.

### 1.5 Amendments

| # | When | Change | Reason |
|---|---|---|---|
| 1 | After the timing pilots, before any eval | Micro-batch 1 × 8 accumulation; N = 36 steps for every arm; max_completion 1,024 | Micro-batch 2 OOMs at the 6.5 GiB cap (Qwen3.5 even with gradient checkpointing). Qwen3.5 runs at 79 s/step in fp32. The longest LIDC train target is 881 tokens. |
| 2 | After the 36-step SmolVLM and LFM evals | Exploratory 150-step tier with the same criteria; seed 1 at 36 steps for LFM and SmolVLM | No arm learned the format at 36 steps. |
| 3 | Before any Qwen3.5 eval | Qwen3.5 frozen base weights in bf16 (LoRA fp32) | The fp32 run hit the cap at step 13. |
| 4 | Before any Qwen3.5 eval | Qwen3.5 added to the 150-step tier (about 98 min, over the 45-min guideline) | bf16 made it 34 s/step. |
| 5 | Post-hoc, exploratory | LoRA on the language model only (vision tower and projector frozen) for Qwen3.5 and SmolVLM at 36 steps | Motivated by the vision-collapse diagnostic. No verdict attached. |

## 2. Results

### 2.1 Test 1: zero-shot box prior

Greedy decoding. Cells are median box IoU [95 % CI], parse rate 100 % unless noted. The prior row is the box of A4's centred disc.

| | RIGA disc, centred | RIGA disc, shifted | LIDC, centred | LIDC, shifted |
|---|---|---|---|---|
| prior (image-blind) | 0.855 | 0.388 | 0.624 | 0.000 |
| Qwen3.5-0.8B, native `bbox_2d` | 0.403 [0.369, 0.425] | **0.414 [0.377, 0.435]** | 0.140 [0.088, 0.186] | 0.072 [0.045, 0.118] |
| LFM2.5-VL-450M, native `bbox` [0, 1] | 0.272 [0.245, 0.302] | 0.315 [0.269, 0.353] | 0.060 | 0.053 |
| SmolVLM-256M, generic (its only option) | 0, parse 0 % | 0, parse 0 % | 0, parse 0 % | 0, parse 0 % |

- **Paired against the prior, shifted crops:**
  - Qwen RIGA +0.033 [−0.029, +0.098]; LFM RIGA −0.095 [−0.150, −0.042].
  - LIDC: Qwen +0.052 [0.005, 0.094], LFM +0.043 [−0.003, +0.090]. Here the prior scores 0, so any box near the nodule wins.
- **Paired Qwen − LFM (post-hoc, `zs_posthoc.py`):** RIGA shifted +0.128 [0.078, 0.177]; LIDC shifted +0.009 [−0.026, +0.046].
- **Generic 0–1000 prompt:** both prior-bearing models ignore the instruction and answer in their native convention. Qwen writes `bbox_2d` in a ```json fence; LFM writes [0, 1] floats.
  - Rescored leniently (post-hoc), the generic prompt gives Qwen RIGA 0.39–0.40 and LFM 0.06–0.08.
  - SmolVLM writes placeholder JSON (`[x1,y1,x2,y2]`).
- **D10 polygon prompt, zero-shot:** 0 % parse for all three (16 items per dataset). SmolVLM and Qwen produce multi-line JSON; LFM produces runaway `[0,0],[0,1],…` rings.
- **Verdict H-A1-1: INCONCLUSIVE.** Qwen3.5's shifted RIGA box IoU is 0.41 < 0.6, so it does not pass. It is not < 0.4 on both datasets, so it is not rejected. The prior that A1 hoped for is weak: real against LFM on RIGA, absent on LIDC, and never better than the image-blind prior.

### 2.2 Test 2: learnability, 36 steps (primary budget, seed 0)

64 LIDC val items × K = 8. Cells are point estimates [95 % CI over items]. Medians are pooled over samples, with illegal samples counted as 0.

| arm | view | parse | legal | median IoU_tc | mean IoU_tc | IoU ≥ 0.5 | mean (model − prior) | dup | σ > floor | tokens/answer |
|---|---|---|---|---|---|---|---|---|---|---|
| prior | centred / shifted | — | — | 0.592 / 0.000 | — | — | — | — | — | — |
| SmolVLM-256M | centred | 0.959 | 0.502 [0.459, 0.543] | 0.000 [0, 0.056] | 0.170 | 0.096 | −0.420 | 0.004 | 1.00 | 122 |
| | shifted | 0.967 | 0.459 [0.416, 0.500] | 0 | 0.022 | 0.002 | −0.045 [−0.072, −0.022] | 0.004 | 1.00 | 124 |
| LFM2.5-VL-450M | centred | 0.951 | 0.588 [0.539, 0.639] | 0.133 [0.056, 0.185] | 0.187 | 0.086 | −0.403 | 0.006 | 1.00 | 81 |
| | shifted | 0.940 | 0.555 [0.522, 0.590] | 0 | 0.029 | 0.002 | −0.038 [−0.064, −0.015] | 0.010 | 1.00 | 80 |
| Qwen3.5-0.8B | centred | 0.971 | **0.830 [0.795, 0.863]** | 0.190 [0.155, 0.225] | 0.230 | 0.135 | −0.360 | 0.004 | 1.00 | 81 |
| | shifted | 0.969 | **0.822 [0.785, 0.857]** | 0 | 0.030 | 0.000 | −0.036 [−0.063, −0.015] | 0.006 | 0.98 | 79 |

**Paired differences, shifted view:**
- Legal:
  - SmolVLM − LFM −0.096 [−0.150, −0.043];
  - SmolVLM − Qwen3.5 −0.363 [−0.420, −0.309];
  - LFM − Qwen3.5 −0.268 [−0.315, −0.219].
- The pooled median IoU_tc difference is 0 [0, 0] for every pair, because all medians are 0.

**Centred view:** median IoU_tc LFM − Qwen3.5 −0.057 [−0.142, +0.005]; mean IoU_tc LFM − Qwen3.5 −0.043 [−0.072, −0.016].

**Seed 1 (LFM and SmolVLM):**

| arm | legal centred / shifted, seed 1 | seed 0 for comparison |
|---|---|---|
| LFM | 0.424 / 0.463 | 0.588 / 0.555 |
| SmolVLM | 0.480 / 0.453 | 0.502 / 0.459 |

- LFM's seed-to-seed change (−16 / −9 points) exceeds its item-bootstrap CI.
- At seed 1, SmolVLM is not below LFM: the legal difference is +0.057 [−0.008, +0.123] centred.
- So **the 36-step LFM vs SmolVLM order is not robust**. Qwen3.5's lead (+24 to +36 points) is larger than any observed seed swing, but Qwen3.5 has a single seed.

### 2.3 Exploratory 150-step tier (seed 0)

| arm | view | parse | legal | median IoU_tc | mean IoU_tc | IoU ≥ 0.5 | mean (model − prior) | dup | tokens/answer |
|---|---|---|---|---|---|---|---|---|---|
| SmolVLM-256M | centred | 0.982 | 0.668 [0.631, 0.705] | 0.283 [0.239, 0.314] | 0.270 | 0.217 | −0.320 | 0.000 | 119 |
| | shifted | 0.982 | 0.668 | 0 | 0.041 | 0.008 | −0.026 [−0.048, −0.008] | 0.000 | 119 |
| LFM2.5-VL-450M | centred | **0.994** | 0.832 [0.799, 0.865] | **0.461 [0.432, 0.497]** | 0.414 | 0.440 | −0.176 | 0.000 | 75 |
| | shifted | 0.986 | 0.709 [0.660, 0.758] | 0 | 0.039 | 0.012 | −0.027 [−0.045, −0.012] | 0.000 | 91 |
| Qwen3.5-0.8B | centred | 0.971 | **0.852 [0.820, 0.885]** | 0.396 [0.376, 0.431] | 0.369 | 0.322 | −0.221 | 0.002 | 94 |
| | shifted | 0.971 | **0.834 [0.807, 0.861]** | 0 | 0.041 | 0.008 | −0.026 [−0.044, −0.011] | 0.002 | 95 |

- **Centred, paired:**
  - Median IoU_tc: LFM − Qwen3.5 +0.065 [0.026, 0.104]; SmolVLM − LFM −0.178 [−0.239, −0.142]; SmolVLM − Qwen3.5 −0.114 [−0.168, −0.081].
  - Legal: LFM − Qwen3.5 −0.020 [−0.065, +0.025].
- **Shifted legal:** LFM − Qwen3.5 −0.125 [−0.178, −0.076]. LFM's legality degrades on out-of-distribution crops; Qwen3.5's does not.
- **SmolVLM's centred and shifted samples are identical** (64 of 64 items, same seed). Its output no longer depends on the image at all (§2.5).

### 2.4 Exploratory: LoRA on the language model only (36 steps, AMEND 5)

| arm | view | parse | legal | median IoU_tc | mean IoU_tc | tokens/answer | s/step |
|---|---|---|---|---|---|---|---|
| Qwen3.5-0.8B, LM-only | centred | 0.973 | 0.822 [0.787, 0.857] | 0.200 [0.160, 0.233] | 0.228 | 76 | 31.7 |
| | shifted | 0.973 | 0.820 | 0 | 0.031 | 79 | |
| SmolVLM-256M, LM-only | centred | 0.955 | 0.451 [0.408, 0.492] | 0.000 | 0.150 | 133 | 8.7 |
| | shifted | 0.939 | 0.471 | 0 | 0.021 | 126 | |

- Same metrics as the all-linear adapter at 36 steps.
- The centre bias stays (95 % of samples closer to the frame centre on shifted crops), because the training crops are still centred.

### 2.5 Diagnostics (post-hoc): does the model look at the image?

**`vision_feats.py`** compares the projector output (what the language model sees) for 4 crops and a black image. **`image_use.py`** measures the next-token TV between the centred crop and a shifted crop, a black image or another item's crop, teacher-forced on the first 24 target tokens of 8 items. Syntax tokens dilute the TV, so compare it across rows, not as an absolute.

| checkpoint | cosine distance between images | feature norm | TV vs shifted | TV vs black | TV vs other item |
|---|---|---|---|---|---|
| Qwen3.5 base | 0.297 | 4.5 | 0.036 | 0.076 | 0.042 |
| Qwen3.5, 36 steps, all-linear | 0.019 | 1.6 | 0.0005 | 0.0004 | 0.0007 |
| Qwen3.5, 150 steps, all-linear | 0.002 | 16.8 | 0.0000 | 0.0000 | 0.0000 |
| Qwen3.5, 36 steps, LM-only | 0.297 (frozen) | 4.5 | 0.008 | 0.015 | 0.009 |
| SmolVLM base | 0.171 | 110 | 0.020 | 0.026 | 0.025 |
| SmolVLM, 36 steps, all-linear | 0.003 | 183 | 0.0003 | 0.0004 | 0.0004 |
| SmolVLM, 150 steps, all-linear | 0.0002 | 155 | 0.0000 | 0.0001 | 0.0001 |
| SmolVLM, 36 steps, LM-only | (frozen) | — | 0.003 | 0.006 | 0.004 |
| LFM base | 0.471 | 130 | 0.006 | 0.023 | 0.011 |
| LFM, 36 steps, all-linear | 0.321 | 310 | 0.004 | 0.006 | 0.004 |
| LFM, 150 steps, all-linear | 0.373 | 289 | 0.034 | 0.031 | 0.034 |

**Centre bias (`centre_bias.py`):** share of parsed samples closer to the frame centre than to the nodule.

| | centred view | shifted view |
|---|---|---|
| All arms, both depths | 0.41–0.65 (the two points nearly coincide) | 0.94–0.98 |
| Median distance to frame centre / nodule centroid, shifted view | — | 2.1–3.4 px / 14.4–14.8 px (median shift 15 px) |

**Vertices (`vertex_stats.py`):** median per parsed sample.

| arm | vertices (36 / 150 steps) | tokens per vertex |
|---|---|---|
| SmolVLM | 11 / 12 | 9.3 |
| LFM | 10 / 10 | 6.5 |
| Qwen3.5 | 8 / 8 | 9.5 |
| targets | 9 | — |

This matches A1's 1.5× tokenizer gap. Qwen3.5 offsets it by drawing fewer vertices, so its answers are about as long in tokens as LFM's.

### 2.6 Test 3: H-A1-8, legality at T = 1

| checkpoint | legal centred / shifted | dup | verdict (legal ≥ 80 %, dup ≤ 10 %, both views) |
|---|---|---|---|
| Qwen3.5, 36 steps | 0.830 / 0.822 | ≤ 0.006 | **pass** |
| Qwen3.5, 150 steps | 0.852 / 0.834 | 0.002 | **pass** |
| Qwen3.5, 36 steps, LM-only | 0.822 / 0.820 | 0.002 | **pass** |
| LFM, 36 steps (seed 0 / 1) | 0.588 / 0.555; 0.424 / 0.463 | ≤ 0.02 | fail |
| LFM, 150 steps | 0.832 / 0.709 | 0.000 | fail (centred passes) |
| SmolVLM, all depths | 0.45–0.67 | ≤ 0.014 | fail |

- **LFM zero-led number tokens** (e.g. `042`): 0 to 0.4 % of samples contain one; 0 to 0.2 % fail with `leading_zero`. Not material.
- **Truncation** at 1,024 tokens: ≤ 0.2 %.
- **σ above the floor:** 97–100 % of groups. Duplicates are not a GRPO risk at these depths.

### 2.7 Cost on the Mac (MPS)

Training is 8 rows per step. "Gen s/item" is K = 8 samples for one item. The Mac's load average varied from 2 to 25 during the runs, so times are indicative.

| arm | s/step (36-step run / 150-step run) | peak MPS memory, training | peak MPS memory, eval | gen s/item |
|---|---|---|---|---|
| SmolVLM-256M, fp32 | 14.4 / 10.5 | 3.9 / 5.0 GB | 2.8–3.9 GB | 6.6–11.3 |
| LFM2.5-VL-450M, fp32 | 15.7 / 12.3 | 3.8 GB | 3.6 GB | 6.6–10.6 |
| Qwen3.5-0.8B, bf16 base | 34.4 / 37.7 | 4.9 / 6.0 GB | 3.6 GB | 14.7–24.4 |
| Qwen3.5-0.8B, fp32 | 79 (pilot) | OOM at 6.5 GiB, step 13 | — | — |

Profiling (`prof_step.py`, transformers source):
- On macOS Qwen3.5 runs its Gated-DeltaNet and causal-conv1d layers on the reference PyTorch implementations. `fla` is Triton-only (NVIDIA, AMD and Intel), and transformers warns "correct but much slower".
- A no-grad fp32 forward of one 329-token LIDC row takes 3.5 s for Qwen3.5, against about 0.8 s forward+backward for LFM.
- **Qwen3.5's Mac cost says nothing about Vast**: install `flash-linear-attention` and `causal-conv1d` there.

## 3. Verdicts

| hypothesis | verdict | basis |
|---|---|---|
| H-A1-1 (Qwen3.5 has a usable zero-shot prior) | **INCONCLUSIVE** | Shifted RIGA disc 0.41 < 0.6, not above the image-blind prior. Beats LFM's prior by +0.13 on RIGA. No LIDC prior for any model. |
| H-A1-2 (LFM ≥ Qwen3.5 > SmolVLM at equal steps) | **INCONCLUSIVE** (pre-registered rule) | All shifted medians are 0 and no arm learns the format, so the ranking statistic ties. Secondary: legality Qwen3.5 ≫ LFM ≈ SmolVLM at 36 steps (LFM vs SmolVLM flips with the seed). At 150 steps, LFM ≈ Qwen3.5 > SmolVLM on legality, and LFM > Qwen3.5 > SmolVLM on centred IoU. |
| A1's SmolVLM drop rule | not triggered as written (shifted medians all 0) | On centred crops SmolVLM is 0.18 below LFM at 150 steps, its vision collapses, and it has no prior. **Drop it from the Vast screen.** |
| "Uses the image" (coordinator) | **FAIL for every arm** | Centring shortcut, plus vision collapse under the all-linear adapter (Qwen3.5, SmolVLM). |
| H-A1-8 (legal sampling without constrained decoding) | **PASS for Qwen3.5** (36 and 150 steps). LFM: pass at 150 steps on centred crops, fail on shifted. SmolVLM: fail. | Duplicates ≤ 2 % everywhere. LFM zero-led tokens are not material. |

## 4. What this changes in the plan

1. **Ladder:**
   - LFM2.5-VL-450M → Qwen3.5-0.8B → Qwen3-VL-2B. Screen LFM and Qwen3.5 side by side in the first Vast session.
   - LFM is the smallest; it is cheapest per token (6.5 vs 9.5 tokens per vertex), keeps its vision under LoRA, and has the best centred IoU at 150 steps.
   - Qwen3.5 learns legal syntax 4× faster in steps, holds legality on shifted crops, and has the better RIGA box prior.
   - SmolVLM is dropped from the Vast screen. It remains, at most, a free control.
2. **Vision tower:**
   - `configs.lora` (all-linear, vision included, lr 5e-4) erases image information for Qwen3.5 and SmolVLM within 36 steps.
   - For LoRA: exclude the vision tower and projector (`exclude_modules`, as `common.lora("lm")`), or give them a separate learning rate ≤ 0.1×.
   - For full fine-tuning (configs: SmolVLM, LFM): test the vision learning rate before trusting it (H-A2-6, now P1).
   - Add a "vision intact" check to every SFT checkpoint: `vision_feats.py`'s cosine distance and `image_use.py`'s TV.
3. **Crop jitter is mandatory in SFT** (B6), and the eval protocol must include shifted crops. Without jitter every arm draws at the frame centre.
4. **Constrained decoding:**
   - Not needed for Qwen3.5 at GRPO start (legal ≥ 82 % from 36 steps).
   - For LFM, decide after the Vast SFT: 71–83 % legal at 150 Mac steps.
5. **Per-model settings:**

| model | settings |
|---|---|
| Qwen3.5-0.8B | Non-thinking template: TRL `qwen3_5_nothink_training`, prefix-consistent with the generation prompt. bf16. Install `flash-linear-attention` and `causal-conv1d`. Sample with T = 1, top_p 1, top_k 0, no presence penalty; this overrides the card's T 0.7 / top_p 0.8 / top_k 20 / presence 1.5. max_completion 1,024 is enough for LIDC (longest train target 876 tokens). |
| LFM2.5-VL-450M | Native template (has generation markers). Override the card's text sampling defaults (T 0.1, min_p 0.15, repetition 1.05) to T = 1 for gates and GRPO. The [0, 1]-float prior is replaced by our 0–1000 format within 36 steps (parse 94–99 %). max_completion 1,024 (longest LIDC target 559). |
| SmolVLM-256M (if kept as a control) | `longest_edge` 512 for LIDC (64 visual tokens). Leading whitespace in completions is stripped by the parser. |

## 5. Limits

- **Not the Vast recipe:**
  - MPS LoRA r16 at 36 or 150 steps, not full fine-tuning for epochs.
  - One learning rate for all arms (5e-4 constant).
  - Qwen3.5 has bf16 base weights while the other arms are fp32 (AMEND 3).
  - LIDC only; 64 val items.
- **Seeds:**
  - Seed 1 exists only for LFM and SmolVLM at 36 steps. It shows run-to-run swings (up to 16 legal points) larger than the item-bootstrap CIs.
  - Qwen3.5 and every 150-step result are single-seed.
- **Post-hoc material:**
  - The 150-step tier and the LM-only arms were added after seeing results (AMEND 2, 5); they are exploratory.
  - The image-use diagnostics are post-hoc and small: 4 to 8 items. The TV is diluted by syntax tokens.
- **The test differs from A1's H-A1-2 design:** LIDC only, all raters, no RIGA rows, K = 8 at T = 1 instead of greedy decoding.
- **Timing:** Mac load varied (load average 2–25), and one Qwen3.5 eval overlapped about 10 s with an accidental second process. Treat s/step as indicative.

## 6. Reproduce

Run from the repo root. Scripts read untrusted `data/`, hence `-I`.

```bash
HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B1.zeroshot --arm qwen3.5-0.8b   # and the other two arms
HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B1.train_eval --arm lfm2.5-vl-450m --steps 36 --seed 0
HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B1.train_eval --arm qwen3.5-0.8b --steps 36 --base-dtype bfloat16 [--lora-scope lm]
.venv/bin/python -I -m experiments.phase2.lab.B1.analyze
.venv/bin/python -I -m experiments.phase2.lab.B1.{zs_posthoc,centre_bias,vertex_stats}
HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase2.lab.B1.{vision_feats,image_use} --arm <arm> --run <run|base>
```

Kept checkpoints: `<scratch>/phase2/B1/checkpoints.json`.

## 7. Sources

- Qwen3.5-0.8B model card (non-thinking default; recommended sampling; grounding): https://huggingface.co/Qwen/Qwen3.5-0.8B
- LFM2.5-VL-450M model card (grounding prompt with [0, 1] coordinates; recommended sampling; image token limits): https://huggingface.co/LiquidAI/LFM2.5-VL-450M
- Qwen3-VL grounding convention (`bbox_2d`, relative 0–1000): https://github.com/QwenLM/Qwen3-VL/blob/main/cookbooks/2d_grounding.ipynb ; A1 §2–3
- flash-linear-attention (Triton kernels; NVIDIA, AMD and Intel): https://github.com/fla-org/flash-linear-attention
- ms-swift Qwen3.5 best practice (fla / causal-conv1d for training speed; `enable_thinking false`): https://swift.readthedocs.io/en/v4.2/BestPractices/Qwen3_5-Best-Practice.html
- Unsloth Qwen3.5 fine-tuning (dependencies; 0.8B LoRA memory): https://unsloth.ai/docs/models/qwen3.5/fine-tune
- transformers 5.20.1 `models/qwen3_5/modeling_qwen3_5.py` (reference PyTorch `torch_chunk_gated_delta_rule` fallback), local source
- TRL 1.16 SFTTrainer (training templates, assistant-only loss): https://huggingface.co/docs/trl/main/en/sft_trainer ; GRPOTrainer HF generation path: https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/grpo_trainer.py
- PEFT LoRA (`target_modules="all-linear"`, `exclude_modules`, `merge_and_unload`): https://huggingface.co/docs/peft/main/en/package_reference/lora
- Image-LoRA, selective adaptation of visual tokens and heads: https://arxiv.org/pdf/2512.19219
- Vision tower frozen vs unfrozen: A2 §4.2 (Qwen fine-tuning README `tune_mm_vision`; Axolotl frozen vision for LoRA; VPP; VectorLLM)
- torch.mps memory API: https://docs.pytorch.org/docs/stable/mps.html
- Repo: A1 (`A1_models_tokens.md` §1.1, §7), A2 §4, A4 §0.1 and `lab/A4/null_prior*.py`, A5 §1 (gate sampling settings), `experiments/phase1/RESULTS.md` §8.
