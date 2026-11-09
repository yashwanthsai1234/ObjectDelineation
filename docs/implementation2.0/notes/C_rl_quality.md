# C. Why the GRPO-tuned Qwen3-VL-2B stops at IoU_tc 0.645, and how to push it further (2026-11-08)

Scope: read only. I changed no repo files, ran no git commands, rented no GPU. Python that read `data/` ran with `-I`.
Scripts and their outputs are in this scratchpad (`ana/join.py`, `ana/stats.py`, `ana/rows.json`).
Tags: **MEASURED** = I computed it from repo files or run logs. **CODE** = I read it in the repo or in the installed
TRL 1.17.1. **LIT** = a published source. **EST** = my estimate; nobody has measured it.
I treated the earlier notes (plateau diagnosis, full RL, final eval) as claims to test. I did not read session
transcripts, the memory folder or `.claude/`.

## 0. Skills used

| Skill | What it added |
|---|---|
| `engineering:debug` | The structure of this file: reproduce, isolate, diagnose, fix, prevent (sections 2 to 6). |
| `engineering:code-review` | A close read of `reward/advantages.py`, `train/fast.py` (`is_ratio`, C10), `cuda/gpucheck.cpp` (IoU_tc) and TRL's IS code. It turned up two findings: the IS masking issue (3.2) and a wrong definition of IoU_tc in the brief (3.8). |
| `engineering:system-design` | Splitting train-time from test-time verification. The checker needs rater masks, so it cannot verify at test time (3.6). Every fix also gets a cost/benefit line. |
| `engineering:architecture` | Each representation change is weighed as an ADR-style option (keep the polygon text, add a box-first prefix, Text4Seg-style descriptors, VLM + SAM hybrid) in 3.7 and section 4. |
| `engineering:documentation` | The bottom line comes first. Every number is tagged. Links point to repo files and are not copied. |
| `huggingface-skills:huggingface-best` | The main thread's live leaderboard lookups (`MAIN_skill_lookups.md`, corrected version). On the card, RefCOCO avg is 88.2 for Qwen3-VL-4B and 84.8 for both Qwen3.5-2B and our Qwen3-VL-2B. A same-size swap to Qwen3.5-2B is therefore not a grounding upgrade. The upgrade candidate is the 4B class (Qwen3-VL-4B / Qwen3.5-4B), in experiment E8. |
| `huggingface-skills:huggingface-papers` | Read arXiv 2606.13156 (Iterative Visual Thinking). Visual self-correction gains vanish under label-free stopping rules, which lowers the rank of multi-turn refinement. |
| `nvidia-skills:nvidia-skill-finder` | The main thread checked the live catalog. No catalog skill covers GRPO training of a VLM. `nv-segment-ct` / `nv-segment-ctmr` (VISTA-3D-style CT segmentation) could serve as a non-VLM reference baseline. Nothing was installed. |
| WebSearch / Exa | Literature in section 4 (URLs inline). |

## 1. Bottom line

1. **Most of the gap is detection, not drawing (MEASURED).** On the val split (shifted view, 1,084 items), items
   with item-mean IoU_tc < 0.1 ("misses", 12.6 % of items) cause **76 % of the gap** to the doctor ceiling. Items
   the model finds score **0.706 against a ceiling of 0.734**. Misses sit off centre: 31.8 % of items with offset
   ≥ 0.3 are misses; under 0.15 it is 0.4 %.
2. **RL itself created the offset penalty (MEASURED, new evidence).** On the 128 probe items, the SFT checkpoint e2
   scores the same at every offset (0.49 / 0.49 / 0.54 / 0.51). After RL: +0.24 on centred items, only +0.05 at
   offset ≥ 0.3, and the miss rate there rose from 8 % to 12 % (n = 24, noisy). RL rows really are centred: offset
   median 0.012, 0.0 % above 0.2, against 60.5 % above 0.2 on val. The fix already built (`train_rl --train-views`)
   targets the right cause.
3. **Trainer and vLLM disagree, and the IS correction throws away gradient (MEASURED + CODE, new).** TRL's default
   `vllm_importance_sampling_mode="sequence_mask"` with cap 3 is active. The mean sequence IS ratio is 0.70-0.86
   after step 10 (0.99 at the start). The minimum is 0 in **98 % of steps** and the maximum is about 2.9 in every
   step. With one update per batch the ratio should be exactly 1, so all of this is numeric mismatch: bf16, the
   merged LoRA in vLLM against the unmerged LoRA in the trainer, different kernels. EST: the effective batch is
   about 40-55 % of the 256 rollouts.
4. **The polygon representation is not the ceiling (MEASURED).** One pixel is 20 grid units (W = 50) or 10
   (W = 100), so quantization error is 0. Targets average 14.2 vertices and 177 characters. Truncation is 0.12 %.
5. **The fastest owner-visible wins need no training (EST).** At test time, keep only legal samples (the legality
   checker needs no labels, so legal → ~100 %). Then pick the medoid of K samples (MBR) or fuse the K masks by pixel
   majority. The checker **cannot** be used as a verifier at test time: IoU and F5 need the rater masks.
6. **A credible "top of field" claim needs a fair ceiling (CODE).** The model's IoU_tc is the best of all k raters.
   The doctor ceiling is the best of the other k − 1. The comparison favours the model, so state the metric before
   making claims.

## 2. Reproduction: the facts I checked

| Fact | Value | Source |
|---|---|---|
| Test (2,200) IoU_tc / F5 / P / R / Dice / legal | 0.645 / 0.816 / 0.695 / 0.748 / 0.701 / 99.4 % | `runs/r5/eval/eval_test/final_eval.json` (MEASURED: read) |
| Val (1,084) IoU_tc, LOO ceiling (plain IoU), F5 | 0.618, 0.723, 0.806 | `ana/stats.py` on `eval_val/items.jsonl` + view targets |
| RL train items / crop sides | 7,341 RL-eligible; W = 50: 6,267, W = 100: 1,068, W = 125: 6 | `data/prepared/lidc/train.jsonl` |
| Object offset (Chebyshev, fraction of side): train / val | median 0.012 / 0.230; share > 0.2: 0.0 % / 60.5 %; > 0.3: 0.0 % / 22.0 % | MEASURED |
| Object area, px (median): train / val | 106 / 85; share filling < 5 % of the crop: 59.5 % / 68.1 % | MEASURED |
| Target vertices: mean / median / p95 / max | train 14.2 / 11 / 32 / 145; val 12.7 / 10 / 28 / 79 | MEASURED |
| Multi-part / holed items | train 1.6 % / 3.7 %; val 0.6 % / 4.0 % | MEASURED |
| Training IoU_tc, steps 1-10 → 700-770 | 0.523 → 0.771; F5 0.717 → 0.943; g_area 0.737 → 0.944 | `runs/rl/train_rl/rl.json` |
| Entropy | 0.361 (1-10) → 0.241 (10-50) → 0.066 (700-770) | rl.json |
| Advantage sd | 1.44 → 0.94 → 0.61 → 0.39 → 0.33 | rl.json |
| KL (logged mean) and beta | 0.29-0.33 after step 50; beta 1e-3 (`configs.py:188`), so the KL term is about 3e-4 of the loss | rl.json, CODE |
| KL / loss spikes | step 156 (KL 1,020), 342 (3,817), 407 (730), 746 (398,284; loss 28.3); grad norm unchanged | rl.json |
| IS ratio (sequence_mask, cap 3) | mean 0.99 (1-10), 0.70 (100-199), 0.80-0.85 after; min = 0 in 98 % of steps; max ≈ 2.9 | rl.json |
| Completion length | 105 tokens mean; clipped 0-0.2 %; truncation errors 0.12 % | rl.json |
| PPO clip active | never (`clip_ratio/*` = 0): one update per generation batch | rl.json |

The plateau note's main claims hold under my own recomputation. Train offsets: 0.012 vs its 0.018. Val share above
0.2: 60.5 % vs its 68 %. I measure the Chebyshev offset; it seems to use a Euclidean one. The miss concentration off
centre holds too. Its entropy and advantage-sd trends also hold.

## 3. Root-cause analysis

### 3.1 Primary: off-centre detection failures that RL never trained on (confirmed)

Val, shifted view (MEASURED, `ana/stats.py`). In the table, gap = LOO ceiling − model, summed over items.

| Offset | n | IoU_tc | ceiling | F5 | miss | share of gap | IoU_tc of found items |
|---|---|---|---|---|---|---|---|
| < 0.15 | 255 | 0.744 | 0.732 | 0.930 | 0.4 % | −3 % | 0.747 |
| 0.15-0.2 | 173 | 0.673 | 0.719 | 0.883 | 5.2 % | 7 % | 0.707 |
| 0.2-0.3 | 417 | 0.603 | 0.727 | 0.799 | 12.2 % | 46 % | 0.685 |
| ≥ 0.3 | 239 | 0.472 | 0.710 | 0.633 | 31.8 % | 50 % | 0.687 |

- Found items score 0.69-0.75 at **every** offset. The offset only changes how often the model misses.
- Small objects miss more often: 20-40 px, 33 % misses; < 1 % of the crop area, 40 % misses. Area matters less once
  the object is found (0.665-0.725).
- **SFT vs RL by offset** (128 probe items, `runs/rl/eval/{ckpt,final}/heldout.json` joined with my offsets):

  | Offset | n | e2 (SFT) | final (RL) | gain | miss e2 → final |
  |---|---|---|---|---|---|
  | < 0.15 | 35 | 0.488 | 0.728 | +0.240 | 0 → 0 % |
  | 0.15-0.2 | 20 | 0.488 | 0.644 | +0.156 | 10 → 5 % |
  | 0.2-0.3 | 49 | 0.536 | 0.665 | +0.128 | 4 → 6 % |
  | ≥ 0.3 | 24 | 0.513 | 0.559 | +0.046 | 8 → 12 % |

  SFT trained on shifted and rescaled windows and has no offset gradient. RL trained only on centred frames
  (`dataset.grpo_rows`, `train_rl.py:51` before `--train-views`), so it improved centred items the most. It may have
  learned a centre prior. This is direct evidence, beyond correlation, that the training distribution causes the
  offset penalty.
- Upper bounds (EST, val): found items average 0.706. If every miss became an average found item, IoU_tc would be
  about 0.70 (+0.085). If every offset bucket matched the centred one, it would be about 0.744 (+0.126).

### 3.2 Secondary, and cheap to fix: train-inference mismatch and sequence-level IS masking

- CODE: with C10 there is one update per batch, so π_old = π_θ and the IS ratio is pure numeric mismatch
  (`fast.py:530-538, 581-597`). TRL 1.17.1 defaults are `vllm_importance_sampling_mode="sequence_mask"` and
  `clip_max=3.0` (`trl/trainer/grpo_config.py:330-347`). The ratio is exp(Σ_t (logp_train − logp_vllm)) per
  sequence. A ratio above 3 sets that sequence to 0. A ratio below 1 is kept but **down-weights** the sequence,
  because there is no lower clip (`fast.py:351-363`; TRL `grpo_trainer.py:2784-2818`).
- MEASURED: the mean ratio is 0.99 in steps 1-10, when the LoRA delta is 0, then 0.70-0.86 as the adapter grows.
  The minimum is 0 in 98 % of steps and the maximum sits just under the cap. A growing adapter delta fits the vLLM
  side: there the LoRA is merged into bf16 weights, W + BA rounded, while the trainer computes Wx + BAx unmerged.
  The bf16 merge/unmerge drift in `full_rl_2026-11-08.md` adds to it. The single-token KL spikes (up to 4e5) are
  tokens that vLLM sampled and the trainer gives about 0 probability.
- EST: if the log sequence ratio is roughly normal with E[r] = 1, a masked mean of 0.73-0.85 means σ ≈ 0.8-1.0 nats.
  About 4-6 % of sequences are dropped, but they carry 15-27 % of the IS weight. The effective sample size falls to
  about exp(−σ²) ≈ 37-53 % of 256. Two things follow. The gradient is noisier. And the signal is biased against the
  sequences the trainer likes more than vLLM does, which are those it has just learned to prefer.
- LIT: the same failure is reported with fixes. FP16 rollouts and training remove most of it
  ([Liu et al. 2510.26788](https://arxiv.org/abs/2510.26788); about 24x lower token KL than BF16). Bitwise-consistent
  vLLM + TorchTitan trains in fewer steps to a higher reward
  ([vLLM blog](https://vllm.ai/blog/2025-11-10-bitwise-consistent-train-inference)). Sequence-level likelihood ratios
  with length normalization are more tolerant of precision gaps ([GSPO 2507.18071](https://arxiv.org/abs/2507.18071)).
- Verdict: this did not cause the plateau, because centred training items saturate anyway. It does make every
  future step worse, and the hard off-centre windows of E1 need clean gradients.

### 3.3 Consequence, not cause: the in-distribution reward saturates and entropy collapses

- Training F5 0.943 and g_area 0.944 are near their tops. F5 is calibrated so that a new doctor scores 0.95
  (`spec.py:36`). Inside the egg-yolk band F5 is flat, so it cannot reward boundary moves smaller than the doctors'
  own spread. Per-group sds fall below the floors (0.05 for IoU, F5, g_area; `spec.py:29`), so z-scores shrink
  (`advantages.py:100`). This explains why the advantage sd falls 1.44 → 0.33 without any dead groups
  (`frac_reward_zero_std` ≤ 3 %).
- Entropy 0.36 → 0.066 per token. LIT: performance trades off against entropy, R = −a·exp(H) + b
  ([2505.22617](https://arxiv.org/abs/2505.22617)). On centred crops the model already beats the LOO ceiling
  (0.744 vs 0.732), so low entropy hurts only where exploration is needed: the off-centre items it never saw.
- Clip-higher (`epsilon_high` 0.28, `configs.py:188`) has **no effect** because the ratio is always 1 in the update.
  None of the DAPO/GSPO clip mechanics act here.

### 3.4 The reward gives no gradient for "all 8 missed"

- When every rollout misses, IoU_tc = 0 for all of them and r_topo = 0, because validity needs IoU ≥ 0.5
  (`advantages.py:73-76`). F5 still varies with distance (a miss item logs F5 0.03), but a spread below the
  0.05 floor is scaled down. The group gives little or no signal. After `--train-views`, off-centre small nodules
  will produce many such groups. EST, from val miss rates: 25-35 % of groups for > 0.3 offsets at first.
- The topology term (c_topo = 2, `spec.py:30`) acts on single-ring LIDC items as a large step at IoU_tc = 0.5. That
  helps a group with some hits and some misses. It does nothing for a group with no hits.
- LIT: Seg-Zero / VisionReasoner add L1 box/point rewards to IoU exactly to shape localization
  ([2503.06520](https://arxiv.org/abs/2503.06520), [2505.12081](https://arxiv.org/abs/2505.12081)). Perception-R1
  finds that reward design sets the ceiling for perception RL ([2504.07954](https://arxiv.org/abs/2504.07954)).

### 3.5 Perception capacity: frozen vision tower, 256 tokens, small objects

- CODE: LoRA r16 on the language model only; the vision tower and projector are excluded
  (`configs.py:90-95`, `VISION` regex). Each crop is upscaled with bicubic filtering to 512 px (`frame.py:74`), which
  gives 256 vision tokens (about 32 rendered px per token, so 3.1 CT px for W = 50 and 6.25 CT px for W = 100). A
  20-40 px nodule (5-7 px across) covers about 2 × 2 tokens in a W = 50 crop, or about one token in a W = 100 crop.
- MEASURED: the miss rate falls steadily with object fraction (40 % below 1 % of the crop, 3.4 % at 10-20 %).
  RESULTS.md (R1 diagnosis) reported MedSAM given the expert box at 0.77 on the shifted view. With the location
  given, a drawing model reaches 0.77; ours reaches 0.706 on found items. Locating is the bigger loss.
- Verdict: a plausible secondary limit for small objects off centre. Test it only after E1/E3 (it costs SFT time and
  slows every step).

### 3.6 No test-time selection, and the verifier exists only at training time

- Eval draws K = 8 and averages them (`final_eval_2026-11-08.md`); nothing selects. Legality (parse + ring checks) is
  **label-free**, so a legal-only filter is deployable. IoU_tc and F5 need the rater masks, so "the checker as an
  exact verifier" holds only during training.
- MEASURED (from item means): partial items (item-mean IoU_tc 0.1-0.5, 12.6 % of val) carry 37 % of the gap. These
  items mix hits and misses, which is where label-free selection (medoid / MBR by mutual IoU, or pixel-majority
  fusion) can help. Misses (mean < 0.1) have at most about one hit in 8, so selection rarely saves them.
- LIT: under label-free stopping, visual self-correction's +2.4 pp disappears; the model cannot tell which of its
  boxes is better (r ≈ 0.22) ([2606.13156](https://arxiv.org/abs/2606.13156)). Repeated sampling helps directly only
  with a verifier; without one, majority vote and reward models plateau
  ([Large Language Monkeys 2407.21787](https://arxiv.org/abs/2407.21787)). So report label-free selection, and give
  the oracle best-of-K only as a labelled upper bound.

### 3.7 Not the ceiling (checked and ruled out, or minor)

| Suspect | Finding | Verdict |
|---|---|---|
| 0-1000 quantization vs pixel size | Targets are pixel-corner polygons. A pixel is 20 (W = 50) or 10 (W = 100) grid units, and the grid points are exact (`simplify.py` `to_grid`). Quantization error is 0. | Not a cause. A pixel grid (0-50) would cut about 1 digit per coordinate (EST −25 % tokens). That is speed, not quality. |
| Vertex budget / token cap | 14.2 vertices mean, p95 32; cap 1,000 vertices (`text.py:59`) / 1,536 tokens (`configs.py:75`); truncation 0.12 % | Not a cause. |
| Rectilinear vs smooth targets | LIDC targets must rasterize to the rater mask pixel for pixel (`build.py:67`), so they are exact for an IoU computed on pixels | Not a cause. Smoothing would lower IoU against pixel masks. |
| SFT target quality / rater mixture | One SFT row per voting rater (`dataset.py:19-23`) gives a rater-mixture policy with high initial entropy, which RL then sharpens | OK. Gives the needed diversity at the start. |
| KL / beta | beta 1e-3 → about 0 effect; the ref pass costs about 3-4 s/step (code review #3) | Set beta = 0 for speed; no quality effect expected. |
| LR / schedule | 5e-5 constant LoRA (`prereg.py:64`, `g_probe.py:86`); grad norm 3-6e-3 and stable | Fine. |
| Group size G = 8 | Enough for centred items; too few for hard items (3.4) | Raise G only for hard items (E3). |
| Data size | 7,341 items × 3 epochs. The train vs centred-val gap is small (0.77 vs 0.744). | Not limiting for centred items. The shifted windows add 2 × 8,208 rows. |
| Legality | 99.4 % per sample at test; training defects about 0.1 % | The legal-only filter at inference fixes it (E0). |

### 3.8 Measurement issues the owner should know

- **IoU_tc definition (CODE, corrects the brief).** IoU_tc is the *plain* IoU against the best rater with the same
  topology signature, or against the vote if no rater matches (`gpucheck.cpp:1503-1521`, `score.py:8-9`).
  Disagreement pixels count as neither only in the **F5** band, not in IoU_tc.
- **The ceiling is not like-for-like.** The model gets the best of k raters; each doctor gets the best of the other
  k − 1. A fair "doctor level" claim scores the model against each set of k − 1 raters and averages over which rater
  is left out. On centred items the model sits 0.012 above the ceiling, which may be this bias alone (EST: the bias
  is 0.01-0.02).
- **The F5 gap is mostly misses.** F5 on found items is 0.896; misses pull the val mean to 0.806.
- **Comparing with the literature.** LIDC 2D work usually uses centred 64-128 px crops and per-rater Dice (for
  example ProSona Dice-max about 0.92, [ISBI 2026](https://doi.org/10.1109/isbi61048.2026.11515434); older papers
  report IoU 0.55-0.74 against inter-observer 0.54-0.72, [iW-Net table](https://www.nature.com/articles/s41598-019-48004-8/tables/3)).
  Our task also includes finding the object in shifted crops. Report one centred-protocol number next to it, or
  "top of field" cannot be checked.

## 4. Literature-backed menu of fixes

Expected effects are EST unless marked. Costs assume one RTX 5090 at $0.45-0.60/h and 23.6 s/step with
`--shared-prefix` (R5).

| # | Fix | Source | What it addresses here | Expected effect (EST) | Cost |
|---|---|---|---|---|---|
| F1 | RL on shifted + rescaled train windows (built: `--train-views`) | in-repo diagnosis; distribution match | 3.1 | +0.04 to +0.07 IoU_tc; misses 12.6 % → 6-8 % | 300 steps ≈ 2 h + eval ≈ $1.6 |
| F2 | Fix the IS mismatch: `vllm_importance_sampling_mode="token_truncate"` (TIS, cap 2-3) or `sequence_truncate` on a length-normalized ratio; then FP16 for trainer + vLLM | [2510.26788](https://arxiv.org/abs/2510.26788); [GSPO 2507.18071](https://arxiv.org/abs/2507.18071); [vLLM bitwise RL](https://vllm.ai/blog/2025-11-10-bitwise-consistent-train-inference) | 3.2 | Effective batch about 2x; +0.01-0.02, faster per step; fewer KL spikes | Config change; FP16 needs a check of TRL/vLLM LoRA FP16 paths. A/B 2 × 150 steps ≈ $1.2 |
| F3 | Label-free test-time selection: legal filter + medoid/MBR over K samples (mutual IoU), or pixel-majority fusion of K masks, then re-polygonize | MBR decoding; [Large Language Monkeys](https://arxiv.org/abs/2407.21787) (verifier vs no verifier); [2606.13156](https://arxiv.org/abs/2606.13156) (report label-free only) | 3.6, legality | legal → ~100 %; +0.01-0.04 IoU_tc; F5 up on partial items | Inference only: K = 8 already drawn. Needs per-sample texts saved: one eval rerun ≈ $0.25 |
| F4 | Dense localization reward: add a z-scored term such as exp(−d_centroid / r_union), or IoU of dilated masks, with its own floor; or a box-first output `{"box":[...],"nodule":[...]}` with a box-IoU/L1 reward | [Seg-Zero](https://arxiv.org/abs/2503.06520), [VisionReasoner](https://arxiv.org/abs/2505.12081), [Perception-R1](https://arxiv.org/abs/2504.07954), [PolyFormer](https://arxiv.org/abs/2302.07387) (box then polygon), [Text4Seg++](https://arxiv.org/abs/2509.06321) (box-wise descriptors) | 3.4 | Turns all-miss groups into signal; misses −30-50 % beyond F1 | Reward change = new reward version and new caches fingerprint (`spec.py:3-4`); 300 steps ≈ $1.6 |
| F5 | Online hard-example mining + curriculum: oversample items whose rollouts miss; skip groups with no spread and refill the batch (DAPO dynamic sampling); G = 16 for hard items | [DAPO 2503.14476](https://arxiv.org/abs/2503.14476); online difficulty filtering | 3.4, 3.3 | Concentrates compute where the gap is; +0.01-0.03 beyond F1 | Sampler code; same run as F4 |
| F6 | Entropy management for the hard windows: KL-Cov / Clip-Cov, or Pass@k advantage | [2505.22617](https://arxiv.org/abs/2505.22617); [Pass@k Training 2508.10751](https://arxiv.org/abs/2508.10751) | 3.3 | Keeps exploration on off-centre items; small on centred | Loss code; tune on a 150-step A/B |
| F7 | Rejection-sampling SFT / expert iteration: sample K on train windows, keep the checker-best sample per item (also when it beats a rater), then SFT | [RAFT 2304.06767](https://arxiv.org/abs/2304.06767); [ReST-EM 2312.06585](https://arxiv.org/abs/2312.06585) | 3.1 if RL on windows stalls; stable and cheap to debug | Similar to F1; stronger when most groups have a hit | ≈ 128k completions (≈ 1.1 h) + SFT ≈ 1 h ≈ $1.3 |
| F8 | Higher input resolution for small objects: 768-1024 px renders (576-1024 vision tokens), or adapt the vision tower (LoRA on vision blocks + merger) in SFT | Perception-R1 (perceptual perplexity); in-repo G-vision gate | 3.5 | Small-object misses −?; uncertain | SFT top-up $2-4; step time 1.5-2.5x |
| F9 | Backbone upgrade to the 4B class: Qwen3-VL-4B or Qwen3.5-4B (fits one 5090 with LoRA) | `MAIN_skill_lookups.md` (RefCOCO avg: Qwen3-VL-4B 88.2 vs 2B models 84.8; Qwen3.5-2B = Qwen3-VL-2B on box grounding) | 3.5 | Better grounding (+3.4 RefCOCO points on natural images); unknown on CT | Full SFT + RL ≈ $8-14 (step time about 2x) |
| F10 | Multi-turn self-refinement: the model sees its own rasterized outline and corrects it | [IVT 2606.13156](https://arxiv.org/abs/2606.13156) (gain vanishes label-free); [RSAgent 2512.24023](https://arxiv.org/abs/2512.24023), [IBISAgent](https://en.papernotes.org/CVPR2026/medical_imaging/ibisagent_reinforcing_pixel-level_visual_reasoning_in_mllms_for_universal_biomed/) (tool-based) | 3.6 | Doubtful without a verifier; 2x tokens + extra images | High (code + SFT traces + RL) |
| F11 | VLM box + MedSAM mask hybrid (Seg-Zero style), as a reference ceiling | [Seg-Zero](https://arxiv.org/abs/2503.06520), [Seg-R1 2506.22624](https://arxiv.org/abs/2506.22624), [SAM-R1 2505.22596](https://arxiv.org/abs/2505.22596) | 3.5 | MedSAM with a GT box: 0.77 (RESULTS.md); sets what drawing is worth once the box is right | Low (inference only); changes the product if adopted |
| F12 | Algorithm swaps: Dr.GRPO (in use), RLOO, REINFORCE++ global normalization, GSPO | [Dr.GRPO 2503.20783](https://arxiv.org/abs/2503.20783), [RLOO 2402.14740](https://arxiv.org/abs/2402.14740), [REINFORCE++ 2501.03262](https://arxiv.org/abs/2501.03262), GSPO | 3.3 | Small. The estimator is not the bottleneck; GSPO's value here is mismatch tolerance (F2) | Config / loss code |
| F13 | Coordinate schemes: pixel-unit integers, location tokens (Florence-2 / Pix2Seq bins), number-token loss in SFT | [Pix2Seq 2109.10852](https://arxiv.org/abs/2109.10852), [Florence-2 2311.06242](https://arxiv.org/abs/2311.06242), [NTL 2411.02083](https://arxiv.org/abs/2411.02083), [VistaLLM 2312.12423](https://arxiv.org/abs/2312.12423) | 3.7 | ≈ 0 on IoU (no quantization loss now); fewer tokens | SFT redo; low priority |

## 5. Ranked experiment plan

Budget note: `rl_plateau_diagnosis` says Phase 2 has $0.60 left of the $17.50 hard stop. Every row below needs the
owner's budget.

| Rank | Experiment | What changes | Expected (EST) | $ (5090) | Go / no-go |
|---|---|---|---|---|---|
| **E0** | Test-time selection on the existing RL model | Rerun `final_eval` on val (and test) with K = 8, **saving every sample's text**. Offline on the Mac: legal filter; medoid by mutual IoU; pixel-majority fusion; oracle best-of-8 (labelled as an upper bound); LOO-fair ceiling (3.8) | legal ≥ 99.9 %; IoU_tc +0.01-0.04 | ≈ $0.25 (≈ 25 min with boot) | **Go** (adopt as the default inference) if medoid or fusion raises val IoU_tc by ≥ 0.015 with a paired-bootstrap CI above 0, and test agrees in sign. Also gives oracle−medoid (how much a verifier would be worth) and the fair ceiling. |
| **E1** | `rl_shifted.yaml` as built (300 steps from rl-full-3ep, `--shared-prefix`) **plus logging only**: share of sequences with IS ratio 0, and share of zero-signal groups per window kind | Train distribution | IoU_tc 0.645 → 0.69-0.71; misses ≤ 7 % | ≈ $1.6 (2.6 h) | **Go** if val IoU_tc ≥ 0.66 and the miss rate at offset ≥ 0.3 falls below 20 % (from 31.8 %), with centred items not below 0.73. **No-go** (switch to E3/E7) if the share of zero-signal groups on shifted windows stays > 25 % after 100 steps. |
| **E2** | IS / precision A/B, 150 steps each from the E1 model: (a) TRL default; (b) `token_truncate` cap 2; (c) FP16 trainer + vLLM, if supported | Gradient quality | Per-step slope of training IoU on windows up ≥ 1.3x; +0.01-0.02 val | ≈ $1.2-1.8 | **Go** for (b) or (c) if the IS mean is ≥ 0.97 with < 0.5 % of sequences masked and the val IoU_tc at 150 steps is not lower than (a) (CI). (c) also needs no fp16 overflow (NaN rate 0). |
| **E3** | Localization term + hard-example mining (F4 + F5): a centroid-distance or dilated-IoU term with floor 0.05; oversample items whose last rollouts missed (×4); G = 16 on them | Reward v5 + sampler | Misses −30-50 % vs E1; +0.02-0.04 | ≈ $1.6 + Mac lab round for reward v5 (caches rebuild) | **Go** if val misses fall by ≥ 25 % relative to E1 and found-item IoU does not drop by more than 0.01. Run reward hacking checks first (for example, does a large outline game the centroid term? g_area should guard against it). |
| **E4** | Rejection-sampling SFT (F7) | Only if E1 shows many groups with no hits | Similar to E1 | ≈ $1.3 | Run only if E1 fails its no-go rule. **Go** if the probe IoU at offset ≥ 0.3 rises ≥ 0.05 over E1. |
| **E5** | Entropy control (F6), KL-Cov k = 2e-3 | Loss | +0.01 on hard windows | ≈ $0.8 (150 steps) | **Go** if entropy stays ≥ 2x baseline and val IoU is not lower. |
| **E6** | Resolution 768/1024 px (+ SFT top-up) or vision-tower LoRA | Perception | Small-object misses −? | $2-4 | Only if misses below 1 % of the crop area are still > 25 % after E1/E3. |
| **E7** | Box + MedSAM hybrid as a reference | Inference only | Sets the drawing ceiling | ≈ $0.3 | A reference point, not product. Run it if the owner wants a "best possible with SAM" bar. |
| **E8** | 4B rung (Qwen3-VL-4B or Qwen3.5-4B; not Qwen3.5-2B, which ties the current model on grounding) | Backbone | Unknown on CT | $8-14 | Only after E0-E3, if the found-item IoU stays < 0.72 or the misses do not move. |

Order: E0 → E1 → E2 → E3; E4-E8 are conditional. E0 and the E1 logging are almost free and decide what comes after.
Cumulative EST if E0-E3 all pass: val/test IoU_tc about 0.70-0.73, F5 about 0.86-0.89 (found-item F5 is already
0.896), legal ~100 %. This is near the (pessimistic) doctor ceiling of 0.72-0.73 on this shifted protocol. That
claim must use the fair ceiling from E0.

## 6. Do not spend on

- More epochs on centred frames: the training reward is saturated and centred val already beats the LOO ceiling.
- Clip-higher / PPO-clip tuning / GSPO clipping for their own sake: the clip never fires with one update per batch.
- Coordinate format or vertex-budget changes for quality: quantization error is 0 and truncation is 0.12 %.
- Multi-turn self-refinement before E0 shows how large the oracle-minus-label-free gap is.
- Raising beta: KL is not the limit. beta = 0 saves 3-4 s/step (code review #3).

## 7. Open points I could not settle

- Per-sample texts are not saved (`items.jsonl` has item means only), so I could not measure best-of-8, medoid or
  the fair ceiling. E0 fixes this.
- The test view is not on the Mac (`data/eval_views` has val only), so my bucket tables use val. The plateau note's
  test tables agree in direction.
- The IS ratio distribution is inferred from the mean, min and max. The masked share and ESS are EST until E1's
  extra logging.
- How much of the mismatch is the bf16 merged LoRA and how much is kernels: an FP16 or merged-weights trainer check
  (E2c) would split them.
