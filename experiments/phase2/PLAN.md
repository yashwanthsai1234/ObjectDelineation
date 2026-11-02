# Phase 2 plan: SFT of the smallest polygon VLM on Vast AI (draft C)

Written 2026-10-20 by plan agent C (wave C of `RESEARCH.md`). Part of [ROADMAP.md](../../ROADMAP.md); follows
[Phase 1 RESULTS §8](../phase1/RESULTS.md#8-phase-2-handoff). Status: **approved (judge round 4, `JUDGE_R4.md`); revision 8** (R4's optional polish, §9.2d) on top of revision 7 (judge round 3, `JUDGE_R3.md`: APPROVE WITH CHANGES; G-prior (a) report-only, 8 SHOULD-FIX; §9.2c) on top of revision 6 (B2: D12 rule R-B, signal panel folded into R1 / R2) and revision 5 (B1: ladder confirmed, vision tower frozen under LoRA, G-vision gate, Qwen3.5 R0 arm) and revision 4 (judge round 2, `JUDGE_R2.md`: 6 MUST-FIX resolved, §9.2b) and revision 3 (owner cost directive, 2026-10-20: one recipe bundle, Vast only for C1 / C2 / G-probe, §6 and §9.1), after the wave-D judge's round 1
(`JUDGE_R1.md`: APPROVE WITH CHANGES, 6 MUST-FIX, 12 SHOULD-FIX; how each was addressed: §9.2). **Wave B is complete** (B6
moved into the bundle and the Vast runs) and no PENDING marker remains. B3's
Addendum A (hole-calibrated topology probe, C2 primary metric and pass rule; `lab/B3/REPORT.md` §Addendum A) is in and
is written into §1, D20, §5 and §9 (revision 2).

How this plan was made:
- **6 research notes** (wave A): `docs/research/notes/Phase 2 SFT/A1_models_tokens.md` … `A6_performance.md`
  (cited as A1 §x … A6 §x).
- **5 lab tests** (wave B): `lab/B3/REPORT.md` (probes, offline augmentation), `lab/B4/REPORT.md` (reward serving),
  `lab/B5/REPORT.md` (environments), `lab/B1/REPORT.md` (model learnability, MPS), `lab/B2/REPORT.md` (SFT depth → GRPO
  signal, interpolation) are done. B6 moved to Vast A/Bs (RESEARCH.md scope change 2026-10-20); revision 3 replaced all per-component
  A/Bs with one evidence-chosen recipe bundle and a single pre-registered fallback (§6).
- Phase 1 evidence: `experiments/phase1/PLAN.md` (P1-D1…D14), `RESULTS.md` (§2 gates, §3 hashes, §7 decisions, §8
  handoff), `impl/ISSUES.md` (I34–I45 open).

Conventions: Phase 2 decisions are **D1…D27**; Phase 1's are cited as **P1-D10** etc. Every number names its source.
Mac numbers were taken under load (load average 8–48), so Mac timings are upper bounds; every timing gate is
re-measured on Vast.

## 1. Goal, scope, claims, non-goals

**Goal.** Teach the smallest VLM that can to emit legal, precise D10 polygons for LIDC nodules (and RIGA disc/cup when
permitted, D16), pick an SFT checkpoint that is a good GRPO start, and pass the revised Phase 2 entry gates (§5) on
Vast GPUs, within the budget in §4.9.

**Claims tested in Phase 2** (each pre-registered, §5.2):
- **C1 (smallest viable).** The smallest rung of the ladder (D1) whose selected checkpoint passes every entry gate
  (§5.1) on the judged datasets. Rungs are screened smallest first and the ladder **stops at the first rung that
  passes**; larger rungs are not run. Qwen3-VL-2B runs only as the next rung or as the conditional capacity control of
  §6.1 (R-ctrl), never as an unconditional control.
- **C2 (A0 vs A0′: holes).** A0 is "A0-multi" (all outer rings, parts kept, `"holes":[]`; D2), so A0 and A0′ differ
  **only in holes**; parts are already expressible in A0. C2 claims: on a synthetic topology probe whose hole sizes
  follow real LIDC holes (D20; B3 Addendum A: median 12 px, 44 % ≤ 10 px), A0′ recovers holes that A0 cannot.
  A0 scores 0 hole recall in every bin by construction (B3 Addendum A4; first probe: non-simple hit 0.50, B3 §3 P2),
  so the paired difference A0′ − A0 equals A0′'s hole recall and C2 is stated on A0′'s hole recall (§5.2).
  C2 is the **single confirmatory test** of Phase 2 (§5.0). A failure is reported and holes are left to GRPO
  (A1 H-A1-6 fail branch). ROADMAP Phase 2 "done when" is rescoped to holes accordingly (2026-10-20).

**In scope (critical path, §6).** Vast operations (image, scripts, guardrails); the library changes in §7; **one
recipe bundle** (§6.0) used in every run; and **only these Vast runs** (§6.1): box checks with folded measurements
(R0), screening smallest first with a stop at the first passing rung (R1), G-probe on the selected checkpoint (R2),
C2 on the chosen rung with 3 A0′ seeds and 1 A0 seed (R3–R5, revision 4), plus at most one bundle-level diagnostic run (§6.3) and at most one C2-outcome
conditional run (×8 seed on a partial, capacity control on a fail; §5.2). Inference-only baselines (image-blind prior, MedSAM, SAM 2.1-tiny) run on the Mac CPU. **Removed or
triggered** (§6.4): every per-component A/B (E1–E10), coordinate tokens, mixing, nnU-Net training, K = 64 / I35,
HD95/ASSD/centroid in `objdel.eval`, the Phase 3 environment smoke (moved to Phase 3 start), RIGA (S6).

**Non-goals.**
- GRPO training runs (Phase 3). Phase 2 runs only the 20-step G-probe on the selected checkpoint (§5.1), which also
  measures the GRPO step time.
- RIGA on any rented host before written permission (D16); publishing any RIGA-derived crop, mask or polygon (A4 §5.2).
- Coordinate tokens, constrained decoding in training, the C++ batch scorer: each only behind its pre-registered
  trigger (§6.4, D13, D14).
- Multi-GPU, H100/H200/B200 (A3 §3.4), 3D, natural images, the CUDA port.

**Done when** every §5.1 gate is decided for every ladder rung that was run, C1 and C2 have verdicts, the spend is
within the §4.9 hard stop, and `make test`, `make lint`, `make cov` (≥ 90 %) pass with the §7 code merged.

## 2. Decisions

| # | Topic | Decision | Rejected (why) | Evidence |
|---|---|---|---|---|
| D1 | **Model ladder** (B1 filled, revision 5) | **LFM2.5-VL-450M → Qwen3.5-0.8B → Qwen3-VL-2B**, smallest first, **stop at the first rung that passes** (owner directive), all three with a vLLM path for Phase 3 (B5 Q2); Qwen3-VL-2B is the third rung (owner approval, §8.2) and otherwise only the conditional capacity control R-ctrl (§6.1). **B1 verdicts** (`lab/B1/REPORT.md` §3): H-A1-2 **INCONCLUSIVE** under its pre-registered rule (every shifted median IoU_tc is 0 and no arm learns the format within 150 Mac steps, so the ranking statistic ties) and H-A1-1 **INCONCLUSIVE** (Qwen3.5 shifted RIGA disc box IoU 0.41: not ≥ 0.6, not < 0.4 on both datasets), so the default stands (B1 rule 4) and rule 2 cannot remove Qwen3.5. **Smallest first, not side by side.** B1 §0.5 recommends screening LFM and Qwen3.5 side by side in the first Vast session because Mac depth could not separate them: Qwen3.5 learns legal syntax ~4× faster in steps (legal 82–83 % at 36 steps vs LFM 42–59 %) and holds legality on shifted crops (83 % vs 71 % at 150 steps), while LFM has the better centred median IoU_tc at 150 steps (+0.065 [0.026, 0.104]), fewer tokens per vertex (6.5 vs 9.5) and kept its vision under the all-linear adapter (B1 §2.2–2.5). The plan takes the cheap part of that advice: **R0 runs a 30-step Qwen3.5 arm** (s/step, peak memory → D6, GDN kernel route → D15, base-weight parity → D13, base G-vision references; ≈ 0.2 GPU-h, ≈ $0.13; §4.7 step 8, §4.9), so a Qwen3.5 R1 after an LFM failure starts in the same V1 session with a measured cost. The arm has **no pick rule**: it decides nothing about order, because C1 is "smallest passing" and rung order is by parameter count whatever B1's ranking. **SmolVLM-256M: control only, never screened on Vast, never C1** (B1: last on every metric at 150 steps, vision collapse under the all-linear adapter, no zero-shot prior; B1 §2.1, §2.3, §2.5, §3; and TRL #6294 open, PRs #6718 / #6719 unmerged, `grpo_config` forces `use_vllm=False`, excluded from TRL's GRPO CI over bf16, A1 §2); B1 on the Mac is its record. Qwen3-VL-4B or InternVL3.5-1B-HF are extra rungs only after the 2B rung fails (§8.2, owner approval). **B1 caveats** (B1 §1.5, §5): LoRA r16 on MPS at 36 / 150 steps with one LR, not the Vast recipe; Qwen3.5 single-seed with a bf16 base (AMEND 3); the 150-step tier and LM-only arms exploratory (AMEND 2, 4, 5); LFM's seed-to-seed swing (16 legal points) exceeds its item CI; B1's shifted view used 0 px clearance and no `other_nodules` check (`lab/B1/common.py`), looser than M1's | **LFM and Qwen3.5 both screened in full in V1** (the literal reading of B1 §0.5): runs Qwen3.5 unconditionally against "stop at the first passing rung"; +1.7–3.9 GPU-h / +$1.1–2.6 on scenario A, ≈ 20 GPU-h / ≈ $13 at MFU 0.05, past the $10 alarm (§4.9). **Short side-by-side arms with a pick or skip rule**: B1's 36-step ranking did not predict its 150-step ranking (LFM legal 55 → 83 %, and LFM overtook Qwen3.5 on centred IoU; B1 §2.2 vs §2.3), so a short screen measures learning speed, not the gate endpoint, and skipping LFM on it would turn C1 into "smallest passing among the rungs run". Qwen3.5 first (larger; C1 is smallest passing). SmolVLM-256M as first rung (no prior, excluded from TRL GRPO CI over bf16, TRL #6294 vLLM bug, 1 visual token per 64 px; CropVLM needed ~3 A100-h for 4-number boxes: A1 §2, §5.1). SmolVLM-500M (P1-D12). Florence-2 (encoder-decoder, no TRL GRPO; baseline only). Ovis2.5, Rex-Omni, Isaac, Moondream 3, Gemma 3n/4, PaliGemma 2, MedGemma (licence, no TRL path, size or no prior: A1 §6) | B1 §0, §2–§4; A1 §2 table, §5.1; Qwen3.5-0.8B: Apache-2.0, RefCOCO 79.3, in TRL GRPO CI; answers natively in `bbox_2d` inside a JSON fence although its card does not document it (B1 §2.1); LFM: RefCOCO-M 81.28; tokens/vertex measured 6.5 vs 9.5 (B1 §2.5; A1 §1.1 estimated 5.7–6.1 vs 8.6–8.9) |
| D2 | Output format | **Keep P1-D10** for every model (compact JSON, ints 0–1000, CW outer / CCW holes, top-most-then-left-most start, `"holes":[]` always). **A0** = same grammar with `holes` emptied, parts kept ("A0-multi"); A0-single reported only as a ceiling. Qwen3.5 in non-thinking mode, trained with TRL's `qwen3_5_nothink_training` template (prefix-consistent with the generation prompt, checked in B1 §1.1); its 4 `<think>\n\n</think>\n\n` tokens inside the SFT loss span are benign (same sequence as the generation prompt; B5 Q2 ‡). Zero-shot, no model emits D10 (0 % parse for all three in B1 §2.1) and the native box priors are weak (B1 §2.1, H-A1-1 inconclusive), so the plan relies on no zero-shot prior; LFM replaces its [0, 1]-float convention with ours within 36 steps (parse 94–99 %, B1 §4.5) | Delta (−25 % tokens, but a ±1 delta-px slip costs 13 IoU pts on LIDC; A1 §1.2). Pixel units (−23 % LIDC tokens; 1 px per last-digit slip, −2.8 IoU; cannot express bevel vertices in 0.34–0.47 % of targets). Coarse grids (A1 §1.5, RIGA: G = 100 costs −2.7 IoU pts and −2.2 pts legality; G = 250 costs −1 IoU pt and −1.5 pts legality). Coordinate tokens by default (SFT-scale ties: GETok 92.9 vs 93.0, Shikra; init-sensitive, HiVG 0.44 vs 0.72): only via E8. Random start (−40 IoU, LPM T3) | A1 §0, §1, §3, §6; every target coordinate in [168, 832] (A1 §1.1) |
| D3 | SFT framework | **TRL 1.16 `SFTTrainer`** via `objdel.train.configs.sft_config`; assistant-only NLL with TRL's training-template swap (stop token trained for all models; a custom collator must pass `get_training_chat_template`) | Unsloth, LLaMA-Factory, ms-swift, Axolotl, NeMo, oumi (model-coverage gaps, second template/mask path); torchtune (unmaintained). Packing / padding-free (refused by TRL 1.16 for vision datasets; ≤ 0.2 % padding left after grouping) | A2 §1.1, §1.3, §2 |
| D4 | Throughput knobs | `train_sampling_strategy="group_by_length"` with a `length` column (completion tokens) from `sft_rows`; `dataloader_num_workers=6`, persistent, prefetch 2; gradient checkpointing **off** for ≤ 0.8B models when peak ≤ 20 GB, on (selective) for Qwen3-VL-2B; SDPA attention; bf16 (`configs._bf16`); `PYTORCH_ALLOC_CONF=expandable_segments:True`. Measured in the R0 timing prelude (§6.2), not a separate A/B | Liger (deprecated path in TRL SFT, no `lfm2_vl`); FA3/FA4 (wrong hardware, sequences ≤ 2.7k); torch.compile (every batch a new shape); FSDP/DeepSpeed (single GPU) | A2 §1.2 (padding 7–32 % → ≤ 0.2 %), §1.4 (collate 110–350 ms per batch of 8, one core), §3; A6 §1a (10–25 ms CPU per row), §2 |
| D5 | Schedule and depth | **WSD**: 3 % warmup, constant LR, 15 % linear decay to 0; one stable run of **2 epochs** (each epoch a fresh augmented view, D10) with **decay branches at 0.5 / 1 / 2 epochs**, all kept; a third stable epoch only if the 2-epoch branch is the D12 pick and still improving on shifted-val IoU_tc by ≥ 0.01 over the 1-epoch branch (owner approval, §4.9). Compute per screening run ≈ **2.5 LIDC epoch-equivalents** (2 stable + 0.075 + 0.15 + 0.3 decay; §4.9 uses this). D12 picks the depth among these branches (B2 found no over-training signal at Mac depth, B2 §3). Peak LR 5e-5 full FT, 5e-4 LoRA (10×); no LR sweep. Effective batch 32 sequences. `seed` and `data_seed` set. **Why fresh views per epoch:** the crop leak (A4 S5–S6, D10), not data uniqueness. Quagmires (2510.01624, abstract re-read 2026-10-20): "SFT training on unique examples for one epoch underperforms training on half examples for two epochs, either after SFT or SFT-then-RL", so 1–2 epochs is the evidence-backed range and the 2-epoch cap loses little | Cosine single-depth run (one depth only). 3 stable epochs + 3-epoch branch by default (revision 2; +1.5 epoch-equivalents per run for a depth the evidence does not favour). Easy-first / short-only curriculum (Quagmires: short-only "often leads to worse outcome after RL") | A2 §4.3 (WSD); SimpleSeg (2601.19228, Tables 4–5, re-read 2026-10-20): **SFT** peak 5e-5 cosine to 2e-6, 3 % warmup; RL constant 2e-6 (GSPO, Muon, 7B), §4.7; Quagmires abstract (held-out loss and pass@large-k predict RL outcome better than SFT score); P1 `configs.py` LRs |
| D6 | Full FT vs LoRA | LFM2.5-VL-450M and SmolVLM: **full FT**. Qwen3.5-0.8B: full FT if the R0 Qwen3.5 arm (§4.7 step 8) shows peak ≤ 20 GB on the 4090 at micro-batch 8, else LoRA r16. Every LoRA adapter is **r16 on the language model's linear layers only, vision tower and projector frozen** (PEFT `exclude_modules`; `configs.lora` changes accordingly, M2; D7), merged after SFT. Qwen3-VL-2B: **LoRA**, on the RTX PRO 6000 (A2 §4.1 caveat: LoRA "is less tolerant of large batch sizes"; effective batch 32 is within the usual range, recorded); full FT (~38 GB) only as a fallback if LoRA misses format at every branch | Qwen3-VL-2B full FT on 24–32 GB (does not fit; 85 min/epoch on a 4090 even if it did). DoRA / rsLoRA (no gain expected at r16) | A2 §4.1 (LoRA Without Regret: all layers, 10× LR), §1.2; A3 §3.4; `estimate_cost.py` (full FT 6 / 10 / 38 GB) |
| D7 | Vision tower | **LoRA (Qwen3-VL-2B; Qwen3.5 if D6 picks LoRA): vision tower and projector frozen** (PEFT `exclude_modules`, B1's `common.lora("lm")` regex). **Full FT (LFM, Qwen3.5): vision tower and projector at 0.1× the LM learning rate** (5e-6), via optimizer param groups. Both checked by the **G-vision** gate on every SFT branch (§5.1) and by an in-training projector check (M8): if at a rolling checkpoint **both** G-vision (a) ratio < 0.5 **and** (b) TV < 0.002 (B1's collapse hit both; (b) costs 8 items × 24 tokens on the GPU), the run is aborted and restarted once with tower and projector frozen (a pre-registered repair, not the §6.3 diagnostic); one of the two alone is warn-only (logged) and the decision is left to `select_ckpt`'s G-vision (revision 7). **Why:** the all-linear adapter at lr 5e-4 (`configs.lora`, vision included) erased image information for Qwen3.5 and SmolVLM within 36 steps (projector cosine distance between images 0.297 → 0.019 and 0.171 → 0.003; next-token TV real vs black or another item ≤ 0.0007 vs 0.02–0.08 at base); LFM kept its vision (0.471 → 0.321–0.373); LM-only LoRA kept Qwen3.5's features with unchanged legality (82 %) and centred IoU (0.20 vs 0.19) (B1 §0.2, §2.4–2.5). Frozen, not 0.1×, for LoRA because frozen is the setting B1 tested; 0.1× on vision LoRA weights is untested. Full FT keeps 0.1× because A2 §4.2's evidence supports unfreezing for localisation (VectorLLM, Visual Position Prompt, VLM4VLA) and a 5e-6 full-FT update is a different regime from B1's 5e-4 LoRA one, so it is checked rather than assumed (B1 §4.2; H-A2-6 now P1). Phase 3 note: vLLM encoder cache must be reset after each weight sync when the tower trains (TRL #7198) | All-linear LoRA incl. vision at the LM LR (revision 4's Qwen3-VL-2B setting; collapses vision, B1 §2.5). 0.1× vision LoRA LR (untested; B1 offers frozen as the tested alternative). Freezing for full FT by default (discards A2 §4.2's localisation evidence; G-vision and the repair cover the risk). A frozen / unfrozen A/B (owner directive; E4 trigger in §6.4). Qwen recipe LR 2e-7–1e-6 frozen vision (tuned for mixed image/video data, not a new format) | B1 §0.2, §2.4–2.5, §4.2; A2 §4.2; A5 §3, H-A5-9 |
| D8 | Loss | **Assistant-only NLL. Number-token loss NOT in the bundle.** The evidence is too weak to adopt without a test: NTL's gains are on math and text regression (math interpolation 64 → 75 %, rJokes RMSE 2.01 → 1.81; A1 §3, A2 §4.5), with no VLM-coordinate follow-up found (A1 §3; A2 §10); the closest polygon evidence, Gaussian label smoothing on coordinate tokens, left IoU flat (83.77 → 83.67) and hurt on WHU (LPM, A1 §3); NTL is designed for single-digit tokenizers while the first rung, LFM, uses one token per 0–999 value (A1 §1.5). Its expected benefit (boundary error) is exactly what G-precision and legality measure, so it is the triggered option if those fail (§6.4; M6 built only then) | Adopting NTL by default (unverified on VLM coordinates; new loss code on the critical path). Label smoothing (deprecated in TRL SFT, extra forward pass, does not repair over-trained checkpoints) **although A4 §4.3 recommends 0.1**. DFT (mode-seeking, entropy collapse; H-A2-8, untested: B2 ran no training, B2 §0). GEM, CHORD (parked) | A1 §3; A2 §1.5 (prototype matches TRL's fused NLL), §4.5, §10 |
| D9 | Targets | **All voting raters as rows** (P1 I27; no consensus/STAPLE). **Per-variant dihedral targets built offline** from transformed masks with **crop-before-trace** (B3 P4: hash `1af819f4…` unchanged, RIGA −23 % per target). Every view has its own `target_set_hash`, the K5/K6/R1 re-check and the per-variant token-budget filter (I5) | Naive vertex transform + canonicalize (breaks the IoU ≥ 0.99 certificate on 7.3 % of RIGA cups; A4 S1). On-the-fly target rebuild (26 ms per RIGA target in the collator; A6 §1b). One random rater per epoch (same expected gradient as all raters, more bookkeeping; A2 §4.6) | B3 P3 (215,208 variant targets, 0 errors, two runs byte-identical, 16.5–18.1 CPU-min); A4 S1; A6 §4 |
| D10 | Augmentation | **Crop shift + scale jitter is mandatory** in every run (the centred frame leaks location and size: the image-blind prior passes RIGA's gates on centred frames and collapses under shift, A4 S5–S6; and every B1 arm trained on centred crops drew at the frame centre: on shifted crops 94–98 % of parsed samples lie closer to the frame centre than to the nodule, median IoU_tc 0, every arm significantly below the image-blind prior, B1 §0.1, §2.5). Shift: uniform, union inside with ≥ 2 px clearance, and no other outlined nodule (`meta.other_nodules`) inside the window (B3-2 rule kept). Scale: margin drawn from [8 px, 60 %] so S moves between divisors (LIDC 50/100/125, RIGA 200/250/500), S ≤ 500, scale 1. **D4 for LIDC; horizontal flip only for RIGA** (fundus orientation is physiological). One op per rater row: op = hash(seed, epoch, item, row). **CT window jitter** L −300 ± 60, W 1400 ± 80 HU on raw HU, noise σ ≤ 10 HU, applied to **int16 HU context crops written from `data/store`** (M1; `data/prepared` images are already lung-windowed uint8, so jitter cannot be applied to them); random CT windowing follows 2510.08116 (liver-tumour CT, a single study), so this component is the weakest-evidenced in the bundle; RIGA brightness/contrast 0.75–1.25, gamma 0.7–1.5, p 0.15–0.3, no hue. Images are cut and transformed at frame resolution **before** the 512-px Pillow upscale | Elastic, non-90° rotation, mixup, cutout, copy-paste, random start vertex (masks resample and need recertification; copy-paste changes the ground truth; A4 §4.2). RIGA vertical flips / rotations (non-physiological; AADG). Transforming the upscaled LIDC image (rotations differ by up to 16 grey levels; A4 S1) | A4 §0.1, §4.2, S1, S5–S6; B1 §2.5 centre bias; ROADMAP's "flips, 90° rotations" for both datasets is superseded for RIGA |
| D11 | Data mixing | All Phase 2 Vast runs are **LIDC-only** (D16; LIDC carries the topology claim). A RIGA joint run, if RIGA is permitted, uses natural proportions by targets (A2 §4.8: RIGA 12 % of rows, ~30 % of answer tokens); the T = 2 mixing A/B (E5) is removed (trigger §6.4). **Non-simple LIDC train targets ×4 in every run** (bundle, §6.0); A0 runs use the same repeated rows with holes emptied | Separate per-dataset models (A2 §4.8). A ×1/×4/×8 oversampling A/B (E1, removed) | A2 §4.8, §7.4; A1 §5.2; P1 RESULTS §7 #9 |
| D12 | Checkpoint selection (B2 filled, revision 6: rule R-B) | **Eligible:** WSD branches (0.5 / 1 / 2 epochs) passing G-format (parse ≥ 98 %, legal ≥ 90 %) and G-signal (duplicates ≤ 10 %, σ-share ≥ 50 %) at **T = 1**, K = 8 on the shifted val, and G-vision (§5.1). **Select the deepest eligible branch.** **Signal veto:** take the next shallower eligible branch instead only if both hold: (i) its geometric reward spread **sd_geom = √(p·Var_legal(R))** on the shifted val (R = the G-probe raw composite, illegal = 0; per item over its K = 8 group, then averaged) is higher, with the paired nodule-cluster CI excluding 0 **and** by ≥ 15 % relative; (ii) its shifted-val median IoU_tc is not significantly lower (paired CI). Applied pairwise downward from the deepest branch; ties stay deep. **Veto statistic: sd_geom, a post-hoc statistic, adopted by this decision and labelled so** (B2 §1 "post-hoc", §4 note, §7). Why not r_std (B2's pre-registered statistic): Var(R) = p(1−p)·μ_L² + p·Var_L(R) exactly, and the legality part rewards a branch for being *less* legal; it carried 55–92 % of the group variance in B2 and the whole 36 → 150-step fall of r_std (B2 §2.1). Inside the eligible set p(1−p) still moves from 0.09 to 0.02 between legal 90 % and 98 %, enough by itself to clear a 15 % margin. **Checked by R2's folded measurement** (§6.2): per-item Spearman of step-0 sd_geom, and of step-0 r_std, with the G-probe gain. That per-item correlation is confounded by legality headroom (items near p ≈ 0.5 have the most room to gain, which mechanically favours r_std), so it can only **retain** sd_geom or **flag** disagreement; switching to r_std needs the branch-level comparison (D12(a) / (b), when a second branch is probed) and is a Phase 3 planning decision, not an automatic rule (revision 7). Both statistics are reported per branch with a "statistics disagree" flag when the veto would fire under one only; the choice is written into `prereg.py` (M8) before R1. **Never select by:** greedy legal or IoU; σ-share, active-group share or distinct share (saturated at 98–100 %, B2 §2.1); total r_std or headroom (fall mechanically as legality rises). **Reported per branch, not used:** r_std, headroom, mixed-legality share, legal-but-wrong share, mean \|A_v4\|, TRL entropy, held-out NLL, pass@8 (B2 §4). The selected branch then runs **G-probe** (§5.1); a G-probe result never climbs the ladder on its own. **D12(c) repair** (selected branch not PASS on G-probe): the **next-earlier eligible branch**, probed once. Base ↔ LoRA-delta interpolation is **not** a repair (B2 H-B2-2 REJECTED: α = 0.75 loses 7–9 legal pts on both views, α = 0.5 reaches 38 % legal and 10.5 % duplicates; B2 §2.2). WiSE-FT between two SFT branches (0.5-epoch ↔ selected, A2's H-A2-7 proper, untested: B2 could only run the base ↔ late form) is scored beside the repair branch when the repair fires, at inference cost only (merge at α = 0.5, shifted val K = 8, ≈ 0.05 GPU-h), report-only, for Phase 3. **D12(a) / (b)** (does the signal branch out-gain the greedy one under GRPO?) are not testable on the Mac (B2 §3, §6): when a second branch is probed, both paired gains are reported side by side. H-A2-3 is INCONCLUSIVE at Mac depth (B2 H-B2-1: LFM's r_std and headroom fall beyond seed noise, all in the legality part; the geometric part is flat or rises), so B2 shows no over-training at ≤ 150 LoRA steps; the veto exists for the 2-epoch full-FT regime where 2606.18487's rank inversion may appear | **R-A**, revisions 1–5's signal-first lexicographic rule (σ-share, duplicates, pass@8, NLL, IoU): its first keys are saturated (σ-share and distinct 98–100 %, duplicates ≤ 1.4 % at every depth), so it would decide on noise, and B2's pre-registered rule keeps R-A only if H-B2-1 is ACCEPTED (B2 §1, §3). **r_std as the veto statistic** (pre-registered, but it penalises legality gains; above). Total r_std or headroom as selection keys. Base ↔ LoRA interpolation as the repair (H-B2-2 REJECT). Selecting by SFT score (rank inversion: pre-RL pass@1 up, peak GRPO pass@10 0.806 → 0.481; 2606.18487) | B2 §0–§7 (`lab/B2/REPORT.md`); A2 §4.10, H-A2-3 / H-A2-7; 2606.18487 (advantage variance p(1−p)(g−1)/g, rank inversion); 2504.10478 (early ↔ late WiSE-FT recovers pass@k); 2510.01624; RESULTS §8 GRPO-signal gate |
| D13 | Eval sampling | **vLLM 0.32.0 offline `LLM.generate` in `/opt/vllm`** (no objdel there), fed by a JSONL bridge from the training env (chat text with `tokenize=False`, image path, `mm_processor_kwargs` = `processor_kwargs` minus `padding_side`); returns token ids and texts; truncation (TRL's rule), parsing and scoring stay in `/opt/venv`. `SamplingParams(n=8, temperature=1, top_p=1, top_k=0, min_p=0, repetition_penalty=1, presence_penalty=0, frequency_penalty=0, max_tokens=MAX_COMPLETION[model], stop_token_ids=rewards.stop_ids, seed)` all explicit, and the same values passed explicitly to HF generate and to G-probe's GRPO: they override the cards' defaults, Qwen3.5 T 0.7 / top_p 0.8 / top_k 20 / presence 1.5 and LFM T 0.1 / min_p 0.15 / repetition 1.05 (B1 §1.4, §4.5). `MAX_COMPLETION` 1,024 for LFM and Qwen3.5 on LIDC (longest train target 559 / 876–881 tokens; truncation ≤ 0.2 % at 1,024; B1 §1.5, §2.6, §4.5). Its numbers count only after H-A5-1 (prompt ids 100 % equal) and H-A5-2 (logprob and distribution parity) pass on CUDA. **`objdel.train.sample` (HF generate) is the reference and fallback.** **GPU sequencing:** sampling is never done by the trainer process. `select_ckpt` is a separate process started only after `train_sft` has exited (branch checkpoints are model-only on disk); it does not initialise CUDA itself and launches `sample_vllm` (asserting `torch.cuda.is_initialized()` is False in the launcher), which checks `torch.cuda.mem_get_info()` free ≥ 90 % before building `LLM(gpu_memory_utilization=0.85)`; HF-path sampling and the G-probe GRPO steps run in their own later processes. **Qwen3.5-0.8B HF-only path, pre-decided:** if its CUDA parity (H-A5-2; CPU parity inconclusive, B5) or its GDN kernels (D15) fail in R0's Qwen3.5 arm (§4.7 steps 3, 7, 8) (vLLM's own GDN bug #34948 was fixed in 0.27.1, so it is not a 0.32 risk), the whole rung samples with HF generate: branch selection on a fixed 256-item shifted-val subset × K = 8, the full 1,084 items only for the selected branch (cost in §4.9). **Gates are unconstrained**; the D10-regex arm is a diagnostic only (H-A5-5). **H-A1-8** (legal ≥ 80 %, duplicates ≤ 10 % at T = 1 on both views): Qwen3.5 passes from 36 Mac steps (legal 82–85 %, B1 §2.6), so it needs no constrained decoding at GRPO start; LFM failed on shifted crops at Mac depth (71 % at 150 steps) and is decided by its Vast SFT, where G-format (legal ≥ 90 %) and G-signal (duplicates ≤ 10 %) on the shifted val imply H-A1-8's pass, so a C1 rung never needs it; the H-A5-5 trigger stays as in §6.4. Never `vllm_structured_outputs_regex` in GRPO without grammar-consistent trainer logprobs. K: 8 for format/geometry/signal gates, 16 for GED16/HM-IoU16, 64 for TV/retention on disagreement items and the I35 audit recalibration at Phase 3 entry (§6.4) | TRL + vLLM colocate for eval (breaks SmolVLM: TRL #6294 reports 855 vs 38 prompt tokens; A5's sandbox reproduction at 512 px gave 89 → 155 tokens and garbage output; needs the conflicting env). Constrained gates (hide format failures: base model gave `{"nodule":[]}` 8/8). transformers continuous batching (no image path). SGLang, TensorRT-LLM, LMDeploy, MLC (coverage) | A5 §0, §1, §3, §6, §8; B5 Q2: prompt ids 32/32 for LFM, Qwen3.5, SmolVLM on both datasets; fp32 logprob parity LFM median 1.0e-7, SmolVLM 4.0e-6 nats; Qwen3.5 inconclusive on CPU |
| D14 | Reward serving | **Preload the RL caches of the datasets that process scores, and only in processes that score** (`select_ckpt`, G-probe, `k3_harness`, Phase 3 GRPO), never in `train_sft` (SFT needs no caches beyond the voting indices stored in the records, M2). LIDC-only: 0.33 GB array bytes (B4 §1); RIGA adds 3.35 GB (Linux RSS ~4.2–4.5 GB) only on RIGA runs. Load into a dict **before model/CUDA init**; `gc.freeze()`; `ProcessPoolExecutor(4, mp_context=fork, initializer=cv2.setNumThreads(0))` warmed with 4 no-op tasks so all workers fork immediately; jobs carry `CacheRef(key)` (pickles 211 KB per RIGA batch vs 96.2 MB with caches); `TRLReward(caches=serve.loader(root), pool=pool)` unchanged. A dead worker raises `BrokenProcessPool` → fail loudly, resume from checkpoint. **C++ batch scorer deferred**; build it (B4 §4 design, ~600–900 lines + differential tests) if, on a quiet Vast box, RIGA-only pool-4 paired wall > 125 ms or speed-up < 3×, or Phase 3 needs larger groups, or < 4 cores are spare beside the dataloader | LRU-64 `cache_loader` (32/32 misses per RIGA step: 4.79 s vs 0.65 s preloaded). Threads (0.8–1.2×). `multiprocessing.Pool` (hangs forever on a killed worker). Forkserver/spawn with a TRL script as `__main__` (re-imports torch: macOS crash 3/3, Linux pool start 49.7 s). GPU scoring (exactness, launch cost) | B4 §1–5; A6 §0.4–0.5, §4 |
| D15 | Environments | **Phase 2: one image, two venvs**: `/opt/venv` (repo pins: torch 2.14.1 cu130, transformers 5.20.1, TRL 1.16.0, numpy 2.5.4, numba 0.68.0) and `/opt/vllm` (`vllm==0.32.0` lock: torch 2.13.0, transformers 5.19.0, numpy 2.3.5, numba 0.65.0, xgrammar 0.2.8). **Qwen3.5's GDN kernels in `/opt/venv`** (without them transformers runs Gated-DeltaNet and its causal conv on the reference PyTorch path, "correct but much slower"; B1 §2.7): transformers 5.20.1 resolves each of `chunk_gated_delta_rule`, `fused_recurrent_gated_delta_rule`, `causal_conv1d_fn`, `causal_conv1d_update` as (1) HF `kernels` hub builds when `use_kernels=True` (`kernels-community/fla`, `kernels-community/mamba-ssm`), else (2) the packages `fla` / `causal_conv1d`, else (3) the reference path (`integrations/hub_kernels.py`, local source read 2026-10-20). The image installs `flash-linear-attention` (Triton, no compile step; uses torch 2.14.1's Triton) and `causal-conv1d` (a compiled CUDA extension: needs a prebuilt wheel for torch 2.14.1 / cu130 / CPython 3.12, because the slim image has no nvcc), with the `kernels` hub route as the fallback for either. **UNVERIFIED for these pins**: B5 tested only vLLM's CPU wheel (B5 Q2) and no lab installed either package; M7's image build records which route each function resolves to, R0 step 3 asserts it on the GPU, and the vLLM sampler has its own GDN kernels. **Phase 3: one re-pinned colocate venv** (TRL 1.16, vLLM 0.32.0, transformers 5.19, torch 2.13.0) with uv overrides **numba 0.68.0 and numpy 2.5.4** written in the lock; it passed imports, the fast suite (666/1, the 1 a 1-ms wall-clock flake) and `xarch_hashes --strict` 26/26. Re-run the fast suite and strict hashes on every re-pin; expect TRL's out-of-range warning until TRL declares vLLM ≥ 0.32 (PR #7799) | vLLM in the pinned `.venv` (torch/transformers/numba pins conflict). vLLM ≤ 0.31 as A3 §6.1 proposed (0.31 caps transformers < 5.18; vision assistant-only loss needs ≥ 5.18). Clean re-pin to numpy 2.3.5 / numba 0.65 (fingerprint and caches change → `prepare_rl` re-run). **Server mode on one GPU** as A5 §8.4(a) proposed (TRL 1.16's client is NCCL-only, no IPC path; same device as two NCCL ranks is unsupported; PR #7799 KeyError on vLLM 0.32's `/inference/v1/generate`). SGLang (no TRL, no SmolVLM) | B5 Q1, Q3; A5 §0.2; A3 §6.1 |
| D16 | Data governance | **LIDC (CC BY 3.0) on any verified host. RIGA on a rented host only with written permission (depositor + MESSIDOR partners) AND only on Secure Cloud** (`datacenter=True`, ISO 27001): encrypted `bundle-riga.tar.age` decrypted into `/dev/shm`; **key delivery:** the age identity lives only in the Mac keychain and is streamed over the direct SSH session's stdin into the decrypt command (`age -d -i /dev/stdin bundle-riga.tar.age \| tar -x -C /dev/shm/riga`), never as an env var, template field, onstart argument or file on the host; the ledger records the recipient (public) key only; exact invocation tested in M7; no Cloud Sync, destroy on completion, ledger row per session; `sync_up` refuses it to any offer without `datacenter=True` (H-GOV). Until permission: RIGA gates are not judged (`phase2_entry(..., datasets=("lidc",))`); Mac MPS remains the last-resort RIGA path. **RIGA-free K3 harness** (`k3_harness`, synthetic RIGA-like items) on community hosts. Never publish RIGA-derived data; Phase 9 release covers LIDC only. Register at ADCIS for MESSIDOR | RIGA on community hosts even encrypted (host root, Vast backups outlive instances; A3 V12, V14). Third-party RIGA mirrors. Option D split (only 290/750 images non-MESSIDOR) | RESULTS I44 (permission AND Secure Cloud); A3 §5; A4 §5; B4 §3 (synthetic RIGA: CPU-time ratio 0.972 [0.935, 1.016]; digest `e31d8f1e…` on arm64 and x86) |
| D17 | Baselines | Same split, same `reward.score_mask` path, and **every row labelled with its view; the gate-relevant rows use the same shifted val as the VLM gate** (and the rescaled view where reported): (1)–(3) are inference-only and run on the **Mac CPU** before V1, and also compute the **per-item centred win share over the prior for MedSAM-oracle and SAM 2.1-tiny** (K = 1 each; the reference for G-prior (a), §5.1) (A4 §3.3: ~3.5 h CPU for all frames; $0). (1) **image-blind prior** (disc at frame centre, train-median diameter per S; mandatory, gating) on centred, shifted and rescaled views; (2) **MedSAM ViT-B re-run**, official min-max preprocessing, on the centred **and shifted** views, box from that view's union (+1 px LIDC, +5 px RIGA: still an oracle box, now on the matching view) plus a jittered-box row (0–20 px at 1024 scale); (3) **SAM 2.1 hiera-tiny**, box, same views (Apache-2.0; H-A4-10); (4) **nnU-Net v2 2D** (training removed from Phase 2, Phase 9 baseline, §6.4), label-sampling and majority-vote, 1 fold, 250 epochs, **trained on the jittered train view (D10 shifts + scale jitter) and evaluated on the shifted val**, never on centred frames (it would inherit the S5 leak); (5) MedSAM2 optional, research-use weights only. Published Dice converted with IoU = D/(2 − D); our LIDC numbers never share a column with Kohl v1/v2 | Phase 0 MedSAM numbers (LIDC 0.33 was a crop/scoring artefact; re-run gives 0.78). SAM 3 / Medical SAM3 (11.9 % external Dice text-only; 3.45 GB), BiomedParse (trained on LIDC: contaminated; NC), LISA/GLaMM/MedPLIB (LIDC Dice ≤ 0.006), MedSAM-Agent (inconsistent) | A4 §0.4, §3, S7 (MedSAM oracle box: LIDC 0.780, disc 0.941, cup 0.847; box-ellipse null 0.602) |
| D18 | Entry gates | Revised in §5.1: geometry on **shifted** val and on a **rescaled** val (shift + S moved to a neighbouring divisor); **paired win over the image-blind prior** on the shifted and rescaled views (the centred win share, G-prior (a), is report-only since revision 7); precision probe re-specified (D19); unanimous-topology gate **replaced** (not added to) in `eval/gates.py` by the synthetic topology probe (D20) with a hole-recall floor; **G-probe**: a 20-step GRPO run must improve the checkpoint; "IoU near MedSAM" reported, not gated (MedSAM gets an oracle box the policy does not; A4 §3.2) | RESULTS §8 gates as written (geometry passable image-blind on centred RIGA frames; topology gate vacuous: 0 unanimous non-simple val items; ≤ 2 px probe looser than the readers) | A4 §2.5, S5–S6; B3 §1–2 |
| D19 | Precision probe | Synthetic shapes, **generator 2** (elongation E = 0.2, amplitude A = 0.15; calibrated: centroid-disc IoU 0.822 vs real 0.806; σ 0.5 → IoU 0.915 / legal 96.4 %; σ 1 → 0.816 / 85.4 %), LIDC diameter distribution, uncentred frames (window + scale jitter), real LIDC **train** lung patches. Protocol per rung (revision 3): **zero-shot** first, the screening checkpoint on 1,000 held-out probe items (inference only); the capacity run (probe-only SFT, 20k items, one-LIDC-epoch budget, ~0.3 GPU-h on a 4090, then the same 1,000 items) only if zero-shot fails; **patient-disjoint backgrounds**: LIDC train patients are split by a seeded patient hash into probe-train patches (the 20k) and probe-test patches (the 1,000), tested in M4 (0 shared patients); the topology probe (D20) uses val patients only; greedy + K = 8 at T = 1. Pass on 6–16 px objects: legal ≥ 95 %, median ASSD ≤ 0.6 px, median IoU ≥ 0.80, median centroid error ≤ 0.5 px, **paired per-item IoU win over the equal-area centroid disc ≥ 60 %**. Topology moves to D20 | RESULTS §8 "≤ 2 px boundary error" (readers' ASSD 0.93 px; σ = 2 px gives LIDC IoU 0.63 and 50.9 % legal; A4 S2–S4). Probe on real LIDC (rater noise dominates). Generator 1 (failed calibration: disc IoU 0.771, σ = 1 legal 69 %) | A4 §6; B3 P2 (exact, brute-force-equal, 1-px shift detected 100 %, discrimination PASS; disc passes every geometric criterion, hence the paired-win rule: σ ≤ 0.5 wins 78–100 %, σ ≥ 0.75 / shift / disc ≤ 49 %) |
| D20 | Topology probe (B3 Addendum A) | **Synthetic**, independent items on real LIDC **val** lung patches: **120 one-hole items per hole-area bin (≤ 10 / 11–30 / > 30 px; 360 total), 160 two-part items (descriptive only: A0 keeps parts), 160 simple controls**. Outlines generator 2 (E = 0.2, A = 0.15). Each hole jointly resamples one real LIDC **train** hole record (holed part's diameter, hole area, centroid offset, wall to the outside); each two-part item resamples one train multi-part record (diameter, part-area ratio, gap); every hole ≥ the size rule (signature (1,) raw and at `min_area(\|GT\|)`; the rule max(5 px, ⌈0.5 %·\|U\|⌉) lives in the reward code, `objdel.checker` `min_area`, P1 PLAN §4.2, not in a note). **Real hole geometry** (train, voting raters, size rule; `hole_census.json`): 285 holes in 241 holed targets; area p10 / p50 / p90 = 5 / **12** / 49 px; bins ≤ 10 / 11–30 / > 30 px = **44 / 36 / 20 %**; hole = 3.6 % of its part (p50); wall p50 1.24 px (p10 0.41); holed part ≈ 24 px across; **75 % of holed targets come from items where only one rater drew a hole; no hole is unanimous**. **Primary metric (C2): per-item hole hit averaged over K = 8 samples at T = 1 and over the A0′ seeds (3 since revision 4; B3 A5 specified 2); pooled hole recall = train-share-weighted mean of the bins (0.439 / 0.361 / 0.200), stratified bootstrap of items within bins** (±3.1 pts at recall 0.5, item difficulty σ_u = 1.5; greedy ±4.9, any-of-8 ±4.4 and saturating at 0.99 for σ = 0.5). Greedy and any-of-8 reported, never deciding. **Generator fidelity (P5.2): failed narrowly** as pre-registered: the final generator (v3e) passed KS p ≥ 0.01 and bin shares ±5 pts on fresh seed 8 but not on fresh seed 7 (lowest p 0.004, hole area); the residual is a ~1 px area under-draw of exact 5-px holes in the smallest bin (probe bin medians 6 / 17 / 48 px vs real 5 / 17 / 49) and a slightly less central low-offset tail (p10 0.25 vs 0.18); wall, fraction, part ratio and gap match on every seed (p ≥ 0.11). Iteration stopped to avoid tuning to the test. **Why it remains fit:** the probe is stratified by area bin, so the bias acts only within the ≤ 10 px bin (making those holes ~1 px larger, i.e. slightly easier), C2 also requires every bin's point estimate ≥ 0.30, and the smallest-bin number is reported with this caveat. Every train non-simple item stays in SFT. The 44-item val non-simple slice is reported descriptively (its holes are smaller, 54 % ≤ 10 px, mostly single-rater); its matching quantity is the K-sample hole rate vs the rater hole fraction (TV), not recall | ~300 held-out real items (H-A1-7: no real design reaches ±6 pts at ICC 0.15–0.19; removes 55–100 % of train non-simple targets). B3's first probe (hole geometry uncalibrated; JUDGE_R1 MUST-FIX 1). Any-of-8 as primary (saturates, rewards random hole stamping). Paired A0′ − A0 test (equals A0′'s recall; adds nothing) | B3 P1 (census, cluster CIs); B3 Addendum A1 (real geometry), A2 (fidelity), A3 (power, `power_c2.json`), A4 (sanity: exact 1.0, A0 0 / 0 / 0, centroid disc 0, σ = 0.5 px mean-over-8 0.65 / 0.81 / 0.90, pooled 0.76 [0.73, 0.78], 0 false non-simple), A5 (metric, seeds) |
| D21 | C++ vs Python | Per module, §3.4. Python stays everywhere in Phase 2 except polycore (P1-D1/D14); numba stays for the simplifier; the only C++ candidate is the D14 batch scorer behind its trigger; an exact FH EDT in C++ only if caches are rebuilt online (H-8) | GPU raster/EDT/CCL (exactness, launch overhead, GPU busy); free-threaded 3.14t (OpenCV abi3 only); piecemeal numba/C++ ports (1.6× serial, no thread scaling) | A6 §0, §4, §6; B4 §4 |
| D22 | GPUs | Five cards (§4.8): RTX 4090 (default screening), RTX 5090 (best DLPerf/$, most Secure Cloud supply: the RIGA card), RTX PRO 6000 Blackwell 96 GB (Qwen3-VL-2B; Phase 3 colocate), A100 SXM4 80 GB (mature fallback), RTX 3090 (cheapest ops/smoke box) | L40S / RTX 6000 Ada / A6000 (1–6 offers, no Secure Cloud); H100/H200/B200 ($2.2–10.8/h, 1–5 offers, wall clock does not bind in Phase 2) | A3 §3 |
| D23 | Rental policy | On-demand, verified, single GPU, reliability > 0.98, ssh `--direct`, image pinned by digest, **destroy never stop**, Mac pulls every artifact, no cloud credentials on hosts, scoped API key, prepaid credit (§4) | Reserved (month prepay), interruptible for decisive runs (outbid → stopped, storage bills), Cloud Sync (credentials on host), volumes (machine-local), Jupyter mode, serverless | A3 §1 (V1–V21), §2, §7, §10 |
| D24 | Tracking, config, determinism | trackio (`report_to="trackio"`) plus run-dir JSON (`trainer_state.json`, gate reports, resolved config, git SHA + dirty-diff hash, image digest, host id, driver, prepared-manifest hash, `pip freeze`); overrides in YAML via `TrlParser`. No `full_determinism=True` (forces `CUDA_LAUNCH_BLOCKING=1`); deterministic mode (`CUBLAS_WORKSPACE_CONFIG=:4096:8`, `use_deterministic_algorithms(True, warn_only=True)`) for the decisive runs only if E7 shows overhead ≤ 15 % | W&B (only if the owner already uses it), Hydra, TensorBoard/MLflow | A2 §5 |
| D25 | Metrics | Add to `objdel.eval` now: **HM-IoU (K = 16)**, **NSD@τ** (τ = 2 px LIDC, 9 px disc, 16 px cup at frame resolution). **Removed** (§6.4): HD95, ASSD and centroid error in `objdel.eval` (the probes keep their own, B3 `probe_lib`). G-format also reports truncation share, mean / p95 completion tokens and held-out NLL per branch. Report GED16 with the sample count, D_max, soft Dice (5 thresholds), human LOO references (LIDC 0.768, disc 0.945, cup 0.815 IoU_tc). `eval_run` passes the prepared caches into `metrics` (1.3–5× faster) | Boundary IoU (equals mask IoU under ~30 px) | A4 §2.1, §2.4, S2; A6 §1c |
| D26 | Image resolution | Keep P1-D12: LIDC upscaled to 512 px (Pillow bicubic, I32); LFM/Qwen processor defaults (≤ 256 image tokens); SmolVLM `longest_edge` 512 LIDC / 1024 RIGA, RIGA frames pre-resized offline if SmolVLM stays (H-7: −40 % collate). RIGA 512-px upscale (H-A1-5) deferred to the joint run | — | A1 §1.3 (250 px RIGA frames get 64 tokens); A2 §4.4; A6 §1a |
| D27 | Repository layout | Library changes in `src/objdel` with tests (§7); Phase 2 scripts in `experiments/phase2/` run as modules with `-I`; Vast scripts and the tested governance/ledger helper in `experiments/phase2/vast/`; lab code under `lab/` stays evidence and is promoted, not imported, by production scripts (the `k3_harness` promotion included) | — | CLAUDE.md; README §Rules |
| D28 | Checker-filtered rejection-sampling SFT (fallback) | **Removed from the critical path (revision 3);** kept as triggered E10 (§6.4): only with the owner's approval, when a rung passes parse but legal stays in [85 %, 90 %) at every branch (the legality knife-edge, A4 S4): sample K = 8 at T = 1 per LIDC train item from the selected checkpoint, keep samples that are legal with IoU_tc ≥ 0.5 by the exact reward path, deduplicate by raster hash, keep ≤ 2 per item, mix 1:1 with rater targets, one short-SFT pass at decay LR; judged by the same gates plus G-signal and G-probe (diversity guard). Records that the Phase 3 start point then differs from plain SFT (Phase 3's best-of-8 control must start from the same checkpoint) | Doing it by default (it trains toward the model's own modes and can cut diversity; it overlaps Phase 3's best-of-8 / DPO controls) | RFT (2308.01825: rejection-sampled correct paths as SFT data, LLaMA-7B GSM8K 35.9 → 49.3 %; more distinct paths help more); ReST-EM (2312.06585: sample, filter by binary feedback, fine-tune, repeat; beats human-data-only fine-tuning); abstracts re-read 2026-10-20; ROADMAP Phase 3 best-of-8 control |

## 3. Architecture

### 3.1 System and data flow (Mac → image / Vast box → checkpoints and eval back to the Mac)

```
 Mac (arm64; owner; .venv pinned)                  GHCR (private)                  Vast host (x86-64, driver ≥ 580, verified)
 ────────────────────────────────                  ──────────────                  ──────────────────────────────────────────
 data/store (37fe7665…) ─► prepare_rl ─► data/prepared (3747acbc…)
                      └─► build_views ─► data/views/<view>/  (per-epoch row specs: dihedral op, shift, S, HU window;
                          (M2, offline)   certified target texts; per-view target_set_hash; context crops once per item)
                      └─► build_eval_views (gate views: fixed-seed shifted val + rescaled val; targets + caches rebuilt)
                      └─► probes: precision (gen. 2, train patches) + topology (480 items, val patches)
            │
            ├─ bundle-lidc.tar + SHA manifest ─────── rsync over direct SSH ────────► /workspace/data/lidc   (verified vs manifest)
            ├─ bundle-riga.tar.age ── only: written permission ∧ datacenter=True ──► /dev/shm/riga          (H-GOV refuses otherwise)
            ├─ git ls-files tarball (no data/, .venv, build/) ── rsync ─────────────► /workspace/src  (uv pip install -e, ~20 s)
            └─ uv pip freeze ─► Dockerfile.phase2 ─buildx amd64─► objdel-gpu:p2-<pinhash> ─pull─►  /opt/venv + /opt/vllm
 offers.sh / launch.sh (scoped key, keychain) ───────── vastai API ─────────────► instance (ssh --direct, onstart watchdog, label run_id)
 pull.sh (launchd, every 15–30 min) ◄──────────────── rsync ◄──────────────── /workspace/runs/<run_id>/
      ▼                                                                         (trainer_state, logs, ckpt last + best + WSD branches,
 runs/<run_id>/ on the Mac ─► eval_run / phase2_entry / reports (-I)             samples JSONL, k3 / gate JSON, manifest)
 ledger.csv ◄── offer JSON, show invoices-v1, destroy confirmation ◄── teardown.sh (final pull → destroy → verify absent)
```

### 3.2 Inside one training box

```
 container (image objdel-gpu:p2-<pinhash>; code at /workspace/src; run under `timeout <cap>` in tmux; heartbeat file)
 GPU owners run strictly one after another (run.sh chains them; each starts only after the previous process exited):

 ① /opt/venv: train_sft (python -I -m experiments.phase2.train_sft --config run.yaml)          [GPU: trainer only]
      start: verify bundle SHA (no cache preload, no reward pool: SFT never scores)
      SFTTrainer ◄── DataLoader (6 persistent workers): row spec → context crop → D4 / shift / HU window at frame
      (WSD, group_by_length)    resolution → Pillow upscale → processor; targets are precomputed texts (D9, D10)
      saves: rolling ckpt (optimizer, every ~10 min, limit 3; G-vision check there, abort only if (a) and (b) both fail, D7)
             + WSD decay-branch ckpts and the stable-phase ckpts at 0.5 / 1 / 2 epochs (model-only, kept) → EXITS
          │
          ▼
 ② /opt/venv: select_ckpt (D12)                                                              [GPU: never touched]
      preload LIDC caches (0.33 GB; RIGA only on RIGA runs) → gc.freeze() → ProcessPoolExecutor(4, fork) warmed
      per branch: export prompts JSONL ─────────────────────────────┐
                  launch ③ as a child (asserts torch.cuda.is_initialized() is False here)
                  ◄── completions JSONL ◄───────────────────────┐   │
      score: reward path + pool (CacheRef) → eval_run → summarize → gates JSON (G-format/geometry/prior/signal/vision)
             + G-vision forward passes on the CPU (B1 protocol) + B2 signal panel per branch → D12 rule R-B
          │                                                     │   ▼  JSONL only (no objdel, no shared pins)
          │      ③ /opt/vllm (vLLM 0.32.0, torch 2.13.0): sample_vllm                 [GPU: sampler only]
          │         asserts mem_get_info() free ≥ 90 % → LLM(merged ckpt, processor saved with our kwargs,
          │         gpu_memory_utilization 0.85, bf16, prefix caching) → n = 8, explicit SamplingParams → EXITS
          ▼             (Qwen3.5 HF-only path: ③ is objdel.train.sample in its own /opt/venv process instead)
 ④ /opt/venv: g_probe (selected branch; the D12(c) repair branch only on failure): preload caches → pool → GRPOTrainer 20 steps, use_vllm=False,
      eval_on_start + eval at step 20 on the fixed G-probe item set → G-probe JSON → EXITS    [GPU: GRPO only]

 CPU budget (≥ 12 effective cores): ① 1 trainer + 6 dataloader workers; ②/④ 1 + 4 reward workers; ③ vLLM
```

### 3.3 Phase 2 run graph (the only Vast work; current at revision 8)

```
 V0  H-WATCH dry run on the cheapest 3090 (cap 10 min, no job; gone at 11 min)                            ≈ $0.06
        ▼
 V1  one 4090 session, LIDC only; processes chained one at a time (§3.2); Mac pulls every 15–30 min
   R0 box checks: image_check ─► xarch_hashes ─► k3_harness (idle) ─► sampler parity on a few items
      ─► timing prelude (4 × 30 steps from R1's init, discarded: workers 0 / 6, deterministic on / off)
      ─► Qwen3.5 arm (30 steps, discarded: s/step, peak memory, GDN kernel route, parity, G-vision base refs)
        │ blocking failure (non-accepted hash diff, K3 fail) ─► stop and fix (C++ scorer trigger, D14)
        ▼
   R1 screening, rung k = smallest remaining (D1): bundle (§6.0), A0′, seed 1, WSD 2 epochs + branches 0.5/1/2
        │ G-vision (a) and (b) both fail at a rolling ckpt ─► abort; D7 repair: rerun frozen (once); one alone = warn
      ─► train_sft EXITS ─► select_ckpt: shifted val K = 8 per branch ─► D12 R-B; selected branch: K = 16,
      rescaled val, precision probe (zero-shot), topology probe; k3_harness once under load ─► §5.1 gates
        │ gate fail ─► first failure in Phase 2 only: R-diag, one leave-one-out run (§6.3)
        │               ├─ passes ─► rung k passes with bundle minus that component (kept for C2)
        │               └─ fails  ─► k := next rung (LFM ─► Qwen3.5-0.8B ─► Qwen3-VL-2B*), full bundle
        ▼ gates pass
   R2 G-probe on the selected checkpoint (probe LR, step-0 double sample = noise floor) ── not PASS ─► D12(c) repair
      (next-earlier eligible branch; WiSE-FT scored, report-only); still not PASS ─► C1 with a G-probe warning (no climb)
        ▼ pass  ══►  C1 = rung k.  STOP the ladder.
   R3–R5 C2 on rung k: A0′ seeds 2 and 3, A0 seed 1 (same ×4 rows, holes emptied), each to the depth of the checkpoint with the final
      G-probe verdict; A0′ seed 1 = that checkpoint ─► topology probe ─► C2 3-seed rule (§5.2)
        │ C2 partial (pooled [0.30, 0.50), or ≥ 0.50 with a failed clause) ─► R-x8: one ×8 A0′ seed*  │ C2 fail (pooled < 0.30) ─► R-ctrl (Qwen3-VL-2B)*
        ▼
   teardown (final pull, destroy, ledger)                    * owner approval; R-ctrl on the PRO 6000
 Mac (no GPU): image-blind prior, MedSAM, SAM 2.1-tiny on centred / shifted / rescaled val (D17); gate reports re-run
                from the pulled samples
 Phase 2 exit: §5.1 per rung run, C1 + C2 verdicts, spend ≤ hard stop. Everything else: §6.4 triggers.
```

### 3.4 C++ / Python split in Phase 2 (D21)

| Module / step | Language | Decision and trigger | Evidence |
|---|---|---|---|
| Raster, legality, signed area | C++ polycore (unchanged) | keep | P1-D1, P1-D14 |
| SFT rows, collation, image ops | Python (TRL collator, Pillow, numpy) in 6 workers | keep; no C++ (CPU 10–25 ms/row is hidden by ≥ 4 workers on every card) | A6 §1a, §2 |
| SmolVLM RIGA Lanczos resize | offline pre-resize (if SmolVLM stays) | −40 % collate (H-7) | A6 §1a (181 of 407 ms) |
| Augmented targets | numba simplifier, offline, crop-before-trace | keep numba; never online | B3 P3–P4; A6 §1b |
| Expert caches (EDT) | scipy, preloaded | keep (deterministic, I30); C++ FH EDT only if caches are rebuilt online (H-8) | A6 §4; B4 §1 |
| Reward scoring | Python + OpenCV + polycore, fork pool of 4 | keep; C++ batch scorer only on the D14 trigger | B4 §2, §4 |
| Eval metrics, probes | Python + numpy | keep (K = 64 full val ≈ 17 s with loaded caches) | A6 §1c |
| `geometry.parse` | regex + json | keep (7 % of reward time) | A6 §1d |
| vLLM sampler | separate venv, JSONL bridge | keep | B5 Q2 |

### 3.5 Package changes (summary; tasks in §7)

```
src/objdel/
  targets/frame.py      + frame(item, cfg, offset=(dx,dy), side=S, window=(L,W))   shift / scale / HU-window views
  targets/build.py      + crop-before-trace in _rings (hash-preserving)
  targets/augment.py    NEW: D4 ops on masks and images at frame resolution; row-spec hashing; view builder helpers
  train/configs.py      + qwen3.5-0.8b entry (non-thinking template, bf16, MAX_COMPLETION 1,024); explicit sampling
                          overrides per model; lora() LM-only (vision tower + projector excluded, was all-linear incl.
                          vision); vision + projector LR groups (0.1×) for full FT; WSD, group_by_length, workers, grad-ckpt
  train/dataset.py      + sft_rows over views (op per rater row, ×k non-simple, A0 transform, mixing), length column,
                          voting indices in the record; cache_loader(maxsize=None) preload mode
  train/serve.py        NEW: preload, CacheRef, fork_pool (warm, gc.freeze, initializer)
  train/sample.py       + prompt export / completion import JSONL contract for the vLLM bridge
  eval/metrics.py       + HM-IoU, NSD@τ (HD95, ASSD, centroid deferred); truncation and token stats in summarize
  eval/signal.py        NEW: B2 signal panel, variance split (sd_geom), D12 rule R-B, per-item Spearman
  eval/vision.py        NEW: G-vision: projector cosine distance between images, image-use TV (from lab/B1)
  eval/probes.py        NEW: precision probe (generator 2) and topology probe: generators, scorers, pass rules
  eval/gates.py         revised phase2_entry: shifted + rescaled views, prior paired win, probes, G-probe;
                        the unanimous-topology gate is REMOVED, not kept beside the probe
experiments/phase2/     train_sft, select_ckpt, sample_vllm, g_probe, k3_harness, probes, baselines, c2_analyze,
                        build_views, build_eval_views (shifted + rescaled val); vast/ (scripts + tested policy.py); Dockerfile.phase2/phase3
```

## 4. Vast operations

### 4.1 Account (owner, before V0; A3 §2)

1. Verify the email (lifts the new-account `spend_rate_limit`; V15). Add the Mac's ed25519 key.
2. Scoped API key for scripts: `misc`, `instance_read`, `instance_write`, `user_read`; no billing/team scopes; macOS
   keychain. Full-power key stays in the browser (V18, G3).
3. Prepaid credit in $25 tranches; low-credit email/webhook at $10 (G1, G2). Card vs no card: owner decision (§8.3).
4. Private GHCR package for the image; read-only `read:packages` PAT used only as Vast `--login` (V10).
5. `.vast-allowlist`: offer filters (§4.2) and the governance tier per bundle (D16).

### 4.2 Instance type and filter

On-demand, verified, single GPU, ssh mode `--direct`, `--cancel-unavail`, container disk sized at creation (cannot
grow; V6): 60–80 GB for ≤ 0.8B models, 150 GB for Qwen3-VL-2B (A3 §6.4).

```
vastai search offers 'num_gpus=1 gpu_name=<GPU> gpu_ram<=<stock+2> verified=true rentable=true reliability>0.98
  cpu_cores_effective>=12 cpu_ram>=32 cuda_vers>=13.0 disk_space>=80 inet_down>=500 direct_port_count>=2
  duration>=1 dph<=<ceiling>' -o 'dph_total' --storage 80 --raw          # + datacenter=True for the RIGA bundle
```

Deviations from A3 §3.3 (recorded in §9): `cpu_cores_effective ≥ 12` (not 8), because a training box runs 1 trainer
+ 6 dataloader workers + 4 reward workers (D4, D14); snapshot medians are 16–34 vCPUs for all five cards, so supply
holds. `cpu_ram ≥ 32` for LIDC-only boxes (caches 0.33 GB, preloaded only by scoring processes, D14); `cpu_ram ≥ 48` for
RIGA boxes (RIGA caches ~4.2–4.5 GB Linux RSS, B4 §1) and the PRO 6000. `gpu_ram ≤ 26` for a stock 4090 (16 of 104 "4090" offers report 48 GB; A3 §3.2); `gpu_ram ≥ 79` for A100 80 GB.
Price ceilings per card: 1.5 × the snapshot on-demand median (§4.8).

### 4.3 The image

`Dockerfile.phase2` (promoted from `lab/B5/Dockerfile.phase2` + `lab/A3/Dockerfile.vast`): base
`python:3.12-slim-trixie` at the bench image's digest; no CUDA base image (the torch wheel carries CUDA 13.0 user-space;
the host injects the driver). `/opt/venv` = the Mac's `uv pip freeze` (same file as `x86.sh`), `/opt/vllm` = the
`vllm==0.32.0` lock; openssh-server, rsync, tmux, procps, build-essential, cmake, ninja. `/opt/venv` adds image-only extras for Qwen3.5 outside the Mac freeze (they cannot install on macOS): `flash-linear-attention`, `causal-conv1d`, `kernels` (D15; UNVERIFIED for torch 2.14.1 / cu130; nothing in `objdel` imports them, so `xarch_hashes` is unaffected). Code is never baked in.
Measured: `/opt/venv` layer 3.61 GB compressed / 6.4 GB on disk, sm_75/80/86/90/100/120 kernels, H-IMG PASS (A3 §6.3);
sampler stage 4.14 GB / 7.9 GB (B5 Q3); full image ≈ 7.5 GB compressed (estimated, B5). `Dockerfile.phase3` (D15) is
built at Phase 3 start only. Rebuilt only when pins change; tag `p2-<pinhash>`; pulled by digest.

### 4.4 Data and code transfer

- `bundle-lidc.tar`: `data/prepared` LIDC (images 1.3 GB, jsonl 8 MB, caches 42 MB; A3 §6.2) + `data/views` LIDC:
  int16 HU context crops (one npz per train and val item; the crop covers every shift and the neighbouring-divisor
  scale: ≈ 2 S_max per side, ≈ 200 px for the 85 % of items at S = 50 and ≈ 250 px at S = 100: ≈ 85 KB raw per item ×
  8,425 items ≈ 0.7 GB raw, **≈ 0.3–0.5 GB compressed**, an estimate M1 replaces; B3's 8-bit PNG crops were 17.6 MB but
  cannot carry window jitter) + target texts ~57 MB for 8 ops (B3 P3) + eval views + probes; SHA-256 manifest checked on
  the box before training (mismatch aborts).
- `bundle-riga.tar.age`: RIGA prepared + views (images 148 MB; RIGA context crops 129 MB; B3 P3); D16 only; decrypted
  with the identity streamed over the SSH session's stdin (D16 key delivery), plaintext only in `/dev/shm`.
- Code: `git ls-files` tarball, as `x86.sh`; `uv pip install -e .` on the box (~20 s on x86, M9).
- Models: HF Hub at a pinned revision SHA (recorded); no HF token for public models.
- Checkpoints: pulled by the Mac every 15–30 min; SmolVLM/LFM checkpoint with optimizer ~3–5 GB; Qwen-2B ~25 GB, so
  pull weights only plus the final optimizer state (A3 §6.2). Resume = rsync the checkpoint to a new box +
  `resume_from_checkpoint` (H-PULL).
- Nothing leaves the host except by the Mac's pull; off-Mac durability is pushed **from the Mac** (private HF repo or B2).

### 4.5 Scripts (`experiments/phase2/vast/`; A3 §6.5)

| Script | Does | Guardrail |
|---|---|---|
| `offers.sh <gpu> [--secure]` | §4.2 search → cheapest 5 → human picks → offer JSON to the ledger | price ceiling; `gpu_ram` bounds; `cuda_vers ≥ 13.0`; `datacenter=True` when the RIGA bundle is requested |
| `launch.sh <offer> <run_id>` | `create instance … --ssh --direct --cancel-unavail --label run_id --onstart-cmd <watchdog>`; poll `show instance --raw` ≤ 15 min | destroy on timeout / `exited` / `unknown` / `offline` (V19); refuses a 2nd concurrent instance without `--allow-parallel` |
| `sync_up.sh` | rsync code + allowed bundle; remote SHA verify | refuses RIGA to non-datacenter offers (H-GOV, decided by tested `policy.py`); refuses anything under `data/raw` |
| onstart watchdog | env → `/etc/environment` (V9); self-destroy via `CONTAINER_API_KEY` | uptime > cap, or GPU idle 30 min and `PULLED_OK` exists |
| `run.sh <config>` | `timeout -k 120 <cap> /opt/venv/bin/python -I -m experiments.phase2.<script> …` in tmux; heartbeat | wall-clock cap per config |
| `pull.sh <run_id>` | touch `PULLED_OK.next`; rsync `runs/<run_id>/` → Mac without `rolling/ decay/ stable/` (full checkpoints); checksums; rename to `PULLED_OK` (the pull's start time) | — |
| `push_hf.sh <run_id> <config>` | refuse a keychain HF token that is not fine-grained and scoped to `svenk029/objdel-sft` write only; ssh `hf_push` with the token on stdin: the selected branch + base LICENSE + model card, one commit on repo branch `<output_dir>`; refused on the box for a public repo, no selected branch or a non-LIDC run | licence (RIGA, LFM §4) |
| `teardown.sh <run_id>` | final pull (3 tries) → always `vastai destroy instance <id> -y` → verify absent → invoices to ledger | — |
| `sweep_orphans.sh` (launchd, hourly) | destroy instances not in `active_runs.json` or past their cap; alert | catches forgotten instances |

### 4.6 Cost guardrails (A3 §7)

G1 prepaid credit is the hard cap ($25 tranches, owner approves each); G2 low-credit alert at $10; G3 scoped key;
G4 destroy never stop; G5 per-run `timeout` + self-destroying watchdog; G6 hourly orphan sweep; G7 price ceiling in
every search, interruptible only with an explicit `--bid_price` (V3) and only for idempotent jobs (eval sampling,
replicates); G8 poll loops with timeout + destroy; G9 ledger row per session (offer JSON, start/stop, invoices,
bundles sent, destroy confirmation; no secrets); **G10 Phase 2 alarm at $10 spent, hard stop at $15** (owner directive);
passing either needs the owner's written approval (§4.9).

### 4.7 The first session (V0 + V1; LIDC only)

V0 (cheapest 3090, ~10 min, ≈ $0.06): H-WATCH, launch with cap = 10 min and no job; pass = `show instances` empty at
11 min and charge ≤ 12 min. V1 (one RTX 4090, or a 3090 if the 4090 is above $0.70/h; 80 GB disk) runs R0 → R1 → R2 →
R3–R5 back to back on one box.

| # | Step (R0 unless stated) | Pass / record | Source |
|---|---|---|---|
| 1 | Preconditions: verified email, key, scoped key, $25 credit, alert at $10 | — | A3 §2 |
| 2 | `offers.sh RTX_4090` → `launch.sh`; record `lscpu`, `nproc`, cgroup `cpu.max`, driver, `nvidia-smi -q` | ledger row | A3 §8.2 |
| 3 | `image_check.sh`: `cuda.is_available()`; bf16 8192³ matmul ≥ 60 % of dense peak (4090: 165 TFLOPS); `objdel._polycore` and `/opt/vllm` import; Qwen3.5 GDN route per function (D15: hub / package / reference) recorded | all true; a reference route is not a stop, the Qwen3.5 arm (step 8) prices it | A3 §8.3; A1 §2 |
| 4 | `sync_up.sh` LIDC bundle; SHA verify 100 %; fixed G-probe item set and rescaled val hashes = Mac build | — | A3 §2; §5.1 |
| 5 | `xarch_hashes --out vast.json`; on the Mac `--compare build/x86/hashes-mac.json vast.json` (expect 22/26 + 4 accepted `splits.*`) and `--strict` vs the Docker x86 run (26/26). AVX-512 difference in `reward.advantages` → rerun with `NPY_DISABLE_CPU_FEATURES="X86_V4"` and decide | non-accepted difference blocks training | RESULTS §8; I34 |
| 6 | `experiments.phase2.k3_harness` idle (LIDC + synthetic RIGA; serial, 8 threads, fork pool 2/4/8; 5 batches); once more during R1 (CPU contention) | K3 ≤ 250 ms, K4 p99 ≤ 5 ms, pool hashes = serial; pool-4 RIGA wall (D14 trigger 125 ms / 3×) | B4 §3–5 |
| 7 | Parity on a few items, base weights of the rung's model: prompt ids on 32 LIDC val items; bf16 logprob parity on 8 items × n = 8; for LFM and Qwen3.5 both, in R0 (B5 §What this changes, item 4) | ids 100 %; median \|Δ\| ≤ 0.02 nats, p99 per-sequence \|ΣΔ\| ≤ 0.5; echoed T = 1, top_p 1, top_k off | A5 §10; B5 Q2 |
| 8 | Timing prelude: 4 × 30 steps from R1's init with the bundle config, discarded: workers 0 vs 6; `use_deterministic_algorithms(True)` on vs off (`CUBLAS_WORKSPACE_CONFIG=:4096:8` set for the session). **Qwen3.5 arm** (revision 5, D1): 1 × 30 steps of Qwen3.5 full FT with the bundle config, workers 6, discarded, capped at 0.3 h. G-vision base references for both models (D7, §5.1) | s/step, data-wait share, GPU util, peak memory (grad-ckpt off fits?), determinism overhead; projected R1 hours for LFM and for Qwen3.5 (D6 full FT vs LoRA from its peak memory) | A2 H-A2-1/9/10; A6 H-5; A3 H-STEP |
| 9 | **Stop rule:** if the projected R1 time exceeds 2× the §4.9 high estimate, stop and re-plan before training; Qwen3.5's projection is checked against its own §4.9 row before its R1 would start | — | §4.9 |
| 10 | R1, R2, R3–R5 per §6.1; R2 records the GRPO step split (generation / scoring / update) with `use_vllm=False` | gates §5.1; G-cost | RESULTS §8 |
| 11 | `pull.sh`, `teardown.sh`, ledger; measured s/step into `estimate_cost.py` (adding the Qwen3.5 entry it lacks) and ROADMAP §Compute (coordinator) | instance absent | A3 §8.8 |

CI on GitHub has already run green (run 39359008526), so I34's CI half is closed; its native-x86 half is steps 5–6.

### 4.8 The five GPUs

Prices: Vast public offer search, **snapshot 2026-10-19 ~24:00 UTC**, verified, ≥ 8 effective cores, reliability > 0.98,
CUDA ≥ 12.8, 100 GB storage included (`lab/A3/vast_snapshot.json`; A3 §3.2). Per-run hours and $ from
`lab/A3/estimate_cost.py` (2 epochs over LIDC 23,469 + RIGA 3,282 rows, MFU 0.10–0.30, +0.25 h eval; full FT except
Qwen-2B on ≤ 32 GB cards = LoRA). These are estimates; the R0 timing prelude replaces them.

| GPU | Memory / arch | Role (Phase 2; Phase 3–4) | On-demand $/h min / median (n) | Secure Cloud $/h (n) | LFM2.5-VL-450M run: h, $ | Qwen3-VL-2B run: h, $ |
|---|---|---|---|---|---|---|
| RTX 4090 | 24 GB, Ada 8.9 | default screening card (LFM, Qwen3.5); RESULTS cost gate; Phase 3 small-model GRPO | 0.51 / 0.66 (23) | 1.57 (1) | 0.66–1.49 h, $0.32–0.98 | LoRA 2.36–6.59 h, $1.14–4.35 |
| RTX 5090 | 32 GB, Blackwell 12.0 | same models with headroom; **the RIGA card** (most Secure Cloud supply); DLPerf 198 vs 97 | 0.63 / 1.02 (37) | 0.80–2.29 (9) | 0.58–1.23 h, $0.32–1.18 | LoRA 1.92–5.25 h, $1.07–5.04 |
| RTX PRO 6000 Blackwell | 96 GB, 12.0 | **Qwen3-VL-2B** screening; Phase 3 Qwen GRPO colocate on one card (H-PRO) | 1.36–1.44 / 1.56–1.60 (17 + 10) | 1.61–2.05 (2) | 0.39–0.66 h, $0.53–1.04 | 0.95–2.35 h, $1.29–3.68 |
| A100 SXM4 80 GB | 80 GB, Ampere 8.0 | fallback if the Blackwell stack misbehaves; Secure Cloud RIGA on 80 GB | 1.09 / 1.19 (7; `gpu_ram ≥ 79` filter; unfiltered A100 SXM4 0.48 / 0.92, n = 14) | 0.69–1.19 (4, 40 + 80 GB mixed) | 0.47–0.91 h, $0.51–1.08 | 1.37–3.61 h, $1.49–4.29 |
| RTX 3090 | 24 GB, Ampere 8.6 | V0 H-WATCH dry run; cheapest fallback box; Phase 8 CUDA dev | 0.20 / 0.32 (21) | 0.39–0.40 (2) | 1.21–3.14 h, $0.29–1.04 | LoRA 5.17–15.0 h, $1.24–4.95 |

Risks per card (A3 §3.4): 4090 modified 48 GB cards (filter); 5090/PRO 6000 youngest stack, FA2 not listed for sm_120
(SDPA used anyway, D4), PRO 6000 BF16 peak derived (UNVERIFIED); A100 no FP8, 40 GB variants mixed in; 3090 slowest.
Price ceilings (G7): 4090 $1.00, 5090 $1.50, PRO 6000 $2.40, A100 $1.80, 3090 $0.50.

### 4.9 Phase 2 budget (revisions 3–8, bottom-up; replaced by the R0 timing prelude)

**Per-run estimate and its uncertainty.** `estimate_cost.py` prices 2 epochs over LIDC + RIGA at MFU 0.10–0.30; for
LFM2.5-VL-450M on a 4090 that is 0.41–1.24 h of training (its total minus 0.25 h eval). A screening run is ≈ 2.5
LIDC-only epoch-equivalents (D5) and LIDC carries 86 % of LFM's tokens (A2 §1.2), so the training part scales by
2.5 × 0.86 / 2 = 1.075, and ×4 hole oversampling adds 3 × 414 rows (+5.3 %): **× 1.13 → 0.46–1.40 h**, plus ≈ 0.3 h of
evaluation (UNVERIFIED) (3 branches × 1,084 shifted-val items × K = 8 with
vLLM, then the selected branch's K = 16, rescaled val and both probes) → **R1(LFM) ≈ 0.76–1.70 h**. The range is ×3
from MFU alone; vLLM decode speed and the GRPO step time (assumed 15–70 s, G-probe 0.2–0.5 h) are UNVERIFIED. The
R0 prelude (§4.7 step 8) replaces the estimate before any long run, with a stop if it projects > 2× the high end.
Qwen3.5-0.8B is taken as 2× LFM: ≈ 1.8× the parameters (0.8B vs 450M), and B1 measured equal answer lengths (75–95
tokens per answer for both; Qwen3.5's 9.5 tokens/vertex vs LFM's 6.5 is offset by 8 vs 10 vertices; B1 §2.2–2.3,
§2.5), so A1's 1.5× tokenizer gap adds no answer-token surcharge. B1's Mac s/step (Qwen3.5 34–38 s vs LFM 12–16 s)
is not used: Qwen3.5 ran its Gated-DeltaNet layers on the reference PyTorch path (B1 §2.7). The 2× stays
**UNVERIFIED** until R0's Qwen3.5 arm measures it. Qwen3-VL-2B scales by 1.06 on a PRO 6000.

Prices: 4090 median $0.66/h (ceiling $1.00), PRO 6000 $1.60/h (§4.8). Baselines run on the Mac CPU ($0).

| Block | Card | GPU-h | $ at median (at ceiling) |
|---|---|---|---|
| V0 H-WATCH | 3090 | 0.2 | 0.06 |
| R0 box checks + timing prelude + Qwen3.5 arm (revision 5: +0.2 h) | 4090 | 0.9 | 0.59 (0.90) |
| R1 LFM screening (C1) | 4090 | 0.76–1.70 | 0.50–1.12 (0.76–1.70) |
| R2 G-probe | 4090 | 0.2–0.5 | 0.13–0.33 (0.2–0.5) |
| R3–R5 C2: A0′ seeds 2–3 + A0 seed 1 to the depth of the checkpoint with the final G-probe verdict (≤ 2.3 epoch-eq each: 0.43–1.30 h with ×4 rows) + probe evals | 4090 | 1.6–4.2 | 1.06–2.77 (1.6–4.2) |
| idle while gates are computed on the box | 4090 | 0.2 | 0.13 (0.2) |
| **Critical path (LFM passes, scenario A)** | | **3.9–7.7** | **2.5–5.0 (≤ 7.6)** |
| Scenario A at MFU 0.05 (no ≤ 2B VLM has a measured s/step, A6 §7; R1 ≈ 3.1 h stays under the R0 stop rule's 3.4 h) | 4090 | ≈ 13 | ≈ 8.6 |
| *Conditional* R-diag (one leave-one-out run, §6.3) | 4090 | 0.5–1.4 | 0.3–0.9 |
| *Conditional* WiSE-FT inference check (D12(c) repair fires; revision 6) | 4090 | ≈ 0.05 | ≈ 0.03 |
| *Conditional* vision repair (D7): the aborted part of a run whose G-vision (a) and (b) both failed; the frozen rerun replaces it | 4090 | ≤ one R1 | ≤ 0.5–1.1 |
| *Conditional* Qwen3.5-0.8B rung, R1 + R2 (2 × 0.76–1.70 + 0.2–0.5) | 4090 | 1.7–3.9 (+0.7–3 if HF-only sampler, D13) | 1.1–2.6 (+0.5–2) |
| *Conditional* C2 on Qwen3.5 instead of LFM (extra) | 4090 | +1.3–3.9 | +0.9–2.6 |
| *Conditional* Qwen3-VL-2B rung, R1 + R2 | PRO 6000 | 1.5–3.3 | 2.4–5.3 |
| *Conditional* C2 partial: R-x8, one ×8 A0′ seed (E1; §5.2) | 4090 | 0.5–1.4 | 0.3–0.9 |
| *Conditional* C2 fail: R-ctrl, one Qwen3-VL-2B A0′ seed | PRO 6000 | 0.8–2.3 | 1.3–3.7 |

Scenarios: **A** (LFM passes) 3.9–7.7 GPU-h, $2.5–5.0. **B** (LFM fails, one diagnostic, Qwen3.5 passes, C2 on
Qwen3.5) 7.2–16.3 GPU-h (+0.7–3 HF-only), $4.7–10.7 (+0.5–2): within $15, but above ~12 GPU-h at the high end.
**C** (Qwen3-VL-2B needed; C2 on it) 8.0–18.0 GPU-h, $8.9–21.4, built up as V0 + R0 + idle, LFM R1 (0.76–1.70) →
R-diag (0.5–1.4) → Qwen3.5 R1 only (1.5–3.4; no R2 after a gate failure) → Qwen3-VL-2B R1 + R2 (1.5–3.3, PRO 6000) →
C2 as 3 runs at the R-ctrl row's cost (2.4–6.9, PRO 6000); the 2B C2's probe evals (≈ 0.3 h on the PRO 6000, ≈ $0.5)
are not in the R-ctrl row. Scenario C can pass the hard stop, so the 2B rung starts only
with the owner's approval of a further tranche. **GPU-hours are likely 1.5–2× low** (JUDGE_R2: MFU 0.10 is an
unmeasured floor; the eval block is unverified): the MFU 0.05 row shows scenario A at ≈ 13 GPU-h / ≈ $8.6: still under the $15 stop, though above the
~12 GPU-h target; the R0 prelude's stop rule (projected R1 > 2× the high estimate) halts before that is exceeded.

**Revision 6 (B2) adds no run and no critical-path cost:** the signal panel is CPU scoring of samples `select_ckpt`
already draws, the stable-phase checkpoints cost disk only, the Spearman check uses G-probe's own samples, and the
only new GPU item is the conditional WiSE-FT check above (≈ 0.05 GPU-h).

**Side by side vs smallest first (revision 5; B1 §0.5, D1).** Screening LFM and Qwen3.5 in full in V1 would replace
the 0.2-GPU-h R0 arm with the conditional Qwen3.5 row (1.7–3.9 GPU-h, $1.1–2.6), unconditionally: scenario A
5.4–11.4 GPU-h, $3.4–7.5, and ≈ 20 GPU-h / ≈ $13 at MFU 0.05, past the $10 alarm, to answer a question C1 does not
ask (C1 stops at the first passing rung). The R0 arm costs ≈ 0.2 GPU-h, $0.13 ($0.20 at the ceiling), and buys the
measurements that make a later Qwen3.5 R1 safe to start: s/step, peak memory, kernel route, parity.

**Caps:** alarm **$10**, hard stop **$15** for all Phase 2 Vast spend (G10; owner directive "≤ ~12 GPU-hours /
≤ ~$15 for the critical path"). Scenario A meets both targets with margin even at price ceilings.

**Old vs new budget.**

| Version | Scope | GPU-h | $ | Caps |
|---|---|---|---|---|
| ROADMAP | 3-model screening + A0 vs A0′ on an A100-class card | ~6 | — | — |
| Revision 1 | all rungs, 9 A/Bs, probes, baselines, RIGA, S7 | ≈ 50 | ≈ 47 | alarm 40, stop 60 |
| Revision 2 | critical path / all blocks | 23–42 / 45–68 | 19–34 / 41–59 | alarm 40, stop 60 |
| Revision 3 | critical path (scenario A) | 3.5–7.2 | 2.3–4.7 | alarm 10, stop 15 |
| Revision 4 | critical path (scenario A; ×4 rows counted; MFU 0.05 row added) | 3.7–7.5 (≈ 13 at MFU 0.05) | 2.3–4.9 (≈ 8.5) | alarm 10, stop 15 |
| **Revision 5** | **critical path (scenario A; + R0 Qwen3.5 arm)** | **3.9–7.7 (≈ 13 at MFU 0.05)** | **2.5–5.0 (≈ 8.6)** | **alarm 10, stop 15** |

Where the saving comes from: stop at the first passing rung (no unconditional 2B control: −2.2–4.6 h PRO 6000 and
−2.4–5.3 h for Qwen3.5); 2.5 instead of 4 epoch-equivalents per screening run; no per-component A/B units (−7 h on the
rev-2 path, −10.4 h deferred); C2 reuses R1 as A0′ seed 1 and trains the other three runs only to the selected depth;
precision capacity runs only on a zero-shot failure; baselines on the Mac CPU; the Phase 3 environment smoke moved to
Phase 3; nnU-Net, K = 64 and RIGA removed from Phase 2.

## 5. Gates and KPIs

### 5.0 Statistics and views (fixed before any run)

- **Per-item statistic** = the mean over seeds and over the K samples of that item (K = 8 at T = 1 unless stated);
  greedy decoding is reported separately and never decides.
- **Bootstrap**: 95 % percentile CI, B = 2,000, resampling **items**; real LIDC val items are resampled by **nodule**
  (adjacent slices correlate, ICC 0.15–0.19; B3 P1); synthetic probe items are independent and are resampled within
  their hole-area bin. Paired comparisons
  resample the same items for both arms.
- **Seeds**: multi-seed arms report the CI **pooled** (per-item mean over seeds) and **per seed**. **Caveat (B3
  Addendum A5):** the item bootstrap is conditional on the trained checkpoints; it does not cover seed-to-seed
  variation of the recipe (coverage of the seed-population mean 95 % with no seed effect, 54 % at σ_v = 0.3 logit, 30 %
  at 0.6), and two seeds cannot estimate that variance. Therefore every multi-seed decision also applies **the per-seed clauses of its decision rule** (C2: ≥ 2 of 3 seeds
  ≥ 0.50, none < 0.30; §5.2), with `prereg.py` (M8) authoritative where it and this text differ, and a **seed-discordant** flag is raised when the per-seed CIs lie on opposite
  sides of the decision threshold. C2 has 3 A0′ seeds by default (§5.2); a fourth seed in deterministic mode (cost
  measured in the R0 prelude) is added only if the 3-seed rule is decided by a single seed's side of 0.50 and the
  owner approves. That seed is part of the C2 verdict: C2 is re-judged on 4 seeds with the per-seed clause restated as
  ≥ 3 of 4 seeds ≥ 0.50 and none < 0.30 (in `prereg.py`), and only then is the outcome (pass / partial / fail, §5.2)
  assigned, so it precedes R-x8 / R-ctrl and two conditionals cannot fire on one result. Single-seed runs (screening, the §6.3 diagnostic) carry no seed claim.
- **Multiplicity**: the only other decisions are pre-registered rules without multiplicity correction and with no
  claim attached: the C1 stop rule (gates), D12 checkpoint selection, and the single §6.3 fallback diagnostic. **C2 is the single confirmatory test** of Phase 2: its primary metric is A0′'s pooled hole
  recall, mean over K = 8 at T = 1 and over the 3 A0′ seeds, with a stratified (within-bin) item bootstrap (D20), judged by the
  §5.2 3-seed rule (95 % CI lower bound > 0.30 plus per-seed clauses). C1 is a gate decision, not a hypothesis test.
- **Views.** *Centred val* = `targets.frame` (P1). *Shifted val* = the same S×S window moved by one fixed-seed uniform
  offset per item, union inside with ≥ 2 px clearance, no other outlined nodule inside, targets and caches rebuilt (A4
  §2.5 proposal 1). *Rescaled val* = one fixed-seed draw per item moving S to a neighbouring divisor (LIDC 50 ↔ 100 ↔
  125; RIGA 200 ↔ 250 ↔ 500) plus a shift, where the invariants allow (union inside with ≥ 2 px clearance, S ≤ 500, no
  other nodule); items with no admissible neighbour keep S and are counted and reported (if > 30 %, the rescaled check
  is noted as weaker than stated here). It removes the size-from-S cue that the
  translation-only view leaves (training already jitters S, D10). Gate sampling is unconstrained (D13).

### 5.1 Phase 2 entry gates (revised; per rung; `phase2_entry` per dataset, B4-5)

| Gate | Definition | Threshold | Evidence for the threshold | Measured how / when |
|---|---|---|---|---|
| **G-format** | parse and legal rates over all K = 8 samples at T = 1, shifted val. **Reported beside it (not gated):** truncation share, mean / p95 completion tokens, held-out NLL | parse ≥ 98 %, legal ≥ 90 % | unchanged from RESULTS §8; legality is the binding precision constraint on LIDC (σ ≥ 0.6 px vertex error drops legal below 90 %; A4 S4); truncation is the runaway-cost failure (A5 §0.6: 1,536 tokens, 5–15× a normal answer) | `select_ckpt` (vLLM after parity, else HF) → `eval_run` → `summarize`; every WSD branch (D12) and the selected checkpoint |
| **G-geometry** | median IoU_tc per structure (illegal = 0) and the share of samples legal with IoU_tc ≥ 0.5 on every object, on the **shifted and the rescaled** val | LIDC ≥ 0.60; disc ≥ 0.85; cup ≥ 0.75; share ≥ 70 %; on both views | RESULTS §8 values; human LOO 0.768 / 0.945 / 0.815 (A4 S2); image-blind prior on shifted frames 0.00 / 0.32 / 0.07 (A4 S6) | same run; selected checkpoint |
| **G-prior** | per-item IoU_tc = the K = 8 item mean with illegal = 0 (§5.0). (a) centred val: share of items where it beats the image-blind prior's, sign test, **report-only for C1** (revision 7); (b) shifted and rescaled val: paired mean (model − prior) IoU_tc CI lower bound > 0, **gate** | (b) gates; (a) reported against 60 % and beside MedSAM-oracle's and SAM 2.1-tiny's centred shares (D17). If the owner wants (a) as a gate, its threshold is set at ≤ SAM 2.1-tiny's share − 5 pts and written into `prereg.py` (M8) **before V1, never after seeing C1** | A4 §2.5 proposal 2 (prior on centred frames passes RIGA, disc 0.904, cup 0.783, and narrowly fails LIDC, 0.588 / 65.3 % vs 0.60 / 70 %; A4 S5); B1 T2 "uses the image" rule: after centred-crop SFT every B1 arm was significantly *below* the prior on shifted crops (paired −0.026 to −0.045, B1 §0.1), so (b) is the check that catches the centring shortcut. (a) is report-only (JUDGE_R3 MUST-FIX 1): the centred prior scores 0.588 with 65 % of items ≥ 0.5, so a model at legal 90 % and shifted median 0.60–0.65 can win < 60 % of centred items while passing (b), G-geometry on both views and G-vision, and the 60 % (A4 §2.5 proposal 2) has no reference model's share behind it; image use is decided by (b), G-geometry and G-vision | prior scored by `reward.score_mask` on all three views (A4 `null_prior.py`); MedSAM-oracle and SAM 2.1-tiny centred shares in the Mac baseline run before V1 (D17); selected checkpoint |
| **G-precision** | precision probe (D19), zero-shot, then the capacity run only on a zero-shot failure, 6–16 px objects, patient-disjoint probe-test patches | legal ≥ 95 %; median ASSD ≤ 0.6 px; median IoU ≥ 0.80; median centroid error ≤ 0.5 px; paired per-item IoU win over the centroid disc ≥ 60 % | A4 §6.2 table (σ = 0.5 → 95.7 % legal; σ ≈ 1 → ASSD 0.61; 1 px shift −23 IoU pts); B3 P2 discrimination (disc fails only via the paired rule: ≤ 49 % wins vs σ ≤ 0.5 at 78–100 %). Legal ≥ 95 % is knife-edge (σ = 0.5: 92.3 % gen. 1, 96.4 % gen. 2) | R1: **zero-shot** on the selected branch (inference only, 1,000 probe-test items, greedy + K = 8). Pass → done. Fail → the 20k-item capacity run (~0.3 GPU-h): capacity pass → rung kept, LIDC SFT reported as the bottleneck (A4 §6.2 reading); capacity fail → the §6.3 rule (R2 component) or the rung is dropped from LIDC |
| **G-topology** | topology probe (D20, B3 Addendum A), selected A0′ checkpoint: (a) simple controls with correct (β0, β1) = (1, 0) — **gate**; (b) **hole recall on one-hole items** (pooled; per hole-area bin) — **report-only for C1** (revision 4): a hole-recall miss never fails a rung or climbs the ladder; it is decided by C2 only, and R1's value is reported early so R-ctrl's trigger can be evaluated | (a) ≥ 95 % (aligned with C2 clause (4), false non-simple ≤ 5 %); (b) reported against 0.30 | replaces "topology-any ≥ 90 % on unanimous items" (val has 0 unanimous non-simple items: vacuous; B3 P1). The σ = 0.5 synthetic model scores 0 false non-simple (B3 Addendum A4), so 95 % is not tight. 0.30 is the C2 fail line (A1 H-A1-6; B3 Addendum A6); per-bin estimates are reliable at n = 120: P(estimate < 0.30) = 0.010 at true recall 0.40, 0.998 at 0.20 (B3 A3) | R1 (selected branch) and R3–R5, K = 8 at T = 1 (greedy reported) |
| **G-signal** | at T = 1, G = 8, shifted val: duplicate share (hash of format/legality flags + packed rasters, B2-6) and share of groups with geometric σ above its floor | duplicates ≤ 10 %; σ-share ≥ 50 % | RESULTS §8; DAPO (equal-reward groups give zero advantage); K12 smoke 0 % / 100 % on a non-collapsed checkpoint | every WSD branch (D12 eligibility; a floor, not a ranking key: σ-share and duplicates were saturated in B2 §2.1) |
| **G-vision** (revision 5) | "the model still looks at the image", B1 §2.5 protocol promoted to `eval/vision.py` (M4): (a) mean pairwise cosine distance of the projector output (`get_image_features`) over a fixed image set (4 centred LIDC val crops + a black image, as `lab/B1/vision_feats.py`), as a ratio to the base model's value measured in R0; (b) teacher-forced next-token TV over the first 24 target tokens of 8 fixed LIDC val items, real crop vs a black image and real crop vs another item's crop (`lab/B1/image_use.py`) | (a) ratio ≥ 0.5; (b) both TVs ≥ 0.002 | B1 §2.5: collapsed checkpoints (Qwen3.5 and SmolVLM, all-linear LoRA, 36 / 150 steps) have ratio 0.001–0.064 and TV ≤ 0.0007; intact ones (LFM all-linear 36 / 150 steps; LM-only Qwen3.5 and SmolVLM) ratio 0.68–1.0 (frozen = 1) and TV ≥ 0.004. (a) sits 1.4× below the lowest intact ratio and 8× above the highest collapsed one; (b) 2× below the lowest intact TV and ~3× above the highest collapsed one. B1's diagnostics are post-hoc and small (4–8 items), and TV is diluted by syntax tokens, so (b) holds only under this exact protocol | every WSD branch in `select_ckpt`, on the CPU as in B1, so ② stays GPU-free (minutes per branch; D12 eligibility); (a) and (b) also in `train_sft` at each rolling checkpoint: both below threshold → abort and D7 repair; one alone → warn-only, logged |
| **G-probe** | "GRPO can improve this checkpoint": 20 GRPO steps (`grpo_config` with a **pre-registered probe LR**: **5e-6 full FT, 5e-5 for LoRA**, constant; G = 8, 32 LIDC train prompts per step with ×4 disagreement oversampling, HF generation, `use_vllm=False`); fixed G-probe set = 128 shifted-val LIDC items, K = 8 at T = 1. **Noise floor at zero training cost:** at step 0 the set is sampled **twice** with different sampler seeds; the per-item difference gives the null CI (half-width h₀). Score per item: the **raw composite** R = 0.3·IoU_tc + 0.3·F5 + 0.1·legal + 0.3·g_area (v4 weights on the raw, legality-gated terms; TRL's aggregator is group-mean-centred, so its mean is ≈ 0) | **PASS**: paired per-item gain (R₂₀ − mean of the two R₀ draws) CI lower bound > 0, legal at T = 1 falls by ≤ 5 pts, σ-share at step 20 ≥ 50 %. **FAIL** (counts): the gain CI lies entirely below the null CI (R degraded beyond sampling noise), or legal falls > 5 pts, or σ-share < 50 %. Otherwise **INCONCLUSIVE**. Reported too: h₀, and the training-log form (mean logged R, steps 11–20 vs 1–5) | G-signal is necessary, not sufficient (2606.18487: advantage variance p(1−p)(g−1)/g can vanish). The probe LR is **the probe's, not Phase 3's**: production `grpo_config` uses 1e-6 (full FT) / 1e-5 (LoRA), which over 20 steps may not move R beyond h₀; SimpleSeg's RL used a constant 2e-6 for ~3k steps at 7B (Tables 4–5); 5× the full-FT production LR stays 10× below the SFT peak (5e-5) and LoRA keeps its 10× ratio (LoRA Without Regret). Thresholds, null procedure and LR pre-registered in `prereg.py` (M8) | ④ in §3.2, on the selected branch. Not PASS → the D12(c) repair (next-earlier eligible branch), probed once. **Still not PASS → no ladder climb**: the rung passes C1 with a "G-probe warning" handed to Phase 3 (which starts with its own LR choice and a longer run), because D12 eligibility already requires G-signal; the ladder climbs only when no branch passes G-signal (a §5.1 gate failure) |
| **G-x86** | K3: 256-rollout batch wall (best of serial / pool), RIGA-only (synthetic) and mixed; K4: p99 per rollout; hashes serial = pool | K3 ≤ 250 ms; K4 ≤ 5 ms; equal hashes | RESULTS §8; B4 (Mac pool-8 51–61 ms; A6 serial RIGA 241 ms at load 9) | `k3_harness` in R0 (idle) and once under R1's load (§4.7 step 6) |
| **G-cost** | measured s/step, SFT and GRPO, per (model, card); projected Phase 3 Track A GRPO cost | measurement recorded; projection ≤ $150 on-demand, else coordinator re-plans (card, arms, seeds) | RESULTS §8 cost gate; A3 §3.4 (~$105 for ~70 GPU-h at $1.5/h) | R0 timing prelude and R1 logs (SFT); R2 (GRPO, HF generation); colocate step time measured at Phase 3 start |
| **G-parity** | vLLM sampler vs HF: prompt ids; per-token \|Δ logprob\|; distribution parity (32 items × K = 8: parse, legal, median IoU_tc, mean length, duplicates) | 100 % ids; median ≤ 0.02 nats bf16, p99 per sequence ≤ 0.5; every metric difference inside its bootstrap 95 % CI | A5 H-A5-1/2; B5 Q2 (fp32 1e-7–4e-6) | R0 / start of each rung's R1 (ids, logprobs on base weights, few items); distribution parity on R1's selected branch. Fail → HF sampler is the gate path for that model (D13) |
| Housekeeping | LFM licence recorded; novelty sentences narrowed (P1 §8); RIGA permission status recorded | done | RESULTS §8 | coordinator |

"Near MedSAM" (A4 §3.2) is **reported, not gated**, with every row labelled by view: model vs MedSAM (oracle box from
that view's union) − 0.10 per structure (centred: LIDC 0.68, disc 0.84, cup 0.74), on the centred **and the shifted**
val, next to the jittered-box MedSAM row, SAM 2.1-tiny and the prior on the same views (D17).

### 5.2 Claim tests

| Claim | Pre-registered procedure | Pass | Evidence for thresholds |
|---|---|---|---|
| **C1 smallest viable** | rungs screened smallest first with the one recipe bundle (§6.0); select per D12; §5.1 on LIDC (and RIGA if permitted); **stop at the first rung that passes** | C1 = the first (smallest) rung passing every §5.1 **gate** (G-prior (a) is report-only; G-topology (b) is report-only, so a hole-recall miss goes to C2's "holes left to GRPO", not up the ladder; G-probe never climbs the ladder on its own, §5.1). No larger rung is run once one passes: passing is sufficient for "smallest viable", and a larger model's score adds nothing to it. Qwen3-VL-2B runs only as the next rung or as R-ctrl (§6.1) | P1-D12; ROADMAP Phase 2; owner directive (revision 3) |
| **C2 A0 vs A0′ (holes; confirmatory)** | C1 rung, bundle; depth = that of the checkpoint that received the final G-probe verdict (the D12(c) repair branch if it was probed, not necessarily the first D12 pick); A0′ = bundle with ×4 non-simple rows, A0 = the same rows with holes emptied; A0′ seed 1 = that checkpoint, **A0′ seeds 2–3 and A0 seed 1** trained to that depth (R3–R5): **3 A0′ seeds, 1 A0 seed** (revision 4: the A0 arm scores 0 hole recall by construction and enters no hole decision, while the A0′ seeds are the scarce resource, B3 A5); topology probe D20 (360 one-hole items, 120 per bin); primary metric and bootstrap as D20 / §5.0. A0 is run under the same protocol; on hole items A0′ − A0 equals A0′'s recall (A0 = 0 by construction), so it is reported for completeness only | **Pass iff all hold** (B3 Addendum A6 adapted to 3 seeds, revision 4): (1) A0′ pooled hole recall (mean over K = 8 and the **3 seeds**, within-bin bootstrap) ≥ 0.50 with 95 % CI lower bound > 0.30; (2) **≥ 2 of the 3 seeds' pooled point estimates ≥ 0.50 and none < 0.30**; (3) no bin's point estimate < 0.30; (4) false non-simple on controls ≤ 5 %; (5) legal ≥ 90 %. Reconciliation with B3 A5's caveat (item CI conditional on the trained checkpoints; 54 % coverage of the seed mean at σ_v = 0.3): the item-bootstrap CI is kept for item noise, and seed variation is handled by the per-seed point-estimate clause; with 3 seeds a majority rule tolerates one weak seed while the < 0.30 floor keeps the fail line; the across-seed t-interval of the 3 seed means (2 df) is reported, not gated. Seed-discordant flag (seeds on both sides of 0.50) is reported. **A0 guard (pre-registered, report-only):** paired (A0′ 3-seed mean − A0) shifted-val median IoU_tc CI lower bound > −0.02, else "hole-awareness costs geometry" is reported. **Outcomes (revision 7), each with at most one conditional run, owner approval required:** **pass** (all of (1)–(5)) → none; **partial** (not pass, pooled A0′ recall in [0.30, 0.50), or ≥ 0.50 with another clause failed) → one ×8 A0′ seed (E1, recipe question), not run when clause (4) false non-simple failed (more hole rows would raise it); **fail** (pooled < 0.30) → R-ctrl (capacity question), "holes left to GRPO" unless it passes. R-ctrl never follows a partial and ×8 never follows a fail. Part recall is descriptive.  | B3 Addendum A3 (±3.1 pts pooled; per-bin floor reliability), A4 (σ = 0.5 px synthetic model: pooled 0.76 [0.73, 0.78], bins 0.65 / 0.81 / 0.90), A5 (seed caveat), A6 (rule); A1 H-A1-6 (≥ 50 % / < 30 %); A2 H-A2-5 (hallucination); JUDGE_R1 MUST-FIX 1, 3 |

### 5.3 Engineering and operations KPIs

| KPI | Definition | Threshold | Evidence | Measured |
|---|---|---|---|---|
| Coverage | `make cov` line coverage of `src/objdel` | ≥ 90 % | K13 (93.03 % at Phase 1 close) | every merge |
| Tests / lint | `make test`, `make lint` | 0 failures, clean | CLAUDE.md | every merge |
| CI | GitHub Actions on `ubuntu-latest` | green | first run green: run 39359008526 (I34 CI half closed) | every push |
| Hash reproducibility | `xarch_hashes`: Vast vs Mac 22/26 + 4 accepted; Vast vs Docker x86 26/26 strict; per-view `target_set_hash` identical across two builds and process counts; identity view = `1af819f4…` subsets | exact | M9; B3 P3.3–P3.4 (`e86e07e6…`) | R0 (V1); every view build |
| Phase 3 env | fast suite + `xarch_hashes --strict` in the re-pinned image | 0 non-flaky failures; 26/26 | B5 Q1 | every re-pin; Phase 3 start |
| Data integrity | bundle SHA verified on box | 100 % | A3 §6.2 | every session |
| Governance | RIGA bundle never sent to a non-datacenter offer; ledger row per session | 0 violations | D16; H-GOV unit test | every session |
| Spend / orphans | ledger total; instances alive outside `active_runs.json` | ≤ $10 alarm, ≤ $15 stop; 0 orphans > 1 h | G6, G10 | hourly sweep, per session |
| GPU utilisation | `nvidia-smi dmon` mean during SFT | ≥ 70 % (else CPU-side fix first) | A2 §6.4; H-A2-10 | R0 prelude, every run log |
| Audit recalibration (I35) | hack-audit thresholds re-fit on the selected checkpoint's K = 64 samples | documented in §9 | I35, I17/I18 | moved (§6.4); required before the first Phase 3 GRPO run |

## 6. The recipe bundle and the only Vast runs (revision 3, owner cost directive)

### 6.0 One recipe bundle, used in every run

Chosen from the evidence; no component gets a standalone A/B. A0 runs use the identical bundle with holes emptied.

| Component | Setting | Evidence that justifies adopting it untested |
|---|---|---|
| Crop shift + scale jitter (D10) | shift uniform, union inside, ≥ 2 px clearance, no other nodule; margin [8 px, 60 %] so S moves between divisors | mandatory: the centred frame leaks location and size: the image-blind prior passes RIGA's centred gates and collapses under shift (0.00 / 0.32 / 0.07; A4 S5–S6), and every B1 arm trained on centred crops drew at the frame centre (94–98 % of shifted samples closer to the frame centre than to the nodule, B1 §2.5); the gates are defined on shifted / rescaled views |
| Offline dihedral targets (D9, D10) | D4 for LIDC (RIGA: hflip only); one op per rater row; crop-before-trace; per-view hash | 215,208 variant targets 100 % certified, LIDC pixel-exact, byte-identical rebuilds (B3 P3–P4); exact transforms are standard for LIDC (A4 §4.1); zero runtime cost |
| CT window jitter (D10) | L −300 ± 60, W 1400 ± 80 HU on **int16 HU context crops** (one npz per LIDC item, written by M1 from `data/store`), noise σ ≤ 10 HU, windowed at collate time | random windowing for CT (liver tumours, 2510.08116; A4 §4.1): one study, the weakest evidence in the bundle, kept because it is applied before windowing (certification untouched), costs ≈ 0 and is the §6.3 rule-3 removal; at zero jitter the HU path reproduces the prepared PNG byte-for-byte (M1 test) |
| All voting raters (D9) | one row per voting rater, every epoch | same expected gradient as random-rater sampling, keeps the rater distribution D5 matches (A2 §4.6; P1 I27) |
| Hole oversampling (D11) | non-simple LIDC train targets ×4 (414 of 23,459 targets → ≈ 6.7 % of rows) | A2 §7.4 and A1 §5.2 recommend ×4–8; P1 used ×4 at 3.9 % prevalence (RESULTS §7 #9); the low end limits hole hallucination (A2 H-A2-5 guard, checked by G-topology (a)) |
| Vision tower + projector (D7) | LoRA: frozen; full FT: 0.1× LR; G-vision on every branch, in-training projector check | all-linear LoRA incl. vision collapsed Qwen3.5's and SmolVLM's image features in 36 steps, LM-only LoRA kept them at equal legality / IoU (B1 §2.4–2.5); A2 §4.2 (unfreeze for localisation; 0.1× its own proposal) for full FT |
| Fine-tuning (D6) | full FT ≤ 0.8B (Qwen3.5 by R0's peak memory); LoRA r16 on the language model for Qwen3-VL-2B | A2 §4.1; B1 §2.4 |
| Per-model settings | Qwen3.5: `qwen3_5_nothink_training` template, bf16; `MAX_COMPLETION` 1,024 (LFM, Qwen3.5); every sampling field explicit (T 1, top_p 1, top_k 0, min_p 0, repetition 1, presence 0), overriding both cards (D13) | B1 §1.1, §2.6, §4.5 |
| Schedule (D5) | WSD, 3 % warmup, stable 2 epochs, 15 % decay branches at 0.5 / 1 / 2 epochs; LR 5e-5 (LoRA 5e-4); effective batch 32 | A2 §4.3; Quagmires (1–2 epochs) |
| Throughput (D4) | `group_by_length` + `length` column; 6 persistent workers; grad-ckpt off when peak ≤ 20 GB; SDPA; bf16 | padding 7–32 % → ≤ 0.2 % (A2 §1.2); collation 10–25 ms/row needs workers (A6 §1a, §2); measured in R0, not A/B-tested |
| Loss (D8) | assistant-only NLL; **no number-token loss** | NTL's evidence is text-only, untested on VLM coordinates, weakest on LFM's 3-digit tokens; it stays a triggered option (D8, §6.4) |
| Not in the bundle | label smoothing, DFT, constrained decoding, coordinate tokens, mixing | D8, D13, D2, D11 |

### 6.1 The Vast runs (complete list)

Ladder (D1): **LFM2.5-VL-450M** → Qwen3.5-0.8B → Qwen3-VL-2B (PRO 6000) → stop and
re-plan (§8.2). All on LIDC; all in one V1 session where possible (§4.7).

| Run | Purpose | Design | Decides / measures | Est. (LFM, 4090) |
|---|---|---|---|---|
| **R0** | box readiness | image check, `xarch_hashes`, `k3_harness` idle, parity on a few items (LFM and Qwen3.5), 4 × 30-step timing prelude and the 30-step Qwen3.5 arm (§6.2, D1) | G-x86, G-parity (ids, logprobs), throughput and determinism cost, projected R1 time for LFM and Qwen3.5 (stop rule), Qwen3.5 kernel route and D6 choice, G-vision base references | 0.9 h |
| **R1(k)** | **C1** screening of rung k | bundle, A0′, seed 1, WSD 2 epochs + branches 0.5 / 1 / 2; `select_ckpt` after `train_sft` exits: shifted val K = 8 per branch; the D12-selected branch also gets K = 16, the rescaled val, the precision probe (zero-shot), the topology probe and the **free attribution check** (identity vs one rotated shifted-val view, inference only, ≈ 0.05 h, report-only: rotated ≈ identity means dihedral training transferred) | §5.1 gates; checkpoint depth (D12 rule R-B, B2 signal panel); SFT s/step. If G-vision (a) and (b) both fail in training, the D7 repair reruns it frozen (§4.9 conditional) | 0.76–1.70 h |
| **R2(k)** | "GRPO can improve it" | G-probe on the selected checkpoint (§5.1: probe LR 5e-6 / 5e-5 LoRA, step-0 double sample as the null); the D12(c) repair branch only if not PASS; never climbs the ladder alone | G-probe; GRPO step time and split (G-cost); per-item Spearman of step-0 sd_geom / r_std vs gain (D12 veto check) | 0.2–0.5 h |
| **Stop** | — | the first rung whose R1 + R2 pass is C1; no larger rung is run | C1 | — |
| **R3–R5** | **C2** (confirmatory) on the C1 rung | **A0′ seeds 2 and 3; A0 seed 1** (same ×4 rows, holes emptied), each trained to the depth of the checkpoint that received the final G-probe verdict (the D12(c) repair branch if it was probed; stable to that depth + one 15 % decay); A0′ seed 1 = that checkpoint. Same run count and budget as revision 3 (A0 seed 2 became A0′ seed 3) | topology probe → §5.2 3-seed rule and outcome (pass / partial / fail); A0 geometry guard | 1.6–4.2 h |
| *R-diag* | bundle-level fallback | at most **one** in Phase 2: §6.3 | whether a single bundle component caused a gate failure | 0.5–1.4 h |
| *R-x8* | conditional recipe check (revision 7) | **only if** C2's outcome is **partial** (§5.2) and clause (4) passed: one A0′ seed with non-simple rows ×8 instead of ×4, at the C2 depth, owner approval. One seed, so its pass is **exploratory**, as R-ctrl's: clauses (1) (on its seed), (3), (4), (5), judged on a topology-probe draw **disjoint** from C2's (same generator and bins, fresh seed; §9.1 revision 2, E1); pass → "more hole rows may lift recall (exploratory, 1 seed)", else "holes left to GRPO" | whether more hole rows lift recall into the pass range | 0.5–1.4 h |
| *R-ctrl* | conditional capacity control | **only if** C2's outcome is **fail** (pooled hole recall < 0.30, §5.2; a partial goes to the ×8 seed instead): one Qwen3-VL-2B A0′ seed at the same depth on the PRO 6000 (owner approval, as for every C2 conditional, §5.2). One seed cannot meet the 3-seed rule, so R-ctrl's pass is **exploratory**: clauses (1) (on its seed), (3), (4), (5); if it passes, report "holes may need ≥ 2B (exploratory, 1 seed)"; else "holes left to GRPO". Also run Qwen3-VL-2B when every smaller rung fails (then it is simply the next rung). Not run otherwise | capacity vs recipe for holes | 0.8–2.3 h (PRO 6000) |

The precision-probe capacity run (~0.3 GPU-h) is run only when the zero-shot G-precision fails (§5.1).

### 6.2 Measurements folded into those runs (no separate experiments)

| Measurement | Folded into | How | Replaces (or source; cost) |
|---|---|---|---|
| Reward K3/K4, x86 hashes | R0, once under R1 load | `k3_harness` (5 batches), `xarch_hashes` compare | first-session steps of rev 2 |
| Loader / throughput settings | R0 timing prelude | 2 × 30 steps from R1's init, workers 0 vs 6 (data order identical); grouping benefit from the logged padding share; peak memory decides grad-ckpt off | E6 (H-A2-1, H-A2-10, A6 H-5) |
| Determinism cost | R0 timing prelude | 2 × 30 steps, deterministic algorithms on vs off (`CUBLAS_WORKSPACE_CONFIG` set for the session); overhead ≤ 15 % → deterministic mode for any fourth C2 seed | E7 cost (bitwise repeat: §6.4) |
| Qwen3.5 cost, memory, kernels | R0 Qwen3.5 arm (D1) | 30 discarded steps; GDN route per function (D15); peak memory → D6 | Qwen3.5's UNVERIFIED 2× cost (§4.9); B1 §0.5's side-by-side screen |
| G-vision base references | R0 | `eval/vision.py` on base LFM and Qwen3.5 | — (needed by G-vision's ratio) |
| Sampler parity | R0 and R1 | prompt ids on 32 items, logprobs on 8 items × n = 8 (base weights, LFM and Qwen3.5 in R0); vLLM vs HF on 32 items × K = 8 on the selected branch; the same comparison gives the HF vs vLLM wall time | H-A5-1, H-A5-2, H-A5-8 |
| SFT step time, GPU utilisation | R1 logs | s/step, data wait, `nvidia-smi dmon` | H-STEP |
| GRPO step time and split | R2 | generation / scoring / update | RESULTS §8 cost gate |
| B2 signal panel per branch (r_std, headroom, mixed, legal-but-wrong, mean \|A_v4\|, sd_geom / legality split, σ-share, duplicates) | R1 `select_ckpt` | CPU scoring (M11) of the shifted-val K = 8 samples it already draws per branch; feeds D12 R-B | B2 §5; 0 GPU |
| Depth curve, 3 points | R1 | the 0.5 / 1 / 2-epoch decay branches plus the stable-phase checkpoints at those points, kept (model-only) | B2 §5; 0 compute, disk only |
| Does the signal predict GRPO gain (H-A2-3 mechanism; D12 veto check) | R2 | per item (128): Spearman of step-0 sd_geom, step-0 r_std and the mixed flag (from the step-0 double sample) against the gain R₂₀ − R₀, bootstrap CIs; retains sd_geom or flags disagreement (confounded by legality headroom); never switches the statistic (D12) | B2 §5–6; 0 (G-probe's own samples) |
| Shallow vs deep under GRPO (D12(a) / (b)) | R2, only when the repair branch is probed anyway | both branches' paired gains side by side, report-only | B2 §5; 0 extra |
| Entropy proxy | R1 logs | TRL `entropy` and held-out NLL per branch (already logged) | B2 §5; 0 |
| WiSE-FT 0.5-epoch ↔ selected | R1 / R2, only when the D12(c) repair fires | merge at α = 0.5, shifted val K = 8, report-only | H-A2-7 proper; ≈ 0.05 GPU-h |
| Watchdog | V0 | 10-min dry run on a 3090 | H-WATCH |

### 6.3 Bundle-level fallback (pre-registered; at most one diagnostic run in Phase 2)

**Trigger:** the selected checkpoint of a rung fails G-format, G-geometry, G-prior (b), G-precision (after its capacity
run) or G-topology (a). G-signal and G-probe failures go to the D12(c) repair, and G-vision failures to the D7 repair (tower and projector frozen, once per rung), not here. Only the first such failure in
Phase 2 gets a diagnostic; later rungs use the bundle as it then stands.

**Component to remove (first matching rule):**
1. ≥ 50 % of the illegal samples carry a hole-related defect (hole outside its outer, hole crossing a ring, wrong hole
   winding), or G-topology (a) fails → **hole oversampling ×4 → ×1**.
2. Else, G-geometry, G-prior (b) or G-precision fails with legal ≥ 90 % → **dihedral variants off** (identity view only):
   the largest change to the training distribution that the val never shows, spending a ≤ 0.8B model's capacity on
   8 orientations.
3. Else (parse, legality or truncation failure) → **CT window jitter off**.

Never removed: crop shift / scale jitter (the gates are defined on shifted and rescaled views; removing it reintroduces
the leak, which every B1 arm learned, B1 §2.5), all voting raters, the WSD branches.

**Run:** R1's design trained to the failing run's selected depth only (one decay; **the deepest branch when no branch was eligible**), same seed, with that component
removed; evaluated on the failed gates plus G-probe. **Pass** → the rung passes with bundle-minus-component, which becomes
the bundle for C2 and any later run (recorded in §9). **Fail** → the bundle is restored and the ladder climbs; no further
diagnostics, no factorial.

### 6.4 Removed or triggered (owner directive, revision 3)

| Item (rev 1–2 ID) | Status | Trigger that brings it back |
|---|---|---|
| E1 hole oversampling ×1 / ×4 / ×8 | removed; ×4 fixed in the bundle | §6.3 rule 1 (drop to ×1); R-x8: one ×8 A0′ seed only on C2's **partial** outcome (§5.2, §6.1) and the owner's approval |
| E2 number-token loss (M6 code) | removed; not in the bundle (D8) | the passing rung's legality is within 3 pts of 90 %, or G-precision fails on two rungs; or Phase 3 planning |
| E3 crop jitter on / off | removed for good (revision 5: B1 showed the leak is learned without jitter, B1 §2.5) | none: jitter is mandatory and the shifted / rescaled gates plus G-prior detect a leak |
| E4 vision tower frozen / unfrozen | removed (LoRA frozen, full FT 0.1× in the bundle; G-vision and the D7 repair cover collapse) | G-precision fails on two consecutive rungs |
| E5 RIGA T = 2 mixing | removed | RIGA permission and a RIGA gate miss at natural proportions |
| E6 loader / grad-ckpt / grouping A/B | folded into R0 (§6.2) | — |
| E7 determinism | cost folded into R0; bitwise repeat removed | seed-discordant flag on C2 |
| E8 coordinate tokens (H-A1-4) | removed | R2 shows generation > 60 % of the GRPO step → Phase 3 planning |
| E9 HU window jitter A/B | removed (jitter in the bundle) | Phase 9 ablation only |
| E10 rejection-sampling SFT (D28) | removed | legal in [85 %, 90 %) at every branch of the last affordable rung, with owner approval |
| Precision-probe capacity runs | conditional | zero-shot G-precision fails |
| Unconditional Qwen3-VL-2B control; screening of larger rungs | removed | a smaller rung fails, or R-ctrl's condition (§6.1) |
| nnU-Net v2 training (D17 row 4) | removed from Phase 2 | Phase 9 paper baselines |
| K = 64 evaluation, I35 audit recalibration | moved | before the first Phase 3 GRPO run |
| HD95 / ASSD / centroid in `objdel.eval` | removed | paper tables (the probes keep their own) |
| H-A5-5 constrained diagnostic, H-A5-4 caching study | removed | Phase 3 planning, or G-format parse < 98 % |
| S7 Phase 3 environment smoke (colocate) | moved | first step of Phase 3 |
| H-PULL resume test | removed | first host failure (resume procedure documented, A3 §6.2) |
| MedSAM / SAM 2.1 GPU runs | moved to the Mac CPU (inference only) | — |
| S6 RIGA runs | removed from the critical path | written RIGA permission (D16) and owner approval of a tranche |

## 7. Milestones and task graph (test-driven: tests first, then promote lab code, then the gate)

| M | Milestone | Owner | Promotes from | Tests (written first) | Depends |
|---|---|---|---|---|---|
| M1 | **Frames and views**: `frame(item, cfg, offset, side, window)`; crop-before-trace in `_rings`; `targets.augment` (D4 on masks + images at frame resolution, row-spec hash); `build_views` (per-epoch row specs, certified targets, per-view hash, token filter; **LIDC context crops as int16 HU npz written from `data/store`**, not 8-bit PNGs; RIGA context crops as RGB PNG); `build_eval_views` (shifted and rescaled val; targets + caches) | targets | `lab/B3/dihedral.py`, `crop_trace.py`; A4 `null_prior_jitter.py`; B1 `common.py` views | offset (0,0) = current frame and PNG hash; invariants (S divides 1000, S ≤ 500, union inside, ≥ 2 px clearance, no other nodule); identity view hash = `1af819f4…` subsets; 100 % certified per op; LIDC variants pixel-exact; two builds byte-identical; collate-time image = offline reference on sampled rows; **byte identity: HU crop + `window=(−1000, 400)` at zero shift / scale / jitter reproduces the prepared LIDC PNG byte-for-byte (and its `image_sha256`)**; HU crops equal the store's int16 pixels on the crop window | — |
| M2 | **Training data and configs**: `sft_rows` over views (op per rater row, ×4 non-simple (bundle), A0 transform, mixing T, `length` column, voting indices stored, I45 manifest-driven `cache_loader`); `configs`: `qwen3.5-0.8b` (MAX_COMPLETION 1,024, `qwen3_5_nothink_training`, bf16 base), explicit per-model sampling overrides (D13), **`lora()` changed to LM-only** (`exclude_modules` = vision tower + projector; today it is all-linear incl. vision, the setting that collapsed vision in B1 §2.5), vision + projector 0.1× param groups for full FT, WSD, grouping, workers, grad-ckpt per model | train | A2 §7; A6 §4; `lab/B1/common.py` (`lora("lm")`, `VISION`) | row counts and repeats; A0 rows have `"holes":[]` and keep parts; op hash deterministic; length column = tokens; config snapshot tests per model; I27 voting-only kept; `lora()` targets no module matching the vision / projector regex on each model's module names (recorded name lists, no download); full-FT param groups put every vision / projector parameter at 0.1× LR; sampling kwargs per model equal T 1, top_p 1, top_k 0, min_p 0, repetition 1, presence 0 | M1 |
| M3 | **Reward serving**: `train/serve.py` (per-dataset preload, `CacheRef`, `fork_pool` warm, `gc.freeze`, initializer), `cache_loader(maxsize=None)`; used only by scoring processes, never by `train_sft` (D14) | train | `lab/B4/serve.py`, `preload.py`, `trl_pattern.py` | hashes pool = serial (all workloads); one dispatch per generation batch; killed worker → `BrokenProcessPool` within 1 s; pool created before CUDA (assert `torch.cuda.is_initialized()` is False); preload digest checks; preload of a LIDC-only run loads no RIGA cache | — |
| M4 | **Eval extensions and probes**: HM-IoU16, NSD@τ (HD95/ASSD/centroid in `objdel.eval` deferred, §6.0); truncation/token stats in `summarize`; `eval_run` passes prepared caches; `eval/probes.py`: precision probe (generator 2, E = 0.2, A = 0.15; probe-train and probe-test patches from patient-disjoint LIDC train halves) and topology probe on val patches per D20 (B3 Addendum A: 120 one-hole items per bin, 160 two-part, 160 simple; holes and parts jointly resampled from real train records, size-rule minimum; stratified bootstrap); `eval/vision.py` (G-vision: projector cosine ratio, image-use TV, fixed image / item sets); revised `phase2_entry` (§5.1: shifted + rescaled views, prior, probes, G-vision, G-probe) **replacing** the unanimous-topology gate | eval | `lab/B3/probe_lib.py`, `precision_probe.py`, `topo_probe.py`, `probe_rule.py`; A4 `probe_proto.py`; `lab/B1/vision_feats.py`, `image_use.py` | H-A4-14: HM-IoU = brute-force assignment (k ≤ 4, K ≤ 8, Hypothesis), NSD = 1 and HD95 = 0 on raters' own targets, NSD = DeepMind reference to 1e-6; probe exactness (ASSD = HD95 = 0 on ground truth); calibration (disc IoU 0.806 ± 0.03); discrimination (exact/σ ≤ 0.5 pass, σ ≥ 1/disc/prior/1-px fail); A0 hole recall 0 in every bin; centroid disc 0; exact answers 1.0; probe geometry reproduces B3 A2's KS and bin-share checks (seeds 7 and 8 of v3e, with the recorded P5.2 residual); stratified bootstrap equals `power_c2.py`'s on a fixture; 0 shared patients between probe-train and probe-test patches; `topo_any_unanimous` no longer read by `phase2_entry`; gate table tests incl. a model dropped from LIDC; G-vision on a tiny random VLM fixture: a model compared with itself gives ratio 1, a projector forced to a constant gives ratio 0 and TV 0 (fail), TV equals a numpy reference and lies in [0, 1], image and item sets fixed by hash | — |
| M5 | **Sampling bridge**: `train.sample` export (chat text, image path, kwargs) / import (ids, texts, finish_reason → `truncated` by TRL's rule); `experiments/phase2/sample_vllm.py` (runs in `/opt/vllm`, no objdel import) | train | `lab/B5/export_prompts.py`, `vllm_parity.py`, `hf_teacher.py` | JSONL schema round trip; truncation rule equals `rewards.truncated`; `processor_kwargs` minus `padding_side`; explicit SamplingParams echoed; sampler module imports nothing from objdel (static check) | M2 |
| M6 | **NTL loss**: **removed** (revision 3); built only if E2's trigger fires (§6.4) | — | A2 `ntl_proto.py` | (if triggered) λ = 0 reproduces TRL's fused NLL to 1e-4; one-hot correct digit → 0; gradients flow; LFM value mapping | M2 |
| M7 | **Vast ops**: scripts (§4.5), `vast/policy.py` (offer filter, governance, ledger), `Dockerfile.phase2/phase3` (with the Qwen3.5 GDN extras, D15), `build_images.sh`, `image_check.sh` (records the GDN route per function) | ops | A3 lab, `lab/B5/Dockerfile.*` | H-GOV: refuses 100 % of recorded `hosting_type` 0 offers for RIGA, allows datacenter; RIGA decrypt reads the age identity from stdin only (no env var, file or onstart field carries it; dry-run check); ceiling and `gpu_ram` filters; ledger has no secret fields; dry-run rsync excludes `data/raw` | — |
| M8 | **Phase 2 scripts**: `train_sft` (TrlParser YAML, WSD branches, run manifest, G-vision check at each rolling checkpoint, abort on both (a) and (b), D7; exits before any sampling), `select_ckpt` (D12 incl. G-vision; separate process; launches the sampler), `g_probe` (§5.1), `k3_harness` (promoted), `probes`, `baselines` (MedSAM/SAM 2.1 on frames, nnU-Net export), `c2_analyze` (stratified bootstrap, seed checks), `prereg.py` (C1 stop rule, D12 rule R-B with the sd_geom veto, G-prior (a) threshold only if the owner gates it before V1, G-probe, C2 outcomes, §6.3 fallback rule, C2 rule) | coordinator + train/eval | `lab/B4/k3_harness.py`, `lab/B1/analyze.py`, A4 `medsam_check.py` | `--help` for every script (`test_experiment_scripts`); a 2-step CPU smoke of `train_sft` → `select_ckpt` on fixtures, asserting that `select_ckpt` launches the sampler only when `torch.cuda.is_initialized()` is False in its own process and that `train_sft` has exited; `sample_vllm` refuses to start when `mem_get_info()` free < 90 % (mocked); `train_sft` aborts with a recorded reason only when the G-vision ratio is < 0.5 and TV < 0.002, and only warns when one is (mocked); `c2_analyze` CI equals a reference bootstrap; `select_ckpt` applies R-B to fixture branch panels (deepest wins, veto fires only per D12) | M1–M5, M11 |
| M11 | **Signal panel** (B2): `objdel.eval.signal` with the per-group metrics r_std, σ-share, mixed, headroom, distinct / duplicates, legal-but-wrong, active and mean \|A_v4\|, the exact variance split Var(R) = p(1−p)·μ_L² + p·Var_L(R) and sd_geom, the paired nodule-cluster bootstrap, the R-B selection and veto, and the per-item Spearman with bootstrap CI; called by `select_ckpt` on the samples it already draws and by `g_probe` | eval | `lab/B2/score.py`, `posthoc.py` (`parts`), `analyze.py` | split identity to 1e-12 on random groups; sd_geom = 0 when all legal samples are equal and = r_std when all are legal; metrics equal B2's `score.py` on a stored-sample fixture; R-B: ties stay deep, the veto needs CI > 0 **and** ≥ 15 % **and** no significant IoU loss; Spearman equals scipy's | M4 |
| M9 | **I45 fixes**: `prepare` checks `cv2.imwrite`; `cache_loader` reads the manifest; `fetch_lidc_pixels` refuses a shared z key and HU truncation; `eval_run` handles single-rater items | data/train/eval | ISSUES I45 | one regression test each | — |
| M10 | **Docs**: ROADMAP Phase 2 / Compute / GPUs updated from R0 / V1, including ROADMAP lines 165–176: the SmolVLM-256M / SmolVLM-500M / Qwen3-VL-2B ladder with "about 1 GPU-hour each" (→ D1), "IoU near MedSAM" as a gate (→ reported, §5.1), SmolVLM as the non-Qwen control (→ Mac control only, D1), "each epoch's target is a random expert's polygon" (→ all voting raters every epoch, D9), "flips, 90° rotations, scale jitter" for both datasets (→ D10: D4 LIDC, hflip RIGA, crop shift / scale, CT window); I41 stale numbers; RESULTS for Phase 2 | coordinator | — | `test_readme` still passes | V1 |

Graph: M1, M3, M4, M7, M9 start in parallel; M2 needs M1; M5 needs M2; M11 needs M4; M8 needs M1–M5 and M11. **V1 needs M1–M5,
M7, M8 and M11** (bundle views, reward pool, probes, signal panel, sampler bridge, scripts, image). Every
milestone ends with `make test`, `make lint`, `make cov` ≥ 90 %.

## 8. Risks, kill/continue criteria, owner questions

### 8.1 Risks

| Risk | Mitigation |
|---|---|
| No rung passes G-geometry on shifted LIDC (the 0.60 gate was set on centred frames; human LOO 0.768) | report per-rung shifted numbers; kill/continue rule below (§8.2: Qwen3-VL-4B or InternVL3.5-1B as extra rungs after the 2B rung) |
| Legality at ≥ 90 % is knife-edge for small nodules (σ ≥ 0.6 px vertex error suffices to fail; A4 S4) | §6.3 fallback (rule 3 / rule 2); E2 / E10 triggers (§6.4); precision probe; D12 selection on T = 1 legality |
| Shifted windows include another outlined nodule (ground truth then misses it) | `meta.other_nodules` check in the shift draw; M1 test |
| Qwen3.5 GDN kernels unresolved for torch 2.14.1 / cu130 (`causal-conv1d` needs a matching prebuilt wheel; D15), leaving the reference PyTorch path (on the Mac 79 s/step fp32, 34 s/step with a bf16 base; not representative of CUDA, B1 §2.7), and the CPU-only vLLM parity gap (B5); vLLM #34948 is fixed since 0.27.1 | `kernels` hub route as fallback; route recorded at image build and in R0 step 3; the R0 Qwen3.5 arm prices whichever route resolves before any Qwen3.5 R1; HF-only sampler path (D13); Qwen3-VL-2B stays as the 2B rung |
| Vision collapse during SFT (B1: all-linear LoRA erased Qwen3.5's and SmolVLM's image features within 36 steps, B1 §2.5) | LoRA on the language model only, full-FT vision at 0.1× (D7); in-training G-vision check (abort only when (a) and (b) both fail); G-vision on every branch; D7 repair (frozen) once per rung |
| LFM's legality drops on shifted crops (71 % vs 83 % centred at 150 Mac steps; B1 §2.3) | crop-jitter training (D10); G-format judged on the shifted val; Qwen3.5 is the next rung, its cost measured in R0 |
| vLLM/HF parity fails for a model | HF sampler is the gate path (slower; A5 H-A5-8) |
| Reward pool or fork misbehaves on Vast (threads at fork, OpenBLAS/libgomp) | pool before CUDA, `cv2.setNumThreads(0)`, `ProcessPoolExecutor` only; C++ scorer trigger (D14) |
| Over-trained SFT collapses diversity (no GRPO signal) | WSD branches; D12 R-B with the sd_geom veto; next-earlier-branch repair; R2's check of the veto statistic (B2 found no over-training at Mac depth, but its checkpoints were under-trained, B2 §6) |
| Holes not learnable at 1.35 % prevalence | ×4 oversampling in the bundle; R-ctrl capacity control; C2 fail branch: holes left to GRPO |
| RIGA permission never arrives | Phase 2 exits on LIDC; RIGA arms (band claim) wait; Mac MPS fallback for RIGA SFT only |
| Phase 3 re-pinned env drifts (vLLM 0.33 moves torch to 2.14.0 exact pin, not our 2.14.1) | re-run fast suite + strict hashes per re-pin; re-evaluate at vLLM 0.33 (B5 §What this changes) |
| Blackwell stack issues (sm_120, CUDA 13, driver ≥ 580) | A100 80 GB fallback; V1 on Ada |
| Host failure mid-run | pulls every 15–30 min; resume tested (H-PULL) |
| Forgotten instance / runaway spend | G1–G10 |
| The probe's smallest-bin holes are ~1 px larger than real ones (P5.2 failed narrowly) | stratified design confines the bias to the ≤ 10 px bin; per-bin floor in C2; the bin is reported with the caveat; the real val slice is reported via TV against the raters' hole fraction |
| G-probe is noisy at 20 steps / 128 items | paired per-item design (§5.0); one re-probe after the D12(c) repair; the training-log form reported alongside |
| One bundle, no per-component A/Bs: a harmful component is not isolated in advance | the §6.3 leave-one-out rule targets the most suspicious component per failed gate; §6.4 triggers bring back the specific A/B only when its failure mode appears |

### 8.2 Kill / continue

- **Precision:** a rung failing zero-shot and capacity G-precision goes to §6.3 (if no diagnostic was used yet), else is
  dropped from LIDC (RIGA-only if permitted).
- **Screening:** if no affordable rung passes G-format on A0′ (ROADMAP kill rule "A0′ < 90 % valid"): (1) E10 or E2 per
  their triggers, with owner approval; (2) Qwen3-VL-4B or InternVL3.5-1B as extra rungs after the 2B rung; (3) E8 coordinate tokens;
  (4) a coarser grid G = 100 for **LIDC only** (lossless on LIDC: A1 §1.5; rejected for RIGA). Mask tokens are rejected
  (no polygon, no certificate; A1 §6). If a rung passes format but not G-geometry on the shifted / rescaled val: continue
  to Phase 3 with it only if it passes G-prior (b), G-signal and G-probe, and record the geometry miss as a Phase 3
  baseline (RL adds +9.7 / +9.8 / +10.5 gIoU over SFT in SimpleSeg Table 3 and +11.9 mIoU in Hi-GAR; A1 §3); otherwise stop and re-plan.
- **After G-probe:** if the selected branch and the D12(c) repair branch are both not PASS, the rung still passes
  C1 with a G-probe warning handed to Phase 3 (G-signal holds by D12 eligibility); the ladder climbs only on a §5.1
  gate failure. Phase 3 then starts with its own LR choice and a longer pilot before committing arms.
- **After C2:** a partial or fail outcome does not block Phase 3 (all Phase 3 arms start from A0′); each outcome has at
  most one conditional run (§5.2: partial → ×8 seed, fail → R-ctrl); holes are then a GRPO question, or a capacity
  question if R-ctrl passes.
- **Budget:** $10 alarm → coordinator review; $15 hard stop; the Qwen3-VL-2B rung, R-ctrl past the stop, E10, RIGA and
  any triggered A/B each need the owner's approval of a tranche.

### 8.3 Open questions for the owner

1. **RIGA permission** (I44, H-A4-16): send the emails now: RIGA depositor (almazroa@med.umich.edu, per the readme) and
   the MESSIDOR partners; register at ADCIS. Until a written yes, RIGA stays off Vast (D16).
2. **Commercial use**: is this project done for a company? RIGA and MESSIDOR forbid "unauthorized commercial use"
   (A4 §5.2); LFM Open 1.0 allows commercial use only under $10M revenue (A1 §2); MedSAM2 weights are research-only.
   The answer can remove RIGA and LFM from the plan.
3. **Vast funding**: approve $25 of prepaid credit under the **$15 hard stop** (alarm $10); card (deletion protection, but auto-charges) vs no card
   (hard cap, destroy at $0) (A3 G1).
4. **GHCR**: create a private GHCR package for the image and a read-only `read:packages` token for Vast's `--login`
   (V10). (CI already ran green on GitHub, run 39359008526; I34's CI half is closed.)
5. Tracking: trackio (default) or an existing W&B account.

## 9. Decisions log and open items

### 9.1 Log

**2026-10-20, draft C.**
- Adopted A1's ladder as the default (D1), pending B1. ROADMAP's ladder (SmolVLM-256M/500M, Qwen3-VL-2B) and P1-D12's
  are superseded once B1 lands.
- Gates revised (D18): shifted val, prior paired win, re-specified precision probe with the centroid-disc paired win
  (B3, exploratory threshold registered here), synthetic topology probe. "Near MedSAM" reported, not gated.
- RIGA needs written permission **and** Secure Cloud (RESULTS I44), stricter than A3/A4's "Secure Cloud by default,
  email in parallel".
- Deviations from the notes, with reasons: label smoothing rejected (A2 over A4); vLLM 0.32 sampler (B5 over A3's
  ≤ 0.31); server mode on one GPU rejected (B5 over A5 §8.4); `cpu_cores_effective ≥ 12`, `cpu_ram ≥ 32/48` (§4.2 over A3
  §3.3); Phase 2 budget (§4.9) replaces ROADMAP's ~6 GPU-h.
- B1 (D1 rule) and B2 (D12 rule) were pending; both rules were written so the result fills them without re-planning
  (resolved in revisions 5 and 6).
  (Revision 4 note: SmolVLM is a Mac control only and never screened on Vast; B2 informs the depth that C2's runs train to.)

**2026-10-20, revision 1 (after JUDGE_R1).** C2 rescoped to holes (ROADMAP Phase 2 wording changed to match);
G-probe added; statistics fixed in §5.0; baselines on matching views; GPU ownership sequenced (train → select →
sample → probe, one process at a time); evidence corrections (Quagmires direction, B3 P1/P2 labels, D5 epoch-equivalent
costing, B1 AMEND 2/3); critical path trimmed with triggers (§6.0); rescaled val view; G-topology (b) replaced; E10
rejection-sampling fallback; CI closed (run 39359008526). New marker: PENDING B3-addendum (closed in revision 2).

**2026-10-20, revision 2 (B3 Addendum A).** PENDING B3-addendum resolved: real hole geometry (median 12 px; bins
≤ 10 / 11–30 / > 30 px = 44 / 36 / 20 %; 75 % single-rater; none unanimous); probe = 120 one-hole items per bin + 160
two-part + 160 simple; primary metric = mean over K = 8 at T = 1 averaged over 2 seeds, items bootstrapped within
bins (pooled ±3.1 pts); seed caveat adopted (item CI conditional on the two checkpoints → both seeds' point estimates
must meet the threshold, discordance flag); C2 rule copied from B3 Addendum A6; generator fidelity P5.2 recorded as
**failed narrowly** (≈ 1 px area under-draw of 5-px holes in the smallest bin; one of two fresh-seed confirmations
failed, iteration stopped) with the reasons the stratified probe remains fit (D20). C2 is now a one-sample threshold
test on A0′'s hole recall; the paired A0′ − A0 test is dropped as uninformative (A0 = 0 by construction). E1 uses a
probe draw disjoint from C2's.

**2026-10-20, revision 3 (owner cost directive).** "Reduce Phase 2's Vast work to the critical questions and test
components combined, not one at a time." Changes: one evidence-chosen recipe bundle in every run (§6.0); number-token
loss left out (D8: text-only evidence, untested on VLM coordinates, weakest on LFM's tokens); hole oversampling fixed
at ×4; WSD reduced to 2 stable epochs with branches 0.5 / 1 / 2 (D5); Vast runs only for C1 (smallest first, **stop
at the first passing rung**), C2 (2 seeds, reusing the screening run as A0′ seed 1) and G-probe on the selected
checkpoint (§6.1); one pre-registered leave-one-out diagnostic as the only fallback (§6.3); throughput, loader,
determinism cost, K3/K4 and sampler parity folded into R0–R2 (§6.2); baselines on the Mac CPU; Qwen3-VL-2B only as the
next rung or the conditional capacity control R-ctrl; budget recomputed: critical path 3.5–7.2 GPU-h, $2.3–4.7
(≤ $7.3 at price ceilings), alarm $10, hard stop $15 (§4.9, with old vs new). **Removed or triggered** (each kept with
its trigger in §6.4): E1–E10, precision capacity runs (conditional), the unconditional 2B control and larger-rung
screening, nnU-Net training, K = 64 / I35 (Phase 3 entry), HD95/ASSD/centroid, H-A5-4/5, the Phase 3 environment smoke
(Phase 3 start), H-PULL, RIGA runs. Judge items in §9.2 that named removed A/Bs (SHOULD-FIX 1, 7, 11) remain addressed:
their tests now sit in §6.4 with triggers.

**2026-10-20, revision 4 (judge round 2).** G-topology (b) report-only for C1; G-probe calibrated (step-0 null,
pre-registered probe LR, no ladder climb on G-probe alone); C2 reallocated to 3 A0′ seeds + 1 A0 seed with a 3-seed
rule; SmolVLM control-only; int16 HU crops for CT window jitter (M1, byte-identity test); vision tower at 0.1× LR;
evidence corrections; budget 3.7–7.5 GPU-h / $2.3–4.9 for scenario A, ≈ 13 GPU-h / ≈ $8.5 at MFU 0.05. Details: §9.2b.

**2026-10-20, revision 5 (B1 report, `lab/B1/REPORT.md`).** PENDING B1 resolved. D1: H-A1-2 and H-A1-1 are
inconclusive under their pre-registered rules, so the default ladder LFM → Qwen3.5 → Qwen3-VL-2B stands; smallest first
is kept over B1's side-by-side recommendation (a full side-by-side costs +1.7–4.2 GPU-h [superseded: 1.7–3.9 since revision 7] and passes the $10 alarm at
MFU 0.05), and a 30-step Qwen3.5 arm in R0 takes its place (+0.2 GPU-h, +$0.13); SmolVLM is a control only. D6 / D7 /
M2: vision tower and projector frozen under LoRA (`configs.lora` becomes LM-only), 0.1× LR under full FT; new
G-vision gate (§5.1, thresholds from B1 §2.5) on every branch, an in-training projector check and the D7 repair;
`eval/vision.py` (M4). D10 / §6: crop-jitter training is mandatory (B1: 94–98 % of shifted samples at the frame
centre); E3 closed for good. D13: explicit sampling overrides for both cards; `MAX_COMPLETION` 1,024 for Qwen3.5
(M2 had 1,536); H-A1-8 passes for Qwen3.5 and is implied by G-format + G-signal for LFM. D15 / §4.3: Qwen3.5 GDN
kernel routes, UNVERIFIED for torch 2.14.1 / cu130. D2: no reliance on a zero-shot prior. §4.9: Qwen3.5 stays 2× LFM
(equal answer lengths; Mac s/step not used), UNVERIFIED until R0; scenario A 3.9–7.7 GPU-h / $2.5–5.0 (≈ 13 GPU-h /
≈ $8.6 at MFU 0.05).

**2026-10-20, revision 6 (B2 report, `lab/B2/REPORT.md`; wave B complete).** PENDING B2 resolved; no PENDING marker
remains. D12: rule R-B (deepest eligible branch; signal veto only when the shallower branch's sd_geom is significantly
and ≥ 15 % higher with IoU not significantly lower), replacing the lexicographic R-A, whose first keys B2 found
saturated. Veto statistic: sd_geom, post-hoc and labelled so, chosen over the pre-registered r_std because r_std's
legality part rewards less legal branches; R2's per-item Spearman confirms or reverses it for later rungs [superseded in revision 7: it can only retain sd_geom or flag]. D12(c):
next-earlier eligible branch; base ↔ LoRA interpolation rejected (H-B2-2); WiSE-FT between two branches only as a
report-only inference check. §6.2: signal panel, stable-phase checkpoints, Spearman, shallow vs deep, entropy proxy.
New M11 (`objdel.eval.signal`, from `lab/B2/score.py` + `posthoc.py`). D8: DFT is H-A2-8 and was not run in B2 (the
plan had called it B2's negative control). §4.9: no run and no critical-path cost added (conditional WiSE-FT check
≈ 0.05 GPU-h). Consistency pass: §6.1 R1 estimate aligned with §4.9 (0.76–1.70 h); G-probe repair wording aligned
in §3, §5.1, §6.1, §8.2; M8 dependencies and script names (`c2_analyze`) aligned with the graph.

**2026-10-20, revision 7 (judge round 3, `JUDGE_R3.md`: APPROVE WITH CHANGES).** MUST-FIX 1: G-prior (a) (centred
win share, uncalibrated 60 %) is report-only for C1; image use is decided by G-prior (b), G-geometry on both views and
G-vision; the Mac baseline run computes MedSAM-oracle's and SAM 2.1-tiny's centred shares, and a gated (a) would need a
threshold ≤ SAM 2.1-tiny's − 5 pts in `prereg.py` before V1. SHOULD-FIX 1–8: §5.0 per-seed clauses with `prereg.py`
authoritative; in-training vision abort only when G-vision (a) and (b) both fail, else warn; R2's Spearman can only
retain sd_geom or flag, a switch is Phase 3 planning; C2 outcomes pass / partial / fail with one conditional each
(new R-x8 on partial, R-ctrl on fail); C2 depth = the checkpoint with the final G-probe verdict; stale names fixed;
M10 lists ROADMAP lines 165–176; Qwen3.5 conditional row 1.7–3.9 h (scenario B 7.2–16.3 GPU-h / $4.7–10.7). §9.2c.

**2026-10-20, revision 8 (judge round 4: APPROVE).** Optional polish only (§9.2d); no GPU cost, no decision change.

### 9.2 Judge round 1 (`JUDGE_R1.md`): how each item was addressed

| Item | Addressed where | How |
|---|---|---|
| MUST-FIX 1 (C2 hole-only, probe uncalibrated) | §1 C2, D20, §5.2 C2, M4, ROADMAP Phase 2 | C2 and ROADMAP "done when" rescoped to holes (A0-multi keeps parts); probe hole area / part ratio drawn from the 414 size-rule non-simple train targets with the size-rule minimum; hole recall per bin (≤ 10, 11–30, > 30 px) with "no bin below 0.30" in C2; recalibration done by B3 Addendum A (120 per bin; P5.2 failed narrowly, recorded in D20); 44-item val slice stays descriptive (TV vs rater hole fraction) |
| MUST-FIX 2 (no "GRPO can improve" KPI) | §5.1 G-probe, D12, §3.2 ④, §8.2 | G-probe: 20 GRPO steps, step-0 vs step-20 paired per item on 128 fixed shifted-val items, CI > 0 on the raw composite, legal drop ≤ 5 pts, σ-share ≥ 50 % at step 20; pre-registered in `prereg.py` (M8); D12(c) repair fires on G-probe failure. Deviation: paired eval on a fixed item set instead of training-step windows 11–20 vs 1–5, because TRL's aggregator is group-mean-centred (mean ≈ 0) and training prompts differ per step; the judge's form is reported |
| MUST-FIX 3 (seed pooling, primary metric) | §5.0, D20, §6.1 | per-item mean over seeds and K, greedy separate; item/nodule bootstrap; pooled + per-seed CIs, decision on pooled, seed-discordant flag and its consequence; E-tests are uncorrected decision rules; C2 the single confirmatory test; topology primary metric fixed by B3 Addendum A5 (mean over K = 8, 2 seeds, within-bin bootstrap) with the seed caveat (both seeds' point estimates must meet the threshold) |
| MUST-FIX 4 (baselines on the gate view) | D17, §5.1 "near MedSAM" | every row labelled with its view; prior, MedSAM, SAM 2.1 on centred, shifted (and rescaled for the prior) views with the box from that view's union; nnU-Net trained on the jittered train view and evaluated on shifted val |
| MUST-FIX 5 (trainer vs vLLM on one GPU) | D13, §3.2, M8 | `train_sft` exits before any sampling; `select_ckpt` is a separate, CUDA-free launcher; the sampler asserts free memory ≥ 90 %; G-probe is a later process; M8 smoke asserts all three |
| MUST-FIX 6a (Quagmires reversed) | D5, inconsistency table | abstract re-read (arXiv 2510.01624): one epoch on unique data **underperforms** two epochs on half; "unique > repeated" dropped as a reason; fresh views per epoch justified by the crop leak only |
| MUST-FIX 6b (B3 P1 vs P2) | §1 C2, D20, §5.1, §5.2 | A0 0.50 [0.45, 0.55] cited to B3 REPORT §3 P2; census and cluster CIs stay B3 P1 |
| MUST-FIX 6c (2 vs ~4 epochs) | D5, §4.9 | revision 3/4: D5 = 2 stable epochs + branches 0.5 / 1 / 2 ≈ 2.5 LIDC epoch-equivalents; §4.9 scales `estimate_cost.py` by 2.5 × 0.86 / 2 = 1.075 and × 1.053 for ×4 rows (1.13), plus branch eval and G-probe |
| MUST-FIX 6d (B1 AMEND 2/3) | D1 | AMEND 1–4 recorded (AMEND 4: Qwen3.5 joins the 150-step tier), with their consequences for the B1 rule; SmolVLM is now control-only (revision 4) |
| SHOULD-FIX 1 (scope) | §6, §3.3, §4.9, §1 | revision 3 (owner directive) supersedes the rev-2 critical path: one bundle; Vast only for R0 (folded measurements), R1 screening with stop at the first passing rung, R2 G-probe, R3–R5 C2; one §6.3 diagnostic; everything else in §6.4 with triggers |
| SHOULD-FIX 2 (G-topology (b) vacuous) | §5.1 | (b) is now hole recall ≥ 0.30 on one-hole items (the C2 fail line) |
| SHOULD-FIX 3 (scale cue in val) | §5.0, §5.1, D18 | rescaled val view (S moved to a neighbouring divisor + shift); G-geometry and G-prior (b) on both shifted and rescaled |
| SHOULD-FIX 4 (unneeded preload) | D14, §3.2, §4.2, M3 | caches preloaded per dataset and only by scoring processes; `train_sft` loads none; LIDC-only boxes `cpu_ram ≥ 32` |
| SHOULD-FIX 5 (RIGA key delivery) | D16, §4.4, M7 | age identity streamed over the SSH session's stdin into `age -d -i /dev/stdin`; never env var, template, onstart field or file; ledger records the public recipient only; tested in M7 |
| SHOULD-FIX 6 (probe patient disjointness) | D19, M4 | probe-train and probe-test backgrounds from patient-disjoint halves of LIDC train (seeded patient hash), tested; topology probe on val patients |
| SHOULD-FIX 7 (rejection-sampling SFT) | D28, E10, §6.0 | researched (RFT 2308.01825, ReST-EM 2312.06585; abstracts read); pre-registered E10, triggered by legal in [85, 90) % at every branch, with a diversity guard (G-signal, G-probe) |
| SHOULD-FIX 8 (D12 shallow-branch bias) | D12 | revision 6: rule R-B selects the deepest eligible branch; a shallower one wins only on a significant and ≥ 15 % higher sd_geom with no significant IoU loss (B2 §4), so noise cannot pick a shallow branch (the revision-1 lexicographic rule with minimum deltas is superseded) |
| SHOULD-FIX 9 (G-format reporting) | §5.1 G-format, D25 | truncation share, mean / p95 completion tokens and held-out NLL per branch reported (not gated) |
| SHOULD-FIX 10 (Qwen3.5 sampler) | D13, §4.9, §3.2 | HF-only path pre-decided on R0 parity/kernel failure: branch selection on 256 items × K = 8, full val only for the selected branch; cost row added (0.7–3 h, ≤ $2, UNVERIFIED decode rate) |
| SHOULD-FIX 11 (S7 scope) | §6.4 | revision 3: the Phase 3 environment smoke (S7) moved to Phase 3 start; E8 removed with its trigger (generation > 60 % of the GRPO step in R2) |
| SHOULD-FIX 12 (minor) | D1, D18, §3.5, M4 | B1's looser shift draw (0 px clearance, no `other_nodules` check) recorded as a caveat; M4 replaces `topo_any_unanimous` in `eval/gates.py` rather than adding beside it |

### 9.2b Judge round 2 (`JUDGE_R2.md`): how each item was addressed (revision 4)

| Item | Where | Resolution |
|---|---|---|
| MUST-FIX 1 (G-topology (b) as a C1 gate) | §5.1 G-topology, §5.2 C1 | (b) hole recall is report-only for C1; decided by C2 only; R1 reports it early for R-ctrl |
| MUST-FIX 2 (G-probe uncalibrated) | §5.1 G-probe, D12, §3.3, §6.1, §8.2 | step-0 double sample gives the null CI at zero training cost; PASS / FAIL / INCONCLUSIVE defined against it; pre-registered probe LR 5e-6 full FT / 5e-5 LoRA (5× production, 10× below the SFT peak; SimpleSeg RL 2e-6 constant verified in its Tables 4–5), recorded as the probe's, not Phase 3's; a non-PASS after the repair hands C1 to Phase 3 with a warning and never climbs the ladder (G-signal is already required by D12) |
| MUST-FIX 3 (A0 arm buys no decision) | §5.2 C2, §5.0, §6.1, §3.3, §4.9 | A0 seed 1 + A0′ seeds 1–3, same 3 runs; 3-seed rule (pooled CI > 0.30 and pooled ≥ 0.50; ≥ 2 of 3 seeds ≥ 0.50, none < 0.30; bins, controls, legality unchanged); reconciled with B3 A5 (item CI for item noise, per-seed clause for seed noise, across-seed t-interval reported); A0 guard pre-registered report-only (paired CI > −0.02). ROADMAP has no Phase 2 seed count, so it is unchanged |
| MUST-FIX 4 (SmolVLM as C1) | D1, §6.1, §4.9 | SmolVLM is a Mac control only, never screened on Vast and never C1 (TRL #6294 open, `use_vllm=False`, excluded from TRL GRPO CI) |
| MUST-FIX 5 (CT jitter needs HU) | D10, §6.0, §4.4, M1 | M1 writes int16 HU context crops (npz) from `data/store`; byte-identity test vs the prepared PNG at zero jitter; bundle size ≈ 0.3–0.5 GB compressed (estimate); 2510.08116 cited as random windowing for liver-tumour CT, the weakest component |
| MUST-FIX 6a–j (evidence) | D13, §8.1, D2, D1, D5, §8.2, G-format, D20, §9.2, §8.3 | (a) vLLM #34948 fixed since 0.27.1; live risks are the CPU parity gap and `fla` / `causal-conv1d` wheels; (b) #6294: 855 vs 38 (issue) and 89 → 155 (A5 sandbox); (c) G = 100: −2.7 IoU / −2.2 legality, G = 250: −1 / −1.5; (d) Qwen3.5 card has no `bbox_2d`: UNVERIFIED, B1 T1; (e) SimpleSeg SFT 5e-5 cosine → 2e-6, RL constant 2e-6, RL gain +9.7 / +9.8 / +10.5 gIoU (Table 3; the paper's KL is 0.01 in the text and 0.1 in Table 5); (f) A5 §0.6; (g) size rule cited to the reward code; (h) §9.2 rows 6c, 6d, SHOULD-FIX 1, 11 rewritten; (i) AMEND 4 in D1; (j) $15 stop in §8.3 |
| SHOULD-FIX 1 (controls 90 vs 95 %) | §5.1 G-topology | (a) raised to ≥ 95 % |
| SHOULD-FIX 2 (§6.3 depth) | §6.3 | deepest branch when no branch is eligible |
| SHOULD-FIX 3 (R-ctrl one seed) | §6.1 R-ctrl | exploratory pass on clauses (1), (3), (4), (5) of one seed |
| SHOULD-FIX 4 (B2 not started) | D12 | B2 done (revision 6): rule R-B; repair = next-earlier eligible branch |
| SHOULD-FIX 5 (vision-tower LR) | D7, §6.0 | 0.1× LR adopted in the bundle, justified (A2 §4.2's own proposal; conservative H-A2-6 arm); E4's trigger now tests 1× vs 0.1×; revision 5 (B1 §2.5): frozen under LoRA, 0.1× kept for full FT (D7) |
| SHOULD-FIX 6 (free attribution check) | §6.1 R1 | identity vs rotated view on the selected branch, inference only, report-only |
| SHOULD-FIX 7 (budget realism) | §4.9, §4.8 | MFU 0.05 row (≈ 13 GPU-h, ≈ $8.5), ×4 rows counted (+5.3 %), eval block marked UNVERIFIED, A100 quote filter stated |
| SHOULD-FIX 8 (rescaled coverage) | §5.0 | share with no admissible neighbour reported; > 30 % flagged as a weaker check |
| SHOULD-FIX 9 (LoRA batch caveat) | D6 | recorded |
| SHOULD-FIX 10 (LIDC prior) | §5.1 G-prior | LIDC centred prior 0.588 / 65.3 % (fails narrowly) cited |

### 9.2c Judge round 3 (`JUDGE_R3.md`): how each item was addressed (revision 7)

| Item | Where | Resolution |
|---|---|---|
| MUST-FIX 1 (G-prior (a) uncalibrated C1 gate) | §5.1 G-prior, §5.2 C1, D17, D18, §6.3, §8.2, M8 | option (i): (a) report-only, computed on the K = 8 item mean with illegal = 0; (b), G-geometry and G-vision decide image use; D17's Mac run computes MedSAM-oracle and SAM 2.1-tiny centred shares; a gated (a) needs a threshold ≤ SAM 2.1-tiny − 5 pts in `prereg.py` before V1, never after C1 |
| SHOULD-FIX 1 (§5.0 "both seeds") | §5.0 | per-seed clauses of the decision rule (C2: ≥ 2 of 3 ≥ 0.50, none < 0.30); `prereg.py` authoritative |
| SHOULD-FIX 2 (in-training abort) | D7, §5.1 G-vision, §3.2, §3.3, §6.1, §4.9, §8.1, M8 | abort only when (a) ratio < 0.5 and (b) TV < 0.002 both hold; one alone is warn-only, left to `select_ckpt` |
| SHOULD-FIX 3 (Spearman as a switch) | D12, §6.2, §9.3, §9.4 | confounded by legality headroom: retains sd_geom or flags only; a switch needs the branch-level comparison and is a Phase 3 planning decision |
| SHOULD-FIX 4 (C2 outcome taxonomy) | §5.2, §6.1 (R-x8, R-ctrl), §6.4 E1, §3.3, §4.9, §8.2, §1 | pass → none; partial ([0.30, 0.50), or ≥ 0.50 with a failed clause) → R-x8 only (not if clause (4) failed); fail (< 0.30) → R-ctrl only |
| SHOULD-FIX 5 (C2 depth after a repair) | §5.2, §6.1, §3.3 | A0′ seed 1 and the C2 depth are the checkpoint that received the final G-probe verdict |
| SHOULD-FIX 6 (stale names) | M10, §5.3, §9.2, §9.3, D13, D1, §8.1, §8.2, §5.1 G-precision | S0 → R0 / V1; S7 → Phase 3 start (§6.4); K = 64 at Phase 3 entry; 4B / InternVL3.5-1B as extra rungs after the 2B rung; G-precision "zero-shot, then capacity run" |
| SHOULD-FIX 7 (ROADMAP M10 scope) | M10 | ROADMAP lines 165–176 listed item by item with the plan decision that replaces each |
| SHOULD-FIX 8 (Qwen3.5 row) | §4.9, D1 | 1.7–3.9 h ($1.1–2.6); scenarios B / C and the side-by-side figure recomputed |

### 9.2d Judge round 4 (`JUDGE_R4.md`): APPROVE (revision 8)

Revision 7 approved with no MUST-FIX or SHOULD-FIX. The six optional polish items were applied, with no GPU cost and
no decision change: (1) §4.9 scenario C 8.0–18.0 GPU-h / $8.9–21.4 with its build-up and the 2B C2 probe-eval note;
(2) §6.1 R-x8 exploratory pass rule (clauses (1) on its seed, (3), (4), (5)) on a disjoint topology-probe draw;
(3) §5.0 fourth seed is part of the C2 verdict (≥ 3 of 4) and precedes R-x8 / R-ctrl; §6.2 "third" → "fourth";
(4) §6.1 R-ctrl approval aligned with §5.2; (5) §3.3 partial label and the §3.3 / §4.9 headers; (6) log hygiene
(superseded figures annotated, revision entries in order).

### 9.3 Inconsistencies between the notes (resolved here)

| Topic | Notes disagree | Resolution |
|---|---|---|
| Quagmires (2510.01624) | A2 §4.10 and §6.1 say "one epoch on unique data beats more epochs on half the data"; the abstract says the opposite ("…underperforms training on half examples for two epochs") | the abstract (D5); A2's paraphrase is wrong (JUDGE_R1 MUST-FIX 6a) |
| SimpleSeg RL settings | A2 §4.3 and revision 1–3 cited "SimpleSeg LR 5e-5" without saying it is the SFT LR; the paper's RL LR is a constant 2e-6, and its KL is 0.01 in the text but 0.1 in Table 5 | SFT 5e-5 → 2e-6 cosine, RL 2e-6 constant (Tables 4–5); the KL inconsistency is the paper's own |
| Label smoothing | A4 §4.3/§7.6 recommends 0.1; A2 §4.5 rejects it (deprecated in TRL SFT, extra forward pass, no repair of over-training) | rejected (D8) |
| vLLM version for Phase 2 sampling | A3 §6.1/H-VLLM: ≤ 0.31.0; A5/B5: 0.32.0 (0.31 caps transformers < 5.18) | 0.32.0 (D15) |
| Phase 3 server mode | A5 §3/§8.4 says TRL pushes weights over `ipc` on a shared GPU; B5 read TRL 1.16's client: NCCL-only, no IPC, separate devices required, PR #7799 KeyError | server mode rejected on one GPU (D15) |
| Model ladder | A3 and A2 plan around SmolVLM-256M → LFM → Qwen3-VL-2B (P1-D12); A1 proposes LFM → Qwen3.5-0.8B → Qwen3-VL-2B; `estimate_cost.py` has no Qwen3.5 entry | A1 default, kept after B1 (inconclusive → default, D1); Qwen3.5 cost marked UNVERIFIED until R0 (§4.9) |
| Screening order | B1 §0.5 / §4.1: screen LFM and Qwen3.5 side by side in the first Vast session; owner directive: smallest first, stop at the first passing rung | smallest first plus a 30-step Qwen3.5 arm in R0 (D1, §4.9) |
| Vision tower under LoRA | revision 4 (D6, D7) and `configs.lora`: LoRA on all linear layers incl. vision, tower at 0.1×; B1 §2.5: that adapter at 5e-4 collapses Qwen3.5's and SmolVLM's image features | LoRA frozen, full FT 0.1×, G-vision gate (D7) |
| Qwen3.5 `MAX_COMPLETION` | M2 (revision 4): 1,536; B1 §4.5: 1,024 (longest LIDC train target 876 tokens; B1 §1.5 AMEND 1 says 881) | 1,024 for LIDC (both figures fit; truncation ≤ 0.2 %, B1 §2.6) |
| D12 veto statistic | B2 pre-registered r_std and headroom; its post-hoc split shows both fall with legality, and its §4 text proposes sd_geom | sd_geom, labelled post-hoc; R2's per-item Spearman can only retain it or flag disagreement; a switch is a Phase 3 planning decision (D12) |
| WiSE-FT / DFT labels | B2's brief called interpolation H-A2-8; A2 and D12 use H-A2-7 = WiSE-FT, H-A2-8 = DFT; B2 tested only a base ↔ late form of H-A2-7 | A2's labels (B2 §0); H-A2-7 proper (early ↔ late) and H-A2-8 stay untested |
| Qwen3.5 grounding prior | A1 §2: native 0–1000 `bbox_2d` prior; the card does not document it (JUDGE_R2 6d) | B1 §2.1: Qwen3.5 does answer in `bbox_2d`, but the prior is weak (shifted RIGA disc 0.41, LIDC ≤ 0.07, never above the image-blind prior); the plan relies on no zero-shot prior (D2) |
| Non-simple counts | A1: 366 train / 48 val items with a non-simple rater, 23 unanimous train; B3 census (size rule, voting raters): 301 / 44, 12 unanimous. A1 H-A1-6 says "517 non-simple train targets" but its own table gives 317 + 200 − 14 = 503 unique (`census.json` `raw_nonsimple_targets` 503); under the size rule 414 | B3's definitions (they match the reward): 414 non-simple train targets, 301 items |
| Qwen3-VL-2B fine-tuning | A2: LoRA r16 (merged); A3: "Qwen3-VL-2B full SFT" on the PRO 6000 | LoRA r16 default on the PRO 6000, full FT as fallback (D6) |
| Prepared set | A6 profiled on prepared manifest `cd66d79f…`; the final is `3747acbc…` (RESULTS §3) | timings only (LIDC upscale change F3); re-measured on Vast anyway |
| Cache preload | A6: RIGA 4.4 GB, ~35 s serial, miss cost ~2 s/step; B4: 3.35 GB array bytes (Linux RSS ~4.2–4.5 GB), 72 s at load 23 with 2 threads, miss cost 4.2 s | B4's measurements; ≥ 4 threads on Vast |
| Prices | `estimate_cost.py` minimums (3090 0.24, 4090 0.48, 5090 0.56) differ from the snapshot table (0.20, 0.51, 0.63) in A3 §3.2 | §4.8 quotes the snapshot table; per-run $ quoted from the script as estimates |
| Augmentation set | ROADMAP: flips and 90° rotations for all; A4: D4 for LIDC, horizontal flip only for RIGA; B3 built all 8 ops for both | A4 (D10); B3's RIGA variants beyond hflip are unused |
| Mixing | A2: natural proportions, RIGA ×2 only if it misses; A4: T = 2 by items | natural default; T = 2 (E5) triggered only (D11, §6.4) |
| GRPO step in Phase 2 | `configs.grpo_config` sets `use_vllm=True` (colocate) for LFM/Qwen, but the Phase 2 training venv has no vLLM (TRL raises ImportError) | Phase 2 measurements pass `use_vllm=False`; colocate measured at Phase 3 start (§6.4) |
| Precision probe gate | RESULTS §8: ≤ 2 px; A4: ASSD ≤ 0.6 etc.; B3: the centroid disc passes A4's geometric rule | D19 (A4 + B3's paired-win rule) |

### 9.4 Open items

**Carried:** I34 (CI half closed by run 39359008526; native-x86 hashes and K3/K4 timing are R0, §4.7 steps 5–6), I13 (R4 lean; TV/JS on disagreement items at K = 64), I15/I16 (hole adversaries; baseline-relative
TV/retention at K = 64 once real samples exist), I35 (audit recalibration, Phase 3 entry), I36 (rim vertices, M1), I38 (Phase 3),
I41 (docs, M10), I42 (D15), I43 (D10, D18), I44 (D16), I45 (M9).

**New in revision 5 (from B1):** Qwen3.5 GDN kernel route for torch 2.14.1 / cu130 (D15, UNVERIFIED; M7 build, R0
step 3); Qwen3.5's Vast s/step (§4.9, R0 arm); G-vision thresholds rest on small post-hoc diagnostics (4–8 items, B1
§5), so R0's base references are checked against B1's base values before the gate is trusted; B1's Qwen3.5 numbers
are single-seed.

**New in revision 6 (from B2):** whether any signal statistic predicts the GRPO gain (R2 Spearman, §6.2); the
over-trained regime of 2606.18487 (B2's checkpoints were under-trained, B2 §6); WiSE-FT between two SFT branches
(report-only when the repair fires); the D12 veto statistic is post-hoc (R2 can retain or flag it, not switch it).

## 10. Sources

**Research notes** (`docs/research/notes/Phase 2 SFT/`): A1 §0–§7 (models, tokens per vertex, A0 ceiling, ladder,
H-A1-1…10); A2 §0–§9 (TRL facts, sequence budget, collation, recipe, checkpoint selection, H-A2-1…10); A3 §0–§9
(Vast facts V1–V21, operations, snapshot, GPUs, governance, image, guardrails, first session); A4 §0–§9 (S1–S7,
gates sanity, baselines, augmentation, licences, probe); A5 §0–§10 (vLLM vs HF, TRL hand-off, grammar, IS mismatch,
parity hypotheses); A6 §0–§5 (profiles, CPU vs GPU budget, per-module decisions).

**Lab reports** (`experiments/phase2/lab/`): B3/REPORT.md §1–3 (topology probe P1, precision probe P2, dihedral P3,
crop-before-trace P4; Addendum A: hole census, hole-calibrated topology probe, C2 power and rule, files
`hole_census.json`, `topo_probe_v3*.json`, `power_c2.json`; prereg `B3/prereg.py` P5); B4/REPORT.md §0–5 (preload, fork pool, synthetic RIGA, C++ decision;
prereg `B4/prereg.py`); B5/REPORT.md Q1–Q3 (resolution matrix, parity, images; prereg `B5/prereg.py`); B1/REPORT.md §0–§7
(zero-shot priors, learnability at 36 / 150 steps, LM-only arms, vision and centre-bias diagnostics, H-A1-8, Mac cost,
per-model settings; prereg `B1/prereg.py` AMEND 1–5; tools `vision_feats.py`, `image_use.py`, `centre_bias.py`); B2/REPORT.md §0–§9 (depth vs GRPO signal, H-B2-1 inconclusive, base interpolation H-B2-2 rejected, variance split,
D12 rule R-B text, folded measurements; prereg `B2/prereg.py` AMEND 1; `score.py`, `posthoc.py`); A3 `vast_snapshot.json`, `estimate_cost.py`; A4 `null_prior*.json`,
`medsam_check.json`, `rater_ceiling.json`, `noise_model.json`.

**Phase 1**: `experiments/phase1/PLAN.md` §2, §7; `RESULTS.md` §2, §3, §7, §8; `impl/ISSUES.md` I34–I45; `ROADMAP.md`
§Phase 2, §Compute, §GPUs, §Kill criteria.

**Key external sources** (full lists in the notes):
- TRL v1.16.0 SFT/GRPO source and docs: https://github.com/huggingface/trl/blob/v1.16.0/trl/trainer/sft_trainer.py ,
  https://huggingface.co/docs/trl/main/en/vllm_integration ; TRL #6294, PR #6718, #7198, PR #7799.
- vLLM v0.32.0 release and requirements: https://github.com/vllm-project/vllm/releases/tag/v0.32.0 ;
  https://pypi.org/pypi/vllm/0.32.0/json ; torch 2.14.1 metadata https://pypi.org/pypi/torch/2.14.1/json
- Vast docs: https://docs.vast.ai/guides/instances/pricing , https://docs.vast.ai/cli/reference/search-offers ,
  https://docs.vast.ai/guides/reference/faq/security , https://docs.vast.ai/guides/instances/docker-environment ,
  ToS https://vast.ai/terms
- Licences: RIGA DataCite https://api.datacite.org/dois/10.7302/Z23R0R29 ; MESSIDOR
  https://www.adcis.net/en/third-party/messidor/ ; TCIA LIDC-IDRI https://www.cancerimagingarchive.net/collection/lidc-idri/
- Models: https://huggingface.co/Qwen/Qwen3.5-0.8B , https://huggingface.co/LiquidAI/LFM2.5-VL-450M ,
  https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct
- Recipe evidence: LoRA Without Regret https://thinkingmachines.ai/blog/lora/ ; WSD https://arxiv.org/abs/2410.05192 ;
  Quagmires https://arxiv.org/abs/2510.01624 ; SFT rank inversion https://arxiv.org/abs/2606.18487 ; WiSE-FT for
  reasoning https://arxiv.org/abs/2504.10478 ; NTL https://arxiv.org/abs/2411.02083 ; SimpleSeg
  https://arxiv.org/abs/2601.19228 ; LPM ordering (Amazon CVPRW 2026, A1 §8)
- Evaluation: Metrics Reloaded https://arxiv.org/abs/2206.01653 ; Nikolov NSD https://arxiv.org/abs/1809.04430 ;
  MedSAM https://github.com/bowang-lab/MedSAM ; LUNA16 box study https://arxiv.org/html/2608.16832
- Rejection-sampling / self-training SFT (D28, E10; abstracts read 2026-10-20): RFT https://arxiv.org/abs/2308.01825 ;
  ReST-EM https://arxiv.org/abs/2312.06585 ; Quagmires abstract (D5 correction) https://arxiv.org/abs/2510.01624
- SimpleSeg training settings (Tables 3–5, PDF read 2026-10-20): https://arxiv.org/abs/2601.19228
- Judge rounds: `experiments/phase2/JUDGE_R1.md`, `JUDGE_R2.md`; B1 amendments `lab/B1/prereg.py` (AMEND 1–5); B3 census
  `lab/B3/census.json`
- Constrained decoding: XGrammar https://arxiv.org/abs/2411.15100 ; GrammarRL https://arxiv.org/html/2609.39869
- Checkpoint selection (D12; B2 §9; abstracts re-read 2026-10-20): SFT rank inversion and p(1−p) advantage variance
  https://arxiv.org/abs/2606.18487 ; WiSE-FT for reasoning (early ↔ late interpolation recovers pass@k)
  https://arxiv.org/abs/2504.10478 ; Quagmires in SFT → RL https://arxiv.org/abs/2510.01624 ; DAPO
  https://arxiv.org/abs/2503.14476
- Process pools: https://docs.python.org/3/library/concurrent.futures.html ; https://docs.python.org/3/library/gc.html
- Qwen3.5 kernels and settings (D15, D13; B1 §7): flash-linear-attention https://github.com/fla-org/flash-linear-attention ;
  causal-conv1d https://github.com/Dao-AILab/causal-conv1d ; HF kernels https://github.com/huggingface/kernels ;
  transformers 5.20.1 `integrations/hub_kernels.py` and `models/qwen3_5/modeling_qwen3_5.py` (local source, read
  2026-10-20); ms-swift Qwen3.5 best practice https://swift.readthedocs.io/en/v4.2/BestPractices/Qwen3_5-Best-Practice.html
- PEFT LoRA `exclude_modules` (D6, D7): https://huggingface.co/docs/peft/main/en/package_reference/lora
