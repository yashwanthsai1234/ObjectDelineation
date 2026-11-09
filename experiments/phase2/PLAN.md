> Frozen record (written 2026-11-06). Versions, budgets and paths were true then. Current facts: docs/README.md.

# Phase 2 plan: SFT of a small polygon VLM on Vast AI

Part of [ROADMAP.md](../../ROADMAP.md). Results: [RESULTS.md](RESULTS.md). Build log: [V0_V1.md](V0_V1.md). Issues:
[impl/ISSUES.md](impl/ISSUES.md). Approved at revision 8 (judge round 4, `JUDGE_R4.md`; rounds 1–3 in `JUDGE_R1–R3.md`).
Research notes: `docs/research/notes/Phase 2 SFT/A1…A6` (cited as A1 §x), lab reports `lab/B1…B5` (cited as B1 §x).

## Revision 9 (2026-11-02): owner decisions after rung 1. These are binding. Where a section below disagrees, this box wins.

1. **Ladder.** Rung 1 (LFM2.5-VL-450M) was trained and failed three gates (RESULTS §1). It is not trained again.
   Qwen3.5-0.8B is skipped. The next rung is **Qwen3-VL-2B**, LM-only LoRA r16, one RTX 4090 (D1, D6, §4.8).
2. **CUDA checker now.** The checker is built in `cuda/` (ROADMAP Phase 8, first part). Lithography stays in Phase 8.
   Phase 2 does not wait for it.
3. **Preflight gate.** `make preflight` must pass before any GPU rental (§4.7 step 0).
4. **Data.** The scale jitter (D10) did not cover the rescaled exam view. The data are rebuilt: `lidc-bundle-v2` (D10, §6.0).
5. **No automatic fallback run.** §6.3 rule 1 (hole oversampling ×4 to ×1) targeted the wrong cause: the control failures
   were 4.3 % illegal, 1.3 % unreadable and 0.7 % false holes. Rules 2 and 3 had no evidence. All three are dropped (§6.3).
6. **G-precision.** Rung 1's verdict is zero-shot only. The capacity run was not run, so that verdict is not final (§5.1).
7. **Money.** Vast credit is $3.18. North America only. RIGA never on a rented host. GPU-hours and dollars (§4.9) go to the
   owner before anything is rented.

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
  passes**; larger rungs are not run. Rung 1 (LFM) failed (RESULTS §1); Qwen3-VL-2B is rung 2 (revision 9).
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
C2 on the chosen rung with 3 A0′ seeds and 1 A0 seed (R3–R5, revision 4), plus at most one C2-outcome
conditional run (×8 seed on a partial; §5.2). Inference-only baselines (image-blind prior, MedSAM, SAM 2.1-tiny) run on the Mac CPU. **Removed or
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
| D1 | **Model ladder** (revision 9) | **Qwen3-VL-2B is the next rung.** Rung 1, LFM2.5-VL-450M, was trained and failed (RESULTS §1); it is not trained again. Qwen3.5-0.8B is skipped (owner, 2026-11-02). **Stop at the first rung that passes.** SmolVLM-256M: control only, never C1 (B1: last on every metric, vision collapse under all-linear LoRA, TRL #6294). Qwen3-VL-4B or InternVL3.5-1B-HF: only after the 2B fails, owner approval. | LFM rerun on fixed data (owner: no repeat). Qwen3.5-0.8B (owner; its GDN layers need extra kernels). Side-by-side screening (B1 §0.5: costs more, answers nothing C1 asks). Florence-2, Ovis2.5, Moondream 3, Gemma, PaliGemma 2, MedGemma (A1 §6) | B1 §0, §2–§4; A1 §2, §5.1; RESULTS §1. Qwen3-VL-2B-Instruct: Apache-2.0, 2.13 B parameters (LM 1.72 B, vision 0.41 B), 28 layers, 8 KV heads, 152k vocabulary, 256 image tokens at 512 px (checked on the model files 2026-11-02) |
| D2 | Output format | **Keep P1-D10** for every model (compact JSON, ints 0–1000, CW outer / CCW holes, top-most-then-left-most start, `"holes":[]` always). **A0** = same grammar with `holes` emptied, parts kept ("A0-multi"); A0-single reported only as a ceiling. Qwen3-VL-2B-Instruct has no thinking mode. `SFTTrainer` swaps in TRL's `qwen3_vl_training` template; the assistant mask then covers the answer and `<|im_end|>` (checked on the real processor, 2026-11-02, `tests/test_p2_train.py`). Zero-shot, no model emits D10 (0 % parse for all three in B1 §2.1) and the native box priors are weak (B1 §2.1, H-A1-1 inconclusive), so the plan relies on no zero-shot prior; LFM replaces its [0, 1]-float convention with ours within 36 steps (parse 94–99 %, B1 §4.5) | Delta (−25 % tokens, but a ±1 delta-px slip costs 13 IoU pts on LIDC; A1 §1.2). Pixel units (−23 % LIDC tokens; 1 px per last-digit slip, −2.8 IoU; cannot express bevel vertices in 0.34–0.47 % of targets). Coarse grids (A1 §1.5, RIGA: G = 100 costs −2.7 IoU pts and −2.2 pts legality; G = 250 costs −1 IoU pt and −1.5 pts legality). Coordinate tokens by default (SFT-scale ties: GETok 92.9 vs 93.0, Shikra; init-sensitive, HiVG 0.44 vs 0.72): only via E8. Random start (−40 IoU, LPM T3) | A1 §0, §1, §3, §6; every target coordinate in [168, 832] (A1 §1.1) |
| D3 | SFT framework | **TRL 1.16 `SFTTrainer`** via `objdel.train.configs.sft_config`; assistant-only NLL with TRL's training-template swap (stop token trained for all models; a custom collator must pass `get_training_chat_template`) | Unsloth, LLaMA-Factory, ms-swift, Axolotl, NeMo, oumi (model-coverage gaps, second template/mask path); torchtune (unmaintained). Packing / padding-free (refused by TRL 1.16 for vision datasets; ≤ 0.2 % padding left after grouping) | A2 §1.1, §1.3, §2 |
| D4 | Throughput knobs | `train_sampling_strategy="group_by_length"` with a `length` column (completion tokens) from `sft_rows`; `dataloader_num_workers=6`, persistent, prefetch 2; gradient checkpointing **off** for LFM when peak ≤ 20 GB, on for Qwen3-VL-2B; SDPA attention; bf16 (`configs._bf16`); `PYTORCH_ALLOC_CONF=expandable_segments:True`. Measured in R0 (§6.2), not a separate A/B | Liger (deprecated path in TRL SFT, no `lfm2_vl`); FA3/FA4 (wrong hardware, sequences ≤ 2.7k); torch.compile (every batch a new shape); FSDP/DeepSpeed (single GPU) | A2 §1.2 (padding 7–32 % → ≤ 0.2 %), §1.4 (collate 110–350 ms per batch of 8, one core), §3; A6 §1a (10–25 ms CPU per row), §2 |
| D5 | Schedule and depth | **WSD**: 3 % warmup, constant LR, 15 % linear decay to 0; one stable run of **2 epochs** (each epoch a fresh augmented view, D10) with **decay branches at 0.5 / 1 / 2 epochs**, all kept; a third stable epoch only if the 2-epoch branch is the D12 pick and still improving on shifted-val IoU_tc by ≥ 0.01 over the 1-epoch branch (owner approval, §4.9). Compute per screening run ≈ **2.5 LIDC epoch-equivalents** (2 stable + 0.075 + 0.15 + 0.3 decay; §4.9 uses this). D12 picks the depth among these branches (B2 found no over-training signal at Mac depth, B2 §3). Peak LR 5e-5 full FT, 5e-4 LoRA (10×); no LR sweep. Effective batch 32 sequences. `seed` and `data_seed` set. **Why fresh views per epoch:** the crop leak (A4 S5–S6, D10), not data uniqueness. Quagmires (2510.01624, abstract re-read 2026-10-20): "SFT training on unique examples for one epoch underperforms training on half examples for two epochs, either after SFT or SFT-then-RL", so 1–2 epochs is the evidence-backed range and the 2-epoch cap loses little | Cosine single-depth run (one depth only). 3 stable epochs + 3-epoch branch by default (revision 2; +1.5 epoch-equivalents per run for a depth the evidence does not favour). Easy-first / short-only curriculum (Quagmires: short-only "often leads to worse outcome after RL") | A2 §4.3 (WSD); SimpleSeg (2601.19228, Tables 4–5, re-read 2026-10-20): **SFT** peak 5e-5 cosine to 2e-6, 3 % warmup; RL constant 2e-6 (GSPO, Muon, 7B), §4.7; Quagmires abstract (held-out loss and pass@large-k predict RL outcome better than SFT score); P1 `configs.py` LRs |
| D6 | Full FT vs LoRA (revision 9) | LFM2.5-VL-450M and SmolVLM: **full FT**. **Qwen3-VL-2B: LoRA r16 on the language model's linear layers; vision tower and projector frozen** (PEFT `exclude_modules`); base in bf16, adapters in fp32; merged after SFT; LR 5e-4. **Memory (2.13 B parameters).** Full FT with AdamW needs 16 bytes per parameter (fp32 weights 4, gradients 4, two moments 8) = about 34 GB before activations: it does not fit a 24 GB 4090 or a 32 GB 5090, and needs a 48 GB card (L40S / A6000 class, few US offers, about twice the price). LoRA needs 4.3 GB of weights, 0.3 GB of adapter state (17.4 M trainable parameters, 16 bytes each), and fits the 4090 with micro-batch 4 and checkpointing. The FLOP cost per step is the same as LFM's full FT (6·N·T with recompute). **Risk:** a frozen tower limits localisation (A2 §4.2). If rung 2 fails G-geometry with legal ≥ 90 %, the first change to consider is training the tower (owner decision, §8.3). | Full FT on 48 GB (about 2 x $/h, no evidence it is needed yet). LoRA on the tower too (B1: collapsed image features at 5e-4). DoRA / rsLoRA | A2 §4.1; A3 §3.4; B1 §2.4; memory arithmetic in `configs.py` and §4.9 |
| D7 | Vision tower | **LoRA (Qwen3-VL-2B): vision tower and projector frozen** (PEFT `exclude_modules`, B1's `common.lora("lm")` regex). **Full FT (LFM, SmolVLM): vision tower and projector at 0.1× the LM learning rate** (5e-6), via optimizer param groups. Both checked by the **G-vision** gate on every SFT branch (§5.1) and by an in-training projector check (M8): if at a rolling checkpoint **both** G-vision (a) ratio < 0.5 **and** (b) TV < 0.002 (B1's collapse hit both; (b) costs 8 items × 24 tokens on the GPU), the run is aborted and restarted once with tower and projector frozen (a pre-registered repair, not the §6.3 diagnostic); one of the two alone is warn-only (logged) and the decision is left to `select_ckpt`'s G-vision (revision 7). **Why:** the all-linear adapter at lr 5e-4 (`configs.lora`, vision included) erased image information for Qwen3.5 and SmolVLM within 36 steps (projector cosine distance between images 0.297 → 0.019 and 0.171 → 0.003; next-token TV real vs black or another item ≤ 0.0007 vs 0.02–0.08 at base); LFM kept its vision (0.471 → 0.321–0.373); LM-only LoRA kept Qwen3.5's features with unchanged legality (82 %) and centred IoU (0.20 vs 0.19) (B1 §0.2, §2.4–2.5). Frozen, not 0.1×, for LoRA because frozen is the setting B1 tested; 0.1× on vision LoRA weights is untested. Full FT keeps 0.1× because A2 §4.2's evidence supports unfreezing for localisation (VectorLLM, Visual Position Prompt, VLM4VLA) and a 5e-6 full-FT update is a different regime from B1's 5e-4 LoRA one, so it is checked rather than assumed (B1 §4.2; H-A2-6 now P1). Phase 3 note: vLLM encoder cache must be reset after each weight sync when the tower trains (TRL #7198) | All-linear LoRA incl. vision at the LM LR (revision 4's Qwen3-VL-2B setting; collapses vision, B1 §2.5). 0.1× vision LoRA LR (untested; B1 offers frozen as the tested alternative). Freezing for full FT by default (discards A2 §4.2's localisation evidence; G-vision and the repair cover the risk). A frozen / unfrozen A/B (owner directive; E4 trigger in §6.4). Qwen recipe LR 2e-7–1e-6 frozen vision (tuned for mixed image/video data, not a new format) | B1 §0.2, §2.4–2.5, §4.2; A2 §4.2; A5 §3, H-A5-9 |
| D8 | Loss | **Assistant-only NLL. Number-token loss NOT in the bundle.** The evidence is too weak to adopt without a test: NTL's gains are on math and text regression (math interpolation 64 → 75 %, rJokes RMSE 2.01 → 1.81; A1 §3, A2 §4.5), with no VLM-coordinate follow-up found (A1 §3; A2 §10); the closest polygon evidence, Gaussian label smoothing on coordinate tokens, left IoU flat (83.77 → 83.67) and hurt on WHU (LPM, A1 §3); NTL is designed for single-digit tokenizers while the first rung, LFM, uses one token per 0–999 value (A1 §1.5). Its expected benefit (boundary error) is exactly what G-precision and legality measure, so it is the triggered option if those fail (§6.4; M6 built only then) | Adopting NTL by default (unverified on VLM coordinates; new loss code on the critical path). Label smoothing (deprecated in TRL SFT, extra forward pass, does not repair over-trained checkpoints) **although A4 §4.3 recommends 0.1**. DFT (mode-seeking, entropy collapse; H-A2-8, untested: B2 ran no training, B2 §0). GEM, CHORD (parked) | A1 §3; A2 §1.5 (prototype matches TRL's fused NLL), §4.5, §10 |
| D9 | Targets | **All voting raters as rows** (P1 I27; no consensus/STAPLE). **Per-variant dihedral targets built offline** from transformed masks with **crop-before-trace** (B3 P4: hash `1af819f4…` unchanged, RIGA −23 % per target). Every view has its own `target_set_hash`, the K5/K6/R1 re-check and the per-variant token-budget filter (I5) | Naive vertex transform + canonicalize (breaks the IoU ≥ 0.99 certificate on 7.3 % of RIGA cups; A4 S1). On-the-fly target rebuild (26 ms per RIGA target in the collator; A6 §1b). One random rater per epoch (same expected gradient as all raters, more bookkeeping; A2 §4.6) | B3 P3 (215,208 variant targets, 0 errors, two runs byte-identical, 16.5–18.1 CPU-min); A4 S1; A6 §4 |
| D10 | Augmentation | **Revision 9 (2026-11-02): the scale jitter did not cover the rescaled gate view.** Rung 1 trained with S drawn from a margin rule; train nodules under 10 % of the crop were 0.5 % of rows, the rescaled exam had 27 %. The new bundle `data/views/lidc-bundle-v2` draws S uniformly from the policy side and its ladder neighbours (LIDC 50 / 100 / 125), which are the sides of the exam views (`build_views.py`, `targets/augment.py`; coverage check `eval/coverage.py`, `tests/test_eval_coverage.py`). **PLACEHOLDER: the final rule, the coverage numbers per exam view and the bundle hashes are copied here when the data work is final; until then `configs/r1_qwen3vl2b.yaml` names the view `lidc-bundle-v2`.** The text below is revision 8. **Crop shift + scale jitter is mandatory** in every run (the centred frame leaks location and size: the image-blind prior passes RIGA's gates on centred frames and collapses under shift, A4 S5–S6; and every B1 arm trained on centred crops drew at the frame centre: on shifted crops 94–98 % of parsed samples lie closer to the frame centre than to the nodule, median IoU_tc 0, every arm significantly below the image-blind prior, B1 §0.1, §2.5). Shift: uniform, union inside with ≥ 2 px clearance, and no other outlined nodule (`meta.other_nodules`) inside the window (B3-2 rule kept). Scale: margin drawn from [8 px, 60 %] so S moves between divisors (LIDC 50/100/125, RIGA 200/250/500), S ≤ 500, scale 1. **D4 for LIDC; horizontal flip only for RIGA** (fundus orientation is physiological). One op per rater row: op = hash(seed, epoch, item, row). **CT window jitter** L −300 ± 60, W 1400 ± 80 HU on raw HU, noise σ ≤ 10 HU, applied to **int16 HU context crops written from `data/store`** (M1; `data/prepared` images are already lung-windowed uint8, so jitter cannot be applied to them); random CT windowing follows 2510.08116 (liver-tumour CT, a single study), so this component is the weakest-evidenced in the bundle; RIGA brightness/contrast 0.75–1.25, gamma 0.7–1.5, p 0.15–0.3, no hue. Images are cut and transformed at frame resolution **before** the 512-px Pillow upscale | Elastic, non-90° rotation, mixup, cutout, copy-paste, random start vertex (masks resample and need recertification; copy-paste changes the ground truth; A4 §4.2). RIGA vertical flips / rotations (non-physiological; AADG). Transforming the upscaled LIDC image (rotations differ by up to 16 grey levels; A4 S1) | A4 §0.1, §4.2, S1, S5–S6; B1 §2.5 centre bias; ROADMAP's "flips, 90° rotations" for both datasets is superseded for RIGA |
| D11 | Data mixing | All Phase 2 Vast runs are **LIDC-only** (D16; LIDC carries the topology claim). A RIGA joint run, if RIGA is permitted, uses natural proportions by targets (A2 §4.8: RIGA 12 % of rows, ~30 % of answer tokens); the T = 2 mixing A/B (E5) is removed (trigger §6.4). **Non-simple LIDC train targets ×4 in every run** (bundle, §6.0); A0 runs use the same repeated rows with holes emptied | Separate per-dataset models (A2 §4.8). A ×1/×4/×8 oversampling A/B (E1, removed) | A2 §4.8, §7.4; A1 §5.2; P1 RESULTS §7 #9 |
| D12 | Checkpoint selection (B2 filled, revision 6: rule R-B) | **Eligible:** WSD branches (0.5 / 1 / 2 epochs; R1 Qwen3-VL-2B also 3, owner 2026-11-05) passing G-format (parse ≥ 98 %, legal ≥ 90 %) and G-signal (duplicates ≤ 10 %, σ-share ≥ 50 %) at **T = 1**, K = 8 on the shifted val, and G-vision (§5.1). **Select the deepest eligible branch.** **Signal veto:** take the next shallower eligible branch instead only if both hold: (i) its geometric reward spread **sd_geom = √(p·Var_legal(R))** on the shifted val (R = the G-probe raw composite, illegal = 0; per item over its K = 8 group, then averaged) is higher, with the paired nodule-cluster CI excluding 0 **and** by ≥ 15 % relative; (ii) its shifted-val median IoU_tc is not significantly lower (paired CI). Applied pairwise downward from the deepest branch; ties stay deep. **Veto statistic: sd_geom, a post-hoc statistic, adopted by this decision and labelled so** (B2 §1 "post-hoc", §4 note, §7). Why not r_std (B2's pre-registered statistic): Var(R) = p(1−p)·μ_L² + p·Var_L(R) exactly, and the legality part rewards a branch for being *less* legal; it carried 55–92 % of the group variance in B2 and the whole 36 → 150-step fall of r_std (B2 §2.1). Inside the eligible set p(1−p) still moves from 0.09 to 0.02 between legal 90 % and 98 %, enough by itself to clear a 15 % margin. **Checked by R2's folded measurement** (§6.2): per-item Spearman of step-0 sd_geom, and of step-0 r_std, with the G-probe gain. That per-item correlation is confounded by legality headroom (items near p ≈ 0.5 have the most room to gain, which mechanically favours r_std), so it can only **retain** sd_geom or **flag** disagreement; switching to r_std needs the branch-level comparison (D12(a) / (b), when a second branch is probed) and is a Phase 3 planning decision, not an automatic rule (revision 7). Both statistics are reported per branch with a "statistics disagree" flag when the veto would fire under one only; the choice is written into `prereg.py` (M8) before R1. **Never select by:** greedy legal or IoU; σ-share, active-group share or distinct share (saturated at 98–100 %, B2 §2.1); total r_std or headroom (fall mechanically as legality rises). **Reported per branch, not used:** r_std, headroom, mixed-legality share, legal-but-wrong share, mean \|A_v4\|, TRL entropy, held-out NLL, pass@8 (B2 §4). The selected branch then runs **G-probe** (§5.1); a G-probe result never climbs the ladder on its own. **D12(c) repair** (selected branch not PASS on G-probe): the **next-earlier eligible branch**, probed once. Base ↔ LoRA-delta interpolation is **not** a repair (B2 H-B2-2 REJECTED: α = 0.75 loses 7–9 legal pts on both views, α = 0.5 reaches 38 % legal and 10.5 % duplicates; B2 §2.2). WiSE-FT between two SFT branches (0.5-epoch ↔ selected, A2's H-A2-7 proper, untested: B2 could only run the base ↔ late form) is scored beside the repair branch when the repair fires, at inference cost only (merge at α = 0.5, shifted val K = 8, ≈ 0.05 GPU-h), report-only, for Phase 3. **D12(a) / (b)** (does the signal branch out-gain the greedy one under GRPO?) are not testable on the Mac (B2 §3, §6): when a second branch is probed, both paired gains are reported side by side. H-A2-3 is INCONCLUSIVE at Mac depth (B2 H-B2-1: LFM's r_std and headroom fall beyond seed noise, all in the legality part; the geometric part is flat or rises), so B2 shows no over-training at ≤ 150 LoRA steps; the veto exists for the 2-epoch full-FT regime where 2606.18487's rank inversion may appear | **R-A**, revisions 1–5's signal-first lexicographic rule (σ-share, duplicates, pass@8, NLL, IoU): its first keys are saturated (σ-share and distinct 98–100 %, duplicates ≤ 1.4 % at every depth), so it would decide on noise, and B2's pre-registered rule keeps R-A only if H-B2-1 is ACCEPTED (B2 §1, §3). **r_std as the veto statistic** (pre-registered, but it penalises legality gains; above). Total r_std or headroom as selection keys. Base ↔ LoRA interpolation as the repair (H-B2-2 REJECT). Selecting by SFT score (rank inversion: pre-RL pass@1 up, peak GRPO pass@10 0.806 → 0.481; 2606.18487) | B2 §0–§7 (`lab/B2/REPORT.md`); A2 §4.10, H-A2-3 / H-A2-7; 2606.18487 (advantage variance p(1−p)(g−1)/g, rank inversion); 2504.10478 (early ↔ late WiSE-FT recovers pass@k); 2510.01624; RESULTS §8 GRPO-signal gate |
| D13 | Eval sampling (revision 9) | **HF `generate` is the sampler for every config now** (`hf_sample`, `sampler: hf`). R0 parity (2026-11-02) failed for LFM (p99 per sequence 1.11 nats, limit 0.5). The vLLM path (`sample_vllm`, `/opt/vllm` 0.32.0, JSONL bridge) stays in the image; R0 repeats the parity for Qwen3-VL-2B, and a pass lets R1 use it. **Revision 10 (owner, 2026-11-05, before any R1 result): R1 Qwen3-VL-2B samples with vLLM.** R0 v3 (2026-11-03): prompt ids 100 %, token median \|Δ\| 4.65e-6 (limit 0.02), token p99 0.109, sequence p99 1.79 (limit 0.5). The sequence p99 becomes report-only: bf16 rounding adds up over ~580-token answers. The question it stood for, whether vLLM's answers score like HF's, is measured directly: R0's 256 base-weight answers of each sampler scored the same (all 0 for format, legal, IoU; tokens 581 vs 568, difference CI [-69, 101]), which is not informative for a trained model, so `select_ckpt.hf_check` re-samples the judged branch with HF on 64 seeded shifted-val items x K 8 and gate **G-hf-check** fails when any vLLM - HF metric difference (format_ok, legal, iou_tc, f5, g_area; per-item means, nodule-cluster bootstrap) has 0 outside its 95 % CI. HF projection 18.9 h > stop 5 h; vLLM (0.09 s per answer, 95 s engine start, measured in R0 v3) plus the check: see `project_r1`. HF batches are 8 items x 8 samples = 64 sequences, left-padded, bf16 (Qwen3-VL-2B's KV cache is 112 KB per token: 13 GB at 1,856 tokens; 16 items would need 27 GB and one runaway answer would overflow 24 GB). Every sampling field is explicit and the same for HF, vLLM and G-probe's GRPO: T 1, top_p 1, top_k 0, min_p 0, repetition 1, presence 0, frequency 0, `max_tokens = MAX_COMPLETION[model]` (1,024 LFM; 1,536 Qwen3-VL-2B and SmolVLM), stop ids from the checkpoint. They override the cards (Qwen3-VL T 0.7 / top_p 0.8 / top_k 20; LFM T 0.1 / min_p 0.15 / repetition 1.05). **GPU sequencing:** the trainer never samples. `select_ckpt` starts after `train_sft` exits, does not touch CUDA itself, and launches the sampler as a child; G-probe runs later in its own process. Branch selection uses a fixed 256-item shifted-val subset x K = 8; the full 1,084 items only for the selected branch. **Gates are unconstrained** (a constrained decoder hides format failures). K: 8 for format / geometry / signal, 16 for HM-IoU16 (report-only), 64 later. | TRL + vLLM colocate for eval (breaks SmolVLM: TRL #6294). Constrained gates. Continuous batching (no image path). SGLang, TensorRT-LLM, LMDeploy, MLC | A5 §0–§8; B5 Q2; R0 parity 2026-11-02 (V0_V1.md) |
| D14 | Reward serving | **Preload the RL caches of the datasets that process scores, and only in processes that score** (`select_ckpt`, G-probe, `k3_harness`, Phase 3 GRPO), never in `train_sft` (SFT needs no caches beyond the voting indices stored in the records, M2). LIDC-only: 0.33 GB array bytes (B4 §1); RIGA adds 3.35 GB (Linux RSS ~4.2–4.5 GB) only on RIGA runs. Load into a dict **before model/CUDA init**; `gc.freeze()`; `ProcessPoolExecutor(4, mp_context=fork, initializer=cv2.setNumThreads(0))` warmed with 4 no-op tasks so all workers fork immediately; jobs carry `CacheRef(key)` (pickles 211 KB per RIGA batch vs 96.2 MB with caches); `TRLReward(caches=serve.loader(root), pool=pool)` unchanged. A dead worker raises `BrokenProcessPool` → fail loudly, resume from checkpoint. **C++ batch scorer deferred**; build it (B4 §4 design, ~600–900 lines + differential tests) if, on a quiet Vast box, RIGA-only pool-4 paired wall > 125 ms or speed-up < 3×, or Phase 3 needs larger groups, or < 4 cores are spare beside the dataloader | LRU-64 `cache_loader` (32/32 misses per RIGA step: 4.79 s vs 0.65 s preloaded). Threads (0.8–1.2×). `multiprocessing.Pool` (hangs forever on a killed worker). Forkserver/spawn with a TRL script as `__main__` (re-imports torch: macOS crash 3/3, Linux pool start 49.7 s). GPU scoring (exactness, launch cost) | B4 §1–5; A6 §0.4–0.5, §4 |
| D15 | Environments (revision 9) | **One image, two venvs**: `/opt/venv` (repo pins: torch 2.14.1 cu130, transformers 5.20.1, TRL 1.16.0, PEFT 0.21.2, numpy 2.5.4, numba 0.68.0) and `/opt/vllm` (`vllm==0.32.0`: torch 2.13.0, transformers 5.19.0). Qwen3-VL-2B runs in `/opt/venv` with no extra kernels (class `qwen3_vl`, SDPA). The image still holds the Qwen3.5 extras (`flash-linear-attention`, `causal-conv1d`, `kernels`); no run needs them now. **Phase 3: one re-pinned colocate venv** (TRL 1.16, vLLM 0.32.0, transformers 5.19, torch 2.13.0) with uv overrides numba 0.68.0 and numpy 2.5.4; it passed imports, the fast suite and `xarch_hashes --strict` 26/26. Re-run both on every re-pin. | vLLM in the pinned `.venv` (pins conflict). vLLM ≤ 0.31 (caps transformers < 5.18; vision assistant-only loss needs ≥ 5.18). Server mode on one GPU (TRL 1.16's client is NCCL-only). SGLang | B5 Q1, Q3; A5 §0.2; A3 §6.1 |
| D16 | Data governance | **LIDC (CC BY 3.0) on any verified host. RIGA on a rented host only with written permission (depositor + MESSIDOR partners) AND only on Secure Cloud** (`datacenter=True`, ISO 27001): encrypted `bundle-riga.tar.age` decrypted into `/dev/shm`; **key delivery:** the age identity lives only in the Mac keychain and is streamed over the direct SSH session's stdin into the decrypt command (`age -d -i /dev/stdin bundle-riga.tar.age \| tar -x -C /dev/shm/riga`), never as an env var, template field, onstart argument or file on the host; the ledger records the recipient (public) key only; exact invocation tested in M7; no Cloud Sync, destroy on completion, ledger row per session; `sync_up` refuses it to any offer without `datacenter=True` (H-GOV). Until permission: RIGA gates are not judged (`phase2_entry(..., datasets=("lidc",))`); Mac MPS remains the last-resort RIGA path. **RIGA-free K3 harness** (`k3_harness`, synthetic RIGA-like items) on community hosts. Never publish RIGA-derived data; Phase 9 release covers LIDC only. Register at ADCIS for MESSIDOR | RIGA on community hosts even encrypted (host root, Vast backups outlive instances; A3 V12, V14). Third-party RIGA mirrors. Option D split (only 290/750 images non-MESSIDOR) | RESULTS I44 (permission AND Secure Cloud); A3 §5; A4 §5; B4 §3 (synthetic RIGA: CPU-time ratio 0.972 [0.935, 1.016]; digest `e31d8f1e…` on arm64 and x86) |
| D17 | Baselines | Same split, same `reward.score_mask` path, and **every row labelled with its view; the gate-relevant rows use the same shifted val as the VLM gate** (and the rescaled view where reported): (1)–(3) are inference-only and run on the **Mac CPU** before V1, and also compute the **per-item centred win share over the prior for MedSAM-oracle and SAM 2.1-tiny** (K = 1 each; the reference for G-prior (a), §5.1) (A4 §3.3: ~3.5 h CPU for all frames; $0). (1) **image-blind prior** (disc at frame centre, train-median diameter per S; mandatory, gating) on centred, shifted and rescaled views; (2) **MedSAM ViT-B re-run**, official min-max preprocessing, on the centred **and shifted** views, box from that view's union (+1 px LIDC, +5 px RIGA: still an oracle box, now on the matching view) plus a jittered-box row (0–20 px at 1024 scale); (3) **SAM 2.1 hiera-tiny**, box, same views (Apache-2.0; H-A4-10); (4) **nnU-Net v2 2D** (training removed from Phase 2, Phase 9 baseline, §6.4), label-sampling and majority-vote, 1 fold, 250 epochs, **trained on the jittered train view (D10 shifts + scale jitter) and evaluated on the shifted val**, never on centred frames (it would inherit the S5 leak); (5) MedSAM2 optional, research-use weights only. Published Dice converted with IoU = D/(2 − D); our LIDC numbers never share a column with Kohl v1/v2 | Phase 0 MedSAM numbers (LIDC 0.33 was a crop/scoring artefact; re-run gives 0.78). SAM 3 / Medical SAM3 (11.9 % external Dice text-only; 3.45 GB), BiomedParse (trained on LIDC: contaminated; NC), LISA/GLaMM/MedPLIB (LIDC Dice ≤ 0.006), MedSAM-Agent (inconsistent) | A4 §0.4, §3, S7 (MedSAM oracle box: LIDC 0.780, disc 0.941, cup 0.847; box-ellipse null 0.602) |
| D18 | Entry gates | Revised in §5.1: geometry on **shifted** val and on a **rescaled** val (shift + S moved to a neighbouring divisor); **paired win over the image-blind prior** on the shifted and rescaled views (the centred win share, G-prior (a), is report-only since revision 7); precision probe re-specified (D19); unanimous-topology gate **replaced** (not added to) in `eval/gates.py` by the synthetic topology probe (D20) with a hole-recall floor; **G-probe**: a 20-step GRPO run must improve the checkpoint; "IoU near MedSAM" reported, not gated (MedSAM gets an oracle box the policy does not; A4 §3.2) | RESULTS §8 gates as written (geometry passable image-blind on centred RIGA frames; topology gate vacuous: 0 unanimous non-simple val items; ≤ 2 px probe looser than the readers) | A4 §2.5, S5–S6; B3 §1–2 |
| D19 | Precision probe | Synthetic shapes, **generator 2** (elongation E = 0.2, amplitude A = 0.15; calibrated: centroid-disc IoU 0.822 vs real 0.806; σ 0.5 → IoU 0.915 / legal 96.4 %; σ 1 → 0.816 / 85.4 %), LIDC diameter distribution, uncentred frames (window + scale jitter), real LIDC lung patches. Protocol per rung (revision 3): **zero-shot** first, the screening checkpoint on 1,000 held-out probe items (inference only); the capacity run (probe-only SFT, 20k items, one-LIDC-epoch budget, ~0.3 GPU-h on a 4090, then the same 1,000 items) only if zero-shot fails; **patient-disjoint backgrounds** (revision 9, 2026-11-02): probe-test patches (the 1,000) come from LIDC **val** patients, never seen in SFT training (40 patches, 19 patients, the topology probe's pool), probe-train patches (the 20k capacity set) from LIDC train patients; the R1 probe used train-patient patches split by a patient hash (43 patients also in SFT training, found by preflight), archived in `data/eval_views_archive/probe-precision-r1`; the topology probe (D20) uses val patients only; greedy + K = 8 at T = 1. Pass on 6–16 px objects: legal ≥ 95 %, median ASSD ≤ 0.6 px, median IoU ≥ 0.80, median centroid error ≤ 0.5 px, **paired per-item IoU win over the equal-area centroid disc ≥ 60 %**. Topology moves to D20 | RESULTS §8 "≤ 2 px boundary error" (readers' ASSD 0.93 px; σ = 2 px gives LIDC IoU 0.63 and 50.9 % legal; A4 S2–S4). Probe on real LIDC (rater noise dominates). Generator 1 (failed calibration: disc IoU 0.771, σ = 1 legal 69 %) | A4 §6; B3 P2 (exact, brute-force-equal, 1-px shift detected 100 %, discrimination PASS; disc passes every geometric criterion, hence the paired-win rule: σ ≤ 0.5 wins 78–100 %, σ ≥ 0.75 / shift / disc ≤ 49 %) |
| D20 | Topology probe (B3 Addendum A) | **Synthetic**, independent items on real LIDC **val** lung patches: **120 one-hole items per hole-area bin (≤ 10 / 11–30 / > 30 px; 360 total), 160 two-part items (descriptive only: A0 keeps parts), 160 simple controls**. Outlines generator 2 (E = 0.2, A = 0.15). Each hole jointly resamples one real LIDC **train** hole record (holed part's diameter, hole area, centroid offset, wall to the outside); each two-part item resamples one train multi-part record (diameter, part-area ratio, gap); every hole ≥ the size rule (signature (1,) raw and at `min_area(\|GT\|)`; the rule max(5 px, ⌈0.5 %·\|U\|⌉) lives in the reward code, `objdel.checker` `min_area`, P1 PLAN §4.2, not in a note). **Real hole geometry** (train, voting raters, size rule; `hole_census.json`): 285 holes in 241 holed targets; area p10 / p50 / p90 = 5 / **12** / 49 px; bins ≤ 10 / 11–30 / > 30 px = **44 / 36 / 20 %**; hole = 3.6 % of its part (p50); wall p50 1.24 px (p10 0.41); holed part ≈ 24 px across; **75 % of holed targets come from items where only one rater drew a hole; no hole is unanimous**. **Primary metric (C2): per-item hole hit averaged over K = 8 samples at T = 1 and over the A0′ seeds (3 since revision 4; B3 A5 specified 2); pooled hole recall = train-share-weighted mean of the bins (0.439 / 0.361 / 0.200), stratified bootstrap of items within bins** (±3.1 pts at recall 0.5, item difficulty σ_u = 1.5; greedy ±4.9, any-of-8 ±4.4 and saturating at 0.99 for σ = 0.5). Greedy and any-of-8 reported, never deciding. **Generator fidelity (P5.2): failed narrowly** as pre-registered: the final generator (v3e) passed KS p ≥ 0.01 and bin shares ±5 pts on fresh seed 8 but not on fresh seed 7 (lowest p 0.004, hole area); the residual is a ~1 px area under-draw of exact 5-px holes in the smallest bin (probe bin medians 6 / 17 / 48 px vs real 5 / 17 / 49) and a slightly less central low-offset tail (p10 0.25 vs 0.18); wall, fraction, part ratio and gap match on every seed (p ≥ 0.11). Iteration stopped to avoid tuning to the test. **Why it remains fit:** the probe is stratified by area bin, so the bias acts only within the ≤ 10 px bin (making those holes ~1 px larger, i.e. slightly easier), C2 also requires every bin's point estimate ≥ 0.30, and the smallest-bin number is reported with this caveat. Every train non-simple item stays in SFT. The 44-item val non-simple slice is reported descriptively (its holes are smaller, 54 % ≤ 10 px, mostly single-rater); its matching quantity is the K-sample hole rate vs the rater hole fraction (TV), not recall | ~300 held-out real items (H-A1-7: no real design reaches ±6 pts at ICC 0.15–0.19; removes 55–100 % of train non-simple targets). B3's first probe (hole geometry uncalibrated; JUDGE_R1 MUST-FIX 1). Any-of-8 as primary (saturates, rewards random hole stamping). Paired A0′ − A0 test (equals A0′'s recall; adds nothing) | B3 P1 (census, cluster CIs); B3 Addendum A1 (real geometry), A2 (fidelity), A3 (power, `power_c2.json`), A4 (sanity: exact 1.0, A0 0 / 0 / 0, centroid disc 0, σ = 0.5 px mean-over-8 0.65 / 0.81 / 0.90, pooled 0.76 [0.73, 0.78], 0 false non-simple), A5 (metric, seeds) |
| D21 | C++ vs Python | Per module, §3.4. Python stays everywhere in Phase 2 except polycore (P1-D1/D14); numba stays for the simplifier; the only C++ candidate is the D14 batch scorer behind its trigger; an exact FH EDT in C++ only if caches are rebuilt online (H-8) | GPU raster/EDT/CCL (exactness, launch overhead, GPU busy); free-threaded 3.14t (OpenCV abi3 only); piecemeal numba/C++ ports (1.6× serial, no thread scaling) | A6 §0, §4, §6; B4 §4 |
| D22 | GPU (revision 9) | **One RTX 4090 (24 GB), US / CA hosts**, about $0.42 per hour (§4.8). Fallbacks: micro-batch 2 x 16, then an RTX 5090. A 48 GB card only for full FT (D6). RTX 3090 for dry runs. | PRO 6000 96 GB and A100 80 GB (2–3 x the price; no need for LoRA); H100 and larger (wall clock does not bind) | §4.8; V0_V1 |
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
          ▼             (with sampler: hf, ③ is hf_sample in its own /opt/venv process instead)
 ④ /opt/venv: g_probe (selected branch; the D12(c) repair branch only on failure): preload caches → pool → GRPOTrainer 20 steps, use_vllm=False,
      eval_on_start + eval at step 20 on the fixed G-probe item set → G-probe JSON → EXITS    [GPU: GRPO only]

 CPU budget (≥ 12 effective cores): ① 1 trainer + 6 dataloader workers; ②/④ 1 + 4 reward workers; ③ vLLM
```

### 3.3 Phase 2 run graph (revision 9)

```
 DONE  V0 watchdog dry run (3090) ─► V1 day 1 and rerun (4090): R0, R1 rung 1 LFM ─► gates FAIL (RESULTS §1)

 Mac, no GPU:  make preflight passes  ─────────────────────────────────────────────►  owner sees GPU-h and $  ─► rent
        ▼
 V2  one 4090 (US / CA), LIDC only; processes chained one at a time (§3.2); Mac pulls every 15–30 min
   R0  image_check ─► xarch_hashes ─► k3_harness idle ─► parity (Qwen3-VL-2B) ─► hf_rate ─► G-vision base reference
       ─► 30-step arm (peak memory, s/step)        │ blocking failure ─► stop and fix
       ▼ Mac: project_r1 (stop rule 6.0 h); sampler: vllm (D13 revision 10)
   R1  Qwen3-VL-2B, lidc-bundle-v3, A0′, seed 1, WSD 3 epochs + branches 0.5/1/2/3 ─► train_sft EXITS
       ─► select_ckpt (D12 R-B; selected branch: rescaled val, probes, K = 16) ─► §5.1 gates
        │ gate fail ─► Mac diagnosis; owner chooses the next run (§6.3, §8.2)
        ▼ gates pass
   R2  G-probe on the selected checkpoint ── not PASS ─► D12(c) repair branch; still not PASS ─► C1 with a warning
        ▼ pass ══► C1 = Qwen3-VL-2B.
   R3–R5 C2: A0′ seeds 2 and 3, A0 seed 1 ─► topology probe ─► §5.2 3-seed rule
        │ partial ─► R-x8 (owner approval)
        ▼
   teardown (final pull, destroy, ledger)
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
  train/configs.py      + qwen3-vl-2b (LoRA, bf16 base, MAX_COMPLETION 1,536); explicit sampling
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
grow; V6): 80 GB is enough for Qwen3-VL-2B with LoRA (base 4.3 GB; four merged branches and three stable copies about 30 GB).

```
vastai search offers 'num_gpus=1 gpu_name=<GPU> gpu_ram<=<stock+2> verified=true rentable=true reliability>0.98
  cpu_cores_effective>=12 cpu_ram>=32 cuda_vers>=13.0 disk_space>=80 inet_down>=500 direct_port_count>=2
  duration>=1 dph<=<ceiling>' -o 'dph_total' --storage 80 --raw          # + datacenter=True for the RIGA bundle
```

Deviations from A3 §3.3 (recorded in §9): `cpu_cores_effective ≥ 12` (not 8), because a training box runs 1 trainer
+ 6 dataloader workers + 4 reward workers (D4, D14); snapshot medians are 16–34 vCPUs for all five cards, so supply
holds. `cpu_ram ≥ 32` for LIDC-only boxes (caches 0.33 GB, preloaded only by scoring processes, D14); `cpu_ram ≥ 48` for
RIGA boxes (RIGA caches ~4.2–4.5 GB Linux RSS, B4 §1). `gpu_ram ≤ 26` for a stock 4090 (16 of 104 "4090" offers report 48 GB; A3 §3.2).
Hosts in the US or Canada only (`policy.reasons()`, a fixed rule). Price ceiling: see `vast/policy.py`.

### 4.3 The image

`Dockerfile.phase2` (promoted from `lab/B5/Dockerfile.phase2` + `lab/A3/Dockerfile.vast`): base
`python:3.12-slim-trixie` at the bench image's digest; no CUDA base image (the torch wheel carries CUDA 13.0 user-space;
the host injects the driver). `/opt/venv` = the Mac's `uv pip freeze` (same file as `x86.sh`), `/opt/vllm` = the
`vllm==0.32.0` lock; openssh-server, rsync, tmux, procps, build-essential, cmake, ninja. `/opt/venv` adds image-only extras that Qwen3.5 needed (`flash-linear-attention`, `causal-conv1d`, `kernels`). Qwen3-VL-2B needs none of them: the pinned transformers 5.20.1 has `qwen3_vl`, so **the image is not rebuilt for this rung**. Removing the extras would cut the pull (owner decision, §8.3). Code is never baked in.
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
- Checkpoints: pulled by the Mac every 15–30 min. A Qwen3-VL-2B merged checkpoint is 4.3 GB (bf16); a LoRA rolling checkpoint is small. `pull.sh` skips `rolling/ decay/ stable/`. Resume = rsync the checkpoint to a new box +
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

### 4.7 The rung-2 session (V0 and the LFM session are done)

One RTX 4090, LIDC only, 80 GB disk. **Step 0, before any rental: `make preflight` passes** (the preflight suite is
`experiments/phase2/preflight.py`; owner decision 2026-11-02). Then R0 → R1 → (R2 → R3–R5 if gates pass).

| # | Step | Pass / record |
|---|---|---|
| 0 | `make preflight`; the GPU-hours and dollars of this plan are shown to the owner | no rental before both |
| 1 | `offers.sh RTX_4090` → `launch.sh` (US / CA only); record `lscpu`, `nproc`, cgroup `cpu.max`, `nvidia-smi -q` | ledger row |
| 2 | `sync_up.sh` LIDC bundle (`lidc-bundle-v2`); SHA verify 100 % | — |
| 3 | `run.sh … configs/r0.yaml`: `image_check.sh`, `xarch_hashes`, `k3_harness` idle, `r0_parity` and `hf_rate` for Qwen3-VL-2B, G-vision base reference, the 30-step arm | image, hashes and K3 pass (`~` steps are recorded, not blocking) |
| 4 | On the Mac: `pull.sh`; `xarch_hashes --compare`; `project_r1 --model qwen3-vl-2b …` | stop rule: STOP above 6.0 h projected (5.0 h until 2026-11-05) |
| 5 | If parity passed: set `sampler: vllm` in `r1_qwen3vl2b.yaml`. Peak memory > 22 GiB: micro-batch 2 x 16 | — |
| 6 | R1, then R2 and R3–R5 per §6.1 | gates §5.1 |
| 7 | `pull.sh`, `teardown.sh`, ledger; measured s/step into ROADMAP §Compute | instance absent |

### 4.8 The GPU (revision 9)

**One RTX 4090 (24 GB), US or Canada hosts only** (`policy.reasons()`). Measured price: $0.40–0.44 per hour (Quebec,
Texas; 2026-11-02). Why this card for Qwen3-VL-2B with LoRA: 4.3 GB of bf16 weights, adapter state 0.3 GB, activations
under checkpointing, logits at micro-batch 4. Peak should stay below 20 GiB (R0 measures). Fallbacks, in order: micro-batch
2 x 16 on the same card; an RTX 5090 (32 GB). A 48 GB card is for full fine-tuning only (D6). The RTX 3090 stays the cheap
box for dry runs. The 2026-10-19 snapshot of five cards (A3 §3.2) is in `lab/A3/vast_snapshot.json`; its prices are out of date.
Rented hosts never get RIGA (D16).

### 4.9 Phase 2 budget for rung 2 (revision 9, bottom-up)

Price: $0.42 per hour (measured). Credit on 2026-11-02: $3.18. Basis, all from rung 1 (RESULTS §1) except where marked:

- **Steps.** 771 per epoch at batch 32 (24,701 rows per epoch). Stable 1,542, decays 58 + 116 + 231 = 405: **1,947 steps**. **Owner 2026-11-05: 3 epochs** on `lidc-bundle-v3` (`build_views --epochs 4`; epochs 0-2 byte-identical to v2): stable 2,313, decays 58 + 116 + 231 + 347 = 752: **3,065 steps** (2.90 h at the R0 3.40 s). With the vLLM exam (D13 revision 10): **5.34 h projected; the stop rule rises from 5.0 to 6.0 h** by the same decision.
- **Step time.** LFM full FT: 0.69 s per step (measured). Qwen3-VL-2B LoRA with checkpointing does 6·N·T work, as LFM did. Per step: 32 samples x about 560 tokens (256 image + 60 prompt + about 240 answer, 1.5 x LFM's tokens per vertex) = 17.9k tokens; LM 1.72 B parameters: 184 TFLOP; vision tower (forward only): 26 TFLOP. LFM ran at about 70 TFLOPS. So **3.1 s per step** (range 2.4–4.5), 1.7 h (1.3–2.4 h). Measured by the R0 arm.
- **HF sampling** (assumes the R0 parity fails). LFM: 0.132 s per sample on base weights, 0.048 trained, 128 sequences per call, 16.5 ms per decode step. Qwen3-VL-2B: 28 layers, about 29 ms per step, 64 sequences per call, longest answer about 550 tokens: **about 0.25 s per sample** (range 0.15–0.40). The R0 rate step measures the base value.
- **vLLM** (if the R0 parity passes): the plan's allowance is 0.3 h (UNVERIFIED), taken as 0.4 h.

Samples in `select_ckpt` as chained today: 62,944 (3 branches x 256 x 8; then shifted 1,084 x 8, rescaled 1,084 x 8, rotated
1,084 x 8, shifted K = 16, both probes 1,680 x 8). Report-only parts: rotated 8,672 and K = 16 17,344 = 26,016. Without
them: 36,928.

| Block | GPU-h central (range) | $ at 0.42 |
|---|---|---|
| Launch, pull, data, teardown | 0.2 | 0.08 |
| R0: checks, parity, rate, reference, arm | 0.5 (0.4–0.7) | 0.21 |
| R1 training | 1.7 (1.3–2.4) | 0.71 |
| R1 evaluation, vLLM path | 0.4 (0.3–0.8) | 0.17 |
| R1 evaluation, HF path with the report-only sets cut | 2.6 (1.5–4.2) | 1.09 |
| R1 evaluation, HF path as chained | 4.4 (2.6–7.0) | 1.85 |
| Branch checks on CPU, idle gates | 0.25 | 0.11 |
| R2 G-probe (only if gates pass; GRPO memory unmeasured) | 0.5 | 0.21 |

| Scenario | GPU-h | $ central | $ high |
|---|---|---|---|
| vLLM path, no R2 | 3.0 | 1.3 | 1.9 |
| HF path, cuts, no R2 | 5.2 | 2.2 | 3.3 |
| HF path as chained, no R2 | 7.0 | 2.9 | 4.5 |

With $3.18 only the first two fit, and the second has no margin at the high end. `project_r1` stops R1 above 5.0 h
(about $2.1; training plus evaluation). Cuts that need a change in `vast/policy.py` / `select_ckpt` (not done here): drop
`--rotated` and K = 16. The cheapest lever is the R0 parity for this model.

**Caps.** Prepaid credit is the cap. The older alarm ($10) and hard stop ($15) of §4.6 apply to the whole of Phase 2; spend
so far is in `ledger.csv`. C2 (3 more trainings and evaluations) costs about 3 x $0.7 for training plus its evaluations: it
needs a new tranche in any scenario.

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
  assigned, so it precedes R-x8 and two conditionals cannot fire on one result. Single-seed runs (screening) carry no seed claim.
- **Multiplicity**: the only other decisions are pre-registered rules without multiplicity correction and with no
  claim attached: the C1 stop rule (gates) and D12 checkpoint selection. **C2 is the single confirmatory test** of Phase 2: its primary metric is A0′'s pooled hole
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
| **G-precision** | precision probe (D19), zero-shot, then the capacity run only on a zero-shot failure, 6–16 px objects, patient-disjoint probe-test patches | legal ≥ 95 %; median ASSD ≤ 0.6 px; median IoU ≥ 0.80; median centroid error ≤ 0.5 px; paired per-item IoU win over the centroid disc ≥ 60 % | A4 §6.2 table (σ = 0.5 → 95.7 % legal; σ ≈ 1 → ASSD 0.61; 1 px shift −23 IoU pts); B3 P2 discrimination (disc fails only via the paired rule: ≤ 49 % wins vs σ ≤ 0.5 at 78–100 %). Legal ≥ 95 % is knife-edge (σ = 0.5: 92.3 % gen. 1, 96.4 % gen. 2) | R1: **zero-shot** on the selected branch (inference only, 1,000 probe-test items, greedy + K = 8). Pass → done. Fail → the 20k-item capacity run (~0.3 GPU-h) decides: capacity pass → rung kept, LIDC SFT reported as the bottleneck (A4 §6.2 reading); capacity fail → the rung is dropped from LIDC. **A zero-shot fail is not final until the capacity run is done. Rung 1's capacity run was not run (RESULTS §1, correction).** |
| **G-topology** | topology probe (D20, B3 Addendum A), selected A0′ checkpoint: (a) simple controls with correct (β0, β1) = (1, 0) — **gate**; (b) **hole recall on one-hole items** (pooled; per hole-area bin) — **report-only for C1** (revision 4): a hole-recall miss never fails a rung or climbs the ladder; it is decided by C2 only, and R1's value is reported early | (a) ≥ 95 % (aligned with C2 clause (4), false non-simple ≤ 5 %); (b) reported against 0.30 | replaces "topology-any ≥ 90 % on unanimous items" (val has 0 unanimous non-simple items: vacuous; B3 P1). The σ = 0.5 synthetic model scores 0 false non-simple (B3 Addendum A4), so 95 % is not tight. 0.30 is the C2 fail line (A1 H-A1-6; B3 Addendum A6); per-bin estimates are reliable at n = 120: P(estimate < 0.30) = 0.010 at true recall 0.40, 0.998 at 0.20 (B3 A3) | R1 (selected branch) and R3–R5, K = 8 at T = 1 (greedy reported) |
| **G-signal** | at T = 1, G = 8, shifted val: duplicate share (hash of format/legality flags + packed rasters, B2-6) and share of groups with geometric σ above its floor | duplicates ≤ 10 %; σ-share ≥ 50 % | RESULTS §8; DAPO (equal-reward groups give zero advantage); K12 smoke 0 % / 100 % on a non-collapsed checkpoint | every WSD branch (D12 eligibility; a floor, not a ranking key: σ-share and duplicates were saturated in B2 §2.1) |
| **G-vision** (revision 5) | "the model still looks at the image", B1 §2.5 protocol promoted to `eval/vision.py` (M4): (a) mean pairwise cosine distance of the projector output (`get_image_features`) over a fixed image set (4 centred LIDC val crops + a black image, as `lab/B1/vision_feats.py`), as a ratio to the base model's value measured in R0; (b) teacher-forced next-token TV over the first 24 target tokens of 8 fixed LIDC val items, real crop vs a black image and real crop vs another item's crop (`lab/B1/image_use.py`) | (a) ratio ≥ 0.5; (b) both TVs ≥ 0.002 | B1 §2.5: collapsed checkpoints (Qwen3.5 and SmolVLM, all-linear LoRA, 36 / 150 steps) have ratio 0.001–0.064 and TV ≤ 0.0007; intact ones (LFM all-linear 36 / 150 steps; LM-only Qwen3.5 and SmolVLM) ratio 0.68–1.0 (frozen = 1) and TV ≥ 0.004. (a) sits 1.4× below the lowest intact ratio and 8× above the highest collapsed one; (b) 2× below the lowest intact TV and ~3× above the highest collapsed one. B1's diagnostics are post-hoc and small (4–8 items), and TV is diluted by syntax tokens, so (b) holds only under this exact protocol | every WSD branch in `select_ckpt`, on the CPU as in B1, so ② stays GPU-free (minutes per branch; D12 eligibility); (a) and (b) also in `train_sft` at each rolling checkpoint: both below threshold → abort and D7 repair; one alone → warn-only, logged |
| **G-probe** | "GRPO can improve this checkpoint": 20 GRPO steps (`grpo_config` with a **pre-registered probe LR**: **5e-6 full FT, 5e-5 for LoRA**, constant; G = 8, 32 LIDC train prompts per step with ×4 disagreement oversampling, HF generation, `use_vllm=False`); fixed G-probe set = 128 shifted-val LIDC items, K = 8 at T = 1. **Noise floor at zero training cost:** at step 0 the set is sampled **twice** with different sampler seeds; the per-item difference gives the null CI (half-width h₀). Score per item: the **raw composite** R = 0.3·IoU_tc + 0.3·F5 + 0.1·legal + 0.3·g_area (v4 weights on the raw, legality-gated terms; TRL's aggregator is group-mean-centred, so its mean is ≈ 0) | **PASS**: paired per-item gain (R₂₀ − mean of the two R₀ draws) CI lower bound > 0, legal at T = 1 falls by ≤ 5 pts, σ-share at step 20 ≥ 50 %. **FAIL** (counts): the gain CI lies entirely below the null CI (R degraded beyond sampling noise), or legal falls > 5 pts, or σ-share < 50 %. Otherwise **INCONCLUSIVE**. Reported too: h₀, and the training-log form (mean logged R, steps 11–20 vs 1–5) | G-signal is necessary, not sufficient (2606.18487: advantage variance p(1−p)(g−1)/g can vanish). The probe LR is **the probe's, not Phase 3's**: production `grpo_config` uses 1e-6 (full FT) / 1e-5 (LoRA), which over 20 steps may not move R beyond h₀; SimpleSeg's RL used a constant 2e-6 for ~3k steps at 7B (Tables 4–5); 5× the full-FT production LR stays 10× below the SFT peak (5e-5) and LoRA keeps its 10× ratio (LoRA Without Regret). Thresholds, null procedure and LR pre-registered in `prereg.py` (M8) | ④ in §3.2, on the selected branch. Not PASS → the D12(c) repair (next-earlier eligible branch), probed once. **Still not PASS → no ladder climb**: the rung passes C1 with a "G-probe warning" handed to Phase 3 (which starts with its own LR choice and a longer run), because D12 eligibility already requires G-signal; the ladder climbs only when no branch passes G-signal (a §5.1 gate failure) |
| **G-x86** | K3: 256-rollout batch wall (best of serial / pool), RIGA-only (synthetic) and mixed; K4: p99 per rollout; hashes serial = pool | K3 ≤ 250 ms; K4 ≤ 5 ms; equal hashes | RESULTS §8; B4 (Mac pool-8 51–61 ms; A6 serial RIGA 241 ms at load 9) | `k3_harness` in R0 (idle) and once under R1's load (§4.7 step 6) |
| **G-cost** | measured s/step, SFT and GRPO, per (model, card); projected Phase 3 Track A GRPO cost | measurement recorded; projection ≤ $150 on-demand, else coordinator re-plans (card, arms, seeds) | RESULTS §8 cost gate; A3 §3.4 (~$105 for ~70 GPU-h at $1.5/h) | R0 timing prelude and R1 logs (SFT); R2 (GRPO, HF generation); colocate step time measured at Phase 3 start |
| **G-parity** | vLLM sampler vs HF: prompt ids; per-token \|Δ logprob\|; distribution parity (32 items × K = 8). Measured per model in R0. LFM failed (p99 per sequence 1.11 nats), so rung 1 used HF | 100 % ids; median ≤ 0.02 nats bf16; every metric difference inside its bootstrap 95 % CI (G-hf-check on the trained branch). p99 per sequence ≤ 0.5: report-only since 2026-11-05 (D13 revision 10) | A5 H-A5-1/2; B5 Q2 (fp32 1e-7–4e-6) | R0 `r0_parity` for Qwen3-VL-2B. Fail → HF is the sampler for that model (D13). Pass → R1 may set `sampler: vllm`. Qwen3-VL-2B: ids and median pass, sequence p99 1.79 (report-only) → R1 vLLM with G-hf-check |
| Housekeeping | LFM licence recorded; novelty sentences narrowed (P1 §8); RIGA permission status recorded | done | RESULTS §8 | coordinator |

"Near MedSAM" (A4 §3.2) is **reported, not gated**, with every row labelled by view: model vs MedSAM (oracle box from
that view's union) − 0.10 per structure (centred: LIDC 0.68, disc 0.84, cup 0.74), on the centred **and the shifted**
val, next to the jittered-box MedSAM row, SAM 2.1-tiny and the prior on the same views (D17).

### 5.2 Claim tests

| Claim | Pre-registered procedure | Pass | Evidence for thresholds |
|---|---|---|---|
| **C1 smallest viable** | rungs screened smallest first with the one recipe bundle (§6.0); select per D12; §5.1 on LIDC (and RIGA if permitted); **stop at the first rung that passes**. Rung 1 (LFM) failed on data with a known gap (D10, revision 9), so its failure is weak evidence about model size, and a pass at 2B does not show that 450M cannot pass | C1 = the first (smallest) rung passing every §5.1 **gate** (G-prior (a) is report-only; G-topology (b) is report-only, so a hole-recall miss goes to C2's "holes left to GRPO", not up the ladder; G-probe never climbs the ladder on its own, §5.1). No larger rung is run once one passes: passing is sufficient for "smallest viable", and a larger model's score adds nothing to it. Qwen3-VL-2B is rung 2 (§6.1) | P1-D12; ROADMAP Phase 2; owner directive (revision 3) |
| **C2 A0 vs A0′ (holes; confirmatory)** | C1 rung, bundle; depth = that of the checkpoint that received the final G-probe verdict (the D12(c) repair branch if it was probed, not necessarily the first D12 pick); A0′ = bundle with ×4 non-simple rows, A0 = the same rows with holes emptied; A0′ seed 1 = that checkpoint, **A0′ seeds 2–3 and A0 seed 1** trained to that depth (R3–R5): **3 A0′ seeds, 1 A0 seed** (revision 4: the A0 arm scores 0 hole recall by construction and enters no hole decision, while the A0′ seeds are the scarce resource, B3 A5); topology probe D20 (360 one-hole items, 120 per bin); primary metric and bootstrap as D20 / §5.0. A0 is run under the same protocol; on hole items A0′ − A0 equals A0′'s recall (A0 = 0 by construction), so it is reported for completeness only | **Pass iff all hold** (B3 Addendum A6 adapted to 3 seeds, revision 4): (1) A0′ pooled hole recall (mean over K = 8 and the **3 seeds**, within-bin bootstrap) ≥ 0.50 with 95 % CI lower bound > 0.30; (2) **≥ 2 of the 3 seeds' pooled point estimates ≥ 0.50 and none < 0.30**; (3) no bin's point estimate < 0.30; (4) false non-simple on controls ≤ 5 %; (5) legal ≥ 90 %. Reconciliation with B3 A5's caveat (item CI conditional on the trained checkpoints; 54 % coverage of the seed mean at σ_v = 0.3): the item-bootstrap CI is kept for item noise, and seed variation is handled by the per-seed point-estimate clause; with 3 seeds a majority rule tolerates one weak seed while the < 0.30 floor keeps the fail line; the across-seed t-interval of the 3 seed means (2 df) is reported, not gated. Seed-discordant flag (seeds on both sides of 0.50) is reported. **A0 guard (pre-registered, report-only):** paired (A0′ 3-seed mean − A0) shifted-val median IoU_tc CI lower bound > −0.02, else "hole-awareness costs geometry" is reported. **Outcomes (revision 7), each with at most one conditional run, owner approval required:** **pass** (all of (1)–(5)) → none; **partial** (not pass, pooled A0′ recall in [0.30, 0.50), or ≥ 0.50 with another clause failed) → one ×8 A0′ seed (E1, recipe question), not run when clause (4) false non-simple failed (more hole rows would raise it); **fail** (pooled < 0.30) → "holes left to GRPO" (R-ctrl is removed: the rung is already the 2B model; a larger model only on the owner's order). ×8 never follows a fail. Part recall is descriptive.  | B3 Addendum A3 (±3.1 pts pooled; per-bin floor reliability), A4 (σ = 0.5 px synthetic model: pooled 0.76 [0.73, 0.78], bins 0.65 / 0.81 / 0.90), A5 (seed caveat), A6 (rule); A1 H-A1-6 (≥ 50 % / < 30 %); A2 H-A2-5 (hallucination); JUDGE_R1 MUST-FIX 1, 3 |

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
| Crop shift + scale jitter (D10, **revised in revision 9**) | shift uniform, union inside, ≥ 2 px clearance, no other nodule; S uniform over the policy side and its ladder neighbours (LIDC 50 / 100 / 125), the sides of the shifted and rescaled exam views | mandatory: the centred frame leaks location and size: the image-blind prior passes RIGA's centred gates and collapses under shift (0.00 / 0.32 / 0.07; A4 S5–S6), and every B1 arm trained on centred crops drew at the frame centre (94–98 % of shifted samples closer to the frame centre than to the nodule, B1 §2.5); the gates are defined on shifted / rescaled views |
| Offline dihedral targets (D9, D10) | D4 for LIDC (RIGA: hflip only); one op per rater row; crop-before-trace; per-view hash | 215,208 variant targets 100 % certified, LIDC pixel-exact, byte-identical rebuilds (B3 P3–P4); exact transforms are standard for LIDC (A4 §4.1); zero runtime cost |
| CT window jitter (D10) | L −300 ± 60, W 1400 ± 80 HU on **int16 HU context crops** (one npz per LIDC item, written by M1 from `data/store`), noise σ ≤ 10 HU, windowed at collate time | random windowing for CT (liver tumours, 2510.08116; A4 §4.1): one study, the weakest evidence in the bundle, kept because it is applied before windowing (certification untouched), costs ≈ 0 and is the §6.3 rule-3 removal; at zero jitter the HU path reproduces the prepared PNG byte-for-byte (M1 test) |
| All voting raters (D9) | one row per voting rater, every epoch | same expected gradient as random-rater sampling, keeps the rater distribution D5 matches (A2 §4.6; P1 I27) |
| Hole oversampling (D11) | non-simple LIDC train targets ×4 (414 of 23,459 targets → ≈ 6.7 % of rows) | A2 §7.4 and A1 §5.2 recommend ×4–8; P1 used ×4 at 3.9 % prevalence (RESULTS §7 #9); the low end limits hole hallucination (A2 H-A2-5 guard, checked by G-topology (a)) |
| Vision tower + projector (D7) | LoRA: frozen; full FT: 0.1× LR; G-vision on every branch, in-training projector check | all-linear LoRA incl. vision collapsed Qwen3.5's and SmolVLM's image features in 36 steps, LM-only LoRA kept them at equal legality / IoU (B1 §2.4–2.5); A2 §4.2 (unfreeze for localisation; 0.1× its own proposal) for full FT |
| Fine-tuning (D6) | full FT for LFM / SmolVLM; LoRA r16 on the language model for Qwen3-VL-2B (full FT needs ~34 GB) | A2 §4.1; B1 §2.4; memory math in D6 |
| Per-model settings | Qwen3-VL-2B: TRL's `qwen3_vl_training` template (swapped in by `SFTTrainer`), bf16 frozen base, `MAX_COMPLETION` 1,536, micro-batch 4 x 8, checkpointing on; every sampling field explicit (T 1, top_p 1, top_k 0, min_p 0, repetition 1, presence 0), overriding the card (D13) | B1 §1.4; checked on the real processor 2026-11-02 |
| Schedule (D5) | WSD, 3 % warmup, stable 2 epochs, 15 % decay branches at 0.5 / 1 / 2 epochs; LR 5e-5 (LoRA 5e-4); effective batch 32 | A2 §4.3; Quagmires (1–2 epochs) |
| Throughput (D4) | `group_by_length` + `length` column; 6 persistent workers; grad-ckpt off when peak ≤ 20 GB; SDPA; bf16 | padding 7–32 % → ≤ 0.2 % (A2 §1.2); collation 10–25 ms/row needs workers (A6 §1a, §2); measured in R0, not A/B-tested |
| Loss (D8) | assistant-only NLL; **no number-token loss** | NTL's evidence is text-only, untested on VLM coordinates, weakest on LFM's 3-digit tokens; it stays a triggered option (D8, §6.4) |
| Not in the bundle | label smoothing, DFT, constrained decoding, coordinate tokens, mixing | D8, D13, D2, D11 |

### 6.1 The Vast runs (revision 9)

Ladder (D1): rung 1 LFM2.5-VL-450M **done, failed** (RESULTS §1) → **rung 2 Qwen3-VL-2B** → stop and re-plan (§8.2).
All on LIDC, one 4090 session (§4.7). Every rental needs `make preflight` to pass first.

| Run | Purpose | Design | Decides / measures | Est. (4090) |
|---|---|---|---|---|
| **R0** | box readiness | `configs/r0.yaml`: image check, x86 hashes, idle K3, parity and HF rate for Qwen3-VL-2B, G-vision base reference, 30-step arm | G-x86, G-parity, peak memory, s/step, projected R1 hours (`project_r1` stop rule) | 0.5 h |
| **R1** | **C1** screening of rung 2 | `configs/r1_qwen3vl2b.yaml`: new bundle (D10 revision), A0′, seed 1, WSD 2 epochs + branches 0.5 / 1 / 2, LoRA r16; then `select_ckpt` (shifted val K = 8 per branch; the selected branch also gets the rescaled val, both probes, K = 16 and the rotated check) | §5.1 gates; D12 selection | 3.0 h vLLM path, 5.2 h HF path with cuts (§4.9) |
| **R2** | "GRPO can improve it" | G-probe on the selected checkpoint (§5.1), only if the gates pass; the D12(c) repair branch only if not PASS | G-probe; GRPO step time | 0.5 h (est.) |
| **Stop** | — | a passing rung is C1; no larger rung is run | C1 | — |
| **R3–R5** | **C2** (confirmatory) | A0′ seeds 2 and 3, A0 seed 1 (same ×4 rows, holes emptied), each to the depth of the checkpoint with the final G-probe verdict; `c2_*.yaml` (already pointed at rung 2) | topology probe → §5.2 3-seed rule | 3 x (R1 training) + evals |
| *R-x8* | conditional | only if C2 is **partial** (§5.2) and clause (4) passed: one A0′ seed with non-simple rows ×8 (`rx8.yaml`), owner approval; exploratory | whether more hole rows lift recall | 1 seed |
| *R-ctrl* | removed | Qwen3-VL-2B is now the rung itself. If C2 fails, "holes left to GRPO" (§5.2); a larger model only on the owner's order | — | — |

The precision capacity run (~0.3 GPU-h) runs only when the zero-shot G-precision fails (§5.1). It was not run for rung 1.

### 6.2 Measurements folded into those runs (no separate experiments)

| Measurement | Folded into | How | Cost |
|---|---|---|---|
| Reward K3/K4, x86 hashes | R0 | `k3_harness` idle, `xarch_hashes` compare | minutes |
| Peak memory and s/step of Qwen3-VL-2B | R0 arm | 30 discarded steps at R1's settings; the longest batch comes first | 0.06 h |
| G-vision base reference | R0 | `train_sft --vision_reference_only` on the base model | CPU |
| Sampler parity for this model | R0 | `r0_parity` (ids, logprobs on 8 items x 8); a pass allows vLLM in R1 | 0.1 h |
| HF sampling rate on base weights | R0 | `hf_rate` (32 items x K = 8) for `project_r1` | 0.15 h |
| SFT s/step, GPU use | R1 logs | s/step, data wait | 0 |
| B2 signal panel per branch | R1 `select_ckpt` | CPU scoring of samples already drawn (feeds D12 R-B) | 0 GPU |
| Depth curve, 3 points | R1 | decay branches 0.5 / 1 / 2 and the stable checkpoints at those points | disk |
| Does the signal predict GRPO gain | R2 | per-item Spearman of step-0 sd_geom / r_std vs gain | 0 |
| GRPO step time and split | R2 | generation / scoring / update | in R2 |
| WiSE-FT 0.5-epoch to selected | R2, only if the D12(c) repair fires | merge at α = 0.5, report-only | 0.05 h |

Done and not repeated (rung 1, RESULTS §1 and V0_V1): the 4 timing arms (workers 6, determinism off), the LFM parity and
HF rate, the watchdog dry run.

### 6.3 After a failed gate (revision 9: the automatic fallback is removed)

Revision 8 had a pre-registered fallback: one leave-one-out run, chosen by a rule. Rule 1 (hole oversampling x4 to x1)
fired on G-topology (a). Rung 1 shows it targets the wrong cause. Of the simple control items, 4.3 % were illegal, 1.3 % were
unreadable and 0.7 % had a false hole (total 6.3 % wrong, gate 5 %). Oversampled holes explain at most the 0.7 %. The owner
rejected the run. Rules 2 and 3 (dihedral off, window jitter off) had no evidence either: the rung-1 diagnosis points to
crop size and detection (RESULTS §1). They are removed too, with their configs.

New rule: a failed gate gets a Mac diagnosis of the same kind as RESULTS §1 (by crop size, by object size, against MedSAM
with the expert box, gallery). The owner then decides the single next run. The bundle's never-removed parts stay: crop
shift and scale jitter, all voting raters, the WSD branches.

### 6.4 Removed or triggered (owner directive, revision 3)

| Item (rev 1–2 ID) | Status | Trigger that brings it back |
|---|---|---|
| E1 hole oversampling ×1 / ×4 / ×8 | removed; ×4 fixed in the bundle | R-x8: one ×8 A0′ seed only on C2's **partial** outcome (§5.2, §6.1) and the owner's approval |
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
| Unconditional Qwen3-VL-2B control; screening of larger rungs | Qwen3-VL-2B is now the screened rung (revision 9) | a larger rung only after the 2B fails (§8.2) |
| nnU-Net v2 training (D17 row 4) | removed from Phase 2 | Phase 9 paper baselines |
| K = 64 evaluation, I35 audit recalibration | moved | before the first Phase 3 GRPO run |
| HD95 / ASSD / centroid in `objdel.eval` | removed | paper tables (the probes keep their own) |
| H-A5-5 constrained diagnostic, H-A5-4 caching study | removed | Phase 3 planning, or G-format parse < 98 % |
| S7 Phase 3 environment smoke (colocate) | moved | first step of Phase 3 |
| H-PULL resume test | removed | first host failure (resume procedure documented, A3 §6.2) |
| MedSAM / SAM 2.1 GPU runs | moved to the Mac CPU (inference only) | — |
| S6 RIGA runs | removed from the critical path | written RIGA permission (D16) and owner approval of a tranche |

## 7. Milestones and task graph (test-driven: tests first, then promote lab code, then the gate)

All milestones are built (Stage 1, 2026-10-20 to 11-02). This table is the record of what was built; Qwen3.5 items in it are history (revision 9).

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
| Qwen3-VL-2B also fails G-geometry (the 0.60 gate was set on centred frames; human LOO 0.768) | report shifted and rescaled numbers; Mac diagnosis first (RESULTS §1 form); §8.2 |
| Legality at ≥ 90 % is knife-edge for small nodules (σ ≥ 0.6 px vertex error suffices to fail; A4 S4) | E2 / E10 triggers (§6.4); precision probe; D12 selection on T = 1 legality |
| Shifted windows include another outlined nodule (ground truth then misses it) | `meta.other_nodules` check in the shift draw; M1 test |
| Qwen3-VL-2B peak memory: 152k-token logits, 28 layers, 256 image tokens per sample | LoRA, micro-batch 4 x 8, checkpointing; the R0 arm sees the longest batch first (length grouping) and measures the peak; fall back to micro-batch 2 x 16. G-probe (GRPO, micro-batch 8 = G, 1,536-token completions) is not measured: set its batch and checkpointing in R0 before R2 |
| HF sampling is slow for 2B (64 sequences per call, ~29 ms per decode step, a batch lasts as long as its longest answer) | R0 parity may allow vLLM (D13); else cut the report-only sets (§4.9) |
| Vision collapse during SFT (B1: all-linear LoRA erased Qwen3.5's and SmolVLM's image features within 36 steps, B1 §2.5; the tower is frozen under LoRA) | LoRA on the language model only, full-FT vision at 0.1× (D7); in-training G-vision check (abort only when (a) and (b) both fail); G-vision on every branch; D7 repair (frozen) once per rung |
| vLLM/HF parity fails for a model | HF sampler is the gate path (slower; A5 H-A5-8) |
| Reward pool or fork misbehaves on Vast (threads at fork, OpenBLAS/libgomp) | pool before CUDA, `cv2.setNumThreads(0)`, `ProcessPoolExecutor` only; C++ scorer trigger (D14) |
| Over-trained SFT collapses diversity (no GRPO signal) | WSD branches; D12 R-B with the sd_geom veto; next-earlier-branch repair; R2's check of the veto statistic (B2 found no over-training at Mac depth, but its checkpoints were under-trained, B2 §6) |
| Holes not learnable at 1.35 % prevalence | ×4 oversampling in the bundle; C2 fail branch: holes left to GRPO |
| RIGA permission never arrives | Phase 2 exits on LIDC; RIGA arms (band claim) wait; Mac MPS fallback for RIGA SFT only |
| Phase 3 re-pinned env drifts (vLLM 0.33 moves torch to 2.14.0 exact pin, not our 2.14.1) | re-run fast suite + strict hashes per re-pin; re-evaluate at vLLM 0.33 (B5 §What this changes) |
| Blackwell stack issues (sm_120, CUDA 13, driver ≥ 580) | A100 80 GB fallback; V1 on Ada |
| Host failure mid-run | pulls every 15–30 min; resume tested (H-PULL) |
| Forgotten instance / runaway spend | G1–G10 |
| The probe's smallest-bin holes are ~1 px larger than real ones (P5.2 failed narrowly) | stratified design confines the bias to the ≤ 10 px bin; per-bin floor in C2; the bin is reported with the caveat; the real val slice is reported via TV against the raters' hole fraction |
| G-probe is noisy at 20 steps / 128 items | paired per-item design (§5.0); one re-probe after the D12(c) repair; the training-log form reported alongside |
| One bundle, no per-component A/Bs: a harmful component is not isolated in advance | a failed gate gets a Mac diagnosis (§6.3); §6.4 triggers bring back the specific A/B only when its failure mode appears |

### 8.2 Kill / continue

- **Rung 2 fails a gate.** Diagnose on the Mac first (crop size, detection against MedSAM, gallery: RESULTS §1). The owner
  then picks one change: unfreeze the vision tower (D6), Qwen3-VL-4B or InternVL3.5-1B-HF, a data change, or stop. No automatic run.
- **Precision.** A rung that fails zero-shot G-precision gets the 20k-item capacity run (~0.3 GPU-h) before it is judged.
- **Format.** If no affordable rung passes G-format on A0′ (ROADMAP kill rule "A0′ < 90 % valid"): E10 or E2 by their
  triggers, then E8 coordinate tokens, then a coarser grid G = 100 for LIDC only (lossless on LIDC: A1 §1.5). Mask tokens
  are rejected (no polygon, no certificate; A1 §6). A rung that passes format but not G-geometry goes on to Phase 3 only
  if it passes G-prior (b), G-signal and G-probe; the miss is then a Phase 3 baseline (RL adds about +10 gIoU over SFT in
  SimpleSeg Table 3; A1 §3).
- **After G-probe.** If the selected branch and the D12(c) repair branch are both not PASS, the rung still passes C1 with a
  G-probe warning for Phase 3.
- **After C2.** A partial or fail does not block Phase 3; holes are then a GRPO question (§5.2).
- **Budget.** Prepaid credit is the cap ($3.18 on 2026-11-02). Every extra tranche, E10, RIGA and each triggered
  A/B need the owner's approval. Alarm and hard stop: §4.6.

### 8.3 Open questions for the owner

1. **Credit.** $3.18 covers one rung-2 session only on the cheap path (§4.9). Add credit, or accept the cuts listed there.
2. **Sampler.** Is 0.5 nats per summed bf16 sequence too strict for parity (it failed LFM at 1.11)? R0 measures Qwen3-VL-2B.
3. **Vision tower.** Keep it frozen (LoRA, 24 GB card) or pay for a 48 GB card and train it (D6)?
4. **RIGA permission** (I44): until a written yes, RIGA stays off Vast (D16). **Commercial use:** RIGA, MESSIDOR and the LFM
   licence limit it (A1 §2, A4 §5.2).
5. **Image diet.** The image still holds Qwen3.5-only extras (D15). Rebuild without them (smaller pull) or leave.

## 9. Decisions log and open items

### 9.1 Log

**2026-10-20, draft C to revision 8.** Eight revisions on one day, each answering a report or a judge round. What stays true:
- Draft C: A1's ladder; gates on shifted val; RIGA needs written permission and Secure Cloud; deviations from the notes
  (label smoothing rejected, vLLM 0.32 sampler, `cpu_cores_effective ≥ 12`).
- Revision 1 (JUDGE_R1): C2 rescoped to holes. Revision 2 (B3 Addendum A): hole-calibrated topology probe, C2 rule.
- Revision 3 (owner cost directive): one recipe bundle (§6.0), no per-component A/Bs, WSD 2 epochs with branches
  0.5 / 1 / 2, Vast runs only for C1, C2 and G-probe, stop at the first passing rung, baselines on the Mac CPU.
- Revision 4 (JUDGE_R2): G-topology (b) report-only; G-probe calibrated; C2 with 3 A0′ seeds and 1 A0 seed.
- Revision 5 (B1): vision tower frozen under LoRA, G-vision gate (§5.1), SmolVLM control only.
- Revision 6 (B2): D12 rule R-B and the signal panel (M11). Revision 7 (JUDGE_R3): G-prior (a) report-only; C2 outcomes
  pass / partial / fail with R-x8. Revision 8 (JUDGE_R4): approved.

**2026-11-02, revision 9 (rung 1 failed; owner decisions).** Details in the box at the top and in RESULTS §1.
- Evidence: LFM2.5-VL-450M fails G-geometry (median IoU_tc shifted 0.495, rescaled 0.092), G-precision (zero-shot only)
  and G-topology (a) (controls 93.7 %). MedSAM with the expert box scores 0.77 / 0.80, so the policy fails at finding
  the nodule in a large crop, not at drawing. Train rows with a nodule under 10 % of the crop: 0.5 %; the rescaled exam: 27 %.
- Decisions: Qwen3-VL-2B is the next rung (LoRA, one 4090). No LFM rerun, no Qwen3.5. CUDA checker now. `make preflight`
  before every rental. Training data rebuilt (`lidc-bundle-v2`). §6.3 rules dropped. Budget re-estimated (§4.9).

### 9.2 Judge rounds

The four judge rounds (`JUDGE_R1.md` ... `JUDGE_R4.md`) and the note inconsistencies they resolved are in those files and in
`RESEARCH.md`. This plan keeps only their results (the revisions in §9.1).

### 9.4 Open items

- I13, I15/I16, I35 (audit recalibration, Phase 3 entry), I38 (Phase 3), I44 (RIGA permission), as before.
- Rung 2 (Qwen3-VL-2B): peak memory and s/step at micro-batch 4 (R0 arm); vLLM parity for this model (R0); the
  frozen vision tower may limit localisation (D6); G-probe memory with a 152k-token vocabulary (§8.1).
- The data fix (D10, revision 9): the coverage check of the new bundle against all three exam views.
- Open since rung 1: G-precision capacity run (§5.1); holes are C2 / Phase 3 work by design (§5.2).

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
- Models: https://huggingface.co/LiquidAI/LFM2.5-VL-450M ,
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
- PEFT LoRA `exclude_modules` (D6, D7): https://huggingface.co/docs/peft/main/en/package_reference/lora
