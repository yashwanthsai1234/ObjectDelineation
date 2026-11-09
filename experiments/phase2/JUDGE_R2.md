> Frozen record (written 2026-10-20). Versions, budgets and paths were true then. Current facts: docs/README.md.

# Judge verdict, round 2: Phase 2 plan (revision 3, 2026-10-20)

Independent review of `experiments/phase2/PLAN.md` revision 3 against `JUDGE_R1.md`, `RESEARCH.md`, notes A1–A6,
`lab/B3/REPORT.md` (incl. Addendum A), `lab/B4/REPORT.md`, `lab/B5/REPORT.md`, `lab/B1/{prereg,analyze,common}.py`
(decision rule only; pilot running; no `lab/B2/` exists yet), `experiments/phase1/RESULTS.md` §8, `ROADMAP.md`
(working-tree diff), the source (`train/configs.py`, `eval/gates.py`, `data/store.py`, `targets/frame.py`), the lab
JSON (`hole_census.json`, `power_c2.json`, `topo_probe_v3e_seed{7,8}.json`, `census.json`, `vast_snapshot.json`,
`estimate_cost.py` re-run) and external sources (arXiv 2510.01624, 2606.18487, 2601.19228, 2510.08116, 2308.01825,
2312.06585; PyPI vllm 0.32.0 / trl 1.16.0 / torch 2.14.1; HF cards Qwen3.5-0.8B, LFM2.5-VL-450M; TRL #6294, #7198,
PR #7799; vLLM #34948; Thinking Machines LoRA post).

Judged under the owner's cost directive: Vast only for C1 screening (stop at the first passing rung), C2 (2 seeds),
G-probe; one recipe bundle, no per-component A/Bs, at most one leave-one-out diagnostic.

## Overall verdict: APPROVE WITH CHANGES

| Dimension | Score | One line |
|---|---|---|
| Scientific validity | 7 | C2 is now a calibrated, pre-registered, single confirmatory test; but G-topology (b) silently became a C1 kill gate, G-probe is uncalibrated at the production LR, and C2 spends two of three runs on an arm that enters no decision |
| Evidence | 8 | Lab numbers check out to the JSON (hole census, power, seed-8 sanity, prices); six external claims need small corrections (vLLM #34948 is fixed, 2510.08116 is liver CT, #6294 numbers, A1 §1.5 grid, SimpleSeg LR, Qwen3.5 card) |
| Engineering | 8 | Sequencing, bridge, preload-per-dataset and key delivery are fixed as asked; one real gap: CT-window jitter needs int16 HU crops and the bundle/M1 still inherit B3's 8-bit PNG crops |
| Feasibility / cost | 8 | $2.3–4.7 scenario A is credible in dollars with the $15 stop; GPU-hours are likely 1.5–2× low (MFU 0.10 floor is unmeasured for any ≤ 2B VLM, A6 §7), which the R0 stop rule absorbs |
| KPIs | 7 | Every gate is measurable; two logic holes (hole-recall as a rung gate; one-seed R-ctrl against a two-seed rule) and one threshold mismatch (controls 90 % vs 95 %) |
| Simplicity | 8 | 3 + 3 runs, one bundle, one diagnostic: the directive is met without losing C1 or C2; §9.2 still carries revision-1 wording in three rows |

## Answers to the four questions

1. **Is the bundle evidence-backed?** Mostly. Crop shift + scale jitter (own data: prior 0.904/0.783 → 0.00/0.32/0.07
   under shift), offline dihedral targets (215,208 certified, byte-identical), all voting raters (A2 §4.6 argument),
   ×4 hole oversampling (P1 precedent, A1/A2 ×4–8), WSD + 5e-5 (A2 §4.3; SimpleSeg's 5e-5 is its SFT LR, verified),
   NTL left out (A2 itself only proposed an A/B; LFM's 0–999 tokens are the worst case): all justified. Two components
   are weaker than the table says: the vision tower "unfrozen at the same LR" (A2 §4.2 says unfreeze, LR undecided,
   0.1× arm in H-A2-6; Qwen's own recipe and Axolotl argue low/frozen) and CT window jitter ("standard" rests on one
   liver-CT paper, 2510.08116, and needs HU plumbing that does not exist yet, MUST-FIX 5).
2. **What is lost by cutting the A/Bs?** Attribution only, plus the ability to see a harm smaller than a gate margin.
   Harms the gates *do* catch: hole hallucination (G-topology (a), C2 (4), IoU_tc = 0 on wrong topology), a location
   / size leak (G-prior, rescaled view), legality loss (G-format, D12 on T = 1 legality), diversity collapse (G-signal,
   G-probe). Harms they *do not* catch: a few IoU points from dihedral capacity on a ≤ 0.8B model, from window jitter
   or from a too-hot vision LR, as long as the rung still clears 0.60. Under the directive that is acceptable; the
   one free insurance is SHOULD-FIX 6 (score the selected checkpoint on the identity vs rotated val views, inference
   only). The §6.3 leave-one-out rule targets the right component per failure mode.
3. **Are the runs sufficient?** C1: yes, once G-topology (b) stops being a rung gate (MUST-FIX 1) and SmolVLM cannot
   be C1 (MUST-FIX 4). C2: yes in design, but the A0 arm enters no decision and the plan's own caveat (54 % coverage
   at σ_v = 0.3) says the A0′ seeds are the scarce resource (MUST-FIX 3). Phase 3 handoff: the checkpoint, GRPO step
   time (HF path, conservative) and the re-pinned image (B5 26/26) are enough; G-probe needs a noise floor before its
   failure can climb the ladder (MUST-FIX 2).
4. **Is the budget credible?** In dollars yes: scenario A $2.3–4.7 at the $0.66 median (snapshot: 4090 on-demand
   0.506 / 0.656, n = 23, verified in `vast_snapshot.json`), ≤ $7.3 at ceilings, hard stop $15; scenario B stays
   under $15 at the high end. In hours it is optimistic: no note has a measured s/step for any ≤ 2B VLM (A6 §7), the
   0.10 MFU floor is an assumption, ×4 oversampling adds ~5 % rows, and the 0.3 h eval block (3 vLLM loads, ~26k
   samples, probes) is unverified. At MFU 0.05 scenario A is ≈ 6–12 GPU-h, still ≈ $4–8. The R0 stop rule (projected
   R1 > 2× high estimate) is the right guard; state the MFU 0.05 row so the owner sees it before V1.

## R1 MUST-FIX status

| R1 item | Status | Check |
|---|---|---|
| 1 C2 hole-only, probe uncalibrated | **Resolved** | B3 Addendum A pre-registered (P5 in `prereg.py`); numbers match `hole_census.json` (285 holes, p50 12 px, bins 125/103/57 = 43.9/36.1/20.0 %, wall p50 1.236, 182/241 single-rater = 75.5 %) and `topo_probe_v3e_seed8.json` (σ = 0.5 mean-8 0.649/0.813/0.904, pooled 0.759 [0.733, 0.784]; A0 0/0/0). P5.2 recorded honestly as "failed narrowly" (seed 7 KS p 0.0037 on area); residual confined to the ≤ 10 px bin. ROADMAP C2 wording changed. Minor: B3 A1 "holed part 23.6 px (all splits)" is the all-splits value; train p50 is 22.7 |
| 2 No "GRPO can improve" KPI | **Resolved in form, not calibrated** | G-probe defined with three thresholds, pre-registered; but see new MUST-FIX 2 |
| 3 Seed pooling / primary metric | **Resolved** | §5.0 and D20 fix per-item mean over seeds and K, within-bin bootstrap (`power_c2.json`: mean-8 ±3.1 at σ_u 1.5, coverage 0.54 / 0.30), both-seeds rule, discordance flag, single confirmatory test |
| 4 Baselines on the gate view | **Resolved** | every row labelled by view; nnU-Net removed from Phase 2 with the jittered-train spec kept for Phase 9 |
| 5 Trainer vs vLLM on one GPU | **Resolved** | `train_sft` exits; `select_ckpt` CUDA-free launcher; sampler asserts free ≥ 90 %; M8 smoke asserts all three |
| 6a Quagmires reversed | **Resolved** | abstract verified: "unique examples for one epoch underperforms … half examples for two epochs"; D5 no longer uses "unique > repeated" |
| 6b B3 P1 vs P2 | **Resolved** | cited to REPORT §3 P2 |
| 6c 2 vs ~4 epochs | **Resolved in D5/§4.9, stale in §9.2** | D5 = 2.525 epoch-eq and §4.9 scales by 1.075 (arithmetic checks: 0.41–1.24 h × 1.075 = 0.44–1.33); §9.2's 6c row still says "≈ 4 LIDC epoch-equivalents … 1.72 / 1.70 / 1.45" |
| 6d B1 AMEND 2/3 | **Partly** | AMEND 1–3 recorded; `lab/B1/prereg.py` has **AMEND 4** (Qwen3.5 at 34 s/step in bf16 joins the 150-step tier), so D1's "clause (2) compares Qwen3.5 with LFM at 36 steps only" is stale; clause (2) can use the 150-step comparison |

## New MUST-FIX (blocking approval)

1. **G-topology (b) is a C1 rung gate; a hole-recall miss would climb the ladder (§5.1, §5.2 C1, §6.3, §8.2).**
   Problem: C1 = "first rung passing every §5.1 gate"; §5.1 lists G-topology (b) hole recall ≥ 0.30. §6.3's trigger
   list names only (a), and §8.2 says "C2 failure does not block Phase 3". So a rung with hole recall 0.25 fails C1
   with no diagnostic, the ladder climbs to Qwen3.5 ($1.1–2.7 + C2), and the C2 fail branch ("holes left to GRPO")
   never fires. Fix: make (b) report-only for C1 (the rung passes on (a)); (b) is decided by C2 only. Keep (b) in R1's
   report so R-ctrl's trigger can be evaluated early.
2. **G-probe is uncalibrated at the production GRPO LR and its failure is the costliest false negative (§5.1, D12, §8.2).**
   Problem: `configs.grpo_config` sets LR 1e-6 (full FT), constant, 20 steps × 32 prompts; the pass rule is a paired
   CI > 0 on 128 items × 8 samples. Nothing says whether 20 steps at 1e-6 moves R by more than the probe's sampling
   noise. The chain "selected → next branch → D12(c) repair → rung fails → ladder climbs" turns probe noise into
   $1–5 of Vast spend and discards a checkpoint that passed every other gate. Fix: (a) **noise floor at zero cost**:
   at step 0 sample the 128-item set twice (different sampler seeds) and report the paired CI; pre-register that a
   G-probe "fail" counts only if the step-20 CI is also outside the null CI; (b) pre-register the probe LR
   (the production LoRA 1e-5, or 5× full-FT) and record it as the probe's, not Phase 3's; (c) consequence: G-probe
   fail after the repair → **flag and hand to Phase 3 with the warning** if G-signal passes; climb the ladder only if
   G-signal also fails. Add the null and the LR to `prereg.py` (M8).
3. **C2's A0 arm has no decision role; two of three C2 runs buy a descriptive guard (§5.2, §6.1 R3–R5, §4.9).**
   Problem: A0 scores 0 hole recall by construction (`topo_probe_v3e_seed8.json`), the paired difference is dropped
   as uninformative, and the only A0 output is the "guard reported" (A0′ shifted IoU_tc ≥ A0 − 0.01, no decision).
   A0 seeds 1 and 2 cost 0.8–2.5 GPU-h, while the plan's own caveat says two A0′ seeds give 54 % coverage of the seed
   mean at σ_v = 0.3. Fix, same budget: train **A0 seed 1 only** (guard, pre-registered as a threshold: paired
   A0′ − A0 shifted-val median IoU_tc CI lower bound > −0.02, else "hole-awareness costs geometry" is reported) and
   use the freed run for **A0′ seed 3**, making the "third seed if discordant" rule the default. C2 rule with three
   seeds: pooled 3-seed CI lower bound > 0.30 and pooled ≥ 0.50; ≥ 2 of 3 seeds ≥ 0.50 and none < 0.30; per-bin,
   control and legality clauses unchanged. If the owner keeps 2 + 2, state in §5.2 that A0 is descriptive and why.
4. **A SmolVLM C1 would hand Phase 3 a model without a vLLM path (D1, §6.1, §5.2 C1).**
   Problem: B1 rule (1) can keep SmolVLM-256M and §6.1 screens it first; with "stop at the first passing rung",
   C1 = SmolVLM is possible. TRL #6294 (SmolVLM + vLLM prompt corruption, 855 vs 38 tokens) is open, PRs #6718/#6719
   unmerged, `grpo_config` already forces `use_vllm=False` for it, and A1 records its bf16 exclusion from TRL's GRPO
   CI. Phase 3's cost projection then rests on HF generation (Mac: 13–17 s generation per step). Fix: SmolVLM is a
   control, never C1; or add a pre-registered Phase 3 viability clause to C1: G-parity passes **or** the HF-path
   G-cost projection ≤ $150, measured in R2.
5. **CT window jitter needs int16 HU crops; the bundle and M1 inherit 8-bit PNG crops (D10, §4.4, M1, §6.0).**
   Problem: `data/prepared` images are lung-windowed uint8 (`targets/frame.py` LUNG_WINDOW (−1000, 400) = L −300,
   W 1400); B3's context crops are PNG bytes (`dihedral.py` `_png_bytes`, 17.6 MB), and §4.4 ships exactly those.
   Raw HU exists only in `data/store` (int16, `store.py` B3-7). A4 §4.2.3 says "this needs `frame` to take the
   window"; M1 adds the parameter but no HU crop. On the box the jitter cannot be applied. Fix: M1 writes one int16 HU
   context crop per LIDC item (npz) from the store into `data/views`, M1 test: `window=(−1000, 400)` reproduces the
   prepared PNG byte-for-byte; §4.4 bundle size recomputed (≈ 7,341 × crop² × 2 bytes, compressed). Also cite
   2510.08116 as "random windowing for CT (liver tumours)", not "standard CT augmentation".
6. **Evidence corrections (text only).**
   (a) vLLM #34948 (Qwen3.5 GDN illegal memory access) was fixed in vLLM 0.27.1 (FlashInfer PR #4553, vLLM #52553):
   D13/§8.1 cite it as a live 0.32 risk; the live risks are the CPU parity gap (B5) and `fla`/`causal-conv1d` wheels
   for torch 2.14.1 in `/opt/venv`. (b) TRL #6294's numbers are 855 vs 38 prompt tokens, not 89 → 155 (D13).
   (c) A1 §1.5: −2.7 IoU is G = 100, −1 and −1.5 % legality is G = 250 (D2). (d) Qwen3.5-0.8B's card does not mention
   `bbox_2d` / 0–1000 (that is Qwen3-VL); mark "UNVERIFIED on the card; B1 T1 tests it" (D1). (e) SimpleSeg: 5e-5 is
   its SFT LR (cosine), RL LR 2e-6, GRPO gain +9.7–10.5 cIoU (D5, §8.2 "+10–12"). (f) A5's 1,536-token runaway is §0.6,
   not §7 (G-format row). (g) The size rule max(5 px, 0.5 %·|U|) is in the reward code, not in any note (D20). (h)
   §9.2 rows 6c, SHOULD-FIX 1 and 11 describe revision 1 (4 epoch-eq; "E3 1-seed sanity, S7 smoke"; "S7 … 1 h"):
   rewrite to revision 3. (i) D1: record AMEND 4. (j) §8.3 Q3 "$25 tranches up to the $60 cap" → $15 stop.

## SHOULD-FIX

1. **Threshold mismatch on controls.** G-topology (a) ≥ 90 % correct (1, 0) on simple controls vs C2 (4) false
   non-simple ≤ 5 %: align (a) to ≥ 95 % (σ = 0.5 synthetic model scores 0 false non-simple, so 95 % is not tight).
2. **§6.3 depth when no branch is eligible.** "Trained to the failing run's selected depth" is undefined when
   G-format fails at every branch; say "the deepest branch".
3. **R-ctrl is one seed against a two-seed rule.** C2 clause (2) cannot be met by one Qwen3-VL-2B seed; define
   R-ctrl pass as clauses (1), (3), (4), (5) on one seed, labelled exploratory, before "holes need ≥ 2B" is written.
4. **B2 has not started** (no `lab/B2/`). Add to D12: "B2 not run before V1 → the rule stands as written and the
   D12(c) repair is the next-earlier branch (WiSE-FT untested)". B1's rule is sound as pre-registered; D1 already
   limits what 36 / 150 LoRA steps can say about a full-FT rung.
5. **Vision-tower LR.** A2 §4.2's evidence is "unfreeze"; the same-LR choice is the plan's. Either adopt 0.1× for the
   tower (the conservative arm of H-A2-6, matching Qwen's and Axolotl's practice) or add "vision LR 0.1×" as the
   first option of §6.3 rule 2 (same cost as "dihedral off", smaller change to the training distribution).
6. **Free attribution check.** Score the selected checkpoint on the identity view and on one rotated val view
   (inference only, ≈ 0.05 h): if rotated ≈ identity, dihedral training transferred; if not, the capacity spent on
   8 orientations is visible without a run. Report, no gate.
7. **Budget table.** Add the MFU 0.05 row for scenario A (≈ 6–12 GPU-h, ≈ $4–8), note the eval block is unverified,
   and that ×4 oversampling adds ≈ 5 % rows; quote §4.8's A100 row (1.09 / 1.19, n = 7) with its filter, since the
   snapshot's unfiltered A100 SXM4 on-demand is 0.477 / 0.923 (n = 14).
8. **Rescaled val coverage.** Report the share of LIDC items with no admissible neighbouring divisor (50 ↔ 100 ↔ 125
   halves or doubles the apparent size); if > 30 %, G-geometry on the rescaled view is a weaker check than §5.0 implies.
9. **A2 §4.1 caveat dropped:** LoRA is less tolerant of large batches (Qwen3-VL-2B at effective batch 32 is fine;
   record it).
10. **A4 S5 LIDC prior** is 0.588 / 65.3 %, i.e. fails the centred LIDC gate narrowly; §5.1 G-prior's evidence line
    should say so rather than cite RIGA only.

## Strengths (new since R1)

- The hole probe is now the right instrument: real hole geometry drawn jointly from train records, size-rule floor,
  stratified bins, powered pooled metric (±3.1), per-bin floor reliability (P < 0.30 at true 0.4 = 0.010), honest
  P5.2 failure with iteration stopped before tuning to the test.
- Revision 3 keeps every removed experiment with a concrete trigger (§6.4) and folds the cheap measurements into R0;
  the critical path is 3 + 3 runs on one box with Mac-side baselines, and the $10 / $15 caps are enforced by G1–G10.
- Environment claims re-verified: vllm 0.32.0 pins torch 2.13.0 / numba 0.65.0 / transformers < 5.20; trl 1.16.0
  declares vllm ≤ 0.31.0 and PR #7799 is still open, so the "expect the warning" note is correct.
- Rank inversion (2606.18487: pass@10 0.806 → 0.481) and Quagmires are now cited correctly, with the right direction.
- The shifted/rescaled views, G-prior and the sequencing fix close the Phase 1 leak (I43) end to end.

## What remains for APPROVE

Apply MUST-FIX 1–6 in the plan text (1, 3, 4, 6 are wording and rule changes; 2 adds a null re-sample and a probe
LR to `prereg.py`; 5 is one M1 code item plus a test) and the ROADMAP line if MUST-FIX 3 changes the seed count.
No new lab run is needed. With those, every Vast run in §6.1 feeds a decision, no gate can send the ladder up on
noise or on a C2-scoped quantity, and the bundle can actually be executed on the box.
