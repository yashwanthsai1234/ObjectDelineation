> Frozen record (written 2026-10-20). Versions, budgets and paths were true then. Current facts: docs/README.md.

# Phase 2 plan: judge round 3

Independent review of `experiments/phase2/PLAN.md` revision 6 against `JUDGE_R1.md`, `JUDGE_R2.md`, `RESEARCH.md`, `lab/B1/REPORT.md`, `lab/B2/REPORT.md` (and their `prereg.py`, plus `analysis.json` / `posthoc.json` in this session's scratchpad), `lab/B3/REPORT.md` Addendum A, `lab/B4`, `lab/B5`, notes A1–A6 (spot checks), `ROADMAP.md` (working-tree diff), `experiments/phase1/RESULTS.md` §8, and the source the thresholds rest on (`src/objdel/train/configs.py`, `src/objdel/eval/gates.py`, `src/objdel/reward/spec.py`, `lab/B1/common.py`). Read-only; no file changed.

## Overall verdict: APPROVE WITH CHANGES (one MUST-FIX, a threshold calibration; the rest is wording)

The plan now does what the owner's brief and the cost directive ask: one evidence-chosen bundle, Vast only for C1 / G-probe / C2, stop at the first passing rung, $10 alarm / $15 stop, every lab verdict carried into the text with its honest label. Every R1 and R2 item is resolved in the text, not only in §9.2. The one remaining defect that can produce a wrong C1 decision is G-prior (a), a gate on the centred view whose 60 % threshold has no calibration and can fail a model that passes every image-use check.

## 1. R1 and R2 items: resolved in the text?

| Item | Status | Where verified |
|---|---|---|
| R1 MF1 C2 hole-only, calibrated probe | resolved | §1 C2, D20 (120/bin, bins 44/36/20, P5.2 "failed narrowly"), §5.2; matches B3 Addendum A1–A6 |
| R1 MF2 "GRPO can improve" KPI | resolved | §5.1 G-probe with PASS / FAIL / INCONCLUSIVE against a step-0 null |
| R1 MF3 seed pooling / primary metric | resolved, one stale sentence | §5.0 (per-item mean over seeds and K, nodule / within-bin bootstrap); "both seeds' point estimates" survives from the 2-seed era (SHOULD-FIX 1) |
| R1 MF4 baselines on the gate view | resolved | D17 (every row labelled; nnU-Net Phase 9) |
| R1 MF5 trainer vs vLLM on one GPU | resolved | D13, §3.2 ①–④, M8 smoke |
| R1 MF6a–d evidence | resolved | D5 (Quagmires direction), D20 (B3 P2), §4.9 (2.5 epoch-eq, ×1.13), D1 (AMEND 1–5) |
| R1 SF1–12 | resolved (SF1 by the directive; SF8 by R-B; SF12: `eval/gates.py` still carries `topo_any_unanimous`, M4 replaces it) | §6, D12, M4 |
| R2 MF1 G-topology (b) as C1 gate | resolved | §5.1 "report-only for C1", §5.2 C1, §6.3 trigger list names only (a) |
| R2 MF2 G-probe uncalibrated | resolved | §5.1: step-0 double sample (h₀), probe LR 5e-6 / 5e-5 (production 1e-6 / 1e-5 confirmed in `configs.py`), no ladder climb on G-probe alone; consistent in §3.3, §6.1, §8.2 |
| R2 MF3 A0 arm buys no decision | resolved | §5.2 / §6.1: 3 A0′ seeds + 1 A0 seed, 3-seed rule, A0 guard report-only; §4.9 row and §3.3 agree |
| R2 MF4 SmolVLM as C1 | resolved | D1 "control only, never screened on Vast, never C1" |
| R2 MF5 CT jitter needs HU | resolved | D10, §4.4 (int16 npz, 0.3–0.5 GB est.), M1 byte-identity test |
| R2 MF6a–j | resolved | D13 / §8.1 (#34948 fixed in 0.27.1), D13 (#6294 855 vs 38), D2 (G = 100 / 250), D1 + §9.3 (`bbox_2d` undocumented, observed in B1), D5 (SimpleSeg SFT 5e-5 / RL 2e-6), G-format (A5 §0.6), D20 (size rule in reward code), §9.2 rows rewritten, D1 AMEND 4, §8.3 Q3 $15 |
| R2 SF1–10 | all resolved | §5.1 (a) 95 %; §6.3 deepest branch; R-ctrl exploratory; B2 done; D7 0.1× / frozen; R1 attribution check; §4.9 MFU 0.05 row, A100 filter; §5.0 rescaled coverage; D6 LoRA batch caveat; G-prior cites 0.588 / 65.3 % |

## 2. B1 and B2 fidelity; the new decisions

**B1 is represented faithfully.** Checked against `lab/B1/REPORT.md` and `analysis.json`: H-A1-1 and H-A1-2 INCONCLUSIVE (not accept), the SmolVLM drop rule "not triggered as written", legality 82–83 % vs 42–59 % at 36 steps (paired shifted LFM − Qwen3.5 −0.268 [−0.315, −0.219]), 83 % vs 71 % shifted at 150 steps (−0.125 [−0.178, −0.076]), LFM centred median +0.065 [0.026, 0.104], tokens/vertex 6.5 vs 9.5, the vision-collapse numbers (0.297 → 0.019 → 0.002; 0.171 → 0.003; LFM 0.471 → 0.321 / 0.373; TV ≤ 0.0007), and the caveats (single-seed Qwen3.5, AMEND 5 post-hoc, looser shift draw in `common._shift`). The G-vision thresholds are derived correctly from those numbers: ratio 0.5 sits 1.4× below the lowest intact ratio (0.68) and 8× above the highest collapsed (0.064); TV 0.002 is 2× below the lowest intact (0.004) and ~3× above the highest collapsed (0.0007). The plan says they rest on 4–8 post-hoc items and makes R0 re-measure the base references; that is the right label.

**B2 is represented faithfully.** H-B2-1 INCONCLUSIVE and H-B2-2 REJECT match `analysis.json` (`D12_rule: R-B`, `D12_repair: next-earlier branch`); the shifted falls (LFM r_std −0.028 [−0.044, −0.014], headroom −0.047 [−0.070, −0.025]), the interpolation losses (α 0.75: −0.086 / −0.070 legal; α 0.5: 0.375 legal, 10.5 % duplicates), and "legality part 55–92 % of group variance" (geom_share 0.08–0.45 in `posthoc.json`) all check out.

**The post-hoc sd_geom veto is handled honestly.** D12 labels it post-hoc, explains why the pre-registered r_std was not kept (its legality term rewards a less legal branch; p(1−p) moves 0.09 → 0.02 across the eligible set, enough to clear 15 % alone), reports both statistics per branch with a "statistics disagree" flag, writes the choice into `prereg.py` before R1, and pre-registers an R2 check. The one weakness is the check itself (SHOULD-FIX 3): a per-item Spearman within one checkpoint confounds "spread predicts gain" with "items at p ≈ 0.5 legality have the most room", which mechanically favours r_std; switching Phase 3's veto statistic on that would be a decision on an artefact. It affects later rungs only.

**Consistency of D1, D6/D7, G-vision, mandatory crop jitter, R-B.** They agree with each other and with the evidence: smallest-first with a 0.2 GPU-h R0 Qwen3.5 arm (no pick rule) is the correct reading of B1 §0.5 under the directive; LoRA frozen tower (the tested setting) vs full-FT 0.1× (A2's own proposal, checked by G-vision and the D7 repair) is stated as a choice, not as B1 evidence; crop jitter mandatory follows from B1 §2.5 and E3 is closed; R-B is B2's own rule with the disclosed statistic change. B2 §4's note ("if the plan prefers only pre-registered statistics, use r_std") is answered in D12.

## 3. Internal consistency

Gates have fixed thresholds (§5.1 every row; `prereg.py` M8). §6.1 runs, §3.3 diagram, §4.9 rows and §4.7 steps agree (R0 0.9 h; R1 0.76–1.70; R2 0.2–0.5; R3–R5 1.6–4.2; R-diag 0.5–1.4; R-ctrl 0.8–2.3). I recomputed §4.9: scenario A 3.9–7.7 GPU-h / $2.5–5.0 / ≤ $7.6 at ceilings, scenario B 7.2–16.6 GPU-h / $4.7–10.9, R3–R5 per-run 0.43–1.30 h, the MFU 0.05 row (R1 ≈ 3.1 h < the 3.4 h stop), and the 2.5 epoch-eq scaling: all arithmetic holds. §7 milestones and the graph agree ("V1 needs M1–M5, M7, M8, M11"); M6 removed.

Stale references that remain (all wording): M10 "updated from S0" and "Depends: S0"; §5.3 Hash row "Measured: S0"; §9.3 GRPO row "colocate measured in S7 (§4.7)" (S7 moved to Phase 3, §6.4); §5.0 "both seeds' point estimates" (3-seed rule now in §5.2); D1 "InternVL3.5-1B … only if LFM and Qwen3.5 both fail" vs the 2B third rung, and §8.1 "Qwen3-VL-4B rung" vs §8.2 "4B or InternVL"; D13 "K … 64 for TV/retention … I35" while §6.4 moves K = 64 to Phase 3; §5.1 G-precision definition column says "capacity run" though the protocol is zero-shot first.

## 4. Risks: wasted Vast money, or passing without the image; the cost directive

**Image use is well guarded.** G-geometry at median 0.60 on the shifted and rescaled views cannot be met by a location/size prior (prior 0.00 on shifted LIDC), G-prior (b) catches the centring shortcut B1 exposed (every B1 arm was below the prior on shifted crops), the rescaled view removes the size-from-S cue, and G-vision catches tower collapse. I found no residual leak in the views or the prompt.

**Money risks, in order:**
1. **G-prior (a) can fail C1 on an uncalibrated threshold** (MUST-FIX 1). The centred image-blind prior scores median 0.588 with 65 % of items ≥ 0.5, i.e. near the gate floor; the per-item statistic is a mean over K = 8 with illegal = 0, so a model at legal 90 % and shifted median 0.60–0.65 has a centred item-mean close to the prior's and may win fewer than 60 % of items while passing (b), G-geometry on both views and G-vision. The consequence is R-diag (0.3–0.9 $) then a ladder climb (+$2–5 and a wrong "LFM not viable"). The 60 % comes from A4 §2.5 proposal 2 with no reference model's win share behind it.
2. In-training vision abort on ratio (a) alone with a full-FT tower (SHOULD-FIX 2): thresholds come from LoRA collapse vs intact; a legitimate full-FT feature shift could trip 0.5 and discard part of an R1 for a frozen rerun. Bounded by rolling-checkpoint cadence and ≤ one R1.
3. Overlapping conditionals after a C2 fail (SHOULD-FIX 4): pooled recall in [0.30, 0.50) triggers both R-ctrl (§6.1) and the ×8 run (§6.4 E1); both need approval, but the plan should say which fires.
4. G-probe cost (0.2–0.5 h) and the eval block (0.3 h) are assumptions; the R0 stop rule and the $15 stop bound them.

**Cost directive:** respected. One bundle, no per-component A/Bs (each kept with a trigger), stop at the first passing rung, no unconditional 2B control, C2 reuses R1, baselines on the Mac, the only added run since revision 3 is the 0.2 GPU-h R0 Qwen3.5 arm with an explicit no-pick rule, and revision 6 adds no GPU cost.

## 5. KPI scores (1–5)

| Area | Score | One line |
|---|---|---|
| Research depth | 5 | Six primary-source notes, corrected across three rounds; every disagreement between notes resolved in §9.3 |
| Evidence / testing rigour | 4 | B1–B5 pre-registered with honest INCONCLUSIVE / REJECT verdicts and post-hoc labels; one gate threshold (G-prior (a)) is still uncalibrated and G-vision rests on 4–8-item diagnostics (acknowledged) |
| Architecture clarity | 5 | §3.1–3.3 drawn, process ownership of the GPU explicit, diagrams agree with §6.1 and §4.9 |
| Algorithmic modernity | 4 | WSD branches, rank-inversion-aware selection with an exact variance split, certified offline augmentation, vLLM 0.32 bridge; NTL / coordinate tokens correctly parked behind triggers |
| Engineering / test-driven | 4 | Tests-first milestones, M8 CPU smoke of the sequencing, hash gates; a handful of stale names (S0, S7) and ROADMAP's Phase 2 text still stale pending M10 |
| Cost discipline | 5 | Directive met in design and in arithmetic; MFU 0.05 row, R0 stop rule, $10 / $15 caps |
| Vast ops | 5 | Five cards with roles, ceilings and fallbacks from a dated snapshot; scripts, G1–G10, tested governance, destroy-never-stop |

## MUST-FIX (blocking)

1. **G-prior (a) is an uncalibrated C1 gate on the centred view (§5.1 G-prior, D18).** Fix, zero GPU cost: compute, in the Mac baseline run (D17, already planned before V1), the per-item win share over the centred prior for MedSAM-oracle and SAM 2.1-tiny, and either (i) make (a) **report-only** (G-prior (b) on shifted + rescaled, G-geometry on both views and G-vision already decide image use), or (ii) keep (a) as a gate with its threshold pre-registered at or below the SAM 2.1-tiny win share minus 5 pts, recorded in `prereg.py` (M8) with the measured values in the "Evidence for the threshold" column. Either way, state that (a) is computed on the K = 8 item mean with illegal = 0, as §5.0 implies.

## SHOULD-FIX

1. **§5.0 "both seeds' point estimates" vs the §5.2 3-seed rule.** Replace with "the per-seed clauses of the decision rule (C2: ≥ 2 of 3 ≥ 0.50, none < 0.30)" and state that `prereg.py` is authoritative.
2. **D7 / §5.1 G-vision in-training abort (full FT).** Require both (a) ratio < 0.5 and (b) TV < 0.002 at the rolling checkpoint before aborting (collapse in B1 hit both; (b) costs 8 items × 24 tokens on the GPU), or make the in-training check warn-only and leave the decision to `select_ckpt`.
3. **D12 / §6.2 R2 Spearman as the switch for the veto statistic.** Record that the per-item correlation is confounded by legality headroom, so it can only *retain* sd_geom or flag disagreement; switching to r_std needs the branch-level comparison (D12(a)/(b), when a second branch is probed) and is a Phase 3 planning decision, not an automatic rule.
4. **§5.2 / §6.1 / §6.4: C2 fail taxonomy.** Define three outcomes (pass; partial, pooled in [0.30, 0.50); fail, < 0.30) and name the single conditional that each may trigger (R-ctrl for < 0.50 only when no ×8 run is approved, or ×8 first), so two conditional runs cannot both fire on one result.
5. **§6.1 R3–R5 depth after a D12(c) repair.** Say that A0′ seed 1 and the C2 depth are those of the checkpoint that received the final G-probe verdict (the repair branch if it was probed), not necessarily the first D12 pick.
6. **Stale names.** M10 and §5.3 "S0" → "R0 / V1"; §9.3 "colocate measured in S7 (§4.7)" → "Phase 3 start (§6.4)"; D13 K = 64 clause → "Phase 3 entry (§6.4)"; D1 InternVL and §8.1 "Qwen3-VL-4B" aligned with §8.2's "4B or InternVL3.5-1B as extra rungs after the 2B rung"; §5.1 G-precision definition "zero-shot, then capacity run".
7. **ROADMAP Phase 2 section (M10 scope).** Lines 165–176 still list SmolVLM-256M / 500M / Qwen3-VL-2B, "a random expert's polygon per epoch" and "flips, 90° rotations" for both datasets; M10 should list these lines explicitly so the handoff doc and the plan do not disagree at Phase 2 close. (The C2 wording change in the working-tree diff is correct.)
8. **§4.9 conditional Qwen3.5 row.** 2 × (0.76–1.70) + (0.2–0.5) = 1.7–3.9 h, not 4.2; harmless (conservative) but worth aligning.

## Strengths since R2

- B1 and B2 were run as pre-registered, reported with their failures (no arm learns the format on the Mac; interpolation rejected; H-B2-1 inconclusive), and the plan took the conservative reading each time (default ladder kept, vision frozen under LoRA, deepest-branch rule) rather than the optimistic one.
- The vision-collapse finding is turned into a gate with thresholds that sit between the two measured populations, plus an in-training check and a bounded repair.
- The G-probe is now decidable before the run: a null at zero cost, a probe LR with a stated rationale, and a failure path that cannot spend money on noise.
- Revision 6 added measurement, not compute: the signal panel, depth curve and Spearman check all reuse samples that `select_ckpt` and `g_probe` already draw.

## What would make me APPROVE

MUST-FIX 1 applied (a report-only (a), or a threshold calibrated on the Mac baselines and written into `prereg.py`), plus SHOULD-FIX 1 and 5 so that C2 is unambiguous before R3 starts. No new lab run is needed.
