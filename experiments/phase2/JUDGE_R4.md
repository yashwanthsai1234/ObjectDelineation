# Phase 2 plan: judge round 4

Confirmation review of `experiments/phase2/PLAN.md` revision 7 against `JUDGE_R3.md`, read in full (lines 1–911), with the working-tree `ROADMAP.md` diff. Read-only; no file changed, no git action.

## Overall verdict: APPROVE

MUST-FIX 1 and SHOULD-FIX 1 and 5 (the three conditions I set in round 3) are applied in the plan body, in every section that carries the item, not only in §9.2c. The remaining SHOULD-FIXes are applied too. The new revision-7 content (C2 outcome taxonomy, R-x8) is sound, costed, approval-gated and inside the owner's cost directive. The §4.9 scenario arithmetic holds for A and B; scenario C is understated by ≈ 0.5 GPU-h / ≈ $0.4–0.9 at the high end, which does not change its conclusion (correct figures below). I found no new inconsistency that would cause a wrong decision, wasted money or an unrunnable step.

## 1. Round-3 items: verified in the body

| Item | Verdict | Where verified (body, not §9.2c) |
|---|---|---|
| MUST-FIX 1 G-prior (a) uncalibrated gate | fixed, option (i) | §5.1 G-prior: "(a) … report-only for C1", "per-item IoU_tc = the K = 8 item mean with illegal = 0", reported beside MedSAM-oracle and SAM 2.1-tiny shares, gating only with a threshold ≤ SAM-tiny − 5 pts in `prereg.py` before V1; §5.2 C1 "(G-prior (a) is report-only …)"; D17 computes the two reference shares in the Mac run (K = 1); D18 "report-only since revision 7"; §6.3 trigger lists G-prior (b) only; §8.2 "passes G-prior (b), G-signal and G-probe"; M8 `prereg.py` "G-prior (a) threshold only if the owner gates it before V1". No section still gates on (a). |
| SHOULD-FIX 1 "both seeds" | fixed | §5.0: "the per-seed clauses of its decision rule (C2: ≥ 2 of 3 seeds ≥ 0.50, none < 0.30; §5.2), with `prereg.py` (M8) authoritative where it and this text differ". The only remaining "both seeds" / "2 seeds" strings are in the §9.1 revision-2/3 log entries and the §9.2 R1 MUST-FIX 3 row, which are history. |
| SHOULD-FIX 2 in-training abort | fixed, consistent | D7 ("both (a) … and (b) …, one of the two alone is warn-only"), §3.2 ①, §3.3 R1 line, §5.1 G-vision "Measured" column, §6.1 R1, §4.9 vision-repair row, §8.1 risk row, M8 tests ("aborts … only when ratio < 0.5 and TV < 0.002, and only warns when one is"). |
| SHOULD-FIX 3 Spearman as switch | fixed | D12 ("can only retain sd_geom or flag disagreement; switching to r_std needs the branch-level comparison … Phase 3 planning decision"), §6.2 row, §9.3 D12 row, §9.4. The §9.1 revision-6 entry still says "confirms or reverses" (log; superseded by the revision-7 entry). |
| SHOULD-FIX 4 C2 taxonomy | fixed, consistent | §5.2 pass / partial / fail with exactly one conditional each and the exclusions ("R-ctrl never follows a partial and ×8 never follows a fail"; ×8 not run when clause (4) failed); §6.1 R-x8 and R-ctrl rows agree; §6.4 E1; §3.3; §4.9 two rows; §8.2 "After C2"; §1. |
| SHOULD-FIX 5 C2 depth after a repair | fixed | §5.2 ("the D12(c) repair branch if it was probed, not necessarily the first D12 pick"), §6.1 R3–R5, §3.3, §4.9 R3–R5 row. |
| SHOULD-FIX 6 stale names | fixed | M10 "from R0 / V1", Depends V1; §5.3 hash row "R0 (V1)"; §9.3 GRPO row "colocate measured at Phase 3 start (§6.4)"; D13 "K = 64 … at Phase 3 entry (§6.4)"; D1 and §8.1 / §8.2 "Qwen3-VL-4B or InternVL3.5-1B as extra rungs after the 2B rung"; §5.1 G-precision "zero-shot, then the capacity run only on a zero-shot failure". The two remaining "S7" strings are A4's section S7 (D17) and the §6.4 / budget-history rows, which are correct. |
| SHOULD-FIX 7 ROADMAP M10 scope | fixed | M10 lists ROADMAP lines 165–176 item by item with the superseding decision. The working-tree ROADMAP diff (C2 wording, A0-multi) is consistent with §1 and D2. |
| SHOULD-FIX 8 Qwen3.5 row | fixed | §4.9 row 1.7–3.9 h / $1.1–2.6; D1 rejected-column uses the same; side-by-side paragraph 5.4–11.4 GPU-h / $3.4–7.5 / ≈ 20 GPU-h at MFU 0.05 all recompute from the rows. |

## 2. Revision-7 content

**C2 outcome taxonomy (§5.2).** Exhaustive (pooled < 0.30 → fail; [0.30, 0.50) or ≥ 0.50 with any of clauses (1)–(5) failing → partial; all five → pass), each with one named conditional, both approval-gated, and the clause-(4) exclusion for ×8 is the right direction (more hole rows raise false non-simple; §6.3 rule 1 drops to ×1 for the same reason).

**R-x8.** Costed in §4.9 (0.5–1.4 GPU-h, $0.3–0.9, 4090), listed in §6.1, §6.4 E1, §3.3, §8.2, §1, owner-approval gated, exploratory (one seed), run only on a partial. Against the directive: it is a triggered single-component check that fires only when its failure mode appears (the §6.4 pattern since revision 3), on the smallest passing model, and scenario A plus R-x8 stays at ≤ 9.1 GPU-h / ≤ $6.0 at the median, under the $10 alarm. The $10 / $15 caps are unchanged (G10, §4.9, §5.3, §8.2). Critical-path cost unchanged.

**§4.9 arithmetic recomputed from the rows.**
- Scenario A: 0.2 + 0.9 + (0.76–1.70) + (0.2–0.5) + (1.6–4.2) + 0.2 = 3.86–7.7 GPU-h; $0.06 + 0.59 + (0.50–1.12) + (0.13–0.33) + (1.06–2.77) + 0.13 = $2.47–5.00; ceilings $7.56. Plan: 3.9–7.7 / $2.5–5.0 / ≤ $7.6. Correct. MFU 0.05 row: R1 3.1 h (< 3.4 h stop), total 13.0 h, $8.6. Correct.
- Scenario B (build-up: A minus LFM's R2, which is not run when LFM fails its gates, + R-diag + Qwen3.5 R1+R2 + C2 extra): (3.86 − 0.2) + 0.5 + 1.72 + 1.29 = 7.2; (7.7 − 0.5) + 1.4 + 3.9 + 3.9 = 16.4 GPU-h; $4.66–10.74. Plan: 7.2–16.3 / $4.7–10.7. Correct to rounding (16.3 vs 16.4).
- Scenario C (build-up: V0 0.2 + R0 0.9 + LFM R1 0.76–1.70 + R-diag 0.5–1.4 + Qwen3.5 R1 only 1.52–3.40 (its R2 is not run when it fails) + 2B R1+R2 1.5–3.3 on the PRO 6000 + C2 on 2B = 3 × the R-ctrl per-run cost 0.8–2.3 = 2.4–6.9 + idle 0.2): **8.0–18.0 GPU-h, $8.9–21.4** ($0.06 + 0.59 + 0.50–1.12 + 0.33–0.92 + 1.00–2.24 + 2.4–5.3 + 3.9–11.1 + 0.13). If the ≈ 0.3 h probe-eval block for the three 2B C2 runs is added (the R-ctrl row is training only): ≈ 8.3–18.3 GPU-h, ≈ $9.4–21.9. The plan's 7.9–17.5 / $8.9–21.0 understates the high end by ≈ 0.5 GPU-h / $0.4–0.9. The stated conclusion (can pass the $15 stop; 2B rung starts only with an approved tranche) is unchanged.

## 3. New inconsistencies or stale references

None that affect a decision, a cost or a runnable step. Items found are wording (listed below).

## Optional (polish; no GPU cost, no decision change)

1. §4.9 scenario C: replace with 8.0–18.0 GPU-h / $8.9–21.4 and state the build-up (LFM R1 → R-diag → Qwen3.5 R1 only → 2B R1+R2 → 3 × R-ctrl-cost C2 runs), as done for A; mention that the 2B C2 probe evals (≈ 0.3 h on the PRO 6000) are not in the R-ctrl row.
2. §6.1 R-x8: give it the same exploratory pass rule as R-ctrl (clauses (1) on its seed, (3), (4), (5)) and say whether it uses C2's topology-probe draw or a disjoint one (§9.1 revision 2 says E1 uses a disjoint draw).
3. §5.0 fourth-seed clause (added on a "decided by a single seed's side of 0.50" result) vs §5.2 "each outcome has at most one conditional run": say the fourth seed, if approved, is part of the C2 verdict (re-judged with the per-seed clause restated for 4 seeds) and precedes R-x8 / R-ctrl, so two conditionals cannot fire on one result. Also §6.2 determinism row "any third C2 seed" → "fourth".
4. §6.1 R-ctrl "(owner approval if past the $15 stop)" vs §5.2 "owner approval required" for every C2 conditional: align to §5.2.
5. §3.3 partial branch label "(pooled [0.30, 0.50))" omits the "≥ 0.50 with a failed clause" case; header "revisions 3–6" and §4.9 header "revisions 3–5" are pre-revision-7.
6. Log hygiene: §9.1 revision-5 entry "+1.7–4.2 GPU-h" (now 1.7–3.9) and revision-6 "confirms or reverses"; the revision-4 entry is placed before revision 3.

## What was checked and how

Full read of PLAN.md revision 7; grep of the body for every item's wording across sections (G-prior (a), abort conditions, R-x8 / R-ctrl / ×8, owner approval, S0 / S7, "both seeds", the 4.2 h figure); scenario totals recomputed row by row; `git diff ROADMAP.md` read for the C2 / A0 wording. No file was edited and no git state was changed.
