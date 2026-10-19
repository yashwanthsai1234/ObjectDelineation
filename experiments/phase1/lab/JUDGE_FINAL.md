# Judge final review of PLAN.md (2026-10-19)

Read: PLAN.md, JUDGE_VERDICT.md, ITERATION3/4.md, INTEGRATION.md, novelty_check.md, t_data/lidc_full/{parity_all,fixture_summary}.json. Code checked: `t_stack/{prod,r2_trl_equiv}.py`, `t_reward/scripts/{inflation_incentive,honest,r4_scale}.py`, raw `iter4_{r4,inflation_incentive,adversaries_meaniou}.txt`, `tests/test_reward_v2.py::test_valid_always_beats_invalid_for_all_reachable_states`, `rewardcore/{reward.area_guard,advantages.reward_v4}`.

## Verdict: APPROVE-WITH-EDITS

The reward is now tuned on 251 real disagreement items with realistic honest rollouts, the K9/R4/K10/K3 evidence is real and consistent with the raw files, and the plan is honest about the hedges it cannot beat. Two gate rows overstate what was run; fix them and the plan can go to the user.

## R1–R10

| R | Status | Honest in plan? |
|---|---|---|
| R1 | Addressed on the fixture set (6,338 targets; fastgeom bit-exact; canonicalization invariance). Not yet "every train target" (full LIDC ≈ 35k masks). | Number honest; "✅" slightly early. M5 runs the full set. |
| R2 | **Partial.** `prod.py` calls `reward_v2`; `r2_trl_equiv.py` compares to `advantages_v2`; rewardcore is the iteration-2 snapshot. v4 (g_area term, no batch z, debias, clip) has never run inside TRL. | **No.** D7 and §7 present 8.7e-8 as if for the frozen reward. |
| R3 | Addressed, then split (see below). | Yes, with one overclaim (local gradient, below). |
| R4 | Addressed: 251 items, both honest modes, CI half-width 0.051/0.056. | Yes. D5 cites the v3 number (+0.009), not v4. |
| R5 | Addressed (per-stratum K8, S_m test, single-rater exclusion; k=2 prior cross-fitted per protocol). | Yes. |
| R6 | Addressed and tested. | Yes. |
| R7 | Addressed (Phase 2 gates, px/token table). Missing the synthetic-shape precision probe INTEGRATION recommends. | Yes. |
| R8 | Addressed (D9 descoped). | Yes. |
| R9, R10 | Pending M1/M8/M9. | Yes, marked pending. |
| K12 | Ran with reward v2 (same snapshot as R2). | **No**: "production reward" is overstated. |

## Re-registrations

**K8b: legitimate.** The feasibility analysis shows the absolute form fails even with a per-item oracle tolerance for any non-step shape (2.1 % of items have a rater outside the 2τ-dilated vote). Ranking is what GRPO consumes. Post-hoc elements the plan must say: the new thresholds (median ≥ 0.90, p10 ≥ 0.70) were set after measurement, and the erode direction (97.1 %, LIDC 95.2 %) is borderline and should be reported next to the 99.9 % dilate number.

**R3 split: legitimate but overclaimed.** The union/hull/core need every rater's mask; their band is 1.0 by construction, so no monotone band can beat them, and "consensus-like answers score well" is partly intended. The split is acceptable **only if** the eval hack audit carries a numeric alarm (area outside [A_min, A_max], |M|/|E|) in Phase 3; the plan says "monitored" without a threshold. The supporting local-gradient claim is overstated: `inflation_incentive.py` moves one rollout in 6 groups per item on ≈1,000 **fixture** items (not the 251 LIDC items) and reports the **mean** ΔA with a CI (all negative, −0.12 to −0.89). It does not compute a per-item/per-rollout rate, so "≥ 95 %" (§7 R3a) and "always" (§9) are unearned; and it is first-order evidence (uniform dilation of one honest copy), not proof of unreachability.

## Reward v4 flaws, none blocking

- **Area guard, zero tolerance.** `area_guard` uses all raters, so jittered copies are in-range by construction; no leave-one-rater-out test of g_area exists for v4 (the iteration-3 probe showed a held-out largest rater exceeds the others' max by 0.17 τ·P median). The guard bought 13–19 points on hedges that still fail, cost 3–9 points on hole adversaries, and moved 1-px neutrality from −0.003 to +0.047 (the τ·P variant stayed at −0.001). Defensible for the 1.5× steeper local gradient, but a LOO g_area test belongs in M6 with a stated switch rule.
- **Holes stamped 64 % vs the self-set 60 %.** A consensus-geometry answer with a rater-drawn topology; needs all raters' holes. Not a learnable hack; M6 spike is the right handling. State the fallback if the spike fails (keep v4, monitor).
- Exact-mask neutrality for v4 is −0.014 [−0.020, −0.008] (CI excludes 0; tiny). Mention it.
- Single-rater exclusion: fine (LIDC items are all k ≥ 2 after the multi-reader filter; RIGA k = 6).

## Architecture and milestones

Sequencing is correct (M1–M3 independent; M6 after M3/M4; M7 after M2/M5/M6). Missing: ExpertCache storage decision (on-the-fly 12 ms vs disk); a group-contiguity-per-process test in M7 (the "batch gather" risk is moot under v4, which has no batch z); the R7 synthetic precision probe as a Phase 2 entry item; numeric hack-audit thresholds in M8.

## Novelty

Correctly narrowed and matches novelty_check §5. Two wording edits (below).

## Required edits to PLAN.md

1. **§2 D7 evidence and §7 R2/K12:** replace "✅ 8.7e-8" / "✅ (from a non-collapsed SFT checkpoint)" with "✅ with reward v2 (`prod.py`, iteration-2 rewardcore); v4 equivalence and K12 re-run are the M7 gate".
2. **§7 R3(a):** replace "dilating or eroding an honest rollout lowers its advantage in ≥ 95 %" with "mean ΔA of a dilated/eroded honest rollout < 0 with 95 % CI excluding 0 for all six moves (0.25τ, 0.5τ, 1 px), ≈1,000 fixture items × 6 groups; per-rollout rate and the 251-item run added in M6". **§9:** "always lowers" → "lowers on average (CI excludes 0)"; "cannot reach by local improvement" → "no first-order incentive".
3. **§2 D5 evidence:** cite v4 realistic neutrality +0.038 [−0.011, +0.091], not +0.009.
4. **§7 K8b:** add erode 97.1 % (LIDC 95.2 %) and "thresholds set after measurement".
5. **§4.3:** add "LOO g_area test: distribution of g_area for each held-out rater; if median < 0.9, switch to the τ·P tolerance (neutral at −0.001)". Add the storage decision for the per-prompt cache.
6. **§4.7 / §9:** give the hack audit numeric alarm lines (share of rollouts with area outside [A_min, A_max]; stamped-hole-like outputs) used by the Phase 3 alarm.
7. **§7 Phase 2 entry:** add the synthetic-shape precision probe (boundary error in crop px per model/setting).
8. **§8:** item 1 prefixed "To our knowledge (targeted search)"; item 2 "Ours:" → "Our instantiation:" and no novelty implied for the LOO estimator.
9. **§10:** "313 train slices" → "251 under the size filter (313 unfiltered)"; "multi-process batch gather untested" → "group contiguity per process untested (no batch z in v4)".
10. **§7 R1:** "✅ 6,338 fixture targets; full train set at M5".
