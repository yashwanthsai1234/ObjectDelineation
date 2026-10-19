# Judge verdict on the Phase 1 decision record

Reviewer role: independent senior ML research engineer. Read in full: DECISIONS.md, KPIS_DRAFT.md, the six lab reports, ITERATION2.md, a_tables/b_tables, the six research notes (takeaways and recommendations), ROADMAP.md, the novelty report, RESULTS_0b.md. Code spot-checked: `common/oracle.py` and its tests, `t_polycore/tests/*` (K1, K2, helpers, compare_libs), `t_reward/rewardcore/{band,advantages,reward,topology,perturb,reference}.py`, `t_reward/scripts/{k9_discrimination,k9_v2,k8_band,k8_adjusted,topo_neutrality}.py`, `t_reward/tests/{test_reward_v2,test_band}.py`, raw `results/iter2/*.txt`, `t_dynamics/{grpo_toy,test_grpo_toy,run_b_alarm,curriculum_sim}.py`, `t_stack/{rewards,t3_grpo}.py` and `t5_checks.json`, `t_targets/targets/builder.py` imports and test list, `t_data/splits/test_splits.py` test list, `lidc_parity/SUMMARY.json`.

Overall: the engineering core (polycore, topology, targets, data) is well built and well tested. The reward design (F5, DM-loo topology, IoU_tc) rests on thin and partly circular evidence, has never been run through TRL, and the record overstates two results. The plan can be written, but four things must become Phase 1 exit gates before any GPU is rented.

---

## Scores by axis

### A. Correctness and rigor of evidence: 3 / 5

What is solid:
- K1/K2 are genuinely adversarial. The oracle is a pure-integer spec implementation; polycore is cross-checked against it on 10k shapes with on-edge centres, half-pixel lattices, degenerate rings, and additionally against third parties where a third party exists (skimage/rasterio off-boundary, GEOS validity and containment, the numba prototype), plus 1M-shape ASan/UBSan fuzz. The on-edge tie rule is only checked against the in-house oracle, which is unavoidable (it is a convention), and the report correctly notes the asymmetry means GT and predictions must share the rasterizer.
- Topology: exact method cross-checked against a scipy brute force and a pure-Python BFS, and the bit-quad χ certificate is independent. 0/10,000 + 0/8,738.
- K8 leave-one-rater-out is not leaky: `loo_item` builds `ExpertCache` from `rest` only, so T5 never sees the held-out rater; k=2 cases use a cross-fitted size prior. Cluster bootstrap by item is the right CI.
- t_dynamics: 20 seeds, CIs, gradient checked against finite differences, a replica of the repo's `group_advantages` verified equal. The collapse mechanism ("rich get richer" under set reward) is demonstrated, not asserted.

What is not solid:
1. **K8 goal-post move.** Pre-registered K8 was pooled |ρ| ≤ 0.10. F5 scored +0.107 [0.088, 0.126]: a clean fail. The record then adopts "partial ρ given rater count", under which F5 passes but so do F2, F3 and F4, i.e. the adjusted metric no longer discriminates between candidates, and F5's selection is justified only by the pooled metric it failed. Within-protocol scale independence is a defensible definition (resolution mixing between CT and fundus is a real confound), but it was chosen after seeing the result. Verdict: legitimate definition, illegitimate process. Must be re-registered as a per-stratum KPI, and a separate discrimination KPI must carry the weight of choosing F5.
2. **F5 self-calibration is partly circular with K8.** T5 is chosen so that a new expert's expected score equals s* = 0.95 on every image, independent of size. The LOO median ≥ 0.9 and ρ ≈ 0 are therefore nearly by construction; K8 tests the Gaussian order-statistic correction S_m/S_(m−1), not whether F5 is a good reward. What matters for a reward (monotone, non-saturating, discriminating beyond rater spread) is only covered by the synthetic monotonicity table in k8_band, which is not a KPI.
3. **K9 construction.** (a) Pairs are not independent (dozens per item), and no CI is reported. (b) Every group in k9_v2 holds 1 wrong + 7 correct fillers, so the DM-loo term always assigns the wrong sample the minimum; K9-v2 largely tests whether geometric terms can outweigh c_topo, not discrimination. (c) The `topo_conditional_iou` fallback (no rater with the sample's topology → max over all raters) means a wrong-topology sample built from a rater's mask keeps IoU_tc ≈ 1 (ITERATION2 §K9, p5 gap −0.19). This is a design flaw that c_topo = 2 papers over. (d) Defects that happen to match a rater's topology are discarded (k9_discrimination.py L104–105), so the one exploit specific to set-valued rewards (punch a hole to match the minority rater) is never scored.
4. **Neutrality claim is overstated.** DECISIONS says "minority answer no longer disadvantaged (0/14 items, was 10/14)". The raw file shows: on exact rater masks IoU_tc and band are constant, so only the DM term varies, and v2 gives a statistically significant pro-minority bias (+0.082 [+0.063, +0.103]) that grows with c_topo. On jittered masks (the realistic case) v2 c_topo=2 is −0.144 [−0.504, +0.167], minority loses in 7/14 items, CI spans ±0.5. That is "underpowered and not neutral", not "fixed". n = 14 items with single-rater minorities cannot support the choice of c_topo.
5. **K12 passed with the wrong reward.** `t_stack/rewards.py` imports `objdel.checker.polygon_io` (PIL rasterizer, the one D1 rejects) and `objdel.reward.AnnotatorSet` (vote IoU, topo_near, the terms D5/D6 reject). The gdpo-floor aggregator in D7 was never implemented in TRL; only `normalize_then_sum` and the controller callback were exercised. "2–4 GRPO steps run on 4 models with per-term logs" is true of the old reward only.
6. **Controller evidence is simulation with invented dynamics** (REPORT caveat) on a 600-step horizon with eval every 10 steps. The alarm's median delay is 270 steps = 27 evals, longer than the 300-step GRPO runs the roadmap budgets. Thresholds were tuned on the same simulator.
7. **Tests calibrated to pass.** `test_no_topology_class_favoured` asserts |diff| < 0.15 with measured +0.08; `test_scale_free_for_uniform_excursion` asserts spread < 0.12. These are regression pins, not acceptance tests.

### B. Optimality / cutting-edge choice per module: 4 / 5

Every decision matches or improves on the research notes' recommendation; I found no clearly better known tool missed:
- D1/D2: exact integer scanline with winding prefix sum, graded pair classification, no repair. Correct; nothing in the notes beats it. numba was 5× slower and still 25× under budget, so the C++ core is over-engineered for Phases 1–3 on its own; it is justified only by the CUDA shared header and dependency hygiene. Accept, but it must ship as a wheel.
- D3: Spaghetti CCL + background labeling with pixel-above assignment; faster than the bit-quad χ and exact. Correct.
- D4: distance-weighted, perimeter-normalized band with rater-derived tolerance follows NSD/Boundary-IoU evidence. The rater-count correction is new and sensible. The notes also recommended a soft boundary-in-band term and a band gate (κ) instead of an IoU gate; neither is in the record.
- D5–D7: match the Q2/Q4 evidence (GDPO + floors + gates; mean-centre the dense term; distribution matching). Missing from the research note's own list: the multiplicative area/part-inflation guard (LocPen/odLength analogue), std-threshold group dropping (described, not implemented in TRL), and a size filter on the topology signature (Phase 0 used 5 px / 100 px; rewardcore uses none).
- D8: TRL 1.16 is the right single stack for SmolVLM; EasyR1/verl would be better throughput for Qwen only. Dr. GRPO + clip-higher + small β is the consensus. Fine.
- D9: the lab correctly demoted typed sampling (its own negative result). Learnability sampling needs a TRL sampler subclass, not a callback; untested.
- D10: compact JSON, 0–1000 relative, canonical order, no caps: matches the output_format note. Fine.
- D11: potrace polygon stage + lossless DSS + exact recertification beats DP as the note predicted. Good. Hybrid with repair is heavier than needed for a first SFT but is built and tested.
- D12/D13: consistent with the notes.
- D14: correct split, but now three implementations of one spec exist (polycore C++, fastgeom numba in t_targets, oracle), plus the repo's PIL path still used by t_stack.

### C. Novelty: 3 / 5

Ranked strongest to weakest:
1. **Set-valued per-instance topology of the policy's own output, against expert rater sets, with Phase 0 evidence that raters disagree.** Still the headline and still unanticipated by anything readable. Weak spot: VectorLLM++ remains unread (the record lists it as a risk; it must be read before the plan freezes claims).
2. **Rater-count-corrected self-calibrated band (F5).** NSD already uses inter-observer tolerance as a metric; the per-image, order-statistic-corrected version as an RL reward is a modest but real contribution. Needs validation on a second multi-rater set (the record admits the Gaussian assumption).
3. **Distribution-matching topology term.** Honest description: `r = 1[valid]·(1 + λ(log f(a) − log π̂_loo(a)))` is FlowRL's reverse-KL-to-target reward with the target set to rater frequencies, estimated by a leave-one-out group count. The "share" variant they tested is literally fitness sharing (Goldberg & Richardson 1987). This is an application of a known family, not a new algorithm; the risk of a reviewer calling it "FlowRL renamed" is real and should be pre-empted by naming it that way ("FlowRL-style rater-frequency matching with a group LOO estimator") and citing FlowRL, Multi-Answer RL, DPH-JS and fitness sharing. The defensible novelty is the target (verified rater histogram) and the empirical finding that coverage bonuses over-spread while frequency matching does not.
4. **IoU_tc.** This is D_max (max-over-annotators Dice, already in the mask_algorithms note's eval list) restricted to same-topology raters. A design detail, not a claim.
5. **Typed alarms + rollback.** Weakest. Components all have precedents (ODF, PODS, NRFR alarms, best-checkpoint early stopping); the lab's own result is that typing does not help sampling; evidence is simulated. Present as engineering, not as a contribution.
6. **Bit-exact shared-spec core.** Solid engineering claim; not ML novelty.

### D. Goal alignment: 3 / 5

Serves the goals: smallest-first screening, CPU reward at ~10–50 ms per 256 rollouts, 4090-class budget, rater-set claim carried by LIDC.

Missing for Phase 2 to start:
- **No quantitative SFT gate.** ROADMAP says "≥ 90 % valid, IoU near MedSAM"; nothing defines parse rate, gate pass rate (IoU_tc ≥ 0.5), duplicate-output rate, or group reward spread. The 4/4 overfit result is not evidence of learnability.
- **No eval harness.** The terms exist in rewardcore, but there is no module or test for the eval protocol the research note requires (D_max, per-image TV/JS to the rater histogram over K samples, pass@k, hack audit: |M|/|E|, part count, duplicate-hash rate). Without it Phase 2/3 cannot be judged.
- **No target→reward round-trip.** Nothing checks that every SFT target text, parsed by the production parser and scored by the production reward, scores at the top. That is the single-spec guarantee in one test, and it does not exist.
- **Topology-disagreement volume.** 14 fixture items; Phase 0 found 326 disagreement slices in 4,087, and the split reports 69.7/10.1/20.2 % of disagreement slices, so ~200+ train items exist. The reward was tuned on 14. The full LIDC build (all 1,018 series; parity checked on 67) has not been run.
- **Small-model precision feasibility.** SmolVLM-256M at `longest_edge=512` sees 64 image tokens, i.e. one token per 64×64 px. Producing 0–1000-grid boundaries on a 500 px crop from that is doubtful, and on LIDC the IoU_tc ≥ 0.5 gate would then silence the topology term for most rollouts (MedSAM reached IoU ≥ 0.5 on only 22 % of nodule slices). This needs a pre-check, not a hope.
- Vast x86 timing is still unmeasured; `ExpertCache` memory (two float32 EDTs per RIGA item ≈ 2 MB × 4,500) needs a storage decision (on-the-fly at 12 ms median is fine; say so).

### E. Engineering robustness: 4 / 5

Strong: test-first in every lab, determinism (K10 across threads/processes/interpreters), Linux Docker on arm64 and emulated amd64, ASan/UBSan fuzz, exceptions never crash, split tests with negative controls, SHA-256 manifests.

Weak:
- Four implementations of the raster/legality spec in play (polycore, fastgeom, oracle, repo PIL). Only one should certify anything in production.
- K13 (coverage ≥ 90 %, CI green) is not measured anywhere.
- No cibuildwheel; Vast will build from source.
- numba becomes a hard dependency through the target builder (offline, acceptable, but pin it).
- OpenCV 5 API wrinkle (`ccltype=`), opencv-python-headless not yet a core dependency.
- TRL 1.16 bug #6294 blocks vLLM for SmolVLM; HF generation at 17.6 s/step on the Mac; the 4090 number is unknown.
- Two floating-point tie-breaks in the target builder: pin and hash (the lab says so; make it a test).
- `hash(tuple)` for signature ids is deterministic for int tuples but should be replaced by a dict index (as topo_neutrality does) to remove any doubt.

### F. KPIs: 3 / 5

Right families, wrong rigor in three places: K8 redefined post hoc; K9 lacks CIs and the set-specific adversaries; K12 passed with the old reward; K13 unmeasured. No Phase 1 exit gates and no Phase 2 entry gates exist. Proposals below.

---

## Verdict per decision

| # | Decision | Verdict | Required change |
|---|---|---|---|
| D1 | polycore raster | **ACCEPT** | Ship as a wheel (cibuildwheel, CI per OS); make polycore the only production rasterizer: targets certify with it (R1), repo PIL path deleted. Cap inputs (≤ 1,000 vertices/ring, ≤ 10 rings) and test the cap. |
| D2 | legality | **ACCEPT-WITH-CHANGES** | Reward v2 uses a hard legality gate, so "graded" is logging-only: say so, and keep per-class counts in logs. Decide and test the two open semantics (wrong_winding on net area; hole crossing its outer gets both flags). Decide whether 8-connected pinches are legal in outputs (targets bevel them; the reward must agree). |
| D3 | mask topology | **ACCEPT** | Add opencv-python-headless to core deps; pin the OpenCV 5 API. Decide the signature size filter (R6). |
| D4 | band F5 | **ACCEPT-WITH-CHANGES** | Re-register K8 per stratum (R5); add a discrimination KPI; a rater-subsample test of the S_m correction on RIGA; a policy for single-rater items (exclude from RL or F4 with prior); validate on IMA++ or QUBIQ when available. |
| D5 | topology term | **ACCEPT-WITH-CHANGES** | Rename and cite as FlowRL-style rater-frequency matching with LOO estimator. Re-tune c_topo and λ on ≥ 200 disagreement items (R4), and treat c_topo as a swept hyperparameter in the first real run with per-image TV as the metric. Report jittered neutrality honestly. |
| D6 | IoU_tc | **ACCEPT-WITH-CHANGES** | Change the fallback for a topology no rater drew from "max over all raters" to IoU vs the vote (or min over raters); re-run K9 at c_topo ∈ {1, 2} (R3). |
| D7 | advantages | **ACCEPT-WITH-CHANGES** | Implement the aggregator inside TRL and prove it bit-identical to `rewardcore.advantages_v2` (R2). Implement "drop low-spread groups" as zero advantages and test it. Add an area/part-inflation guard or show the adversary suite does not need it (R3). |
| D8 | RL loss | **ACCEPT** | Keep GSPO as an ablation; note that `beta=1e-3` is a workaround for the full-FT crash, not a tuned value. |
| D9 | controller | **ACCEPT-WITH-CHANGES** | Descope for Phases 2–3 to typed logging + best-validation-checkpoint rollback + manual guard switch. Re-express alarm thresholds in eval units and re-tune on real curves; the simulated 270-step lag is unusable at 300 steps. Learnability sampling needs a TRL sampler subclass and a test, or is dropped. |
| D10 | output format | **ACCEPT** | Fix the grammar (`"holes":[]` always); truncation at `max_new_tokens` counts as a format failure; add the canonicalization-invariance test (start vertex / direction changes must not change the reward). |
| D11 | SFT targets | **ACCEPT-WITH-CHANGES** | Certify with polycore, not fastgeom; keep fastgeom only as a differential test (R1). Pin numba and library versions; hash the target set. Lossless-vs-hybrid SFT ablation stays in Phase 2. |
| D12 | model shortlist | **ACCEPT-WITH-CHANGES** | Define the SFT gate numerically (R7) and add an image-token precision feasibility check before spending on SmolVLM-256M. Confirm the LFM1.0 licence before any run. Keep Qwen3-VL-2B as the non-smallest control for the spurious-reward comparison. |
| D13 | data | **ACCEPT-WITH-CHANGES** | Run the full LIDC build with the four parser fixes; parity on ≥ 200 scans; report disagreement item counts per split (R4). Document the outline rule. Storage: compute ExpertCache on the fly (12 ms) or cache on disk; decide. |
| D14 | C++ vs Python | **ACCEPT** | One production implementation per operation; prototypes retired to `prototypes/`. |

---

## Top 5 risks, ranked

1. **The headline reward is tuned on 14 items.** c_topo = 2, the IoU_tc fallback, and the neutrality claim all rest on 14 LIDC fixture items with single-rater minorities. The jittered neutrality CI spans ±0.5. If the full-data re-run flips the sign, the Phase 3 arms were designed on noise.
2. **Reward v2 has never run in TRL.** The smoke test used the old PIL rasterizer and vote-IoU reward; the aggregator is unimplemented. Integration bugs (group boundaries inside reward functions, `scale_rewards="none"` semantics, weight-0 loggers) will surface on rented GPU time, and single-spec drift already exists in the repo.
3. **The smallest model may produce no gradient signal.** At 64 image tokens the IoU_tc ≥ 0.5 gate may silence topology on most LIDC rollouts, and near-constant groups will be dropped. No SFT gate or feasibility check exists, so "smallest viable" is undefined.
4. **F5/K8 evidence is partly circular and post hoc.** The band is calibrated to give 0.95 by construction and the KPI was redefined after a fail. If the Gaussian order-statistic correction is wrong on 2-rater LIDC items, the band is systematically lenient or harsh on the exact subset that carries the topology claim.
5. **Novelty framing.** The DM term is FlowRL/fitness sharing with a rater-frequency target; VectorLLM++ is still unread; the controller's contribution is simulated and lagging. A reviewer will raise all three. The claims must be narrowed before the plan hard-codes them.

---

## Required changes before the plan is written

Each item says whether it needs a new test and describes the test precisely.

**R1. Single-spec round-trip (new test, Phase 1 exit gate).** For every SFT target in RIGA and LIDC train (all raters, both lossless and hybrid modes): serialize → parse with the production parser → rasterize with polycore at the crop size → score with rewardcore against that item's rater set. Assert, over 100 % of targets: legal; per-instance signature equals the source rater's signature; IoU_tc ≥ 0.99 (= 1.0 in lossless mode); F5 ≥ 0.95 for ≥ 99 % of targets. Second assertion: polycore raster of each target equals the fastgeom raster bit for bit; after it passes, fastgeom is demoted to `prototypes/`. Third assertion: applying a canonicalization-preserving change (rotate ring start, reverse-and-fix winding) leaves every reward term unchanged.

**R2. TRL aggregator equivalence (new test, Phase 1 exit gate).** Build a synthetic batch of B = 4 prompts × G = 8 with hand-set term vectors: one group with one illegal rollout, one group all-identical (must be dropped → zero advantages), one group with a lone minority topology, one normal. Run it through the exact TRL path (weight-0 per-term reward functions + aggregator reward function, `scale_rewards="none"`, `sum_then_normalize`, `num_generations=8`) with a dummy model, capture the advantages the trainer computes (hook or logged `advantages`), and assert equality with `rewardcore.batch.advantages_v2` to 1e-6. Also assert the per-term logs match the term vectors.

**R3. IoU_tc fallback fix + adversary suite (code change + new test).** Change the fallback for an unmatched topology to IoU vs the vote. Re-run K9 at c_topo ∈ {1, 2} and report K9 with a cluster bootstrap CI over items. Extend `wrong_variants` with set-specific adversaries scored as groups: union dilated by 2τ; convex hull of the union; vote with every rater's holes stamped in; vote plus a speck inside the band but detached; vote with the minority rater's hole punched (topology matches a rater: must not beat that rater's own jittered mask); ring with reversed winding. KPI: each adversary's advantage is below the median honest jittered rollout in ≥ 95 % of groups, with CI. If the dilated-union or hull adversary wins, add the area-inflation guard.

**R4. Disagreement-scale re-validation (data build + re-run, Phase 1 exit gate).** Run the full LIDC parser with the four fixes over all 1,018 series; parity vs pylidc on ≥ 200 scans; report per-split counts of items with per-instance rater disagreement. Re-run topo_gap, topo_neutrality and K9-v2 on every train disagreement item (expected ≥ 200). Acceptance: jittered neutrality 95 % CI half-width ≤ 0.10 and the interval contains 0 for the chosen c_topo; otherwise re-choose c_topo/λ and record the choice. Freeze c_topo only after this.

**R5. K8 re-registration + band discrimination (KPI change + new test).** K8 := within each stratum (riga/disc, riga/cup, lidc k=2, k=3, k=4) |Spearman ρ(F5, area)| ≤ 0.10, pooled ρ reported but not gated. New KPI K8b (discrimination): on ≥ 95 % of items, F5 of the vote mask moved by a uniform 2τ excursion ≤ 0.50 and F5 of every held-out rater ≥ 0.85; monotone non-increasing in |excursion| on ≥ 98 % of items. New test for the rater-count correction: on RIGA (6 raters) subsample m ∈ {2, 3, 4, 5} raters per image, compute the LOO F5 of the held-out expert; the median must stay in [0.92, 0.98] for every m. Decide the single-rater policy (exclude from RL, or F4 with a documented prior) and test that the cache builds deterministically for k = 1.

**R6. Signature size filter (decision + new test).** Pick one `min_area` rule (absolute px by dataset, or a fraction of |U|) and use it in the Phase 0 census re-run, the targets, the reward and the eval. Test: a 1-px speck added to a rater mask below the filter does not change the item's topology set; a component above the filter does; the reward's `topo_any` for a prediction is invariant to sub-filter specks.

**R7. SFT gate and precision feasibility (definition + Phase 2 entry gate).** Before the first SFT run, for each candidate model record px per image token at the chosen crop size. Phase 2 entry gate on the validation split after SFT: parse ≥ 98 %; legal ≥ 90 %; median IoU_tc ≥ 0.85 (RIGA disc), ≥ 0.75 (cup), ≥ 0.60 (LIDC); share of rollouts passing IoU_tc ≥ 0.5 ≥ 70 %; topology-any ≥ 90 % on unanimous items; duplicate-output rate ≤ 10 %; at temperature 1 with G = 8, ≥ 50 % of groups have geometric σ above the floor (otherwise GRPO has no signal and the model is not viable). The smallest model that passes is "smallest viable"; the next size up is the control.

**R8. Controller descoping (plan change + one test).** Phases 2–3 ship typed per-term logging, held-out typed error curves, and best-validation-checkpoint rollback only. The alarm is re-expressed in eval counts and re-tuned on real curves in Phase 5. If learnability sampling is kept, add a test that the TRL sampler subclass changes per-item visit counts as intended across 2 steps.

**R9. Portability and CI (Phase 1 exit gate).** Linux wheel via cibuildwheel (x86-64 and arm64); CI runs pytest + sanitizer fuzz per OS; K13 coverage actually measured; K3/K4 re-measured on a Vast x86 host with 8 cores; pinned versions for TRL, transformers, torch, numba, OpenCV, with the target-set SHA-256 recomputed under the pins.

**R10. Eval harness (new module + tests, Phase 1 exit gate).** One module computing: IoU/Dice vs majority, D_max, F5, topology-any, per-image TV and JS between the K-sample topology histogram and the rater histogram, minimum-mode retention, pass@k for k ≤ 64, and the hack audit (|M|/|E|, component count, duplicate-hash rate, NRFR-style newly-rewarded failures). Tests: metrics are zero or one on the raters' own masks; TV = 0 when samples reproduce the rater histogram exactly; the audit flags the R3 adversaries.

---

## Proposed KPI additions

Phase 1 exit gates (all must pass before renting a GPU): R1 round-trip 100 %; R2 TRL equivalence 1e-6; R3 K9 with CI ≥ 95 % and adversary suite ≥ 95 %; R4 neutrality CI half-width ≤ 0.10 on ≥ 200 items; R5 K8 per-stratum and K8b; R6 size-filter tests; R9 wheel + CI + Vast timing (K3 ≤ 250 ms, K4 ≤ 5 ms on x86); R10 eval harness tests green; K13 ≥ 90 % measured.

Phase 2 entry gates: R7 thresholds; VectorLLM++ read and the novelty sentences updated; LFM licence resolved; per-step wall time on the 4090 measured for the smallest model with G = 8 and 1,536 max tokens, and the per-run GPU-hour estimate in the roadmap replaced by the measurement.

---

## Go / no-go

**GO, conditional.** Write the plan. The plan must list R1–R6, R9 and R10 as Phase 1 exit gates and R7 as the Phase 2 entry gate, must restate the novelty claims in the narrowed form above, and must mark c_topo, λ and the IoU_tc fallback as provisional until R3/R4 are done. No GPU is to be rented before R1, R2 and R4 pass.
