# M6: reward v4 (`objdel.reward`)

Plan: D4–D7, §3.2, §4.3, gates K3, K4, K8, K8b, K9, K10, R3, R4, R5. Promoted from `lab/t_reward` (reward v4 as
frozen in `results/ITERATION4.md`) and `lab/t_stack/prod.py` (`score_text`, low-spread drop).

## Result

| Gate | Requirement | Measured (gate tables below; Apple M3, load average 9) | Verdict |
|---|---|---|---|
| K3 | 256 completions ≤ 250 ms | serial: mixed 48.6 ms (max 51.2), RIGA-only 103.5 ms (max 127.7); text → parse → polycore → terms → advantages | PASS (Mac; x86 at M9) |
| K4 | p99 per rollout ≤ 5 ms | p99 0.53 (mixed) / 0.60 ms (RIGA); sandbox at load 14: RIGA p99 3.9 ms, max 9.8 ms | PASS |
| K8 | per-stratum abs ρ(F5, area) ≤ 0.10 | RIGA disc +0.008 [−0.028, +0.046], cup −0.040 [−0.077, −0.002], LIDC k=3 +0.012 [−0.022, +0.045], k=4 −0.023 [−0.048, +0.006] (equal to the lab) | PASS |
| K8b | monotone ≥ 98 %; 2τ vote < median held-out ≥ 95 %; median ≥ 0.90, p10 ≥ 0.70 | 100 %; 100 % (eroded 99.2 %); 0.975 / 0.824 | PASS |
| K9 | ≥ 95 % with CI | fixtures + composites 99.9 % [99.8, 100.0] (8,018 pairs); full-LIDC 98.9 % [98.0, 99.6] (789 pairs, 229 items; lab 98.86 [98.00, 99.62]) | PASS |
| K10 | identical SHA-256 | serial ×2, 8 threads, 3 processes, OpenCV 8 threads, fresh `-I` interpreter (gate script and unit test) | PASS |
| R3 (a) | dilated union, in-band speck below median honest ≥ 95 % | 100.0 % [100.0, 100.0]; 96.4 % [94.9, 97.5] (lab 96.5 [95.2, 97.7]) | PASS (speck CI low 94.9 %) |
| R3 (a) | local gradient on average | every move's mean ΔA < 0 with CI < 0 (251 items and fixtures) | PASS |
| R3 (a), M6 | per-rollout local-gradient rate ≥ 95 % | **81.6 %** (251 items, 26,350 moves), 84.1 % (fixtures); 97.6 % / 98.3 % for moves away from the median rater area, 64 % / 69 % for moves toward it | **FAIL** |
| R3 (b) | reported | adversary beats the median honest rollout: hull 34.3 %, union 50.2 %, core 44.2 %; holes stamped win 62.7 % [58.2, 67.4] (re-registered limit 60 %: FAIL as in the lab), minority hole 54.5 % [48.7, 60.1] (PASS) | reported |
| R4 | ≥ 200 items, CI contains 0, half-width ≤ 0.10 | 251 items. 1-px: +0.007 [−0.047, +0.061] ACCEPT. Realistic: +0.056 [+0.014, +0.100] reject in the gate run; 5 more pool draws: +0.037, +0.066, +0.031, +0.030, +0.037 (4 of 5 accept) | **borderline** (1-px PASS; realistic passes in 4 of 6 draws) |
| R5 | medians in [0.92, 0.98] for m = 2–5; single-rater policy | 0.961, 0.971, 0.976, 0.980 (m = 5 at the edge, as in the lab); k = 1 refused by `build` | PASS |

**Superseded by `fixes/F2.md`** (break-it B2): sub-filter raters abstain from the topology vote (B2-4), so the
251-item gate set is 239; validity gated on `iou_min` (B2-7). Re-run: K9 full-LIDC 99.3 % [98.5, 99.9] (228 items);
R3(a) speck 96.3 %; holes stamped 61.6 %; R4 realistic +0.033 [−0.013, +0.078] ACCEPT, 1 px +0.050 [+0.004, +0.095]
reject (still borderline); local gradient away 97.3 %. The gate tables below are M6's run.

**A/B decisions.**
- **AB2 (Phase 0 / v2 / v4): v4 confirmed.** Phase 0 loses the in-band speck (52.8 %), is strongly majority-biased
  (R4 −0.563 / −0.458 / −0.178) and ranks topology worse on full LIDC (K9 90.7 %). v2 fails realistic neutrality
  (+0.111 [+0.062, +0.161]) and K9 is lower (96.2 %); it lets the hull, union and core win more often (below
  median 51 / 37 / 40 % vs v4 66 / 50 / 56 %). v2 is slightly better on the hole adversaries (holes stamped
  60.4 % vs 62.7 %), the known cost of the guard (lab ITERATION4).
- **AB3 (hole-matching term): rejected.** Topology-conditional consensus holes make the hole adversaries win
  more (66.4 % / 58.4 %) and break neutrality (1-px −0.130, realistic −0.054); all-rater consensus holes leave
  holes-stamped unchanged (62.7 %) and break neutrality badly (−0.23 to −0.29). Kept only in
  `tests/oracles/reward_variants.py`.
- **AB4 (g_area zero vs τ·P): keep zero tolerance.** Leave-one-rater-out medians 1.000 (RIGA), 0.991 (LIDC
  fixtures), 0.971 (full-LIDC disagreement) ≥ 0.9, so the switch rule does not fire; p10 is 0.895 / 0.757 / 0.637
  (held-out raters outside the others' area range pay up to ~35 %). τ·P passes 97–100 % of held-out raters but
  also every hedge (hull below median 51 % vs 66 %).

## What was built

`src/objdel/reward/` (≈ 470 lines including docstrings) implements reward v4 only:

| File | Contents |
|---|---|
| `spec.py` | `REWARD_V4` (frozen dataclass of every constant of plan §4.3), `TERMS`, `fingerprint()` (SHA-256 of REWARD_V4, the size rule probed at fixed areas, numpy and OpenCV versions; + scipy since F2, B2-8) |
| `cache.py` | `ExpertCache` (`build`, `save`, `load`, pickles as masks + τ + T5), `rl_eligible`, `boundary`, exact `edt`, `band_depth`, `excursion_moment` (S_m), `rater_tau`, `f5_tolerance` (T5) |
| `score.py` | `ObjectTerms`, `score_mask` (IoU_tc, F5, g_area, signature, f), `area_guard`, `sig_id`, `score_completion`, `score_batch` |
| `advantages.py` | `advantages_v4`, `topology_reward` (r_topo, for logging), `loo_bias` (b(f, G), exact binomial sum) |

Deleted: `src/objdel/reward/annotators.py` (Phase 0 `AnnotatorSet`, `group_advantages`) and its test. Its logic
survives only as the `phase0` arm of `tests/oracles/reward_variants.py`, checked equal to the original in the
sandbox (`scratchpad/impl/M6/ab2_phase0_check.py`: 41 items × 8 rollouts, per-rollout terms equal, max advantage
difference 0).

Tests and oracles: `tests/test_reward.py` (44 tests, no data, ~10–20 s), `tests/test_reward_gates.py` (data:
9 fast + 1 strict xfail + 3 slow), `tests/oracles/reward_reference.py` (independent slow reference),
`tests/oracles/reward_harness.py` (honest rollouts, adversaries, K9 pairs, gate computations, ported from the lab
scripts), `tests/oracles/reward_variants.py` (Phase 0, v2, τ·P guard, hole-matching arms for AB2–AB4).
Gate script: `experiments/phase1/reward_gates.py` (writes the tables at the end of this note).

## Contract (EXECUTION.md §3, unchanged; details M7/M8 depend on)

```python
from objdel.reward import REWARD_V4, ExpertCache, score_completion, score_batch, advantages_v4, topology_reward

cache = ExpertCache.build(stack)            # stack bool (k >= 2, H, W); ValueError for k < 2 (rl_eligible)
cache.save(path) -> sha256 ; ExpertCache.load(path)   # ValueError on SHA mismatch; rebuilds if the fingerprint changed
score_completion(text: str, truncated: bool, names: Sequence[str], caches: dict[str, ExpertCache], W: int, H: int) -> dict
  # format_ok, legal, iou_tc, f5, g_area, topo_any, topo_freq, sig_id (floats); error (str); defects (dict[str, int])
score_batch(jobs: Sequence[tuple], pool: Executor | None = None, chunk: int = 8) -> dict
  # jobs = tuples of score_completion's arguments; numeric keys -> float64 arrays (N,), "error"/"defects" -> lists
advantages_v4(terms: dict[str, array (B*G,)], group: int = 8) -> np.ndarray (B*G,)
  # reads legal, iou_tc, f5, g_area, topo_any, topo_freq, sig_id; prompt-contiguous groups; ValueError otherwise
topology_reward(terms, group=8) -> np.ndarray (B*G,)   # r_topo after the legality and validity gates (logging)
```
*Superseded by F2 (EXECUTION.md §3 is current):* `ExpertCache.load(path, expected_digest=None, rebuild=False)`
raises `StaleCacheError` on a fingerprint change unless `rebuild=True` (B2-8); `score_completion` adds `iou_min`
(B2-7); `rater_answers` / `voters` (B2-4, B2-9); `score_objects` added for `objdel.eval`.

Semantics decided here:
- **Multi-object completions** (RIGA disc + cup): mean of the per-object IoU_tc, F5 and g_area; legal = every
  object legal; the topology answer is joint: `topo_freq` = share of raters whose signatures equal the
  prediction's for every object (rater r of one object's cache must be rater r of the others), `sig_id` =
  CRC32 of the tuple of per-object signatures. This is the lab rule (`t_stack/prod.py score_text`) and reduces to
  the per-object terms for LIDC.
- **Failures**: unparseable or truncated → `format_ok = legal = 0`, every term 0, `sig_id = -1`, `error` = the M2
  parser code. Parsed but illegal → `legal = 0`, terms computed and logged, `defects` = summed polycore counts;
  `advantages_v4` zeroes every term of a rollout with `legal < 0.5` (hard gate, as in the lab).
- **Empty object list** (M2: legal empty mask): scored as a mask; g_area = 0 for an empty mask (lab rule), IoU_tc vs
  the vote unless a rater drew nothing (no fixture rater mask is empty).
- **Low-spread drop**: the lab `prod.py` rule extended with g_area: a group whose four ungated terms all span
  < 1e-2 and whose rollouts share one `sig_id` gets advantage 0.
- K10 also covers the r_topo bias b(f, G): computed by an exact binomial sum instead of `scipy.stats.binom`
  (differs from the lab by ≤ 1e-13 in advantages).

## Research (sources read for this milestone)

| Topic | Finding used | Source |
|---|---|---|
| GDPO | per-reward group z, sum, then batch z; official TRL code divides by (std + 1e-4), zero-std reward → 0. v4 keeps per-term group z but replaces the epsilon by σ floors and drops the batch z (lab dark-room and batch-composition evidence) | https://arxiv.org/abs/2601.05242, https://github.com/NVlabs/GDPO |
| TRL | `normalize_then_sum` = GDPO incl. batch z and ignores `scale_rewards`; `sum_then_normalize` with `scale_rewards="none"` = mean-centring only, which is what the M7 aggregator needs | https://huggingface.co/docs/trl/main/en/grpo_trainer, https://github.com/huggingface/trl/blob/main/trl/trainer/grpo_trainer.py |
| Dr. GRPO | advantage R − mean(R), no std division (difficulty bias), constant length normaliser | https://arxiv.org/abs/2503.20783 |
| DAPO | dynamic sampling drops all-equal groups; we drop by a spread threshold because rewards are continuous | https://arxiv.org/abs/2503.14476 |
| MaD-RL | group-frequency target matching; reverse-KL reward log P*(z) − log P̂π(z) (Eq. 10), gated form R = R_corr (1 + λ R_div) (Eq. 17); notes the finite-group bias of nonlinear rewards. r_topo is this reverse-KL with P* = rater histogram, a leave-one-out P̂π and an exact bias correction | https://arxiv.org/abs/2609.31644 |
| LOO estimators | RLOO baseline; Mohri et al. 2026 show Jensen bias of plug-in log rewards and give U-statistic / polynomial estimators. Our b(f, G) is the exact bias at π = f for the LOO estimator (enumeration test: equal expected r_topo per class to 1e-9); their minimax estimator is a possible future replacement (proposal only) | https://arxiv.org/abs/2402.14740, https://arxiv.org/abs/2605.02435 |
| FlowRL | reward-distribution matching via a learned partition function; cited for the matching idea | https://arxiv.org/abs/2509.15207 |
| Inter-rater tolerance | Nikolov et al. set the surface-Dice tolerance per organ as the P95 of inter-rater distances (our τ); STAPLE (EM over rater performance) and the egg–yolk model (core/union band) | https://arxiv.org/abs/1809.04430, https://github.com/google-deepmind/surface-distance, https://doi.org/10.1109/TMI.2004.828354 |
| Exact EDT (superseded by S2/M9 F1: now `scipy.ndimage.distance_transform_edt` cast to float32, bit-identical on the Mac, deterministic on x86; OpenCV's IPP path is not, opencv#24082) | OpenCV DIST_L2 + DIST_MASK_PRECISE runs Felzenszwalb–Huttenlocher (exact; 3×3 / 5×5 masks are chamfer approximations). Pitfall: with no zero pixel it returns 65536 (we return 1e6, explicitly). scipy `distance_transform_edt` is Maurer et al. 2003 and returns meaningless finite values without background | https://theoryofcomputing.org/articles/v008a019/, https://github.com/opencv/opencv/blob/4.x/modules/imgproc/include/opencv2/imgproc.hpp, https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.distance_transform_edt.html |
| Segmentation RL hacks | Seg-R1 (S-measure hacked by all-black masks), SAM-R1 (over-segmentation), VLM-R1 (over-prediction; count guard min(1, L_gt/L_pred) — the analogue of g_area), Dr. Seg (count consistency). No RL paper rewarding against multiple annotators was found | https://arxiv.org/abs/2506.22624, https://arxiv.org/abs/2505.22596, https://arxiv.org/abs/2504.07615, https://arxiv.org/abs/2603.00152 |
| Hungarian | `linear_sum_assignment` (modified Jonker–Volgenant, Crouse 2016), rectangular, `maximize=True` | https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html |
| Threads | OpenCV releases the GIL in every wrapped call (ERRWRAP2); `setNumThreads` is process-global and on the macOS GCD backend only values ≤ 0 take effect (`setNumThreads(1)` is a no-op there, `getNumThreads()` stays 8 — verified locally); numpy releases the GIL only for loops > 500 elements | https://github.com/opencv/opencv/blob/4.x/modules/python/src2/cv2_util.hpp, https://github.com/opencv/opencv/blob/4.x/modules/core/include/opencv2/core/utility.hpp |

## Sandbox (scratch `impl/M6/`)

- **Lab differential** (`diff_lab.py`): on 251 fixture items × 8 rollouts, production τ, T5 and every per-object
  term (IoU_tc, F5, g_area, f, signature) are bit-identical to `lab/t_reward/rewardcore`; v4 advantages agree to
  1.0e-13 (the binomial bias sum). Cache build 3.0 ms per item on average for both.
- **Scorer variants** (`opt_score.py`, `opt_iou.py`, µs per object, loaded machine): rater intersections from a
  per-pixel rater bit-code + bincount: RIGA −18 %, LIDC +27 %; packed bits + `np.bitwise_count`: RIGA 56 → 33 µs,
  LIDC 6 → 5 µs (≈ 12 ms per all-RIGA batch); a gather-based band: slower. None adopted (gains ≤ 10 % of the
  batch time for extra code).
- **Execution** (`bench_k3.py`, `bench_k3b.py`, load average 5–14): per-rollout work is ~10 small numpy/OpenCV
  calls, so Python threads contend for the GIL: 4–8 threads were slower than serial on mixed batches
  (serial 54–65 ms vs 8 threads 65–91 ms) and helped only RIGA-only batches when OpenCV was made sequential
  (`setNumThreads(0)`: 185 → 139 ms). A 3-process pool helped only when every worker held all caches (lab style;
  47 vs 65 ms mixed) and lost once caches travel with the jobs (RIGA-only 392 vs 258 ms): a RIGA cache is ~3 MB
  derived state. Decision: serial by default, any `Executor` accepted; caches pickle as masks + τ + T5 and
  rederive exactly.

## Tests

`tests/test_reward.py` (no data): frozen constants; k = 1 refused; cache vs the independent reference on four
synthetic stacks (incl. a border-touching one); τ and T5 floors for identical raters; S_m vs Monte Carlo and the reference integral; EDT exactness and
the no-background case; in-band masks score F5 = 1; save/load round trip, tamper → ValueError, stale fingerprint
→ rebuild; pickle rederivation; cache digest equal in a fresh `-I` interpreter; per-object terms vs the reference
(100 random masks incl. noise and empty); IoU_tc same-topology / vote fallback; sub-filter speck invariance and
filter-size speck; area guard values; band monotone in excursion; shape check; completion contract keys;
truncated / incomplete / syntax / key-order failures; illegal bow-tie (terms ignored by the gate); empty object
list; joint topology over disc + cup; completion terms vs the reference; `score_batch` order and thread
equality; advantages vs the reference (300 Hypothesis batches, ≤ 1e-9); low-spread drop; group independence (no
batch z) and centring; debiased r_topo favours no class at π = f (exact enumeration); `loo_bias` vs scipy and the
lab values; valid > invalid at every reachable state (r_topo and advantage); dark room (1e-3 noise on a constant
term moves advantages ≤ 0.006, legality ≤ 4e-4); K10 (serial, 8 threads, 2 processes, fresh interpreter: equal
SHA-256).

`tests/test_reward_gates.py` (`@pytest.mark.data`): cache/terms vs reference on real items; cache store round
trip; K9 on the 251 items (CI low ≥ 95 %) and on LIDC fixtures + composites; local gradient on average (each move's
mean ΔA CI < 0); the per-rollout ≥ 95 % gate as a strict xfail (it fails, below); R4 on the 251 items (real and
1-px; half-width ≤ 0.10 and |mean| ≤ 0.07); K8 LIDC strata, K8b; slow: R3 reachable hacks on every fixture
item, K8 RIGA strata and R5 on all RIGA structures.

**Test-first caveat.** The production modules were written as the sandbox prototype before the repo tests.
To show the tests are not vacuous, `scratchpad/impl/M6/mutate.py` applies 13 single mutations to
`src/objdel/reward` and runs `tests/test_reward.py`: all 13 mutations are killed (debias, clip, σ floors, legality gate, low-spread drop, IoU gate on validity, IoU_tc fallback, deflation guard, joint topology, all-objects legality, rater-count correction, strict vote, τ floor). The τ-floor mutation survived the first run; `test_identical_raters_hit_the_tau_and_t5_floors` was added for it.

## The two gates that do not pass cleanly

**Per-rollout local gradient (M6 gate, target ≥ 95 %): 81.6 % on the 251 items, 84.1 % on the fixtures.**
The mean effect is strongly negative for every move (ΔA −0.25 to −1.76, CIs exclude 0), so there is no
systematic inflation or deflation incentive, as the lab found. Per rollout, however, an honest rollout is a
rater mask with a smooth perturbation of up to τ; when it is smaller than the raters, dilating it moves it
toward them and its advantage rises. Split by direction: moves away from the median rater area lower the
advantage in 97.6 % (251 items) and 98.3 % (fixtures) of cases; moves toward it only in 64–69 %. The original
≥ 95 % per-rollout criterion cannot hold for any reward that rewards closeness to the raters. Kept in the
suite as a strict `xfail`. **Proposal for the coordinator (a re-registration after measurement, so flagged as
such):** gate the per-rollout rate on moves away from the rater consensus (measured 97.6 %).

**R4 realistic neutrality: borderline.** The plan's lab value was +0.038 [−0.011, +0.091]. Across six independent
draws of the honest pools the mean is +0.030 to +0.066 (pooled ≈ +0.043) with half-width ≈ 0.045; the CI excludes
0 in 2 of 6 draws (+0.056 [+0.014, +0.100] in the gate run). The 1-px mode passes in every draw (+0.007 to +0.049).
So v4 carries a small pro-minority lean (≈ 0.04 advantage units, ~2 % of c_topo) under realistic rollouts. The
bootstrap CI resamples items only, not the pool construction, which is why single draws disagree.
`test_r4_neutrality_on_the_lidc_disagreement_items` (239 items since F2) pins half-width ≤ 0.10 and |mean| ≤ 0.07 on one fixed draw. No
definition change is proposed (the frozen constants were chosen on this set); the lean is reported for the
coordinator and for M8's per-image TV monitoring.

## Proposals (not adopted: the definition is frozen)

1. Low-spread drop: decide "one topology answer" from the gated r_topo spread instead of `sig_id` equality. Today a
   group with one signature whose IoU_tc straddles 0.5 within < 1e-2 is dropped although validity (r_topo 0 vs ≥ 0.1)
   differs. Rare; not measured on data.
2. r_topo estimator: Mohri et al. 2026 give minimax / U-statistic estimators for log rewards from small groups
   (https://arxiv.org/abs/2605.02435); the exact b(f, G) correction already removes the bias at π = f, so a switch
   would only matter away from the target.
3. Local-gradient KPI re-registration (above).

## Open issues

- **Phase 0 scripts are broken by the deletion:** `experiments/phase0/{riga_outlines, riga_census, baseline_vlm,
  baseline_medsam}.py` import `objdel.reward.AnnotatorSet` (frozen records; coordinator decision). *(Decided:
  frozen record, run against commit `faaacc4`; `experiments/phase0/README.md`.)*
- `betti` users after the deletion: `checker/certify.py`, `experiments/phase1/census_sizerule.py`,
  `experiments/phase0/lidc_census.py`, `tests/test_checker.py`, `tests/test_topology.py`. Not used by the reward.
  (Also `experiments/phase1/bench/xarch_hashes.py`, M9.)
- **OpenCV threads:** `cv2.setNumThreads(1)` is a no-op on the macOS GCD backend (`getNumThreads()` stays 8);
  `setNumThreads(0)` makes OpenCV sequential on every backend. M3's note and M7 should use 0.
- **Threads do not scale on the Mac** (GIL-bound small calls); serial meets K3 with margin here. M9 should time
  `score_batch` serially and with `ThreadPoolExecutor(8)` on the x86 box.
- R3 (b) holes stamped stays above the re-registered limit (62.7 % vs ≤ 60 %), unchanged from the lab.
- The in-band speck's CI lower bound is 94.9 % (point 96.4 %).
- R5 at m = 5 is 0.980, the upper edge (lab 0.97996).
- Slow tests (`-m slow`): R3 on every fixture item, R5 and RIGA K8 strata on all RIGA structures (subsamples are too
  noisy for these gates: RIGA cup ρ reached −0.14 on 1/8 of the images).

## Gate tables (written by `experiments/phase1/reward_gates.py`; 941 s with 3 processes)


<!-- gates:begin -->
Fixture items with ≥ 2 raters: 1002; full-LIDC train disagreement items: 251.

### K8, K8b, R5 (band)

| stratum (leave one rater out) | n | Spearman ρ(F5, area) [95 % CI] | median F5 | K8 (abs ρ ≤ 0.10) |
|---|---|---|---|---|
| riga/disc | 1200 | +0.008 [-0.028, +0.046] | 0.974 | PASS |
| riga/cup | 1200 | -0.040 [-0.077, -0.002] | 0.987 | PASS |
| lidc k=3 | 423 | +0.012 [-0.022, +0.045] | 0.961 | PASS |
| lidc k=4 | 1248 | -0.023 [-0.048, +0.006] | 0.968 | PASS |

Pooled ρ (reported, not gated): +0.070. LIDC k = 2 items have no leave-one-out case (one remaining rater is not RL-eligible; production has no single-rater prior).

| K8b criterion | result | target | verdict |
|---|---|---|---|
| (i) F5 monotone in excursion (0.25–2 τ, dilate and erode), all 1002 items | 100.0 % | ≥ 98 % | PASS |
| (ii) 2τ-dilated vote below the median held-out rater (853 items, k ≥ 3) | 100.0 % (eroded: 99.2 %) | ≥ 95 % | PASS |
| (iii) held-out F5 median / p10 (4071 cases) | 0.975 / 0.824 | ≥ 0.90 / ≥ 0.70 | PASS |

| R5: m reference raters (RIGA) | median F5 of a held-out expert | p10 | verdict ([0.92, 0.98]) |
|---|---|---|---|
| 2 | 0.961 | 0.784 | PASS |
| 3 | 0.971 | 0.805 | PASS |
| 4 | 0.976 | 0.823 | PASS |
| 5 | 0.980 | 0.836 | PASS |

Single-rater items: `ExpertCache.build` refuses k = 1 (`rl_eligible`); tested.

### K9 (topology ranking; v4 gated ≥ 95 %)

| item set | pairs | phase0 | v2 | v4 | v4_tauP | v4_hole_tc | v4_hole_all |
|---|---|---|---|---|---|---|---|
| fixtures + 2-part composites | 8018 pairs, 1327 items | 97.0 % [96.4 %, 97.7 %]; IoU ≥ 0.5: 99.2 % | 99.4 % [99.2 %, 99.6 %]; IoU ≥ 0.5: 99.7 % | 99.9 % [99.8 %, 100.0 %]; IoU ≥ 0.5: 100.0 % | 99.9 % [99.9 %, 100.0 %]; IoU ≥ 0.5: 100.0 % | 100.0 % [99.9 %, 100.0 %]; IoU ≥ 0.5: 100.0 % | 100.0 % [99.9 %, 100.0 %]; IoU ≥ 0.5: 100.0 % |
| full-LIDC train disagreement | 789 pairs, 229 items | 90.7 % [87.8 %, 93.7 %]; IoU ≥ 0.5: 99.3 % | 96.2 % [94.5 %, 97.7 %]; IoU ≥ 0.5: 98.0 % | 98.9 % [98.0 %, 99.6 %]; IoU ≥ 0.5: 99.7 % | 98.9 % [97.9 %, 99.6 %]; IoU ≥ 0.5: 99.6 % | 99.1 % [98.4 %, 99.7 %]; IoU ≥ 0.5: 99.7 % | 99.1 % [98.4 %, 99.7 %]; IoU ≥ 0.5: 99.7 % |

### R3 adversaries (dilated union and in-band speck gated ≥ 95 %; others reported) and AB2/AB3/AB4


Honest rollouts: **real**. Cell = share of groups where the adversary's advantage is below the median honest advantage [95 % cluster-bootstrap CI]; hole rows: win rate against the median honest rollout of the same topology class (re-registered KPI ≤ 60 %, CI upper ≤ 65 %), then the original below-median share.

| adversary | items | phase0 | v2 | v4 | v4_tauP | v4_hole_tc | v4_hole_all |
|---|---|---|---|---|---|---|---|
| dilated_union | 1002 | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] |
| band_speck | 281 | 52.8 % [48.6 %, 56.8 %] | 94.9 % [93.4 %, 96.1 %] | 96.4 % [94.9 %, 97.5 %] | 96.4 % [95.0 %, 97.7 %] | 96.4 % [94.9 %, 97.6 %] | 96.4 % [94.9 %, 97.6 %] |
| hull | 1002 | 73.4 % [71.8 %, 75.0 %] | 51.1 % [49.0 %, 53.1 %] | 65.7 % [63.8 %, 67.5 %] | 51.0 % [49.1 %, 53.1 %] | 65.7 % [63.8 %, 67.5 %] | 65.7 % [63.9 %, 67.5 %] |
| union | 1002 | 55.9 % [54.1 %, 57.7 %] | 36.5 % [34.8 %, 38.4 %] | 49.8 % [47.6 %, 51.8 %] | 36.5 % [34.6 %, 38.4 %] | 49.8 % [47.8 %, 51.8 %] | 49.8 % [47.6 %, 51.9 %] |
| core | 1001 | 25.7 % [24.0 %, 27.4 %] | 40.2 % [38.2 %, 42.1 %] | 55.8 % [53.7 %, 58.0 %] | 40.2 % [38.1 %, 42.2 %] | 55.8 % [53.6 %, 57.8 %] | 55.8 % [53.7 %, 58.0 %] |
| holes_stamped | 101 | win 99.9 % [99.6 %, 100.0 %]; below 0.4 % | win 60.4 % [55.6 %, 65.2 %]; below 51.9 % | win 62.7 % [58.2 %, 67.4 %]; below 54.4 % | win 60.4 % [55.8 %, 64.9 %]; below 57.7 % | win 66.4 % [62.0 %, 70.8 %]; below 54.3 % | win 62.7 % [58.1 %, 67.1 %]; below 64.4 % |
| minority_hole | 70 | win 99.6 % [98.9 %, 100.0 %]; below 0.6 % | win 52.1 % [46.6 %, 57.5 %]; below 58.4 % | win 54.5 % [48.7 %, 60.1 %]; below 60.4 % | win 52.1 % [46.6 %, 57.7 %]; below 64.8 % | win 58.4 % [52.9 %, 63.9 %]; below 59.1 % | win 54.5 % [48.9 %, 59.9 %]; below 75.4 % |

Honest rollouts: **jit1px**. Cell = share of groups where the adversary's advantage is below the median honest advantage [95 % cluster-bootstrap CI]; hole rows: win rate against the median honest rollout of the same topology class (re-registered KPI ≤ 60 %, CI upper ≤ 65 %), then the original below-median share.

| adversary | items | phase0 | v2 | v4 | v4_tauP | v4_hole_tc | v4_hole_all |
|---|---|---|---|---|---|---|---|
| dilated_union | 1002 | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] | 100.0 % [100.0 %, 100.0 %] |
| band_speck | 281 | 67.3 % [63.6 %, 70.7 %] | 99.7 % [99.3 %, 99.9 %] | 99.9 % [99.8 %, 100.0 %] | 99.9 % [99.8 %, 100.0 %] | 99.9 % [99.8 %, 100.0 %] | 99.9 % [99.8 %, 100.0 %] |
| hull | 1002 | 78.8 % [76.8 %, 80.9 %] | 60.0 % [57.4 %, 62.8 %] | 71.2 % [68.8 %, 73.4 %] | 58.7 % [56.0 %, 61.4 %] | 71.2 % [68.7 %, 73.5 %] | 71.2 % [68.8 %, 73.5 %] |
| union | 1002 | 51.8 % [49.3 %, 54.0 %] | 38.8 % [36.1 %, 41.7 %] | 49.6 % [46.8 %, 52.4 %] | 37.9 % [35.0 %, 40.7 %] | 49.6 % [46.8 %, 52.4 %] | 49.6 % [46.8 %, 52.4 %] |
| core | 1001 | 26.2 % [24.2 %, 28.3 %] | 46.5 % [43.6 %, 49.2 %] | 56.7 % [54.0 %, 59.5 %] | 45.7 % [42.7 %, 48.8 %] | 56.7 % [53.9 %, 59.6 %] | 56.7 % [53.9 %, 59.3 %] |
| holes_stamped | 101 | win 100.0 % [100.0 %, 100.0 %]; below 0.1 % | win 55.2 % [47.8 %, 62.7 %]; below 55.3 % | win 63.2 % [56.4 %, 69.9 %]; below 57.3 % | win 55.2 % [47.7 %, 62.6 %]; below 61.5 % | win 73.5 % [67.8 %, 79.2 %]; below 55.1 % | win 63.2 % [56.6 %, 69.8 %]; below 67.5 % |
| minority_hole | 70 | win 100.0 % [100.0 %, 100.0 %]; below 0.2 % | win 49.9 % [40.6 %, 59.6 %]; below 61.1 % | win 58.9 % [49.9 %, 67.6 %]; below 65.9 % | win 49.9 % [40.3 %, 59.5 %]; below 70.3 % | win 71.2 % [63.6 %, 78.4 %]; below 61.9 % | win 58.9 % [50.1 %, 68.1 %]; below 79.8 % |

| R4 neutrality (minority − mode advantage), 251 items, 300 batches | phase0 | v2 | v4 | v4_tauP | v4_hole_tc | v4_hole_all |
|---|---|---|---|---|---|---|
| exact | -0.563 [-0.725, -0.401] reject | +0.056 [+0.049, +0.063] reject | -0.009 [-0.015, -0.004] reject | -0.015 [-0.021, -0.009] reject | -0.004 [-0.013, +0.005] ACCEPT | -0.294 [-0.345, -0.243] reject |
| jit1px | -0.458 [-0.559, -0.349] reject | +0.014 [-0.062, +0.090] ACCEPT | +0.007 [-0.047, +0.061] ACCEPT | -0.028 [-0.064, +0.008] ACCEPT | -0.130 [-0.198, -0.064] reject | -0.275 [-0.355, -0.195] reject |
| real | -0.178 [-0.240, -0.118] reject | +0.111 [+0.062, +0.161] reject | +0.056 [+0.014, +0.100] reject | +0.048 [+0.007, +0.089] reject | -0.054 [-0.104, -0.003] reject | -0.225 [-0.292, -0.160] reject |

### Local gradient (dilating or eroding one honest rollout, others fixed)

| local gradient, 251 full-LIDC items × 20 groups | moved rollouts | ΔA < 0 (per rollout) | mean ΔA [95 % CI over items] |
|---|---|---|---|
| dilate 1px | 4869 | 82.5 % | -0.522 [-0.591, -0.461] |
| dilate 0.25tau | 3573 | 83.4 % | -0.662 [-0.739, -0.586] |
| dilate 0.5tau | 4733 | 87.9 % | -1.168 [-1.243, -1.094] |
| erode 1px | 4869 | 77.0 % | -0.848 [-0.937, -0.768] |
| erode 0.25tau | 3573 | 72.7 % | -0.972 [-1.060, -0.884] |
| erode 0.5tau | 4733 | 84.5 % | -1.761 [-1.842, -1.685] |
| all moves | 26350 | **81.6 %** (target ≥ 95 %) | |
| moves away from the median rater area (dilating a rollout ≥ it, eroding one ≤ it) | 13728 | 97.6 % | |
| moves toward the median rater area | 12622 | 64.2 % | |

| local gradient, fixture items × 6 groups | moved rollouts | ΔA < 0 (per rollout) | mean ΔA [95 % CI over items] |
|---|---|---|---|
| dilate 1px | 6002 | 78.4 % | -0.456 [-0.492, -0.421] |
| dilate 0.25tau | 2912 | 81.8 % | -0.245 [-0.267, -0.223] |
| dilate 0.5tau | 4826 | 86.4 % | -0.666 [-0.697, -0.636] |
| erode 1px | 6002 | 84.9 % | -0.777 [-0.828, -0.728] |
| erode 0.25tau | 2912 | 82.7 % | -0.390 [-0.421, -0.360] |
| erode 0.5tau | 4826 | 90.3 % | -1.106 [-1.144, -1.068] |
| all moves | 27480 | **84.1 %** (target ≥ 95 %) | |
| moves away from the median rater area (dilating a rollout ≥ it, eroding one ≤ it) | 14083 | 98.3 % | |
| moves toward the median rater area | 13397 | 69.2 % | |

### AB4: leave-one-rater-out g_area (held-out rater vs the other raters' area range)

| item set | held-out cases | zero tolerance: median / p10 / share < 1 | τ·P tolerance: median / p10 / share < 1 |
|---|---|---|---|
| RIGA fixtures | 2400 | 1.000 / 0.895 / 33.3% | 1.000 / 1.000 / 0.4% |
| LIDC fixtures | 1671 | 0.991 / 0.757 / 51.8% | 1.000 / 1.000 / 2.0% |
| full-LIDC disagreement | 732 | 0.971 / 0.637 / 55.9% | 1.000 / 1.000 / 3.0% |

Switch rule (median zero-tolerance g_area < 0.9 → τ·P): medians 1.000, 0.991, 0.971 → **keep zero tolerance**.

### K3, K4, K10 (completion text → reward → advantages)

| K3 / K4 workload (256 completions = 32 prompts × 8) | executor | batch wall median / max (ms) | per-completion p50 / p99 (ms) |
|---|---|---|---|
| mixed (40 % RIGA disc + cup, 60 % LIDC) | serial | 48.6 / 51.2 | 0.06 / 0.53 |
| mixed (40 % RIGA disc + cup, 60 % LIDC) | 8 threads | 55.3 / 59.1 | 0.06 / 0.53 |
| RIGA only | serial | 103.5 / 127.7 | 0.39 / 0.60 |
| RIGA only | 8 threads | 121.7 / 165.2 | 0.39 / 0.60 |

Load average during the run: 9.1 (8-core Apple M3, shared with other agents).

K10: SHA-256 of 256 completions' terms + advantages identical across serial run 1, serial run 2, 8 threads, 3 processes, OpenCV 8 threads, fresh interpreter: **PASS** (62c3ea691b426131…).

Total run time 941 s, 3 worker processes.
<!-- gates:end -->
