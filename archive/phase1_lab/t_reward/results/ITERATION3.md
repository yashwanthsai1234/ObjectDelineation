# t_reward iteration 3: judge items R3–R6 and the bracket-clip addendum

The judge's verdict (`../JUDGE_VERDICT.md`) asked for R3–R6 and the addendum. All runs used ≤ 4 processes; worker comparisons in K3 and K10 are therefore 1 vs 4 this time.

- **Raw outputs:** `results/iter3_*.txt`.
- **New scripts:** `k9_v3.py`, `adversaries.py`, `hole_adversaries.py`, `hole_weights.py`, `guard_probe.py`, `k8_r5.py`, `k8b_feasibility.py`, `r4_scale.py`, `dynamics_clip.py`.
- **New tests:** `tests/test_iter3.py`, plus additions to `tests/test_reward_v2.py`.
- **Status:** 73 tests pass and ruff is clean.

## Frozen reward v3

Given a group of G = 8 rollouts for one item with k ≥ 2 raters:

```
min_area   = max(5 px, ceil(0.005·|U|))                    (one size filter: census, targets, reward, eval)
sig(x)     = sorted holes per component of x, ignoring components and holes smaller than min_area
IoU_tc     = max IoU over the raters with sig(rater) = sig(x); if no rater drew sig(x): IoU vs the strict-majority vote
F5         = 1 / (1 + (2D/P) / T5²)                        (global rater band, as in iteration 1)
valid      = [sig(x) ∈ rater sigs] · [IoU_tc ≥ 0.5]
r_topo     = valid · max(0.1, 1 + λ(log f(a) − b(f(a), G) − log π̂_loo(a)))
             π̂_loo = max(n(a)−1, 0.5)/(G−1)
             b(f, G) = log f − E_{N~Bin(G−1,f)} log(max(N, 0.5)/(G−1))
A          = Σ_k w_k · z_k  +  c_topo · (r_topo·1[legal] − mean_group)
             z_k = per-group z-score of term k ∈ {IoU_tc, F5, legal}, with σ_k ← max(σ_k, floor_k)
             no batch-level z-score; every term is zeroed when the rollout is illegal (hard gate)
```

Constants:

| constant | value |
|---|---|
| w (IoU_tc, F5, legal) | 0.3 / 0.3 / 0.1 |
| σ floors | 0.05 / 0.05 / 0.25 |
| λ | 0.5 |
| c_topo | 2 |
| G | 8 |
| ε (bracket clip) | 0.1 |
| s\* | 0.95 |
| τ floor | 1 px |

**Single-rater items (k = 1) are excluded from RL** (`rl_eligible`).

The topology term is FlowRL-style rater-frequency matching with a group leave-one-out estimator. It is an instance of MaD-RL's gated reverse-KL (Kulkarni et al. 2026) with the rater topology histogram as the target. Cite FlowRL, MaD-RL, Multi-Answer RL and fitness sharing.

Code locations:

- `rewardcore/reward.py`: `signature_min_area`, `topo_conditional_iou`, `rl_eligible`.
- `rewardcore/advantages.py`: `dm_loo_bias`, `dm_loo_topology`, `reward_v3`.
- `rewardcore/batch.py`: `advantages_v3`.
- `rewardcore/topology.py`: `filtered_signature`.

What changed since iteration 2:

1. The IoU_tc fallback for a topology no rater drew is now the vote (R3).
2. The size filter is added (R6).
3. The batch-level z-score of the geometric part is **removed**.
   - With it, the topology weight depended on batch composition; this is what made the band-speck adversary win.
   - Without it, the geometric advantage of one rollout is bounded by Σw·√(G−1) = 1.59, so c_topo = 2 dominates it structurally.
4. **The DM-loo term is debiased** by the exact binomial bias b(f, G) (R4).
5. The bracket is clipped at ε = 0.1 (addendum).

## R3: IoU_tc fallback + adversary suite

**K9 (`iter3_k9.txt`, `iter3_k9_full_lidc.txt`).**

- Pairs:
  - Defects are sized at and above the size filter (A and 2A px).
  - Groups mix in 0–3 extra wrong-topology fillers; 75 % of groups have at least one.
  - 95 % CIs are a cluster bootstrap over items.

| item set | IoU alone | frozen v3, c_topo = 1 / 1.5 / 2 / 3 | v1 |
|---|---|---|---|
| fixtures + 2-part composites (7,963 pairs, 1,328 items) | 52.0 % | **99.94 % [99.87, 99.99]** at every c_topo; IoU ≥ 0.5: **100.00 %** | 98.91 % |
| full-LIDC train disagreement items (789 pairs, 229 items) | 57.7 % | **98.86 % [98.03, 99.60]**; IoU ≥ 0.5: 99.57 % [98.98, 100] | 95.56 % |

- **Verdict: PASS.**
- With the vote fallback but before the batch-z removal, c_topo = 1 gave 95.48 % [94.68, 96.23] and c_topo = 2 gave 99.79 %. The vote fallback alone did not make K9 insensitive to c_topo; removing the batch z did.

**Adversaries (`iter3_adversaries_final.txt`, `iter3_hole_adversaries.txt`).**

- Each adversary sits in a group with 7 honest rollouts: jittered rater masks, with the size-filtered topology preserved.
- KPI: the adversary's advantage is below the median honest rollout in ≥ 95 % of groups, with the 95 % CI lower bound used as the pass criterion.
- Columns: **1 px** = honest rollouts are 1-px-jittered rater masks (as specified). **0.5 τ** = honest rollouts are noisier (stress test).

| adversary | items | frozen v3, 1 px | v3, 0.5 τ (stress) | v1, 1 px | verdict |
|---|---|---|---|---|---|
| union dilated by 2τ | 1,002 | **100 % [100, 100]** | 99.9 % | 100 % | PASS |
| detached speck (min_area px) inside the band | 281 | **100 % [99.9, 100]** | 94.0 % [91.6, 96.1] | 88.4 % | PASS (stress: borderline) |
| convex hull of the union | 1,002 | 60.2 % [57.4, 62.9]: disc 99.7 %, cup 87.9 %, **LIDC 37.9 %** | 12.2 % | 79.8 % | **FAIL** |
| vote with every rater's holes stamped | 108 | 53.1 % [46.9, 58.7] | 12.3 % | 5.0 % | **FAIL** |
| vote + minority rater's hole, below median honest | 70 | 72.3 % [65.7, 78.6] | 21.4 % | – | **FAIL** |
| vote + minority rater's hole, below that rater's own jittered mask | 70 | 52.7 % [44.1, 61.3] | 13.6 % | 0.1 % (combined) | **FAIL** |
| reversed-winding ring | – | skipped | – | – | n/a: a mask has no orientation; winding is a polygon-legality matter (t_polycore) |

The hull and hole adversaries fail. Here is why, and the guard result.

- **Area-inflation guard (added as instructed, then rejected).** g = 1/(1 + (x/0.3)²), with x = max(0, area − max rater area)/(T5·P_U), multiplies IoU_tc and F5.
  - It raises the hull to 86.3 % (1 px; LIDC 78 %) and 72.9 % (stress). It still fails.
  - It lowers the band speck under stress from 94.0 % to 88.9 %.
  - It hits real experts. The probe (`iter3_guard_probe.txt`) shows a held-out largest rater (k ≥ 3) exceeds the other raters' maximum area by x = 0.169 (median) and 0.63 (p90). The union's x is 0.128 (median) and 0.32 (p90).
  - A legitimately large expert is therefore at least as "inflated" as the union. The guard would cost such an expert a median factor of 0.76 on the geometric terms.
  - **No guard is in the frozen reward.**
  - The union/hull win is concentrated on tiny LIDC nodules. There, a 1-px jitter is about 0.5τ, so the honest rollouts are worse than the union, which is the rater egg.
  - Recommendation: monitor inflation in the eval hack audit (|M|/|E|, R10) rather than shape the reward.
- **Hole adversaries.**
  - These masks have consensus geometry plus a topology that some rater drew, so the reward treats them as valid.
  - The pure-geometry test (`iter3_hole_weights.txt`) compares the adversary with the minority rater's own jittered mask, where the DM term is identical. IoU_tc favours the honest mask in only 60.2 % of pairs (median −0.027); F5 always favours the adversary (+0.015).
  - **No IoU/band weighting reaches 95 %**: the scores are 53 % at 0.3/0.3, 57 % at 0.45/0.15 and 60 % with IoU only.
  - A topology-conditional band (each class banded by its own raters) raised this test only to 59.5 % and broke neutrality (below), so it is not used.
  - Under the frozen reward these outputs are about as good as a 1-px-jittered copy of the minority rater.
  - **FAIL as specified.** Making them lose would require penalising consensus geometry, which conflicts with the band's purpose. The coordinator should decide whether this adversary describes a defect or a legitimate answer.

## R5: K8 per stratum, K8b, rater-count correction, single-rater policy

`iter3_r5.txt`, `iter3_k8b_feasibility.txt`.

**K8: |Spearman ρ(F5, area)| ≤ 0.10 per stratum.** Pooled ρ = +0.107 (reported, not gated).

| stratum | n | F5 ρ | verdict |
|---|---|---|---|
| riga/disc | 1,200 | +0.008 | PASS |
| riga/cup | 1,200 | −0.040 | PASS |
| lidc k=3 | 423 | +0.012 | PASS |
| lidc k=4 | 1,248 | −0.023 | PASS |
| lidc k=2, single reference, prior fitted on all items | 298 | +0.129 | FAIL |
| lidc k=2, LIDC-only σ prior, 2-fold cross-fitted (c ≈ 0.48, α ≈ 0.40) | 298 | **−0.063** | PASS |

- **Fix for k = 2:** the prior is now fitted on the item's own protocol. This is an honest correction (the pooled prior mixed fundus and CT), not a target change.
- **Single-rater policy:** k = 1 items are **excluded from RL**. Their band needs a prior, and their topology set is a single answer. They stay usable for SFT and for eval (F4/F5 with the per-dataset prior).
- **Tests:** a k = 1 cache builds bit-identically twice, after pickling, and in a fresh interpreter (`test_single_rater_cache_is_deterministic`).

**K8b: band discrimination.**

| criterion | target | measured | verdict |
|---|---|---|---|
| F5 non-increasing as \|excursion\| grows (0.25–2τ, dilate and erode) | ≥ 98 % of items | **100 %** | PASS |
| F5 of the vote dilated by 2τ ≤ 0.50 | ≥ 95 % of items | 89.0 % (RIGA 100 %, LIDC 81.7 %); median 0.231 | **FAIL** |
| F5 of every held-out rater ≥ 0.85 | ≥ 95 % of items | **48.0 %** (median of the item minimum 0.845); per held-out rater 87.8 % | **FAIL** |

**Why K8b cannot be passed honestly (feasibility analysis).**

- In 2.1 % of items, the worst held-out rater is farther out than the 2τ-dilated vote. No monotone band can pass those.
- **Global calibration.** Choosing one tolerance scale for all items (equivalent to choosing s\*) and changing the shape gives at best:

  | shape | best share of items |
  |---|---|
  | 1/(1+x) | 52.8 % |
  | exp(−x) | 59.0 % |
  | 1/(1+x²) | 71.6 % |
  | 1/(1+x³) | 77.4 % |

- **Per-item oracle.** Even an oracle tolerance set separately for every item caps at:

  | shape | share of items |
  |---|---|
  | 1/(1+x) | 84.4 % |
  | 1/(1+x²) | 94.8 % |
  | 1/(1+x³) | 96.4 % |

  Only near-step shapes reach 95 %, and a step-like band is poor for GRPO.
- **Root cause:** real outlier raters exist. Their excursion beyond the other raters' band is often within a factor of about 2.3 (5th percentile) of a 2τ consensus dilation.
- F5 is left unchanged; there is no target change and no step-like reward.
- **K8b FAIL as specified.** I recommend re-registering it with a per-rater criterion (87.8 % of held-out raters are ≥ 0.85) or a percentile criterion. That decision is the coordinator's.

**S_m correction test.** RIGA, m random reference raters plus 1 held-out expert, 3 draws per image.

| m | median F5 | disc | cup | verdict |
|---|---|---|---|---|
| 2 | 0.961 | 0.957 | 0.967 | PASS |
| 3 | 0.971 | – | – | PASS |
| 4 | 0.976 | – | – | PASS |
| 5 | 0.97996 | 0.972 | 0.988 | PASS, at the upper edge |

The median is inside [0.92, 0.98] for every m. It drifts upward with m (+0.019 from m = 2 to 5), so the i.i.d. correction slightly over-tolerates large rater sets.

## R6: signature size filter

**One rule:** `min_area = max(5 px, ⌈0.5 %·|U|⌉)` components and holes count. This is 5 px on every LIDC item and about 100 px on a RIGA disc, the Phase 0 values. It is implemented once (`signature_min_area`, `filtered_signature`); the filter is exact, with filled holes merging islands.

**Census.**

- The filter changes 9 of 4,369 fixture rater masks.
- Fixture disagreement items drop from 14 to 5.
- Of the 18 extra parts and holes in LIDC fixture raters, 12 are 1–3 px annotation noise.
- Full LIDC under the filter: **251 train, 37 val and 87 test disagreement items**.

**Tests (`tests/test_iter3.py`):**

- A sub-filter speck or pinhole added to a rater mask leaves the item's topology set unchanged.
- A component or hole at the filter size changes it.
- `topo_any` and `sig` of a prediction are invariant to 50 random sub-filter specks and pinholes (still logged in `sig_raw`).
- A filter-size speck flips `topo_any` to 0.
- The filter equals objdel's `size_filter` followed by topology on 500 random masks; the size filter itself matches on 2,000 masks.

**Verdict: PASS.**

## R4: re-validation at scale (251 full-LIDC train disagreement items)

`iter3_r4.txt`.

**Per-topology gap**, production emulation, 239 items with both classes:

| reference | mean gap [95 % CI] |
|---|---|
| IoU vs vote (v1) | +0.130 [+0.108, +0.153] |
| **IoU_tc** | **−0.004 [−0.013, +0.005]** (median \|gap\| 0.026) |

**Neutrality (jittered):** minority minus mode advantage, 300 batches.

| setting | mean [95 % CI] | half-width | verdict |
|---|---|---|---|
| **v3, c_topo = 2, λ = 0.5** | **+0.009 [−0.026, +0.042]** | **0.034** | **ACCEPT** |
| v3, every point of c_topo ∈ {1, 1.5, 2, 3} × λ ∈ {0.25, 0.5, 1} | between +0.002 and +0.015 | ≤ 0.045 | all ACCEPT |
| without the debias | +0.077 [+0.044, +0.110] | – | reject |
| with the topology-conditional band | +0.038 [+0.006, +0.070] | – | reject |
| v1 | −0.450 [−0.590, −0.316] | – | reject |

**Neutrality (exact rater masks).**

| setting | mean [95 % CI] | note |
|---|---|---|
| v3 | −0.012 [−0.018, −0.006] | tiny; residual from group-composition effects on the group mean |
| without the debias | +0.056 | – |

- **The debias, b(f, G) exact.** Exhaustive enumeration over all groups of 8 drawn at π = f gives equal expected r_topo per class to 1e-9 (test).
- Values: b = 0.145 at f = 1/4, 0.088 at f = 1/2, 0.028 at f = 3/4, 0 at f = 1.
- **The class band is rejected.** It is the in-sample source of the +0.038: jittered copies sit inside their own class band, and a singleton class gets the larger m = 1 tolerance.
- **Choice: c_topo = 2, λ = 0.5.** All grid points pass neutrality, and K9 is flat in c_topo.
  - c_topo = 2 is the smallest value that dominates the geometric bound under stress: band speck 45 % at c_topo = 1 vs 94 % at c_topo = 2.
  - λ = 0.5 follows t_dynamics. At λ = 0.5 the ε clip never binds for k ≤ 4. At f = 1/6 (k = 6) the clip binds only when n = 8 (minimum bracket 0.068).

**c_topo and λ are now frozen.**

## Addendum: a valid rollout must never score below an invalid one

`iter3_dynamics_clip.txt`.

- The bracket is clipped at ε = 0.1. Before the debias, the minimum unclipped bracket over reachable states (k = 2..6, n = 1..8) was 0.552 at λ = 0.25, 0.104 at λ = 0.5 and −0.792 at λ = 1.
- **Property test (PASS):** `test_valid_always_beats_invalid_for_all_reachable_states`, for λ ∈ {0.25, 0.5, 1, 2}.
  - It covers every rater share j/k with k = 2..6 and every n(a) from 1 to 7.
  - The valid r_topo is strictly greater than the invalid one (0), and so is its mean-centred contribution.
- **Toy re-run** (t_dynamics `grpo_toy.py`, read-only; 20 seeds; suite plus a "hole 1/6" prompt).

  | variant (λ = 0.5) | TV, topology only | TV, IoU gap 0.1 | clip binding |
  |---|---|---|---|
  | original | 0.064 | 0.055 | – |
  | clipped | 0.064 | 0.055 | never |
  | clipped + debiased (frozen) | 0.060 | 0.061 | never |

  - Min-mode retention, clipped + debiased: 0.853 (topology only) and 0.795 (IoU gap 0.1).
  - At λ = 1 the clip binds on 1.6–3.1 % of valid rollouts, and TV stays at 0.079–0.083.

## Also re-checked with v3

- **Dark room** (`iter3_advantage_replay.txt`):
  - Noise on a near-constant term changes advantages by at most 0.004.
  - Ordering stays at 98.3–98.6 %; the remaining misses are pairs where the "worse" output really is closer to some rater.
  - Converged groups: median max|A| 0.051 (GRPO 1.92, GDPO 2.00).
- **K3/K4** (`iter3_k3.txt`, 4 workers, load average about 7):

  | workload | serial batch | 4 processes |
  |---|---|---|
  | mixed | 32.7 ms | 14.4 ms |
  | RIGA-only | 75.5 ms | 35.3 ms |

  Wall p99 ≤ 4.2 ms everywhere. **PASS.**
- **K10** (`iter3_k10.txt`): terms plus v3 advantages are bit-identical across serial ×2, 4 threads, 1 and 4 processes, OpenCV with 4 threads, and fresh interpreters. Signature ids use CRC32 instead of `hash`. **PASS.**

## Remaining risks

- **Three adversaries fail by design or by data** and need a decision on whether they describe defects:
  - the hull (LIDC);
  - vote + minority hole;
  - vote with all holes stamped.
- **K8b fails.** It needs re-registration; no honest F5 variant reaches it.
- **Neutrality was tested in-sample.** Samples are jittered copies of the raters; a new-expert test of singleton classes is impossible on this data.
- **The S_m correction drifts upward with m**, and m = 5 sits at the 0.98 edge.
- **The band speck under 0.5τ stress is 94.0 % [91.6, 96.1]**, borderline at c_topo = 2.
