# t_reward iteration 2: anti-collapse topology term + topology-neutral geometry

This iteration adds the t_dynamics anti-collapse topology term to the reward core and re-tests every KPI.

The outputs are in `results/iter2/`. The new scripts are:

- `scripts/topo_gap.py`
- `scripts/topo_neutrality.py`
- `scripts/k9_v2.py`
- `scripts/k8_adjusted.py`

The updated scripts are `advantage_replay.py`, `k3_throughput.py` and `k10_determinism.py`. New tests are in `tests/test_reward_v2.py`, and `tests/test_determinism.py` now covers the v2 advantages.

## Reward v2

`R = 1[legal] · [w_iou·IoU_tc + w_band·F5 + w_topo·r_topo]`, combined at the advantage level as follows:

- **Geometric part, `[IoU_tc, F5, legal]`.**
  - Hard legality gate.
  - Per-term group z-score with σ floors: 0.05 for IoU and band, 0.25 for legality.
  - Weights 0.3 / 0.3 / 0.1.
  - Batch-level z-score of the weighted sum.
  - This is `gdpo_floor` on the gated terms.
- **IoU_tc.** The maximum IoU over the raters whose per-instance topology equals the sample's. If no rater has that topology, the maximum is taken over all raters.
- **Topology part.**
  - `r_topo = 1[valid]·(1 + λ(log f(a) − log π̂_loo(a)))`, with λ = 0.5.
  - `valid` = the topology matches a rater and IoU_tc ≥ 0.5.
  - `π̂_loo = max(n(a)−1, 0.5)/(G−1)`, where n(a) counts the group members with answer a, and G = 8.
  - The term is legality-gated, only **mean-centred** per group, and added with `c_topo = 2` in advantage units.
- **Code.**
  - `rewardcore/advantages.py`: `dm_loo_topology`, `reward_v2`.
  - `rewardcore/reward.py`: `topo_conditional_iou`.
  - `rewardcore/batch.py`: `advantages_v2`.

## 1. Real per-topology IoU gap (`results/iter2/topo_gap.txt`)

- **Where disagreement exists.** Raters disagree on per-instance topology in 14 of 602 LIDC items and 0 of 400 RIGA items.
- **Which minorities occur.** Every minority topology is a single rater: 7 one-hole items, 1 four-hole item, 2 two-hole items, and 4 items with an extra part.
- **Gap definition.** Gap = mean IoU of samples whose topology matches the vote's topology − the mean of those whose topology does not, averaged over items. CIs are a bootstrap over items.

| reference | rater masks, leave-own-out | production emulation (jittered rater masks vs all raters) |
|---|---|---|
| (c) vote, all raters (iteration 1) | **+0.154 [+0.066, +0.262]**, median +0.137, max +0.684 | **+0.127 [+0.057, +0.209]**, median \|gap\| 0.110 |
| (c') vote, other raters | +0.055 [+0.007, +0.107] | – |
| (a) consensus of same-topology raters | +0.079 [+0.020, +0.137] (singleton minority falls back to vote) | +0.018 [−0.085, +0.164], median \|gap\| 0.082, range up to +0.80 (unstable) |
| (b) max over same-topology raters | +0.067 [+0.015, +0.121] (same fallback) | **+0.007 [−0.014, +0.026]**, median \|gap\| 0.023 |

- **The gap is real, and above the ~0.1 the dynamics lab flagged.** The iteration-1 vote IoU shows a mean gap of +0.13 to +0.15.
- **(b) is the only topology-neutral reference, so v2 uses it.** Option (a) is unstable because a singleton class's consensus is that single rater's mask.
- **The leave-own-out numbers cannot rate (a) or (b).** With singleton minorities, removing the sample's own rater forces a fallback to the vote.

## 2. Reward v2 re-tests

### K9 (`results/iter2/k9_v2.txt`)

Each pair is placed inside a group of 8 with 6 topology-correct fillers, because the topology term depends on the group. Batches have 256 groups. The 12,077 pairs are the ones used in iteration 1.

| scheme | all | IoU ≥ 0.5 | pairs matched on IoU_tc (n = 7,087) |
|---|---|---|---|
| v1 (`gdpo_floor_gated`) | 98.72 % | 99.97 % | 99.90 % |
| v2, c_topo = 1 | 93.72 % **FAIL** | 94.62 % | 99.99 % |
| **v2, c_topo = 2** | **98.93 % PASS** | **99.36 %** | **100.00 %** |
| v2, c_topo = 3 | 99.64 % | 99.91 % | 100.00 % |

The c_topo = 1 failures come from how the pairs were matched:

- Pairs were matched on IoU vs the vote, but v2 scores IoU_tc.
- When the base is a rater's mask, "rater + speck" keeps IoU_tc ≈ 1. Its correct partner can be up to 0.19 lower (5th percentile).
- With c_topo = 2 the topology term outweighs that gap.

### Dark room (`results/iter2/advantage_replay_v2.txt`)

The topology term is mean-centred, never z-scored. As a result:

| scenario | v2 effect | ordering |
|---|---|---|
| D1: legality noise 1e-3 | max Δadvantage 0.000 | 98.4 % |
| D2: one output with legality 0.99 | 0.006 | 98.5 % |
| D3: topology noise 1e-3 | 0.001 | 98.4 % |
| D4: converged groups | median max\|A\| 0.13 (GRPO/GDPO: 1.9–2.0) | – |

Ordering on clean groups drops from 100 % to 98.4 %. Every miss is a group where the "known-worse" output (the consensus grown by 2 px or shifted) is closer to some individual rater than the consensus is: IoU_tc is higher in 9.0 % and 3.5 % of groups respectively. That is intended under any-rater semantics.

### Topology-class neutrality (`results/iter2/topo_neutrality.txt`)

Setup:

- Groups of 8 are drawn i.i.d. from each disagreement item's raters.
- There are 400 batches, each containing 14 groups.
- Minority vs mode advantage is computed per item, then averaged with items weighted equally.
- 95 % CIs are a bootstrap over items.

| samples | v1 (iteration 1) | v2, c_topo = 2 | v2 with vote-IoU |
|---|---|---|---|
| exact rater masks | **−0.512 [−1.314, +0.395]**; minority loses in 10/14 items | **+0.082 [+0.063, +0.103]**; minority loses in 0/14 items | −0.579 |
| jittered rater masks | **−0.597 [−1.195, +0.032]**; per-class means: mode +0.141, minority −0.456 [−0.881, −0.016] | **−0.144 [−0.504, +0.167]**, not significant; 7/14 items | −0.653 [−1.255, −0.004] |

- **Exact rater masks.** For v2, IoU_tc and band are constant (each sample matches itself), so only the DM-loo term varies. It leans slightly towards rare answers: c_topo = 1 gives +0.041, c_topo = 3 gives +0.123.
- **Jittered masks.** The residual splits into IoU_tc −0.18 and band −0.12, both not significant.
  - The IoU_tc part is a class-size effect: the max is taken over 3 raters for the mode class but over 1 for a minority.
  - The band part is consistent with minority raters being outliers at the band edge.
- **Both fixes are needed.** v2 with vote-IoU is as biased as v1.

### K3 / K4 (`results/iter2/k3_v2.txt`)

These runs include IoU_tc and the v2 group advantages. Load average was about 18 (iteration 1 ran at 60–190, which inflated wall time roughly 15×).

| workload | backend | batch wall | per-rollout wall p99 |
|---|---|---|---|
| mixed | serial | **22.5 ms** | 0.50 ms |
| mixed | 8 threads | 28.9 ms | 2.95 ms |
| mixed | 8 processes | **10.7 ms** | 0.76 ms |
| RIGA-only | serial | **53.2 ms** | 0.96 ms |
| RIGA-only | 8 threads | 46.1 ms | 4.46 ms |
| RIGA-only | 8 processes | **30.2 ms** | 1.67 ms |

- **K3 PASS:** every batch is ≤ 250 ms.
- **K4 PASS:** wall p99 ≤ 5 ms in every configuration.
- **Cache build per prompt:** median 0.6 ms (mixed) and 11.7 ms (RIGA); p99 88 ms.

### K10 (`results/iter2/k10_v2.txt`)

- **PASS.** The SHA-256 of per-rollout terms plus v2 advantages is identical across 8 settings, for 2 batches of 256.
- The settings are serial ×2, 8 threads, 1 process, 8 processes, OpenCV with 8 threads, and fresh interpreters.
- Topology ids use tuple-of-int hashes, which are not salted.

## 3. K8 adjusted for rater count (`results/iter2/k8_adjusted.txt`)

| formula | pooled ρ | partial ρ given k [95 % CI] | stratified ρ [95 % CI] | per-stratum ρ (k = 2 / 3 / 4 / 6) | median | verdict |
|---|---|---|---|---|---|---|
| F0 | +0.531 | +0.105 | +0.141 | +0.18 / +0.07 / +0.06 / +0.19 | 0.970 | FAIL |
| F1 | +0.333 | −0.122 | −0.077 | | 0.922 | FAIL |
| F2 / F3 | +0.391 | −0.096 | −0.053 | | 0.973 | PASS |
| F4 | +0.209 | −0.075 | −0.053 | | 0.967 | PASS |
| **F5** | +0.107 | **−0.089 [−0.108, −0.069]** | **−0.074 [−0.095, −0.053]** | +0.13 / +0.01 / −0.02 / −0.14 | 0.973 | **PASS** |

- **F5 passes K8 under the rater-count-adjusted definition.** The upper end of the partial-ρ CI (0.108) touches the limit.
- **The k = 6 stratum mixes RIGA discs and cups.** Separately, they give ρ = +0.01 and −0.04.
- **The adjusted ρ is less discriminating.** F2–F4 pass it too. F5 remains the best formula on pooled ρ and on the k = 3 and k = 4 strata.

## Recommendation and remaining risks

**Recommendation.** Adopt reward v2 with c_topo = 2, λ = 0.5, G = 8, IoU_tc and band F5. Report K8 as rater-count-adjusted ρ.

**Remaining risks.**

- **Thin evidence on topology disagreement.** Only 14 items (all LIDC) have it, and every minority class is a single rater.
- **v2 changes IoU semantics** from "agree with the consensus" to "match the closest rater of your topology". This is more permissive geometrically. The band term (F5) remains the consensus-anchored check.
- **c_topo trades two things off.** It balances K9 (needs ≥ 2) against a pro-rare lean on exact rater masks (+0.04 per unit of c_topo).
- **The topology term needs real groups.** It was validated in toy dynamics (t_dynamics) and in these replays, but not in a GRPO run.
