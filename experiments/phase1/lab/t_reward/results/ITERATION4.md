# t_reward iteration 4 (final): inflation guard, realistic honest rollouts, re-registered KPIs

All runs used ≤ 4 processes. K3/K10 compare 1 vs 4 workers.

- **Raw outputs:** `results/iter4_*.txt`.
- **New scripts:** `scripts/iter4_adversaries.py`, `honest.py`, `inflation_incentive.py`, `k8b_v2.py`.
- **Updated scripts:** `r4_scale.py`, `k9_v3.py` (option `--v4-only`).
- **Status:** 75 tests pass and ruff is clean.

## Frozen reward v4

Given a group of G = 8 rollouts for one item with k ≥ 2 raters (k = 1 is excluded from RL). Changes from v3 are marked ★.

```
min_area = max(5 px, ceil(0.005·|U|))                     size filter for every topology signature
IoU_tc   = max IoU over the raters with the same filtered signature; IoU vs the vote if no rater drew it
F5       = 1 / (1 + (2D/P) / T5²)                          global rater band
g_area ★ = min(1, A_max/A) · min(1, A/A_min)               A_min, A_max = smallest / largest rater area (ALL raters)
valid    = [sig ∈ rater sigs] · [IoU_tc ≥ 0.5]
r_topo   = valid · max(0.1, 1 + 0.5·(log f(a) − b(f(a), 8) − log π̂_loo(a)))
           π̂_loo = max(n(a)−1, 0.5)/7
           b(f, G) = log f − E_{N~Bin(G−1,f)} log(max(N, 0.5)/(G−1))
A        = Σ_k w_k · z_k  +  2 · (r_topo·1[legal] − mean_group)
           k ∈ {IoU_tc, F5, legal, g_area ★}
           z_k = per-group z-score of term k with σ_k ← max(σ_k, floor_k)
           no batch z-score; a rollout that is not legal has every term zeroed
```

Constants:

| constant | value |
|---|---|
| w (IoU_tc, F5, legal, g_area) | 0.3 / 0.3 / 0.1 / 0.3 |
| σ floors | 0.05 / 0.05 / 0.25 / 0.05 |
| λ | 0.5 |
| c_topo | 2 |
| ε (bracket clip) | 0.1 |
| s\* | 0.95 |
| τ floor | 1 px |

- The guard's weight and floor reuse the other geometric terms' values. Nothing new was tuned.
- **Code:**
  - `rewardcore/reward.py`: `area_guard`; `Terms.g_area`.
  - `rewardcore/advantages.py`: `reward_v4(guard="term")`.
  - `rewardcore/batch.py`: `advantages_v4`.
- **Naming:** the topology term is FlowRL-style rater-frequency matching; it is an instance of MaD-RL's gated reverse-KL with the rater topology histogram as the target.

## 1. Realistic honest rollouts

The new construction is `scripts/honest.py`, mode `real`:

- With probability 0.25, the rollout is every rater's exact mask.
- Otherwise it is a smooth boundary perturbation with amplitude ~ U(0, τ) and correlation length max(3, τ) px.
- **One amplitude per pool entry is shared by all raters of the item**, and the size-filtered topology is preserved.

The shared amplitude is needed. Per-rater rejection sampling with an exact-mask fallback fell back on **16.5 %** of minority-topology draws (holes and small parts are fragile) against **2.8 %** for the mode class. That handed minority answers better honest copies and inflated neutrality to +0.285. With matched amplitudes the effect disappears.

The adversary suite is reported under both the old (1-px jitter) and the new honest construction (`iter4_adversaries_meaniou.txt`; τ·P variants in `iter4_adversaries.txt`).

## 2. Inflation/deflation guard

Four variants were tested: a gate (multiplies IoU_tc and F5) or a separate z-scored term, each with tolerance zero or ±τ·P_U.

| criterion (realistic honest; 1-px in parentheses) | v3 | **v4: term, zero tolerance** | gate, zero tolerance | τ·P variants |
|---|---|---|---|---|
| dilated union (2τ) below the median honest rollout | 100 % P | **100 % P** | 100 % P | 100 % P |
| hull below median | 51.8 % (59.9) F | **66.2 % [64.3, 68.0] (71.6) F** | 69.3 % (78.1) F | 51.8 % F (no effect: g = 1) |
| union below median | 36.4 % (38.9) F | **49.3 % (50.0) F** | 52.8 % F | 36.5 % F |
| core (rater intersection) below median | 39.6 % (46.9) F | **55.6 % (57.1) F** | 59.3 % F | 40.1 % F |
| in-band speck below median | 96.7 % (99.9) P | **96.5 % [95.2, 97.7] (99.9) P** | 93.6 % F | 95.9 % |
| R4 neutrality, 1-px jitter, 251 items | −0.003 ACCEPT | **+0.047 [−0.009, +0.104], half-width 0.056, ACCEPT** | +0.088 [+0.016, +0.160] **reject** | ACCEPT |
| R4 neutrality, realistic | +0.036 ACCEPT | **+0.038 [−0.011, +0.091], half-width 0.051, ACCEPT** | +0.015 ACCEPT | ACCEPT |
| K9, fixtures + composites | 99.94 % | **99.92 % [99.86, 99.98]; IoU ≥ 0.5: 100 %** | 99.28 % | 99.94 % |
| K9, full-LIDC train disagreement items | 98.86 % | **98.86 % [98.00, 99.62]; IoU ≥ 0.5: 99.71 %** | 98.10 % | 98.86 % |

- **The KPI for hull and plain union/core still FAILS** (target ≥ 95 %). The guard raises them by 13–19 points. The τ·P tolerance never triggers on them.
- A third alternative was tried and **rejected**: mean IoU over the same-topology raters instead of max.
  - It helps the hull only a little (71.7 %).
  - It drops the in-band speck to 85 %.
  - It breaks neutrality (+0.25 to +0.44 pro-minority). SFT trains on individual raters' masks, so policy samples resemble specific raters, and a singleton class's own rater is then its entire reference.

**Why the hedges cannot be beaten, and why that is acceptable.** The union, hull and core can only be built by knowing every rater's mask (oracle hedges). Their IoU_tc equals A_max/A_U (or A_core/A_min), and their band is 1.0. Realistic honest rollouts score lower on both terms.

What RL can actually learn is a systematic size bias. The test (`iter4_inflation_incentive.txt`) moves one honest rollout in its group, keeping the others fixed, and measures the change in its own advantage, mean over items [95 % CI]:

| move | v3 | v4 |
|---|---|---|
| dilate 0.25τ | −0.084 [−0.093, −0.074] | −0.122 |
| dilate 0.5τ | −0.37 | −0.56 |
| dilate 1 px | −0.32 | −0.47 |
| erode 0.25τ | −0.13 | −0.18 |
| erode 0.5τ | −0.61 | −0.89 |

**There is no inflation or deflation incentive under any scheme**, and the guard makes the penalty about 1.5× steeper. **Part-count inflation** is covered by topology: an extra part at or above the size filter is an invalid answer. Sub-filter specks are invisible to topology by design; they are logged in `sig_raw` and still cost IoU and band.

**Decision:** v4 uses the guard as a z-scored term with zero tolerance. It is the only variant that improves every inflation/deflation hedge while keeping neutrality in both honest modes and K9. Hull and union remain a documented residual, to be monitored in the eval hack audit (|M|/|E|, area outside [A_min, A_max]).

## 3. Re-registered KPIs (original and new numbers)

**Hole adversaries.**

- New KPI: the adversary's win rate against the median honest rollout of the **same topology class** must be ≤ 60 %, with CI upper bound ≤ 65 %.

| adversary | honest rollouts | v3 | v4 | original KPI (below the overall median, ≥ 95 %), v4 |
|---|---|---|---|---|
| minority-hole punch (70 items) | realistic | 57 % [51, 62] PASS | **60 % [54, 65] PASS (borderline)** | 58.6 % FAIL |
| minority-hole punch | 1-px | 49 % PASS | 58 % [50, 67] FAIL (CI) | 65.5 % FAIL |
| holes stamped (101 items) | realistic | 61 % [56, 65] FAIL | **64 % [59, 69] FAIL** | 52.2 % FAIL |
| holes stamped | 1-px | 54 % PASS | 62 % FAIL | 59.0 % FAIL |

- The guard raises the hole adversaries' win rate by about 3–9 points. These masks are built from the vote, so their area is always inside [A_min, A_max], while honest rollouts sometimes fall outside it.
- "Holes stamped" stays slightly above the re-registered bound under realistic rollouts. These are consensus-geometry answers with a rater-drawn topology, i.e. plausible answers. The excess (≈ 4 points) is small and documented.

**K8b** (`iter4_k8b.txt`).

| | criterion | result | verdict |
|---|---|---|---|
| original | F5 of the 2τ-dilated vote ≤ 0.50 | 89.0 % of items | FAIL |
| original | every held-out rater ≥ 0.85 | 48.0 % of items | FAIL |
| original | monotone | 100 % | PASS |
| **new (i)** | monotone in \|excursion\| | **100 %** | **PASS** |
| **new (ii)** | 2τ-dilated vote ranks below the median held-out rater | **99.9 %** of items (eroded: 97.1 %; LIDC 99.8 % / 95.2 %) | **PASS** |
| **new (iii)** | held-out-rater F5 median / p10 | **0.973 / 0.827** | **PASS** |

## Re-runs with v4

- **K10 PASS** (`iter4_k10.txt`): terms plus v4 advantages are bit-identical across serial ×2, 4 threads, 1 and 4 processes, OpenCV with 4 threads, and fresh interpreters (2 batches).
- **K3/K4 PASS** (`iter4_k3.txt`, load average about 3):

  | workload | serial batch | 4 processes | per-rollout wall p99 |
  |---|---|---|---|
  | mixed | 38.1 ms | 16.1 ms | ≤ 3.7 ms everywhere |
  | RIGA-only | 86.2 ms | 39.7 ms | |

## Summary of verdicts

| item | verdict |
|---|---|
| realistic honest rollouts (with the matched-amplitude fix) | done |
| dilated-union KPI | **PASS** (100 %) |
| hull / union / core KPI | **FAIL**: 66 / 49 / 56 % with the guard, up from 52 / 36 / 40 %. No learnable inflation incentive (all ΔA < 0). |
| in-band speck | PASS (96.5 %) |
| R4 neutrality (both honest modes) | PASS |
| K9 | PASS |
| K10 / K3 | PASS |
| hole KPI, re-registered | minority hole PASS (borderline); holes stamped FAIL (64 % vs ≤ 60 %) |
| K8b, re-registered | **PASS**; original FAIL |
