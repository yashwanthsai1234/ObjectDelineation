# ADR-0002: Score against the set of rater outlines, not one mask

**Status:** Accepted
**Date:** 2026-10-19 (commit `1c54bc4`, reward v4 in `src/objdel/reward`; census 2026-10-18, commit `29868b9`)
**Deciders:** Repository owner

## Context

LIDC-IDRI has up to four raters per nodule. RIGA has six per photograph. The raters disagree.
The Phase 0 census measured how often they disagree on topology (`docs/research/records/phase0/RESULTS.md`):

- 114 of 1,877 nodules (6.1 %) have a topology disagreement after the size filter.
- 65 nodules have a hole in at least one rater's outline. In all 65, at least one other rater drew no hole.
- 38 of those 65 were outlined by all four raters.

A single consensus mask would erase this. A model trained against it could never be right about a hole that half the raters drew.

## Decision

Reward v4 treats the raters as a set of valid answers:

1. **IoU_tc** is the best IoU against the raters who drew the same topology as the answer.
   When no rater did, it is the IoU against the strict-majority vote.
2. **F5** scores the boundary against a band between the core (all raters) and the union (any rater). Pixels inside the band are neither right nor wrong.
3. **topology** is rewarded by rater frequency. An answer with a topology that some rater drew is valid. A common topology is worth more than a rare one, in proportion.
4. **g_area** guards against inflation and deflation with the smallest and largest rater areas.

Full formulas are in the [reward reference](../../reference/reward.md).

## Options considered

### Option A: Set-valued scoring (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | High. Needs an expert table per object |
| Cost | Tables are built once on the CPU and stay on the device |
| Correctness | Respects disagreement. No rater is wrong by definition |

**Pros:** a hole that any rater drew counts. The reward matches the rater histogram.
**Cons:** raters must be independent draws. Copied outlines count twice, so the LIDC parser merges identical masks (rule P5).

### Option B: One consensus mask (strict majority)

**Pros:** simple. One IoU.
**Cons:** erases minority topologies. A rater-drawn hole would always cost reward.

### Option C: Mean IoU against each rater

The lab tested this as a variant. Its result files carry the label "meanIoU" (`archive/phase1_lab/t_reward/results/iter4_adversaries_meaniou.txt`).
The frozen reward uses IoU_tc instead. The lab's own reasons are in `archive/phase1_lab/t_reward/results/ITERATION4.md`. This record does not restate them.

## Trade-off analysis

Option A costs one expert table per object. The table is the price of exactness: it holds the band, the vote and the signatures.
In exchange the reward has no wrong answer inside the rater set.
The topology term is only mean-centred, never divided by its group standard deviation, so one outlier cannot dominate the geometric terms (lab ITERATION3).

## Consequences

- Easier: the model can learn holes and multiple parts without being punished for a rater's choice.
- Easier: the evaluation uses the same terms as training, with no second metric definition.
- Harder: items with fewer than two raters cannot enter RL (`rl_eligible`). Their band has no width.
- Harder: a change of any `REWARD_V4` value invalidates every stored expert table (`spec.fingerprint`).
- Revisit: abstaining raters (all components below the size filter) leave the topology vote.
  They stay in the band, the vote and the area range.

## Action items

1. [x] Freeze the constants (`REWARD_V4`).
2. [x] Store tables with a SHA-256 and the reward fingerprint (`ExpertCache.save`).
3. [ ] Test the reward on RIGA RL. The RL path reads LIDC only today.
