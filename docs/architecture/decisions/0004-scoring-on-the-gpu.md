# ADR-0004: Score on the GPU. Masks never return to the host.

**Status:** Accepted
**Date:** 2026-11-07 (commit `9c5dd2e`, "Scoring on the GPU"; first GPU checker 2026-11-04, commit `b905e7e`)
**Deciders:** Repository owner

## Context

After ADR-0001 the checker drew masks on the GPU.
The reward then copied masks to the host and computed IoU, boundary score, area and topology there.
That put a device-to-host copy and a CPU pool between every step's answers and their rewards.
The owner's rule is that the checker stays on the GPU (`docs/research/notes/gpu_scoring.md`).

## Decision

Every check and every score runs in the CUDA checker: polygon, mask, defects, IoU_tc, F5, g_area and the size-filtered signature.

- The expert tables of each item upload once (`ExpertCache.device()` calls `add_tables`) and stay resident until the cache is garbage-collected.
- One call, `score_polygons`, scores all parsed objects of a batch against their tables.
- Only per-object numbers return: defects, area, three float64 terms, the topology frequency and the signature.
- The host keeps only the combination of per-object numbers into per-answer terms (mean over objects, joint topology, `sig_id`).
- The scorer is bit-identical to the host scorer it replaced. The replaced scorer lives on as a test oracle (`tests/oracles/score_numpy.py`).

## Options considered

### Option A: Score on the device (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | High. Labelling and sums run in kernels |
| Cost | About 140 MB of resident tables for all 7,341 LIDC train items |
| Latency | Rewards cost 0.01 s per step in the full run (`full_rl_2026-11-08.md`) |

### Option B: Draw masks on the GPU, score on the CPU

**Pros:** easy to write. **Cons:** copies a mask per answer to the host every step and needs a CPU pool sized to the box. This was the state before the decision.

### Option C: Score in PyTorch on the GPU

Not evaluated in this repository. The exactness argument below is stated for hand-written kernels with fixed summation order.

## Trade-off analysis

The exactness argument carries Option A. Counts are integers. Distance sums are exact multiples of 2^-23. They add per row in int64 and across rows in int128. So their float64 value equals any other summation order while the sum is below 2^30.
The remaining terms are the same IEEE operations in the same order as the Python code. The build forbids FMA contraction.
Option C cannot offer that guarantee.

## Consequences

- Easier: a reward call is one checker call. Scoring overlaps the trainer on a background thread because the call releases the GIL.
- Easier: no CPU scoring pool for training.
- Harder: the scorer needs a GPU to be tested on its target. The CPU build runs the same code for the tests.
- Harder: tables are tied to a frame. An object frame must equal its table frame, or the call raises.
- Not covered: exam metrics (pairwise IoU, NSD, GED, HM-IoU) still read masks on the host. They run a few times per run, not per step.

## Evidence

| Check | Result | Source |
|---|---|---|
| Equivalence with the CPU checker's recorded results on a 5090 | 0 mismatches over 144,133 LIDC polygons and 37,093 LIDC masks | `docs/research/notes/issues_solved_2026-11-07.md` section 8 |
| Sanitizers | memcheck and initcheck: 0 errors on 52 tests | same, and `docs/research/notes/full_rl_2026-11-08.md` |
| Step time | 56.6 s to 36.3 s (R2 to R3). Most of the gain is image opening, not scoring | `issues_solved_2026-11-07.md` |

## Action items

1. [x] `add_tables`, `score_polygons`, `score_masks` in the checker.
2. [x] Archive polycore and keep its results as goldens.
3. [ ] Move the exam metrics to the device if they start to cost wall time.
