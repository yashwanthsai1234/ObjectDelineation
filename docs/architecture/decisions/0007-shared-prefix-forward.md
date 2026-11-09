# ADR-0007: Run each prompt once per group of 8 answers in the training forwards

**Status:** Accepted
**Date:** 2026-11-08 (commit `be52b1c`)
**Deciders:** Repository owner

## Context

TRL scores each answer as its own row, `[prompt | answer]`.
The prompt has about 325 real tokens, mostly image tokens. The 8 answers of a group share it.
So the prompt ran 8 times per group in every forward and backward: the reference pass and the update.
After earlier speed work, the R4 profile showed the cost in the update (12.46 s of 29.94 s).

## Decision

Pack the group into one row (`FastForward.packed_chunk`, enabled by `train_rl --shared-prefix`):

```text
[prompt minus its last token] [p_last a1 ...] [p_last a2 ...] ... [p_last a8 ...]
```

- The last prompt token repeats at the head of each answer block, with its own position. Its hidden state equals the state in TRL's row.
- A 4-D boolean mask lets each block attend to the prompt and to itself, never to another answer.
- Positions are each row's own 3-D rope positions from `RopeCache`.
- The vision tower gets one image per group.
- Update micro-batches keep whole groups together (`by_group_length`).
- A cost model decides per chunk (`packing_wins`): cost of n tokens = n + n^2 / 12,000. SDPA's masked kernel skips no tiles, so packed attention is dense. With 8 answers at the 1,536-token cap packing would cost 1.5x, so it falls back to TRL's rows.
- Off with `--compile`: packed widths are not bucketed and would recompile.

## Options considered

### Option A: Shared-prefix packing (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | High. A custom mask, positions and scatter-back |
| Cost | 0.36 of the tokens at the full run's sizes |
| Correctness | Log-probabilities equal TRL's per-row forward (atol 1e-5), gradients equal (atol 1e-4) on a tiny Qwen3-VL |

### Option B: TRL's per-row forward (default)

Simple. Pays the prompt 8 times.

### Option C: KV-cache reuse of the prompt

Not evaluated in the files read.

## Trade-off analysis

Packing trades code complexity for tokens. The measured result was 31.04 s to 23.65 s per step (-24 %) on the same card and seed, with step 1 identical (same answers, same gradient norm).
The cost model keeps the flag safe: it never packs when packing loses.

## Consequences

- Easier: 24 % faster steps at the full run's sizes.
- Harder: a model without the `RopeCache` hooks (a separate reference model) takes the per-row path. Flash attention is refused.
- Harder: with `top_entropy_quantile < 1` TRL's own shuffle stays, so no packing.
- The flag is off by default. `rl_shifted.yaml` turns it on.
- Revisit: after packing, fixed cost per micro-batch (launches and host syncs) dominates the update. `docs/research/notes/cuda_streams_graphs_review_2026-11-08.md` lists the avoidable syncs.

## Evidence

| Check | Result | Source |
|---|---|---|
| Time per step | 31.04 s base, 23.65 s shared (-24 %) | `docs/research/notes/shared_prefix_2026-11-08.md` |
| Token ratio | 0.361 measured, 0.37 predicted | same |
| Equivalence | max differences 5e-7 (log-probs), 1.1e-6 (gradients) in a fresh-context review | same |

## Action items

1. [x] `packed_chunk` and `packing_wins`.
2. [x] Flag on `train_rl` and `g_probe`.
3. [ ] Remove avoidable host syncs in the update (review file, section 1).
