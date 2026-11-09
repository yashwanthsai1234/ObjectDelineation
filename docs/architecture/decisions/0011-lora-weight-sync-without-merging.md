# ADR-0011: Push LoRA weights to vLLM as W + delta in a new tensor, never by merging in place

**Status:** Accepted
**Date:** 2026-11-08 (commit `1870328`)
**Deciders:** Repository owner

## Context

After each optimizer step, the new LoRA weights must reach vLLM.
TRL 1.17.1 `VLLMGeneration` merges the adapter into the bf16 base weights in place, pushes them, then unmerges (`trl/generation/vllm_generation.py`).
In bf16, `W + d - d` is not `W`. The base weights drift at every sync.

The full RL run found it after the fact (`docs/research/notes/full_rl_2026-11-08.md`):

- A 2048 x 2048 bf16 matrix after 770 merge and unmerge cycles with a fixed rank-16 delta had 6.4 % of entries changed. The largest change was 4.9e-4.
- The trainer and vLLM both see the drifted weights, so the importance ratio does not flag it.
- TRL tracks it as issue #6688 (maintainer-confirmed) and #7423, with fix PR #7427 open at the time.

## Decision

`FastGRPO.sync_without_merge` replaces vLLM's named-parameter iterator with `lora_merged_params`.
For every LoRA layer it computes W + delta in W's dtype, as PEFT's merge does, into a new tensor, and yields it under vLLM's parameter name.
No parameter is written. It applies on one process. FSDP and ZeRO-3 keep TRL's gathered path.
It refuses DoRA and LoRA bias (`NotImplementedError`).

## Options considered

### Option A: Add the delta in a new tensor (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | Low. One generator function |
| Cost | A temporary tensor per layer during the sync |
| Result | Same pushed values as PEFT's merge. The base stays bit-identical |

### Option B: Keep TRL's merge and unmerge

Zero code. The frozen base drifts at every step.

### Option C: Wait for TRL's fix (PR #7427)

Costs nothing now. The fix was open when the decision was made. A pinned stack could not wait.

## Trade-off analysis

Option A adds a patch that Option C would remove later. Its value is reproducibility: the saved model is the base plus the trained adapter and nothing else.
The tests demand it: the pushed weights equal TRL's merge exactly, and 20 syncs leave every parameter bit-identical while TRL's path changes them (`tests/test_train_fast.py`).

## Consequences

- Easier: the base model stays equal to the SFT checkpoint plus the adapter.
- Easier: the learning is unchanged. In R4, step 1 matched the full run's and steps 3 to 15 matched to three decimals (`r4_speed_2026-11-08.md`).
- Harder: the patch depends on `vg._iter_named_params`, a private name in TRL.
- Known limit: the held-out numbers of the full run (IoU_tc 0.659) were measured on the model as trained with the old path, drift included. They are real measurements of that saved model.

## Action items

1. [x] `lora_merged_params` and `sync_without_merge`.
2. [x] Tests for equality with TRL's merge and for no drift.
3. [ ] Remove the patch when TRL ships its fix.
