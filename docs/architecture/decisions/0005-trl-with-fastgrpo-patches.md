# ADR-0005: Keep TRL and patch it (FastGRPO) instead of changing the RL framework

**Status:** Accepted
**Date:** 2026-11-06 (commit `f68b32d`, "RL-ready: fast GRPO trainer"; TRL choice 2026-10-19, `experiments/phase1/PLAN.md` D8)
**Deciders:** Repository owner

## Context

Phase 1 chose TRL's `GRPOTrainer` with `loss_type="dr_grpo"` (`experiments/phase1/PLAN.md`, D8).
The first profile of TRL 1.17 on an RTX 5090 (P1, 2026-11-06) found two problems:

- **Speed.** TRL ran at 80.3 s per step. TRL pads all 256 answers of a step to the longest one and runs every forward at that width.
- **Correctness.** vLLM sampled off-policy for Qwen3-VL. The token log-probability gap reached 25 nats and the importance ratio mean was 0.65 to 0.81.
  TRL sends vLLM processor-expanded prompt ids with the image, and vLLM's own processor expands them again (`experiments/phase2/perf/p1_results.md`).

## Decision

Keep TRL 1.17.1 (pinned in `pyproject.toml`). Subclass `GRPOTrainer` as `FastGRPO` in `src/objdel/train/fast.py`.
Patch only the places that cost time or break correctness.

| Patch | What it does |
|---|---|
| `FastForward._get_per_token_logps_and_entropies` | Cut every forward chunk to its own window, rounded to 64 columns (or to compile buckets) |
| `VisionCache` | Run the frozen vision tower once per distinct image per step |
| `RopeCache` | Compute Qwen3-VL's 3-D rope positions once per row per step |
| `ImageMemo`, `OpenOnce` | Open and prepare each image once per step (32 images, not 4,096 opens) |
| `one_pad_per_image` plus a per-prompt check | vLLM samples from exactly the image features the trainer scores |
| `skip_old_logps`, `set_is_ratio` | Drop the old-log-probability pass. The update's own detached log-probabilities set the importance ratio |
| `by_length`, `by_group_length` | Cut update micro-batches from the length-sorted batch |
| `no_lora_input_cast` | Remove the fp32 copy of each LoRA input under bf16 autocast (bit-identical) |
| `prefetch_rewards` | Score on a background thread while the trainer computes reference log-probabilities |
| `lora_merged_params` | Push W + delta to vLLM without writing the bf16 base ([ADR-0011](0011-lora-weight-sync-without-merging.md)) |
| `packed_chunk` | Shared-prefix forward ([ADR-0007](0007-shared-prefix-forward.md)) |

Every patch is a method override or a mixin. No file of TRL is copied or edited.

## Options considered

### Option A: TRL plus FastGRPO (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | Medium. Overrides depend on TRL internals |
| Cost | No migration. Reward functions, config and callbacks stay |
| Team familiarity | High |

### Option B: TRL as shipped

80.3 s per step and off-policy sampling (above). Not acceptable.

### Option C: Another framework

`experiments/phase1/PLAN.md` D8 lists the options rejected in Phase 1: verl and EasyR1 (no SmolVLM support then), LLaMA-Factory (no GRPO), and RLOO (called the most hack-prone).
No migration was tried after P1.

## Trade-off analysis

Option A keeps the reward integration (`TRLReward` as weight-0 loggers plus one aggregator) and the `Controller` callback.
The price is coupling to TRL 1.17.1 internals. The pin and the tests carry that risk: `tests/test_train_fast.py`, `tests/test_train_trl.py`.
Each patch was measured on the card before it stayed (`experiments/phase2/profile.py`).

## Consequences

- Easier: one pinned stack. 44.3 s per step (1.81x) with on-policy sampling: the log-probability gap mean is 0.009 nats (P1).
- Easier: each patch has a test that demands the same log-probabilities, entropies or gradients as TRL (tolerances in `tests/test_train_fast.py`, `docs/research/notes/shared_prefix_2026-11-08.md`).
- Harder: a TRL upgrade needs a review of every override. Several patches read private names (`_get_per_token_logps_and_entropies`, `vllm_generation._iter_named_params`).
- Harder: `OpenOnce` needs `dataloader_num_workers = 0`. `get_train_dataloader` raises otherwise.
- Revisit: TRL issues #6688 and #7423 track the merge drift. When TRL fixes it, `lora_merged_params` can go.

## Evidence

| Setting | s per step | Source |
|---|---|---|
| TRL 1.17 as shipped | 80.3 | `experiments/phase2/perf/p1_results.md` |
| FastGRPO, 16 x 16 | 44.3 | same |
| Full run mean, 770 steps | 27.5 | `docs/research/notes/full_rl_2026-11-08.md` |
| With `--shared-prefix`, 15 steps | 23.65 (base 31.04) | `docs/research/notes/shared_prefix_2026-11-08.md` |

## Action items

1. [x] Pin TRL 1.17.1 and add the overrides.
2. [x] Test each patch against TRL's own result.
3. [ ] Remove `lora_merged_params` once TRL ships its fix.
