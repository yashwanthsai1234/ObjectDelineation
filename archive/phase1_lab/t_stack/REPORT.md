# t_stack: TRL GRPO/SFT and small VLMs on the Mac

Saved by the coordinator from the agent's final message; the agent could not write this file.
Scripts are `t1`–`t5`, with `results/*.json`, `logs/` and `versions.txt`.

Versions: TRL 1.16.0, transformers 5.20.1, torch 2.14.1, peft 0.21.2, MPS.
A concurrent lab was swapping, so wall times are upper bounds.

## Results

| # | Test | Result |
|---|---|---|
| 1 | Tokens per polygon | Qwen, Qwen3.5 and SmolVLM split every digit: minified 0–1000 integer JSON costs ~8.2 tokens per vertex. LFM2.5 has single tokens for 3-digit numbers: ~5.2. LFM's native 0–1 floats cost about 2×. A spaced JSON dump costs about 22 % more. Image tokens for a 448² crop: 196 (LFM, Qwen); 1088 for SmolVLM by default, 64 with `longest_edge=512`. |
| 2 | Zero-shot disc (5 crops; fixed centre box scores IoU 0.08) | No model draws usable polygons. Box IoU: Qwen3-VL-2B 0.76 (5/5), LFM 0.43, Qwen3.5-0.8B 0.29, SmolVLM 0. The nested holes schema gets 0 parses zero-shot. **SFT is required for every model.** |
| 3 | GRPO (no vLLM) | Runs on SmolVLM-256M (full fine-tune, 17.6 s/step), LFM-450M (LoRA), Qwen3.5-0.8B (LoRA) and Qwen3-VL-2B (LoRA, G=2, 65 s/step). Each reward is logged separately. Details below. |
| 4 | SFT overfit, 4 examples, assistant-only loss | SmolVLM-256M full fine-tune: loss 1.14 → 0.0003. LFM LoRA: 1.53 → 0.0001. 4/4 outputs valid, IoU 0.977 (same as the targets). The checkpoint loads into GRPO, and 4 GRPO steps run with the real rewards. |
| 5 | Issue #6294 | No symptoms without vLLM. Untested whether single-tile SmolVLM avoids the vLLM bug. |

GRPO details (test 3):
- `normalize_then_sum` gives no NaN with a constant term. But it blows up near-constant terms: a 0.001 legality difference flips the best sample's advantage from +1.32 to +0.15. The ε is a hard-coded 1e-4, there is no std floor, and `scale_rewards` is ignored in this mode.
- A controller `TrainerCallback` that changes `trainer.beta` or `trainer.reward_weights` takes effect at the next generation.
- β = 0 at start with full fine-tuning crashes when β is raised later (no reference model). With LoRA it works.

## Shortlist (smallest first)

1. **SmolVLM-256M:** cheapest, no coordinate prior, use HF generation (no vLLM) until PR #6718 merges.
2. **LFM2.5-VL-450M:** about 37 % fewer tokens per vertex, weak box prior, use LoRA; check the LFM1.0 licence.
3. **Qwen3-VL-2B:** strongest prior.
4. **Qwen3.5-0.8B:** stretch only.

SmolVLM-500M adds nothing over the 256M version.

## Settings (TRL 1.16)

**SFT**
- `assistant_only_loss=True`, `max_length=None`, no packing.
- Full fine-tune lr ~5e-5, or LoRA r16 on all linear layers at lr ~5e-4.

**GRPO**
- `loss_type="dr_grpo"`, `epsilon=0.2`, `epsilon_high=0.28`, `beta=1e-3` (non-zero at init).
- `mask_truncated_completions=True`, constant learning-rate schedule.
- Left padding; SmolVLM `longest_edge=512`.

**Arms**
- GDPO with σ floor and gates: pass the objdel terms as weight-0 reward functions (logging only), plus one aggregator reward function that does the per-group floor and gating. Use `sum_then_normalize` with `scale_rewards="none"`. Alternatively, subclass the advantage block.

**Controller**
- A `TrainerCallback` that mutates `trainer.beta` and `trainer.reward_weights`. Every reward function can see `trainer_state`.

## Blockers

- **Mac memory:** full fine-tune at 450M needs 12 GB and swaps.
- **Speed:** 13–65 s per step, enough for smoke tests only.
- **TRL gaps:** no GDPO std floor; β cannot go from 0 to >0 with full fine-tuning; #6294 still open.
