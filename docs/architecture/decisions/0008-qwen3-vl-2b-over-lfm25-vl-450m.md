# ADR-0008: Use Qwen3-VL-2B with LoRA. Do not train LFM2.5-VL-450M again.

**Status:** Accepted
**Date:** 2026-11-02 (commits `f8aa4ef` and `734b1ec`; `experiments/phase2/PLAN.md` revision 9)
**Deciders:** Repository owner

## Context

Phase 2 ran a model ladder and a rule: stop at the first rung that passes the entry gates (`experiments/phase2/PLAN.md` D1).
Rung 1 was LFM2.5-VL-450M with full fine-tuning. It was trained, selected and judged on 2026-11-02.
It failed the entry gates (`experiments/phase2/RESULTS.md` section 1):

| Gate | Result |
|---|---|
| G-geometry | Fail. Median IoU_tc 0.495 (shifted), 0.092 (rescaled). Needs 0.60 |
| G-precision | Fail on the zero-shot probe. Probe IoU 0.61 |
| G-topology-a | Fail. Simple controls correct 93.7 % (needs 95 %). Hole recall 0.5 % |
| G-format, G-prior-b, G-signal, G-vision | Pass |

The correction note in the same file says the failure was partly data-side.
Train rows with a nodule under 10 % of the frame were 0.5 %. The rescaled exam view had 27 %.
The training data were rebuilt (`lidc-bundle-v2`).

## Decision

- Do not train LFM2.5-VL-450M again.
- Skip Qwen3.5-0.8B (owner decision).
- The next rung is Qwen3-VL-2B-Instruct (Apache-2.0, 2.13 B parameters: language model 1.72 B, vision 0.41 B), trained with LoRA rank 16 on the language model's linear layers.
  The vision tower and projector stay frozen (`configs.lora`, `VISION` regex).
- SmolVLM-256M stays a control only.

## Options considered

### Option A: Qwen3-VL-2B with LoRA (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | Low. Same grammar, same reward |
| Cost | LoRA needs about 4.3 GB of weights and 17.4 M trainable parameters. Fits a 24 GB card with micro-batch 4 and checkpointing |
| Licence | Apache-2.0 |

### Option B: Full fine-tune of Qwen3-VL-2B

AdamW needs 16 bytes per parameter: about 34 GB before activations. It fits neither a 24 GB RTX 4090 nor a 32 GB RTX 5090 (`configs.py`, PLAN D6).
A 48 GB card would cost about twice as much per hour.

### Option C: Rerun LFM2.5-VL-450M on the rebuilt data

Declined by the owner. The correction note records "LFM is not trained again".

### Option D: Larger or other models

Qwen3-VL-4B or InternVL3.5-1B-HF only after the 2B fails, with owner approval (PLAN D1).

## Trade-off analysis

The honest reading is this. Rung 1's failure had two causes that the data cannot separate: model capacity and a data gap.
The owner chose to move up one rung and fix the data at the same time.
The decision is a bet on the larger model, and it is reversible: the grammar and the reward do not depend on the model.
The risk named in PLAN D6: a frozen vision tower limits localisation. If rung 2 fails G-geometry, the first change to consider is to train the tower.

## Consequences

- Easier: GRPO uses `lora_alpha = 32`, learning rate 1e-5 (`grpo_config`) or 5e-5 (G-probe rate). The LoRA path merges into a plain model folder at the end (`train_rl`).
- Easier: the model folder is small (4.0 GB merged).
- Harder: LoRA only. No full fine-tune on the 24 GB and 32 GB cards that the project rents.
- Harder: Qwen3-VL's 152k-token vocabulary makes logits large. SFT uses micro-batch 4 with 8 accumulation steps (`r1_qwen3vl2b.yaml`).
- Outcome so far: SFT checkpoint e2 scored IoU_tc 0.511 on the 128-item shifted-validation probe set. After 770 RL steps the same set scores 0.659 (`docs/research/notes/full_rl_2026-11-08.md`).
  The final model scores IoU_tc 0.645 on the test split (`final_eval_2026-11-08.md`).
- Not measured: whether rung 2 passes all entry gates. The R1 entry-gate report for Qwen3-VL-2B is not in `RESULTS.md`, which freezes only rung 1.

## Action items

1. [x] Config `r1_qwen3vl2b.yaml`, LoRA in `configs.lora`.
2. [x] Rebuild training data with the scale ladder (`lidc-bundle-v2`, then v3).
3. [ ] Record the rung 2 entry-gate verdict in `RESULTS.md` as a new frozen section.
