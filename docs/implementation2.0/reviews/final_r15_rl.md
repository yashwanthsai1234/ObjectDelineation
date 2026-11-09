# Final review round 15: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read in the sources: B.md B.1 (items 5, 8), B1.5, B1.9, B2.0-B2.4, B3.1-B3.10, B.6.2 (Gates B1-B3), B.7,
B.8.1, B.8.3, B.9 rows U38, U45, U49, U51, U68, U70, U73, U76, I22b and the U51 sketch; C.md C5, C6, C7, C8, C9,
C11, C.6; 00_front sections 1, 1.1, 3 and 6.3 (V3 table, V3 amount rule, V3 launch lines); 99_back decisions 13-19;
A.md A7; `reviews/final_r13_rl*.md`, `reviews/final_r14_rl*.md`; review-log.md "Review log (final round 14)" of
every phase. The two round-14 RL changes (rl-2: B3 soft train step, step 3b `--merge-last`, V3 pull job; rl-3:
`--opsm-start-delta` and `start_delta`) were checked against the code and hold (below).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL GRPO, OPSM, RLOO or vLLM LoRA serving. No item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:351-362` (`is_ratio`), `:524-616` (C10 skip, `_compute_loss`, `set_is_ratio`, `log`), `experiments/phase2/train_rl.py:1-60` (constant LR via `make_trainer`; `final/` only after `trainer.train`), `src/objdel/train/configs.py:150-196` (RL: `lr_scheduler_type` "constant", warmup 0), `experiments/phase2/vast/policy.py:70-82, 180-193, 525-536`. Confirmed that the fp16 arm (60 steps) and B1 (260 steps) share the LR at steps 41-60 (no schedule confound), so B2.3's paired comparison holds. |
| engineering:debug | Traced the round-14 rl-2 path: crash in step 3 -> `~` soft step logs and continues -> step 3b merges the last `checkpoint-N` on `--ckpt {rd}/ckpt` (`rl-full-3ep`) -> step 4 evals `{rd}/b3-v5/final`. Traced rl-3: preset 0.5 -> `--opsm-start-delta X` after the preset -> TRL copies it once at init (`grpo_trainer.py:882`) and reads it per loss call (`:3190-3203`) -> step-20 callback overwrites it. Both hold. |
| engineering:system-design | Compared the B3 A/B against what B3.9/B3.10 credit, under the default and under the owner's beta-0 alternative (item 1). Checked the V3 item 4 minute sum against the B1 re-eval line of the same row (item 2). |
| engineering:architecture | Checked front section 3 against the HF API (below) and the code (`configs.py:73-77`): no change to the model decision. |
| engineering:documentation | Each item names the source file and heading, quotes the doc line and the code line or measurement, and gives an exact fix. |
| huggingface-skills:huggingface-best | HF API `models?author=Qwen&sort=createdAt` (2026-11-09): newest small image-text model is `Qwen/Qwen-Drive-1.0-4B` (2026-08-27, 4,539,265,536 params, Apache-2.0, `base_model: Qwen/Qwen3.5-4B`, architecture `QwenDriveForPlanning`, a driving fine-tune with BEV and planning heads). It is the parked Qwen3.5 hybrid family and above 4.0 B, so front section 3 needs no change (round 14's "no new Qwen VLM up to 4B" is imprecise but has no effect). |
| huggingface-skills:trl-training | Read TRL 1.17.1 in `.venv`: `grpo_config.py:268` (`mask_truncated_completions` default False; ours sets True at `configs.py:188`), `grpo_trainer.py:2603-2611` (truncation mask), `:2776-2824` (IS modes), `:3065-3089` (OPSM; denominator `mask.sum`), `:3184-3290` (loss: OPSM mask, IS weight on the policy term, `beta * per_token_kl` added after, dr_grpo normaliser `size(0) x max_completion_length`). Gave item 1's code line. |
| huggingface-skills:huggingface-papers | No new paper claim in scope; the RLOO identity, the Adam (g / r, eps / r) argument and the masked-truncation bias identity are algebra and were re-derived (they hold). |
| Web (WebSearch/WebFetch loaded; raw GitHub at tag v0.32.0) | `vllm/entrypoints/llm.py:418-442` at v0.32.0: `sampling_params: SamplingParams | Sequence[SamplingParams]`, "When it is a list, the list must have the same length as the prompts and it is paired one by one with the prompt" (item 4). `vllm/lora/worker_manager.py:289-325` at v0.32.0: with `load_inplace=True` the LRU manager reloads, `remove_adapter`s, re-adds and `activate_adapter`s, so C5's "Re-adding after each wake also restores the LoRA slots that level-2 sleep discards" holds. |

## Verified, no item

- Round-14 rl-2 fix: B.md B3.9 now reads "Step 3, the `train_rl` step, writes `--out {rd}/b3-v5` with `--save-steps 50`
  ... It is a soft step (first word `~`) ... Step 3b is the hard step `train_rl --merge-last --model qwen3-vl-2b --ckpt
  {rd}/ckpt --out {rd}/b3-v5`"; `{rd}/ckpt` is `rl-full-3ep` (front 6.3 V3 item 0), the B3 start model, so the merge
  base is right. Front 6.3 V3 launch lines load and unload the pull job. I22b asserts the `~`, `--save-steps 50` and
  step 3b. Holds against `policy.py:484-490` and `train_rl.py:55-60`.
- Round-14 rl-3 fix: B2.4 rule 4 "Until step 20, delta = `start_delta`"; B2.2 `--opsm-start-delta` is the only
  path; U45, U51 (sketch lines 3001-3005), U68 test the number. If `cf/opsm_share@0.5` > 2 %, more than 2 % of
  sequences have `seq_kl` > 0.5, so q98 > 0.5 and the rule-3 q99.9 is always above 0.5 (the rule always raises).
- B3.4 RLOO: R_i - (S - R_i)/(G - 1) = G/(G - 1) (R_i - S/G). Worked check 0.1 x (0.62 - 0.374/7) = +0.057. Holds.
- B3.1 popcount column sum: masks 0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 0xFF00FF00, 0xFFFF0000 select bit 0-4 of the
  position, weights 1, 2, 4, 8, 16. Holds.
- B3.4 Adam scaling: weight decay 0, grad norm 0.0018-0.022 under clip 1.0, KL weight beta / r: AdamW with
  (g / r, eps / r) takes the steps of (g, eps). Holds. RMS 0.0048 / sqrt(17,432,576) = 1.15e-6 = 115 eps. Holds.
- B2.2 F comparator: PEFT `set_adapter("b1")` flips `requires_grad` to the active adapter and back on
  `set_adapter("default")`, and the restore deletes `b1` before step 1, so the trainable set at step 1 is the fresh
  `default` adapter. The F sync order under level-2 sleep is consistent with vLLM 0.32.0's LRU LoRA manager (above).
- RL LR is constant with no warmup (`configs.py:189`), so the `--steps 60` fp16 arm and the 260-step B1 run train
  at the same LR over the compared steps 41-60 (B2.3).

## Items

### 1. MINOR (B, C, 99_back): under the owner's beta-0 alternative for B3, the fixed ADR-0014 credit sentence and B3.9's "keeps v4's balance" are false

- Doc lines:
  - 99_back "RL settings" decision 15: "B3: `1e-3 / r` (B3.7 rule 6)]. Alternative: beta 0.0 for B3 only."
  - C.md "C8" Files: "only if the owner picks 0.0 in decision 15, set `--grpo '{"beta": 0.0}'` in
    `experiments/phase2/configs/b3_v5.yaml`."
  - B.md "B3.9 B3 run": "Against B1, the reward changes, plus the truncation mask (part of v5, ADR-0014), plus `beta`
    and `adam_epsilon` scaled by the same r (which keeps v4's balance and Adam's step)."
  - B.md "B3.10": "What is credited: "v5 + unmasked truncation on the B1 preset `bf16-seqmask`". ADR-0014 writes this
    sentence." B.6.2 Gate B3: "so the gate credits "v5 + unmasked truncation"".
- Code: TRL adds the KL term to the loss only when beta is non-zero: `if self.beta != 0.0: per_token_loss =
  per_token_loss + self.beta * per_token_kl` [code: .venv trl/trainer/grpo_trainer.py, `_compute_loss`, after the
  IS weight line `:3268-3269`]. With beta 0 the B3 arm drops the KL pull that B1 had (beta 1e-3), so the A/B changes
  reward, truncation mask and the KL term.
- Why minor: it is an owner-chosen alternative, the Gate B3 pass rules are outcome rules that do not read the credit
  sentence, and the default (`1e-3 / r`) is consistent.
- Fix: B.md B3.10 "What is credited": "v5 + unmasked truncation on the B1 preset `bf16-seqmask` (plus 'KL removed'
  when the owner chose beta 0.0 in decision 15)". B3.9: "... scaled by the same r (which keeps v4's balance and
  Adam's step); under the beta-0 alternative only `adam_epsilon` is scaled and the KL term is a third change".
  B.6.2 Gate B3 credit bullet: same parenthesis.

### 2. MINOR (00_front): V3 item 4 prices the B3 val eval without the 2.6 min of untimed imports that the same row adds to the B1 re-eval

- Doc line: 00_front "6.3" V3 table, item 4 basis: "3 + 5 [estimate: as V1 smoke] + 13.0 (205 s + 2 x 1,084 x
  0.193 s + 2.6 imports) + 6 + 260 x step time + 10.4". The B3 val eval is the same one-process, two-view
  `final_eval` form (B3.9: "one `final_eval` step in the B1 form (B1.5 step 4)"), so it has the same untimed imports
  (00_front section 6 inputs: "Process imports before the load are not timed").
- Measurement: 2 x 1,084 x 0.193 s + 205 s = 623 s = 10.4 min; with the imports 13.0 min. Item 4 and the B-S2 total
  are low by 2.6 min ($0.03 at r_s). `BudgetStop` and the V3 amount are unaffected (the amount is set at the end of
  step 260, and the $0.35 tail = 30.9 min >= 13.0 + 2 + 3).
- Fix: item 4: "then B3 on `lidc-val-shifted` and `lidc-val-edge` in one process (13.0 min, imports included)";
  minutes 159.2 / 142.5; recompute the B-S2 and "Items 0, 1, 3 and 4" totals and the 99_back totals that cite them.

### 3. MINOR (B): B.8.1 lists the optional V3 rows as "2a, 2b and 5", while front 6.3's one definition lists 1, 2a, 2b, 3 and 5

- B.md "B.8.1": "B-S2 = rental V3 items 0 and 4 ... Items 2a, 2b and 5 (the test eval, when the owner has closed
  the 4B gate) are optional V3 rows, approved by name; front 6.3 holds the one definition". 00_front "6.3": "One
  definition: B-S2 (B.8.1) is V3 items 0 and 4, the only required V3 items. Items 1, 2a, 2b, 3 and 5 are optional
  rows."
- Fix: B.8.1: "The other V3 rows (front 6.3 "One definition") are optional and approved by name." (no list copy).

### 4. MINOR (C): C11 step 3 marks the per-row `SamplingParams` list "unverified for 0.32.0"; the pinned source confirms it

- C.md "C11" step 3: "`LLM.generate` with a `SamplingParams` list of the same length as the prompts gives each
  answer its own `max_tokens` [web: docs.vllm.ai LLM.generate docstring, latest; unverified for 0.32.0]."
- Source: vLLM tag v0.32.0 `vllm/entrypoints/llm.py:421` `sampling_params: SamplingParams |
  Sequence[SamplingParams] | None = None` and `:438-442` "When it is a list, the list must have the same length as
  the prompts and it is paired one by one with the prompt."
- Fix: replace the citation with "[web: vllm v0.32.0 entrypoints/llm.py:421, 438-442]". G-RC0 stays the box check.

### 5. MINOR (A, B, C): earlier RL minors are still neither applied nor logged

- review-log.md "Review log (final round 14)" for Phase B lists only rl-2 and rl-3; no row applies or rejects the
  other items. Still in the sources:
  - `final_r12_rl.md` items 1-3: B.md "B2.2" F step 2.4 and C.md "C5" F sync step 4 still say "the prefix-cache key
    holds the LoRA name, not its weights" (B.md line 1342, C.md line 824); `_llm_weights_sleeping` and
    `off_policy_mask_threshold` are absent from A.md "A7" (grep: 0 hits in A.md) while B2.4 rule 4 says "The
    attribute name joins the A7 start-up check"; `ess` has no unit.
  - `final_r13_rl.md` item 1 (verdict minor): truncated-row IS weight under `sequence_mask` not logged.
  - `final_r14_rl.md` items 4-6: B2.1 `seq_kl` still divides by `max(1, sum_t ok)` while TRL divides by
    `mask.sum(...).clamp(min=1.0)` [code: grpo_trainer.py:3086-3087]; B2.2 F step 4 still "seed s = 0" with no
    per-request seed list (grep "seed = 1000": 0 hits); no `seq_gaps.f32` truncation on `--resume` (grep
    "truncates `seq_gaps": 0 hits). The verdict-downgraded `final_r14_rl.md` item 1 (`stopped_by`) is also absent.
- Fix: apply those items as written in their files, or add one row per item to the Phase A/B/C round-15 logs that
  rejects it with a reason.

## Counts

blocking 0, major 0, minor 5 (items 1-5).
