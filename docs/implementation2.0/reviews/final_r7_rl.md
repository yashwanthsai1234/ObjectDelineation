# Final review round 7: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read: B.md B2.0-B2.4 (incl. the new F step 4b), B3.1-B3.10; C.md C5-C11, C.6, C.7, C.8; 00_front section 3;
99_back "RL settings" and "Model"; review-log.md rows of rounds 4-6 for B, C, front and back (so nothing below
repeats a fixed item; round-6 rl-3..rl-7 minors are not re-filed).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL/vLLM LoRA serving or GRPO; no item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:488-620, 696-717` (`sync_without_merge` patches `_iter_named_params`, so TRL's own `sync_weights` still runs today; `lora_merged_params` merges the active adapters). Led to item 5. |
| engineering:debug | Traced the engine state through B2.2 F step 4b in mode B, call by call, against TRL's sleep calls (`vllm_generation.py:376-379, 500-503, 602-603, 771-773`). Found the hash after `sleep(level=2)` (item 1). |
| engineering:system-design | Checked what each comparator number measures. Found that the output-level ratio has a d = 0 floor (item 2) and that the timing ratio depends on answer length (item 3). |
| engineering:architecture | Checked that the rule texts agree: B2.2 step 4b.5, C5 step 5(c) and C.6 G-EQ5 all put the base hash "after the restore" in both modes (item 1). |
| engineering:documentation | Each item quotes the source file, the heading and the line or measurement that proves it. |
| huggingface-skills:huggingface-best | HF `config.json` read for Qwen3-VL-2B/4B-Instruct: vision depth 24, width 1024, DeepStack [5, 11, 17] in both; merger out 2,048 vs 2,560; text layers 28 vs 36. The front section 3 table is confirmed. |
| huggingface-skills:trl-training | Read TRL 1.17.1 `grpo_trainer.py:2777-2824` (IS modes, ratio from `old - sampling`), `3066-3089` (OPSM, `sampling - per_token_logps.detach()`), `3188-3200` (OPSM input), `3262-3286` (order OPSM, IS, `+ beta*kl`, dr_grpo constant), `1838-1841` (sync once per global step), and `vllm_generation.py` (start-up `sleep(level=2)`, `_llm_weights_sleeping`). |
| huggingface-skills:huggingface-papers | The HF paper page for 2601.13633 returned 404. I read the arXiv HTML instead. Table: 2B Instruct avg 84.2, 4B Instruct 87.3; Thinking 83.6 vs 87.2; + EGM 89.6 vs 91.0. Front section 3 reason 2 is confirmed. |
| Web | vLLM sleep-mode docs (docs.vllm.ai/en/latest/features/sleep_mode): "Level 2 sleep will discard the model parameters and the KV cache"; "the parameters are not backed up"; after a level-2 wake the weights must be reloaded. arXiv HTML 2601.13633. |

## Verified, no item

- EGM numbers in front section 3 reason 2 [paper: arXiv 2601.13633 HTML table]. Qwen3-VL 2B/4B config facts [measured:
  HF config.json, 2026-11-08].
- C7 and B.1 item 8 bias sign: masking the truncated rows removes E[1{trunc}(0 - b) grad log pi] = -b grad P(trunc), so
  the masked update is grad J + b grad P(trunc). This matches "-b grad P(not truncated)".
- B3.2 monotonicity: dR/dIoU >= 1 - alpha_loc / miss_iou >= 0 when alpha_loc <= miss_iou (loc_raw <= 1, g_area <= 1).
- B3.4 worked check: 0.1 x (0.62 - (0.29 + 0.084) / 7) = 0.0567.
- B2.4 KL/eps scaling: KL enters `per_token_loss` after the IS weight and before the dr_grpo constant
  [code: grpo_trainer.py:3271-3286], so (beta / r, eps / r) keeps Adam's step.
- C8: `if self.beta != 0.0` guards the reference pass [code: grpo_trainer.py:2827].
- F base push into LoRA-wrapped vLLM modules: round 5 already verified it in the 0.32.0 source. I did not re-check it.

## Items

### 1. BLOCKING (B, C): in mode B the comparator hashes the vLLM base tensors after `sleep(level=2)` has discarded them

- Sources: B.md "B2.2 Numerics presets", Option F plumbing step 4b.5; C.md "C5. Weight push" step 5(c); C.md
  "C.6. Acceptance gates" row G-EQ5.
- Doc lines:
  - B2.2 step 4b.5: "Restore, in both modes: `push_weights(vg, model, values="base")`, `remove_lora(2)`,
    `model.set_adapter("default")`, `model.delete_adapter("b1")`, `reset_prefix_cache()`, then in mode B
    `sleep(level=2)`. Then hash every vLLM base tensor on the host against `p.data`, one tensor at a time (the C5 step
    5(b) method). A mismatch fails the F check before step 1."
  - C5 step 5(c): "In both modes, after the comparator's base restore, every vLLM base tensor hashes equal to
    `p.data`, as in (b)."
  - G-EQ5: "in both modes, after the F comparator's base restore before step 1 (B2.2 F step 4b), every base tensor
    hashes equal to `p.data`".
- Proof:
  - TRL: "Sleep level 2 discards the weights; track it so that generate() knows it must re-push them"
    [code: .venv trl/generation/vllm_generation.py:376-379].
  - vLLM docs: "Level 2 sleep will discard the model parameters and the KV cache ... the parameters are not backed
    up." After a level-2 wake, the weights must be reloaded [web: docs.vllm.ai features/sleep_mode].
  - So in mode B the hash reads unmapped or undefined device memory. The likely result is a CUDA illegal-address
    error, which kills the process. The other result is a hash mismatch. Either way the F check fails before step 1
    on every mode-B run.
- Why mode B is a real path: B2.2 F step 2.2 says "Mode B (sleep on, today)". Under (b'), C5 mode A never lands. C.6
  also says G-MEM "can switch mode A to mode B".
- Result: rule 0 rejects F for a plumbing reason. ADR-0013 records a false result. The run restarts under the
  fallback preset and pays one start-up (about 6 min) plus the comparator.
- Fix (B2.2 F step 4b.5; mirror it in C5 step 5(c) and G-EQ5):
  - Mode A: `push_weights(values="base")`, then the host hash against `p.data`, then `remove_lora(2)`,
    `set_adapter("default")`, `delete_adapter("b1")`, `reset_prefix_cache()`.
  - Mode B: no restore push and no hash here. Step 1's own `sync_lora` re-pushes the base after the wake (F sync step
    2.2), and the C5 step 5(c) mode-B base hash runs after that sync. Then `remove_lora(2)`, `set_adapter("default")`,
    `delete_adapter("b1")`, `reset_prefix_cache()`, `sleep(level=2)`.
  - U76 gains a planted case: a fake engine whose tensors raise when read while asleep. Mode B must not read them.

### 2. BLOCKING (B): the output-level comparator's 0.5 bar is below what a perfect F can reach at B1's drift, so F fails even when it removes all merge rounding

- Sources: B.md "B2.2 Numerics presets", Option F plumbing step 4; B.md "B2.4 Decide" rule 0 (restart choice reads
  the result); U76 and G05 (`f_check.json`).
- Doc line (B2.2 step 4): "(iii) The trainer scores both answer sets with way (b) of B2.3b. Pass: mean |f_lora - b|
  <= 0.5 x mean |f_merged - b|. This is the number F exists to fix."
- Proof: the bf16 engine-vs-trainer gap does not go to 0 when merge rounding goes away.
  - At step 1 the adapter delta is exactly 0, so there is no merge error. The logp difference there is still 0.00956
    on the plain path (R5 base arm, step 1) and 0.00952 on the shared-prefix path (R5 shared arm, step 1).
    The full run's step 1 gives 0.00956 [measured: runs/r5/base/rl.json, runs/r5/shared_rl.json,
    runs/rl/train_rl/rl.json, 2026-11-08].
  - At the drift of a 260-step adapter, the full run's median `sampling/sampling_logp_difference/mean` over steps
    240-280 is 0.0176 [measured: runs/rl/train_rl/rl.json]. Other ranges: 101-300 gives 0.0236, 301-500 gives
    0.0185.
  - Best case for F: f_lora reaches the d = 0 floor. Then the ratio is 0.0096 / 0.0176 = 0.54, which is above
    0.5. The B1 adapter has 260 steps (B1.5) on top of `rl-full-3ep`, the same order of drift.
  - So the rule rejects an F that works. The restart then costs about 6 min plus 15 steps of GPU time, and ADR-0013
    records "F failed".
- Fix (B2.2 step 4, B2.4 rule 0, U76):
  - Measure the floor in the same comparator at no extra generation cost. The ABAB timing (step 4b.4) already
    generates the batch with `lora_request=None` on base W. Score those answers with the trainer under
    `with model.disable_adapter():`. The result is f0 = mean |f_base - b_base| (d = 0, so no merge error).
  - Pass if (mean |f_lora - b| - f0) <= 0.5 x (mean |f_merged - b| - f0).
  - If mean |f_merged - b| - f0 < 0.002 [estimate threshold], record "comparator uninformative" in `f_check.json`.
    The output-level check then passes, and F is decided by the tensor, NaN and timing checks.
  - Write f0, both gaps and the rule into `f_check.json` and ADR-0013.
  - U76 gains two cases: (f0 0.0096, merged 0.0176, lora 0.0096) passes; (f0 0.0096, merged 0.0176, lora 0.0140)
    fails.

### 3. MAJOR (B): the in-process LoRA / base timing ratio compares two different answer sets, so answer length drives it, not LoRA kernel cost

- Sources: B.md "B2.4 Decide" rule 0 (timing rule); B.md "B2.2 Numerics presets" step 4b.4; U51 (`restart_preset`
  ratio 1.15 / 1.151).
- Doc lines:
  - B2.4 rule 0: "generate the same 32 x 8 batch twice, with the comparator's `lora_request` and with
    `lora_request=None` (the base, same engine), in the order ABAB, seeded; pass if the median LoRA / base
    generate-time ratio is <= 1.15".
  - B2.2 step 4b.4: "(the LoRA kernel cost depends on rank and shapes, not on the adapter values)".
- Proof:
  - The kernel cost does not depend on the adapter values, but the answers do. Arm A samples from W + d_B1 and arm B
    from W (`rl-full-3ep`), so they write different answers with different lengths.
  - Decode time follows the longest answer: "Generate time grows 4.3 ms per token of the longest answer" (C.md
    "C11", [measured fit: B_gpu_audit.md:40]).
  - In the full run, 25.5 % of steps 101-300 had a 1,536-token answer (B.md "B.1" item 8 table [measured: rl.json]).
    So in about a quarter of batches one arm has a runaway and the other may not.
  - Example: longest answer 1,536 against 800 tokens means 736 x 4.3 ms = 3.2 s more in one arm, on an 11.1 s
    generate phase (R5). The ratio moves by about 0.29 for that reason alone.
  - The repeats are seeded, so each arm writes the same answers twice. The median of the ABAB pairs does not average
    this out.
  - Result: the 1.15 bar can fail a cheap F, or pass a slow one. The F check then restarts or keeps the wrong preset.
- Fix (B2.4 rule 0 timing rule and B2.2 step 4b.4):
  - Time both arms with `SamplingParams(max_tokens=512, ignore_eos=True, seed=s)`, so every row decodes exactly 512
    tokens in both arms and the ratio is kernel cost only.
  - Run three ABAB pairs with seeds s = 0, 1, 2 and take the median of the 6 ratios.
  - Keep the comparator's quality generations (arms (i) and (ii)) on the normal `SamplingParams`.
  - U51 and U76: the ratio input is the fixed-length timing, and a planted length difference does not move it.

### 4. MINOR (B): step 4b does not say when (iii) scoring runs, and "FastGRPO patches off" has no switch inside a running FastGRPO process

- Source: B.md "B2.2" step 4 (iii) and step 4b.4-4b.5.
- Step 4b.4 orders (ii), the timing and (i), and step 4b.5 then deletes adapter `"b1"`. Scoring (iii) needs `"b1"`
  active in the trainer, and no step places it. Way (b) is "the plain TRL forward (FastGRPO patches off)", but the
  process has VisionCache, RopeCache and the shared-prefix patches installed, and the doc names no way to turn them off.
  Also, "Arm (ii) runs first, on the start-up base W" is wrong in mode B, where the start-up base has been discarded
  and `sync_lora` re-pushes `p.data`. The values are the same, but the wording misleads.
- Fix:
  - Add step 4b.4a: "(iii) score both answer sets, and the floor set of item 2, with `"b1"` active, before step 4b.5."
  - Score with the run's own forward (way (a)). That is the forward the IS ratio uses. R5 step 1 shows (a) and (b)
    give the same floor (0.00952 vs 0.00956).
  - Reword the arm (ii) text as "on base W (mode A: the start-up load; mode B: `sync_lora`'s base push)".

### 5. MINOR (C): the C5 `install` replaces `vg.sync_weights`, but the doc does not name TRL's `_llm_weights_sleeping` flag, which TRL's generate reads

- Source: C.md "C5. Weight push", Interface (`push_weights`, `install`) and the U-SY8 test.
- Today `sync_without_merge` patches only `_iter_named_params`, so TRL's own `sync_weights` still clears the flag
  [code: fast.py:495-504; vllm_generation.py:503]. Once `install` replaces `sync_weights`, TRL's generate calls
  `self.sync_weights()` again whenever `_llm_weights_sleeping` is True [code: vllm_generation.py:602-603]. That gives
  two pushes per step and two `add_lora` calls under F, which fails the F "load count per RL step is 1" check. The doc
  cites "TRL's order, vllm_generation.py:500-503", and those lines include the reset, so a careful engineer gets it
  right. But U-SY8 counts per sync, not per step, so nothing catches the mistake.
- Fix:
  - In the C5 interface, state that `push_weights` and `sync_lora` wake unconditionally when sleep mode is on and set
    `vg._llm_weights_sleeping = False`, as TRL does.
  - Add `_llm_weights_sleeping` to the A7 private-name check.
  - Add one U-SY8 case: one step through TRL's generate path on the fake gives exactly one sync.

### 6. MINOR (B): `completions/truncated_grad_share` leaves out the IS weight, and the IS weight hits long rows hardest

- Source: B.md "B3.5 Wire v5 into the TRL adapter", bullet "Token weight of truncated rows".
- Doc line: "`completions/truncated_grad_share` = sum over truncated rows of |A| x tokens / sum over all rows of |A| x
  tokens ... If the median exceeds 10 % over steps 50-260, the owner may switch on C11's option".
- B3 keeps `bf16-seqmask`, where "A sequence ratio is exp(sum of token gaps), so long answers pass the cap first" (B.1
  item 5). A 1,536-token truncated row is the most likely to be masked (ratio > 3) or down-weighted to near 0. Its real
  gradient mass is |A| x sum_t(w_IS x mask), not |A| x tokens. So the 10 % owner trigger reads an overstated share.
  It may also hide that the unmasked truncation signal B3 credits is largely removed by the IS mask.
- Fix: compute the share with |A| x sum_t(w_IS,t x mask_t), using the same `importance_sampling_ratio` that
  `set_is_ratio` writes. Log the unweighted share beside it. Add both to the B3 run report and to Gate B3's
  clipped-share line.

## Counts

blocking 2, major 1, minor 3.
