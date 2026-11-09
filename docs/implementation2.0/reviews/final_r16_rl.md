# Final review round 16: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read in the sources: B.md B2 (B2.0-B2.4), B3 (B3.1-B3.10), B.6.2 (Gates B1-B3); C.md C5, C6, C7, C8, C11,
C.6, C.7; 00_front sections 1 (order, "End state under (b')"), 3 (model) and 6.3 (V3 table, amount rule, launch
lines); 99_back decisions 15-19; E.md E4 (for the cross-phase eval path). Earlier RL reviews `final_r13_rl*`,
`final_r14_rl*`, `final_r15_rl.md`, and every "Review log (final round 15)" block of review-log.md. Round 15 changed
nothing in B2/B3/C (its applied items are data-2, arch-1 and cuda-1), so this round looks for items earlier rounds
did not raise.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL GRPO, OPSM, RLOO or vLLM LoRA serving. No item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:351-362` (`is_ratio`), `:500-616` (C10 skip, `_compute_loss`, `set_is_ratio`, `log`), `:670-717` (`patch_vllm_prompts`, `lora_merged_params`); `experiments/phase2/train_rl.py:1-75` (`--resume` restarts at step 0 when no checkpoint exists, `final/` after `trainer.train`); `experiments/phase2/g_probe.py:75-100` (`logging_steps` 1, `--grpo` JSON applied last); `src/objdel/train/configs.py:170-196`; `experiments/phase2/vast/policy.py:80, 484-490`. Gave items 2 and 4. |
| engineering:debug | Traced a resume after a crash between step 20 and the first checkpoint (step 50) through B2.4 rule 4's `opsm_delta.json` rule (item 4); traced a post-E `final_eval --adapter` through E4's `VllmSampler` (item 2). |
| engineering:system-design | Checked the cross-phase order (front "End state under (b')": `after-raise` merges after E0 by a path map) against E4's behaviour change of `final_eval` (item 2). |
| engineering:architecture | Re-checked front section 3 (keep Qwen3-VL-2B) against the HF Hub today and the EGM paper; no change. |
| engineering:documentation | Each item names the source file and heading, quotes the doc line and the code line or measurement, and gives an exact fix. |
| huggingface-skills:huggingface-best | HF API `models?pipeline_tag=image-text-to-text&sort=trendingScore` and `models?author=Qwen&search=Qwen3.8` (2026-11-09): the new Qwen3.8 family starts at 27B (no model up to 4B); `LiquidAI/d1-2B` (1,918,181,888 params, `base_model: LiquidAI/LFM3-VL-2B`, `license:other`) is a fine-tune of the parked LFM3-VL-2B (front section 3 row 5). No change to the model decision. |
| huggingface-skills:trl-training | Read TRL 1.17.1 in `.venv`: `grpo_trainer.py:882` (threshold copied at init), `:1031-1037` (no ref model at beta 0 or with PEFT), `:2776-2824` (IS modes), `:2880-2907` (`sum_then_normalize`, `scale_rewards="none"`: advantages = rewards - group mean, so RLOO's zero-mean A passes unchanged, B3.4 holds), `:3066-3089` (OPSM), `:3184-3290` (loss; KL added after the IS weight; dr_grpo normaliser `size(0) x max_completion_length` at :3286); `grpo_config.py:340-360` (`use_bias_correction_kl` default True); `generation/vllm_generation.py:43, 356-379, 495-528, 600-603, 675-697, 700-775`. All B2/B3 code citations checked hold. |
| huggingface-skills:huggingface-papers | arXiv 2601.13633 (EGM, fetched from arxiv.org/html; not indexed on hf.co/papers): Table 4 Avg 84.2 (2B Instruct), 87.3 (4B Instruct), 83.6 / 87.2 (Thinking), 89.6 / 91.0 after SFT+RL; Table 9 the same. Front section 3 reason 2 holds exactly. |
| Web (WebSearch/WebFetch; Exa) | SkyRL `OffPolicyCorrectionConfig.token_tis_ratio_clip_high: float = 2.0` [web: github.com/NovaSky-AI/SkyRL skyrl/train/config/config.py; docs.skyrl.ai off_policy_correction], so B2.2's `_TIS` "C = 2.0, SkyRL's token default [web]" holds. |

## Verified, no item

- Full-run measurements quoted in B2.2 and B2.4 [measured now: runs/rl/train_rl/rl.json, 770 logged steps]:
  step-1 logp gap 0.009556 (doc 0.00956); R5 base 0.009563, shared 0.009515 (doc 0.00956 / 0.00952); median gap
  steps 1-20 0.0110, 21-100 0.0278, 101-300 0.0236, 240-280 0.0176, 301-500 0.0185 (all as the doc); IS ratio min
  = 0 in 98.05 % of steps (doc 98.1 %); grad norm 0.0018 / median 0.00475 / max 0.0221 (doc 0.0018-0.022, median
  0.0048); run `lr` 5e-5 (B2.0). All 16 checkpoints of B2.0's input list exist (50-770).
- C5 adapter count: r 16 x (4,096 + 3,072 + 3,072 + 4,096 + 3 x 8,192) = 622,592 per layer, x 28 = 17,432,576; the
  196 touched tensors 50.33 M x 28 = 1.409 B; tied embedding 151,936 x 2,048 = 0.311 B. Hold.
- B2.1 `cf/opsm_share` sign and rule equal TRL's `get_off_policy_mask` (drop iff A < 0 and mean over valid tokens of
  `sampling - logp` > delta). Under C10 the trainer logps that `set_is_ratio` sees are the same detached logps TRL's
  OPSM uses, so B1's counterfactual is the quantity the later OPSM run applies.
- B3.4 Adam argument: total gradient (g_pg + beta g_kl) / r with eps / r gives the steps of (g, eps); clip 1.0 is
  never reached (max grad norm 0.022). Holds; `use_bias_correction_kl` True does not change it (it scales the KL
  term by `coef_1`, value 1 under C10, in both arms).
- B2.2 fp16 plumbing: `from vllm import LLM` at module level (`vllm_generation.py:43`) and `self.llm = LLM(...)`
  resolve the global at call time, so the `functools.partial` patch works.
- C6: TRL discards vLLM's prompt ids in colocate (`_, completion_ids, logprobs, _, sampling_mask =
  self.vllm_generation.generate(...)`, `grpo_trainer.py:1844`), and `max_num_seqs = 16 x 1 x 16 = 256`
  (`grpo_trainer.py:1186-1188`). Hold.
- B3 unmasked truncation and trainer memory: today TRL zeroes the attention mask of truncated rows, so `window`
  trims them away; under B3 they are 1,536 tokens wide. Chunks of that width already ran (naturally ended answers to
  1,493 tokens in R5, steps with `max_length` >= 1,536 in 23.6 % of full-run steps), so no new OOM case; the added
  trainer time sits inside the 27.5 s planning value. No item.

## Items

### 1. MINOR (B): "the mismatch grows as the adapter drifts" is not what the full run shows; the gap peaks at steps 21-100 and then falls

- Doc lines: B.md "B2. Number match": "The mismatch is near zero for a fresh adapter and grows as the adapter drifts
  (B.1 item 5)." and "Merge rounding ... is a systematic gap that grows with the adapter". B.md "B2.4" rule 4: "a
  fresh adapter has the smallest gap in steps 1-20 (median logp diff 0.011 against 0.028 in steps 21-100)".
- Measurement [runs/rl/train_rl/rl.json, `sampling/sampling_logp_difference/mean`, medians]: steps 1-20 0.0110,
  21-100 0.0278, 101-300 0.0236, 301-500 0.0185, 501-770 0.0158. The gap rises fast and then falls by 43 % while the
  adapter keeps drifting. This fits merge rounding: once |d| is well above ulp(W)/2, the rounding error stays
  bounded and its share of |d| falls, which is also why B2.0 reads the step-100 relative error.
- Why minor: no rule changes. Rule 0 reads step 100, rules 2-3 read steps 50-260 (median 0.0267, near the peak),
  and the rule-4 guard reason (steps 1-20 lowest) still holds.
- Fix: B.md B2 intro: "The mismatch is near zero for a fresh adapter (0.011, steps 1-20), peaks at 0.028 (steps
  21-100) and falls to 0.016 by steps 501-770 [measured: runs/rl/train_rl/rl.json]. Merge rounding fits this: its
  error is bounded by ulp(W)/2, so its share of d falls as d grows (B2.0)." Replace "grows with the adapter" with
  "is largest early in the drift".

### 2. MINOR (B, E, 00_front): `final_eval --adapter` is specified on the pre-E trainer path, and E4 removes that path without carrying the adapter

- Doc lines: B.md "B2.4 Decide", eval rule: "With `--adapter`, `final_eval` builds the engine with
  `fast.vllm_kwargs(enable_lora=True, max_lora_rank=16)` (B2.2 F step 1) and passes `LoRARequest("eval", 1, <dir>)`
  to every generate, through `fast.vllm_generate_checked`". E.md "E4" step 3: "`final_eval` calls `VllmSampler` +
  `evaluate` ... It no longer calls `build_trainer` or `trainer._generate`", with `VllmSampler.__init__(self, ckpt,
  model, gpu_memory_utilization, max_model_len)` (step 1) and no adapter argument; grep of E.md for `enable_lora`,
  `--adapter`, `LoRARequest`: 0 hits. 00_front "End state under (b')": "If E0 has started, the owner rebases the
  branch onto post-E paths, with E's move tables as the path map".
- Why it matters: an F-trained run (the first run after Gate B3 when B2.0's error is >= 0.2) can only happen after a
  raise, which may come after E. The rebase is described as a path map, but here the target function changed:
  `VllmSampler` builds text prompts (`build_prompt`), not TRL's processor-expanded ids, so `vllm_generate_checked`
  has no role there, and the adapter must become a sampler argument.
- Why minor: B's I48 (`final_eval --adapter` with a fake engine) fails on the Mac at $0 after the rebase and names
  the missing path (caught-before-spend rule).
- Fix: E.md E4 step 1: `VllmSampler.__init__(self, ckpt, model, gpu_memory_utilization, max_model_len, adapter:
  Path | None = None)`; with `adapter`, the engine gets `enable_lora=True, max_lora_rank=16` and every generate
  passes `LoRARequest("eval", 1, str(adapter))`; `final_eval --adapter` passes it through; `tests/policy/
  test_sampler.py` takes over I48's assertions on a fake `LLM`. B.md B2.4 eval rule: add "After E4 the same rule
  lives in `VllmSampler(adapter=...)` (E4 step 1)". 00_front "End state under (b')": add "B2.4's `final_eval
  --adapter` maps to E4's `VllmSampler(adapter=...)`, not to a path".

### 3. MINOR (B): rule 3's fallback veto reads the ESS at OPSM delta 0.1, a delta the run never uses

- Doc lines: B.md "B2.4" rule 3: "fall back to `bf16-seqmask` only if the median `cf/ess_E@0.1` < 0.5"; rule 4:
  "Until step 20, delta = `start_delta`: 0.5 ... else the rule-3 q99.9 value", then the step-20 q99.5. B2.1 logs
  `cf/ess_E@{0.05,0.1,0.5}`.
- Code: OPSM keeps a sequence when `avg_seq_kl <= off_policy_threshold` [code: .venv trl/trainer/grpo_trainer.py:
  3086-3089], so a smaller delta drops more sequences and lowers the ESS. The veto at 0.1 measures a stricter mask
  than the run's start delta (0.5 or larger) and can drop `bf16-tis-opsm` for a mask it will not apply.
- Fix: rule 3: "fall back ... only if the median `cf/ess_E@0.5` < 0.5 (the start delta's ESS; when rule 3's first
  veto raised the start delta above 0.5, `cf/ess_E@0.5` is a lower bound)". Keep `cf/ess_E@0.1` as a report line.

### 4. MINOR (B): rule 4's resume rule reads `opsm_delta.json` even when the resume restarts before step 20, and the guard state is not saved

- Doc line: B.md "B2.4" rule 4: "a `TrainerCallback.on_step_end` at step 20 sets `trainer.off_policy_mask_threshold
  = delta` ... and writes `delta` to `<out>/opsm_delta.json` beside the step-20 checkpoint. A resume reads that
  file". Guard: "recompute delta once ... and never change it again".
- Code: B runs save every 50 steps (`--save-steps 50`, B1.5/B3.9), so no checkpoint exists at step 20; `train_rl
  --resume` with no checkpoint prints "starting at step 0" and trains from step 0 [code:
  experiments/phase2/train_rl.py:42-45], while `opsm_delta.json` from the failed attempt stays in `<out>`. Read as
  written, the restarted run applies the frozen delta from step 1, not `start_delta`. The guard's ring buffer and its
  "fired" flag live only in memory, so a resume after the guard fired can fire it a second time. Also, HF passes
  `model`, `optimizer` etc. to callbacks but not the trainer [code: .venv transformers/trainer_callback.py:543-556],
  so the callback needs a trainer handle at construction (U68 forces this).
- Fix: rule 4: "`opsm_delta.json` holds `{"delta", "set_at_step": 20, "guard_fired", "guard_step"}` and is rewritten
  when the guard fires. A resume reads it only when the resume checkpoint's step >= 20; a resume from an earlier
  step (or from no checkpoint) deletes it and uses `start_delta`. The callback is built with the trainer
  (`OpsmDelta(trainer)`)." U68: add a resume from step 0 with a stale file (delta = `start_delta`) and a resume after
  the guard fired (no second change).

### 5. MINOR (C, 00_front): "G-EQ1, G-EQ4 and G-EQ5 re-run in the first run that uses it" names no command for a `train_rl` run

- Doc lines: C.md "C.6" protocol: "If B2.4 later picks `fp16-seqmask` or `lora-seqmask` for the runs after Gate B3,
  G-EQ1, G-EQ4 and G-EQ5 re-run in the first run that uses it, because the dtype or the push changes." 00_front
  "6.3" V3 item 2b repeats it. C.6 table: G-EQ1 test `c_speed_check --pixels`, G-EQ4 `c_speed_check --encoder`,
  G-EQ5 `c_speed_check --sync` (under F `train_rl --f-check`).
- Problem: G-EQ1 and G-EQ4 are separate `c_speed_check` processes, not part of `train_rl`, and `c_speed_check` has no
  `--numerics` flag in C (grep of C.md: 0 hits for a preset option on it). Under `fp16-seqmask` G-EQ5 is the merged
  check, which only `c_speed_check --sync` runs. So "in the first run" has no command and no dtype input.
- Why minor: that run has no plan in this doc yet; its plan is written after Gate B3.
- Fix: C.6: "re-run as one `c_speed_check --pixels --encoder --sync --numerics <preset>` step before `train_rl` in the
  plan of the first run that uses the preset (under `lora-seqmask` the G-EQ5 part is the F check inside `train_rl
  --f-check`, so only `--pixels --encoder` run there); `c_speed_check` gains `--numerics`, which loads the model and
  the vLLM engine in that preset's dtype." 00_front 6.3 item 2b: point to that C.6 sentence.

### 6. MINOR (A, B, C): the round-15 RL minors are still neither applied nor logged

- review-log.md "Review log (final round 15)" lists only arch-1, data-2 and cuda-1. Checked in the sources today:
  C.md C11 step 3 still says "unverified for 0.32.0" (2 hits in C.md); B.md B.8.1 still lists "Items 2a, 2b and 5"
  (line 2140); B3.10's credit sentence has no beta-0 clause (grep "KL removed": 0 hits); B.md B2.2 F step 2.4 and C.md
  C5 F sync step 4 still say "the prefix-cache key holds the LoRA name, not its weights" (`final_r12_rl.md` item 1).
- Fix: apply `final_r15_rl.md` items 1-4 as written there, or add one row per item to the round-16 logs of B and C
  that rejects it with a reason.

## Counts

blocking 0, major 0, minor 6 (items 1-6).
