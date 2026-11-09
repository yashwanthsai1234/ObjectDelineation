# Final review round 8: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read: B.md B.1 items 5 and 8, B2.0-B2.4 (with the round-7 floor arm and fixed-length timing), B3.1-B3.10,
Gate B2 and Gate B3 (B.6.2); C.md C5-C8, C11, C.6; 00_front sections 2.5 and 3; 99_back "RL settings" and "Model";
review-log.md rounds 4-7 for front, B, C and back. Round-7 minors rl-4, rl-5 and rl-6 are not filed again.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL/vLLM LoRA serving or GRPO. No item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:515-620` (`skip_old_logps`, `set_is_ratio`, `log`) and `:696-717` (`lora_merged_params`). Confirmed that the IS ratio and `sampling_logp_difference/mean` come from the update forward, as B2.1 says. Read `configs.py:60-100, 177-196`: LoRA `all-linear` minus the vision tower, `vllm_gpu_memory_utilization` 0.35. |
| engineering:debug | Traced the comparator in modes A and B against TRL's sync trigger. `_last_loaded_step = -1` means step 0 always syncs [code: .venv trl/trainer/grpo_trainer.py:1203, 1838-1841], so step 1 never meets an engine that has no "policy" adapter. No item. |
| engineering:system-design | Sized the fixed-length timing batch against the vLLM KV budget. Item 1. |
| engineering:architecture | Checked that B2.2 step 4/4b, B2.4 rule 0, Gate B2, C5 step 5(c) and G-EQ5 agree after the round-7 fixes. They agree. |
| engineering:documentation | Each item quotes the source file and heading, plus the code line or measurement. |
| huggingface-skills:huggingface-best | Re-read the front section 3 ranking against the owner rule ("any model up to 4B"). I also checked the 4B memory plan. In mode B, level-2 sleep frees vLLM's weights and KV during the update. So 0.55 x 31.84 GiB for vLLM beside 8.9 GB of 4B trainer weights fits under the 30 GB gate during generation. No item. |
| huggingface-skills:trl-training | Read TRL 1.17.1 `grpo_trainer.py`: 2603-2611 (truncation mask), 2777-2824 (IS modes), 2880-2907 (`sum_then_normalize`, `scale_rewards="none"`: centring only), 2975-2992 (`sampling_logp_difference/mean` = mean \|old - sampling\| over valid tokens), 3066-3089 (OPSM), 3262-3286 (OPSM, IS, `+ beta*kl`, dr_grpo constant). Also read `vllm_generation.py:352-379` (LLM kwargs: `gpu_memory_utilization`, `max_num_batched_tokens=4096`, sleep at start-up). B3.4's "RLOO passes through TRL unchanged" holds. So do B2.4's "(beta/r, eps/r) keeps Adam's step" and B.1 item 8's bias sign. |
| huggingface-skills:huggingface-papers | No new paper claim in scope. Front section 3 EGM numbers (2601.13633) were verified in round 7. I did not re-fetch them. |
| Web | Searched for vLLM LoRA-vs-base decode overhead numbers (Punica 8-12 % on Llama-7B, second-hand; Fireworks 10-30 % prefill). None is a merged-vs-unmerged decode figure for a 2B model, so item 3 stays minor. |

## Verified, no item

- Comparator noise vs margin. Per-step SD of `sampling/sampling_logp_difference/mean` is 0.00107 over steps 240-280
  (consecutive-difference estimate 0.00097). It is 0.00010 over steps 1-5 [measured: runs/rl/train_rl/rl.json,
  2026-11-08]. The floor-subtracted gap there is about 0.0176 - 0.0096 = 0.008. The pass bar at 0.5x gives a
  0.004 margin, about 4 SD of the noise. So the round-7 rule can separate a working F from a failed one.
- B.1 item 5 table re-measured: medians 0.011 (1-20), 0.0278 (21-100); IS min = 0 in 0.9805 of steps; median step-mean
  ratio 0.823; steps 101-300 with a 1,536-token answer 0.255 [measured: rl.json].
- B.1 item 8 bias sign: E[1{nt}(R - b) grad log pi] = grad J - b grad P(nt) when R = 0 on truncated rows.
  TRL zeroes only `completion_mask` [code: grpo_trainer.py:2603-2611]. Truncated rows return 0 terms, not None
  [code: src/objdel/reward/score.py:25], so TRL's `unscorable_mask` (NaN rows get A = 0) does not catch them, and
  B3's "A = -b" holds.
- KL placement: `+ beta * per_token_kl` comes after the IS weight and before the dr_grpo constant [code:
  grpo_trainer.py:3271-3286], so the B2.4/B3.4 r-scaling argument holds. Clipping never binds (grad norm median 0.0048
  vs `max_grad_norm` 1.0).
- Prompt length per row is 325 tokens (num_tokens delta minus 256 x mean length, R5 steps 2-5) [measured:
  runs/r5/shared_rl.json]. Item 1 uses it.

## Items

### 1. MINOR (B): the fixed-length timing batch needs more KV than vLLM's whole memory budget, so it runs with preemption, not "at batch 256"

- Source: B.md "B2.2 Numerics presets", Option F plumbing step 4b.4 and step 4 (cost line); B.md "B2.4 Decide" rule 0
  timing rule.
- Doc lines:
  - B2.2 step 4b.4: "both timed arms use `SamplingParams(max_tokens=512, min_tokens=512, ignore_eos=True, seed=s)`
    ... every row then decodes exactly 512 tokens in both arms ... Limit: a full-batch 512-token decode measures the
    kernel cost at batch 256".
  - B2.2 step 4: "(14 generations of about 7.8 s each [estimate: 5.61 s + 4.3 ms x 512 tokens ...]), about 5 min".
- Proof (arithmetic from measured inputs):
  - KV per token for Qwen3-VL-2B: 2 x 28 layers x 8 KV heads x 128 x 2 B = 114,688 B [measured:
    runs/r3/ckpt/config.json text_config].
  - Tokens held at the end of one timed generate: 32 shared prompts x 325 + 256 x 512 = 141,472 tokens
    = 15.1 GiB of KV [prompt length measured: runs/r5/shared_rl.json].
  - vLLM's whole budget is 0.35 x 31.84 GiB = 11.1 GiB, before weights (4.15 GiB) and activations
    [code: src/objdel/train/configs.py:192 `vllm_gpu_memory_utilization: 0.35`; .venv
    trl/generation/vllm_generation.py:358 passes it to `LLM`; 32,607 MiB and "Model loading took 4.15 GiB" measured
    in runs/*/env_check_vllm.log]. At 0.5 the engine had 83,584 KV tokens (runs/r4/env_check_vllm.log), so at 0.35 it
    has roughly 40-47k [estimate]. That is about 3x short.
  - So vLLM V1 preempts and recomputes rows. The 512-token count still holds, but the measured regime is about
    50-60 concurrent rows plus recompute prefill, not "batch 256". The 7.8 s per generate is also too low. The
    comparator takes longer than "about 5 min" [estimate].
  - Both arms share one engine and one schedule, so the ratio is not confounded by length. That is why this item is
    minor. But the doc's stated limit is false. Recompute prefill, where LoRA costs relatively less, also dilutes
    the ratio toward 1.
- Fix (B2.2 step 4b.4, B2.4 rule 0, U76):
  - Size the timed length from the live engine, not as a literal:
    `kv = llm.llm_engine.vllm_config.cache_config.num_gpu_blocks * block_size`;
    `T = min(512, 16 * floor((0.9 * kv - 32 * prompt_len) / (256 * 16)))`; require `T >= 64`, else the timing check
    records `kv_too_small` and fails.
  - Use `max_tokens=T, min_tokens=T`. At about 45k KV tokens, T is about 112, close to the mean answer of 131.
  - Record `T`, `kv` and vLLM's preemption count delta (must be 0) in `f_check.json`.
  - Change "kernel cost at batch 256" to "kernel cost at batch 256 and T tokens, with no preemption".
  - Recompute the cost line as 14 x (5.61 s + 4.3 ms x T) [estimate].
  - U76: a fake engine with a small `num_gpu_blocks` gives the hand T. A preemption count of 1 fails the timing.

### 2. MINOR (B, C): the 1.15 timing bar sits inside the doc's own 10-30 % LoRA-decode estimate, so F can be dropped for speed alone

- Sources: B.md "B2.4 Decide" rule 0 ("pass if the median of the 6 LoRA / base generate-time ratios is <= 1.15
  [estimate threshold]"); C.md "C5" ("LoRA kernels may slow decode [estimate 10-30 %, ...B_gpu_audit §4.6]") and
  "C.6" ("vLLM's LoRA kernels may slow decode by 10-30 %").
- By the doc's own estimate, the timing check fails F in a good part of the range. A failure restarts under fp16 or
  `bf16-tis-opsm`, so a speed check throws away the numerics fix. The owner's binding rule ranks quality first
  (BRIEF section 3). No measurement exists for this model, so the item is minor.
- Fix: in B2.4 rule 0, make a timing failure alone (all other F flags true) restart under F, and write the ratio and
  the step-time cost into ADR-0013 as an owner decision line. If the owner wants a hard bar, they set it before the
  run, next to the G-C budget (C.6), so both bars come from one step-time budget.

### 3. MINOR (B): B.1 item 5 treats an IS ratio of 0 as "masked", but phase A says a 0 can also be an fp32 underflow

- Source: B.md "B.1 Goal and why" item 5: "`importance_sampling_ratio/min` is 0 in 98.1 % of 770 steps, so at least 1
  of 256 sequences was masked in those steps." A.md "Skills used" row trl-training: "in `sequence_mask` a zero ratio
  is either '> 3' or an fp32 underflow. A6 splits the two."
- `exp(sum gap)` is 0 in fp32 when the summed gap is below about -103.97 [code: .venv grpo_trainer.py:2791
  `torch.exp(logps_diff)`; src/objdel/train/fast.py:358 `torch.exp(diff)` on the FastGRPO path], which a long runaway can reach. The two parts disagree on wording only. The B2.1
  counterfactual counts are unaffected.
- Fix: B.1 item 5 reads "so at least 1 of 256 sequences was masked (ratio > 3) or underflowed to 0 (A6 counts the
  two apart)".

## Counts

blocking 0, major 0, minor 3.
