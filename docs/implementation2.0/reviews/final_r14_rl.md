# Final review round 14: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read in the sources: B.md B.1 (items 5, 8), B1.5, B2.0-B2.4, B3.1-B3.10, B.6.2 (Gates B2, B3), B.9 U38, U51
(and its sketch), U62, U68, U73, I16, I22b; C.md C5, C6, C7, C8, C9, C10, C11, C.6, C.9.1 (U-RC1..RC6); 00_front
sections 1.1, 3 and 6.2/6.3 (V1 and V3 launch lines, V3 table, V3 amount rule); 99_back decision 15;
`reviews/final_r13_rl.md` and its verdict; review-log.md (no round-13 rows exist for B or C, so B and C are
unchanged since round 13).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL GRPO, OPSM or vLLM LoRA serving. No item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:351-362` (`is_ratio`), `:524-602` (C10 skip, `_compute_loss`, `set_is_ratio`), `src/objdel/train/configs.py:108-195`, `experiments/phase2/train_rl.py:19-68`, `experiments/phase2/g_probe.py:75-100`, `experiments/phase2/vast/policy.py:484-494`. Gave items 1-2 (plan-step chaining, `final/` saved after any stop). |
| engineering:debug | Traced a LocGuard stop and a `train_rl` crash through the B3 plan: stop -> `should_save` -> `trainer.train` returns -> `final/` saved -> step 4 eval -> Gate B3 with no stop rule (item 1); crash -> `&&` chain ends -> no B3 eval (item 2). Traced the OPSM delta from the preset literal through TRL's `__init__` (line 882) to the loss (line 3190) (item 3). |
| engineering:system-design | Compared the B1 and B3 run plans side by side: B1 has a soft train step, `--merge-last`, a 20-min pull job and a "stopped early" decision row; B3 has none of the four (items 1-2). |
| engineering:architecture | Checked front section 3 (model decision) against the code (`configs.py:70-77`) and the HF API: no new Qwen VLM up to 4B since round 12 (Qwen3-VL Embedding/Reranker 2B are retrieval models; Qwen3.5 small models are still the parked hybrid-attention family). No model-choice item. |
| engineering:documentation | Each item names the source file and heading, quotes the doc line and the code line, and gives an exact fix. |
| huggingface-skills:huggingface-best | Queried `huggingface.co/api/models?author=Qwen&search=VL` and `search=Qwen3.5`, sorted by date (2026-11-09): newest VL generators are the 2025-10-31 GGUF re-uploads; nothing changes the table in front section 3. |
| huggingface-skills:trl-training | Read TRL 1.17.1 in `.venv`: `grpo_config.py:791-971` (defaults: `sequence_mask`, cap 3.0, clip_min None, OPSM None, `mask_truncated_completions` False), `grpo_trainer.py:882` (OPSM threshold read from args at init), `:1365-1377` (RepeatSampler), `:2603-2611` (truncation mask), `:2776-2824` (IS modes), `:3066-3089` (OPSM: denominator `mask.sum`), `:3184-3290` (loss; IS weight on the policy term; dr_grpo normaliser), `generation/vllm_generation.py:356-378, 690-760` (engine-level seed, no per-request seed). Gave items 3 and 4-5. |
| huggingface-skills:huggingface-papers | Re-checked arXiv 2601.13633 (EGM) exists with the title the front cites; the table numbers were verified in earlier rounds and are not re-raised. |
| Web (WebSearch/WebFetch loaded) | arXiv abstract fetch for 2601.13633; HF API calls above. |

## Verified, no item

- R5 step-1 reference values of G-EQ6 (IoU_tc 0.5467, legal 0.8984, mean length 117.89, grad norm 0.003264, logp
  difference 0.009515), C7 "6 of 15 R5 steps had a naturally ended answer > 1,200 (max 1,493)" and C11 "5 of R5
  steps 3-15 hit 1,536" [measured: runs/r5/shared_rl.json, 2026-11-09].
- Full run `lr` 5e-5 and 770 steps [measured: runs/rl/train_rl/rl.json]; B2.0 and B3.9 use 5e-5.
- B.1 item 8 bias identity: E[1{nt}(R - b) grad log pi] = grad J - b grad P(nt) when R = 0 on truncation. Correct.
- B3.4 RLOO identity, the worked check (+0.057), U73's 0.626 and U-RC1's [214, 1,264, 1,536]: arithmetic holds.
- B3.4 Adam argument: with `weight_decay` 0 (TrainingArguments default) and `max_grad_norm` 1.0 never reached
  (grad norm 0.0018-0.022), Adam with (g / r, eps / r) and KL weight beta / r takes exactly the steps of (g, eps,
  beta). Holds.
- B3.7 rule 9 / C7 step 5: TRL's soft overlong penalty formula matches `trl/rewards/other_rewards.py:83-95`.
- TRL has no native vLLM LoRA path (no `enable_lora` or `LoRARequest` in `trl/`), so B2.2's `vllm_kwargs` wrapper
  is needed; TRL builds `LLM(...)` without `dtype` (`vllm_generation.py:356-372`). Holds.

## Items

### 1. MAJOR (B, 00_front): a B3 run stopped early by LocGuard (or cut by BudgetStop) has no Gate B3 rule, so a run that LocGuard flagged as gaming the location term can still pass Gate B3 and make v5 the only reward

- Doc lines:
  - B.md "B3.5b Location guard (run guard)": "It watches whether the model learns to move outlines toward the nodule
    without hitting it. ... Rule: stop (and save) only when, in **two disjoint consecutive windows** of 50 logged
    steps, both hold: (a) the least-squares slope of the location mean has t-statistic > 3; (b) the upper 95 % bound
    of the IoU_tc slope is < 0."
  - B.md "B3.9 B3 run": "A lower value cuts the 260 steps, and Gate B3 then compares unequal runs." and the step
    order "0 CUDA build, 1 G08b smoke, 2 B1 re-eval, 3 `train_rl`, 4 B3 val eval, 5 G-RC0 (soft)".
  - B.md "B.6.2 Run gates", Gate B3: the after-run rules are miss share, both-found IoU, edge IoU and legal share
    only; no line reads a LocGuard stop or the step count.
  - B.md "B3.10 Clean after gate B3": "Pass: v5 becomes the only reward. Move `RewardV4` ... to
    `archive/reward_v4/`".
  - 00_front "1.1 After V1: decision table" has a row "B1 stopped early: `rl.json` has `partial: true` ... a pass
    does not lead to B-S2 until the owner reviews it". No B3 row exists, and the Pass row's final candidate is "Best
    of B1 and B3 on val-select by the Gate B3 rules".
- Code: a guard stop saves and ends training normally (B.9 I16: "`BudgetStop` and `LocGuard` set
  `should_training_stop` and `should_save`"); `train_rl` then merges and saves `final/` after `trainer.train`
  returns [code: experiments/phase2/train_rl.py:55-60], so step 4 evaluates the stopped model and Gate B3 is
  computed on it like a full run. I22b asserts this exact step list for `b3_v5.yaml`, so no $0 check flags it.
- Why it changes the outcome: LocGuard's stop is the plan's own evidence that the location term is being gamed
  (B.7 risk row "Location term gamed"). A run stopped at step 100-150 has trained less than B1 (260) yet can still
  meet the four Gate B3 rules on val-select (for example through the miss-share rule, which a policy that moves
  outlines toward the nodule can lower without raising IoU). B3.10 would then archive v4 and make the gamed reward
  the only reward, and the 1.1 Pass row could pick B3 as the one model that reads test. A BudgetStop cut (slow host:
  the V3 amount rule assumes 27.5 s/step) gives the "unequal runs" that B3.9 names but nobody rules on.
- Fix:
  - B.md B.6.2 Gate B3, add a first after-run bullet: "**Stop rule.** If LocGuard stopped the run (`rl.json`
    `stopped_by: "LocGuard"`), Gate B3 fails whatever the val numbers; B3.10 'Fail' applies and ADR-0014 records the
    LocGuard windows. If BudgetStop or a crash ended the run before step 260 (`partial: true` or `steps` < 260),
    Gate B3 is computed and labelled 'partial run, N steps'; a pass makes v5 the only reward only after the owner's
    written review."
  - B.md B3.5b: "On a stop, LocGuard writes `stopped_by: "LocGuard"`, the step and both window slopes into
    `rl.json`." BudgetStop writes `stopped_by: "BudgetStop"`. `train_rl` copies `stopped_by` into `rl.json`.
  - 00_front 1.1: add a note under the table: "Gate B3 has the same 'stopped early' rule (B.6.2 Gate B3 'Stop
    rule'); a LocGuard stop is a Gate B3 fail."
  - Tests: extend I16 (`rl.json` gets `stopped_by` for each guard) and U49 (planted `stopped_by: "LocGuard"` with
    passing val rules -> Gate B3 fail; planted `partial: true` -> "partial run, N steps", no automatic pass).

### 2. MAJOR (B, 00_front): the B3 run lacks the crash protection that B1 has (soft train step, `--merge-last`, 20-min pull job), so a crash or lost host late in the 260 steps gives no Gate B3 result after about $1.8 of V3

- Doc lines:
  - B.md "B1.5" step 3b: "Why: `train_rl` writes `final/` only after `trainer.train` returns ..., and the `&&` join
    skips every later step after a failure. Without step 3b a crash at step 210 leaves `trainer/checkpoint-200` but
    no B1 eval."
  - 00_front 6.2 "Why the pull job (final round 11, fresh-arch-2)": "Without it, nothing moves checkpoints off the box
    during the 2-3 h run, and a lost host loses every checkpoint."
  - B.md "B3.9 B3 run": "So the order of `b3_v5.yaml` is: 0 CUDA build, 1 G08b smoke, 2 B1 re-eval, 3 `train_rl`,
    4 B3 val eval, 5 G-RC0 (soft)." Step 3 has no `~` and there is no `--merge-last` step.
  - 00_front 6.3 "V3 launch lines": `run.sh <run> experiments/phase2/configs/b3_v5.yaml # item 4; wait for EXIT`;
    no `launchctl load` of the pull job (the V1 lines in 6.2 have it).
- Code: a non-soft step joins with `&&` [code: experiments/phase2/vast/policy.py:484-490, `jobs = " && ".join(...)`],
  so a `train_rl` failure skips steps 4-5; `final/` is written only after `trainer.train` returns [code:
  experiments/phase2/train_rl.py:55-60].
- Why it matters: V3 item 4 is 156.6 min; the B3 `train_rl` starts at about minute 66 and ends at 185.2 (front 6.3
  "V3 amount rule"). The plan rates a crash or lost host as Medium likelihood for the same kind of run (00_front FR16),
  and fixed it for B1 in round 11. The same failure in B3 spends the session and leaves Gate B3 open; without the
  pull job a lost host also loses the checkpoints, so not even a later $0 merge is possible. No planned $0 check
  covers it: I22b asserts the current step list ("then the `train_rl` step ... one `final_eval` step after
  `train_rl`").
- Fix:
  - B.md B3.9: make step 3 soft (`[~, -m, experiments.phase2.train_rl, ...]`) and add step 3b, the hard step
    `train_rl --merge-last --model qwen3-vl-2b --ckpt {rd}/ckpt --out {rd}/b3-v5`, exactly as B1.5 step 3b. New
    order: "0 CUDA build, 1 G08b smoke, 2 B1 re-eval, 3 `train_rl` (soft), 3b `--merge-last`, 4 B3 val eval,
    5 G-RC0 (soft)". The `partial: true` record then feeds item 1's "partial run" rule.
  - 00_front 6.3 V3 launch lines: load the pull job (the 6.2 `sed ... && launchctl load ...` line) before `run.sh
    <run> experiments/phase2/configs/b3_v5.yaml`, and `launchctl unload` after `teardown.sh`, with a pointer to the
    6.2 "Why the pull job" note (the 30-min idle destroy applies the same way).
  - B.9 I22b: assert the soft first word of the `train_rl` step and the `--merge-last` step 3b; I47-style argv parse
    of `b3_v5.yaml`'s `--merge-last` step. FR16: "during the V1 or V3 run".

### 3. MAJOR (B): B2.4 rule 3 raises the OPSM starting delta, but rule 4, the preset literal, U38 and U68 fix delta at 0.5 for steps 1-20, and nothing carries the raised value into the run

- Doc lines:
  - B.md "B2.4 Decide" rule 3: "if `cf/opsm_share@0.5` > 2 % of sequences, raise the starting delta to the q99.9 of
    `seq_kl` over the per-sequence values of B1 steps 50-260 (read from the B1 run's `seq_gaps.f32`, B2.1)".
  - B.md "B2.4 Decide" rule 4: "Until step 20, delta = 0.5 (TRL docstring)."
  - B.md "B2.2 Numerics presets": `"bf16-tis-opsm": {**_TIS, "off_policy_mask_threshold": 0.5},  # ADR option E;
    delta replaced live (B2.4)`.
  - B.md B.9 U38: "bf16-tis-opsm = bf16-tis + `off_policy_mask_threshold` 0.5 (the starting value; `None`
    refused)"; U68: "delta = 0.5 for steps 1-20; at step 20 it is set once to q99.5 ...".
  - B.md B.9 U51 and its sketch: `assert d.preset == "bf16-tis-opsm" and d.raise_delta_to == "q99.9"`. No other
    line reads `raise_delta_to`; `train_rl --numerics` takes only a preset name (B2.2: "`train_rl --numerics
    {bf16-seqmask,...}`").
- Code: TRL copies the threshold from the config once in `__init__` (`self.off_policy_mask_threshold =
  args.off_policy_mask_threshold` [code: .venv trl/trainer/grpo_trainer.py:882]) and applies it at every loss call
  (`:3190-3203`), so the preset's 0.5 is the delta of steps 1-20 unless something overrides it.
- Why it changes the outcome: rule 3 fires exactly when B1's counterfactual says delta 0.5 would drop more than 2 %
  of sequences (all negative-advantage, so the kept advantages lean positive; B2.1 "OPSM drops only
  negative-advantage sequences, so the kept advantages no longer sum to 0 per group"). With the tests as written the
  engineer builds 0.5 for steps 1-20 and rule 3's output is computed, tested and silently unused; an engineer who
  follows rule 3 instead breaks U38/U68. Both readings pass every planned $0 test (U51 tests one function, U68 the
  other), so the conflict reaches the paid run after Gate B3. It also decides the 20 steps whose q99.5 then freezes
  the delta for the rest of the run (rule 4), so the starting value feeds the frozen value.
- Fix:
  - B.md B2.4 rule 4: "Until step 20, delta = `start_delta`: 0.5 (TRL docstring), or the rule-3 value when rule 3
    fired. `start_preset` / `restart_preset` return `start_delta` as a number (the q99.9 computed from the B1
    `seq_gaps.f32`), not the string 'q99.9'."
  - B.md B2.2: `train_rl --opsm-start-delta X` (only valid with `bf16-tis-opsm`; argparse error otherwise) sets
    `off_policy_mask_threshold` after the preset; the run plan writes the ADR-0013 value; `rl.json` records it.
  - B.9: U51 asserts `d.start_delta == pytest.approx(np.quantile(planted_steps_50_260, 0.999))` when `opsm05 =
    0.021` and `== 0.5` when 0.020; U68 runs with a planted `start_delta` 0.7 and asserts 0.7 for steps 1-20; U38
    keeps 0.5 as the preset default only.

### 4. MINOR (B): B2.1's counterfactual `seq_kl` divides by the non-NaN token count, TRL's OPSM divides by the whole mask, while B2.1 calls it "TRL's formula"

- B.md "B2.1": "`seq_kl = sum_t((sampling - logp) * ok) / max(1, sum_t ok)`, ... with NaN sampling logprobs set to
  0 as TRL does." TRL: `kl_div = torch.nan_to_num(kl_div, nan=0.0)`; `avg_seq_kl = seq_kl_sum /
  mask.sum(dim=1, keepdim=True).clamp(min=1.0)` with `mask = completion_mask` (or times `tool_mask`) [code: .venv
  trl/trainer/grpo_trainer.py:3076-3087, 3148]. With any NaN token the counterfactual `seq_kl` is larger than the
  active OPSM's value. Today vLLM returns no NaN logprob at T = 1 with no filter, so the effect is near zero, and
  U62 compares with TRL's function (it catches the gap only if its 200 random tensors hold a NaN).
- Fix: B2.1: denominator `max(1, sum_t mask)` (the completion mask, NaN tokens included), "exactly TRL's
  `get_off_policy_mask`"; U62: plant NaN in 10 % of `sampling` entries of half of the 200 tensors.

### 5. MINOR (B): the F comparator's "fixed seed" is not defined; a per-request `seed=0` on the 8 copies of each prompt makes them identical

- B.md "B2.2" F step 4: "each sent 8 times in the run's own request form (256 answers ...), seed s = 0 for arms (0),
  (i) and (ii)", and step 4b.4: "`SamplingParams(max_tokens=512, min_tokens=512, ignore_eos=True, seed=s)`".
  TRL's own request form has no per-request seed: it seeds the engine (`seed=accelerator.process_index // tp`)
  [code: .venv trl/generation/vllm_generation.py:365-366, 696]. A vLLM `SamplingParams.seed` seeds each request's
  own generator, so 8 requests with the same prompt and the same seed sample the same answer: the comparator then
  scores 32 distinct answers, not 256, and the timing arms decode 8 copies of one row. The pass rule still works
  (the mean over about 32 x 130 tokens), so it is minor.
- Fix: B2.2 F step 4: "Seed per request: `seed = 1000 x s + row index` (row index 0-255 of the 32 x 8 batch), the
  same list in arms (0), (i), (ii) and in both timing arms"; U76 asserts 256 distinct seeds and equal lists across
  arms.

### 6. MINOR (B): `seq_gaps.f32` gets duplicate blocks after a resume

- B.md "B2.1": "It appends them to `<out>/seq_gaps.f32` (... one 256 x 3 block per logged step, in step order ...)";
  B1.5 "Resume" restarts from the last checkpoint (every 50 steps). Steps logged after that checkpoint and before
  the crash stay in the file, so a resume appends their blocks again and every later block is off by up to 49 steps;
  rules 3 and 4 then read the wrong steps. A resume needs a new rental (00_front FR16), so it is rare.
- Fix: B2.1: "On `--resume` from `checkpoint-N`, `train_rl` truncates `seq_gaps.f32` to N x 256 x 3 x 4 bytes and
  rewrites `seq_gaps.json` before training." Extend R01b: a 4-step run stopped after step 3 and resumed from
  checkpoint-2 gives a file of exactly 4 blocks equal to an unbroken run's.

### 7. MINOR (A, B, C): the round-12 and round-13 RL minors are still neither applied nor logged

- review-log.md has no round-13 row for B or C, and no `final_r12_rl_verdict.md` exists. Still open:
  `final_r12_rl.md` items 1-3 (prefix-cache key wording in B2.2 F step 2.4 and C5 F sync step 4; `_llm_weights_sleeping`
  missing from the A7 list; `ess` unit) and `final_r13_rl.md` item 1 (verdict: minor; the truncated-row IS weight
  under `sequence_mask` is not logged).
- Fix: apply them as written in those files, or add a round-14 log row in B and C that rejects each with a reason.

## Counts

blocking 0, major 3 (items 1, 2, 3), minor 4 (items 4-7).
