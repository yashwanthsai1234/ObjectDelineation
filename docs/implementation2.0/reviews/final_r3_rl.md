# Final review, round 3: RL for LLM/VLM (B2, B3, C, model choice)

Reviewer role: RL research scientist (GRPO, IS correction, reward design, TRL/vLLM internals).
Target: `REPO/docs/implementation2.0/Implementation.md` (14,506 lines, round 3). Earlier files read:
`reviews/final_r1_rl.md`, `reviews/final_r2_rl.md` (18 items), and the logs B.R3, C.R4 and the front "final round 2"
log. All 18 round-2 RL items are handled in the text. None is repeated here. Where a round-2 fix created a new defect,
the item says so.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers GRPO numerics, LoRA serving or IS correction. No item comes from it. Nothing installed. |
| engineering:code-review | I read TRL 1.17.1 in `.venv`: sequence IS ratio = exp of the SUMMED token gap (`grpo_trainer.py:2783-2790`); OPSM uses `sampling - trainer` (`:3077-3089`); dr_grpo normaliser (`:3286`); truncation mask (`:2603-2611`); `self.off_policy_mask_threshold` read per loss call (`:882, :3190`); sleep(level=2) after every generate and the re-push on wake (`vllm_generation.py:376-379, 500-528, 602-603, 771-773`). I read `fast.py:350-363, 495-717` and `configs.py:60-197` (`vllm_enable_sleep_mode: True`, `exclude_modules=VISION`, `adam` defaults). Items 1, 2, 4, 5, 6 come from this. |
| engineering:debug | I traced option F end to end through vLLM 0.32.0 source (the box's pinned version). Two failure paths found: level-2 sleep discards the base weights that F no longer pushes (item 1), and `load_inplace` on the generate request reloads the adapter from disk at every engine step (item 2). |
| engineering:system-design | I traced which variables change between B1 and B3 when rule 0 picks F (item 8) and how the F arm's timing gate compares two rentals (item 7). |
| engineering:architecture | ADR-0013 couples the serving fix (F) with the IS mode of the control (item 3). ADR-0014's "AdamW is scale-invariant" premise uses v4 gradient sizes (item 6). ADR-0022's title contradicts ADR-0013 F (item 11). |
| engineering:documentation | Each item has a severity, a section, the problem and an exact fix. |
| huggingface-skills:huggingface-best | HF API, trending `image-text-to-text`, 2026-11-08: the new small entries are above 4.0 B (`microsoft/Mage-VL-Cowork` 4,741,793,792; `maritaca-ai/S5-nano` 4,659,865,088 parameters [measured: HF API safetensors.total]). The section 3 candidate list stands. Item 9 is about the 4B SFT recipe, not the list. |
| huggingface-skills:trl-training | TRL facts behind items 3-5 and 14: sequence-mode ratios sum over all completion tokens; `mask_truncated_completions=False` puts truncated rows into the IS sums; OPSM is applied to `per_token_loss` before the dr_grpo constant normaliser. |
| huggingface-skills:huggingface-papers | I re-read the OPSM definition (DeepSeek-V3.2, as implemented in TRL) and Dr. GRPO (2503.20783) for the normaliser argument in item 14. |
| Web | vLLM v0.32.0 source on GitHub (raw files, 2026-11-08): `vllm/lora/request.py` (`load_inplace` field: "forces reloading the adapter even if one with the same lora_int_id already exists"); `vllm/lora/worker_manager.py:274-323` (`_apply_adapters` calls `add_adapter` for every active request; `load_inplace` -> `_load_adapter` from disk); `vllm/v1/worker/gpu_model_runner.py:2234-2241` (`set_active_loras` in every `_prepare_inputs`); `vllm/lora/lora_weights.py:37-43, 284-291` (`optimize()` folds the scaling into `lora_b`); `vllm/config/lora.py:45, 138-139` (`lora_dtype` "auto" = the model dtype); `vllm/v1/engine/llm_engine.py:469-471` (`add_lora`); `vllm/entrypoints/llm.py:836-858` (level 2 discards weights and KV). vLLM LoRA docs "In-Place LoRA Reloading" (server endpoint example). |

## Items

### Blocking

**1. [blocking] B (B2.2 "Option F plumbing" step 2, U70, G05) and C (C5 "Under option F", G-EQ5). Under F the plan
stops pushing the base weights, but level-2 sleep discards them every step.**
Problem: `grpo_config` sets `vllm_enable_sleep_mode: True` [code: src/objdel/train/configs.py:193]. TRL then calls
`llm.sleep(level=2)` after every generate [code: .venv trl/generation/vllm_generation.py:771-773]. Level 2 "discards
all GPU memory (weights + KV cache)" [web: vllm v0.32.0 entrypoints/llm.py:850-853]. TRL's only way back is the full
re-push in `sync_weights` on the next generate [code: vllm_generation.py:500-528, 602-603]. B2.2 F step 2 says "The
merged push (`lora_merged_params`) is skipped under F". U70 asserts the sync "never" calls `lora_merged_params`. C5
says "The push is the LoRA adapter only". So in mode B (today's mode, and C5's "always valid" fallback), vLLM decodes
from freed memory from step 2 on. The LoRA slot buffers are model memory too, so they are discarded as well, while
vLLM's adapter manager still lists id 1 as loaded. The 15-step F arm then fails (NaN or reward collapse). Rule 0
falls to A for a plumbing bug, not for a numerics reason.
Exact fix:
- B2.2 F step 2, replace with: "Each sync under F, in this order: (1) if sleeping, `empty_cache()` and
  `wake_up(tags=["weights"])`; (2) mode B only: push the UNMERGED base weights with the C5 one-call push, values
  `p.data` (no delta; `lora_merged_params(model, with_delta=False)`); mode A: no base push after the first sync;
  (3) write the adapter and register it once (item 2); (4) `llm.reset_prefix_cache()`. The prefix-cache key holds
  the LoRA name, not its weights, so blocks computed with the old adapter would be reused."
- C5: replace "The push is the LoRA adapter only" with "mode A: adapter only; mode B: base W without delta plus the
  adapter". C5 step 1 then reads "values `p.data + delta` under merge presets, `p.data` under F (mode B)".
- U70: "under F in mode B, `load_weights` receives every base name with values `torch.equal` to `p.data` (no delta);
  in mode A it receives nothing after the first sync; `reset_prefix_cache` runs once per sync".
- G05 and G-EQ5 under F: run the equality check after one sleep/wake cycle, not only after the first sync.

**2. [blocking] B (B2.2 F step 2, ADR-0013 F text, U70). `load_inplace=True` on the generate request reloads the
adapter from disk at every engine step.**
Problem: the plan passes `lora_request=LoRARequest("policy", 1, <tmpfs dir>, load_inplace=True)` to every
`generate` call. In vLLM 0.32.0 the model runner calls `set_active_loras` in every `_prepare_inputs`, that is, once
per scheduler step [web: v1/worker/gpu_model_runner.py:2234-2241]. That calls `_apply_adapters`, which calls
`add_adapter` for each active request. With `load_inplace` True it runs `_load_adapter` (a safetensors read of the
whole adapter), `remove_adapter` and `add_adapter` [web: lora/worker_manager.py:274-323]. One generate call has
about as many engine steps as the longest answer (mean 1,271 tokens [measured: runs/r5/shared_rl.json]). So F
re-reads a 70 MB adapter (17,432,576 fp32 parameters [measured: B2.2]) about 1,300 times per RL step. The F arm then
fails the 1.15 x step-time rule, and F is rejected for a call-site bug.
Exact fix: in B2.2 F step 2 and the ADR-0013 F text write: "After each sync, call
`llm.llm_engine.add_lora(LoRARequest("policy", 1, <dir>, load_inplace=True))` once [web: v1/engine/llm_engine.py:
469-471]. The `patch_vllm_prompts` wrapper passes `lora_request=LoRARequest("policy", 1, <dir>)` with the default
`load_inplace=False`." Add `llm_engine.add_lora` to the A7 private-name start-up check. U70: "`add_lora` is called
exactly once per sync with `load_inplace=True`; every `generate` call carries a `LoRARequest` with
`load_inplace=False`". Add a box assertion to the F arm: the vLLM LoRA load count per RL step is 1.

**3. [blocking] B (B2.2 F step 4, B2.4 rule 0, Gate B2, U70, G05) and C (G-EQ5). "The vLLM LoRA tensors equal the
trainer's adapter tensors" can never pass, so rule 0 can never pick F.**
Problem: three vLLM facts make bit equality false by construction.
- `lora_dtype` "auto" casts the adapter to the model dtype, bf16 [web: vllm config/lora.py:45, 138-139]. The trainer
  keeps the adapter in fp32 (B2.2 fp16 step 4; `fast.py` `no_lora_input_cast` docstring).
- `optimize()` multiplies `lora_b` by the scaling and sets the scaling to 1 [web: lora/lora_weights.py:37-43]. Here
  the scaling is alpha / r = 32 / 16 = 2 [code: configs.py `LoraConfig(r=16, lora_alpha=32)`].
- vLLM stacks q, k, v into `qkv_proj` and gate, up into `gate_up_proj` (`PackedLoRALayerWeights`, same file
  :284-291). So the names and shapes differ from PEFT's.
Gate B2 and rule 0 require "LoRA tensors equal". The check fails on every run, and B3 falls to A or E.
Exact fix:
- Save the adapter already cast to bf16 (vLLM casts to bf16 anyway, so no information is lost):
  `model.save_pretrained(dir, safe_serialization=True)` after a bf16 cast of a state-dict copy, or write the
  safetensors file directly from `{k: v.to(torch.bfloat16)}`.
- Define the check: for each layer and module, the vLLM slice of the packed `lora_a` equals `A.to(bf16)`, and the
  slice of `lora_b` equals `(2 * B.to(bf16))` with `torch.equal`. Multiplying by 2 is exact in bf16 except at
  overflow. The module map is q/k/v -> `qkv_proj` slices 0/1/2 and gate/up -> `gate_up_proj` slices 0/1.
- Add an output-level check that does not depend on vLLM internals: way (f) against way (b) on the 256-answer batch
  (B2.3b) under F must show a mean |gap| at most half of the B1 run's (f) against (b). This is the number F exists
  to fix.
- Rewrite G05, U70, Gate B2 ("LoRA tensors equal under the packed map, bf16, scaling folded") and G-EQ5 the same way.

### Major

**4. [major] B (B2.1 `cf/opsm_share`, the quantile buffer, B2.4 rule 4, U54, U68). The sign of the OPSM gap is not
defined, and B2.1's own gap has the opposite sign.**
Problem: B2.1 defines "the signed sequence log gap Σ_t (logp - sampling)", which is trainer minus vLLM. TRL's OPSM
uses `kl_div = sampling_per_token_logps - per_token_logps` and keeps a sequence when its mean is `<= delta`
[code: .venv trl/trainer/grpo_trainer.py:3077-3089]. Only B.2 states the TRL sign. `cf/opsm_share@δ`, the buffer of
"per-sequence mean gaps", rule 4's "q99.5 of the per-sequence mean gap" and U54's "mean gap 0.2" give no sign. An
engineer who reuses B2.1's gap computes the counterfactual and the live delta from the wrong tail. The frozen delta is
then wrong for the whole B3 run.
Exact fix: in B2.1 add a named quantity `seq_kl = sum_t((sampling - logp) * ok) / max(1, sum_t ok)`, with NaN
sampling logprobs set to 0 as TRL does (`ok` = valid token). State that `cf/opsm_share@δ`, the quantile buffer and
rule 4 use `seq_kl` and no other gap. Keep `sampling/seq_logdiff_*` for the trainer-minus-vLLM sum. U54: add a
negative-advantage sequence with trainer-minus-vLLM mean +0.2 (so `seq_kl` = -0.2). It is kept by OPSM at 0.1.
The mirror sequence (`seq_kl` = +0.2) is dropped.

**5. [major] B (B2.4 rule 4, U68). The live OPSM delta is fixed where the mismatch is smallest, nothing watches it
after step 20, and the patch point is not named.**
Problem: B3 starts from `rl-full-3ep` with a fresh adapter (B3.9; C8 "a fresh LoRA"). B.1 item 5 shows the mismatch
is smallest in steps 1-20 (median logp diff 0.011) and 2.5 times larger in steps 21-100 (0.028) [measured:
runs/rl/train_rl/rl.json]. Rule 4 freezes delta at the q99.5 of steps 1-20. Under E (no F), a larger share of
negative-advantage sequences can cross that delta later. OPSM drops only negative-advantage sequences, so the kept
advantages no longer sum to 0 per group. The update then leans toward raising the logprob of sampled answers. No
run guard or gate reads the live masked share. Also, rule 4 does not say how the frozen value reaches TRL. TRL reads
`self.off_policy_mask_threshold` at each loss call [code: grpo_trainer.py:882, 3190]. U68 says a resume "reads it
from `rl.json`", but `rl.json` is written at the end of the run (the full run's file holds TRL's log and timing
[measured: docs/research/notes/full_rl_2026-11-08.md]).
Exact fix:
- B2.1: log `sampling/opsm_masked_share` (share of sequences OPSM drops in the active run) and
  `sampling/opsm_masked_adv_sum` (sum of the dropped advantages).
- B2.4 rule 4: "a `TrainerCallback.on_step_end` at step 20 sets `trainer.off_policy_mask_threshold = delta` and
  writes `delta` to `<out>/opsm_delta.json` beside the step-20 checkpoint. A resume reads that file." Add the attribute
  name to the A7 start-up check.
- Add a guard: if `opsm_masked_share` exceeds 2 % in 20 consecutive steps, recompute delta once as the q99.5 of the
  last 20 steps, log the change, and never change it again. Test the guard (U68b) on a synthetic rising gap.
- Price the rule on the B1 run at $0: apply the steps 1-20 q99.5 of the B1 buffer to steps 50-300 as a
  counterfactual, and report the share it would mask.

**6. [major] B (ADR-0014, B3.4 "Why no scale constant", B3.7 rule 6). Removing the std division makes the gradient
about r times smaller, so AdamW's `eps` stops being negligible.**
Problem: the premise is "per-parameter gradients are far above eps = 1e-8 [estimate: 0.0048 / sqrt(trainable LoRA
parameters)]". That uses v4's median grad norm 0.0048 [measured: rl.json]. 0.0048 / sqrt(17,432,576) = 1.15e-6 root
mean square per element [estimate: arithmetic], about 115 eps. B3.7 rule 6 itself predicts that v5 advantages are r
times smaller than v4's, and the plan scales beta by 1/r for that reason. At r = 10 [estimate], the RMS is 1.15e-7,
about 11.5 eps. Many LoRA elements sit well below the RMS, and for an element at 0.1 x RMS, eps equals the gradient
and halves the step. So v5 lowers the effective learning rate unevenly across parameters. That is a second change
in the B3 A/B besides the reward.
Exact fix: Adam with gradient g / r and eps / r takes exactly the same steps as with g and eps. So in B2.4 and B3.9
set `adam_epsilon = 1e-8 / r` beside `beta = 1e-3 / r`, from the same r of B3.7 rule 6. Record both in ADR-0014.
Replace the sentence "neither is active" with this rule. Add U-test: two AdamW runs on a toy model, (g, eps) and
(g / r, eps / r), give equal parameters after 5 steps (`torch.allclose`, rtol 1e-6).

**7. [major] B (B2.4 rule 0, Gate B2, B.8.1, B.8.3 item 9) and C (C.8.2). The F arm compares step times across two
rentals, and its cost has two values.**
Problem: rule 0 accepts F if "its s/step is at most 1.15 x the B1 run's s/step over B1 steps 1-15". The F arm runs in
B-S2 on another host. B1 runs in V1. The host moves the step: the CPU phase is 3.0-3.2 s of the R5 step (BRIEF
section 2), and it scales with the host's CPU. Also, B1 steps 1-15 and the F arm's 15 steps sample different
answers, so the decode tail (4.3 ms per token of the longest answer [measured fit: B_gpu_audit.md:40]) differs too.
A 15 % rule cannot be read across two hosts and two samples. B.8.1 prices the arm at 6.9 min = $0.07 with no start-up. C.8.2 prices
the same arm at 13 min = $0.14, with 6 min start-up per new process. An F engine needs `enable_lora` at
construction, so it is a new process.
Exact fix:
- Rule 0's timing part: measure the LoRA cost inside the F process, on one host. After the first F sync, generate
  the same 32 x 8 batch twice, with `lora_request` and with `lora_request=None` (the base, same engine). Order them
  ABAB, seeded. Pass if the median LoRA / base generate-time ratio is <= 1.15. Keep the cross-run s/step as a report
  line only.
- Make the F arm the first 15 steps of the B3 run itself (v5, unmasked rows, F). If it fails, restart B3 under A or E.
  That removes one process start-up.
- One price in B.8.1, B.8.3 item 9 and C.8.2: 15 steps + 6 min start-up (12.9 min at 27.5 s = $0.14), or $0 extra
  when the arm is B3's own start.

**8. [major] B (B2.4 rule 0, B3.9, Gate B3, B3.10). With F near-certain, the B3 A/B changes four things, and B3.10
credits all of them to v5.**
Problem: rule 0 picks F when the step-100 bf16 relative error is >= 0.2. The seed run gives 0.531 [measured: B.1
item 5], so F is the expected B3 preset. Against B1, the B3 run then changes the reward (v5), the truncation mask
(B3.5), the serving path (F: vLLM samples from the whole update instead of one with 46-69 % of its delta elements
erased), and beta (1e-3 / r). B2.4 says only that ADR-0013/0014 "state that the B3-vs-B1 comparison includes the
numerics change". B3.10 still makes v5 "the only reward" on a Gate B3 pass, and archives v5 on a fail. F alone can
move the miss share, so Gate B3 cannot tell the two apart.
Exact fix: pick one rule and write it into B2.4, Gate B3 and B3.10.
- (a) Default if the V1 cap allows: F joins B1. Its Mac fakes and the A7 check pass before V1. V1 starts with a
  15-step F check (item 7's same-process timing). If F fails, B1 runs on `bf16-seqmask`. Gate B1's
  difference-in-differences rule (far gain minus near gain) already removes a global gain such as a serving fix. B3
  then changes the reward and the truncation mask only.
- (b) Otherwise B3 keeps B1's preset. F is adopted only after Gate B3, for later runs.
In both cases, B3.10 says what was credited: "v5 + unmasked truncation" (a) or "v5 + unmasked truncation on the B1
preset" (b).

**9. [major] 00_front (section 3 "If the gate fires" step 0, 99_back model decisions). The 4B SFT recipe differs
from the 2B recipe, so a 4B gain is not a model-size gain.**
Problem: step 0 says "Use `sft_config` on the B1 rows' targets (LoRA r 16, alpha 32, same schedule)". The 2B SFT
checkpoint behind `rl-full-3ep` came from `r1_qwen3vl2b.yaml`: view `lidc-bundle-v3` (all items, single-doctor items
included), WSD over 3 epochs, decay branches 0.5/1/2/3, branch e2 chosen on val, seed 1, effective batch 32
[code: experiments/phase2/configs/r1_qwen3vl2b.yaml; docs/research/notes/full_rl_2026-11-08.md "starting from the
SFT checkpoint e2"]. The B1 rows hold only the about 9,600 scheduled windows of RL-eligible items (2 or more doctors,
B.3.2). Those are the same windows that the 4B RL step then scores. "Same schedule" is undefined for a different
dataset, and the branch choice is missing. So a 4B-vs-2B difference mixes model size, SFT data and branch choice.
Exact fix: step 0 reads: "SFT the 4B with `r1_qwen3vl2b.yaml`, changing only `model: qwen3-vl-4b` and the output
dir: same view `lidc-bundle-v3`, same WSD, seed 1. Train only the e2 branch (decay from the 2-epoch checkpoint), so
no branch selection is needed and val stays unread. Price it from the 2B SFT ledger rows times the 4B/2B step-time
ratio of a 50-step timing probe." Add the SFT minutes and dollars to the gate text, and record the recipe in 99_back
model decisions.

### Minor

**10. [minor] B (B2 intro, B.7 row "Merge rounding", B2.4 rule 0). "No IS mode can undo merge rounding" is not
correct.**
Problem: the IS ratio uses vLLM's own sampling logprobs, which come from the rounded policy μ. IS is built to correct
exactly this kind of π ≠ μ gap: the untruncated sequence ratio is unbiased for any μ with known logprobs. The real
case for F is that the gap is systematic and grows with the adapter. That makes truncated or masked IS biased and
noisy (IS min = 0 in 98.1 % of steps), so removing the gap at its source is better than weighting it.
Exact fix: replace the sentence in the B2 intro, B.7 and rule 0 with "merge rounding is a systematic gap that grows
with the adapter; truncated and masked IS modes then trade bias for variance on every step; F removes the gap at
its source".

**11. [minor] 00_front (section 2.5 diagram and sources; section 4 "RL trainer" row trigger; ADR 0022 title), C
(C5 paragraph "ADR-0013 ... records it as option F", C10 ADR-0022 text), 99_back (decisions 16 and 16b). Stale
option-F text.**
Problem: B now picks F by rule 0 from the B2.0 table. Several places still say otherwise:
- 99_back 16: "Option F [only if B2.3b points to merge rounding]". 16b and front section 4: trigger "B2.3b points to
  merge rounding". B2.0 already shows it.
- C5: "RL review round 2 item 1 proposes ... <= 15 % ... B.md has not changed yet". B's rule is 1.15 x, and B has
  changed.
- ADR-0022 title "(LoRA serving not adopted)" contradicts ADR-0013 F.
- Section 2.5 draws only "push merged weights W + BA".
Exact fix: 99_back 16, "Option F [picked by B2.4 rule 0 when the B2.0 step-100 relative error is >= 0.2; its arm per
item 7]". 16b and the front trigger, "the F arm fails twice on a TRL or vLLM bump". C5, cite B2.4 rule 0 and delete
the "has not changed" sentence. ADR-0022 title, "GPU pixel backend, K3, one-call weight push (LoRA serving decided by
ADR-0013)". Section 2.5, add "under F: base W (mode B only) + `add_lora` in place".

**12. [minor] C (C8 "Files"). C8 tells engineers to put beta 0 into "the B run configs", but B1 fixes beta 1e-3.**
Problem: C8 "Files: the B run configs `--grpo '{"beta": 0.0}'`". B.8.3 item 5 and 99_back decision 15 fix B1 at
1e-3 and default B3 to 1e-3 / r. An engineer who follows C8 changes `b1_rows.yaml` and breaks Gate B1's "data only"
claim.
Exact fix: C8 "Files: only if the owner picks 0.0 in decision 15, set `--grpo '{"beta": 0.0}'` in `b3_v5.yaml`.
`b1_rows.yaml` keeps 1e-3." Add a config test: `b1_rows.yaml` has beta 1e-3.

**13. [minor] B (B3.5b LocGuard tests). The false-stop rate is tested only on independent noise.**
Problem: the unit test uses "seeded noise with flat true trends (per-step sd as in R5)". Per-step means in a real run
are autocorrelated (the policy moves slowly, and the 20-step blocks repeat the cell mix). An OLS slope t-statistic
over 50 steps is then too large, so the < 1 % false-stop rate does not hold.
Exact fix: compute the slope's standard error with Newey-West (lag 5) in LocGuard. Add a test that fits AR(1) to the
full run's per-step `rewards/iou_tc` [measured source: runs/rl/train_rl/rl.json], simulates flat series with that
coefficient, and asserts < 1 % stops in 300-step runs.

**14. [minor] B (B3.5, B.1 item 8, run reports). The token weight of an unmasked truncated row is not reported.**
Problem: dr_grpo divides by a constant (`size(0) x max_completion_length` [code: grpo_trainer.py:3286]). A
1,536-token truncated row with A = -b gets about 12 times the gradient mass of a mean 131-token answer
[measured lengths: runs/r5/shared_rl.json]. The negative signal then falls on all its tokens, including a valid
outline prefix. This is unbiased for the sequence objective, but it can dominate a step's update.
Exact fix: B3.5 logs `completions/truncated_grad_share = sum over truncated rows of |A| x tokens / sum over all rows
of |A| x tokens`. The run report prints its median and max. If the median exceeds 10 % over steps 50-300, the owner
may switch on the C11 option "mask the first min(len, T_i) tokens" for truncated rows too.

**15. [minor] B (B2.0 test). The erased share has no stated denominator for zero deltas.**
Problem: the planted case has two zero entries. If zeros count, the erased share is 3/4, not the stated 0.5. The test
"d = 0 gives share 0" implies that zeros are excluded, but the definition does not say so.
Exact fix: B2.0 definition, "share of NON-ZERO delta elements with |d| < ulp(W)/2; 0 when no element is non-zero".
Keep the test values.

## Counts

blocking 3, major 6, minor 6.
