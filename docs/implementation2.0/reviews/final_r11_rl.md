# Final review round 11: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read in the sources: B.md B.1 (items 5, 8), B2.0-B2.4 (F steps 1-5, 4b.1-4b.5, eval rule, beta), B3.1-B3.10,
B.6.2 (Gates B2, B3), B.9 rows U51, U70, U76, I48, G05; C.md C5-C11, C.6, C.9.1 "Weight push (C5)" (U-SY1..U-SY8);
00_front sections 2.5, 3 and 6.3 (V3 table, V3 amount rule); 99_back "RL settings" and "Model". Review log read
(review-log.md: front, B, C, back up to round 10) and `reviews/final_r10_rl.md` with its verdict. Round-10 items
(rl-1 refuted; rl-2..rl-5 minor) and the carried round-9 minors are not repeated.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL/vLLM sleep mode, LoRA serving or GRPO IS correction. No item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:351-362` (`is_ratio`), `:495-504` (`sync_without_merge`), `:524-616` (C10 skip, `_compute_loss`, `set_is_ratio`, `log`), `:618-690` (`_generate`, `time_vllm_phases`, `patch_vllm_prompts`), `:693-717` (`lora_merged_params`); `configs.py:177-195`. This gave item 1 (today's code only swaps `_iter_named_params`, so TRL's own `sync_weights` still clears its flag; C5 replaces `sync_weights` itself). |
| engineering:debug | Traced the mode-B engine state and TRL's `_llm_weights_sleeping` flag through start-up, the F comparator, the restore and step 1, for each way an engineer can read "if sleeping". Only one of the three readings works (item 1). |
| engineering:system-design | Checked who owns the "asleep" state: TRL's flag (set at `vllm_generation.py:379, 503, 773`) versus the replaced `sync_weights` and the comparator's raw `llm.sleep(level=2)`. Two owners, no rule (item 1). |
| engineering:architecture | Checked that C5 "Interface", C5 F sync step 1, B2.2 F step 2.1, 4b.3, 4b.5, U-SY3, U-SY8, U70 and U76 agree on the sleep/wake state. None names the flag (item 1). |
| engineering:documentation | Each item names its source file and heading, quotes the doc line, and cites the code line that proves it. |
| huggingface-skills:huggingface-best | Re-read front section 3 (ranking, size limit, 4B gate, 4B SFT step 0). The parameter counts and licences were re-checked in round 10 with the HF API; nothing in section 3 changed since. No model-choice item. |
| huggingface-skills:trl-training | Read TRL 1.17.1 in the repo .venv: `generation/vllm_generation.py:356-379` (start-up sleep, flag True), `:494-528` (`sync_weights`: wake, flag False at :503, per-parameter push), `:599-603` (generate re-pushes when the flag is True), `:718-719`, `:771-773` (sleep, flag True); `trainer/grpo_trainer.py:1186-1188` (`max_num_seqs` 256), `:1372-1379` (RepeatSampler), `:1838-1841` (one `sync_weights` per global step), `:2603-2611` (truncation mask), `:2776-2824` (IS modes), `:3066-3089` (OPSM), `:3190-3290` (loss, KL, dr_grpo normaliser); `trainer/grpo_config.py:322-360, 1160-1185` (IS keys and checks; `use_bias_correction_kl` default True). |
| huggingface-skills:huggingface-papers | No new paper claim in scope. The EGM numbers (2601.13633) were verified in round 7. |
| Web (WebSearch/WebFetch loaded) | Not needed for a new item: every claim in item 1 is in the pinned TRL source in the repo .venv. |

## Verified, no item

- IS modes: TRL's `sequence_mask` with `clip_min` None zeroes only ratios above 3 and keeps small ratios as small
  weights [code: grpo_trainer.py:2803-2818]; `fast.is_ratio` copies it [code: fast.py:351-362]. B.1 item 5 holds.
- OPSM sign and formula in B2.1 match `get_off_policy_mask` (`kl_div = sampling - per_token`, NaN to 0, keep if
  `adv >= 0` or mean KL <= delta) [code: grpo_trainer.py:3078-3088]. TRL reads `self.off_policy_mask_threshold` at
  each loss call (set at :882, read at :3190), so the live delta of B2.4 rule 4 works.
- KL term: with `use_bias_correction_kl` True (TRL default), TRL multiplies `per_token_kl` by `coef_1`, the policy
  ratio, not by the vLLM IS weight [code: grpo_trainer.py:3218-3225, 3268-3272]. Under C10 `coef_1` is 1 in value.
  The term stays linear in beta, so the B3 scaling `beta = 1e-3 / r`, `adam_epsilon = 1e-8 / r` keeps Adam's step
  (B2.4, B3.4, B3.7 rule 6).
- RLOO: `A_i = R_i - (sum R - R_i)/(G-1) = G/(G-1) (R_i - mean)`; group mean 0; `scale_rewards="none"` only centres.
  The B3.4 worked check (+0.057) recomputes: 0.1 x (0.62 - (0.29 + 6 x 0.014)/7) = 0.0567.
- Reward v5 monotonicity: d/dIoU of IoU + alpha (1 - IoU/m) g loc = 1 - alpha g loc / m >= 0 when alpha <= m (B3.3
  constructor rule); continuous at IoU = m.
- `max_num_seqs` = 16 x 1 x 16 = 256 [code: grpo_trainer.py:1186-1188] (C6 step 2).
- dr_grpo normaliser `size(0) x max_completion_length` [code: grpo_trainer.py:3286] (B3.5, C7 step 2).
- Model choice (front 3): ranking, size rule, 4B gate on val and the 4B SFT step 0 are consistent with 99_back
  decisions 17-19.

## Items

### 1. MAJOR (B, C): no rule for TRL's `_llm_weights_sleeping` flag once C5 replaces `sync_weights`; one of three readings works

- Sources: C.md "C5. Weight push" (Interface; "Under option F" step 1); B.md "B2.2 Numerics presets" (Option F step
  2.1; step 4b.3; step 4b.5).
- Doc lines:
  - C5 Interface: "`push_weights(vg, model, values="merged")` ...: if sleeping, `empty_cache()` and
    `wake_up(tags=["weights"])` (TRL's order, vllm_generation.py:500-503)" and "`install(trainer) -> bool`: replaces
    `vg.sync_weights` on the instance with `push_weights` (merged presets) or `sync_lora` (`lora-seqmask`)".
  - B2.2 F step 2.1: "If sleeping: `empty_cache()` and `wake_up(tags=["weights"])` (as TRL does [code:
    vllm_generation.py:500-503])."
  - B2.2 F step 4b.5: "... `reset_prefix_cache()`; then, in mode B only, `sleep(level=2)` last."
  - The doc never says how "if sleeping" is decided, and never names TRL's flag (grep of all sources for
    `_llm_weights_sleeping`, `is_sleeping`: 0 hits).
- Code that proves it (TRL 1.17.1, repo .venv):
  - `vllm_generation.py:376-379`: start-up `self.llm.sleep(level=2)`, then `self._llm_weights_sleeping =
    self.enable_sleep_mode`.
  - `:500-503`: TRL's own `sync_weights` wakes the weights and sets `self._llm_weights_sleeping = False` (line 503).
  - `:602-603`: `generate()` calls `self.sync_weights()` again when `self._llm_weights_sleeping` is True. Through
    `install`, that call is the replaced `push_weights` / `sync_lora`.
  - `:771-773`: after each generate, `sleep(level=2)` and the flag back to True.
  - `grpo_trainer.py:1838-1841`: the trainer already calls `sync_weights()` once per global step before generate.
  - Today this works only because `sync_without_merge` swaps `_iter_named_params`, not `sync_weights`, so TRL's own
    line 503 still runs [code: fast.py:495-504].
- What happens in mode B (sleep on: the C fallback, V2 item 3c `c_speed_b.yaml`, and the F run), by reading:
  - (a) The replaced sync does not touch the flag. The flag stays True from the last sleep. Each step then runs
    the sync twice: once from the trainer (:1840), once from `generate` (:603). Under F, `add_lora` runs twice per
    RL step, so the F check's "LoRA load count per RL step is 1" fails at step 1. `train_rl` raises, and B2.4 rule 0
    restarts under the fallback: one start-up (about 6 min) plus GPU time for nothing, and F, the expected preset,
    is dropped with no evidence against it. Under the merged mode-B arm, each step pushes twice (up to 0.2 s more
    by G-SP2, plus a second wake), which eats into the 2.7 s mode-B margin of G-SP1. `time/sync_s` is a mean per
    call, so G-SP2 does not show it.
  - (b) The replaced sync copies TRL's lines 500-503 as cited, including line 503, and reads "if sleeping" from the
    flag. Arm (ii)'s `sync_lora` then sets the flag False. The restore's raw `sleep(level=2)` (4b.5) does not set it
    back. At step 1 the trainer's sync sees the flag False, skips `wake_up(tags=["weights"])`, and `load_weights`
    writes into the slept (unmapped) weight memory. `generate` also skips its re-sync (flag False). Result: a CUDA
    fault or garbage weights at step 1, and the F check fails or the run crashes.
  - (c) The sync reads vLLM's own sleep state and also clears TRL's flag: this works. The doc does not say it.
- Caught before spend? No planned $0 check covers it. U-SY3 and U-SY8 test one sync call in isolation. U70 counts
  `add_lora` "once per sync", not once per RL step. U76's "2 `add_lora` calls in one RL step" is a planted failure
  that tests the check, not the real path. U76's fake "base tensors raise when read while asleep", but reading (b)
  writes, not reads. None of them runs TRL's `generate()` with its flag.
- Fix:
  - C5 Interface, after "if sleeping": "TRL's flag `vg._llm_weights_sleeping` is the one sleep state. `push_weights`
    and `sync_lora` wake only when it is True, and set it False after the wake (TRL line 503). Otherwise TRL's
    `generate` re-syncs (vllm_generation.py:602-603), and every step pushes twice."
  - B2.2 F step 4b.5: "then, in mode B only, `sleep(level=2)` last, and set `vg._llm_weights_sleeping = True` (as TRL
    does at vllm_generation.py:772-773). Step 1's sync then wakes the weights before it pushes `p.data`."
  - B2.2 F step 4b.3: "Before each comparator generation, wake the weights only through `sync_lora` /
    `push_weights` or with the flag update above. Never call a raw `wake_up(tags=["weights"])` without setting the
    flag."
  - Tests: U-SY3 and U-SY8 assert the flag is False after a mode-B push. New case in U70: wrap a real TRL
    `VLLMGeneration` around the fake `LLM` (patched `trl.generation.vllm_generation.LLM`). Run 3 trainer steps in
    mode B. Assert exactly 1 `load_weights` and 1 `add_lora` per step. A planted sync that leaves the flag True gives
    2 and fails. U76 mode B: after the comparator restore, `vg._llm_weights_sleeping` is True. The fake engine raises
    on `load_weights` while asleep, and the good run passes. C.9.1 G-SP2 note: also log `time/sync_calls` per step
    (must be 1).

### 2. MINOR (B): `completions/truncated_grad_share` leaves out the IS weight, so it overstates what truncated rows carry

- Source: B.md "B3.5 Wire v5 into the TRL adapter" ("`completions/truncated_grad_share` = sum over truncated rows of
  |A| x tokens / sum over all rows of |A| x tokens ... If the median exceeds 10 % over steps 50-260, the owner may
  switch on C11's option").
- Problem: B3 keeps `bf16-seqmask`. Its weight is exp(sum of token gaps), and B.1 item 5 says "long answers pass the
  cap first". A summed gap of -3 gives weight 0.05. So the 1,536-token truncated rows are the rows most likely to be
  zeroed or down-weighted. The share ignores this and can report 10 % while the weighted share is near 0. This
  feeds an optional owner switch and a report line, not a gate. No wrong gate result: minor.
- Fix: define the share with the final weight: sum over truncated rows of |A| x w x tokens / the same over all rows
  (w = `importance_sampling_ratio`, broadcast for sequence modes). Also log `completions/truncated_is_masked_share`
  (truncated rows with w = 0) beside it. The run report prints both.

### 3. MINOR (B): the OPSM delta callback has no trainer handle

- Source: B.md "B2.4 Decide", rule 4 ("a `TrainerCallback.on_step_end` at step 20 sets
  `trainer.off_policy_mask_threshold = delta`").
- Problem: the HF Trainer passes `model`, `optimizer`, `lr_scheduler`, dataloaders and `processing_class` to
  callbacks, not the trainer. So the callback needs a reference given at build time. An engineer finds this at $0
  when U51's live-delta test runs.
- Fix: "`OpsmDelta(trainer)` holds the trainer reference, set in `make_trainer` after the trainer is built
  (`trainer.add_callback(OpsmDelta(trainer))`)."

### Carried (not re-counted)

Round-10 rl-2 (eval LoRA id 1), rl-3 (seed s = 0 for all 8 copies), rl-4 (dead G-SP1 "Under F" margins: still in
C.md "C.6" G-SP1 row) and the round-9 minors rl-3..rl-8 are still in the text as filed. Their fixes are in
`final_r10_rl.md` and `final_r9_rl.md`.

## Counts

blocking 0, major 1 (item 1), minor 2 (items 2-3).
