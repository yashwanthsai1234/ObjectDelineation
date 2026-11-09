# Final review round 12: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read in the sources: B.md B.1 items 5 and 8, B1.5, B2.0-B2.4 (F steps 1-5, 4b.1-4b.5, eval rule, beta),
B3.1-B3.10, B.6.2 (Gates B2, B3), B.9 rows U51, U68, U70, U73, U76, I48, G05; C.md C5-C11, C.6, C.7, C.8.1, C.9.1
"Weight push (C5)" (U-SY1..U-SY9); A.md A7 (contract list and the "later phases" rule); 00_front sections 2.5 and 3;
99_back "RL settings" and "Model". Review log read (review-log.md, front/B/C/back, up to round 11) and
`reviews/final_r11_rl.md` with its verdict. Round-11 rl-1 (sleep flag) is fixed in C5 "Sleep state", B2.2 F 2.1 /
4b.3 / 4b.5, U-SY3, U-SY8, U-SY9, U70, U76 and G-SP2; I checked the fix against TRL's code and it holds (below).
Round-11 rl-2/rl-3 and the carried round-9/10 minors are not repeated.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL GRPO, vLLM LoRA serving or sleep mode. No item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:351-362` (`is_ratio`), `:495-504` (`sync_without_merge`), `:524-616` (C10 skip, `_compute_loss`, `set_is_ratio`, `log`), `:693-717` (`lora_merged_params`); `configs.py:1-30, 90-196`; `reward/advantages.py` (v4, legality gate on `r_topo`); `reward/score.py:24-27`; `train/rewards.py:26-28, 68-76` (truncated = format failure). Confirmed the B3 reward/RLOO claims and the C10/OPSM interplay. |
| engineering:debug | Traced the round-11 sleep-flag fix through TRL's start-up sleep (vllm_generation.py:376-379), the step-start sync (grpo_trainer.py:1838-1841), `generate`'s re-sync (:602-603) and the end sleep (:771-773), in mode A and mode B, for merged presets and F. The fixed rule gives exactly one sync per step in both modes. |
| engineering:system-design | Checked one owner per state: TRL's flag is now the one sleep state; the F comparator restore sets it; `push_weights` / `sync_lora` clear it. Checked one home per name for the private-name contract (item 1). |
| engineering:architecture | Checked that the B2.4 rules (start/restart preset), C.6 "Numerics first", 2.5's diagram and 99_back decisions 15-16b agree: B1 and B3 on `bf16-seqmask`, F only after Gate B3. They do. |
| engineering:documentation | Each item names its source file and heading, quotes the doc line and the code or web line. |
| huggingface-skills:huggingface-best | Re-checked front section 3 against the HF API on 2026-11-09: Qwen3-VL-2B and -4B configs (vision depth 24, width 1024, DeepStack 5/11/17, merger out 2,048 vs 2,560; 4B 4,437,815,808 params, Apache-2.0) match the table. Listed recent image-text-to-text models and recent Qwen uploads: the new small ones are derivatives of parked families (`LiquidAI/d1-2B` from LFM3-VL-2B; `Qwen/Qwen-Drive-1.0-4B` from Qwen3.5-4B, 4.54 B), and no new Qwen VLM up to 4B exists (Qwen3.6/3.8 are 27B+). No model-choice item. |
| huggingface-skills:trl-training | Read TRL 1.17.1 in the repo .venv: `grpo_trainer.py:2603-2611` (truncation mask), `:2776-2824` (IS modes), `:2870-2930` (`sum_then_normalize`, `scale_rewards="none"`, unscorable rows), `:3066-3089` (OPSM), `:3180-3290` (OPSM call, KL with `use_bias_correction_kl`, IS weight applied before the KL term, dr_grpo normaliser); `generation/vllm_generation.py:356-382, 490-530, 595-610, 700-780`. |
| huggingface-skills:huggingface-papers | No new paper claim in scope; EGM (2601.13633) was verified in round 7. |
| Web (WebSearch/WebFetch loaded; raw GitHub at tag v0.32.0) | `vllm/lora/worker_manager.py:274-324` (LRU `add_adapter`: `load_inplace` reloads and re-activates the slot; `_apply_adapters` calls `add_adapter` per request, so a request with `load_inplace=True` reloads every engine step) confirms B2.2 F step 3 and C5 step 3. `vllm/v1/core/kv_cache_utils.py:578-594` gave item 2. |

## Verified, no item

- Round-11 fix: TRL's own `sync_weights` wakes whenever sleep mode is on and sets the flag False
  (vllm_generation.py:500-503); `generate` re-syncs only when the flag is True (:602-603); the end sleep sets it True
  (:771-773). C5 "Sleep state" copies this, and B2.2 F 4b.5 sets the flag after the comparator's raw sleep. U-SY9
  runs TRL's real `VLLMGeneration` with planted failures in both directions. One sync per step in both modes.
- RLOO through TRL: with `sum_then_normalize` and `scale_rewards="none"`, TRL computes `rewards - group mean` and
  never divides (grpo_trainer.py:2880-2907); the weight-0 loggers drop out through `nansum`. A v5 aggregator output
  with group mean 0 passes unchanged (B3.4).
- `r_topo` is gated by `legal` (advantages.py `topology_valid` -> `legal_gate`), so B3.4's
  `R = legal * (...) + c_topo * r_topo` gives an illegal answer R = 0 (T-tests hold).
- Truncated answers are format failures in the reward (rewards.py:75-76, score.py:24-27), so B3.5's "A = -b" for an
  unmasked truncated row and C11's "capped row gets R = 0" hold.
- B.1 item 8: E[1{nt}(R - b) grad log pi] = grad J - b grad P(nt), using R = 0 when truncated. Correct.
- KL: TRL multiplies `per_token_loss` by the vLLM IS weight before it adds `beta * per_token_kl`
  (grpo_trainer.py:3268-3272), and the KL is weighted only by `coef_1` (1 in value under C10). So the B3 scaling
  `beta = 1e-3 / r`, `adam_epsilon = 1e-8 / r` keeps v4's policy-to-KL balance and Adam's step (U71). AdamW weight
  decay is decoupled and 0 by default, so it does not break the argument.
- OPSM reads `self.off_policy_mask_threshold` at each loss call (grpo_trainer.py:3190) and prefers
  `sampling_per_token_logps`; under C10 `old_per_token_logps = per_token_logps.detach()` (:3187-3188). B2.1's sign
  and B2.4 rule 4 match.
- dr_grpo normaliser `per_token_loss.size(0) * max_completion_length` (:3286): B3.5's "12x gradient mass"
  (1,536 / 131 = 11.7) and C7 step 2's 1,536 / L scale hold.
- Adapter size 622,592 x 28 = 17,432,576 (C5), the 4b.4 timing arithmetic (736 x 4.3 ms / 11 s = 0.29; 5.61 s +
  4.3 ms x 512 = 7.8 s), and the B3.4 RMS arithmetic (0.0048 / sqrt(17,432,576) = 1.15e-6) recompute.
- vLLM 0.32.0 LoRA: `add_lora(..., load_inplace=True)` re-loads id 1 and re-activates its GPU slot, which also
  restores the slot after level-2 sleep (worker_manager.py:296-324). The comparator's ids 1/2 fit the default
  `max_loras` 1 (LRU eviction; only one id is active per batch).
- `shuffle_dataset=False` (B1.5) makes "rows 0-31 of the `--train` root" the step-1 prompts (B2.2 F step 4).

## Items

### 1. MINOR (A, B, C): `add_lora` has two homes in the private-name contract, and the round-11 flag has none

- Sources: A.md "A7" ("B2.2 option F calls `llm.llm_engine.add_lora` once per sync. That is vLLM, not TRL, and vLLM
  is not on the Mac. So B adds it to the vLLM list of `env_check.sh` on the box, not to `ALL_ITEMS`."); C.md "C5"
  Interface ("Add `LLMEngine.add_lora` to the A7 start-up check"); B.md "B2.2" F step 3 ("The patch points (`LLM`
  kwargs, `patch_vllm_prompts`, `llm_engine.add_lora`) go into the Phase A7 private-name start-up check").
- Problem: (a) A7 says `add_lora` goes to `env_check.sh`, not `ALL_ITEMS`; C5 and B2.2 say the A7 check. If an
  engineer puts it in `ALL_ITEMS`, `check()` imports vLLM on the Mac and `make test` fails at $0, so this is caught
  before spend. (b) The round-11 fix makes C5 read and write `VLLMGeneration._llm_weights_sleeping`, a private TRL
  attribute (vllm_generation.py:379, 503, 773). A7's rule ("the commit that starts to use a new TRL name adds it to
  the contract in the same commit") covers it, but A7's list of planned names and C5's file list do not name it.
  On a TRL bump that renames it, `push_weights` would write a new attribute that TRL never reads, and every step
  would sync twice; U-SY9 catches that only if it runs against the new TRL.
- Fix: C5 Interface: "Add `LLMEngine.add_lora` to the vLLM list of `env_check.sh` (A7 rule), and add
  `_llm_weights_sleeping` to `trl_adapter.PRIVATE_ATTRS` (substring `self._llm_weights_sleeping =` in the
  VLLMGeneration source) in the C5 commit." B2.2 F step 3: "`LLM` kwargs and `patch_vllm_prompts` go into the A7
  contract; `llm_engine.add_lora` goes into `env_check.sh` (A7)." A7 planned-names list: add a third bullet for
  `_llm_weights_sleeping` (C5, PRIVATE_ATTRS).

### 2. MINOR (B, C): the stated reason for `reset_prefix_cache()` after `add_lora` is incomplete for vLLM 0.32.0

- Sources: B.md "B2.2" F step 2.4 ("`llm.reset_prefix_cache()`: the prefix-cache key holds the LoRA name, not its
  weights"); C.md "C5" F sync step 4 (same words).
- Web: vLLM v0.32.0 `vllm/v1/core/kv_cache_utils.py:578-594`, `_gen_lora_extra_hash_keys` returns
  `[("lora", lora_request.lora_name, lora_request.lora_path)]` ("The adapter path is included so that re-pointing a
  LoRA name at a different adapter does not reuse KV computed with the previous one").
- Problem: the key holds the name AND the path. The reset is still needed, because every F sync writes the new
  adapter into the same `<tmpfs dir>`, so name and path stay equal while the weights change. The conclusion holds;
  only the reason is incomplete, and a reader could think a new folder per sync removes the need (it would, but the
  plan does not do that). No wrong result.
- Fix: B2.2 F step 2.4 and C5 F step 4: "`llm.reset_prefix_cache()`: the prefix-cache key holds the LoRA name and
  path, not its weights [web: vllm v0.32.0 v1/core/kv_cache_utils.py:578-594]; every sync reuses the same
  `<tmpfs dir>`, so blocks from the old adapter would match."

### 3. MINOR (B): the active `ess` has no stated unit, so it cannot be compared with `cf/ess_C` of the same mode

- Source: B.md "B2.1" ("`ess`: effective sample size `(sum w)^2 / (N sum w^2)` of the active weights"; "counterfactual
  ESS on one unit, token level (a sequence weight is broadcast to its tokens): `cf/ess_C` (sequence_mask, cap 3)").
- Problem: in the B1 run the active mode is `sequence_mask` cap 3, the same mode as `cf/ess_C`. Today's
  `set_is_ratio` counts sequence-mode ratios per sequence (`rk` is all-ones of shape (B, 1) [code: fast.py:594]), so
  the natural build gives a sequence-level `ess` and a token-level `cf/ess_C`, two different numbers for one mode.
  Both are report lines (rule 3 reads only `cf/ess_E@0.1`), so no gate changes.
- Fix: B2.1: "`ess` uses the same token-level unit as the `cf/ess_*` keys (a sequence weight broadcast to its valid
  tokens). Under `bf16-seqmask`, `ess == cf/ess_C` (test in `tests/test_train_fast.py`)."

## Counts

blocking 0, major 0, minor 3 (items 1-3).
