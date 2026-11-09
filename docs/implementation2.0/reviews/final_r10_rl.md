# Final review round 10: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read in the sources: B.md B.1 item 5, B1.5, B2.0-B2.4 (with the round-9 text of F step 4, 4b.1-4b.5 and the
`lora-seqmask` eval rule), B3.1-B3.10, B.6.2 (Gates B2, B3), B.9 rows U51, U70, U76, I48, G05; C.md C5-C8, C11, C.6,
C.8.1; 00_front sections 2.5 and 3; 99_back "RL settings" and "Model". Review log read: B, C, front and back rounds 8
and 9, and `reviews/final_r9_rl.md` with its verdict. Items fixed in round 9 (rl-1, rl-2) are not repeated.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL/vLLM sleep mode, LoRA serving or GRPO IS correction. No item depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:351-362` (`is_ratio`), `:495-504` (`sync_without_merge`), `:524-616` (`skip_old_logps`, `_compute_loss`, `set_is_ratio`, `log`), `:618-626` (`_generate`), `:670-690` (`patch_vllm_prompts`), `:693-717` (`lora_merged_params`); `configs.py:70-195`; `train_rl.py:19-68`; `final_eval.py:46-85` (it samples through `g_probe.make_trainer` and `trainer._generate`, so TRL's vLLM path). These reads gave items 1 and 2. |
| engineering:debug | Traced the engine state of a mode-B F run from `LLM(...)` to the first comparator generate. TRL puts the engine to level-2 sleep right after it builds it, and the comparator only wakes it. No step re-pushes the base before arm (ii). This is item 1. |
| engineering:system-design | Checked the LoRA slot ids across the three users of `enable_lora`: the F sync (`"policy"`, id 1), the comparator (`"b1"`, id 2) and `final_eval --adapter` (`"eval"`, id 1). That gave item 2. |
| engineering:architecture | Checked that B2.2 F step 4b, C5 step 5(c), G-EQ5, U70 and U76 agree on comparator order and adapter names. They agree. None of them pushes the base before the first arm in mode B (item 1). |
| engineering:documentation | Each item names its source file and heading, quotes the doc line, and cites the code line, web source or measurement that proves it. |
| huggingface-skills:huggingface-best | Re-checked front section 3 with the HF API. Qwen3-VL-4B-Instruct 4,437,815,808, Qwen3-VL-2B-Instruct 2,127,532,032, Qwen3.5-4B 4,659,865,088 parameters, all Apache-2.0. 4B vision: depth 24, width 1024, DeepStack 5/11/17, `out_hidden_size` 2560 (2B: 2048). Qwen3.5-4B `layer_types` has `linear_attention` 3 of every 4 layers. Every front claim holds: no model-choice item. |
| huggingface-skills:trl-training | Read TRL 1.17.1 in the repo .venv: `generation/vllm_generation.py:356-379` (start-up `sleep(level=2)` and the `_llm_weights_sleeping` flag), `:497-503`, `:599-603` ("waking up isn't enough: they must be re-pushed"), `:718-719`, `:771-773`; `trainer/grpo_trainer.py:1029-1036`, `:1332-1377` (RepeatSampler: step 1 = rows 0-31 with `shuffle=False`), `:2760-2824` (IS modes), `:2880-2907` (`sum_then_normalize` with `scale_rewards="none"` only centres), `:3066-3089` (OPSM), `:3195-3290` (IS weight, then the unweighted `beta * per_token_kl`, then the dr_grpo normaliser). |
| huggingface-skills:huggingface-papers | No new paper claim in scope. The EGM numbers (2601.13633) were verified in round 7; I did not fetch them again. |
| Web | vLLM sleep-mode docs (docs.vllm.ai/en/latest/features/sleep_mode): "Level 2 sleep will discard the model parameters and the KV cache"; "The content of the parameters and KV cache is forgotten"; after waking, restore them before you generate. This gives item 1. |

## Verified, no item

- Measured values recomputed from `runs/rl/train_rl/rl.json` (770 logged steps): step-1 logp difference 0.00956
  (also 0.00956 in `runs/r5/base/rl.json`, 0.00952 in `runs/r5/shared_rl.json`); median over steps 240-280 is
  0.0176, over 101-300 0.0236, over 301-500 0.0185, over 1-20 0.0110, over 21-100 0.0278. IS min = 0 in 98.05 % of
  steps. Median step-mean ratio 0.823. grad_norm 0.0018-0.0221, median 0.00475. All match B2.2 F step 4, B2.4
  rule 4 and B3.4.
- RLOO through TRL: `rewards - mean_grouped`, with no division when `scale_rewards == "none"`
  [code: grpo_trainer.py:2904-2906]. So the group-mean-0 RLOO advantage passes through unchanged (B3.4).
- The beta and eps scaling (B2.4, B3.4, B3.7 rule 6) holds. TRL adds `beta * per_token_kl` after the IS weight,
  and that term carries no weight [code: grpo_trainer.py:3268-3272]. So (PG/r + (beta/r) KL) with eps/r gives Adam
  the step of (PG + beta KL, eps). AdamW `weight_decay` is 0 by default, so no other term scales.
- The step-1 prompts of the F comparator (round-9 fix) are right. `RepeatSampler(batch_size=256 // 8 = 32,
  shuffle=False)` [code: grpo_trainer.py:1372-1377] gives rows 0-31 of the file at step 1.
- Block rule: the C arms (20), the fp16 arm (60), B1 and B3 (260) all divide by the 20-step block (B1.5).
- Masked-truncation bias: E[1(nt)(R - b) grad log pi] = grad E[R] - b grad P(nt), because R = 0 when truncated. This
  is the formula in B3.7 rule 9 and C7 "Facts".

## Items

### 1. BLOCKING (B; C follows): in mode B, comparator arms (ii) and (0) and the timing arms decode from weights that TRL's start-up level-2 sleep discarded. No base push comes before them.

- Sources: B.md "B2.2 Numerics presets", Option F step 4b.3 and 4b.4; B.md B.9.1 U76; C.md "C5" step 5(c) and C.6
  G-EQ5 (they point to 4b).
- Doc lines:
  - 4b.3: "Engine state. Before each comparator generation: `wake_up(tags=["weights"])` if asleep, then
    `wake_up(tags=["kv_cache"])`. In mode B, TRL sleeps the engine at start-up and wakes the KV cache only inside its
    own generate".
  - 4b.4: "Order. Arm (ii) runs first, on the start-up base W. The floor arm (0) runs next, on the same base W. The
    B2.4 rule 0 timing runs after them ... Arm (i) runs last: its merged push writes `W + d_B1` into vLLM's base
    tensors".
  - F sync step 2.2: "Mode B (sleep on, today): push the UNMERGED base weights ... The one exception is the
    comparator (step 4b): it pushes `W + d` once, so it restores the base once, before step 1, in both modes." The
    doc pushes nothing before arm (ii).
- Proof:
  - TRL builds the engine, then at once calls `if self.enable_sleep_mode: self.llm.sleep(level=2)`, with the comment
    "Sleep level 2 discards the weights; track it so that generate() knows it must re-push them"
    [code: .venv trl/generation/vllm_generation.py:376-379]. `grpo_config` sets `vllm_enable_sleep_mode: True`
    [code: src/objdel/train/configs.py:193].
  - TRL's own generate does not rely on a wake. "Sleep level 2 discards the weights, so waking up isn't enough:
    they must be re-pushed from the training model", and then `if ... self._llm_weights_sleeping:
    self.sync_weights()` [code: vllm_generation.py:599-603].
  - vLLM docs: "Level 2 sleep will discard the model parameters and the KV cache"; "The content of the parameters
    and KV cache is forgotten"; restore them after the wake [web: docs.vllm.ai/en/latest/features/sleep_mode].
  - The comparator runs before step 1, so no `sync_weights` has run yet. It calls `llm.generate` (through
    `vllm_generate_checked`), not TRL's generate, so TRL's re-push guard never fires.
  - Effect: in mode B, arm (ii) (`"b1"` LoRA on base W), the floor arm (0) (base W) and the timing runs sample
    from weight memory that holds no model. f_lora and f0 are then garbage. Expect NaN, or a comparator ratio far
    above 0.5, so the F check fails.
  - B2.4 rule 0 then restarts the run under the fallback preset. That costs one start-up (about 6 min) plus 15 steps
    of GPU money for nothing. It also drops F, the expected preset (seed error 0.531), with no real evidence against
    it.
  - Only arm (i) works, because its merged push writes every tensor.
  - Mode B is today's default and the fallback whenever G-MEM fails (C5 step 3). Mode A (sleep off) is not
    affected.
  - U76's fake models the sleep only as "base tensors raise when read while asleep". It has no "woken from level 2
    but never pushed" state, so the test cannot catch this.
- Fix:
  - B2.2 F step 4b.3: add "In mode B the engine has been at level-2 sleep since start-up, so its weights are gone
    (vllm_generation.py:376-379). Before arm (ii), the comparator wakes the weights and pushes the unmerged base
    once: `push_weights(vg, model, values="base")` (the F sync step 2.2 call, values `p.data`). Then it hashes every
    vLLM base tensor on the host against `p.data` (the C5 step 5(b) method). Then it wakes the KV cache. Mode A needs
    no push."
  - 4b.4: replace "on the start-up base W" with "on the base W that step 4b.3 pushed (mode B) or loaded at start-up
    (mode A)".
  - F sync step 2.2: "The comparator pushes the base twice in mode B (before arm (ii) and at the restore). It pushes
    the base once in mode A (at the restore)."
  - U76 and U70: the fake engine's base tensors read as NaN after a level-2 wake until `load_weights` has written
    them. In mode B, the first comparator `generate` must come after one `load_weights` call whose values are
    `torch.equal` to `p.data`. A planted run that skips the push fails with `comparator.ok` false.
  - C5 step 5(c) and G-EQ5: add "mode B: also before arm (ii)" beside the restore hash.

### 2. MINOR (B): `final_eval --adapter` reuses LoRA int id 1, the id of the F sync's `"policy"` slot

- Sources: B.md "B2.4 Decide", eval rule ("passes `LoRARequest("eval", 1, <dir>)` to every generate"); B.md "B2.2" F
  step 2.3 (`LoRARequest("policy", 1, <tmpfs dir>, load_inplace=True)`); I48.
- Problem: vLLM keys a loaded adapter by `lora_int_id`. A request with a new name and an id that is already loaded
  is served by the adapter already in that slot. `final_eval` builds a trainer through `g_probe.make_trainer`
  [code: experiments/phase2/final_eval.py:72-75]. If an engineer builds it with the `lora-seqmask` sync installed,
  `sync_lora` first adds the trainer's fresh, zero adapter as id 1. Every eval generate then samples the base
  model and gives no error. With today's merged sync this does not happen. The doc does not say which sync
  `final_eval --adapter` installs, so there is no proof of a wrong result: minor.
- Fix: eval rule: "`LoRARequest("eval", 3, <dir>)` (id 3: ids 1 and 2 belong to the F sync and the comparator);
  `final_eval --adapter` keeps the merged sync, whose push of a fresh adapter writes base W exactly". I48 asserts id
  3.

### 3. MINOR (B): "seed s = 0" for all 256 comparator rows gives 8 identical answers per prompt

- Source: B.md "B2.2" F step 4 ("each sent 8 times in the run's own request form (256 answers ...), seed s = 0 for
  arms (0), (i) and (ii)") and 4b.4 (`SamplingParams(..., seed=s)`).
- Problem: one per-request seed on 8 copies of the same prompt draws the same sample 8 times. So the 32 x 8 batch
  holds 32 distinct answers, not 256. The comparator means stay unbiased; only their noise grows. The timing arms
  are not affected (fixed length).
- Fix: "seed = 1000 x s + j for copy j = 0..7 of each prompt, the same in every arm (common random numbers across
  arms, 256 distinct draws)".

### 4. MINOR (C): G-SP1's "Under F" margins can never be applied

- Source: C.md "C.6" G-SP1 row ("Under F: `- 2.7 s` (A) or `- 2.1 s` (B) ... under F `- 1.5 s` (A) or `- 1.25 s` (B)")
  against the same section: "V2 has no F arm", "Under F (runs after Gate B3) no per-parameter TRL path exists, so
  there is no F baseline arm", and "Under F that run reports its step time against the B3 run as a report line;
  G-SP1 is not re-gated".
- Fix: delete the F margins from G-SP1, and the "under F" rows of C.8.1 that feed only them. Keep one sentence: "Under
  F, step time is a report line against the B3 run (C.6 Numerics first)".

### 5. Carried round-9 minors (still open in the text; not re-counted)

Round-9 rl-3 to rl-8 were confirmed minor, and the round-9 fixer did not process them. They are still in the text:
B.md:134 ("so at least 1 of 256 sequences was masked", rl-6); B.md:1315 ("measures the kernel cost at batch 256",
rl-4); B.md:1460 ("Until step 20, delta = 0.5", rl-8); and rl-3, rl-5 and rl-7 as written in `final_r9_rl.md`. The
fixes are those of `final_r9_rl.md`.

## Counts

blocking 1, major 0, minor 3 (items 2-4; the carried round-9 minors are not counted again).
