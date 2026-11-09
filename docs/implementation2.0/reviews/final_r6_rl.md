# Final review round 6: RL (B2, B3, C, model choice)

Reviewer role: RL for LLM/VLM research scientist (GRPO, IS correction, reward design, TRL 1.17.1 and vLLM internals).
Scope read: B.md B1.5, B2.0-B2.4, B3.1-B3.10; C.md C5-C11, C.6; 00_front section 3; 99_back "RL settings" and
"Model"; review-log.md (all RL rows of rounds 1-5, so nothing below repeats a fixed item).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL/vLLM LoRA serving; no item depends on it. |
| engineering:code-review | Read `fast.py:488-620, 696-717` (C10 skip, `_compute_loss`, `set_is_ratio`, `lora_merged_params`) against the B2.2 F text; found that `lora_merged_params` uses the active adapters and that `get_peft_model_state_dict` reads `"default"` (item 2). |
| engineering:debug | Traced the F comparator step by step against TRL's sleep calls (`vllm_generation.py:376-377, 500-503, 718-719, 771-773`) and C5 mode A; found the doubled delta (item 1) and the asleep KV cache (item 2). |
| engineering:system-design | Checked which engine state each F check step needs (base W, adapter slot, KV awake) in mode A and mode B. |
| engineering:architecture | Checked ADR-0013 rule (b) scope against rules 0-4 and decision 15 (items 3, 5). |
| engineering:documentation | Each item names the source file + heading and the proving line. |
| huggingface-skills:huggingface-best | HF API `safetensors.total`: Qwen3-VL-2B 2,127,532,032; Qwen3-VL-4B 4,437,815,808; Qwen3.5-4B 4,659,865,088. Front section 3 numbers confirmed. |
| huggingface-skills:trl-training | Read TRL 1.17.1 `grpo_trainer.py:2780-2824` (vLLM IS modes), `3066-3089` (OPSM sign), `3150-3290` (loss order: IS ratio applied before `+ beta * kl`; dr_grpo constant), `2880-2905` (`sum_then_normalize` with `scale_rewards="none"`), `1025-1036` (no ref copy with PEFT). B3.4 and C7/C8 claims confirmed. |
| huggingface-skills:huggingface-papers | arXiv 2601.13633 (EGM) abstract read; Tables 4/9 not re-read (no item depends on it). |
| Web | arXiv abs page for 2601.13633. |

## Verified, no item

- B3.4: RLOO `A_i = R_i - (sum R - R_i)/(G-1)` has group mean 0, so TRL's `sum_then_normalize` + `scale_rewards="none"`
  returns it unchanged [code: .venv trl/trainer/grpo_trainer.py:2880-2905; configs.py:190-191].
- B2.4 KL: `beta = 1e-3/r` with `adam_epsilon = 1e-8/r` is consistent. KL enters `per_token_loss` before the dr_grpo
  constant [code: grpo_trainer.py:3271-3286], so the policy and KL parts scale together, and Adam is invariant to
  (g/r, eps/r) with weight decay 0.
- B2.1 OPSM sign: `seq_kl = mean(sampling - logp)` matches TRL [code: grpo_trainer.py:3078-3086].
- B2.2/C5: adapter count 17,432,576 = 28 x 622,592 [arithmetic]; checkpoint-100 header 392 keys,
  `base_model.model.model.language_model.layers.*.lora_[AB].weight`, F32, `r` 16, `lora_alpha` 32 [measured: header
  read 2026-11-08].
- C8: with PEFT TRL keeps no reference copy [code: grpo_trainer.py:1031-1036].

## Items

### 1. BLOCKING (B, C): under C5 mode A the F comparator leaves W + d in vLLM, so every later F step samples from W + 2d

- Sources: B.md "B2.2 Numerics presets", Option F plumbing step 2.2 and step 4; C.md "C5. Weight push", F sync step 2.
- Doc lines:
  - B2.2 F step 4: "(i) It pushes the adapter merged (`values="merged"`, C5), generates the 32 x 8 batch ... (ii) It
    runs `sync_lora` with the same adapter, generates the same batch".
  - B2.2 F step 2.2: "Mode A (vLLM awake, C5): no base push after the first sync."
  - C5 F sync step 2: "Mode A: no base push; vLLM loaded the same base checkpoint at start-up, and it never sleeps."
    C5 interface: `sync_lora` "calls `push_weights(..., values="base")` in mode B only."
  - B.8.4 "Faster step" row: "Under `lora-seqmask` (runs after Gate B3), C5's push is: mode A the adapter only".
    C.6: "Fast arm = C1-C5, mode A (vLLM awake) unless G-MEM fails". So mode A + F is a planned combination.
- Code: `values="merged"` pushes `lora_merged_params`, which yields `p.data + d.to(p.dtype)` for every LoRA layer
  [code: src/objdel/train/fast.py:716-717] into vLLM's `load_weights`, which overwrites vLLM's base tensors. In
  mode A nothing writes `p.data` back afterwards.
- Result: (ii) and every training step after it decode with base W + d_B1 plus the LoRA slot. During steps 1-15 the
  slot holds the run's own adapter, so the sampler is W + d_B1 + d_run, not the trainer's W + d_run. Two outcomes:
  the output-level check fails for a plumbing reason, so rule 0 wrongly rejects F, ADR-0013 records a false result,
  and the restart costs about 6 min + 15 steps; or it passes and the run trains with a wrong sampler, and the IS
  mask hides the gap. The G-EQ5 F check does not catch it, because in mode A it compares only the LoRA slices
  ("in mode B every base tensor hashes equal to `p.data`").
- Fix (B2.2 F step 4 and C5): after the comparator's merged arm, always restore the base, in both modes: call
  `push_weights(vg, model, values="base")` before (ii). Better still, run the merged arm last, then push the base
  before step 1. In mode A, also hash every base tensor against `p.data` once after the comparator (C5 step 5(c)
  and G-EQ5: "in both modes after the comparator, every base tensor hashes equal to `p.data`"). U70/U76: add a mode-A
  fake case where the comparator runs and the next `load_weights` call before step 1 carries `p.data`.

### 2. MAJOR (B, C): the F comparator and ABAB timing have no plumbing for the adapter, the request or the sleep state

- Sources: B.md "B2.2 Numerics presets" F steps 2.3, 4 and 5; B.md "B2.4 Decide" rule 0 (timing rule); B.md B.9 U70;
  C.md "C5" step 5(c).
- Doc lines:
  - B2.2 F step 4: "It loads the B1 final adapter (`final_adapter/` of HF revision `b1-rows`, B1.5) on `rl-full-3ep`.
    ... (ii) It runs `sync_lora` with the same adapter ... (iii) The trainer scores both answer sets with way (b)".
  - B2.2 F step 2.3: the adapter file comes "from `peft.get_peft_model_state_dict(model)`" and
    "`model.peft_config["default"].save_pretrained(<tmpfs dir>)`".
  - B2.2: "The `patch_vllm_prompts` wrapper ... passes `lora_request=LoRARequest("policy", 1, <tmpfs dir>)` ... to
    every `generate` call." U70: "every `generate` call carries `LoRARequest("policy", 1, <dir>)`".
  - Rule 0 timing: "after the first F sync, generate the same 32 x 8 batch twice, with the `lora_request` and with
    `lora_request=None`".
- Code and why the engineer must guess:
  - (a) Adapter. The training process holds the run's own fresh "default" adapter. `get_peft_model_state_dict` reads
    adapter `"default"` by default [code: .venv peft/utils/save_and_load.py], and `lora_merged_params` merges the
    *active* adapters [code: fast.py:703-707]. If B1 is loaded as a second adapter, `sync_lora` writes the fresh
    adapter (B = 0) while (i) and (iii) use B1. The check then compares base against B1 and fails. If B1 is loaded
    into "default", the run starts training from the B1 adapter instead of its own start state. The doc says
    neither how to load the adapter nor how to remove it before step 1.
  - (b) Sleep state. In mode B TRL puts the engine to sleep at start-up [code: .venv
    trl/generation/vllm_generation.py:376-377]. It wakes the KV cache only inside its own generate path (`wake_up(tags=["kv_cache"])`, :718-719)
    and sleeps at level 2 after each generate (:771-773). `push_weights` wakes only the weights (C5 interface). So the
    comparator's direct generations before step 1 run with the KV cache asleep. The ABAB timing "after the first F
    sync" has no stated place before TRL's level-2 sleep discards the weights. After that sleep it would time
    generation from freed memory.
  - (c) Request. Arm (i) must generate with `lora_request=None`, but the doc says the wrapper adds the request to
    "every" call, and U70 asserts this. Before the first `sync_lora` the tmpfs folder does not exist yet.
  - (d) Base. The comparator says "on `rl-full-3ep`". The first run after Gate B3 has no `--ckpt` stated anywhere. If
    it starts from another model, the trainer and vLLM in that process hold a different base.
- Fix (B2.2 F step 4, new sub-steps; mirror in C5 step 5(c)):
  1. Load B1 with `model.load_adapter(<b1 dir>, adapter_name="b1")` and `model.set_adapter("b1")`. Write the
     comparator's adapter folder with `get_peft_model_state_dict(model, adapter_name="b1")` and
     `peft_config["b1"]`. `sync_lora` gains an `adapter_name` argument (default `"default"`).
  2. Before each comparator generation: `wake_up(tags=["weights"])` if needed, then `wake_up(tags=["kv_cache"])`.
     Generate through a direct `llm.generate` that bypasses the wrapper. Pass `lora_request=None` for (i) and the
     comparator request for (ii). After the comparator: `set_adapter("default")`, `delete_adapter("b1")`, push the
     base (item 1), `remove_lora(1)`, then `sleep(level=2)` in mode B. This gives step 1 the same engine state as a
     run without `--f-check`.
  3. Run the ABAB timing inside the first F sync, right after `add_lora` and before TRL's generate, with the KV
     cache woken. Run it on the run's own adapter, and log its place in `f_check.json`.
  4. State "comparator base = the run's own `--ckpt`". If that is not `rl-full-3ep`, use the run's start adapter, or
     skip the output-level check with a recorded reason.
  5. U70: "every generate call *from TRL's generate path* carries the LoRARequest". U76: planted cases for a wrong
     adapter name and for a comparator that leaves "b1" active at step 1.

### 3. MINOR (B): rule 3 and rule 4 both set the starting OPSM delta

- Source: B.md "B2.4 Decide" rules 3 and 4.
- Rule 3: "if `cf/opsm_share@0.5` > 2 % of sequences, raise the starting delta to the q99.9 of `seq_kl` ... of B1 steps
  50-260". Rule 4: "Until step 20, delta = 0.5 (TRL docstring)".
- Fix: rule 4 should read "Until step 20, delta = 0.5, or rule 3's q99.9 when its veto fired". `bf16-tis-opsm`'s
  literal 0.5 in B2.2 is then replaced at start-up in the same way. Add one U51 case.

### 4. MINOR (B): the rule-4 callback cannot reach the trainer through HF's callback arguments

- Source: B.md "B2.4 Decide" rule 4: "a `TrainerCallback.on_step_end` at step 20 sets `trainer.off_policy_mask_threshold
  = delta`".
- HF passes `model`, `optimizer`, `lr_scheduler`, `processing_class` and the dataloaders to callbacks, not the trainer.
  `off_policy_mask_threshold` is a trainer attribute [code: grpo_trainer.py:882].
- Fix: "the callback is built with the trainer (`OpsmDelta(trainer)`) after `FastGRPO(...)`, then
  `trainer.add_callback(cb)`". The same applies to the device lists that `set_is_ratio` fills.

### 5. MINOR (C, 99_back): G-C's beta condition uses a value that no B3 default has

- C.md "C.6" G-C: "If the owner keeps beta 1e-3 (C8), G-C is reported as a computed number, not gated". 99_back
  decision 15: B3 default is "`1e-3 / r`" and "If beta stays above 0, C-G-C (18 s) is reported". G-C is measured
  in the first run after Gate B3, and no decision sets that run's beta.
- Fix: in C.6, write "If beta > 0 in the run that measures G-C". In decision 15, add "the runs after Gate B3 keep
  B3's beta unless the owner changes it".

### 6. MINOR (B): the B2.0 seed rounds differently from the code it models

- B.md "B2.0" promotes `notes/r2rl/merge_rounding.py`. That script computes `merged = (wb + d).to(bf16)` with d in
  fp32. `lora_merged_params` pushes `p.data + d.to(p.dtype)`, which rounds d to bf16 first [code: fast.py:717].
  The effect is about 2^-9 relative on d, so the 0.531 does not change. But the promoted script should mirror the
  code.
- Fix: in B2.0, add "the relative error uses `(W_bf16 + d.to(bf16))` in bf16, as `lora_merged_params` does". Keep
  the fp32-sum line as a second column.

### 7. MINOR (B): option F fallback restart needs a plan step that the doc does not name

- B.md "B2.2" F step 5: "On a failure it raises after writing the file, so the plan step fails and the run restarts
  under the B2.4 rule 0 fallback". B2.4: "No code edit inside a rental. A run picks its preset with `train_rl
  --numerics <name>`". Plan steps are `&&`-chained (B1.5 step 2), so a failing step stops the chain, and the doc
  names no fallback plan file.
- Fix: the plan of the first run after Gate B3 ships a second plan, `<run>_fallback.yaml`, with the rule-0 fallback
  preset. The main agent launches it after the `f_check.json` read. The cost line keeps "6 min + 15 steps".

## Counts

blocking 1, major 1, minor 5.
