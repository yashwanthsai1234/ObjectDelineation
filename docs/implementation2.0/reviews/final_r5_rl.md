# Final review round 5: RL (GRPO, IS correction, reward, TRL 1.17.1 / vLLM 0.32.0 internals)

Scope: B2, B3, C (C5-C11 in depth) and the model choice (front section 3); blocking items anywhere.
The doc was read by sections from `docs/implementation2.0/Implementation.md`. Every claim below was checked against the
sources in `docs/implementation2.0/`, the repo code, the TRL 1.17.1 install in `.venv` and the vLLM 0.32.0 source in the
local uv cache (`~/.cache/uv/archive-v0/LSXZLRiG9Ee6shbS/vllm`, the `+cpu` wheel; the Python is the same as the box
wheel). Items already in `review-log.md` (rounds 1-4) are not repeated.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers TRL/vLLM LoRA serving or GRPO numerics; no finding depends on it. |
| engineering:code-review | Read the code each doc claim names: `fast.py:495-617` (C10 path, `set_is_ratio`, `log`), `configs.py:177-196`, TRL `grpo_trainer.py` (2600-2611, 2760-2824, 2880-2905, 3066-3089, 3180-3200, 3270-3290, 3416-3422), `vllm_generation.py:356-372, 495-528, 595-610, 700-775`. Items 1, 3, 6 come from it. |
| engineering:debug | Traced the vLLM F path in the 0.32.0 source: `lora/worker_manager.py:288-321` (remove then re-add on `load_inplace`), `lora/layers/base.py:36-46` (wrapper passes `load_weights` to the base layer, so the F base push works: no finding), `lora/utils.py:167-220` and `lora/peft_helper.py:176-181` (item 1). |
| engineering:system-design | Checked which process and which rental each gate runs in (items 2, 5, 8). |
| engineering:architecture | Checked the owner-decision texts against the rules they feed (item 4, item 9). |
| engineering:documentation | Each item quotes the source file and heading and the proving line. |
| huggingface-skills:huggingface-best | HF API re-check of the model table: Qwen3-VL-4B 4,437,815,808, Qwen3-VL-2B 2,127,532,032, Qwen3.5-4B 4,659,865,088, MedMO-4B 4,437,815,808 parameters, all Apache-2.0; config.json: vision depth 24 / width 1024 / DeepStack [5, 11, 17] in both, merger out 2,048 vs 2,560, text 28 vs 36 layers. The front section 3 numbers match. |
| huggingface-skills:trl-training | Confirmed the TRL 1.17.1 semantics the doc relies on: `sequence_mask` = exp(sum of token gaps), masked only above `clip_max` when `clip_min` is None; OPSM keeps `adv >= 0 or seq_kl <= delta` with `seq_kl` = vLLM minus trainer; `sum_then_normalize` + `scale_rewards="none"` = R minus group mean only; `dr_grpo` divides by `B x max_completion_length`; `prediction_step` calls `compute_loss` (item 3). |
| huggingface-skills:huggingface-papers | Confirmed arXiv 2601.13633 is "EGM: Efficient Visual Grounding Language Models" (the cited source of the 2B-vs-4B gap); OPSM is DeepSeek-V3.2 (2512.02556), as the TRL docstring says. |

Web: one search on vLLM LoRA `base_layer` loading ([docs.vllm.ai LoRA base layer](https://docs.vllm.ai/en/v0.29.0/api/vllm/lora/layers/base/)); it pointed to the `base.py` passthrough, which I then read in the 0.32.0 source.

Verified and correct (no item): B2.0 planted 2x2 case (share 0.5, ulp(1) bf16 = 2^-7); B2.1 OPSM sign and formula;
B3.4 RLOO identity, the 8/7 factor, the "unchanged by `sum_then_normalize`" claim, the Adam (g/r, eps/r) invariance
and the 1.15e-6 RMS; B3.4 worked check (+0.0567); C5 1.41 B touched parameters; C6 `max_num_seqs` 256; C7 1,536 / L
scaling; C8 beta-0 path; F base push through LoRA-wrapped modules; F re-add after a level-2 wake (remove, add,
activate); `LLM` patch point at module level (`vllm_generation.py:42-43`).

Counts: blocking 0, major 5, minor 7.

## Items

### 1. major. F adapter file: format not named; the natural reading makes `add_lora` raise on the box (B, C)

- Doc, B.md "B2.2 Numerics presets", Option F plumbing step 2.3 (B.md:1133-1134): "Write the adapter, cast to bf16,
  to `<tmpfs dir>` (`safetensors` from `{k: v.to(torch.bfloat16)}` of a state-dict copy; 17,432,576 parameters".
  C.md "C5. Weight push", F sync step 3 (C.md:759): "Write the adapter as bf16 safetensors to a tmpfs folder".
  U70 / U-SY8 check only "the adapter file holds bf16 adapter tensors only". Nothing names the key format or
  `adapter_config.json`.
- Proof. A PEFT model's state dict names its tensors `...q_proj.lora_A.default.weight` (adapter name inside the key).
  vLLM 0.32.0 accepts only `...lora_A.weight` / `...lora_B.weight`: `parse_fine_tuned_lora_name` checks
  `parts[-2] in ["lora_A", "lora_B"] and parts[-1] == "weight"` and otherwise ends in
  `raise ValueError(f"{name} is unsupported LoRA weight")` (`vllm/lora/utils.py:207-220`). vLLM also reads
  `os.path.join(lora_path, "adapter_config.json")` (`vllm/lora/peft_helper.py:176-181`), which a bare safetensors
  write does not create. The saved checkpoints show the format vLLM needs: `checkpoint-100/adapter_model.safetensors`
  has 392 keys of the form `base_model.model.model.language_model.layers.0.mlp.down_proj.lora_A.weight` (no
  `.default.`), F32 [measured: header read, 2026-11-08]; PEFT strips the adapter name in
  `get_peft_model_state_dict` (`peft/utils/save_and_load.py:148`).
- Effect. An engineer who copies "a state-dict copy" gets `ValueError` at the first `add_lora` of the first run after
  Gate B3. B2.4 rule 0 then says "If the F check fails, the run restarts under `fp16-seqmask` ... or
  `bf16-tis-opsm`": F is rejected for a file-format bug, and the restart costs about 6 min plus 15 steps.
- Fix (B2.2 F step 2.3, C5 F sync step 3, ADR-0013 F text): "Write `peft.get_peft_model_state_dict(model)` (keys
  without the adapter name), each tensor cast to bf16, as `adapter_model.safetensors`, and
  `model.peft_config['default'].save_pretrained(<tmpfs dir>)` (writes `adapter_config.json`), before `add_lora`."
  Add to U70 and U-SY8: every key matches `^base_model\.model\..*\.lora_[AB]\.weight$` (no `.default.`), the key
  count is 392, and `adapter_config.json` exists with `r` 16 and `lora_alpha` 32. Add to G05: `add_lora` returns
  True.

### 2. major. B2.4 rule 1 gates fp16 on a cross-rental s/step ratio that rule 0 calls report-only (B)

- Doc, B.md "B2.4 Decide", rule 1 (B.md:1241-1243): "`fp16-seqmask` passes if ... AND its s/step is at most 1.05 x
  the B1 run's". The same section, rule 0 (B.md:1237-1239): "A cross-run s/step ratio is a report line only: the
  host moves the CPU phase (3.0-3.2 s of the R5 step, BRIEF section 2), and two runs sample different answers, so
  the decode tail differs too."
- Proof. The fp16 arm is V3 item 2b; the B1 run is V1 (front 6.3 and 6.2): two rentals, so two hosts in general.
  The front's own planning inputs carry a 16 % host and data spread for the same code (step time 27.5 s vs 23.65 s,
  front section 6 table), and the decode-tail fit is 4.3 ms per token of the longest answer (C.1 table), so two runs
  with different longest answers differ by seconds. A 1.05 band is about 1.2-1.4 s at 23.65-27.5 s. The doc fixed
  the same cross-host defect for D-G7 (review-log Phase D, revision 2 item 4 "D-G7 compared two hosts").
- Effect. Rule 1 can drop or pass fp16 on host noise, and rule 1 feeds rule 0's fallback ("restarts under
  `fp16-seqmask` (if rule 1 passed)").
- Fix (B2.4 rule 1 and Gate B2): delete "AND its s/step is at most 1.05 x the B1 run's"; add "s/step of both runs is
  a report line". If a speed bar is wanted, it must be in-process like rule 0 (not possible across dtypes in one
  engine), so the honest rule is report-only.

### 3. major. Held-out eval rows enter the B2.1 buffers and break `seq_gaps.f32` (B)

- Doc, B.md "B2.1 Log what is lost", (B.md:1034-1052): `set_is_ratio` "appends the micro-batch's per-sequence
  `seq_kl` values and advantages to a device list"; `log` writes them to `<out>/seq_gaps.f32` as "one 256 x 3 block
  per logged step, in step order" with `{"rows_per_step": 256, ...}`; rules 3 and 4 price their thresholds from that
  file. B.md "B1.5" (B.md:818-829) offers the held-out curve (`eval_strategy="steps", eval_steps=50,
  num_generations_eval=8`, 64 items x 8 = 512 answers) and requires "the tiny TRL loop test ... runs one held-out
  eval with FastGRPO's patches on" before the `v1-b1` tag.
- Proof. TRL's eval calls the loss: `prediction_step` -> `self.compute_loss(model, inputs)`
  (`grpo_trainer.py:3416-3420`). FastGRPO's `_compute_loss` sets `self.loss_inputs` whenever
  `"importance_sampling_ratio" not in inputs and self.use_vllm and self.vllm_importance_sampling_correction`, with no
  training-mode check (`fast.py:567-571`); C10's skip also runs in eval (`fast.py:524-538`), so TRL leaves the ratio
  out and `set_is_ratio` runs on every eval micro-batch (`fast.py:573-577`). `Trainer.evaluate` then calls `log`
  with `model.training` False (`grpo_trainer.py:3424-3425`).
- Effect. With the held-out option on, one eval adds a 512-row block of val-shape rows to the file (and to the
  callback lists if eval falls in steps 1-20). A reader that reshapes by `rows_per_step` 256 shifts every later
  step by two, so the steps 50-260 windows of rules 3 and 4 read the wrong steps and include val rows; the doc's own
  rule is that val-shape only informs, never trains or prices the run.
- Fix (B2.1): "`set_is_ratio` appends to `seq_gaps` and the rule-4 callback lists only when `self.model.training`;
  `log` writes `seq_gaps.f32` only in train mode." Add to U68: one eval between two train steps leaves the file at
  2 x 256 rows and the quantiles unchanged.

### 4. major. Owner decision 15 offers "beta 0.0 from the B1 run", which C8, B.8.3 and U-CO5 forbid (99_back, C, B)

- Doc, 99_back.md "RL settings", decision 15 (99_back.md:147-148): "Alternative: beta 0.0 from the B1 run (RL review
  item 4)." C.md "C8. KL weight 0" (C.md:935): "`b1_rows.yaml` keeps 1e-3 in every case: changing it would break
  Gate B1's 'data only' claim"; U-CO5 "asserts that `b1_rows.yaml` has beta 1e-3". B.md "B.8.3" item 5: "B1 beta 1e-3
  (fixed, as R5)".
- Proof. The owner approves decision 15 from this menu before V1 (6.1 step 0). If the owner picks the alternative
  as written, the B1.5 plan must change, U-CO5 fails, `make test` fails, and the `v1-b1` tag cannot be made (order
  rule 2: the tests pass on that tree). The round-3 C fix (review-log C, item 11) changed C8 but left this
  alternative text.
- Fix (decision 15): "Alternative: beta 0.0 for B3 only (B1 keeps 1e-3 in every case, C8). It saves the reference
  pass of the B3 run ..."

### 5. major. G05 is a separate pytest process "inside" the run's first 15 steps (B)

- Doc, B.md "B.9.5 GPU tests", command table (B.md:2038): "| G05 | the first run that uses `fp16-seqmask` or
  `lora-seqmask` (after Gate B3), inside its first 15 steps | `python -I -m pytest -m gpu tests/test_train_configs.py
  -k G05` |". B2.4 rule 0 (B.md:1232-1233): "The F check is the first 15 steps of that run itself, in one process on
  one host". B.8.1: "(0 extra minutes)".
- Proof. A pytest command is its own process with its own vLLM engine. The training process already holds vLLM
  plus the trainer, about 28 of 32 GB (C5 step 2), so a second engine cannot start beside it; run before the training
  process, it costs one more start-up (about 6 min, front section 6 inputs), not 0 minutes, and it checks another
  engine than the one the F check is about. G05's own content (packed `lora_a` / `lora_b` equality after the first
  sync and after a sleep/wake cycle, "the vLLM LoRA load count per RL step is 1") needs the training run's engine.
  Also `-k G05` contradicts B1.5's rule "Node ids, not `-k`".
- Fix (B.9.5 G05 and B2.4 rule 0): "The F check runs inside `train_rl` (`--f-check`, on by default under
  `lora-seqmask` and `fp16-seqmask`): at sync 1 and after the first sleep/wake cycle it runs the G05 / G-EQ5 tensor
  checks, the comparator and the ABAB timing, writes `<out>/f_check.json`, and raises on a failure. G05 on the box
  is a node-id pytest that reads `f_check.json` after step 15; the fake-engine part stays U70."

### 6. minor. C5 `push_weights` must clear TRL's sleeping flag, or TRL syncs twice per step (C)

- C.md "C5" interface (C.md:777): "`push_weights` ...: if sleeping, `empty_cache()` and `wake_up(tags=["weights"])`
  (TRL's order, vllm_generation.py:500-503)". `install` replaces `vg.sync_weights` on the instance.
- TRL's generate re-syncs when `self._llm_weights_sleeping` is True (`vllm_generation.py:602-603`), and only TRL's
  own `sync_weights` sets it False (`vllm_generation.py:503`). The cited range includes that line, so a careful
  engineer copies it; the interface text and U-SY tests do not say it. Missing it doubles the push under mode B (and
  the base push, adapter write and `add_lora` under F), which inflates G-SP1's fast arm.
- Fix: add "then sets `vg._llm_weights_sleeping = False`" to `push_weights` and `sync_lora`; add the name to the A7
  private-name list; U-SY test: one TRL `generate` after a sleep calls the replaced sync exactly once.

### 7. minor. `completions/truncated_grad_share` ignores the IS weight, which drops long rows first (B)

- B.md "B3.5" (B.md:1379): "sum over truncated rows of |A| x tokens / sum over all rows of |A| x tokens".
- Under `bf16-seqmask` a row's weight is exp(sum of token gaps), masked above 3 (`fast.py:355-362`); the sum grows
  with length, and B.1 item 5 says "long answers pass the cap first". 1,536-token truncated rows are therefore masked
  more often than 131-token rows, and the metric overstates their gradient share; it drives the C11 masking option
  (median > 10 %).
- Fix: weight each row by its IS weight (sequence ratio after masking, or the mean token weight under token modes):
  sum over truncated rows of w_i |A_i| T_i / sum over all rows of w_i |A_i| T_i; U73 adds a masked truncated row.

### 8. minor. "The first run after Gate B3" has no start checkpoint or LoRA state (B, C, 00_front)

- It carries the F check, C-G-C, G-RC2, C7's L and rule 4's OPSM delta. B2.2 F step 4 loads "the B1 final adapter
  ... on `rl-full-3ep`" in that run's process, which needs that run's base to be `rl-full-3ep`
  (`vllm_generation.py:356`, `LLM(model=model.name_or_path, ...)`). Rule 4's reason ("a fresh adapter has the
  smallest gap in steps 1-20") holds only if that run starts a fresh LoRA. B.8.1 says "The front prices it when that
  rental is planned", so it is out of the current spend; one line avoids a later guess.
- Fix (front 6.3 or B2.4 "Scope"): "The runs after Gate B3 start, like B1 and B3, from `rl-full-3ep` with a fresh
  LoRA, unless the owner names another start; then the F comparator loads the adapter on that run's base."

### 9. minor. B.8.3 item 11 says $2.20 means 260 steps; the front says 220 (B)

- B.md "B.8.3" item 11 (B.md:1777-1778): "Keeping the older $2.20 means a 260-step B1 run". Front 6.2 notes: "If the
  owner keeps $2.20 anyway: ... the run plans 220 steps". The front wins by its own rule, so not blocking.
- Fix: "Keeping the older $2.20 means a 220-step B1 run (front 6.2 notes)."

### 10. minor. 4B recipe order: step 0 runs `train_sft` before step 1 adds the model name (00_front)

- Front section 3, "If the gate fires", step 0 (00_front.md:443): "Change exactly five lines: `model:
  qwen3-vl-4b`..." and runs `train_sft` on the box; step 1 then adds `"qwen3-vl-4b"` to `MODELS`.
  `train_sft.py:188` and `:308` read `configs.MODELS[args.model]`, so in the written order step 0 raises `KeyError`.
  The block is headed "code changes, not now", so a careful reader lands steps 1-2 first.
- Fix: put step 1 (model map) and step 2 (per-model `vllm_gpu_memory_utilization`) before step 0, or add to step 0
  "after steps 1-2 are merged"; add a Mac test that `sft_qwen3vl4b.yaml`'s `model` is a key of `MODELS`.

### 11. minor. Rule 3 "starting delta" against rule 4's step-20 freeze (B)

- B2.4 rule 3 raises "the starting delta to the q99.9 of `seq_kl`" over B1 steps 50-260; rule 4 says "Until step 20,
  delta = 0.5 (TRL docstring)" and then freezes q99.5 of steps 1-20. So rule 3 changes only steps 1-20; the text
  does not say so.
- Fix (rule 3): "raise the steps 1-20 delta (rule 4's starting value) from 0.5 to that q99.9; rule 4 still sets the
  frozen delta at step 20."

### 12. minor. The rule-4 callback needs a trainer reference (B)

- B2.4 rule 4: "a `TrainerCallback.on_step_end` at step 20 sets `trainer.off_policy_mask_threshold = delta`". HF
  callbacks receive `args, state, control, model, ...`, not the trainer.
- Fix: "the callback is built with the trainer (`OpsmDelta(trainer)`) and added with `trainer.add_callback`"; U-test
  with a fake trainer.
