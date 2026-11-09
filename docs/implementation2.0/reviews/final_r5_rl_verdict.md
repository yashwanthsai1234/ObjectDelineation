# Verdict on final_r5_rl.md (round 5, RL review)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers vLLM LoRA loading or TRL eval paths; no verdict depends on it. |
| engineering:code-review | Read `fast.py:524-602` (C10 skip, `_compute_loss`, `set_is_ratio`, `log`) and TRL `grpo_trainer.py:2750-2782, 3416-3425`; confirmed rl-3. |
| engineering:debug | Traced vLLM 0.32.0 (local uv cache): `lora/utils.py:167-220` (key parse raises on `.default.`), `lora/peft_helper.py:176-181` (reads `adapter_config.json`), `lora/lora_model.py:239-240` (reads `adapter_model.safetensors`); header read of `checkpoint-100/adapter_model.safetensors`: 392 keys, 0 with `.default.`. Confirmed rl-1. |
| engineering:system-design | Checked which rental/process each gate runs in: fp16 arm V3 item 2b vs B1 run V1 (rl-2); G05 pytest vs the training process holding vLLM 0.35 x 32 GB + trainer ~17 GB (rl-5). |
| engineering:architecture | Owner decision 15 menu against C8 / B.8.3 item 5 / U-CO5 (rl-4). |
| engineering:documentation | Each verdict names source file + heading and the proving line. |
| huggingface-skills:huggingface-best | Not used for any verdict (model list not in dispute). |
| Web | docs.vllm.ai LoRA page: `adapter_config.json` named in the resolver section, `load_inplace` semantics; the file names were then confirmed in the 0.32.0 source. |

## Verdicts

- rl-1: CONFIRMED major. B.md "B2.2 Numerics presets" F step 2.3 ("`safetensors` from `{k: v.to(torch.bfloat16)}` of a state-dict copy") and C.md "C5" F sync step 3 ("Write the adapter as bf16 safetensors to a tmpfs folder") name no file name, key format or config file, while vLLM 0.32.0 reads `os.path.join(lora_dir, "adapter_model.safetensors")` (`lora_model.py:239`) and `adapter_config.json` (`peft_helper.py:176-181`) and raises `"{name} is unsupported LoRA weight"` unless `parts[-2] in ["lora_A","lora_B"]` (`lora/utils.py:207-220`); a PEFT `state_dict()` key carries `.default.`, so the literal reading fails at the first `add_lora` and U70 (fakes) cannot catch it.
- rl-2: CONFIRMED major. B.md "B2.4 Decide" rule 1 gates on "its s/step is at most 1.05 x the B1 run's", while rule 0 in the same section says "A cross-run s/step ratio is a report line only" (host moves the CPU phase, different answers move the decode tail); the fp16 arm is V3 item 2b and the B1 run is V1 (00_front 6.3 table), two rentals, so the 1.05 bar is an untestable gate that can flip the rule-0 fallback on host noise (the reviewer's "16 % host spread" compares different code, 27.5 s non-shared vs 23.65 s shared-prefix, but rule 0's own text proves the point).
- rl-3: CONFIRMED major. B.md "B2.1" says `set_is_ratio` "appends the micro-batch's per-sequence `seq_kl` values" and `log` writes "one 256 x 3 block per logged step" with no train-mode condition, and B1.5 offers the held-out eval (owner decision B.8.3 item 15); in eval TRL's `prediction_step` calls `compute_loss` (`grpo_trainer.py:3416-3420`), the C10 skip has no training check (`fast.py:524-538`, `skip_old_logps`), and `_compute_loss` then routes every eval micro-batch to `set_is_ratio` (`fast.py:567-577`), so with the held-out curve on, 512 val-shape rows enter `seq_gaps.f32` at steps 50, 100, ... and shift the rules 3/4 windows.
- rl-4: CONFIRMED major. 99_back.md "RL settings" decision 15 offers "Alternative: beta 0.0 from the B1 run", but C.md "C8" says "`b1_rows.yaml` keeps 1e-3 in every case", U-CO5 asserts it, and B.md "B.8.3" item 5 says "B1 beta 1e-3 (fixed, as R5)"; an owner who picks the offered alternative breaks Gate B1's data-only claim and fails `make test` (U-CO5), blocking the `v1-b1` tag.
- rl-5: CONFIRMED major. B.md "B.9.5 GPU tests" command table runs G05 as `python -I -m pytest -m gpu tests/test_train_configs.py -k G05` "inside its first 15 steps" of a training run, but a pytest is its own process and cannot start a second vLLM engine beside the run (`vllm_gpu_memory_utilization: 0.35` at `configs.py:192` plus a ~17 GB trainer peak, C.md C5), G05's "vLLM LoRA load count per RL step is 1" needs the run's engine, B2.4 rule 0 says the F check is "in one process on one host" with "0 extra minutes" (B.8.1), and B.md B1.5 step 2 forbids `-k` ("an empty selection exits 5 and stops the `&&` chain", `policy.py:490`); the engineer must guess where the check lives.
- rl-6: minor (as filed).
- rl-7: minor (as filed).
- rl-8: minor (as filed).
- rl-9: minor (as filed; front wins by its own rule).
- rl-10: minor (as filed; block is headed "code changes, not now").
- rl-11: minor (as filed).
- rl-12: minor (as filed).

Counts: confirmed 5 (blocking 0, major 5), refuted 0, downgraded 0.
