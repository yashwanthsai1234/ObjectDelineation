# Verdict on final_r6_rl.md (fresh verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill applies to TRL/vLLM LoRA sync; no verdict depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:488-504, 696-717`: `lora_merged_params` yields `p.data + d.to(p.dtype)` over the ACTIVE adapters (item 1, item 2a). |
| engineering:debug | Traced engine state through the comparator: merged push (i) -> `sync_lora` (ii) -> step 1 under C5 mode A (item 1). |
| engineering:system-design | Checked which engine state (base W, LoRA slot, KV awake) each comparator arm needs in mode A and B (item 2b). |
| engineering:architecture | Checked rule 0 timing vs U70 "every generate" (item 2c). |
| engineering:documentation | Each line names the proving doc text and code line. |
| huggingface-skills:huggingface-best | No model choice is under review here; no change. |

## Verdicts

- rl-1 (BLOCKING): CONFIRMED blocking. C.md "C5" interface says `sync_lora` "calls `push_weights(..., values="base")` in mode B only" and F step 2 says "Mode A: no base push" (U-SY8: "mode A has no `load_weights` after start-up"), while B.md "B2.2" F step 4(i) pushes `values="merged"` = `p.data + d.to(p.dtype)` [code: fast.py:717] into vLLM's base tensors; nothing restores W, so in mode A (planned: C.6 "Fast arm = C1-C5, mode A") arm (ii) and every step decode from W + d_B1 + slot; B.md F step 2.2 "no base push after the first sync" contradicts C5 (it would push once and hide the bug), so the two rule texts also disagree.
- rl-2 (MAJOR): CONFIRMED major. B.md "B2.2" F step 4 says only "loads the B1 final adapter ... on `rl-full-3ep`" while `sync_lora` writes `get_peft_model_state_dict(model)` whose default is `adapter_name="default"` [code: .venv peft/utils/save_and_load.py:148-151] and `lora_merged_params` merges active adapters [code: fast.py:703-707]; no step says how B1 coexists with the run's own adapter, how the KV cache is woken in mode B before step 1 (TRL sleeps at start-up and wakes `kv_cache` only in its own generate [code: vllm_generation.py:376-377, 718-719]), or how arm (i) and the rule-0 base timing run with `lora_request=None` when B2.2 and U70 put the request on "every `generate` call" -- an engineer must guess, and a wrong guess fails or falsifies the F check.
- rl-3 (minor): as filed, minor (rule 3/rule 4 start delta overlap; wording).
- rl-4 (minor): as filed, minor (callback construction detail).
- rl-5 (minor): as filed, minor (G-C beta condition wording).
- rl-6 (minor): as filed, minor (seed script rounding differs by ~2^-9 of d; 0.531 unchanged).
- rl-7 (minor): as filed, minor (fallback plan file name).

Counts: confirmed 2 (1 blocking, 1 major); refuted 0; downgraded 0; minor as filed 5.
