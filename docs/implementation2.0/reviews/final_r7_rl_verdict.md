# Verdict on final_r7_rl.md (fresh verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL/vLLM sleep mode or GRPO gates; no verdict depends on it. |
| engineering:code-review | Read TRL 1.17.1 `vllm_generation.py` in the repo .venv (lines 376-379, 495-528, 596-606, 712-722, 766-776) to check the sleep/wake claims of rl-1. |
| engineering:debug | Traced B2.2 F step 4b.5 call by call in mode B: the hash runs after `sleep(level=2)`, so it reads unmapped memory (rl-1). |
| engineering:system-design | Checked what each comparator number measures: the 0.5 bar ignores the d = 0 floor (rl-2); the timing ratio compares two answer sets of different length (rl-3). |
| engineering:architecture | Checked that B2.2 4b.5, C5 step 5(c) and G-EQ5 all place the base hash "after the restore", so all three need the order fix (rl-1). |
| engineering:documentation | One line per item, each with the doc line and the code line or measurement. |
| huggingface-skills:huggingface-best | No model choice in scope; not used for any verdict. |
| Web tools | Loaded (WebSearch, WebFetch). Not needed: the TRL source comments in the .venv state the level-2 sleep behaviour directly. |

## Proof read for the verdicts

- rl-1: B.md "B2.2 Numerics presets", F step 4b.5: "... `reset_prefix_cache()`, then in mode B `sleep(level=2)`. Then
  hash every vLLM base tensor on the host against `p.data`". TRL [code: .venv trl/generation/vllm_generation.py:378]
  "Sleep level 2 discards the weights"; [:497-499] "Without this, load_weights() writes to freed/unmapped memory when
  sleep mode is active". So a host read of the base tensors after `sleep(level=2)` reads unmapped memory. C.md C5 step
  5(c) ("after the comparator's base restore, every vLLM base tensor hashes equal to `p.data`") and G-EQ5 inherit the
  order. Simplest fix: hash before `sleep(level=2)` (mode B), or the reviewer's variant (mode B: no restore push, the
  hash runs after step 1's own sync).
- rl-2: B.md "B2.2", F step 4 (iii): "Pass: mean |f_lora - b| <= 0.5 x mean |f_merged - b|". Measured (re-run,
  `runs/*/rl.json`, `sampling/sampling_logp_difference/mean`): step 1 (fresh LoRA, delta 0) 0.009563 (r5 base),
  0.009515 (r5 shared), 0.009556 (full run); full-run medians 240-280: 0.01764, 101-300: 0.02361, 301-500: 0.01853.
  A perfect F (at the floor) gives 0.54 / 0.41 / 0.52: the gate outcome depends on the drift range, not on F. Not
  blocking (101-300 would pass a perfect F), but the gate cannot separate a working F from a failed one: major.
- rl-3: B.md "B2.4 Decide" rule 0: "generate the same 32 x 8 batch twice, with the comparator's `lora_request` and
  with `lora_request=None` ... pass if the median LoRA / base generate-time ratio is <= 1.15"; B2.2 4b.4 claims
  kernel cost only. Arm A samples W + d_B1, arm B samples W, so answers and their maximum length differ. Measured
  slope: generate = 5.61 s + 4.3 ms x longest answer, r = 0.92 [notes/B_gpu_audit.md:40, 84]. A 736-token gap in the
  longest answer moves the ratio by about 0.29 on an 11 s generate, twice the 0.15 margin. Seeded ABAB repeats the
  same answers, so the median does not cancel it. Major (the gate outcome depends on a variable the doc does not fix).

## Verdicts

- rl-1, confirmed (blocking): B.md F step 4b.5 hashes the vLLM base after `sleep(level=2)`, which TRL's own code says discards and unmaps the weights, so every mode-B F check fails or crashes before step 1 and triggers a paid restart.
- rl-2, confirmed (severity lowered to major, not blocking): the 0.5 bar ignores the measured d = 0 floor (0.0096 vs 0.0176-0.0236 at comparable drift), so a perfect F passes or fails by the drift range, an untestable gate; the floor-subtracted rule with f0 fixes it.
- rl-3, confirmed (major): the ABAB timing compares two policies whose answers differ in length, and decode time follows the longest answer at 4.3 ms per token, so the 1.15 ratio measures length, not LoRA kernel cost; fixed-length timing (`max_tokens=512, ignore_eos=True`) fixes it.
- rl-4, minor (as filed): scoring placement of (iii) and the "patches off" switch can be inferred by a careful engineer; wording only.
- rl-5, minor (as filed): the reviewer agrees the cited TRL lines include the flag reset; a test gap, not a wrong result.
- rl-6, minor (as filed): the truncated-share metric only feeds an optional owner trigger; a metric refinement.

Counts: confirmed 3 (1 blocking, 2 major), refuted 0, downgraded to minor 0 (rl-2 lowered from blocking to major).
