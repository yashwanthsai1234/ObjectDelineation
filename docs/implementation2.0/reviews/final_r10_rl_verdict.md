# Verdict on final review round 10: RL (final_r10_rl.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers TRL/vLLM sleep mode or LoRA sync. No verdict depends on it. |
| engineering:code-review | Read the TRL 1.17.1 start-up and generate code in the repo .venv (vllm_generation.py:356-379, 492-528, 599-603, 718-719, 771-773). It confirms the start-up `sleep(level=2)`. |
| engineering:debug | Traced the mode-B engine state through the comparator, arm by arm. Arm (ii) calls `sync_lora`, and `sync_lora` pushes the base in mode B. This refutes item 1. |
| engineering:system-design | Checked the one sync home (C5 `sync_lora`, steps 1-4) against the comparator calls in B2.2 F step 4b.2. |
| engineering:architecture | Checked that B2.2 F step 2, F step 4 (ii), 4b.2, C5 "Interface" and U70 agree that every F sync pushes `p.data` in mode B. |
| engineering:documentation | Each verdict names the source file, heading and line that decides it. |
| huggingface-skills:huggingface-best | No model-choice item in this review. Not used for a verdict. |
| Web (WebSearch/WebFetch loaded) | The review's vLLM claim ("level 2 discards the weights") is not disputed. TRL code says the same (vllm_generation.py:378). No fetch was needed. |

## Verdicts

- rl-1 (blocking): REFUTED. Arm (ii) "runs `sync_lora` with the same adapter" (B.md "B2.2", F step 4 (ii)), and the comparator calls it as `sync_lora(adapter_name="b1", slot=2)` (B.md F step 4b.2). `sync_lora` is "the F sync, steps 1-4 above. It calls `push_weights(..., values="base")` in mode B only" (C.md "C5", Interface). F sync steps 2.1-2.2 wake the weights, then push the unmerged base `p.data` in mode B (B.md F step 2). U70 asserts this push for every mode-B sync. So in mode B, arm (ii)'s own sync pushes base W after the start-up level-2 sleep (TRL vllm_generation.py:376-379). Arm (0) and the timing arms run after it on the same awake engine, through direct `llm.generate`, so nothing sleeps the engine between them. The garbage-weights scenario does not happen. Only the phrase "on the start-up base W" in 4b.4 is loose wording in mode B, which is minor.
- rl-2 (minor): not re-graded. The reviewer graded it minor and gave no proof of a wrong result (`final_eval --adapter` sync not fixed by the doc).
- rl-3 (minor): not re-graded. The comparator means stay unbiased, so only the noise grows.
- rl-4 (minor): not re-graded. A dead margin in G-SP1 that the same section says is never applied.
- rl-5 (carried round-9 minors): not re-graded. These were already confirmed minor in final_r9_rl_verdict.md.

Counts: confirmed blocking/major 0; refuted 1 (rl-1); downgraded 0. Items 2-5 were minor as filed.
