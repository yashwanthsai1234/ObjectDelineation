# Verdict on final_r11_rl.md

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers TRL/vLLM sleep-flag handling. No verdict depends on it. |
| engineering:code-review | Re-read the TRL 1.17.1 code paths cited by item 1 in the repo .venv (vllm_generation.py:376-379, 500-503, 600-603, 718-719, 771-773; grpo_trainer.py:1203, 1838-1841) and fast.py:495-504. Every code quote in item 1 matches. |
| engineering:debug | Traced the flag through start-up, comparator (B2.2 F 4b.3/4b.5), and step 1 for each reading of "if sleeping". Found `_last_loaded_step = -1` (grpo_trainer.py:1203), so step 1 always calls the replaced sync; under a flag-based reading the flag is False after the comparator, so the wake is skipped. |
| engineering:system-design | Confirmed two owners of the "asleep" state (TRL's flag vs the doc's raw `sleep(level=2)` and "if sleeping"), with no rule linking them. |
| engineering:architecture | Checked C5 Interface, C5 F sync step 1, B2.2 F 2.1, 4b.3, 4b.5, U70, U76, U-SY3, U-SY8: none names the flag (grep `_llm_weights_sleeping`/`is_sleeping` in all phase sources: 0 hits). |
| engineering:documentation | Each verdict quotes source + heading and code line. |
| huggingface-skills:huggingface-best | No model-choice item in this review; not used for a verdict. |
| Web tools | Loaded (WebSearch, WebFetch). Not needed: all evidence is in the pinned TRL source in the repo .venv. |

## Verdicts

- rl-1 (MAJOR) — CONFIRMED, major. C.md "C5" Interface says "`push_weights` ...: if sleeping, `empty_cache()` and `wake_up(tags=["weights"])` (TRL's order, vllm_generation.py:500-503)" and B.md "B2.2" F step 2.1 / 4b.3 ("`wake_up(tags=["weights"])` if asleep") / 4b.5 ("in mode B only, `sleep(level=2)` last") never define the sleep state, while TRL keys both the re-sync and its own state on `self._llm_weights_sleeping` (vllm_generation.py:379, 503, 602-603, 772-773) and the trainer always syncs at step 1 (`_last_loaded_step = -1`, grpo_trainer.py:1203, 1838-1841); a flag-based "if sleeping" leaves the flag False after the comparator's raw sleep, so step 1 writes into slept weight memory, and a vLLM-state-based check never clears the flag, so every step syncs twice and the F check's "LoRA load count per RL step is 1" fails on the box; U70 counts `add_lora` "once per sync", U-SY3/U-SY8 test one call, U76's mode-B fake raises on reads only, so no planned $0 check covers it and the failure costs a rental start-up.
- rl-2 (minor) — filed as minor; agree (feeds an optional owner switch and a report line, not a gate; U73 defines the metric as written).
- rl-3 (minor) — filed as minor; agree (callback needs a trainer reference; found at $0 by U68's live-delta test).

Counts: confirmed 1 (major), refuted 0, downgraded 0, minor as filed 2.
