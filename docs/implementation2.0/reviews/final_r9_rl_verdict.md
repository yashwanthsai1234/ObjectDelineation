# Verdict on final_r9_rl.md (fresh skeptical verifier)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog lane covers TRL/vLLM LoRA serving or the comparator prompt path; no verdict depends on it. |
| engineering:code-review | Read `src/objdel/train/fast.py:351-362` (`is_ratio`), `:407-415` (`one_pad_per_image`), `:618-626` (`_generate` installs the wrapper), `:670-690` (`patch_vllm_prompts`); `experiments/phase2/train_rl.py:58-59` (merge before save); `final_eval.py:48-58` (only `--ckpt`, no adapter flag). |
| engineering:debug | Traced the comparator before step 1: the wrapper is installed only inside `_generate`, so a direct `llm.generate` is vLLM's own; TRL's own path sends processor-expanded ids (`.venv trl/generation/vllm_generation.py:722-743`). Confirms item 1's failure path. |
| engineering:system-design | Traced the F-trained adapter's life: unmerged in training, merged in bf16 at `final/`, evaluated merged through `final_eval --ckpt`. Confirms item 2. |
| engineering:architecture | Checked whether B2.3b, C5, G-EQ5, U70 or U76 define the comparator's prompts or their pad form: none do (B2.3b's batch is saved answers for scoring, not the comparator's generate input). |
| engineering:documentation | Each confirmed item names its source file and heading and the code line that proves it. |
| huggingface-skills:huggingface-best | No model-choice claim in the major items; no verdict depends on it. |

## Verdicts

- rl-1: CONFIRMED (major). B.md "B2.2 Numerics presets" F step 4 uses "the 32 x 8 batch" (line 1245) without defining it, and step 4b.3 says the comparator "calls `llm.generate` directly, not through the `patch_vllm_prompts` wrapper"; the wrapper is installed only in `_generate` (fast.py:618-622) and is the only code that cuts pads with `one_pad_per_image` (fast.py:680-688), while TRL's prompts are processor-expanded ids (vllm_generation.py:722-743), so the obvious input double-expands the image and the F decision rests on a corrupted comparator.
- rl-2: CONFIRMED (major). B.md "B1.5 Train on the ordered rows" names `final/` (merged) as what `final_eval --ckpt` reads, and no line in B, C, front or back gives an eval path for a `lora-seqmask` run; `train_rl.py:58-59` runs `merge_and_unload()` before `save_pretrained(final)` and `final_eval.py:48` has only `--ckpt` (0 `enable_lora` hits in src/experiments), so an F run is scored on round(W+d), the very policy (rel. error 0.531 at step 100, B2.0) F exists to avoid.
- rl-3: confirmed as MINOR (no measurement of truncated rows masked by the sequence IS mask; reviewer's own grading, not blocking or major).
- rl-4: confirmed as MINOR (KV-budget estimate shows preemption, but the timing is a ratio of two arms under the same preemption; the doc's "batch 256" wording is the defect).
- rl-5: confirmed as MINOR (the 1.15 bar may fail on speed alone, but no benchmark proves it; estimate only).
- rl-6: confirmed as MINOR (fp32 `torch.exp` underflow to 0 is real at fast.py:358; wording fix in B.1 item 5).
- rl-7: confirmed as MINOR (G-C estimate omits F cost; no action follows from a G-C fail, so no outcome changes).
- rl-8: confirmed as MINOR (rules 3 and 4 overlap only for steps 1-20; small, bounded ambiguity).

Counts: confirmed major 2, refuted 0, downgraded 0 (items 3-8 were already minor).
