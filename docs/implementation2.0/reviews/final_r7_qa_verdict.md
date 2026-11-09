# Verdict on final_r7_qa.md (fresh verifier, 2026-11-08)

## Skills used
| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product claim at stake; pointed me to the box toolchain path (`run_box.sh --gpu-check`, `D_GPU_TESTS`) used to check qa-2. |
| engineering:code-review | Read `experiments/phase2/vast/policy.py` `cmd_run` (lines 465-479) to test whether a bare pytest can go through `run.sh` (qa-2). |
| engineering:debug | Traced the restart path: `--f-check` fails -> `train_rl` raises -> `restart_preset` decides; checked which inputs the U51/B.9.8 seam has (qa-1). |
| engineering:system-design | Compared B2.4 rule 0, B2.2 F step 5 / U76 (the `f_check.json` fields, `--f-check` default on under fp16-seqmask too) and U51/B.9.8 side by side (qa-1). |
| engineering:architecture | Checked the V3 launch lines in 00_front 6.3 against D8's `d_gpu_test.yaml` (one step only) for qa-2. |
| engineering:documentation | Checked "one rule, one place" for the restart rule and the G02b/G03 run place. |
| huggingface-skills:huggingface-best | No model choice involved; no change. |
| Web tools | Loaded (WebSearch, WebFetch); not needed: both items turn on repo text and code only. |

## Verdicts
- qa-1: CONFIRMED (major) - B "#### B2.4" rule 0 says a failed check restarts under fp16-seqmask or rules 2-3, and a run that started under fp16-seqmask restarts under rules 2-3; but B "### B.9.1" U51 defines "F passed (NaN 0, LoRA tensors equal ..., ratio <= 1.15)" with no comparator/load-count/base-hash term, and the B.9.8 sketch `F = lambda ok, ratio: dict(nan=..., lora_equal=ok, lora_over_base=ratio)` and `restart_preset(CF(), f_check=f, fp16=fp, b1=B1)` have no start preset and no `passed` flag, while U76 says `--f-check` is on by default under fp16-seqmask too; a seam built to U51 returns lora-seqmask after a comparator-only failure (paid restart into the failed preset).
- qa-2: CONFIRMED (major) - A "### A.9.6" G02b/G03 say "as its own command after D's `make gpu-check` ... the front prices it inside R-D2", but D "D8" `d_gpu_test.yaml` has one step (`run_box.sh --gpu-check`, `pytest $D_GPU_TESTS -m gpu`) and `D_GPU_TESTS` excludes `tests/test_a_gpu.py`; front 6.3 V3 item 1 and its launch line name only `d_gpu_test.yaml`; and `cmd_run` dies on anything but a module or yaml (`die("module must be experiments.phase2.<script> or experiments/phase2/configs/<name>.yaml")`, policy.py:479), so A-G2's G02b and the GPU-rule check G03 have no runnable step.
- qa-3: minor (kept minor) - GP5 `group_id` source; a reader can take it from A3 answers.jsonl; no outcome change proved.
- qa-4: minor (kept minor) - noise-floor print line is a nice-to-have; no gate change.
- qa-5: minor (kept minor) - seed/split pin in U51; E.8 step 3 already writes `--split val --seed 11`.
- qa-6: minor (kept minor) - naming only; the doc already says "at the 90 % level".
- qa-7: minor (kept minor) - runner-table file name for G01/G06.
- qa-8: minor (kept minor) - 9 vs 10 min wording; the front 6.3 sum is the priced one.
- qa-9: minor (kept minor) - wording "before training starts".

Count: 2 confirmed (2 major), 0 refuted, 0 downgraded (7 items were already minor).
