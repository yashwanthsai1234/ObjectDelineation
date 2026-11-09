# Phase 1 one-off scripts, archived 2026-11-09

**Nothing here is built, imported or tested.** Four Phase 1 scripts with no live caller.

## What it was

`census_sizerule.py`, `eval_demo.py`, `sample_eval.py`, `smoke_grpo.py`. Their tests moved to `archive/phase1_scripts/tests/test_scripts.py` and are not collected:

- `tests/test_experiment_scripts.py::test_sample_eval_uses_bf16_only_on_cuda_cards_that_support_it`
- `tests/test_experiment_scripts.py::test_smoke_grpo_group_stats_use_the_gate_definitions`
- `tests/test_eval.py::test_eval_demo_builds_targets_with_the_prepare_size_rule`

(The three functions no longer exist in those live files. The `census_sizerule` help check left the live parametrized lists.)

## What replaced it

`experiments/phase1/{build_data,build_targets,fetch_lidc_pixels,prepare_rl,eval_run,reward_gates}.py` and `src/objdel` cover what these did.

## Retiring commit

Code last present in commit `207685bf1b57` (tag `pre-clean-2026-11-08`). The commit that adds this README records the move.

## How to restore

`git restore --source=pre-clean-2026-11-08 -- experiments/phase1/census_sizerule.py experiments/phase1/eval_demo.py experiments/phase1/sample_eval.py experiments/phase1/smoke_grpo.py` (owner).
