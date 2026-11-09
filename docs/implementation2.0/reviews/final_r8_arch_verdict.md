# Verdict on final review round 8, architecture (final_r8_arch.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes an NVIDIA product claim. No catalog lookup or install needed. |
| engineering:code-review | I opened every cited code line: `policy.py:422-495` (`chain`, `cmd_run`), `train_rl.py:20-52`, `g_probe.py:125-160`, `select_ckpt.py:365-392`, `train_sft.py:95-103, 255-261`, `final_eval.py:69`, `tests/test_p2_scripts.py:300-390`. |
| engineering:debug | I traced `run.sh <run> <yaml>` to `cmd_run`: it calls `chain(ROOT / a.module, ...)` on the Mac after launch, sync and fetch (item 1). I grepped every live `fork_pool`/`pool_size` caller (item 3). |
| engineering:system-design | I walked the V4 rental order to price a path guard in `chain()` (item 1). |
| engineering:architecture | I compared front 2.2/2.3 with E5 (`policy/callbacks.py`), E8 step 1 (no `infra/paths.py`) and E2.5 (`core/provenance.py`) (item 4). |
| engineering:documentation | I re-ran a check 14 style grep for scratch paths with no DOC1 step 0 row: 18 lines in 15 files with my narrower pattern. The defect exists, but every hit sits in a copied review record (item 6). |
| huggingface-skills:huggingface-best | No model choice. Nothing changed. |
| WebSearch/WebFetch | Loaded. Not needed: every claim is local text or code. |

## Verdicts

- arch-1: confirmed (blocking). E.md I29 says "`chain()` reads only `specs/*.yaml`", and E-G9 says "the box chain reads only `specs/*.yaml`". But E6 step 8b, E.8, U51 and front 6.3 put the live V4 `RunSpec` at `experiments/phase2/configs/e_v4.yaml`, and `cmd_run` calls `chain()` only after the paid launch, sync and fetch [policy.py:469-470]. No Mac test calls `chain()` on that file, so if you build the gate as written, V4 dies on a paid box.
- arch-2: confirmed (major). U52 asserts dests == stage fields minus kind/name/soft/dropped. But `TrainerFlags` gives `RlStage` the fields `compile` and `plain`, and `train_rl.py:23-36` has neither flag (no phase adds them). It gives `ProbeStage` an `lr` field that `g_probe.py:131-145` has no flag for. `SftStage` drops `stages`, `sampler` and `rows_plan`, which `TrlParser([RunArgs])` still exposes [train_sft.py:100-101, 260]. So E-G9's U52 fails on E's own models, and the engineer must pick an interface change.
- arch-3: confirmed (major). D7 step 6 removes the last live `fork_pool` callers (`final_eval.py:69`, `g_probe.py:157`, `select_ckpt.py:390`; CLEAN2 removes `train_rl.py:50`) before E0. The R-2 branch in E3/E.3.3 then archives `pool_size`. Yet `select_ckpt.py:374` and `g_probe.py:139` still call `default=serve.pool_size()`, E6.1 still types `procs`, and CLEAN16 step 1 claims the deletion. Three texts disagree on who deletes `--procs`, and the preflight dry run [test_p2_scripts.py:369, 384] breaks.
- arch-4: confirmed (major). Front 2.3 "Target tree" lists `rl/ ... callbacks.py`, `infra/ paths.py` and `core/ spec.py` only. E5 puts `TraceStep` in `policy/callbacks.py`, because `train_sft.py:390` and `g_probe.py:93` both use it. E8 step 1 says "There is no `infra/paths.py`". E2.5 creates `core/provenance.py`. A tree built from the front breaks the rl|sft import contract.
- arch-5: downgraded (minor). CLEAN8 step 1 already says to give the select stage "the same override and a test". The field or flag name is a design detail that does not change the outcome. The step's purpose (a token cap that reaches `select_ckpt`) shows that the Qwen case keeps the select call.
- arch-6: downgraded (minor). The finding is real, but every hit is a scratch path inside a copied review record (for example `reviews/final_r4_data.md:27`, `reviews/A_phase_review.md:10`). Any replacement makes check 14 green, so no outcome changes, and review-record text is never major.
- arch-7: minor (as filed).
- arch-8: minor (as filed).
- arch-9: minor (as filed).
- arch-10: minor (as filed).
- arch-11: minor (as filed).
- arch-12: minor (as filed).
