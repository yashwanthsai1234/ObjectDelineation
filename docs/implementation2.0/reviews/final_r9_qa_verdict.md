# Verdict on final review round 9, QA (final_r9_qa.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice is in scope. It pointed me at the box build chain (one `_gpucheck` build per venv, `make cuda`), which decides qa-2. |
| engineering:code-review | I read `policy.py:375-380, 400-440, 475-492`, `train_sft.py:71, 82`, `gpucheck_py.cpp:113`, `CMakeLists.txt:15`, `Makefile:41-45`, `tests/test_cuda_checker_gpu.py`, `tests/test_gpu_scoring.py` against each claim. |
| engineering:debug | I traced the G08b `-m` step: cwd `/workspace/src`, only `OBJDEL_RUN_DIR` exported, bundle at `/workspace/data/lidc` (qa-1). I recomputed the R5 step-time SD (qa-3). |
| engineering:system-design | I matched each C.9.5 arm against the V2 item rows of front 6.3 (qa-4). |
| engineering:architecture | One owner per rule: G-SP1 vs E-SP1 (qa-3); A.9.6 skip rule vs D S13 fail rule (qa-6). |
| engineering:documentation | I checked whether front 6.1 step 4d's content is already owned by named tests (I20, I22a, I22b, U51); it is, so qa-5 is a pointer fix. |
| huggingface-skills:huggingface-best | No model choice is in scope. No change. |

Measurement made: R5 `step_time` over steps 3-15 (13 values each): base mean 30.78 s, SD 1.887 s; shared mean 23.43 s, SD 2.110 s [measured: `runs/r5/base_rl.json`, `runs/r5/shared_rl.json`, `statistics.stdev`, 2026-11-09].

## Verdicts (blocking and major items)

- qa-1 | confirmed (major, down from blocking) | B3.9 says the smoke reads `{data}/checker_eval/masks.npz`, but the `-m` step argv carries no data path and only `OBJDEL_RUN_DIR` is exported [code: policy.py:486, 492], and D8's only stated root rule (`OBJDEL_DATA_ROOT`, Mac default `data`) resolves to `/workspace/src/data/...`, which the code sync drops [code: policy.py:379]; not blocking because the repo already has a box fallback (`BOX_DATA = /workspace/data/lidc` when it is a dir [code: experiments/phase2/train_sft.py:71, 82]) that an engineer might reuse, so the doc forces a guess rather than a certain failure.
- qa-2 | confirmed (major) | B.9.5 G01 ("CPU build and CUDA build give bit-equal results") and the B.9.8 sketch `test_sums_cpu_gpu_bit_equal(cpu_build, cuda_build)` need two builds in one process, but `build()` returns one value [code: cuda/gpucheck_py.cpp:113] and `make cuda` reinstalls in place [code: Makefile:41-45]; A.9.6 G02b already states this and uses a Mac golden, so the G01/G08b reference is undefined and G01 is named in two files (B.9.5 runner row vs B.9.8 sketch).
- qa-3 | confirmed (blocking) | C.9.5 E-SP1 "std of steps 3-15 at most 0.3 s in each arm, else "noise, rerun"" fails on every real arm (measured SD 1.89 s and 2.11 s), so E-SP1 always orders an unpriced 43 min / $0.49 rerun while C.6 G-SP1 says "Anything above: not adopted" with no rerun: two rules for one decision, and following E-SP1 spends GPU money for nothing.
- qa-4 | confirmed (major) | C.9.5 asks for a second baseline arm (E-EQ6 control), a second fast arm (E-REP, mapped to G-EQ6 and G-SP1 in C.9.7), a forced-1,536 run (E-MEM) and a B = 128 arm (E-SP4), but front 6.3 V2 item 3 prices only "12 + 8 + 7 + 4 ... + 2 x 6 start-up" (one baseline, one fast+check arm) and lists no such optional rows, so the engineer must guess whether to run unpriced arms or drop tests that gates depend on.
- qa-5 | downgraded (minor) | The content of front 6.1 step 4d is already owned and specified: A I20 (`a3_measure.yaml` views, no `--split test`, build-if-missing first), B I22a (`b1_rows.yaml` step 0 build, views exactly val-shifted and val-edge, no `--split test`), B I22b (`b2_fp16.yaml`, `b3_v5.yaml` step 0 = `b1_rows.yaml` step 0), C `c_speed.yaml` step 0, E U51 (`e_v4.yaml` stage 0 = `b1_rows.yaml` step 0); since every plan builds itself, the launch-order check adds nothing, and the fix is to cite those tests instead of a "grep test in `tests/test_repo_hygiene.py`".
- qa-6 | confirmed (major) | D.9.6 S13 (`data` marker, "a missing file fails, never skips", reads `runs/<V3 run>/...`) contradicts A.9.6 ("A test whose `runs/...` input is absent skips with a reason naming the file") and A.6 "All" (`make test-data` passes); under default (b') V3 never runs, so `make test-data` stays red from D7 on, and `<V3 run>` has no rule to resolve it.

## Minor items (7-14)

Not re-verified; they stay minor as filed by the reviewer.

Totals: 5 confirmed (qa-3 blocking; qa-1, qa-2, qa-4, qa-6 major), 0 refuted, 1 downgraded to minor (qa-5).
