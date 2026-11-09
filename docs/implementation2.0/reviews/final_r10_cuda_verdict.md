# Verdict: final review round 10, CUDA (final_r10_cuda.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Review only, no catalog install; no NVIDIA-tool claim needed a new check for the blocking/major items. |
| engineering:code-review | Read run_box.sh:22-24, Dockerfile.vllm_native:20, Dockerfile.phase2:37 and policy.py:464-480 against the doc lines. |
| engineering:debug | Traced the post-E6 `run.sh <run> configs/...yaml` path (cmd_run -> chain locally after launch) and the post-E8 dry-run import. |
| engineering:system-design | Walked the V3 Mac-side pre-launch steps (front 6.1, D8 Steps, D9 runbook) to test item 4. |
| engineering:architecture | Checked E6 step 2, R-4 and E8's name table against the plan files that front 6.3 names. |
| engineering:documentation | One line per item, each with the source heading and the code line. |
| huggingface-skills:huggingface-best | No model choice in scope; nothing changed. |

Web: not needed; every blocking/major item is settled by the repo code and the doc sources.

## Verdicts

- cuda-1 | CONFIRMED (blocking) | 00_front "6.3" says after E6 only the c_speed/b3_v5 lines change ("Nothing else in this section changes"), yet the V3 lines keep `run.sh <run> experiments/phase2/configs/d_gpu_test.yaml` and `.../a3_test.yaml` (and V2 keeps `a3_measure.yaml`), while E I29/R-4 say `chain()` reads only `specs/*.yaml`, E6 step 2's live/finished list omits `a3_test`, no source names `specs/d_gpu_test.yaml` or `specs/a3_test.yaml` (grep: 0 hits), and `cmd_run` calls `chain(ROOT / a.module, ...)` only when `run.sh` runs after launch [code: experiments/phase2/vast/policy.py:471]; two rules disagree and the box bills while plans are ported by hand.
- cuda-2 | DOWNGRADED (minor) | The code part is true (`CONSTRAINTS` is set only from `/opt/requirements.txt` [code: run_box.sh:22-23], which only Dockerfile.phase2:37 writes; Dockerfile.vllm_native:20 copies other files), so the "(so both venvs load the same numpy)" sentence is wrong, but no measurement shows that a numpy version difference moves a D-G7 median by the 5 % margin; the gate flip is speculative, so this is a wording fix plus a nice-to-have version check.
- cuda-3 | CONFIRMED (major) | D "D8" `--gpu-check --dry` asserts through `experiments.phase2.vast.policy.code_files()` once the tar exists (the V3 pre-launch gate), while E "E8" step 1 splits policy.py into `objdel.infra.vast.*` (name table: `code_files` -> `infra/vast/remote.py`) and rewrites only the `vast/*.sh` shims; `run_box.sh` is not named and U8 scans only `objdel.` literals in Python, so on a post-E8 tree (allowed by front "End state under (b')") the D8 gate fails with ModuleNotFoundError and no source says what replaces it (same class as round-9 cuda-2, confirmed major).
- cuda-4 | DOWNGRADED (minor) | The pre-launch step is already in the doc: D "D8" Steps ("On the Mac before the V3 launch ($0 ...): `make old-engine-tar`, then `make gpu-check DRY=1`") and the D9 `docs/runbooks/gpu-run.md` list both, and the box step 0 stops in seconds; what remains is the dead preflight-hook sentence (stamped V3 config is `b3_launch.yaml`, front 6.1 step 4) and a missing cross-reference in front 6.1, which D's own log records as a residual.
- cuda-5 | kept minor | Wording ("no forked worker opens its own CUDA context" vs one worker that does); no outcome change.
- cuda-6 | kept minor | Test-design detail for D-G15's growth check; an engineer sees the test fail on the Mac at $0.
- cuda-7 | kept minor | Sanitizer coverage wording; D-G5 runs all four tools on the port anyway.
- cuda-8 | kept minor | Nice-to-have sanitizer flag and lane note.
- cuda-9 | kept minor | Wording about where nvcc sits; conclusion unchanged.

Counts: confirmed 2 (blocking 1, major 1), refuted 0, downgraded 2, kept minor 5.
