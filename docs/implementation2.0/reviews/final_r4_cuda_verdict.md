# Verdict on final_r4_cuda.md (round 4, CUDA and GPU performance review)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers rental timing, CMake build options or test payloads; no verdict depends on it. Nothing installed. |
| engineering:code-review | Read the repo lines each item cites: `CMakeLists.txt:15` (`OBJDEL_CUDA ... OFF`), `policy.py:399` (editable install, no CUDA define), `policy.py:432-452` (`chain`: `steps:` returned literally; build only for scoring stages), `r5_shared.yaml` (no build step), `runs/r5/env_check.json` (`nvcc present false`), `checker_gpu.yaml` header (code-only sync ships no data), `data/checker_eval/` (masks.npz, polygons.npz, manifest only). Confirmed items 2 and 4. |
| engineering:debug | Reproduced each failure path on paper: CT-72 at the end of Wave 2 after D1/D2 edit `cuda/` (item 1); A.9 hook skips `gpu` tests on a CPU build so G02/G08 exit 0 (item 2); D-G8 replay reads `data/prepared/caches/lidc`, absent from R-D2 (item 4). |
| engineering:system-design | Recomputed the V1 timeline from `ledger.csv` launch and `run` rows (17.7 / 24.1 / 38.2 min) against the 30-min item 1 and the uncounted 6-min `train_rl` start-up (item 3). |
| engineering:architecture | One-rule check: CLEAN.spec C1 and CT-72/CT-55 (diff against the Wave 1 tag) against D1/D2/D8 edits of the same paths (item 1). |
| engineering:documentation | Each verdict names the source file and heading and the proving code line or measurement. |
| huggingface-skills:huggingface-best | Not used for any verdict (model choice not in dispute). |

## Verdicts

- cuda-1 | confirmed (blocking) | CLEAN.md L7 CT-72 ("run at the end of each wave": "`git diff --name-only pre-clean-2026-11-08 -- cuda CMakeLists.txt` is empty (CLEAN.spec constraint C1)") and L5 CT-55 (CLEAN16.1, Wave 3, gate CLEAN-G2: `git diff --stat pre-clean-2026-11-08 -- cuda/ tests/test_cuda_checker_gpu.py tests/test_gpu_scoring.py` empty) diff against the tag made once in CLEAN0 step 6 (Wave 1, CLEAN.md:332), while D1 changes `cuda/gpucheck.cpp` and `CMakeLists.txt` (D.md:238, 245), D2 deletes engine symbols (D.md:317) and D8 edits `tests/test_cuda_checker_gpu.py` (D.md:674-676) before Wave 2 ("after phase D3") and Wave 3 ("after phases D and E"); no WAVE guard exists, so both checks fail on a correct tree.
- cuda-2 | confirmed (blocking) | Front 6.2 "Plan files" and B1.5 list `b1_rows.yaml` as env_check, G02/G08 smoke, train_rl, final_eval with no build step; `sync_up` installs with `uv pip install ... -e .` (policy.py:399), `option(OBJDEL_CUDA ... OFF)` (CMakeLists.txt:15), `chain` returns `steps:` literally (policy.py:432-433) and adds `run_box.sh --build-only` only for scoring stages (policy.py:450-452); R5's `steps:` plan had no build and its box reports `"nvcc": {"present": false}` (runs/r5/env_check.json); so V1 rewards and evals score on the CPU build (owner rule "CUDA ... GPU everywhere") and A.9's hook (A.md:1236-1238) skips G02/G08 so the smoke passes by construction. Same defect as qa-2.
- cuda-3 | confirmed (major) | Front 6.2 item 1 puts launch, sync, downloads, env_check, smoke AND pull/push/teardown at 30 min, and the note says the stop "trips only above about 28.5 s per step [estimate: (172.3 - 30) x 60 / 300]", but the front's own inputs table says the ~6 min `train_rl` start-up is "not in the step minutes" and is in no V1 row, and ledger.csv measures launch to first `run` row at 17.7 (rl), 24.1 (r4) and 38.2 min (r5); pre-step time is then about 30-50 min, so `BudgetStop` ($2.05, session minute 172.3) ends the B1 run at about 266-297 steps instead of the planned 300 the owner approves; money is safe, so major, not blocking.
- cuda-4 | confirmed (major) | D1 defines R-D2's payload only as code-only sync plus `data/checker_eval` (D.md:266-267; repo `data/checker_eval/` holds only masks.npz, polygons.npz, manifest), but D-G8 replays on caches from `data/prepared/caches/lidc` (D.md:292-294, 786), D-G10's box half reads 1,000 LIDC caches (D.md:788, 592), and GPU02 needs "npz golden saved on the Mac" with no path (D.md:1024) and cannot ride the tracked-code sync if over 1 MB (CLEAN.spec C8); D8 step 3 fails on any skip (D.md:714) and `--gpu-check --dry` only collects, so V3 item 1 money can be spent with an undefined input.
- cuda-5 | minor (as filed) | Wording of "no nvcc" and the wheel swap; no wrong result for an engineer who follows the doc.
- cuda-6 | minor (as filed) | `nsys status -e` exit code is undocumented; the fallback to torch.profiler still gives a result.
- cuda-7 | minor (as filed) | Range-name list in a diagram.
- cuda-8 | minor (as filed) | One-cent rounding.
- cuda-9 | minor (as filed) | Stale size term; the bound still holds.
- cuda-10 | minor (as filed) | Pinned-index detail; G-PX2 would expose it on the Mac.

Counts: confirmed 4 (blocking 2, major 2), refuted 0, downgraded 0.
