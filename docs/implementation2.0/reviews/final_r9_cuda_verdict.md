# Verdict: final_r9_cuda.md

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill needed; CUDA facts checked against repo code instead. |
| engineering:code-review | Read run_box.sh:10-30, policy.py:374-435 and 480-495, Makefile lint target to test each claim. |
| engineering:debug | Traced the bench's input paths from cwd /workspace/src to the bundle at /workspace/data/lidc (item 1). |
| engineering:system-design | Checked the front time line (items 5-7, "End state under (b')") for V3 running after E (item 2). |
| engineering:architecture | Checked E3a's importer rewrite and TID251 per-file-ignores against D-G7's pre-E old venv (item 2). |
| engineering:documentation | One line per item with the proof. |
| huggingface-skills:huggingface-best | No model choice in scope; nothing changed. |

## Verdicts

- cuda-1 | CONFIRMED (blocking) | D.md "D1" step 2 names the bench inputs as `data/checker_eval/masks.npz`, `data/prepared/lidc/val.jsonl`, `data/prepared/caches/lidc`, and the D1 file list for `bench_checker.py` names no data-root flag or variable; the only `OBJDEL_DATA_ROOT` readers in D are the pytest payload check and `box_input()` (D.md:805, 826, 867); on the box cwd is `/workspace/src` [code: run_box.sh:16; policy.py:489 `cd /workspace/src`], `data/` is dropped from the code sync [code: policy.py:379], and the bundle unpacks at `/workspace/data/lidc` [code: policy.py:406-433], so step 6's `--dump-replay` and crops workload fail after the build, pytest and sanitizer minutes, and `--dry` (collect-only, S08, S10 `--synthetic`, tar check) never runs the real inputs; D-G7/D-G8 are lost on a paid rental.
- cuda-2 | CONFIRMED (major) | E.md E3a says "rewrite every direct importer to `objdel.geom.engine` ... re-grepped on the E0 commit" with TID251 banning `objdel._gpucheck` and per-file-ignores listing no `experiments/phase2/perf/bench_checker.py` (E.md:486-498, 630-632), while `make lint` covers `experiments` [code: Makefile:33]; the front lets V3 run on a post-E tree ("If E0 has started, the owner rebases the branch onto post-E paths", 00_front.md "End state under (b')"), and D-G7's old venv is the pre-E `pre-d2` tree (today `src/objdel` has no `geom` [measured: ls src/objdel, 2026-11-09]), so `--engine old` would raise ModuleNotFoundError; no D or E text names this exception.
- cuda-3 | minor (as filed; not promoted) | Alternation wording between D-G7 and D8 step 6 is a measurement-quality nit, not a wrong result.
- cuda-4 | minor (as filed; not promoted) | Counter-key timing in D.4 is a clarity gap; the compare's key set is a small implementation detail.

Counts: confirmed 2 (blocking 1, major 1), refuted 0, downgraded 0.
