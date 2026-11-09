# Verdict on final_r4_arch.md (fresh skeptical verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item needs an NVIDIA product or catalog skill; no lookup, nothing recommended. |
| engineering:code-review | Checked each cited code line in the repo: serve.py:50-125 (STORE AnnAssign, _POOLS guard), dataset.py:65,131,172-175, select_ckpt.py:52-54,187-220,271,354,397-403, g_probe.py:44-48,81-93, train_sft.py:390, preflight.py:215-276. |
| engineering:debug | Traced item 1 through import-linter's exhaustive rule (docs fetched: an unlisted container child fails the contract) and item 7 through preflight `main` (a crash sets `pass` False). |
| engineering:system-design | Judged items 6 and 9 by whether the R-3 gate or an on-disk record removes the guess (both do), which lowers them. |
| engineering:architecture | Checked items 2 and 3 against the `rl | sft` independence and the src-never-imports-experiments rule (U5, E-G8). |
| engineering:documentation | Looked for the same fact handled elsewhere (E-G14, E.8, front 6.3, CLEAN16, U19, B.md "R5 = rl-full-3ep"). |
| huggingface-skills:huggingface-best | No model choice in these items; no change. |
| WebFetch | import-linter layers docs: "If, say, a module existed called `mypackage.foo.extra`, the contract will fail as it is not listed as a layer." |

## Verdicts

- arch-1: CONFIRMED (blocking). E.md E8 step 4 and U1 fix `exhaustive_ignores = ["_gpucheck"]`, but E4 step 3 and E-G14 keep `src/objdel/train/` (sample.py) until the optional post-E V4 follow-up, and import-linter fails an exhaustive contract on any unlisted container child, so `make lint` (R-3) fails at E8: two rules that both claim to hold.
- arch-2: CONFIRMED (major). E.md E2 step 6 moves "the D12 selection logic of select_ckpt.py" into `sft/select.py` with no name table, while that logic calls `experiments.phase2.prereg` (select_ckpt.py:397-403), `baselines` (271-290), `train_sft.vision_check` (354) and uses a `__file__` ROOT (54), which U5/E-G8 ban in src; the engineer must guess which names move and where prereg/baselines live.
- arch-3: CONFIRMED (major). E.md E5 table puts `TraceStep` in `rl/callbacks.py` and `stop_ids` in `rl/reward.py`, but train_sft.py:390 (moving to `sft/train.py`) uses `callbacks.TraceStep` and select_ckpt.py:187-200 (moving to `sft/select.py`) calls `rewards.stop_ids`, so `sft -> rl` breaks the stated `rl | sft` independence and the engineer must pick a new home.
- arch-4: CONFIRMED (major). E.md E4 step 3 calls `sample()` "the only external `_generate` caller" and deletes it in the V4 follow-up, but g_probe.py:45-48 (a live CLI, chain stage `probe`) calls `sample.sample(trainer, rows, 8)` lazily, so g_probe breaks at run time and `export`, `write`, `import_completions` get no stated home.
- arch-5: DOWNGRADED (minor). `STORE` is indeed missing from the `verify/store.py` row and the AnnAssign is invisible to U10, but its home beside `preload`/`CacheRef` is obvious (they fail with NameError otherwise, caught by the serve tests), and CLEAN16 step 1 deletes only `_Pool, fork_pool, pool_size, cpus, quota`, so `_POOLS` stays and `preload` does not break.
- arch-6: DOWNGRADED (minor). The change-before-move order for `prepare(prompts=...)` is unstated, but a move-first attempt fails at import time (dataset.py:131 PROMPTS vs dataset.py:175 `_pixels_sha`) under the R-3 gate on the Mac, so the engineer reorders with no wrong result and no spend.
- arch-7: CONFIRMED (major). CLEAN.md CLEAN8 step 4 archives `project_r1.py` in the E6.6 commit, but preflight.py:217 `projection` imports it first thing and preflight.py:269-276 turns the crash into `pass: False`; no source removes the `projection` check, so every later preflight fails and no rental can launch until the engineer finds and fixes the missing step.
- arch-8: CONFIRMED (major). CLEAN.md CLEAN0 step 4 and CLEAN-G13 hard-code `make preflight CONFIG=experiments/phase2/configs/r1_qwen3vl2b.yaml`, re-run at Wave 3 (CLEAN0 step 5), but E.md E6 step 6 archives that file and E6 step 8 makes preflight take a RunSpec; the Wave 3 gate cannot run as written (untestable gate).
- arch-9: DOWNGRADED (minor). V4 has no named spec file, but every input is pinned: E4/E.8 say "same checkpoint, view, K and seed" as `runs/r5/eval/eval_val/final_eval.json` (ckpt R5 final = `rl-full-3ep` per front 6.3/B.md, view `lidc-val-shifted`, K 8, seed 11 [measured]), and E.9 gives the pytest list; writing the spec is mechanical.
- arch-10 to arch-18: minor as filed (no change).

Counts: confirmed 6 (1 blocking, 5 major); refuted 0; downgraded 3 (items 5, 6, 9).
