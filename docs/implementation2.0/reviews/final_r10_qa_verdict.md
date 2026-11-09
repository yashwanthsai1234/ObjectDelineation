# Verdict on final_r10_qa.md (fresh verifier, 2026-11-09)

## Skills used
| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in scope; no catalog skill recommended. Pointed the check at the `score_masks` device-table binding (item 1). |
| engineering:code-review | Read `cuda/gpucheck_py.cpp:104-106, 213-224`, `cuda/gpucheck.cpp:1174-1185`, `src/objdel/reward/score.py:67-77`, `experiments/phase2/vast/policy.py:422-497` against each claim. |
| engineering:debug | Traced `chain` -> `cmd_run`: one `&&` job under one `timeout`; no per-job kill command exists (only `teardown`), so an in-plan stop rule cannot fire (item 2). |
| engineering:system-design | Checked V2/V3 tables, launch lines and owner totals (items 2-4). |
| engineering:architecture | One owner per rule: C9 home (C9 Files, C0, C.8.2, G-CF1 runner) and B-S2 (front table, front bullet, B.8.1). |
| engineering:documentation | Checked the V2 Mac test path against the V3/V4 env-var rule (item 5). |
| huggingface-skills:huggingface-best | No model choice in scope. No change. |
Web tools loaded (WebSearch, WebFetch); no external claim needed checking.

## Verdicts (blocking/major items)
- qa-1, confirmed (major): B.md B.9.5 G01 says "Every other field ... equal tests/goldens/b_g01.json" and the B.9.8 sketch calls `G.score_masks(m[None])  # ... (adapt)`, but `score_masks(bits, hw, start, tables)` needs one table handle per object in the object's own frame [code: cuda/gpucheck_py.cpp:213-224; cuda/gpucheck.cpp:1181-1182], and no source text says which rater stack builds those tables, so the golden-compared fields depend on an engineer's guess (the `np.nonzero` sum part is fine).
- qa-2, confirmed (major): front 6.3 V2 item 3 stop rule "Step-1 rewards differ: stop the fast arm" and item 3c "only if mode A fails C-G-MEM" sit inside one `c_speed.yaml` run, but `chain` returns the steps literally and `cmd_run` joins them with `&&` in one `bash -c` [code: policy.py:432-433, 487-494] and policy.py has no per-job stop command (only teardown), so neither condition can act and the engineer must guess what goes in the plan.
- qa-3, confirmed (major): C.md C9 "Files: `experiments/phase2/configs/c_speed.yaml` arms" and C0 "each 17-step arm (baseline, 3c, 3d, C9)" contradict C.8.2 "C7, C8, C9 and C11 ... ride inside the Phase B runs" and the G-CF1 runner row "inside the C9 arm (a B run)"; no V2 or V3 row prices C9, so the two readings give different V2 plans and spend.
- qa-4, confirmed (major): B-S2 has three item sets (front 6.3 table "items 0 + 4", front bullet "items 0, 2a and 4", B.md B.8.1 "items 0, 2a, 2b and 4"); option (a) prices B-S2 at $2.17 (items 0 + 4) while the V3 launch lines run items 1 and 3 with no approval condition, and 6.4 rule 5 then puts the B3 BudgetStop where the run ends well before 260 steps, so Gate B3 compares unequal runs.
- qa-5, confirmed (major), with a corrected fix: G-RP2/G-RM1 run "after `pull.sh` of V2" but no text gives the Mac run folder (V3/V4 have `OBJDEL_V3_RUN`/`OBJDEL_V4_RUN`; V2 has none) and the check arm's `--record-batches DIR` is never named; the review's quote is partly wrong (C1 already names `{rd}/c_check/batches/b{0..4}.pt` for the `c_speed_check` batches), and its fix would put the check-arm batches in that same folder, which `--record-batches` refuses (non-empty dir -> FileExistsError, C.md C.9 item 11), so the fix must use a separate folder (for example `{rd}/c_check/arm_batches`) plus an `OBJDEL_V2_RUN` rule.

## Minor items (not re-graded upward)
- qa-6 .. qa-11: stay minor as filed (wording, stale numbers, unmeasured sanitizer time, robustness of a record file).

Totals: 5 confirmed (all major), 0 refuted, 0 downgraded.
