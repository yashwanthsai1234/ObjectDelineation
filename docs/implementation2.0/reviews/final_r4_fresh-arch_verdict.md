# Verdict on final_r4_fresh-arch.md (fresh skeptical verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog lookup needed: no item needs a new NVIDIA tool. Item 10 (sanitizer pin) stays minor. |
| engineering:code-review | Each blocking or major item checked against the code lines: policy.py:434-490 (chain, `&&` join, soft `~`), policy.py:450 (default stages), configs.py:176-179 (num_generations_eval guard), rewards.py:132-135, bundle.py:47-60, train_rl.py:55-60. |
| engineering:debug | Traced item 5 end to end: deadline = started + cap + PULL_S, run.sh cap = deadline - now - MARGIN_S (both 600 s), so the second run.sh ends at session minute 201; found the slack is still positive with the measured 47.62 s save. |
| engineering:system-design | Judged item 4 and item 5 as guarded paths (4c gate, V2 item 1b fallback), which lowers their severity. |
| engineering:architecture | Checked item 2 against the decision-17 record: front and back each claim to hold the 4B recipe and disagree. |
| engineering:documentation | Looked for the same rule in two places (G08 row vs G01 row vs order rule 2; front 3 vs 99_back 17; A3 header vs front 6.1 4c). |
| huggingface-skills:huggingface-best | No model query needed; item 2 is a recipe/runner seam, not a model choice. |

## Verdicts

- fresh-arch-1: CONFIRMED (blocking). B.md B.9.5 row G08 includes "G01 parity" and G01 requires `sum_x`, `sum_y` (B3.1), while front order rule 2 says `v1-b1` holds no B3 code and `b1_rows.yaml` step 2 runs G08 on that tree; `grep sum_x cuda src` finds nothing and policy.py:490 joins steps with `&&`, so a failing G08 stops `train_rl`; G01's own command row puts it in V3 (contradiction).
- fresh-arch-2: CONFIRMED (blocking). Front section 3 step 0 says `branches: [2.0]` means "no branch selection runs and val stays unread", but policy.py:450 adds select and probe (val views, policy.py:441-448) whenever `branches` is set and `stages` is absent, and the 2B `vision.reference` (r1_qwen3vl2b.yaml:9) trips the run guard (policy.py:471-475) or gives a wrong G-vision baseline; 99_back decision 17 says the line must change and front was not changed.
- fresh-arch-3: CONFIRMED (major). B.md B1.5 and CLEAN.md CLEAN6 4a set `num_generations_eval=4`, but configs.py `grpo_config` raises for any value other than None or 8 and TRLReward needs whole groups of 8; the engineer must choose 4 vs 8 and the item count, which changes cost and the guard.
- fresh-arch-4: DOWNGRADED (minor). The A3 header does put the test view into data/eval_views before bundling (bundle.py globs every eval_views/lidc-*), but front 6.1 step 4c is a pre-spend gate that catches it and names the remedy; either remedy gives the same result, so no money is lost and no outcome-changing guess is needed.
- fresh-arch-5: DOWNGRADED (minor). The 47.62 s save is real [runs/r5/r5_shared.log] and uncounted, but the reviewer's own arithmetic still fits (29.8 <= 30.9 min in the worst offer; about 6.7 min slack at 27.5 s/step), the risk exists only when BudgetStop trips (step > 28.5 s vs measured 23.6 s), and V2 item 1b is the documented fallback; this is a margin refinement.
- fresh-arch-6: minor as filed (no change).
- fresh-arch-7: minor as filed (no change).
- fresh-arch-8: minor as filed (no change); runs/r5/EXIT = 137 verified.
- fresh-arch-9: minor as filed (no change).
- fresh-arch-10: minor as filed (no change).
- fresh-arch-11: minor as filed (no change).
- fresh-arch-12: minor as filed (no change).

Counts: confirmed 3 (2 blocking, 1 major); refuted 0; downgraded 2 (items 4, 5).
