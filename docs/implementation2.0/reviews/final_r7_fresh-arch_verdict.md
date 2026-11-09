# Verdict: final round 7, fresh-arch review

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA tool. I did no catalog lookup and installed nothing. |
| engineering:code-review | I read `experiments/phase2/train_rl.py:20-70`. `--train-views` is at line 26 and `roots = [a.train, *a.train_views]` at line 47, which feeds `dataset.grpo_rows_views`. A grep shows `--train-views` only in `rl_shifted.yaml:26`, `train_rl.py`, `g_probe.py` and two tests. |
| engineering:debug | I traced the "Pass" row: V1, Gate B1, B1.8 on `main`, a raise, the branch merges, V2. C's arms still pass `--train-views` (C.md C2 tests, C.8.2; front 6.3), but B1.8 removes it. B1.8 also removes the path without a schedule, and B.md B1.3 refuses `--steps` that is not a multiple of `block_steps` (B.md:560, 825), so 17 and 32 fail too. |
| engineering:system-design | I checked the end state of each decision-table row. Row 3 ("no data effect") has no end condition for "B1.8 waits" (front 1.1, B1.8, 99_back decision 2b). |
| engineering:architecture | For items 2 and 3 I looked for the rule that governs. E's own "When E starts" and R-4 define the order, and the front "End state under (b')" calls itself "the one rule". D.4 says D2-D8 do not wait for a GPU run. The section 5 kernel rule quotes K3 DESIGN section 9, and that table covers K3's parts only. |
| engineering:documentation | I checked every quote against its source file and heading. |
| huggingface-skills:huggingface-best | No item is about choosing a model. No change. |
| WebSearch/WebFetch (loaded) | Every item is about the doc agreeing with itself or with the repo, so I made no web query. |

I graded only the blocking and major items (1-4). The review already rates items 5-10 as minor, so I did not grade them again.

## Verdicts

- fresh-arch-1: CONFIRMED (blocking). The front, section 6.3 "Own batches", says "The check arm of item 3 reads the same two views (`--train-views`)". `c_launch.yaml` "has no `rows_plan`: the C arms read the bundle's train views". C.md (C2 tests and C.8.2) trains the check and fast arms on `--train-views`. B.md "B1.8 Clean after gate B1" deletes `--train-views` [code: train_rl.py:26,47] and "the 'without schedule' branch" on the "Pass" and "far-miss only" rows. U42 and CT-33 even assert argparse exit 2 after B1.8. V2 runs only after a raise merges `after-raise` into `main`, which by then holds B1.8. So item 3 of V2 dies after items 1-2 are paid ($0.58 + $0.15). B1.8's grep cannot catch `c_speed.yaml`, because that file sits on the branch.
- fresh-arch-2: DOWNGRADED to minor. The "no post-E exception" bullet in E.8 and the E5 rows without "if landed" are stale text. The same phase gives the rule correctly in "When E starts" and in R-4, and the front "End state under (b')" calls itself "the one rule". `old_names.json` is recorded on the E0 commit, so a source file that does not exist adds no names, and `test_move_complete` does not fail. Whether `count_syncs` lands on `main` changes no result under (b').
- fresh-arch-3: DOWNGRADED to minor. Front section 5 quotes K3 DESIGN section 9, and that decision table covers only K3's parts ((a), (b), (c)). D.4 says plainly: "D2-D8 use the CPU build as the oracle, so they do not wait for R-D1. The R-D1 GPU numbers gate only speed claims (D-G7 ...)". The front end state lands D1-D9 with their box gates "open until a rental". The rule that governs is clear. Making section 5's sentence name its scope is a wording fix.
- fresh-arch-4: CONFIRMED (major). On row 3 of the decision table (front 1.1, B.md B1.8, 99_back decision 2b), "B1.8 waits; the rows code stays" has no trigger and no end. After B1.5 archives `rl_shifted.yaml` (the only config with `--train-views`, `rl_shifted.yaml:26`), no live plan on `main` uses `train_rl.py:26,47`, `dataset.grpo_rows_views`, `serve.VIEW_SEP` or `g_probe.py:90`. They stay beside the rows path through E and CLEAN, which conflicts with the owner's archive and dead-code rule. The engineer must guess whether and when to delete them. The fix should go with fresh-arch-1: once C's arms use the rows root, the old path has no user after B1.5.
- fresh-arch-5..10: not graded again (the review rates them minor).

Counts: confirmed 2, refuted 0, downgraded 2.
