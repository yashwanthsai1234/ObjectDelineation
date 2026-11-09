# Verdict on final_r5_data.md (fresh skeptical verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in these items; no catalog lookup needed. Nothing changed. |
| engineering:code-review | Read `experiments/phase2/hf_pull.py` (no `--subdir` flag today; `snapshot_download(..., local_dir=a.out)`), `policy.py` `cmd_fetch_hf` and `{data}` = `/workspace/data/lidc`. Confirmed item 2. |
| engineering:debug | Reproduced item 2 from the library source: `.venv` `huggingface_hub/_snapshot_download.py:152` "If `local_dir` is provided, the file structure from the repo will be replicated in this location". |
| engineering:system-design | Traced the V1 rows path (B1.4 upload to `rows/lidc-train-b1/` -> `fetch_hf --rows` -> `train_rl --train {data}/rows/lidc-train-b1`). |
| engineering:architecture | Compared the three after-V1 tables (front 1.1, 99_back 2b, B.4.1) with B.6.2 and A3 step 5 (item 1). |
| engineering:documentation | Checked whether any other section defines a Gate B3 arm for B1.9 (none; item 3). |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |

## Verdicts

- data-1 (blocking): CONFIRMED. 00_front "1.1 After V1: decision table" row 4 and 99_back decision 2b say "view identity fails: Stop. Owner review ... R5". B.6.2 "View identity (from A-G3)" says "If this correlation rule fails, Gate B1 still pairs with A3's items ... only the R5 report line is dropped", and A3 step 5 says the same. B.9.8 `test_gate_b1_needs_all_hard_rules` asserts that `("view_corr", 0.74)` makes `gate_b1(...).ok` false. B.4.1 row 4 has only "edge guard or legal share" and no row for the offset < 0.15 guard. U49 lists that guard as a hard rule. So two texts that each claim to be the rule give different outcomes: fund B3, or stop with R5 as the final candidate.
- data-2 (blocking): CONFIRMED. B "B1.5 Rows-root download" says `hf_pull --subdir rows/lidc-train-b1 --out /workspace/data/lidc/rows/lidc-train-b1` with `snapshot_download(..., allow_patterns="rows/lidc-train-b1/**")`, and I40 fixes the output folder at that path. `--subdir` exists nowhere in today's code [code: hf_pull.py:18-22] and the doc never defines it. `snapshot_download` with `local_dir` copies the repo paths [code: .venv huggingface_hub/_snapshot_download.py:152]. So the literal reading puts the files at `.../rows/lidc-train-b1/rows/lidc-train-b1/`. Then `train_rl --train {data}/rows/lidc-train-b1` [code: policy.py:433] finds no manifest at b1_rows.yaml step 3, after paid box time. The Mac test (I40, Hub mocked) does not catch this.
- data-3 (major): CONFIRMED. B1.9 says "apply every Gate B1 and Gate B3 rule ... Print P(pass) per rule and jointly", and "joint P(pass) < 0.5" starts an owner rule-change decision. But Gate B3 pairs a B3 arm against B1 on val-select (B.6.2), and B1.9 simulates only a B1 arm against R5. No other section (grep `gate_power`, "B3 power") and no test (E10, E10b) defines a B3 arm, so engineers would compute the joint value differently.
- data-4 (minor): stays minor (reviewer label). It is a wording conflict between CR1 and the percentile interval in B1.9. It is not re-escalated.
- data-5 (minor): stays minor. It is advice to the owner. The doc already has the branch "joint P(pass) < 0.5: the owner decides before the rental".
- data-6 (minor): stays minor. It is a stated-share accuracy issue in DG14. DG14 is a report gate (hard=False, U19).
- data-7 (minor): stays minor (reviewer label). The B1.0 difference branch contradicts DG15. The reviewer measured this branch as unlikely: 40 of 40 rebuilt items match.
- data-8 (minor): stays minor. Two codings differ, but the outcome does not change (1,183 < 1,436).
- data-9 (minor): stays minor. It is a stale number in a diagram.
- data-10 (minor): stays minor. It is item-list wording. The skippable items are explained elsewhere.

Counts: 3 confirmed (2 blocking, 1 major), 0 refuted, 0 downgraded. The 7 minor items keep their minor label.
