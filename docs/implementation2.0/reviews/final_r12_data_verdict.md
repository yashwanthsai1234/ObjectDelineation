# Verdict on final_r12_data.md (fresh verifier, 2026-11-09)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA product choice in scope; no catalog lookup, nothing installed.
- engineering:code-review: read B1.1 `pool` signature and steps, B1.3 DG3 / `neardup_review.json`, B1.4 `build_rl_rows` flags and steps, B1.7 Gate DG17, and the U09 / U24 / I17 / I24 / I30 / E04 rows against each other.
- engineering:debug: reproduced the 27-item count instead of trusting the review (train.parquet: 380786 -> 15, 228428 -> 5, 179144 -> 4, 295628 -> 3).
- engineering:system-design: traced store -> `pool` -> schedule -> `window_record` -> rows root and looked for any hop that reads a drop decision; none does.
- engineering:architecture: checked whether another file (front, A, 99_back, review-log) owns a train-drop mechanism; A.md and 99_back defer to "B1.7 DG17 owns the rule"; the review log has no DG17 train-drop entry.
- engineering:documentation: each verdict cites the source line and the code or measurement.
- huggingface-skills:huggingface-best: no model choice in scope.

## Verdicts
- data-1: CONFIRMED (major). B.md B1.7 "Gate DG17: 0 unresolved hits. A train hit is resolved by mapping by `imageSOP_UID` or by dropping the items from train; the report says which", but B1.1 `pool(items, draws: int, seed: int = 0)` has no exclusion input and `dropped` has only "not_rl"/"shifted_none"/"rescaled_none"/"padded", B1.4 `build_rl_rows --store data/store ...` step 1 "`pool(train items, draws=n)`" reads the whole store, DG1b checks only split/group/image hash, and I17 passes "when the report marks it resolved" while I24 is a fixture with no rows root; the 4 train hit series hold 27 train items today [measured: data/store/lidc/train.parquet, 15+5+4+3]; no planned $0 check fails, so the engineer must guess the drop mechanism and a report-only "dropped" trains suspect labels in V1 (owner rule 1); not in any review log; 0.37 % of train, so major not blocking.
- data-2: minor as filed: front 6.1 step 4c is a named $0 pre-rental check that fails on a `lidc-test-*` folder [code: experiments/phase2/vast/bundle.py:53-62].
- data-3: minor as filed: label of `item_sd` affects only B1.9's record-only re-run, no gate or spend.
- data-4: minor as filed: B.spec.md drift is caught by U49 / U02 at $0 (same ruling as r11 data-2).
- data-5: minor as filed: U49 pins the absolute 20.5 %; the readings differ only near the point.
- data-6: minor as filed: carry-over list from r4-r11, not re-argued.

Counts: confirmed 1 (major), refuted 0, downgraded 0.
