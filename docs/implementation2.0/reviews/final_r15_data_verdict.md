# Verdict: final review round 15, data

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Loaded. No item depends on an NVIDIA product; no catalog lookup. |
| engineering:code-review | Read `experiments/phase2/train_rl.py:51-68` (normal return writes `final/` and `rl.json` with `steps` and no `partial`) and `experiments/phase2/vast/policy.py:374-383` (`drop = re.compile(r"^(data|docs)/|...")`). |
| engineering:debug | Traced both paths: BudgetStop trip -> `trainer.train` returns -> plain `rl.json`; `--heldout` on the box -> `val_half` -> `docs/reference/val_halves.json` -> not in `code_files()`. |
| engineering:system-design | Checked which files a box process can read (code sync + bundle) against the B1.0 storage path of the val-half table. |
| engineering:architecture | Checked one owner per rule: B1.6 keep regex, CLEAN1 `kept_docs`, D1 step 0 `keep` all name only `frozen_views.json` (plus kernel drafts). |
| engineering:documentation | Each verdict quotes the source line and the code line. |
| huggingface-skills:huggingface-best | No model choice in scope. |

## Verdicts

- data-1, downgraded (minor): the gap is real (`train_rl.py:61-68` writes `"steps": done` with no `partial` key after a `BudgetStop` return, and 00_front 1.1 keys the row on "`rl.json` has `partial: true`"), but it is the same label gap r14 fresh-arch-8 filed as minor ("The run then ends without `partial: true`, so no 1.1 row labels the shorter run"), and r14 rl-1's verdict calibrated the missing "partial run" label as "a consistency nicety, not a wrong result": Gate B1/B3 are held-out val rules, B-S2 already waits for "the owner's budget decision" (1.1 Pass row), and the final candidate chosen on val-select is not wrong; only the v5-vs-B1 attribution gets an unstated step confound.
- data-2, confirmed (major): B1.0 stores the val-half table at `docs/reference/val_halves.json` and `val_half` "is then a lookup in that table"; B1.5's optional held-out curve ("The eval set is 64 items of **val-shape** only") and I44 ("a planted val-select id is refused") need `val_half` inside `train_rl` on the box, but `code_files()` drops `^(data|docs)/` [code: policy.py:379-380] and every keep rule (B1.6 `keep = re.compile(r"^docs/reference/frozen_views\.json$")`, CLEAN.md:1986 `kept_docs`, D.md:331) keeps only `frozen_views.json` and kernel drafts; I44 and E08 run on the Mac where `docs/` exists and I50 checks only `frozen_views.json`, so with owner item 15 on (B.8.3, $0.02) the V1 `train_rl` fails at start, step 3b finds no checkpoint and the chain stops after about 53 min of paid start-up, and no planned $0 check catches it.
- data-3: minor as filed (no change).
- data-4: minor as filed (no change).
- data-5: minor as filed (open items re-listed, not re-graded).

## Counts

confirmed 1 (data-2), refuted 0, downgraded 1 (data-1).
