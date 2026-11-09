# Verdict on final_r6_data.md (fresh skeptical verifier, 2026-11-08)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in these items; no catalog lookup needed. Nothing changed. |
| engineering:code-review | Read `dataset._prepare:114` (`rec["disagreement"] = len(rater_answers(caches, names)) > 1`, per record), `dataset._check/_verified:170-180`, `build_eval_views.py:93` (manifest `files` hashes only the jsonl), and the `.venv` `huggingface_hub` `upload_folder` docstring (multi-commit "(part N)"). Confirmed items 1 and 3. |
| engineering:debug | Re-ran item 1's measurement on `data/prepared/lidc/train.jsonl` (`python -I`): 1,247 nodules, 100 with a disagreement slice, 100 of those mixed, 1,394 items in them, 289 / 7,341 flagged items. Counted 14,684 files in `data/views/lidc-train-shifted`. |
| engineering:system-design | Traced the rows root Mac build -> HF upload -> box download -> `train_rl` start-up: the only checks are the local-folder preflight and the box `manifest.json` SHA; image files are first checked in `_verified` after about 53 min / $0.70 (front 6.2 item 1). |
| engineering:architecture | Checked the order: front time line step 7 puts E after the last A-D commit, and B1.8 runs on two of the four section-1.1 rows, so E must handle both the deleted and the kept multi-root loader; E's move table and U33 handle only the kept case. |
| engineering:documentation | Searched E for a "re-grep / absent on E0" rule that covers `VIEW_SEP`/`_tagged`/`train_views` (exists only for `_gpucheck` importers and the fork pool, R-2); searched B/front for a remote-revision check (none). |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |

## Verdicts

- data-1 (major): CONFIRMED major. B "B.9.1 Unit tests" U11 says "`Window.disagreement` is True iff the doctors' topologies disagree (the repo's rule for prepared records); equal across all windows of one nodule", but the repo rule is per slice [code: src/objdel/train/dataset.py:114] and 100 of 100 train nodules with a disagreement slice also have agreeing slices (1,394 items) [measured: python -I on data/prepared/lidc/train.jsonl], so the test spec contradicts itself and a nodule-level reading moves the B1.2 quota from 14.1 % to about 48 %.
- data-2 (major): CONFIRMED major. B "B1.8 Clean after gate B1" deletes `serve.VIEW_SEP`, `_tagged`, `grpo_rows_views` and `--train-views`, while E "E3b" moves `VIEW_SEP`, `_tagged` into `verify/store.py`, E.9.1 U33 asserts `VIEW_SEP == "::"` and `loader(r1, r2)`, and E6 `RlStage` carries `train_views` "(train_rl --train-views)"; E runs after Gate B1 (front time line step 7) and B1.8 runs on two of four rows, so on those rows the two texts give different code and E has no rule for the absent names (CLEAN CT-33/CT-35 staleness alone would be minor, since `make test` exposes it).
- data-3 (major): CONFIRMED major. B "B1.4" step 5 says "Record the HF commit as `rows_revision`", but `upload_folder` commits a large folder "in several batches ... (part N)" [code: .venv huggingface_hub/hf_api.py upload_folder docstring] and the root is about 16,640 files (today's 7,341-window view: 14,684 files [measured: find | wc -l]); B1.3 change 2 preflight reads only the local folder, the B1.5 box check reads only `manifest.json`, and the manifest `files` field hashes only `lidc/train.jsonl` [code: build_eval_views.py:93], so a wrong or partial pinned revision is first caught by `_verified` [code: dataset.py:176-180] after the ~53 min / $0.70 front 6.2 item 1 block.
- data-4 (minor): stays minor as filed (test-row key count wording).
- data-5 (minor): stays minor as filed (test sentence in the wrong section; R03 tests DG10 on the schedule).
- data-6 (minor): stays minor as filed (stale step count in an owner-decision copy).
- data-7 (minor): stays minor as filed (A3 wording; B.6.2 holds the rule).
- data-8 (minor): stays minor as filed (table precedence wording).
- data-9 (minor): stays minor as filed (report-only label-noise line).
- data-10 (minor): stays minor as filed (download time estimate; nice-to-have).
- data-11 (minor): stays minor as filed (still-open wording conflict from r5 data-4).
- data-12 (minor): stays minor as filed (U59 already runs `final_eval` on a view fixture and would catch it).

Counts: 3 confirmed (0 blocking, 3 major), 0 refuted, 0 downgraded. The 9 minor items keep their minor label.
