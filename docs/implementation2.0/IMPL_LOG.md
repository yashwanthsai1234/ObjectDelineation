# Implementation log (deviations, flags, open items)

One line per deviation from Implementation.md: `- <task>: <what changed> (why)`.

## Flags kept (defaults = B1 control)

## W0
- DOC1: step 0 copies were already committed (ab86df3); `git mv` used for docs/plan because research2.0 was already tracked, not untracked (why: plan assumed untracked)
- DOC1: test_no_reference_to_the_old_plan_path skips docs/implementation2.0/ (why: the plan names the old path as history)
- DOC1: scratch-path scan test skips docs/implementation2.0/notes/research/ (why: other agents write it in parallel); 5 `SP/` hits in notes/B3.md, C.v1.md and 3 reviews became repo-relative cites
- DOC1: owner-run pre-publication scan (b), RIGA figure OK and `git add` list not run by the agent; baseline in runs/doc/baseline.txt is collect-only (1497 tests), `make test` summary skipped (why: only touched tests run in W0)
- DOC2: TRL link/text replaced in 8 files incl. Dockerfile.phase2; lines rewrapped to stay under ruff E501 (why: longer replacement text)
- DOC2: banner dates are each file's last git commit date (why: plan says YYYY-MM-DD, no source given)
- CLEAN0: tag pre-clean-2026-11-08 created locally at the DOC2 commit, not pushed (owner pushes); before.json also has vulture_lines and target_set_sha256; tests/test_clean_baseline.py added (skips when runs/ is absent)
- CLEAN0: `make test` run twice, 1460/1463 passed (3 CLEAN0 tests red before the records existed), 6 skipped, 31 deselected; baseline in runs/clean/baseline.txt; preflight r1_qwen3vl2b passed (10/10 checks)

## L1 eval

## L2 data

## L3 cuda

## L4 train

## L5 clean

## W1b

## W2

## E

## Break-fix

## GPU rental

## DOC and CLEAN final
