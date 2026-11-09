# Verdict: final_r9_data.md (fresh verifier, 2026-11-09)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA product in scope (headline publication and data statistics); no catalog lookup; nothing changed.
- engineering:code-review: read `experiments/phase2/final_eval.py:41` (writes bare `ci95`) and grepped every `ci95` reader in src/experiments/tests (only final_eval.py and tests/test_final_eval.py:24-26), as A1 step 3 claims.
- engineering:debug: reproduced the path claim with `ls runs/r5`: the box writes `runs/r5/eval_test/` (only `trainer/` left there), and the pulled files sit at `runs/r5/eval/eval_test/final_eval.json`, so the `eval/` level was added by hand after the pull.
- engineering:system-design: traced A1 schema 2 -> A3 `--out {rd}/eval_{val,test}` -> DOC6 `headline.py` (`runs/<ID>/eval/eval_{test,val}`, `ci95`) -> DOC11 card.
- engineering:architecture: checked ownership: DOC "Docs per phase" row A and DOC Dependencies ("Phase A provides ... a new `headline.json` (via `make headline`, `ci_unit: patient`)") give A the job, but A9 lists no headline step and DOC never mentions schema 2.
- engineering:documentation: one line per item with its proof.
- huggingface-skills:huggingface-best: no model choice in scope; nothing changed.

## Verdicts
- data-1: CONFIRMED (major). DOC.md "DOC6" says `headline.py` "Keeps the 6 metrics with `mean` and `ci95`" and I04 refuses a file with a "metric key missing", while A.md "A1" step 3 says "The bare key `ci95` is dropped in schema 2" and DOC has no word on schema 2 / `gate_ci`; also A3 writes `{rd}/eval_val` and `{rd}/eval_test` (A.md:497, 513) while DOC6 reads `runs/<ID>/eval/eval_{test,val}` (that layout exists for r5 only because the files were moved by hand, `ls runs/r5`), so the DOC row-A deliverable "`headline.json` via `make headline` (`ci_unit: patient`)" cannot be produced without guessing the key and the paths. Not blocking: no GPU spend, no gate depends on it, and r5 (schema 1) still works. One small error in the item: A.md does have a "Words" block (A.md:11), so only the DOC6 pointer text needs a sharper anchor, not a missing-section fix. Fix touches DOC (DOC6, I03/I04, owner table) and A (A1 step 3 grep sentence / optional `--resummarize`, A9).
- data-2: minor (as filed). BudgetStop mid-block is an estimate-only path, and the reviewer states no gate result changes.
- data-3: minor (as filed). Interval-method switch at 50 patients is a consistency nit; percentile vs cluster_t at 51-52 patients does not change a stated outcome.
- data-4: minor (as filed). Replacing "not yet measured" with a capacity-bound estimate is a nice-to-have; B1.4 already measures it in the run.
- data-5: minor (as filed). Carried-over list from rounds 4-8; not re-argued.

Counts: confirmed 1, refuted 0, downgraded 0.
