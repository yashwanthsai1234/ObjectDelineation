# Verdict: final round 8, fresh-arch review

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Item 1 is about the data bundle, not an NVIDIA tool. I did no catalog lookup and installed nothing. |
| engineering:code-review | I read `experiments/phase2/vast/bundle.py:26, 50-70`. Line 54 globs every `data/eval_views/lidc-*`. Only `views` has the "train only when named" rule (lines 55-56). The review's citation `bundle.py:51` is off by three lines. The claim itself is correct. |
| engineering:debug | I traced the V1 path: A3 step 1 builds `lidc-test-shifted`, then the A3 header bundles with no move step. Front 6.1 step 4c is a pre-spend check that says "no test folder" and names the remedy. The leak is caught at $0. |
| engineering:system-design | I checked whether either remedy in step 4c ("moved out" or "a bundle `--view` filter") changes a result. Both keep the test view off the V1 box. Neither changes a score, a hash or the spend. |
| engineering:architecture | I compared this with round 4. Fresh-arch-4 in round 4 graded the same issue minor, for the same reason. The new evidence (A.3.1 arrow label, B.8.1 note gone) is wording drift, not a new failure path. |
| engineering:documentation | I checked every quote against `00_front.md` (6.1 step 4c, line 636; 6.4 rule 7, line 961), `A.md` (A.3.1 line 173; A3 lines 507-518) and `B.md` (B.8.1 has no "bundle"; the review log row 10 cites it). The quotes are accurate. |
| huggingface-skills:huggingface-best | No item is about choosing a model. No change. |
| WebSearch/WebFetch | Not used. The item is about the doc agreeing with itself and with `bundle.py`. |

I graded only the blocking and major items. The review rates items 2-9 minor, so I did not grade them again.

## Verdicts

- fresh-arch-1: DOWNGRADED to minor. The quotes and the code are right (`bundle.py:54` ships every `eval_views/lidc-*`, and A3 and the A.3.1 arrow put `lidc-test-shifted` in the V1 bundle). But front 6.1 step 4c is a $0 check before the rental that requires "no `lidc-test-*` folder" and names the remedy. Either remedy gives the same result (the test view stays off the box), and 6.4 rule 7 already says the view is moved out only for a rental whose table does not name it. So no money is lost and no guess changes the outcome. Round 4 (fresh-arch-4) graded the same issue minor for this reason. The A3/A.3.1 wording and the stale B review-log reference are wording fixes. The review's `--eval-view` filter is a good minor follow-up.
- fresh-arch-2..9: not graded again (the review rates them minor).

Counts: confirmed 0, refuted 0, downgraded 1.
