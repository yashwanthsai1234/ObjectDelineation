# Verdict: final round 5, fresh-arch review

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA tool; no catalog lookup or install needed for this verdict. |
| engineering:code-review | Checked each quoted code line myself: policy.py:480 (cap formula), policy.py:61 (INET_CEIL 0.03), policy.py:190-193 (accrued). |
| engineering:debug | Re-ran the money path under the default (b'): money left vs V2 item 1, 1b, 2 costs. |
| engineering:system-design | Judged whether a V2 subset under (b') yields any gate result (none does). |
| engineering:architecture | Checked that decision 1 (b') and decision 3b / E.8 both claim to be the rule for V4. |
| engineering:documentation | Located every source line that repeats "V2 with the money left" or "V4 inside the $3". |
| huggingface-skills:huggingface-best | No model-choice item; no change. |

Only the blocking and major items need verdicts (item 1). Items 2-8 are minor in the review; I did not re-grade them.

## Verdicts

- fresh-arch-1: CONFIRMED (blocking). Front 6.3 "Money left after V1" gives $0.46 / $0.65, and the 6.3 V2 table prices item 1 at $0.55, item 1b at $0.12 and item 2 at $0.15. So under the default decision 1 [(b')] ("V2 runs with the money left", 99_back decision 1; front 6.3 totals) no V2 item that gives a gate result fits ($0.46 < $0.55; $0.55 + $0.15 = $0.70 > $0.65; $0.55 + $0.12 = $0.67 > $0.65). The cap formula `cap = int(min(a.cap_hours * 3600, run["deadline"] - time.time() - MARGIN_S))` [code: experiments/phase2/vast/policy.py:480] ends the box inside item 1 or 2: GPU money spent for nothing. Separately, front 6.3 "V4" ("inside the $3 unless the owner raises the cap"), 99_back decision 3b and E.md E.8 contradict the front 6.3 totals and 99_back decision 1 ("V3 and V4 wait"), and V4 ($0.54) plus V2 cannot both fit $0.46-0.65. The fix must touch 00_front (6.3 heading, totals, V4 paragraph, line 175 order rule, 2.4 diagram "V2 (money left)"), 99_back (decisions 1, 3, 3b), E.md (E.8), A.md (A-G6 and diagram "open until V2 item 3 (money left)"), and D.md (lines 13, 308, 867 "drops it first when the money left is too small").
- fresh-arch-2..8: not re-graded (already minor in the review).

Counts: confirmed 1, refuted 0, downgraded 0.
