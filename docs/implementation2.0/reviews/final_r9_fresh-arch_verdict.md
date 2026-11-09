# Verdict: final round 9, fresh-arch review

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA tool. I did no catalog lookup and installed nothing. |
| engineering:code-review | I read `experiments/phase2/vast/policy.py:422-435` (`chain`, the `steps:` branch with `str(x).format(rd=rd, data=data)`). This confirmed item 1. |
| engineering:debug | I copied the `a3_measure.yaml` block from A.md (A3, lines 492-497) unchanged and ran `policy.chain(Path(file), 'r')` with the repo venv. It raised `KeyError(' bash experiments/cuda_checker/run_box')` [measured 2026-11-09; file: reviews/r9v_a3.yaml]. |
| engineering:system-design | For item 3 I traced what the default (b') leaves on `main`: E4 step 1, E.8, front 6.3 "V4" and 99_back decision 3b all say the old samplers stay until V4, and the owner gets that choice in decision 3b. |
| engineering:architecture | Item 3: I weighed keeping an untested-path fallback against the archive rule. The doc states the trade-off and names who decides, so it is not a hidden gap. |
| engineering:documentation | I checked every quote against its source file and line. |
| huggingface-skills:huggingface-best | No item is about choosing a model. No change. |
| WebSearch/WebFetch (loaded) | Every item is about the doc against the repo, so I made no web query. |

I graded only the blocking and major items (1-3). The review already rates items 4-9 as minor, so I did not grade them again.

## Verdicts

- fresh-arch-1: CONFIRMED (blocking). A.md A3 "Plan file `a3_measure.yaml`" step 2 holds a bash group `{ bash ... ; }`, and `chain` runs `str(x).format(rd=rd, data=data)` on every `steps:` entry [code: policy.py:434-435]. Running `policy.chain` on the doc's exact block raises `KeyError(' bash experiments/cuda_checker/run_box')` [measured 2026-11-09]. So `run.sh <run> a3_measure.yaml` (V1 items 4-5, "Not skippable") fails on the Mac. That kills the Gate B1 paired baseline and Gate A, and V2 item 1b fails the same way. No other plan in the doc has `{` inside a step (grep of phases/*.md). Fix in A: use a subshell `( ... )` or `{{ ... ; }}`, and make I20 call `policy.chain` on both A3 plans. The B guard (I47 or a new test that calls `chain` on every `steps:` plan named in the launch lines) costs $0.
- fresh-arch-2: DOWNGRADED to minor. The stale "CR1" word in B1.9 contradicts B.6.2 (percentile patient interval for Gate B1, 78 patients). But B1.9's own earlier bullet says "the A1 decision interval", and the same bullet says "the same as the gate's", so the rule in force is clear. The choice also barely changes the result. On R5 val (1,084 items, 78 patients), the CR1 interval is [0.5837, 0.6527] (width 0.0690). The percentile patient bootstrap (10,000 draws) is [0.5846, 0.6527] (width 0.0682) [measured 2026-11-09: runs/r5/eval/eval_val/items.jsonl, patient = series uid]. A 1 % width difference cannot move joint P(pass) enough to matter. This is a wording fix.
- fresh-arch-3: DOWNGRADED to minor. The doc says the same thing everywhere it names the trigger. E4 step 1: "Until then they stay". Front 6.3: "Until V4 passes, E4 keeps the old samplers". 99_back decision 3b gives the owner that exact consequence with "[no default: owner decision]". So the engineer has nothing to guess: under (b') the two scripts stay, and the owner knows why (a box-tested vLLM fallback). As modules run as `__main__`, they also do not trip CLEAN-G9 (`make dead`). Three points remain and are wording fixes: the front's "no old code sits beside new code" sentence is too absolute, E.8's "E4 is done" disagrees with the E4 gate's "not done until E-G10's box part passes", and the end state under (b') could say "until V4 or never". The preferred "fix" would archive the only box-tested sampler before the new one has run on a box, which is a trade-off for the owner, not a defect.
- fresh-arch-4..9: not graded again (the review rates them minor).

Counts: confirmed 1, refuted 0, downgraded 2.
