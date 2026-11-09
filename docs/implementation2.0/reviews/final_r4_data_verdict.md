# Verdict on final_r4_data.md (round 4, data and statistics)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product in scope; no catalog skill applies. Nothing changed. |
| engineering:code-review | Read the argparse code in `train_rl.py:23-25` and `final_eval.py:50` against the B1.5 commands (item 5). |
| engineering:debug | Reproduced the 128-item retest (`r4data/retest.py`, same output). Traced the seed path in `final_eval` and TRL to judge whether the two samplings are independent (item 2). |
| engineering:system-design | Traced `run.sh` plan expansion (`policy.py:471-490`): no flag is injected into a plan step, so missing required flags reach the box (item 5). |
| engineering:architecture | Checked which text owns `b1_rows.yaml` (B1.5, front 6.2) against the test I22 (item 1). |
| engineering:documentation | Checked whether each gap is handled elsewhere (U49 for item 4, E08 for item 5). |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |

Web tools loaded; no web fact was needed for these verdicts.

## Verdicts (blocking and major items)

- data-1: confirmed (blocking). B.md B.9.3 row I22 asserts `--budget-usd 1.90` and a fourth step "R5 on val edge", while B.md B1.5 step 3 says `--budget-usd 2.05` and "R5 on val edge is not a `b1_rows.yaml` step", and 00_front.md 6.2 and B.md line 2767 agree with B1.5; a test and the rule contradict, and `make test` must pass.
- data-2: confirmed (major). B.md B1.9 "Noise model" says "The baseline values stay fixed (they exist)" and calls the one-arm noise "an exact upper bound", but B.md B.6.2 "Paired baseline" makes A3's fresh K = 8 sampling the primary baseline (it does not exist before the rental), so the simulated gate omits the baseline arm; the reproduced retest (var of difference 0.00324 = 7.7 % of the two-arm bound, `r4data/retest.py` on `runs/rl/eval/final/heldout.json` vs `runs/r5/eval/eval_val/items.jsonl`) shows the printed P(pass) that drives the owner's "< 0.5: change a rule" decision is far off. Caveat for the fix: both runs used seed 11 and different code paths (TRL seeds the colocated vLLM engine with `process_index`, `trl/generation/vllm_generation.py:366`), so independence is not proven; the fix must label c = 0.077 as one measured pair and keep c = 1 as the stress line.
- data-3: confirmed (major). B.md B1.2 step 3-4 draw with weight 4 and "Reject a window ... already used in the run", yet claim "about 4p / (1 + 3p) = 14 % of rows ... so B1 changes only the geometry against R5"; weighted draws without replacement deplete the 289 disagreement items (of 7,341, `data/prepared/lidc/train.jsonl`) early, so the share is below 14 % overall and drifts by block (`r4data/disagree.py`), while R5 repeated rows 4x and shuffled (`dataset.py:81`, `grpo_rows` oversample); the engineer cannot meet the stated default and DG14's "expected about 14 %".
- data-4: downgraded (minor). B.md B.6.2 Gate B1 states the far-miss point as a fixed number "<= 20.5 %", and test U49 (B.md B.9.1) asserts exactly "far-miss share <= 20.5 %"; the parenthetical "25 % below R5's 27.3 %" is an explanation, so an engineer does not have to guess; making it relative is a wording/design choice, not a contradiction.
- data-5: confirmed (major). B.md B1.5 step 3 has no `--model` and step 4 has no `--train`, but `train_rl.py:23` has `ap.add_argument("--model", required=True, ...)` and `final_eval.py:50` has `ap.add_argument("--train", type=Path, required=True, ...)`; `run.sh` adds only `-I` (`policy.py:486`), and no Mac test (I22, E08) parses the plan argv, so both steps exit 2 on the paid box.

## Minor items (not re-assessed; stay minor as filed)

- data-6: minor as filed (B1.9 names CR1 where B.6.2 uses the percentile interval; wording).
- data-7: minor as filed (B1.1 pool test for the S0 share belongs to the schedule).
- data-8: minor as filed (third-party duplicate-z patients in val).
- data-9: minor as filed (A3 TOST power sentence).
- data-10: minor as filed (+0.018 decomposition arithmetic).
- data-11: minor as filed (A3 header comment bundle order; front 6.1 check 4c catches it).

Counts: 4 confirmed (1 blocking, 3 major), 0 refuted, 1 downgraded.
