# Final evaluation of the RL model, 2026-11-08

Model: the full RL run's final model (HF `svenk029/objdel-sft`, revision `rl-full-3ep`). Script:
`experiments/phase2/final_eval.py`. Box: R5's RTX 5090, right after the shared-prefix test.

Every RL-eligible item of each split, K = 8 samples per item at T = 1, seed 11, scored by the checker. The split
views are the **shifted** crops (object off centre, `build_eval_views`), the harder view the G-probe uses; the test
view was built with `build_eval_views --split test` (seed 0, same rules as the val view). Per item the mean of its 8
samples, then the mean over items; 95 % intervals from 2,000 item bootstrap draws. Illegal answers score 0.

| Metric | Validation (1,084 items) | Test (2,200 items) |
|---|---|---|
| Legal answers | 99.5 % [99.2, 99.7] | 99.4 % [99.1, 99.6] |
| IoU_tc (overlap with the doctors' band) | **0.618** [0.601, 0.635] | **0.645** [0.633, 0.656] |
| IoU_tc, median item | 0.737 | 0.756 |
| Boundary score f5 (a new doctor scores 0.95) | 0.806 [0.789, 0.823] | 0.816 [0.805, 0.829] |
| Precision vs the doctors' majority outline | 0.689 [0.670, 0.708] | 0.695 [0.682, 0.708] |
| Recall vs the majority outline | 0.716 [0.696, 0.735] | 0.748 [0.735, 0.760] |
| Dice vs the majority outline | 0.680 [0.661, 0.697] | 0.701 [0.689, 0.713] |
| Sampling time | 415 s | 634 s |

- Precision and recall are new in `objdel.eval.metrics` (`precision_vote`, `recall_vote`): answer pixels inside the
  strict-majority outline / answer pixels, and majority pixels covered / majority pixels; 0 for failed or illegal
  answers, as the other geometric terms. Test: `tests/test_eval.py::test_precision_and_recall_vs_the_majority_outline`.
- Check against the earlier held-out number: the 128 G-probe items inside this validation run score 0.667 (0.659 in
  `full_rl_2026-11-08.md`, a different sampling order), so the two scripts agree.
- The test split had never been scored before this run. It scores slightly higher than validation on every metric.
- Recall is above precision on both splits: the model's outlines are a little larger than the majority outline.
- Training-time f5 (0.945 at the end of RL) is on training images; on unseen images it is 0.81-0.82.
- Not measured: the SFT checkpoint e2 on the full splits (only on the 128-item probe: IoU_tc 0.511), so the RL gain
  on the whole test split is not measured.

Results: `runs/r5/eval/eval_{val,test}/` (`final_eval.json`, `items.jsonl`, git-ignored).
