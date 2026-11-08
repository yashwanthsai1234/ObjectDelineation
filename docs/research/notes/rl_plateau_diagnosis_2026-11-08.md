# Why RL stopped improving, 2026-11-08

Question: held-out IoU_tc rose 0.511 (e2) -> 0.634 (step 250) -> 0.659 (step 770), test 0.645. Why so little after
the first epoch? Loops (seconds, deterministic): `curve.py` over `runs/rl/train_rl/rl.json` (50-step means) and a
per-item join of `runs/r5/eval/eval_test/items.jsonl` with the test view's rater caches.

## Training curve

| Steps | IoU_tc | F5 | legal | entropy | advantage sd | groups with no spread |
|---|---|---|---|---|---|---|
| 1-50 | 0.641 | 0.833 | 0.968 | 0.271 | 1.04 | 0.0 % |
| 51-100 | 0.733 | 0.920 | 0.996 | 0.127 | 0.61 | 0.0 % |
| 251-300 | 0.753 | 0.934 | 0.996 | 0.099 | 0.46 | 0.4 % |
| 501-550 | 0.759 | 0.934 | 0.986 | 0.068 | 0.46 | 2.1 % |
| 751-770 | 0.771 | 0.943 | 0.996 | 0.062 | 0.34 | 2.7 % |

Slope of training IoU_tc: +0.48 per 100 steps over steps 1-50, then +0.005 to +0.008 per 100 steps. Learning rate
constant (5e-5), gradient norm steady (4-7e-3), median KL 0.27-0.40, PPO clip never active (one update per batch).

## Root cause: RL trains on centred crops, the evaluation uses shifted crops

Object centre's distance from the frame centre (fraction of the side):

| Data | median | share > 0.2 | share > 0.3 |
|---|---|---|---|
| RL training rows (`data/prepared`, `grpo_rows`) | 0.018 | 0.1 % | 0.0 % |
| val, shifted view | 0.259 | 68.1 % | 37.5 % |
| test, shifted view | 0.265 | 69.9 % | 37.5 % |

SFT trained on shifted and rescaled windows (`dataset.py` view rows, `augment.render`); RL's `grpo_rows` reads only
the centred policy frames. Test IoU_tc by offset:

| Offset | n | model | doctor ceiling | misses (IoU < 0.1) |
|---|---|---|---|---|
| < 0.15 | 382 | 0.76-0.77 | 0.71-0.74 | 0-0.5 % |
| 0.15-0.2 | 280 | 0.710 | 0.729 | 1.8 % |
| 0.2-0.3 | 713 | 0.670 | 0.742 | 5.9 % |
| > 0.3 | 825 | 0.545 | 0.724 | 23.0 % |

Doctor ceiling = mean over doctors of the best IoU against the other doctors (leave one out; it ignores topology, and
the model's IoU_tc takes the best of all k doctors, so the ceiling is slightly pessimistic). Whole test split: model
0.645, ceiling 0.732. 10.8 % of items are misses; without them the model scores 0.722 against a ceiling of 0.740.

Size only matters off centre: centred objects of 40-80 px score 0.746 (1.2 % misses); off-centre (> 0.3) 0-40 px
score 0.326 (47 % misses), 40-80 px 0.496 (30 %), > 80 px 0.619 (14 %).

## Reading

1. The training reward is close to saturated: on centred crops the model is at the doctors' level (training IoU_tc
   0.77; centred test items 0.77 against a ceiling of about 0.73). Entropy fell from 0.27 to 0.06 and the advantage
   spread by two thirds, so each step carries less signal. More steps on the same data add little.
2. The held-out gap is a distribution shift that training never shows: off-centre objects, mostly small ones, are
   missed. RL improved them only by transfer.
3. Ruled out: memorisation (the training/held-out gap is already there at step 1, when every training image is new:
   0.641 vs 0.511), a dead signal (groups without spread at most 2.7 %), learning-rate decay (constant), divergence
   (KL and gradient norm steady).

## Fix (built, not yet run on a GPU)

- `train_rl --train-views DIR ...`: RL rows from the centred frames plus windows of the same train items, each scored
  on its own window's rater masks. `dataset.grpo_rows_views` tags the window rows' keys `<i>::<item id>`;
  `serve.loader(*roots)` resolves a tagged key in the i-th root (one root: unchanged).
- Windows: `build_eval_views --split train --out data/views` (2 min 11 s on the Mac), the same draws as the exam views.
  Object offset, median / share > 0.2: train shifted 0.247 / 65.9 %, train rescaled 0.312 / 76.6 % (exams: 0.26 /
  68-70 %). Rows: 8,208 centred + 8,208 shifted + 8,208 rescaled = 24,624.
- Tests: `tests/test_p2_views.py::test_rl_rows_mix_the_centred_frames_with_shifted_windows` (keys, counts, the shifted
  window's caches differ from the centred ones), `::test_one_root_loader_is_unchanged_and_unknown_view_tags_fail`,
  `tests/test_p2_train_rl.py::test_full_rl_trains_on_shifted_windows_beside_the_centred_frames` (TRL's real loop, 2
  steps), `tests/test_profile.py::test_shifted_rl_chain_trains_on_the_windows_and_scores_both_splits`.
- Run plan: `experiments/phase2/configs/rl_shifted.yaml`: 300 steps from the full run's final model with
  `--shared-prefix`, then the whole val and test splits; ~2.6 h, ~$1.6. Needs the owner's budget (Phase 2 has $0.60
  left of the $17.50 hard stop).
