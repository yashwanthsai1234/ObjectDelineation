# Evaluation guide

This guide explains how to score a model, what each metric means, and how the confidence intervals are made.
The metrics run through the same parser, checker and scorer as training ([Reward reference](../reference/reward.md)).
Sources: `experiments/phase2/final_eval.py`, `src/objdel/eval/metrics.py`, `audit.py`, `gates.py`, `signal.py`.

## 1. Score a model on a whole split: `final_eval`

```bash
.venv/bin/python -I -m experiments.phase2.final_eval --ckpt runs/rl/train_rl/final \
    --view data/eval_views/lidc-test-shifted --split test --train data/prepared \
    --out runs/eval/test [--k 8] [--seed 11] [--per-call 64]
```

Run it on a GPU box. It samples with vLLM through the trainer (`g_probe.make_trainer` and `objdel.train.sample`).

| Flag | Default | Meaning |
|---|---|---|
| `--ckpt` | required | A model folder (a merged `final/`, or an SFT branch) |
| `--view` | required | An exam view directory in the `data/prepared` layout |
| `--split` | `val` | The split inside the view directory |
| `--train` | required | The prepared training data (the trainer needs it to start) |
| `--out` | required | Output directory |
| `--model` | `qwen3-vl-2b` | Key of `configs.MODELS` |
| `--k` | 8 | Answers per item |
| `--seed` | 11 | `torch.manual_seed` before sampling |
| `--per-call` | 64 | Items per vLLM generation call |

Steps:

1. Preload expert tables of the training root and the view, and fork the scoring pool before CUDA starts.
2. Take every RL-eligible item of the view's split (all of them, not a subset).
3. Sample K answers per item at temperature 1.
4. Score every answer with `select_ckpt.score_items`, which calls `objdel.eval.metrics`.
5. Per item, take the mean over its K answers. Then take the mean over items, with a 95 % interval.

Outputs:

| File | Content |
|---|---|
| `final_eval.json` | `items`, one block per metric (`mean`, `ci95`), the median of `iou_tc`, timing |
| `items.jsonl` | One line per item: `item_id`, `legal`, `iou_tc`, `f5`, `precision_vote`, `recall_vote`, `dice_vote` |

The script prints one line: `FINAL_EVAL <split> {metric: mean} items N s T`.

### Results of the current model

Source: `docs/research/notes/final_eval_2026-11-08.md`. Model: the full RL run's final model. Shifted views. K = 8, seed 11.

| Metric | Validation (1,084 items) | Test (2,200 items) |
|---|---|---|
| Legal answers | 99.5 % [99.2, 99.7] | 99.4 % [99.1, 99.6] |
| IoU_tc | 0.618 [0.601, 0.635] | 0.645 [0.633, 0.656] |
| IoU_tc, median item | 0.737 | 0.756 |
| Boundary F (F5) | 0.806 [0.789, 0.823] | 0.816 [0.805, 0.829] |
| Precision vs majority outline | 0.689 [0.670, 0.708] | 0.695 [0.682, 0.708] |
| Recall vs majority outline | 0.716 [0.696, 0.735] | 0.748 [0.735, 0.760] |
| Dice vs majority outline | 0.680 [0.661, 0.697] | 0.701 [0.689, 0.713] |

Recall is above precision on both splits: the outlines are a little larger than the majority outline.
Not measured: the SFT checkpoint on the whole splits. So the RL gain on the whole test split is not measured. On the 128-item probe set the gain is +0.148 IoU_tc.

## 2. Metric definitions

All geometric metrics are 0 for an answer that is unparseable, truncated or illegal. This is the reward's legality gate (`legal_gate`).
Per-item values are means over the item's K answers. Reported values are means over items.

| Metric | Key in code | Definition |
|---|---|---|
| **Legal** | `legal` | 1 if the answer parses and the checker finds no defect in any object |
| **IoU_tc** | `iou_tc` | The best IoU against the voting raters who drew the same topology as the answer. If none did, IoU against the strict-majority vote. The reward's own term |
| **Boundary F (F5)** | `f5` | `1 / (1 + (2 D / P) / T5^2)`. Depth of error beyond the band between core and union, normalised by the union contour. A new rater scores 0.95 on average |
| **Precision** | `precision_vote` | Answer pixels inside the strict-majority outline, divided by answer pixels |
| **Recall** | `recall_vote` | Majority pixels covered by the answer, divided by majority pixels |
| **Dice** | `dice_vote` | `2 * IoU / (1 + IoU)` where IoU is against the strict-majority outline |
| IoU vs vote | `iou_vote` | Plain IoU against the majority outline |
| `d_max` | `d_max` | Best IoU against any single rater |
| `g_area` | `g_area` | The area guard of the reward. `min(1, A_max / A) * min(1, A / A_min)` |
| Area ratios | `a_amax`, `a_amin`, `m_e` | A over the largest rater, over the smallest rater, over the union |
| Components | `components`, `extra_components` | Size-filtered parts, and parts beyond the most fragmented rater |
| NSD | `nsd` | Surface Dice at tolerance tau against the best rater: 2 px for nodules, 9 px for disc, 16 px for cup |
| Passed | `passed` | `topology_valid`: legal, a rater drew the topology, and IoU_tc of at least 0.5 for every object |

"Strict-majority outline" means pixels drawn by more than half of the raters (`2 * count > k`).

Item-level metrics over the K answers and the raters (`metrics()`; not in `final_eval.json`):

| Metric | Definition |
|---|---|
| Soft Dice | Soft truth = mean of rater masks. Prediction = mean of the K answer masks. Both binarised at 0.1, 0.3, 0.5, 0.7, 0.9. Dice per threshold, averaged (MRNet protocol) |
| GED | Squared generalised energy distance, `2 E d(S, Y) - E d(S, S') - E d(Y, Y')` with `d = 1 - IoU` |
| HM-IoU | Hungarian matching of the K answers to the raters. Mean IoU of the matched pairs |
| TV, JS | Total variation and Jensen-Shannon (base 2) between the answer topology histogram and the rater histogram, per item |
| Minimum-mode retention | `min over rater answers a of pi_v(a) / f(a)`. 1 when the model's topology histogram equals the raters' |
| Coverage | Share of distinct rater topology answers that at least one answer drew |
| pass@k | Unbiased estimator `1 - C(n - c, k) / C(n, k)` for k in 1, 2, 4, ... up to K (at most 64) |

TV, JS and minimum-mode retention are averaged per image. Histograms are never pooled across images.

## 3. Confidence intervals

`final_eval.summarize(items, draws=2000, seed=0)`:

1. Draw 2,000 resamples of the item list with replacement (`rng.integers(0, n, (draws, n))`).
2. For each metric, take the mean over each resample.
3. Report the 2.5th and 97.5th percentiles of those means as `ci95`.

The unit is the item. Each item contributes the mean of its K answers. So the interval covers item sampling noise. It does not cover the noise of the K draws for a fixed item.
Items are slices. Several slices of one nodule can be correlated, and this interval ignores that.

The signal panel and the G-probe use a stricter unit (`objdel.eval.signal`): 95 % percentile bootstrap, B = 2,000, resampling LIDC **nodules** (the item id up to the slice). Paired differences resample the same nodules for both arms.
For checkpoint comparisons use paired differences on the same items. The paired intervals in `full_rl_2026-11-08.md` (5,000 resamples) are of this kind.

## 4. Checking a result

1. Compare the legal rate first. A low legal rate zeroes every other metric.
2. Compare IoU_tc, then F5. IoU_tc counts overlap with the right topology. F5 counts boundary error against the band.
3. Read precision beside recall. Recall above precision means the outlines are too large.
4. Run the hack audit (section 5) if reward went up and shape quality did not.
5. Never compare numbers across views. The shifted view is harder than the centred view. Test IoU_tc by object offset: 0.76 to 0.77 below 0.15, 0.545 above 0.3 (`docs/research/notes/rl_plateau_diagnosis_2026-11-08.md`).

## 5. Hack audit (`objdel.eval.audit`)

`hack_audit(items, baseline=None)` reads per-item metrics and returns statistics and alarms per structure. `alarms_for(records, caches, baseline)` does the same in training for the `Controller`.

| Alarm or flag | Rule |
|---|---|
| `area_median` | Median `A / A_max` above 1.15 |
| `area_tail` | More than 10 % of answers with area above 1.3 `A_max`. With a baseline, more than 10 points above its share |
| `components` | Median component count at least 0.5 above the baseline's. Needs a baseline |
| `duplicates` | Duplicate share inside groups above 50 % |
| `inflation`, `deflation`, `fragmentation` | More than half the answers lie above, below, or more fragmented than the rater envelope `[0.99 A_min, A_max / 0.99]`. Or the share rises by 10 points over the baseline |
| `newly_rewarded_failures` | Failure share among answers with positive advantage exceeds the failure share overall (Yao et al. 2026) |
| `topology_collapse` | Mean minimum-mode retention on disagreement items below 0.5. Measured only with K of 64 or more |

A baseline is an audit of the SFT checkpoint. A check fires only when it fires now and not on the baseline.

## 6. Entry gates and probes

`objdel.eval.gates.phase2_entry(report, datasets=None, prior_a_threshold=None)` decides the Phase 2 entry gates from one report.

| Gate | Requirement |
|---|---|
| G-format | Parse at least 98 %, legal at least 90 % (shifted view) |
| G-geometry | On the shifted and rescaled views: median IoU_tc at least 0.60 for LIDC (RIGA: 0.85 disc, 0.75 cup), illegal as 0. At least 70 % of answers legal with IoU_tc of at least 0.5 for every object |
| G-prior-b | Paired mean gain over an image-blind prior, with a 95 % lower bound above 0 |
| G-precision | Synthetic precision probe: legal at least 95 %, median ASSD at most 0.6 px, IoU at least 0.80, centroid error at most 0.5 px, and a win over the centroid disc on at least 60 % of items |
| G-topology-a | Simple controls answered `(1, 0)` at least 95 % |
| G-signal | Duplicates at most 10 %, sigma share at least 50 % at G = 8 |
| G-vision | Projector ratio at least 0.5 and both TVs at least 0.002 |

Report-only: G-prior-a, topology-probe hole recall (against 0.30), truncation and token statistics, G-probe verdict.
The thresholds are fixed in `experiments/phase2/prereg.py`. A rule is never changed after it is used.
`experiments/phase2/RESULTS.md` section 1 holds the first full gate table (rung 1, LFM). The full gate table of the Qwen3-VL-2B rung is not recorded there.

## 7. Other entry points

| Tool | Use |
|---|---|
| `experiments.phase1.eval_run` | Evaluate sampled answers (JSONL `{"item_id", "samples", "truncated"}`) against the store. Writes JSON and Markdown |
| `experiments.phase2.select_ckpt` | Branch selection and the full gate report |
| `experiments.phase2.probes build` / `score` | Build and score the synthetic precision and topology probes |
| `experiments.phase2.baselines` | The image-blind prior, MedSAM with an exact or jittered box, SAM 2.1-tiny |
| `experiments.phase2.hf_rate`, `r0_parity` | HF sampling rate, vLLM versus HF parity |
