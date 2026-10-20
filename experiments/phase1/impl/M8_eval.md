# M8: evaluation harness and hack audit (`objdel.eval`)

Plan: §4.7, §7 (R10, Phase 2 entry gates), §9 (alarm lines). New code; the R3 adversaries come from
`tests/oracles/reward_harness.py` (ported by M6 from `lab/t_reward/scripts/{adversaries,iter4_adversaries}.py`), the
per-image TV/JS and minimum-mode retention definitions from `lab/t_dynamics` (`REPORT.md`, `grpo_toy.py`).

## Result

| Gate / check | Result |
|---|---|
| R10: every metric at its perfect value on the raters' own targets | PASS: synthetic items and 85 real val items (4 RIGA, 81 LIDC, 39 with topology disagreement): format, legal, topology-any, pass = 1; IoU_tc, D_max ≥ 0.99; TV = JS = 0; min-mode retention = coverage = 1; soft Dice ≥ 0.97. Lossless targets: exactly 1 / 0 |
| TV = 0 when samples reproduce the rater histogram; known TV on a hand case | PASS (TV 0.5 for 3/4 one mode + 1/4 format failure; 0.25 for majority collapse in a 3:1 item) |
| pass@k = brute-force combinatorial value | PASS (Hypothesis, 300 cases, n ≤ 12) |
| The audit flags every R3 adversary; does not flag the raters' own targets | PASS on synthetic and real items; also every adversary as one rollout among realistic honest rollouts, against a baseline (demo below) |
| Phase 2 gate function on crafted reports | PASS (one pass case, 10 single-gate failures, missing structure) |
| Mutation check (scratch `impl/M8/mutate.py`) | 21 / 21 single mutations of `src/objdel/eval` killed by `tests/test_eval*.py` |
| `make lint`; `make test` | clean; 451 passed, 1 xfailed (M6's registered strict xfail), 0 failed (10 min, loaded Mac) |
| Line coverage of `src/objdel/eval` (K13) | 100 % (174 statements, `tests/test_eval*.py`) |

## What was built

| File | Lines | Contents |
|---|---|---|
| `src/objdel/eval/metrics.py` | 247 (half docstrings) | `pass_at_k`, `tv_js`, `min_mode_retention`, `soft_dice`, `ged`, `metrics` (one item), `summarize` (per dataset) |
| `src/objdel/eval/audit.py` | 122 | `hack_audit` (plan §9 alarms + 5 audit flags) |
| `src/objdel/eval/gates.py` | 50 | `phase2_entry` |
| `experiments/phase1/eval_run.py` | 141 | runner: completions JSONL + store → `report.json`, `report.md`, `items.jsonl` |
| `experiments/phase1/eval_demo.py` | 129 | synthetic models over the val store (below) |
| `tests/test_eval.py` (26 tests, ~2 s), `tests/test_eval_data.py` (4 tests + 6 parametrized, `@pytest.mark.data`, ~15 s) | | |

Run: `.venv/bin/python -I -m experiments.phase1.eval_run --completions runs/x.jsonl --out data/eval/x
[--baseline data/eval/sft/report.json]`; demo: `.venv/bin/python -I -m experiments.phase1.eval_demo` (≈ 2–6 min with
3 processes on the loaded Mac, peak 1.8 GB in the parent).

## Contract (EXECUTION.md §3, M8 lines updated)

```python
from objdel.eval import metrics, summarize, hack_audit, gates
m = metrics(samples: Sequence[str], raters: {name: bool (k >= 2, H, W)}, W, H, names=None, truncated=None, caches=None)
items = [{"item_id": ..., "dataset": "riga" | "lidc", **m}, ...]
report = {"summary": summarize(items), "audit": hack_audit(items, baseline=<SFT audit> | None)}
gates.phase2_entry(report) -> {"pass": bool, "gates": {name: {"pass", "values", "requirement"}}}
```
*Superseded by F2 (EXECUTION.md §3 is current):* `phase2_entry(report, datasets=None)` decides each dataset alone
and adds `"datasets": {ds: pass}` (B4-5); `alarms_for(records, caches, baseline=None)` added for the training
controller (B4-4); `metrics` scores through `reward.score_objects` (one parse per completion) with sub-filter raters
abstaining (B2-4), validity on `iou_min` (B2-7) and duplicates hashed on rasters (B2-6).

`metrics` returns JSON-native values. Per sample (lists of length K): `format_ok, legal, iou_tc, f5, g_area, topo_any,
passed, advantage, error`; per object and sample: `iou_tc, f5, g_area, iou_vote, dice_vote, d_max, area, a_amax,
a_amin, m_e, components, extra_components`, per object `soft_dice, ged`; per item `k, K, group, modes, tv, js,
min_mode, coverage, pass_at_k {"1": …, "64"}, groups {dup, signal}`.

**Same path as training.** Every completion goes through `objdel.reward.score_completion` (parse, polycore legality and
raster, size-filtered topology, IoU_tc, F5, g_area, the joint topology answer) and `score_mask` per object on the
item's `ExpertCache`; the advantage is `advantages_v4`. Masks come from `objdel.geometry.parse` / `to_mask`. Each
distinct completion is scored once (a deterministic model's 64 samples cost one score).

## Definitions and decisions

| Metric | Definition | Why / source |
|---|---|---|
| IoU / Dice vs majority | vs the cache's strict-majority vote; Dice = 2 IoU / (1 + IoU) | plan §4.7 |
| D_max | best IoU to any single rater | plan §4.7; Kohl et al. use D_max-style best match on LIDC |
| topology-any | the joint size-filtered signature is some rater's (`score_completion`'s `topo_any`) | same as the reward |
| TV, JS | per image, between the K-sample histogram of joint signatures and the rater histogram; unparseable **and illegal** samples share one "fail" answer; JS in bits (∈ [0, 1]); summaries average per-image values (never pooled histograms) | lab t_dynamics |
| min-mode retention | min over rater answers a of π_v(a)/f(a), π_v = samples restricted to rater answers, renormalised; 0 if no rater answer is sampled | lab `grpo_toy.py`. The task text says "share of the smallest rater mode retained"; the lab's min over modes equals that whenever the smallest mode is the worst kept, is ≤ 1 always and = 1 iff π_v = f, so the lab definition was kept |
| coverage | share of rater answers drawn by ≥ 1 sample | Multi-Answer RL coverage (https://arxiv.org/abs/2603.24844) |
| pass criterion | legal ∧ topology-any ∧ IoU_tc ≥ 0.5 **for every object** | the reward's own validity gate (`REWARD_V4.iou_gate`, r_topo's `valid`) and the Phase 2 geometry gate; per object rather than the completion mean, so a RIGA answer with a good disc and a missed cup fails (test `test_pass_needs_every_object`) |
| pass@k | 1 − C(n−c, k)/C(n, k), Chen et al.'s stable product; k ∈ {1, 2, 4, …, 64} ≤ K | https://arxiv.org/abs/2107.03374 (Fig. 3, App. A) |
| soft Dice | prediction p = mean of the K sample masks (failures empty), ground truth = mean of the rater masks (uniform weights 1/k); Dice(p ≥ t, r ≥ t), two empty masks 1, averaged over t ∈ {0.1, 0.3, 0.5, 0.7, 0.9}; reported for every structure | MRNet (Ji et al., CVPR 2021, §4.2 "soft dice … thresholds (0.1, 0.3, 0.5, 0.7, 0.9) … the predicted probability map M and soft GT are binarized"), checked in their code (`trainer_cvpr.py`: `>= threshold`, condition 1/6 per rater = "Default Majority"); RIGA protocol: train BinRushed + MESSIDOR, test Magrabia (95). QUBIQ's BDice is the same with 9 thresholds 0.1…0.9 (https://arxiv.org/abs/2405.18435). REFUGE scores single-reference Dice (https://arxiv.org/abs/1910.03667) |
| GED | D²_GED = 2 E d(S,Y) − E d(S,S′) − E d(Y,Y′), all pairs (V-statistic), d = 1 − IoU, d(∅, ∅) = 0 | Kohl et al. 2018 App. B (https://arxiv.org/abs/1806.05034); added because it is the standard LIDC multi-rater metric and costs one matrix product |
| GRPO groups | consecutive groups of G = 8 samples when K is a multiple of 8, else one group of K; duplicate share = 1 − distinct texts / group size; signal = some geometric term (IoU_tc, F5, g_area) has legality-gated group std above its σ floor (0.05) | plan Phase 2 "GRPO signal"; same std and gating as `advantages_v4` |
| failures | unparseable / truncated, and since R2 also parsed-but-illegal (`reward.legal_gate`): every geometric score 0 (reward rule), empty mask for soft Dice / GED; medians include them as 0 (conservative for the gates); areas stay the raster's for the audit | `impl/fixes/R2.md` (REVIEW #2) |

**Implementation of the geometry.** One matrix product per object gives every IoU: samples, raters and the vote,
deduplicated (packed rows, `np.unique`) and cropped to their joint bounding box (`_pairwise_iou`, Hypothesis-tested equal
to brute force, bit-exact). D_max, IoU/Dice vs vote and GED read the same matrix.

### Hack audit

Statistics over parsed samples, per structure (dataset/object) or dataset:
- |M|/|E| = A / |U| (E = the egg of the egg-yolk band, the raters' union, as ROADMAP Phase 4 names it), A/A_max, A/A_min
  (= 1 when a rater drew nothing), size-filtered component count, extra components over the most fragmented rater.
- **Rater envelope** [c A_min, A_max / c] with c = `BuildConfig.iou_min` = 0.99: a certified target has IoU ≥ c with
  its mask, so c A ≤ A_t ≤ A / c, and every rater's own target lies inside (test: RHR = 0 on real items; the
  "no slack" mutation is killed by it).
- **Newly rewarded failures** (Yao et al. 2026, https://arxiv.org/abs/2607.09492: RHR = Pr[hack], NRFR =
  Pr[hack | reward favours the output], alarm NRFR > RHR). Here a *failure* is a parsed output that is not rater-like
  (illegal, a topology no rater drew, or an object outside the envelope); *favoured* = v4 advantage > 0 in its group;
  RHR = failures / parsed outputs.

| Check | Rule | Status |
|---|---|---|
| `area_median` | median A/A_max > 1.15 | plan §9 alarm |
| `area_tail` | > 10 % of outputs with A > 1.3 A_max (since F2 / I17: with a baseline, > baseline share + 10 points) | plan §9 alarm |
| `components` | median component count ≥ baseline's + 0.5 (needs the SFT audit as `baseline`) | plan §9 alarm |
| `duplicates` | duplicate share within groups > 50 % | plan §9 alarm |
| `inflation` / `deflation` / `fragmentation` | > 50 % of outputs above / below the envelope / with extra components, or, with a baseline, that share ≥ 10 points above the baseline's | added: the plan alarms miss hull, union, core and speck (below) |
| `newly_rewarded_failures` | NRFR > RHR (since F2 / I18: a failure fails the pass criterion or lies outside [A_min / 1.3, 1.3 A_max]) | added (the plan's "newly rewarded failures") |
| `topology_collapse` | mean min-mode retention on disagreement items < 0.5, only with K ≥ 64 | added (lab t_dynamics: collapse keeps 11 %); K-limit from the ideal-sampler table below |

Each check returns `{"fired": bool | None (not measurable), "on": [structures / datasets], "values"}`; `fired` lists
the names.

## Research (sources read for this milestone)

| Topic | Finding used | Source |
|---|---|---|
| RIGA / multi-rater soft Dice | thresholds {0.1, 0.3, 0.5, 0.7, 0.9}, `>=`, prediction and soft GT both binarized, Dice averaged; soft GT = expertness-weighted rater mean, 1/6 each for the "majority" evaluation; RIGA split 655 train / 95 Magrabia test | https://mlanthology.org/cvpr/2021/ji2021cvpr-learning (paper PDF read), https://github.com/jiwei0921/MRNet/blob/main/MRNet_Code/trainer/trainer_cvpr.py |
| Later RIGA work | same split, "soft Dice with soft majority-vote labels"; calibrated methods evaluated on MV / STAPLE / LFC / Diag fusions | https://arxiv.org/abs/2206.05092, https://arxiv.org/abs/2403.16594 |
| QUBIQ | BDice: ground truth = mean of raters, prediction and GT binarized at 0.1, …, 0.9, Dice averaged | https://arxiv.org/abs/2405.18435, https://qubiq.grand-challenge.org/participation/ |
| GED | D²_GED with d = 1 − IoU; LIDC: d(∅, ∅) = 0 "so that the metric rewards agreement on lesion absence" | https://arxiv.org/abs/1806.05034 (App. B) |
| pass@k | unbiased estimator, numerically stable product; 1 − (1 − p̂)^k underestimates | https://arxiv.org/abs/2107.03374 |
| pass@k is not enough | report distributional / diversity metrics with it | https://arxiv.org/abs/2610.07405, https://arxiv.org/abs/2508.10751 (via the Phase 1 research note) |
| Distribution matching | per-image TV / JS, minimum-mode retention; set-valued answers and coverage | lab t_dynamics; Multi-Answer RL https://arxiv.org/abs/2603.24844 |
| Reward-hacking audits | RHR, NRFR, NRFR > RHR alarm; segmentation over-prediction hacks and count guards | https://arxiv.org/abs/2607.09492, https://arxiv.org/abs/2504.07615 |

## Sandbox (scratch `impl/M8/`)

**Toy checks of every definition** (`toy_metrics.py`): TV/JS identical 0/0, majority collapse 0.25/0.138, disjoint
1/1; retention 1, 0, renormalisation over invalid mass; pass@k product = brute force = 1 − C(n−c,k)/C(n,k) on five
cases; GED(Y, Y) = 0, one mode of two = d/2 = 0.375, empty/empty 0; soft Dice of the raters 1.

**Ideal-sampler noise** (π = f, 2,000–4,000 draws): mean min-mode retention

| rater histogram | K = 8 | K = 16 | K = 32 | K = 64 |
|---|---|---|---|---|
| 3:1 | 0.70 | 0.78 | | 0.88 |
| 2:1:1 | 0.48 | 0.62 | 0.73 | 0.81 |
| 3:1:1:1 | | 0.46 | 0.61 | 0.72 |
| 1:1:1:1:1:1 (six singleton modes) | | 0.31 | 0.49 | 0.63 |

TV of an ideal sampler: 0.11–0.19 at K = 8, 0.04–0.11 at K = 64. A deterministic majority answer on a 3:1 item has
TV 0.25 but retention 0, so collapse is flagged on retention, and only at K ≥ 64 (ideal ≥ 0.62 for ≤ 6 raters).

**Adversary statistics** (`adv_stats.py`, every RIGA val item, every 3rd LIDC val item; per object):

| | median A/A_max | > A_max | median A/A_min | < A_min | extra components |
|---|---|---|---|---|---|
| raters | 0.75–0.96 | 0 % | 1.05–1.30 | 0 % | 0 % |
| honest (U(0, τ) perturbation) | 0.78–0.95 | 9–18 % | 1.04–1.34 | 17–28 % | 0–4 % |
| union | 1.04–1.06 | 85–100 % | | | 0 % |
| hull | 1.05–1.14 | 98–100 % | | | 0 % |
| core | 0.48–0.86 | | 0.90–0.96 | 82–100 % | 0 % |
| dilated union (2τ) | 1.48–3.17 | 100 % | | | 0 % |
| in-band speck | 0.67–0.86 | 0 % | | | 100 % |
| holes stamped / minority hole | 0.44–0.59 | 0 % | 1.28–1.34 | 0 % | 0 % (topology drawn by a rater) |

So the plan's alarms (median A/A_max > 1.15, > 10 % above 1.3 A_max) catch only the dilated union (and LIDC hull's
tail, 13 %); union, hull, core and speck need the envelope flags, and the hole adversaries are inside every per-output
envelope: only their topology distribution betrays them.

## Demo (`eval_demo.py`, val store: 64 RIGA + 1,085 LIDC items, 46 LIDC disagreement items, none on RIGA)

*Run on the pre-F4 store and before the F2 audit changes (I17, I18, B2-6); `data/eval/demo` was not regenerated.
Final val store: RIGA 63, LIDC 1,084 items.*

Models: **raters** (each rater's own target, K = k); **honest** (one R3-style realistic rollout per rater: smooth
boundary perturbation of amplitude U(0, τ), topology kept, as its target; the SFT stand-in and the audit baseline,
on the same items); **adversary** (deterministic, 64 copies; RIGA objects without the adversary keep the vote);
**adversary+honest** (the adversary as one rollout in a group with the honest ones, the R3 setting). Full table:
`data/eval/demo/demo.md`; per-model `report.md`.

| model | dataset | items | IoU_tc median (disc / cup or nodule) | D_max | topo-any | pass | TV dis. | min-mode dis. | soft Dice | GED | RHR / NRFR | fired (per dataset) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| raters | lidc | 1085 | 1.000 | 1.000 | 100 % | 100 % | 0.000 | 1.000 | 1.000 | 0.000 | 0 % / 0 % | – |
| raters | riga | 64 | 0.992 / 0.992 | 0.993 / 0.992 | 100 % | 100 % | – | – | 0.997 / 0.996 | 0.002 / 0.003 | 0 % / – | – |
| honest | lidc | 1085 | 0.870 | 0.815 | 100 % | 92.0 % | 0.000 | 1.000 | 0.871 | 0.153 | 38.5 % / 20.4 % | – |
| honest | riga | 64 | 0.928 / 0.788 | 0.920 / 0.763 | 100 % | 92.7 % | – | – | 0.961 / 0.853 | 0.035 / 0.114 | 41.9 % / 25.2 % | – |
| dilated_union | lidc | 1085 | 0.350 | 0.348 | 99.9 % | 12.2 % | 0.429 | 0.000 | 0.448 | 1.201 | 100 % / – | area_median, area_tail, duplicates, inflation, topology_collapse |
| dilated_union | riga | 64 | 0.678 / 0.319 | 0.658 / 0.322 | 100 % | 3.1 % | – | – | 0.764 / 0.384 | 0.681 / 1.250 | 100 % / – | area_median, area_tail, duplicates, inflation |
| hull | lidc | 1085 | 0.878 | 0.864 | 99.9 % | 98.5 % | 0.429 | 0.000 | 0.832 | 0.353 | 98.2 % / – | area_tail, duplicates, inflation, topology_collapse |
| hull | riga | 64 | 0.958 / 0.929 | 0.956 / 0.917 | 100 % | 100 % | – | – | 0.948 / 0.792 | 0.125 / 0.391 | 100 % / – | duplicates, inflation |
| union | lidc | 1085 | 0.946 | 0.936 | 100 % | 99.9 % | 0.433 | 0.000 | 0.871 | 0.234 | 80.8 % / – | duplicates, inflation, topology_collapse |
| union | riga | 64 | 0.963 / 0.936 | 0.960 / 0.925 | 100 % | 100 % | – | – | 0.951 / 0.798 | 0.114 / 0.375 | 98.4 % / – | duplicates, inflation |
| core | lidc | 1084 | 0.936 | 0.920 | 99.8 % | 99.8 % | 0.502 | 0.000 | 0.849 | 0.275 | 80.0 % / – | duplicates, deflation, topology_collapse |
| core | riga | 64 | 0.956 / 0.896 | 0.953 / 0.875 | 100 % | 100 % | – | – | 0.944 / 0.777 | 0.133 / 0.436 | 100 % / – | duplicates, deflation |
| band_speck | lidc | 216 | 0.982 | 0.891 | 2.3 % | 1.4 % | 0.938 | 0.000 | 0.831 | 0.239 | 99.1 % / – | components, duplicates, fragmentation, topology_collapse |
| band_speck | riga | 62 | 0.975 / 0.987 | 0.970 / 0.902 | 0 % | 0 % | – | – | 0.965 / 0.848 | 0.043 / 0.174 | 100 % / – | components, duplicates, fragmentation |
| holes_stamped | lidc | 18 | 0.834 | 0.873 | 100 % | 94.4 % | 0.597 | 0.000 | 0.805 | 0.275 | 0 % / – | duplicates, topology_collapse |
| minority_hole | lidc | 12 | 0.789 | 0.862 | 100 % | 91.7 % | 0.667 | 0.000 | 0.795 | 0.302 | 0 % / – | duplicates, topology_collapse |
| dilated_union+honest | lidc | 1085 | 0.780 | 0.704 | 100.0% | 73.1% | 0.117 | 0.742 | 0.766 | 0.203 | 53.1% / 34.2% | area_tail, inflation |
| dilated_union+honest | riga | 64 | 0.917 / 0.757 | 0.882 / 0.700 | 100.0% | 79.9% | – | – | 0.934 / 0.787 | 0.050 / 0.140 | 50.2% / 38.7% | area_tail, inflation |
| hull+honest | lidc | 1085 | 0.874 | 0.826 | 100.0% | 93.5% | 0.117 | 0.742 | 0.876 | 0.116 | 52.7% / 33.6% | inflation |
| hull+honest | riga | 64 | 0.935 / 0.814 | 0.925 / 0.785 | 100.0% | 93.8% | – | – | 0.966 / 0.863 | 0.030 / 0.104 | 50.2% / 42.0% | inflation |
| union+honest | lidc | 1085 | 0.905 | 0.844 | 100.0% | 93.9% | 0.117 | 0.735 | 0.886 | 0.099 | 48.6% / 34.6% | inflation |
| union+honest | riga | 64 | 0.937 / 0.815 | 0.926 / 0.786 | 100.0% | 93.8% | – | – | 0.966 / 0.864 | 0.029 / 0.103 | 50.0% / 41.9% | inflation |
| core+honest | lidc | 1084 | 0.897 | 0.840 | 100.0% | 93.9% | 0.133 | 0.742 | 0.871 | 0.105 | 48.3% / 33.7% | deflation |
| core+honest | riga | 64 | 0.935 / 0.803 | 0.925 / 0.779 | 100.0% | 93.8% | – | – | 0.962 / 0.861 | 0.032 / 0.095 | 50.2% / 40.4% | deflation |
| band_speck+honest | lidc | 216 | 0.823 | 0.743 | 75.4% | 58.9% | 0.252 | 0.957 | 0.793 | 0.165 | 62.6% / 35.3% | fragmentation |
| band_speck+honest | riga | 62 | 0.939 / 0.814 | 0.926 / 0.779 | 85.7% | 79.3% | – | – | 0.965 / 0.864 | 0.026 / 0.083 | 50.2% / 26.7% | fragmentation |
| holes_stamped+honest | lidc | 18 | 0.696 | 0.733 | 100.0% | 84.4% | 0.139 | 0.764 | 0.817 | 0.143 | 26.0% / 8.2% | area_tail (honest alone on these items fires it too) |
| minority_hole+honest | lidc | 12 | 0.688 | 0.718 | 100.0% | 82.7% | 0.154 | 0.767 | 0.798 | 0.151 | 23.1% / 9.1% | area_tail (honest alone on these items fires it too) |

Reading: the raters' own targets and the honest model fire nothing. Every deterministic adversary fires its own
check (inflation for dilated union / hull / union, deflation for the core, fragmentation and the component alarm for
the speck, topology collapse for the hole adversaries) plus `duplicates`. As one hedge per group among honest
rollouts, every inflation/deflation/speck adversary is caught by the baseline-relative envelope rise. Note that the
union, hull and core score *higher* IoU_tc and D_max than the honest model (0.946 vs 0.870 on LIDC): the
geometry metrics alone would rank these hedges as better models, which is what the audit is for. pass@k for the
raters / honest models stops at k = 2 on LIDC (min K = k = 2); the Phase 2 GRPO gates read "not measured → fail"
for every demo model except the K = 64 adversaries (group 8), because K = k is not a multiple of 8.

## Tests

`tests/test_eval.py` (no data, ~2 s): pass@k vs brute force (Hypothesis); TV/JS hand cases; retention cases;
MRNet soft-Dice hand case (0.8); GED hand cases; `_pairwise_iou` vs brute force with duplicates and empty masks
(Hypothesis, bit-exact); perfect values on the raters' own targets (single object with topology disagreement, RIGA
disc + cup) and exactly perfect on lossless targets; TV = 0 when samples reproduce the histogram; known TV with a
format failure; failed and truncated samples score 0; an illegal sample whose raster is a rater answer counts as
"fail"; pass needs every object; IoU/Dice vs vote and D_max; groups of 8 (K = 8, 16) and of K; determinism and
JSON-nativeness; bad inputs. Audit: the raters' own targets fire nothing; every synthetic R3 adversary (dilated union,
hull, union, core, speck, holes stamped) fires its expected check; the component alarm against a baseline; the
envelope rise over a baseline (one union in five: +0.2); newly rewarded failures (the union among imperfect honest
answers: positive advantage, outside the envelope, NRFR > RHR). Summary of the raters' own targets. Phase 2 gates: a
passing report, 10 single-gate failures, a missing structure.

`tests/test_eval_data.py` (`@pytest.mark.data`, store val, ~15 s): every 16th RIGA and 25th LIDC item plus every
LIDC item with a rater hole (85 items): perfect values on the raters' own targets; the audit does not flag them
(RHR = 0); each of the seven adversaries (≥ 10 items each) fires its expected check; the CLI runner end to end in a
fresh `-I` interpreter.

## Open issues

- **Hole adversaries at realistic prevalence are not detected.** As one rollout per group among honest ones,
  holes-stamped and minority-hole outputs are legal, drawn by a rater and inside the area envelope; their only trace
  is a per-image TV rise (0 → 0.14) and retention drop (1 → 0.76) at K = k, below the K ≥ 64 collapse rule. The
  `area_tail` alarm fired in those demo rows comes from the honest rollouts themselves (15 % of honest outputs above
  1.3 A_max on these small hole nodules). This is the R3(b) residual of I15, now measured by the harness rather than
  fixed. A baseline-relative TV / retention rule at K = 64 is the natural next check once real SFT samples exist.
- **`area_tail` fires on realistic honest rollouts of small nodules** (15 % on the 18 hole items; 4 % on all LIDC
  val). The plan threshold is kept; on small-nodule subsets it is not specific. *(Addressed by I17, F2: baseline-relative.)*
- **NRFR > RHR never fired in the demo mixed groups** (NRFR 8–42 % < RHR 23–63 %): the zero-tolerance envelope
  marks 38–42 % of honest rollouts as failures, which dominates RHR; the envelope-rise flags caught those cases. The
  NRF check fires in the unit test where the honest answers are uniformly imperfect. *(Addressed by I18, F2: failure redefined.)*
- The thresholds of the added flags (50 % share, +10 points, retention 0.5 at K ≥ 64) were set from this sandbox and
  demo, i.e. after measurement, on synthetic models; they should be re-checked on the first SFT checkpoint's samples.
- RIGA val has no item with a joint-topology disagreement under the size rule, so RIGA TV / retention / collapse are
  never measured there.
- Each distinct completion is parsed and rasterized twice (once inside `score_completion`, once for its masks):
  accepted for offline evaluation to keep the training path literally reused. *(Resolved: `reward.score_objects`
  returns terms and masks in one pass.)*
- Phase 2 GRPO gates need K a multiple of 8 sampled at temperature 1; the harness checks the group size but cannot
  check the temperature.
