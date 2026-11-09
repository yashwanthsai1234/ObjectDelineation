# Final review, round 2: data quality and statistics (Implementation.md)

Reviewer role: medical imaging data scientist and statistician (LIDC, multi-rater, leakage, bootstrap).
Scope: the whole doc, with A and B1 first. I checked each r1 data item (`final_r1_data.md`, 27 items) against the
revised text. Items 1-27 are fixed or answered, so they are not repeated here. Exceptions: data 21 and data 26 are
only partly fixed (see items 3, 16 and 17).
Line numbers refer to `docs/implementation2.0/Implementation.md` as of 2026-11-08.
Counts: **1 blocking, 6 major, 11 minor (18 items).**

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers data audit, bootstrap or LIDC QA. The doc's existing recommendation (`dicom-series-preflight` for B1.7) still fits. Item 8 shows that local sidecar files answer part of B1.7 at $0, before any DICOM download. Nothing recommended or installed. |
| engineering:code-review | Checked claims against code: `augment.shift/draw/neighbour_side/row_rng`, `build_eval_views._record`, `frame` padding, `final_eval` (seed 11, `prompts()` = RL items only), `splits.deal/lidc_stratum`, `lidc.slice_index/ct_slices`, `fetch_lidc_pixels.fetch` sidecars, `git log` of the render path. |
| engineering:debug | Reproduced the numbers before judging. The B.1 table matches exactly (val 0.6182 / 12.64 %, far 406 / 0.503 / 27.3 %; test 0.6448 / 10.82 %, far 807 / 0.552 / 22.2 %). The trainlaw "Today" values do not use the adopted size axis (item 9). The held-out curve contradicts the "about 0" claim (item 4). |
| engineering:system-design | Traced each Gate B1 input back to the rental item that produces it. Found that the edge-view baselines and the two-seed floor have no producer in V1 (item 1). Found that the paired val baseline has no hash link to R5 (item 3). |
| engineering:architecture | Treated the CI decision rule as one decision with three conflicting texts (item 2). Treated the rescaled-kind target as a trade-off that must be written in ADR-0014 (item 10). |
| engineering:documentation | Each fix names the section, the line, the file and the exact replacement text. |
| huggingface-skills:huggingface-best | No model choice in scope. Used only to check that the 4B gate stays on val (it does) and that test is read once (item 14 finds a second read). |
| huggingface-skills:huggingface-datasets | Checked the rows-root path to the HF repo (pinned `rows_revision` plus manifest SHA, fixed in r1). No new item. |

Web checks: the TCIA LIDC-IDRI collection page, "Annotation and Markup Issues/Comments" items 1-8. Item 2: reader
order is not consistent across XML files. Item 6: LIDC-IDRI-0396 (139.xml) had a wrong SOP Instance UID at
position 1420, fixed 2018-06-28. Item 7: LIDC-IDRI-0510 has internalStructure = 5 (invalid). Item 8: 8 patients
have two timepoints. Also listed there: the 0101 XML was replaced on 2012-03-21, and 044.xml and 191.xml lack
ratings (2018). Licence: CC BY 3.0 [web: https://www.cancerimagingarchive.net/collection/lidc-idri/, fetched
2026-11-08].

New measurements for this review (read-only, `python3 -I`, script
`docs/implementation2.0/reviews/final_r2_data_check.py` run on `docs/implementation2.0/notes/b1_out/geom.json`):

- Today's off-centre rows against the train shifted law, with the doc's own apparent size axis: worst 0.557,
  TV 0.253, KS 0.109. The doc reports 0.759 / 0.108 / 0.109 (raw `ext_px`). The full mix gives 0.375 / 0.350 / 0.331
  (doc: 0.507 / 0.266 / 0.331). 77.5 % of rescaled windows change cell when the apparent size is used.
- Null distribution: random train-patient subsamples of the train shifted law against the full law (500 draws,
  seed 0). With 78 patients: TV p50 0.060, p95 0.100; KS p50 0.028, p95 0.052; worst-ratio p5 0.299. With
  38 patients: TV p50 0.092, p95 0.159; KS p95 0.086. The whole val gives TV 0.077, KS 0.061, worst 0.376.
- Size-bin shares of the target law: 0.149 / 0.471 / 0.305 / 0.075. Apparent size of the rescaled draw-0 windows:
  0.620 / 0.234 / 0.067 / 0.080. Items that can reach each apparent size bin with either neighbour side:
  4,549 / 1,718 / 924 / 1,074 (of 7,341).
- Padded windows (beyond the slice): test shifted 0.0000, val shifted 0.0046, train shifted 0.0038,
  val rescaled 0.0480, train rescaled 0.0341. Fallback and kept-side: 0 everywhere.
- Fetch sidecars `data/lidc_slices/*.json` (999 files): 4 series have kept slices != annotated SOP UIDs (item 8).
- Held-out IoU_tc on the 128-item set: e2 0.5113 [measured: runs/rl/eval/ckpt/heldout.json], step 250 0.634
  [docs/research/notes/full_rl_2026-11-08.md:37], final 0.6590 [runs/rl/eval/final/heldout.json]. `rl.json`
  `log` has no eval key in any of its 771 entries.

---

## Items

### 1. [blocking] [B, 00_front] Gate B1's hard edge-view rule and the "two-seed" noise floor have no producer in the paid plan

- Problem: Gate B1 has a hard rule: "val IoU on the `edge` view: CI lower bound of (B1 - R5) > -0.02 (R5 value from
  A3)" (line 3722). B.8.4 says A3 supplies "R5 evals on the val edge view ... R5 val eval with a second sampling
  seed". B.6.2 (line 3707), B.7 (line 3753), B1.9 (line 3307) and B.R2 #49 say "A3 evaluates R5 on val twice with
  two sampling seeds". But A3 does not do this. `a3_measure.yaml` runs only `lidc-val-shifted` and
  `lidc-test-shifted` (A3, Plan file). A3 step 7 says a second val sampling "is not planned". Front 6.2 says "A3 is
  not needed for Gate B1". The V1 table item 3 (line 470) evaluates the B1 model on one val view only (6.9 min), so
  the B1 model's edge eval is not there either. A3 (V1 item 4) is also skippable ("Skip if the ledger is above
  $1.90"). As planned, a hard Gate B1 rule cannot be computed after V1. Adding the two edge evals (2 x 6.9 min)
  makes V1 216 min = $2.31 at 27.5 s, which is above the $2.20 cap [estimate: arithmetic on front 6.2 inputs].
- Fix:
  - Front 6.2, V1 item 3: "Val evals of the B1 model: `lidc-val-shifted` and `lidc-val-edge` (2 x 6.9 min)". New
    item 3b (not skippable): "R5 (`rl-full-3ep`) on `lidc-val-edge`, 6.9 min". Recompute the totals and the cap
    (216 min = $2.31 at 27.5 s; 197 min = $2.10 at 23.65 s). Owner decision 2: either raise the V1 cap to $2.35, or
    keep $2.20 with `BudgetStop` at $1.75 for the B1 run (about 270 steps at 27.5 s, a multiple of 20 = 260).
  - Add the edge view to the A3 bundle and to `a3_measure.yaml` only if A3 is moved before item 3. Otherwise item 3b
    is its own non-skippable item.
  - Delete every "two sampling seeds" / "A3 two-seed" phrase (lines 3307, 3707, 3753, 4693; B.8.4 row A3). Replace
    it with the split-half floor of A3 step 7 (see item 5). In B.8.4, change the A3 row to "R5 eval on the val edge
    view (V1 item 3b)".
  - Front 6.2: replace "A3 is not needed for Gate B1" with "A3 is not needed for Gate B1 except as the paired
    baseline (item 3); the R5 edge baseline is V1 item 3b".

### 2. [major] [A, B] Three different interval rules decide the same gates

- Problem: A1 step 5 says the percentile patient interval is the decision value on the whole val (78) and on test
  (158). It also says a gate on fewer than 50 patients uses `signal.cluster_t` (t interval of the patient means,
  n - 1 df). B.6.2 says "The gate decision uses the BCa interval". B.8.4 asks A1 for "BCa for the decision". No task
  builds a BCa function (A1's files list `draws` and `cluster_t` only). Gate B3 runs on val-select (about
  39 patients after the re-deal): A says cluster t, B says BCa. The rules are "fixed before the run", but the doc
  fixes two rules. A borderline gate can pass under one rule and fail under the other.
- Fix: keep A1's rule, which is built and tested. B.6.2 "Pairing and CI": replace "The gate decision uses the BCa
  interval; the percentile interval is reported beside it" with "The gate decision uses the A1 gate CI: the
  percentile patient interval when the gate set has >= 50 patients (Gate B1, whole val), and `signal.cluster_t` on
  the per-patient mean differences when it has < 50 (Gate B3, val-select). The other interval is printed beside it
  with the cluster count." Change B.8.4 row A1 to the same words. Change B.R2 #2 "BCa for the decision" to "A1 gate
  CI". For paired subset rules (far misses, near IoU), apply `cluster_t` to per-patient mean differences over the
  patients that have items in the subset. Print that patient count (far subset: 35 patients on val-select, r1
  measurement).

### 3. [major] [A, B] The paired baseline of Gate B1 has no proof that it used the same images

- Problem: Gate B1 pairs the B1 model's val answers with `runs/r5/eval/eval_val/items.jsonl` by `item_id`. R5 scored
  `/workspace/data/lidc/eval_views/lidc-val-shifted`, a view built on the box. It has no hash [measured:
  runs/r5/eval/eval_val/final_eval.json `view`]. The doc freezes only `lidc-test-shifted` and `lidc-val-edge`
  (B1.6). The val view is rebuilt in B1.0/B1.6 and overwrites the only local copy,
  `data/eval_views/lidc-val-shifted`, built 2026-10-20 [measured: ls -la; manifest `target_set_sha256` 2da645bf...].
  So the B1-vs-R5 pairs rest on "same code" only. A-G3 cannot catch a different view. Its pass rule is a TOST on the
  mean difference, and two different shift draws give about the same mean IoU. Only the per-item correlation
  separates them, and the doc makes it a "second line". Also, B1.6 line 3263 says A3's "`view_manifest_sha256` ...
  must equal the frozen entry". The frozen entry holds the content hash, and A3 step 1 says
  `view_content_sha256`. The two texts disagree.
- Fix:
  - B1.0, before any rebuild: run `final_eval --hash-view data/eval_views/lidc-val-shifted --out
    runs/a3/val_pre_rebuild.json`. Rebuild val into a new folder. Require an equal content hash and
    `target_set_sha256`, then replace the old folder. Add `lidc-val-shifted` as a third entry in
    `docs/reference/frozen_views.json` and in DG15.
  - A3: evaluate `rl-full-3ep` on the frozen val view (already in the plan). In Gate B1 (B.6.2), make **A3's val
    items the paired baseline**: the same frozen view, the same code and the same box as the B1 eval. Keep R5's
    items as a second line. Then a view or host difference cannot enter the B1 - baseline difference. For R5
    0.618, write "A3 value; R5 0.618 for reference".
  - A-G3 (A.6) and A3 step 5: add a hard view-identity rule for val and test. The Pearson correlation of per-item
    IoU_tc between A3 and R5 must be >= `r_sh`. `r_sh` is the Spearman-Brown K = 8 reliability from A3's split
    halves: `r_sh = 2 r_h / (1 + r_h)`, with `r_h` the correlation of the answer 0-3 and 4-7 means. Use the 95 %
    patient-bootstrap lower bound of `r_sh`. A lower correlation means a different view: stop and pair only on
    rebuilt views.
  - B1.6 line 3263: replace "`view_manifest_sha256`" with "`view_content_sha256`".
  - CLEAN10 step 4: "keep every `data/eval_views/lidc-val-*` folder until B1.0 has recorded its content hash".

### 4. [major] [B] The extra-steps confound is about as large as the expected data effect, and its cited evidence does not exist

- Problem: ADR-0014 "Confound" says "Evidence that more centred steps add about 0: the held-out IoU plateau over the
  last epochs". Gate B1 reports "the slope of R5's held-out IoU over the last 300 of its 770 steps
  (`runs/rl/train_rl/rl.json` eval series)". `rl.json` has no eval series: no eval key in 771 log entries
  [measured]. The held-out set (128 items) has three points: 0.511 (e2), 0.634 (step 250), 0.659 (step 770)
  [measured: heldout.json files; full_rl_2026-11-08.md:37]. That is +0.025 over the last 520 steps. A linear
  extrapolation gives +0.014 for 300 more steps [estimate]. Gate B1's own decomposition expects +0.018 from the
  data fix. So rule 1 (CI lower bound of B1 - R5 > 0) can pass on extra training alone. The near-centre control is
  report-only.
- Fix:
  - ADR-0014 "Confound": replace the "about 0" sentence with "held-out IoU_tc rose 0.634 -> 0.659 from step 250 to
    770 (128 items) [measured]; a linear extrapolation gives up to +0.014 for 300 steps; the confound is the same
    size as the expected data effect".
  - Gate B1: add a hard difference-in-differences rule. With the same patient draws, the CI lower bound of
    (far gain - near gain) > 0. Far gain = mean (B1 - baseline) over offset > 0.3 items; near gain = the same over
    offset < 0.15. Extra centred training lifts both subsets; the data fix lifts mostly the far one. Add this rule to
    B1.9's power simulation.
  - Delete the "slope ... `rl.json` eval series" bullet. Replace it with "held-out 128-item points at steps 250 and
    770 (two points; no series exists)".

### 5. [major] [A, B] The noise floor is undefined, and it does not measure what the text says

- Problem: A1 says the CIs exclude training-seed variance, so "A3 reports a sampling-noise floor, so a later gain can
  be compared with it". A3 step 7 reports "the paired difference of the two K = 4 means with its patient CI, and its
  standard deviation times 1/sqrt 2". It does not say whether this is the per-item SD or the SD of the mean. B.6.2
  says "A point gain must exceed it", but this rule is in no gate bullet. A per-item SD (about 0.1) would block any
  realistic gain. The SD of the mean (about 0.003) is already inside the paired patient CI, because each B1
  per-item value is a fresh K = 8 sample. And neither measures training-seed variance, which is the risk A.7 row 7
  names.
- Fix:
  - A3 step 7: define two keys. `noise_floor.item_sd = sd_i(d_i) / sqrt 2`, with `d_i` = (mean of answers 0-3) -
    (mean of answers 4-7). `noise_floor.mean_sd` = the patient-bootstrap SD of `mean_i(d_i)`, divided by sqrt 2.
    State: "These are K = 8 sampling noise. The paired patient CI of every gate already contains them."
  - B.6.2 "Noise floor": replace "A point gain must exceed it" with "Reported beside each gain: `mean_sd`.
    Training-seed variance is not measured; a gain is a single-run result." Add the same sentence to DOC11 Known
    issues.
  - A.7 row 7 and B.7 row "Gains are seed noise": the mitigation is "stated as a limitation; a second-seed B1 run is
    the only measurement, owner decision, not in $3". Do not cite the sampling floor as the mitigation.

### 6. [major] [B] B1.9 cannot get its inputs before the spend it is meant to price

- Problem: B1.9 runs "$0, before any spend", but its item noise comes from "the A3 two-seed R5 evals" or, "without
  A3, the K = 8 answer spread / sqrt 8". A3 runs inside V1, the same rental as the B1 run, so it does not exist
  before V1. There is no two-seed eval (item 1). `runs/r5/eval/eval_val/items.jsonl` holds per-item K-means only, so
  the answer spread is not available either [code: experiments/phase2/final_eval.py:26-30]. The owner would price
  V1 with an undefined simulation.
- Fix: B1.9 noise model (Mac, R5 items only). For each item, the new B1 value is the simulated mean plus noise with
  SD `sqrt(mu_i (1 - mu_i) / 8)`. This is an exact upper bound, because Var(X) <= mu (1 - mu) for X in [0, 1], and
  `mu_i` is the per-item K-mean. The baseline values stay fixed (existing items). Run the simulation twice, with
  that bound and with noise 0, and report P(pass) for both. Rules are fixed from the bound run. After V1, re-run
  with A3's split-half `item_sd` for the record only. Delete "plus the A3 two-seed R5 evals when they exist".
  Update the B1.9 test: with noise 0 and no effect, P(pass) of the IoU rule = 0.

### 7. [major] [A, B] The test hash is frozen before the contour-to-slice check, and test hits have no rule

- Problem: A3 step 1 builds the test view and writes its content hash before the rental. B1.6 and DG15 require later
  builds to reproduce it. B1.7 (DG17) checks every series, but it resolves hits only "by dropping the items from
  train or val". If a test series has a z hit, two outcomes are possible. If the fix changes masks, it changes
  `target_set_sha256` and DG15 fails, with no rule for what follows. If the hit is ignored, test labels are known
  to be wrong. The time line (front section 1) puts A's Mac work and B1 in parallel, so B1.7 is not ordered before
  A3 step 1.
- Fix:
  - Front section 1, step 2, and A3 step 1: "B1.7 (`check_series_z.py`) runs before A3 step 1 builds and hashes any
    eval view."
  - B1.7 Gate DG17, add: "Test hits: the test view keeps its items and labels (R5 pairing). Every final test report
    adds a sensitivity line that excludes the affected items (ids listed in `frozen_views.json` under
    `known_label_issues`). A test fix that changes labels is applied only before A3 step 1. If one is needed after
    it, the owner decides in writing, and the frozen entry is re-issued with the reason."

### 8. [minor] [B] B1.7 ignores the free evidence already on the Mac

- Problem: B1.7 says no DICOM is on the Mac, so the duplicate-z check needs a TCIA download. But the fetch step wrote
  one sidecar per series from every DICOM header (`slices`, `annotated_sops`, and for 5 series `annotated_z`,
  `matched_by_z`) [code: experiments/phase1/fetch_lidc_pixels.py:115-118]. A duplicate z, or a skipped ambiguous
  slice, shows as `slices != annotated_sops`. 4 of 998 series show it: uid suffix `...295628` (6 vs 8), `...380786`
  (27 vs 28), `...228428` (5 vs 10), `...074972` (37 vs 48) [measured: data/lidc_slices/*.json]. `lidc.ct_slices`
  drops slice keys that are not within 0.01 mm of a pylidc z [code: src/objdel/data/lidc.py:303-312]. So contours on
  a dropped slice lose their image, and contours on a shared z may sit on the wrong image.
- Fix: B1.7 step 1 adds check (d): read every `data/lidc_slices/<series>.json` (`python -I`) and list the series with
  `slices != annotated_sops` or `matched_by_z > 0`. These go to step 2 (DICOM header check) with the (a)-(c) hits.
  DG17 counts them. I17 and I35 add a fixture sidecar with `slices < annotated_sops`.

### 9. [minor] [B] The "Today" coverage numbers use the size axis that the doc rejects

- Problem: ADR-0014 and B1.1 set the cell size to the apparent size `extent_px x S0 / side`. The "Today" values of
  DG7 (0.759), DG8 (0.108), DG9b (0.507 / 0.266 / 0.331) and B.1 item 3 come from `trainlaw.py`. That script's
  `cell()` uses raw `ext_px` for rescaled rows. With the adopted axis they are 0.557 / 0.253 / 0.109 (off-centre)
  and 0.375 / 0.350 / 0.331 (full mix). 77.5 % of rescaled windows change cell [measured: this review].
- Fix: B.6.1, DG7/DG8/DG9b "Today": 0.557 / 0.253 / 0.109 and 0.375 / 0.350 / 0.331, tagged
  `[measured: apparent size, final_r2_data_check.py]`. B.1 item 3: the same numbers. Add a U-test that `cell_coverage`
  and `trainlaw`-style baselines call `Window.apparent_px`.

### 10. [minor] [B] Matching the rescaled kind to the shifted law removes most of what rescaling adds

- Problem: of the rescaled draw-0 windows, 62.0 % have an apparent size under 10 px. The target gives that size bin
  14.9 % of the slots. Bins 2 and 3 (20 px and up) get 38.0 % of the slots, but only 924 and 1,074 items can reach
  them [measured]. So the scheduler fills rescaled slots mostly from the few large-nodule items, over and over, and
  drops the small-apparent windows. Those windows are the scale shift the rescaled kind was added for. ADR-0014
  does not state this trade-off.
- Fix: ADR-0014 "Size axis": add the trade-off with these numbers. Option (i), the default: keep the shifted target
  for both kinds and report in DG14 the number of distinct source items per (rescaled, cell) and the maximum reuse
  of one item. Option (ii): schedule the rescaled slots against the rescaled kind's own draw-0 apparent law, with
  DG7-DG9 per kind. Add this as B.8.3 item 3b (owner picks).

### 11. [minor] [B] DG16 has no reference for "how far is noise"

- Problem: DG16 shows the train law against the whole val as 0.755 / 0.077 / 0.061. That is below the SD3 limits
  (0.8 / 0.10 / 0.06) on two counts, and the reader cannot tell whether this is noise. Measured null for 78 random
  train patients: TV p95 0.100, KS p95 0.052, worst-ratio p5 0.299 [measured]. So the val TV is noise, but the val
  offset KS (0.061) is above the 95th percentile. Val (and test, median 0.262) offsets run slightly larger than
  train (0.247).
- Fix: `schedule_report` DG16 adds the null quantiles from `n_val_patients` random train-patient subsamples (500
  draws, seed 0) beside each distance. The DG16 row in B.6.1 shows "TV 0.077 (null p95 0.100); KS 0.061 (null p95
  0.052)". It stays a report.

### 12. [minor] [B] Two offset definitions are mixed in B.1

- Problem: B.1 item 1 cites "test offset median 0.265 and 37.5 % above 0.3". These come from the diagnosis note,
  whose table has 825 test items above 0.3 [rl_plateau_diagnosis_2026-11-08.md:24-38]. The B.1 table, under "one
  definition for all of Phase B", uses `geom.json`: median 0.262, 36.7 % (807 items) [measured]. B1.6 also says
  0.262. The test row shows patients as "-", but it is 158. No task says where Gate B1 gets each val item's offset
  (`geom.json` is a session note).
- Fix: B.1 item 1: "median 0.262, 36.7 % above 0.3 (807 of 2,200) [measured: geom.json, union-box definition]; the
  diagnosis note uses another centre and gives 0.265 / 37.5 %". Front section 1 and the BRIEF line get the same
  numbers. Test row: "158 / 2,200". B1.6: `build_eval_views` writes `view.offset_frac = windows.offset_of(...)` in
  each record, and the gate scripts read it from the frozen view, not from `geom.json`.

### 13. [minor] [B] The val halves reuse the split seed and deal disagreement patients by count

- Problem: B1.0 deals val patients with `splits.deal(groups, ..., seed=0)` inside each `lidc_stratum`. `lidc_split`
  chose the val patients with the same `_key(0, gid)` hash order [code: src/objdel/data/splits.py:59-70, 91-101].
  The halves are then a fixed function of the original deal order. Also, the `dis=1` stratum is dealt by
  `deal_weighted` on disagreement count in `lidc_split`, but by patient count here. The aim "balance
  disagreement" is not met when a few patients hold most disagreement slices (44 items on val).
- Fix: B1.0: `seed="val_half:0"` (a new string seed), and for the `dis=1` stratum use `splits.deal_weighted(...,
  weight=disagree)`, as `lidc_split` does. The test in `tests/test_splits.py` adds: the disagreement-item counts of
  the two halves differ by at most the largest single-patient count.

### 14. [minor] [B, 00_front] Test can be read for two candidates

- Problem: B-S2 ends with "the test eval of the final candidate" (B.4 step 4, B.8.1). The 4B gate (front section 3,
  owner decision 17) is decided after the 2B B runs, and "test is read once, on the final candidate of whichever model
  wins". If the 4B arm runs, test has already been read for the 2B model.
- Fix: B.4 step 4 and B.8.1 B-S2: drop the test eval. Add: "Test is read once, after the owner closes the 4B gate
  (front section 3), on the single final candidate. B-S2 ends with the val evals." Move the 10.6 min test cost to
  the session that follows the 4B decision.

### 15. [minor] [B] The edge view's fallback and padding rules are not defined

- Problem: B1.6 defines the edge offset as "the largest legal shift in a random direction", with crowding checks. It
  does not say what happens when the boundary point is crowded, or when no direction is admissible. (The shifted
  view falls back to the policy frame and counts it [code: build_eval_views.py:60-61].) Boundary windows are the
  ones most likely to extend beyond the slice. The pool drops padded windows (DG18 = 0), so training never shows
  them. A hard guard gate on `edge` then also tests padding, a second factor [padded share: val shifted 0.46 %,
  val rescaled 4.8 %, measured].
- Fix: B1.6: "Edge draw: up to 64 directions from `row_rng(1, "edge", item_id)`; take the first whose boundary
  window is not crowded. If none: fallback = the policy frame, counted in the manifest `fallback`. Report the padded
  share of the edge view; Gate B1's edge rule uses non-padded items, and the padded items are a report line." Test:
  edge fallback is counted, and every edge offset lies on the admissible boundary.

### 16. [minor] [A, B] "Code identity" names only the generator files

- Problem: `frozen_views.json` will say "equal to the R5 box view by code identity (generator files unchanged since
  commit 734b1ec)". The content hash includes `image_sha256` and `target_set_sha256`, which come from
  `dataset._prepare` (render, certified targets). `src/objdel/train/dataset.py` changed on 2026-11-02, 11-07 and 11-08
  (be52b1c) [measured: git log]. Whether R5's box view used the same render code depends on the commit R5 ran at,
  and the doc does not record it.
- Fix: A3 step 1 and B1.6: record `render_files` (`targets/augment.py`, `targets/frame.py`, `targets/*` used by
  `build_target`, `train/dataset.py`, `experiments/phase2/build_eval_views.py`) with their last commits, plus the
  R5 code commit from `runs/r5` (ledger or `CODE_SHA` once A3's sync writes it). Keep the statement "not measured".
  The view-identity rule of item 3 is the measured check.

### 17. [minor] [DOC] DOC11 known issues 4 and 5 are still marked unverified, but they are verified

- Problem: DOC11 "Known issues" items 4 (reader identity not consistent across scans) and 5 (0396, 0510) say
  "[unverified]". The TCIA page states both (issues 2, 6 and 7). B1.7 (c) says only "the 2018-06-28 SOP UID
  correction", without the file and position.
- Fix: DOC11 item 4: "[web: TCIA LIDC-IDRI, Annotation and Markup Issues item 2]". Item 5: "LIDC-IDRI-0396
  (139.xml): wrong SOP Instance UID at position 1420, fixed 2018-06-28; LIDC-IDRI-0510: internalStructure = 5,
  invalid [web: same page, items 6-7]". Add item 5b: "0101 XML replaced 2012-03-21; 044.xml and 191.xml lack ratings
  (items 4-5); ratings are not used here". B1.7 (c): "139.xml carries the corrected SOP Instance UID at position
  1420".

### 18. [minor] [B] Val-select stops being clean if B1 is rebuilt after Gate B1

- Problem: Gate B1 uses the whole val, which includes val-select. B.7 allows a rebuild after a B1 result ("set
  rescaled 0 in a rebuild"). A rebuild decision based on whole-val numbers lets val-select shape B1 before it judges
  Gate B3.
- Fix: B.7 row "Rescaled share hurts" and B.6.2: "Any change to the B1 rows after Gate B1 is decided on val-shape
  numbers and training logs only. Val-select numbers from the B1 run are reported, not used to change B1." Add the
  same sentence to ADR-0014.
