# Final review, round 3: data quality and statistics (Implementation.md)

Reviewer role: medical imaging data scientist and statistician (LIDC, multi-rater, leakage, bootstrap).
Scope: the whole doc, A and B1 first. I checked the 18 items of `final_r2_data.md` against the round-3 text (front
log "final round 2", A "Review log, part 3", B.R3 items 1-43). All 18 are fixed in the section that owns them, so
none is repeated here. The new items below come from cross-section checks and from new measurements.
Line numbers refer to `docs/implementation2.0/Implementation.md` as of 2026-11-08 (14,506 lines).
Counts: **1 blocking, 5 major, 6 minor (12 items).**

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No catalog skill covers bootstrap, split design or LIDC label QA. The doc's existing `dicom-series-preflight` recommendation (B1.7 step 2) still fits item 6. Nothing installed or newly recommended. |
| engineering:code-review | Checked claims against code: `splits.deal/deal_weighted/lidc_stratum` (item 3), `lidc.parse` / `slice_index` / `sop_z` / `ct_slices` (items 6, 11), `build_eval_views._record` and `augment.draw` (item 2 de-risk), the b1_pixels seed script (item 10). |
| engineering:debug | Reproduced before judging: re-ran the B1.0 val-half rule as written (item 3); recomputed all 1,084 val shifted windows with today's code (item 2); joined the B1.7 sidecar hits and the third-party duplicate-z list to the store splits (items 6, 8). |
| engineering:system-design | Traced each V1 item to the plan file that runs it; found three different V1 lists and a plan file that would exceed the cap (item 1). |
| engineering:architecture | Treated "which baseline does Gate B1 pair with" as one decision with two texts (item 2), and the gate interval as one estimand (item 4). |
| engineering:documentation | Every fix names the section, the text to change and the exact replacement. |
| huggingface-skills:huggingface-best | No model choice in scope. Used only to confirm that model selection stays on val and test is read once (it does; item 7 adds that val already chose R5). |
| huggingface-skills:huggingface-datasets | Checked the HF dataset card that B1.7 cites for duplicate-z series (`YongchengYAO/LIDC-IDRI-Lite`); two of its five patients are in our val split (item 6). |

Web checks (2026-11-08): the TCIA LIDC-IDRI collection page (items 2, 4-8 of "Annotation and Markup Issues") for
items 6, 8 and 11; the HF card `YongchengYAO/LIDC-IDRI-Lite` for the duplicate-z list (0085, 0146, 0418, 0572,
0979). Cluster-robust (CR1) variance of a ratio mean: standard sandwich form (Cameron and Miller 2015, "A
practitioner's guide to cluster-robust inference", J. Human Resources 50(2):317-372, doi 10.3368/jhr.50.2.317;
the paper states that few clusters need a small-sample correction and t critical values with G - 1 df) for item 4
[web: econpapers.repec.org/RePEc:uwp:jhriss:v:50:y:2015:i:2:p:317-372, fetched 2026-11-08].

New measurements for this review (read-only, `.venv/bin/python -I`, nothing written under `data/` or the repo):

- Val: 78 patients, items per patient min 2, median 8, max 78; the 5 largest patients hold 24.9 % of the 1,084
  items. Item-weighted mean IoU_tc of R5 0.6182; patient-weighted mean 0.6472 [measured: runs/r5/eval/eval_val/
  items.jsonl joined with data/store/lidc/val.parquet `group_id`].
- The held-out set of the full RL run is 128 items of `lidc-val-shifted` from 51 of the 78 val patients
  [measured: runs/rl/eval/final/heldout.json `per_item_iou` keys joined with val.parquet].
- B1.0 val-half rule exactly as written (`deal` with seed `"val_half:0"` in the `dis=0` strata, `deal_weighted` by
  `disagree` in `dis=1`, item counts as `Patient.slices`): 38 / 40 patients, **454 / 630 items**, 22 / 22
  disagreement items; the `dis=1` stratum splits 4 / 8 patients. The same `deal` with seed `"val_half:78"` (the
  first of 0-199 that meets |delta items| <= 4 % and |delta disagreement items| <= 4): 533 / 551 items, 24 / 20.
- All 1,084 val shifted windows recomputed with today's `augment.draw(bare, row_rng(0, "shifted", id), side=S0)`
  equal the offsets and sides in the local `data/eval_views/lidc-val-shifted` (built 2026-10-20): 1,084 / 1,084.
- B1.7 sidecar hits by split: train 27 items (LIDC-IDRI-0013: 15, -0017: 5, -0365: 4, -0566: 3); val 0; **test
  4 items (LIDC-IDRI-0659, series suffix 074972)**. Third-party duplicate-z patients: val LIDC-IDRI-0085 (2 items)
  and **LIDC-IDRI-0146 (15 items)**; train LIDC-IDRI-0572 (12); 0418 and 0979 not in the store
  [measured: data/store/lidc/{train,val,test}.parquet].
- Trainer checkpoint `runs/rl/train_rl/trainer/checkpoint-250` exists on the Mac (LoRA + optimizer) [measured: ls].
- Split leakage re-checked: patient and series overlap 0 for train-val, train-test, val-test; 557 / 78 / 158
  patients; 7,341 train items, all `rl`, 289 with `disagreement` [measured: store parquet, prepared train.jsonl].

---

## Items

### 1. [blocking] [00_front, A, B] Three different V1 item lists; A's plan file as written duplicates an eval, reads test in V1 and breaks the $2.20 cap

- Problem: the owner approves one V1 list (99_back decision 2), but the doc holds three.
  - Front 6.2 (lines 511-534) and 99_back decision 2: 5 items. Item 4 = R5 on `lidc-val-edge` (not skippable);
    item 5 = A3-val, R5 on `lidc-val-shifted` only, skippable; "No nsys window, no A3 test half in V1"; total
    195.1 min = $2.08 at 27.5 s; cap $2.20 or $2.30.
  - B.8.1 (lines 4371-4381): 6 items. Item 6 = A3 val part **plus "the A6 nsys window (about 10.5)"**; total
    205.6 min = $2.20 (at the cap); cap option $2.35. A says V1 has no nsys step (A.4, A3 plan text).
  - A3 (lines 1121-1176), A.4, A-G3 and A.8: `a3_measure.yaml` runs R5 on val shifted, **val edge** and **test**
    (24.4 min); "the two val steps are required"; A.8: "The front's V1 table must list these three steps".
  If the main agent runs `a3_measure.yaml` as written after front item 4, R5 on val edge runs twice, the test split
  is read in V1 (front 6.2 moved A3's test half to V3), and V1 at 27.5 s is 30 + 137.5 + 13.8 + 6.9 + 24.4 =
  212.6 min = $2.27 > the $2.20 cap [estimate: arithmetic on front 6.2 inputs]. A paid rental with an item list the
  owner did not approve breaks the owner rule "owner approves the amount and the item list" (6.1 step 6).
- Fix (front 6.2 is the one list; 99_back decision 2 already copies it):
  - A3: `a3_measure.yaml` keeps one `final_eval` step: `lidc-val-shifted` (seed 11, saved answers). Delete its
    val-edge step (front item 4 is a step of `b1_rows.yaml`: `final_eval --ckpt rl-full-3ep --view
    {data}/eval_views/lidc-val-edge --split val --out {rd}/eval_val_edge_r5`, on the A3 code, so its answers are
    saved too). Move the test step to a new `experiments/phase2/configs/a3_test.yaml` for V3 item 3. Replace "The
    two val steps are required" with "The val-shifted step is front 6.2 item 5 and follows its skip rule". A.8:
    "The front's V1 table lists one A3 step (val shifted); val edge is front item 4; test is V3 item 3". A-G3 val
    part reads the val-edge answers from `eval_val_edge_r5`.
  - B.8.1: delete the nsys part of row 6 (the A6 capture is V2 item 3, as A and the front say). Totals become
    195.1 min = $2.08 (27.5 s) and 175.9 min = $1.88 (23.65 s); cap options $2.20 or $2.30 (not $2.35). Change
    B.8.3 item 11 to the same two numbers.
  - Add a front 6.1 check 4d: "the plan files run in V1 (`b1_rows.yaml` steps, `a3_measure.yaml`) contain exactly
    the `final_eval` views of the 6.2 table and no `--split test`" (grep test in `tests/test_repo_hygiene.py`).

### 2. [major] [00_front, A, B] Gate B1's primary baseline is named two ways, and at the base step time it never exists

- Problem:
  - Front 6.2 notes: "The main rule pairs item 3 against R5's `runs/r5/eval/eval_val/items.jsonl`; item 5, when it
    runs, is a same-view R5 baseline." B.6.2 and A.4 say the opposite: A3's R5 val items are **the** baseline and
    R5's file is a second line.
  - Front 6.2 also says item 5 "is skipped at 27.5 s" (90 % rule: session spend at its start is (30 + 137.5 +
    20.7) min = $2.01 > $1.98). So at the base step time the primary baseline of B.6.2 is never made in V1, and the
    fallback ("R5's items only when the B1.0 pre-rebuild hash equals the frozen entry") becomes the plan. V2 item 3b
    runs after Gate B1 is due, so Gate B1 either waits for V2 or uses the fallback; no text says which.
  - A-G3 is "hard" (TOST and view identity), but no text says what Gate B1 does when A-G3 fails.
  - B1.0 step 2: "A difference stops B1.0" has no next step.
  - De-risk fact for B1.0: today's code reproduces the geometry of all 1,084 local val shifted windows (offset and
    side equal) [measured: this review]. A hash difference can then come only from rendering or targets.
- Fix:
  - Front 6.2 item 5 skip rule: replace "skipped when session spend is above $1.98" with "runs when session spend
    + $0.08 (its 6.9 min + imports) <= the cap". At 27.5 s: $2.01 + $0.08 = $2.09 <= $2.20, so it runs. Mark it
    "skippable only when the cap cannot cover it". Same words in 99_back decision 2 and B.8.1 row 6. 6.4 rule 2
    gets: "an item whose full cost fits under the cap is not stopped by the 90 % rule".
  - Front 6.2 notes: replace the "main rule pairs ... R5's items.jsonl" sentence with "the main rule pairs item 3
    with item 5 (A3's R5 val items, same frozen view, same box). R5's `eval_val/items.jsonl` is a second line, and
    the fallback baseline only under B.6.2's hash condition."
  - B.6.2 "Paired baseline": add "If A-G3's view-identity rule fails, Gate B1 still pairs with A3's items (same
    frozen view); only the R5 second line is dropped. If A-G3's TOST fails, Gate B1 is decided on A3's items and
    the report states the A3-R5 difference."
  - B1.0 step 2: "A difference: keep both folders, freeze the rebuilt view, record both hashes, and mark A3-val
    (front item 5) as the only Gate B1 baseline (not skippable). Add to B1.0 'Today': geometry of all 1,084 windows
    reproduces with current code [measured: final_r3_data review]."

### 3. [major] [B] The new val-half rule does not balance items; its own test fails

- Problem: B1.0 says the stratified deal "balances items and disagreement slices" and replaces the parity rule
  because 449 / 635 items is unbalanced. Re-running the rule exactly as written gives **454 / 630 items** (38 / 40
  patients), and the `dis=1` stratum splits **4 / 8 patients** [measured: this review]. `deal` balances patient
  counts, not items, and items per patient range from 2 to 78. `deal_weighted(..., weight=disagree)` balances
  disagreement items (22 / 22) but cannot keep "each stratum's two halves differ in size by at most 1 patient",
  which the B1.0 test asserts. So the test fails, and val-select would again hold 58 % of the items and the harder
  patients (the old halves gave IoU 0.645 / 0.599).
- Fix: B1.0 "Val halves", replace the deal rule with: "Deal with `splits.deal` inside each `lidc_stratum` stratum
  (all strata, `dis=1` included). Try seeds `"val_half:<s>"` for s = 0, 1, 2, ... and keep the first that meets both
  (a) |items(shape) - items(select)| <= 4 % of val items and (b) |disagreement items(shape) - (select)| <= 4. The
  criterion uses counts only, never a model score. Record `s` and both counts in `docs/reference/val_halves.json`."
  Today's val gives s = 78: 533 / 551 items, 24 / 20 disagreement items [measured: this review]. Test: replace
  "each stratum's two halves differ by at most 1 patient" with "the stored halves meet (a) and (b), the stored
  seed is the first that does, and no score field is read".

### 4. [major] [A, B] `cluster_t` on per-patient means estimates a different number than the gate's point value

- Problem: A1 step 5 and B.6.2 use `signal.cluster_t` below 50 patients: "the per-patient means, then the t interval
  of their mean". That interval is for the **patient-weighted** mean. Every point rule (Gate B3 "miss share falls
  by >= 25 %", the item means in `final_eval`, the percentile patient bootstrap at >= 50 patients) uses the
  **item-weighted** mean. On val they differ by 0.029: 0.618 item-weighted, 0.647 patient-weighted; items per
  patient run 2-78 and 5 patients hold 24.9 % of items [measured: this review]. On val-select (about 39 patients) a
  rule can pass its point on one estimand and its CI on the other, and the CI is not centred on the reported gain.
- Fix: A1 step 5, rename to `signal.cluster_t(x, groups, level)` = CR1 cluster-robust t interval of the
  item-weighted mean. With per-patient sums `S_p`, counts `N_p`, P patients: `m = sum S_p / sum N_p`;
  `se^2 = P / (P - 1) * sum_p (S_p - m N_p)^2 / (sum N_p)^2`; interval `m +- t_{P-1, (1+level)/2} se`. Paired rules
  apply it to per-item differences. Subset rules apply it to the items of the subset, P = patients with items in
  it. Tests (A1 step 6): with equal `N_p` it equals the t interval of the patient means (the current hand test);
  with unequal `N_p` its centre equals the item mean exactly. B.6.2 "Pairing and CI" and B.8.4 row A1: replace
  "on the per-patient mean differences" with "CR1 t interval of the item-weighted mean difference (A1 step 5)".

### 5. [major] [B] The difference-in-differences rule cannot separate extra training from the data fix, because near items sit at the ceiling

- Problem: Gate B1's hard rule "CI lower bound of (far gain - near gain) > 0" assumes extra steps lift near and far
  items alike. Near val items (offset < 0.15) score 0.751, above the legacy val doctor ceiling 0.734 (B.1 item 2
  table; A.1 item 2). So a generic gain cannot show on near items, and the DiD rule reduces to "far gain > 0". The
  doc's own evidence shows the generic gain does reach shifted items: the held-out set that rose 0.634 -> 0.659
  between steps 250 and 770 of centred-only training is 128 **shifted** val items [measured: heldout.json keys;
  full_rl note "128 items of eval_views/lidc-val-shifted"]. So ADR-0015's confound is not controlled.
- Fix:
  - ADR-0015 "Confound" and B.6.2 Gate B1: add "Near items are at the doctor ceiling (0.751 vs 0.734 legacy), so the
    DiD rule cannot detect a generic gain. A2 reports the near-subset headroom = mean(fair doctor score - R5) on
    offset < 0.15 items, with its patient CI. If the headroom upper bound is below 0.014 (the extrapolated generic
    gain), the DiD rule is reported as uninformative and the B1 result is stated as 'B1 vs R5', not 'data effect'."
  - Owner option (B.8.3 new item 12, $0.08): add to V1 one eval of `runs/rl/train_rl/trainer/checkpoint-250`
    (exists on the Mac) on `lidc-val-shifted`. Per subset, (R5 - ckpt250) is the measured gain of 520 centred
    steps. Gate B1 control rule then becomes: CI lower bound of [far gain - (300/520) x far gain(250 -> 770)] > 0,
    same patient draws. B1.9 simulates it with the 250 -> 770 per-item values once measured.

### 6. [major] [B, A] A val label hit has no rule that fits the frozen val view; two flagged val patients exist

- Problem: DG17 resolves a train or val hit "by mapping by `imageSOP_UID` or by dropping the items from train or val",
  and "a val change moves the val baselines: re-run `valsplit.py`". But B1.0 freezes `lidc-val-shifted` (DG15:
  content hash must equal the pre-rebuild hash), A3 and R5 pair on its 1,084 items, and val halves are stored.
  Dropping or remapping val items breaks DG15, the R5 pairing and the halves at once; the two hard gates
  (DG15, DG17) then cannot both pass. This is not hypothetical: the card that B1.7 cites lists LIDC-IDRI-0085 and
  -0146 as duplicate-z series, and both are in our val split (2 + 15 = 17 items, 1.6 % of val) [measured: this
  review; web: huggingface.co/datasets/YongchengYAO/LIDC-IDRI-Lite].
- Fix: B1.7 DG17, change the val rule to the test rule: "Val hits keep their items and labels in the frozen view
  (R5 and A3 pairing). Their ids go to `frozen_views.json` under `known_label_issues`. Every Gate B1 and Gate B3
  rule is computed on the paired items without them (the exclusion is the same in both arms), and the report
  prints the line with them. A val label fix is applied only before B1.0 hashes the view." Remove "re-run
  `valsplit.py`" from DG17. Add I-test: a planted val hit leaves the content hash unchanged and is excluded from
  both arms of a Gate B1 rule.

### 7. [minor] [B, A] Val already chose the R5 checkpoint; the doc does not say so

- Problem: the full-run note picks the final model over checkpoint-250 on 128 val shifted items ("the final model is
  the one to keep"), and the G-probe used the same set. These items come from 51 of the 78 val patients
  [measured: this review]. So R5 is mildly selected on val. B.7 says "No val data reaches training" (true), but the
  Limitations list (B.6.2, DOC11) does not say that val chose R5. The bias favours R5, so Gate B1 is conservative;
  Gate B3 (B3 vs B1) is not affected.
- Fix: B.6.2 "Limitations" and DOC11 item 8: add "128 val shifted items (51 patients) chose R5 over checkpoint-250
  and served the G-probe [docs/research/notes/full_rl_2026-11-08.md]; R5's val numbers are slightly optimistic on
  them. Gate B1 prints a sensitivity line without these 128 items." List the ids in `docs/reference/
  val_halves.json` under `r5_selection_items`.

### 8. [minor] [B, A] B1.7 already has a test hit; the decision is due before A3 step 1

- Problem: B1.7 "Today" lists 5 sidecar hits by series suffix only. One is a test series: LIDC-IDRI-0659 (suffix
  074972, 37 of 48 slices kept, 26 matched by z), 4 test items. The other four are train (27 items in LIDC-IDRI-0013,
  -0017, -0365, -0566) [measured: this review]. A3 step 0 lets a test-label fix happen only before A3 step 1, so the
  owner must choose now, but no text says the hit exists.
- Fix: B1.7 "Today": add the split, patient and item count of each hit (the numbers above). A3 step 0: "B1.7 already
  flags 4 test items (LIDC-IDRI-0659). Before step 1, the owner chooses: fix (DICOM header check of B1.7 step 2,
  then rebuild) or keep (`known_label_issues`, sensitivity line)." Add it to 99_back "Data and quality".

### 9. [minor] [A] The gate CI has no golden test

- Problem: A1 creates `tests/goldens/r5_test_items.jsonl.gz` with `item_id` and `iou_tc` only, and pins `ci95_item`
  and `ci95_nodule`. The patient CI is the gate CI for every gate, and A.1 already measured it
  ([0.6252288200704239, 0.6644106494958667]), but the golden has no `group_id`, so no test can pin it. Item ids hold
  the series, not the patient.
- Fix: A1 Files: add `group_id` to the golden rows (public LIDC patient ids). A1 step 6 and A-G1: "`ci95_patient` of
  IoU_tc equals [0.6252288200704239, 0.6644106494958667] exactly on the recorded numpy and machine, else within
  1e-12". Re-measure the gzip size.

### 10. [minor] [B] The near-duplicate rule is not implementable as written

- Problem: DG3 and `quality.near_duplicates` say "cosine >= 0.98 on centred 32x32 HU crops". The seed that produced
  the "Today" numbers subtracts each crop's mean before the cosine (a Pearson correlation) and pads edge crops with
  -1,000 HU [code: docs/implementation2.0/seeds/b1_scripts/b1_pixels.py:39-41, 67]. A raw cosine on HU crops is near 1
  for almost any two lung crops, because every value sits near -800 HU [estimate], so an engineer following the
  text would flag most pairs.
- Fix: B1.3 `near_duplicates` and DG3: "crop = 32x32 int16 HU at the nodule centre, padded with -1,000 to 32x32;
  subtract the crop mean, divide by its norm (+1e-6); flag pairs with dot product >= 0.98 (the b1_pixels.py rule)".
  Test: two crops that differ by a constant HU offset have score 1; a crop and its negative have -1.

### 11. [minor] [DOC] DOC11 Known issues: one wrong claim and one wording drift

- Problem: item 5 says the 0396 SOP UID fix does not affect this data because "the parser places a contour by
  `imageZposition`, not by SOP UID". True for the parser, but the CT fetch picks slices by SOP UID (`lidc.sop_z`
  feeds `fetch_lidc_pixels` [code: src/objdel/data/lidc.py:272-279]), so a wrong UID can drop a slice. That is why
  B1.7 (c) checks 139.xml. Item 8 says "every design number uses val only"; B.6.2 says "train or val only" (the
  target law is the train shifted law).
- Fix: DOC11 item 5: "The parser places contours by z; the CT fetch picks slices by SOP UID (`lidc.sop_z`). B1.7 (c)
  checks that the local 139.xml carries the corrected UID; the result is written here." Item 8: "train or val only".

### 12. [minor] [00_front] Stop rule 7 does not list the eval views that leave the Mac

- Problem: 6.4 rule 7 says nothing from `data/` leaves the Mac "except the B1 rows root and the K1 mask bundle". The
  V1 bundle ships every `data/eval_views/lidc-*` folder [code: experiments/phase2/vast/bundle.py:53-56], so the
  frozen val views (and the test view, once A3 step 1 builds it) leave too. All are LIDC, so the licence is fine, but
  the rule as written is false, and the test view reaches V1 although test is not read there.
- Fix: rule 7: "... except LIDC-only artifacts: the B1 rows root, the K1 mask bundle and the LIDC eval views named in
  the rental's table". 6.1 step 4c: "the V1 bundle holds `lidc-val-shifted` and `lidc-val-edge` and no
  `lidc-test-*` folder" (bundle `--view` filter or a pre-bundle move).
