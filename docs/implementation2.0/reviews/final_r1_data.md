# Final review, round 1: data quality and statistics (Implementation.md)

Reviewer role: medical imaging data scientist and statistician (LIDC, multi-rater, leakage, bootstrap).
Scope: whole doc, A and B1 first. No earlier `final_r*_*.md` file exists, so nothing here repeats a fixed item.
Counts: **3 blocking, 10 major, 14 minor (27 items).**

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA catalog skill covers data audit or bootstrap. The doc's own catalog result (`dicom-series-preflight`) is relevant to item 12 only. Nothing recommended beyond that. |
| engineering:code-review | Each claim checked against code: `final_eval.per_item` (K-mean, item 3), `augment.draw/shift/row_rng`, `build_eval_views._record`, `dataset.prompts`, `cache.voters`, `lidc.py` P1-P6, `splits.py`. |
| engineering:debug | Reproduced the key numbers before judging: miss rate 0.10818 is the K-mean rule; R5 val-select CIs by nodule vs patient; val-shape vs val-select vs train cell laws (script run on `docs/implementation2.0/notes/b1_out/geom.json`). |
| engineering:system-design | Traced the data path artifact by artifact (store -> pool -> schedule -> rows root -> HF -> box). Found the checks that run on the store but not on the artifact trained on (item 7) and the prefix problem (item 6). |
| engineering:architecture | Target-law choice treated as a decision with options (val-shape, train law, pooled); item 1 picks the train law with measured reasons. |
| engineering:documentation | Each fix names the section, file, function and gate text to change. |
| huggingface-skills:huggingface-best | Not applicable to data review (no model choice); used only to flag that the 4B gate uses test (item 4). |
| huggingface-skills:huggingface-datasets | Checked the HF upload path for the rows root (item 25) and the third-party `YongchengYAO/LIDC-IDRI-Lite` dataset card that is the source of the duplicate-z list (item 12). |

Web checks: TCIA LIDC-IDRI collection page and wiki (licence CC BY 3.0; known issues: readers not consistent across
scans, LIDC-IDRI-0396 SOP UID fix 2018-06-28, LIDC-IDRI-0510 invalid attribute, 8 two-timepoint patients)
[web: https://www.cancerimagingarchive.net/collection/lidc-idri/ ; https://wiki.cancerimagingarchive.net/pages/viewpage.action?pageId=1966254];
`YongchengYAO/LIDC-IDRI-Lite` card (5 duplicate-z series: 0085, 0146, 0418, 0572, 0979)
[web: https://huggingface.co/datasets/YongchengYAO/LIDC-IDRI-Lite].

New measurements made for this review (script inline, read-only, `python3 -I`):
- Cell-law distances, cells (SIZE_EDGES 10/20/40 px x OFFSET_EDGES 0.15/0.2/0.3), from `geom.json`
  [measured: this review]:

  | Set | items | patients | offset < 0.15 | offset > 0.3 | extent < 10 px | extent >= 40 px |
  |---|---|---|---|---|---|---|
  | val-shape | 449 | 40 | 0.156 | 0.421 | 0.241 | 0.031 |
  | val-select | 635 | 38 | 0.214 | 0.342 | 0.134 | 0.063 |
  | train shifted (draw 0) | 7,341 | 557 | 0.199 | 0.325 | 0.149 | 0.075 |

  Target val-shape vs val-select: TV 0.169, offset KS 0.099. Target train-shifted vs val-select: TV 0.063, KS 0.037.
- R5 val-select bootstrap (B = 4,000, seed 0) [measured: this review, `runs/r5/eval/eval_val/items.jsonl` x geom.json]:
  IoU 0.5993; CI by item [0.576, 0.622] (width 0.046), by nodule (106) [0.564, 0.634] (0.069), by patient (38)
  [0.556, 0.641] (0.085). Far misses (offset > 0.3, 217 items, 86 nodules, 35 patients): 0.318, nodule CI width
  0.128, patient CI width 0.141. Near IoU (offset < 0.15, 136 items, 28 patients): patient CI width 0.098.

---

## Items

### 1. [blocking] [B] The coverage target is a 40-patient sample that differs from the half that judges the gates

- Problem: B1.0, B1.2, B1.4 step 1 and D7-D9 take the target cell shares from val-shape (449 items, 40 patients,
  94 nodules). The offset law is the generator's (`augment.shift`: uniform over the admissible rectangle), so the
  only thing a patient sample adds is the size mix, and 40 patients estimate it badly. Measured: val-shape has 24.1 %
  of objects under 10 px against 13.4 % in val-select and 14.9 % in train; 42.1 % at offset > 0.3 against 34.2 % and
  32.5 %. Distance to the held-out half: val-shape TV 0.169 / KS 0.099; the train-split shifted law (7,341 items,
  557 patients, same generator) TV 0.063 / KS 0.037 [measured: this review]. So the plan tunes the training mix to the
  noise of 40 patients, over-weights small far nodules (Bae et al., B.2 row 3: harder-than-needed data loses
  accuracy), and uses doctor-derived geometry of val patients to shape training. It also risks `ValueError` on target
  cells that the train geometry cannot fill.
- Fix:
  - B1.4 step 1: target = `cell_shares` of the draw-0 shifted windows of the train items (the D16 law; geometry from
    the train store, no val data). Rename the CLI flag `--target` to `--target-law train-shifted` and record
    `target_law`, `target_sha256` in `manifest.schedule`.
  - B1.3 `schedule_report`: D7-D9 against the train law (implementation checks, unchanged limits); new report row
    D16 = coverage of the schedule against val-shape and val-select (worst ratio, TV, KS), report only.
  - B.6.1 "Today" for D7-D9: recompute against the train law; replace the D16 "Today" value.
  - B.6.2: since val no longer shapes B1, Gate B1 uses the whole val (78 patients, 200 nodules, 1,084 items; R5 0.618,
    far misses 27.3 %) for more power. Keep the val halves only for B3 (B3.7 tunes on val-shape, Gate B3 on
    val-select). Update the B.3.1 diagram ("val-SHAPE half" box -> "train shifted law") and B.7 row 1.

### 2. [blocking] [A, B] Gates resample nodules, but the independent unit in LIDC is the patient (scan)

- Problem: A1 adds a patient CI but every gate uses the nodule CI: A-G3, A-G4 (A.4 A4, A.6), Gate B1 and Gate B3
  ("the Phase A1 nodule-level bootstrap", B.6.2), and E G10 / GP5 (90 % nodule-level TOST). In LIDC all nodules of a
  scan share the same four reading sessions, and reader 1 in one XML is not the same radiologist as reader 1 in
  another [web: TCIA LIDC-IDRI known issue 2]; scanner, slice thickness and reconstruction are also per scan; 8
  patients have two scans [code: src/objdel/data/lidc.py:17]. Nodules of one patient are therefore correlated in
  both the image and the reference. Measured on R5 val-select: patient CI width 0.085 vs nodule 0.069 (+23 %);
  far-miss 0.141 vs 0.128 [measured: this review]. The two val halves differ by 0.046 in mean IoU (0.645 vs 0.599
  [measured: valsplit.json]), which shows large patient-level variance. Nodule CIs make the gates anti-conservative,
  and a false pass unlocks spend.
- Fix: in A1, make `ci95_patient` (key `group_id`) the gate interval and keep `ci95_nodule` as a secondary line.
  Change the gate text of A-G3, A-G4, Gate B1, Gate B3 and E G10/GP5 from "nodule" to "patient (group_id), paired".
  With 35-38 clusters, state the cluster count beside each CI and use the BCa interval (or a cluster t interval with
  n_clusters - 1 df) for the gate decision; keep the percentile interval as the reported value. Add an A.9 test:
  two patients with two nodules each, perfectly correlated within patient, give a patient CI wider than the nodule CI.

### 3. [blocking] [00_front] "miss" has three meanings; section 0 contradicts every measured number

- Problem: section 0 defines a miss as "an image where the best answer scores IoU below 0.1". Every measured number
  uses the K-sample mean: `final_eval.per_item` stores the mean IoU_tc of the 8 answers
  [code: experiments/phase2/final_eval.py:26-30], and the 10.8 % (0.10818181818181818) is the share of items with
  that mean below 0.1 [measured: this review]; A5 step 4 also uses the K-mean. B3.2 uses a third meaning (per answer,
  `iou_tc < miss_iou`). An implementer who follows section 0 (best of 8) gets a lower miss rate and silently moves
  the Gate B1 baseline (31.8 %) and Gate B3. This breaks the doc's own "one word for each thing" rule.
- Fix: section 0: "miss: an image whose mean IoU_tc over its K = 8 answers is below 0.1"; "found image: not a miss".
  Add a new word "missed answer: one answer with IoU_tc below 0.1" and use it in B3.2, B3.4, B3.7 and B3.9. Add
  "the K-mean rule" to A5 step 4 and to `docs/glossary.md` (A9).

### 4. [major] [00_front] The 4B gate reads the test split

- Problem: section 3 "Gate for the 4B arm": "test IoU on shifted views is still below about 0.70 and the miss rate is
  still above about 6 %". A model choice made on test makes the final test number a selected value. B.6.2 says
  "test once, final candidate only".
- Fix: gate on val (whole val after item 1, or val-select): "val IoU_tc below about 0.68 and val miss rate above
  about 7 %" (owner sets the numbers; val R5 is 0.618 / 12.6 % against test 0.645 / 10.8 %, so the test thresholds do
  not transfer one to one). Owner decision 12 changes the same way.

### 5. [major] [A] A-G4 lets the test split veto the answer-choice method

- Problem: A-G4 requires "The test paired difference of the chosen method ... has the same sign as on val". If it
  fails, the default does not change, so test chooses. A.4 A4 says "It never chooses the method"; the gate text
  contradicts it.
- Fix: A-G4: "The test paired difference is reported with its patient CI as a direction line. It is not a pass
  condition." Delete "has the same sign as on val" from A-G4 and from the A4 Gate A paragraph.

### 6. [major] [B] The trained prefix of the schedule is never checked for coverage

- Problem: B1.2 step 2 assigns cell slots to positions "by a shuffle" over the whole 300-step list. The plan then
  trains on prefixes: the B2 probe uses the first 640 rows (20 steps, B2.3); section 6.3 says the $3 budget allows
  only about 160-190 steps with `BudgetStop`. A random prefix meets the target only on average; D7-D9 and D11 are
  computed on the full 9,600 rows, not on what is trained.
- Fix: B1.2 step 2: assign cell slots per block of 20 steps (640 rows) by largest remainder, so every block meets the
  target counts exactly; shuffle positions only inside a block. Build with `--steps` equal to the budgeted steps
  (owner decision 1), and record `planned_steps` in `manifest.schedule`. B1.3: run D7-D9 and D11 on every block
  prefix (20, 40, ...); B1.5 refuses `--steps` that is not a block multiple. Add a property test: every 20-step
  prefix has exact per-kind cell counts.

### 7. [major] [B] Leakage gates run on the store, not on the rows root that is trained

- Problem: D1-D3 call `quality.leakage(store_root)`. A builder bug (wrong `--split`, a stale view folder) would put
  val or test items in `data/views/lidc-train-b1` and no gate would see it. The rows root, not the store, goes to HF
  and to the box.
- Fix: add `quality.rows_leakage(rows_root, store_root) -> dict` and hard gate D1b: every `source_item_id` is in the
  store train split; its `group_id` is not a val or test patient; no row `image_sha256` equals any val or test view
  image hash, or the frozen test view hashes (`docs/reference/frozen_views.json`); 0 of each. Call it in
  `build_rl_rows` before the manifest is written and again in `preflight`. Test: a planted val item in a 30-row
  fixture fails D1b.

### 8. [major] [B] No label-noise audit of the training rows

- Problem: B1 audits geometry and leakage but not the reference quality of the rows it trains on. B1.md Q3 (keep or
  down-weight the 4.1 % low-contrast train items, contrast < 5 grey levels) is left open; D14 only reports the share.
  Rows where the doctors disagree strongly give a near-random reward. A2 already builds the leave-one-out doctor
  score, which is the per-item agreement measure, and it is free on the Mac.
- Fix: B1.3 add `quality.agreement(rows_root) -> dict`: per scheduled window, the mean doctor LOO IoU_tc from
  `objdel.eval.ceiling.loo_caches` + `score_mask` (doctor side only; CPU build). Report the distribution per cell and
  the share with doctor LOO IoU_tc < 0.3, joint with low contrast. Add D14b (report) and owner decision B.8.3 item 4:
  "keep all / drop items with doctor LOO IoU_tc < X / weight 0.5", decided before B1.4 runs. Default: keep all
  (honest exam), but the number is shown.

### 9. [major] [B] The plan trains on the exam's exact generator and deletes every out-of-law check

- Problem: B1 trains on many draws of `augment.draw(side=S0)` with `row_rng(0, kind, ...)`, the same generator that
  builds the test view [code: experiments/phase2/build_eval_views.py:46-56]. A3 deletes `lidc-test-rescaled` and
  `lidc-test-rotated`; B1.0 and B1.6 build `--kinds shifted` only. A test gain then measures fit to one synthetic
  shift law, not robustness to framing. Clinical crops will not follow `augment.shift`.
- Fix: keep the local val views `lidc-val-rescaled` and `lidc-val-rotated` (they exist [measured: ls data/eval_views])
  and add a report line to Gate B1 and Gate B3: paired patient-CI difference vs R5 on those two views, with a
  non-inferiority note at -0.02 (report, not a gate). B1.0 rebuilds val with all kinds (`--kinds` default all). DOC11
  "Limitations": "the exam is a synthetic shift law; training matches it".

### 10. [major] [B] No power check for Gate B1 and Gate B3 before spend

- Problem: the gates hinge on small held-out sets: far misses in val-select are 217 items from 35 patients (patient
  CI width 0.141 around 0.318), near IoU 136 items from 28 patients (width 0.098) [measured: this review]. The point
  rule asks for -0.079 on far misses and the near rule for a lower bound above -0.02. Nobody has estimated the chance
  that a real effect passes. B.7 row 2 then lets "the owner accept a near-miss in writing", which is not a
  pre-registered rule.
- Fix: new $0 task B1.9 `experiments/phase2/gate_power.py`: from R5 per-item values, simulate a B1 effect (convert
  25 % of far misses, chosen by patient, into hits with IoU drawn from the found-image distribution of the same cell;
  leave other items with item-level noise from two K = 8 resamples of A3 answers), run the paired patient bootstrap
  of each Gate B1 and Gate B3 rule 1,000 times, print P(pass) per rule and jointly. Put P(pass) in B.8.3 item 1 so
  the owner prices the run with it. Replace "the owner may accept a near-miss in writing" with fixed rules decided
  before the run.

### 11. [major] [00_front, B] The model decision rests on the ceiling that A2 calls unfair

- Problem: section 3 reason 1, its "Free check" and B.1 item 7 use "found images 0.722 against 0.740". The 0.740 is
  the old ceiling (doctor vs k - 1 with plain IoU, model vs k with IoU_tc), which A.1 item 2 says favours the model
  [docs/research/notes/rl_plateau_diagnosis_2026-11-08.md: ceiling definition]. The found subset is also selected by
  the model's own success. "A bigger model cannot gain there" says nothing about misses, which are the loss.
- Fix: section 3 reason 1 and B.1 item 7: "[to be confirmed by A2]". After A2 (doctor side free before the rental,
  model side from A3 answers), recompute on val: fair `diff` on found images and on all images, by patient CI. Free
  check text: "this bounds the gain on found images only; misses are not covered". Section 3 decision stays "keep
  2B", but its reason cites the A2 number.

### 12. [major] [B] The contour-to-slice check covers 4 series from a third-party list

- Problem: B1.7 checks only 0085, 0146, 0572, 0510. The duplicate-z list comes from one third-party dataset card
  (it also names 0418 and 0979) [web: YongchengYAO/LIDC-IDRI-Lite]. `lidc.slice_index` maps every ROI by nearest z
  [code: src/objdel/data/lidc.py:162-164], so any series with ambiguous z can mis-assign a contour. The XML check is
  free and runs in seconds on all series.
- Fix: `check_series_z.py` runs on every series of the store: (a) z values shared by different `imageSOP_UID`
  (`lidc.sop_z`); (b) ROIs whose `|imageZposition - nearest pylidc zval|` exceeds half the slice spacing; (c) the
  local XML for LIDC-IDRI-0396 carries the 2018-06-28 SOP UID correction [web: TCIA known issue 6]. New hard gate
  D17: 0 unresolved hits (each resolved by SOP mapping or by dropping the items, written in the report).

### 13. [major] [00_front] The view bank has no task and targets data that B1 deletes

- Problem: sections 2.3, 2.4 and 4 make `data/viewbank.py` the data format ("~2.2 GB for 22,023 views"), and owner
  decision 11 archives "the 51,400 PNG view files after the view bank passes its check". No phase task builds it
  (grep: `viewbank` appears only in lines 140, 170, 296). The 22,023 views are the fixed 3-view set that B1.8
  deletes; B1 trains on a new 9,600-window rows root of PNGs. The only check named is "RGB equals gray on 300 + 1
  sampled" (open question 6), but every record carries `image_sha256` of the RGB pixels [code:
  src/objdel/train/dataset.py:87-88,108-110].
- Fix: either remove the view bank from sections 2.3, 2.4, 4 and owner decision 11 until a phase owns it, or add a
  task (Phase C or E) that builds it from the B1 rows root and verifies, for every image, that the expanded RGB
  equals the PNG pixels by `_pixels_sha` (all images, not a sample). No PNG is archived before that check passes.

### 14. [minor] [A] A-G3 uses "CI contains 0", which is not an equivalence test

- Problem: a noisier run passes more easily. E G10 already uses TOST correctly.
- Fix: A-G3: "the 90 % patient CI of the mean per-item IoU_tc difference lies inside ±0.01 (TOST); a wider CI is
  'not shown', never a pass". Drop "contains 0".

### 15. [minor] [B] D13 checks 150 of 9,600 windows

- Problem: a 1 % defect rate is found with probability 1 - 0.99^150 = 0.78. Materializing all windows takes about a
  minute (B.8.2).
- Fix: D13 re-renders every scheduled window (box, side, `image_sha256`); limit 0 mismatches.

### 16. [minor] [B] D3 omits a measured hit and the review step

- Problem: `neardup.json` has test->val 1 hit at cos >= 0.98 besides test->train 1 [measured:
  docs/implementation2.0/notes/b1_out/neardup.json]; B1.md G4 asked that each hit be reviewed.
- Fix: D3 Today: "1 test->train, 1 test->val, 0 val->train". Rule: every hit is listed with both item ids and a
  side-by-side PNG in the audit output, and marked reviewed before the gate passes.

### 17. [minor] [B] Disagreement weight 4 moves the training mix away from the exam

- Problem: 289 of 7,341 train items are disagreement items (3.9 %) [measured: data/prepared manifest]; weight 4 inside
  each cell gives about 14 % of rows [estimate: 4p / (1 + 3p)], against 44 of 1,084 (4.1 %) in val. This contradicts
  the "match the exam" rule of B1 and is not argued.
- Fix: B1.2: state the expected share, report it in D14, and add to B.8.3 item 3 "disagreement weight 4 (inherited
  from `DISAGREEMENT_OVERSAMPLE`) or 1 (exam law)".

### 18. [minor] [B] Coverage is checked on off-centre rows only

- Problem: D7-D9 skip the 15 % centred rows, which sit at offset about 0.01, so the offset law of the full training
  mix is never shown.
- Fix: add D9b (report): offset KS and cell TV of the full mix (all kinds) against the target law.

### 19. [minor] [B] The size axis ignores the scale of rescaled windows

- Problem: `cell` uses `extent_px` (slice pixels). For a rescaled window S changes (85 % are S = 100 [measured:
  B1.md F4]), so the object appears half as large as its cell says, and rescaled rows are matched to shifted cells.
- Fix: define the size axis as apparent size `extent_px x S0 / S` for every kind (equal to `extent_px` for centred and
  shifted), or schedule the rescaled kind against its own draw-0 law. Write the choice in ADR-0014.

### 20. [minor] [B] The text says "test is read once", but test already shaped the design

- Problem: B.1 item 3 and SD3 measure the training mix against test (TV 0.299, worst ratio 0.77); the B.1 table and
  the diagnosis use test results per offset. The final test number is therefore not an untouched estimate.
- Fix: B.1 item 3: replace the test comparisons with val-select and train-law numbers (item 1 table). Add to B.6.2
  and to DOC11 Limitations: "the test split informed the diagnosis (per-offset results of R5); from 2026-11-08 every
  design number uses val only".

### 21. [minor] [A, B] The R5 test view is not on the Mac and was never hashed

- Problem: R5 scored `/workspace/data/lidc/eval_views/lidc-test-shifted`, built on the box; `final_eval.json` has no
  manifest SHA and no local copy exists [measured: runs/r5/eval/eval_test/final_eval.json; ls data/eval_views].
  A-G3 pairs a rebuilt view with R5 by `item_id`, and B1.6 freezes a rebuilt view. Equality is inferred from "same
  code" (the generator files last changed in commit 734b1ec, 2026-11-02 [measured: git log]), not measured.
- Fix: A3 step 1 computes the B1.6 content hash of the rebuilt test view and writes it to the pre-rental record;
  B1.6 must reproduce that hash after the B1.0 rebuild. `frozen_views.json` records `generator_commit` and states
  "equal to the R5 box view by code identity, not measured".

### 22. [minor] [A] A2 stratifies by raters but leaves out voters

- Problem: strata `k=2`, `k>=3` come from rater counts (281 / 803 [notes/A1.md F3]); the left-out set is
  `cache.voters`, which can be smaller than the raters [code: src/objdel/reward/cache.py:108-111].
- Fix: A2: "k = len(cache.raters)"; `ceiling.json` reports `n_voters_lt_raters` per stratum, and the doctor side
  averages over voters only (state it).

### 23. [minor] [B] The B2 rule uses a t bound over 10 dependent steps

- Problem: per-step differences of two diverging trajectories are autocorrelated; the one-sided t lower bound assumes
  independence.
- Fix: B2.3 `b2_summary.py`: paired bootstrap over the 320 prompts of steps 11-20 (same prompts in each arm),
  resampling items, as the decision number; keep the t bound as a second line.

### 24. [minor] [B] Gate B1 credits data for the effect of 300 more RL steps

- Problem: the B1 run continues `rl-full-3ep` for 300 steps; R5 is the start model. Part of any gain is extra
  training (the diagnosis measured +0.005 to +0.008 training IoU_tc per 100 steps late in the run
  [docs/research/notes/rl_plateau_diagnosis_2026-11-08.md]).
- Fix: Gate B1 report adds the near-centre (offset < 0.15) gain as the control: a data effect shows as far gain much
  larger than near gain. ADR-0014 states the confound.

### 25. [minor] [B, 00_front] Rows-root storage and pinning are inconsistent

- Problem: section 4 and owner decision 10 plan a separate private HF dataset repo (owner OK needed); B1.4 step 5
  uploads to the model repo `svenk029/objdel-sft` under `rows/`. The run config pins no HF commit, so the box could
  load a different upload.
- Fix: pick one place (owner decision 10), record `rows_revision` (HF commit) and the manifest SHA-256 in
  `b1_rows.yaml`; `train_rl` asserts the SHA-256 of the downloaded `manifest.json`. The rows root is LIDC only (CC BY
  3.0, attribution in the repo README [web: TCIA]); RIGA stays on the Mac (section 6.4 rule 7).

### 26. [minor] [DOC] The data card content is not specified, and one line contradicts A1

- Problem: DOC11 adds datasheet headings with "no new numbers" but names no known issue; `docs/reproducibility.md`
  cites "item bootstrap CIs in `final_eval`", which A1 replaces.
- Fix: DOC11 "Known issues" lists: population = slices with >= 2 distinct reader masks (P5 copies count once); masks
  about 1 px generous (P2); readers not consistent across scans; 8 two-scan patients; LIDC-IDRI-0396 and -0510 TCIA
  notes; duplicate-z result of item 12; the exam is a synthetic shift law; test informed the diagnosis (item 20).
  Reproducibility line: "patient-cluster bootstrap CIs (A1)".

### 27. [minor] [B] The val halves are unbalanced

- Problem: the parity-of-first-byte rule gives 449 / 635 items and IoU 0.645 / 0.599 [measured: valsplit.json]. After
  item 1 the halves serve only B3, but B3.7 still tunes on a half that differs from the judging half.
- Fix: B1.0 `val_half` deals patients with `splits.deal` stratified by `lidc_stratum` (multi-reader slice count),
  seeded, so the halves balance items and disagreement. Re-run `valsplit.py` and update the section 1 table.
