# Final review round 13: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).

Scope. I read the sources by section, not the built file:
- A: header, A.1, A.3.1, A.4 header, A1 (steps 1-6, 5b, output keys), A2, A3 (files, plans, header block, steps
  0-8), A8, A.6, A.7, A.8.
- B: header, B.1, B.3.4 (ADR-0015), B.4, B1.0-B1.9 (whole), B3.9, B.6.1, B.6.2, B.7, B.8.3, B.8.4, the test rows
  U09, U24, U55, P04, I17, I24, I30, E04, the gate and task maps, the round 10-12 logs.
- Front: 6.1 (whole), 6.2 (launch config and plan list), time line step 2. 99_back "Data and quality".
- E, D, CLEAN: grep for `build_rl_rows`, `train_exclusions`, `frozen_views`, `data/eval_views` (no data item found).
- `review-log.md`: rounds 10-12 of front, A and B, so handled items are not repeated. The round-12 data-1 fix (the
  exclusion list) is in B.md and I re-checked it against the code; it holds.

Scripts (all `.venv/bin/python -I`, repo path as argument, 2026-11-09):
`(planning probe, not kept: r13data)` (dry run of B1.7 checks (a) and (b) on every store series), plus the inline checks
quoted below.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in the data scope. B1.7's DICOM helper stays "recommend only". No catalog lookup, nothing installed. |
| engineering:code-review | Read the code behind the data claims: `fetch_lidc_pixels.{fetch,incomplete}` (item 2), `lidc.{sop_z,slice_index,ct_slices,scans}`, `build_data.main` and `store.write` (the "mapped" rebuild path), `splits.{deal,lidc_stratum}`, `build_eval_views._record`, `augment.{shift,draw,neighbour_side}`, `dataset.{prompts,grpo_rows}`, `cache.rl_eligible`. |
| engineering:debug | Reproduced each number before leaning on it: store counts per split, the val-half search (s = 78), CR1 on R5 val, the five sidecar hits, the far / near patient counts. A crash of my own sidecar loop (KeyError `slices`) found item 3. |
| engineering:system-design | Traced XML -> sidecar -> store -> `pool(exclude)` -> rows root -> DG1b, and asked which planned check could still move the train hit count. Running checks (a) and (b) now answered it (item 1). |
| engineering:architecture | Checked which section owns each order rule: B.4 step 1 and B1.0 ("do first") against A3's header and B1.7's val-fix branch (item 5). |
| engineering:documentation | Each item quotes the source line and the code line or measurement behind it. |
| huggingface-skills:huggingface-best | No model choice in scope. Front section 3 unchanged. |
| huggingface-skills:huggingface-datasets | Read the YongchengYAO/LIDC-IDRI-Lite card (the third-party duplicate-z list B1.7 cites). No Dataset Viewer endpoint is on the data path; the rows root goes up by `HfApi.upload_folder` (B1.4 step 5). |

Web checks (2026-11-09):
- HF card YongchengYAO/LIDC-IDRI-Lite [web: https://huggingface.co/datasets/YongchengYAO/LIDC-IDRI-Lite]: lists
  LIDC-IDRI-0085, -0146, -0418, -0572, -0979, "two or more DICOM slices sharing one `ImagePositionPatient` z",
  dropped from that repo. Matches B1.7's citation. Item 1 measures what this means for our store.
- TCIA LIDC-IDRI "Annotation and Markup Issues" (cited by B1.7 (c)): re-read in round 12; unchanged.

## Claims verified (no item)

- Store: train 7,341 items / 1,247 nodules / 562 series / 557 patients; val 1,084 / 200 / 78 / 78; test 2,200 / 401 /
  160 / 158 [measured: `data/store/lidc/*.parquet`]. Matches A.1 item 1 and Gate B1's "78 patients, 200 nodules,
  1,084 items". No store item has a null image (0 of 10,625), so a contour on a slice the fetch dropped is already
  absent from the store.
- A1 step 5 / B.6.2: on R5 val, item mean 0.61819, patient mean 0.64722, CR1 [0.58371, 0.65267], top 5 patients
  24.9 % of items, 2-78 items per patient [measured: `runs/r5/eval/eval_val/items.jsonl` with `group_id`].
- B1.0: strata 41 / 23 / 2 / 12; first passing s = 78; 40 / 38 patients; 533 / 551 items; 24 / 20 disagreement items
  [measured: `splits.deal` on `data/eval_views/lidc-val-shifted/lidc/val.jsonl`].
- Gate B1 subset rules on the whole val use the percentile interval: the far subset (offset > 0.3) has 406 items in
  75 patients, the near subset (< 0.15) 206 items in 52 patients, 49 patients in both; far misses 111 / 406 =
  27.3 % [measured: R5 items with `(planning probe, not kept: r11data)`]. Both subsets are >= 50 patients, so A1 step 5's
  "< 50 patients -> `cluster_t`" switch does not fire for Gate B1 (it does for Gate B3 on val-select).
- B1.7 check (d): 5 sidecar hits (074972, 179144, 228428, 295628, 380786); every one of the 800 store series has a
  sidecar [measured: `data/lidc_slices/*.json`].
- The "mapped" rebuild path is byte-deterministic (`store.write`, ascending item ids, fixed write options), so the
  val / test parquet SHA-256 check in B1.7 is a real test, not a guaranteed failure.
- `pool` step 1 `rl_eligible(item.raters)`: `rl_eligible` is `len(stack) >= 2` [code: src/objdel/reward/cache.py:103-105],
  so a list of rater names works.
- B1.5 reads a schedule root with `oversample=1`, so the B1.2 disagreement quota is not multiplied by
  `DISAGREEMENT_OVERSAMPLE` a second time [code: src/objdel/train/dataset.py:217-227].

## Items

No blocking and no major item. Every item below is minor: either a planned $0 check (DG17, I17, the DG17 re-run)
fails first and points at it, or it has no effect on a result.

### 1. MINOR. B1.7 checks (a) and (b) now measured: they add no series. The train hit set is final at 4 series / 27 items, and the val "duplicate-z" patients are not hits

Files: B (B1.7 "Exclusion list" bullet "steps (a) and (b) have not run yet, so the count can grow"; the paragraph
after the hit table; B.6.1 DG17 "Today"; B.7 row "A val label issue breaks the frozen val view"), A (A3 step 0),
99_back (8c).

Doc lines:
- B.md B1.7: "Today that is 27 items in 4 series ... steps (a) and (b) have not run yet, so the count can grow."
- B.md B.6.1 DG17 Today: "the third-party duplicate-z list names 2 val patients (17 items)".
- B.md B.7: "B1.7 hit on a val series (the duplicate-z list names LIDC-IDRI-0085 and -0146)".
- 99_back 8c: "The third-party duplicate-z list also names two val patients, LIDC-IDRI-0085 (2 items) and -0146 (15
  items) ... [keep, for test and val]".

Measurement [`(planning probe, not kept: r13data)`, 2026-11-09; `lidc.sop_z` on `data/raw/lidc_xml`, pylidc z values, all
800 store series]:
- (a) z shared by different `imageSOP_UID`: 3 series, LIDC-IDRI-0013, -0017 (train, 20 items) and -0659 (test, 4).
- (b) ROI z farther than half the slice spacing from the nearest scan z: 1 series, LIDC-IDRI-0566 (train, 3 items).
- Both lists are inside the (d) list. So the train hit set is 0013 (15), 0017 (5), 0365 (4), 0566 (3) = 27 items;
  the test hit is 0659; there is no val hit.
- LIDC-IDRI-0085, -0146, -0418, -0572, -0979: pylidc z values have 0 duplicates (267, 309, 226, 447, 249 slices) and
  check (a) finds no shared z. The fetch keeps the DICOM whose SOP UID the contour cites
  [code: experiments/phase1/fetch_lidc_pixels.py:8-12], so the image of each item is the cited one. The duplicate
  DICOM of the card is never read. This closes the r4 data-8 carry-over with a measurement.

Fix: B1.7: replace "steps (a) and (b) have not run yet, so the count can grow" with "checks (a) and (b), dry-run on
2026-11-09, hit only series already in the (d) list (a: 0013, 0017, 0659; b: 0566) [measured:
`check_series_z` dry run]; the count grows only if the XML or the store changes". B.6.1 DG17, B.7 val row and 99_back
8c: "the third-party duplicate-z patients 0085 and 0146 are not hits: the XML cites one SOP UID per z and pylidc's z
list has no duplicate; no val item is in `known_label_issues` today".

### 2. MINOR. The "mapped" resolution cannot clear any of the four known train hits; only "dropped" can

Files: B (B1.7 "mapped" bullet).

Doc line: B.md B1.7: "**mapped**: fix the cause (contours mapped by `imageSOP_UID`; for a (d) hit, re-fetch the series
with `experiments.phase1.fetch_lidc_pixels --refetch-incomplete`, which rewrites its sidecar ...)".

Code and measurement:
- `--refetch-incomplete` re-fetches only "series with fewer slices than annotated z positions"
  [code: experiments/phase1/fetch_lidc_pixels.py:128-136, `incomplete`]. Sidecars: 380786 slices 27 / annotated_z
  27; 228428 5 / 5; 179144 9 / 9; 295628 6 / 8 [measured: `data/lidc_slices/*.json`]. So the flag re-fetches only
  295628, and it reads the same TCIA DICOMs, so the sidecar and check (d) come back the same.
- The fetch docstring records why: "LIDC-IDRI-0365's XML cites no SOP UID of its TCIA series, and readers 2 and 4 of
  LIDC-IDRI-0017 cite five absent UIDs" [code: fetch_lidc_pixels.py:12-14]. These are source-data facts.
- A global change of `lidc.slice_index` to SOP mapping would also touch val and test series, and B1.7 then stops on
  the val / test parquet SHA check.

Why minor: DG17 reads the re-run, which still lists the hit, so the gate fails at $0 and the engineer then drops.
Mac time only. Fix: add to B1.7: "For the four (d) train hits the only resolution is 'dropped' (the re-fetch reads
the same TCIA files). 'Mapped' applies only to an (a) or (b) hit whose fix is limited to that series (a per-series
SOP map), so the val and test parquet SHA-256 stay equal."

### 3. MINOR. Check (d) "reads every sidecar", but the folder also holds `selection.json`; the count is 998 sidecars, not 999

Files: B (B1.7 step 1 (d)).

Doc line: "read every `data/lidc_slices/<series>.json` ... Today 5 of 999 sidecars hit (d)".

Measurement: `data/lidc_slices/` holds 999 `*.json` files: 998 series sidecars and `selection.json` with keys
`control`, `disagree` (Phase 0 selection lists). A loop over `*.json` that reads `d["slices"]` raises
`KeyError: 'slices'` on it [measured: my first loop, 2026-11-09]. Caught at $0 on the first run. Fix: "read the
sidecar of every store series (`data/lidc_slices/<series uid>.json`, one per store series; 800 today) ... Today 5 of
998 series sidecars hit (d)."

### 4. MINOR. Check (c) names "139.xml", but five files have that name, and LIDC-IDRI-0396 is not in the store

Files: B (B1.7 step 1 (c)), I17.

Doc line: "(c) the local XML for LIDC-IDRI-0396 (139.xml) carries the corrected SOP Instance UID at position 1420".

Measurement: `data/raw/lidc_xml/tcia-lidc-xml/{185,186,187,188,189}/139.xml` all exist; the store holds 793 patients
and LIDC-IDRI-0396 is not one of them [measured: store `group_id`]. So (c) cannot change any item. Fix: "(c), for the
record only (0396 is not in the store): the XML whose `SeriesInstanceUid` is 0396's series ...". I17 keeps its
fixture.

### 5. MINOR. B1.0 hashes and freezes the val view before B1.7 runs, but two lines say no view is hashed before B1.7

Files: A (A3 header block, A3 step 0), B (B.4 step 1, B1.0, B1.7 "Val hits").

Doc lines:
- A3 header: "`check_series_z` # B1.7 first: no eval view is hashed before it".
- B1.7 and A3 step 0: "A val label fix is applied only before B1.0 hashes the view."
- B1.0: "(do first; one rebuild)", step 1 runs `final_eval --hash-view data/eval_views/lidc-val-shifted`; B.4 step 1
  orders "B1.0, B1.1 to B1.7" and puts B1.7 only before A3 step 1 and B1.4.

Why minor: item 1 shows no val hit, so the val-fix branch is empty today. Fix: B.4 step 1: "B1.7 step 1 runs first,
before B1.0 step 1 hashes `lidc-val-shifted`" (the XML check takes seconds), or delete the val-fix sentence and keep
only the `known_label_issues` rule.

### 6. MINOR. DG17's "a mapped series that the latest report still lists" can be read as failing every mapped series

Files: B (B1.7 "mapped" bullet and "Gate DG17").

Doc lines: "the re-run names it under `mapped` for the record" and "a 'mapped' series that the latest report still
lists are unresolved". I17 pins the meant rule ("passes when the re-run ... lists no hit for it"), so a literal build
fails I17 at $0. Fix: "a 'mapped' series that the latest report still lists as a hit".

### 7. MINOR. DG14b label-noise decision: its input exists only after B1.4, but the owner decides before B1.4; "drop" and "weight 0.5" have no mechanism

Files: B (B1.3 `agreement`, B.8.3 item 4), 99_back (8).

Doc lines: B1.3: "`agreement(rows_root) -> dict` (DG14b, report): per scheduled window ... The keep / drop /
down-weight rule is owner decision B.8.3 item 4, made before B1.4 runs." B.8.3 item 4: "keep all (default ...), or
drop items with doctor LOO IoU_tc < X, or weight them 0.5." The exclusion list takes only reasons `"DG17 ..."` and
`"DG3 ..."` (B1.7), and `schedule` has no weight input (B1.2).

Why minor: the default (keep all, 99_back 8) needs no mechanism. Fix: `agreement(prepared_root)` on the centred train
root (one value per item; the doctor masks do not change with the window shift), run before B1.4; "drop" writes
`"DG14b <value>"` entries into `train_exclusions.json`; remove "weight 0.5" or name where the weight enters.

### 8. MINOR. Two stale numbers in B

- B.7: "Gain comes from 300 more steps, not data". B1 is 260 steps (B1.5 step 3, ADR-0015 "260 more steps"). Fix:
  260.
- B.6.2 common rules: "Gate B3, val-select, about 39 patients". B1.0's dealt halves give 38 [measured above; B.7
  says 38]. Fix: "38 patients (B1.0, s = 78)".

### 9. MINOR (still open from earlier rounds; not re-argued). Listed so they are not lost

Re-read 2026-11-09; each text is unchanged in the sources:
- r12 data-2: A3 still says "`lidc-test-shifted` reaches the box when it exists on the Mac at bundle time"; the
  header block has no move of the test folder before `bundle` (front 6.1 step 4c requires it).
- r12 data-3: A3 step 7 still labels `item_sd = sd_i(d_i) / sqrt 2` "per-item K = 8 noise" (it is the SD of a
  two-sampling difference).
- r12 data-4: the B.md header still has no "this phase section is the rule" line; `B.spec.md` is still revision 2.
- r12 data-5: Gate B1 "point <= 20.5 % (25 % below R5's 27.3 %)" still does not say it is absolute.
- r11 data-3: A3 step 5 "The power of this test is not measured".
- r10 data-2, r10 data-3, r10 data-4, r9 data-4, r9 data-5 list, r5 data-5, r7 data-5: as listed in round 12.
- r4 data-8 is closed by item 1 (measured: no hit).

## Counts

Blocking 0. Major 0. Minor 9 (items 1-9).
