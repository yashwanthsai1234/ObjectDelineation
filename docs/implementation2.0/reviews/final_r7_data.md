# Final review round 7: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).
Scope: the whole doc, read by section from the sources. Read in full: A.1, A1-A5, A8, A.6-A.8; B.1, B.3.4 (ADR-0015),
B.4, B1.0-B1.9, B.6.1, B.6.2, B.8.4, the B1 rows of B.9.1, B.9.3, B.9.7 and the B.9.8 sketches for
`test_data_windows.py` and `test_reward_spec.py`; front sections 0, 1.1, 3, 6.1, 6.2; 99_back "Data and quality";
DOC11 known issues; E4 E-G10 and GP5. Round-6 fixes in `review-log.md` (data-1, data-2, data-3) are not repeated.
Script for this round: `(planning probe, not kept: r7data)`.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in the data scope (B1.7's DICOM helper stays "recommend only"). No catalog lookup, no install. Nothing changed. |
| engineering:code-review | Read the code behind each data claim: `splits.{deal,lidc_stratum,lidc_split}`, `experiments/phase1/build_data.py:122-128` (how the split built `Patient`), `dataset.{_prepare,prompts,grpo_rows}`, `build_eval_views._record`, `augment.{row_rng,shift,draw,neighbour_side}`, `frame.{frame,crowded,pixels}`, `bundle.members`. Found items 1, 2 and 7. |
| engineering:debug | Reproduced numbers before using them: the val halves search (first s = 78; 40 / 38 patients, 533 / 551 items, 24 / 20 disagreement items), the seed-0 halves of the sketch (item 1), the CR1 interval on R5 val ([0.5837, 0.6527] around 0.6182; patient-weighted mean 0.6472; top 5 patients 24.9 %), the val k strata (281 / 330 / 473), train targets (23,469 in all three roots), padded windows per split and kind (item 3). |
| engineering:system-design | Traced the rows path pool -> schedule -> materialize -> manifest -> preflight -> box, and which gate reads which split (item 2). Checked the bundle filter: `lidc-train-*` views ship only when named [code: bundle.py:56-57], so the rows root never rides the bundle. |
| engineering:architecture | Checked which text owns each rule: B1.0 owns the val halves, U02 and the B.9.8 sketch both claim to test it (item 1); ADR-0015 assigns the near-subset headroom to A2, whose spec has no such output (item 5). |
| engineering:documentation | Checked each test row and sketch against the task text it claims to test (items 1-4, 7) and listed the round-6 minors that are still open (item 8). |
| huggingface-skills:huggingface-best | No model choice in scope; front section 3 left as is. |
| huggingface-skills:huggingface-datasets | Used for the web check of LIDC dataset cards on the Hub (the MedOtter LIDC-IDRI card names the eight two-scan patients and the five duplicate-z series), which cross-checks B1.7 and DOC11 known issue 3. No finding. |

Web checks: TCIA LIDC-IDRI page ("Annotation and Markup Issues" items 1-8) [web:
https://www.cancerimagingarchive.net/collection/lidc-idri/]; Armato et al. 2011 (1,018 scans from 1,010 patients;
eight patients have two scans under one patient ID) [web: https://pmc.ncbi.nlm.nih.gov/articles/PMC3041807/];
MedOtter/LIDC-IDRI Hub card (two-scan patients 0132/0151/0315/0332/0355/0365/0442/0484; duplicate-z series 0085,
0146, 0418, 0572, 0979) [web: https://huggingface.co/datasets/MedOtter/LIDC-IDRI]. No documented cross-patient
duplicate scan exists, so grouping by `group_id` (the TCIA patient id) is the right leakage unit. The store's
series surplus (train 562 - 557 = 5, test 160 - 158 = 2) fits seven of the eight two-scan patients.

## Claims verified (no item)

- Val records: 1,084 items, 78 patients, k = 2 / 3 / 4 for 281 / 330 / 473 items (so A2's 281 / 803), 44
  disagreement items; train 7,341 items, 557 patients, 289 disagreement items; every store item is `rl`
  [measured: `data/eval_views/lidc-val-shifted/lidc/val.jsonl`, `data/prepared/lidc/{train,val}.jsonl`, `python -I`].
- A8's centred, shifted and rescaled train roots each hold 7,341 records and 23,469 targets, with 0 kept-side and 0
  fallback windows [measured: same].
- A1 step 5 / B.6.2: CR1 on R5 val = [0.58371, 0.65267] around 0.61819; patient-weighted mean 0.64722; 2-78 items
  per patient; top 5 patients 24.9 % of items; miss share 12.64 % [measured: `runs/r5/eval/eval_val/items.jsonl`
  joined with `group_id`]. The equal-N_p reduction of CR1 to the t interval of patient means holds algebraically.
- B1.0: the first s that meets (a) and (b) is 78, with 40 / 38 patients, 533 / 551 items, 24 / 20 disagreement
  items; the revision-4 rule gives 454 / 630 [measured: `(planning probe, not kept: b_r4)` on the val view].
- DG7-DG9 do hold by construction on every 20-step prefix under the B1.2 rounding: worst ratio >= 0.960, TV <= 0.006,
  offset KS p95 0.036 on one block, 0 of 200 simulated schedules above the 0.06 KS limit; all 16 cells of the train
  shifted law are >= 1 % (smallest 1.19 %) [measured: `(planning probe, not kept: r7data)` on `notes/b1_out/geom.json`;
  an approximation that draws each cell's offsets from its draw-0 windows].
- B.1 item 8 (masked truncation bias): E[1{not truncated}(R - b) grad log pi] = grad J - b grad P(not truncated),
  since R = 0 when truncated. Correct.

## Items

### 1. MAJOR. The `val_half` pytest sketch deals with `seed=0`; B1.0 and U02 deal with `"val_half:<s>"`. Seed 0 halves fail the B1.0 balance rules on the real val

Files: B (B1.0, B.9.1 U02, B.9.8 sketch `tests/test_reward_spec.py`).

Doc lines:
- B, "B1.0" Val halves: "deal the val patients with `splits.deal(groups, (("shape", 0.5), ("select", 0.5)),
  seed=f"val_half:{s}")` inside each `splits.lidc_stratum` stratum ... Try s = 0, 1, 2, ... and keep the first s
  that meets both: (a) |items(shape) - items(select)| <= 4 % of val items, and (b) |disagreement items(shape) -
  disagreement items(select)| <= 4 ... Today's val gives s = 78".
- B, "B.9.1" U02: "seeds `"val_half:<s>"` from s = 0 ... every s below the stored one fails (a) or (b) ... an unknown
  `group_id` raises `KeyError`."
- B, "B.9.8 Pytest sketches" header: "keep the assertions, change the call sites if the real seams differ". Sketch
  `test_val_half_is_the_stratified_deal`: `want |= splits.deal(gids, (("shape", 0.5), ("select", 0.5)), seed=0)`
  then `assert {g: val_half(g) for g in want} == want`.
- Same sketch, `test_val_half_independent_of_hashseed`: `print(val_half('lidc/LIDC-IDRI-0001'))` in three
  subprocesses, then `assert len(outs) == 1`.

Proof:
- Seed 0 on the real val (strata from `lidc_stratum`, the record-based `Patient` that B1.0's numbers use) gives
  40 / 38 patients, **567 / 517 items** (difference 50 > 43.4 = 4 % of 1,084) and **14 / 30 disagreement items**
  (difference 16 > 4) [measured: `python -I`, `splits.deal(..., seed=0)` per stratum on
  `data/eval_views/lidc-val-shifted/lidc/val.jsonl`, 2026-11-08]. It fails both (a) and (b). The B1.0 rule gives
  s = 78 and 533 / 551, 24 / 20 [measured: `(planning probe, not kept: b_r4)`].
- Following the header ("keep the assertions") forces `val_half` to be the seed-0 deal. Val-select then holds 30 of
  the 44 disagreement items and val-shape 14. B3.7 tunes on val-shape and Gate B3 judges on val-select (B.6.2), so
  the two halves differ in their hardest items. Following B1.0 and U02 fails the sketch. The two test specs give
  different `val_halves.json` files, so the engineer must guess.
- The sketch also calls `val_half` on synthetic ids from `tmp_store_val_patients`, while U02 says `val_half` is a
  lookup in `docs/reference/val_halves.json` and an unknown id raises `KeyError`. So the sketch cannot pass as
  written under either rule.
- The hash-seed sketch is vacuous. Group ids have no `lidc/` prefix (`"group_id": "LIDC-IDRI-0026"` [measured: val
  view record]), and LIDC-IDRI-0001 is a train patient (8 train items, 0 val items [measured: grep]). So all three
  subprocesses raise `KeyError`, stdout is empty in each, and `len(outs) == 1` passes.
- Related missing input: B1.0 never says which `Patient(slices, disagree)` feeds `lidc_stratum` for the halves. The
  split built it from the XML (`n + len(sl)` slices and `disagree_5px` flags) [code:
  experiments/phase1/build_data.py:122-127]. B1.0's measured s = 78 comes from per-patient record counts
  (`Patient(items, disagreement items)` in `valhalf.py`). The two sources can put patients in different strata.

Fix:
1. Replace the sketch body with the U02 rule: build a synthetic store, run `splits.val_half_table(records, max_s)`
   (a new pure function that `val_half`'s table is written from), and assert: the returned s is the first s that meets
   (a) and (b); every smaller s fails one of them; the table equals `deal(..., seed=f"val_half:{s}")` per stratum.
   Test `val_half` itself on the committed table only: one real val id, for example `LIDC-IDRI-0026`, returns
   "shape" or "select", and `val_half("LIDC-IDRI-0001")` raises `KeyError`.
2. Hash-seed sketch: print `val_half('LIDC-IDRI-0026')` (a val id, no prefix) and add `assert outs == {"shape\n"}`
   or `{"select\n"}`, the value from the table, so a `KeyError` fails the test.
3. B1.0: add "Patient for `lidc_stratum` = `Patient(slices=items of the patient in the frozen val view,
   disagree=its items with `disagreement` true)`, the form that gives s = 78 [measured: valhalf.py]". Record that
   form in `val_halves.json` beside `s`.

### 2. MINOR. I08 says `build_rl_rows` never opens val or test pixels; B1.4 step 4 runs every gate there, and DG2 and DG3 read val and test pixels

Files: B (B1.4 step 4, B1.3 `leakage` and `near_duplicates`, B.9.3 I08).

B1.4 step 4: "Run `quality.rows_leakage` and all other gates. Write the manifest only after `quality.failures(...)`
is empty". B1.3: `leakage(store_root)` hashes "every int16 slice ... cross-patient and cross-split" and
`near_duplicates` builds "the 32x32 int16 HU patch" of val and test items (DG2, DG3). I08: "`build_rl_rows` never
opens any val or test pixels, masks or view records (... only the frozen-hash file and val/test image hashes for
DG1b are read)". Reading val pixels for a leakage check is not leakage, and preflight runs all gates again, so the
outcome does not change; but the literal code fails I08. Fix: B1.4 step 4 runs the rows gates only (DG1b, DG4-DG14,
DG18); `audit_data.py` and preflight run the store gates DG1-DG3, DG15-DG17 too. I08 keeps its wording.

### 3. MINOR. The draw-0 equality tests ignore padded windows, which the pool drops

Files: B (B1.1 tests, B.9.3 I02). B1.1: "A window that extends beyond the slice is dropped and counted under
`"padded"`"; its tests and I02 say "draw 0 equals `build_eval_views` `(dx, dy, side)` for shifted on 200 val
items and for rescaled where the test view is not a kept-side window". `build_eval_views._record` keeps padded
windows [code: build_eval_views.py:55-62]. Padded draw-0 windows: val shifted 5 / 1,084, val rescaled 52 / 1,084,
train shifted 28, train rescaled 250 [measured: `padded` column of `notes/b1_out/geom.json`]. On 200 val items about
10 rescaled windows have no pool window, so I02 fails on correct code. Fix: I02 and the B1.1 test add "and the view
window is not padded; for a padded one the pool has no draw-0 window of that kind and `dropped["padded"]` counts
it".

### 4. MINOR. I20 omits the disagreement-reuse rule of B1.4 step 2

Files: B (B1.4 step 2, I20). B1.4 step 2: take the smallest `draws` that passes DG7-DG9, does not raise "and uses no
disagreement window more than once. If no value meets the last condition, take the smallest value that passes the
rest". I20: "picks the smallest value that passes DG7-DG9 on every block prefix and does not raise". Code that
ignores reuse passes I20 and may pick a smaller `draws` with up to 4 repeats. Fix: I20 adds a planted pool where
draws 1 passes with reuse 2 and draws 2 passes with reuse 1; the picked value is 2, and
`manifest.schedule.disagreement_max_reuse` is 1.

### 5. MINOR. ADR-0015 asks A2 for a "near-subset headroom"; A2 does not define it

Files: B (B.3.4 ADR-0015 "Confound", B.6.2 Gate B1 DiD rule), A (A2). ADR-0015: "A2 reports the near-subset headroom
= mean(fair doctor score - R5) on offset < 0.15 items, with its patient CI. If the headroom's upper bound is below
0.013 ... Gate B1 reports the rule as uninformative". A2's strata are `all`, `k>=3`, `k=2`, `found`, `miss`; A2 reads
no offsets, and "R5" is not said to be `model` (against k - 1 doctors) or `model_all_k`. The check changes only the
report wording, so minor. Fix: A2 `summarize` gains an optional `offsets` map (from the frozen view's
`view.offset_frac`) and a stratum `near` (offset < 0.15) with `headroom = doctor - model` (both against the same
k - 1 doctors), patient CI; ADR-0015 and B.6.2 cite that key.

### 6. MINOR. E-G10's trainer arm stays R5's old `eval_val/items.jsonl` even when B1.0 makes that file a report line

Files: E (E4 E-G10, GP5), B (B1.0 step 2). B1.0: on a content-hash difference "mark A3-val ... as the only Gate B1
baseline. R5's `eval_val/items.jsonl` is then a report line only". E-G10 always reads the trainer arm from
`runs/r5/eval/eval_val/items.jsonl`, while the V4 sampler arm runs on the frozen view. In that branch the TOST
compares two views. Fix: E-G10 and GP5 read A3's `eval_val/items.jsonl` (trainer path, `rl-full-3ep`, frozen view)
when it exists, and R5's file only when B1.0 recorded equal hashes; the report names the arm used.

### 7. MINOR. Reused disagreement windows are materialized once per scheduled row

Files: B (B1.4 step 3, B1.4 tests). B1.2 step 4 lets a disagreement window appear up to 4 times. B1.4 step 3
"Materialize only scheduled windows ... in a `multiprocessing.Pool(procs)`" with `window_id` as the item id, so a
reused window is prepared up to 4 times, possibly in parallel to one PNG and one cache folder. The test "every
window's PNG path and cache folder is unique" then fails or must be read per distinct id. Fix: "materialize each
distinct `window_id` once; `lidc/train.jsonl` repeats its record in schedule order"; the test reads "unique per
distinct `window_id`", and `files_count` counts distinct files.

### 8. MINOR (still open from round 6; minor items, not fixed). Listed so they are not lost

- r6 data-4: U09 "`dropped` has exactly the 3 int keys"; B1.1 defines 4.
- r6 data-5: B1.1 tests "S0 windows are >= 0.75 of off-centre windows" (the pool gives about 0.5; the rule is DG10
  on the schedule).
- r6 data-6: B.8.3 item 11 "older $2.20 means a 260-step B1 run"; front 6.2 gives 220. Item 14 and B.7 row "Gate B1
  fails" still say "decision table B.4.1" (front 1.1 holds it).
- r6 data-7: A3 step 5 "A lower correlation means a different view: stop"; the next bullet and B.6.2 say Gate B1
  still pairs with A3.
- r6 data-8: front 1.1 has no precedence when two rows match.
- r6 data-11: B1.9 "The decision interval of the simulation is the CR1 interval" (Gate B1 uses the percentile
  interval at 78 patients).
- r6 data-12: B1.5 "Add columns `kind` and `cell` ... when the record has `view`" (eval-view records have `view`
  without `kind`).

## Counts

Blocking 0. Major 1 (item 1). Minor 7 (items 2-8).
