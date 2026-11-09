# Final review round 9: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).
Scope: the whole doc, read by section from the sources. I read all of these: A.1, A.3.2, A1-A5, A8, A.6-A.8, A.9.1 (U01-U18),
A.9.6 (G04, G05); B.1, B.3.1-B.3.4, B1.0-B1.9, B.6.1, B.6.2, B.7, B.8, B.9.1 (U01-U21) and the B.9.8 sketch
`tests/test_data_windows.py`. From the front I read 0, 1.1, 3, 6.1 and 6.2. From the back I read "Data and quality". From DOC I read DOC6 and DOC11 (card).
I checked the round-8 fixes in `review-log.md` (data-1: the val legacy reference 0.723, the found-stratum line 0.734 / 0.740, ADR-0015 0.731).
They are in place in A, B, DOC and 99_back, so I do not repeat them. Scripts for this round are in `(planning probe, not kept: r9data)`
(`ci_check.py`, `subset_patients.py`, `pool_sim.py`, `capacity.py`, `capacity2.py`). All of them ran with `.venv/bin/python -I`
on 2026-11-09, and each script reads the repo path as an argument.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Nothing in the data scope needs an NVIDIA product choice. B1.7's DICOM helper stays "recommend only". I ran no catalog lookup and installed nothing. Nothing changed. |
| engineering:code-review | I read the code behind each data claim: `signal.draws/_ci/ci`, `splits.{deal,lidc_stratum,lidc_split}`, `dataset._prepare/prompts/records/_verified`, `augment.{row_rng,shift,draw,neighbour_side}`, `frame.{union_box,centred,frame,policy_side}`, `build_eval_views._record`, `policy.chain` paths, and `huggingface_hub` `upload_folder` (2.2.0). This confirmed that `rl` is `rl_eligible(raters)` and does not depend on the window masks, so a rescaled window can never drop out of `prompts()`. |
| engineering:debug | I reproduced each number before I wrote it: patient CIs and CR1 on val and test, the near/far subset patient counts, a full window pool for draws 1-8 on the train store, and per-cell capacity against the 260-step schedule (items 3 and 4). |
| engineering:system-design | I traced where each published interval comes from: A1 schema 2, then `final_eval.json`, then DOC6 `headline.py`, then `headline.json`, then README and HF card (item 1). |
| engineering:architecture | I checked who owns each rule. The DOC owner table gives A the job of moving the headline to `ci_unit: patient`, but no A task does it, and A1's schema change breaks the DOC6 reader (item 1). |
| engineering:documentation | Each item quotes the doc line and the code line or measurement that proves it. |
| huggingface-skills:huggingface-best | No model choice is in scope. Front section 3 stays as it is. |
| huggingface-skills:huggingface-datasets | I checked the rows-root Hub path: `upload_folder` splits large folders into " (part N)" commits and returns the last commit (`_upload_pipeline.py:476, 679-692`), so B1.4's `CommitInfo.oid` rule is right. No finding. |

Web checks (2026-11-09):
- TCIA LIDC-IDRI page. 0396 SOP UID fix of 2018-06-28: confirmed. 8 two-timepoint patients: confirmed, and `splits.py` keeps them together. CC BY 3.0: confirmed.
- TCIA 044.xml / 191.xml errata. I read the local `tcia-lidc-xml/187/044.xml` comment: it concerns missing characteristic ratings, not contours, so it has no effect on masks.
- YongchengYAO/LIDC-IDRI-Lite duplicate-z list (0085, 0146, 0418, 0572, 0979): confirmed.

## Claims verified (no item)

- **Patient CI on test.** [0.6252288200704239, 0.6644106494958667] (B 2,000, seed 0, 158 patients): reproduced exactly. Test CR1 is
  [0.6251, 0.6645].
- **Val, CR1 and percentile.** Val CR1 [0.5837, 0.6527] around 0.6182, and val percentile [0.5851, 0.6526]. The patient-weighted mean is 0.6472.
  Items per patient run from 2 to 78, and the top 5 patients hold 24.9 %. These match A1 step 5 and B.6.2 [measured: `r9data/ci_check.py`].
- **Counts.** Val has 200 nodules, 78 series and 78 patients. Test has 401 nodules, 160 series and 158 patients. Miss shares are 12.64 % (val) and 10.82 % (test).
- **Patient overlap.** Overlap across the store splits is 0 / 0 / 0, with 557 / 78 / 158 patients.
- **Train items.** All 7,341 train items are `rl` (rater counts 3,379 x 4, 2,029 x 3, 1,933 x 2), and 289 of them are disagreement items. The 23,469 targets equal
  sum(k).
- **Rescaled windows.** 62.0 % of the draw-0 rescaled windows have an apparent size under 10 px (7,341 windows; 61.1 % of the 7,091 non-padded ones).
  This matches ADR-0015 [measured: `r9data/pool_sim.py`].

## Items

### 1. MAJOR. The plan's only way to publish a patient-CI headline cannot run. DOC6 reads the bare `ci95` that A1 removes, and a run layout that A3 does not write

Files: DOC (DOC6, the "Docs per phase" owner table, DOC11 `card_sections`), A (A1 step 3, A3 Files, A9).

Doc lines:
- DOC, "DOC6": "`headline.py` ... Reads `runs/<ID>/eval/eval_test/final_eval.json` and
  `runs/<ID>/eval/eval_val/final_eval.json` ... Keeps the 6 metrics with `mean` and `ci95` unrounded." and
  "`ci_unit` is `"item"` because `final_eval` uses an item bootstrap ...; phase A writes `"patient"` when it switches
  the gate CI to the patient-cluster bootstrap ... [A.md, "Words", gate CI]". A.md has no "Words" section.
- DOC, owner table "Docs per phase", row A: "`results/headline.json` via `make headline` (`ci_unit: patient`)".
- DOC, "DOC11": `card_sections` writes "the 6 test metrics with 95 % CI, view, items, K, seed, `ci_unit`" from `headline.json`
  into the HF card.
- A, "A1" step 3: "The bare key `ci95` is dropped in schema 2, so an old reader fails loudly ... Only `final_eval.py`
  and `tests/test_final_eval.py` read `ci95` today". A, "A.1" item 1: "The intervals are too narrow."
- A, "A3" Files and plan files: `--out {rd}/eval_val` (V1, `a3_measure.yaml`) and `--out {rd}/eval_test` (V3,
  `a3_test.yaml`). A9 lists no headline step.

Proof:
- After A1, every new `final_eval.json` is schema 2. Its metrics hold `ci95_nodule/ci95_patient/ci95_item` and no `ci95`
  (A1 step 3, U08 "no key `ci95` anywhere"). DOC I04 makes `headline.main` refuse a file with a missing metric key. So
  `make headline` on any A-era run exits non-zero.
- The paths do not match. R5's files are at `runs/r5/eval/eval_{val,test}/final_eval.json` [measured: ls]. A3 writes
  `runs/<run>/eval_val` in the V1 run and `runs/<run>/eval_test` in the V3 run [code: `policy.chain` puts `{rd}` =
  `/workspace/runs/<run>`, policy.py:433; `pull.sh` copies `runs/<run>/`]. So no single `RUN=<id>` has
  `eval/eval_test` and `eval/eval_val`.
- Effect when following the doc: `headline.json`, the README check and the HF card keep the R5 item CI
  [0.6334, 0.6562] (width 0.023) as "95 % CI". A1 calls that interval too narrow, and the gate CI is
  [0.6252, 0.6644] (width 0.039) [measured: `r9data/ci_check.py`, reproduced to the digit]. The A row of the DOC table
  cannot be done, and no task says how `headline.py` reads schema 2. The engineer must guess which key to copy and
  which run folders to read.

Fix:
1. DOC6 `headline.py`: replace `--run ID` with `--test PATH --val PATH` (two `final_eval.json` files, each under
   `runs/`), and keep `--run ID` as a short form for `runs/<ID>/eval/eval_{test,val}`. Accept `schema` 1 and 2.
   Schema 1 gives `ci95 := ci95`, `ci_unit := "item"`. Schema 2 gives `ci95 := metrics[m][gate_ci]` and `ci_unit := gate_ci`
   without the `ci95_` prefix (`"patient"`). Refuse a test/val pair with different schemas. Add I03/I04 cases for
   schema 2.
2. A1 step 3: add `final_eval --resummarize DIR`. It reads `DIR/items.jsonl`, joins `group_id` from the store,
   writes a schema-2 `final_eval.json` beside it, and loads no model and does no scoring. Run it on `runs/r5/eval/eval_{val,test}` after A1. Then
   `make headline RUN=r5` publishes the patient CI [0.6252, 0.6644] with `ci_unit: patient`. Update A1's grep
   sentence to name `experiments/phase2/headline.py` when DOC6 has landed first (DOC6 may run before A, DOC "any
   order after DOC1").
3. A9: add "`docs/results/headline.json` regenerated with `ci_unit: patient` (DOC6)". In DOC6, change "[A.md,
   "Words", gate CI]" to "[A.md, A1 step 5]".

### 2. MINOR. `BudgetStop` can end the B1 run inside a block, but B1.2 and B.7 say only whole blocks are trained

Files: B (B1.2, B1.5 "Run guard", B.7 row "A budget cut ..."), 00_front (6.2 notes).

Doc lines: B1.2: "So every prefix of whole blocks has exact per-kind cell counts". B1.5: "Also require `--steps <=
schedule.steps` and `--steps % schedule.block_steps == 0` (only whole blocks are trained)". B1.5 "Run guard": "it
stops training and saves when fixed_usd + (now - session_start) x rate reaches the cap". B.7: "A budget cut trains a
prefix that misses the target | ... | 20-step blocks; only whole blocks are trained". Front 6.2 notes: "With the
53-min start the stop trips only above about 27.5 s per step ... If it trips, the adapter is saved".

Proof: at 28.5 s per step with the 53-min start, the stop trips at (172.3 - 53) x 60 / 28.5 = 251 steps, which is 11 of
20 steps of block 13 [estimate: the front's own formula]. The rows-root download is "not in the measured rows"
(front 6.2), so a later start is plausible. The block shuffle is uniform, so the partial block is near the target and no
gate result changes. The doc's claim and the B.7 mitigation do not hold on this path. Fix: `BudgetStop` decides only
at block ends. At step s with `s % block_steps == 0`, it stops if spend at s plus `block_steps` x (mean step time so
far) x rate would pass the cap. `rl.json` records `trained_steps`. If a stop still lands inside a block (a crash),
the Gate B1 report runs DG7-DG9 on the trained prefix and prints the result.

### 3. MINOR. Gate B1's near rule sits at 52 patients, two above the 50-patient switch of the interval type

Files: B (B.6.2 "Pairing and CI", Gate B1 rule 3), A (A1 step 5).

B.6.2: "For a subset rule ... apply the same choice to the per-item differences of the subset, with P = the patients
that have items in it". [measured: `r9data/subset_patients.py`, union-box offset as `r8data/near_ceiling.py`]: the
near subset (offset < 0.15) has 206 items from **52** patients. Without the two B1.7 val patients (0085, 0146) it has 203
items from **51** patients. The far subset has 406 items from 75 patients (399 / 73). If B1.7 adds one more val hit with
near items, rule 3 switches from the percentile interval to `cluster_t`. A label-check outcome then changes the
interval method of one rule but not the others, while the difference-in-differences rule uses the whole-val draws. Fix: B.6.2
"Gate B1 rules (1)-(4) all use the percentile patient interval of the whole-val draws (78 patients), whatever the
subset's P. The `cluster_t` line is printed beside each subset rule". B1.9 and U49 use the same wording.

### 4. MINOR. B1.4 step 2 says the feasible `draws` "is not yet measured"; it can be measured now

Files: B (B1.4 step 2, B.7 row "`draws` needed for DG7 is large").

B1.4 step 2: "(Today one draw fails [measured: SD3 s5]; the value that passes is not yet measured.)". I built the pool of
B1.1 (geometry only, draw-0 rngs as `build_eval_views`, `side=S0` for shifted, padded windows dropped) for draws 1-8.
I checked the per-(kind, cell) capacity against the 260-step need (13 blocks; 96 / 480 / 64 kind slots; target =
draw-0 shifted law; B1.2 step 3 quota per block, w = 4) [measured: `r9data/pool_sim.py` (13 s on 6 processes) and
`capacity.py` / `capacity2.py`, 2026-11-09]:
- With reuse <= 4: draws = 1 fails one set, (shifted, cell (1, 1)): 13 disagreement slots against 3 windows x 4 = 12.
  draws = 2, 4 and 8 pass.
- With no disagreement reuse: draws = 1, 2 and 4 fail (17, 12 and 1 sets); draws = 8 passes.
- Draw-0 drops: 28 padded shifted windows and 250 padded rescaled windows. No `shifted_none` or `rescaled_none`.

So under B1.4's rule ("smallest value that ... uses no disagreement window more than once"), the build takes
draws = 8 (if the full scheduler, with its one-item-per-step rule, agrees). This is a capacity bound only. Fix: replace the
parenthesis with these numbers. Add "the B1.4 run confirms them with the real `schedule`". The B.7 row then reads
"draws 8 expected; only scheduled windows are materialized".

### 5. MINOR (still open from earlier rounds; not re-argued). Listed so they are not lost

- r4 data-7 / r6 data-5: B1.1 tests "S0 windows are >= 0.75 of off-centre windows" on the pool. The real draw-0 pool
  gives 0.508 [measured: `r9data/pool_sim.py`], so the test fails on correct code. DG10 is a schedule rule (0.75 / 0.85
  = 0.88 there). Move the assertion to the schedule test (R03).
- r6 data-4: U09 still says "`dropped` has exactly the 3 int keys". B1.1 defines 4 (`padded` included).
- r5 data-4 / r6 data-11: B1.9 still says "The decision interval of the simulation is the CR1 interval of the
  item-weighted mean (B.6.2)". B.6.2 and U49 use the percentile interval at >= 50 patients.
- r5 data-7: DG15 "the val shifted entry equals the B1.0 pre-rebuild hash" cannot pass on B1.0 step 2's difference
  branch ("freeze the rebuilt view").
- r5 data-8: B1.3 `tokens()` counts "certified target tokens (+1 end token)" and then gates "`max + 1 < limit - 100`".
  The +1 is counted twice (A8 step 5: "no '+1' is added again"). The DG12 row "max 1,181 + 1" is right.
- r6 data-7: A3 step 5 still says "A lower correlation means a different view: stop". The next bullet and B.6.2 say
  that a failure drops only the R5 line.
- r7 data-2..5 and r8 data-2..4 (wording and test-text items): unchanged.

## Counts

Blocking 0. Major 1 (item 1). Minor 4 (items 2-5).
