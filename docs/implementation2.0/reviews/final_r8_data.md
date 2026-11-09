# Final review round 8: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).
Scope: the whole doc, read by section from the sources. Read in full: A.1, A1, A2, A3, A8, A.6-A.8, A.9.6 (G04, G05);
B.1, B.3.4 (ADR-0015), B.4, B1.0-B1.9, B.6.1, B.6.2, the B.9.8 sketch `tests/test_reward_spec.py`; front 1.1, 3,
6.1, 6.2; 99_back "Data and quality". Round-7 fixes in `review-log.md` (data-1: val-halves `Patient` form,
`val_half_table`, new sketches) are not repeated. Scripts for this round: `docs/implementation2.0/notes/r8data/{legacy.py,
legacy_found.py, near_ceiling.py}` (all run with `.venv/bin/python -I`, 2026-11-08).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in the data scope (B1.7's DICOM helper stays "recommend only"). No catalog lookup, no install. Nothing changed. |
| engineering:code-review | Read the code behind each data claim: `splits.{deal,lidc_stratum,lidc_split}`, `build_data.py:123-127`, `dataset.{_prepare,records,_verified,prompts,grpo_rows}`, `augment.{shift,draw,neighbour_side}`, `build_eval_views._record`, `final_eval.main`, `cache.{voters,rl_eligible}`. Confirmed that repeated rows survive `records`/`prompts` (no de-duplication), so the schedule order holds. |
| engineering:debug | Reproduced before writing: legacy LOO d_max ceiling on the store, both pooling orders, all and found items (item 1); near/far subset ceiling with the B.1 offset definition (206 / 406 items, 27.34 % far misses, matching B.1 item 2); rater counts 281 / 330 / 473 (val), 633 / 502 / 1,065 (test); B1.7 sidecar hits (5 of 999); `-I` and `PYTHONHASHSEED` (item 3). |
| engineering:system-design | Traced which number each gate reads: A-G2 / G04 read the whole-val legacy ceiling; the 0.734 reference comes from the found-image row of Research2.0 (item 1). |
| engineering:architecture | Checked rule ownership: A2 step 4, A-G2, G04 and 99_back decision 12 all state the same 0.734 reference, so one wrong number propagates to four places; ADR-0015 cites it as the near-item bound. |
| engineering:documentation | One item per defect, each with the doc line and the code line or measurement. |
| huggingface-skills:huggingface-best | No model choice in scope; front section 3 left as is. |
| huggingface-skills:huggingface-datasets | No Hub dataset is read by A or B1 (LIDC comes from TCIA into the local store). No finding. |

Web checks: none needed for the items below; every claim is checked against the repo code and local data. The
B1.7 TCIA and third-party references were checked in round 7 and are unchanged.

## Claims verified (no item)

- Val rater counts k = 2 / 3 / 4: 281 / 330 / 473 items (A2's 281 / 803); test 633 / 502 / 1,065 [measured:
  `store.load("lidc", split)`].
- B.1 item 2 val subsets with the union-box offset: near (< 0.15) 206 items, R5 IoU 0.7511, misses 0.49 %; far
  (> 0.3) 406 items, IoU 0.5032, misses 27.34 % [measured: `r8data/near_ceiling.py`].
- Test legacy ceiling (item pooling) 0.7319, so the test half of A-G2 (0.732 ± 0.005) passes; test found-image
  ceiling 0.7403 (1,962 items) and model 0.7219 match B.1 item 7 and front 3 ("0.722 against 0.740")
  [measured: `r8data/legacy_found.py`].
- B1.7: 5 of 999 sidecars hit (d), the five series and counts as in the B1.7 table [measured: `data/lidc_slices`].
- B1.4 / B1.5: `dataset.records` and `prompts` keep repeated records in file order (no dict by id) [code:
  src/objdel/train/dataset.py:159-161, 217-232], so a reused disagreement window stays in schedule order and the
  `len(rows)` / `ids_sha256` checks hold.

## Items

### 1. BLOCKING. A-G2 / G04 compare the whole-val legacy ceiling with 0.734, which is the found-image value. A correct A2 fails its gate and `make test-data` fails

Files: A (A.1 item 2, A2 step 4, A.6 A-G2, A.9.6 G04, A.6 row "All"), 99_back (decision 12), B (B.3.4 ADR-0015
"Confound").

Doc lines:
- A, "A2" step 4: "Legacy check (A-G2): if neither `legacy_loo_dmax_item` nor `legacy_loo_dmax_pooled` is within
  ±0.005 of 0.734 (val) / 0.732 (test), the gate fails and A2 stops until the owner accepts the difference."
- A, "A2" Definition: `legacy_loo_dmax_item`: "per item, the mean over voting r of the max plain IoU of rater r with
  any other rater; then the mean over items"; `legacy_loo_dmax_pooled`: "the same max IoU values pooled over all
  (item, r) pairs, then one mean". Strata `all` on val = 1,084 items.
- A, "A.9.6" G04 (marker `data`): "`min(abs(item - 0.734), abs(pooled - 0.734)) <= 0.005` on val and `0.732` on
  test; otherwise fail".
- A, "A.6" row All: "`make test-data` passes on the Mac with `data/`"; A.9.0: `test-data` runs `-m "data and not
  slow"`.
- 99_back "Data and quality" 12: "Neither legacy line is within ±0.005 of 0.734 (val) / 0.732 (test)."

Proof:
- Whole val, the A2 definition (store masks, all 1,084 items): item pooling **0.7232**, pair pooling **0.7419**;
  with voting raters only (signature with `min_area`, `cache.voters`): 0.7232 / 0.7419 [measured:
  `docs/implementation2.0/notes/r8data/legacy.py`, 2026-11-08]. |0.7232 - 0.734| = 0.0108 and |0.7419 - 0.734| = 0.0079, both above
  0.005. So G04 and A-G2 fail on correct code.
- The 0.734 is the found-image value: val items with K-mean IoU_tc >= 0.1 (947 of 1,084) give item pooling
  **0.7339** (model 0.7056) [measured: `docs/implementation2.0/notes/r8data/legacy_found.py` joined with
  `runs/r5/eval/eval_val/items.jsonl`]. That is the pair "0.706 against a ceiling of 0.734" of
  `docs/implementation2.0/notes/C_rl_quality.md` section 1. Its offset table weights to (255 x 0.732 + 173 x 0.719 + 417 x 0.727 + 239 x
  0.710) / 1,084 = 0.7232 for all items, the same as my value. The source of the error is
  `docs/research/research2.0/Research2.0.md:196-197`, which prints 0.734 for both "All images" and "found".
- The test reference is the all-item value (0.7319 measured, 0.732 in the doc), so val and test use two different
  estimands in one rule.
- Effect when following the doc: A2 stops for an owner review caused by a wrong reference, and G04 fails inside
  `make test-data` for good, so the A.6 "All" gate can never pass as written.

Fix:
1. A2 step 4, A-G2, G04 and 99_back 12: val reference = **0.723** (whole val, item pooling; measured 0.7232),
   test reference 0.732 (unchanged). Write: "both references are the item-pooled legacy value on all items;
   measured 0.7232 (val) and 0.7319 (test) on the store masks, 2026-11-08 [`r8data/legacy.py`]".
2. Add a second G04 line after A3 (found stratum, K-mean rule on R5's items): `legacy_loo_dmax_item` on found
   items within ±0.005 of 0.734 (val, 947 items) and 0.740 (test, 1,962 items).
3. A.1 item 2: "0.732 (test) and 0.723 (val; 0.734 on found val images)".
4. B.3.4 ADR-0015 "Confound": replace "above the legacy val doctor ceiling 0.734" with "above the legacy doctor
   ceiling of the same near items, 0.731 (206 items; whole val 0.723) [measured: r8data/near_ceiling.py]". The
   conclusion of that sentence does not change.
5. DOC phase: correct the "All images" row of `Research2.0.md:196` to 0.723 (val), so no later reader copies 0.734.

### 2. MINOR. The hash-seed sketch sets `PYTHONHASHSEED` but runs `python -I`, which ignores it

Files: B (B.9.8 sketch `tests/test_reward_spec.py`, `test_val_half_independent_of_hashseed`).

The sketch runs `[sys.executable, "-I", "-c", code]` with `env={"PYTHONHASHSEED": s}` for s in 0, 1, 2. `-I`
implies `-E`, which ignores every `PYTHON*` variable. Measured: `env -i PYTHONHASHSEED=0 python3 -I -c
"print(hash('abc'))"` printed three different values in three runs; without `-I` the value is the same each time
[measured, 2026-11-08]. The test still passes on correct code (each run gets a random seed), so no outcome changes,
but it does not test the seeds it names. The script reads only `docs/reference/val_halves.json`, not `data/`. Fix:
run `[sys.executable, "-s", "-P", "-c", code]` (no user site, no script dir on the path; `PYTHONHASHSEED` is then
honoured), or keep `-I` and drop the `env` loop, with the comment "three random hash seeds".

### 3. MINOR. B1.2 step 3 "largest-remainder rounding of q x slots" per single slot set does not say what is rounded together

Files: B (B1.2 step 3, DG11, DG14).

B1.2 step 3: "Per block and per (kind, cell) slot set ..., the disagreement slot count is the largest-remainder
rounding of `q x slots`". Largest remainder needs a set of quantities with a fixed total; one slot set has one
number. With 16 cells at about 6 shifted slots per block each (smallest cell 1.19 %, 480 shifted slots), plain
rounding of q x slots per cell can move the block's disagreement total by several rows from sum(q x slots), so
DG14's "14.1 % overall" is not what B1.2 says it keeps "in every block". Two engineers can write two DG11 oracles.
Fix: "For each kind in a block, the disagreement counts of its (kind, cell) slot sets are the largest-remainder
rounding of the vector (q_c x slots_c) over the kind's cells, so the kind's block total is round(sum_c q_c x
slots_c)". DG11 and P06 assert that vector rule.

### 4. MINOR. A.8 says B1 audits A8's shifted and rescaled train views; B1.4 audits the new rows root

Files: A (A.8 Dependencies, A8 step 5), B (B1.3 `tokens`, B1.4).

A.8: "B1 uses ... the A8 length rows of the shifted and rescaled views (same manifest SHAs); it audits those
views." A8 step 5: "B1 then reuses these views and only audits them." B1.4 builds `data/views/lidc-train-b1` with
draws d >= 1 windows that are not in `lidc-train-{shifted,rescaled}`, and DG12 runs `quality.tokens(rows_root,
...)` on that root. No outcome changes (DG12 reads the right root). Fix: A.8 and A8 step 5: "B1's DG12 counts
the tokens of the rows root itself; A8's view rows are the pre-B1 record".

### 5. MINOR (still open from rounds 6 and 7; not re-argued). Listed so they are not lost

- r6 data-4: U09 "`dropped` has exactly the 3 int keys"; B1.1 defines 4.
- r6 data-5: B1.1 tests "S0 windows are >= 0.75 of off-centre windows" (the pool gives about 0.5; DG10 is a
  schedule rule).
- r6 data-11: B1.9 "The decision interval of the simulation is the CR1 interval" (Gate B1 on 78 patients uses the
  percentile interval).
- r7 data-2: I08 "`build_rl_rows` never opens any val or test pixels" vs B1.4 step 4 "all other gates" (DG2, DG3
  read val and test pixels).
- r7 data-3: the draw-0 equality tests ignore padded windows, which the pool drops.
- r7 data-4: I20 omits the disagreement-reuse rule of B1.4 step 2.
- r7 data-5: ADR-0015's "near-subset headroom" has no A2 output key.

## Counts

Blocking 1 (item 1). Major 0. Minor 4 (items 2-5).
