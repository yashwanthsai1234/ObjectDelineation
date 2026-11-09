# Final review round 12: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).

Scope. I read the sources by section, not the built file:
- A: header (incl. the round-11 spec rule), A.1, A1 (steps 1-6, 5b, output keys), A2, A3 (files, plans, steps 0-8),
  A4, A5, A8, A.6.
- B: header, B.1, B.2, B.3.1-B.3.4 (ADR-0015, ADR-0013, ADR-0014), B.4 task order, B1.0-B1.9, B3.7, B.6.1, B.6.2,
  the DG3 / DG17 test rows (U24, I17, I24, I30), B.9.0, B.9.12, the B round 7-11 logs.
- Front: 1.1, 6, 6.1 (steps 0-6), 6.2 (whole), 8. 99_back "Data and quality". DOC6 headline schema handling.
- `phases/B.spec.md` header and §§1-5 (Gate B1 text).
- `review-log.md`: rounds 8-11 for A and B, so handled items are not repeated.

Scripts (all `.venv/bin/python -I`, repo path as argument, 2026-11-09): `(planning probe, not kept: r12data)`, plus
the inline `-c` checks quoted below. I reused `(planning probe, not kept: r11data)` (union-box offsets of the 1,084 val
items, B.1 item 2 definition).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in the data scope. B1.7's DICOM helper stays "recommend only". No catalog lookup, nothing installed. |
| engineering:code-review | Read the code behind the data claims: `augment.{row_rng,shift,draw,neighbour_side}`, `build_eval_views._record/build`, `dataset.{_prepare,prompts,grpo_rows}`, `cache.rl_eligible`, `spec.fingerprint`, `lidc.series_files` (the 0101 correction file), the `b1_pixels.py` seed (DG3 score). |
| engineering:debug | Reproduced each number I lean on before using it (train 7,341 items / 557 patients / 23,469 targets / 289 disagreement items, val strata 281 / 803, the 27 train items in the 4 B1.7 hit series, the disagreement flag on the two existing train views). |
| engineering:system-design | Traced the data path store -> `pool` -> `schedule` -> `window_record` -> rows root -> `quality` gates -> `train_rl`, and asked at each hop which input can remove an item. That found item 1 (no hop removes a "dropped" DG17 item). |
| engineering:architecture | Checked which file owns each rule: B.md versus `B.spec.md` (the round-11 precedence rule landed in A only, item 4); A3's header versus front 6.1 step 4c (item 2). |
| engineering:documentation | Each item quotes the source line and the code line or measurement that proves it. |
| huggingface-skills:huggingface-best | No model choice in scope. Front section 3 unchanged. |
| huggingface-skills:huggingface-datasets | No Hub dataset is read on the data path. The rows root goes to the private model repo through `HfApi.upload_folder` (B1.4 step 5), not through the Dataset Viewer, so no Viewer endpoint applies. |

Web checks (2026-11-09):
- TCIA LIDC-IDRI page, "Annotation and Markup Issues" [web: https://www.cancerimagingarchive.net/collection/lidc-idri/]:
  8 items. Item 4 (LIDC-IDRI-0101 XML replaced, 2012-03-21) is handled: `data/raw/lidc_xml` holds
  `161-resubmitted-correction-3-9-12.xml`, and `lidc.series_files` takes the file whose name holds "correction"
  [code: src/objdel/data/lidc.py:263-269]. Item 6 (0396) is B1.7 check (c). Item 7 (0510, `internalStructure` = 5)
  does not touch outlines. Item 8 (8 patients with two timepoints) is safe because the split key is the patient
  id (`group_id`), so both scans of one patient fall in one split (B.1 item 4: patient overlap 0). LIDC-IDRI-0101 and
  -0510 are train patients, -0396 is not in the store [measured: store `group_id`].
- A search for "same subject under two LIDC patient ids" found no source; nothing in the plan depends on it.

## Claims verified (no item)

- Train store: 7,341 items, all with >= 2 doctors (`rl` true for 7,341 of 7,341), 557 patients, 289 disagreement
  items, 23,469 non-empty targets [measured: `data/prepared/lidc/train.jsonl`, `data/store/lidc/train.parquet`].
  So `pool`'s `"not_rl"` count is 0 on today's store, and 7,341 - 289 + 4 x 289 = 8,208 holds.
- Disagreement per nodule (B1.1): 1,247 train nodules, 100 with a disagreement slice, 0 of them all-disagreement,
  1,394 items in those nodules [measured: inline check]. Matches B1.1.
- The `disagreement` flag that `_prepare` recomputes on a window equals the source item's flag on every one of the
  7,341 records of `lidc-train-shifted` and of `lidc-train-rescaled` (0 differ) [measured: inline check]. So B1.1's
  "every window copies it" is consistent with what `window_record` will write.
- Val: 1,084 items, k = 2 / 3 / 4 doctors on 281 / 330 / 473 items (A2's 281 / 803 strata), 44 disagreement items
  [measured: `data/eval_views/lidc-val-shifted/lidc/val.jsonl`].
- A3 step 7 / B.6.2 noise floor: `d_i` = mean(answers 0-3) - mean(answers 4-7) has variance sigma^2 / 2, and the
  difference of two independent K = 8 means has variance sigma^2 / 4. So `sd(d_i) / sqrt 2` is the SD of a
  two-sampling difference, directly comparable with the measured 0.0569 pair. The formula is right (see item 3 for
  the label). Spearman-Brown `r_sh = 2 r_h / (1 + r_h)` is the right K = 8 reliability for two K = 4 halves.
- B1.9 arithmetic: 0.077 x 0.0422 = 0.00325 (measured 0.00324); per-arm 0.00162; 1 / 0.077 = 13.
- DG3 score: `b1_pixels.py:39-41, 67` subtracts the crop mean and divides by its norm (a Pearson correlation), as
  B1.3 says. Its null rate is at the cross-split rate (within-train, other patients: 4 pairs >= 0.98 over 7,341
  items; test -> train: 1 over 2,200) [measured: `docs/implementation2.0/notes/b1_out/neardup.json`], so DG3 shows no
  sign of a duplicated patient across splits.

## Items

### 1. MAJOR. DG17 counts a train hit as "resolved" from the report alone. No step removes a "dropped" train item from the pool or the rows root, and no test checks it

Files: B (B1.1 `pool`, B1.4, B1.7 "Gate DG17", B.6.1 DG17 and DG1b, B.9 I17 and I24).

Doc lines:
- B.md, B1.7: "Gate DG17: 0 unresolved hits. A train hit is resolved by mapping by `imageSOP_UID` or by dropping the
  items from train; the report says which."
- B.md, B.6.1 DG17: "unresolved train hits ... | 0".
- B.md, B1.1: "`pool(items, draws: int, seed: int = 0) -> WindowPool`" and "`dropped: Mapping[str, int]` (keys
  `"not_rl"`, `"shifted_none"`, `"rescaled_none"`, `"padded"`)". The pool has no exclusion input.
- B.md, B1.4: `build_rl_rows.py --store data/store ...`, step 1 "`pool(train items, draws=n)`". No exclusion flag.
- B.md, B.6.1 DG1b checks only that each `source_item_id` is a train item, its `group_id` is not val or test, and no
  image hash equals a val or test image.
- B.md, I17: "DG17 fails on an unresolved hit and passes when the report marks it resolved". I24: "items from an
  affected series are dropped or mapped as the decision says" (a fixture of `check_series_z`; it names no rows root).
- Contrast DG3, where the same word is backed by a test: U24 "marked 'drop train item' (and the train item absent from
  the pool)" and I30 "marking it 'drop train item' removes the train item from the pool". DG17 has no such case.

Proof:
- The four train hit series hold 27 train items in today's store: LIDC-IDRI-0013 (15), -0017 (5), -0365 (4), -0566
  (3) [measured: `data/store/lidc/train.parquet`, series suffixes 380786, 228428, 179144, 295628, 2026-11-09]. Steps
  (a) and (b) of B1.7 have not run yet, so the count can grow.
- An engineer who follows the doc writes "dropped" in the B1.7 report. DG17 then reads 0 unresolved and passes; I17 and
  I24 pass; `pool` still reads all 7,341 train items from `data/store`, and the scheduler can put windows of those
  items in the rows root. Every planned $0 check passes. The owner's first rule ("dataset quality is the highest
  priority") is then met on paper only: suspect labels are trained in V1. The effect on Gate B1 is small (27 of
  7,341 items, 0.37 %), so this is major (a gate that does not measure what it claims, and a missing step), not
  blocking.
- The "mapping" branch has the same gap: it changes `lidc.slice_index` or the store build, but no step says the store
  and the prepared roots are rebuilt and B1.7 is re-run before B1.4.

Fix:
1. One exclusion list, `docs/reference/train_exclusions.json`: `{item_id: {"reason": "DG17 <series>" | "DG3 <pair>",
   "date": ...}}`. B1.7 step 2 writes the DG17 drops into it; DG3's "drop train item" writes into the same file.
2. B1.1: `pool(items, draws, seed=0, exclude=frozenset())`; an excluded item gives no window and is counted under
   `dropped["excluded"]` (five keys; U09 follows). B1.4: `build_rl_rows --exclude docs/reference/train_exclusions.json`
   (required flag; an empty list is written as `{}`); the manifest `schedule` records the list's SHA-256. The target
   law is then computed without the excluded items; say so in B1.4 step 1.
3. DG17 rule: a train hit counts as resolved only if (mapped) the store and prepared roots were rebuilt after the fix
   and a re-run of `check_series_z` lists no hit for that series, or (dropped) every store item of the series is in
   `train_exclusions.json`. DG1b gains: "no `source_item_id` of the rows root is in `train_exclusions.json`".
4. Tests: I17 adds "a hit marked dropped whose items are not in `train_exclusions.json` fails DG17"; I24 adds "a rows
   root built from a store with a planted hit series holds no window of it when the series is dropped"; E04 adds a
   planted excluded item in the rows root (DG1b fails).

### 2. MINOR. A3's header comment bundles with the test view still in `data/eval_views` and says it "reaches the box"; front 6.1 step 4c requires no `lidc-test-*` folder in the V1 bundle

Files: A (A3, the paragraph before the `a3_measure.yaml` header comments, and the header comment block), front (6.1
step 4c).

Doc lines:
- A.md, A3: "So `lidc-val-edge` (built by B1.6) reaches the V1 box, and `lidc-test-shifted` reaches the box when it
  exists on the Mac at bundle time." The header block then runs `build_eval_views --split test ...`, `final_eval
  --hash-view .../lidc-test-shifted`, and `vast.bundle lidc --view lidc-bundle-v3`, with no move of the test folder.
- Front 6.1 step 4c: "The V1 bundle holds the frozen `lidc-val-edge` view (B1.6) and `lidc-val-shifted`, and no
  `lidc-test-*` folder ... so the test view is moved out before bundling (or a bundle `--view` filter is added)."
- Code: the bundle ships every `data/eval_views/lidc-*` folder [code: experiments/phase2/vast/bundle.py:53-62].

Why minor: step 4c is a named $0 check before the rental, and it fails and points at the folder. Fix: in A3, replace
the "reaches the box" sentence with "`lidc-test-shifted` is moved out of `data/eval_views` before the V1 bundle (front
6.1 step 4c) and comes back for V3", and add that move line to the header block between `--hash-view` and `bundle`.

### 3. MINOR. `noise_floor.item_sd` is labelled "per-item K = 8 noise", but the formula gives the SD of a two-sampling difference; B1.9's record re-run can double the per-arm variance

Files: A (A3 step 7), B (B1.9 noise model, B.6.2 "Noise floor").

Doc lines:
- A3 step 7: "`noise_floor.item_sd = sd_i(d_i) / sqrt 2` (per-item K = 8 noise)".
- B1.9: "Item i gets independent noise on each arm, with variance `c x mu_i (1 - mu_i) / 8`" and "After B-S1, re-run
  with A3's split-half `item_sd` (A3 step 7), for the record only."

Proof: Var(d_i) = 2 x sigma^2 / 4 = sigma^2 / 2, so `sd(d_i) / sqrt 2` = sigma / 2. One K = 8 mean has SD
sigma / sqrt 8 = sigma / 2.83; the difference of two independent K = 8 means has SD sigma / 2. So `item_sd` is the
two-sampling difference SD (the quantity of the measured 0.0569 pair), and the per-arm K = 8 SD is
`item_sd / sqrt 2`. An engineer who reads "per-item K = 8 noise" and puts `item_sd^2` in as the per-arm variance
doubles the noise of the record run. Record only, so minor.

Fix: A3 step 7: "`noise_floor.item_sd = sd_i(d_i) / sqrt 2`: the per-item SD of the difference of two independent
K = 8 samplings (comparable with the measured 0.0569); the per-arm K = 8 SD is `item_sd / sqrt 2`." B1.9: "c from
A3 = mean_i(item_sd_i^2) / mean_i(2 mu_i (1 - mu_i) / 8)", the same form as the measured 0.077.

### 4. MINOR. The round-11 spec precedence rule landed in A only; `B.spec.md` is still revision 2 and contradicts Gate B1, the val halves and the T-list

Files: B (header, B3.6, B.6.2 Gate B3, B.9.0), `phases/B.spec.md` (copied by DOC1 step 0 to `docs/implementation2.0/specs/`).

Doc lines:
- `B.spec.md` header: "Revision 2 (2026-11-08) ... This file is the contract for test writers."
- `B.spec.md` §2: "A val item is in val-shape when the first byte of SHA-256(group_id) is even, else in val-select"
  (B1.0: stratified `deal` with the first passing seed `"val_half:<s>"`, s = 78 today).
- `B.spec.md` Gate B1: "CI lower bound of (B1 - R5) IoU_tc > 0 (R5 0.599); miss share at offset > 0.3 <= 23.9 % (R5
  31.8 %)", on val-select with "a nodule-level" paired CI (B.6.2: whole val, 0.618, 20.5 %, percentile patient CI,
  DiD and edge rules).
- B.9.0: the test sets "were written from `B.spec.md` and the repo"; B3.6 and Gate B3 still point to `B.spec.md` for
  "T1-T12, T14 and T15", which it does not list.
- A.md header now says "Where `A.spec.md` and this phase section differ, this phase section is the rule". B.md has no
  such line. The r11 verdict kept r11 data-2 minor because "its fix (header precedence rule) rides with data-1"; the
  round-11 A log shows the rule was added to A.md only.

Why minor: U49 (B.md) pins the whole-val Gate B1 and would fail a spec-built evaluator at $0; U02 pins the halves.
Fix: add the A.md precedence sentence to the B.md header (and to C, D, E, DOC, CLEAN, whose specs are also older),
list T1-T12, T14, T15 in B3.6 itself, and either update `B.spec.md` §§2, 5 to B1.0 / B.6 or cut its §5 to a pointer
to B.6 before DOC1 copies it.

### 5. MINOR. Gate B1's far-miss point is an absolute 20.5 % while the baseline is A3's re-sampling; "25 % below R5's 27.3 %" reads as a relative rule

Files: B (B.6.2 Gate B1), front (1.1 row 2 reads it).

Doc lines: B.6.2: "'R5' below means the paired baseline of the common rules: A3's value, with R5's 0.618 for
reference" and "val misses at offset > 0.3: point <= 20.5 % (25 % below R5's 27.3 %)". U49 pins 20.5 %.

Measurement: A3 is a fresh K = 8 sampling of the same model. With the measured noise (c = 0.077 on R5's per-item
means), A3's far-miss share has mean 0.268, SD 0.0046, 90 % range 0.261-0.276; 0.75 x A3 then ranges 0.196-0.207
against the fixed 0.205 (c = 1: up to 0.214) [measured: `(planning probe, not kept: r12data)`, 2,000 draws, 406 far
items, 2026-11-09]. The two readings give different decisions only near the point, and U49 fixes the absolute form,
so minor. Fix: write "point <= 20.5 % (an absolute value: 0.75 x R5's 27.3 %, fixed before V1; it does not move with
A3's far-miss share, which the report prints beside it)".

### 6. MINOR (still open from earlier rounds; not re-argued). Listed so they are not lost

Re-read 2026-11-09; each text is unchanged in the sources:
- r11 data-3: A3 step 5 still says "The power of this test is not measured"; the measured pair bounds it (90 %
  patient half-width about 0.0027 at item SD 0.0569).
- r10 data-2: A8 step 5 still says "B1 then reuses these views and only audits them" (DG12 counts the rows root).
- r10 data-3: front 6.2 "Fallback" still says "must equal the one in the view manifest that R5's eval scored" (no
  such manifest; B.6.2 states the measurable B1.0 hash chain).
- r10 data-4: B1.6 still gives two edge laws ("uniformly on the boundary" and "the largest legal shift in a random
  direction").
- r9 data-4: B1.4 step 2 still says the feasible `draws` "is not yet measured" (r9 pool simulation: 8 without
  disagreement reuse; 2 with reuse <= 4).
- r9 data-5 list: B1.1 test "S0 windows are >= 0.75 of off-centre windows" on the pool (0.508 measured); B1.9 "The
  decision interval of the simulation is the CR1 interval" (B.6.2 uses the percentile interval at >= 50 patients);
  B1.3 `tokens` "(+1 end token)" plus DG12 "max + 1"; A3 step 5 "A lower correlation means a different view: stop"
  (the next bullet and B.6.2 say it drops only the R5 line).
- r5 data-5: the far-miss point sits at the 25 % design effect, so B1.9's joint P(pass) at 25 % is near 0.5 by
  construction (B1.9 says this for Gate B3, not for Gate B1).
- r4 data-8: LIDC-IDRI-0085 and -0146 (17 val items) are on the third-party duplicate-z list but are not B1.7 step 1
  (a)-(d) hits, so step 2 never checks their headers, while 99_back 8c says "[keep, for test and val]".
- r7 data-5: ADR-0015's "near-subset headroom" has no A2 output key.

## Counts

Blocking 0. Major 1 (item 1). Minor 5 (items 2-6).
