# Final review round 10: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).

Scope. I read the doc from its sources by section:
- A: header "Words", A.1, A.3, A1 (steps 1-6 and 5b), A2, A3 (plans, steps 0-8), A4, A5, A8, A9, A.6, A.7, A.8.
- B: B1.0-B1.9, B.6.1, B.6.2, B.7.
- Front: 0, 1 (including 1.1), 6, 6.1, 6.2.
- DOC: DOC6, DOC11.
- 99_back: "Data and quality".

I checked the round-9 fixes in `review-log.md`: data-1 (DOC6 reads both schemas; A1 step 5b `--resummarize`; A9 headline line) and fresh-arch-1 (subshell in `a3_measure.yaml`). They are in place in A and DOC, so I do not repeat them.

My scripts are in `(planning probe, not kept: r10data)` (`dupz.py`, `counts.py`). Both ran with `.venv/bin/python -I` on 2026-11-09. Each script takes the repo path as an argument.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Nothing in the data scope needs an NVIDIA product choice. B1.7's DICOM helper stays "recommend only". I ran no catalog lookup and installed nothing. |
| engineering:code-review | I read the code behind each data claim: `final_eval.py` (`summarize`, `main`), `dataset._prepare/records/prompts`, `build_eval_views._record`, `augment.{row_rng,shift,draw,neighbour_side}`, `lidc.{parse,slice_index,sop_z,ct_slices,clusters,slices}` and `bundle.members`. |
| engineering:debug | I traced every Gate B1 input through the V1 failure paths (front 6.2 "Fallback", 6.3 item 1b, B.6.2 "Paired baseline", U49). That trace found item 1. I re-ran the B1.7 check (a) on the local XML. |
| engineering:system-design | I traced the data flow of each Gate B1 rule: frozen view, then V1 item 3/4/5 outputs, then the evaluator. The edge rule's baseline exists in only one place, item 4 (item 1). |
| engineering:architecture | I checked which section owns each rule. Front 6.2 owns the fallback, B.6.2 owns the gate rules, and U49 tests them. All three are silent on the edge rule under the fallback. |
| engineering:documentation | Each item quotes the doc line, then the code line or measurement that proves it. |
| huggingface-skills:huggingface-best | No model choice is in scope. Front section 3 stays as it is. |
| huggingface-skills:huggingface-datasets | No Hub dataset is read in this scope. The rows-root Hub path was verified in round 9, and nothing in it has changed. |

Web checks (2026-11-09):
- Cameron and Miller 2015, J. Human Resources 50(2):317-372, DOI 10.3368/jhr.50.2.317: the citation is confirmed [web: https://econpapers.repec.org/RePEc:uwp:jhriss:v:50:y:2015:i:2:p:317-372].
- The CR1 factor for a mean, `G/(G-1) x (N-1)/(N-K)` with K = 1, reduces to `P/(P-1)`. So the formula in A1 step 5 and B.6.2 is exact for the item-weighted mean.

## Claims verified (no item)

- **Counts.** Train has 7,341 items, all `rl`, and 289 are disagreement items. The val shifted view has 1,084 items, all `rl`, and 44 are disagreement items. The 1,084 R5 val item ids all occur in the local `lidc-val-shifted` view [measured: `r10data/counts.py`].
- **R5 runs.** `runs/r5/eval/` holds only `eval_val` and `eval_test`. No R5 edge-view file exists [measured: `os.listdir`]. B1.9 says the same.
- **B1.7 check (a).** On the local XML it flags 3 store series: LIDC-IDRI-0013 (train, 1 shared z), -0017 (train, 5) and -0659 (test, 11). It flags no val series, and the pylidc `zvals` hold no duplicate z for any store series [measured: `r10data/dupz.py`]. This matches round 4 data-8, which is still open (item 5). LIDC-IDRI-0085 and -0146 are not hits under steps 1(a) and 1(d).
- **Patient key.** `group_id = scan.patient` [code: src/objdel/data/lidc.py:316]. So the two-scan patients stay in one cluster for the patient bootstrap and in one split.
- **A2 legacy check.** It runs on view caches at frame resolution. Each view crop holds the whole doctor union [code: augment.shift keeps the union inside with CLEARANCE]. So the legacy IoUs on the view equal those on the store masks, and the ±0.005 check against 0.7232 / 0.7319 is sound.

## Items

### 1. MAJOR. In the V1 fallback, the Gate B1 edge rule has no baseline. No section says what the evaluator does then

Files: 00_front (6.2 "Fallback", 6.3 item 1b, 1.1), B (B.6.2 "Paired baseline" and the Gate B1 edge rule, B.9.1 U49).

Doc lines:
- 00_front, 6.2 notes: "Fallback, only if the R5 process (items 4-5) gives no output: Gate B1 pairs with R5's
  `items.jsonl`. The view-identity rule is then "hash only, correlation not run" ... B.6.2 and U49 follow this rule."
- 00_front, 6.2 table, item 4: "R5 (`rl-full-3ep`) on `lidc-val-edge`: the Gate B1 edge baseline (process start +
  first view)". Items 4 and 5 run in one process (`a3_measure.yaml`).
- 00_front, 6.3 item 1b: "Fails again: Gate B1 keeps the 6.2 fallback". Section 1, "End state under (b')":
  "Under (b') no rental follows V1."
- B, B.6.2 Gate B1: "guard law: val IoU on the `edge` view, non-padded and non-fallback items: CI lower bound of
  (B1 - R5) > -0.02 (R5 on the edge view is front 6.2 item 4)". Also: ""Pass" needs every point rule and every CI rule."
- B, B.6.2 "Paired baseline", fallback bullet: it replaces only the shifted baseline ("R5's `eval_val/items.jsonl` is
  the baseline"). It says "Gate B1 waits for A3-val in V2 item 1b" only when the hash chain is broken.
- B, U49: "**item 5 gave no output** (front 6.2 fallback): the evaluator uses R5's items with the hash-chain identity
  ...". The test has no case for item 4 missing.
- 00_front, 1.1 row 4: "Fails ... the edge guard ... | Stop. Owner review".

Proof:
- No older R5 edge values exist. `runs/r5/eval/` holds `eval_test`, `eval_test.log`, `eval_val` and `eval_val.log` only
  [measured: `r10data/counts.py`, 2026-11-09]. B1.9 says the same: "no R5 edge-view values exist before V1:
  `runs/r5/eval/` holds `eval_val` and `eval_test` only".
- The only R5 edge source is item 4. Item 4 runs in the same `final_eval` process as item 5, edge first [code: A3
  plan `a3_measure.yaml`, one `-m experiments.phase2.final_eval` step with two `--view`]. So the fallback state "the R5
  process gave no output" always lacks the edge baseline. Under (b'), V2 item 1b never runs. Under a raise, item 1b
  can fail again. Either way the edge rule cannot be computed.
- What an engineer would do. The fallback text says "Gate B1 pairs with R5's `items.jsonl`" as if the gate can be
  decided, but the edge rule has nothing to pair with. The engineer must guess one of three readings:
  - skip the rule: Gate B1 can then pass with no guard-law check, against B.6.2's "every ... CI rule";
  - count it as a failure: decision table row 4, "Stop. Owner review";
  - mark it open: no row applies.
  These readings pick different rows of the decision table in 1.1. That table is approved before V1 "so the result
  cannot bend it". This is a missing rule that changes the outcome, so it is major and not blocking: it applies only
  on a failure path, and no money is spent for nothing.

Fix:
1. 00_front 6.2 "Fallback": add this sentence. "If item 4 gave no output, the edge rule has no baseline. Gate B1 is
   then "open: edge baseline missing" and is never a pass. Under a raise, V2 item 1b gives the baseline. Under (b'),
   or when item 1b fails again, the gate takes row 4 of section 1.1 (stop, owner review), with the reason "edge
   baseline missing". The other Gate B1 rules are computed and printed."
2. B.6.2, "Paired baseline": add the same bullet, and make the edge rule's text say "(R5 on the edge view is front 6.2
   item 4; without it, see the fallback bullet)".
3. U49: add two cases. (a) Item 4 and item 5 both missing: the result is "open: edge baseline missing", never a pass,
   and the decision row is 4 when the planted flag `raise=False` is set. (b) Item 4 present and item 5 missing: the edge
   rule pairs with item 4, and the shifted rules use the R5 hash-chain fallback.
4. 1.1 row 4: add "or the edge baseline is missing after the fallback (6.2)".

### 2. MINOR. A8 and A.8 say B1 "reuses" the A8 train views and "only audits them". B1.4 builds a new rows root, and DG12 counts on that root

Files: A (A8 step 5, A.8 "Dependencies").

Doc lines:
- A8 step 5: "B1 then reuses these views and only audits them."
- A.8: "B1 uses ... the A8 length rows of the shifted and rescaled views (same manifest SHAs); it audits those views."

Proof:
- B1.4 step 3 materializes new windows for `draws` up to 8 into `data/views/lidc-train-b1`.
- B1.1 drops the `kept_side` fallback windows that the `lidc-train-rescaled` view contains [code:
  build_eval_views.py:57-62].
- B1.3 `tokens(rows_root, ...)` (DG12) counts on the rows root, not on the A8 views.

So A8's draw-0 counts are context only. No gate result changes, because DG12 reads the right root.

Fix: in A8 step 5, write "B1's DG12 counts the rows root itself (B1.3 `tokens`). The A8 rows of the draw-0 train
views are context for the limit choice." In A.8, drop "it audits those views".

### 3. MINOR. Front 6.2 states the fallback view-identity check in a form that cannot be measured. B.6.2 states the measurable form

Files: 00_front (6.2 "Fallback"), B (B.6.2, B1.0).

Doc lines:
- Front 6.2: "the frozen `lidc-val-shifted` `target_set_sha256` must equal the one in the view manifest that R5's
  eval scored".
- B.6.2: "Its view identity is then the measured hash chain of B1.0 (the R5 launch stamp's data hash equals the
  recomputed one ...), recorded in `runs/a3/val_pre_rebuild.json`".
- B1.0: "R5 scored a val view built on the box, with no hash".

Proof: `runs/r5/eval/eval_val/final_eval.json` holds no view manifest or view hash. Its `view` field is
`/workspace/data/lidc/eval_views/lidc-val-shifted` [measured: 2026-11-09]. So "the view manifest that R5's eval scored"
exists only through B1.0's stamp chain.

Fix: in front 6.2, replace the clause with "the B1.0 hash chain recorded in `runs/a3/val_pre_rebuild.json` holds
(B.6.2)".

### 4. MINOR. B1.6 defines the edge offset law in two ways that do not give the same distribution

Files: B (B1.6, I32).

Doc line: "draws the offset uniformly on the **boundary** of the admissible rectangle (the largest legal shift in a
random direction)" and "try up to 64 directions from `row_rng(1, "edge", item_id)`".

Problem:
- A point drawn uniformly on the rectangle's perimeter is one law.
- A ray cast at a uniform angle to the boundary is another: its density per unit length falls toward the corners.

The frozen hash pins whichever form the engineer picks, so no gate breaks. But the view's offset law is then not
defined by the doc.

Fix: pick one form. The "direction" form fits the 64-tries rule: "draw an angle uniformly in [0, 2 pi); move from the
centred corner along it to the admissible rectangle's boundary; round to the nearest admissible integer offset on the
boundary". In I32, state that the angles, not the boundary positions, are uniform.

### 5. MINOR (still open from earlier rounds; not re-argued). Listed so they are not lost

- r4 data-8: B1.7 and DG17 still name LIDC-IDRI-0085 and -0146 (17 val items). The doc does not say whether they go to
  `known_label_issues`. My re-run today again shows that check (a) flags no val series [measured: `r10data/dupz.py`].
- r9 data-2, data-3, data-4 and the r9 data-5 list: unchanged in B (BudgetStop inside a block; the 50-patient
  interval switch on the near subset; "draws not yet measured"; the pool S0-share test; U09 `dropped` keys; the B1.9
  CR1 wording; the DG15 val-shifted branch; the DG12 "+1").
