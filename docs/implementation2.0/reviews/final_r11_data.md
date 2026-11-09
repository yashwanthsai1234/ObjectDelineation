# Final review round 11: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).

Scope. I read the sources by section, not the built file:
- A: header ("Words", revisions), A.4 rental lines, A1 (steps 1-6, 5b), A2, A3 (files, plans, steps 0-8), A4, A5,
  A8, A9, A.6, A.7, A.8, A.9 header, A.9.0, U08-U12, U51, I07, I22, G01, A.9.8.
- B: header, B.1, B.3.2, B.3.4 (ADR-0015), B.4 task order, B1.0-B1.9, B.6.1, B.6.2, U07, U09, U49, U60, U63, I07, I32.
- Front: 0, 1.1, 6, 6.1, 6.2 (incl. the round-10 "edge rule under the fallback"), 6.3 item 1b.
- DOC6 (headline schema handling), CLEAN10, 99_back "Data and quality".
- The test-writer contracts that the phases name as rules: `phases/A.spec.md`, `phases/B.spec.md` (DOC1 copies
  them unchanged to `docs/implementation2.0/specs/`; `diff` of the B pair: identical).

I checked `review-log.md` round 10: data-1 (edge baseline under the fallback) is in front 6.2, 1.1 row 4, V2 item 1b,
B.6.2 "Edge baseline" and U49. I do not repeat it.

Scripts (all `.venv/bin/python -I`, repo path as argument, 2026-11-09): `(planning probe, not kept: r11data)`, `offsets.py`,
`rebuild_val.py` (output `rebuild_val.out`), plus inline `-c` checks quoted below.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in the data scope (B1.7's DICOM helper stays "recommend only"). No catalog lookup, nothing installed. |
| engineering:code-review | Read the code behind each claim: `final_eval.py` (`summarize`, `main`), `signal.{cluster,draws,_ci,ci}`, `build_eval_views._record/build`, `augment.{row_rng,shift,draw,neighbour_side}`, `dataset.{_prepare,records,prompts,grpo_rows}`, `cache.rl_eligible`, `store.load`, `tests/test_gpu_scoring.py` (exact CPU/CUDA equality, so A-G3's "to the bit" rescore is sound). |
| engineering:debug | Reproduced every load-bearing number before relying on it (patient CI, CR1, Gate B1 baselines, noise pair, golden size), and rebuilt `lidc-val-shifted` with today's code to test B1.0's equal-hash branch. |
| engineering:system-design | Traced the summary-key flow A1 -> A4 `paired` -> A-G4 / U51 `gate_a` and found that the key list A1 points to (A.spec 3.2) gives a different CI unit (item 1). |
| engineering:architecture | Checked which file owns each rule: A.md / B.md own the gates, but A1, B3.6 and Gate B3 point to the spec files for keys and T-lists, and the spec files still carry revision-2 gates (items 1, 2). |
| engineering:documentation | Each item quotes the source line and the code line or measurement that proves it. |
| huggingface-skills:huggingface-best | No model choice in scope. Front section 3 unchanged. |
| huggingface-skills:huggingface-datasets | No Hub dataset is read in the data path. I used the Hub only to verify the third-party duplicate-z card (web check below). |

Web checks (2026-11-09):
- HF card `YongchengYAO/LIDC-IDRI-Lite`: names LIDC-IDRI-0085, -0146, -0418, -0572, -0979, "two or more DICOM slices
  sharing one `ImagePositionPatient` z", excluded (1,018 -> 1,013 CTs) [web:
  https://huggingface.co/datasets/YongchengYAO/LIDC-IDRI-Lite]. B1.7's list is quoted correctly.

## Claims verified (no item)

- **Patient CI and CR1.** Val: 1,084 items, 78 patients, 200 nodules, 0 items without `group_id` in the store; patient
  percentile CI [0.5851, 0.6526] around 0.6182; CR1 [0.5837, 0.6527]; patient-weighted mean 0.6472; 2-78 items per
  patient. Test: 2,200 items, 158 patients, 401 nodules; patient CI [0.6252288200704239, 0.6644106494958667] exactly
  (A1 step 6, A-G1) [measured: `r11data/ci.py`].
- **Golden size.** `r5_test_items.jsonl.gz` = 33,761 bytes with `json.dumps` default separators, level 9, mtime 0
  (A1 Files) [measured: inline check].
- **Gate B1 baselines.** Val offset median 0.257; offset > 0.3: 406 items, 75 patients, miss share 27.3 %; offset <
  0.15: 206 items, 52 patients, IoU 0.751; found IoU 0.706; R5 val legal 0.9952 [measured: `r11data/offsets.py` with
  the union-box / window-box definition of B.1 item 2; `final_eval.json`].
- **Noise pair (B.6.2, B1.9).** 128 items: per-item difference SD 0.0569, variance 0.00324, Pearson 0.971, mean
  -0.0076, c = 0.0768 [measured: inline check of `runs/rl/eval/final/heldout.json` against R5 val items].
- **B1.0 equal-hash branch holds today.** Rebuilding `lidc-val-shifted` with today's code (`build_eval_views`, shifted
  kind only, seed 0) gives the same B1.6 content hash (`405b3b97...`), the same `target_set_sha256` (`2da645bf...`),
  0 differing `image_sha256`, boxes, targets, raters, `disagreement`, `rl` and offsets on all 1,084 items [measured:
  `r11data/rebuild_val.py`, 2026-11-09]. So B1.0 step 2's difference branch is not expected; the fingerprint split and
  `view.offset_frac` do not enter the content hash.
- **A-G3 rescore to the bit.** `tests/test_gpu_scoring.py:124,139,165` assert `device_terms(m, cache) ==
  host.score_mask(m, cache)` (exact), so the Mac CPU rescore of box answers can match bit for bit.

## Items

### 1. MAJOR. A1 delegates the output keys to `A.spec.md` §3.2, which is still revision 2: it gives Gate A's paired difference a nodule CI, while A4 and A-G4 decide on the patient CI. No test pins the unit

Files: A (A1 "Config keys" line; A4 selection block; A.6 A-G4; A.9 U11, U51, I22). Also `phases/A.spec.md` (copied by
DOC1 to `docs/implementation2.0/specs/A.spec.md`).

Doc lines:
- A.md, A1, last line: "Config keys: none (CLI only). Output keys: listed in `docs/implementation2.0/specs/A.spec.md`
  §3.2, plus `gate_ci`."
- A.md header: "The test writers' contract is `docs/implementation2.0/specs/A.spec.md`."
- A.md, A4: "Paired differences `sel_medoid - sel_mean_of_k` and `sel_vote - sel_mean_of_k`: `signal.ci` of the
  per-item differences under patient draws (the gate CI), with the nodule CI as a second line."
- A.md, A-G4: "Val: `sel_medoid - sel_mean_of_k` (and `sel_vote - sel_mean_of_k`) reported with paired patient CIs
  (78 patients). The default changes only if the val mean difference is >= 0.015 and its patient CI lower end is > 0."
- A.md, A5 step 4: "`iou_found` (the mean IoU_tc over found images, with the patient CI on that subset; nodule CI as
  a second line)".
- A.spec.md header: "Revision 2 (2026-11-08): matches `A.md` revision 2." (A.md is at revision 5 plus rounds 4-10.)
- A.spec.md §3.2: "`paired`: `{"medoid_minus_mean": [d, lo, hi], "vote_minus_mean": [d, lo, hi] or null}`. These are
  the mean per-item differences with nodule-cluster CIs." and "`iou_found`: the mean IoU_tc over found items, with a
  nodule CI on that subset."
- A.spec.md §5 item 8 (Gate A): "the medoid (or vote) becomes the default only if the val mean ... is >= 0.015 and
  its nodule CI lower end is > 0." Item 7 (A-G3): "the nodule CI of the per-item IoU_tc difference contains 0 and the
  mean difference is within ±0.01" (A.md A-G3: TOST, 90 % patient CI inside ±0.01, plus view identity and the
  frozen-view hash). §3.2's `final_eval.json` key list also has no `view_content_sha256`, which A3 writes and Gate
  B1's frozen-view rule reads.

Proof that no planned $0 check catches the wrong unit:
- U11 checks only `medoid_minus_mean[0] == 0.0125` and `lo <= d <= hi` (no cluster unit).
- U51 `final_eval.gate_a(val_summary, ...)` "reads the patient CI (`gate_ci`)" but takes `medoid_minus_mean` as one
  `[d, lo, hi]` triple. Fed a nodule triple, its hand cases pass unchanged.
- I22 checks only that `paired.medoid_minus_mean` reaches `final_eval.json`.
- The unit changes the decision. On R5 val-select the patient CI is 23 % wider than the nodule CI (0.085 vs 0.069)
  [measured: B.6.2, `docs/implementation2.0/notes/b_final/cluster_ci.json`]. Today `signal.draws` has only the nodule
  key [code: src/objdel/eval/signal.py:102-107], so the nodule CI is the default path. A medoid gain whose nodule
  lower end is just above 0 and whose patient lower end is below 0 changes the report default under the spec and
  not under A-G4.
- Two parts both claim to be the rule (the A1 pointer plus the "test writers' contract", against A4 and A-G4). The
  engineer must guess. Gate A is decided on the Mac after V1, so no money is lost, but the answer-choice default for
  every later report can differ. Major, not blocking.

Fix:
1. A1: replace the pointer with the key list itself: "Output keys (schema 2): top level `schema, items, nodules,
   patients, ci_level, boot, gate_ci`; per METRICS and selection key `{mean, ci95_nodule, ci95_patient, ci95_item,
   higher_is_better}`; `iou_found` and `miss_rate` with the same CI keys (subset draws); `paired.{medoid_minus_mean,
   vote_minus_mean}` = `{mean, ci95_patient, ci95_nodule}` (patient draws are the gate CI); run keys `ckpt, view,
   split, K, seed, seconds, view_manifest_sha256, view_content_sha256, code_sha, checker_build, partial`."
2. U11: recompute the paired CI with `signal.draws(ids, key=groups.__getitem__)` and assert `ci95_patient` equals
   it (1e-12), with `ci95_nodule` as the second key. U10: the same for `iou_found`. U51: read
   `paired.medoid_minus_mean[summary["gate_ci"]]`, and add a hand case where the nodule lower end is > 0 and the patient
   lower end is < 0 (expect no change).
3. `phases/A.spec.md`: bring §3.2 (paired, `iou_found`, `view_content_sha256`, `gate_ci`, `cluster_t`) and §5 items
   7-8 (TOST 90 % patient CI, view identity, frozen-view hash; Gate A on the patient CI) in line with A.md, and set
   its revision line to the A.md revision. Add one rule to the header of every `*.spec.md`: "Where this spec and the
   phase section differ, the phase section is the rule." DOC1 copies the corrected files.

### 2. MINOR. `B.spec.md` (revision 2) still states the old data gates and Gate B1, and B3.6 / Gate B3 point to it for a T-list it does not hold

Files: B (header, B3.6, B.6.2 Gate B3). Also `phases/B.spec.md`.

Doc lines and proof:
- B.md header: "`docs/implementation2.0/specs/` (`B.spec.md`)" as the contract; B3.6: "List in
  `docs/implementation2.0/specs/B.spec.md` (T1-T12, T14; T13 is removed with `adv_scale`; T15 is new in revision 4)";
  Gate B3: "T1-T12, T14 and T15 pass (`B.spec.md`)". B.spec.md §5.2 lists T1-T13 and has no T14 or T15.
- B.spec.md §5.1: D7-D9 measure against "val-shape target cell shares" and "val-shape offsets" (B1.4 step 1 and
  ADR-0015: the train shifted law; "No val data shapes training", B.3.2 item 2); D3 "rate <= within-train rate"
  (DG3: 0 unreviewed pairs); D13 "on 150 sampled windows" (DG13: every window); D15 test view only (DG15: three
  frozen views). §5.3 Gate B1 is on val-select with "R5 0.599", "<= 23.9 % (R5 31.8 %)", "R5 0.732", no
  difference-in-differences rule and no edge rule (B.6.2: whole val, 0.618, 20.5 %, 0.751, DiD and edge).
- Why minor: each difference fails a named B test on the Mac at $0 and points at it: I07 (target bit-identical when
  any val item changes), U49 (Gate B1 rules on the whole val), the B1.3 DG3 test ("fails DG3 until the pair is marked
  reviewed"), DG13's every-window test, and U63 / the T14 sketch.

Fix: B3.6 and Gate B3: list T1-T12, T14, T15 in B.md itself and drop "List in B.spec.md". Update `phases/B.spec.md`
§5.1-5.3 to B.6.1 / B.6.2 (train shifted law, every window, three frozen views, whole-val Gate B1 with DiD and edge),
or cut §5 to a pointer to B.6, under the same header rule as item 1.

### 3. MINOR. A3 step 5 says the TOST power "is not measured"; the measured noise pair already bounds it

Files: A (A3 step 5, A-G3).

Doc line: "TOST: the 90 % patient CI (paired, 78 patients) of the mean per-item IoU_tc difference must lie inside
±0.01 ... The power of this test is not measured [estimate: the paired CI is narrower than the unpaired one, by how
much is unknown]."

Measurement: with independent per-item differences of SD 0.0569 (the measured two-box pair), the 90 % patient-bootstrap
half-width on the 1,084 val items of 78 patients is 0.0027 (max 0.0032 over 50 simulations; B 2,000) [measured: inline
simulation on the R5 val ids and store `group_id`, 2026-11-09]. So the TOST passes only if the systematic A3 - R5
offset is below about 0.0073. The one measured two-box pair had a mean difference of -0.0076 (128 items, SE about
0.005). A failure does not stop Gate B1 (B.6.2), so this changes no spend.

Fix: replace the bracket with "[estimate: 90 % patient half-width about 0.0027 under the measured per-item SD 0.0569;
the TOST passes only if the box-to-box offset is below about 0.007; the one measured two-box pair gave -0.0076 on 128
items]". Add this line to the A-G3 report.

### 4. MINOR (still open from earlier rounds; not re-argued). Listed so they are not lost

- r10 data-2: A8 step 5 still says "B1 then reuses these views and only audits them", and A.8 "it audits those views".
  DG12 counts the rows root (B1.3 `tokens`).
- r10 data-3: front 6.2 "Fallback" still says "the frozen `lidc-val-shifted` `target_set_sha256` must equal the one in
  the view manifest that R5's eval scored" (no such manifest exists; B.6.2 states the measurable B1.0 hash chain).
- r10 data-4: B1.6 still defines the edge law as "uniformly on the boundary" and as "the largest legal shift in a random
  direction" (two different laws).
- r4 data-8: LIDC-IDRI-0085 and -0146 (17 val items) are on the third-party duplicate-z list (web check above) but are
  no hit under B1.7 step 1 (a)-(d), so step 2's DICOM header check never runs for them, while 99_back 8c says "[keep,
  for test and val]". Fix: add step 1 (e) "series on the third-party list" so step 2 checks their headers, and state
  whether they enter `known_label_issues`.
- r9 data-2 (BudgetStop inside a block), r9 data-3 (near rule at 52 patients, 51 without 0085/0146), r9 data-4 (the
  `draws` capacity numbers can replace "not yet measured").
- r9 data-5 list: B1.1 test "S0 windows are >= 0.75 of off-centre windows" on the pool (0.508 measured); U09 "exactly
  the 3 int keys" (B1.1 defines 4, `padded` included); B1.9 "The decision interval of the simulation is the CR1
  interval" (B.6.2 and U49 use the percentile interval at >= 50 patients; also, CR1 is defined for one item-weighted
  mean, not for rule (4)'s far-minus-near difference); DG15 "equals the B1.0 pre-rebuild hash" on B1.0's difference
  branch (today's rebuild gives the equal hash, see "Claims verified", so this branch is unlikely); B1.3 `tokens`
  counts "+1 end token" and DG12 adds "+1" again; A3 step 5 "A lower correlation means a different view: stop".
- r7 data-5: ADR-0015's "near-subset headroom" has no A2 output key (A2 strata are `all`, `k>=3`, `k=2`, `found`,
  `miss`).

## Counts

Blocking 0. Major 1 (item 1). Minor 3 (items 2-4).
