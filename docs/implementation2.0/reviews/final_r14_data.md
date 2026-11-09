# Final review round 14: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice is in the data scope. B1.7's DICOM helper stays "recommend only". No catalog lookup and no install. |
| engineering:code-review | I read the code behind each data claim I re-checked: `dataset._prepare` (per-window `disagreement`), `dataset.prompts`, `cache.rl_eligible`, `ExpertCache.__init__` (vote, voters), `score.py` (IoU_tc falls back to the vote), `build_eval_views._record`, `augment.{draw,shift,neighbour_side}`, `preflight.hashes`, `bundle.members`, `final_eval.main` defaults, the `b1_pixels.py` near-duplicate seed. |
| engineering:debug | I reproduced each number before using it (scripts listed below). Two suspicions died this way: (1) per-window disagreement could differ from the source item's (it does not: 0 of 7,341 in both train views); (2) the near-duplicate crop pads on the wrong side at a top or left edge (true in the seed, but no store item is within 16 px of an edge). |
| engineering:system-design | I traced store -> pool -> schedule -> rows root -> DG1b, and val view -> A3 / B1 evals -> Gate B1 / B3 intervals, and asked which planned check still covers each step. |
| engineering:architecture | I checked one owner per rule: B1.3 (DG3 resolution) against B1.7 (exclusion list), B.6.2 (interval choice) against A1 step 5. |
| engineering:documentation | Each item quotes the source line and the measurement. |
| huggingface-skills:huggingface-best | No model choice in the data scope. Front section 3 unchanged. |
| huggingface-skills:huggingface-datasets | Re-checked the rows-root upload claim against the installed `huggingface_hub` 2.2.0: `upload_folder` commits a large folder in batches with " (part N)" messages when `hf_xet` is installed [code: .venv huggingface_hub/hf_api.py:5778-5783], so B1.4 step 5 is right to pin the returned `CommitInfo.oid`. No Dataset Viewer endpoint is on the data path. |

Web tools loaded (WebSearch, WebFetch). Web check, 2026-11-09: MacKinnon and Webb, "Wild bootstrap inference for
wildly different cluster sizes", J. Applied Econometrics 32(2):233-254, 2017 [web:
https://econpapers.repec.org/RePEc:wly:japmet:v:32:y:2017:i:2:p:233-254; working paper
https://ageconsearch.umn.edu/record/274639/files/qed_wp_1314.pdf]: cluster-robust t tests over-reject when cluster
sizes are very unequal, also at about 50 clusters; the wild cluster bootstrap mostly fixes it (item 2).

## Scope and state

- The A and B sources have not changed since round 13 (`phases/A.md` and `phases/B.md` last written 02:40 and
  02:43; the round-13 data review was written 02:54). Round 13 applied fixes only to D and CLEAN
  (`review-log.md` "final round 13" sections). So I did not repeat round 13's A/B reading line by line. I looked for
  what earlier rounds did not measure, and I read the round-13 D and CLEAN changes for data effects.
- Read: A header, A.1, A1, A2, A3, A4, A5, A8, A.6-A.8; B header, B.1, B.3.4 ADR-0015, B.4, B1.0-B1.9, B.6.1,
  B.6.2, B.8.3; front 1.1, 6.1, 6.2; 99_back "Data and quality"; CLEAN5 (round-13 RIGA `.gitignore` fix) and CLEAN10;
  E grep for the B1 data modules; DOC11 datasheet list.
- Scripts (all `.venv/bin/python -I` with the repo path as argument, 2026-11-09), in
  `docs/implementation2.0/`: `r14_dis.py` (window vs item disagreement), `r14_edge.py` (nodule centres near a slice edge),
  `r14_cr1.py` (coverage of the gate intervals with the real val cluster sizes).

## Claims verified (no item)

- Doctors per item: val 281 / 330 / 473 items with k = 2 / 3 / 4; train 1,933 / 2,029 / 3,379 [measured: view and
  prepared records]. Matches A2 and A-G2 (281 + 803 = 1,084).
- Every train item is RL-eligible: prepared `rl` true for 7,341 of 7,341, and the val view manifest has
  `rl_items` 1,084 of 1,084. So B1.1's `dropped["not_rl"]` is 0 today, and Gate B1 scores no single-doctor item.
- B1.1 "every window copies the source item's `disagreement`": in `data/views/lidc-train-shifted` and
  `lidc-train-rescaled`, the per-window `disagreement` that `_prepare` recomputes on the window's own masks equals the
  prepared item's for 7,341 of 7,341 records each (289 true in all three) [measured: `r14_dis.py`]. So the B1.2
  quota (which uses the copied flag) and the materialized records (which `window_record` recomputes) agree today.
- Val disagreement items 44 of 1,084 (B1.2) [measured: view manifest].
- A2's k = 2 rule "IoU_tc equals plain IoU": with the stack `[other, other]` the vote is `other`, and IoU_tc falls back
  to the vote when no doctor shares the answer's signature [code: src/objdel/reward/cache.py:226; src/objdel/reward/score.py:7-8].
  So the rule holds for any topology, not only for matched signatures.
- A3 and B1.5 evals both use seed 11: `final_eval --seed` defaults to 11 [code: experiments/phase2/final_eval.py:53],
  and R5's val eval used seed 11 [measured: runs/r5/eval/eval_val/final_eval.json].
- The offset definition behind Gate B1's subset counts (406 far items, 206 near items) is B1.1's `offset_of`
  (union-box centre to window centre over the side) [code: (planning probe, not kept: r11data)].
- B1.0 rebuild of `lidc-val-shifted` with today's code: same content hash `405b3b97...` and same
  `target_set_sha256` `2da645bf...` [measured: (planning probe, not kept: r11data)].
- Near-duplicate crops: no store item has its union-box centre within 16 px of a slice edge (minimum distance:
  train 16, val 26, test 41 px) [measured: `r14_edge.py`]. The seed pads a clipped crop only at the bottom and right
  (`np.pad(crop, ((0, 32-h), (0, 32-w)))` [code: docs/implementation2.0/seeds/b1_scripts/b1_pixels.py:40-41]), which
  would misalign a crop clipped at the top or left, but no item is affected today.
- CLEAN10 keeps every `data/eval_views/lidc-val-*` folder until B1.0 records the val hash; CLEAN5's new RIGA block
  touches no LIDC data path. The bundle never ships `data/views/lidc-train-*` unless named
  [code: experiments/phase2/vast/bundle.py:52-54], so the rows root `lidc-train-b1` leaves the Mac only by the B1.4
  HF upload.

## Items

No blocking and no major item. Each item below is minor: it either has no effect on a result today, or it changes a
reported number by a small, measured amount, and no proof shows a wrong gate decision.

### 1. MINOR. DG3: the 0.98 threshold sits inside the look-alike tail, and a "drop train item" removes one slice, not the nodule

Files: B (B1.3 `near_duplicates`; B.6.1 DG3; B.1 item 4).

Doc lines:
- B1.3: "pairs (train item, val or test item) and (val item, test item) with score >= 0.98 ... A 'drop train item'
  pair counts as reviewed only when its train item id is also a key of `docs/reference/train_exclusions.json`
  (reason `"DG3 <pair>"`)".
- B.1 item 4: "Two near-duplicate pairs (cos >= 0.98) cross splits".

Measurement [docs/implementation2.0/notes/b1_out/neardup.json, the seed's output]: the best-match score of each test item
against all train items has median 0.847 and p99 0.962; 97 test items score >= 0.95, 1 scores >= 0.98, 0 score >=
0.995. Val -> train: p99 0.964, 37 >= 0.95, 0 >= 0.98, 0 >= 0.995. Train against other train patients: 231 >= 0.95,
4 >= 0.98. So 0.98 is the far tail of a smooth look-alike distribution: a 32x32 mean-subtracted HU crop centred on a
round nodule correlates highly with any other round nodule of similar size. A true duplicate scan (the same CT under
two patient ids) would give scores near 1.0 on every slice of the nodule; today no pair reaches 0.995.

Two consequences:
- The owner reviews look-alike pairs, which is cheap; but the review rule does not ask for the evidence that tells a
  duplicate from a look-alike.
- If a pair is judged a real duplicate, "drop train item" puts only that one slice into the exclusion list. The other
  slices of the same train nodule (about 5.9 per nodule: 7,341 items / 1,247 nodules) stay in the rows root, and
  DG1b passes, because it checks only listed ids.

Why minor: no pair reaches 0.995 today, and the owner reviews each pair with a PNG.

Fix (B1.3, DG3 row): (a) report, per listed pair, the best score of every other slice of the same two nodules and the
count of pairs >= 0.995; (b) "A pair the owner judges a duplicate drops every train item of its nodule
(`item_id.rsplit("/", 1)[0]`), each with reason `"DG3 <pair>"`; a look-alike pair is resolved by the owner's OK";
(c) B.1 item 4: "score (mean-subtracted correlation) >= 0.98", not "cos" (the seed subtracts the crop mean
[code: b1_pixels.py:67]).

### 2. MINOR. With the real val cluster sizes, the gate intervals cover the miss share about 91.5 % of the time, not 95 %

Files: A (A1 step 5), B (B.6.2 "Pairing and CI"; B1.9).

Doc lines:
- A1 step 5: "The percentile interval is the decision value on the whole val split (78 patients) ... A gate on fewer
  than 50 patients (for example a val half) uses a cluster t interval with n_patients - 1 degrees of freedom".
- B.6.2: "the percentile patient interval when the gate set has >= 50 patients (Gate B1, whole val, 78 patients), and
  `signal.cluster_t` when it has < 50 (Gate B3, val-select, about 39 patients)".
- B1.9 test: "with no effect, P(pass) of the IoU rule is <= 0.05 + 0.03" (a synthetic set, not the real clusters).

Measurement [`(planning probe, not kept: r14_cr1.py)`, 2026-11-09]: population = the 78 val patients with R5's per-item values (2 to 78
items per patient). Each repeat draws G patients with replacement and builds the interval of the item-weighted mean;
coverage = share of repeats whose interval holds the population value.

| Interval | G | IoU_tc coverage | miss share coverage |
|---|---|---|---|
| CR1 t (`cluster_t`) | 38 | 0.934 (4,000 repeats) | 0.916 |
| CR1 t | 78 | 0.944 | 0.930 |
| percentile patient bootstrap, B = 1,000 | 78 | 0.938 (1,000 repeats) | 0.915 |

The literature gives the cause: cluster-robust tests over-reject when cluster sizes differ widely, also near 50
clusters [web: MacKinnon and Webb 2017, J. Applied Econometrics 32(2):233-254]. Here the top 5 patients hold 24.9 %
of the val items (A1 step 5).

Why minor: these are one-sample means. The gates use paired differences, and every Gate B1 rule must pass together,
so the joint false-pass rate is lower than one rule's rate. The miss rules (Gate B1 far-miss, Gate B3 miss) are the
ones closest to the edge: a one-sided error of about 4 % instead of 2.5 %. No gate result exists yet.

Fix (no change to any rule): B1.9 adds a null line. With effect 0 and c = 0.077, on the real val clusters (and on
the val-select patients for Gate B3), print the false-pass rate of each rule that reads item values. If a rule's rate
is above 0.05, the report says so beside the gate. Optional second line in A1 step 5: a wild cluster bootstrap
(Rademacher or Webb weights) interval, printed and never deciding.

### 3. MINOR (still open from round 13; not re-argued). The round-13 data items were not applied

Round 13 changed only D and CLEAN. Each text below is unchanged in `phases/B.md` and `phases/A.md` (re-read
2026-11-09):
- r13 data-1: B1.7 still says "steps (a) and (b) have not run yet, so the count can grow"; B.6.1 DG17 "Today", the
  B.7 val row and 99_back 8c still name 0085 and 0146 as val hits. Round 13 measured that (a) and (b) add no series
  and that 0085 and 0146 are not hits.
- r13 data-2: the "mapped" resolution cannot clear the four (d) train hits.
- r13 data-3: "5 of 999 sidecars" (998 sidecars plus `selection.json`).
- r13 data-4: check (c) names "139.xml" (five files have that name; 0396 is not in the store).
- r13 data-5: B1.0 hashes the val view before B1.7 runs, against the A3 header line "no eval view is hashed before it".
- r13 data-6: DG17 "a 'mapped' series that the latest report still lists".
- r13 data-7: DG14b input timing and the "drop" / "weight 0.5" mechanism.
- r13 data-8: B.7 "300 more steps" (260); B.6.2 "about 39 patients" (38).
- r13 data-9 list (r12 data-2, -3, -4, -5; r11 data-3; r10 data-2, -3, -4; r9 data-4, -5; r5 data-5; r7 data-5),
  and the B1.9 CR1 sentence (r5 data-4, r6 data-11).

## Counts

Blocking 0. Major 0. Minor 3 (items 1-3).
