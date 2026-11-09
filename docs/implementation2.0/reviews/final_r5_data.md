# Final review round 5: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).
Scope: the whole assembled doc, read by section, with A1-A3, A8, B1.0-B1.9, B.6.1, B.6.2, B.4.1 and front 1.1 / 6 in
full. Items already fixed in `docs/implementation2.0/review-log.md` (round 4: data-1, data-2, data-3, data-5) are not
repeated.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product decision in the data scope (the B1.7 DICOM helper is already "recommend only"). No catalog install. Nothing changed. |
| engineering:code-review | Read the code each data claim cites: `augment.py` (`draw`, `shift`, `neighbour_side`), `build_eval_views._record`, `splits.py`, `signal.draws`, `hf_pull.py`, `bundle.members`. Found the rows-root path problem (item 2) by reading `hf_pull.py` against the B1.5 text. |
| engineering:debug | Reproduced numbers instead of trusting them: A1 golden CIs, CR1 interval, val strata, far-miss share, a 40-item rebuild of the val views, the disagreement quota share, and a far-miss point-rule simulation (items 5-7). |
| engineering:system-design | Traced the V1 data path end to end (B1.4 build -> HF upload -> `fetch_hf.sh --rows` -> `train_rl --train`) and found where the layout breaks (item 2). |
| engineering:architecture | Checked which text owns the after-V1 decision table: front 1.1, B.4.1 and 99_back 2b all claim it, and they disagree with B.6.2 (item 1). |
| engineering:documentation | Checked each gate for one definition and one owner; listed the stale copies (items 4, 9, 10). |
| huggingface-skills:huggingface-best | No model choice in scope; front section 3 left as is. |
| huggingface-skills:huggingface-datasets | Used for the HF side: read the HF download guide (`snapshot_download` `local_dir` keeps repo paths, item 2) and the `YongchengYAO/LIDC-IDRI-Lite` card (duplicate-z list 0085, 0146, 0418, 0572, 0979: confirmed). |

Web checks: HF download guide (https://huggingface.co/docs/huggingface_hub/en/guides/download: "The downloaded files
will maintain their original file structure within the specified folder"); LIDC-IDRI-Lite card (B1.7 list confirmed).

## Claims verified (no item)

- A1 goldens on `runs/r5/eval/eval_test/items.jsonl` + test store `group_id`: nodule CI [0.6279017988891661,
  0.6620867640136667], patient CI [0.6252288200704239, 0.6644106494958667] (158 patients), item CI [0.6334010210,
  0.6561989222], miss rate 0.10818181818181818. All equal the doc.
- A1 step 5 / B.6.2: CR1 on R5 val = [0.58371, 0.65267] around 0.61819; patient-weighted mean 0.64722; top-5
  patients 24.9 % of items; 2-78 items per patient, median 8. All equal.
- Val: 1,084 items, 78 patients, 200 nodules, k=2 281 / k>=3 803, 44 disagreement items. Store: train 7,341 / val
  1,084 / test 2,200 items, all k >= 2; 557 / 78 / 158 patients. B1.2: 289 disagreement items in train.
- B1.0 rebuild risk: rebuilding the first 40 val items with today's code (`build_eval_views --limit 40`, into
  `(planning probe, not kept: r5data)`) gives 0 differences in `image_sha256`, `view`, `targets`, `box`, `raters`, `caches` for the
  shifted and the rescaled view.
- 99_back 7b: rescaled draw-0 windows under 10 px = 62.0 %; size bins >= 20 px get 38.0 % of the slots; 1,074 items
  reach them. Equal.

## Items

### 1. BLOCKING. Two after-V1 decision tables, and the gate tests, disagree with B.6.2 on view identity and the near guard

Files: 00_front, 99_back, B.

Doc lines:
- 00_front, "1.1 After V1: decision table", row 4: "| Fails the offset < 0.15 guard, the edge guard, legal or view
  identity | Stop. Owner review | B1.8 waits | R5 |". Same text in 99_back "Owner decisions needed", decision 2b:
  "The offset < 0.15 guard, the edge guard, legal or view identity fails: stop, owner review; B1.8 waits; R5." 2b
  says it approves "(front section 1.1, B.4.1)" as if they were one table.
- B, "B.4.1 After V1: decision table", row 4: "| Fail on the edge guard or on legal share | Stop; owner review |
  B1.8 waits | R5 |". No row for the offset < 0.15 guard, no row for view identity.
- B, "B.6.2 Run gates", common rule "View identity (from A-G3)": "If this correlation rule fails, Gate B1 still
  pairs with A3's items (same frozen view, same box); only the R5 report line is dropped." A, "A.6" A-G3: "A val
  failure does not block Gate B1".
- B, "B.9.1" U49 lists "when item 5 ran, view-identity correlation >= the `r_sh` lower bound" among "all hard rules",
  and in the same row "item 5 ran but the correlation rule fails: the gate still pairs with A3's items and drops the
  R5 line". The sketch in "B.9.8" (`tests/test_phase_b_decisions.py::test_gate_b1_needs_all_hard_rules`) asserts
  `("view_corr", 0.74)` (below `view_corr_min=0.75`) makes `gate_b1(...).ok` false.

Why blocking: if B1 passes every B.6.2 rule but the A3-R5 correlation is below the `r_sh` bound, front 1.1 / 99_back
say stop, final candidate R5; B.6.2 / A-G3 say Gate B1 passes on A3's items, so B-S2 is funded. If only the near guard
fails, B.4.1 has no row. The owner approves "the" table before V1, and `make test` must pass: an engineer who codes
`gate_b1` per B.6.2 fails the B.9.8 sketch; one who codes the sketch breaks B.6.2. The outcome (fund B3 or stop) and
the final candidate change.

Fix:
1. One table, in front 1.1 only. B.4.1 and 99_back 2b keep a pointer ("front 1.1 owns the table").
2. Row 4 text, in front 1.1: "Fails the offset < 0.15 guard, the edge guard or legal share | Stop. Owner review |
   B1.8 waits | R5". Delete "view identity" there and in 99_back 2b. Add one sentence under the table: "View identity
   is not a Gate B1 rule: a failed correlation drops only the R5 report line (B.6.2)."
3. U49: delete "when item 5 ran, view-identity correlation >= the `r_sh` lower bound" from the hard-rule list.
4. B.9.8 sketch: remove `view_corr`/`view_corr_min` from `good` and the `("view_corr", 0.74)` case; add
   `assert gate_b1(**good, view_corr=0.74, view_corr_min=0.75).ok and gate_b1(...).r5_line_dropped`.

### 2. BLOCKING. The rows-root download puts the root one level too deep; `train_rl --train {data}/rows/lidc-train-b1` finds no rows on the paid box

Files: B (B1.5, B.9.3 I40), 00_front (6.2 launch lines).

Doc lines:
- B, "B1.5 Train on the ordered rows", "Rows-root download": "runs `hf_pull --subdir rows/lidc-train-b1 --out
  /workspace/data/lidc/rows/lidc-train-b1` with `snapshot_download(..., allow_patterns="rows/lidc-train-b1/**",
  revision=...)`".
- B, "B.9.3" I40: "`snapshot_download` gets `allow_patterns="rows/lidc-train-b1/**"` and the pinned revision, the
  output folder is `/workspace/data/lidc/rows/lidc-train-b1`" (Hub mocked).
- B1.5 step 3: `--train {data}/rows/lidc-train-b1`; `{data}` = `/workspace/data/lidc` [code: policy.py:433].

Proof: `hf_pull.py` calls `snapshot_download(repo_id=a.repo, revision=a.revision, local_dir=a.out, token=token)`
[code: experiments/phase2/hf_pull.py:23]. The HF guide: "The downloaded files will maintain their original file
structure within the specified folder. For example, if `filename="data/train.csv"` and
`local_dir="path/to/folder"`, the resulting filepath will be `"path/to/folder/data/train.csv"`" [web:
huggingface.co/docs/huggingface_hub/en/guides/download]. So the files land in
`/workspace/data/lidc/rows/lidc-train-b1/rows/lidc-train-b1/manifest.json`. `train_rl` reads
`{data}/rows/lidc-train-b1/manifest.json` and `lidc/train.jsonl`: missing. The doc never defines what `--subdir`
does, and I40 mocks the Hub, so no Mac test sees the layout. The failure comes at `b1_rows.yaml` step 3, after
launch, download, build, `env_check` and smoke (item 1 of 6.2, about 53 min, $0.70), and the `&&` chain then loses
the B1 run and its evals (items 2-3).

Fix:
1. B1.5: drop `--subdir`. `cmd_fetch_hf --rows REV ROOT` runs `hf_pull --repo svenk029/objdel-sft --revision REV
   --out /workspace/data/lidc --allow rows/ROOT/**`, i.e. `snapshot_download(local_dir="/workspace/data/lidc",
   allow_patterns=f"rows/{ROOT}/**", revision=REV)`, so the repo path `rows/lidc-train-b1/` becomes
   `{data}/rows/lidc-train-b1/`.
2. `cmd_fetch_hf --rows` ends with a box check over the same ssh call: `test -f
   /workspace/data/lidc/rows/<ROOT>/manifest.json && sha256sum` of it equals the plan's `--rows-manifest-sha256`;
   else exit non-zero (the main agent tears down before `run.sh`).
3. I40: assert `local_dir == "/workspace/data/lidc"` and `allow_patterns == "rows/lidc-train-b1/**"`; add I40b: a
   fake `snapshot_download` that writes the repo-relative paths under `local_dir` (as the real one does), then
   `train_rl`'s root resolution finds `manifest.json` at `{data}/rows/lidc-train-b1`.

### 3. MAJOR. B1.9 asks for Gate B3 pass probabilities but defines no B3 arm

Files: B (B1.9, B.7, B.9.4 E10).

Doc lines: B, "B1.9 Gate power check": "Simulated B1 effect: convert 25 % of far misses ..." and "For each
simulated run, apply every Gate B1 and Gate B3 rule exactly as B.6.2 writes it ... Print P(pass) per rule and
jointly"; "If joint P(pass) < 0.5 for the 25 % effect in the c = 0.077 run, the owner decides before the rental
whether to change a rule." B.7: "Small gate sets (78 patients for B1, 38 for B3) | ... | B1.9 power check before the
rental".

Problem: Gate B3 pairs B3 against B1 on val-select (B.6.2 "Gate B3"). B1.9 simulates only a B1 arm against R5. No B3
effect, no B3 noise and no val-select split are defined, so the Gate B3 rows cannot be computed as written, and
"jointly" is undefined (Gate B1 only, or B1 and B3). The joint value drives the owner's "change a rule" decision, so
two engineers print different numbers. Measurement that shows the gap: E10/E10b (B.9.4, B.9.16) test only Gate B1
rules; no test names a B3 arm.

Fix: in B1.9 write "apply every Gate B1 rule"; "jointly" = all Gate B1 hard rules. Add one optional line: "Gate B3
power (report): on val-select, B3 arm = the simulated B1 values with a further 25 % relative cut of the remaining
misses (chosen by patient, IoU from the same cell's found distribution), both arms noisy with the same c; print P(pass)
per Gate B3 run rule; not used for the V1 decision." Update E10 to assert the joint value covers Gate B1 rules only.

### 4. MINOR (still open from round 4, data-6). B1.9 names CR1 as the Gate B1 interval; the gate uses the percentile interval

Files: B (B1.9, B.6.2, B.8.4).

B1.9: "The decision interval of the simulation is the CR1 interval of the item-weighted mean (B.6.2), the same as
the gate's." B.6.2: "the percentile patient interval when the gate set has >= 50 patients (Gate B1, whole val, 78
patients)". The DiD rule (far gain - near gain) has no CR1 form in A1 step 5. Fix: "The decision interval is the A1
gate CI of each rule (percentile patient interval for Gate B1, 78 patients; `cluster_t` below 50), the same draws as
the gate."

### 5. MINOR. The far-miss point rule sits at the design effect, so the B1.9 outcome is already known

Files: B (B.6.2 Gate B1, B1.9).

B.6.2: "val misses at offset > 0.3: point <= 20.5 % (25 % below R5's 27.3 %)"; B1.9 simulates exactly a 25 % cut.
[measured: `(planning probe, not kept: r5data)` on R5 val items, offset approximated as |offset| / side: 404 far
items, 112 far misses (27.7 %); 2,000 runs, conversion by patient] P(point rule passes) = 0.43 (noise 0), 0.63
(c = 0.077), 0.22 (c = 1), before the CI part and the other rules. So joint P(pass) at the 25 % effect is near or
below 0.5 by construction, and the "owner decides whether to change a rule" branch is close to certain. Fix: say this
in B1.9 now, and offer the owner the choice before B1.9 runs: (a) keep the point rule and expect the
smallest effect with P >= 0.8 to be above 25 %; or (b) make the point value a report line and keep only the CI rule
"CI upper bound of (B1 - R5) far-miss share < 0".

### 6. MINOR. DG14 "14.1 % overall" is not what the per-cell quota gives

Files: B (B1.2, B.6.1 DG14).

B1.2: "The quota with w = 4 gives q = 4p / (1 + 3p) = 14.1 % at p = 3.9 % in every block"; DG14: "14.1 % overall at
weight 4". q is computed per (kind, cell) from that cell's p, and disagreement items are large (share of
disagreement windows by size bin, train shifted draw 0: 0.0 / 0.8 / 7.0 / 19.3 %). [measured: geom.json joined with
`data/prepared/lidc/train.jsonl` `disagreement`, draw-0 pools, target = train shifted law] Overall share = 13.7 %;
the shifted kind 12.1 %, the rescaled kind 24.9 %, centred 14.1 %. Fix: DG14 text "about 13.7 % overall at draw 0
[measured]; per kind printed"; B1.2 "14.1 % in the centred set; per cell q(p_cell)".

### 7. MINOR. DG15 cannot pass on the B1.0 "rebuild differs" branch

Files: B (B1.0 step 2, B.6.1 DG15).

B1.0 step 2: "On a difference: keep both folders, freeze the rebuilt view ... mark A3-val ... as the only Gate B1
baseline." DG15 (hard): "the val shifted entry equals the B1.0 pre-rebuild hash". On that branch DG15 fails forever
and blocks V1. Low risk: 40 of 40 rebuilt items match today (see "Claims verified"). Fix: DG15 "equals the B1.0
pre-rebuild hash, or, on the B1.0 step 2 difference branch, the rebuilt hash recorded with the reason in
`frozen_views.json`".

### 8. MINOR. DG12 adds the end token twice

Files: B (B1.3 `tokens`), A (A8 step 5).

B1.3: "`tokens(rows_root, tokenizer, limit)`: certified target tokens (+1 end token) ... Rule (DG12): `max + 1 <
limit - 100`". A8 step 5: "1,182 includes the end token, so no "+1" is added again." No outcome change (1,183 <
1,436), but two codings differ. Fix B1.3: "Rule (DG12): `max < limit - 100`, where max counts the end token (A8)".

### 9. MINOR. Stale "9,600 rows"

File: 00_front, "2.4 GPU-everywhere data path" diagram line: "-> B1 rows root (9,600 rows)". V1 builds 8,320 rows
(B1.2: "V1 budgets 260 steps ..., so the B1 rows root holds 8,320 rows"). Fix: "(8,320 rows)".

### 10. MINOR. B-S2 has three item lists

Files: 00_front (6.3), B (B.4, B.8.1).

Front 6.3 table: "**B-S2 = items 0 + 4**"; front 6.3 bullet: "B-S2 in B.8.1 is V3 items 0, 2a and 4"; B.4 step 4 and
B.8.1: "rental V3 items 0, 2a, 2b and 4". Fix: front 6.3 bullet "B-S2 = V3 items 0 and 4, plus the skippable items 2a
and 2b when they run"; B.8.1 and B.4 copy that sentence.

## Counts

Blocking 2 (items 1, 2). Major 1 (item 3). Minor 7 (items 4-10).
