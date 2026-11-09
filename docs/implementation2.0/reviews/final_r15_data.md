# Final review round 15: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice is in the data scope. B1.7's DICOM helper stays "recommend only". No catalog lookup, no install. |
| engineering:code-review | I read the code behind each claim I used: `train_rl.py:19-68` (what `rl.json` holds after a normal return), `vast/policy.py:374-383` (`code_files` drop rule), `vast/bundle.py:48-73` (what the bundle ships; `lidc-train-*` views only when named), `build_eval_views._record`, `augment.{shift,draw,neighbour_side}`, `dataset._prepare/prompts/grpo_rows`, TRL `RepeatSampler` call (`grpo_trainer.py:1372-1377`) and `shuffle_dataset` (`grpo_config.py:509`). |
| engineering:debug | I traced two paid-run paths to their end: (1) `BudgetStop` trips before step 260 -> `trainer.train` returns -> `final/` and `rl.json` without `partial` -> front 1.1 row; (2) `--heldout` on the box -> `val_half` -> `docs/reference/val_halves.json` -> dropped by the code sync. |
| engineering:system-design | I followed the files that a box process reads (code sync, bundle, HF rows root) against the files B1 writes under `docs/reference/`. |
| engineering:architecture | One owner per rule: front 1.1 "B1 stopped early" vs B1.5 step 3b vs B3.9 "same `--steps`"; B1.6 keep rule vs every `docs/reference/*.json` that a box step reads. |
| engineering:documentation | Each item quotes the source line, the code line or the measurement. |
| huggingface-skills:huggingface-best | No model choice in the data scope. Front section 3 unchanged. |
| huggingface-skills:huggingface-datasets | Checked the rows-root download path: `snapshot_download(local_dir=/workspace/data/lidc, allow_patterns="rows/lidc-train-b1/**")` keeps its `.cache/huggingface/` folder at the `local_dir` root, outside `rows/lidc-train-b1/`, so the box file-count check of B1.5 is not polluted. No Dataset Viewer endpoint is on the data path. |

Web tools loaded (WebSearch, WebFetch). Web check, 2026-11-09: TCIA LIDC-IDRI collection page
[web: https://www.cancerimagingarchive.net/collection/lidc-idri/]: license "CC BY 3.0" (B1.4 step 5 is right);
"Annotation and Markup Issues" item 6 is the LIDC-IDRI-0396 (139.xml) SOP Instance UID fix of 2018-06-28 (B1.7
check (c) cites it correctly); no TCIA item mentions duplicate z positions (B1.7 names a third-party card for that,
correctly).

## Scope and state

- `phases/A.md` is unchanged since 02:40 (before the round-13 and round-14 data reviews). `phases/B.md` changed at
  03:28 for round-14 items rl-2 (B3.9 soft `train_rl`, step 3b) and rl-3 (OPSM start delta). I read the round-14
  changes for data effects, then looked for what no earlier data round traced: what a paid box process reads, and
  how a run that stops early is labelled.
- Read: A1, A2, A3, A8; B.3.4 ADR-0015, B.4, B1.0-B1.9, B3.5b, B3.9, B.6.1, B.6.2, B.7, B.8.3; front 1.1, 6, 6.1,
  6.2, 6.4; 99_back "Data and quality"; E and CLEAN greps for the B1 data modules and the keep rule.
- Scripts (`.venv/bin/python -I`, repo path as argument, 2026-11-09) in `(planning probe, not kept: r15data)`: `counts.py` (items,
  patients, nodules, series per split), `far.py` (val far/near subsets, miss share, found IoU).

## Claims verified (no item)

- Store counts: train 7,341 items / 557 patients / 1,247 nodules / 562 series; val 1,084 / 78 / 200 / 78; test
  2,200 / 158 / 401 / 160 [measured: `counts.py` on `data/store/lidc/*.parquet`]. They match ADR-0015, Gate B1 and
  A1. Five train patients and two test patients have two series; the split is by patient, so no physical nodule
  crosses splits.
- Gate B1 numbers: 406 val items with offset > 0.3, far miss share 0.2734, 206 items with offset < 0.15 (52
  patients) at IoU 0.7511, found IoU 0.7056, val IoU 0.6182 [measured: `far.py`, offsets from
  `lidc-val-shifted` record `box` and `augment.union_box`, values from `runs/r5/eval/eval_val/items.jsonl`]. So the
  decomposition 0.25 x 0.273 x 406 = 28 and 28 x 0.706 / 1,084 = +0.018 holds, and r9 data-3's 52-patient near
  subset holds.
- The test view and the train views draw shifted windows with `side=S0` (`build_eval_views._record` passes
  `side=S0` to `augment.draw` [code: build_eval_views.py:56-58]), so ADR-0015's "the train law equals the test
  view's generator law" holds for the side.
- B1.1 bin edges: `np.searchsorted((0.15, 0.2, 0.3), 0.15, side="right") = 1` and apparent 15 px -> size bin 1, as
  the B1.1 tests state.
- The bundle never ships the rows root: `members()` drops `lidc-train-*` views unless named [code: bundle.py:55-56],
  and the A3 header bundles with `--view lidc-bundle-v3`. So B1.4 step 5's HF upload is the only path.
- B1.5's unshuffled order: `shuffle_dataset` is a public `GRPOConfig` field [code: .venv trl/trainer/grpo_config.py:509]
  and TRL passes it to `RepeatSampler(shuffle=...)` [code: grpo_trainer.py:1372-1377]. 8,320 rows / 32 = 260 whole
  generation batches.

## Items

### 1. MAJOR. A `BudgetStop` that ends the B1 run before step 260 writes no `partial: true`, so front 1.1's "B1 stopped early" row never fires

Files: 00_front (1.1; 6.2 notes), B (B1.5 step 3b and the `train_rl` bullets; B3.9).

Doc lines:
- 00_front "1.1 After V1: decision table", row: "B1 stopped early: `rl.json` has `partial: true` (a crash, then
  `train_rl --merge-last`, 6.2) and N < 260 steps | Gate B1 is computed on the merged checkpoint-N and labelled
  "partial run, N steps". The row the rules give applies, but a pass does not lead to B-S2 until the owner reviews
  it".
- 00_front "6.2", notes, "Why 260 steps": "With the 53-min start the stop trips only above about 27.5 s per step
  ... If it trips, the adapter is saved and items 3-5 still run." The planning step time is itself 27.5 s
  (section 6 inputs: "Step time, base | 27.5 s ... Shifted rows may lengthen the decode tail"), and "The rows-root
  download is not in the measured rows". So a trip is a planned case, not a corner case.
- B "B1.5", step 3b: "If `{rd}/b1-rows/final/` exists (the run ended normally, or `BudgetStop` ended it), it does
  nothing and exits 0."
- B "B1.5", `train_rl` bullets: "Also require `--steps <= schedule.steps` and `--steps % schedule.block_steps == 0`
  (only whole blocks are trained)."
- B "B3.9": "the same `--steps` as the B1 run (260, front 6.2)".

Code: `BudgetStop` sets `should_training_stop` (U46), so `trainer.train` returns normally. Then `train_rl` saves
`final/` and writes `rl.json` with `"steps": done` and no `partial` key [code: experiments/phase2/train_rl.py:55-68].
Only `--merge-last` writes `partial: true` (B1.5 step 3b, I49). I43 tests the row only with a record that has
`partial: true` (B.9.3 I43).

Failure: the B1 run trips `BudgetStop` at, for example, step 237 (inside a 20-step block). `rl.json` has
`steps: 237` and no `partial`. The front 1.1 row needs `partial: true`, so the evaluator applies the plain row. A
pass funds B-S2 with no owner review, and the result is not labelled "partial run". B-S2 then trains B3 for 260
steps (B3.9 "the same `--steps` as the B1 run (260)") and Gate B3 pairs B3 (260 steps) with B1 (237 steps): 23
extra steps that Gate B3 credits to v5 (B3.9 itself says "A lower value cuts the 260 steps, and Gate B3 then
compares unequal runs"). No planned $0 test covers this path: I43 and I49 use `partial: true`, U46 tests only the
guard flags.

Fix:
- B1.5 `train_rl`: on every exit path write `partial = (global_step < --steps)` and `stopped_by` (`"budget"`,
  `"loc_guard"`, `"crash"` from `--merge-last`, or null) into `rl.json`. U46 gains a case: a fake clock that trips the
  guard at step 3 of a 4-step stub run gives `rl.json` `partial: true, stopped_by: "budget", steps: 3`.
- Front 1.1 row: "B1 stopped early: `rl.json` `steps` < 260 (a crash and `--merge-last`, or `BudgetStop`)". I43
  gains the case `steps: 237, partial: true, stopped_by: "budget"`.
- B1.5 "only whole blocks are trained": add "except a `BudgetStop` trip; the label above then applies" (this also
  closes r9 data-2).
- B3.9: "`--steps` = the B1 run's `rl.json` `steps` rounded down to a multiple of 20 when B1 stopped early (owner
  review of the 1.1 row decides), else 260".

### 2. MAJOR. The optional held-out curve reads `val_halves.json` on the box, but the code sync drops it; a B1 run with `--heldout` fails at start

Files: B (B1.0; B1.5 "Held-out curve"; B1.6 "Ship the record to the box"; I44), CLEAN (CLEAN1 keep regex).

Doc lines:
- B "B1.0", val halves: "Its output is stored as `docs/reference/val_halves.json` ... `val_half` is then a lookup in
  that table, needs no store at run time, and raises `KeyError` on an id that is not a val patient."
- B "B1.5", "Held-out curve (optional, owner decision B.8.3 item 15 ...)": "`train_rl --heldout <view>
  --heldout-items 64 --heldout-every 50` ... The eval set is 64 items of **val-shape** only, so val-select stays
  clean."
- B.9.3 I44: "the eval set holds val-shape items only (a planted val-select id is refused)". So `train_rl` on the
  box must call `val_half`, which reads `docs/reference/val_halves.json`.
- B "B1.6", "Ship the record to the box": "The code sync drops all of `docs/` ... So B1.6 adds the one keep rule of
  `code_files()` ...: `keep = re.compile(r"^docs/reference/frozen_views\.json$")`". CLEAN's form
  (`CLEAN.md` line 1986): `kept_docs = re.compile(r"^docs/reference/frozen_views\.json$|^docs/implementation2.0/kernels/K[123]/")`.

Code: `drop = re.compile(r"^(data|docs)/|...")` and `files = [f for f in files if not drop.search(f)]`
[code: experiments/phase2/vast/policy.py:379-380]. `docs/reference/` today holds no JSON at all [measured: ls].

Failure: the owner turns on item 15 for V1 (B.8.3 says it costs $0.02). `train_rl` starts on the box, builds the
held-out set, calls `val_half`, and the file is missing (`FileNotFoundError`). Step 3 is soft, step 3b finds no
checkpoint and exits non-zero, the `&&` chain stops [code: policy.py:490]: no B1 model, no Gate B1 eval, and about
53 min ($0.70, front 6.2 item 1) plus the start-up are spent for nothing. I44 runs on the Mac, where `docs/`
exists, and I50 checks only `frozen_views.json`, so no planned $0 check catches it.

Fix (pick one, state it in B1.5):
- (a) B1.6 keep rule becomes `^docs/reference/(frozen_views|val_halves)\.json$` (CLEAN1's `kept_docs` and D1 step 0
  widen the same regex); I50 also finds `val_half("LIDC-IDRI-0026")` in a child `python -I` on the `code_files()`
  copy. The file holds patient ids and item ids, no pixels, so front 6.4 rule 7 still holds; or
- (b) the Mac resolves the 64 held-out ids before the rental (`experiments/phase2/heldout_ids.py`, val-shape only,
  patient-stratified, seed 0) into the plan step as `--heldout-ids <file under the rows root>`, and `train_rl` on the
  box reads no `docs/` file. I44 then checks the id file, not `val_half`.

### 3. MINOR. B1.3 says `preflight` re-runs `rows_leakage` "on the box"; no plan step runs preflight on the box, and the box has no store

File: B (B1.3 `rows_leakage`).

Doc line: "Called by `build_rl_rows` before the manifest is written ..., and again by `preflight` on the downloaded
root. The code sync drops `docs/` ..., so on the box `preflight` passes the manifest's `schedule.excluded_ids`".
`make preflight` runs on the Mac (front 6.1 step 4); `b1_rows.yaml` steps 0-4 have no preflight step; the bundle
ships no `data/store` and no test view [code: bundle.py:48-73], and `rows_leakage(rows_root, store_root, ...)` needs
both. The box integrity check is the manifest SHA and file count of B1.5 "Rows-root download", which already pins
the Mac-checked root.

Fix: replace "and again by `preflight` on the downloaded root ... on the box `preflight` passes the manifest's
`schedule.excluded_ids`" by "and again by `preflight` on the Mac on the local root (front 6.1 step 4). On the box,
the manifest SHA-256 and file count (B1.5) prove that the root is the one the Mac checked." Keep `excluded_ids` in
the manifest as the record.

### 4. MINOR. B1.0 says "Nothing in training reads val", but the optional held-out curve evaluates val-shape items during the B1 run

File: B (B1.0 last bullets; B1.5 "Held-out curve").

Doc lines: B1.0 "Val-shape tunes B3 (B3.7). Val-select judges Gate B3. Gate B1 uses the whole val. Nothing in
training reads val." B1.5: "The eval set is 64 items of **val-shape** only ... Report only, never a stop rule."

Why minor: the curve changes no weight and no stop rule. But Gate B1 is computed on the whole val, and the B-S2 step
count can be chosen from this curve (B1.5 "without a curve the B-S2 step count is chosen blind"), so val-shape then
informs a training choice.

Fix: B1.0: "No gradient step reads val. The optional held-out curve (B1.5) reads 64 val-shape items as a report; if
the owner uses it to choose the B-S2 step count, the Gate B3 report says so."

### 5. MINOR (still open; not re-argued). Earlier data items not yet applied

Re-checked in `phases/B.md` and `phases/A.md` on 2026-11-09; the text is unchanged:
- r14 data-1 (DG3 drops one slice, not the nodule; "cos" wording), r14 data-2 (real-cluster coverage of the gate
  intervals, about 91.5 % for the miss share).
- r13 data-1 to data-8 (B1.7 "count can grow", 0085 / 0146 named as val hits in B1.7, DG17 "Today" and 99_back 8c;
  "mapped" cannot clear the four (d) hits; 998 sidecars; "139.xml" ambiguity; B1.0 hashes before B1.7; DG17 mapped
  wording; DG14b input timing; "about 39 patients" is 38).
- B1.1 pool test "S0 windows are >= 0.75 of off-centre windows" (r4 data-7, r6 data-5): the pool holds d shifted
  (S0) and up to d rescaled (non-S0) windows per item, so the share is near 0.5; it is a schedule property (DG10).
  The test fails on the Mac at $0.
- A8 step 5 "B1 then reuses these views and only audits them" and the `rl_shifted.yaml:10` build command (r10
  data-2, r8 data-4): B1.4 builds a new rows root, and B1.5 archives `rl_shifted.yaml`.
- DG12 end token (r5 data-8): B1.3 `tokens()` already adds the end token, and DG12 adds "+ 1" again; A8 step 5 says
  1,182 already includes it.

## Counts

Blocking 0. Major 2 (items 1, 2). Minor 3 (items 3, 4, 5).
