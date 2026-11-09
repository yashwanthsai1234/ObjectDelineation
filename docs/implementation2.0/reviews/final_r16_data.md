# Final review round 16: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in the data scope. B1.7's DICOM helper stays "recommend only". No catalog lookup, no install. |
| engineering:code-review | I read the code behind each claim I used: `build_eval_views._record` and `build` (draw order, `kept_side`, `fallback`), `augment.{shift,draw,neighbour_side}`, `dataset._prepare/_verified/prompts/cache_loader` (duplicate window ids, per-window `disagreement`), `splits.{deal,lidc_stratum}`, `vast/bundle.members` (which eval views ship), `vast/policy.code_files`, `final_eval.main` defaults (K 8, seed 11), TRL `RepeatSampler` call. |
| engineering:debug | I traced the Mac order "bundle, then `make preflight`" (A3 header) into the DG15 hard gate that preflight runs through `quality.failures`, and the E01 fixture call into the B1.2 block rule. |
| engineering:system-design | I checked which folders each gate reads (rows root, store, frozen views) against the folders that 6.1 step 4c moves before bundling. |
| engineering:architecture | One owner per rule: B1.4 manifest keys vs the E01 key list; 6.1 step 4c vs DG15. |
| engineering:documentation | Each item quotes the source line and the code line or measurement. |
| huggingface-skills:huggingface-best | No model choice in the data scope. Front section 3 unchanged. |
| huggingface-skills:huggingface-datasets | Web check of the Hub rate limits for the rows-root download (free user: 5,000 resolver requests per 5-min fixed window [web: https://huggingface.co/docs/hub/rate-limits, table "in September '25"]). This is the r12 fresh-arch-3 item, downgraded there; I do not re-raise it. No Dataset Viewer endpoint is on the data path. |

Web tools loaded (WebSearch, WebFetch). Web check, 2026-11-09: Hugging Face Hub rate limits page (above).

## Scope and state

- `phases/A.md` is unchanged since 02:40 (already read by rounds 13-15). `phases/B.md` and `00_front.md` changed
  at 03:47 for the round-15 fixes (data-2 keep regex, `val_half` path; arch-1). I read the round-15 changes for
  data effects, then the parts no earlier data round traced: the B.9 end-to-end tests for B1.4, and the folders
  that the Mac gates read after 6.1 step 4c.
- Read: A1, A2, A3, A4, A5, A.6-A.8; B.3.4, B.4, B1.0-B1.9, B.6.1, B.6.2, B.7, B.8.3, B.9.1-B.9.4 (B1 rows), B.9.8
  `test_data_windows.py` sketch; front 1, 1.1, 3, 6, 6.1, 6.2, 6.4, 7; 99_back "Data and quality"; greps of C, D,
  E, DOC, CLEAN for the B1 data files and the keep rule.
- Scripts (`.venv/bin/python -I`, 2026-11-09) in `(planning probe, not kept: r16data)`: `dg7.py` (train shifted target law and the
  per-block largest-remainder counts), `resc.py` (per-cell pool size against the 260-step need), a 40-item val
  rebuild into `r16data/views/`.

## Claims verified (no item)

- Val view rebuild is byte-stable: rebuilding the first 40 `lidc-val-shifted` items today gives equal offsets,
  `image_sha256`, targets and cache digests for 40 of 40 [measured: `build_eval_views --split val --limit 40` into
  `r16data/views`, compared with `data/eval_views/lidc-val-shifted/lidc/val.jsonl`]. So B1.0 step 2's "equal
  content hash" branch is the likely one.
- Val halves: s = 78 gives 40 / 38 patients, 533 / 551 items, 24 / 20 disagreement items; strata 41 / 23 / 2 / 12
  [measured: seed search re-run with `splits.deal` and `lidc_stratum` on the val view records]. As B1.0 states.
- Disagreement: 100 of 1,247 train nodules have a disagreement slice, all 100 also have agreeing slices, 1,394
  items in them; val 44 of 1,084; val items per patient 2-78 (median 8), top 5 patients 24.9 % [measured]. Per-window
  `disagreement` equals the source item's flag on all 7,341 windows of both `lidc-train-shifted` and
  `lidc-train-rescaled` [measured], so B1.1's "every window copies it" holds after `_prepare` recomputes it.
- Raters: every LIDC item has k >= 2 (train 1,933 / 2,029 / 3,379 for k = 2 / 3 / 4; val 281 / 330 / 473; test
  633 / 502 / 1,065) [measured: store `raters`]. So `dropped["not_rl"]` is 0 on the real store, and A2's k - 1 >= 1
  holds on val and test.
- Far / near subsets on whole val: 406 items in 75 patients (far), 206 items in 52 patients (near); far miss share
  0.2734; near IoU 0.7511; found IoU 0.7056 [measured: `r15data/far.py`]. Both subsets have >= 50 patients, so the
  percentile interval applies to both sides of the difference-in-differences rule.
- Feasibility of the schedule (B1.4 step 2), against the train shifted target law (7,341 draw-0 windows): per-block
  counts (640 rows) are 480 shifted and 64 rescaled; with `draws` 1 the only short cell is (rescaled, size 2,
  offset 0): 76 non-padded windows against 13 x 6 = 78 needed; every shifted cell has >= 1.10 x the need
  [measured: `r16data/resc.py` on `docs/implementation2.0/notes/b1_out/geom.json`]. So `draws` 2 is very likely the B1.4 value.
  The combined off-centre worst cell ratio after per-block largest-remainder rounding is 0.961 (TV 0.006)
  [measured: `r16data/dg7.py`], so B.6.1's "DG7-DG9 pass by construction" holds.
- Retest noise pair: Pearson 0.9712, mean difference -0.0076, SD 0.0569, variance ratio 0.0768 [measured:
  `r4data/retest.py` re-run]. B1.9's c = 0.077 and B.6.2's numbers hold.
- Offsets: test median 0.262, 36.7 % (807) above 0.3; val 0.257, 37.5 % (406) [measured: geom.json]. Front
  section 1 and B1.6 hold.
- `final_eval` defaults K = 8, seed 11 [code: experiments/phase2/final_eval.py:54-55], so `a3_measure.yaml` needs no
  `--seed` to match "K = 8, seed 11".

## Items

### 1. MINOR. E01 calls `build_rl_rows` with `steps 3`, which B1.2 refuses, and asserts the stale manifest key `target_half "shape"`

Files: B (B.9.4 E01; B1.2; B1.4 step 4).

Doc lines:
- B.9.4 E01: "`build_rl_rows` on a 30-40 item tiny store (steps 3, per-step 4 or 8, draws 2): exit 0; ... manifest
  keys exact, `target_half "shape"`, `ids_sha256 == sha256("\n".join(ids))` ...".
- B1.2: "`steps` must be a multiple of `block_steps`; else `ValueError`." The `build_rl_rows` command line of B1.4
  has no `--block-steps` flag, and B1.4's test list says "`--steps 30` with `block_steps` 20 is refused".
- B1.4 step 4: `schedule = {steps, planned_steps, block_steps, per_step, seed, draws, shares, target_law:
  "train-shifted", target_sha256, ids_sha256, ...}`. No `target_half`. The only other `target_half` in the sources is
  `B.spec.md:171-172` (revision 2: "`target_half = "shape"`", the old val-shape target).

Why minor: E01 fails on the Mac at $0 and names the problem (a `ValueError` on the block rule, or an extra key).
But an engineer who "fixes" E01 by adding `target_half: "shape"` writes the old val-shape target name into the
record of a run whose target is the train shifted law.

Fix:
- B1.4: add `--block-steps N` (default 20, recorded as `schedule.block_steps`) to `build_rl_rows`, used by test
  fixtures only; `b1_rows.yaml` and I22a keep the default.
- E01: "(steps 3, `--block-steps 1`, per-step 4 or 8, draws 2)"; replace "`target_half "shape"`" by
  "`target_law "train-shifted"`, and the schedule keys exactly as B1.4 step 4 lists them".

### 2. MINOR. DG15 needs the frozen test view folder at `make preflight`, but 6.1 step 4c moves that folder out before bundling, which comes first; no gate names where it reads the frozen views

Files: 00_front (6.1 step 4c), A (A3 header comment block), B (B1.3 `rows_leakage`, `audit_data`, preflight change
2; B.6.1 DG15).

Doc lines:
- 00_front 6.1 step 4c: "The bundle ships every `data/eval_views/lidc-*` folder ..., so the test view is moved out
  before bundling (or a bundle `--view` filter is added)." No target folder is named.
- A3 header order (A.md lines 565-570): `check_series_z`, `build_eval_views --split test`, `final_eval --hash-view`,
  "B1.6 builds lidc-val-edge", `vast.bundle lidc --view lidc-bundle-v3`, then `make preflight
  CONFIG=experiments/phase2/configs/b1_launch.yaml`.
- B1.3 preflight change 2: "`preflight` ... calls `quality.failures` on the root. A hard failure ... is a failed
  check ("ROWS root")". U19: DG15 is hard (only DG9b, DG14, DG14b, DG16 are `hard=False`).
- B.6.1 DG15: "Frozen views (test shifted, val edge, val shifted): content hash ... and target-set hash equal
  `docs/reference/frozen_views.json`". The content hash needs each view's records (B1.6: tuples `(item_id,
  view.offset, view.side, image_sha256)`).
- B1.3: `audit_data.py --rows <root> --store data/store --target-law train-shifted [--val-view ...] [--json ...]`:
  no argument for the frozen views. B1.3 `rows_leakage`: "no row `image_sha256` equals any val or test view image
  hash or a frozen test view hash (`docs/reference/frozen_views.json`)". That file holds view-level
  `content_sha256` and `target_set_sha256`, never per-image hashes (B1.6 field list), so that comparison can never
  match.

Code: `bundle.members` adds every `data/eval_views/lidc-*` folder for LIDC [code:
experiments/phase2/vast/bundle.py:53-54, 65]; it has no exclude filter today.

Why minor: preflight fails on the Mac at $0 and names DG15 (caught before spend). The engineer must still guess
where the test view went and where DG15 looks.

Fix:
- 6.1 step 4c: name one hold folder, `data/eval_views_hold/` (the test view is moved there before bundling and
  stays there), or add `bundle --exclude-eval 'lidc-test-*'` so nothing moves (preferred: no folder moves).
- B1.3: `audit_data` and the preflight "ROWS root" check read the frozen views from `--frozen-root` (default
  `data/eval_views`, then `data/eval_views_hold` when the hold folder is used); a missing frozen view is a DG15
  failure that names the folders searched. B.4 task order adds "B1.6 before B1.4" (DG15 needs all three entries
  when B1.4 runs its gates).
- B1.3 `rows_leakage`: replace "or a frozen test view hash (`docs/reference/frozen_views.json`)" by "(the test view
  image hashes come from the frozen test view folder; `frozen_views.json` holds view-level hashes only)".

### 3. MINOR. B1.4 step 2 still says the passing `draws` "is not yet measured"; the pool numbers now bound it

File: B (B1.4 step 2). Re-raises r9 data-4 with new numbers, not re-argued.

Doc line: "(Today one draw fails [measured: SD3 s5]; the value that passes is not yet measured.)"

Measurement (above): with `draws` 1 only (rescaled, cell (2, 0)) is short, 76 against 78; all shifted cells have
>= 1.10 x the 260-step need; the off-centre DG7 worst ratio after per-block rounding is 0.961 [measured:
`r16data/resc.py`, `r16data/dg7.py`]. The disagreement-reuse condition of step 2 is not measured here.

Fix: "Today `draws` 1 fails only at (rescaled, size 2, offset 0): 76 windows against 78 needed; `draws` 2 is the
expected value [measured: final_r16_data]; the step-2 table decides."

### 4. MINOR (still open; not re-argued). Earlier data items not yet applied

Re-checked in the sources on 2026-11-09; the text is unchanged:
- r15 data-1 (downgraded to minor): front 1.1 row "B1 stopped early" still keys on "`rl.json` has `partial: true`",
  which a `BudgetStop` return never writes.
- r15 data-3: B1.3 still says `rows_leakage` runs "again by `preflight` on the downloaded root" on the box.
- r15 data-4: B1.0 still says "Nothing in training reads val" beside the optional val-shape held-out curve.
- r14 data-1, data-2; r13 data-1 to data-8 (B1.7 "so the count can grow" is still at B.md line 1086; DG14b
  timing and the missing "drop" / "weight 0.5" mechanism of B.8.3 item 4; "about 39 patients" is 38 [measured:
  item above]); U24 still says "cosine 0.98".
- r12 data-5: Gate B1 "point <= 20.5 % (25 % below R5's 27.3 %)" still does not say it is absolute.
- B1.1 pool test "S0 windows are >= 0.75 of off-centre windows" (r4 data-7, r6 data-5).
- A8 step 5 "B1 then reuses these views and only audits them" (A.md line 1137) and the `rl_shifted.yaml:10` build
  command (r10 data-2, r8 data-4).
- DG12 end token counted twice (r5 data-8).

No blocking or major item found in the other phases in the parts I read (front 6.x, 99_back data decisions, the
B1 data-file rows of C, D, E, DOC, CLEAN).

## Counts

Blocking 0. Major 0. Minor 4 (items 1-4; items 1 and 2 are new).
