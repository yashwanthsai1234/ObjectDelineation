# Final review round 6: data and statistics (A and B1 first)

Reviewer role: medical imaging data scientist and statistician (LIDC multi-rater, leakage, cluster bootstrap).
Scope: the whole assembled doc, read by section from the sources. I read A.1, A1-A5, A8, A.6-A.8, B.4, B1.0-B1.9,
B.6.1, B.6.2, B.8.3, the B.9 unit rows for B1.1/B1.2, front 1 and 1.1, front 6.2, 99_back "Owner decisions", and
the E/CLEAN parts that touch B1's code. Round-5 fixes in `review-log.md` (data-1, data-2, data-3) are not repeated.
Scripts for this round: `(planning probe, not kept: r6data) offlaw.py}`.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product decision in the data scope (the B1.7 DICOM helper is "recommend only"). No catalog lookup or install. Nothing changed. |
| engineering:code-review | Read the code that each data claim cites: `augment.{row_rng,shift,draw,neighbour_side}`, `build_eval_views._record`, `dataset._prepare/_verified/prompts`, `serve.loader/_tagged/VIEW_SEP`, `train_rl.py:26,47`, `g_probe.py:90`, `splits.{deal,lidc_stratum}`, `frame.crowded`, `lidc.py` P1/P6, `configs.grpo_config`, `hf_pull.py`. Found items 2 and 12. |
| engineering:debug | Reproduced numbers before using them: val far/near subsets from target union boxes (406 far items, 75 patients, 27.3 % misses; 206 near items, 52 patients, IoU 0.751: all equal the doc); patient / series disjointness of the store; disagreement flags across centred, shifted and rescaled roots (0 differences); disagreement per nodule (item 1). |
| engineering:system-design | Traced the rows-root path Mac build -> HF upload -> box download -> `train_rl` and found where a partial revision is caught only on the paid box (item 3). |
| engineering:architecture | Checked which text owns the post-Gate-B1 cleanup (B1.8) against the later phases that still name the deleted code (CLEAN2 CT-33/CT-35, E3b move row, E-G4 U33, E6 `train_views`) (item 2). |
| engineering:documentation | Checked each test row against the rule it claims to test (U09, U11, the B1.1 test list) and each owner-decision copy against front 6.2 (items 1, 4, 5, 6). |
| huggingface-skills:huggingface-best | No model choice in scope; front section 3 left as is. |
| huggingface-skills:huggingface-datasets | Used for the HF side of the rows root: read the Hub storage-limit guide (files per repo/folder, "merge data into fewer files", Parquet/WebDataset) and the installed `huggingface_hub` 2.2.0 `upload_folder` docstring (large folders are split into several commits). Items 3 and 10. |

Web checks: TCIA LIDC-IDRI page (items 6-8 of "Annotation and Markup Issues": 0396 SOP UID fix 2018-06-28, 0510
invalid `internalStructure`, 8 two-timepoint patients; no duplicate-patient item) [web:
https://www.cancerimagingarchive.net/collection/lidc-idri/]; HF storage limits [web:
https://huggingface.co/docs/hub/storage-limits]; a search for LIDC cross-patient duplicate scans found none
documented.

## Claims verified (no item)

- DG1: store patients 557 / 78 / 158, series 562 / 78 / 160, items 7,341 / 1,084 / 2,200; 0 patient and 0 series
  overlap between any two splits [measured: `data/store/lidc/*.parquet`, `python -I`]. Every store item has an
  image and k >= 2.
- Gate B1 subsets on the frozen-to-be val view (offset from the doctors' union box in target coordinates):
  far 406 items / 75 patients / 27.34 % misses; near 206 items / 52 patients / IoU 0.7511. Equal to B.1 and B.6.2.
  Both subsets have >= 50 patients, so the percentile interval applies to them under the A1 rule.
- Disagreement flags are identical in `data/prepared`, `lidc-train-shifted` and `lidc-train-rescaled` (289 each,
  0 differences), so a window inherits its source item's flag (supports B1.2's p).
- Rescaled supply under 99_back 7b option (i): every cell with demand >= 0.5 row has enough draw-0 windows at
  draws >= 2 for 260 steps (tightest: cell 20, need 76.4 vs 77 draw-0 windows) [measured: `r6data/resc_feas.py` on
  `notes/b1_out/geom.json`].
- Offset law: train shifted median 0.247, > 0.3 32.5 %; val 0.257 / 37.5 %; test 0.262 / 36.7 % (train nodules are
  larger: median extent 16 vs 14 px) [measured: `r6data/offlaw.py`]. Consistent with DG16's note.
- A1/U53/U55/U20 hand values recomputed (CR1 equal-N_p form = sd/sqrt(P); item mean 0.6875; d mean 0.05;
  legacy item mean 0.7): all correct.
- `grpo_config`: 16 x 16 completions = 32 prompts x 8 per step, so `per_step=32` matches TRL's sampler formula.

## Items

### 1. MAJOR. U11 defines `Window.disagreement` per nodule; the rule and every share in B1.2 are per slice

Files: B (B1.1, B1.2, B.9.1 U11).

Doc lines:
- B, "B.9.1 Unit tests", U11: "`Window.disagreement` is True iff the doctors' topologies disagree (the repo's rule
  for prepared records); equal across all windows of one nodule."
- B, "B1.1": "`@dataclass(frozen=True) class Window`: `item_id, group_id, nodule, kind, draw, ...,
  disagreement`" (the nodule is a separate field from the slice item).
- B, "B1.2" Disagreement weight: "289 of 7,341 train items are disagreement items (3.9 %) ... The quota with w = 4
  gives q = 4p / (1 + 3p) = 14.1 %".

Proof: the repo's rule is per record (slice): `rec["disagreement"] = len(rater_answers(caches, names)) > 1`
[code: src/objdel/train/dataset.py `_prepare`]. On `data/prepared/lidc/train.jsonl`, 100 of 1,247 nodules have a
disagreement slice, and all 100 also have agreeing slices; those nodules hold 1,394 items [measured: `python -I`,
grouping `item_id.rsplit("/", 1)[0]`, 2026-11-08]. So "equal across all windows of one nodule" is false on the real
data. An engineer who makes U11 pass with a nodule-level flag marks 1,394 items (19.0 %), not 289 (3.9 %); the w = 4
quota is then 4 x 0.19 / (1 + 3 x 0.19) = 48 % disagreement rows instead of 14.1 % [arithmetic], and DG11, DG14 and
"B1 changes only the geometry against R5" all change.

Fix: U11 text: "`Window.disagreement` equals the source item's `disagreement` field (per slice, `dataset._prepare`
rule); equal across all windows (kinds, draws) of one **item**. Planted case: two slices of one nodule, one with
disagreeing topologies, give one True and one False." Add to B1.1: "`disagreement` is per item (slice), copied from
the source item's prepared record; 100 of 100 train nodules with a disagreement slice also have agreeing slices."

### 2. MAJOR. B1.8 deletes the multi-root path, but CLEAN wave 1 and Phase E still test and move it

Files: B (B1.8), CLEAN (CLEAN2 CT-33, CT-35, sketch F), E (E3b move row, U33 / E-G4, E6 `RlStage`).

Doc lines:
- B, "B1.8 Clean after gate B1": "delete `dataset.grpo_rows_views`, `serve.VIEW_SEP` and `_tagged` ... (keep
  single-root `loader(root)`), `--train-views` ... Update `tests/test_p2_views.py:610`,
  `tests/test_p2_train_rl.py:76-90`, `tests/test_profile.py:227`." It runs on the front 1.1 rows "Pass" and "Fails
  only the far-miss rule".
- CLEAN, "CT-33 --procs is an argparse error" (CLEAN2.1, Wave 1a, on the tag `v1-b1`): "`--model --ckpt --train
  --train-views --out --steps --save-steps` present"; sketch F (`tests/test_p2_train_rl.py`)
  `test_help_has_no_procs_and_keeps_the_other_options` asserts `"--train-views" in out`. CT-35: "Run with
  `--train-views a b` (3 preloads)".
- E, "E3b" move table: "`verify/store.py` | `train/serve.py`: `preload`, `CacheRef`, `VIEW_SEP`, `loader`,
  `_tagged`"; E.9.1 U33 (E3b, gate E-G4): "`verify.store` exposes `preload, loader, VIEW_SEP == "::"`;
  `loader(r1, r2)` takes the first root holding the key". E6 `RlStage`: "`train_views: tuple[str, ...] = ()   # extra
  view roots (train_rl --train-views)`".

Proof: today `--train-views` and the tagged loader exist [code: experiments/phase2/train_rl.py:26, 47;
src/objdel/train/serve.py:107-121; experiments/phase2/g_probe.py:90]. CLEAN2 lands before V1 (CLEAN "Wave 1 is part
of the owner tag `v1-b1`"), so after Gate B1 passes, B1.8's own list leaves CT-33 and CT-35 failing, and `make
test` must pass before the B1.8 commit (repo rule). E runs after Gate B1 (front time line step 7), so on the Pass
row E-G4's U33 tests names that no longer exist, and E6 ports a flag the CLI no longer has. The engineer must guess
whether to keep the multi-root loader (to satisfy E) or delete it (B1.8); the two texts give different code.

Fix:
1. B1.8: replace the line-number list with "update every test that names `--train-views`, `VIEW_SEP`, `_tagged` or
   `grpo_rows_views` (grep on the B1.8 commit; today `tests/test_p2_views.py`, `tests/test_p2_train_rl.py`
   including CLEAN2's CT-33 option list and CT-35, `tests/test_profile.py`)". CT-33's option list becomes "the
   pre-CLEAN literal minus `--procs` (and minus `--train-views` after B1.8)"; CT-35 uses two `--train` roots only if
   B1.8 has not run, else it checks one preload.
2. E3b row and U33: "On the E0 commit: if B1.8 ran, `VIEW_SEP` and `_tagged` are absent and `verify.store.loader`
   takes one root; U33 asserts the one-root form. If B1.8 did not run, the row and U33 stay as written." E6:
   "`train_views` exists only if `train_rl --train-views` exists on the E0 commit (E2 rule: fields come from the
   CLI on the E0 commit)."

### 3. MAJOR. No check before spend that the pinned HF revision holds the whole rows root

Files: B (B1.4 step 5, B1.3 preflight change 2, B1.5 "Rows-root download"), 00_front (6.1, 6.2 item 1).

Doc lines:
- B, "B1.4" step 5: "Upload the rows root before any rental ... Record the HF commit as `rows_revision` and the
  SHA-256 of `manifest.json` as `rows_manifest_sha256`".
- B, "B1.3" change 2: preflight "requires the local root's `manifest.json` SHA-256 to equal the flag, and calls
  `quality.failures` on the root" (the local folder only).
- B, "B1.5" Rows-root download: "The same ssh call then checks the box: `manifest.json` exists under that folder and
  its SHA-256 equals the plan's `--rows-manifest-sha256`".

Proof: the root is many small files: today's 7,341-window train view holds 14,684 files, 1.2 GB [measured: `find
data/views/lidc-train-shifted -type f | wc -l`], so 8,320 windows give about 16,640 files (one PNG and one cache
per window [code: dataset.py `_prepare`]). The installed `huggingface_hub` 2.2.0 `upload_folder` says "large
folders are automatically committed in several batches to stay below server limits (follow-up commits get a
` (part N)` suffix)" [code: .venv huggingface_hub/hf_api.py:5779-5784]; the HF guide recommends < 100 files per
commit [web: huggingface.co/docs/hub/storage-limits]. So "the HF commit" is one of many. If `rows_revision` names a
part commit, or an interrupted upload is pinned, nothing on the Mac sees it: preflight reads the local folder, and
the box check reads only `manifest.json`. The gap shows when `train_rl` starts and `_verified` opens a missing
image [code: dataset.py:178-182], after launch, sync, build, `env_check` and smoke (front 6.2 item 1: 53 min,
$0.70 of the $3).

Fix:
1. B1.4 step 5: upload with one `HfApi.upload_folder(...)` call and record `CommitInfo.oid` of its return value
   (the last commit). Write the local file list and count into `manifest.schedule.files_sha256` (SHA-256 of the
   sorted relative paths) and `files_count`.
2. New pre-rental Mac check in `make preflight` (B1.3 change 2b, front 6.1 step 4): `HfApi.list_repo_files(repo,
   revision=rows_revision)` under `rows/lidc-train-b1/` equals the local sorted path list (count and SHA), and
   `hf_hub_download(..., "rows/lidc-train-b1/manifest.json", revision=rows_revision)` has SHA-256
   `rows_manifest_sha256`. A mismatch is a failed "ROWS root" check. Token from `objdel.secrets`, $0.
3. Box check (B1.5): after the download, also compare the file count under `{data}/rows/lidc-train-b1` with
   `files_count`. Test: I40c with a fake Hub whose revision lacks one PNG fails the Mac check.

### 4. MINOR. U09 says `dropped` has 3 keys; B1.1 defines 4

Files: B (B1.1, U09). B1.1: "`dropped: Mapping[str, int]` (keys `"not_rl"`, `"shifted_none"`, `"rescaled_none"`,
`"padded"`)". U09: "`dropped` has exactly the 3 int keys even at 0". One of them fails. Fix U09: "exactly the 4 int
keys (`not_rl`, `shifted_none`, `rescaled_none`, `padded`) even at 0".

### 5. MINOR. B1.1 asks the pool for an S0 share it cannot have

Files: B (B1.1 tests). "Tests `tests/test_data_windows.py` ...: S0 windows are >= 0.75 of off-centre windows". The
pool makes one shifted (S0) and one rescaled (not S0) window per item per draw, so its S0 share is about 0.5
(train: 7,341 shifted and 7,341 rescaled draw-0 windows [measured: geom.json]). The 0.75 rule is DG10 on the
schedule (R03 tests it). Fix: move the sentence to B1.2's tests ("the schedule's S0 share of off-centre rows is
>= 0.75"); B1.1 tests "shifted windows have `side == S0`" (as U10 does).

### 6. MINOR. B.8.3 item 11 gives 260 steps for $2.20; front 6.2 and 99_back decision 2 give 220

Files: B (B.8.3 item 11). "Keeping the older $2.20 means a 260-step B1 run (front 6.2 notes)." Front 6.2 notes:
"If the owner keeps $2.20 anyway: ... the run plans 220 steps"; 99_back decision 2: "The older $2.20 would cut the
B1 run to 220 steps". Fix: "Keeping the older $2.20 means a 220-step B1 run (front 6.2 notes)". Also item 14: "The
decision table B.4.1" -> "front section 1.1".

### 7. MINOR. A3 step 5 still says "stop" on a failed view-identity rule

Files: A (A3 step 5). "A lower correlation means a different view: stop, and pair only on rebuilt, hashed views
(final r2 data item 3)." The next bullet and B.6.2 say Gate B1 still pairs with A3's items and drops only the R5
line. Fix: "A lower correlation means R5's box view differed from the frozen view: drop the R5 `items.jsonl` second
line (B.6.2). Gate B1 pairs with A3's items."

### 8. MINOR. Front 1.1 has no precedence when two rows match

Files: 00_front (1.1). Example: the IoU rule and the edge guard both fail: row 3 says "R5, unless the 4B arm wins
on val", row 4 says "Stop. Owner review ... R5". Fix, one line under the table: "If rules of two rows fail, the
lowest row wins (row 4 over row 3 over row 2)."

### 9. MINOR. Shifted windows can show nodules the doctors marked only as points

Files: B (B1.1, B1.3 `window_checks`, DG5). `crowded` refuses only boxes in `meta.other_nodules`, the nodules
>= 3 mm that a reader outlined [code: src/objdel/targets/frame.py:124-127; src/objdel/data/lidc.py P6]. Nodules
< 3 mm and >= 3 mm non-nodules are single points in the XML [code: lidc.py:3-4] and are not refused. A shifted
window shows more context than the centred frame, so it more often shows such a mark with no outline in the
target. Train and test views share this, so it is a label-noise report, not a leak. Fix: DG14 adds a report line
from the XML ($0): the share of scheduled windows (and of the val/test views) that contain a small-nodule or
non-nodule point of any reader.

### 10. MINOR. The rows-root download time is not measured, and 16,640 small files is the slow form

Files: 00_front (6.2 notes), B (B1.4 step 5). Front 6.2: "The rows-root download is not in the measured rows."
260 steps end at minute 172.2 against the 172.3 stop at 27.5 s/step, so any download time can cut steps (at
23.65 s there are about 16.8 min of slack). The HF guide says "merge data into fewer files ... Parquet or
WebDataset" [web: huggingface.co/docs/hub/storage-limits]. Fix: B1.4 step 5 uploads the root as a few tar shards
(about 200 MB each) plus `manifest.json`, and the box check untars them; or keep the folder form and time
`snapshot_download` of the pinned revision once on the preflight host (front 6.1), and add that time to item 1.

### 11. MINOR (still open from round 5, data-4). B1.9 still names CR1 as the Gate B1 interval

Files: B (B1.9). The rule list says "the A1 decision interval", the closing paragraph says "The decision interval
of the simulation is the CR1 interval of the item-weighted mean (B.6.2), the same as the gate's". Gate B1 (78
patients; far 75, near 52) uses the percentile interval (B.6.2). Fix: delete the CR1 sentence; keep "the A1
decision interval (percentile at >= 50 patients, `cluster_t` below)".

### 12. MINOR. B1.5 adds `kind` and `cell` columns "when the record has `view`", but eval-view records have neither

Files: B (B1.5). "Add columns `kind` and `cell` (string "ij") to prompt rows when the record has `view`." Eval-view
records also have `view` (keys `fallback, kept_side, offset, op, side`, no `kind`, no top-level `cell`)
[measured: data/eval_views/lidc-val-shifted/lidc/val.jsonl; code: build_eval_views.py `_record`], and `final_eval`
builds its rows with `dataset.prompts(view, ...)` [code: final_eval.py:71]. A literal `r["view"]["kind"]` raises
`KeyError` in A3 and in the B1 evals. Fix: "when the record has `view.kind` (rows roots only); eval views get no
`kind`/`cell` column". U59 already runs `final_eval` on a view fixture and would catch it; add the eval-view case to
the B1.5 test list.

## Counts

Blocking 0. Major 3 (items 1, 2, 3). Minor 9 (items 4-12).
