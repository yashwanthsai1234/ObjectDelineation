# Final review round 15: fresh outside architect

Reviewer role: a fresh architect with no prior opinion. Scope: the whole plan `docs/implementation2.0/Implementation.md`
(round 15), read through its sources in `phases/`. The focus is whether this is the best end-to-end architecture
up to Phase E, with GPU everywhere, and what a top lab would do differently. I read `docs/implementation2.0/review-log.md`
(round 13-14 entries of every phase) and the r14 fresh-arch review and verdict first. I do not re-raise r14
fresh-arch-1 (fixed), -2 (downgraded), or -3..-8 (minor, still open; item 1 below only points at -8 where its fix
touches the same line).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Loaded. No item depends on an NVIDIA product or SDK choice; the front already lists the catalog skills for the 3D step. No catalog lookup was needed and nothing was installed. |
| engineering:code-review | Read `experiments/phase2/train_rl.py` (whole file), `experiments/phase2/hf_pull.py`, `policy.py` (`CARDS` :51, constants :58-73, `accrued`/`spend` :190-199, `cmd_launch` refusal :321-329, `chain` :422-452, `cmd_run` :464-495, `cmd_push_hf`/`cmd_fetch_hf` :511-536, sweep :588-592), TRL `RepeatSampler` call (grpo_trainer.py:1372-1379), and huggingface_hub 2.2.0 `snapshot_download` / `hf_hub_download` / retry code. Found item 1. |
| engineering:debug | Traced the V1 time line from launch through `fetch_hf.sh --rows` to step 1 and to `BudgetStop`, with the HF request limits applied (item 1). |
| engineering:system-design | Re-ran the V1 minute arithmetic of front 6.2 with the download time added (item 1). Checked that a3_measure's separate plan is a deliberate failure boundary (it still runs when step 3b fails), so no item. |
| engineering:architecture | Wrote the overall judgement below as an ADR-style review: seams, layers, data-first order, (b') end state. |
| engineering:documentation | Every major item quotes the source file and heading, plus the code line or web source. |
| huggingface-skills:huggingface-best | No model choice is in question (Qwen3-VL-2B now; 4B gate on val). Used the HF Hub side: account tier limits for item 1. |
| system-design-skills:system-design | Ran the loop on "data onto the box": requirement (pinned, hashed rows root), estimate (16,640 objects, about 1.4 GB), constraint (free-tier resolver quota), failure mode (silent throttling eats training minutes). Found item 1. |
| engineering:tech-debt | Checked what stays on `main` under (b') after E and CLEAN: no new duplicate found after the r14 sampler fix. Item 2 is a small consistency debt. |

Web checks:
- Hugging Face Hub rate limits: "All quotas are calculated over 5-minute fixed windows"; Free user: API 1,000,
  Resolvers 5,000; PRO 12,000. Resolvers are "all the URLs that contain a `/resolve/` segment in their path"
  [web: https://huggingface.co/docs/hub/rate-limits, table "Rate limit Tiers (in September '25)"].
- `huggingface_hub` 1.2.0+ waits for the `RateLimit` reset on a 429 and then retries
  [web: same page, "Smart rate limit handling"]. The box lock pins `huggingface-hub==2.2.0`
  [code: experiments/phase2/docker/locks/train.txt:84], so the download waits, it does not fail.

## Overall judgement

The architecture is sound, and close to what a small top-lab team would build under a $3 cap:
- Data first answers the measured diagnosis ("trained centred, tested shifted"). The rows root is built and gated
  on the Mac (21 B data gates), pinned by HF revision and manifest SHA-256, and checked again on the box.
- One CUDA engine source for the CPU build and nvcc, with the CPU build as the oracle, and
  `tests/test_gpu_scoring.py` as step 0 of every paid plan. A bad GPU build stops in minutes.
- The three seams (`geom.polygonize`, `policy.sampler`, `rl.trl_adapter`) are at the right places. B1.4a puts the
  polygonizer seam in front of today's engine before the targets are built, so the owner rule "CUDA for the
  polygonizer" holds from V1 on, and D2 swaps K1 in behind it later.
- The (b') end state is honest: the GPU pixel path and reward v5 stay off `main` until a raise funds their gates.
- Failure boundaries in V1 are explicit: a soft `train_rl` step, a hard `--merge-last` step, a separate
  `a3_measure.yaml` plan so the R5 baseline is made even when B1 fails, a pull job for checkpoints, and a
  session-spend `BudgetStop` that leaves the eval tail paid for.

What a top lab would do differently (not defects unless an item says so):
1. Ship every dataset to the box as a few immutable shards (one object per root, or a handful), not as 16,640
   loose files. The repo already does this for the bundle (`sync_up` copies one tar and runs `tar -xf`
   after a check [code: experiments/phase2/vast/policy.py:386, 417]). The rows root is the exception, and it collides with the HF free-tier
   request quota (item 1).
2. Treat every way a training run can end short (crash, `BudgetStop`, cap) as one labelled state in the result
   file, and key the decision table on that label. Today only the crash path writes `partial: true` (r14
   fresh-arch-8, still open; item 1's fix adds the one field it needs).
3. Long term, sample windows on the device from the source slices (ADR-0015 option (b)) instead of materialising
   PNG views. At the 3D voxel step, materialising every scheduled sub-volume will not scale. ADR-0015 already names
   the trigger; no item.
4. Training-seed variance is not measured, so a Gate B1 pass is a single-run result. Under $3 that cannot be
   bought. The plan says so (B.6.2 "Noise floor", DOC11). No item.

## Items

### 1. MAJOR. The V1 rows-root download is about 16,640 HF file requests on a free account (5,000 per 5 min). The 10-15 min wait is in no item's minutes, and at the base step time it makes `BudgetStop` cut the 260-step B1 run at about step 228-239, with no label

Sources: B (B1.4 step 5; B1.5 "Rows-root download"), 00_front (6.2 item 1, notes "Why 260 steps"; 6.3 V2 item 1
and V3 item 0; 4 "Data versioning" row).

Problem.
- B1.4 step 5: "the root holds about 16,640 files (one image and one cache per window; today's 7,341-window train
  view holds 14,684 files [measured ...])". I re-counted today's view: 7,341 `png` + 7,341 `npz` + 2 = 14,684,
  1.2 GB [measured: `find data/views/lidc-train-shifted -type f`, 2026-11-09].
- B1.5 "Rows-root download": "It runs `hf_pull --repo svenk029/objdel-sft --revision REVISION --out
  /workspace/data/lidc --allow "rows/ROOT/**"` ... So the call is `snapshot_download(local_dir=..., allow_patterns=
  "rows/lidc-train-b1/**", revision=...)`".
- Code: `snapshot_download` lists the tree, then calls `hf_hub_download` once per file in a thread map of 8 workers
  [code: .venv huggingface_hub/_snapshot_download.py:519-545]. Each `hf_hub_download` first asks for the file
  metadata at the file's `/resolve/` URL [code: file_download.py `get_hf_file_metadata`, retry on (408, 429) at
  :398]. The xet read token is fetched once per revision and cached, so it does not add per-file calls [code:
  utils/_xet.py:145-194]. So one rows-root fetch makes at least 16,640 resolver requests.
- HF limit: a free user has 5,000 resolver requests per 5-minute fixed window [web: huggingface.co/docs/hub/
  rate-limits]. The repo belongs to a free account: "the owner's free HF account has no card; storage only"
  [code: experiments/phase2/hf_push.py:8-9]. 16,640 / 5,000 = 3.3 windows, so the call waits for at least three window
  resets: about 10-15 min plus the transfer [estimate: (4 - 1) windows x 5 min, minus the unknown position inside
  the first window]. huggingface_hub 2.2.0 waits on each 429 and retries, so the call does not fail; it is slow.
- 00_front 6.2 item 1: "launch to first `run` up to 38.2 min [measured] + build 3 + `env_check` 1 + smoke 5 +
  start-up 6 [estimate]" = 53 min. The notes say: "The rows-root download is not in the measured rows. At 27.5 s and
  a 53-min start, 260 steps end at minute 53 + 119.2 = 172.2, before the stop." The stop is at minute 172.3. So
  the plan has 0.1 min of slack at the base step time, and the download is the one known unpriced step before
  step 1.
- With the download: step 1 starts at minute 63-68. `BudgetStop` at $2.05 still trips at minute 172.3, so the run
  ends at step (172.3 - 63..68) x 60 / 27.5 = 228-239 instead of 260 [estimate: section 6 inputs]. At the r4
  start (39 min + 10-15) it ends at about step 258-260. At 23.65 s it still fits (53 + 15 + 102.5 = 170.5).
- A `BudgetStop` end is a normal end: `train_rl` writes `final/` after `trainer.train` returns [code:
  experiments/phase2/train_rl.py:55-60], and step 3b "does nothing and exits 0" when `final/` exists (B1.5). No
  `partial: true` is written, so 00_front 1.1 row "B1 stopped early" does not fire, and the Gate B1 report labels a
  ~230-step model as the planned 260-step run. The cut is also mid-block (228-239 is not a multiple of 20), while
  B1.2/B1.4 gate coverage (DG7-DG9) only "on every block prefix".
- The same fetch is in V2 item 1 and V3 item 0 ("the rows-root download is not in the measured rows (as in V1)"),
  and V3 adds two more `fetch_hf.sh` calls in the same window.
- No planned check measures it: 6.1 step 4 lists the HF tree on the Mac (API bucket, about 17 paged calls), but no
  step downloads the root. So it is first seen on the paid box. This passes every $0 check, so the caught-before-
  spend rule does not apply.

Result. An engineer who follows the doc exactly, at the plan's own base step time and the worst measured start, gets
a B1 run 8-12 % shorter than the owner approved (260 steps), cut inside a block, with no label in `rl.json` and no
matching row in the decision table. Money is not lost (the eval tail is still paid for), but the one paid
experiment does not test what the approved table says it tests, and about $0.11-0.17 of box time goes to waiting.

Fix (B, 00_front):
1. B1.4 step 5: after the gates pass, pack the rows root into ONE uncompressed tar
   (`rows/lidc-train-b1.tar`; PNG and npz are already compressed, about 1.4 GB) with sorted member names, mode 0644,
   uid/gid 0 and mtime 0, so the tar bytes are reproducible. Upload the tar and `manifest.json` (as
   `rows/lidc-train-b1.manifest.json`) with one `upload_folder` call. Record the tar's SHA-256 in the manifest
   (`schedule.tar_sha256`) and in `b1_rows.yaml` beside `--rows-manifest-sha256`. This matches how the bundle
   already reaches the box (one tar, `tar -xf` after a check).
2. B1.5 "Rows-root download": `hf_pull --allow "rows/lidc-train-b1.*"` (2 resolver requests). The same ssh call
   checks the tar SHA-256, then extracts with Python `tarfile` and `filter="data"` into `{data}/rows/lidc-train-b1/`
   (the box treats the tar as data), then runs the existing manifest SHA-256, `files_count` and `files_sha256`
   checks. 6.1 step 4 compares the HF LFS `sha256` of the tar at `rows_revision` (from `get_paths_info`) with the
   local tar, instead of the 16,640-path list. Tests: the I-test of the rows fetch plants a tar member with `..` and
   expects a refusal; a second build of the tar from the same root has the same SHA-256.
3. If the owner prefers loose files: add a 6.1 step 4e ($0, on the Mac with the same account): time a full
   `hf_pull --allow "rows/lidc-train-b1/**"` into a temp folder, add the measured minutes to 6.2 item 1 (and V2
   item 1, V3 item 0), and recompute `--steps` as the largest multiple of 20 that ends before the stop at 27.5 s.
4. Either way, close the label gap on the same line (r14 fresh-arch-8): `BudgetStop` writes
   `budget_stop_step` to `rl.json` when it trips, and 00_front 1.1 row "B1 stopped early" also matches
   `budget_stop_step < 260`.

### 2. MINOR. 00_front 6.2 item 1 says its 53 min include the rows-root download, while its notes say the download is not in the measured rows

Sources: 00_front (6.2 item 1 row; notes "Why 260 steps"), 00_front (6.3 V2 item 1, V3 item 0).

The item 1 row text is "Before step 1: launch, sync, checkpoint and rows-root download from HF (pinned revisions),
CUDA build, ... | 53 |". Its basis column and the notes say "The rows-root download is not in the measured rows".
A reader of the table alone believes the download is priced. This does not change money by itself (item 1 above
does). Fix: in the item 1 row, write "rows-root download (not priced; see notes)" until item 1's fix lands, then
"rows-root tar download (2 requests, about 1 min [estimate])" with the minute added to the 53.

### 3. MINOR. The V1 run plan has 0.1 min of slack between the last planned step and `BudgetStop` at the base step time

Sources: 00_front (6.2 notes "Why 260 steps").

"At 27.5 s and a 53-min start, 260 steps end at minute 53 + 119.2 = 172.2, before the stop" (172.3). Any minute
not in the estimate (item 1's download, a slower checkpoint save, a slower host) cuts the run. The plan already
accepts this ("If it trips, the adapter is saved and items 3-5 still run"), so this is a choice, not a defect.
A top-lab plan would keep one block (20 steps, about 9 min) of slack, or state the expected step count as a range
in the owner message of 6.1 step 0 (for example "240-260 steps at 27.5 s, 260 at 23.65 s"). Fix: add that range
to the 6.1 step 0 owner message and to 99_back decision 2.

## Counts

blocking 0, major 1 (item 1), minor 2 (items 2-3).
