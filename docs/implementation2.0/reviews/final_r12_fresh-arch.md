# Final review round 12: fresh outside architect

Reviewer role: a fresh architect with no prior opinion. Scope: the whole plan `docs/implementation2.0/Implementation.md`
(round 12), read by its sources in `phases/`. I focused on the one paid path that the default (b') funds (V1, front
6.2), on the GPU-everywhere data path, and on the layer design up to Phase E. I opened the code each item cites.
I did not re-raise items that `docs/implementation2.0/review-log.md` marks as handled (r10 fresh-arch-5 "branch work
after E", r10 fresh-arch-6 "rental tables from YAML", r11 fresh-arch-1..7 were checked before writing).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Loaded. No item depends on an NVIDIA product or SDK beyond the repo's own CUDA engine, so I did no catalog lookup and recommend nothing. The front's catalog recommendations (sanitizer skill, `nv-segment-ct`, `dicom-series-to-volume`) are unchanged. |
| engineering:code-review | Read the code behind each V1 claim: `policy.py:372-383` (`code_files()` drops `^(data|docs)/`), `bundle.py:26-75` (bundle members: `prepared/`, `views/`, `eval_views/`, `checker_eval/`; no `docs/`), `policy.py:480-500` (`cmd_run`: `rm -f {rd}/EXIT` at the start of each plan, heartbeat killed at the end), `policy.py:500-536` (`cmd_pull`, `cmd_push_hf`, `cmd_fetch_hf`), `watchdog.sh:8, 36-40`, `com.objdel.vast.pull.plist`, `train_rl.py:1-72`, `final_eval.py:1-110`, `certify.py`, `gpucheck_py.cpp:140-158`. This found item 1. |
| engineering:debug | Traced the V1 box from `sync_up` to the last eval: which files exist on the box, which step reads them, which step refuses. Item 1 is the trace "B1.6 frozen-view check -> E05 refusal on a missing file -> `code_files()` drop regex". Item 2 is the trace "pull job -> PULLED_OK after the idle start -> watchdog destroy". |
| engineering:system-design | Checked the V1 time budget as a capacity plan: item 1's 53 min has no term for the rows-root download, and the 260-step plan has 0.1 min of slack at the base step time (item 3). Checked the eval cost model against `final_eval.json`: scoring is 0.7 s of 415 s, so the pool-size-1 rule costs nothing. |
| engineering:architecture | Judged the end-to-end design (section "Overall judgement" below) as an ADR reviewer: three seams, eight layers, (b') end state, `after-raise` branch. Items 5 and 6 are design trade-offs, so they stay minor. |
| engineering:documentation | Every blocking and major item quotes the source file and heading and the code line or measurement. |
| huggingface-skills:huggingface-best | No model choice is in question; the 2B-first decision and the 4B gate stand. The HF part of this review is the Hub transfer path (item 3): I read the Hub rate-limit page and `huggingface_hub` 2.2.0 (`train.txt:84`, `vllm_native.txt:236`) `snapshot_download` code. |
| system-design-skills:system-design | Ran the loop on V1 only: requirements (Gate B1 needs items 3-5), estimate (minutes per item), failure modes (what destroys the box, what refuses an eval). The failure-mode pass found items 1 and 2. |
| engineering:tech-debt | Scored the long-lived `after-raise` branch and the size of the plan as debt; both were raised before (r10), so only item 6 notes the residual. |

Web checks: Hugging Face Hub rate limits (https://huggingface.co/docs/hub/rate-limits, table "current rate limits (in
September '25)": free user 5,000 resolver requests per 5-minute fixed window; `huggingface_hub` >= 1.2.0 waits for the
window reset on a 429).

## Overall judgement

The design is sound where it matters most. Data first (B1 with coverage gates before any spend) is the right call
for a diagnosis that says "trained centred, tested shifted". One CUDA engine source for CPU and nvcc, with the CPU
build as the oracle, is the right way to keep the checker exact on a $0 Mac. The seams (`geom.polygonize`,
`policy.sampler`, `rl.trl_adapter`) are at the right places. A top lab would do four things differently, and none is
a defect of this plan: (a) ship training data as a few large shards, not 16,640 small files (item 3); (b) run one
plan file per rental so no idle gap exists between paid steps (item 2); (c) use one generation path for training and
evaluation instead of a second `VllmSampler` that E4 must prove equal by TOST (item 5); (d) budget at least one
seed replicate before calling a +0.01-0.02 gain a data effect (the plan states this limitation, 99_back 8b). The
one serious gap is item 1: as written, every V1 eval refuses on the box, so the only funded rental produces no
Gate B1 result.

## Items

### 1. BLOCKING. Every V1 eval refuses on the box: `docs/reference/frozen_views.json` never reaches it

Sources: B (B1.6 "Freeze the test view; build the guard views"; B.9.4 E05), A (A3 "Steps" preamble), 00_front (6.2
"Plan files", 6.1 step 4c). Also CLEAN (CLEAN1 step 1) and D (section 5 / D1 step 0 keep rule), which keep `docs/`
out of the code sync.

Problem.
- B.md, B1.6: "`experiments/phase2/final_eval.py`: for a frozen view, load `frozen_views.json`, recompute the content
  hash and the target-set hash, and stop with a clear message on a mismatch." The file is "`docs/reference/frozen_views.json`,
  one entry per frozen view (`lidc-test-shifted`, `lidc-val-edge`, `lidc-val-shifted`)".
- B.md, B.9.4 E05: "`final_eval` refuses a changed test view ... **missing `frozen_views.json` gives a clear refusal**".
- A.md, A3: "The stamp's data hash covers the val views in `preflight.EXAMS`, not the test or edge view
  [code: preflight.py:28]. Those two are checked by `final_eval` against `frozen_views.json` (B1.6)."
- 00_front 6.2: item 3 runs `final_eval` on `lidc-val-shifted` and `lidc-val-edge` (B1 final); items 4-5 run
  `final_eval` on `lidc-val-edge`, then `lidc-val-shifted` (R5). All three are on the box.
- Code: the code sync drops every `docs/` path: `drop = re.compile(r"^(data|docs)/|^experiments/[^/]+/(lab|results)/|..."`
  [code: experiments/phase2/vast/policy.py:379]. CLEAN1 step 1 keeps the drop (`^(data|docs|archive|\.claude)/`),
  and D1 step 0 adds a keep rule only for `docs/implementation2.0/kernels/K1|K2|K3`. The bundle holds only `prepared/`,
  `views/`, `eval_views/` and `checker_eval/` [code: experiments/phase2/vast/bundle.py:49-75]. So `/workspace/src/docs/`
  does not exist on the box.

Failure path. Item 3's process reaches `lidc-val-edge` (a frozen view) and refuses (E05), so no B1 eval is written.
Item 4 is the first view of the R5 process, so items 4-5 refuse too. Gate B1 then has no box input: no B1 val
numbers, no A3-val baseline, no edge baseline. Under (b') no V2 item 1b follows, so front 1.1 ends in row 4 ("edge
baseline missing", stop) after $2.40 of spend. The same refusal hits V2 item 1b, V3 (`a3_test.yaml`, B3.9 step 3, the
B3 eval, the one test eval) and V4.

No planned $0 check catches it. 6.1 step 4c compares the bundle's view hashes with `frozen_views.json` on the Mac,
where `docs/` exists. CLEAN1 step 5 (`test_bundle_imports_standalone`) copies `code_files()` but only imports the
modules; the file is read inside `main()`. E05 and I05 run in the Mac repo. I47 parses argv only.

Fix (B1.6, with one line each in A3, front 6.2 "Plan files" and 6.4 rule 7):
1. Ship the record with the code: in `code_files()` add one keep line beside D1 step 0's kernel rule, so
   `docs/reference/frozen_views.json` ships (CLEAN1 keeps its drop regex; the keep rule runs after it). Or move the
   record to a shipped path (`experiments/phase2/configs/frozen_views.json`) and update every reference (B1.0 step 3,
   B1.3 `rows_leakage`, B1.6, A3, DG15, I05, E05, DOC "reference" list).
2. Add a $0 test in `tests/test_vast.py` next to CLEAN1's standalone test: `docs/reference/frozen_views.json` (or the
   new path) is in `P.code_files()`, and in the `code_files()` copy of `test_bundle_imports_standalone`, the
   function that `final_eval` uses to look up a frozen entry finds `lidc-val-edge` (child `python -I`, `cwd=tmp_path`).
3. Name the test in front 6.1 step 4 and in order rule 2's list of checks that pass on the `v1-b1` tree.

### 2. MAJOR. The new pull job lets the watchdog destroy the V1 box between `b1_rows` and `a3_measure`

Sources: 00_front (6.2 "Launch lines" and "Why the pull job (final round 11, fresh-arch-2)"), B (B1.5 "V1 launch
lines").

Problem.
- 00_front 6.2: "`run.sh <run> experiments/phase2/configs/b1_rows.yaml       # items 2-3; wait for EXIT`" then
  "`run.sh <run> experiments/phase2/configs/a3_measure.yaml    # items 4-5; start right after the b1_rows EXIT`".
- 00_front 6.2 "Why the pull job": "One side effect: the watchdog destroys an idle box 30 min after a pulled `EXIT`
  ... so `a3_measure.yaml` must start within minutes of the `b1_rows` `EXIT`, **as today**."
- Code: the watchdog destroys when the GPU is idle for `IDLE_S` (default 1,800 s), `EXIT` exists, and `PULLED_OK` is
  newer than the idle start [code: experiments/phase2/vast/watchdog.sh:8, 36-40; policy.py:603 `--idle-minutes`
  default 30]. `cmd_run` writes `EXIT` at the end of a plan and removes it only when the next plan starts [code:
  policy.py:492]. The pull job runs `pull.sh` every 1,200 s [code: com.objdel.vast.pull.plist:11], and `cmd_pull` writes
  `PULLED_OK` on the box after every pull [code: policy.py:503-508].

"As today" is not true. Today `PULLED_OK` appears only when the main agent itself pulls, so an agent that checks for
`EXIT` with `ssh test -f` and is slow only pays for idle time; `a3_measure` still runs later. With the job loaded,
a pull after the idle start comes within 20 min without the agent, so the box is destroyed 30 min after `EXIT`
whatever the agent does. The plan never says how the agent waits for `EXIT` (poll interval, tool), and the `EXIT`
time is not known in advance (260 steps, `BudgetStop`, or a crash plus step 3b). If the agent starts `a3_measure`
more than 30 min late, items 4-5 (not skippable; the paired baseline and the only edge baseline) are lost, and under
(b') Gate B1 takes row 4 of 1.1. No $0 check covers the gap.

Fix (00_front 6.2 and B1.5; A3 keeps `a3_measure.yaml` for V2 item 1b):
- Preferred: remove the gap. Append the two `a3_measure.yaml` steps to `b1_rows.yaml` as steps 5-6, and make step 3b
  and step 4 soft (`~`), so a failed B1 eval or merge does not skip the R5 process (`&&` join, policy.py:490). V1
  then has one `run.sh` call for items 2-5. 6.1 step 4d and I47 check that steps 5-6 equal `a3_measure.yaml`'s steps
  byte for byte. Re-state "Why items 3-5 always fit" for one plan (the tail is unchanged: 26 min).
- Or, if the two plans stay: pass `--idle-minutes 120` on the V1 launch line (the worst case is set by the cap,
  6.2 "Worst case", so it does not change), and write the wait rule: the main agent polls `test -f {rd}/EXIT` over
  ssh every 2 min from step 260's expected time and starts `a3_measure` within 5 min of it.
- In both cases delete "as today".

### 3. MAJOR. The rows-root download is unpriced, and the 260-step plan has 0.1 min of slack at the base step time

Sources: 00_front (6.2 item 1 and "Why 260 steps"; 6.3 V2 item 1 and V3 item 0), B (B1.4 step 5, B1.5 "Rows-root
download").

Problem.
- 00_front 6.2 item 1: "53 | launch to first `run` up to 38.2 min [measured] + build 3 + `env_check` 1 + smoke 5 +
  start-up 6 [estimate]". "Why 260 steps": "The rows-root download is not in the measured rows. At 27.5 s and a
  53-min start, 260 steps end at minute 53 + 119.2 = 172.2, before the stop" (`BudgetStop` at minute 172.3).
  So the plan has 0.1 min of slack at the base step time, and the download time is set to zero.
- B1.4 step 5: the root "holds about 16,640 files (one image and one cache per window)", 1.2-1.6 GB.
- Code: `fetch_hf --rows` is `snapshot_download(... allow_patterns="rows/lidc-train-b1/**")` (B1.5), which calls
  `hf_hub_download` once per file with 8 workers [code: .venv huggingface_hub/_snapshot_download.py:519-546]. For a
  non-Xet file, each call makes a `HEAD /resolve` request first [code: .venv huggingface_hub/file_download.py:1363-1366,
  "otherwise => HEAD /resolve call"]. The box image pins `huggingface-hub==2.2.0` [code:
  experiments/phase2/docker/locks/train.txt:84].
- Web: a free account gets 5,000 resolver requests per 5-minute fixed window [web: huggingface.co/docs/hub/rate-limits];
  `huggingface_hub` >= 1.2.0 sleeps to the window reset on a 429.

The download time is not known. If the files are Xet-stored, per-file latency at 8 workers still costs minutes. If
they are not, 16,640 resolver requests need at least 4 windows: more than 15 min. Each minute at 27.5 s is 2.2
steps. At 15 min the stop cuts the run at about (172.3 - 68) x 60 / 27.5 = 227 steps, not 260, and the cut is
mid-block (B1.2 trains whole 20-step blocks; `BudgetStop` does not). The same unpriced minutes sit in V2 item 1 and
V3 item 0. A $0 measurement exists and is not planned: download the uploaded root to the Mac with the same
`hf_pull --allow` call.

Fix (B1.4 step 5, B1.5, 00_front 6.2):
1. Ship the rows root as one file: B1.4 step 5 writes `lidc-train-b1.tar` (uncompressed; PNGs are already
   compressed) plus its SHA-256 and uploads it with `HfApi.upload_file` to `rows/lidc-train-b1.tar`. `fetch_hf --rows`
   downloads that one file (1-2 resolver requests), checks the tar SHA-256, untars to `{data}/rows/lidc-train-b1/`,
   then runs today's manifest SHA-256 and file-count check. The preflight "ROWS revision" check compares the tar's
   SHA-256 at `rows_revision` instead of a 16,640-path listing. This is not a stateful system, so section 4's
   "Not adopted" row is unaffected.
2. Before V1 ($0, Mac): time `fetch_hf --rows` against the pinned revision into a scratch folder, record the minutes
   as `[measured]` in item 1, and recompute "Why 260 steps". If slack at 27.5 s is below 5 min, plan 240 steps (a
   multiple of 20) or raise the amount through owner decision 2.

### 4. MINOR. `push_hf.sh <run> b1_launch.yaml` has no stated way to select `hf_push --rl`

Sources: 00_front (6.2 "The launch lines need five code changes"), B (B1.5 "B1 adapters to V3").

B1.5: "B1.5 adds `hf_push --rl <train_rl out>` ... `push_hf.sh <run> b1_launch.yaml` runs it on the box". But
`cmd_push_hf` builds `hf_push --run .../{output_dir} --repo ... --revision {output_dir}` with no mode flag [code:
policy.py:511-524], and I41 tests `hf_push --rl` directly, not the command `cmd_push_hf` sends. If the engineer adds
only `hf_push --rl`, the box call takes the SFT path and fails, and the `&&` in "`pull.sh <run> && push_hf.sh ... &&
teardown.sh`" skips teardown (the box then idles until the watchdog or the cap). The model is still on the Mac via
`pull.sh`, so nothing is lost. Fix: `cmd_push_hf` passes `--rl` when the config has `rows_plan`; add a `DRY_RUN=1` test
in `tests/test_vast.py` that asserts the sent command for `b1_launch.yaml` contains `--rl`; put `teardown.sh` on its
own line after the push (as `sweep_orphans.sh` is) so a push failure never skips it.

### 5. MINOR. Two generation paths for one model (top-lab note)

Sources: 00_front (2.2 seam table: "`policy.sampler` ... Evaluation only; training generation stays in TRL's
`VLLMGeneration`"), E (E4).

Training samples through TRL's `VLLMGeneration`; evaluation, after E4, through a new `VllmSampler`. The plan then
needs V4 (a budget raise) and a TOST to show the two agree (E-G10), and until V4 passes E4 keeps the old samplers.
A top lab keeps one sampler: one thin wrapper over the vLLM engine that both TRL's generation hook and `final_eval`
call, so prompt ids and `SamplingParams` are equal by construction and E-G10's box half is a unit test. Not a defect:
the plan's TOST is a valid guard. Suggested only for ADR 0020.

### 6. MINOR. Residual debt of the (b') shape (top-lab note)

Sources: 00_front (1 "End state under (b')"), E (E.1, E.4 R-4).

Under (b') the plan builds C1-C11 and B2.2-B3 code on `after-raise` (pre-E paths), then E moves every file and
replaces the YAML plans, so a later raise needs the E R-4 port (more than a path map after E6). r10 fresh-arch-5
already raised the order; I add one point: Phase E restructures a 2D pipeline right before the stated 3D voxel
direction (front section 4, last rows). A top lab would land E0's guardrails (import-linter contract, goldens) now
and defer the E5/E6 splits until the first 3D interface is known (ADR 0020 "3D extension"), so the layers are cut
once. Owner choice; no outcome changes.

## Counts

blocking 1 (item 1), major 2 (items 2, 3), minor 3 (items 4, 5, 6).
