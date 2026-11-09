# Final review round 13: fresh outside architect

Reviewer role: a fresh architect with no prior opinion. Scope: the whole plan `docs/implementation2.0/Implementation.md`
(round 13), read by sections, with the sources in `phases/` as the files to fix. Focus: is this the best end-to-end
architecture up to Phase E with GPU everywhere, and what would a top lab do differently. I opened the code each item
cites. I checked `docs/implementation2.0/review-log.md` (front r6-r12, A r9-r12, B r9-r12, E r1-r3) and the r12
fresh-arch review and verdict before writing, and I do not re-raise r12 fresh-arch-2 (pull job gap), -3 (rows-root
download), -5 (two generation paths) or -6 (E before 3D).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Loaded. No item depends on an NVIDIA SDK beyond the repo's own CUDA engine, so no catalog lookup and no install recommendation. I did one web check of the plan's claim "no GPU PNG decoder" (section 2.1): NVIDIA lists only nvJPEG, nvJPEG2000 and nvTIFF as GPU codecs, PNG falls back to the CPU, so the claim stands and is not an item. |
| engineering:code-review | Read `policy.py` (`code_files` :374-383, `cmd_run` :465-495, `cmd_pull`, `cmd_push_hf`, `cmd_fetch_hf`), `bundle.py` :26-75, `train_rl.py`, `final_eval.py`, `serve.py` pool code, `tests/test_p2_breakit_ops.py:297-300`, `tests/test_gpu_scoring.py`, `.gitignore`. This found item 1. |
| engineering:debug | Traced the CLEAN5 move from `mv` to the owner's `git add -A` against the real `.gitignore` with `git check-ignore --no-index` (item 1), and traced the default-budget path from V1 to "the one test read" (item 2). |
| engineering:system-design | Checked the V1 money and time arithmetic for the item 2 fix (tail, cap, worst case, `launch_refusal`, reserve) at r_s = $0.01133/min. |
| engineering:architecture | Wrote the overall judgement below as an ADR review: three seams, eight layers, (b') end state, data-first order. Design trade-offs stay minor. |
| engineering:documentation | Every blocking and major item quotes the source file and heading and the code line or measurement. Item 4 is a stale header. |
| huggingface-skills:huggingface-best | No model choice is in question (2B first, 4B gate on val). The HF part of this review is the test-read path of item 2 (answers saved on the box, scored on the Mac). |
| system-design-skills:system-design | Ran the loop on "the default end state": requirement (one test number for the final candidate), constraint ($3, (b')), failure mode (no path to the number). Found item 2. |
| engineering:tech-debt | Scored the CLEAN wave-1 moves as debt removal with a hidden liability (ignored files that change status when moved). Found item 1; items 3-5 are small debts. |

Web checks: nvImageCodec codec support (https://developer.nvidia.com/nvimagecodec: GPU codecs nvJPEG, nvJPEG2000,
nvTIFF, CPU fallbacks; https://arikpoz.github.io/posts/2025-04-07-fast-image-loading-with-nvidia-nvimagecodec: PNG
falls back to CPU decoding).

## Overall judgement

The architecture is sound and close to what a small top-lab team would build under a $3 cap. Data first (B1 with 21
data gates on the Mac before any spend) is the right answer to the measured diagnosis "trained centred, tested
shifted". One CUDA engine source for CPU and nvcc, with the CPU build as the oracle and `tests/test_gpu_scoring.py` as
step 0 of every paid plan, keeps the checker, scorer and polygonizer exact on a $0 Mac and stops a bad GPU build in
the first minutes of a rental. The three seams (`geom.polygonize`, `policy.sampler`, `rl.trl_adapter`) sit at the
right places, and the layer order core < geom < verify < data < policy < {rl | sft} < infra matches the import
graph the plan measured.

What a top lab would do differently (none of these is a defect except where an item says so):
1. Under a fixed budget, produce the goal metric in the one funded rental: generate the final candidate's test
   answers inside V1, keep them unscored, and score them on the Mac only after the val-only decision closes
   (item 2). The plan today defers the test read to a V3 that (b') never funds.
2. Treat every move of a results tree as a licence event: re-check ignore status after the move (item 1).
3. Ship data as a few large shards and run one plan per rental (r12 fresh-arch-2/-3, already logged).
4. Keep one generation path for training and evaluation (r12 fresh-arch-5, already logged).
5. Budget a seed replicate before reading a +0.01-0.02 gain as a data effect (stated limitation, 99_back).
6. Keep the GPU pixel path decision visible: under (b') `main` keeps "today's CPU pixel path only" (00_front "End
   state under (b')"), so the owner rule "GPU everywhere" is met for checker, scorer and polygonizer but not for the
   image path. The plan tells the owner this in decision 1; that is the right handling.

## Items

### 1. BLOCKING. CLEAN5's moves and `git add -A` commit RIGA per-case results that are git-ignored only at their old paths

Sources: CLEAN (CLEAN5 "Files to move", step 8; CLEAN.9 CT-87), with 00_front 6.4 rule 7 ("RIGA never leaves") and
DOC8 (records under `docs/research/records/`).

Problem.
- CLEAN.md, CLEAN5 "Files to move (plain `mv`)": "`experiments/phase0/*.md` and `experiments/phase0/results/` ->
  `docs/research/records/phase0/`" and "`experiments/phase1/lab/` (all remaining files) -> `archive/phase1_lab/`".
- CLEAN.md, CLEAN5 step 8: "Owner runs ...: `git add -A experiments archive docs/research/records tests
  pyproject.toml src` ... then commits."
- 00_front 6.4 rule 7: "RIGA never leaves." CLEAN.2 row "RIGA terms": "copy, redistribution ... prohibited";
  "RIGA-derived tables already left the public MIT repo. So no tracked golden holds RIGA masks or per-case RIGA
  results (CLEAN4)."
- Code: `.gitignore` keeps these files out of git by their OLD paths only, under the comment "RIGA-derived results
  stay local: no permission to share RIGA data (repo public, owner 2026-11-06)":
  `/experiments/phase0/results/medsam_riga.csv`, `/experiments/phase0/results/riga/`,
  `/experiments/phase0/results/riga_outlines/`, `/experiments/phase0/results/vlm_box/`,
  `/experiments/phase0/results/vlm_polygon/`, `/experiments/phase1/lab/t_reward/results/h5_items.csv`,
  `/experiments/phase1/lab/t_data/riga/out_maxrg_check/`, `/experiments/phase1/lab/t_data/riga/maxrg_vs_prime.json`,
  `/experiments/phase1/lab/t_data/splits/out/riga_images.csv`,
  `/experiments/phase1/lab/t_reward/results/advantage_groups.json`,
  `/experiments/phase1/lab/t_reward/results/iter3_r5_raw.json`, `/experiments/phase1/lab/t_stack/results/t[24]_*.json`
  (the commit that added them: `e9f23d7 MIT licence; RIGA-derived result tables leave the public repo`).
- Measured (2026-11-09, read-only): the files exist and are untracked-ignored. `git status --porcelain --ignored`
  lists `!! experiments/phase0/results/medsam_riga.csv`, `!! .../results/riga/`, `!! .../riga_outlines/`,
  `!! .../vlm_box/`, `!! .../vlm_polygon/`, and the six lab paths above; `riga/masks.csv` holds per-image, per-expert
  RIGA rows (`BinRushed/BinRushed1-Corrected/image1,1,0.3454...`), `medsam_riga.csv` and `vlm_box/scores.csv` hold
  per-image RIGA IoU rows (`BinRushed_BinRushed2,43,disc,...`). After the move the new paths are NOT ignored:
  `git check-ignore -v --no-index docs/research/records/phase0/results/medsam_riga.csv
  docs/research/records/phase0/results/riga/x.png archive/phase1_lab/t_reward/results/h5_items.csv` returns exit 1
  (no rule matches).

Failure path. The engineer follows CLEAN5 exactly: `mv` carries the ignored files along (a plain `mv` moves the whole
folder), and the owner's `git add -A ... archive docs/research/records` stages them because no ignore rule matches
the new paths. The CLEAN5 commit then holds per-case RIGA results in a public MIT repo, which breaks the RIGA terms
and the owner decision of 2026-11-06. No planned $0 check catches it: the only RIGA guard,
`test_code_snapshot_has_no_riga_derived_results` [code: tests/test_p2_breakit_ops.py:297-300], reads
`P.code_files()`, which drops `docs/` and (after CLEAN1) `archive/`; CLEAN-G16 reads only `tests/goldens/*.npz`;
CT-87 counts the moved `results/` files and so passes when they are all there; the CLEAN.9 `git ls-files | grep -iE
'riga.*\.(png|tif|npy|npz)'` check misses `.csv` and `.json`, and `masks.csv`, `images.csv`, `scores.csv` have no
"riga" in the file name.

Fix (CLEAN5, plus one line in DOC8 and front 6.4 rule 7):
1. New CLEAN5 step 0, before any move: list every ignored file under the trees that move:
   `git ls-files --others --ignored --exclude-standard -- experiments/phase0 experiments/phase1/lab
   experiments/phase2/lab` (read-only). Move every hit that is not a `*.log` or `__pycache__` file to a git-ignored
   local folder outside the moved trees, `data/records/riga_local/<old path>` (`/data/` is ignored, `.gitignore:13`),
   and write that list into `archive/phase1_lab/README.md` and `docs/research/records/phase0/README.md` as "kept
   local, not in git (RIGA terms)". Only then run the moves. Delete the old `.gitignore` lines in the same commit,
   because their paths no longer exist.
2. Step 8: replace `git add -A ... archive docs/research/records ...` with `git add -A experiments tests pyproject.toml
   src` plus `git add archive docs/research/records` only after `git status --porcelain --ignored -- archive
   docs/research/records` shows no `??` path from the step 0 list.
3. New $0 test in `tests/test_p2_breakit_ops.py` next to the snapshot test: `test_no_riga_derived_file_is_trackable`.
   For every path in `git ls-files -co --exclude-standard` (tracked plus untracked-not-ignored, the set `git add -A`
   stages), fail when the path has a `riga` part in any folder or file name and ends in `.csv`, `.json` or `.jsonl`,
   or when it is one of the step 0 paths under its new prefix. Planted case: a temp repo with
   `docs/research/records/phase0/results/riga/masks.csv` fails; the same file under `data/` passes. Add it to CLEAN-G5
   and to CLEAN5 step 9's gates, and run it before the owner's `git add`.

### 2. MAJOR. Under the default (b'), the final candidate never gets its one test read; the plan has no step for it

Sources: 00_front (1.1 decision table header and rows; 6.3 V3 item 5 and the (b') rule; "End state under (b')"), 99_back
(decisions 1, 2b, 17), A (A3 `--rescore`), B (B1.5 `b1_rows.yaml` step 4).

Problem.
- 00_front 1.1 header: "Final candidate (the one model that reads test)"; row "Pass": "Best of B1 and B3 on
  val-select ... A tie goes to B1". Time line step 6: "The test eval runs once, after the owner closes the 4B gate
  (section 3)."
- 00_front 6.3 V3 table, item 5: "The one test eval of the final candidate, only if the owner has closed the 4B gate".
  It is a V3 row. 00_front 6.3: "Rule: under the default (b') (owner decision 1) no rental follows V1. V2, V3 and V4
  each wait for a budget raise". 99_back decision 1 marks (b') as the default.
- 00_front "End state under (b')" lists what lands on `main` and what waits on `after-raise`, and says nothing about
  the test read. No plan file or task samples the test view for a new candidate: the only plan file that names
  `lidc-test-shifted` is `a3_test.yaml` (R5, V3 item 3; A.md:251 "No V1 plan file reads the test split"), and the
  candidate's test eval exists only as the V3 item 5 row [measured: grep of `phases/00_front.md`, `A.md`, `B.md` for
  `lidc-test-shifted`]; B I28 pins that only `final_eval` scores the test view.
- 6.1 step 4c: the V1 bundle holds "no `lidc-test-*` folder".

Result. On the default path, Gate B1 passes (or fails only the far-miss rule), B3 never runs, so B1 is the final
candidate, and the plan's goal number (B1 on test, to compare with R5's 0.645) cannot be made without a rental that
(b') forbids. The engineer must guess: rent again (breaks the budget rule), read test from something else, or report
no test number. Each guess changes the outcome, and the owner approves (b') in decision 1 without being told that it
means "no test number for any new model".

A cheap fix fits the $3 and keeps "test is read once, after the 4B gate". Generating answers is not reading them; A3
already has the parts (answers written before scoring, `--rescore DIR` later [A.md A3 "write order", `--rescore`]).

Fix:
1. A3: new `final_eval --answers-only` mode. It writes `answers.jsonl` and the partial `final_eval.json` (step 1 of
   A3's write order), plus `answers.sha256`, then exits 0 without `score_items`, without `items.jsonl` and without the
   `FINAL_EVAL` print [code: final_eval.py:75-86 is the part it skips]. Test: no score key and no summary in the
   output folder.
2. B1.5 `b1_rows.yaml` step 4: append a third view group to the same process, `--view
   {data}/eval_views/lidc-test-shifted --split test --out {rd}/sealed_b1_test --answers-only`. Cost: 2,200 x 0.193 s
   = 7.1 min = $0.08 in the same process (no new fixed load) [estimate: section 6 inputs].
3. 00_front 6.2: item 3 becomes 17.5 min; the eval tail is 26 + 7.1 = 33.1 min. Keep `--budget-usd 2.05` and the
   260 steps; the V1 amount becomes $2.46 (tail $0.41 buys at least 36.2 min >= 33.1 + 3 min of cap rounding), cap
   floor((2.46 - 0.10) / 0.680) = 3.45 h = 207 min >= 172.3 + 33.1 = 205.4; worst case $2.46 + $0.11 = $2.57;
   `launch_refusal`: $16.90 + $2.57 = $19.47 <= $19.90; with the reserve $2.82 <= $3.00 [estimate: arithmetic at r_s].
   6.1 step 4c: the V1 bundle holds `lidc-test-shifted` (frozen, hash-checked) and no other `lidc-test-*` folder; 6.4
   rule 7 already allows a named LIDC view. 6.1 step 4d's grep test allows `--split test` only together with
   `--answers-only`.
4. 00_front 1.1 and "End state under (b')", 99_back decision 17: "Under (b') the test read is `final_eval --rescore
   runs/<V1>/sealed_b1_test` on the Mac, after the owner closes the 4B gate, and only if B1 is the final candidate;
   when R5 is the final candidate, R5's measured test stays the number and the sealed answers stay unscored."
   The Mac scoring uses the CPU build of the same engine, which `tests/test_gpu_scoring.py` holds bit-identical to
   the CUDA build on every term. Under a raise, V3 item 5 still runs for a B3 candidate.
5. If the owner rejects the extra $0.06 of amount, write the consequence instead, in decision 1 (b'): "No test number
   for B1 or any later model until a raise; README keeps R5's 0.645."

### 3. MINOR. "pool size 1 ... so no forked worker opens its own CUDA context" is not what pool size 1 does

Sources: A (A3 `final_eval` change bullet; A2 step 2), 00_front (6.2 "Why step 0").

A3: "on a CUDA build ... `final_eval` caps the pool at 1 worker ..., so no forked worker opens its own CUDA context."
With pool size 1 there is still one forked worker, and `score_items(pool, ...)` sends every scoring job to it
[code: experiments/phase2/final_eval.py:69, 76]; A2 step 2 says each worker uploads its own tables
(`ExpertCache.device` is per process). So exactly one worker opens a CUDA context. This is safe (the fork happens
before the parent touches CUDA, `serve.fork_pool` guard [code: src/objdel/train/serve.py:167-169]), and it costs one
context of memory, so nothing breaks. Fix the wording in A3 and 6.2: "at most one forked worker opens a CUDA context
(one context, created after the fork); D7 step 6 removes it." D-G17's "no new context" rule applies only after D7.

### 4. MINOR. The front header still says round 3

Sources: 00_front (title block).

"Revision: final round 3 plus the round-4 consistency pass." The doc is at round 13 and its review log runs to final
round 12. A reader of the built doc cannot tell which revision they hold. Fix: have `assemble.py` write the round
number into the header, or edit the line to "Revision: final round 13. Review logs: `docs/implementation2.0/review-log.md`."

### 5. MINOR. E3/E5 move table still lists `grpo_rows_views`, which B1.5 deletes before E0

Sources: E (E5 move table row `rl/rows.py`), B (B1.5 "Delete the `--train-views` path in the same commit").

E: "`rl/rows.py` | `dataset.py`: `prompt_rows`, `grpo_rows`, `grpo_rows_views`, `OpenOnce`". B1.5: "delete
`dataset.grpo_rows_views`, `serve.VIEW_SEP` and `_tagged`" before the tag `v1-b1`. The `verify/store.py` row already
notes the B1.5 deletion of `VIEW_SEP` and `_tagged`. E's AST test that the table names every top-level name of
`dataset.py` on the E0 commit fails at $0 on the extra name, so this is found during implementation. Fix: drop
`grpo_rows_views` from the row and add "(`grpo_rows_views` deleted by B1.5)".

## Counts

blocking 1 (item 1), major 1 (item 2), minor 3 (items 3, 4, 5).
