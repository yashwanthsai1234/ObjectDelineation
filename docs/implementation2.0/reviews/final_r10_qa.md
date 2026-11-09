# Final review round 10: QA and technical writing (Tests sections, gates, testability, cross-section consistency)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice is in scope, so I recommend no catalog skill. It pointed me at the CUDA scorer binding: `score_masks` needs a device table per object, which is the cause of item 1. |
| engineering:code-review | I read `cuda/gpucheck_py.cpp:104-106, 213-224`, `cuda/gpucheck.cpp:1174-1185, 1808-1833`, `src/objdel/reward/score.py:67-77`, `tests/test_gpu_scoring.py`, `experiments/phase2/vast/policy.py:420-510`, `experiments/phase2/vast/bundle.py:48-73` and the `Makefile` against each test claim. |
| engineering:debug | I traced the V2 plan path. `policy.chain` returns a fixed list of steps, and `cmd_run` joins them with `&&` into one tmux job. So no step can depend on a gate result, and the main agent cannot act inside one `run.sh` call (item 2). |
| engineering:system-design | I matched the V2 and V3 launch lines against the priced rows and the owner totals (items 2, 3, 4). |
| engineering:architecture | One owner per rule: I compared the C9 arm home (C9, C0, C.8.2, front 6.3) and the B-S2 definition (front 6.3 table, front 6.3 bullet, B.8.1). |
| engineering:documentation | I checked stale header and baseline text (items 7, 8) and the cross-references in front 6.1 step 4d (item 9). |
| huggingface-skills:huggingface-best | No model choice is in scope. No change. |
| engineering:testing-strategy | I checked each gate for a test with a defined input, location and pass rule. Items 1 and 5 have gates whose inputs are not defined. |
| mattpocock-skills:tdd | I checked that the pytest sketches call the public interface with its real signature. The G01 sketch does not (item 1). |

Web tools: loaded (WebSearch, WebFetch, Exa). This round needed no external claim check: every item rests on repo code or a local measurement.

Measurements made this round ($0, Mac):
- Coverage with `data` deselected (the A and D Makefile change): `.venv/bin/python -I -m pytest -q -m "not slow and not data" --cov=objdel` gives TOTAL 3,866 statements, 253 missed, 93 %; 1,388 passed, 6 skipped, 90 deselected, 257 s [measured: `reviews/r10qa/cov_nodata.txt`, 2026-11-09]. So the A result note gets 93 %, and the D-G12 floor stays 90. This is not a finding. E8's coverage gate still applies: (522 + 33) missed of 4,608 = 88.0 % < 90 [computed from the E8 numbers and this run].
- `data/checker_eval/masks.npz` (read with `-I`): bit-packed `bits` plus `hw` and `start`. 37,093 crops; the first 1,000 have 155 distinct frames, from 20 x 20 up to 76 x 91. The file holds no doctor masks and no expert tables.
- `docs/implementation2.0/Implementation.md` is 1,665,519 bytes; `review-log.md` is 342,286 bytes (`ls -l`, 2026-11-09).

## Items

### 1. major: G01 and G08b score crops against expert tables that the doc never defines (B)

- Doc (B.md, B.9.5 G01): "Inputs: 200 random masks, 10,000 fuzz polygons ..., and the first 1,000 crops of
  `checker_eval/masks.npz` (through `box_input`). `sum_x`, `sum_y`: equal the `np.nonzero` sums of each mask ...
  Every other field: dtype, shape and SHA-256 of `tobytes()` equal `tests/goldens/b_g01.json`".
- Doc (B.md, B.9.5 G08b): "`test_g08b_smoke`: the G01 checks (the `np.nonzero` sums and the `windows_1000` entry
  of `tests/goldens/b_g01.json`) on the first 1,000 crops of `checker_eval/masks.npz`".
- Doc (B.md, B.9.8 sketch `tests/test_scorer_sums.py`): `r = G.score_masks(m[None])  # call shape per
  tests/test_cuda_checker.py (adapt)`.
- Code: `score_masks(bits, hw, start, tables)` needs one table handle per object [code: cuda/gpucheck_py.cpp:213-224,
  104-106]. Each handle must have the object's own frame: `"an object's frame must be its tables' frame"` [code:
  cuda/gpucheck.cpp:1181-1182]. The only Python caller builds the handle from an `ExpertCache` of the same shape
  [code: src/objdel/reward/score.py:67-77]. The existing tests build that cache from a stack of rater masks
  (`ExpertCache.build(stack)`) [code: tests/test_gpu_scoring.py:120, 135, 151].
- Measurement: `masks.npz` has no rater masks and no tables. Its first 1,000 crops have 155 different frames
  (20 x 20 to 76 x 91). The 200 random masks of the sketch have random frames from 8 to 512.
- Problem: the doc does not say which rater stack builds each table. The `iou_tc`, `f5`, `g_area`, `sig` and
  `topo_freq` fields, and so the golden SHA-256 values, depend on it. One engineer may use the mask itself as a
  one-rater cache. Another may use a fixed disk set. The two goldens differ. Also, the sketch call fails as written.
  The `np.nonzero` sum check does not depend on the tables, but "every other field" does.
- Fix (B.md B.9.5 G01, G08b and B.9.8 sketch): add one rule. "Tables for G01 and G08b: for each mask m, the
  table is `ExpertCache.build(np.stack([m, np.roll(m, 1, axis=1), np.roll(m, 1, axis=0)]))` (three raters in the
  mask's own frame; the same code in `make_b_goldens.py` and in the test). The crops are passed to `score_masks`
  as stored (`bits`, `hw`, `start` slices of the first 1,000 entries), with one handle per crop." Change the sketch
  to build these caches and call `G.score_masks(bits, hw, start, np.array([c.device() for c in caches], np.int64))`.
  Add the same rule to the `make_b_goldens.py` text. Then the Mac golden and the box test use the same input.

### 2. major: the V2 stop rule and the conditional V2 arms cannot run inside one `run.sh` call (00_front, C)

- Doc (00_front.md, 6.3 V2 table): item 3, stop rule "Step-1 rewards differ: stop the fast arm"; item 3c "Optional
  mode B arm, only if mode A fails C-G-MEM"; item 3d "Optional C6 arm ... Skippable (6.4 rule 2)".
- Doc (00_front.md, 6.3 V2 launch lines): `run.sh <run> experiments/phase2/configs/c_speed.yaml  # items 2-3 (3c,
  3d only if approved); wait for EXIT`, and "one `run.sh` call per item, in table order; between two calls the main
  agent applies 6.4 rule 2".
- Doc (00_front.md, 6.4 rule 2): "`policy.chain` has no spend-based skip (a leading `~` only lets a step fail ...).
  So the main agent applies this rule on the Mac, between two `run.sh` calls".
- Doc (C.md, C.9.5 E-EQ6): "a step-1 mismatch is a G-EQ6 fail: Gate C-adopt fails ... and front 6.3's stop rule
  stops the fast arm".
- Code: `chain` returns the `steps:` list as written, with only `{rd}`/`{data}` filled in [code:
  experiments/phase2/vast/policy.py:432-433]. `cmd_run` joins every step with `&&` into one `bash -c` job under one
  `timeout` [code: policy.py:487-494]. No step can read a gate result of an earlier step and skip itself.
- Problem: items 2, 3, 3c and 3d are one plan file and one `run.sh` call. So (a) nothing stops the fast arm after a
  step-1 mismatch. It and its check steps run to `--steps 40`, about 19 min = $0.22 after C is already rejected
  [estimate: C.8.2 8 + 7 + 4 min]. (b) 3c "only if mode A fails C-G-MEM" needs a result from the same plan. (c) 3d
  "only if approved" has no way to leave the plan. An engineer must guess: put 3c and 3d in the plan (unpriced
  spend when they are not needed), or leave them out (no mode B result when G-MEM fails).
- Fix (00_front 6.3 V2, C.md C0 files and C.9.5): split the plan. `c_speed.yaml` holds items 2 and 3 only.
  Add `c_speed_b.yaml` (item 3c) and `c_speed_c6.yaml` (item 3d), each with its own `run.sh` line marked "only if
  ...", so 6.4 rule 2 and the G-MEM condition apply between calls. For the step-1 stop: the fast arm's `train_rl`
  gets `--expect-step1 {rd}/c0_base/rl.json`. After step 1 it compares the three G-EQ6 reward numbers with the
  baseline's step 1 and exits non-zero on a mismatch, which stops the `&&` chain. Add a Mac unit test (fake
  `rl.json`, planted mismatch, exit code non-zero before step 2). U-CF4 then also parses the two new plans.

### 3. major: the C9 arms have two homes, and V2 does not price them (C, 00_front)

- Doc (C.md, C9): "Files: `experiments/phase2/configs/c_speed.yaml` arms; `train_rl --no-ckpt-layers K`".
- Doc (C.md, C0 files): "each 17-step arm (baseline, 3c, 3d, C9) runs `--steps 20`".
- Doc (C.md, C.8.2): "C7, C8, C9 and C11 need no rental of their own: they ride inside the Phase B runs (C9 and C11
  arms only if B's budget allows; each 17-step arm ... about 8 min plus 6 min start-up)".
- Doc (00_front.md, 6.3 V2 table): rows 1, 1b, 2, 3, 3c, 3d, 4, 5 and no C9 row. V3 has no C9 row either.
- Problem: the engineer who writes `c_speed.yaml` from C9 and C0 adds three arms (k = 4, k = 8, 32 x 8): about
  3 x 14 = 42 min = $0.48 that no V2 row holds [estimate: C.8.2 rates at r_s]. The cap is computed from the approved
  amount, so the cap then ends the plan job inside an arm, or the arms use the money of items 4-5. The engineer who
  follows C.8.2 leaves them out. The two readings give different V2 costs and results.
- Fix (C.md C9 and C0; front 6.3): pick one home. Proposed: C.8.2 wins. C9 "Files" says "a separate plan file
  `c_mem.yaml`, run only if the owner adds C9 rows to a rental table (front 6.3); not in `c_speed.yaml`". C0 drops
  "C9" from the arm list. E-C9 in C.9.5 gets "Optional, not in V2" like E-REP.

### 4. major: three definitions of B-S2; the V3 amount does not cover the V3 launch lines (00_front, B)

- Doc (00_front.md, 6.3 V3 table): "**B-S2 = items 0 + 4** | **182.6 = $2.17** | **165.9 = $1.98**".
- Doc (00_front.md, 6.3 bullet under it): "B-S2 in B.8.1 is V3 items 0, 2a and 4 (plus 5 when the gate is closed);
  B.8.1 points here."
- Doc (B.md, B.8.1): "**B-S2 = rental V3 items 0, 2a, 2b and 4, front 6.3** (plus item 5 ...)".
- Doc (00_front.md, 6.3 "Totals for owner decision 1"): "(a) B1 then B3: V1 + B-S2 + reserve = $2.29 + $2.17 +
  $0.25 = $4.71".
- Doc (00_front.md, 6.3 V3 launch lines): one `run.sh` call each for item 1 (`d_gpu_test.yaml`), 2b, item 3
  (`a3_test.yaml`) and item 4 (`b3_v5.yaml`); "skipped items have no call". Items 1 and 3 have no approval
  condition. 6.4 rule 2 skips a skippable item only once spend reaches 90 % of the amount.
- Doc (00_front.md, 6.4 rule 5): the `BudgetStop` value is "approved amount - $0.35 ... - the $ of every later item
  in the same rental". `BudgetStop` counts session spend from launch [code: policy.py:190-193, as the front cites].
- Problem: the front table, the front bullet and B.8.1 name three different item sets for one label. The owner
  price for option (a) ($2.17) covers items 0 and 4 only. The V3 launch lines also run items 1 and 3, $0.46-0.67
  before item 4 [front table: 0.31-0.52 + 0.15]. With an approved V3 amount of $2.17, the B3 `BudgetStop` sits at
  $1.82. Items 0, 1 and 3 use about $1.0-1.2 first, so the B3 run stops near step 85-110, not 260 [estimate:
  (1.82 - 1.22) / 0.01133 = 53 min, minus 14 min of build, smoke and start-up, at 27.5 s per step]. Gate B3 then
  compares a short B3 run with the 260-step B1 run. Also, no line states the B3 `--budget-usd` value: B3.9 says
  only "BudgetStop (session spend) on".
- Fix: (00_front 6.3) one definition: "B-S2 = V3 items 0 and 4; items 1, 2a, 2b, 3 and 5 are separate optional
  rows, each priced in this table". Add a V3 total row "items 0, 1, 2a, 3, 4 (default V3)" and use it in the
  totals for option (a), or mark items 1 and 3 "only if approved" in the launch lines. Change the bullet to match.
  (B.md B.8.1) "B-S2 = V3 items 0 and 4 (front 6.3); 2a, 2b and 5 are optional rows". (B.md B3.9) add "`train_rl
  --budget-usd <V3 approved amount - 0.35 - $ of item 5 when planned>` (front 6.4 rule 5), written in `b3_v5.yaml`
  before preflight".

### 5. major: the Mac half of G-EQ2 and G-EQ3 has no V2 run-folder rule and no batch folder (C)

- Doc (C.md, C.9.4 runner table): "G-RP2, G-RM1 (equality part) | Mac, after `pull.sh` of V2 | `make test-data` on
  the pulled box batches plus the A3 Mac batches (marker `data`, not `gpu`)".
- Doc (C.md, C.9.5 G-RP2): "Runs on the Mac (marker `data`) after V2 is pulled. At least 50 batches (at least 15
  from the box: the check-arm `--record-batches` files and the 5 `c_speed_check` batches)".
- Doc (C.md, C.6 check arm): "`--record-batches DIR`". No line names `DIR` or the folder of the 5 `c_speed_check`
  batches.
- Contrast: S13 reads `runs/$OBJDEL_V3_RUN` and skips without the variable (D.md D.9.6 S13). GP5 reads
  `OBJDEL_V4_RUN` (E.md E.9). The round 9 fix gave V3 and V4 this rule; V2 has none.
- Code: `pull.sh` copies the box folder `/workspace/runs/<run_id>/` to `runs/<run_id>/` [code:
  experiments/phase2/vast/policy.py:498-507]. So the Mac path depends on a run id that no test knows.
- Problem: the engineer must guess where G-RP2 and G-RM1 look and what they do with no V2 run. Under (b') V2 never
  runs. Then `make test-data` (A.6 "All" gate) either fails on "at least 15 box batches" or skips with no stated
  reason. Gate C-adopt reads this half (C.6: "the Mac half of G-EQ2/G-EQ3 passes on A3's answers before V3").
- Fix (C.md C0 files, C.6, C.9.4, C.9.5): name the folders: `--record-batches {rd}/c_check/batches` and
  `c_speed_check --batches-out {rd}/c_speed_check/batches`. G-RP2 and G-RM1 read `runs/$OBJDEL_V2_RUN/c_check/batches`
  and `runs/$OBJDEL_V2_RUN/c_speed_check/batches`. Without `OBJDEL_V2_RUN` they skip with "no V2 run pulled
  (OBJDEL_V2_RUN unset)" (A.9.6 rule). With it, fewer than 15 box batches fails. Front 6.3 V2: after `pull.sh`, add
  `OBJDEL_V2_RUN=<run> make test-data`, as for V3.

### 6. minor: a later `final_eval` can overwrite the D-G17 record, and S13 does not check which eval wrote it (D, 00_front)

- Doc (D.md, D7 step 5): `final_eval` "writes `$OBJDEL_RUN_DIR/d/topology.json`: per view, the method, the count
  before scoring ...". The thread runs whenever `OBJDEL_RUN_DIR` is set and the build is CUDA.
- Doc (00_front.md, 6.3 V3 lines): `a3_test.yaml` (item 3, a `final_eval` on `lidc-test-shifted`) runs before
  `b3_v5.yaml` (item 4) in the same run folder.
- Doc (D.md, D.9.6 S13): "`test_b06` reads `runs/$OBJDEL_V3_RUN/d/topology.json` (D-G17 (a) and (c))".
- Problem: if the B3 val eval does not run (an error before it), S13 reads the A3 test eval's record and can pass
  D-G17 on the wrong process.
- Fix (D.md D7 step 5 and S13): write the file per output folder (`<--out>/topology.json`), or keep one file and
  key it by view. `test_b06` asserts that the record holds the view `lidc-val-shifted` from `eval_b3_val_shifted`.

### 7. minor: E-MEM and E-SP6 test rules that C.6 does not hold (C)

- Doc (C.md, C.9.5 E-MEM): "peak at most 30 GB on every V2 arm; drift under 100 MB across steps 3-15". C.6 G-MEM:
  "Peak GPU memory <= 30 GB on every arm" (no drift rule).
- Doc (C.md, C.9.5 E-SP6): "At least 70 % over steps 16-17". C.6: "G-SP2..6 and G-MEM are reported".
- Fix: E-MEM reports the drift as a report line ("drift reported, not gated"), or C.6 G-MEM adds it. E-SP6 says
  "reported against the 70 % target; never a failure of Gate C-adopt".

### 8. minor: `make test-data` has two filters (A, D)

- Doc (A.md, A.9 conventions): "new `test-data: $(EXT_STAMP)` runs `$(PY) -I -m pytest -q -m "data and not slow"`".
- Doc (D.md, D8 files): "`make test-data` (`$(PY) -I -m pytest -q -m data`)".
- D adds the target only "if no earlier phase did", so A's form wins, but D's text misstates it.
- Fix (D.md D8): "`make test-data` (A.9: `-I -m "data and not slow"`; `test-all-data` runs the slow ones)".

### 9. minor: front 6.1 step 4d still cites a grep test (00_front)

- Doc (00_front.md, 6.1 step 4d): "A grep test in `tests/test_repo_hygiene.py` checks it."
- The round 9 verdict (qa-5, minor) found that A I20, B I22a, B I22b, C `c_speed.yaml` step 0 and E U51 already
  own these checks. The text is unchanged.
- Fix: replace that sentence with "A I20, B I22a and I22b, C U-CF4 and E U51 check it".

### 10. minor: G01 joins the sanitizer corpus with no size limit (B, D)

- Doc (D.md, D8): every other sanitizer case has a cap (for example "`test_san_check_masks` (<= 200 masks ...)").
  "Once B3.1 lands, its `tests/test_scorer_sums.py` joins this corpus; its GPU cases carry `gpu` and `sanitizer`."
  Each tool has a 300 s `timeout`.
- Doc (B.md, B.9.5 G01): 200 masks, 10,000 fuzz polygons and 1,000 crops.
- Problem: under memcheck or racecheck, 11,200 cases may pass 300 s, and R-D2 stops on "Any failure". This is not
  measured, so it stays minor.
- Fix (B.md G07, D.md D8): the `sanitizer` case of `test_scorer_sums.py` uses the 200 random masks only. The
  10,000 polygons and 1,000 crops carry `gpu` only.

### 11. minor: stale header and baseline (00_front, DOC)

- Doc (00_front.md, line 4): "Revision: final round 3 plus the round-4 consistency pass." This is round 10.
- Doc (DOC.md, DOC1 and DOC-G13 baseline): "1,480,096 bytes [measured: `ls -l`, 2026-11-08]". Measured today:
  1,665,519 bytes (83 % of the 2 MB limit); `review-log.md` 342,286 bytes.
- Fix: update the revision line; update the baseline cell. Optional: DOC-G13 names `review-log.md` with the same
  1 MB rule.

## Not repeated

The round 9 fixes (qa-1 G08b `bash -c` with `OBJDEL_DATA_ROOT={data}`, qa-2 one G01 file with a Mac golden, qa-3
E-SP1 report lines, qa-4 optional C.9.5 arms, qa-6 `OBJDEL_V3_RUN`) are in place. I checked them against
`policy.py:432-494`, `bundle.py:66-70` and the B, C and D texts. Item 1 is a new defect in the qa-2 golden design
(the tables), not a repeat.

Totals: 0 blocking, 5 major (1-5), 6 minor (6-11).
