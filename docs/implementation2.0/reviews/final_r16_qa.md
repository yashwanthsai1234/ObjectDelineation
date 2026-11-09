# Final review round 16: QA and technical writing (Tests sections, gates, testability, cross-section consistency)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | This round makes no NVIDIA product choice, so I recommend no catalog skill. It sent me back to the round 15 engine-counter change (D.4 "Engine counters", D-G14, P03, I19, S12, S13, B04). I read it against `cuda/gpucheck.cpp`. Result: items 4 and 5. |
| engineering:code-review | I read the code that the cited lines name. `Makefile:1-40` (today's `test`, `test-all`, `cov`; there is no `test-data` target yet). `pyproject.toml` markers (`slow`, `data`). `experiments/phase2/vast/policy.py:373-420` (`code_files`, `cmd_sync_up`: the code lands at `/workspace/src/`), `:498-506` (`cmd_pull`), `:598-606` (`--idle-minutes` 30). `experiments/phase2/train_rl.py:50-68` (`rl.json` is written in `--out`). `src/objdel/train/configs.py:178-188`. TRL `grpo_config.py:935-942`. |
| engineering:debug | I traced three paths step by step. (1) Phase A's `make headline TEST=... VAL=...` through the DOC6 recipe, DOC.9.5 F3 and I05 (item 1). (2) A `slow` + `data` test (D-R08) through each `make` target that the plan defines (item 2). (3) CLEAN-G17's grep paths after E's moves (item 3). |
| engineering:system-design | One owner per fact. Two make targets have more than one definition: `headline` (DOC6 against F3/I05) and `test-data` (A.9 against B.9.0, D8 and CLEAN.4). Items 1 and 2. |
| engineering:architecture | Round 15's `val_halves.json` path, checked on the box layout. `code_files()` is rsynced to `/workspace/src/` [code: policy.py:398]. So `parents[3]` of `/workspace/src/src/objdel/data/splits.py` is `/workspace/src`, which holds `docs/reference/val_halves.json`. I50 (d) asserts that `__file__` lies under the copy. No finding. |
| engineering:documentation | Carry-over check. I re-read the round 13-15 qa minors and checked each one with grep. All 13 are still open (list at the end). |
| huggingface-skills:huggingface-best | This round makes no model choice. I checked one HF claim only. U-RC6 downloads `processor_config.json` and `chat_template.jinja`. Saved checkpoints hold these exact file names [measured: `runs/rl/train_rl/final/`, `runs/r3/ckpt/`]. No finding. |
| engineering:testing-strategy | For each gate, I asked two questions. Which make target runs its tests? Is that target defined the same way everywhere? A `slow` + `data` test runs under no target that a gate names (item 2). DOC-G6 never checks `ci_unit`, so the item CI can be published while every test passes (item 1). |
| mattpocock-skills:tdd | A gate must be able to go red. D-G14's `live_device_bytes <= arena_bytes + table_bytes` cannot go red when the counter undercounts (item 5). P03's `eval_sets` cap of 4,096 B conflicts with its own "raises ValueError" clause (item 4). |

Web tools: loaded (WebSearch, WebFetch). No new external claim needed a web check this round. The one HF file-name
claim was checked against local checkpoints (above).

Measurements this round ($0, Mac, 2026-11-09):
- `docs/implementation2.0/Implementation.md` is 1,765,919 bytes. `review-log.md` is 385,693 bytes (`ls -l`). Both are
  under the DOC-G13 limits.
- No match of `/[p]rivate/tmp` or `/[U]sers/` in either file (`grep -cE`: 0 and 0).
- Structural id checks (scripts in `reviews/r16qa/`). Every `A-G`, `D-G`, `E-G`, `DOC-G` and `CLEAN-G` id that the
  doc cites has its own gate row. In D, the undefined test-like ids are only task ids, note ids and E's `U8`.
- 47 `python` fences: 45 parse with `ast.parse`. Two do not: E.4 E4 step 1 `params_for` and the E2.1 `load_spec`
  stubs. Both are signature sketches, not test sketches. No item.
- `runs/r5/eval/eval_test/final_eval.json` has no `schema` key, and `iou_tc` holds `mean`, `ci95`, `median`
  (schema 1, item CI). `runs/r5/eval_a1/` does not exist yet.

## Items

### 1. major: `make headline` has three forms that disagree. A's command fails under F3/I05, and the only form F3 allows publishes the item CI where A9 needs the patient CI (A, DOC)

The doc lines, in order:
- A.md, A1 step 5b: "Then `make headline TEST=runs/r5/eval_a1/eval_test/final_eval.json
  VAL=runs/r5/eval_a1/eval_val/final_eval.json` (DOC6) publishes the patient CI of test IoU_tc, [0.6252, 0.6644],
  with `ci_unit: patient` (A9)."
- A.md, A9: "`docs/results/headline.json` regenerated with `ci_unit: patient` (DOC6): `make headline TEST=... VAL=...`".
  The same form appears in DOC.md "Docs per phase" row A and DOC.8.
- DOC.md, DOC6, `Makefile` bullet. The recipe is `$(if $(RUN),--run $(RUN),--test $(TEST) --val $(VAL)) --model-ref
  $(MODEL_REF)`. It accepts both forms.
- DOC.md, DOC.9.5 F3 ("Findings that change the spec"): "Recipe first validates `RUN` and `MODEL_REF`
  (`case "$(RUN)" in ""\|*[!A-Za-z0-9_.-]*) ... exit 2`)". An empty `RUN` therefore exits 2.
- DOC.md, DOC.9.3 I05: "`make headline` without `RUN`, with `RUN=` or empty `MODEL_REF` exits non-zero".
- DOC.md, DOC.3.2 step 5: "After a new evaluation, the owner runs `make headline RUN=<id>` ... It copies the metrics
  of `runs/<id>/eval/eval_test/final_eval.json` and `eval_val/final_eval.json`".
- DOC.md, DOC6, `headline.py` bullet: "Why two paths: A3 writes `{rd}/eval_val` in the V1 run and `{rd}/eval_test` in
  the V3 run ..., so no one run folder holds both."
- DOC.md, DOC.9.3 I10: "If `runs/r5` exists, regenerating with `main` gives the same `metrics` as the tracked file".
- DOC.md, DOC.6 DOC-G6: "README contains all strings from `headline_strings(headline.json)`". It never checks `ci_unit`.

What goes wrong:
- A's two commands give no `MODEL_REF` and no `RUN`. Under F3 and I05, both exit 2.
- DOC.3.2 step 5 tells the owner to use `RUN=<id>` after every new evaluation. That form reads `runs/<id>/eval/...`,
  which is only the R5 layout. DOC6 itself says that no A3 run folder holds both files.
- So the one form that F3 allows is `make headline RUN=r5`. It reads the schema-1 file, which holds the item CI
  [measured: `runs/r5/eval/eval_test/final_eval.json` has no `schema` key; `iou_tc` keys are `mean`, `ci95`,
  `median`]. The result is `ci_unit: "item"`, not the patient CI that A1 step 5b and A9 require.
- Every planned check still passes. DOC-G6 compares only strings. I10's regeneration rule (`main` on `runs/r5`) then
  agrees with the item-CI file. And if the engineer builds the correct schema-2 file from `runs/r5/eval_a1/`, I10
  goes red, because its regeneration does not name a source.
- No test covers the `TEST`/`VAL` form. I03 and I05 test only `RUN=r5 MODEL_REF=x`.

Fix (DOC.md, plus two A.md lines):
- F3 and I05: validate `RUN` only when it is set. The recipe exits 2 when (a) `RUN` and `TEST`/`VAL` are both given,
  (b) neither is given, or (c) `MODEL_REF` is empty. `TEST` and `VAL` must match `^runs/[A-Za-z0-9_./-]+\.json$`
  with no `..`.
- I05: add the case `make -n headline TEST=runs/a/t.json VAL=runs/b/v.json MODEL_REF=x`. It prints `--test`, `--val`
  and `-I`.
- DOC.3.2 step 5 and the DOC.3.1 box label: "`make headline TEST=<eval_test file> VAL=<eval_val file>
  MODEL_REF=<repo@revision>`. `RUN=<id>` is only for the R5 layout."
- A1 step 5b and A9: append `MODEL_REF=svenk029/objdel-sft@rl-full-3ep` (the I03 value) to both commands.
- I10: "If `runs/r5/eval_a1` exists, regenerating with `main --test ... --val ...` from the `source` paths that are
  recorded in the file gives the same `metrics`."
- DOC-G6 (or a new I10 line): "`ci_unit == "patient"` once A1 has landed".

### 2. major: `make test-data` has two definitions. Under A's definition, which lands first, the `slow` + `data` test D-R08 runs under no target that a gate names (A, B, D, CLEAN)

The doc lines:
- A.md, A.9 (the Makefile edit; A owns it because A lands first): "new `test-data: $(EXT_STAMP)` runs `$(PY) -I -m
  pytest -q -m "data and not slow"`; new `test-all-data` ... `-m data` (the slow `data` tests too)".
- B.md, B.9.0: "runs only through `make test-data`, which is `$(PY) -I -m pytest -q -m data`".
- D.md, D8: "`make test-data` (`$(PY) -I -m pytest -q -m data`) ... and `test-all-data` runs `$(PY) -I -m pytest -q
  -m data`". This gives the same command to two targets.
- CLEAN.md, CLEAN.4 task rules: "`make test-data` (`$(PY) -I -m pytest -q -m data`)".
- D.md, D.9.4 D-R08: "`data/checker_eval/masks.npz` read with `-I`: 37,093 crops ... every crop equals the oracle,
  parts equal cv2, `slow`, `data`, CPU build on the Mac". Its gate is D-G3. D8 says that `data` tests "run only
  through `make test-data`".

Proof of the gap:
- Today's Makefile has no `test-data` target. Its `test` target is `-m "not slow"` [code: Makefile:19-20].
- A.9 changes `test` to `-m "not slow and not data"`, and `test-data` to `-m "data and not slow"`. So D-R08 is
  deselected by `make test`, by `make test-data` and by `make cov`.
- Only `test-all-data` runs D-R08, and no gate or run step in the doc names `test-all-data` (grep: its only uses are
  the A.9 and D8 definitions).
- A deselected test is not a failure, so D-G3's Mac half can close without D-R08 having run. Under (b') the box half
  of D-G3 (R-D2, V3) never runs. Then the 37,093-crop oracle and cv2 check of the new polygonizer has no planned
  run. D-R01 still compares target hashes, but it does not check the oracle or the cv2 parts.
- E.9 X2 (`data`, `slow`) has the same exposure. Its gate file (E-G3, `tests/core/test_move_invariance.py`) does not
  name it.

Fix:
- One definition, A.9's. B.9.0, D8 and CLEAN.4 replace their restated commands with "`make test-data` and
  `make test-all-data` as defined in A.9".
- D8 deletes "and `test-all-data` runs `$(PY) -I -m pytest -q -m data`" from its own list, or keeps it as the only
  copy and drops the `test-data` one.
- D-G3 (Mac half) adds "and `make test-all-data -- -k "d_r08"` (or `make test-all-data`) passes on the Mac". The D-R08
  row says "run by `make test-all-data` (slow)".
- E.9 X2 names its run command in the same way (`make test-all-data`).

### 3. minor: CLEAN-G17 greps pre-E paths, but Wave 3 runs after E, which moves `TRLReward` to `rl/reward.py` (CLEAN)

- CLEAN.md, CLEAN.6: "CLEAN-G17 | `grep -n "f.trl_reward = self" src/objdel/train/rewards.py`; `grep -n "gc.freeze()"
  experiments/phase2/train_rl.py` | 1 hit each". CT-09: "Source gate: `grep -n "f.trl_reward = self"` gives 1 hit".
- E.md, E5 move table: "`rl/reward.py` | `rewards.py`: `TRLReward`, ...". E.md, step 4 rewrites `train_rl.py` as a
  thin CLI. CLEAN.8: "Wave 3 needs E2, E6 and D7 step 6".
- After E the first grep names a file that no longer exists. CLEAN-G13 already switches its path by wave, but
  CLEAN-G17 does not.
- Severity: the grep fails on the Mac at $0 and names the missing file.
- Fix: in CLEAN-G17 and CT-09, give the path by wave: Waves 1-2 use `src/objdel/train/rewards.py`; Wave 3 uses
  `src/objdel/rl/reward.py`. For `gc.freeze()`, use the file that holds the RL run body after E2 (E's move table),
  or replace both greps with the AST check that CT-09 and CT-35 already run.

### 4. minor: P03's `eval_sets` cap of 1<<12 bytes conflicts with its own "too small raises ValueError" clause and with "peak > 0" (D)

- D.md, D.9.2 P03: "a case too small for one mask raises ValueError, never wrong data; same for `eval_sets` caps
  {1<<12, 1<<16, 1<<28}, where `gc.stats()["eval_scratch_peak"]` after `reset_stats()` is > 0 and <= the cap".
- D.md, D.4: `eval_scratch_peak` is "the largest scratch that one `eval_sets` chunk asked for ... 0 until the first
  `eval_sets` call".
- The doc gives no P03 mask size. One 512 x 512 bit mask is 32 KiB, which is more than 4,096 B. Then the 1<<12 case
  raises. The doc does not say whether a refused chunk raises `eval_scratch_peak`, so "> 0" may hold or fail.
- Severity: the test is Mac-only and goes red at $0.
- Fix: P03 names the `eval_sets` fixture size (for example 40 masks of at most 32 x 32). It states: "a cap below one
  set's scratch raises ValueError and leaves `eval_scratch_peak` unchanged; the peak assertion applies only to caps
  that succeed".

### 5. minor: D-G14, B04 and S13 gate `live_device_bytes <= arena_bytes + table_bytes`, but D.4 defines equality, so an undercounting counter passes (D)

- D.md, D.4: "After D5c every device block is an arena or a table, so `live_device_bytes == arena_bytes +
  table_bytes` in a correct engine." I19 (Mac) asserts `==`.
- D.md, D.6 D-G14, D.9.6 S13 and D.9.7 B04: "`live_device_bytes` <= `arena_bytes` + `table_bytes` (no hidden device
  allocation)".
- With `<=`, a `live_device_bytes` that misses a block (for example a pool block outside the helpers) still passes.
  That is the defect the record exists to find.
- Fix: use `==` in D-G14, S13 `test_b04` and B04, as in I19.

### 6. minor: the V1 cap arithmetic reads 3.39 h; (2.40 - 0.10) / 0.680 = 3.38 h (00_front)

- 00_front.md, 6.2: "At the planning values the cap is 3.35 h = 201 min [estimate: (2.40 - 0.10) / 0.680 = 3.39 h,
  rounded down]."
- 2.30 / 0.680 = 3.382 [arithmetic]. The floor to 0.05 (3.35 h) and the $2.49 worst case (0.10 + 3.5167 x 0.680 =
  2.491) are correct.
- Fix: "= 3.38 h".

## Carry-overs still open (minor; checked with grep this round; not re-argued)

- r15 qa-1: A3's header builds `lidc-test-shifted`, then bundles, and no line moves the test view out (A.md A3 header;
  00_front 6.1 step 4c still offers "(or a bundle `--view` filter is added)").
- r15 qa-2: D-R01's always-run golden names no pre-D source (D.md D8 payload item 3, D-R01).
- r15 qa-3: `rl.json` `opsm_start_delta` is 0.5 when the flag is absent, also under presets with OPSM off (B.md B2.2,
  U45).
- r14 qa-2: A2, A8, C2 and C11 still read `runs/a3/eval_val/answers.jsonl`. There is no `OBJDEL_V1_RUN` (A.md:428,
  1109; C.md:651, 1093).
- r14 qa-3: E8's list of hard-coded test paths does not name `tests/variants.py`.
- r14 qa-4: D5c sizes `in` by "names per answer" (D.md:695). There is no `NAMES_MAX`.
- r14 qa-5: `GC_ALIAS_SCAN_FOR_TEST` does not name the aliased `scan` call (D.md:430).
- r13 qa-3: the B1.3 DG1b row still says "(box: the manifest's `excluded_ids`)" (B.md:1942).
- r13 qa-4: C.9.4 still has no `limit_mm_per_prompt={"image": 0}`.
- r13 qa-5: C.md:1968 still says "15.4 min < 30.9 min".
- r13 qa-7: the C10 keep list still omits `--encode-check` and `--checks-from-step`.
- r13 qa-8: E.md:645 and E.md:747 still say "moved by D3".
- r13 qa-9: 00_front.md line 4 still says "Revision: final round 3 plus the round-4 consistency pass".

## Totals

Blocking 0, major 2, minor 4 (new), plus 13 carry-over minors.
