# Final review, round 9: QA lead and technical writer (Tests sections, gates, testability, cross-section consistency)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice is in scope. It pointed me at the box CUDA-build and sanitizer chain, so I checked which build sits in `/opt/venv` when each box test runs (items 1, 2). |
| engineering:code-review | I read `policy.py:379, 417, 484-492`, `run_box.sh`, `Makefile:41-45`, `CMakeLists.txt:15`, `bundle.py:66-70`, `tests/test_gpu_scoring.py:1-20` and `pyproject.toml:66-71` to check each box-test claim against the code. |
| engineering:debug | I traced V3 item 4 step by step: sync drops `data/`, the `-m pytest` step gets no data root, the `&&` chain stops (item 1). |
| engineering:system-design | I matched each Tests-section arm against the V2 price rows of front 6.3; three test arms have no row (item 4). |
| engineering:architecture | One owner per rule: G-SP1 (C.6) and E-SP1 (C.9.5) both state the decision rule and disagree (item 3). |
| engineering:documentation | I checked that each gate names its inputs and its test. The 6.1 step 4d test has no id and no owner (item 5). |
| huggingface-skills:huggingface-best | No model choice is in scope. No change. |
| engineering:testing-strategy | I checked that each box test has a reference that exists on the box. A CPU-build reference cannot exist there (item 2). I checked skip rules against `make test-data` (item 6). |
| mattpocock-skills:tdd | I checked that the expected values come from an independent source. G08b and G01 need a CPU-build answer that the box cannot make; the fix uses `np.nonzero` sums plus a Mac golden (item 2). |
| Web (WebSearch, WebFetch loaded) | I checked the Nsight Systems user guide: `--cuda-memory-usage` "may cause significant runtime overhead", and the `--cuda-graph-trace` default is `graph`. Both match C0 step 3 and step 4. No item comes from this. |

Measurements made for this review:
- R5 per-step time over steps 3-15 (`step_time` in `runs/r5/base_rl.json` and `runs/r5/shared_rl.json`): base mean 30.78 s, SD 1.89 s; shared mean 23.43 s, SD 2.11 s [measured: 2026-11-09, `statistics.stdev` of 13 values].
- R5 found items (K-mean IoU_tc >= 0.1): val 947 of 1,084, test 1,962 of 2,200; this matches A-G2 and G04 [measured: `runs/r5/eval/eval_{val,test}/items.jsonl`].
- `data/checker_eval/masks.npz` exists on the Mac (keys `bits`, `hw`, `start`, ...; 37,093 masks of different sizes) [measured: `np.load`].

Sources: phases 00_front (6.1, 6.2, 6.3), A (A3, A.6, A.9.4, A.9.6), B (B1.5, B2.3, B3.1, B3.9, B.6, B.9.5, B.9.8), C (C0, C.6, C.8.2, C.9.4, C.9.5), D (D7, D8, D.6, D.9.6), E (E.9 U51, U52); review log rounds 1-8.

## Items

### 1. blocking | sources: B (also D) | The G08b smoke step cannot find its input on the box, so the paid B3 run never starts

Problem.
- B "#### B3.9 B3 run": "Step 1 is the G08b smoke, before `train_rl` ... `[-m, pytest, -q, -p, no:cacheprovider, -m, gpu, tests/test_train_fast.py::test_g08b_smoke]`. Its 1,000 real windows are the first 1,000 crops of `{data}/checker_eval/masks.npz`". B "### B.9.5" runner row says the same.
- `{data}` is a plan placeholder. It is filled only in the step argv [code: policy.py:432-433], and this argv holds no data path. The test itself cannot see `{data}`.
- D "### D8" fixes the convention for every `gpu` test that opens a file: "A `gpu` test opens a file only through a helper `box_input(name)` ... it FAILS (never skips) when the file is missing". The data root is `OBJDEL_DATA_ROOT`, "on the Mac the default is `data`". `d_gpu_test.yaml` sets it explicitly in a `bash -c` step: "`OBJDEL_DATA_ROOT={data} /opt/venv/bin/python -I -m pytest ...`".
- The G08b step is a `-m` step. It runs as `/opt/venv/bin/python -I -m pytest ...` [code: policy.py:486], so it cannot set an environment variable. Only `OBJDEL_RUN_DIR` is exported [code: policy.py:492].
- With the default root, the test reads `/workspace/src/data/checker_eval/masks.npz`. The code sync drops `data/` [code: policy.py:379 `drop = re.compile(r"^(data|docs)/...")`]. The bundle unpacks to `/workspace/data/lidc` [code: policy.py:417].
- Result: the test fails (D8 rule). The `&&` chain stops [code: policy.py:490]. V3 item 4 pays launch, sync, build and smoke (front 6.3: 39 + 3 + 5 min, about $0.6), and `train_rl` never runs. If the engineer instead lets it skip, the smoke passes with nothing checked.

Fix.
- B3.9, B.9.5 runner row, and I22b: make step 1 a `bash -c` step: `[bash, -c, "OBJDEL_DATA_ROOT={data} /opt/venv/bin/python -I -m pytest -q -p no:cacheprovider -m gpu tests/test_train_fast.py::test_g08b_smoke"]`.
- B.9.5 G08b: the test reads the file through `box_input("checker_eval/masks.npz")`, and `tests/test_train_fast.py` defines `BOX_INPUTS = ("checker_eval/masks.npz",)`.
- Add `tests/test_train_fast.py` to D's S08 module list. Then the $0 box-input check covers it.
- I22b: assert that the step string contains `OBJDEL_DATA_ROOT={data}`.

### 2. major | sources: B (also A, D) | G01 and G08b compare "CPU build and CUDA build" in one process, but a box venv holds only one build

Problem.
- B.9.5 G01: "CPU build and CUDA build give bit-equal results on every field". G08b: "`test_g08b_smoke`: G01 parity on 1,000 windows ... with `sum_x`, `sum_y`".
- B.9.8 sketch `tests/test_scorer_sums.py`: `def test_sums_cpu_gpu_bit_equal(cpu_build, cuda_build)` with `@pytest.mark.gpu`. The fixtures `cpu_build` and `cuda_build` are listed as assumed seams, with no definition.
- A "### A.9.6" G02b states the fact: "One venv holds one build: `build()` returns one value [code: cuda/gpucheck_py.cpp:113], and the build is a CMake option at install [code: CMakeLists.txt:15]. So the CPU side is a Mac golden". `make cuda` reinstalls the package with `OBJDEL_CUDA=ON` in place [code: Makefile:41-45].
- So `cpu_build` cannot exist on the box. In R-D2, D8 step 3 "fails if `skipped` > 0", and D8 says `tests/test_scorer_sums.py` joins the box corpus after B3.1. In V3 item 4, G08b has no reference at all.
- Also, B.9.5 runner table puts G01 in `tests/test_gpu_scoring.py` and `tests/test_cuda_checker_gpu.py`, but the B.9.8 sketch puts G01 in `tests/test_scorer_sums.py`. These are two files for one test.
- `test_g08b_fp16_dtype` ("the `vllm_dtype` engine check") names no model path. The test step gets no argv. The engineer must guess where the checkpoint is (`$OBJDEL_RUN_DIR/ckpt`?).

Fix.
- B.9.5 G01 and G08b: on the box, the reference for `sum_x`, `sum_y` is `np.nonzero` sums of each decoded mask (independent of any build). All other fields are compared with a Mac golden `tests/goldens/b_g08b.json`. `make_a_goldens.py` or a new `make_b_goldens.py` writes it with the CPU build: input SHA-256, then dtype, shape and SHA-256 of `tobytes()` per field (the GPU02 form, D8). A missing golden fails.
- B.9.8: delete the `cpu_build, cuda_build` fixtures and the sketch `test_sums_cpu_gpu_bit_equal`. Replace it with the golden form above.
- Pick one file for G01 (`tests/test_scorer_sums.py`) and name it in the B.9.5 runner row and in D's `D_GPU_TESTS`.
- `test_g08b_fp16_dtype`: the model path is `os.environ["OBJDEL_RUN_DIR"] + "/ckpt"`, the folder that `fetch_hf.sh <run> rl-full-3ep` fills (front 6.3 V3 item 0). The test fails, never skips, when the folder is missing.

### 3. blocking | sources: C | E-SP1 adds a rerun rule that real step times always trip, and it conflicts with G-SP1

Problem.
- C "### C.9.5" E-SP1: "Report median and std; std of steps 3-15 at most 0.3 s in each arm, else "noise, rerun"."
- C "## C.6" G-SP1 states the gate with no such rule: "pass if `fast_mean <= base_mean - 3.4 s` (mode A) ... Anything above: not adopted". Both parts claim to be the decision rule for C adoption.
- Measurement: R5 per-step SD over steps 3-15 is 1.89 s (base arm) and 2.11 s (shared arm) [measured: `runs/r5/base_rl.json`, `runs/r5/shared_rl.json` `step_time`]. The decode tail changes from step to step (38 % of R5 steps hit 1,536 tokens, BRIEF). So every arm fails the 0.3 s rule. E-SP1 then says "rerun" every time.
- A rerun of item 3 is 43 min, $0.49 (front 6.3), and no V2 row or stop rule prices it. An engineer who follows E-SP1 either reruns with no money or cannot close G-SP1.

Fix.
- C.9.5 E-SP1: delete "std of steps 3-15 at most 0.3 s in each arm, else "noise, rerun"". Replace it with: "Report the per-arm SD of steps 3-15 and the Welch 95 % interval of `base_mean - fast_mean` as report lines. They never change the G-SP1 result, and no arm is rerun."
- Add the reference values: R5 SD 1.89 s (base) and 2.11 s (shared). With 13 steps per arm, the SE of the difference is about 0.8 s, so the 3.4 s margin is about 4 SE [estimate: sqrt(1.89^2/13 + 2.11^2/13) = 0.79 s].

### 4. major | sources: C (also 00_front) | Four C.9.5 tests need box arms that front 6.3 does not price

Problem.
- Front 6.3 is "the only place that lists rental items, their order, minutes". V2 item 3 basis: "12 + 8 + 7 + 4 ... + 2 x 6 start-up": one baseline arm, one fast arm (with the check steps). C.8.2 gives the same arms.
- C.9.5 asks for more:
  - E-EQ6: "Control: run the baseline twice first; if it differs from itself report "baseline noise"". That is a second baseline arm, about 18 min.
  - E-REP: "The fast arm twice with the same seed: step times within 5 %". That is a second fast arm. C.9.7 lists E-REP under G-EQ6 and G-SP1.
  - E-MEM: "also with all 256 answers forced to 1,536 tokens". That is another run.
  - E-SP4: "count does not scale per row (B = 128 vs 256 within 10 %)". That is a B = 128 arm.
- The engineer must guess: run these arms (unpriced, above the owner amount; front 6.4 rule 2 cannot stop an item it does not know), or skip them (then G-EQ6 and G-SP1 lose the tests that C.9.7 maps to them).

Fix.
- C.9.5: mark E-EQ6's control, E-REP, the forced-1,536 part of E-MEM and the B = 128 part of E-SP4 as "not in V2; optional arms, each priced in front 6.3 only if the owner adds it". Remove E-REP from the C.9.7 rows for G-EQ6 and G-SP1.
- E-EQ6: without the control arm, a step-1 mismatch is bisected by flag on the box only if the owner funds it; otherwise the result is "C not adopted" (G-SP1 rule).

### 5. major | sources: 00_front (also A) | The 6.1 step 4d launch-order test has no test id, no owner and no spec

Problem.
- 00_front "### 6.1" step 4d: "A grep test in `tests/test_repo_hygiene.py` checks it. ... The same grep test reads the V1-V4 launch lines of 6.2-6.3 in order and checks the rule for each `run.sh` line." Pass column: "test green".
- `tests/test_repo_hygiene.py` is A's file. A "### A.9.4" lists S01-S12 and none of them reads the launch lines or the plan views. A I20 checks only `a3_measure.yaml` and `a3_test.yaml`. B I22a/I22b check single config files. No phase lists the order test in its coverage map.
- The test must parse plan steps and decide "builds the CUDA checker". It must also handle V3 items that are skippable and V4's `specs/e_v4.yaml`, which is a `RunSpec`, not a `steps:` file (front 6.3 "V4 launch config"). None of this is specified. A preflight pass condition that no phase writes cannot be checked before V1.

Fix.
- A.9.4: add S13 "plan views and build order" with this rule. (a) `b1_rows.yaml` and `a3_measure.yaml` name only the `final_eval` views of the 6.2 table and no `--split test`. (b) For each block of launch lines in front 6.2-6.3 (fenced block after "Launch lines"), take the `run.sh` lines in order. For every subset that drops lines marked "only if" or "skippable", the first plan either has a step equal to `b1_rows.yaml` step 0, or has A3's build-if-missing step, or (V4) has `to_argv` of stage 0 equal to `b1_rows.yaml` step 0 (E U51). (c) A planted plan without a build step fails.
- Add S13 to A.9.7 coverage for front 6.1 step 4d, and change 6.1 step 4d to cite "A S13".

### 6. major | sources: D (also A) | S13 fails `make test-data` until a V3 run exists, and it cannot locate the V3 run

Problem.
- D "### D.9.6" S13: "`tests/test_d_records.py`, marker `data`, reads `runs/` ... `test_b04` reads `runs/<V3 run>/b3-v5/rl.json` ... a missing file fails, never skips".
- A "### A.9.6": "A test whose `runs/...` input is absent skips with a reason naming the file." A "## A.6" row "All": "`make test-data` passes on the Mac with `data/`". `make test-data` runs every `data` test (A.9 Makefile line 2265; D8 line "make test-data (`$(PY) -I -m pytest -q -m data`)").
- Under the default (b'), V3 never runs (front 6.3 rule). So from the D7 commit on, `make test-data` is red for good. Any later phase that runs it (B data gates, CLEAN L4, D7 step 4) sees one failure that is not a defect.
- `<V3 run>` is a run id chosen at launch. No rule says how the test finds it (no env var, no index file).

Fix.
- D.9.6 S13: the run folder comes from `OBJDEL_V3_RUN` (the main agent sets it after `pull.sh`). Without that variable, S13 skips with "no V3 run pulled". With it set, a missing file fails. Name the eval folder that `test_b04` reads (the `--out` of B3's `lidc-val-shifted` eval step in `b3_v5.yaml`).
- Front 6.3 V3 note: after `pull.sh`, run `OBJDEL_V3_RUN=<run> make test-data`.

### 7. minor | sources: 00_front | Two arithmetic slips in 6.2

- 6.2 notes: "(2.05 - 0.10) / 0.01133 = 172.3". The value is 172.1 (1.95 / 0.011333 = 172.06). The derived step counts move by under 1 step. Fix: write 172.1 and re-derive "(172.1 - 53) x 60 / 27.5 = 260".
- 6.2: "(2.40 - 0.10) / 0.680 = 3.39 h". The value is 3.38 h. The floor to 0.05 (3.35 h) does not change. Fix: write 3.38.

### 8. minor | sources: A, C | A says no V2 plan step runs pytest; C adds one

- A.9.6 G06: "no V2 plan step runs pytest on the box". C.9.4: "pytest" = `"$PY" -I -m pytest -m gpu tests/test_c_gpu.py`, a step of `c_speed.yaml` after `c_speed_check` in V2 item 2. Fix: in G06, write "the V2 pytest step (C.9.4) runs only `tests/test_c_gpu.py`".

### 9. minor | sources: C | C.9.4 says box tests "never skip", but G-RM1 skips, and the batch path is not named

- C.9.4 intro: "its `gpu` tests run and never skip". Row G-RM1 (sync and time part): "pytest, on the 5 `c_speed_check` batches (skips with a reason if absent)". The test gets no path to the batches. Fix: name the folder (`$OBJDEL_RUN_DIR/c_speed_check/batches/`), and fail when it is missing.

### 10. minor | sources: A, D | Two definitions of `make test-data`

- A (line "new `test-data: $(EXT_STAMP)` runs `$(PY) -I -m pytest -q -m "data and not slow"`", plus `test-all-data`). D8: "`make test-data` (`$(PY) -I -m pytest -q -m data`)". Fix: D8 cites A's target and does not restate it.

### 11. minor | sources: D | The preflight branch for `d_gpu_test.yaml` never runs

- D8 Makefile item: "`experiments/phase2/preflight.py` runs `make gpu-check DRY=1` when the config is `d_gpu_test.yaml`". Front 6.1 step 4: the V3 stamp is `b3_launch.yaml`, and "The plan files ... never follow `preflight CONFIG=`". So the branch is dead code (owner rule: remove dead code). The round-8 log already calls it "Residual". Fix: run `make gpu-check DRY=1` when the config is `b3_launch.yaml`, or delete the branch and keep the D8 Steps rule only.

### 12. minor | sources: B | The G05 record test has no path to `f_check.json`, and its skip makes it vacuous

- B.9.5 G05: "`test_g05_f_check_record` reads `<out>/f_check.json` ... (skips with a reason when the file is absent)". `<out>` is not given to the test. Fix: read `$OBJDEL_RUN_DIR/<out name>/f_check.json` with the out name given as a test constant, and fail when the file is missing on a box (`build() == "cuda"`).

### 13. minor | sources: B, C | `DEFAULT_NUMERICS` has no stated value; the C arms depend on it

- B2.2: "`grpo_config(..., numerics=DEFAULT_NUMERICS)`". U38: "default equals `DEFAULT_NUMERICS`". The value is not written. C0 step 3's baseline argv has no `--numerics`, and C.6 requires "both arms run under `bf16-seqmask`". Fix: B2.2 states `DEFAULT_NUMERICS = "bf16-seqmask"` until B2.4's change. U-CF4 asserts that every C arm's `train_rl` command either has no `--numerics` or has `--numerics bf16-seqmask`.

### 14. minor | sources: A, B, D, E | Test ids repeat across phases without a prefix

- S12 is A's "`reward` never imports `train`" and D's "V3 records are written". I29 is B's sums test and E's `chain()` test. G01 is A's answers test and B's parity test. A-G6 says "S12 passes" with no phase. Fix: in every cross-reference, prefix the test id with its phase (`A S12`, `D S12`), as the gate ids already do.

## Totals

2 blocking (items 1, 3), 4 major (items 2, 4, 5, 6), 8 minor (items 7-14).
