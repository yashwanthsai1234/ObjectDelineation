# Final review round 11: QA and technical writing (Tests sections, gates, testability, cross-section consistency)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice is in scope, so I recommend no catalog skill. It pointed me at the nsys wrapper: `nsys_capture.sh` turns any child failure into `exit 4`, and that hides the new `--expect-step1` exit 3 (item 5). |
| engineering:code-review | I read `experiments/phase2/vast/policy.py:422-508` (`chain`, `cmd_run`, `cmd_pull`, `PULL_EXCLUDE`), `experiments/phase2/train_rl.py:11-67` (`rl.json` holds `log`), `src/objdel/train/configs.py:70-77` and `tests/test_gpu_scoring.py` (host and device scorer compared with `==`) against the test claims. |
| engineering:debug | I traced the V2 arm order in one `&&` chain. A mode A OOM stops the chain before step 18, so the check arm never writes its batches (item 3). I also traced which input the new `--expect-step1` check needs (item 1). |
| engineering:system-design | I matched the V3 amount rule and the decision 1 (a) total against the `launch_refusal` guard (item 6). |
| engineering:architecture | One owner per rule: I checked that every box test has a plan step in a rental (item 2: G-RC0, G-GR1). |
| engineering:documentation | I checked stale text (the header revision line, the doc size) and the 6.1 step 4 config names after E6 (items 9, 10). |
| huggingface-skills:huggingface-best | No model choice is in scope. No change. |
| engineering:testing-strategy | I checked each gate for a defined input, a run place and a pass rule. Items 1, 2 and 3 are gates or stop rules whose input or run place is missing. |
| mattpocock-skills:tdd | I checked that the sketches call real interfaces and carry the markers the runner selects. The G01 sketch has no `sanitizer` marker, which D8's `-m "gpu and sanitizer"` selects (item 4). |

Web tools: loaded (WebSearch, WebFetch). One web check: whether `nsys profile` returns the exit code of the target app. The search found no official statement, so item 5 depends only on the wrapper's `|| exit 4`, which is in the doc.

Measurements this round ($0, Mac):
- `docs/implementation2.0/Implementation.md` is 1,689,731 bytes and `review-log.md` is 353,579 bytes (`ls -l`, 2026-11-09).
- `PULL_EXCLUDE = ("rolling/", "decay/", "stable/", "PULLED_OK", "PULLED_OK.next")` [code: policy.py:80]. So `pull.sh` brings back the `.pt` batch folders that G-RP2 and G-RM1 read. This is not a finding.

## Items

### 1. major: `--expect-step1` must know the G-EN0 result, but no input gives it to `train_rl` (C, 00_front)

- Doc (C.md, C0 files): "flag `--expect-step1 FILE` ... After step 1 it reads step 1 from the TRL log in FILE
  (another arm's `rl.json`) and applies the E-EQ6 comparison, the same function that the gate uses. On a mismatch it
  writes `step1_mismatch.json` into `--out` and exits with code 3 before step 2."
- Doc (C.md, C.9.5 E-EQ6): "`grad_norm` and `sampling/sampling_logp_difference/mean` equal to 4 significant digits
  ... If G-EN0 showed different encode bits, `grad_norm` and the logp difference use relative difference <= 1e-3".
  C.9.9 rule 9 says the same.
- Doc (C.md, C.9.4 G-EN0): "the result (bits equal or not, max abs feature difference) and the noise floor ... are
  written to the check log." No file or key holds the result in a form a program can read.
- Doc (C.md, C.6 G-EQ6): "R5's base and shared arms differ by 6e-4 in `grad_norm`: 0.003266 vs 0.003264
  [measured: runs/r5]". At 4 significant digits, 0.003266 and 0.003264 differ.
- Doc (C.md, U-CF5): it plants only a `rewards/iou_tc/mean` mismatch. No case tests the `grad_norm` branch.
- Problem: the comparison has two rules, and G-EN0 (V2 item 2, `c_speed_check --encoder`, the same `&&` chain)
  picks the rule. The flag gets only FILE. An engineer who uses the strict 4-digit rule kills the fast arm at step 1
  when G-EN0 shows different bits and `grad_norm` moves by about 6e-4. That is the case the 1e-3 rule exists for. Then
  Gate C-adopt fails on a false mismatch after V2's $1.21 is spent. No planned $0 check covers it.
- Fix (C.md C0 files, C3 step 0, C.9.4 G-EN0, U-CF5; 00_front 6.3 V2 item 3):
  - `c_speed_check --encoder` writes `{rd}/c_check/encoder.json` with `{"bits_equal": bool, "max_abs": float,
    "floor_max": float, "floor_p99": float}`.
  - The fast arm passes `--expect-step1 {rd}/c0_base/rl.json --encode-check {rd}/c_check/encoder.json`. A missing
    or unreadable file is an error at start-up, before the model loads. The comparison function takes `bits_equal`
    as an argument.
  - U-CF5 adds two cases: `grad_norm` 0.003266 against 0.003264 with `bits_equal` false goes on to step 2; the same
    values with `bits_equal` true exit with code 3.
  - U-CF4 asserts the new flag on the fast arm.

### 2. major: G-RC0 is priced in V3 item 0, but no V3 plan runs it (C, B, 00_front)

- Doc (00_front.md, 6.3 V3 table, item 0): "... C11's 1-min vLLM check (G-RC0); plus b | 39".
- Doc (C.md, C11 step 3): "vLLM check (box, rental V3 item 0, front 6.3; 1 min, before the B3 run ...".
- Doc (C.md, C.9.4): "| G-RC0 | V3, before the B3 run | pytest `-k rc0`, tiny model on the pinned vLLM (C11 step 3) |".
  The same section defines "pytest" only for V2: "`"$PY" -I -m pytest -m gpu tests/test_c_gpu.py`, a step of
  `c_speed.yaml`".
- Doc (00_front.md, 6.3 V3 launch lines): `d_gpu_test.yaml` (item 1), `b2_fp16.yaml` (item 2b), `a3_test.yaml`
  (item 3), `b3_v5.yaml` (item 4). Item 0 is launch, sync and fetch only, with no `run.sh` line.
- Doc (B.md, B3.9 and I22b): `b3_v5.yaml` = step 0 (CUDA build), step 1 (G08b smoke), the `train_rl` step and one
  `final_eval` step. No G-RC0 step. `grep -n rc0` over all phase sources finds only C.md:1845.
- Code: `chain` runs a `steps:` plan exactly as written [code: experiments/phase2/vast/policy.py:432-433]. So a test
  that no plan names never runs on the box.
- Rule: B.9.5 says "`pytest -m gpu` without a file list is never run on a box". `pytest -k rc0` names no file.
- Problem: the one box check that vLLM 0.32.0 honours a per-row `SamplingParams` list never runs. The C11 arm in the
  first run after Gate B3 then spends money on an unchecked premise. The C11 text marks that premise "unverified for
  0.32.0". G-GR1 has the same gap: "V2 optional C6 arm | pytest `-k gr1`", but C0 files and C6 do not list it as a
  step of `c_speed_c6.yaml`.
- Fix:
  - B.md B3.9 and I22b: `b3_v5.yaml` step 2 (after G08b, before `train_rl`) is `[-m, pytest, -q, -p,
    no:cacheprovider, -m, gpu, tests/test_c_gpu.py::test_g_rc0_per_row_caps]`. I22b asserts the step, and
    `pytest --collect-only` on the node id lists exactly 1 item.
  - C.md C.9.4: the G-RC0 row names that node id and "`b3_v5.yaml` step 2". The G-GR1 row names
    `tests/test_c_gpu.py::test_g_gr1_requests_independent` as the first step of `c_speed_c6.yaml` (C0 files).
    U-CF4 checks both.
  - 00_front 6.3 item 0: "G-RC0 runs as `b3_v5.yaml` step 2 (its minute stays in item 0)".

### 3. major: if mode A goes out of memory before step 18, no check arm runs, and mode B (3c) cannot pass Gate C-adopt (C, 00_front)

- Doc (C.md, C.6): "Check arm = steps 18-32 of the fast arm's process (`--steps 40 --checks-from-step 18 ...
  --record-batches {rd}/c_check/arm_batches`)". "Fast arm = C1-C5, mode A (vLLM awake) unless G-MEM fails, then
  mode B. Mode B is V2 item 3c, its own plan `c_speed_b.yaml`".
- Doc (00_front.md, 6.3 V2 item 3c): "Optional mode B arm, `c_speed_b.yaml` ..., only if mode A fails C-G-MEM (17
  steps read, `--steps 20`)". So it has no check steps. U-CF4 asks for `--expect-step1` and `--record-batches` on
  "the fast arm" only, which is the mode A arm.
- Doc (C.md, C.6 G-EQ2, G-EQ3): ">= 50 batches, of which >= 15 come from the box". Doc (C.md, C.9.4): "with it,
  fewer than 15 box batches fails". The 5 `c_speed_check` batches alone are fewer than 15.
- Doc (C.md, C5 step 2): mode A "vLLM about 11.2 GB resident beside a trainer peak of about 17 GB = about 28 GB of
  32 GB [estimate]". So an out-of-memory crash is the likely way for mode A to fail G-MEM.
- Code: one plan is one `&&` job [code: policy.py:487-494]. A crash of the fast arm process ends `c_speed.yaml`
  before step 18, so `c_check/arm_batches` is empty.
- Problem: in exactly the case that funds item 3c, G-EQ2 and G-EQ3 cannot reach 15 box batches. Gate C-adopt then
  fails whatever mode B shows, and the $0.18 of item 3c is spent for nothing. The engineer must guess whether to add
  the check steps to `c_speed_b.yaml`, and that guess changes whether C can be adopted.
- Fix (C.md C0 files, C.6, C.9.4, U-CF4; 00_front 6.3 item 3c):
  - `c_speed_b.yaml`'s arm runs `--steps 40 --checks-from-step 18 --check-rope-every 1 --check-meta-every 1
    --record-batches {rd}/c_check/arm_batches_b --expect-step1 {rd}/c0_base/rl.json` (plus item 1's
    `--encode-check`).
  - G-RP2 and G-RM1 read `arm_batches` or `arm_batches_b`, whichever holds 15 files.
  - Item 3c minutes: 8 x 40/17 + 6 = about 25 min = $0.28 [estimate: C.8.2 rate], or keep 16 min when the mode A
    arm reached step 32. U-CF4 asserts the flags on both fast arms.

### 4. minor: the G01 sketch has no `sanitizer` marker, and the G01 inputs may pass D-G5's 300 s limit (B, D)

- Doc (D.md, D8): "its `tests/test_scorer_sums.py` joins this corpus; its GPU cases carry `gpu` and `sanitizer`".
  D-G5 runs `-m "gpu and sanitizer"` "each tool inside `timeout 300`".
- Doc (B.md, B.9.8 sketch): `test_g01_golden` has `@pytest.mark.gpu` only. If the engineer follows the sketch,
  `-m "gpu and sanitizer"` deselects it, and G07 passes with no scorer case.
- Doc (B.md, G01): "`b_sets` ... returns the five sets ... in 100.9 s" on the Mac CPU build, and
  "`ExpertCache.build` ... the builds take 8.6 s, 0.6 s and 90.6 s". Under memcheck or racecheck the 10,203 table
  uploads and score calls also run slower. R10 qa-10 raised the size limit as minor, and no log marks it handled.
  The 100.9 s figure is new evidence.
- Fix (B.md G07 and the sketch; D.md D8): add `test_g07_sanitizer_random` with `@pytest.mark.gpu` and
  `@pytest.mark.sanitizer`. It runs the `random` set only (203 masks, `random_masks(np.random.default_rng(1))`) with
  the `np.nonzero` sums check and the golden `random` entry. `test_g01_golden` keeps `gpu` only.

### 5. minor: `nsys_capture.sh` maps the step-1 stop to `exit 4`, which A6 step 10 treats as a profiler failure (A, C)

- Doc (A.md, A6 step 8): `"$NSYS" profile ... "$@" --profile-window "$S" "$E" || exit 4`.
- Doc (A.md, A6 step 10): "Fallback, if `nsys_capture.sh` exits 3 or 4: run the same arm command with `--trace ...`".
- Doc (C.md, C0 files): on a mismatch `train_rl` "exits with code 3". The fast arm runs inside `nsys_capture.sh`
  (C0 step 3).
- Problem: the main agent sees `exit 4` for a G-EQ6 stop and may rerun the rejected arm as the A6 fallback (about 6
  min of start-up plus 1 step [estimate]).
- Fix (C.md C0 files; A.md A6 step 10): "Before the A6 fallback, check `{rd}/c_fast/step1_mismatch.json`. If it
  exists, the arm stopped on G-EQ6: no fallback."

### 6. minor: the decision 1 (a) total uses the B-S2 price, not the V3 approved amount, so the guard can refuse V3 (00_front, 99_back)

- Doc (00_front.md, 6.3 totals): "(a) B1 then B3: V1 + B-S2 + reserve = $2.29 + $2.17 + $0.25 = $4.71 ... the V3
  approved amount is $2.40". "V3 amount rule ... Worst case = amount + $0.11". "Any raise also raises `HARD_STOP` by
  the same amount".
- Arithmetic at planning values: a raise to $4.71 gives `HARD_STOP` = $16.90 + $4.71 = $21.61. At the V3 launch,
  the committed spend is $16.90 + $2.29 = $19.19, and the worst case to the deadline is $2.49 to $2.51, which gives
  $21.68 to $21.70 > $21.61. So `launch_refusal` refuses V3. V1 itself uses the stricter form: "$2.49 + $0.25 = $2.74
  <= $3.00".
- The 6.1 step 4b check would show the refusal at $0, so this is minor (caught before spend). The owner would still
  need a second approval message.
- Fix (00_front 6.3 totals; 99_back decision 1): state (a) as V1 + V3 worst case + reserve = $2.29 + $2.51 + $0.25
  = $5.05 (with items 1 and 3: $2.29 + $3.19 + $0.25 = $5.73). The raise is that total minus $3.00.

### 7. minor: under (b'), a broken R5 hash chain leaves Gate B1 waiting with no end rule (B)

- Doc (B.md, B.6.2 Paired baseline): "If that chain is broken, Gate B1 waits for A3-val in V2 item 1b (front 6.3)."
- Under (b') V2 never runs. The edge bullet has an end rule ("Under (b') ... Gate B1 takes front 1.1 row 4"), but
  the shifted bullet has none.
- Fix: add "Under (b'), or when item 1b fails again, Gate B1 takes front 1.1 row 4 with the reason 'shifted baseline
  missing'." Add the case to U49.

### 8. minor: G-RC0 and the C11 arm depend on the `after-raise` merge, but the V3 item 0 text does not say so (00_front)

- Doc (00_front.md, 1 "End state under (b')"): C1-C11 code waits on `after-raise`. V3 item 0 lists G-RC0 with no
  tree condition.
- Fix: with item 2, add "(on the merged `after-raise` tree; C11 code)" to the item 0 cell.

### 9. minor: front 6.1 step 4 names only the pre-E6 launch configs (00_front)

- Doc (00_front.md, 6.1 step 4): "Before V2: `CONFIG=experiments/phase2/configs/c_launch.yaml`; before V3:
  `CONFIG=experiments/phase2/configs/b3_launch.yaml`". Front 6.3 "After E6" says these become `specs/c_speed.yaml`
  and `specs/b3_v5.yaml`.
- Fix: add "(after E6: `specs/c_speed.yaml`, `specs/b3_v5.yaml`; front 6.3 'After E6')".

### 10. minor: stale header, and the doc keeps growing toward its 2 MB limit (00_front, DOC)

- Doc (00_front.md, line 4): "Revision: final round 3 plus the round-4 consistency pass." This is round 11. R10
  qa-11 raised it, and the line is unchanged.
- Measured: 1,689,731 bytes (84 % of DOC-G13's 2 MB). It was 1,665,519 bytes in round 10.
- Fix: update the revision line and the DOC1 baseline cell. Keep new review text in `review-log.md` only.

## Not repeated

These round 10 fixes are in place, and I checked each against the code:
- qa-1: G01 tables rule, `b_tables`, the 4-argument call. Host and device scorers compare with `==` in
  `tests/test_gpu_scoring.py:124, 139`, so a bit-level golden is sound.
- qa-2: three V2 plans and `--expect-step1`. `rl.json` holds `log` [code: train_rl.py:63-67].
- qa-3: C9 has one home.
- qa-4: one B-S2 definition and the V3 amount rule. I re-did the arithmetic: 39 + 14 = 53, then 172.2 min gives
  $2.05 / $2.40; 112.6 then 231.8 min gives $2.73 / $3.08.
- qa-5: `OBJDEL_V2_RUN` and `arm_batches`.

Items 1 and 3 are new defects in the round 10 qa-2 design, not repeats.

Totals: 0 blocking, 3 major (1-3), 7 minor (4-10).
