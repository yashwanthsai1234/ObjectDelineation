# Final review round 12: QA and technical writing (Tests sections, gates, testability, cross-section consistency)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice is in scope, so I recommend no catalog skill. It made me check the GPU-only parts of the D-G3 box gate (the CUB scan runs only in the CUDA build), which led to item 1. |
| engineering:code-review | I read `experiments/cuda_checker/compare.py:21, 165-169` (`mask_reference` calls `objdel.checker.polygonize`), `tests/test_cuda_checker.py:13, 35-38`, `tests/test_cuda_checker_gpu.py:46-51`, `src/objdel/checker/__init__.py:6`, `certify.py:18`, `targets/build.py:53, 93`, `experiments/phase2/vast/policy.py:422-536` and `experiments/phase2/train_rl.py:50-68` against the round 11 text. |
| engineering:debug | I traced what each oracle import resolves to after B1.4a, and I ran a probe of Python's submodule binding (item 1). I traced a G-RC0 failure through the `&&` chain of `b3_v5.yaml` (item 2). |
| engineering:system-design | I checked which tree V3 runs under decision 1 (a) against the time line rule 5 (item 3), and the V3 amount rule against its own minutes (item 5). |
| engineering:architecture | One owner per rule: the C10 keep list against the flags that C0 makes required (item 4); E's typed spec against the new `--merge-last` step (item 7). |
| engineering:documentation | I checked the header revision line, the doc size, and which round 11 minor items have no log entry (items 9-13). |
| huggingface-skills:huggingface-best | No model choice is in scope. No change. |
| engineering:testing-strategy | I checked each new round 11 gate for an independent oracle, a run place and a $0 check. Item 1 is a gate whose oracle became the code under test. Item 2 is a box test with no $0 check that gates an unrelated paid run. |
| mattpocock-skills:tdd | I looked for tautological tests (the expected value comes from the code under test). Item 1 is one. U10 and U11 call `signal.ci`/`signal.draws` as the expected value, but they also pin literal intervals, so they are not tautological (I re-computed the literals, below). |

Web tools: loaded (WebSearch, WebFetch). No external claim in my items needed a web check; every proof is in the repo.

Measurements this round ($0, Mac):
- `docs/implementation2.0/Implementation.md` is 1,714,302 bytes, `review-log.md` 363,778 bytes (`ls -l`, 2026-11-09).
  Round 11 measured 1,689,731 bytes.
- Python submodule probe (scratch package `pkg.chk` with `__init__.py` = `from .rings import polygonize`, and files
  `rings.py` and `polygonize.py`, each with a function `polygonize`): `from pkg.chk import polygonize` gives the
  `rings` function; after `from pkg.chk.polygonize import polygonize`, `from pkg.chk import polygonize` gives
  `<module 'pkg.chk.polygonize'>`. The import system sets the parent's attribute to the submodule.
- U10/U11 literals re-computed with A1 step 1's `key` written out (boot 2000, seed 0): found nodule [0.5, 0.62],
  patient [0.525, 0.6]; differences nodule upper 0.1333, patient upper 0.15. They match the doc. Not a finding.
- R5 `grad_norm` pair and `PULL_EXCLUDE` were checked in round 11; not repeated.

## Items

### 1. major: after B1.4a the engine's mask-ring parity tests compare the engine with itself, so D2's box gate D-G3 cannot fail on ring fields (B, D)

- Doc (B.md, B1.4a): "Point the production callers at it: `checker/__init__.py:6` and `checker/certify.py:18`
  import `polygonize` from `.rings`." and "Tests that need the oracle import
  `objdel.checker.polygonize.polygonize` by module path." B1.4a names no oracle caller to change.
- Doc (D.md, D3 Files): only D3 changes "the oracle callers: `tests/test_checker.py`, `tests/test_targets.py` ...,
  `tests/test_cuda_checker.py`, `tests/test_cuda_checker_gpu.py:46`, ..., `experiments/cuda_checker/compare.py`
  (:21, 167), `experiments/cuda_checker/eval.py` (:42, 247), ... `docs/implementation2.0/kernels/K1/cross_numpy.py:16`".
  D3 runs after D2 ("1. Land D2 first").
- Doc (D.md, D2 Files): "change `tests/test_cuda_checker.py`: ... Compare against the numpy oracle (D3 moves it to
  `tests/oracles/`)" and "change `tests/test_cuda_checker_gpu.py`: a 37,093-crop equality test against the oracle".
- Doc (D.md, D.6 D-G3): on the RTX 5090, `pytest -m gpu` on `D_GPU_TESTS` "passes with 0 skipped tests and 0
  mismatches on the R3 corpus (... 6,024 fuzz masks, ... 37,093 LIDC masks". D-G1 (K1 on the GPU) is "R-D1
  (optional; else D-G3 in R-D2)", and R-D1 is dropped on the time line (front 6.3 item 4).
- Code: `experiments/cuda_checker/compare.py:21` is `from objdel.checker import bitquads, polygonize`, and
  `mask_reference` builds the reference rings with `polygonize(mask)` (`compare.py:167`). The mask tests use it:
  `tests/test_cuda_checker.py:38` (`C.mask_diff(out, i, C.mask_reference(...))`) and
  `tests/test_cuda_checker_gpu.py:50`. After B1.4a, `objdel.checker.polygonize` is `rings.polygonize`, which calls
  `_gpucheck.check_masks`. So `ring_sha`, `n_rings`, `b0` and `b1` of the reference come from the engine under test.
  Only the quad fields stay independent (`bitquads`).
- Problem: from B1.4a through D2 (the K1 port into `check_masks`), the existing mask tests and the D2 "against the
  oracle" tests that reuse `C.mask_reference` cannot see a ring bug. On the box this is the only mask-ring check:
  `tests/test_rings.py` is in `D_GPU_TESTS`, but its cases are "CPU build" (D3) and `-m gpu` deselects them. A
  CUDA-only defect (the CUB `DeviceScan` path exists only in the CUDA build, D2 Files) passes D-G3, and R-D2 money
  buys a gate that cannot fail. No planned $0 check flags a test that always passes.
- Second effect (same root): B1.4a keeps a submodule `objdel.checker.polygonize` beside a package attribute
  `polygonize`. Once any test imports the oracle by module path, a later `from objdel.checker import polygonize`
  returns the module (probe above), and `targets/build.py:93` then raises `TypeError: 'module' object is not
  callable`, in an order-dependent way. This half is loud (minor), but the fix is the same.
- Fix (B.md B1.4a; D.md D2, D3; CLEAN11 owner line if it names the move):
  - B1.4a moves the oracle out of the package now: `git mv src/objdel/checker/polygonize.py
    tests/oracles/polygonize_oracle.py` (D3 keeps only the `bitquad` move). `Ring` stays in `rings.py`.
  - B1.4a changes every caller that uses `polygonize` as a reference to import `tests.oracles.polygonize_oracle`:
    `experiments/cuda_checker/compare.py:21, 167` (so `mask_reference` stays numpy), `experiments/cuda_checker/eval.py:42`,
    `docs/implementation2.0/kernels/K1/cross_numpy.py:16`, `experiments/phase1/bench/xarch_hashes.py:45`,
    `tests/test_checker.py`, `tests/test_targets.py` (:14, 28, 32, 314) where they check rings. D3's list keeps only
    the `bitquad` callers.
  - Add to `tests/test_rings.py` a guard: `inspect.getsourcefile(compare.mask_reference)`'s module imports
    `polygonize` from `tests.oracles` (AST check), and `objdel.checker.polygonize` is a function after
    `import tests.oracles.polygonize_oracle`.
  - D2's GPU test names its oracle: "the 37,093-crop equality test compares `check_masks` with
    `tests.oracles.polygonize_oracle.polygonize`, never with `objdel.checker.polygonize`".

### 2. major: G-RC0 is a hard step of `b3_v5.yaml` before the B3 run, but its "tiny model on vLLM" has no defined load path and no $0 check, so a test-harness failure stops the B3 run (B, C, 00_front)

- Doc (B.md, B3.9): "Step 2 is C's G-RC0 box check by node id ...: `[-m, pytest, -q, -p, no:cacheprovider, -m, gpu,
  tests/test_c_gpu.py::test_g_rc0_per_row_caps]`". Step 3 (the B1 re-eval) and the `train_rl` step follow it.
- Doc (C.md, C.9.4): "**G-RC0 per-row caps on the pinned vLLM (rev 3)**. ... Tiny model, vLLM 0.32.0, a
  `SamplingParams` list with caps [4, 8, 16]". C.9.0 defines the tiny model: "tiny Qwen3-VL from
  `tests/test_train_fast.py::tiny_model`".
- Code: `tests/test_train_fast.py:24-35` builds `Qwen3VLForConditionalGeneration(cfg)` in memory, with
  `vocab_size=128`, no tokenizer and no processor files. `vllm.LLM` loads a model from a path or a Hub id; the doc
  names no save step, no tokenizer, and no `skip_tokenizer_init`/token-prompt form.
- Doc (C.md, C11 step 3): the premise is "[web: docs.vllm.ai LLM.generate docstring, latest; unverified for
  0.32.0]". C11's arm applies only "from the first run after Gate B3" (B3.9 last paragraph), so G-RC0 serves a later
  run, not B3.
- Code: a step without a leading `~` joins the chain with `&&` [code: experiments/phase2/vast/policy.py:487-490].
- Problem: the first time this test touches vLLM is on the paid V3 box. If the tiny model cannot be loaded (no
  tokenizer or processor for Qwen3-VL's multimodal profiling), pytest fails, the chain stops before step 3 and the
  `train_rl` step, and V3's required item (B-S2, $2.32 at 27.5 s) gives no Gate B3 result. The $0.54 of item 0 plus
  the build and smoke minutes are spent for nothing. I22b only checks `--collect-only` on the Mac, where the `gpu`
  test skips, so no planned $0 check runs it.
- Fix:
  - B.md B3.9 and I22b: step 2 becomes a soft step, `[~, -m, pytest, ..., tests/test_c_gpu.py::test_g_rc0_per_row_caps]`
    ("a G-RC0 failure is logged; it blocks only C11's arm, never B3"). I22b asserts the leading `~`. Better still,
    move it after the B3 val eval as the last soft step, so it can never delay `train_rl`.
  - C.md C.9.4 (G-RC0 and G-GR1): name the vLLM model. Either the real `{rd}/ckpt` (already on the box; text-only
    prompts, `max_model_len=256`, `gpu_memory_utilization=0.3`; about 2 min [estimate]) or a saved tiny model with
    its tokenizer and processor files in `tests/fixtures/tiny_qwen3vl/`. Add a Mac test that the chosen directory
    loads with `AutoConfig`, `AutoTokenizer` and `AutoProcessor`.
  - 00_front 6.3 item 0: keep the minute; add "soft step: a failure does not stop B3".

### 3. minor: under decision 1 (a) (V3 without V2), the time line does not put C's test file on the V3 tree, but `b3_v5.yaml` step 2 needs it (B, C, 00_front)

- Doc (00_front.md, time line 5): "C1-C5 and B2/B3 code lands on `main` only when a raise funds its gate (V2, V3)".
  Totals (a): "B1 then B3: V1 + B-S2 + reserve" and "The V2 items for C adoption add $1.21".
- Doc (C.md, C.9.0): "C's pytest box tests are all in `tests/test_c_gpu.py`". C11 code is C code
  (`vllm_requests.py`, C6 module).
- Problem: under (a) the V3 tree may hold no `tests/test_c_gpu.py`. I22b's `--collect-only` on the G-RC0 node id
  then fails on the Mac, so it is caught at $0. But the doc does not say whether to bring C's test file and C11 code
  onto the V3 tree or to drop step 2.
- Fix (B3.9; 00_front 6.3 item 0): "Step 2 exists only when C11 code is on the V3 tree (V2 funded and Gate C-adopt
  passed). Otherwise `b3_v5.yaml` has no step 2, and I22b asserts its absence."

### 4. minor: the C10 keep list and U-CL1 omit `--encode-check` and `--checks-from-step`, which C0 makes required (C)

- Doc (C.md, C0 Files): "`--expect-step1` requires `--encode-check FILE`" and the new `--checks-from-step S`.
- Doc (C.md, C10 and U-CL1): keep list "`--allow-cpu-pixels`, `--check-*-every`, `--count-syncs`,
  `--record-batches`, `--expect-step1`, `--group-requests`, `--no-ckpt-layers`, `--runaway-cap`". U-CF3: "After C10
  only the kept flags exist".
- Problem: U-CF3 then fails while `--encode-check` exists, or the engineer deletes it and U-CF5 fails. Both are $0
  checks, so it is minor.
- Fix (C10, U-CL1, U-CF3): add `--encode-check` and `--checks-from-step` to the keep list.

### 5. minor: V3 item 0 rounds 38.2 + 1 down to 39 and leaves out the B1 model download, so the V3 amount rule's `BudgetStop` cuts step 260 at 27.5 s (00_front)

- Doc (00_front.md, 6.3 V3 item 0): "39 | ... launch to first `run` up to 38.2 min [measured] + 1 + b; the rows-root
  download is not in the measured rows". Item 0 also runs `fetch_hf.sh <run> b1-rows --as b1` (a second full merged
  model plus `final_adapter/`), which the measured 38.2 min (one checkpoint) does not hold.
- Doc (00_front.md, V3 amount rule): "step 260 ends at minute 185.2, so `BudgetStop` = $0.10 + 185.2 x $0.01133 =
  $2.20 (rounded up)" and "6.4 rule 5 then gives a `BudgetStop` that does not cut the 260 steps at 27.5 s".
- Arithmetic: (2.20 - 0.10) / 0.01133 = 185.35 min of headroom. With item 0 at 39.2 min, step 260 ends at 185.37
  min, so the stop trips inside step 260 at the planning step time. Any download minute widens the gap by 2.2 steps.
  Gate B3 then compares unequal runs (B3.9: "A lower value cuts the 260 steps").
- Fix (6.3 item 0, V3 amount rule, B-S2 row, totals (a); B.md I22b constant): item 0 = 38.2 + 1 + 1 (B1 model
  download [estimate]) = 40.2, rounded up to 41 min. Recompute: step 1 starts at 68, step 260 ends at 187.2,
  `BudgetStop` $0.10 + 187.2 x $0.01133 = $2.23 (rounded up), amount $2.58. Apply the same 2-min shift to the
  "items 0, 1, 3 and 4" row.

### 6. minor: no test pins that Gate B3 pairs with the V3 B1 re-eval, not V1's B1 eval (B)

- Doc (B.md, B.6.2 Gate B3): "paired against the B1 re-eval of the same rental (B3.9 step 3, `eval_b1v3_val_shifted`
  and `eval_b1v3_val_edge` ...). V1's `eval_b1_val_*` is a report line only".
- Doc (B.md, U49): the B3 part tests only the rule thresholds. No case names the input folders. Both folders exist on
  the Mac after the V3 pull.
- Fix (U49): add "the Gate B3 evaluator, given a run folder with both `eval_b1_val_*` and `eval_b1v3_val_*`, reads
  `eval_b1v3_val_*` for every rule and prints `eval_b1_val_*` as the report line; without `eval_b1v3_val_*` it returns
  'open: B1 re-eval missing', never a pass".

### 7. minor: E's typed spec does not model the new `train_rl --merge-last` step (E)

- Doc (B.md, B1.5 step 3b): "`train_rl --merge-last --model qwen3-vl-2b --ckpt {rd}/ckpt --out {rd}/b1-rows` ...
  `--steps` and `--train` are not required there (argparse checks them only without the flag)".
- Doc (E.md, E6.1): "`FLAG_ORDER[stage.kind]`. For `rl` the tuple starts `model, ckpt, train, out, steps, ...`" and
  "B1.5 writes `b1_rows.yaml` in `FLAG_ORDER`"; E ports `b1_rows.yaml` (golden 5). `grep merge-last` finds no line in
  E.md.
- Problem: an `RlStage` with required `train` and `steps` refuses step 3b, so golden 5 fails at E6 (Mac, $0).
- Fix (E6.1, U49 of E): add `merge_last: bool = False` to `RlStage` (`train` and `steps` optional only when it is
  True, by a model validator), its place in `FLAG_ORDER["rl"]`, and a golden-5 case for step 3b.

### 8. minor: `--merge-last` trusts any existing `final/`, but `train_rl` writes `final/` in two saves before `rl.json` (B)

- Doc (B.md, B1.5 step 3b): "If `{rd}/b1-rows/final/` exists (the run ended normally, or `BudgetStop` ended it), it
  does nothing and exits 0."
- Code: `experiments/phase2/train_rl.py:59-63` saves the model to `final/`, then the processor, then writes
  `rl.json`. A crash or kill between them leaves a `final/` without a processor and no `rl.json`.
- Fix (B1.5 step 3b, I49): "does nothing only when `{rd}/b1-rows/rl.json` exists and `final/` holds the processor
  files; otherwise it re-merges from the last checkpoint into `final/`". I49 adds the case "a `final/` with only
  `model.safetensors` and no `rl.json` is rebuilt".

### 9. minor (round 11 qa-4, no log entry): the G01 sketch still has no `sanitizer` marker (B, D)

- Doc (B.md, B.9.8 sketch): `@pytest.mark.gpu` then `def test_g01_golden():`. D-G5 runs `-m "gpu and sanitizer"`
  and names `tests/test_scorer_sums.py`, so the scorer file contributes no case to D-G5.
- Fix: as round 11 qa-4 (a `test_g07_sanitizer_random` with both markers on the `random` set only).

### 10. minor (round 11 qa-5, no log entry): A6 step 10 still treats `nsys_capture.sh` exit 3 or 4 as a profiler failure (A, C)

- Doc (A.md, A6 step 10): "Fallback, if `nsys_capture.sh` exits 3 or 4: run the same arm command with `--trace`".
  A G-EQ6 stop (`step1_mismatch.json`, exit 3 inside the wrapper) reaches the main agent as exit 4.
- Fix: as round 11 qa-5 ("check `{rd}/c_fast/step1_mismatch.json` first; if it exists, no fallback").

### 11. minor (round 11 qa-6, no log entry): the decision 1 (a) total still uses the B-S2 price, not the V3 amount plus `PULL_S` (00_front, 99_back)

- Doc (00_front.md, totals): "(a) B1 then B3: V1 + B-S2 + reserve = $2.29 + $2.32 + $0.25 = $4.86 ... the V3
  approved amount is $2.55".
- Arithmetic: a raise of $1.86 gives `HARD_STOP` $21.76. At the V3 launch: $16.90 + $2.29 + ($2.55 + $0.11) =
  $21.85 > $21.76, so `launch_refusal` refuses. 6.1 step 4b shows it at $0.
- Fix: (a) = V1 + V3 worst case + reserve = $2.29 + $2.66 + $0.25 = $5.20 (after item 5: re-compute with the new
  amount); the raise is that total minus $3.00. Same form for the other totals and 99_back decision 1.

### 12. minor (round 11 qa-7, no log entry): under (b') a broken R5 hash chain still leaves Gate B1 waiting with no end rule (B)

- Doc (B.md, B.6.2): "If that chain is broken, Gate B1 waits for A3-val in V2 item 1b (front 6.3)." U49: "a broken
  chain returns 'wait for V2 item 1b' (never a pass)". Under (b') V2 never runs.
- Fix (B.6.2, U49): "Under (b'), or when item 1b fails again, Gate B1 takes front 1.1 row 4 with the reason
  'shifted baseline missing'", as the edge bullet already says.

### 13. minor: stale header and doc size (00_front, DOC)

- Doc (00_front.md, line 4): "Revision: final round 3 plus the round-4 consistency pass." This is round 12.
- Doc (DOC.md, DOC1 and DOC-G13): "`Implementation.md` 1,480,096 bytes". Measured 1,714,302 bytes (86 % of 2 MB),
  up 24,571 bytes since round 11.
- Fix: update the revision line and the DOC1/DOC-G13 cells; keep new review text in `review-log.md` only.

## Not repeated

These round 11 fixes are in place, and I checked each against the code or the doc:
- qa-1: `--encode-check` and `{rd}/c_check/encoder.json` (C0 Files, C3 step 0, U-CF5 with the 0.003266/0.003264
  cases, E-EQ6).
- qa-2: G-RC0 as `b3_v5.yaml` step 2 and G-GR1 as `c_speed_c6.yaml` step 0, both by node id; the V2 pytest step
  deselects both. Item 2 above is a new defect in that placement (hard step, no load path), not a repeat.
- qa-3: `c_speed_b.yaml` carries check steps 18-32 with `arm_batches_b`; G-RP2 reads the adopted arm's folder; 3c
  priced at 25 min.
- rl-1: TRL's `_llm_weights_sleeping` lines match the doc [code: .venv trl/generation/vllm_generation.py:379, 503,
  602, 773].
- fresh-arch-2: soft steps exist today (`job()` in `cmd_run`, `policy.py:483-486`); the pull plist has `__REPO__`
  and `StartInterval` 1200.
- fresh-arch-3: B1.4a's code references (`__init__.py:6`, `certify.py:18`, `build.py:53`, `gpucheck_py.cpp:140-158`)
  are correct. Item 1 is a new defect in its oracle rule.
- data-1: U10/U11 literals re-computed (Measurements).

Totals: 0 blocking, 2 major (1, 2), 11 minor (3-13).
