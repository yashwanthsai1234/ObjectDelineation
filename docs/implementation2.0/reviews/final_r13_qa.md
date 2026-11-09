# Final review round 13: QA and technical writing (Tests sections, gates, testability, cross-section consistency)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice is in scope, so I recommend no catalog skill. It made me check how the CPU-build test variants of the CUDA engine (`K1_ORDER`, `GC_POISON_SCAN`, `K1_MUTATE_FG4`) get loaded, which led to item 1. |
| engineering:code-review | I read `src/objdel/checker/__init__.py:6`, `checker/certify.py:18`, `targets/build.py:53`, `experiments/cuda_checker/compare.py:21, 167`, `eval.py:42, 247`, `xarch_hashes.py:45, 255`, `K1/cross_numpy.py:16, 98`, `tests/test_checker.py:4, 76, 120, 132`, `tests/test_targets.py:14, 28, 32, 314` against B1.4a (all match), `experiments/phase2/vast/policy.py:374-383, 420-492, 527-536`, `cuda/gpucheck.cpp:136-141`, `K1/k1poly.cpp:138-146`, and the editable finder `.venv/lib/python3.12/site-packages/_editable_skbc_objdel.py`. |
| engineering:debug | I ran a probe of what a child Python imports as `objdel._gpucheck` when a variant directory is first on `PYTHONPATH` or `sys.path` (item 1). I traced the soft G-RC0 step through `policy.cmd_run.job` (`~` -> `{ cmd || echo ...; }`). |
| engineering:system-design | I re-did the V3 minute and dollar arithmetic of front 6.3 (item 0 39, item 4 156.6, B-S2 195.6 = $2.32, step 260 at minute 185.2, `BudgetStop` $2.20) and the C.9.4 tail sum (item 5). |
| engineering:architecture | One owner per rule: who checks a frozen val view (B1.6 and I50 against E05, item 2); where the exclusion list is read on the box (item 3). |
| engineering:documentation | Header revision line, stale "moved by D3" lines after the round 12 move (items 7, 8). |
| huggingface-skills:huggingface-best | No model choice is in scope. I checked only that `runs/rl/train_rl/final` (the local copy of `rl-full-3ep`) holds `config.json`, `generation_config.json`, `tokenizer*.json`, `processor_config.json`, `chat_template.jinja` and `model.safetensors` (U-RC6's file list loads). |
| engineering:testing-strategy | For each new round 12 check I asked: does it run in a named gate, and can it fail? Item 1 is a check that cannot fail as specified and that no gate runs. |
| mattpocock-skills:tdd | Red before green: a poison or mutant run must first show red on a planted defect. Item 1's runs stay green whatever the variant does, because the variant is never loaded. |

Web tools: loaded (WebSearch, WebFetch). One web check (item 4: vLLM multimodal profiling with a small `max_model_len`).

Measurements this round ($0, Mac, 2026-11-09):
- `docs/implementation2.0/Implementation.md` 1,731,404 bytes, `review-log.md` 370,597 bytes (`ls -l`).
- Train store items per excluded series (`python -I`, `data/store/lidc/train.parquet`): 380786 -> 15, 228428 -> 5,
  179144 -> 4, 295628 -> 3 (27), 7,341 items. Matches B1.7.
- Variant-load probe (scratch dir `r13qa/mut/objdel/{__init__.py,_gpucheck.py}`, run from `/`):
  `PYTHONPATH=mut .venv/bin/python -I -c "import objdel._gpucheck as g; print(g.__file__)"`,
  the same without `-I`, and `-I` with `sys.path.insert(0, mut)` all print
  `.venv/lib/python3.12/site-packages/objdel/_gpucheck.abi3.so`. The variant is never imported.
- V3 arithmetic: 3 + 5 + 13.0 + 6 + 260 x 27.5 / 60 + 10.4 = 156.57 min; 39 + 156.6 = 195.6 min; $0.10 + 195.6 x
  $0.01133 = $2.316. Step 260 ends at 39 + 14 + 13 + 119.17 = 185.2 min; $0.10 + 185.2 x $0.01133 = $2.198. All match.

## Items

### 1. major: the CPU-build test variants (`GC_POISON_SCAN`, `K1_ORDER`, `K1_MUTATE_FG4`) cannot be loaded by the plan's mechanism, and `make test-order` is never defined or gated, so the round 12 scan-aliasing check passes without testing anything (D)

- Doc (D.md, D2 step 7): "Keep the CPU-build test hooks as compile definitions of a test-only CPU build target:
  `K1_ORDER` 0/1/2 (... `make test-order` runs `tests/test_cuda_checker.py -k masks` three times), `GC_POISON_SCAN`
  (... `make test-order` runs `tests/test_cuda_checker.py -k masks` a fourth time under it, so a `scan` call given
  `L.dev` or `L.score` fails on the Mac at $0, step 1)".
- Doc (D.md, D2 step 1): "the CPU `scan` ignores the block ... so no plain CPU test sees it (step 7 adds the poison
  build that does)".
- Doc (D.md, D.9.9 "Open test-support needs" item 1): "the mutant is built into its own directory and the test runs
  a subprocess with that directory first on `PYTHONPATH`". Sketch S3 (D.9.10):
  `env = dict(os.environ, PYTHONPATH=path + ...)` and `subprocess.run([sys.executable, "-I", "-c", code], env=env, ...)`.
  This is the only load mechanism the doc gives for a variant build.
- Doc (D.md, D.6 D-G3): the Mac half lists `tests/test_cuda_checker.py`, `test_rings.py`, ... "pass on the Mac CPU
  build". No gate, no D8 Makefile change (D8 adds only `gpu-check` and `old-engine-tar`) and no CI job names
  `make test-order`; the Makefile has no such target today [code: Makefile targets setup, build, test, test-all,
  cov, preflight, lint, fmt, wheels, cuda].
- Code: the editable install puts its finder first: `sys.meta_path.insert(0, ScikitBuildInplaceFinder(...))`, with
  the static map `{'objdel._gpucheck': 'objdel/_gpucheck.abi3.so'}` [code:
  .venv/lib/python3.12/site-packages/_editable_skbc_objdel.py, `install(...)` at the end of the file;
  `find_spec`: "if fullname in self.known_wheel_files: ... return self._make_spec(...)"]. And `-I` ignores
  `PYTHONPATH` (CLEAN1 step 5 (a) says so too).
- Measurement: the probe above. With the variant first on `PYTHONPATH` (with or without `-I`) or first on
  `sys.path`, `objdel._gpucheck` is the installed `site-packages/objdel/_gpucheck.abi3.so`. CLEAN1 step 5 (b)
  measured the same for the package: "with `-I` only, 33 of 44 first-party modules came from the repo".
- Result: an engineer who builds the `GC_POISON_SCAN` variant and loads it the way D.9.9 says gets the normal build.
  The fourth `-k masks` run then passes whether or not D2 maps `tmp` to a live block. So the check that round 12
  (cuda-1) added to catch scan aliasing "on the Mac at $0" cannot fail. An aliased map passes every planned $0 check,
  and shows only in D-G3 on the RTX 5090 (R-D2 = V3 item 1, $0.31-0.52, "Any failure: stop; bisect on the Mac"), or
  never when the owner does not name item 1. `K1_ORDER` runs pass the same way. (The `K1_MUTATE_FG4` test I10/D-G16
  fails red at $0, because it asserts `False` and gets the normal `True`; that part alone would be minor.)
- Fix (D.md D2 Files and step 7, D.9.9 item 1, sketch S3, D.6 D-G3 and D-G16, D8 Makefile list; 00_front 6.1 step 3
  or step 1):
  1. One build mechanism for every CPU test variant: a separate venv per variant under `build/variants/<name>/`
     (`uv venv`, then a NON-editable `uv pip install --no-deps . -Ccmake.define.<DEFINE>=ON
     -Cbuild-dir=build/variants/<name>/cmake`), so no editable finder exists in that venv. The variant tests run
     with that venv's python (`build/variants/<name>/bin/python -I -m pytest ...`), never through `PYTHONPATH`.
  2. A load guard that turns a wrong load red: D2 adds a test-only binding `_gpucheck.test_hooks() -> tuple[str, ...]`
     (the compile definitions present, empty in shipped builds). Every variant run starts with
     `assert set(g.test_hooks()) == {"GC_POISON_SCAN"}` (or `K1_ORDER=1`, `K1_MUTATE_FG4`), and the normal suite
     asserts `g.test_hooks() == ()`.
  3. Red first: a planted aliasing (a test-only define `GC_ALIAS_SCAN_FOR_TEST` that passes `L.dev` to `scan` in
     one call) must make the poison run fail once, on the Mac, before D2 merges (record the red run in D2's note).
  4. D8 adds the Makefile target `test-order` (the three `K1_ORDER` runs, the `GC_POISON_SCAN` run, and I10/D-G16 on
     the `K1_MUTATE_FG4` venv). D-G3's Mac half and D-G16 name `make test-order`; D-G12 / the CI job runs it.
  5. Sketch S3: replace the `PYTHONPATH` env with `[str(mutant_python), "-I", "-c", code]`, where the fixture
     `mutant_python` is `build/variants/k1_mutate_fg4/bin/python`.

### 2. minor: E05 says final_eval checks only the test view, while B1.6 and I50 make it check the frozen val views (B)

- Doc (B.md, B1.6): "`experiments/phase2/final_eval.py`: for a frozen view, load `frozen_views.json`, recompute the
  content hash and the target-set hash, and stop with a clear message on a mismatch". The frozen views are
  `lidc-test-shifted`, `lidc-val-edge`, `lidc-val-shifted`.
- Doc (B.md, I50): "gets an entry with `content_sha256` from `frozen_entry("lidc-val-edge")`. A planted copy without
  the file gives `None`, and `final_eval` refuses (E05)."
- Doc (B.md, E05): "missing `frozen_views.json` gives a clear refusal ...; a non-test view is not checked."
- Problem: two $0 tests that contradict each other on `lidc-val-edge` (E05 says no check, I50 says refuse). One of
  them fails on the Mac, so it is caught at $0. `final_eval.json` always carries `view_content_sha256` (A.md run
  keys), so Gate B1's hash rule works either way.
- Fix (E05): "a frozen view (test shifted, val edge, val shifted) is checked; a view that is not a key of
  `frozen_views.json` (val rescaled, val rotated, probes) is not checked; a missing file refuses only for a view
  name in `FROZEN_VIEWS` (a constant beside `frozen_entry`)". I50: name `FROZEN_VIEWS` as the source of the view set.

### 3. minor: B1.3 and DG1b describe a "box preflight" that does not exist (B)

- Doc (B.md, B1.3 `rows_leakage`): "Called by `build_rl_rows` ... and again by `preflight` on the downloaded root.
  The code sync drops `docs/` ..., so on the box `preflight` passes the manifest's `schedule.excluded_ids`".
  DG1b row (B.6.1): "no `source_item_id` in `train_exclusions.json` (box: the manifest's `excluded_ids`)".
- Code: preflight runs only on the Mac: "Preflight: everything that can be proven on the Mac, before any GPU is
  rented" [code: experiments/phase2/preflight.py:1]; B1.3 change 2 checks the local root `data/views/<name>`.
- Fix: B1.3 "again by `preflight` on the local root that B1.4 uploaded (Mac; the file is present there)". Remove
  "on the box" and "(box: the manifest's `excluded_ids`)", or name the box caller if one is meant (none is planned).

### 4. minor: G-RC0 and G-GR1 load Qwen3-VL in vLLM with `max_model_len=256` and no `limit_mm_per_prompt` (C)

- Doc (C.md, C.9.4 "Model of the two vLLM box tests"): "`vllm.LLM(model=..., max_model_len=256,
  gpu_memory_utilization=0.3, enforce_eager=True)`, text-only prompts".
- Web: vLLM warns that a short profiling length cannot hold the worst-case multimodal embeddings and suggests a
  larger `max_model_len` or fewer multimodal items; for text-only use of a Qwen-VL model, the vLLM forum answer is
  `limit_mm_per_prompt={"image": 0}` [web: https://github.com/vllm-project/vllm/issues/20123;
  https://discuss.vllm.ai/t/problem-with-qwen2-5-vl-7b/1065]. Not verified for vLLM 0.32.0.
- Impact: G-RC0 is a soft last step; G-GR1 is step 0 of an optional V2 arm. A start-up error costs about a minute.
- Fix (C.9.4): add `limit_mm_per_prompt={"image": 0}` to the `vllm.LLM` call of both tests (text-only prompts need no
  vision profiling); U-RC6 stays as is.

### 5. minor: the C.9.4 tail sum leaves out the B3 val eval's process imports (C, 00_front)

- Doc (C.md, C.9.4): "tail allowance (B3 val eval 10.4 + G-RC0 2 + cap rounding 3 = 15.4 min < 30.9 min)".
- Doc (00_front 6.2 notes): each eval process adds untimed imports, "about 5.2 min for ... the two eval processes"
  (2.6 min each); 6.3 item 4 adds "2.6 imports" to the B1 re-eval but not to the B3 val eval's 10.4.
- With 2.6 min: 18.0 min < 30.9 min. The rule still holds. Fix: write "B3 val eval 10.4 + its imports 2.6 + G-RC0 2 +
  cap rounding 3 = 18.0 min < 30.9 min".

### 6. minor (carry-over, round 12 qa-3, still open): under decision 1 (a) the V3 tree may not hold `tests/test_c_gpu.py` (B, C, 00_front)

- Doc (00_front, time line 5 and "End state under (b')"): C1-C11 code and C0's V2 tools wait on `after-raise`; order
  rule 3: "B2 and B3 use C code only after C-G-EQ1..6 ... and C-G-SP1 pass". Yet 6.1 step 4 needs U-RC6
  (`tests/test_c_gpu.py::test_vllm_check_model_loads`) to pass before V3, and `b3_v5.yaml` step 5 runs
  `tests/test_c_gpu.py::test_g_rc0_per_row_caps`.
- Caught at $0 (I22b collect-only, U-RC6 pass rule). Fix (B3.9, C.9.0): "`tests/test_c_gpu.py` (the two vLLM tests
  and U-RC6 only) lands with B3.9 on the V3 tree whatever decision 1 picks; it imports no C1-C11 module."

### 7. minor (carry-over, round 12 qa-4, still open): the C10 keep list and U-CL1 omit `--encode-check` and `--checks-from-step` (C)

- Doc (C.md, C10): "Keep `--allow-cpu-pixels`, `--check-*-every`, `--count-syncs`, `--record-batches`,
  `--expect-step1`, `--group-requests`, `--no-ckpt-layers`, `--runaway-cap`". U-CF3: "After C10 only the kept flags
  exist". C0 makes `--encode-check` required with `--expect-step1` (U-CF5).
- Fix: add `--encode-check` and `--checks-from-step` to the C10 keep list and U-CL1.

### 8. minor: stale "D3 moves the polygonizer" lines after the round 12 move to B1.4a (E, CLEAN)

- E.md E3 move table: "E-G16 lists them as "moved by D3""; E-G16 row: "minus ... D3's moves to `tests/oracles/`";
  E.5: "`tests/oracles/polygonize_oracle.py` after D3". CLEAN.6 CLEAN-G14 fail column: "numpy `polygonize.py` in
  `src/` (D3 not landed)".
- Fix: "moved by B1.4a (`polygonize.py`) and D3 (`bitquad.py`)"; CLEAN-G14 fail: "`polygonize.py` in `src/` (B1.4a
  not landed) or `bitquad.py` in `src/` (D3 not landed)".

### 9. minor: header revision line is stale (00_front)

- Doc (00_front, line 4): "Revision: final round 3 plus the round-4 consistency pass." The doc is at round 13.
- Fix: "Revision: final round 13. The review logs are in `docs/implementation2.0/review-log.md`."

## Checked, no finding

- B1.4a's line references and the oracle move: every cited line matches the code; `tests.oracles` imports under `-I`
  through the repo-root `.pth`, and the code sync ships untracked `tests/oracles/*` (`git ls-files -co`).
- I50 (c): `experiments` and `experiments/phase2` are regular packages, so a child with the copy first on `sys.path`
  imports `final_eval` from the copy, and `parents[2]` is the copy root. The test can pass and can fail.
- The soft G-RC0 step: `job()` wraps a `~` step as `{ cmd || echo ...; }`, so a failure leaves `EXIT` 0 and B3's
  outputs intact [code: policy.py:484-487].
- B.9 and C.9 test ids: every id that the coverage maps cite is defined (script check; the only "undefined" names
  are gate ids such as D-G14 and the inline U-SP1).
