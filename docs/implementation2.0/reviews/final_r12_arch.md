# Final review round 12: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc reviewed: the built
`docs/implementation2.0/Implementation.md` (16,750 lines [measured: `wc -l`, 2026-11-09]). Findings are written against
its sources in `impl/phases/`. I read `review-log.md` (rounds 9-11 for front, E, DOC, CLEAN, back), the round-11
reviews `final_r11_arch.md` and `final_r11_fresh-arch_verdict.md`. This round focuses on the text that round 11
added (B1.4a polygonizer seam, `train_rl --merge-last`, the pull job) and how it meets E, DOC and CLEAN. No item
repeats an item that a log marks as handled. The round-11 arch minors (items 1-10) are still open in the sources; I do
not re-file them. I edited no repo file and ran no git command that changes anything. Probe files are in
`reviews/r12arch_probe/`.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA product claim. No catalog lookup, no install. |
| engineering:code-review | I read B1.4a against the real importers of `polygonize` (`grep -rn polygonize src experiments tests`): `checker/__init__.py:6`, `certify.py:18`, `targets/build.py:53`, `compare.py:21,167`, `eval.py:42`, `xarch_hashes.py:45`, `K1/cross_numpy.py:16`, `tests/test_targets.py:14`, `tests/test_checker.py:4` (items 1, 2). |
| engineering:debug | I reproduced the package-attribute shadowing in a scratch package (item 1, `r12arch_probe/probe.py`): after a submodule import by path, `from pkg import polygonize` returns a module and the call raises `TypeError`. |
| engineering:system-design | I walked the order B1.4a -> D2 -> D3 and asked which parity check is real in each window (item 2). |
| engineering:architecture | Seam check: B1.4a keeps a numpy module named `polygonize` inside the package whose `__init__` re-exports a function of the same name. E's own R-1 names this hazard (`eval/__init__.py` re-exports `metrics`). |
| engineering:documentation | DOC9 owner table vs B1.4a's production switch (item 5); stale "D3 puts the polygonizer in production" text (item 4). |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |
| system-design-skills:system-design | Same order walk, framed as "which gate first sees a K1 bug". |
| system-design-skills:service-decomposition | One package with import-linter seams (a modular monolith). Right form; no split proposed. |
| mattpocock-skills:codebase-design | Deletion test on the numpy `checker/polygonize.py` after B1.4a: no `src/` caller is left, so it is a test oracle and belongs outside the package now, not at D3 (items 1, 2). Interface check of `RlStage` vs the new `--merge-last` mode (item 3). |
| WebFetch | Python reference, "The import system", section 5.4.2 "Submodules": "When a submodule is loaded using any mechanism ... a binding is placed in the parent module's namespace to the submodule object" [web: https://docs.python.org/3/reference/import.html]. The probe confirms it on the Mac Python (item 1). |

## Items

### 1. minor: B1.4a keeps the numpy module `checker/polygonize.py` beside the re-exported function `checker.polygonize`; importing the oracle by module path rebinds the package name to the module.

Sources: B (B1.4a), D (D3).

- B.md "#### B1.4a": "Point the production callers at it: `checker/__init__.py:6` and `checker/certify.py:18` import
  `polygonize` from `.rings`." and "The numpy `checker/polygonize.py` stays only as the test oracle. Tests that need the
  oracle import `objdel.checker.polygonize.polygonize` by module path. D3 moves it to `tests/oracles/` (D3 step 4)."
- Python binds an imported submodule as an attribute of its package. So once a test runs
  `from objdel.checker.polygonize import polygonize`, the name `objdel.checker.polygonize` is the module, not the
  function. Every later `from objdel.checker import polygonize` gets a module [measured: `reviews/r12arch_probe/probe.py`:
  "before oracle import: function engine", "after oracle import: module", "TypeError: 'module' object is not
  callable"].
- In one pytest process `tests/test_rings.py` (B1.4a, imports the oracle by path) is collected before
  `tests/test_targets.py`, whose line 14 is `from objdel.checker import instance_topology, min_area, polygonize, ...` and
  whose lines 28, 32, 314 call `polygonize(...)` [code]. Those tests then fail. `targets/build.py:53` has the same
  import form, so it breaks whenever it is first imported after the oracle.
- E.md R-1 already names this hazard for `eval/__init__.py` ("re-exports the function `metrics`, which hides the
  submodule of the same name").
- Severity: `make test` fails loudly on the Mac at $0 before the `v1-b1` tag, so minor. But the B1.4a text leads to it.

Fix (B.md B1.4a, D.md D3, CLEAN.md CLEAN11 owner paragraph): B1.4a does D3's oracle move for the polygonizer now:
`git mv src/objdel/checker/polygonize.py tests/oracles/polygonize_oracle.py`. The oracle imports `Ring` from
`objdel.checker.rings`. B1.4a re-points the polygonize oracle callers in D3's list (`tests/test_checker.py`,
`tests/test_targets.py`, `tests/test_cuda_checker.py`, `tests/test_cuda_checker_gpu.py:46`,
`experiments/cuda_checker/compare.py:21,167`, `experiments/cuda_checker/eval.py:42,247`,
`experiments/phase1/bench/xarch_hashes.py:45`, `docs/implementation2.0/kernels/K1/cross_numpy.py:16`) to
`tests.oracles.polygonize_oracle`. D3 then moves only `bitquad.py`, because `certify` still uses bit-quads until D3's
leg 1. Replace "Tests that need the oracle import `objdel.checker.polygonize.polygonize` by module path" with "Tests
import the oracle from `tests.oracles.polygonize_oracle`". D3's file list and CLEAN11's owner paragraph say "B1.4a
moved `polygonize.py`; D3 moves `bitquad.py`". This also meets the owner rule "no old code beside new code" one phase
earlier.

### 2. minor: from B1.4a to D3 the mask-parity tests and scripts compare the engine with itself.

Sources: B (B1.4a), D (D2, D3).

- `experiments/cuda_checker/compare.py:21` is `from objdel.checker import bitquads, polygonize`, and `mask_reference`
  (:165-169) builds the reference rings with `polygonize(mask)` [code]. `tests/test_cuda_checker.py:131-174`
  (`test_fuzz_masks_match_polygonize_and_bitquads`, `test_phase0_lidc_masks`) and
  `tests/test_cuda_checker_gpu.py:46-52` (`test_large_mask_batch_matches_polygonize`) use `C.mask_reference`
  [code]. After B1.4a, `objdel.checker.polygonize` is `rings.polygonize`, that is `_gpucheck.check_masks`. So the ring
  half of these tests compares `check_masks` with `check_masks` (only the bit-quad half still uses numpy). The same
  holds for `experiments/cuda_checker/eval.py:42,247`, `experiments/phase1/bench/xarch_hashes.py:45` and
  `K1/cross_numpy.py:16`, the Mac half of the K1 runbook.
- B1.4a re-points only `src/` callers ("After B1.4a no `src/` module calls it"). D3 re-points the oracle callers
  ("change the oracle callers: ... `experiments/cuda_checker/compare.py` (:21, 167)"). D2 comes between them and says:
  "change `tests/test_cuda_checker.py`: port K1-draft's generators ... Compare against the numpy oracle". An engineer
  who adds that family through the file's existing `mask_mismatches` helper compares K1 with K1.
- Not blocking or major: B1.4a's `tests/test_rings.py` compares `polygonize` with the oracle on the fuzz families and
  the 37,093 crops (Mac, CPU build). D3 lands right after D2 and fixes the callers before R-D2 (D-G3 "Mac + R-D2"). V1
  runs no mask-parity test on the box.

Fix: the same change as item 1 (B1.4a re-points the polygonize oracle callers). If item 1 is not taken, add one line
to B1.4a: "the same commit re-points `compare.py:21,167`, `eval.py:42,247`, `xarch_hashes.py:45` and
`K1/cross_numpy.py:16` to `objdel.checker.polygonize.polygonize` by module path, so parity keeps a numpy reference
until D3". Add one line to D2: "the new fuzz family goes through `compare.mask_reference`, which reads the oracle".

### 3. minor: `train_rl --merge-last` (B1.5 step 3b, round 11) has no `RlStage` form; `RlStage` requires `train` and `steps`.

Sources: E (E2.1, E2.3, E6 step 0), B (B1.5).

- B.md "#### B1.5", V1 step 3b: "`train_rl --merge-last --model qwen3-vl-2b --ckpt {rd}/ckpt --out {rd}/b1-rows` (a
  hard step). ... In `--merge-last` mode `train_rl` reads no rows and trains nothing, so `--steps` and `--train` are not
  required there (argparse checks them only without the flag)."
- E.md "### E2" step 1: `RlStage` has `train: str` and `steps: int = Field(gt=0)` with no default, and
  `@model_validator ... def _save(self): ... save_steps <= steps`. E2 step 3: "`run_rl(stage: RlStage, rd: Path)`
  (the body of `train_rl.main` after argument parsing: resume from the last checkpoint, train, merge LoRA, save
  `<out>/final`, write `rl.json`)". E.md has no hit for `merge-last` or `merge_last` [measured: grep].
- Golden 5 records `b1_rows.yaml` (E0), and E6 step 0 must parse "every golden-5 argv list of the `rl` steps" into
  `RlStage`. Step 3b's argv has no `--train` and no `--steps`, so the spike and E-G9 fail on the Mac at $0. So the
  item is minor. But the engineer must still choose the form: a new field, a new stage kind, or a `CmdStage`.

Fix (E.md E2 step 1 and step 3, E6 step 1): add `merge_last: bool = False` to `RlStage`. Make `train: str | None =
None` and `steps: int | None = Field(None, gt=0)`. The `_save` validator becomes: without `merge_last`, `train` and
`steps` are required and `save_steps <= steps`. With `merge_last`, `train`, `steps`, `resume` and `rows_*` must be
unset. `FLAG_ORDER["rl"]` puts `merge_last` first, before `model`, because B1.5 step 3b writes `--merge-last` first.
Golden 5 then pins step 3b string for string. `run_rl` starts with `if stage.merge_last:
return merge_last(stage, rd)` (no `session_budget`, no rows). Add one U-case: step 3b's argv round-trips through
`cli_parse`/`to_argv`.

### 4. minor: stale "D3 puts the polygonizer in production" text after round 11.

Sources: E, 00_front, 99_back.

- E.md "### E.3.3", `data/*` row: "D3 put the polygonizer in production on the engine before E". After round 11,
  B1.4a does that and D3 only adds `trace`. E.md "## E.0" (system-design row): "one polygonizer seam (D3's
  `rings.polygonize`)". E3a import-rewrite paragraph: "D adds importers, for example `checker/rings.py`". B1.4a adds it.
- 99_back.md "Skills used" (architecture row): "0015 polygonizer exception". The front ADR table says "0015 | Scheduled
  RL rows" (00_front.md section 4).
- 00_front.md section 1, D bullet: "D. CUDA core. The CUDA polygonizer in production, ...".

Fix: E.3.3 row: "B1.4a put the polygonizer on the engine before B1.4; D3 added `trace`". E.0 row: "one polygonizer
seam (`rings.polygonize`, B1.4a, extended by D3)". E3a: "B1.4a and D add importers (`checker/rings.py`,
`eval/device.py`)". 99_back skills row: "0015 scheduled RL rows". Front D bullet: "K1 inside the engine behind the
B1.4a seam (D2-D3)".

### 5. minor: no phase owns the doc update when B1.4a switches the production polygonizer.

Sources: DOC (DOC9), B.

- DOC.md "### DOC9": overview section 5 gets a polygonizer row "(write the row that is true on the commit DOC9 lands
  on)". The owner table's B1 row updates only "`reference/data.md` (view schema, counts), `guides/training.md:106,151`,
  `reference/configs.md`". The D row owns `architecture/overview.md` section 5 and ADR 0016 "CUDA polygonizer in
  production" (D9).
- If DOC9 lands before B1.4a, its row says "CPU, numpy". B1.4a then changes production, and no row in the owner table
  tells B1 to update the overview row or `reference/package.md` (new `checker/rings.py`). ADR 0016 is written at D9,
  one phase after the switch. DOC10's rule says "a pull request that changes ... a CUDA contract adds or supersedes an
  ADR". B1.4a changes no engine contract, so the rule does not fire.
- Docs only, so minor.

Fix: DOC9 owner table, B1 row "Updates": add "`architecture/overview.md` section 5 polygonizer row and
`reference/package.md` checker section (B1.4a: `checker/rings.py`)". D row: "ADR 0016 records B1.4a as the switch and
D2-D3 as the K1 port".

### 6. minor: `LIVE_CLIS` does not list B1.4a's new live script `verify_target_hashes`.

Sources: E (E8 step 2).

B1.4a creates `experiments/phase2/verify_target_hashes.py` (`-I`), and D3 step 3 re-runs it. E8 step 2's live list
and `LIVE_CLIS` do not name it, so E8's `-m experiments.*` string gate fails on its usage line at $0. Round-11 arch
item 5 covers the same class of gap and is still open. Fix: when that item is fixed, add `verify_target_hashes` to the
E8 step 2 list.

## Counts

blocking 0, major 0, minor 6 (items 1-6).
