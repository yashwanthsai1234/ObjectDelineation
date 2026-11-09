# CLEAN phase: adversarial, property, performance and failure-mode tests (set 2)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA catalog skill is needed for repo-hygiene tests; GPU tests only reuse the existing `objdel._gpucheck` extension. No change.
- engineering:testing-strategy: pyramid split (many unit/property tests, few integration, one-off gates as scripts).
- mattpocock-skills:tdd: tests sit at public seams only (`code_files`, `main`, `boxes_meet`, `grpo_config`, archive layout), assert behaviour not internals, and expected values come from independent sources (brute-force rasterisation, shapely, literals captured before CLEAN).
- Not invoked (no effect on a test list): code-review, debug, system-design, architecture, documentation, huggingface-best.

Sources read: `CLEAN.spec.md`; `experiments/phase2/vast/policy.py` (`code_files`, line 374, uses `git ls-files -z -co --exclude-standard` plus one regex), `tests/test_vast.py` (fake-binary style, `test_code_files_keep_paths_with_spaces_whole`), `src/objdel/data/lidc.py:257` and `src/objdel/targets/frame.py:124` (`crowded`), `src/objdel/train/configs.py` (`processor_kwargs`, `grpo_config`, `MODELS`), `src/objdel/train/rewards.py:107` (`trl_reward`), `experiments/phase2/train_rl.py` (`--procs` line 35, `fork_pool` line 50), `Makefile`, `pyproject.toml` pytest table (`testpaths = ["tests"]`). Facts: [code: file:line] as cited. Environment facts [measured: this session]: `hypothesis 6.170.2` is installed; `vulture` and `shapely` are NOT installed in `.venv`; 802 files are tracked/untracked-visible by `git ls-files`.

Conventions. Level = unit / property / integration / GPU / e2e. IDs `T01...`. "Pre-CLEAN literal" = a value the owner captures with a script on the commit tagged `pre-clean-2026-11-08` before any change, and pastes as a literal into the test (never recomputed by the code under test).

## Spec gaps found while writing the tests (the implementer must resolve these)
1. B1 "`riga*.csv|json|jsonl`" is ambiguous: today's regex `riga[^/]*\.(csv|jsonl?)$` matches `riga` anywhere in the basename (`xriga.csv`), and ignores `.CSV` only through `re.I`. Tests T05 pin: basename containing `riga` anywhere, any case, extension csv/json/jsonl. `src/objdel/data/riga.py` must stay.
2. B1 "`/lab/` or `/results/` directly under `experiments/<phase>/`": nested `experiments/p/sub/lab/x.py` is not covered. T05 pins "kept" (spec literal). Confirm with owner.
3. `code_files()` uses `git ls-files -c`, which also lists tracked files that were moved/deleted in the working tree but not yet `git rm`/added. During CLEAN, before the owner commits, `rsync --files-from` would hit missing files. Spec silent. T04 asserts every returned path exists on disk.
4. B5 requires backticked paths in `README.md` and `docs/**` to exist, but C9 says the owner's uncommitted `README.md`/`docs/` stay out of CLEAN commits, and `docs/` is untracked today, so a fresh clone has none (test is vacuous in CI). T26 states the rule and tests the checker on synthetic markdown; the live run is a local gate only.
5. B5 hygiene paths such as `data/...` exist locally but not in CI (data is git-ignored). Exempt `data/` or the test fails in CI.
6. B9 / B22 / B14 string scans: the test file that forbids `smolvlm`/`lfm`/`experiments/phase1/lab` would contain the literal and fail itself. Build the needle by concatenation (`"smol" + "vlm"`). Also "lfm" is a substring of ordinary words (`selfmade`, `shelfmark`) and can occur inside hashes in JSON fixtures. Use `(?<![a-z])lfm` on text source (`.py .yaml .md .toml`), never on `.json/.npz`, and never scan `__pycache__` (use `git ls-files -co --exclude-standard`).
7. B7 caps arrays at `ndim <= 1` and non-bool: a flattened `uint8` mask (262144 values) passes. T14 adds a size cap [estimate: 4096 elements]. A 1-D per-case RIGA vector also passes "aggregate only"; T14 requires each RIGA array to have size 1.
8. B7 "at least 99 %": integer arithmetic (`100*match >= 99*N`), and N must equal the golden's recorded case count or an empty/truncated fixture passes vacuously.
9. B13 "prints nothing": vulture's exit codes are 0 (clean), 3 (dead code), 1 (bad input). Test must assert both `returncode == 0` and empty stdout+stderr, plus a canary proving the config still detects dead code.
10. B19 version check: `vastai --version` output format is not in the spec. T37 uses exact-token comparison so `1.11.0` does not match `1.11.01` or `11.11.0`.
11. B1a: under `python -I` neither cwd nor `PYTHONPATH` is on `sys.path`, and `objdel` is editable-installed from REPO. The test must insert the copy first and then assert on every module's `__file__`; absence of an `ImportError` proves nothing because the editable finder would silently serve REPO files.

---------------------------------------------------------------------------------------------------------------

## A. B1 / B1a: bundle manifest and bundle import

### T01 `code_files_never_returns_a_banned_root` (property, unit)
- Setup: tmp git repo (pattern of `test_code_files_keep_paths_with_spaces_whole`), `monkeypatch.setattr(P, "ROOT", tmp)`. Create, tracked and untracked mixed: `tests/oracles/o.py`, `.claude/settings.json`, `.claude/worktrees/w1/src/a.py`, `archive/phase0/x.py`, `archive/smolvlm/README.md`, `data/raw/a.csv`, `docs/a.md`, `src/objdel/__init__.py`, `experiments/phase2/train_rl.py`.
- Input: `P.code_files()`.
- Expected: exactly `["experiments/phase2/train_rl.py", "src/objdel/__init__.py", "tests/oracles/o.py"]` (sorted).
- Why: pins the whole B1 contract at once; the old regex dropped only `data|docs`, so `.claude/` (2 files) and `archive/` (10 files, G6) leaked.

### T02 `code_files_does_not_over_match_prefixes` (adversarial unit)
- Setup: tmp repo with these files: `src/objdel/data/lidc.py`, `src/objdel/data/riga.py`, `src/docs/x.py`, `tests/data/loader.py`, `archive_tools/x.py`, `docs_tools/x.py`, `data_utils.py`, `.claudeignore`, `.claude_x/a.py`, `experiments/phase2/lab_utils.py`, `experiments/phase2/results_table.py`, `experiments/phase2/labs/a.py`, `experiments/phase2/resultsets/a.py`, `tests/oracles/o.py`.
- Expected: ALL of them returned (14 files). Anchors are directory boundaries at repo root only.
- Why: a regex like `data/` unanchored, or `archive` without a slash, would drop live library code (`src/objdel/data/lidc.py`) from a bundle and break the GPU box after money is spent. The spec says `data/` roots only at top level.

### T03 `code_files_riga_table_rule` (property, parametrised)
- Input names (all tracked): drop: `experiments/phase1/riga_cases.csv`, `experiments/phase1/RIGA_Cases.CSV`, `experiments/x/xriga.json`, `tests/goldens/riga-agg.jsonl`, `a/b/Riga1.JsonL`. Keep: `src/objdel/data/riga.py`, `experiments/phase1/riga_notes.txt`, `riga/readme.py`, `tests/test_data_riga.py`, `experiments/phase1/rigaX.csv.py`.
- Expected: dropped set equals the five; kept set all present.
- Why: C8 forbids RIGA masks or per-case results leaving the machine (D16); case-insensitivity and the "anywhere in basename" reading are the places a regex rewrite regresses. Extension must be the end of the name (`rigaX.csv.py` kept).

### T04 `code_files_lists_only_existing_files_no_duplicates` (property, adversarial)
- Setup: tmp repo; `git add -A`, commit nothing needed; then `os.remove("src/old.py")` (tracked, deleted in working tree), `shutil.move("experiments/phase1/a.py", "archive/phase1_scripts/a.py")` without `git add`.
- Expected: every returned path satisfies `(tmp/p).is_file()`; `len(set(files)) == len(files)`; `archive/phase1_scripts/a.py` is absent (banned root), `experiments/phase1/a.py` absent (file gone).
- Why: `git ls-files -c` reports index entries even when the file was deleted, so the mid-CLEAN working tree makes `rsync --files-from` abort with "No such file". Expected to FAIL on today's code: the implementer must filter by existence (spec gap 3).

### T05 `code_files_nested_lab_results_and_worktree_nested_repo` (adversarial unit)
- Setup: files `experiments/phase1/lab/a.py`, `experiments/phase2/results/r.json`, `experiments/phase2/sub/lab/b.py`, `lab/c.py`, and `.claude/worktrees/w1/.git` (a file `gitdir: ../x`) with `.claude/worktrees/w1/src/z.py`.
- Expected: first two dropped; `experiments/phase2/sub/lab/b.py` and `lab/c.py` kept (spec literal, gap 2); nothing under `.claude/` returned even when git treats `w1` as a nested repo.
- Why: nested git dirs make `ls-files -o` print a directory entry (`.claude/worktrees/w1/`) rather than files; a filter on file suffix would let a directory path reach rsync.

### T06 `code_files_without_oracles_stops_with_error` (failure mode)
- Setup: tmp repo with only `src/objdel/__init__.py`; and a second case with `tests/oracles_old/x.py` only; and a third with `archive/tests/oracles/o.py` only.
- Expected: `SystemExit` (policy `die`) in all three, message contains `tests/oracles`.
- Why: the existence check must run AFTER filtering. If it runs before, the archived oracle satisfies it and a bundle without the K3 harness ships (import fails on the GPU box, `IMPORT_CHECK` in sync_up).

### T07 `code_files_untracked_not_ignored_in_ignored_out` (unit)
- Setup: tmp repo with `.gitignore` containing `build/` and `*.pyc`; untracked `src/new.py`, ignored `build/x.so`, `src/__pycache__/n.pyc`.
- Expected: `src/new.py` returned; the other two absent.
- Why: spec says untracked, not-ignored files still ship (a new file not yet `git add`ed must reach the box).

### T08 `gitignore_covers_claude_worktrees_only` (unit)
- Setup: tmp repo with the repo's real `.gitignore` copied in.
- Input: `git check-ignore -q <path>` for `.claude/worktrees/a/b.py` (exit 0), `.claude/worktrees` (0), `.claude/settings.json` (exit 1), `foo/.claude/worktrees/x` (exit 1, pattern is root-anchored), `.claude/worktrees_x/a` (1).
- Why: the anchor `/` and the trailing-slash semantics are the failure; ignoring all of `.claude/` would hide the owner's settings and is not what B1 asks.

### T09 `bundle_copy_imports_from_the_copy` (integration, no GPU; B1a, G14) [pytest sketch S1]
- Setup: `dest = tmp_path/"bundle"`; copy exactly `P.code_files()` (real repo) into `dest`; locate `objdel/_gpucheck*.so` via `importlib.util.find_spec("objdel._gpucheck").origin` and copy it to `dest/src/objdel/`. Run `sys.executable -I -c SCRIPT dest`.
- SCRIPT: `sys.path.insert(0, f"{dest}/src"); sys.path.insert(0, dest)`; import `objdel, objdel._gpucheck, experiments.phase1.reward_gates, tests.oracles.reward_harness, experiments.phase2.{train_rl,g_probe,select_ckpt,final_eval}`; then print JSON of `[(name, m.__file__) for name, m in sys.modules.items() if name.split(".")[0] in {"objdel","experiments","tests"}]`.
- Expected: exit 0; list non-empty (>= 12 modules); every `__file__` starts with `str(dest)`; none starts with the REPO root; `data/`, `docs/`, `archive/` absent from the copy.
- Why: the real bundle was once broken by a missing file that Python silently satisfied from the editable install on the dev box (spec gap 11). Only the `__file__` assertion catches it. Catches a banned root that holds a needed module (e.g. a live import reaching `archive/`).

### T10 `bundle_import_canary_detects_leak` (meta test of T09)
- Setup: as T09 but delete `dest/experiments/phase1/reward_gates.py` from the copy first.
- Expected: the T09 assertion helper fails (either ImportError or a REPO-rooted `__file__`); the test asserts the helper raises `AssertionError`.
- Why: proves T09 is not vacuous.

### T11 `bundle_size_and_count_gate` (performance/regression)
- Input: `len(P.code_files())`, `sum(size)`.
- Expected: count <= 802 - 2 (.claude) - 10 (archive) - 426 (lab/phase0 moved, G5, but they live under `archive/`, so already dropped) ... concretely: count after CLEAN <= pre-CLEAN count and total bytes <= pre-CLEAN bytes + 50 KB (new fixtures, goldens, README files). Capture both pre-CLEAN literals [measure before].
- Why: the bundle is rsync'd per rental; accidental inclusion of `archive/` or goldens bloats upload time on a metered box.

### T12 `code_files_speed_gate` (performance)
- Setup: tmp repo with 5000 tracked files in 200 dirs.
- Expected: `code_files()` returns in < 3 s [estimate: one git call + one regex pass]; time is O(n) (10000 files takes < 2.5x of 5000).
- Why: catches an accidental per-file `git` subprocess or an O(n^2) existence check introduced for T04.

---------------------------------------------------------------------------------------------------------------

## B. B2: `train_rl` without `--procs`

### T13 `train_rl_rejects_procs_in_every_spelling` (adversarial, failure mode) [pytest sketch S2]
- Input: `train_rl.main(argv)` for argv = required options plus each of `["--procs","2"]`, `["--procs=2"]`, `["--proc","2"]` (abbreviation), `["--procs"]`, `["-procs","2"]`.
- Expected: `SystemExit` with `code == 2`; captured stderr contains `unrecognized arguments` (or `ignored explicit argument`); `serve.preload`, `serve.fork_pool` and `gc.freeze` mocks are never called.
- Why: with `allow_abbrev` default True, a leftover option sharing the prefix `--proc` would swallow the abbreviation silently; the exit must happen before any data preload (a bad flag must not cost a 5-minute preload on a rented box).

### T14 `train_rl_help_has_no_procs_and_keeps_other_options` (unit)
- Input: `main(["--help"])` under capsys (`SystemExit(0)`).
- Expected: `--procs` not in stdout; `--model --ckpt --train --train-views --out --steps --save-steps` all in stdout (the full option list is a pre-CLEAN literal captured by `ap._actions`; values/defaults equal the literal except the removed one).
- Why: "all other options unchanged" can be broken by deleting one line too many.

### T15 `train_rl_creates_no_process_pool` (adversarial integration, stubbed) 
- Setup: reuse the stub harness of `tests/test_p2_train_rl.py`; monkeypatch to raise: `multiprocessing.Pool`, `multiprocessing.pool.Pool`, `multiprocessing.get_context`, `concurrent.futures.ProcessPoolExecutor`, `os.fork`, `objdel.train.serve.fork_pool`.
- Input: run `main([... valid args ...])` up to the first training call (stubbed).
- Expected: main returns 0; no patched object was called. Plus a static test: AST of `train_rl.py` contains no `Name`/`Attribute` called `fork_pool`, `pool_size`, `Pool`, and no string `"--procs"`.
- Why: B2 says "creates no process pool"; a forked pool after CUDA init corrupts the process (D14).

### T16 `train_rl_gc_freeze_once_after_last_preload` (unit, ordering) [pytest sketch S3]
- Setup: shared `calls: list[str]`; monkeypatch `serve.preload` -> appends `"preload"`; `gc.freeze` -> appends `"freeze"`; stub the heavy parts; run with `--train-views a b` (two extra roots, so `preload` is called more than once) and again with no views.
- Expected: `calls.count("freeze") == 1`; `calls.index("freeze") > len(calls) - 1 - calls[::-1].index("preload")` (freeze after the last preload); `calls[-1]` may be anything after the freeze but no `"preload"` follows it.
- Why: freezing before the last preload leaves the final cache unfrozen (copy-on-write page growth); freezing twice or inside a loop is the typical regression. Multiple preload roots make "last" meaningful.

### T17 `train_rl_gc_freeze_not_reached_on_bad_args_and_not_on_failure` (failure mode)
- Input: (a) `--steps` missing, (b) preload raising `FileNotFoundError`.
- Expected: (a) exit 2, no freeze call; (b) exception propagates, `gc.freeze` never called, no leftover monkeypatch side effect (`gc.get_freeze_count()` equals its value before the test).
- Why: a freeze in a `finally` or before validation would hide the real error behind a crash in the stub, and would leak frozen objects into the rest of the pytest session.

### T18 `rl_json_keys_unchanged` (regression)
- Setup: the stub run of T15 writes `out/rl.json`.
- Expected: `sorted(json.load(...))` equals the pre-CLEAN literal list (captured by running the old code under the same stub); no key `procs` anywhere (recursive).
- Why: downstream (`select_ckpt`, `g_probe`) reads these keys; the spec says unchanged.

---------------------------------------------------------------------------------------------------------------

## C. B10: `boxes_meet`

### T19 `boxes_meet_matches_pixel_overlap_oracle` (property, independent oracle) [pytest sketch S4]
- Setup: `hypothesis` strategies: `x0 < x1`, `y0 < y1` in `[0, 12]`, 0..4 other boxes with `a < c`, `b < d` in `[0, 12]`.
- Oracle (independent of the formula): rasterise each box into a 12x12 bool grid (`grid[y0:y1, x0:x1] = True`); expected = `any((g_box & g_other).any())`.
- Expected: `boxes_meet(box, others) == oracle`, for tuples, lists, and one-shot generators (`iter(others)`), and for numpy `np.int64` coordinates.
- Why: the formula `a < x1 and x0 < c ...` is easy to rewrite as `<=` ; the half-open rule has an independent meaning (shared pixels). 600+ examples cover every touching/abutting configuration.

### T20 `boxes_meet_edge_cases_literal` (unit, parametrised)
| box | others | result | reason |
|---|---|---|---|
| (0,0,10,10) | [] | False | empty iterable |
| (0,0,10,10) | [(10,0,20,10)] | False | edge touching, half-open |
| (0,0,10,10) | [(9,0,20,10)] | True | one pixel column shared |
| (0,0,10,10) | [(0,10,10,20)] | False | touching below |
| (0,0,10,10) | [(-5,-5,1,1)] | True | corner pixel (0,0) |
| (0,0,10,10) | [(2,2,8,8)] | True | contained |
| (2,2,8,8) | [(0,0,10,10)] | True | containing |
| (0,0,10,10) | [(20,20,30,30), (5,5,6,6)] | True | hit is last |
| (0,0,10,10) | [(5,5,5,9)] | True | zero-width other, literal formula (documented) |
| (5,5,5,5) | [(0,0,10,10)] | True | zero-area box, literal formula (documented) |
| (0,0,10,10) | np.empty((0,4),int) | False | empty array |
- Why: pins today's behaviour, including degenerate boxes the formula accepts; if the implementer "fixes" degenerate boxes, `crowded` outputs (and the dataset it filters) change silently.

### T21 `boxes_meet_algebraic_properties` (property)
- Properties (hypothesis, random boxes incl. degenerate): symmetry `meets(A,[B]) == meets(B,[A])`; union `meets(A, X+Y) == meets(A,X) or meets(A,Y)`; monotone (adding a box never turns True to False); translation invariance `meets(A+t, X+t) == meets(A,X)`; order independence (shuffle `others`).
- Why: independent of the formula; a swapped variable (`b < y0`) breaks symmetry or translation at once.

### T22 `crowded_wrappers_equal_boxes_meet` (property + regression)
- Setup: fake `Slice` (`SimpleNamespace(others=[...])`) and fake `Item` (`SimpleNamespace(meta={...})`).
- Expected: `lidc.crowded(s, box) == boxes_meet(box, s.others)`; `frame.crowded(item, box)` equals `boxes_meet(box, item.meta.get("other_nodules", ()))` for `meta == {}` (False), `{"other_nodules": []}` (False), key present with one overlapping box (True). `frame.crowded` with `{"other_nodules": None}` keeps today's behaviour (raises `TypeError`); the test pins whatever the pre-CLEAN code does [measure before].
- Plus data test (marker `data`): `crowded` over all LIDC fixture slices x 5 windows equals the pre-CLEAN list of results, hashed (SHA-256 of the bool bytes), pre-CLEAN literal.
- Why: two copies of the rule were merged; the dataset split by crowding must not change (C1: target-set hashes).

### T23 `boxes_meet_performance_gate` (performance, relative)
- Setup: 20 000 `box`es, 40 `others` each; reference = the old inline generator expression defined in the test.
- Expected: new function time <= 1.5 x reference time (same process, best of 3); results equal.
- Why: `crowded` is called per window in dataset builders; wrapping must not add overhead. Relative gate is machine independent.

---------------------------------------------------------------------------------------------------------------

## D. B9 (SmolVLM) and config behaviour

### T24 `smolvlm_gone_from_configs` (unit)
- Expected: `"smolvlm-256m" not in MODELS` and not in `MAX_COMPLETION`; `processor_kwargs(m, d) == {"padding_side": "left"}` for every `m in MODELS` and `d in ("lidc","riga")` (exact dict equality, no extra keys such as `size`); `processor_kwargs("qwen3-vl-2b", "anything")` also works (the old dict lookup on `dataset` must not remain).
- Why: the old code did `{"riga":1024,"lidc":512}[dataset]` only for smolvlm; a half-edit leaves a `KeyError` for unknown datasets or a stray `size` key that changes tokenisation for Qwen.

### T25 `grpo_config_use_vllm_override_matrix` (property, parametrised)
- Cases: no override -> `True`; `use_vllm=False` -> `False`; `use_vllm=True` -> `True`; `use_vllm=None`: pin pre-CLEAN behaviour; `grpo_config("qwen3-vl-2b", out, use_cpu=True, use_vllm=False)` constructs (no exception) with `use_vllm is False`.
- Also: with the `use_vllm` default the full config equals the pre-CLEAN literal dict for the model (`dataclasses.asdict` subset: `loss_type, epsilon, epsilon_high, beta, num_generations, per_device_train_batch_size, max_completion_length, vllm_gpu_memory_utilization, vllm_enable_sleep_mode, vllm_max_model_length, learning_rate, reward_weights`).
- Why: overrides are merged with `|` after the default; reordering that expression makes the default win over the caller and a CPU test starts vLLM.

### T26 `forbidden_words_scan_is_exact_and_not_self_tripping` (unit + integration) [part of hygiene B14]
- Part 1 (checker unit): the scan helper over synthetic texts: `selfmade`, `shelfmark`, a base64 blob `aGVsZm1Gx...` (contains `lfm`), and a `.json` fixture containing `lfm` must not match; `LFM2.5-VL`, `lfm2.5-vl-450m`, `# uses LFM` must match; `SmolVLM`, `smolvlm_256m`, `SMOLVLM` match.
- Part 2 (repo): over `git ls-files -co --exclude-standard` text sources in `src/`, `experiments/phase2/*.py`, `tests/**/*.py`: zero `smolvlm` matches (case-insensitive). Needles are written as `"smol" + "vlm"` and `r"(?<![a-z])" + "l" + "fm"`.
- Why: gap 6; a self-matching scan test fails forever, and an over-broad `lfm` scan blocks Wave 3 on false positives.

---------------------------------------------------------------------------------------------------------------

## E. B3 / B4 / B5 / B6 / B15: moves, archive, imports, docs

### T27 `retired_paths_hold_nothing_on_disk_or_in_index` (unit, git read-only)
- Input: for each path in B3: `git ls-files -- <path>` empty AND `not (ROOT/path).exists()` (stale `__pycache__` leaves the directory, and an empty directory still imports as a namespace package).
- Plus subprocess `python -I -c "importlib.import_module(...)"` with `sys.path` set to REPO for each of the six importable phase1 modules and `experiments.phase1.bench.xarch_hashes`; and the negative: `import experiments.phase0` and `import experiments.phase1.lab` raise `ModuleNotFoundError`.
- Why: the negative import is the real behavioural meaning of "removed".

### T28 `archive_layout_and_readme_contract` (property over the folder)
- For each direct subfolder `D` of `archive/`: `D/README.md` exists; headings `## What it was`, `## What replaced it`, `## Retiring commit`, `## How to restore` each appear exactly once, in that order, as whole lines, each followed by a non-empty body before the next heading; the restore body contains `pre-clean-2026-11-08`; `archive/__init__.py` absent and no `__init__.py` anywhere under `archive/`; expected set of folders is a superset of `{phase0, phase1_lab, phase1_scripts, phase2_lab, smolvlm}` plus `controller` iff `not hasattr(callbacks, "Controller")` is false (i.e. B8 default taken). No stray file at `archive/` root other than folders and an optional top README.
- Why: heading drift (`## How to Restore`) and a missing tag name defeat the "nothing deleted without a copy" promise (C2).

### T29 `archive_is_never_collected_linted_or_scanned` (integration)
- Checks: (a) `pytest --collect-only -q` (subprocess) lists no node id starting with `archive/` even if an archived `test_x.py` exists (the test plants `archive/_probe/test_probe.py` in a tmp copy of pyproject + rootdir, or relies on `testpaths = ["tests"]`; pin it); (b) `tomllib` reads `[tool.ruff] extend-exclude == ["archive"]`; (c) `ruff check archive --force-exclude` reports no files checked / exit 0 even with a deliberate syntax error planted in a tmp archive file; (d) `[tool.vulture].paths` contains no `archive`.
- Why: archived legacy code has syntax and import errors by design; one collected file turns `make test` red.

### T30 `archived_tests_are_listed_by_node_id` (property, C5)
- For every `def test_*` in `archive/**/test_*.py` (AST), the string `<archive-relative path>::<name>` appears in that folder's README. One-time gate (script, G2): `pytest --collect-only -q` ids at tag `pre-clean-2026-11-08` (saved list) minus ids now collected must each appear in some archive README or be in an explicit rename map file in the PR.
- Why: removed tests are the cheapest way to make `make test` green falsely.

### T31 `tests_do_not_touch_lab_or_archive` (property, B5/G7)
- AST walk of every `tests/**/*.py`: no `Import`/`ImportFrom` whose module starts with `archive` or contains `.lab.` / equals `lab`; no call `sys.path.insert/append` with an argument whose source contains `lab` or `archive`; no `Constant` string (outside docstrings) containing `experiments/phase1/lab`, `experiments/phase2/lab`, or `archive/` used with `open/read_text/Path(...)`. Docstrings and comments citing `archive/...` are allowed (test must use the AST, not grep, so docstrings are exempt).
- Also `src/objdel/**`: no import of `experiments` or `tests` (including `importlib.import_module("experiments....")` string literals, and lazy imports inside functions).
- Why: AST resolves aliasing (`import archive.x as a`), grep does not; lazy in-function imports are where `src` leaks to `experiments`.

### T32 `repo_paths_in_docs_exist_checker` (unit on synthetic markdown, plus local gate)
- Checker unit cases (tmp `README.md`): backticked `src/objdel/geometry/box.py` (exists -> ok); `src/objdel/gone.py` (missing -> reported); `src/objdel/train/configs.py::grpo_config` (strip `::sym`, ok); `experiments/phase2/configs/*.yaml` (glob, must match >= 1 file); `<phase>` placeholders and `data/...` skipped; fenced code blocks skipped; URLs skipped; `make test`, `.venv/bin/python`, `python -I -m x.y` not treated as paths; trailing punctuation `` `a/b.py`, `` handled; paths inside exempt files (`docs/changelog.md`, `docs/research/**`, `experiments/*/impl/**`, frozen phase `.md`) never reported.
- Local gate: the same checker over the real `README.md` and `docs/**/*.md` returns `[]` (skipped with a clear reason, not silently, when `docs/` is not tracked; see gap 4).
- Why: a naive regex flags commands and globs (flaky test that gets disabled) or misses paths (useless).

### T33 `fixtures_moved_and_constants_inline` (unit)
- Expected: `tests/fixtures/vast_snapshot.json` and `tests/fixtures/topo_probe_v3e_seed*.json` exist (at least one `seed*`), parse with `json.loads(..., parse_constant=raise)`; `tests/test_vast.py` has no `Path`/string containing `lab`; the Dockerfile digest check compares against a `sha256:` literal defined in the test (regex `[0-9a-f]{64}`) and equals the digest in `experiments/phase2/docker/Dockerfile`.
- Why: B6; the pre-move test pointed at `experiments/phase2/lab/A3/vast_snapshot.json` (current `SNAP`), a path the archive step deletes.

---------------------------------------------------------------------------------------------------------------

## F. B7 / B7a: goldens

### T34 `goldens_hold_no_masks_and_riga_only_aggregates` (property, adversarial) [pytest sketch S5]
- Helper `check_goldens(dir)` runs over `tests/goldens/*.npz` except `polycore.npz`; `np.load(p, allow_pickle=False)` must succeed (no object arrays).
- Expected for `targets_lab.npz`: every array `ndim <= 1`, `dtype != bool`, `size <= 4096` [estimate: a 512x512 mask flattened is 262144, a crop at least 64x64 = 4096], file size < 1 MB; keys starting `riga` are exactly three and each has `size == 1`; no key matches `riga/.+/.+` (per-case); LIDC keys all match `^lidc/[^/]+/[^/]+$`; number of distinct LIDC items >= the case count recorded in a literal.
- Meta (helper test on synthetic npz in tmp): a 2-D `uint8`, a 1-D bool, an object array, a 1-D `uint8` of 100 000 values, and a `riga/case1/r1` key must each make the helper raise; a clean file passes.
- Why: C8 is a hard privacy/legal rule; `ndim <= 1` alone is defeated by `.ravel()`.

### T35 `differential_vs_lab_golden_is_not_vacuous` (integration, marker `data`)
- Properties of `test_fixtures_differential_vs_lab_golden`: N = number of LIDC cases read from `data/fixtures` == N recorded in golden (assert before the 99 % rule); N > 0; match rule uses integers: `100*match >= 99*N`; RIGA: `prod_total <= lab_total`; recorded crop hashes checked before the builder runs (fail loudly on a crop mismatch, never "fix" the golden).
- Mutation check (data marker, slow-ish): monkeypatch the production builder's simplification tolerance (or drop the last vertex of one case) and assert the test fails: one vertex fewer in one case must fail the "vertex count equals lab count per case" check.
- Companion non-data test: golden keys count equals the literal case count, so a truncated golden fails even on a machine without data.
- Why: the old differential test imported the lab builder; the golden version can silently turn into a no-op that always passes (empty loop, wrong key regex).

### T36 `p2_lab_refs_golden_properties` (property + mutation)
- Expected: `tests/goldens/p2_lab_refs.json` parses with `parse_constant` raising (no NaN/Infinity), every leaf float finite; the three named tests exist (`inspect`/AST) and use `rtol/atol` equal to the pre-CLEAN literals (values from the old tests, written down at capture time); none imports a lab module (AST, same walker as T31).
- Mutation: copy the JSON to tmp, scale one value by `1 + 3*rtol`, run the comparison helper against live code: must fail; scale by `1 + 0.1*rtol`: must pass. Seeds in the test are literal; regenerating inputs twice gives identical arrays (determinism).
- Why: a golden that is compared with a looser tolerance than before, or that holds a stale value, makes the live comparison a tautology.

---------------------------------------------------------------------------------------------------------------

## G. B8 (Controller), B16 (volume)

### T37 `controller_gone_trl_reward_intact` (unit + regression)
- Expected: `callbacks.TraceStep` exists; `Controller`, `audit_alarm`, `FIRED` absent; `objdel.eval` has `hack_audit, gates, metrics, summarize, vision` and no `alarms_for`; `TRLReward(...)` has no `recorded` attribute (`vars(r)` has no key `recorded`); for two instances `r1, r2`: `all(f.trl_reward is r1 for f in r1.funcs()[0])` and the same for `r2` (no cross-wiring); the aggregator and every logger have the attribute; `copy.copy(f).trl_reward` still resolves; `fast.py`'s prefetch path (`next(f.trl_reward ...)`) returns the reward when given `SimpleNamespace(reward_funcs=r.funcs()[0])`.
- Regression of numbers: on a seeded synthetic batch (the one in `tests/test_train_rewards.py`), the 16 returned reward floats `==` a pre-CLEAN literal list (bit-equal, `np.array_equal`), and `fingerprint()` equals the pre-CLEAN literal hex string.
- Why: `recorded` fed the Controller; deleting it may drop a side effect that also normalised something; G16 requires the `trl_reward` assertion explicitly.

### T38 `no_controller_words_in_src_or_phase2` (property)
- Word-boundary scan (`\bController\b`, `\balarms_for\b`) over `src/` and `experiments/phase2/*.py` (G10) through the shared scan helper; the test file builds the needle by concatenation.
- Why: dangling references in docstrings/imports surface as vulture noise or stale docs.

### T39 `breakit_wiring_runs_end_to_end_without_controller` (e2e on CPU)
- Run the existing `tests/test_breakit_wiring.py` node (store -> reward -> TRL -> eval -> gates) with `monkeypatch.setattr(callbacks, "__getattr__", ...)` guard: any access to a removed name raises. Expected: passes, and the gates output equals pre-CLEAN literal (gate names + booleans).
- Why: the spec requires the test to still exist; this adds a tripwire so it does not silently pass through a stubbed path.

### T40 `topology_any_behaviour_and_docstring` (property + unit)
- Property (hypothesis, 6x6x6 random bool volumes, 0..3 readers): `topology_any(pred, readers) == (betti3d(pred) in {betti3d(r) for r in readers})` and `== any(betti3d(pred) == betti3d(r))`; empty readers -> `False`; single reader equals pred -> True; a hollow cube (b2 = 1) vs a solid cube -> False; a torus-like ring (b1 = 1) vs a ball -> False; readers supplied as a generator work (consumed once).
- Docstring: `"reward" not in (volume.__doc__ or "").lower()` is too coarse; test that no sentence of the module docstring matches `(?i)\breward\b.*\b(uses|calls|computes)\b.*topology_any` and `topology_any` name is still in `volume.__dict__`.
- Why: B16 says "behaves as before" and not claim the reward uses it (it does not, G16 style honesty).

---------------------------------------------------------------------------------------------------------------

## H. B11 / B12 / B13 / B14 / B17 / B18 / B19: tooling

### T41 `preflight_default_config_loads` (integration, no GPU)
- Setup: `make -n preflight` (dry run, no stamp rebuild needed because `-n`) in a subprocess with `CONFIG` unset in the environment.
- Expected: output contains `--config experiments/phase2/configs/r1_qwen3vl2b.yaml`; does not contain `r1_lfm.yaml` or `rl_shifted.yaml`; `make -n preflight CONFIG=x.yaml` contains `x.yaml`; the config file exists and `experiments.phase2.preflight` loads it (its loader function, e.g. `PF.load(cfg)`; use the same entry as `tests/test_preflight.py`); `rl_shifted.yaml` fails to load as a preflight config with a clear message (it is steps-only) rather than a traceback.
- Gate G12 (one-time script): `make preflight CONFIG=<r1>` JSON, each check's `pass` value equal pre- vs post-CLEAN (diff of sorted `(name, pass)` pairs).
- Why: a Makefile typo makes the pre-rental check (mandatory before money is spent) die on a missing file or, worse, pass on the wrong config; `CONFIG ?=` is overridden by an environment variable `CONFIG` of the same name, so the test must clear it.

### T42 `train_sft_load_items_rename` (unit)
- Expected: `hasattr(train_sft, "load_items")`, `not hasattr(train_sft, "preload")`; `inspect.signature(load_items)` parameter names, order, defaults equal the pre-CLEAN `preload` literal (`["cfg", "model"]` today, annotations `dict`, `str | None`); callers: AST grep of `experiments/` and `tests/` finds no `train_sft.preload` / `from ...train_sft import preload`; `serve.preload` must remain (train_rl uses it, T16): assert `hasattr(serve, "preload")`.
- Why: the rename is trivial to half-do; `serve.preload` has the same name and must not be touched.

### T43 `vulture_gate_exact_and_not_vacuous` (integration; fail, never skip) [part of B14]
- (a) `importlib.metadata.version("vulture") == "2.16"` and `pyproject` dev extra has `vulture==2.16` (parsed with `tomllib`, normalised name). `ModuleNotFoundError` -> test fails with an install hint (use `pytest.fail`, never `importorskip`).
- (b) `subprocess.run([sys.executable, "-m", "vulture"], cwd=ROOT)` (uses the `[tool.vulture]` table; no extra args) -> `returncode == 0` and `stdout + stderr == ""`.
- (c) Canary: in tmp, write a project with a dead function `def _never_called(): ...` and copy of the repo's `[tool.vulture]` table with `paths = ["pkg"]`; the same invocation exits 3 and mentions `_never_called` (proves `min_confidence = 60` and the ignore lists do not mask real dead code).
- (d) `[tool.vulture]` structure: `paths == ["src","experiments","tests"]`, `min_confidence == 60`; each of `ignore_names`/`ignore_decorators` group is preceded by a `#` comment line in the raw TOML (line scan); no file matching `*whitelist*` or `_vulture*` tracked; `make dead` (`make -s dead`) exits 0 with empty output.
- (e) Test fakes keep real parameter names: AST check that no function in `tests/` has a parameter named `_`-prefixed solely to dodge vulture where the production signature differs (heuristic: for fakes assigned via `monkeypatch.setattr(mod, "f", fake)`, `inspect.signature(fake).parameters` names are a superset of the real function's names).
- Why: the gate is the only protection against dead code returning; an over-wide `ignore_names` (like `.*`) makes it pass forever (c catches this). `-m vulture` honours `pyproject.toml` only from the cwd.

### T44 `vulture_runtime_gate` (performance)
- `vulture` over the three roots completes in < 60 s [estimate; measure before]; the whole hygiene test file in < 90 s, no data, no GPU, no network (monkeypatch `socket.socket` to raise).
- Why: `make test` must stay usable on the Mac; hygiene tests that need the network or data break the "needs no data and no GPU" requirement.

### T45 `adr_library_keep_list_is_consistent_with_pyproject` (property)
- Setup: find ADR files `docs/architecture/decisions/NNNN-*.md`; numbers are unique, contiguous from 0001 (0001..0011 plus the new one), the keep-list ADR number = max + 1 at the time of writing; headings of the template (copy the exact heading set from `docs/architecture/decisions/README.md`/the newest ADR) all present.
- Content: sections `KEEP`, `DEFERRED`, `NOT ADOPTED` exist; every DEFERRED row names an owning phase (regex `phase [A-Z]\d*`); parse library names (first backticked token per bullet), normalise per PEP 503 (`lower`, `[-_.]+` -> `-`): every KEEP library is declared in `pyproject.toml` (any extra) or is stdlib; every NOT ADOPTED library is absent from every dependency list (incl. `vastai`); a library appears in exactly one section.
- Why: the ADR is a decision record; if it disagrees with pyproject it misleads the next phase (B19 vastai and B20 shapely are exactly such calls). Normalisation trap: `PyYAML` vs `pyyaml`.

### T46 `vastai_is_isolated_and_pinned` (unit, optional B19)
- Setup: fake `vastai` executables on PATH (pattern of `tests/test_vast.py` `FAKE`).
- Expected: no dependency in any `pyproject.toml` list whose normalised name is `vastai` or `vast-ai`; `pyyaml` still in `dev`; `docs/runbooks/gpu-run.md` contains the pinned version string and `vastai` isolated-install command; preflight version check: fake printing `1.11.0` -> pass; `1.11.1`, `1.11.01`, `11.11.0`, `1.11`, empty output, non-zero exit, and no binary on PATH -> each FAIL with a message (never an exception traceback); `policy.vast("search offers ...")` invokes `subprocess.run` with `argv[0] == "vastai"` (monkeypatched), and `"vastai" not in sys.modules` after the call.
- Why: exact-token matching, and failure (not crash) when the binary is absent, are the two mistakes a `startswith`/`in` check makes; the dependency test protects C10 (Pillow, cryptography pins).

### T47 `environment_pins_unchanged` (regression, C10, one-time gate)
- Script: `uv pip freeze --python .venv/bin/python` hash before/after CLEAN, excluding `vulture`/`shapely`/`pyyaml` additions; `Pillow` version equals literal; `git diff --name-only pre-clean-2026-11-08 -- cuda CMakeLists.txt` empty (C1).
- Also: `fingerprint()` equals the pre-CLEAN literal (G3); `target_set_sha256` literal in `tests/test_p2_views.py` unchanged (`git diff pre-clean-2026-11-08 -- tests/test_p2_views.py | grep target_set_sha256` empty) and the live manifest value (data marker) equal (G4).
- Why: C1 and C10 are binary promises that a big refactor can break indirectly (a dev-extra change re-resolving Pillow).

### T48 `local_data_retirement_B17` (e2e, local only, marker `data`, run with `-I`)
- Expected: the three old dirs absent; `data/views/lidc-bundle-v3` and `data/eval_views/lidc-val-rotated` present and their `rows.jsonl` SHA-256 equal to pre-move literals; `data/_retired/2026-11-08/LIST.txt` has one line per moved dir `path<TAB>size<TAB>sha256`, and each listed sha256 equals `sha256(moved/rows.jsonl)`, each size equals recomputed `du` bytes; moved dir byte count equals pre-move literal (a `mv` across volumes that fails midway leaves a partial copy); no live config/yaml/py references `lidc-bundle` NOT followed by `-v3` (regex `lidc-bundle(?![-\w])|lidc-bundle-v2|probe-precision-r1`, excluding `archive/`, `docs/changelog.md`, `docs/research/**`).
- Why: the plain-prefix trap (`lidc-bundle` is a prefix of `lidc-bundle-v3`) would delete or flag the kept dataset.

---------------------------------------------------------------------------------------------------------------

## I. B20 (optional shapely cross-check)

### T49 `shapely_crosscheck_defect_classes_and_raster` (property/differential; fail if shapely missing, it is a dev dependency)
- Setup: `shapely==2.2.*` in dev extra (installed in `.venv` is NOT today: [measured: ImportError]). Generators: random simple polygons (star-shaped around a centre, 5..40 vertices, int coordinates 0..1000), plus hand-built defects: bowtie (self-crossing), repeated vertex, zero-area spike, hole outside shell, ring touching itself at a vertex, collinear degenerate, <3 vertices.
- Expected: for every mapped defect class the CUDA checker's class (host build of `_gpucheck`) is the one mapped to `shapely.is_valid_reason(poly)` text (`Self-intersection`, `Ring Self-intersection`, `Too few points`, `Hole lies outside shell`, ...); unmapped reasons are asserted to be `not legal` only. For legal polygons: the raster mask equals `shapely.contains_xy(poly, xs+0.5, ys+0.5)` at pixel centres, EXCLUDING pixels whose centre lies exactly on an edge (use integer vertices doubled and shifted by 0.25 so no centre is on a boundary, or filter by `poly.boundary.distance(pt) < 1e-9`). `shapely.is_valid(poly)` of a legal polygon is True; area by shoelace equals `poly.area` exactly (integers).
- Why: independent oracle (GEOS) for the validity and raster rules; `validity_oracle.py` and `raster_oracle.py` were written by the same team as the checker. Boundary pixels are where a convention mismatch (half-open vs closed) hides.

---------------------------------------------------------------------------------------------------------------

## J. Wave 2 and Wave 3 (post-dependency)

### T50 `polygonizer_moved_to_oracle_exactly` (regression + property)
- Expected: `src/objdel/checker/polygonize.py` absent; `not hasattr(objdel.checker, "polygonize")`; `tests/oracles/polygonize_oracle.py` imports nothing from `objdel.checker.polygonize` (AST) and `polygonize(mask, offset=(0,0))` works.
- One-time gate: `ast.dump` of the oracle equals `ast.dump` of `git show pre-clean-2026-11-08:src/objdel/checker/polygonize.py` apart from the module docstring (a move, not a rewrite).
- Differential: for 300 seeded random masks (sizes 1..40, density 0.1..0.9) outputs equal a pre-CLEAN hash (SHA-256 of `repr` of the rings), pre-CLEAN literal.
- Properties: empty mask -> `[]`; 1x1 `True` -> one ring with 4 vertices; full HxW -> one rectangle ring with 4 vertices at `(0,0),(W,0),(W,H),(0,H)` in some fixed orientation; mask with a hole -> 2 rings, signed areas sum to the mask area (`np.count_nonzero`); `offset=(7,-3)` equals unshifted output translated; diagonal-touching pixels (`[[1,0],[0,1]]`) follow today's connectivity (pin pre-CLEAN result, 2 rings); mask touching all four borders; non-bool input (`uint8` 0/1, int64) equals bool input; 512x512 random mask completes in < 5 s [estimate; measure before].
- Why: after the move nothing in `src` tests the oracle; this keeps it equivalent to what production last used, and `-0` / ring orientation flips are classic differences.

### T51 `checker_modules_removal_and_no_dangling_imports` (property, B21)
- Expected: `certify.py`, `bitquad.py` absent only after phase D's device check lands (marker/condition on `objdel.checker` exporting the device check; else the test is an explicit `xfail(strict=True)`, not a silent skip); `not hasattr(objdel.checker, name)` for `polygonize, bitquads, certify`; AST scan of `src, tests, experiments` (including `experiments/phase1/bench/xarch_hashes.py`, `experiments/cuda_checker/`) finds no `from objdel.checker import <those>`, `import objdel.checker.certify`, `from objdel.checker.bitquad import`, relative `from .certify import` in `src/objdel/checker`, and no `importlib.import_module("objdel.checker.certify")` string.
- Subprocess smoke: `python -I -c "import objdel.checker, experiments.phase1.bench.xarch_hashes"` exits 0 (no dangling import at runtime; `__all__` in `checker/__init__.py` lists only existing names, verified by `all(hasattr(m, n) for n in m.__all__)`).
- Why: `from X import *`-style `__all__` entries pointing at deleted names only fail at import of the package.

### T52 `lfm_and_full_finetune_gone` (property, B22 Wave 3)
- `set(MODELS) == {"qwen3-vl-2b"}`; no `lfm` (T26 needle) in `src/`, `experiments/phase2/*.py`, `tests/`; AST: no `param_groups` identifier and no string constant `"float32"`/`torch.float32` branch in `train/configs.py` or its phase E home; `grpo_config`/`sft_config` never read `dtype()` for a full-FT path; `experiments/phase2/configs/r1_lfm.yaml` absent and no reference to it (T41).
- Slow smoke (marker `slow`): tiny random-weight Qwen3-VL built from a config (2 layers, hidden 32, vocab 1000), `max_completion_length` cap 8, 2 steps on CPU: finishes < 120 s [estimate], loss is finite, no download (`HF_HUB_OFFLINE=1`, assert no `socket.connect`).
- Why: dead branches survive as "supporting code" unless a negative test names them.

### T53 `cpu_pool_gone_and_metrics_on_gpu_equal_cpu` (Wave 3: B23; GPU + unit)
- Unit: `select_ckpt`, `g_probe`, `final_eval` `main(["--procs","2", ...])` -> `SystemExit(2)`; `serve` has no `fork_pool`/`pool_size` (`hasattr`); AST finds no `multiprocessing` import in `experiments/phase2/*.py`.
- GPU (marker `gpu`, skipped without CUDA, run in R-test rental only): for 200 seeded (pred polygon set, reader masks) pairs, GPU metrics (IoU, F, area, betti) vs the CPU numpy oracle `tests/oracles/score_numpy.py`: integer quantities (intersection, union, area counts, Betti numbers) bit-equal (`==`); float ratios identical when computed from the same integer pair in float64 on host (assert `abs(a - b) <= 1 ulp`, `np.nextafter` check) [estimate: ratio in float32 may differ; the spec asks for exact scores, so require equality of the int pair first].
- Failure mode: a device-side OOM or kernel error surfaces as an exception (not a silent zero score): feed a 0-sized mask and an all-ones 512x512 mask; expected scores equal oracle (0 vs 1; IoU of two empty sets pinned to pre-CLEAN value).
- Why: moving metrics onto the GPU and deleting the pool leaves no CPU fallback; bit-exact integers are the only check that does not need a tolerance.

---------------------------------------------------------------------------------------------------------------

## K. Cross-cutting gates (one-time scripts for the owner run; C3 forbids agents from changing git)

### T54 `commit_hygiene_read_only` (e2e, read-only git)
- `git diff --cached --name-only` before the owner's commit: no `README.md`, nothing under `docs/` except as allowed by CLEAN tasks (C9; note gap 4), nothing under `data/`, `.claude/`; `git log -1 --format=%B` has no `Co-Authored-By`, no `Generated with`, no `claude.ai/code/session` line (project rule + C7); `git ls-files | grep -iE 'riga.*\.(png|tif|npy|npz)'` empty and no tracked `.npz`/`.npy`/image over 1 MB (C8).
- Why: three owner rules in one cheap script; nothing in `make test` guards the commit message.

### T55 `full_suite_accounting` (regression, G1/G2/G8)
- `make lint` exit 0; `make test` 0 failed; `passed == 1447 - removed + new` where `removed` = node ids listed in archive READMEs (T30) and `new` = collected new ids; `skipped == 6` unless a documented reason; deselected 31 unchanged; `make dead` exit 0 and 0 output lines; wall time of `make test` <= 1.25 x the pre-CLEAN time [estimate: capture pre-CLEAN time; hygiene tests add vulture ~ seconds].
- Why: guards "green by deletion" and test-time creep.

---------------------------------------------------------------------------------------------------------------

## Pytest sketches (the five most important)

Ranking: S1 bundle import (T09, a broken bundle wastes rental money), S2 `--procs` rejection (T13), S3 freeze ordering (T16), S4 `boxes_meet` oracle (T19), S5 golden privacy (T34). Others (T04 existence filter, T43 vulture canary) are close behind.

### S1 (T09)
```python
import json, shutil, subprocess, sys, importlib.util
from pathlib import Path
from experiments.phase2.vast import policy as P

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = """
import sys, json
dest = sys.argv[1]
sys.path[:0] = [dest + "/src", dest]
import objdel, objdel._gpucheck
import experiments.phase1.reward_gates, tests.oracles.reward_harness
import experiments.phase2.train_rl, experiments.phase2.g_probe
import experiments.phase2.select_ckpt, experiments.phase2.final_eval
mods = {n: getattr(m, "__file__", None) for n, m in sys.modules.items()
        if n.split(".")[0] in {"objdel", "experiments", "tests"}}
print(json.dumps(mods))
"""

def build_copy(dest: Path, files: list[str]) -> None:
    for f in files:
        (dest / f).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / f, dest / f)
    ext = Path(importlib.util.find_spec("objdel._gpucheck").origin)
    shutil.copy2(ext, dest / "src" / "objdel" / ext.name)

def foreign_modules(dest: Path) -> list[str]:
    out = subprocess.run([sys.executable, "-I", "-c", SCRIPT, str(dest)],
                         capture_output=True, text=True, timeout=120, cwd=dest.parent)
    assert out.returncode == 0, out.stderr
    mods = json.loads(out.stdout.splitlines()[-1])
    assert len(mods) >= 12
    return [n for n, f in mods.items() if f is None or not f.startswith(str(dest))]

def test_bundle_copy_imports_from_the_copy(tmp_path):
    files = P.code_files()
    assert not [f for f in files if f.startswith((".claude/", "archive/", "data/", "docs/"))]
    dest = tmp_path / "bundle"
    build_copy(dest, files)
    assert foreign_modules(dest) == []

def test_bundle_import_canary_detects_a_missing_file(tmp_path):
    files = [f for f in P.code_files() if f != "experiments/phase1/reward_gates.py"]
    dest = tmp_path / "bundle"
    build_copy(dest, files)
    try:
        leaked = foreign_modules(dest)
    except AssertionError:      # ImportError in the subprocess: also a detection
        return
    assert leaked, "a file missing from the copy was served from outside it"
```

### S2 (T13)
```python
import pytest
from experiments.phase2 import train_rl
from objdel.train import serve

REQUIRED = ["--model", "qwen3-vl-2b", "--ckpt", "c", "--train", "t", "--out", "o", "--steps", "1"]

@pytest.mark.parametrize("extra", [["--procs", "2"], ["--procs=2"], ["--proc", "2"], ["--procs"], ["-procs", "2"]])
def test_procs_is_an_argparse_error_before_any_work(extra, monkeypatch, capsys):
    called = []
    monkeypatch.setattr(serve, "preload", lambda *a, **k: called.append("preload"))
    monkeypatch.setattr(serve, "fork_pool", lambda *a, **k: called.append("pool"), raising=False)
    monkeypatch.setattr("gc.freeze", lambda: called.append("freeze"))
    with pytest.raises(SystemExit) as e:
        train_rl.main(REQUIRED + extra)
    assert e.value.code == 2
    assert called == []
    assert "--procs" not in capsys.readouterr().out

def test_help_lists_the_other_options(capsys):
    with pytest.raises(SystemExit) as e:
        train_rl.main(["--help"])
    out = capsys.readouterr().out
    assert e.value.code == 0 and "--procs" not in out
    for opt in ("--model", "--ckpt", "--train", "--train-views", "--out", "--steps", "--save-steps"):
        assert opt in out
```

### S3 (T16)
```python
def test_gc_freeze_once_after_the_last_preload(monkeypatch, tmp_path, stub_train_rl):
    # stub_train_rl: fixture from tests/test_p2_train_rl.py that stubs model, trainer, CUDA and returns main's argv
    calls: list[str] = []
    monkeypatch.setattr(serve, "preload", lambda *a, **k: calls.append("preload") or {})
    monkeypatch.setattr("gc.freeze", lambda: calls.append("freeze"))
    argv = stub_train_rl(tmp_path, train_views=["v1", "v2"])
    assert train_rl.main(argv) == 0
    assert calls.count("preload") >= 2
    assert calls.count("freeze") == 1
    last_preload = len(calls) - 1 - calls[::-1].index("preload")
    assert calls.index("freeze") > last_preload
```

### S4 (T19)
```python
import numpy as np
from hypothesis import given, strategies as st
from objdel.geometry.box import boxes_meet

G = 12
coord = st.integers(0, G)
box = st.tuples(coord, coord, coord, coord).filter(lambda b: b[0] < b[2] and b[1] < b[3])

def raster(b):
    g = np.zeros((G, G), bool)
    g[b[1]:b[3], b[0]:b[2]] = True
    return g

@given(box, st.lists(box, max_size=4))
def test_boxes_meet_equals_shared_pixels(b, others):
    expected = any((raster(b) & raster(o)).any() for o in others)
    assert boxes_meet(b, others) == expected
    assert boxes_meet(b, iter(others)) == expected
    assert boxes_meet(tuple(np.int64(v) for v in b), [tuple(np.int64(v) for v in o) for o in others]) == expected
```

### S5 (T34)
```python
import re
from pathlib import Path
import numpy as np
import pytest

GOLDENS = Path(__file__).parent / "goldens"
MAX_ELEMENTS = 4096  # a 512x512 mask is 262144; a 64x64 crop is already 4096

def check_goldens(directory: Path) -> list[str]:
    bad = []
    for p in sorted(directory.glob("*.npz")):
        if p.name == "polycore.npz":
            continue
        with np.load(p, allow_pickle=False) as z:
            for k in z.files:
                a = z[k]
                if a.ndim > 1 or a.dtype == bool or a.size > MAX_ELEMENTS:
                    bad.append(f"{p.name}:{k} ndim={a.ndim} dtype={a.dtype} size={a.size}")
                if re.match(r"(?i)^riga/.+/.+", k):
                    bad.append(f"{p.name}:{k} per-case RIGA key")
                if k.lower().startswith("riga") and a.size != 1:
                    bad.append(f"{p.name}:{k} RIGA array is not an aggregate")
    return bad

def test_goldens_hold_no_masks():
    assert check_goldens(GOLDENS) == []

@pytest.mark.parametrize("name,arr", [
    ("2d", np.zeros((4, 4), np.uint8)), ("bool", np.zeros(8, bool)),
    ("flat_mask", np.zeros(100_000, np.uint8)), ("riga/c1/r1", np.zeros(1, np.int64)),
    ("riga_vec", np.arange(50)),
])
def test_checker_rejects_masks_and_per_case_riga(tmp_path, name, arr):
    np.savez(tmp_path / "x.npz", **{name: arr})
    assert check_goldens(tmp_path)

def test_checker_accepts_a_clean_golden(tmp_path):
    np.savez(tmp_path / "x.npz", **{"lidc/a/1": np.arange(30), "riga_total_vertices": np.array([7]),
                                    "riga_case_count": np.array([3]), "riga_hash": np.array("ab12")})
    assert check_goldens(tmp_path) == []
```
