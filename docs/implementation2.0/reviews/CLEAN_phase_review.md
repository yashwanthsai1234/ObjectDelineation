# Review of Phase CLEAN (docs/implementation2.0/CLEAN.md)

Date: 2026-11-08. Reviewer: fresh critical reviewer. No REPO file edited, no git write, no GPU, no money spent.
Every claim below was checked against the code or a command run on 2026-11-08. Tags: [code: file:line],
[measured: command], [web: url], [pypi: package].

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Loaded. CLEAN picks no NVIDIA product. I checked the one NVIDIA claim (cuCIM 26.10): the docs page returns HTTP 200 and `cucim-cu12`/`cu13` are 26.10.0 on PyPI. No finding. |
| engineering:code-review | Correctness check on each "delete this line" step. It found item 2: `rewards.py:107` is live code that FastGRPO reads, not a docstring. |
| engineering:debug | Reproduce-then-isolate on the CLEAN0 command. I ran it and found the `tracked` count is wrong (1,038 instead of 802) because paths contain spaces (item 5). |
| engineering:system-design | Checked dependencies between the waves and phases B/C. It found CLEAN2/6/7 edit the same files as B3 and C (item 15). |
| engineering:architecture | Layer check on `boxes_meet` (item 18) and the SkyPilot keep-list evidence (item 20). |
| engineering:documentation | Stale-path sweep: 13 docs files and about 30 src/experiments docstrings still name the moved folders (item 9). |
| huggingface-skills:huggingface-best | Loaded. CLEAN picks no model. No finding. |
| engineering:testing-strategy | Checked each archived test: did it test frozen code or live code? Three tests compare LIVE code against lab code (item 3). Also: no $0 gate proves the shipped bundle still imports (item 14). |
| engineering:tech-debt | Ranked findings by impact x risk. Data-licence and silent-regression items come first. |
| Web tools (WebSearch, WebFetch, Exa) | PyPI JSON for vulture, vastai, vastai-sdk, shapely, import-linter, pydantic (all versions in the plan confirmed). vastai 1.11.0 `requires_dist` shows exact pins on pillow and cryptography (item 4). rsync man page confirms `--files-from` does not imply `-r`. vulture README confirms the Python API. SkyPilot issue #6446 and the FAQ show Vast region filtering (item 20). |

Confirmed correct (no action): 802 tracked files, 276/130/20 per folder; 34 of 39 vendor files byte-identical
(5 rewardcore files differ); `code_files()` lists `.claude/settings.json`, the worktree folder and 10 `archive/`
files; `train_rl` creates the pool and only shuts it down; `MODELS`/`MAX_COMPLETION` hold the SmolVLM and LFM keys;
`Makefile:28` defaults to `r1_lfm.yaml`; test baseline `1447 passed, 6 skipped, 31 deselected`; `testpaths =
["tests"]` so `archive/**/tests` is not collected; data view sizes 270 M / 543 M / 155 M / 188 M.

---

## Findings

### 1. blocking: CLEAN4 golden would publish RIGA rater masks

- Location: CLEAN4, step 2 (and step 3 size rule).
- Problem: the golden stores the cropped rater masks, bit-packed, for RIGA and LIDC, in the tracked file
  `tests/goldens/targets_lab.npz`. RIGA terms forbid copying and redistribution
  [docs/research/notes/Phase 1 stack and algorithms/data_engineering.md:191]. The repo is public under MIT, and
  RIGA-derived tables were already removed from it [docs/changelog.md:41]. `tests/` also ships in every GPU bundle,
  while RIGA itself only travels age-encrypted. The current guard checks only `.csv/.json/.jsonl`
  [code: tests/test_p2_breakit_ops.py:297-299], so an `.npz` would pass unnoticed.
- Fix: store no masks. The test is already `@pytest.mark.data` and reads the masks from `data/fixtures`
  [code: tests/test_targets.py:381, tests/conftest.py:26-29]. For each case key, store only the lab outputs: for LIDC,
  the SHA-256 of the lab text and `vertices`; for RIGA, the vertex sum. Always port the `pad_div` crop to
  `tests/oracles/crops_oracle.py` and delete the "size rule" branch. Add a check to
  `test_code_snapshot_has_no_riga_derived_results`: every array in each `tests/goldens/*.npz` has `ndim <= 1` or comes
  from the polycore golden.

### 2. blocking: CLEAN6 step 3 deletes code that FastGRPO prefetch uses

- Location: CLEAN6, step 3 ("the docstring lines 41-43 and 107").
- Problem: `rewards.py:107` is code, not a docstring: `f.trl_reward = self`. `FastGRPO` finds the reward through this
  line to start background scoring [code: src/objdel/train/fast.py:633]. If you delete it, prefetch turns off without
  any error and the GPU wait grows. The tests stay green because `tests/test_train_fast.py:295,328` set
  `trl_reward` by hand.
- Fix: keep line 107. Change only its comment to "lets FastGRPO.prefetch (fast.py:633) find the reward". Add an
  assertion to an existing reward test: `all(f.trl_reward is r for f in r.funcs()[0])`. Gate: `grep -n "f.trl_reward =
  self" src/objdel/train/rewards.py` returns 1 hit.

### 3. blocking: CLEAN5 step 2 archives tests of LIVE code as if they tested "frozen analyses"

- Location: CLEAN5, step 2.
- Problem: three tests compare live production code against lab reference code:
  - `tests/test_p2_signal.py::test_live_b2_functions_and_group_scoring` tests `objdel.eval.signal` `per_item` and
    `group` [code: tests/test_p2_signal.py:91-120].
  - `test_c2_stratified_bootstrap_equals_power_c2` tests `c2_analyze.pooled_ci` [code: :189-196].
  - `tests/test_p2_probes.py::test_stratified_bootstrap_equals_power_c2` tests `objdel.eval.probes.stratified_ci`
    [code: tests/test_p2_probes.py:174-180].
  If these move to the archive, the bootstrap-CI code that phase A depends on loses its reference check. Gate G2 would
  still accept the change, because the dropped tests appear on the README list.
- Fix: use the CLEAN4 pattern. Before the move, a recorder inside `experiments/phase2/lab/` runs the lab functions on
  the same seeded inputs (rng 3, 11, 1/5) and writes `tests/goldens/p2_lab_refs.json`. That file holds
  `A.group_metrics(v) | PH.parts(v)`, the `score_view` raw columns and both `power_c2.pooled` results. Rewrite the
  three tests to compare live outputs against the golden with the same tolerances (`atol=1e-12`). Only then move the
  lab. Gate: the old and new tests both pass on the same commit.

### 4. blocking: CLEAN13 pin would downgrade Pillow in the project environment

- Location: CLEAN13, step 1 (`infra = ["vastai==1.11.0", ...]`).
- Problem: vastai 1.11.0 pins `pillow==12.2.0` and `cryptography==49.0.0` exactly [pypi: vastai 1.11.0 requires_dist,
  checked 2026-11-08]. The `.venv` has Pillow 13.0.0 [measured: `PIL.__version__`]. Pillow bicubic makes the model
  images, and the plan itself keeps it because it is bit-identical [code: src/objdel/targets/frame.py:19-21]. If you
  install the extra, the resolver downgrades Pillow, and image bytes and view hashes can change. The local CLI also
  lives outside the project venv today (`~/venvs/myenv/bin/vastai`) [measured: `which vastai`].
- Fix: do not add vastai to `pyproject.toml`. Install it as an isolated tool: `uv tool install vastai==1.11.0` (or
  `pipx install vastai==1.11.0`). Record the version in `docs/runbooks/gpu-run.md:15`. Make `preflight` refuse when
  `vastai --version` differs from the pinned string. Leave `pyyaml` in `dev` unchanged.

### 5. major: CLEAN0 baseline command counts tracked files wrongly

- Location: CLEAN0, step 2.
- Problem: `stdout.split()` splits paths that contain spaces (for example `docs/research/notes/Phase 1 stack and
  algorithms/...`). It gives `tracked = 1038`; the true count is 802 [measured: `len(t.split())` vs
  `len(t.splitlines())`]. The `git ls-files` subprocess also depends on the cwd, and the `lab` count leaves out
  `experiments/phase0` (which gate G5 counts).
- Fix: in the command, use `subprocess.run(['git','-C',ROOT,'ls-files','-z'],...).stdout.split('\0')[:-1]` with
  `ROOT = p.ROOT`. Add `'phase0': sum(x.startswith('experiments/phase0/') for x in t)`. Also record the vulture count
  (item 13) and the per-check preflight result (item 21).

### 6. major: CLEAN2 also removes `gc.freeze()` from the training process

- Location: CLEAN2, steps 1-2.
- Problem: `fork_pool` does more than fork. It imports `SCORERS` and calls `gc.freeze()` on the preloaded store
  [code: src/objdel/train/serve.py:170-172]. Without the pool, the training process loses `gc.freeze()`, so every
  full GC pass scans the whole store during the 23.4 s steps. The cost is not measured, and the plan's "no result
  changes" check does not cover speed. `docs/runbooks/troubleshooting.md:34` also still says to "Reduce `--procs`".
- Fix: replace line 50 with `gc.freeze()  # preloaded store: leave it out of the GC scans (was inside fork_pool)`
  rather than deleting it. Update the troubleshooting row to name only `select_ckpt`/`g_probe`/`final_eval`. Gate: a
  test that monkeypatches `gc.freeze` and asserts that `train_rl.main` calls it once after the last `serve.preload`.

### 7. major: CLEAN6 test list is incomplete and conflicts with CLEAN7

- Location: CLEAN6, step 4.
- Problem: tests also construct `Controller` outside the named "190-400 region". It appears at
  `tests/test_train_trl.py:481,488,511`, `tests/test_breakit_wiring.py:51` (the end-to-end reward->TRL wiring test)
  and in the `tests/test_breakit_reward.py:59` section [measured: grep]. CLEAN7 step 1 edits the same
  `test_breakit_wiring.py` test.
- Fix: build the list with `grep -rn "Controller\|audit_alarm\|alarms_for" tests` and write the node ids into the task.
  In `test_breakit_wiring.py`, keep the reward/TRL wiring assertions and remove only the Controller part, so the
  integration test survives. Order: CLEAN6 before CLEAN7.

### 8. major: CLEAN5 step 3 misses an assertion and gives the wrong line range

- Location: CLEAN5, step 3.
- Problem: `tests/test_experiment_scripts.py:52` asserts that `"census_sizerule"` is in `SCRIPTS`, so the test fails
  after the move. The `smoke_grpo` test runs from line 77 to at least line 100, not 67-90 [code].
- Fix: change line 52 to `{"eval_run", "prepare_rl", "build_data"} <= set(SCRIPTS)`. Remove the
  `test_sample_eval_*` and `test_smoke_grpo_*` functions by name, not by line range.

### 9. major: stale paths in docs and docstrings after the move

- Location: CLEAN5, steps 4-5.
- Problem: only README, `getting-started.md`, `frame.py:49` and the test oracles get updated. Thirteen files under
  `docs/` still name `experiments/phase0`, `experiments/phase1/lab` or `experiments/phase2/lab`: `contributing.md`,
  `runbooks/troubleshooting.md`, ADR 0001/0002/0003, `reference/reward.md`, `changelog.md`, and more. `ROADMAP.md`
  does too. About 30 src/experiments docstrings cite `lab/B3`, `lab/t_stack` and similar (for example
  `src/objdel/eval/probes.py:8,298,328,452,552` and `src/objdel/train/configs.py:58,78`) [measured: grep].
- Fix: add a step that runs `grep -rln "experiments/phase0\|experiments/phase[12]/lab\|lab/[A-Z][0-9]\|lab/t_" src
  experiments tests docs README.md ROADMAP.md`. Rewrite the hits to `archive/phase1_lab/...`,
  `archive/phase2_lab/...` or `docs/research/records/phase0/...`. Frozen records (`changelog.md`,
  `experiments/*/impl/`) are exempt. Add a hygiene test: every backticked repo path in `README.md` and `docs/**/*.md`
  (records exempt) exists.

### 10. major: the CLEAN5 `git add -A` sweeps in the owner's uncommitted docs work

- Location: CLEAN5, step 8.
- Problem: today `README.md` is modified and `docs/getting-started.md` (with most of `docs/`) is untracked in the
  owner's tree [gitStatus at session start]. `git add -A ... README.md docs/getting-started.md` would commit that
  unrelated work inside the cleanup commit.
- Fix: add a precondition: `git status --porcelain README.md docs/` is empty, or the owner commits the docs work first.
  Give the owner `git add -p README.md docs/getting-started.md` for those two files.

### 11. major: the archive README text and the tag are not accurate

- Location: CLEAN0, step 4; CLEAN5, step 7.
- Problem: "Retiring commit: parent tagged `archive/pre-clean-2026-11-08`" is false, because CLEAN1-4 each commit
  before CLEAN5. The tag exists only on the Mac. A tag named `archive/...` next to a folder named `archive/` invites
  confusion when someone mistypes a path.
- Fix: name the tag `pre-clean-2026-11-08`. Owner commands: `git tag pre-clean-2026-11-08` and then `git push origin
  pre-clean-2026-11-08`. README text: "Removed after tag `pre-clean-2026-11-08`. Restore: `git restore
  --source=pre-clean-2026-11-08 -- <old path>`."

### 12. major: the vulture whitelist will fail `make lint`

- Location: CLEAN12, steps 2 and 4.
- Problem: vulture whitelist lines are bare expressions such as `FastGRPO._iter_named_params`. Ruff selects `F` and
  `B` [code: pyproject.toml:88-90], so these lines raise F821 (undefined name) and B018 (useless expression) in
  `tests/`. The lint target checks `tests` [code: Makefile:32-33]. Also, `test_no_dead_code` "skips if vulture is not
  installed", and the local `.venv` has no vulture today [measured: ModuleNotFoundError]. So G9 passes silently on
  the Mac.
- Fix: add `"tests/vulture_whitelist.py" = ["F821", "B018"]` to `[tool.ruff.lint.per-file-ignores]`, or import the
  classes and write `_ = FastGRPO._iter_named_params`. Make the test fail, not skip, when vulture is missing: it is
  in the `dev` extra, which CI installs [code: .github/workflows/tests.yml:30].

### 13. major: the G9 baseline "6" uses a different scope from the gate

- Location: G9; CLEAN12, step 2.
- Problem: `v_src_withTests.txt` scanned src and tests only. The gate scans src, experiments and tests. The whitelist
  in step 2 names `same_sig`, `_iter_named_params` and `floating_point_ops`, which appear only in the no-tests run
  (`v_src_allexp_notests.txt`). Four of the six "today" items are Controller `on_*` methods that CLEAN6 archives.
- Fix: at CLEAN0, run the exact gate command (`vulture src experiments tests --min-confidence 60`) and store the
  count in `before.json`. Make the whitelist from that output with `vulture ... --make-whitelist`, after CLEAN6.

### 14. major: no $0 check that the shipped bundle still imports

- Location: CLEAN1 and CLEAN5; risk table row 3.
- Problem: the only guard against "the box misses a file" is `IMPORT_CHECK` at `sync_up`, and that runs on a rented
  GPU. CLEAN1 changes the bundle filter, and CLEAN5 moves 426 files.
- Fix: add `tests/test_repo_hygiene.py::test_bundle_imports_standalone`. It copies `policy.code_files()` into
  `tmp_path` and then runs `python -I -c "import objdel, experiments.phase2.train_rl, experiments.phase2.g_probe,
  experiments.phase2.select_ckpt, experiments.phase2.final_eval, experiments.phase1.reward_gates"` with
  `PYTHONPATH=tmp_path/src:tmp_path`, plus `pytest --collect-only -q tests` in that folder. Mark it `slow` if it takes
  more than 30 s.

### 15. major: hidden merge conflicts with phases B and C

- Location: section 8, "CLEAN does not block phases A, B or C".
- Problem: CLEAN2 (`train_rl.py`), CLEAN6 (`rewards.py`, `callbacks.py`) and CLEAN7 (`configs.py`) edit the same
  files as B2/B3 (reward v5, trainer-vLLM match) and C (pixels on GPU, outline limit, n=8 requests).
- Fix: put CLEAN2, CLEAN6 and CLEAN7 (about 4.5 h) first in Wave 1 and land them before any B3 or C branch starts. Say
  so under Dependencies. The rest of Wave 1 can run in parallel.

### 16. major: CLEAN11 import list is incomplete

- Location: CLEAN11, step 3.
- Problem: more files import from the checker: `tests/test_topology.py:11` (`bitquads`),
  `tests/test_breakit_geometry.py:17` (`objdel.checker.bitquad`) and `tests/test_cuda_checker.py:15` (`certify`).
  `checker/__init__.py:4-5,11` also re-exports `bitquads` and `certify` [code]. If you remove only `polygonize`, the
  CLEAN11 step 2 move breaks the imports.
- Fix: list all of these files. Remove `bitquads` and `certify` from `checker/__init__.py` and `__all__` when they
  move. Gate: `grep -rn "checker.bitquad\|checker.certify\|import.*\bcertify\b\|bitquads" src` returns 0.

### 17. major: CLEAN1 treats the symptom, not the cause

- Location: CLEAN1.
- Problem: `.claude/worktrees/` is untracked, so `ls-files -co --exclude-standard` lists it only because it is not
  ignored. `git status` stays noisy, and other tools (ruff, vulture, pytest rootdir scans) may walk 1.5 GB.
- Fix: also add `/.claude/worktrees/` to `.gitignore`. Keep the regex drop for `.claude/` (`settings.json` is tracked
  [measured: git ls-files .claude]) and for `archive/`.

### 18. minor: `boxes_meet` placed in a dataset loader

- Location: CLEAN9, step 1.
- Problem: a pure geometry predicate in `data/lidc.py` makes `targets/frame.py` import an LIDC-specific module (with
  sqlite3 and xml). That edge runs against phase E's layers (geom below data and targets). "2,000 random box lists"
  is not the hypothesis default (100 examples).
- Fix: put `boxes_meet` in `src/objdel/geometry/box.py` (or `geometry/shape.py`) and import it in both places. Use
  `@settings(max_examples=2000)`.

### 19. minor: CLEAN13 step 1 wording contradicts itself

- Location: CLEAN13, step 1 ("Move pyyaml from dev to infra and keep it in dev too").
- Problem and fix: the step becomes moot with item 4. Delete it.

### 20. minor: wrong SkyPilot evidence in the keep-list ADR

- Location: CLEAN14 ("SkyPilot has no documented US/CA filter").
- Problem: SkyPilot can restrict Vast regions with `resources.any_of` and entries such as `infra: vast/Florida, US,
  NA` [web: https://github.com/skypilot-org/skypilot/issues/6446; https://docs.skypilot.ai/en/stable/reference/faq.html].
- Fix: give the real reason to keep `vast/policy.py`: the repo's governance (scoped key, image allowlist, ledger,
  cap hours, destroy-and-verify, preflight stamp, North-America rule inside `policy.reasons()`), which SkyPilot does
  not provide.

### 21. minor: gate G13 can pass on FAIL = FAIL

- Location: G13; CLEAN9, step 3.
- Problem: `data/eval_views` has no `lidc-test-shifted` today [measured: ls], and `rl_shifted.yaml`'s own comments
  run preflight on `r1_qwen3vl2b.yaml` [code: rl_shifted.yaml:17-19]. "Same verdict" can mean two failures.
- Fix: at CLEAN0, store the preflight per-check result list (name -> pass/fail) and compare each check. Fix the
  `rl_shifted.yaml` comment so that it uses its own file before you make it the `Makefile` default.

### 22. minor: archived `t_stack/prod.py` will not run

- Location: CLEAN3.
- Problem: `experiments/phase1/lab/t_stack/prod.py` references `vendor/`. After the twins are deleted, the archived
  copy cannot run.
- Fix: in the README under "How to restore", say "run it from the tag, not from `archive/`".

### 23. minor: the same test is written twice

- Location: CLEAN1 step 3 and CLEAN12 `test_bundle_never_ships_agent_or_archive_files`.
- Fix: keep CLEAN1's temp-repo test and add one real-repo assertion to it. Drop the CLEAN12 copy.

### 24. minor: run order of CLEAN9 and CLEAN10

- Location: CLEAN10.
- Problem: the `Makefile` default config (`r1_lfm.yaml`) points at `views/lidc-bundle`, which CLEAN10 moves away.
- Fix: run CLEAN9 step 3 before CLEAN10. Hash the `rows.jsonl` files with `.venv/bin/python -I`.

### 25. minor: CI skips the CLEAN4 gate test

- Location: CLEAN4, step 6.
- Problem: the test is `@pytest.mark.data` and skips without `data/fixtures`, which CI does not have. It runs only on
  the Mac, where the fixtures exist [measured: ls data/fixtures].
- Fix: gate command: `.venv/bin/python -m pytest -q -rs tests/test_targets.py -k lab`. The output must show
  "1 passed" and no skip.

### 26. minor: CLEAN7 file list misses a file that its own gate greps

- Location: CLEAN7, Files and step 5.
- Problem: `experiments/phase2/sample_vllm.py:10` mentions SmolVLM, and the gate greps `experiments/phase2/*.py`.
- Fix: add the file and edit the docstring to name Qwen3-VL's stop id only.

### 27. minor: GPU versus CPU

- Location: section 3.4.
- Problem: none beyond item 6. CLEAN correctly defers the numpy polygonizer (CLEAN11 waits for phase D) and the CPU
  pool (CLEAN16 waits for D5). The CLEAN4 recorder is a one-time CPU run of frozen CPU code.
- Fix: none.
