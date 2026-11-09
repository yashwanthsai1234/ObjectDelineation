# Final review round 7: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc: `docs/implementation2.0/Implementation.md`
(built 23:07 from `docs/implementation2.0/*.md`). Repo at be52b1c. I checked every item against the round 1-6 rows in
`docs/implementation2.0/review-log.md` and `reviews/final_r6_arch*.md`. No item below repeats a fixed one. Round 6
fixes arch-1 (tuple `EvalStage`), arch-2 (`CliSuppress` defaults, `cli_parse`, `same_stage`), arch-4 (V4 launch
lines) and arch-5 (bracketed patterns in the phase pages) are in place; I re-ran the fixer's probe
(`(planning probe, not kept: r6fixE)/probe_eval.py`): all round trips print `True`.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes a new NVIDIA product claim (no CUDA, NIM or driver change in this round). I recommend no catalog skill and install none. |
| engineering:code-review | I opened every code line I quote: `policy.py:420-536` (`chain`, `job`, `cmd_push_hf`, `cmd_fetch_hf`), `train_sft.py:75-111`, `final_eval.py:44-56`, `select_ckpt.py:362-379`, `g_probe.py:131-145`, `preflight.py:26-35`, `tests/test_p2_scripts.py:340-390`, `tests/test_p2_train_rl.py`. |
| engineering:debug | I traced three paths by hand and by probe: the E2 commit on a clean Mac venv (item 1); the preflight dry run after E6.4 (item 2); the V4 chain step 2 through `cmd_run.job` (item 3). |
| engineering:system-design | I walked the V4 rental as a run (launch, fetch, chain steps 0-3, pull, teardown). Each chain step is a process that `job()` starts; that found item 3. |
| engineering:architecture | I checked the stage seam against its CLIs: E2.1's rule ("every flag ... becomes a field") against E6.1's field lists. That found item 2. I checked where the `--hash-view` mode lives once `final_eval` parses only `EvalStage` (item 4). |
| engineering:documentation | I ran the check 14 regex on the built `Implementation.md`, `review-log.md` (0 hits each) and on `docs/implementation2.0/kernels/**` (16 hits in 7 files), which check 14 also covers. That found item 5. |
| huggingface-skills:huggingface-best | This review makes no model choice, so nothing changed. |
| system-design-skills:system-design | Stress test of each chain step's failure path: a non-soft step that exits non-zero stops the chain, so a broken step 2 loses the V4 eval (item 3). |
| system-design-skills:service-decomposition | I re-confirmed the modular-monolith call (layers in one package, optional E7 HTTP face). No item. |
| mattpocock-skills:codebase-design | I tested the stage as the interface of each CLI: one model, two adapters (CLI and spec). A CLI flag with no field is a flag the interface silently drops (items 2 and 4). |
| WebFetch | pydantic-settings source (`pydantic_settings/sources/providers/cli.py`, main): a one-letter field gets one dash (`flag_prefix[: 1 if len(name) == 1 else None]`), and `CliSuppress` hides help only; the argument stays parsable. Both E6 claims hold. |
| Local probe | `(planning probe, not kept: r7arch)/probe_select.py` (pydantic 2.14.3, pydantic-settings 2.16.2, the r6 venv): E6.1's `SelectStage` with the one parser form exits 2 on `--procs 1` and on `--no-vision`; with a `procs` field it passes; a `timeout_s` field accepts `--timeout-s` and exits 2 on `--timeout_s`. |

## Counts

blocking 0, major 5, minor 4.

## Items

### 1. major. E2 makes `core/spec.py` import pydantic-settings, but the dependency is added only in E6

Sources: E.

- E.md "### E2", Libraries line: "Libraries: pydantic, base dependency as a range `pydantic>=2.14,<3` ... So: (a) add
  `import pydantic` to `IMPORT_CHECK` and to `image_check.sh`".
- E.md "### E2", step 1 code (round 6 fix): `class StageBase(...)`: `kind: CliSuppress[str]`, `name: CliSuppress[str]
  = ""`, `soft: CliSuppress[bool] = False`. `CliSuppress` comes from `pydantic_settings`.
- E.md "### E6", Libraries line: "pydantic-settings `>=2.16,<3` (2.16.2 in the box lock ...; base dependency beside
  pydantic ...)". E.8 owner decision 4 lists both, with no commit. E-G13 and the E.7 risk row name `pydantic` only.
- Measured: `.venv/bin/python -I -c "import pydantic_settings"` gives `ModuleNotFoundError` (the Mac venv has neither
  package today). `pyproject.toml` has no `pydantic` line.

Problem: the E2 commit adds `pydantic` only. Its new `core/spec.py` imports `pydantic_settings`, so every test that
imports `objdel.core.spec` fails at import on a fresh `pip install -e .`, and E2's gate (R-3, `make test`) is red.
The box check also misses it: `IMPORT_CHECK` tests `pydantic` only, and the phase2 image (`--no-deps` install,
`policy.py:399`) has neither.

Fix (E2 Libraries line, E.7 pydantic row, E-G13): "Libraries: pydantic `>=2.14,<3` and pydantic-settings
`>=2.16,<3`, both added to base deps in the E2 commit (E2.1 uses `CliSuppress`); E6 adds no dependency."
`IMPORT_CHECK` and `image_check.sh` get `import pydantic, pydantic_settings`. In E6's Libraries line, replace "base
dependency beside pydantic" with "already a base dependency since E2". Add to U6 (or the E2 tests): `import
objdel.core.spec` in a fresh process succeeds.

### 2. major. E6.1's `SelectStage` and `ProbeStage` drop live CLI flags, so the preflight dry run exits 2 after E6.4

Sources: E (CLEAN16 for the later deletion).

- E.md "### E2", step 1 text: "every flag that A-D add to `train_rl` or `g_probe` ... becomes a `TrainerFlags` or
  `RlStage` field ... The rule covers every CLI that a stage maps to ...: `train_rl`, `g_probe`, `final_eval`,
  `select_ckpt` and `train_sft`, each read on the E0 commit".
- E.md "### E6", step 1 code: `class SelectStage`: `run, prepared, probes, sampler, items, shifted, rescaled,
  rotated, vision (None -> --no-vision), only`; `class ProbeStage`: `model, train, probe, out, ckpt, select, items`
  plus `TrainerFlags`. Neither has `procs`; `SelectStage` also lacks `k`, `seed`, `timeout_s`, `out`, `model`, and has
  no field for `--no-vision`.
- E.md "### E6", step 4: "CLIs take `--spec PATH --stage NAME` or the old flags. The old flags go through
  `cli_parse`".
- Code: `select_ckpt.py:369-379` has `--model`, `--k`, `--seed`, `--procs`, `--timeout_s`, `--no-vision`;
  `g_probe.py:139` has `--procs`. CLEAN16 step 1 deletes `--procs` only in wave 3, after E.
- Code: the preflight dry run `DRY_RUN = "tests/test_p2_scripts.py::test_cpu_smoke_train_then_select"`
  [preflight.py:35] calls `select_ckpt ... "--procs", "1"` [tests/test_p2_scripts.py:369] and `g_probe ... "--procs",
  "1"` [:384]. E6.4 rewrites only the train call of this test.
- Measured (`(planning probe, not kept: r7arch)/probe_select.py`): `cli_parse(SelectStage, [... "--procs", "1"])` with E6.1's fields exits 2,
  "unrecognized arguments: --procs 1". `--no-vision` exits 2 the same way. With a `procs` field it passes. A
  `timeout_s` field reads `--timeout-s`; `--timeout_s` exits 2.

Problem: two parts that both claim to be the rule disagree: E2.1 says every E0 flag of these CLIs becomes a field,
and E6.1's code lists fewer fields. An engineer who copies E6.1 gets a red dry run at the E6.4 commit, and preflight
refuses every rental until it is fixed. The `--no-vision` path cannot round-trip at all.

Fix (E6 step 1 code, step 4, step 7, U22):
- `SelectStage` adds `model: str | None = None`, `k: int = 8`, `seed: int = 0`, `procs: int | None = None` (None ->
  `serve.pool_size()` in the CLI, today's default), `timeout_s: float | None = None`, `out: str | None = None`, and
  `no_vision: bool = False` (validator: exactly one of `vision` and `no_vision`; `to_argv` writes `--no-vision` as
  a bare flag). `ProbeStage` adds `procs: int | None = None`. FLAG_ORDER entries follow `select_ckpt.py:362-379`
  and `g_probe.py:131-145`.
- E6 step 4 adds: "`--timeout_s` becomes `--timeout-s` (kebab rule); no plan or test writes it today".
- E6 step 4, preflight bullet: the dry run's `select_ckpt` and `g_probe` calls keep `--procs 1` (now a field)
  until CLEAN16 deletes the field and the flag together.
- New test beside U22 (`tests/core/test_spec.py`): for each of the five CLIs, the set of argparse `dest` names on the
  E0 commit (recorded in a golden file) equals the stage fields minus `kind`, `name`, `soft`, minus a listed
  dropped set (`train` of `final_eval`, E4).

### 3. major. V4 step 2 names the box tests as `pytest -m gpu ...`, which the chain runs as `python -I pytest`

Sources: E (front 6.3 points to E.8).

- E.md "## E.8", V4 list: "2. the E.9 box test list by file (`pytest -m gpu tests/verify/test_verifier_gpu.py
  tests/infra/test_box_smoke.py tests/policy/test_sampler_box.py`)". Step 1 of the same list starts with `bash`.
  E6 step 8b: "steps 0-2 of E.8 as `CmdStage`s with the same argv". U51 checks "the E.9 box test list by file (no
  bare `pytest -m gpu`, no `test_sampler_tost.py`)" but not the leading tokens.
- Code: `cmd_run.job` [policy.py, `def job(s)`]: `cmd = shlex.join(s[soft:]) if s[soft] == "bash" else
  f"/opt/venv/bin/python -I {shlex.join(s[soft:])}"`, and the steps are joined with `&&`.
- B.md I22a pins the same kind of step for V1: "step 2 the smoke by the two node ids, starting with `-m, pytest`".

Problem: an argv that starts with `pytest` becomes `/opt/venv/bin/python -I pytest -m gpu ...`, which fails with
"can't open file '/workspace/src/pytest'". The step is not soft, so the `&&` chain stops, and step 3 (the sampler
arm, the only reason for V4) never runs. The rental is paid (about $0.71, front 6.3) and gives no E-G10 result.
Preflight does not read `CmdStage.argv` (E6.8 table), so nothing catches it before the spend.

Fix (E.8 V4 step 2, E.9 conventions "Each box test names ...", U51): write the argv:
`CmdStage(name="box_tests", argv=("-m", "pytest", "-q", "-p", "no:cacheprovider", "-m", "gpu",
"tests/verify/test_verifier_gpu.py", "tests/infra/test_box_smoke.py", "tests/policy/test_sampler_box.py"))`.
U51 adds: "every `CmdStage.argv` of `e_v4.yaml` starts with `bash` or `-m`; step 2 starts with `-m, pytest`".

### 4. major. `--hash-view` has no home once `final_eval` parses only `EvalStage`

Sources: E, A.

- E.md "### E6", step 1 code: `class EvalStage(...)`: `ckpt: str` (required), `view: ... Field(min_length=1)`,
  `out: ... Field(min_length=1)`, `extra="forbid"`; comment: "`--hash-view` stays out (a Mac-only mode, not a
  stage)". E6 step 4: "The old flags go through `cli_parse`".
- A.md "### A3": "new flag `--hash-view VIEW --out FILE`: write `{view, items, content_sha256, ...}` and exit (no
  model, no scoring)". A.md test "U54 view content hash": "`--hash-view` writes `{view, items, content_sha256,
  target_set_sha256, generator_commit, note}` and exits 0 without loading a model" (in `tests/test_final_eval.py`,
  A.md coverage table). B.md B1.0 step 1 and I36 run `final_eval --hash-view`.
- Measured (same parser form, `(planning probe, not kept: r7arch)/probe_select.py` pattern): an unknown flag exits 2 ("unrecognized
  arguments"), and a stage without `--ckpt`/`--view` exits 2 ("required").

Problem: after E6.4, `final_eval --hash-view VIEW --out FILE` exits 2, so U54 turns red at the E6.4 commit
(`make test` must pass) and the Mac view-hash tool that A3 and B1.0 use is gone. The doc says the mode "stays out"
but not where it goes. The engineer must choose between a pre-parse branch in `main`, a new CLI, or an
`EvalStage` field, and each choice changes A's and B's commands.

Fix (E6 step 4, new bullet; A's U54 and B's I36 follow): "`final_eval.main` keeps one branch before `cli_parse`:
when `argv[0] == "--hash-view"`, a two-flag argparse (`--hash-view VIEW --out FILE`) runs `view_content_sha256`
and exits; every other argv goes to `cli_parse(EvalStage, argv)`. A3's command form stays the same." Add a test to
E6.7: `final_eval.main(["--hash-view", v, "--out", f])` writes the A3 keys, and
`final_eval.main(["--hash-view", v])` exits 2. Also state that `--rescore DIR` still needs `--ckpt`, `--view` and
`--out` (all required fields), or make `ckpt` optional when `rescore` is set (validator), and check A's `--rescore`
tests against the chosen form.

### 5. major. DOC1 leaves scratchpad paths in the kernel folders, which check 14 scans

Sources: DOC.

- DOC.md "### DOC4", check 14: "0 matches of `/[p]rivate/tmp`, `/[U]sers/`, `S[P]/impl` or `S[P]/[A-Z]` in
  `doc_files` ... and in `docs/implementation2.0/**/*.md` and `docs/implementation2.0/kernels/**` text files".
- DOC.md "### DOC1", step 0: "in every copied file and in `Implementation.md`, replace each absolute scratchpad
  prefix and each short scratchpad prefix" (the kernel folders are not in the copy table); kernel bullet: "Replace
  the absolute paths in `K1/GPU_RUN.md` (2), `K2/GPU_RUN.md` (1), `K3/GPU_RUN.md` (1) and `K1/results/build_cpu.txt`
  ... [measured: `grep -c -E "/[p]rivate/tmp|/[U]sers/"`]"; re-scan: "run the check 14 scan on the built
  `Implementation.md` and `review-log.md` (and the copied files) until it prints nothing".
- Measured (`grep -rnE` with the check 14 patterns on `docs/implementation2.0/kernels`, 16 lines in 7 files):
  `K1/GPU_RUN.md:18,19,20,25,31,126` (`S[P]=...`, `REPO=...`, `K1=docs/implementation2.0/kernels/K1`, `$REPO/data/k1_bundle/lidc` twice,
  `docs/implementation2.0/kernels/K1/results/`), `K1/DESIGN.md:3,19,242,304` (`docs/implementation2.0/kernels/K1/`, `docs/implementation2.0/notes/H1_geometry_algorithms.md`,
  `docs/implementation2.0/notes/B_gpu_audit.md`), `K2/DESIGN.md:3`, `K2/GPU_RUN.md:33`, `K3/DESIGN.md:3,9` (`docs/implementation2.0/notes/doc/nv_catalog.json`),
  `K3/GPU_RUN.md:15`, `K1/results/build_cpu.txt:1`.

Problem: the kernel bullet fixes only the absolute paths (its pattern has no `S[P]/`), and the re-scan does not
include the kernel folders. After DOC1 as written, at least 11 lines remain (the `$S[P]/...` lines of `K1/GPU_RUN.md`
and the four `DESIGN.md` files), so check 14 and `make test` are red at the DOC4 commit, the same failure class as
round 6 arch-5.

Fix (DOC1 kernel bullet and re-scan sentence): replace `docs/implementation2.0/kernels/K<n>` with `$K<n>` (=
`$REPO/docs/implementation2.0/kernels/K<n>`), `$REPO/data/k1_bundle` with `$REPO/data/k1_bundle`, `docs/implementation2.0/kernels/K<n>/`
with `docs/implementation2.0/kernels/K<n>/`, `docs/implementation2.0/notes/doc/nv_catalog.json` with `docs/implementation2.0/notes/doc/nv_catalog.json`
(copied in step 0), and `S[P]/<name>.md` with `docs/implementation2.0/notes/<name>.md`; drop the `S[P]=` line. Change the
measure to `grep -rnE "/[p]rivate/tmp|/[U]sers/|S[P]/impl|S[P]/[A-Z]" docs/implementation2.0/kernels` (16 lines in 7
files on 2026-11-08). Re-scan sentence: "... on the built `Implementation.md`, `review-log.md`, the copied files
and `docs/implementation2.0/kernels/**`".

### 6. minor. DOC1's measured numbers are stale

Source: DOC (DOC1 step 1 (a), step 2). "Today (a) prints 15 lines: the absolute paths in `Implementation.md` (11)
and the kernel files (5 lines in 4 files)" and "one Markdown text file of 1,480,096 bytes". Measured now: 0 lines in
`Implementation.md`, 16 lines in 7 kernel files (item 5), and `Implementation.md` is 1,577,738 bytes. The 2 MB
limit still holds. Fix: write the counts as "re-measured at hand-off" with the commands, or update them to the
measured values above.

### 7. minor. V4's last line skips teardown when the pull fails

Source: 00_front (6.3 V4 launch lines). `pull.sh <run> && teardown.sh <run> && sweep_orphans.sh`. A failed pull
stops the line before `teardown.sh`; the watchdog cap still ends the box, so the cost is bounded by the cap. Fix:
add one sentence: "If `pull.sh` fails, retry it once; then run `teardown.sh <run>` in any case (6.4 stop rules)".

### 8. minor. Round 6 minors still open

Sources: E, 00_front. Unchanged since round 6 (`final_r6_arch.md` items 6-10):
- arch-6: `SelectStage.vision` path form and `--vision` JSON spacing (E6.1).
- arch-7: `SftStage.data_root`'s default reads `BOX_DATA`, which E2.6 puts in `sft/train.py` (E6.1, E2.6 table).
- arch-8: `vision: False` is dropped by "a `False` boolean never" (E6.1 SFT argv); write it as `--vision off`.
- arch-9: E.3.1 still says `train.py run_sft(stage)` (E.md:209) and E.3.2 still shows
  `train_rl.py --spec ... --stage train` (E.md:256), against E6.4.
- arch-10 (r5 minors): 00_front 2.2/2.3 still list `callbacks` and `infra/ paths.py` (00_front.md:278, 333, 336);
  E.8 decision 1 still reads "Keep `Controller` + `eval/audit.py`" (E.md:1719).
Apply the fixes as written in `final_r6_arch.md`.

### 9. minor. E-G13 and GP3 name the import check without pydantic-settings

Source: E (E-G13 row: "`IMPORT_CHECK` (with `pydantic`, `vllm`, `objdel.policy.sampler`)"). After item 1 the
check also needs `pydantic_settings`. Fix: add it to the E-G13 row and GP3, in the same edit as item 1.
