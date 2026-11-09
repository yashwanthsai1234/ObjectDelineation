# Final review round 9: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. The built doc is
`docs/implementation2.0/Implementation.md` (1,632,532 bytes, 16,124 lines [measured: `ls -l`, `wc -l`, 2026-11-09]). I
read the sources in `impl/phases/` by heading. Items already fixed in rounds 1-8 (review-log.md) are not repeated.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes an NVIDIA product claim. I did no catalog lookup and installed nothing. |
| engineering:code-review | I opened every cited code line: `policy.py:300-330, 420-500` (`cmd_launch`, `chain`, `vision_refs`, `cmd_run`, `job`), `train_rl.py:15-70`, `g_probe.py:125-165`, `train_sft.py:76-112`, `configs.py:60-200`, `dataset.py:52-158`, `secrets.py:14`, `preflight.py:24-36`, `.github/workflows/tests.yml:34-44`, `Makefile:20-30`, `pyproject.toml:40-80`. |
| engineering:debug | I traced one RL run on the post-E6 tree: `run.sh` exports the session clock, `chain()` gives `to_argv`, and `train_rl` parses with `cli_parse`. A probe showed that the clock never reaches the stage (item 2). |
| engineering:system-design | I walked the "raise after E0" branch from section 1 through E6.8 and E6.8b. I checked which reader takes `c_launch.yaml`, `b3_launch.yaml`, `c_speed.yaml` and `b3_v5.yaml` after E6 (item 1). |
| engineering:architecture | I checked the layer seams E3-E8 build against the CI coverage gate. I measured the coverage of the infra code that E8.1 moves into `src/objdel` (item 3). |
| engineering:documentation | Each blocking or major item quotes the source file, the heading and the code line or measurement. |
| huggingface-skills:huggingface-best | No model choice is in scope. Nothing changed. |
| system-design-skills:system-design | I used its "make interfaces concrete" check on the stage API (`RlStage` fields against the env-var flag defaults, item 2). |
| system-design-skills:service-decomposition | I checked E7 (`verify.api` in process, HTTP as an optional extra). It stays a modular monolith with one owner per piece of data. No finding. |
| mattpocock-skills:codebase-design | I applied the deletion test and the one-adapter rule to `policy.sampler`, `rl.trl_adapter` and `core.spec`. The seams hold. Item 6 is a list owner question (`LIVE_CLIS`). |
| WebSearch/WebFetch | import-linter layers docs: `exhaustive` and `exhaustive_ignores` exist as E0 and E8.4 use them [web: https://import-linter.readthedocs.io/en/stable/contract_types/layers/]. The pydantic-settings env behaviour comes from a local probe, not from a web page (item 2). |

## Counts

blocking 0, major 3, minor 6.

## Items

### 1. major. A raise after E6 merges `after-raise`, and that branch holds plan and launch-config files that no post-E6 reader takes. The doc calls this merge a mechanical path map.

Sources: 00_front, E (B and C own the files).

- 00_front.md "## 1", "End state under (b')": "If E0 has started, the owner rebases the branch onto post-E paths, with
  E's move tables as the path map (E R-4). This is the one exception to order rule 1." Order rule 1: "The one
  exception is the `after-raise` branch merged after E0 ... it maps its paths through E's move tables". The same
  section puts "`c_speed_check`, `c_speed.yaml`, `c_launch.yaml`, ... B2.2-B2.4 ..., B3.1-B3.10 (... configs)" on
  that branch.
- 00_front.md "## 0. Words", launch config: "A `RunArgs` YAML ... V2: `c_launch.yaml`; V3: `b3_launch.yaml`". 6.3
  V2/V3 launch lines: `make preflight CONFIG=experiments/phase2/configs/c_launch.yaml`, `push_hf.sh <run>
  b3_launch.yaml`.
- E.md "### E.4 Rules", R-4: the `after-raise` branch "use[s] the move tables below as the path map (a mechanical
  change)". E.md "## E.8", dependencies: "B2 ... and B3 ... land before E on pre-E paths (front order rule 1; no
  post-E exception)". This contradicts front order rule 1.
- E.md "### E6" step 8: "The E6.8 commit deletes `preflight.load`, the `rows_plan` reader and the plan-file
  `TypeError` path of B1.3". Step 8b: "After step 8 no reader takes a RunArgs file". I29: "`chain()` reads only
  `specs/*.yaml`: `chain()` on a path outside `specs/` raises". Step 2 names `b2_fp16`, `b3_v5` and `c_speed` as
  plans that "are marked live or finished by the owner on the E0 commit". Under (b') those plans are not on `main`
  at E0 (front: they wait on `after-raise`).
- E.md "### E2" step 1: "E2 takes the fields from the CLI on the E0 commit". At E6.4, `train_rl` stops using
  argparse (`cli_parse(RlStage, argv)`, E6 intro). B2.2's `--numerics`, B3's `--reward v5` and C's
  `--check-rope-every`/`--record-batches` live on the branch. So no `RlStage` field exists for them, and with
  `extra="forbid"` a spec that sets one fails to load.
- Code: `cmd_launch` calls `preflight.refusal(Path(a.config))` before every rental [code:
  experiments/phase2/vast/policy.py:315]. After E6.8 that call is `load_spec`, so `c_launch.yaml` is refused.

Problem: an engineer who follows front order rule 1 and E R-4 after a raise that comes after E6 cannot launch V2 or
V3. Preflight refuses both launch configs. `chain()` refuses both plan files. The branch's `train_rl`, `g_probe` and
`final_eval` flag edits target argparse tables that E6.4 removed. The move tables map module paths only. They do not
say how to port a `steps:` plan or a RunArgs launch config to a `RunSpec`, or how to add fields for flags that land
after E0. It is caught at $0 (preflight refuses), and only on a non-default path. That is why I rate it major, not
blocking. But the doc gives no way forward, and two texts claim the rule ("mechanical" against "chain reads only
specs").

Fix:
- E R-4 and front order rule 1, one added sentence each: "If E6 has landed, the merge is not only a path map. The
  same merge commit (i) ports each branch plan (`c_speed.yaml`, `b2_fp16.yaml`, `b3_v5.yaml`) to `specs/<name>.yaml`
  (`CmdStage` for `bash`/nsys steps, `RlStage`/`EvalStage` for the rest). (ii) It folds each launch config into its
  spec: `bundle_view` = its `view`, and the rl stage `out` gives the push revision (E6.8b). (iii) It adds each new
  flag as an `RlStage`/`TrainerFlags` field, with its `FLAG_ORDER` entry and a golden-5 row (E2.1 rule). (iv) It
  re-points C's U-CF4 and B's I22b to the specs."
- Front 6.3 V2/V3 launch lines: add "after E6: `CONFIG=specs/c_speed.yaml` / `specs/b3_v5.yaml`; `run.sh <run>
  specs/<name>.yaml`".
- E.8 dependencies: delete "no post-E exception" and point to front order rule 1.

### 2. major. After E6.4, `train_rl` loses the env-var defaults of `--session-start`, `--rate-usd-per-h` and `--fixed-usd`. The ported live plan `b1_rows` then has a budget cap with no session clock.

Sources: E (B owns B1.5).

- B.md "B1.5" (V1 step 3 and the run-guard bullet): "`train_rl --budget-usd X --session-start <unix epoch>
  --rate-usd-per-h Y --fixed-usd Z`. `run.sh` exports ... `OBJDEL_SESSION_START` = `started`,
  `OBJDEL_RATE_USD_PER_H` = ..., `OBJDEL_FIXED_USD` = .... `train_rl` uses them as the flag defaults. `--budget-usd`
  without a session start (flag or environment) is an argparse error." `b1_rows.yaml` step 3 writes
  `--budget-usd 2.05` and no `--session-start`.
- E.md "### E6", intro: "One parser form for every CLI: `CliApp.run(Stage, cli_args=argv,
  cli_settings_source=CliSettingsSource(Stage, cli_kebab_case=True, cli_implicit_flags=True))`". Stages are
  `BaseModel` (E2 step 1: `class StageBase(BaseModel, frozen=True, extra="forbid")`). E2 step 1: "every flag that
  A-D add to `train_rl` ... becomes a `TrainerFlags` or `RlStage` field". E6 step 2 ports `b1_rows.yaml` as the live
  RL plan. A grep of E.md for `BudgetStop`, `budget`, `OBJDEL_SESSION_START` and `OBJDEL_RATE` finds nothing.
- Measurement: `CliApp.run` with that exact parser form on a frozen `BaseModel` with `budget_usd` and `session_start`
  fields, `OBJDEL_SESSION_START=123.0` and `SESSION_START=123.0` set, and argv `--budget-usd 2.05`, gives
  `budget_usd=2.05 session_start=None` [measured: pydantic-settings 2.16.2 (the box-lock version),
  `impl/reviews/r9arch_probe/env_probe.py`, 2026-11-09]. A `BaseModel` stage reads no environment.
- E handled the same pattern for SFT: "Fields left at their defaults (`data_root`, `run_dir`) are not emitted, so the
  box keeps today's defaults (`OBJDEL_DATA_ROOT`, ..., `OBJDEL_RUN_DIR`)". `RunArgs` reads them through
  `default_factory=lambda: os.environ.get(...)` [code: experiments/phase2/train_sft.py:81-85]. No such rule exists
  for the B1.5 flags.

Problem: on the post-E6 tree, every RL run from a spec (the ported `specs/b1_rows.yaml`, B-S2/B3 after a raise) gets
`session_start=None`. The engineer must guess. (a) If B1.5's "budget without start is an error" becomes an `RlStage`
validator, then `load_spec(specs/b1_rows.yaml)` raises, E-G9 ("every file in `specs/` loads") fails, and golden 5
cannot hold. (b) If it is dropped, `BudgetStop` is not built, and the session-spend guard (front 6.4 rule 5) is gone
on a paid run. (c) If the env is read in a validator, `same_stage` (spec on the Mac against CLI on the box) and
`to_argv` drift. Each choice changes the guard or the gate, so this is major.

Fix (E2 step 1 and E6 step 4, plus a U-test):
- "`RlStage.session_start`, `rate_usd_per_h` and `fixed_usd` take `Field(default_factory=...)` that reads
  `OBJDEL_SESSION_START`, `OBJDEL_RATE_USD_PER_H` and `OBJDEL_FIXED_USD` (`None` when unset), as `SftStage.data_root`
  does. A default is not in `model_fields_set`, so `to_argv` and golden 5 do not change."
- "The rule 'budget needs a session start' moves from argparse to `rl.run.run_rl` (box runtime). It is not a model
  validator, so the Mac `load_spec` passes."
- Test: with the three env vars set, `cli_parse(RlStage, to_argv(b1_rows rl stage)[2:])` gives the env values, and
  `run_rl` builds `BudgetStop` with them. With them unset and `budget_usd` set, `run_rl` raises before
  `build_trainer`.

### 3. major. E8.1 moves about 740 statements of half-covered infra code into `src/objdel`, and the CI coverage gate (`--cov=objdel --cov-fail-under=90`) then fails. No E step handles it.

Sources: E (DOC8 states the coverage source).

- E.md "### E8" step 1: `src/objdel/infra/`: "`vast/{offers,ledger,remote,chain,__main__}.py` (split of
  `experiments/phase2/vast/policy.py`, 628 lines), `bundle.py`, `preflight.py`, `hf.py` (hf_push + hf_pull)". The
  thin-CLI rule moves the library bodies.
- E.md "### E7" step 5 knows the gate: "the CI train job fails below 90 % line coverage [code: tests.yml:44]". It
  handles it only for `verify/http.py`. DOC.md "DOC8": "Coverage source is `objdel` only."
- Code: `.github/workflows/tests.yml:44` runs `pytest -q -m "not slow" --cov=objdel ... --cov-fail-under=90`.
  `Makefile:25-26` `cov` does the same (D-G12 binds it).
- Measurement: the full fast suite on be52b1c with coverage on `objdel` and on the files E8.1 moves gives: `src/objdel`
  3,866 statements, 220 missed, 94.3 %. `vast/policy.py` 379/191 missed (50 %), `preflight.py` 179/92 (49 %),
  `vast/bundle.py` 91/13 (86 %), `hf_push.py` 79/5, `hf_pull.py` 14/1. Combined: 4,608 statements, 522 missed,
  88.7 % [measured: `pytest -q -m "not slow" --cov=objdel --cov=experiments.phase2.vast --cov=experiments.phase2.preflight
  --cov=experiments.phase2.hf_push --cov=experiments.phase2.hf_pull`, Mac, 1447 passed, 574 s,
  `impl/reviews/r9arch_cov.txt`, 2026-11-09]. CI skips the `data` tests, so the CI number is not higher
  [estimate]. A-D change both counts before E0 [estimate].

Problem: the E8.1 move commit follows the doc exactly and turns the CI `train` job red (88.7 % < 90 % on today's
code). `make cov` and D-G12 fail too. The commit rule is "`make test` and `make lint` must pass", so the engineer
must guess: lower the gate, omit `objdel/infra` from coverage, or write tests for the Vast code. Each changes a repo
gate.

Fix (E8 step 1, new bullet "Coverage"): "On the E0 commit, record the coverage of the files E8.1 moves
(`make cov` with `--cov=experiments.phase2.vast --cov=experiments.phase2.preflight --cov=experiments.phase2.hf_push
--cov=experiments.phase2.hf_pull`) in `runs/E/baseline.txt`. If `objdel` plus those files is below 90 %, the change
commit before the E8.1 move adds unit tests for `ledger`, `offers` and `preflight.local`/`CHECKS` (fakes, no ssh),
until the combined figure is >= 90 %. Do not add `objdel/infra` to coverage `omit`: the spend guard (`spend`,
`HARD_STOP`) lives there." Add a gate row to E-G15 or E-G1: "`make cov` passes on the E8.1 commit".

### 4. minor. E.8 V4 step 3 still says "pool size 1", but no pool exists after D7 step 6 and E3b

E.md "## E.8": "`final_eval` with K 8 and R5's seed, on the CUDA build with pool size 1 (as A3)". E.3.3 row: "D7
step 6 ... deletes the `--procs` flags, so no `fork_pool` caller is left on the E0 commit ... E3b archives the pool
names". Fix: delete "with pool size 1 (as A3)".

### 5. minor. The `vision_refs` guard message names a plan that E6.6 archives

Code: `cmd_run` writes "missing G-vision reference {r}: run.sh {run} experiments/phase2/configs/r0.yaml first"
[code: experiments/phase2/vast/policy.py:473-474]. E6.6 archives `r0.yaml`, and E6.8b ports `vision_refs` without
the message. Fix (E6.8b vision-guard bullet): "the refusal text names the spec stage that writes the reference (the
SFT stage with `vision_reference_only`), not `r0.yaml`. Test: the message holds no `experiments/phase2/configs/`
path."

### 6. minor. `LIVE_CLIS` has two owners and an incomplete list, while the E8 gate requires every `-m experiments.*` string in it

- CLEAN.md "### CLEAN1" step 5 creates `P.LIVE_CLIS` with 4 box CLIs in wave 1, before E. E.md "### E8" step 1
  calls it "a new constant `LIVE_CLIS`". E8 step 2 lists 14 CLIs. The E8 gate: "every `-m experiments.*` module
  string in `src/`, `specs/`, `Makefile` and live `experiments/` resolves ... and is in `LIVE_CLIS`".
- The A-D and DOC texts add CLIs that E8.2 does not list: `experiments.phase2.headline` (DOC6 Makefile target),
  `ceiling` (A2), `check_series_z` (B1.7), `outline_lengths` (A8), `fp16_check`, `b2_attribution`,
  `perf.bench_checker` (D) [measured: grep of the phase sources].
- CLEAN1's `test_bundle_imports_standalone` imports every `LIVE_CLIS` module from a `code_files()` copy. Growing the
  tuple to the Mac data builders makes that box-import test cover Mac-only CLIs.

Fix (E8 step 1): "E8.1 extends CLEAN1's `LIVE_CLIS` (it does not create it). `LIVE_CLIS` = every module named by
`-m experiments.*` in `Makefile`, `specs/` and live `experiments/`, found on the E0 commit by grep. A second tuple
`BOX_CLIS` (CLEAN1's four, plus any `-m` target of a spec stage) feeds `test_bundle_imports_standalone`."

### 7. minor. `VllmSampler(gpu_memory_utilization, max_model_len)` has no stated source for `final_eval`

E.md "### E4" step 1: `class VllmSampler: def __init__(self, ckpt, model, gpu_memory_utilization, max_model_len)`.
`EvalStage` has neither field (E6.1). The TOST trainer arm ran with `vllm_max_model_length: 4096` and
`vllm_gpu_memory_utilization: 0.35` [code: src/objdel/train/configs.py:192-193], and `sample_vllm.py` uses 0.85 and
its own `--max-model-len` [code: sample_vllm.py:15, 94-95]. Fix (E4 step 3): "`final_eval` builds
`VllmSampler(..., gpu_memory_utilization=0.85, max_model_len=4096)`: the trainer arm's context length, and the
standalone sampler's memory share (no trainer on the card). Constants live in `policy/sampler.py`."

### 8. minor. E.3.2 data-flow diagram still differs from the E2 signatures

E.md "### E.3.2": `rl.run.run_rl(stage: RlStage) -> rl.build.build_trainer(model, ckpt, rows, out, roots, flags)`.
E2 step 2: `build_trainer(model, ckpt, rows, out, roots, steps, flags, **over)`. E2 step 3: `run_rl(stage: RlStage,
rd: Path)`. The `--spec ... --stage` box line was round 8 minor 11 and is not repeated here. Fix: add `steps` and
`rd` to the diagram.

### 9. minor. E2.1's mixin rule for A6's `--profile-window` misses one U52 detail

E.md "### E2" step 1: "An A-D flag that `build_trainer` reads (for example A6's `--profile-window` ...) goes into
`TrainerFlags`, so the other CLI gains it too and it joins that CLI's added set". U52 then lists the `added` sets "on
be52b1c plus A6" as `{compile, plain}` and `{lr}`. If `--profile-window` is on `train_rl` only (A6 step 3), the
`g_probe` set must hold `profile_window`. Fix (U52): write the E0-commit sets explicitly once A6 lands, for example
"`g_probe` added = `{lr, profile_window}` if A6 adds the flag to `train_rl` only".

## What looks good

- One rule for each seam: `geom.engine` is the only `_gpucheck` importer, `rl.trl_adapter` holds every private TRL
  name, and `core.provenance.repo_root()` is the one root (front 2.2/2.3 now agree with E).
- The round 8 fixes hold. V4 lives in `specs/` (I29). The pool has one owner (D7 step 6, E3b, CLEAN16 check). U52 now
  has `dropped` and `added` sets.
- DOC4 check 16 holds on the current build: 0 tables with a "$" header outside front section 6 [measured: probe of
  `Implementation.md`, 2026-11-09]. `Implementation.md` (1.63 MB) is under the 2 MB DOC-G13 limit. The 1,480,096
  figure in DOC1 step 2 and DOC-G13 is stale (round 8 minor 9).
