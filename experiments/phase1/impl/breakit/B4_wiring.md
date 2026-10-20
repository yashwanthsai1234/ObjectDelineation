# Break-it B4: end-to-end wiring, install, CI, docs

Scope: the PLAN §3.1 / §3.2 chain as code (store → targets → geometry → reward → TRL → controller → eval → gates),
the EXECUTION §3 contracts, the build (wheel and editable), the CI workflow, and README / plan references.
Reproducers: `tests/test_breakit_wiring.py` (9 strict xfails `B4-<n>`, 12 passing checks; the no-data end-to-end
test takes ~12 s). Scratch: `scratchpad/breakit/B4/` (`e2e.py`, `edges.py`, fresh venvs, build logs).
No production code, data, repo `.venv` or other test file was changed.

**Totals: 39 checks tried, 30 held, 9 broke** (2 critical, 4 medium, 3 low). Not repeated here: B2's findings (B2-1
rollback ordering is next to B4-4), I19, I20.

## Findings

| id | severity | area | one line |
|---|---|---|---|
| B4-1 | **critical** | CI | `astral-sh/setup-uv@v10` does not exist (no major tags since v8.0.0): the `python` and `gates` jobs cannot start |
| B4-2 | **critical** | CI / packaging | CI installs `.[dev]` only; `objdel.data` imports pyarrow (`[data]` extra), so 10 test modules fail at collection |
| B4-3 | medium | wiring (train → eval) | TRL's completions table has no item key and no truncation flag: training / held-out completions cannot reach `eval_run` or `hack_audit` |
| B4-4 | medium | wiring (eval → controller), deepens I21 | `Controller.alarm(metrics)` gets only TRL's scalar means; no hack-audit statistic exists during training, so `hack_audit` cannot be plugged in as the alarm |
| B4-8 | medium | CI / K13 | No CI job installs `[train]`: R2, TRL, dataset and callback tests always skip in CI; CI-like coverage is 82 % (train/ 0–4 %), no `--cov-fail-under` |
| B4-6 | medium | docs | README quick start imports `objdel.reward.AnnotatorSet` (deleted); layout, status and make targets are pre-Phase-1 |
| B4-5 | low | eval.gates | `phase2_entry` requires both datasets, but PLAN §7 drops a model that fails the precision probe from LIDC: such a model can never pass |
| B4-7 | low | experiments | `census_sizerule` and `fetch_lidc_pixels` take `sys.argv` positionally: `--help` is read as a data path / raises TypeError |
| B4-9 | low | docs | PLAN §3.1 / §3.4 name `data.ingest` / `data.export` (no such modules); M7 note says the hack audit can read the completions table (it cannot, B4-3); CI comment cites the public-repo arm64 changelog for a private repo |

### B4-1 (critical) CI: setup-uv tag does not resolve
- **Reproduction:** `test_b4_1_setup_uv_is_pinned_to_an_existing_ref`; shell:
  `gh api repos/astral-sh/setup-uv/git/ref/tags/v10` → HTTP 404, while `v10.0.0 … v10.3.0` exist
  (`gh api 'repos/astral-sh/setup-uv/git/matching-refs/tags/v'`: last major tag `v7`, then only `v8.0.0`, `v9.0.0`,
  `v10.0.0`, …).
- **Expected vs actual:** every `uses:` ref resolves. Actual: `.github/workflows/tests.yml` uses `@v10` in the `python`
  (3 OS) and `gates` jobs; GitHub fails such a job with "Unable to resolve action" before any step runs. Only
  `sanitizers` and `wheels` can run. actionlint 1.7.12 does not catch it (it checks syntax, not tag existence).
- **Root cause:** setup-uv v8.0.0 release notes: "No more major and minor tags … You won't be able to use `@v8` or
  `@v8.0` any longer … Use the immutable tag … or even better the githash"
  (https://github.com/astral-sh/setup-uv/releases/tag/v8.0.0).
- **Fix:** `astral-sh/setup-uv@v10.3.0` (or its commit SHA). Pinning by SHA is GitHub's own hardening advice
  (https://docs.github.com/en/actions/security-for-github-actions/security-guides/security-hardening-for-github-actions#using-third-party-actions).

### B4-2 (critical) CI: pyarrow missing from the CI install
- **Reproduction:** `test_b4_2_ci_installs_what_objdel_data_imports`; shell (what CI does, from a wheel):
  ```
  uv build --python 3.12 -o $S/dist                      # sdist -> wheel, objdel-0.1.0-cp312-abi3-*.whl
  uv venv --python 3.12 $S/venv && uv pip install --python $S/venv/bin/python "objdel[dev] @ file://$S/dist/<wheel>"
  cp -R tests pyproject.toml $S/wt/ && cd $S/wt && OBJDEL_FIXTURES=/x OBJDEL_STORE=/x ../venv/bin/python -m pytest -q -m "not slow"
  ```
  → `10 errors during collection` (`ModuleNotFoundError: No module named 'pyarrow'`): test_breakit_geometry,
  test_breakit_reward, test_data_{lidc,riga,splits,store}, test_eval, test_eval_data, test_targets, test_targets_frame.
  The editable CI install (`uv pip install -e ".[dev]"`) imports the same way.
- **Expected vs actual:** `make test`'s suite runs in CI with data tests skipped. With pyarrow added to the same venv:
  447 passed, 39 skipped, 25 xfailed (95 s), so pyarrow is the only missing piece.
- **Root cause:** `objdel/data/store.py` imports pyarrow at module level; `objdel.data.__init__` imports `store`;
  `objdel.targets.frame` imports `objdel.data.Item`, and `objdel.eval` tests / `eval_run` import `objdel.targets`.
  pyarrow is only in the `[data]` extra; local `.venv` has every extra, so `make test` never shows it.
- **Fix (pick one):** move `pyarrow` to `[project].dependencies` (the item store is core: targets and eval need it);
  or move `Item` out of `store.py` and import pyarrow lazily inside `write`/`load`; or install `.[dev,data]` in CI
  (keeps the hole for users of the base package).

### B4-3 (medium) training completions cannot reach the eval harness
- **Reproduction:** `test_b4_3_logged_completions_carry_item_key_and_truncation` (real GRPOTrainer, production
  `grpo_config`, `log_completions=True`): `completions_00002.parquet` columns are `step, prompt, completion, format_ok,
  legal, iou_tc, f5, g_area, topo_any, r_topo, advantage_v4, error, advantage`.
- **Expected vs actual:** `eval_run` needs `{"item_id", "samples", "truncated"}` and `hack_audit` needs per-item
  metrics. The table has neither the item key (the prompt is `PROMPTS[dataset]`, identical for every item of a dataset;
  images are stripped) nor the truncation flag (TRL computes it but does not log it). M7's note says "the hack audit
  (M8) can read it".
- **Root cause:** `TRLReward._terms` calls `log_extra("error", …)` only. TRL writes `log_extra` columns into the table
  (`grpo_trainer.py` `log`, `_logs["extra"]`), so the fix is local.
- **Fix:** in `_terms`, also `log_extra("key", keys)` and `log_extra("truncated", [not ids or ids[-1] not in
  self.stop …])`; a 20-line exporter (completions parquet → eval JSONL grouped by key, in order) then closes the arrow.
  Source: TRL GRPO docs, "log_extra" (https://huggingface.co/docs/trl/main/en/grpo_trainer).

### B4-4 (medium, deepens I21) the alarm hook cannot carry the hack audit
- **Reproduction:** `test_b4_4_alarm_receives_hack_audit_statistics`: the dict the Controller passes to `alarm` holds
  `eval_rewards/<term>/{mean,std}`, `eval_error/*`, `eval_defect/*`, loss/KL/length stats; no area ratio, component
  count or duplicate share.
- **Expected vs actual:** I21 reads as "plug hack_audit into `Controller(alarm=...)`". That cannot work as is:
  `hack_audit` takes `eval.metrics` items (per-sample areas vs A_max, size-filtered components, duplicate hashes),
  which need the completions per item (B4-3) and the rater masks; nothing computes them during training.
- **Fix:** an alarm adapter in `objdel.eval` or `objdel.train` that, in `on_evaluate`, reads the last held-out
  completions (with key, B4-3), runs `eval.metrics` on the prepared caches (`metrics(..., caches=cache_loader(key))`)
  and `hack_audit(items, baseline)`, and returns `", ".join(audit["fired"])`. Evaluate it before updating `best`
  (B2-1).

### B4-8 (medium) CI never runs the TRL integration; K13 not enforced
- **Reproduction:** `test_b4_8_ci_runs_the_trl_integration_tests`; shell: the B4-2 venv plus pyarrow, then
  `pytest -q -m "not slow" --cov=objdel` → TOTAL 82 %; `train/callbacks.py` 4 %, `train/configs.py` 0 %,
  `train/dataset.py` 0 %; `test_train_trl.py`, `test_train_configs.py`, `test_train_dataset.py` and 4 tests of
  `test_breakit_reward.py` skip (`importorskip`). Also the end-to-end test of this file skips in CI.
- **Expected:** PLAN K13 "line coverage of src/objdel ≥ 90 %; CI green"; R2 "TRL advantages equal the reference"
  checked continuously. Actual: CI measures neither; no `--cov-fail-under`.
- **Fix:** one Linux job with `.[dev,data,train]` (CPU torch: `uv pip install --torch-backend cpu`,
  https://docs.astral.sh/uv/guides/integration/pytorch/) running `-m "not slow and not data" --cov-fail-under=90`.

### B4-6 (medium) README is pre-Phase-1
- **Reproduction:** `test_b4_6_readme_quick_start_imports_exist` (`objdel.reward.AnnotatorSet` missing).
- **Also stale:** Layout lists `checker/ reward/ models/` only (no `geometry/ targets/ data/ train/ eval/`, no `cpp/`,
  no `_polycore`); "checker … (numpy/scipy only)" (now OpenCV and the C++ extension); Status "Next: Phase 1: finish
  the reward terms, prepare data"; Quick start lacks `make build`, `test-all`, `cov`, `fuzz`, `wheels`; no
  `python -I -m experiments.phase1.<script>` convention (CLAUDE.md has it).
- **Fix:** rewrite Layout from PLAN §3.4, replace the snippet with `score_completion` / `ExpertCache.build`.

### B4-5 (low) Phase 2 gates cannot pass a RIGA-only model
- **Reproduction:** `test_b4_5_gates_accept_a_model_dropped_from_lidc`: a perfect RIGA-only summary → `pass` False
  (every LIDC value None). `test_gates_pass_a_perfect_two_dataset_report` holds.
- **Contradiction:** PLAN §7 Phase 2 "Precision probe … or the model is dropped from LIDC" and "the smallest model
  passing them is smallest viable"; `gates.py`: "Every dataset must pass; a … dataset missing from the report fails".
- **Fix:** `phase2_entry(report, datasets=("riga", "lidc"))`, the caller passing the datasets the probe kept, and the
  report recording the dropped ones.

### B4-7 (low) two scripts have no CLI
- **Reproduction:** `test_b4_7_scripts_without_argparse[census_sizerule|fetch_lidc_pixels]`:
  `.venv/bin/python -I -m experiments.phase1.census_sizerule --help` → `FileNotFoundError: '--help/fixtures/…'`;
  `fetch_lidc_pixels --help` → `TypeError: main() missing … 'out'`. The other 8 scripts print usage
  (now `tests/test_experiment_scripts.py::test_phase1_script_help`, every script; S4 removed the slow
  `test_every_phase1_script_has_help` duplicate).
- **Fix:** argparse with the current positional defaults (as `prepare_rl` / `eval_run` do).

### B4-9 (low) stale references
- PLAN §3.1 boxes `data.ingest` and `data.export: TRL chat rows`, §3.4 `data/ … export`: the code is
  `objdel.data.{riga,lidc}` + `experiments/phase1/build_data.py`, and `objdel.train.dataset.{sft_rows,grpo_rows}`
  (EXECUTION §2 records the move; the diagrams were not updated). `ls src/objdel/data` → `lidc riga splits store`.
- PLAN §4.5 cites `build_splits.py` (a lab file; production is `objdel/data/splits.py`).
- `impl/M7_train.md` (log_completions row): "the hack audit (M8) can read it" — not without B4-3.
- `.github/workflows/tests.yml` comment links the public-repo arm64 announcement; the repo is private, where
  `ubuntu-24.04-arm` works since 2026-01-29 (https://github.blog/changelog/2026-01-29-arm64-standard-runners-are-now-available-in-private-repositories)
  with 2 vCPUs instead of 4. Not a failure; the link should be the private-repo one.

## Diagram arrows with no code

| Diagram | Arrow | Status |
|---|---|---|
| §3.1 | eval hack audit → train.callbacks alarm / rollback | **no code** (I21; B4-4: the hook's input cannot carry the audit) |
| §3.1 | GRPO completions → eval (same checker, metrics, hack audit) | **no code** (B4-3: completions are logged without item key / truncation; no exporter) |
| §3.1 | checkpoint → K sampled completions per val item → `eval_run` JSONL | **no code** (Phase 2 sampler; only `eval_demo` writes synthetic JSONL) |
| §3.1 | "alarm automated in Phase 5" | deferred by design (D9) |
| §3.1 | data.ingest → splits → store → targets / reward.cache → chat rows | code exists (`data.riga/lidc`, `data.splits`, `data.store`, `targets`, `reward.cache`, `train.dataset`); names stale (B4-9) |
| §3.1 | SFT (SFTTrainer) → checkpoint | code exists (`sft_rows`, `sft_config`); exercised only by `smoke_grpo.py` |
| §3.2 | every arrow (parse → validate → rasterize → topology → cache → terms → advantages → aggregator + loggers → typed log) | code exists, run end to end by `test_e2e_store_to_trl_to_controller_to_eval_and_gates` |

## Checks that held

Wiring (end to end, `test_e2e_store_to_trl_to_controller_to_eval_and_gates`, no data, ~12 s):
- synthetic store → `prepare` (train + val) → `grpo_rows` for RIGA and LIDC concatenated → production `TRLReward` +
  `cache_loader` + `stop_ids` in a real `GRPOTrainer` built from `configs.grpo_config("smolvlm-256m")` (dr_grpo,
  sum_then_normalize, mask_truncated_completions, log_completions) with a scripted policy;
- R2 on the production path: every train and eval call's TRL advantages = `advantages_v4(score_batch(...))` of the
  same texts on the prepared caches (≤ 1e-5); RIGA and LIDC groups mixed in one reward call;
- a group of 8 junk completions through TRL: advantages exactly 0; a rater target: legal, IoU_tc ≥ 0.99, topology-any;
  a completion without EOS: `error/truncated`, format 0; `error/incomplete` logged; no NaN in any log;
- `Controller` receives TRL's held-out evaluation, records `best` = `eval_rewards/iou_tc/mean`, calls `alarm` once
  per evaluation;
- the same held-out completions through `eval_run.load_frames` (store frames) + `evaluate` + `build_report`: per-item
  advantages = the ones TRL trained on (≤ 1e-5); report JSON-native; the six Phase 2 gates present; Markdown renders.

Real data (`test_real_lidc_sample_prepared_equals_store_and_eval_frames`, `@pytest.mark.data`, 2 s): 4 LIDC val items
and the first disagreement item: `frame(store item)` W/H/box/scale = the prepared record; `ExpertCache.build` digests
= the prepared caches; `build_target` texts = the prepared targets; every target scores legal, topology-any,
IoU_tc ≥ 0.99 on both cache sets; `eval_run` frames = the prepared caches' rater masks. `data/prepared/manifest.json`
pins the current store manifest SHA-256 and the current reward fingerprint (7,619 / 1,085 LIDC, 591 / 64 RIGA items;
4 LIDC train items skipped without pixels; 311 + 46 disagreement items).

Edges (`test_eval_runner_edges`, scratch `edges.py`): empty JSONL → failing report, no crash; unknown id → KeyError
naming it; an item without pixels is skipped by `prepare` and still scorable by `eval_run`; a dataset missing from a
report fails its gates (by design; B4-5); a single-rater item makes `metrics` raise, but the store has 0 single-rater
items in any split (checked), so it is unreachable.

Single spec: one parser (`geometry.parse`, used by reward and eval), one completion rasterizer (polycore via
`geometry.to_mask`), one topology (`checker.topology`, used by targets, reward cache / score and eval). Other
rasterizing / labeling code is for other inputs: `data/lidc.draw` (PIL, reader contours, pylidc parity P2),
`data/riga.py` scipy FG8/BG4 (stroke extraction; local structures, not the removed checker exports),
`targets/build._rings` (cv2 8-connected components to split parts before polygonizing), `reward_gates.py`
(findContours to *generate* completions).

Contracts: every EXECUTION §3 name exists and its `inspect.signature` matches (polycore, geometry, topology, store,
targets, reward, train.{dataset,rewards,callbacks,configs}, eval, gates); every `__all__` name of each package exists
and is used; no import of `polygon_io`, `annotators`, `AnnotatorSet`, `betti3d`, or checker FG8/BG4 outside
experiments/phase0 and the lab (README's snippet is B4-6).

Install: `uv build` (sdist → wheel) gives `objdel-0.1.0-cp312-abi3`, no `.cpp`, no repo `.pth`; the wheel's tests
pass once pyarrow is present (B4-2). Fresh editable venv (`.[models,train,data,dev]`) on a copy of the tree:
installs, `_objdel_repo.pth` makes `python -I -m experiments.phase1.*` work, a new module / subpackage imports without
reinstall (scikit-build-core's docs say redirect mode needs a reinstall for new files,
https://scikit-build-core.readthedocs.io/en/latest/configuration/editable.html; it works here through the editable
`.pth`'s `src` path), 114 tests of the train / eval / targets / wiring files pass. 8 / 10 phase1 scripts print `--help`.

CI: lint runs `ruff check src tests experiments` like `make lint`, same ruff (0.17.1; `ruff>=0.6` is unpinned, so a
new ruff release can turn CI red); `-m "not slow"` as in `make test`; the slow tests outside
`test_polycore_oracle.py` are data-backed and would skip anyway; actionlint 1.7.12 reports 0 errors; `checkout@v6`,
`upload-artifact@v7`, `cibuildwheel@v4.3.1` exist; cibuildwheel's `pytest` (entry point, not `python -m`) resolves
`from oracles import …` (50 passed).
