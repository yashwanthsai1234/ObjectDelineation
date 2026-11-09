# Spec: Phase CLEAN (archive old code, remove dead code, prefer libraries)

REPO = the ObjectDelineation repository root. All paths are relative to REPO.

## Goal

Only live code stays in `src/`, `experiments/` and `tests/`. Retired code is in `archive/<name>/` with a README.
No hashed result changes. `make test` and `make lint` pass.

## Required behaviours

### Wave 1 (no dependency on other phases)

B1. Bundle manifest. `experiments.phase2.vast.policy.code_files() -> list[str]` never returns a path that starts with
    `.claude/`, `archive/`, `data/` or `docs/`, or that contains `/lab/` or `/results/` directly under
    `experiments/<phase>/`, or that matches `riga*.csv|json|jsonl` (case-insensitive). It still returns untracked,
    not-ignored files outside those roots. It still returns `tests/oracles/reward_reference.py`,
    `experiments/phase1/reward_gates.py` and `src/objdel/__init__.py`, and still stops with an error when no
    `tests/oracles/` file is present. `.gitignore` ignores `/.claude/worktrees/`.

B1a. Bundle imports. A copy of exactly the files `code_files()` returns, plus the built `objdel._gpucheck`
    extension, imports `objdel`, `objdel._gpucheck`, `experiments.phase1.reward_gates`,
    `tests.oracles.reward_harness`, `experiments.phase2.{train_rl,g_probe,select_ckpt,final_eval}` under `python -I`,
    and every first-party module loaded comes from the copy. A test checks this with no GPU.

B2. `python -m experiments.phase2.train_rl` has no `--procs` option: passing `--procs` is an argparse error (exit
    code 2). `--help` output does not contain `--procs`. `train_rl` creates no process pool. `train_rl.main` calls
    `gc.freeze()` exactly once, after the last `serve.preload`. All other options and the `rl.json` output keys are
    unchanged.

B3. These paths hold no tracked file: `experiments/phase0/`, `experiments/phase1/lab/`, `experiments/phase2/lab/`,
    `experiments/phase1/census_sizerule.py`, `experiments/phase1/eval_demo.py`, `experiments/phase1/sample_eval.py`,
    `experiments/phase1/smoke_grpo.py`. These stay importable as modules:
    `experiments.phase1.{build_data,build_targets,fetch_lidc_pixels,prepare_rl,eval_run,reward_gates}` and
    `experiments.phase1.bench.xarch_hashes`.

B4. `archive/` subfolders exist: `phase0/`, `phase1_lab/`, `phase1_scripts/`, `phase2_lab/`, `smolvlm/`, and
    `controller/` (unless the owner keeps the Controller, see B8). Each direct subfolder of `archive/` has a
    `README.md` containing the headings `## What it was`, `## What replaced it`, `## Retiring commit`,
    `## How to restore`. The restore text names the tag `pre-clean-2026-11-08`. `archive/` has no `__init__.py` at
    its root.

B5. No file under `tests/` imports a module from `archive` or a module path that contains `.lab.`, inserts a lab or
    archive folder into `sys.path`, or contains the string `experiments/phase1/lab` or `experiments/phase2/lab` in
    code (docstrings may cite `archive/...`). No file under `src/objdel` imports `experiments` or `tests`.
    No live file (outside `docs/changelog.md`, `docs/research/**`, `experiments/*/impl/**` and the frozen phase
    `.md` records) names `experiments/phase0`, `experiments/phase1/lab` or `experiments/phase2/lab`. Every
    backticked repo-relative path in `README.md` and `docs/**/*.md` (same exemptions) exists.

B6. Test fixtures that came from the phase 2 lab live in `tests/fixtures/`: `vast_snapshot.json` and
    `topo_probe_v3e_seed*.json`. The Dockerfile digest check in `tests/test_vast.py` uses a constant in the test,
    not a lab file.

B7. Golden for the target builder. `tests/goldens/targets_lab.npz` exists and holds no masks: every array has
    `ndim <= 1` and a non-bool dtype; RIGA appears only as aggregate values (one lab vertex total, one case count, one
    hash), never per case. `tests/oracles/crops_oracle.py` provides the `pad_div` crop. A test named
    `test_fixtures_differential_vs_lab_golden` in `tests/test_targets.py` (marker `data`) reads the masks from
    `data/fixtures`, crops them with the oracle, checks the recorded crop hashes, runs the production builder, and
    asserts: per LIDC case the production nodule vertex count equals the recorded lab count; the production LIDC text
    hash equals the recorded lab text hash in at least 99 % of LIDC cases; the RIGA production vertex total is at most
    the recorded lab total. The old test `test_fixtures_differential_vs_lab_builder` does not exist. A test fails
    when any `tests/goldens/*.npz` other than `polycore.npz` holds an array with `ndim > 1` or bool dtype.

B7a. Golden for the phase 2 lab references. `tests/goldens/p2_lab_refs.json` holds the lab outputs on seeded
    synthetic inputs. These tests still exist, pass, and compare live code to it with the old tolerances:
    `tests/test_p2_signal.py::test_live_b2_functions_and_group_scoring` (`objdel.eval.signal.per_item`, `group`),
    `tests/test_p2_signal.py::test_c2_stratified_bootstrap_equals_power_c2` (`c2_analyze.pooled_ci`, `WEIGHTS`),
    `tests/test_p2_probes.py::test_stratified_bootstrap_equals_power_c2` (`objdel.eval.probes.stratified_ci`).
    None imports a lab module.

B8. Controller retired (default; skip if the owner decides to keep it):
    - `objdel.train.callbacks` exports `TraceStep` and has no `Controller`, `audit_alarm`, `FIRED`.
    - `objdel.eval` exports `hack_audit`, `gates`, `metrics`, `summarize`, `vision`; it has no `alarms_for`.
    - `objdel.train.rewards.TRLReward` instances have no `recorded` attribute; reward values are unchanged.
    - Every function in `TRLReward.funcs()[0]` still has `trl_reward` set to its `TRLReward` (FastGRPO prefetch reads it).
    - The end-to-end wiring test in `tests/test_breakit_wiring.py` (store -> reward -> TRL -> eval -> gates) still
      exists and runs without the Controller.

B9. SmolVLM retired:
    - `objdel.train.configs.MODELS` and `MAX_COMPLETION` have no key `"smolvlm-256m"`.
    - `objdel.train.configs.processor_kwargs(model, dataset)` returns `{"padding_side": "left"}` for every model key.
    - `objdel.train.configs.grpo_config(model, output_dir, **overrides)` sets `use_vllm=True` unless
      `overrides["use_vllm"]` says otherwise; `grpo_config("qwen3-vl-2b", out, use_cpu=True, use_vllm=False)` gives a
      config with `use_vllm == False`.
    - No file in `src/`, `experiments/phase2/*.py` or `tests/` contains `smolvlm` (case-insensitive).

B10. One box-overlap rule. `objdel.geometry.box.boxes_meet(box, others) -> bool`: `box = (x0, y0, x1, y1)` half-open;
     `others` an iterable of `(a, b, c, d)`; returns True when some `(a, b, c, d)` has `a < x1 and x0 < c and b < y1
     and y0 < d`; False for an empty iterable. `objdel.data.lidc.crowded(s, box) == boxes_meet(box, s.others)`.
     `objdel.targets.frame.crowded(item, box) == boxes_meet(box, item.meta.get("other_nodules", ()))`.

B11. `make preflight` with no `CONFIG` uses `experiments/phase2/configs/r1_qwen3vl2b.yaml` (a config that
     `preflight` can load; `rl_shifted.yaml` is a steps-only run file).

B12. `experiments.phase2.train_sft` has a function `load_items` and no function `preload`.

B13. Dead-code gate. `pyproject.toml` dev extra pins `vulture==2.16`. `pyproject.toml` has a `[tool.vulture]` table:
     `paths = ["src", "experiments", "tests"]`, `min_confidence = 60`, and `ignore_names` / `ignore_decorators` for
     accepted false positives, each group with a comment. `make dead` runs vulture with that table, exits 0 and prints
     nothing. No whitelist module exists. Parameters of test fakes keep their real names.

B14. `tests/test_repo_hygiene.py` exists, needs no data and no GPU, and contains tests that fail when: `make dead`'s
     vulture run reports anything, or vulture is not installed (fail, not skip); B5 is broken (imports, lab strings,
     stale paths, missing doc paths); `smolvlm` appears in `src/` (and, after Wave 3, `lfm`). The `.claude/` and
     `archive/` bundle check is in `tests/test_vast.py` only (B1).

B15. `pyproject.toml` `[tool.ruff] extend-exclude` is `["archive"]`.

B16. `objdel.checker.volume.topology_any` exists and behaves as before; its module docstring does not claim the reward
     uses it.

B17. Local data (git-ignored): `data/views/lidc-bundle`, `data/views/lidc-bundle-v2` and
     `data/eval_views_archive/probe-precision-r1` are moved under `data/_retired/2026-11-08/` with a `LIST.txt`
     (path, size, `rows.jsonl` SHA-256). `data/views/lidc-bundle-v3` and `data/eval_views/lidc-val-rotated` stay.

B18. A new ADR `docs/architecture/decisions/NNNN-library-keep-list.md` (next free number) has the ADR template
     headings and lists KEEP, DEFERRED (with owning phase) and NOT ADOPTED libraries.

B19. Optional (only if the read-only offer query with vastai 1.11.0 returns every offer field `policy.py` reads):
     `vastai` is not a dependency in `pyproject.toml` (its exact pins on Pillow and cryptography would change the
     project environment). `docs/runbooks/gpu-run.md` names the pinned version and the isolated install. `preflight`
     has a check that fails when `vastai --version` differs from the pinned version. `policy.vast()` stays a
     subprocess call to `vastai`. `pyyaml` stays in `dev`.

B20. Optional: `tests/test_checker_shapely_crosscheck.py` (shapely 2.2 in dev extra) checks the CUDA checker's
     defect class against `shapely.is_valid_reason` for the mapped classes, and the raster mask against
     `shapely.contains_xy` at pixel centres on legal polygons. `tests/oracles/validity_oracle.py` and
     `raster_oracle.py` remain.

### Wave 2 (after phase D puts the CUDA polygonizer in production)

B21. `src/objdel/checker/polygonize.py` does not exist; `objdel.checker` does not export `polygonize`. The numpy
     polygonizer is available to tests as `tests/oracles/polygonize_oracle.py` with the same function
     `polygonize(mask, offset=(0, 0))`. `checker/certify.py` and `checker/bitquad.py` are not in `src/objdel` once
     phase D's device check covers their callers; `objdel.checker` then exports neither `bitquads` nor `certify`, and
     no file in `src/`, `tests/` or `experiments/` (including `experiments/phase1/bench/xarch_hashes.py` and
     `experiments/cuda_checker/`) imports them from `objdel.checker`.

### Wave 3 (after phases D5, E2, E6)

B22. `objdel.train.configs.MODELS` (or its phase E home) has only Qwen3-VL entries (one today: `"qwen3-vl-2b"`); no
     `lfm` string in `src/`, `experiments/phase2/*.py`, `tests/`. No full-fine-tune code path (`param_groups`, a
     float32 `dtype` branch). The slow CPU smoke test uses a tiny random Qwen3-VL with a completion cap.

B23. The CPU process pool (`fork_pool`, `pool_size`, `--procs` in `select_ckpt`, `g_probe`, `final_eval`) is gone
     after evaluation metrics run on the GPU.

## Interfaces

| Name | Signature | Notes |
|---|---|---|
| `experiments.phase2.vast.policy.code_files` | `() -> list[str]` | B1 |
| `experiments.phase2.train_rl.main` | `(argv: list[str] \| None = None)` | B2; no `--procs` |
| `objdel.geometry.box.boxes_meet` | `(box: tuple[int, int, int, int], others: Iterable[tuple[int, int, int, int]]) -> bool` | B10 |
| `objdel.data.lidc.crowded` | `(s: Slice, box) -> bool` | unchanged result |
| `objdel.targets.frame.crowded` | `(item: Item, box) -> bool` | unchanged result |
| `objdel.train.configs.grpo_config` | `(model: str, output_dir: str, **overrides) -> GRPOConfig` | B9 |
| `objdel.train.configs.processor_kwargs` | `(model: str, dataset: str) -> dict` | B9 |
| `experiments.phase2.train_sft.load_items` | same parameters as the old `preload` | B12 |
| `make dead` | Makefile target; reads `[tool.vulture]` | B13 |
| `tests/goldens/targets_lab.npz` | npz of 1-D arrays; LIDC keys per case `"lidc/{item}/{rater}"`; RIGA aggregate keys only | B7 |
| `tests/goldens/p2_lab_refs.json` | JSON of lab outputs on seeded synthetic inputs | B7a |
| `tests/oracles/crops_oracle.py` | `pad_div` crop of the phase 1 lab, ported | B7 |

## Constraints

C1. No change to: `objdel.reward.spec.fingerprint()` output; the target-set hashes (`target_set_sha256` in
    `data/targets` manifests and the value checked by `tests/test_p2_views.py`); any file under `cuda/`.
C2. Nothing is deleted without a copy: the owner tags the pre-clean commit `pre-clean-2026-11-08` and pushes the
    tag; byte-identical duplicates and the SmolVLM branches are deleted, all else is moved.
C3. Agents run no git command that changes anything; the owner runs tag, add and commit.
C4. Python that reads `data/` runs with `-I`.
C5. `make test` and `make lint` pass after each task. Every removed test is listed by node id in the README of the
    archive folder that holds it.
C6. No GPU time and no money.
C7. No Claude or AI attribution in commits or pull requests.
C8. No tracked file holds RIGA masks or per-case RIGA results.
C9. No CLEAN commit contains the owner's unrelated uncommitted work (`README.md`, `docs/`).
C10. The project environment's pinned packages (Pillow included) do not change.
C11. Order: Wave 1a (B2, B8, B9) lands before any phase B2, B3 or C branch starts.

## Acceptance gates

| # | Gate | Target |
|---|---|---|
| G1 | `make lint` | exit 0 |
| G2 | `make test` | 0 failed; passed = 1447 minus listed removed tests plus new tests (baseline 1447 passed, 6 skipped, 31 deselected) |
| G3 | `fingerprint()` | equal before and after |
| G4 | target-set hashes | equal before and after |
| G5 | tracked files under `experiments/phase0`, `experiments/phase1/lab`, `experiments/phase2/lab` | 0 (from 426) |
| G6 | `code_files()` under `.claude/` or `archive/` | 0 (from 2 and 10) |
| G7 | tests importing, path-inserting or reading lab or archive code | 0 (from 3 importing and 2 reading files) |
| G8 | `make dead` | exit 0, 0 output lines (from 56 at the gate scope) |
| G9 | `smolvlm` in src, experiments/phase2/*.py, tests | 0 (Wave 1); `lfm` 0 (Wave 3) |
| G10 | `Controller`, `alarms_for` in src and experiments/phase2 | 0 (default branch) |
| G11 | archive READMEs with the four headings | every direct subfolder |
| G12 | `make preflight CONFIG=experiments/phase2/configs/r1_qwen3vl2b.yaml` | each check's pass value equal to before CLEAN |
| G13 | Wave 2: `src/objdel/checker/polygonize.py` | absent, with G4 holding |
| G14 | bundle import test (B1a) | passes |
| G15 | `tests/goldens/*.npz` except `polycore.npz` | no array with `ndim > 1` or bool dtype |
| G16 | B7a tests and `trl_reward` assertion (B8) | pass |
