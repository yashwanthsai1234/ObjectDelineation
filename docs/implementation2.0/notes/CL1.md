# CL1: dead and legacy code audit (CLEAN phase)

Repo: ObjectDelineation, branch gpu-scorer-open-once, 802 tracked files. No REPO file was edited. No git change was made.
Tool: vulture 2.16 in a temp venv (`docs/implementation2.0/notes/venv`), not in REPO/.venv. Raw output: `docs/implementation2.0/notes/v_src_*.txt`.

## Skills used
- nvidia-skills:nvidia-skill-finder: loaded; no catalog lookup needed (cleanup task, no NVIDIA product choice). No change.
- engineering:code-review: used its correctness/maintainability lens on the callers of each candidate before marking it dead.
- engineering:debug: used reproduce-isolate method for "who still imports this" (grep + vulture, then pytest collection).
- engineering:system-design, engineering:architecture: ADR shape used for the archive mechanism decision (options table, trade-offs).
- engineering:documentation: shaped the archive README requirement (what was it, what replaced it, how to run again).
- huggingface-skills:huggingface-best: loaded; no model search needed (no model change in this task).
- engineering:tech-debt: Impact/Risk/Effort scoring used for the ranking in Recommendations.
- mattpocock-skills:codebase-design: "deletion test" and "one adapter = hypothetical seam" used to judge pass-through modules and the numpy polygonizer (oracle, not a seam).
- Web tools loaded (ToolSearch). No web research was needed: this is a repo-internal static analysis.

## Method (evidence)
1. `vulture src` alone: 90 hits; nearly all are functions called only from experiments/tests (false positives for a library).
2. `vulture src + non-lab experiments`: 9 real candidates in src (list below). Adding tests removes only 6 (TrainerCallback hooks, false positives).
3. grep import graph from tests, src, vast policy, Makefile, CI to every phase0/phase1/lab path.
4. md5 over all tracked .py/.cpp/.hpp/.h for exact duplicates.
5. `make test` (not slow) started as baseline; see "Baseline" at the end.

## Findings

### F1. Size of the legacy layer [measured: git ls-files, du]
- experiments/phase1/lab: 276 tracked files, ~6.4 MB (t_data 34, t_dynamics 18, t_polycore 26, t_reward 95, t_stack 70, t_targets 24, common 2, tools 1, plus 6 md).
- experiments/phase2/lab: 130 tracked files (A1,A3,A4,A5,A6,B1..B5), ~1.5 MB.
- experiments/phase0: 3.0 MB (7 scripts + results csv). experiments/phase1 total 7.6 MB. phase2 3.2 MB.
- phase1/lab is already ruff-excluded (`pyproject.toml: [tool.ruff] extend-exclude`) and is dropped from the Vast bundle (`experiments/phase2/vast/policy.py:379` regex `^experiments/[^/]+/(lab|results)/`). It is frozen evidence, not code.

### F2. Exact duplicate files [measured: md5 of tracked files]
28 files exist twice inside phase1/lab, byte-identical:
- `t_stack/vendor/polycore_src/**` (16 files) == `t_polycore/**` (cpp/hpp, _polycore.cpp, python/polycore/__init__.py, tests/*, bench/bench_k3k4.py).
- `t_stack/vendor/targets/*` (5) == `t_targets/targets/*` (simplify, fastgeom, crops, builder, __init__).
- `t_stack/vendor/rewardcore/*` (6) == `t_reward/rewardcore/*`.
- `t_stack/vendor/oracle.py` == `common/oracle.py`.
`vendor/` is a copy for running the stack lab standalone. Safe to delete (originals exist).
Also: `archive/cpu_checker/cpp/polycore/polycore.hpp` and `lab/t_polycore/cpp/polycore/polycore.hpp` differ (not duplicates; the archive one is the promoted, later version).

### F3. Dead or tests-only code inside src/objdel [code: vulture; grep]
Only 9 candidates survive when experiments (non-lab) count as callers:
| Symbol | Where | Who uses it | Verdict basis |
|---|---|---|---|
| `topology_any` | checker/volume.py:62 | tests/test_volume.py only; docstring says reward compares with it, reward does not | dead for production |
| `nsd`, `ged` (public wrappers) | eval/metrics.py:215,244 | nobody; `_nsd`, `_ged`, `_nsd_best` are used by `metrics()` | dead wrappers (check test_eval first) |
| `INVERSE` | targets/augment.py:38 | nobody (`DATASET_OPS` is used by build_views) | dead constant; may be wanted for a rotated-view check |
| `ExpertCache.same_sig`, `.atexit` | reward/cache.py:235,256 | vulture noise (attributes set for weakref/atexit) | keep |
| `FastGRPO._iter_named_params`, `.floating_point_ops` | train/fast.py:503,692 | TRL/Trainer hooks | keep |
| `checker.topology` `writeable` | topology.py:63 | numpy flag assignment | keep |
Larger tests-only cluster (used by tests and doc strings, not by train_rl / g_probe / select_ckpt / final_eval):
- `train/callbacks.py: Controller, audit_alarm` and its dependency `eval/audit.py` (`alarms_for`), plus `rewards.py` "recorded" audit hook: the production RL (`experiments/phase2/train_rl.py`, `g_probe.make_trainer`) uses FastGRPO and never installs Controller. Only `tests/test_train_trl.py` and the Phase 1 `smoke_grpo.py` doc mention it. [code: grep Controller src experiments/phase2 tests]
- `train/callbacks.py: TraceStep` is used (g_probe, train_sft): keep.
- `train/configs.py`: `MODELS` has `smolvlm-256m` and `lfm2.5-vl-450m`; `MAX_COMPLETION` and `sampling` carry LFM/SmolVLM branches; `use_vllm: model != "smolvlm-256m"` (configs.py:71-75,110,136,192). Production model is qwen3-vl-2b only.
- `geometry/text.py:19` and `train/sample.py:22` docstrings quote LFM/SmolVLM token costs (comments only).
- `eval/probes.py` (616 lines), `eval/signal.py`, `eval/vision.py`, `eval/coverage.py`, `eval/gates.py`: used by select_ckpt/g_probe/prereg (alive), but tied to the R1 two-model selection protocol (G-vision, G-signal). Keep until Phase B/C decides the exam; not dead.

### F4. Worker pools [code: grep ProcessPool|Pool|fork_pool]
- `train/serve.py: fork_pool/_Pool/pool_size` : a CPU ProcessPoolExecutor used only for exam item metrics (`select_ckpt._score` via `pool.map`, `g_probe.score`, `final_eval`). `train_rl.py:50` creates `pool = serve.fork_pool(a.procs)` and only calls `pool.shutdown()` at line 61: it is never submitted a job in RL training. That pool and `--procs` in train_rl are unused (idle forked processes holding STORE). Remove from train_rl.
- Scoring itself is on the GPU (`reward.score_batch` -> `_gpucheck.score_polygons`), so the pool only covers `objdel.eval.metrics` (numpy/scipy per-item metrics: nsd, hm_iou, ged, soft_dice). Owner rule "GPU everywhere": the pool is the CPU leftover; it disappears when metrics move to the GPU (phase D).
- `train/dataset.py: prepare` uses `multiprocessing.Pool(3)` (offline data prep); `experiments/phase2/build_views.py:101` and `build_eval_views.py:73` use `Pool(procs)` (offline view building). Offline, CPU, numba/numpy target builder. Alive; replace when the polygonizer and mask ops move to CUDA.
- `reward/score.py: score_completion, score_mask, score_many` are the CPU-shaped API over a GPU call; `score_mask` is used by `eval/metrics.py:263` and `experiments/phase2/baselines.py`. Keep, they call the device.

### F5. The numpy polygonizer [code: grep polygonize]
- `src/objdel/checker/polygonize.py` (152 lines; boundary_edges, link, ring_labels, list_rank, polygonize) is NOT oracle-only yet. Production callers:
  1. `targets/build.py:53,93` (SFT target builder: `polygonize(sub, offset=(x0, y0))`).
  2. `checker/certify.py:18,46` (used by `data/riga.py:226` when extracting RIGA masks).
  3. `checker/__init__.py` re-export.
- Test and experiment callers: `tests/test_checker.py`, `tests/test_targets.py:14,28,32,314`, `tests/test_cuda_checker.py:131,158` (CUDA vs polygonize fuzz), `tests/test_cuda_checker_gpu.py:46`, `experiments/cuda_checker/compare.py:21,167`, `eval.py:42,247`.
- The CUDA port exists: `cuda/gpucheck.cpp:17` "port of objdel.checker.polygonize + bitquad", `_gpucheck.check_masks(bits, hw, start)` returns rings in polygonize's order (`gpucheck.h:106`), pinned equal by tests/test_cuda_checker.py.
- Therefore: it can become an oracle only after build.py and certify.py call `_gpucheck.check_masks` (or a thin `objdel.geometry.rings` wrapper). Also `bitquad.py`, `topology.py` (scipy/OpenCV Spaghetti), `certify.py` are the other two methods of the same certificate; `topology.py` is also used by reward cache and targets, so it stays in src.

### F6. Phase 0 scripts (experiments/phase0) [code: grep]
| Script | Referenced by | Note |
|---|---|---|
| baseline_medsam.py | `targets/frame.py:49` comment (LUNG_WINDOW) | MedSAM baseline, finished; results in phase0/results |
| baseline_vlm.py | nothing | zero-shot VLM baseline, finished |
| lidc_census.py, riga_census.py | `phase1/census_sizerule.py`, README | census, finished; numbers cited in README/RESULTS |
| lidc_fetch.py | nothing | LIDC download helper, superseded by `phase1/fetch_lidc_pixels.py` |
| riga_outlines.py | nothing | RIGA outline stats, finished |
| summarize_0b.py | nothing | summary writer for 0b |
| results/*.csv/json (3.0 MB) | README links, test_p2_breakit_ops.py:299 (string check only) | published evidence |
No test imports a phase0 module. `phase0/__init__.py` exists only so `-m experiments.phase0.*` works.

### F7. Phase 1 scripts (experiments/phase1, non-lab) [code: grep, tests]
| Script | Live users | Verdict basis |
|---|---|---|
| build_data.py | tests/test_data_build.py, test_vast.py (R0 module); README | data store builder; alive |
| build_targets.py | tests/test_eval.py; README | target set builder; alive |
| fetch_lidc_pixels.py | `data/lidc.py` docs, tests (test_p2_views, test_data_lidc, test_breakit_wiring) | alive |
| prepare_rl.py | `train/dataset.py`, `configs.py`, README, tests | alive (prepared caches) |
| eval_run.py | `eval/gates.py`, `train/sample.py`, build_eval_views, select_ckpt, 6 test files | alive, core exam harness |
| reward_gates.py | `phase2/k3_harness.py`, `vast/policy.py`, tests/oracles/reward_harness.py, test_p2_serve | alive (imported by phase 2) |
| bench/xarch_hashes.py | `configs/r0.yaml`, policy.py, test_vast.py | cross-arch hash gate for R0; alive for preflight |
| census_sizerule.py | only tests/test_experiment_scripts.py, test_breakit_wiring.py (--help runs) | legacy census; no pipeline user |
| eval_demo.py | only tests/test_eval.py | demo of eval on rater targets; legacy |
| sample_eval.py | only tests/test_experiment_scripts.py (bf16 check) | HF-generate sampling on a Mac; superseded by phase2 sample path |
| smoke_grpo.py | only tests/test_experiment_scripts.py (group_stats) | Mac SmolVLM-256M K12 smoke; SmolVLM-only path |
phase1 `.md` (PLAN, EXECUTION, RESULTS, impl/**: 38 md files) are history; README links PLAN/RESULTS.

### F8. LFM / SmolVLM / small-model only paths [code: grep -i lfm|smolvlm|liquid]
- src: `train/configs.py` (MODELS, MAX_COMPLETION, `sampling` LFM branch, `lora` handling, smolvlm size/use_vllm), docstrings in `geometry/text.py`, `train/sample.py`.
- experiments/phase2: `configs/r1_lfm.yaml`, `r0.yaml` (parity arm), `r0_qwen3vl_arm.yaml`; `preflight.py:35 HF_RATE` map; `project_r1.py:33 R1` table (LFM stop hours); `Makefile:28 CONFIG ?= experiments/phase2/configs/r1_lfm.yaml` (default points at a retired model); `hf_push.py` docstring + LFM license text (§4(a)); `train_sft.py:340` comment on LFM key conversion; hf_rate/hf_sample/r0_parity/g_probe usage strings.
- `vast/ledger.csv` mentions r1_lfm_resume.yaml and _r1_lfm_select.yaml (history, keep as ledger).
- Memory/owner decision (2026-11-02): next run Qwen3-VL-2B, no LFM retrain. Brief section 3: open to any model up to 4B. So the SmolVLM arm (HF-generate fallback, longest_edge sizing) and the LFM sampling/limits are legacy; they also cost tests (`tests/test_train_configs.py`, `test_train_trl.py` param cases).
- Phase 1 lab: `t_stack/t3_grpo.py, t4_sft.py, results/t3_*` for SmolVLM-256M, LFM2.5-VL-450M, Qwen3.5-0.8B (frozen evidence).

### F9. Experiment scripts under phase2 [code: grep]
Alive (library-like, tested): g_probe, train_rl, train_sft, select_ckpt, final_eval, build_views, build_eval_views, preflight, parity, probes, profile, baselines, hf_pull/push/rate/sample, sample_vllm, prereg, k3_harness.
One-offs for finished decisions (no config/Makefile/docs users): `c2_analyze.py` (tests=5, docs=0), `project_r1.py` (R1 hour projection; tests=2, docs=0), `r0_parity.py` (R0 arm), `hf_rate.py`/`hf_sample.py` (R0 HF parity), `k3_harness.py` (K3 throughput gate; docs=0 but vast policy runs it, keep).
Configs (16): finished run plans `c2_a0_seed1, c2_a0p_seed2, c2_a0p_seed3, p1_rl, r0, r0_qwen3vl_arm, r1_lfm, r1_qwen3vl2b, r2_qwen3vl2b, r3_gpu_test, r4_speed, rx8, checker_gpu`; live: `r5_shared`, `rl_full`, `rl_shifted` (plan). Preflight tests read some of them.
phase2 docs: PLAN.md, EXECUTION.md, JUDGE_R1..R4, RESEARCH.md, V0_V1.md, RESULTS.md, impl/ (13), perf/ (5): history; owner said ignore phase-2 PLAN.md.

### F10. Phase 2 lab [code: tests grep]
Used by tests (cannot move without a change): `lab/B2/{analyze,posthoc,score}.py` (tests/test_p2_signal.py:95-97), `lab/B3/power_c2.py` (test_p2_signal.py:190, test_p2_probes.py:175 importorskip), `lab/B3` path (test_p2_probes.py:17), `lab/A3/vast_snapshot.json` (test_vast.py:31), `lab/B5/Dockerfile.phase2` (test_vast.py:789, digest pin). Everything else in phase2/lab (A1, A4, A5, A6, B1, B4, most of B2/B3) is frozen evidence; `src/objdel/eval/{signal,vision}.py` docstrings say they were promoted from B1/B2.

### F11. Tests that read legacy paths [code: grep]
- `tests/test_targets.py:383-414` `test_fixtures_differential_vs_lab_builder` imports `experiments/phase1/lab/t_targets/targets/{builder,crops}` via `sys.path.insert` (`data` marker). The only test that executes lab code.
- `tests/oracles/*` are docs-referenced ports of lab oracles; they do not import lab. Already the right "oracle" home (raster_oracle, validity_oracle, topology_oracle, reward_reference, reward_variants, score_numpy, reward_harness).
- `tests/test_experiment_scripts.py` runs `--help` on phase1 scripts.
- `tests/test_vast.py` hard-codes `experiments.phase1.bench.xarch_hashes`, `experiments.phase1.reward_gates`, `experiments.phase1.build_data` (policy chain). Moving these needs policy.py and test updates together.
- `experiments/phase2/vast/policy.py:code_files()` uses `git ls-files -co`: includes UNTRACKED files. `code_files.txt` today lists `.claude/worktrees/agent-.../` and `.claude/settings.json` (302 lines): a bundle would rsync agent worktrees. Also an `archive/` folder is NOT dropped by the regex (`^(data|docs)/|^experiments/[^/]+/(lab|results)/`), so it would ship. [code: policy.py:377-380]

### F12. Duplicate helpers [code: grep, vulture]
- Thin wrappers, not duplicates: `select_ckpt.view_rows` (wraps `dataset.prompts`), `select_ckpt.stop_ids` (wraps `rewards.stop_ids` for a path). Keep.
- Real near-duplicates: `crowded` in `data/lidc.py:257` and `targets/frame.py:124` and `phase1/build_data.py:66`; `lidc_series` in phase1/build_data.py:83 and phase0/baseline_medsam.py:131; `per_image`, `write_csv` in phase0 riga_census and riga_outlines; `window` mask helper in phase0/baseline_medsam; `_ring` x3 (shape.py:95, text.py:98, plus tests); `records` in `train/dataset.py:159` and `eval/probes.py:467`. Low value: the phase0 ones disappear with the archive; fix `crowded` x2 in src only if kept.
- `train_sft.preload` vs `train/serve.preload`: different jobs, same name (confusing); rename `train_sft.preload` -> `load_items`.

### F13. Static quality [measured: ruff, vulture]
- `ruff --select F401,F841,F811` on the two lab trees: 4 findings (2 F401, 2 F841): negligible.
- tests: 100%-confidence unused variables (`rid` in test_vast/test_hf_pull/test_p2_hf_push, `has_tool_images`/`multimodal_fields` in test_breakit_wiring and test_train_trl, `return_tensors` test_p2_v1_select, `hlo/hup` in phase1/reward_gates.py:138): trivially removable.

## Recommendations (numbered; Priority = (Impact+Risk) x (6-Effort))

Archive mechanism (decision, ADR-style)
- Context: repo already uses `archive/cpu_checker/` with a README table ("Was / replaced by / run again"). `.venv`-level imports must stay green; Vast bundle ships tracked files; CI lints `src tests experiments`.
- Option A: `archive/<name>/` folder + README (existing precedent). Pros: greppable, code visible, one clone. Cons: lint/bundle/test discovery must be told to skip it; tree gets bigger.
- Option B: git tag only (e.g. `archive/phase1-lab-2026-11-08`), delete from tree. Pros: tree is clean, no CI rules. Cons: hard to browse, evidence links in README/docs break, tests that read lab need replacement first.
- Option C (recommended, hybrid): FIRST tag `archive/pre-clean-2026-11-08` on the current commit (owner runs `git tag`; no agent runs it). THEN `git mv` to `archive/phase1_lab/`, `archive/phase0/`, etc., each with a README using the cpu_checker template. DELETE pure duplicates (F2) and dead code (F3) with no archive copy (tag is the copy).
- Rules that keep tests green: (a) `archive/` has no `__init__.py`; (b) pyproject `[tool.pytest.ini_options] norecursedirs = ["archive"]` and `testpaths = ["tests"]` already scopes tests; (c) `[tool.ruff] extend-exclude = ["archive", ...]` replacing the phase1/lab entry; (d) `.github/workflows` unchanged (ruff paths `src tests experiments` do not include archive); (e) `policy.py:379` regex add `|^archive/`, and drop the `-co` untracked inclusion (or add `.claude/`); (f) docs links: README/docs table update (see 12).

1. Fix the Vast bundle manifest first (smallest, highest risk). `experiments/phase2/vast/policy.py:code_files()`: use tracked files only (`git ls-files -z`, no `-co`), add `archive/` and `.claude/` to `drop`; update `tests/test_vast.py` (lines ~528-580, 1189-1198) and regenerate `vast/code_files.txt` (it is a generated state file: add to .gitignore). Effort 1.
2. Remove the idle pool from `train_rl.py` (lines 35, 50, 61: `--procs`, `serve.fork_pool`, `pool.shutdown`); update `tests/test_p2_train_rl.py`. Effort 1. Keep `fork_pool` for select_ckpt/g_probe/final_eval until their metrics are on the GPU (phase D), then delete `train/serve.py: _Pool, fork_pool, pool_size, cpus, quota` and the `SCORERS` import dance.
3. Delete exact duplicates: `experiments/phase1/lab/t_stack/vendor/` (28 files). Zero test impact (grep: no test reads vendor). Effort 1.
4. Archive Phase 1 lab (216 remaining files, 6 MB) to `archive/phase1_lab/` after step 6 (the differential test). Delete from the tree: `t_polycore/` (the promoted polycore is already in `archive/cpu_checker`; its tests are duplicated there), `t_dynamics/figs/*.png` unless README links them. Keep README table in `archive/phase1_lab/README.md`. Remove ruff extend-exclude entry.
5. Archive Phase 0 (`experiments/phase0/*`, 3 MB) to `archive/phase0/` except the results the README cites (`results/lidc/clusters.csv`, `summary.json`, RESULTS.md, RESULTS_0b.md): move those evidence files to `docs/research/evidence/phase0/`. `phase1/census_sizerule.py` reads phase0 census modules: archive census_sizerule with them (item 7). Update README lines 30, 98-99, 113, 116, 184 links.
6. Retire the only lab-executing test: `tests/test_targets.py::test_fixtures_differential_vs_lab_builder`. Replace with a golden: record (once, before the move) the lab builder's per-item texts hash / vertex counts on the fixture crops into `tests/goldens/targets_lab.json` (same pattern as `archive/cpu_checker/record_goldens.py`), and assert prod_v <= recorded lab_v and same >= 0.99. Do this BEFORE item 4.
7. Archive phase1 one-offs: `census_sizerule.py`, `eval_demo.py`, `sample_eval.py`, `smoke_grpo.py` -> `archive/phase1_scripts/`. Keep `build_data, build_targets, fetch_lidc_pixels, prepare_rl, eval_run, reward_gates, bench/xarch_hashes` (live imports) and do not move them yet: they are referenced by `vast/policy.py`, `configs/r0.yaml`, 6+ tests. Later (phase E layering) promote them to `src/objdel/...` or `experiments/data/` in one commit with policy.py/test_vast updates. Update `tests/test_experiment_scripts.py` (drop census_sizerule/sample_eval/smoke_grpo cases; keep `--help` for live ones). Effort 2.
8. Polygonizer to oracle (needs code first): (a) add a wrapper `objdel.geometry.rings(masks)` over `_gpucheck.check_masks` returning `Ring` objects; (b) switch `targets/build.py:93` and `checker/certify.py:46` to it; (c) `git mv src/objdel/checker/polygonize.py tests/oracles/polygonize_oracle.py`; (d) update imports in tests/test_checker.py, test_targets.py, test_cuda_checker.py, test_cuda_checker_gpu.py, experiments/cuda_checker/{compare,eval}.py; (e) remove `polygonize` from `checker/__init__.py`. Target-set hash must not change: gate by `target_set_hash` equality on fixtures (rings are pinned equal by tests/test_cuda_checker.py). Then `bitquad.py`, `certify.py` follow (certify becomes `tests/oracles/certify_oracle.py`; riga.py:226 uses GPU masks check). Effort 3.
9. Delete dead src (after confirming no test needs it): `checker/volume.py: topology_any` (+ its test in tests/test_volume.py); `eval/metrics.py: nsd, ged` public wrappers (use `_nsd_best`, `_ged` internally); `targets/augment.py: INVERSE`. Effort 1.
10. Decide Controller: if the owner keeps FastGRPO as the only trainer, archive `train/callbacks.py: Controller, audit_alarm`, `eval/audit.py`, the `recorded` hook in `train/rewards.py`, tests in `tests/test_train_trl.py` (Controller cases) and `test_breakit_*` bits to `archive/controller/`. Keep `TraceStep`, move it to `train/trace.py`. Alternative: wire it into train_rl as the held-out hack alarm. This is an owner decision (see questions). Effort 2.
11. Drop LFM/SmolVLM paths: `train/configs.py` MODELS/MAX_COMPLETION keep only `qwen3-vl-2b` (and new <=4B entries as added); remove smolvlm `size`/`use_vllm` branches and LFM sampling; strip LFM text in `geometry/text.py`, `train/sample.py`; `Makefile:28` default CONFIG -> `experiments/phase2/configs/rl_shifted.yaml`; `preflight.HF_RATE`, `project_r1.R1`; archive `configs/{r1_lfm,r0,r0_qwen3vl_arm,r1_qwen3vl2b,r2_qwen3vl2b,c2_*,p1_rl,rx8,checker_gpu,r3_gpu_test,r4_speed}.yaml` to `archive/phase2_runs/configs/` (keep the ledger as is), keeping `r5_shared, rl_full, rl_shifted`. Check `tests/test_train_configs.py`, `test_preflight.py`, `test_vast.py` for config names first (grep before move). Effort 3.
12. Phase 2 one-off scripts to archive after item 11: `project_r1.py`, `c2_analyze.py`, `r0_parity.py`, `hf_rate.py`, `hf_sample.py` (keep `k3_harness.py`, `parity.py`, `probes.py`, `profile.py`). Tests: `tests/test_p2_hf_eval.py`, `test_p2_scripts.py`, `test_p2_probes.py` need case removal. Effort 3.
13. Phase 2 lab: promote what tests use, archive the rest. Move `lab/B2/{analyze,posthoc,score}.py` and `lab/B3/power_c2.py` into `src/objdel/eval/` (signal already holds promoted parts) or `tests/helpers/`, and make `test_p2_signal.py`, `test_p2_probes.py` import them there; move `lab/A3/vast_snapshot.json` to `tests/fixtures/` and the digest source `lab/B5/Dockerfile.phase2` to a pinned constant in `experiments/phase2/docker/` (test_vast.py:31,789). Then `git mv experiments/phase2/lab archive/phase2_lab`. Effort 3.
14. Docs of history: move phase1/phase2 PLAN/EXECUTION/JUDGE/impl/perf md (about 70 files) to `archive/<phase>_docs/`; keep RESULTS.md files (README cites them) in `docs/research/evidence/`.
15. Tidy: rename `train_sft.preload` -> `load_items`; remove trivial unused vars listed in F13; add `vulture` to dev extras with a whitelist (`tools/vulture_whitelist.py` for TrainerCallback hooks, nanobind-registered names) and a `make dead` target (min-confidence 80; a 60% scan is advisory) so dead code does not return. Effort 1.

Suggested order: 1, 2, 3, 6, 9, then 7, 4, 5, then 11, 12, 13, 14 (config and script archive after the new run plan is fixed), then 8 and 10 in phase D/E.
Expected result: tracked files 802 -> about 330 outside archive (phase1/lab 276 + phase0 ~25 + phase2/lab 130 + phase1/phase2 md ~70 + one-off scripts); src stays 7.2k lines minus about 600 (Controller, audit, wrappers).

## Open questions (owner)
1. Keep `Controller` (held-out hack alarm and best-checkpoint) and `eval/audit.py`, or archive them because FastGRPO does not use them?
2. Archive mechanism: OK with the hybrid (tag before the move, then `archive/` folders with README, delete pure duplicates)? Tags are a git change: owner must run `git tag archive/pre-clean-2026-11-08`.
3. Is LFM2.5-VL-450M permanently retired (item 11)? The HF repo `svenk029/objdel-sft` holds `r1-lfm`; keep its `hf_push` license text or not?
4. Is the CPU exam-metrics pool (`objdel.eval.metrics`) to be ported to CUDA in phase D? If yes, `serve.py` pool code goes then, not now.
5. Should the numpy polygonizer stay importable for the targets builder on a Mac (no GPU) until the CPU build of `check_masks` is accepted as the fallback? The CPU build of `_gpucheck` is the same source, so (a) should hold; confirm that `make test` on macOS stays the test path.
6. Phase 0 evidence (3 MB csv/json): keep under `docs/research/evidence` or `archive/phase0`?

## Baseline
`make test` equivalent (`pytest -q -m "not slow"`) was started to confirm green before any change; it is slow (about 14 % after 3 minutes, 0 failures so far). Result file: `docs/implementation2.0/notes/CL1_baseline_pytest.txt`. Not complete when this note was written; re-run before step 1.
