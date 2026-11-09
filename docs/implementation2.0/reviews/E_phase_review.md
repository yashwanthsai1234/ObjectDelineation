# Review of Phase E plan (docs/implementation2.0/E.md)

Reviewer: fresh critical review, 2026-11-08. Base commit read: be52b1c. REPO = the ObjectDelineation repo.
Tags: [code: file:line] read in source this session; [measured: cmd] run this session; [web: url] read this session;
[estimate] not measured.

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | Phase E is about package layout. No NVIDIA catalog skill applies, so I recommend no install. I did not run the catalog CLI again because the plan already cites a snapshot (nv_catalog.json). |
| engineering:code-review | I checked each move table against the real imports. That found items 1, 2, 3, 14 and 18. |
| engineering:debug | I traced what each move breaks at runtime: paths built from `__file__` (item 9), module names written as strings (item 14), hard-coded test paths (item 8), and the preflight loader (item 4). |
| engineering:system-design | I checked whether the layer order holds for every caller (items 1, 2, 15, 27) and whether the infra code belongs inside the wheel (item 32). |
| engineering:architecture | ADR-level choices. I propose pydantic-settings, which is already in the image lock, in place of hand-written argparse and `to_argv` (item 26). I name the missing stage types (item 12). |
| engineering:documentation | Fixes must be exact and checkable. Wording fix for serve.py (item 27). Suggest the `broken_contract_guidance` option that points to the layer docs (item 16). |
| huggingface-skills:huggingface-best | No model choice in Phase E. I checked only that the per-model sampling values stay in one place (item 10). I ran no leaderboard query. |
| engineering:testing-strategy | Golden tests that cannot see the change they guard (item 6). Goldens that cannot match across CI platforms (item 7). Weak statistical gate (item 11). A coverage gate that can drop (item 29). |
| engineering:tech-debt | Duplicate sources the plan leaves or adds: regions (item 13), two vLLM sampling paths (item 10), the interim lint gaps (item 17). |
| WebSearch / WebFetch / PyPI JSON | import-linter versions and release notes (item 16), layers syntax (item 17), ruff banned-api page (item 3), pydantic-settings CLI support (item 26), pydantic latest version (2.14.3 is confirmed). |

## Verified as correct (no action)

- 33 of 63 `tests/test_*.py` import `experiments` [measured: grep]. 41 host-sync call sites in src [measured: grep -roE].
- LR sources: configs.py:189 (1e-5 / 1e-6), prereg.py:63-64, g_probe.py:86 [code]. Reward weights: spec.py:28 and signal.py:31 [code].
- train_rl.py: `--procs` at line 35, `fork_pool` at line 50, `pool.shutdown()` at line 61, and the pool gets no job [code].
- The TRL 1.17.1 method lines (1482, 1591, 1648, 2286, 2429, 3141) match. `split_tensor_dict`, `split_pixel_values_by_grid` and `unsplit_pixel_values_by_grid` are in `trl/trainer/utils.py` (1052, 1158, 1224) [measured: grep .venv].
- IMPORT_CHECK at policy.py:69 and `chain()` at policy.py:422 [code]. The baseline is 1447 passed / 6 skipped / 31 deselected in 523 s [measured: CL1_baseline_pytest.txt].
- pydantic 2.14.2 is in vllm_native.txt. fastapi 0.136.3, uvicorn 0.54.0 and httpx 0.28.1 are in the lock [code: locks/vllm_native.txt]. PyPI latest pydantic is 2.14.3 [measured: PyPI JSON].

## Findings (most severe first)

1. **blocking** | E5 move table, `rl/rows.py` row (E.md lines 395-400) |
   **Problem:** the plan puts `prompts` (and with it `records`, `_verified`, `cache_loader`) in `rl/rows.py`. But `select_ckpt` calls `D.prompts` [code: select_ckpt.py:184], and the plan moves select_ckpt to `sft/select.py`. `final_eval` calls `dataset.prompts` [code: final_eval.py:72], and the plan moves that logic to `policy/evaluate.py`. `eval/audit.py` uses `cache_loader`. So sft would import rl, which breaks the `rl | sft` contract. policy would import rl, which breaks the layer order. The table also leaves functions with no place: `view_rows`, `_no_holes`, `_contexts`, `_inside` (SFT views, used by train_sft and build_views), `prompt_rows`, `_tokens`, `_caches`, `_voting`, `_image`, `_check`, `_pixels_sha`, `_jobs` and `DISAGREEMENT_OVERSAMPLE` (used by final_eval, g_probe, train_rl) [code: dataset.py:86-392].
   **Fix:** add `policy/rows.py` with `prompts`, `records`, `_verified`, `_image`, `_check`, `_caches`, `_voting`, `cache_loader` and `DISAGREEMENT_OVERSAMPLE`. policy may import data and `policy.prompt`. Keep only `grpo_rows`, `grpo_rows_views`, `prompt_rows`, `OpenOnce` and `_ViewRows` in `rl/rows.py`. `sft/rows.py` gets `sft_rows`, `view_rows`, `_no_holes`, `_contexts`, `_inside` and `_tokens`. `data/prepare.py` gets `prepare`, `_prepare`, `_jobs` and `_pixels_sha`. List every top-level name of dataset.py in the table. Add a test: `set(dir(old_module))` before the move equals the union of the new modules' names.

2. **blocking** | E3b, `verify/parity.py` row (E.md line 362) |
   **Problem:** `experiments/cuda_checker/compare.py` adds `tests/` to `sys.path` [code: compare.py:24-26]. It then imports `oracles.raster_oracle`, `oracles.validity_oracle` [code: compare.py:139-140, 253-254] and `test_topology.random_masks` [code: compare.py:278]. If this code moves into `src/objdel/verify/parity.py`, the library depends on the tests. `test_src_never_imports_tests` checks only names that start with `tests`, so it misses `oracles.*` and `test_topology`. import-linter cannot see them either.
   **Fix:** do not move compare.py into src. Move it to `tests/oracles/cuda_parity.py`, because it is a test harness. Its 5 test users and the box gate already run under pytest. Add an AST test: 0 `sys.path` changes and 0 imports of `oracles` or `test_*` under `src/objdel`.

3. **blocking** | E0 ruff banned-api table and E3a (E.md lines 263-279, 348-350) |
   **Problem:** `make lint` cannot pass with the planned exceptions.
   (a) `tests/test_train_fast.py` imports `trl.trainer.grpo_trainer` at lines 427, 670 and 753, and `trl.generation.vllm_generation` at line 483 [code]. E0 bans both modules.
   (b) 13 test files and `tests/oracles/{raster,validity}_oracle.py` import `objdel._gpucheck` directly [measured: grep -l, 13 files]. Library files that do the same: `checker/__init__.py` and `geometry/shape.py:31`. E3a rewrites only `reward/cache.py:88`, `reward/score.py:41` and `train/rewards.py:60`.
   (c) In ruff, per-file-ignores work per rule code. A file freed from TID251 for TRL is also freed from the `experiments` ban, so TID251 cannot carry rule G8 for test files.
   **Fix:**
   (a) In the E3a move commit, rewrite all 13 tests, the 2 oracles and compare.py to import from `objdel.geom.engine`.
   (b) Add per-file-ignores for TID251: `tests/test_train_fast.py` (comment: "removed in E5"), then `tests/rl/test_trl_adapter*.py`.
   (c) Enforce "tests never import experiments" with an AST test in `tests/core/test_layers.py`, not with TID251.
   (d) Run `ruff check` once in E0 and paste the 0-error output into the commit message.

4. **blocking** | E1.5 and the E6 gate (E.md lines 326-328, 572) |
   **Problem:** `preflight.load` runs `TS.RunArgs(**raw)` [code: preflight.py:49-56]. `RunArgs` has no `steps` field [code: train_sft.py:77-101], so a `steps:` plan such as rl_shifted.yaml raises TypeError. The rl_shifted.yaml header itself says to preflight `r1_qwen3vl2b.yaml`. So E1.5 breaks the default `make preflight`. The E6 gate (`make preflight CONFIG=specs/rl_shifted.yaml`) has no task that rewrites `preflight.load`, `hashes` or `refusal`. `refusal` gates every launch [code: policy.py:313-316].
   **Fix:** drop E1.5, or set the default to `experiments/phase2/configs/r1_qwen3vl2b.yaml`. Add a task E6.8: `preflight` takes a `RunSpec`. Data and view hashes come from the spec stages. The stamp is keyed on the spec path. `policy.launch` calls `refusal(spec_path)`. Tests: preflight passes on `specs/rl_shifted.yaml` with fixture data and refuses an edited spec.

5. **major** | E1.1 (E.md lines 308-316) |
   **Problem:** the plan deletes `lr_full` and `lr_lora` from `prereg.G_PROBE`. The prereg.py header says: "a rule is never changed after a run sees data (a later re-specification gets a new name and an AMEND note)" [code: prereg.py:1-4]. The deletion also breaks `tests/test_p2_breakit_runners.py:60`, which asserts both values. The plan does not list that test.
   **Fix:** keep G_PROBE unchanged, because it is the pre-registered record. Add `rl_lr(model) -> float` in configs (5e-5 LoRA / 5e-6 full) and use it in `grpo_config` and `make_trainer`. Add `tests/experiments/test_prereg_lr.py`: `rl_lr(m) == G_PROBE["lr_lora" if lora(m) else "lr_full"]` for every model. G7 then counts code sources, not the prereg record.

6. **major** | Golden 4 and the E1.1 result check (E.md lines 290-291, 316) |
   **Problem:** the loss and the gradient norm of one update do not depend on the learning rate. So golden 4 cannot see the LR change it is meant to guard ("golden 4 unchanged"). It also pins the E2 override order only where that order changes the loss.
   **Fix:** add golden 6, the resolved `GRPOConfig.to_dict()` (without output_dir and logging_dir) for the rl_shifted stage and for a G-probe. Add to golden 4 the L2 norm of the parameter change after the step.

7. **major** | Goldens 1, 3 and 4 against CI (E.md lines 281-292; .github/workflows/tests.yml) |
   **Problem:** CI runs on 3 operating systems plus a Linux torch-CPU job [code: tests.yml]. Golden 1 hashes the numpy, OpenCV and scipy versions [code: spec.py:50]. numpy and scipy are not pinned (`numpy>=2.0`, `scipy>=1.13` [code: pyproject.toml]), so every resolver update changes golden 1. Golden 3 stores float hex from the CPU build, and clang and gcc on x86 and arm can give different bits [estimate]. Golden 4 depends on the CPU kernels of torch. So bit-equal goldens recorded on the Mac fail in CI.
   **Fix:** for golden 1, store `sha256(repr(astuple(REWARD_V4), min_area probes))` without the versions, and skip the full fingerprint check when the versions differ from the recorded ones. Record goldens 3 and 4 per platform (`invariance.<sys.platform>-<machine>.json`), with the recording script in `tests/goldens/phase_e/README`. Skip on a platform with no recorded file. The Mac run is the binding gate.

8. **major** | E8.3 test move against hard-coded test paths (E.md lines 589-591, risk row 1) |
   **Problem:** the plan counts only `experiments.phase*` strings (32). Test paths are also hard-coded in these places:
   - `pyproject.toml:58` (cibuildwheel test-command)
   - `tests.yml:54` (gates job)
   - `preflight.py:29-34` (CHECKER_TESTS, QUICK_TESTS, DRY_RUN) and `:212`
   - the live plan `rl_full.yaml:28,36-37`
   - `checker_gpu.yaml`, `p1_rl.yaml` and `r3_gpu_test.yaml`
   - `policy.chain` docstring (tests/test_gpu_scoring.py)

   Also, CI runs `ruff check` directly, not `make lint` [code: tests.yml:31], so lint-imports would never run in CI.
   **Fix:** E0 adds a `.venv/bin/lint-imports` step to tests.yml. E8 adds a test: every `tests/[\w/]+\.py(::\w+)?` string in pyproject.toml, .github/, experiments/, specs/ and Makefile names a file that exists. Add this item to the E8 grep gate.

9. **major** | E8.1 move of `secrets.py` and `vast/policy.py` (E.md lines 576-581) |
   **Problem:** `secrets.py` sets `ROOT = Path(__file__).resolve().parents[2]` [code: secrets.py:14], and `policy.py` sets `ROOT = HERE.parents[2]` [code: policy.py:46-47]. After the move to `src/objdel/infra/…`, these resolve to `src/` or `src/objdel/`. Then `.secrets`, `git ls-files`, rsync roots and `.vast-allowlist` resolve to the wrong place. `tests/conftest.py:16-17` forces `OBJDEL_SECRETS` to a fake path, so no test sees the change.
   **Fix:** one `objdel.infra.paths.repo_root()` (env `OBJDEL_REPO`, else `git rev-parse --show-toplevel`, else an error). Replace every `__file__.parents[n]` in moved files with it. Add a test that calls `repo_root()` and checks for `pyproject.toml`, and a test that `secrets.path()` is `repo_root()/".secrets"` when the env var is unset.

10. **major** | E4.1 `SamplingParams`, `VllmSampler` and duplicate samplers (E.md lines 495-513) |
    **Problem:**
    (a) The dataclass drops `min_p`, `repetition_penalty`, `presence_penalty` and `frequency_penalty`. `configs.sampling` sets all of them on purpose, to override the model card defaults (D13) [code: configs.py:108-114]. So the G10 check "params equal" cannot pass.
    (b) The plan does not list what VllmSampler must carry over from sample_vllm.py: `mm_processor_kwargs`, `request_seed` per item, `VLLM_USE_FLASHINFER_SAMPLER=0` set before the vLLM import, `check_memory`, `limit_mm_per_prompt` and `enable_prefix_caching` [code: sample_vllm.py:39-101].
    (c) The plan keeps sample_vllm.py, hf_sample.py and select_ckpt's child-process bridge beside the new sampler. That makes two vLLM paths, against the owner's archive rule.
    (d) sample_vllm.py "never imports objdel": it runs in the phase2 image's `/opt/vllm` venv (vllm.txt lock: torch 2.13.0, numba 0.65.0). An in-process sampler works only in the vllm_native image.
    **Fix:** `params_for` returns `configs.sampling(model, "vllm") | {"n", "seed", "stop_token_ids"}` as a plain dict, with no vLLM import. Move `requests`, `request_seed`, `ENV` and `check_memory` into `policy/sampler.py`. Make `python -I -m objdel.policy.sampler` the child-process entry point that select_ckpt keeps (D13 isolation, timeouts, resume). Archive sample_vllm.py and hf_sample.py in the same commit. State "vllm_native image only" in ADR 0015, and add `import vllm, objdel.policy.sampler` to IMPORT_CHECK.

11. **major** | E4 box gate and §8 cost (E.md lines 528-531, 702-703) |
    **Problem:** "the CI of the difference contains 0" is a weak test, because a small or noisy run passes it. It also needs both the trainer path and the sampler path on the same checkpoint, which is an extra eval. That contradicts §8, which says E4 "adds no extra eval".
    **Fix:** make the main gate deterministic. vLLM `prompt_token_ids` (sample_vllm's `--logprobs` path already records them) must equal the TRL-path prompt ids on 64 rows, and the `SamplingParams` repr must equal the one TRL builds. For the statistical part, use an equivalence test (TOST) with margin ±0.01 IoU on val-shifted, with nodule-level bootstrap, and state n. Add the extra eval time to §8. One val eval is about 9 min [estimate: half of the "~18 min of evals" in the rl_shifted.yaml header].

12. **major** | E6 Stage union (E.md lines 539-560) and E2.4 g_probe |
    **Problem:** the ported `r1_qwen3vl2b.yaml` uses the `stages:` schema (train, select, probe, with `sampler`, `branches` and `vision`) [code: policy.py:429-440; r1_qwen3vl2b.yaml]. The union has no SelectStage and no ProbeStage. g_probe has `--probe`, `--items` and `--select`, which RlStage does not have [code: g_probe.py:131-145]. So "g_probe calls build_trainer" needs a fake RlStage.
    **Fix:** add `SelectStage` (fields from `chain()`: sampler, items, vision, rotated, rescaled) and `ProbeStage` (probe, items, select | ckpt, compile, shared_prefix, plain, trace). Give `build_trainer` explicit typed arguments (model, out, roots, flags). RlStage and ProbeStage then map onto them. Add r1_qwen3vl2b to golden 5.

13. **major** | E6 `ComputeSpec.regions` (E.md line 548) |
    **Problem:** BRIEF §3 says North America hosts only. The rule today is a fixed constant, `REGION = ("US", "CA")` [code: policy.py:75]. A spec field adds a second source and a setting the owner rule does not allow. It also goes against the "single source" goal of E1.
    **Fix:** remove `regions` from ComputeSpec, so `extra="forbid"` rejects it. Make `infra.vast.offers.REGION` the only source. Add a test that a spec with a `regions` key fails.

14. **major** | Rule R-1 "move = git mv + import rewrite" and serve.SCORERS (E.md lines 207-209, 361) |
    **Problem:** some module names are strings, so an import rewrite does not change them. `serve.SCORERS = ("objdel.reward", "objdel.eval.metrics", "objdel.eval.probes")` [code: serve.py:52] is passed to `importlib.import_module`. After E3, verify/serve.py would name `objdel.data.probes`. That is an upward dependency (verify -> data) that import-linter cannot see. Other strings: train_sft.py:201 and test_p2_train.py:68,81 (monkeypatch targets).
    **Fix:** SCORERS becomes a parameter. Callers pass `modules=` (fork_pool already accepts it [code: serve.py:165]). E0 adds a test: every string literal that matches `objdel(\.\w+)+` in src, experiments and configs resolves with `importlib.util.find_spec`.

15. **major** | E8.1-2 IMPORT_CHECK and k3_harness (E.md lines 580-584) |
    **Problem:** k3_harness is a live CLI. Its import chain is `experiments.phase1.reward_gates` [code: k3_harness.py:33], then `tests/oracles/reward_harness.py` [code: reward_gates.py:9]. `code_files()` stops when tests/oracles is missing [code: policy.py:375-382]. E8 changes the IMPORT_CHECK string but keeps this chain, and phase1 is not in the E8 group list. So G8's "IMPORT_CHECK imports 0 tests" hides a real dependency on the tests.
    **Fix:** move `digest` and `make_jobs` into `verify` (or move k3_harness to `tests/experiments/`, since it is a test harness). Then remove the tests/oracles guard in `code_files`. Add a G8 test: run `python -I -X importtime -m experiments.rl.k3_harness --help` and check that no `tests.` module is loaded.

16. **major** | §2 and E0 import-linter pin (E.md lines 87-92, 232) |
    **Problem:** the plan says 2.11 (2026-03-06) is the current release and that 2.12 is not confirmed. That is wrong. PyPI lists 2.13 (2026-07-03), 2.14 (2026-08-28), 2.15 (2026-09-04) and 2.16 (2026-11-06) [measured: pypi.org/pypi/import-linter/json]. 2.16 makes validation stricter: unknown contract options and repeated modules are errors, and an empty config fails. 2.14 adds `broken_contract_guidance` [web: https://import-linter.readthedocs.io/en/latest/release_notes/].
    **Fix:** pin `import-linter==2.16`. Add `broken_contract_guidance = "see docs/architecture/layers.md"` to each contract. Run the E0 negative tests again on 2.16.

17. **major** | E0 contract "valid on day 1" (E.md lines 238-261) |
    **Problem:** every layer is optional and none exists at E0, so the layers contract checks nothing. The current cycle (eval <-> train) is not guarded. The `pure-core` source modules do not exist yet, so the check passes with nothing to check or fails as a missing module (not tested). The docs show `(foo)` and `a | b`, but not the two together in `(rl) | (sft)` [web: https://import-linter.readthedocs.io/en/stable/contract_types/layers/].
    **Fix:** E0 adds interim contracts on the current packages:
    - forbidden: `objdel.{checker,geometry,reward,targets,data}` must not import the ML stack.
    - layers: `["train", "eval", "reward | targets", "checker | geometry"]`, with `ignore_imports = ["objdel.eval.vision -> objdel.train.dataset"]`. E3 removes that ignore.

    Add `pure-core` sources only as each package appears. In E0, test `(x) | (y)` on a scratch package. If the syntax is rejected, use plain `rl | sft` from E5.

18. **minor** | E3a `geom/engine.py` (E.md line 342) |
    **Problem:** `from objdel._gpucheck import *` fails ruff F403 and F405, because "F" is selected [code: pyproject.toml `[tool.ruff.lint]`].
    **Fix:** import the names explicitly and set `__all__ = ("DEFECTS", "score_polygons", "check_polygons", "check_masks", "build", "device_count", ...)`. List every name the 13 tests use.

19. **minor** | E1 gate grep (E.md lines 330-331) |
    **Problem:** `grep "5e-5\|5e-6"` also matches the SFT literal (configs.py:166) and docstrings (configs.py:32 and 156, g_probe.py:6) [measured: grep]. So "exactly one place" cannot be measured.
    **Fix:** replace the grep with the test in item 5 plus an AST check: `G_PROBE["lr_*"]` is not read anywhere in src or experiments.

20. **minor** | E1.1 step "read configs.grpo_config(...).learning_rate before the overrides" |
    **Problem:** this builds a second GRPOConfig (TrainingArguments `__post_init__`) just to read one number.
    **Fix:** use `configs.rl_lr(model)` (item 5).

21. **minor** | E5.2 REQUIRED list and E5.3 `_prepare_inputs` copy (E.md lines 410-435) |
    **Problem:** the copied body reads private trainer state that the check does not list: `self._step` (grpo_trainer.py:995, 1582), `self._buffered_inputs` (998) and `self.num_iterations` [measured: grep]. The plan describes only the train branch. The copy must keep TRL's eval `else` branch and our fallback guard [code: fast.py:553-556].
    **Fix:** add `_step`, `_buffered_inputs` and `num_iterations` to REQUIRED (an attribute check on a built tiny trainer) and to the PRIVATE scan set. In the equivalence test, check eval mode too.

22. **minor** | `test_no_private_trl_names_outside_adapter` (E.md lines 295-298) |
    **Problem:** a scan of all attribute names flags `_metrics`, `_generate` and `_step` anywhere in src, including our own classes.
    **Fix:** flag only attribute access on `self` in classes that inherit GRPOTrainer or VLLMGeneration, plus `trainer.<name>` and `super().<name>`. Keep TID251 for module imports.

23. **minor** | G11 host-sync count (E.md lines 186-187, 299-300) |
    **Problem:** the regex also counts numpy `.tolist()` and `.numpy()` in verify and E7 code, where no GPU sync happens. It would block JSON conversion in `verify/api.py`.
    **Fix:** count only in modules that import torch (policy, rl, policy.caches), recorded per file. Store the per-file baseline in `tests/goldens/phase_e/host_syncs.json`.

24. **minor** | E2.1 `RlStage` defaults (E.md lines 454-470) |
    **Problem:** `save_steps: int | None = None`, but the train_rl CLI default is 50 [code: train_rl.py:31]. A spec without the field behaves differently. On a frozen model, a `dict` field makes `hash()` raise.
    **Fix:** set `save_steps: int = Field(50, gt=0)`, and add a rule `save_steps <= steps` when set. Type `grpo` as `Mapping[str, JsonValue]` and never hash the stage, or store it as a JSON string that is validated on load.

25. **minor** | pydantic as a base dependency (E.md lines 450, 690, 728) |
    **Problem:** the box installs objdel with `--no-deps` [code: policy.py:399]. The phase2 image lock (train.txt) has no pydantic [measured: grep -c = 0]. If `.vast-allowlist` names that image, `import objdel.core.spec` fails on a paid box. An exact `==` pin in base deps also clashes with any future vLLM bump.
    **Fix:** use `pydantic>=2.14,<3` in base deps and the exact version only in the locks. Add `import pydantic` to IMPORT_CHECK and image_check.sh. Write down that live runs use the vllm_native image, or archive the phase2 image in CLEAN.

26. **minor (better library)** | E6 CLI and `to_argv` (E.md lines 535-537, 554, 562) |
    **Problem:** the plan rejects tyro because it is in no lock. But pydantic-settings 2.16.2 is already in vllm_native.txt [code: locks/vllm_native.txt:651]. It builds a CLI from a pydantic model (`cli_parse_args`, `CliApp.run`, `cli_kebab_case`, subcommands) [web: https://pydantic.dev/docs/validation/latest/concepts/pydantic_settings/].
    **Fix:** build each CLI with `CliApp.run(RlStage)`, so there is no second argparse table. Generate `to_argv` from `model_dump` with one kebab-case rule. Golden 5 still pins the strings. Pin `pydantic-settings==2.16.2` in the same extra.

27. **minor** | E3b `verify/serve.py` note "deleted by Phase D" (E.md line 361; §3.3 row) |
    **Problem:** serve.py also holds `preload`, `loader`, STORE and VIEW_SEP. The training reward reads these [code: g_probe.py:91; train_rl.py:49; dataset.py:244]. Only the fork pool can go.
    **Fix:** split into `verify/store.py` (preload, loader, VIEW_SEP) and `verify/pool.py` (fork_pool, pool_size). Phase D deletes only pool.py.

28. **minor** | E8.1 `code_files()` move |
    **Problem:** `git ls-files -co` sends untracked files to rented boxes, and the current list starts with `.claude/worktrees/…` and `.claude/settings.json` [code: policy.py:377; code_files.txt lines 1-3].
    **Fix:** when the function moves to `infra/vast/remote.py`, add `^\.claude/|^archive/` to the drop pattern. Add a test that no `.claude/` path is listed.

29. **minor** | E7 HTTP face and the coverage gate |
    **Problem:** the CI train job fails below 90 % coverage [code: tests.yml train job]. `verify/http.py` is skipped when the `serve` extra is missing, so its lines are not run and coverage drops.
    **Fix:** install `.[dev,data,train,serve]` in the train job, or add `omit = ["*/verify/http.py"]` with a comment.

30. **minor** | E5.2 env_check integration and G13 (E.md lines 425-427, 675) |
    **Problem:** a separate `python -I -c` line does not write into `env_check.json`. Each check inside the heredoc does [code: env_check.sh: the `check()` wrapper].
    **Fix:** add `@check("trl_private_names")` to the existing heredoc, so G13's evidence file holds the result.

31. **minor (GPU first)** | E3a/E3c polygonizer seam (E.md lines 345, 375) |
    **Problem:** `targets/build.py` imports `polygonize` from `objdel.checker` [code: build.py:53]. After E3 it would import the numpy module by name. The owner rule says the polygonizer must be CUDA. Phase D then edits every caller.
    **Fix:** in E3a, expose one `geom.polygonize(mask)` function that calls the numpy code today. Phase D changes only that function to `engine.check_masks` when `engine.build() == "cuda"`, behind a parity test.

32. **minor (architecture)** | E8.1 infra inside the wheel |
    **Problem:** `wheel.packages = ["src/objdel"]` [code: pyproject.toml]. So the published objdel wheel would carry the Vast rental, ssh and secrets code.
    **Fix:** add `wheel.exclude = ["**.cpp", "objdel/infra/**"]`, or keep infra as a separate top-level package outside the wheel. Record the choice in ADR 0012.
