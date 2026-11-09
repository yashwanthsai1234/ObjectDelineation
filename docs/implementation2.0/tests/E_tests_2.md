# Phase E test cases, set 2 (adversarial, property, performance, failure modes)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA skill needed; Phase E has no GPU work (spec constraint).
- engineering:testing-strategy: pyramid split; most cases are unit/property on CPU, few integration, GPU cases only skip-marked.
- mattpocock-skills:tdd: tests sit at public seams only (load_spec, to_argv, Verifier.score, check_private_names, evaluate, lint-imports); expected values come from literals, goldens or an independent oracle, never from re-running the code under test.
- engineering:code-review, debug, system-design, architecture, documentation, huggingface-skills:huggingface-best: read for adversarial thinking only; they changed no case.

Source: only phases/E.spec.md (B1-B14, G1-G16) and the existing tests/ style (plain pytest, `_box`/`_compact` helpers, conftest tiny offline TRL model, fixtures skipped when data is absent).
Conventions: "CPU" = darwin-arm64 CPU build of `objdel._gpucheck`. `hx = float.hex`. Fixtures named `goldens` = `tests/goldens/phase_e/*`. Mark: U unit, P property, I integration, G GPU (skipped without GPU), E2E.

---------------------------------------------------------------------------------------------------------------
## A. Layers and imports (B1, B2, B3, B13, B14)

1. **layers_contract_is_exhaustive_and_clean** (I). Setup: parse `pyproject.toml`/`.importlinter`. Input: run `lint-imports` as a subprocess. Expect: exit 0; contract type `layers`, container `objdel`, layers exactly `infra > {rl|sft} > policy > data > verify > geom > core` (rl and sft as sibling layer: `rl | sft`), `exhaustive = true`, exhaustive_ignores == `["_gpucheck"]`, `ignore_imports` empty or absent, every contract has non-empty `broken_contract_guidance`. Why: a loosened contract (an added ignore) silently allows a layer breach.

2. **layers_contract_catches_planted_violation** (I). Setup: copy `src/objdel` to tmp_path, append `import objdel.rl` to `objdel/core/__init__.py`, run `lint-imports` against the copy. Expect: non-zero exit; output names `objdel.core` and `objdel.rl` and prints the guidance text. Repeat for (a) `sft` importing `rl`, (b) `verify` importing `data`, (c) a new top-level package `objdel.newpkg` (exhaustive must fail). Why: proves the linter really guards; a contract that passes on everything passes vacuously.

3. **old_packages_raise_ModuleNotFoundError** (U). Input: `importlib.import_module(n)` for each of `objdel.checker, geometry, reward, eval, targets, train, train.fast, train.dataset`. Expect: `ModuleNotFoundError` for each; also `objdel.__path__` dirs contain no `*.pyc`-only leftover (no `__pycache__` dir of a deleted package, since a stale `.pyc` in a package dir without source is not importable but stale dirs shadow namespace packages: assert `not (SRC/"objdel"/n).exists()`). Why: a leftover empty dir becomes a namespace package and `import objdel.train` would succeed silently.

4. **no_ml_stack_below_policy_subprocess** (I). For each `m` in submodules found by `pkgutil.walk_packages` of core, geom, verify, data: `python -I -c "import m, sys; bad=[x for x in ('torch','trl','transformers','vllm','peft','datasets','accelerate') if x in sys.modules]; sys.exit(bool(bad))"`. Expect exit 0 for every module, one process each (not one process importing all, so order cannot hide a leak). Also run with `sys.modules` poisoned: a meta-path finder that raises on those names; import must still succeed. Why: lazy imports inside functions are fine, top-level ones are not; the poisoned finder also catches indirect imports through a try/except ImportError.

5. **no_ml_stack_property_over_new_modules** (P). Hypothesis-free enumeration: for every `.py` under core/geom/verify/data, AST-walk all `Import/ImportFrom` at module level and inside `if TYPE_CHECKING` excluded; assert top-level name not in banned set. Why: catches a function-level import that someone later hoists.

6. **verify_names_no_higher_layer_in_strings** (U). AST + regex over every string constant in `objdel/verify/**`: no match of `objdel\.(policy|rl|sft|infra)`. Also in imports. Why: B13 explicitly includes strings (the old `fork_pool` built-in module list).

7. **only_engine_imports_gpucheck** (U). grep AST of src/tests/experiments: `_gpucheck` imported only in `objdel/geom/engine.py`. Plus set equality: every attribute `_gpucheck.<x>` / `from objdel._gpucheck import x` used anywhere in src, tests, experiments is in `objdel.geom.engine.__all__`, and every name in `__all__` exists on `objdel._gpucheck`. Why: dead or missing re-export.

8. **geom_engine_all_is_superset_after_mutation** (U, adversarial). Monkeypatch-free: temp copy of tests tree with `from objdel._gpucheck import score_polygons` added to a test; the checker script from case 7 must report it. Why: proves the checker itself detects.

9. **direction_no_experiments_tests_oracles_imports** (U). AST over `src/objdel`: no import of root `experiments|tests|oracles|test_*`; no `sys.path` mutation (`sys.path.insert/append/extend`, `sys.path[...] =`, `site.addsitedir`). Over `tests/**` excluding `tests/experiments/`: no `import experiments`/`from experiments`. Also string-form `importlib.import_module("experiments...")`. Why: string imports bypass the AST walk of imports.

10. **IMPORT_CHECK_string_is_clean** (U). `exec`-free: parse `IMPORT_CHECK` with `ast.parse` after splitting into statements; list imported modules. Expect: none starts with `tests` or `experiments`; includes `pydantic`, `vllm`, `objdel.policy.sampler`. Then run it on the Mac with `vllm` stubbed via a tmp `vllm/__init__.py` on PYTHONPATH (no -I here, document why): exit 0. Why: the box check must not need tests that are not shipped.

11. **cli_help_loads_no_test_modules** (I). For every live script under `experiments/**` (not `archive/`): `python -I -m experiments.x --help`; run with `-X importtime`-free approach: `python -I -c "import runpy,sys; sys.argv=[m,'--help']; try: runpy.run_module(m, run_name='__main__') except SystemExit as e: code=e.code; print(sorted(k for k in sys.modules if k.split('.')[0] in ('tests','oracles')))"`. Expect exit code 0 (`SystemExit(0)`), printed list `[]`. Timeout 20 s each. Why: argparse built at import time that pulls oracles breaks on boxes without tests.

12. **fork_pool_imports_exactly_given_modules** (U). `fork_pool(2, ["json"])` then submit `lambda`-free worker `sys.modules` probe; compare to `fork_pool(2)`: difference of loaded modules == {"json"} family only (assert `"objdel.policy" not in` workers' sys.modules when not passed). `fork_pool(1, ["no.such.module"])` raises `ModuleNotFoundError` in the parent at creation, not a hang in workers. Why: no hidden built-in list; failure mode of a bad name must be loud.

13. **pool_size_cgroup_parsing** (P). `pool_size(cgroup=tmp)` with files: `cpu.max` = "max 100000" -> os.cpu_count(); "200000 100000" -> 2; "150000 100000" -> 1 (floor, min 1); "50000 100000" -> 1; garbage / missing file / empty -> os.cpu_count() (no crash); cgroup v1 `cpu.cfs_quota_us=-1`. Property: result is int, `1 <= r <= os.cpu_count()` (or documented cap) for 200 random quota/period pairs (hypothesis `integers(1,10**9)`). Why: container sizing bug oversubscribes fork workers.

14. **every_name_has_a_home** (U). Load `tests/goldens/phase_e/names.json` (old module -> sorted top-level names, recorded at base commit; recorded not computed). For each old module compute the union of top-level names (ast, `Assign/FunctionDef/ClassDef/Import` aliases excluded) across its new modules from the rename map; expect `old == union | archived | renamed_old_names` and `union - old` listed in `added`. Also each `archived` name exists under `archive/<name>/`. Why: silent drops of helper functions in a split.

---------------------------------------------------------------------------------------------------------------
## B. TRL adapter (B4, G5, G6)

15. **check_private_names_passes_on_pin_fast** (U). Input: `check_private_names()` in a fresh subprocess timing the whole `import trl` + call. Expect: returns None, wall < 10 s (G5), and `trl.__version__ == "1.17.1"`. Skip with message if trl absent. Why: perf gate on the boot path of every run.

16. **check_private_names_negative_cases** (U). Each case uses `monkeypatch` on the real TRL objects and expects `TrlDriftError` whose `str` contains the missing/changed name AND `"1.17.1"`:
   a) `monkeypatch.setattr(trl, "__version__", "1.17.2")` (also "1.17.1.dev0", "1.17.10", "", None -> still TrlDriftError, never TypeError).
   b) `delattr(GRPOTrainer, name)` for each of the six overridden methods (6 parametrised cases; use `raising=True`).
   c) replace a method with `def f(self, other_name): ...` (changed params) for each of the six; also add an extra defaulted param; also reorder two params.
   d) replace a method with a function with identical signature but different body (so the SHA-256 differs). Also whitespace-only change to the source (SHA over `inspect.getsource` must differ; documents that formatting counts). 
   e) built trainer (tiny model) with `del trainer._metrics` (and `_step`, `_buffered_inputs`, `num_iterations` each).
   f) `inspect.getsource` raising `OSError` (frozen/pyc-only install) -> `TrlDriftError`, not a raw OSError.
   Why: the guard is only worth its negative cases; (f) is a real deployment failure mode.

17. **check_private_names_message_is_actionable** (U). For case 16b, assert message contains method name, "1.17.1", and the installed version string. Why: operators read this on a rented box.

18. **fastgrpo_init_checks_before_trl_init** (U). Monkeypatch `GRPOTrainer.__init__` with a sentinel that sets `called=True`, and `check_private_names` with a function raising `TrlDriftError`. `FastGRPO(...)` must raise `TrlDriftError` with `called is False`. Order of calls recorded in a list equals `["check", "trl_init"]` in the passing case. Why: otherwise TRL init crashes with an unreadable error first.

19. **private_names_only_in_adapter** (U). AST over `src/objdel/**` except `rl/trl_adapter.py`: no `Attribute` node with attr in the 11 names whose value is `self` in a class deriving `GRPOTrainer|VLLMGeneration`, `trainer`, `super()` call, or a name bound to a `trl` module (`trl.trainer.grpo_trainer.shuffle_sequence_dict`). Allowed: unrelated classes using `_step`. Fixture: a hand-made snippet where an unrelated class uses `self._step` (must pass) and one where `trainer._metrics` is used (must fail). Why: the checker rule has a false-positive edge; test the checker with both snippets.

20. **no_assignment_to_shuffle_sequence_dict** (U). After importing `objdel.rl` fully (all submodules, run `FastGRPO` construction and one `_prepare_inputs` on the tiny model), `trl.trainer.grpo_trainer.shuffle_sequence_dict is ORIGINAL` (captured via `inspect` before import in a subprocess). AST: zero `Assign/AugAssign/setattr` targeting that name. Also run an interleaved case: build `FastGRPO` then a plain `GRPOTrainer` in one process; the plain one's `_prepare_inputs` must use the original shuffle (compare micro-batch order to the TRL reference with a seeded RNG). Why: the global patch used to leak between trainers.

21. **fast_ordering_golden_microbatches_and_loss** (I). Tiny offline model, CPU, seed 0, `shared_prefix in (False, True)`, dr_grpo. Expect: micro-batch index lists (list of lists of sample ids) equal `goldens/phase_e/microbatches.json`; loss `hx` equals golden; shared-prefix on vs off loss within the recorded tolerance (golden stores both; assert each equals its own golden exactly). Why: B10.4 plus G6.

22. **fast_fallbacks_equal_trl_prepare_inputs** (I, parametrised). Cases: eval mode; loss_type != "dr_grpo"; entropy cut on (`top_entropy_quantile < 1`); `gradient_accumulation_steps % steps_per_generation != 0` (e.g. gas=3, spg=2). For each, run the FastGRPO `_prepare_inputs` and TRL's `GRPOTrainer._prepare_inputs` on identical inputs with the same torch seed (seed set directly before each call) and assert tensors `torch.equal` key by key and equal key sets. Adversarial extra: the boundary gas == spg and gas = 2*spg must NOT fall back (assert the Fast path was taken through a call counter on the adapter's ordering helper). Why: a wrong modulo test sends production into fallback (slow) or Fast path in the unsupported mode (wrong grads).

23. **fast_ordering_property_permutation** (P, hypothesis). Random lengths `lists(integers(1,200), min=8, max=64)` for a group structure; the Fast ordering of micro-batches must be (a) a permutation of all sample indices (each exactly once), (b) every group of 8 answers to one prompt stays in a micro-batch set that preserves loss-sum (loss computed with dr_grpo normalisation equals the unsorted loss to 1e-12 relative on a toy linear model), (c) deterministic for equal seeds. Why: sorting by length is where off-by-one drops samples.

24. **rl_loss_grad_param_delta_golden** (I). Tiny model, CPU, seed 0, one FastGRPO update: `hx(loss)`, `hx(grad_l2)`, `hx(param_delta_l2)` equal `goldens/phase_e/invariance.darwin-arm64.json`. On a platform without the file: `pytest.skip("no golden for <platform>")`. Why: B10.4 move-invariance; exactly bit-equal on the binding platform.

---------------------------------------------------------------------------------------------------------------
## C. Single sources and config (B5, B6)

25. **rl_lr_all_model_keys_match_prereg** (P). For every key in the model registry (`objdel.policy.models` keys, enumerated, not hardcoded): `rl_lr(m) == G_PROBE["lr_lora" if is_lora(m) else "lr_full"]`, with literals `5e-5` / `5e-6` asserted independently in the test. Unknown key -> raises `KeyError`/`ValueError` (assert it is not silently 5e-5). Why: a new model key silently inherits a wrong lr.

26. **no_code_reads_G_PROBE_lr_outside_prereg** (U). Regex over src/experiments/tests (except prereg.py and its own test): `G_PROBE\[["']lr_(full|lora)["']\]` has zero matches; also `G_PROBE.get("lr_`. Why: second source of truth.

27. **grpo_config_uses_rl_lr_when_none_given** (U). `grpo_config(m, out)` for each model: `.learning_rate == rl_lr(m)`; `grpo_config(m, out, learning_rate=1e-3).learning_rate == 1e-3`; `build_trainer(..., flags=TrainerFlags(lr=None))` -> `Built.lr == rl_lr(m)`; `lr=0.0`: documented behaviour must be "0.0 is honoured, not treated as None" (assert `Built.lr == 0.0`). Why: classic `lr or default` falsy bug.

28. **weights_single_source** (U). `WEIGHTS == dict(zip(TERMS, REWARD_V4.weights))` and `len(TERMS) == len(REWARD_V4.weights)`; `sum` equals recorded literal; no second dict literal: AST scan of src for dict displays whose keys are exactly the TERMS names (>=3 of them) outside `reward.py`. Why: drift between stats and reward.

29. **region_rule_constant** (U). `REGION == ("US","CA")` exactly; a `RunSpec` dict containing `compute.regions` raises `ValidationError`; `offers` filtering on a fixture list with hosts in "US", "CA", "DE", "us" (lower), "" and missing: only exact "US"/"CA" kept (document the lower-case expectation as the spec decides; the test pins current literal behaviour). Why: policy rule is money-relevant (memory: North America only).

30. **no_worker_pool_in_rl_training** (U). AST over `experiments/**` live RL scripts and `objdel/rl/**`: no `ProcessPoolExecutor|multiprocessing.Pool|fork_pool` calls. Runtime variant on tiny model: `run_rl` with `monkeypatch` of `concurrent.futures.ProcessPoolExecutor.__init__` raising; run completes. Why: B5 last bullet.

31. **build_trainer_signature_and_flags** (U). `inspect.signature(build_trainer)`: no parameter annotated `argparse.Namespace`, params named exactly (model, ckpt, rows, out, roots, steps, flags, **over); passing a Namespace as `flags` raises `TypeError`/`ValidationError`. `flags.plain=True` -> `type(trainer) is GRPOTrainer` (exact, not subclass); `False` -> `type(trainer) is FastGRPO`. `grep make_trainer` over src, experiments, tests: 0 hits (excluding this assertion's own file: build the pattern with string concat). Why: B6 and G14.

32. **grpo_config_goldens** (I). `build_trainer` resolved `GRPOConfig.to_dict()` minus keys {output_dir, logging_dir, run_name, `*_time*`, `hub_*` path keys; list recorded in golden} for `rl_shifted` train_rl stage and default G-probe `qwen3-vl-2b` equals `goldens/phase_e/grpo_config_{rl_shifted,gprobe}.json` exactly (json round trip; floats compared by `hx`). Adversarial: assert the golden's key-removal list is itself identical in test and golden (so the test cannot hide a changed key by widening the removal set). Why: removal-list creep.

---------------------------------------------------------------------------------------------------------------
## D. Run spec (B8)

33. **load_spec_every_file_loads** (U). Parametrise `specs/*.yaml` (glob at collection; assert glob non-empty and includes rl_shifted, rl_full, r5_shared, r1_qwen3vl2b). `load_spec(p)` returns `RunSpec`; `spec.__class__.model_config["frozen"]`; assigning `spec.name = "x"` raises `ValidationError`; hashing a stage is not required (`stages` contain dicts: assert `hash(stage)` raises `TypeError` OR is documented unsupported; pin whichever the spec says: "never hashed", so just assert no code path hashes: skip). Why: broken port of a YAML.

34. **load_spec_rejects_unknown_keys_everywhere** (P). Take each valid spec; for every nested mapping path (spec root, compute, each stage, `grpo`-sibling levels except `grpo` which is a free Mapping) inject key `"zz_unknown": 1` -> `pydantic.ValidationError` whose error loc ends with `zz_unknown`. Exception: inside `grpo` an arbitrary key is accepted (assert that, and that a non-JSON value such as a YAML timestamp/`set` raises). `compute.regions` raises. Why: extra="forbid" missing on one nested model is the typical hole; walking every path finds it.

35. **load_spec_structural_failures** (U, parametrised). Expect `ValidationError` (or the documented `ValueError` wrapper; one type, pinned) for: `schema_version: 2`, `0`, `"1"` (string; strict Literal int must reject? pin per spec: Literal[1] rejects "1" in strict, but pydantic lax mode may coerce: test documents actual behaviour and fails if it flips), missing schema_version; unknown `kind: "rlx"`; duplicate stage names (same name, different kinds too); empty `stages: []` (pin: allowed or not, per implementation; spec silent -> record as open question); `steps: 0`, `-1`; `save_steps: 0`; `save_steps > steps` (steps=10,save_steps=11) raises, `save_steps == steps` passes; RlStage default `save_steps == 50` when steps=10 raises (default 50 > steps 10): this is an adversarial case from the spec text; the test pins whichever the implementation chose and flags it. `ProbeStage` with both ckpt and select, and neither; `cap_hours: 0`, `-1`, `nan`, `inf`; empty file; YAML that is a list; non-UTF8 bytes; YAML anchors producing a duplicate key (duplicate mapping key must raise, not last-wins; PyYAML's default is last-wins, so the loader must be strict). Why: each is a past-run footgun; duplicate keys in YAML silently change a run.

36. **to_argv_goldens_string_for_string** (I). For rl_shifted, rl_full, r5_shared: `to_argv(stage, rd, data)` == `goldens/phase_e/chain_argv.json[spec][stage]` (list equality, not set), using the recorded rename map `{old_module: new_module}`; apply the map to the golden's `-m` token only. For r1_qwen3vl2b: select, probe, cmd stages equal; the SftStage fields equal `dataclasses.asdict(train_sft.RunArgs(...))` loaded from the golden JSON (the old class is gone; the golden holds the dict). Adversarial: assert the rename map has no entry that maps a module to itself and no two old modules to the same new one unless listed. Why: argv drift changes a paid run.

37. **to_argv_cmd_literal_and_substitution** (P). `CmdStage(argv=("echo","{rd}/a","{data}",  "{x}", "$HOME", "a b", ""))` with rd="R D", data="d": result == `["echo","R D/a","d","{x}","$HOME","a b",""]` (only `{rd}` and `{data}` replaced; unknown braces, spaces, shell metachar, empty string preserved; no shlex splitting; no env expansion). Hypothesis: for random text lacking the two tokens, `to_argv` is the identity. Text containing `{rd}{rd}` -> both replaced; `{{rd}}` -> pin behaviour (literal replacement gives `{R D}`; test records it). Why: injection/quoting bugs.

38. **cli_old_flags_equal_spec_stage** (I). For each live CLI: parse old flags -> stage S1; `--spec P --stage NAME` -> S2; `S1 == S2` (pydantic equality). Use the spec stage's own `to_argv` as the old-flag input. `--stage nope` -> exit code 2 with message listing available stage names; `--spec missing.yaml` -> exit != 0 with a clean message, no traceback of pydantic internals; `--spec` without `--stage` for a multi-stage spec -> error. Why: two entry paths drift.

39. **preflight_stamp_keyed_on_spec_sha** (I). Fixture data tmp dir. Preflight on `specs/rl_shifted.yaml` -> pass and stamp file written with key containing the path and SHA-256. Append a space/newline to a COPY of the spec (same path) -> launch refusal raised. Revert -> accepted again. Touch mtime only (no content change) -> still accepted. Copy the identical content to a different path -> refused (stamp keyed on the path too). Stamp file corrupted (truncated JSON) -> refusal, not crash. Why: money guard.

40. **run_spec_json_written** (I). A dry/tiny run writes `runs/<id>/run_spec.json` with keys spec, git_sha, trl, torch, transformers, reward_fingerprint; `git_sha` is 40 hex or `"unknown"` when no git (run with `PATH` lacking git; must not crash); `reward_fingerprint == fingerprint()`; the `spec` round-trips: `RunSpec.model_validate(json["spec"]) == load_spec(path)`. Run twice with same id -> second does not silently overwrite a different spec (raise or new id; pin). Why: reproducibility record.

41. **save_json_rejects_nan_and_handles_numpy** (U). `save_json({"a": float("nan")}, p)` raises `ValueError` and leaves NO partial file at `p` (atomic write: assert p absent and no `*.tmp` residue); `inf`, `-inf` same; nested nan in list/np array same; `np.float32(1.5)`, `np.int64(3)`, `np.bool_(True)`, `np.arange(3)` serialise to `1.5, 3, true, [0,1,2]`; `np.float64("nan")` raises; a pre-existing valid file is intact after a failed write. Why: crashes mid-run lose results; NaN in JSON is non-standard.

---------------------------------------------------------------------------------------------------------------
## E. Sampler and evaluation (B7, G10)

42. **params_for_equals_sampling_union** (P). For each model key, n in {1,8}, seed in {0, 2**31-1}, stop ids in {(), (1,), (1,2)}: `params_for(...) == sampling(m,"vllm") | {"n":n,"seed":seed,"stop_token_ids":list(s)}` and contains `min_p, repetition_penalty, presence_penalty, frequency_penalty`; result `stop_token_ids` is a `list` not tuple; input `sampling()` dict is not mutated (deepcopy before/after equality). Why: shared-dict aliasing bugs.

43. **request_seed_properties** (P, hypothesis). `request_seed(seed, item_id)`: deterministic; in `[0, 2**63)` (pin range from implementation, must be valid for vLLM); different item ids give different seeds over 10,000 random ids with collision count 0 (64-bit space; document tolerance); different base seeds give different seeds for the same id; `item_id` containing unicode/empty string works; stable across processes and `PYTHONHASHSEED` values (run in two subprocesses with `PYTHONHASHSEED=1,2`, equal output). Why: Python `hash()` use would break reproducibility across runs.

44. **requests_carry_seed_and_mm_kwargs** (U). Rows with and without `processor_kwargs`: each request dict has `seed == request_seed(seed, id)`; `mm_processor_kwargs` present iff the line has it; order of requests equals order of rows; duplicate item ids in rows -> raises or documented. Empty rows -> `[]`. Why: B7.

45. **ENV_set_before_vllm_import** (I). Subprocess with a stub `vllm` package that records `os.environ.get("VLLM_USE_FLASHINFER_SAMPLER")` at its import time: `import objdel.policy.sampler; VllmSampler(...)` constructs; recorded value "0". Case: user has `VLLM_USE_FLASHINFER_SAMPLER=1` already -> pin whether setdefault or override (spec says it is set to "0": override; assert "0"). Why: ordering bug is invisible until a box.

46. **sampler_import_is_light** (I). Subprocess: `import objdel.policy.sampler, objdel.policy.evaluate`; assert `trl.trainer`, `vllm` not in `sys.modules` (torch allowed). Why: B7 last bullet.

47. **hf_sampler_shape_invariants** (I, tiny offline model; CPU). For rows of k in {1,3} and n in {1,4}: one `Sample` per row, in row order, `len(texts) == len(token_ids) == len(truncated) == n`; `truncated[i] == truncated(token_ids[i], stop)`; `texts[i] == tokenizer.decode(token_ids[i], skip_special_tokens=True)` (independent oracle via the tokenizer); with `max_tokens=1` all rows truncated; determinism: same seed twice -> identical token ids; different seeds -> at least one difference over 8 samples at temperature > 0 (flaky guard: use temperature 1.5, n=8; allowed to skip if tiny model degenerate). Empty rows -> `[]`. Why: B7.

48. **truncated_function_table** (U/P). `truncated([], s)` is True? Pin by spec semantics "no stop token seen at the end = truncated": cases `([5,6,7],[7])->False`, `([5,6,8],[7])->True`, `([7,5],[7])->True` (stop not last), `([],[7])->True`, `([5],[])->True`, `([7,7],[7,8])->False`, tuple vs list input. Property: `truncated(ids+[s], stops)` is False for any `s in stops`. Why: this decides which rollouts are penalised as cut off; off-by-one at last position is the classic.

49. **check_memory_failure_mode** (U). `check_memory(lambda: (free, total), need)`: free >= need*total passes; free just below raises with a message naming free/total GB; `mem_get_info` raising `RuntimeError` -> clean error, not hang; `total == 0` -> no ZeroDivisionError. Why: pre-flight before vLLM load.

50. **sampler_main_child_process_cli** (I). `python -I -m objdel.policy.sampler --kind hf --ckpt <tiny> --rows rows.jsonl --out out.json` writes JSON with one entry per row; `--kind bogus` -> exit 2; missing rows file -> exit != 0 and nothing in out; `select_ckpt.sampler_cmd(...)` returns an argv whose parse via the sampler's own parser gives the same args (round trip through `argparse` of `main`). Child must not inherit `sys.path` entries from cwd (`-I`): run from a dir containing a malicious `json.py` and assert not imported. Why: planted-module defence per repo rules.

51. **evaluate_builds_no_trainer** (I). `monkeypatch` `GRPOTrainer.__init__`, `FastGRPO.__init__`, `VLLMGeneration._generate`/`GRPOTrainer._generate` to raise `AssertionError`; `evaluate(fake_sampler, rows, params, view_root)` completes with a fake Sampler returning scripted `Sample`s. Expect result dict whose per-item metrics equal values hand-computed from known box polygons (independent oracle: IoU of two overlapping 20x20 squares = 100/700 etc.). Why: G10.

52. **evaluate_adversarial_samples** (U). Fake sampler returns: all-empty texts; truncated True items (must score as the trainer path does: reward 0 / flagged per `ScoreItem.truncated`); a Sample whose `item_id` is not in the rows (raise); fewer Samples than rows (raise, no silent zip truncation); unparseable text (counted in errors, no crash); unknown view key (VerifyError `unknown_key`). Why: silent zip truncation shortens the eval set and inflates means.

53. **build_prompt_parity_with_trl_path** (I, data fixture, skip without). For 64 val rows (fixtures) with `qwen3-vl-2b` processor: `build_prompt(row, processor, model)` token ids equal TRL's path ids (`prepare_multimodal_messages` + `apply_chat_template(add_generation_prompt=True)`) as lists, exact. Property variant on synthetic rows: image sizes 1x1, 512x512, 4000x3, 0-length instruction text, text with `<|im_end|>` literal: both paths equal or both raise. Why: prompt drift between train and eval is the headline risk (G13 binding).

54. **archive_only_for_old_samplers** (U). `sample_vllm.py`, `hf_sample.py` exist only under `archive/`; `archive/` dirs each contain README.md; `archive` not on any `testpaths`, not in ruff `src`, not matched by `code_files()`. Why: G14.

---------------------------------------------------------------------------------------------------------------
## F. Verifier API (B9) including bit-equality and edge geometry

55. **verifier_bit_equal_to_score_batch** (I). Fixed corpus of 64 polygon texts (goldens). `Verifier("reward_v4", view_root).score(items)` vs `score_batch`: for each item `hx(reward)` equal, each term `hx` equal, defects equal. Run CPU only (GPU variant case 66). Also order: results in input order for shuffled input (shuffle 20 seeds; result[i].id == items[i].id). Why: B9, B10.3.

56. **verifier_golden_corpus_hex** (I). `score_batch` over the 64 texts equals `goldens/phase_e/invariance.darwin-arm64.json` per float `hx`; platform file missing -> skip with message. Corpus must include (check by test that the corpus contains these categories, so goldens cannot be vacuous): empty outline, single-pixel polygon (a 1x1 box in the 0-1000 grid), polygon touching all four borders (0,0)-(1000,1000), degenerate zero-area (collinear points), self-intersecting bow-tie, duplicate vertices, polygon with a hole, hole touching the outer ring, hole outside the outer ring, two nodules, nested outers, coordinates 1001 / -1 / float / string, truncated JSON, non-JSON, 1.5 MB text, 10,000 vertices, RIGA disc+cup with cup not inside disc. Expect the 13 defect classes are each hit at least once (assert union of defect keys across corpus covers all 13 names read from the checker). Why: vacuous goldens catch nothing.

57. **verifier_stamp_fields** (U). `stamp.api_version == API_VERSION == "1.0"`, `reward_name == "reward_v4"`, `reward_fingerprint == fingerprint()`, `checker_build in {"cuda","cpu"}` (on the Mac == "cpu"), `code_version` non-empty str; `Stamp` frozen (assignment raises `FrozenInstanceError`). `Verifier("reward_v9", root)` -> `VerifyError(code="spec_mismatch")` (or the documented code); `Verifier("reward_v4", missing_dir)` -> clear error. Why: stamps are the audit trail.

58. **verifier_unknown_key_whole_batch_fails** (U). Items: 3 valid keys + 1 unknown at position 2 -> `VerifyError` with `code == "unknown_key"` raised; NO partial results returned; message names the key. Key with the wrong VIEW_SEP form (`"a:b"`, `"a::b::c"`, empty key, key with newline) -> `unknown_key`, never an IndexError or path traversal: key `"../../etc/passwd"` and `"/abs/path"` must be `unknown_key` and must not open files outside the store (monkeypatch `open` / `Path.read_bytes` guard to fail on out-of-root paths). Why: store key is untrusted input (also HTTP).

59. **verifier_unparseable_is_result_not_exception** (U). 20 bad texts (empty, "null", "[]", "{", binary garbage, lone surrogate `"\ud800"`, 5 MB of `[`, NaN literal, deeply nested JSON depth 100,000 -> must not RecursionError) each give `ScoreResult(error=<non-empty str>, reward is None or 0.0 per pin)` and neighbours in the same batch still score normally (equal to scoring them alone). Why: one bad answer must not kill a 128-answer reward call.

60. **verifier_batch_independence_property** (P, hypothesis). For random batch of texts from the corpus, `score(batch)[i] == score([batch[i]])[0]` bitwise (`hx`) and invariant to permutation and to duplicates inside the batch; empty batch -> `[]`. Why: a batched kernel with cross-item state (shared buffer reuse, stale mask) would break this; this is the best adversarial test of the CUDA engine contract.

61. **verifier_score_mask_symmetries** (P). Using synthetic square nodule against a store with fixed reference masks: translating the polygon by (dx,dy) and the reference by the same moves nothing in IoU (bit-equal when integer shifts stay inside); IoU of a polygon vs itself == 1.0 exactly; swapping hole and outer yields a defect; reward in [0,1] (or the REWARD_V4 range) for 1,000 random boxes; monotonic: shrinking a box that fully covers the reference never increases IoU. Why: independent mathematical oracles, not recomputed from code.

62. **verifier_max_size_and_borders** (U). Outline covering the entire 512x512 image; polygon with 10,000 vertices; text of exactly the largest allowed length +/- 1 char (read limit from checker, `max_completion_length 1536` tokens equivalent); coordinate 0 and 1000 exactly (valid), 1001 (defect). Expect no crash; defect flags exactly those the checker documents. Why: border off-by-one in rasterising.

63. **verifier_thread_safety** (P). `ThreadPoolExecutor(8)` calling `score` 200 times concurrently with different batches; each result equals the single-thread result bitwise. Why: reward calls come from threaded callers (cf. existing test_reward_one_checker_call).

64. **http_parity_and_contract** (I; skip if fastapi/httpx absent, but CI coverage job must NOT skip: separate test `test_ci_installs_serve_extra` asserts `importlib.util.find_spec("fastapi")` when env `CI=true`). With `TestClient(app)`: `POST /v1/score` floats equal in-process (`hx` after `float(json)`; JSON roundtrip of Python floats is exact via repr; assert); every response (200, 4xx, 413, 422, 500) has `api_version`; error body shape exactly `{"error": {"code","message","retryable"}}` (no extra keys); `retryable` True only for `device_error` (parametrise over all four codes; force `device_error` by monkeypatching Verifier.score to raise `VerifyError("device_error")`). Batch of `MAX_BATCH` -> 200; `MAX_BATCH + 1` -> 413 `batch_too_large` and the verifier's score is NOT called (call counter). Malformed JSON body -> 422 or 400 with the same error shape. Missing `items`; items not list; item missing `truncated` (pin 422); `reward` unknown -> `spec_mismatch`. Body of 50 MB -> rejected without reading it all (use `Content-Length` check if implemented; otherwise record). `NaN` in JSON response never appears (reward None -> `null`). `GET /v1/spec`, `/livez`, `/readyz` -> 200 and `/readyz` -> 503 when the store is not loaded (monkeypatch). Unicode ids round-trip. Why: B9 HTTP surface.

65. **verify_api_has_no_torch_pydantic** (U). Subprocess `python -I -c "import objdel.verify.api,sys; assert not {'torch','pydantic','fastapi'} & set(sys.modules)"`. Why: B9 constraint.

66. **cuda_cpu_bit_equality** (G; skip without GPU). Same 64 corpus + 2,000 random polygons (seeded) scored by `checker_build="cuda"` vs `"cpu"` (same source, `gpucheck.cpp`): every reward, term, defect equal by `hx`; skip message when no GPU. Polygons include holes, touching borders, single pixel. Exactness is the contract of "one source for CPU and nvcc"; allow no tolerance unless the spec records one (open question to flag: float reductions on device may reorder; if bit-equality fails the test must print the max ulp difference). Why: the heart of the engine.

67. **device_error_is_retryable_failure_mode** (G/U). With a fake `_gpucheck` raising `RuntimeError("CUDA error")`, `Verifier.score` raises `VerifyError(code="device_error")` (not a bare RuntimeError); a second call after the failure still works when the fake recovers. Why: HTTP retryable semantics.

---------------------------------------------------------------------------------------------------------------
## G. polygonize and store (B13)

68. **polygonize_property_roundtrip** (P, hypothesis). Random binary masks (random rectangles unions, 8x8..64x64, with holes): `polygonize(mask)` then rasterise via `objdel.geom` rasteriser == mask exactly for simply-connected 4/8-connected shapes where the engine defines it; `polygonize(mask, offset=(dx,dy))` == `polygonize(mask)` with every vertex shifted by (dx,dy) (exact int equality); empty mask -> `[]`; single pixel -> one ring of 4 vertices (pin) with area 1; full mask touches all borders; mask with a hole -> one outer, one hole, hole orientation opposite; touching diagonal pixels -> pin the 4/8-connectivity choice; negative offset; non-bool dtype (uint8 0/255) either accepted equally or raises, never silently different; shape (0,0) -> `[]`. Production callers: AST scan shows no other function with polygonize behaviour is called (names `polygonize*`, `mask_to_poly*` other than the one) in src/experiments. Why: B13 single entry.

69. **store_loader_and_preload** (U). `loader(root1, root2)` resolves key `"ds/id::view"` from first root having it; unknown -> `KeyError`/`VerifyError` per layer; `VIEW_SEP == "::"`; ids containing `::` -> pin (split on last or first; adversarial case `a::b::c`). `preload(root, ["ds"], threads=0)` -> raises or uses 1 thread (pin), `threads=4` returns count == number of files; empty dataset -> 0; repeat preload idempotent; corrupt file in the store -> clean error naming the file. Why: threads=0 and corrupt files are failure modes.

---------------------------------------------------------------------------------------------------------------
## H. No stale strings, paths and root (B11, B12)

70. **host_syncs_not_increased** (U). For each module in `src/objdel` importing torch (AST): count `re.findall(r"\.(item|cpu|tolist|numpy)\(\)", source)` <= `host_syncs.json[module]` after mapping moved files; module absent from golden -> fails (new torch module must record a count; adversarial: new file with 1 sync and no golden entry must fail). Negative control test: feed the counter a source with `x.item()`, `y.cpu() `, `z.tolist()` and `w.numpy()` -> 4; `x.item(0)` and `x.items()` -> 0 (document what the regex misses: `.to("cpu")`, `.item(` with args; mention as known gap). Also strings/comments containing `.item()` count (conservative). Why: G11.

71. **objdel_strings_resolve** (U). Over string literals (via AST, so docstrings excluded or included? include all str constants; allowlist file for prose) in src, live experiments, tests matching `objdel(\.\w+)+`: resolve statically to a module file or `objdel._gpucheck`, and optional trailing top-level name defined by AST in it. Negative control: a fixture string `"objdel.train.fast.FastGRPO"` must fail. Strings like `"objdel.egg-info"` are covered by an explicit allowlist with a reason, and the allowlist must not contain an entry that now resolves (no dead allowlist). Why: B12, mock.patch strings go stale silently.

72. **test_path_strings_exist** (U). Regex `tests/[\w/]+\.py(::\w+)?` over pyproject, .github, live experiments, specs, src, Makefile: file exists; `::name` is a top-level `def`/class in the file (AST). Include `tests/test_gpu_scoring.py` in `code_files()`. Negative control with a temp file. Why: CI gates pointing at a renamed test pass vacuously.

73. **repo_root_and_secrets** (U). `OBJDEL_REPO=tmp` with a `pyproject.toml` -> `repo_root() == tmp`; `OBJDEL_REPO` pointing at a dir without pyproject -> raises (not silent); unset -> git toplevel (monkeypatch `subprocess.run`) ; not a git repo and no env -> clear error; calling from `chdir(other)` gives the same root. `OBJDEL_SECRETS` unset -> `secrets.path() == repo_root()/".secrets"`; set -> that path. AST: no file in `src/objdel/infra/` uses `__file__` to derive a root (`Path(__file__).parents`). Why: B12 relocatable.

74. **code_files_snapshot** (U). `code_files()` excludes any path with parts `.claude` or `archive` (create fixture tree with `.claude/worktrees/x.py`, `archive/old/y.py`, `src/objdel/z.py`, `tests/test_gpu_scoring.py`, `.secrets`, `data/x.npy`, `__pycache__/a.pyc`, `.venv/..`); result contains the last-but-one fixture; contains NO `.secrets`, NO `data/`, NO `.pyc`; symlink escaping the repo root excluded. Why: secrets must never ship to a rented box.

---------------------------------------------------------------------------------------------------------------
## I. Regression of today's numbers and performance gates

75. **reward_fingerprint_golden** (U). golden 1: `sha256(repr((astuple(REWARD_V4), [min_area(a) for a in (0,999,1001,10**4,10**6)])))` equals literal; full `fingerprint()` hex equals literal only if `numpy`, `cv2`, `scipy` versions equal the recorded ones else skipped with message that names the versions. Property: `fingerprint()` changes when one weight is perturbed by 1 ulp (monkeypatch). Why: reward identity.

76. **target_set_hash_golden** (U). golden 2 for fixture items; skip without fixtures. Why: targets unchanged by the move.

77. **tests_count_gate** (E2E, run by CI script not pytest). `pytest -q -m "not slow"` summary: failed 0; passed >= 1447 - removed (list read from commit file `tests/goldens/phase_e/removed_tests.txt`); skipped <= 6 + allowed new skips; `deselected` = 31 unless commit says otherwise. Adversarial: the number of removed tests listed must equal the diff of collected node ids (`pytest --collect-only -q` against `base_ids.txt`), so deleting a test without listing fails. Why: G1 cannot be gamed by silent deletion.

78. **perf_check_private_names** (I). 5 cold subprocess runs, median < 10 s, max < 15 s (Mac). Also 100 warm repeated calls in one process: total < 1 s (source hashing cached, no re-reading TRL each call). Why: G5, boot latency.

79. **perf_lint_imports** (I). `lint-imports` wall < 60 s on the Mac and `make lint` < 120 s; fails on regression > 2x of recorded `perf.json`. Why: CI cost; also guards against caching hiding a failure (run with `--no-cache`).

80. **perf_verifier_cpu_throughput** (I). Scoring 128 answers (the GRPO batch) on CPU build: median of 5 runs < recorded bound x 1.25 (`perf.json`, per platform; skip without file). Equal-size batch via the API vs raw `score_batch`: API overhead < 5 % + 2 ms (the API is a thin wrapper). Memory: peak RSS growth over 50 repeated calls < 50 MB (no leak in the store cache). Why: Phase E must not slow rewards (step time 23.6 s is dominated by training, so a 2x reward slowdown matters).

81. **perf_http_overhead** (I). `TestClient` 128-item request: server-side time minus in-process time < 50 ms; 413 path answers in < 100 ms for MAX_BATCH+1 items of 1 KB each. Why: DoS guard.

82. **perf_gpu_scoring_gate** (G). Existing `tests/test_gpu_scoring.py` numbers unchanged: the GPU scorer scores 128 answers no slower than recorded x 1.25 (skip without GPU). Also no mask copy to host: patch `torch.Tensor.cpu`/`.numpy` counters during `Verifier.score` on GPU: 0 calls with mask-sized tensors (only score arrays). Why: memory rule "pixels never return to CPU".

83. **import_time_budget** (I). `python -X importtime -I -c "import objdel.core, objdel.geom, objdel.verify"` total < 2 s (Mac) and no module > 0.5 s self time that is not numpy/cv2/scipy. Why: ML stack leakage shows up as seconds.

---------------------------------------------------------------------------------------------------------------
## J. End-to-end

84. **e2e_tiny_rl_run_from_spec** (E2E, slow-marked). `load_spec` of a tiny spec (fixture YAML with steps 2, save_steps 1, tiny model, CPU) -> `run_rl` -> metrics file valid JSON (no NaN), checkpoint dirs at steps 1 and 2, `run_spec.json` present, rewards equal the golden floats `hx` for step 1 (seed 0). Resume: `resume=True` after step 1 continues to step 2 with identical step-2 loss (`hx`) as the uninterrupted run. Kill during step (SIGTERM via fake callback raising) -> no half-written checkpoint directory (checkpoint written to tmp then renamed; assert no dir without `adapter_model.safetensors`). Why: ties B6, B8, B4 and persistence together.

85. **e2e_eval_path_vs_reward_path_agree** (E2E). Same fake rollouts through `evaluate()` and through the reward function used by trainer (`score_batch` via the rl reward callback): mean reward equal `hx`. Why: G13 claim that evaluation and trainer use the same scoring.

---------------------------------------------------------------------------------------------------------------
## Pytest sketches for the 5 most important tests

Selection reasons: (1) bit-equality of the verifier (core invariant, B9/B10.3), (2) batch independence property (catches kernel state leaks), (3) check_private_names negative cases (G5 gate), (4) YAML strictness incl. duplicate keys (money-relevant spec safety), (5) lint-imports planted violation (guards the guard).

```python
# tests/test_phase_e_adversarial.py
import copy, json, shutil, subprocess, sys, inspect
from pathlib import Path
import pytest
from hypothesis import given, settings, strategies as st

ROOT = Path(__file__).resolve().parents[1]
GOLD = ROOT / "tests" / "goldens" / "phase_e"


def hx(x):  # None stays None, floats become exact hex
    return None if x is None else float(x).hex()


# ---- 1. Verifier bit-equal to score_batch, on the recorded corpus (case 55/56) ---------------------------------
@pytest.fixture(scope="module")
def corpus():
    return json.loads((GOLD / "corpus64.json").read_text())      # [{"id","key","text","truncated","job"}...]


def test_verifier_bit_equal_to_score_batch(corpus, view_root):
    from objdel.verify.api import Verifier, ScoreItem
    from objdel.verify.reward import score_batch

    items = [ScoreItem(c["id"], c["key"], c["text"], c["truncated"]) for c in corpus]
    got = Verifier("reward_v4", view_root).score(items)
    ref = score_batch([tuple(c["job"]) for c in corpus])          # independent entry point
    assert [r.id for r in got] == [c["id"] for c in corpus]
    for i, r in enumerate(got):
        assert hx(r.reward) == hx(ref["reward"][i]), (i, corpus[i]["text"][:80])
        for k, v in r.terms.items():
            assert hx(v) == hx(ref[k][i]), (i, k)
        assert dict(r.defects) == ref["defects"][i]
    plat = f"{sys.platform}-{__import__('platform').machine()}"
    gold = GOLD / f"invariance.{plat}.json"
    if not gold.exists():
        pytest.skip(f"no golden for {plat}")
    assert [hx(r.reward) for r in got] == json.loads(gold.read_text())["score_batch_rewards"]


# ---- 2. Batch independence / permutation invariance (case 60) ---------------------------------------------------
@settings(max_examples=50, deadline=None)
@given(st.data())
def test_score_is_batch_independent(data, corpus, view_root):
    from objdel.verify.api import Verifier, ScoreItem

    v = Verifier("reward_v4", view_root)
    idx = data.draw(st.lists(st.integers(0, len(corpus) - 1), min_size=1, max_size=40))   # duplicates allowed
    items = [ScoreItem(f"{n}", corpus[i]["key"], corpus[i]["text"], corpus[i]["truncated"]) for n, i in enumerate(idx)]
    together = v.score(items)
    for it, r in zip(items, together):
        alone = v.score([it])[0]
        assert hx(r.reward) == hx(alone.reward) and dict(r.terms) == dict(alone.terms)
    perm = data.draw(st.permutations(range(len(items))))
    shuffled = v.score([items[p] for p in perm])
    assert [hx(shuffled[j].reward) for j in range(len(perm))] == [hx(together[p].reward) for p in perm]


# ---- 3. check_private_names negative cases (case 16) ------------------------------------------------------------
SIX = ["_get_per_token_logps_and_entropies", "_generate_and_score_completions", "_prepare_inputs",
       "_compute_loss", "_generate", "_calculate_rewards"]   # the real six come from the adapter's recorded table


def test_check_passes_on_pin():
    from objdel.rl import trl_adapter
    assert trl_adapter.check_private_names() is None


@pytest.mark.parametrize("name", SIX)
@pytest.mark.parametrize("mode", ["missing", "signature", "source"])
def test_check_detects_drift(monkeypatch, name, mode):
    from trl import GRPOTrainer
    from objdel.rl import trl_adapter as A

    orig = getattr(GRPOTrainer, name)
    if mode == "missing":
        monkeypatch.delattr(GRPOTrainer, name)
    elif mode == "signature":
        ns = {}
        exec("def f(self, renamed_param): pass", ns)
        monkeypatch.setattr(GRPOTrainer, name, ns["f"])
    else:  # same signature, different body -> SHA-256 differs
        ns = {}
        params = ", ".join(inspect.signature(orig).parameters)
        exec(f"def f({params}):\n    return None", ns)
        monkeypatch.setattr(GRPOTrainer, name, ns["f"])
    with pytest.raises(A.TrlDriftError) as e:
        A.check_private_names()
    assert name in str(e.value) and "1.17.1" in str(e.value)


@pytest.mark.parametrize("ver", ["1.17.2", "1.17.1.dev0", "", None])
def test_version_drift(monkeypatch, ver):
    import trl
    from objdel.rl import trl_adapter as A
    monkeypatch.setattr(trl, "__version__", ver)
    with pytest.raises(A.TrlDriftError, match="1.17.1"):
        A.check_private_names()


# ---- 4. Spec strictness incl. YAML duplicate keys and unknown keys at every level (cases 34, 35) -----------------
SPECS = sorted((ROOT / "specs").glob("*.yaml"))


def _paths(node, prefix=()):
    if isinstance(node, dict):
        yield prefix
        for k, v in node.items():
            if k != "grpo":                      # free-form mapping by design
                yield from _paths(v, prefix + (k,))
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from _paths(v, prefix + (i,))


@pytest.mark.parametrize("spec", SPECS, ids=lambda p: p.name)
def test_unknown_key_rejected_at_every_level(spec, tmp_path):
    import yaml, pydantic
    from objdel.core.spec import load_spec
    base = yaml.safe_load(spec.read_text())
    load_spec(spec)                                                    # control: valid
    for path in _paths(base):
        doc = copy.deepcopy(base)
        node = doc
        for p in path:
            node = node[p]
        node["zz_unknown"] = 1
        f = tmp_path / "s.yaml"
        f.write_text(yaml.safe_dump(doc))
        with pytest.raises(pydantic.ValidationError):
            load_spec(f)


def test_duplicate_yaml_key_is_an_error(tmp_path):
    from objdel.core.spec import load_spec
    text = (ROOT / "specs" / "rl_shifted.yaml").read_text()
    f = tmp_path / "dup.yaml"
    f.write_text(text + "\nschema_version: 1\n")                       # second top-level key of the same name
    with pytest.raises(Exception) as e:                                # type pinned by the implementation
        load_spec(f)
    assert "schema_version" in str(e.value) or "duplicate" in str(e.value).lower()


# ---- 5. lint-imports guards the guard: planted violations must fail (case 2) -------------------------------------
@pytest.mark.parametrize("target, line", [
    ("core/__init__.py", "import objdel.rl"),
    ("verify/__init__.py", "import objdel.data"),
    ("sft/__init__.py", "import objdel.rl"),
])
def test_planted_layer_violation_fails(tmp_path, target, line):
    dst = tmp_path / "repo"
    shutil.copytree(ROOT / "src", dst / "src", ignore=shutil.ignore_patterns("__pycache__", "*.so"))
    shutil.copy(ROOT / "pyproject.toml", dst)
    for extra in (".importlinter",):
        if (ROOT / extra).exists():
            shutil.copy(ROOT / extra, dst)
    with open(dst / "src" / "objdel" / target, "a") as fh:
        fh.write(f"\n{line}\n")
    r = subprocess.run(["lint-imports", "--no-cache"], cwd=dst, capture_output=True, text=True, timeout=120)
    assert r.returncode != 0, r.stdout
    assert "BROKEN" in r.stdout.upper() and "objdel" in r.stdout
    clean = subprocess.run(["lint-imports", "--no-cache"], cwd=ROOT, capture_output=True, text=True, timeout=120)
    assert clean.returncode == 0, clean.stdout
```

## Open questions the tests expose (spec is silent; tests pin the implementation's answer)
- RlStage default `save_steps=50` with `steps < 50` raises (spec literal); is a short probe-length RL stage meant to need an explicit save_steps?
- Is `stages: []` valid? Is the YAML loader strict on duplicate keys (PyYAML default is last-wins)?
- Case 66 demands CPU/GPU bit-equality; the spec only says "exactness CPU vs GPU where relevant". If device reductions reorder floats, a tolerance must be written into the spec.
- `check_private_names` behaviour when `inspect.getsource` fails (installed pyc-only).
- Case-sensitivity of REGION ("us").
