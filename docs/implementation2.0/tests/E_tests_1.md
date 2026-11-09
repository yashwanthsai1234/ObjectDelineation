# Phase E test cases, set 1 (written from E.spec.md only)

## Skills used
- engineering:testing-strategy: test plan shape (what, level, example cases, gaps).
- mattpocock-skills:tdd: tests at public seams, independent expected values (literals/goldens, never recomputed from the code under test), vertical slices.
- nvidia-skills:nvidia-skill-finder: loaded; Phase E needs no GPU, so no NVIDIA skill applies beyond the GPU-marked skip rule.
- Other listed skills (code-review, debug, system-design, architecture, documentation, huggingface-best): not loaded; they do not change test design for a pure restructure.

Conventions: style of REPO/tests (plain pytest, hypothesis, `tmp_path`, `monkeypatch`, autouse secrets fixture, `pytest.skip` when data is missing, `-I` for data). Spec refs are B1..B14 and G1..G16. "Level" is unit/property/integration/GPU/end-to-end. Test files named in the sketches are proposals under `tests/phase_e/`.

---
## A. Layers, imports, direction (B1, B2, B3, B12, B14)

**1. layers_contract_passes** (integration). Setup: run `lint-imports` via subprocess in REPO. Input: repo config. Expect: exit 0; config has a `layers` contract over container `objdel`, layers exactly `infra > {rl | sft} > policy > data > verify > geom > core` (rl and sft as independent siblings), `exhaustive = true`, exhaustive ignore only `_gpucheck`, `ignore_imports` empty, every contract has `broken_contract_guidance`. Parse `pyproject.toml`/`.importlinter` with tomllib/configparser and assert each. Why: a loosened config (extra ignore, missing guidance, non-exhaustive) silently disables the architecture.

**2. layers_contract_catches_violation** (integration, adversarial). Setup: copy `src/objdel` to tmp, append `import objdel.rl` to `objdel/core/x.py`, run lint-imports with that package path. Expect: non-zero exit, message names `core -> rl`. Repeat for `sft -> rl`, `rl -> sft`, `verify -> policy`, `geom -> data`. Why: proves the contract is enforcing, not vacuous (e.g. wrong container name passes forever).

**3. layers_unlisted_package_fails_exhaustive** (integration). Add `objdel/newpkg/__init__.py` in the tmp copy. Expect lint-imports fails (exhaustive). Why: B1 exhaustive.

**4. old_packages_gone** (unit). For each of `objdel.checker, geometry, reward, eval, targets, train, train.fast` assert `importlib.import_module` raises `ModuleNotFoundError`; assert directories absent under `src/objdel`; no `objdel.secrets` flat module if moved (only if spec lists it; skip otherwise). Assert no file named `fast.py`. Why: a leftover shim or stale `__pycache__`-only package (namespace-package ghost) would keep old imports alive. Also check `importlib.util.find_spec` is None (catches namespace package from a stray empty directory).

**5. new_packages_exist** (unit). Import each of `core geom verify data policy rl sft infra` and `objdel._gpucheck`. Expect success on Mac CPU.

**6. only_geom_engine_imports_gpucheck** (unit, AST). Walk `src/`, `tests/`, `experiments/` (exclude `archive/`); parse AST; collect `import objdel._gpucheck`, `from objdel import _gpucheck`, `from objdel._gpucheck import ...`, `importlib.import_module("objdel._gpucheck")`, string `"objdel._gpucheck"` used in `__import__`. Expect: in `src/` only `objdel/geom/engine.py`. Also: every name accessed as `_gpucheck.<n>` / imported from `_gpucheck` anywhere is in `geom.engine.__all__`; every name in `__all__` is an attribute of the real extension and `geom.engine.<n> is _gpucheck.<n>`. Why: `from objdel import _gpucheck` form is easy to miss with a regex; `__all__` drift breaks callers at import.

**7. no_ml_stack_below_policy** (integration, subprocess). For each pkg in `core geom verify data` and for every submodule discovered with `pkgutil.walk_packages`: `python -I -c "import <mod>, sys; bad={'torch','trl','transformers','vllm','peft','datasets','accelerate'}; print(sorted(bad & {m.split('.')[0] for m in sys.modules}))"`. Expect `[]`. Fresh process per package (not per session; sys.modules pollution hides leaks). Also `import objdel.verify.api` alone and `import objdel.policy.sampler`, `objdel.policy.evaluate` must not load `trl.trainer` or `vllm` (B7). Why: one top-level `import torch` in a helper breaks CPU-only verifier boxes.

**8. dependency_direction_ast** (unit). AST over `src/objdel`: no import of `experiments`, `tests`, `oracles`, or module whose top part starts with `test_`; no string literal passed to import_module with those; no `sys.path` mutation (`sys.path.insert/append/extend`, `sys.path =`, `site.addsitedir`). Plus tests outside `tests/experiments/` import no `experiments` (AST, including `from experiments.x import y`). Negative self-test: feed the checker a synthetic source with each violation and expect it flagged (guards the checker itself).

**9. import_check_string** (unit). `from objdel.infra.vast import IMPORT_CHECK`; parse as Python; assert imported modules include `pydantic`, `vllm`, `objdel.policy.sampler`; none start with `tests` or `experiments`. Run it in a subprocess on Mac: expect failure only on `vllm` ImportError (not a tests/experiments error) when vllm is absent.

**10. live_cli_help_loads_no_test_modules** (integration). For every script under `experiments/**` outside `archive/`: `python -I -m <mod> --help` with `-X importtime` or a wrapper `-c "runpy; print modules"`; expect exit 0 and no `sys.modules` key starting `tests.`/`oracles.`. Why: a CLI importing a test helper breaks on boxes that ship no tests.

**11. harness_locations** (unit). `tests/oracles/` contains the CUDA parity harness; `tests/harness/` the K3 harness; `src/` contains no `score_numpy` or `oracles`.

**12. stale_module_strings** (unit, B12). Regex `objdel(\.\w+)+` over string literals (AST Constant str, plus docstrings excluded or included consistently) in src, live experiments, tests; resolve statically: longest prefix must be a module file/package dir; remainder at most one name, which must be a top-level def/class/assign in that file (AST). Expect none unresolved. Negative self-test: `"objdel.train.fast"` and `"objdel.geom.nope_fn"` flagged. Edge: `objdel._gpucheck` accepted; attribute-chained text like `objdel.geom.engine.__all__` has name `__all__` (assign, ok); file-name-like `objdel.egg-info` not matched by `\w+`.

**13. stale_test_paths** (unit). Regex `tests/[\w/]+\.py(::\w+)?` over pyproject.toml, `.github/**`, live experiments, specs, src, Makefile. Each file exists; each `::name` is a top-level def in that file (AST). Must include `tests/test_gpu_scoring.py`. Edge: parametrised ids `::test_x[abc]` are not matched by regex (acceptable) but a rename of the function is caught.

**14. every_name_has_a_home** (unit, B14). Golden `tests/goldens/phase_e/names.json`: for each of `train/dataset, configs, fast, rewards, callbacks, serve` the top-level name set recorded on the base commit, the new module list, and the archived/renamed map. Expect `old == union(new modules' top-level names) - archived - renamed_old + renamed_new`. Computed by AST (no import). Negative: delete one function from a moved module in tmp and expect failure. Edge: names from `from x import y` re-exports do not count as definitions; `__all__` aside.

---
## B. TRL adapter (B4, G5, G6)

**15. trl_imports_only_in_adapter** (unit, AST). `trl.trainer.grpo_trainer` and `trl.generation.vllm_generation` imported (import, from-import, string import) only in `src/objdel/rl/trl_adapter.py` (and test files of the adapter). Why: spreading private imports defeats drift detection.

**16. private_names_only_in_adapter** (unit, AST). For the 11 private names, find attribute access `self.<n>`, `trainer.<n>`, `super().<n>`, `<trl module alias>.<n>` in src outside the adapter, but only inside classes deriving GRPOTrainer/VLLMGeneration or on those variables. Expect none. Edge: `self._metrics` in an unrelated class (e.g. a Logger) must NOT be flagged (spec says other classes may). The test needs its own class-context logic; add a self-test fixture with a non-TRL class using `self._step`.

**17. no_assignment_to_shuffle_sequence_dict** (unit, AST + grep). No `Assign/AugAssign/setattr/monkeypatch.setattr` targeting `shuffle_sequence_dict` on a TRL module anywhere in src; also assert at runtime: after building a FastGRPO on the tiny model and running a step, `trl.trainer.grpo_trainer.shuffle_sequence_dict is ORIGINAL` captured before. Why: the old code swapped the global during `_prepare_inputs` and restored it; an exception between swap and restore would leak. Adversarial runtime test: make FastGRPO's `_prepare_inputs` raise midway (monkeypatch a callee) and assert the global is still the original.

**18. check_private_names_passes** (integration). On pinned trl 1.17.1: returns `None`. Time < 10 s in a fresh subprocess including `import trl` (G5). Skip with a message if installed trl != 1.17.1 (then see test 19a instead).

**19. check_private_names_negative_cases** (unit, monkeypatch). Each must raise `TrlDriftError` whose message contains the offending name and `1.17.1`:
 a. `trl.__version__ = "1.18.0"` (also "1.17.0", "1.17.1.dev0", "1.17.10": exact match only).
 b. delete one of the six overridden methods from `GRPOTrainer` (parametrise all six).
 c. replace a method with a function of different parameter names (rename `inputs` to `batch`; also add a parameter; also remove one).
 d. replace a method with a function of the same signature but different source (SHA-256 mismatch), e.g. add a `pass` line. Also a comment-only source change (must raise if SHA is over the raw source; document the decision: expected behaviour = raise).
 e. built trainer lacking each of `_metrics`, `_step`, `_buffered_inputs`, `num_iterations` (parametrise four; delete attr from a built tiny trainer or use a stub object per the adapter's check hook).
 Why: a drift guard that never fires is worse than none. Also assert the passing case still passes after each monkeypatch is undone (no cross-test leak).

**20. fastgrpo_checks_before_super_init** (unit). Monkeypatch `check_private_names` to raise a sentinel and `GRPOTrainer.__init__` to record a call; construct FastGRPO; expect sentinel raised and the recorded init NOT called. Why: ordering requirement of B4; a call after TRL init wastes a model load and fails later with a confusing error.

**21. env_check_records_trl_private_names** (integration). Run `env_check.sh` in a mode that writes `env_check.json` (or its python part) into tmp; expect key `trl_private_names` present with passed status. With check_private_names forced to raise: key present with failed status and the message, and script exit non-zero. Skip if script needs box-only tools; then test the python helper it calls.

**22. fast_microbatches_and_loss_golden** (integration, CPU). On the tiny offline model (tests/conftest `tiny`), FastGRPO with dr_grpo, shared_prefix in {False, True}, seed 0: micro-batch composition (list of index lists / lengths) and loss (`float.hex`) equal golden recorded on pre-change commit. Why: this is the core move-invariance claim for the removal of the global patch.

**23. fast_fallback_equals_trl** (integration, parametrised). Cases: eval mode (`model.eval()`/ not training); loss_type != dr_grpo; entropy cut on (`top_entropy_quantile < 1`); `gradient_accumulation_steps % steps_per_generation != 0`. For each: `FastGRPO._prepare_inputs(batch)` output equals `GRPOTrainer._prepare_inputs` run with the same RNG state: tensor-by-tensor `torch.equal`, same list structure/order. Plus the positive control: with none of the conditions, output differs in ORDER from TRL's (length ordering is active) but is a permutation of the same rows (multiset equality). Why: proves the fallback condition is exactly as specified, and that the positive path is not silently the fallback.

**24. length_order_property** (property, hypothesis). Random lengths (1..1536), group size 8, random micro-batch counts: ordering function output is a permutation of input indices; whole groups never split when shared-prefix is on; micro-batch max length sum <= that of random order (never worse); deterministic given seed; None and 0-d tensors in the sequence dict pass through untouched (spec carries over the field rules of TRL's shuffle). Edge: one prompt, all equal lengths, length 1, 1536 (max).

---
## C. Single sources (B5)

**25. rl_lr_values** (unit). `rl_lr(m)` is `5e-5` for every LoRA model key, `5e-6` for full fine-tune, exact `==`; iterate over all model keys in `objdel.policy.models`. Unknown key raises (KeyError/ValueError; whichever, not a silent default). Expect `rl_lr(m) == G_PROBE["lr_lora" if LoRA else "lr_full"]` for every key. `G_PROBE` keeps `lr_full == 5e-6`, `lr_lora == 5e-5`.

**26. lr_default_flows_into_config_and_trainer** (integration). `grpo_config(m, out)` has `learning_rate == rl_lr(m)`; with `learning_rate=1e-4` override it is 1e-4; `build_trainer(..., flags=TrainerFlags(lr=None))` returns `Built.lr == rl_lr(m)`; with `lr=2e-5` returns 2e-5 (not 0 or falsy-replaced: also `lr=0.0` must be kept literally or rejected, never replaced by default; `is None` check not truthiness). Why: classic `lr or default` bug.

**27. no_other_reader_of_G_PROBE_lr** (unit, AST). `G_PROBE["lr_full"]`/`["lr_lora"]` subscripts or `.get("lr_*")` appear only in `prereg.py` (and the equality test).

**28. weights_single_source** (unit). `verify.stats.WEIGHTS == dict(zip(reward.TERMS, REWARD_V4.weights))`, key order identical. AST search: no dict literal with the same term names and numeric weights elsewhere in src (heuristic: any dict whose keys are a superset of TERMS with float values, outside reward.py).

**29. region_constant** (unit). `objdel.infra.vast.offers.REGION == ("US", "CA")` (tuple, not list). No `regions` field in `ComputeSpec.model_fields`, in `RunSpec`, or in any spec YAML. Offer filtering: offers with country codes "US", "CA", "DE", "CN", "" and lowercase "us": only US and CA pass (decide: lowercase is excluded unless code upper-cases; test documents actual spec: exact tuple membership). Why: policy constant memory rule "North America only".

**30. no_worker_pool_in_rl_training** (unit/integration). AST: no `fork_pool`, `ProcessPoolExecutor`, `multiprocessing.Pool` in experiments RL training scripts. Runtime: run `run_rl` on the tiny model with `fork_pool` monkeypatched to raise; passes.

---
## D. Typed run spec (B6, B8)

**31. load_spec_every_yaml** (integration). Every `specs/*.yaml` loads to a frozen `RunSpec`; assigning a field raises `ValidationError` (frozen); nested stages frozen; `hash(stage)` is not relied on (spec says never hashed: assert no `__hash__` use in loader).

**32. unknown_keys_rejected_at_every_level** (property/unit). Take a valid spec dict; inject `extra="x"` at: top level, `compute` (incl. `regions`), each stage kind, `grpo` mapping stays free-form (JsonValue; extra keys inside `grpo` ARE allowed: assert this to avoid over-rejecting). Expect `pydantic.ValidationError` for all but `grpo`. Parametrise by walking the model tree so a new model without `extra="forbid"` fails.

**33. rejections** (unit, parametrised). Each raises ValidationError: unknown `kind: "foo"`; missing `kind`; `schema_version: 2`, `0`, `"1"` (string; Literal[1] strict? expect raise in strict, document), `1.0` (document); duplicate stage names (two stages named "a", also case variations are NOT duplicates); `steps: 0`, `-5`; `save_steps: 0`; `save_steps > steps` (51 vs 50); `cap_hours: 0`, `-1`, `nan`; `ProbeStage` with both `ckpt` and `select`, with neither. Passing boundary cases: `save_steps == steps` ok; omitted `save_steps` default 50 with `steps=10` -> expect raise (50 > 10)! That is a real trap: default 50 and `steps < 50`. Test pins the behaviour: either spec accepts (validator applies only to explicit value) or raises; whichever the code does must be deliberate and documented. Flag as open question if the spec's "save_steps > steps raises" is ambiguous for the default. Also `steps=50` default OK.

**34. sft_stage_fields_equal_run_args** (unit). `set(SftStage.model_fields) - {"kind","name"} == set(dataclasses.fields(train_sft.RunArgs))` names; defaults equal; for `specs/r1_qwen3vl2b.yaml` the SftStage field values equal `RunArgs` built from the old YAML (golden).

**35. to_argv_goldens** (integration). `to_argv(stage, rd, data)` for every stage of `rl_shifted`, `rl_full`, `r5_shared` equals golden 5 lists exactly (string equality, order, no extra whitespace); for `r1_qwen3vl2b`, select/probe/cmd match; SFT by field values. Compare after the recorded module rename map only. Parametrise rd/data with paths containing spaces, `{`, and unicode to confirm substitution is literal and single-pass (e.g. `data="/a/{rd}"` must not re-expand).

**36. cmd_stage_literal** (unit). `CmdStage(argv=("echo", "{rd}/x", "{data}", "{other}", "$HOME", "*"))` -> `["echo", "<rd>/x", "<data>", "{other}", "$HOME", "*"]`. Unknown braces untouched; no shell interpretation; empty `argv=()` raises or yields empty (document; prefer raising since exec of nothing is a bug: flag).

**37. flags_to_argv_property** (property). For random `TrainerFlags` combos, argv contains `--shared-prefix` iff flag set, `--plain` iff plain, `--lr X` iff lr not None (lr `0.0` present), `--compile mode` iff set, `--trace` iff set; `grpo` mapping serialised deterministically (sorted keys) and round-trips; False flags absent.

**38. cli_old_flags_equal_spec** (integration). For each live CLI: parse old flags -> stage object `A`; parse `--spec specs/X.yaml --stage name` -> stage `B`; `A == B` for flags that mirror the YAML. Missing `--stage` with `--spec` -> argparse error exit 2; unknown stage name -> error naming available stages; `--spec` plus conflicting old flag -> define behaviour (error preferred) and test it.

**39. build_trainer_signature** (unit). `inspect.signature(build_trainer)` has typed params (model, ckpt, rows, out, roots, steps, flags); passing an `argparse.Namespace` as `flags` raises TypeError/ValidationError (not attribute fallback). No function named `make_trainer` in src/experiments/tests (AST).

**40. build_trainer_plain_flag** (integration, CPU tiny). `flags.plain=True` -> `type(trainer) is GRPOTrainer`; False -> `FastGRPO`. `Built.lr` consistent with config.

**41. grpo_config_golden** (integration). `grpo_config` dict (`to_dict()` minus path/time fields: output_dir, logging_dir, run_name?, any key containing "time", "dir", "path") for `train_rl` stage of `rl_shifted` and for a default G-probe on `qwen3-vl-2b` equals golden 6 via canonical JSON compare. Failure message prints the diff of keys. Why: silent hyperparameter drift (e.g. a flag default moved) changes training with no failing unit test.

**42. preflight_with_spec** (integration). Preflight passes on `specs/rl_shifted.yaml` with fixture data. Then append a newline/space change to a tmp copy: launch refusal fires; revert bytes: passes again (stamp keyed on path AND sha256). Same content at a different path: refused (path in key). Missing spec file: clear error. Edge: spec edited between preflight and launch (race) is refused since the stamp is recomputed at launch.

**43. run_spec_json_written** (integration). `run_rl` on tiny model writes `runs/<id>/run_spec.json` with keys: spec (round-trips to equal `RunSpec`), git sha (40 hex or "unknown" in non-git tmp), trl/torch/transformers versions equal installed, reward fingerprint equal `fingerprint()`. Written before training begins (kill after first step still leaves the file). No NaN in json (strict loads).

**44. save_json_contract** (unit/property). `save_json` rejects `nan`, `inf`, `-inf` (ValueError) including nested and inside lists; serialises `np.float32`, `np.int64`, `np.bool_`, `np.ndarray` (nested, 0-d, empty) to plain JSON; result readable with `json.loads(parse_constant=raise)`. Atomic: if serialisation fails halfway, no partial file (target absent or unchanged). Edge: numpy scalar NaN inside array.

---
## E. Sampler and evaluation (B7)

**45. params_for_exact** (unit). For every model key: `params_for(m, n=8, seed=3, stop_token_ids=(1,2)) == sampling(m, "vllm") | {"n":8,"seed":3,"stop_token_ids":[1,2]}`; includes `min_p`, `repetition_penalty`, `presence_penalty`, `frequency_penalty` keys (each asserted present); `stop_token_ids` is a `list` not tuple; empty stop ids -> `[]`; `n=0` and negative raise (or documented). Input mapping not mutated (deepcopy compare).

**46. request_seed** (property). Deterministic for same `(seed,item_id)`; differs across different item ids over 10k random ids with collision rate ~0 (< 1 in 10^4 for 31-bit space); within `[0, 2**31)` or vLLM's valid range; independent of Python hash randomisation (run in two subprocesses with different `PYTHONHASHSEED`, equal values: catches use of builtin `hash()`). Golden literal for 3 inputs, recorded on base commit.

**47. requests_shape** (unit). `requests(...)` per row: seed equals `request_seed`; `mm_processor_kwargs` present iff line has processor kwargs (both branches); prompt text from `build_prompt`; multimodal data list non-empty for image rows. Row ordering preserved; empty rows list -> empty list. Duplicate item ids: define (raise preferred).

**48. env_set_before_vllm_import** (integration, subprocess). `python -I -c "import objdel.policy.sampler, os; print(os.environ.get('VLLM_USE_FLASHINFER_SAMPLER'))"` prints `0` and `'vllm' not in sys.modules`. For `VllmSampler.__init__` with stubbed `vllm` module: env set before the stub's import statement executes (stub records env at import time). If env is pre-set to "1" by user: behaviour pinned (override to "0"; spec says set).

**49. build_prompt_equals_trl_ids** (integration; data+tokenizer required, skip otherwise). 64 validation rows of `qwen3-vl-2b`: `build_prompt` token ids == TRL path (`prepare_multimodal_messages` + `apply_chat_template(add_generation_prompt=True)`) token-for-token. Binding at the box (G13); offline tiny-tokenizer variant on CPU with synthetic rows to run everywhere. Edge rows: no image, very long text, unicode, empty system prompt.

**50. hf_sampler_contract** (integration, CPU tiny llama). `HfSampler.sample(rows, params)` returns one `Sample` per row in input order, `len(texts)==len(token_ids)==len(truncated)==n` (n=1, 4, 8); `truncated[i] == truncated(token_ids[i], stop)`. Empty rows -> `[]`. Greedy/seeded determinism: same seed same ids twice; different `request_seed` for different items.

**51. truncated_fn** (property+unit). `truncated([], stops)` False (empty completion: pin); last id in stops -> False; no stop id anywhere -> True; stop id present in middle but not at the end -> define (hypothesis check vs an independent reference: "truncated iff no stop id is the final token" or "iff none in ids", whichever the spec of vLLM finish_reason uses; test must be written against vLLM semantics: finish_reason=="length" ⇒ last token not stop). Empty `stop_token_ids` -> True for any non-empty ids. Max length 1536 ids ok.

**52. check_memory** (unit). `check_memory(lambda: (free, total), need)`: free/total ratio below `need` raises with numbers in message; exactly at `need` passes (boundary, `>=`); `total == 0` raises, not ZeroDivisionError; `mem_get_info` raising propagates.

**53. sampler_main_child** (integration). `python -I -m objdel.policy.sampler --help` exit 0; `--kind` outside vllm|hf -> exit 2; `select_ckpt.sampler_cmd(...)` returns an argv whose first elems are `[sys.executable or python, "-I", "-m", "objdel.policy.sampler", "--kind", ...]` and that argparse in `main` accepts (feed the argv, minus interpreter, to `main(argv)` with the sampler stubbed). Adds round trip so the builder and parser cannot drift.

**54. evaluate_builds_no_trainer** (integration). `evaluate(StubSampler, rows, params, view_root)` on fixture rows with trl `GRPOTrainer.__init__`, `VLLMGeneration._generate` and `_generate` monkeypatched to raise: passes; stub sampler called exactly once with `(rows, params)`. Output dict contains per-item and aggregate numbers; with `truncated=True` samples, they are scored as defined (pin against `train/rewards` golden). Edge: zero rows -> well-defined empty result (no NaN mean/ZeroDivision); n=1; a sample with empty text -> scored as parse failure, not exception.

**55. sample_scripts_archived** (unit). `experiments/phase2/sample_vllm.py` and `hf_sample.py` absent outside `archive/`; present under `archive/<name>/` with README; `archive/` not matched by pytest collection (`testpaths`), ruff, lint-imports config, `code_files()`.

---
## F. Verifier API (B9, B13)

**56. verifier_bit_equal_score_batch** (integration, CPU build). Fixed corpus (64 polygon texts of golden 3 plus adversarial: empty string, `"[]"`, 1 vertex, 2 vertices, self-intersecting bowtie, duplicate vertices, vertices at 0 and 1000, hole ring, 3 nested rings, truncated JSON, non-JSON, NaN text, 100 kB text). `[r.reward for r in Verifier("reward_v4", root).score(items)]` equals `score_batch(...)` with `float.hex` equality, each position; for items with `error` set, reward is `None` (spec) while `score_batch` has its own convention: assert `error` is not None exactly where parse fails. Order preserved; ids echoed; batch of 0 -> `[]`; batch of 1; duplicate ids allowed.

**57. verifier_batch_split_invariance** (property, hypothesis). Scoring items in one batch == scoring in two halves == one at a time (bit-equal). Permuting the batch permutes results identically. Why: catches batch-position dependent state (shared buffers, stream reuse).

**58. verifier_stamp** (unit). `stamp.api_version == API_VERSION == "1.0"`; `reward_name == "reward_v4"`; `reward_fingerprint == fingerprint()`; `checker_build in {"cuda","cpu"}` (== "cpu" on Mac); `code_version` non-empty str. Stamp frozen. `Verifier("nope", root)` raises `VerifyError` (`spec_mismatch`) or ValueError: pin; `Verifier("reward_v4", missing_dir)` raises at init, not at first score.

**59. unknown_key_whole_batch** (unit). Batch of 5 with the 4th key unknown: raises `VerifyError` with `code == "unknown_key"` and NO partial results; message names the key; other 4 not scored (assert scorer not called via spy). Key with wrong `VIEW_SEP` form (`"a:b"` vs `"a::b"`), empty key, key with unicode: unknown_key. Key format `"id::view"` resolves.

**60. parse_failure_is_a_result** (unit). Texts: `""`, `"not json"`, `"[[1,2]]"`, huge ints, floats out of 0-1000 grid, nested wrong depth. Each returns `ScoreResult` with `error` non-empty, `reward is None`, `terms == {}` or zeros (pin), no exception; the next item in the batch scores normally.

**61. verify_api_no_pydantic_no_torch** (integration, subprocess). Import `objdel.verify.api` with `-I`: `pydantic`, `torch`, `fastapi` not in `sys.modules`. Dataclasses frozen: mutate raises `FrozenInstanceError`.

**62. http_score_equals_inprocess** (integration, skip if no `serve` extra). Use `fastapi.testclient`/`httpx` on `objdel.verify.http.app`. `POST /v1/score` floats equal in-process (`float.hex` of JSON round-trip: JSON float repr must round-trip exactly; compare `float(x).hex()`). Every response body (200, 4xx, 413) has `api_version`. `GET /v1/spec` equals `Verifier.spec()`; `/livez` 200 always; `/readyz` 200 when ready and 503 before the Verifier is loaded (startup lifespan).

**63. http_errors** (integration, skip if no serve). Batch `MAX_BATCH + 1` -> 413, `error.code == "batch_too_large"`, `retryable == false`; batch exactly `MAX_BATCH` -> 200 (boundary; use cheap items); empty `items` -> 200 `[]` or 422 (pin); unknown key -> 4xx with `unknown_key`, retryable false; malformed JSON body -> 422 with the same error envelope shape `{"error":{"code","message","retryable"}}` (not FastAPI default `detail`; catches missing exception handler); wrong `reward` name -> `spec_mismatch`; simulated `device_error` (monkeypatch scorer to raise) -> 5xx with `retryable == true`; ALL other codes `retryable == false` (parametrise all four codes). Extra field in an item -> rejected (forbid) or ignored (pin).

**64. http_coverage_job_installs_serve** (unit). `.github/workflows/*` coverage job installs extra `serve` (grep `\[.*serve.*\]`); `lint-imports` step exists in CI (B1).

**65. store_and_pool_seams** (unit). `objdel.verify.store` exposes `preload, loader, VIEW_SEP == "::"`; `objdel.verify.pool` exposes `fork_pool, pool_size`. `loader(root)("missing::view")` raises KeyError-like with key; `preload(root, [], 4) == 0`; `preload` with threads=1 vs 8 returns equal count. `fork_pool(2, modules=["json"])` workers have `json` imported and `fork_pool(1)` imports nothing extra: check in the worker `sys.modules` difference against a baseline set; `fork_pool(1, ["no.such.module"])` raises at creation (not silently). `pool_size(cgroup)` with fixture cgroup files: `max` -> cpu_count; quota 200000/period 100000 -> 2; quota 150000/100000 -> 1 (floor, minimum 1); quota 10/100000 -> 1; missing file -> cpu_count; garbage content -> cpu_count; cgroup v1 layout if supported.

**66. verify_names_no_higher_layer** (unit, AST). In `objdel.verify` source: imports and ALL string literals contain no `objdel.policy|rl|sft|infra`. Self-test with a synthetic string.

**67. polygonize_single_entry** (unit/property). `objdel.geom.polygonize(mask, offset=(0,0))`: empty mask -> `[]`; single pixel -> one ring with 4 vertices enclosing area 1; solid rectangle -> 4 vertices; ring with a hole -> outer + hole (hole role marked per existing format); two blobs -> two rings; mask touching all four borders; full-image mask 512x512; `offset=(dy,dx)` shifts every vertex by exactly that amount and `offset` argument order is (row, col) vs (x, y): pin with an asymmetric mask. Property: rasterise(polygonize(mask)) == mask for random blobs (round trip through `to_masks`); dtype bool, uint8 0/1, and 0/255 inputs behave as documented (non-bool 255 treated as True). AST: production code calls only `geom.polygonize`, not the old `checker.polygonize` internals (no import of `objdel.geom.<private polygonize module>` outside geom).

---
## G. Move invariance goldens and host syncs (B10, B11)

**68. golden1_reward_spec** (unit). `sha256(repr((astuple(REWARD_V4), [min_area(a) for a in (0, 999, 1001, 10**4, 10**6)])))` equals the recorded hex string. Full `fingerprint()` hex equal only if numpy/OpenCV/scipy versions equal recorded versions; otherwise skip with message naming the versions. Never recompute the golden in the test.

**69. golden2_target_set_hash** (unit). `target_set_hash` of the fixture items equals recorded (per fixture item id).

**70. golden3_score_batch_corpus** (integration). 64-text corpus, CPU build, each float `float.hex` equals `invariance.<sys.platform>-<machine>.json`; missing file -> `pytest.skip("no golden for darwin-x86_64")`; darwin-arm64 file must exist (test fails, not skips, when `sys.platform=="darwin" and machine=="arm64"` and file missing: binding platform). Also checks key set equal (no corpus item dropped).

**71. golden4_fastgrpo_update** (integration, CPU, tiny model, seed 0). Loss, grad L2 norm, parameter-change L2 norm as `float.hex` equal golden per platform; run the update twice in one process: identical (determinism precondition; if not, the golden is meaningless). Variation: shared_prefix False and True each have a golden.

**72. golden5_chain_argv** (see 35). Include all four specs; `rl_shifted` etc. under the rename map: `experiments.phase2.X` -> new module path. Test that the rename map itself is total: every `-m <module>` token in goldens maps to an existing module (`importlib.util.find_spec`).

**73. golden6_grpo_config** (see 41).

**74. host_sync_budget** (unit). For each module in src importing torch: count `\.(item|cpu|tolist|numpy)\(\)` regex matches (raw text, as spec says; include comments/docstrings consistently with the recorded counts) <= recorded in `tests/goldens/phase_e/host_syncs.json`; moved files via old->new map. New torch-importing module with no record -> FAIL (count must be recorded as 0 or explicit). Recorded entry for a missing module -> fail (stale). Negative self-test: synthetic source with `x.item()` over budget flagged. Edge: `.item( )` with space not matched by spec regex; `.cpu().numpy()` counts 2.

**75. paths_and_secrets** (unit). `repo_root()` is the dir containing `pyproject.toml`; `OBJDEL_REPO=tmp` takes priority (tmp with pyproject) ; `OBJDEL_REPO` pointing to a dir with no pyproject -> raises or falls back (pin); from a different cwd gives the same result. `OBJDEL_SECRETS` unset -> `secrets.path() == repo_root()/".secrets"`; set -> that value (the autouse fixture sets it; use `monkeypatch.delenv`). AST: no module under `src/objdel/infra/` uses `__file__`. `code_files()` lists no path starting `.claude/` or `archive/` (also nested `x/.claude/y`), contains `tests/test_gpu_scoring.py`; deterministic sorted order; no path outside repo, no `data/`, no `.secrets` (secret leakage guard).

---
## H. Whole-suite, performance, GPU

**76. fast_suite_gate** (end-to-end, CI). `make test`: 0 failed; passed >= 1447 - removed (the commit's list file must name each removed test; the gate script diffs collected node ids of base vs head: every missing id is listed or fail). `make lint` exit 0. `pytest --collect-only` yields no import errors (collection errors count as failures, a classic stale-import symptom).

**77. check_private_names_time** (performance). < 10 s on Mac in a fresh subprocess (cold import of trl). Use wall time with a generous margin in CI (retry once); also assert it imports no vLLM.

**78. import_time_floor** (performance). `python -I -X importtime -c "import objdel.verify"` total < budget recorded (e.g. 1.5 s on Mac; budget [estimate], recorded at implementation) and loads no torch. Guards accidental heavy imports.

**79. score_throughput_no_regression** (performance, slow). `Verifier.score` on 1,024 corpus items: wall time <= 1.10x recorded base-commit time for `score_batch` on the same machine (relative, measured in same session by both calls to avoid machine noise; assert verifier/score_batch ratio < 1.10). Verifier must add no per-call stamp/fingerprint recomputation (fingerprint cached: assert `fingerprint` called once across 10 score calls via spy).

**80. gpu_marked** (GPU, skipped without a GPU). `checker_build == "cuda"`, CUDA vs CPU `Verifier` rewards bit-equal (`float.hex`) on the corpus (existing test_gpu_scoring parity extended to the new import path); `device=1` on a one-GPU box -> `VerifyError(device_error)`; box smoke for G13 (`trl_private_names`, `IMPORT_CHECK` exit 0). Mark `@pytest.mark.gpu`; skip when `torch.cuda.is_available()` false or `_gpucheck` built without CUDA.

**81. regression_today_numbers** (integration, data/ckpt-gated, skip when absent). Pinned results from the repo history: test IoU 0.645, F 0.816, precision 0.695, recall 0.748 (final_eval on `rl-full-3ep`) recomputed with the new `evaluate` path through `Verifier` on the same samples must equal within bootstrap tolerance 0 for per-item rewards (bit-equal on the same sampled completions). Only the scoring step; no sampling (sampling needs a GPU). Mark slow.

**82. cross_cutting_negative_self_tests** (unit). Every AST/regex guard above (tests 6, 8, 12, 13, 15, 16, 17, 27, 66, 74) ships a tiny synthetic-source test that the guard FAILS on a planted violation. Why: static guard tests are the most likely to be vacuously green (wrong path glob, empty file list). Also assert each guard scanned >= N files (non-empty glob sanity: e.g. src has > 30 modules).

---
## Pytest sketches for the 5 most important tests

Priority reasoning: (1) move invariance of reward (core claim, silent failure), (2) no-ML-stack isolation, (3) drift guard actually fires, (4) fallback equals TRL / no global patch, (5) spec to argv goldens.

```python
# tests/phase_e/test_moves.py
import json, subprocess, sys, platform
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[2]
GOLD = ROOT / "tests" / "goldens" / "phase_e"

# 1. verifier == score_batch == recorded golden, bit-equal
def test_verifier_bit_equal_to_score_batch_and_golden(corpus_items, view_root):
    from objdel.verify.api import Verifier
    from objdel.verify.reward import score_batch
    got = [r.reward for r in Verifier("reward_v4", view_root).score(corpus_items)]
    ref = score_batch([(i.key, i.text, i.truncated) for i in corpus_items])  # adapt to real signature
    ok = [(i, r) for i, r in zip(got, ref) if i is not None]
    assert [g.hex() for g, _ in ok] == [float(r).hex() for _, r in ok]
    f = GOLD / f"invariance.{sys.platform}-{platform.machine()}.json"
    if not f.exists():
        if (sys.platform, platform.machine()) == ("darwin", "arm64"):
            pytest.fail("binding platform has no golden")
        pytest.skip(f"no golden for {f.name}")
    assert [x.hex() for x in got if x is not None] == json.loads(f.read_text())["score_batch"]

# 2. no ML stack below policy, fresh process per package
LOW = ["objdel.core", "objdel.geom", "objdel.verify", "objdel.data"]
BAD = {"torch", "trl", "transformers", "vllm", "peft", "datasets", "accelerate"}
@pytest.mark.parametrize("pkg", LOW)
def test_no_ml_stack(pkg):
    import importlib, pkgutil
    code = (
        "import importlib, pkgutil, sys\n"
        f"p = importlib.import_module({pkg!r})\n"
        "for m in pkgutil.walk_packages(p.__path__, p.__name__ + '.'):\n"
        "    importlib.import_module(m.name)\n"
        "print(sorted({k.split('.')[0] for k in sys.modules} & " + repr(BAD) + "))\n"
    )
    out = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "[]", out.stdout

# 3. drift guard fires in every negative case
@pytest.mark.parametrize("name", ["_get_per_token_logps_and_entropies", "_generate_and_score_completions",
                                  "_prepare_inputs", "_compute_loss", "_generate", "_calculate_rewards"])
def test_missing_method_raises(name, monkeypatch):
    from objdel.rl import trl_adapter as A
    from trl.trainer.grpo_trainer import GRPOTrainer
    owner = GRPOTrainer if hasattr(GRPOTrainer, name) else None
    if owner is None:
        pytest.skip("method lives on another TRL class; covered by the vllm-generation case")
    monkeypatch.delattr(GRPOTrainer, name)
    with pytest.raises(A.TrlDriftError) as e:
        A.check_private_names()
    assert name in str(e.value) and "1.17.1" in str(e.value)

def test_version_changed_and_source_changed_raise(monkeypatch):
    from objdel.rl import trl_adapter as A
    import trl
    monkeypatch.setattr(trl, "__version__", "1.17.10")
    with pytest.raises(A.TrlDriftError, match="1.17.1"):
        A.check_private_names()

# 4. fallback == TRL output; no global patch (tiny model fixture from conftest)
@pytest.mark.parametrize("fallback", ["eval", "loss_type", "entropy", "accum"])
def test_fallback_equals_trl_and_global_untouched(fallback, tiny, tmp_path):
    import torch, trl.trainer.grpo_trainer as g
    original = g.shuffle_sequence_dict
    fast = build_tiny_trainer(tiny, tmp_path, plain=False, fallback=fallback)
    plain = build_tiny_trainer(tiny, tmp_path, plain=True, fallback=fallback)
    batch = one_generation_batch(tiny)
    torch.manual_seed(0); a = fast._prepare_inputs(batch)
    torch.manual_seed(0); b = plain._prepare_inputs(batch)
    assert_same_structure_and_tensors(a, b)       # torch.equal on every tensor, same order
    assert g.shuffle_sequence_dict is original    # never swapped, even transiently left behind

# 5. spec -> argv golden, string for string
@pytest.mark.parametrize("spec", ["rl_shifted", "rl_full", "r5_shared"])
def test_to_argv_matches_pre_phase_e_chain(spec):
    from objdel.core.spec import load_spec, to_argv
    rs = load_spec(ROOT / "specs" / f"{spec}.yaml")
    gold = json.loads((GOLD / "chain_argv.json").read_text())[spec]   # recorded on base commit, after rename map
    got = {s.name: to_argv(s, "RD", "DATA") for s in rs.stages}
    assert got == gold
```

Notes for the implementer of the test files: helpers named above (`corpus_items`, `view_root`, `build_tiny_trainer`, `one_generation_batch`, `assert_same_structure_and_tensors`) go in `tests/phase_e/conftest.py`, reusing `tests/conftest.py` `tiny`. Goldens are recorded on the pre-Phase-E commit only, never regenerated from new code.

## Open questions found while writing (spec ambiguities)
1. `save_steps` default 50 vs `steps < 50`: does the default trip the `save_steps > steps` rule? (test 33)
2. `truncated(ids, stops)` exact semantics for stop id not at the end (test 51).
3. Empty `CmdStage.argv` and empty `items` in HTTP: accept or reject (tests 36, 63).
4. Does the SHA-256 of a TRL method include comments/whitespace? Test 19d assumes yes.
5. `ScoreResult.terms` content when `error` is set (test 60).
