# Python tooling the plan adopts: exact usage at the pinned versions

Date: 2026-11-09. Tags: [run] = I ran it in a scratch venv this session (Python 3.13, macOS) with the version named;
[doc] = read in official docs (URL given); [src] = read in installed source; [unverified] = not checked.
Scratch work lived under the session scratchpad, nothing in the repo changed.
PyPI latest on 2026-11-09 [run: pypi.org/pypi/<pkg>/json]: import-linter 2.16 (2026-11-06), pydantic 2.14.3,
pydantic-settings 2.16.2, zensical 0.1.0, mkdocs-material 9.8.1, vulture 2.16, ruff 0.17.4, huggingface-hub 2.4.0,
pytest 9.2.0, mkdocstrings-python 2.0.10, uv 0.13.7. Plan pins differ from latest in 3 places: ruff 0.17.1 (latest
0.17.4), pydantic 2.14.2 in box locks (latest 2.14.3; behaviour below tested on 2.14.3), huggingface-hub 2.2.0 in
`locks/train.txt` and `vllm_native.txt` but 2.4.0 in `locks/vllm.txt`. All HF signatures below are identical in 2.2.0
and 2.4.0 [run: inspect.signature on both].

## 1. import-linter 2.16 (plan: E0, E3a-c, E5, E8; CLEAN)

Requires Python >=3.11, deps `grimp>=3.19`, `click`, `rich` [run: PyPI 2.16 JSON].

### Run
- CLI `lint-imports` (entry point `importlinter.cli:lint_imports_command`); options: `--config`, `--contract ID`
  (repeatable), `--cache-dir`, `--no-cache`, `--debug`, `--no-logo`, `--show-timings`, `--verbose` [run: --help].
- Exit code: 0 all kept, 1 any broken or config error [run]. Config is read from `pyproject.toml`
  `[tool.importlinter]` (or `.importlinter` / `setup.cfg`) in the cwd [doc: https://import-linter.readthedocs.io/en/stable/usage/ - page 404 via fetch; confirmed by run].
- **There is no `python -m importlinter`** (`No module named importlinter.__main__`) [run]. From pytest use either
  `subprocess.run(["lint-imports", "--no-cache", "--no-logo"], cwd=ROOT)` (resolve the script next to
  `sys.executable`, not via PATH), or in-process `from importlinter.cli import lint_imports; rc = lint_imports(no_cache=True, no_logo=True)`
  which returns int 0/1 [run; src: importlinter/cli.py:188]. In-process inserts `os.getcwd()` into `sys.path`
  (cli.py:219): run it with cwd = repo root and expect the side effect. Prefer subprocess in tests for isolation.
- Use `--no-cache` in negative tests that edit files in a scratch tree; the default cache dir is `.import_linter_cache`
  in the cwd (add to `.gitignore`) [unverified name; the flag `--cache-dir` is verified].
- `root_package` must be importable (editable install, or cwd on path). `include_external_packages = true` is
  required if any contract names an external package (forbidden `torch`) [doc: https://import-linter.readthedocs.io/en/stable/contract_types/forbidden/].

### Layers contract [doc: https://import-linter.readthedocs.io/en/stable/contract_types/layers/; every claim re-run]
- Fields: `layers` (high to low, names relative to `containers`; no wildcards), `containers` (wildcards allowed),
  `exhaustive` (default false), `exhaustive_ignores`, `ignore_imports`, `unmatched_ignore_imports_alerting`,
  `broken_contract_guidance`.
- Siblings: `"a | b"` independent (may not import each other); `"a : b"` may import each other; mixing separators on one
  line is invalid. Optional layer: wrap in parentheses `"(medium)"`; a non-optional layer that does not exist as a
  module is an error ("Missing layer in container 'objdel': module objdel.nope does not exist") [run, exit 1].
- `"(rl) | (sft) | (nope)"` (all parenthesised, independent) is accepted and KEPT [run]. This answers plan E0 negative test (c):
  the `(x) | (y)` form works at 2.16, so E5 may use it for the interim list while one package is missing.
- Indirect imports count (chains via other modules) [doc].
- `exhaustive = true` needs `containers`; every child module of the container not named as a layer fails unless listed in
  `exhaustive_ignores` (`"_gpucheck"`, `"secrets"` as bare child names) [run: planted `objdel/extra.py` gives
  "The following modules are not listed as layers: - objdel.extra"]. A top-level `.py` module counts as a child.
- Upward import is reported as `objdel.verify is not allowed to import objdel.rl: - objdel.verify.up -> objdel.rl (l.1)` [run].
- Missing non-optional layer: exit code 1 [run, re-checked without a pipe]. (My first capture showed 0 because of `| head`.)

### Forbidden contract [doc: https://import-linter.readthedocs.io/en/stable/contract_types/forbidden/]
- Fields: `source_modules`, `forbidden_modules` (both accept wildcards), `allow_indirect_imports` (default false:
  an indirect chain `objdel.data -> objdel.policy -> torch` also breaks), `as_packages` (default true: a listed module
  covers its descendants).
- External packages: only root-level names (`torch`, not `torch.nn`) and the top-level `include_external_packages = true`.
- **Function-level imports ARE detected** [run: `def f(): import torch` in `objdel/data/lazy.py` gave
  `objdel.data.lazy -> torch (l.2)`]. So plan negative test (d) passes; keep E-G4 anyway for `importlib.import_module("torch")`
  (string imports are invisible to grimp).
- Because indirect imports count, `pure-core` source `objdel.verify` fails if any `verify` module imports a module that
  imports torch, even through `objdel.policy`. Set `allow_indirect_imports = true` only if you want direct edges only.

### ignore_imports [doc: https://import-linter.readthedocs.io/en/stable/contract_types/]
- Form `"pkg.a.importer -> pkg.b.imported"`. `*` = one module name (no subpackages), `**` = subpackages too;
  the wildcard must replace a whole module name: `pkg.foo*` is invalid.
- **An ignore that matches no import is an error by default** [run: "No matches for ignored import objdel.data.zz ->
  objdel.verify.zz", exit 1]. `unmatched_ignore_imports_alerting = "error" | "warn" | "none"` changes it. Plan R: each move
  commit must drop the ignore its move removes. If the E0 CI needs to survive a missing ignore, set `"warn"` on that
  one contract; do not set `"none"`.
- `broken_contract_guidance` is printed under the failure text [run: "see docs" appeared].

### Plan checks
- E0 TOML in the plan is valid for 2.16 except: `exhaustive = true` only at E8 (needs containers; E0 has them).
  The E0 contract with `"reward | targets"` etc. is fine.
- A contract with `source_modules` naming a module that does not exist errors (plan says so; consistent with the missing-layer behaviour) [unverified for forbidden].
- Add to `Makefile` `lint:` the call `.venv/bin/lint-imports`; the CI step needs the project installed (editable) so grimp finds `objdel`.

## 2. pydantic 2.14 and PyYAML duplicate keys (plan: E2, E6)

Docs: https://docs.pydantic.dev/latest/concepts/models/ , /concepts/unions/#discriminated-unions , /concepts/strict_mode/ ,
/api/config/ . Every behaviour below [run, pydantic 2.14.3, PyYAML 6.0.3].

- `class Base(BaseModel, frozen=True, extra="forbid")` (class kwargs) is valid; `ConfigDict(extra="forbid", strict=True)` equivalent.
- Discriminated union: `Stage = Annotated[Rl | Sft, Field(discriminator="kind")]` with `kind: Literal["rl"]` in each member.
  Error types: unknown tag `union_tag_invalid` at loc `('stages', 0)`; missing tag `union_tag_not_found`;
  unknown extra key inside a member `extra_forbidden` at loc `('stages', 0, 'rl', 'typo')` (the tag name is part of loc).
- `extra="forbid"` is NOT inherited by field-typed submodels you forgot to declare: each nested `BaseModel` needs it
  (hence plan P2's tree walk). Inheriting from one `Base` solves it.
- **Lax mode coerces YAML scalars**: `n: '3'` and `n: 3.0` validate to `3` [run]. For typed specs use
  `Model.model_validate(data, strict=True)` or `ConfigDict(strict=True)`; strict rejects `'3'`, `3.0`, and `True` for `int`
  (`int_type`) [run]. `StrictInt` also exists per-field.
- **YAML 1.1 pitfalls with PyYAML** [run: `yaml.safe_load`]: `lr: 1e-4` loads as the string `'1e-4'` (PyYAML float
  needs a dot: write `1.0e-4`); `yes`/`no`/`on`/`off` load as bool; `010` loads as octal 8. Under `strict=True` the string `1e-4`
  fails with `float_type` for a float field; in lax mode it is silently parsed. Pick strict and write `1.0e-4` in specs.
- Frozen model hashing: `hash()` of a frozen model with a `list` field raises `TypeError: unhashable type: 'list'`;
  use `tuple[...]` fields to hash [run]. Matches plan E2 note on `Mapping` fields.
- `model_validate_json` runs JSON-mode validation (different coercion rules from Python-mode); use `model_validate(yaml_dict)` for YAML.

### Duplicate YAML keys (PyYAML ignores them silently: last wins)
Custom loader, tested [run]: it raises on duplicate scalar keys, and handles merge keys.
```python
import yaml
class UniqueKeyLoader(yaml.SafeLoader): pass
def _construct_mapping(loader, node, deep=False):
    seen = set()
    for key_node, _ in node.value:
        if key_node.tag == "tag:yaml.org,2002:merge":   # `<<: *anchor`: skip, else construct_object raises
            continue
        key = loader.construct_object(key_node, deep=True)
        if key in seen:
            raise yaml.constructor.ConstructorError(None, None, f"duplicate key {key!r}", key_node.start_mark)
        seen.add(key)
    return yaml.SafeLoader.construct_mapping(loader, node, deep)
UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_mapping)
data = yaml.load(text, Loader=UniqueKeyLoader)     # error message includes line and column
```
Pitfall: without the merge-tag skip, `<<: *a` raises "could not determine a constructor for the tag ...:merge" [run].
Unhashable keys (a list as key) would fail at `key in seen`; specs never use them. Do not use `yaml.BaseLoader` for specs
(everything becomes str); plan DOC4 uses BaseLoader only for `nav` page-path reading, which is fine.

## 3. pydantic-settings 2.16.2 (plan: E6 `CliApp.run`)

Docs: https://docs.pydantic.dev/latest/concepts/pydantic_settings/ (section "Command Line Support"). [run unless noted]
- `from pydantic_settings import BaseSettings, SettingsConfigDict, CliApp`.
  `CliApp.run(Model, cli_args=[...])` parses the list (not `sys.argv`) and returns the model instance.
- `SettingsConfigDict(cli_kebab_case=True, cli_exit_on_error=False, cli_implicit_flags=True, env_prefix="OBJDEL_", extra="forbid", frozen=True)` accepted.
- With `cli_kebab_case=True`, `--max-steps 5` works; `--max_steps 5` is rejected: `SettingsError: error parsing CLI: unrecognized arguments: --max_steps 5`. An unknown flag raises the same `SettingsError` (with `cli_exit_on_error=False`; otherwise argparse exits 2).
- `list[float]` flag takes JSON text `--lr-list "[1.0,2.0]"`; a nested model takes JSON `--sub '{"a":4}'`; `cli_implicit_flags=True` makes `--dry-run` a bool flag.
- Env vars (`OBJDEL_MAX_STEPS=9`) apply both to `S()` and to `CliApp.run(S, cli_args=[])`. Precedence default: init/CLI args > env > defaults.
- `BaseSettings` with a `float | None` field and an env var set: plan E6 already measured the `OBJDEL_SESSION_START=123.0` case; I did not re-run it.

## 4. Zensical 0.1.0 and Material fallback (plan: DOC3)

- Python >=3.11. Commands: `zensical build`, `zensical serve`, `zensical new` [run: --help]. `build` options:
  `-f/--config-file PATH`, `-c/--clean` (clean cache), `-s/--strict` ("abort the build on warnings") [run].
  `--site-dir` is not supported; use `site_dir` in the config [doc: https://zensical.org/compatibility/cli].
- **`--strict` works at 0.1.0** (a web search snippet said Zensical ignores it; that is stale): a page with `[bad](nope.md)`
  printed "Warning: page does not exist ... 1 issue found", then `RuntimeError: Aborted because --strict flag is set`, exit code 1.
  Without `--strict` the same build exits 0 [run]. The traceback is noisy but the exit code is the gate.
- Reads `mkdocs.yml` directly (FAQ: reads it indefinitely [doc: https://zensical.org/docs/community/faqs/ via search summary; page fetch 404]). Native config is `zensical.toml`.
  Config keys: `site_dir` (relative to config file), `docs_dir` (cannot be `.`), `nav` supported; `hooks` NOT supported; `theme.name`,
  `validation`, `strict` keys not addressed by docs [doc: https://zensical.org/compatibility/configuration/, https://zensical.org/docs/setup/basics/].
- Plugins [doc: https://zensical.org/compatibility/plugins]: `search` (built in, on even if unlisted), `mkdocstrings` (since 0.0.11; install
  `mkdocstrings-python` separately), `api-autonav` (since 0.0.66; name is `api-autonav`), `awesome-nav`, `macros` supported;
  `git-revision-date` not. **A plugin that is not supported is silently ignored, not an error** [doc; run: `some-unknown-plugin` listed, build passed].
  Consequence: a typo in a plugin name will not fail CI; DOC4's test should assert the plugin list in `mkdocs.yml` itself.
  Supported plugins reject unknown settings with a config error [doc].
- Build writes a `.cache/` directory in the project root [run]; plan DOC already adds `/site/` to `.gitignore`: also add `/.cache/`.
- Material fallback: `mkdocs-material` 9.8.1 prints a MkDocs 2.0 warning; silence with `NO_MKDOCS_2_WARNING=1` (warning began at 9.7.2);
  Material critical maintenance runs through 2027-05-05; 9.7.5+ limits `mkdocs<2` [doc: https://squidfunk.github.io/mkdocs-material/blog/2026/02/18/mkdocs-2.0].
  Fallback build is `mkdocs build --strict` with `mkdocs-api-autonav` and `mkdocstrings-python` [plan; not re-run].
- Keep `mkdocs.yml` plain: no `hooks`, no plugin Zensical lacks, no `!!python/name` tags in markdown_extensions (a YAML loader that is not unsafe cannot read them; `yaml.BaseLoader` in tests also rejects them) [unverified for Zensical].

## 5. ruff 0.17.x docstrings and TID251 (plan: E0, DOC5)

Docs: https://docs.astral.sh/ruff/rules/banned-api/ , https://docs.astral.sh/ruff/settings/#lint_pydocstyle_convention . [run: ruff 0.17.1]
- Banned-api table syntax is valid exactly as the plan writes it:
  `[tool.ruff.lint.flake8-tidy-imports.banned-api]` then `"trl.trainer.grpo_trainer".msg = "..."`. TID251 must also be in `select`.
- What TID251 flags [run]: `from trl.trainer import grpo_trainer` (col of the imported name), `from trl.trainer.grpo_trainer import X`,
  `import trl.trainer.grpo_trainer as g`, `import trl.generation.vllm_generation`, `from objdel import _gpucheck`,
  `from objdel._gpucheck import score`, `import objdel._gpucheck`, and **function-level** imports. Not flagged: `from trl import trainer`,
  `importlib.import_module("trl.trainer.grpo_trainer")` (strings are invisible).
- `per-file-ignores = {"src/objdel/ok.py": ["TID251"]}` frees the file for TID251 only; `"tests/**" = ["D"]` frees all D rules [run]. Plan E0 negative test passes: a tests file with the banned import was reported.
- Docstrings: `select = [..., "D100","D101","D102","D103","D104","D106"]` with `[tool.ruff.lint.pydocstyle] convention = "google"`
  gives D101 (class), D102 (method), D103 (function) on missing docstrings [run]. D105 (magic) and D107 (`__init__`) are not selected by the plan: fine.
  Convention `google` switches off D203, D204, D213, D215, D400, D401, D404, D406, D407, D408, D409, D413 only; it does not touch D1xx. Selecting full `D`
  would add D205/D209 churn (plan measured 444 hits).
- `ruff check --isolated --select D100,... path` ignores `pyproject.toml` (plan's measuring command is valid).
- Pin note: `ruff==0.17.1` in the plan; 0.17.4 is latest [run: PyPI]. Rule counts in the plan are for 0.17.1: do not bump ruff inside E/DOC.

## 6. vulture 2.16 (plan: CLEAN12)

Source: https://github.com/jendrikseipp/vulture ; [run: vulture 2.16, and [src]].
- Config in `[tool.vulture]`: `paths`, `min_confidence`, `ignore_names`, `ignore_decorators`, `exclude`, `sort_by_size`, `make_whitelist`. `python -m vulture` with no args reads `paths` from pyproject [run].
- **Exit codes: 0 clean, 3 dead code found** [run]. Output lines `path:line: unused function 'x' (60% confidence)`.
- `ignore_names = ["on_*", "cmd_*"]` is a glob (fnmatch) on the bare name: `on_step` and `cmd_go` were not reported, `other` was [run].
  `ignore_decorators = ["@pytest.fixture"]` takes `@`-prefixed patterns [doc].
- A class whose only members are ignored is still reported as `unused class 'CB'` (60%) [run]: if the plan's Trainer subclasses show up, list the class name.
- **`exclude` semantics trap** [src: vulture/core.py `prepare_pattern`]: a pattern WITHOUT `*?[` is turned into `*pattern*` (substring, case-insensitive);
  a pattern WITH a wildcard is matched with `fnmatch` against the full path as scanned [run: `"src/skip*"` did NOT exclude `src/skipme.py`; `"*/skipme.py"` and `"src/skipme.py"` did].
  The plan's `--exclude "experiments/phase1/lab,..."` (no wildcard) therefore works as a substring; never write a wildcard pattern without a leading `*`.
- `# noqa` comments: bare `# noqa` (or `# noqa: vulture`, which parses as bare) silences all vulture findings on the line; codes are `V104` (unused import) / `V107` (unused local) only [src: vulture/noqa.py].
  Do not add `# noqa: vulture` comments in the repo: ruff PGH/RUF100 may flag them if ever selected.
- Whitelist file: pass as an extra path (`vulture src wl.py`), lines like `unused_fn  # noqa` or `_.attr`; `--make-whitelist` prints one [run]. The plan rejects a whitelist file (ruff F821/B018): `ignore_names` covers it.

## 7. pytest strict markers (plan: Phase A, DOC4)

Docs: https://docs.pytest.org/en/stable/how-to/mark.html [run: pytest 9.2.0].
- Registration: `[tool.pytest.ini_options] markers = ["gpu: needs a CUDA device", "slow: ..."]`; strictness: `addopts = "--strict-markers"`.
  With a typo `@pytest.mark.gpuu` collection fails: `'gpuu' not found in 'markers' configuration option`, exit 2 (whole file errors) [run].
- pytest 9 also accepts a native table `[tool.pytest]` with `strict_markers = true` [run]. The repo pins `pytest>=8.0`:
  keep `[tool.pytest.ini_options]`, which works on 8 and 9.
- Markers named only in `-m "expr"` are not validated by `--strict-markers` [unverified]; register `docs_build`, `post_commit`, `gpu`, `slow`, `data` in the same commit that adds `--strict-markers`.
- `pytest --markers` lists them with descriptions [run].

## 8. huggingface_hub 2.2.0 / 2.4.0 (plan: B1.5, DOC11, `hf_push`)

Docs: https://huggingface.co/docs/huggingface_hub/package_reference/hf_api ; signatures [run: inspect, identical in 2.2.0 and 2.4.0].
- `HfApi.upload_folder(*, repo_id, folder_path, path_in_repo=None, commit_message=None, commit_description=None, token=None, repo_type=None, revision=None, create_pr=None, parent_commit=None, allow_patterns=None, ignore_patterns=None, delete_patterns=None, run_as_future=False)`.
  `path_in_repo` = target directory in the repo (default root), e.g. `"rows/lidc-train-b1"`. `ignore_patterns` always gets the default `.git` ignores added [src].
  `delete_patterns` deletes remote files matching the patterns in the same commit (never `.gitattributes`): use it only on purpose; a wrong pattern deletes a checkpoint.
  `create_pr=True` together with `revision` is an error ("Don't set `revision` when `create_pr=True`") [src].
  Large folders are committed in several batches with " (part N)" messages (plan 4428); `parent_commit` applies to the first only [doc string].
- `parent_commit` = optimistic concurrency: with `create_pr=False` the commit fails if `revision` no longer points at `parent_commit`; 7-char shorthand allowed [src docstring]. Use it in `hf_push` to stop a second writer.
- `revision` for writes names a branch ("git revision to commit from", default `main`). Whether a raw commit SHA is accepted as a write target is [unverified]: pin by `parent_commit`, not by passing a SHA as `revision`.
  Create a pinned reference with `create_tag(repo_id, *, tag, tag_message=None, revision=None, repo_type=None, exist_ok=False)` or `create_branch(repo_id, *, branch, revision=None, ...)`.
- `HfApi.list_repo_commits(repo_id, *, repo_type=None, token=None, revision=None, formatted=False) -> list[GitCommitInfo]`, newest first;
  `GitCommitInfo(commit_id, authors, created_at, title, message, formatted_title, formatted_message)`. Index `[0]` is the head, `[-1]` the initial system commit [src docstring]. No page-size argument; it lists all commits.
- `snapshot_download(repo_id, *, repo_type=None, revision=None, cache_dir=None, local_dir=None, etag_timeout=10, force_download=False, token=None, local_files_only=False, allow_patterns=None, ignore_patterns=None, max_workers=8, tqdm_class=None, dry_run=False)`.
  `revision` accepts a branch, tag or full commit SHA for reads (pin runs this way). It returns `str(os.path.realpath(local_dir))` when `local_dir` is set. `dry_run=True` downloads nothing and returns `list[DryRunFileInfo]` (usable for the preflight "does the pattern match" test).
  Pattern matching is `filter_repo_objects` [run]: `"rows/lidc-train-b1/**"`, `["rows/lidc-train-b1/*"]` and `"rows/lidc-train-b1/"` (trailing slash is expanded) all selected exactly the 2 files under that folder and not `rows/other/x`. Matching uses fnmatch, so `*` crosses `/`.
  Repo-relative paths are kept under `local_dir` (plan, verified there against `_snapshot_download.py`).
- `upload_file(*, path_or_fileobj, path_in_repo, repo_id, token=None, repo_type=None, revision=None, commit_message=None, commit_description=None, create_pr=None, parent_commit=None, run_as_future=False)`: all keyword-only; plan's fake `FakeApi.upload_file(**k)` reads `k["path_in_repo"]`, `k["revision"]` - consistent.
- `create_commit(repo_id, operations, *, commit_message, ..., revision=None, parent_commit=None, num_threads=5)` for multi-file atomic commits.
- Token: pass `token=` explicitly or rely on `HF_TOKEN`; plan forbids argv. [doc: https://huggingface.co/docs/huggingface_hub/quick-start#authentication]
- Hub CLI is `hf` (the old `huggingface-cli` is deprecated) [skill: huggingface-skills:hf-cli, description]. Not used by the plan's code.

## 9. uv tool install with a pinned version (plan: CLEAN13)

Docs: https://docs.astral.sh/uv/concepts/tools/ [doc, quoted]. uv 0.13.7 latest.
- `uv tool install vastai==1.11.0` is the documented form (`uv tool install ruff==0.5.0` is the doc example); `pkg@1.11.0` is also accepted by `uv tool install`.
- Each tool gets an isolated virtualenv, "always run isolated from the project": vastai's own pins (e.g. Pillow) never touch `.venv`.
- `uv tool upgrade` respects the constraint given at install ("==1.11.0" stays); to change the version, run `uv tool install` again with the new pin.
  Add `--force` if replacing a different build of the same version [unverified flag name for this case; `uv tool install --force` exists].
- Executables go to the uv tool bin dir; `uv tool update-shell` adds it to PATH. `--python 3.12` picks the interpreter.
- Uninstall (plan's rollback) `uv tool uninstall vastai` and `uv tool list` were not on the fetched page [unverified here; standard uv subcommands].
- `uvx vulture==2.16 ...` as written in the plan (CLEAN baseline step): the page documents `uvx pkg@latest`; whether `uvx pkg==ver` parses is [unverified]. Safer: `uvx vulture@2.16`, or use `.venv/bin/python -m vulture` after the dev extra is installed.

## 10. CITATION.cff schema 1.2.0 (plan: DOC11)

Schema: https://github.com/citation-file-format/citation-file-format/blob/1.2.0/schema.json (fetched raw, tag 1.2.0) [run].
- Required top-level keys: `authors`, `cff-version`, `message`, `title`. `cff-version` must match `^1\.2\.0$`.
  `type` enum `software|dataset` (default software). Optional: `abstract`, `commit`, `contact`, `date-released`, `doi`, `identifiers`, `keywords`, `license`, `license-url`, `preferred-citation`, `references`, `repository`, `repository-artifact`, `repository-code`, `url`, `version`.
  Top level has `additionalProperties: false`: an unknown key fails validation.
- `authors` items: person (`family-names`, `given-names`, `orcid`, `email`, `affiliation`, ...) or entity (`name`). A person has no required key in the schema.
- `license` must be an SPDX identifier: `MIT` passed, `mit-license` failed [run via `cffconvert --validate`]. Match the repo `LICENSE`.
- `version` string should equal `pyproject.toml` `[project].version` (plan DOC4 test 9); `date-released` is `YYYY-MM-DD`; do not add `date-released` unless a test keeps it fresh.
- Validate in CI/tests: `pip install cffconvert` then `cffconvert -i CITATION.cff --validate` prints "Citation metadata are valid according to schema version 1.2.0." and exits 0 [run]. DOC4 uses PyYAML only (no schema): either accept that, or add `cffconvert` to the dev extra and call it from the test [decision for DOC11 owner].
- Plan item "no AI tool in authors": satisfied; never list a model as author (also matches CLAUDE.md git rule).

## Pitfall summary (for the implementer)
1. import-linter unmatched `ignore_imports` is an error; drop ignores in the same commit as the move.
2. No `python -m importlinter`; use the `lint-imports` script or `importlinter.cli.lint_imports()`.
3. Forbidden contracts see function-level imports and indirect chains; string imports are invisible to import-linter and ruff.
4. YAML: `1e-4` is a string, `yes` is bool, duplicates pass silently without the custom loader (merge keys need the tag skip).
5. pydantic lax mode coerces `'3'` to `3`; use strict validation for specs; frozen + list field is unhashable.
6. Zensical ignores unsupported plugins silently; `--strict` is real (exit 1 on warning). Add `.cache/` to `.gitignore`.
7. vulture `exclude` with a wildcard needs a leading `*`; no wildcard means substring.
8. Pin the HF write by `parent_commit`; pin reads by `revision=<sha>`; `allow_patterns` with a trailing `/` works.
9. Locks carry different huggingface-hub versions (2.2.0 and 2.4.0); APIs used here are identical.
