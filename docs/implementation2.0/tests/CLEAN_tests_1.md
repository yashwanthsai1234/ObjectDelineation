# CLEAN phase: test cases, set 1 (written from CLEAN.spec.md only, no design context)

## Skills used
- mattpocock-skills:tdd: tests only at spec interfaces (seams); expected values are literals, never recomputed the way the code computes them; one behaviour per test.
- engineering:testing-strategy: pyramid split (many unit/property, few integration, GPU/e2e only where the spec needs them), explicit gates.
- Not invoked (no effect on a test plan written from a spec): nvidia-skill-finder, code-review, debug, system-design, architecture, documentation, huggingface-best. Web tools not needed (no external facts claimed).

Conventions. REPO = repo root. Style follows tests/ (plain pytest, `P = experiments.phase2.vast.policy`, tmp_path, monkeypatch). `data` marker = needs data/fixtures. Facts marked [spec: Bn] or [code: file:line]. Baseline numbers are from the spec (1447 passed, 6 skipped, 31 deselected; 426 tracked files; 56 vulture findings).
Test ids are CT-n. "Top 5" tests have a pytest sketch at the end.

---------------------------------------------------------------------------------------------------------
## A. B1 bundle manifest: `code_files()`

CT-1 code_files drops every forbidden root [TOP5]
- level: unit. setup: tmp git repo; monkeypatch `P.ROOT`; create and `git add` (or leave untracked) one file per case: `.claude/worktrees/x/a.py`, `.claude/settings.json`, `archive/phase0/a.py`, `data/x.npz`, `docs/a.md`, `experiments/phase1/lab/a.py`, `experiments/phase2/results/r.json`, `experiments/phase3/lab/b.py`, `experiments/phase2/riga_x.csv`, `experiments/phase2/RIGA_X.JSON`, `experiments/phase2/Riga-cases.jsonl`, plus keepers `src/objdel/__init__.py`, `tests/oracles/reward_reference.py`, `experiments/phase1/reward_gates.py`, `experiments/phase2/train_rl.py`.
- expected: result equals exactly sorted([the 4 keepers]); no element starts with `.claude/ archive/ data/ docs/`; none contains `/lab/` or `/results/` directly under `experiments/<phase>/`; none matches `riga*.csv|json|jsonl` case-insensitively.
- catches: a regex anchored wrong (old regex had `^(data|docs)/` only; `.claude/` and `archive/` leak today: G6 from 2 and 10); case-insensitivity lost; `riga` rule dropped when the dir is nested.

CT-2 untracked, not-ignored files outside the roots are kept
- unit. setup: tmp git repo with `.gitignore` containing `ignored.py`; untracked `experiments/phase2/new_script.py`, tracked `src/a.py`, ignored `ignored.py`, plus `tests/oracles/o.py`.
- expected: `new_script.py` in result; `ignored.py` not in result.
- catches: switching from `ls-files -co --exclude-standard` to tracked-only (a new script then misses the box) or dropping `--exclude-standard`.

CT-3 boundary of the filters: look-alikes are kept (property/adversarial)
- unit. input paths (each created): `experiments/phase2/labs/x.py`, `experiments/phase2/lab.py`, `experiments/phase2/resultsx/y.py`, `experiments/phase2/vast/lab/z.py` (lab at depth 3, NOT directly under `experiments/<phase>/`), `src/objdel/archive_util.py`, `tests/data_helpers.py`, `docs_tools/x.py`, `dataset/x.py`, `src/riga_loader.py`, `experiments/p/riga.py`, `experiments/p/riga.csv.py`, `tests/oracles/o.py`.
- expected: all kept EXCEPT none dropped; i.e. result contains every listed path (`riga.py` and `.csv.py` do not match `riga*.csv|json|jsonl` because the extension differs). The `vast/lab/z.py` case is kept because the spec says "directly under experiments/<phase>/" (flag in the review if the implementer drops it: the spec text is the contract).
- catches: unanchored substring regexes (`lab` matching `labs`, `data` matching `dataset`), `riga` rule matching on stem only.

CT-4 paths with spaces and non-ASCII stay whole
- unit. setup: `experiments/my dir/a b.py`, `experiments/é/x.py`, `tests/oracles/o.py`. expected: all three returned verbatim, no quoting, no octal escapes (git `-z`).
- catches: dropping `-z` (git quotes non-ASCII as "\303\251"). Regression for existing test_code_files_keep_paths_with_spaces_whole (tests/test_vast.py:1189), which must still pass.

CT-5 missing tests/oracles stops with an error
- unit. setup: tmp repo with only `src/a.py`. expected: `P.die` fires (SystemExit or the module's die type; assert via `pytest.raises(SystemExit)` matching how `die` is raised today) with message containing `tests/oracles`. Also: repo whose only oracles file is under `archive/tests/oracles/o.py` (must die: archive is filtered first, so the check sees none); repo whose only oracle is `docs/tests/oracles/x.py` likewise.
- catches: the oracle check running before the filter (a filtered-out oracle satisfying it).

CT-6 empty and single-file repos
- unit. `git ls-files` returns nothing: dies with the oracles message, not IndexError. A repo with exactly one file `tests/oracles/o.py` returns `["tests/oracles/o.py"]` (the `[:-1]` trailing-NUL strip must not eat a real entry; also check output with no trailing NUL if git changes).

CT-7 real repo: gates G6 and keepers [integration, needs git only]
- setup: the real REPO, after CLEAN. expected: `code_files()` has no entry starting with `.claude/` or `archive/`; contains `tests/oracles/reward_reference.py`, `experiments/phase1/reward_gates.py`, `src/objdel/__init__.py`; every returned path exists on disk (`os.path.isfile`), including no deleted-but-tracked file (git ls-files lists tracked-deleted files: catches CLEAN's moved files before `git add`; decide per owner flow, mark xfail-free by filtering with isfile in the implementation or documenting).
- also `.gitignore` contains the line `/.claude/worktrees/`; assert `git check-ignore .claude/worktrees/x/y` exits 0 in a tmp repo seeded with the repo's `.gitignore`.
- catches: the 2 `.claude/` files and 10 `archive/` files shipping (spec G6); anchoring `/` dropped (would ignore any nested `.claude/worktrees`, harmless but wrong).

CT-8 existing callers keep passing
- regression. tests/test_p2_breakit_ops.py:298-299 (no riga table, no `experiments/phase0/results/vlm_*`), tests/test_vast.py:528-533, 576, 1082-1083, 1189-1198 pass unchanged. A mock-free assertion that `code_files.txt` written by `cmd_sync_up` equals `"\n".join(code_files())`.

## B. B1a bundle import (no GPU) [TOP5]

CT-9 copy of exactly the code_files() set imports and loads only copy modules
- level: integration, no GPU, no data. setup: `dest = tmp_path/"bundle"`; copy each path from `code_files()` (use REPO as ROOT) to `dest`, preserving structure; copy the built extension file(s): `objdel/_gpucheck*.so|.pyd` from `objdel._gpucheck.__file__` into `dest/src/objdel/` ONLY IF src layout; else locate per the editable install (the test must work out where `objdel` lives in the copy: spec says the copy plus the built extension). Run `subprocess.run([sys.executable, "-I", "-c", PROBE], cwd=dest, env={"PYTHONPATH": ...minimal...})` where PROBE imports: `objdel, objdel._gpucheck, experiments.phase1.reward_gates, tests.oracles.reward_harness, experiments.phase2.train_rl, experiments.phase2.g_probe, experiments.phase2.select_ckpt, experiments.phase2.final_eval`, then prints `json.dumps([m.__file__ for m in sys.modules.values() if getattr(m,"__name__","").split(".")[0] in {"objdel","experiments","tests"} and getattr(m,"__file__",None)])`.
- expected: return code 0; every printed file `Path(f).resolve().is_relative_to(dest.resolve())` (no module from REPO leaks through the editable install or `sys.path`); the list is non-empty and includes all 8 modules' files.
- caveat: `-I` ignores PYTHONPATH and cwd; the test must make `dest` importable some other way: write a `.pth` into a temp venv site dir, or `-c "import sys; sys.path[:0]=[...]"` inside PROBE (explicit insert is allowed: it runs inside the child). The editable install of REPO must be absent from the child: use `sys.flags.isolated` plus `-S`-free but strip the repo's `__editable__` finder by running the child with a venv without objdel, or assert after the fact via the file check above (the file check is the real gate).
- skip rules: skips only if trl/vllm-free imports fail in the parent for an optional dep; NOT allowed to skip because `.claude/` exists. If an import needs torch/trl and it is missing in CI, fail explicitly with a message (the spec says "a test checks this with no GPU").
- catches: a file that train_rl imports but `code_files()` drops (e.g. a helper under `experiments/phase2/lab/` or `docs/`), a test oracle that imports from `archive/`, a hidden dependency on `.claude/`, and silent import from REPO via editable install (the Vast box failure mode: imports fine locally, ImportError on the box).

CT-10 negative control for CT-9
- integration. Delete one file the bundle needs (`tests/oracles/reward_harness.py`) from the copy; the same probe must fail with returncode != 0. Catches a vacuous probe (for example one that imports from REPO and so always passes).

CT-11 IMPORT_CHECK covers the spec list
- unit. `P.IMPORT_CHECK` (policy.py:69) imports at least `objdel._gpucheck`, `experiments.phase1.reward_gates`, `tests.oracles.reward_harness` (today's value); and is executed under `python -I` in cmd_sync_up (assert the remote command contains `python -I -c`). Catches weakening the remote check while changing the manifest.

## C. B2 `--procs` removed, gc.freeze

CT-12 `--procs` is an argparse error [TOP5]
- level: unit (subprocess or in-process argparse). Run `python -I -m experiments.phase2.train_rl --procs 1 --model qwen3-vl-2b --ckpt x --train y --out z`. expected: returncode 2; stderr contains `unrecognized arguments: --procs` (argparse text; prefix-abbreviation hazard below).
- adversarial: `--proc 1` and `--pro 1` must ALSO fail with code 2 or be a documented abbreviation of another option. argparse `allow_abbrev=True` accepts a unique prefix: if another option begins with `--pro` (e.g. `--profile`) then `--pro` would silently bind to it. Test: for each prefix `--p, --pr, --pro, --proc` assert either exit 2 or that the parsed value lands in an option other than a pool size (inspect `ap.parse_known_args`). Catches an abbreviation silently swallowing `--procs`.
- `--help`: stdout does not contain the substring `procs` nor `scoring workers`. expected exit 0.
- catches: option left in as hidden (`help=argparse.SUPPRESS`), or kept as ignored.

CT-13 no process pool is created
- unit. monkeypatch `multiprocessing.get_context`, `multiprocessing.Pool`, `concurrent.futures.ProcessPoolExecutor.__init__`, `os.fork` with functions that raise AssertionError("pool"); run `train_rl.main([...minimal args...])` with `serve`, model loading and trainer stubbed (reuse fakes from tests/test_p2_train_rl.py / test_p2_scripts.py). expected: no assertion fires; also `hasattr(train_rl, "pool")` is False and source of `experiments/phase2/train_rl.py` contains neither `fork_pool` nor `pool_size`. (Source grep is a weak guard but cheap; the monkeypatch is the real test.)
- catches: a pool created but unused (leaks CPUs and forks before CUDA).

CT-14 gc.freeze exactly once, after the last serve.preload
- unit. record call order: monkeypatch `serve.preload` to append `("preload", args)`, `gc.freeze` (patch `train_rl.gc.freeze`, or `gc.freeze` itself) to append `("freeze",)`, stub rest of main. Run with `--train` of two datasets / main code path that calls preload k times (k=1 today at train_rl.py:49; also parametrise a config that makes it 2 if reachable). expected: `calls.count(("freeze",)) == 1`; the index of freeze > index of the last preload; and freeze is called before the trainer starts (stub `trainer.train` appends ("train",), assert freeze index < train index).
- failure paths: if main exits early with SystemExit before preload (bad args), freeze is not called (count == 0, no exception). If main is called twice in one process, each call freezes once (count 2): documents the contract.
- catches: freeze before preload (the preloaded objects then are not frozen and gen-2 GC cost remains), freeze twice, freeze in a try/finally that runs on error.

CT-15 other options and rl.json keys unchanged [TOP-adjacent regression]
- integration (slow marker as in test_full_rl_saves_resumes...). Update the old test to drop `"--procs", "1"`; assert `set(json.loads(rl.json))` equals a literal key set captured BEFORE CLEAN (capture it now from the current test run and paste as a literal set in the test; do not recompute from the code). Also `ap._actions` option strings minus `--procs` equals the literal pre-CLEAN list.
- catches: an accidental rename/removal of another option or output key during the edit.

## D. B3/B4/B5/B6 file layout

CT-16 removed paths hold no tracked file
- unit (needs git). `git ls-files <path>` empty for each B3 path (the 7 listed), AND `git ls-files experiments/phase0 experiments/phase1/lab experiments/phase2/lab | wc -l == 0` (G5, from 426). Untracked leftovers: also assert the folders do not exist on disk (a leftover untracked `__pycache__` is allowed: use `ls-files -co --exclude-standard`).
- catches: files removed from index but left; or moved copy left behind in both places.

CT-17 modules that must stay importable
- unit. `importlib.import_module` for `experiments.phase1.{build_data,build_targets,fetch_lidc_pixels,prepare_rl,eval_run,reward_gates}` and `experiments.phase1.bench.xarch_hashes`. expected: no exception, no `.lab.` module name in `sys.modules` afterwards (`not any(".lab." in n or n.endswith(".lab") for n in sys.modules)`).
- catches: a live module that imported a lab helper (the phase 1 scripts did `from lab...`).

CT-18 archive READMEs [TOP-adjacent]
- unit. For each of `phase0 phase1_lab phase1_scripts phase2_lab smolvlm` (+ `controller` unless owner keeps it: read the decision from a module constant or env; default requires it): dir exists; for EVERY direct subfolder `d` of `archive/` (use `iterdir`, not the fixed list; this also catches an extra folder with no README): `d/README.md` exists and contains the four headings at line start, in this order: `## What it was`, `## What replaced it`, `## Retiring commit`, `## How to restore`. The restore section (text after the last heading) contains `pre-clean-2026-11-08`. `archive/__init__.py` does not exist; no `__init__.py` anywhere makes `archive` a package root (`archive/*/__init__.py` allowed only if a subfolder keeps a package layout: assert `archive/__init__.py` only).
- edge: heading in a code fence or inside a sentence does not count (match with `re.M` and `^## What it was$`).
- C5 check: every test node id removed by CLEAN appears in the README of the archive folder holding it. Implement as: parse a list `tests/…` ids recorded in a test constant? No: spec gives no list. Test instead that any `archive/*/README.md` mentioning `tests/` ids uses the `path::name` form (regex) and that each cited path under tests/ no longer defines that name (`ast` scan of the cited file; if the file still exists, the function must be absent).
- catches: missing headings, wrong tag name (typo `pre-clean-2026-11-08`), archive turned into an importable package (pytest collection and ruff then pick it up).

CT-19 B5 hygiene scans (these run as tests/test_repo_hygiene.py, B14)
- unit, no data, no GPU. Implement with `ast`, not regex, for imports:
  a. for each `tests/**/*.py`: no `Import/ImportFrom` whose module starts with `archive` or contains a `lab` segment (`.lab.` or `.lab` at end); no `sys.path.insert/append` call whose argument source text contains `lab` or `archive`; no string constant (code, not docstring: skip `Expr(Constant)` first statements of module/class/function) containing `experiments/phase1/lab` or `experiments/phase2/lab`.
  b. for each `src/objdel/**/*.py`: no import of top-level `experiments` or `tests`.
  c. stale paths: scan live text files (git ls-files minus the exemptions `docs/changelog.md`, `docs/research/**`, `experiments/*/impl/**`, frozen phase `.md` records) for `experiments/phase0`, `experiments/phase1/lab`, `experiments/phase2/lab`.
  d. every backticked repo-relative path in `README.md` and `docs/**/*.md` (minus exemptions) exists: regex ``` `([\w./-]+/[\w./*-]+)` ```, keep those whose first segment is a top-level repo dir (src, tests, experiments, docs, cuda, archive, data excluded since git-ignored?). Decision needed on `data/...` paths (git-ignored, absent in CI): exempt `data/`, `build/`, `.venv/`, `.claude/`; glob patterns (`*`) checked with `glob`. Also trailing `:line` or `::name` stripped, trailing `/` dir ok.
- negative controls (important, they prove the scanner is not vacuous): run the scanner functions on a tmp tree seeded with one violation of each kind and assert each is reported; run on a clean tmp tree and assert nothing.
- adversarial: `from . import lab` and `import importlib; importlib.import_module("lab...")`: the second is not detectable by AST import scan; document as known limit (a runtime import check CT-17 covers it). A docstring mentioning `archive/...` must pass. `sys.path.insert(0, str(ROOT / "archive" / "x"))` must be reported.
- catches: G7 (3 importing and 2 reading lab files today).

CT-20 fixtures moved (B6)
- unit. `tests/fixtures/vast_snapshot.json` exists and parses as JSON; `tests/fixtures/topo_probe_v3e_seed*.json` glob has >= 1 match and each parses; no file under `tests/` or `src/` names `experiments/phase2/lab`; the Dockerfile digest check in tests/test_vast.py reads no path outside `tests/` (ast: the test function has a string constant matching `sha256:[0-9a-f]{64}` and no `open/read_text` on a path containing `lab`). Run test_vast.py and test_p2_* fully: still pass.
- catches: moved fixtures with old path hard-coded (FileNotFoundError hidden by `skipif exists`): scan for `skipif`/`exists()` guards around these fixtures and fail if the fixture path is conditional (a skip would hide the loss).

## E. B7 golden for target builder [TOP5: the no-mask golden gate]

CT-21 golden holds no masks and no per-case RIGA [TOP5]
- unit, no data (the file is tracked). setup: `z = np.load("tests/goldens/targets_lab.npz", allow_pickle=False)`. expected: for every key `z[k].ndim <= 1` and `z[k].dtype != bool` (also reject `object` dtype: pickle risk, load with allow_pickle=False must succeed); all LIDC keys match `^lidc/[^/]+/[^/]+$`; no key starts with `riga/` followed by a case id: allowed RIGA keys are exactly the aggregate set (assert `[k for k in z.files if k.startswith("riga")]` has exactly 3 entries: lab vertex total, case count, hash; names fixed in the test as literals once the implementer picks them; until then assert count == 3 and each array has size == 1 or (hash) a bytes/uint8 vector of 32).
- size guard: `Path(...).stat().st_size < 2_000_000` (a mask leak is megabytes). Also no key contains `mask`.
- C8 sibling: `git ls-files | grep -i riga` shows no `.npz/.npy/.png` mask and no per-case RIGA result.
- catches: a 2-D crop array sneaking in, bool dtype masks, per-case RIGA data (owner rule C8), pickled objects.

CT-22 generic golden rule (G15)
- unit. For every `tests/goldens/*.npz` except `polycore.npz`: no array with `ndim > 1` or bool dtype. Negative control: write a tmp dir with one npz having a `(2,2)` array and one with a bool vector and call the same checker function; both reported. Also verify polycore.npz is the only exemption by name (a renamed `polycore2.npz` must NOT be exempt).
- catches: exemption by prefix/glob; checker skipping when the dir is empty (assert `polycore.npz` exists so the glob is non-vacuous).

CT-23 test_fixtures_differential_vs_lab_golden [data marker]
- level: integration, `@pytest.mark.data`. In tests/test_targets.py. Steps: read masks from `fixtures_dir/{riga,lidc}_experts.npz` (same `step`s as the old test, 20 for riga, 6 for lidc, so the sample is the same cases as the recorded golden); crop with `tests.oracles.crops_oracle.pad_div`; check crop hash recorded in golden per case equals sha256 of the cropped mask bytes (fail with the case id); run production builder `check(m, names)`; assert (1) per LIDC case: production nodule vertex count == `z[f"lidc/{item}/{rater}"]` lab vertex count: exact, per case, report all mismatches; (2) fraction of LIDC cases with production text hash == recorded lab text hash >= 0.99 (computed over ALL cases, denominator printed; 3/1969 deviations today = 99.85 %, so the gate has room of 19 cases); (3) RIGA production vertex total <= recorded lab total.
- boundary: if a case key in the golden is missing from fixtures (or the reverse), the test FAILS, not skips (assert key sets equal; the comparison must not iterate over the intersection only). If fixtures are absent, standard `data` skip.
- the old test `test_fixtures_differential_vs_lab_builder` must not exist (ast scan of tests/test_targets.py); `sys.path` untouched (CT-19a).
- mutation checks (run once by hand, listed here as expected behaviour): (a) flip one hash byte in a copy of the golden -> test fails on step 1; (b) perturb one LIDC vertex count -> step (1) fails with that case id; (c) set RIGA lab total to production total - 1 -> step (3) fails.
- catches: a golden recorded from the wrong crop (hash check), a vacuous comparison over zero cases (assert `lidc >= 300` approximate size: choose from the golden: `len(lidc keys) > 0` and equals recorded case count constant), tolerance creep.

CT-24 crops_oracle
- unit (synthetic, no data). `pad_div` on a hand-built mask: single pixel at (0,0) in 512x512; mask touching the right border; empty mask (all False); full mask (all True); an off-centre 3x3 blob; non-square input. Expected values: hand-derived literals from the lab's rule (read `tests/oracles/crops_oracle.py` after it exists and derive the expected shape/offset independently on paper from the spec name `pad_div`: output side divisible by the divisor; crop contains the whole object; padded with False). Property: `out.sum() == mask.sum()` (no pixel lost) for random masks, output dtype bool, output H and W both divisible by the same divisor, deterministic (same input, same output; hash twice).
- empty mask: must either return a defined empty crop or raise a specific error: assert the behaviour matches the lab reference (record it as a literal), not "no exception".
- catches: pixel loss at borders (boundary-touching objects), off-by-one in padding, a port that changed the lab's behaviour.

## F. B7a golden for lab references

CT-25 p2_lab_refs.json and the 3 named tests [TOP-adjacent]
- unit. (a) the three test node ids exist (`pytest --collect-only -q` lists exactly these: `tests/test_p2_signal.py::test_live_b2_functions_and_group_scoring`, `::test_c2_stratified_bootstrap_equals_power_c2`, `tests/test_p2_probes.py::test_stratified_bootstrap_equals_power_c2`) and pass. (b) AST scan: none of the 3 imports a lab module or reads a path containing `lab`. (c) JSON parses; contains literal numbers (floats) not NaN (`json.loads(..., parse_constant=raise)`), and is generated from seeded inputs: check that the seed values are stored in the JSON. (d) tolerances unchanged: AST-extract `rel`/`abs` of `pytest.approx` / `assert_allclose` calls in these 3 tests and compare to literals captured from the pre-CLEAN versions (record them now from git HEAD: copy the numbers into the test plan before implementing). (e) mutation: change one number in the JSON by 1e-3 relative -> the matching test must fail (proves the golden is actually compared; run by hand).
- edge: seeded synthetic inputs must be independent of numpy version: ensure the inputs come from the JSON/`default_rng(seed)` not legacy global state (a changed numpy then changes inputs, not only outputs; catches silent drift).
- catches: golden regenerated from live code (tautology: live == golden by construction); weakened tolerances.

## G. B8 Controller retired

CT-26 callbacks and eval surface
- unit. `objdel.train.callbacks` has `TraceStep`; `not hasattr(callbacks, n)` for `Controller, audit_alarm, FIRED`; `objdel.eval` exports (`hasattr`) `hack_audit, gates, metrics, summarize, vision`, `not hasattr(objdel.eval, "alarms_for")`. `grep -rE "Controller|alarms_for|audit_alarm|\bFIRED\b" src experiments/phase2 tests` returns 0 lines outside archive (G10); the README of `archive/controller` may contain them.
- catches: dangling import (`from objdel.eval import alarms_for` in a script only run on the box).

CT-27 TRLReward has no `recorded`, values unchanged [TOP5]
- unit/regression. Build `TRLReward` as tests/test_train_rewards.py does (reuse its fixtures, a seeded synthetic group of 16 completions, mix of valid, empty string, malformed JSON, truncated). expected: `not hasattr(r, "recorded")` for a fresh instance AND after calls; reward list returned by `r(prompts, completions, ...)` equals a literal list of floats captured BEFORE CLEAN (paste exact values, compare with `==` or `math.isclose(rel_tol=0)`: bit-equal, because spec says "reward values are unchanged"). Also `objdel.reward.spec.fingerprint()` equals a literal string captured before CLEAN (G3). Memory: calling 1000 times does not grow the instance `__dict__` size (no hidden accumulator): `len(vars(r))` equal before/after.
- catches: removing `recorded` while changing the reward path (e.g. a lazy fetch that no longer caches), fingerprint drift via an edited spec constant.

CT-28 trl_reward back-reference stays
- unit. `funcs, *_ = r.funcs()` (or `r.funcs()[0]`); assert `len(funcs) > 0` and for every `f in funcs`: `f.trl_reward is r` (identity). Then simulate the FastGRPO prefetch: call the code path in `objdel.train.fast` that reads `.trl_reward` (find via grep `trl_reward` in src/objdel/train/fast.py; use test_reward_prefetch.py's harness) and assert one checker call per group (existing test_reward_one_checker_call.py stays green).
- negative control: a function copy without the attribute makes the prefetch path fall back/raise exactly as before (document literal).
- catches: the `recorded`-removal edit also stripping the attribute assignment (G16), which silently disables prefetch: slower, same rewards, so no reward test sees it.

CT-29 breakit wiring end to end without Controller
- integration. tests/test_breakit_wiring.py: the test exists (name captured now: keep id), runs store -> reward -> TRL -> eval -> gates, passes; assert it no longer references `Controller` (ast) and it still asserts at least: one reward per completion, gates verdict object present. Mutation by hand: break `gates` import -> test fails (proves the stage runs).
- catches: test deleted with the Controller (spec says must still exist).

## H. B9 SmolVLM retired

CT-30 configs surface
- unit. `set(configs.MODELS) == {"qwen3-vl-2b", "lfm2.5-vl-450m"}` in Wave 1 (and `{"qwen3-vl-2b"}` after Wave 3: parametrise on a wave constant); `"smolvlm-256m" not in MODELS and not in MAX_COMPLETION`; `set(MODELS) == set(MAX_COMPLETION)` (every model has a budget: a key in one dict only crashes `grpo_config` with KeyError).
- `processor_kwargs(m, d) == {"padding_side": "left"}` for every `m in MODELS`, `d in ("riga", "lidc")` (exact dict equality, so no stray `size` key). `processor_kwargs("smolvlm-256m", "lidc")` raises KeyError or ValueError (assert any Exception; do not pin the type).
- catches: `size` kwarg leaking to Qwen (changes image tokens, so changes every prompt).

CT-31 grpo_config use_vllm [TOP5]
- unit (needs trl, `importorskip`). For each model: `c = grpo_config(model, tmp)`; `c.use_vllm is True`. `grpo_config("qwen3-vl-2b", out, use_cpu=True, use_vllm=False).use_vllm is False`. Explicit `use_vllm=True` stays True. All the plan-D8 assertions of test_grpo_config_is_plan_d8 stay (loss_type dr_grpo, epsilon 0.2/0.28, beta 1e-3, 1536, groups) with the smolvlm line replaced by `c.use_vllm is True`. `budget` branch for lfm 1024 stays until Wave 3.
- adversarial: an override `use_vllm=None` or `0` (falsy, not False): spec says "unless overrides['use_vllm'] says otherwise": define: any explicit override wins; assert `grpo_config(..., use_vllm=0).use_vllm` is falsy. Other overrides (`num_generations=8`) still raise ValueError (unchanged).
- the overrides dict merge order: `|` overrides last (config.py:195), so no hidden second `use_vllm` key (TypeError duplicate keyword would show as a crash: covered by the call).
- catches: `use_vllm` hardcoded True ignoring override (the CPU smoke tests need False); leftover `model != ...` expression.

CT-32 no smolvlm text in code (G9)
- unit. For `src/**`, `experiments/phase2/*.py`, `tests/**/*.py` (git ls-files -co): read text, `re.search("smolvlm", text, re.I)` yields nothing. The test file itself must not contain the word: build the pattern as `"smol" + "vlm"`. `archive/smolvlm/**` is exempt. Also `HuggingFaceTB` string absent.
- catches: leftover branches (processor size rule, `_GRAD_CKPT`, comments, model hub id in test fixtures such as tests/probe_models.py).

## I. B10 boxes_meet [TOP5]

CT-33 boxes_meet truth table [TOP5]
- unit + property. `from objdel.geometry.box import boxes_meet`. Literals (box = (0,0,10,10)):
  - empty iterable -> False; empty tuple, empty list, empty generator, `iter(())` all False.
  - overlap interior (5,5,15,15) -> True; contained (2,2,3,3) -> True; containing (-5,-5,20,20) -> True; identical (0,0,10,10) -> True.
  - touching edges, half-open so NOT meeting: (10,0,20,10) False; (-10,0,0,10) False; (0,10,10,20) False; (0,-10,10,0) False; corner touch (10,10,20,20) False.
  - one pixel overlap (9,9,20,20) -> True; (9,0,20,10) True.
  - overlap on x but not y (5,20,15,30) False; and y but not x False.
  - degenerate other (empty extent) (5,5,5,5): a<x1 and x0<c and b<y1 and y0<d = 5<10 and 0<5 and 5<10 and 0<5 -> True (spec formula, even though zero area); (0,0,0,0): 0<10 and 0<0 False. These two pin the literal formula so an "improved" version that rejects zero-area boxes is caught.
  - degenerate query box (5,5,5,5) vs other (0,0,10,10): 0<5 and 5<10 and 0<5 and 5<10 -> True.
  - negative coordinates, large coordinates (2**31, 2**40: Python ints, no overflow; numpy int64 values and numpy scalars as inputs work and return a Python `bool`: `type(result) is bool`, not `np.bool_`).
  - `others` as list, tuple, set of tuples, generator, `dict.values()`; a generator is consumed at most once and the function short-circuits: pass a generator that raises after the first matching box: result True without raising (early exit; spec says "some ... has"). A generator whose first item does not match and which has a bad tuple of length 3 raises (ValueError from unpack): unchanged behaviour, assert raises.
- property (hypothesis if installed, else seeded loop of 20,000 random boxes in [-50,50]): `boxes_meet(b, [o]) == (not (o.x1 <= b.x0 or b.x1 <= o.x0 or o.y1 <= b.y0 or b.y1 <= o.y0))` (independent De Morgan form; ensure ranges include x0 > x1 inverted boxes: for inverted box formula must still equal the spec expression, so compare against the spec expression literally copied as a lambda in the test, in a form different from the implementation: write it with `max/min`: overlap iff `max(x0,a) < min(x1,c) and max(y0,b) < min(y1,d)` holds only for valid boxes: restrict the property to x0<x1,y0<y1,a<c,b<d and test inverted separately against the literal spec lambda). Symmetry: `boxes_meet(b,[o]) == boxes_meet(o,[b])` for valid boxes. Monotone: adding a box never turns True to False; `boxes_meet(b, A+B) == boxes_meet(b,A) or boxes_meet(b,B)`.
- catches: `<=` vs `<` (the half-open bug), an `all` for `any`, True on empty, a dropped axis.

CT-34 crowded wrappers equal boxes_meet [TOP-adjacent]
- unit. `objdel.data.lidc.crowded(s, box) == boxes_meet(box, s.others)` and `objdel.targets.frame.crowded(item, box) == boxes_meet(box, item.meta.get("other_nodules", ()))` over the CT-33 literals (use a tiny `types.SimpleNamespace(others=...)` for the slice, an `Item`-like with `meta` dict). Cases: `meta` has no `other_nodules` key -> False; `other_nodules` key present with `[]` -> False; `meta["other_nodules"]` is None (what does the old code do? `any(... in None)` raises TypeError: the spec keeps `get(..., ())`, so None raises; pin as a literal: `pytest.raises(TypeError)`), values given as lists `[a,b,c,d]` not tuples work.
- regression, data-backed variant [data]: over real store items (`store` fixture, first 2,000 items x 5 windows) the new `crowded` equals the OLD expression evaluated inline in the test (copy the one-line `any(a < x1 and x0 < c ...)`: an independent restatement, acceptable because it is the spec's formula) -> identical bool per call. Count of True equals a literal captured before CLEAN.
- AST check: neither `crowded` body contains a `Compare` with `<` chain any more (the rule lives in one place): count of nodes `Compare` with 2+ ops across `src/` for the pattern `a < x1 and x0 < c` equals 1 (only box.py).
- catches: wrappers drifting; a changed default for the missing meta key; wrapper returning a truthy non-bool.

## J. B11-B15 configuration and gates

CT-35 preflight default config (B11)
- unit. Parse Makefile (`make -n preflight` in REPO, `-n` prints, runs nothing): output contains `experiments/phase2/configs/r1_qwen3vl2b.yaml` and not `r1_lfm.yaml`/`rl_shifted.yaml`. `experiments/phase2/configs/r1_qwen3vl2b.yaml` exists; `load` it with the preflight module's own config loader (find via tests/test_preflight.py) without error; `rl_shifted.yaml` is NOT used as a default (steps-only run file: the loader on it fails or lacks model keys: assert whichever preflight does today as a literal). G12: with CONFIG override `make -n preflight CONFIG=x.yaml` uses x.yaml.
- catches: Makefile default pointing at a deleted LFM config (broken `make preflight` with no args).

CT-36 train_sft.load_items (B12)
- unit. `hasattr(train_sft, "load_items")` and `not hasattr(train_sft, "preload")`; `inspect.signature(load_items)` parameters equal (names, defaults) a literal list captured from the old `preload` signature before the rename; every caller updated: `grep -rn "train_sft.preload\|train_sft import preload" src experiments tests` returns nothing. Behaviour: one existing test that called `preload` (tests/test_p2_train.py or test_p2_scripts.py; find by grep) now calls `load_items` and passes with identical expectations. Edge: the name `preload` is also a different function `serve.preload` (train_rl): ensure that one is NOT renamed (`hasattr(serve, "preload")` stays True, train_rl.py:49 still works). 
- catches: rename that breaks the unrelated `serve.preload` through a global search-and-replace.

CT-37 dead code gate (B13, G8) [TOP-adjacent]
- unit/integration. `pyproject.toml` parsed with `tomllib`: `"vulture==2.16"` in `project["optional-dependencies"]["dev"]` (exact pin string); `[tool.vulture]` has `paths == ["src","experiments","tests"]`, `min_confidence == 60`, and `ignore_names`/`ignore_decorators` lists. Each of the two groups is preceded by a comment: read raw text, assert the line immediately above (or on the same line as) `ignore_names` and `ignore_decorators` starts with `#` (or a comment block is within 3 lines above).
- `subprocess.run([sys.executable,"-m","vulture"], cwd=ROOT, capture_output=True)`: returncode 0 and `stdout == "" and stderr == ""`. `make dead` equals the same (`make -s dead`, no echo of the command: assert Makefile recipe line starts with `@` or `make -s`, so "prints nothing" holds; check `stdout.strip() == ""`). Whitelist: `git ls-files | grep -i whitelist` and `*_whitelist*.py` is empty; no file named `whitelist.py`.
- negative control (run on a tmp copy): add `def _unused_zzz(): pass` to a file under src -> `make dead`/vulture exits 3 with that name in stdout. Proves the gate bites. If vulture is not installed, test FAILS (B14), not skips: `importlib.util.find_spec("vulture") is not None` asserted with a message.
- test fakes: AST scan tests/ for classes named `Fake*`/`_Fake*` with methods whose params were renamed to `_x`, `_unused`, `_` : assert param names do not start with `_` unless the real function has them (spec: "keep real names"). Weak heuristic: flag params named `_`, `_a`, `_arg`, `unused*` in tests.
- catches: gate passing because paths are empty (e.g. `paths` relative to cwd), whitelist crutch, vulture version drift.

CT-38 ruff exclude (B15) and lint
- unit. `tomllib.load(pyproject)["tool"]["ruff"]["extend-exclude"] == ["archive"]` (exact list). Negative control in a tmp dir: `archive/x/bad.py` with an unused import: `ruff check` there exits 0 with config, and a copy of the same file under `src/` exits non-zero. `make lint` exits 0 (G1).
- catches: archive code failing lint, or `extend-exclude` wiping the default excludes (it must be `extend-exclude`, not `exclude`: assert key name; `.venv` still excluded by default: `ruff check --show-files` lists nothing under `.venv`).

CT-39 hygiene file shape (B14)
- unit. `tests/test_repo_hygiene.py` exists; AST: it contains no `pytest.mark.data`, no `gpu` marker, no `skip`/`skipif`/`importorskip` calls (vulture missing must fail); `torch` not imported. It runs in < 30 s (performance gate: `pytest --durations` shows each test under 20 s; vulture itself ~ seconds). The "bundle" check is NOT in this file (`.claude/`, `archive/` strings absent from its code) and is in tests/test_vast.py.
- after Wave 3: `lfm` absent from `src/` (case-insensitive, word fragments `lfm2.5` included): parametrise by a WAVE constant defaulting to 1 so Wave 1 passes.

CT-40 checker topology_any (B16)
- unit. `objdel.checker.volume.topology_any` exists; behaves as before: capture before CLEAN (run now) 6 literal cases on tiny volumes: empty volume, single voxel, two voxels diagonal-touching (6- vs 26-connectivity answer as today), a hollow cube (cavity, Betti-2), a ring/torus (Betti-1), volume touching border. Pin the output literals (type and values) in the test; tests/test_volume.py stays green. `topology_any.__module__ == "objdel.checker.volume"`. Module docstring: `"reward" not in doc.lower() or "not" precedes it`: simplest exact check: the docstring does not contain the phrases `the reward uses`, `used by the reward`, `reward uses` (case-insensitive) and does not contain `reward` followed within 40 chars by `topology_any` unless negated; require instead a positive statement test: docstring contains `not` near `reward` is too fuzzy, so restrict to the negative phrases above and review by hand.
- catches: the cleanup removing or renaming the function, or its behaviour on 26-connectivity.

## K. B17 data moves, B18 ADR

CT-41 local data moves (B17) [local, marker `data`, needs data/]
- unit. Using `-I`. After the move: `data/views/lidc-bundle`, `data/views/lidc-bundle-v2`, `data/eval_views_archive/probe-precision-r1` do not exist; they exist under `data/_retired/2026-11-08/` (same relative layout or flat: pin the layout chosen); `data/views/lidc-bundle-v3` and `data/eval_views/lidc-val-rotated` still exist. `LIST.txt` has one line per moved folder with three fields (path, size in bytes, sha256 hex 64 chars); recompute: each size equals `du`-style recursive byte total (sum of `st_size`) and each sha256 equals hashlib of the moved `rows.jsonl` (independent recomputation: allowed since the golden is the file content). A folder with no `rows.jsonl`: LIST.txt shows a literal `-` / none (decide; test pins it); fail if the field is empty. No file lost: moved folder file count == pre-move count (record counts before moving in a sidecar and compare).
- safety: `data/` is untrusted: open files only with `-I` python and never execute anything inside.
- catches: a copy instead of a move (both places exist, space doubled), wrong checksum, deleted-not-moved (C2).

CT-42 ADR (B18)
- unit. `docs/architecture/decisions/` listing: the new file matches `^(\d{4})-library-keep-list\.md$`, its number is exactly max(other numbers)+1 (so no collision, no gap); the other ADR files' numbers are unique. It has every heading of the ADR template: read the template headings from an existing ADR in the same directory (the two most recent ones; intersect their `##` headings) and assert all are in the new file. Contains the three section words KEEP, DEFERRED, NOT ADOPTED (headings or table columns); DEFERRED entries each name an owning phase: regex per DEFERRED line `phase [A-Z]\d*` or a phase letter; at least one library in each section; `vulture`, `shapely` appear (spec B13, B20); `vastai` appears in NOT ADOPTED/dependency context if B19 is taken.
- Also docs link check (CT-19d) covers any backticked path in the ADR.
- catches: number collision, missing template headings, DEFERRED with no owner.

## L. B19/B20 optional

CT-43 vastai pin (only if B19 taken)
- unit. `tomllib`: no dependency starting `vastai` in any of `dependencies` and every `optional-dependencies` group; `pyyaml` in `dev`. `docs/runbooks/gpu-run.md` contains the pinned version string `1.11.0` and an install line using an isolated tool (`uv tool install vastai==1.11.0` or `pipx`: assert regex `(uv tool install|pipx install)\s+vastai==1\.11\.0`). Preflight: monkeypatch `subprocess.run` for `vastai --version` returning `vastai 1.11.0` -> check passes; returning `1.10.0`, `1.11.1`, empty output, FileNotFoundError (not installed), nonzero exit -> check FAILS with a message containing both versions (not an exception). `P.vast()` still calls subprocess with argv[0] == "vastai" (monkeypatch subprocess.run; assert no `import vastai`). Existing test_preflight.py and test_vast.py stay green.
- catches: unpinned/ignored version drift, a crash instead of a failed check on a missing binary, Pillow pin change (C10: `uv pip freeze` hash of project env unchanged: compare `pip list` text of pinned packages Pillow, cryptography to pre-CLEAN literals).

CT-44 shapely cross-check (only if B20 taken; tests/test_checker_shapely_crosscheck.py)
- property + integration. Dev extra pin `shapely==2.2.*` (string match). (a) defect class: for 2,000 seeded random polygons (self-touching, bow-ties, repeated vertex, spike, collinear, tiny-area, ring with hole as 1 ring) the checker's class, for the mapped classes only (mapping table lives in the test, with literal reason strings: "Self-intersection" -> self-intersection class, "Ring Self-intersection", "Too few points", ...), agrees with `shapely.is_valid_reason(Polygon(pts))` prefix. Unmapped classes are counted and the count printed; assert at least 5 of the 13 classes are mapped and each mapped class occurs >= 20 times in the sample (else vacuous). (b) raster: for legal polygons on a 64x64 grid, the checker mask equals `shapely.contains_xy(poly, xs+0.5, ys+0.5)` exactly (pixel centres); points exactly on an edge are a boundary case: use integer-vertex polygons where pixel centres at .5 can lie on diagonal edges (e.g. triangle (0,0),(8,0),(0,8): centre (3.5,4.5) lies on x+y=8): `contains_xy` is False on the boundary; the checker's rule must be pinned to a literal for that pixel (inclusion convention): if they differ the test must show the one pixel and name the convention rather than loosen to a tolerance. (c) `tests/oracles/validity_oracle.py` and `raster_oracle.py` still exist (file check) and are still imported by tests/test_checker_oracle.py.
- catches: oracles replaced by shapely (spec says they remain), vacuous class coverage.

## M. Wave 2 and 3 (guarded by a WAVE constant: xfail(strict=False) not allowed; use skipif WAVE < n with a visible reason)

CT-45 polygonize move (B21, G13)
- unit. `not Path("src/objdel/checker/polygonize.py").exists()`; `not hasattr(objdel.checker, "polygonize")`; `tests/oracles/polygonize_oracle.py` exists and `polygonize(mask, offset=(0,0))` works: single pixel at (3,4) with offset (10,20) -> literal ring (capture today from the live function BEFORE the move and paste: bit-for-bit vertex list and orientation); empty mask -> literal today's output (likely `[]`); full-image mask; mask with a hole (outer + hole ring); two components; mask touching all four borders; 1x1 mask; non-contiguous array and Fortran-order input same result as C-order (property); input dtype bool vs uint8 same result. Equivalence to the pre-move implementation: sha256 of `repr(polygonize(m))` over 200 seeded random masks equals a literal list captured before the move (G4 holds: target-set hashes equal).
- `certify.py` and `bitquad.py` absent from src (only after phase D device check); `objdel.checker` has no `bitquads`, `certify`; scan `src tests experiments` with ast for `from objdel.checker import ...bitquads|certify` and `objdel.checker.certify` / `.bitquad`: zero hits, including `experiments/phase1/bench/xarch_hashes.py` (CT-17 shows it imports) and `experiments/cuda_checker/`.
- catches: oracle copy drifting, an importer missed (xarch_hashes is exactly the kind: imports only at run time).

CT-46 Wave 3 (B22, B23)
- unit. `set(MODELS) == {"qwen3-vl-2b"}`; `grep -ri lfm` 0 hits in src, experiments/phase2/*.py, tests (build the word as `"l"+"fm"`); `param_groups` absent from `src/objdel/train/*.py`; no float32 `dtype` branch (AST: no `torch.float32` compare in configs/sft code); `fork_pool`, `pool_size` absent; `--procs` absent from `select_ckpt`, `g_probe`, `final_eval` (argparse exit 2 for each, like CT-12) and `--help` free of `procs`. Slow CPU smoke: tiny random Qwen3-VL, `max_completion_length` small (<= 32), finishes in < 120 s on CPU, saves a checkpoint, and loss is finite.

## N. Cross-cutting gates

CT-47 acceptance gates in one runner (marker: ci)
- integration. G1: `make lint` exit 0. G2: `make test` -> parse the pytest summary line: `0 failed`; passed == 1447 - R + N where R = removed tests (list kept in a constant by node id, each must be absent now) and N = new tests (constant; collected count equality via `pytest --collect-only -q` count compare). Exact equality (not >=) so a silently lost test fails. Skipped count <= 6 + new data/gpu tests (explicit), no new skip in hygiene. G3 fingerprint literal. G4: `target_set_sha256` of manifests unchanged: literal captured before (tests/test_p2_views.py). C1: `git diff --stat pre-clean-2026-11-08 -- cuda/` empty (read-only git diff; skip if the tag is absent locally with an explicit message rather than pass).
- catches: tests deleted to make the suite green.

CT-48 no hashed result changed (regression of today's numbers)
- GPU-free. Re-run, with `-I`, today's pinned numbers: fingerprint literal; target-set hash literal; `polycore.npz` golden comparison (tests/test_geometry*.py stay green); reward oracle agreement (`tests/oracles/reward_reference.py` vs live `reward` over the seeded group of CT-27, exact equality); cuda checker CPU path (`test_cuda_checker.py`) green. Numbers to hard-code: test IoU 0.645 / F 0.816 / P 0.695 / R 0.748 are GPU-run results [measured: README/memory], not reproducible here: not tested; only the deterministic hashes are.
- catches: a refactor of `crowded`/targets changing even one case's crop window (hash flips).

CT-49 CPU vs GPU bit-equality (GPU marker; must not run in CLEAN CI; only if the owner later runs on a GPU, C6 forbids spending now)
- GPU. `tests/test_cuda_checker_gpu.py` and `test_gpu_scoring.py`: unchanged by CLEAN, so the test is "files untouched": `git diff --stat pre-clean-2026-11-08 -- tests/test_cuda_checker_gpu.py tests/test_gpu_scoring.py cuda/` shows 0 changes, OR the same tests pass unchanged on an RTX box (not executed in this phase). Rationale: CLEAN must not change CPU-vs-GPU equality; the cheapest guard is that neither the checker source nor those tests moved. Also the CPU build of the extension is rebuilt from unchanged `cuda/` (stamp unchanged: `build/.gpucheck-stamp` not newer than sources).
- catches: moved/deleted GPU test files while claiming no change.

CT-50 performance gates
- hygiene suite total wall time < 60 s; `code_files()` < 5 s on the real repo (one subprocess); `make dead` < 60 s; `import experiments.phase2.train_rl` time unchanged ±20 % (no pool creation at import): measure `python -X importtime -I -c "import experiments.phase2.train_rl"` before and after and record; gc.freeze present means no per-step regression: not measurable on CPU, state as GPU-run item (R-step time 23.6 s/step with shared-prefix [measured: memory], untouched by CLEAN).
- bundle size: total bytes of `code_files()` after CLEAN <= the bytes before (never grows), and the 426 removed lab files plus 12 `.claude/archive` entries are gone (count `len(code_files())` drop >= 12). Catches: bundle growth and upload time (V0: slow uploads cost real rental money).

## O. Failure modes list (each has a test above)
- bundle imports locally but not on box: CT-9/10.
- `--procs` silently accepted by abbreviation: CT-12.
- freeze order wrong: CT-14.
- trl_reward attr dropped, prefetch silently off: CT-28.
- reward drift: CT-27, CT-48.
- golden holds masks / RIGA per-case: CT-21/22.
- golden regenerated from live code: CT-25 (mutation) and CT-23 (hash check on crops).
- half-open box bug: CT-33.
- archive importable / linted: CT-18, CT-38.
- missing test hidden by deletion: CT-47.
- gate passes because empty path / skipped: CT-37, CT-39 (no skip), CT-23 (key-set equality).

---------------------------------------------------------------------------------------------------------
## P. Pytest sketches (5 most important)

Ranking by bug-catching value and cost: (1) CT-9 bundle import (the box failure mode), (2) CT-1 code_files filters, (3) CT-33 boxes_meet (plus wrappers CT-34), (4) CT-27/28 reward and trl_reward, (5) CT-21/22 golden no-mask gate. CT-12/14 follow closely (sketched as bonus).

```python
# tests/test_vast.py  (CT-1)
import subprocess
import pytest
from experiments.phase2.vast import policy as P

def _repo(tmp_path, tracked, untracked=()):
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    for rel in (*tracked, *untracked):
        f = tmp_path / rel
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text("x")
    subprocess.run(["git", "add", "--", *tracked], cwd=tmp_path, check=True)
    return tmp_path

KEEP = ["src/objdel/__init__.py", "tests/oracles/reward_reference.py",
        "experiments/phase1/reward_gates.py", "experiments/phase2/train_rl.py"]
DROP = [".claude/worktrees/x/a.py", ".claude/settings.json", "archive/phase0/a.py", "data/x.npz", "docs/a.md",
        "experiments/phase1/lab/a.py", "experiments/phase2/results/r.json", "experiments/phase2/riga_x.csv",
        "experiments/phase2/RIGA_X.JSON", "experiments/phase2/Riga-cases.jsonl"]

def test_code_files_drops_forbidden_roots_and_keeps_the_rest(tmp_path, monkeypatch):
    root = _repo(tmp_path, tracked=KEEP + DROP[:5], untracked=DROP[5:] + ["experiments/phase2/new_script.py"])
    monkeypatch.setattr(P, "ROOT", root)
    got = P.code_files()
    assert sorted(got) == sorted(KEEP + ["experiments/phase2/new_script.py"])  # untracked kept, all DROP gone

@pytest.mark.parametrize("keep", ["experiments/phase2/labs/x.py", "experiments/phase2/vast/lab/z.py",
                                  "dataset/x.py", "docs_tools/x.py", "src/riga_loader.py", "experiments/p/riga.py"])
def test_code_files_filters_are_anchored(tmp_path, monkeypatch, keep):
    monkeypatch.setattr(P, "ROOT", _repo(tmp_path, [keep, "tests/oracles/o.py"]))
    assert keep in P.code_files()

def test_code_files_dies_without_oracles_even_if_one_is_archived(tmp_path, monkeypatch):
    monkeypatch.setattr(P, "ROOT", _repo(tmp_path, ["src/a.py", "archive/tests/oracles/o.py"]))
    with pytest.raises(SystemExit):  # adapt to how P.die raises
        P.code_files()
```

```python
# tests/test_vast.py  (CT-9, CT-10: bundle import, no GPU)
import json, shutil, subprocess, sys
from pathlib import Path
import pytest
import objdel._gpucheck as gpucheck
from experiments.phase2.vast import policy as P

MODS = ["objdel", "objdel._gpucheck", "experiments.phase1.reward_gates", "tests.oracles.reward_harness",
        "experiments.phase2.train_rl", "experiments.phase2.g_probe", "experiments.phase2.select_ckpt",
        "experiments.phase2.final_eval"]
PROBE = """
import sys, json
sys.path[:0] = {paths!r}                      # the copy only; -I ignores PYTHONPATH and cwd
import importlib
for m in {mods!r}: importlib.import_module(m)
print(json.dumps([getattr(m, "__file__", None) for n, m in list(sys.modules.items())
                  if n.split(".")[0] in ("objdel", "experiments", "tests") and getattr(m, "__file__", None)]))
"""

def _bundle(tmp_path):
    dest = tmp_path / "bundle"
    for rel in P.code_files():
        (dest / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(P.ROOT / rel, dest / rel)
    ext = Path(gpucheck.__file__)
    shutil.copy2(ext, dest / "src" / "objdel" / ext.name)   # the built extension goes beside the copied package
    return dest

def _run(dest, mods=MODS):
    code = PROBE.format(paths=[str(dest / "src"), str(dest)], mods=mods)
    return subprocess.run([sys.executable, "-I", "-c", code], cwd=dest, capture_output=True, text=True)

def test_bundle_copy_imports_and_loads_only_its_own_modules(tmp_path):
    dest = _bundle(tmp_path)
    r = _run(dest)
    assert r.returncode == 0, r.stderr[-3000:]
    files = [Path(f).resolve() for f in json.loads(r.stdout)]
    assert files and all(f.is_relative_to(dest.resolve()) for f in files), [str(f) for f in files if not f.is_relative_to(dest.resolve())]

def test_bundle_probe_is_not_vacuous(tmp_path):
    dest = _bundle(tmp_path)
    (dest / "tests/oracles/reward_harness.py").unlink()
    assert _run(dest).returncode != 0
```
(The editable install of REPO in the interpreter's site-packages stays importable inside the child: the sys.path prepend wins, and the file-location assertion is the real gate. A stricter variant runs the child with a venv that lacks objdel.)

```python
# tests/test_geometry_box.py  (CT-33, CT-34)
import itertools, random
import pytest
from types import SimpleNamespace
from objdel.geometry.box import boxes_meet
from objdel.data import lidc
from objdel.targets import frame

B = (0, 0, 10, 10)

@pytest.mark.parametrize("other,want", [
    ((5, 5, 15, 15), True), ((2, 2, 3, 3), True), ((-5, -5, 20, 20), True), (B, True), ((9, 9, 20, 20), True),
    ((10, 0, 20, 10), False), ((-10, 0, 0, 10), False), ((0, 10, 10, 20), False), ((0, -10, 10, 0), False),
    ((10, 10, 20, 20), False), ((5, 20, 15, 30), False), ((20, 5, 30, 15), False),
    ((5, 5, 5, 5), True), ((0, 0, 0, 0), False),     # zero-size boxes: the literal spec formula
    ((2**40, 2**40, 2**41, 2**41), False),
])
def test_boxes_meet_half_open_literals(other, want):
    got = boxes_meet(B, [other])
    assert got is want and type(got) is bool

def test_boxes_meet_empty_and_lazy():
    for empty in ([], (), iter(()), (x for x in ())):
        assert boxes_meet(B, empty) is False
    def gen():
        yield (5, 5, 15, 15)
        raise AssertionError("not short-circuited")
    assert boxes_meet(B, gen()) is True

def test_boxes_meet_matches_spec_formula_randomly_and_is_symmetric():
    rng = random.Random(0)
    spec = lambda b, o: o[0] < b[2] and b[0] < o[2] and o[1] < b[3] and b[1] < o[3]
    def box():
        x, y = rng.randint(-50, 50), rng.randint(-50, 50)
        return (x, y, x + rng.randint(1, 30), y + rng.randint(1, 30))
    for _ in range(20000):
        b, o = box(), box()
        assert boxes_meet(b, [o]) == spec(b, o) == boxes_meet(o, [b])

@pytest.mark.parametrize("others", [[], [(5, 5, 15, 15)], [(10, 0, 20, 10)], [(10, 0, 20, 10), (1, 1, 2, 2)]])
def test_crowded_wrappers_equal_boxes_meet(others):
    s = SimpleNamespace(others=others)
    item = SimpleNamespace(meta={"other_nodules": others})
    assert lidc.crowded(s, B) == frame.crowded(item, B) == boxes_meet(B, others)
    assert frame.crowded(SimpleNamespace(meta={}), B) is False
```

```python
# tests/test_train_rewards.py  (CT-27, CT-28)
GOLDEN_REWARDS = [...]   # literal floats captured from the pre-CLEAN code on the seeded group; paste here
GOLDEN_FINGERPRINT = "..."  # literal, captured before CLEAN

def test_trlreward_has_no_recorded_and_values_are_unchanged(make_trl_reward, seeded_group):
    r = make_trl_reward()
    out = r(**seeded_group)                  # same call shape as the neighbouring tests
    assert out == GOLDEN_REWARDS             # exact equality: bit-for-bit
    assert not hasattr(r, "recorded")
    n = len(vars(r)); r(**seeded_group); assert len(vars(r)) == n   # nothing accumulates

def test_every_reward_function_points_back_to_its_trlreward(make_trl_reward):
    r = make_trl_reward()
    funcs = r.funcs()[0]
    assert funcs and all(f.trl_reward is r for f in funcs)   # FastGRPO prefetch reads this

def test_fingerprint_unchanged():
    from objdel.reward import spec
    assert spec.fingerprint() == GOLDEN_FINGERPRINT
```

```python
# tests/test_repo_hygiene.py  (CT-21, CT-22: golden gate; same file hosts CT-19, CT-37, CT-32)
from pathlib import Path
import numpy as np
import pytest

GOLD = Path(__file__).parent / "goldens"
EXEMPT = {"polycore.npz"}

def bad_arrays(path):
    z = np.load(path, allow_pickle=False)        # pickled objects: ValueError, a test failure
    return [k for k in z.files if z[k].ndim > 1 or z[k].dtype == bool]

def test_goldens_hold_no_masks():
    files = sorted(GOLD.glob("*.npz"))
    assert "polycore.npz" in {p.name for p in files} and "targets_lab.npz" in {p.name for p in files}  # non-vacuous
    assert {p.name: bad_arrays(p) for p in files if p.name not in EXEMPT} == \
           {p.name: [] for p in files if p.name not in EXEMPT}

def test_targets_lab_golden_has_only_aggregate_riga_values():
    z = np.load(GOLD / "targets_lab.npz", allow_pickle=False)
    riga = [k for k in z.files if k.lower().startswith("riga")]
    assert len(riga) == 3 and all(z[k].size <= 32 for k in riga)    # vertex total, case count, hash
    assert all(__import__("re").fullmatch(r"lidc/[^/]+/[^/]+", k) for k in z.files if k not in riga)
    assert (GOLD / "targets_lab.npz").stat().st_size < 2_000_000

def test_golden_checker_catches_a_mask(tmp_path):
    np.savez(tmp_path / "a.npz", m=np.zeros((2, 2), bool))
    np.savez(tmp_path / "b.npz", v=np.zeros(4, bool))
    assert bad_arrays(tmp_path / "a.npz") == ["m"] and bad_arrays(tmp_path / "b.npz") == ["v"]
```

Bonus sketch (CT-12, CT-14):
```python
def test_train_rl_rejects_procs_and_help_hides_it():
    base = [sys.executable, "-I", "-m", "experiments.phase2.train_rl"]
    for flag in ("--procs", "--proc", "--pro"):
        r = subprocess.run([*base, flag, "1"], capture_output=True, text=True)
        assert r.returncode == 2, (flag, r.stderr[-300:])
    h = subprocess.run([*base, "--help"], capture_output=True, text=True)
    assert h.returncode == 0 and "procs" not in h.stdout

def test_gc_freeze_once_after_last_preload(monkeypatch, train_rl_stubbed):   # fixture stubs model, trainer
    calls = []
    monkeypatch.setattr(train_rl.serve, "preload", lambda *a, **k: calls.append("preload"))
    monkeypatch.setattr(train_rl.gc, "freeze", lambda: calls.append("freeze"))
    train_rl.main([...minimal args...])
    assert calls.count("freeze") == 1 and calls.index("freeze") > len(calls) - 1 - calls[::-1].index("preload")
```

## Q. Open points for the implementer (flag in review, not guessed here)
1. `P.die` exception type (CT-5): read it from policy.py.
2. Spec says `/lab/` "directly under experiments/<phase>/"; CT-3 pins that deeper `lab/` dirs are kept. If the owner wants deeper dirs dropped, change the spec first.
3. Where the built extension sits in the copied bundle (CT-9): depends on src layout; the test locates it from `objdel._gpucheck.__file__`.
4. `data/` paths in docs are git-ignored: CT-19d exempts them (decision to confirm).
5. Pre-CLEAN literals (rewards, fingerprint, topology_any, polygonize hashes, signatures, tolerances) must be captured BEFORE the code changes, from git HEAD; this plan does not invent them.
