# Final review, round 3: architecture, code organization, layering, APIs, DOC, CLEAN

Reviewer role: principal software architect for production ML systems.
Document: `REPO/docs/implementation2.0/Implementation.md` (14,506 lines, round 3 text).
Earlier files read: `reviews/final_r1_arch.md`, `reviews/final_r2_arch.md`. All 20 round-2 items are fixed in the
text (checked by grep and by reading the changed tasks): ADR numbers in B and C (only DOC lines name 0012), one
archive folder `archive/phase2_runs/configs/`, `b1_rows` everywhere in E, Vast `home()`, `StageBase.soft`,
`RlStage.rows_*`, the E8.1 name table, CLEAN8 on post-E names, DOC4's `specs/*.yaml` glob, the DOC9 E row, the B1.3
no-transformers precondition. None of them is repeated below.
Clean context kept: no transcripts, no memory folder, no `REPO/.claude/`, no git write, no rental, no REPO edit. One
probe script was written in the scratchpad (`reviews/r3arch_probe/field_order.py`).

## Skills used

| Skill | What it changed in this review |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | No NVIDIA catalog skill covers package layering, run specs or docs drift. No install, no recommendation. The kernel build home (E.5) is unchanged since round 2; no item. |
| `engineering:code-review` | Every claim was checked in the working tree: `preflight.py:26,49-69`, `policy.py:300-340,420-456`, `train_sft.py:70,196,225-233,297`, `hf_push.py:60,75`, `serve.py` imports, `dataset.py:86-160`, `eval/__init__.py:11-13`, the three RL plans, both launchd plists. This found items 1, 4, 6 and 9. |
| `engineering:debug` | I traced what fails at run time, not at import: the V1 preflight `TypeError` (item 1), a code hash over 0 files after E2.5 (item 4), a licence file that is not found after E8 (item 9), and a list rename that breaks lint once the ignores go (item 5). |
| `engineering:system-design` | I followed V1 from launch to the first `train_rl` step. The rows root has no way onto the box, and the manifest SHA has no way into `train_rl` (item 2). |
| `engineering:architecture` | Layer order per commit: the interim contract must reorder `verify` below `data` at E3b (item 5). The `truncated` rule has the wrong home (item 13). |
| `engineering:documentation` | One rule, one test: CLEAN12 and DOC4 both test doc paths, with different scopes (item 11). DOC accepts red tests until DOC13 (item 3). V4 has two prices (item 16). |
| `huggingface-skills:huggingface-best` | No leaderboard query: model choice is outside this review. I checked only that the HF parts touched by E (rows root on `svenk029/objdel-sft`, the `hf_push` card licence) still work after the moves (items 2, 9). |
| `system-design-skills:system-design` | The phase order is the main hypothesis. I checked each E commit against the gates it must pass (R-3: `make test` and `make lint` green after every commit). Items 3 and 7 show two places where the plan makes them red by design. |
| `system-design-skills:service-decomposition` | A modular monolith with enforced layers is still right; no service split is asked for. Item 4 puts provenance (`repo_root`, `code_hashes`) in `core`, so every layer can use one function without an upward import. |
| `mattpocock-skills:codebase-design` | Interface check: `to_argv` is the seam between a spec and a box chain, and its interface must fix the flag order (item 8). Package `__init__` re-exports are part of each package's interface and have no move rows (item 6). |
| Web (WebFetch) and a local probe | import-linter layers page: a listed layer that does not exist breaks the contract; `exhaustive` needs `containers` [web: import-linter.readthedocs.io/en/stable/contract_types/layers/]. pydantic field order with the `TrainerFlags` mixin, probed with pydantic 2.11.7: `['lr', 'shared_prefix', 'grpo', 'kind', 'name', 'soft', 'model', 'ckpt', 'steps']` [measured: reviews/r3arch_probe/field_order.py] (item 8). ruff 0.17.1 docstring count on the files E moves into `src/` (item 7). |

## Counts

blocking 3, major 8, minor 6 (17 items).

## Items

### 1. blocking | 00_front (also B, A) | V1 cannot get a preflight stamp: both configs in front 6.1 step 4 are `steps:` plans

Problem. Front 6.1 step 4 and order rule 2 require `make preflight CONFIG=experiments/phase2/configs/b1_rows.yaml`
("the V1 launch config") and `make preflight CONFIG=experiments/phase2/configs/a3_measure.yaml`. Front 6.2 says "the
launch config is `b1_rows.yaml`". `preflight.load` builds `train_sft.RunArgs(**raw)` [code:
experiments/phase2/preflight.py:49-56]. A `steps:` plan raises `TypeError: ... unexpected keyword argument 'steps'`
(measured by CLEAN: `CLEAN.0` debug row, line 12158). `b1_rows.yaml` is a `steps:` plan of `train_rl` argv (B1.5,
line 3616: `--ckpt rl-full-3ep`, `--train ...`, `--steps ...`). `a3_measure.yaml` is one too (A3, line 1161). A3
already says that a plan file cannot get a stamp (line 1191) and points back to "the V1 launch config (front 6.1
step 4)", which is `b1_rows.yaml`. So the pointer is circular. `cmd_launch` refuses without a stamp of its
`--config` [code: experiments/phase2/vast/policy.py:313-316]. CLEAN flagged this (lines 12749-12753, 14191) and left
it to "B and front". Neither changed. As written, V1 cannot launch.

Fix. Front 6.1 step 4: "On the tag `v1-b1`: `make preflight CONFIG=experiments/phase2/configs/b1_launch.yaml`
(the V1 launch config; a `RunArgs` file like `r1_qwen3vl2b.yaml`). `b1_rows.yaml` and `a3_measure.yaml` are plan
files: `run.sh` runs them after `launch.sh --config b1_launch.yaml`. A plan file never follows `preflight CONFIG=`."
Order rule 2 and front 6.2: the same change. B1.5 creates `b1_launch.yaml` (fields `model: qwen3-vl-2b`,
`output_dir`, `data_root`, `view` = the local B1 rows root, no `branches`, so `chain()` gives only a train stage that
is never run). B1.3's preflight wiring adds the SHA-256 of the rows-root `manifest.json` to the stamp's data hash
(`preflight.hashes`), so an edited rows root refuses the launch. B.9 adds a test: `preflight.load(b1_launch.yaml)`
passes; `preflight.load(b1_rows.yaml)` raises the pinned `TypeError`. 99_back decision 2 names `b1_launch.yaml`.

### 2. major | B (also 00_front) | The rows root has no way onto the box, and its manifest SHA has no way into `train_rl`

Problem. B1.4 step 5 uploads the rows root to `svenk029/objdel-sft` under `rows/lidc-train-b1/` and pins
`rows_revision` and `rows_manifest_sha256` "in `b1_rows.yaml`" (line 3605). Front 6.2 item 1 counts a "rows-root
download from HF (pinned revision)". But:
- No task downloads it. `fetch_hf.sh` / `cmd_fetch_hf` pulls a checkpoint revision into `runs/<id>/ckpt`
  [code: policy.py:528-536]. No B, A or CLEAN step extends it, and no `b1_rows.yaml` step does it. The HF token
  reaches the box only on the stdin of an ssh command from the Mac [code: policy.py:524, 536], so a plan step on the
  box cannot fetch a private repo.
- `train_rl` gets only argv. B1.5 says "Assert SHA-256 ... equals `rows_manifest_sha256` from the config" (line
  3625), but no `train_rl` flag carries the value (grep: no `--rows-manifest-sha256` or `--rows-revision` in the
  doc). E2.1 says "`to_argv` emits both as the `train_rl` flags that B1.5 reads; B's flag spelling wins", but B has
  no spelling.
B E08 (the dry run) builds the config with a stub trainer and does not cover either path.

Fix. B1.5:
(a) `train_rl` flags `--rows-revision HEX40` and `--rows-manifest-sha256 HEX64`; the SHA check reads the flag;
`b1_rows.yaml` passes both. E2.1's `to_argv` uses these spellings.
(b) `cmd_fetch_hf --rows REVISION ROOT` (and `fetch_hf.sh <run> --rows <rows_revision> lidc-train-b1`): it streams
the token as today and runs `hf_pull --subdir rows/lidc-train-b1 --out /workspace/data/lidc/rows/lidc-train-b1`
with `snapshot_download(..., allow_patterns="rows/lidc-train-b1/**", revision=...)`.
(c) The V1 launch line in front 6.2 note and in A3's header becomes: `sync_up.sh <run> && fetch_hf.sh <run>
rl-full-3ep && fetch_hf.sh <run> --rows <rows_revision> lidc-train-b1 && run.sh <run> .../b1_rows.yaml`.
(d) Tests: U61 also passes the flag; a new test with the Hub mocked checks the `allow_patterns`, the revision and
that the token goes only to stdin; E08 reads the rows root from the same path layout.

### 3. blocking | DOC (also E, CLEAN) | DOC accepts a red `make test` from the first E move until DOC13

Problem. DOC.7 says: "Phase E moves `experiments/phase2` | path and link checks fail in many docs | Expected; DOC13
fixes them in one pass" (line 11576). `tests/test_docs.py` runs inside `make test` (DOC.3.2 step 3). So:
- Check 2 (`test_backtick_paths_exist`) reads every backtick `src/`, `experiments/`, `tests/` path in `docs/**`
  outside frozen and plan files. The engineering pages name `src/objdel/train/fast.py`, `checker/...`,
  `eval/...` paths. Every E3/E5/E8 move commit breaks it.
- Check 6 fails on the E6 commit unless `docs/reference/configs.md` names each new `specs/*.yaml`; check 2 fails on
  the 12 archived plan paths.
- Check 8, U18 and U19 import `POLICY_MODULE = "experiments.phase2.vast.policy"` (line 11198). E8.1 deletes that
  module (its names move to `infra/vast/*`). U18 says a missing module is "a clear error, never a silent pass", so
  the suite fails. DOC13 step 3 changes the constant only after CLEAN wave 3.
This breaks the repo rule "`make test` and `make lint` must pass before a commit" (CLAUDE.md), E's R-3 and E-G1.
CLEAN archive commits (CLEAN5, CLEAN8) have the same exposure.

Fix. Add to E R-1: "A move commit also rewrites, with the same old -> new path map, every backtick repo path in the
DOC4 `doc_files` set (non-frozen, non-plan docs, `README.md`, `experiments/README.md`). The E6 switch commit adds
each `specs/*.yaml` to `docs/reference/configs.md`. The E8.1 move commit sets `POLICY_MODULE =
"objdel.infra.vast.ledger"` in `tests/docs_check.py` (a module string, R-1)." Add the same rule to the CLEAN task
rules for archive commits. Change the DOC.7 row to "Expected: each E or CLEAN commit fixes its own doc paths (E R-1);
DOC13 only re-runs the checks and rewrites `package.md`". Change DOC13 step 3 from "Set `POLICY_MODULE`" to "Check
that `POLICY_MODULE` imports (X9)".

### 4. blocking | E | E2.5 moves `code_hashes` into `sft/train.py`; the preflight code hash then covers 0 files

Problem. `train_sft.code_hashes()` globs `CODE = ("src/objdel/**/*.py", "cuda/*.cpp", ..., "experiments/phase2/*.py",
...)` under `ROOT = Path(__file__).resolve().parents[2]` [code: experiments/phase2/train_sft.py:70,225-233].
`preflight.hashes()` uses it for the stamp's `code` hash [code: preflight.py:60,69], and `refusal()` refuses a
launch when that hash changes. Run manifests on the box use it too [code: train_sft.py:196,297]. E2.5 moves "the
library part of `train_sft.py`" into `src/objdel/sft/train.py` (line 9518) and names no home for `code_hashes`. In
`src/objdel/sft/train.py`, `parents[2]` is `src/`. `ROOT.glob("src/objdel/**/*.py")` then matches nothing, the hash
is a constant, and the launch gate no longer sees code edits. Nothing catches it:
- E-G15's AST test looks for `parents[` only under `src/objdel/infra/` (E8.1, U34).
- `infra.paths.repo_root()` sits in `infra`, which `sft` and `rl` must not import (layer rule).
- `repo_root()` calls `git rev-parse`. The box has no `.git`: `code_files()` ships tracked files only [code:
  policy.py:377], and the box installs with `pip -e .` [code: policy.py:399].
`rl.run` (E6.5: `run_spec.json` with the git sha) needs the same function on the box.

Fix. Create `src/objdel/core/provenance.py` (stdlib only, so `core` still imports nothing from objdel):
`repo_root()` = env `OBJDEL_REPO`, else the first parent of `Path(objdel.__file__)` that holds `pyproject.toml` and
`src/objdel` (works on the Mac and on the box with the editable install), else raise; `CODE`; `code_hashes()`
(git sha `None` when there is no `.git`, as today). Move them there in E2 (a move commit for `CODE`/`code_hashes`, a
change commit for `repo_root`). `infra.paths.repo_root` re-exports nothing: E8 callers import
`objdel.core.provenance.repo_root`. Tests: `code_hashes()["files"] >= 100` (be52b1c count re-measured on the E0
commit); editing one `src/objdel` file in a tmp copy changes the hash; AST: 0 `__file__`-based roots under all of
`src/objdel/` (widen U34 and E-G15).

### 5. major | E | The interim contract edit at E3b is a rename that keeps `verify` above `data`; dropping the ignores then fails

Problem. E0 says each move commit "changes its name in the list (... E3b: `reward` -> `verify`; ...)" and E3b drops
the `reward` and `eval.ceiling` ignores (line 9137). A rename gives `["train", "eval", "verify | targets", "data",
"geom"]`. `verify` is then above `data`, so B1's edges `data.* -> verify.reward` (`rl_eligible`) and
`data.quality -> verify.ceiling` point up, and the ignores are gone: `lint-imports` fails on the E3b commit. The same
E0 text says ceiling moves "to `verify`, below `data`", so the plan contradicts itself. E3c, E5 and E2 (which
creates `core`) also need their exact lists.

Fix. Replace the rename sentence with the exact `layers` list per commit (E0 re-derives the edges, but the order is
fixed here):
- E0: `["train", "eval", "reward | targets", "data", "checker | geometry"]`
- E3a: `["train", "eval", "reward | targets", "data", "geom"]`
- E3b: `["train", "eval", "targets", "data", "verify", "geom"]` (ignores left: `eval.vision`, the two `targets`)
- E3c: `["train", "eval", "data", "verify", "geom"]` (ignore left: `eval.vision`)
- E5: `["rl | sft", "policy", "data", "verify", "geom"]` (no ignore)
- E2: add `"core"` last; E8: the target form with `infra` first, `exhaustive = true`.
Add one E0 negative test per list: a planted `objdel.verify -> objdel.data` import fails from E3b on.

### 6. major | E | Package `__init__` re-exports have no move rows; 54 `from objdel.eval import ...` lines depend on them

Problem. `eval/__init__.py` re-exports `gates`, `vision`, `alarms_for`, `hack_audit`, `metrics` and `summarize`
[code: src/objdel/eval/__init__.py:11-13]. `metrics` there is the function, which shadows the submodule name. The
working tree has 54 `from objdel.eval import ...` lines (`vision` 14, `metrics` 13, `gates` 12, `summarize` 8,
`signal` 8, ...) and package-level imports of `objdel.checker` (20), `objdel.data` (13), `objdel.targets` (10)
[measured: grep of src, experiments, tests, lab excluded]. E3a defines `geom/__init__.py` exports. E3b, E3c and E5
define none for `verify`, `data`, `policy`, `rl`, `sft`. After E3b, `objdel.eval` still exists (vision, probes). Its
`__init__` would re-export `metrics` from `verify`, which is a shim (R-2 forbids it), or the 13 callers break. A
mechanical rewrite of `from objdel.eval import metrics` to `from objdel.verify import metrics` imports the module,
not the function. `old_names.json` does not say whether `__init__` re-exports count as names.

Fix. Add to E R-1: "Package `__init__` files define no logic. A move commit rewrites every caller of a re-exported
name to import it from its defining module (`from objdel.verify.metrics import metrics, summarize`) and deletes the
re-export from the old `__init__`. No new `__init__` re-exports, except `geom/__init__.py`'s `polygonize` (the
seam)." `test_move_complete` records definitions only (`FunctionDef`, `ClassDef`, top-level `Assign`), not
`ImportFrom`. Add a test: after E3b, `objdel/eval/__init__.py` holds no `from objdel.verify` import.

### 7. major | E (also DOC) | DOC5's docstring gate makes the E move commits that bring experiment code into `src/` fail lint

Problem. DOC5 adds `D100-D104, D106` to ruff `select` for `src/objdel` and exempts `experiments/**` (line 11262). It
lands before B1's preflight (line 11254), so it is active during E. E2 (`train_sft`, `select_ckpt`, `g_probe`,
`train_rl`), E4 (`sample_vllm`, `hf_sample`, `final_eval`) and E8 (`vast/policy.py`, `vast/bundle.py`,
`preflight.py`, `hf_push.py`, `hf_pull.py`) move code from `experiments/` into `src/objdel/`. Those files have 59
public names with no docstring: 58 D103 + 1 D101; `policy.py` 18, `preflight.py` 10, `train_sft.py` 6,
`select_ckpt.py` 6, `bundle.py` 5, `hf_push.py` 4, `sample_vllm.py` 4 [measured: `ruff check --isolated --select
D100,D101,D102,D103,D104,D106` on those 12 files, ruff 0.17.1, 2026-11-08]. R-1 allows no edit beyond the import
rewrite in a move commit, and R-3 needs `make lint` green. So each of these move commits either fails or breaks R-1.

Fix. Add E rule R-6: "Before a move that brings `experiments/` code into `src/objdel/`, one change commit adds
one-line docstrings (DOC5 style) to every name that `ruff check --isolated --select D100,D101,D102,D103,D104,D106
<files>` lists (59 at be52b1c; re-count on the E0 commit). The move commit then passes lint." List it in E2.5, E4.1
and E8.1. Add to E.8 Mac time: about 1 h.

### 8. major | E | `to_argv` has no flag order, so golden 5's string-for-string test cannot pass; sketch 8 ignores stage order

Problem. Golden 5 is "string for string (list equality)" against today's `chain()` (I13). `to_argv` is "the ONE
definition of each CLI flag" (E6.1), but no rule fixes the order of the flags. A field-order loop emits the mixin
fields first: with `class RlStage(StageBase, TrainerFlags)`, pydantic orders the fields `lr, shared_prefix, grpo,
kind, name, soft, model, ckpt, steps` [measured: pydantic 2.11.7, reviews/r3arch_probe/field_order.py]. Today's
plans put `--shared-prefix` and `--resume` last [code: configs/r5_shared.yaml step 4; rl_shifted.yaml step 2;
rl_full.yaml last step]. The spike (E6.0) tests only that `CliApp` parses the argv; it does not test the order. Also
sketch 8 builds `got = {s.name: to_argv(...)}` and asserts `got == gold`. A dict comparison ignores the stage
order, and the old `steps:` plans have no stage names, so the golden's keys are undefined.

Fix. In E6.1, add `FLAG_ORDER: dict[str, tuple[str, ...]]` beside `to_argv`, one tuple per `kind`, in the order of the
be52b1c plans (for `rl`: `model, ckpt, train, train_views, out, steps, save_steps, resume, shared_prefix`, then
the A-D flags in the order B1.5 writes them). `to_argv` emits fields in that order and skips defaults that today's
argv omits. A test asserts that each tuple covers every field of its stage exactly once. Golden 5 is recorded as an
ordered list of argv lists (stage index = position). Sketch 8 becomes `got = [to_argv(s, "RD", "DATA") for s in
rs.stages]; assert got == gold`. B1.5 writes `b1_rows.yaml` in `FLAG_ORDER`.

### 9. major | E | E8 moves `hf_push`'s body into `objdel/infra/hf.py`; its licence file path then points to nothing

Problem. `APACHE = Path(__file__).with_name("licenses") / "Apache-2.0.txt"` [code: experiments/phase2/hf_push.py:60];
the card upload reads it [code: hf_push.py:75]. The file is `experiments/phase2/licenses/Apache-2.0.txt`
[measured: ls]. E8.1's thin-CLI rule moves the library body to `objdel.infra.hf`. `__file__` is then
`src/objdel/infra/hf.py`, and `src/objdel/infra/licenses/` does not exist. `hf_push` runs on the box after a paid
run [code: policy.py:522-524]. E-G15 and U34 look only for `parents[`, not `with_name` or `parent /`. CLEAN8 step 4
edits this licence text (line 12712) but does not move it.

Fix. In E8.1: move `experiments/phase2/licenses/` to `src/objdel/infra/licenses/` in the same move commit; read it
with `importlib.resources.files("objdel.infra").joinpath("licenses/Apache-2.0.txt").read_bytes()`; add the folder
to the wheel's package data if `wheel.exclude` ever drops infra (ADR 0017 note). Test: the card builder returns
bytes whose SHA-256 starts with `cfc7749b` (the value in the `hf_push.py:60` comment). Widen the E-G15 AST test to
any `__file__`-relative path (`with_name`, `.parent /`, `parents[`) under `src/objdel/`.

### 10. major | E | `EvalStage` leaves out live `final_eval` flags, and E2's "new flags become fields" rule covers only `train_rl` and `g_probe`

Problem. `EvalStage` has `ckpt, view, split, train, out, k` (E6.1). Today's `final_eval` also takes `--model`
(default `qwen3-vl-2b`), `--seed` (default 11) and `--per-call` (default 64) [code:
experiments/phase2/final_eval.py:53-56]. A3 adds `--rescore` and `--hash-view` (lines 1150-1152, 1181). The E4 TOST
pairs the sampler arm with the trainer arm "on the same checkpoint, view, K and seed" (seed 11 [measured:
runs/r5/eval/eval_val/final_eval.json]), but the stage cannot carry a seed. A 4B eval (front section 3) cannot name
its model. After E4, `train` is dead (no trainer). E2.1's rule ("every flag that A-D add to `train_rl` or `g_probe`
becomes a field") leaves out `final_eval`, `select_ckpt` and `train_sft`. With `extra="forbid"`, a ported plan with
one of these flags fails to load.

Fix. `EvalStage`: add `model: str = "qwen3-vl-2b"`, `seed: int = 11`, `per_call: int = Field(64, gt=0)`, and A3's
`rescore: str | None = None`; keep `hash_view` out (a Mac-only mode, not a chain stage); delete `train` in the E4
change commit. Widen E2.1's rule to "every flag of every CLI that a stage maps to (`train_rl`, `g_probe`,
`final_eval`, `select_ckpt`, `train_sft`) on the E0 commit". Widen the E6.0 spike to parse the `eval`, `select`
and `probe` argv of golden 5 as well.

### 11. major | CLEAN (also DOC) | CLEAN12 adds a second doc-path test with another scope; it fails on the plan file on day one

Problem. CLEAN12 step 4 adds `test_doc_paths_exist` over `README.md` and `docs/**/*.md`, exempting only
`docs/changelog.md` and `docs/research/**` (line 12872). DOC4 check 2 tests the same rule with frozen and plan
exemptions and the archive-README exemption. CLEAN12's scope includes `docs/implementation2.0/**`. DOC measured that 285
of the plan's 1,240 backtick paths do not exist yet (line 11140). So CLEAN12's test fails as soon as it lands, and two
tests own one rule with two scopes. CLEAN12 itself says "The number of broken paths today is not measured".

Fix. Delete `test_doc_paths_exist` from CLEAN12 step 4. Write instead: "Doc paths: DOC4 check 2
(`tests/test_docs.py::test_backtick_paths_exist`) is the only test; CLEAN commits keep it green (E R-1 rule, item
3)." Remove it from CLEAN-G9 and the CLEAN.9 rows that list it.

### 12. minor | 00_front, E | `data/views.py` is in the target tree and the layer boxes, but no task creates it; "thin CLIs only" has no plan for the data builders

Problem. Front 2.2, front 2.3 (line 249) and E.3.1 list `views` in `data`. No E3c row, no other task and no move
table fills it (grep: only line 249 names `views.py`). E8.2 says "`experiments/` keeps thin CLIs only" and lists
`build_data`, `build_targets`, `fetch_lidc_pixels`, `prepare_rl`, `build_views`, `build_eval_views` as live CLIs.
They hold 907 lines of logic [measured: wc], and no table moves it. No gate checks thinness.

Fix. Drop `views.py` from front 2.2, 2.3 and E.3.1. In E8.2 write: "Thin CLIs: `train_rl`, `g_probe`,
`final_eval`, `train_sft`, `select_ckpt`, `preflight`, `hf_push`, `hf_pull` (their bodies move in E2, E4, E8). The
data builders in `experiments/phase1/` and `build_views.py` / `build_eval_views.py` stay as they are in this plan; a
later task may move them with a name table."

### 13. minor | E | E4 puts the truncation rule in `policy/sampler.py`, so the RL reward imports the sampler module

Problem. E4.1 moves `rl.reward.truncated` to `policy/sampler.py` "and `rl.reward` imports it". `truncated(ids,
stop_ids)` is a pure rule about token ids. The sampler module holds `VllmSampler`, the child-process entry and the
vLLM environment set-up. The reward then depends on the sampler's interface for one pure function.

Fix. Put `truncated(ids, stop_ids)` in `policy/models.py`, beside `sampling` (the stop ids come from there).
`policy/sampler.py` and `rl/reward.py` both import it from `policy.models`.

### 14. minor | E | `rl/__init__.py` is created "in a separate commit" after the move commit that needs it

Problem. E5 step 1: "In a separate commit, add `src/objdel/rl/__init__.py` with a docstring only". The move commit
that creates `rl/*.py` needs the package file for import-linter and the test layout. DOC5's D104 gate needs a
docstring in every public package `__init__`.

Fix. "The move commit creates `src/objdel/{policy,rl,sft}/__init__.py`, each with a one-line docstring and no
imports (`import objdel.rl.rows` must not import TRL; U40 checks it)."

### 15. minor | E | R-1 names no import-rewrite tool

Problem. R-1 says "`git mv` plus an import rewrite". The rewrite touches 54 `from objdel.eval import` lines and
dozens of package-level imports (item 6), monkeypatch strings and docstrings. The owner rule prefers a library over
own code, and the tool decides how reviewable each move commit is.

Fix. Add to R-1: "The import rewrite is one libcst codemod, `tests/goldens/phase_e/rewrite.py` (input: the move
table as JSON), run in every move commit. The same JSON feeds `old_names.json` and the doc-path rewrite of item 3."
(libcst is a library; `rope`'s move refactoring is the fallback.)

### 16. minor | 00_front, E | V4 has two prices

Problem. Front 6.3 gives V4 "about 39 min = $0.42" and calls E.8's "standalone 15 min" wrong. E.8 gives "26-40 min,
about $0.28-0.43", and E owner decision 7 repeats it. 99_back decision 3b quotes both. One planned rental needs one
number.

Fix. E.8 and E owner decision 7: "about 39 min = $0.42 at $0.641/h [estimate: front 6.3]". Keep the 17-31 min
launch range as the reason. 99_back 3b drops the E.8 range.

### 17. minor | DOC | `FROZEN_GLOBS` freezes folders that hold live code

Problem. `FROZEN_GLOBS` includes `experiments/phase1/**` and `experiments/cuda_checker/**` (line 11135). Frozen
files skip checks 2 and 7. But `experiments/phase1/` holds the live data CLIs that E8.2 keeps (E3c edits
`prepare_rl.py`). `experiments/cuda_checker/run_box.sh` is live: front section 6 changes its toolchain line, and the
chain runs it [code: policy.py:451-453]. A stale path or a wrong version string in those files is not reported.

Fix. Replace the two globs with `experiments/phase1/lab/**`, `experiments/phase1/impl/**`,
`experiments/cuda_checker/results/**` and `experiments/cuda_checker/RESULTS.md`. Run checks 2 and 7 once on the
newly checked files and fix the hits in DOC4's commit.
