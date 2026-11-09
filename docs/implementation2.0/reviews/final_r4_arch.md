# Final review round 4: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc reviewed:
`docs/implementation2.0/Implementation.md` (assembled 2026-11-08 21:08, 14,787 lines, 1,457,224 bytes). Code checked at
be52b1c (working tree). The review log (`docs/implementation2.0/review-log.md`) was read for E, DOC and CLEAN (rounds
1-3 and the round-4 prep pass). No item below repeats a logged item (grep of the log for `STORE`, `stop_ids`,
`project_r1`, `prereg.select_rb`, `sample.sample`, `import_completions`: 0 hits each; `TraceStep` and
`exhaustive_ignores` hits are about other questions).

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice is in E, DOC or CLEAN. No catalog lookup needed beyond the front matter's; nothing recommended or installed. |
| engineering:code-review | Read every module the E move tables name (`train/{serve,dataset,configs,fast,rewards,callbacks,sample}.py`, `eval/*`, `reward/*`, the experiment scripts that move into `src/`). Found the unlisted `STORE` and the `rl`-only homes of `TraceStep` and `stop_ids`. |
| engineering:debug | Traced each E commit through `make lint` (import-linter, exhaustive contract) and through `make preflight`: found the E8 exhaustive failure (item 1), the E3c circular import (item 6) and the preflight crash after `project_r1.py` is archived (item 7). |
| engineering:system-design | Checked the run path after E end to end (spec -> chain -> CLI -> library). Found that V4 has no spec, so it cannot pass preflight (item 9). |
| engineering:architecture | Checked each layer edge against the target contract `infra > rl | sft > policy > data > verify > geom > core`. Two shared helpers sit in `rl` but `sft` needs them (item 3). |
| engineering:documentation | Checked that front 2.2/2.3, 99_back and CLEAN agree with E and CLEAN.9.1 (items 10, 11). |
| huggingface-skills:huggingface-best | No model choice in these phases. No change. |
| system-design-skills:system-design | Treated the layer plan as a hypothesis and stress-tested each seam with the real import graph (AST scan). |
| system-design-skills:service-decomposition | One process, a modular monolith: the right choice. Checked that no seam needs a shared mutable global across a layer: `STORE` + `_POOLS` split across `verify/store.py` and `verify/pool.py` is such a seam (item 5). |
| mattpocock-skills:codebase-design | Deletion test on `train/sample.py`, `select_ckpt.py`, `train_sft.py`: their names reappear in 3-4 callers, so they need name tables like `vast/policy.py` has (items 2, 4). |
| Web tools (WebSearch, WebFetch) | import-linter layers docs: "If a layer is listed in the contract, the contract will be broken if the layer doesn't exist"; exhaustive contracts fail on unlisted container children unless they are in `exhaustive_ignores` [web: import-linter.readthedocs.io/en/stable/contract_types/layers/]. PyPI JSON, 2026-11-08: import-linter 2.16 (2026-11-06), libcst 1.9.0, pydantic-settings 2.16.2, zensical 0.1.0 (2026-11-05), mkdocs-material 9.8.1 (2026-11-07), mkdocstrings-python 2.0.10 (2026-10-19). All match the doc. |

## Counts

blocking 1, major 8, minor 9.

## Items

### 1. blocking | E | E8 makes the contract exhaustive while `objdel.train` (and `objdel.eval`) still exist; `make lint` fails at E8

Doc (E.md, E8 step 4): "Contract: the `layers` target form of E0 with `exhaustive = true`, `exhaustive_ignores =
["_gpucheck"]`, all layers non-optional, `ignore_imports` empty". U1 (E.9.1) asserts `exhaustive_ignores ==
["_gpucheck"]`. E-G2: "at the end of E8 `exhaustive = true` and 0 `ignore_imports`".

Doc (E.md, E4 step 3): "The V4 follow-up commit (step 1) then deletes `train/sample.py`'s `sample()` ... Until that
commit, `train/sample.py` is the only file left in `train/`." E.8: V4 is "optional, after E; owner approval needed",
and 99_back E-Q7 has "no default". The order is E4 -> E6 -> E8, so at the E8 commit `src/objdel/train/__init__.py`
and `train/sample.py` exist. Also, E3b keeps `objdel/eval/__init__.py` ("the docstring stays"), and no step deletes it
after E5 moves `vision`.

Proof: import-linter: an exhaustive layers contract fails when a child of the container is not a layer and not in
`exhaustive_ignores` [web: import-linter layers docs, "Exhaustive contracts"]. So `lint-imports` reports
`objdel.train` (and `objdel.eval`) at the E8.4 commit. R-3 needs `make lint` green after every commit (owner rule,
BRIEF section 3). The two rules (E8.4 + U1 vs E4.3 + E-G14 "listed exceptions") cannot both hold. Also the TID251
table still says `"src/objdel/train/sample.py" = ["TID251"]  # removed in E4 (file deleted)`, and the E5 table header
says "the package is gone after E5"; both contradict E4.3.

Fix (E.md):
- E8 step 4: `exhaustive_ignores = ["_gpucheck", "train"]` with the comment "`train` leaves with the E4 V4 follow-up
  commit". The V4 follow-up commit deletes `src/objdel/train/` and drops `"train"` from the list.
- U1: `exhaustive_ignores == ["_gpucheck"]` if `src/objdel/train/` is absent, else `["_gpucheck", "train"]`.
- E5 step 1: "The commit that empties a package deletes its `__init__.py`": `eval/__init__.py` goes in the E5 move
  commit (vision was its last module); `checker/__init__.py` and `geometry/__init__.py` go in the E3a move commit.
- E5 table header: "the package is gone after E5, except `train/sample.py` (E4)". TID251 comment: "removed by the
  E4 V4 follow-up commit".

### 2. major | E | `select_ckpt`, `train_sft`, `g_probe`, `train_rl`, `final_eval` move into `src/` with no name table; the D12 logic imports `experiments.*`

Doc (E.md, E2 step 6): "`src/objdel/sft/select.py` takes the D12 selection logic of `select_ckpt.py` (426 lines
[measured: wc]); it reads rows through `policy.rows.prompts` (not rl), so sft never imports rl." E2 step 6 also moves
"the library part of `experiments/phase2/train_sft.py` (440 lines)". U10/E-G16 (`test_move_complete`) cover
`dataset, configs, fast, rewards, callbacks, serve`, the A-D modules and the E8.1 infra sources only, not these
scripts.

Proof: the D12 logic calls experiment modules:
- `select_ckpt.py:397` `from experiments.phase2 import prereg`; `name, why = prereg.select_rb(branches)`;
  `gate = prereg.entry(report, ["lidc"])` (in `main`).
- `select_ckpt.py:271` `from experiments.phase2 import baselines`; `baselines.fit_prior`, `baselines.score_prior`
  (in `judge`).
- `select_ckpt.py:354` `from experiments.phase2.train_sft import vision_check as check`.
- `select_ckpt.py:54` `ROOT = Path(__file__).resolve().parents[2]` (a `__file__` root that E-G15 bans in `src/`).
In `src/objdel/sft/select.py` these imports break CLEAN12's `test_src_never_imports_experiments_or_tests` (E0 relies
on it) and E-G8 ("src imports of `experiments` ... : 0"). The doc does not say whether `prereg` (160 lines) and
`baselines` (213 lines) move into `src` (and where), or whether `main`/`judge` stay in the CLI. That choice changes
the layer graph and which tests must move.

Fix (E.md, E2 step 6): add a name table for `select_ckpt.py`, `train_sft.py`, `g_probe.py`, `train_rl.py`,
`final_eval.py` in the style of the E8.1 `vast/policy.py` table, and add the five files to `old_names.json` (U10).
Suggested homes: `prereg.select_rb`, `prereg.entry`, `G_PROBE`, `g_probe` stay in `experiments/phase2/prereg.py`
(the frozen record), and the CLI `select_ckpt.main` keeps the calls to `prereg` and `baselines`; `sft/select.py`
takes only `score_view`, `sample_view`, `run_sampler`, `run_group`, `complete`, `finished`, `exited`,
`cuda_untouched`, `hf_check`, `prior_ci`, `vision_config`, `TRIES`, `TIMEOUT_S`, `HF_CHECK_*` and `SUBSET_SEED`,
with `judge` taking the prior functions as arguments (dependency injection, as E3c does for `prompts`). `ROOT` is
replaced by `core.provenance.repo_root()`.

### 3. major | E | Two helpers that `sft` needs are placed in `rl`; `rl | sft` independence breaks at E2/E4

Doc (E.md, E5 table): "`rl/callbacks.py` | `callbacks.py`: `TraceStep`"; "`rl/reward.py` | `rewards.py` (TRLReward,
LOGGED, stop_ids, truncated)". E4 step 1 moves only `truncated` to `policy/models.py`. E5 gate: "`lint-imports`
shows `rl` and `sft` as independent."

Proof:
- `experiments/phase2/train_sft.py:390`: `calls = [Rolling(keep_at, stop_at)] + ([callbacks.TraceStep(args.profile,
  WARM)] if args.profile else [])`. E2 step 6 moves this body into `sft/train.py`, so `sft.train -> rl.callbacks`.
- `experiments/phase2/select_ckpt.py:187-201`: `def stop_ids(path, proc)` calls `rewards.stop_ids(...)`;
  `sample_view` (line 216) uses it. Whether it lands in `sft/select.py` or `policy/evaluate.py`, it imports
  `rl.reward` (`sft -> rl` or `policy -> rl`, both forbidden).
- `g_probe.py:93` also uses `TraceStep`, so it is shared by both upper layers.

Fix (E.md, E5 table and E4 step 1): `TraceStep` goes to `policy/marks.py` (the profiling module, beside A6's
NVTX steps) or a new `policy/callbacks.py`. `stop_ids` goes to `policy/models.py` beside `truncated`;
`rl/reward.py` and the sampler import both from `policy.models`. Add U-test: `rl.callbacks` and `rl.reward` are not
imported by any `sft.*` or `policy.*` module (the contract does this, but a named planted-violation case per
helper makes the failure readable).

### 4. major | E | `train/sample.py` has no name table; `g_probe` still calls `sample.sample()` after E4.3 deletes it

Doc (E.md, E5 table): "(E4) `policy/sampler.py` | `sample.py`". E4 step 3: "The V4 follow-up commit (step 1) then
deletes `train/sample.py`'s `sample()` (the only external `_generate` caller [code: sample.py:64])". E4 step 1 moves
only "`build_prompt` ... `sample.export_prompts` logic".

Proof: `g_probe.py:44-47`: `from objdel.train import sample` ... `sample.sample(trainer, rows, 8)` inside `score()`,
the G-probe's K = 8 draw from the trainer in memory. g_probe is a live CLI (`ProbeStage`; the chain's `probe` stage
[code: policy.py:447]; front section 3 step 3 uses an RL probe for the 4B gate). `select_ckpt.py:212-221`
(`sample_view`) uses `sample.export_prompts` and `sample.import_completions`. `sample.py` defines `sample, export,
write, _path, export_prompts, import_completions` [code: sample.py:49-137]; only `export_prompts` gets a home. So
after the follow-up commit, `g_probe` fails at import time, and `import_completions`, `export`, `write` have no
home. `sample()` calls the private `trainer._generate`, so by E-G5 its only legal home is `rl/trl_adapter.py`.

Fix (E.md, E4): add a name table for `train/sample.py`: `sample` -> `rl/trl_adapter.py` (as `generate_k(trainer,
rows, k, items)`, the private-name owner; `g_probe` calls it); `export_prompts` -> `policy/sampler.py`
(`build_prompt` + writer); `import_completions`, `export`, `write`, `_path` -> `policy/sampler.py`. Add `train/sample.py`
to `old_names.json`. The V4 follow-up commit then deletes the empty `train/` package (item 1).

### 5. major | E (also CLEAN) | `serve.STORE` has no home, the name test cannot see it, and `preload` depends on the pool file that CLEAN16 deletes

Doc (E.md, E3b): `verify/store.py` gets "`preload`, `CacheRef`, `VIEW_SEP`, `loader`, `_tagged`"; `verify/pool.py`
gets "`_Pool`, `_POOLS`, `CGROUP`, `POOL_CAP`, `quota`, `cpus`, `pool_size`, `fork_pool`". E5: `test_move_complete`
"records definitions only (`FunctionDef`, `ClassDef`, top-level `Assign`)". E.3.1/E.5: "CLEAN16 step 1 deletes this
file" (`verify/pool.py`); R-2 alternative: "archive these names instead of moving them".

Proof: `src/objdel/train/serve.py:50` `STORE: dict[str, dict[str, dict[str, ExpertCache]]] = {}` is an `AnnAssign`,
the only top-level `AnnAssign` in all moved modules [measured: AST scan of `src/objdel/**`, `experiments/phase2/*.py`,
`experiments/phase2/vast/*.py`, this session]. It is in neither row, and the `Assign`-only test does not record it.
`preload` reads `_POOLS` (`serve.py:59-61`) and writes `STORE` (line 80); `CacheRef` reads `STORE`. So
`verify/store.py` must import `_POOLS` from `verify/pool.py`; when CLEAN16 step 1 (or R-2's archive branch) removes
the pool, `preload` breaks. CLEAN16 step 1 itself lists only `_Pool, fork_pool, pool_size, cpus, quota` and leaves
`_POOLS`, `CGROUP`, `POOL_CAP`, `SCORERS`, which contradicts "deletes this file".

Fix (E.md E3b, E5; CLEAN.md CLEAN16):
- Add `STORE` to the `verify/store.py` row. Make `test_move_complete` record `AnnAssign` targets too (U10 text).
- Put `_POOLS` in `verify/store.py` (it guards the store), and let `verify/pool.py` import it from there; the
  dependency then points from the temporary file to the permanent one.
- CLEAN16 step 1: "delete `verify/pool.py` (all its names) and the `_POOLS` guard in `preload`".

### 6. major | E | E3c: the `prepare(prompts=...)` change must come before the move; the order is not stated, and move-first gives a circular import

Doc (E.md, E3c): "`data/prepare.py` | `train/dataset.py`: `prepare`, `_prepare`, ... Change (separate commit):
`prepare(..., prompts: Mapping[str, str])` takes the prompt table as an argument (dependency injection), because the
manifest stores `PROMPTS` [code: dataset.py:131]." R-1: "A commit is either a move or a change, never both."

Proof: `dataset.py:131` `prepare` reads the module global `PROMPTS` (`"prompts": PROMPTS`). If the move commit comes
first, the codemod must add `from objdel.train.dataset import PROMPTS` to `data/prepare.py`, while `train/dataset.py`
imports `_pixels_sha` from `data.prepare` for `_check` (`dataset.py:172-175`): a circular import at load time. It
also adds `data -> train` with no ignore in the E3c row (`["train", "eval", "data", "verify", "geom"]`, ignores
`eval.vision` only), and `pure-core` (source `objdel.data`, forbidden `datasets`) fails through `train.dataset`'s
top-level `from datasets import Dataset, Image` (`dataset.py:65`). The E3c commit then fails R-3.

Fix (E.md, E3c): "Order: (1) change commit on today's path: `train/dataset.prepare(..., prompts)` with the CLI
`prepare_rl.py` passing `dataset.PROMPTS` (U39 green); (2) move commit." Same rule for any moved name that reads a
global that stays behind: the E0 move-table JSON lists, per moved name, the globals it reads, and the codemod refuses
an upward import.

### 7. major | E, CLEAN | E6.6 archives `project_r1.py`; the preflight `projection` check then crashes and every preflight fails

Doc (CLEAN.md, CLEAN8 step 4): "`project_r1.py`, `r0_parity.py`, `hf_rate.py`, `c2_analyze.py` and their tests move
to `archive/phase2_runs/` ... This move goes in the same commit as E6.6." E6 step 8 rewrites preflight to read a
`RunSpec` but says nothing about its check list. Front 6.4 rule 1: "Preflight fails: no rental."

Proof: `experiments/phase2/preflight.py:250` `CHECKS` holds `("CODE projected GPU hours and dollars", projection)`;
`projection` starts with `from experiments.phase2 import project_r1 as P` (line 217) before any config test, and
`main` turns a crash into `"check crashed: ..."` and `stamp["pass"] = False` (lines 271-276). So after E6.6,
`make preflight` fails on every config, `policy.launch` refuses (`refusal`, policy.py:313-316), and V4 (after E)
cannot launch. `tests/test_p2_hf_eval.py:273-330` (the project_r1 tests) are archived with it, so `make test` does not
flag the dead import.

Fix (E.md E6.8, CLEAN.md CLEAN8 step 4): the E6.6 commit (which archives `project_r1.py`) also deletes `projection`,
`HF_RATE` and `ceiling` from `preflight.py` and its row from `CHECKS`. Add a test: every `CHECKS` function's lazy
imports resolve (`find_spec`) on the post-E6 tree, and a preflight run on the fixture spec has no "check crashed"
problem.

### 8. major | CLEAN (also E) | CLEAN-G13 cannot run in Wave 3: its config is archived by E6.6, and preflight reads only specs after E6.8

Doc (CLEAN.md, CLEAN.6): "CLEAN-G13 | `make preflight CONFIG=experiments/phase2/configs/r1_qwen3vl2b.yaml` ... each
entry of the stamp's `checks` map has the same `pass` value as in `before.json`". CLEAN0 step 5: "Re-run steps 2-5
at the start of Wave 2 and Wave 3 into `runs/clean/before_w2.json` and `before_w3.json`" (step 4 is this preflight).
Wave 3 runs "after phases D and E" (CLEAN16, CLEAN8).

Proof: E6 step 6 archives "the 4 ported originals (their specs replace them)", which include `r1_qwen3vl2b.yaml`.
E6 step 8: "`preflight` takes a `RunSpec` (via `load_spec`)"; the Makefile default becomes `specs/b1_rows.yaml`.
So in Wave 3 the gate's command names a missing file, and a `RunArgs` YAML is no longer a valid input. The stamp's
`checks` map also loses the `projection` entry (item 7), so "same `pass` value per entry" cannot be met even with the
right path. The gate is untestable as written.

Fix (CLEAN.md CLEAN0 step 4, CLEAN-G13): "Wave 1 and Wave 2: `CONFIG=experiments/phase2/configs/r1_qwen3vl2b.yaml`.
Wave 3 (after E6): `CONFIG=specs/r1_qwen3vl2b.yaml`; record `before_w3.json` from that command; compare only the
check names present in both maps, and list removed checks (the E6.6 `projection` removal) by name."

### 9. major | E, 00_front | V4 has no run spec or launch config, so it cannot pass preflight or launch

Doc (00_front, 6.3): "**V4** (optional, after E ...): E-G10 box part (prompt ids and SamplingParams equal; TOST
sampler arm on val), E4 deterministic check and E-G13 ... about 39 min = $0.54". Front 6.4 rule 1: "Preflight fails:
no rental." E6 step 8: "`preflight` takes a `RunSpec` ... `policy.launch` calls `refusal(spec_path)`". E.9:
"All five run at V4 item 'E checks' (E.8) with one file list: `pytest -m gpu tests/verify/test_verifier_gpu.py
tests/infra/test_box_smoke.py tests/policy/test_sampler_box.py`."

Proof: grep of E.md, 00_front.md and 99_back.md for a V4 spec or plan (`specs/v4`, `v4.yaml`, `e_checks`): 0 hits.
V1 has `b1_launch.yaml` + `b1_rows.yaml`, V3 has `a3_test.yaml` etc.; V4 has nothing. The inputs that V4 needs are
not named anywhere as one list: the checkpoint fetch (`fetch_hf.sh <run> rl-full-3ep`, the trainer arm's `ckpt`),
the shipped view (`lidc-val-shifted`; front 6.4 rule 7 moves every other eval view out first), the `EvalStage`
values (K 8, seed 11, `per_call` 64), the `CmdStage` for `env_check.sh` and the pytest file list. The engineer must
guess them, and a wrong K or seed makes the TOST compare two different samplings.

Fix (E.md E.8, front 6.3): add `specs/v4_e_checks.yaml` to E6 step 2 (ported and preflighted with the others):
`compute: {gpu: RTX_5090, cap_hours: 0.75}`; stages: `CmdStage` env_check, `CmdStage` fetch of `rl-full-3ep`,
`EvalStage {ckpt: <fetched>, view: lidc-val-shifted, split: val, k: 8, seed: 11, per_call: 64}`, `CmdStage` pytest
with the E.9 file list. Front 6.3 names it as the V4 launch spec; front 6.4 rule 7's view list for V4 is
`lidc-val-shifted` only.

### 10. minor | 00_front | Section 2.2 and 2.3 disagree with E on `core`, `infra/paths.py`, `verify/pool.py`, `tests/harness/` and the DOC task of `mkdocs.yml`

Doc (00_front 2.2): "`infra   paths, vast (...)`" and "`core    spec (RunSpec, pydantic)`". 2.3: "`core/ spec.py`",
"`infra/ paths.py vast/ ...`", "`mkdocs.yml  docs site config (DOC10)`". E (E8.1): "There is no `infra/paths.py` and
no re-export"; E.3.1: "`core  spec.py (RunSpec, pydantic), provenance.py (repo_root, CODE, code_hashes; stdlib
only)`"; E8.2 creates `tests/harness/`; E.3.1 lists `verify/pool.py` until CLEAN16. DOC3 owns the site
(`mkdocs.yml`), DOC10 owns the ADR. Front order rule 1 makes E's tables the only source of post-E paths, so this is not
a rule conflict, but a reader of section 2 gets a wrong tree.

Fix (00_front): 2.2 `core spec, provenance`; drop `paths` from `infra`; 2.3 `core/ spec.py provenance.py`, `infra/
vast/ bundle.py preflight.py hf.py (+ licenses/) secrets.py`, `tests/ ... tests/harness/ (k3_harness)`,
`mkdocs.yml (DOC3)`.

### 11. minor | 99_back | Decision 32 still states rules for Q-T3 and Q-T4, which CLEAN.9.1 retired

Doc (99_back, decision 32): "Q-T3 exempt `data/`, `build/`, `.venv/`, `.claude/` (CT-15); Q-T4 the live doc-path
check skips with a reason while `docs/` is untracked." CLEAN.9.1: "Q-T3. Retired in revision 5 ... CLEAN has no
doc-path test." "Q-T4. Retired in revision 5 for the same reason." CT-84 still says "(local gate, `docs/` is
untracked, Q-T4)".

Fix: 99_back decision 32: "Q-T1 and Q-T2 have one rule each; Q-T3 and Q-T4 are retired (DOC4 check 2 owns doc
paths)". CLEAN CT-84: drop "`docs/` is untracked, Q-T4" (DOC1 tracks `docs/` before CLEAN).

### 12. minor | E | E.3.2 shows the box running `--spec ... --stage`, but E6.4 says the box chain uses the flag form

Doc (E.3.2): "box: experiments/phase2/train_rl.py --spec ... --stage train". E6 step 4: "The box chain uses the flag
form (`to_argv`), which golden 5 pins."

Fix (E.3.2): "box: `experiments/phase2/train_rl.py <flags from to_argv(stage)>` (golden 5); `--spec/--stage` is the
Mac form".

### 13. minor | E | "core imports nothing that imports torch, even lazily" has no test; U7 allows function-level imports in `core`

Doc (E.3.1): "`core` ... Nothing that imports torch, even lazily." U7: "AST over `core geom verify data`: no
module-level import of torch, ... (function-level and `TYPE_CHECKING` imports are allowed)". E-G4 is a `sys.modules`
test, which cannot see a function-level import.

Fix (U7): for `core` only, ban function-level imports of the ML stack too (AST walk of every node, not only the
module body).

### 14. minor | E | E5 table leaves 13 top-level names without a row (the U10 test catches them, so the engineer guesses at commit time)

Proof [measured: AST of `train/configs.py`, `fast.py`, `rewards.py`]: `configs.py` `VISION`, `WARMUP`, `DECAY`,
`_RENAME`, `_bf16`; `fast.py` `FINE`, `ROW_KEYS`, `OTHER_VISION`, `PAIR`, `INT`; `rewards.py` `GATED`, `AGGREGATOR`,
`_text`. `train_sft.py:117` imports `DECAY, WARMUP` (an `sft` home); `as_names` needs `_RENAME` (policy).

Fix (E5 table): `WARMUP`, `DECAY` -> `sft/config.py`; `VISION`, `_RENAME`, `_bf16` -> `policy/models.py`; `FINE`,
`ROW_KEYS`, `OTHER_VISION`, `PAIR`, `INT` -> with their users (`rl/forward.py` or `policy/caches.py`, by AST);
`GATED`, `AGGREGATOR`, `_text` -> `rl/reward.py`.

### 15. minor | E | The codemod lives in `tests/goldens/phase_e/rewrite.py`

Doc (E R-1): "one libcst codemod, `tests/goldens/phase_e/rewrite.py`". `tests/goldens/` holds recorded data
(`polycore.npz`, `invariance*.json`); a tool there is hard to find and is a module under a data folder.

Fix: `tests/harness/phase_e_rewrite.py` (E8.2 already creates `tests/harness/` for tools that are not tests), run as
`python -I -m tests.harness.phase_e_rewrite <move_table.json>`.

### 16. minor | E | `infra/vast/__main__.py` holds the nine `cmd_*` bodies

Doc (E8.1 table): "`infra/vast/__main__.py` | `cmd_offers`, `cmd_launch`, ..., `FLAG`, `CLI`, `main`". Tests then
import a module named `__main__`, and `python -m objdel.infra.vast` runs the file as `__main__` while a test import
loads it as `objdel.infra.vast.__main__`: two module objects for one file.

Fix: `infra/vast/cli.py` gets the `cmd_*`, `FLAG`, `CLI`, `main`; `__main__.py` is two lines
(`from objdel.infra.vast.cli import main; raise SystemExit(main())`).

### 17. minor | E | `tests/test_secrets.py` reads `secrets.ROOT`, which E8.1 deletes

Proof: `tests/test_secrets.py:8` `ROOT = secrets.ROOT`; line 38 uses it. E8.1: "Replace every `__file__`-based root
in moved files ... with [`repo_root()`]". The test edit is not listed (the codemod rewrites imports, not attribute
reads).

Fix (E8.1 "Repo root, one function"): "`tests/test_secrets.py` uses `core.provenance.repo_root()` in place of
`secrets.ROOT` (same move commit)".

### 18. minor | DOC | DOC1 step 2 quotes a stale size for the plan file

Doc (DOC.md, DOC1 step 2): "That file is one Markdown text file of 1,480,096 bytes [measured: `ls -l`,
2026-11-08]". Measured now: 1,457,224 bytes (`ls -l docs/implementation2.0/Implementation.md`, after the assembler moved
the review logs to `review-log.md`, 250,537 bytes). The 2 MB limit still holds.

Fix: quote the size as "about 1.46 MB (assembled plan, review logs in `review-log.md`)" or drop the number and keep
only the 2 MB limit, which DOC-G13 tests.
