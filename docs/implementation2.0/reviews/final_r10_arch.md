# Final review round 10: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc reviewed: the built
`docs/implementation2.0/Implementation.md` (16,352 lines, 1,665,519 B). Findings are written against its sources in
`impl/phases/`. I read `review-log.md` round 9 (E, DOC, CLEAN, front) and `final_r9_arch_verdict.md` and
`final_r8_arch_verdict.md`. No item below repeats an item that was fixed. I edited no repo file and ran no git
command that changes anything.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes an NVIDIA product claim. No catalog lookup and no install. |
| engineering:code-review | I opened `Makefile:20-27`, `.github/workflows/tests.yml:30-44`, `pyproject.toml:40-80`, `src/objdel/train/configs.py:70-75`, the imports of `preflight.py`, `vast/policy.py`, `bundle.py`, `hf_push.py`, `hf_pull.py`, and the top-level names of `select_ckpt.py`, `g_probe.py:30-36`, `final_eval.py:60-66`. |
| engineering:debug | Two probes. (1) vulture 2.16 on a pydantic stage model with a `model_validator` (item 1). (2) The full fast suite with `--cov` on the scripts that E2 and E4 move into `src/objdel` (item 2), then per-function coverage from `coverage json`. |
| engineering:system-design | I walked the E commit order (E0 -> E1 -> E3 -> E5 -> E2 -> E4 -> E6 -> E8) against each repo-wide gate that runs on every commit: `make test` (with CLEAN12's vulture test inside), `make cov` (E-G1), `make lint`. Items 1 and 2 came from this walk. |
| engineering:architecture | I checked the ownership of the per-model table (`MODELS`, `MAX_COMPLETION`) across C7, front section 3 (4B step 1) and CLEAN8 (item 3). |
| engineering:documentation | I re-ran DOC4 checks 15 and 16 on the built doc: both pass. I checked DOC1 step 0's copy table against the evidence files that E cites (item 7). |
| huggingface-skills:huggingface-best | No model choice is in my scope. Nothing changed. |
| system-design-skills:system-design | Same commit-order walk as engineering:system-design, framed as "which gate breaks first". |
| system-design-skills:service-decomposition | One process, one package: the plan is a modular monolith with import-linter seams. That is the right form. No service split is proposed, and I propose none. |
| mattpocock-skills:codebase-design | Deletion test on `policy/callbacks.py`: it now holds `TraceStep`, `BudgetStop`, `LocGuard`. `LocGuard` is RL-only (it reads `kind/*/loc`), so a seam in `policy` is shallow. The outcome does not change, so I did not file it. |
| WebFetch | vulture README: unused attribute, class, function, method, property and variable all get 60 % confidence; `--ignore-decorators` is the documented escape. |

## Items

### 1. major: CLEAN12's dead-code test turns `make test` red on the first E commit that adds a pydantic validator. No phase adds the ignore.

Sources: CLEAN, E.

- CLEAN.md "### CLEAN12" step 2: `min_confidence = 60`, `ignore_decorators = ["@pytest.fixture"]`, and "Re-run ...
  and trim any name that no longer matches (keep the list minimal)". Step 4: `test_no_dead_code` asserts "exit code 0
  and empty stdout", and it "FAILS" if vulture is missing. It sits in `tests/test_repo_hygiene.py`, which `make test`
  runs on every commit from Wave 1 on.
- E.md "### E2" step 1 adds `@model_validator(mode="after") def _save(self): ...` to `RlStage`. E6 step 1 adds
  `RunSpec._rules`, `EvalStage._groups` and before-validators. E R-3 requires `make test` with 0 failed after each
  commit. E.md has 0 hits for `vulture`, `make dead` or `ignore_decorators` (grep).
- Proof [measured: `reviews/r10arch_probe/` (`src/spec.py`, `pyproject.toml`, `out.txt`), `uvx vulture==2.16`,
  2026-11-09]. A pydantic stage with a `model_validator` gives exit 3 and
  `src/spec.py:11: unused method '_save' (60% confidence)`. It also gives `unused variable` for fields that are
  read only through `model_dump` or a generic loop. The vulture README gives 60 % for an unused method, function and
  variable [web: github.com/jendrikseipp/vulture]. The repo has 0 pydantic validators today (grep of `src`,
  `experiments/phase2/*.py`, `tests/*.py`), so CLEAN12's measured "4 findings" does not cover them.
- Effect: the E2.1 change commit fails R-3. CLEAN12 says to keep the list minimal and has no rule for new library
  idioms. So the engineer must guess the fix: an ignore entry, a whitelist file (CLEAN12 rejects a whitelist file:
  ruff F821/B018), or a rename. E7's FastAPI route functions (`verify/http.py`) hit the same rule.

Fix: in E.md E2 step 1, add one bullet. The E2.1 change commit adds `"@model_validator"`, `"@field_validator"` and
`"@field_serializer"` to `[tool.vulture] ignore_decorators`, with a `#` comment line ("pydantic calls these by
decorator"), as CLEAN12 step 2's comment rule asks. It also adds any field that only `to_argv`/`model_dump` reads
(for example `schema_version`) to `ignore_names` in the same commit. E7 step 5 adds `"@app.*"` (FastAPI routes) in its
commit. In CLEAN.md CLEAN12 step 2, add one line: "Phase E adds the pydantic and FastAPI decorator ignores (E2.1, E7);
they are not stale names". CLEAN-G9 stays "exit 0, 0 lines".

### 2. major: the coverage gate covers only the E8.1 move. The E2 and E4 moves already push `objdel` below 90 %, and E-G1 binds on every E commit.

Sources: E.

- E.md "## E.6" E-G1: "Line coverage of `objdel` >= the `--cov-fail-under` value of the E0 commit (90 at be52b1c
  ...) on every E commit". The code gate is `--cov=objdel ... --cov-fail-under=90` [code: Makefile:26;
  .github/workflows/tests.yml:44].
- E.md "### E8" step 1, "Coverage gate": "Today the files this step moves are outside `objdel`, so they are not
  measured. After the move they are." Its E0 record measures only `--cov=experiments.phase2.vast
  --cov=experiments.phase2.preflight --cov=experiments.phase2.hf_push --cov=experiments.phase2.hf_pull`. The "add
  tests first" rule runs only "before the E8.1 move". U55 is scoped "On the E8.1 move commit and every later E
  commit".
- But E2 (step 6 name tables) and E4 (step 2) move the library bodies of `train_sft.py`, `select_ckpt.py`,
  `g_probe.make_trainer`, `train_rl.main` (as `run_rl`) and `final_eval.per_item/summarize` into `src/objdel`. That
  happens four tasks earlier. The same sentence applies: today these files are outside `objdel`.
- Proof [measured: `reviews/r10arch_cov.txt`, full fast suite on the Mac at be52b1c, 1447 passed, 2026-11-09]:
  `train_sft.py` 293 statements / 166 missed (43 %; `train` 97/74, its Trainer subclass and `Rolling` callback 0 %),
  `train_rl.py` 45/24, `g_probe.py` 121/40 (`make_trainer` 20/20), `select_ckpt.py` 269/35, `final_eval.py` 57/25.
  `objdel` itself is 3,866/220. Only the names in the E2/E4 tables move. Then `objdel` after E2+E4 is about
  4,289-4,349 statements with 452 missed, which is **89.5-89.6 %** [computed: `reviews/r10arch_probe/e2_move_estimate.txt`,
  per-function data from `coverage json`]. So `make cov` fails on the E2 move commit, and E-G1 fails there.
- Second effect: the E8 rule's E0 figure leaves out these E2/E4 files. So "tests until it reaches that value" at E0
  does not make `make cov` pass at E8.1, where `objdel` already holds the E2 code.
- A-D change both counts, as the E8 bullet says. But the plan gives the engineer no rule for the E2 commit, and
  under the E8 rule's own reasons (no `omit`, no lower floor) the only legal fix is new tests. The engineer must find
  this out at E2.

Fix: in E.md E8 step 1 "Coverage gate", widen the E0 record. It covers every file that E2, E4 and E8 move into
`src/objdel`: add `--cov=experiments.phase2.train_sft --cov=experiments.phase2.select_ckpt
--cov=experiments.phase2.g_probe --cov=experiments.phase2.train_rl --cov=experiments.phase2.final_eval
--cov=experiments.phase2.sample_vllm --cov=experiments.phase2.hf_sample`. Then record two figures in
`runs/E/baseline.txt`: (a) `objdel` plus the E2/E4 files, and (b) (a) plus the E8.1 files. If (a) is below the
floor, one change commit before the E2 move (after its R-6 docstring commit) adds fake-only tests in `tests/sft/`
and `tests/rl/`. First targets: `train_sft.train` with a stub `Trainer` and a 2-row fixture view, `Rolling.on_save`
and `on_step_end`, `vision_reference`, `replace_dir`, `latest`, `write_json`, `preload`, `make_trainer` on the tiny
offline model, and `select_ckpt.score_view`/`prior_ci`/`score_items`. The same rules apply (no `omit`, no
`pragma: no cover`, no lower floor). Move U55's scope to "the E2 move commit and every later E commit". Add E2.6 to
U55's row in the E.9.8 coverage map.

### 3. major: CLEAN8 step 3 writes literal per-model values, which overwrite values that C7 and the 4B step may have set by owner decision.

Sources: CLEAN, front (section 3), C.

- CLEAN.md "### CLEAN8" step 3: "`MODELS = {"qwen3-vl-2b": ...}`; `MAX_COMPLETION = {"qwen3-vl-2b": 1536}`".
  CLEAN8 is Wave 3, after E.
- C.md "### C7. Outline limit (owner decision, from the first run after Gate B3)", Files: "change
  `src/objdel/train/configs.py` `MAX_COMPLETION` only after the owner decides". C7 step 1: the owner picks L
  ("Direction recommends about 1,200").
- 00_front.md "## 3. Model decision", 4B step 1: "Add `"qwen3-vl-4b"` ... to the model map and to `MAX_COMPLETION`,
  `_LORA`, `_GRAD_CKPT`: in `policy/models.py` after E ... CLEAN8 keeps `_LORA` and `_GRAD_CKPT` (CLEAN must
  follow)." The front makes CLEAN keep two names. It says nothing about the `MODELS`/`MAX_COMPLETION` rows.
- Code today: `MAX_COMPLETION = {"smolvlm-256m": 1536, "lfm2.5-vl-450m": 1024, "qwen3-vl-2b": 1536}`
  [code: src/objdel/train/configs.py:75]. CLEAN8's job is to delete the LFM row. Its literal also resets the
  `qwen3-vl-2b` value and drops any row that is not `qwen3-vl-2b`.
- Effect: say the owner set L (C7) or added the 4B row before Wave 3. An engineer who types CLEAN8 step 3 exactly
  then reverts the outline limit to 1,536, or deletes the 4B model. The silent change of `max_completion_length`
  changes the Dr. GRPO normaliser of every later run (C.md section 2: "the outline limit rescales the loss (C7)").
  CLEAN-G10's literal subset test (CT-06) pins the pre-CLEAN values, so it would fail. CT-06 then contradicts step 3.

Fix: in CLEAN.md CLEAN8 step 3, replace the two literals with: "delete the `lfm2.5-vl-450m` keys from `MODELS`,
`MAX_COMPLETION` and every other per-model table. Keep every other key and value as it is on the tree (C7's L and
front section 3's 4B row included)." Add a test: after CLEAN8, `MAX_COMPLETION` minus the LFM key equals its
pre-CLEAN8 value (the same literal-subset pattern as CT-06).

### 4. minor: B's I42 still pins "flag defaults from the environment". The E6.4 commit re-points U45 only.

Sources: E, B.

E.md E2 step 1 ("Session-spend flags after E6.4") says: "The E6.4 commit re-points U45's 'all given build
`BudgetStop`' case to `session_budget`". B.md I42 asserts "`train_rl` takes them as the flag defaults". After E6.4
that statement is false by design: `cli_parse` reads no environment. I42 is not named, so the E6.4 commit goes red
on it. The engineer would re-point it the obvious way. Fix: add "and B's I42 (env -> `session_budget`, not
`cli_parse`)" to the same E2.1 bullet and to E6.4 "Env defaults that leave the parser". Add I42 to U54's row in the
E.9.8 coverage map.

### 5. minor: E2.3 and E6.4 disagree on which commit creates `rl.run.session_budget`.

Sources: E.

E2 step 3 says `run_rl`'s "first call is `session_budget(stage)`". E2 runs before E6. E6.4 says "in the E6.4 commit
these three reads move to `rl.run.session_budget`". Both orders work, because between E2 and E6.4 the argparse env
defaults fill the fields first. But the text should name one commit. Fix: E2 step 3 creates `session_budget`
(the env fallback is then a no-op, because the argparse defaults already filled the fields). E6.4 only deletes the
argparse env defaults. Write that sentence in both places.

### 6. minor: R-4 calls a merge "mechanical" before E6.4. But from E2 on, a new `train_rl` flag also needs an `RlStage` field.

Sources: E.

E.md R-4: "Before E6.4 this is a mechanical change." E2 step 4 makes `train_rl` "argparse -> `RlStage` -> `run_rl`".
`RlStage` is `extra="forbid"`, and E2 step 1 says each new flag "becomes a `TrainerFlags` or `RlStage` field in the
same commit". So an `after-raise` merge between E2 and E6.4 that brings `--numerics` or `--reward v5` must also add
the field, its `FLAG_ORDER` entry and a golden-5 row (R-4 (iii)). Else `RlStage(...)` raises at `train_rl` start.
U52 lands only in E6, so no Mac test catches it before E6. Fix: R-4 "Before E2 this is a mechanical change; from E2
on, step (iii) applies too".

### 7. minor: DOC1 step 0 does not copy the probe and measurement files that E cites as proof.

Sources: DOC, E.

DOC.md "### DOC1" step 0 copies `impl/reviews/*.md` only. E cites `reviews/r9arch_probe/env_probe.py` (E2.1),
`reviews/r9arch_cov.txt` (E8 coverage gate) and `docs/implementation2.0/reviews/r3arch_probe/field_order.py` (E6 flag
order). This round adds `reviews/r10arch_cov.txt` and `reviews/r10arch_probe/`. These are `.py`/`.txt` files and
subfolders, so they are lost when the scratchpad is cleared. DOC4 check 2 skips plan files, so no test sees the dead
citation. Fix: add one row to the step 0 table: `impl/reviews/*.{py,txt,yaml}` and `impl/reviews/r*_probe/` ->
`docs/implementation2.0/reviews/` (skip `docrev/`, `.venv*`, `.coverage*`).

### 8. minor: DOC1 step 2 and DOC-G13 quote a stale size.

Sources: DOC.

DOC1 step 2 says "1,480,096 bytes". The built file is 1,665,519 B today [measured: `ls -l`, 2026-11-09]. It is still
under the 2 MB limit. Fix: "about 1.7 MB at round 10; DOC1 re-measures with `ls -l` before the commit".

### 9. minor: CLEAN8 timing words disagree with E.

Sources: CLEAN, E.

The CLEAN8 heading says "(Wave 3, after phase E2)", but the body says "CLEAN8 runs after E5, E2 and E6". Step 4
says "`HF_RATE` ... in `objdel.infra.preflight` if E8 has landed" and "This move goes in the same commit as E6.6".
E6.6 always comes before E8, and E6 step 8 deletes `HF_RATE` in that same E6.6 commit. Fix: heading "(Wave 3; step 4
inside the E6.6 commit, steps 1-3 after E6)". In step 4, drop the `HF_RATE` clause ("E6.6 deletes `HF_RATE` with
`projection`").

## Counts

blocking 0, major 3 (items 1-3), minor 6 (items 4-9).
