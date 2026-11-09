# Final review round 11: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc reviewed: the built
`docs/implementation2.0/Implementation.md` (16,550 lines, 1,689,731 B [measured: `wc -l`, `ls -l`, 2026-11-09]).
Findings are written against its sources in `impl/phases/`. I read `review-log.md` (E, DOC, CLEAN rounds 6-10)
and `final_r10_arch.md` with its verdicts. No item repeats an item that a log marks as handled, and no item
repeats a round-10 minor item. I edited no repo file and ran no git command that changes anything. Probe files are in
`reviews/r11arch_probe/`.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item makes an NVIDIA product claim. No catalog lookup, no install. |
| engineering:code-review | I read the E5 and E3b move tables against the real top-level names of `train/configs.py`, `fast.py`, `serve.py` and `dataset.py` (AST listing, item 1), and the CLEAN8 steps against `grep -rni lfm` of the tree (item 2). |
| engineering:debug | I checked each claim before filing: the `NO_SYNC` and `_bf16` users (`configs.py:140-193`), the `_POOLS` read inside `preload` (`serve.py:53-61`), the docstring `-m experiments.*` strings of live scripts (item 5), and the Mac venv contents (item 6). |
| engineering:system-design | I walked the order front time line step 8 -> DOC13 -> CLEAN waves 2-3 against DOC.4 and DOC13 (item 4), and the E0 coverage record against the E3b archive (item 7). |
| engineering:architecture | Layer check of the E5 table: one shared constant (`NO_SYNC`) placed in `rl` while `sft` reads it breaks the `rl \| sft` independence the plan itself asks for (item 1). |
| engineering:documentation | Ownership of `docs/architecture/layers.md` across E.8, E8 step 5 and DOC9 (item 3); DOC13 position (item 4). I re-ran DOC4 checks 14 and 16 on the built doc: 0 hits each [measured: `reviews/r11arch_probe/chk.py`]. |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |
| system-design-skills:system-design | Same order walk, framed as "which gate breaks first". |
| system-design-skills:service-decomposition | One package, import-linter seams: a modular monolith. Right form; no split proposed. |
| mattpocock-skills:codebase-design | Deletion test on the shared trainer-argument helpers (`NO_SYNC`, `_bf16`): both SFT and RL configs use them, so they belong to the seam below both (`policy`), not inside `rl` (item 1). |
| WebFetch / PyPI JSON | import-linter layers docs: "Exhaustive contracts are only supported for layers that define containers" (the plan sets `containers = ["objdel"]`: OK); `unmatched_ignore_imports_alerting` default is `error` (matches E0). PyPI: import-linter 2.16 (2026-11-06), pydantic-settings 2.16.2, libcst 1.9.0, vulture 2.16: all match the plan. |

## Items

### 1. minor: the E5 move table puts `NO_SYNC` in `rl/config.py`, but `sft_config` reads it; five `configs.py` names, five `fast.py` names and `serve.STORE` have no row.

Sources: E.

- E.md "### E5", move table: "`rl/config.py` | `configs.py`: `rl_lr` (E1), `grpo_config`, `NO_SYNC`" and
  "`sft/config.py` | `configs.py`: `sft_config`". E5 gate: "`lint-imports` shows `rl` and `sft` as independent".
  E.md (after the table): "The table also covers `configs.py`, `fast.py`, `rewards.py`, `callbacks.py`, `sample.py`,
  `serve.py` (E3b)".
- Code: `sft_config` spreads `**NO_SYNC` [code: src/objdel/train/configs.py:172] and calls `_bf16(overrides)` (:165);
  `grpo_config` does the same (:183, :193). So `sft/config.py` must import `NO_SYNC` from `objdel.rl.config`, an
  `sft -> rl` edge that the contract forbids. `_bf16` (:140), `VISION` (:79, read by `lora` and `param_groups`),
  `WARMUP, DECAY` (:80, read by `sft_config` and `train_sft.py:117,288`) and `_RENAME` (:117, read by `as_names`)
  have no row. `fast.py` names with no row: `FINE` (:31), `ROW_KEYS` (:32), `OTHER_VISION` (:33), `PAIR` (:81), `INT`
  (:242). `serve.py` `STORE` (:51) has no row in E3b, and E3b's `preload` reads `_POOLS` (:53, :60-61), a name that
  E3b archives, so a pure move leaves an undefined name.
- Severity: lint-imports (E-G2), `test_move_complete` (E-G16) and ruff F821 all fail on the Mac on the E5/E3b move
  commits, so the engineer finds it at $0. But the table names the wrong home for `NO_SYNC`, and the claim "the table
  also covers" is false.

Fix: in E.md E5 table, move `NO_SYNC` and `_bf16` to the `policy/models.py` row (shared trainer-argument helpers; the
reason: "read by both `sft_config` and `grpo_config`"), add `VISION` and `_RENAME` to the same row, add
`WARMUP, DECAY` to the `sft/config.py` row, and add `FINE, ROW_KEYS, OTHER_VISION, PAIR` to the `rl/forward.py` row
and `INT` to the `policy/caches.py` row (`VisionCache` reads it [code: fast.py:242, 245-268]). In E3b's
`verify/store.py` row add `STORE`, and add one sentence: "a change commit before the E3b move deletes the `_POOLS`
check from `preload` (no pool exists after D7 step 6)".

### 2. minor: CLEAN8 step 2 replaces only `lfm2.5-vl-450m`, but its gate and CT-07 need 0 `lfm` matches; several LFM strings have no step.

Sources: CLEAN.

- CLEAN.md "### CLEAN8" step 2: "Find the files with `grep -rl "lfm2.5-vl-450m" tests`". Step 5: "Gate: `grep -rni
  "lfm\|liquidai" src experiments/phase2/*.py tests` = 0 hits outside frozen records". CT-07: "After Wave 3 (`WAVE
  >= 3`): 0 `lfm` and 0 `liquidai`" over `src/`, `experiments/phase2/*.py`, `tests/**/*.py`.
- Strings that no CLEAN8 step removes [measured: grep, be52b1c]: `src/objdel/geometry/text.py:19` (a measured fact,
  "2.14 characters per token (LFM2.5-VL ...)"), `src/objdel/train/sample.py:22` ("64 / 64 LFM items"),
  `experiments/phase2/g_probe.py:3` (`--model lfm2.5-vl-450m` usage line), `train_sft.py:3` and `preflight.py:3`
  (`r1_lfm.yaml` usage lines), `hf_push.py:102` (`e.g. r1-lfm`), `tests/test_p2_breakit_scoring.py:450-475`
  (`("lfm", ...)` labels), `tests/test_p2_signal.py:18,290`, `tests/test_p2_hf_push.py:66-106` (`r1-lfm`,
  `"lfm1.0"`). The `r1-lfm` HF revision is a kept record (99_back decision 29, CLEAN.Q2).
- CT-07 fails on the Mac at $0, so this is minor. But the engineer must guess whether to delete a measured fact or a
  record name.

Fix: CLEAN8 step 2: "Find the files with `grep -rliI "lfm\|liquidai" src experiments/phase2/*.py tests` (the gate
scope)". Add one rule: a measured fact that names LFM is reworded to name the tokenizer by its HF id in a
`docs/research/notes/` reference, not deleted; test labels `"lfm"` become `"model_a"`; the kept HF revision name is
built at run time (`"r1-" + "lfm"`) like CT-07's needle.

### 3. minor: E and DOC each say the other writes `docs/architecture/layers.md`; DOC13 does not list it.

Sources: E, DOC.

- E.md "## E.8", dependencies: "After E8, DOC writes the layer page and the README layout from E8.5's facts". E8 step
  5: "Hand-off to DOC (DOC writes these pages, E gives the facts): `docs/architecture/layers.md`".
- DOC.md "### DOC9" owner table, row E: "Adds `docs/architecture/layers.md` (links `[tool.importlinter]`)", and
  "Each phase writes its own pages in its own pull request; DOC only owns the table and the checks". DOC13 steps 1-4
  do not name the page.
- E0's contract text points to it from the first E commit: `broken_contract_guidance = "Layer rules:
  docs/architecture/layers.md. ..."` (E0 toml). No DOC4 check reads `pyproject.toml`, so a missing page is not caught.

Fix: one owner. E.md E8 step 5 and E.8: "the E8.4 commit (exhaustive contract) writes `docs/architecture/layers.md`
and adds it to `nav` (DOC4 check 13); DOC13 step 2 only checks it". Keep DOC9's row as is.

### 4. minor: DOC.md puts DOC13 after CLEAN; the front and CLEAN put it before CLEAN waves 2-3.

Sources: DOC, 00_front, CLEAN.

- 00_front.md "## 1", time line step 8: "DOC final sweep, then CLEAN waves 2-3." CLEAN.md "## CLEAN.8": "the DOC final
  sweep (DOC13) runs after E, then CLEAN Waves 2-3."
- DOC.md "## DOC.4": "DOC13 runs last, after phase E and CLEAN have moved code." "### DOC13. Final sweep after phase
  E and CLEAN": "Runs once, after E has moved packages and CLEAN has archived code." DOC.8 lines "DOC13 runs after
  it" (twice) and E05 "DOC13 re-run after phase E and CLEAN".
- Also CLEAN.8 says "Waves 2-3 change no doc text except the archive READMEs and the lines they must fix", while
  DOC9's row CLEAN updates `reference/package.md` and `glossary.md` ("retired" rows): CLEAN8 (Wave 3) retires LFM.
- No result changes (each CLEAN commit keeps docs green), so minor.

Fix: DOC.md: "DOC13 runs after E and CLEAN wave 1, before CLEAN waves 2-3 (front time line step 8). CLEAN waves 2-3
update `reference/package.md` and the glossary "retired" rows in their own commits (DOC9 row CLEAN)." Change the
DOC13 heading to "Final sweep after phase E". CLEAN.8: add "and the `package.md` / glossary rows of DOC9".

### 5. minor: `LIVE_CLIS` misses live scripts that the plan keeps, so E8's string gate fails on their own usage lines.

Sources: E.

- E.md "### E8" step 2: the live CLIs are "build_data, build_targets, fetch_lidc_pixels, prepare_rl, build_views,
  build_eval_views, train_sft, select_ckpt, train_rl, g_probe, final_eval, preflight, hf_push/pull; `LIVE_CLIS`
  (step 1) lists their module names". E8 gate: "every `-m experiments.*` module string in `src/`, `specs/`, `Makefile`
  and live `experiments/` resolves (`find_spec`) and is in `LIVE_CLIS`"; E-G15 the same.
- Scripts the plan keeps live, with a `-m experiments.*` usage line in their docstring [measured: grep]:
  `experiments/phase2/profile.py:4` (E2 step 4 rewrites it as a thin caller), `experiments/phase1/eval_run.py:3`
  (E3b: a torch-free caller of `verify/store.py`; CLEAN6 keeps `hack_audit` for it), `experiments/phase2/baselines.py:3`
  (E2.6: "serves only `judge` and its own CLI"), plus `probes`, `prereg`, `parity`, `phase1/reward_gates`,
  `phase1/bench/xarch_hashes`, `cuda_checker/eval`, A2's `ceiling.py` (CLEAN16 checks it) and D1's
  `perf/bench_checker.py`. CLEAN5 archives only `census_sizerule, eval_demo, sample_eval, smoke_grpo`.
- The gate fails on the Mac at $0. But "live" is not defined, and the engineer must pick between growing
  `LIVE_CLIS` (then `test_cli_smoke`, `test_cli_imports` and `test_bundle_imports_standalone` run them) and archiving.

Fix: E8 step 2: add `profile`, `eval_run`, `baselines`, `probes`, `ceiling` and `reward_gates` to the live list (each
is called by a live CLI, a test or a runbook on the E0 commit; re-check with grep there), and add one rule: "a script
in `experiments/` outside `archive/` and outside the frozen globs is live; E0 lists the others for archive in E8.2".
Scope the string gate to strings outside the script's own module docstring, or require each such script in
`LIVE_CLIS`.

### 6. minor: no step installs the new tools in the Mac venv; Mac `make cov` loses `verify/http.py` at E7.

Sources: E.

- E0 adds `import-linter==2.16` and libcst to the dev extra; E2 adds `pydantic>=2.14,<3`; E6 adds
  `pydantic-settings>=2.16,<3`; E7 adds the `serve` extra. The Mac `.venv` holds none of pydantic, pydantic_settings,
  import_linter, libcst or fastapi [measured: `ls .venv/lib/python*/site-packages`, 2026-11-09; E.md section 2 says
  the same for pydantic]. The Makefile's reinstall rule uses `--no-deps` [code: Makefile:16], so `make test` never
  pulls them. CLEAN12 names its reinstall (`uv pip install -e ".[dev]"`); E names none.
- E7 step 5 adds `serve` to the CI train job only. U55 binds `make cov` "on the E2 move commit and every later E
  commit". On the Mac `tests/verify/test_http.py` skips, so `verify/http.py` is uncovered there, and E2's tests were
  sized to reach the floor, not to keep a margin.
- Loud at $0, so minor.

Fix: E0, E2, E6 and E7 each end with the owner command `uv pip install --python .venv/bin/python -e ".[dev,data,train]"`
(E7: `.[dev,data,train,serve]`). E7 step 5: "the Mac venv installs `serve` too, so `make cov` runs `test_http.py`".

### 7. minor: E0's coverage figure (a) counts the pool code that E3b archives before the E2 move.

Sources: E.

E8 step 1 "Coverage gate": figure (a) = "`objdel` plus the names the E2/E4 tables move", and tests are added "until
(a) reaches that value". E3b archives `_Pool`, `_POOLS`, `CGROUP`, `POOL_CAP`, `quota`, `cpus`, `pool_size`,
`fork_pool` with their tests, before E2. That code is fully covered today (`serve.py` 98 statements, 0 missed
[measured: `reviews/r10arch_cov.txt`]; about 41 statements are pool code [estimate: AST count,
`reviews/r11arch_probe`]). So a figure (a) at exactly 90.0 % is about 89.9 % on the E2 move commit. `make cov` (U55)
catches it on the Mac. Fix: figure (a) also subtracts the names that E1, E3b and E5 delete or archive (E3b pool row,
E5 step 3's `shuffle_sequence_dict` swap), and the test commit targets (a) + 0.3 points.

### 8. minor: Wave 3 CLEAN steps run pre-E commands and names.

Sources: CLEAN.

CLEAN0 step 5: "Re-run steps 2-5 at the start of Wave 2 and Wave 3". Step 2's command imports
`objdel.reward.spec` and `experiments.phase2.vast.policy` and reads `p.ROOT`. CLEAN16 step 3: "Re-run CLEAN0's
command into `runs/clean/after.json`". Wave 3 runs after E, where these are `objdel.verify.reward.spec`,
`objdel.infra.vast.remote.code_files` and `core.provenance.repo_root()` (E8.1 deletes `ROOT`). CT-63 also checks
"`serve` has no `_Pool` ...", but `serve` does not exist after E3b. All fail loudly on the Mac. Fix: CLEAN0 step 2
gives a second command for Wave 3 with the post-E names (as CLEAN0 step 4 already does for preflight); CT-63 checks
`objdel.verify.store` and the absence of the pool names in `src/objdel/**`.

### 9. minor: E still names reward v4 as fixed, but under a raise B3.10 may archive it before E0.

Sources: E, B.

E "When E starts": "Under a raise that funds V2 or V3 before E, E0 also waits for C10 and B3.10". B.md B3.10 pass:
"Move `RewardV4`, `advantages_v4`, their tests and the v4 TRL wiring to `archive/reward_v4/`". E still writes golden 1
as `astuple(REWARD_V4)`, E3b "`advantages_v4` stays here", E7 "`reward="reward_v4"` maps to `REWARD_V4`", and
`grpo_config` reads `REWARD_V4.group`. The golden script fails on the Mac at $0. Fix: one sentence in E0 and E7: "the
live reward on the E0 commit (`REWARD_V4`, or `REWARD_V5` after a B3.10 pass; the fail case archives v5 and the E3b
`advantages_v5` note drops)".

### 10. minor: stale text.

Sources: E, 00_front.

- E.md E.3.1 box and 00_front.md 2.2 box: "callbacks.py (TraceStep: profiler, used by rl and sft)". The E5 table
  also puts `BudgetStop` and `LocGuard` there (final round 9). Fix: "callbacks.py (TraceStep, BudgetStop, LocGuard)".
- E.md E.8 V4 step 3: "on the CUDA build with pool size 1 (as A3)". D7 step 6 removes the pool and `--procs` before
  E0. Fix: drop "with pool size 1 (as A3)".

## Counts

blocking 0, major 0, minor 10 (items 1-10).
