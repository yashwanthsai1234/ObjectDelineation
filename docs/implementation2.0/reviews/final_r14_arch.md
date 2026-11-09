# Final review round 14: architecture (E, code organization, layering, APIs, DOC, CLEAN)

Reviewer role: principal software architect for production ML systems. Doc reviewed: the built
`docs/implementation2.0/Implementation.md` (17,050 lines, 1,752,215 bytes [measured: `wc -l`, `ls -l`, 2026-11-09]).
Findings are written against the sources in `impl/phases/`. I read `review-log.md` (round 13 for D and CLEAN, round 12
for the rest), `final_r13_arch.md` and its verdict, and the round-13 changes as they meet E, DOC and CLEAN: CLEAN5
step 0 / step 9 RIGA gate / CT-89 / CT-73 / CLEAN.Q6, and D2 step 7 / D8 (`tests/variants.py`, `test_hooks()`,
`make test-order`). No item repeats an item that a log marks as handled. The round-13 arch minors that are still open
(r13 items 2-3: "moved by D3" wording in E3a/E.5; `final_eval.py` name table has no home for `view_content_sha256`
and `frozen_entry`) are not filed again. I edited no repo file and ran no git command that changes anything (only
`git ls-files`, `git check-ignore`). Probe files: `reviews/r14arch_probe/` (`build.txt`, `build_clean.txt`).

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No item depends on an NVIDIA product claim (the round-13 changes are CPU-build test plumbing and git hygiene). No catalog lookup, no install. |
| engineering:code-review | Re-checked the round-13 CLEAN5 fix against the real `.gitignore` (lines 37-49) and the ignored files on disk (31 of 31, as the fix says); checked CT-89's name rule against today's committable set (only the 2 allowed tables match). Read E3a's `_gpucheck` rewrite rule against D2's new direct users of the extension (item 2). |
| engineering:debug | Reproduce-then-isolate for item 1: built a copy of today's `docs/` + `src/` with DOC3's exact `mkdocs.yml` and Zensical 0.1.0; the strict build aborts with 3 issues; mapped each line back to its phase source. |
| engineering:system-design | Walked the order D2/D8 (adds `tests/variants.py`, CI step) -> E3a (TID251 ban on `objdel._gpucheck`) -> E8.3 (test files move into layer folders) and asked which planned check sees a broken helper path (items 2-3). |
| engineering:architecture | Seam check of the extension seam: `geom.engine` is "the ONLY importer of `objdel._gpucheck`", but the D2 load guard must look at the extension module itself, not at the facade (item 2). |
| engineering:documentation | DOC3's "No issues found" claim is from 2026-11-08 and predates `review-log.md` and the current assembled file (item 1); stale size, owner-task labels and status text (items 4-6). |
| huggingface-skills:huggingface-best | No model choice in scope. Nothing changed. |
| system-design-skills:system-design | Same order walk as engineering:system-design, framed as "which artifact (helper module, hard-coded test path, docs page) does each later move break, and who guards it". |
| system-design-skills:service-decomposition | One package with an exhaustive import-linter `layers` contract (a modular monolith) is still the right form; no split proposed. |
| mattpocock-skills:codebase-design | Deletion test on `geom.engine` as a facade: the load guard of a test variant is a caller that must see past the facade, so it needs an explicit, named exception (item 2). |
| WebSearch | import-linter docs, layers contract: "exhaustive" fails on a module in a container that is not a listed layer, `exhaustive_ignores` skips named modules, and "exhaustive contracts are only supported for layers that define containers" [web: https://import-linter.readthedocs.io/en/latest/contract_types/layers/]. E0/E8.4's contract (container `objdel`, `exhaustive = true`, `exhaustive_ignores = ["_gpucheck"]`) fits that. No item. |

## Checked and found correct (no item)

- CLEAN5 step 0b block: `.gitignore:38-49` today holds exactly the 12 anchored lines CLEAN5 lists; `git ls-files -o -i
  --exclude-standard experiments/phase0 experiments/phase1/lab experiments/phase2/lab` minus `__pycache__` and `*.log`
  = 31 [measured 2026-11-09]. The step 9 `sed` map matches the move list (`experiments/phase0/results/` ->
  `docs/research/records/phase0/results/`). CT-89's name rule over today's `git ls-files -co --exclude-standard`
  hits only `experiments/phase1/lab/t_data/riga/{inventory,repo_parity}.json`, the two literals CLEAN.Q6 allows.
- Other ignored files outside the CLEAN5 trees are Vast state under `experiments/phase2/vast/` (ledger, offers,
  invoices, known_hosts), which E8.1 keeps in place ("Vast state stays where it is", `.gitignore:25-35` unchanged).
- `python -I -m tests.variants` resolves: `tests/` has no `__init__.py`, but the editable install puts the repo root
  on `sys.path`, so `tests` imports as a namespace package under `-I` [measured: `.venv/bin/python -I -c "import
  tests"` from `/tmp` prints the repo `tests` path].

## Items

### 1. minor: the strict docs build fails on today's assembled files (3 unresolved autorefs), and two `docs/implementation2.0` pages are missing from DOC3's nav.

Sources: DOC (DOC3), A (A.0 Skills used, A review-log row "data 4"), B (B.1).

- DOC.md "### DOC3": "With it excluded, the tree of 2026-11-08 plus `docs/implementation2.0/Implementation.md` builds with
  'No issues found', exit 0, in 3.4 s [measured: docs/implementation2.0/notes/doc/r1build/build.txt, 2026-11-08]."
- Measured 2026-11-09: a copy of today's `docs/` and `src/` with DOC3's `mkdocs.yml` (same plugins, `strict: true`,
  `validation.links.* : warn`, the full `exclude_docs` block) and `zensical==0.1.0`, `mkdocstrings-python==2.0.10`:
  "3 issues found ... RuntimeError: Aborted because --strict flag is set", exit 1
  [measured: `reviews/r14arch_probe/build_clean.txt`]. The three warnings are "unresolved autoref" at:
  - `implementation/Implementation.md:1179` <- A.md line 84 (A.0 Skills used, engineering:debug row): "...CI is
    [0.6252288200704239, 0.6644106494958667] [measured: scratchpad script, 2026-11-08]."
  - `implementation/Implementation.md:3788` <- B.md line 178: "0.740 [measured: diagnosis note] [to be confirmed by A2]".
  - `implementation/review-log.md:414` <- A.md line 2642 (review-log row "data 4"): "CR1 interval [0.5837, 0.6527]
    [measured: 2026-11-08]".
  Two bracket groups separated by one space parse as a reference link; mkdocstrings' autorefs then reports the label.
- `review-log.md` and `Architecture-explained.md` did not exist when DOC3's trial ran. DOC3's nav "Plan" lists only
  `implementation/README.md` and `implementation/Implementation.md`, so DOC4 check 13 reports both pages.
- Caught-before-spend: DOC3 step 2 (`make docs`, expect "No issues found") and step 3 (check 13) fail on the Mac at
  $0 and point at the lines. So minor. But the plan's own measured claim is wrong for the file it will build.

Fix:
1. A.md line 84: write "... CI is 0.6252-0.6644 (bootstrap) [measured: scratchpad script, 2026-11-08]" (or
   "CI = (0.6252288200704239, 0.6644106494958667)"). B.md line 178: "0.740 [measured: diagnosis note; to be confirmed
   by A2]". A.md line 2642: "CR1 interval (0.5837, 0.6527) [measured: 2026-11-08]". Re-assemble.
2. DOC3 "Files", `nav:` Plan: add `implementation/review-log.md` and `implementation/Architecture-explained.md`.
3. DOC3 step 3: after "re-runs the build on the final file", add "and on `review-log.md`. Any `] [` pair outside a
   code span is rewritten as in the fix above; record the clean build in the DOC phase report." Replace the
   2026-11-08 sentence with the result of that re-run.

### 2. minor: E3a's `_gpucheck` ban and rewrite rule do not say how D2's new direct users of the extension are handled.

Sources: E (E0 toml, E3a), D (D2 step 7).

- E.md E3a: "Then add the TID251 ban on `objdel._gpucheck`. In the same commit, rewrite every direct importer to
  `objdel.geom.engine`." The `geom/engine.py` row lists the names `DEFECTS`, `score_polygons`, `score_masks`,
  `check_polygons`, `check_masks`, `build`, `device_count`, `resident`, `add_tables`, `drop_tables`; "This is the only
  importer of the extension." The E0 toml has one non-`src` exception, `experiments/phase2/perf/bench_checker.py`.
- D.md D2 step 7 (round 13): "Load guard: the child asserts `objdel._gpucheck.__file__` is under
  `build/variants/<name>/site` and `set(g.test_hooks()) == <defs>` ... The normal suite asserts `g.test_hooks() == ()`
  (a test in `tests/test_cuda_checker.py`)." So after D, `tests/variants.py` and `tests/test_cuda_checker.py` use the
  new binding `test_hooks()`, and the guard's point is the path of the compiled `.so`.
- A mechanical E3a rewrite gives `engine.__file__`, which is `geom/engine.py` inside the variant `site`: the guard
  still passes, but it no longer proves which `.so` loaded. `test_hooks` is not in the engine name list.
- Caught-before-spend: E3a re-greps importers on the E0 commit, its `__all__` coverage test fails on the missing
  `test_hooks`, and `make lint` (TID251) fails on `tests/variants.py` at $0. The engineer must still choose how the
  guard sees the extension. So minor.

Fix (E.md E3a and the E0 toml): add `test_hooks` to the `geom/engine.py` name list. Add a second non-`src` exception:
`"tests/variants.py" = ["TID251"]  # added in E3a: the D2 load guard checks objdel._gpucheck.__file__ (the .so),
not the facade`. One sentence in E3a: "The D2 load guard keeps `import objdel._gpucheck`; every other test, including
the `test_hooks() == ()` test, goes through `geom.engine`."

### 3. minor: E8.3 moves the test files that `tests/variants.py` names by path, and no E check reads that helper.

Sources: E (R-3, E8 step 3), D (D2 step 7, D8).

- D.md D2 step 7: "`order1` ... and `order2` ...: `tests/test_cuda_checker.py -k masks` exits 0; `poison` ...: the same
  run exits 0; `alias` ...: the same run MUST exit 1"; `mutfg4` runs I10 (in `tests/test_phase_d_core.py`, D.9.10).
  D8: the CI lint + pytest job runs `.venv/bin/python -I -m tests.variants`.
- E.md E8 step 3: "Tests mirror the layers: `tests/{core,geom,verify,data,policy,rl,sft,infra}/`" and
  "`tests/infra/test_test_paths.py`: every string matching `tests/[\w/]+\.py(::\w+)?` in `pyproject.toml`,
  `.github/`, `experiments/` (live), `specs/`, `src/` and `Makefile` names a file ... that exists." `tests/*.py`
  helpers are not in that scope. E.md R-3's per-commit gate is `make test`, `make lint` and the move-invariance test,
  not `make test-order`.
- After the move, pytest gets a missing path and exits 4 (not 0, not 1), so `order1` fails and the CI step goes red.
  That is a $0 check D names (D-G12: "`make test`, `make lint` and `make test-order` pass"), so minor. But the move
  commit itself passes R-3, and the red shows one push later, away from the commit that caused it.

Fix (E.md E8 step 3 and R-3): widen `test_test_paths.py` to every non-`test_*.py` module under `tests/` (today
`tests/variants.py`, `tests/docs_check.py`, `tests/oracles/*`, `tests/harness/*`). Add `make test-order` to R-3 for
each E8.3 move commit. List `tests/variants.py` beside `preflight.py:29-34` in E8 step 3's "Hard-coded test paths"
list.

### 4. minor: DOC1 step 2 gives a stale size for `Implementation.md`, and the 2 MB limit has 248 KB of headroom left.

Sources: DOC (DOC1 step 2, DOC-G13).

- DOC.md DOC1 step 2: "That file is one Markdown text file of 1,480,096 bytes [measured: `ls -l`, 2026-11-08]; its own
  limit is 2 MB (DOC-G13)." Today it is 1,752,215 bytes [measured: `ls -l`, 2026-11-09 03:07], 272 KB more in about
  one day of review rounds. DOC-G13 fails at $0 if it passes 2 MB, so minor.

Fix: DOC1 step 2: "its size at the DOC1 commit (1,752,215 bytes on 2026-11-09; re-measure)". Add to the assembler
(or DOC4) a warning at 1.9 MB that names the largest phase source, so the owner can split before DOC-G13 fails.

### 5. minor: front 2.3 says DOC10 owns `mkdocs.yml`; DOC3 creates it (DOC10 is the ADR process).

Sources: 00_front (2.3).

- 00_front.md 2.3 tree: "`mkdocs.yml                 docs site config (DOC10)`". DOC.md "### DOC3 ... Files to create
  or change: `mkdocs.yml` (repo root)"; "### DOC10. ADR process".

Fix: "(DOC3)".

### 6. minor: `docs/implementation2.0/README.md` disagrees with DOC3 and with the review state.

Sources: DOC (DOC1 step 0, DOC3); the file itself is not a phase source, it is written at hand-off.

- DOC.md DOC3: "`implementation/README.md` (new, DOC1 step 0) names the copied folders as backtick paths, not as links,
  because those folders are not site pages." Today's file links them: "`[kernels/](kernels/)`, `[specs/](specs/)`,
  `[tests/](tests/)`, `[reviews/](reviews/)`, `[notes/](notes/)`, `[seeds/](seeds/)`" and says "phase reviews, then
  final rounds 1-10", "Ten final review rounds ran" [measured: `docs/implementation2.0/README.md`]. Zensical 0.1.0 did not
  warn on the folder links in item 1's build, so only the text is wrong.

Fix (DOC1 step 0): rewrite the six folder links as backtick paths (the DOC3 rule), and set the status lines from the
review log at hand-off ("Fourteen final review rounds ran; confirmed blocking plus major per round: ...").

## Counts

blocking 0, major 0, minor 6 (items 1-6).
