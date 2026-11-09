# DOC1: audit of REPO/docs/ and README.md, and the final docs and repo tree

## Skills used
- nvidia-skills:nvidia-skill-finder: loaded. Catalog check not run (no NVIDIA product question in a docs audit; no install). No change to output.
- engineering:code-review: used the Correctness and Maintainability lenses on docs-vs-code (stale constants, dead links, unchecked claims).
- engineering:debug: reproduce/isolate method for each stale fact (grep the constant, find the source of truth, name the owner file).
- engineering:system-design: shaped the "docs tree by reader and by lifecycle" proposal and the per-phase doc impact table.
- engineering:architecture: ADR template check; proposed ADR-0012..0017 as the new decisions.
- engineering:documentation: README / runbook / architecture / onboarding doc types, "link, don't duplicate".
- huggingface-skills:huggingface-best: loaded; no model comparison needed for DOC1 (model choice belongs to other tasks).
- mattpocock-skills:domain-modeling: glossary rules (no implementation detail in the glossary; ADRs only if hard to reverse + surprising + real trade-off). Used to judge which new items get a glossary row versus an ADR.
- Web tools loaded (WebSearch, WebFetch). Not needed: DOC1 is an internal consistency audit; no external fact is asserted.

## Method
Read docs/README.md, changelog, package.md, ADR index, runbooks; ran a link/path checker over docs/**/*.md and README.md (my script, scratchpad/linkcheck.py); grepped constants against code; compared configs listed in docs with experiments/phase2/configs/. No REPO file edited, no git write.

## Findings

### F1. Docs are largely accurate (what is good)
- Every `src|experiments|tests|cuda|docs|archive/...` path quoted in backticks in docs/**/*.md and README.md exists [measured: linkcheck.py, BADPATHS = []].
- Every relative markdown link in docs/ resolves. The 5 README report links flagged are only percent-encoded spaces (`%20`); the files exist [measured: ls docs/research/reports].
- All 16 YAMLs in experiments/phase2/configs/ appear in docs/reference/configs.md [measured: loop over configs].
- Hard stop in docs (17.5) matches code [code: experiments/phase2/vast/policy.py:58 `ALARM, HARD_STOP = 10.0, 17.5`]; docs/runbooks/gpu-run.md:116,121 and ADR-0009:74,81 already say the docstring is stale.

### F2. Stale facts in code and comments that the docs point at (docs are right, code is wrong)
1. policy.py module docstring says "$10 alarm, $15 stop" [code: policy.py:16] vs constant 17.5 [code: policy.py:58]. Open action item ADR-0009 #3.
2. Dockerfile.phase2 header says "TRL 1.16.0" [code: experiments/phase2/docker/Dockerfile.phase2:6]; pins are trl==1.17.1 [code: pyproject.toml train extra]. Need to check the actual `pip` line and locks in the same file (train.txt is rewritten by build_images.sh --lock per perf/implementation.md:91).
3. train_sft.py docstring cites "TRL 1.16 SFTTrainer" [code: experiments/phase2/train_sft.py:49]. `profile.py:7` already says 1.17.
4. configs/r2_qwen3vl2b.yaml:15 says "traced against TRL 1.17.1 but first run on a GPU here" (history; r2 has since run: PASS per brief). Stale comment.
5. PLAN.md (D3, D15) and JUDGE_R1/R2 say TRL 1.16.0 and vLLM 0.32.0 pins. These are frozen research records. Do not edit; add a one-line banner "superseded by TRL 1.17.1 (perf/implementation.md B1)" at the top of PLAN.md only.
6. Hard-stop wording in frozen records ($15): PLAN.md:286, JUDGE_R2/R3/R4, V0_V1.md:109 (`ALARM, HARD_STOP = 10.0, 15.0`). Same treatment: banner, not edit.
7. tests/test_vast.py has hard-stop tests named `..._at_15`, `_at_16_60`, `_at_17`, `_at_17_50` [measured: .pytest_cache nodeids]. Names are history of the raises; after the $3 budget change they must be rewritten around a budget argument, not another constant bump.

### F3. Gaps in docs/ (things a new engineer cannot find)
1. No doc for the NEW target system. docs/ describes what exists; Architecture.md in docs/research/research2.0 is "not yet built". Phases A-E need a living "current vs planned" page (see R5).
2. No page on the `--train-views` data pipeline end to end (build_views -> train_rl rows 24,624 = 8,208 centred + shifted + rescaled). Only a flag row in guides/training.md:106 and the offset sentence at :151. Data quality gates (coverage audit, leakage checks) have no page. B1 is the highest owner priority and has no home. [code: src/objdel/eval/coverage.py, audit.py; experiments/phase2/build_views.py]
3. No profiling/NVTX doc (A, C): docs/testing.md has sanitizers (:112-113) but no performance-measurement guide, no gate numbers (s/step budget, GPU idle share 40 %).
4. No doc of the TRL private-name start-up check (A) although ADR-0005:59 admits coupling to TRL 1.17.1 internals. Needs a "TRL coupling table": each private name FastGRPO touches, the test that fails if it moves.
5. No runbook for the HF repo (svenk029/objdel-sft, token, push/pull): hf_push.py, hf_pull.py, hf_rate.py exist; only mentioned in passing. Needs docs/runbooks/hf-checkpoints.md (no token values, only names).
6. No doc for CI: .github/workflows/{tests.yml,image.yml}; docs/testing.md mentions CI/gates job but not image.yml, no "GPU CI job" (D) section.
7. No doc for the 3D voxel direction beyond figures/voxels.svg and checker/volume.py (package.md:`checker.volume`). Fine until after E; reserve a stub in the tree, do not write CAD.
8. No ARCHIVE policy doc: archive/cpu_checker/README.md is good but contributing.md:61 is one line. Needs the rules (R7).
9. No per-run record index. 44 GB of runs/ (git-ignored) hold r1..r5, rl, v1-v3, k1, p1, preflight [measured: du -sh runs]; results live in docs/research/notes/*.md with ad hoc names (explanations_2026-11-07.md, issues_solved_2026-11-07.md ...). No index of "run id -> config -> commit -> note -> HF tag".
10. docs/research/ is a dumping area: notes/ mixes 6 loose .md-less folders ("GPU delineation RL mix and match", "NVIDIA CUDA performance for VLM training", "RL ready implementation", "RL speed stack on RTX 5090", "Phase 1 stack and algorithms", "Phase 2 SFT", "GPU delineation use case map") that duplicate names in reports/, plus 2212.07020v1.pdf (untracked) in docs/research/ root [measured: ls docs/research/notes, git status].
11. docs/figures holds 8 figures; no list of which doc uses which, nor of how to regenerate (hero/pipeline/gpu_timeline SVGs hand-made?). Open question.
12. docs are not tested: docs/README.md:60 "These documents are not checked by a test." Only README.md has tests/test_readme.py (quick-start exec + make targets).

### F4. Facts to re-verify when numbers change (rot risk)
- Hardcoded run numbers duplicated across README.md, docs/changelog.md, guides, ADRs: test IoU 0.645/F 0.816/P 0.695/R 0.748 [measured: runs/r5/eval/eval_test/final_eval.json], 23.6 vs 31.0 s/step, $4.56, 770 steps, "1,484 tests / 1,453 default / 31 slow" (docs/testing.md:5). Each phase will change several. docs/README.md says "a number that changes lives in a results file" but quotes it anyway.
- docs/testing.md:5 suite size and package.md "Package version 0.1.0" will go stale at the first test or version change.
- Offset medians: docs quote 0.018 (train rows) and 0.259 (shifted VAL view) [docs/guides/training.md:151, changelog]; BRIEF quotes 0.265 for TEST. Both right (different split); docs must always say which split.
- reference/reward.md describes REWARD_V4 [package.md:236 fingerprint]; reward v5 (B3) changes the term table, the fingerprint, advantages (RLOO, no std division), and the "rewards hash" claims in pyproject comments ("topology results are hashed").
- docs/runbooks/gpu-run.md:115-118 and ADR-0009:40,67 describe ALARM 10 / HARD_STOP 17.5 and "$15.59 of $17.50 committed after R4". The owner's new budget is $3 for the new GPU scripts. That is a different accounting window: it needs a new constant and ledger rule, not a fifth bump (see R9).
- README.md layout block says `docs/  research reports and notes, figures` (README.md ~line 177); docs/ now has 11 engineering areas. Stale.
- README.md "TRL 1.17" badge and speed table are right today [code: pyproject.toml].
- ROADMAP.md was "Written 2026-10-18 ... Revised the same day" and still shows phase table rows with "—" status [code: ROADMAP.md:1-30]. It predates Phases 1-2 results, the CUDA checker, and A-E. docs/README.md:51 calls ROADMAP.md "the plan", but Direction.md is the live plan.

### F5. Naming and structure issues
- Two vocabularies for phases: "Phase 0/1/2" (experiments/ dirs, frozen) vs A-E/DOC/CLEAN (Direction.md). `experiments/phase2/` holds ~35 live scripts including everything RL-related; naming it "phase2" will be wrong the moment A-E starts.
- `experiments/phase1/{lab,impl}`, `phase2/{lab,impl}` (6.6 MB + 1.5 MB + 0.5 MB), `phase2/JUDGE_R*.md`, `V0_V1.md`, `RESEARCH.md`, `EXECUTION.md`: research record mixed beside runnable scripts.
- `experiments/cuda_checker/` is a phase-less experiment (eval, compare, run_box.sh).
- Two package-style trees for tests: 71 files; 18 named test_p2_* and 11 test_breakit_*/test_p2_breakit_* mixed with unit tests [measured: ls tests | grep -c test_p2 = 18]. File names carry history (p2), not subject.
- docs/research/research2.0/ has a space file name "previous vs new.md": links need `<...>` (docs/README.md already does this). Rename to previous-vs-new.md in the move.
- docs/changelog.md duplicates git log by design; fine, but it needs an owner rule (append on each merged phase).
- Glossary: one word per thing is respected for image/outline (grep for scan/photo only hits dataset-specific sense: LIDC "series/scan", RIGA "photographs" in data.md:15). OK but add "slice", "window", "view", "crop" definitions if not present (views are central to B1).

## Recommendations (numbered, concrete)

### Immediate fixes (small, do in the DOC/CLEAN phase)
1. Fix policy.py docstring to read "$10 alarm, hard stop HARD_STOP" without a number (refer to the constant). File: experiments/phase2/vast/policy.py:16-19. Close ADR-0009 action item #3 and delete the two "docstring still says $15" sentences in docs/runbooks/gpu-run.md:121 and ADR-0009:74.
2. Fix TRL version strings: Dockerfile.phase2:6 and train_sft.py:49 to 1.17.1, r2_qwen3vl2b.yaml:15 comment. Better: stop writing the version in comments at all; say "see pyproject.toml [train]". Add a test `tests/test_pins.py` that parses pyproject.toml and asserts every `TRL x.y.z` string in src/, experiments/phase2/*.py, docker/Dockerfile.* equals the pin (one regex, 15 lines).
3. Add the banner "Frozen record. TRL 1.16.0 and the $15 stop were true at the time; now see ADR-0005, ADR-0009" to experiments/phase2/PLAN.md, JUDGE_R1..R4.md, V0_V1.md (or move them, see R7).
4. Fix README.md layout block (docs/ line, add archive/ line stays) and add a link to docs/README.md at the top of the Quick start. Keep tests/test_readme.py; extend it (R6).
5. docs/testing.md: replace the hard-coded "1,484 / 1,453 / 31" with "run `pytest --collect-only -q`" plus a CI-printed number, or have `make docs-check` print and compare it (R6).
6. Move 2212.07020v1.pdf out of docs/research/ root into docs/research/papers/ or leave out of git (licence: ask owner).

### Final docs tree (proposal)
```
docs/
  README.md                       index by reader (keep; add a "status" column: Current / Planned / Frozen)
  getting-started.md
  contributing.md                 add archive rules, docs-update rule, phase-gate checklist
  glossary.md                     add: view, window, centred/shifted/rescaled, miss, ceiling (fair), medoid, RLOO, GED, HM-IoU, delta sync
  testing.md
  changelog.md
  architecture/
    overview.md                   CURRENT system only (as built)
    target.md                     NEW: the A-E target (moved from research2.0/Architecture.md, trimmed); each section tagged Built / Planned(phase)
    layers.md                     NEW (phase E): geom / verify / data / policy / rl / sft / infra, allowed imports, enforced by a test
    decisions/                    ADR-0001..0011 + 0012..
  reference/
    package.md  cuda-engine.md  reward.md  output-format.md  data.md  configs.md
    run-spec.md                   NEW (E): typed run spec fields replacing the YAML schema text in configs.md
    trl-coupling.md               NEW (A): private names touched, pinned version, failing test
    metrics.md                    NEW (A): ceiling, CIs, GED, HM-IoU, Dice definitions with code refs (split from guides/evaluation.md)
  guides/
    training.md  evaluation.md
    data-quality.md               NEW (B1): views, coverage audit, leakage checks, gates, how to read the report
    profiling.md                  NEW (A/C): NVTX ranges, nsys command, step budget, gate numbers
    cuda-kernels.md               NEW (D): how to add a kernel: source, test oracle, sanitizer run, benchmark, ADR
  runbooks/
    gpu-run.md  troubleshooting.md
    hf-checkpoints.md             NEW: push/pull, names only, no tokens
    ci.md                         NEW (D): tests.yml, image.yml, GPU CI job
  runs/
    INDEX.md                      NEW: run id -> config -> commit -> note -> HF tag -> $ (one table, appended per run)
  research/                       FROZEN evidence (tag: no edits)
    reports/ notes/ papers/ research2.0/ (rename "previous vs new.md" -> previous-vs-new.md)
    notes/ folder dirs without .md: merge into reports/ (they duplicate it)
  figures/                        add figures/README.md (which doc uses which, how generated)
```

### Final repo structure (proposal)
```
ObjectDelineation/
  README.md  ROADMAP.md (see R8)  CLAUDE.md  LICENSE  Makefile  pyproject.toml  CMakeLists.txt
  cuda/                 kernels + gpucheck.* ; add cuda/kernels/ per custom kernel (D)
  src/objdel/           layers (E): geom/ verify/ data/ policy/ rl/ sft/ infra/  (+ _gpucheck extension)
  tests/                unit/ integration/ gpu/ oracles/ goldens/   (rename test_p2_* by subject, drop "p2")
  experiments/
    core/               scripts still run (train_sft, train_rl, final_eval, preflight, hf_*, views, configs/, vast/, docker/)  <- from phase2/
    cuda_checker/       keep
    records/            FROZEN: phase0/, phase1/, phase2 plans/judges/impl/lab (moved, read-only, banner)
  docs/                 (tree above)
  archive/              retired CODE only: cpu_checker/ (exists), + CLEAN outputs: numpy_polygonizer/ (when D lands), old reward versions, LFM path
  runs/, data/          git-ignored
```
archive/ goes at the REPO ROOT (it already does). Rules: not built, not imported, not tested, excluded from ruff and pytest testpaths; each subfolder has a README saying what it was, what replaced it, the removal commit, and how to restore. Research records go to experiments/records/ or docs/research/, NOT archive/ (archive = code).
Compatibility note: moving experiments/phase2 breaks `python -I -m experiments.phase2.*` used by Makefile (CONFIG default), CLAUDE.md, docs, README, tests, vast bundle code_files.txt, Dockerfiles and the owner's memory. Do the move in ONE commit with a grep-and-fix and keep `experiments/phase2` only if the owner prefers; open question Q1.

### What each phase must add or update in docs
| Phase | Add | Update |
|---|---|---|
| A Measure | guides/profiling.md; reference/metrics.md (fair ceiling, nodule-level bootstrap CI, medoid selection, GED, HM-IoU, Dice); reference/trl-coupling.md; guide for "true outline-length count" (record the number, then ADR if the 1,536 limit changes); runs/INDEX.md; saved-answers format in reference/output-format.md or data.md | guides/evaluation.md (CIs are per nodule, not per row); testing.md (NVTX no-op on CPU); glossary; changelog |
| B1 Data | guides/data-quality.md (views, coverage audit, leakage checks, gates, thresholds); data.md (view schema, row counts) | guides/training.md (:106 flag row, :151 offset paragraph move to data-quality); reference/configs.md (rl_shifted.yaml); README results caveat that RL data changed |
| B2 Number match | ADR-0012 (match trainer and vLLM numerics: e.g. dtype/kernels chosen, IS mode) with before/after |ratio min=0 98.1 % of steps, per-token 0.018| | ADR-0005, ADR-0006; guides/training.md; troubleshooting.md |
| B3 Reward v5 | ADR-0013 (RLOO, no std division, location term) | reference/reward.md (full rewrite of the term table; REWARD_V5 fingerprint), ADR-0002 note "extended by", package.md `reward.*`, glossary |
| C GPU speed | guides/profiling.md step budget table; ADR-0014 (pixels on GPU / delta sync / vLLM awake / outline limit / KL 0, one ADR per hard-to-reverse choice only) | ADR-0006 (awake), ADR-0007, ADR-0011 (delta sync), architecture/overview.md CPU/GPU split table, README speed table, docs/figures/gpu_timeline.svg |
| D CUDA core | guides/cuda-kernels.md; ADR-0015 (CUDA polygonizer in production, ring tracing, run-based CCL); runbooks/ci.md (GPU CI job); testing.md (new oracle for EDT/GED/HM-IoU) | reference/cuda-engine.md (new entry points, memory, exactness); ADR-0001, ADR-0004; package.md `checker.polygonize` -> archived numpy path; figures/warp_ballot.svg |
| E Structure | architecture/layers.md; reference/run-spec.md; ADR-0016 (layer rules + typed run spec), ADR-0017 (TRL adapter boundary) | package.md (rewrite per layer), contributing.md (import rules), configs.md, getting-started.md (tour), docs/README.md |
| DOC | the tree above; `make docs-check` | everything in F2 |
| CLEAN | archive/<name>/README.md for each retired thing; changelog row | package.md, glossary rows "retired" (like polycore, glossary.md:110), contributing.md |
| After E (3D) | stub architecture/target.md section "3D voxels" only | ADR-0003 scope note |

### Doc quality gates
6. Add `tests/test_docs.py` (fast): (a) every relative link in docs/**/*.md (excluding docs/research) and README.md resolves, with `urllib.parse.unquote` and `<...>` support; (b) every backticked repo path exists; (c) every `make X` in docs is a Makefile target; (d) every ADR file is in decisions/README.md and has the nine template headings; (e) every config in experiments/**/configs/*.yaml is in reference/configs.md; (f) no "TRL 1.16" outside frozen records; (g) HARD_STOP and ALARM in docs equal policy.py constants (import policy, regex docs). These caught-nothing-today checks make F2 impossible to repeat. Add `make docs-check`, run in tests.yml `python` job.
7. Single source for numbers: docs quote only `runs/INDEX.md` or a named results file; put headline metrics in `docs/runs/INDEX.md` and have README link there. Add a test that README headline numbers equal the values in `runs/r5/eval/eval_test/final_eval.json` is not possible in CI (runs/ is git-ignored); instead copy the 8 headline numbers into a tracked `docs/runs/headline.json` and test README against it.
8. ROADMAP.md: either replace its body with a pointer to docs/research/research2.0/Direction.md plus a status table for A-E, or archive the old text to docs/research/ROADMAP-2026-10-18.md. CLAUDE.md and docs/README.md name ROADMAP.md as the plan, so keep the file and make it the short live plan. Update CLAUDE.md "Plan:" line accordingly (needs owner OK).
9. Budget docs: the $3 GPU budget (new accounting window). Document in runbooks/gpu-run.md and ADR-0009 addendum; code side (not DOC1): replace constants with `BUDGET_USD` per "window" read from a tracked file (e.g. experiments/core/vast/budget.yaml with window start date) and make the ledger sum only rows since the window. Rewrite test_vast hard-stop tests to parametrize on the budget.
10. Docs for experiment moves: add experiments/README.md (what is core, what is frozen, how to run as module) so the `-I -m` rule is stated in one place; remove duplicated statements in README, CLAUDE.md, getting-started.md.
11. Template for per-phase doc deliverable in contributing.md: "Definition of done for a phase = code + tests + `make test lint docs-check` + the doc rows from the table above + changelog row + runs/INDEX.md row."
12. Glossary additions (no implementation detail): view, window, centred/shifted/rescaled view, miss (IoU<0.1), fair ceiling, medoid, GED, HM-IoU, delta sync, RLOO. Keep one word per thing (image, outline).

## Open questions
Q1. Rename experiments/phase2 to experiments/core (breaks every `-m experiments.phase2.*` command and the vast bundle)? Or keep the names and only add experiments/README.md? Owner decision.
Q2. Is 2212.07020v1.pdf redistributable? If not, keep it untracked and cite by arXiv id only. [paper: arXiv 2212.07020, not read by me]
Q3. How are docs/figures/*.svg produced (hand-written or by a script)? If a script exists it must be tracked; if not, note in figures/README.md.
Q4. Should ROADMAP.md stay the live plan or yield to Direction.md? CLAUDE.md and README point at ROADMAP.md.
Q5. Dockerfile.phase2 and locks: does the image still install TRL 1.16.0 anywhere beyond the comment? I did not open the lock files (train.txt hash is rewritten by build_images.sh --lock).
Q6. Where do frozen JUDGE_*/impl/lab files go: banner in place (cheap, no path breakage) or move to experiments/records/ (clean, breaks links in docs and tests such as test_experiment_scripts.py)? I recommend banner first, move in CLEAN.
Q7. docs/changelog.md says 87 commits; verify when adding rows (not re-counted).
