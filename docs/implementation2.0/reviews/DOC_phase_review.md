# Review of Phase DOC (docs/implementation2.0/DOC.md, with DOC.spec.md)

Reviewer: fresh critical reviewer, 2026-11-08. No REPO file edited, no git write, no GPU, no spend. Probes ran only in
SP (`docs/implementation2.0/reviews/docrev/`, and the existing scratch tree `docs/implementation2.0/notes/doc/proj`, restored after each probe).

## Skills used

| Skill | What it changed in this review |
|---|---|
| nvidia-skills:nvidia-skill-finder | Checked the cached catalog `docs/implementation2.0/notes/doc/nv_catalog.json`: `compute-sanitizer-diagnose-cuda-bugs` exists (the DOC9 recommendation is real). No NVTX/Nsight docs skill, so no extra recommendation. |
| engineering:code-review | Correctness lens: every file:line the plan names was opened. Found the hf_push/RL gap (item 4), the dropped ruff `ignore` (item 17), stale `blob/v1.16.0` URLs (item 14). |
| engineering:debug | Reproduce-then-isolate: re-ran the Zensical strict build, seeded a broken link, a bad anchor and an omitted nav page to test what the build really catches (items 7, 8, 21). |
| engineering:system-design | Looked at where the docs live and how they flow (live plan placed under the frozen folder, item 3; card publication path, item 4). |
| engineering:architecture | ADR process: found that DOC2 edits the body of an accepted ADR, which DOC10's new rule forbids (item 5). |
| engineering:documentation | "Link, do not duplicate" and runbook rules: item 27 (CPU/GPU table has no first author), item 31 (strict build output). |
| huggingface-skills:huggingface-best | No model choice in DOC. Used only to check the HF card path against `hf_push.py` (item 4). |
| engineering:testing-strategy | Ran a prototype of drift checks 1-3 on today's tree: check 2 fails on day one (item 2); check 6 has false positives (item 1); gates need a recorded baseline (item 30). |
| engineering:tech-debt | Sequencing: docstrings for code that CLEAN deletes (item 16); unpinned transitive docs dependencies (item 13). |
| Web tools (WebSearch, WebFetch, PyPI JSON, GitHub API) | Confirmed: zensical 0.1.0 (2026-11-05), mkdocs-material 9.8.1 (2026-11-07), mkdocstrings-python 2.0.10, mkdocs-api-autonav 0.4.0, import-linter 2.16, lychee-action v2.9.0 (2026-07-09), CVE-2024-48908 fixed in 2.0.2; Zensical plugin page: mkdocstrings and api-autonav supported, gen-files not; lychee-action README needs `actions/cache` for `--cache`; lychee has `--offline` and `--include-fragments`. |

## Verified as correct (no action)

- 37 files under `docs/` outside `docs/research/`, 8 tracked; research2.0 0 of 4 tracked [measured: find / git ls-files].
- 63 missing docstrings (25 D102, 23 D103, 15 D101); 444 hits with all of `D` [measured: ruff 0.17.1].
- 8 "TRL 1.16" strings in the listed scope; pin `trl==1.17.1` [code: pyproject.toml:21; locks/train.txt:335].
- `policy.py` docstring says $15 [code: experiments/phase2/vast/policy.py:16,19]; constants 10.0 / 17.5 [code: :58].
- r5 test numbers: iou_tc 0.6448, f5 0.8164, precision 0.6946, recall 0.7478, dice 0.7007, legal 0.9938, items 2200,
  K 8 [measured: runs/r5/eval/eval_test/final_eval.json].
- Zensical 0.1.0 strict build of the trial tree: "No issues found", exit 0, 2.36 s; API pages generated, no
  `_gpucheck` page [measured: docs/implementation2.0/notes/doc/proj].
- All 16 config files are named in `docs/reference/configs.md` (check 5 passes today) [measured].
- ADR files 0001 and 0011 have the exact `**Status:**` / `**Date:**` / `**Deciders:**` lines and the six headings.

## Findings

1. **blocking** - DOC4 check 6 (`test_version_strings_match_pins`), spec B5.7.
   Problem: it compares every torch/transformers/trl mention with the `pyproject.toml` pin, but legitimate mentions
   name other environments. `experiments/phase2/docker/Dockerfile.vllm_native:4` says "torch 2.13.0 ...
   transformers 5.19.0" (the RL image: `locks/vllm_native.txt:838,878`); `Dockerfile.phase2:4` names the `/opt/vllm`
   venv (torch 2.13.0, transformers 5.19.0, `locks/vllm.txt:599,630`); `experiments/phase2/lab/B5/resolve_matrix.py:21-43`
   holds `trl==1.16.0` as experiment input data; `lab/B5/prereg.py`, `lab/B5/Dockerfile.phase2/phase3` and
   `experiments/phase2/impl/*.md` quote 1.16 as history. None carries the banner, and a Markdown banner cannot go
   into `.py` data. G3 fails on day one, or the agent rewrites frozen evidence.
   Fix: in `tests/docs_check.py` define `FROZEN_GLOBS = ["docs/research/**", "experiments/phase0/**",
   "experiments/phase1/**", "experiments/**/lab/**", "experiments/**/impl/**", "experiments/*/PLAN.md",
   "experiments/*/JUDGE_*.md", "experiments/*/V0_V1.md", "experiments/*/RESULTS*.md", "experiments/*/EXECUTION.md",
   "experiments/*/RESEARCH.md"]` and exempt by path, not by banner. Make the pin lock-aware with one map:
   `LOCK_FOR = {"experiments/phase2/docker/Dockerfile.vllm_native": "locks/vllm_native.txt"}` and, for
   `Dockerfile.phase2`, compare the `/opt/vllm` line with `locks/vllm.txt`; every other file uses `pyproject.toml`.
   Add a negative test: `TRL 1.16` in a tmp `src/x.py` is reported; `torch 2.13.0` in a tmp `Dockerfile.vllm_native`
   is not.

2. **blocking** - DOC4 check 2 (`test_backtick_paths_exist`) with scope `archive/**/README.md`, DOC8.
   Problem: an archive README must name retired paths. `archive/cpu_checker/README.md` names
   `src/objdel/_polycore.cpp`, `src/objdel/checker/polycore.py` and `experiments/prepare_data.py`, which do not exist
   [measured: prototype probe `docs/implementation2.0/reviews/docrev/probe.py`: 0 broken links, 3 missing backtick paths, all
   three in that file]. G3 fails on day one, and every future CLEAN folder fails the same way.
   Fix: drop `archive/**/README.md` from check 2 (keep it in checks 1 and 11), or skip spans under the
   "What it was" heading. Add a negative test that a retired path in a tmp archive README is not reported.

3. **major** - Section 1 item 5, DOC1, DOC3 nav, spec Constraints.
   Problem: the plan calls `docs/research/research2.0/Direction.md` the live plan, but puts it in
   "Research (frozen)", outside drift-check scope, under a spec rule "docs/research/** gets at most a one-line
   banner". DOC1 itself must edit `Research2.0.md:10` (it names `previous vs new.md`), which breaks that rule.
   Fix: in DOC1, before the owner's first `git add`, move the folder: `docs/research/research2.0/` ->
   `docs/plan/` (it is untracked, so a plain `mv`, no history lost), rename `previous vs new.md` ->
   `previous-vs-new.md`, update `docs/README.md:46` and `Research2.0.md:10`. Add `docs/plan/**` to `doc_files`
   (the probe found 0 broken links and 0 missing paths in these 4 files today [measured]). Nav section "Plan"
   above "Research (frozen)". Update G1 to "docs/plan: 4 of 4 tracked".

4. **blocking (for DOC11)** - DOC11 `hf_push.card(..., headline)`, owner command, Risks "card push".
   Problem: the headline numbers belong to the RL model `runs/r5/final_model`. `hf_push.main` needs
   `<run>/select/selected.json` and exits when no branch is selected [code: experiments/phase2/hf_push.py:106-108];
   `runs/r5` has no `select/` [measured: ls runs/r5]. `hf_push` also runs on the rented box [code: hf_push.py:1;
   policy.py:511-512 `cmd_push_hf`], so "re-publish the card" would need a rental, against the $3 budget.
   Fix: add `card_update(repo, revision, text)` to `hf_push.py` that uploads only `README.md` from the Mac:
   `HfApi(token=sys.stdin.readline().strip()).upload_file(path_or_fileobj=text.encode(), path_in_repo="README.md",
   repo_id=repo, revision=revision, commit_message="card: headline from docs/results/headline.json")`, token on
   stdin as today, `token_problems` checked first. Make `card()` accept an RL run (`sel=None` -> "RL from SFT
   branch X" line from the run config). CLI: `--card-only --revision REV --headline PATH`. Test: tmp run with no
   `select/` builds a card; a fake `HfApi` receives exactly one `README.md`. Owner runs the push.

5. **major** - DOC4 check 7 (`test_cost_constants_in_docs`), DOC2 ADR-0009 edits, DOC10 rule.
   Problem (a): history lines name old stops: `docs/runbooks/gpu-run.md:121` "The hard stop rose ...: $15, then
   $16.60, $17.00 and $17.50" and ADR-0009:72 "rose from $15 to $16.60 ...". "Every hard stop $X" is not a regex;
   a loose one fails on these lines. Problem (b): DOC2 deletes a sentence in ADR-0009 line 72 and ticks action item 3
   (line 81), while DOC10 adds "an accepted ADR changes only its Status".
   Fix (a): spec the regex exactly: `r"(?i)\b(hard stop|alarm)\s+(?:is|=|of|at)\s+\$(\d+(?:\.\d+)?)"`; history
   sentences stay legal. Negative test: "the hard stop is $15" in tmp docs is reported; "rose from $15" is not.
   Fix (b): word the DOC10 rule as "an accepted ADR changes only its Status line and its action-item check boxes";
   in ADR-0009 tick item 3 and leave the sentence (it is history), or supersede.

6. **major** - DOC7 ROADMAP rewrite, DOC9, owner rule "$3 budget" (BRIEF section 3).
   Problem: `Direction.md` section 4 quotes phase costs ("Phase A ... about $0.25", "Phase B ... $4.4-5.0")
   [code: Direction.md:106,123], and `policy.py` enforces 10.0 / 17.5 [code: policy.py:58]. A ROADMAP table copied
   from Direction would publish a cost above the owner's $3 budget, and no doc states the $3 rule.
   Fix: the ROADMAP table gets no cost column; one line under it: "GPU budget for phases A-E: $3 total on Vast AI,
   North America hosts, preflight before every rental (owner decision 2026-11-08)". Add the same line to
   `docs/runbooks/gpu-run.md` top. Put the mismatch (HARD_STOP 17.5 is cumulative Phase 2 spend, not the $3 budget)
   in the DOC phase report as an owner decision; the agent does not change the constant.

7. **major** - DOC3 step 3 "fix any file not in nav (Zensical lists them)".
   Problem: false. With `nav:` listing only 2 pages and `validation.nav.omitted_files: warn`, Zensical 0.1.0
   `build --strict -c` printed "No issues found", exit 0 [measured: docs/implementation2.0/notes/doc/proj, config restored]. With an
   explicit nav, every new phase page (DOC9 owner table) can vanish from the site with no failure.
   Fix: add check 13 `test_nav_complete` in `tests/test_docs.py`: every `docs/**/*.md` not excluded by
   `exclude_docs` and not under `reference/objdel/` appears in `mkdocs.yml` `nav` (read with
   `yaml.load(..., Loader=yaml.BaseLoader)`); negative test with a tmp page not in nav. Add it to G3 and spec B5.

8. **major** - DOC3 `mkdocs.yml`, DOC4 check 1 ("drop #anchor").
   Problem: no layer checks anchors, yet docs use them (`docs/runbooks/gpu-run.md:18` ->
   `../reference/data.md#7-bundle-experimentsphase2vastbundle`). A heading rename breaks it silently.
   Fix: add `validation: anchors: warn` to `mkdocs.yml`. Measured in SP: the tree still builds clean, and a seeded
   `[x](guides/training.md#no-such-anchor)` gives "1 issue found", exit 1. Re-check once with `-c` (my clean run
   used the cache). For the four files outside `docs/`, see item 21.

9. **major** - DOC6 `headline.json` schema, spec B7.
   Problem: the schema drops `view` (r5 is `lidc-test-shifted` [measured: final_eval.json `view`]) and the CI unit.
   Phase A switches to nodule-level bootstrap CIs and a fair ceiling, so a later file with different CI semantics
   would silently replace today's image-level CIs.
   Fix: schema adds `"view": "lidc-test-shifted"`, `"ckpt"`, `"seed"`, `"ci_unit": "image"` (phase A writes
   `"nodule"`), and `headline.py` refuses an input without `view`. README string gets "(shifted crops)" from `view`.
   Test: a `ci_unit` change makes `test_readme_headline_matches` fail until the README names it.

10. **major** - DOC1 publish step; `.gitignore` comment "repo public, owner 2026-11-06".
    Problem: DOC1 puts 33 untracked files into a public repo with no pre-publication check. My scan found no tokens,
    IPs or instance ids [measured: grep over the 33 files], but the plan has no such step, and RIGA-derived outputs
    must stay local (`.gitignore` RIGA block).
    Fix: DOC1 step 1b, run before the owner's `git add`: `git ls-files -z --others --exclude-standard docs |
    xargs -0 grep -nE "hf_[A-Za-z0-9]{20,}|ghp_[A-Za-z0-9]{20,}|ssh -p [0-9]+ root@|([0-9]{1,3}\.){3}[0-9]{1,3}"`
    must print nothing; list every figure or table derived from RIGA per-image results and get owner OK; the owner
    reads `git diff --cached --stat` before the commit.

11. **major** - DOC9 "owner table" row C and D; `architecture/overview.md` CPU/GPU table.
    Problem: the table is only "updated" by C and D; no task writes its first version. The owner rule is "CUDA for
    checker, scorer AND polygonizer"; today production uses the numpy polygonizer and builds pixel_values on the CPU
    (BRIEF section 2). Docs that say nothing about this read as "all on GPU".
    Fix: DOC9 "write now": the CPU/GPU table in `docs/architecture/overview.md` with today's facts and sources:
    polygonizer numpy [code: src/objdel/checker/polygonize.py], pixel_values fp32 built twice on CPU and rope
    positions in a Python loop [note: docs/implementation2.0/notes/B_gpu_audit.md], checker and scorer on GPU [ADR-0004]. Columns: part,
    device today, target device, phase.

12. **major** - DOC3 `docs` extra, "2 packages only".
    Problem: `zensical==0.1.0` and `mkdocstrings-python==2.0.10` pull unpinned mkdocs 1.6.1, mkdocstrings 1.0.8,
    griffelib 2.4.1, mkdocs-autorefs 1.4.5, pymdown-extensions 13.0, zensical_ui 0.0.41 [measured:
    docs/implementation2.0/notes/doc/.venv2 site-packages]. A transitive release can turn the strict build red with no repo change; the
    repo pins everything else (pyproject comments, docker locks).
    Fix: add `experiments/phase2/docker/locks/docs.txt` (or `docs/requirements.lock`) made by
    `uv pip compile pyproject.toml --extra docs -o <lock>`; CI and `make docs` install with `-c <lock>`; DOC4 check 6
    reads `zensical` from it.

13. **major** - DOC2 "Keep the URLs".
    Problem: the TRL source links still point at `blob/v1.16.0` in `src/objdel/train/configs.py:48`,
    `sample.py:30,33`, `rewards.py:50`, `callbacks.py:46`. Check 6 as specified (`TRL 1.16`, `trl==1.16.0`) does not
    match `v1.16.0` in a URL, so stale links survive G5.
    Fix: the names these docstrings cite exist in the installed TRL 1.17.1 (`_generate_and_score_completions`,
    `_calculate_rewards`, `def _generate(`, `RepeatSampler`, `prepare_multimodal_messages` found in
    `trl/trainer/grpo_trainer.py` / `trl/generation/` [measured: .venv trl 1.17.1]). Change the URLs to
    `blob/v1.17.1/`. Extend `version_mentions` to match `/blob/v(\d+\.\d+\.\d+)/` for trl and transformers.

14. **major** - DOC2 and DOC5 edit `src/objdel`; Risks table; order "DOC2-DOC12 beside phases A-E".
    Problem: the preflight stamp holds `objdel_and_script_sha256` and refuses a launch when it changes
    [code: experiments/phase2/preflight.py:69,79-87; train_sft.py:233]. A docstring merge between B1's preflight and
    its launch refuses the launch inside a paid or booked window.
    Fix: add a risk row "src text edit after preflight -> launch refused"; rule: DOC2 and DOC5 land before B1's
    preflight, or no `src/` edit from preflight to teardown. `docs/results/INDEX.md` records `git_head` and the code
    hash per run (both are in the run manifest).

15. **minor** - DOC5 sequencing with CLEAN.
    Problem: CLEAN removes dead code (CL1 lists LFM branches in `train/configs.py`, `checker/volume.py:62
    topology_any`); 11 of the 63 missing-docstring symbols are also flagged by vulture (60 % confidence, some are TRL
    callback hooks) [measured: overlap of ruff D101-D103 list and docs/implementation2.0/notes/vulture_src.txt].
    Fix: write the docstrings after CLEAN's removal list is final; enable the ruff gate in the same pull request.

16. **minor** - DOC5 `pyproject.toml` snippet.
    Problem: the snippet shows `[tool.ruff.lint]` with only `select`; the file also has `ignore = ["E731", "E741"]`
    and `[tool.ruff.lint.isort]` [code: pyproject.toml]. Copied literally, E731/E741 become errors and `make lint`
    fails.
    Fix: write it as a diff: add the six D codes to `select`; keep `ignore`; per-file-ignores become
    `"experiments/**" = ["E501", "B023", "B905", "D"]` and `"tests/**" = ["D"]`; add `[tool.ruff.lint.pydocstyle]`.

17. **minor** - DOC2 `policy.py` docstring vs G5 vs spec B9.
    Problem: DOC2 replaces "$10 alarm, $15 stop" and "passes $15" but line 19 also says "reaches $10". G5 checks only
    "$15"; B9 says "no dollar amount". The two gates disagree.
    Fix: replace "reaches $10" too; G5 text = B9 text: `re.search(r"\$\d", policy.__doc__) is None`.

18. **minor** - DOC3 CI job `docs`, fallback "install only `.[docs]` without `-e .`".
    Problem: the test step then has no pytest (`dev` extra) and cannot import `experiments.phase2.vast.policy`
    (repo root reaches `sys.path` only via the editable install [code: pyproject.toml scikit-build comment]).
    Fix: keep `uv pip install -e ".[dev,docs]"`; if slow, split: build step `.[docs]`, test step runs in the existing
    `python` job (it already runs `make test`, which collects `tests/test_docs.py`).

19. **minor** - DOC3 rollback to Material; G2.
    Problem: Material is built once in SP; no CI exercises the rollback, and its extra is a comment.
    Fix: weekly `links.yml` also runs `uv pip install mkdocs-material==9.8.1 mkdocs-api-autonav==0.4.0
    mkdocstrings-python==2.0.10 && NO_MKDOCS_2_WARNING=1 mkdocs build --strict` (non-blocking).

20. **minor** - DOC3 lychee job.
    Problem: `--cache --max-cache-age 1d` keeps `.lycheecache` only inside one runner; the lychee-action README
    restores and saves it with `actions/cache` [web: github.com/lycheeverse/lychee-action README].
    Fix: add `actions/cache/restore` and `actions/cache/save` (SHA-pinned, `path: .lycheecache`,
    `key: cache-lychee-${{ github.sha }}`, `restore-keys: cache-lychee-`, save `if: always()`); add `.lycheecache`
    and `/.cache/` (Zensical writes `.cache/` in the project root: 286 entries in the SP trial [measured]) to
    `.gitignore` in DOC1.

21. **minor** - DOC4 checks 1 (`relative_links`, `broken_links`); owner rule "prefer a good library".
    Problem: Zensical strict already fails on broken internal links in `docs/` (seeded `[x](does-not-exist.md)` ->
    exit 1 [measured]). The custom checker duplicates it for 29 of its 32 files.
    Fix: limit check 1 to the files outside `docs_dir` (`README.md`, `ROADMAP.md`, `experiments/README.md`,
    `archive/**/README.md`), or replace it with `lychee --offline --include-fragments` on those files in the `docs`
    CI job and a `make links-local` target. Keep the negative test.

22. **minor** - DOC7 "phases A to I from `Direction.md` section 3".
    Problem: the phases are in section 4 ("## 4. The phases", line 104); section 3 is "The plan in one picture"
    [code: docs/research/research2.0/Direction.md:80,104]. CAD is "Parked" between H and I.
    Fix: cite section 4; add CAD as one row with status "Parked (owner 2026-11-08)", no gate.

23. **minor** - DOC3 step 4 subpackage list.
    Problem: the build also makes `reference/objdel/secrets/` (module `src/objdel/secrets.py`) [measured]. The step
    as written fails its own check.
    Fix: list `secrets` and check that the page shows no secret value (it documents keychain names only).

24. **minor** - DOC6 rounding.
    Problem: means are rounded to 4 places, then formatted with `.3f`. Double rounding can move the third digit
    (for example 0.74449 -> 0.7445 -> "0.745", direct "0.744").
    Fix: keep full precision in `headline.json`; round only inside `headline_strings`.

25. **minor** - Goal "Numbers in docs come from one tracked results file"; README validation numbers.
    Problem: README line 121 also quotes validation (IoU 0.618, F 0.806, precision 0.689, recall 0.716); those stay
    hand-copied. `runs/r5/eval/eval_val` exists [measured].
    Fix: `headline.json` gets `"test"` and `"val"` blocks; `headline_strings` returns both sets; or narrow the goal to
    "test headline".

26. **minor** - DOC6 file name `docs/results/INDEX.md`.
    Problem: MkDocs/Zensical treat only `index.md` and `README.md` as folder index pages; `INDEX.md` becomes
    `/results/INDEX/`, and on the Mac's case-insensitive file system it collides with any later `index.md`.
    Fix: name it `docs/results/README.md` (or `runs.md`) and update DOC7, DOC9, DOC11 and the definition of done.

27. **minor** - Section 3.3 "about 120 Markdown files".
    Problem: unverified. `docs/` has 98 `.md` files, 29 outside `docs/research/`; check scope is about 32
    [measured: find].
    Fix: write "32 files in scope (98 under docs/)" with the command.

28. **minor** - DOC3 `Makefile`.
    Problem: `.PHONY` lists `setup build test test-all cov lint fmt wheels cuda preflight` only [code: Makefile:6].
    Fix: add `docs docs-serve headline` to `.PHONY`; `headline` uses `$(PY) -I -m experiments.phase2.headline
    --run $(RUN)` and fails when `RUN` is empty.

29. **minor** - G7 "collected test count = before + new docs tests".
    Problem: "before" is not recorded; `docs/testing.md:5` is the only number and DOC4 deletes it.
    Fix: record `.venv/bin/python -m pytest --collect-only -q | tail -1` in the DOC phase report before DOC1.

30. **minor** - DOC9 `docs/runbooks/ci.md`.
    Problem: a Zensical strict failure ends in a Python traceback ("RuntimeError: Aborted because --strict flag is
    set"); the real issue is printed in a box above it [measured]. A reader stops at the traceback.
    Fix: one line in `ci.md`: "On a docs failure, scroll up to the box that says `N issue(s) found`."

## Verdict

Request changes. Fix items 1, 2 and 4 before DOC4 and DOC11 start, and settle item 3 (where the live plan lives)
in DOC1, because DOC1 decides the paths every later check reads. The tool choice (Zensical with a Material rollback)
and most numbers hold up under re-measurement.
