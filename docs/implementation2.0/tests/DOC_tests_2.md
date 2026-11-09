# DOC phase tests, set 2 (adversarial, property, performance, failure modes)

## Skills used
- engineering:testing-strategy: shaped the plan by level (unit, property, integration, e2e) and a coverage target per area.
- mattpocock-skills:tdd: tests sit only at the seams the spec names (`tests/docs_check.py` pure functions, `headline.main`,
  `hf_push.card_update`, `make` targets); expected values are literals, never recomputed the way the code computes them.
- Not invoked (no effect on a docs-only test plan): the other skills in the task list. This file was written from DOC.spec.md
  and the repo's tests/ style only (`tests/test_p2_hf_push.py`, `tests/conftest.py`), with zero design context.

Conventions. REPO = repo root. Level "GPU" does not occur: the spec forbids GPU use, so GPU/CPU bit-equality applies
only as "G12: nothing in the checker/golden path changed" (T30-T31). All tests are stdlib + pytest + pyyaml +
hypothesis (already in the dev extra; allowed because it is already pinned in dev; if the owner reads "no new test
dependency" as "not even hypothesis in these files", use the seeded `random.Random(0)` loop variants noted in each
property test). Fixture tree = `tmp_path` with a fake repo (helper `mkrepo(tmp_path, **files)` writes files from a dict).
Numbering: T1..T40. Priority sketches for T1, T9, T14, T21, T27 are at the end.

---------------------------------------------------------------------------------------------------------------------
## A. Link, path and Make-target checks (B5.1, B5.2, B5.3)

**T1. Link extractor: hostile Markdown shapes** (property + unit; docs_check.relative_links)
- Setup: none.
- Input (each a separate parametrize case) and expected list:
  - `[a](x.md)` -> `["x.md"]`
  - `[a](<my file.md>)` -> `["my file.md"]` (angle form kept, brackets stripped)
  - `[a](my%20file.md)` -> `["my%20file.md"]` or decoded; broken_links must treat both as the file `my file.md` (see T2)
  - `[a](x.md#sec)` -> link is reported but its anchor part is dropped by broken_links (T2)
  - `[a](#only-anchor)`, `[a](https://e.com/x.md)`, `[a](http://e.com)`, `[a](mailto:a@b.c)` -> `[]`
  - `![img](fig.png)` -> `["fig.png"]` (images are links too; a missing figure is a real bug)
  - A link inside a fenced block ```` ``` ```` or inline code `` `[a](nope.md)` `` -> `[]` (code is not a link)
  - `[a](x.md "title")` -> `["x.md"]` (title dropped)
  - `[a](x(1).md)` (balanced parentheses) -> `["x(1).md"]`
  - Two links on one line, one reference-style `[a][ref]` with `[ref]: y.md` -> includes `y.md`
  - `[a]()` (empty target) -> `[]` and no exception
  - `[a](HTTPS://E.COM)` (upper-case scheme) -> `[]`
  - `[a](//cdn.example.com/x)` (scheme-relative) -> `[]`
  - `[a](/abs/path.md)` (root-absolute) -> expected: reported (a doc-scope file must not rely on the host root); the test fixes the
    chosen behaviour so it cannot flip silently.
- Why: a regex that treats code spans as links gives false failures that teach the owner to ignore the check; one that skips
  images or reference links lets a broken figure through.

**T2. broken_links: path resolution edge cases** (unit; tmp_path)
- Setup: tree `README.md`, `docs/a.md`, `docs/sub/b.md`, `docs/my file.md`, `docs/x(1).md`.
- Input and expected:
  - `docs/sub/b.md` links `../a.md` -> `[]`.
  - `docs/a.md` links `sub/b.md#heading` -> `[]` (anchor dropped, file exists).
  - `docs/a.md` links `my%20file.md` and `<my file.md>` -> `[]`.
  - `docs/a.md` links `missing.md` -> `["missing.md"]` (exact string, so a message can name it).
  - `docs/a.md` links `sub` (a directory) -> expected: reported unless the directory has `index.md`/`README.md`; test fixes behaviour.
  - Link escapes root: `docs/a.md` links `../../etc/passwd` -> reported (never `exists()`-true outside root is not asserted;
    the function must not follow it out of `root`). Create `tmp_path.parent/"outside.md"` and link `../../outside.md`: reported.
  - Case: on a case-insensitive file system (Mac) `Docs/A.md` exists for `docs/a.md`. A link `A.md` to `a.md` must be REPORTED
    (CI is Linux and fails it). Skip the case if the tmp file system is case-sensitive but assert on case-insensitive: use
    `os.listdir` name comparison. This is the classic "passes on the Mac, fails in CI" bug.
  - Symlink loop `docs/loop -> docs`; link `loop/loop/loop/a.md` -> no hang (timeout 5 s), result is either [] or reported, never an exception.
  - Link with trailing `?query` (`a.md?plain=1`) -> treated as `a.md`.
  - CRLF file and a file with a UTF-8 BOM: links still found.
  - A file that is not valid UTF-8 -> raises a clear error naming the file, or is reported; never a silent pass.
- Why: each case is a distinct way the check can say "fine" for a broken link or crash the whole suite.

**T3. Doc-scope selection (doc_files, is_frozen, outside_docs_dir)** (unit + property)
- Setup: tree with one file at every frozen glob and one near-miss for each:
  `docs/research/x.md` (frozen) vs `docs/researchX/x.md` (doc scope); `experiments/phase0/a.md` (frozen) vs
  `experiments/phase05/a.md`; `experiments/phase2/PLAN.md` frozen vs `experiments/phase2/PLAN_NOTES.md` frozen (glob `PLAN*.md`)
  vs `experiments/phase2/MYPLAN.md` (doc scope); `experiments/phase2/lab/deep/x.md` frozen; `experiments/phase2/impl/a.md` frozen;
  `experiments/cuda_checker/README.md` frozen; `experiments/README.md` doc scope; `archive/foo/README.md` doc scope;
  `archive/foo/NOTES.md` NOT doc scope (only README); `docs/plan/Direction.md` doc scope (plan is not research);
  `node_modules/x/README.md`, `.venv/lib/x/README.md`, `build/x.md` not doc scope.
- Expected: `is_frozen` True/False exactly per the comments; `doc_files(root)` returns sorted, de-duplicated, only `.md`;
  `outside_docs_dir` returns exactly those not under `docs/` (README.md, ROADMAP.md, experiments/README.md, archive READMEs).
- Property: for 200 random relative paths from a small alphabet of segments, `is_frozen(p) == any(fnmatch-like match of FROZEN_GLOBS)`
  AND `is_frozen` is invariant to a leading `./` and to Windows separators; `doc_files` never contains a frozen path.
- Why: a too-wide frozen glob silently removes a file from every check; a too-narrow one fails the build on files the owner must not edit.

**T4. backtick_paths: what counts as a path** (unit/property)
- Input and expected:
  - `` `src/objdel/checker` `` -> reported name; `` `src/objdel/` `` (trailing slash) -> reported (existing dir passes in broken-check).
  - `` `src/objdel/x.py:42` `` (line suffix) and `` `src/objdel/x.py::fn` `` -> spec says span names a path: test fixes that a `:NN`
    suffix is stripped before the existence test (otherwise every code reference fails). Expect 0 failures when file exists.
  - `` `experiments/<phase>/x.py` ``, `` `docs/*.md` ``, `` `src/{a,b}` ``, `` `docs/a b.md` `` -> not returned (space, `*`, `<`, `{`).
  - `` `docs` `` (no slash), `` `./docs/a.md` ``, `` `/docs/a.md` ``, `` `my/docs/a.md` `` -> not returned (must START with a prefix + `/`).
  - `` `tests/` `` alone and `` `archive/` `` alone -> returned and exist.
  - `` `docs/a.md` `` inside a fenced code block -> not returned (shell transcripts are not claims).
  - Double-backtick span ``` ``docs/a.md`` ``` -> returned.
  - Two spans on one line; a span that crosses a newline -> not returned.
  - `` `src/objdel/gone.py` `` in `archive/x/README.md` -> still returned by the function; exemption is in the test (B6 legal case),
    so assert the exemption by running the full check on that tree and getting 0 reports.
- Why: path claims in prose are the main doc-rot signal; the three false-positive families (globs, placeholders, line suffixes)
  decide whether the owner keeps the check on.

**T5. make_targets / make_mentions: Makefile corner cases** (unit)
- Setup Makefile text: `.PHONY: a b\n\na: ## doc\n\tx\n\nb c: d\n\nVAR := x\n$(EXT_STAMP): $(EXT_SRC)\n%.o: %.c\n.SUFFIXES:\nCONFIG ?= foo\n\ttab-indented line: with colon\n`.
- Expected `make_targets`: `{"a","b","c"}` (plus whatever the code decides for `.PHONY`; fix: dot-targets and `%`/`$(` patterns
  EXCLUDED; variables excluded; tab-indented lines excluded; `.PHONY: a b` alone is not enough to define a target, but a target
  listed only in `.PHONY` should still count only if rule exists: fix behaviour).
- `make_mentions` input and expected:
  - "run `make docs` now" -> `{"docs"}`; "make docs-serve" -> `{"docs-serve"}` (hyphen kept).
  - "make sure", "make it", "make a copy", "to make test cases" -> must produce NO false mention for plain English. Decide and fix:
    only `` `make X` `` code spans and lines starting `make X`/`$ make X` count. "make test" in running prose without backticks: not counted.
  - "`make headline RUN=r5 MODEL_REF=x`" -> `{"headline"}` (variables are not targets).
  - "`make -C sub docs`", "`make -j4 test`" -> flags skipped, `{"docs"}` / `{"test"}`.
  - "`make test && make lint`" -> `{"test","lint"}`.
  - "`make`" alone -> `set()` (default target, not a mention).
  - "`make docs`" inside a table cell and inside fenced block -> counted (commands in fences ARE commands; document this choice).
  - Upper-case "`Make docs`" -> not counted.
- Why: both directions matter: false positives on English prose make the check noisy; missing flags makes it blind.

**T6. Make targets exist and behave** (integration; `make -n`, no execution)
- Setup: REPO root. Run `make -n docs`, `make -n docs-serve`, `make -n headline RUN=r5 MODEL_REF=x`.
- Expected: exit 0 each; output of `-n docs` contains `zensical build --strict`; `make -n headline` (RUN empty) exits non-zero with
  "RUN" in stderr and prints no `python` line; `make -pn | grep '^.PHONY'` contains `docs`, `docs-serve`, `headline`.
- Failure variant: `make -n headline RUN= MODEL_REF=x` and `RUN="r 5"` (space) and `RUN=../x` -> second and third must not reach
  python with an unquoted value (`-n` output shows RUN quoted or the make fails). `RUN='a; rm -rf /'`: output must never contain a bare `; rm`.
- Why: spec only says "fails when RUN is empty"; shell injection through a Make variable is the adversarial case.

---------------------------------------------------------------------------------------------------------------------
## B. Version, cost, citation, ADR, config, nav checks (B5.4-B5.13)

**T7. version_mentions: positives and the exempt units** (unit/property)
- For name in {trl, transformers, torch, vllm}. Expected returns the version string; None for exempt.
  - "trl==1.17.1" -> ["1.17.1"]; "TRL 1.16" -> ["1.16"]; "TRL  1.16.0" (two spaces) -> ["1.16.0"];
    "trl 1.17.1." (sentence period) -> ["1.17.1"] (trailing dot stripped); "torch==2.14.1+cu130" -> ["2.14.1"] (local tag dropped).
  - ".../huggingface/trl/blob/v1.16.0/x.py" -> ["1.16.0"] for trl only; for torch -> [].
  - "vLLM 0.09 s", "vllm 0.5 ms", "torch 3 GB", "torch 100 %", "trl 2 x", "vllm 0.5s" (no space) -> []. "vllm 0.5 seconds" -> decide:
    unit list is `s, ms, GB, MB, %, x`; "seconds" is NOT in the list, so it is a mention; the test pins this.
  - "trl 1.17.1 s" -> [] ? The unit rule says a unit after the version exempts; pin it.
  - "the trl library", "trlx 1.0", "ptrl 1.2", "atorch 2.0" -> [] (word boundary).
  - "TRL: 1.16" (colon) -> [] (not "name then whitespace/==").
  - "transformers 5" (no dot) -> [] or ["5"]; spec: two-part mention is legal if prefix matches. Pin: one-part mention is ignored.
  - "torch>=2.14" -> [] (spec lists only `==` and whitespace).
  - Mixed case in one line: "trl==1.17.1 and transformers==5.20.1" with name=trl -> only trl's.
- Property (200 random strings from the grammar `name sep version [unit]`): result non-empty iff no unit and sep in {`==`, space}.
- Why: version drift is the whole point; an over-eager unit rule silently exempts a real mention ("trl 1.16 x" typo), an
  under-eager one fails on "vLLM 0.09 s" timings that exist in research notes.

**T8. allowed_versions and the two-part rule** (unit; tmp_path)
- Setup: `pyproject.toml` with `"trl==1.17.1",`, `"torch==2.14.1"`, `"ruff==0.17.1"`, `"pytest>=8.0"`, `"zensical==0.1.0"`;
  `locks/train.txt` with `trl==1.17.1`, `torch==2.14.1+cu130`, `# via trl==9.9.9` (comment), and `vllm.txt` with `vllm==0.21.0 \`
  (hash continuation). `vllm_native.override` is NOT `*.txt` so it is not read.
- Expected: `allowed_versions(root,"trl") == {"1.17.1"}` (comment line ignored); torch == {"2.14.1"} or {"2.14.1+cu130","2.14.1"}
  (pin the choice; a mention "torch==2.14.1" must pass); vllm == {"0.21.0"}; `>=` pins never allowed; name match is exact and
  case-insensitive with `_`/`-` normalised (`Torch`).
- Check 7 two-part rule: allowed {"1.17.1"}: "trl 1.17" legal, "trl 1.1" ILLEGAL (prefix "1.1." does not start "1.17.1"),
  "trl 1.1" vs allowed "1.17.1": must fail. "trl 1.17.10" illegal (not a prefix match on a dot boundary).
- Failure mode: no pin for the package anywhere -> allowed set empty; every mention reported (not a crash, not a silent pass).
- Why: "1.1" matching "1.17" by string prefix is the classic off-by-boundary bug.

**T9. Version check on the real repo: today's counts and the regression numbers** (integration; regression of today's numbers)
- Run `version_mentions` for trl across doc scope + the listed code globs on the REAL repo before the fix commit and after.
- Expected before DOC fixes (record, not asserted forever): 8 "TRL 1.16" strings and 5 `blob/v1.16.0` links outside frozen
  paths (spec G5). Expected after: 0 for both; `grep -rn "TRL 1.16" . --include='*.md' ` outside frozen paths returns nothing.
- Also assert in the test: every mention found in a frozen path (`docs/research/**`) is NOT in the check's report (skip works), and
  the count of such frozen mentions is > 0 (proves the skip is exercised, not that the files vanished).
- Why: tests that only run on seeded input cannot tell that the real repo is clean; a count of frozen mentions > 0 guards a
  vacuous pass (the check scanning nothing).

**T10. cost_mentions and the policy tie** (unit + integration)
- Input and expected:
  - "| Alarm | $10 |" -> [("alarm", 10.0)]; "| Hard stop | $17.50 |" -> [("hard stop", 17.5)];
  - "| **Alarm** | $10 |" (bold) and "| alarm | $10 |" (case) -> decide; pin: both match.
  - "| Hard stop | $1,000 |" -> [("hard stop", 1000.0)] (comma); "| Alarm | $10.00 |" -> 10.0 (not "10.0.0"); "| Alarm | USD 10 |" -> [] (no `$`)
    and in check 8 must be flagged only for the runbook row presence, not value.
  - "The hard stop is $17.5." (sentence period) -> 17.5, not "17.5." parse error.
  - "alarm at $10" , "Hard stop = $17.50", "hard stop of $17.5", "hard stop is $17.5", -> one tuple each.
  - History: "rose from $15 to $17.5" -> only if the sentence is "hard stop is $17.5"; "hard stop rose from $15" -> [].
  - "| Budget | $3 |", "| Alarm clock | $10 |" ("Alarm clock": first cell is not exactly "Alarm") -> [].
  - Row with two dollar values "| Alarm | $10 (was $7) |" -> first value 10.0 only, or both; pin one and require the policy compare to pass.
  - `$` in LaTeX or shell (`$PATH`, `$1`) -> no tuple.
- Integration: real `docs/runbooks/gpu-run.md` has exactly one Alarm row and one Hard stop row; both values equal
  `policy.ALARM`, `policy.HARD_STOP` (today 10.0 and 17.5, [code: experiments/phase2/vast/policy.py]).
- Mutation test: monkeypatch `policy.HARD_STOP = 17.0` and run check 8 on the real docs -> reports the runbook (and every other doc row).
  Also set `POLICY_MODULE` to a nonexistent module -> check fails with a clear error (not a silent pass, not a skip).
- Why: a doc that disagrees with the $ stop is a money bug; the check must not be satisfiable by deleting the row (spec requires
  presence, T10 deletes the Hard stop row in a copy and expects a report).

**T11. Policy docstring has no dollar amount (B9)** (unit)
- `re.search(r"\$\d", policy.__doc__) is None`. Mutation: append "$10 alarm" to a copy of the docstring -> regex finds it.
- Also: `policy.ALARM` and `policy.HARD_STOP` unchanged (10.0 and 17.5): constants are owned by code, [code: policy.py ALARM, HARD_STOP].
- Why: the docstring is the second source of truth that went stale ("$10 alarm, $15 stop" while the constant is 17.5, observed today).

**T12. ADR checks (B5.4, B5.5)** (unit; tmp_path)
- Setup per case: `decisions/0001-a.md`, `0002-b.md`, `_template.md`, `README.md` with a table.
- Cases:
  - File `0003-c.md` not in the index -> reported; index row `0004-d.md` with no file -> reported; both at once -> two reports.
  - `_template.md` absent from index -> NOT reported; `README.md` is never an ADR; `0001-a.MD` (upper ext) and `001-a.md` (3 digits)
    and `0001_a.md` (underscore) -> not ADR files by name pattern (`adr_files` ignores them); decide whether a malformed name is itself a
    report (recommend yes: a "NNNN-" typo hides an ADR). Pin it.
  - Index links in all shapes: `[a](0001-a.md)`, `[a](./0001-a.md)`, `[a](0001-a.md#x)`, bare text `0001-a.md` -> `adr_index` finds each.
    A row for the same file twice -> set (no duplicate report).
  - Headings: all six present -> pass; `## Context ` (trailing space) -> pass; `### Context` (level 3) -> FAIL; `## context` (lower case) -> FAIL;
    a heading only inside a fenced block -> FAIL (not a real heading); `**Status:** Accepted` indented by 2 spaces -> FAIL (must start the line);
    `**Status:**` present twice -> pass; missing `**Deciders:**` -> report names exactly the missing item.
  - ADR 0012 (new) must itself pass all of check 5, and `_template.md` must contain every required heading (so a new ADR copied from it passes).
- Why: ADR index drift is silent; the template is the source of future ADRs, so it is tested as a producer.

**T13. Config check (B5.6) and nav check (B5.13)** (unit + integration)
- Config: `experiments/phase2/configs/a.yaml`, `experiments/phase3/configs/b.yml` (extension .yml: spec says `*.yaml`, pin that
  `.yml` is not collected), `experiments/phase2/configs/sub/c.yaml` (nested: `configs/*.yaml` is not recursive -> not collected),
  and `docs/reference/configs.md` containing `a.yaml` only in a link to `../x/a.yaml.bak` (substring trap: must be reported, since the
  name is `a.yaml.bak`), then with `a.yaml` in backticks (passes). Name mentioned as a substring of another name
  (`ra.yaml` mentions `a.yaml`) -> must be reported. Two configs with the same file name in different phases: one mention of the bare
  name does not cover both -> require the path, or the bare name only if unique. Pin it.
- Real repo: all 16 files listed in `experiments/phase2/configs/*.yaml` today are named in `docs/reference/configs.md`.
- Nav: `nav_pages` on a mkdocs.yml with nested lists, `- Title: file.md`, `- file.md` (bare), `- Title: <url>` (external, ignored),
  and a section with a `index.md`. Pages excluded by `exclude_docs` (gitignore-style: `*.pdf`, `research/`, `!research/keep.md`)
  must not be required. A page under `docs/reference/objdel/` is not required. Add a docs page -> reported; remove it from nav -> reported;
  rename in nav with wrong case -> reported. YAML with the Material `!!python/name:` tags (mkdocs.yml may contain
  `emoji_index: !!python/name:material.extensions.emoji.twemoji`) must be parsed with a loader that ignores tags, or `nav_pages` crashes: test with such a tag present.
- Why: the tag crash only appears when the extension keys get added; the substring trap makes a stale config pass.

**T14. Negative-test completeness (B6 itself)** (meta; the most important one against a vacuous suite)
- For each of checks 1,2,3,4,7,8,11,13 there is a seeded-bad `tmp_path` repo; the test asserts `len(report) == expected` (not just
  `> 0`) and that the message contains the offending path/name. Each legal case in B6 is asserted to yield 0 reports on the same
  seeded tree (one tree contains BOTH the bad and the legal inputs, so a check that over-reports or under-reports is caught in one run).
- Mutation harness (in the test): for each check function in the module, monkeypatch its core predicate to `lambda *a, **k: True`/`[]`
  and assert the corresponding negative test fails. (Implement as: replace the function with a no-op returning an empty list; run
  the seeded test body; expect AssertionError.) Guarantees no negative test is tautological.
- Why: the spec's biggest risk is a check that always passes; this makes the passing state falsifiable.

**T15. Archive README check (B5.11)** (unit)
- Cases: `archive/old/README.md` with all four parts (what it was, what replaced it, the retiring commit, how to restore) -> pass.
  Missing "retiring commit" -> reported naming the missing part. `archive/old/` with no README -> reported. `archive/README.md`
  (the index) is not a subfolder README: not required to contain the four parts. `archive/.hidden/` -> ignored or required (pin).
  A commit SHA that is not 7-40 hex characters ("TBD", "abc") -> reported. A nested `archive/old/sub/` without README -> not required (direct subfolders only).
  Part detection must not be satisfied by the words appearing in a code fence or one shared heading `## Restore` containing the
  text "TODO" -> reported when the section body is empty (heading with no text below it).
- Why: empty headings pass naive substring checks.

**T16. Attribution strings (B5.12) and CITATION (B5.10)** (unit)
- Files: PR template, contributing, CITATION.cff. Cases: "Co-Authored-By: Claude" (exact), "co-authored-by: claude" (case), "Co-authored-by:  Claude <noreply@anthropic.com>"
  (extra space), "Generated with [Claude Code]" -> all reported? Spec literal strings are `Co-Authored-By: Claude` and `Generated with Claude`;
  the third and fourth are the adversarial inputs: pin that the check is case-insensitive and tolerates Markdown (`[Claude Code]`, `🤖 Generated with`).
  Clean text: "claude" used in a sentence about the model Qwen-vs-Claude -> not reported.
- CITATION: `version: 0.1.0` vs pyproject `version = "0.1.0"` pass; quoted YAML `version: "0.1.0"` pass; `version: 0.1` fail; `version: 0.1.0.dev1` fail;
  CITATION missing `version` key -> reported (KeyError is a failure, not a pass); pyproject with dynamic version -> clear error.
  A CITATION that parses as YAML float (`version: 0.10` -> 0.1 float) must compare as text "0.10": the test reads it as string.
- Why: YAML float coercion makes 0.10 equal 0.1, a real trap.

---------------------------------------------------------------------------------------------------------------------
## C. Headline file, `headline.py`, strings (B5.9, B7, G6)

**T17. headline_strings: exact output, order and rounding** (unit)
- Input: dict with test iou_tc mean 0.64549999, f5 0.8159, precision_vote 0.6954, recall_vote 0.7476, dice_vote 0.70051, legal 0.99357; val
  iou_tc 0.6184, f5 0.8063, precision 0.6891, recall 0.7164; `ci_unit="item"`, test view `lidc-test-shifted`.
- Expected exactly (list equality, order matters):
  `["IoU 0.645","boundary F 0.816","precision 0.695","recall 0.748","Dice 0.701","legal 99.4 %","shifted crops","IoU 0.618","F 0.806","precision 0.689","recall 0.716"]`
  (0.64549999 -> "0.645"; 0.99357 -> "99.4"). Length 11 for `ci_unit="item"`, 12 with `ci_unit="nodule"` ending `"nodule-level"`.
- Rounding adversarial: mean 0.0005 -> "0.001"? (Python `format` rounds half-even on the binary value: 0.0005 -> "0.001" or "0.000"; pin the actual
  Python output with a literal computed offline and written in the test, not recomputed.) mean 1.0 -> "1.000"; mean 0 -> "0.000"; legal 1.0 -> "100.0 %".
- View words: `lidc-test-shifted` -> "shifted crops"; test with each other known view from `experiments/phase2/build_eval_views.py` outputs
  (read the real names with glob on `data/eval_views`? NO: tests must not read data; hard-code `lidc-val-shifted`, `lidc-test-centred`, an unknown name `foo-bar`).
  Unknown view -> pin: raises ValueError with the view name (not silently a wrong sentence).
- Missing key (`dice_vote` absent) -> KeyError naming the key; extra metric key -> ignored.
- Why: G6 numbers in the README are the project's public claim; one swapped key or a rounding mode difference ships a wrong number.

**T18. README contains every headline string, and is checked after rounding trap** (integration; real repo)
- Read `docs/results/headline.json` and `README.md`. Expected: every string from T17's function is a substring of README.
  For the r5 values: all 11 literal strings from G6 appear. Negative: mutate one number in a copy of headline.json (0.645 -> 0.646) -> check 9 reports exactly "IoU 0.646".
- Adversarial: README contains "IoU 0.6450" (more digits) -> "IoU 0.645" is a substring so it passes; pin that this is accepted OR require a word boundary
  (recommended: boundary regex `(?<![\d.])IoU 0\.645(?![\d])` so "IoU 0.6451" fails). Add the boundary case "IoU 0.6451" -> reported.
- "IoU 0.645" appearing only in the VAL position while the test number differs: not detected by substrings (documented limitation); add one test with the
  test row first and the val row second to make the check at least order-independent but present.
- Why: substring checks pass on the wrong table row; boundary guard prevents prefix-collision.

**T19. headline.py refusal matrix** (integration; tmp_path, `main(argv)`)
- Setup helper writes `runs/ID/eval/eval_test/final_eval.json` and `eval_val/final_eval.json` from a minimal valid fixture
  (6 metrics each with `mean` and `ci95`, `split`, `view`, `items`, `k`, `seed`); `--out tmp_path/"h.json"`; run with cwd=tmp_path or
  a `--runs` root if present (use monkeypatch.chdir).
- Refusals (each: `main(...)` returns non-zero, `out` file does not exist, no `.tmp` sibling left behind, and stderr/exit text names the cause):
  1. test file missing; 2. val file missing; 3. both missing; 4. test file `split == "val"` (swapped); 5. val file `split == "test"`;
  6. test file lacks `view`; 7. `view` is `""`; 8. `view` is `null`; 9. file is empty; 10. file is truncated JSON; 11. file is a JSON list not object;
  12. a metric key missing (e.g. no `dice_vote`); 13. `ci95` has one element or NaN; 14. `mean` is a string "0.645"; 15. `--run ""`; 16. `--run ../x`
  and `--run a/b` (path traversal -> refuse); 17. `--model-ref` empty; 18. `final_eval.json` is a symlink to outside runs -> refuse or document.
- Failure atomicity: when out exists from an earlier run and a new run is refused, the old file is byte-identical afterwards.
- Why: a half-written or wrong-split headline file becomes the "one source of truth" for the README and the HF card.

**T20. headline.py success: exact copy, schema and determinism** (integration + property)
- Input: fixture with means such as `0.6454999999999999` and CI `[0.6312345678901234, 0.6591234567890123]`.
- Expected: output JSON parses; keys exactly `run, model_ref, ckpt, seed, K, ci_unit, test, val`; each split has exactly keys `source, view, items, metrics`;
  metrics keys exactly the 6 (extra input metric `foo` is dropped, `n_samples` too); `mean` and `ci95` floats are `==` to input (compare parsed floats with `==`,
  and compare text with `repr`; unrounded); `view` is last path part: input `data/eval_views/lidc-test-shifted/` (trailing slash) and `data\\x\\v` and
  `lidc-test-shifted` -> `"lidc-test-shifted"`; `test.source == "runs/ID/eval/eval_test/final_eval.json"` (posix, relative);
  `ci_unit == "item"`; `seed`, `K` from the inputs (name where they come from: pin `k` in file -> `K`).
- Determinism: run twice -> byte-identical files (sorted or fixed key order, trailing newline, no timestamp); the output ends with "\n".
- Property (Hypothesis `floats(allow_nan=False, allow_infinity=False)` in [0,1] for means; or 200 seeded randoms): `json.loads(out)[...]["mean"] == input mean` exactly,
  i.e. `float(repr(x)) == x` round trip.
- Also: items count equals input `items` (an int); `items` as float 2200.0 -> 2200 int or refusal (pin).
- Why: "copied unrounded" is a spec promise; json writers that format `%.6g` violate it quietly.

**T21. Real r5 numbers regression** (integration; skipped if `runs/r5` is absent, never reads `data/`)
- Input: `docs/results/headline.json` (tracked) after creation, source `runs/r5/eval/eval_test/final_eval.json` if present.
- Expected: `test.metrics.iou_tc.mean` rounds to 0.645, `f5` 0.816, `precision_vote` 0.695, `recall_vote` 0.748, `dice_vote` 0.701, `legal` 0.994
  (legal in test: 0.994 per [measured: runs/r5/eval/eval_test/final_eval.json]); `items == 2200`; `val.iou_tc` 0.618, `f5` 0.806, `precision` 0.689, `recall` 0.716.
  Doctor ceiling 0.732 is NOT in headline.json (it is a different file; test asserts no key `ceiling`).
- If runs/ is present: regenerating with `main(["--run","r5","--model-ref","r5","--out",tmp])` gives the same `metrics` as the tracked file (no drift).
  If runs/ is absent: only the tracked file's own numbers are checked (the test must not fail in CI, must not read `runs/`; use `skipif`).
- Why: ties the public numbers to the measured run; the tracked file is the only copy CI can see.

---------------------------------------------------------------------------------------------------------------------
## D. HF card (B8)

**T22. card unchanged (G12)** (regression; unit)
- Input: fixed `base="Qwen/Qwen3-VL-2B-Instruct"`, meta `{"license":"apache-2.0"}`, man (config seed 7, branches [...], model, lora true), sel, run `Path("runs/x")`.
- Expected: output equals a golden string written to `tests/goldens/hf_card_sft.md` BEFORE the DOC phase edits `hf_push.py` (record with current code, commit it first).
  Compare bytes. Run with `lora=False` and with `gates_failed=["g1"]` too (three goldens).
- Why: B8 says "returns exactly what it returns today"; only a golden recorded from the old code proves it.

**T23. card_update: idempotence, preservation, replacement** (unit + property)
- Old card variants (parametrize):
  a. no Evaluation/Limitations sections; b. both present with stale numbers; c. Evaluation only; d. Limitations only; e. heading order reversed;
  f. empty string; g. card with YAML front matter (`---` block) and a `## Evaluation notes` heading (NOT the same as `## Evaluation`);
  h. an "Evaluation" word inside a code fence or in text; i. `## Evaluation` as the LAST line with no trailing newline; j. CRLF line endings;
  k. Evaluation section containing a subsection `### Per-view` (must the subsection go with its parent? pin: yes, until next `## `);
  l. Evaluation heading with trailing `#`s or different case (`## EVALUATION`) -> pin.
- Expected, for each: (1) `card_update(old, h)` contains each of: 6 test means with 95 % CI text, view, items, K, seed, `ci_unit`, 4 validation means;
  exactly one `## Evaluation` and exactly one `## Limitations` heading; (2) every line of `old` outside those two sections appears unchanged and in the same order
  (compute with `difflib` on lines excluding the two sections: line-level equality); (3) `card_update(card_update(old,h),h) == card_update(old,h)`;
  (4) YAML front matter is byte-identical; (5) final newline policy: ends with exactly one "\n".
- Property (random old text made from a grammar of headings, paragraphs, fences, 100 seeds): idempotent AND `card_update(old, h1)` followed by `card_update(.., h2)`
  equals `card_update(old, h2)` (the last headline wins, no residue of h1: assert no string unique to h1, e.g. its run name, remains).
- Case g adversarial: `## Evaluation notes` stays; `## Data` stays.
- Why: a regex `## Evaluation.*?(?=##|$)` eats "## Evaluation notes" or a code fence; the card is published to the private repo and cannot be unsent.

**T24. card_sections content, no secrets, no leakage** (unit)
- Output of `card_sections(h)` contains no substring from the secret list (`hf_`, `ghp_`, `github_pat_`, `AGE-SECRET-KEY-`, `Bearer `); contains
  no absolute local path (`/[U]sers/`, `/private/`, `C:\`); contains no RIGA-derived text (headline is LIDC: add `"dataset": "riga"` in a seeded headline -> `card_sections`/`main` refuses, mirroring
  hf_push's RIGA rule [code: hf_push.py "RIGA's licence forbids redistribution"]); CI text format "[0.631, 0.659]" with 3 decimals.
- `run` and `model_ref` containing markdown/HTML (`<script>`, `](http://evil)`) are escaped or refused (pin: the card is Markdown on a private repo, but a
  `model_ref` with a newline injects headings: input `"r5\n## Hacked"` -> refused or flattened to a single line).
- Why: the card is the one HF write path for headline data; injection by newline is the realistic bug.

**T25. `--card-only` upload gate** (integration; fake hub as in tests/test_p2_hf_push.py)
- Setup: monkeypatch `H.api` to a fake hub recording calls (`repo_info`, `create_commit`, `create_branch`, `upload_file`); `sys.stdin` = StringIO("tok\n"); `headline` file via tmp_path;
  `whoami` patched to return `token()` shapes from the existing test helper.
- Expected on the happy path: exactly ONE `create_commit` (or one upload) with exactly one operation whose path is `README.md`, `revision == V`, `repo_type == "model"`;
  no `create_branch` call unless the revision does not exist (pin: card-only must not create branches; a missing revision is an error). The README body equals `card_update(<fetched old card>, headline)`.
- Gate failures (each: return code non-zero or SystemExit, `create_commit` NOT called, zero uploads): token with over-broad role (reuse the 6 parametrized tokens of
  `test_an_over_broad_token_is_refused`); repo is public; empty stdin; token on stdin with trailing spaces only; `--run` given together with `--card-only` (argparse error);
  `--headline` missing file; `--headline` is invalid JSON; headline with `dataset: riga`; repo name not in the allowlist/`token_problems` repo mismatch.
- Order property: `token_problems` and `repo_info(...).private` are evaluated BEFORE any network read of the old card (assert call order list: whoami, repo_info, then download).
- Token never appears in stdout/stderr (capsys) or in the commit message, even on failure: assert `"tok" not in captured.out + captured.err` using a distinctive token `"hf_SECRETXYZ"`.
- Why: this path writes to Hugging Face with the owner's write token; the gates must be proven to run first and to block.

---------------------------------------------------------------------------------------------------------------------
## E. Repo-level invariants, git tracking, secrets scan, build

**T26. Tracked-files invariant (G1, B4, B10)** (integration; uses `git ls-files`, read-only)
- Expected: `git ls-files docs | grep -v '^docs/research/' | wc -l` equals `find docs -type f -not -path 'docs/research/*'` (excluding `.DS_Store` which `.gitignore` must list);
  `docs/plan/` has exactly Direction.md, Architecture.md, Research2.0.md, previous-vs-new.md, all tracked; `docs/research/research2.0` absent;
  ignored patterns present via `git check-ignore`: `site/x`, `.cache/x`, `.lycheecache`, `docs/research/a.pdf` -> exit 0 each.
- Adversarial: a new file `docs/new.md` created in a temp git repo copy (`git init` in tmp_path with the same tree) is untracked -> invariant reports it. Test the comparison function on the
  temp repo, not on the real one, so the test is stable before and after the owner's commit. The real-repo variant is marked `integration` and is only run by the DOC gate, since before
  the owner's first `git add` it must FAIL (that is its purpose); make it `xfail(strict=False)` is not allowed: use a pytest marker `post_commit`, excluded from `make test`.
- Never invoke any mutating git command (CLAUDE.md): `git init` in `tmp_path` is allowed because it does not touch the repo; no `git add` in the real repo.
- Also assert `tests/test_docs.py` source contains no string `runs/` or `data/` as a path it opens (`ast` walk of `open`, `read_text`, `Path(...)` string literals): B4.
- Why: the project's failure today is 14 untracked doc paths; the test must be stable when the commit has not happened yet.

**T27. Secrets scan, B12** (unit; tmp_path tree; the scanner is a function the owner runs, not a gate)
- Seeded files in `docs/`: `a.md` with `hf_abcdefghijklmnopqrstuvwxyz0123456789ABCD` (HF token shape: `hf_` + 34 alphanumerics), `b.md` with
  `ghp_` + 36 alnum, `c.md` with `github_pat_11AAAAAAA0abcdefghijklmnopqrstuvwxyz`, `d.md` with `ssh -p 40022 root@ssh5.vast.ai`, `e.md` with `203.0.113.7`,
  `f.md` with `AGE-SECRET-KEY-1ABCDEF...`, `g.md` with `-----BEGIN PRIVATE KEY-----`.
- Expected: each of the seven files reported with line number and rule name; the report prints the rule and path but NOT the secret value (assert the secret string is not in the output).
- Clean (must print nothing): version strings `1.17.1`, `0.1.0`, `5.20.1` (look like partial IPs); `127.0.0.1`? (loopback is an IPv4: pin as reported or allowed-list `127.0.0.1`, `0.0.0.0`; document in test),
  `1.2.3.4.5`, a hex hash (64 hex chars, a SHA-256), the words "hf_push", "hf_pull", `objdel-hf-write` (keychain NAME, not a token), "ghp_" alone, "ssh -p N root@" placeholder (literal N).
- Edge: token split across a line wrap; base64-looking 40-char strings (not flagged); binary file (PNG) in `docs/figures` -> skipped without crash; file > 5 MB -> scanned in chunks;
  `docs/research/2212.07020v1.pdf` (a PDF in the untracked list) -> skipped as binary, not decoded with errors.
- Real repo smoke (read-only): scan the untracked docs files list from `git status --porcelain` and assert zero findings; it is allowed to be the owner's pre-`git add` check.
- Why: the scan is the only control before a public push of 14 untracked paths; both false negatives (leak) and false positives on version numbers (mutes the control) matter.

**T28. Build gates (B1, B2, G2)** (integration; requires zensical and mkdocstrings installed, skipped otherwise; marked `docs_build`)
- Run `zensical build --strict -c` (or `make docs`) in a COPY of the repo tree in tmp_path (so `site/` and `.cache/` do not dirty the repo).
- Expected: exit 0; stdout contains "No issues found"; `site/reference/objdel/` has a page for each of checker, data, eval, geometry, reward, secrets, targets, train;
  NO `site/reference/objdel/_gpucheck*`; `site/` contains no string matching the secret regexes of T27 (grep) and no real secret value: for the `secrets` page, set env `OBJDEL_SECRETS`
  to a file containing `hf_TESTSECRET123` and assert the string is absent from the whole site (documents that `objdel.secrets` is not executed to render).
- Needs no torch/TRL/GPU: run the build with a `PYTHONPATH`/env that hides torch (`python -c "import torch"` fails inside the venv used; use `sys.modules['torch']=None` sitecustomize) -> still exit 0.
- Negative variants in the copy: (a) add `[x](missing.md)` to `docs/a.md` -> build exits non-zero; (b) add `[x](a.md#nope)` -> anchors are `warn` and strict true -> exit non-zero;
  (c) add a docstring with an invalid cross-reference to a nonexistent module -> non-zero or warn pinned; (d) a docs page with a bad YAML front matter -> non-zero;
  (e) nav entry to a missing file -> non-zero.
- Material parity: `mkdocs build --strict` (Material 9.8.1) on the same `mkdocs.yml` exit 0, in a separate venv; skipped when not installed.
- Why: B1's exact semantics (anchors warn + strict) is the only part of the contract that can be wrong silently; variant (b) proves it.

**T29. Build performance and size gates** (performance)
- `zensical build --strict` on the full docs tree: wall time <= 60 s on the Mac and <= 120 s in CI (assertion generous; [estimate: no measurement yet]. Record the first
  measured value in `docs/results/runs.md`? NO, that file is for runs. Record it in the test as a constant `BUILD_BUDGET_S` that the author sets after measuring once.)
- `tests/test_docs.py` full run <= 5 s (it reads only text files; guards against an accidental full-repo `rglob("*")` including `.venv`, `data` or `runs`).
  Implement with `time.perf_counter()` around the whole checks and assert < 5.0.
- `doc_files(root)` on a tree with 20,000 files (create 20,000 tiny `.py` files and 100 `.md`) returns in < 2 s, and does not descend into `.venv`, `node_modules`, `.git`, `data`, `runs`, `build`,
  `site`: assert by counting `os.scandir` calls with a monkeypatch wrapper, or by placing an unreadable (mode 000) directory named `data/` and expecting no PermissionError.
- Site size gate: `du -s site` < 50 MB (guards against accidentally publishing `runs/` or `docs/research/*.pdf`; the PDF must not appear: `find site -name '*.pdf'` is empty).
- Why: a performance gate catches a scan of `.venv` (tens of thousands of files) that makes `make test` slow and could also read untrusted `data/`.

**T30. Nothing in library behaviour changed (B11, G12)** (regression)
- `git diff --stat origin/main -- src/objdel cuda tests/goldens tests/oracles` (read-only) shows only docstring-line changes for `src/objdel`: assert via `ast`: for every module
  changed, `ast.dump` of each function with docstrings stripped (`ast.get_docstring` removal via a NodeTransformer) is identical between HEAD (git show) and the working tree.
  Zero diff in `cuda/` and `tests/goldens/`.
- The full golden test modules (`test_geometry*.py`, `test_checker*.py`, `test_cuda_checker.py`, `tests/test_targets*.py` if present, reward goldens) pass unchanged: run
  `pytest -q tests/test_checker_oracle.py tests/test_geometry_oracle.py -m "not slow"` and expect pass count equal to the baseline recorded before DOC1.
- Test count: collected `pytest --collect-only -q | tail -1` is >= 1,484 + new tests (spec G7), and the set of node ids before DOC is a subset of the set after (no test removed or renamed): store baseline ids in `tests/goldens/` ? NO (adds a file);
  keep the baseline in the PR description and compute at gate time.
- Why: adding docstrings is the only allowed source edit; an AST comparison proves "no behaviour change" better than reading the diff.

**T31. CPU-vs-GPU exactness is untouched** (GPU gate; run only by the owner on a rented box, not in DOC; listed so nobody claims it was verified)
- DOC makes no GPU change; the expected result is that the CUDA-vs-CPU bit-equality tests (`tests/test_cuda_checker_gpu.py`, `tests/test_gpu_scoring.py`) have the same
  pass/skip counts as before. In DOC these tests must SKIP cleanly on the Mac (no `cuda` build): assert the skip reason is "no CUDA" and not an import error caused by the docstring edits
  (`python -I -c "import objdel; import objdel._gpucheck as g; print(g.build())"` prints `cpu`).
- Also: `src/objdel/**/*.py` still imports with docstrings (an `import` of every module with `importlib` in a loop, skipping modules that need torch: those are imported only when torch is present).
- Why: a docstring containing a stray triple quote or backslash escape (`"\d"` -> SyntaxWarning, an error in Python 3.12+ with `-W error`) breaks import. Run `python -W error -c "import compileall..."`
  on every changed file: `python -I -W error -m compileall -q src/objdel` exits 0.

**T32. Lint rule scope (B3, G4)** (integration; runs ruff)
- `ruff check src/objdel --select D100,D101,D102,D103,D104,D106` -> 0 errors after DOC (63 before, [measured: spec G4]); record and assert nothing less than 0.
- Seed a temp copy: add a module with each of a missing module docstring (D100), class (D101), method (D102), function (D103), package `__init__` (D104), nested class (D106) -> each rule reports exactly once
  (`ruff check --output-format json` and count by code). A private function `_x` (no docstring) -> not reported (D103 ignores private); `__init__` method (D107) -> not reported (not enabled);
  D105 magic methods -> not reported; D200-series (docstring style) -> not enabled (a multi-line docstring with the summary on the second line passes); a test file `tests/test_x.py` without docstrings
  -> not reported; `experiments/phase2/x.py` without docstring -> not reported; `archive/old/x.py` -> excluded (`extend-exclude`).
- Config retention: parse `pyproject.toml` with `tomllib`: `ignore == ["E731","E741"]`, `select` still contains E, F, I, B, UP plus the D codes; `extend-exclude` contains both `experiments/phase1/lab` and `archive`;
  isort `known-first-party` unchanged; per-file-ignores for `experiments/**` still has E501, B023, B905 and now `D`; `tests/**` has `D`.
- `make lint` run on the repo exits 0 and does not scan `site/` or `.cache/` (add them to extend-exclude if ruff would find python there: create `site/x.py` with a syntax error and expect `make lint` still 0).
- Why: enabling D can silently disable other rules if `select` is replaced by a new list instead of extended.

---------------------------------------------------------------------------------------------------------------------
## F. Failure modes and adversarial repo-level cases

**T33. Frozen-file protection** (integration; git diff read-only)
- For every file matching `FROZEN_GLOBS` in `git diff --name-only origin/main`: the diff is at most one added line at the top of a `.md` (lines added == 1, deleted == 0, line index 0),
  and none for `.py` or `Dockerfile*`. Implement as a function over `git diff -U0 --numstat` and `git diff -U0` hunks: hunk header must be `@@ -0,0 +1 @@`.
- Adversarial: banner added in the middle -> fails; trailing whitespace change -> fails; whole file rewritten (line-ending change CRLF->LF) -> fails (numstat > 1).
  Move of `docs/research/research2.0/*` to `docs/plan/*` is the only allowed rename out of a frozen path; the test accepts exactly those four renames (`git diff -M --name-status` shows R100 or R9x).
  R-score below 90% for `previous vs new.md` -> fails if content changed beyond link updates.
- Why: frozen records are evidence; a mass link-fix script that touches them is the likely accident.

**T34. Link update after the move (B10)** (integration)
- Search all doc scope + frozen files for the old path strings `research/research2.0`, `research2.0/`, `previous vs new.md`; expected: no occurrences in doc-scope; occurrences in frozen paths
  are allowed only if they are historical (pin: zero in `docs/plan/**`, `README.md`, `ROADMAP.md`, `docs/**` outside research). Also, a URL-encoded form `previous%20vs%20new.md` -> none.
- Every file in `docs/plan/` has at least one inbound link from the nav or an index. Anchors: links like `docs/plan/Direction.md#...` resolve to an existing heading slug (the slugger must match
  zensical's: lower-case, spaces to hyphens, punctuation dropped: test one heading with `&`, one with a backticked word, one duplicate heading -> `-1` suffix). Pin the slug function against the built HTML id attributes
  (compare with `site/` output when T28 ran).
- Why: the move is the one change in this phase that can break many inbound links at once.

**T35. Docs CI lock (B13) and extras** (unit)
- Parse `experiments/phase2/docker/locks/docs.txt`: every non-comment, non-continuation line has `name==version`; no `>=`, `~=`, no URLs, no unpinned lines; `zensical==0.1.0` and `mkdocstrings-python==2.0.10`
  present; the set of top-level packages named in `pyproject.toml` `docs` extra is a subset of the lock; the lock contains no `torch`, `trl`, `vllm`, `transformers` (B2: build needs none);
  every transitive dependency is pinned: if `uv pip compile` is available, `uv pip compile --universal pyproject docs extra` output equals the lock (skip otherwise);
  `pip install --dry-run --no-deps -r docs.txt` is not used (network).
- The `.github/workflows/tests.yml` job `docs` text contains `locks/docs.txt` and `zensical build --strict`, and the `python` job runs `tests/test_docs.py` (YAML-parse the workflow with `yaml.safe_load`; the key `on` parses as boolean True in YAML 1.1: use `data[True]` or quote); `links.yml` actions are pinned
  `uses: owner/name@<40 hex>` (regex `@[0-9a-f]{40}\b`), version comment `v2.9.0` present for lychee; the `actions/cache` pin is a SHA; weekly cron present; the weekly job has `continue-on-error: true`;
  PR trigger has `paths: ['**/*.md']`.
- Why: unpinned actions and loose lock files are the supply-chain failures the spec's pins aim at.

**T36. Hostile Markdown inputs do not crash or hang the checks** (property/fuzz; time-bounded)
- Generate 300 random Markdown texts (seeded `random.Random(0)`): random mixes of `[`, `]`, `(`, `)`, backticks, `<`, `>`, `%`, `#`, newlines, long runs of `[` (10,000 chars), 1 MB of backticks,
  deeply nested `[[[[...](x)`.
- Expected for `relative_links`, `backtick_paths`, `make_mentions`, `version_mentions`, `cost_mentions`: each returns within 0.5 s per input (catastrophic regex backtracking is the bug: input
  `"[" * 20000 + "](" ` and `"`" * 100000 + "src/"` and `"trl " + "1." * 50000` must complete), and returns a list/set of str (no exception), and is deterministic (same output on second call).
  Idempotence property: `relative_links(text + "\n" + text) == relative_links(text) * 2` order-preserving for the list-returning ones (or set equality for the set ones).
- Use `signal.alarm(5)` or `pytest-timeout` if already in dev (it is not: use `signal.alarm`, Unix only, skip on Windows).
- Why: regex ReDoS in a test helper would make CI hang and be mistaken for infrastructure flakiness.

**T37. Order independence and no hidden state** (property)
- Run all 13 checks twice, once with `doc_files` returning reversed order, and once with `PYTHONHASHSEED` 0 and 1 (subprocess): reports are equal as sorted lists; reports are sorted before
  comparison and the failure message lists items in sorted order (so CI logs diff cleanly).
- `tests/test_docs.py` does not write any file in the repo: record `git status --porcelain` before and after (read-only) -> equal; and `tmp_path` is the only place written.
- Why: set iteration order leaking into messages causes flaky diffs; writing into the repo from a test breaks the "tracked only" promise.

**T38. Platform and encoding differences (Mac vs Linux CI)** (integration)
- Case-sensitivity (see T2); `Path.exists` with a trailing `/` on a file (`docs/a.md/` should be a miss); unicode in file names (NFC vs NFD: `é` composed vs decomposed in a link
  vs on APFS, which normalises, vs Linux, which does not): link written in NFC, file in NFD -> reported only when they differ byte-wise on a case/normalisation-sensitive FS; pin as "the check uses
  `unicodedata.normalize("NFC")` on both sides and ALSO requires the exact bytes to exist in `os.listdir`", so the Mac matches CI.
- `.DS_Store` in `docs/` is ignored by doc_files and by the nav check; Windows line endings in `Makefile` (`\r\n`) -> `make_targets` still works.
- Why: the check runs on the Mac but CI is Linux; every Mac-only pass is a CI failure later.

**T39. Hypothesis-style invariants for the pure functions as a set** (property; 100 seeds each)
- `relative_links(text)` never contains a string starting with `http://`, `https://`, `mailto:`, `#`.
- `backtick_paths(text)` results all start with one of the six prefixes and contain none of ` *<{`.
- `make_targets(Makefile)` is a subset of tokens appearing before a `:` at column 0.
- `adr_files(d)` results all match `^\d{4}-.+\.md$`.
- `version_mentions(t, n)` results all match `^\d+(\.\d+)+`.
- `cost_mentions(t)` results have kind in {"alarm","hard stop"} and dollars >= 0 and finite.
- `headline_strings(h)` length in {11, 12}; all strings non-empty; round trip: each numeric token appears in `format(x, ".3f")` of the matching mean.
- Why: the invariants hold regardless of the implementation's regexes and catch the case where a fix for one bug breaks another.

**T40. Failure mode: the checks cannot be disabled by accident** (meta)
- `tests/test_docs.py` has no `pytest.skip`, `skipif`, `xfail` or `importorskip` on the 13 repo-level checks (AST scan): the only allowed skips are for the optional build tests (T28) and runs/-dependent T21.
- Each repo-level check asserts it examined > 0 items (e.g. check 1 saw >= 1 link, check 2 saw >= 10 backtick paths, check 3 saw >= 1 make mention, check 6 saw all 16 configs, check 7 saw
  >= 1 frozen-skipped mention): a count-of-examined assertion makes a vacuous pass (`doc_files` returned `[]`) a failure.
- Collected test count for `tests/test_docs.py` and `tests/test_headline.py` is >= 1 per check + 8 negative tests (B6) + the legal-input tests (4); `pytest --collect-only -q tests/test_docs.py | tail -1`.
- Why: a green suite over an empty file list is the commonest docs-test bug.

---------------------------------------------------------------------------------------------------------------------
## Pytest sketches for the 5 most important tests

Priority order: T14 (negative tests prove checks are not vacuous), T7+T8 (version check, money-relevant drift), T19 (headline refusal matrix), T23 (card_update idempotence and preservation), T25 (HF upload gate).

```python
# --- tests/test_docs.py (sketch of T14 + T7/T8) -------------------------------------------------
import pytest
from tests import docs_check as D

BAD_AND_LEGAL = {
    "README.md": "see [x](docs/gone.md) and `make nope`\n",
    "docs/a.md": "TRL 1.16 and trl==1.17.1 and vLLM 0.09 s\n",
    "docs/research/frozen.md": "trl==1.16.0\n",
    "archive/old/README.md": "was `src/objdel/retired.py`\n",
    "Makefile": ".PHONY: test\ntest:\n\ttrue\n",
    "pyproject.toml": '[project]\nversion="0.1.0"\ndependencies=["trl==1.17.1"]\n',
}

def mkrepo(root, files):
    for rel, text in files.items():
        p = root / rel; p.parent.mkdir(parents=True, exist_ok=True); p.write_text(text)
    return root

def test_check_reports_exactly_the_seeded_faults_and_not_the_legal_inputs(tmp_path):
    root = mkrepo(tmp_path, BAD_AND_LEGAL)
    assert [m for m in D.broken_links(root / "README.md", root)] == ["docs/gone.md"]
    assert D.make_mentions((root / "README.md").read_text()) - D.make_targets(root / "Makefile") == {"nope"}
    allowed = D.allowed_versions(root, "trl")
    bad = [v for v in D.version_mentions((root / "docs/a.md").read_text(), "trl") if v not in allowed]
    assert bad == ["1.16"]                       # "trl==1.17.1" and "vLLM 0.09 s" are legal, only 1.16 is not
    assert D.is_frozen("docs/research/frozen.md")  # frozen path skipped by the caller
    assert D.is_frozen("archive/old/README.md") is False  # archive README is exempt from check 2 by the test, not by is_frozen

@pytest.mark.parametrize("check", ["broken_links", "make_mentions", "version_mentions"])
def test_a_check_that_always_passes_is_caught(monkeypatch, tmp_path, check):
    monkeypatch.setattr(D, check, lambda *a, **k: [] if check != "make_mentions" else set())
    root = mkrepo(tmp_path, BAD_AND_LEGAL)
    with pytest.raises(AssertionError):
        test_check_reports_exactly_the_seeded_faults_and_not_the_legal_inputs(root)

@pytest.mark.parametrize("text, allowed, bad", [
    ("trl 1.17", {"1.17.1"}, []), ("trl 1.1", {"1.17.1"}, ["1.1"]), ("trl 1.17.10", {"1.17.1"}, ["1.17.10"]),
    ("torch 100 %", {"2.14.1"}, []), ("vllm 0.5 ms", {"0.21.0"}, []), ("trlx 1.0", {"1.17.1"}, []),
])
def test_version_two_part_rule_and_units(text, allowed, bad):
    name = text.split()[0] if text.split()[0] != "trlx" else "trl"
    got = [v for v in D.version_mentions(text, name) if not any(a == v or a.startswith(v + ".") for a in allowed)]
    assert got == bad
```

```python
# --- tests/test_headline.py (sketch of T19 + T20) ----------------------------------------------
import json
import pytest
from experiments.phase2 import headline as Hd

METRICS = {k: {"mean": 0.5, "ci95": [0.4, 0.6]} for k in ("iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote", "legal")}

def write_run(root, run="r5", test=None, val=None):
    for split, over in (("test", test), ("val", val)):
        d = root / "runs" / run / "eval" / f"eval_{split}"; d.mkdir(parents=True)
        doc = {"split": split, "view": f"data/eval_views/lidc-{split}-shifted", "items": 2200, "k": 8, "seed": 11,
               "ci_unit": "item", "metrics": METRICS} | (over or {})
        if over == "missing": continue
        (d / "final_eval.json").write_text(json.dumps(doc))

@pytest.mark.parametrize("test, val", [
    ({"split": "val"}, None), (None, {"split": "test"}), ({"view": None}, None), ({"view": ""}, None),
])
def test_refusal_writes_nothing_and_keeps_the_old_file(tmp_path, monkeypatch, test, val):
    monkeypatch.chdir(tmp_path)
    write_run(tmp_path, test=test, val=val)
    out = tmp_path / "h.json"; out.write_text("OLD")
    assert Hd.main(["--run", "r5", "--model-ref", "x", "--out", str(out)]) != 0
    assert out.read_text() == "OLD"
    assert not list(tmp_path.glob("h.json.*"))

def test_missing_input_refuses_and_leaves_no_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "runs/r5/eval/eval_test").mkdir(parents=True)
    out = tmp_path / "h.json"
    assert Hd.main(["--run", "r5", "--model-ref", "x", "--out", str(out)]) != 0
    assert not out.exists()

def test_means_are_copied_unrounded_and_output_is_deterministic(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    m = dict(METRICS); m["iou_tc"] = {"mean": 0.6454999999999999, "ci95": [0.6312345678901234, 0.6591234567890123]}
    write_run(tmp_path, test={"metrics": m})
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    for o in (a, b):
        assert Hd.main(["--run", "r5", "--model-ref", "x", "--out", str(o)]) == 0
    assert a.read_bytes() == b.read_bytes() and a.read_text().endswith("\n")
    got = json.loads(a.read_text())["test"]["metrics"]["iou_tc"]
    assert got == {"mean": 0.6454999999999999, "ci95": [0.6312345678901234, 0.6591234567890123]}
    assert set(json.loads(a.read_text())["test"]["metrics"]) == set(METRICS)
    assert json.loads(a.read_text())["test"]["view"] == "lidc-test-shifted"

@pytest.mark.parametrize("bad", ["", "../x", "a/b", "a b"])
def test_run_id_that_could_escape_runs_is_refused(tmp_path, monkeypatch, bad):
    monkeypatch.chdir(tmp_path)
    assert Hd.main(["--run", bad, "--model-ref", "x", "--out", str(tmp_path / "h.json")]) != 0
```

```python
# --- tests/test_p2_hf_push.py additions (sketch of T23 + T25) ------------------------------------
import io
import pytest
from experiments.phase2 import hf_push as H

HEAD = {"run": "r5", "model_ref": "r5", "ckpt": "c", "seed": 11, "K": 8, "ci_unit": "item",
        "test": {"view": "lidc-test-shifted", "items": 2200, "metrics": {k: {"mean": 0.645, "ci95": [0.63, 0.66]}
                 for k in ("iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote", "legal")}},
        "val": {"view": "lidc-val-shifted", "items": 128, "metrics": {k: {"mean": 0.618, "ci95": [0.6, 0.63]}
                for k in ("iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote", "legal")}}}

OLDS = [
    "", "---\nlicense: apache-2.0\n---\n\n# t\n\nbody\n",
    "# t\n\n## Evaluation\n\nstale 0.1\n\n## Data\n\nkeep me\n",
    "# t\n\n## Evaluation notes\n\nkeep me too\n\n## Limitations\n\nold\n",
    "# t\n\n```\n## Evaluation\n```\n\ntext\n", "# t\n\n## Evaluation",
]

@pytest.mark.parametrize("old", OLDS)
def test_card_update_is_idempotent_and_keeps_every_other_line(old):
    once = H.card_update(old, HEAD)
    assert H.card_update(once, HEAD) == once
    assert once.count("\n## Evaluation\n") + once.startswith("## Evaluation\n") == 1
    assert once.count("## Limitations") == 1
    kept = [l for l in old.splitlines() if l.strip() and l != "stale 0.1" and l != "old"]
    for line in kept:
        if line == "## Evaluation" and "```" not in old: continue   # replaced, not kept
        assert line in once.splitlines()
    assert "0.645" in once and "0.618" in once and "[0.63, 0.66]" in once

def test_a_second_headline_leaves_no_residue_of_the_first():
    h2 = {**HEAD, "run": "r6"}
    assert "r6" in H.card_update(H.card_update(OLDS[2], HEAD), h2)
    assert "r5" not in H.card_update(H.card_update(OLDS[2], HEAD), h2).replace("r5-extra", "")

def test_card_only_uploads_one_file_after_the_gates(monkeypatch, tmp_path, capsys):
    calls = []
    class Hub:
        def repo_info(self, repo, repo_type): calls.append("repo_info"); return type("I", (), {"private": True})()
        def create_commit(self, repo, operations, commit_message, revision, repo_type):
            calls.append(("commit", [o.path_in_repo for o in operations], revision))
        def create_branch(self, *a, **k): calls.append("branch")
    monkeypatch.setattr(H, "api", lambda token: Hub())
    monkeypatch.setattr(H, "whoami", lambda token: __import__("tests.test_p2_hf_push", fromlist=["token"]).token())
    monkeypatch.setattr(H, "fetch_old_card", lambda hub, repo, rev: "# t\n")      # name pinned by the implementer
    p = tmp_path / "h.json"; p.write_text(__import__("json").dumps(HEAD))
    monkeypatch.setattr("sys.stdin", io.StringIO("hf_SECRETXYZ\n"))
    assert H.main(["--card-only", "--repo", "svenk029/objdel-sft", "--revision", "rl-full-3ep", "--headline", str(p)]) == 0
    assert calls == ["repo_info", ("commit", ["README.md"], "rl-full-3ep")]       # exactly one file, no branch
    out = capsys.readouterr(); assert "hf_SECRETXYZ" not in out.out + out.err

@pytest.mark.parametrize("private", [False])
def test_a_public_repo_blocks_the_upload(monkeypatch, private):
    calls = []
    Hub = type("Hub", (), {"repo_info": lambda s, *a, **k: type("I", (), {"private": private})(),
                           "create_commit": lambda s, *a, **k: calls.append("commit")})
    monkeypatch.setattr(H, "api", lambda t: Hub())
    monkeypatch.setattr("sys.stdin", io.StringIO("tok\n"))
    with pytest.raises(SystemExit):
        H.main(["--card-only", "--repo", "svenk029/objdel-sft", "--revision", "v", "--headline", "/nonexistent"])
    assert calls == []
```

## Coverage map (spec to tests)
B1 T28, T6; B2 T28, T35; B3 T32; B4 T26, T29, T40; B5.1 T1, T2, T3; B5.2 T4, T3; B5.3 T5, T6; B5.4/5 T12; B5.6/13 T13; B5.7 T7, T8, T9;
B5.8 T10, T11; B5.9 T17, T18; B5.10 T16; B5.11 T15; B5.12 T16; B6 T14; B7 T19, T20, T21, T6; B8 T22, T23, T24, T25; B9 T9, T11;
B10 T26, T34; B11 T30; B12 T27; B13 T35. Gates: G1 T26, G2 T28, G3 T14, G4 T32, G5 T9, G6 T17-T18-T21, G7 T30, G8 not testable offline (lychee run by CI),
G9 T16, G10 T16, G11 T12, G12 T22-T30-T31, G13 T27.

## Open questions the implementer must answer (each test above says "pin")
1. Is a root-absolute link or a directory link a break? 2. Are `.yml` configs and nested configs collected? 3. Does one bare config file name cover two phases?
4. Is a malformed ADR file name itself a report? 5. Is "seconds" a unit exemption? 6. Does `card_update` treat `## Evaluation notes` as a different section (tests say yes)?
7. `--card-only` with a missing revision: error or create (tests say error; the spec's "uploads exactly one file" implies no branch creation)?
8. Is `hypothesis` allowed in `test_docs.py` (spec: stdlib + pyyaml; hypothesis is already pinned in dev, seeded-random fallbacks are given)?
