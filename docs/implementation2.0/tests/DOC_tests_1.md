# DOC phase test cases, set 1 (written from DOC.spec.md only)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA product in scope (docs/lint/CPU only); no catalog skill needed, nothing installed.
- engineering:testing-strategy: pyramid split (many unit on pure functions in `tests/docs_check.py`, few integration on the real repo, one end-to-end on `headline`/`card_update`); gate list.
- mattpocock-skills:tdd: tests only at public seams (`docs_check` functions, `headline.main`, `hf_push.card_update/main`, `make` targets); expected values are literals from the spec, never recomputed with the code under test.
- engineering:code-review, engineering:debug: used to hunt regex edge cases (units, two-part versions, frozen paths) and failure modes (partial writes, non-idempotent update).
- engineering:documentation: doc-scope and link rules read as a documentation contract.
- engineering:system-design, engineering:architecture, huggingface-skills:huggingface-best: loaded per the brief; no design input (by design: zero design context); HF only matters for the card/upload seam.

## Conventions
Style follows `tests/test_p2_hf_push.py` (module docstring, `parametrize`, `tmp_path`, `# fmt: skip`). Markers: unit = pure function on a string or a `tmp_path` tree; property = `hypothesis` (already a dev dep); integration = real repo tree or real `make`; GPU = none in this phase (stated, section Z). `R` = repo root. `tree(tmp_path, {rel: text})` is a helper that writes files.
IDs: T01.. . "Why" = the real bug the test can catch.

---
## A. Frozen paths and doc scope (`is_frozen`, `doc_files`, `outside_docs_dir`)

**T01 is_frozen: every glob in the definition, positive** (unit, parametrize)
Input -> expected True: `docs/research/x.md`, `docs/research/reports/a/b.md`, `experiments/phase0/a.py`, `experiments/phase1/lab/z.md`, `experiments/phase2/lab/x.md`, `experiments/phase2/impl/y.md`, `experiments/cuda_checker/README.md`, `experiments/phase2/PLAN.md`, `experiments/phase2/PLAN_v2.md`, `experiments/phase2/JUDGE_1.md`, `experiments/phase2/V0_V1.md`, `experiments/phase2/RESULTS.md`, `experiments/phase2/RESULTS_r5.md`, `experiments/phase2/EXECUTION.md`, `experiments/phase2/RESEARCH.md`.
Why: a missed glob makes frozen history fail the link/version checks (red CI on files nobody may edit).

**T02 is_frozen: near misses, negative** (unit)
Expected False: `docs/researchX/a.md`, `docs/plan/Direction.md`, `docs/reference/configs.md`, `experiments/phase2/README.md`, `experiments/README.md`, `experiments/phase2/docker/locks/train.txt`, `experiments/phase2/headline.py`, `experiments/phase2/sub/PLAN.md` (`experiments/*/PLAN*.md` is one level only), `src/objdel/x.py`, `README.md`, `archive/a/README.md`, `docs/research` (no trailing slash and no file).
Why: over-broad globs (`**` instead of `*`, prefix match without `/`) silently exempt live docs from every check.

**T03 is_frozen: Windows/odd separators and leading `./`** (unit)
`is_frozen("./docs/research/a.md")` and `is_frozen("docs\\research\\a.md")`: pick one documented behaviour (normalise or raise); assert it is the same for both, never True for one and False for the other.
Why: `Path.relative_to` vs `str` forms differ between Mac and CI.

**T04 doc_files returns exactly doc scope** (unit, tmp tree)
Tree: `README.md`, `ROADMAP.md`, `docs/a.md`, `docs/sub/b.md`, `docs/plan/Direction.md`, `docs/research/r.md`, `docs/x.txt`, `docs/pic.png`, `experiments/README.md`, `experiments/phase2/README.md`, `experiments/phase2/PLAN.md`, `archive/old/README.md`, `archive/old/NOTES.md`, `src/objdel/README.md`, `tests/README.md`, `.github/pull_request_template.md`, `node_modules/x/README.md`.
Expected set (relative, sorted): `README.md, ROADMAP.md, docs/a.md, docs/sub/b.md, docs/plan/Direction.md, experiments/README.md, archive/old/README.md`. `experiments/phase2/README.md` is NOT in scope (spec lists only `experiments/README.md`). Pin this literally; if the owner wants it in, the spec must change.
Why: scope drift is the most likely silent weakening of all 13 checks.

**T05 doc_files on missing files and empty tree** (unit)
Empty `tmp_path` -> `[]` (no exception). Tree with only `docs/` (no README.md) -> only docs files. `docs` as a file not a dir -> `[]` or clear error, not a traceback from `rglob`.
Why: negative tests in B6 run on tiny trees.

**T06 outside_docs_dir** (unit)
Input files `[README.md, docs/a.md, experiments/README.md, archive/x/README.md, docs/plan/p.md]` -> `[README.md, experiments/README.md, archive/x/README.md]`, order preserved. Input with absolute paths under root gives the same answer. A file named `docsfoo/a.md` is outside.
Why: B5.1 applies link checking only to these; `startswith("docs")` bug.

**T07 doc_files ignores symlink loops and `.venv`/`site`/`.cache`** (unit)
Tree has `docs/loop -> docs` symlink and `site/index.md`, `.venv/lib/README.md`. Expected: terminates, no entries from `site/`, `.venv/`, `.cache/`.
Why: a locally built `site/` and a venv would otherwise be scanned, failing on third-party text.

---
## B. Link checks (`relative_links`, `broken_links`, B5.1)

**T08 relative_links extraction, parametrize** (unit). Expected list per input:
- `[a](b.md)` -> `["b.md"]`
- `[a](<my file.md>)` -> `["my file.md"]` (angle form kept as text inside)
- `[a](my%20file.md)` -> `["my%20file.md"]` (decoding happens in `broken_links`; pin whichever, but `broken_links` must find the file `my file.md`)
- `[a](https://x.org)`, `[a](http://x)`, `[a](mailto:a@b.c)`, `[a](#sec)` -> `[]`
- `[a](b.md#sec)` -> `["b.md#sec"]` (anchor dropping is later)
- `![img](../fig/a.png)` -> `["../fig/a.png"]` (images are links)
- `[a](b.md "title")` -> `["b.md"]`
- `[a](b.md) and [c](d.md)` on one line -> 2 items
- link inside a fenced code block -> `[]`; inside inline backticks -> `[]`
- `[a][ref]` with `[ref]: x.md` definition -> `["x.md"]`
- `[a](dir/(paren).md)` -> `["dir/(paren).md"]` (balanced parens)
- `[a]()` -> `[]` and no crash.
Why: each is a known false-positive or false-negative source (fenced examples in contributing guide).

**T09 broken_links: existing, missing, anchor, URL-encoded** (unit, tmp_path)
Tree: `README.md` containing `[a](docs/x.md#top) [b](docs/missing.md) [c](docs/my%20file.md) [d](<docs/my file.md>) [e](#local) [f](docs/)` ; files `docs/x.md`, `docs/my file.md`.
Expected: `["docs/missing.md"]` only (`docs/` directory target exists -> not broken; decide and pin). The anchor `#top` is NOT validated (anchors are the site build's job, B1).
Why: B5.1 spec says exactly: decode %20, accept `<...>`, drop anchor.

**T10 broken_links: relative to the file's own folder, not root** (unit)
`experiments/README.md` has `[a](phase2/README.md)` and `[b](../docs/a.md)`; both files exist -> `[]`. `[c](docs/a.md)` from `experiments/README.md` (wrong base) -> reported.
Why: resolving against root hides real breakage.

**T11 broken_links: escapes the repo root** (unit, adversarial)
`README.md` has `[a](../../etc/passwd)` and a file exists at `tmp_path/../etc/passwd` -> still reported broken (target outside root counts as broken) OR documented; must never read outside root.
Why: link checker must not leak a pass for paths outside tracked tree.

**T12 broken_links: case sensitivity** (unit)
Link `Docs/A.md` when `docs/a.md` exists -> reported (Linux CI is case sensitive; Mac FS is not). Implement via listing names, not `exists()`.
Why: passes on Mac, fails in CI and on the site.

**T13 broken_links on the real repo** (integration)
For every file in `outside_docs_dir(doc_files(R))`, `broken_links(f, R) == []`. Also include `archive/**/README.md` here (links in archive README must resolve, only backtick paths are exempt).
Why: this is check 1 itself (G3).

**T14 property: relative_links never returns an absolute URL** (property, hypothesis)
Random markdown from alphabet `[]()<>#!:/ .%a-z` ; invariant: no result starts with `http:`, `https:`, `mailto:`, `#`; function never raises.
Why: regex catastrophic input; crash in CI.

**T15 relative_links performance** (performance gate)
1 MB text of `[` repeated 10^6 times, and 10^5 `](`: finishes < 1 s (no catastrophic backtracking).
Why: nested-quantifier regex DoS makes `make test` hang.

---
## C. Backtick paths (`backtick_paths`, B5.2)

**T16 backtick_paths extraction, parametrize** (unit)
- `` `src/objdel/geometry/mask.py` `` -> returned
- `` `experiments/phase2/configs/r1_lfm.yaml` `` -> returned
- `` `tests/test_docs.py` ``, `` `cuda/gpucheck.cpp` ``, `` `docs/plan/` ``, `` `archive/old/` `` -> returned
- `` `src/objdel/*.py` `` (star), `` `src/<pkg>/x.py` `` (`<`), `` `src/{a,b}.py` `` (`{`), `` `src/objdel x` `` (space) -> NOT returned
- `` `make test` ``, `` `objdel/geometry` ``, `` `data/lidc` ``, `` `runs/x/final_eval.json` `` -> NOT returned (not a listed prefix)
- `` `source/x` ``, `` `srcs/x` `` -> NOT returned (prefix must be followed by `/`)
- `` `src` `` alone -> NOT returned
- double-backtick span ``` ``src/a b`` ``` -> not returned; fenced block containing `` `src/x` `` -> decide; pin (recommended: fenced content ignored).
- `` `src/objdel/x.py:12` `` and `` `src/objdel/x.py::fn` `` -> pin: either stripped to `src/objdel/x.py` or reported; must be consistent and documented. (Real docs use `file:line` heavily; this is the highest false-positive risk.)
- trailing punctuation inside span `` `docs/a.md.` `` -> reported as nonexistent (it is not stripped).
Why: B5.2's four exclusions are exact; each wrong rule produces a false alarm on real docs or hides a rotted path.

**T17 check 2 negative: seeded rotted path** (unit, tmp_path)
`README.md` contains `` `src/objdel/gone.py` `` and `` `docs/gone/` ``; neither exists. Check reports both, with file name. Creating them makes the report empty.
Why: B6 requires a seeded fault test for check 2.

**T18 check 2 legal: archive README may name retired paths** (unit)
`archive/old/README.md` with `` `experiments/phase0/gone.py` `` -> not reported. Same text in `docs/guide.md` -> reported. Same in `archive/old/NOTES.md` (not in doc scope) -> not read.
Why: B6 explicit legal input.

**T19 check 2: directory vs file, trailing slash** (unit)
`` `docs/guides/` `` exists as dir -> OK; `` `docs/guides` `` -> OK; `` `docs/guides/x.md` `` missing -> reported. A path that exists only as git-ignored (e.g. `docs/research/foo.pdf` untracked): the check reads tracked files only (B4). Test with a tree that has an untracked file: use `git ls-files` list injection or assert the function takes an existing-paths provider, so that "exists on disk but not tracked" is reported. (If implementation uses disk `exists()`, this test documents the gap; owner to decide.)
Why: B4 states reads only tracked files; ignored artefacts (`site/`) would make a path pass locally and fail on a clean clone.

**T20 check 2 on the real repo** (integration)
All backtick paths in doc scope (except `archive/**/README.md`) exist. Today's expected failure count is the baseline recorded in DOC1; after DOC the count is 0.
Why: G3.

**T21 property: backtick_paths output obeys the filter** (property)
Random spans; every returned item starts with one of the 6 prefixes and contains none of space, `*`, `<`, `{`.
Why: invariant independent of the regex.

---
## D. Makefile targets and mentions (B5.3)

**T22 make_targets parsing** (unit, tmp_path Makefile)
Makefile text:
```
PY := x
.PHONY: a b
a:  ## doc
	echo
$(STAMP): dep
b c: dep
VAR ?= 1
d::
%.o: %.c
```
Expected `{a, b, c, d}` plus `.PHONY` excluded. `$(STAMP)`, `%.o`, `PY`, `VAR` excluded. A recipe line (tab-indented) containing `x: y` is not a target.
Why: variable assignments with `:=` look like targets.

**T23 make_targets on the real Makefile** (integration)
Contains `docs`, `docs-serve`, `headline`, `test`, `lint`, `build`, `preflight`, `cuda`, `wheels`, `cov`, `fmt`, `setup`, `test-all`. All three new targets are listed in `.PHONY` (parse `.PHONY` lines, including continuation `\`).
Why: B-interfaces says all three in .PHONY; a file named `docs` (the directory!) makes a non-phony `docs` target silently do nothing ("make: 'docs' is up to date").

**T24 make_mentions, parametrize** (unit)
- `` `make test` `` -> `{test}`; `make docs-serve` -> `{docs-serve}`; `make test-all.` -> `{test-all}` (trailing period stripped)
- `make headline RUN=r5 MODEL_REF=x` -> `{headline}` (variables not targets)
- `make -C dir docs`, `make -j4 test` -> pin: flags skipped, target found
- "to make it work", "make sure", "make a plan", "we make progress" -> `{}` (English; this is the largest false-positive risk. Spec says "a `make X` mention": pin that only backtick/code-fenced/`$ make` forms count, or an allowlist of English words; test must assert "make sure" is never reported).
- `makefile`, `cmake`, `Makefile` -> `{}`; `make` alone -> `{}`
- `make test && make lint` -> `{test, lint}`
Why: check 3 would fail on prose, or miss `make docs-serve` after rename.

**T25 check 3 negative** (unit): doc says `` `make deploy` ``, Makefile has no `deploy` -> reported; add the target -> clean.
**T26 check 3 real repo** (integration): every mention in doc scope is a target (G3).

---
## E. ADR checks (`adr_files`, `adr_index`, `adr_headings`, B5.4, B5.5)

**T27 adr_files** (unit)
Dir: `0001-a.md`, `0012-docs-site-with-zensical.md`, `_template.md`, `README.md`, `notes.md`, `001-short.md`, `0003-x.txt`, `0004-.md`(?), subdir `sub/0005-z.md`.
Expected: `0001-a.md`, `0012-...md` only. 3-digit and non-md excluded; subdirectories not recursed. Sorted.
Why: `_template.md` would wrongly demand an index row.

**T28 adr_index** (unit)
README table:
```
| # | Title | Status |
|---|---|---|
| 1 | [Use X](0001-a.md) | Accepted |
| 2 | [Y](./0002-b.md) | Accepted |
| 3 | [Z](0003-c.md#ctx) | Proposed |
```
plus prose line `see [x](0009-prose.md)` outside the table. Expected `{0001-a.md, 0002-b.md, 0003-c.md}`; the prose link is NOT an index row (pin). Row linking `<0004 d.md>` angle form parsed. Rows without link ignored. Empty README -> `set()`.
Why: `./` prefix and anchors break naive string compare.

**T29 check 4 negative both directions** (unit)
(a) file `0007-new.md` exists, no row -> reported "missing from index". (b) row for `0008-ghost.md`, no file -> reported "row without file". (c) both clean -> empty. (d) duplicate row for same file -> pin (report or accept); (e) file exists in index with different case -> reported.
Why: B6 negative for check 4 + both directions in G11.

**T30 adr_headings** (unit)
File with `## Context`, `## Decision`, ... returned as a set of heading texts. Edge: heading inside fenced code block not counted; `### Context` (level 3) not counted as `## Context`; `##Context` (no space) not counted; trailing spaces `## Decision ` counted; `## Options considered` vs `## Options Considered` (case) -> pin exact-case match (spec text is lowercase c); `## Trade-off analysis` with a non-ASCII hyphen not accepted.
Why: check 5 headings are exact strings.

**T31 check 5 each requirement missing, parametrize over 9 items** (unit)
Base valid ADR containing `**Status:** Accepted`, `**Date:** 2026-11-08`, `**Deciders:** owner`, and 6 headings. Remove exactly one at a time -> reported naming that item (9 cases). Line `Status: Accepted` (no bold) -> reported. `**Status:**` not at line start (indented ` **Status:**` or in a blockquote `> **Status:**`) -> reported (spec: "lines starting"). Valid base -> clean.
Why: B6 does not list check 5 negative, so this fills a gap; each of 9 items needs its own detector.

**T32 check 5 on real ADRs** (integration): every `docs/architecture/decisions/NNNN-*.md` passes, including new `0012-docs-site-with-zensical.md`; `_template.md` itself contains all 9 items (so copies pass).
**T33 ADR immutability** (integration, uses `git diff main -- docs/architecture/decisions/00*.md`, skipped if no git): for ADRs existing on main, the only changed lines match `^[-+]\*\*Status:\*\*` or `^[-+]- \[[ x]\]`. (Spec constraint: accepted ADR changes only Status and check boxes.) Mark as owner-run/CI, not in `make test` if git is unavailable.
Why: silent edits to decision records.

---
## F. Config coverage (`config_files`, B5.6)

**T34 config_files** (unit)
Tree: `experiments/phase2/configs/a.yaml`, `experiments/phase2/configs/b.yml`, `experiments/phase2/configs/sub/c.yaml`, `experiments/phase3/configs/d.yaml`, `experiments/phase2/other/e.yaml`, `experiments/configs/f.yaml` (no phase dir), `experiments/phase2/lab/configs/g.yaml`.
Expected under the glob `experiments/**/configs/*.yaml`: `a.yaml`, `d.yaml`, `f.yaml` (if `**` matches zero dirs, pin), `g.yaml` (frozen paths are NOT exempt from this check in B5.6: pin this), not `b.yml`, `c.yaml`, `e.yaml`. Sorted, relative.
Why: `**` zero-dir semantics differ between `glob` and `fnmatch`.

**T35 check 6 negative** (unit): `docs/reference/configs.md` names `a.yaml` only; `d.yaml` exists -> reported by name. Mention only as a substring of another name (`rr_a.yaml` vs `a.yaml`) -> reported (whole-name or full-path match, pin which, and test substring trap).
**T36 check 6 real repo** (integration): every `experiments/**/configs/*.yaml` is named in `docs/reference/configs.md`. New config added by a later phase (e.g. `rl_shifted.yaml`) is covered the same way.

---
## G. Version mentions (`allowed_versions`, `version_mentions`, B5.7)

**T37 allowed_versions** (unit, tmp_path)
`pyproject.toml` has `"trl==1.17.1"`, `"transformers==5.9.0"`, `"torch>=2.0"` (range, not exact, ignored), `"vllm==0.20.1; sys_platform=='linux'"`; `locks/train.txt` has `trl==1.17.1`, `torch==2.13.0+cu130`, `vllm==0.20.1`; `locks/vllm.txt` has `torch==2.11.0`; `locks/vllm_native.override` is NOT read (only `*.txt`); `locks/docs.txt` exists too (read; all 'name==X' pins).
Expected `allowed_versions(root,"trl") == {"1.17.1"}`, torch `{"2.13.0+cu130","2.11.0"}` (pin whether the local `+cu130` suffix is kept; a mention `torch 2.13.0` must be legal when the lock says `2.13.0+cu130`: B6 requires `torch 2.13.0` in `Dockerfile.vllm_native` to NOT be reported, so either normalise the local label away or have the vllm_native lock pin `2.13.0`). Name matching is case-insensitive and treats `-`/`_` alike (`Torch==`), but `torchvision==0.99` is not torch.
Why: prefix name collisions (`torch` vs `torchvision`, `trl` vs `trlx`).

**T38 version_mentions, positives** (unit, parametrize; name -> expected list)
- `"TRL 1.16"` (trl) -> `["1.16"]`; `"TRL 1.16.0"` -> `["1.16.0"]`; `"trl==1.16.0"` -> `["1.16.0"]`
- `"https://github.com/huggingface/trl/blob/v1.16.0/trl/x.py"` -> `["1.16.0"]` (v stripped)
- same URL for `transformers` with `huggingface/transformers/blob/v5.9.0/` -> `["5.9.0"]`; the trl URL asked with name `torch` -> `[]`
- `"torch==2.13.0+cu130"` -> `["2.13.0+cu130"]`; `"torch 2.13.0"` -> `["2.13.0"]`
- `"vllm==0.20.1"`, `"vLLM 0.20.1"` (case-insensitive name), `"vllm\t0.20.1"` (whitespace)
- two mentions in one line -> both, in order.
- `"transformers>=5.0"`: not `==` or whitespace-then-version -> `[]` (pin; `>=` is a range not a mention).
- trl in a pip extra `"trl==1.17.1"` -> `["1.17.1"]`.
**T39 version_mentions, units and non-mentions** (unit)
`"vLLM 0.09 s"`, `"vllm 0.09s"`, `"torch 2 ms"`? (no: "2" with space-unit), `"TRL 1.6 GB"`, `"torch 2.1 MB"`, `"trl 50 %"`, `"trl 2x"`, `"vllm 0.5x"`, `"torch 2.13.0 x"` -> all `[]`. `"torch 2.13.0s"` -> `[]`. But `"torch 2.13.0 and"` -> `["2.13.0"]`, `"torch 2.13.0."` (sentence end) -> `["2.13.0"]` (trailing dot not part of version), `"torch 2.13.0,"` -> `["2.13.0"]`.
Non-mentions: `"torch is fast"`, `"torch.nn"`, `"torchvision 0.9"` for name `torch`, `"pytorch 2.1"`, `"trl-lib 3.2"`, `"trlx 1.0"`, `"vllm_native 1.0"`, `"the trl 1 case"`? (single integer: pin; recommend not a mention, version needs a dot), `"in torch 5 minutes"`, `"torch 3 times"`: pin: a bare integer after the name is not a mention unless `==`. `"torch 2.13 0"`. Digits-only dotted dates `"trl 2026.11.08"` is a mention (reported as non-allowed; document).
Why: the spec's unit list is `s, ms, GB, MB, %, x`; the `s` rule must not eat "s" of a word: `"torch 2.13.0 sources"` -> mention `["2.13.0"]` (unit must be a whole token), whereas `"vLLM 0.09 s"` -> `[]`. Test both sides of that.

**T40 two-part legality** (unit)
allowed `{1.17.1}`: mention `1.17` legal (allowed version starts with `1.17.`); `1.1` NOT legal (`1.17.1` starts with `1.1` but not `1.1.`); `1` not a mention; `1.17.1` legal; `1.17.10` NOT legal; `1.17.2` not legal; `1.170` not legal. With allowed `{2.13.0+cu130}` mention `2.13` legal, `2.13.0` legal (see T37), `2.1` not.
Why: prefix-matching without the dot lets `1.1` and `1.17.10` pass.

**T41 frozen and exempt contexts for check 7** (unit)
`docs/research/old.md` with `trl==1.16.0` -> not scanned. `experiments/phase0/a.py` with `trl==1.16.0` -> not scanned. `experiments/phase2/lab/x.py` -> not scanned. `experiments/phase2/RESULTS.md` -> not scanned. `experiments/phase2/docker/Dockerfile.vllm_native` with `torch 2.13.0` -> NOT reported (given T37 lock). `experiments/phase2/docker/Dockerfile.train` with `ARG TRL=1.16.0`? (no whitespace or `==`: not a mention). `src/objdel/train/x.py` with comment `# trl 1.16` -> reported. `tests/test_x.py` containing `"trl==1.16.0"` as a test literal -> reported (tests are in scan list); BUT `tests/test_docs.py` itself must not contain seeded bad versions as literals unless built at runtime (otherwise the check fails on its own source). Assert that: scan `tests/test_docs.py` -> `[]` mentions of non-allowed versions.
Why: self-referential failure is the most likely first CI break.

**T42 check 7 negative seeded** (unit): tmp tree with pyproject pinning trl 1.17.1 and `docs/a.md` "TRL 1.16" -> reported with file and line; fix to "TRL 1.17" -> clean.
**T43 check 7, B9 regression of today's numbers** (integration, before/after)
Before DOC: 8 strings "TRL 1.16" and 5 `blob/v1.16.0` links outside frozen paths (G5). After DOC: `grep -rn "TRL 1.16" . --include=*.md --include=*.py ...` outside frozen = 0 and `huggingface/trl/blob/v1.16.0/` = 0. Record the pre-count as a documented literal in the test docstring (8 and 5) so a reviewer sees the starting state. The test asserts only the after-state.
**T44 property: version_mentions is monotone under unrelated text** (property)
For random surrounding text with no package name, mentions == []; inserting `"trl==1.2.3"` into random text yields exactly that one mention (surrounding chars drawn from non-alnum to avoid creating `xtrl==`). Also `"xtrl==1.2.3"` and `"trl_==1.2.3"` -> no mention (word boundary).

---
## H. Cost checks (`cost_mentions`, B5.8)

**T45 cost_mentions, parametrize** (unit; expected list)
- `"| Alarm | $10.00 |"` -> `[("alarm", 10.0)]`
- `"| Hard stop | $17.50 |"` -> `[("hard stop", 17.5)]`
- `"| alarm | $10 |"` (case-insens first cell), `"| **Alarm** | $10 |"` (bold) -> pin: bold stripped -> `[("alarm",10.0)]`
- `"| Alarm | ~$10 |"`, `"| Alarm | USD 10 |"` -> pin (recommend: `$` only per spec) `[]` for USD; `~$10` -> 10.0
- `"| Alarm | $1,000.50 |"` -> `[("alarm", 1000.5)]`
- `"| Alarm | $10 (was $15) |"` -> first dollar only? Spec: "a dollar value in a ... row"; history wording is legal: pin `[("alarm",10.0)]` plus decide `$15` (recommend: all `$` values in the row are checked, so this row is reported; document that history must be in prose, not in the Alarm row). Test pins the chosen behaviour.
- `"The hard stop is $17.50."` -> `[("hard stop",17.5)]`; `"alarm = $10"`, `"hard stop of $17.5"`, `"Alarm at $10.00"` -> one each; `"hard stop (is) $17.5"` -> []
- `"rose from $15 to $17.50"` -> `[]` (B6 legal) ; `"the alarm rose from $15"` -> `[]`
- `"alarm is $10 and hard stop is $17.50"` -> both, in order
- `"| Cost | $10 |"`, `"| Budget | $3 |"` -> `[]` (first cell not Alarm/Hard stop)
- `"| Alarm hours | $10 |"` -> `[]` (first cell must be exactly Alarm)
- `"| Alarmed | $10 |"` -> `[]`
- `"alarms are $10"` -> `[]` (plural not in form)
- Multi-line: `"hard stop is\n$17.50"` -> pin: sentences are single-line -> `[]`; ensure no crash.
- `"| Alarm | 10 |"` (no `$`) -> `[]` and the check treats the runbook as "has row" (T47) with an unparsable value: pin as failure (value missing) rather than pass.
Why: this check pins money limits; false negatives let a wrong stop value ship.

**T46 check 8 vs policy** (integration + unit with fake module)
`POLICY_MODULE = "experiments.phase2.vast.policy"`; import it and read `ALARM=10.0`, `HARD_STOP=17.5`. Docs row `$10.00` equals 10.0; `$17.50` equals 17.5; `$17.5` equal; `$17.51` differs -> reported (float compare exact via `Decimal` or `==` after parse; assert `$10.00 == 10.0`). Substitute a fake module via monkeypatching the constant to `ALARM=12.0` and assert today's docs are now reported (proves the check reads the module, not a literal).
Why: the check must follow the constant, not today's 10/17.5.

**T47 runbook completeness** (unit + integration)
`docs/runbooks/gpu-run.md` missing the Alarm row -> reported; missing the Hard stop row -> reported; row present but in a different file -> still reported; rows in the wrong table (first cell "Alarm" in a prose list) -> reported. Real repo: both rows exist and equal policy.
**T48 B9 policy docstring** (unit): `re.search(r"\$\d", policy.__doc__) is None` (today fails: 4 amounts on 2 lines, `$10 alarm, $15 stop` and `$15`). Also assert the docstring no longer says `$15` anywhere. Note today's docstring says stop `$15` and alarm `$10` while code says HARD_STOP `17.5`: a stale-number bug this check class exists for. Regression: the docstring may mention `HARD_STOP` by name.
**T49 No doc publishes a plan phase cost above $3** (integration, soft): scan `docs/plan/*.md` and `docs/runbooks/*.md` for sentences with "phase" and `$N` where N>3 not marked history; expected: none, or each in the allowlist. (Spec constraint line.) Mark xfail-strict only if the owner confirms the wording is scannable; otherwise leave as a review checklist item.

---
## I. Drift checks 9 to 13

**T50 headline_strings exact output for r5** (unit; the G6 literal)
Input: the r5 headline with means `iou_tc 0.6449?`. NOTE: use the real file; the expected list is exactly (order matters):
`["IoU 0.645","boundary F 0.816","precision 0.695","recall 0.748","Dice 0.701","legal 99.4 %","shifted crops","IoU 0.618","F 0.806","precision 0.689","recall 0.716"]` (11 strings; no `ci_unit` string since `item`).
**T51 headline_strings rounding and format edges** (unit, parametrize)
Means (test.iou_tc): `0.6454` -> `IoU 0.645`; `0.6456` -> `IoU 0.646`; `0.0` -> `IoU 0.000`; `1.0` -> `IoU 1.000`; `0.0005` -> pin exact `f"{0.0005:.3f}"` result (`0.001`? binary repr gives `0.001`; use Python's own format as the contract: expected literal computed by hand, `"0.001"`), `0.9995` -> `"1.000"`? hand-check: `f"{0.9995:.3f}"` is `'0.999'` or `'1.000'` depending on binary repr; the test uses the literal observed on CPython 3.12 and a comment, to flag if the implementation switches to `round()` or `Decimal` (which would change it).
legal: `0.994` -> `legal 99.4 %`; `1.0` -> `legal 100.0 %`; `0.0` -> `legal 0.0 %`; `0.9995` -> `legal 100.0 %` or `99.9 %` (hand literal from `f"{0.9995*100:.1f}"`); `0.99449999` -> `99.4`. Note `0.994*100` = `99.4` exactly? (`99.39999999999999`; formats to `99.4`). Assert that.
**T52 view words** (unit): `lidc-test-shifted` -> includes `"shifted crops"`; `lidc-test` (centred) -> pin the word (spec only defines the shifted mapping; an unknown view must either map to a documented word or raise `KeyError`/`ValueError`, never produce an empty string that `"" in README` always satisfies). Add assertion: no returned string is empty or whitespace.
Why: an empty required string makes check 9 pass trivially.
**T53 ci_unit string** (unit): `ci_unit="item"` -> list has 11 items, no `"-level"`; `ci_unit="patient"` -> 12 items, last is `"patient-level"`; `ci_unit=""` or missing -> pin (raise).
**T54 headline_strings does not mutate input; deterministic; missing key -> KeyError** (unit): deep-copied input equal after the call; a headline lacking `test.metrics.dice_vote` raises KeyError (not silently skips).
**T55 check 9 negative** (unit): README with all 11 strings -> clean; README missing `"Dice 0.701"` -> reported; README containing `IoU 0.6451` ("IoU 0.645" is a prefix) -> counts as present (substring); README with `IoU 0.640` -> reported with the missing string. Check 9 runs against the real `docs/results/headline.json` and real README (G6). Also: `README.md` and the HF card use the same strings (T66).
**T56 check 10 CITATION version** (unit, tmp): `pyproject` `version = "0.1.0"` and CITATION `version: 0.1.0` / `version: "0.1.0"` / `version: '0.1.0'` -> clean (strip quotes); `0.1.1` -> reported; CITATION missing `version` -> reported; `cff-version: 1.2.0` must not be mistaken for `version` (line starting `cff-version`). pyproject `[project].version` must be read from the `[project]` table, not `[tool.x] version`. Use `tomllib`. Dynamic version -> clear failure.
**T57 check 11 archive READMEs** (unit): tree `archive/a/README.md` with the four parts (what it was, what replaced it, retiring commit, how to restore), `archive/b/` without README -> reported; `archive/c/README.md` missing "how to restore" -> reported; `archive/README.md` itself not required to have the four parts; `archive/.hidden/` ignored; a file directly in `archive/` (not a folder) ignored. How the four parts are detected (headings vs keywords) is not in the spec: pin from implementation, but test that a README with a 40-hex commit string passes the "commit" part and one with the word "commit" but no hash does not (recommended). Real repo: every `archive/*/` complies.
**T58 check 12 attribution** (unit, parametrize): each of `Co-Authored-By: Claude`, `co-authored-by: claude` (case), `Co-Authored-By: Claude Opus <noreply@anthropic.com>`, `Generated with Claude Code`, `Generated with [Claude Code]`(markdown link; pin whether caught; the spec says `Generated with Claude`, so the link form would slip through; test documents), in `.github/pull_request_template.md`, `docs/contributing.md`, `CITATION.cff` -> reported. Negative: a contributing guide line that says "never add Co-Authored-By trailers" (no `Claude`) -> clean; the repo owner's rule text "never add a `Co-Authored-By: Claude …` trailer" in `docs/contributing.md` would be reported (!): pin that the guide must phrase the rule without the literal string, e.g. test asserts the real contributing guide states the rule and passes. Missing file -> clean (not crash) except pr template required? pin.
**T59 check 13 nav coverage** (unit, tmp tree + mkdocs.yml)
`nav_pages(mkdocs.yml)`: yaml with nested nav
```
nav:
  - Home: README.md
  - Guides:
    - Start: getting-started.md
    - Sub:
      - X: guides/x.md
  - plan/Direction.md
  - External: https://example.org
```
Expected `{"README.md","getting-started.md","guides/x.md","plan/Direction.md"}` (relative to `docs_dir`, URL excluded). `exclude_docs` (gitignore syntax: `research/`, `*.draft.md`) honoured: an excluded file need not be in nav; a file under `reference/objdel/` is exempt; frozen `docs/research/**` is exempt via `exclude_docs`; a non-excluded `docs/orphan.md` -> reported (B6). Also: nav entry pointing at a file that does not exist -> pin (not check 13's job, but `zensical build --strict` catches it; assert the build catches it in T69). YAML with `!!python/name` tag in mkdocs.yml (Material uses `!!python/name:material.extensions.emoji.twemoji` in `pymdownx.emoji`): `nav_pages` must parse with a loader that tolerates unknown tags (custom `SafeLoader` subclass with multi-constructor), else `yaml.safe_load` raises `ConstructorError` on real mkdocs.yml. HIGH VALUE: include a test that parses a tmp mkdocs.yml with `!!python/name:` and `!ENV` tags.
**T60 check 13 on real repo** (integration): all of `docs/**/*.md` outside excluded/generated is in nav, includes `docs/plan/**` (4 files), `docs/results/runs.md`, `docs/reproducibility.md`.

---
## J. Headline writer (`experiments/phase2/headline.py`, B7)

Fixture helper `final_eval(split, view="x/y/lidc-test-shifted", means=...)` writes `runs/ID/eval/eval_{test,val}/final_eval.json` in `tmp_path` and calls `main` with a monkeypatched repo root / `--out` into tmp. (If `main` hard-codes `runs/`, change cwd via `monkeypatch.chdir`.) The real file schema must be read from `src/objdel/eval/*final_eval*` before the sketch is finalised; input shape in tests is whatever `final_eval` emits.

**T61 happy path writes the exact schema** (integration)
Both inputs valid; `main(["--run","r5","--model-ref","svenk029/objdel-sft@rl-full-3ep","--out",str(out)])` returns 0; JSON has keys `run, model_ref, ckpt, seed, K, ci_unit, test, val`; `test.source == "runs/r5/eval/eval_test/final_eval.json"`; `test.view == "lidc-test-shifted"` (last path part of `x/y/lidc-test-shifted`); `test.items` equals input items; metrics keys exactly `{iou_tc,f5,precision_vote,recall_vote,dice_vote,legal}` (extra keys in input dropped, `len==6`); each has `mean` and `ci95` list of 2; `ci_unit == "item"`.
**T62 numbers copied unrounded and bit-exact** (property)
Random floats in [0,1] with 17 significant digits as means and CI bounds: output `float` equals input `float` exactly (`==`, not `approx`); JSON round trip via `repr` preserves. Also NaN/inf in input -> refuse (non-zero, no file), since JSON `NaN` is invalid.
Why: B7 "copied unrounded"; `round(x,3)` bug would break the README string only on boundary cases.
**T63 refusals write nothing** (parametrize, unit)
(a) test file missing; (b) val file missing; (c) both missing; (d) test file `split == "val"` (swapped); (e) val file `split == "test"`; (f) `split` key absent; (g) file with no `view`; (h) `view` null or `""`; (i) invalid JSON / truncated file; (j) empty file; (k) metric key missing (e.g. no `dice_vote`); (l) `ci95` missing. Each: `main(...)` returns non-zero (or raises SystemExit non-zero) and `--out` path does NOT exist AND a pre-existing `--out` file is left byte-identical (not truncated).
Why: B7 states exit non-zero and writes nothing; partial or truncated headline silently weakens the README check.
**T64 atomic write and permissions** (unit): output written via temp + rename (monkeypatch `os.replace` to raise -> original file intact, no `.tmp` litter). Output ends with newline, `indent=2`, keys in schema order (stable diffs).
**T65 make headline** (integration, uses `make -n`): `make -n headline RUN=r5 MODEL_REF=x` prints a command containing `experiments.phase2.headline`, `--run r5`, `--model-ref x`, and `-I`; `make headline` (no RUN) and `make headline RUN= MODEL_REF=x` exit non-zero with a message and run nothing; `MODEL_REF` empty -> pin (recommend fail). Value with a space or shell metachar `RUN='a b'` or `RUN='r5; touch pwned'`: command must quote (no file `pwned` created). Use `make -n` and assert the quoted form.
Why: shell injection through make variables; owner runs it by hand.
**T66 headline round trip vs README check** (end-to-end, tmp)
Generate headline.json with `main` from fixture inputs; build a README containing all `headline_strings(h)`; check 9 passes; change one input mean by 0.0005 across a rounding boundary, regenerate, check 9 now fails until README is updated. Proves writer, formatter and checker agree on rounding (same function).
**T67 regression of today's numbers (r5)** (integration, skipped if `docs/results/headline.json` absent)
Load the tracked headline.json; assert the 11 literal strings of T50 equal `headline_strings`; assert `test.metrics.iou_tc.mean` within [0.6445, 0.6455), ci95 lower < mean < upper, `legal.mean` in [0.9935,0.9945); `val` view differs from `test` view; `test.items` is an int > 0; `model_ref` and `ckpt` non-empty; `seed` int; `K` int. These bounds come from the spec (G6) not from re-running.

---
## K. HF card (`card`, `card_sections`, `card_update`, `--card-only`, B8)

**T68 card() unchanged** (regression, golden)
Run `H.card(...)` on the existing test's fixture inputs (reuse from `tests/test_p2_hf_push.py`) before the change; save as `tests/goldens/hf_card_before.md`; after: byte-equality. Existing hf_push tests unchanged and pass (G12).
**T69 card_update on a card with no Evaluation/Limitations** (unit)
old = `"# T\n\nintro\n\n## Usage\n\ntext\n"`. Output: every line of `old` present in the same relative order; plus `## Evaluation` and `## Limitations` appended; contains 6 test means with CI text, view, items, K, seed, `ci_unit`, 4 val means. Literal checks: `"0.645"`, CI bounds formatted, `"99.4"`.
**T70 idempotent** (property + unit): `card_update(card_update(old,h),h) == card_update(old,h)` for: old empty `""`, old without trailing newline, old with Windows line endings `\r\n`, old with sections already present, old with Evaluation in the middle followed by other sections, old with both sections reversed, old with `## Evaluation` heading inside a fenced code block (must not be treated as the section), old with `### Evaluation results` (level 3, different heading, not replaced), old containing `## Evaluation` twice (pin: replace first, or all; second pass stable either way), hypothesis-generated random markdown lines for old (alphabet incl. `#`, newline, backtick) with the invariant plus "every non-section line of old kept".
**T71 other lines kept exactly** (property): remove the two sections from `card_update(old,h)` (by heading span) and compare to `old` minus its own old sections: equal line for line. Sections replaced in place keep position; a section ends at next `## ` heading or EOF, not at a blank line (multi-paragraph and table inside the section are replaced fully; no leftover old table rows).
**T72 section content from headline, not from old** (unit): old card's Evaluation says `IoU 0.500`; after update it says `0.645` and `0.500` is gone. Stale numbers never survive. 95 % CI rendered for the 6 test metrics (12 bounds) and ONLY means for val (4) (spec: "4 validation means"): assert val CI bounds do not appear (guard against over-publishing).
**T73 card_update does not leak RIGA / local paths / run internals** (unit): output contains no `runs/`, no `/[U]sers/`, no IPv4, no `hf_`-prefixed token pattern; `test.source` path is NOT required in card (pin).
**T74 card_sections(headline)** (unit): returns text beginning `## Evaluation` and containing `## Limitations`; deterministic; no trailing whitespace lines beyond one newline; the same function used by `card_update` (call both and assert `card_sections(h)` is a substring of `card_update("", h)` modulo leading newline).
**T75 --card-only CLI: token and privacy gates** (integration with `FakeApi`, style of test_p2_hf_push.py)
Cases (each asserts upload call list):
(a) token via stdin, `token_problems==[]`, repo private -> exactly one upload, path_in_repo `README.md`, revision `V`, content equals `card_update(<fetched old card>, headline)`;
(b) repo public -> no upload, non-zero exit;
(c) over-broad token (`role=write`) -> no upload, non-zero;
(d) empty stdin / whitespace-only token -> no upload, non-zero, token not echoed;
(e) `--run` given together with `--card-only` -> pin (error or ignored), no run files read;
(f) `--headline` file missing/invalid -> no upload, non-zero, BEFORE any network call (assert FakeApi.calls empty);
(g) `--revision` missing -> parser error;
(h) token must never appear in stdout/stderr or in any exception message (capsys: assert token string not in output, for success and failure);
(i) upload failure from API -> non-zero, no retry loop, no second upload;
(j) argv contains the token -> refused (token only from stdin).
(k) fetched old card absent (404) -> pin (create from `card_sections` only, or fail); never overwrites with an empty base silently.
(l) upload uses `create_commit` with single `CommitOperationAdd("README.md")` (assert exactly one operation, not a folder upload; no `LICENSE`, `*.safetensors`, `adapter*` touched).
Why: B8 sentence is the whole security contract of this mode.
**T76 token_problems unchanged** (regression): existing parametrised test passes unmodified.

---
## L. Docs build (B1, B2, B13), CI

**T77 build passes** (integration; skipped if `zensical` not installed; CI `docs` job runs it)
`zensical build --strict -c` exit 0 and output contains "No issues found". Run in a clean env without torch/trl/vllm/objdel importable (use `python -I -S`-like env with only docs lock): proves B2 "needs no torch, TRL, GPU or installed objdel". Assert `import objdel` fails inside that env (`pytest.importorskip` inverse) so the claim is tested, not assumed.
**T78 seeded build faults each break `--strict`** (integration on a tmp copy of docs): (a) broken internal link `[x](nope.md)`; (b) link to missing anchor `a.md#nope` (spec: `anchors: warn` + strict -> fail); (c) nav entry for a missing file; (d) file not in nav -> (only INFO unless configured; pin) ; (e) bad YAML in `mkdocs.yml`. Expected non-zero in a,b,c,e.
Why: B1 says the build fails on broken link or missing anchor; config typo `validation.links.anchors: ignore` would silently drop it.
**T79 mkdocs.yml static assertions** (unit, yaml with tolerant loader): `strict: true`; `validation.links.anchors == "warn"`; `docs_dir` is `docs`; `exclude_docs` contains `research/` (or equivalent), and does not contain `plan/`; nav includes `plan/Direction.md`, `plan/Architecture.md`, `plan/Research2.0.md`, `plan/previous-vs-new.md`; `mkdocstrings` handler lists `paths: [src]` and does not import (static analysis; `handlers.python.options` has no `load_external_modules` requiring import); no entry referencing `_gpucheck`.
**T80 generated API pages (B2)** (integration, after T77): site has `reference/objdel/{checker,data,eval,geometry,reward,secrets,targets,train}/index.html` (or the layout the config yields); NO file path containing `_gpucheck`; the string `_gpucheck` does not appear in `reference/objdel/*` HTML at all except none (search the site's `search/*.json` too).
**T81 secrets page shows no secret value** (integration): HTML text of the `secrets` page, run the repo's secret scanner regexes (HF `hf_[A-Za-z0-9]{30,}`, GitHub `gh[pous]_[A-Za-z0-9]{36,}`, `AKIA[0-9A-Z]{16}`, `-----BEGIN`, IPv4, `ssh -p \d+ root@`) -> no match. Additionally set a canary: monkeypatch a module-level constant in a tmp copy of `secrets` with `"hf_CANARY..."`; build; canary absent from the page iff the page does not render values (docstring only). If mkdocstrings renders module attributes with values, canary appears -> test fails -> config needs `show_if_no_docstring: false` / `members` filter.
**T82 Material 9.8.1 builds same mkdocs.yml** (integration, weekly/optional, skipped without it): `mkdocs build --strict` exit 0. Only Zensical-specific keys under a namespace Material ignores; test fails if the config uses a key Material rejects.
**T83 docs lock completeness (B13)** (unit): `experiments/phase2/docker/locks/docs.txt`: every non-comment line matches `^[A-Za-z0-9_.\-\[\]]+==[^=\s;]+` (exact pin, with `--hash` optional); includes `zensical==0.1.0` and `mkdocstrings-python==2.0.10`; contains transitive packages (`griffe`, `mkdocstrings`, `markdown`, `pymdown-extensions`, `jinja2`, `pyyaml` present; assert count >= 15 as a sanity floor and every `Requires-Dist` of the two top packages (read from the built env in CI only) is present: in CI job, `pip install -r docs.txt && pip check` exit 0 and `pip freeze` is equal as a set to the lock). `pyproject.toml` `docs` extra equals the two pins exactly. No unpinned `>=` lines. CI workflow file installs from the lock path (grep `docs.txt` in `.github/workflows/tests.yml` `docs` job), not from `.[docs]`.
**T84 workflow security** (unit on YAML): `links.yml`: every `uses:` for `lycheeverse/lychee-action` and `actions/cache` is a 40-hex SHA (comment may show `v2.9.0`); `on.pull_request.paths` includes `**/*.md`; has `schedule`; `permissions` is read-only (`contents: read`); no `pull_request_target`; weekly job `continue-on-error` for Material. `tests.yml` `docs` job has no secrets. Test uses a tolerant YAML loader (`on:` parses as `True` key in YAML 1.1: handle `data[True]`; that is a classic pitfall).
**T85 `.gitignore`** (unit): `git check-ignore` (or parse) says ignored: `site/index.html`, `.cache/x`, `.lycheecache`, `docs/research/paper.pdf`; NOT ignored: `docs/README.md`, `docs/results/headline.json`, `docs/plan/Direction.md`, `docs/figures/README.md`. `/site/` is anchored: `docs/site/x.md` NOT ignored.

---
## M. Repo state gates (G1, B10, B12, G7, G10, G11)

**T86 every doc file is tracked (G1)** (integration, `git ls-files`, skipped outside git): `git ls-files docs | grep -v ^docs/research/ | wc -l` equals `find docs -type f -not -path 'docs/research/*' | wc -l`; set equality of names, not only counts (counts can match with one extra and one missing). `docs/plan/` has exactly 4 tracked files named as spec; `docs/research/research2.0` does not exist; `docs/plan/previous vs new.md` (with spaces) does not exist. Today this fails (many `??` untracked): that is the point. Gate for owner, not for `make test` of a dirty working tree: mark as CI-only (`GITHUB_ACTIONS`) or `pytest.mark.skipif(dirty)`. Caveat: `tests/test_docs.py` itself must read only tracked files (B4): implement via `git ls-files` or a fixed tree walk limited to doc scope and test it with an untracked decoy `docs/zz_untracked.md` containing a broken link: check 1 must NOT report it when the decoy is untracked (and runs/ or data/ decoys never read: assert via monkeypatched `open` that no path under `runs/` or `data/` is touched; use `strace`-free approach: patch `Path.read_text`/`open` and record).
**T87 no moved-path references remain (B10)** (integration): grep over doc scope, `Makefile`, `README.md`, `pyproject.toml`, `.github/` for `research2.0` and `previous vs new` -> 0 hits outside frozen paths (frozen files may keep old links; check 1 does not run on them).
**T88 secret and address scan of untracked docs (B12)** (integration, owner-run script): over `docs/**` untracked files regex: `hf_[A-Za-z0-9]{20,}`, `gh[pousr]_[A-Za-z0-9]{30,}`, `github_pat_`, `ssh -p \d+ root@`, `\b(?:\d{1,3}\.){3}\d{1,3}\b` (allow `127.0.0.1`, `0.0.0.0`, version-like within code spans? pin an allowlist), `AGE-SECRET-KEY-`, `sk-[A-Za-z0-9]{20,}`. Expected no output. Seed test: a tmp file with each pattern is detected (9 detections), and `1.17.1`, `0.20.1`, `2026.11.08` are NOT flagged as IPv4 (they have 3 dots? `1.17.1` has 2; `10.0.0.1` has 3 flagged; `v1.2.3.4` ambiguity: flagged, document).
**T89 lint docstring gate (B3, G4)** (integration): `ruff check src/objdel --select D100,D101,D102,D103,D104,D106` -> 0 errors after DOC (63 today); `ruff check src tests experiments` clean; `pyproject` has `D` rules only for those 6 codes (parse `[tool.ruff.lint] select` and `extend-select`: assert no `D2xx`, `D4xx`, no bare `"D"`), `per-file-ignores` exempts `tests/**` and `experiments/**` from `D`, `extend-exclude` contains both `experiments/phase1/lab` and `archive`; `ignore` still has `E731`, `E741`; isort `known-first-party` retained. Seeded: a tmp package with a public class lacking a docstring, a nested public class, a public method, a `__init__.py` w/o docstring -> D101/D106/D102/D104 reported; a `_private` function and a `__init__` method (D107 not enabled) not reported; `tests/x.py` with undocumented function not reported (ruff run with the repo config on a tmp copy, `--config`).
Why: D105 (magic) and D107 (`__init__`) are deliberately NOT in the list; enabling `D` wholesale would flag thousands.
**T90 public/private boundary of docstring gate** (unit): `src/objdel/_gpucheck*` stubs or `.pyi`: not flagged; `__all__`-only re-exports in `__init__.py` need a docstring (D104).
**T91 test count (G7)** (CI step, not pytest): `pytest --collect-only -q` count >= 1,484 + new tests; no test id from the baseline list is missing (store baseline ids in the PR description or compare to `git show main:` list: diff of node ids shows only additions).
**T92 G12 no behaviour change** (integration): golden-file tests (`tests/goldens`), target-hash tests, reward and checker tests pass unchanged; `git diff --stat main -- src/objdel tests/goldens` shows only docstring-only changes in `src` (assert via AST: `ast.dump` of every changed `src/objdel/**/*.py` with docstrings stripped equals the main version; the strongest guard that "library behaviour does not change"). Implementation: for each changed file, parse both versions, remove `Expr(Constant(str))` first statements in module/class/def, compare dumps. Also assert no `src/objdel` edit happened between a preflight stamp and its teardown: preflight stamp hash of `src/objdel` equals current hash? (docstrings change the hash: so after DOC the stamp is stale; the test documents that a new preflight is required before the next rental: assert the preflight refuses a stale stamp).
**T93 G10 no AI attribution** (unit): scan the three files + `git log --format=%B main..HEAD` for `Co-Authored-By`, `Generated with`, `Claude`, `noreply@anthropic.com`, session URL lines (owner rule: none in repo files; note the harness attribution line `Claude-Session:` is separate and not covered by CLAUDE.md for commits; the test scans files, not commit trailers, to avoid contradicting the owner's rule).
**T94 G11 ADR index both ways on the real repo** (integration): same as T29 on real tree; includes `0012` in index with status `Proposed` or `Accepted`, with Date in ISO format `YYYY-MM-DD`.
**T95 G8 lychee** (CI-only, skipped locally): `.lychee.toml` parse: has `exclude` with the localhost/example/private HF URLs; `max_retries` set; `accept` includes 200,429 (rate limits) and not 404; `exclude_path` contains `docs/research/` only if frozen files may keep dead links (pin: the spec says 0 broken outside the exclude list). Seeded: a tmp md with `https://example.invalid/x` fails lychee (offline: use `--offline` and a relative broken link to be hermetic).

---
## N. Cross-cutting adversarial, performance, determinism

**T96 CRLF, BOM, tabs, non-UTF8** (unit): every `docs_check` function given text with `\r\n`, a leading `﻿`, tab-indented code, and a file containing bytes `\xff` (non-UTF8) via `Path`: either decode with `errors="replace"` or raise a clear error naming the file; never `UnicodeDecodeError` without file name. Link `[a](b.md)\r\n` -> `b.md` (no `\r`).
**T97 unicode paths and names** (unit): `` `docs/guides/café.md` `` exists (NFC) while disk has NFD (macOS) -> pin: normalise to NFC both sides; link `%C3%A9` decodes.
**T98 huge doc** (performance gate): `doc_files + all checks` on the real repo complete within 5 s CPU (target: well under 2 s) on the Mac; `tests/test_docs.py` whole module < 10 s; a 5 MB synthetic md processed by each of the 7 pure functions < 1 s each. The build gate: `zensical build --strict` < 60 s on CI (no hard assertion, log only; fail > 300 s as a hang guard).
**T99 determinism** (property): every function returning a list returns a sorted or document-order list, stable across two calls and across `PYTHONHASHSEED=0,1` (subprocess run twice, compare outputs). Sets are compared as sets.
**T100 failure messages name file and line** (unit): each check's failure message contains the relative path and the offending string; asserts via `pytest.raises(AssertionError, match=...)` or the returned list. Tests that a failing run prints ALL offenders (not just the first): seed 3 faults across 3 files, expect 3 reported.
**T101 checks do not mutate the tree** (unit): hash tmp tree before/after running every check; equal.
**T102 checks are hermetic** (unit): run `tests/test_docs.py` with network disabled (monkeypatch `socket.socket` to raise) and with `runs/` and `data/` replaced by directories that raise on access (`chmod 000` on POSIX): all pass. Satisfies B4.
**T103 GPU/CUDA none** : this phase has no GPU/bit-exact CPU vs GPU test (B11: no golden, reward or checker result changes). The single bit-equality-style gate is T92 (AST equality except docstrings) plus the unchanged goldens: `pytest tests/test_checker_oracle.py tests/test_gpu_scoring.py tests/test_cuda_checker.py -q` pass; the GPU variant `test_cuda_checker_gpu.py` stays skipped on the Mac (assert it still collects).

---
## Pytest sketches for the 5 most important tests

Ranking: (1) check 7 version regex (T38-T41) as it has the most exact edge rules; (2) card_update idempotence and line preservation (T69-T71); (3) headline refusal writes nothing (T63); (4) link and backtick negatives with legal inputs (T09, T16-T18, B6); (5) check 8 policy cost vs docs (T45-T48). Helper `tree` shown once.

```python
"""test_docs sketches (set 1). Pure functions come from tests/docs_check.py."""
import json
import re
from pathlib import Path

import pytest
from hypothesis import given, strategies as st

from tests import docs_check as D          # exact import style: follow tests/conftest.py
from experiments.phase2 import hf_push as H
from experiments.phase2 import headline as HL

def tree(root: Path, files: dict[str, str]) -> Path:
    for rel, text in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return root

# ---- 1. version mentions: spec literals, both sides of every rule (T38-T41)
@pytest.mark.parametrize("name, text, expected", [
    ("trl", "TRL 1.16", ["1.16"]),
    ("trl", "pin trl==1.16.0 now", ["1.16.0"]),
    ("trl", "https://github.com/huggingface/trl/blob/v1.16.0/trl/x.py", ["1.16.0"]),
    ("torch", "torch 2.13.0 and more", ["2.13.0"]),
    ("torch", "torch 2.13.0.", ["2.13.0"]),
    ("vllm", "vLLM 0.09 s", []),                  # B6 legal: unit follows
    ("vllm", "vllm 0.09s", []),
    ("torch", "torch 2.1 GB", []),
    ("trl", "trl 50 %", []),
    ("torch", "torch 2.13.0 sources", ["2.13.0"]),  # 's' is a unit only as a whole token
    ("torch", "torchvision 0.99", []),
    ("trl", "xtrl==1.2.3", []),
    ("trl", "trl is fast", []),
])
def test_version_mentions(name, text, expected):
    assert D.version_mentions(text, name) == expected

def legal(mention: str, allowed: set[str]) -> bool:   # mirrors B5.7 wording, written independently
    return mention in allowed or any(a.startswith(mention + ".") for a in allowed)

@pytest.mark.parametrize("m, ok", [("1.17.1", True), ("1.17", True), ("1.1", False),
                                   ("1.17.10", False), ("1.16", False), ("1.170", False)])
def test_two_part_legality(m, ok):
    assert legal(m, {"1.17.1"}) is ok

def test_frozen_and_dockerfile_inputs_not_reported(tmp_path):
    root = tree(tmp_path, {
        "pyproject.toml": '[project]\nname="x"\nversion="0.1.0"\ndependencies=["trl==1.17.1"]\n',
        "experiments/phase2/docker/locks/vllm_native.txt": "torch==2.13.0\nvllm==0.20.1\n",
        "docs/research/old.md": "trl==1.16.0\n",
        "experiments/phase2/docker/Dockerfile.vllm_native": "# torch 2.13.0\n",
        "docs/a.md": "TRL 1.16\n",
    })
    bad = D.version_offenders(root)               # name per implementation; returns [(path, line, text)]
    assert [b[0] for b in bad] == ["docs/a.md"]

# ---- 2. card_update: kept lines, idempotent, stale numbers gone (T69-T72)
HEAD = {"test": {"view": "a/b/lidc-test-shifted", "items": 10,
        "metrics": {k: {"mean": 0.645, "ci95": [0.6, 0.7]}
                    for k in ["iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote", "legal"]}},
        "val": {"view": "a/b/lidc-val", "items": 5,
                "metrics": {k: {"mean": 0.618, "ci95": [0.5, 0.7]}
                            for k in ["iou_tc", "f5", "precision_vote", "recall_vote", "dice_vote", "legal"]}},
        "run": "r5", "model_ref": "m", "ckpt": "c", "seed": 0, "K": 8, "ci_unit": "item"}

OLDS = ["", "# T\n\nintro\n", "# T\n\n## Evaluation\nIoU 0.500\n\n| a | b |\n|--|--|\n| 1 | 2 |\n\n## Usage\nu\n",
        "# T\r\n## Limitations\r\nold\r\n", "```\n## Evaluation\n```\n", "no newline"]

@pytest.mark.parametrize("old", OLDS)
def test_card_update_idempotent_and_stale_free(old):
    once = H.card_update(old, HEAD)
    assert H.card_update(once, HEAD) == once
    assert "IoU 0.500" not in once and "| 1 | 2 |" not in once
    kept = [l for l in old.splitlines() if l.strip()]
    assert all(l in once.splitlines() for l in kept if "0.500" not in l and l not in ("old", "| a | b |", "|--|--|", "| 1 | 2 |"))

@given(st.lists(st.text(alphabet="ab #`\n", max_size=20), max_size=12).map("\n".join))
def test_card_update_property(old):
    once = H.card_update(old, HEAD)
    assert H.card_update(once, HEAD) == once

# ---- 3. headline refusals write nothing (T63)
def write_eval(root, run, split, **over):
    d = root / "runs" / run / "eval" / f"eval_{split}"
    d.mkdir(parents=True)
    doc = {"split": split, "view": "x/y/lidc-test-shifted", "items": 3, "metrics": {}}  # use real final_eval shape
    doc.update(over)
    (d / "final_eval.json").write_text(json.dumps(doc))

@pytest.mark.parametrize("break_it", [
    "missing_test", "missing_val", "swapped_split", "no_split", "no_view", "bad_json", "empty"])
def test_headline_refuses_and_writes_nothing(tmp_path, monkeypatch, break_it):
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "out.json"
    out.write_text("PRE")                        # pre-existing file must survive byte-identical
    if break_it != "missing_test":
        write_eval(tmp_path, "r", "test", **({"split": "val"} if break_it == "swapped_split" else {}))
    if break_it != "missing_val":
        write_eval(tmp_path, "r", "val")
    # ... apply remaining breakages by rewriting the file (no_split pops the key, no_view pops 'view', ...)
    rc = HL.main(["--run", "r", "--model-ref", "m", "--out", str(out)])
    assert rc != 0
    assert out.read_text() == "PRE"
    assert list(tmp_path.glob("*.tmp")) == []

# ---- 4. link/backtick negatives and B6 legal inputs (T09, T16-T18)
def test_broken_links_exact(tmp_path):
    root = tree(tmp_path, {
        "README.md": "[a](docs/x.md#top) [b](docs/missing.md) [c](docs/my%20file.md) "
                     "[d](<docs/my file.md>) [e](#local) [f](https://x.org) [g](mailto:a@b.c)\n",
        "docs/x.md": "x", "docs/my file.md": "y"})
    assert D.broken_links(root / "README.md", root) == ["docs/missing.md"]

def test_retired_path_legal_only_in_archive_readme(tmp_path):
    root = tree(tmp_path, {"archive/old/README.md": "see `experiments/phase0/gone.py`\n",
                           "docs/g.md": "see `experiments/phase0/gone.py`\n"})
    spans = lambda p: [s for s in D.backtick_paths((root / p).read_text()) if not (root / s).exists()]
    assert spans("docs/g.md") == ["experiments/phase0/gone.py"]
    # archive README text is returned by backtick_paths too; the CHECK skips archive/**/README.md:
    assert D.check_backtick_paths(root) == ["docs/g.md: experiments/phase0/gone.py"]

@pytest.mark.parametrize("span, returned", [
    ("`src/objdel/a.py`", True), ("`docs/plan/`", True), ("`src/*.py`", False), ("`src/<p>/a`", False),
    ("`src/{a,b}`", False), ("`src/a b`", False), ("`make test`", False), ("`runs/x.json`", False),
    ("`src`", False), ("`srcs/x`", False)])
def test_backtick_filter(span, returned):
    assert bool(D.backtick_paths(span)) is returned

# ---- 5. cost check follows the policy constants (T45-T48)
@pytest.mark.parametrize("text, expected", [
    ("| Alarm | $10.00 |", [("alarm", 10.0)]),
    ("| Hard stop | $17.50 |", [("hard stop", 17.5)]),
    ("The hard stop is $17.50.", [("hard stop", 17.5)]),
    ("alarm = $10", [("alarm", 10.0)]),
    ("rose from $15 to $17.50", []),
    ("| Budget | $3 |", []),
    ("| Alarm hours | $10 |", [])])
def test_cost_mentions(text, expected):
    assert D.cost_mentions(text) == expected

def test_policy_docstring_has_no_dollar_amount():
    from experiments.phase2.vast import policy
    assert re.search(r"\$\d", policy.__doc__) is None

def test_runbook_amounts_equal_policy(monkeypatch):
    from experiments.phase2.vast import policy
    text = Path("docs/runbooks/gpu-run.md").read_text()
    got = dict(D.cost_mentions(text))
    assert got["alarm"] == policy.ALARM and got["hard stop"] == policy.HARD_STOP
    monkeypatch.setattr(policy, "ALARM", 12.0)           # the check must follow the module
    assert got["alarm"] != policy.ALARM
```

---
## Z. Notes for the owner and the implementer (open points the spec leaves undecided)
1. `make X` mentions in English prose ("make sure"): spec silent; T24 forces a decision.
2. `file:line` suffix in backtick paths (T16): spec silent; real docs use it.
3. `+cu130` local version labels vs `torch 2.13.0` (T37): B6 legal case needs normalisation.
4. Value lists in an Alarm row with a history clause (T45).
5. Views other than `lidc-test-shifted` (T52): the empty-string trap makes check 9 vacuous.
6. YAML `!!python/name` and `!ENV` tags in `mkdocs.yml` crash `yaml.safe_load` (T59) and `on:` parses as `True` (T84).
7. `tests/test_docs.py` contains seeded bad literals (T41) and the contributing guide states the attribution rule (T58): both can trip their own check.
8. G1/G7/G8 are repo-state or CI gates; run them as owner or CI steps, not inside the dirty-tree `make test`.
