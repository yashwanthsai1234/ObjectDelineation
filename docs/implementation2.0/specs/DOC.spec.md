# Spec: Phase DOC (documentation and repo structuring)

REPO = the ObjectDelineation repository root. All paths are relative to REPO. Revision 2 (after
`docs/implementation2.0/reviews/DOC_phase_review.md`).

## Goal

Every engineering document is tracked in git, builds into one site with no warning, and is checked by tests against
the code. The headline test and validation numbers come from one tracked results file. Library behaviour does not
change.

## Definitions

- Frozen paths: `docs/research/**`, `experiments/phase0/**`, `experiments/phase1/**`, `experiments/**/lab/**`,
  `experiments/**/impl/**`, `experiments/cuda_checker/**`, `experiments/*/PLAN*.md`, `experiments/*/JUDGE_*.md`,
  `experiments/*/V0_V1.md`, `experiments/*/RESULTS*.md`, `experiments/*/EXECUTION.md`, `experiments/*/RESEARCH.md`.
  A file is frozen by its path, not by a banner.
- Doc scope: `README.md`, `ROADMAP.md`, `docs/**/*.md` except frozen paths (includes `docs/plan/**`),
  `experiments/README.md`, `archive/**/README.md`.
- Allowed versions of a package: every exact `name==X` pin in `pyproject.toml` and in
  `experiments/phase2/docker/locks/*.txt`.

## Required behaviours

B1. `make docs` runs `zensical build --strict` with `mkdocs.yml` at REPO root and exits 0. The build fails on a broken
    internal link or a missing anchor under `docs/` (`validation.links.anchors: warn` with `strict: true`). The same
    `mkdocs.yml` also builds with Material for MkDocs 9.8.1 `mkdocs build --strict`, exit 0.
B2. The built site has generated API pages under `reference/objdel/` for each public subpackage and module of
    `src/objdel` (checker, data, eval, geometry, reward, secrets, targets, train), and no page for the private
    extension `objdel._gpucheck`. The `secrets` page shows no secret value. The build needs no torch, TRL, GPU or
    installed `objdel`.
B3. `make lint` fails when a public module, class, method, function or nested class in `src/objdel` has no docstring
    (ruff D100, D101, D102, D103, D104, D106). `tests/` and `experiments/` are exempt from all `D` rules. Other `D`
    rules are not enabled. Existing ruff settings (`ignore`, `extend-exclude`, isort) are kept; `archive` is added to
    `extend-exclude`.
B4. `make test` includes `tests/test_docs.py`. It reads only tracked files. It never reads `runs/` or `data/`.
B5. `tests/test_docs.py` fails when any of these is true:
    1. In a doc-scope file outside `docs/`, a relative Markdown link target does not exist (after URL-decoding `%20`,
       accepting `<...>` targets, and dropping `#anchor`).
    2. In a doc-scope file other than `archive/**/README.md`, a backtick span that starts with `src/`, `experiments/`,
       `tests/`, `cuda/`, `docs/` or `archive/` and has no space, `*`, `<` or `{` names a path that does not exist.
    3. A `make X` mention in doc scope names no target in `Makefile`.
    4. A file `docs/architecture/decisions/NNNN-*.md` is missing from the index table in
       `docs/architecture/decisions/README.md`, or an index row has no file. Files starting with `_` are ignored.
    5. An ADR lacks any of: lines starting `**Status:**`, `**Date:**`, `**Deciders:**`; headings `## Context`,
       `## Decision`, `## Options considered`, `## Trade-off analysis`, `## Consequences`, `## Action items`.
    6. A file `experiments/**/configs/*.yaml` is not named in `docs/reference/configs.md`.
    7. A version mention for `trl`, `transformers`, `torch` or `vllm` in `src/`, `tests/`, `experiments/**/*.py`,
       `experiments/**/Dockerfile*`, `experiments/**/configs/*.yaml` or doc scope is not an allowed version (a
       two-part mention is legal when an allowed version starts with it plus "."). Frozen paths are skipped. A mention
       is the package name followed by `==` or whitespace and a version, not followed by a unit (`s`, `ms`, `GB`,
       `MB`, `%`, `x`); for trl and transformers also a `github.com/huggingface/<name>/blob/vX.Y.Z/` URL.
    8. In doc scope, a dollar value in a Markdown table row whose first cell is "Alarm" or "Hard stop", or in a
       sentence of the form "hard stop|alarm (is|=|of|at) $X", differs from `ALARM` or `HARD_STOP` in the policy
       module named by the constant `POLICY_MODULE` (today `experiments.phase2.vast.policy`); or
       `docs/runbooks/gpu-run.md` has no Alarm row or no Hard stop row. History wording ("rose from $15") is legal.
    9. `README.md` lacks any string returned by `headline_strings(headline)` for `docs/results/headline.json`.
    10. `CITATION.cff` `version` differs from `pyproject.toml` `[project].version`.
    11. A direct subfolder of `archive/` has no `README.md`, or its README lacks the four parts: what it was, what
        replaced it, the retiring commit, how to restore.
    12. `.github/pull_request_template.md`, `docs/contributing.md` or `CITATION.cff` contains
        `Co-Authored-By: Claude` or `Generated with Claude`.
    13. A doc-scope file under `docs/` that is not excluded by `exclude_docs` and not under `docs/reference/objdel/`
        is missing from `nav` in `mkdocs.yml`.
B6. Checks 1, 2, 3, 4, 7, 8, 11 and 13 each have a negative test on seeded bad input in `tmp_path` that the check
    reports. These inputs must NOT be reported: a retired path in an `archive/**/README.md` (check 2);
    `torch 2.13.0` in a `Dockerfile.vllm_native`, `vLLM 0.09 s`, and `trl==1.16.0` under a frozen path (check 7);
    "rose from $15" (check 8).
B7. `python -I -m experiments.phase2.headline --run ID --model-ref REF [--out PATH]` reads
    `runs/ID/eval/eval_test/final_eval.json` and `runs/ID/eval/eval_val/final_eval.json` and writes
    `docs/results/headline.json` (or PATH). It exits non-zero and writes nothing if either file is missing, a file's
    `split` is not the expected one (`test`, `val`), or a file has no `view`. `make headline RUN=ID MODEL_REF=REF`
    calls it and fails when `RUN` is empty.
B8. `hf_push.card(...)` returns exactly what it returns today. `hf_push.card_update(old, headline)` returns `old` with
    an "Evaluation" section (6 test means and 95 % CIs, view, items, K, seed, `ci_unit`, 4 validation means) and a
    "Limitations" section replaced or appended; every other line of `old` is kept; applying it twice gives the same
    text. `hf_push` with `--card-only --repo R --revision V --headline PATH` reads the token from stdin, uploads
    nothing unless `token_problems` is empty and the repo is private, and then uploads exactly one file,
    `README.md`, to revision V.
B9. No "TRL 1.16" and no `huggingface/trl/blob/v1.16.0/` link remains outside frozen paths.
    `re.search(r"\$\d", policy.__doc__)` is `None` for the policy module.
B10. Every file under `docs/` outside `docs/research/` is tracked in git. The plan lives in `docs/plan/`
     (`Direction.md`, `Architecture.md`, `Research2.0.md`, `previous-vs-new.md`), all tracked;
     `docs/research/research2.0/` no longer exists; every link to the old paths is updated. `/site/`, `/.cache/`,
     `.lycheecache` and `docs/research/*.pdf` are git-ignored.
B11. No golden-file, target-hash, reward or checker result changes.
B12. Before the owner's first `git add` of `docs/`, a scan of the untracked `docs/` files for HF tokens, GitHub
     tokens, `ssh -p N root@` lines and IPv4 addresses prints nothing, and each RIGA-derived figure or table has the
     owner's OK.
B13. Docs CI installs docs packages from `experiments/phase2/docker/locks/docs.txt` (all transitive packages pinned).

## Interfaces

Files created: `mkdocs.yml`, `CITATION.cff`, `.lychee.toml`, `.github/workflows/links.yml`,
`.github/pull_request_template.md`, `tests/docs_check.py`, `tests/test_docs.py`, `tests/test_headline.py`,
`experiments/phase2/headline.py`, `experiments/phase2/docker/locks/docs.txt`, `experiments/README.md`,
`archive/README.md`, `docs/results/headline.json`, `docs/results/runs.md`, `docs/reproducibility.md`,
`docs/runbooks/hf-checkpoints.md`, `docs/runbooks/ci.md`, `docs/architecture/status.md`, `docs/figures/README.md`,
`docs/architecture/decisions/_template.md`, `docs/architecture/decisions/0012-docs-site-with-zensical.md`,
`docs/research/roadmap-2026-10-18.md`.

Files moved: `docs/research/research2.0/*` -> `docs/plan/*` (`previous vs new.md` -> `previous-vs-new.md`).

Make targets: `docs`, `docs-serve`, `headline` (variables `RUN`, `MODEL_REF`); all three in `.PHONY`.

`pyproject.toml` extra: `docs = ["zensical==0.1.0", "mkdocstrings-python==2.0.10"]`.

CI: job `docs` in `.github/workflows/tests.yml` installs from the docs lock and runs the strict build;
`tests/test_docs.py` runs in the existing `python` job. Workflow `links.yml` runs lychee-action v2.9.0 pinned by
commit SHA, with `.lycheecache` restored and saved by SHA-pinned `actions/cache`, on pull requests touching
`**/*.md` and weekly; the weekly run blocks nothing and also runs a non-blocking Material 9.8.1 strict build.

`tests/docs_check.py` (pure functions; `Path` arguments are absolute or relative to REPO):

```python
FROZEN_GLOBS: list[str]
def is_frozen(rel: str) -> bool: ...
def doc_files(root: Path) -> list[Path]: ...
def outside_docs_dir(files: list[Path]) -> list[Path]: ...
def relative_links(text: str) -> list[str]: ...          # excludes http(s):, mailto:, pure #anchor
def broken_links(md: Path, root: Path) -> list[str]: ...
def backtick_paths(text: str) -> list[str]: ...
def make_targets(makefile: Path) -> set[str]: ...
def make_mentions(text: str) -> set[str]: ...
def adr_files(decisions_dir: Path) -> list[Path]: ...     # NNNN-*.md only, not _template.md, not README.md
def adr_index(readme: Path) -> set[str]: ...              # file names linked in the index table
def adr_headings(md: Path) -> set[str]: ...
def config_files(root: Path) -> list[Path]: ...
def allowed_versions(root: Path, name: str) -> set[str]: ...
def version_mentions(text: str, name: str) -> list[str]: ...  # "TRL 1.16", "trl==1.16.0", ".../trl/blob/v1.16.0/"
def cost_mentions(text: str) -> list[tuple[str, float]]: ...   # ("alarm" | "hard stop", dollars)
def headline_strings(headline: dict) -> list[str]: ...
def nav_pages(mkdocs_yml: Path) -> set[str]: ...
```

`headline_strings(h)` returns, in order, with `t = h["test"]["metrics"]`, `v = h["val"]["metrics"]`:
`f"IoU {t['iou_tc']['mean']:.3f}"`, `f"boundary F {t['f5']['mean']:.3f}"`,
`f"precision {t['precision_vote']['mean']:.3f}"`, `f"recall {t['recall_vote']['mean']:.3f}"`,
`f"Dice {t['dice_vote']['mean']:.3f}"`, `f"legal {t['legal']['mean'] * 100:.1f} %"`, the view words of
`h["test"]["view"]` (`lidc-test-shifted` -> `"shifted crops"`), `f"IoU {v['iou_tc']['mean']:.3f}"`,
`f"F {v['f5']['mean']:.3f}"`, `f"precision {v['precision_vote']['mean']:.3f}"`,
`f"recall {v['recall_vote']['mean']:.3f}"`, and `f"{h['ci_unit']}-level"` only when `ci_unit != "item"`.

`docs/results/headline.json` schema:

```json
{"run": "str", "model_ref": "str", "ckpt": "str", "seed": 0, "K": 0, "ci_unit": "item",
 "test": {"source": "runs/<run>/eval/eval_test/final_eval.json", "view": "str", "items": 0,
          "metrics": {"iou_tc": {"mean": 0.0, "ci95": [0.0, 0.0]}, "f5": {...}, "precision_vote": {...},
                      "recall_vote": {...}, "dice_vote": {...}, "legal": {...}}},
 "val":  {"source": "runs/<run>/eval/eval_val/final_eval.json", "view": "str", "items": 0,
          "metrics": {"same 6 keys": "..."}}}
```

Means and CI bounds are copied from the input unrounded. Exactly these 6 metric keys per split. `view` is the last
path part of the input's `view`. `ci_unit` is `"item"` for `final_eval`'s item bootstrap.

`experiments/phase2/headline.py`: `main(argv: list[str] | None = None) -> int` (0 on success, non-zero on refusal).

`experiments/phase2/hf_push.py`:
`card(base: str, meta: dict, man: dict, sel: dict, run: Path) -> str` (unchanged);
`card_sections(headline: dict) -> str`; `card_update(old: str, headline: dict) -> str`;
`main()` accepts `--card-only --repo R --revision V --headline PATH` (no `--run` in this mode).

## Constraints

- No GPU, no rental, no spend. All work runs on the CPU (Mac and CI); the card update runs on the Mac.
- The agent runs no git or HF write. The owner runs every commit and the card update.
- No `Co-Authored-By` or "Generated with" lines in commits, PR template or docs.
- Python that reads `runs/` or `data/` runs with `-I`. `tests/test_docs.py` reads neither.
- Frozen Markdown records get at most a one-line banner at the top; no other edit. No banner in `.py` or
  `Dockerfile` files.
- An accepted ADR changes only its Status line and its action-item check boxes.
- No edit under `src/` between a preflight and its teardown (the preflight stamp hashes `src/objdel`).
- Docs never write a number that a results file or a code constant owns without naming the source. No doc publishes
  a phase cost above the owner's $3 GPU budget as a plan.
- Stdlib plus already-pinned dev packages (`pyyaml`) for tests; no new test dependency.
- `make test` and `make lint` pass after every task.

## Acceptance gates

| # | Gate |
|---|---|
| G1 | `git ls-files docs \| grep -v ^docs/research/ \| wc -l` equals `find docs -type f -not -path 'docs/research/*' \| wc -l`; `docs/plan/` has 4 of 4 files tracked; `docs/research/research2.0/` does not exist. |
| G2 | `zensical build --strict -c` exit 0 with "No issues found" and anchors validated; `mkdocs build --strict` (Material 9.8.1) exit 0 on the same `mkdocs.yml`. |
| G3 | All 13 checks in B5 pass on the repo; each negative test in B6 detects its seeded fault and each listed legal input is not reported. |
| G4 | `ruff check src/objdel --select D100,D101,D102,D103,D104,D106` reports 0 errors (63 today). |
| G5 | B9 holds (today: 8 "TRL 1.16" strings, 5 `blob/v1.16.0` links, 4 dollar amounts on 2 lines of the `policy.py` docstring). |
| G6 | README contains all `headline_strings` of `docs/results/headline.json` (for r5: "IoU 0.645", "boundary F 0.816", "precision 0.695", "recall 0.748", "Dice 0.701", "legal 99.4 %", "shifted crops", "IoU 0.618", "F 0.806", "precision 0.689", "recall 0.716"). |
| G7 | `make test` and `make lint` pass; collected tests = baseline recorded before DOC1 (1,484 on 2026-11-08) + new tests; no existing test removed. |
| G8 | lychee on `main` reports 0 broken links outside the exclude list. |
| G9 | `CITATION.cff` version equals `pyproject.toml` version (0.1.0 today). |
| G10 | 0 AI-attribution strings in the PR template, contributing guide and CITATION. |
| G11 | 0 ADR files missing from the index; 0 index rows without a file. |
| G12 | Golden-file and target-hash tests pass unchanged; `hf_push.card()` output unchanged. |
| G13 | B12 holds. |
