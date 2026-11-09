# Contributing

These rules come from `CLAUDE.md`, the `Makefile`, `pyproject.toml`, the CI workflow and the owner's instructions for coding assistants.
Read [Getting started](getting-started.md) first.

## 1. The rules in short

1. `make test` and `make lint` must pass before a commit.
2. Library code lives in `src/objdel` and needs tests in `tests/`.
3. Experiment scripts live in `experiments/<phase>/` and run as modules in an isolated interpreter.
4. Python that reads `data/` runs with `-I`. Downloaded data is untrusted.
5. Commits are authored by the repository owner only. No `Co-Authored-By` trailer. No AI attribution line.
6. A pull request description has no "Generated with Claude Code" line or similar.
7. Git commands that change anything, and every `gh` command that changes a pull request or an issue, are run by the owner. Assistants give the exact commands and wait.

## 2. Workflow

```text
  branch  ->  change  ->  test  ->  lint  ->  commit  ->  pull request  ->  CI  ->  merge
  (owner)     src +       make       make     (owner)     (owner)           green
              tests       test       lint
```

1. Work on a branch. The main branch is `main`.
2. Change code and tests together. A behaviour change changes its [document](README.md) in the same pull request.
3. Run:

   ```bash
   make test
   make lint
   ```

   `make fmt` fixes lint issues that ruff can fix.
4. Run `make cov` when you add a module. The fast suite must keep at least 90 % line coverage of `src/objdel`. CI fails under 90 %.
5. If you changed `cuda/` or `CMakeLists.txt`, `make test` rebuilds the extension. For a GPU change also run `make cuda` on a CUDA box and the GPU test group ([Testing](testing.md)).
6. The owner commits and pushes. CI runs the jobs `python` (three operating systems), `train`, `gates` and `wheels`.

### Git rules that apply to everyone

- Do not run `commit`, `push`, `merge`, `rebase`, `reset`, `checkout`, `stash`, `tag` or branch creation unless the owner says to.
- Read-only commands are fine without asking: `status`, `log`, `diff`, `show`, `fetch`, `gh pr view`, `gh search`.
- "Fix the issues" is not permission to push or to change GitHub.
- Commit messages in this repository are one descriptive line with the numbers that matter. Example: `Shared-prefix forward: each prompt runs once per group of 8 answers (R5 on an RTX 5090: 31.0 -> 23.6 s/step, same rewards) ...`.

## 3. Where things go

| What | Where |
|---|---|
| Library code | `src/objdel/<package>/` and tests in `tests/test_<area>.py` |
| Slow exact references | `tests/oracles/` |
| Recorded results of retired code | `tests/goldens/` |
| Experiment scripts | `experiments/<phase>/<script>.py`. Run as `.venv/bin/python -I -m experiments.<phase>.<script>` |
| Run configs | `experiments/phase2/configs/*.yaml` |
| Results of a phase | `experiments/<phase>/RESULTS.md`. Frozen sections are never edited. Corrections are dated notes below them |
| Research reports | `docs/research/reports/` |
| Research notes | `docs/research/notes/` |
| Engineering documentation | `docs/` (this set) |
| Plan | `ROADMAP.md` |
| Layout and conventions | `README.md` |
| Downloaded data | `data/` (git-ignored) |
| Retired code | `archive/` (not built, not imported, not tested) |

Phase scripts under `experiments/` are frozen records of what produced published numbers. They are linted lightly (`E501`, `B023`, `B905` are ignored there).
`experiments/phase1/lab` is frozen evidence and is excluded from ruff.

## 4. Code conventions

- Python 3.12. Line length 120. Ruff rules `E`, `F`, `I`, `B`, `UP`. `E731` and `E741` are ignored. `objdel` and `experiments` are first-party imports.
- Pins are exact in `pyproject.toml` on purpose. Reward and target hashes are recorded against them. Do not loosen a pin without re-running the hash checks.
- Docstrings state the contract: inputs, outputs, errors, and the reasons for non-obvious choices. Many name an issue id (`B2-7`, `I14`) that points into `experiments/*/impl/ISSUES.md`.
- Design guidance (not enforced by a tool): modules with a small interface and the behaviour behind it. Add a seam only where two adapters exist, as the checker has (CPU build and nvcc).
- Never raise from the parser. `parse` returns a code.
- Do not copy masks to the host for scoring. See [ADR-0004](architecture/decisions/0004-scoring-on-the-gpu.md).
- In the checker, write algorithms as functors over independent items. Keep integer arithmetic. Keep the no-FMA build flags. See [CUDA engine](reference/cuda-engine.md).
- Do not change a `REWARD_V4` value. It is a new reward, needs a new lab round, and invalidates stored expert tables.
- Do not change a rule in `experiments/phase2/prereg.py` after you have seen the result it decides.

## 5. Data and secrets

- Downloaded files and extracted archives are untrusted. Put each in its own new, empty directory. Keep scripts elsewhere.
  Run every script that reads them with `python -I`, and pass paths as arguments. A planted `json.py` would otherwise run on `import json`.
- Never commit data. No image or annotation of LIDC-IDRI or RIGA belongs in the repository. RIGA-derived tables stay on the owner's machine.
- Secrets are read by `objdel.secrets` from `<repo>/.secrets` (mode 600, git-ignored) or the macOS keychain. Never print a value. Never write one to a log, a ledger, a config or a command line.
- RIGA goes only to Secure Cloud hosts, with recorded written permission. See [Data](reference/data.md).

## 6. Pull request checklist

- [ ] `make test` and `make lint` pass.
- [ ] New library code has tests. Exact code compares with an oracle or a golden, not a tolerance.
- [ ] The matching document in `docs/` is updated.
- [ ] No data, secret, or large result file is added.
- [ ] A GPU-side change names the box run that tested it, or says it is untested on a GPU.
- [ ] A change to the reward, the grammar, the frame policy or the checker states its effect on stored hashes and expert tables.
- [ ] The commit message has no AI attribution. The description has no "Generated with" line.

## 7. Writing documentation

- Use simple technical English. Sentences have at most 20 words, with active voice and present tense. One idea per sentence.
- Use one word for each thing. Define a new term once in the [glossary](glossary.md) and use only that word afterwards.
- Cite the source of every number: a code path, a config, a test or a results file. Do not write a fact you have not checked.
- Use ASCII boxes and arrows for diagrams. Give a worked example with real numbers.
- Record a decision as an ADR only if it is hard to reverse, surprising without context, and a real trade-off ([decision records](architecture/decisions/README.md)).
