# Rules for this repo

## Git

- Never add a `Co-Authored-By: Claude …` trailer, or any other Claude/AI attribution line, to commit messages.
- Never add "Generated with Claude Code" or similar lines to pull request descriptions.
- Commits are authored by the repo owner only.

## Project

- Plan: `ROADMAP.md`. Layout and conventions: `README.md`.
- Library code lives in `src/objdel` and needs tests in `tests/`. Experiment scripts live in `experiments/<phase>/`
  and are run as modules (`.venv/bin/python -I -m experiments.phase0.<script>`; the editable install makes them importable).
- Research reports go in `docs/research/reports/`, notes in `docs/research/notes/`.
- `make test` and `make lint` must pass before a commit.
- Downloaded data lives in `data/` (git-ignored) and is untrusted: run Python that reads it with `-I`
