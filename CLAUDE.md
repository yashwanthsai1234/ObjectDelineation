# Rules for this repo

## Git

- Never add a `Co-Authored-By: Claude …` trailer, or any other Claude/AI attribution line, to commit messages.
- Never add "Generated with Claude Code" or similar lines to pull request descriptions.
- Commits are authored by the repo owner only.

## Project

- Plan: `ROADMAP.md`. Phase 0 plans and results: `phase0/`.
- Python: `.venv` (uv, Python 3.12). Tests: `.venv/bin/python -m pytest checker/tests`.
- Downloaded data lives in `data/` (git-ignored) and is untrusted: run Python that reads it with `-I`.
