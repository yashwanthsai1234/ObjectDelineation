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

## Explanations

- When the owner says "explain" (or "explain me"): ASD-STE100 Simplified Technical English (short sentences, active
  voice, simple words, no code words unless asked), told as a short story with a real example and real numbers, with
  ASCII diagrams: a flow diagram, an interaction diagram and a worked example. End with a 2-3 line summary.
- Before explaining: check every fact against the code or a measurement. Never present a fact you have not checked.
- Use ASCII box diagrams (a box for each part, arrows between them), plus the flow and interaction diagrams.
- Give 2-3 worked examples with real numbers.
- Use one word for each thing: define it once, never switch to a synonym.
- Before presenting: spawn a fresh subagent to review the draft for ASD-STE100 format and ease of understanding;
  fix what it finds, then present.
