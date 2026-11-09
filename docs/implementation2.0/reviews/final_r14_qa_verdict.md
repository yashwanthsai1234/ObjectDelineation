# Verdict: final round 14 QA review (final_r14_qa.md)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product choice in scope; no catalog skill recommended. No change to verdicts. |
| engineering:code-review | Read the cited code: .gitignore:37-49, tests/conftest.py:44-56, tests/test_p2_breakit_ops.py:297-300, policy.py code_files() (drop regex), tests/goldens/ contents. |
| engineering:debug | Traced where a D-R01 golden would go: tests/goldens/ -> git (public repo) and code_files() -> box code sync; checked which guard would stop it. |
| engineering:system-design | Checked the one-owner rule for RIGA: D.3.2 and D.4 say RIGA never leaves the Mac; D-R01 contradicts that. |
| engineering:architecture | Licence boundary: LIDC goldens allowed, RIGA per-case data not (owner rule, .gitignore:37; CLEAN.2 RIGA terms). |
| engineering:documentation | Checked whether any other doc part already handles D-R01 (no review-log entry covers it). |
| huggingface-skills:huggingface-best | No model choice in scope. No change. |

Web tools loaded (WebSearch, WebFetch); no web check was needed: every claim is decidable from the repo and the sources.

## Items

- qa-1: CONFIRMED (blocking). D.md D.9.4 row D-R01 says "Always-run copy: 50 fixed LIDC and RIGA items in `tests/goldens/` equal the stored pre-D text"; an always-run test has no `data/`, so the golden must hold the RIGA rater inputs and/or per-case RIGA target texts, which breaks the owner rule "RIGA-derived results stay local: no permission to share RIGA data (repo public, owner 2026-11-06)" [code: .gitignore:37] and contradicts D.3.2 "RIGA data never leaves the Mac"; no planned check catches it, because CLEAN-G16/CT-11 read only `tests/goldens/*.npz`, `test_code_snapshot_has_no_riga_derived_results` and `code_files()` match only file names containing "riga" (`riga[^/]*\.(csv|jsonl?)$` in policy.py code_files; tests/test_p2_breakit_ops.py:298), and `code_files()` ships `tests/` to the box; a name like `tests/goldens/d_r01.json` passes every check and is published at the D3 commit (irreversible). Fix touches D.md only (D-R01 LIDC-only golden plus a "no riga id" assert; D-R02 wording "synthetic disc/cup/rim shapes").
- qa-2: not re-judged (already minor in the review; caught at $0 by file-not-found).
- qa-3: not re-judged (already minor; `make test-order` goes red at $0).
- qa-4: not re-judged (already minor).
- qa-5: not re-judged (already minor; both failure cases are visible at $0).
- Carry-overs (round 13 qa-3..qa-9): minor, not re-argued.

Totals: confirmed 1, refuted 0, downgraded 0.
