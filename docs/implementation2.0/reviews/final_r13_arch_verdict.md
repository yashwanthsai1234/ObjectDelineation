# Verdict on final_r13_arch.md (fresh verifier, 2026-11-09)

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | No NVIDIA product claim in the items; no catalog lookup, no install. |
| engineering:code-review | Read CLEAN5 steps 4, 8, 9, CT-73, CT-87 and the CLEAN12/CLEAN4 guards against the real `.gitignore` and `code_files()`. |
| engineering:debug (listed twice in the task; run once) | Reproduced item 1: listed the ignored files in the moving trees and ran `git check-ignore` at old and new paths. |
| engineering:system-design | Walked the order DOC1 -> CLEAN5 -> commit and asked which planned check sees an ignored file after it moves. |
| engineering:architecture | Seam check: the only RIGA leak test reads `code_files()`, which drops `docs/` (and `archive/` after CLEAN1). |
| engineering:documentation | Checked whether any other phase (DOC, E, front, back) re-handles the phase 0 results move. None does. |
| huggingface-skills:huggingface-best | No model choice in scope; nothing changed. |
| WebFetch | git-scm.com/docs/gitignore: "If there is a separator at the beginning or middle ... the pattern is relative to the directory level of the particular .gitignore file itself." |

## Items

| id | verdict | reason |
|---|---|---|
| arch-1 | confirmed (blocking) | Verified: `.gitignore:37-49` anchors 12 RIGA patterns to `/experiments/phase0/results/...` and `/experiments/phase1/lab/...`; `git ls-files -o -i --exclude-standard experiments/phase0 experiments/phase1/lab` lists exactly 31 RIGA-derived files (minus `__pycache__`, `*.log`); `git check-ignore -q` exits 0 at `experiments/phase0/results/medsam_riga.csv` and 1 at `docs/research/records/phase0/results/medsam_riga.csv` and `archive/phase1_lab/t_data/riga/maxrg_vs_prime.json`; CLEAN.md "### CLEAN5" moves the folders with plain `mv` and step 8 runs `git add -A experiments archive docs/research/records ...` with no `.gitignore` edit; no planned check catches it (`test_code_snapshot_has_no_riga_derived_results` reads `code_files()`, whose drop regex `^(data|docs)/` [policy.py:379] hides `docs/`, and CLEAN1 adds `archive/`; CT-73 greps only `png|tif|npy|npz`; CT-87 counts files; DOC1's RIGA owner check runs before CLEAN5 and DOC1 step 3 excludes `docs/research`); many moved names (`vlm_box/scores.csv`, `t_stack/results/t2_*.json`, `riga/images.csv`) have no `riga` in the basename, so even CLEAN1's `riga[^/]*\.(csv|jsonl?)$` misses them; result: RIGA-derived tables pushed to the public repo, against the RIGA terms CLEAN.2 quotes and the owner's `.gitignore` rule; not in any review log. Note for the fixer: the proposed CLEAN12 test (path contains `riga` or `records/phase0/results/vlm_`) misses `archive/phase1_lab/t_stack/results/t[24]_*.json`; the `.gitignore` rewrite (or the move into `data/`) plus the step 9 mapped-path `check-ignore` gate is what closes the hole. |
| arch-2 | minor (as filed) | Stale "moved by D3" wording for the numpy oracle; `old_names.json` is recorded on E0 after both moves, so no outcome changes. |
| arch-3 | minor (as filed) | U10 and the E2.5 AST test both fail on the Mac at $0 and point at the gap (caught-before-spend rule). |

Counts: confirmed 1 (blocking), refuted 0, downgraded 0 (items 2-3 were already minor).
