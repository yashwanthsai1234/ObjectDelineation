# Verdict: final_r13_fresh-arch

## Skills used
- nvidia-skills:nvidia-skill-finder: loaded; no item depends on an NVIDIA product, so no catalog lookup.
- engineering:code-review: read `.gitignore` (RIGA block), `tests/test_p2_breakit_ops.py:297-300`, CLEAN5 steps 1-9, CT-87, the CLEAN.9 RIGA greps.
- engineering:debug: reproduced item 1 read-only (`git status --porcelain --ignored`, `git check-ignore -v --no-index` on the post-move paths: exit 1).
- engineering:system-design: traced item 2 through time line step 6, 6.3 V3 item 5, the (b') rule and decision 1.
- engineering:architecture: item 2 is a decision-input gap under a clear "no rental after V1" rule, not a contradiction.
- engineering:documentation: every confirmed item quotes the source heading and the code line or measurement.
- huggingface-skills:huggingface-best: no model choice in question.

## Verdicts
- fresh-arch-1: CONFIRMED (blocking). CLEAN.md CLEAN5 "Files to move (plain `mv`)" moves `experiments/phase0/results/` to `docs/research/records/phase0/` and `experiments/phase1/lab/` to `archive/phase1_lab/`, and step 8 has the owner run `git add -A experiments archive docs/research/records ...`; `.gitignore` ignores the RIGA per-case files only at their old paths ("RIGA-derived results stay local ... owner 2026-11-06": `/experiments/phase0/results/medsam_riga.csv`, `/riga/`, `/vlm_box/`, ... `/experiments/phase1/lab/t_reward/results/h5_items.csv` ...); measured 2026-11-09: these files exist as `!!` ignored (medsam_riga.csv holds per-image RIGA IoU rows `BinRushed_BinRushed2,1,disc,6,0.9485,...`), and `git check-ignore -v --no-index docs/research/records/phase0/results/medsam_riga.csv archive/phase1_lab/t_reward/results/h5_items.csv` returns exit 1, so the step 8 `git add -A` stages them; no planned check catches it (`test_code_snapshot_has_no_riga_derived_results` reads `P.code_files()`, which drops `docs/` and `archive/` [CLEAN1 regex, CLEAN.md:368]; CT-87 only counts files; the CLEAN.9 grep covers png/tif/npy/npz only; DOC-G1 excludes `docs/research`); no review log marks it handled.
- fresh-arch-2: DOWNGRADED (minor). The plan is consistent, so no guess changes the outcome: 00_front 6.3 places the candidate's test eval in V3 (item 5), 6.3 and 99_back decision 1 (b') say "No rental follows V1: V2 ..., V3 ... and V4 each wait for a budget raise", and time line step 6 says "Under the default (b') neither runs"; the engineer therefore produces no new test number and spends nothing; the gap is only that decision 1 does not spell out "no test number for B1 under (b')", a one-line wording fix (reviewer fix 5), while the proposed V1 `--answers-only` change would reopen the V1 bundle rule (6.1 step 4c "no `lidc-test-*` folder"), the 4d grep test and the approved V1 amount.
- fresh-arch-3: minor as filed (no change).
- fresh-arch-4: minor as filed (no change).
- fresh-arch-5: minor as filed (no change).

Counts: confirmed 1, refuted 0, downgraded 1.
