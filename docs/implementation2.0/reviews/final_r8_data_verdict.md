# Verdict: final_r8_data.md (fresh verifier, 2026-11-08)

## Skills used
- nvidia-skills:nvidia-skill-finder: no NVIDIA product in scope (data statistics); no catalog lookup; nothing changed.
- engineering:code-review: read r8data/legacy.py and legacy_found.py against the A2 definition (A.md:365-367) before trusting their numbers.
- engineering:debug: reproduced item 1 myself (both scripts on data/store and runs/r5/eval) and the item 2 hash-seed claim.
- engineering:system-design: traced the 0.734 reference to Research2.0.md:196-197 (one number used for "All images" and "found").
- engineering:architecture: checked that A2 step 4, A.6 A-G2, A.9.6 G04 and 99_back 12 all own the same wrong reference.
- engineering:documentation: one line per item, each with its proof.
- huggingface-skills:huggingface-best: no model choice in scope; nothing changed.

## Verdicts
- data-1: CONFIRMED (blocking). My own run of `.venv/bin/python -I r8data/legacy.py data/store` gives val item 0.7232 / pooled 0.7419 (voters the same) and test 0.7319 / 0.7486, and legacy_found.py gives val found(947) 0.7339; so the A2 step 4 / A-G2 (A.md:420, 1122) / G04 (A.md:1906) / 99_back 12 reference "0.734 (val)" is the found-image value (Research2.0.md:196 prints 0.734 for "All images" too), and a correct A2 fails both legacy lines by 0.0108 and 0.0079 > 0.005, so A2 stops and `make test-data` cannot pass. Fix touches A, 99_back, B (B.md:341 ADR-0015 "Confound"), DOC (Research2.0.md:196); A.spec.md:315 carries the same text and should follow.
- data-2: minor (as filed). Reproduced: `env -i PYTHONHASHSEED=0 python3 -I` prints different hashes each run, but the test still passes on correct code; no outcome change.
- data-3: minor (as filed). Wording of the rounding rule; not escalated.
- data-4: minor (as filed). DG12 already reads the right root; wording only.
- data-5: minor (as filed). Carried-over list from rounds 6-7; not re-argued.

Counts: confirmed 1, refuted 0, downgraded 0.
