# Implementation plan (phases A to E, DOC, CLEAN)

Status: plan only (2026-11-09). No code in `src/`, `experiments/` or `cuda/` was changed to write it. No GPU was rented.

| File or folder | What it holds |
|---|---|
| [Implementation.md](Implementation.md) | The plan: front matter (words, order of work, target architecture, model, tech stack, kernels, the only GPU rental plan, risks), then phases A, B, C, D, E, DOC, CLEAN, each with tasks, gates, costs and tests; then owner decisions and sources |
| [Architecture-explained.md](Architecture-explained.md) | The same architecture in ASD-STE100 plain English, with diagrams and worked examples. Read this first |
| [review-log.md](review-log.md) | What each review round found and what changed (phase reviews, then final rounds 1-15) |
| [kernels/](kernels/) | CUDA kernel drafts K1 (polygonizer), K2 (scoring), K3 (RL-loop helpers): `DESIGN.md`, `GPU_RUN.md`, source, CPU tests. Not built into the package |
| [specs/](specs/) | The short spec of each phase that the zero-context test writers read |
| [tests/](tests/) | The two independent test plans per phase, as written (merged into each phase's section 9) |
| [reviews/](reviews/) | Phase reviews, final review rounds and their verdict files |
| [notes/](notes/) | Research notes, measurements and probe scripts that the plan cites |
| [seeds/](seeds/) | Seed scripts that phase B promotes into the code |

## Review status

Fifteen final review rounds ran (six expert reviewers each, rounds 4-15 with a refute-first verifier per review).
Confirmed blocking plus major items per round: 29, 22, 20, 14, 14, 14, 11, 8, 5, 4, 4, 4 (rounds 4-15); each round's
items were fixed. The owner stopped the loop after round 15; its items were fixed but not re-reviewed, and the round-15
QA and RL reviews have no verdict file. Expect detail-level gaps (exact commands, config fields) to surface
during implementation; fix them in the phase section and the review log.

## Not kept

Files the plan cites that are not in the repo (over 1 MB):

- not kept: notes/b1_out/geom.json, 9,108,536 B, over 1 MB; rebuild with seeds/b1_scripts
- not kept: notes/b1_out/pix.json, 3,622,228 B, over 1 MB; rebuild with seeds/b1_scripts
