# Task graph and dependency graph (implementation of Implementation.md)

The plan is `Implementation.md`. This file only says **who does what, in which order, on which files**. If the
plan and this file differ on content, the plan wins. If they differ on order or ownership, this file wins.

## Execution decisions (owner, 2026-11-09)

1. **Everything is implemented**, including the work the plan parks on the `after-raise` branch under (b')
   (C0-C11, B2.2-B2.4, B3.1-B3.10). It lands on `impl2` on pre-E paths, before E, as the plan's "raise" case
   ("if E0 has not started, the branch merges on pre-E paths"). Gates that need a paid run (C-G-SP1, Gate B2,
   Gate B3) are not read. So the plan's temporary flags stay, and **their defaults keep the B1 control**: reward v4,
   `bf16-seqmask`, CPU pixel path, merged weight push. The owner turns them on in a run spec. C10 and B3.10 do only
   their docs and ADR parts. Every flag is listed in `docs/implementation2.0/IMPL_LOG.md`.
2. Owner decisions take the bracketed defaults of "Owner decisions needed" in the plan.
3. GPU: one test rental of at most $3 (Vast, US/CA), after E. It tests GPU scripts and runs a short RL smoke.
   It is not the B1 run. The owner starts the RL run.
4. Git: the main agent and lane agents may commit. Author = repo owner. **No `Co-Authored-By`, no "Generated with"
   line, no AI attribution.** Never push, never force, never rewrite history. Only the main agent merges.
5. No HF push and no Vast call by lane agents.

## Line ranges in Implementation.md

Front matter 1-1106 (0 words 20-72; 1 order 73-261; 2 architecture 262-437; 6 rentals 604-1060).
Phase A 1107-3621 (tasks 1347-2272, gates 2284, tests 2367). Phase B 3622-6823 (design 3849, tasks 4028-5535,
gates 5552, tests 5854). Phase C 6824-9150 (tasks 7125-7949, gates 7967, tests 8171). Phase D 9151-10749 (tasks
9348-10256, gates 10269, tests 10361). Phase E 10750-13248 (design 10939, tasks 11060-12525, gates 12556, tests
12696). DOC 13249-14652 (tasks 13431-14156, tests 14256). CLEAN 14653-16851 (tasks 14913-15673, tests 15792).
Owner decisions 16879-17082.

Tasks: A1 1363, A2 1476, A3 1573, A4 1792, A5 1844, A6 1886, A7 2112, A8 2206, A9 2247 | B1.0 4064, B1.1 4136,
B1.2 4195, B1.3 4248, B1.4a 4348, B1.4 4394, B1.5 4445, B1.6 4598, B1.7 4659, B1.8 4732, B1.9 4748, B2.0 4804,
B2.1 4817, B2.2 4896, B2.3 5104, B2.3b 5125, B2.4 5159, B3.1 5262, B3.2 5280, B3.3 5299, B3.4 5312, B3.5 5350,
B3.5b 5370, B3.6 5387, B3.7 5399, B3.9 5466, B3.10 5521 | C0 7142, C1 7282, C2 7406, C3 7483, C4 7548, C5 7606,
C6 7734, C6b 7763, C7 7773, C8 7818, C9 7840, C10 7857, C11 7889 | D1 9414, D2 9562, D3 9664, D4 9722, D4b 9760,
D5 9784, D6 9875, D7 9900, D7b 9994, D8 10038, D9 10231 | E0 11137, E1 11334, E3 11378, E5 11455, E2 11569,
E4 11756, E6 11889, E8 12246, E7 12480 | DOC1 13458 ... DOC13 14136 | CLEAN0 14968 ... CLEAN16 15659.
Short specs: `specs/<phase>.spec.md`. Test plans: `tests/`. Kernel drafts: `kernels/K1..K3`.
Research notes for builders: `notes/research/*.md` (written in wave R).

## Dependency graph (waves)

```
 W0  DOC1 DOC2 CLEAN0  (impl2)          R  4 web research notes (parallel with W0)
        |                                   |
        +-----------------+-----------------+
                          v
 W1  five lanes, each in its own worktree and branch (parallel)
     L1 eval   A1 > A2 > A4 > A5 > A3(code) > A8 > B1.9
     L2 data   B1.0 > B1.7 > B1.1 > B1.2 > B1.3 > B1.6
     L3 cuda   B1.4a > D1 > D2 > D3 > D4 > D4b > D5 > D6 > D7b > B3.1
     L4 train  A7 > A6 > B2.0 > B2.1 > B2.3b(Mac) > B1.5 > CLEAN2 > CLEAN6 > CLEAN7
     L5 clean  CLEAN1 > CLEAN3 > CLEAN4 > CLEAN5 > CLEAN9 > CLEAN10 > CLEAN12 > CLEAN13 > CLEAN14 > CLEAN15
                          |  merge M1 (main agent), full suite green
                          v
 W1b L2b rows  B1.4 (build rows root on Mac data) > all 21 B data gates > B1.8 rows path > b1 dry run + preflight
     L3b cuda  D7 > D8 > D9
                          |  merge M2, full suite green
                          v
 W2  L6 speed  C0 > C1 > C2 > C3 > C4 > C5 > C6 > C7 > C8 > C9 > C11 > C10(docs/ADR only)
     L7 reward B3.2 > B3.3 > B3.4 > B3.5 > B3.5b > B3.6 > B3.7 (lab, $0) > B3.9 (configs) > B3.10 (docs only)
     L8 numer. B2.2 > B2.3 (configs) > B2.4 (decision code)
                          |  merge M3, full suite green
                          v
 W3  E (serial, on impl2): E0 > E1 > E3 > E5 > E2 > CLEAN8 > E4 > E6 > E8 > E7
                          v
 W4  break-it round 1 (4 agents) > fix fleet > GPU test rental (<= $3) > break-it round 2 > fix fleet
                          v
 W5  DOC3-DOC13 (+A9 content) > CLEAN11 > CLEAN16 final sweep > full suite + lint + docs build
```

Why this order: B1.0 changes the cache fingerprint (one rebuild, first). B1.7 runs before any eval view is hashed.
B1.4 needs B1.4a's polygonizer seam (L3) and B1.1-B1.3 (L2). D7 edits `final_eval.py` and `fast.py`, so it waits for
L1 and L4. C, B2.2+ and B3.2+ edit `fast.py`, `configs.py`, `rewards.py`, so they wait for W1. E moves files, so it
waits for every code lane. CLEAN8 needs E2. DOC final sweep and CLEAN16 come last.

## File ownership in W1 (a lane edits only its own files; shared files get the smallest possible hunk)

| Lane | Owns |
|---|---|
| L1 | `experiments/phase2/final_eval.py`, `select_ckpt.py`, `ceiling.py`, `headline.py`, `outline_lengths.py`, `gate_power.py`, `build_eval_views.py` (A3 part), `src/objdel/eval/{metrics,signal,select,ceiling}.py`, `src/objdel/train/sample.py` |
| L2 | `src/objdel/data/{windows,coverage,quality}.py`, `src/objdel/reward/spec.py` (B1.0), `src/objdel/eval/{coverage,gates,audit}.py` (B1.3), `experiments/phase2/{audit_data,check_series_z}.py`, `docs/reference/{frozen_views,val_halves,train_exclusions,neardup_review}.json`, `train/dataset.py` (B1.6) |
| L3 | `cuda/*`, `src/objdel/checker/*`, `src/objdel/targets/build.py`, `src/objdel/reward/cache.py`, `experiments/cuda_checker/*`, `experiments/phase2/perf/bench_checker.py`, `verify_target_hashes.py`, oracles |
| L4 | `src/objdel/train/{fast,callbacks,rewards,configs,trl_adapter,nvtx}.py`, `experiments/phase2/{train_rl,train_sft,g_probe,profile,preflight,merge_rounding,b2_attribution,nsys_baseline}.py` |
| L5 | `experiments/phase2/vast/*`, `archive/*`, `experiments/phase0/*`, `experiments/phase1/lab/*`, `tests/test_repo_hygiene.py`, `tests/vulture_whitelist.py`, ADR 0021 |

`preflight.py` and `vast/policy.py` are shared: L4 owns `preflight.py`, L5 owns `vast/policy.py`; other lanes add
one small hunk each and say so in their report. The merge resolves the rest.

## Rules for every lane agent

- Work only in your worktree. Run Python from its `.venv`. Data reads use `python -I`.
- TDD per task: write the test from the plan's section 9 sketch, run it, see it fail for the right reason, write the
  simplest code that passes, refactor. A proposal is tested (scratch script, Docker, or a test) before it lands.
- Karpathy rules: least code that does the job; no speculative options; no dead or stale code; replace, do not
  keep old beside new; read the code you change first.
- During work run only the tests you touch. At lane end run `.venv/bin/python -m pytest -q -m "not slow" -x` once
  and `.venv/bin/python -m ruff check src tests experiments`. Both must pass.
- Commit after each task: `git commit -m "<task id>: <what>"`. No AI attribution of any kind.
- Plan gaps or plan-vs-code conflicts: choose the simplest fix that keeps the plan's gate, and write one line in
  `docs/implementation2.0/IMPL_LOG.md` under your lane (`- <task>: <deviation> (why)`).
- No GPU rental, no Vast CLI, no HF upload, no `git push`.
