# Phase 2 lab, archived 2026-11-09

**Nothing here is built, imported or tested.** The sandboxed prototypes and reports (folders `A1`..`A6`, `B1`..`B5`) that chose the phase 2 design. Some scripts import modules that no longer exist (`polycore`, `objdel.checker.bitquad`).

## What it was

Prototypes, reports and result files of the phase 2 lab. Production code was promoted from it into `src/objdel` and `experiments/phase2`.

Three tests compared live code against lab functions. They now compare against `tests/goldens/p2_lab_refs.json`, recorded by `record_p2_lab_refs.py` in this folder (it runs only from a checkout of the tag). Fixtures that tests read moved to `tests/fixtures/` (`vast_snapshot.json`, `topo_probe_v3e_seed*.json`).

## What replaced it

`src/objdel/eval/{signal,probes,vision}.py`, `experiments/phase2/*.py`, `tests/goldens/p2_lab_refs.json`, `tests/fixtures/`.

## Retiring commit

Code last present in commit `207685bf1b57` (tag `pre-clean-2026-11-08`). The commit that adds this README records the move.

## How to restore

`git restore --source=pre-clean-2026-11-08 -- experiments/phase2/lab` (owner).
