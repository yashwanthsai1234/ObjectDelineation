# Phase 0 scripts, archived 2026-11-09

**Nothing here is built, imported or tested.** The Phase 0 census and baseline scripts. Their results and notes are in `docs/research/records/phase0/`.

## What it was

Seven scripts and an `__init__.py` that produced the Phase 0 numbers: `lidc_fetch.py`, `lidc_census.py`, `riga_census.py`, `riga_outlines.py`, `baseline_vlm.py`, `baseline_medsam.py`, `summarize_0b.py`. Five of them import library code that Phase 1 removed. The record README (`docs/research/records/phase0/README.md`) explains how to rerun them from commit `faaacc4`.

## What replaced it

The records moved to `docs/research/records/phase0/` (`RESULTS.md`, `RESULTS_0b.md`, `PLAN*.md`, `results/`). The library modules replaced them in `src/objdel` (Phase 1). `archive/phase0/baseline_medsam.py` is still cited for `LUNG_WINDOW` in `src/objdel/targets/frame.py`.

## Retiring commit

Code last present in commit `207685bf1b57` (tag `pre-clean-2026-11-08`). The commit that adds this README records the move.

## How to restore

`git restore --source=pre-clean-2026-11-08 -- experiments/phase0` (owner).
