# Phase 1 lab, archived 2026-11-09

**Nothing here is built, imported or tested.** The evidence behind `experiments/phase1/PLAN.md`. The RIGA-derived results of this folder stay on the owner's disk only and are git-ignored (RIGA terms; `.gitignore` RIGA block).

## What it was

Sandboxed prototypes and tests that chose each Phase 1 module. Data, fixtures, model caches and large outputs were not copied (licences, size). Online research notes: `docs/research/notes/Phase 1 stack and algorithms/`.

| Folder | What it tested | Report |
|---|---|---|
| `common/` | exact pure-Python oracle (raster spec, legality labels) + its tests | `common/test_oracle.py` |
| `t_polycore/` | C++ polygon core: rasterizer + legality, nanobind, Linux Docker | `REPORT.md` |
| `t_reward/` | reward core on real expert masks: topology, band F5, K9, advantages, adversaries, neutrality | `REPORT.md`, `results/ITERATION2.md`–`ITERATION4.md` |
| `t_targets/` | SFT target builder (mask → canonical polygon JSON) | `REPORT.md` |
| `t_stack/` | TRL SFT/GRPO on 4 small VLMs; integration R1/R2/K12; resolution and licence | `REPORT.md`, `results/INTEGRATION.md` |
| `t_data/` | RIGA staged extractor, near-duplicates, LIDC parity vs pylidc, splits, storage | `REPORT.md` |
| `t_dynamics/` | GRPO toy: set-valued topology collapse and fixes; Level 2 controller simulation | `REPORT.md`, `results/*_tables.md` |
| `tools/` | fixture builder | — |
| `DECISIONS.md`, `KPIS_DRAFT.md` | decision record and draft KPIs reviewed by the judge | — |
| `JUDGE_VERDICT.md` | independent judge review (required changes R1–R10) | — |

`t_stack/prod.py` puts `vendor/` on `sys.path`. 34 vendor files were deleted as byte-identical twins (path, md5 and twin path of each: `archive/phase1_lab/vendor_deleted.tsv`), so run it from the tag, not from `archive/`.

`t_targets/record_target_goldens.py` recorded `tests/goldens/targets_lab.npz`; it imports the lab, so it also runs only from the tag.

## What replaced it

`src/objdel` with its tests in `tests/` (promoted during Phase 1 milestones M1 to M8), and the goldens `tests/goldens/polycore.npz` and `tests/goldens/targets_lab.npz`.

## Retiring commit

Code last present in commit `207685bf1b57` (tag `pre-clean-2026-11-08`). The commit that adds this README records the move.

## How to restore

`git restore --source=pre-clean-2026-11-08 -- experiments/phase1/lab` (owner).
