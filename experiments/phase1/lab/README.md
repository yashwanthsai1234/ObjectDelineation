# Phase 1 lab: the evidence behind `../PLAN.md`

Sandboxed prototypes and tests that chose each Phase 1 module. They are frozen; production code is promoted into `src/objdel` with
its tests during Phase 1 milestones M1–M8. Data, fixtures, model caches and large outputs were not copied (licences, size).

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

Online research notes: `docs/research/notes/Phase 1 stack and algorithms/`.
