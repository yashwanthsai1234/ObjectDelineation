# Phase 0: frozen record

Phase 0 is done. Its results are committed in this folder ([RESULTS.md](RESULTS.md), [RESULTS_0b.md](RESULTS_0b.md),
`results/`). The scripts are kept unchanged in `archive/phase0/` as the record of what produced those numbers.

The RIGA-derived results of this folder stay on the owner's disk only and are git-ignored (RIGA terms; `.gitignore` RIGA block).

## Run the scripts against commit `faaacc4`

Phase 1 replaced the library code that five of the scripts import. Running them against the current tree fails at
import:

| Script | Imports removed in Phase 1 | Phase 1 replacement |
|---|---|---|
| `lidc_census.py` | `objdel.checker.betti3d` | none yet (3D is Phase 7) |
| `riga_census.py` | `objdel.reward.AnnotatorSet` | `objdel.reward` (reward v4); the Phase 0 scheme survives as `phase0` in `tests/oracles/reward_variants.py` (A/B test AB2) |
| `riga_outlines.py` | `objdel.checker.topology.FG8`/`BG4`, `AnnotatorSet` | `objdel.checker.topology` (OpenCV labeling) |
| `baseline_vlm.py` | `objdel.checker.polygon_io`, `AnnotatorSet` (and `betti3d` through `lidc_census`) | `objdel.geometry` (parser, `to_mask`) |
| `baseline_medsam.py` | `AnnotatorSet` (and `betti3d` through `lidc_census`) | `objdel.reward` |

`lidc_fetch.py` and `summarize_0b.py` do not import the library and still run on the current tree
(Phase 1's `fetch_lidc_pixels.py` reuses `lidc_fetch.sop_to_z`).

Commit `faaacc4` is the first commit with the `src/` + `experiments/` layout. Its `src/` and `experiments/phase0/`
are identical to those of `ab1b0b9`, the last commit before Phase 1 code (`git diff faaacc4 ab1b0b9 --
src experiments/phase0` is empty). To rerun a script, check that commit out in a separate worktree and point it at
the same (git-ignored) data:

```bash
git worktree add ../od-phase0 faaacc4
cd ../od-phase0
make setup                       # .venv with the Phase 0 extras (models, data, dev)
ln -s ../ObjectDelineation/data data
.venv/bin/python -I -m experiments.phase0.<script> ...   # the command at the top of each script
git worktree remove ../od-phase0 # when done
```

Checked on 2026-10-19 (a `git archive faaacc4` copy, `[dev]` extra): its test suite passes (20 tests); the five
scripts above import (the two baselines need the `models` extra); `summarize_0b` regenerates
`results/summary_0b.json` byte for byte.

Sources: `git worktree` (https://git-scm.com/docs/git-worktree).
