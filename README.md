# ObjectDelineation

An exact, self-certifying geometry checker for masks and polygons. It is used as a **verifiable reward**
for training vision-language models to segment, through RL with GRPO.

The checker counts parts and holes exactly (bit-quads, Betti numbers) and turns masks into polygons.
Every result comes with a certificate. Model outputs are scored against the **set** of answers that expert
annotators gave, not against a single ground truth. A shape that any expert drew counts as correct, and pixels
where experts disagree count as neither right nor wrong.

## Status

| | |
|---|---|
| Phase 0 (census and baselines) | done: [results](experiments/phase0/RESULTS.md), [model baselines](experiments/phase0/RESULTS_0b.md) |
| Next | Phase 1: finish the reward terms, prepare data |
| Full plan | [ROADMAP.md](ROADMAP.md) |

## Layout

```
src/objdel/        installable library
  checker/         exact geometry: bit-quads, Betti numbers, polygons, certificates (numpy/scipy only)
  reward/          scoring against an annotator set; GRPO group advantages
  models/          model wrappers (MedSAM, …); needs the `models` extra
tests/             pytest suite for the library
experiments/       one folder per phase: plan, scripts, results; frozen once the phase is done
cuda/              the original C++/CUDA engine (CMake, GoogleTest)
docs/research/     research reports and notes behind the direction
data/              downloaded datasets (git-ignored, never committed)
```

**Rules**
- Reusable code goes in `src/objdel`, with tests.
- Experiment scripts only wire library code to data.
- Results are committed next to the script that produced them.
- Data is never committed.

## Quick start

```bash
make setup        # .venv with Python 3.12 (uv), package installed editable with all extras
make test         # Python tests
make lint
make cuda         # build the C++/CUDA engine
```

```python
import numpy as np
from objdel.checker import certify
from objdel.reward import AnnotatorSet

report = certify(mask)                      # parts, holes, area, perimeter, rings and a certificate
score = AnnotatorSet(expert_masks).score(prediction)   # IoU, band, topology against all experts
```

To reproduce Phase 0, see the commands at the top of each script in `experiments/phase0/`.
The datasets are LIDC-IDRI (TCIA) and RIGA. Each is downloaded separately under its own licence.
