# ObjectDelineation

An exact, self-certifying geometry checker for masks and polygons. It is used as a **verifiable reward**
for training vision-language models to segment, through RL with GRPO.

The checker counts parts and holes exactly and validates and rasterizes polygons with integer arithmetic (C++).
Model outputs are scored against the **set** of answers that expert annotators gave, not against a single ground
truth. A shape that any expert drew counts as correct, and pixels where experts disagree count as neither right nor
wrong.

## Status

| | |
|---|---|
| Phase 0 (census and baselines) | done: [results](experiments/phase0/RESULTS.md), [model baselines](experiments/phase0/RESULTS_0b.md); scripts run against commit `faaacc4` ([note](experiments/phase0/README.md)) |
| Phase 1 (reward and data) | done: [results and gates](experiments/phase1/RESULTS.md), [plan](experiments/phase1/PLAN.md), [execution](experiments/phase1/EXECUTION.md), [notes and issues](experiments/phase1/impl/) |
| Phase 2 (SFT on Vast AI) | in progress: rung 1 (LFM2.5-VL-450M) failed the entry gates; next rung Qwen3-VL-2B (LoRA): [results](experiments/phase2/RESULTS.md), [plan](experiments/phase2/PLAN.md), [build log](experiments/phase2/V0_V1.md), [issues](experiments/phase2/impl/ISSUES.md) |
| Full plan | [ROADMAP.md](ROADMAP.md) |

## Layout

```
cpp/polycore/      header-only C++17 core: polygon rasterizer, legality checks, signed area (+ native fuzz tests)
src/objdel/        installable library
  _polycore        nanobind extension built from cpp/polycore (scikit-build-core)
  checker/         bit-quads, Betti numbers, instance topology + size rule, polygons, certificates, polycore wrapper
  geometry/        the model output spec: strict parser, serializer, canonical rings, polygon -> mask
  targets/         model frame (crop/scale) and certified SFT polygon targets
  reward/          reward v4: per-item expert caches, per-rollout terms, group advantages
  data/            RIGA extractor, LIDC parser, group-aware splits, Parquet item store
  train/           TRL integration: chat rows, reward functions, controller callback, configs, sampler
  eval/            metrics harness, hack audit, Phase 2 entry gates
tests/             pytest suite (oracles in tests/oracles/)
experiments/       one folder per phase: plan, scripts, results; frozen once the phase is done
cuda/              the C++/CUDA checker engine (CMake, GoogleTest); its eval is in experiments/cuda_checker/
docs/research/     research reports and notes behind the direction
data/              downloaded datasets (git-ignored, never committed)
```

**Rules**
- Reusable code goes in `src/objdel`, with tests.
- Experiment scripts only wire library code to data. They run as modules in an isolated interpreter, since data is
  untrusted: `.venv/bin/python -I -m experiments.phase1.<script> --help`.
- Gate scripts (`reward_gates`, `eval_demo`) reuse the test harness in `tests/oracles/`, so the gates measured on
  data and in pytest share one implementation. They run from a checkout: the editable install puts the repo root on
  `sys.path`.
- Results are committed next to the script that produced them.
- Data is never committed.

## Quick start

```bash
make setup        # .venv with Python 3.12 (uv), package editable with the train, data and dev extras; builds the C++ extension
make test         # fast test suite (rebuilds the extension first if cpp/ changed); make test-all adds the slow gates
make lint
make cov          # fast suite with line coverage of src/objdel (fails below 90 %, K13)
make build        # rebuild the C++ extension only
make fuzz         # native polycore fuzz under ASan/UBSan
make wheels       # manylinux wheels in Docker (cibuildwheel)
make cuda         # build the C++/CUDA engine (the CUDA checker)
make preflight    # prove data, checker and code on the Mac before any GPU is rented (CONFIG=...)
```

```python
import numpy as np

from objdel.checker import certify
from objdel.reward import ExpertCache, score_completion
from objdel.targets import build_target

y, x = np.ogrid[:100, :100]
raters = np.stack([(y - 50) ** 2 + (x - 50) ** 2 <= r * r for r in (20, 22, 24)])  # 3 raters, one nodule

report = certify(raters[0])                         # b0 parts, b1 holes, area, perimeter and a certificate
text = build_target({"nodule": raters[1]}, ("nodule",)).text   # rater 2 as canonical polygon JSON (grid 0-1000)
caches = {"nodule": ExpertCache.build(raters)}      # built once per item, offline
score = score_completion(text, False, ("nodule",), caches, 100, 100)
print(score["format_ok"], score["legal"], score["iou_tc"], score["topo_any"])   # 1.0 1.0 1.0 1.0
```

The datasets are LIDC-IDRI (TCIA) and RIGA. Each is downloaded separately under its own licence.

Rebuilding `data/` (store, prepared RL set, SFT targets) from the raw downloads uses the Phase 1 scripts in order
(`fetch_lidc_pixels`, `build_data`, `prepare_rl`, `build_targets`). The commands, timings and expected hashes are
in [RESULTS.md §8](experiments/phase1/RESULTS.md#8-phase-2-handoff).
