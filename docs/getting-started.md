# Getting started

This guide takes you from a clean checkout to a first scored outline.
It takes about 10 minutes on a laptop. No GPU and no dataset are needed.

## Prerequisites

| Need | Version | Why |
|---|---|---|
| Python | 3.12 | `requires-python = ">=3.12"` in `pyproject.toml` |
| `uv` | any recent | `make setup` uses `uv venv` and `uv pip` |
| C++17 compiler, CMake 3.26+ | any | The checker builds as a Python extension |
| CUDA toolkit with `nvcc` | 13.0 | Only for `make cuda`. Not needed for the tests |
| NVIDIA GPU | sm_80 or newer | Only for the GPU tests and for training |

The repository was developed on an Apple-silicon Mac (CPU build of the checker) and on rented RTX 5090 boxes.
CI runs on Linux x86-64, Linux arm64 and macOS arm64 (`.github/workflows/tests.yml`).

## Step 1. Install

```bash
make setup
```

`make setup` runs these commands (see `Makefile`):

1. `uv venv --python 3.12` creates `.venv`.
2. `uv pip install -e ".[train,data,dev]"` installs the package in editable mode with all extras.
3. The install builds the C++ extension `objdel._gpucheck` with scikit-build-core and nanobind.

The extras are `train` (torch, transformers, trl, peft), `data` (pydicom, pylidc) and `dev` (pytest, ruff, hypothesis).
The pins are exact on purpose: reward and target hashes are recorded against these versions.

Check the build:

```bash
.venv/bin/python -I -c "from objdel import _gpucheck as g; print(g.build(), g.device_count())"
```

On a Mac this prints `cpu 0`. On a CUDA box after `make cuda` it prints `cuda 1` (or more devices).

## Step 2. Run the fast tests

```bash
make test        # pytest -m "not slow"
make lint        # ruff check src tests experiments
```

`make test` rebuilds the extension first when a file under `cuda/` or `CMakeLists.txt` changed.
Tests that need `data/` skip themselves when the data is absent.
Tests that need a GPU skip themselves when `nvidia-smi -L` lists no device.
`make test-all` adds the slow exhaustive oracle runs.
See [Testing](testing.md).

## Step 3. Build for CUDA (optional)

```bash
make cuda        # needs nvcc in $CUDACXX
CUDA_ARCHS="80;90" make cuda    # cross-compile for other cards
```

`make cuda` compiles the same `cuda/gpucheck.cpp` with nvcc and then asserts `g.build() == "cuda"`.
See [CUDA engine](reference/cuda-engine.md).

## Step 4. Score a first outline

Run this in a Python shell started with `.venv/bin/python -I`.
Three raters drew discs of radius 20, 22 and 24 pixels. A model answer is a 16-sided polygon.

```python
import math
import numpy as np
from objdel.geometry import canonicalize, serialize
from objdel.reward import ExpertCache, score_completion

y, x = np.ogrid[:100, :100]
raters = np.stack([(y - 50) ** 2 + (x - 50) ** 2 <= r * r for r in (20, 22, 24)])
caches = {"nodule": ExpertCache.build(raters)}          # built once per item, on the CPU

pts = [[int(500 + 220 * math.cos(2 * math.pi * i / 16)),
        int(500 + 220 * math.sin(2 * math.pi * i / 16))] for i in range(16)]
text = serialize(canonicalize({"nodule": [(np.array(pts), [])]}), ("nodule",))
print(text)
print(score_completion(text, False, ("nodule",), caches, 100, 100))
```

The real output of this code is:

```text
{"nodule":[{"outer":[[499,280],[584,296],[655,344],[703,415],[720,500],[703,584],[655,655],[584,703],[500,720],[415,703],[344,655],[296,584],[280,500],[296,415],[344,344],[415,296]],"holes":[]}]}
{'format_ok': 1.0, 'legal': 1.0, 'iou_tc': 0.9519, 'iou_min': 0.9519, 'f5': 1.0, 'g_area': 1.0,
 'topo_any': 1.0, 'topo_freq': 1.0, 'sig_id': 2367359347.0, 'error': '', 'defects': {}}
```

How to read it:

- The polygon has radius 220 grid units. The frame is 100 pixels wide, so one pixel is 10 grid units. The radius is 22 pixels.
- `format_ok` and `legal` are 1: the text parses and the checker finds no defect.
- `iou_tc` is 0.9519: the answer overlaps the best rater of the same topology by 95 %.
- `topo_any` is 1: every rater drew one part without a hole, and so did the answer.

Now break the outline on purpose. Reverse the point order so the ring runs counter-clockwise:

```python
bad = serialize({"nodule": [(np.array(pts)[::-1], [])]}, ("nodule",))
print(score_completion(bad, False, ("nodule",), caches, 100, 100)["defects"])
# {'wrong_winding': 1}
```

An outer ring must run clockwise on screen. The checker reports `wrong_winding` and `legal` becomes 0.
In training, the legality gate then sets every term of this answer to 0.

## Step 5. Tour the repository

```text
cuda/               gpucheck.h, gpucheck.cpp, gpucheck_py.cpp: the checker and its Python binding
archive/            retired code, not built: the CPU checker polycore
src/objdel/
  checker/          bit quads, Betti numbers, topology, polygonize, certify, 3-D volume
  geometry/         output grammar (parse, serialize), canonical form, polygon -> mask
  targets/          frame (crop policy), build (certified targets), augment (views), simplify
  reward/           ExpertCache, per-answer terms, group advantages, REWARD_V4
  data/             LIDC parser, RIGA extractor, splits, Parquet store
  train/            TRL configs, dataset rows, reward functions, FastGRPO, callbacks, vLLM serving
  eval/             metrics, hack audit, gates, probes, signal panel, vision check
tests/              pytest suite; oracles/ holds slow exact references; goldens/ holds recorded results
experiments/
  phase0/ phase1/   census, data build, reward gates (frozen records)
  phase2/           train_sft, train_rl, final_eval, preflight, configs/, vast/, docker/
  cuda_checker/     checker-versus-reference eval and the box script run_box.sh
docs/               these documents, figures/, research/
data/               downloaded data and built sets (git-ignored; treat as untrusted)
```

The layers depend on each other in one direction:

```text
  data  ->  targets  ->  reward  ->  train  ->  experiments
    \          \           ^            ^
     \          +-> geometry (parse, canonical form, mask)
      +-> checker (topology, certificates)       eval reads reward and geometry
```

Details are in the [architecture overview](architecture/overview.md) and the [package reference](reference/package.md).

## Step 6. Run an experiment script

Experiment scripts run as modules in an isolated interpreter.
The editable install puts the repository root on `sys.path`, so `experiments.*` imports work.

```bash
.venv/bin/python -I -m experiments.phase2.train_rl --help
.venv/bin/python -I -m experiments.phase2.final_eval --help
```

Always pass `-I` when a script reads `data/`. The data comes from the internet and is untrusted.
`-I` stops Python from importing code that sits next to the data.

## Next steps

| Goal | Read |
|---|---|
| Train a model | [Training guide](guides/training.md) |
| Rent a GPU | [GPU run runbook](runbooks/gpu-run.md) |
| Change the reward | [Reward reference](reference/reward.md). A change of `REWARD_V4` is a new reward and invalidates stored expert tables |
| Change the checker | [CUDA engine](reference/cuda-engine.md) and [Testing](testing.md) |
| Open a pull request | [Contributing](contributing.md) |
