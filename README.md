<div align="center">

# ObjectDelineation

**Vision-language models that draw outlines you can trust.**

A small VLM outlines objects as exact polygons. A GPU-native geometry engine checks every outline,
and reinforcement learning teaches the model from those checks.

![License: MIT](https://img.shields.io/badge/license-MIT-blue)
![Python 3.12](https://img.shields.io/badge/python-3.12-3776AB)
![PyTorch](https://img.shields.io/badge/PyTorch-2.14-EE4C2C)
![CUDA](https://img.shields.io/badge/CUDA-13.0-76B900)

<img src="assets/hero.svg" alt="Image to VLM to polygon outline to checker reward" width="900">

</div>

## The problem

Segmentation models paint pixels. Vision-language models describe what they see. Neither gives you an outline
you can verify: a shape with the right number of parts and holes, valid geometry, and a measurable error against
what experts actually drew. In medical imaging, where experts often disagree, that gap decides whether a model
can be trusted.

## What we built

- **A model that writes geometry.** Qwen3-VL-2B, fine-tuned to return polygons (parts, rings and holes) as
  structured output, not masks or prose.
- **A checker that proves it.** A custom CUDA engine validates every polygon exactly and scores it against the
  full set of expert answers, including topology: how many parts, how many holes.
- **Training that rewards correctness.** Supervised fine-tuning, then GRPO reinforcement learning, where the
  checker is the reward. Generation, training and checking share one GPU.

## Results

Lung nodule outlining on LIDC-IDRI CT scans, scored on images the model never saw:

| | First RL run | RL run on shifted views (headline) |
|---|---|---|
| Overlap with experts (IoU, test, 2,200 images) | 0.645 | **0.735** |
| Off-centre nodules (IoU, edge validation view) | 0.540 | **0.692** |
| Valid outlines (test) | 99.4 % | **99.8 %** |

- **Test set (2,200 images):** IoU 0.735, boundary F-score 0.902, 99.8 % valid outlines.
- **Exact by construction:** 0 mismatches between the GPU engine and the reference checker on 195,809 cases.
- **Fast:** 7.7 M polygons checked per second on one RTX 5090, 141x the first GPU build.
- **Cheap:** the latest RL run trained 780 steps in 2 h 17 min on one GPU; the whole rental cost about $2.40.

<img src="assets/rl_gallery.png" alt="Held-out CT crops: expert outlines in green, the model's outline in red, before and after RL" width="900">

## How it works

<img src="assets/pipeline.svg" alt="Data, supervised fine-tuning, then the reinforcement-learning loop" width="900">

1. **Data.** Expert annotations become certified polygon targets, with patient-level train, validation and test splits.
2. **Fine-tune.** The model learns the output format and a first estimate of each shape.
3. **Reinforce.** The model writes several outlines per image. The checker scores each one, and the model learns
   from the differences.
4. **Evaluate.** Every claim is reported with patient-level confidence intervals on held-out data.

## Quick start

```bash
make setup     # Python 3.12 environment, builds the checker
make test      # test suite
```

```python
import numpy as np
from objdel.data.targets.build import build_target, min_areas
from objdel.verify.reward.cache import ExpertCache
from objdel.verify.reward.score import score_completion

y, x = np.ogrid[:100, :100]
experts = np.stack([(y - 50) ** 2 + (x - 50) ** 2 <= r * r for r in (20, 22, 24)])  # three expert outlines

outline = build_target({"nodule": experts[1]}, ("nodule",), None, min_areas({"nodule": experts})).text
score = score_completion(outline, False, ("nodule",), {"nodule": ExpertCache.build(experts)}, 100, 100)
print(score["legal"], score["iou_tc"], score["topo_any"])  # 1.0 1.0 1.0
```

## Repository layout

```text
cuda/            the geometry engine: one C++ source for the CPU build and the nvcc build
src/objdel/      the library, in layers; each imports only the layers below it
  infra/         rented GPU boxes, preflight, the data bundle, the private model repo, secrets
  rl/  sft/      reinforcement learning (GRPO) | supervised fine-tuning and checkpoint selection
  policy/        models, prompts, rows, caches, the sampler, evaluation
  data/          datasets, splits, the item store, targets, training windows
  verify/        reward, metrics, gates, the verifier API (no PyTorch)
  geom/          output grammar, canonical form, rings, topology, the engine wrapper (no PyTorch)
  core/          the typed run spec and provenance
specs/           run specs: what a GPU run does, checked before any rental
experiments/     thin command-line scripts, by phase
tests/           the test suite, oracles and goldens
reference/       records the code reads: frozen view hashes, validation halves, training exclusions, headline numbers
assets/          the images in this README
```

## Built with

PyTorch, Hugging Face Transformers and TRL, vLLM, and CUDA C++ with a single source for CPU and GPU builds.

## Citation

If you use this work, cite it as described in [CITATION.cff](CITATION.cff).

## License and data

The code is released under the [MIT License](LICENSE). Datasets are not included: LIDC-IDRI (TCIA, CC BY 3.0) and
RIGA are downloaded separately under their own licenses, and no images or annotations from either are in this
repository.
