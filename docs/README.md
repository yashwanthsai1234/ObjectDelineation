# Documentation index

This folder holds the engineering documentation of ObjectDelineation.
The project trains a small vision-language model to draw outlines of lung nodules and optic discs.
A custom CUDA checker proves that each outline is legal and scores it.

Read the [glossary](glossary.md) first if a word is new to you.
The glossary fixes one word for each thing.
The other documents use those words only.

## Who reads what

| You are | Start here | Then read |
|---|---|---|
| A new engineer | [Getting started](getting-started.md) | [Architecture overview](architecture/overview.md), [Contributing](contributing.md) |
| A model trainer | [Training guide](guides/training.md) | [Reward](reference/reward.md), [Configs](reference/configs.md), [Evaluation](guides/evaluation.md) |
| A CUDA engineer | [CUDA engine](reference/cuda-engine.md) | [ADR-0001](architecture/decisions/0001-single-source-cuda-checker.md), [ADR-0004](architecture/decisions/0004-scoring-on-the-gpu.md), [Testing](testing.md) |
| The person who rents GPUs | [GPU run runbook](runbooks/gpu-run.md) | [Troubleshooting](runbooks/troubleshooting.md), [ADR-0009](architecture/decisions/0009-vast-ai-north-america-and-preflight.md) |
| A data engineer | [Data](reference/data.md) | [Output format](reference/output-format.md), [Package reference](reference/package.md) |
| A reviewer | [Architecture overview](architecture/overview.md) | [Decision records](architecture/decisions/README.md), [Changelog](changelog.md) |

## Map of all documents

| Document | Type | Purpose |
|---|---|---|
| [getting-started.md](getting-started.md) | Onboarding guide | Install, test, build, score a first outline, tour the repository |
| [architecture/overview.md](architecture/overview.md) | Architecture | System context, components, data flow of one RL step, CPU and GPU split |
| [architecture/decisions/](architecture/decisions/README.md) | ADRs | Eleven accepted decisions with options and consequences |
| [reference/package.md](reference/package.md) | Reference | Every module of `src/objdel` and its public API |
| [reference/cuda-engine.md](reference/cuda-engine.md) | Reference | The `gpucheck` engine: Python API, kernels, memory, exactness, build |
| [reference/reward.md](reference/reward.md) | Reference | Every reward term, gate and the way terms combine |
| [reference/output-format.md](reference/output-format.md) | Reference | The polygon JSON grammar, parser rules and legality rules |
| [reference/data.md](reference/data.md) | Reference | Datasets, store layout, splits, views, licences |
| [reference/configs.md](reference/configs.md) | Reference | The experiment YAML files and their schema |
| [guides/training.md](guides/training.md) | Guide | SFT then RL, flags, checkpoints, resume |
| [guides/evaluation.md](guides/evaluation.md) | Guide | `final_eval`, metric definitions, confidence intervals |
| [runbooks/gpu-run.md](runbooks/gpu-run.md) | Runbook | Preflight to teardown on Vast AI, cost caps, recovery |
| [runbooks/troubleshooting.md](runbooks/troubleshooting.md) | Runbook | Known failures and fixes |
| [testing.md](testing.md) | Strategy | Fast and slow tests, oracles, goldens, GPU tests, sanitizers, CI |
| [contributing.md](contributing.md) | Process | Workflow rules for code, commits and data |
| [glossary.md](glossary.md) | Glossary | One definition for each domain term |
| [changelog.md](changelog.md) | History | Milestones from the git log, newest first |
| [research/research2.0/Research2.0.md](research/research2.0/Research2.0.md) | Research | What the Research 2.0 study asked, how, and what it found |
| [research/research2.0/Direction.md](research/research2.0/Direction.md) | Plan | The phased plan with gates, costs and owner decisions |
| [research/research2.0/Architecture.md](research/research2.0/Architecture.md) | Architecture | The proposed target system (not yet built) |
| [research/research2.0/previous vs new.md](<research/research2.0/previous vs new.md>) | Comparison | The present system against the proposed system |
| [implementation2.0/Architecture-explained.md](implementation2.0/Architecture-explained.md) | Plan | The implementation architecture in plain English, with diagrams (read first) |
| [implementation2.0/Implementation.md](implementation2.0/Implementation.md) | Plan | The implementation plan for phases A-E, DOC and CLEAN: tasks, gates, costs, tests |

## What these documents do not replace

- `README.md` in the repository root is the public front page. It holds the figures and the results tables.
- `ROADMAP.md` is the plan. These documents describe what exists now.
- `experiments/*/PLAN.md`, `RESULTS.md` and `docs/research/` are the research record. They are frozen evidence.
  These documents link to them and do not copy them.

## How these documents stay true

- Every fact comes from code, a config, a test or a results file. Each document names its source.
- A number that changes with a new run lives in a results file. Documents quote it with its source path.
- When you change behaviour, change the matching document in the same pull request.
- The root `README.md` is checked by `tests/test_readme.py`. These documents are not checked by a test.
