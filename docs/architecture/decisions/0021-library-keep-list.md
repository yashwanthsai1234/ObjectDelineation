# ADR-0021: Keep our own parser, store, frame, metrics and Vast policy. Adopt no new framework

**Status:** Accepted
**Date:** 2026-11-09
**Deciders:** Repository owner

## Context

Implementation 2.0 searched for libraries that could replace repository code. Most hits were stateful systems or covered only part of the job.
This record lists each decision once, with its evidence and the task that owns a later change.

## Decision

### KEEP

- `geometry/text.py`: the strict grammar and error codes. Its regex also feeds vLLM structured output. No library gives both.
- `geometry/shape.py`: the hashed ring order (row-major start rule).
- `targets/simplify.py`: no library does lossless minimum-vertex simplification. Output is hash-pinned.
- `data/store.py`: byte-deterministic Parquet and a SHA-256 manifest.
- `Pillow`: bicubic resize in `targets/frame.py`. It is bit-identical across CPUs. Version pins stay as they are.
- `eval/signal.py`: patient-cluster bootstrap. SciPy has no cluster argument.
- `eval/metrics.py`: GED, HM-IoU and NSD. No maintained library implements multi-rater GED and HM-IoU.
- `K2 EDT`: our separable distance transform on the device (phase D). cuCIM is at most a second oracle in one box test.
- `vast/policy.py`: used with the `vastai` CLI 1.11.0 as an isolated `uv tool` (CLEAN13). The reason is the repo's governance: scoped key, image allowlist, spend ledger, cap hours, destroy-and-verify teardown, preflight stamp and the North America rule. SkyPilot can restrict Vast regions, so region filtering is not the reason.

### DEFERRED

- `trackio`: `report_to` live plots. Trigger: a run needs live plots.
- `events.jsonl`: run event log. Trigger: a run needs more than the ledger and the step log.
- `uv.lock`: `uv sync --locked`. Trigger: a second machine or CI.
- `pyright`: type checking. Trigger: the owner asks for it.
- `pre-commit`: with `ruff format` and a no-AI-attribution hook. Trigger: the owner asks for hooks.
- `pytest-xdist`: D8 runs without `-n`. Trigger: the Mac suite is too slow to use.
- `view-bank`: a cache of decoded views. Trigger: PNG decode above 0.2 s per step (C1).
- `pydantic`: run spec, owned by phase E (E6).
- `import-linter`: layer contracts, owned by phase E (E0).

### NOT ADOPTED

- `hydra-core`: no config sweep need. `omegaconf` goes with it.
- `structlog`: logs stay plain.
- `monai`: its metrics miss multi-rater GED.
- `surface-distance`: NSD is kept in `eval/metrics.py`.
- `medpy`: same reason.
- `skypilot`: no ledger, no preflight stamp.
- `dstack`: same reason.
- `modal`: same reason.
- `vastai-sdk`: deprecated in favour of `vastai`. Method return types are not documented.
- `cucim`: not for the production EDT.
- `nvimagecodec`: GPU image decode has no measured need.
- `nvidia-dali`: same reason.
- `nvcomp`: no measured need.
- `lance`: stateful store, no measured need.
- `webdataset`: no measured need.
- `dvc`: stateful system, no measured need.
- `lakefs`: stateful system, no measured need.
- `ray`: stateful system, no measured need.
- `triton`: stateful system, no measured need.
- `airflow`: stateful system, no measured need.
- `celery`: stateful system, no measured need.
- `redis`: stateful system, no measured need.
- `nats`: stateful system, no measured need.
- `prometheus-client`: the Prometheus stack has no measured need.

## Options considered

### Option A: Keep the code, list the decisions (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | None added |
| Risk | None: hashes and bytes stay as recorded |

### Option B: Replace parts with libraries

Each replacement changes bytes or drops a check the repo relies on (see KEEP). Not chosen.

## Trade-off analysis

A library helps only if it keeps the pinned hashes and the strict rules. None of the listed ones does today. Each DEFERRED row names the event that reopens it.

## Consequences

- Easier: the next phase reads one list instead of repeating the search.
- Harder: this list must change when a trigger fires.
- Revisit: when a DEFERRED trigger fires, or when a NOT ADOPTED library gains a documented need.

## Action items

1. [x] Record the list in this file and in the index.
