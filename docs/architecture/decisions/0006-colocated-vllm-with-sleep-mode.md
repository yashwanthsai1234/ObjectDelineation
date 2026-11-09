# ADR-0006: Generate with vLLM on the same GPU, and let it sleep during the update

**Status:** Accepted
**Date:** 2026-11-06 (commit `96ed7b7`, "TRL 1.17.1; vllm_native image ...; perf fixes F1-F3")
**Deciders:** Repository owner

## Context

RL needs 256 sampled answers per step. Before this decision the G-probe sampled with Hugging Face `generate`.
`experiments/phase2/perf/issues.md` (decision D-a) measured the two choices:

- HF `generate`: about 1 s per answer, about 2.3 h for the G-probe, and out of memory at 256 sequences (20 to 51 GiB of KV cache).
- vLLM: about 0.04 to 0.09 s per answer, with a paged KV cache.

One RTX 5090 has 32 GB. vLLM and the trainer must share it.

## Decision

Use vLLM in colocate mode: the same process and GPU as the trainer.

- `vllm_mode="colocate"`, `vllm_gpu_memory_utilization = 0.35`, `vllm_max_model_length = 4096` (`configs.grpo_config`).
- `vllm_enable_sleep_mode = True`. vLLM frees its memory during the update and wakes for the next generation.
- Training and vLLM live in one virtual environment, in the image `vllm_native`. The box venv re-pins torch 2.13.0, transformers 5.19.0 and vLLM 0.32.0. The Mac pins torch 2.14.1 and transformers 5.20.1 (`Dockerfile.vllm_native`, `pyproject.toml`).
- SmolVLM-256M is the exception: it uses HF `generate` because TRL re-expands its image tokens under vLLM (TRL issue #6294) (`configs.py`).

## Options considered

### Option A: Colocated vLLM with sleep mode (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | Medium. One image with both stacks. Handover code in TRL |
| Cost | One GPU for the whole loop |
| Memory | Trainer peak 17 of 32 GB in the full run |

### Option B: HF `generate` in chunks

Too slow and too large (above). The G-probe's earlier plan used it (`experiments/phase2/PLAN.md` D13).

### Option C: A separate vLLM server on a second GPU

Not evaluated in the files read. The project rents one card at a time (`policy.cmd_launch` refuses a second instance without `--allow-parallel`).

### Option D: Colocated vLLM kept awake

Tested in R4 (`docs/research/notes/r4_speed_2026-11-08.md`). Step time 29.89 s awake against 29.94 s asleep: no gain.

## Trade-off analysis

Sleep mode costs 1.6 s per step in handovers: sync 1.28 s, wake 0.12 s, sleep 0.25 s (R4).
Awake mode saves 0.05 s per step, which is within noise. Sleep mode stays on because it frees memory.
The GPU use of 67 to 73 % comes from inside the generate and batch phases (CPU-side image processing and scheduling), not from the handovers.

## Consequences

- Easier: 256 answers in 7 to 12 s on one card.
- Easier: the same sampler path serves RL, `g_probe` and `final_eval` (`objdel.train.sample`).
- Harder: CUDA Graphs through `torch.compile(mode="reduce-overhead")` crash, because per-layer graph replay overwrites outputs and conflicts with vLLM's sleep allocator (`README.md`, P1 table).
- Harder: a harmless line `CUDA Error: invalid argument at cumem_allocator.cpp:249` prints at exit. It is vLLM's allocator shutting down. The exit code is 0 (`full_rl_2026-11-08.md`).
- Revisit: the idle GPU seconds between phases. See `docs/research/notes/grpo_utilization_research_2026-11-08.md`.

## Action items

1. [x] vLLM colocate in `configs.grpo_config`.
2. [x] Time sync, wake and sleep (`time/sync_s`, `time/wake_s`, `time/sleep_s`).
3. [ ] Re-measure when vLLM or TRL change their sleep path.
