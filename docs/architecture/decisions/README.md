# Architecture decision records

Each record follows one template: Status, Date, Deciders, Context, Decision, Options considered, Trade-off analysis, Consequences, Action items.
All records are **Accepted**. They describe decisions already taken and visible in code, configs, results files and the git log.
The date is the date of the commit that put the decision into the repository.
Deciders is the repository owner. The owner's name appears in commit messages as "owner" or "owner decision".

To add a record, copy the template of any file here, use the next number, and add a row below.
Record a decision only if it is hard to reverse, surprising without context, and a real trade-off.

| ADR | Title | Date | Status |
|---|---|---|---|
| [0001](0001-single-source-cuda-checker.md) | One source file builds the checker for the CPU and for nvcc | 2026-11-02 | Accepted |
| [0002](0002-set-valued-multi-rater-reward.md) | Score against the set of rater outlines, not one mask | 2026-10-19 | Accepted |
| [0003](0003-exact-integer-geometry.md) | Exact integer geometry on a 0-1000 grid | 2026-10-19 | Accepted |
| [0004](0004-scoring-on-the-gpu.md) | Score on the GPU. Masks never return to the host | 2026-11-07 | Accepted |
| [0005](0005-trl-with-fastgrpo-patches.md) | Keep TRL and patch it (FastGRPO) | 2026-11-06 | Accepted |
| [0006](0006-colocated-vllm-with-sleep-mode.md) | Generate with vLLM on the same GPU, with sleep mode | 2026-11-06 | Accepted |
| [0007](0007-shared-prefix-forward.md) | Run each prompt once per group of 8 answers | 2026-11-08 | Accepted |
| [0008](0008-qwen3-vl-2b-over-lfm25-vl-450m.md) | Use Qwen3-VL-2B with LoRA. Do not train LFM2.5-VL-450M again | 2026-11-02 | Accepted |
| [0009](0009-vast-ai-north-america-and-preflight.md) | Rent from North American Vast AI hosts only. Preflight first | 2026-11-01 | Accepted |
| [0010](0010-compact-json-grammar-and-strict-parser.md) | Compact JSON on a 0-1000 grid with a strict parser | 2026-10-19 | Accepted |
| [0011](0011-lora-weight-sync-without-merging.md) | Push LoRA weights to vLLM as W + delta, never merged in place | 2026-11-08 | Accepted |
| [0021](0021-library-keep-list.md) | Keep our own parser, store, frame, metrics and Vast policy. Adopt no new framework | 2026-11-09 | Accepted |

## How the decisions fit together

```text
  0003 exact integer geometry
        |
        v
  0001 one source, CPU + nvcc  ---->  0004 scoring on the GPU
        |                                   |
        v                                   v
  0002 set-valued reward  <-----------  0005 TRL + FastGRPO  --> 0006 colocated vLLM
        ^                                   |                       |
        |                                   +--> 0007 shared prefix |
  0010 output grammar                       +--> 0011 LoRA sync     |
        |                                                           |
  0008 Qwen3-VL-2B (LoRA)  --------------------------------------> 0009 Vast AI + preflight
```
