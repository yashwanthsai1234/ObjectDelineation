# Phase 2 planning: research → think → test → plan

Phase 2 is SFT: teach the smallest vision-language model that can, to emit legal polygons (A0 outer ring vs A0′
hole- and part-aware), on Vast GPUs, and pass the Phase 2 entry gates (`experiments/phase1/RESULTS.md` §8,
`objdel.eval.gates.phase2_entry`). This file coordinates the planning agents. The plan itself goes in `PLAN.md`.

Effort split (owner's brief): 40 % online research, 25 % thinking / pattern-finding, 25 % testing hypotheses with
rigorous tests, 10 % writing the plan. Every proposal that reaches the plan is evidenced by a sandbox test.

## Inputs every agent reads first

- `ROADMAP.md` (Phase 2, Compute, GPUs, Kill criteria; the model ladder there is stale: Phase 1 D12 chose
  SmolVLM-256M → LFM2.5-VL-450M → Qwen3-VL-2B).
- `experiments/phase1/RESULTS.md` (what Phase 1 delivered; §8 Phase 2 handoff, entry gates, open issues).
- `experiments/phase1/PLAN.md` §2 decisions (D1–D13: output format D10, models D12), §3 architecture.
- `experiments/phase1/EXECUTION.md` §3 (library API: `objdel.train.{dataset,configs,rewards,callbacks,sample}`,
  `objdel.eval`, `objdel.targets`).
- `CLAUDE.md` (repo rules). Data in `data/` is untrusted: Python that reads it runs with `.venv/bin/python -I`.

## Task graph

```
 wave A: research (parallel, fresh contexts, web)            wave B: lab tests              wave C       wave D
 ┌──────────────────────────────┐
 │ A1 models + output tokens    │──┐
 ├──────────────────────────────┤  │   ┌───────────────────────┐
 │ A2 SFT methods + stack       │──┼──►│ B1..Bn: one sandbox   │   ┌──────────────┐  ┌───────────────┐
 ├──────────────────────────────┤  │   │ test per selected     │──►│ C plan draft │─►│ D Fable judge │─┐
 │ A3 Vast ops + GPUs           │──┤   │ hypothesis (Mac/MPS,  │   │ + arch + KPIs│  │ verdict       │ │
 ├──────────────────────────────┤  │   │ Docker amd64, tiny    │   └──────────────┘  └───────────────┘ │
 │ A4 eval, baselines, data,    │──┤   │ SFT runs)             │          ▲                            │
 │    licence, augmentation     │  │   └───────────────────────┘          └──── loop until KPIs pass ──┘
 ├──────────────────────────────┤  │
 │ A5 decoding + serving        │──┤
 ├──────────────────────────────┤  │
 │ A6 performance: C++ vs Python│──┘
 └──────────────────────────────┘
```

| Task | Question | Notes file (`docs/research/notes/Phase 2 SFT/`) |
|---|---|---|
| A1 | Which small VLMs and which coordinate/polygon output representation (text digits, special/bin tokens, quantization, ordering) give the best legal-polygon SFT? What do 2025–2026 VLM segmentation systems do? | `A1_models_tokens.md` |
| A2 | SFT recipe and stack: losses (coordinate-aware, number-token losses), multi-rater targets, curriculum, LoRA vs full, packing, kernels (Liger, FlashAttention), frameworks (TRL, Unsloth, LLaMA-Factory, ms-swift, Axolotl, torchtune) | `A2_sft_stack.md` |
| A3 | Vast.ai: how we use it (CLI/SDK, templates, Docker images, on-demand vs interruptible, storage, data transfer, security, checkpointing), and 5 GPUs for Phase 2 with memory/throughput/price evidence | `A3_vast_gpus.md` |
| A4 | Evaluation and data: multi-rater metrics, MedSAM/SAM-family baselines on our splits, RIGA licence vs rented machines, augmentation that keeps targets certified, precision probe | `A4_eval_data.md` |
| A5 | Decoding and serving for SFT eval / sampling: constrained decoding (grammar engines), vLLM / SGLang support for our models, K-sample throughput | `A5_decoding_serving.md` |
| A6 | Performance: where Phase 2's CPU time goes (collation, augmentation re-targeting, scoring eval samples); which modules justify C++ (nanobind, like polycore) | `A6_performance.md` |

Each research note ends with: ranked recommendations, rejected alternatives with reasons, **hypotheses to test**
(each with a concrete sandbox test design and pass criterion), and a source list (URL per claim).

## Rules

- Research online for every non-trivial claim (WebSearch/WebFetch via ToolSearch `select:WebSearch,WebFetch`);
  prefer primary sources (papers, official docs, source code, release notes) dated 2025–2026; cite URLs.
- Sandbox: `<scratchpad>/phase2/<task>/`; lab code that becomes evidence goes to `experiments/phase2/lab/<task>/`.
  Never edit `src/`, `tests/`, `experiments/phase1/` (Phase 1 is closing). No git state changes.
- Mac: 8 cores, 16 GB shared by all agents: ≤ 2 processes and ≤ 3 GB per agent (an MPS training test: ≤ 7 GB,
  one at a time, coordinate through the coordinator). Docker (linux/amd64 emulation) is available; no CUDA locally.
- Wait loops must not match their own command line (wait on a PID, or `pgrep -f "[p]attern"`).

## Wave B: hypothesis tests

Wave A notes (`docs/research/notes/Phase 2 SFT/A1–A6`) proposed ~45 hypotheses. Tests that decide the plan run
now; the rest become first-Vast-session A/Bs in the plan. The Mac has one GPU (MPS): MPS training runs one at a
time, in this order; CPU/Docker tasks run in parallel.

```
 MPS slot (serial):  B1 model learnability ──► B2 SFT depth → GRPO signal (+ interpolation) ──► B6 A0′ hole
                     (H-A1-1, H-A1-2)          (H-A2-3, H-A2-7, H-A2-8, H-A1-8)                 oversampling + NTL
                                                                                                (H-A2-5/H-A1-6, H-A2-2)
 CPU / Docker:       B3 probes + offline augmentation (H-A1-7, A4 precision probe, H-A6-6, H-A6-9)
                     B4 reward speed (H-A6-2 forked pool, H-A6-4 preload, H-K3SYN, H-A6-1 harness)
                     B5 environments (H-A5-7 Phase 3 env, H-A5-1/2 vLLM parity, Vast images from A3)
```

| Task | Owner paths (lab) | Report |
|---|---|---|
| B1 | `experiments/phase2/lab/B1/` | `experiments/phase2/lab/B1/REPORT.md` |
| B2 | `experiments/phase2/lab/B2/` | `…/B2/REPORT.md` |
| B3 | `experiments/phase2/lab/B3/` | `…/B3/REPORT.md` |
| B4 | `experiments/phase2/lab/B4/` | `…/B4/REPORT.md` |
| B5 | `experiments/phase2/lab/B5/` | `…/B5/REPORT.md` |
| B6 | `experiments/phase2/lab/B6/` | `…/B6/REPORT.md` |

Every report: hypothesis, test design, rigorous pass criterion fixed before running, results with confidence
intervals where sampling is involved, verdict (accept / reject / inconclusive), what it changes in the plan, and
sources. Lab code is ruff-clean and run as modules with `-I` when it reads `data/`.

**Scope change (2026-10-20, coordinator):** B6 (A0′ hole oversampling, number-token loss, crop-jitter training)
moves to pre-registered A/Bs in the first Vast sessions (plan §Experiments): under the shared Mac's load an MPS
pilot ran only 36 steps per arm (B1), too few to separate these effects. B2 (SFT depth → GRPO signal) still runs on
MPS after B1 because it sets the checkpoint-selection rule. Wave C (plan draft) started in parallel with B1/B2;
its model and recipe sections are filled from B1/B2.

**Status (2026-10-20):** B1 done (`lab/B1/REPORT.md`), written into PLAN revision 5: ladder kept, vision tower
frozen under LoRA plus a G-vision gate, and crop-jitter training mandatory in the recipe bundle (not an A/B: every B1
arm drew at the frame centre). B2 done (`lab/B2/REPORT.md`), written into PLAN revision 6: D12 rule R-B (deepest
eligible branch, sd_geom veto), base interpolation rejected as a repair, signal panel folded into R1 / R2. Wave B is
complete.
