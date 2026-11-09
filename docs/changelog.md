# Changelog

Milestones from the git log, newest first. Dates are commit dates. Numbers come from the commit messages and the results files named in each row.
The log has 87 commits up to `be52b1c` (2026-11-08). This page groups them. It does not list every commit.
Commits `a927ce9` (2025-03-20) and `793d926` (2025-03-21) are the initial repository commits.

## 2026-11-08: full RL run, final evaluation, speed work

| Commit | Change | Numbers | Source |
|---|---|---|---|
| `be52b1c` | Shared-prefix forward: each prompt runs once per group of 8 answers (`--shared-prefix`). RL plateau diagnosis. `train_rl --train-views` adds shifted and rescaled train windows. The bundle ships train windows only when named. CUDA streams and graphs review | R5 on an RTX 5090: 31.0 to 23.6 s per step (-24 %), same rewards. Plateau: RL rows have median object offset 0.018 of the side, the shifted exam 0.259 | `shared_prefix_2026-11-08.md`, `rl_plateau_diagnosis_2026-11-08.md` |
| `302ca3a` | Precision and recall against the majority outline. `final_eval` scores a checkpoint on a whole split: 8 samples per item, bootstrap 95 % intervals | Test (2,200 items): IoU_tc 0.645, boundary F 0.816, precision 0.695, recall 0.748 | `final_eval_2026-11-08.md` |
| `8851ea3` | README: before and after RL gallery of the model's polygons on held-out LIDC crops, picked by percentile | | `README.md` |
| `3fd83d8` | R4 speed test: 64-wide chunks and the no-merge vLLM sync | 8 % faster per step (32.66 to 29.94 s), identical rewards. vLLM sleep costs 1.6 s per step. Awake gains 0.05 s, within noise | `r4_speed_2026-11-08.md` |
| `255cb67` | Chunks round to 64 wide unless the text layers are compiled. Timers `time/sync_s`, `wake_s`, `sleep_s`. R4 config | | |
| `1870328` | Full RL run, 3 epochs. vLLM sync no longer writes the bf16 base weights (TRL #6688). Code review and utilisation research | 770 steps in 5 h 52 min, 27.5 s per step, about $4.56. Held-out IoU_tc 0.511 to 0.659 (+0.148, 95 % CI +0.115 to +0.178, better on 112 and worse on 14 of 128 items). Legal 94.5 % to 98.8 % | `full_rl_2026-11-08.md` |

## 2026-11-07: scoring moves to the GPU

| Commit | Change | Numbers | Source |
|---|---|---|---|
| `9f83c07`, `0a4b93f` | Vast image `vllm_native` in four layers. The one 4.14 GB layer timed out on 4 hosts | Largest layer 1.54 GB | `Dockerfile.vllm_native` |
| `edf3129`, `74a0135`, `cf2f507` | Phase 2 hard stop raised for the full RL run (owner): $15, then $16.60, $17.00, $17.50. Teardown of a box with no SSH skips the pull | Hard stop `HARD_STOP = 17.5` | `policy.py` |
| `0073f2d`, `c110022`, `2316212` | The full RL run: `train_rl`, resumable checkpoints every `--save-steps`, config `rl_full.yaml` | 770 steps = 8,208 rows / 32 x 3 epochs | `rl_full.yaml` |
| `a300197` | `ImageMemo` copies as the plain image processor. TRL deep-copies the processor on every save, which crashed the first save | | `fast.py` |
| `16b24da`, `ca44340` | R3 GPU test of the GPU scorer, `OpenOnce` and `ImageMemo` (RTX 5090, billed $0.440) | 40 GPU tests pass, 0 sanitizer errors, 56.6 to 36.3 s per step | `issues_solved_2026-11-07.md` |
| `6161da8` | `OpenOnce`: one stack of 256 rows per step. `ImageMemo`: one image preparation per distinct image | 4,096 opens to 32 (12.9 s to 0.14 s per step on the Mac). 2.12 s to 0.78 s per step for the processor | `issues_solved_2026-11-07.md` |
| `9c5dd2e` | Scoring on the GPU: IoU, F-score, area and topology on the device. Expert tables resident. polycore archived | Bit-identical to the host scorer | [ADR-0004](architecture/decisions/0004-scoring-on-the-gpu.md) |
| `5759224` (11-06), `791a52c`, `ed74aca` | R2 chain: the pre-registered G-probe on e2 | PASS: held-out reward 0.694 / 0.679 to 0.769, legal 94.5 % / 93.5 % to 98.2 %. Billed $0.550 | `issues_solved_2026-11-07.md` |
| `cc722aa`, `e617a6c` | Review fixes. Length-sorted update batches. No LoRA input cast | | `README.md` |
| `2fe5c15`, `e0fc6c4`, `69251d7` | K1: the CUDA checker under compute-sanitizer. Pin 13.4.92 | memcheck and initcheck: 0 errors on 40 GPU tests (RTX 3090) | `perf/k1_results.md` |
| `37158b9`, `aa45438`, `01e9162`, `b004ca3` | C10: skip TRL's old-log-probability pass. C5: `RopeCache`. 3-D voxel topology (`betti3d`, `topology_any`) | | `fast.py`, `checker/volume.py` |

## 2026-11-06: RL-ready

| Commit | Change | Numbers | Source |
|---|---|---|---|
| `647ec01` | P1 on an RTX 5090. FastGRPO is 1.81x TRL as shipped and on-policy. CUDA checker exact and fast | 80.3 to 44.3 s per step. Checker 7.74 M polygons per second (141x the first CUDA build, 87x the CPU loop). Billed $1.253 | `perf/p1_results.md` |
| `f68b32d` | `FastGRPO`, the GPU checker's host side (arena, pinned staging, two streams), profile arms, box chain | | [ADR-0005](architecture/decisions/0005-trl-with-fastgrpo-patches.md) |
| `96ed7b7` | TRL 1.17.1. The `vllm_native` image: training and vLLM in one venv. Perf fixes F1 to F3 | | [ADR-0006](architecture/decisions/0006-colocated-vllm-with-sleep-mode.md) |
| `e9f23d7` | MIT licence. RIGA-derived result tables leave the public repository | | `LICENSE` |

## 2026-11-02 to 2026-11-05: Phase 2 SFT and the CUDA checker

| Commit | Change | Numbers | Source |
|---|---|---|---|
| `b905e7e` (11-04) | Scoring on the CUDA checker only: batched `to_masks`, GPU build and check before every scoring chain | | |
| `87c0fb3` (11-03) | CUDA checker scored on an RTX 4090 (R0 v3). CMake arch fix | 195,809 of 195,809 cases identical to the CPU checker. First CUDA build 54,710 polygons per second at 16 vertices, batch 1,024 | `experiments/cuda_checker/RESULTS.md` |
| `999969b` (11-03) | GPU check before the bundle upload. Machine 28360 blocked | | `policy.py` |
| `734b1ec` (11-02) | Phase 2 reset after R1: single-source CUDA checker, preflight gate, coverage-fixed data, Qwen3-VL-2B | | [ADR-0001](architecture/decisions/0001-single-source-cuda-checker.md), [ADR-0008](architecture/decisions/0008-qwen3-vl-2b-over-lfm25-vl-450m.md) |
| `f8aa4ef` (11-02) | Phase 2 RESULTS: rung 1 (LFM2.5-VL-450M) frozen | Failed entry gates: median IoU_tc 0.495 shifted, 0.092 rescaled. Hole recall 0.5 % | `experiments/phase2/RESULTS.md` |
| `5d003e4`, `6faec26` (11-02) | `hf_push`: the selected branch to the private HF repo | | `hf_push.py` |
| `787c096` (11-01) | Vast offers default to North America | | [ADR-0009](architecture/decisions/0009-vast-ai-north-america-and-preflight.md) |
| `b37917b` (11-01) | R1 and C2 sample with HF generate. V1 R0 vLLM parity failed | | |

## 2026-10-18 to 2026-11-01: Phases 0 to 2 built

| Commit | Change | Numbers | Source |
|---|---|---|---|
| `1210dde`, `4f6159a`, `a6d8fb3` (10-20 to 11-01) | Phase 2 V1: R0 chain, stop rule, parity driver, Vast offer filter, image built from locks on GitHub Actions | | `experiments/phase2/V0_V1.md` |
| `d8891ce` (10-20) | Phase 2 stage 1: augmented views, SFT runner, vLLM bridge, eval gates, Vast tooling | | |
| `7516c9f` (10-20) | Phase 2 plan, judge reviews, lab evidence | | `experiments/phase2/PLAN.md` |
| `1c54bc4` (10-19) | Phase 1: polygon core, reward v4, data, TRL wiring, evaluation | K1: 0 of 849,802,293 pixels differ on 10,000 shapes. K2: 0 of 10,000 rings differ | `experiments/phase1/RESULTS.md` |
| `ab1b0b9` (10-19) | Phase 1 plan, lab evidence, research notes | | `experiments/phase1/PLAN.md` |
| `29868b9` (10-18) | Certified geometry checker, Phase 0 study, roadmap | 0 certificate failures on 26,820 RIGA masks and 40,189 LIDC slices. 114 of 1,877 nodules (6.1 %) with topology disagreement | `experiments/phase0/RESULTS.md` |
| `faaacc4` (10-18) | Repository restructured into library, experiments, cuda and docs | | |

## Totals that matter

| Quantity | Value | Source |
|---|---|---|
| Phase 2 committed spend | $15.59 of $17.50 after R4. R5 and later rows are in the ledger | `r4_speed_2026-11-08.md`, `vast/ledger.csv` |
| RL speed, TRL as shipped to current | 80.3 s per step to 23.6 s with `--shared-prefix` (different sampling lengths, see each source) | `perf/p1_results.md`, `shared_prefix_2026-11-08.md` |
| Tests | 1,484 collected (1,453 fast, 31 slow) at `be52b1c` | `pytest --collect-only` |
