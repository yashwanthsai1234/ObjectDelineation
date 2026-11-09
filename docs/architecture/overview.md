# Architecture overview

Status: current as of branch `gpu-scorer-open-once`, commit `be52b1c` (2026-11-08).
Decisions behind this design are in the [decision records](decisions/README.md).

## 1. Purpose and context

A model that outlines a lung nodule must get the shape right: how many parts, and whether there is a hole.
The system trains Qwen3-VL-2B (LoRA) to write polygon JSON.
It rewards each answer with an exact geometry checker.
The checker compares the answer with the set of outlines that expert raters drew.

Goals:

- The reward is exact. The same answer always gets the same score, on the CPU and on the GPU.
- The reward respects rater disagreement. A shape that any rater drew counts as correct.
- RL is fast enough to run 770 steps on one rented RTX 5090 for about $4.56 (`docs/research/notes/full_rl_2026-11-08.md`).

Non-goals:

- No serving stack. The project produces model folders, not an inference service.
- No multi-GPU training. FSDP and ZeRO paths are left to TRL and are not used.
- No redistribution of RIGA data (its licence forbids it).

## 2. System context

```text
   LIDC-IDRI (TCIA)     RIGA (Univ. of Michigan)
        |                      |                       untrusted downloads in data/raw
        v                      v
   +-----------------------------------+
   |  data: parse, extract, split      |   CPU, Mac
   |  -> Parquet store (exact masks)   |
   +-----------------+-----------------+
                     v
   +-----------------------------------+
   |  targets + reward tables          |   CPU, Mac
   |  frames, certified SFT targets,   |
   |  ExpertCache per object           |
   +-----------------+-----------------+
                     v   bundle (tar, SHA-256)
   +------------------------------------------------------------+
   |  rented box: RTX 5090, Vast AI, US or CA host              |
   |  +-----------+   +--------------+   +------------------+   |
   |  | SFT       |-->| RL (FastGRPO)|-->| final_eval       |   |
   |  +-----------+   +------+-------+   +------------------+   |
   |                         | answers (text)                   |
   |                         v                                  |
   |                 +---------------+                          |
   |                 | CUDA checker  | one call per step        |
   |                 +---------------+                          |
   +----------------------------+-------------------------------+
                                v   pull.sh (rsync + checksum)
                       model folder, logs on the Mac; HF private repo
```

External systems: TCIA and RIGA (data), Vast AI (boxes), GitHub Container Registry (image), Hugging Face Hub (model storage, private repo `svenk029/objdel-sft`).

## 3. Components

```text
  +-------------------------------- src/objdel ---------------------------------+
  |                                                                             |
  |  data ---------> targets ---------> reward <--------- geometry              |
  |  (store, lidc,   (frame, build,     (ExpertCache,     (text grammar,        |
  |   riga, splits)   augment,           score, advantages, canonical form,     |
  |                   simplify)          spec)              polygon -> mask)    |
  |       |               |                 |  ^                 |              |
  |       |               |                 |  |                 |              |
  |       v               v                 v  |                 v              |
  |    checker (bit quads, topology, polygonize, certify, volume)               |
  |                                         |                                   |
  |                                         v                                   |
  |                              objdel._gpucheck  <---- cuda/gpucheck.cpp      |
  |                                                                             |
  |  train (configs, dataset, rewards, fast, serve, callbacks, sample)          |
  |  eval  (metrics, audit, gates, probes, signal, vision, coverage)            |
  +-----------------------------------------------------------------------------+
  experiments/phase2: train_sft, train_rl, final_eval, preflight, configs, vast
```

| Component | Responsibility | Interface (deep modules) |
|---|---|---|
| `objdel._gpucheck` | Validate outlines, draw masks, trace rings, score terms. Exact integer and float64 arithmetic. | 9 functions: `check_polygons`, `check_masks`, `add_tables`, `drop_tables`, `resident`, `score_polygons`, `score_masks`, `build`, `device_count` |
| `objdel.geometry` | The single spec of model output text. Parse strictly. Canonicalise. Turn parts into masks. | `parse`, `serialize`, `canonicalize`, `to_masks` |
| `objdel.reward` | Reward v4. Hides the expert tables, the device handles and the gates. | `ExpertCache.build`, `score_batch`, `advantages_v4` |
| `objdel.targets` | Frames and certified SFT targets. | `frame`, `build_target` |
| `objdel.train.fast` | Speed patches of TRL's `GRPOTrainer`. | `FastGRPO` |
| `objdel.eval` | Exam metrics, hack audit, entry gates. | `metrics`, `hack_audit`, `gates.phase2_entry` |

The checker is a deep module in the sense of the `codebase-design` vocabulary.
A small interface (nine functions that take packed arrays) hides about 1,800 lines of kernels, memory management and exactness arguments.
The seam sits at `cuda/gpucheck.h`.
Two adapters sit behind it: the C++ compiler (CPU build) and nvcc (GPU build). A real seam has two adapters ([ADR-0001](decisions/0001-single-source-cuda-checker.md)).

Rules between components:

- `reward` and `eval` never rasterise outlines themselves. They call the checker.
- `geometry.text` is the only place that knows the output grammar. SFT targets, the reward and evaluation all import it.
- `train` modules import each other only through the names in `objdel/train/__init__.py`. Importing one does not import the whole training stack.
- Nothing under `src/objdel` imports from `experiments/`.

## 4. Data flow of one RL step

One step handles 32 images with 8 answers each: 256 answers.
Numbers are from `rl_full.yaml` and `fast.py`.
The time split is from R5 (`docs/research/notes/cuda_streams_graphs_review_2026-11-08.md`): 23.43 s per step, of which generation 11.11 s, scoring and forwards 4.65 s, update 7.67 s.

```text
  1. rows          OpenOnce hands TRL 256 rows (32 images x 8) and opens each image once
        |
        v
  2. generate      vLLM (colocated) writes 8 answers per image, T = 1, max 1,536 tokens
        |          vLLM gets one image pad per image; ids are checked against the trainer's
        |
        +----------------------------+
        v                            v
  3a. score (background thread)   3b. reference log-probs (trainer, GPU)
      parse 256 texts (CPU)           per-chunk windows, vision cache, rope cache
      pack rings (CPU)                shared-prefix packing when it wins
      ONE score_polygons call
      against resident tables (GPU)
        |                            |
        +-------------+--------------+
                      v
  4. advantages    advantages_v4 on 256 term vectors, group by group (CPU, NumPy)
                      |
                      v
  5. update        16 micro-batches of 16 rows, sorted by length, whole groups together
                   each chunk trimmed to its own width; the forward also sets the importance ratio
                      |
                      v
  6. sync          LoRA delta added in a new tensor and pushed to vLLM; vLLM sleeps, then wakes
```

Interaction order across processes and devices:

```text
  Trainer(py)      vLLM(GPU)         Reward thread(py)    CUDA checker(GPU)     Trainer(GPU)
     |  prompts ------>|                    |                    |                   |
     |<--- 256 answers-|                    |                    |                   |
     |  prefetch ------------------------->|                    |                   |
     |                                      |-- score_polygons ->|                   |
     |  reference forward ------------------------------------------------------->|
     |                                      |<-- per-object terms|                   |
     |<--- terms ---------------------------|                    |                   |
     |  advantages (CPU)                    |                    |                   |
     |  update forwards + backward ------------------------------------------------>|
     |  sync LoRA delta ------>|            |                    |                   |
```

Facts about this flow, with their source:

- The checker call releases the Python GIL, so the scoring thread overlaps the trainer (`gpucheck_py.cpp`, `nb::gil_scoped_release`).
- The reward wait is 0.01 s per step in the full run (`docs/research/notes/full_rl_2026-11-08.md`).
- `FastGRPO.skip_old_logps` drops TRL's old-log-probability pass. With one generation per optimizer step the weights do not move before the update, so the update's own detached log-probabilities are the old ones (`fast.py`).
- The weight sync uses `lora_merged_params`. It never writes the bf16 base weights ([ADR-0005](decisions/0005-trl-with-fastgrpo-patches.md)).

## 5. Where each part runs

| Part | CPU | GPU | Notes |
|---|---|---|---|
| Parse and build the store | yes | | `pylidc`, OpenCV, NumPy, Parquet |
| `ExpertCache.build` (EDT, tau, T5) | yes | | scipy `distance_transform_edt`. Offline, once per object |
| `ExpertCache.device()` | upload | resident | tables stay on the device until the cache is garbage-collected |
| Targets and frames | yes | | numba kernels, Pillow resize |
| Output parse (`parse`) | yes | | one regex plus `json.loads` |
| Pack rings into arrays | yes | | into a pinned staging buffer |
| Raster, validate, signature, IoU, F5, area | | yes | one warp per pixel row for the raster |
| Group advantages | yes | | NumPy, 256 values |
| vLLM generation | | yes | colocated, sleep mode |
| Model forward and backward | | yes | bf16 autocast, LoRA |
| Item-level exam metrics (pairwise IoU, NSD, GED) | yes | | read masks from `to_masks`; run a few times per run, not per step |
| Evaluation sampling | | yes | `final_eval` through the trainer's vLLM |

The rule "pixels never return to the host for scoring" is [ADR-0004](decisions/0004-scoring-on-the-gpu.md).
Only the exam metrics and the target builder copy masks to the host.

## 6. Boundaries and trust

| Boundary | Rule | Enforced by |
|---|---|---|
| Downloaded data | Untrusted. Python that reads it runs with `-I`. | `CLAUDE.md`, command lines in docs |
| RIGA data and rented hosts | RIGA never goes to a rented host without recorded permission and a datacenter offer. | `policy.governance`, `run_box.sh` refuses RIGA files |
| Secrets | Read from `<repo>/.secrets` (mode 600) or the macOS keychain. Never printed, never shipped to a box. | `objdel.secrets`, `sync_up` ships the `git ls-files` list, which skips ignored files |
| HF token on a box | Fine-grained, write access to one repo, passed on stdin only. | `hf_push.token_problems`, `policy.cmd_fetch_hf` |
| Money | A launch is refused above the hard stop. A launch above the alarm needs `--owner-approved`. | `policy.cmd_launch`, constants `ALARM = 10.0`, `HARD_STOP = 17.5` |
| Rental | No box without a passing, fresh preflight stamp. | `preflight.refusal` |
| Reward definition | `REWARD_V4` is frozen. A new value needs a new lab round and invalidates stored expert tables. | `spec.fingerprint`, `StaleCacheError` |

## 7. Failure modes and recovery

| Failure | Effect | Recovery |
|---|---|---|
| Answer does not parse or is truncated | `format_ok = legal = 0`, every term 0, `sig_id = -1` | Nothing. This is the reward signal |
| Answer parses but is illegal | Terms are kept for logs. The legality gate zeroes them in the advantage | Nothing |
| A scoring worker dies (OOM kill) | `BrokenProcessPool` within a second | The run fails loudly. Resume from the last checkpoint |
| Box lost during RL | Up to 50 steps lost | Same chain on a new box, `--resume` continues from the last `checkpoint-N` |
| Expert table built under another fingerprint | `StaleCacheError` | Re-run `prepare_rl`, or load with `rebuild=True` |
| Hacked reward (hedging, inflation) | Audit alarm | `Controller` rolls back to the best held-out checkpoint |

See [Troubleshooting](../runbooks/troubleshooting.md) for operational failures.

## 8. Scale and limits

| Quantity | Value | Source |
|---|---|---|
| Answers per step | 256 | `configs.grpo_config`: 16 x 16 micro-batches, G = 8 |
| Max answer length | 1,536 tokens (Qwen3-VL) | `configs.MAX_COMPLETION` |
| Checker caps | W, H <= 4,096; <= 10 rings per object; <= 1,000 vertices per ring; 2 to 32 raters per table | `gpucheck.h` |
| Resident tables | about 140 MB for all 7,341 LIDC train items | `docs/research/notes/gpu_scoring.md` |
| Step time | 27.5 s mean over the full run; 23.6 s with `--shared-prefix` (R5, 15 steps) | `full_rl_2026-11-08.md`, `shared_prefix_2026-11-08.md` |
| Checker throughput | 7.74 M polygons per second at 16 vertices, batch 1,024, RTX 5090 | `experiments/phase2/perf/p1_results.md` |

What to revisit as the system grows:

- A second GPU. The checker serialises calls with one mutex and keeps two lanes on one device.
- RIGA training. The RL path reads LIDC only (`train_rl` preloads `lidc`).
- Larger frames. The target builder refuses frames above 500 pixels (`MAX_SIDE`).
- RL on shifted windows. The plan exists in `rl_shifted.yaml` and is not yet run (see [Training guide](../guides/training.md)).
