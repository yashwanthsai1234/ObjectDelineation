# SD2: Production readiness and evolvability (system design, up to Phase E)

Author: fresh system designer, no prior opinion. Date 2026-11-08. Lens: production readiness and evolvability.
Scope: module boundaries, versioned verifier API, typed run specs, reproducibility, observability, CI with a GPU job,
the path to 2+ GPUs (async RL) and to 3D voxels, with GPU everywhere. Every fact is tagged
[measured: file], [code: file:line], [web: URL], [paper: arXiv id] or [estimate]. REPO = the project root.

---

## 0. Skills used (and what each changed)

| Skill | What it changed in this design |
|---|---|
| system-design-skills:design / system-design (orchestrator) | Forced the loop: requirements, numbers, API, trade-offs, failure, evolution; the design doc layout below; the final score with the quality bar. |
| system-design-skills:service-decomposition | Kept a **modular monolith**. Only one seam is real: the verifier, and only once a second process needs it. No gateway, no mesh. |
| system-design-skills:api-design | Concrete verifier contract: per-item error envelope with `retryable`, additive versioning, batch caps, natural idempotency (pure function, no key needed), float encoding rule. |
| system-design-skills:consistency-coordination | Framed async RL as a consistency problem: every sample carries `policy_version`; the buffer bounds staleness; a reward carries `reward_id`. Resume semantics chosen (drop buffer, exactly-once per step). |
| system-design-skills:scaling-evolution | Climb one rung at a time: 1 GPU in turns, then 2 GPUs async, then N. Named the breaking point of each rung and the gate to climb. |
| system-design-skills:resilience-failure | Inverted the usual degradation story: in RL a wrong reward is worse than a stopped run, so the verifier **fails closed**. Retries only for eval and uploads. |
| system-design-skills:observability | Training SLOs (step time, GPU busy, IS share at 0, cost per run) and "watch the watcher" (heartbeat, NVML sampler). Symptom guards stop the run; cause metrics are logged only. |
| system-design-skills:distributed-logging | YAGNI for a pipeline: one host, structured JSONL per process with `run_id/step/component`, pulled with the run. Revisit only at multi-host. |
| system-design-skills:back-of-the-envelope | Sized the verifier load (about 11 answers/s now, about 1,300/s at 64 GPUs), the doctor-store memory for 3D, and the cost of a GPU CI job. |
| system-design-skills:architecture-diagram | Rendered the component view to `SD2_production_diagram.html` (same folder). |
| engineering:system-design, engineering:architecture | ADR-style decision rows (section 9) and the "what to revisit as it grows" list. |
| engineering:documentation | Doc layout: most useful facts first; link, do not copy (pointers to REPO files instead of repeating them). |
| engineering:code-review | Read the reward, cache, config and run-record code as a reviewer: found the provenance gaps (section 5.1) and the layer leaks (section 3.2). |
| engineering:debug | Reproduce-isolate-diagnose framing for the failure table (section 7): each failure has a detector, not only a fix. |
| huggingface-skills:huggingface-best | No model change proposed here; kept the model as a field of the run spec so a 4B swap is a spec change, not a code change. |
| nvidia-skills:nvidia-skill-finder | Checked relevance (CUDA, NeMo-RL). No catalog lookup run and no skill installed; NeMo-RL facts below come from NVIDIA docs, not a skill. |

Web research (2025-2026 RL post-training systems) is in section 2.

---

## 1. Requirements

### 1.1 Functional (in scope, Phase A to E)
1. Train a VLM (now Qwen3-VL-2B + LoRA, any model up to 4B later) to write outlines as text: SFT, then GRPO.
2. Score every answer exactly against the SET of doctors (legality, IoU_tc, F5, area, topology) on the GPU.
3. Evaluate a checkpoint on a whole split with fair metrics and nodule-level CIs; save the answers so a new reward
   can re-score them without new generation.
4. Offer the verifier as a reusable, versioned library (in process) and, when needed, a service (out of process).
5. Run every GPU job from one typed run spec, with full provenance in every artefact.
6. Grow to 2+ GPUs (async RL) and to a 3D voxel shape kind without a change to reward maths.

### 1.2 Non-functional
| Quality | Target | Source |
|---|---|---|
| Exact | CPU and GPU builds give identical bits for every term | Architecture.md §1.2; golden files |
| GPU first | GPU busy >= 80 % of a step; CPU only where unavoidable (PNG decode, text parse) | owner rule; Architecture.md §1.2 |
| Reproducible | Every artefact names code, checker build, reward, data, view, model, image, spec | owner rule; today only SFT does it (5.1) |
| Cheap and safe | $3 GPU budget for scripts; US/CA hosts; preflight + owner approval before rental | BRIEF §3 |
| Evolvable | Trainer framework swap (TRL -> verl/NeMo-RL) touches only the `rl` adapter; 3D adds a shape kind, not a new API | this design |
| Simple | Prefer a library over own code; archive old code | owner rule |

### 1.3 Out of scope
- An inference/serving product (model folders only) [code: REPO/docs/architecture/overview.md non-goals].
- CAD solids (parked). Physical AI assets after Phase E: only keep the shape-kind seam open.
- Multi-region, HA, autoscaling: one rented host at a time. A verifier "outage" stops a run; it is not a user outage.

### 1.4 Load-bearing questions (asked here; defaults assumed so the design can proceed)
| # | Question for the owner | Default I assumed | What changes if the answer differs |
|---|---|---|---|
| Q1 | Does any program other than our trainer need the verifier before Phase E? | No | If yes: build the HTTP adapter in Phase E, not at stage 2. |
| Q2 | 2 GPUs on one host (2x5090) or two hosts? | One host | Two hosts: NCCL over network, verifier per host, RIGA rule blocks remote scoring of RIGA. |
| Q3 | Max policy staleness in async RL? | 1 step | Larger: needs a stronger off-policy correction (CISPO/GSPO) and a measured reward-curve match. |
| Q4 | May run metrics leave the host (private HF dataset / trackio Space)? | Yes, metrics only, never data | No: JSONL only, pulled at the end. |

---

## 2. What top labs do (web research, 2025-2026) and what it means here

| System | Structure (fact) | Lesson for this project |
|---|---|---|
| verl v0.7 | Single controller (MPMD) orchestrates; SPMD workers inside; **rollout server mode by default**, SPMD rollout removed [web: https://verl.readthedocs.io/en/latest/blog/v0.7.html]. Fully async: Rollouter -> MessageQueue -> Trainer, NCCL param sync, `async_training.staleness_threshold`, partial rollout [web: https://verl.readthedocs.io/en/latest/advance/fully_async.html]. Reward: `RewardLoopManager` with `RewardWorker`s, streaming reward as each sample finishes; custom reward signature `compute_score(data_source, solution_str, ground_truth, extra_info)` [web: https://verl.readthedocs.io/en/latest/advance/reward_loop.html ; https://verl.readthedocs.io/en/latest/preparation/reward_function.html] | Our verifier must work **per sample, streaming**, keyed by a data id (verl passes `ground_truth`/`extra_info`, so our `key` goes in `extra_info`). Batch-joining is an optimisation inside the verifier, not a requirement on the caller. |
| NeMo-RL + NeMo Gym | Async GRPO: background trajectory collector, replay buffer with weight-version tags, `max_trajectory_age_steps`, IS correction **required**, vLLM async needs non-colocated GPUs [web: https://docs.nvidia.com/nemo/rl/latest/guides/async-grpo.html]. Verifiers live in FastAPI "Resources Servers" with `/verify`; Gym is CPU-only [web: https://docs.nvidia.com/nemo/gym/about/architecture ; https://docs.nvidia.com/nemo/rl/latest/design-docs/nemo-gym-integration.html] | Version-tag every sample. Our verifier is GPU-bound, so it cannot sit in a CPU-only Gym server: it runs as a GPU actor with a thin HTTP face. |
| TRL AsyncGRPOTrainer | Rollout worker is a separate process with `CUDA_VISIBLE_DEVICES=""`; GPU rewards must be served over HTTP; **images not supported yet**; LoRA sync since v1.14 [web: https://huggingface.co/docs/trl/async_grpo_trainer ; https://huggingface.co/blog/asyncgrpo-lora-hfjobs] | Confirms Architecture.md D1: TRL async cannot do our VLM task. It also shows the general rule: once async, the GPU verifier becomes a service. |
| AReaL | Fully async, interruptible rollout workers, staleness bound eta enforced by throttling generation requests, "parallel reward service", up to 2.77x speed-up [paper: arXiv 2505.24298] | Staleness is a controller rule, not a trainer rule. Reward is a separate service. |
| PipelineRL | In-flight weight updates; generation engine needs only 3 HTTP endpoints [paper: arXiv 2509.19128] | Keep the sampler behind a small protocol (generate, init weight group, update weights). |
| slime | Megatron + SGLang + **Data Buffer** as the one bridge; reward/verifier runs inside rollout; `staleness/mean`, `staleness/max` metrics; decomposition `log p_vt - log q_vs` = parameter staleness + same-version engine mismatch [web: https://thudm.github.io/slime/advanced/rl-systems.html ; https://github.com/thudm/slime] | Log staleness and engine mismatch **separately**. Our measured mismatch (IS min = 0 in 98.1 % of steps) is pure engine mismatch at age 0 [measured: REPO/runs/rl/train_rl/rl.json via BRIEF §2]; async adds the other term on top. |
| OpenEnv (Meta + HF) | Gymnasium-style `reset/step/state` over HTTP in Docker; "protocol layer, not a reward framework"; rewards computed inside the env via Rubrics [web: https://huggingface.co/blog/openenv-agentic-rl ; https://github.com/huggingface/OpenEnv/blob/main/rfcs/002-env-spec.md] | If an external user wants our verifier, an OpenEnv-compatible single-step env wrapper is the cheapest interop. Do not invent our own env protocol. |
| Tinker (Thinking Machines) | Client loop on CPU; primitives `forward_backward`, `optim_step`, `sample`, `save_state`; LoRA only; RL losses incl. `importance_sampling`, `ppo`, `cispo` [web: https://tinker-docs.thinkingmachines.ai/tinker/sdk-cheatsheet/] | The clean seam in 2026 is: **loop owns algorithm; engine owns forward/backward/sample**. Our `rl` layer should look like that so the engine (TRL now, verl later) is swappable. |

Common pattern across all of them: (1) rollout, train and reward are separate roles even when colocated;
(2) samples carry a policy version and staleness is bounded by the controller; (3) reward/verifier is a service
behind a narrow, typed request; (4) typed config (Hydra/OmegaConf YAML in verl and NeMo-RL).

---

## 3. Module boundaries (library vs experiments)

### 3.1 Proposed layout (modular monolith, two distributions from one repo)

```
REPO/
  src/objdel/
    core/      NEW. Types only, no deps: RunSpec, VersionStamp, ArtifactRef, ids, errors.   (pydantic v2)
    geom/      engine binding (_gpucheck), parse, raster, polygonize, labels, topology.    (numpy only)
    verify/    reward specs (v4, v5), DoctorStore (device tables), Verifier.score(),
               metrics (GPU), stats (CIs), http.py (thin FastAPI face, optional extra).     NO torch / TRL
      build/   doctor-store BUILDER (EDT, tau, T5). GPU EDT target; scipy only as oracle.
    data/      store, LIDC/RIGA readers, splits, views, targets, quality gates, leakage checks
    policy/    model registry, prompt, processor settings, Sampler protocol (vLLM, HF)
    rl/        algorithm (advantage, loss config), TrainerAdapter protocol, adapters/trl.py
               (+ fast.py patches), later adapters/verl.py
    sft/       supervised training
    runs/      NEW. spec resolve/lock, provenance stamping, run record, guards, CLI entry points
    infra/     Vast host policy, preflight, bundle, cost ledger (moved from experiments/phase2/vast)
  experiments/<phase>/   research scripts and lab notes only; frozen after the phase; archived to
                         archive/<phase>/ when superseded. Runbooks never call experiments/.
  specs/                 NEW. typed run specs (YAML validated by objdel.core.RunSpec), one per run
  archive/               old code, read-only, excluded from lint/test/coverage
```

Distributions (uv workspace or two `pyproject` targets, same cibuildwheel job that exists today
[code: REPO/.github/workflows/tests.yml wheels job]):
- `objdel-verify` = core + geom + verify. Deps: numpy + the extension. Installable in a verl/NeMo Ray worker.
- `objdel` = everything, extras `[train]`, `[data]`, `[serve]`.

Layer rule (one direction): core <- geom <- verify <- data <- policy <- {rl, sft} <- runs <- infra.
Enforced with **import-linter** "layers" contract in `make lint` (a library, not our own checker).
Note: `data` sits above `verify` because view building needs the verifier (certified targets, per-view doctor masks).

### 3.2 Why change today's layout (findings)
- The reward package (future `verify`) imports `cv2` and `scipy` at module level for the cache build
  [code: REPO/src/objdel/reward/spec.py:17-19; REPO/src/objdel/reward/cache.py imports]. The run-time scorer needs
  neither. Splitting `verify.build` keeps `objdel-verify` light and isolates the CPU EDT that the owner wants on GPU.
- Eval metrics use `scipy.ndimage` and `linear_sum_assignment` on the CPU and a fork pool
  [code: REPO/src/objdel/eval/metrics.py:94-95; REPO/src/objdel/train/serve.py docstring]. GPU-first means these move to
  `verify.metrics` on the engine; the Hungarian match on <= 32 doctors can stay host-side (a few hundred numbers).
- `advantages_v4` lives in the reward package [code: REPO/src/objdel/reward/__init__.py:7]. Group statistics are
  algorithm, not verification; verl and NeMo-RL compute their own advantages from a per-sample score. Move to `rl`.
- Run plans are argv lists in YAML [code: REPO/experiments/phase2/configs/rl_shifted.yaml `steps:`], not typed. The
  ops scripts that cost money (`vast/policy.py`, `preflight.py`) live under `experiments/` beside lab notes and
  invoices [measured: `ls REPO/experiments/phase2`]. Production code that spends money needs library status and tests.

---

## 4. Estimates (numbers that decide the design)

| Quantity | Value | Basis |
|---|---|---|
| Answers per step | 32 prompts x 8 = 256 | [code: REPO/src/objdel/train/configs.py grpo_config] |
| Step time now | 23.6 s (R5, shared prefix) | [measured: BRIEF §2 / R5] |
| Verifier load now | 256 / 23.6 = about 11 answers/s; about 4 KB/s at about 400 B/answer | [estimate: 131 tokens mean x ~3 B/token] |
| Verifier load at 64 GPUs async | about 1,300 answers/s, about 1.3 MB/s | Architecture.md §4.2 worked example 2 |
| Conclusion | Plain HTTP/JSON or an in-process call is enough at every stage up to Phase E. No gRPC, no streaming protocol, no queue product. | api-design: escalate only when a number forces it |
| Doctor tables per 2D object | k*H*W + 3*H*W bytes (u8 planes) + 4*H*W (f32 distance); S <= 500 | [code: REPO/cuda/gpucheck_py.cpp:163-192; REPO/src/objdel/targets/frame.py:6-7] |
| 2D worst case, k=4, S=500 | 11 x 250,000 B = about 2.75 MB/object; 8,208 train items x 3 views = 24,624 rows -> up to about 68 GB if all were S=500 and all resident | [estimate]; most LIDC frames are far smaller than 500, but the store has **no byte budget** today (lazy upload, resident until the Python object dies) [code: REPO/src/objdel/reward/cache.py:241-256] |
| 3D object, 128^3 voxels, k=4 | 2.1 M voxels x 11 B = about 23 MB/object; 1,000 objects = about 23 GB | [estimate] -> a byte-budgeted, LRU device store is **required** before 3D |
| 3D answer length | about 10 slices x 131 tokens = about 1,300 tokens per answer | [estimate]; near today's 1,536 cap, so decode tail grows and generation dominates: this is what forces async |
| GPU CI cost | GitHub GPU larger runner (T4, 4-core) about $0.07/min = $4.20/h, needs a Team/Enterprise plan [web: https://runs-on.com/features/gpu-runners ; https://docs.github.com/en/enterprise-cloud@latest/billing/reference/actions-runner-pricing]; a Vast RTX 5090 about $0.60/h [code: REPO/experiments/phase2/configs/rl_shifted.yaml header] | GPU CI = gated, path-filtered, on a Vast box with the existing preflight, about 10 min = about $0.10 per run [estimate] |

---

## 5. Data and version model (reproducibility)

### 5.1 What exists today (code review)
| Artefact | Carries versions? | Evidence |
|---|---|---|
| SFT output manifest | Yes: config, view hash, code SHA-256, step plan | [code: REPO/experiments/phase2/train_sft.py:28, 287-307] |
| Item store | Yes: SHA-256 per Parquet file, numpy/pyarrow versions | [code: REPO/src/objdel/data/store.py:164-167] |
| Expert cache | Yes: reward fingerprint (REWARD_V4 constants, size rule, numpy/cv2/scipy versions) + content SHA-256 | [code: REPO/src/objdel/reward/spec.py fingerprint(); cache.py docstring] |
| RL run record | **No**: keys `steps, lr, ckpt, resumed_from, final, timing_s, log` | [measured: REPO/runs/r5/shared_rl.json] |
| Final eval | **No**: keys `ckpt, view, split, K, seed, items, legal, iou_tc, f5, precision_vote, recall_vote, dice_vote, seconds`; `ckpt` is a host path `/workspace/runs/r5/final_model` | [measured: REPO/runs/r5/eval/eval_test/final_eval.json] |

So the headline result (test IoU 0.645) cannot be tied by its own file to a code commit, a reward version or a
data digest. This is the most important production gap.

### 5.2 Identifiers (content-addressed, short = first 12 hex of SHA-256)
```
code_id     = git commit + dirty flag + sha256(src tree)         (dirty => run refused on a rented host)
checker_id  = _gpucheck build: "cuda"|"cpu" + CUDA archs + sha256(cuda/*.cpp,*.h) + nvcc version
reward_id   = sha256(canonical JSON of the reward spec constants + checker ABI version)     e.g. reward_v5@3f2a..
store_id    = sha256(item-store manifest.json)                                              (data/store.py)
view_id     = sha256(view manifest: store_id + view params + per-item doctor-mask digests)
tables_id   = sha256(doctor tables: view_id + reward_id + build-library versions)          (today: fingerprint)
policy_id   = base model HF repo + revision SHA + adapter safetensors sha256
env_id      = container image digest (GHCR) + uv.lock sha256
spec_id     = sha256(locked RunSpec)
answers_id  = sha256(answers.parquet)     answers = (key, view_id, policy_id, sampler, seed, text, truncated)
```

### 5.3 Lineage (every arrow is a recorded id)
```
 raw LIDC/RIGA --> store_id --> view_id --+--> tables_id (reward_id) --+
                                          |                            |
 base model rev + SFT spec --> policy_id(sft) --> RL spec --> policy_id(rl) --> answers_id --> score(reward_id) --> eval.json
                                                     ^                              ^
                                        code_id, checker_id, env_id, spec_id  (stamped on every node)
```
Rule: **generation and scoring are separate artefacts.** `answers.parquet` is the expensive GPU product; scores are a
pure function `score(answers, tables, reward_id)` and can be recomputed on any GPU (or the CPU oracle) later. This makes
reward v5 vs v4 comparisons free of new generation, and makes Phase A's "save answers" a structural rule.

### 5.4 Run record (one per run, written by `objdel.runs`, never by a script)
```
runs/<run_id>/
  spec.locked.yaml       RunSpec with every *_id resolved (the "lockfile")
  stamp.json             VersionStamp {code_id, checker_id, reward_id, store_id, view_ids, tables_id,
                                       policy_id_in, env_id, spec_id, host{gpu, driver, cuda}, seed}
  events.jsonl           structured events (5.6), append-only
  metrics.jsonl          per-step numbers (section 8)
  answers/<split>.parquet
  eval/<split>.json      metrics + CIs + {stamp, answers_id}
  ckpt/step-N/           adapter + optimizer + dataloader position + RNG + stamp
```
The HF model repo commit message and model card carry `stamp.json`. A missing id stops the run at start (Architecture.md
§8.2 already says this; today nothing enforces it).

### 5.5 Typed run spec (contract)
```python
# objdel/core/spec.py  (pydantic v2, extra="forbid", frozen)
class ModelRef(BaseModel):   base: str; revision: str; adapter: ArtifactRef | None; lora: LoraSpec | None
class DataRef(BaseModel):    store_id: str; views: list[ViewRef]; split: Literal["train","val","test"]
class RewardRef(BaseModel):  name: Literal["reward_v4","reward_v5"]; reward_id: str | None   # None => resolve and lock
class SamplerSpec(BaseModel): n: int = 8; temperature: float = 1.0; max_tokens: int; seed: int
class AlgoSpec(BaseModel):   loss: Literal["dr_grpo","cispo","gspo"]; beta: float; lr: float; group: int = 8
                             is_correction: Literal["sequence_mask","token_trunc","none"]; is_cap: float
                             max_staleness: int = 0                                   # 0 = on-policy
class ComputeSpec(BaseModel): gpus: int; gpu_type: str; region: Literal["US","CA"]; cap_hours: float; cap_usd: float
class Stage(BaseModel):      kind: Literal["sft","rl","generate","score","eval"]; args: dict   # validated per kind
class RunSpec(BaseModel):    run_id: str; model: ModelRef; data: DataRef; reward: RewardRef
                             sampler: SamplerSpec; algo: AlgoSpec | None; compute: ComputeSpec
                             engine: Literal["trl","verl","nemo_rl"] = "trl"; stages: list[Stage]
```
- `objdel.runs.lock(spec) -> LockedSpec` resolves every id from local manifests (no GPU), refuses a dirty tree, and
  checks invariants (e.g. `sampler.n == algo.group == reward group`, today checked inside `grpo_config`
  [code: REPO/src/objdel/train/configs.py grpo_config]).
- `preflight(locked)` runs the existing data/checker/code gates on exactly that locked spec; the launch requires the
  preflight stamp for the same `spec_id`.
- Engine adapters map the spec to their native config: TRL `GRPOConfig` now; verl Hydra overrides later. The spec is
  the stable surface; framework configs are generated, never hand-edited.

### 5.6 Event schema (logs)
```
{"ts":"2026-11-08T12:00:00.123Z","run_id":"rl2","step":41,"component":"verifier","level":"INFO",
 "event":"score_batch","n":256,"ms":38.2,"reward_id":"reward_v5@3f2a","policy_version":41}
```
Bounded labels: `component`, `level`, `event`. High-cardinality values (`key`, text) stay in the body or in Parquet.

---

## 6. Components and sequences

### 6.1 Component view (stage 1, one GPU)
```
                    owner (Mac)                                   rented GPU host (US/CA)
 +-------------------------------------------+      +--------------------------------------------------------+
 | specs/rl2.yaml --> runs.lock --> preflight |      |  runs.launch(locked spec)                              |
 |        (typed)     (ids)      (Mac gates)  |--+-->|    |                                                   |
 | infra: host policy, cost cap, bundle       |  |   |    v                                                   |
 +-------------------------------------------+  |   |  rl.Loop  (algorithm: advantage, loss cfg, guards)      |
                                                |   |    |            |                   |                   |
           bundle (tar + SHA-256 manifest)  ----+   |    v            v                   v                   |
                                                    | policy.Sampler  verify.Verifier     rl.adapters.trl     |
                                                    | (vLLM colocate) (in process,        (TRL 1.17.1 +       |
                                                    |                  low-prio stream)    FastGRPO patches)  |
                                                    |    |            |  DoctorStore (GPU, byte budget, LRU) |
                                                    |    +------------+----------- geom._gpucheck (CUDA) -----+
                                                    |  runs.record: stamp.json, events.jsonl, metrics.jsonl  |
                                                    |  answers/*.parquet, ckpt/ --> pull.sh, HF private repo |
                                                    +--------------------------------------------------------+
```

### 6.2 Sequence: one stage-1 RL step (in turns on one GPU)
```
rl.Loop      Sampler(vLLM)      Verifier        TRL adapter        runs.record
   |  batch(32 keys, v=41)          |                |                  |
   |----------->|                   |                |                  |
   |            | 256 answers + logp|                |                  |
   |<-----------|                   |                |                  |
   |  score(256 items, reward_id)   |                |                  |
   |------------------------------->| one engine call|                  |
   |<-------------------------------| terms, reward, stamp              |
   |  advantages (RLOO, in rl)      |                |                  |
   |  update(batch, adv, logp_old)  |                |                  |
   |------------------------------------------------>| fwd/bwd/opt      |
   |<------------------------------------------------| loss, IS stats   |
   |  push LoRA delta (v=42)        |                |                  |
   |----------->|                   |                |                  |
   |  metrics(step 41: phase times, IS share at 0, staleness=0, ...)    |
   |------------------------------------------------------------------->|
```

### 6.3 Sequence: stage 2, two GPUs, async (verl first choice)
```
Controller      Rollouter(GPU0,vLLM)   Verifier actor(GPU1, 0.1 GPU)   Buffer        Trainer(GPU1)
    |  submit prompts while  |                    |                    |               |
    |  floor((Nr-1)/B) <= v+eta (AReaL rule)      |                    |               |
    |----------------------->| sample done (v_s)  |                    |               |
    |                        |------------------->| score (streaming)  |               |
    |                        |                    |------------------->| put(sample,   |
    |                        |                    |                    |  v_s, reward) |
    |                        |                    |                    |<--- take B ---|
    |                        |                    |                    | drop if        |
    |                        |                    |                    | v_t - v_s>eta  |
    |                        |                    |                    |               | update -> v_t+1
    |                        |<------------- LoRA delta (NCCL / CUDA IPC), version v_t+1 -----------|
    |  log staleness/mean,max; engine mismatch at age 0; dropped-stale count                        |
```
Consistency model: bounded staleness eta (default 1). Each sample is consumed at most once. On resume, the buffer is
dropped and regenerated (at-least-once generation, exactly-once consumption per optimizer step). A trajectory that
spans two versions (partial rollout) is tagged with its oldest version.

### 6.4 Sequence: run lifecycle and provenance
```
Owner        runs.lock        preflight(Mac)     infra/host policy     GPU host           HF / Mac
  | spec.yaml   |                  |                    |                   |                  |
  |------------>| resolve ids,     |                    |                   |                  |
  |             | refuse dirty tree|                    |                   |                  |
  |             |----------------->| data coverage, leakage, checker goldens, TRL names, dry-run
  |<---------------- cost, GPU-h, spec_id, preflight stamp                    |                  |
  | approve --------------------------------------------->| US/CA, cap $ / h  |                  |
  |                                                       |------------------>| env_check (GPU job
  |                                                       |                   |  tests, image digest)
  |                                                       |                   | run stages       |
  |                                                       |                   | ckpt every N --->| push (retry)
  |                                                       |                   | stamp + record ->| pull.sh (checksum)
  |                                                       |<-- teardown verified (API says gone) |
```

---

## 7. API contracts

### 7.1 Verifier, in process (the primary interface)
```python
# objdel.verify (no torch, no TRL)
@dataclass(frozen=True)
class ScoreItem:
    id: str                    # caller's id, echoed back
    key: str                   # "<dataset>/<split>/<item>::<view>", resolved against the loaded view
    text: str
    truncated: bool
    policy_version: int | None = None   # passthrough for async bookkeeping
    expect_W: int | None = None         # optional asserts; the store owns W, H, names
    expect_H: int | None = None

@dataclass(frozen=True)
class ScoreResult:
    id: str
    reward: float              # per-answer scalar, a pure function of (reward_id, tables, text)
    terms: dict[str, float]    # format_ok, legal, iou_tc, iou_min, f5, g_area, location, topo_any, topo_freq
    defects: dict[str, int]
    error: str | None          # parse code, or None

class Verifier:
    def __init__(self, reward: str, tables: DoctorStore): ...
    stamp: VersionStamp        # reward_id, tables_id, checker_id
    def score(self, items: Sequence[ScoreItem]) -> list[ScoreResult]: ...      # one engine call
    def score_masks(self, masks, keys) -> list[ScoreResult]: ...                 # baselines, gates
    def spec(self) -> RewardSpec: ...

class DoctorStore:
    @classmethod
    def open(cls, view_root, reward, device=0, budget_bytes=...) -> "DoctorStore": ...   # verifies tables_id
    def ensure(self, keys) -> None    # upload missing tables (one batched H2D), evict LRU over budget
    resident_bytes: int
```
Rules: (1) `reward` is per answer; group terms (RLOO baseline, dead-group drop) live in `rl`. (2) The store, not the
caller, knows W, H and the object names of a key: the client cannot send a mismatched schema. (3) The result is exact
and deterministic: same inputs give the same bits on any GPU and the CPU build.

### 7.2 Verifier, out of process (only when a second process needs it: stage 2 or an external user)
```
POST /v1/score                                  (JSON, max 4,096 items, max 8 MB body)
{ "reward": "reward_v5", "reward_id": "3f2a91c0d4e1",        # reward_id optional; if given, must match
  "items": [ {"id":"a_000123","key":"lidc/train/0042::shift","text":"{...}","truncated":false,
              "policy_version":41} ] }

200
{ "stamp": {"reward_id":"3f2a91c0d4e1","tables_id":"9b0c..","checker_id":"cuda-sm120-7d1e..","code_id":"be52b1c"},
  "results": [ {"id":"a_000123","reward":0.73,
                "terms":{"format_ok":1,"legal":1,"iou_tc":0.71,"f5":0.90,"location":0.0,"topo_any":1},
                "defects":{}, "error":null} ] }

errors (whole batch):
{ "error": {"code":"spec_mismatch|unknown_key|batch_too_large|device_error|not_ready",
            "message":"...", "request_id":"r_..", "retryable": false} }
GET /v1/spec   -> reward spec + stamp        GET /livez (process)   GET /readyz (tables loaded, device ok)
GET /metrics   -> Prometheus text (requests, items, ms p50/p99, resident_bytes, errors by code)
```
- Floats: terms are float32 on device, sent as float64 shortest repr; a round-trip test proves bit equality.
- Idempotency: scoring is a pure read, so a retry is safe without a key. Only `device_error` is retried, once, by the
  client; then the run stops (fail closed, section 8).
- Versioning: additive fields only within `/v1`; a new reward is a new `reward` name, never a changed meaning.
- Unknown key is a batch error (a data bug), not a per-item 0 reward.
- Implementation: FastAPI in `verify/http.py` as the `[serve]` extra; a verl `compute_score` wrapper and an
  OpenEnv single-step env wrapper are thin clients over the same in-process `Verifier`.

### 7.3 Trainer adapter and sampler (the swappable seams)
```python
class Sampler(Protocol):                       # vLLM colocate / vLLM server / HF
    def generate(self, prompts: Batch, n: int, params: SamplerSpec) -> Rollouts   # ids, text, logp, truncated
    def load_adapter(self, delta: LoraDelta, version: int) -> None

class TrainerAdapter(Protocol):               # TRL now, verl later (Tinker-like primitives)
    def forward_backward(self, batch: TrainBatch, loss: AlgoSpec) -> StepStats    # IS ratio stats inside
    def optim_step(self) -> int                                                   # new policy version
    def export_adapter(self) -> LoraDelta
    def save(self, path) -> None; def load(self, path) -> None
```
TRL private names stay in one file with a start-up contract test [Architecture.md §6.2]. When verl arrives, the
adapter is replaced and `rl.Loop` either drives it or (more likely) verl's own controller drives our Verifier
through `compute_score`; in that case `rl` keeps only spec mapping, advantage settings and guards.

---

## 8. Failure modes and recovery

Principle: in training, **a wrong reward is worse than a stopped run**. The verifier fails closed. There is no
"stale default reward". Degrade only on the paths that do not touch the gradient (uploads, metrics sink, eval retries).

| # | Failure | Detector | Response / recovery |
|---|---|---|---|
| F1 | Host lost or preempted (Vast) | heartbeat file stale; cost ledger | Checkpoint every N steps (50 today [code: rl_shifted.yaml]) with dataloader position + RNG; resume on a new host from the last pushed checkpoint; buffer dropped. |
| F2 | Verifier device error / OOM | exception from engine | One retry after freeing the pool; then stop the run with the batch dumped to `events.jsonl`. Never return 0. |
| F3 | Doctor tables built under another reward/library | `tables_id` mismatch (today `StaleCacheError`) [code: cache.py docstring] | Refuse at start (preflight), not mid-run. |
| F4 | TRL private API drift on upgrade | start-up contract test | Run refuses to start; TRL stays pinned. |
| F5 | Trainer-vLLM mismatch (IS min = 0 in 98.1 % steps) | per-step IS share at 0 > 0.5 % alarm | Phase B2 fix before stage 2; guard stops a run whose share stays high for K steps. |
| F6 | Async staleness grows (rollouter faster than trainer) | `staleness/max` > eta | Controller throttles submissions (AReaL rule); samples older than eta are dropped and counted. |
| F7 | Runaway answers hit the token cap (38 % of R5 steps) | share at limit > 5 % | Masked today; the cap is set from a true outline-length count (Phase A). |
| F8 | Dead groups (no spread) | share > 25 % | Logged; curriculum or view mix adjusted by a spec change, not mid-run. |
| F9 | Data leakage across splits via views | leakage gate in preflight (group-aware keys) | Refuse launch. |
| F10 | Doctor store exceeds GPU memory (3D) | `resident_bytes` vs budget | LRU eviction with batched re-upload; budget is a spec field; preflight computes the worst case. |
| F11 | HF Hub or network outage at push | push error | Local copy kept; retry with backoff + jitter, 3 attempts; pull.sh is the second path. |
| F12 | Cost overrun / forgotten box | cap_hours, cap_usd, teardown verify | Hard stop; teardown confirmed through the provider API. |
| F13 | CPU vs GPU drift in the engine | golden files on both builds; GPU CI job | Block the merge (CUDA paths) or the launch (preflight env_check). |
| F14 | Silent provenance loss | `stamp.json` validator | A run with a missing id refuses to start; an eval without `answers_id` is not accepted as a result. |
| F15 | Logging fills the host disk | bytes written per hour | Bounded per-step logs; full completion tables sampled (1 group in 8) except in eval. |

---

## 9. Trade-offs (solves / worsens / when to change)

| Decision | Solves | Worsens | Change it when |
|---|---|---|---|
| Modular monolith, verifier library first; HTTP face optional | No network failure modes; one deploy; exact in-process calls | The verifier service is untested until stage 2 | A second process needs scoring (stage 2) or an external user appears (Q1). |
| Two distributions (`objdel-verify`, `objdel`) | verl/NeMo workers install the verifier without torch/TRL | Two version numbers to keep in step | Never merge back; split further only if 3D adds heavy deps. |
| Pydantic RunSpec + generated framework configs | One validated surface; ids locked before money is spent | A mapping layer per engine; Hydra users lose overrides by CLI | verl becomes the only engine: consider emitting Hydra directly from the spec, still validated. |
| Content-addressed ids + run record | Every result reproducible and auditable | Hashing time at lock (seconds); discipline in every stage | Data grows to many GB: hash manifests, not files (already the store's pattern). |
| Separate generate and score artefacts | Re-score under a new reward without GPU generation | Storage for answers (2,200 x 8 x ~400 B = about 7 MB per eval [estimate]) | Never; it is cheap. |
| Fail closed on verifier errors | No silent reward corruption | A transient fault costs a resume | Fault rate makes resumes expensive: fix the fault, do not add defaults. |
| Bounded staleness eta = 1 at stage 2 | Overlap of generation and training; gates off-policy risk | Needs IS correction that already misbehaves at age 0 | Reward curve matches stage 1 with eta = 1: try eta = 2 with CISPO/GSPO. |
| verl as first async engine (NeMo-RL second) | VLM + server-mode rollout + streaming reward + staleness control exist | FastGRPO patches (shared prefix 31.0 -> 23.6 s/step) do not port; re-measure from verl's baseline | verl cannot run Qwen3-VL + LoRA + image prompts at our step budget: try NeMo-RL. |
| JSONL logs + optional trackio, no log pipeline | Zero ops; works offline on a rented box | No live cross-run search | Multi-host runs or >1 concurrent run. |
| Byte-budgeted LRU doctor store | 3D and many views fit one GPU | Re-upload cost on misses | Working set fits: budget = all; set it from preflight numbers. |
| GPU CI on a gated Vast box (not GitHub GPU runners) | ~$0.10 per run; same GPU class as training | Manual approval; not on every push | The repo moves to a plan with GPU runners, or CUDA changes become daily. |
| Hosted compile-only nvcc job on every CUDA change | Catches build breaks for sm_80/90/120 for free | Does not execute kernels | Never removed; it is the cheap first gate. |

---

## 10. Observability, logging and CI

### 10.1 Training SLOs (symptoms that stop or flag a run)
| SLI | SLO | Action |
|---|---|---|
| Step time p50 | <= 18 s stage 1 target (Architecture.md §1.2) | Flag; profile with NVTX + Nsight Systems |
| GPU busy share (NVML sampler, 1 Hz) | >= 80 % | Flag; the per-phase NVTX ranges show the idle phase |
| IS ratio share at 0 | <= 0.5 % | Stop after K steps (guard) once B2 lands |
| Legal share | >= 0.99 (test today 0.994 [measured: final_eval.json]) | Stop on collapse |
| Staleness max (stage 2) | <= eta | Throttle (F6) |
| Cost per run | <= approved $ | Hard stop (F12) |
Cause metrics (CPU prepare time, host syncs, H2D bytes, resident bytes) are logged, never paged.
"Watch the watcher": heartbeat age, metrics lines per minute, log bytes per hour.

### 10.2 Logging
Per process `events.jsonl` + `metrics.jsonl` with the schema in 5.6; pulled with the run; merged by `(run_id, step, ts)`.
No ELK/Loki (YAGNI: one host, < 1 MB/min [estimate]). Optional trackio sink for dashboards (HF library, local-first).

### 10.3 CI
```
on every push (GitHub-hosted, free)                    on CUDA paths or release (gated, owner-approved)
+------------------------------------------+           +-----------------------------------------------+
| lint: ruff + import-linter layers        |           | gpu: Vast US/CA box, existing preflight       |
| test: CPU build, oracles, goldens (fast) |           |  - make cuda (nvcc, sm_120)                   |
| train: CPU torch, coverage >= 90 %       |           |  - test_cuda_checker_gpu + goldens on device  |
| gates: slow K1/K2 exhaustive oracles     |           |  - compute-sanitizer memcheck/initcheck       |
| wheels: cibuildwheel (exists)            |           |  - verifier service smoke (/readyz, /v1/score)|
| NEW nvcc-compile: nvidia/cuda devel      |           |  - 20-step RL smoke on a tiny model (optional)|
|   image, builds sm_80;90;120, no run     |           |  result stamp uploaded as a CI artefact       |
| NEW spec: every specs/*.yaml validates   |           +-----------------------------------------------+
|   and locks against fixture manifests    |
| NEW contract: TRL private names exist    |
+------------------------------------------+
```
Existing CI has lint/pytest on three OSes, a CPU-torch coverage job, slow oracle gates and wheels, but no GPU job;
GPU tests are skipped where `nvidia-smi -L` lists no GPU [code: REPO/.github/workflows/tests.yml;
REPO/tests/test_cuda_checker_gpu.py:27].

---

## 11. Scale evolution

| Rung | Trigger to climb | What changes | Breaking point |
|---|---|---|---|
| 1. One GPU in turns (now) | - | Phase C speed work; verifier in process | Generation > 30 % of step after Phase C, or 3D answers (~1,300 tokens) make the decode tail dominate |
| 2. Two GPUs async, one host | Gate: Phase B2 IS match fixed; same reward curve as rung 1 over 50 steps; FastGRPO gains re-measured on verl | verl fully async (eta = 1), Verifier as a GPU Ray actor with HTTP face, LoRA delta sync | Trainer GPU saturated or one rollout GPU not enough |
| 3. N GPUs (8+) | A 4B model or 3D needs it | More rollout replicas; trainer FSDP2; one Verifier per node (tables per node) | Doctor tables per node exceed the budget; multi-host logs need a pipeline |
| 3D shape kind (orthogonal) | Phase E done | `volume3d` kind in geom + verify; DoctorStore byte budget; 3D goldens + oracle | Answer length; table memory (23 GB per 1,000 objects [estimate]) |

The verifier API does not change across rungs: `task`/kind is implied by the key's view; the reward name selects
the maths. Only placement (in process -> actor -> per-node service) changes.

---

## 12. What I would change in Architecture.md

1. **§0 / intro "ready for production use"**: remove the claim. RL and eval artefacts carry no versions today
   [measured: runs/r5/shared_rl.json, final_eval.json]. Say "designed for" until §1.2 "Reproducible" is enforced.
2. **§1.2 Reproducible**: name the enforcement point: `runs.lock` + `stamp.json` validator; define the ids (5.2 here).
   Add `checker_id`, `env_id`, `policy_id`, `answers_id`, not only "code, reward, data".
3. **§2 layers**: add a bottom `core` layer for shared types (RunSpec, VersionStamp); split `verify.build` (EDT, tau,
   T5) from run-time `verify`; move group advantages from verify to `rl`; move `infra` code out of experiments.
   Name the enforcement tool: import-linter layers contract.
4. **§2.1**: state that `verify` ships as its own wheel (`objdel-verify`) without torch/TRL; that is what makes D2 real.
5. **§4.2 contract**: drop client-sent `W`, `H`, `names` (the store owns them; keep as optional asserts); add
   `reward_id` hash, `policy_version` passthrough, `retryable` in the error envelope, `/livez` + `/readyz` instead of
   one `/healthz`, batch and body caps, the float64 round-trip rule, and "unknown key = batch error".
6. **§4.2 rules**: add "fail closed: no default reward on error" explicitly.
7. **§4.3 shape kinds**: add a device-memory budget and LRU for the doctor store; give the 3D table size estimate.
8. **§6.3**: state that the reward is a per-answer pure function; RLOO baseline, dead-group handling and any group
   term are algorithm settings in `rl`. Check that reward v5's "topology" term is per-answer (v4's topology reward was
   in advantage units and group-centred [code: REPO/src/objdel/reward/spec.py:30 c_topo]).
9. **§6.4 stage gates**: add to the stage-2 gate "B2 trainer-vLLM match fixed" and "FastGRPO gains re-measured on
   verl's baseline"; cite that TRL AsyncGRPO has no image support and runs rollout with no GPU, so the GPU verifier
   must be a service there [web: https://huggingface.co/docs/trl/async_grpo_trainer]. Name NeMo-RL's
   `max_trajectory_age_steps` as the fallback engine's staleness knob.
10. **§7.2 evaluation**: add "generation and scoring are separate artefacts; answers saved with ids".
11. **§8**: define how "GPU busy" is measured (NVML sampler + NVTX), log staleness and engine mismatch separately
    (slime decomposition), and name the metrics sink (JSONL source of truth; trackio optional).
12. **§9 tests**: make the GPU job concrete: free hosted nvcc compile job on every CUDA change, plus a gated Vast GPU
    job (goldens, sanitizer, verifier smoke). GitHub GPU runners need a paid plan and cost about 7x a Vast 5090 hour.
13. **§10 operations**: insert `lock` before preflight; the launch requires a preflight stamp with the same `spec_id`;
    teardown is verified through the provider API.
14. **§11 D8**: keep "one modular program"; add "the HTTP face is built only when a second process needs it" (YAGNI).
15. **Missing section**: archive policy (experiments frozen per phase, moved to `archive/`, excluded from lint/test).

---

## 13. Ten checks any implementation plan must pass

1. **Layer check**: import-linter passes with the order core <- geom <- verify <- data <- policy <- {rl, sft} <- runs
   <- infra, and `objdel-verify` installs and scores in a fresh venv with **no torch, no TRL, no scipy, no cv2**.
2. **Provenance check**: every artefact the plan produces (checkpoint, answers, eval JSON, HF commit) contains a
   complete `stamp.json`; a test deletes one id and the run refuses to start.
3. **Spec check**: every GPU run is one `specs/*.yaml` that validates as `RunSpec` and locks; no argv-list plans; the
   launch fails without a preflight stamp for the same `spec_id`.
4. **Exactness check**: CPU and GPU builds give identical bits on the golden files (195,809 cases per Architecture.md
   §9) and on the verifier's float64 JSON round trip; the GPU job runs this before any change to `cuda/` merges.
5. **Fail-closed check**: a fault injected into the engine (device error, unknown key, tables_id mismatch) stops the
   run with a dumped batch; no code path returns a default reward.
6. **Re-score check**: an eval can be recomputed from saved `answers.parquet` under reward v4 and v5 with no
   generation, and the v4 numbers equal the original eval's numbers bit for bit.
7. **GPU-everywhere check**: the plan lists every remaining CPU step with a reason (PNG decode, text parse, Hungarian
   on <= 32 doctors); doctor-table build (EDT), polygonizer and eval metrics run on the GPU, with the CPU code kept
   only as test oracles.
8. **Swap check**: the TRL adapter can be replaced by a stub adapter in a test without editing `verify`, `data` or
   `policy`; a verl `compute_score` wrapper scores a fixture batch with the same bits as the in-process call.
9. **Async readiness check**: before any 2-GPU spend, the plan shows (a) IS share at 0 <= 0.5 % at stage 1,
   (b) per-sample `policy_version` carried end to end, (c) staleness and engine mismatch logged separately,
   (d) a stage-1 vs stage-2 reward-curve comparison as the gate.
10. **3D readiness check**: the doctor store has a byte budget with LRU and a preflight worst-case number; a
    `volume3d` kind adds goldens and an oracle; the verifier HTTP and in-process contracts do not change.

---

## 14. Score and diagnosis (quality bar, 0-5)

| Dimension | Score | What would raise it |
|---|---|---|
| Fundamentals (consistency, coordination) | 4 | A measured stage-2 staleness/IS curve, not only a rule |
| Building blocks under stress | 4 | Real verl + Qwen3-VL + LoRA + images test; today it is an unverified assumption |
| Requirements | 4 | Owner answers to Q1-Q4 |
| Trade-offs | 5 | - |
| Numbers | 3 | **Weakest.** Doctor-table sizes use worst-case frames; the real frame-size histogram and 3D slice counts are unknown; verl's baseline step time on our task is unknown |
| Failure | 4 | Fault-injection tests that prove F2/F3/F14 |
| Collaboration | 3 | Q1-Q4 are asked but not yet answered; this is a one-pass design |

Quick diagnostic: requirements stated; estimates present but partly worst-case; API concrete; trade-offs complete;
reads (score) vs writes (checkpoints, answers) handled; SPOFs covered (host, verifier, Hub); 10x/100x known per rung.
Coverage sweep: blob storage = HF private repo + local Parquet (no new store needed); ids = content hashes (no
sequencer needed); DNS/LB/CDN/search/sharded counters = not relevant (single host, no users); scheduling = owner-driven
launches (no scheduler needed); logs/SLOs = section 10.

**Current bottleneck**: not architecture but provenance and the trainer-vLLM mismatch. The headline result has no
version stamp, and about half of each step's signal is lost to IS masking. **At the next order of magnitude**
(10x tokens per answer from 3D, or 8+ GPUs), the bottleneck moves to generation and to doctor-table memory: that is
when async with bounded staleness, a per-node verifier actor and a byte-budgeted store stop being optional. This design
is a hypothesis; the verl-on-VLM assumption is the part most likely to break.
