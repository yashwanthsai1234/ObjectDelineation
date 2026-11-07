# Implementation plan: fast RL (R2) on an RTX 5090 (2026-11-06)

Inputs:
- `issues.md` (what is slow and why) and `fix_test.md` (the first fix list);
- research report `docs/research/reports/RL speed stack on RTX 5090.md`, with its notes;
- earlier research report `NVIDIA CUDA performance for VLM training.md`.

Owner decisions:
- RL answers come from **vLLM** (colocate), not HF `generate`;
- the card is an **RTX 5090**;
- no more SFT runs.

Rules:
- every step has its tests;
- old code that a step replaces is deleted;
- `make test` and `make lint` pass before a commit, and the owner runs git;
- nothing is rented before `make preflight` passes and the owner sees the hours and dollars.

```
Stage A  config fixes F1-F3 ............................ DONE on the Mac (tests/test_p2_perf.py)
Stage B  vllm_native image (training + vLLM in one venv): name, TRL 1.17.1, lock → owner builds → allowlist
Stage C  Mac code: C1 trim · C2 checker · C3 scoring · C4 5090 · C5 positions · C6 vision cache · C7 async copies
Stage D  box chain + profile script (GRPO arms, smoke checks)
Stage E  profiling run on the 5090 (~1 h, quote first) ═══► pick settings
Stage F  update step: eager vs compile vs reduce-overhead vs manual CUDA Graphs ═► keep the fastest that passes
Stage G  R2 G-probe with the chosen settings
Stage H  after R2: evaluation fixes, checker kernel stores
```

---

## Stage A: config fixes (done)

| Fix | Change | Test (`tests/test_p2_perf.py`) |
|---|---|---|
| F1 | `g_probe` no longer forces HF: vLLM colocate, `vllm_gpu_memory_utilization` 0.35, sleep mode on | `test_g_probe_generates_with_vllm_colocate`, `test_vllm_kv_budget_holds_one_runaway_group` (24 and 32 GB) |
| F2 | GRPO micro-batch 16 × 16 accumulation (was 8 × 32), still 256 answers per step | `test_grpo_step_keeps_256_answers_in_larger_micro_batches` |
| F3 | `logging_nan_inf_filter` off, non-blocking copies, no FLOP counting (SFT and GRPO) | `test_configs_turn_off_the_per_micro_batch_syncs`, `test_trainers_skip_flop_counting` |

---

## Stage B: the `vllm_native` image

### B0. Name: the `vllm_native` image (owner, 2026-11-06)

The single-venv image (training + vLLM in one `/opt/venv`) is renamed from "phase3" to **`vllm_native`**. The old name
is removed everywhere; nothing is kept as an alias.

| Today | New |
|---|---|
| `docker/Dockerfile.phase3` | `docker/Dockerfile.vllm_native` |
| `docker/locks/p3.txt`, `p3.override` | `docker/locks/vllm_native.txt`, `vllm_native.override` |
| `build_images.sh … phase3` (tag `p3-<hash>`) | `build_images.sh … vllm_native` (tag `vllm_native-<hash>`) |
| `image.yml` input `which: [phase2, phase3]` | `which: [phase2, vllm_native]` |
| `tests/test_p2_image.py` | `tests/test_images.py` (both images: no "p2" in the name) |
| `tests/test_p2_vast.py`, `tests/p2_vast_*.json` | `tests/test_vast.py`, `tests/vast_*.json` (rental tools for every image) |

Files to change:
- `experiments/phase2/docker/build_images.sh` (the lock branch and the build loop; the tag is set by name, no longer
  by `${w/phase/p}`);
- `.github/workflows/image.yml`;
- `experiments/phase2/vast/code_files.txt`;
- the two renamed test files and `experiments/phase2/preflight.py:32` (its test list);
- the "Phase 3 single-venv image" wording in PLAN D15 and `impl/M7.md` (with a dated note);
- this plan. `lab/B5` (history) stays as it is.

Tests:
- `test_image_names`: `build_images.sh` accepts exactly `phase2`, `sampler` and `vllm_native`, and its tags are
  `p2-<hash>`, `sampler-<hash>` and `vllm_native-<hash>`;
- `test_no_phase3_image_name_left`: `phase3` and `p3.` appear in no file under `experiments/phase2/docker`,
  `.github` or `experiments/phase2/vast`;
- the workflow test checks `options: [phase2, vllm_native]`.

The build command becomes:
```
gh workflow run image.yml -f which=vllm_native
```

### B1. TRL 1.17.1 (owner decision, 2026-11-06)

TRL 1.16.0 declares `vllm<=0.31.0`. TRL 1.17.0 (2026-10-29) added vLLM 0.32 support. We move to **TRL 1.17.1**.

Steps (Mac, in this order, because `build_images.sh --lock` copies the versions from the Mac `.venv`):
1. `pyproject.toml`: `trl==1.17.1`. Install it in the Mac `.venv`.
2. Read the TRL 1.17.0 / 1.17.1 release notes for GRPO and SFT changes. Fix our code where they break it, deleting the
   old code paths. Known change: TRL 1.17 refuses `top_p < 1` or `min_p > 0` without `top_k` when importance sampling
   is on. We use `top_p` 1 and `min_p` 0, so we are not affected.
3. `make test` and `make lint` pass.
4. Check that TRL 1.17.1 accepts the `vllm_native` pins: transformers 5.19.0, torch 2.13.0. If not, stop and report
   to the owner.
5. Run `build_images.sh --lock`. This rewrites `train.txt` (TRL 1.17.1 for the phase2 image as well, with a new hash)
   and `vllm_native.txt`.

Tests:
- `tests/test_images.py` (existing: locks committed, workflow manual and pinned, image choice);
- `test_vllm_native_lock_pins_a_trl_that_supports_its_vllm`: reads `locks/vllm_native.txt` and asserts
  `trl==1.17.1` and `vllm==0.32.0`;
- the full suite on TRL 1.17.1 (step 3).

### B2. Build (owner runs it: GitHub)

The image holds only packages; code and data are never baked in. So the build depends only on B0 and B1, and it runs
on GitHub **while Stage C is written**. The Mac has 12 GiB free and the build needs ~35 GB, so it runs on GitHub
Actions.

```
git add pyproject.toml experiments/phase2/docker experiments/phase2/preflight.py experiments/phase2/vast/code_files.txt \
        .github/workflows/image.yml tests/test_images.py tests/test_vast.py tests/vast_charges_v0.json \
        tests/vast_offer_3090.json tests/test_p2_image.py tests/test_p2_vast.py \
        tests/p2_vast_charges_v0.json tests/p2_vast_offer_3090.json
git commit -m "TRL 1.17.1; vllm_native image; image workflow builds phase2 or vllm_native; image and Vast tests renamed"
git push
gh workflow run image.yml -f which=vllm_native
gh run watch
```

### B3. Allowlist

`.vast-allowlist` holds **one** image, and `launch.sh` uses it for every rental (`policy.py:332`). Replace its
`image` line with the printed `image: ghcr.io/…/objdel-gpu:vllm_native-<hash>@sha256:…`, with the layer sizes as for p2
(ISSUES R9-14). The p2 image stays in the registry. Going back is a one-line change.

Tests: the existing allowlist and policy tests (`tests/test_vast.py`).

### B4. Our scripts on a one-venv image

Two places name the second venv by path:
- `select_ckpt.py:55`: `VLLM_PYTHON = "/opt/vllm/bin/python"`;
- `image_check.sh:34`: the `/opt/vllm` import check.

`vllm_native` has only `/opt/venv`, with vLLM inside it.

Change: the sampler Python is `/opt/vllm/bin/python` when it exists, else the running interpreter. `image_check.sh`
checks whichever venv holds vLLM. Both images then run every script. R2 needs neither path (G-probe runs vLLM inside the
trainer), but later exams do.

Tests:
- `test_sampler_python_uses_the_vllm_venv_when_present` and `..._else_the_running_python` (`tmp_path` fake
  `/opt/vllm`);
- `test_image_check_handles_one_and_two_venvs` (reads the script's branches).

Notes for the `vllm_native` image:
- **The checker build.** It installs `cuda-toolkit[nvcc,cccl]` at run time, as `run_box.sh` already does. sm_120 is
  built with `native`.
- **The allocator setting.** vLLM sleep mode does not work with `PYTORCH_ALLOC_CONF=expandable_segments`. Only
  `train_sft` sets it, so G-probe must not set it. Test: `test_g_probe_does_not_set_expandable_segments` (reads the
  G-probe source and the box chain config).

---

## Stage C: Mac code

### C1. Trim each micro-batch to its own length (new issue RL-6)

TRL pads the 256 answers of a step to the **longest** answer, then shuffles and splits them
(`grpo_trainer.py:1541-1545`, checked). It does not trim after the split:

```
generation batch (256 answers), padded to the longest = 1,536 tokens
micro-batch 1: [████░░░░░░░░░░░░░░░░]  real ~300, computed 1,536
micro-batch 2: [██████░░░░░░░░░░░░░░]  real ~450, computed 1,536
...
after C1:      [████]  [██████]  ...    each micro-batch computed to its own longest answer
```

The waste is in the update forward and backward. The reference forward and the old-logprob forward are computed over
the whole generation batch, before the split, so C1 does not trim them.

Change: in `objdel.train`, add one small function `trim_completions(inputs) -> inputs`. It:
1. finds `L = completion_mask.sum(1).max()`;
2. rounds `L` up to the next bucket of `BUCKETS = (128, 256, 384, 512, 768, 1024, 1536)`, 7 shapes, below torch's
   recompile limit of 8 (the buckets also serve Stage F);
3. slices every per-completion-token tensor to `L`: `completion_ids`, `completion_mask`, `old_per_token_logps`,
   `ref_per_token_logps`, `tool_mask`, and the vLLM sampling logprobs. The exact key names are read from the pinned
   TRL's `_generate_and_score_completions` when C1 is written; they are not verified yet;
4. removes the prompt columns that are left-padding in every row.

The G-probe trainer (a small `GRPOTrainer` subclass in `objdel.train`) calls it on the result of `_prepare_inputs`.

Tests (`tests/test_train_trim.py`):
- `test_trim_cuts_to_the_longest_real_answer_rounded_to_a_bucket`, on a synthetic inputs dict;
- `test_trim_keeps_every_masked_token_and_value`: the loss inputs under the mask are unchanged;
- `test_trim_cuts_every_completion_length_tensor`: every tensor whose width is the completion width is cut. A new
  TRL key of that width fails the test, so a TRL upgrade cannot pass silently;
- `test_loss_is_unchanged_by_trim`: a tiny GRPO model (the CPU fixture of `tests/test_train_trl.py`) computes
  `_compute_loss` on the same inputs with and without trim. The loss must be equal within 1e-6.

### C2. CUDA checker host side (issues CK-1 to CK-3)

Change (`cuda/gpucheck.cpp`, `cuda/gpucheck_py.cpp`, `src/objdel/geometry/shape.py`), following the research notes:
- **One grow-only device arena** (`cudaMalloc`), made at the first call. It grows to the next power of two and is
  never freed per call. Delete the per-call `Buf` (`cudaMallocManaged` + `cudaMemset` + sync) and `cudaFree`.
- **Pinned staging buffers** (`cudaMallocHost`), one for input and one for output, kept with the arena. Then one
  `cudaMemcpyAsync` in, the kernels on one stream, one `cudaMemcpyAsync` out, and one `cudaStreamSynchronize`.
- **Output: the even-odd mask only, bit-packed**, shape `(N, H, ceil(W / 32))` uint32, LSB = leftmost pixel. Each warp
  builds one 32-bit word with `__ballot_sync`, so the stores are coalesced. Python unpacks with
  `np.unpackbits(view(uint8), bitorder="little")[..., :W]`. Delete the nonzero raster and the mask memsets.
- **A mutex around the GPU section**, and the GIL released around it (`nb::gil_scoped_release`). The return array is
  created after the GIL is taken back.
- **One CUDA process.** A process that forked after CUDA started cannot use CUDA, so pool workers must never call
  `_gpucheck` (see C3).

Tests:
- Mac:
  - `test_unpack_matches_a_reference_mask` (widths 1, 31, 32, 33, 63, 64, 65, 512), with a fake packed output;
  - `test_to_masks_unchanged_on_fixtures`: the bit-exact CPU fixtures stay green;
  - `bash cuda/docker_check.sh`: builds archs 80, 86, 89, 90 and 120, so sm_120 compiles without a GPU.
- Box (Stage D):
  - `run_box.sh`: 0 mismatches against polycore, and the timing sweep;
  - a new `test_arena_reuse_is_exact` (GPU only, skipped on the Mac): 200 calls with random batch sizes 1–1024 and
    random raster sizes give the same masks as one-at-a-time calls;
  - `compute-sanitizer --tool memcheck` on the checker tests.
- Done when: the fixed cost is ≤ 0.2 ms per call, and the GPU beats the polycore loop at batch ≥ 64.

### C3. One checker call per RL step (issues RL-5, CK-5)

`TRLReward._terms` already scores the whole generation batch (256 answers) in one `score_batch` call. Today that call
sends one group per pool job, and the job calls the checker inside the worker.

Change (`src/objdel/reward/score.py`, `src/objdel/train/rewards.py`):
- the main process parses the 256 answers and calls `to_masks` **once**;
- it sends `(masks, cache key)` jobs to the fork pool;
- the workers run `score_mask` (CPU only).

Delete the per-group checker call in the worker.

Tests:
- `test_reward_calls_the_checker_once_per_reward_call`: a fake `check_polygons` counts its calls (must be 1 for 256
  answers);
- `test_rewards_equal_the_previous_path`: exact equality on the synthetic fixtures, against values recorded from the
  current code before the change;
- `test_pool_jobs_never_touch_cuda`: the worker entry point imports no `_gpucheck` (checked through `sys.modules` in a
  spawned test worker);
- the existing reward tests stay green.

### C4. RTX 5090 support

Change:
- `preflight.DPH_CEILING` comes from `policy.CARDS[gpu]["ceiling"]` ($1.50 for RTX_5090). Delete the constant.
- The offer query already has `cuda_vers>=13.0`, which implies driver ≥ 580 (the CUDA 13 minimum).
- A start-of-rental check script, `experiments/phase2/docker/env_check.sh`, writes `env_check.json`. Each check fails
  loudly:
  - `nvidia-smi` (name, driver, PCIe generation and width, power limit);
  - `torch.cuda.get_device_capability() == (12, 0)` and `sm_120` in `torch.cuda.get_arch_list()`;
  - a bf16 matmul TFLOPS test (`image_check.sh` peak 209.5);
  - a Triton kernel compile (torch 2.13 cu13 no longer bundles `ptxas`);
  - `torch.compile` of a small MLP, forward and backward;
  - a vLLM engine start for Qwen3-VL-2B with CUDA Graphs, `VLLM_USE_FLASHINFER_SAMPLER=0` and the attention backend
    logged;
  - `nvcc --list-gpu-arch` after the toolkit install.

Tests:
- `test_preflight_ceiling_follows_the_card` (RTX_5090 → 1.50, RTX_4090 → 1.00);
- `test_env_check_lists_every_check`: the script names each check, and its JSON keys match the list;
- the existing `test_vast.py` stays green.

### C5. Precomputed position ids (removes syncs; needed for Stage F)

Qwen3-VL computes its 3-D rope positions inside `forward` when `position_ids` is None (`modeling_qwen3_vl.py:1228`).
That code has `.tolist()` and `.item()` calls (`:911-913, :985, :1044`): host syncs in every forward, and graph breaks
for compile.

Change:
- Compute `position_ids` **once per generation batch**, with the model's own `compute_3d_position_ids`, on the trimmed
  `input_ids`.
- Pass them in the update, reference and old-logprob forwards. This is one override in the G-probe trainer subclass
  (C1).

Tests:
- `test_precomputed_positions_equal_the_model_positions`: a tiny Qwen3-VL config (CPU, random weights), the same ids
  with and without precompute: identical `position_ids` and identical logits;
- `test_forward_does_not_recompute_positions`: a spy on `compute_3d_position_ids` is called 0 times in the forward
  when the ids are given.

### C6. Vision features once per image, not once per answer (new issue RL-7, estimate)

Each prompt's image is repeated for its G = 8 answers, so one step holds 256 copies of 32 images. The frozen vision tower
runs on every copy, in every pass: the reference forward, the update forward and, with checkpointing, the recompute.

```
today : 32 images × 8 copies × ~3 passes = ~768 vision forwards per step
cache : 32 images × 1                     =   32 vision forwards per step   (frozen tower: same output)
```

*Estimate*: ~0.7 TFLOP per image → ~180 TFLOP per pass at 256 copies, about 1–2 s per pass. The profiling run
measures the real share (time in `model.visual`).

Change:
- The G-probe trainer subclass wraps `model.visual` (frozen, LoRA-free) with a per-step cache. The key is each image's
  `image_grid_thw` row plus a hash of its `pixel_values` slice.
- The cache returns the stored `(image_embeds, deepstack_embeds)` for an identical image.
- The cache is cleared at every new generation batch, so the memory holds 32 images only.
- The tower runs under `no_grad` (frozen).

Tests:
- `test_vision_cache_returns_identical_features`: a tiny Qwen3-VL config; 2 images × 4 copies give exactly the logits
  that the uncached model gives;
- `test_vision_cache_runs_the_tower_once_per_image`: a counter on the tower's forward equals the number of distinct
  images;
- `test_vision_cache_clears_per_generation_batch`;
- `test_vision_tower_has_no_trainable_parameters` (the cache is valid only for a frozen tower).

### C7. Asynchronous transfers (CUDA streams and double buffering), where data crosses PCIe

Where the bytes move in one RL step (checked in TRL 1.16):

```
generation batch (once per step):  CPU ──.to(device), blocking──► GPU   ids ~4 MB, pixel_values ~1.6 GB fp32
update micro-batches (16 per step): GPU slices of that batch            no PCIe transfer at all
checker (once per step):           CPU ──polygons──► GPU ──masks──► CPU
```

- **Update micro-batches:** they are slices of tensors that are already on the GPU (`grpo_trainer.py:1541-1545`,
  2443-2462). Double buffering cannot save time there.
- **Generation-batch move:** pin the `pixel_values` and copy them with `non_blocking=True` on a side stream. The
  reference forward can start on the token ids while the images arrive; the stream waits on an event before the
  vision tower. After C6, only the 32 distinct images must cross PCIe (~200 MB, not 1.6 GB).
- **Checker (C2):** cut a large batch into 2 chunks with 2 streams and 2 pinned buffer pairs:

```
stream 1: [copy in A][kernels A][copy out A]
stream 2:            [copy in B][kernels B][copy out B]
             ▲ chunk B's copy overlaps chunk A's kernels
```

Tests:
- `test_pixel_values_are_pinned_and_copied_non_blocking` (a spy on `Tensor.to` / `pin_memory` in the subclass);
- `test_side_stream_waits_before_use`: on CPU the stream is a no-op, and the order of the event record and wait is
  checked with a fake stream object;
- checker on the box: `test_two_stream_chunks_equal_one_call` (exact masks), and an nsys timeline where the copies
  overlap the kernels.
- Done when: the profile shows the H2D copy hidden under compute (no idle GPU gap at the start of the update).

---

## Stage D: the box chain and the profile script

### D1. `experiments/phase2/profile.py` (built 2026-11-06; replaces `profile_sft.py`, deleted)

One g_probe process per arm on R1's e2: 4 optimizer steps of 32 prompts × 8 answers, 8 probe items; step 1 (vLLM,
allocator and compile warm-up) is left out of s/step. nvidia-smi is sampled every 500 ms beside each arm.

| Arm | g_probe flags | Micro-batch × accumulation | What it measures |
|---|---|---|---|
| trl | `--plain` | 16 × 16 | TRL 1.17 as shipped: the reference |
| fast | `--trace DIR` | 16 × 16 | `FastGRPO` (C1 trim, C6 vision cache, C9 overlap); torch.profiler on step 3 |
| mb32 | | 32 × 8 | fewer, larger launches |
| no_ckpt | | 16 × 16, no checkpointing | may run out of memory (recorded) |
| compile | `--compile default` | 16 × 16 | regional torch.compile of the text layers (F) |
| graphs | `--compile reduce-overhead` | 16 × 16 | the same with CUDA Graphs (F) |

Per arm, in `profile.json`: exit code, s/step, the split **generate / rewards / old+ref logprobs / update** (from
`FastGRPO`'s `time/*_s` metrics and TRL's `step_time`), mean answer length, minimum importance-sampling ratio,
nvidia-smi peak MiB, mean utilization and power; for the traced arm the step's GPU busy share, kernels, CUDA API
time and host syncs. The torch.profiler callback is shared with train_sft (`callbacks.TraceStep`).

Tests (`tests/test_profile.py`): every arm keeps 256 answers per step; the flags per arm; the split from a
synthetic g_probe log; the trace parse; nvidia-smi parse; the p1 chain.

### C9. Rewards overlap the old and reference logprobs (built 2026-11-06)

TRL 1.17 computes, per generation batch: generate → old logprobs (GPU) → reference logprobs (GPU) → rewards (CPU and
the checker). `FastGRPO._generate` hands the decoded answers to `TRLReward.prefetch` as soon as vLLM returns; a
background thread parses them, runs the one checker call (own CUDA streams, GIL released) and the pool scoring. The
reward call takes that result only when the answers, truncation flags, keys and sizes are identical; otherwise it
scores as before. The checker and the training both stay on the GPU. Tests: `tests/test_reward_prefetch.py`
(equal to the recorded values, one checker call, off the main thread; a mismatched prefetch is ignored) and
`test_fast_grpo_starts_scoring_when_vllm_returns`.

### Fresh-context review (2026-11-06): fixed before the box run

Two read-only reviewers (one on the CUDA checker, one on the trainer) read the code without our conclusions.
- **Fixed:** `FastGRPO._generate` crashed at g_probe's step-0 sampling (no generation batch yet): `rows` is None
  outside a batch, and the prefetch is skipped. Test: `test_fast_grpo_prefetch_waits_for_a_generation_batch`.
- **Fixed:** `profile.py` died on the `trl` arm (TRL as shipped has no `time/*` metrics) and on any malformed arm
  output; the split is now optional, parse errors are recorded, and the next arm runs.
- **Fixed:** the traced arm's profiler slowed 2 of its 3 timed steps: tracing is its own last arm (`trace`), untimed.
- **Fixed (wrong premise):** TRL 1.17 trains on its own processor-expanded ids, not vLLM's (`grpo_trainer.py:1844`),
  so the "TRL dropped the images" check could never fire; it is gone. The real risk is vLLM expanding the one pad to
  another count (sampling from other image features, silently off-policy): the vLLM wrapper now compares each
  returned prompt to TRL's ids and raises.
- **Fixed:** the compile arms would pass Dynamo's recompile limit (8 widths × grad on/off × adapter on/off > 8) and
  run eager with a warning: limits raised in `compile_text_layers`; `TORCH_LOGS=recompiles` in every arm's log.
- **Fixed (checker):** a failed CUDA call (an OOM, say) left the runtime's last error set and failed the next good
  call; every call clears it first. Arena rounding is capped at 64 MiB (one large `check_masks` call kept twice its
  need, never freed). Checker streams are created at the lowest priority, so the trainer's kernels go first.
- **Open, decided on the box:** `__int128` division in device code (nvcc must accept it: the build step shows it);
  occupancy of the 16 KiB shared window; one-thread witness tests on illegal shapes with many vertices (tail
  latency); GIL contention of the prefetch parse with the main thread's launches (the `fast` vs `trl` arms show it).
- **C10 (built 2026-11-06 after P1; tests on CPU, speed measured on the next rental):** with `steps_per_generation = gradient_accumulation_steps` and one iteration, the old
  logprobs equal the update's own (the weights do not change in between); TRL computes them only for the vLLM
  importance-sampling ratio. Moving that ratio into the loss drops one of three forwards over 256 rows. Built only if
  the box's `logprobs` share justifies changing TRL's loss. P1 measured 12.3 s of 44.3 for old + reference
  logprobs, so it was built: `FastGRPO.skip_old_logps` turns TRL's IS flag off for the generation batch (TRL then
  skips the old pass), and the update's forward sets `inputs["importance_sampling_ratio"]` with `fast.is_ratio`
  (TRL's four modes, worked-example test) from its own detached logprobs, plus TRL's sampling metrics per
  micro-batch. With two iterations or misaligned steps, TRL's old pass runs as before.

### C5 built (2026-11-06, after P1)

P1's trace counted 7,233 `.item()` reads per step, matching Qwen3-VL's `get_rope_index` over 768 row computations
(256 rows x 3 passes). `fast.RopeCache` keeps each row's (3, width) positions for the generation batch, keyed by an
exact hash of the row's ids and mask (the update's micro-batches are shuffled subsets), and hands each chunk its
column slice through a forward pre-hook on the Qwen-VL base model. Test: `test_rope_cache_gives_trls_logps_and_
computes_each_row_once` (equal to TRL's untrimmed logprobs; 5 rows x 3 passes compute 5 rows). Speed on the GPU is
measured on the next rental, beside C10.

### D2. The chain config `experiments/phase2/configs/p1_rl.yaml` (`steps:` form, as `r0.yaml`)

Before the chain: `sync_up.sh p1` (code and the LIDC bundle) and `fetch_hf.sh p1 r1-qwen3vl2b` (e2 from the private
HF repo to `{rd}/ckpt`; the token goes over ssh stdin only, as push_hf). No HF generation anywhere: vLLM generates
in RL and the G-probe (owner, 2026-11-06). The vLLM/TRL image check (TRL #6294) is the `fast` arm itself:
`FastGRPO` sends vLLM one image pad per image and raises if TRL drops the images.

1. `~` box facts (lscpu, cgroup CPU limit, nvidia-smi -q, /dev/shm).
2. `env_check.sh {rd}/ckpt` (C4). A failure stops the chain.
3. The checker: `run_box.sh` (build for the card, GPU tests incl. arena reuse and two-stream chunks, 0-mismatch
   eval, timing sweep), then `tests/test_gpu_scoring.py`. A failure stops the chain.
4. `~` compute-sanitizer memcheck on the checker GPU tests, when the toolkit ships it.
5. The arms (D1).

Cap 1.7 h (RTX 5090 at ≤ $0.79/h: ≤ $1.35, inside the owner's $1.50 for this test).

## Stage E: profiling run (box)

**Done 2026-11-06 (P1, $1.253):** results, decisions and the next steps in `p1_results.md`.

1. `make preflight` (C4 ceiling).
2. `offers.sh RTX_5090` (read-only quote, US/CA).
3. The hours (≤ 1.7 h cap) and dollars (quote × cap) stay inside the owner's $1.50 test allocation (2026-11-06).
4. `make preflight CONFIG=experiments/phase2/configs/r1_qwen3vl2b.yaml` (the gate checks R1's data, checker and code), then `launch.sh OFFER p1 --config experiments/phase2/configs/r1_qwen3vl2b.yaml --cap-hours 1.7` and `run.sh p1 experiments/phase2/configs/p1_rl.yaml`.

**Pass:** every smoke step passes.
**Choose:** the arm with the lowest s/step that fits in memory, with no more than 5 % of answers masked by importance
sampling. Write the choice into `grpo_config` with its date and the measured numbers. Its test asserts the chosen
values.

---

## Stage F: torch.compile and CUDA Graphs for the update step (measured, off by default)

Why this stage exists:
- vLLM already uses CUDA Graphs for **generation**;
- the **update** step (forward + backward) is still eager and launch-bound.

What the research says:
- `mode="reduce-overhead"` records the forward and the backward as CUDA Graphs. It records again for **every new
  input shape**, and it skips graphs that contain CPU ops or `.item()`.
- The recompile limit is 8 (torch 2.13/2.14 source).
- No measured result exists for Qwen-VL training with graphs.
- One vendor benchmark shows compile 9 % **slower** on a 5090 for a small model. So we measure; we do not assume.

How we make it possible:

```
shapes per micro-batch:  [prompt 330 + completion L]  → C1 rounds L up to BUCKETS
                          BUCKETS: 128, 256, 384, 512, 768, 1024, 1536 → 7 shapes (≤ recompile limit 8)
host syncs inside layers: C5 removes the rope .tolist(); vision runs eagerly under no_grad (frozen)
compiled region:          each TEXT decoder layer only (regional compilation: 28 identical layers, one compile)
```

| Step | Change (flag in `grpo_config` overrides, default off) | Test (Mac, CPU) | Box measure |
|---|---|---|---|
| F1 compile | `layer.compile(dynamic=False)` on each text decoder layer, after `GRPOTrainer` wraps the PEFT model; vision and LM head stay eager | `test_compiled_layers_match_eager`: tiny Qwen3-VL config, loss and grads equal within 1e-4 (fp32, CPU inductor) | arm g_compile: s/step, compile time, `TORCH_LOGS=recompiles` count |
| F2 graphs | `mode="reduce-overhead"` on the same layers, `torch.compiler.cudagraph_mark_step_begin()` at each micro-batch, the 7 BUCKETS (no limit change) | `test_bucket_list_bounds_recompiles`: the bucket set is ≤ the limit | arm g_graphs: s/step, peak GiB (graph pools), `cudagraph` skip reasons (must be 0 for the layers) |
| F3 (only if F2 is skipped) | find the remaining break with `torch._dynamo.explain` / tlparse, fix it, measure again | — | — |
| F4 manual graphs | `torch.cuda.make_graphed_callables` captures forward **and** backward of the text decoder layers, without torch.compile. One capture per bucket (7), with shared memory pools (`pool=`). Layers 0–2 are captured one by one (the deepstack step between them is eager: boolean-mask scatter); layers 3–27 are captured as **one** block. Static inputs per bucket: `hidden_states`, rotary `cos/sin`, the attention mask, and `position_ids` (C5). Checkpointing is off inside the captured region, and autocast is entered before the capture with `cache_enabled=False`. The LoRA parameters keep their addresses (AdamW updates them in place outside the graph). | `test_graphed_block_matches_eager` (CUDA only, runs on the box: loss and LoRA grads equal eager within bf16 noise for every bucket); `test_capture_inputs_are_static_per_bucket` (Mac: the buffer shapes per bucket) | arm g_manual_graphs: s/step, capture time, peak GiB, GPU busy share |

The arms compare four ways for the update step: **eager**, **compile**, **compile + reduce-overhead**, and **manual
graphs**. We keep the fastest one that passes the rule below. torch.compile is one option; the manual graphs need no
compiler and avoid its 5090 risk.

```
manual graphs, one micro-batch in bucket 512:
CPU : [layers 0-2: 3 graph launches + 2 eager deepstack steps][block 3-27: 1 launch][head]
GPU :   [■■■][■][■■■][■][■■■][■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■][■■■]
         ~1,100 forward orders → ~10 launches; backward the same
```

Adopt a flag only when:
- it is **≥ 10 % faster** in s/step on the 5090;
- the step-1 loss equals eager within bf16 noise (|Δ| ≤ 1e-3);
- peak memory stays ≤ 28 GiB.

Otherwise the code is deleted, not left behind a flag.

---

## Stage G: R2 G-probe

- `g_probe` on R1's e2, with the Stage E/F settings, through the existing chain (`train_sft` is not run).
- Projection: the measured s/step × 20 steps + 3 probe draws.
- Before launch: preflight, then the quote shown to the owner.
- Tests: the existing G-probe tests; `test_g_probe_settings_are_the_chosen_ones` (Stage E).

---

## Stage H: after R2 (evaluation and kernel work)

| Fix | Issue (`issues.md`) | Test / check |
|---|---|---|
| One vLLM engine start per checkpoint in `select_ckpt` (all views in one child) | EV-1 | prompt-id prefix per view; outputs equal to the per-view path |
| HF check: larger batch and a static KV cache, tested on 5.19 | EV-2 | metric equality on 32 items; timing |
| A/B of per-request seeds in vLLM | EV-3 | same prompts with and without seeds: throughput |
| Scoring overlaps the next branch's sampling; pool size from `nproc`; sampler timeout ~3× the measured time | EV-4, EV-5 | timing log per phase |
| Checker kernels: row-tile coalesced stores, edges in shared memory, 32-bit math (bit-exact) | CK-4 | `run_box.sh` 0 mismatches; Nsight Compute sectors per request |
| Flash-attention-2 for the vision tower (image change); after C6, the gain is small | RL-4 | loss parity |
| Trim the reference and old-logprob forwards per chunk: **done by C1** (`FastForward` covers every pass) | RL-6 | logprob equality |

## Coverage: every issue and its stage

| Issue | Stage |
|---|---|
| RL-1 out of memory with 256 sequences in HF generate, RL-2 slow HF generation | A (F1: vLLM) + B (`vllm_native` image) |
| RL-3 small micro-batches | A (F2) + E (8 / 16 / 32 measured) |
| RL-4 syncs | A (F3) + C5 (rope positions) |
| RL-5 small checker calls from 4 processes | C3 |
| RL-6 padding to the longest answer | C1 (update); H (ref / old forwards) |
| RL-7 vision tower on every copy of an image | C6 |
| Transfers | C7 |
| Rewards wait for the logprob passes | C9 |
| CK-1, CK-2, CK-3 checker memory | C2 |
| CK-4 checker store pattern | H |
| CK-5 batch size | C3 |
| Launch-bound update step (CUDA Graphs) | C1 buckets + C5 + F (compile, reduce-overhead, manual graphs) |
| RTX 5090 | C4 + D (env checks) |
| vLLM + TRL version range, Qwen3-VL support | B1 + D2 (`fast` arm: one pad per image, loud image-drop check) |
| EV-1 to EV-6 evaluation | H (EV-6 is solved by F1) |
| SFT-1 to SFT-4 | not planned (no SFT run); train_sft keeps `profile:` (`callbacks.TraceStep`) |

## What is not in this plan

- **bitsandbytes, QLoRA, FP8:** memory is not our limit, and they are slower or have no proven gain.
- **CUDA streams:** they do not reduce the CPU cost of writing orders.
- **Whole-step manual CUDA Graph capture:** dynamic shapes, the optimizer and the data path make it fragile. The
  per-layer graphs in Stage F cover most of the launch cost.
- **SFT speed work:** there are no more SFT runs. `profile.py` keeps the SFT arms for later.
