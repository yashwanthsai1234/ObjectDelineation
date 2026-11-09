# Phase C: GPU speed (Implementation doc)

Date 2026-11-08. Author: lead engineer, Phase C. REPO = `$REPO`.
SP = `(planning scratch, not kept)`.
Test-writer contract: `docs/implementation2.0/specs/C.spec.md`.

Evidence tags: [measured: file] from a repo log or a run I did now; [code: file:line] read in code now; [paper: arXiv id];
[web: url]; [estimate] my arithmetic or judgment (the inputs are named); [unverified] not checkable on the Mac.
TRL, transformers and torch paths are in `REPO/.venv/lib/python3.*/site-packages/` (trl 1.17.1, transformers 5.20.1,
torch 2.14.1, peft 0.21.2 [code: REPO/pyproject.toml:19-25]). vLLM 0.32.0 is on the box image only
[code: REPO/experiments/phase2/docker/Dockerfile.vllm_native:5; measured: REPO/runs/r5/r5_shared.log:39].

## Skills used

| Skill | What it changed in this phase |
|---|---|
| nvidia-skills:nvidia-skill-finder | NVIDIA catalog route checked through the K3 note (catalog saved at `docs/implementation2.0/notes/doc/nv_catalog.json`). Only `compute-sanitizer-diagnose-cuda-bugs` fits Phase C, and only if the K3 fallback kernel is built (task C1 step 6). Nothing installed (owner approval needed). Nsight Systems + NVTX stay the profiling tool (C0). |
| engineering:code-review | Each task has a correctness rule before speed: "same bits or a stated statistical gate". It found two wrong claims in the notes: the `deepcopy` cost (C1 note S20: TRL copies 32 images, not 256, see §4 C1) and "keep vLLM awake = -1.5 s" (R4 measured -0.05 s). |
| engineering:debug | C0 is reproduce-first: one nsys trace and one sync count before any change. Every task has a bisect order if a gate fails (§7). |
| engineering:system-design | The §3 box diagram, the GPU/CPU placement table with a reason for every CPU item, and the load numbers (bytes per step, syncs per step). |
| engineering:architecture | ADR-shaped decisions inside tasks: pixel backend order (HF on CUDA, then plain torch, then K3 kernel), sleep vs awake, one `load_weights` call vs vLLM LoRA serving, `n=8` request shape. |
| engineering:documentation | Most useful first (goal, gates), exact paths, a runbook-style GPU session in §8. |
| huggingface-skills:huggingface-best | No model change in Phase C (the model stays Qwen3-VL-2B; M1-M3 notes own model choice). Leaderboard workflow not run. |
| engineering:testing-strategy | Each task names its Mac unit test, its recorded-batch equivalence test and its box test; equality tests are bit tests (`torch.equal`), not tolerances, except where the note says why. |
| mattpocock-skills:codebase-design | New code is three deep modules with small interfaces (`gpu_pixels.patchify`, `rope_positions`, `vllm_sync.push_weights`) at the seams E1 already plans (`rl/caches.py`, `rl/sync.py`). No seam with one adapter is added (the K3 kernel is a fallback behind `patchify`, not a second public interface). |
| huggingface-skills:trl-training | GRPO colocate facts used: `n=1` per request in colocate, `generation_kwargs.update` override point, `sync_weights` loop, `beta == 0` skips the reference pass, dr_grpo normaliser. |

Web research done for this phase (beyond the notes): vLLM V1 parallel sampling splits an `n>1` request into `n` child
requests with `n=1`, and prefix caching removes the repeated prefill [web: https://github.com/vllm-project/vllm/pull/10980;
web: https://github.com/vllm-project/vllm/issues/39647]. vLLM model `load_weights` takes an iterable of `(name, tensor)`, and
`AutoWeightsLoader` builds its parameter maps once per call [web: https://github.com/vllm-project/vllm/blob/main/vllm/model_executor/models/utils.py;
web: https://github.com/vllm-project/vllm/blob/7cc302dd/vllm/model_executor/models/qwen3_vl.py] (main branch, not the
0.32.0 tag: [unverified] for 0.32.0, checked on the box by the C5 parameter-equality test).

---

## 1. Goal and why

**Goal.** Cut the R5 step (23.43 s) by removing host work, host waits and wasted transfers, with the same numbers.
No task in C1-C6 may change one sampled token or one gradient bit. C7 and C8 change the objective; they are owner
decisions and run inside Phase B.

**Why (evidence).**

| Fact | Tag |
|---|---|
| R5 step 23.43 s; generate 11.11 s; batch 15.76 s; weight sync 1.27 s, wake 0.11 s, sleep 0.29 s (steps 3-15) | [measured: REPO/runs/r5/shared_rl.json, recomputed now] |
| GPU idle about 40 % of the step; CPU phase 3.0-3.2 s; about 1,250 host syncs per step | [measured/estimate: BRIEF §2, docs/implementation2.0/notes/B_gpu_audit.md]; my own tally of trainer-side syncs is 900-1,100 [estimate: docs/implementation2.0/notes/C1.md §2.3] |
| TRL runs the image processor twice per step; each pass builds 256 x 1,024 x 1,536 fp32 = 1.61 GB on the CPU; pass 2 copies it pageable and blocking | [code: fast.py:458-475; trl grpo_trainer.py:2645; transformers trainer.py:2345-2356 (C1 note)]; size = 256 x 1,024 x 1,536 x 4 B [estimate: arithmetic] |
| One processor pass costs 0.75-1.31 s on the Mac | [measured: REPO/docs/research/notes/cpu_gpu_trace.md row 2.3, cited by docs/implementation2.0/kernels/K3/DESIGN.md §1] |
| Rope positions: HF `get_rope_index` runs as a Python loop per row on host copies | [code: fast.py:336-348; modeling_qwen3_vl.py:925] |
| Weight push: TRL calls vLLM `load_weights` once per parameter (about 700 calls), all weights, every step | [code: trl/generation/vllm_generation.py:520-522; fast.py:696-717] |
| vLLM awake vs asleep: 29.89 vs 29.94 s per step (-0.05 s); awake sync still 1.01 s | [measured: REPO/docs/research/notes/r4_speed_2026-11-08.md] |
| Decode: generate = 5.61 s + 4.3 ms x (longest answer in the step), r = 0.92 | [measured fit: docs/implementation2.0/notes/B_gpu_audit.md:40, :84] |
| 5 of steps 3-15 hit the 1,536 limit; mean longest answer 1,271 tokens; mean answer 131 tokens | [measured: runs/r5/shared_rl.json, recomputed now] |
| Reference pass (beta = 1e-3) costs 1.4-1.7 s per step | [estimate: two-point fit, docs/implementation2.0/notes/B_gpu_audit.md:90] |
| 32 x 8 micro-batches were slower than 16 x 16 on older code (59.1 vs 44.3 s/step) | [measured: perf/p1_results.md via configs.py:184 comment] |

The owner rule "GPU everywhere; CPU only where unavoidable" [BRIEF §3] makes the CPU pixel build and the CPU rope loop
defects, not only slow parts.

**Corrections to Direction.md Phase C (verified here).**
1. Direction C1 says "keep vLLM awake, about 1.5 s". R4 measured -0.05 s for awake alone [measured: r4 note],
   because the awake sync still pushed every weight (1.01 s). The cost is the per-parameter push loop. Awake matters
   as the precondition of a delta push: TRL sleeps at level 2 (vllm_generation.py:377, 772), which discards the
   weights [web: sleep_mode docs], so only a full push is valid after a level-2 wake. C5 pairs "delta push" with
   "awake", and falls back to "sleep + one-call full push" if memory does not allow awake.
2. Direction C3 (32 x 8, -1.3 to -2.5 s) conflicts with the only measurement (slower) [measured: configs.py:184]. It
   moves to C9 as a measured arm after C3 frees memory.
3. Direction C6 (outline limit, 0.3-4 s): with the measured decode fit, a 1,200 limit saves 178 tokens per step on the
   R5 steps (0.77 s) and 81 tokens per step over the 770-step full run (0.35 s) [estimate: 4.3 ms x mean(max_len -
   min(max_len, 1200)); inputs measured from runs/r5/shared_rl.json and runs/rl/train_rl/rl.json]. Not 3 s.
4. Direction's Gate C (18 s with C1-C5) is the optimistic end of my estimate. §6 sets the must-pass gate from the
   estimate's midpoint and keeps 18 s for C1-C8.

## 2. Research basis

| Source | Use in Phase C |
|---|---|
| TRL 1.17.1 `VLLMGeneration.sync_weights` [code: vllm_generation.py:492-528] | The loop C5 replaces; `reset_prefix_cache()` after the push is kept. |
| vLLM `AutoWeightsLoader.load_weights(Iterable)` [web: vllm utils.py, qwen3_vl.py, main branch] | One call with a generator replaces about 700 calls (C5). |
| vLLM sleep mode docs: level 2 frees weights and KV; wake weights, load, wake KV [web: https://docs.vllm.ai/en/latest/features/sleep_mode/] | Sleep level 2 forces a full push each step unless the weights are preserved. C5 keeps TRL's wake order. |
| vLLM V1 parallel sampling [web: https://github.com/vllm-project/vllm/pull/10980] | `n=8` is split into 8 child requests; prefix cache already dedupes prefill. So C6 can save only host-side request work. |
| vLLM automatic prefix caching (hash includes LoRA id, mm hash) [web: https://docs.vllm.ai/en/stable/design/v1/prefix_caching.html] | Prefix cache must be reset after each push (kept). |
| transformers 5.20.1 Qwen3-VL `mm_encoder_outputs` [code: modeling_qwen3_vl.py:1155-1190] | Public forward argument: precomputed vision outputs, so the trainer never needs `pixel_values` after one encode (C3). |
| HF fast image processors accept `device`; torchvision backend [web: https://huggingface.co/docs/transformers/main_classes/image_processor] | First pixel backend for C1. |
| K3 CPU reference == HF processor, 41 images, 0 of 64.5 M values differ [measured: docs/implementation2.0/kernels/K3/build/test_hf_cpu.log] | Proves the fused formula `RN(RN(x - m') / s')` and the patch order; the "digit image" index test is reused in C1 tests. |
| K3 rope template + extrapolation == HF `get_rope_index`, Lp 325, 16 rows x 1,868 [measured: K3/build/test_hf_cpu.log] | C2 closed form. |
| Nsight Systems + NVTX practice [web: A2 note F5, https://dev-discuss.pytorch.org/t/using-nsight-systems-to-profile-gpu-workload/59] | C0 and all speed gates. |
| `torch.cuda.set_sync_debug_mode("warn")` [code: torch 2.14.1 API present, checked now] | Library way to count synchronizing ops per step (C0 gate), no custom hooks. |
| Dr. GRPO normaliser [paper: arXiv 2503.20783]; TRL dr_grpo = sum / (B x max_completion_length) [code: grpo_trainer.py:3286] | Micro-batch shape does not change the update (C9); the outline limit rescales the loss (C7). |
| Training-inference mismatch, TIS/MIS, FP16 [paper: arXiv 2510.26788] | Boundary: Phase B2 owns numerics. Phase C must work under `sequence_mask` and `token_truncate`. |
| APRIL partial rollouts [paper: arXiv 2509.18521] | Rejected for Phase C: off-policy tokens inside an answer. |
| CCE [paper: arXiv 2411.09009], Liger [paper: arXiv 2410.10989] | Not adopted: TRL's fused head already avoids full logits (C3 note F4); loss kernels conflict with packed rows. |
| Speculative decoding at large batch [paper: arXiv 2601.11580] | Not adopted now (C2 note F8): batch is large, outputs are novel numbers. |
| XGrammar [paper: arXiv 2411.15100] | Deferred: a grammar bounds runaways but changes the legal-rate signal. |

## 3. Design

### 3.1 Box diagram (one RTX 5090; `==>` PCIe, `-->` on device)

```
 HOST (CPU)                                   | GPU (RTX 5090, 32 GB)
                                              |
 +-----------------------------+              |
 | OpenOnce / ViewBank         |  32 x uint8  |   +-------------------------------+
 | 32 distinct images per step |==25 MB pinned=>  | C1 patchify (HF on CUDA /     |
 | (PNG decode: CPU, see 3.3)  |   (async)    |   |  torch / K3), 32 images only  |
 +-----------------------------+              |   +---------------+---------------+
                                              |                   | 201 MB fp32 patches (32 images)
 +-----------------------------+              |                   v
 | C2 prompt-ids memo          |              |   +-------------------------------+
 | (text, image size) -> ids   |              |   | C3 vision tower, once per     |
 | once per prompt layout      |              |   | generation batch, no grad     |
 +--------------+--------------+              |   | -> mm_encoder_outputs[32]     |
                | prompt ids (lists)          |   +---------------+---------------+
                v                             |                   |
 +-----------------------------+              |   +---------------v---------------+
 | vLLM host loop (scheduler,  |--requests--> |   | vLLM decode (CUDA graphs)     |
 | sampler) C6: 32 x n=8       |<--token ids--|   | awake (C5 mode A) or asleep   |
 +--------------+--------------+              |   +---------------^---------------+
                | completion ids, lengths     |                   | C5 push: 1 load_weights call,
                v                             |                   | 196 LoRA-touched tensors
 +-----------------------------+              |   +---------------+---------------+
 | C4 RowMeta (numpy)          |--plan------> |   | trainer forward / backward    |
 | prompt_len, completion_len, |  (pinned     |   | C2 rope_positions (torch)     |
 | image_id, group, order,     |   int64      |   | C3 encoder outputs by row id  |
 | windows, packing indices    |   staging)   |   | C4 loss metrics on device     |
 +-----------------------------+              |   +---------------+---------------+
                                              |                   |
 +-----------------------------+              |   +---------------v---------------+
 | answer text parse (ms)      |<--text------ |   | verifier (CUDA checker, Phase |
 | log: ONE device read / step |<--7+N floats-|   | D): rewards on device         |
 +-----------------------------+              |   +-------------------------------+
```

### 3.2 Data flow of one step (after Phase C)

1. OpenOnce hands 32 PIL images (or, after T1, 32 uint8 planes from the view bank) for 256 rows (8 per image).
2. C2: `_tokenize_prompts` returns prompt ids from a memo keyed by (prompt text, image size). No pixels are built.
3. C5: if the step changed, push the 196 LoRA-touched weights in one `load_weights` call (mode A, vLLM awake; full
   push in one call in mode B, sleep on; full push also at step 0 and
   after resume), then `reset_prefix_cache()`.
4. vLLM generates 256 answers (C6 arm: 32 requests with `n=8`). Completion ids and lengths come back as host lists.
5. C4: `RowMeta` is built from those host lists: prompt length, completion length (0 when truncated and masked),
   image id, group id. Sort order, chunk windows, packed-row indices come from `RowMeta`, not from device reads.
6. TRL pass 2 calls the processor: C1 `ImageMemo` returns CUDA `pixel_values` (built from 32 images by one gather), so
   TRL's `.to(device)` is a no-op. After C3, the batch carries a per-row image id instead of `pixel_values`.
7. C3: the vision tower encodes the 32 images once (no grad); every forward chunk gets `mm_encoder_outputs` gathered
   by image id. C2: positions come from one rope template plus a closed form, on the device.
8. Reference pass (absent when beta = 0, C8), update micro-batches, backward, optimizer. C4: loss metrics stay on the
   device and are read once in `FastGRPO.log`.
9. Rewards from the CUDA checker (Phase D owns it) overlap the forwards (existing `prefetch_rewards`, fast.py:629).

### 3.3 Placement: GPU vs CPU (GPU everywhere; each CPU item justified)

| Work | Today | After C | Why CPU, when CPU |
|---|---|---|---|
| PNG decode, 32 per step | CPU | CPU (or zero after T1 view bank preload) | No GPU PNG decoder: nvImageCodec decodes PNG on the CPU only [web: github.com/NVIDIA/nvImageCodec, T1 note]. 32 images x 3 ms [measured-Mac: T1 F5]. |
| H2D of images | none (fp32 1.61 GB pageable later) | 25.2 MB pinned uint8, async | Unavoidable transfer; 64 x smaller [estimate: 32 x 512 x 512 x 3 B]. |
| Rescale, normalize, patchify | CPU, twice, 1.61 GB fp32 | GPU (C1) | - |
| Chat template, tokenization | CPU, 512 calls per step | CPU, once per prompt layout (C2 memo) | Text work; the tokenizer is a CPU library. One call per layout per run. |
| Rope positions | CPU Python loop + 9.8 MB H2D | GPU closed form (C2) | Template: one HF call per prompt layout, CPU, once per run. |
| Vision encode | GPU, keyed by a pixel hash | GPU, once per batch by image id (C3) | - |
| Row metadata (lengths, order, windows, packing plan) | read back from the device (`.tolist()`) | numpy on the host (C4) | The host already holds these numbers (vLLM returns lists). Reading them from the device is the defect. |
| Packed-row indices | Python lists -> pageable `torch.tensor` | numpy into one pinned staging tensor, async H2D (C4) | Index plan is host knowledge; one small async copy per chunk. |
| vLLM scheduler and sampler loop | CPU | CPU | Inside vLLM; per-decode-step host work is inherent. |
| Weight push orchestration | ~700 Python calls | 1 call; tensors stay on the device (C5) | Python drives the loader; no tensor crosses PCIe. |
| Answer text parse | CPU, ms, overlapped | CPU | vLLM delivers token ids to the host; a GPU parser saves about 0 s [SP: B_gpu_audit §3.1, SD1 §6.1]. |
| Loss metrics | 7 `.item()` per micro-batch | device accumulators, 1 read per step (C4) | - |
| Logging, checkpoints | CPU | CPU | I/O. |

---

## 4. Tasks

Order: C0 first (measure), then C1-C5 on the Mac with equivalence tests, one GPU session (§8), then C6-C9 as arms.
Every task runs behind a temporary `train_rl`/`g_probe` flag until its gate passes; then the old path is deleted
(C10) and the flag removed. Library code in `src/objdel/train/` (E1 later moves it to `rl/caches.py`, `rl/sync.py`,
`rl/forward.py` [SP/impl/notes/E1.md R1]); tests in `REPO/tests/`. Run Python that reads `data/` with `-I`.

Pinned stack (no change in Phase C): torch 2.14.1, torchvision 0.29.1, transformers 5.20.1, trl 1.17.1, peft 0.21.2,
accelerate 1.15.0 (Mac and box); vllm 0.32.0 (box image `vllm_native`). New dependency: none. Box tool: Nsight
Systems CLI `nsys` (presence on the image [unverified]; C0 has a fallback).

New private TRL names touched by Phase C (add each to the Phase A7 start-up check, `src/objdel/train/trl_adapter.py`):
`GRPOTrainer._tokenize_prompts` (C2), `VLLMGeneration.sync_weights` (C5, replaced on the instance),
`GRPOTrainer._compute_loss` source hash (C4 step 6), `trl.trainer.utils.split_pixel_values_by_grid` /
`unsplit_pixel_values_by_grid` (C3).

### C0. Baseline profile and sync count (needs Phase A6 NVTX ranges)

Files:
- change `src/objdel/train/nvtx.py` (created by A6): add `count_syncs()`.
- change `src/objdel/train/fast.py`: `FastGRPO.timed` (fast.py:643-650) logs `perf/host_syncs` per phase when
  counting is on.
- create `experiments/phase2/configs/c_speed.yaml` (the GPU session, §8).
- change `experiments/phase2/train_rl.py` and `g_probe.py`: flags `--count-syncs` (bool), `--profile-steps A-B`
  (calls `torch.cuda.profiler.start()/stop()` around steps A..B for `nsys --capture-range=cudaProfilerApi`).

Steps:
1. `count_syncs()` is a context manager: `torch.cuda.set_sync_debug_mode("warn")`, `warnings.catch_warnings(record=True)`,
   count records whose message says "synchroniz", restore mode 0 on exit. No-op without CUDA. It counts `.item()`,
   `.tolist()`, `nonzero`, D2H and pageable copies; explicit `torch.cuda.synchronize()` may not be counted
   [unverified on 2.14.1; check on the box against nsys `cudaStreamSynchronize` + `cudaDeviceSynchronize` counts].
2. `timed` already syncs at phase end (fast.py:647). With `OBJDEL_NVTX_SYNC=0` (A6) it must record CUDA events and
   read them once in `log`, so the trace shows real overlap.
3. Capture on the box: `nsys profile -t cuda,nvtx,osrt --capture-range=cudaProfilerApi --capture-range-end=stop
   --cuda-memory-usage=true -o c0 python -I -m experiments.phase2.train_rl ... --profile-steps 4-5`. Fallback when
   `nsys` is absent: `g_probe --trace` (existing `callbacks.TraceStep`, callbacks.py:210) plus `--count-syncs`.
4. Report per NVTX range: wall time, GPU-active share, sync count, H2D and D2H bytes. Save `c0_summary.json`
   (`nsys stats --report cuda_gpu_trace,nvtx_sum,cuda_api_sum`).

Gate: the baseline numbers exist (they replace the two-point fit). No speed target in C0.

### C1. Pixels on the GPU from 32 images

Files:
- create `src/objdel/train/gpu_pixels.py` (deep module, one public function plus one error type).
- change `src/objdel/train/fast.py`: `ImageMemo.__call__` (fast.py:458-475) uses `gpu_pixels.patchify`; add
  `ImageMemo.hits`, `ImageMemo.misses` counters; start-up self-check in `FastGRPO.__init__` (fast.py:488).
- create `tests/test_gpu_pixels.py`; extend `tests/test_image_memo.py`.
- create `experiments/phase2/c_speed_check.py` (box checks for C1, C3, C5; module-run).

Interface (`gpu_pixels.py`):
- `class UnsupportedImage(ValueError)`.
- `patchify(images, image_processor, *, device=None, out_dtype=torch.float32, backend="auto") -> (pixel_values, image_grid_thw)`.
  `images`: `torch.uint8` tensor `(N, H, W, 3)` or a sequence of PIL RGB images / `(H, W, 3)` uint8 arrays.
  Output equals `image_processor(list_of_images, return_tensors="pt")` bit for bit (`pixel_values` `(sum t*h*w, 1536)`,
  grid `(N, 3)` int64). Raises `UnsupportedImage` when `smart_resize(H, W) != (H, W)` or the image is not RGB uint8,
  or when the processor's settings differ from the ones checked (resize on with another size, `do_convert_rgb`
  unexpected, patch size not 16, merge not 2, temporal patch not 2).
- `backend`: `"hf"` = the processor's own fast path with `device="cuda"` on a tensor input; `"torch"` = the same ops in
  plain torch (`x.float().sub_(mean').div_(std')`, then HF's patchify reshape/permute
  `(0,2,5,3,6,1,4,7)` and temporal expand [code: image_processing_qwen2_vl.py `patchify`]); `"k3"` = the K3
  `preprocess` extension (built only on a box); `"auto"` = the first backend that passed the start-up self-check.
  `mean'` and `std'` come from HF's own fused-factor function [code: transformers/image_processing_backends.py:300-340],
  never from constants in our code.

Steps:
1. Stack: `np.asarray(pil)` for each distinct image -> one pinned `torch.empty((N,H,W,3), uint8, pin_memory=True)`
   -> `.to(device, non_blocking=True)`. N = 32 per step, 25.2 MB [estimate: arithmetic].
2. `ImageMemo.__call__`: keep the memo by `(id(image), kwargs)` (fast.py:466); compute misses with `patchify` on
   the device; assemble the 256-row output with one `torch.index_select` on the per-image rows (row -> image index
   built on the host). Output `pixel_values` is a CUDA tensor. `BatchFeature(..., tensor_type="pt")` keeps it there.
3. On `UnsupportedImage` call the inner processor (the library path, as today for other processor families).
4. Start-up self-check (once, in `FastGRPO.__init__` when CUDA is present): take the first training image, run
   `patchify(backend=b)` for b in ("hf", "torch") and the inner CPU processor; the first `b` with `torch.equal` on
   both outputs becomes the backend. If none is equal, log `gpu_pixels/backend=cpu` and use the inner processor
   (old speed, same bits). Raise nothing.
5. Pass 1 (`_tokenize_prompts`) no longer calls the image processor at all after C2, so the pass-1/pass-2 id question
   (C1 note §2.3: `load_image` returns new objects) is moot; `hits`/`misses` prove it (expected per step: 32 misses,
   224 hits).
6. K3 fallback only if step 4 fails for both "hf" and "torch" on the box: build `docs/implementation2.0/kernels/K3/k3_torch.cpp`
   per `K3/GPU_RUN.md` steps 1-4, copy `k3_common.cuh`, `k3_kernels.{cu,h}`, `k3_torch.cpp` into
   `REPO/cuda/k3/` with a `torch.utils.cpp_extension.load` wrapper in `gpu_pixels.py`, run compute-sanitizer
   memcheck/racecheck (GPU_RUN step 6). Otherwise K3 (a) is not copied into REPO.
7. Do not optimize the `copy.deepcopy(prompts)` in TRL `_generate` (grpo_trainer.py:2291): Python's deepcopy memo
   copies each distinct PIL object once, so it copies 32 images (25 MB), not 256 [code: grpo_trainer.py:2291;
   Python `copy.deepcopy` memo semantics]. Revisit only if C0 shows it above 50 ms.

Tests:
- Mac (`tests/test_gpu_pixels.py`, CPU device): `patchify(backend="torch", device="cpu")` equals the HF processor on
  (a) 32 random images, (b) all-0, all-255, (c) the 256-value x 3-channel ramp, (d) three "digit images" (pixel =
  byte k of the source element's flat index; equality on these proves the index map, K3 DESIGN §6), (e) a grey image
  as RGB. Shapes 512 x 512 and 96 x 64. `UnsupportedImage` for 500 x 500 and for RGBA.
- Mac: `ImageMemo` returns the same `BatchFeature` as the inner processor for 256 rows of 32 images; counters as above.
- Box (`c_speed_check --pixels`): `torch.equal` of CUDA output vs HF CPU on >= 2,000 views drawn from centred,
  `lidc-train-shifted`, `lidc-train-rescaled`, and on all 32 images of 5 recorded steps; time per call.

Config keys: `train_rl --gpu-pixels {auto,off}` (temporary; default `auto` after the gate).

### C2. Prompt ids and rope positions once per prompt layout

Files:
- change `src/objdel/train/fast.py`: add `FastGRPO._tokenize_prompts` override; replace `RopeCache.keys`
  (fast.py:321-327) and the miss path of `RopeCache.positions` (fast.py:336-348) by `rope_template` +
  `rope_positions`.
- extend `tests/test_train_fast.py`.

Interfaces:
- `rope_template(base, input_ids, mm_token_type_ids, image_grid_thw) -> RopeTemplate` with fields `pos (3, Lp) int64`
  (device), `end: int`, `ids (Lp,) int64`, `types (Lp,) int64`. Built by one HF `get_rope_index` call on one bare
  prompt row. Cached per key `(prompt ids tuple, grid tuple)`.
- `rope_positions(tmpl, input_ids, attention_mask, mm_token_type_ids) -> (positions (3, B, W) int64, bad (B,) bool)`.
  For real token r (mask != 0, counted by `cumsum(mask) - 1`): `pos = tmpl.pos[:, r]` for r < Lp, else
  `tmpl.end + (r - Lp)` on all three axes; masked columns 0 (HF starts from zeros [code: modeling_qwen3_vl.py:954-957,
  C1 note verdict]). `bad[i]` is True when row i's first Lp real tokens differ from `tmpl.ids`/`tmpl.types`, when a
  real token after Lp is not a text token, or when the row has fewer than Lp real tokens. No host read inside.

Steps:
1. `_tokenize_prompts(prompts)`: group rows by `(rendered prompt text, image (H, W))`; for each new key call TRL's
   own `super()._tokenize_prompts([prompt])` once and keep its `prompt_ids[0]`; return `(prompt_ids, images, {})` with
   `images` extracted exactly as TRL does (grpo_trainer.py:1759-1771). The vLLM path drops `multimodal_fields`
   (C1 note S1), so `{}` is safe there; for the non-vLLM path call `super()`.
   The existing `one_pad` check (fast.py:680-686) stays as the guard that vLLM's expansion equals the trainer's ids.
2. `RopeCache.positions(...)`: get the row's template by its prompt layout (from `RowMeta` after C4; before C4 from
   the first Lp ids), call `rope_positions`, accumulate `bad.any()` into a device flag read once in `log`.
   Rows with `bad` take HF `get_rope_index` (old path kept only as the fallback, logged as `rope/fallback_rows`).
3. `--check-rope-every N` (default 50, 0 = off): on step 1 and every N-th step compare with HF on all rows and raise
   on a difference.
4. Pick the plain-torch form (one `cumsum`, one `gather`, one `where`) and not the K3 `rope_positions` kernel: same
   exactness, no build [K3 DESIGN §9 decision "A by default"].

Tests:
- Mac: `torch.equal(rope_positions(...)[0], get_rope_index(...)[0])` on the tiny Qwen3-VL test model and on the real
  processor's prompt (Lp 325, grid (1, 32, 32)), rows with answers of 0, 1, 131, 1,536 tokens, left padding, a mask
  hole, a truncated-and-masked row; Hypothesis property test with random lengths. Three adversarial rows set `bad`.
- Mac: `_tokenize_prompts` returns the same `prompt_ids` as TRL's for 256 rows (tiny model).

Config keys: `train_rl --check-rope-every N`.

### C3. Vision features once per batch; no `pixel_values` in the trainer batch

Files:
- change `src/objdel/train/fast.py`: `VisionCache` (fast.py:245-291) becomes `EncoderCache` with
  `encode(images_patches, grid) -> None` (once per generation batch) and `outputs(image_ids) -> BaseModelOutputWithDeepstackFeatures`;
  `FastForward._get_per_token_logps_and_entropies` (fast.py:102) and `packed_chunk` (fast.py:164-239) pass
  `mm_encoder_outputs={"image": ...}` and no `pixel_values`; `prompt_groups` (fast.py:62-78) drops the pixel compare;
  `per_row` (fast.py:50) removed.
- change `src/objdel/train/fast.py` `FastGRPO._generate_and_score_completions` (fast.py:524): after TRL builds the
  batch, replace `pixel_values` by `image_ids` (int64 per row) and keep `image_grid_thw`.
- patch (module-level swap, same pattern as `shuffle_sequence_dict`, fast.py:561-565)
  `trl.trainer.utils.split_pixel_values_by_grid` / `unsplit_pixel_values_by_grid` to pass the batch through when it has
  `image_ids` and no `pixel_values`.

Steps:
1. Encode the 32 distinct images once per generation batch, before the first forward, under `torch.no_grad()`,
   with the frozen tower (LoRA excludes it [code: configs.py:90-95]). Store per image: `pooler_output` slice
   (rows = t*h*w / merge^2) and each DeepStack level's slice.
2. Per forward chunk: gather the per-row image's tensors in row order into a tuple (`pooler_output`) and a list per
   level (`deepstack_features`), as `Qwen3VLModel.forward` expects (`torch.cat(image_embeds)`,
   modeling_qwen3_vl.py:1183-1186). In the packed path the vision input is one image per group (same as today).
3. Pass `mm_encoder_outputs={"image": out}` and omit `pixel_values` (the model raises if both are given,
   modeling_qwen3_vl.py:1161-1163). `position_ids` always comes from `RopeCache`, so `image_grid_thw` is not needed
   by the model; keep it in the batch for the rope fallback.
4. VRAM: the 1.61 GB `pixel_values` and TRL's unsplit copy (C1 note S4) leave the batch [estimate: frees 1.6-3.2 GB].
   Log `torch.cuda.max_memory_allocated()` after the update.

Tests:
- Mac (tiny Qwen3-VL with DeepStack, existing fixtures in `tests/test_train_fast.py`): per-token logprobs and
  entropies `torch.equal` between the `pixel_values` path and the `mm_encoder_outputs` path, packed and unpacked,
  batch 8 and 4, SDPA and eager; gradients of every LoRA parameter `torch.equal` for one micro-batch.
- Box (`c_speed_check --encoder`): same equality on one recorded generation batch of the real model, reference pass
  and all 16 update chunks.

Config keys: `train_rl --encoder-once` (temporary).

### C4. Host-side metadata; fewer host syncs

Files:
- change `src/objdel/train/fast.py`: add `RowMeta` + `row_meta()`; change `window` (fast.py:36-47), `prompt_groups`
  (fast.py:62-78), `by_length` / `by_group_length` (fast.py:375-404), `packed_chunk` (fast.py:164-239),
  `FastForward` entry (fast.py:108-115), `FastGRPO._generate` (fast.py:618-627).
- optional step 6: `FastGRPO._compute_loss` override for `loss_type == "dr_grpo"`.
- extend `tests/test_train_fast.py`, `tests/test_train_trl.py`.

Interface:
- `@dataclass(frozen=True) RowMeta: prompt_len, completion_len, image_id, group (np.ndarray[int32], shape (B,))`.
- `row_meta(prompt_ids_list, completion_ids_list, image_ids, G, mask_truncated, max_completion_length) -> RowMeta`.
  `completion_len[i] = 0` when the answer is truncated and `mask_truncated` is True (mirrors
  `mask_truncated_completions`, grpo_trainer.py:2604-2611).
- `window_from_meta(meta, rows, logits_to_keep, buckets) -> (start, keep)` equals `window(attention_mask, ...)`.
- `group_order(meta, G) -> np.ndarray[int64]` equals the order `by_group_length` produces; `length_order(meta)`
  equals `by_length`.

Steps:
1. In `FastGRPO._generate` capture `prompt_ids_list` and `completion_ids_list` (host lists TRL already has) and the
   OpenOnce image identity per row; build `RowMeta`; store on `self.meta`.
2. Put `meta` row-aligned into the batch as a Python list field (`host_meta`) so TRL's `split_tensor_dict` /
   `shuffle` and our `move` (fast.py:383-386) carry it with the rows [code: trl/trainer/utils.py:1076-1078].
3. Replace each device read with the host value: `window` start/keep; `prompt_groups` runs = equal `(image_id,
   prompt_len)`; `lengths.tolist()` (fast.py:125); `order.tolist()` (fast.py:381, 397); `grid.prod(-1).tolist()`
   (fast.py:111) from a constant per image size.
4. `packed_chunk`: build `seg`, `idx`, label offsets with numpy into one pinned int64 staging tensor and copy it
   with `non_blocking=True`; build `tok`/`lab`/`mmt`/`posp` with `index_select` from the row tensors instead of the
   per-row `torch.full`/`torch.cat` loop.
5. Debug mode `--check-meta-every N` (default 50): on step 1 and every N-th step run the old device derivation and
   assert equal integers.
6. Loss metrics (decision rule): copy TRL's `_compute_loss` (grpo_trainer.py:3141-3405) into `FastGRPO` for
   `dr_grpo` only, with the seven `.item()` calls (grpo_trainer.py:3376-3405) replaced by device accumulators read once
   in `log` (fast.py:604-616). Guard: SHA-256 of `inspect.getsource(GRPOTrainer._compute_loss)` equals the pinned value
   (A7 mechanism), else raise at start-up. Adopt only if the box A/B shows >= 0.1 s per step saved; else drop the copy.
7. Optional, last (C4b): precompute the fused-head label positions (trl/trainer/utils.py:1694-1697) and the DeepStack
   positions (modeling_qwen3_vl.py:853-854) once per chunk from `RowMeta`, via a version-gated monkeypatch. Adopt
   only after C0/C4 traces show these syncs on the critical path.

Tests:
- Mac: on 50 recorded batches (save them with a `--record-batches DIR` flag on a tiny-model run), the host path and
  the device path give equal `(start, keep)`, runs, order, `seg`, `idx` (integers, exact).
- Mac: `loss` `torch.equal` old vs new; logged metric values equal; a counting patch on `torch.Tensor.item` shows 0
  calls inside `_compute_loss` (step 6).

Config keys: `train_rl --host-meta` (temporary), `--check-meta-every N`.

### C5. Weight push: one call, LoRA-touched weights only

Files:
- create `src/objdel/train/vllm_sync.py`.
- change `src/objdel/train/fast.py`: `FastGRPO.sync_without_merge` (fast.py:495-504) calls `vllm_sync.install`;
  `lora_merged_params` (fast.py:696-717) gains `names: set[str] | None` (None = every weight, as today).
- create `tests/test_vllm_sync.py`.

Interface (`vllm_sync.py`):
- `lora_touched_names(model) -> set[str]`: vLLM-side names of the weights a LoRA layer changes (the
  `*.base_layer.weight` of each `LoraLayer`, mapped with the same rules as `lora_merged_params`). Expected 196
  (28 layers x 7 linear modules, LoRA r16 all-linear, vision excluded [code: configs.py:90-95]) [estimate: count].
- `push_weights(vg, model, *, full: bool) -> PushStats(n_tensors, bytes, seconds)`: if sleeping, `empty_cache()` and
  `wake_up(tags=["weights"])` (TRL's order, vllm_generation.py:500-503); ONE call
  `model_runner.model.load_weights(lora_merged_params(model, names=None if full else touched))`; then
  `llm.reset_prefix_cache()`. A delta push (`full=False`) while the engine sleeps at level 2 raises `ValueError`
  (the untouched weights were discarded).
- `install(trainer) -> bool`: replaces `vg.sync_weights` on the instance (one process, no FSDP; same guards as
  fast.py:500-503). Full push on the first call, after `--resume`, every `--full-push-every N` steps (default 50),
  and on every call while sleep mode is on.

Steps:
1. Keep the values exactly `p.data + delta.to(p.dtype)` (fast.py:717), so vLLM receives the same bf16 W + d as today.
2. Mode A (target): vLLM awake (`vllm_enable_sleep_mode: False`) + delta push. Removes the ~700-call loop, the
   wake (0.11 s) and the sleep (0.29 s) [measured: R5]. Memory: vLLM about 11.2 GB stays resident beside a trainer
   peak of about 17 GB = about 28 GB of 32 GB [doc: Architecture.md §3.4, estimate]; C3 frees 1.6-3.2 GB [estimate].
   Adopt only if G-MEM (<= 30 GB) passes on the box.
3. Mode B (fallback, always valid): sleep mode on (configs.py:193) + full push in ONE `load_weights` call. Keeps
   wake/sleep; removes only the per-call overhead (awake full push measured 1.01 s with ~700 calls [measured: r4
   note]; one call [estimate: well below that, measured on the box]).
4. Do not patch TRL's `sleep(level=2)` to level 1 in Phase C: level 1 restores 4.3 GB of weights from host RAM on each
   wake, which is a large H2D per step and against G-SP5 [estimate: 2.13 B params x 2 B].
5. Box precondition test (`c_speed_check --sync`): in mode A, after 3 delta pushes, copy every vLLM parameter
   (int16 views) and `torch.equal` against a full push of the same adapter state; in mode B, the same check after a
   sleep/wake cycle. Also hash the untouched parameters before and after a delta push (must be equal).
6. Keep the `time/sync_s` timer (fast.py:652-668).
5. Rejected: vLLM native LoRA serving (no merge, 35 MB push) because LoRA kernels slow decode [estimate 10-30 %,
   docs/implementation2.0/notes/B_gpu_audit.md §4.6] and change numerics vs the trainer's merged W + d. Revisit at two GPUs (Phase F).

Tests:
- Mac: a fake vLLM model whose `load_weights` records calls: exactly 1 call per push; names = touched set (delta push)
  or all names (full push); values `torch.equal` to `lora_merged_params`; no parameter of the trainer model changes;
  a delta push on a fake engine marked as sleeping raises `ValueError`.
- Box: the step-3 parameter equality, plus step-1 identity in the A/B (§6 G-EQ5).

Config keys: `train_rl --fast-sync`, `--full-push-every N`; GRPOConfig `vllm_enable_sleep_mode` (False = mode A,
True = mode B), via `--grpo '{"vllm_enable_sleep_mode": false}'`.

### C6. vLLM request shape: 32 requests with `n=8` (measured arm, default off)

Files: create `src/objdel/train/vllm_requests.py` (`group_requests`, `expand_outputs`); change
`FastGRPO.patch_vllm_prompts` (fast.py:670-690) to call them when the flag is on; `tests/test_vllm_requests.py`.

Interface:
- `group_requests(prompts: list[dict], G: int) -> (unique: list[dict], index: list[int])`: consecutive runs of G
  identical prompt dicts (same `prompt_token_ids`, same image object) become one; raises if a run is not G long.
- `expand_outputs(outputs, index, G) -> list`: 256 single-completion outputs in the original row order (child j of
  request k goes to row k*G + j).

Steps:
1. Send `SamplingParams(n=G, ...)` copied from TRL's kwargs (vllm_generation.py:676-697) with `n` replaced.
2. Print `max_num_seqs` and `max_num_batched_tokens` of the engine at start-up (not in the R5 log [measured:
   r5_shared.log:39 has no max_num_seqs]).
3. Expected saving is host-side only (vLLM V1 splits `n` into child requests; prefix caching already dedupes prefill
   [web: vllm PR 10980]): multimodal input processing for 32 instead of 256 requests [unverified size].
4. Not bit-equal (per-request RNG streams change). Gate is statistical (§6 G-ST1).

Config keys: `train_rl --group-requests` (stays a flag unless the gate passes).

### C7. Outline limit (owner decision, run inside Phase B)

Facts: LIDC train doctor outlines: max 1,182 tokens, 0 above 1,200, 13 at or above 800 of 23,469 [measured:
docs/implementation2.0/notes/A1.md F9, counted with the Qwen3-VL tokenizer]. RIGA train: 1 above 1,200 and 1 above 1,536 [same].
In R5, 6 of 15 steps had a naturally ended answer longer than 1,200 tokens (max 1,493); in the full run 23 of 770
steps [measured: runs/r5/shared_rl.json, runs/rl/train_rl/rl.json, counted now]. These answers are longer than any
LIDC doctor outline, but a 1,200 limit turns them into truncated, masked rows.

Files: `experiments/phase2/configs/rl_shifted.yaml` (or the B config) adds `--grpo '{"max_completion_length": L}'`;
change `src/objdel/train/configs.py` `MAX_COMPLETION` only after the owner decides (docstring configs.py:16-21 updated).

Steps:
1. Phase A8 publishes the voting-target length table (A1 R8). The owner picks L (Direction recommends about 1,200).
2. dr_grpo divides by `max_completion_length` (grpo_trainer.py:3286), so the loss scales by 1,536 / L. AdamW is
   invariant to a constant gradient scale except through `eps` and gradient clipping; record this in the run note.
3. Log `completions/clipped_ratio` and the share of answers ended in (L, 1,536] in the previous run's distribution.
4. LIDC-only runs only; a RIGA run keeps 1,536 until its own count is decided.

Expected saving: 0.35-0.77 s per step [estimate, §1 correction 3].

### C8. KL weight 0 (owner decision, run inside Phase B2)

Facts: with `beta == 0` TRL loads no reference adapter and skips the reference pass [code: grpo_trainer.py:526,
1031-1032, 2827]. Our RL starts from a merged SFT checkpoint with a fresh LoRA (`--ckpt` is "a plain model folder",
train_rl.py:24), so no "ref" adapter copy is involved. configs.py:6-7 warns that beta 0 at init blocks raising it later
under full fine-tuning; with LoRA the reference is the adapter-off model, so that warning does not apply
[code: grpo_trainer.py:1033-1036].

Files: B2 config `--grpo '{"beta": 0.0}'`; no code change. Test: `tests/test_train_configs.py` asserts that
`grpo_config(..., beta=0.0)` builds and that a FastGRPO on the tiny model runs one step with no reference forward
(count calls of `_get_per_token_logps_and_entropies`: old and ref passes absent).

Expected saving: 1.4-1.7 s per step [estimate: docs/implementation2.0/notes/B_gpu_audit.md:90].

### C9. Memory-headroom arms (measured, one at a time, after C1-C5 pass)

Arms (each 15 steps, same seed, §8 protocol): (a) partial gradient checkpointing: turn checkpointing off on k = 4, 8
of 28 decoder layers (`GradientCheckpointingLayer.gradient_checkpointing = False` on those layers), keep 16 x 16;
(b) micro-batch 32 x 8 (`per_device_train_batch_size: 32, gradient_accumulation_steps: 8`; same dr_grpo update
[code: grpo_trainer.py:3286]); (c) is C5 mode A itself (vLLM awake), measured in the main fast arm, not here.

Files: `experiments/phase2/configs/c_speed.yaml` arms; `train_rl --no-ckpt-layers K` (new flag) in
`experiments/phase2/train_rl.py` -> `g_probe.make_trainer`.

Rule: adopt an arm only if step time falls by >= 0.5 s and peak memory (`nvidia-smi` max used and
`torch.cuda.max_memory_allocated`) stays <= 30 GB. Today's peaks: vLLM about 11.2 GB, trainer about 17 GB
[doc: REPO/docs/research/research2.0/Architecture.md §3.4]. Estimates: (a) -1 to -2.5 s [estimate: docs/implementation2.0/notes/B_gpu_audit.md
via C3 note F5]; (b) unknown (old measurement slower). Arms (a) and (b) run in mode A only if its peak leaves room.

### C10. Archive, delete, document

Files and actions (after each gate passes):
- delete from `src/objdel/train/fast.py`: `VisionCache.keys` and the `INT` table (C3), `per_row` (C3), pixel compare in
  `prompt_groups` (C3), `RopeCache.keys` (C2), device reads replaced in C4, `FastGRPO.sync_without_merge`'s
  `_iter_named_params` override (C5). Git history plus the owner's CLEAN tag is the copy; add one row per removed
  piece to `REPO/archive/README.md` ("Was / replaced by / run again", CL1 note archive decision).
- remove the temporary flags `--gpu-pixels`, `--encoder-once`, `--host-meta`, `--fast-sync` once their path is the
  default. Keep `--check-*-every`, `--count-syncs`, `--profile-steps`, `--group-requests`, `--no-ckpt-layers`.
- K3: kernels not chosen (lse, is_ratio, rope; preprocess unless C1 step 6 fired) stay in `docs/implementation2.0/kernels/K3/` and are
  not copied into REPO. Record the decision in `REPO/docs/decisions/` (one ADR: "GPU pixel backend and K3").
- update `REPO/docs/architecture/overview.md` (step data flow), `REPO/docs/runbooks/` (nsys capture), module docstring of
  `fast.py` (no "1.16" references; E1 F5), and `src/objdel/train/configs.py` docstring (C7/C8 decisions).
- `make test` and `make lint` pass.

---

## 5. Custom kernels

Folder: `docs/implementation2.0/kernels/K3/` (design `K3/DESIGN.md`, runbook `K3/GPU_RUN.md`). Note: unrelated to
`REPO/experiments/phase2/k3_harness.py` (a reward-speed gate with the same short name).

| K3 piece | Phase C use | Why |
|---|---|---|
| (a) `preprocess` (`k3_kernels.cu`, `k3_torch.cpp` `preprocess`) | Fallback #3 behind `gpu_pixels.patchify(backend="k3")`, only if "hf" and "torch" fail the bit test on the box (C1 step 6) | Owner prefers a library; the HF processor on CUDA and plain torch use the same IEEE sub/div, so they are expected bit-equal [estimate; K3 DESIGN §6]. K3's CPU reference and digit-image test are the proof harness either way. |
| (c) `rope_positions` | Not adopted; plain torch form in C2 | Same exactness, no build [K3 DESIGN §9]. |
| (b) `lse_logprob` | Not in Phase C | Only consumer would be the no-grad reference pass; it disappears with beta = 0 (C8), and the old-logprob pass is already skipped (fast.py:540-547). Revisit only if beta stays non-zero and the K3 bench shows it no slower than TRL. |
| (b) `is_ratio` | Not in Phase C | Saves launches only (< 1 ms [estimate: K3 DESIGN §7]); B2 may change the IS mode. |

K3 status: CPU tests pass and every kernel compiles for sm_120 with 0 spills [measured: K3/build/test_cpu.log,
K3/build/docker_nvcc_check2.log]; GPU bit-equality and timings are not run yet [K3 DESIGN §2].

---

## 6. Acceptance gates

Protocol: one GPU session (§8), R5 settings (`--shared-prefix`, 16 x 16, beta 1e-3, sleep on), same seed, 15 steps
per arm, means over steps 3-15 (as R4/R5). Baseline arm = today's code on the same host.

Equality gates (must pass before any speed number counts):

| Id | Gate | Test |
|---|---|---|
| G-EQ1 | GPU `pixel_values` and grid `torch.equal` to the HF CPU processor on >= 2,000 views (centred, shifted, rescaled) and on all 32 images of 5 recorded steps | `c_speed_check --pixels` |
| G-EQ2 | `rope_positions` `torch.equal` to HF `get_rope_index` on all 256 rows of >= 50 recorded batches; 0 `bad` rows on LIDC | Mac test + `--check-rope-every` on box |
| G-EQ3 | Host `RowMeta` derivations equal the device derivations (integers) on >= 50 recorded batches | Mac test + `--check-meta-every` |
| G-EQ4 | Per-token logprobs and entropies `torch.equal` between the `pixel_values` path and the `mm_encoder_outputs` path, reference pass and 16 update chunks of one recorded batch; one micro-batch's LoRA gradients `torch.equal` | `c_speed_check --encoder` |
| G-EQ5 | Mode A: after 3 delta pushes, every vLLM parameter `torch.equal` to a full push of the same adapter state; untouched parameters unchanged (hash). Mode B: same after a sleep/wake cycle | `c_speed_check --sync` |
| G-EQ6 | Step 1 of the fast arm equals the baseline arm: reward, legal, IoU_tc, mean answer length, grad norm (same digits as R5: 0.5467 / 0.8984 / 117.89 / 0.0033 on e2) | A/B log |

Speed and resource gates (fast arm with C1-C5):

| Id | Gate | Basis |
|---|---|---|
| G-SP1 | Mean step time <= 19.5 s in mode A (vLLM awake) or <= 20.5 s in mode B (sleep). Fail (re-profile; Phase B keeps today's code) if > 20.7 s (A) or > 21.9 s (B) | Estimated saving 2.7-5.6 s (A) or 1.55-5.0 s (B) for C1-C6; gate = 23.43 s minus the midpoint, rounded up to 0.5 s; fail = 23.43 s minus the low end [estimate, §8.1] |
| G-SP2 | Mode A: `time/sync_s` <= 0.2 s and no wake/sleep (from 1.27 + 0.11 + 0.29 s). Mode B: `time/sync_s` < 1.01 s (the measured awake full push with ~700 calls) | C5 [estimate: C1 note R-4 "~0.1 s"; measured: r4 note] |
| G-SP3 | NVTX `prepare` range (TRL pass 2 + H2D) <= 0.5 s | Architecture.md §8.2 alarm "prepare above 0.5 s" |
| G-SP4 | Host syncs per step < 200 (from 900-1,250) | C1 note G0 target |
| G-SP5 | H2D bytes per step < 50 MB (from >= 1.61 GB) | C1 note G0 target; images are 25.2 MB |
| G-SP6 | GPU-active share of the step >= 70 % (from about 60 %) | BRIEF §2 idle about 40 %; Architecture.md target 80 % stays the end goal after decode-tail work |
| G-MEM | Peak GPU memory <= 30 GB on every arm | SD1 §7.1 rule |
| G-ST1 (C6 only) | Over steps 3-15 vs baseline: mean reward, legal, IoU_tc within the baseline's step-to-step spread (2 standard deviations), answer-length mean within 5 %; and >= 0.2 s per step saved, else drop C6 | statistical (not bit-equal) |
| G-C | With C7 and C8 added (owner decisions, inside Phase B): mean step <= 18.0 s | Direction.md Gate C value; estimate 15.4-19.0 s (A) or 16.0-20.1 s (B) (§8.1) |

Code gates: `make test` and `make lint` pass on the Mac; no test needs a GPU unless marked `cuda`.

---

## 7. Risks and rollback

| Risk | Signal | Mitigation | Rollback |
|---|---|---|---|
| HF/torchvision drift changes the pixel formula or patch order | start-up self-check unequal | backend falls to CPU automatically; digit-image test in CI | `--gpu-pixels off` |
| A training view is not 512 x 512 (rescaled views) or not RGB | `UnsupportedImage` | per-image CPU fallback; log count | none needed |
| Mixed LIDC + RIGA prompts in one batch | `bad` rows | one template per layout; HF fallback per bad row | `--check-rope-every 1` to find it |
| `mm_encoder_outputs` order wrong in packed chunks (DeepStack) | G-EQ4 fails | debug order: one image, unpacked, then packed, then DeepStack levels one by one | `--encoder-once` off (old VisionCache) |
| Host metadata out of step with rows after TRL shuffles | G-EQ3 / periodic check raises | `host_meta` travels inside the batch dict | `--host-meta` off |
| Mode A (awake) does not fit in memory | G-MEM fails | C3 frees 1.6-3.2 GB first; else mode B (sleep + one-call full push) | `vllm_enable_sleep_mode: true` |
| Delta push leaves a vLLM weight stale | G-EQ5 fails | full push every 50 steps and on resume; parameter hash check | `--fast-sync` off |
| Stale prefix cache after push | logprob gap rises (`sampling/sampling_logp_difference/mean` vs baseline 0.011-0.018 [measured: r4 note, BRIEF]) | `reset_prefix_cache()` kept after every push | as above |
| TRL private names change | A7 start-up check fails | pin trl 1.17.1 | none (run does not start) |
| `_compute_loss` copy diverges from TRL | source-hash guard | only adopted if >= 0.1 s saved | drop the copy |
| `n=8` changes the sampling stream | G-ST1 | default off | flag off |
| Outline limit cuts answers that end between L and 1,536 tokens | share logged (6/15 R5 steps had one) | owner decides L from A8 table | restore 1,536 |
| beta 0 lets the policy drift | KL no longer logged; watch entropy and legal rate | owner decision inside B2 with its gates | beta 1e-3 |
| Memory: awake vLLM or larger micro-batch OOM | G-MEM | arms one at a time | config back |

Debug order when G-SP1 fails: read the C0-style trace of the fast arm; compare each NVTX range with the baseline;
the range that did not shrink names the task; then its sync count (`perf/host_syncs`) tells whether a host read came
back.

Rollback of code: each task is one commit (owner commits); revert = the owner runs `git revert <sha>`. During A/B the
flags give a no-code rollback.

---

## 8. Cost, time and dependencies

### 8.1 Step-time estimate (from R5 23.43 s)

| Task | Saving per step | Tag |
|---|---|---|
| C1 pixels on GPU + C2 prompt-ids memo | 0.75-1.6 s | [estimate: C1 note R-1 parts: 2 x 1.61 GB CPU copies 0.3-0.5, pass-1 processing 0.2-0.5, pageable H2D 0.15-0.3, tokenization 0.1-0.3] |
| C2 rope | 0.1-0.3 s | [estimate: C1 note R-3] |
| C3 encoder outputs once | 0.1-0.3 s | [estimate: C1 note R-1 level B] |
| C4 host metadata, syncs | 0.3-0.8 s | [estimate: C1 note R-2] |
| C5 weight push | mode A: 1.45-1.6 s (sync 1.27 + wake 0.11 + sleep 0.29 measured, minus a new sync of 0.05-0.2 s); mode B: 0.3-1.0 s [unverified: one-call full push cost] | [estimate] |
| C6 `n=8` | 0-1 s | [unverified; expected low end, §4 C6 step 3] |
| **C1-C6** | **mode A 2.7-5.6 s -> 17.8-20.7 s; mode B 1.55-5.0 s -> 18.4-21.9 s** | [estimate: sum] |
| C7 limit 1,200 | 0.35-0.77 s | [estimate, §1] |
| C8 beta 0 | 1.4-1.7 s | [estimate: two-point fit] |
| **C1-C8** | **mode A 4.45-8.1 s -> 15.4-19.0 s; mode B 3.3-7.5 s -> 16.0-20.1 s** | [estimate: sum] |
| C9 arms | (a) 1-2.5 s; (b) unknown; (c) <= 0.4 s | [estimate] |

### 8.2 GPU session (one rental; owner approval and preflight first; North America host)

| Item | GPU minutes | Basis |
|---|---|---|
| Box setup (image pull, sync, model fetch, env check) | about 40-50 | [estimate: R4 cost $0.65 at $0.593/h = 66 min for two 15-step arms of about 7.5 min each, measured: r4 note] |
| `c_speed_check --pixels --encoder --sync` | about 5 | [estimate: 2,000 views, one batch, one push] |
| C0 baseline arm with nsys capture on steps 4-5 | about 8 | [measured: R5 base arm 7 min 45 s for 15 steps] + capture overhead [estimate] |
| Fast arm (C1-C5), 15 steps | about 6 | [measured: R5 shared arm 5 min 54 s] scaled down [estimate] |
| Optional C6 arm, 15 steps | about 6 | [estimate] |
| Optional K3 GPU_RUN (only if C1 step 6 fires) | <= 20 | [estimate: K3/GPU_RUN.md] |
| **Total** | **about 65-95 min = $0.65-0.95 at $0.593/h** | [estimate; rate measured: r4/r5 notes] |

C7, C8, C9 cost no extra rental: they ride inside the Phase B runs (C9 arms only if B's budget allows; each 15-step
arm about 6 min = about $0.06 [estimate]).

Payback: a 300-step B run at 23.6 s/step [measured: R5 shared] takes 118 min; at 19.5 s/step (G-SP1, mode A) it takes
98 min (-20 min, about -$0.20 per 300 steps at $0.593/h) [estimate]. Over Phase B's three runs this returns about the session's cost.

Budget note: BRIEF §3 sets $3 for GPU scripts; Direction.md §8 shows $0.60 left under the $17.50 limit. The owner
must approve this session against the B budget (Direction §6 decision 1).

Mac time (engineer days, [estimate]): C0 0.5; C1 1.5; C2 1; C3 1.5; C4 2; C5 1; C6 0.5; C7-C9 config only; C10 0.5.

### 8.3 Dependencies

| Depends on | What | Why |
|---|---|---|
| Phase A6 | NVTX ranges in `FastGRPO.timed`, event-based timing option | C0 and every speed gate |
| Phase A7 | TRL private-name start-up check (`trl_adapter.py`) | C2, C4, C5, C3 touch more private names (§4 list) |
| Phase A8 | voting-target outline-length table | C7 limit |
| Phase B1 data | `lidc-train-shifted`, `lidc-train-rescaled` views (exist locally: T1 F3) | G-EQ1 covers every view kind B1 trains on |
| Phase B2 | IS mode and beta decision | C8 runs inside B2; C must pass under `sequence_mask` and `token_truncate` |
| T1 view bank (optional) | uint8 planes instead of PNG | `patchify` accepts `(N,H,W,3)` uint8; grey -> RGB expand on device |
| Phase D | CUDA checker / polygonizer | independent; rewards already on GPU |
| Phase E | file moves (`rl/caches.py`, `rl/sync.py`, `rl/forward.py`) | land C first, then E moves the new modules unchanged |

Phase C blocks: the Phase B runs should use the C code only after G-EQ1..6 and G-SP1 pass; otherwise B runs on
today's code (correct, slower).
