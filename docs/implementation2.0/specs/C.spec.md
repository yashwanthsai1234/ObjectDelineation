# Phase C spec: GPU speed (revision 2)

## 1. Goal

Make one GRPO step faster with the same numbers. Tasks C1-C5 must not change any sampled token, log-prob, loss or
gradient bit. C6 (and the C6b investigation) has a statistical gate. C7 and C8 are config changes chosen by the owner.

Stacks (fixed). Mac: torch 2.14.1, torchvision 0.29.1, transformers 5.20.1, trl 1.17.1, peft 0.21.2, accelerate 1.15.0.
GPU box image `vllm_native`: torch 2.13.0, torchvision 0.28.0, transformers 5.19.0, trl 1.17.1, peft 0.21.2,
accelerate 1.15.0, vllm 0.32.0. Every requirement below must hold on both stacks; bit-equality gates are judged on the
box against the box's own HF CPU path. Model: Qwen3-VL-2B with LoRA r16 on the language model (vision tower frozen).
Step shape: 32 prompts x 8 answers = 256 rows; one 512 x 512 image per prompt (RGB for LIDC, grey `L` for RIGA).

## 2. Modules and interfaces

All library code is in `src/objdel/train/`; tests are in `tests/`. Tests marked `cuda` run only on a GPU box; every
other test runs on a Mac CPU.

### 2.1 `objdel.train.gpu_pixels`

```python
class UnsupportedImage(ValueError): ...

def patchify(images, image_processor, *, device=None, out_dtype=torch.float32,
             backend: str = "auto") -> tuple[torch.Tensor, torch.Tensor]: ...
```

- `images`: a `torch.uint8` tensor of shape `(N, H, W, 3)` or `(N, H, W)` (grey), or a sequence of N PIL images in
  mode `RGB` or `L`, or `(H, W, 3)` / `(H, W)` uint8 numpy arrays.
- `image_processor`: the Qwen2-VL/Qwen3-VL image processor loaded by `AutoProcessor` for the model.
- Returns `(pixel_values, image_grid_thw)`: `pixel_values` shape `(sum over images of t*h*w, 1536)`, dtype
  `out_dtype`, on `device` (default: CUDA when available, else CPU); `image_grid_thw` shape `(N, 3)`, int64, ALWAYS on
  the CPU whatever `device` is.
- Required: both outputs are bit-equal (`torch.equal`, after moving to one device) to
  `image_processor(list_of_images, return_tensors="pt")` for the same images, for every backend that returns without
  raising. A grey input gives the same result as the same image converted with `PIL.Image.convert("RGB")`.
- `backend` in `{"auto", "hf", "torch", "k3"}`; `"torch"` must work on CPU and CUDA; an unknown name raises
  `ValueError`.
- Raises `UnsupportedImage` when the processor would resize the image (`smart_resize(H, W) != (H, W)`), when the image
  mode is not `RGB` or `L` (for example RGBA, P, I;16) or the array has another channel count, or when the processor
  settings are not the checked ones (patch size 16, merge size 2, temporal patch size 2).
- Valid sizes include 512 x 512 and 96 x 64. 500 x 500 raises `UnsupportedImage`.
- Normalisation constants come from the processor object, never from literals in our code.

### 2.2 `objdel.train.fast.ImageMemo` (changed behaviour)

- `ImageMemo.__call__(images, **kwargs) -> BatchFeature`.
- Attribute `per_image: bool` (default False).
  - `per_image == False`: returns the same keys and the same values (`torch.equal`) as the wrapped processor for the
    same call.
  - `per_image == True`: `pixel_values` holds each distinct image object's patches exactly once, in first-appearance
    order, and equals the wrapped processor's `pixel_values` for the list of distinct images; `image_grid_thw` has one
    row per input image (as the wrapped processor's).
- When CUDA is available, `pixel_values` in the result is a CUDA tensor; `image_grid_thw` is always a CPU tensor.
- `prefill(images) -> None` processes the distinct image objects not yet in the memo.
- Each distinct image object is processed at most once between two `clear()` calls.
- Counters `hits: int` and `misses: int` count per-row lookups (and `prefill` counts one miss per newly processed
  image). For 256 rows over 32 distinct images in one call after `clear()`: `misses == 32`, `hits == 224`. After
  `prefill` of those 32 images, the same call adds `hits == 256` and no misses.
- An `UnsupportedImage` falls back to the wrapped processor for the whole call; the result still equals the wrapped
  processor's output.

### 2.3 GPU pixel start-up check (`objdel.train.fast.FastGRPO`)

- When CUDA is present, start-up compares `patchify` (backends `"hf"` then `"torch"`) with the wrapped CPU processor
  on 8 images (training views of each kind present, a grey image when RIGA rows are present, one generated digit
  image) and selects the first backend equal on all 8.
- If no backend passes: `RuntimeError` at start-up, unless `--allow-cpu-pixels` is given; then the CPU processor is
  used and `gpu_pixels/backend=cpu` is logged.

### 2.4 Rope positions (`objdel.train.fast`)

```python
@dataclass(frozen=True)
class RopeTemplate:
    pos: torch.Tensor    # (3, Lp) int64
    end: int
    ids: torch.Tensor    # (Lp,) int64
    types: torch.Tensor  # (Lp,) int64

def rope_template(base, input_ids, mm_token_type_ids, image_grid_thw) -> RopeTemplate: ...
def rope_positions(tmpl: RopeTemplate, input_ids, attention_mask, mm_token_type_ids
                   ) -> tuple[torch.Tensor, torch.Tensor]: ...
```

- `rope_template` takes one prompt row (`input_ids` shape `(1, Lp)`) and the model's Qwen-VL base (the module with
  `get_rope_index`).
- `rope_positions` returns `(positions, bad)`: `positions` shape `(3, B, W)` int64 on the input device; `bad` shape
  `(B,)` bool.
- Required: for every row whose real tokens (attention mask 1) start with exactly the template's `Lp` ids and token
  types and continue with text tokens only, `positions[:, i]` equals `base.get_rope_index(...)[0][:, i]` bit for bit,
  including left padding, masked holes, 0-token, 1-token, 131-token and 1,536-token answers. Masked columns are 0.
- `bad[i]` is True when the row's first `Lp` real tokens differ from the template (ids or types), when a real token
  after `Lp` is not a text token, or when the row has fewer than `Lp` real tokens.
- No host synchronisation inside `rope_positions` (no `.item()`, `.tolist()`, `.cpu()`, boolean-mask indexing).

### 2.5 Prompt ids (`objdel.train.fast.FastGRPO._tokenize_prompts`)

- On the vLLM path, returns `(prompt_ids, images, {})` where `prompt_ids` equals TRL 1.17.1
  `GRPOTrainer._tokenize_prompts(prompts)[0]` and `images` equals its `[1]` (same objects, same order).
- The image processor is not called by this method.
- TRL's tokenization runs at most once per distinct (rendered prompt text, image size) per run.
- On the non-vLLM path, the result equals TRL's method.

### 2.6 `objdel.train.fast.PromptMemo`

```python
class PromptMemo:
    @classmethod
    def install(cls, processing_class, image_memo) -> "PromptMemo | None": ...
```

- `install` returns None and changes nothing for a processor family other than Qwen-VL.
- After `install`, `isinstance(processing_class, T)` is still True for every type `T` it was an instance of, and
  attribute access is unchanged.
- For a call `processing_class(images=..., text=..., padding=True, return_tensors="pt")`, the result has the same keys
  and values (`torch.equal`) as the processor before `install`, including mixed image sizes and mixed prompt texts in
  one call.
- The original processor runs at most once per distinct (text, image (H, W)) per run; `pixel_values` and
  `image_grid_thw` come from `image_memo`.
- Any other call signature behaves exactly as before `install`.

### 2.7 Vision encoder outputs once per generation batch (`objdel.train.fast`)

- The frozen vision tower is called exactly once per generation batch (one call over the distinct images), covering
  the reference pass and every update micro-batch.
- The trainer batch returned by `FastGRPO._generate_and_score_completions` has no `pixel_values` key; it has
  `host_meta` (section 2.8) and keeps `image_grid_thw`.
- Every model forward in the reference pass and the update passes `mm_encoder_outputs={"image": ...}` and no
  `pixel_values`.
- After TRL shuffles and splits the batch, each micro-batch's per-row image ids match its rows.
- If the per-row image ids are missing while this path is on, the forward raises (no silent fallback).
- Required: per-token log-probs and entropies are `torch.equal` to the `pixel_values` path, packed (`shared_prefix`)
  and unpacked, SDPA and eager, chunk sizes 8 and 4; LoRA gradients of one micro-batch are `torch.equal`.
- With this path on, the generation batch never holds a `(256 x 1024, 1536)` pixel tensor.

### 2.8 Row metadata (`objdel.train.fast`)

```python
@dataclass(frozen=True)
class RowMeta:
    prompt_len: np.ndarray      # (B,) int32
    completion_len: np.ndarray  # (B,) int32
    image_id: np.ndarray        # (B,) int32
    group: np.ndarray           # (B,) int32

def row_meta(prompt_ids_list, completion_ids_list, image_ids, G: int, mask_truncated: bool,
             max_completion_length: int) -> RowMeta: ...
def window_from_meta(meta: RowMeta, rows, logits_to_keep: int, buckets) -> tuple[int, int]: ...
def length_order(meta: RowMeta) -> np.ndarray: ...
def group_order(meta: RowMeta, G: int) -> np.ndarray: ...
```

- `completion_len[i] == 0` when answer i reached `max_completion_length` without an end token and `mask_truncated`
  is True; otherwise the answer's token count. `group[i] == i // G`. Rows with the same image object share an
  `image_id`.
- `window_from_meta(...)` equals `window(attention_mask, logits_to_keep, buckets)` for the same rows.
- `length_order(meta)` equals the row order `by_length` produces; `group_order(meta, G)` equals the order
  `by_group_length(batch, G)` produces; whole groups stay together in `group_order`.
- The packed-chunk indices (`seg`, `idx`) built from `RowMeta` equal the ones built from device tensors.
- These functions read no device tensor.
- The trainer batch carries the metadata as a row-aligned Python list field `host_meta`.

### 2.9 Loss metrics (only if adopted, see gate rule in section 5)

- `FastGRPO._compute_loss` for `loss_type == "dr_grpo"` calls `torch.Tensor.item` 0 times.
- The loss is `torch.equal` to TRL 1.17.1's `_compute_loss` on the same inputs.
- After `FastGRPO.log`, the logged values of `entropy`, `kl` (when beta != 0), `clip_ratio/low_mean`,
  `clip_ratio/high_mean`, `clip_ratio/region_mean`, `clip_ratio/low_min`, `clip_ratio/high_max` equal TRL's.
- At start-up, if `inspect.getsource(GRPOTrainer._compute_loss)` does not hash to the pinned SHA-256, the trainer raises
  `RuntimeError` naming the method.

### 2.10 `objdel.train.vllm_sync`

```python
@dataclass(frozen=True)
class PushStats:
    n_tensors: int
    bytes: int
    seconds: float

def push_weights(vg, model) -> PushStats: ...
def install(trainer) -> bool: ...
```

- `push_weights` calls the vLLM model's `load_weights` exactly once per push, with an iterable holding every
  `(name, tensor)` that `lora_merged_params(model)` yields, with equal values.
- When the engine uses sleep mode, `push_weights` calls `wake_up(tags=["weights"])` before `load_weights`.
- After loading, `reset_prefix_cache()` is called once.
- No parameter of the trainer model changes.
- `PushStats.seconds` is the duration of the `load_weights` call only; `n_tensors` and `bytes` describe the iterable.
- `install(trainer)` replaces `trainer.vllm_generation.sync_weights` on that instance and returns True; it returns
  False (and changes nothing) when there is no vLLM generation, no PEFT model, FSDP, or more than one process. Every
  installed sync is a full push.

### 2.11 `objdel.train.vllm_requests`

```python
def group_requests(prompts: list[dict], G: int) -> list[dict]: ...
```

- Returns `prompts[::G]` when the list is consecutive runs of exactly G identical prompt dicts (equal
  `prompt_token_ids`, same image object); raises `ValueError` if a run is not exactly G long or the length is not a
  multiple of G.
- With `train_rl --group-requests`, vLLM receives `len(prompts) / G` requests with `n == G`, and the trainer receives
  `len(prompts)` completions in the original row order (completion j of request k is row `k * G + j`).
- Default off.

### 2.12 Measurement (`objdel.train.nvtx`, created by Phase A6)

```python
@contextlib.contextmanager
def count_syncs() -> Iterator[Counter]: ...
```

- Without CUDA: a no-op that yields a counter with value 0 and raises nothing.
- With CUDA: counts every synchronizing operation reported by `torch.cuda.set_sync_debug_mode("warn")` inside the block,
  including repeats from the same code line (1,000 `.item()` calls in a loop count 1,000), and restores the previous
  mode on exit (also on an exception).
- `FastGRPO.timed(phase)` logs `perf/host_syncs` per phase when `--count-syncs` is on, and still logs
  `time/<phase>_s` as today. Syncs inside the `generate` phase are not in `perf/host_syncs`; they are logged as
  `perf/host_syncs_vllm`.
- The weights wake and the KV-cache wake are timed apart: `time/wake_weights_s`, `time/wake_kv_s`.

### 2.13 GPU-active share (`experiments/phase2/gpu_active.py`, module-run)

- GPU-active share = (length of the union of all kernel and memcpy/memset intervals on the device, over all streams,
  inside a given time range) / (length of that range).
- Accepts an nsys `cuda_gpu_trace` CSV and a torch.profiler Chrome trace JSON; both give the same share for the same
  intervals.
- Unit test (Mac): a synthetic trace with overlapping intervals on two streams gives the exact expected share.

## 3. Command-line flags (`experiments/phase2/train_rl.py` and `g_probe.py`)

| Flag | Values | Behaviour |
|---|---|---|
| `--gpu-pixels` | `auto` (default), `off` | `off` keeps the CPU processor path |
| `--allow-cpu-pixels` | bool, default off | allows start-up to continue on the CPU processor when no GPU backend passes (section 2.3) |
| `--prompt-memo` | bool | installs section 2.6 |
| `--encoder-once` | bool | section 2.7 path (sets `ImageMemo.per_image`) |
| `--host-meta` | bool | section 2.8 path |
| `--fast-sync` | bool | installs section 2.10 |
| `--group-requests` | bool, default off | section 2.11 |
| `--check-rope-every` | int, default 50; 0 = off | compares rope positions with HF on step 1 and every N-th step; raises on a difference |
| `--check-meta-every` | int, default 50; 0 = off | compares host metadata with the device derivation; raises on a difference |
| `--count-syncs` | bool | section 2.12 |
| `--profile-steps` | `A-B` | `torch.cuda.profiler.start()` before step A, `stop()` after step B |
| `--record-batches` | directory | saves each generation batch's ids, masks, grids, metadata and positions |
| `--no-ckpt-layers` | int K, default 0 | turns gradient checkpointing off on K decoder layers |

GRPOConfig keys used (through `--grpo` JSON): `vllm_enable_sleep_mode`, `max_completion_length`, `beta`,
`per_device_train_batch_size`, `gradient_accumulation_steps`.

## 4. Constraints

- No new Python dependency. TRL stays 1.17.1. The box image gains the Nsight Systems CLI (`nsys`).
- No test reads `data/` without `python -I`. No test needs a GPU unless marked `cuda`.
- Bit tests use `torch.equal` (or exact integer equality), not tolerances, except the C6 statistical gate and the two
  4-significant-digit values in G-EQ6.
- The trainer and vLLM must work with `vllm_importance_sampling_mode` set to `sequence_mask` and to `token_truncate`.
- Every private TRL name the new code overrides or patches (`_tokenize_prompts`, `_compute_loss`,
  `VLLMGeneration.sync_weights`) is in the Phase A7 start-up check.
- Each new wrapper (rope, encoder, prompt memo) installs only for Qwen-VL-style processors/models and otherwise leaves
  TRL's path unchanged.
- `make test` and `make lint` pass.

## 5. Acceptance gates

Protocol: one GPU box session, same host and seed for all arms, R5 settings (`--shared-prefix`, 16 x 16, beta 1e-3).
Timed arms run 17 steps; times are the mean of steps 3-15; nsys captures steps 16-17 only. A separate 15-step check arm
(checks every step, batches recorded) supplies box batches; its times are not used.

Equality (GPU box, real model; Mac tests cover the same with the tiny model or the real processor where possible):

| Id | Gate |
|---|---|
| G-EQ1 | `patchify` output (pixels, CPU grid) `torch.equal` to the box's HF CPU processor on >= 2,000 views (centred, shifted, rescaled) and all 32 images of 5 recorded steps; 0 syncs around one pass-2 processor call |
| G-EQ2 | `rope_positions` `torch.equal` to HF `get_rope_index` on all 256 rows of >= 50 batches, >= 15 of them from the box; 0 `bad` rows on LIDC |
| G-EQ3 | `RowMeta`-based windows, orders, runs and packed indices equal the device-based ones on >= 50 batches, >= 15 of them from the box |
| G-EQ4 | Log-probs and entropies `torch.equal` between the `pixel_values` path and the `mm_encoder_outputs` path (reference pass and 16 update chunks of one batch); LoRA gradients of one micro-batch `torch.equal`; vision tower called once per generation batch |
| G-EQ5 | After one `push_weights`, every vLLM parameter `torch.equal` to the state after TRL's per-parameter push of the same adapter state (vLLM awake: 3 pushes on 3 adapter states; sleep mode: after a sleep/wake cycle) |
| G-EQ6 | Step 1 of the fast arm vs the baseline arm: `rewards/iou_tc/mean`, `rewards/legal/mean`, `completions/mean_length` equal; `grad_norm` and `sampling/sampling_logp_difference/mean` equal to 4 significant digits; `reward` not compared |

Speed and resources (fast arm = C1-C5 only):

| Id | Gate |
|---|---|
| G-SP1 | Mean step time <= 20.0 s with vLLM awake, or <= 20.7 s with sleep mode (baseline 23.43 s); fail above 21.2 s (awake) or 21.9 s (sleep). Thresholds are recomputed from the measured one-call push time before the fast arm runs |
| G-SP2 | `PushStats.seconds` <= 0.2 s; vLLM awake: no wake or sleep phase; sleep mode: `time/sync_s` <= `PushStats.seconds` + baseline `time/wake_weights_s` + 0.1 s |
| G-SP3 | NVTX `prepare` range <= 0.5 s per step |
| G-SP4 | Trainer-side host synchronisations per step (outside `generate`) < 200 and <= 25 % of the baseline's trainer-side count |
| G-SP5 | Trainer-side host-to-device bytes per step < 50 MB (vLLM's own copies reported, not gated) |
| G-SP6 | GPU-active share (section 2.13) of the step >= 70 % |
| G-MEM | Peak GPU memory <= 30 GB on every arm |
| G-ST1 | `--group-requests` arm: mean reward, legal and IoU_tc within 2 standard deviations of the baseline's step-to-step values; mean answer length within 5 %; >= 0.2 s per step saved |
| G-C | With the owner's outline limit and beta = 0 added: mean step time <= 18.0 s |
| Loss-metric copy (2.9) | kept only if it saves >= 0.1 s per step in the box A/B |
| `PromptMemo` (2.6) | kept only if the baseline trace shows the pass-2 text work at >= 50 ms per step |
