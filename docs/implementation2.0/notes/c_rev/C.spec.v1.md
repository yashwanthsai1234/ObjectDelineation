# Phase C spec: GPU speed

## 1. Goal

Make one GRPO step faster with the same numbers. Tasks C1-C6 must not change any sampled token, log-prob, loss or
gradient bit (C6 is the one exception and has a statistical gate). C7 and C8 are config changes chosen by the owner.

Stack (fixed): torch 2.14.1, torchvision 0.29.1, transformers 5.20.1, trl 1.17.1, peft 0.21.2, accelerate 1.15.0;
vllm 0.32.0 on the GPU box only. Model: Qwen3-VL-2B with LoRA r16 on the language model (vision tower frozen).
Step shape: 32 prompts x 8 answers = 256 rows; one 512 x 512 RGB image per prompt.

## 2. Modules and interfaces

All library code is in `src/objdel/train/`; tests are in `tests/`. Tests marked `cuda` run only on a GPU box; every
other test runs on a Mac CPU.

### 2.1 `objdel.train.gpu_pixels`

```python
class UnsupportedImage(ValueError): ...

def patchify(images, image_processor, *, device=None, out_dtype=torch.float32,
             backend: str = "auto") -> tuple[torch.Tensor, torch.Tensor]: ...
```

- `images`: a `torch.uint8` tensor of shape `(N, H, W, 3)`, or a sequence of N PIL RGB images or `(H, W, 3)` uint8
  numpy arrays.
- `image_processor`: the Qwen2-VL/Qwen3-VL image processor loaded by `AutoProcessor` for the model.
- Returns `(pixel_values, image_grid_thw)`: `pixel_values` shape `(sum over images of t*h*w, 1536)`, dtype
  `out_dtype`, on `device` (default: CUDA when available, else CPU); `image_grid_thw` shape `(N, 3)`, int64.
- Required: both outputs are bit-equal (`torch.equal`) to `image_processor(list_of_images, return_tensors="pt")`
  for the same images, for every backend that returns without raising.
- `backend` in `{"auto", "hf", "torch", "k3"}`; `"torch"` must work on CPU and CUDA; an unknown name raises
  `ValueError`.
- Raises `UnsupportedImage` when the processor would resize the image (`smart_resize(H, W) != (H, W)`), when the image
  is not 3-channel uint8 (for example RGBA or grey with 1 channel), or when the processor settings are not the
  checked ones (patch size 16, merge size 2, temporal patch size 2).
- Valid sizes include 512 x 512 and 96 x 64. 500 x 500 raises `UnsupportedImage`.

### 2.2 `objdel.train.fast.ImageMemo` (changed behaviour)

- `ImageMemo.__call__(images, **kwargs) -> BatchFeature` returns the same keys and the same values (`torch.equal`)
  as the wrapped processor for the same call.
- When CUDA is available, `pixel_values` in the result is a CUDA tensor.
- Each distinct image object is processed at most once between two `clear()` calls.
- Counters `hits: int` and `misses: int` exist and count per-row lookups; for 256 rows over 32 distinct images in one
  call after `clear()`: `misses == 32`, `hits == 224`.
- An `UnsupportedImage` falls back to the wrapped processor for the whole call; the result still equals the wrapped
  processor's output.

### 2.3 Rope positions (`objdel.train.fast`)

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

### 2.4 Prompt ids (`objdel.train.fast.FastGRPO._tokenize_prompts`)

- On the vLLM path, returns `(prompt_ids, images, {})` where `prompt_ids` equals TRL 1.17.1
  `GRPOTrainer._tokenize_prompts(prompts)[0]` and `images` equals its `[1]` (same objects, same order).
- The image processor is not called by this method.
- The chat template runs at most once per distinct (rendered prompt text, image size) per run.
- On the non-vLLM path, the result equals TRL's method.

### 2.5 Vision encoder outputs once per generation batch (`objdel.train.fast`)

- After this change, the trainer batch built by `FastGRPO._generate_and_score_completions` has no `pixel_values` key;
  it has `image_ids` (int64, one per row; rows of the same image share an id) and keeps `image_grid_thw`.
- The frozen vision tower runs once per distinct image per generation batch.
- Forward chunks call the model with `mm_encoder_outputs={"image": ...}` and without `pixel_values`.
- Required: per-token log-probs and entropies are `torch.equal` to the `pixel_values` path, packed (`shared_prefix`)
  and unpacked, SDPA and eager, chunk sizes 8 and 4; LoRA gradients of one micro-batch are `torch.equal`.

### 2.6 Row metadata (`objdel.train.fast`)

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
  is True; otherwise the answer's token count. `group[i] == i // G`.
- `window_from_meta(...)` equals `window(attention_mask, logits_to_keep, buckets)` for the same rows.
- `length_order(meta)` equals the row order `by_length` produces; `group_order(meta, G)` equals the order
  `by_group_length(batch, G)` produces; whole groups stay together in `group_order`.
- The packed-chunk indices (`seg`, `idx`) built from `RowMeta` equal the ones built from device tensors.
- These functions read no device tensor.

### 2.7 Loss metrics

- `FastGRPO._compute_loss` for `loss_type == "dr_grpo"` calls `torch.Tensor.item` 0 times.
- The loss is `torch.equal` to TRL 1.17.1's `_compute_loss` on the same inputs.
- After `FastGRPO.log`, the logged values of `entropy`, `kl` (when beta != 0), `clip_ratio/low_mean`,
  `clip_ratio/high_mean`, `clip_ratio/region_mean`, `clip_ratio/low_min`, `clip_ratio/high_max` equal TRL's.
- At start-up, if `inspect.getsource(GRPOTrainer._compute_loss)` does not hash to the pinned SHA-256, the trainer raises
  `RuntimeError` naming the method.

### 2.8 `objdel.train.vllm_sync`

```python
@dataclass(frozen=True)
class PushStats:
    n_tensors: int
    bytes: int
    seconds: float

def lora_touched_names(model) -> set[str]: ...
def push_weights(vg, model, *, full: bool) -> PushStats: ...
def install(trainer) -> bool: ...
```

- `lora_touched_names` returns the vLLM-side names (same naming as `lora_merged_params`) of every base weight that a
  LoRA layer changes. For the Qwen3-VL-2B LoRA config it has 196 names; no vision-tower name is in it.
- `push_weights` calls the vLLM model's `load_weights` exactly once per push.
  - `full=True`: the iterable holds every weight `lora_merged_params(model)` yields, with equal values.
  - `full=False`: the iterable holds exactly the `lora_touched_names(model)` weights, with values equal to
    `lora_merged_params`'s for those names.
  - When the engine uses sleep mode, `push_weights` wakes the weights first (`wake_up(tags=["weights"])`); a call with
    `full=False` while sleep mode is on raises `ValueError`.
  - After loading, `reset_prefix_cache()` is called once.
  - No parameter of the trainer model changes.
- `install(trainer)` replaces `trainer.vllm_generation.sync_weights` on that instance; returns False (and changes
  nothing) when there is no vLLM generation, no PEFT model, FSDP, or more than one process. The installed sync does a
  full push on its first call, after resume, every `full_push_every` calls (default 50), and on every call while sleep
  mode is on; otherwise a delta push.

### 2.9 `objdel.train.vllm_requests`

```python
def group_requests(prompts: list[dict], G: int) -> tuple[list[dict], list[int]]: ...
def expand_outputs(outputs: list, index: list[int], G: int) -> list: ...
```

- `group_requests` merges each run of G consecutive identical prompt dicts (equal `prompt_token_ids`, same image
  object) into one; raises `ValueError` if a run is not exactly G long or the length is not a multiple of G.
- `expand_outputs` returns `len(index) * G` single-completion outputs in the original row order: completion j of
  request k goes to row `k * G + j`.
- Used only when `train_rl --group-requests` is given (default off).

### 2.10 Measurement (`objdel.train.nvtx`, created by Phase A6)

```python
@contextlib.contextmanager
def count_syncs() -> Iterator[Counter]: ...
```

- Without CUDA: a no-op that yields a counter with value 0 and raises nothing.
- With CUDA: counts synchronizing operations reported by `torch.cuda.set_sync_debug_mode("warn")` inside the block
  and restores the previous mode on exit (also on an exception).
- `FastGRPO.timed(phase)` logs `perf/host_syncs` per phase when `--count-syncs` is on, and still logs
  `time/<phase>_s` as today.

## 3. Command-line flags (`experiments/phase2/train_rl.py` and `g_probe.py`)

| Flag | Values | Behaviour |
|---|---|---|
| `--gpu-pixels` | `auto` (default), `off` | `off` keeps the CPU processor path |
| `--encoder-once` | bool | section 2.5 path |
| `--host-meta` | bool | section 2.6 path |
| `--fast-sync` | bool | installs section 2.8 |
| `--full-push-every` | int, default 50 | full weight push period |
| `--group-requests` | bool, default off | section 2.9 |
| `--check-rope-every` | int, default 50; 0 = off | compares rope positions with HF on step 1 and every N-th step; raises on a difference |
| `--check-meta-every` | int, default 50; 0 = off | compares host metadata with the device derivation; raises on a difference |
| `--count-syncs` | bool | section 2.10 |
| `--profile-steps` | `A-B` | `torch.cuda.profiler.start()` before step A, `stop()` after step B |
| `--no-ckpt-layers` | int K, default 0 | turns gradient checkpointing off on K decoder layers |

GRPOConfig keys used (through `--grpo` JSON): `vllm_enable_sleep_mode`, `max_completion_length`, `beta`,
`per_device_train_batch_size`, `gradient_accumulation_steps`.

## 4. Constraints

- No new Python dependency. TRL stays 1.17.1.
- No test reads `data/` without `python -I`. No test needs a GPU unless marked `cuda`.
- Bit tests use `torch.equal` (or exact integer equality), not tolerances, except the C6 statistical gate.
- The trainer and vLLM must work with `vllm_importance_sampling_mode` set to `sequence_mask` and to `token_truncate`.
- The image processor's normalisation constants come from the processor object, never from literals in our code.
- Every private TRL name the new code overrides or patches (`_tokenize_prompts`, `_compute_loss`,
  `VLLMGeneration.sync_weights`, `split_pixel_values_by_grid`, `unsplit_pixel_values_by_grid`) is in the Phase A7
  start-up check.
- `make test` and `make lint` pass.

## 5. Acceptance gates

Equality (GPU box, real model; Mac tests cover the same with the tiny model where possible):

| Id | Gate |
|---|---|
| G-EQ1 | `patchify` output `torch.equal` to the HF CPU processor on >= 2,000 views (centred, shifted, rescaled) and all 32 images of 5 recorded steps |
| G-EQ2 | `rope_positions` `torch.equal` to HF `get_rope_index` on all 256 rows of >= 50 recorded batches; 0 `bad` rows on LIDC |
| G-EQ3 | `RowMeta`-based windows, orders, runs and packed indices equal the device-based ones on >= 50 recorded batches |
| G-EQ4 | Log-probs and entropies `torch.equal` between the `pixel_values` path and the `mm_encoder_outputs` path (reference pass and 16 update chunks of one batch); LoRA gradients of one micro-batch `torch.equal` |
| G-EQ5 | vLLM awake: after 3 delta pushes every vLLM parameter `torch.equal` to a full push of the same adapter state, untouched parameters unchanged; sleep mode: same after a sleep/wake cycle with full pushes |
| G-EQ6 | Step 1 of the fast run equals the baseline run (same seed, same host): reward, legal, IoU_tc, mean answer length, grad norm |

Speed and resources (15 steps per arm, mean of steps 3-15, R5 settings, same host as the baseline):

| Id | Gate |
|---|---|
| G-SP1 | Mean step time <= 19.5 s with vLLM awake, or <= 20.5 s with sleep mode (baseline 23.43 s) |
| G-SP2 | vLLM awake: `time/sync_s` <= 0.2 s and no wake/sleep phases; sleep mode: `time/sync_s` < 1.01 s |
| G-SP3 | NVTX `prepare` range <= 0.5 s per step |
| G-SP4 | Host synchronisations per step < 200 |
| G-SP5 | Host-to-device bytes per step < 50 MB (vLLM-internal copies excluded) |
| G-SP6 | GPU-active share of the step >= 70 % |
| G-MEM | Peak GPU memory <= 30 GB on every arm |
| G-ST1 | `--group-requests` arm: mean reward, legal and IoU_tc within 2 standard deviations of the baseline's step-to-step values; mean answer length within 5 %; >= 0.2 s per step saved |
| G-C | With the owner's outline limit and beta = 0 added: mean step time <= 18.0 s |
