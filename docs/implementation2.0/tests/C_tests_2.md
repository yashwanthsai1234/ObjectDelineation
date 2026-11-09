# Phase C test cases, set 2 (adversarial, property, performance, failure-mode)

Written from `phases/C.spec.md` plus the repo code (`src/objdel/train/fast.py`, `tests/test_image_memo.py`,
`tests/test_train_fast.py`). No design notes were read.

## Skills used
- nvidia-skills:nvidia-skill-finder: loaded; no catalog lookup done (no install wanted). Changed nothing.
- engineering:testing-strategy: pyramid split. Mac CPU unit/property tests first, `cuda` tests only where a GPU is the point.
- mattpocock-skills:tdd: tests assert behaviour at the spec's public seams, with expected values from HF/TRL (an
  independent oracle), never recomputed by the code under test.
- Not invoked (no effect on a test-case file): code-review, debug, system-design, architecture, documentation,
  huggingface-best.

## Conventions
- Oracle for all bit tests: the unmodified HF processor / HF `get_rope_index` / TRL 1.17.1 method / device-derived
  metadata. Compare with `torch.equal` or exact integer equality. No `allclose`.
- Mac processor: `Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)` as in `test_image_memo.py`.
  Tiny model: `tiny_model()` from `test_train_fast.py`.
- `cuda` marker is NOT registered in `pyproject.toml` (only `slow`, `data`). Add `"cuda: needs a GPU box"` and a
  `conftest` skip when `not torch.cuda.is_available()`, else `--strict-markers` or warnings hit.
- Hypothesis is not a listed dependency ("no new Python dependency"). Property tests below use seeded
  `numpy.random.default_rng` loops (>= 200 cases) with the seed in the failure message.
- Level key: U unit, P property, I integration, G GPU, E end-to-end (box session).

---------------------------------------------------------------------------------------------------------------

## A. patchify (2.1)

**T01 patchify_equals_hf_random_sizes** (P, Mac)
- Setup: processor above; rng seed 0; 200 cases. Sizes drawn from multiples of 32 (h,w in 32..256) so
  smart_resize is identity (factor 32 = patch 16 x merge 2). N in 1..5 images per call, content uint8 uniform.
- Input: list of PIL RGB, and the same as `(N,H,W,3)` uint8 tensor only when sizes are equal.
- Expected: `torch.equal` on pixels and grid vs `ip(list, return_tensors="pt")`, backend `"torch"` and `"hf"`.
  grid dtype int64, device CPU.
- Why: catches wrong patch ordering for non-square or non-512 shapes (e.g. swapped h/w in the reshape/permute) that
  a single 512x512 test misses.

**T02 patchify_extreme_pixel_values** (U)
- Input: 64x64 images all 0, all 255, a checkerboard 0/255, a ramp 0..255, one with 1-pixel value 255 at corner
  (0,0) and another at (H-1,W-1) on zero background.
- Expected: `torch.equal` to HF for each.
- Why: float rounding of `(x/255 - mean)/std`. Order of ops (rescale then normalise; fused vs two-step) differs in
  the last bit exactly at values like 0, 255 and 1/255 steps. Single-pixel images catch off-by-one in the corner patch.

**T03 patchify_all_256_gray_levels** (U)
- Input: a 16x16 patch tiled to 32x32, one image per level v in 0..255 (256 images, one call). Also the 3 channels
  set to different levels (r=v, g=255-v, b=(v*7)%256).
- Expected: bit-equal to HF.
- Why: exhaustive over the uint8 domain, so no input value can hit a rounding tie the sampled tests skipped.
  Catches channel-order (RGB vs BGR) bugs because channels differ.

**T04 patchify_grey_equals_rgb_convert** (U)
- Input: 64x64 `L` PIL image, `(H,W)` uint8 numpy, `(1,H,W)` uint8 tensor, vs `Image.convert("RGB")` fed to HF.
- Expected: all three bit-equal to HF on the converted image. Also `(N,H,W)` with N=3 equals the HF list.
- Why: ambiguity trap. A `(3,H,W)` tensor with H=3 could be read as `(H,W,3)`/channels. Test `(3, 64, 64)` grey batch of
  3 and a `(64, 64, 3)` single RGB: each must go to the correct reading (spec: `(N,H,W,3)` or `(N,H,W)`; a
  `(H,W,3)` numpy array is one RGB image). Also `(64,64,3)` as torch tensor: spec lists only 4-D/3-D tensors, so
  3-D tensor means `(N,H,W)` grey; with last dim 3 and H=W=64 that reads as N=64 grey images of 64x3, which
  raises `UnsupportedImage` (too small? 64x3 not multiple of 32) rather than silently mis-batching. Pin whichever
  behaviour the spec resolves; flag spec ambiguity to the owner.

**T05 patchify_unsupported_raises** (U)
- Input table, each expects `UnsupportedImage` (and `isinstance(e, ValueError)`):
  500x500 RGB; 511x512; 512x513; 1x1 RGB; 0-size `(0,64,64,3)` (see T06); RGBA PIL; `P` PIL; `I;16` PIL; `1` PIL;
  `F` PIL; numpy `(64,64,4)`; numpy `(64,64,2)`; float32 tensor; int64 tensor; uint8 `(2,64,64,1)`.
  Also processor with `patch_size=14`, with `merge_size=1`, with `temporal_patch_size=1`.
  Expect also: one bad image among 31 good ones in a list raises for the whole call (no partial output).
- Why: silent acceptance of a wrong mode would produce plausible-looking tensors that differ from HF. Float/int64
  tensors must not be silently cast.

**T06 patchify_empty_and_single_pixel_patch** (U)
- Input: `(0,64,64,3)` tensor; empty list; one 32x32 image (smallest valid: 2x2 patches = 1 merged token, grid
  `[1,2,2]`, 4 rows x 1536).
- Expected: smallest valid bit-equal to HF; empty input behaves like HF for the same call (HF raises) or raises
  `UnsupportedImage`/`ValueError` — never returns a `(0,1536)` tensor with a `(0,3)` grid unless HF does. Record
  HF's actual behaviour in the test with a `pytest.raises` mirror.
- Why: boundary N=0 breaks `torch.cat` / reshape(-1, ...) code paths.

**T07 patchify_max_size** (U + G)
- Input: 2048x2048 RGB (smart_resize identity if inside processor max_pixels; if processor max_pixels forces a resize
  the spec requires `UnsupportedImage`). Read the processor's `max_pixels` and derive the largest multiple-of-32
  square that is unchanged, plus 32 more.
- Expected: largest unchanged size bit-equal; one step above (resize required) raises `UnsupportedImage`.
- Why: boundary of `smart_resize` is computed from the processor, not a literal; catches hard-coded 512.

**T08 patchify_does_not_mutate_input** (U)
- Input: uint8 tensor and numpy array, copy before; call; compare after. Also a non-contiguous view
  (`x[:, ::1, :, :].permute(...)`/ a transposed numpy array `np.asfortranarray`) gives the same output as the
  contiguous copy.
- Expected: input bytes unchanged; non-contiguous result bit-equal to contiguous.
- Why: in-place normalise on the caller's buffer would corrupt the shared ImageMemo images (which are reused by
  later passes); strided input exposes `.view` vs `.reshape` bugs.

**T09 patchify_backends_and_names** (U)
- `backend="auto"|"hf"|"torch"` equal to HF on CPU; `"k3"` on CPU either equals HF or raises (spec: every backend
  that returns without raising must be bit-equal). `backend="cuda"`, `""`, `"HF"`, `None` raise `ValueError`
  (not `UnsupportedImage` — check `type(e) is ValueError` is NOT required, but message names the backend).
- Why: `UnsupportedImage` is a `ValueError` subclass, so a lazy implementation that raises UnsupportedImage for a bad
  name would make `ImageMemo` fall back silently. Test that for a bad name the ImageMemo path propagates (see T20).

**T10 patchify_out_dtype_and_device** (U + G)
- `out_dtype=torch.bfloat16`, `float16`, `float32`; `device="cpu"`; on box `device="cuda"`, `"cuda:0"`.
- Expected: dtype/device as asked; grid always CPU int64; bf16 result equals `hf_fp32.to(bf16)` bit for bit (the
  cast happens last, once).
- Why: normalising in bf16 gives different bits than casting fp32 output.

**T11 patchify_cpu_vs_cuda_bit_equal** (G)
- Setup: box. 2,000 random views (centred, shifted by random offsets, rescaled with bilinear) built from 20 real
  recorded images (or noise if data absent); batches of 1, 8, 64, 256.
- Expected: `patchify(..., device="cuda")` moved to CPU `torch.equal` `patchify(..., device="cpu")` and HF. Run each
  backend. Repeat 3 times: identical bits each time (no atomics / nondeterminism).
- Why: GPU fused multiply-add (`fma`) contraction changes the last bit versus CPU. Gate G-EQ1.

**T12 patchify_no_sync** (G)
- `with count_syncs() as c: patchify(cuda_uint8_batch, ip, device="cuda")` on a batch already on device.
- Expected: `c == 0` after warm-up call (warm-up excluded). Also assert no D2H copy of pixel_values
  (`torch.cuda.memory` trace) and grid is built on CPU from shapes only.
- Why: gate G-EQ1 "0 syncs"; a hidden `.item()` for the grid or `.tolist()` fails it.

**T13 patchify_normalisation_constants_from_processor** (U)
- Build a processor with `image_mean=[0.1,0.2,0.3]`, `image_std=[0.5,0.6,0.7]`, and `rescale_factor=1/200`; compare to
  HF with the same processor. Then a processor with `do_normalize=False` and one with `do_rescale=False`.
- Expected: bit-equal to HF; if the implementation does not support the toggles it raises `UnsupportedImage`,
  never returns standard-constant output.
- Why: spec "never literals in our code". Grep test: `inspect.getsource(gpu_pixels)` contains none of
  `0.5, 0.5, 0.5`, `0.48145466`, `0.26862954`, `255`-as-divisor-literal other than via `rescale_factor`.

**T14 patchify_pad_tail_for_batch_of_odd_temporal** (U)
- Single image: Qwen duplicates the frame to fill `temporal_patch_size=2`. Check both temporal slices of every patch
  row are identical: `pv.view(-1, 3, 2, 16, 16)[:,:,0] == [:,:,1]`.
- Why: independent structural invariant (not derived from HF). Catches a temporal-axis reshape bug even if both
  sides share a wrong helper.

---------------------------------------------------------------------------------------------------------------

## B. ImageMemo (2.2)

**T20 memo_counters_exact** (U)
- Setup: 32 distinct 64x64 images, row list `[img[i//8] for i in range(256)]`, fresh memo, `clear()`.
- Expected: after call: `misses == 32`, `hits == 224`. Call again: `hits == 480`, `misses == 32`. `clear()`; `prefill(32
  images)`: `misses` +32; then call: `hits += 256`, `misses` unchanged. `prefill` twice: second adds 0.
- Why: spec numbers are exact. Catches counting per distinct image instead of per row, and double-count in prefill.

**T21 memo_per_image_true_values** (U)
- Rows `[a,a,b,c,a,b,c,a]` (c another size). `per_image=True`.
- Expected: `pixel_values == hf([a,b,c]).pixel_values` (first appearance order, `torch.equal`); `image_grid_thw` shape
  `(8,3)` equal to `hf(rows).image_grid_thw`. Total patch rows = sum of distinct only.
- Adversarial: rows `[c,a,c,b,a]` -> order `c,a,b`. Two equal-content but distinct objects `a1`, `a2` (same pixels)
  -> two blocks (keyed by object, not content). A single image repeated 256 times -> one block.

**T22 memo_object_identity_hazard** (U)
- Mutate an image in place (`img.paste(...)` or `np_array[...] = 0`) between two calls without `clear()`.
- Expected (pin the contract): memo returns the OLD patches (documented "keyed by object"); with `clear()` returns the
  new. If implemented otherwise, the test records it. Another case: `del img` then allocate a new image that reuses
  the same `id()`; because the memo holds a reference the id cannot be reused — assert `memo.memo` keeps the image alive
  (`weakref` stays alive until `clear()`).
- Why: id() reuse would serve wrong pixels to a different image: silent training corruption.

**T23 memo_kwargs_keying** (U)
- Same image with `return_tensors="pt"` then `return_tensors=None`, `do_resize=False`, `size=...`.
- Expected: each result equals the wrapped processor for that exact kwargs; never a cross-served result of other kwargs.

**T24 memo_cuda_placement** (G; Mac part: monkeypatch `torch.cuda.is_available` -> False gives CPU tensors)
- On box: `pixel_values.is_cuda`, `image_grid_thw.device.type == "cpu"`, both for per_image True and False and for the
  fallback path (T25). Values `torch.equal` after `.cpu()` to the wrapped processor.

**T25 memo_unsupported_falls_back_whole_call** (U)
- Rows: 31 valid 64x64 images and one 500x500 (needs resize).
- Expected: result `torch.equal` wrapped processor for the whole call (all rows, order kept); no exception; the valid
  images are not left half-memoised in a state that changes later results (a following call with only valid images
  equals HF). Counters stay consistent: hits+misses == rows seen on the GPU path only (pin the choice).
- Also: a processor that raises `RuntimeError` (not `UnsupportedImage`) inside the wrapped call propagates — only
  `UnsupportedImage` is caught.
- Why: the catch-all `except ValueError` trap (UnsupportedImage subclasses ValueError, `ValueError("bad backend")` must NOT trigger fallback).

**T26 memo_property_random_row_lists** (P)
- 300 seeded cases: pool of 1..12 images of 3 sizes; row list length 1..64 drawn with repetition and random order;
  random `per_image` flag, random `prefill` before or not, random `clear()` between calls, 1..3 calls.
- Invariants: (i) `per_image=False` equals HF on rows; (ii) `per_image=True` equals HF on distinct list; (iii) wrapped
  processor total image count over all calls <= number of distinct objects per clear-epoch; (iv) `hits+misses ==
  total rows looked up`; (v) `sum(grid.prod(-1)) == pixel_values.size(0)` when `per_image=False`; when
  `per_image=True`, `pixel_values.size(0) == sum over distinct`.
- Why: combinatorial interplay of prefill, clear, mixed sizes.

**T27 memo_thread_safety_and_reentrancy** (U)
- Two threads (reward prefetch uses threads in this repo) call the memo on overlapping image sets, 200 iterations.
- Expected: no exception, results equal HF, each distinct image processed at most once per epoch (counting wrapper
  with a lock). If the implementation is not thread-safe the spec must say so; this test pins that decision.

**T28 memo_deepcopy_and_getattr_preserved** (U)
- `copy.deepcopy(memo)` returns the plain processor (regression of existing behaviour); `memo.patch_size == 16`;
  `ImageMemo.__getattr__` on a bare instance (`object.__new__`) raises `AttributeError`, no infinite recursion; new
  attributes `per_image`, `hits`, `misses` must survive `copy.copy` and are not forwarded to the inner processor.
- Why: new attributes added to `__init__` after `__getattr__` can recurse during pickling.

**T29 memo_memory_bound** (G)
- 256 rows over 32 images, 512x512: GPU memory held by the memo after the call <= 32 x 1024 x 1536 x 4 B = 201 MB
  (+10 %); with `per_image=True` the result holds exactly that, with False the result holds up to 256 x ... only if
  a copy is made (spec 2.7 forbids the 256x1024x1536 tensor on the encoder-once path: assert
  `max_memory_allocated` delta < 400 MB on that path).
- Why: 256 x 1024 x 1536 x 4 B = 1.61 GB; gate for 2.7 last bullet.

---------------------------------------------------------------------------------------------------------------

## C. Start-up pixel check (2.3)

**T30 startup_selects_first_equal_backend** (U with fakes)
- Monkeypatch `patchify` with a fake: `"hf"` returns output with one flipped bit in image 7 of 8; `"torch"` exact.
- Expected: selected backend `"torch"`; log line names it. If `"hf"` exact: selected `"hf"` and `"torch"` not run.
- Adversarial: fake `"hf"` wrong only on the GREY image (RIGA) or only on the digit image: still rejected (all 8 must
  match). Fake that raises `UnsupportedImage` on one of the 8 -> that backend counts as failed.

**T31 startup_no_backend_raises_unless_flag** (U)
- Both fakes wrong. Expected: `RuntimeError` at start-up; with `--allow-cpu-pixels`: continues, CPU processor used,
  log contains `gpu_pixels/backend=cpu`. With `--gpu-pixels off`: no check, no error.
- Failure mode: flag parse `--allow-cpu-pixels` default False: assert via argparse defaults.

**T32 startup_check_image_set** (U)
- Expected: 8 images; contains each training-view kind present (centred, shifted, rescaled), a grey image iff RIGA
  rows exist, exactly one generated digit image. No-CUDA machine: the check is skipped (no error).

---------------------------------------------------------------------------------------------------------------

## D. Rope positions (2.4)

Oracle: `base.get_rope_index(input_ids, image_grid_thw, None, attention_mask)` on the tiny model (and on the real
processor's token layout with Lp from a rendered prompt). Build prompts as `[START, IMG*n, END, text...]` with
`mm_token_type_ids` 1 on image tokens.

**T40 rope_equal_all_lengths_and_padding** (P)
- 300 seeded cases: B in 1..16, Lp from template (e.g. 3 + 16 image tokens + 8 text), left pad 0..40 per row,
  answer length from {0, 1, 2, 7, 131, 1536} plus random 0..200, optional masked holes inside the answer region
  (attention 0 for random answer columns; the spec says holes "masked").
- Expected: `positions[:, i]` `torch.equal` `get_rope_index(...)[0][:, i]` for every good row; masked columns are 0
  (check `positions.masked_select(~mask) == 0`); `bad` all False.
- Why: Rope value after a hole: HF counts only attended tokens (cumulative of mask) — a naive `arange` would shift
  everything after the hole.

**T41 rope_all_masked_prefix_rows** (U)
- Rows with 0 real tokens (all mask 0) or fewer than Lp real tokens.
- Expected: `bad[i] == True`; no exception, no NaN/out-of-range gather; positions for other rows unaffected.

**T42 rope_bad_detection_each_cause** (U)
- Case table (each flips exactly one thing, others in same batch stay good):
  1. one id in the template region changed; 2. one token type changed but ids same; 3. an image-type token in the
  answer region (`mm_token_type_ids=1` after Lp); 4. a different image count (second `START..END` block in the
  answer); 5. real tokens = Lp-1; 6. the row's prefix is the template shifted by one (misaligned left pad: pad count
  differs but real tokens identical -> must be GOOD); 7. a prefix equal to the template only after
  skipping a masked hole (hole inside template region) -> BAD (spec: first Lp real tokens; holes within the prompt
  region mean the real-token sequence still equals the template; define: real tokens are taken by mask, so a hole
  inside the prompt does NOT make it bad but HF's positions for it differ in structure — compare to HF; whatever
  `bad` says, if `bad[i]` False then `positions` must equal HF.)
- Master invariant for all cases: `not bad[i]` implies `positions[:,i] == hf[:,i]`. Never a false negative.
- Why: the safety property of the whole optimisation; a false-good row trains on wrong positions silently.

**T43 rope_no_false_good_fuzz** (P)
- 500 random mutations of valid rows (flip random id / type / mask bit anywhere, truncate, insert image token,
  swap two tokens). Compute HF positions per row. Assert for every row where `~bad`: equal to HF.
- Also assert at least 5 % of mutations are flagged bad (the test has teeth) and that mutations in the answer text
  region that stay text are NOT flagged (no pointless fallback; efficiency).

**T44 rope_no_host_sync** (U + G)
- Mac: monkeypatch `torch.Tensor.item`, `.tolist`, `.cpu`, `.numpy`, `torch.Tensor.__bool__`, `torch.nonzero` to
  raise during `rope_positions` and run; boolean-mask indexing check via `torch.Tensor.__getitem__` wrapper that
  raises when the index is a `bool` tensor. Expected: no raise. 
- Box: `count_syncs()` == 0 around `rope_positions` on CUDA tensors (after warm-up).
- Why: the stated no-sync requirement. A masked `index` silently syncs.

**T45 rope_dtype_device_shape** (U)
- Output `positions` is int64, shape `(3,B,W)`, on input device; `bad` bool `(B,)`. Input as int32 ids or attention mask
  given as bool/int64: result same (HF may take int64 only; pin).
- B=1, W=Lp (no answer); B=256, W=2000 (max size) completes < 50 ms CPU.

**T46 rope_template_validation** (U)
- `rope_template` with `input_ids` of shape `(2, Lp)` raises (ValueError); with no image tokens works (text-only
  prompt) and `rope_positions` equals HF; with `image_grid_thw` row count 2 but one image in prompt -> raises or equals
  HF's error. `RopeTemplate` frozen: assigning `tmpl.end = 0` raises `FrozenInstanceError`.
- `tmpl.pos.shape == (3, Lp)`; `tmpl.end` equals `pos.max() + 1` (independent structural check of HF's meaning).

**T47 rope_real_stack_positions** (G, box, gate G-EQ2)
- 256 rows x 50 recorded batches (15 from box). Expected: equal; `bad.sum() == 0` on LIDC. Includes batches where an
  answer reached 1536 tokens (truncated) and 0-token answers.

---------------------------------------------------------------------------------------------------------------

## E. Prompt ids and PromptMemo (2.5, 2.6)

**T50 tokenize_prompts_equals_trl** (I, tiny processor + fake vLLM flag)
- Use a real Qwen3-VL tokenizer/processor from the local HF cache if present (skip otherwise) with tiny model.
- Rows: 32 prompts x 8 identical copies, plus a variant with 2 image sizes and 2 prompt texts.
- Expected: `prompt_ids == TRL._tokenize_prompts(prompts)[0]` (list equality), `images is`-equal element-wise (`a is b`)
  to `[1]`, third item `{}`; `image_processor` call count == 0 (counting wrapper). Non-vLLM path: whole tuple equals TRL.

**T51 tokenize_at_most_once_per_distinct_key** (U)
- Count TRL tokenizer-call invocations over: two calls (same 256 rows); a call after `clear()`; same text with a
  different image size (must re-run), same image size and different text (must re-run).
- Expected: calls == number of distinct `(text, (H,W))` seen in the run; the key must include image size. Adversarial:
  the same text with the same size but different image CONTENT -> same ids (spec: key is text and size) — pin this
  and confirm token count depends only on size (true for Qwen pad count).
- Does the memo ever go stale across `clear()`? Spec: "per run" so text memo survives `clear()`; assert so.

**T52 promptmemo_install_family_guard** (U)
- Non-Qwen processor (a fake `CLIPProcessor`-like with `image_processor`) -> returns None and `processing_class.__class__`,
  `__dict__` identical to before (deep compare of `vars`).
- Qwen processor: second `install` call returns the same memo or an equal wrapper, not double-wrapped (calls still once).

**T53 promptmemo_isinstance_and_attrs** (U)
- After install: `isinstance(p, orig_cls)` and for every class in `type(p).__mro__` before install; `p.tokenizer is
  orig.tokenizer`; `p.chat_template`; `p.save_pretrained(tmp)` works and writes the same files as before; `copy.deepcopy(p)`
  works; `pickle.dumps(p)` same as before when it was picklable.

**T54 promptmemo_values_equal_mixed_batches** (P)
- 150 seeded cases: rows with random (text from 4 templates, image size from 3); `padding=True`, `return_tensors="pt"`.
- Expected: every key (`input_ids, attention_mask, mm_token_type_ids, pixel_values, image_grid_thw`, whatever the
  original returns) `torch.equal` to the pre-install result (keep an unmodified second processor as oracle). Includes
  left vs right padding side as configured, and batch of 1.
- Why: assembling per-row results with padding by hand is where off-by-one pad sides appear; mixed sizes break
  stacking assumptions.

**T55 promptmemo_other_signatures_untouched** (U)
- Calls: text only (no images); `images` only; with `videos=`; `padding="max_length", max_length=...`; with
  `return_tensors=None`; extra kwargs `add_special_tokens=False`; `processor.apply_chat_template(...)`. Each result equals
  the unmodified processor's. Also the call counter on the inner processor increments (the memo is bypassed, not
  used).

**T56 promptmemo_runs_original_once_per_key_count** (U)
- 256 rows, 1 text, 1 size: original processor text/ids call counts == 1 for the run (second `clear()`+call also 0
  new). 2 texts x 2 sizes: 4.
- Failure mode: an exception inside the original call must not poison the memo: after the error, the same call
  retried reruns the original and succeeds (assert call count 2, no cached exception).

---------------------------------------------------------------------------------------------------------------

## F. Vision encoder once per generation batch (2.7)

**T60 encoder_once_equal_logps_matrix** (I, tiny model, CPU fp32)
- Matrix: `shared_prefix` in {False, True} x attn in {sdpa, eager} x chunk in {8, 4}. Batch 16 rows = 2 prompts x 8
  answers (and a 4-image variant), answer lengths varied incl. 0 and 1.
- Expected: per-token log-probs and entropies `torch.equal` between `pixel_values` path and
  `mm_encoder_outputs` path; tower call count (hook on `visual.forward`) == 1 per generation batch for the new path
  (and covers reference + every update chunk), versus > 1 for the old path (control, proves the counter works).
- Why: gate G-EQ4. The control prevents a counter that always reads 1.

**T61 encoder_once_gradients_equal** (I)
- One micro-batch: backward with LoRA on the tiny model; every LoRA param grad `torch.equal` between paths. Also the
  frozen tower params have `grad is None` on both and no autograd graph reaches the tower output on the new path
  (`mm_encoder_outputs["image"].requires_grad is False`).

**T62 encoder_once_batch_has_no_pixels_key** (U)
- `_generate_and_score_completions` output: `"pixel_values" not in batch`, `"host_meta" in batch`, `image_grid_thw`
  present; `host_meta` is a list of length 256 (row aligned); every forward spy sees `mm_encoder_outputs` and no
  `pixel_values` kwarg.

**T63 encoder_once_shuffle_alignment** (P)
- 100 seeded permutations: TRL shuffle (`shuffle_sequence_dict`), `by_length`, `by_group_length`, then split into
  chunks of 8 and 16 (and the last short chunk). Expected: for each micro-batch row the image id travels with it:
  `image_ids[chunk] == original_image_id[perm][chunk]` and the encoder outputs gathered for that row `torch.equal` the
  tower output of that row's own image computed standalone.
- Adversarial: groups where two prompts share one image object (`image_id` equal across groups); a micro-batch that
  mixes images of different grid sizes (`(1,32,32)` and `(1,16,48)`) so the per-image output rows have unequal length
  -> gather must use per-image offsets, not equal-size reshape. Expected equal to the pixel path.

**T64 encoder_once_missing_ids_raises** (U)
- Remove the per-row image id field from a micro-batch with the path on. Expected: forward raises (specific exception
  type pinned: `RuntimeError`/`KeyError` containing "image"), and the tower call count stays 0 afterwards (no silent
  re-encode). With the path off and ids missing: normal forward works.
- Also: ids of wrong length (B-1), ids out of range (>= distinct count), negative ids -> raise, not a modulo wrap.

**T65 encoder_once_stale_cache_across_generation_batches** (I)
- Two generation batches with different images; between them weights unchanged. Batch 2 outputs must come from batch
  2's images: `torch.equal` to the pixel path on batch 2. Spoof: reuse the same image objects with changed content
  (`clear()` is called each batch) -> outputs updated. A test that skips `clear()` documents the failure.
- Another hazard: LoRA update between batches does not change the tower (frozen, assert `torch.equal` of tower
  outputs before and after an optimizer step on LoRA params); if the tower had a trainable param (config error),
  startup must raise/disable the path. Test with `requires_grad=True` on one tower param: path refuses.

**T66 encoder_once_memory_gate** (G)
- Peak memory of the generation batch with the path on: no tensor of shape `(256*1024, 1536)` is alive
  (scan `torch.cuda.memory_snapshot()` or `gc` for numel >= 256*1024*1536 -> none). Peak <= 30 GB (G-MEM).

**T67 encoder_once_token_count_mismatch** (U)
- `image_grid_thw` says an image yields 256 merged tokens but the prompt has 255 image tokens (tokenizer drift).
  Expected: raise a clear error (HF raises "Image features and image tokens do not match"); not silently truncated.

---------------------------------------------------------------------------------------------------------------

## G. Row metadata (2.8)

Oracle: the existing device functions `window`, `by_length`, `by_group_length`, and TRL masks.

**T70 rowmeta_completion_len_rules** (U)
- G=8, `max_completion_length=16`. Rows: ended with EOS at 5 (len 5 incl. EOS: pin count of tokens given); length 0
  (EOS first token / empty); length 16 with no EOS and `mask_truncated=True` -> 0; same with False -> 16; length 16
  ending EOS exactly at position 16 -> 16 (not truncated!); length 15 no EOS -> 15 (not hit limit).
- Expected: exact int values, dtype int32 for all four arrays, shape `(B,)`.
- Why: the truncated-vs-EOS-at-limit boundary decides whether a row gets masked from the loss.

**T71 rowmeta_group_and_image_id** (U + P)
- `group == arange(B)//G`. image ids: rows sharing object -> same id, distinct objects -> distinct, ids dense from
  0 in first-appearance order (pin) even if the same image object appears in non-adjacent groups. B not divisible by
  G raises `ValueError`. B=0 -> empty arrays (no crash).

**T72 rowmeta_window_equals_device_window** (P)
- 500 seeded cases: random left-padded prompts, answers lengths 0..L, `logits_to_keep`, buckets in
  {BUCKETS, FINE, custom}, subsets `rows` (contiguous chunks and random subsets), including all rows with zero-length
  answers and one row with `completion_len == logits_to_keep`.
- Expected: `window_from_meta(...) == window(attention_mask[rows], ...)` exactly (tuple of ints). Include prompt-only
  rows where all attention is on the prompt and the `keep >= 1` clamp applies.
- Why: gate G-EQ3; the `max(longest,1)` clamp and the bucket-overflow fallback `W` are the likely mismatches.

**T73 rowmeta_orders_equal_device_orders** (P)
- 500 cases incl. ties (many equal lengths), all-equal lengths, strictly decreasing, B=8 (one group), B=256.
- `length_order(meta)` equals `torch.sort(n, descending=True, stable=True).indices`; `group_order(meta, G)` equals
  the `order` inside `by_group_length`. Check that stable ties keep original order (all-equal lengths -> identity).
- Invariants: `group_order` is a permutation; groups contiguous: `order.reshape(-1,G) // G` has constant rows; each
  group's internal order is ascending.
- Adversarial: the TRL-masked row lengths: `n` uses real masks, not `completion_len` after mask_truncated? Spec says
  `completion_len` is 0 for truncated rows when mask_truncated, yet `by_length` sorts by `completion_mask.sum`
  (which is 0 for such rows too in TRL: the mask is zeroed). Test with a truncated row: equal both ways. If they
  differ, report spec bug.

**T74 rowmeta_packed_indices_equal** (P)
- `seg`/`idx` built from RowMeta vs from device tensors for the same chunk, 200 cases with packed chunks
  (`packing_wins` both true and false), groups of 8 with varied lengths. `torch.equal` / `np.array_equal`.

**T75 rowmeta_reads_no_device_tensor** (U + G)
- Pass numpy / list inputs only; monkeypatch `torch.Tensor.__array__`, `.item`, `.tolist`, `.cpu` to raise. Functions
  work. On box: `count_syncs()` == 0 across all four functions with host inputs. A tensor argument to `row_meta`
  (a CUDA tensor in a prompt list) raises `TypeError` or is converted without sync — pin one.

**T76 host_meta_batch_field_alignment** (I)
- After TRL shuffle + split of the batch dict, `host_meta` (python list field) is reordered with the rows
  (`move` handles lists) -> for each chunk `host_meta[i].prompt_len == prompt_mask[i].sum()`. Scalar `0-d` and `None` fields
  untouched. A `host_meta` of the wrong length raises at split time.

**T77 check_meta_every_raises_on_difference** (U)
- Inject a corrupted `completion_len` by 1 in one row; with `--check-meta-every 1` the step raises (message names the
  field and row); with 0 no check; with 50 only steps 1, 51, ... check (count calls over 101 steps == 3).
- Same for `--check-rope-every` (flip one position). Default value 50 asserted via argparse.

---------------------------------------------------------------------------------------------------------------

## H. Loss metrics (2.9, only if adopted)

**T80 loss_equal_to_trl_and_zero_item** (I)
- Tiny model batch, `loss_type="dr_grpo"`. Count `torch.Tensor.item` calls during `_compute_loss` via monkeypatch ->
  0. Loss `torch.equal` to `GRPOTrainer._compute_loss` on the same inputs; also for `beta=0` and `beta=1e-3`, both
  `vllm_importance_sampling_mode` values (`sequence_mask`, `token_truncate`), and a batch where all advantages are 0.

**T81 logged_metrics_equal_after_log** (I)
- Exact equality of `entropy`, `kl`, `clip_ratio/{low_mean,high_mean,region_mean,low_min,high_max}` against TRL on the
  same two micro-batches; `kl` absent when beta == 0. Edge: all rows masked (zero tokens) -> TRL's value (0 or nan)
  reproduced exactly (`torch.equal` with nan handled via `isnan` both).

**T82 loss_source_hash_guard** (U)
- Monkeypatch `inspect.getsource` to return the source plus a space -> start-up raises `RuntimeError` whose message
  contains `_compute_loss`. Unmodified source passes (pins the SHA to trl 1.17.1).
- Same guard for `_tokenize_prompts` and `VLLMGeneration.sync_weights` (A7 list).

---------------------------------------------------------------------------------------------------------------

## I. vllm_sync (2.10)

Use a fake `vg` (records calls in order) and the tiny PEFT model.

**T90 push_weights_exact_payload** (U)
- Expected: `load_weights` called once; iterable materialised equals `list(lora_merged_params(model))` names and
  `torch.equal` values; `PushStats.n_tensors == len`, `bytes == sum(numel*element_size)`; `seconds` measured only
  around `load_weights` (fake load sleeps 0.05 s, the generator creation sleeps 0.3 s: `0.05 <= seconds < 0.15`).
- Hazard: the iterable is a generator consumed lazily by `load_weights`; stats must still count after consumption
  (n_tensors correct although computed during the load) and a fake `load_weights` that consumes it twice must not
  see an empty second pass (spec "iterable" - pin generator vs list).

**T91 push_order_sleep_mode** (U)
- Calls order for sleep mode: `wake_up(tags=["weights"])`, `load_weights`, `reset_prefix_cache` (once). Awake: no
  `wake_up`, still `load_weights` then `reset_prefix_cache`. Exception inside `load_weights`: propagates; the engine
  is not left with an unreset prefix cache? (pin: `reset_prefix_cache` still called or the error documented).

**T92 push_does_not_change_trainer_model** (U)
- Hash all parameters and buffers (name, dtype, `torch.equal` copy) before/after: unchanged, including LoRA A/B (no merge-
  and-unmerge rounding drift), `requires_grad` flags, `training` flag, and adapter enabled state. 3 pushes in a row
  give bit-identical payloads.
- Adversarial: model in `train()` mode before push stays in train mode; pushing under `torch.no_grad` not required.

**T93 install_returns_false_cases** (U)
- Table: no `vllm_generation`; model not PEFT; FSDP flagged; `world_size == 2`; each returns False and
  `trainer.vllm_generation.sync_weights` is the same object as before (`is`). True case: replaced on the instance only,
  `type(vg).sync_weights` unchanged (a second trainer's vg is not affected). Installing twice does not stack wrappers.
- Installed sync called with the TRL signature (no args) performs a full push.

**T94 push_equals_trl_per_parameter_push** (G, gate G-EQ5)
- Real vLLM 0.32: 3 adapter states (random LoRA perturbations); compare every vLLM param `torch.equal` after fast push vs
  TRL's own per-parameter push; sleep variant after sleep/wake. Include a parameter that vLLM fuses (qkv, gate_up) and
  a tied embedding/lm_head. Time `PushStats.seconds <= 0.2 s` (G-SP2).

**T95 push_repeatable_no_memory_growth** (G)
- 50 pushes: `torch.cuda.memory_allocated` after each push within 1 MB of the first (no leak from merged copies);
  `seconds` stable.

---------------------------------------------------------------------------------------------------------------

## J. group_requests (2.11)

**T100 group_requests_valid_and_invalid** (U)
- `[p]*8 + [q]*8` with `G=8` -> `[p, q]` (same objects, `is`). Invalid -> `ValueError`: length 15; runs 8 + 7 + 1;
  run of 9; two different prompt dicts interleaved `p q p q...`; same ids but different image object (identity check:
  equal `prompt_token_ids`, different `multi_modal_data["image"]` objects -> ValueError); equal image object but
  different ids -> ValueError; empty list -> `[]` (pin) ; `G=1` returns the list; `G=0` raises.
- Adversarial: two adjacent groups with IDENTICAL prompt dicts (same prompt sampled twice, 16 equal rows, G=8) ->
  valid as two runs? Ambiguous: expected `[p, p]` (runs of exactly 8 are by position). Pin and flag.

**T101 group_requests_order_restoration** (I, fake vLLM)
- Fake engine returns for request k, n=G completions tagged `(k, j)`. Expected trainer completions list, in order,
  `[(0,0),(0,1),...,(31,7)]`, length 256, completion j of request k is row `k*G+j`. vLLM received 32 requests, each
  `n == G` (sampling params). Default off: 256 requests, n=1. Property over G in {1,2,4,8}, K in {1,3,32}.
- Failure: engine returns fewer than G completions for a request -> raise, not shift rows.

---------------------------------------------------------------------------------------------------------------

## K. count_syncs (2.12) and gpu_active (2.13)

**T110 count_syncs_no_cuda** (U, Mac)
- Yields a Counter-like with value 0 (`c["syncs"] == 0` or `sum(c.values()) == 0`; pin the key), body runs, exception
  inside the body propagates and nothing else raised. Nested use works.

**T111 count_syncs_exact** (G)
- 1,000 `x.item()` in a loop on one line -> count == 1000 (not 1: warnings dedupe by line unless the implementation
  hooks differently — this is the repeat trap). 1 `x.cpu()` of a CUDA tensor -> 1; `torch.cuda.synchronize()` -> 1;
  async kernel only -> 0; `x.tolist()` -> 1.
- Mode restore: `torch.cuda.get_sync_debug_mode()` before == after, set `"error"` before to check restoration to
  "error" (and exception in body still restores). Re-entrancy: nested blocks count independently (outer >= inner).
- Warnings filter: user `warnings.simplefilter("ignore")` outside must not suppress counting.

**T112 timed_logs_syncs_per_phase** (I)
- `--count-syncs`: `perf/host_syncs` logged per phase; `time/<phase>_s` still logged; `generate` syncs only in
  `perf/host_syncs_vllm`, absent from `perf/host_syncs` (assert key sets). Flag off: no `perf/` keys and no mode
  change (`get_sync_debug_mode() == 0`). `time/wake_weights_s` and `time/wake_kv_s` both logged under sleep mode.

**T113 gpu_active_share_exact** (U, Mac)
- Synthetic intervals (ns): stream A `[0,40]`, `[30,60]`; stream B `[50,80]`, `[90,100]`; memcpy `[100,110]`; range
  `[0,200]`. Union = `[0,80] + [90,110]` = 100 -> share 0.5 exactly. Range `[10,50]` -> union `[10,50]` -> 1.0.
  Range `[85,95]` -> `[90,95]` -> 0.5. Range with no kernels -> 0.0. Zero-length range -> raises `ValueError` (not
  ZeroDivision) (pin).
- Adversarial: intervals fully nested (`[0,100]` containing `[10,20]`), touching (`[0,10]`,`[10,20]` -> 20 not 20+
  double count), duplicated rows, unsorted input, intervals crossing range edges (clipped), negative-duration or
  zero-duration rows ignored, 1e6 intervals runs in < 2 s.
- Share must satisfy 0 <= s <= 1 and sum-of-lengths-without-union >= union (property over 300 random sets; compare to a
  brute-force per-nanosecond boolean array for small coordinates, independent oracle).

**T114 gpu_active_csv_equals_json** (U)
- Write the same interval set as an nsys `cuda_gpu_trace` CSV (columns `Start (ns)`, `Duration (ns)`; mem ops included)
  and as a Chrome-trace JSON (`ph:"X"`, `ts`/`dur` in microseconds, `cat` kernel/gpu_memcpy/gpu_memset; CPU
  events and `cuda_runtime` events present as distractors that must be ignored). Expected: both shares equal (floats
  equal to 1e-12 — the one place a tolerance is acceptable, due to us->ns conversion; prefer using integer-us
  values so equality is exact).
- Failure: empty trace, wrong header, malformed JSON -> clear exception, non-zero exit via `-m` run.

---------------------------------------------------------------------------------------------------------------

## L. Cross-cutting, regression and failure-mode

**T120 regression_existing_suites_unchanged** (I)
- `tests/test_image_memo.py`, `test_train_fast.py`, `test_open_once.py`, `test_p2_train_rl.py` pass unchanged
  (with `per_image` default False). `make test` and `make lint` clean.

**T121 regression_todays_numbers** (E, box)
- Baseline (spec): 23.43 s/step R5. Fast arm step 1 vs baseline arm step 1: `rewards/iou_tc/mean`, `rewards/legal/mean`,
  `completions/mean_length` equal (exact float); `grad_norm`, `sampling/sampling_logp_difference/mean` equal to 4
  significant digits (`f"{a:.3e}" == f"{b:.3e}"` is NOT 4 s.f.: use `math.isclose(rel_tol=5e-4)` after rounding
  both to 4 s.f. with `float(f"{x:.4g}")`; flag boundary rollovers like 0.99995 vs 1.0000 with a `rel_tol=1e-3` assert
  as the fallback check). `reward` ignored.
- Same seed, two runs of the SAME arm: identical (determinism control; if not, G-EQ6 cannot discriminate and the
  test reports "baseline noise" instead of fail).

**T122 each_flag_off_equals_baseline** (I, tiny)
- With every new flag False (`--gpu-pixels off`, no memo, no encoder-once, no host-meta, no fast-sync, no group-requests,
  checks 0): a 2-step tiny GRPO run is bit-identical (losses, grads, logged metrics) to a run built from today's
  `FastGRPO` at HEAD (saved golden of losses from a seeded run: store as small JSON in `tests/goldens/`). Each flag
  turned on alone (7 runs) gives the same logps/loss (these are C1-C5, bit-exact) except `--group-requests` (statistical).
- Why: a flag whose "off" path changed behaviour; and flags that interact.

**T123 flag_pairs_interaction** (I)
- All 2^5 combos of {gpu-pixels, prompt-memo, encoder-once, host-meta, fast-sync} on the tiny trainer for 1 step: no
  exception, and loss `torch.equal` across combos. Combos with encoder-once but prompt-memo off, and host-meta off while
  encoder-once on (needs `host_meta`? spec: batch has `host_meta` only with host-meta path — pin behaviour; either
  works or raises a clear flag-dependency error at start-up).

**T124 importance_sampling_modes** (I)
- `vllm_importance_sampling_mode` in {`sequence_mask`, `token_truncate`} with each fast path on: step runs, `is_ratio`
  masks equal the baseline's; constraint from spec section 4.

**T125 non_qwen_untouched** (U)
- Each wrapper (rope, encoder-once, PromptMemo) with a non-Qwen model/processor stub installs nothing and the
  trainer's `_get_per_token_logps_and_entropies`/`_tokenize_prompts` are TRL's functions (`is`).

**T126 private_trl_name_drift** (U)
- Delete or rename `GRPOTrainer._tokenize_prompts` / `VLLMGeneration.sync_weights` on a stub: start-up check raises with
  the missing name in the message; does not fall through to a silent no-op.

**T127 failure_cuda_oom_and_partial_state** (G)
- Force OOM inside `patchify` (batch of 4096 images 512x512) -> exception propagates; `ImageMemo` state not
  corrupted (next normal call equals HF); `torch.cuda.memory_allocated` returns to the pre-call level after
  `gc.collect()` + `empty_cache()` (+ 1 MB).

**T128 failure_untrusted_input** (U)
- Image with a malicious huge declared size (PIL decompression bomb of 40000x40000 `L` via `Image.new` is memory heavy:
  use a lazily opened PNG header file in tmp) -> `UnsupportedImage` (or PIL's error) raised BEFORE allocating; no
  tensor over 100 MB created. NaN/inf can't occur in uint8; float tensor raises (T05).

**T129 determinism_repeat** (G)
- 10 repeated runs of `patchify`, `rope_positions`, `window_from_meta` on the same input: identical bytes (hash).
  Different CUDA stream and different `torch.backends.cudnn.benchmark` setting: identical.

---------------------------------------------------------------------------------------------------------------

## M. Performance gates (box, `cuda`, `slow`)

Protocol (spec 5): same host and seed, 17 steps, mean of steps 3-15, baseline 23.43 s.

**T140 perf_step_time** (E): fast arm mean <= 20.0 s awake / <= 20.7 s sleep; FAIL above 21.2 / 21.9 s; thresholds
recomputed from the measured one-call push time before the fast arm runs (assert the recompute step ran and logged).
Report the median and std; assert std of steps 3-15 < 0.3 s, else the A/B is noise (rerun, do not pass/fail).

**T141 perf_push** (G): `PushStats.seconds <= 0.2 s` (median of 20 pushes); awake: no `time/wake_*` or sleep phase keys;
sleep: `time/sync_s <= PushStats.seconds + baseline wake_weights_s + 0.1 s`.

**T142 perf_prepare_range** (G): NVTX `prepare` <= 0.5 s per step (from nsys steps 16-17); microbench: `patchify` of 256
rows of 512x512 on GPU < 20 ms warm; `rope_positions` for 256 rows < 5 ms; `row_meta` < 3 ms; `window_from_meta`,
`length_order`, `group_order` each < 1 ms. Each must be at least 5x faster than its device/CPU predecessor
(compare to `ip(...)` CPU: ~ seconds per 256 rows) or the optimisation is not kept.

**T143 perf_syncs** (G): trainer-side host syncs per step (outside `generate`) < 200 and <= 25 % of the baseline's measured
count; assert baseline count > 0 (so the ratio is meaningful) and that the count scales (not per row: the same
at B=128 and B=256 within +10 %).

**T144 perf_h2d_bytes** (G): trainer-side H2D bytes per step < 50 MB (from profiler memcpy HtoD events outside vLLM's
stream/thread); regression check that `patchify` takes uint8 H2D (256 x 512 x 512 x 3 = 201 MB if per row! so it
must send only 32 distinct images = 25 MB; test asserts H2D of the pass <= 32x512x512x3 + 1 MB = 25.2 MB).

**T145 perf_gpu_active** (G): `gpu_active.py` over steps 16-17 >= 70 %; sanity: script on a deliberately idle
trace (sleep 1 s between kernels) < 30 %.

**T146 perf_memory** (G): peak `max_memory_allocated` <= 30 GB on every arm (also with the 1,536-token answers, i.e. the
worst case window: force all 256 answers to max length via `max_tokens` logit bias); no growth across steps 3-15
(< 100 MB drift).

**T147 perf_group_requests** (E, G-ST1): reward/legal/IoU_tc within 2 standard deviations of the baseline's step-to-step
values (compute SD from baseline steps 3-15, require >= 10 samples); mean answer length within 5 %; >= 0.2 s per step
saved; both comparisons use the same seeds. The test fails closed if baseline SD is 0 (degenerate).

**T148 perf_scaling_property** (P, Mac, cheap): `row_meta` and `group_order` runtimes grow linearly: time at B=256 vs
B=2560 ratio < 15 (no quadratic list scans); `ImageMemo` on 256 rows with 32 images does 32 dictionary builds not 256
(count via counting wrapper, not wall time).

**T149 perf_no_regression_on_mac_cpu** (U): the Mac CPU path of every new function must not be > 3x slower than the
old path on the tiny shapes (guards accidental O(n^2) on CPU tests; `slow` marker not needed).

---------------------------------------------------------------------------------------------------------------

## Pytest sketches for the 5 most important tests

Importance ranking: (1) T43 rope false-good fuzz (silent wrong training signal if broken); (2) T11/T01 patchify
bit-equality; (3) T63 encoder-once alignment after shuffle; (4) T73 order/window equals device; (5) T113 gpu_active
exact union (the gate G-SP6 measuring tool itself; a wrong tool passes bad runs). T64/T25 are close follow-ups.

```python
# tests/test_c_rope.py
import numpy as np, pytest
torch = pytest.importorskip("torch")
from objdel.train import fast
from tests.test_train_fast import tiny_model, IMG, START, END   # helpers from the existing tests

def make_prompt(n_img_tokens=16, text=(5, 6, 7)):
    ids = [START] + [IMG] * n_img_tokens + [END] + list(text)
    types = [0] + [1] * n_img_tokens + [0] * (1 + len(text))
    return ids, types

def hf_positions(base, ids, mask, types, grid):
    pos, _ = base.get_rope_index(input_ids=ids, image_grid_thw=grid, attention_mask=mask,
                                 mm_token_type_ids=types)
    return pos

def test_rope_never_calls_a_mutated_row_good_unless_it_equals_hf():
    model = tiny_model()
    base = model.model
    p_ids, p_types = make_prompt()
    grid = torch.tensor([[1, 8, 8]])                    # 64 patches -> 16 merged tokens
    tmpl = fast.rope_template(base, torch.tensor([p_ids]), torch.tensor([p_types]), grid)
    rng = np.random.default_rng(0)
    flagged = 0
    for case in range(500):
        B, W = 4, 80
        ids = torch.zeros(B, W, dtype=torch.long); mask = torch.zeros(B, W, dtype=torch.long)
        types = torch.zeros(B, W, dtype=torch.long)
        for i in range(B):
            pad = int(rng.integers(0, 20)); ans = int(rng.integers(0, 30))
            seq, ty = list(p_ids) + [int(x) for x in rng.integers(0, 90, ans)], list(p_types) + [0] * ans
            ids[i, pad:pad + len(seq)] = torch.tensor(seq); types[i, pad:pad + len(seq)] = torch.tensor(ty)
            mask[i, pad:pad + len(seq)] = 1
            for _ in range(int(rng.integers(0, 3))):          # mutations
                j = int(rng.integers(0, W)); kind = int(rng.integers(0, 3))
                if kind == 0: ids[i, j] = int(rng.integers(0, 90))
                elif kind == 1: types[i, j] ^= 1
                else: mask[i, j] ^= 1
        pos, bad = fast.rope_positions(tmpl, ids, mask, types)
        flagged += int(bad.sum())
        for i in range(B):
            if not bad[i]:
                want = hf_positions(base, ids[i:i+1], mask[i:i+1], types[i:i+1], grid)
                assert torch.equal(pos[:, i:i+1], want), f"case {case} row {i}"
        assert pos.dtype == torch.int64 and bad.dtype == torch.bool
    assert flagged > 25                                      # the fuzz has teeth
```

```python
# tests/test_c_patchify.py
import numpy as np, PIL.Image, pytest
torch = pytest.importorskip("torch")
from transformers import Qwen2VLImageProcessor
from objdel.train import gpu_pixels

@pytest.fixture(scope="module")
def ip():
    return Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)

@pytest.mark.parametrize("backend", ["hf", "torch"])
def test_patchify_is_bit_equal_to_hf_for_random_valid_sizes(ip, backend):
    rng = np.random.default_rng(0)
    for case in range(200):
        n = int(rng.integers(1, 6))
        imgs = [PIL.Image.fromarray(rng.integers(0, 256, (int(h), int(w), 3), dtype=np.uint8))
                for h, w in zip(32 * rng.integers(1, 9, n), 32 * rng.integers(1, 9, n))]
        want = ip(imgs, return_tensors="pt")
        pv, grid = gpu_pixels.patchify(imgs, ip, device="cpu", backend=backend)
        assert torch.equal(pv, want["pixel_values"]), (case, backend)
        assert torch.equal(grid, want["image_grid_thw"]) and grid.device.type == "cpu" and grid.dtype == torch.int64

def test_patchify_rejects_what_the_processor_would_resize_or_cannot_read(ip):
    bad = [PIL.Image.new("RGB", (500, 500)), PIL.Image.new("RGBA", (64, 64)), PIL.Image.new("P", (64, 64))]
    for b in bad:
        with pytest.raises(gpu_pixels.UnsupportedImage):
            gpu_pixels.patchify([b], ip, device="cpu")
    with pytest.raises(ValueError) as e:
        gpu_pixels.patchify([PIL.Image.new("RGB", (64, 64))], ip, backend="nope")
    assert not isinstance(e.value, gpu_pixels.UnsupportedImage)   # must not trigger ImageMemo's fallback
```

```python
# tests/test_c_encoder_once.py
import numpy as np, pytest
torch = pytest.importorskip("torch")
from objdel.train import fast

def test_image_ids_follow_rows_through_shuffle_order_and_split():
    rng = np.random.default_rng(1)
    G, K = 8, 32
    for case in range(100):
        base_ids = np.repeat(np.arange(K), G)
        if case % 3 == 0:
            base_ids[G:2 * G] = base_ids[:G]                  # two groups share one image object
        plen = rng.integers(10, 40, K * G).astype(np.int32)
        clen = rng.integers(0, 200, K * G).astype(np.int32)
        meta = fast.RowMeta(plen, clen, base_ids.astype(np.int32), (np.arange(K * G) // G).astype(np.int32))
        for order in (rng.permutation(K * G), fast.length_order(meta), fast.group_order(meta, G)):
            ids = base_ids[order]
            for chunk in (8, 16, 4):
                for s in range(0, K * G, chunk):
                    rows = order[s:s + chunk]
                    assert np.array_equal(ids[s:s + chunk], base_ids[rows])
        # whole groups stay together under group_order
        go = fast.group_order(meta, G).reshape(K, G)
        assert (go // G == (go // G)[:, :1]).all()
```
(The sketch above tests the index bookkeeping at the unit level. The matching real test, T63 second half, runs the
tiny model with unequal grids and asserts `torch.equal` of per-row encoder rows against standalone tower outputs.)

```python
# tests/test_c_meta_equals_device.py
import numpy as np, pytest
torch = pytest.importorskip("torch")
pytest.importorskip("trl")
from objdel.train import fast

def build(rng, B=16, G=8, L=64, P=24):
    plen = rng.integers(4, P + 1, B); clen = rng.integers(0, L + 1, B)
    am = torch.zeros(B, P + L, dtype=torch.long)
    for i in range(B):
        am[i, P - plen[i]:P] = 1; am[i, P:P + clen[i]] = 1
    return am, plen, clen

@pytest.mark.parametrize("buckets", [fast.BUCKETS, fast.FINE, (32, 64, 128)])
def test_window_and_orders_equal_the_device_versions(buckets):
    rng = np.random.default_rng(2)
    for case in range(500):
        am, plen, clen = build(rng)
        L = 64
        meta = fast.row_meta([[1] * int(p) for p in plen], [[2] * int(c) for c in clen],
                             np.repeat([0, 1], 8), 8, mask_truncated=False, max_completion_length=L)
        for rows in (slice(0, 8), slice(8, 16), slice(0, 16), [3, 9, 12]):
            r = np.arange(16)[rows] if isinstance(rows, slice) else np.array(rows)
            assert fast.window_from_meta(meta, r, L, buckets) == fast.window(am[r], L, buckets), case
        n = am.sum(1)
        assert np.array_equal(fast.length_order(meta), torch.sort(n, descending=True, stable=True).indices.numpy())
        batch = {"prompt_mask": am[:, :24], "completion_mask": am[:, 24:], "x": torch.arange(16)}
        assert np.array_equal(fast.group_order(meta, 8), fast.by_group_length(batch, 8)["x"].numpy())
```

```python
# tests/test_gpu_active.py   (experiments/phase2/gpu_active.py)
import json, pytest
from experiments.phase2 import gpu_active

INTERVALS = [(0, 40), (30, 60), (50, 80), (90, 100), (100, 110)]      # ns; 2 streams, one memcpy
def brute(iv, lo, hi):
    on = [False] * (hi - lo)
    for s, e in iv:
        for t in range(max(s, lo), min(e, hi)): on[t - lo] = True
    return sum(on) / (hi - lo)

@pytest.mark.parametrize("rng_", [(0, 200), (10, 50), (85, 95), (150, 200), (0, 110)])
def test_union_share_is_exact_and_csv_equals_json(tmp_path, rng_):
    lo, hi = rng_
    csv = tmp_path / "t.csv"; csv.write_text("Start (ns),Duration (ns),Name\n" +
        "".join(f"{s},{e - s},k\n" for s, e in INTERVALS))
    js = tmp_path / "t.json"; js.write_text(json.dumps({"traceEvents": [
        {"ph": "X", "cat": "kernel", "ts": s / 1000, "dur": (e - s) / 1000, "name": "k", "pid": 0, "tid": i % 2}
        for i, (s, e) in enumerate(INTERVALS)] + [{"ph": "X", "cat": "cpu_op", "ts": 0, "dur": 500, "name": "cpu"}]}))
    want = brute(INTERVALS, lo, hi)
    assert gpu_active.share(csv, lo, hi) == pytest.approx(want, abs=1e-12)
    assert gpu_active.share(js, lo, hi) == pytest.approx(want, abs=1e-12)
```
Note: `gpu_active.share(path, lo, hi)` is an assumed name; the spec fixes only the module path and the behaviour. The
test author adapts the call to the real function.

## Open spec questions found while writing the tests
1. 2.1: a 3-D uint8 tensor `(H,W,3)` vs `(N,H,W)` is ambiguous; the spec lists numpy `(H,W,3)` only (T04).
2. 2.2: counters under the `UnsupportedImage` fallback are undefined (T25).
3. 2.8: `completion_len` is 0 for truncated rows, while `by_length` sorts by the real mask sum; equal only if TRL zeroes the
   mask for truncated rows (T73).
4. 2.11: two adjacent equal prompts are indistinguishable from one run of 2G (T100).
5. 2.12: the counter key name is not fixed (T110).
6. `UnsupportedImage` subclasses `ValueError` while unknown backend must raise plain `ValueError` (T09, T25).
7. G-EQ6 "4 significant digits" is rollover-sensitive (T121).
