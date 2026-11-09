# Phase C test cases, set 1 (written from C.spec.md only, no design notes)

## Skills used
- engineering:testing-strategy: sorted tests by pyramid level (unit/property/integration/GPU/e2e) and added perf gates.
- mattpocock-skills:tdd: every test goes through a public seam named in the spec; expected values come from an independent
  source (the HF processor, HF `get_rope_index`, TRL 1.17.1 methods, fixed literals), never recomputed with our code path.
- Not invoked (no effect on a test-only document): nvidia-skill-finder, code-review, debug, system-design, architecture,
  documentation, huggingface-best.

## Conventions
- Level `unit` = Mac CPU, no model download (tiny Qwen3-VL as in tests/test_train_fast.py, processor built as in
  tests/test_image_memo.py: `Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)`).
- Level `GPU` = marked `cuda`, runs on the box. `E2E` = box, real Qwen3-VL-2B.
- Oracle for every pixel test: `image_processor(list_of_images, return_tensors="pt")` of the same stack. For rope: `base.get_rope_index`.
  For prompts: TRL 1.17.1 `GRPOTrainer._tokenize_prompts`. For loss: TRL 1.17.1 `_compute_loss`.
- Equality means `torch.equal` / exact integer equality. No `allclose` anywhere except G-EQ6's 4-digit values and C6.
- Images are random uint8 from `np.random.default_rng(seed)`, never smooth gradients (smooth data hides rounding-order bugs).
  Add hostile pixel values: all 0, all 255, alternating 0/255, values 127/128 (rounding of x/255 - mean)/std).
- Test ids: C-<area><n>. Area codes: PX pixels, IM ImageMemo, ST start-up, RP rope, PI prompt ids, PM PromptMemo,
  EN encoder-once, RM row meta, LM loss metrics, SY sync, GR group requests, MS measurement, GA gpu_active, CF flags.

---------------------------------------------------------------------------------------------------------------------

## A. gpu_pixels.patchify

### C-PX1 patchify equals the processor, 512x512 RGB (unit, backend torch, CPU)
- Setup: processor as above. Images: 3 random 512x512x3 uint8 as PIL RGB.
- Input: `patchify(imgs, ip, device="cpu", backend="torch")`.
- Expected: `pv.shape == (3*1*32*32, 1536)`, `pv.dtype == float32`, `grid.tolist() == [[1,32,32]]*3`, grid dtype int64, grid.device cpu;
  `torch.equal(pv, ref["pixel_values"])`, `torch.equal(grid, ref["image_grid_thw"])`.
- Why: base case. The patch order (merge-block-major, then C, temporal, py, px) is easy to get transposed; 1536 columns
  all in the wrong order still have the right shape.

### C-PX2 hostile pixel values (unit)
- Input: four 512x512 images: all 0; all 255; checkerboard 0/255; every pixel in {127,128}. Plus one image containing every
  value 0..255 (tile `arange(256)` over the plane, channel-shifted).
- Expected: `torch.equal` to processor for each, individually and as one batch.
- Why: normalisation order. `(x*(1/255) - mean)/std` vs `(x/255 - mean)/std` vs fused `x*a+b` differ by 1 ulp on some of the 256
  values; only an exhaustive value sweep sees it. This is the most likely bit-equality killer.

### C-PX3 property: random sizes and contents (property, hypothesis-free loop, 200 cases, fixed seeds)
- Input: H, W each drawn from multiples of 32 in [32, 640]; pick N in 1..4; images random; modes RGB or L.
- Expected: equal to processor for every case that does not raise. If `smart_resize(H,W) != (H,W)` (min/max pixels reached),
  `UnsupportedImage` is raised, never a wrong result.
- Why: finds grid-shape bugs (h != w, h > w, h < w) that 512x512 hides, e.g. swapped gh/gw in the reshape.

### C-PX4 non-square: 96x64 and 64x96 (unit)
- Expected: `grid == [[1,6,4]]` and `[[1,4,6]]`; pixels equal the processor's. Batch of both in one call equals the processor on the list
  (rows concatenated in input order).
- Why: spec names 96x64 as valid. A swapped h/w gives right shape (384 x 1536 vs 24 x 1536 patches) but wrong order.

### C-PX5 minimum size, single patch block (unit)
- Input: 32x32 image (smallest multiple of patch*merge). If processor `min_pixels` rejects it, expected is `UnsupportedImage` and the
  test asserts that instead; whichever the processor does, patchify must agree (compute `smart_resize` in the test from the
  processor's own settings).
- Why: edge of the resize rule; guards off-by-one in the "would resize" test.

### C-PX6 500x500 raises (unit)
- Input: 500x500 RGB. Expected: `pytest.raises(UnsupportedImage)`; also `isinstance(e, ValueError)`. Also 512x500, 500x512, 513x512, 511x512.
- Why: spec literal; a wrong implementation would crop or pad silently.

### C-PX7 over max_pixels raises (unit)
- Input: an image whose area exceeds the processor's `max_pixels` but is a multiple of 32 (derive H=W from the processor's value).
- Expected: `UnsupportedImage`. Why: the resize rule has two sides; test the upper side too.

### C-PX8 grey equals converted RGB (unit)
- Input: random 512x512 `L` PIL image, random `(512,512)` uint8 tensor, `(512,512)` numpy.
- Expected: all three results `torch.equal` to `image_processor([pil.convert("RGB")])`, and to each other.
- Why: channel replication must happen before normalisation (the three channels have different mean/std; normalising a 1-channel
  tensor with channel-0 constants would be wrong for channels 1, 2).

### C-PX9 every accepted input type gives one answer (unit)
- Input: the same 4 images as (a) list of PIL RGB, (b) uint8 tensor (N,H,W,3), (c) list of numpy (H,W,3), (d) uint8 tensor on the
  same device.
- Expected: all four outputs `torch.equal`. Also a single `(H,W,3)` numpy array and a single `(H,W)` array: N == 1, grid `[[1,32,32]]`.
- Why: layout confusion between `(N,H,W,3)` and `(N,3,H,W)`; ambiguity when W == 3 (use a 32x... no: a (3,32,3)-shaped input must
  be rejected, see C-PX10).

### C-PX10 unsupported modes and channel counts (unit)
- Input, each must raise `UnsupportedImage`: PIL RGBA, PIL P, PIL I;16, PIL 1, PIL CMYK; numpy `(512,512,4)`, `(512,512,1)`,
  `(512,512,2)`; uint8 tensor `(N,512,512,4)`; a `(3,512,512)` channel-first numpy array (3 is not a size multiple of 32 ... must raise either
  `UnsupportedImage` or be read as H=3: assert it raises, not that it silently works).
- Also float32 and uint16 tensors: expected `UnsupportedImage` or `TypeError` (spec silent; test pins "never a wrong result" by comparing
  to the processor if it does not raise, and asserting a raise otherwise).
- Why: RGBA given to a "ignore alpha" implementation silently equals the processor only by luck; processor itself converts RGBA to RGB.
  The spec wants a raise, which makes `ImageMemo` fall back to the processor.

### C-PX11 unchecked processor settings raise (unit)
- Input: processors with `patch_size=14`, `merge_size=1`, `temporal_patch_size=1`, and `merge_size=4` (each with a size valid for it).
- Expected: `UnsupportedImage` for each. Also a processor whose `image_mean` is changed to (0.1,0.2,0.3) and `image_std` to (0.5,0.4,0.3): result
  must equal THAT processor's output (constants come from the object).
- Also `do_rescale=False`, `do_normalize=False`, `do_convert_rgb=False`, `rescale_factor=1/128`, `resample` changed: expected
  equal to processor if no raise, else `UnsupportedImage`. Never a wrong result.
- Why: spec says constants never from literals. A hard-coded mean passes C-PX1 and fails here. The settings test catches code that
  assumes `do_normalize=True`.

### C-PX12 out_dtype and device (unit + GPU)
- Input: `out_dtype=torch.bfloat16` and `torch.float16`, device cpu.
- Expected: `pv.dtype` matches; `torch.equal(pv, ref.to(out_dtype))` (cast of the fp32 reference, one rounding only). grid stays int64 CPU.
- GPU variant: `device="cuda"`: `pv.device.type == "cuda"`, `grid.device.type == "cpu"`; `torch.equal(pv.cpu(), ref)`.
  Default device: cuda when available else cpu.
- Why: casting before normalisation (double rounding) would differ from cast-after.

### C-PX13 backends (unit/GPU)
- `backend="auto"`, `"hf"`, `"torch"` on CPU; `"k3"` on GPU. Each equal to the processor if it returns; `"k3"` on CPU may raise but not silently
  differ. `backend="bogus"` and `backend=""` and `backend=None`: `ValueError` (and the error is raised before any work: patch a counter on the
  processor and assert it is not called).
- Why: "for every backend that returns without raising" is the spec's equality contract.

### C-PX14 empty input and zero images (unit)
- Input: `[]`, `torch.empty(0,512,512,3,dtype=uint8)`. Spec silent. Pin the behaviour: either `ValueError`/`UnsupportedImage`, or
  `pv.shape == (0,1536)` and `grid.shape == (0,3)`. Test passes with either, fails on a crash with an index error or a hang.
- Why: ImageMemo calls patchify with the missed-image list; `prefill([])` must not reach it with an empty list and break.

### C-PX15 inputs are not modified (unit)
- Input: a uint8 tensor and a numpy array. Copy before; call; assert unchanged bytes. Also the result does not alias the input
  (mutate the output, input stays).
- Why: an in-place normalise on a uint8->float view or a shared-memory `torch.from_numpy` would corrupt the caller's image that
  OpenOnce hands to the other rows.

### C-PX16 non-contiguous and read-only inputs (unit)
- Input: a transposed/sliced uint8 tensor (stride != contiguous), a read-only numpy array (`setflags(write=False)`), a PIL image from `Image.frombytes`.
- Expected: equal to the processor on the contiguous copy; no warning "non-writable".
- Why: `torch.from_numpy` of a read-only array warns; `.view` on non-contiguous raises.

### C-PX17 property: same image twice in a batch, order (property)
- Input: `[a, b, a, c]`. Expected: `pv` equals processor on `[a,b,a,c]`: rows for `a` appear twice, in order.
- Why: guards a dedup inside patchify that would reorder (dedup belongs to ImageMemo).

### C-PX18 CUDA vs CPU bit equality (GPU)
- Input: 64 random + 8 hostile images, backend `torch` on cuda and on cpu, and `hf`, `k3` where they exist.
- Expected: all pairwise `torch.equal`. Run under `torch.backends.cuda.matmul.allow_tf32=True` and `torch.set_float32_matmul_precision("high")`.
- Why: a GPU path using FMA contraction (`x*a+b` fused) changes the last bit versus CPU; TF32 flags must not matter (no matmul used).
  Also run with `CUDA_LAUNCH_BLOCKING=1` to catch race bugs.

### C-PX19 G-EQ1 scale test (GPU/E2E)
- Input: >= 2,000 views from the real training views code (centred, shifted, rescaled) and all 32 images of 5 recorded steps; box HF CPU processor as oracle.
- Expected: 100 % `torch.equal`. Report mismatch count and first mismatch index if any.
- Sync gate: `count_syncs()` around ONE pass-2 `image_processor` call (the ImageMemo path) gives 0 (see C-MS3).

### C-PX20 views with odd content (unit, uses repo data only with `python -I`; skip if no fixtures)
- Input: real LIDC crops and one RIGA grey crop from `data/fixtures` via the project's loader.
- Expected: equal to the processor. Why: real images are smooth with long constant runs; catches shortcuts that exploit nothing in random data.

### C-PX21 performance (GPU, gate)
- Input: 32 images 512x512. Expected: patchify (torch backend, cuda, images already on device) < 5 ms after warm-up (median of 50);
  `ImageMemo` pass-2 call < 10 ms. The CPU HF processor takes the baseline (record it in the test output; assert patchify is >= 5x faster).
- Why: the whole point of C1 is speed; a regression to a CPU round trip (`.cpu().numpy()`) would pass equality and fail speed.
  Mark it as a soft gate (skip on shared CI) but always run on the box.

### C-PX22 no host synchronisation inside patchify when input is on device (GPU)
- Input: uint8 CUDA tensor `(32,512,512,3)`. `count_syncs()` around the call. Expected: 0. The returned grid is CPU and built from shapes
  (no `.cpu()` of a device tensor). Why: spec G-EQ1 "0 syncs".

---------------------------------------------------------------------------------------------------------------------

## B. ImageMemo

### C-IM1 per_image False equals the wrapped processor (unit)
- Input: rows `[a,a,b,c,a,b,c,a]` (c is 96x64), call as in tests/test_image_memo.py, both flat and nested `[[x]]`, with and without
  `return_tensors="pt"`.
- Expected: same keys, `torch.equal` values. (Regression of today's test; keep it as is, run with GPU pixels on and off.)

### C-IM2 per_image True: patches once, first-appearance order (unit)
- Input: `memo.per_image = True`; rows `[b, a, b, c, a, a, c, b]` (distinct order b, a, c).
- Expected: `pixel_values` equals `ip([b, a, c])["pixel_values"]`, `shape[0] == sum of prod(grid of distinct)`; `image_grid_thw` has 8 rows equal to
  `ip(rows)["image_grid_thw"]`. A distinct-image sort by id (`id()` order) must not pass: use images allocated in a loop and shuffled.
- Why: ordering by dict of ids vs first appearance; the encoder-once path indexes by this order.

### C-IM3 per_image True with mixed sizes (unit)
- Input: rows `[c, a, c, a]` with c 96x64 and a 64x64. Expected: pixel_values = `ip([c, a])`, grid rows 4.
- Why: row counts differ per image; a stride computed from the first image breaks.

### C-IM4 equal-content but different objects are distinct (unit)
- Input: `a1 = image(0)`, `a2 = image(0)` (same pixels, two objects). Rows `[a1, a2, a1]`.
- Expected: per_image True: `pixel_values` holds 2 images (a1, a2). misses == 2, hits == 1. Spec says "distinct image object".
- Why: keying by content hash would be correct but is not the spec; keying by `id` after the object is garbage collected would alias.
  Variant: delete `a1` and allocate a new image; memo keeps the old object alive so `id` is not reused (assert `memo.memo` holds a reference).

### C-IM5 counters (unit)
- Input: 256 rows over 32 distinct images (OpenOnce-like: 8 identical object refs each), after `clear()`.
- Expected: `misses == 32`, `hits == 224`. After `clear()` + `prefill(32 images)`: counter after prefill is `misses == 32, hits == 0`; the
  same call adds `hits == 256`, `misses` unchanged (32). Also: `clear()` resets counters? Spec silent: pin it by asserting that
  after `clear()` a new call counts 32 new misses (cumulative or reset both fine; assert the DELTA).
- Why: exact spec numbers.

### C-IM6 prefill (unit)
- `prefill([a,b,a])`: processes 2 images, calls the wrapped processor once with 2 images (or patchify once); `prefill([a,b])` again: 0 calls, 0
  new misses. `prefill([])`: no call, no error. `prefill` then a call with different kwargs (`return_tensors` absent vs "pt"): still no
  second preparation. (Spec: "at most once between two clear() calls".)
- Counting: wrap `ip.__class__.__call__` as the existing test does; assert the sum of images seen == distinct count.

### C-IM7 at most once even if image objects interleave across calls (property)
- Input: random sequences of calls over a pool of 10 images with random kwargs from {none, return_tensors="pt"}; random `clear()`.
- Invariant: total images processed between clears == number of distinct objects seen; output always equals the wrapped processor's.

### C-IM8 UnsupportedImage falls back for the whole call (unit)
- Input: rows `[a(512), bad(500x500), a]` with the GPU path on.
- Expected: result `torch.equal` to wrapped processor over all rows (including `a`); no exception. After the call, the memo may hold
  entries but the result must still be right on a following call `[a, a]`. In `per_image` mode: the fallback result still equals the processor's
  output for the whole call (spec: "the result still equals the wrapped processor's output"). Pin: grid has 3 rows.
- Why: fallback that returns only `bad` via the CPU and `a` via GPU in a different order, or half-poisoned memo.

### C-IM9 RGBA / palette image falls back (unit)
- Input: RGBA PIL image among RGB. Expected equal to processor (processor converts RGBA to RGB itself). No raise.

### C-IM10 device and dtype of the result (GPU + unit)
- GPU: `pixel_values.device.type == "cuda"`, `image_grid_thw.device.type == "cpu"` in both per_image modes and in the fallback path.
  Mac: both on CPU. dtype float32.
- Why: spec; a grid on CUDA makes `grid.prod(-1).tolist()` a sync in TRL.

### C-IM11 clear releases memory (unit)
- After `clear()`: `memo.memo == {}`; image objects are collectible (weakref to a PIL image dies after `del` + `gc.collect()`).
  Counters follow C-IM5.

### C-IM12 deepcopy and getattr (unit)
- Keep the existing tests: `deepcopy(memo)` returns the plain processor (no cached pixels); `memo.merge_size == 2`. Add: `per_image`, `hits`,
  `misses` are attributes of the memo not forwarded; `copy.copy(bare)` before init raises AttributeError, not recursion.

### C-IM13 thread safety is not required, but re-entrancy is (unit)
- Reward prefetch runs on a thread while the generation uses the memo? Spec silent: skip; instead test that calling the memo from inside a
  `PromptMemo` call (C-PM) does not double count.

---------------------------------------------------------------------------------------------------------------------

## C. Start-up check (FastGRPO)

### C-ST1 first passing backend is chosen (unit, patchify monkeypatched)
- Setup: fake `patchify` where `hf` returns correct, `torch` correct. Expected: backend `"hf"` chosen (order hf then torch); the check used 8 images:
  3 training-view kinds present (centred, shifted, rescaled), 1 grey if a RIGA row exists, 1 generated digit image.
- Variant: `hf` returns a tensor differing in ONE element of ONE image (e.g. image 7, last column): `hf` fails, `torch` selected.
- Variant: no RIGA rows: 7 images or fewer, no grey image.
- Why: a check that compares only `shape` or only the first image would accept a broken backend.

### C-ST2 no backend passes (unit)
- Expected: `RuntimeError` naming the problem. With `allow_cpu_pixels=True`: no raise, CPU processor used, logged `gpu_pixels/backend=cpu`
  (value check on the metrics dict). Without CUDA (Mac): the check is skipped and no error (spec: "When CUDA is present").
- Also: `--gpu-pixels off`: no check, no patchify call at all (monkeypatched `patchify` raises if called).

### C-ST3 a backend that raises at start-up counts as failing (unit)
- `hf` raises `UnsupportedImage` or `RuntimeError`; `torch` is correct: `torch` selected, no crash. Both raise: RuntimeError (or CPU with the flag).

### C-ST4 the generated digit image is truly different (unit)
- The 8 images are not all the same object and not all equal; assert pairwise-distinct bytes. Why: a start-up check on 8 copies of one image is blind.

### C-ST5 the check images follow the dataset's sizes (unit)
- With a dataset whose views are all 512x512, the check images are 512x512 (valid). A dataset view that would be resized (an invalid size) does not make
  the start-up fail: it is skipped (UnsupportedImage is a fallback, not a start-up failure).

---------------------------------------------------------------------------------------------------------------------

## D. Rope positions

Shared setup: tiny Qwen3-VL from tests/test_train_fast.py (`tiny_model()`), `base = model.model`. Prompt template: `[text..., START, IMG*n, END, text...]` with
`n = t*h*w/4` for a grid. Use grid (1,4,4) for 4 pads, (1,8,6) for 12 pads, and the real 512x512 grid (1,32,32): 256 pads.
`mm_token_type_ids`: 1 on image pads, 0 elsewhere (check against the tiny model's expectation of the type ids).

### C-RP1 template + positions equal HF for plain left-padded batch (unit)
- Input: prompt Lp = 20 (with 4 pads), B = 6 rows, answers of lengths [0, 1, 5, 12, 12, 3] padded right; left padding on prompts of different real length is NOT needed
  when all prompts equal the template but add 3 left-pad columns to some rows (mask 0).
- Expected: `positions.shape == (3, B, W)`, int64, device same as input; for each row `torch.equal(positions[:, i], base.get_rope_index(...)[0][:, i])`
  (oracle called on the same full batch). `bad.tolist()` all False. Masked columns are 0.

### C-RP2 answer lengths 0, 1, 131, 1536 (unit)
- Input: rows with answer length 0, 1, 131, 1536 in one batch (W = Lp + 1536).
- Expected: bit-equal to HF for all; for 0-token rows, positions for real tokens stop at Lp-1.
- Why: spec names these lengths; positions after the last image continue as `max+1` increments, 1536 stresses int range/arange sizes.

### C-RP3 masked holes (unit)
- Input: attention mask with zeros inside the answer region (hole at column Lp+3 and Lp+4, real tokens after). Expected: equal to HF; HF's position
  for tokens after a hole equals cumulative count of real tokens (check literal: tokens after the hole continue without gaps) and masked columns 0.
- Why: spec says "masked holes" — a position = arange computed ignoring the mask would be off by the hole size.

### C-RP4 left padding of varying width (unit)
- Rows with 0, 1, 7, 100 left-pad columns; prompts identical otherwise. Expected: equal to HF; position of the first real token is 0 in every row.

### C-RP5 bad flags (unit)
- (a) first Lp real tokens differ in one id (text token changed) -> `bad[i]`; (b) same ids but a different `mm_token_type_ids` at one place -> bad;
  (c) a real token after Lp with type != 0 (an image pad in the answer) -> bad; (d) row with fewer than Lp real tokens (truncated prompt, e.g.
  only Lp-1 real) -> bad; (e) a row where the prompt differs only in masked (mask 0) columns -> NOT bad.
- Expected: `bad` exactly equal to the intended literal vector, bad rows' positions are unspecified (do not compare).
- Also: a row with a DIFFERENT image grid but same ids length (the pad count equal: (1,8,2) vs (1,4,4) both 4 pads). Spec: bad only on ids/types; the
  template has fixed positions from its own grid, so this row passes `bad == False` but positions would differ from HF. Pin the behaviour: the doc
  says "ids and types". Test asserts bad is False AND documents (xfail strict) that positions differ, so the owner sees the hole; flag it in the review.
  (This is the one case where spec-literal behaviour can yield a wrong rope: the check-rope-every flag is the guard, test C-CF3.)

### C-RP6 answer containing the end-of-text/special ids (unit)
- Answer tokens equal to IMG id (100) with type 0? Real vLLM answers cannot contain image pads but ids may equal START/END ids as text. Input: answer token = START (102)
  with type 0. Expected: HF, which keys on types, gives the same as the template (positions of text). `bad == False`. Why: an implementation keying on ids would misfire.

### C-RP7 single image vs multiple images (unit)
- Template with 2 images (grids (1,4,4) and (1,2,2)): equal to HF. 0 images in prompt (text-only template): equal to HF. 
- Why: the template walk has per-image offset logic.

### C-RP8 B = 1 and W = Lp (no answer columns) (unit)
- W == Lp (completion width 0): positions `(3,1,Lp)` and equal to HF. Why: slicing `[:, :, Lp:]` of width 0.

### C-RP9 no host sync (unit via monkeypatch + GPU)
- Unit: patch `torch.Tensor.item`, `.tolist`, `.cpu`, `.numpy`, `torch.nonzero`, `Tensor.__getitem__` with a bool mask (check by wrapping `__getitem__` and
  raising if the index dtype is bool) to raise; call `rope_positions`. Expected: no raise.
- GPU: `count_syncs()` around the call == 0, also with `torch.cuda.set_sync_debug_mode("error")` (raises on any sync).
- Why: the spec bans these. Mac patches catch the obvious; the GPU debug mode catches `torch.where(mask)` forms.

### C-RP10 property: random batches (property, 300 batches)
- Random B in 1..16, left pad 0..30, random answer lengths 0..64 with random holes, random bad rows (10 %).
- Invariant: for every non-bad row, equal to HF (oracle on the full batch, with bad rows replaced by template-valid ones so HF runs); `bad` equals a numpy reference
  computed independently from ids/types/mask.

### C-RP11 recorded batches (G-EQ2 gate, GPU/E2E)
- >= 50 batches (>= 15 from the box) x 256 rows; 0 bad rows on LIDC; every row `torch.equal`. Runs through `--record-batches` outputs: ids, masks, grids, positions.

### C-RP12 dtype/device invariants (unit + GPU)
- `positions.dtype == int64`, `bad.dtype == bool`, same device as `input_ids`. Template tensors live on CPU or device but results work with either (call with
  cuda inputs and a CPU template, and vice versa) without a host sync.

### C-RP13 template immutability (unit)
- `RopeTemplate` is frozen: assigning `tmpl.end = 1` raises `FrozenInstanceError`; `rope_positions` does not modify `tmpl.pos` in place (hash tensor before/after).

### C-RP14 interplay with RopeCache (unit, regression)
- Keep tests/test_train_fast.py's RopeCache tests green. If the new path is enabled, the model's forward gets `position_ids` equal to HF's
  (hook injects). `--check-rope-every 0` skips comparisons; 1 compares every step and raises `RuntimeError` on a deliberately corrupted position (flip one entry).

---------------------------------------------------------------------------------------------------------------------

## E. Prompt ids and PromptMemo

### C-PI1 vLLM path equals TRL (integration, tiny processor with a fake chat template)
- Input: 16 prompts (2 distinct images x 8 rows, OpenOnce-style same object) via `FastGRPO._tokenize_prompts(prompts)` with `use_vllm=True`.
- Expected: `prompt_ids == TRL._tokenize_prompts(prompts)[0]` (list of lists, exact), `images` is TRL's `[1]` with the same objects (`is`) in order, third element `{}`.
  Image processor NOT called (counting processor: 0 calls).
- Mixed: prompts with different text (two different system prompts) and image sizes (512 and 96x64): equal.

### C-PI2 once per distinct (text, size) (unit)
- Count TRL tokenizer-side invocations (patch the processor `__call__`). 256 rows: 32 images x 8 rows of identical text: <= 1 call per distinct (text, (H,W)) -> 1 call when
  all images are 512x512 and the text is the same. Two sizes: 2 calls. A second call to `_tokenize_prompts` (next step) with new images of the same size and the same text: 0 calls.
- The pad-count dependence: ids for (512,512) have 256 pads; for (96,64) 6 pads. A cache keyed only by text returns 256 pads for the small image: test asserts ids equal TRL's for both.

### C-PI3 non-vLLM path equals TRL (unit)
- `use_vllm=False`: result equals TRL's exactly, including the images object, and the processor IS called (as TRL does).

### C-PI4 text with the same rendered string but different image content (unit)
- Same text, two different 512x512 images: ids equal (same pad count); cached once. Different content must not matter. But same text, same size, image of mode L vs RGB: ids equal TRL's.

### C-PI5 cache lifetime (unit)
- Cache is per run: a new trainer instance starts empty; ids still right if the processor's tokenizer settings change between runs (new instance, new settings).

### C-PM1 install on non-Qwen processor (unit)
- A text-only processing class and a processor from another family (fake class without `image_processor`/qwen name): `install` returns None; `type(p)` and
  `p.__dict__` unchanged (compare a deep snapshot of attribute ids).

### C-PM2 install keeps isinstance and attribute access (unit)
- After `install`: `isinstance(proc, T)` True for every `T in type(proc).__mro__`; `proc.tokenizer is before_tokenizer`; `proc.chat_template` same; `proc.save_pretrained(tmp)`
  works and the saved files equal the pre-install files (hash). `copy.deepcopy(proc)` works. `pickle` if the original pickles.
- Why: TRL calls `isinstance(processing_class, ProcessorMixin)`, saves the processor at each checkpoint.

### C-PM3 output equality (unit)
- Input: `proc(images=[..], text=[..], padding=True, return_tensors="pt")` with 16 rows, mixed sizes (512 + 96x64) and mixed texts (3 different).
- Expected: same keys and `torch.equal` for each value: `input_ids`, `attention_mask`, `mm_token_type_ids` (if present), `pixel_values`, `image_grid_thw`. Pad side (left)
  identical, including when padding_side is overridden by a kwarg.

### C-PM4 original runs once per distinct (text, H, W) (unit)
- Count original `__call__`: 256 rows with 32 distinct images of one size, one text -> 1 call (per run, across repeated invocations within the run). Mixed 2 sizes x 2 texts -> 4.
  `pixel_values`/`image_grid_thw` come from `image_memo` (memo.hits increases; no second image-processor pass: counting image processor == 0 for the memo-covered images).

### C-PM5 other call signatures behave as before (unit)
- `proc(text=["a"], return_tensors="pt")` (no images); `proc(images=..., text=..., padding="max_length", max_length=300)`; `return_tensors=None`; `proc.apply_chat_template(...)`.
- Expected: identical to the pre-install processor (compare against a second, never-installed copy). No cache entries created (assert the cache length).

### C-PM6 text is a single string vs list; images nested vs flat (unit)
- `text="..."`, `images=img` (single, not list); nested `[[img]]`. Equal to original.
- Why: cache key construction on a non-list text fails with a TypeError; spec says equal output.

### C-PM7 input_ids lengths differ per row: padding stays right (unit)
- Rows whose texts have different lengths; with `padding=True` the batch width is the max; per-row cached ids (unpadded) must be re-padded the way the original does (left padding, pad id). Assert
  `attention_mask` pattern equals the original's. Why: caching padded rows from one call would paste the wrong width into another call.

### C-PM8 does not mutate returned cached tensors (unit)
- Mutate the result's `input_ids` in place after a call; call again with the same inputs; second result still equals the original. Why: aliasing cache storage.

### C-PM9 `image_memo` None / missing prefill (unit)
- With the memo cleared between calls: still correct. With `memo.per_image=True` set by the encoder-once path, PromptMemo must still return the full-row `pixel_values`
  expected by TRL's tokenize step OR the call is the vLLM path where pixel values are dropped: pin: results equal to the original when `per_image` False; when True, `image_grid_thw` has one row per row.

### C-PM10 TRL tokenizer parity, real processor (GPU box / Mac if cached)
- With the real Qwen3-VL-2B processor (skip if not cached): 8 real prompts built by the project's `prompt` code + the digit image; PromptMemo output `torch.equal` to the un-installed processor.

---------------------------------------------------------------------------------------------------------------------

## F. Vision encoder once per generation batch

### C-EN1 log-probs and entropies equal, packed and unpacked (unit, tiny model)
- Setup: tests/test_train_fast.py harness, 32-prompt-like batch scaled down to 4 prompts x 4 answers (B=16), `Plain` (pixel_values path) vs `Fast` with `mm_encoder_outputs`.
- Matrix: `shared_prefix` in {True, False} x attention in {sdpa, eager} x chunk size in {8, 4}.
- Expected: `torch.equal` for per-token log-probs and entropies in all 8 cells.
- Why: spec G-EQ4; the deepstack features must be passed with the pooled output, else the first decoder layers differ only in text rows touching deepstack indexes.

### C-EN2 LoRA gradients of one micro-batch equal (unit, tiny model with LoRA)
- Compute the loss on one micro-batch via both paths; `torch.equal` for every LoRA A/B grad (`p.grad`). Also the vision tower has `requires_grad False` and gets no grad on either path.

### C-EN3 vision tower called exactly once per generation batch (unit)
- Hook `model.visual.forward` (counter). Run reference pass (2 chunks) + 4 update micro-batches. Expected: counter == 1 per generation batch; after the next `_generate_and_score_completions`, counter == 2.
- Variant: the old-logps pass (when `skip_old_logps` False) also uses the cached output: still 1.
- Called with the DISTINCT images only: input patch rows == sum of distinct (e.g. 32 x 1024), not 256 x 1024; assert `visual` input shape.

### C-EN4 batch dict has no pixel_values, has host_meta and grid (unit)
- After `_generate_and_score_completions` (mocked generation): `"pixel_values" not in batch`, `"host_meta" in batch`, `batch["image_grid_thw"]` present with `B` rows.
  `len(batch["host_meta"]) == B` (row-aligned list).

### C-EN5 every forward uses mm_encoder_outputs and no pixel_values (unit)
- Spy on `model.forward` kwargs for reference, old, and update passes: `"mm_encoder_outputs" in kwargs and kwargs["mm_encoder_outputs"].keys() == {"image"}` and `"pixel_values" not in kwargs`.
- Also `--encoder-once` off: spy sees pixel_values and no mm_encoder_outputs (regression).

### C-EN6 shuffling and splitting keep image ids with rows (unit)
- Take a batch of 32 rows with image ids [0..3 x 8], shuffle with TRL's `shuffle_sequence_dict` (seeded), split in micro-batches of 8 or 4; for each micro-batch assert that the per-row image ids equal the ids of its rows (by tagging each row with its original index and image id).
- Expected: exact.
- Mutation check: permute the image-id field but not the rows -> the log-probs differ from the `pixel_values` path (test asserts the equality test would FAIL, to show the oracle has teeth: use two images with very different pixel data).

### C-EN7 missing image ids raises (unit)
- Delete the per-row image-id field from a micro-batch with this path on. Expected: forward raises (RuntimeError/KeyError/ValueError), never falls back to zero images or to pixel_values.
- Also: the id field present but out of range (id 99 with 32 images): raises, not wraparound/negative index; negative id -1 raises (guards Python negative indexing silently selecting the last image).

### C-EN8 rows without images (unit)
- A text-only row mixed into a batch (if the data ever has one): spec silent; test that an all-text batch with the path on does not call the encoder and equals TRL's.

### C-EN9 image order after dedup (unit)
- Batch rows reference images in order [2, 0, 1, 0, 2]; the `per_image` ImageMemo order is first-appearance [2, 0, 1]; each row's features equal the row's own image's `pixel_values` features. Use 3 images with different sizes (grids), since feature row counts differ (offsets).

### C-EN10 memory gate (GPU)
- During the generation batch with this path, `torch.cuda.max_memory_allocated` increase attributable to pixels < 100 MB; no tensor of shape `(256*1024, 1536)` exists (scan `gc.get_objects()` for CUDA tensors with that shape) -> none. G-MEM: peak <= 30 GB.

### C-EN11 bf16 real model equality (GPU/E2E, G-EQ4)
- Reference pass and 16 update chunks of one real batch. `torch.equal` log-probs, entropies; LoRA grads of one micro-batch equal; vision tower call count 1.
  Run twice (determinism): the same `pixel_values` path twice gives equal results first (otherwise the comparison is meaningless; if not, mark the test non-deterministic and report).

### C-EN12 the `ImageMemo.per_image` flag is set only with the flag (unit)
- `--encoder-once` off: `per_image` False. On: True and `FastGRPO` clears it before each generation. Another trainer instance does not share it.

---------------------------------------------------------------------------------------------------------------------

## G. Row metadata

Setup: python lists of token ids; compare each against the device derivation using tests/test_train_fast.py's batches and TRL-style tensors. No device tensor is read: run
with CUDA tensors "poisoned" (a tensor subclass raising on any data access) passed anywhere reachable; all functions take only numpy / lists.

### C-RM1 row_meta basics (unit)
- Input: G = 4, 8 rows, prompt_ids_list lens [20]*8, completion lens [0,1,5,3,7,2,4,6], image_ids [0,0,0,0,1,1,1,1], `mask_truncated=False`, max_completion_length 10.
- Expected: `prompt_len == [20]*8` int32, `completion_len == [0,1,5,3,7,2,4,6]` int32, `image_id == [0,0,0,0,1,1,1,1]` int32, `group == [0,0,0,0,1,1,1,1]` int32 (`i//G`). Shapes `(8,)`.

### C-RM2 truncated answers (unit)
- A row with completion len == max_completion_length and no end token: `mask_truncated=True` -> `completion_len == 0`; `mask_truncated=False` -> the count.
  A row with len == max and an end token as last token: kept (full length) under both. A row with len > max: spec silent: assert it raises or equals the count; pin whatever TRL does.
- Why: this mirrors TRL's `completion_mask * ~is_truncated`; the "ended with eos at the last column" case is the classic off-by-one.
- Edge: `max_completion_length = 1`, answers of len 1 with and without eos.
- Edge: eos id given how? The signature has no eos: the lists must carry the info (spec). Test uses lists where an end token id appears; if the function cannot know the eos, this exposes a spec gap: flag it (review note).

### C-RM3 image_id equality by object (unit)
- `image_ids` come from the caller; test the helper that builds them (or the caller): rows with the same image object share an id; distinct objects with equal content do not. dtype int32.

### C-RM4 window_from_meta == window (property, 300 random batches)
- Build random attention masks with left padding and right padded answers, W up to 2048, `logits_to_keep` in 1..1536, buckets in {BUCKETS, FINE}, random row subsets (the micro-batch rows).
- Invariant: `window_from_meta(meta, rows, ltk, buckets) == fast.window(attention_mask[rows], ltk, buckets)` exactly (tuple of ints).
- Edges: all answers length 0 (keep >= 1); one row; the longest answer equals logits_to_keep; the prompt has no left pad; the bucket exceeds W (min with W); rows with `completion_len == 0` due to truncation masking but the answer's columns are still present in the attention mask? (TRL masks truncated rows' completion_mask to 0 -> consistent: the device window sees the masked mask).
- Why: any difference changes chunk widths and thus shapes -> numbers.

### C-RM5 length_order (property)
- Equals the row order of `fast.by_length` for random batches (get the order by tagging rows with their index). Ties: stable (equal lengths keep index order). All equal lengths -> identity. Dtype int64 index array, a permutation (sorted == arange).

### C-RM6 group_order (property)
- Equals `by_group_length(batch, G)`'s order; whole groups stay contiguous (`order.reshape(-1,G) // G` constant per row, and within a group the original member order `j` ascending); ties between groups stable; the group key is the max member length. G in {1, 2, 8}. G = 1 equals `length_order`.
- Edge: B not multiple of G -> `ValueError` (the device version fails at `view`); pin that it raises.

### C-RM7 packed-chunk indices (property)
- `seg`, `idx` from RowMeta equal the ones `packed_chunk` builds from device tensors (call the helper on both inputs, or extract by instrumenting `packed_chunk`): `torch.equal` / `np.array_equal` on 300 random chunks with groups of 8 and 4, answer lengths 0..40, including 0-length answers.

### C-RM8 no device reads (unit)
- Pass CUDA-less "tripwire" tensors where possible; monkeypatch `torch.Tensor.item/tolist/cpu/numpy/__array__` to raise; call the four functions. Expected no raise.

### C-RM9 host_meta list alignment (unit)
- `batch["host_meta"]` is a Python list of length B; after `shuffle_sequence_dict` + split, each micro-batch's `host_meta` rows correspond to its rows (tag by index). A shared RowMeta object in each row would be shared state: the test asserts per-row independence (mutating one row's entry does not change another's).

### C-RM10 --check-meta-every (unit)
- With flag 1 and a corrupted meta (`completion_len[3] += 1`): trainer raises on that step; with flag 0: no comparison, no raise; default 50: compares on step 1 and step 51 but not step 2 (count comparisons with a spy over 52 steps = 2).

### C-RM11 G-EQ3 recorded (GPU/E2E)
- >= 50 batches (>= 15 box): windows, orders, runs, packed indices equal.

---------------------------------------------------------------------------------------------------------------------

## H. Loss metrics (only if adopted)

### C-LM1 zero item() calls (unit, tiny model)
- Patch `torch.Tensor.item` with a counting wrapper (also `.tolist`, `float(tensor)`); run `_compute_loss` for `loss_type == "dr_grpo"`; counter == 0.

### C-LM2 loss equals TRL's (unit)
- Same inputs as TRL 1.17.1 `GRPOTrainer._compute_loss`; `torch.equal`. Cases: beta 0 and 1e-3, `sequence_mask` and `token_truncate` importance sampling, all-zero mask rows
  (completion_len 0), all rows truncated, advantages zero, one-row batch, extreme ratios (logp diff +-20) to test clipping.
- Also a non-`dr_grpo` loss type: falls back to TRL's behaviour (item calls allowed).

### C-LM3 logged metrics equal TRL's (unit)
- After `FastGRPO.log`, `entropy`, `kl` (beta != 0), `clip_ratio/low_mean`, `clip_ratio/high_mean`, `clip_ratio/region_mean`, `clip_ratio/low_min`, `clip_ratio/high_max` equal TRL's
  values (exact, same reduction order). Cases: no clipped tokens (min over empty set: TRL's value for empty), all tokens clipped, mask all zero (division by zero -> same as TRL, nan vs 0).
- Also metrics accumulate across micro-batches and are logged once with the mean over micro-batches.

### C-LM4 pinned hash (unit)
- Replace `inspect.getsource(GRPOTrainer._compute_loss)` via monkeypatch with a modified source: start-up raises `RuntimeError` whose message contains `_compute_loss`. With the real source (TRL 1.17.1): no raise.
  Whitespace/comment-only change also raises (hash is of exact source): assert it, as spec says "does not hash to the pinned SHA-256".
- Why: silent drift when TRL is upgraded.

### C-LM5 flag off keeps TRL's method (unit)
- With the copy off: `FastGRPO._compute_loss` is the previous one (the current behaviour, tests/test_train_fast.py unchanged).

---------------------------------------------------------------------------------------------------------------------

## I. vllm_sync.push_weights / install

Use a fake `vg` with fake `llm` recording calls: `wake_up(tags=...)`, `load_weights` receives an iterable (record by consuming), `reset_prefix_cache()`.
Tiny PEFT model (LoRA r2 on a tiny Qwen3-VL language model), `lora_merged_params(model)`.

### C-SY1 load_weights once with every tensor (unit)
- Expected: `load_weights` called exactly 1 time; consumed `(name, tensor)` list has the same names, in any order, as `list(lora_merged_params(model))`, and `torch.equal` values (compare dict by name);
  `n_tensors == len(list)`, `bytes == sum(t.numel()*t.element_size())`. Note the iterable may be a generator: it is consumable once; a test that consumes it twice must fail the fake.
- Merged check: with non-zero LoRA B, the pushed weight for a target layer equals `W + delta` (compute with `m.get_delta_weight`) and differs from W; for a non-target layer equals W.

### C-SY2 no parameter changes (unit)
- Snapshot `{n: p.clone()}` of model params (including LoRA weights, bf16 base) before; after `push_weights` all `torch.equal`; `layer.merged` flags unchanged (False); repeating 100 pushes: still equal (no drift; TRL's merge/unmerge would drift in bf16 — build the bf16 case where `W + d - d != W` to prove the oracle: assert that merge/unmerge changes W while push_weights doesn't).

### C-SY3 sleep mode ordering (unit)
- Engine with `sleep mode` enabled: call order == `["wake_up(tags=['weights'])", "load_weights", "reset_prefix_cache"]`; without sleep mode: no `wake_up` at all. `reset_prefix_cache` called exactly once per push (not zero, not per tensor).
- `wake_up` raising: the exception propagates and `load_weights` is not called.

### C-SY4 PushStats.seconds measures only load_weights (unit)
- Fake `wake_up` sleeps 0.3 s, fake `load_weights` sleeps 0.05 s, reset 0.1 s. Expected: `0.05 <= stats.seconds < 0.15`. Also if the iterable is lazy, generating tensors happens inside `load_weights` consumption: spec "describes the iterable" — `n_tensors` and `bytes` are still right (counted while iterating).

### C-SY5 empty / odd models (unit)
- Model with no LoRA layers (deltas empty): pushes base weights unchanged, n_tensors == number of params. Model with `modules_to_save`, tied embeddings (name dedupe: spec says whatever `lora_merged_params` yields; assert equality of the name list with the generator, not our expectation).
- DoRA / LoRA bias: `NotImplementedError` from the generator surfaces from `push_weights` (no swallowing).

### C-SY6 install conditions (unit)
- Returns True and replaces `trainer.vllm_generation.sync_weights` (instance attribute; the class method is untouched) when: vLLM generation present, PEFT model, single process, no FSDP.
  Returns False and `sync_weights` is the same object as before when: no `vllm_generation` (None), non-PEFT model, FSDP (`accelerator.state.fsdp_plugin` set), `num_processes > 1`.
  Calling the installed `sync_weights()` performs a full push (spy `push_weights` called once); two consecutive syncs = two full pushes (no "only changed" logic).
- Idempotence: calling `install` twice returns True and does not wrap twice (one push per sync).

### C-SY7 G-EQ5 (GPU/E2E)
- After one `push_weights`, every vLLM parameter `torch.equal` to the state after TRL's per-parameter push of the same adapter state. Awake: 3 pushes on 3 adapter states (perturb LoRA B randomly). Sleep mode: after a sleep/wake cycle.
- vLLM stacked params (`qkv_proj`, `gate_up_proj`) are loaded via the `load_weights` mapping: compare after load. Check generation determinism: with temperature 0, 16 prompts, generated token ids equal in both arms.

### C-SY8 G-SP2 performance (GPU gate)
- `PushStats.seconds <= 0.2` for the real model; in sleep mode `time/sync_s <= PushStats.seconds + baseline wake_weights_s + 0.1`.

---------------------------------------------------------------------------------------------------------------------

## J. group_requests

### C-GR1 happy path (unit)
- 4 distinct prompt dicts each repeated 8 times (same object refs for images) -> returns `prompts[::8]` (4 dicts, `is` the same objects as inputs 0, 8, 16, 24).

### C-GR2 errors (unit)
- A run of 7 then 9; a run of 8 followed by a run broken in the middle (`[a]*8 + [b]*4 + [c]*4`); length not multiple of G (25 rows, G=8) -> `ValueError` each.
- `a, a, ..., a` (8) then the same prompt again (`a` x 16): run of 16, G = 8 -> `ValueError` (spec: exactly G long) — or two groups of the same prompt? Spec: "consecutive runs of exactly G identical". Pin: 16 identical rows with G=8: `ValueError` (it's one run of 16). Flag in review as ambiguity: the dataset never does this when each prompt is unique, but two consecutive samples with the same image + text would be legal in a dataset with repeats.
- Equality semantics: "same image object": two rows with equal `prompt_token_ids` but different image objects of equal content -> `ValueError` (not identical).
- G = 1: returns all prompts; empty list, G = 8: returns `[]` (0 is a multiple); G = 0 -> `ValueError` (no ZeroDivisionError).

### C-GR3 completions re-expanded in original order (integration, fake vLLM)
- Fake `llm.generate(prompts, sampling_params)` where `sampling_params.n == G` and returns for request k outputs with token ids `[k*100 + j]`.
- Expected: trainer gets `len(prompts)` completions; row `k*G + j` has the ids `[k*100 + j]`; `vllm` received `len(prompts)/G` requests; `n == G`; the `logprobs`, `finish_reason`, and `prompt_token_ids` alignments hold for every row. Default off: vLLM receives all rows with n == 1 (regression).
- With `one_pad_per_image` patch (`patch_vllm_prompts`) active: the length check of prompt ids per request still passes.

### C-GR4 seeds (integration)
- With `n == G`, the G answers of a request are different (not copies): fake engine can't test; the GPU test (C-GR5) asserts the G answers for a prompt are not all identical at temperature 1.0 and top-p 1 (at least 2 distinct among 8 on average across 32 prompts).

### C-GR5 statistical gate G-ST1 (GPU/E2E)
- Mean reward, legal and IoU_tc within 2 stds of the baseline's step-to-step values; mean length within 5 %; >= 0.2 s per step saved. Report the effect sizes, not just pass/fail.

---------------------------------------------------------------------------------------------------------------------

## K. Measurement (count_syncs, timed, gpu_active)

### C-MS1 count_syncs without CUDA (unit, Mac)
- Inside the block `counter` value 0, no exception; nested use works; exception inside the block propagates and nothing is left patched.

### C-MS2 count_syncs counting rules (GPU)
- 1000 `.item()` calls in a loop on a CUDA scalar -> counter value == 1000 (not 1: repeats from the same line count). `torch.cuda.synchronize()` x 5 -> 5 (or note
  that explicit sync is counted by the debug mode? Spec: "every synchronizing operation reported by `set_sync_debug_mode('warn')`" — assert items count >= 1000 for `.item()`
  and exact 1000 for pure item loops). A block with no sync -> 0. Warnings are not printed to stderr in a flood (use `warnings.catch_warnings` with a handler; spec: counts).
- Mode restore: set the mode "error" before, enter/exit, the mode after == "error". Also restore on exception inside the block. The initial mode "default" restored too.
- Warning filter state restored (the global `warnings.filters` list equal before/after).

### C-MS3 0 syncs around one pass-2 processor call (GPU, G-EQ1)
- `count_syncs()` around `ImageMemo.__call__` for 256 rows after `prefill`: 0.

### C-MS4 timed logs perf/host_syncs only when --count-syncs on (unit/GPU)
- Flag off: `_metrics` has `time/<phase>_s` and no `perf/host_syncs`. On: both, per phase. `generate` phase: its syncs go in `perf/host_syncs_vllm` and NOT in
  `perf/host_syncs` for that phase (assert the key absent or 0). Nested phases (a phase inside a phase): spec silent: assert the inner counts are not double counted in the outer (or document).
- `time/<phase>_s` still logged as today (regression: tests/test_train_fast.py).
- On Mac (no CUDA): flag on logs 0 and does not raise.

### C-MS5 wake timers (unit with fake llm)
- `time/wake_weights_s` and `time/wake_kv_s` logged separately: `wake_up(tags=["weights"])` charges only the first, `wake_up(tags=["kv_cache"])` only the second; `wake_up()` with no tags: both or the legacy `wake_s`? Pin: no tags -> charged to a documented key and never lost (sum of all wake keys == total).

### C-GA1 gpu_active exact share (unit, Mac)
- Synthetic trace (module-run script, tested through its function): range [0, 100] us. Stream A kernels: [0,30], [20,50]; stream B memcpy [40,60], memset [90,95]. Union = [0,60] + [90,95] = 65 -> share 0.65 exactly.
- Variants: no intervals -> 0.0; one interval covering the range -> 1.0; intervals partially outside the range are clipped (an interval [-10, 10] contributes 10); an interval fully outside contributes 0; identical duplicated intervals counted once; touching intervals [0,10],[10,20] union 20; zero-length interval 0; range of length 0 -> `ValueError` (no ZeroDivisionError).
- Non-GPU events (CPU ops, `cudaLaunchKernel` runtime calls, NVTX ranges) in the trace are ignored.

### C-GA2 two formats agree (unit)
- The same intervals written as an nsys `cuda_gpu_trace` CSV (columns Start (ns), Duration (ns), ...; units ns) and as a torch.profiler Chrome trace JSON (ts, dur in us; categories `kernel`, `gpu_memcpy`, `gpu_memset`): `share_csv == share_json` (exact for integer-us intervals: convert units exactly, use integer arithmetic or `fractions`).
- Unit-mix error: CSV in ns vs JSON in us; test with intervals not multiples of 1000 ns? JSON ts is a float in us: use values with exact binary fractions (e.g. .5 us).
- Malformed file (missing columns, empty file): raises with a clear error, not a 0 share.
- Run as a module (`python -I -m experiments.phase2.gpu_active trace.csv --range 0 100`): exit code 0 and printed share parsed equals the function's.

---------------------------------------------------------------------------------------------------------------------

## L. Flags and wiring

### C-CF1 train_rl argument parsing (unit)
- Defaults: `--gpu-pixels auto`, `--allow-cpu-pixels` False, `--group-requests` False, `--check-rope-every 50`, `--check-meta-every 50`, `--no-ckpt-layers 0`. `--gpu-pixels bogus` -> SystemExit(2).
  `--profile-steps 3-5` parsed as (3, 5); `3`, `5-3`, `a-b`, `-1-2` -> error. `--no-ckpt-layers -1` error; `K` larger than the layer count (28) -> error (no silent clamp).
- `g_probe.py` accepts the same flags (spec section 3 title says both).

### C-CF2 flags independent (unit)
- Each of `--prompt-memo --encoder-once --host-meta --fast-sync --group-requests` alone builds a trainer config (tiny) that runs one step; all together: runs one step. Equality step 1 of every
  subset vs the baseline on the tiny model: `rewards`/`completions/mean_length` equal. 2^5 = 32 subsets is cheap on a tiny model: run all.

### C-CF3 --check-rope-every raises on a difference (unit)
- Corrupt a position; step 1 raises (default 50 -> step 1 and step 51 compared). 0 -> no raise ever. Flag N=3: compared on steps 1, 4, 7 (spy count over 7 steps = 3).

### C-CF4 --record-batches (unit)
- Directory created if missing; each generation batch saves ids, masks, grids, metadata, positions in files readable back (`torch.load(weights_only=True)`) and equal to the batch; the
  file count == steps; an existing non-empty directory: spec silent -> does not overwrite silently? Pin: raise or append with distinct names, never overwrite.

### C-CF5 --no-ckpt-layers K (unit/GPU)
- K layers have checkpointing off: count modules with `_gradient_checkpointing` toggled == K (the first K? spec silent; assert exactly K are off); K=0 unchanged (regression).
  Log-probs and grads `torch.equal` with K = 0 and K = 4 (checkpointing is numerically exact only if dropout is off and the RNG is not used; the test asserts bit-equality because
  LoRA dropout must be 0: if the config has `lora_dropout > 0` the test sets it 0 and notes it).
- Memory: peak memory with K = 4 increases but <= 30 GB (G-MEM).

### C-CF6 importance sampling modes (unit/GPU)
- Run steps of tiny trainer with `vllm_importance_sampling_mode` = `sequence_mask` and `token_truncate` with every new flag on: no exception; shapes right; for `sequence_mask`, the mask row behaviour equals the baseline.

### C-CF7 every private TRL name is in the A7 start-up check (unit)
- Monkeypatch each of `_tokenize_prompts`, `_compute_loss`, `VLLMGeneration.sync_weights` to a function with a different source: start-up raises naming the method. With TRL 1.17.1 as installed: passes.
  Each new wrapper installs only for Qwen-VL: for a non-Qwen processor/model (a tiny Llama), `FastGRPO` leaves TRL's `_get_per_token_logps` path (spy) and does not hook rope/encoder.

---------------------------------------------------------------------------------------------------------------------

## M. End-to-end and performance gates (GPU box)

### C-E1 G-EQ6 step 1 fast vs baseline (E2E)
- Same host, seed, data order. `rewards/iou_tc/mean`, `rewards/legal/mean`, `completions/mean_length` exactly equal; `grad_norm` and `sampling/sampling_logp_difference/mean` equal to 4 significant digits
  (compare `f"{x:.3e}"` strings); `reward` not compared. If the sampling differs at step 1, bisect by flag (run each of C1-C5 alone) and report which flag breaks equality.
- Nondeterminism control: run the BASELINE twice first; if the baseline differs from itself, the gate is meaningless -> report as "vLLM not deterministic" rather than fail the fast arm.

### C-E2 regression of today's numbers (E2E/unit)
- From REPO docs/research and memory (not read here): the R5 numbers (31.0 -> 23.6 s/step, baseline 23.43 s of the spec) are the baseline; the test asserts the fast arm's mean of steps 3-15 against spec thresholds:
  <= 20.0 s awake, <= 20.7 s sleep; hard fail > 21.2 / > 21.9; thresholds recomputed from the measured push time. The test prints the thresholds it used.
- Unit side regression: all existing tests (`make test`) pass unchanged: tests/test_train_fast.py, tests/test_image_memo.py, tests/test_open_once.py, tests/test_p2_train_rl.py.

### C-E3 G-SP3..SP6 and G-MEM from the profile (E2E)
- `prepare` NVTX range <= 0.5 s per step (mean over steps 3-15); trainer-side syncs per step < 200 and <= 25 % of the baseline's; host-to-device bytes < 50 MB per step; GPU-active share >= 70 % over the step
  (steps 16-17 nsys); peak memory <= 30 GB on every arm. Each gate has its own test id so a single failure does not mask the others.

### C-E4 repeatability (E2E)
- Run the fast arm twice with the same seed: step times within 5 %, step-1 numbers identical. Report the spread, so a noisy gate (20.0 vs 20.7) is not decided by one run.

### C-E5 failure modes
- Kill vLLM mid-push (fake): the error propagates; the trainer does not proceed with half-updated weights. A `load_weights` that raises: `reset_prefix_cache` not called; the exception message includes the parameter name if known.
- CUDA OOM during `patchify` of a 1000-image batch: not required; skip.
- `UnsupportedImage` for one image in a step with GPU pixels: that call falls back to the processor; training continues; the log has `gpu_pixels/fallbacks` (if present) counting 1. Equality of that step to the baseline.

---------------------------------------------------------------------------------------------------------------------

## Pytest sketches for the 5 most important tests

Ranking rule: (1) breaks numbers silently, (2) breaks the equality contract of C1-C5, (3) hard to see in a normal run.

### 1. C-PX2 hostile pixel values (patchify equals the processor on all 256 values and both extremes)
```python
import numpy as np, PIL.Image, pytest
torch = pytest.importorskip("torch")
from transformers import Qwen2VLImageProcessor
from objdel.train import gpu_pixels

def ip():
    return Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)

def hostile():
    ramp = (np.arange(512 * 512 * 3) % 256).astype(np.uint8).reshape(512, 512, 3)
    check = np.indices((512, 512)).sum(0) % 2 * 255
    return [np.zeros((512, 512, 3), np.uint8), np.full((512, 512, 3), 255, np.uint8),
            np.repeat(check[..., None], 3, 2).astype(np.uint8),
            np.where(np.indices((512, 512)).sum(0)[..., None] % 2, 127, 128).astype(np.uint8).repeat(3, 2),
            ramp, np.roll(ramp, 1, axis=2)]

@pytest.mark.parametrize("backend", ["hf", "torch"])
def test_patchify_is_bit_equal_on_hostile_pixels(backend):
    proc, imgs = ip(), [PIL.Image.fromarray(a) for a in hostile()]
    want = proc(imgs, return_tensors="pt")
    pv, grid = gpu_pixels.patchify(imgs, proc, device="cpu", backend=backend)
    assert pv.dtype == torch.float32 and pv.shape == (len(imgs) * 1024, 1536)
    assert grid.dtype == torch.int64 and grid.device.type == "cpu"
    assert torch.equal(grid, want["image_grid_thw"])
    assert torch.equal(pv, want["pixel_values"])
```

### 2. C-PX11 normalisation constants come from the processor
```python
def test_constants_are_read_from_the_processor():
    proc = ip()
    proc.image_mean, proc.image_std = [0.1, 0.2, 0.3], [0.5, 0.4, 0.3]
    img = PIL.Image.fromarray(np.random.default_rng(0).integers(0, 256, (512, 512, 3), dtype=np.uint8))
    pv, _ = gpu_pixels.patchify([img], proc, device="cpu", backend="torch")
    assert torch.equal(pv, proc([img], return_tensors="pt")["pixel_values"])
    for bad in (dict(patch_size=14), dict(merge_size=1), dict(temporal_patch_size=1)):
        with pytest.raises(gpu_pixels.UnsupportedImage):
            gpu_pixels.patchify([img], Qwen2VLImageProcessor(**{**dict(patch_size=16, temporal_patch_size=2, merge_size=2), **bad}),
                                device="cpu", backend="torch")
```

### 3. C-RP10 rope_positions equals HF get_rope_index on random batches with holes, padding, 0/1/131-token answers
```python
from objdel.train import fast
# tiny_model(), IMG, START, END from tests/test_train_fast.py; build() makes (ids, mask, types, grid) with a fixed 4-pad prompt.
@pytest.mark.parametrize("seed", range(20))
def test_rope_positions_equal_hf_on_random_batches(seed):
    rng = np.random.default_rng(seed)
    base = tiny_model().model
    ids, mask, types, grid = build_batch(rng, Lp=20, B=8, answer_lens=rng.choice([0, 1, 5, 131], 8),
                                         left_pad=rng.integers(0, 9, 8), holes=True)
    tmpl = fast.rope_template(base, ids[:1, first_real(mask[0]):Lp_end], types[:1, ...], grid[:1])
    pos, bad = fast.rope_positions(tmpl, ids, mask, types)
    want, _ = base.get_rope_index(ids, mm_token_type_ids=types, image_grid_thw=grid, attention_mask=mask)
    assert pos.dtype == torch.int64 and pos.shape == want.shape and not bad.any()
    assert torch.equal(pos, want)
    assert (pos.permute(1, 0, 2)[mask == 0] == 0).all()  # masked columns are 0

def test_rope_bad_flags_are_exact():
    ...  # C-RP5: five corrupted rows -> bad == [F, T, T, T, T, F(masked-only diff)]
```

### 4. C-EN1 log-probs and entropies equal, pixel_values path vs mm_encoder_outputs path
```python
@pytest.mark.parametrize("shared_prefix", [False, True])
@pytest.mark.parametrize("attn", ["sdpa", "eager"])
@pytest.mark.parametrize("chunk", [8, 4])
def test_encoder_once_path_is_bit_equal(shared_prefix, attn, chunk):
    model = tiny_model(); model.config._attn_implementation = attn
    batch = make_batch(model, images=4, G=4)                 # distinct images of different pixel data
    ref = Plain_logps(model, batch, chunk)                   # TRL path with pixel_values
    t = trainer(Fast, shared_prefix=shared_prefix, chunk=chunk)
    out = Fast_logps_encoder_once(t, model, batch, chunk)    # mm_encoder_outputs + per-row image ids
    assert torch.equal(out.logps, ref.logps)
    assert torch.equal(out.entropies, ref.entropies)
    assert visual_calls(model) == 1                          # frozen tower ran once over the distinct images

def test_missing_or_negative_image_ids_raise():
    ...  # C-EN7: ids removed, id 99, id -1 -> raises, never falls back
```

### 5. C-IM2 / C-IM5 / C-IM8 ImageMemo per_image order, counters, fallback
```python
def test_per_image_order_counters_and_fallback(counted):
    ip, seen = counted
    a, b, c = image(0), image(1), image(2, 96, 64)
    memo = fast.ImageMemo(ip); memo.per_image = True
    rows = [b, a, b, c, a, a, c, b]
    got = memo(rows, return_tensors="pt")
    plain = Qwen2VLImageProcessor(patch_size=16, temporal_patch_size=2, merge_size=2)
    assert torch.equal(got["pixel_values"], plain([b, a, c], return_tensors="pt")["pixel_values"])
    assert torch.equal(got["image_grid_thw"], plain(rows, return_tensors="pt")["image_grid_thw"])
    assert (memo.misses, memo.hits) == (3, 5)
    memo.clear(); big = [image(i) for i in range(32)]
    d0 = (memo.misses, memo.hits)
    memo.per_image = False; memo(sum(([x] * 8 for x in big), []), return_tensors="pt")
    assert (memo.misses - d0[0], memo.hits - d0[1]) == (32, 224)
    bad = PIL.Image.fromarray(np.zeros((500, 500, 3), np.uint8))          # UnsupportedImage -> wrapped processor
    mixed = [a, bad, a]
    out = fast.ImageMemo(ip)(mixed, return_tensors="pt")
    assert torch.equal(out["pixel_values"], plain(mixed, return_tensors="pt")["pixel_values"])
```

---------------------------------------------------------------------------------------------------------------------

## Spec ambiguities and holes the tests exposed (for the owner)
1. `row_meta` takes no end-token id, yet needs to know whether an answer of max length "reached it without an end token" (C-RM2). The signature cannot
   tell; the lists must carry the end token or the caller must pass it.
2. `rope_positions.bad` checks ids and token types only, not image grids: two images with equal pad counts but different (h, w) share ids but not positions (C-RP5 last case). Safe on the project data (one size per prompt template) but only `--check-rope-every` guards it.
3. `group_requests` on 16 identical consecutive rows with G = 8 raises, though a dataset with a repeated sample could produce it (C-GR2).
4. `patchify` behaviour for an empty list and for float or uint16 tensors is unspecified (C-PX10, C-PX14).
5. Counter reset on `ImageMemo.clear()` is unspecified (C-IM5 asserts deltas only).
6. `count_syncs` nested use and explicit `torch.cuda.synchronize()` counting are unspecified (C-MS2).
7. Spec says "every private TRL name ... is in the Phase A7 start-up check" but A7 is outside this spec; C-CF7 tests the hash/guard behaviour only through the three names listed.
