# A6: where Phase 2's CPU time goes, and what is worth moving to C++

Written 2026-10-19 by research agent A6 (`experiments/phase2/RESEARCH.md`).

**Question.** Where does CPU time go in Phase 2? The candidates are:
- SFT data loading and collation;
- augmentation that rebuilds certified targets;
- scoring K eval samples per item;
- the GRPO-readiness gate of 256 rollouts in ≤ 250 ms on x86.

And which modules are worth moving to C++, numba, vectorised numpy or the GPU?

**Evidence.**
- Lab scripts: `experiments/phase2/lab/A6/` (run commands are in each module docstring).
- Raw results: `experiments/phase2/lab/A6/results/*.json`.
- Data: `data/prepared` (manifest `cd66d79f…`), read with `.venv/bin/python -I`.
- Stack: the pinned stack of `experiments/phase1/RESULTS.md` §3.
- Profiles: cProfile plus stage timers. I did not use py-spy (it needs root on macOS to sample) or scalene. Stage timers replay the library calls and are checked against the library outputs.

**Measurement caveat.** The 8-core Apple M3 was shared with other agents throughout: load average 8–9 for the first runs and 15–23 for the later ones (each JSON records it).
- Wall times are therefore upper bounds.
- Ratios come from interleaved A/B runs within one process.
- Run to run, the same RIGA-only batch varied 149–652 ms at load 18–20 (`control_k3.json`).
- Every absolute number here must be re-measured on the Vast box (H-1).

## 0. Bottom line

1. **SFT (Phase 2 proper) is GPU-bound once the dataloader has workers.**
   - Measured CPU cost per SFT row on one loaded core: 10–25 ms (PNG decode, then the processor).
   - Estimated GPU time per row: 5–220 ms, depending on model and card (§2).
   - TRL's default `dataloader_num_workers=0` would serialize them; on an H100 that would make the LIDC/SmolVLM step up to ~3× slower.
   - **Decision:** `dataloader_num_workers ≥ 4` and persistent workers. No C++ is needed here.
2. **Augmentation should be precomputed, not done on the fly.**
   - Rebuilding certified targets after a flip or rotation costs 26 ms per RIGA target (172 ms per 6-rater item) and 0.9 ms per LIDC target.
   - The cost does not depend on the transform (§1b).
   - Precomputing all 8 dihedral variants of the train set takes ≈ 15 CPU-minutes. That is cheaper and safer (certified once, hashable) than doing it on the fly.
   - **Keep numba** for the simplifier.
3. **Eval metrics are cheap. Keep them in Python.**
   - K = 64: 145 ms per RIGA item and 8 ms per LIDC item when the loaded caches are passed in.
   - Whole val split: ≈ 17 s.
   - `eval_run` today rebuilds every ExpertCache inside `metrics`, which costs 1.3–5× more. Passing the prepared caches is a one-line fix.
4. **The GRPO reward is the one place where C++ may pay off. It is also an x86 entry-gate risk.**
   - RIGA-only 256-rollout batches on the production (prepared) caches measured **241 ms serial at load 9**, and 300–650 ms at load 15–20.
   - The gate limit is 250 ms. M6's 103.5 ms was a quiet-machine best case.
   - Python threads do not help: 0.8–1.2× (§1d, as in M6).
   - Two remedies, in order:
     1. a fork-after-preload process pool (cheap; H-2);
     2. a single nanobind call per batch with C++ threads (H-3). polycore's own `score_batch` scales 9× at 8 threads (109 → 12 ms).
   - Mixed batches (113 ms at load 9) and LIDC-only batches (37 ms) pass with margin.
5. **New finding: the cache LRU miss costs more than the whole reward.**
   - `ExpertCache.load` costs 32 ms per RIGA object (64 ms per item), 90 % of it in scipy EDTs. LIDC is 1.1 ms.
   - TRL's sampler permutes the dataset, and `cache_loader`'s LRU holds 64 items. A RIGA step therefore pays ≈ 32 misses × 64 ms ≈ 2 s on the main process, about 8× the reward budget.
   - **Decision:** preload every cache at start-up (RIGA train 4.4 GB, LIDC train 0.4 GB; ≈ 35 s serial, and less with threads, since scipy's EDT releases the GIL: 1.6–2.7× on 2 threads).
6. **No GPU port and no free-threaded Python for now.**
   - Free-threading is blocked by OpenCV: it ships abi3 wheels only, and free-threaded CPython supports neither the limited API nor the stable ABI.
   - A GPU path is not worth it for 50–500 px masks: kernel launch and sync overhead, and tie-rule exactness.

## 1. Profiles (measured)

### (a) SFT row building and collation

`prof_collate.py`, load ≈ 9. Batches are drawn at random from `sft_rows(…, "train")`.
- *fetch* = `dict(ds[i])`, which includes the lazy PNG decode by `datasets.Image`.
- *collate* = TRL 1.16 `DataCollatorForVisionLanguageModeling(processor, assistant_only_loss=True)`.
- Processor kwargs come from `configs.processor_kwargs`: SmolVLM `longest_edge` 1024 (RIGA) / 512 (LIDC).

| Model / dataset (bs 16) | fetch ms | collate ms (p95) | tokens / batch | pixel_values | CPU ms / row |
|---|---|---|---|---|---|
| SmolVLM-256M / RIGA | 83 | 323 (388) | 11,776 | 16 × 5 tiles × 3 × 512² | **25** |
| SmolVLM-256M / LIDC | 101 | 57 (78) | 4,339 | 16 × 1 × 3 × 512² | 10 |
| LFM2.5-VL-450M / RIGA | 95 | 138 (153) | 7,952 | 16 × 1024 × 768 | 15 |
| LFM2.5-VL-450M / LIDC | 90 | 109 (123) | 6,699 | 16 × 1024 × 768 | 12 |
| Qwen3-VL-2B / RIGA | 78 | 98 (115) | 9,309 | 15,616 × 1536 | 11 |
| Qwen3-VL-2B / LIDC | 105 | 96 (170) | 7,364 | 16,384 × 1536 | 13 |

Hot spots (cProfile of one 16-row collate):
- **SmolVLM / RIGA:** `torch._C._nn._upsample_lanczos2d_aa` takes 181 of 407 ms. That is the 500 → 1024 px Lanczos resize before tiling, the single largest CPU item in Phase 2's SFT path. The Idefics3 `_preprocess` itself takes 66 ms, tensor `.to` 46 ms.
- **LFM / Qwen:** `apply_chat_template` plus the tokenizer's `flatten` / tensor `unbind` loops (≈ 9k Python calls per batch) take 20–35 ms. Bilinear/bicubic resizes take 17 ms.
- **PNG decode** (`decode.py`, load 18): RIGA 500² frames (267 kB) and LIDC 512² frames (169 kB), median ms:

  | Decoder | RIGA | LIDC |
  |---|---|---|
  | PIL | 9.4 | 6.9 |
  | OpenCV | 7.9 | 5.8 |
  | torchvision `decode_png` | 7.7 | 5.7 |
  | copy of a RAM-resident array | 0.09 | 0.08 |

  The three decoders are within 20 % of each other because all are zlib-bound. Published libpng/libspng speeds of 200–300 MP/s on an idle Zen 4 core predict ≈ 1–1.5 ms per 0.26 MP frame (A6 web notes, §Sources). Ours is ×5 because of load.

`sft_rows` (one-off per run):

| Split | Rows | Wall time | Time per row |
|---|---|---|---|
| RIGA train | 3,281 | 23 s | 7.1 ms |
| LIDC train | 23,459 | 33 s | 1.4 ms |

On RIGA, 14.2 of 22 s go to scipy's `euclidean_feature_transform`. `_voting` loads full ExpertCaches only to read `rater_sigs`; the second-largest item is the PNG hash check. This is a harmless start-up cost (a fix is listed in §3).

### (b) Target rebuild after flip / rotate (`prof_targets.py`, load 15)

Masks are the frame-size rater stacks stored in the caches. Check: the identity transform reproduces **144 / 144 RIGA and 974 / 974 LIDC** stored target texts.

| | RIGA (24 items × 6 raters) | LIDC (300 items, 3.2 raters on average) |
|---|---|---|
| `build_target` per target, median / p95 / max | 25.5–26.2 / 35–36 / 40–186 ms | 0.91–0.94 / 1.7–1.8 / 2.5–3.0 ms |
| per item (all raters), median | 168–178 ms (the same for identity, hflip, vflip, rot90, rot180, transpose) | 2.9–3.1 ms |
| `ExpertCache.build` of the augmented item (GRPO would need it) | **280 ms** (p95 473) | 1.5 ms |
| tolerance reached after augmentation | τ = 3 on 52 %, down to 0.75 | τ = 3 on 100 % |

cProfile, RIGA, per item (207 ms under the profiler):

| Function | ms per item | Share |
|---|---|---|
| numba `simplify` (≈ 33 calls; the ladder re-simplifies an object at each rung) | 86 | 42 % |
| `np.nonzero` inside `polygonize.boundary_edges`, on the **full 500² frame** 4× per ring set | 40 | 19 % |
| OpenCV CCL (`instance_topology` in `_kept` / `signature`) | 23 + 6 | 14 % |
| everything else (polycore, canonicalize, …) | < 10 each | |

Cold start (`cold_start.py`, load 19):
- The first `build_target` with an empty numba cache takes **20.2 s** of JIT. With a warm cache, the first call takes 1.3 s and later calls 0.08 s.
- Every dataloader worker on a fresh box would pay this unless the cache is warmed once in the parent before workers start.

### (c) Eval metrics (`prof_eval.py`, load 11)

Samples come from the K3 generator: rater masks jittered, shifted, grown or exact, then traced. Distinct shares: RIGA 0.95 (K = 8) and 0.75 (K = 64); LIDC 0.85 and 0.59.

| `metrics()` per item, median (p95) | K = 8, loaded caches | K = 8, caches built (eval_run path) | K = 64, loaded caches | K = 64, caches built |
|---|---|---|---|---|
| RIGA (10 val items) | 48 (141) ms | 253 (351) ms | 145 (179) ms | 339 (453) ms |
| LIDC (60 val items) | 1.6 (2.8) ms | 3.5 (7.4) ms | 8.4 (13.6) ms | 10.7 (17.5) ms |

- Hot spots at K = 64 (RIGA): `_pairwise_iou` takes 52 ms per item: `np.unique` over packed rasters (argsort 13 ms plus structured-dtype `_promote_fields` 9 ms), then a float64 matmul. The per-sample `score_objects` path makes up the rest.
- Whole val split (58 RIGA + 1,084 LIDC items) at K = 64: ≈ 8 s + 9 s with loaded caches, ≈ 20 s + 12 s via eval_run.
- This is small next to generating 73k samples.

### (d) Reward scoring of 256 rollouts (`prof_reward.py`, load 8.9; prepared caches; `cv2.setNumThreads(0)`)

Wall time per batch of 32 prompts × 8 completions, `advantages_v4(score_batch(jobs, pool))`, median (max):

| Workload | Characters per completion | serial | 2 threads | 4 threads | 8 threads | Per completion p50 / p99 / max | C++ floor `polycore.score_batch` (raster + validity + vote intersection), 1 / 8 threads |
|---|---|---|---|---|---|---|---|
| mixed (40 % RIGA) | 422 | **113 (118)** | 149 | 142 | 146 | 0.14 / 1.69 / 4.9 ms | 19.6 / 5.7 ms (360 objects) |
| RIGA only | 851 | **241 (257)** | 239 | 239 | 286 | 0.89 / 2.23 / 15.3 ms | 109 / 12.0 ms (512 objects) |
| LIDC only | 119 | 37 (81) | 38 | 52 | 51 | 0.15 / 1.39 / 18.1 ms | 2.0 / 1.3 ms |

- K4 (p99 per rollout ≤ 5 ms) passes everywhere.
- In CPU time (`cpu_k3.py`, load 15–17), RIGA-only costs 296 ms on the prepared caches and 299 ms on M6's fixture caches. The data is not the cause of the gap to M6's 103.5 ms; machine load is (E-cores and shared caches inflate even CPU time).

Stage split, replaying `score_objects` stage by stage (ms per RIGA-only batch; results identical to the library's):

| Stage | ms | Who holds the GIL |
|---|---|---|
| `signature` (two OpenCV Spaghetti labelings + crop/pad + numpy glue; the slow size-filter path in 0.2–0.8 % of objects) | 88 | CCL releases it; the glue does not |
| IoU against the same-signature raters (`&` + `count_nonzero` on the ROI, ≈ 6 per object) | 57 | numpy (mostly held: many small calls) |
| `band_depth` (boolean-mask gathers + float64 sums) | 39 | numpy |
| ROI: `cv2.boundingRect` on the full frame + slicing + area | 30 | cv2 releases it; 28 µs per call |
| `polycore.rasterize` (C++) | 17 | released |
| `parse` (regex fullmatch + `json.loads`) | 17 | held |
| `polycore.pack` + `validate` (Python wrapper `_as_i32` ≈ 10) | 18 | partly |
| `rater_answers` (recomputed for every completion) | 3 | held |

**Prototype H1 (`proto_fused.py`, `micro_score_mask.py`).** I fused score_mask's pixel work (bbox, area, per-rater and vote intersections, band depth) into one numba `nogil` kernel as a proxy for C++.
- **Exactness:** 0 mismatched NUMERIC terms over 6 batches × 256 completions. The band sums are exact in float64 in any order: every distance is a float32 ≥ 1 with fewer than 2¹⁸ terms below 2¹⁰ each, all multiples of 2⁻²³. So a C++ port can be bit-identical (K10 hashes unchanged).
- **Speed, interleaved per object at load ≈ 20, median µs:**

  | Dataset | lib `score_mask` | ROI kernel | `signature` on the ROI | → kernel + signature |
  |---|---|---|---|---|
  | RIGA | 573 | 171 | 181 | **352 (1.6×)** |
  | LIDC | 92 | 13 | 29 | **42 (2.2×)** |

  A naive whole-frame scan without an ROI was *slower* than numpy (scalar loop, 4× the pixels).
- **Threads in this hybrid** (Python glue + GIL-free kernels): 4 threads give 1.3–1.9× over the hybrid's own serial, but no better than the library's serial at that load. Amdahl's law limits it: parse, `rater_answers`, the wrappers and the signature glue keep the GIL.
- **Conclusion:** the gain needs the *whole* per-completion path in one GIL-free call, which polycore's 9× thread scaling shows to be achievable. Moving pieces piecemeal is not enough.

**Cache load (`prof_cache.py`, load 17).** `ExpertCache.load` per object, median (p95):

| Dataset | Median (p95) | Stages | Memory per object |
|---|---|---|---|
| RIGA (6 × 500 × 500) | 32 ms (90) | npz inflate 1.8, SHA-256 0.26, `_band` 33 (two scipy EDTs on the full frame: 8 and 20 ms); OpenCV's exact EDT on the same mask: 2.4 ms | 4.05 MB |
| LIDC | 1.1 ms | | 0.05 MB |

## 2. CPU budget against the GPU step (derived, to be replaced by the Vast measurement)

Training FLOPs per row ≈ 6·N·tokens for trained parameters, plus 6·N_vis·patches for the vision tower (4·N for LoRA: forward plus input gradients).

Inputs:
- token and patch counts measured above;
- SmolVLM-256M: 93M SigLIP + 135M LM; 1,024 patches per 512 tile (SmolVLM paper);
- LFM2.5-VL-450M: 86M + 350M (model card);
- Qwen3-VL-2B: LoRA, LM ≈ 1.7B, vision tower ≈ 0.4B (**assumed**);
- peak dense BF16: 4090 165 TFLOPS (Ada whitepaper), A100 312 (datasheet), H100 SXM 989;
- MFU 15 %, a cautious value for sub-0.5B models (lmms-engine reports 30–40 % for 7–8B VLMs with packing).

| Model / dataset | FLOPs per row | GPU ms per row, 4090 @ 15 % | GPU ms per row, H100 @ 15 % | CPU ms per row, measured (no augmentation) | Workers needed to hide CPU on H100 |
|---|---|---|---|---|---|
| SmolVLM / RIGA | 3.5 T | 141 | 24 | 25 (+26 with on-the-fly targets) | 2 (3) |
| SmolVLM / LIDC | 0.8 T | 32 | 5.3 | 10 | 2–4 |
| LFM / RIGA | 1.6 T | 63 | 11 | 15 | 2 |
| LFM / LIDC | 1.4 T | 57 | 9.5 | 12 | 2 |
| Qwen3-VL LoRA / RIGA | 5.5 T | 222 | 37 | 11 | 1 |
| Qwen3-VL LoRA / LIDC | 4.8 T | 194 | 32 | 13 | 1 |

- The CPU side is 0.05–2× the GPU side per row. With `num_workers = 0` (the TRL default) it adds to every step. With ≥ 4 workers it is hidden on every listed card.
- SmolVLM's vision tower dominates its FLOPs on RIGA (2.9 of 3.5 T): the 1,024 px setting buys 5 tiles.
- **GRPO:** TRL calls the reward synchronously on each process after generation (async only for `async def` rewards). 0.1–0.25 s per step against multi-second generation of 256 × ≤ 1,536 tokens is a few percent. The 250 ms gate is therefore a guard against the reward becoming visible, not the step's bottleneck. The 2 s of cache misses (§0.5) would be visible.

## 3. Research per hot spot (2025–2026, primary sources)

| Hot spot | State of the art / facts | Bearing on us |
|---|---|---|
| Image processors | transformers v5 has torchvision backends, "up to 33× faster" than PIL on tensor batches, with `device=` for GPU preprocessing. Idefics3/SmolVLM fall back to PIL for Lanczos unless torchvision > 0.27 (`backend="torchvision"`). Qwen3-VL ships `Qwen2VLImageProcessorFast` | We already run torch's Lanczos (the profile shows `_upsample_lanczos2d_aa`). The remaining cost is the resize itself, so remove it offline (H-7) or run the processor on the GPU |
| TRL data path | `dataloader_num_workers` defaults to 0. Image datasets skip `prepare_dataset` and preprocess in the collator. `remove_unused_columns=False` is needed with custom columns | Set workers ≥ 4, `dataloader_persistent_workers`, `dataloader_pin_memory` |
| Image decode | PNG decoders run at 200–375 MP/s on Zen 4 (≈ 1 ms per frame). nvJPEG/`decode_jpeg(device="cuda")` is JPEG only. DALI 2.3 (2026-08). pillow-simd is sdist-only and lags Pillow 13 | PNG must stay (lossless, `image_sha256`). Decode once into RAM (RIGA 0.4 GB; LIDC 1.9 GB as single-channel) if workers are scarce. Not DALI |
| Dataloader frameworks | torchdata DataPipes were removed in 0.10, and `torchdata.nodes` is the successor (0.11, Feb 2025). WebDataset 1.0.2 and MosaicML Streaming 0.13 target multi-node streaming | 1.5 GB and ≈ 10k images: a map-style dataset in RAM. None of these frameworks |
| EDT | scipy (Maurer) releases the GIL (`NPY_BEGIN_THREADS` in `ni_morphology.c`; measured 1.6–2.7× on 2 threads). OpenCV's exact FH path is 3–8× faster here but nondeterministic on x86 IPP (opencv#24082, M9). cupyx PBA+ is exact on the GPU. seung-lab `edt` 3.1 is multi-threaded with cp314t wheels | For cache build/load: preload first. If needed, a ~100-line exact Felzenszwalb–Huttenlocher EDT in `polycore.hpp` (integer squared distances, one sqrt per pixel, correctly rounded to float32) would be deterministic on every CPU (H-8) |
| CCL | OpenCV Spaghetti / Spaghetti4C (Bolelli et al. 2020). The parallel path only runs at ≥ 2 × nCPU rows, so our masks run serially anyway. BKE on the GPU takes ≈ 0.45 ms per small image excluding transfers (Allegretti et al.). cc3d is 3-D oriented | Keep OpenCV on the CPU. In a C++ batch scorer, link the same algorithm or a two-pass union-find with hole counting, checked against `instance_topology` |
| Polygon simplification | numba ≈ C (M5: about 300× pure Python). `numba.pycc` is deprecated with no free-threaded AOT. numba 0.65+ supports 3.14t | Keep numba. Precompute instead of porting |
| JSON parse | orjson is 4–6× stdlib on large documents but holds the GIL and has no free-threaded wheels. msgspec ≈ orjson | parse is 17 ms of 241 (7 %). Not worth a dependency; inside a C++ scorer, a hand parser of the fixed grammar |
| Rasterization on the GPU | nvdiffrast is triangle-only, point-sampled with OpenGL conventions | It would break the exact top-left rule (K1). Rejected |
| GPU kernels for tiny masks | 3.4–9.6 µs launch + sync per kernel (CUDA Graphs blog), against tens of µs of CPU work per op | Only a single batched kernel over all 512 masks could pay, and the GPU is busy training. Rejected for Phase 2–3 |
| nanobind | `nb::call_guard<nb::gil_scoped_release>` (polycore already releases the GIL: `_polycore.cpp` lines 86/110/157). `FREE_THREADED` flag since 2.2. ~3× lower call overhead than pybind11 | A batch API (one call per 256 rollouts, internal `std::thread` pool) is the proven pattern: 9× at 8 threads |
| Free-threaded CPython | PEP 779: officially supported, optional in 3.14 (~10 % single-thread cost on Linux). numpy, scipy, torch, numba and pillow have cp314t wheels; **OpenCV has abi3 only**, and free-threaded builds cannot load abi3 | Blocked until OpenCV ships cp314t (or CCL/boundingRect leave cv2) |
| RL frameworks | verl's RewardLoopManager runs rule rewards in `reward.num_workers` workers, per sample as rollouts finish. OpenRLHF uses remote reward servers via Ray. TRL calls reward functions inline, awaiting `async def` ones concurrently | The process-pool pattern (H-2) is what verl does. Our TRL path can do the same inside `TRLReward` via `pool=` |

## 4. Decisions per module

| Module / step | Today | Decision | Expected effect | Engineering cost | Blocks the x86 gate? |
|---|---|---|---|---|---|
| `train.dataset` rows + TRL collator | Python, `num_workers` 0 | **Keep Python.** `dataloader_num_workers=4–8`, persistent workers, pin memory; warm the numba cache in the parent | CPU hidden; up to ~3× step time on H100 / LIDC versus 0 workers (derived) | config only | no |
| SmolVLM RIGA Lanczos resize (181 of 407 ms per 16 rows) | torch CPU, per batch | **Pre-resize offline** to the processor's longest edge (H-7), or run the processor with `device="cuda"` in the main process | −40 % SmolVLM/RIGA collate | small (prepare step plus a hash) | no |
| `sft_rows._voting` (14 of 22 s on EDTs) | loads full caches | Store voting rater indices in the prepared record, or read `rater_sigs` without `_band` | −60 % start-up | small, in a Phase 2 code change | no |
| Augmentation targets (`targets.build_target`, numba) | none (I7) | **Precompute the dihedral variants offline** (8 × train ≈ 15 CPU-minutes), each variant with its own `target_set_hash` and K5/K6/R1 check. Keep numba. Crop to the bbox before `polygonize` (H-9) | 0 runtime cost; −15–20 % build | small | no |
| On-the-fly augmentation (if A4 wants continuous transforms) | — | Allowed only in workers: 26 ms per RIGA row plus ≥ 1 extra worker | — | — | no |
| GRPO augmentation (caches) | — | **Not online.** `ExpertCache.build` takes 280 ms per RIGA item; precompute or skip | — | — | — |
| `eval.metrics` | Python | **Keep.** Pass the prepared caches into `metrics` from `eval_run` / `sample_eval` | 1.3–5× faster eval (RIGA 339 → 145 ms per item at K = 64) | one line, outside Phase 1 code | no |
| `reward.cache_loader` (LRU 64) | load on a miss: 64 ms per RIGA item | **Preload all caches at start-up** (`maxsize ≥ dataset`; thread pool, since the EDT releases the GIL) | removes ≈ 2 s per RIGA step | small | **yes**: it is per-step reward latency |
| `reward.score_batch` (Python + numpy + OpenCV + polycore) | serial | Phase 2: **keep**, and measure on Vast (H-1). If RIGA-only is > 250 ms there: first **fork-after-preload process pool** (H-2); if that fails or risks CUDA/fork problems, a **C++ batch scorer** in `cpp/polycore` (H-3) | H-2: ≈ 3× on 4 workers (expected); H-3: ≥ 1.6× serial (measured proxy), ≈ 6–9× at 8 threads (polycore's scaling) | H-2 ≈ 50 lines; H-3 ≈ 600–900 lines C++ (parser fast path, CCL with holes + size filter, fused terms) plus differential tests against the Python reference: 1–2 weeks | **yes, the main risk** |
| `reward.cache` EDT (scipy) | 8–20 ms per 500² EDT | Keep scipy, which is deterministic (M9 F1). A C++ EDT only if caches are rebuilt online (H-8) | ≥ 3× on build/load | ≈ 100 lines + tests | no (after preload) |
| `geometry.parse` | regex + json | Keep (7 % of the reward) | — | — | no |
| GPU (rasterizer / EDT / CCL) | — | **Rejected** (exactness, launch overhead, GPU busy) | — | high | — |
| Free-threaded 3.14t | — | **Rejected for now** (cv2 abi3) | — | — | — |

## 5. Hypotheses to test

Each hypothesis has a sandbox design and a pass criterion. "Vast" means the first rented box, before any training step, with nothing else running. Report the CPU model and `nproc`.

| # | Hypothesis | Test design | Pass criterion |
|---|---|---|---|
| H-1 | The x86 reward gate holds on production caches | On the Vast box: `python -I -m experiments.phase2.lab.A6.prof_reward --out …` and `cpu_k3` (mixed, RIGA-only, LIDC-only; serial and 2/4/8 threads), 3 runs | RIGA-only **and** mixed serial median ≤ 250 ms, p99 per rollout ≤ 5 ms. If RIGA-only fails, run H-2/H-3 before GRPO |
| H-2 | A process pool forked after preloading caches (copy-on-write, jobs carry only text + key) scales the reward | Prototype `TRLReward(pool=…)` with a `fork`-context pool of 4 workers created before CUDA init. Workers look caches up by key from the parent's preloaded dict. Ten 256-rollout RIGA-only batches | ≥ 3× over serial wall, the SHA-256 of terms + advantages identical to serial (K10), RSS growth per worker < 200 MB |
| H-3 | A C++ batch scorer (parse fast path → raster → fused terms → CCL signature, internal threads, one nanobind call per batch) beats Python | Sandbox header beside polycore. Differential test against `score_completion` on 10 batches × 3 workloads plus the M6 fixture adversaries | 0 mismatched NUMERIC terms; serial ≥ 1.5× and 8 threads ≥ 5× faster than Python serial on 256 RIGA-only rollouts |
| H-4 | Preloading every cache at start-up is cheap and removes the per-step miss cost | Load all train caches with a `ThreadPoolExecutor(4)`; then time `TRLReward._terms` on 32 cold prompts vs preloaded | Preload ≤ 60 s, RSS ≤ 6 GB; reward call on cold prompts ≥ 1.5 s slower than preloaded (confirms the cost); preloaded within 10 % of K3 |
| H-5 | Dataloader workers hide the SFT CPU path | Vast 4090: SmolVLM LIDC and RIGA, 50 SFT steps, `dataloader_num_workers` ∈ {0, 2, 4, 8}; log data-wait time per step (time from batch request to arrival) | With 4 workers: data wait < 5 % of step time, step ≥ 1.2× faster than with 0 workers (H100: ≥ 2×) |
| H-6 | Offline dihedral augmentation is cheap and keeps every certification | Build 8 variants of all train targets (6 processes); run the R1 checks (parse, legal, signature, IoU ≥ 0.99, F5 ≥ 0.95) per variant | ≤ 20 CPU-minutes; 100 % certified per variant; identity variant hash = `1af819f4…`'s train subset |
| H-7 | Pre-resizing RIGA frames to SmolVLM's 1024 px edge removes the Lanczos cost | Resize offline with the processor's own resample, store PNGs, compare `pixel_values` and collate time | pixel_values identical (or max abs diff ≤ 1/255 if the processor still resamples); collate ≥ 40 % faster |
| H-8 | An exact FH EDT in C++ matches scipy bit-for-bit and is faster | `edt_fh` in a sandbox header; compare with `objdel.reward.cache.edt` on 103k fixture masks plus random masks, arm64 and x86 Docker | 0 differing float32 values on both architectures; ≥ 3× faster on 500² frames |
| H-9 | Cropping masks to their bbox before `polygonize` / `connectedComponents` in `build_target` cuts target cost without changing texts | Patch in the sandbox, rebuild all train/val targets | target-set hash unchanged (`1af819f4…`); RIGA median per target ≥ 15 % faster |
| H-10 | Eval stays negligible | `eval_run` with prepared caches on full val at K = 64 (synthetic samples) | ≤ 60 s wall single-process; outputs identical to the cache-building path |

## 6. Rejected alternatives (with reasons)

- **Python threads for the reward.** Measured at 0.8–1.2× here and in M6. The GIL is held by the numpy/regex glue between GIL-free kernels.
- **Porting `score_mask` piecemeal to a GIL-free kernel.** Only 1.6× serial on RIGA (prototype) and no thread scaling without the parse and signature glue. Only worth doing as part of H-3.
- **Porting the numba simplifier to C++.** numba is already compiled (≈ 300× pure Python, M5). The real cost is the ladder's repeated rungs and the full-frame polygonize scan, which are fixed in Python (H-9) or made irrelevant by precomputation.
- **orjson / msgspec for parse.** 7 % of the reward, a new dependency, and no free-threaded wheel for orjson.
- **DALI, WebDataset, Streaming, torchdata.nodes.** The dataset fits in RAM, and PNG decode is ≈ 1 ms on an idle core.
- **GPU rasterization / EDT / CCL.** Exactness (top-left rule), launch and sync costs at 50–500 px, and the GPU is the busy resource.
- **Free-threaded Python 3.14t.** OpenCV wheels are abi3-only.

## 7. Sources

Repo evidence:
- `experiments/phase1/RESULTS.md` §2 (K3/K4), §6 (x86), §8 (entry gates);
- `experiments/phase1/impl/M1_polycore.md` (score_batch scaling, GIL release);
- `M5_targets.md` (numba vs pure Python, build time);
- `M6_reward.md` (K3 serial vs threads, process pool with caches);
- `src/objdel/_polycore.cpp` (`gil_scoped_release` at lines 86, 110, 157);
- `src/objdel/train/dataset.py` (`cache_loader` LRU 64, `_voting`).

Web (one URL per claim):
- SmolVLM processor defaults, `longest_edge` 2048, splitting, 64 tokens per tile: https://huggingface.co/HuggingFaceTB/SmolVLM-256M-Instruct/raw/main/preprocessor_config.json , https://huggingface.co/HuggingFaceTB/SmolVLM-256M-Instruct/raw/main/processor_config.json , https://huggingface.co/docs/transformers/main/en/model_doc/smolvlm , https://arxiv.org/pdf/2504.05299
- Qwen3-VL-2B patch 16 / merge 2, fast image processor: https://huggingface.co/Qwen/Qwen3-VL-2B-Instruct/raw/main/preprocessor_config.json
- LFM2.5-VL-450M (SigLIP2 NaFlex 86M, ≤ 256 image tokens): https://huggingface.co/LiquidAI/LFM2.5-VL-450M , https://www.liquid.ai/blog/lfm2-vl-efficient-vision-language-models
- MFU of 7–8B VLM training, 30–40 % with packing: https://lmms-engine.readthedocs.io/en/latest/_sources/reference/mfu.md.txt ; padding removal: https://lmms-engine.readthedocs.io/en/latest/models/qwenvl.html
- nanoVLM, ~6 h on one H100: https://huggingface.co/blog/nanovlm
- RTX 4090 165.2 TFLOPS dense BF16 (FP32 accumulate): https://images.nvidia.com/aem-dam/Solutions/geforce/ada/nvidia-ada-gpu-architecture.pdf (summarised at https://aiwiki.ai/wiki/nvidia_geforce_rtx_4090) ; A100 312 TFLOPS BF16 dense: https://www.nvidia.com/content/dam/en-zz/Solutions/Data-Center/a100/pdf/nvidia-a100-datasheet-us-nvidia-1758950-r4-web.pdf ; H100: https://www.nvidia.com/en-us/data-center/h100/
- TRL SFT defaults (`dataloader_num_workers`, `skip_prepare_dataset` for images, packing / padding-free, `max_length=None` for VLMs): https://huggingface.co/docs/trl/main/en/sft_trainer
- TRL GRPO reward-function call pattern, async rewards: https://huggingface.co/docs/trl/main/en/grpo_trainer
- transformers image-processor backends (torchvision default, "up to 33×", `device=`, Idefics3 Lanczos exception): https://huggingface.co/docs/transformers/main/image_processors
- JPEG decode benchmarks: https://arxiv.org/html/2501.13131v1 ; multi-worker caveat: https://arxiv.org/pdf/2605.08731 ; nvJPEG via torchvision: https://docs.pytorch.org/vision/stable/generated/torchvision.io.decode_jpeg.html ; A100 hardware JPEG: https://developer.nvidia.com/blog/leveraging-hardware-jpeg-decoder-and-nvjpeg-on-a100/
- PNG decoder throughput on Zen 4 (libpng 209, libspng 300 MP/s; forum benchmark, search snippet only): https://redlib.hackliberty.org/r/rust/comments/1ha7uyi/comment/m1978ve
- DALI 2.x, free-threaded support: https://docs.nvidia.com/deeplearning/dali/archives/dali_2_0_0/release-notes/index.html ; pillow-simd / Pillow / DALI / WebDataset / Streaming / torchdata release state: https://pypi.org/project/pillow-simd/ , https://pypi.org/project/nvidia-dali-cuda120/ , https://pypi.org/project/webdataset/ , https://pypi.org/project/mosaicml-streaming/ , https://pypi.org/project/torchdata/
- torchdata.nodes: https://docs.pytorch.org/data/beta/getting_started_with_torchdata_nodes.html ; Streaming design: https://www.databricks.com/en/blog/mosaicml-streamingdataset
- scipy EDT releases the GIL (`NPY_BEGIN_THREADS`): https://raw.githubusercontent.com/scipy/scipy/main/scipy/ndimage/src/ni_morphology.c ; nogil labeling: https://raw.githubusercontent.com/scipy/scipy/main/scipy/ndimage/src/_ni_label.pyx
- OpenCV CCL / EDT flags and the parallel-CCL row threshold: https://raw.githubusercontent.com/opencv/opencv/4.x/modules/imgproc/include/opencv2/imgproc.hpp ; GIL release in the bindings: https://raw.githubusercontent.com/opencv/opencv/4.x/modules/python/src2/cv2_util.hpp ; IPP EDT nondeterminism: https://github.com/opencv/opencv/issues/24082
- cupyx exact EDT (PBA+): https://docs.cupy.dev/en/latest/reference/generated/cupyx.scipy.ndimage.distance_transform_edt.html ; cupyx label sync: https://docs.cupy.dev/en/latest/reference/generated/cupyx.scipy.ndimage.label.html
- GPU BKE timings: https://iris.unimore.it/bitstream/11380/1179616/1/2018_TPDS_Optimized_Block_Based_Algorithms_to_Label_Connected_Components_on_GPUs.pdf ; YACCLAB: https://github.com/prittt/YACCLAB ; cuCIM: https://developer.nvidia.com/blog/cucim-rapid-n-dimensional-image-processing-and-i-o-on-gpus/
- Kernel launch overhead: https://developer.nvidia.com/blog/cuda-graphs/
- nvdiffrast conventions: https://nvlabs.github.io/nvdiffrast/
- seung-lab edt / cc3d: https://github.com/seung-lab/euclidean-distance-transform-3d , https://pypi.org/project/edt/ , https://seunglab.org/connected-components-3d/benchmarks/
- numba pycc deprecation: https://numba.readthedocs.io/en/latest/reference/deprecation.html ; free-threading in numba: https://numba.readthedocs.io/en/0.63.0/release/0.63.0-notes.html , https://numba.readthedocs.io/en/latest/release/0.65.0-notes.html
- orjson: https://github.com/ijl/orjson ; msgspec: https://msgspec.dev/benchmarks.html
- nanobind GIL management: https://nanobind.readthedocs.io/en/latest/functions.html ; free-threaded: https://nanobind.readthedocs.io/en/latest/free_threaded.html ; ndarray: https://nanobind.readthedocs.io/en/latest/ndarray.html ; benchmarks: https://nanobind.readthedocs.io/en/latest/benchmark.html
- PEP 779 (free-threaded supported in 3.14): https://peps.python.org/pep-0779/ ; free-threaded builds cannot use the limited API / stable ABI: https://docs.python.org/3.14/howto/free-threading-extensions.html ; opencv-python-headless wheels (abi3): https://pypi.org/project/opencv-python-headless/
- verl reward loop: https://verl.readthedocs.io/en/latest/advance/reward_loop.html ; OpenRLHF remote rewards: https://openrlhf.readthedocs.io/en/latest/agent_training.html

Not found: a published tokens/s or step time for SFT of any ≤ 2B VLM on a 4090/A100/H100 (hence the derived §2 table, which the Phase 2 "Cost" entry gate replaces), and a head-to-head of a Python thread pool against an internal C++ pool driving GIL-free code (our polycore numbers are the evidence).
