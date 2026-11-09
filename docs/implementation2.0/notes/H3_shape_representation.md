# H3: How a language model should write a shape as tokens

Date: 2026-11-08. Read-only on the repo. CPU only, no spend. Every number marked "measured" comes from the scripts
in `scratchpad/h3exp/` run on the repo's own LIDC data. Literature facts carry a source link.

## Skills used

| Skill | What it contributed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Catalog check (main thread used the live NVIDIA catalog). Medical matches exist (nv-segment-ct, dicom-series-to-volume, VISTA-3D style CT segmentation), but none covers token representations of shapes. No install. Relevant only to the 3D section as a source of volume masks for a 3D verifier. |
| engineering:code-review | Used to review my own experiment scripts for correctness: found and fixed an output path that wrote into the repo root (file moved out, repo left clean), a numpy 2 `ptp` break, and the unequal perturbation units in `perturb.py` (now stated in the table). |
| engineering:debug | Used its reproduce/isolate/diagnose loop on "why do small nodules fail": isolate the representation (lossless, so not the cause) from localization (1-px sensitivity) and from model shape errors (best-shift test). |
| engineering:system-design | Requirements then trade-offs for the representation choice: token length, legality by construction, error propagation, verifier fit. Section 5 uses its options/trade-off format. |
| engineering:architecture | ADR-style layout for the recommendation (context, options, consequences) in Section 5. |
| engineering:documentation | Report structure: most useful information first (TL;DR), measured tables, paths to reproduce. |
| huggingface-skills:huggingface-best | Main-thread lookup of small VLMs (Qwen3.5-0.8B/2B/4B). Used here to check that the Qwen3.5 tokenizer splits digits the same way as Qwen3-VL (it does: identical token counts). |
| huggingface-skills:huggingface-papers | Paper lookup method (arXiv IDs) for VectorLLM 2507.04664, cadrille 2505.22914, OctGPT 2504.09975, BPT (CVPR 2025), Text4Seg (ICLR 2025). |
| huggingface-skills:hf-cli | Confirmed tokenizers are in the local HF cache (Qwen3-VL-2B-Instruct, Qwen3.5-2B); no downloads needed. |

## TL;DR

1. **The current representation is lossless.** All 3,444 LIDC val targets rasterize back to the rater mask exactly
   (S divides 1000, so every pixel corner is a grid point; 99.7 % of vertices sit on the pixel lattice). It costs
   nothing in IoU. The 0-1000 grid is 10-20x finer than the pixel lattice (1 px = 20 grid units at S = 50), so
   it spends digits on precision that does not exist in the data.
2. **The binding ceilings are not the representation.** Label noise: mean inter-rater IoU 0.683 (0.363 for objects
   under 5 px, 0.644 for 5-9 px). Pixel sensitivity: a 1-px shift of the true mask alone gives IoU 0.455 (under 5 px)
   and 0.665 (5-9 px). Closed form for a disc of diameter d shifted by s px: IoU ≈ 1 - 8s/(πd).
3. **Tokens: 115 per answer, 8.0 per vertex** (fit: 8.01 x vertices + 11.9). Qwen3-VL and Qwen3.5 tokenizers split
   every digit into its own token. Lossless re-spellings cut this 22-69 %. Smooth low-order curves (radial function
   K = 12-16, B-spline with 6-8 control points, Fourier with 2-4 harmonics) need 35-47 tokens and **lose at most
   0.007 IoU against the other raters**. In contrast, coarser min-# polygons hurt small objects badly.
4. **First test: an anchored radial (star) function with K = 16 radii in half pixels.** It needs 47 tokens (-60 %).
   It is legal by construction (100 % measured), canonical by construction (no start vertex, no orientation),
   and splits position (2 numbers) from shape. Its rater-agreement ceiling is 0.680 vs 0.683 lossless, and a
   one-unit token error costs 3x less IoU than in the absolute polygon. Cheap fallback: keep the polygon but
   write pixel units (0..S), which is lossless and -22 % (JSON) to -34 % (flat).
5. **Avoid delta coding without closure.** One wrong delta costs 0.18 IoU on 5-9 px objects (drift), vs 0.06 for an
   absolute coordinate.

## 0. Data and method (measured)

- LIDC val: 3,444 rater masks (S = 50: 3,101; S = 100: 343). Size = bounding-box longest side d in frame px:
  d < 5: 21; 5-9: 911; 10-19: 1,627; 20-39: 774; >= 40: 111. Holes in 1.4 % of masks, more than one part in 0.17 %.
- The truth is the rater mask from `data/prepared/caches/lidc/.../nodule.npz`. Every candidate representation is
  decoded to integer grid polygons and rasterized by the repo's checker (`objdel.geometry.to_masks`, same
  even-odd at pixel centres as the reward). Bitmap-type codes are compared in pixel space directly.
- Tokens: `Qwen/Qwen3-VL-2B-Instruct` tokenizer from the local cache, no special tokens. Qwen3.5-2B gives the
  same count for the current text (115.33 mean for both).
- Model errors: SFT checkpoint completions on the shifted val view (`runs/r1/r1-qwen3vl2b/select/e3/shifted`,
  8,495 parsed samples) scored by mean IoU against each rater mask.
- Scripts (run with `cd <repo> && .venv/bin/python -I <script>`):
  `scratchpad/h3exp/ceilings.py` (main table, about 50 s), `pred_errors.py`, `rows.py`, `quadtree.py`,
  `perturb.py`. Outputs: `ceil_val.json`, `vsother.txt`.

## 1. Error budget

### 1.1 Quantization of the 0-1000 grid

A vertex on a G-level grid over a crop of S px has step q = S/G px. At G = 1000: q = 0.05 px (S = 50) or 0.1 px
(S = 100). The frame policy picks S dividing 1000, so every pixel corner is exactly on the grid and the round trip
mask -> polygon -> mask has zero error (`frame.py`). Measured: all targets reproduce the mask; IoU ceiling = 1.000
in every size bin. **The current representation imposes no IoU ceiling on small objects.**

Coarser grids (measured, rounding the exact target vertices):

| Grid | Step at S = 50 / 100 | IoU vs own mask | IoU vs other raters | Tokens (JSON) |
|---|---|---|---|---|
| 0-1000 (now) | 0.05 / 0.1 px | 1.000 | 0.683 | 115.3 |
| 0-100 | 0.5 / 1 px | 1.000 (legal 99.7 %) | 0.683 | 80.0 |
| 0-50 | 1 / 2 px | 0.990 (0.921 at d >= 40) | 0.682 | 80.0 |

The meaningful lattice is the pixel lattice: S+1 values per coordinate, log2(51) = 5.7 bits at S = 50. The grid
offers 1,001 values, so 95 % of the value space is off-lattice. The SFT model already learned this: only 2.5 % of
its vertices are off-lattice, and snapping them changes mean IoU by less than 0.001 (measured).

General rounding bound: a uniform rounding error δ ~ U(-q/2, q/2) per coordinate gives E|δ| = q/4. For a disc of
diameter d, a mean normal boundary error ē gives IoU ≈ 1 - 4ē/d. So rounding costs about q/d in IoU: 0.005 at q =
0.05 px and d = 10 px, and it is even smaller after rasterization at pixel centres.

### 1.2 Pixel sensitivity (the real budget for small objects)

| Perturbation of the true mask (measured) | all | d < 5 | 5-9 | 10-19 | 20-39 | >= 40 |
|---|---|---|---|---|---|---|
| Inter-rater IoU (label noise) | 0.683 | 0.363 | 0.644 | 0.695 | 0.709 | 0.708 |
| Shift 1 px (x) | 0.767 | 0.455 | 0.665 | 0.775 | 0.861 | 0.895 |
| Shift 1 px (diagonal) | 0.698 | 0.312 | 0.581 | 0.706 | 0.808 | 0.854 |
| Dilate 1 px | 0.729 | 0.447 | 0.630 | 0.734 | 0.822 | 0.861 |
| Erode 1 px | 0.664 | 0.172 | 0.520 | 0.676 | 0.794 | 0.841 |

Closed forms (disc, diameter d px): shift by s: IoU ≈ 1 - 8s/(πd), which gives 0.49 at d = 5 and 0.68 at d = 8
(measured bins 0.455 and 0.665). Uniform boundary error ē: IoU ≈ 1 - 4ē/d. **For IoU >= 0.8 on an 8-px nodule the
boundary must be within 0.4 px.** For comparison, the Qwen3-VL vision grid is patch 16 px and merge 2
(`preprocessor_config.json`). A 512-px image therefore becomes 16 x 16 = 256 visual tokens, and at S = 50 each
visual token covers 3.1 frame px. So the model must place a boundary to about 1/8 of a visual token. That is
the small-object bottleneck, and no output format removes it.

### 1.3 Where the model's error actually is (measured, SFT e3, shifted val)

| bin | n | IoU | after lattice snap | best integer shift within ±3 px | miss (IoU = 0) | median centroid offset (px) | median area ratio | tokens |
|---|---|---|---|---|---|---|---|---|
| all | 8,495 | 0.505 | 0.505 | 0.552 | 7 % | 1.31 | 0.94 | 121 |
| 5-9 | 1,520 | 0.437 | 0.437 | 0.495 | 14 % | 1.06 | 1.07 | 86 |
| 10-19 | 4,341 | 0.543 | 0.543 | 0.592 | 5 % | 1.08 | 0.96 | 98 |
| 20-39 | 2,228 | 0.499 | 0.499 | 0.537 | 5 % | 2.11 | 0.87 | 162 |
| >= 40 | 399 | 0.386 | 0.386 | 0.421 | 9 % | 4.80 | 0.71 | 274 |

A perfect re-centring (oracle shift) recovers only +0.047 IoU. Most of the gap is shape and extent (large nodules
are under-drawn: area ratio 0.71) plus outright misses (14 % at 5-9 px). An anchor-relative format alone therefore
cannot close the gap. Its value lies in shorter answers, structural legality and decoupled error, not in a
higher ceiling.

### 1.4 Vertex budget vs IoU (rate-distortion for polygons)

Min-# approximation (Imai-Iri: shortest path in the shortcut graph, O(n²) after Chan-Chin; min-ε is its dual
by parametric search) is what `simplify.py` already does, potrace-style, with an L∞ tolerance τ. Measured, with
the repo's simplifier at τ in pixels and JSON-1000 spelling:

| τ (px) | mean vertices | IoU vs own mask | IoU vs other raters | d < 5 vs others | 5-9 vs others | tokens | legal |
|---|---|---|---|---|---|---|---|
| exact target | 12.9 | 1.000 | 0.683 | 0.363 | 0.644 | 115.3 | 1.000 |
| 0.5 | 10.2 | 0.972 | 0.679 | 0.332 | 0.636 | 93.5 | 0.999 |
| 1.0 | 5.8 | 0.870 | 0.654 | 0.192 | 0.578 | 58.4 | 0.998 |
| 1.5 | 4.7 | 0.780 | 0.616 | 0.104 | 0.493 | 49.3 | 0.985 |
| 2.0 | 4.0 | 0.682 | 0.564 | 0.000 | 0.418 | 43.7 | 0.961 |
| 3.0 | 3.3 | 0.506 | 0.442 | 0.000 | 0.178 | 38.2 | 0.792 |

On staircase pixel boundaries a min-# polygon with τ >= 1 px collapses small objects (triangles that miss pixel
centres). **Fewer polygon vertices is the wrong way to save tokens for small nodules.** Smooth parametric curves
(Section 2.3) reach the same token count with almost no loss.

## 2. Token efficiency

### 2.1 Tokenizer facts (measured)

- Qwen3-VL-2B and Qwen3.5-2B: one token per digit; `,` is 1 token, `],[` is 1 token, `]],` is 1 token. A vertex
  `[420,380],` costs 8 tokens; `1000` costs 4. Vocabulary 151,669, but digit tokens carry at most log2(10) = 3.32
  bits and the punctuation carries 0 bits (fixed by the grammar).
- BPE merges runs of letters (`EEESSWWN` -> `EE`,`ESS`,`WW`,`N`; `..##` rows merge as well), so letter codes are
  short but have irregular token boundaries that depend on context.

### 2.2 Lossless spellings (measured mean tokens; all bins, then by size)

| Representation | all | d < 5 | 5-9 | 10-19 | 20-39 | >= 40 | vs now |
|---|---|---|---|---|---|---|---|
| JSON 0-1000 (now) | 115.3 | 53.5 | 70.6 | 97.5 | 178.9 | 312.2 | - |
| JSON, pixel units 0..S | 89.5 | 43.1 | 55.9 | 76.1 | 137.2 | 237.3 | -22 % |
| flat `x y` pixel units | 76.5 | 30.1 | 42.9 | 63.1 | 124.2 | 224.1 | -34 % |
| anchor (bbox corner) + relative `x y` | 67.7 | 26.8 | 35.3 | 54.9 | 113.3 | 209.3 | -41 % |
| delta `x0 y0 +dx+dy` (px) | 52.9 | 21.8 | 30.3 | 43.9 | 85.1 | 152.8 | -54 % |
| crack chain code `ESWN` (pixel boundary) | 35.6 | 12.7 | 20.5 | 30.7 | 56.0 | 95.2 | -69 % |
| crack chain RLE `E3S2` | 79.5 | 22.1 | 42.9 | 67.3 | 129.5 | 220.2 | -31 % |
| mask RLE in bbox (COCO-style) | 74.4 | 18.8 | 34.4 | 60.3 | 129.8 | 234.6 | -35 % |
| bbox bitmap rows `.#` | 53.6 | 12.0 | 23.0 | 42.8 | 95.2 | 179.9 | -54 % |
| quadtree E/F/S in bbox | 64.3 | | | | | | -44 % |
| special coordinate tokens (1 token per coordinate, estimated 2/vertex + 4) | ~29.8 | | | | | | -74 % |
| row intervals, delta (lossy: 80 % of masks are row-convex; IoU 0.992) | 60.5 | | | | | | -48 % |

### 2.3 Lossy smooth curves (measured; IoU vs own mask / vs other raters / tokens)

| Representation | vs own | vs others (lossless 0.683) | 5-9 vs others (0.644) | tokens | legal |
|---|---|---|---|---|---|
| radial star K = 8, r in 0.5 px | 0.847 | 0.659 | 0.616 | 25.9 | 1.000 |
| radial K = 12 | 0.910 | 0.676 | 0.632 | 36.3 | 1.000 |
| **radial K = 16** | 0.939 | **0.680** | 0.638 | **46.6** | **1.000** |
| radial K = 24 | 0.967 | 0.682 | 0.644 | 67.4 | 1.000 |
| cubic periodic B-spline, 4 int-px control points | 0.858 | 0.674 | 0.638 | 22.9 | 1.000 |
| B-spline, 6 control points | 0.896 | 0.679 | 0.637 | 35.0 | 0.9997 |
| B-spline, 8 control points | 0.920 | 0.681 | 0.638 | 47.0 | 0.9985 |
| Fourier (EFD-equivalent) N = 2, 10 reals at 0.1 px | 0.900 | 0.683 | 0.647 | 41.0 | 0.9997 |
| Fourier N = 4 | 0.958 | 0.684 | 0.645 | 73.0 | 0.9985 |
| frame grid 16 x 16 labels (Text4Seg-like) | 0.709 | | (d < 5: 0.237 vs own) | 38.6 | n/a |
| frame grid 32 x 32 labels | 0.807 | | | 90.3 | n/a |
| bbox SDF 8 x 8, 0.5 px, bilinear | 0.966 | | | 151.6 | n/a |

Reading: against the label distribution, a smooth 10-16-number curve is as good as the exact staircase. Fourier
N = 3-4 even agrees slightly *better* with the other raters (0.684), because smoothing removes one rater's pixel
noise. Fourier tokens are inflated by decimals and signs (`-12.3` = 5 tokens); with half-pixel integers they would
cost about 30 % less (estimate, not measured). Frame-level coarse grids (Text4Seg, PaliGemma-like) fail on small
objects: a 16 x 16 frame grid has a 3-px cell at S = 50 (d < 5 IoU 0.24). SDF grids are accurate but cost many
number tokens.

Literature anchors: Text4Seg (ICLR 2025) uses 16 x 16 patch labels and row-wise RLE, cutting 583 to 154 tokens
(-74 %) "without compromising performance" (https://arxiv.org/abs/2410.09855). PaliGemma writes a mask as 4 `<loc>`
tokens (1,024 bins) + 16 `<seg>` tokens from a 128-entry VQ-VAE codebook, decoded to 64 x 64 inside the box
(Google blog, https://developers.googleblog.com/gemma-explained-paligemma-architecture/). That is about 16 x 7 =
112 bits for any mask, but it needs a learned decoder, so legality and exactness depend on a network rather than
on geometry. VectorLLM (Qwen3-0.6B/1.7B) adds `[x0]..[x(W-1)]`, `[y0]..` tokens and reports that one token per
coordinate beats digit strings in accuracy, length and compute (https://arxiv.org/abs/2507.04664).

### 2.4 Information-theoretic bounds vs spend (measured on val)

- Lossless, boundary view: order-3 context entropy of the relative crack code (straight/left/right) = 1.07
  bits/step. The mean perimeter is 59.8 steps, so about 64 bits, plus position log2(51²) = 11.3 bits: **about 75
  bits per mask**.
- Lossless, polygon view: empirical joint entropy of vertex deltas = 6.61 bits/vertex x 12.9 + 11.3 = about 96
  bits.
- Spend: 115 tokens, about 0.65-0.83 bits per token. A digits-only code needs at least 75/3.32 = 23 tokens; a
  coordinate-token code (5.7 bits/token) needs at least 13.
- Combinatorial view: the number of digital convex shapes in an m x m box grows like exp(Θ(m^(2/3)))
  (Acketa-Žunić), so the *shape* information of compact blobs grows sub-linearly in size; position and
  rater-specific staircase detail dominate the lossless bit count.
- Lossy view (the one that matters): the rater disagreement (IoU 0.68) is the distortion floor. Codes of 10-16
  numbers (radial, B-spline, Fourier) hit that floor at 35-47 tokens. Everything beyond that encodes one rater's
  staircase.

### 2.5 Error propagation (measured, `perturb.py`): change ONE number by one unit, IoU vs unperturbed shape

| bin | absolute polygon (±1 px on one coordinate) | delta polygon (±1 px on one delta) | radial K = 16 (±0.5 px on one radius) |
|---|---|---|---|
| all | 0.970 | 0.879 | 0.990 |
| d < 5 | 0.844 | 0.644 | 0.981 |
| 5-9 | 0.942 | 0.821 | 0.984 |
| 10-19 | 0.974 | 0.887 | 0.990 |
| >= 40 | 0.997 | 0.948 | 0.996 |

(The radial unit is 0.5 px, half the others. A ±1 px radius error is roughly twice its loss, so about 0.97 at 5-9,
still well above delta's 0.82.) Delta and chain codes drift: one error moves every later vertex and breaks closure.
They are compact but fragile unless closure is enforced during decoding.

## 3. Topology-safe representations (legality by construction)

| Representation | Simple ring guaranteed? | Constraint type | Decoding enforcement |
|---|---|---|---|
| Absolute polygon (now) | No (learned; the checker certifies) | Self-intersection is geometric, not regular | Regex covers syntax only. A custom logits processor could test each candidate lattice value (≤ S+1 per coordinate) against the k existing edges, O(k) per candidate. |
| Radial star r(θ_k), fixed angles, r > 0 | **Yes**: angles strictly increase, so the polygon is star-shaped about the anchor and simple | Fixed length, bounded integers: **regular** | Plain regex in xgrammar or vLLM; the length is fixed, so no early-stop or truncation failure. Canonical by construction (no start vertex or orientation choice). Cannot draw non-star shapes, but the measured cost is small (vs-others 0.680). |
| Row intervals (x-monotone, one run per row, consecutive runs overlap) | **Yes**: simply connected | Local (row i vs row i-1) | Small logits processor; regular if bounded |
| Crack chain code with closure + self-avoidance | **Yes**: a closed self-avoiding lattice path is a Jordan curve | Closure needs two counters (#E = #W, #N = #S): not context-free in general, but **finite-state for bounded S** ((dx, dy) ∈ [-S, S]², about 10⁴ states). Self-avoidance needs a visited bitmap plus a "can still close within the remaining budget" check (BFS on ≤ 51² nodes) | Stateful logits processor, O(1)-O(S²) per step |
| Bitmap / RLE in bbox | No "legality" concept: every bitmap is a valid mask. Topology (holes, parts) is not controlled | Width-bounded rows: regular | Regex per width. Errors stay local (no drift) |
| Fourier / B-spline | No (measured 99.6-99.97 % legal); smooth | Not enforceable cheaply | Post-hoc check |

Implication for constrained decoding: today's regex (`text.py: pattern`) already goes to xgrammar. Radial and
bitmap codes keep that property (syntax = semantics). Polygon and chain codes need a stateful processor to give
legality. That is possible (vLLM custom logits processors), but it costs per-step CPU work in the RL rollout loop.

## 4. 3D and CAD

- **LIDC is already 3D.** A nodule spans several slices. The clinical standard (DICOM RTSTRUCT) is a stack of
  planar closed contours per slice, which is exactly the repo's 2D grammar repeated with a z index. Pairing: the
  same CUDA checker per slice plus the repo's 3D Betti check (`objdel.checker.volume.betti3d`) gives an exact
  verifier with no new geometry. A radial code per slice (16 numbers) plus a shared z-anchor extends cleanly.
  Inter-slice smoothness could use a 2D Fourier/spherical-harmonic code, but star-shapedness in 3D is a stronger
  assumption.
- **Voxels:** RLE per row or octree/Morton codes. OctGPT serializes octree split bits in z-order per depth, plus
  VQ binary codes at the finest level, and reports 13x faster training than prior AR 3D models
  (https://arxiv.org/abs/2504.09975). For small blobs my 2D quadtree measurement says hierarchy does not pay
  (64 tokens vs 36 for chain code), because blobs are boundary-dominated. Octrees pay at large resolutions.
- **Meshes:** PolyGen (vertices then faces via pointers, 8-bit coordinates). MeshGPT (learned residual-VQ face
  tokens). MeshAnything V2 (adjacent mesh tokenization) and EdgeRunner reach about 46-47 % of the vanilla length;
  BPT reaches about 26 % (block-wise indexing gives ≤ 2 tokens per vertex, patch aggregation halves faces again;
  CVPR 2025). BPT states that topology is *not* guaranteed. LLaMA-Mesh writes OBJ text with quantized integer
  coordinates. Validity (manifold, watertight, no self-intersections) is learned, never structural, and the
  exact verifiers (exact-predicate self-intersection tests, winding numbers) are heavier than a 2D rasterizer.
- **CAD:** DeepCAD sketch-and-extrude command sequences with 8-bit quantized parameters. CAD-Recode and cadrille
  write CadQuery Python. cadrille (Qwen2-VL, SFT then online RL "Dr. CPPO" with a programmatic IoU reward) cuts
  invalid programs to 0.0-0.2 % and reaches IoU 90.2 % on DeepCAD point clouds (https://arxiv.org/abs/2505.22914,
  ICLR 2026 oral). That is the same recipe as this repo's (SFT then GRPO-family RL with an exact checker), and the
  closest published precedent.
- **Best pairing with an exact verifier:** a sketch-extrude DSL on an integer grid. Closed loops are the same
  closure constraint as polygons. Extruding a simple planar profile gives a valid solid, and Booleans of valid
  solids stay valid in exact arithmetic. So most legality can be structural, and IoU can be computed exactly
  (rational coordinates). Free-form CadQuery is more expressive but verifiable only through the OpenCascade
  kernel's tolerances. For medical volumes, the contour stack (or radial per slice) is the natural choice.

## 5. Recommendation (ADR-style)

**Context.** The answer format is lossless today, and the model already lands on the pixel lattice. Small-object
failures come from the 1-px budget (≈ 0.4 px boundary accuracy for IoU 0.8 at d = 8; one visual token is 3.1 frame
px) and from shape/extent errors and misses, not from the format. A format can still help in four ways: shorter
answers (cheaper RL rollouts, fewer chances to err), structural legality, canonical form, and errors that do not
propagate.

**Options.**

| Option | Tokens | Ceiling vs other raters | Legality | Error propagation | Change size |
|---|---|---|---|---|---|
| A. Keep JSON 0-1000 | 115 | 0.683 | learned | 0.942 (5-9) | none |
| B. Same polygon, pixel units 0..S (S in prompt), flat syntax | 77-90 | 0.683 | learned | same as A | small: prompt + serializer; reward multiplies by 1000/S |
| **C. Anchor + radial K = 16 (half-px radii)** | **47** | **0.680** | **by construction** | **0.984 (5-9, 0.5-px unit)** | medium: new serializer and decoder to a grid polygon; the checker and reward are unchanged |
| D. B-spline 6-8 integer control points | 35-47 | 0.679-0.681 | 99.9 % | not measured | medium |
| E. One special token per coordinate (VectorLLM / PaliGemma style) | ~30 | 0.683 | learned | as A | large: new embeddings, tokenizer change, vLLM config |
| F. Bitmap rows in bbox | 54 | 0.683 | n/a (any bitmap is a mask) | local | medium; irregular BPE boundaries |
| Reject: min-# with τ ≥ 1 px, frame-level 16x16/32x32 grids, SDF grids, quadtree, chain RLE, delta without closure | | | | | |

**Decision (to test first).** Run C against B as two SFT arms at equal steps. Both decode to integer grid
polygons, so `objdel.geometry.to_masks`, the reward and the eval stay as they are. Anchor = polygon area
centroid in half-pixels; angles θ_k = 2πk/16 starting at +x (screen y down); radius = farthest boundary
crossing, in half-pixels. Prompt: "crop is S px; answer cx cy r0..r15". Grammar: a fixed-length regex of 18
bounded integers. For non-star or multi-part masks (0.17 % multi-part; 1.4 % with holes), fall back to the
polygon. A one-character mode prefix such as `R`/`P` keeps it regular.

**Expected gain.** Certain: -60 % answer tokens (115 -> 47; 5-9 px objects 71 -> 38) and 100 % legal outputs. On
the decode-bound share of each RL step this means shorter rollouts; I have not measured the step-time effect.
IoU: no ceiling loss (-0.003 vs other raters). A plausible 0 to +3 points comes from canonical targets, decoupled
position, and a 3x smaller IoU cost per token error. That range is an estimate, not a measurement. Because
oracle re-centring is worth only +0.047, do not expect the format to fix misses; those need visual resolution
(e.g. a larger upscale so a 6-px nodule covers more than 2 visual tokens) or the shifted-view training already
planned.

**Consequences.** Easier: grammar-exact decoding, no orientation or start canonicalization, fixed-length answers
(no truncation failure). Harder: deep concavities are lost (36 % of staircase boundaries are strictly star-shaped
under a 360-ray test, yet IoU vs others drops only 0.003). Multi-part and holed masks need the fallback. The model
must infer the centroid. Revisit with K = 24 if large nodules (d >= 40) under-fit (their vs-others score is 0.710
at K = 16, so not yet).

### Small offline experiment (no GPU), ready to run

Already built and run on val (results above). To finish the ceiling study before any rental:

1. `ceilings.py train` and the test split: confirm the val tables hold (about 6 min for train on CPU).
2. Extend `perturb.py` to equal physical units (±1 px for radial too) and to two- and three-error draws, so the
   comparison of token-error robustness across formats is fair.
3. Teacher-forced likelihood probe (CPU feasible for a 2B model on about 200 examples, slow but free): score the
   *current SFT checkpoint's* per-token negative log-likelihood on targets rewritten in format B. This shows how far
   the pretrained prior is from each spelling before any new SFT.
4. Gate: adopt C for the SFT A/B only if (1) holds within ±0.005 on train and test, and the fallback rate is under
   2 %.

Paths: report `docs/implementation2.0/notes/H3_shape_representation.md`;
scripts and outputs in `(planning scratch, not kept: h3exp)/`.

Caveats: the size bins use the bounding-box longest side of each rater mask, not the union and not area; the
"0-40 px objects missed 47 %" figure from the brief was not re-derived here. Model errors come from the SFT
checkpoint on the shifted view, not from the final RL model (no saved completions found for it). The d < 5 bin
has only 21 masks (7 model samples).
