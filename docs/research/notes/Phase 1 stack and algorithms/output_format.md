# VLM segmentation output format and SFT target construction (2023–2026)

Scope: how VLMs are trained to emit segmentation geometry (polygons, contours, mask tokens), what that implies for our
output format, and how we should build SFT targets for SmolVLM-256M/500M or Qwen3-VL-2B before GRPO with the exact
geometry checker. Research date 2026-10-18.

Evidence labels used below:
- **[paper]**: read in the primary PDF this session (text extracted with pdftotext and read or grepped).
- **[secondary]**: from a search snippet or a third-party page only. Not checked against the primary source.
- **[local]**: measured in this repo on 2026-10-18 (method is described inline; the scripts were scratch files and were not committed).
- **[memory]**: background knowledge, not re-fetched this session. Treat as unverified.

---

## 1. Survey of output schemes: tokenization, ordering, holes and multi-part, vertex counts, simplification, accuracy

### Takeaway
Two families exist:
- **Polygon-as-text.** SimpleSeg, VistaLLM, PolyFormer, VectorLLM, Florence-2, Rex-Omni, FloorplanVLM and Groundbench all use it. They converge on a few choices: quantize coordinates to about 1000 bins relative to the image or crop, trace clockwise, start from the top-left-most vertex, and use roughly 16–40 vertices per object.
- **Mask tokens.** PaliGemma VQ-VAE, HiMTok, UFO, Text4Seg/GeoGround descriptors and LISA's [SEG] token belong here. These handle holes for free but need a decoder or a new vocabulary.

None of the polygon-as-text VLM papers I found represents holes explicitly. SimpleSeg lists "object with holes" as a failure case. Our `outer`/`holes` schema therefore has no direct precedent.

### Cited Findings

**SimpleSeg (Moonshot AI + Nanjing U., arXiv 2601.19228, Jan 2026) [paper]** — [arXiv](https://arxiv.org/abs/2601.19228), project/code: https://simpleseg.github.io/
- Format: `<mask> : [[x1, y1], [x2, y2], ..., [xV, yV]]`. Coordinates are normalized to [0, 1] and printed with 3 decimals, for example `[[[0.819, 0.453], [0.779, 0.453], ...]]`. V is variable. The paper says it uses "a minimal JSON-like grammar to reduce decoding entropy". Points are `[x, y]` and boxes are `[x1, y1, x2, y2]`. — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- Target construction: Suzuki–Abe contour tracing (OpenCV `findContours`), "enforce a consistent traversal order (clockwise), and optionally apply a tolerance-based sparsification" with tolerance ε. The paper does not name the simplifier. It is presumably Douglas–Peucker via `approxPolyDP`, but this is unverified. The ablation sweeps ε = 0.01 and 0.001, which suggests ε is in normalized-image units, i.e. 10 or 1 units on a 0–1000 scale (my inference). — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- Vertex density versus accuracy (Fig. 4, SFT only) is unimodal in token length:
  - too few points: 35.6 cIoU at 78 tokens;
  - too many: 72.5 cIoU at 859 tokens, attributed to "long-horizon decoding errors";
  - best: about 221 tokens.
  - Fig. 1 example caption: "around 70 points". — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- Ordering ablation (Fig. 7, qualitative only): "without enforcing clockwise ordering, the coordinate sequence fails to form a valid polygon"; "alternative orders [random, row-major] confuse the model and yield chaotic or repeated points". — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- Training:
  - Models: Qwen2.5-VL-7B and Kimi-VL (2.8B activated MoE).
  - SFT: 800k RefCOCO-family samples, 1 epoch, global batch 256 (about 3.1k steps), LR 5e-5 cosine to 2e-6, Muon optimizer.
  - RL: GSPO on 400k prompts, 2 epochs, 8 responses per group, rollout temperature 0.8, constant LR 2e-6.
  - KL coefficient is reported as 0.01 in the text but 0.1 in Table 5 (internal inconsistency). — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- Rewards:
  - mask IoU, thresholded at τ (stated range "[0.0, 0.1]", probably a typo);
  - negative centroid MSE distance;
  - format reward that zeroes the reward when the output does not parse.
  - Reward ablation (gIoU, RefC/RefC+/RefCg):
    - IoU only: 76.9 / 71.9 / 72.9
    - plus distance: 77.1 / 72.2 / 73.1
    - plus length penalty: 66.7 / 62.4 / 62.0. The paper says "imposing hard length constraints degrades the performance".
  - This is the source of the "caps hurt" claim. Strictly, the finding is about a length-penalty reward, not a parser-level vertex cap. — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- During RL "the model can adaptively adjust its output length". With dense targets (ε=0.001) the length drops. With sparse targets (ε=0.01) it rises slightly. RL "lets the model discover alternative yet valid trajectories (e.g., different starting points)". — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- Training stages (gIoU on RefCOCO val):
  - SFT only: 65.5
  - SFT + RL: 75.2 (+9.7)
  - web pretraining, no SFT: 25.3
  - pretraining + SFT: 70.1
  - pretraining + SFT + RL: 78.5. — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- Main results, RefCOCO val cIoU:
  - SimpleSeg-Qwen2.5-VL-7B: 75.6
  - SimpleSeg-Kimi-VL: 76.9
  - with pretraining: 80.9 and 80.0
  - UFO (LLaVA-1.5-7B): 77.2
  - Text4Seg (InternVL2-8B): 74.7. — [SimpleSeg](https://arxiv.org/abs/2601.19228). Columns were extracted from the PDF table; some testA values look out of order, so recheck before quoting.
- Multi-part: several polygons are written as separate `[[...]]` lists (Fig. 11). Holes are listed as a failure case: Fig. 12 shows a donut whose hole had to be prompted separately and then combined. — [SimpleSeg](https://arxiv.org/abs/2601.19228)

**Groundbench (arXiv 2609.27821, Aug 2026) [paper]** — [arXiv](https://arxiv.org/abs/2609.27821), [project](https://co-minder.github.io/Groundbench)
- Task: 1,500 RefCOCO-family referents, each re-targeted to exact-N polygons with N ∈ {8, 16, 24, 32, 64}. The prompt asks for "one line of 2N integers in {0, ..., 999}", which are mapped to pixels as x·W/999. — [Groundbench](https://arxiv.org/abs/2609.27821)
- Target construction (Algorithm 1):
  - To reach exactly N vertices from a sparse contour, repeatedly insert the midpoint of the longest edge.
  - From a dense contour, apply Visvalingam–Whyatt (repeatedly remove the vertex with the smallest neighbour-triangle area). Ties take the lowest index.
  - Start point and direction are unrestricted in scoring, because IoU is invariant to them.
  - Multi-part annotations keep only the part with the most vertices. This affects 148/1,500 referents; the kept part covers 71.76% of the union on average.
  - Construction fidelity (IoU of the N-vertex target against the source mask): 88.6 / 96.2 / 98.3 / 99.0 / 99.4 for N = 8 / 16 / 24 / 32 / 64. — [Groundbench](https://arxiv.org/abs/2609.27821)
- Legality rule: vertices must be distinct, in range and integer, and the polygon must be simple with positive area. Unparseable outputs score 0. — [Groundbench](https://arxiv.org/abs/2609.27821)
- Accuracy is non-monotone in N and collapses at 64. Fixed-denominator IoU for N = 8 / 16 / 24 / 32 / 64:
  - Qwen 3.7 Plus (max thinking): 48.2 / 49.7 / 49.1 / 48.6 / 38.5
  - Gemini 3.5 Flash: 62.6 / 64.7 / 65.9 / 57.8 / 37.7
  - At N=64 Gemini's legal-polygon coverage and legal IoU fall to 7.5 and 4.2.
  - The decline persists "when extra vertices only subdivide edges without changing the filled target". — [Groundbench](https://arxiv.org/abs/2609.27821)
- Box versus polygon: best box IoU 88.2 (Acc@.5 97.1) against direct polygons 57.7 (Acc@.5 69.2). A predicted box resampled to N points and scored against the same contour gets 57.3. In other words, zero-shot frontier polygons are no better than boxes. Qwen with thinking off falls to 20.5 / 14.7 / 17.0 IoU at N = 8 / 16 / 24. — [Groundbench](https://arxiv.org/abs/2609.27821)
- Groundbench evaluates hosted models zero-shot only; it does no fine-tuning. — [Groundbench](https://arxiv.org/abs/2609.27821)

**PolyFormer (CVPR 2023, arXiv 2302.07387) [paper]** — [arXiv](https://arxiv.org/abs/2302.07387)
- This is a specialist seq2seq model, not an LLM. It regresses floating-point coordinates (with bilinear-interpolated coordinate embeddings) instead of quantizing.
- Vertices run clockwise. The start is the vertex closest to the image's top-left corner.
- Multiple polygons are separated by `<SEP>` and ordered by the distance of their start points to the image origin. No holes. Sequences are wrapped in `<BOS>`/`<EOS>`. — [PolyFormer](https://arxiv.org/abs/2302.07387)
- Ablation (RefCOCO mIoU):
  - randomly ordered polygons: 55.92;
  - canonical ordering: +12.43;
  - "polygon augmentation" (densify the contour, then uniformly down-sample at a random interval to produce polygons of varying granularity): +3.72;
  - multitask with boxes: +3.07;
  - `<SEP>` token: +0.82. — [PolyFormer](https://arxiv.org/abs/2302.07387)
- Prior work SeqTR was limited to a single 18-vertex polygon. — [PolyFormer](https://arxiv.org/abs/2302.07387)

**VistaLLM (CVPR 2024, arXiv 2312.12423) [paper]** — [arXiv](https://arxiv.org/abs/2312.12423)
- 32 points per mask, quantized to 1000 bins. New tokens "10–999" were added because Vicuna's vocabulary only has the digits 0–9.
- Points are chosen by "gradient-aware adaptive sampling": keep the N points with the highest turning angle. This replaces uniform clockwise sampling. — [VistaLLM](https://arxiv.org/abs/2312.12423)
- Upper bound on RefCOCO val with 32 points: uniform sampling 94.70 mIoU, adaptive 97.26. Adaptive sampling improves trained results by 3–4 mIoU. Going from 16 to 32 points gives +3.6 (Ref) and +4.5 (Ref+). — [VistaLLM](https://arxiv.org/abs/2312.12423)
- Multiple objects: the prompt templates ask for "32 points along each object's contour", with separators between objects. No holes. — [VistaLLM](https://arxiv.org/abs/2312.12423)

**VectorLLM (arXiv 2507.04664, 2025) [paper]** — [arXiv](https://arxiv.org/abs/2507.04664)
- This is the most relevant small-model precedent. The LLM is **Qwen3-0.6B** or Qwen3-1.7B with a RADIO-L ViT and learnable position embeddings.
- Input is a single building cropped by enlarging the detector box 1.1–1.5×. Output is the building's corner points. — [VectorLLM](https://arxiv.org/abs/2507.04664)
- Coordinates are crop-pixel special tokens `[x0]..[x(W-1)]`, `[y0]..[y(H-1)]`. Example output: `[x85][y32][x160][y63]...[x85][y32]`. The sequence starts at the top-left-most corner, runs clockwise and repeats the start point at the end. — [VectorLLM](https://arxiv.org/abs/2507.04664)
- Start point: randomly shuffling the start point helped during *pretraining*. In *SFT*, a fixed start rule "brings slight performance improvements and stability". — [VectorLLM](https://arxiv.org/abs/2507.04664)
- Coordinate format ablation (Table 8, WHU): one special token per value scores 88.1 mAP. Splitting each value into hundreds/tens/ones special tokens scores 86.4 mAP and makes sequences longer. — [VectorLLM](https://arxiv.org/abs/2507.04664)
- Model size (WHU AP): Qwen3-0.6B 88.1, Qwen3-1.7B 90.2. Training used SFT followed by preference optimization (DPO-style chosen/rejected contours). — [VectorLLM](https://arxiv.org/abs/2507.04664)
- Building polygons with holes (courtyards) are not addressed. The model handles one instance per crop. — [VectorLLM](https://arxiv.org/abs/2507.04664)

**VectorLLM++ [secondary]** — [Research Square preprint](https://www.researchsquare.com/article/rs-10949232/v1)
- According to the search snippet: a "vector tokenizer that discretizes coordinates into a 1024-token vocabulary", a mask tokenizer with a codebook, 1.2M polygon annotations, 62.9M pretraining samples, 7.2M SFT instructions, and RL with verifiable geometric rewards.
- I did not read the full text. Ordering, simplification, holes and base model are unknown.

**FloorplanVLM (arXiv 2602.06507, 2026) [paper HTML via fetch]** — [arXiv](https://arxiv.org/abs/2602.06507)
- Base model: Qwen2.5-VL-3B.
- Coordinates are plain-text integers on a [0, 1024] grid after scaling the longer image edge to 1024.
- Output is hierarchical JSON: walls (start, end, thickness, curvature, openings), then rooms as closed cycles of wall indices. — [FloorplanVLM](https://arxiv.org/html/2602.06507)
- Format ablation:
  - JSON: 3,095 tokens, validity 90.20%, external IoU 0.8567
  - Python code: 2,038 tokens, validity 88.15%, external IoU 0.8555
  - Authors' conclusion: "pre-training alignment outweighs token efficiency".
  - 1,391 custom special tokens for JSON syntax and keys cut length by about 25%. — [FloorplanVLM](https://arxiv.org/html/2602.06507)
- Training:
  - SFT on Floorplan-2M for 2 epochs, then Floorplan-HQ-300K for 10 epochs ("quality annealing");
  - then GRPO with G=8 and β=0.01;
  - rewards: validity 0.1, external IoU 0.5, internal structure 0.4.
  - Final validity 96.10%, external IoU 0.9252. Holes are not addressed. — [FloorplanVLM](https://arxiv.org/html/2602.06507)

**Rex-Omni (arXiv 2510.12798, Oct 2025) [paper]** — [arXiv](https://arxiv.org/abs/2510.12798)
- Base model: Qwen2.5-VL-3B. The last 1,000 vocabulary entries are repurposed as coordinate tokens `<0>..<999>` (relative, quantized).
- A box costs 4 tokens versus 15 "atomic" tokens for digit text. Polygons (OCR) are written as `[x0, y0, x1, y1, ...]` using the same tokens. — [Rex-Omni](https://arxiv.org/abs/2510.12798)
- Data: 22M SFT samples, then GRPO on a 66K subset.
- SFT produced duplicate and degenerate predictions: 20.5% oversized "large boxes" on Dense200, falling to 3.5% after GRPO. The authors attribute this to teacher forcing never exposing the model to its own errors. This matches our zero-shot "repeated coordinate runs". — [Rex-Omni](https://arxiv.org/html/2510.12798v1)
- I found no clean ablation of special tokens versus digits in the parts I read.

**Florence-2 (CVPR 2024, arXiv 2311.06242) [paper + secondary]** — [arXiv](https://arxiv.org/abs/2311.06242)
- Adds location tokens for 1,000 quantized bins. Polygons are `(x0, y0, ..., xn, yn)` as location tokens, "in clockwise order". — [Florence-2](https://arxiv.org/abs/2311.06242)
- Output wraps each instance in `<ins>` and each part in `<poly>`, so one object can have several parts. — [Florence-2 docs (secondary)](https://docs.curiosity.ai/florence2-sharp/guides/phrase-grounding)
- One third-party page reports that `REFERRING_EXPRESSION_SEGMENTATION` "produces malformed location tokens" in practice. — [fal.ai (secondary, unverified)](https://fal.ai/models/fal-ai/florence-2-large/referring-expression-segmentation/llms.txt)

**PaliGemma / PaliGemma 2 (arXiv 2407.07726) [paper]** — [arXiv](https://arxiv.org/abs/2407.07726)
- Adds 1,024 `<loc0000>`–`<loc1023>` tokens and 128 `<seg000>`–`<seg127>` VQ-VAE codebook tokens.
- A referring mask is 4 loc tokens (box) plus 16 seg tokens describing a mask inside the box. A VQ-VAE decoder is needed to recover pixels. — [PaliGemma](https://arxiv.org/abs/2407.07726); [Google blog](https://developers.googleblog.com/en/gemma-explained-paligemma-architecture/)

**Text4Seg (ICLR 2025, arXiv 2410.09855) and GeoGround (arXiv 2411.11904) [paper]**
- Text4Seg writes "semantic descriptors": one text label per patch on a 16×16 grid (256 entries).
- Row-wise RLE (R-RLE) shortens this by 74% and speeds inference 3×. An optional SAM refiner is used.
- Holes and multiple parts are free because the output is a raster, but resolution is coarse. — [Text4Seg](https://arxiv.org/abs/2410.09855)
- GeoGround applies the same idea ("Text-Mask") to remote sensing at 32×32 with R-RLE. — [GeoGround](https://arxiv.org/abs/2411.11904)

**HiMTok (arXiv 2503.13026) [paper]**
- Represents a mask with up to 32 hierarchical, coarse-to-fine mask tokens. It needs a de-tokenizer and a 3-stage training recipe. — [HiMTok](https://arxiv.org/abs/2503.13026)

**UFO (arXiv 2503.01342) [paper]**
- Emits `<MASK>` tokens whose embeddings are dot-producted with image features ("embedding retrieval"). Predicting several mask tokens upsamples the mask: "any shape using only 16 tokens". — [UFO](https://arxiv.org/abs/2503.01342)

**LISA-style [SEG] [memory]**
- The LLM emits a `[SEG]` token whose hidden state prompts a SAM decoder, so no geometry appears in text. — [LISA](https://arxiv.org/abs/2308.00692) (not re-read this session)
- Recent RL variants such as Seg-R1 and Dr. Seg use GRPO to output boxes or points that prompt SAM2. Dr. Seg trained on about 10k samples. — [Dr. Seg](https://arxiv.org/abs/2603.00152), [Seg-R1](https://arxiv.org/html/2506.22624)

**Pix2Poly (arXiv 2412.07899) [paper]**
- Pix2Seq-style discretized building-corner tokens (vocabulary = N_B bins + start/end/pad). Polygon connectivity comes from a separate optimal-matching head that predicts a permutation matrix, so it is not pure text. — [Pix2Poly](https://arxiv.org/abs/2412.07899)

**Hi-Token (arXiv 2608.03471, 2026) [paper]** — [arXiv](https://arxiv.org/abs/2608.03471)
- Compares a "flat" vocabulary of 1,000 coordinate tokens with axis-specific hundreds/tens/ones digit tokens (60 token types, 12 tokens per box).
- Data: 80k RefCOCO samples on Qwen2.5-VL-3B.
- Digit tokens beat flat tokens by +14.1 mIoU, +14.9 P@0.5 and +8.7 P@0.95. The gain carries over to Qwen3-VL-2B (+9.0 P@0.95) and Qwen3-VL-4B (+2.6).
- Mechanism: with flat tokens, 985 of the 1,000 types appear fewer than 1,000 times in 80k samples. Digit tokens get about 50× denser supervision per type. — [Hi-Token](https://arxiv.org/abs/2608.03471)
- **Conflict:** VectorLLM found the opposite (single token beats multi-digit, 88.1 vs 86.4 mAP). VectorLLM had large building-pretraining data, while Hi-Token is the data-limited case.

**DenseMLLM (arXiv 2602.14134, 2026) [paper abstract]**
- Supervises vision-token logits for dense labels without a decoder. This alternative to polygon text is noted only for completeness. — [DenseMLLM](https://arxiv.org/abs/2602.14134)

### Inferences
- The polygon-as-text convention most papers converge on:
  - integer coordinates on about 1000 bins relative to the image or crop;
  - one closed ring per part;
  - clockwise traversal from the vertex nearest the top-left corner;
  - parts separated and ordered by start point;
  - an adaptive vertex count of roughly 16–40.
- Our draft (0–1000 crop-relative integers) fits that convention.
- Holes are an open problem in polygon-as-text VLMs. The most precedent-compatible way to handle them is our explicit `holes` list, plus a canonical orientation rule for each ring.
- For RIGA the rim's hole can be avoided entirely: emit the disc and cup as two simple polygons and compute rim = disc − cup in the checker.

### Gaps
- I did not inspect SimpleSeg's released code, so its exact simplifier and the units of ε are unverified.
- I could not obtain the VectorLLM++ full text.
- I did not find a "UFO" polygon variant or any explicit "Groundbench-style" fine-tuning study.
- I found no paper with a quantitative ablation of hole encodings in VLM text output.

---

## 2. Evidence on fixed vs adaptive vertex count, caps, ordering, relative vs absolute, special tokens, JSON vs compact, small-model SFT

### Takeaway
- Use an adaptive vertex count with a moderate token budget (SimpleSeg's optimum was about 221 tokens).
- Do not penalize length or cap it hard in the reward.
- Use a deterministic canonical order in SFT. This is the strongest single effect measured: +12.4 mIoU in PolyFormer, and SimpleSeg's model fails without it.
- Use relative quantized coordinates.
- Keep plain digit text rather than adding untrained `<loc>` tokens when data is small (Hi-Token).
- Use compact JSON: JSON was more valid than a terser code format (FloorplanVLM).

### Cited Findings
- **Fixed-N is fragile; adaptive is better:**
  - Groundbench's exact-N results collapse at N=64 even though the targets are more faithful. — [Groundbench](https://arxiv.org/abs/2609.27821)
  - VistaLLM's turning-angle sampling beats uniform sampling for a fixed N=32 (upper bound 97.26 vs 94.70). — [VistaLLM](https://arxiv.org/abs/2312.12423)
  - SimpleSeg uses tolerance-based (adaptive) sparsification, and during RL the length drifts to an "accuracy-efficiency balance". — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- **Caps and length penalties:** adding a length penalty to SimpleSeg's GSPO reward cost about 10 gIoU (76.9 → 66.7 on RefCOCO). — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- **Ordering:**
  - PolyFormer: random order 55.92 mIoU, canonical order +12.43. — [PolyFormer](https://arxiv.org/abs/2302.07387)
  - SimpleSeg: without clockwise order "no segmentation mask can be derived". — [SimpleSeg](https://arxiv.org/abs/2601.19228)
  - VectorLLM: a fixed start in SFT is slightly better and more stable; random start helps only in pretraining. — [VectorLLM](https://arxiv.org/abs/2507.04664)
- **Relative vs absolute coordinates:**
  - Rex-Omni chose relative 0–999 because it "reduces learning complexity" (bounded classes) and is shorter (4 vs 15 tokens per box). This was a design argument, not an ablation that I found. — [Rex-Omni](https://arxiv.org/abs/2510.12798)
  - Qwen3-VL itself moved from Qwen2.5-VL's absolute pixels to 0–1000 relative coordinates (see §3). — [Qwen3-VL issue #1616](https://github.com/QwenLM/Qwen3-VL/issues/1616)
- **Special tokens vs digits:** the evidence conflicts (Hi-Token favours digit tokens at 80k samples; VectorLLM favours single tokens). See §1. — [Hi-Token](https://arxiv.org/abs/2608.03471); [VectorLLM](https://arxiv.org/abs/2507.04664)
- **JSON vs compact:**
  - FloorplanVLM: JSON beat Python on validity (90.20% vs 88.15%) despite being 52% longer. — [FloorplanVLM](https://arxiv.org/html/2602.06507)
  - SimpleSeg uses a JSON-like nested list "to reduce decoding entropy". — [SimpleSeg](https://arxiv.org/abs/2601.19228)
  - Groundbench prompts for a flat line of integers, and frontier models still fail at large N. — [Groundbench](https://arxiv.org/abs/2609.27821)
- **Small models and SFT data scale:**
  - VectorLLM with Qwen3-0.6B reaches 88.1 AP on WHU single-building crops (1.7B: 90.2). — [VectorLLM](https://arxiv.org/abs/2507.04664)
  - Hi-Token SFT with 80k samples on Qwen3-VL-2B gives +9.0 P@0.95 over flat tokens. — [Hi-Token](https://arxiv.org/abs/2608.03471)
  - SimpleSeg uses 800k SFT samples for 1 epoch (about 3.1k steps at batch 256) on 7B models. — [SimpleSeg](https://arxiv.org/abs/2601.19228)
  - FloorplanVLM (3B) uses 2M + 10 epochs on 300K. — [FloorplanVLM](https://arxiv.org/html/2602.06507)
  - SFT alone leaves degenerate duplicates that GRPO removes (Rex-Omni: 20.5% → 3.5%). — [Rex-Omni](https://arxiv.org/html/2510.12798v1)
- **Our own vertex/IoU/token measurements on RIGA [local]:**
  - Setup: 40 random RIGA expert masks (`data/riga_true`, 127 = rim, 254 = cup). Each crop is a square of 1.5× the disc bbox (median 332 px), with coordinates on 0–1000 relative to the crop.
  - Ground-truth rings are the exact pixel-edge ("crack") rings from `objdel.checker.polygonize`.
  - Simplification is either Douglas–Peucker (`cv2.approxPolyDP`, ε in 0–1000 units) or uniform arc-length resampling to N points. After rounding to integers, polygons are rasterized with an exact pixel-centre even-odd test.
  - Tokens are Qwen3-VL-2B counts of the draft compact JSON. Medians:

| object | simplifier | vertices | IoU median | IoU p5 | tokens |
|---|---|---|---|---|---|
| disc | DP ε=0.5–1 | 488 | 1.0000 | 0.9996 | 3,915 |
| disc | DP ε=3 | 87 | 0.9961 | 0.9951 | 707 |
| disc | DP ε=5 | 26 | 0.9929 | 0.9915 | 223 |
| disc | DP ε=8 | 18 | 0.9889 | 0.9862 | 155 |
| disc | arc N=16 / 24 / 32 | 16 / 24 / 32 | 0.9744 / 0.9877 / 0.9918 | 0.9712 / 0.9845 / 0.9896 | 139 / 203 / 267 |
| cup | DP ε=3 | 44 | 0.9928 | 0.9890 | 359 |
| cup | DP ε=5 | 18 | 0.9873 | 0.9826 | 155 |
| cup | arc N=16 / 24 / 32 | 16 / 24 / 32 | 0.9739 / 0.9859 / 0.9892 | 0.9679 / 0.9808 / 0.9862 | 139 / 203 / 267 |
| rim (outer+hole) | DP ε=5 | 44 | 0.9862 | 0.9772 | 364 |
| rim (outer+hole) | DP ε=8 | 32 | 0.9787 | 0.9657 | 272 |
| rim (outer+hole) | arc N=24+24 | 48 | 0.9789 | 0.9678 | 396 |

  - At equal vertex count, DP beats arc-length resampling. For example, disc DP ε=5 gives 26 vertices at 0.9929, against arc N=24 at 0.9877.
  - The 0–1000 quantization of a crop of about 332 px is lossless (DP ε=0.5 gives IoU 1.0000).
  - One caveat: with the contour taken through *pixel centres* (`cv2.findContours` + 0.5) instead of crack edges, the IoU ceiling drops to about 0.99 for the disc and 0.98 for the cup, and small DP ε then appears *worse* than ε=3. That is a measurement artifact of the contour convention.

### Inferences
- RIGA's expert ambiguity band is large: median ambiguous band 17% for the disc and 55% for the cup (from `experiments/phase0/RESULTS_0b.md`). DP ε≈5 targets lose only about 1% IoU. Simplification error is therefore an order of magnitude below annotator disagreement, so a tolerance of about 5 units costs nothing that matters for the reward.
- Converting SimpleSeg's best SFT budget (about 221 tokens) at our measured 8 tokens per vertex gives about 27 vertices per object. That matches DP ε=5 (disc 26, cup 18). The SimpleSeg budget used decimals, so treat this as an approximate cross-check.
- Use each of the 6 RIGA experts (and the 2–4 LIDC readers) as a separate SFT target for the same image. This is a natural form of PolyFormer-style target diversity. It may also teach the model to hedge toward a typical expert.

### Gaps
- I found no paper that runs polygon SFT on models under 1B *VLM* parameters (SmolVLM class). VectorLLM-0.6B is the closest, and its vision encoder is RADIO-L, not SigLIP-base.
- No paper reports a steps-to-validity curve. Our own measurement (hypothesis H3 below) is required.

---

## 3. Native coordinate conventions of Qwen2.5-VL / Qwen3-VL / SmolVLM, and tokenizer cost

### Takeaway
- Qwen3-VL natively uses 0–1000 relative integer coordinates. Qwen2.5-VL used absolute pixels of the resized input.
- SmolVLM has no documented grounding or coordinate convention.
- Both tokenizers split numbers into single digits. A compact `[[x,y],...]` vertex costs about 8 tokens; default `json.dumps` spacing costs 10; SimpleSeg-style decimals cost 13.6.

### Cited Findings
- Qwen maintainers recommend absolute coordinates for Qwen2.5-VL and 0–1000 normalized coordinates for Qwen3-VL. The formula given is `round(bbox / [w,h,w,h] * 1000)`. — [QwenLM/Qwen3-VL #1616](https://github.com/QwenLM/Qwen3-VL/issues/1616)
- Some third-party runners describe the range as 0–999 and convert with x·w/999. Groundbench also maps with /999. This is an off-by-one ambiguity of about 0.1% scale. — [Groundbench](https://arxiv.org/abs/2609.27821); [search summary (secondary)](https://lmstudio.ai/ceveyne/analyse-image/files/src/fastvlm_server/qwen3_vl_runner.py)
- Rex-Omni describes Qwen2.5-VL as absolute coordinates whose digits are tokenized one by one ("1921" becomes 1, 9, 2, 1). — [Rex-Omni](https://arxiv.org/abs/2510.12798)
- The SmolVLM technical report (arXiv 2504.05299) mentions grounding only when describing Qwen-VL in related work. I found no grounding or coordinate data in its training mix. The training emphasis is documents (25%) and captioning (18%). — [SmolVLM](https://arxiv.org/abs/2504.05299); [search summary](https://huggingface.co/HuggingFaceTB/SmolVLM-Instruct/blame/77f80ed86046795eab54cfb60d9442c2ad813f57/README.md)
- **Tokenizer measurements [local]** (Qwen/Qwen3-VL-2B-Instruct and HuggingFaceTB/SmolVLM-256M-Instruct / SmolVLM2-500M-Video-Instruct tokenizers):
  - Every digit is a separate token.
    - `[[512,487],[1000,7]]` becomes `[[ 5 1 2 , 4 8 7 ],[ 1 0 0 0 , 7 ]]`.
    - A space before a number is its own token (`Ġ`).
    - `[[0.512, 0.487]]` becomes `[[ 0 . 5 1 2 , Ġ 0 . 4 8 7 ]]`.
  - Tokens per vertex, 32-vertex ring (identical for Qwen3-VL and SmolVLM-256M):
    - compact JSON `[[x,y],...]` 0–1000: **8.03**
    - flat `x,y,x,y`: 7.97
    - default `json.dumps` with spaces: 10.0
    - absolute pixels on a 2376-px image: 8.97
    - SimpleSeg decimals `[[0.xxx, 0.yyy], ...]`: 13.59
  - Full draft schema `{"rim":[{"outer":[...],"holes":[[...]]}]}`:
    - compact: 8.1–8.4 tokens per vertex
    - `json.dumps` default: 10.1–10.4
    - pretty-printed (indent=1): 13–17
  - The fixed key overhead is only about 10 tokens.

### Inferences
- With Qwen3-VL-2B, 0–1000 crop-relative integers match the pretraining prior. Expect better zero-shot behaviour and faster SFT than with any new token scheme.
- For SmolVLM there is no coordinate prior to preserve, so the choice should be driven by data efficiency. Plain digits reuse 10 well-trained token types, which Hi-Token's density argument favours at our data scale (thousands, not millions, of masks). Adding 1,000 `<loc>` tokens would leave most of them under-trained.
- Use the same format for both models so that comparisons between them are clean.
- Emit compact JSON with no spaces (about 20% fewer tokens than default `json.dumps`), and train and parse with exactly that serialization. One parser-side caution: the current `extract_json` in `src/objdel/checker/polygon_io.py` accepts any JSON, so whitespace variants still parse (correction, 2026-10-19: `polygon_io.py` was deleted in Phase 1 (M2); the production path is `objdel.geometry` (strict `parse`, polycore `to_mask`, pixel-centre even-odd); whitespace inside the text is a parse error there).

### Gaps
- I did not find primary documentation (model card or cookbook) for Qwen3-VL *polygon* output. The issue above covers boxes only.
- I could not confirm whether Qwen3-VL was pretrained on any polygon or contour data.

---

## 4. Mask-from-polygon semantics and pitfalls (fill rule, pixel-centre convention, self-intersection, winding)

### Takeaway
Most papers do not state their rasterizer. Groundbench is the only one that freezes a policy ("clipped cyclic fill", with legality checked separately). Our measurements show the convention matters a lot for small LIDC nodules. Mixing contour and rasterizer conventions costs 3–35% IoU on objects of radius 3–12 px. **Our repo's current `rasterize` (PIL `polygon(fill=1, outline=1)`) inflates crack-edge polygons by 12–56% in area for nodules of radius 3–12 px.**

### Cited Findings
- Groundbench calls rasterization a set of "explicit choices about disconnected parts, self-intersections, clipping, and discrete integration". Its choices:
  - freeze one policy;
  - score filled IoU and legality separately (legality: distinct vertices, simple, positive area);
  - check robustness against a continuous-area integrator, which moved means by ≤1.1 points on RefCOCO-size objects. — [Groundbench](https://arxiv.org/abs/2609.27821)
- SimpleSeg traces contours with Suzuki–Abe / OpenCV. It does not state how it rasterizes back to a mask for the IoU reward. — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- Common toolchains differ: COCO uses `pycocotools frPyObjects`; supervision uses `cv2.fillPoly`. — [supervision docs](https://supervision.roboflow.com/0.24.0/detection/utils/)
- **Measurement [local]** (synthetic disks and annuli with jittered centres, 20 per size; crop = 4r+4 px; round trip through 0–1000 integers; table shows mean IoU against the source mask, with area ratio in brackets):
  - Contour = exact crack-edge rings (`objdel.checker.polygonize`), rasterized with a pixel-centre even-odd test: **IoU 1.000, area 1.00, for every size r=3–20, disks and annuli.** The round trip is exact.
  - Crack-edge rings rasterized with PIL fill+outline (the repo's then-current `polygon_io.rasterize`; correction, 2026-10-19: `polygon_io.py` was deleted in Phase 1 (M2); the production path is `objdel.geometry` (strict `parse`, polycore `to_mask`, pixel-centre even-odd)):

| radius r (px) | 3 | 5 | 8 | 12 | 20 |
|---|---|---|---|---|---|
| IoU | 0.643 | 0.789 | 0.897 | 0.878 | 0.969 |
| area ratio | 1.56 | 1.27 | 1.12 | 1.14 | 1.03 |

  - Pixel-centre contours (`cv2.findContours`) rasterized with a pixel-centre fill: area 0.71 / 0.89 / 0.90 / 0.93 / 0.93 for r = 3 / 5 / 8 / 12 / 20. Annuli are worse (0.42 at r=3), because the hole contour hugs the foreground.
  - My `cv2.fillPoly` runs were confounded by a half-pixel offset in my own script (OpenCV puts pixel centres at integer coordinates). Treat those numbers as a demonstration of how fragile conventions are, not as a benchmark.
- Winding: the repo's `polygonize` emits clockwise-on-screen outer rings (positive shoelace area, y down) and negative-area holes, starting from a canonical vertex (minimum (y, x) edge key). — `src/objdel/checker/polygonize.py` [local]

### Inferences
- Use one convention end to end. Vertices sit on pixel *corners*, in 0–1000 units relative to the crop, with px = v/1000·S. Fill by testing whether each pixel *centre* (i+0.5, j+0.5) is inside, using even-odd across the outer ring and its holes. Then the checker's IoU of a perfect target is exactly 1.0 at every nodule size.
- With the current PIL outline=1 fill, the reward systematically favours polygons drawn about 0.5 px *inside* the true boundary. That bias is largest exactly where LIDC nodules live (r ≈ 3–10 px).
- Even-odd versus nonzero fill only differs for self-intersecting or overlapping rings. Reject those through legality (`ring_problems` already flags `self_crossing` and `self_touching`) rather than relying on a fill rule.
- Additional legality checks are needed:
  - holes lie inside their outer ring;
  - holes are disjoint from each other;
  - parts are disjoint from each other.
- Douglas–Peucker can create self-intersections, especially on thin rims and at a hole near its outer ring. Target construction must validate every simplified ring. When a ring fails, retry with a smaller ε, or use topology-preserving simplification.

### Gaps
- No paper reports how many of its SFT targets were invalid after simplification.
- I did not measure DP self-intersection rates on LIDC.

---

## 5. Data augmentation and small-model SFT recipe evidence

### Takeaway
- Published polygon-VLM recipes augment the *targets* rather than the images:
  - PolyFormer uses random vertex granularity (+3.72 mIoU);
  - VectorLLM uses random start points (pretraining only) and crop jitter (1.1–1.5× box).
- I found no paper that reports geometric flips or rotations with polygon transforms for VLM polygon SFT.

### Cited Findings
- PolyFormer polygon augmentation: densify, then down-sample at a random interval. Gain +3.72 mIoU. — [PolyFormer](https://arxiv.org/abs/2302.07387)
- VectorLLM:
  - random scaling of the crop box by 1.1–1.5×;
  - random start point helps pretraining;
  - fixed start is better for SFT. — [VectorLLM](https://arxiv.org/abs/2507.04664)
- SimpleSeg adds task multiplexing over the tuple [text, point, box, mask]: a box or point prompt maps to a mask, and weak labels are recombined, "multiplying supervision sources". — [SimpleSeg](https://arxiv.org/abs/2601.19228)
- FloorplanVLM does not mention flips or rotations. — [FloorplanVLM](https://arxiv.org/html/2602.06507)

### Inferences
- Horizontal or vertical flips and 90° rotations are exact on the pixel grid and preserve topology. They reverse winding (flips) or move the canonical start (both). After transforming, re-canonicalize: re-orient to clockwise and re-pick the top-most-then-left-most start. Otherwise the augmentation teaches inconsistent order, which is the very thing PolyFormer and SimpleSeg show is harmful.
- Fundus images have a left/right-eye asymmetry, so horizontal flips are anatomically plausible for RIGA. Check whether 90° rotations hurt, because the optic disc's vertical elongation is diagnostic (vertical cup-to-disc ratio).
- Varying ε per sample (for example ε ∈ [3, 8]) is the cheap PolyFormer-style analogue. It may also reduce over-commitment to exact vertex positions before GRPO.

### Gaps
- There are no published numbers on SFT sample counts or steps for VLMs under 1B.

---

## Candidate hypotheses to test locally

Each test is cheap and measurable, and most run on the M3 without training.

1. **H1 – Token budget per real object (no training).** Use the RIGA DP-ε table above. Repeat it on LIDC nodule crops: crack-edge rings, ε ∈ {2, 3, 5, 8} in 0–1000 units, with crop sizes as planned for Phase 1. Measure vertices, tokens (Qwen3-VL-2B and SmolVLM tokenizers), and the IoU loss compared with the inter-reader band. Expected: under 150 tokens per nodule at ε=5.
2. **H2 – Rasterizer bias (no training).** Replace `polygon_io.rasterize` (PIL fill+outline) with a pixel-centre even-odd fill on corner-convention vertices (done in Phase 1: plan D1, `objdel.checker.polycore`; `polygon_io.py` deleted). Re-score the Phase-0 MedSAM and VLM outputs and the expert-vs-expert IoUs. Measure how much IoU and band membership change for LIDC slices with nodule radius ≤ 8 px. Expected: shifts of several IoU points on small nodules.
3. **H3 – Steps to validity for small models.** SFT SmolVLM-256M, SmolVLM-500M and Qwen3-VL-2B (LoRA vs full) on RIGA disc+cup targets at ε=5. Every N steps, record:
   - parse rate;
   - legal-polygon rate (`ring_problems`);
   - duplicate-run rate (the zero-shot failure mode);
   - IoU against the nearest expert.
   This answers "how many examples or steps" directly.
4. **H4 – Rim as a hole vs disc+cup as two simple rings.** Train two otherwise identical SFTs:
   - (a) `{"rim":[{"outer":...,"holes":[...]}]}`
   - (b) `{"disc":[...],"cup":[...]}` with rim derived.
   Compare rim IoU, topology accuracy (1 part, 1 hole) and legality.
5. **H5 – Adaptive DP vs fixed-N arc-length targets.** Compare ε=5 against N=24 at matched mean tokens. Measure post-SFT IoU and validity. The prior (VistaLLM, Groundbench) predicts adaptive is better.
6. **H6 – Canonical order vs random start in SFT.** Use the fixed top-left-most start, clockwise, against random start plus random direction. The prior (PolyFormer +12.4, VectorLLM) predicts canonical wins for SFT.
7. **H7 – Digits vs added location tokens (SmolVLM only).** Compare plain 0–1000 digits with 1,001 new `<loc>` tokens at our data scale. Hi-Token predicts digits win when data is limited.
8. **H8 – Compact vs spaced JSON vs flat list.** Measure validity after SFT at equal steps. Expected to be a small effect; FloorplanVLM suggests keeping JSON.
9. **H9 – Expert-as-target diversity.** Train on 1 random expert per image per epoch versus a fixed single expert. Measure band-membership rate and the GRPO starting reward.
10. **H10 – Flip/rot90 augmentation with re-canonicalization**, against no augmentation and against flips *without* re-canonicalization (a negative control).
11. **H11 – DP legality.** Measure the fraction of simplified LIDC and RIGA rings that self-intersect, or whose hole touches its outer ring, at each ε, and the cost of the fallback.

---

## Recommendation

1. **Coordinates.** Use integers 0–1000 relative to the square crop, with x before y and vertices on pixel corners (px = v·S/1000). This is Qwen3-VL's native convention and matches Florence-2, Rex-Omni and Groundbench-style quantization. Use plain digit text with no new tokens, for both Qwen3-VL-2B and SmolVLM. Document the /1000 divisor. Groundbench and some Qwen runners use /999, which is a 0.1% difference.
2. **Serialization.** Keep JSON (evidence: FloorplanVLM validity; SimpleSeg's JSON-like grammar), but serialize it compactly: `separators=(",",":")`, no spaces, no closing-vertex repeat. That costs about 8 tokens per vertex. Recommended schemas:
   - **RIGA:** `{"disc":[[x,y],...],"cup":[[x,y],...]}`. Two simple rings, with rim = disc − cup computed by the checker. The rim's hole then comes for free, and the topology check reduces to "cup inside disc".
   - **LIDC:** `{"nodule":[{"outer":[[x,y],...],"holes":[[[x,y],...]]}]}` as drafted. Always emit `"holes":[]` when there are none, so the grammar is fixed.
3. **Canonical order.**
   - Outer rings run clockwise on screen (y down); holes run counter-clockwise. This matches `polygonize` signed areas.
   - Each ring starts at its top-most, then left-most vertex.
   - Parts are sorted by that start vertex (row-major); holes within a part likewise.
   - Re-canonicalize after every augmentation.
   - The checker's reward should stay invariant to start and direction (IoU already is), so GRPO can explore, as SimpleSeg notes.
4. **Target construction.**
   - Start from exact crack-edge rings (`objdel.checker.polygonize`), not `cv2.findContours` pixel centres.
   - Simplify with Douglas–Peucker at ε ≈ 5 on the 0–1000 crop scale. Expect about 26 disc, 18 cup and roughly 150–220 tokens per object, at IoU ≥ 0.987 against the source mask.
   - Validate every ring (simple; holes inside and disjoint). If validation fails, retry with smaller ε.
   - For SFT diversity, sample ε uniformly from [3, 8] and the expert uniformly from the 6 RIGA experts or 2–4 LIDC readers.
   - Do not use fixed-N targets.
5. **No hard caps in the reward.** Only a generous `max_new_tokens` safety limit, about 3× the 95th-percentile target length, and truncation counts as a format failure. SimpleSeg's length penalty cost about 10 gIoU.
6. **Fix the rasterizer before any RL.** Use a pixel-centre even-odd fill on corner-convention vertices. The current PIL `outline=1` fill biases IoU and band rewards on small nodules.
7. **Augmentation.**
   - Crop jitter of 1.1–1.5× the object box (VectorLLM).
   - Horizontal flips, with re-canonicalized order.
   - ε jitter.
   - Hold off on 90° rotations for RIGA until H10 shows they don't hurt.
8. **Expect SFT to leave degenerate repeats that GRPO must clean up.** Rex-Omni went from 20.5% to 3.5% degenerate outputs after GRPO, and SimpleSeg's SFT-only result was +9.7 gIoU below its post-RL result. Track duplicate-run rate as a first-class SFT metric (H3).
