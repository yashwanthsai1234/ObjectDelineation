# Geometry-in-the-loop ML (round 2): an exact GPU polygonizer and bit-quad topology counts inside RL fine-tuning, training, data generation and evaluation (2025–2026)

Scope and method notes:
- "The engine" means the user's CUDA operator. It takes a binary bitmap, classifies 2x2 bit-quads, and emits closed rings with holes. The same pass gives χ, area, perimeter and vertex count.
- This round builds on round 1 (`ml_geometry_representation.md`, `algorithms_ml.md`, `ai_infra_data.md`, `semiconductor_eda.md`) and does not repeat it.
- **Coverage caveat:** the session-wide web-search budget ran out about a third of the way through this task. After that, only direct fetches of known arXiv or GitHub URLs were possible. As a result, these were **not re-verified this round**: UniVG-R1, Visual-RFT, VLM-R1 and Perception-R1 reward details (only the Perception-R1 abstract was fetched); floorplan or CAD generators with validity rewards; PolyFormer's multi-polygon handling; and 3D-printing/CNC manufacturability-aware ML. They are listed under Gaps.

---

## Q1. RL for segmentation and grounding (2025–2026): what rewards are used, how are they computed, and is reward computation a reported bottleneck?

### Takeaway
Nearly all 2025–2026 RL segmentation work follows one pattern: the VLM emits boxes and points, a frozen SAM/SAM2 makes the mask, and the reward is box IoU, point L1 or mask IoU, plus format terms. Only SimpleSeg (polygon points emitted as text, IoU after rasterizing) and ALTo (mask tokens with a length penalty) reward the geometry the model emits directly. **None of the VLM-segmentation papers I fetched uses topology, validity, hole or vertex-budget rewards, and none reports reward computation as a bottleneck.** The bottleneck case appears in **lithography RL**. LithoGRPO (ICML 2026) found that the standard shot-count metric took about 34–148 s per mask, which made GRPO sampling impractical. The authors had to write a fast CPU approximation that runs in 0.16–0.30 s per mask.

### Cited Findings
**VLM segmentation and grounding with RL**
- **Seg-Zero** (2025): a decoupled MLLM plus SAM, trained only with GRPO. The rewards are format rewards plus accuracy rewards "based on IoU of segmentation and L1 distance metrics of the prompt". It scores 57.5 zero-shot on ReasonSeg, 18% above LISA-7B, after training on 9,000 RefCOCOg samples. — [Seg-Zero (alphaxiv 2503.06520)](https://www.alphaxiv.org/abs/2503.06520.md); [Reasoning Segmentation survey arXiv 2505.18816](https://arxiv.org/pdf/2505.18816)
- **SAM-R1** (NeurIPS 2025): SAM runs "directly into the reinforcement learning training loop as a reward signal generator".
  - The reward is a tiered mask IoU: 4 if IoU>0.80, 3 if 0.70–0.80, 2 if 0.50–0.70, otherwise 0. Reasoning-format and segmentation-format rewards are added.
  - Setup: 8 responses per question, 3,000 training samples, 8×A100.
  - The paper gives no reward latency figures. — [arXiv 2505.22596](https://arxiv.org/html/2505.22596); [NeurIPS 2025](https://proceedings.neurips.cc/paper_files/paper/2025/hash/c9c54ac0dd5e942b99b2b51c297544fd-Abstract-Conference.html)
- **VisionReasoner** (ICLR 2026) states explicitly: "We use bounding boxes and points to calculate the rewards rather than binary masks."
  - Rewards: box IoU>0.5, box L1<10 px and point L1<30 px, each adding 1/max{N,K}.
  - Multi-object Hungarian matching with batching takes 5×10⁻⁴ s for 30 objects, against 2×10⁻³ s without batching, a 4× speedup.
  - SAM2 produces masks only at inference. — [arXiv 2505.12081v5](https://arxiv.org/html/2505.12081v5); [ICLR 2026 poster](https://www.iclr.cc/virtual/2026/poster/10009554)
- **LENS** (2025): R = λ₁·format + λ₂·box IoU + λ₃·segment IoU, where segment IoU is the mask IoU of SAM2's output against ground truth.
  - Setup: GRPO with G=8, batch 64, 16 RL epochs on 16×L40S.
  - Inference latency is dominated by the MLLM at about 1.5 s; the connector takes 0.002 s. — [arXiv 2508.14153](https://arxiv.org/html/2508.14153)
- **Seg-R1** (2025): GRPO on point and box prompts to a frozen SAM2, trained only on COD/SOD foreground image-mask pairs. It reaches 0.873 S-measure on COD10K. — [arXiv 2506.22624](https://arxiv.org/abs/2506.22624)
- **GenSeg-R1** (2026): Qwen3-VL emits boxes plus interior keypoints for a frozen SAM 2. One variant uses "a SAM 2 in-the-loop reward that directly optimizes mask quality". It reports 0.7127 cIoU on RefCOCOg val. — [arXiv 2602.09701](https://arxiv.org/abs/2602.09701)
- Other GRPO variants found but not fetched:
  - GRPO-SEG: box and point prompts for SAM2, with re-prompting of hard samples. — [Springer TVC 2026](https://link.springer.com/article/10.1007/s00371-026-04679-5)
  - GRASP: geospatial pixel reasoning. — [arXiv 2508.17102](https://arxiv.org/html/2508.17102)
  - FineRS: small objects. — [arXiv 2510.21311](https://arxiv.org/pdf/2510.21311)
- **Perception-R1**: rule-based GRPO on Qwen2.5-VL-3B. It reports +4.2% on RefCOCO+, +17.9% on PixMo-Count and 31.9 AP on COCO2017 val, and finds that "perceptual complexity is a major factor in determining the effectiveness of RL". Its reward details were not in the fetched abstract. — [arXiv 2504.07954](https://arxiv.org/abs/2504.07954)

**Policies that emit the geometry directly**
- **SimpleSeg** (arXiv 2601.19228, Kimi-VL team) emits polygon points as text.
  - RL uses **GSPO** with three rewards: (a) "Mask IoU reward: the direct IoU between the predicted and ground-truth mask", obtained by rasterizing the predicted contour; (b) a negative centroid MSE distance; (c) a format reward that returns zero on malformed output.
  - Setup: 8 responses per group, global batch 256, 400k RL prompts per epoch, 2 epochs, 32 GPUs.
  - Imposing hard length constraints "degrades the performance".
  - The paper does not describe the cost of reward computation or how rasterization is implemented. — [alphaxiv summary of 2601.19228](https://www.alphaxiv.org/abs/2601.19228.md)
- **ALTo** (NeurIPS 2025) uses adaptive-length mask tokens.
  - The GRPO reward is R = 𝕀_format + IoU − α·L, with mask IoU taken after the mask de-tokenizer.
  - Tokens per mask range from 1 to 32, averaging 17.5 on gRefCOCO val against a fixed 32.
  - The RL stage took 5 h on 8×A100. — [arXiv 2505.16495](https://arxiv.org/html/2505.16495v1)

**Render-based rewards for vector output (closest analogue for polygon rewards)**
- **RLRF** (SVG generation) renders every rollout with **CairoSVG**.
  - Rewards: L2, Canny-edge L2, DreamSim, CLIP, and a quadratic penalty on SVG length.
  - Each step has 32 images × 64 rollouts = **2,048 renders**. Im2SVG RLRF took about 3 days on 4 nodes × 8×H100 for 500 steps.
  - The fetched portion does not describe how invalid SVGs are handled. — [arXiv 2505.20793v2](https://arxiv.org/html/2505.20793v2)

**Lithography RL: geometric reward computation is a documented bottleneck**
- **LithoGRPO** (Lai, …, Bei Yu, Ping Luo; ICML 2026; arXiv 2606.00228) combines flow matching with GRPO fine-tuning for ILT.
  - The reward is R = −Σₖ k(x)/k₀ over four metrics: L2, PVB, EPE and **shot count**, where shot count is the minimum number of rectangles for e-beam writing.
  - The traditional shot count is the LithoBench implementation of NP-hard non-overlapping rectangular partitioning. It took **34.01 s (MetalSet), 64.82 s (StdMetal), 44.81 s (ViaSet) and 148.35 s (StdContact) per mask**.
  - The new "fast shot-counting" method works in four steps:
    1. Enumerate maximal rectangles by histogram sweep, O(N²) at N=512.
    2. Prune contained rectangles in Numba.
    3. Formulate set cover as an ILP.
    4. Solve it with PuLP.
  - It allows overlapping shots, which the authors say "aligns with modern multi-beam mask writing practice". It runs in **0.16–0.30 s**, which is **134.6× to 491.3×** faster, with R²=0.994 against the traditional count.
  - The authors note that group-wise GRPO normalization "cancels any constant offset", so only within-group ranking matters.
  - It does **not** extract contours or polygons.
  - Hardware: 4× RTX 3090, each stage under 8 h, group size G=6, 1,000 RL steps. — [arXiv 2606.00228](https://arxiv.org/html/2606.00228); [ICML 2026 poster](https://icml.cc/virtual/2026/poster/61189)
- **Generative RL ILT** (Haoyu Yang and Haoxing Ren, NVIDIA; arXiv 2602.19027) uses a WGAN generator fine-tuned with GRPO. It samples mask candidates, runs batched ILT refinement and selects on litho metrics such as EPE and process window. It reports a 3× speedup over SOTA numerical ILT and 20% better EPE on ICCAD13. The abstract does not say whether manufacturability (shot count or MRC) enters the reward. — [arXiv 2602.19027](https://arxiv.org/abs/2602.19027)

### Inferences
- **Where the engine plugs into RL rewards.** There are two distinct cases.
  - *Polygon-emitting policies* (SimpleSeg, Groundbench-style models). The reward needs polygon→mask rasterization, which is the reverse of what the engine does, and that step is cheap. The engine adds value through **canonicalization**: rasterize the emitted ring with a fill rule, then re-polygonize. That yields exact validity, β₀/β₁, hole correctness and vertex-count terms, none of which appear in any reward found.
  - *Mask-emitting policies* (SAM-in-loop, ALTo, flow or diffusion mask generators such as LithoGRPO). Any polygon-level reward requires mask→polygon on every rollout. That is exactly the engine's direction, and it is where a CPU cost was found to be prohibitive (LithoGRPO).
- **A reward from bit-quad statistics.** One GPU bit-quad pass gives the rectilinear vertex count n = n1 + n3 + 2nD, plus hole count and χ. LithoGRPO shows that a cheap proxy only needs to preserve ranking within a group. So a microsecond-scale, rank-faithful proxy built from bit-quad counts is plausible.
  - Example: a classical formula for the minimum rectangle partition of a rectilinear polygon is roughly n/2 + h − g − 1 (h holes, g non-intersecting chords). Bit-quads give n and h exactly; g needs the rings. I recalled this formula from memory and did not fetch a source, so verify it before relying on it.
  - The engine would replace a 0.16–0.30 s CPU ILP per 512² mask in the inner loop, about 0.2 s × G=6 per target.
- **Scale estimate.** For VLM segmentation RL at SimpleSeg scale (about 400k prompts × 8 rollouts × 2 epochs ≈ 6.4M reward evaluations), a CPU polygon rasterization plus IoU at about 1 ms each comes to about 1.8 CPU-hours. That is negligible next to VLM decoding (about 1.5 s per sample for LENS). So in VLM RL the engine's value is **new reward types** (topology, validity, holes, vertex budget, per-instance β), not speed. This estimate is my own arithmetic, not a measurement.
- **Dense multi-instance rewards** (SAM 3-style "all instances of a concept", about 100 objects per image) would multiply per-instance polygon or topology checks. That case could make CPU polygonization noticeable, but no paper measures it.

### Gaps
- Reward details for UniVG-R1, Visual-RFT, VLM-R1 and Perception-R1 were not fetched because the search budget was exhausted. From the general literature they are believed to use box IoU and similar rewards, but this was not verified this round.
- No paper reports wall-clock reward computation time for mask-IoU rewards in VLM RL.
- I found no RL paper using a polygon validity, topology or hole reward for VLM segmentation. This absence was established with a truncated search, so it is moderately confident, not certain.

---

## Q2. Straight-through or differentiable approximations for non-differentiable geometric post-processing, and topology losses with runtimes

### Takeaway
I found **no published "exact-forward, soft-backward" polygonizer** (a straight-through estimator that runs exact contour extraction forward and uses surrogate gradients backward). The field takes four other routes:
1. Differentiable *polygon→raster* (BoundaryFormer, DDSL, soft winding number), where the network predicts vertices.
2. Soft surrogates for bit-quad and Minkowski statistics (Fast-χ, Minkowski Image Loss, ECC CUDA kernels), which produce counts but no rings.
3. Learning auxiliary fields so that a non-differentiable polygonizer behaves well (Frame Field Learning).
4. **RL around the non-differentiable geometry** (LithoGRPO's shot count).

Topology losses cost tens to hundreds of milliseconds per loss call, and Betti matching runs on CPU in C++. That leaves room for an exact GPU β₀/β₁ pass used as gating or weighting.

### Cited Findings
**Topology losses and their cost**
- **Topograph** (ICLR 2025) reports average loss-computation time and training time in its Table 2:

  | Loss | Avg. loss calc. | Train time |
  |------|-----------------|------------|
  | clDice | 9.88 ms | 1h27m |
  | Mosin | 233.00 ms | 2h29m |
  | Hu et al. (PH) | 602.99 ms | 3h05m |
  | Betti matching | 327.64 ms | 2h46m |
  | Topograph | 95.94 ms | 1h53m |

  It claims a "3-6 fold improvement compared to the PH methods" and saves about 1 h per run on Platelet. The fetched text does not state the hardware, image size, or whether the computation is on CPU or GPU. — [arXiv 2411.03228v2](https://arxiv.org/html/2411.03228v2); [ICLR 2025](https://iclr.cc/virtual/2025/poster/29716)
- **Efficient Betti Matching 3D** (Stucki, Bürgin, Paetzold, Bauer, 2024) is a "highly optimized implementation of Betti matching, implemented in C++ together with a python interface", with "significant speedups" over Cubical Ripser. That makes it a CPU implementation. — [arXiv 2407.04683](https://arxiv.org/abs/2407.04683)
- **Skeleton Recall Loss** (ECCV 2024) gets its speed by "circumventing intensive GPU-based calculations with inexpensive CPU operations" (CPU skeletonization). It cuts overhead by more than 90% against prior topology losses and is the first such loss to handle multi-class segmentation. — [arXiv 2404.03010](https://arxiv.org/abs/2404.03010)
- **Fast-χ** (arXiv 2507.23763) computes the Euler characteristic from bit-quads as a deep-learning-compatible pipeline, 72.7× faster than persistent homology. **Minkowski Image Loss** (arXiv 2604.11422) uses sigmoid-relaxed area, perimeter and χ. **ECC CUDA kernels** (arXiv 2510.20271) claim 16–2000× speedups. All three were covered in round 1. — [arXiv 2507.23763](https://arxiv.org/abs/2507.23763); [arXiv 2604.11422](https://arxiv.org/abs/2604.11422); [arXiv 2510.20271](https://arxiv.org/abs/2510.20271)

**Differentiable polygon pipelines (all in the polygon→raster direction)**
- BoundaryFormer predicts polygons and is trained with mask supervision through a differentiable rasterizer. — [CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Lazarow_Instance_Segmentation_With_Mask-Supervised_Polygonal_Boundary_Transformers_CVPR_2022_paper.pdf)
- DDSL is a differentiable simplex layer that rasterizes simplices and polygons. — [arXiv 1901.11082](https://arxiv.org/pdf/1901.11082)
- BoxSnake uses differentiable polygon rendering via a soft winding number. — [arXiv 2303.11630](https://arxiv.org/pdf/2303.11630)
- A real-time polygon instance segmentation method (CenterPoly successor) adds region-based and order losses to vertex regression. It runs at 0.045 s per 2048×1024 frame on an RTX 2070. — [arXiv 2305.05490](https://arxiv.org/abs/2305.05490)
- Frame Field Learning adds a frame-field output "to a deep segmentation model" so that a separate polygonization step yields regular building polygons. The polygonization itself sits outside the gradient path. — [arXiv 2004.14875](https://arxiv.org/abs/2004.14875)

**Manufacturability (from round 1 and Q1)**
- CurvyILT enforces curvilinear MRC (minimum shape area and minimum shape distance) with **differentiable morphology**. It stops at a pixel mask: 2.11 s per clip, 0.6 GB. — [arXiv 2411.07311](https://arxiv.org/html/2411.07311v1)
- LithoGRPO handles the non-differentiable shot count through RL rather than through any surrogate gradient. — [arXiv 2606.00228](https://arxiv.org/html/2606.00228)

### Inferences
- **H2's design is open, and the field has partly answered it with RL.** For polygon-level objectives (vertex count, shot count, Manhattan-ness, β₁), 2026 practice either (a) predicts vertices and back-propagates through a soft rasterizer, or (b) treats the polygon metric as a black-box reward. A straight-through layer would be a third option: an exact polygon forward pass, with the backward pass routed through the soft bit-quad statistics (Fast-χ, Minkowski loss) or through the per-pixel locations of the offending quads.
  - The engine gives *localization* that soft χ cannot. It identifies which ring is a spurious hole or component and which quads produce the excess vertices. That localization can turn a global count into a pixel weight map, similar to Topograph's critical regions but exact and cheaper.
  - This design is untested and is my inference.
- **Cost argument.** Per-loss topology costs are 96–600 ms (Topograph table), and Betti matching runs on CPU. One GPU pass that gives exact per-instance β₀/β₁ (outer-ring and hole-ring counts) should take milliseconds at most. Its natural uses are (i) a gate that skips the expensive topology loss when β already matches, (ii) a topology metric logged every step, and (iii) a GRPO reward. This is an inference; I have no measurement of the engine against these numbers.

### Gaps
- I did not fetch 2025–2026 papers on "differentiable vectorization" (for example SuperSVG or Bézier splatting) or "PolyGNN" (3D building reconstruction). Search was exhausted.
- Topograph's hardware and image sizes are unknown, so its ms figures cannot be compared directly with a GPU bit-quad pass.
- Betti matching 3D's per-volume runtime numbers were not in the fetched abstract.

---

## Q3. Manufacturability-aware ML (lithography, chips, PCB, fabrication, floorplans)

### Takeaway
Lithography is the clearest case where geometry-in-the-loop is both *needed* and *documented as too slow*. LithoGRPO's shot-count reward needed a 130–490× speedup to make RL practical, and its "fast" version is still a CPU ILP at about 0.2 s per mask. CurvyILT handles MRC with differentiable morphology but never produces polygons. I could not verify, this round, any floorplan, CAD, PCB or 3D-printing generators with polygon-validity rewards.

### Cited Findings
- **LithoGRPO**: reward L2 + PVB + EPE + shot count. Traditional shot count takes 34–148 s per mask; the fast version takes 0.16–0.30 s with R²=0.994. Mask resolution is 512 px. — [arXiv 2606.00228](https://arxiv.org/html/2606.00228)
- **NVIDIA generative-RL ILT** (WGAN + GRPO plus batched ILT refinement) reduces violations and doubles throughput against numerical baselines on LithoBench. Its manufacturability terms are not stated in the abstract. — [arXiv 2602.19027](https://arxiv.org/abs/2602.19027)
- **CurvyILT**: differentiable morphology for MRC, binarization and morphology at the end, no polygon output. — [arXiv 2411.07311](https://arxiv.org/html/2411.07311v1)
- **OpenLithoHub** (2026) is open source and computes shot count and MRC/DRC gates with Manhattan and curvilinear contour extraction. It runs on CPU reference hardware, with GPU timing listed as "pending". This is from round 1. — [GitHub OpenLithoHub](https://github.com/OpenLithoHub/OpenLithoHub)
- **PRISM** (photonics-informed ILT for manufacturable inverse-designed PICs) and **OptiClear** (differentiable curvilinear design-rule legalization) show the same issue in photonics. Only titles and snippets were seen. — [arXiv 2602.15762](https://arxiv.org/pdf/2602.15762); [arXiv 2607.03632](https://arxiv.org/pdf/2607.03632)
- Floorplan generators (HouseDiffusion, GSDiff, DiffPlanner) work in vector space and avoid rasterize-then-vectorize. DiffPlanner cites "information loss and precision issues" as the reason. This is from round 1. — [arXiv 2508.13738](https://arxiv.org/html/2508.13738v1)

### Inferences
- **Lithography is the strongest beachhead for "geometry in the RL loop".** The policy outputs pixels, the deployment cost is polygonal (shots, vertices, OASIS bytes, MRC), RL rewards are computed per rollout, and the current fast proxy is a 0.2 s CPU ILP that ignores polygon structure.
  - A GPU engine could provide exact Manhattan rings, vertex and jog counts, hole counts, and a rectangle-partition bound per rollout, batched over G×batch masks on the device that already holds them.
  - It would also unlock rewards LithoGRPO does not use: vertex count, hole and island count (an MRC proxy), and OASIS size.
- For floorplans and CAD, the vector-native trend means the engine's role is a *validity checker* (rasterize the generated layout, then re-polygonize to detect overlaps, gaps and holes), not a generator component. This is unverified.

### Gaps
- I did not obtain LithoBench's original shot-count algorithm or source.
- RL floorplan or layout generators with validity rewards, PCB ML with DRC rewards, and 3D-printing/CNC slice-level ML with manufacturability rewards were all not researched this round because the search budget ran out.

---

## Q4. VLM polygon tokenization pipelines: who generates polygon targets at scale, and how are holes and multi-part objects handled?

### Takeaway
Polygon targets for VLMs are made on the fly from SAM-style RLE masks with **OpenCV Suzuki-Abe `findContours`** plus simplification, using one outer ring traversed clockwise. **Holes are not represented, and multi-part objects are dropped or truncated.**
- Groundbench keeps only the part with the most vertices. That affects 148 of 1,500 referents (9.9%), and the kept part covers on average only 71.76% of the object's union.
- SimpleSeg lists "object with holes" as a failure case.
- Every large mask source (SA-1B, 1.1B masks; SA-Co, 4M concepts) ships masks, typically as RLE. So mask→polygon conversion is an unavoidable, under-specified preprocessing step.

The evidence points to a quality and representation problem, not a throughput one.

### Cited Findings
- **SimpleSeg** pipeline: Grounding-DINO → SAM masks → "we extract polygonal contours using the Suzuki–Abe algorithm (OpenCV), enforce a consistent traversal order (clockwise), and optionally apply a tolerance-based sparsification".
  - Pre-training data is LAION and Coyo auto-annotated with this pipeline; the count is not stated. SFT uses about 800k RefCOCO-series samples, and RL uses 400k prompts.
  - ε sweep: 35.6 cIoU at 78 tokens, best at about 221 tokens, 72.5 cIoU at 859 tokens.
  - Figure 12 shows failure cases for "object with holes and texture confusion", with no explicit handling of holes or multi-part objects.
  - SFT gives +9.7 on refCOCO, and SFT+RL gives +13.0 over the pre-trained baseline. — [alphaxiv summary of 2601.19228](https://www.alphaxiv.org/abs/2601.19228.md)
- **Groundbench** (arXiv 2609.27821, Aug 2026) re-targets 1,500 RefCOCO/+/g referents to exact-N polygons at five budgets, with N up to 64.
  - **Multi-part:** "When an annotation contains multiple polygon parts, construction selects the part with the most listed vertices. This affects 148/1,500 UIDs; within them, the selected part covers 71.76% of the all-parts union on average."
  - **Target construction:** the longest edge is subdivided to reach N vertices, or Visvalingam–Whyatt simplification reduces to exactly N, deterministically.
  - **Legality:** "requires distinct in-range integer vertices and a simple positive-area polygon". Audits also check for duplicates, zero area and self-intersection.
  - At N=64, Gemini's filled IoU is 37.7, while legal-polygon coverage and legal IoU drop to **7.5 and 4.2**.
  - Headline scores: 88.2 box IoU against 57.7 direct-polygon IoU, and 57.7 against 57.3 when predicted boxes are rasterized against the same contour target.
  - The paper notes that rasterizing an ordered boundary "involves explicit choices about disconnected parts, self-intersections, clipping, and discrete integration". — [arXiv 2609.27821](https://arxiv.org/pdf/2609.27821)
- **SA-1B**: 11M images, **1.1B masks**, about 100 per image, average 1500×2250 px. "Masks are provided in the COCO run-length encoding (RLE) annotation format." — [Meta SA-1B](https://ai.meta.com/datasets/segment-anything/)
- **SAM 3 / SA-Co**: automatically annotated over 4M unique concepts. The SA-Co benchmark has 270K unique concepts, with instance masks and IDs for every matching object. The format was not stated in the fetched page; round 1 found Roboflow serving defaults to RLE. — [GitHub facebookresearch/sam3](https://github.com/facebookresearch/sam3); [arXiv 2511.16719](https://arxiv.org/abs/2511.16719)
- **PolyFormer** (CVPR 2023) generates polygon vertices autoregressively, with a regression decoder that predicts float coordinates "without any coordinate quantization error". Its data scripts (`create_pretraining_data.py`, `create_finetuning_data.py`) were not inspected. — [arXiv 2302.07387](https://arxiv.org/abs/2302.07387); [GitHub polygon-transformer](https://github.com/amazon-science/polygon-transformer)
- **Text4Seg** (ICLR 2025) avoids polygons altogether by giving each of 16×16 patches a semantic descriptor (256 tokens). Row-wise RLE cuts length by 74% and speeds inference 3×. — [arXiv 2410.09855](https://arxiv.org/abs/2410.09855)
- **ALTo** uses 1–32 learned mask tokens, averaging about 17.5. — [arXiv 2505.16495](https://arxiv.org/html/2505.16495v1)
- **VistaLLM**: 32 adaptive contour points give a 97.26 mIoU ceiling on RefCOCO. From round 1. — [arXiv 2312.12423](https://arxiv.org/pdf/2312.12423)

### Inferences
- **Quality problems an exact, hole-aware engine fixes directly.**
  1. Multi-part objects: about 10% of RefCOCO referents in Groundbench's cohort, losing about 28% of their area when truncated.
  2. Holes, which are currently invisible to the policy.
  3. Non-deterministic start vertex and orientation.
  4. Validity: rings from a bit-quad tracer are simple by construction when 4/8-connectivity is resolved consistently.

  A tokenization with ring delimiters (outer/hole) and part separators would be new. I found no paper that uses one.
- **Throughput is probably not the bottleneck.** At an assumed 0.5–2 ms per mask with CPU `findContours`, 1.1B SA-1B masks take about 150–600 CPU-hours, which is embarrassingly parallel. That is cheap next to VLM training. The engine's advantage is **on-the-fly** generation inside GPU data loaders and RL rollouts, with variable ε or vertex budgets per sample (curricula, multiple budgets per target), plus deterministic canonical targets. This is my own estimate, not a measurement.
- Text4Seg and ALTo show the competing path, mask-like tokens. Polygon tokens win on interpretability and composability, but Groundbench shows they lose on legality at high N. That argues for generating hole-aware polygon targets *and* adding legality rewards in RL.

### Gaps
- No source quantifies the fraction of SA-1B or SA-Co masks that have holes or multiple components. A fast engine could measure this directly, and that would be a publishable data-statistics contribution.
- PolyFormer's and VistaLLM's handling of multi-polygons and separators was not verified.
- Gemini's and Qwen3-VL's native polygon formats for multi-part objects were not re-verified this round.

---

## Q5. Evaluation metrics that need polygons, and their cost at scale

### Takeaway
The most widely adopted boundary metric does **not** need polygons. The official Boundary IoU computes the boundary band by `cv2.erode` pixel morphology, with a dilation ratio of 0.02 of the diagonal. That corrects an implication in round 1. Polygon-specific metrics (PoLiS, C-IoU, N-ratio, max tangent angle, legality, vertex-budget compliance) matter for polygon-output models: building extraction and VLM polygon grounding (Groundbench). There, evaluation requires rasterizing predicted polygons, plus validity audits. For mask-output models compared against polygon models, it requires polygonizing masks. I found no published runtime figures for these metrics at scale.

### Cited Findings
- Boundary IoU (CVPR 2021) uses d = 2% of the diagonal and is adopted in COCO/LVIS boundary AP and Boundary PQ. — [arXiv 2103.16562](https://arxiv.org/abs/2103.16562)
- The official implementation's `mask_to_boundary` uses `dilation = int(round(dilation_ratio * img_diag))`, `cv2.erode(new_mask, kernel, iterations=dilation)` with a 3×3 kernel, then `mask - mask_erode`. It is pixel morphology, not contours. — [boundary-iou-api boundary_utils.py](https://raw.githubusercontent.com/bowenc0221/boundary-iou-api/master/boundary_iou/utils/boundary_utils.py)
- Groundbench evaluation:
  - Filled scoring rasterizes "clipped cyclic fill".
  - Legality checks for distinct in-range vertices and a simple positive-area polygon. Audits cover duplicate points, zero area and self-intersection.
  - Rescoring uses alternate reference masks, eight start phases and continuous-area intersection. — [arXiv 2609.27821](https://arxiv.org/pdf/2609.27821)
- Topology metrics are typically Betti errors; Topograph also proposes a strict homotopy-based metric. Costs are roughly 96–603 ms per loss call for PH- or graph-based methods. — [arXiv 2411.03228v2](https://arxiv.org/html/2411.03228v2)

### Inferences
- Note that the erosion with iterations = 2% of the diagonal is O(d) passes per mask: about 54 iterations for a 1500×2250 SA-1B image. A ring-based band (distance-to-ring) would be an alternative, but this is not a strong selling point.
- The stronger evaluation role is **topology and validity auditing at scale**: per-instance β₀/β₁ error, hole recall, legality, vertex count against budget, and Manhattan-ness, across millions of predictions (SA-Co-scale benchmarks, Groundbench-style polygon suites). All of these are exact by-products of one bit-quad plus ring pass. This is an inference.

### Gaps
- I found no runtime numbers for PoLiS, C-IoU, N-ratio or MTA on large benchmarks, and none for Boundary IoU at LVIS/SA-Co scale.

---

## Q6. Verdicts on H1–H3

### Takeaway
- **H1: partially confirmed, and sharper than stated.** RL with geometric rewards on *emitted* polygons exists (SimpleSeg's rasterized IoU, Groundbench's legality axis as an evaluation). No VLM-segmentation RL reward uses validity, β₀/β₁, holes or vertex budget, which makes it under-tooled. The speed argument holds only in lithography, where LithoGRPO documents a 34–148 s shot-count bottleneck.
- **H2: confirmed as absent.** I found no exact-forward, soft-backward polygonizer. The literature routes around it with soft rasterizers, soft bit-quad statistics, or RL.
- **H3: quality confirmed, throughput not.** Holes and multi-part objects are mishandled in every pipeline inspected. I found no evidence that conversion speed limits data generation.

### Cited Findings
- H1 evidence:
  - SimpleSeg uses GSPO with a mask IoU reward on the rasterized polygon, a centroid reward and a format reward, and no validity or topology terms. — [alphaxiv 2601.19228](https://www.alphaxiv.org/abs/2601.19228.md)
  - SAM-R1, Seg-Zero, LENS and VisionReasoner use mask or box IoU and L1 rewards; VisionReasoner explicitly avoids masks in its rewards. — [2505.22596](https://arxiv.org/html/2505.22596); [2505.12081v5](https://arxiv.org/html/2505.12081v5); [2508.14153](https://arxiv.org/html/2508.14153)
  - Groundbench: legal coverage of 7.5% at N=64. — [2609.27821](https://arxiv.org/pdf/2609.27821)
  - LithoGRPO: shot count of 34–148 s fell to 0.16–0.30 s with a CPU ILP. — [2606.00228](https://arxiv.org/html/2606.00228)
- H2 evidence:
  - BoundaryFormer, DDSL and BoxSnake are polygon→raster. — [CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Lazarow_Instance_Segmentation_With_Mask-Supervised_Polygonal_Boundary_Transformers_CVPR_2022_paper.pdf); [1901.11082](https://arxiv.org/pdf/1901.11082); [2303.11630](https://arxiv.org/pdf/2303.11630)
  - Fast-χ and the Minkowski loss are soft counts with no rings. — [2507.23763](https://arxiv.org/abs/2507.23763); [2604.11422](https://arxiv.org/abs/2604.11422)
  - LithoGRPO puts the non-differentiable metric into RL. — [2606.00228](https://arxiv.org/html/2606.00228)
- H3 evidence:
  - SimpleSeg: OpenCV Suzuki-Abe, clockwise, holes as a failure case. — [alphaxiv 2601.19228](https://www.alphaxiv.org/abs/2601.19228.md)
  - Groundbench: largest part only, 148/1,500, 71.76% coverage. — [2609.27821](https://arxiv.org/pdf/2609.27821)
  - SA-1B: 1.1B RLE masks. — [Meta SA-1B](https://ai.meta.com/datasets/segment-anything/)

### Inferences
- **Most defensible product and research angle, ranked:**
  1. A **GPU geometric-reward library for GRPO**: per-rollout validity, β₀/β₁ per instance, hole IoU, vertex and shot-count proxies, and Manhattan-ness, batched over G×B masks or polygons. The lithography (LithoGRPO-style) shot-count proxy is the speed-critical case; VLM segmentation is the reward-novelty case.
  2. A **hole- and multi-part-aware polygon tokenizer and data engine** for VLM training, with a measured dataset statistic such as the share of SA-1B or SA-Co masks with holes or multiple parts.
  3. A **straight-through polygon layer** as a research contribution, novel but unproven.
- **Remaining uncertainty:**
  - The search budget ran out mid-task, so absence claims (no topology rewards in VLM RL, no STE polygonizer, no RL floorplan validity rewards) rest on round-1 plus round-2 searches, not an exhaustive sweep.
  - Whether VLM RL researchers *want* topology rewards is not demonstrated. Groundbench shows the failure mode, but no paper yet trains on legality.
  - The engine's actual speed against `findContours`, CairoSVG and LithoGRPO's fast count is unmeasured, and those benchmarks should be run.

### Gaps
- Not verified this round: UniVG-R1, Visual-RFT and VLM-R1 rewards; Gemini and Qwen3-VL polygon formats; RL floorplan, CAD and PCB validity rewards; 2025–2026 differentiable vectorization papers.
- No measured statistic exists on the share of holes or multi-part objects in SA-1B or SA-Co.
