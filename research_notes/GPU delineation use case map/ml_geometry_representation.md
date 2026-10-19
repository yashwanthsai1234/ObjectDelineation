# ML side of geometry as a representation: where a fast, exact raster→polygon operator sits inside ML loops (2023–2026)

Scope note: this covers the broader ML representation questions. It deliberately skips parallel CCL algorithms and polygon-native segmentation networks (HiSup, FFL, Pix2Poly, etc.), which another researcher covers. "The engine" below means the user's CUDA/C++ operator: binary bitmap in, 2x2 bit-quad classification in one scanline pass, closed rings (outer plus holes) out.

Overall finding: the field converts between raster and vector all the time, but the conversion is almost always a CPU afterthought (OpenCV `findContours` + Douglas-Peucker, marching squares, or "rasterize the vector and compare pixels"). The 3D world has a mature family of differentiable, GPU-native iso-surfacers (DMTet, FlexiCubes, DiffMC, TetWeave, the `diso` CUDA library). I found no 2D equivalent ("FlexiSquares") that is topology-exact for holes and multi-component masks. Bit-quad counting, which the engine already does, has just been revived as a fast differentiable topology loss: the Fast-χ paper (2025) builds Gray's formula into fixed conv kernels and runs **72.7× faster than a persistent-homology loss**. So one GPU pass could yield both a topology signal and the exact contour.

---

## Q1. Implicit → explicit geometry (neural fields, SDFs, Gaussian splatting → mesh/polygon); is there a 2D "FlexiSquares"?

### Takeaway
In 3D, differentiable iso-surface extraction is an established research line with CUDA libraries: MeshSDF, DMTet, FlexiCubes (SIGGRAPH 2023), DiffMC, NMC, GOF's marching tetrahedra for Gaussian splats, and TetWeave (2025). In 2D, the work that turns learned fields into vectors (fonts, SVG) relies on bespoke vectorization steps. I found no general, differentiable, topology-exact 2D polygonizer that serves as a library primitive. That is a real gap.

### Cited Findings
- FlexiCubes (SIGGRAPH 2023, NVIDIA) is "a high quality isosurface representation specifically designed for gradient-based mesh optimization with respect to geometric, visual, or even physical objectives". Extra per-cell parameters (21 per grid cell, 3 per grid point) are optimized together with the scalar field through autodiff. Extraction builds on Dual Marching Cubes for better topology. — [NVIDIA FlexiCubes](https://research.nvidia.com/labs/toronto-ai/publication/2023_siggraph_flexicubes); [TetWeave arXiv 2505.04590](https://arxiv.org/pdf/2505.04590)
- DMTet and FlexiCubes both jointly optimize the implicit field and the grid structure that the mesh is extracted from. DMTet deforms a precomputed tetrahedral grid. — [TetWeave (2025)](https://arxiv.org/html/2505.04590v1)
- The `diso` repo is "a variety of differentiable iso-surface extraction algorithms implemented in cuda" (DiffMC, DiffDMC). This is a 3D GPU library precedent. — [diso README](https://github.com/SarahWeiii/diso/blob/main/README.md)
- DiffMC splits space into a deformable grid, takes a scalar field on its vertices and produces a mesh through differentiable marching cubes. MeshSDF gives a differentiable surface extraction for DeepSDF. Neural Marching Cubes learns vertex positions and local topologies from training meshes to avoid topological errors. — [Neural Mesh Reconstruction survey (Eurographics digital library)](https://diglib.eg.org/items/ecd5cf15-2aef-4c4e-8a5c-96c45b4ede60/full); [NeuManifold WACV 2025 supp.](https://openaccess.thecvf.com/content/WACV2025/supplemental/Wei_NeuManifold_Neural_Watertight_WACV_2025_supplemental.pdf)
- Gaussian Opacity Fields (GOF) extracts geometry from 3D Gaussians "by identifying its levelset, without resorting to Poisson reconstruction or TSDF fusion", using Marching Tetrahedra on grids induced from the Gaussians. Training takes about 24 min for 6 TNT scenes and about 45 min for Mip-NeRF360 "bicycle". 2DGS uses surfels to improve mesh extraction. MVG-Splatting uses plain Marching Cubes. — [GOF arXiv 2404.10772](https://arxiv.org/html/2404.10772v2); [GOF GitHub](https://github.com/autonomousvision/gaussian-opacity-fields)
- 2D fonts as implicits: Multi-Implicit fonts (NeurIPS 2021) represent glyphs as a permutation-invariant set of occupancy functions "without losing features like edges and corners". — [arXiv 2106.06866](https://arxiv.org/pdf/2106.06866)
- Joint Implicit Neural Representation for vector fonts (ICCV 2023) pairs an SDF with a probabilistic **corner field**, plus a custom vectorization step to get compact vector fonts. Standard contouring would produce redundant segments and round off corners. — [Chen et al. ICCV 2023](https://openaccess.thecvf.com/content/ICCV2023/papers/Chen_Joint_Implicit_Neural_Representation_for_High-fidelity_and_Compact_Vector_Fonts_ICCV_2023_paper.pdf)

### Inferences
- Marching squares and 2x2 bit-quad classification work on the same cell (four pixels or samples). The engine is essentially a binary, topology-exact marching squares that resolves the saddle (diagonal) case and links segments into rings in the same pass. Thresholding a 2D field (occupancy net, SDF slice, 2D Gaussian splat opacity, attention or heat map) and then running the engine gives an exact explicit polygon of the level set at pixel resolution.
- The 3D precedent suggests what a 2D "FlexiSquares" would need: (a) sub-pixel vertex placement by interpolating field values along cell edges (the engine produces pixel-corner vertices from a binary input, so this needs the scalar field as well); (b) per-cell learnable offsets; (c) gradients from vertex positions back to field values. The engine's topology pass could supply the discrete connectivity. Positions are then a differentiable function of the field, which is how DMTet and FlexiCubes separate topology from geometry.
- The font papers show why generic contouring is not enough: they need corner and compactness priors. A GPU polygonizer would be the fast first stage. Corner-aware simplification is a separate step.

### Gaps
- I found no published 2D analog of FlexiCubes/DMTet (a differentiable marching-squares library with learnable per-cell parameters). The search turned up only 3D methods. This looks like an open niche, but a negative search result is not proof that nothing exists.
- I found no runtime benchmarks of 2D contour extraction inside neural-field training loops.
- Planar slicing of 3D Gaussian splats or neural fields into 2D polygons (floor plans from scene reconstructions, for example) came up only indirectly. I found no dedicated GPU pipeline.

---

## Q2. Generative models of vector geometry (floorplans, CAD, fonts, SVG, chip layout): how do they bridge raster loss and vector output?

### Takeaway
There are three bridging strategies: (1) generate vectors directly and supervise in coordinate space (HouseDiffusion, GSDiff, DiffPlanner, CAD-Recode, OmniSVG); (2) render vectors to raster and compare pixels, via differentiable rasterizers (DiffVG/LIVE) or non-differentiable rendering plus an RL reward (RLRF); (3) work in raster and vectorize afterwards. Recent papers call this a source of "information loss and precision issues" and are moving away from it. A fast exact raster→vector operator is most useful as (a) the data engine that turns raster datasets into vector training targets, (b) the reverse check in evaluation or reward loops, and (c) the cleanup step for raster-domain generators.

### Cited Findings
- HouseDiffusion (CVPR 2023) denoises 2D coordinates of room and door corners. Floorplans are "1D polygonal loops". It has continuous and discrete denoising targets to enforce parallelism, orthogonality and corner-sharing, and it can control the exact number of corners per room. — [CVPR 2023](https://openaccess.thecvf.com/content/CVPR2023/html/Shabani_HouseDiffusion_Vector_Floorplan_Generation_via_a_Diffusion_Model_With_Discrete_CVPR_2023_paper.html)
- GSDiff (AAAI 2025) treats floorplans as a structural graph (wall junctions as nodes, walls as edges). Diffusion generates the nodes, a Transformer predicts the edges, and an alignment loss improves geometric consistency. — [arXiv 2408.16258](https://arxiv.org/pdf/2408.16258); [AAAI](https://ojs.aaai.org/index.php/AAAI/article/view/33904)
- DiffPlanner (2025) generates floor plans "entirely within vector space". It explicitly criticizes prior methods that "rely on a complex and redundant rasterization-to-vectorization pipeline, which often leads to information loss and precision issues". — [arXiv 2508.13738](https://arxiv.org/html/2508.13738v1)
- CAD-Recode (ICCV 2025) turns point clouds into executable CadQuery Python (sketch-extrude). It is trained on 1M procedurally generated sequences. At inference it samples 10 candidates, executes each and keeps the one with minimum Chamfer distance to the input. On DeepCAD it reports about 10× better mean Chamfer distance and IoU from 77.6% to 92.0%. — [arXiv 2412.14042](https://arxiv.org/abs/2412.14042v2); [ICCV PDF](https://orbilu.uni.lu/bitstream/10993/67036/1/CADRecode_ICCV_2025.pdf)
- DualVector (CVPR 2023) models glyphs as pairs of closed "positive" and "negative" paths combined by boolean operations. It learns from glyph **images** and then runs a contour-refinement step that aligns the contour to the image representation. — [CVPR 2023](https://openaccess.thecvf.com/content/CVPR2023/html/Liu_DualVector_Unsupervised_Vector_Font_Synthesis_With_Dual-Part_Representation_CVPR_2023_paper.html)
- OmniSVG tokenizes SVG commands and coordinates into discrete sequences for a pretrained VLM, and introduces the MMSVG-2M dataset (2M SVGs). StarVector feeds CLIP image features to a code LLM. — [OmniSVG arXiv 2504.06263](https://arxiv.org/html/2504.06263v3); [HF OmniSVG](https://huggingface.co/OmniSVG/OmniSVG)
- RLRF (Reinforcement Learning from Rendering Feedback, 2025) renders the SVG rollouts, compares them with the target image and uses that as the RL reward. This is a raster loss on vector output, with no differentiable rasterizer needed. — [arXiv 2505.20793](https://arxiv.org/html/2505.20793v2)
- VFIG (2026) vectorizes complex figures with VLMs, using a rubric VLM judge that scores rendered predictions on Presence, Layout and Connectivity. — [arXiv 2603.24575](https://arxiv.org/pdf/2603.24575)
- Optimization-based vectorization (LIVE plus the DiffVG differentiable rasterizer) is slow. At 1080×1920 on an RTX 3090 Ti, one 200-iteration layer took 212 s and 28 layers took about 2 h. At a 256 px max side, 32 layers took 22 min. — [Image Vectorization: a Review, arXiv 2306.06441](https://arxiv.org/pdf/2306.06441)
- DREAMPlace (chip placement) recasts analytical placement as neural-network training in PyTorch. It rasterizes cell density into a grid of rectangular bins and treats density as an electrostatic system solved by FFT-based Poisson. It reports over 30× speedup over RePlAce. — [DREAMPlace (UT Austin)](https://www.cerc.utexas.edu/utda/publications/C252.pdf); [IEEE CEDA](https://ieee-ceda.org/media/dreamplace-deep-learning-toolkit-enabled-gpu-acceleration-modern-vlsi-placement); [NVIDIA AutoDMP blog](https://developer.nvidia.com/blog/autodmp-optimizes-macro-placement-for-chip-design-with-ai-and-gpus/)

### Inferences
- **Data engine:** OmniSVG-style models need millions of vector targets (MMSVG-2M). Floorplan, building and map generators are trained on vector datasets that are often derived from raster masks. A GPU polygonizer that is exact with holes can turn large raster or segmentation corpora (SA-1B-style masks, rendered maps, game-level occupancy grids) into ring-structured vector training data at throughput comparable to the data loader.
- **Reverse validity check for vector generators:** generated polygons can self-intersect or have invalid holes (see Groundbench's "legality" metric in Q3). Rasterizing the generated geometry and re-polygonizing it is a cheap canonicalization that yields a valid, hole-aware ring set. A GPU version would allow this per RL rollout. RLRF and CAD-Recode already use render-and-compare or execute-and-compare loops.
- **Placement and layout:** DREAMPlace works on rasterized density bins. Turning density or overflow maps, or occupancy of macro regions, into polygons (for legalization, blockage export or GNN features) is a natural raster→vector hand-off. I found no source showing anyone does this inside the loop today (inference only).
- Raster-domain generators (diffusion that outputs masks or images) still need vectorization for CAD/GIS/font export. The trend in DiffPlanner and GSDiff is to avoid it, which suggests a vectorizer's value here lies in data generation and evaluation rather than in the generator's output head.

### Gaps
- I did not verify whether MaskPLAN, SketchGraphs/Vitruvion, Text2CAD or DeepVecFont-v2 use rasterize/vectorize steps internally. No sources were fetched for them in this pass.
- I found no source that measures vectorization cost as a bottleneck in SVG dataset construction such as MMSVG-2M.
- I found no source on AlphaChip's internal representation for this pass. It is known to place macros on a grid canvas, but that is not verified here.
- I found no map or road-network generator details in this pass.

---

## Q3. Polygon tokenization for transformers and VLMs: how do current models emit geometry, and how much do polygons save over RLE or masks?

### Takeaway
Polygon-as-text-tokens is now a mainstream way for MLLMs to output segmentation (VistaLLM, SimpleSeg/Kimi-VL, Qwen2.5-VL fine-tunes). Training data comes from **OpenCV `findContours` + Douglas-Peucker** applied to SAM masks. Production APIs are inconsistent: SAM 3 returns RLE or polygons, while Gemini 2.5 returns base64 PNG probability masks inside boxes. VLMs still ground with polygons far worse than with boxes (Groundbench: 57.7 vs 88.2 IoU), and polygon "legality" is a measured failure mode. Fast, exact, hole-aware, consistently oriented contour extraction is the hidden preprocessing step behind all of this.

### Cited Findings
- VistaLLM (CVPR 2024) turns binary masks into point sequences with a "gradient-aware adaptive contour sampling". With 32 points, the maximum achievable mIoU on RefCOCO val is 97.26 (adaptive) vs 94.70 (uniform). Going from 16 to 32 points improves results by +3.6 on Ref and +4.5 on Ref+. — [arXiv 2312.12423](https://arxiv.org/pdf/2312.12423); [CVPR 2024](https://openaccess.thecvf.com/content/CVPR2024/papers/Pramanick_Jack_of_All_Tasks_Master_of_Many_Designing_General-Purpose_Coarse-to-Fine_CVPR_2024_paper.pdf)
- SimpleSeg ("Simple Point Prediction for Pixel-Level VLM Perception", arXiv 2601.19228, Moonshot/Kimi-VL team) outputs boundary coordinates as plain text tokens with no mask decoder. Its data pipeline is Grounding-DINO → SAM mask → contour extraction "enforc[ing] a consistent clockwise traversal order" with optional sparsification. The best tolerance ε=0.005 gives about 221 tokens per object; 78 tokens (too sparse) and 859 tokens (too dense) both perform worse. It reports 87.2 cIoU on refCOCO, and SFT→RL adds +9.7 to 10.5 gIoU. — [hyper.ai summary of 2601.19228](https://hyper.ai/ja/papers/2601.19228)
- A Qwen2.5-VL polygon fine-tuning recipe uses OpenCV `findContours` + Douglas-Peucker, quantizes vertices to [0,255], and wraps them in `<seg>` delimiters. The same source says "a thin line segment requiring ~200 RLE tokens can be represented by ~10 polygon coordinate tokens", a 5 to 10× reduction, and that "clockwise point ordering is critical; unordered sequences cause chaotic outputs". — [Search-result summary attributing these to Qwen2.5-VL polygon work / SimpleSeg (2601.19228)](https://hyper.ai/fr/papers/2601.19228). Caution: the 200 vs 10 token figure came from a search snippet, and I could not confirm which paper it comes from.
- Groundbench (arXiv 2609.27821, 2026) retargets 1,500 referring triples to exact polygons at 5 vertex budgets and scores filled-region IoU and **polygon legality** separately. The strongest configuration reaches 88.2 box IoU (97.1% Acc@0.5) but only 57.7 IoU (69.2%) with direct polygons, and performance "collapses at the densest budget". — [arXiv 2609.27821](https://arxiv.org/pdf/2609.27821)
- ALTo (NeurIPS 2025) is the competing approach: an adaptive-length *mask tokenizer* with a token-length predictor and differentiable token chunking, using GRPO to trade quality against token cost. — [arXiv 2505.16495](https://arxiv.org/html/2505.16495v1); [NeurIPS 2025](https://neurips.cc/virtual/2025/poster/116479)
- SAM 3 (Meta, Nov 2025) supports open-vocabulary "promptable concept segmentation" of all instances. Through the Roboflow integration it can output `rle` (default) or `polygons`. It takes about 30 ms per image with 100 objects on an H200. — [Roboflow SAM3 docs](https://docs.roboflow.com/models/supported-models/sam3); [Ultralytics SAM 3](https://docs.ultralytics.com/models/sam-3); [InfoQ](https://www.infoq.com/news/2025/11/meta-sam3)
- Gemini 2.5 returns segmentation as JSON entries with `box_2d` plus `mask`, a base64 PNG probability map (0 to 255) cropped to the box. The format varies across model versions. — [Simon Willison, Apr 2025](https://simonwillison.net/2025/Apr/18/gemini-image-segmentation/); [Google AI forum](https://discuss.ai.google.dev/t/image-understanding-and-segmentation-mask-support/146559)

### Inferences
- **Where the engine fits:** (1) Offline tokenizer for training data: millions of SAM/SAM 3 masks → rings → simplified, quantized token sequences. Today this is CPU OpenCV and is part of the data-prep cost. (2) Online in RL fine-tuning (SimpleSeg's RL stage, ALTo's GRPO): rewards need fast rasterize↔polygonize round-trips per rollout. (3) In serving: converting SAM 3 masks or Gemini PNG probability maps into polygons for downstream LLM reasoning or GIS export, at around 100 objects per frame at 30 ms per frame, where a CPU `findContours` per object becomes a noticeable fraction of latency. That last point is an inference that should be benchmarked.
- **Holes and multi-part objects are an unaddressed problem:** every pipeline above serializes one clockwise outer contour. The engine emits outer plus hole rings with consistent orientation, which suggests a representation (ring delimiter tokens, hole tokens) that current polygon-token MLLMs lack. This matters for donut-shaped objects, building courtyards and lesions with cavities.
- **Is boundary vectorization a cheap way to give LLMs geometry?** Yes for compactness: roughly 100 to 250 tokens per object at good fidelity per SimpleSeg, against 5 to 10× more for RLE per the snippet. But Groundbench shows that *generating* legal polygons is hard for VLMs. A practical hybrid is to have the VLM emit a mask (or ALTo tokens) and let a deterministic polygonizer produce legal geometry. That puts a fast exact polygonizer in the inference path.

### Gaps
- I found no apples-to-apples published table comparing polygon vs RLE vs mask-token counts on one dataset, beyond the SimpleSeg ε sweep and the unattributed "200 vs 10" snippet.
- Molmo pointing, Florence-2 polygon output and PolyFormer token budgets were not verified in this pass.
- I found no paper on hole-aware or multi-ring polygon tokenization for MLLMs. This appears open.

---

## Q4. ML that consumes polygons (shape encoders, GNNs on planar partitions): how do they get the polygons?

### Takeaway
Polygon encoders are maturing: Poly2Vec at ICML 2025 is a Fourier-based unified encoder for points, lines and polygons, and PolygonGNN at KDD 2024 encodes multipolygons with holes through heterogeneous visibility graphs. They assume clean vector input from GIS databases (OSM, building footprints). Polygons extracted from model masks reach them only through ad hoc conversion.

### Cited Findings
- Poly2Vec (ICML 2025) is a "polymorphic Fourier-based encoding" that handles points, polylines and polygons in one framework. Magnitude encodes extent and phase encodes direction, combined by a learned fusion module. It beats object-specific baselines at preserving topology, direction and distance relationships, and improves population prediction and land-use inference. — [PMLR v267](https://proceedings.mlr.press/v267/siampou25a.html); [arXiv 2408.14806](https://arxiv.org/html/2408.14806v2)
- PolygonGNN (KDD 2024) is built on a heterogeneous visibility graph. Inner edges follow each ring, both the exterior and **holes**. Visibility edges connect vertices of different polygon parts that can see each other. Spanning-tree sampling keeps it efficient. Tasks include shape coding, building-pattern classification and geographic QA. — [arXiv 2407.00742](https://arxiv.org/pdf/2407.00742); [HF model](https://huggingface.co/Dzy6/PolygonGNN)

### Inferences
- PolygonGNN takes exactly the structure the engine outputs: rings grouped as outer plus holes, in vertex order. Poly2Vec's Fourier transform of a polygon has a closed form over edges (NUFT-style), so it also consumes rings directly. A GPU polygonizer can therefore act as a **mask → shape-embedding adapter**: segmentation output → rings → Poly2Vec/PolygonGNN features → downstream reasoning (building typology, land-use, cell morphology), all on the GPU without a CPU round-trip.
- Region adjacency graphs and planar partitions (panoptic or land-cover maps) would need shared-boundary extraction between labels. That is a multi-label extension the binary engine does not do as described, but a bit-quad approach could extend to it.

### Gaps
- I found no paper that chains segmentation network → polygonization → polygon encoder end-to-end (trained jointly). This looks like an open opportunity, since such a chain would need a differentiable polygonizer or a straight-through estimator.
- I did not verify Space2Vec or NUFT polygon encoder details in this pass.

---

## Q5. Topology and morphology as learning signals and metrics: can a GPU bit-quad pass give topology metrics or losses for free?

### Takeaway
Yes, and there is now direct evidence for it. "Topology Optimization in Medical Image Segmentation with Fast χ Euler Characteristic" (arXiv 2507.23763, 2025) uses exactly Gray's bit-quad formula as fixed conv kernels, makes it differentiable with soft probabilities, and reaches **0.0134 s/item vs 0.9748 s/item for a persistent-homology loss (72.7×)**, roughly the cost of clDice (0.0151 s). Betti matching remains much slower even after C++ optimization. The engine's pass already classifies every 2x2 quad, so χ, area and perimeter come out as by-products of polygonization.

### Cited Findings
- The Fast-χ paper computes 2D χ with bit-quads: **χ = ¼(#Q_A − #Q_B − 2#Q_C)**, implemented as fixed convolution kernels. The 3D version uses 2×2×2 "bit-octets", with 256 patterns reduced to 22 symmetry groups. It also produces **topological violation maps** that localize χ errors and feed a correction network. — [arXiv 2507.23763 (HTML)](https://arxiv.org/html/2507.23763)
- Runtimes on CREMI: Fast-χ 0.0134 s/item; PH loss 0.9748 s/item (72.7× slower); Boundary loss 0.0168; Warp loss 0.0151; clDice 0.0151 s/item. — [arXiv 2507.23763](https://arxiv.org/html/2507.23763)
- Betti error (L1 on Betti numbers) falls against a Swin-UNet baseline: CREMI 7.127→4.941; ISBI13 5.643→2.330; FIVES 40.124→26.562; dHCP 3.852→2.032; TopCoW (3D) 6.667→2.267. Dice is preserved, and all results have p<0.05. — [arXiv 2507.23763](https://arxiv.org/html/2507.23763)
- An earlier MICCAI 2023 Euler-characteristic topology paper by the same group (Li et al.) has code available. — [GitHub smilell/Topology-aware-Segmentation-using-Euler-Characteristic](https://github.com/smilell/Topology-aware-Segmentation-using-Euler-Characteristic); [MICCAI 2023](https://conferences.miccai.org/2023/papers/556-Paper0354.html)
- Efficient Betti Matching (2024) is a C++ implementation for arbitrary dimensions, up to about 5× faster than prior PH approaches, which first made 3D Betti-matching training feasible. The authors note it is still CPU-bound and suggest GPU work in the style of Ripser++. — [arXiv 2407.04683](https://arxiv.org/pdf/2407.04683)
- GPU parallelization can cut PH run-time scaling from O(N⁴) to O(N³) or O(N²). — [arXiv 2203.02527](https://arxiv.org/pdf/2203.02527)
- Boundary IoU (CVPR 2021) restricts IoU to a band of width d (typically 2% of the image diagonal) around each mask's contour. It is the basis of Boundary AP and Boundary PQ, and it is more sensitive than Mask IoU for large objects. — [arXiv 2103.16562](https://arxiv.org/abs/2103.16562)
- HyperTopo-Adapters (2026) is another recent geometry- and topology-aware segmentation adapter (leaf lesions) that uses differentiable topology surrogates. — [arXiv 2601.06067](https://arxiv.org/pdf/2601.06067)

### Inferences
- **"For free":** the engine already computes the per-pixel 2x2 pattern class that Gray's χ formula counts, so summing over pattern classes gives χ (components minus holes) at essentially zero extra cost. It is exact at 8-connected or 4-connected foreground depending on the formula variant, and the rings make the counts explicit: #outer rings and #hole rings give β₀ and β₁ *separately*. That is strictly more than χ, and it is what the Fast-χ paper cannot get from χ alone (χ = β₀ − β₁ cannot tell 2 components + 1 hole from 1 component + 0 holes). Area (pixel count) and perimeter (bit-quad edge counts, Gray's estimators) come from the same classification. Together these are the 2D Minkowski functionals (area, perimeter, χ).
- **Losses:** a forward pass on the GPU that counts β₀ and β₁ per image, plus a localized violation map (which rings are spurious), could supply targets or weights for a soft loss. The non-differentiable exact counts go into the metric or weighting, and a soft bit-quad surrogate carries the gradient. The suggested design is to use the engine for exact topology and localization, and Fast-χ-style soft kernels for gradients.
- **Metrics that need polygons:** Boundary IoU, PoLiS, C-IoU, N-ratio, MTA and complexity-aware building metrics all need vector boundaries or contour bands. Large-scale evaluation (millions of predicted instances, SAM-scale) would benefit from GPU polygon extraction feeding GPU metric kernels.

### Gaps
- I did not fetch exact runtime numbers for Betti matching per 2D image (only "up to 5× faster"), nor runtime for computing Boundary IoU at scale.
- Fast-χ reports χ only. I found no paper using separate β₀/β₁ counts from contour tracing as a training signal. It appears open but is not confirmed.
- PoLiS, C-IoU and MTA definitions were not re-fetched here. The polygon-segmentation researcher likely covers them.

---

## Q6. Agents, world models and RL on spatial layouts: do they need vectorized geometry?

### Takeaway
LLM GIS agents in 2025–2026 operate through tool calls on vector and raster data (GeoBenchX, GeoAgentBench). Raster→vector conversion is a standard GIS tool step inside those workflows. I found no evidence that these benchmarks stress conversion speed, but agents chaining "segment → polygonize → spatial join" would call such an operator repeatedly.

### Cited Findings
- GeoBenchX (2025) benchmarks 8 commercial LLMs as tool-calling agents with 23 geospatial functions, over tasks spanning tabular, **vector and raster** data, multi-stage spatial analysis and hallucination rejection. — [arXiv 2503.18129](https://arxiv.org/html/2503.18129v2)
- GeoAgentBench (2026) is an execution sandbox with **117 atomic GIS tools** and 53 spatial-analysis tasks across 6 GIS domains. It criticizes static text or code matching and evaluates dynamic runtime feedback and multimodal spatial outputs. — [arXiv 2604.13888](https://arxiv.org/pdf/2604.13888)

### Inferences
- When the agent's perception is a segmentation model (SAM 3, a VLM) and its reasoning is in GIS terms (area, adjacency, containment, buffer), a polygonizer is the bridge. A GPU, topology-exact version would let an agent's tool return clean vector layers in milliseconds instead of a GDAL `Polygonize` CPU step. That last point is an inference; I did not benchmark GDAL here.
- RL agents on occupancy grids (robot maps, game levels, layout design) could use polygonized free space and obstacles as compact observations. This is speculative and has no source.

### Gaps
- I found no source documenting which raster→vector tool GeoAgentBench or GeoBenchX expose, or its latency.
- I found no world-model or RL paper that explicitly uses polygonized observations in 2023–2026 in this pass.

---

## Q7. Compression and streaming of segmentation outputs (edge → cloud contours)

### Takeaway
Sending contours instead of masks is an old idea, found in patents and in MPEG-4 shape coding. Today's mask APIs use RLE (SAM 3 default) or PNG (Gemini). I found no 2023–2026 ML paper that benchmarks polygon streaming of segmentation for edge-to-cloud links. That is either a gap or simply too engineering-level to publish.

### Cited Findings
- A patent describes segmenting on the edge and transmitting reduced data to the cloud because image transfer is a bottleneck under constrained bandwidth. — [US11336741B2](https://patents.google.com/patent/US11336741)
- A patent on contour display notes that "a typical contour algorithm can generate tens of thousands of polygons that require a lot of memory to store, high levels of bandwidth to transmit", which motivates point-reduction simplification. — [USPTO 10753848 "Efficient contours and gating"](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10753848)
- Piecewise linear approximation is used at the edge to cut transmitted data. — [Chalmers](https://research.chalmers.se/en/publication/?id=510345)
- SAM 3 serving defaults to RLE "which is more memory-efficient", with polygons as an option. — [Roboflow SAM3](https://docs.roboflow.com/models/supported-models/sam3)

### Inferences
- For compact objects, a simplified polygon costs O(perimeter vertices) against O(rows × runs) for RLE. The SimpleSeg numbers (about 221 coordinates per object at near-SOTA fidelity) suggest large savings over raw or RLE masks for typical objects. On-device GPU polygonization (Jetson-class) right after the segmenter would make it practical for drones and satellites, where downlink is scarce. The engine also preserves holes exactly, which simplified-outer-only pipelines lose.

### Gaps
- I found no measured bandwidth comparison (polygon vs RLE vs PNG) for modern segmentation outputs in 2023–2026 literature.
- I found no on-board satellite or drone ML paper that reports vectorization latency.

---

## Cross-cutting GAP list (what does not exist yet, as far as this search found)
1. **2D differentiable iso-contouring library** ("FlexiSquares/DiffMS") with learnable per-cell parameters and exact hole topology. The 3D analogs (FlexiCubes, DMTet, `diso` CUDA DiffMC) exist.
2. **Hole-aware polygon tokenization for MLLMs.** Current pipelines (SimpleSeg, VistaLLM, Qwen fine-tunes) serialize a single clockwise outer contour from OpenCV.
3. **Separate β₀/β₁ (not just χ) GPU topology signals** from one pass. Fast-χ (2025) proves the bit-quad route is about 70× cheaper than PH. Ring-level counting would add localization and hole/component separation.
4. **End-to-end mask → polygon → polygon encoder** (Poly2Vec/PolygonGNN) pipelines on the GPU.
5. **Benchmarks of polygonization cost** in VLM data prep, RL reward loops (RLRF, SimpleSeg RL, ALTo GRPO), SAM-scale evaluation, and edge streaming.
6. **Legal-polygon guarantee for VLM outputs:** Groundbench shows legality and density collapse. A deterministic mask→polygon back end is a plausible fix but has not been evaluated.
