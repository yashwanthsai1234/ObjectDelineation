# Algorithms + ML-Methods Frontier for GPU Mask-to-Polygon Conversion (2025–2026)

Scope: (a) parallel algorithms that could make bitmap → closed rectilinear rings (with hole nesting) fully GPU-native; (b) ML methods whose output is polygons, and where an exact, fast, topology-correct GPU polygonizer would plug in. Research date: 2026-10-07. Some full texts were paywalled (MDPI 403, Springer login); those items are flagged.

---

## Q1. GPU connected-component labeling, contour tracing, marching squares, and parallel ring linking: what is fastest, and does an open-source GPU polygonizer with hole nesting exist?

### Takeaway
GPU connected-component labeling (CCL) is a solved, benchmarked problem: block-based union-find (BKE for 2D, BUF for 3D) is fastest, and BKE ships in OpenCV CUDA (≥4.6). GPU contour/border tracing with hole hierarchy is not solved in the same way. The work that exists is academic (block-wise Suzuki variants, "triad" stitching, a 2025 GPU polygon-reconstruction paper with holes). None of NVIDIA's or the mainstream GPU libraries I checked exposes a GPU `findContours` that returns nested rings. The serial CPU ring-linking step the user still has is exactly the gap the field has left open.

### Cited Findings
**CCL state of the art**
- The 2024 IEEE TPDS survey by Bolelli, Allegretti, Lumetti and Grana ("A State-of-the-Art Review with Code about Connected Components Labeling on GPUs") benchmarks about 20 GPU CCL algorithms. These include OLE, UF, LBUF, KE, BE, DLP, STAVA, 8DLS/M8DLS, Rasmusson, HA4/HA8, the decision-tree ports C_SAUF/C_BBDT/C_DRAG, BUF and BKE. The tests span several OSes, compilers and GPUs (Quadro P1000 up to A100 80GB). Its conclusion: "block-based oriented algorithms outperform all the other algorithmic solutions on both 2D images and 3D volumes, regardless of the selected environment." — [Bolelli et al., TPDS 2024 PDF](https://federicobolelli.it/media/publications/pdfs/2024tpds.pdf)
- BKE (Block-based Komura Equivalence) "is always the best-performing algorithm on 2D datasets" (text, fingerprints, medical, surveillance). BUF (Block-based Union-Find) is usually better on 3D volumes, except on the A100, where BKE stays competitive. Both use "inline compression" (path compression during find). That saves traversal but adds writes. — [Bolelli et al. 2024](https://federicobolelli.it/media/publications/pdfs/2024tpds.pdf)
- BKE "has been available in OpenCV since version 4.6.0" (cv::cuda::connectedComponents) and is the only GPU CCL in OpenCV. — [Bolelli et al. 2024](https://federicobolelli.it/media/publications/pdfs/2024tpds.pdf)
- Example magnitudes from the survey tables (Quadro P1000, Fingerprints dataset, average ms per image): BKE ≈ 0.44, BUF ≈ 0.45, C_DRAG ≈ 0.47, HA8 ≈ 0.48, KE ≈ 0.51, 8DLS ≈ 3.0, M8DLS ≈ 13.1. On 2D natural images the top block-based algorithms sit within about 10% of each other. — [Bolelli et al. 2024, Table 3](https://federicobolelli.it/media/publications/pdfs/2024tpds.pdf)
- HA8 (Hardware-Accelerated 8-connectivity) uses CUDA intrinsics to find runs of foreground pixels, then union-find merges the runs. It is near state of the art, especially for horizontally extended components. — [Bolelli et al. 2024](https://federicobolelli.it/media/publications/pdfs/2024tpds.pdf)
- Decision-tree algorithms (C_DRAG etc.) blow up at around 50% density with granularity 1 because of branch divergence. They do well at <10% or >90% density, which is typical of real masks. — [Bolelli et al. 2024](https://federicobolelli.it/media/publications/pdfs/2024tpds.pdf)
- All code is in YACCLAB (TPDS2024 release, v3.4-alpha). The authors had to reimplement many published algorithms because the original code was not available. — [YACCLAB releases](https://github.com/prittt/YACCLAB/releases); [Bolelli et al. 2024](https://federicobolelli.it/media/publications/pdfs/2024tpds.pdf)

**GPU contour tracing / border following**
- "Parallel border tracking in binary images using GPUs" (J. Supercomputing, 2021) is a block-wise parallel border tracker. — [Springer s11227-021-04260-y](https://link.springer.com/article/10.1007/s11227-021-04260-y)
- "Accelerated border tracking in binary images with GPUs" (J. Supercomputing, 2026) splits the image into rectangular blocks. Each block extracts "triads" (three interconnected, ordered points) in parallel. Triads are then connected within blocks and across block boundaries into closed contours, using three CUDA kernels. Its stated main goal is "to avoid costly data transfers between the CPU and GPU, while maintaining performance at a level similar to that of the CPU." So the headline gain is avoiding host-device transfers, not raw speedup over CPU tracing. — [Springer s11227-026-08488-4](https://link.springer.com/10.1007/s11227-026-08488-4)
- An earlier CUDA Suzuki-style parallel contour tracer splits the image into small rectangles and uses block searching. — [IJSSST Vol-17 paper](https://ijssst.info/Vol-17/No-19/paper1.pdf); [ResearchGate](https://www.researchgate.net/publication/309767994_A_Parallel_Approach_for_Contour_Extraction_Based_on_CUDA_Platform)
- "GPU-Accelerated Algorithm for Polygon Reconstruction" (Applied Sciences 15(3):1111, 2025). Reported average speedups over CPU: 46.1× on non-Manhattan test sets without holes, 45.2× with holes, 17.9× on Manhattan (rectilinear) sets without holes, 16.8× with holes. It explicitly matches inner and outer contours to rebuild polygons with holes. This looks like the closest prior art to the user's engine. Full text returned 403, so the domain (likely IC layout), GPU model and code availability are unverified. — [doi:10.3390/app15031111](https://doi.org/10.3390/app15031111) (search-snippet level)
- "Real-Time Contour Image Vectorization on GPU" (PKU Graphics Lab) is another GPU vectorization pipeline. — [PKU PDF](https://www.icst.pku.edu.cn/GraphicsLab/docs/20200416023556073458.pdf)

**GPU marching squares + stitching**
- The standard GPU pattern for marching-squares/cubes output is two kernels: count the vertices each cell emits, run an exclusive prefix sum to get offsets, then write. This is "embarrassingly parallel." — [WebGPU Marching Cubes, W. Usher 2024](https://www.willusher.io/graphics/2024/04/22/webgpu-marching-cubes/); [Dyken & Ziegler, HistoPyramid MC, GTC 2010](https://www.nvidia.com/content/gtc-2010/pdfs/gtc10_session2020_hp5mc_dyken_ziegler.pdf)
- A USPTO patent ("Determining contour edges for an image") describes GPU marching squares to extract resist-contour edges for OPC (optical proximity correction). The extracted edges go into GPU hierarchical search structures, and the search results are used to stitch contours. This shows that computational lithography already needs GPU contour stitching. — [USPTO 12315131](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12315131)
- NVIDIA cuLitho claims 57× end-to-end acceleration and "reinvents ... computational geometry" primitives. The public abstract gives no details of its contour algorithms. — [arXiv 2602.15036](https://arxiv.org/abs/2602.15036)

**Library coverage**
- OpenCV has no CUDA counterpart of cv::findContours. The forum consensus is that tracing is "not very parallelizable for large objects." — [OpenCV Q&A](https://answers.opencv.org/question/84942/isnt-there-a-opencv-cuda-function-similar-to-findcontours/); [NVIDIA Jetson forum](https://forums.developer.nvidia.com/t/findcontours-call-using-multi-threaded-cuda-application/53486)

### Inferences
- **A fully GPU-native rectilinear ring linker looks feasible with standard primitives (my design inference, not found in the literature):**
  1. Classify each 2×2 quad (the user already does this). Vertex count = n(Q1) + n(Q3) + 2·n(QD).
  2. Compact the vertices with a prefix sum.
  3. On any horizontal grid line, boundary vertices sorted by x pair up as (1,2), (3,4), …, giving each vertex its horizontal-edge partner. Sorting by (x, y) gives the vertical partner the same way. The diagonal (QD) quads need a deterministic split rule (4- vs 8-connectivity) that emits two vertices.
  4. This makes every vertex a node of degree 2 in an H/V-alternating graph. Each ring is a cycle, so ring ID, ring order and per-ring compaction are list-ranking / pointer-jumping problems: Wyllie gives O(n log n) work; sparse-ruling-set or Helman–JáJá gives O(n) work.
  5. Shoelace area by segmented reduction gives orientation, which separates outer rings from holes.
  
  Every step is a radix sort, a scan or pointer jumping, all of which CUB/Thrust already provide.
- **Hole nesting without point-in-polygon tests (inference):**
  - Run BKE on the foreground (8-connected) and on the background (4-connected, or the dual).
  - Each foreground component is exactly one outer ring. Each bounded background component is exactly one hole.
  - A hole's parent is the foreground label of any pixel adjacent to it.
  - This gives the full outer/inner hierarchy as a by-product of two CCL runs (Suzuki-Abe-equivalent topology) at BKE speed (sub-millisecond for typical images per the survey).
  - The parent of a hole ring is found by reading the CCL label at one of its vertices.
- The 2026 triad paper reaching only CPU parity suggests that general (non-rectilinear, 8-direction chain) border following on GPU is still hard. The user's rectilinear setting is easier: vertices have fixed successors determined by local pattern plus a sorted-neighbor lookup. So a big GPU win is plausible here even where general tracing struggles.

### Gaps
- I could not confirm whether cuCIM (cucim.skimage.measure), CV-CUDA, VPI, NPP, RAPIDS cuSpatial, Kornia or TorchVision provide a GPU `find_contours` or GPU polygonizer with hole nesting. Searches returned nothing, and nothing found suggests any of them does. This needs a direct docs check, but it looks like a real ecosystem gap.
- I found no published GPU algorithm that frames rectilinear ring linking explicitly as sort + pair + list ranking. GPU list ranking itself is established (Wyllie; Rehman/Kothapalli/Narayanan ICS 2009), but I did not verify a source this session.
- The Applied Sciences 2025 paper's method details (linking strategy, GPU, code) were not retrievable because of the 403.
- No benchmark compares GPU contour extractors against each other, the way YACCLAB does for CCL. A "YACCLAB for polygonization" does not exist.

---

## Q2. Polygon simplification on GPU, topology-preserving simplification with shared edges, and multi-label (label map → planar partition) polygonization

### Takeaway
Douglas-Peucker has been parallelized on CPUs and GPUs at the level of one feature per thread. Topology-preserving and shared-boundary simplification exists only on CPUs and is mostly sequential (constrained DP, TopoJSON-style arc sharing). The newest ML work (Sep 2026) explicitly cites independent per-instance raster-to-vector conversion as the cause of gaps and overlaps. It fixes this with a common planar graph, but not on the GPU.

### Cited Findings
- GPU-accelerated Ramer-Douglas-Peucker is reported to be at least an order of magnitude faster than multi-CPU versions. Parallelism comes from splitting a polyline into independent segments, or from one polyline per thread. — [ResearchGate: real-time parallel DP](https://www.researchgate.net/publication/251961281_A_real-time_parallel_implementation_of_Douglas-Peucker_polyline_simplification_algorithm_on_shared_memory_multi-core_processor_computers) (search-snippet level)
- Contextual / constrained Douglas-Peucker uses other geometries as constraints. It "incrementally rewinds" to the original polyline to resolve conflicts and avoids self-intersections and topology changes. It is CPU-only. — [Contextual Douglas-Peucker Simplification](https://www.researchgate.net/publication/283800669_Contextual_Douglas-Peucker_Simplification); [topology-preserving line simplification](https://proceedings.ictinnovations.org/attachment/paper/578/a-topology-preserving-line-simplification-algorithm.pdf)
- GPU topology-preserving simplification does exist for triangle meshes: CUDA edge contraction with custom half-edge structures and lock-free skip priority queues. — [PLOS ONE 2021, High-performance simplification of triangular surfaces using a GPU](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0255832)
- Jiao, Zuo and Persello (arXiv, Sept 7 2026) do topology-aware agricultural parcel polygonization. They rebuild "parcel faces from a common planar graph, enabling adjacent predicted parcels to reuse shared boundaries and avoid mutual interior intrusion." They report zero intrusion ratio and the highest shared-edge recall. They criticize prior methods that "rely on heuristic raster-to-vector conversion" or "reconstruct parcels independently." — [arXiv 2609.07520](https://arxiv.org/abs/2609.07520)

### Inferences
- A label map (multi-class or instance IDs) can be turned into a planar partition with the same 2×2 machinery:
  - Any 2×2 quad with two or more distinct labels marks a boundary. Quads with three or more labels are junction (degree-3/4) nodes.
  - Edges between junctions are arcs shared by exactly two faces, the TopoJSON model.
  - Simplifying each arc once (GPU DP, one arc per thread) and then assembling the faces gives gap-free and overlap-free coverage by construction.
  
  This is the GPU-native version of what the 2026 parcel paper does with heuristics. I found no existing GPU implementation, so it is a clear gap.
- For rectilinear output, simplification can be done during linking: drop collinear vertices, or apply stair-step removal at a tolerance, all in a single segmented scan.

### Gaps
- I found no GPU Visvalingam-Whyatt implementation and no GPU arc-sharing (TopoJSON-style) polygonizer.
- I found no benchmark numbers for GPU DP on large polygon sets newer than about 2015.
- GDAL Polygonize and rasterio.features.shapes are the de-facto CPU label-map polygonizers. I did not fetch their docs to confirm their algorithm or complexity this session.

---

## Q3. Polygon-native deep learning: how do these methods do mask→polygon, and is post-processing a bottleneck or error source?

### Takeaway
The field has split three ways:
- **Mask + vectorization post-process:** FFL, HiSup, most SAM pipelines.
- **Mask + learned vertices / attraction:** HiSup, PolyWorld-style vertex graph matching.
- **Fully end-to-end sequence or primitive prediction:** Pix2Poly (WACV 2025), P2PFormer, BoundaryFormer, PolyFormer, E2EC.

The end-to-end camp's main stated motive is that raster→vector post-processing is slow and error-prone. That is direct evidence the step matters. A fast, exact, topology-correct polygonizer would mainly help the first camp, and would also serve as an evaluation and supervision tool for the others.

### Cited Findings
- **Frame Field Learning (FFL)**, Girard et al., CVPR 2021, predicts a mask plus a frame field and then runs a polygonization post-process (ASM-style contour optimization guided by the frame field, then simplification). — [CVPR 2021 paper](https://openaccess.thecvf.com/content/CVPR2021/papers/Girard_Polygonal_Building_Extraction_by_Frame_Field_Learning_CVPR_2021_paper.pdf); [arXiv 2004.14875](https://arxiv.org/pdf/2004.14875)
- **HiSup** (ISPRS J. 2023) states that FFL-type methods "incur low efficiency during inference due to their reliance on post-processing methods to vectorize polygons and reduce redundant vertices, which heavily limits their practical application." HiSup's own post-process is simple: boundary pixels are attracted to the nearest predicted vertex (MaV-Attr). — [ScienceDirect HiSup](https://www.sciencedirect.com/science/article/pii/S0924271623000667); [Spyrosoft summary](https://spyro-soft.com/blog/geospatial/closing-the-geometry-gap-in-building-footprint-extraction-our-first-look-at-hisup)
- HiSup on CrowdAI: AP 79.4, which is +12.4 pp over FFL and +5.6 pp over the prior best joint vertex-edge model; AP_boundary 66.5 (+16.5 pp). — [Spyrosoft summary](https://spyro-soft.com/blog/geospatial/closing-the-geometry-gap-in-building-footprint-extraction-our-first-look-at-hisup) (secondary source; verify against the paper)
- **Pix2Poly** (Adimoolam, Poullis et al., WACV 2025) uses ViT patch embeddings and an autoregressive transformer decoder that emits vertex tokens. Connectivity comes from an optimal-matching network, giving ring graphs end to end. It claims to avoid "complicated, computationally intensive raster loss functions" and handles both building footprints and road networks. — [arXiv 2412.07899](https://arxiv.org/abs/2412.07899); [GitHub](https://github.com/yeshwanth95/pix2poly); [WACV PDF](https://openaccess.thecvf.com/content/WACV2025/papers/Adimoolam_Pix2Poly_A_Sequence_Prediction_Method_for_End-to-End_Polygonal_Building_Footprint_WACV_2025_paper.pdf)
- **P2PFormer** (2024) is a primitive-to-polygon method that "does not rely on any postprocessing." It beats PolyWorld by 2.7 AP and 6.5 AP75. — [arXiv 2406.02930](https://arxiv.org/pdf/2406.02930)
- **BoundaryFormer** (Lazarow, Xu, Tu, CVPR 2022) regresses polygon boundaries with a transformer. Supervision is mask-only, through a differentiable rasterizer. It matches or beats Mask R-CNN on COCO and Cityscapes and transfers better across datasets. — [CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Lazarow_Instance_Segmentation_With_Mask-Supervised_Polygonal_Boundary_Transformers_CVPR_2022_paper.pdf)
- **P³ dataset** (2025, "Pixels, Points and Polygons for Multimodal Building Vectorization") is a new benchmark for building vectorization. — [arXiv 2505.15379](https://arxiv.org/pdf/2505.15379)
- An "Edge Length Loss for Improved Building Polygon Extraction" (ACM, 2025) shows that polygon-geometry losses remain an active topic. — [ACM DL](https://dl.acm.org/doi/10.1145/3764919.3770877)
- **Agricultural parcels (2026):** semantic-guided diffusion plus topology-aware polygonization. Independent per-instance polygonization is the named failure mode. — [arXiv 2609.07520](https://arxiv.org/abs/2609.07520)

### Inferences
- **Where an exact GPU polygonizer unlocks things:**
  1. As a drop-in replacement for the CPU vectorization step in FFL/HiSup-style and SAM-mask pipelines. The stated efficiency complaint is about that step.
  2. As the ground-truth generator for sequence models like Pix2Poly. Training labels are rings, so fast, exact, canonical ring extraction from masks would allow on-the-fly augmentation (rotate/warp the mask, then re-polygonize per batch on GPU) instead of transforming precomputed polygons.
  3. As a metric engine. C-IoU, N-ratio, PoLiS, boundary AP and max-tangent-angle all need polygons, and evaluation currently polygonizes on CPU.
- "SAM to polygon" in practice (SAMGeo, Label Studio, CVAT plugins) is SAM mask + OpenCV findContours + approxPolyDP on CPU. I did not verify a specific tool's source this session. This is exactly the step a GPU polygonizer would replace for SAM-scale batch annotation.

### Gaps
- I did not retrieve exact numbers for PolyWorld, PolyR-CNN, GraphMapper, PolyMapper, E2EC, PolyFormer (CVPR 2023 referring segmentation as vertex sequences), PolyBuilding or Pix2Poly AP/C-IoU from primary sources this session. Note only their architecture class.
- I found no paper that measures the polygonization post-process wall-time share of total inference. HiSup asserts FFL's post-process is slow but gives no number in the material I could retrieve.

---

## Q4. Differentiable rasterization/vectorization and topology-aware losses: can an exact polygonizer be a differentiable or straight-through layer? How do bit-quads connect to Euler characteristic and Minkowski functionals?

### Takeaway
Bit-quad counting, the same 2×2 classification the user's CUDA kernel already does, gives area, perimeter and Euler characteristic (the three 2D Minkowski functionals) in one pass. 2025–2026 papers now make bit-quad Euler characteristic deep-learning-compatible and GPU-fast, as fast topology losses and violation maps. This is the strongest intellectual bridge: the user's classifier already emits topological invariants. A straight-through polygonizer (exact forward, soft-rasterizer backward) is plausible, but I found no published example.

### Cited Findings
- Liu Li, Qiang Ma, Cheng Ouyang, Johannes C. Paetzold, Daniel Rueckert and Bernhard Kainz (arXiv 2507.23763, IEEE TMI) give a fast Euler characteristic χ computation in 2D and 3D. It is described as a "deep learning-compatible pipeline for bit-quads-based 2D Euler characteristic calculation," with bit-octets in 3D. They use χ error as a metric, generate a "topological violation map," and refine with a correction network. It is presented as a fast alternative to persistent homology, which has polynomial complexity. — [arXiv 2507.23763](https://arxiv.org/abs/2507.23763); [ResearchGate](https://www.researchgate.net/publication/394175162_Topology_Optimization_in_Medical_Image_Segmentation_with_Fast_Euler_Characteristic)
- Quarenghi, Cotsakis and Beucler (arXiv 2604.11422, April 2026) propose the "Minkowski Image Loss," a differentiable surrogate for area, perimeter and Euler characteristic. It uses temperature-controlled sigmoids, continuous logical operators and spectrally-normalized Lipschitz CNNs, applied to precipitation super-resolution (EUMETNET OPERA). It eliminated geometric violations, but Lipschitz regularization over-smoothed gradients. — [arXiv 2604.11422](https://arxiv.org/abs/2604.11422)
- Udit Saxena (arXiv 2510.20271, NeurReps @ NeurIPS 2025) gives CUDA kernels for Euler Characteristic Curves with 128B-coalesced access and hierarchical shared-memory accumulation (Ampere). Claimed speedups are 16–2000× over prior GPU implementations on synthetic grids, plus a PyTorch differentiable layer through DECT-style sigmoid relaxation. — [arXiv 2510.20271](https://arxiv.org/abs/2510.20271)
- Minkowski functionals in 2D are area, half-perimeter and Euler characteristic, where χ = #components − #holes. — [arXiv 2604.11422](https://arxiv.org/pdf/2604.11422)
- clDice (Shit et al., CVPR 2021) is a soft-skeleton topology-preserving loss for tubular structures. — [CVPR 2021](https://openaccess.thecvf.com/content/CVPR2021/papers/Shit_clDice_-_A_Novel_Topology-Preserving_Loss_Function_for_Tubular_Structure_CVPR_2021_paper.pdf)
- BoundaryFormer already trains a polygon predictor through a differentiable rasterizer with mask supervision. This is the inverse direction (polygon → mask), and it shows that soft-rasterizer gradients are good enough to train polygon outputs. — [CVPR 2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Lazarow_Instance_Segmentation_With_Mask-Supervised_Polygonal_Boundary_Transformers_CVPR_2022_paper.pdf)

### Inferences
- **The bit-quad identity, which the user's classifier already computes:**
  - Let n1 = quads with one foreground pixel (convex corners), n3 = quads with three (concave corners), nD = diagonal quads.
  - Gray's formulas: χ₄ = (n1 − n3 + 2nD)/4 and χ₈ = (n1 − n3 − 2nD)/4.
  - Rectilinear vertex count = n1 + n3 + 2nD.
  - The geometric reason: each outer ring turns +4 net quarter-turns and each hole ring turns −4.
  
  So the same GPU pass gives χ, the vertex count, area (popcount) and exact pixel-edge perimeter (count of 0/1 transitions) at no extra cost. The vertex count can pre-size output buffers before linking. The formulas are classic (Gray 1971); I did not fetch a primary source for them this session.
- **A "topology certificate" for free (inference):**
  1. After linking, (#outer rings − #hole rings) must equal χ.
  2. #outer rings must equal the foreground CCL count.
  
  Checking both is a cheap correctness test for the GPU linker and a debugging oracle.
- **Straight-through polygonization layer (speculative design):**
  - Forward: an exact GPU polygonizer on the thresholded mask.
  - Backward: gradients via a soft rasterizer (DiffVG/BoundaryFormer-style) or via the sigmoid-relaxed bit-quad statistics (2604.11422, 2510.20271).
  
  This would let polygon-space losses (vertex count, right-angle regularity, χ) train a mask network. That kind of loss is otherwise unavailable without a learned polygon head.
- The soft bit-quad χ losses do soft counting only and never produce rings. An exact polygonizer adds ring identity, so a violation can be localized to a specific spurious hole or component rather than a pixel heat-map.

### Gaps
- I did not fetch DiffVG (Li et al., SIGGRAPH Asia 2020) or LIVE (Ma et al., CVPR 2022) primary pages this session. They are foundational differentiable vector-graphics rasterizers; cite them from their primary papers when writing.
- I found no published work that uses an exact mask→polygon converter as a straight-through estimator inside training.
- I did not research persistent-homology / Betti-matching losses (Hu et al. 2019; Stucki et al. Betti matching 2023) or their 2025–2026 speedups this session.

---

## Q5. Image/vector tracing state of the art (Potrace, VTracer, StarVector, OmniSVG, LLM SVG)

### Takeaway
Raster→SVG has moved to multimodal LLMs that emit SVG code: StarVector (CVPR 2025), OmniSVG (NeurIPS 2025) and VFIG (2026, diagrams). These optimize perceptual or semantic fidelity, not exactness or topology. Classical tracers such as Potrace and VTracer remain the exact-ish baselines. There is a clear niche for "exact, topology-guaranteed, GPU-speed" tracing that the LLM approaches do not target.

### Cited Findings
- StarVector (Rodriguez et al., CVPR 2025) is a multimodal LLM that generates SVG code. It is trained on SVG-Stack (2M samples) and introduces SVG-Bench (10 datasets; tasks are image vectorization, text-to-SVG and diagram generation). — [CVPR 2025](https://openaccess.thecvf.com/content/CVPR2025/html/Rodriguez_StarVector_Generating_Scalable_Vector_Graphics_Code_from_Images_and_Text_CVPR_2025_paper.html); [arXiv 2312.11556](https://arxiv.org/html/2312.11556)
- OmniSVG (Yang et al., NeurIPS 2025) is built on native VLMs. It introduces the MMSVG-2M dataset (2M assets: icons, illustrations, anime characters) and MMSVG-Bench for text-to-SVG and image-to-SVG. It reports that StarVector-8B fails on complex MMSVG-Character inputs. — [NeurIPS 2025 paper](https://papers.nips.cc/paper_files/paper/2025/file/a510f05a574d4203ef3952973672fe2f-Paper-Conference.pdf); [arXiv 2504.06263](https://arxiv.org/html/2504.06263v1)
- VFIG (arXiv 2603.24575, 2026) vectorizes complex diagrams to SVG with vision-language models. — [arXiv 2603.24575](https://arxiv.org/pdf/2603.24575)

### Inferences
- LLM tracers produce token-length-limited, approximate geometry with no guarantee of watertightness or hole correctness. An exact GPU polygonizer can be a pre-tokenizer: it compresses an image into exact rings that an LLM then simplifies or labels. It can also be a verifier that rasterizes the LLM output and compares it with exact rings.

### Gaps
- I did not fetch Potrace (Selinger 2003) or VTracer (visioncortex) primary docs this session. Potrace is the classic CPU curve-fitting tracer for bitmaps; VTracer is a Rust tracer that handles color images. Neither is GPU-accelerated as far as I know, but I did not verify this.
- I found no GPU implementation of Potrace-style optimal polygon / Bezier fitting.

---

## Q6. Compression: polygon/chain code vs RLE as mask representations; edge-to-cloud transmission

### Takeaway
COCO RLE is O(√n) for simple shapes, as are polygons and chain codes; all three scale with boundary length. Chain codes with entropy coding beat compressed raster formats for binary masks in the literature. A GPU polygonizer that emits rectilinear rings is effectively a vertex-only chain-code encoder: each vertex is one turn, so a ring is a sequence of run lengths alternating H/V. That makes it a natural fast mask codec for edge-to-cloud transmission.

### Cited Findings
- COCO RLE size is "proportional to the number of boundaries pixels." It is O(√n) for simple shapes, area/union/intersection can be computed directly on RLE, and LEB128-style variable-bitrate encoding adds compression. Ground-truth COCO masks are mostly stored as polygons and converted to RLE when needed. — [pycocotools mask.py](https://github.com/cocodataset/cocoapi/blob/master/PythonAPI/pycocotools/mask.py)
- Chain-code methods for segmentation contours achieved over 27% bits-per-pixel savings versus prior state of the art. "Even in uncompressed form, chain codes can be more space efficient than a compressed raster representation." — [IEEE 7952399, Efficient representation of segmentation contours using chain codes](https://ieeexplore.ieee.org/document/7952399/); [ResearchGate: compressed storage of large binary images using chain codes](https://www.researchgate.net/publication/383878169_Efficient_compressed_storage_and_fast_reconstruction_of_large_binary_images_using_chain_codes)
- Chain-code compression with modified interpolative coding. — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0045790618327666)
- Contour-based chain-code serialization for lossless compression of voxelized 3D objects (Entropy, 2026). — [doi:10.3390/e28070774](https://doi.org/10.3390/e28070774)
- TorchVision has an open issue about supporting compact/encoded RLE masks in its COCO wrapper, which shows friction in current mask formats. — [pytorch/vision #8351](https://github.com/pytorch/vision/issues/8351)

### Inferences
- A rectilinear ring stores only its corner vertices. Corners alternate H/V, so each vertex needs only one new coordinate delta. This is a "vertex chain code" with about one varint per corner, against roughly two RLE runs per scanline crossing. For blocky or Manhattan masks (buildings, IC layouts, document regions) the polygon form should be much smaller than RLE. For noisy masks the corner count explodes and RLE may win. I found no published benchmark of this trade-off.
- Because χ, area and vertex count come from the classification pass, an edge device can choose the codec (polygon vs RLE) per mask before linking, using the vertex count as an exact size predictor.

### Gaps
- I found no 2024–2026 paper that benchmarks polygon vs RLE vs chain code for ML mask transmission (edge-to-cloud) with GPU encoding times.
- I found no GPU chain-code encoder.

---

## Cross-cutting gaps (things that appear not to exist yet)
1. A fully GPU-resident polygonizer that outputs closed, oriented rings with an outer/hole hierarchy, open-source and library-grade (cuCIM/CV-CUDA/OpenCV CUDA level). The closest items are an Applied Sciences 2025 paper (code status unknown) and a 2026 triad tracker that only reaches CPU parity.
2. A GPU planar-partition (multi-label, shared-arc) polygonizer with topology-preserving simplification.
3. A YACCLAB-style benchmark for contour extraction and polygonization.
4. A published "sort + pair + list-rank" rectilinear ring-linking algorithm that combines CCL-derived hole parentage with a χ-based correctness certificate. This is the user's opportunity.
5. An exact polygonizer used as a straight-through layer in mask-network training.
6. A benchmark of polygon/chain-code vs RLE for GPU-side mask compression in edge-to-cloud ML.
