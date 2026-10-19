# Gap-Closing: Prior Art and Engineering Feasibility for a Fully GPU-Native Raster→Polygon Engine (and RL Reward Engine)

Research date: 2026-10-12. This closes the open gaps from `research_notes/GPU delineation use case map/r2_gpu_polygonizer_prior_art.md`, which ran out of search budget.

**Evidence labels** (same scheme as r2):
- **[SRC]**: I read source code, a header, or a raw repository file.
- **[DOC]**: I read the official documentation text.
- **[PAPER]**: I extracted and read the relevant part of the full text.
- **[META]**: abstract or bibliographic metadata only.
- **[RUN]**: I ran a local script on an Apple M3, single thread, CPU only.
- **[SNIPPET]**: search-snippet level only. Treat as unverified.

**Local artifacts** (scratchpad of this session; not committed):
- NPP 12.4 headers taken from NVIDIA's PyPI wheel `nvidia_npp_cu12-12.4.1.87` (win_amd64): `.../scratchpad/nppwheel/nppi_filtering_functions.h`, `nppdefs.h`.
- GEOS benchmarks: `.../scratchpad/scripts/bench_geos.py` and `bench_rect.py`. Both use a venv with Shapely 2.2.0, GEOS 3.14.1 and NumPy 2.5.4.

---

## Q1. NVIDIA NPP contour functions: how are holes, inner contours and nesting handled?

### Takeaway
NPP's contour API produces **one boundary-pixel list per connected-region label ID**.
- It has **no outer/hole flag, no parent/child link and no hierarchy field** in any of its data structures. I verified this in the header.
- Labeling covers *every* constant-value region, background included. Enclosed background (a hole) therefore gets its own label and its own contour, but nothing links it to the enclosing object.
- NVIDIA's own sample admits that "a few contour pixels can escape insertion into the ordered geometry list". NPP is therefore not a topology-exact ring extractor.
- The functions are documented as far back as the CUDA 11.4 docs.

### Cited Findings
**Labeling semantics: background regions get labels too [SRC]**
- `LabelMarkersUF` header text: "A connected region is any pixel region where all pixels in the region have the same pixel value."
- Connectivity: "nppiNormInf will use 8 way connectivity and nppiNormL1 will use 4 way connectivity."
- Limits and algorithm: ROI "up to 4 gigapixels"; based on "An Optimized Union-Find Algorithm for Connected Components Labeling Using GPUs" by Jun Chen et al.
- Source: `nppi_filtering_functions.h` (NPP 12.4.1.87, from [nvidia-npp-cu12 wheel on PyPI](https://pypi.org/project/nvidia-npp-cu12/)).

**Per-label info is a single record [SRC]**
- `NppiCompressedMarkerLabelsInfo` holds exactly: `nMarkerLabelPixelCount`, `nContourPixelCount`, `nContourPixelsFound`, `oContourFirstPixelLocation` and `oMarkerLabelBoundingBox`.
- That is one contour pixel count and one first pixel per label, with no ring count, hole count or parent. — `nppdefs.h` (same wheel)

**Contour pixel structs: still no hierarchy [SRC]**
- `NppiContourPixelGeometryInfo` holds:
  - ordered, previous, center and next pixel locations;
  - `nOrderIndex` / `nReverseOrderIndex`;
  - `nFirstIndex` / `nLastIndex`, documented as "index of first/last ordered contour pixel **in this subgroup**";
  - next/prev indices, three used/linked/output flags, and `nContourInteriorDirection`.
- `NppiContourPixelDirectionInfo.nMarkerLabelID` is documented as the "MarkerLabelID of contour interior connected region".
- `NppiContourBlockSegment` = {`nMarkerLabelID`, `nContourPixelCount`, `nContourStartingPixelOffset`, `nSegmentNum`} ("relative block segment number within this particular contour ID").
- No struct has any parent, hole, outer or nesting field. — `nppdefs.h`

**Geometry-list generation [SRC]**
- `nppiCompressedMarkerLabelsUFContoursGenerateGeometryLists_C1R_Ctx(...)` takes:
  - `bOutputInCounterclockwiseOrder`: "if nonzero then output geometry list for each contour in counterclockwise order, otherwise clockwise";
  - `nFirstContourGeometryListID` / `nLastContourGeometryListID`.
- Output goes to a **host** list. The header says: "1 channel connected region contours image to generate contours geometry info list in host memory."
- Header warnings, verbatim:
  - the "256K pixels … bypassed" default;
  - "contours approaching a million pixels can take many minutes";
  - "contour ID 0 can contain ALL contours in the image."
- Source: `nppi_filtering_functions.h`.

**NVIDIA sample `NPP+/findContour` [SRC]**
- It labels with `nppiNormInf` (8-way) and calls the generator with counterclockwise output.
- It loops `nID` over label IDs, one geometry list per ID.
- Its comment block adds: "**In general a few contour pixels can escape insertion into the ordered geometry list.** When geometry output is in clockwise order the extra pixels will be at the end of the geometry list otherwise they will be at the start of the list."
- It also shows a second call, `nppiCompressedMarkerLabelsUFContoursOutputGeometryLists_C1R`, that "can be recalled any number of times" after generation.
- Source: [findContourNPPPlus.cpp](https://github.com/NVIDIA/CUDALibrarySamples/blob/master/NPP%2B/findContour/findContourNPPPlus.cpp); [README](https://github.com/NVIDIA/CUDALibrarySamples/tree/master/NPP%2B/findContour).
- The README lists SM 7.0–8.0, CUDA 12.0+, and NPP+ as a separate download.

**Release history [DOC]**
- The group page `group__image__filter__compresssed__marker__labels__info.html`, including `GeometryLists`, `pContoursDirectionImage` and `ContoursImage`, exists in the CUDA **11.4.4, 11.5.2, 11.6.0, 11.6.2** and 11.8 archived NPP docs.
- The same URL returns 404 for 11.2.2 and 11.3.1.
- Sources: [NPP 11.6.2 archive](https://docs.nvidia.com/cuda/archive/11.6.2/npp/group__image__filter__compresssed__marker__labels__info.html); [NPP 11.8 archive](https://docs.nvidia.com/cuda/archive/11.8.0/npp/group__image__filter__compresssed__marker__labels__info.html).
- The 11.8 archive page shows the same `bOutputInCounterclockwiseOrder` parameter and the 256K note.

**Current docs (v13.x) [DOC]**
- The single-page `image_filtering_functions.html` I fetched documents `nppiContoursImageMarchingSquaresInterpolation_{32f,64f}` in full, with all the same caveats, including the new note that bounding boxes "may vary slightly between runs … due to the parallel reduction".
- It mentions `nppiCompressedMarkerLabelsUFInfo_32u_C1R_Ctx` and `nppiCompressedMarkerLabelsUFContoursPixelGeometryLists_C1R_Ctx` only inside parameter descriptions; I found no standalone entry for either.
- The function name differs from the 12.4 header (`...ContoursGenerateGeometryLists...`).
- Source: [NPP image filtering functions](https://docs.nvidia.com/cuda/npp/image_filtering_functions.html).

**Third-party clean-room reimplementation [SRC; not authoritative]**
- `CompilerFans/npp` (0 stars, pushed 2026-09-10) reimplements the geometry lists.
- Its contour pixel test is "has neighbor with different label".
- It traces each label once on the host, starting from `oContourFirstPixelLocation`, using a Moore-style 8-direction search.
- A trace that starts from one first pixel cannot reach a label's separate inner boundary.
- Source: [nppi_compressed_markerLabels.cu](https://github.com/CompilerFans/npp/blob/main/src/nppi/nppi_filtering_functions/nppi_compressed_markerLabels.cu).

### Inferences
- **How NPP represents a donut** (inferred from the semantics above; not run, since no GPU is available):
  - **Foreground ring label A.** Its "contour pixels" are A-pixels adjacent to a different label, which includes **both** its outer and inner boundary pixels. They sit under one contour ID with one first pixel and one count.
  - **Hole label H.** The enclosed background has its own label, and its contour is the ring of H-pixels touching A.
  - **"Contour ID 0 can contain ALL contours."** This is consistent with ID 0 being the unbounded background region, which touches every object.
  - So holes exist only **implicitly**, as separate background-label contours. Recovering nesting requires extra work, such as adjacency to the enclosing label plus an "is this region unbounded" test.
  - The undocumented "subgroup" fields (`nFirstIndex`/`nLastIndex`) and `NppiContourBlockSegment.nSegmentNum` suggest one ID's list can be split into several subgroups, which could hold the outer and inner chains. This is unconfirmed.
- **Differences from the proposed engine** (now verified rather than suspected):
  1. NPP outputs pixel-center 8-chains, not pixel-edge rectilinear rings.
  2. It has no hole hierarchy and no parent links.
  3. Its ordered output is host-side, and NVIDIA admits it may be incomplete ("a few contour pixels can escape").
  4. It skips contours over 256K pixels by default.
  5. It emits nothing in GeoArrow/WKB form.

  The engine's guarantees (exact rings, CCL parents, the χ certificate) are a real differentiator. NPP is still the baseline to benchmark against for the label → contour step.

### Gaps
- I did not run a donut test, because there is no CUDA GPU on this machine. NPP's actual output for (a) a donut and (b) nested islands-in-holes remains unverified. The single most valuable experiment: run `findContour` on a 3-level nested image, then dump the `nID` lists and their `nFirstIndex`/`nLastIndex` subgroups.
- Whether 11.4 is truly the first release is not certain. A 404 on 11.2/11.3 could reflect a different URL scheme rather than absence. I did not check the CUDA release notes.
- The v13 docs page lacks full entries for the Info/Generate functions. Whether they were renamed (`...PixelGeometryLists`), moved to NPP+ only, or simply dropped from the page is unresolved.

---

## Q2. Modern GPU list ranking, pointer jumping and Euler tour: what throughput should we expect for 10^8 nodes?

### Takeaway
I found **no list-ranking benchmark on A100/H100-class hardware from 2015–2026**.
- A June 2026 list-ranking paper's related-work section still cites only the 2009–2012 GPU papers for GPUs.
- The most recent GPU Euler-tour work (IPDPS 2021, GTX 980) uses **Wei–JáJá sparse ruling sets** because they "perform much better than the classical pointer jumping". It also reports that array scan is about 7–8× faster than list ranking.
- The best modern proxy I found is GPU CCL on an A100: BKE reaches 26.6 MPixel/ms.
- My estimate for 10^8 ring corners is a few tens of milliseconds for a sparse-ruling-set or active-set-compacted pointer-doubling pass on an A100/H100. This is an inference.

### Cited Findings
**Polak, Siwiec and Stobierski, "Euler Meets GPU: Practical Graph Algorithms with Theoretical Guarantees" (IPDPS 2021) [PAPER]**
- Hardware: experiments ran on an NVIDIA GeForce GTX 980.
- Method: "For list-ranking we implement the GPU-optimized Wei-JaJa algorithm [64], which performs much better than the classical pointer jumping technique."
- On cost: "on GPU the scan operation (i.e. array prefix sum) is much faster then list-ranking (see e.g. [64], who report 7-8x difference)."
- Design pattern: "call a single list-ranking on the list, and … create an array of (directed) edges in the Euler tour ordering … all the following prefix sum calculations … by using a fast scan primitive."
- They use moderngpu for sort and scan.
- Source: [arXiv 2103.15217](https://arxiv.org/abs/2103.15217).

**Sanders, Schimek, Uhl and Weidmann, "Engineering Scalable Distributed List Ranking" (arXiv 2606.09318, 2026) [PAPER]**
- Their related work says: "The problem has also been examined on GPUs [13, 23]." [13] is Rehman et al., ICS 2009; [23] is Wei & JáJá, Parallel Process. Lett. 2012.
- "The sparse ruling-set variants are up to 15× faster than pointer doubling."
- They rank "over 100 billion elements" in 6.7 s on 24,576 CPU cores (distributed CPU, not GPU).
- Source: [arXiv 2606.09318](https://arxiv.org/abs/2606.09318).

**Sahu et al., "Fully Dynamic Rooted Spanning Tree on GPU" (arXiv 2607.20211, July 2026) [PAPER]**
- Runs on an NVIDIA A100 80 GB.
- Uses Euler-tour ranking plus "parallel pointer jumping on parent[·]" inside dynamic spanning-forest updates.
- Reports update throughput (2M insertions/s, 1.4M deletions/s). It gives **no standalone list-ranking timings**.
- Source: [arXiv 2607.20211](https://arxiv.org/abs/2607.20211).

**Wei & JáJá (2010 tech report; already in r2) [PAPER]**
- 64M-node random list in about 297 ms on a Tesla C1060 (102 GB/s).
- Linear scaling from 1M to 256M nodes.
- Source: [UMIACS-TR-2010-08](http://drum.lib.umd.edu/bitstream/1903/10600/1/UMIACS-TR-2010-08.pdf).

**Rehman, Kothapalli and Narayanan (ICS 2009) [SNIPPET]**
- Snippet: "rank a random list of 32 million elements in about a second," about 8–9× over CPU, on the GPUs of the time.
- Source: [ResearchGate](https://www.researchgate.net/publication/221235489_Fast_and_Scalable_List_Ranking_on_the_GPU).

**Modern GPU CCL proxy: Bolelli, Allegretti, Lumetti and Grana, "A State-of-the-Art Review with Code about Connected Components Labeling on GPUs" (IEEE TPDS 2024) [PAPER]**
- Setup B is an **A100 80GB PCIe** with NVCC 12.0.
- On XDOCS (1,677 images, average 4,853 × 3,387, about 16.4 MP; "the hardest dataset"):
  - BKE: **0.621 ms/image, 26.55 MPixel/ms**;
  - BUF: 0.740 ms (22.28 MPixel/ms);
  - UF: 1.041 ms.
- BKE "is always the best-performing algorithm on 2D datasets" and is "the only GPU-based algorithm available in the OpenCV library."
- Source: [TPDS 2024 PDF](https://federicobolelli.it/media/publications/pdfs/2024tpds.pdf).

**Library coverage [SRC/search]**
- Gunrock, cuGraph and ParlayLib were not found to expose a GPU list-ranking primitive. My searches returned none. Treat this as "not found," not "verified absent."

### Inferences
- **Pointer doubling vs. sparse ruling set in the engine.** The literature consistently finds sparse ruling set (Wei–JáJá / Sibeyn / Helman–JáJá) beats Wyllie pointer doubling: up to 15× in distributed settings, and "much better" on GPU. But the engine's lists are **many disjoint cycles** (rings) whose lengths follow the object-size distribution, not one long list.
  - Recommended hybrid:
    1. Run one or two min-label doubling rounds for all corners.
    2. **Compact** to corners whose ring is still unresolved (rings of length ≤ 2^k finish after k rounds).
    3. Run a ruling-set pass only on the few long rings.
  - For building or SAM-mask workloads, most rings are short, so active-set compaction shrinks the work geometrically.
- **Throughput envelope for V = 10^8 corners** (inference; bandwidth B is a parameter):
  - Each doubling round does about 3 random gathers per active node, at about one 32 B sector each. That is about 10 GB of sector traffic per full round, roughly 5–10 ms per full round on an A100/H100-class part.
  - With active-set compaction and a ruling-set tail, a total of **~20–80 ms** is plausible.
  - Wei–JáJá scaled by bandwidth (297 ms × 102 GB/s ÷ ~1.5–3 TB/s) gives about 10–20 ms per 64M nodes, which is consistent.
  - CCL at 26.6 MPixel/ms implies about 0.4 s for 10^10 pixels, and about 4 ms for 10^8 pixels, if throughput holds at scale. The XDOCS images are small, so this extrapolation is uncertain.
- **Exploit locality.** Corners are emitted in row-major order. The horizontal successor is the adjacent index; only the vertical partner jumps. Gathers are therefore far more cache-friendly than in a random-list benchmark, so the random-list numbers are pessimistic.
- **Use the Euler-tour trick.** List-rank once, scatter ring vertices into a contiguous (ring_id, rank) array, then do all later per-ring work with segmented scans, which are 7–8× cheaper. Per-ring work includes shoelace area, the canonical start (argmin) and WKB sizes.

### Gaps
- There is no measured A100/H100 list-ranking or pointer-jumping number for 10^8 nodes. The estimates above must be replaced by a microbenchmark: a CUB-based doubling kernel versus a Wei–JáJá port on target hardware, using synthetic ring-length distributions.
- I did not open the Euler-Meets-GPU timing figures for raw list-ranking milliseconds; they are only plotted.
- I did not verify whether CCCL/CUB, cuGraph or Gunrock ship list ranking (search found none).

---

## Q3. Does any of Open3D, Halide, Taichi, Earth2Studio, GDAL (GPU), PyTorch3D, Kaolin or nvdiffrast have a GPU polygonizer? And which GPU polygon *rasterizers* could compute rewards when a VLM emits polygons?

### Takeaway
**None of the named libraries contains a raster→polygon (contour-with-holes) GPU implementation.**
- The 3D/graphics libraries (nvdiffrast, Kaolin, PyTorch3D) are **triangle** rasterizers.
- Taichi ships only a marching-squares *example*.
- Earth2Studio, Open3D and Halide had no contour or polygonize code.
- I found no GDAL GPU RFC.

For rewards in the reverse direction (polygon→mask):
- **nvdiffrast**, **Kaolin** and **PyTorch3D** can fill polygons only after **triangulation**. Triangulation is ill-defined for self-intersecting VLM output.
- **torch_contour** gives a dense, differentiable winding-number fill for single rings.
- Purpose-built CUDA polygon scanline rasterizers exist in GIS and lithography (up to 290× vs. CPU on H100 for Manhattan masks), but none is an open PyTorch drop-in.

### Cited Findings
**GitHub code search, repo-scoped [SRC]**
- Queries: `findContours`, `polygonize`, `find_contours`, `marching_squares`.
- Zero hits: Open3D (all four); Halide (all four); PyTorch3D (first three); Earth2Studio (`findContours`, `polygonize`).
- Taichi: only `python/taichi/examples/algorithm/marching_squares.py` and `examples/real_func/algorithm/marching_squares.py`.
- The rest were rate-limited, so I listed repository trees instead:
  - **Earth2Studio**: no file paths matching contour, polygon, raster or vectoriz.
  - **Open3D**: only `SelectionPolygon*` (a UI selection tool).
  - **Kaolin**: `kaolin/csrc/render/mesh/rasterization_cuda.cu` (mesh rasterizer).
  - **nvdiffrast**: `csrc/common/cudaraster/*` (triangle rasterizer).
  - **PyTorch3D**: `csrc/rasterize_meshes`, `rasterize_points`, `rasterize_coarse`.
- Sources: [Open3D](https://github.com/isl-org/Open3D), [earth2studio](https://github.com/NVIDIA/earth2studio), [taichi examples](https://github.com/taichi-dev/taichi/tree/master/python/taichi/examples/algorithm), [kaolin](https://github.com/NVIDIAGameWorks/kaolin), [nvdiffrast](https://github.com/NVlabs/nvdiffrast), [pytorch3d](https://github.com/facebookresearch/pytorch3d).

**nvdiffrast [DOC]**
- "The rasterization operation takes as inputs a tensor of vertex positions and a tensor of vertex index triplets that specify the triangles. Vertex positions are specified in clip space."
- Output is "a 4-channel float32 image with tuple (u, v, z/w, triangle_id) in each pixel."
- It has a `RasterizeCudaContext` (CUDA rasterizer, not only OpenGL), minibatch, instanced and range modes, and depth peeling.
- It has no polygon or fill-rule API. — [nvdiffrast docs](https://nvlabs.github.io/nvdiffrast/)

**torch_contour [SRC: README]**
- Contours are a tensor `(B, N, K, 2)` in [0,1], and outputs are `(B, N, H, H)`.
- The method is "based on the winding-number property … approximated with a smooth tanh so that gradients exist everywhere."
- It uses a fixed K per polygon, takes a single ring per (B, N) slot with no hole argument, and is a dense O(B·N·K·H²) computation.
- Source: [torch_contour](https://github.com/antoinehabis/torch_contour).

**Raster Join (Zacharatou et al., PVLDB 11(3), 2017) [PAPER, partial]**
- GPU (OpenGL) spatial aggregation over arbitrary polygons.
- Polygons are "first triangulated" (constrained Delaunay).
- The "accurate raster join" variant does point-in-polygon tests only for boundary pixels.
- Source: [VLDB PDF](http://www.vldb.org/pvldb/vol11/p352-zacharatou.pdf).

**Siemens Calibre blog (2026-04-13) [DOC, vendor blog]**
- A CUDA polygon rasterizer for computational lithography.
- It has a "five-stage parallel architecture" with "bounding box pruning, fine-grained thread allocation".
- It computes "fractional pixel coverage" via "the trapezoidal area formed by the edge and pixel boundary".
- Speedups on H100 vs. optimized CPU: "up to 290x" for Manhattan, "up to 45x" for curvilinear.
- Error: "less than 1% absolute error".
- Closed source. — [Siemens blog](https://blogs.sw.siemens.com/calibre/2026/04/13/why-a-gpu-rasterizer-matters-for-computational-lithography-in-both-performance-and-precision-at-scale/)

**GIS GPU polygon rasterization [SNIPPET]**
- "Speeding up large-scale geospatial polygon rasterization on GPGPUs" (ACM HPDGIS 2011): ">20×" over CPU for tens of thousands of polygons. — [ACM DL](https://dl.acm.org/doi/10.1145/2070770.2070772)
- "A parallel scheme for large-scale polygon rasterization on CUDA-enabled GPUs": about 1.3×10^8 polygons, 36.91× speedup. — [MetaTOC](https://metatoc.com/papers/74662-a-parallel-scheme-for-large-scale-polygon-rasterization-on-cuda-enabled-gpus) (page timed out on fetch; snippet only)

**OxiGDAL GPU [DOC]**
- `oxigdal-gpu` 0.1.7 (2026-09-30) is a wgpu backend. It offers element-wise ops, statistics, resampling, convolution, band math, reprojection and "tiled processing".
- It has no polygonize, raster-to-vector or contour operation. — [docs.rs/oxigdal-gpu](https://docs.rs/oxigdal-gpu)

**GDAL [search]**
- No GDAL RFC or branch for GPU polygonize was found. The current docs describe the CPU `gdal_polygonize` / `gdal raster polygonize`. — [gdal raster polygonize](https://gdal.org/en/stable/programs/gdal_raster_polygonize.html)

### Inferences
- **Reward rasterizer options for VLM-emitted polygons, ranked for RL:**
  1. **A custom CUDA even-odd/nonzero scanline or tile kernel.** One block per polygon bbox tile; each thread computes crossing parity or winding for its pixel over the polygon's edges, with edges staged in shared memory.
     - It handles **self-intersecting** and holed polygons directly via the fill rule, with no triangulation or `make_valid`.
     - Rough cost: work ≈ Σ (bbox pixels × edges). For 10^4 polygons × 64 edges × 32×32 bbox, that is about 6.5×10^8 edge tests, which should be ≲ 1 ms on an A100. This is an estimate, not measured.
  2. **nvdiffrast range mode with `triangle_id`.** It is fast and differentiable via antialias. It requires GPU-side triangulation (holes are fine; self-intersections are not), and overlapping polygons occlude each other unless rendered as separate range images.
  3. **torch_contour.** Differentiable, but dense and memory-heavy, and holes are not supported.
- **Pair the reward rasterizer with the polygonizer's own primitives.** Bit-quad χ, area and perimeter of the rasterized prediction give β0/β1 or χ rewards without polygonizing back. Corner-pair validity checks (Q6) give a validity reward without GEOS.

### Gaps
- I could not run any GPU rasterizer here (no CUDA), so I have **no measured runtime** for thousands of small polygons in nvdiffrast, Kaolin or torch_contour.
- Code search was rate-limited partway. For Kaolin, nvdiffrast, GDAL and cuCIM I relied on file-tree names, not content search. A function inside a generically named file could have been missed.
- The Siemens rasterizer, the 2011 HPDGIS paper and the CUDA-scheme paper were not read in full.

---

## Q4. GitHub search: "cuda polygonize", "gpu findContours", "cuda contour tracing holes", "gpu raster to vector"

### Takeaway
Authenticated GitHub **code** search, which closes r2's gap, found **no open-source GPU raster→polygon tracer with holes**.
- CUDA hits for `findContours` are files in application repositories (projects such as parking-spot detection, pose estimation and saliency maps). Judging by repo and file names they are not contour libraries, but I opened only a few of them.
- `polygonize` in CUDA means 3D marching cubes.
- The only CUDA "contour + holes" hit is a marching-squares *line renderer*.
- The NPP contour API is used in NVIDIA's sample and wrapped in C# and Rust bindings.

### Cited Findings
**`findContours language:Cuda` (30 hits) [SRC: result list]**
- Repos include qis/horus, Roboy/poseestimator, Maxerience-Electrosystems/CUDA, Parking-Spot-Detection, Cuda-FT-Saliency and compsygroup/3DI.
- None of these repos is a GPU contour library. — [GitHub code search](https://github.com/search?q=findContours+language%3ACuda&type=code)

**`polygonize language:Cuda` (9 hits) [SRC: result list]**
- Marching cubes (maya-implicit-skinning, SparseSurfaceConstraints, VoxelMapping, PRayGround, implicit_skinning) and `gauenk/bist/border_edges.cu`.
- BIST is "Bayesian-Inspired Space-Time Superpixels"; its kernel marks superpixel border edges in a point-cloud graph. — [bist border_edges.cu](https://github.com/gauenk/bist/blob/main/bist/csrc/border_edges.cu)

**`marching squares contour holes language:Cuda` (1 hit) [SRC]**
- `FahrenheitResearch/rustwx/.../contour_lines.cu` is a "GPU port of `rustwx_render::render::draw_contours` — the marching-squares" renderer. It draws anti-aliased isoline segments per cell. "Holes" there means invalid-data cells, and there is no segment linking.
- Source: [contour_lines.cu](https://github.com/FahrenheitResearch/rustwx/blob/main/rustwx-production-cuda/kernels/render/contour_lines.cu)

**`nppiCompressedMarkerLabelsUFContoursGenerateGeometryLists` (11 hits) [SRC: result list]**
- NVIDIA/CUDALibrarySamples `NPP+/findContour`;
- kunzmi/managedCuda (C#);
- vansweej/npp-rs, tinrab/singe (Rust);
- of78ojax/NPPContourTest;
- CompilerFans/npp (open reimplementation).
- Source: [GitHub code search](https://github.com/search?q=nppiCompressedMarkerLabelsUFContoursGenerateGeometryLists&type=code)

**Web search `site:github.com cuda polygonize raster holes` [search]**
- It returned only xarray-spatial polygonize (hybrid; CPU tracing) and GDAL/rasterio issues. — [xarray-spatial polygonize source](https://xarray-spatial.readthedocs.io/en/latest/_modules/xrspatial/polygonize.html)

### Inferences
- The ecosystem gap is confirmed at the code level. Beyond NPP's partial API, no open GPU implementation produces ordered rings with holes.
- The only GPU work on ordered contour extraction remains the 2022/2026 J. Supercomputing border-tracking papers (r2 Q2), which have no public code found.

### Gaps
- GitHub code search covers only indexed default branches and caps results. Private, non-default-branch or GitLab/Gitee code is not covered.
- "gpu raster to vector" and "cuda contour tracing" were run as repository searches in r2, not code searches. I did not re-run them as code searches because of the rate limit.

---

## Q5. Out-of-core, streaming and distributed polygonization: how do existing systems stitch seams?

### Takeaway
I found **no published GPU or out-of-core polygonization paper with exact seam stitching**. Practice falls into four patterns:
1. **Scanline streaming with no tiles** (GDAL two-arm chains; r2).
2. **Unit-edge cancellation across chunks** (xarray-spatial dask; r2).
3. **Strip rectangles + value dissolve** (R `vectra`): exact but heavy.
4. **Accept split polygons and reconcile later.** The 2026 global WSF settlement study did this with PostGIS on 50k × 50k tiles and a **non-parallelizable** cross-tile merge.

Port-pairing via sort (r2 Q5) remains the right GPU design. Its novelty claim at the out-of-core level is unrefuted.

### Cited Findings
**Schorcht et al., "Settlement percolation: global maps of Critical Distances" (arXiv 2603.04439v2, Sept 2026) [PAPER]**
- Input: global WSF 2019 = 5,138 GeoTIFF tiles of 22,488 × 22,487. These were "merged into a seamless global raster" to remove overlaps, then "partitioned into tiles of 50,000 × 50,000 pixels."
- Vectorization: "converted into vector geometries using the PostGIS function ST_DumpAsPolygons … directly adjacent raster pixels connected through a four-neighborhood were merged into single polygons, reducing approximately 8.8 billion settlement pixels to approximately 524 million polygon features … approximately 300 GB."
- Parallelism: "Python's multiprocessing pool … tiled partitions."
- Seams: clustering runs per tile, then "clusters located on adjacent tile boundaries were combined whenever their buffered geometries intersected, as determined by ST_Intersects(geometry). This merging step, which cannot be parallelized…"
- Software: PostgreSQL 15 / PostGIS 3.4.0 / GDAL 3.2.2.
- Source: [arXiv 2603.04439](https://arxiv.org/abs/2603.04439)

**R `vectra::polygonize` [DOC]**
- Reads "one tile-row strip at a time"; "within each strip equal-valued cells along a row collapse to a rectangle"; with `dissolve = TRUE` rectangles are "merged by value through `spatial_dissolve()`" using a sort/partition tier.
- No hole or seam semantics are documented. — [rdrr.io vectra polygonize](https://rdrr.io/cran/vectra/man/polygonize.html)

**Jiao, Zuo and Persello (arXiv 2609.07520, 2026) [SNIPPET]**
- They stitch "intermediate primitives in tile coordinates," then run "topology-aware polygonization once on the full tile," instead of merging independently polygonized patches. — [arXiv 2609.07520](https://arxiv.org/pdf/2609.07520)

**Google Earth Engine `reduceToVectors` [DOC, snippet-level]**
- Exposes `tileScale` ("scaling factor used to reduce aggregation tile size") and `maxPixels`. Seam handling across internal tiles is not documented. — [ee.Image.reduceToVectors](https://developers.google.com/earth-engine/apidocs/ee-image-reducetovectors)

**Carried over from r2 [SRC]**
- xarray-spatial: unit-edge cancellation plus union-find grouping of boundary polygons.
- GDAL: two-line scanline streaming.
- Sources: [xrspatial/polygonize.py](https://github.com/makepath/xarray-spatial/blob/master/xrspatial/polygonize.py); [gdal polygonize.cpp](https://github.com/OSGeo/gdal/blob/master/alg/polygonize.cpp)

**OxiGDAL GPU [DOC]**
- Has "tiled processing" for datasets exceeding VRAM, but no polygonize. — [docs.rs/oxigdal-gpu](https://docs.rs/oxigdal-gpu)

### Inferences
- **The WSF case is a concrete, citable motivating workload:**
  - 8.8 × 10^9 foreground pixels → 5.24 × 10^8 polygons, i.e. about 17 pixels per polygon on average;
  - CPU/PostGIS per-tile vectorization;
  - a serial seam merge.

  A GPU engine with sort-based port pairing would handle this as about 30 tiles of 50k² (2.5 × 10^9 px each), plus one global sort of seam ports. That replaces the "cannot be parallelized" step with a parallel sort.
- **Comparing the stitch strategies:**
  - Edge cancellation and dissolve are O(perimeter) or O(area) at seams.
  - Port pairing is O(#seam crossings).
  - Per-tile "accept split polygons" approaches (WSF, Earth Engine-style) give wrong β0 and polygon counts near seams. That matters directly for an RL reward, which must be seam-invariant.

### Gaps
- I did not find an academic paper titled along the lines of "streaming polygonization" or "out-of-core raster vectorization" with exact seams. The search was moderate (two web queries plus targeted fetches). PostGIS `ST_DumpAsPolygons` and GRASS `r.to.vect` internals were not read.
- WSF compute time for the vectorization step was not reported in the parts I read; only clustering times (4–28 h) are given.

---

## Q6. Shapely/GEOS validity checking cost at scale (is_valid, make_valid): baselines for reward computation

### Takeaway
I found no published polygons-per-second benchmark, so I measured one. Single-thread, Apple M3, Shapely 2.2.0 / GEOS 3.14.1:
- **`is_valid` on valid rectilinear polygonizer output costs about 150 ns per vertex**, scaling linearly. That is about 5.5 µs for a median 17-vertex polygon, or about 6.5M vertices/s per core.
- `is_valid` is cheap even on garbage (2–25 µs), because it exits early.
- **`make_valid` on invalid, self-intersecting polygons (typical noisy VLM output) is catastrophic and heavy-tailed.** With the default linework method:
  - 32 vertices: 9.3 ms;
  - 64 vertices: 71 ms;
  - 128 vertices: 687 ms per polygon.
- A GEOS IoU reward (make_valid + intersection + union) costs 0.09 / 0.88 / 14.9 ms per pair at 16 / 64 / 256 vertices.

So CPU GEOS is fine for *checking* polygonizer output, but unusable as an RL reward on raw VLM polygons at batch scale.

### Cited Findings
**Benchmark 1 [RUN]** (`bench_geos.py`; synthetic shapes; times are per polygon)

| Case | is_valid | make_valid (linework) | make_valid (structure) | area |
|---|---|---|---|---|
| star simple ring n=8 | 3.99 µs | 2.59 µs | 6.47 µs | 0.08 µs |
| star n=32 | 3.25 µs | 3.29 µs | 11.56 µs | 0.11 µs |
| star n=128 | 20.0 µs | 20.7 µs | 45.2 µs | 0.26 µs |
| star n=512 | 156 µs | 153 µs | 312 µs | 0.74 µs |
| star n=2048 | 1,739 µs | 1,536 µs | 3,110 µs | 2.27 µs |
| polygon + 1 hole, 32+32 | 8.2 µs | 8.1 µs | 50.1 µs | 0.45 µs |
| polygon + 1 hole, 512+512 | 295 µs | 301 µs | 1,221 µs | 3.39 µs |
| random self-intersecting n=8 (98% invalid) | 2.4 µs | 187 µs | 27 µs | 0.26 µs |
| random self-intersecting n=32 (100% invalid) | 5.9 µs | 9,297 µs | 517 µs | 1.8 µs |
| random self-intersecting n=64 | 12.2 µs | 71,153 µs | 2,035 µs | 4.2 µs |
| random self-intersecting n=128 | 24.9 µs | 686,892 µs | 9,163 µs | 10.8 µs |

- **IoU reward** (`make_valid` → `intersection` + `union` areas, against a star GT, prediction = GT vertices + N(0, 4 px)): n=16: 93 µs/pair; n=64: 877 µs/pair; n=256: 14.9 ms/pair.
- Note: the "rectilinear stair" generator in this script produced self-overlapping (invalid) rings by construction. Those rows (make_valid linework up to 818 ms at about 4,096 vertices) show invalid-rectilinear repair cost, not polygonizer output.

**Benchmark 2 [RUN]** (`bench_rect.py`: valid rectilinear polygons from thresholded Gaussian-noise blob masks, built as union of pixel squares, collinear vertices removed, with holes)

| Mask | Polygons | Median / max vertices | Holes | is_valid per polygon | ns/vertex | make_valid |
|---|---|---|---|---|---|---|
| 512², σ=2 | 885 | 17 / 721 | 73 | 5.5 µs | 143 | 5.3 µs |
| 512², σ=4 | 208 | 37 / 1,266 | 26 | 12.2 µs | 147 | 12.2 µs |
| 1024², σ=8 | 169 | 73 / 3,572 | 25 | 30.9 µs | 152 | 31.2 µs |
| 1024², σ=16 | 46 | 119 / 4,470 | 3 | 55.7 µs | 157 | 55.0 µs |

- All 100% valid.
- The largest polygons (about 3,700–4,500 vertices) took 0.56–0.74 ms each.

**Shapely docs [DOC, snippet-level]**
- `make_valid` requires GEOS ≥ 3.8.
- The "linework" method nodes all rings and rebuilds polygons; the "structure" method makes rings valid, then merges shells and subtracts holes.
- Source: [shapely.make_valid](https://shapely.readthedocs.io/en/2.1.1/reference/shapely.make_valid.html)

**Published benchmarks [search]**
- None found for is_valid/make_valid throughput. Searches returned only API docs and an unrelated project issue.

### Inferences
- **Reward-computation budget example** (inference from [RUN] numbers): an RL step with 4,096 rollouts × 20 polygons × about 64 vertices.
  - `is_valid` only: about 82k × 12 µs ≈ **1 s single-core** (about 30 ms on 32 cores).
  - GEOS IoU with `make_valid` on noisy outputs: about 82k × 0.9 ms ≈ **74 s single-core**, with a heavy tail (one 128-vertex garbage polygon costs about 0.7 s with linework).
  - Without a GPU reward engine, the reward step would dominate the RL step.
- **Policy recommendation for the reward engine:**
  1. Do not call `make_valid` in the reward loop. Score invalidity directly, e.g. reward = 0 or a penalty, using a GPU self-intersection / validity check.
  2. Compute IoU by **GPU fill-rule rasterization** (Q3), which is well-defined for invalid polygons.
  3. Compute β0/β1 from bit-quads / CCL on the rasterized mask.
  4. If GEOS-compatible repair semantics are needed, `method='structure'` is 10–75× cheaper than the default `linework` on self-intersecting input.
- **GPU validity-check design** (inference):
  - Rectilinear polygonizer output is valid by construction, and the χ / #ring certificate verifies it.
  - For VLM polygons, validity is (a) no proper segment intersections and (b) correct ring nesting. (a) is a batched segment-pair test over each polygon's edges: O(n²) per polygon, trivially parallel for n ≤ 256, or a sort-based sweep for larger n. (b) is point-in-polygon per hole.
  - This is a GPU analogue of GEOS `IsValidOp`. I found no GPU library that implements OGC validity (cuSpatial was archived; SedonaDB libgpuspatial does joins only, per r2).

### Gaps
- The measurements are single-threaded on an Apple M3 laptop CPU. Server x86 per-core numbers may differ by about ±30%, and multi-process scaling was not measured.
- The synthetic "noisy VLM" polygons (Gaussian vertex jitter, uniform random vertices) may over- or under-state real VLM failure modes. Real model outputs should be benchmarked.
- `is_valid` on non-rectilinear star polygons scaled superlinearly (128→2,048 vertices: 20 µs→1.7 ms, about 87× for 16× vertices). I did not investigate the cause; it may be GEOS's ring self-intersection / monotone-chain behaviour on highly non-monotone rings.
- No GPU validity implementation was found to compare against.
