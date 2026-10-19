# Round 2 Deep Dive: Prior Art and Feasibility of a Fully GPU-Native, Topology-Exact Raster→Polygon Engine

Research date: 2026-10-07. Builds on round-1 notes (`algorithms_ml.md`).

**Evidence labels.** Each finding is tagged with how it was established:
- **[SRC]**: I read the source code or the raw repository file.
- **[DOC]**: I read the official documentation text.
- **[PAPER]**: I extracted the full text of the paper.
- **[META]**: abstract or bibliographic metadata only (Semantic Scholar / arXiv API).
- **[RUN]**: verified by executing a local script (pure-Python CPU reference; no GPU is available on this Mac).
- **[SNIPPET]**: search-snippet level only.

**Process limits.** The session's WebSearch budget ran out partway through. Everything after that was done by fetching known primary URLs directly: GitHub raw/API, NVIDIA docs, arXiv, DOI metadata and open-access PDFs. Some breadth items could therefore not be checked; they are flagged under Gaps.

Local verification scripts:
- `/private/tmp/claude-501/r2scripts/verify_bitquad3.py`: bit-quad identities.
- `/private/tmp/claude-501/r2scripts/proto_sortpair.py`: sort-pair ring-linking prototype.

---

## Q1. Does any library or open-source project already provide GPU contour/polygon extraction with rings, holes and nesting?

### Takeaway
No GPU-resident library I checked outputs closed rectilinear rings with a hole hierarchy. Round 1 missed one important item: **NVIDIA NPP does ship a GPU contour API.** The chain is `nppiLabelMarkersUF` → `nppiCompressMarkerLabelsUF` → `nppiCompressedMarkerLabelsUFInfo` → `nppiCompressedMarkerLabelsUFContoursPixelGeometryLists`, with optional marching-squares interpolation. Its limits:
- It emits **ordered boundary-pixel chains per label**, not pixel-edge polygons.
- It needs host-side buffers.
- By default it skips contours larger than 256K pixels. NVIDIA's own docs say million-pixel contours "can take many minutes."
- Its hole/hierarchy semantics are undocumented.

Every other ecosystem tool either:
- has CCL only (CV-CUDA, ArrayFire, Kornia, cuCIM, VPI has neither), or
- runs CCL on GPU and traces on CPU (xarray-spatial), or
- is CPU-only (GDAL, and SedonaDB, which calls GDAL).

### Cited Findings
**NVIDIA NPP (v13.4 docs) [DOC]**
- **API surface.** NPP lists:
  - `nppiLabelMarkersUF_{8u32u,16u32u,32u}_C1R_Ctx` and the batched `nppiLabelMarkersUFBatch_*`;
  - `nppiCompressMarkerLabelsUF_32u_C1IR_Ctx`;
  - `nppiCompressedMarkerLabelsUFInfo_32u_C1R_Ctx`, `nppiCompressedMarkerLabelsUFGetGeometryListsSize` and `nppiCompressedMarkerLabelsUFContoursPixelGeometryLists_C1R_Ctx`;
  - `nppiContoursImageMarchingSquaresInterpolation_{32f,64f}_C1R_Ctx`;
  - `nppiSegmentWatershed_*`.
  
  — [NPP Image Filtering Functions](https://docs.nvidia.com/cuda/npp/image_filtering_functions.html)
- **LabelMarkersUF algorithm.** It is "based on the one described in 'An Optimized Union-Find Algorithm for Connected Components Labeling Using GPUs' by Jun Chen and others." `eNorm` selects 8-connectivity (`nppiNormInf`) or 4-connectivity (`nppiNormL1`). Batch versions support ROIs "up to 4 gigapixels." — [NPP docs](https://docs.nvidia.com/cuda/npp/image_filtering_functions.html)
- **Deprecated predecessor.** The older "Bound Segments" functions are deprecated: "Use the nppiLabelMarkerUF, nppiCompressMarkerLabelsUF, and nppiCompressedMarkerLabelsUFInfo functions to generate connected pixel region boundaries (contours)." — [NPP docs](https://docs.nvidia.com/cuda/npp/image_filtering_functions.html)
- **Documented limitations of the contour geometry-list output:**
  - "by default a contour that contains more than 256K pixels will be bypassed when generating the output geometry list."
  - It can be forced per contour ID, but "doing so for contours approaching a million pixels can take many minutes."
  - "contour ID 0 can contain ALL contours in the image."
  - Bounding boxes "may vary slightly between runs on the same input due to the parallel reduction."
  - Parameters include host pointers (`pContoursPixelGeometryListsHost`, `pContoursPixelsFoundListHost`, `pContoursPixelsStartingOffsetHost`), so the pipeline is not purely device-resident.
  
  — [NPP docs](https://docs.nvidia.com/cuda/npp/image_filtering_functions.html)
- **Data structures:**
  - `NppiContourPixelGeometryInfo` "provides contour (boundary) geometry info of uniquely labeled pixel regions ... in host memory in counterclockwise order relative to contour interiors." Its fields include prev/center/next pixel locations, `nOrderIndex`/`nReverseOrderIndex` and `nContourInteriorDirection`, with 8 compass directions.
  - `NppiCompressedMarkerLabelsInfo` holds, per label, a single `nContourPixelCount`, `oContourFirstPixelLocation` and `oMarkerLabelBoundingBox`.
  - `NppiContourTotalsInfo` holds `nTotalImagePixelContourCount` and `nLongestImageContourPixelCount`.
  
  — [NPP nppdefs](https://docs.nvidia.com/cuda/npp/nppdefs.html)

**cuCIM (RAPIDS) [SRC]**
- `cucim.skimage.measure` exports `label`, `regionprops`, `perimeter`, `perimeter_crofton`, `euler_number`, `approximate_polygon`, `subdivide_polygon` and moments. There is **no `find_contours`** and the module directory has no contour file. — [cuCIM measure/__init__.pyi](https://github.com/rapidsai/cucim/blob/main/python/cucim/src/cucim/skimage/measure/__init__.pyi)
- cuCIM `approximate_polygon` (Douglas-Peucker) is a host-driven loop. It calls `cp.asnumpy(coords[start, :])` per segment, so it is GPU arrays with serial control flow, not a parallel DP. — [cuCIM _polygon.py](https://github.com/rapidsai/cucim/blob/main/python/cucim/src/cucim/skimage/measure/_polygon.py)

**CV-CUDA [SRC]**
- In the public operator headers, the only CCL-related operator is `OpLabel`. There is no contour, findContours or polygon operator. — [CV-CUDA include/cvcuda](https://github.com/CVCUDA/CV-CUDA/tree/main/src/cvcuda/include/cvcuda)

**VPI 4.1 [DOC]**
- The algorithm table lists filters, pyramids, morphology, FFT, remap, stereo, KLT, Harris, optical flow, histograms, background subtraction, Canny, ORB, template matching, trackers and similar. It has **no CCL and no contour extraction**. — [VPI Algorithms](https://docs.nvidia.com/vpi/algorithms.html)

**Kornia, TorchVision, torch_contour, ArrayFire [SRC/DOC]**
- **Kornia** `contrib` exports `connected_components` and `connected_components_union_find`, with no contour function. — [kornia/contrib/__init__.py](https://github.com/kornia/kornia/blob/main/kornia/contrib/__init__.py)
- **TorchVision** ops: `masks_to_boxes`, NMS, box ops, RoI ops and losses. No contours. — [torchvision/ops/__init__.py](https://github.com/pytorch/vision/blob/main/torchvision/ops/__init__.py)
- **torch_contour** (v1.4.5) works in the opposite direction: differentiable contour → mask/distance-map/SDF. — [PyPI torch-contour](https://pypi.org/project/torch-contour/)
- **ArrayFire** `image.h` exposes `regions()` (CCL), `moments()` and `confidenceCC()`. No contour API. — [arrayfire image.h](https://github.com/arrayfire/arrayfire/blob/master/include/af/image.h)

**xarray-spatial [SRC]**
- **Algorithm.** `polygonize.py` does CCL, then boundary following ("identifies where polygon exteriors and holes start, and follows each of these boundaries"). Exteriors are anticlockwise, holes clockwise, and holes are grouped with their enclosing exterior. — [xrspatial/polygonize.py](https://github.com/makepath/xarray-spatial/blob/master/xrspatial/polygonize.py)
- **CuPy path is hybrid.** `_polygonize_cupy` is documented as "Hybrid GPU/CPU: GPU CCL for integer data, CPU CCL for float data, CPU boundary tracing in either case."
  - It runs `cupyx.scipy.ndimage.label` per unique value.
  - It then does `cupy.asnumpy(regions_gpu)` and traces on CPU (`_scan`).
  
  — [xrspatial/polygonize.py](https://github.com/makepath/xarray-spatial/blob/master/xrspatial/polygonize.py)
- **Contours.** `contour.py` uses marching squares. `_contours_cupy` copies to host (`cupy.asnumpy(data)`) before processing. — [xrspatial/contour.py](https://github.com/makepath/xarray-spatial/blob/master/xrspatial/contour.py)

**GDAL [SRC]**
- **Design.**
  - `GDALPolygonize` does two scanline passes. The first is a `GDALRasterPolygonEnumerator` pass with `ProcessLine`, which builds CCL equivalences. The second pass feeds a `Polygonizer` using `TwoArm` arrays sized `nXSize + 2` for the last and current lines.
  - `polygonize_polygonizer.h` states it "Implements The Two-Arm Chains EdgeTracing Algorithm" (kikitte.lee, 2023).
  - Options: `8CONNECTED=8`, otherwise 4-connectedness.
  
  — [gdal/alg/polygonize.cpp](https://github.com/OSGeo/gdal/blob/master/alg/polygonize.cpp); [polygonize_polygonizer.h](https://github.com/OSGeo/gdal/blob/master/alg/polygonize_polygonizer.h)
- **Stated behaviour.** GDAL says the algorithm "attempts to minimize memory use so that very large rasters can be processed" but produces "very dense polygon geometries, with edges that follow exactly on pixel boundaries." It is CPU-only. — [polygonize.cpp](https://github.com/OSGeo/gdal/blob/master/alg/polygonize.cpp)

**Apache SedonaDB [SRC/DOC]**
- **GPU library.** It contains `libgpuspatial`, a "GPU-accelerated Geospatial Processing Library." It currently supports spatial joins (DE-9IM) on WKB ArrowArrays (Point through MultiPolygon). — [sedona-db libgpuspatial README](https://github.com/apache/sedona-db/tree/main/c/sedona-libgpuspatial/libgpuspatial)
- **Polygonize is CPU.** SedonaDB's raster polygonize calls `GDALPolygonize`/`GDALFPolygonize`. — [sedona-gdal polygonize.rs](https://github.com/apache/sedona-db/blob/main/c/sedona-gdal/src/raster/polygonize.rs)

**RAPIDS cuSpatial [SRC]**
- The GitHub repo is **archived** (last push 2025-07-28). Its README says it "represents data in GeoArrow format." — [rapidsai/cuspatial](https://github.com/rapidsai/cuspatial)

**GitHub repository search [SRC]**
- I queried the GitHub repository-search API for "cuda polygonize", "gpu polygonize", "gpu findContours", "cuda contour tracing", "gpu marching squares contour", "cuda vectorize raster" and "gpu raster to vector". There were no raster→polygon GPU polygonizers among the results, which were voxelizers, polygon clippers (GH-CUDA), a QGIS SAM plugin, and NVIDIA/cuspa for transcript-to-cell-polygon assignment. "gpu findContours" and "cuda contour tracing" returned zero repositories. — [GitHub search API](https://api.github.com/search/repositories?q=cuda+polygonize)

### Inferences
- **NPP is now the closest "productized" prior art. It differs from the user's engine in four ways:**
  1. It outputs 8-direction boundary-**pixel** chains (pixel centers), not exact pixel-edge rectilinear rings.
  2. Its ordering stage involves host buffers, and per its own docs it scales badly for long contours (256K default cap; "many minutes" near 1M pixels). That points to a serial or poorly parallel ordering step.
  3. Per-label info holds a single first-pixel location, so outer/hole separation and nesting are not exposed as a hierarchy. This is unverified (see Gaps).
  4. It has no GeoArrow/WKB output.
  
  A fully device-resident, log-depth, exact rectilinear-ring engine with an explicit hole hierarchy remains novel relative to NPP. Any novelty claim must still cite NPP and be benchmarked against it.
- **The common ecosystem pattern is "GPU CCL + CPU trace."** It shows up in xarray-spatial and round-1 OpenCV. GPU CCL is everywhere (NPP, CV-CUDA, cuCIM, Kornia, ArrayFire, OpenCV-BKE). Ordered-ring extraction is the missing layer.

### Gaps
- **NPP hole/nesting semantics.** I could not determine from the docs whether NPP emits separate hole contours per label or how it nests them. This needs a run test with a donut image. I also did not find which CUDA/NPP release introduced the contour functions.
- **Not checked** (search budget exhausted): Open3D (my tree query failed), Halide/Taichi examples, Earth2Studio, NVIDIA cuLitho internals (closed), and any GDAL GPU branch. I found no evidence any of them provides raster→polygon with holes, but this is not verified.
- **GitHub code search** (as opposed to repository search) needs authentication and was not run. Small un-starred CUDA polygonizers inside larger repos could exist.

---

## Q2. What exactly do the closest papers do? (Applied Sciences 2025; J. Supercomputing 2022 and 2026)

### Takeaway
- **Ji, Niu and Chen (Applied Sciences 2025)** is **not raster→polygon.** It is **segment→polygon reconstruction** for IC-layout Boolean operations: oriented line segments go into a GPU hash table, then each thread walks one start. Holes are matched by bounding box plus ray casting. The 17× Manhattan speedup is versus Boost.Polygon on CPU, including PCIe transfers.
- **The 2026 "triad" paper (Alonso-Jordá et al.)** is fully on-GPU binary border tracking. Its contour-identification step is asynchronous **pointer jumping over a successor permutation**, essentially the user's proposed ring-ID step. It computes no hole hierarchy, is only tested on small images, and reaches roughly CPU parity (1.13× step speedup).

Neither paper does the user's combination: exact rectilinear rings plus sort-based pairing plus log-depth ring IDs plus CCL-based nesting.

### Cited Findings
**Ji, R.; Niu, Z.; Chen, L., "GPU-Accelerated Algorithm for Polygon Reconstruction," Appl. Sci. 2025, 15, 1111 (Institute of Microelectronics, CAS) [PAPER]**
- **Problem and output.**
  - Input: "a set of line segments," which "must be oriented in the correct direction." Polygon reconstruction is "the process of connecting a set of unordered line segments to form a closed polygon," used inside polygon Boolean operations for EDA/OPC/ILT.
  - Output: outer contours CCW and inner contours CW, plus outer↔inner association.
  
  — [doi:10.3390/app15031111](https://doi.org/10.3390/app15031111); full text via [MDPI-res PDF](https://mdpi-res.com/d_attachment/applsci/applsci-15-01111/article_deploy/applsci-15-01111.pdf)
- **Method:**
  1. GPU hash table with key = segment start and value = segment end. Hash H(x,y) = (x·length + y) % capacity, with capacity about 1.5× the segment count (load factor ~0.7).
  2. Each thread takes a start vertex and follows the hash links until it returns. To avoid duplicates, a walk is valid only if its start is the rightmost (then upper-most) vertex; a thread that meets a vertex further upper-right aborts. There are two rounds: count, then allocate from a memory pool and write.
  3. Orientation is taken at the extreme vertex, CW = hole.
  4. Hole matching: sort headers by bbox left x. Each hole thread scans backward for an outer whose bbox encloses it, then confirms by ray casting.
  5. A final sort places each hole after its outer.
  
  Touching polygons need a multimap and a policy choice. — [Ji et al. 2025](https://mdpi-res.com/d_attachment/applsci/applsci-15-01111/article_deploy/applsci-15-01111.pdf)
- **Complexity, as claimed by the authors:**
  - The construction phase costs O(p) per thread, where p is the largest polygon's vertex count.
  - Hole matching is O(n) worst case.
  - Total "O(n)" (their accounting).
  
  The validity-check scheme means every vertex's thread walks until it hits an upper-right vertex. Total work can therefore exceed O(n) for long rings. — [Ji et al. 2025](https://mdpi-res.com/d_attachment/applsci/applsci-15-01111/article_deploy/applsci-15-01111.pdf)
- **Experimental setup.** GPU: RTX 2080 Ti. CPU baseline: Boost.Polygon (`polygon_90_set_data` for Manhattan) on an Intel Ultra7 155H. The test polygons came from OpenILT model outputs converted "using edge detection techniques," then replicated. — [Ji et al. 2025](https://mdpi-res.com/d_attachment/applsci/applsci-15-01111/article_deploy/applsci-15-01111.pdf)
- **Numbers (average speedups):**
  - Non-Manhattan: 46.1× without holes, 45.2× with holes.
  - Manhattan: 17.9× without holes, 16.8× with holes.
  - Largest Manhattan run: 608,400 polygons / 10.95M segments in 150.9 ms GPU vs 2638.3 ms CPU.
  - Manhattan with holes: 966,575 polygons / 6.22M segments in 109.7 ms vs 1750.4 ms.
  
  — [Ji et al. 2025, Tables 1–4](https://mdpi-res.com/d_attachment/applsci/applsci-15-01111/article_deploy/applsci-15-01111.pdf)
- **Bottlenecks.**
  - The PCIe link was "only around 4 GB/s," and data transfer is "a critical factor limiting the algorithm's performance," especially for Manhattan sets (about 18 or 6.4 segments per polygon).
  - Pinned-memory allocation added 40% to the core step on Manhattan data.
  - Figure 6 shows the transfer and pin-allocation shares; segment connection is about 40–68% of time depending on the set.
  
  — [Ji et al. 2025](https://mdpi-res.com/d_attachment/applsci/applsci-15-01111/article_deploy/applsci-15-01111.pdf)
- **Code availability.** No code is released. "Data ... available by the authors on request." — [Ji et al. 2025](https://mdpi-res.com/d_attachment/applsci/applsci-15-01111/article_deploy/applsci-15-01111.pdf)

**Alonso-Jordá, Díaz-Cano Lozano, Quintana-Ortí, Folch, "Accelerated border tracking in binary images with GPUs," J. Supercomputing 2026 (CC-BY) [PAPER via Research Square preprint; META for the published version]**
- **Method.** The image is split into blocks. Kernel 1 extracts "triads" (prev, current, next contour pixels) using a 256-entry lookup table over the 8-neighbourhood; this beat an automaton by 1.37–1.57×. Kernel 2 connects triads within and across blocks. Kernel 3 performs "contour identification." — [Research Square preprint rs-7694436](https://www.researchsquare.com/article/rs-7694436/latest.pdf); [doi:10.1007/s11227-026-08488-4](https://doi.org/10.1007/s11227-026-08488-4)
- **Contour identification:**
  - It takes the "next index" vector v (a permutation whose cycles are contours).
  - Algorithm 3 has each thread follow v until it returns to itself or drops below its own index. This was "clearly disappointing."
  - Algorithm 4 is an **asynchronous pointer-jumping** variant: `x(t) ← x(x(t))` while `x(t) > t`, without global synchronization. The authors argue it is correct without barriers.
  - The authors state that the cost "depends on the order of the indices": worst case n (sequential), best case log2 n.
  
  — [preprint](https://www.researchsquare.com/article/rs-7694436/latest.pdf)
- **Numbers** (Quadro RTX 5000; test images with 1,277–10,332 contour points):
  - Triad extraction: 43.6–102 µs.
  - Connection: 2.2–8.3 µs.
  - Identification: improved from 100–1024 µs to 11.7–75.4 µs.
  - In an industrial machine-vision workflow (RTX A4000, 1,000 images), border tracking went from 5.87 to 5.21 ms (1.13×). "Process Image" sped up 2.16× by avoiding transfers, for 1.20× globally.
  
  — [preprint](https://www.researchsquare.com/article/rs-7694436/latest.pdf)
- **No hierarchy.** I found no hole/parent/hierarchy computation in the text: no matches for "inner," "outer," "parent" or "nest." — [preprint](https://www.researchsquare.com/article/rs-7694436/latest.pdf)

**Garcia-Molla, Alonso-Jordá, García-Laguía, "Parallel border tracking in binary images using GPUs," J. Supercomputing 2022 [META]**
- This is a parallel Suzuki variant: one thread per small rectangle, then cross-rectangle connection. On the fastest GPU it was ">10 times faster than the sequential CPU routine running in the fastest CPU." — [doi:10.1007/s11227-021-04260-y](https://doi.org/10.1007/s11227-021-04260-y)

### Inferences
- **The "17× on rectilinear with holes" number is not a raster-polygonizer benchmark.** It measures linking pre-made oriented segments, including about 4 GB/s PCIe transfers, against Boost.Polygon. The user's engine produces its own segments from the bitmap on device, so transfers vanish.
- **Ji et al.'s design choices to avoid:**
  - **Per-thread walks with "rightmost-start" dedup** have load imbalance that grows with ring length. The worst case is one ring with p vertices, so p serial steps.
  - **Bbox + ray-cast hole matching** is O(n) worst case per hole.
  
  The proposed design replaces them with pointer jumping (log depth) and CCL-derived parents (no geometry tests).
- **The triad paper is direct prior art for "ring ID via pointer jumping."** Its variant (jump only while `x(t) > t`) has a permutation-dependent worst case of n rounds. The fix is Wyllie-style **min-label propagation with pointer doubling** (label ← min(label, label[succ]); succ ← succ[succ]). It is guaranteed ⌈log2 L⌉ rounds for ring length L, and I prototyped it ([RUN], Q3). Alternatively, use Helman–JáJá sparse ruling sets (Q3).
- **Novelty statement that survives scrutiny:** "first fully device-resident raster→polygon engine that outputs exact pixel-edge rectilinear rings with outer/hole classification and parent links, using only sort/scan/pointer-jumping primitives with log-depth ring identification." Cite NPP, the 2026 triad paper and Ji et al. 2025 as nearest prior art.

### Gaps
- I read the 2026 triad paper via its Research Square preprint, not the Springer version of record. Numbers could differ slightly.
- The 2022 paper's full text was not read (abstract only). Its GPU models and exact speedups are unverified.
- I found no published comparison of NPP contours versus OpenCV findContours versus the triad method.

---

## Q3. Is the proposed sort → pair → pointer-jump → classify → CCL-parent design correct, and how fast could it be?

### Takeaway
I implemented the full proposed pipeline as a pure-Python reference using only map / sort / pair / pointer-jump / min-reduce / CCL lookups. It was **topologically exact on 3,000 random images, both 4- and 8-connectivity**:
- every outer ring matched a foreground component;
- every hole matched a bounded background component;
- every hole-parent assignment was valid;
- every sort-pair assertion held.

Throughput is bounded by a few streaming passes plus one 17-bit key sort, plus ⌈log2 L_max⌉ random-gather rounds over the corner array. On a modern datacenter GPU this plausibly lands in the tens to low hundreds of milliseconds for ~10^8 corners. That is an inference, not a measurement.

### Cited Findings
- **Prototype results [RUN].** `proto_sortpair.py` did the following:
  1. Emits corners (x, y, h-dir, v-dir) from each 2×2 quad. Diagonal quads emit two corners, paired by the connectivity rule: fg-8 wraps the background pixels, fg-4 wraps the foreground pixels.
  2. Sorts by (y, x, h) and pairs consecutive entries to get the horizontal partner.
  3. Sorts by (x, y, v) and pairs to get the vertical partner.
  4. Computes ring IDs by min-label pointer doubling over succ = vpartner∘hpartner.
  5. Classifies each ring as outer or hole by the pixel inside its lexicographically-min corner.
  6. Assigns a hole's parent as the fg-CCL label of the pixel above that corner. Each parent was validated by checking that the hole stays enclosed when only the parent component is foreground.
  
  Results:
  - fg8: outer==#fgCCL 1500/1500, holes==#bounded-bgCCL 1500/1500, parent valid 1500/1500.
  - fg4: 1500/1500 on all three.
  - The pairing assertions (pair = R-then-L on the same row; D-then-U on the same column) never failed.
  
  — local script `/private/tmp/claude-501/r2scripts/proto_sortpair.py`
- **GPU sort throughput [META].** Onesweep (Adinets & Merrill, NVIDIA) reaches **29.4 GKey/s sorting 256M random 32-bit keys on A100**, about 1.5× CUB's previous LSD radix sort. It needs about 2n global memory operations per digit pass versus about 3n. — [arXiv 2206.01784](https://arxiv.org/abs/2206.01784)
- **GPU list ranking [PAPER].** Wei & JáJá adapt Helman–JáJá (random splitters → sequential sublists → splitter prefix → fix-up):
  - About **297 ms for a 64M-node random list** on a Tesla C1060 (30 SMs, **102 GB/s** peak), versus 271 ms for an ordered list. Ordered, stride and random layouts differ by at most about 10%.
  - Linear scaling from 1M to 256M nodes.
  - Best configuration: 30 blocks × 64 threads with s = 61,440 sublists.
  
  — [UMIACS-TR-2010-08](http://drum.lib.umd.edu/bitstream/1903/10600/1/UMIACS-TR-2010-08.pdf)
- **Other list-ranking prior work [META].** Rehman, Kothapalli and Narayanan, "Fast and scalable list ranking on the GPU," ICS 2009. — [doi:10.1145/1542275.1542311](https://doi.org/10.1145/1542275.1542311)
- **Triad paper caveat [PAPER].** Its asynchronous pointer jumping is n rounds in the worst case for adverse index orders. — [preprint](https://www.researchsquare.com/article/rs-7694436/latest.pdf)
- **Comparison point: Ji et al. 2025** linked 6.2M–11M Manhattan segments in about 110–151 ms on an RTX 2080 Ti, including ~4 GB/s PCIe transfers. — [Ji et al. 2025](https://mdpi-res.com/d_attachment/applsci/applsci-15-01111/article_deploy/applsci-15-01111.pdf)

### Inferences
- **Algorithmic simplifications found while prototyping (design constraints):**
  1. **Horizontal pairing needs no sort.** Compact corners in row-major order (one prefix sum). The result is already sorted by (y, x). At a diagonal quad, emit the "points-left" copy before the "points-right" copy, and partners are then consecutive pairs.
  2. **Vertical pairing** needs only a **stable radix sort on x** of the row-major corner list, which keeps y order. For a 100k-wide image the key is 17 bits, so 2–3 Onesweep digit passes (8-bit digits). The v-direction tie-break at diagonal quads comes from emission order.
  3. **Ring IDs:** use min-label pointer doubling over succ = vpartner∘hpartner, which takes every second corner, then fold in the H-partner. Rounds = ⌈log2(L_max/2)⌉. Do not use the triad paper's order-dependent variant.
     - A work-efficient alternative is Helman–JáJá / Wei–JáJá sparse ruling sets: O(n) work, and practical per the 2010 numbers.
     - Ring order (rank within ring) comes from the same pass with sum instead of min, or from a second list-rank. Output compaction is then a sort by (ring_id, rank) or a direct scatter.
  4. **Outer vs hole:** a segmented argmin per ring of (y, x). A hole is a ring whose min corner has background at its bottom-right pixel. The shoelace sign is equivalent once the ring is directed, but the min-corner test needs no area accumulation.
  5. **Hole parent:** the fg-CCL label of the pixel directly above the hole's min corner (fallback: left), with one CCL on the foreground only. A background CCL is not needed for parenting; it is only needed for verification.
  6. **Topology certificate:** #outer = #fg components; #holes = #bounded bg components; (#outer − #holes) = χ from bit-quads (Q7). This is a free runtime self-check.
- **Back-of-envelope cost for a 100k × 100k mask (10^10 pixels).** This is my inference. Device bandwidth B is a parameter, not a cited spec.
  - Quad classification plus corner count plus compaction: about 2–3 streaming passes over the mask. Bit-packed that is 1.25 GB per pass; byte masks are 10 GB per pass. At multi-TB/s, that is tens of ms or less.
  - Corner count V is content-dependent (= n1 + n3 + 2nD).
  - Pointer doubling costs about ⌈log2 L_max⌉ × V random 8–16 B gathers. For V = 10^8 and L_max ~ 10^6 that is about 20 rounds × 10^8 gathers.
  - A naive bandwidth scaling of Wei–JáJá's 2010 number (64M random nodes in 297 ms at 102 GB/s) to a ~30× higher-bandwidth GPU suggests ~10 ms per 64M-node list-ranking pass. Pointer doubling is a few times more expensive than Helman–JáJá. So **~10^8 corners in roughly 30–200 ms on an A100/H100-class GPU** is a plausible envelope, and CCL is in the sub-second range at this scale.
  - Memory: V × (~24–32 B) = 2.4–3.2 GB for 10^8 corners, plus the mask and labels. This fits on 40–80 GB cards. For larger masks, tile (Q5).
  - These are **unmeasured estimates**. Benchmark against NPP and against the user's current CPU linker.

### Gaps
- I verified no modern (Ampere/Hopper) GPU list-ranking or pointer-jumping benchmark. The search budget was exhausted. Wei & JáJá (2010) is the only measured number in these notes.
- I did not verify CUB `DeviceRadixSort` pair-sort throughput for 64-bit keys + 32-bit values on H100/A100/Orin. Only the Onesweep 32-bit-key A100 figure is verified.
- The prototype is a CPU reference for correctness, not a GPU timing. A CUDA implementation and benchmark remain to be done.
- I did not quantify the real corner density ρ = V/N for target datasets (building masks, SAM masks, IC layouts). This determines everything after classification.

---

## Q4. Multi-label planar partitions and topology-preserving simplification on GPU: any prior work?

### Takeaway
I found **no GPU work** on multi-label (label-map → shared-edge planar partition) polygonization or on shared-arc, topology-preserving line simplification. The only GPU "simplification" in a mainstream library (cuCIM `approximate_polygon`) is a host-driven serial DP loop. Every CPU tool that produces gap-free multi-label output (GDAL, xarray-spatial) does so by tracing exact pixel edges, so shared boundaries are bit-identical integers. The sort-pair design extends to label maps with a modest change to the sort key.

### Cited Findings
- **cuCIM `approximate_polygon`** is serial DP with per-segment host syncs (`cp.asnumpy`). — [cuCIM _polygon.py](https://github.com/rapidsai/cucim/blob/main/python/cucim/src/cucim/skimage/measure/_polygon.py) [SRC]
- **GDAL** polygonizes label rasters region by region, with edges "exactly on pixel boundaries." Neighbouring polygons therefore share identical vertex coordinates. — [GDAL polygonize.cpp](https://github.com/OSGeo/gdal/blob/master/alg/polygonize.cpp) [SRC]
- **xarray-spatial** polygonizes all values of a label raster. Integer rasters group by strict equality, floats by atol/rtol. — [xrspatial/polygonize.py](https://github.com/makepath/xarray-spatial/blob/master/xrspatial/polygonize.py) [SRC]
- **Round-1 evidence:**
  - Jiao, Zuo and Persello (arXiv 2609.07520, Sept 2026) blame independent per-instance raster→vector conversion for gaps and overlaps, and rebuild faces from a common planar graph (CPU).
  - Constrained/contextual DP is CPU-only.
  - GPU topology-preserving simplification exists only for triangle meshes (PLOS ONE 2021).
  
  — see `algorithms_ml.md` Q2: [arXiv 2609.07520](https://arxiv.org/abs/2609.07520); [PLOS ONE 2021](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0255832)

### Inferences
- **Multi-label extension of the sort-pair design:**
  - Corners become (label, x, y, h, v). A quad emits a corner for label ℓ wherever ℓ's region turns at that vertex.
  - With 4 distinct labels at a vertex, it emits up to 4 corners, one per label region touching it.
  - Pairing sorts by (ℓ, y, x) and (ℓ, x, y), so it runs independently per label in one batched sort.
  - Adjacent faces then share exact integer vertices and edges, so there are no gaps or overlaps by construction. This is the GPU analogue of what GDAL does on CPU.
- **Shared-arc (TopoJSON-style) output for topology-preserving simplification:**
  - Junction vertices are quads with ≥3 distinct labels, plus 2-label diagonal quads.
  - Every unit boundary edge has a (left label, right label) pair. Arcs are maximal edge chains between junctions with constant label pair, which is again a list-ranking problem over edges.
  - Simplify each arc once (parallel DP / Visvalingam, one arc per thread or warp, with junctions pinned).
  - Then reassemble faces by referencing arcs. Adjacent faces simplify identically, so no gaps appear.
  - Self-intersection avoidance between arcs would still need a constraint check (e.g., a GPU segment-intersection test after simplification, with rewind on conflict). I found no GPU precedent for this.
- **Rectilinear-preserving "simplification" is nearly free:** removing collinear vertices is inherent, since only corners are emitted. Stair-step collapse at a tolerance can be done per arc with a segmented scan.

### Gaps
- No GPU Visvalingam-Whyatt, GPU shared-arc extraction, or GPU planar-partition polygonizer was found. Search coverage this round was limited by the exhausted search budget, so treat this as "not found," not "proven absent."
- I could not check TopoJSON/mapshaper internals or QGIS/GRASS `r.to.vect` algorithms this round.

---

## Q5. Exact tile-seam stitching for out-of-core or multi-GPU polygonization

### Takeaway
The only working open-source exact stitcher I verified (xarray-spatial's dask path) uses **unit-edge cancellation**:
1. Split boundary-touching rings into unit directed edges.
2. Cancel opposite edges across chunk seams.
3. Re-trace the rest.
4. Re-group exteriors and holes.

It is exact because coordinates are integers, but it runs on CPU in Python/numba. GDAL avoids tiling by streaming scanlines with two line buffers. A GPU design can do exact stitching with the same sort-pair primitive applied only to seam-crossing chain endpoints.

### Cited Findings
- **xarray-spatial `_merge_polygon_rings`** "Uses edge cancellation: splits all rings into unit-length directed edges, cancels opposing edges (which occur at chunk boundaries where the same value continues across), and traces the remaining edges into closed rings." `connectivity_8` selects "the degree-4 vertex pairing rule" for figure-8 rings at chunk corners. — [xrspatial/polygonize.py](https://github.com/makepath/xarray-spatial/blob/master/xrspatial/polygonize.py) [SRC]
- **Seam grouping.** xarray-spatial groups boundary polygons by "spatial topology + value closeness union-find." For 4-connectivity, polygons must share a unit edge; for 8-connectivity, "sharing any integer-grid vertex is enough." It also fixes float-tolerance chain bugs across chunks (issues #2171/#2173/#2583 referenced in code). — [xrspatial/polygonize.py](https://github.com/makepath/xarray-spatial/blob/master/xrspatial/polygonize.py) [SRC]
- **Interior vs boundary results.** Chunks return interior polygons, which are final, separately from boundary polygons, which need merging. Batches are scheduled with `dask.compute`. — [xrspatial/polygonize.py](https://github.com/makepath/xarray-spatial/blob/master/xrspatial/polygonize.py) [SRC]
- **xarray-spatial contours** add "a NaN 1-cell halo" to dask chunks for marching squares, then deduplicate segments and lines across chunks. — [xrspatial/contour.py](https://github.com/makepath/xarray-spatial/blob/master/xrspatial/contour.py) [SRC]
- **GDAL** processes the raster as a scanline stream (`ProcessLine` over the last and current lines; `TwoArm` arrays of width nXSize+2) to "minimize memory use so that very large rasters can be processed." — [GDAL polygonize.cpp](https://github.com/OSGeo/gdal/blob/master/alg/polygonize.cpp) [SRC]
- **Triad paper.** The 2026 paper connects triads "between adjacent blocks" on GPU. This is intra-image block stitching, not out-of-core. — [preprint](https://www.researchsquare.com/article/rs-7694436/latest.pdf) [PAPER]

### Inferences
- **A GPU-native exact seam protocol:**
  1. **Ownership.** Each tile owns the quads whose top-left pixel lies in the tile, and reads a 1-pixel halo, as the user's engine already does. Every global grid vertex is then emitted by exactly one tile, so there are no duplicate corners.
  2. **Local linking.** Rings fully inside a tile close locally and are final, as in xarray-spatial's interior set. Chains that leave the tile end at a seam edge-crossing, so emit a "port" record (seam id, integer coordinate along the seam, direction, local chain id).
  3. **Global stitch.** Sort all ports by (seam, coordinate). Partners are adjacent pairs, the same pairing rule as Q3. This yields a small chain-graph whose cycles are found with the same pointer doubling.
  4. **Labels across seams.** Merge CCL labels across seams with a union-find over seam-adjacent label pairs (the standard tiled-CCL approach). Hole parents then resolve to global labels.
- **Edge cancellation versus port pairing.** Edge cancellation is simpler but materializes unit edges, so it costs O(perimeter in pixels). Port pairing touches only seam crossings, O(#crossings), which is preferable on GPU and multi-GPU. Both are exact because coordinates are integer pixel-grid coordinates. Any float georeferencing must be applied **after** stitching.
- **Connectivity at seam corners.** The diagonal-quad pairing rule must be identical in every tile (the same fg8/fg4 rule as in the prototype). Otherwise figure-8 configurations at seam corners split or merge inconsistently. xarray-spatial explicitly handles this case.

### Gaps
- I found no published GPU or multi-GPU out-of-core polygonization paper this round (search budget exhausted). "Streaming polygonization" or "seamless vectorization of large rasters" literature could not be searched.
- GDAL's multi-threaded or tiled behaviour, if any, was not checked beyond the core polygonize source.

---

## Q6. GPU-native output formats: GeoArrow and WKB

### Takeaway
GeoArrow's Polygon layout is exactly what a GPU linker naturally produces: coordinates plus ring offsets plus polygon offsets, all computable with prefix sums. WKB is also a fixed-size, scan-computable byte layout. The GPU geospatial consumers I verified take these formats: cuSpatial (archived) used GeoArrow, and SedonaDB's libgpuspatial ingests WKB ArrowArrays. So a GPU polygonizer can hand off zero-copy, with no CPU in the loop.

### Cited Findings
- **GeoArrow Polygon layout.** An array of Polygons is `List<List<Coordinate>>`: "each element of the array (Polygon) is a list of rings (the first ring is the exterior ring, optional subsequent rings are interior rings), and each ring is a list of coordinate vertices." Child names are "rings" and "vertices." "The first coordinate and the last coordinate of a ring must be identical (i.e., rings must be closed)." — [GeoArrow format.md](https://github.com/geoarrow/geoarrow/blob/main/format.md) [DOC]
- **GeoArrow coordinates.** They can be separated `Struct<x: double, y: double>` (recommended) or interleaved `FixedSizeList<double>[n_dim]`. — [GeoArrow format.md](https://github.com/geoarrow/geoarrow/blob/main/format.md) [DOC]
- **cuSpatial** represented data in GeoArrow format and has a `cuspatial/io/pygeoarrow.py` module. The repo is now archived. — [rapidsai/cuspatial](https://github.com/rapidsai/cuspatial) [SRC]
- **SedonaDB libgpuspatial** "takes two inputs ... ArrowArrays containing geometries in WKB format," and supports Point through MultiPolygon. — [libgpuspatial README](https://github.com/apache/sedona-db/tree/main/c/sedona-libgpuspatial/libgpuspatial) [DOC]
- **Ji et al.'s output layout** is a header array (index, vertex count, hole flag, memory offset, bbox), then free space, then vertices. Outers are CCW, holes CW, and each hole is stored after its outer. This is a custom GPU layout comparable to GeoArrow offsets. — [Ji et al. 2025](https://mdpi-res.com/d_attachment/applsci/applsci-15-01111/article_deploy/applsci-15-01111.pdf) [PAPER]

### Inferences
- **GeoArrow on GPU, build steps:**
  1. Ring vertex counts (+1 for closure) give ring offsets by exclusive scan.
  2. Sort rings by (parent_component, is_hole, ring_id) so each outer is followed by its holes.
  3. Rings-per-polygon (1 + #holes, a segmented count) gives polygon offsets by exclusive scan.
  4. Scatter coordinates.
  
  This is 2 scans + 1 sort + 1 scatter, and it maps onto cuDF list-of-list columns. Orientation convention: GeoArrow and OGC do not mandate winding. Pick CCW-outer/CW-hole (matching Ji et al. and xarray-spatial) or follow RFC 7946 for GeoJSON consumers.
- **WKB on GPU.** Per polygon the size is 1 (byte order) + 4 (type) + 4 (numRings) + Σ_rings [4 + 16·(n_vertices+1)] bytes for 2D doubles. A scan over per-polygon sizes gives byte offsets, and each thread writes its own polygon. The layout is fully deterministic and needs no CPU. Integer pixel coordinates should be transformed by the affine geotransform during the write, after exact stitching.
- **Compact interchange alternative:** a rectilinear ring only needs alternating H/V coordinates. That gives about 50% coordinate savings over GeoArrow before compression, but it is non-standard. Offer it as an internal or ML format only.

### Gaps
- I did not verify any GPU WKB *writer* implementation in cuDF or cuSpatial. I only verified that libgpuspatial ingests WKB; whether its WKB parsing runs on device was not checked.
- I did not check the current cuDF/libcudf status of GeoArrow extension types after cuSpatial's archival.

---

## Q7. Bit-quad (Gray 1971) formulas for Euler number, area, perimeter and corner count; the Fast-χ paper

### Takeaway
I verified all the bit-quad identities the engine relies on by brute force against independent flood-fill labeling ([RUN], 2,000 random images each, 100% agreement):
- **χ₄ = (n1 − n3 + 2nD)/4** for 4-connected foreground (holes as 8-connected background).
- **χ₈ = (n1 − n3 − 2nD)/4** for 8-connected foreground (holes as 4-connected background).
- **Area = (n1 + 2n2 + 3n3 + 4n4 + 2nD)/4**.
- **Exact pixel-edge perimeter = n1 + n2 + n3 + 2nD**.
- **Rectilinear corner count = n1 + n3 + 2nD**.

Here n2 counts only edge-adjacent pairs and nD the two diagonal patterns, over the zero-padded image. Fast-χ (Li et al., IEEE TMI) uses the minus-sign (8-connectivity) Gray form inside a CNN.

### Cited Findings
- **Executable verification [RUN].** All five identities above held 2000/2000 against flood-fill CCL; holes are counted as bounded background components in the complementary connectivity. My first guess for perimeter (n1 + 2n2 + n3 + 2nD) failed at 277/2000, which confirms the test is discriminating. — local script `/private/tmp/claude-501/r2scripts/verify_bitquad3.py`
- **Primary reference [META].** S. B. Gray, "Local Properties of Binary Images in Two Dimensions," IEEE Transactions on Computers, 1971. Closed access; the full text was not read. — [doi:10.1109/T-C.1971.223289](https://doi.org/10.1109/T-C.1971.223289)
- **Fast χ [PAPER].** Li, Ma, Ouyang, Paetzold, Rueckert and Kainz, "Topology Optimization in Medical Image Segmentation with Fast χ Euler Characteristic" (arXiv 2507.23763; IEEE TMI):
  - **Gray method in 2D:** "χ can be computed from the occurrences of only 10 bit-quads #Qi as ... χ = ¼(#QA − #QB − 2#QC)," where QA = {Q2..Q5} (single-pixel quads), QB = {Q12..Q15} (three-pixel quads) and QC = {Q10, Q11} (diagonals).
  - **3D extension:** bit-octets, with 22 coefficient groups verified empirically.
  - **χ-Net:** reproduces the bit-quads with fixed conv kernels K_A, K_B ∈ {−1,1}^{2×2×4} and K_C ∈ {−1,1}^{2×2×2}, plus pooling to get local χ maps. These feed a "topological violation map" and a correction network.
  - **Code:** github.com/smilell/Topology-aware-Segmentation-using-Euler-Characteristic.
  - The paper text does not name the connectivity. Its minus sign matches the χ₈ form verified above.
  
  — [arXiv 2507.23763](https://arxiv.org/abs/2507.23763); [PDF](https://arxiv.org/pdf/2507.23763)
- **Ohser LUT formulation in skimage/cuCIM.** scikit-image/cuCIM `euler_number` uses LUT coefficients from Ohser, Nagel & Schladitz (2002), "The Euler Number of Discretized Sets — On the Choice of Adjacency in Homogeneous Lattices":
  - `EULER_COEFS2D_4 = [0,1,0,0,0,0,0,-1,0,1,0,0,0,0,0,0]` and `EULER_COEFS2D_8 = [0,0,0,0,0,0,-1,0,1,0,0,0,0,0,-1,0]`, indexed by 2×2 configuration code.
  - This is an integral-geometry variant that weights only some orientations, so it is not Gray's symmetric form, but it is equivalent in total.
  - cuCIM runs this on GPU.
  
  — [cuCIM _regionprops_utils.py](https://github.com/rapidsai/cucim/blob/main/python/cucim/src/cucim/skimage/measure/_regionprops_utils.py) [SRC]
- **cuCIM `perimeter`** follows Benkrid & Crookes. It uses weights of 1, √2 and (1+√2)/2 on 3×3 neighbourhood codes, which is a smoothed estimate, not the exact pixel-edge count. — [cuCIM _regionprops_utils.py](https://github.com/rapidsai/cucim/blob/main/python/cucim/src/cucim/skimage/measure/_regionprops_utils.py) [SRC]

### Inferences
- **Free counts from the classification kernel.** One quad histogram (n0..n4, nD) gives:
  - exact output size: corners V = n1 + n3 + 2nD;
  - exact area and exact edge perimeter;
  - χ for either connectivity.
  
  Uses:
  - V pre-sizes every downstream buffer exactly, with no Ji-style "n/3 headroom."
  - χ plus the CCL counts gives a runtime correctness certificate: #outer − #holes = χ (Q3).
- **Connectivity convention must be consistent end-to-end.** The diagonal-quad corner pairing (fg8 vs fg4), the CCL connectivity used for parents, and the χ sign must all agree. The prototype shows the fg8/bg4 and fg4/bg8 pairs are each internally consistent.

### Gaps
- Gray's 1971 paper itself was not read (closed access). The formulas are verified computationally and via Fast-χ and the Ohser LUT code, not from the primary text. Pratt's *Digital Image Processing* is the usual secondary citation but was not fetched.
- Fast-χ runtime numbers versus persistent homology were not located in the text extracted.
