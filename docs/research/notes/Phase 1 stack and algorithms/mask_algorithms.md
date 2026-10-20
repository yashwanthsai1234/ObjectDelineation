# Mask-side algorithms for the exact geometry reward checker (CCL, Euler/Betti, matching, band/boundary metrics, 2D and 3D)

Scope: the mask half of the reward. A predicted binary mask (rasterized VLM polygon) is scored against 2–6 expert masks. Convention: 8-connected foreground, 4-connected background (26/6 in 3D). Budget: about 256 masks per GRPO step, well under 1 s on 8 CPU cores. Date: 2026-10-18.

**Local evidence.** Timings below come from a micro-benchmark run on this machine: Apple M3, single thread, OMP_NUM_THREADS=1, cv2.setNumThreads(1) (correction, 2026-10-19: on the macOS OpenCV build `setNumThreads(1)` leaves OpenCV at 8 threads and only `setNumThreads(0)` gives 1, so these OpenCV timings used 8 threads; see `experiments/phase1/impl/fixes/F1.md`, I14). Library versions: numpy 2.x, scipy 1.18.1, scikit-image 0.26.0, connected-components-3d 5.0.0, opencv-python-headless 5.0.0, edt. Scripts are in `docs/research/notes/Phase 1 stack and algorithms/prototypes/` (`bench_mask.py`, `holes_fix.py`, `band_demo.py`), cited below as [local bench](prototypes/bench_mask.py) and [local band demo](prototypes/band_demo.py). Synthetic masks come in two kinds. "blob" is one smooth object with occasional holes, like a polygon. "noisy" is smoothed noise thresholded at 0, which gives hundreds of components and holes and serves as a stress test. They are not real LIDC/RIGA masks.

## Q1. Connected-component labeling: what is fastest for small crops in Python?

### Takeaway
Use OpenCV's `connectedComponentsWithAlgorithm(..., ccltype=CCL_SPAGHETTI)` (or CCL_BBDT) for 2D. Measured locally it is about 3x faster than cc3d and 5–7x faster than scipy.ndimage.label on 64–256 px crops (2.1 µs vs 6.0 µs vs 13.8 µs at 64²), and about 4x faster than scipy at 1024². Use cc3d for 3D 26/6, where it is 6–8x faster than scipy. GPU CCL is unnecessary at this scale.

### Cited Findings
- YACCLAB is the reference CCL benchmark. It ships 2D/3D CPU algorithms (SAUF 2009, BBDT 2010, CTB, PRED 2016, DRAG 2018, Spaghetti 2019, Tagliatelle 2021, PRED++ 2021, BRTS/BMRS 2021, EPDT 3D) and GPU algorithms (UF, Komura Equivalence 2018, BUF/BKE 2019, C-SAUF/C-BBDT/C-DRAG, HA4/HA8, …). The README publishes no ranking table: results are regenerated per machine. — [YACCLAB repo](https://github.com/prittt/YACCLAB)
- Spaghetti labeling (Bolelli, Allegretti, Baraldi, Grana, IEEE TIP 29:1999–2012, 2020) reportedly surpasses all compared approaches in all settings, on synthetic and real datasets. — [IRIS record](https://iris.unimore.it/handle/11380/1181437); [PDF](https://iris.unimore.it/retrieve/e31e124d-c7be-987f-e053-3705fe0a095a/2018_TIP_Spaghetti%20Labeling%20Directed%20Acyclic%20Graphs%20for%20Block-Based%20Connected%20Components%20Labeling.pdf)
- The OpenCV header (4.x) maps the CCL flags as follows. `CCL_DEFAULT` and `CCL_BOLELLI`/`CCL_SPAGHETTI` select Spaghetti for 8-way and Spaghetti4C (Bolelli 2021) for 4-way. `CCL_WU`/`CCL_SAUF` select SAUF for both. `CCL_GRANA`/`CCL_BBDT` select BBDT for 8-way and SAUF for 4-way. Parallel versions exist (Bolelli 2017). — [opencv imgproc.hpp](https://raw.githubusercontent.com/opencv/opencv/4.x/modules/imgproc/include/opencv2/imgproc.hpp)
- cc3d (Seung lab) supports 2D 4/8 and 3D 6/18/26 connectivity, plus periodic boundaries. Its core is a two-pass SAUF with decision trees after Wu/Otoo/Suzuki, extended to 3D, and union-find with path compression but no union-by-rank. It adds "phantom labeling" for 4/6-connectivity. It also provides `statistics()` (counts, bounding boxes, centroids), `dust()` (size filter), `largest_k()`, `contacts()`, `region_graph()` and `voxel_connectivity_graph()`. The installed 5.0.0 has no Euler-number function (checked with `dir(cc3d)`). — [cc3d repo](https://github.com/seung-lab/connected-components-3d); [local bench](prototypes/bench_mask.py)
- Local timings, single thread, µs per call, 8-connected labeling of one mask. Columns: scipy `ndi.label`, cv2 Spaghetti, cv2 SAUF, cv2 BBDT, cc3d, skimage `label`. — [local bench](prototypes/bench_mask.py)
  - 64² blob: 13.8 / 2.1 / 3.2 / 2.4 / 6.0 / 19.4
  - 128² blob: 39.9 / 8.2 / 12.5 / 8.1 / 16.0 / 45.1
  - 256² blob: 139 / 38 / 48 / 30 / 53 / 140
  - 512² blob: 552 / 138 / 177 / 119 / 182 / 570
  - 1024² blob: 1946 / 564 / 671 / 451 / 671 / 1968
  - 1024² noisy: 4243 / 1342 / 1754 / 1401 / 1320 / 4200
  - skimage `label` is no faster than scipy; at these sizes it costs about the same.
- Correctness: on 400 random masks (32–100 px, blob and noisy), the labelings from cv2 Spaghetti, cv2 SAUF, cc3d and skimage matched scipy's 8-connected labeling up to permutation, with 0 mismatches. — [local bench](prototypes/bench_mask.py)
- API pitfall: opencv-python-headless 5.0.0 rejects the `ccltype=` keyword on `cv2.connectedComponents`. Use `cv2.connectedComponentsWithAlgorithm(img_u8, 8, cv2.CV_32S, cv2.CCL_SPAGHETTI)`. — [local bench](prototypes/bench_mask.py) (observed error: "unexpected keyword argument 'ccltype'")
- GPU: Allegretti, Bolelli and Grana, "Optimized Block-Based Algorithms to Label Connected Components on GPUs", IEEE TPDS 31(2):423–438. BKE/BUF work on 2×2 (2D) and 2×2×2 (3D) blocks to cut memory accesses and atomics. A survey with code exists. A `bkeccl` PyPI package exists (not inspected). — [TPDS PDF](https://iris.unimore.it/bitstream/11380/1179616/1/2018_TPDS_Optimized_Block_Based_Algorithms_to_Label_Connected_Components_on_GPUs.pdf); [GPU CCL review](https://iris.unimore.it/retrieve/a031cde1-6b22-43db-b3c7-45a5a254c791/A_State-of-the-Art_Review_With_Code_About_Connected_Components_Labeling_on_GPUs.pdf); [bkeccl](https://pypi.org/project/bkeccl/)

### Inferences
- At 64–256 px, Python call overhead dominates. cv2 makes about 2 µs per call; cc3d's fixed overhead is about 4 µs and scipy's about 12 µs. Even scipy is affordable in absolute terms: 256 masks × 2 labelings × 140 µs ≈ 70 ms on one core at 256². CCL is not the bottleneck. EDT and per-component bookkeeping are (see Q2 and Q5).
- BBDT and Spaghetti are within noise of each other on blob masks, consistent with YACCLAB's general finding that block-based trees win. Spaghetti4C covers the 4-connected background pass, so one library handles both passes.
- GPU CCL (Komura, Playne–Hawick, BKE) only pays off when masks already live on the GPU in large batches. At 256 masks of 64–1000 px, the transfer and launch overhead would exceed CPU time. Treat it as a Phase-later CUDA item, not a requirement.

### Gaps
- I did not retrieve YACCLAB's published per-dataset millisecond tables (from the Spaghetti/Tagliatelle/PRED++ papers), so I cannot quote their numbers here. Tagliatelle, PRED++ and BMRS are not exposed in OpenCV, and I did not benchmark them.
- I did not verify Playne–Hawick (2018 TPDS) performance claims.
- All timings are from synthetic masks on an M3. Linux x86 (the training host) may rank cc3d and cv2 differently, so benchmark there.

## Q2. Euler number and Betti numbers: per-component β1, and is persistent homology needed?

### Takeaway
For 2D, compute β0 by 8-connected CCL and per-component β1 (holes per instance) by 4-connected labeling of the padded background. Assign each hole to its enclosing component with a "pixel-above-the-hole's-first-raster-pixel" lookup. This is exact, verified against brute force, and faster than the numpy bit-quad Euler. Keep the bit-quad χ only as an independent certificate check. `skimage.measure.euler_number(m, connectivity=2)` uses the project's 8/4 convention in 2D, and `connectivity=3` gives 26/6 in 3D (verified). Persistent homology is overkill for exact binary-mask topology.

### Cited Findings
- `skimage.measure.euler_number(image, connectivity)` semantics, read from the 0.26 source. `connectivity` is the "maximum number of orthogonal hops": in 2D 1→4-neighborhood and 2→8; in 3D 1→6 and 3→26, and 2 raises NotImplementedError. The docstring states "If object is 4-connected, then background is 8-connected, and conversely". The image is binarized (`>0`) and zero-padded. Configurations are built with `ndi.convolve` over 2×2 (16 bins) or 2×2×2 (256 bins) windows, and an Ohser et al. (2002) coefficient LUT is applied, with the 3D result multiplied by 0.125. — [skimage source `_regionprops_utils.py`](https://github.com/scikit-image/scikit-image/blob/main/skimage/measure/_regionprops_utils.py)
- Empirical semantics check: [local bench](prototypes/bench_mask.py)
  - 2D diagonal pair gives χ=1 with conn=2 and 2 with conn=1, so conn=2 means 8-connected foreground.
  - A 4-pixel diamond ring (8-connected ring around a 4-connected hole) gives χ=0 with conn=2. This matches the bit-quad χ8 and our convention.
  - 3D corner-touching voxels give 1 with conn=3 and 2 with conn=1, so conn=3 means 26-connected foreground.
  - A cube with a 1-voxel cavity gives 2 (β0+β2). A square torus gives 0.
  - On 400 random 2D masks, `euler_number(m, 2)` and the numpy bit-quad χ8 both equal β0−β1 from 8/4 CCL.
- Per-component holes. Label the background of `pad(m,1)` with 4-connectivity; labels other than the outside label (`bl[0,0]`) are holes. For each hole take its first raster pixel (via `np.unique(..., return_index)` over raster-ordered hole pixels). The pixel directly above it is foreground and belongs to the enclosing 8-connected component. The reason: the pixel above the topmost hole pixel cannot be background without being 4-adjacent and therefore part of the hole, and the border of a 4-connected background component is 8-connected. A `bincount` then gives holes per instance. This matched per-component brute force (label each component alone and count holes) on 400 masks with 0 mismatches, after fixing a padding bug (pad `m` and then invert, not `~m`). — [local bench](prototypes/bench_mask.py); [holes_fix.py](prototypes/holes_fix.py)
- Timings in µs, single thread. Columns: numpy bit-quad χ8, skimage `euler_number`, cv2-based β0 + per-component β1 (two Spaghetti labelings plus bookkeeping). — [local bench](prototypes/bench_mask.py)
  - 64²: 30 / 51 / 21
  - 256² blob: 243 / 532 / 160
  - 512² blob: 1126 / 1985 / 546
  - 1024² blob: 4298 / 8424 / 2216
  - 1024² noisy: 4006 / 7765 / 8720 (np.unique cost grows with hole pixels)
- 3D timings, nodule-like volumes. cc3d 26-connected foreground vs 6-connected padded background vs scipy 26-connected vs skimage 3D Euler: [local bench](prototypes/bench_mask.py)
  - 32³: 30 / 69 / 235 / 390 µs
  - 64³: 171 / 361 / 1418 / 2556 µs
  - 128³: 1.5 / 2.4 / 9.8 / 19.3 ms
- Betti matching (Stucki, Bürgin, Paetzold, Bauer 2024). Betti matching uses induced matchings of persistence barcodes to match spatially corresponding features. Betti-matching-3D is C++ with a pybind interface and is MIT licensed. It needs five barcode computations per pair. It is optimized relative to Cubical Ripser, with implementation optimizations up to more than 3x, but a VesSAP label volume was a counter-example where Cubical Ripser was faster. GPU is named only as future work (a Ripser++-style approach). — [arXiv 2407.04683](https://arxiv.org/abs/2407.04683); [Betti-Matching-3D repo](https://github.com/nstucki/Betti-Matching-3D)
- Topograph (Lux et al., ICLR 2025) builds a component graph of prediction and label. Its strict topological metric checks homotopy equivalence between the union and the intersection of the prediction/label pair, and it reports up to 5x faster loss computation than persistent-homology methods. — [arXiv 2411.03228](https://arxiv.org/abs/2411.03228v1); [ICLR poster](https://iclr.cc/virtual/2025/poster/29716)

### Inferences
- For binary masks, every persistence pair is born and dies at the same threshold. PH therefore adds nothing over CCL for Betti numbers. Its value is spatial matching of features (Betti matching) and soft or likelihood inputs. Our hole-matching needs (match predicted holes to consensus holes) can be met with CCL hole labels plus IoU matching (Q4) at about 100x lower cost than barcode computation.
- Topograph's "union vs intersection homotopy" metric is conceptually close to our egg–yolk band (intersection = yolk, union = egg). A cheap proxy is to require β(pred ∩ I) and β(pred ∪ U) to be consistent. This is worth noting as a related-work anchor.
- 3D: β0 comes from cc3d 26-connectivity, β2 = (#6-connected background components of the padded volume − 1), and β1 = β0 + β2 − χ. χ comes from the 2×2×2 (256-entry) LUT. The skimage LUT is correct with conn=3, but at 2.5 ms per 64³ it is about 10x slower than the two labelings. A shift-and-OR uint8 code plus `np.bincount` over 256 bins, with skimage's `EULER_COEFS3D_26` coefficients, should be much faster (hypothesis H3).
- Speck robustness (no persistence needed): apply `cc3d.dust` or an area threshold to predicted components and holes below a scale-aware minimum (e.g. max(4 px, 1% of union area)) before computing (β0, β1). Report the unfiltered counts as a separate penalty so the policy cannot hide specks. This is the binary analogue of a persistence threshold.

### Gaps
- I did not retrieve the primary Lin et al. paper on the 3D Euler LUT. skimage cites Ohser/Nagel/Schladitz (2002) and Rivollier (2010) for the coefficients, and the verified skimage LUT is a usable reference implementation.
- I did not extract Betti-matching-3D's absolute runtime numbers (Table II/III values), only the qualitative result.
- No published per-instance (β0, β1)-equal-any-annotator metric was found. It appears novel to this project.

## Q3. Topology-aware segmentation metrics (2023–2026) and speck robustness

### Takeaway
The current family is Betti error, Betti matching error, clDice, Skeleton Recall and Topograph's strict metric. They target tubular or graph-like structures and use PH or skeletons. For blob-like medical objects (nodules, optic disc/cup) with 0–2 holes, an exact per-instance (β0, β1) check with size filtering is sufficient and far cheaper. Use Betti matching error only as an offline eval metric.

### Cited Findings
- Metrics Reloaded recommends clDice for tubular structures (FP2.3/FP3.3) and notes clDice "can be complemented by application-specific connectivity" metrics. — [Metrics Reloaded arXiv 2206.01653](https://arxiv.org/abs/2206.01653)
- Skeleton Recall Loss (Kirchhoff et al., ECCV 2024) uses cheap CPU operations instead of GPU-heavy clDice-style computations. It cuts compute overhead by more than 90%, matches or beats clDice on 4 of 5 datasets, and handles multi-class labels. — [arXiv 2404.03010](https://alphaxiv.org/abs/2404.03010); [ECCV poster](https://eccv2024.ecva.net/virtual/2024/poster/1449)
- Betti matching error counts only unmatched features between images, refining the plain Betti-number error. — [Betti-Matching-3D repo](https://github.com/nstucki/Betti-Matching-3D)
- Topograph's strict topological metric uses union/intersection homotopy equivalence and is up to 5x faster than PH losses. — [arXiv 2411.03228](https://arxiv.org/abs/2411.03228v1)

### Inferences
- For a reward with 256 masks per step, Betti matching (five cubical barcodes per pair, times 2–6 annotators) is the wrong cost class. Per-instance CCL topology costs about 20–550 µs per mask in total (Q2).
- Speck robustness, in descending order of preference:
  1. Area filters scaled to the object (minimum component and hole area = max(a_min, ρ·|U|)).
  2. "Topology at tolerance": compute (β0, β1) on the prediction after a morphological opening/closing of radius τ_px, and also on the raw prediction. Reward the tolerant version and penalize the raw–tolerant difference lightly.
  3. Hole validity requires the hole to overlap the consensus hole region, or to lie outside every annotator's mask (from hole matching, Q4).
  These are hypotheses; no paper establishes thresholds for this setting.

### Gaps
- I found no 2025–2026 paper defining an instance-wise multi-rater topology score. I did not search the newest MICCAI 2026 proceedings.

## Q4. Instance and hole matching: Hungarian vs unique-IoU>0.5 vs greedy

### Takeaway
Build a sparse component-IoU matrix from a joint label histogram (`np.bincount(lab_p * (n_g+1) + lab_g)`). Match with the PQ rule (IoU > 0.5 is unique, no assignment solver) for the "matched" count. Run `scipy.optimize.linear_sum_assignment` only for soft scores at low IoU: it costs 0.5–2.4 µs for k ≤ 20. lap/lapx are unnecessary.

### Cited Findings
- PQ's unique matching theorem (Kirillov et al.): if IoU(p, g) > 0.5 then no other non-overlapping prediction can have IoU > 0.5 with g, so matches are unique and need no assignment step. — [Panoptic Segmentation arXiv 1801.00868](https://arxiv.org/abs/1801.00868)
- Metrics Reloaded (FP5.4) recommends an IoU > 0.5 assignment when predictions cannot overlap and there are no class scores. Otherwise it recommends greedy-by-score or greedy-by-localization, or optimal Hungarian matching. For small structures (FP3.1/FP3.2) it recommends lower localization thresholds and size stratification. — [Metrics Reloaded](https://arxiv.org/abs/2206.01653)
- `linear_sum_assignment` timings, single thread: 0.5 µs (2×2), 0.8 µs (6×6), 2.4 µs (20×20), 80 µs (100×100). — [local bench](prototypes/bench_mask.py)

### Inferences
- The components of one predicted mask cannot overlap each other, so the PQ uniqueness theorem applies to pred-vs-annotator component matching. Hungarian only matters for partial-credit matching below 0.5, which is useful for tiny components where 1 px shifts break IoU > 0.5. Following Metrics Reloaded FP3.1, use a lower threshold such as 0.3, or a boundary-tolerant IoU, for components under ~200 px. Below 0.5, uniqueness is lost and Hungarian is required.
- Holes: treat consensus holes (holes present in a majority of annotators, matched across annotators by IoU) as instances. Match predicted holes by IoU with Hungarian, since there are few. For tiny holes, use "hole centroid inside a consensus hole dilated by τ" as the localization criterion; Metrics Reloaded lists point-inside-mask style criteria among its localization options.
- Joint-histogram IoU (one `bincount` over n_p × n_g bins) is O(pixels) and avoids per-pair masks.

### Gaps
- I did not benchmark lap/lapx. At k ≤ 20, scipy's sub-3 µs leaves nothing to gain.

## Q5. Boundary metrics, tolerance, distance transforms, and scale-aware band formulas

### Takeaway
The current band (1 − (miss_core + spill)/|U|) is scale-biased. A local demo shows it is harsh on tiny objects and lenient on large ones: a 3 px spill scores 0.00 at r=4 and 0.93 at r=80. The evidence (Boundary IoU, Metrics Reloaded, Nikolov surface DSC) says annotator boundary disagreement is roughly fixed in pixels or mm, not proportional to size. So normalize errors by perimeter × tolerance (F1), or weight them by distance and normalize by perimeter (F2). Set τ from the per-image inter-annotator surface-distance spread. Do not use a √area-scaled tolerance (F4).

### Cited Findings
- Boundary IoU (Cheng et al., CVPR 2021) is |(G_d ∩ G) ∩ (P_d ∩ P)| / |(G_d ∩ G) ∪ (P_d ∩ P)|, where G_d and P_d are pixels within distance d of each contour. It is symmetric, and per the authors it does not over-penalize small objects. Mask IoU, by contrast, "tolerate[s] a larger number of misclassified pixels per each unit of contour length for a larger object". — [Boundary IoU arXiv 2103.16562](https://arxiv.org/abs/2103.16562)
- From the same paper: the distance between two annotators' contours "seldomly exceeds 1% of the image diagonal, irrespective of the object size", confirmed on LVIS double annotations, where "boundary ambiguity is fixed and independent of objects area". The authors choose d so that median expert-vs-expert Boundary IoU exceeds 0.9: d = 2% of the image diagonal (about 15 px) on COCO/ADE20k and 0.5% on Cityscapes. Their advice is that "annotation consistency sets the lower bound on d". They also show that F-measure (contours within d count as matched) is "nearly discontinuous", which makes rewards high-variance. — [Boundary IoU](https://arxiv.org/abs/2103.16562)
- Metrics Reloaded: NSD should be used when tolerance for imprecise annotations is wanted, with τ "set according to the inter-rater variability or, if not available, heuristics". Boundary IoU suits cases where contour errors matter, and its d can also be set from inter-rater variability. Mask IoU "over-penalizes small structures". Boundary IoU can reach 1.0 for an imperfect prediction. — [Metrics Reloaded](https://arxiv.org/abs/2206.01653)
- Surface DSC (Nikolov et al., DeepMind) uses organ-specific τ "defined as the 95th percentile of the distances" between surfaces of independent expert segmentations. It counts both directed surface distances within τ and normalizes by combined surface area. — [arXiv 1809.04430](https://arxiv.org/abs/1809.04430); [JMIR 2021 version](https://research.google/pubs/clinically-applicable-segmentation-of-head-and-neck-anatomy-for-radiotherapy-deep-learning-algorithm-development-and-validation-study/)
- Distance-transform cost, µs single thread. Columns: scipy `distance_transform_edt`, seung-lab `edt` (parallel=1), `cv2.distanceTransform(DIST_L2, DIST_MASK_PRECISE)`. — [local bench](prototypes/bench_mask.py)
  - 64²: 73 / 57 / 20
  - 256² blob: 1117 / 501 / 341
  - 512² blob: 5878 / 2090 / 1488
  - 1024² blob: 31851 / 9691 / 7410
  - 1024² noisy: 43690 / 12523 / 12536
- `edt` uses Felzenszwalb–Huttenlocher parabolic envelopes with Meijster/Saito–Toriwaki axis decomposition and supports anisotropy, `parallel`, and `sdf`/`sdfsq` signed distances. — [edt repo](https://github.com/seung-lab/euclidean-distance-transform-3d)
- Local band-formula demo. Three annotators are concentric disks of radii r−1, r, r+1, so the band is 2 px wide. The prediction spills or misses by δ px beyond U or inside I. τ = 2 px. Scores are listed for r = 4 / 10 / 40 / 80. — [local band demo](prototypes/band_demo.py)

  | Formula | spill δ=1 | spill δ=3 | miss δ=3 |
  |---|---|---|---|
  | F0 current, 1 − (miss+spill)/\|U\| | 0.60 / 0.83 / 0.95 / 0.98 | 0.00 / 0.40 / 0.85 / 0.93 | — / 0.62 / 0.86 / 0.93 |
  | F1 perimeter-normalized, 1 − (miss+spill)/(τ·\|∂U\|) | 0.56 / 0.62 / 0.60 / 0.62 | 0.00 at every r | — / 0.14 / 0.00 / 0.00 |
  | F2 distance-weighted, Σ min(1, d/τ) over error pixels / \|∂U\| | 0.56 / 0.62 / 0.58 / 0.61 | 0.00 at every r | — / 0.00 / 0.00 / 0.00 |
  | F3 NSD-like (pred boundary within τ of the band) | 1.00 at every r | 0.00 / 0.07 / 0.09 / 0.08 | 0.00 at every r |
  | F4 tolerance k·√\|U\| (k=0.1) | 0.60 / 1.00 / 1.00 / 1.00 | 0.00 / 0.64 / 1.00 / 1.00 | — / 0.80 / 1.00 / 1.00 |

  - F1 and F2 are close to scale-invariant for r ≥ 6–10 px; at r=4, discretization adds noise.
  - F3 is step-like.
  - F4 inverts the bias: tolerance grows with size, so large-object errors vanish.

### Inferences
- **Candidate band formulas, ranked:**
  - **(A) Perimeter-normalized band (F1).** Formula: band = 1 − (|pred \ U| + |I \ pred|) / (τ · |∂U|). Error area divided by perimeter is the mean normal excursion in px, so the score is about 1 − mean_excursion/τ. This is scale-free by construction. Support: Boundary IoU's fixed-pixel ambiguity finding and Metrics Reloaded's small-structure guidance. Cheapest option: one 3×3 erosion for |∂U| and no EDT.
  - **(B) Distance-weighted band (F2).** Each spill pixel is weighted min(1, d(x, U)/τ) and each core miss min(1, d(x, ∂I)/τ), normalized by |∂U|. This forgives sub-τ excursions partly, so it is smoother near zero error, and it penalizes far spill harder in shape. Needs two EDTs per annotator set. These can be precomputed once per image, because U and I are fixed across the GRPO group. Per-rollout cost is then only a gather, which removes EDT from the hot loop. Support: surface-DSC tolerance logic (Nikolov).
  - **(C) Boundary-in-tolerance-band (NSD-like, F3)** as a separate, low-weight term. Use the fraction of predicted boundary pixels within τ of the egg–yolk band, plus the fraction of band "centerline" covered. It is step-like (Boundary IoU paper critique), so do not use it alone as a GRPO reward. A soft version, mean of exp(−d/τ) over predicted boundary pixels, fixes the step.
  - **(D) Rejected: tolerance ∝ √area (F4).** It contradicts the evidence that inter-annotator boundary disagreement does not grow with object size. Use it only if a dataset shows otherwise (testable, H5).
- **Setting τ:** per image, τ = max(τ_min, P95 of symmetric surface distances between annotator pairs), following Nikolov. Add a floor τ_min of about 1–1.5 px, because tiny LIDC nodules with near-identical annotations would otherwise give τ≈0. Optionally cap τ at a dataset-level P95. This makes the egg–yolk band and the tolerance consistent: the band already encodes the disagreement, and τ only adds slack beyond it.
- Normalizing by |∂U| fixes the original symptom (band ≈ 0 for ~100 px nodules). With r≈5.6 px, a 1 px all-round spill is about 35 px of error, which F0 divides by |U|≈100 and F1 divides by τ·|∂U|≈2·35. F1 thus reads "1 px of excursion out of a 2 px tolerance" at every scale.
- Hot-loop budget at 256² per mask: CCL β0/β1 ≈ 0.16 ms, IoU ≈ 0.05 ms, band F1 ≈ 0.05 ms (erosion of U is precomputed per image), plus F2 gathers. With EDTs precomputed per image the total is under about 0.3 ms per mask, or about 77 ms per 256 masks on one core and about 10 ms on 8 cores. Even 1024² blob crops at about 3–5 ms per mask stay under 0.2 s on 8 cores. Rollout-side EDT (needed only for HD95/NSD from the prediction's side) is the most expensive item; use cv2.distanceTransform for it.

### Gaps
- Not yet measured: LIDC/RIGA inter-annotator P95 surface distances as a function of nodule or disc size. This is the decisive data for choosing (A)/(B) vs (D), and it is local work (H5).
- I did not fetch MONAI's or DeepMind's `surface-distance` implementation details (crack/edge surface elements vs pixel centers). Surface DSC "rounds the tolerance to the nearest neighbour" of quantized distances (Nikolov appendix), which matters at τ≈1 px.

## Q6. Multi-rater evaluation metrics for ambiguous segmentation

### Takeaway
For the eval set, keep the following:
- IoU/Dice vs majority.
- QUBIQ-style threshold-averaged Dice vs the soft consensus.
- Max-over-annotators Dice and its mean over annotators (D_max-like), so a prediction that matches any expert is credited.
- Combined sensitivity vs the union.
- The new band score (F1/F2).
- Per-instance topology-match-any.

GED and HM-IoU are distribution metrics for samplers. They are only relevant if several rollouts per image are treated as a distribution: GED rewards diversity regardless of fidelity. STAPLE is a consensus-building tool, not a metric.

### Cited Findings
- QUBIQ binarizes the averaged expert map and the prediction at 0.1, 0.2, …, 0.9 and averages the Dice scores over thresholds. — [QUBIQ arXiv 2405.18435](https://arxiv.org/abs/2405.18435)
- QUBIQ found that models trained on averaged annotation maps beat those trained on majority votes or random annotator selection. — [QUBIQ](https://arxiv.org/abs/2405.18435)
- GED: D²_GED = 2E[d(S,Y)] − E[d(S,S′)] − E[d(Y,Y′)] with d = 1 − IoU. Rahman et al. (CIMD, CVPR 2023) state that GED "reward[s] sample diversity regardless of the generated samples' fidelity". They propose the CI score, combining combined sensitivity S_c (TP/(TP+FN) between unions of predictions and of ground truths), D_max (mean over ground truths of the max Dice to any prediction) and diversity agreement D_a. On LIDC: PHi-Seg GED 0.270 / CI 0.736 / D_max 0.904; CIMD GED 0.321 / CI 0.759 / D_max 0.915. — [arXiv 2304.04745](https://arxiv.org/html/2304.04745). The CI formula as extracted (3·S_c·D_max·D_a/(S_c+D_max+D_a)) is not a true harmonic mean; check it against the paper PDF before implementing.
- HM-IoU builds a C×C IoU matrix between predicted samples and annotations and uses Hungarian matching to maximize total IoU. — [ambiguous segmentation literature, e.g. arXiv 2511.06857](https://arxiv.org/pdf/2511.06857)
- LIDC-IDRI: four radiologists segmented nodules ≥3 mm. One cited study reports inter-observer Dice of 0.77 ± 0.09, and another found "very high inter-observer variation" in size metrics for 127 nodules marked by all four readers (2.03–29.4 mm). — [Cornell VIA APRAMB07](https://www.via.cornell.edu/ece5470/home/publications/APRAMB07.pdf). These figures came from search snippets; I did not read the full source.

### Inferences
- For a single-output VLM policy, the right analogues are "match-any" metrics: max-over-annotators Dice/IoU, and topology equal to any annotator. Within-group GRPO rollouts could be scored by HM-IoU against the annotators as a diagnostic of whether the policy's diversity tracks expert diversity. That is an eval metric, not a reward.
- STAPLE gives a probabilistic consensus. It could replace majority voting as the IoU target, but it adds an EM fit per image (precomputable) and offers no clear advantage over mean-map QUBIQ Dice. I did not research STAPLE further; flagged as an assumption.

### Gaps
- I did not verify the full HM-IoU protocol (sample count, handling of empty masks) from the original paper.

## Q7. 3D: bit-cube Euler LUT, voxel surfaces, cc3d 26/6

### Takeaway
3D is feasible on CPU. cc3d 26-connected plus padded 6-connected background labeling costs about 0.5 ms at 64³ and 4 ms at 128³. χ comes from a 256-entry 2×2×2 LUT (skimage's `EULER_COEFS3D_26`, verified semantics). β1 (tunnels) = β0 + β2 − χ.

### Cited Findings
- skimage implements 3D Euler with 2×2×2 configurations (256 bins) and Ohser coefficients `EULER_COEFS3D_26`. For conn=1 (6-connectivity) it reverses the LUT, and it does not support conn=2 (18). Verified locally: cube with cavity χ=2, square torus χ=0, corner-touching voxels χ=1 with conn=3. — [skimage source](https://github.com/scikit-image/scikit-image/blob/main/skimage/measure/_regionprops_utils.py); [local bench](prototypes/bench_mask.py)
- cc3d supports 6/18/26 connectivity and `contacts()` for surface contact area (anisotropy-aware). — [cc3d repo](https://github.com/seung-lab/connected-components-3d)
- Timings: cc3d 26-connected 30 / 171 / 1489 µs at 32³ / 64³ / 128³, vs scipy 26-connected 235 / 1418 / 9774 µs. — [local bench](prototypes/bench_mask.py)

### Inferences
- For a certificate, a voxel surface can be extracted exactly as the set of boundary faces between foreground and 6-adjacent background voxels. The face count gives surface area, and it is a closed 2-manifold only for well-composed volumes. A certificate can cross-check χ(surface)/2 against the bit-cube χ where valid. This is a hypothesis; I did not research marching-cubes or face-mesh libraries in this pass.

### Gaps
- I did not retrieve the primary 2×2×2 LUT literature (Lin/Sha/Lee-Kashyap-Chu), the Betti-matching-3D CUDA status beyond "future work", or voxel-surface library choices.

## Candidate hypotheses to test locally

- **H1 (CCL choice).** On 10k real crops (LIDC slices, RIGA disc/cup, VLM-rasterized polygons; 64–1024 px), cv2 Spaghetti ≤ cc3d ≤ scipy in wall time per mask on the Linux training host, with 0 partition mismatches against a brute-force flood fill (BFS with explicit 8/4 neighborhoods). Also measure `multiprocessing` vs a thread pool. cv2 releases the GIL, so 8 threads may suffice without process overhead.
- **H2 (per-component holes).** The raster-first-pixel hole assignment equals per-component brute force on 10k masks, including adversarial cases: holes touching the crop edge, nested islands inside holes, 1-px diagonal-gap rings, and checkerboards. Bit-quad χ8 = Σ(1 − holes_k) must hold as a certificate invariant.
- **H3 (3D Euler speed).** A shift-and-OR uint8 2×2×2 code with `np.bincount` and the skimage LUT matches `euler_number(v, 3)` exactly and is ≥5x faster than skimage at 64³. β1 = β0 + β2 − χ matches known tori and handles.
- **H4 (band scale invariance on real data).** On LIDC (4 readers) and RIGA (6 readers), apply controlled perturbations (dilate/erode by 1–3 px, translate 1–2 px, polygon simplification) to an annotator mask. Score it against the remaining annotators. Regress score on log(object area). Prediction: F1/F2 slope ≈ 0, F0 slope > 0 (biased toward large objects), F4 slope > 0. Also check the GRPO-relevant property: within-group score variance should not collapse for small objects.
- **H5 (what τ should scale with).** Compute leave-one-annotator-out P95 symmetric surface distance on LIDC/RIGA vs object size. If P95 is flat in size (as in Boundary IoU/LVIS), use a fixed-px or per-image τ (A/B). If it grows like √area, revisit F4.
- **H6 (topology speck filter).** Sweep the minimum component/hole area ρ·|U| over ρ ∈ {0, 0.5%, 1%, 2%} with floors of 2–6 px. Measure the fraction of annotator masks whose (β0, β1) changes; it should be close to 0 at the chosen ρ. Measure the false "topology match" rate on perturbed predictions.
- **H7 (matching threshold for tiny components).** For components under 200 px, compare IoU>0.5 unique matching, Hungarian at IoU>0.3, and centroid-in-dilated-mask matching by agreement with human-judged correspondence on 200 LIDC multi-nodule crops.
- **H8 (end-to-end budget).** Score 256 masks at the realistic size mix in under 100 ms on 8 cores, with per-image precomputation (U, I, ∂U, EDTs of U and I, consensus holes) cached across the GRPO group.

## Recommendation

- **2D CCL:** `cv2.connectedComponentsWithAlgorithm(m_u8, 8, CV_32S, CCL_SPAGHETTI)` for the foreground and connectivity 4 (Spaghetti4C) for the padded background. Keep scipy as the brute-force-adjacent reference in tests. Note the OpenCV 5.0 Python API: `ccltype=` on `connectedComponents` fails. Add opencv-python-headless to core dependencies; it is currently only under `models`.
- **Per-instance topology:** exact (β0, β1-per-component) via the background-label, pixel-above-first-hole-pixel assignment ([holes_fix.py](prototypes/holes_fix.py)), with area-based speck filters scaled to |U|. Keep the bit-quad χ8 identity as the certificate check. Skip persistent homology and Betti matching in the reward; Betti matching error is optional for offline eval.
- **Matching:** joint-histogram IoU matrix, PQ unique matching at IoU>0.5 for counts, and `scipy.optimize.linear_sum_assignment` for soft or low-threshold matching and for holes. No lap/lapx.
- **Band:** replace F0 with (A), the perimeter-normalized band 1 − (spill + core_miss)/(τ·|∂U|) clipped to [0, 1], or (B), its distance-weighted variant. τ is per image: the P95 inter-annotator surface distance with a floor of about 1–1.5 px. Precompute U, I, ∂U and the EDTs once per image per GRPO group. Add (C), a soft NSD-in-band term, only at low weight. Validate with H4/H5 before freezing.
- **Distance transforms:** `cv2.distanceTransform(DIST_L2, DIST_MASK_PRECISE)` for 2D (3–4x faster than scipy). Use `edt` for 3D or anisotropic data and for signed distances.
- **Eval set:** IoU/Dice vs majority, QUBIQ threshold-averaged Dice vs the mean map, max-over-annotators Dice (D_max-like), combined sensitivity vs the union, band F1/F2, per-instance topology-match-any, Boundary IoU with d set from inter-annotator agreement, and HD95 for reporting. Use GED and HM-IoU only for multi-sample diagnostics.
- **3D:** cc3d 26-connected foreground and 6-connected padded background, a fast 256-entry 2×2×2 LUT for χ, and β1 = β0 + β2 − χ. GPU CCL (BKE) is a later CUDA milestone, not needed for the CPU budget.
