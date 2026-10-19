# Polygon side of the exact geometry reward checker: algorithms, libraries, Python/C++ trade-off

Scope: parse, then validity grading (never repair), then exact rasterization with explicit fill rule and pixel-centre convention, plus the reverse path for SFT targets (mask, then exact pixel-edge rings, then topology-preserving simplification). Integer vertices on a 0–1000 grid are mapped onto crops of S ≈ 100–1000 px. Budget: about 256 rollouts per GRPO step, 1–10 rings each, about 10–200 vertices per ring, under 1 s total on 8 cores.

Evidence comes from two places. The first is source and doc reading, cited inline. The second is **local experiments run for this note** on an Apple M3 (8 cores, 4P+4E) with Python 3.13, numpy 2.5.4, Pillow 13.0.0, shapely 2.2.0 / GEOS 3.14.1, scikit-image 0.26.0, OpenCV 5.0.0, matplotlib 3.11.2, rasterio 1.5.2 / GDAL 3.12.2, numba 0.68.0 and potracer (pure-Python potrace port). The scripts are saved next to this note in `prototypes/poly_*.py`. They are labelled **[local experiment]**, and they are first-party measurements, not literature.

Exact frame used throughout the experiments. A vertex v on the 0–1000 grid maps to pixel coordinate x = v·S/1000. The centre of pixel (i, j) is at (i+½, j+½). Multiplying everything by 2000 makes every quantity an integer: vertex V = 2·S·v, centre C = 1000·(2i+1). With S ≤ 1000 every product in the predicates stays below about 10¹³. That is far inside int64, and also inside the exactly representable float64 integers (2⁵³ ≈ 9·10¹⁵).

---

## 1. Exact polygon rasterization with fill rules (even-odd, nonzero, pixel-centre convention)

### Takeaway
No off-the-shelf Python rasterizer gives a **nonzero** mask. None of them documents its tie-break for pixel centres that lie exactly on an edge.
- **Exact at strictly interior/exterior centres:** skimage, matplotlib, shapely and rasterio/GDAL are exact **even-odd** pixel-centre samplers there. They had 0 off-boundary mismatches on 1,200 polygons.
- **Not centre samplers:** **PIL and OpenCV** over-fill boundary-touched pixels. They had about 135 to 175 wrong pixels per polygon at S=250.
- **Speed:** the exact ones are either slow (4–45 ms per polygon) or are GDAL (about 0.2–0.3 ms).

A roughly 30-line integer scanline does better. It accumulates per row a signed difference array of edge crossings (the font-rs style "accumulate then prefix-sum", but point-sampled instead of area coverage), and its prefix sum is the winding number per pixel. That gives provably exact even-odd and nonzero masks in a single pass. In numba it runs in about 75–80 µs per polygon at S=512 and matched brute force bit for bit.

### Cited Findings
- **Pillow `ImageDraw.polygon`:**
  - The C fill (`polygon_generic`) samples scanlines at **integer y with no +0.5 offset**. It computes x-intersections in float, fills spans pairwise between sorted intersections (even-odd style pairing), and rounds span ends with asymmetric `ROUND_UP` / `ROUND_DOWN` macros. Horizontal edges are drawn directly as hlines, and consecutive horizontal segments are merged. — [Pillow Draw.c](https://github.com/python-pillow/Pillow/blob/main/src/libImaging/Draw.c)
  - [local experiment] Pillow 13.0.0 produced **identical masks with and without `outline=1`**. That means the fill already includes boundary pixels, so the repo's `rasterize(..., outline=1)` call in `src/objdel/checker/polygon_io.py` is redundant.
  - [local experiment] Against an exact pixel-centre reference, PIL differed **off the boundary on 295/296 non-self-overlapping polygons** (S=250), averaging about 119 px per polygon. At S=1000 it was about 694 px per polygon. For star polygons the split was 27,020 extra pixels against 1,400 missing pixels over 200 polygons, so it is mostly over-fill.
  - [local experiment] On self-overlapping polygons its error against even-odd was smaller than against nonzero (202,723 vs 276,281 px), which is consistent with even-odd pairing.
  - [local experiment] Speed: 28–104 µs per polygon at S=512 for n=16–200.
- **OpenCV `fillPoly`:**
  - It uses an edge table sorted by (y0, x, dx), integer scanlines, and fixed-point x with `XY_SHIFT` fractional bits rounded with `delta = XY_ONE-1`. Spans are inclusive of the right end (`xr+1`). The fill toggles `draw ^= 1` per crossing. The fetch summary called this "non-zero", but a parity toggle is even-odd. — [OpenCV drawing.cpp](https://github.com/opencv/opencv/blob/4.x/modules/imgproc/src/drawing.cpp)
  - [local experiment] Even with `shift=8` sub-pixel input and a −0.5 centre shift, it differed off-boundary on 295/296 polygons at S=250. Star polygons got 35,180 extra pixels against 165 missing over 200 polygons, so it is a conservative over-fill.
  - [local experiment] Fastest of all at 11–28 µs per polygon, but not centre-exact.
- **scikit-image `draw.polygon` / `polygon2mask`:**
  - It loops over every integer (r, c) in the bounding box and calls `point_in_polygon` (O'Rourke's crossing test with separate left and right ray counts). Pixel centres are therefore at **integer coordinates**. It returns INSIDE, EDGE or VERTEX (non-zero values), so **centres on edges or vertices are counted as inside**, and the vertex test uses a 1e-12 tolerance. The work is O(bbox pixels × n). — [skimage _draw.pyx](https://github.com/scikit-image/scikit-image/blob/main/src/_skimage2/draw/_draw.pyx), [skimage geometry.pyx](https://github.com/scikit-image/scikit-image/blob/main/src/_skimage2/_shared/geometry.pyx)
  - [local experiment] 0 off-boundary mismatches vs exact even-odd on 1,200 polygons (S=250 and S=1000). Speed: 4.4 ms (n=16), 10.5 ms (n=64) and 22.8 ms (n=200) per polygon at S=512.
- **matplotlib `Path.contains_points`:**
  - `point_in_path_impl` is an even-odd crossing test (`subpath_flag[i] ^= 1`). It has **no nonzero option**. A non-zero `radius` contours or strokes the path through Agg, which makes it inexact. — [matplotlib _path.h](https://github.com/matplotlib/matplotlib/blob/main/src/_path.h)
  - [local experiment] 0 off-boundary mismatches vs even-odd. 4–45 ms per polygon at S=512.
- **GDAL / rasterio `rasterize`:**
  - The default burns pixels whose **centre** is inside. `ALL_TOUCHED` burns every touched pixel. The docs say behaviour is "unspecified" when the polygon just touches a pixel centre. — [gdal_rasterize docs](https://gdal.org/programs/gdal_rasterize.html)
  - [local experiment] 0 off-boundary mismatches vs even-odd, including on self-intersecting input. 190–310 µs per polygon at S=512, the fastest of the exact libraries.
- **shapely 2.x `contains_xy` on pixel centres:**
  - [local experiment] 0 off-boundary mismatches vs even-odd (also on invalid polygons, empirically). 5–9.6 ms per polygon at S=512 because it is a point query over S² centres. `contains` excludes the boundary by definition.
- **pycocotools `frPoly`:**
  - `rleFrPoly` upsamples vertices ×5 with `(int)(scale*x+.5)` rounding, digitizes the boundary as DDA lines, and keeps y-boundary crossings at pixel centres (`(xd+.5)/scale-.5`, integer centres). It then builds the RLE by alternating runs, which is parity, so even-odd. Edges are approximated at 1/5 px, so the result is **not exact**. — [cocoapi maskApi.c](https://github.com/cocodataset/cocoapi/blob/master/common/maskApi.c)
- **Coverage (anti-aliasing) rasterizers:** tiny-skia exposes `FillRule::{Winding, EvenOdd}`, and AA is a per-paint flag. — [tiny-skia docs](https://docs.rs/tiny-skia/)
- **Vello sparse strips:** vello_cpu / vello_hybrid use a "sparse strips" CPU/hybrid design based on Levien's "Potato" design, described in Laurenz Stampfl's master's thesis. — [GeekNews summary](https://news.hada.io/topic?id=24300), thesis PDF at [github.com/LaurenzV/master-thesis](https://github.com/LaurenzV/master-thesis/blob/main/main.pdf), [vello_cpu source](https://docs.rs/crate/vello_cpu/0.1.0/source/src/lib.rs). Vello has commits simplifying both NonZero and EvenOdd strip logic. — [vello commit](https://skia.googlesource.com/external/github.com/linebender/vello/+/refs/heads/velo_common%5E%21/)
- **Exact integer scanline written for this note** [local experiment]:
  - **How it works:** for each non-horizontal edge, consider every row j whose centre Y = 1000(2j+1) satisfies min(y) ≤ Y < max(y) (half-open). With num = (Y−y0)·dx + x0·dy and dy > 0, the edge lies to the right of all pixels i ≤ i_max, where i_max = ⌊(⌊(num−1)/(1000·dy)⌋ − 1)/2⌋. The routine adds ±1 (edge direction) at column 0 and subtracts it at i_max+1, then takes a row prefix sum. The result is the winding number w per pixel: even-odd is w odd, nonzero is w ≠ 0.
  - **Validation:** **0 mismatches** against an O(S²·n) brute-force winding reference (same half-open rule) on 500 polygons at S=250. The set mixed star, random self-intersecting and rectilinear polygons.
  - **Speed:** 75–80 µs per polygon at S=512, nearly independent of n because the S² prefix sum dominates. It can be limited to the bounding box.
- [local experiment] **On-edge pixel centres are common** on a 0–1000 grid: about 31 per polygon at S=250 (18,674 over 600 polygons) and about 47 per polygon at S=1000 (27,982 over 600). Every library broke these ties differently: on-boundary mismatch counts against the reference were 2,892–17,747 px. A tie-break rule is therefore part of the spec, not a corner case.

### Inferences
- **Our reference's tie-break is the "top-left" rule.** It counts a crossing only when the centre is strictly left of the edge, with half-open y. So a centre exactly on a left edge or a top edge (y-down) is in, and on a right or bottom edge is out. This is the rasterization rule GPUs use. Two polygons sharing an edge then partition the pixels with no double counting or gaps, which is the property you want for a reward. The D3D/GL top-left rule docs were not fetched in this session.
- **The nonzero XOR even-odd self-overlap flag is incomplete.** eo ≠ nz only where |w| ≥ 2. A bow-tie or figure-8 has lobes with w = +1 and w = −1, so eo == nz everywhere and it would pass. The single-pass winding buffer also gives `any(w<0) & any(w>0)` within one ring, which flags inverted lobes. Use both flags.
- **Float64 vs int64:** float64 is exact for this frame, but int64 integer division makes the i_max computation exact without reasoning about rounding. In CUDA, use int64, or int32 if S ≤ about 200.
- **Anti-aliasing rasterizers are the wrong tool for the reward mask.** Area coverage thresholded at 0.5 is not centre sampling, and fixed-point edge snapping (26.6 / 16.16 formats in Skia and FreeType) cannot represent every v·S/1000 crossing exactly. Their accumulate-then-prefix-sum structure is still the right pattern to borrow. The specific fixed-point formats of Skia and FreeType were not re-verified this session.
- **Matching an existing library:** if a library must be matched, for example for a GDAL-based evaluation, rasterio is the only fast, exact even-odd one, and its on-edge behaviour is unspecified.

### Gaps
- PIL's sampling semantics depend on the Pillow version. Polygon fill changed several times in 9.x–10.x, and only 13.0.0 was tested.
- The non-AA paths of tiny-skia, Blend2D and zeno were not benchmarked. No Python bindings were checked for them.
- The Vello thesis numbers were not read. The thesis PDF was not fetched, only the summary page.
- The font-rs blog (Medium) returned 403. Its accumulation-buffer description is cited from domain knowledge, not fetched.

---

## 2. Self-intersection / validity detection (graded, never repaired)

### Takeaway
On a 0–1000 integer grid every orientation determinant is at most about 2·10⁶. Plain int64 (even int32) arithmetic is therefore exact, and Shewchuk-style adaptive predicates are unnecessary.

GEOS/shapely is fast (2–21 µs per polygon) but is the wrong grader:
- it reports **only the first error**;
- it **ignores winding**;
- it accepts repeated vertices;
- it treats a hole touching the shell at a single point as valid.

A custom single pass over all edges of all rings, with ring ids, classifies every pair: proper crossing, touch, collinear overlap, adjacent-edge spike, ring–ring crossing or touch. Each defect becomes its own graded flag. Then do hole containment and winding by exact point-in-polygon and shoelace sign. For n ≤ about 2,000 edges, compiled O(n²) with bounding-box rejection costs about 0.7 µs (n=16) to 117 µs (n=200, worst case with no crossings). The current pure-Python O(n²) test costs 43 ms at n=200, which alone could take about 30 s per GRPO step.

### Cited Findings
- **GEOS IsValidOp:** follows OGC SFS. `setSelfTouchingRingFormingHoleValid()` (default false) makes inverted shells and exverted holes valid. The docs do not state the algorithm or its complexity. — [GEOS IsValidOp doxygen](https://libgeos.org/doxygen/classgeos_1_1operation_1_1valid_1_1IsValidOp.html)
- [local experiment] **What `shapely.is_valid_reason` (GEOS 3.14.1) reports** for each defect class:

  | Case | GEOS result |
  |---|---|
  | bow-tie | `Self-intersection[5 5]` |
  | vertex pinch (figure-8) | `Ring Self-intersection[5 5]` |
  | vertex touching an edge interior | `Ring Self-intersection` |
  | collinear spike (backtrack) | `Ring Self-intersection` |
  | non-adjacent collinear overlap | `Ring Self-intersection` |
  | all-collinear (zero area) | `Self-intersection` |
  | hole outside shell | `Hole lies outside shell` |
  | hole crossing shell | `Self-intersection` |
  | hole sharing an edge with shell | `Self-intersection` |
  | nested holes | `Holes are nested` |
  | hole cutting the interior | `Interior is disconnected` |
  | two crossing holes | `Self-intersection` |
  | inverted shell | `Ring Self-intersection` |
  | **clockwise shell** | **valid** |
  | **consecutive repeated vertex** | **valid** |
  | **hole touching shell at one point** | **valid** |

  The returned **area of an invalid bow-tie is 0.0**. Area is meaningless on invalid input.
- [local experiment] **Speed:** `shapely.is_valid` takes 2.6 µs (n=16), 7.0 µs (n=64) and 20.3 µs (n=200). `is_valid_reason` is about the same.
- **shapely `make_valid`:** two methods. "linework" nodes the rings and rebuilds polygons with an **even-odd** strategy. "structure" assumes shells and holes are labelled, then unions shells and subtracts holes. Both are repair semantics and may change dimension. — [shapely make_valid docs](https://shapely.readthedocs.io/en/stable/reference/shapely.make_valid.html)
- **Boost.Geometry `is_valid(geom, message)`:** OGC validity. Polygon and multipolygon checks are **worst-case quadratic**. The documentation example message is "Geometry has disconnected interior". — [Boost.Geometry is_valid docs](https://www.boost.org/doc/libs/release/libs/geometry/doc/html/geometry/reference/algorithms/is_valid/is_valid_2_with_message.html)
- **Clipper2** (BSL-1.0, int64 `Paths64` or double `PathsD`, NonZero / EvenOdd fill rules, `SimplifyPaths`) is a clipping and offsetting library. It is a polygon validator only indirectly, through union-and-compare. — [Clipper2 GitHub](https://github.com/AngusJohnson/Clipper2). The `PointInPolygon` page confirms the result is IsInside, IsOutside or IsOn but does not state the rule. — [Clipper2 PointInPolygon](https://angusj.com/clipper2/Docs/Units/Clipper/Functions/PointInPolygon.htm)
- **Classic algorithms** (bibliographic, from domain knowledge, not re-fetched this session):
  - Shamos & Hoey 1976, "Geometric intersection problems": an O(n log n) sweep that **decides** whether any pair of n segments intersects. — [doi:10.1109/SFCS.1976.16](https://doi.org/10.1109/SFCS.1976.16)
  - Bentley & Ottmann 1979: **reports** all k intersections in O((n+k) log n). — [doi:10.1109/TC.1979.1675432](https://doi.org/10.1109/TC.1979.1675432)
  - Shewchuk 1997: adaptive-precision exact orientation and incircle predicates for floating-point input. — [doi:10.1007/PL00009321](https://doi.org/10.1007/PL00009321), [Robust predicates page](https://www.cs.cmu.edu/~quake/robust.html)
- [local experiment] **Exact O(n²) int64 segment classifier in numba** (proper crossing via strict sign tests; touch or collinear overlap via zero orientation plus on-segment test; adjacent edges checked only for collinear backtracking): 0.7 µs (n=16), 3.6 µs (n=64), 110–117 µs (n=200). On spikes, closing-vertex spikes, 8-connected pinch rings and redundant collinear vertices it agrees with the repo's `ring_problems`.
- [local experiment] **The repo's pure-Python `ring_problems`** takes 0.25 ms (n=16), 4.3–4.7 ms (n=64) and **41–43 ms (n=200)**. At 256 rollouts × 3 rings × n=200 that is about 33 s per step serially.

### Inferences
- **One pass covers every defect class.** Put all edges of all rings into one array tagged (ring id, edge index). For each pair that survives bounding-box or sweep pruning, compute the four orientation signs and classify:

  | Pair | Classes |
  |---|---|
  | same ring, adjacent | collinear backtrack (spike) |
  | same ring, non-adjacent | proper cross / vertex-touch / collinear overlap |
  | different rings | ring-crossing / ring-touch / ring-collinear-overlap |

  Then, only if no ring–ring crossing exists, test containment with one exact point-in-polygon per ring pair. Pick a vertex of A that is not on B's boundary (or a doubled-coordinate edge midpoint). This yields:
  - "hole outside its outer ring";
  - nesting depth, where an odd depth means it should be a hole;
  - "holes nested" and "outer inside outer";
  - winding as the sign of the exact shoelace sum, checked against the declared role;
  - zero area as a shoelace of 0.

  Raster-level flags (eo≠nz, mixed-sign w) come free from §1's winding buffer.
- **Algorithm choice:** for at most about 2,000 edges per rollout, brute O(n²) with bounding-box rejection is simpler and fast enough (about 0.1 ms at n=200 in compiled code). Bentley–Ottmann's constant factors and degeneracy handling (collinear overlaps, vertical segments) only pay off for much larger n. If a sweep is used, Shamos–Hoey-style early exit cannot be used for grading, because the grader needs counts per class, not a yes/no answer.
- **Degeneracy handling:** a sweep for all pairs must handle coincident endpoints and collinear overlaps explicitly. Brute force with exact integer signs has no such pitfalls. A uniform grid (bucket edges by cells of the 1000-grid) gives near-linear expected time with trivially correct degeneracy handling, if n ever grows.
- **GEOS as a cross-check, not the grader.** Use it as an oracle in tests: any rollout we call fully valid must be `is_valid` in GEOS, after removing repeated vertices and orienting rings. GEOS cannot produce graded, multi-label, partial-credit signals.
- **Policy decision: 8-connected pinches.** Pixel-exact rings from the repo's polygonizer (8-connected foreground, which joins diagonal pixels at a saddle vertex) contain **vertex pinches**, which both GEOS and the repo's checker call self-touching. [local experiment] 3% (S=256) to 8% (S=1000) of synthetic masks gave rings that are not simple. Either SFT targets must remove pinches during simplification (§3), or the grader must explicitly allow vertex-touch at the same vertex, equivalent to GEOS `setSelfTouchingRingFormingHoleValid`.

### Gaps
- The GEOS IsValidOp algorithm (noding with monotone chains) and its complexity are not documented on the fetched page. JTS / Martin Davis blog posts on the 2021 IsValidOp rewrite were not found by search.
- CGAL `is_simple_2` / `Polygon_with_holes_2::is_valid` and Boost.Geometry's full `validity_failure_type` enum were not fetched.
- No Bentley–Ottmann implementation was benchmarked.

---

## 3. Mask → compact polygon: topology-preserving simplification with ≤ 1 px deviation

### Takeaway
Douglas–Peucker on pixel-corner staircase rings is a poor fit:
- at ε=0.5 it barely removes vertices (−19%), because staircase corners sit about 0.5–0.7 px off the ideal line;
- at ε=1.0 it removes 81% but drops the mean IoU to 0.976 and the minimum to 0.949, which **fails the 0.98 target**.

**Potrace's polygon stage** is minimum number of segments, where a segment is admissible if every skipped pixel corner lies within L∞ ½ of it (a digital-straightness criterion), with a deviation penalty as tie-break. It gave **fewer vertices than DP@1.0, IoU mean 0.982–0.992, 100% topology kept and 100% simple and non-crossing rings** in local tests.

Minimum-vertex simplification of a whole subdivision under topology constraints is NP-hard to approximate within n^(1/5−δ) (Estkowski–Mitchell). So the practical design is per-ring optimal min-# (Imai–Iri / potrace / DSS-style shortest path in a DAG of admissible shortcuts), followed by exact verification and local tightening.

### Cited Findings
- **Potrace (Selinger 2003), polygon stage:**
  - A path is "straight" if some segment passes within max-distance (L∞) ½ of every path vertex and not all four directions occur.
  - All straight subpaths of a closed path are found in **O(n²)** worst case.
  - Optimality is **lexicographic (number of segments, then penalty)**, where the penalty measures Euclidean deviation of the subpath from the segment.
  - Solving it is an optimal cycle in a DAG, O(n·m) with m the longest admissible segment.
  - Vertices are then "adjusted" off the integer corners. — [Selinger, Potrace: a polygon-based tracing algorithm (PDF)](https://potrace.sourceforge.net/potrace.pdf)
  - Licence: potrace itself is **GPL**. potracer is a pure-Python port. Licence status of the port was not verified this session.
- **Kolesnikov & Fränti:**
  - min-ε (fixed segment count M) and min-# formulations, solved by dynamic programming / A* / shortest path in O(N²)–O(N³).
  - "Reduced-search DP" gives a near-optimal solution much faster.
  - Multiresolution variants exist. — [Kolesnikov & Fränti paper (UEF)](https://cs.uef.fi/pub/franti/papers/Kolesnikov/Kolesnikov_Paper4.pdf), [UEF report A-2003-2](https://cs.uef.fi/pub/Reports/A-2003-2.pdf)
- **Estkowski & Mitchell, "Simplifying a polygonal subdivision while keeping it simple":**
  - Minimum-complexity simplification under an error bound, without crossings and without Steiner points, is **MIN PB-complete and NP-hard to approximate within n^(1/5−δ)**.
  - They propose a "Simple Detours" heuristic: standard chain simplification, then untangle. — [Stony Brook record](https://researchconnect.stonybrook.edu/en/publications/simplifying-a-polygonal-subdivision-while-keeping-it-simple/), [PDF](https://www.uni-trier.de/fileadmin/fb4/prof/INF/DEA/Seminar0708/Estkowski.pdf)
- **Guibas, Hershberger, Mitchell, Snoeyink:** minimum-link approximation of polygons and subdivisions, including hardness results for subdivisions. — [UBC TR-92-05](https://www.cs.ubc.ca/tr/1992/tr-92-05)
- **van Kreveld, Löffler, Wiratma:** optimal polyline simplification under Hausdorff and Fréchet distance. Found in search, not read. — [arXiv:1803.03550](https://arxiv.org/abs/1803.03550)
- **GEOS TopologyPreservingSimplifier:**
  - a DP-like "maximum distance difference" algorithm;
  - the result has the **same number of shells and holes, in the same order**;
  - rings touch at no more points than in the input;
  - the result is valid with the same dimension and component count;
  - all vertices lie within the tolerance of the original. — [GEOS TPS doxygen](https://libgeos.org/doxygen/classgeos_1_1simplify_1_1TopologyPreservingSimplifier.html)
- **DGtal (C++):** provides arithmetical DSS recognition (`ArithmeticalDSSComputer`) and `GreedySegmentation` of 4-connected Freeman chains into maximal DSSs. — [DGtal greedy DSS example](https://dgtal.org/doc/2.1/geometry_2curves_2greedy-dss-decomposition_8cpp-example.html). DGtal's licence (LGPL) is from domain knowledge and was not re-verified.
- **Clipper2** includes `SimplifyPaths`. — [Clipper2 GitHub](https://github.com/AngusJohnson/Clipper2)
- **Classic simplification references** (bibliographic, domain knowledge, not re-fetched):
  - Douglas & Peucker 1973 — [doi:10.3138/FM57-6770-U75U-7727](https://doi.org/10.3138/FM57-6770-U75U-7727)
  - Visvalingam & Whyatt 1993 — [doi:10.1179/000870493786962263](https://doi.org/10.1179/000870493786962263)
  - Saalfeld 1999, topologically consistent DP: a shortcut is safe if no foreign vertex lies in the region between shortcut and subpath, tested via convex-hull/winding parity — [doi:10.1559/152304099782424901](https://doi.org/10.1559/152304099782424901)
  - Debled-Rennesson & Reveillès 1995, linear-time arithmetical DSS recognition — [doi:10.1142/S0218001495000249](https://doi.org/10.1142/S0218001495000249)
  - Not linked and unverified: Imai & Iri 1986/1988 (min-# as shortest path in the shortcut graph); Chan & Chin 1996 (O(n²) min-# under the tolerance-zone criterion); Kovalevsky 1990 (DSS recognition).
- **[local experiment] Simplification study** on synthetic smooth blobs (Gaussian-filtered noise, thresholded, largest 8-connected component, holes allowed). Rings came from the repo's bit-quad `polygonize`. Every output was **rounded to the 0–1000 grid**, re-rasterized with the exact even-odd scanline of §1, and compared with the source mask (IoU, Betti numbers b0/b1, `shapely.is_simple` per ring, no ring–ring intersections).

  S=256, 40 masks:

  | Method | mean vertices | IoU min / mean | topology kept | simple & non-crossing | time |
  |---|---|---|---|---|---|
  | exact pixel-edge rings | 543.0 | 1 / 1 | 100% | **97%** (pinches) | 0.8 ms |
  | DP ε=0.5 (preserve_topology T/F) | 438.8 | 0.9991 / 0.9998 | 100% | 97% | 0.24 / 0.03 ms |
  | DP ε=1.0 (preserve_topology T/F) | 101.6 | **0.9486 / 0.9763** | 97% | 97% | 0.26 / 0.02 ms |
  | potrace polygon stage (alphamax=0, opticurve off, turnpolicy=black) | **94.5** | 0.9775 / **0.9918** | **100%** | **100%** | 6.8 ms (pure Python) |

  S=1000, 12 masks:

  | Method | mean vertices | IoU min / mean | topology kept | simple & non-crossing | time |
  |---|---|---|---|---|---|
  | exact rings | 2405.5 | 1 / 1 | 100% | 92% | 9.2 ms |
  | DP ε=0.5 | 1928.5 | 1.0 / 1.0 | 100% | 92% | — |
  | DP ε=1.0, preserve_topology=True | 428.8 | 0.9672 / 0.9769 | 100% | 92% | — |
  | DP ε=1.0, preserve_topology=False | 428.7 | 0.9672 / 0.9769 | **92%** | 92% | — |
  | potrace | **391.8** | 0.9719 / 0.9824 | 100% | 100% | 35 ms |

### Inferences
- **For staircase boundaries, the right admissibility test is digital straightness, not Euclidean DP distance to pixel corners.** Potrace's L∞-½ criterion and arithmetical DSS recognition (DGtal, Debled-Rennesson) both encode it. DP measured to corners penalises the inherent ±0.5 px staircase, so it either does nothing (ε=0.5) or overshoots (ε=1.0).
- **Recommended design for SFT targets** (all integer and deterministic):
  1. Exact rings from the bit-quad polygonizer.
  2. Per ring, build the admissible-shortcut graph. An edge (i→j) is admissible if (a) the subpath is digitally straight (potrace triple test, or DSS recognition, O(n²) worst case and near-linear for typical boundaries), and (b) the shortcut's endpoints are pixel corners, which keeps everything on integers.
  3. Find the min-# cycle with potrace's lexicographic (count, penalty).
  4. Map to the 0–1000 grid.
  5. Run the §2 validity pass plus the §1 re-raster, and check parts/holes and IoU.
  6. If a check fails (crossing with another ring, pinch, or IoU < threshold), forbid the offending shortcuts near the failure and re-solve that ring (Estkowski–Mitchell "detour"-style repair), or fall back to ε=0.5.

  Saalfeld's "no foreign vertex inside the shortcut triangle/region" test can be folded into admissibility (b) to prevent crossings up front.
- **Do not use potrace's vertex adjustment.** It moves vertices off pixel corners. After rounding to the 0–1000 grid that adds up to 0.5·S/1000 px of extra error, about 0.5 px at S=1000, and is a likely cause of the lower IoU at S=1000. Corner-anchored vertices plus crops with S dividing 1000 (125, 200, 250, 500, 1000) make every pixel corner an exact grid point. The round trip mask → rings → simplified → re-raster is then free of quantization error. Untested hypothesis, see H6.
- **The IoU ≥ 0.98 target is shape-dependent.** Boundary error is roughly 0.25–0.5 px × perimeter, so thin or small objects cannot reach 0.98 at 1 px deviation regardless of algorithm. A bound on symmetric difference relative to perimeter (e.g. ≤ 0.5 px × perimeter), or an adaptive ε per ring, is a more consistent acceptance rule. This is an inference from the measured IoU spread.
- **Pinch removal** is a by-product of shortcutting through the pinch vertex, and the min-# solution did so in 100% of cases here. A dedicated rule (split or bevel 8-connected pinches) is still needed for ε→0.
- **GEOS TPS** guarantees ring counts and validity but uses DP distance, so it inherits DP's staircase problem. Good as a fallback or oracle, not the primary tool.

### Gaps
- A true DSS-based (DGtal) segmentation was not benchmarked. Neither was Kolesnikov's optimal min-ε, nor Visvalingam–Whyatt.
- The experiment used synthetic blobs, not the project's real RIGA / LIDC masks. The repo has only CSV summaries locally, not the masks.
- The potracer port is pure Python (7–35 ms). The C potrace (GPL) speed for this workload was not measured.
- Whether potrace guarantees non-crossing between different paths is not stated in the parts of the paper read. Empirically there were 0 crossings in 52 masks.

---

## 4. Python ↔ C++ bridge in 2026, and parallelism for the reward pool

### Takeaway
Per-call binding overhead does not matter here: the kernels take 10–100 µs, while binding overheads are well under 1 µs. **What matters is that all per-vertex and per-pixel loops are compiled and release the GIL.**
- [local experiment] A numba `nogil` exact reward (3 rings × up to 120 vertices, eo+nz raster + IoU + overlap flag at S=512) for **256 rollouts takes 48 ms on 1 thread and 10.5 ms on 8 threads** with a plain `ThreadPoolExecutor` on standard (GIL) CPython. That is about 100× under the 1 s budget, so neither free-threading nor multiprocessing is needed.
- For production, a header-only int64 C++ core shared with the CUDA engine, bound with **nanobind + scikit-build-core + cibuildwheel**, is the cleanest. numba is a fine prototype and reference.

### Cited Findings
- **nanobind vs pybind11:** about 3× lower overhead for simple calls and about 10× lower when passing objects. 2.7–4.4× faster compiles. Binaries 3–5× smaller than pybind11 and 3–12× smaller than Cython. Runtime is comparable to Cython. Measured on Ryzen 9 7950X, Python 3.10, nanobind 1.2.0. Free-threading support is covered on its own docs page. — [nanobind benchmarks](https://nanobind.readthedocs.io/en/latest/benchmark.html), [nanobind "why"](https://nanobind.readthedocs.io/en/latest/why.html)
- **Free-threaded CPython:** experimental in 3.13 (PEP 703). In 3.14 it is "no longer considered experimental" (PEP 779) but is still not the default build. — [py-free-threading guide](https://py-free-threading.github.io/)
- **Free-threaded overhead:** single-thread overhead is "about 1% on macOS aarch64 to 8% on x86-64 Linux" on pyperformance. C extensions not marked free-thread-safe **re-enable the GIL at import**, with a warning. — [Python 3.14 free-threading HOWTO](https://docs.python.org/3.14/howto/free-threading-python.html)
- [local experiment] **Single-threaded per-polygon costs** at S=512 (median):

  | Operation | Cost |
  |---|---|
  | exact numba scan | 75–80 µs |
  | numba O(n²) validity | 0.7–117 µs |
  | shapely.is_valid | 2.6–20 µs |
  | rasterio | 190–310 µs |
  | skimage / matplotlib / shapely-points | 4–45 ms |
  | repo pure-Python validity | 0.25–43 ms |

- [local experiment] **Thread scaling:** 256-rollout batch takes 48.3 ms (1 thread), 25.1 ms (2), 13.4 ms (4) and 10.5 ms (8) on M3 (4P+4E). Scaling is near-linear to 4 threads, and E-cores limit the gains past that.

### Inferences
- **Where C++ (or compiled code) beats Python:**
  - all O(n²) segment tests;
  - all scanline and prefix-sum work;
  - polygonize and simplify (potrace O(n²) straightness);
  - the shortcut DAG.

  Pure Python is 300–400× slower on the validity test (43 ms vs 0.11 ms at n=200). Parsing JSON or text to int arrays is fine in Python: one `json.loads` plus `np.asarray` per rollout, which was not measured but is typically tens of µs.
- **Recommended stack:**
  - **C++17/20 header-only core** with an int64 frame. The same code compiles in nvcc for the CUDA engine and in the host extension.
  - **nanobind** bindings taking `nb::ndarray<int32/int64>`, with the GIL released inside (`nb::gil_scoped_release`).
  - **scikit-build-core** as the build backend, **cibuildwheel** for manylinux, macOS and Windows wheels.
  - Optionally declare free-threaded support with the module flag, so 3.14t users do not silently get the GIL back.
- **Cython** is equivalent in speed (per the nanobind benchmarks) but duplicates logic that the CUDA side needs in C++ anyway.
- **numba** suits the research phase. Its downsides are JIT warm-up, LLVM version pinning, and the fact that its code cannot be shared with CUDA C++.
- **mypyc** does not help numeric inner loops much. This is inference; no benchmark was found.
- **Pool choice:** threads plus GIL-releasing kernels beat a process pool, because there is no pickling and no per-process import. Multiprocessing overhead was not measured. Free-threaded 3.14 adds nothing if the kernels already release the GIL, and costs 1–8% single-thread overhead plus possible GIL re-enabling by unmarked extensions (numba / shapely free-threading status was not checked).

### Gaps
- No direct measurement of nanobind or pybind11 call overhead in ns on this machine.
- No multiprocessing-pool benchmark.
- The status of numba, shapely, OpenCV and rasterio on free-threaded 3.14 was not checked. See the [compatibility tracker](https://py-free-threading.github.io/tracking/).
- cibuildwheel and scikit-build-core docs were not fetched this session. Recommending them is standard practice, not newly sourced evidence.

---

## 5. Existing libraries we could just use (licence, speed, exactness)

### Takeaway
None of these is the core grader or rasterizer:

| Library | Licence | Role |
|---|---|---|
| GEOS / shapely 2.x | LGPL / BSD | fast first-error validity, test oracle, TPS fallback |
| rasterio / GDAL | MIT / BSD | fast exact even-odd centre raster; on-edge ties unspecified; no nonzero |
| Clipper2 | BSL-1.0 | int64 boolean ops with NonZero / EvenOdd; useful to build union-vs-raw self-overlap oracles; no binding listed |
| skimage / matplotlib | BSD | exact even-odd but 100–500× too slow |
| PIL / OpenCV / pycocotools | various | not centre-exact |
| potrace | GPL | best simplification quality found, but GPL; re-implement the polygon stage (about 300 lines) |
| DGtal | LGPL (unverified) | DSS recognition reference |

### Cited Findings
- Clipper2: BSL-1.0, int64 and double, NonZero and EvenOdd, `SimplifyPaths`. Ports to Rust, Go, Java and others. **No Python binding listed** in the README. — [Clipper2 GitHub](https://github.com/AngusJohnson/Clipper2)
- shapely `make_valid` (linework = even-odd reconstruction) is usable as an even-odd area oracle for self-intersecting rollouts. — [shapely make_valid](https://shapely.readthedocs.io/en/stable/reference/shapely.make_valid.html)
- GEOS TPS guarantees on rings and validity. — [GEOS TPS](https://libgeos.org/doxygen/classgeos_1_1simplify_1_1TopologyPreservingSimplifier.html)
- Boost.Geometry `is_valid` with message is worst-case quadratic for polygons. — [Boost.Geometry is_valid](https://www.boost.org/doc/libs/release/libs/geometry/doc/html/geometry/reference/algorithms/is_valid/is_valid_2_with_message.html)
- pycocotools polygon→RLE is approximate (×5 upsampling, rounding). — [maskApi.c](https://github.com/cocodataset/cocoapi/blob/master/common/maskApi.c)
- GDAL centre-sampled rasterization, with on-centre behaviour unspecified. — [gdal_rasterize](https://gdal.org/programs/gdal_rasterize.html)
- tiny-skia FillRule Winding / EvenOdd. — [tiny-skia docs](https://docs.rs/tiny-skia/)

### Inferences
- **Licence-safe production path:** our own C++ code (validity, scanline, simplification), with GEOS, Clipper2 and rasterio used only in tests as oracles.
- **Avoid shipping potrace** (GPL) inside a wheel unless the project is GPL-compatible.

### Gaps
- Not evaluated: CGAL (GPL/LGPL, heavy), manifold, earcut, polylabel, geo (Rust) via Python, `pyclipr` / Clipper2 Python wrappers. They were considered out of the critical path because none provides exact centre-sampled nonzero rasterization or multi-label validity grading. This is an inference and was not checked.

---

## Candidate hypotheses to test locally

The H1 and H2 harnesses already exist in `prototypes/poly_fillcheck.py` and `poly_fillcheck2.py`. H4 is in `poly_timing.py`, H5 in `poly_simp.py`, and H7 in `poly_pool.py`.

- **H1, rasterizer exactness.**
  - Test: the integer top-left-rule scanline equals the O(S²·n) brute-force winding reference bit for bit, for both rules, on ≥10,000 polygons. Cover star, random self-intersecting, rectilinear, collinear-heavy and duplicate-vertex polygons, plus multi-ring inputs with holes, S ∈ {100, 125, 250, 384, 500, 1000}, and the CUDA port.
  - Pass: 0 mismatching pixels. Measured so far: 0 on 500 polygons at S=250, single ring.
- **H2, library tie-breaks.**
  - Test: for skimage, matplotlib, shapely and rasterio, off-boundary mismatches against exact even-odd are 0 (so far 0 on 1,200 polygons). Tabulate their on-edge behaviour per edge orientation to document each convention.
  - Expected: PIL and OpenCV are rejected (>100 px per polygon off-boundary). Re-test PIL across Pillow 9.5 / 10.x / 11.x / 13.0 to see whether semantics drift.
- **H3, self-overlap flags.**
  - Test: on figure-8s, bow-ties, doubly-wound loops and spiral overlaps, check that `eo XOR nz` misses opposite-orientation lobes and that `any(w>0) & any(w<0)` catches them. Measure recall of the combined raster flag against the exact segment classifier's "self-crossing" label on 10k random rings.
- **H4, validity cost and agreement.**
  - Test: the exact int64 one-pass classifier (O(n²) with bounding-box rejection, then a uniform-grid variant) agrees with GEOS `is_valid` on every input where GEOS semantics coincide (after dropping repeated vertices, ignoring winding, and accepting single-point hole touches). Inputs: 100k random multi-ring rollouts.
  - Report: per-class confusion. Time: under 150 µs at 10 rings × 200 vertices (2,000 edges, O(n²) about 2M pairs). If not, switch to the grid.
- **H5, simplification vertex count.**
  - Test: on real RIGA / LIDC / target-domain masks (not synthetic), compare mean vertices at matched constraints (re-raster IoU ≥ 0.98, b0/b1 preserved, rings simple and non-crossing after grid rounding) for:
    - (a) DP ε ∈ {0.5, 0.75, 1.0};
    - (b) GEOS TPS;
    - (c) Visvalingam–Whyatt;
    - (d) potrace polygon stage with and without vertex adjustment;
    - (e) DGtal greedy maximal-DSS segmentation;
    - (f) our corner-anchored min-# with Saalfeld admissibility.
  - Prediction: (f) ≤ (d) < (a@1.0) in vertices, with 100% topology. Synthetic results so far: potrace 94.5 vs DP@1.0 101.6 at S=256.
- **H6, grid-divisor crops.**
  - Test: with S ∈ divisors of 1000 and corner-anchored vertices, the round trip mask → exact rings → identity (no simplification) → grid → exact re-raster reproduces the mask exactly (IoU = 1, same hash), for 10k masks. With S=384 it should not (quantization).
  - Decides whether crops should be snapped to divisor sizes.
- **H7, throughput.**
  - Test: the full reward (parse + validity + eo/nz raster of all rings + IoU vs GT + flags) for 256 rollouts × up to 10 rings × 200 vertices at S=1000 stays under 100 ms on 8 cores (threads, GIL released).
  - Variants: numba, a nanobind C++ build, and a process pool. Measured so far: 10.5 ms for 3 rings × ≤120 vertices at S=512.
- **H8, binding overhead.**
  - Test: nanobind vs pybind11 vs Cython per-call overhead with an `int64[n,2]` ndarray argument and GIL release, 1M calls.
  - Prediction: under 1 µs each, so irrelevant next to 10–100 µs kernels.

## Recommendation

1. **Rasterize with our own exact integer scanline.** Use one pass per rollout over all rings, a signed crossing difference array and a row prefix sum to get the winding w. Derive even-odd = w odd and nonzero = w ≠ 0 from the same buffer, plus the flags eo≠nz (|w| ≥ 2) and mixed-sign w (inverted lobes).
   - Spec: pixel centres at (i+½, j+½) and the **top-left tie-break rule**, written into the docs.
   - Drop PIL entirely: it is not centre-exact, and `outline=1` adds nothing.
   - Keep rasterio and skimage only as test oracles for even-odd.
2. **Validity: one exact int64 pass over all edges of all rings** with per-pair classification into graded flags:
   - self-crossing, self-touch, collinear overlap, spike;
   - ring-crossing, ring-touch;
   - hole outside its outer ring, nesting errors;
   - wrong winding by shoelace sign, zero area, repeated vertices.

   Start with O(n²) plus bounding-box rejection. Move to a uniform grid only if n grows. Never call `make_valid`. Use GEOS `is_valid` as a regression oracle only. Replace the pure-Python `ring_problems` now: at 41 ms for n=200 it alone can blow the step budget.
3. **SFT simplification: re-implement potrace's polygon stage** (straight-path test plus min-# with penalty tie-break), anchored to pixel corners, with no vertex adjustment, which avoids the GPL code. Add Saalfeld-style foreign-vertex admissibility and an exact re-raster verification loop that tightens locally on failure.
   - Treat 8-connected pinches explicitly.
   - Prefer crop sizes that divide 1000 so the grid round trip is exact.
   - Revisit the IoU ≥ 0.98 acceptance rule for thin or small objects; a perimeter-normalised symmetric difference bound is a candidate.
4. **Implementation: header-only C++ int64 core** shared with CUDA, exposed through nanobind (GIL released) and built with scikit-build-core and cibuildwheel. Run a `ThreadPoolExecutor` reward pool on standard CPython. Free-threaded 3.14 and multiprocessing are not needed at the measured about 10 ms per 256 rollouts. Keep the numba versions in `prototypes/` as the executable reference spec and the differential-test oracle for the C++/CUDA ports.
