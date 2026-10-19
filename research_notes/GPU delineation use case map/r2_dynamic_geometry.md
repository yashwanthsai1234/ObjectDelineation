# Dynamic / Streaming Geometry: incremental polygonization, temporally stable polygon IDs, geometry diffs (2025–2026)

Scope: whether a deterministic, seam-free GPU engine (binary bitmap with halos → 2x2 bit-quad classification → closed integer rings with holes) enables (a) re-polygonizing only the tiles that changed, (b) polygon IDs that stay stable over time, and (c) cheap, exact geometry diffs ("git for geometry from AI").

Method note: the session's shared web-search budget ran out partway through, so some threads stop at search-snippet level (marked below). Primary measurements were made directly. DuckDB was run over Overture's public S3 changelogs (2026-07, 08, 09) and the releases (2026-08-19.0 vs 2026-09-23.0) to count change types and to test whether "geometry changed" buildings really changed shape. Round-1 facts reused here keep their original citations.

---

## 1. Algorithms: dynamic/incremental CCL, contour maintenance, streaming polygonization, temporal coherence of contours in video

### Takeaway
Single-pass streaming CCL and raster-to-vector conversion are mature on FPGA and in patents: 4K60 hardware, two-scanline memory. GPU CCL research targets static images. "Dynamic" CC on GPU exists only for graph streams. I found no published GPU method that maintains contours or polygons incrementally as a mask changes. In video, temporal coherence is handled at the mask/model level (SAM 3 memory) or by point-correspondence polygon trackers. Annotation tools still need polygons whose vertex order is consistent from frame to frame.

### Cited Findings
- **Streaming CCL in hardware is mature.** A single-pass CCL FPGA implementation handles 4K/UHD (3840x2160) at 60 fps at 2 or 4 pixels per clock. An older single-pass design handles 640x480 at 30 fps. — [Springer J. Signal Proc. Systems s11265-021-01636-4](https://link.springer.com/article/10.1007/s11265-021-01636-4); [SPIE 11914 streaming CCL on FPGA](https://nanophotonics.spiedigitallibrary.org/conference-proceedings-of-spie/11914/119141I/Connected-component-labeling-algorithm-in-streaming-image-processing-with-FPGAs/10.1117/12.2613546.full) (search-snippet level)
- **LinkRunCCA** is an open Verilog core for single-pass connected-component analysis. It runs at one pixel per cycle with very low memory, and it emits bounding boxes, not contours. — [OpenCores LinkRunCCA](https://opencores.org/projects/linkruncca)
- **Memory-bounded raster-to-vector.** One patented method keeps only two scanlines in memory and scans the raster once. It notes that the usual approach (tile, vectorize, then stitch at the end) "scales up badly as the image grows". — [US8731309B1](https://patents.google.com/patent/US8731309)
- A scan-line method that builds region boundaries dynamically from "Pixel Directed Edges" in linked lists converts raster regions to vectors in one pass, with no intermediate data — [JCAD (jcad.cn)](https://jcad.cn/en/article/id/32388d86-a7f6-4b21-961b-5bca61ffa5a9) (search-snippet level)
- **SweepLSD** (Y. Shimizu, arXiv 2608.22086, Aug 22 2026) is a one-pass line-segment detector with O(width) memory, an integer-only streaming core and a CCL stage. It runs about 11 ms per Full-HD frame on one CPU thread, 4.6x/5.2x/25x faster than ELSED/EDLines/LSD. Its FPGA version is **bit-exact** with the software and runs 1080p30 without a frame buffer. — [arXiv 2608.22086](https://arxiv.org/abs/2608.22086)
- **GPU CCL is a static-image problem in the literature.** The work covers iterative label equivalence, block-based union-find and direct algorithms (DLP-GPU). — [DASIP18 CCL GPU (LIP6)](https://largo.lip6.fr/~lacas/Publications/DASIP18_CCL_GPU.pdf); [GTC17 DLP](https://largo.lip6.fr/~lacas/Publications/GTC17_DLP.pdf). "Incremental Connected Component Detection for Graph Streams on GPU" (Electronics 12(6):1465, 2023) handles streaming graph edges, not image masks. The full text returned 403, so whether it supports deletions was not verified — [MDPI](https://www.mdpi.com/2079-9292/12/6/1465)
- **No GPU findContours.** OpenCV CUDA has no findContours, so users copy GPU masks back to the CPU to extract contours (round 1) — [OpenCV Q&A](https://answers.opencv.org/question/181859/why-no-cuda-version-of-findcontours/)
- **Polygon tracking in video.** *Polygonal Point Set Tracking* (CVPR 2021) carries a set of contour points from frame to frame, giving both the contour and point-to-point correspondence for VFX (motion tracking, deformation, texture mapping) — [arXiv 2105.14584](https://ar5iv.arxiv.org/html/2105.14584). *Coherent Parametric Contours* models object boundaries as parametric curves with spatio-temporal optimization for temporally stable boundaries — [PapersWithCode](https://astro.paperswithcode.com/paper/coherent-parametric-contours-for-interactive). A 2025 vectorized video representation uses "spatio-temporally consistent proxy embedding" — [arXiv 2510.12256](https://arxiv.org/pdf/2510.12256)
- **CVAT polygon tracks.** A polygon track is stored as keyframes, and each vertex is interpolated independently between them. The annotator must keep the **same starting point and direction** on each new keyframe polygon. **Masks are not interpolated** and must be set on every frame. — [CVAT docs: track mode with polygons](https://docs.cvat.ai/docs/annotation/manual-annotation/shapes/annotation-with-polygons/track-mode-with-polygons/); [CVAT manual annotation](https://www.mintlify.com/cvat-ai/cvat/annotation/manual-annotation)
- **SAM 3 video tracking.** Identity comes from per-video model memory. The Roboflow SAM3 Video Tracker block outputs masks with tracker IDs, and its docs mention no polygon conversion. It needs a persistent WebRTC session because on stateless or multi-replica runtimes "the state resets between calls and the output is meaningless for tracking". It requires CUDA. — [Roboflow SAM3 Video Tracker block](https://docs.roboflow.com/workflows/blocks/blocks/video-processing/sam3-video-tracker); [Roboflow changelog, June 2026](https://docs.roboflow.com/changelog/explore-by-month/june-2026/sam-3-point-prompting-and-video-tracking-in-workflows). From round 1: SAM 3 outputs masks, boxes and track IDs but no polygons, and SAM 3.1 tracks 16 objects at 32 fps with about a 7x speedup at 128 objects — [Ultralytics SAM 3 docs](https://docs.ultralytics.com/models/sam-3)
- **Online HD maps.** RTMap (arXiv 2507.00980, Jul 2025) sorts prior-map vector elements into matched, outdated and newly observed. It associates elements across traversals with point-wise matching and **union-find over vertices**, then fuses them with MAP estimation. It reports no latency or bandwidth numbers. — [arXiv 2507.00980](https://arxiv.org/html/2507.00980v1). ArgoTweak (ICCV 2025) breaks map changes into **atomic, element-level changes** through a bijective mapping, so unchanged elements are preserved — [ICCV 2025 open access](https://openaccess.thecvf.com/content/ICCV2025/html/Wild_ArgoTweak_Towards_Self-Updating_HD_Maps_through_Structured_Priors_ICCV_2025_paper.html) (search-snippet level; the page returned 403)
- **Incremental navmesh.** The SE(2) Navigation Mesh (arXiv 2607.01454, Jul 2026) builds polygonal navmeshes from per-yaw footprint masks and updates them online from streaming point clouds (round 1) — [arXiv 2607.01454](https://arxiv.org/abs/2607.01454)

### Inferences
- With halos, re-polygonizing only the dirty tiles is computationally trivial: the work is local, and the engine is already deterministic and seam-free. The hard part is everything after that:
  1. re-stitching rings that cross a dirty/clean tile boundary;
  2. keeping ring IDs and canonical vertex order stable when a ring is regenerated;
  3. emitting a diff record.
  None of the found literature solves this for raster-derived polygons on GPU, so this is a real but narrow gap.
- CVAT's "same start point, same direction" requirement is a manual workaround for missing **canonical ring ordering**. An engine that fixes orientation (for example, CCW outer and CW holes) and a deterministic start vertex (for example, the lexicographically smallest integer vertex) produces polygon tracks that interpolate and diff cleanly.
- Hardware CCL work (bit-exact FPGA, two-scanline memory) shows that people who stream line-scan data value determinism and O(width) memory. A GPU engine that accepts row bands with a halo and emits closed rings as each band completes would match line-scan sensors (NIROPS IR scanners, pushbroom satellites, web inspection). This is speculation: I found no buyer asking for it.

### Gaps
- I found no paper on "incremental marching squares" or dynamic contour/polygon maintenance under pixel flips, on CPU or GPU. The search budget ran out before a second pass with alternative terms (e.g., "dynamic isocontour", "contour update local edit").
- Not reached: Deep Snake video variants, rotoscoping spline-tracking products (e.g., Mocha, Nuke), and DynLP (arXiv 2604.06596, GPU dynamic batch label propagation, title only).

---

## 2. Object identity over time: tracking by overlap, merge/split event graphs, topology-change events. Raster or polygon? GPU?

### Takeaway
Scientific trackers assign identity by **raster overlap or point distance**, not polygon geometry:
- TempestExtremes StitchBlobs and py-eddy-tracker's AreaTracker use overlap.
- tobac merge/split uses a minimum spanning tree over feature points.
- The topology community uses merge trees with optimal transport.

Identity for map features (Overture) is assigned by **polygon IoU** (≥0.5). None of the scientific trackers found has a GPU path, except NVIDIA Earth2Studio's torch TC tracker, which is distance-based. Merge and split are tracked as events in science, but Overture's changelog has no merge/split type.

### Cited Findings
- **TempestExtremes StitchBlobs** tracks areal features with an overlap criterion. A blob takes the ID of the earlier blob it overlaps most. On a split, the child with the largest overlap keeps the ID and the others get new IDs. Merge events are also distinguished. — [Ullrich TempestExtremes poster 2019](https://climate.ucdavis.edu/Ullrich_TempestPoster2019_v2.pdf); [UC Davis TempestExtremes](https://climate.ucdavis.edu/tempestextremes.php). From round 1: DetectBlobs over 41 years of ERA5 took 34 min 42 s, 66% of it file input, with no GPU path — [GMD preprint](https://gmd.copernicus.org/preprints/gmd-2020-303/gmd-2020-303-manuscript-version4.pdf)
- **tobac v1.5** (GMD 2024) added 3-D tracking, splitting and merging, and periodic boundaries. The `merge_split` post-processing uses a **minimum Euclidean spanning tree** on feature positions, with a default distance threshold of **25x the x/y grid spacing**. That is point-based, not shape-overlap. — [tobac.merge_split docs](https://tobac.readthedocs.io/en/latest/api/generated/tobac.merge_split.html); [tobac GMD 2024](https://gmd.copernicus.org/articles/17/6035/2024/gmd-17-6035-2024.html). No GPU support was mentioned.
- **py-eddy-tracker** finds eddies as closed SSH contours (2 mm steps, about 5,936 eddies per day, about 59.4M since 1993) and has an **AreaTracker** that links eddies by contour overlap. It is written in Python (round 1). — [ESSD 2022](https://essd.copernicus.org/articles/14/1087/2022/essd-14-1087-2022.html); [PET AreaTracker](https://py-eddy-tracker.readthedocs.io/en/stable/_autosummary/py_eddy_tracker.featured_tracking.area_tracker.AreaTracker.html)
- **ProbSevere** matches storm objects between consecutive MRMS images about every 2 minutes and outputs GeoJSON polygons. Testbed latency was 2–5 minutes (round 1). — [CIMSS](https://cimss.ssec.wisc.edu/severe_conv/training/slides/radar_data.html); [NSSL EWP blog](https://inside.nssl.noaa.gov/ewp/?p=19142). **Threats-in-Motion** warning polygons are meant to update every minute and move with the storm — [NOAA GSL TIM](https://gsl.noaa.gov/tim/)
- **Topology-based tracking (TVCG 2025).** "Flexible and Probabilistic Topology Tracking with Partial Optimal Transport" (Li, Yan, Yan, Needham, Wang) compares merge trees across time steps using partial optimal transport. It produces **probabilistic tracking graphs** instead of explicit overlap matching. No GPU or performance figures were given. — [IEEE VIS 2025](https://aa.ieeevis.org/year/2025/program/paper_0ffcd3a5-a31e-4532-9adb-8306bbb1518f.html); [arXiv 2302.02895](https://arxiv.org/abs/2302.02895). Earlier merge-tree region tracking builds a DAG from local overlap and similarity and solves shortest paths — [Saikia & Weinkauf](https://people.kth.se/~weinkauf/publications/abssaikia20.html)
- **NVIDIA Earth2Studio** has a GPU tracker, TCTrackerWuDuan: threshold vorticity, run CCL, erode, then link paths within a 300 km radius. Linking is distance-based, not shape-based (round 1). — [Earth2Studio docs](https://nvidia.github.io/earth2studio/main/modules/generated/models/px/WeatherNext2Cyclones/)
- **Polygon identity in maps.** Overture's building matching workflow uses **IoU ≥ 0.5**, with a fallback of IoU ≥ 0.3 plus the Overture centroid falling inside the other polygon. About 13% of matches in the demo used the fallback. Buildings matched 350/412 (85%) in one 1°x1° cell. The page lists fragmentation (building vs. building parts, one land polygon overlapping several) as a pain point. — [Overture Labs workshop: buildings matching](https://labs.overturemaps.org/workshop/7-buildings-matching.html)
- GERS IDs are "intended to be stable (within a reasonable tolerance of error)" across versions — [Overture GERS guide](https://docs.overturemaps.org/guides/gers/)

### Inferences
- Across these domains, identity comes from a **label-overlap matrix** between time t and t+1: counts of pixel pairs (label_t, label_t+1). On GPU that is one histogram pass over two label maps. Polygons are not needed to compute it, and polygon-polygon IoU is the slower way to get the same answer. The engine's edge for identity is therefore less about polygons than about producing **label maps plus per-component stats** (area, bbox, Euler number from the bit-quad counts it already computes) in the same pass, then emitting a birth/death/continue/merge/split event graph.
- Bit-quad classification yields the Euler number per component for free (Gray's bit-quad formula). A change in a component's Euler number or hole count between versions is a cheap, exact topology-change event (hole opened or closed). None of the trackers found expose this. This is an inference, not documented demand.
- A GPU tracker for StitchBlobs-, tobac- or PET-style workloads is unclaimed space. The binding constraint measured so far is I/O (66%), so the value is in-line fusion with ML inference, not faster tracking alone.

### Gaps
- I found no GPU implementation of overlap-based merge/split tracking for gridded fields. Euler-number or Betti-number change detection has been formalized in TDA, but no operational product uses it.
- TRACE and newer eddy and cyclone trackers were not reviewed in round 2 (search budget).

---

## 3. Geometry versioning and diffs: GERS, OSM, Kart/GeoGig, Iceberg v3, change-polygon products

### Takeaway
Production "diffs" today come in three forms:
1. **ID-joined record comparison**: Overture's changelog holds only id, bbox and columns_changed, no geometry diff.
2. **Raster differencing then vectorization**: Planet change feeds.
3. **Row-level version control**: Kart stores one file per feature, and Iceberg offers snapshot versioning.

None of them computes exact geometric change polygons between versions of AI-derived layers. The documented pain is **ID instability** (Overture had a full ID break in June 2025 and a 100% GERS churn event in the base theme in Sept 2026). Diff compute cost is not the documented pain. Measured monthly churn for buildings is 0.6–2.4% of features in normal months. In September 2026, 27.5% (696.7M) of buildings were flagged geometry-changed, and in sampled windows about 97% of those were identical shapes with a rotated or reordered ring. Overture's diff does not canonicalize rings, so false diffs dominate.

### Cited Findings
- **GERS.** Since the June 2025 release, all GERS IDs are UUID v4 strings that Overture keeps stable across releases. June 2025 was a one-time breaking change, shipped with May→June ID mapping Parquet files, and no changelog was produced that month. Places example: of 5,727 new places in Salt Lake City in May, only 1,546 (27%) persisted into June, which Overture calls "a known quality issue with our source datasets". — [Overture blog: Using the GERS system](https://docs.overturemaps.org/blog/2025/06/25/getting-started-gers/); [2025-06-25 release notes](https://docs.overturemaps.org/blog/2025/06/25/release-notes/)
- **Sept 2026 ID churn.** The 2026-09-23 release had "100% GERS churn" in the base theme because of "a change in how UUIDs are generated". Patch 2026-09-23.1 fixed "a bug in the new GERS ID generation for base and building_parts features". In the same release: 2,533,842,612 buildings (+0.17%), Microsoft ML Buildings −3.41%, Google Open Buildings +2.91%, Esri −6.65%. — [Overture 2026-09-23 release notes](https://docs.overturemaps.org/blog/2026/09/23/release-notes/)
- **Changelog design.** Each release has a changelog of Parquet files partitioned by theme/type/change_type. Change types are `added`, `removed`, `data_changed` and `unchanged`. Columns: `id`, `bbox`, `columns_changed`. Geometry is compared "on its normalized text form rather than its raw bytes". `sources` and `confidence` are excluded because they "change on nearly every feature in most releases". There is no change type for splits, merges or restored IDs; these appear as `added`. — [Overture GERS changelog docs](https://docs.overturemaps.org/gers/changelog/)
- **Measured (this research, DuckDB over the public S3 buildings changelogs, `type=building`):**

  | release | added | removed | data_changed | unchanged | release total |
  |---|---|---|---|---|---|
  | 2026-07-22.0 | 6,118,173 | 4,746,103 | 5,793,051 | 2,533,251,510 | 2,545,162,734 |
  | 2026-08-19.0 | 18,695,740 | 34,275,861 | 7,178,501 | 2,503,708,372 | 2,529,582,613 |
  | 2026-09-23.0 | 39,201,698 | 34,941,699 | **696,708,335** | 1,797,932,579 | 2,533,842,612 |

  Each release total equals added + changed + unchanged and matches the published counts (the Sept total matches the release notes exactly). File sizes for 2026-08-19.0: added 0.77 GB, removed 1.39 GB, data_changed 0.30 GB, unchanged 82.2 GB. The 2026-09-23.1 buildings prefix holds only an `unchanged` partition. — source: `s3://overturemaps-us-west-2/changelog/<release>/theme=buildings/type=building/` (public bucket, documented in [Overture changelog docs](https://docs.overturemaps.org/gers/changelog/))
- **What "data_changed" means (measured).** I sampled the `columns_changed` arrays in 8 data_changed files per release. **Geometry** was the changed column in 7,860,070 of 7,862,875 sampled Sept rows (99.96%) and 286,113 of 287,241 Aug rows (99.6%). Attribute edits (class, subtype, names, num_floors) are a rounding error.
- **The Sept 2026 "geometry changes" are almost entirely the same shapes re-encoded (measured).** I joined data_changed buildings to their 2026-08-19.0 and 2026-09-23.0 geometries in four windows (0.03–0.05°), using DuckDB spatial:

  | window | sample | `ST_Equals` (same shape) |
  |---|---|---|
  | France (2.72E, 44.90N) | 177 | 78.5%, with median IoU 1.0 and median centroid shift about 1e-10° |
  | Illinois (−87.88, 42.30) | 114 (103 Microsoft ML + 11 OSM) | 100%, also identical after `ST_Normalize` |
  | Philippines (124.41, 6.41) | 1,173 Microsoft ML | 100% |
  | Philippines, same window | 2 Google | 0% (real change) |

  Pooled: **about 1,426 of 1,466 (≈97%) "changed" geometries are the same polygon with a different ring start vertex or order.** Example, same building in both releases: Aug `POLYGON((2.6754781 44.8812556, 2.675493 44.881284, …))` vs Sept `POLYGON((2.6754192 44.8812712, 2.6754781 44.8812556, …))`. Same vertices, rotated start. The source in the new release is Microsoft ML Buildings `version 2026-08-11`.

  Caveat: four non-random windows. The measurements are mine, not an Overture-published analysis.
- **Mapbox MTS incremental updates.** These require "stable, unique feature IDs" that the user assigns and maintains. They are recommended only when **under 15%** of source data changes per update. They are explicitly **not recommended for "datasets without stable unique identifiers (weather data)"**. Changesets are line-delimited GeoJSON with adds, modifies and deletes, up to 10 files and 50 GB. — [Mapbox MTS incremental updates](https://docs.mapbox.com/mapbox-tiling-service/guides/incremental-updates/)
- **Kart.** Kart is built on Git. "Every database table row is stored in its own file" (MessagePack), geometry is GeoPackageBinary/WKB, and paths are sharded by base64 integer PK or a SHA-256 hash of the PK. Diffs are therefore row/feature-level. Kart supports Postgres, SQL Server, MySQL and GeoPackage working copies, plus raster and point-cloud datasets through Git LFS. — [Kart table v3 spec](https://docs.kartproject.org/en/latest/pages/development/table_v3.html); [Kart docs](https://docs.kartproject.org/en/latest/); [Koordinates Kart](https://koordinates.com/products/kart/)
- **Iceberg v3 / Parquet** add native `geometry` and `geography` types, CRS in the schema, and bounding-box column statistics for file skipping (Parquet 2.11 GEOMETRY logical type). Adoption: AWS Glue at re:Invent 2025, Snowflake as the first major engine to ship v3 geo types (2026), Dremio GA, and Sedona reading and writing Iceberg geo through its pre-v3 Havasu encoding. The Parquet spec PR drew 400+ comments and the Iceberg geo spec 240+. The end-to-end Iceberg Parquet-geo read/write path was still under review at the time of the post. — [Wherobots: Iceberg geospatial types v3](https://wherobots.com/blog/iceberg-geospatial-types-v3/); [Snowflake Iceberg v3 blog](https://snowflake.com/en/blog/apache-iceberg-v3-table-spec-oss-shared-success). I found no Iceberg feature for geometric (shape-level) diffs between snapshots.
- **Planet Road & Building Change Detection** (Analytic Feeds):
  - Weekly or monthly cadence, output as vector GeoJSON.
  - Building F1 0.81, roads F1 0.72, change-detection precision 0.94.
  - Change polygons are mapped to a "standard 8x8 pixel gridcell" with a z-score threshold of 10.5.
  - Pipeline: a modified U-Net per image, then weekly/monthly aggregation, then a predictive time-series model that flags change.
  
  — [Planet techspec](https://docs.planet.com/data/analytic-feeds/road-building-change-detection/techspec/). Planet engineers describe a temporal inference model that denoises monthly rasters into a "consistent state feed". "Successive rasters" are then differenced into a change raster and run through "vector simplification to serve change vectors". — [ODSC / Planet](https://opendatascience.com/working-towards-planetary-scale-location-insights/)
- **OSM.** About 4 million map changes per day as of 2021-Q4 ([OSM wiki Stats](https://wiki.openstreetmap.org/wiki/Stats), search-snippet level). In Japan, changesets grew from 144,625 (2022) to 236,203 (2025). The share of AI-assisted (RapiD/MapWithAI) changesets flagged as suspect fell from 64% to 16%, but there are only about 1,700–1,800 per year — [ISPRS Archives L-4-W1-2026-275](https://isprs-archives.copernicus.org/articles/L-4-W1-2026/275/2026/isprs-archives-L-4-W1-2026-275-2026.pdf). Over 50% of AI-generated buildings in OSM are in the US and 75% are in five countries — [Masaryk Univ. publication](https://ics.muni.cz/en/research/publications/2482277)
- **Google Open Buildings 2.5D Temporal** is distributed as annual 2016–2023 **rasters** of presence, count and height over about 58M km², not as polygons (round 1). — [Earth Engine catalog](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_Research_open-buildings-temporal_v1?hl=en)
- **Incremental tiling exists downstream.** Besides Mapbox MTS, SuperMap has tile-update tooling — [SuperMap iDesktopX Tile Update](https://help.supermap.com/iDesktopX/en/Tutorial/MapTiles/UpdateTiles.html)

### Inferences
- **Calculation: ordinary months are small.**
  - July: about 0.6% of features touched (16.7M).
  - August: about 2.4% (60.1M of about 2.545B).
  
  A geometry-carrying delta in such months would be roughly 40–150x smaller than a full re-publish. These are calculations, not reported figures.
- **Calculation: September was not ordinary.** 27.5% of all buildings (696.7M) were flagged `data_changed`. The samples suggest about 97% of those are **spurious diffs**: identical polygons whose ring start vertex or order changed when Microsoft ML Buildings refreshed (version 2026-08-11). Extrapolating, the order of **~650–680M false "geometry changed" records in one month**. That is a large-scale, documented case of non-canonical ring encoding poisoning a production geometry diff. Overture's "normalized text form" comparison evidently does not canonicalize the ring start vertex. A polygonizer that emits **canonical rings** (fixed orientation, deterministic start vertex, integer grid) removes this whole class of false change. Equivalently, any diff tool can fix it by comparing on `ST_Normalize`d geometry; this is cheap, so the value is in the convention, not in GPU compute.
- **Inference: identity turnover from source swaps.** In August, `removed` (34.3M) was nearly 5x `data_changed` (7.2M). In September, 39.2M were added and 34.9M removed, while source counts swung about ±3% per month (MS −3.41%, Google +2.91%). This suggests heavy **identity turnover from source swaps or conflation**. Content-addressed canonical geometry would help separate "same shape, new ID" from real change. Overture does not publish such a breakdown.
- Planet's pipeline (difference rasters, then vectorize change) is the natural pattern for a deterministic engine. Given two aligned bitmaps, polygonize(A AND NOT B) and polygonize(B AND NOT A) give **exact signed change polygons** with no floating-point overlay, no slivers and no tolerance tuning. Doing the same with vector overlay (ST_Difference on independently vectorized polygons) creates slivers because of vertex jitter. This advantage depends on both versions sharing the same pixel grid, which holds for repeat-pass rasters on a fixed tiling but not for vectors coming from different sources.
- Mapbox's rule ("no stable IDs, as with weather data → don't use incremental updates") is direct evidence that raster-derived layers currently **cannot use incremental vector infrastructure**, because they lack stable IDs. That is the clearest documented demand signal for (b) plus (c).
- Kart (row-level) and Iceberg (snapshot-level) version **records**, not **shapes**. A "git for geometry" layer would sit above them: canonical rings, content hashes, an ID-linking graph and signed change polygons, all written into GeoParquet/Iceberg. I found no product doing this for AI-derived geometry.

### Gaps
- I found no published Overture analysis of the Sept 2026 mass `data_changed` event or of why building removals outnumber changes. The ~97% spurious figure comes from 4 non-random windows (1,466 buildings); a random sample across all 512 files would firm it up.
- GeoGig's current status (believed dormant) and Dolt or other "git for data" geo support were not checked because the search budget ran out.
- Esri change-detection outputs and Google Dynamic World change products were not reviewed in round 2.
- I found no source on the compute cost of producing Overture's changelog or conflation.

---

## 4. Bandwidth: edge, satellite and vehicle devices sending geometry deltas instead of imagery

### Takeaway
The best-documented cases are Mobileye REM (about 10 KB/km of extracted landmarks, 15M+ km mapped per day) and NDS.Live binary tile-delta streaming (up to 75% less transfer). Planet has shown onboard detection on Pelican-4 (Jetson Orin, aircraft found in about 0.5 s) with the stated aim of downlinking results instead of imagery. I found no published compression ratio specifically for **polygon deltas** of AI-derived masks.

### Cited Findings
- **Mobileye REM** vehicles upload anonymized landmark and road data at about **10 KB per km** — [Mobileye REM](https://www.mobileye.com/technology/rem/); [Mobileye blog](https://www.mobileye.com/blog/rem-mapping-avs/). From round 1: 8M+ vehicles, 54.7B km per year, 15M+ km mapped daily, and Stellantis integration from 2027 — [Mobileye REM](https://mobileye.com/technology/rem); [TimesLive Jul 2026](https://www.timeslive.co.za/motoring/2026-07-21-mobileye-to-supply-stellantis-with-its-cloud-based-driver-assistance-tech/)
- **NDS.Live delta updates.** The server diffs old and new tile versions at the **binary** level and the vehicle downloads only the delta, "cutting transfer volumes by as much as 75%" with no quality difference. Deltas are applied on the fly in the cache database without slowing streaming. Layers update at different rates (road layers monthly, speed limits weekly). — [NDS Association: delta updates](https://nds-association.org/nds-delta-updates-of-map-tile-service/); [NNG blog](https://nng.com/blog/low-cost-data-streaming-cloud-navigation)
- **Planet Pelican-4** (Mar 25 2026) ran an onboard NVIDIA Jetson Orin model that found more than a dozen aircraft in about 0.5 s, 500 km above Alice Springs. Planet frames this as sending "already processed and much lighter" products and plans to extend it to the Owl constellation for wildfire and disaster monitoring. — [SatNow](https://www.satnow.com/news/details/4986-planet-labs-demonstrates-ai-powered-onboard-object-detection-on-pelican-4-satellite); [Via Satellite, Apr 7 2026](https://www.satellitetoday.com/imagery-and-sensing/2026/04/07/planet-details-ai-driven-object-detection-onboard-pelican-4-satellite); [BusinessWire](https://secure.businesswire.com/news/home/20260407165913/en/Planet-Successfully-Runs-AI-in-Space)
- RTMap's crowdsourced multi-traversal map fusion gives no bandwidth or upload specification — [arXiv 2507.00980](https://arxiv.org/html/2507.00980v1)

### Inferences
- The NDS approach shows what deltas get you with **no** semantic stability: binary tile diffs save up to 75%. If tile contents were deterministic (unchanged tiles byte-identical), delta size would scale with the area that actually changed. Using the Overture measurement (about 2.4% of features touched per month), semantic deltas could be roughly 40x smaller than full tiles, compared with the 4x NDS reports. This is an extrapolation and has not been measured.
- Onboard polygonization is a precondition for "downlink vectors, not pixels" (Pelican, drones, wildfire aircraft). Deterministic integer output with no floating-point overlay is well suited to bit-exact reproducibility between edge and ground, which is the same property SweepLSD advertises for its FPGA/software parity.

### Gaps
- No Tesla fleet-map upload figures were found (the search budget ran out).
- I found no published comparison of imagery bytes vs. vector-delta bytes for drone, satellite or wildfire perimeter use cases.

---

## 5. Live digital twins: update-latency requirements for geometry, and whether polygon extraction is a stated bottleneck

### Takeaway
Twin vendors (NVIDIA Omniverse/Mega with KION and Accenture) describe real-time asset tracking and sub-second simulation but say nothing about polygon extraction. The clearest geometry-latency requirements come from public-safety and operations domains:
- **NWS Threats-in-Motion**: about 1-minute polygon updates.
- **ProbSevere**: 2-minute cycle, 2–5 minute latency.
- **Wildfire perimeters**: still manually digitized; the NIROPS daily progression dataset covers 12,705 perimeters across 737 fires, delivered up to 12 hours after the flight.

### Cited Findings
- **KION/NVIDIA Mega warehouse twins.** They map assets "in real time" and run simulations "in fractions of a second" to route AMRs and forklifts. The materials mention no raster-to-polygon step. — [Linde MH news](https://www.linde-mh.com/en/technical/News-Detail_5305536.html); [Omniverse warehouse digital twins docs](https://docs.omniverse.nvidia.com/digital-twins/latest/building-warehouse-digital-twins.html); [ERP Today](https://erp.today/kion-teams-with-nvidia-and-accenture-to-optimize-supply-chains-with-ai-powered-robots-and-digital-twins/)
- **Wildfire progression data.** A daily/sub-daily perimeter dataset for the Western US (June 2020–Nov 2024) holds 12,705 observations over 737 fires. "Infrared perimeters were manually digitized by trained infrared interpreters" from NIROPS flights. — [Mendeley Data 95rj5d379g](https://data.mendeley.com/datasets/95rj5d379g). From round 1: products reach crews "up to 12 hours after the original flight" — [USGS/AK Geo Summit](https://dggs.alaska.gov/webpubs/dggs/ago/doc)
- **NWS warning polygons.** TIM and FACETs/PHI aim for continuously updating warning polygons, roughly every minute (round 1) — [NOAA GSL TIM](https://gsl.noaa.gov/tim/); [CIRA TIM](https://www.cira.colostate.edu/?p=14346). ProbSevere outputs polygons about every 2 minutes with 2–5 minute observed latency — [NSSL EWP](https://inside.nssl.noaa.gov/ewp/?p=19142)
- **Robotics.** The SE(2) NavMesh updates polygonal regions online, and nvblox does GPU incremental mapping (177x faster surface reconstruction than CPU) but outputs SDF/mesh, not polygons (round 1) — [arXiv 2607.01454](https://arxiv.org/abs/2607.01454); [arXiv 2311.00626](https://arxiv.org/abs/2311.00626)
- Round 1 noted RAPIDS stopped publishing cuSpatial at v25.06 (June 2025), leaving no maintained GPU vector-geometry library — [RAPIDS RSN 45](https://docs.rapids.ai/notices/rsn0045/)

### Inferences
- Twin platforms store geometry as USD meshes and asset poses. Their latency problems sit in sensing and simulation, not in 2D polygonization. The workflows that need **minute-scale polygon versioning** are operational GIS (warnings, wildfire, flood), robotics navmeshes and annotation/VFX. The hypothesis should target those, not "digital twins" in general.
- Wildfire is the strongest single case for streaming plus versioning together: line-scan input, repeated passes, a perimeter that grows over time, and the need for change polygons (new growth since the last pass) with stable fire IDs. Today all of it is manual.

### Gaps
- I found no source naming polygon extraction as a digital-twin bottleneck, and no stated latency requirement for city-twin geometry updates.
- The operational status of TIM/PHI in 2026 is unconfirmed (testbed or deployed).

---

## 6. Verdict on the hypothesis ("git for geometry from AI" is under-served, and deterministic seam-free GPU polygonization enables it)

### Takeaway
**Partially supported, with one strong new piece of evidence.** In Overture's Sept 2026 release, 696.7M buildings (27.5%) were flagged geometry-changed. In sampled windows about 97% were identical shapes with a different ring start vertex. So **non-canonical polygon encoding already causes huge false diffs in the largest open AI-derived vector dataset**. Canonical, deterministic ring output is directly valuable. That value comes from a convention any tool could adopt (or a `ST_Normalize` step), not from GPU speed.
- **(c) Exact diffs / change polygons**: the strongest of the three claims, but only for same-grid raster sources. The real advantage is "XOR the bitmaps, then polygonize", not faster polygonization.
- **(b) Stable IDs**: demand is real and documented, but it is mainly a matching and conflation problem. Determinism helps (content hashes, no spurious changes) but does not solve it.
- **(a) Incremental re-polygonization of changed tiles**: easy and real, but weakly demanded on its own. The savings matter mainly because they make unchanged output byte-identical, which unlocks the downstream incremental infrastructure (Mapbox MTS, NDS deltas, changelogs) that today requires stable IDs.

### Cited Findings
- **Demand evidence for stable IDs and changelogs:**
  - Overture publishes a changelog every month and has suffered two ID breaks: June 2025 (format change) and Sept 2026 (100% base-theme churn plus a building_parts ID bug) — [Overture GERS blog](https://docs.overturemaps.org/blog/2025/06/25/getting-started-gers/); [Overture 2026-09-23 notes](https://docs.overturemaps.org/blog/2026/09/23/release-notes/)
  - Mapbox excludes data without stable IDs from incremental updates — [Mapbox MTS](https://docs.mapbox.com/mapbox-tiling-service/guides/incremental-updates/)
  - CVAT needs a consistent polygon start vertex and direction across frames — [CVAT](https://docs.cvat.ai/docs/annotation/manual-annotation/shapes/annotation-with-polygons/track-mode-with-polygons/)
- **How diffs are done today:**
  - Overture: an ID-keyed comparison of normalized geometry text, with no geometry in the changelog and no merge/split type — [Overture changelog docs](https://docs.overturemaps.org/gers/changelog/)
  - Planet: raster differencing, then vectorization — [ODSC/Planet](https://opendatascience.com/working-towards-planetary-scale-location-insights/)
  - Kart: row-level — [Kart table v3](https://docs.kartproject.org/en/latest/pages/development/table_v3.html)
- **Measured monthly building churn (Overture changelogs, DuckDB; see Section 3):**
  - July: 0.6% of features touched.
  - August: 2.4%, with removals about 5x in-place changes.
  - September: 27.5% flagged `data_changed`; in sampled windows about 97% of these were the same polygon with a rotated or reordered ring (spurious diffs).
- **Identity is computed by overlap or distance on rasters or points** in science trackers, and by polygon IoU ≥ 0.5 in Overture, with no GPU implementations found — [TempestExtremes poster](https://climate.ucdavis.edu/Ullrich_TempestPoster2019_v2.pdf); [tobac merge_split](https://tobac.readthedocs.io/en/latest/api/generated/tobac.merge_split.html); [Overture matching](https://labs.overturemaps.org/workshop/7-buildings-matching.html)

### Inferences
- **Evidence-backed:**
  1. Stable IDs and changelogs are a live, painful requirement at the largest open vector dataset.
  2. Downstream incremental infrastructure exists but excludes raster-derived layers that lack stable IDs.
  3. No GPU tooling exists for overlap tracking or merge/split event graphs.
  4. Change-polygon products are built by differencing rasters and then vectorizing.
- **Speculative (unvalidated product thesis):** a GPU "geometry versioning kernel" made of:
  1. per-tile deterministic canonical rings plus content hashes, so unchanged tiles are byte-identical;
  2. label-overlap matrices between versions, giving birth/death/continue/merge/split events and Euler-number topology events;
  3. signed change polygons from bitmap XOR;
  4. output as GeoParquet/Iceberg v3 rows with ID lineage.
  
  This would fill a gap no found product covers. The most credible first buyers are repeat-pass EO change feeds (Planet-style, Google 2.5D Temporal rasters), wildfire progression, and weather-object tracking (TempestExtremes/tobac users moving to GPU ML ensembles).
- **Risks:**
  - Identity across *different* sources or models (Overture's main churn driver) cannot be fixed by deterministic polygonization, because shapes genuinely differ.
  - Exact diffs need a shared pixel grid.
  - Competing in-model identity (SAM 3 memory) may make geometry-based identity redundant for video.
  - I found no buyer quote naming polygon-diff compute as a bottleneck.

### Gaps
- No direct customer or engineering-blog evidence that geometry diff compute (as opposed to ID stability) is a cost center.
- No public measurements of how much of Overture's building churn comes from geometry jitter vs. source swaps vs. real change. This could be measured by joining the changelog to two releases and computing IoU of removed/added pairs, which is a recommended next step.
- The round-2 search budget was exhausted, so Esri change detection, Dynamic World, GeoGig/Dolt, Tesla fleet maps, Deep Snake video and rotoscoping tools were not covered.
