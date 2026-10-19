# Dynamic / Temporal Geometry: closing the round-2 gaps (October 2026)

Scope: this round covers the items the round-2 notes (`research_notes/GPU delineation use case map/r2_dynamic_geometry.md`) did not reach. The question is whether a deterministic, canonical GPU polygonizer opens real uses in three areas:
- stable ring hashes;
- exact change polygons (polygonize the XOR / AND-NOT of two bitmaps);
- temporal-consistency rewards for video segmentation and tracking RL (the same object gets the same canonical outline across frames).

Round-2 facts are not repeated except where needed: Overture's ~97% spurious "geometry changed" diffs from ring rotation, Planet's change pipeline (difference rasters, then vectorize), Kart's row-per-file storage, and Mobileye's 10 KB/km.

Method note: about 25 web searches and fetches. Some items stop at search-snippet level and are marked that way. The GeoGig, Dolt, DETER-minimum-area and incremental-contour items rest on fewer sources than the rest.

---

## 1. Change-detection products: how they produce change polygons, and whether they ship vectors

### Takeaway
Operational change products are **pixel-first**. GLAD, RADD and the GFW integrated alerts, Dynamic World, and the Esri difference rasters all work on pixels. Polygons appear in three ways:
- an optional vectorization step: the Esri Change Detection Wizard's "polygon feature class" output, GFW's .gpkg export, Planet's "vector simplification";
- coarse grid cells: Planet's 8x8-pixel cells;
- in one case, a polygon-native alert product: INPE DETER's daily polygons.

None advertises exactness or reproducibility of the polygons. Minimum mapping units and event grouping are coarse:
- DETER-B: >3 ha;
- Planet: 8x8-pixel grid cells;
- GFW: 180-day event merging.

### Cited Findings
- **ArcGIS Pro.** The Change Detection Wizard takes two categorical or continuous rasters, or a multidimensional time series. It outputs a difference raster, a raster function template, **or a polygon feature class**. Other change methods (Compute Change, CCDC, LandTrendr through "Detect Change Using Change Analysis Raster") output rasters. The page I fetched does not mention deep-learning change detection. — [Esri: Change detection in ArcGIS Pro](https://doc.esri.com/en/arcgis-pro/latest/help/analysis/image-analyst/change-detection-in-arcgis-pro.html)
- **GFW / Global Nature Watch Integrated Disturbance Alerts** (post dated Jan 13 2026):
  - Merges GLAD-L (30 m), GLAD-S2 (10 m), RADD (10 m) and DIST-ALERT (30 m) into a layer resampled to **10 m raster**. The post mentions no vector polygons.
  - Optical detection runs "every one-four days" in the dry season. Radar detection "may take up to 12 days".
  - The integrated layer detects alerts "on average nine days earlier" than single systems.
  - Alerts within **180 days** of each other count as **one event**.
  
  — [Global Nature Watch blog](https://globalnaturewatch.org/blog/data-and-tools/integrated-deforestation-alerts/)
- **GLAD / RADD encoding.** Each alert is a pixel. One integer holds confidence (2 = low, 3 = high) and days since 2014-12-31. Alerts can be exported as **.gpkg** with geometry and properties. GLAD lag ranges from 3 days to several months depending on cloud. Sentinel-1 revisit for RADD is 6–12 days. — [GFW Open Data: RADD](https://data.globalforestwatch.org/datasets/gfw::deforestation-alerts-radd/about); [GFW integrated alerts dataset](https://data.globalforestwatch.org/datasets/gfw::integrated-deforestation-alerts/about); [GNW: GLAD alerts explained](https://globalnaturewatch.org/blog/data-and-tools/glad-deforestation-alerts/) (search-snippet level)
- **INPE DETER** is polygon-native. It produces **daily sets of polygons** ("alerts") that outline deforestation or degradation in the Amazon and Cerrado, distributed as shapefile/CSV through TerraBrasilis. Each alert has an image date and a change type. Minimum area: DETER-A (MODIS, 250 m, 2004–2017) >25 ha; DETER-B (CBERS-4 WFI, 64 m, 2014–present) **>3 ha**. — [datazoom.amazonia README (CRAN)](https://cran.r-project.org/web/packages/datazoom.amazonia/readme/README.html) (secondary source); [TerraBrasilis FOSS4G 2024 talk](https://talks.osgeo.org/foss4g-2024/talk/VUYRLW/). I did not verify from a primary INPE source whether DETER polygons are drawn by visual interpretation or by automatic delineation.
- **MapBiomas Alerta** validates and refines alerts from DETER, GLAD and others into polygon reports. The method page was found but not fetched. — [MapBiomas Alerta method](https://alerta.mapbiomas.org/en/metodo-mapbiomas-alerta/) (search-snippet level)
- **Google Dynamic World**:
  - Near-real-time 10 m land cover from Sentinel-2: nine class probabilities, updated with each sufficiently cloud-free scene (up to every ~5 days).
  - Change detection is left to the user, by rules on class-probability changes in Earth Engine.
  - I found **no official Google change-polygon product**. Third-party code exists: `rw-dynamicworld-cd` on PyPI and SkyTruth tooling.
  
  — [Earth Engine tutorial: Dynamic World pt 3](https://developers.google.com/earth-engine/tutorials/community/introduction-to-dynamic-world-pt-3); [Land & Carbon Lab: Dynamic World](https://landcarbonlab.org/data/dynamic-world); [PyPI rw-dynamicworld-cd](https://pypi.org/project/rw-dynamicworld-cd/); [SkyTruth Dynamic World](https://skytruth.org/section/dynamic-world/)
- **Planet** (round 2, repeated for completeness): Road & Building Change Detection outputs GeoJSON vectors on a "standard 8x8 pixel gridcell" with a z-score threshold of 10.5. The pipeline differences successive rasters and then applies "vector simplification". — [Planet techspec](https://docs.planet.com/data/analytic-feeds/road-building-change-detection/techspec/); [ODSC/Planet](https://opendatascience.com/working-towards-planetary-scale-location-insights/)

### Inferences
- Every product found follows the same pattern: decide change per pixel, then optionally vectorize. That pattern is exactly "polygonize(A AND NOT B)". So the canonical polygonizer drops into the existing pipeline position with no workflow change. The advantage it offers (no slivers, deterministic output, byte-stable reruns) is something **no vendor markets or measures**.
- Products deliberately coarsen geometry: 3 ha minimum (DETER-B), 8x8-pixel cells (Planet), 180-day event merging (GFW). For these users, **reproducibility and identity** (the same alert polygon on reprocessing, a stable event ID) is a more plausible pitch than **sub-pixel exactness**.
- GFW ships its core product as raster and offers vector only as an export. So deforestation buyers treat polygons as a convenience format, not the system of record. DETER is the exception and is the most polygon-dependent public alert system.

### Gaps
- Pricing and latency for Esri change detection, and for Planet Analytic Feeds beyond weekly/monthly cadence: not found.
- How DETER polygons are delineated (manual vs. automatic) and how they are reconciled day to day: not verified from INPE primary docs.
- Whether EUDR (EU Deforestation Regulation) due-diligence workflows require reproducible or auditable change polygons: not searched in this round. It is a plausible demand driver and should be checked next.

---

## 2. Geometry versioning tools: GeoGig, Dolt, Kart, Iceberg. How diffs are computed

### Takeaway
All of them diff **records**, not **shapes**:
- **GeoGig** is effectively dead. LocationTech listed it for termination review for inactivity in Aug 2026.
- **Dolt** stores MySQL spatial types and diffs at row/cell level.
- **Kart** diffs per feature by primary key, emitting +/- (insert/update/delete) of whole hexWKB geometries.
- **Iceberg v3** versions snapshots (round 2).

None computes geometric change polygons between versions, and none canonicalizes rings. A ring rotated between versions is therefore an "update" in all of them, the same failure mode measured in Overture in round 2.

### Cited Findings
- **GeoGig.** In a LocationTech PMC thread (Jody Garnett, **Aug 10 2026**), GeoGig is listed with SFCurve, Spatial4j and uDig as facing "Termination review due to no activity". The thread records no final decision. Latest GitHub release seen: 1.2.0. — [locationtech-pmc message](https://www.eclipse.org/lists/locationtech-pmc/msg01252.html); [LocationTech GeoGig governance](https://projects.eclipse.org/projects/locationtech.geogig/governance)
- **Dolt.** Spatial types arrived in beta in Feb 2022 (Point, LineString, Polygon). From v0.50.6 Dolt supports every MySQL type, including the Multi* geometries. Dolt's diff model is row/cell-level, which is general Dolt behavior; I found no geometry-specific diff. — [DoltHub: Introducing Spatial Types](https://dolthub.com/blog/2022-02-09-spatial-types); [DoltHub: Dolt supports every type](https://dolthub.com/blog/2022-11-09-dolt-supports-every-type)
- **Kart diff.** Version control is at "row and cell level". `diff -o json` is organized as `<dataset>/<feature|meta>/<pk>/<+|->/<feature>`. `diff -o geojson` emits a flat array with IDs `<+/-><U|I|D>::<pk>`. Geometries are hexWKB because GeoJSON cannot hold curved types. — [Koordinates: Kart](https://koordinates.com/products/kart/); [kart issue #71 (sno-era diff output design)](https://github.com/koordinates/kart/issues/71) (search-snippet level, sno = Kart's former name)
- **Iceberg v3 / Parquet GEOMETRY** gives snapshot versioning with bbox stats. I found no shape-level diff feature (round 2). — [Wherobots: Iceberg geospatial v3](https://wherobots.com/blog/iceberg-geospatial-types-v3/)

### Inferences
- Every "git for geo" tool keys diffs on a primary key and compares bytes or values. A **content hash of a canonical ring** (fixed orientation, lexicographically smallest start vertex, integer grid) would give these tools free "unchanged" detection, and would let them pair "deleted + inserted" rows into "same shape, new ID". That is a convention plus a hash, not GPU work. The GPU engine's role is emitting canonical output at the source, so downstream tools never see rotated rings.
- No tool offers "show me the area that changed". Signed change polygons (A AND NOT B, B AND NOT A) would be a new **geometry-level diff column** on top of Kart, Dolt or Iceberg. I found no competitor.
- GeoGig's dormancy and Kart being the only active geo-VCS suggest a small market for geo version control as such. The value has to ride inside existing pipelines (GeoParquet/Iceberg, Overture-style changelogs), not a new VCS.

### Gaps
- Kart's current (v0.12+) diff output schema was not fetched directly. The details above come from a sno-era issue snippet.
- I could not confirm whether Dolt does anything geometry-aware in diffs, such as spatial indexes in diff queries.

---

## 3. Fleet-crowdsourced maps: Tesla and Mobileye REM geometry deltas

### Takeaway
**Tesla has published nothing** on geometry deltas. The only public detail comes from reverse-engineering: speed limits and real-time speeds are uploaded about every 2 minutes, and map corrections are downloaded. **Mobileye REM** is better documented:
- about 10 KB/km of Road Segment Data, compared on the cloud side against the existing Roadbook, with updates only "where changes are detected";
- no public spec of the delta format or of how geometry identity is matched.

### Cited Findings
- **Tesla.** Per hacker @greentheonly (reported Apr 3 2023), vehicles "send back detected speed limits and real-time speeds to the Tesla maps service every two minutes". The car fetches "map corrections, such as the locations of stop signs, crosswalks". The article gives no format, size or geometry-delta details. — [Not a Tesla App, Apr 2023](https://www.notateslaapp.com/news/1311/teslahacker-reveals-how-maps-are-augmented-with-fleet-data-between-updates). Tesla's 2017 "Vector Maps" added lane details on the display side — [Electrek 2017](https://electrek.co/2017/11/12/tesla-maps-smoother-detailed-engine/). I found no Tesla engineering publication on fleet geometry deltas.
- **Mobileye REM** (blog, Mar 11 2021):
  - Road Segment Data is "only 10 kilobytes … per kilometer", sent over existing cellular modems (3G is enough).
  - Roadway "signatures" are "compared automatically to the existing map, and where changes are detected, the Roadbook is updated", instead of remapping whole sections.
  - All of Japan (25,000 km) was mapped "in only 24 hours". The database is **400 MB**, which works out to about 16 KB/km for the distributed map.
  - At the time, 65M+ vehicles carried Mobileye tech.
  
  — [Mobileye blog: REM mapping](https://www.mobileye.com/blog/rem-mapping-avs/). Roadbook tiles are distributed in small pieces and "partially updated in real-time" — [Mapbox blog: Mobileye Roadbook](https://blog.mapbox.com/mobileye-roadbook-high-precision-hd-maps-distributed-at-scale-9b4e692d29f) (search-snippet level)

### Inferences
- Both fleets upload **pre-extracted semantic primitives** (landmarks, lane splines, speed values), not polygons of raster masks. So the "exact change polygon" idea has **little purchase in AV mapping**. Their change unit is a matched element: lanes and signs, the same atomic-element diffs as ArgoTweak and RTMap in round 2. Canonical polygonization would matter only for area features such as drivable area, parking lots and construction zones. I found no evidence these are crowd-updated as polygons.
- Tesla's map stack is opaque. Any claim about Tesla geometry deltas would be speculation.

### Gaps
- No public REM delta format, no change-detection thresholds, no figure for time-to-reflect-reality.
- No Tesla primary source.

---

## 4. Video: contour methods, rotoscoping tools, temporal-consistency metrics, and RL rewards for temporal consistency

### Takeaway
- **Production rotoscoping has split into two camps.**
  - Mask-propagation tools that output mattes: Resolve Magic Mask 2, AE Roto Brush 3, SAM 2/3-style trackers.
  - **Spline-propagation tools** that output editable vector shapes: Foundry SmartRoto (Jul 2026, $599/yr) and Mocha Pro's Object Brush (ML mask converted to editable splines, then tracked).
  
  The spline camp exists because artists need vertex-level continuity across frames. That is precisely the canonical-outline property.
- **Temporal-consistency metrics exist.** The label-free, optical-flow-warped mTC metric dates from 2020.
- **No published RL reward measures a model's own cross-frame consistency of segmentation masks or contours.** The 2025–2026 RL-for-VOS papers reward either IoU against ground truth propagated across frames (Veason-R1's "consistency" reward) or temporal localization (keyframes/time ranges).

### Cited Findings
**Contour methods**
- **Deep Snake** (Jan 2020) is image-only. It deforms an initial contour with circular convolution and runs at 32.3 fps at 512x512 on a 1080Ti. Searches for 2024–2025 video/temporal Deep Snake variants found nothing. — [arXiv 2001.01629](https://arxiv.org/abs/2001.01629)
- Polygonal Point Set Tracking (CVPR 2021) and Coherent Parametric Contours were covered in round 2. — [arXiv 2105.14584](https://ar5iv.arxiv.org/html/2105.14584)

**Rotoscoping tools**
- **Foundry SmartRoto** (released **Jul 21 2026**):
  - The artist draws roto shapes on initial frames and keys frames where significant changes occur. The model then "automatically moves and deforms existing roto splines to follow the subject".
  - Output is editable splines inside Nuke's Roto node, not mattes. It uses "the relationship between splines as context", and artists can key non-sequentially.
  - $599/year on top of Nuke 16.1+. Runs locally. Trained on licensed data over four years.
  
  — [CG Channel, Jul 2026](https://www.cgchannel.com/2026/07/foundry-releases-new-ai-powered-smartroto-add-on-for-nuke/); [Foundry SmartRoto](https://www.foundry.com/products/nuke-family/smartroto). Foundry markets it as keeping "complete vector control over every single point, handle, and curve" instead of "a flat black-and-white image mask" (search-snippet level).
- **Community matte-to-spline bridge.** AlphaToRoto is an open-source Nuke tool that traces alpha mattes into animated Roto shapes using embedded potrace. — [GitHub petermercell/AlphaToRoto](https://github.com/petermercell/AlphaToRoto)
- **Mocha Pro 2025 Object Brush** isolates an object with a click, then converts the selection into **editable splines** that are tracked with planar or PowerMesh trackers. — [Boris FX press](https://blog.borisfx.com/press/mocha-pro-unveils-ai-powered-rotoscoping-and-masking); [CG Channel, Dec 2024](https://www.cgchannel.com/2024/12/boris-fx-releases-mocha-pro-2025/)
- **DaVinci Resolve 20.1 Magic Mask v2** (Aug 2025) is point-prompted and AI-propagated. Forum threads report tracking failures and Fusion feature gaps. — [Digital Production, Aug 2025](https://digitalproduction.com/2025/08/07/davinci-resolve-20-1-apple-vision-pro-magic-mask-v2-and-a-few-more-reasons-to-backup/); [Blackmagic forum](https://forum.blackmagicdesign.com/viewtopic.php?f=21&t=222790)
- **AE Roto Brush 3** produces matte-based output with a refine-edge tool. — [Boris FX: top AI roto tools](https://borisfx.com/blog/top-6-ai-rotoscoping-tools-free-and-paid/) (search-snippet level)
- Runway: no specifics found.

**Temporal-consistency metrics**
- **mTC** (Varghese et al., CVPRW 2020) is an unsupervised metric. TC_t is the mIoU between the prediction at t and the prediction at t-1 warped by optical flow (FlowNet2). It needs no labels. — [CVPRW 2020 paper](https://openaccess.thecvf.com/content_CVPRW_2020/papers/w20/Varghese_Unsupervised_Temporal_Consistency_Metric_for_Video_Segmentation_in_Highly-Automated_Driving_CVPRW_2020_paper.pdf). A related "perceptual consistency" metric is in [WACV 2022](https://openaccess.thecvf.com:443/content/WACV2022/papers/Zhang_Perceptual_Consistency_in_Video_Segmentation_WACV_2022_paper.pdf), and a 2026 spatio-temporal attention VSS paper targets consistency in driving — [arXiv 2602.10052](https://arxiv.org/pdf/2602.10052)

**RL for video segmentation and tracking**
- **Veason-R1** (Aug 2025; CVPR 2026 poster) uses Qwen2.5-VL plus GRPO. Its reward is R = R_format + R_keyframe + R_spatial (box IoU, Hungarian) + R_u, all weights 1.0.
  - The "unified consistency" reward R_u = (1/T) Σ IoU(m_t, m_t^gt). Masks are produced by **frozen SAM2** from the predicted keyframe box, so it is **ground-truth IoU over propagated frames, not self-consistency**.
  - Results: ReasonVOS 59.9 J&F (+10.0 over GLUS-7B); MeViS 52.2 zero-shot.
  
  — [arXiv 2508.11538](https://arxiv.org/html/2508.11538); [CVPR 2026 poster](https://cvpr.thecvf.com/virtual/2026/poster/38001)
- **VideoSeg-R1** (Nov 2025) uses Qwen2.5-VL with SAM2 and XMem propagation. Rewards are format, mask IoU, point/box L1, non-repetition and a soft length penalty. There is **no explicit temporal-consistency term**; coherence comes from XMem propagation. Results: Ref-YouTube-VOS 81.3, Ref-DAVIS17 79.8, MeViS 55.3, ReVOS 61.1 J&F. — [arXiv 2511.16077](https://arxiv.org/html/2511.16077)
- **VideoSEG-O3** (Jun 5 2026) uses Qwen3-VL plus SAM2 with multi-turn RL.
  - Its "temporal" reward is keyframe membership plus time-range precision |P∩G|/|P|.
  - Its "progressive" reward rewards IoU gains over the historical maximum across turns, not across frames.
  - Results: ReVOS 70.3, MeViS 60.0, Ref-SAV 65.5 J&F (4B model).
  
  — [arXiv 2606.06819](https://arxiv.org/html/2606.06819)
- **ReVSeg** (Dec 2025) splits VOS into interpretation, temporal evidence selection and spatial grounding, and rewards intermediate steps. — [arXiv 2512.02835](https://arxiv.org/html/2512.02835) (search-snippet level)
- **SAM2RL** (Jul 2025) uses RL to control SAM 2's memory-bank updates. Per-video overfit agents beat heuristics by more than 3x the heuristics' gains. The abstract does not state the reward. — [arXiv 2507.08548](https://arxiv.org/abs/2507.08548)
- **Video generation.** RL here *does* use temporal-geometry consistency rewards. VGGRPO (Mar 2026) adds a camera-motion-smoothness reward and a geometry reprojection-consistency reward against "geometric drift, temporal flicker". — [arXiv 2603.26599](https://arxiv.org/html/2603.26599v2)

### Inferences
- **Gap confirmed.** The label-free consistency metric (mTC, 2020) and RL-for-VOS (2025–2026) both exist, but I found no paper that combines them into a reward. Nor did I find any reward on **outline/contour identity**: same canonical ring, same vertex count, same hash for an unchanged object, or Hausdorff/turning-function stability of the polygon. The nearest analogues:
  - Veason-R1's R_u, which is GT-based;
  - VGGRPO's geometry consistency, which is for generation, not segmentation.
- **Where a canonical polygonizer actually helps a reward.**
  - A mask-level consistency reward (flow-warped IoU, i.e. mTC) needs no polygons at all.
  - Polygons add **contour-shape stability** (boundary jitter, vertex churn, topology flips such as a hole appearing or disappearing via the Euler number) and **cheap exact change regions** (XOR of warped vs. current mask, polygonized for a "where did it flicker" signal).
  - Determinism matters because a reward must be bit-reproducible across GRPO rollouts. Canonical rings make "unchanged" a hash equality instead of a thresholded IoU.
  - This is a plausible, unclaimed niche, but it is a contribution to **reward design**, not a GPU throughput story. Mask-to-ring cost on GPU is small next to SAM2/VLM rollout cost.
- **Demand for vector continuity in VFX** is real and commercial (SmartRoto at $599/yr; Mocha's mask-to-spline step). But those tools want **artist-editable, sparse splines with point correspondence**, not dense integer pixel-edge rings. A canonical pixel-edge ring would need simplification plus correspondence (as in CVAT's "same start vertex, same direction" rule from round 2) to serve this market.

### Gaps
- I could not read SAM2RL's reward definition (abstract only).
- No Runway, AE Roto Brush 3 or Resolve documentation on how temporal consistency is handled internally.
- No temporal-contour Deep Snake successor found. "Contour-based video instance segmentation 2024–2026" deserves one more targeted search.
- I did not check whether SAM 3 / SAM 3.1 papers report a temporal-stability metric.

---

## 5. Incremental / dynamic marching squares and dynamic connected components (especially GPU)

### Takeaway
The literature has three kinds of incremental work:
- **graph** connectivity on GPU (incremental, insert-only, very fast);
- **contour-tree** maintenance on CPU (theory);
- **patents on regenerating marching-squares/cubes output only in changed blocks**.

I found **no GPU algorithm that maintains image connected components or polygon rings under pixel flips** (insertions and deletions), and no "incremental marching squares" paper. The practical method is block-local recompute, which the halo-tiled engine already does.

### Cited Findings
- **GPU incremental graph connectivity**: up to **48.23 billion edges/s** on a Tesla V100, 1.85–13.36x over a 72-core CPU. This is for edge insertions; deletions (decremental) are the hard case. — [arXiv 2008.11839 (Exploring the Design Space of Static and Incremental Graph Connectivity Algorithms on GPUs)](https://ar5iv.labs.arxiv.org/html/2008.11839)
- **Incremental CC for graph streams on GPU** (Electronics 2023), round 2 — [MDPI 12(6):1465](https://www.mdpi.com/2079-9292/12/6/1465)
- **DynLP** (arXiv 2604.06596, Apr 8 2026) does GPU dynamic batch label propagation for graph semi-supervised learning: 13x average and up to 102x speedup. It is **not** image CCL. — [arXiv 2604.06596](https://arxiv.org/abs/2604.06596)
- **Contour trees of dynamic terrains.** This work maintains contour trees of time-varying 2-manifolds under vertex insertion, deletion and edge flip, with O(log n) certificate processing. — [arXiv 1406.4005](https://arxiv.org/abs/1406.4005)
- **Selective regeneration patent.** It keeps contour segments for voxel/cell subsets whose SDF did not cross the threshold, and re-runs Marching Squares only for changed subsets. — [US10366534 "Selective surface mesh regeneration for 3-dimensional renderings"](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10366534) (search-snippet level; assignee not verified). A similar real-time contour-line patent for dynamic map displays — [US7006106B2](https://patents.google.com/patent/US7006106)
- **Quadtree incremental CC.** Incremental maintenance of connected components for binary images in **quadtree** form, visiting fewer nodes than batch recompute (search-snippet level; the exact source in the result list was not verified).
- GPU CCL itself remains static-image work. Examples: Allegretti et al. (IPAS 2018) and LIP6 CCA (ICASSP 2021). — [IPAS 2018 GPU CCL](https://iris.unimore.it/bitstream/11380/1167089/4/2018_IPAS_Optimizing_GPU_Based_Connected_Components_Labeling_Algorithms.pdf); [ICASSP21 CCA GPU](https://largo.lip6.fr/~lacas/Publications/ICASSP21_CCA_GPU.pdf)

### Inferences
- **"Dirty-tile recompute with halos" is the state of practice**, and it is what the selective-regeneration patents describe. The engine needs no new algorithm for incremental *extraction*.
- **The unsolved part is global**: re-stitching rings that cross dirty/clean boundaries and handling deletions, i.e. component splits under pixel removal. This is the classic hard decremental-connectivity problem. A GPU implementation that handles splits for image components, and keeps canonical ring IDs while doing so, would be novel. Even so, the case for doing it rather than re-running a full-frame GPU pass (already milliseconds) is weak except for very large rasters such as continental mosaics.

### Gaps
- No benchmark comparing full-frame GPU re-polygonization with incremental update at video or continental scale.
- Decremental (split-handling) image CCL on GPU: no paper found.

---

## 6. Verdicts: (1) demand for exact change polygons; (2) temporal-consistency rewards in RL

### Takeaway
**(1) Demand for change polygons is real. Demand for *exact* change polygons is not documented.** Products ship coarse, pixel-derived, optionally vectorized change. The value of a canonical polygonizer is **reproducibility and identity**, not precision. **(2) RL rewards for temporal consistency of segmentation outputs essentially do not exist** in self-supervised form. The closest is Veason-R1's GT-IoU-across-frames "consistency" reward. This is an open, defensible niche.

### Cited Findings
- **Change polygons are produced and consumed:** DETER daily polygons ([CRAN datazoom README](https://cran.r-project.org/web/packages/datazoom.amazonia/readme/README.html)), Esri wizard polygon output ([Esri](https://doc.esri.com/en/arcgis-pro/latest/help/analysis/image-analyst/change-detection-in-arcgis-pro.html)), Planet GeoJSON change vectors ([Planet techspec](https://docs.planet.com/data/analytic-feeds/road-building-change-detection/techspec/)), GFW .gpkg export ([GFW RADD](https://data.globalforestwatch.org/datasets/gfw::deforestation-alerts-radd/about)).
- **But at coarse resolution:** >3 ha (DETER-B), 8x8-pixel cells (Planet), 180-day event merging (GFW). The flagship integrated product is raster-only ([GNW](https://globalnaturewatch.org/blog/data-and-tools/integrated-deforestation-alerts/)).
- **Versioning tools diff records by key, not shapes:** Kart ([Koordinates](https://koordinates.com/products/kart/)), Dolt ([DoltHub](https://dolthub.com/blog/2022-11-09-dolt-supports-every-type)). GeoGig is under termination review ([LocationTech PMC, Aug 2026](https://www.eclipse.org/lists/locationtech-pmc/msg01252.html)).
- **RL-for-VOS rewards** are GT-based or about temporal localization: Veason-R1 ([arXiv 2508.11538](https://arxiv.org/html/2508.11538)), VideoSeg-R1 ([arXiv 2511.16077](https://arxiv.org/html/2511.16077)), VideoSEG-O3 ([arXiv 2606.06819](https://arxiv.org/html/2606.06819)). A label-free consistency metric exists ([Varghese CVPRW 2020](https://openaccess.thecvf.com/content_CVPRW_2020/papers/w20/Varghese_Unsupervised_Temporal_Consistency_Metric_for_Video_Segmentation_in_Highly-Automated_Driving_CVPRW_2020_paper.pdf)). Video *generation* RL already uses geometric temporal-consistency rewards ([VGGRPO](https://arxiv.org/html/2603.26599v2)).
- **Commercial demand for vector continuity across frames:** SmartRoto ([CG Channel](https://www.cgchannel.com/2026/07/foundry-releases-new-ai-powered-smartroto-add-on-for-nuke/)) and Mocha Object Brush ([Boris FX](https://blog.borisfx.com/press/mocha-pro-unveils-ai-powered-rotoscoping-and-masking)).

### Inferences
- **Verdict (1): weak-to-moderate for "exact"; moderate for "deterministic and identity-stable".**
  - Nobody asks for exact change polygons by name.
  - The concrete pains found across rounds are spurious diffs (Overture's ~97% ring-rotation false changes in round 2) and ID churn.
  - A canonical engine fixes both cheaply.
  - Best pitch: "byte-stable, canonical, sliver-free change polygons and ring hashes for repeat-pass rasters". Target DETER/MapBiomas-style alert vectorization, Planet-style feeds and Overture-style changelogs. EUDR auditability is a hypothesis to test.
- **Verdict (2): the reward does not exist yet, so this is a research opportunity rather than a market.**
  - Proposed reward: a label-free, deterministic temporal-consistency reward combining:
    - flow-warped mask IoU (mTC);
    - canonical-ring hash equality or contour distance for objects judged static;
    - Euler-number/topology flip penalties;
    - polygonized XOR area as a "flicker area" penalty.
  - It could plug into GRPO pipelines like Veason-R1's or VideoSEG-O3's, which already run SAM2 inside the reward loop.
  - Risks:
    - (a) reward hacking: a model that outputs static or over-smoothed masks scores perfectly, so the term must be paired with GT IoU or a motion-aware warp;
    - (b) mask-level mTC already captures most of the signal, so the polygon contribution must be shown empirically (e.g., boundary-jitter reduction that mTC misses);
    - (c) rollout cost is dominated by the VLM and SAM2, so GPU polygonization speed is not the selling point. Determinism and exact topology are.

### Gaps
- No user or buyer quote asking for "exact" change polygons.
- No paper testing whether contour-level consistency adds signal beyond mask-level mTC.
- EUDR / compliance-driven reproducibility requirements not researched.
