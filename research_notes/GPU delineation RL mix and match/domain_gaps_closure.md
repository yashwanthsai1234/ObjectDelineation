# Domain Gap Closure: items round 1 left uncovered (GPU mask→polygon engine; RL/ML angle)

Scope: these items were listed as "not researched" in the Gaps sections of `GPU delineation use case map/{fabrication_physical, geospatial_eo, biomedical, robotics_av}.md`. This is a breadth pass of about 28 tool calls, so depth per item is uneven. Each item covers how rasters become polygons, scale, tools, pain points, and the ML/RL loop with a possible place for geometric rewards. Engine shorthand: **rings+holes**, **seams** (tile-seam-free), **batch-θ** (batched thresholds), **topo** (component/hole counts), **canon** (canonical encoding).

---

## Q1. Fabrication / physical items

### Takeaway
The best new fits are **mining fragmentation** and **P&ID line graphs**. In fragmentation, touching-fragment topology directly changes the output metric, a size distribution. In P&ID digitization, ML already trains on graph-level targets. Oil spill already has an operational *polygon* product (CleanSeaNet), but it is latency-bound at tens of minutes and does not run at GPU scale. Power-line and wind-blade work mostly stays in point clouds or bounding boxes. SDF fonts are the *reverse* direction (vector→SDF), except for neural-font vectorization.

### Cited Findings

**Oil spill (SAR)**
- EMSA CleanSeaNet (since 2007) delivers oil-spill products whose attributes include a *polygon geometry describing the spill and associated slicks*, center lat/lon, area in m², timestamp, and an oil-confidence class. Nominal delivery is about 20 min for SAR and 40 min for optical. Shapefile and alert delivery is 30 min per 400 km of scene, plus 5 min per extra 200 km — [Wikipedia: CleanSeaNet](https://en.wikipedia.org/wiki/CleanSeaNet) (secondary; consistent with [EMSA/Interspill deck](https://www.interspill.org/downloads/archive/CleanSeaNet-10-years-S-Djavidnia-European-Maritime-Safety-Agency-EMSA.pdf), not fetched).
- Research models are semantic-segmentation masks. DeepSegFusion (SegNet + DeepLabV3+ with attention fusion, arXiv Jan 2026) reports 94.85% accuracy but only **IoU 0.5685** on SAR oil-spill sets including ALOS PALSAR — [arXiv 2601.12015](https://arxiv.org/abs/2601.12015v1). Large-scale detection and categorization with DL dates from 2020 — [arXiv 2006.13575](https://arxiv.org/pdf/2006.13575).

**Power-line vegetation encroachment**
- The dominant pipeline is LiDAR point-cloud semantic segmentation (e.g., RandLA-Net for vegetation, poles and wires), followed by clearance distance checks. The stated challenge is the volume of point data — [ISPRS Annals 2026](https://isprs-annals.copernicus.org/articles/X-3-W4-2025/227/2026/isprs-annals-X-3-W4-2025-227-2026.pdf).
- An image-based variant uses a hybrid U-Net + SegFormer combined via meta-learning for powerline/vegetation detection and encroachment — [ETS Montréal](https://pure.etsmtl.ca/en/publications/a-novel-hybrid-ml-approach-for-powerline-vegetation-encroachment-).
- A patent, "System and method of intelligent vegetation management", identifies *encroachment zones* and *boundaries of identified vegetation* from mixed satellite, aircraft and drone time series, then scores criticality by distance to lines — [US 11842538](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11842538) (assignee not verified).
- The GIS deliverable is per-tree *crown polygons* and treetop points, attributed with the worst clearance violation. Cage-clearance results are reported by distance zone — [NV5 TVMP flyer](https://www.nv5.com/wp-content/uploads/2021/02/NV5-Geospatial_TVMP-flyer.pdf) (via search snippet; not fetched).

**Wind turbine blade inspection**
- The public data is mostly **bounding boxes**. The DTU drone dataset has about 13k RGB images at 586×371 with YOLO boxes (damage, dirt). The DTU Risø 2024 set has 29 DJI Mavic 2 videos. YOLOv8–v11 and Faster R-CNN dominate — [Boxes2Pixels arXiv 2604.11162](https://arxiv.org/pdf/2604.11162); [Mendeley DTU Risø](https://data.mendeley.com/datasets/6nzbdvjn87).
- *Boxes2Pixels* (2026) learns defect **segmentation from noisy SAM masks** prompted by boxes. Turning boxes into masks is an active ML problem — [arXiv 2604.11162](https://arxiv.org/pdf/2604.11162).

**Mining: rock fragmentation and ore sorting**
- Post-blast fragmentation: YOLO12l-seg was fine-tuned on >500 images and reaches mask mAP@0.5 ≈ 0.80 at **~15 FPS**. Masks are converted to 3D coordinates for spatial descriptors (PCA, KDE hotspots, Delaunay edge stats) — [arXiv 2507.20126](https://arxiv.org/html/2507.20126v1).
- Stated pain point: quarry scenes contain "thousands of irregular, densely packed fragments with touching boundaries". Small boundary errors **merge or split fragments** and so change fragment counts and the particle size distribution (P80, oversize, fines) — [RIT structure-aware segmentation](https://www.rit.edu/dirs/research/structure-aware-segmentation-post-blast-rock-fragmentation-analysis); [Roboflow](https://roboflow.com/ai/rock-fragmentation-analysis).
- FragSAM (2026) claims near-real-time fragment segmentation that works across blasting, TBM muck and coastal rock without retraining — [IJMST](https://journal.hep.com.cn/ijmst/EN/10.1016/j.ijmst.2026.03.017).
- Ore sorting: a CNN does semantic segmentation of particles in flight and maps the results onto **pneumatic ejection valves** in real time — [Fraunhofer IOSB](https://iosb.fraunhofer.de/en/competences/optronics/visual-inspection-systems/Sorting-for-bulk-goods-and-recycling.html); [Fraunhofer publica](https://publica.fraunhofer.de/handle/publica/479145). TOMRA OBTAIN applies deep learning for single-particle detection on XRT sorters — [TOMRA](https://cm.tomra.com/en/mining); [IM Mining](https://im-mining.com/?p=101009).

**P&ID / engineering drawing digitization**
- A Microsoft ISE production-style pipeline: CNN symbol detection, OCR, then **Zhang-Suen thinning** followed by the **standard Hough transform** for lines (contour tracing was considered but not chosen). Graph assembly uses NetworkX + Shapely into SQL Graph. It takes about **20 s per 7168×4561 image** and detects about 80% of assets and connections. Pain points include dashed vs solid lines, small symbol/line misalignments (fixed by extending lines with a buffer), and arrow direction heuristics — [MS ISE blog, Feb 2024](https://devblogs.microsoft.com/ise/engineering-document-pid-digitization/).
- Relationformer extracts the symbol+connection **graph** end-to-end. It improves edge detection by more than 25% over the modular baseline and ships the first public graph-annotated P&ID benchmark — [arXiv 2411.13929](https://arxiv.org/abs/2411.13929v3).

**SDF fonts / level-set interfaces**
- msdfgen goes from **vector to SDF** (shapes are closed contours of linear, quadratic and cubic edges, with edge coloring for sharp corners). It is the opposite direction from a polygonizer — [msdfgen](https://github.com/Chlumsky/msdfgen).
- The reverse direction does appear in neural vector fonts. VecFontSDF (CVPR 2023) learns SDF primitives bounded by parabolic curves that convert exactly to quadratic Béziers — [arXiv 2303.12675](https://arxiv.org/html/2303.12675v3). Chen et al. (ICCV 2023) learn SDF plus a probabilistic *corner field* and design a dedicated vectorization step to get compact outlines — [ICCV 2023](https://openaccess.thecvf.com/content/ICCV2023/papers/Chen_Joint_Implicit_Neural_Representation_for_High-fidelity_and_Compact_Vector_Fonts_ICCV_2023_paper.pdf).
- GPU marching squares for floating-point contour edges is patented prior art — [US 12315131](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12315131). GPU dual contouring and marching cubes variants exist — [Eurographics DL](https://diglib.eg.org/items/c94e2e75-fc6a-4c55-8e77-c6ef52460288). Sparse GPU level sets are covered in [TU/e](https://research.tue.nl/en/publications/fast-sparse-level-sets-on-graphics-hardware/).

### Inferences
- **Mining fragmentation is the strongest new topo/reward fit.** The business metric (fragment count and size distribution) is a direct function of connected-component topology. A GPU engine that emits per-fragment rings and counts at 15+ FPS could serve as an eval metric and as a training reward (e.g., penalize merges and splits against ground-truth counts).
- **Oil spill**: the product already *is* polygons with area attributes, so the engine's value would be batch-θ over confidence thresholds (spill area vs threshold curves) and cutting the 20–30 min delivery. A low IoU (~0.57) means uncertainty-aware outlines matter more than speed.
- **P&ID** is skeleton/graph extraction rather than area polygons. The engine fits only if it also offers thin-line/centerline tracing. It is a weak fit for rings+holes, though symbol region masks could use it.
- **Wind blades and power lines** are weak fits today, because the outputs are boxes and point clouds. The exceptions are crown polygons and satellite "encroachment zone" boundaries.
- Neural SDF fonts are a niche case where zero-level-set → ring with exact corner handling is the actual bottleneck step. They are small scale (glyphs), so the value is in training-loop differentiability or rewards, not throughput.

### Gaps
- No RL usage was found in any of these subdomains. The only ML training loop where geometric/topological losses plausibly plug in is fragmentation, and that is my inference.
- No latency budget was found for TOMRA ejection, and no confirmation that sorters ever form polygons. They likely act on per-pixel or per-blob masks.
- Vendor formats for wind-blade defects (e.g., Aerones, SkySpecs) and AiDash/Overstory outputs were not found.
- No level-set *multiphase CFD* interface-extraction performance numbers were found.

---

## Q2. Geospatial items

### Takeaway
National- and global-scale field-boundary polygonization is now real and quantified. Examples are **3.17 B field polygons globally** and **Ukraine (603k km²) in 5.4 h on one consumer GPU desktop**. Pipelines still do overlap merge → raster mosaic → CPU polygonize → topology fix → custom simplification, which is exactly the engine's scope (seams+canon+topology). PostGIS, WhiteboxTools, tippecanoe and Martin are either CPU raster→vector with no published benchmarks, or vector-only consumers. Varuna had no findable technical presence.

### Cited Findings

**Field boundaries (EU IACS/LPIS, global, Varuna)**
- A global 10 m field-boundary map has **3.17 B polygons** (1.62 B in 2024, 1.55 B in 2025) across 241 countries. It comes from a U-Net trained on Fields of The World and applied to Sentinel-2 mosaics, with mean pixel recall 0.85, F1 of 0.89/0.88/0.74 for Austria/Latvia/Finland, and a 500 m confidence layer. The abstract gives no polygonization method or cost — [arXiv 2605.11055](https://arxiv.org/abs/2605.11055).
- **DelAnyFlow** (YOLOv11-based instance segmentation trained on FBIS-22M: 672,909 images and 22.9 M instances at 0.25–10 m) — [arXiv 2511.13417](https://arxiv.org/pdf/2511.13417). Post-processing steps:
  1. per-instance erosion → keep the largest component → dilation;
  2. quality and min-area filtering;
  3. **cross-tile merge of overlapping instances by IoU/area thresholds**;
  4. mosaic into a single raster with unique IDs;
  5. **raster→vector**, topological validation and **sliver removal**;
  6. simplification by rasterizing edges, finding "anchor" vertices shared between adjacent fields, then running RDP between anchors (shared-edge-preserving simplification).
- DelAnyFlow throughput and counts: Ukraine was produced in **5.4 h** on a Ryzen 9 9900X + **RTX 5070 Ti 16 GB** + 64 GB RAM (5 m data). It found 3.75 M fields at 5 m and 5.15 M at 2.5 m, vs Sinergise's 2.66 M. The authors argue *instance IoU* rather than boundary IoU is the operationally aligned metric — same source.
- EU context: AI4Boundaries labels come from member-state GSAA/LPIS data (Austria, Catalonia, France, Luxembourg, Netherlands, Slovenia, Sweden) — [FAO AGRIS](https://agris.fao.org/search/fr/records/68879c937fd4d06c32a9879c). Mask R-CNN field delineation was produced under the H2020 NIVA "New IACS Vision in Action" project — [GitHub field-delineation](https://github.com/ManuelHuber-Github/field-delineation). Planet/Sinergise ship field delineation as an Area Monitoring "marker" — [Planet docs](https://docs.planet.com/data/area-monitoring/markers/field-delineation/).
- **Varuna**: a search for its field-boundary tech returned nothing about the company. One result was a USPTO patent on field-boundary detection from coarse-resolution satellite data plus GIS. Its assignee was not checked, so it should not be attributed to Varuna — [USPTO 12340516](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12340516).

**PostGIS ST_DumpAsPolygons**
- PostGIS 3.3.0 **disabled validation and fixing in ST_DumpAsPolygons for performance**, and the docs warn the output may be invalid. Version 3.7.0 made polygonization honor interrupts and statement timeouts, which implies long-running calls are a real occurrence — [PostGIS docs](https://postgis.net/docs/RT_ST_DumpAsPolygons.html); [3.7 docs](https://postgis.net/docs/manual-3.7/es/RT_ST_DumpAsPolygons.html).
- An older ticket records ST_MakeValid making dump-as-polygons much slower — [trac #1586](https://trac.osgeo.org/postgis/ticket/1586) (via snippet).

**WhiteboxTools**
- `raster_to_vector_polygons` exists (R/Python wrappers), but no published performance data was found — [R docs](https://search.r-project.org/CRAN/refmans/whitebox/html/wbt_raster_to_vector_polygons.html). FME users have an open idea requesting better raster→polygon support (RasterToPolygonCoercer) — [Safe community](https://community.safe.com/ideas/improve-support-for-raster-to-vector-rastertopolygoncoercer-29818) (content not fetched).

**tippecanoe / Martin (downstream consumers)**
- tippecanoe ingests GeoJSON/FlatGeobuf and writes PMTiles natively since v2.17. Protomaps recommends ogr2ogr→PMTiles only for small data because tippecanoe builds better overviews — [Protomaps docs](https://docs.protomaps.com/pmtiles/create). The typical raster workflow is `gdal_polygonize` (one polygon per connected same-value region) → tippecanoe — [GDAL](https://gdal.org/en/stable/programs/gdal_polygonize.html).
- Martin serves MVT on the fly from PostGIS. Each ST_AsMVT query gathers rows, **simplifies to tile resolution, clips to the tile (wagyu), and encodes**. PostGIS 3 parallelized ST_AsMVT and integrated wagyu, which reduced invalid clipped output — [Crunchy Data](https://info.crunchydata.com/blog/waiting-for-postgis-3-st_asmvt-performance); [Martin](https://github.com/maplibre/martin).

**Carbon MRV (forest carbon project boundaries)**
- An open database holds boundaries of **575 nature-based forest carbon projects in 55 countries** from six registries including Verra, obtained from registries, project managers and **manual georeferencing** — [PMC11972284](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11972284/) (via search snippet; full text blocked by CAPTCHA). Distributed on [source.coop/cecil](https://source.coop/cecil/forest-carbon-boundaries).
- Verra REDD+ GIS deliverables include a "single clean polygon that is topologically correct" project boundary, forest stratification and leakage belts — [Workana job post](https://www.workana.com/job/especialista-em-sig-gis-para-projeto-de-carbono-redd-padrao-verra) (weak source; practitioner job ad). Verra has issued RFPs for jurisdictional activity data and forest-cover benchmark maps — [Verra RFP](https://verra.org/documents/request-for-proposals-development-of-jurisdictional-activity-data-and-forest-cover-benchmark-maps-for-vcs-avoiding-unplanned-deforestation-projects).

**Parametric insurance**
- Flood parametric covers use satellite-derived or modeled **flood footprints** and pay on predefined thresholds — [Swiss Re](https://corporatesolutions.swissre.com/alternative-risk-transfer/parametric-solutions/flood.html). Swiss Re tested SAR flood data from ICEYE — [Beinsure](https://beinsure.com/news/swiss-re-iceye-tests-parametric-flood-insurance/). An NYC community pilot uses satellite, sensors and social media to decide trigger breaches — [Floodbase blog](https://www.floodbase.com/blog-posts) (via snippet).

### Inferences
- **DelAnyFlow is the clearest "exact problem statement" found.** It needs a cross-tile instance merge, then mosaic → polygonize → validate → remove slivers → shared-edge simplification. The engine's seams + canonical encoding + topology-preserving shared-edge output would replace steps 3–6. The polygonization share of the 5.4 h is unknown.
- With a 3.17 B-polygon global product on annual refresh, batch-θ (multiple confidence cuts) and canonical encoding (deterministic diffing between years) have concrete value.
- PostGIS trading validity for speed in 3.3 is documentary evidence that a fast *and valid-by-construction* polygonizer is wanted in-database (e.g., a GPU-backed raster extension).
- Carbon MRV and parametric insurance are mainly **polygon-in-raster consumers**: zonal stats of forest-loss or flood rasters inside fixed boundaries. Polygonization matters for *flood footprints* and *deforestation-event polygons*, not for project boundaries. The ML loop is forest/flood segmentation, where topology/area rewards (area agreement inside the boundary) could serve as an evaluation metric.
- RL angle: none found. Field delineation trains with instance-mask losses, and an instance-IoU or count-based geometric reward is a plausible addition (inference).

### Gaps
- No runtime split (inference vs vectorization) for DelAnyFlow or the 3.17 B global map.
- Varuna's pipeline and EU paying-agency (IACS) operational ML systems were not found. NIVA outputs were not fetched.
- No WhiteboxTools or ST_DumpAsPolygons throughput benchmarks exist publicly as far as found.
- No volume numbers for parametric triggers or carbon-MRV polygon processing.

---

## Q3. Biomedical items

### Takeaway
Real-time smart microscopy has sub-200 ms segmentation budgets and *does* use RL and bandits for acquisition decisions. That makes it the one domain in this pass with a live RL loop, where shape and topology signals from masks could feed rewards or state. Connectomics meshing (3D) has hard cost numbers. It is out of 2D scope but is the volumetric analogue of seams via a one-voxel overlap.

### Cited Findings
- **EAP4EMSIG** (event-driven microfluidic live-cell imaging, 2025) evaluated 11 DL segmentation methods. Cellpose reached PQ 93.36%, and the fastest (distance-based) took **121 ms** at PQ 93.02%. MLP autofocus inference starts at 87 ms. The authors stress that delayed reactions can irreversibly ruin experiments — [arXiv 2504.00047](https://arxiv.org/html/2504.00047v2).
- **RL and bandits in microscopy:**
  - Restless multi-process multi-armed bandits schedule which of thousands of ROIs to image in HCS, under photobleaching and compute budgets — [arXiv 2512.14930](https://arxiv.org/pdf/2512.14930).
  - An RL agent controls an LED-array illumination pattern with a reward that balances classification confidence against acquisition cost — [Duke](https://dukespace.lib.duke.edu/items/682d6b59-5ba6-426e-9763-23c452578fba).
  - Deep-RL autofocus also exists — [Microsystems & Nanoengineering 2024](https://doaj.org/article/7e38985063be490993303e4532dd76da).
- **Ultrack** (ECCV 2024, CZ Biohub) builds a *hierarchy of segmentation hypotheses* (ultrametric contour maps) per frame and selects disjoint segments by maximizing overlap across frames via ILP. It scales to **millions of instances in TB-scale 3D+t** data — [arXiv 2308.04526](https://arxiv.org/abs/2308.04526v1); [ECCV](https://www.ecva.net/papers/eccv_2024/papers_ECCV/papers/10114.pdf).
- **Connectomics meshing (Igneous/zmesh)** — [Frontiers Neural Circuits 2022](https://www.frontiersin.org/journals/neural-circuits/articles/10.3389/fncir.2022.977700/full):
  - zmesh does multi-label marching cubes plus simplification;
  - S1 throughput was ~**18.2 gigavoxels/core-hour**;
  - the MICrONS-scale job processed **74 teravoxels in 99,800 core-hours** plus 2,450 core-hours for merging;
  - chunks use a **one-voxel overlap** so fragments stitch trivially;
  - simplification is topology-guarded Garland-Heckbert, 100× at 40 nm max error;
  - pain point: hundreds of millions of fragment files, hence sharded formats.

### Inferences
- **Smart microscopy** is the best RL hook found in this pass. The agents (bandits, illumination RL) consume per-ROI state, and per-cell polygon features computed on GPU in under 10 ms would fit the 100–200 ms loop. Examples are area, perimeter, hole count and component count, which flag division events (1→2 components) or vacuoles (holes). Rewards could include "topology event detected" signals.
- Ultrack's multi-hypothesis hierarchy is effectively **batched thresholds over a contour map**. A batch-θ ring extractor with topology counts per threshold maps directly onto generating its candidate segments, in 2D per slice at least.
- Connectomics shows the seam technique (1-voxel overlap) and cost scale for 3D. A 2D engine does not address it directly.

### Gaps
- No HCS vendor (Harmony, CellProfiler) polygon export or latency figures were found. No organoid contour-tracking specifics.
- No RL paper was found that uses *shape/topology* rewards in microscopy. Existing rewards are classification confidence or acquisition cost.
- FlyWire/CAVE remesh latency was not found.

---

## Q4. Robotics items

### Takeaway
Open-vocabulary scene graphs keep objects as **point clouds**. Room footprints come out as **watershed 2D masks, not polygons**, and the room step is called too slow for real-time. That is a direct mask→polygon gap. Search-and-rescue planners (POMCP/RL) operate on probability grids with area partitioning, where searched-area polygons are a plausible but undocumented need.

### Cited Findings
- **ConceptGraphs** associates class-agnostic 2D masks across RGB-D frames into object nodes. Each node is a **3D point cloud + semantic feature**, and LLM/LVLM add captions and edges, with no training — [arXiv 2309.16650](https://arxiv.org/pdf/2309.16650).
- **HOV-SG** (RSS 2024) finds floors from height-histogram peaks. Rooms are derived as follows: BEV projection → thresholded **wall skeleton mask** → dilate → Euclidean distance field → seeds → **Watershed → 2D region masks** (no polygon output). Object nodes are point clouds. The representation is 75% smaller than dense open-vocabulary maps, and the authors call the method "time-consuming, rendering it unsuitable for real-time mapping" — [arXiv 2403.17846](https://arxiv.org/html/2403.17846v2).
- **Search-and-rescue:**
  - Shrinking POMCP plans UAV search over belief/probability maps in real time — [arXiv 2411.12967](https://arxiv.org/pdf/2411.12967).
  - Multi-UAV wilderness SAR uses agent-based probability models with area partitioning (Voronoi) across UAVs — [arXiv 2411.10148](https://arxiv.org/pdf/2411.10148).
  - Land-coverage-aware multi-UAV path planning uses segmentation of the scene for coverage planning — [arXiv 2505.08060](https://arxiv.org/html/2505.08060v1).

### Inferences
- A GPU ring builder could turn HOV-SG/ConceptGraphs BEV room masks into **room polygons with holes** (pillars, islands) every frame, along with canonical IDs for graph diffing. The room layer is currently offline.
- In SAR, thresholding the belief map at several probability levels (batch-θ) gives nested "probable area" polygons for multi-UAV partitioning and for reporting the searched area to incident command. The planners are POMCP/RL, so polygon-derived quantities (uncovered area, number of disconnected unsearched pockets) are candidate shaped rewards. I found no paper doing this, so it is an inference.

### Gaps
- No runtime numbers for HOV-SG's room segmentation step, and no ConceptGraphs footprint-extraction variant.
- No SAR incident-management software (e.g., SARTopo/CalTopo) documentation on auto-generated searched-area polygons was found.
- No ROS package doing GPU mask→room-polygon was found (not specifically searched this round).

---

## Q5. Cross-item ranking for an RL/ML-reward angle (synthesis)

### Takeaway
Ranked by "an ML or RL loop already exists, and topology or ring output directly affects its objective": (1) smart microscopy RL/bandits, (2) rock fragmentation counts and size distribution, (3) national field-boundary instance pipelines (DelAnyFlow), (4) SAR belief-map planners, (5) neural SDF font vectorization. The rest are mainly throughput/plumbing wins (oil spill, PostGIS, IACS) or poor fits (wind blades, power-line LiDAR, tippecanoe/Martin).

### Cited Findings
- Live RL/bandit acquisition loops with ~100–200 ms segmentation budgets — [arXiv 2512.14930](https://arxiv.org/pdf/2512.14930); [arXiv 2504.00047](https://arxiv.org/html/2504.00047v2).
- Merge/split errors directly change fragmentation metrics — [RIT](https://www.rit.edu/dirs/research/structure-aware-segmentation-post-blast-rock-fragmentation-analysis).
- Field pipelines need cross-tile merge, validation, slivers and shared-edge simplification, and argue for instance-level metrics — [arXiv 2511.13417](https://arxiv.org/pdf/2511.13417).
- Belief-map POMCP planners for UAV SAR — [arXiv 2411.12967](https://arxiv.org/pdf/2411.12967).

### Inferences
- The common reward primitives are per-instance component count, hole count, area, perimeter and shared-edge consistency. They are cheap on GPU once rings exist, and they line up with the metrics practitioners already state (fragment count/PSD, instance IoU, division events).

### Gaps
- No source anywhere in this pass uses an explicit *geometric/topological RL reward* derived from polygonized masks. That remains an open opportunity claim, not a documented practice.
