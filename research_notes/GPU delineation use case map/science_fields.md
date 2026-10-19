# Science & Simulation Fields: Thresholding Gridded Fields into Tracked Objects, Regions and Polygons (2023–2026)

Scope: where gridded fields from physics simulations, ML surrogates/emulators and sensors get thresholded, turned into connected regions or contours, and tracked. This is the target market for a CUDA engine that turns binary bitmap → closed polygon rings (scanline, 2x2 classification, halo tiling). Research was done in Oct 2026 with about 30 search/fetch calls. Wherever a number was not found, it is listed under Gaps rather than guessed.

---

## 1. Weather & climate: ML weather models → extreme-event objects (TCs, ARs, heatwaves, blobs)

### Takeaway
ML weather models (GenCast, AIFS, FourCastNet/SFNO, Pangu) now make ensembles of thousands of members cheaply. The feature-extraction step that follows them is still mostly CPU-bound, file-I/O-bound C++/Python (TempestExtremes, tobac, py-eddy-tracker). NVIDIA's Earth2Studio has started moving trackers into PyTorch on the GPU, using threshold → connected-component labeling → erosion. Nobody publishes a GPU-native "threshold → polygon → track" stage that scales to 1000+ ensemble members. That is the clearest gap for this engine.

### Cited Findings
- **TempestExtremes** is the de facto C++ CLI toolkit for detecting and tracking TCs, extratropical cyclones, monsoon depressions, atmospheric blocks, ARs and MCSs in large climate datasets — [BEAR apps / TempestExtremes](https://bear-apps.bham.ac.uk/applications/TempestExtremes/). It has two families:
  - **DetectNodes/StitchNodes** for point features (e.g. MSLP minima with a closed-contour criterion)
  - **DetectBlobs/StitchBlobs** for areal features, tracked through time with an overlap criterion — [UC Davis TempestExtremes page](https://climate.ucdavis.edu/tempestextremes.php); [TempestExtremes v2.1 GMD preprint](https://gmd.copernicus.org/preprints/gmd-2020-303/gmd-2020-303-manuscript-version4.pdf)
- **TempestExtremes performance:** DetectBlobs over ERA5 at 6-hourly resolution from Jan 1979 to Feb 2020 (494 monthly files) took **34 min 42 s**. **66% of runtime was file input** and ~13% was the Laplacian operator. It is parallelised over files with MPI and has no GPU path — [TempestExtremes v2.1 manuscript](https://gmd.copernicus.org/preprints/gmd-2020-303/gmd-2020-303-manuscript-version4.pdf)
- **ML models still use TempestExtremes for TC tracking:**
  - **TCBench** (2026) benchmarks TIGGE NWP plus AIFS, Pangu-Weather, FourCastNet v2 and GenCast. Where tracks are not supplied, it derives them from model output with TempestExtremes — [TCBench arXiv 2601.23268](https://arxiv.org/html/2601.23268v1)
  - The original GenCast paper used TempestExtremes v2.1 — [same](https://arxiv.org/html/2601.23268v1)
  - DetectNodes finds MSL local minima, then requires a "sufficiently compact closed low-pressure structure", i.e. a closed-contour test — [same](https://arxiv.org/html/2601.23268v1)
- **Huge Ensembles (HENS)** (Mahesh et al., LBNL/NVIDIA, SFNO): **7,424 members per initialization day**, initialized every day from 1 Jun to 31 Aug 2023. It samples events 4σ from the ensemble mean, and was built to study TCs, heatwaves and ARs — [HENS Part II arXiv 2408.01581](https://arxiv.org/abs/2408.01581); [HENS Part I arXiv 2408.03100](https://arxiv.org/pdf/2408.03100)
- Later HENS work: surface temperature extremes from the summer-2023 huge hindcasts — [arXiv 2604.09754](https://arxiv.org/pdf/2604.09754); anomalous diffusion of TCs in huge ensembles — [arXiv 2607.21954](https://arxiv.org/pdf/2607.21954)
- AI models have produced **2000-member TC track forecasts** for the first time — [EGU25-2160](https://meetingorganizer.copernicus.org/EGU25/EGU25-2160.html) (search snippet; the full abstract was not read)
- **NVIDIA Earth2Studio** has "huge ensemble" recipes that combine multiple checkpoints, distributed inference, and post-processing with cyclone trackers — [ECMWF HPC2025 Weissenberger talk](https://events.ecmwf.int/event/460/contributions/5523/attachments/3182/5317/HPC2025_Weissenberger.pdf) (search-snippet level; the PDF did not parse). It ships TC tracker "diagnostic models":
  - **TCTrackerVitart** (Vitart 1997 criteria)
  - **TCTrackerWuDuan**, which thresholds vorticity into a **binary image**, runs **connected component labeling**, then erodes to get inner-core TC seeds. It is a `torch.nn.Module` (GPU-capable). It links paths with a 300 km search radius and outputs a `[batch, path_id, steps, variable]` tensor using NaN padding — [Earth2Studio TCTrackerWuDuan docs](https://nvidia.github.io/earth2studio/main/modules/generated/models/dx/TCTrackerWuDuan/); [Cyclone tracking example](https://nvidia.github.io/earth2studio/main/examples/02_medium_range/05_cyclone_tracking/)
  - Earth2Studio also lists **WeatherNext2Cyclones** prognostic models, i.e. direct ML TC-track prediction — [Earth2Studio docs](https://nvidia.github.io/earth2studio/main/modules/generated/models/px/WeatherNext2Cyclones/)
- **ClimateNet:** an expert-labeled dataset of TCs and ARs in CAM5.1 high-res output, used to train a DL **pixel-level segmentation** model (TC/AR/background). The masks then condition precipitation statistics — [ClimateNet GMD 2021](https://gmd.copernicus.org/articles/14/107/2021/). The output is a mask, not polygons.
- **tobac** (Tracking and Object-Based Analysis of Clouds) is the community Python package with three stages:
  1. feature detection (multi-threshold)
  2. **watershed segmentation**, which flows out from feature markers above a threshold, independently per timestep, and can split adjacent objects
  3. linking

  It is used on OLR, updraft and condensate fields from convection-permitting models — [tobac GMD 2024 (v1.5)](https://gmd.copernicus.org/articles/17/6035/2024/gmd-17-6035-2024.html); [tobac docs](https://tobac.readthedocs.io/en/latest/examples/Example_OLR_Tracking_model/Example_OLR_Tracking_model.html). No GPU acceleration is documented.
- **ByteStorm** (Dec 2025) is a multi-step data-driven TC detection and tracking approach — [arXiv 2512.07885](https://arxiv.org/pdf/2512.07885) (title only). Adversarial attacks on downstream TC trajectory models show that the trackers are themselves part of the ML pipeline — [arXiv 2510.10140](https://arxiv.org/pdf/2510.10140) (title only).

### Inferences
- **Ensemble load estimate:**
  - 7,424 members × 92 init days × a multi-day rollout (e.g. 60 steps at 6 h) is roughly 40M global 0.25° fields (721×1440 ≈ 1.04 Mpx each) per variable per threshold.
  - That is about 4×10^13 pixels to classify per variable/threshold.
  - At a GPU scan rate of ~10 Gpx/s (comparable to the ECC kernel in §8), the engine-equivalent work is ~1 GPU-hour per variable/threshold. The CPU/file-bound TempestExtremes path would be orders of magnitude slower.
  - This is an estimate, not a sourced figure.
- **In-line placement:** the binding constraint is I/O (66% in TempestExtremes). The winning design runs **in-line with ML inference on the GPU**: threshold the tensor in device memory, emit polygons and object attributes, and never write full fields. This matches the in-situ argument in §3.
- **WuDuan pattern:** "Binary image → CCL → erosion → seed" is exactly the pattern where a fast bitmap→polygon engine fits. Polygons also give area, perimeter, centroid and shape metrics (e.g. compactness tests like TempestExtremes' closed-contour criterion) for free.
- **Probabilistic polygons:** probability-of-exceedance polygons from N members need either per-member polygons with overlap/union (a polygon boolean problem) or thresholding the ensemble-mean probability field (one bitmap per probability level). The second is cheap and maps directly onto the engine. Doing the first at N=1000+ is a gap.

### Gaps
- No published wall-clock or GPU-hours figure was found for the TC/heatwave extraction step of HENS. The "25 PB" figure often quoted for HENS could not be verified in what was read.
- No GPU port of TempestExtremes, TECA or tobac was found. TECA (LBNL Toolkit for Extreme Climate Analysis) did not appear in 2023–2026 results, so its current status is unknown.
- No source found for how Aurora, GraphCast or AIFS operational products extract objects such as heatwave regions or fronts. Front-detection tooling was not found in this pass.
- No numbers on the computational cost of turning 50–1000 members into probabilistic polygons.

---

## 2. Operational nowcasting & warning polygons (NWS, FACETs/PHI, ProbSevere, Europe)

### Takeaway
Operational severe-weather systems already emit **polygons every ~2 minutes**. Examples are ProbSevere storm-object GeoJSON, DWD NowCastMIX warning polygons, and FACETs PHI objects on ≤2 km grids. They are built from watershed/threshold segmentation of radar mosaics on CPU (WDSS-II). Observed end-to-end latency is 2–5 minutes, which is a concrete pain point.

### Cited Findings
- **ProbSevere (NOAA/CIMSS)** identifies storms in **MRMS** merged reflectivity using **WDSS-II's enhanced watershed** algorithm and tracks them by matching objects between consecutive images — [CIMSS ProbSevere training](https://cimss.ssec.wisc.edu/severe_conv/training/slides/radar_data.html); [WAF paper](https://cimss.ssec.wisc.edu/severe_conv/training/waf-d-17-0099.1.pdf)
  - About **every 2 minutes** it outputs a GeoJSON feature collection. Each feature carries a **polygon** plus ProbSevere/ProbTor/ProbHail/ProbWind and environmental attributes — [NSSL EWP testbed blog](https://inside.nssl.noaa.gov/ewp/?p=19142)
  - Forecasters in the HWT testbed observed **2–5 min latency (~2 radar scans)** — [NSSL EWP blog](https://inside.nssl.noaa.gov/ewp/tag/glm/page/4)
- **FACETs** (Forecasting a Continuum of Environmental Threats) replaces deterministic warning polygons with rapidly updating **Probabilistic Hazard Information (PHI)** on **≤2 km grids** — [NSSL FACETs](https://www.nssl.noaa.gov:443/projects/facets/FACETs%20SSIP_1.pdf); [Eos project update](https://eos.org/project-updates/next-generation-severe-weather-forecasting-communication)
  - **Warn-on-Forecast** ensembles are the source of the probabilistic storm-scale output, and FACETs is the "delivery mechanism" — [same](https://eos.org/project-updates/next-generation-severe-weather-forecasting-communication)
  - NSSL's PHI prototype tracks storms as **geospatial "objects"** — [same](https://eos.org/project-updates/next-generation-severe-weather-forecasting-communication)
  - Threats-in-Motion (TIM) makes warning polygons move continuously with the storm — [CIRA TIM](https://www.cira.colostate.edu/?p=14346)
- **Europe:**
  - **DWD NowCastMIX** produces forecast warning polygons smoothed in space and time. Forecasters can edit them, and they can be projected onto administrative borders — [GMD 2025 cell-tracking nowcast evaluation](https://gmd.copernicus.org/articles/18/1851/2025/)
  - A 2025 GMD framework evaluates nowcasts (advection, S-PROG, LINDA, **L-CNN**) by **cell tracking** on Swiss radar. L-CNN reproduced convective-cell existence best — [GMD 18, 1851 (2025)](https://gmd.copernicus.org/articles/18/1851/2025/)
  - OPERA-composite cell tracking plus random-forest severity models, trained on 2018–2022 reports, are also in use — [same results](https://gmd.copernicus.org/articles/18/1851/2025/)

### Inferences
- **Latency budget:** radar volume scans arrive every ~2–5 min, and MRMS mosaics are ~1 km CONUS (~7000×3500 px; known grid size, not sourced here). Polygonization itself is milliseconds on GPU. The 2–5 min latency is dominated by the pipeline (mosaicking, CPU segmentation, I/O), so a GPU engine helps mainly by being fused in-line.
- **ML nowcasts:** models like MetNet and L-CNN output probability fields. Turning them into PHI-style polygons at several probability levels × lead times × hazards is a natural batched-threshold workload.
- **Smoothing and simplification:** operational polygons must be smoothed and simplified (NowCastMIX smooths in space and time, NWS polygons have vertex limits). A raw pixel-edge ring output needs a simplification and smoothing stage after it.

### Gaps
- No source on Google MetNet-3 or WeatherNext products emitting polygons. MetNet outputs appear to be gridded probabilities only, which was not verified.
- No public numbers on ProbSevere's per-cycle compute time or on PHI grid-to-object conversion cost.

---

## 3. Ocean: eddies, sea ice floes (marine heatwaves / oil spills not covered)

### Takeaway
Ocean eddy detection is literally **closed-contour extraction on SSH at 2 mm steps**, run daily on a global 1/4° grid (~6000 eddies/day) in Python (py-eddy-tracker). It is a textbook "many-threshold contour → shape test → track" workload with no GPU path. Sea-ice floe segmentation is moving to foundation models (SAM) plus instance segmentation → floe polygons → floe-size distributions.

### Cited Findings
- **py-eddy-tracker (PET)** powers AVISO+ **META3.1exp** (daily global Mesoscale Eddy Trajectory Atlas) — [ESSD 2022, Pegliasco et al.](https://essd.copernicus.org/articles/14/1087/2022/essd-14-1087-2022.html). Detection details:
  - Scans **closed SSH contours from maxima downward (anticyclones) / minima upward (cyclones) with a 2 mm step** (older methods used 1 cm)
  - **Shape-error test:** deviation from the best-fit circle must be below **70%**
  - **Minimum amplitude 0.4 cm** and **≥5 pixels**
  - **1/4° grid**
  - **~5,936 eddies/day globally**, about **59.4 million** detections over 1993 onward — [same](https://essd.copernicus.org/articles/14/1087/2022/essd-14-1087-2022.html)
- PET is in Python (GPL v3) and is used on global daily SLA/ADT with a 500 km high-pass Bessel filter — [ESSD 2025 vol. 17](https://essd.copernicus.org/articles/17/7119/2025/essd-17-7119-2025.html). PET also has an **AreaTracker** for tracking by contour overlap — [PET docs](https://py-eddy-tracker.readthedocs.io/en/stable/_autosummary/py_eddy_tracker.featured_tracking.area_tracker.AreaTracker.html)
- **Eddy work on ML ocean models:** eddy demographics in eddy-permitting models vs reanalysis — [arXiv 2604.22524](https://arxiv.org/pdf/2604.22524); ML global ocean forecasting on unstructured meshes (OceanLight) — [arXiv 2608.16070](https://arxiv.org/pdf/2608.16070); a deep-learning eddy dataset (PhyEddy-Net 2010–2021) — [Mendeley Data](https://data.mendeley.com/datasets/rh6wtx3hyr). These were title/snippet level only.
- **Sea ice floes:**
  - Floe-size distribution (FSD) needs instance segmentation of individual floes from airborne or satellite imagery
  - DL semantic segmentation has been applied to Sentinel-2 in the marginal ice zone — [The Cryosphere 17, 5519 (2023)](https://tc.copernicus.org/articles/17/5519/2023/)
  - GVF active contours vs **SAM** vs hybrid: prompt-driven SAM gave the best accuracy/efficiency trade-off — [arXiv 2409.06641](https://arxiv.org/html/2409.06641v3)
  - Seasonal FSD from Sentinel-1 in Nares Strait — [EGU25-1736](https://meetingorganizer.copernicus.org/EGU25/EGU25-1736.html)

### Inferences
- **Fit with PET:** PET's core loop (hundreds of contour levels per day × closed-contour test × shape-error vs circle × amplitude) is a strong fit for a batched GPU engine. Thresholding at each 2 mm level produces a bitmap per level, and closed rings come straight out. The nesting of rings across levels needs a containment/merge-tree step, i.e. the same structure as a merge tree (§4).
- **Scale:** a 1/4° global grid is small (~1M px), but ML ocean ensembles and multi-decade reanalyses (~12k days × many levels × many members) multiply it quickly.
- **Floe polygons:** floe FSD needs polygons with holes and touching floes separated. A bitmap→ring engine needs instance-aware input (label maps, or edges from SAM).

### Gaps
- No PET runtime numbers found. Marine heatwave region extraction and oil-spill polygon pipelines were not researched in this pass, so no sourced data is available.

---

## 4. Turbulence / CFD / combustion: in-situ feature extraction, topology (TTK, merge trees), data reduction

### Takeaway
Exascale workflows (ECP ALPINE: Ascent, ParaView Catalyst, VisIt) push feature extraction in situ to avoid I/O. **Segmented merge trees** encode all threshold-based features (flame kernels, ignition regions) compactly. TTK provides merge, contour and Morse-Smale trees, mostly on CPU; GPU merge-tree work exists (SYCL, VTK-m), but it is research code. Topology-preserving compression (TopoSZ, 2024; a general framework, 2025) is the current answer to "store less than the full field".

### Cited Findings
- **ECP ALPINE:**
  - Made ParaView and VisIt exascale-ready and built the lightweight in-situ infrastructure **Ascent**
  - Its algorithms target automated feature focus and **data reduction to meet exascale I/O limits**
  - [ALPINE paper (NSF PAR)](https://par.nsf.gov/biblio/10548487-ecp-alpine-project-situ-post-hoc-visualization-infrastructure-analysis-capabilities-exascale); [ECP DAV BOF slides](https://www.exascaleproject.org/wp-content/uploads/2021/01/talk10-2.3.4.16-ECP-DAV-BOF-2021-ALPINE.pdf)
- **In-situ feature extraction for combustion via segmented merge trees** (Sandia/Utah, Bremer/Bennett et al.) encodes "a wide range of threshold-based features". It is the canonical pattern for flame and ignition-kernel tracking without storing fields — [Sandia publication](https://sandia.gov/research/publications/details/in-situ-feature-extraction-of-large-scale-combustion-simulations-using-segm-2014-01-16) (2014, foundational, not new)
- ECP also published feature extraction and visualization for **multiphase flow** (droplet/bubble tracking) — [ECP publication](https://www.exascaleproject.org/publication/feature-extraction-and-visualization-algorithm-improves-functional-memory-and-research-outcomes-in-multiphase-flow-simulation-analysis/)
- **TTK (Topology ToolKit)** implements critical points, persistence diagrams, **merge trees, contour trees, Morse-Smale complexes**, Reeb spaces and more, integrated with ParaView — [TTK arXiv 1805.09110](https://arxiv.org/pdf/1805.09110)
- 2024 merge-tree advances:
  - **ExTreeM** (scalable augmented merge trees via extremum graphs, IEEE TVCG 2024) — [Will lab page](https://vis.cs.rptu.de/team/will/)
  - Merge tree geodesics/barycenters with path mappings (TVCG 2024) — [same](https://vis.cs.rptu.de/team/will/)
- **GPU merge trees:** the triplet merge-tree algorithm has been ported to **SYCL** for GPUs and compared against **VTK-m's** implementation. The motivation is that analysis should run where the simulation runs, i.e. on GPUs — [arXiv 2301.10838](https://arxiv.org/pdf/2301.10838)
- **Topology-preserving compression:**
  - **TopoSZ** (IEEE TVCG 30(1), Jan 2024) keeps extrema and contour-tree level-set relations under SZ1.4 error-bounded compression — [arXiv 2304.11768](https://arxiv.org/pdf/2304.11768)
  - A **2025 general framework** adds topological guarantees to SZ3, TTHRESH, ZFP and neural compressors (Neurcomp) — [arXiv 2502.14022](https://arxiv.org/html/2502.14022v1)

### Inferences
- **Two data-reduction routes:** "store contours instead of fields" competes with:
  1. topology-preserving compression, which keeps the field
  2. merge trees, which keep all thresholds implicitly
- **Where the engine fits:** a fast bitmap→polygon engine is the right tool when the threshold is known in advance (flame front at T = T_iso, vortex region at Q > 0, droplet interface at VOF = 0.5). It emits explicit boundary geometry per timestep, typically KB–MB versus GB–TB fields.
- **2D vs 3D:** in 3D CFD the analog is isosurfaces or marching cubes, and 2D ring output applies to slices, wall-shear maps and 2D simulations. That limits direct applicability.
- **Merge-tree connection:** a scanline 2×2-classification engine is structurally close to the "sweep" in merge-tree construction (union-find over sorted values). A multi-threshold variant could produce nested polygons that equal the merge-tree segmentation at sampled levels.

### Gaps
- No 2023–2026 numbers found for in-situ contour output volumes vs raw field volumes, or for Ascent/Catalyst feature-extraction overheads.
- TTK GPU support status as of 2025–2026 was not confirmed; TTK appears mostly CPU/OpenMP.
- Vortex region tracking (λ2/Q-criterion objects) and flame-front tracking specifics for 2023–2026 were not found.

---

## 5. Materials science: microstructure segmentation, grain boundary networks

### Takeaway
The materials workflow has moved to DL semantic segmentation of grain boundaries (SEM images trained on EBSD labels). The key pain point is **topology/connectivity of the grain-boundary network**: broken boundaries merge grains. The output that matters downstream is a grain-boundary **polygon network** (grain size distributions, triple junctions), but tooling for robust mask→polygon-network conversion is not prominent in the literature.

### Cited Findings
- **PNNL (Jan 2025):** automated grain-boundary segmentation of 347H stainless steel. SEM images are the input and **EBSD maps serve as pixel-wise labels**. **Adding topological information during training improved boundary-network connectivity** and segmentation quality — [PNNL publication](https://www.pnnl.gov/publications/automated-grain-boundary-gb-segmentation-and-microstructural-analysis-347h-stainless)
- **PNNL (Mar 2026):** generalizable segmentation with minimal ground truth, using correlated SEM+EBSD on friction-stir-processed 316L — [PNNL publication](https://www.pnnl.gov/publications/generalizable-image-segmentation-microstructure-characterization-through-integrated)
- DL phase segmentation in EBSD maps of dual-phase steels — [arXiv 2112.03072](https://arxiv.org/pdf/2112.03072). NMF/PCA-based approaches have been used to segment low-angle grain boundaries — [same results](https://hal.univ-lorraine.fr/hal-03480629v2/document)

### Inferences
- **Topology-aware conversion:** a grain map is a **planar partition**, a label image where every pixel belongs to a grain. Polygonizing it needs **shared edges and topology preservation**, not independent per-object rings.
  - A binary-occupancy engine run per grain duplicates every shared boundary.
  - A label-map (multi-class) mode that emits a planar graph (edges + triple junctions) would be a differentiator.
  - The same need appears in land-cover and cadastral domains (other branches).
- **Phase-field:** phase-field simulations (GPU-native, often 2D) produce order-parameter fields that are thresholded into grains every timestep. Grain tracking across timesteps is a direct fit, but no 2023–2026 tooling source was found.

### Gaps
- No sourced numbers for EBSD map sizes, micro-CT volume sizes or pore-segmentation throughput in 2023–2026. No source on GPU grain-tracking in phase-field codes. No source on "digital twin of materials" polygon outputs.

---

## 6. Astronomy: source extraction (Rubin/LSST Footprints, SKA, radio cubes)

### Takeaway
Astronomy pipelines are **threshold → connected region ("Footprint") → peaks → deblend**. Rubin's pipeline represents detections as Footprints and HeavyFootprints. Radio source finding for the SKA era is described (2026) as among the most computationally demanding problems in astronomy, with finders needing adaptation to petabyte scale. Existing finders (SoFiA-2, PyBDSF, Aegean) are CPU/OpenMP/MPI.

### Cited Findings
- **Rubin/LSST Science Pipelines:**
  - `SourceDetectionTask` "finds above-threshold regions and peaks within them (**Footprints**)" and creates a parent source per Footprint
  - `MultibandDeblendTask` creates a child per peak in multi-peak Footprints, each with a **HeavyFootprint** (pixel region + fractional pixel values)
  - Deblending uses **scarlet / scarlet lite** (gradient-descent multiband models)
  - [LSST pipelines scarlet docs](https://pipelines.lsst.io/v/v27_0_0/modules/lsst.meas.extensions.scarlet/index.html); [DP1 deblender footprints tutorial](https://dp1.lsst.io/tutorials/notebook/206/notebook-206-2.html); [Rubin 2023 PCW detection & deblending slides](https://project.lsst.org/meetings/rubin2023/sites/default/files/PCW-detction_and_deblending.pdf)
- **SKAO source finding (2026 review):** source finding "has evolved into one of the most computationally demanding challenges in modern astronomy". Traditional finders (Aegean, PyBDSF) "must now be adapted to handle petabyte-scale datasets" — [arXiv 2607.03736](https://arxiv.org/pdf/2607.03736)
- **SKA Science Data Challenge 2** evaluated 3D HI source finding methods — [arXiv 2303.07943](https://arxiv.org/pdf/2303.07943)
- **SoFiA 2** (WALLABY/ASKAP HI survey):
  - Multithreaded with OpenMP, plus **SoFiA-X** to split cubes across nodes
  - Processes "**hundreds of GB of imaging data per epoch in real time**"
  - Near-100% completeness/reliability above integrated S/N ~5–6
  - [arXiv 2106.15789](https://arxiv.org/abs/2106.15789)
- **Aegean 2.0** was designed for SKA-precursor source finding — [arXiv 1801.05548](https://arxiv.org/pdf/1801.05548)

### Inferences
- **Footprint representation:** in the LSST stack, a Footprint is internally a **SpanSet**, a run-length list of (y, x0, x1) spans. This is from my knowledge of the `afw` library and was not re-verified in this pass. It is directly analogous to the scanline intervals a scanline polygonizer builds internally. An engine that can emit **spans as well as rings** would plug straight into LSST-style data models.
- **Polygons in astronomy:** they are used for survey masks (bright-star masks, footprint masks) and for "region" files (DS9/MOC). Per-source representation is mostly spans or masks, not polygons.
- **Radio:** SKA image cubes are 3D (RA, Dec, frequency). 2D polygonization applies to continuum images and per-channel moment maps. 3D sources need 3D labeling.

### Gaps
- No Rubin per-visit detection timing or Footprint counts were found. No Euclid pipeline (MER) segmentation numbers. No GPU source-finder benchmarks for SKA (e.g. GPU-based SoFiA) were found.

---

## 7. Particle physics / detectors: GPU connected-component clustering in triggers

### Takeaway
HEP trigger systems are the most mature GPU connected-component deployments in science, running at **30–40 MHz collision rates and Tbit/s–TB/s data rates** on hundreds to thousands of GPUs. Their output is **cluster centroids and charges, not polygons**. They are proof that GPU CCL at extreme throughput is solved for sparse, zero-suppressed hits, which is relevant competition and inspiration rather than a market for polygons.

### Cited Findings
- **LHCb Allen:** a GPU first-level trigger for Run 3. It processes the **40 Tbit/s** detector output at the **30 MHz** collision rate on about **500 GPUs**. VELO clustering finds **8-connected groups of pixels** and stores average (x, y), a variant of CCA — [vom Bruch Allen seminar](https://indico.in2p3.fr/event/20264/attachments/57067/76015/vom_Bruch_Allen_LPNHE_seminar.pdf); [arXiv 1912.09161](https://arxiv.org/abs/1912.09161v1)
- **CMS Patatrack:** heterogeneous pixel local reconstruction plus track and vertex reconstruction. On an **NVIDIA T4** it matches the performance of **52–56 Xeon cores** — [search summary of Patatrack literature](https://www.epj-conferences.org/10.1051/epjconf/202429502023) (number from snippet; not re-verified in the full text)
- **ALICE Run 3:**
  - Records **50 kHz Pb-Pb** in continuous readout
  - Detectors produce **3.5 TB/s**, reduced by **zero suppression 3.3 TB/s → 900 GB/s**, then by GPU TPC clusterization, tracking and compression to **170 GB/s**
  - The **EPN farm has 350 servers × 8 AMD GPUs (32 GB each)**, 2×32-core CPUs and 512 GB RAM per server
  - [arXiv 2411.02206](https://arxiv.org/html/2411.02206v1); [Rohr arXiv 2502.09138](https://arxiv.org/html/2502.09138v1); [CERN Courier](https://cerncourier.com/a/alice-ups-its-game-for-sustainable-computing/)
- **Hybrid pixel detectors** (Timepix-class; Dec 2024):
  - GPU **union-find with path compression** over zero-suppressed data reaches **up to 300 million hits/s**, about **100× faster than CPU**, on an RTX 4070 Ti Super
  - Detector rates can exceed **40M hits/s**
  - Cluster sizes range from 2.46 px (gammas) to 7,258 px (lead ions)
  - Uses 256×256 auxiliary matrices, pinned memory and CUDA streams
  - [arXiv 2412.11809](https://arxiv.org/html/2412.11809v1)

### Inferences
- **What carries over:** HEP CCL works on **sparse lists** (hits), while the user's engine works on **dense bitmaps with halos**. Calorimeter clustering (e.g. HGCAL CLUE-style density clustering) and dense pixel frames are closer.
- **Throughput benchmark:** 300M hits/s (sparse) and 5–10 Gvox/s (dense ECC, §8) are useful reference points for claims about the engine's own speed.
- **Polygons in HEP:** few, apart from event display and possibly HGCAL cluster boundary shapes. HEP is mainly a validation domain, not a buyer.

### Gaps
- No precise Patatrack pixel-clustering kernel latency (µs/event) was found in 2023–2026 sources, and no ALICE clusterizer-only throughput figure.

---

## 8. Topological summaries: Euler characteristic curves and Minkowski functionals on GPU

### Takeaway
Euler characteristic curves (ECC) and Minkowski functionals are computed with **local voxel/pixel-configuration lookups**, the same 2×2 (bit-quad) neighborhood logic as the user's engine. GPU implementations reach **5–10 Gvox/s**, and recent tensor/vectorized ECT methods (pyECT, SoCG 2026) are GPU-native. This domain consumes **scalar summaries, not polygons**, but it is the closest algorithmic sibling. The engine's 2×2 classification gives Euler number, perimeter and area essentially for free.

### Cited Findings
- **GPU ECC for imaging data (Stony Brook)**, RTX 2070 8 GB:
  - **5–10 GVox/s** kernel throughput
  - Real-time capability of **512³ at 30 Hz** and **8K (7680×4320) at 120 Hz**
  - 4–6× faster than the CPU tool CHUNKYEuler on files over 0.5 GB
  - Streams **4096³ (0.25 TB)** volumes with 8 GB VRAM by chunking
  - **Disk I/O is 75–99.7% of runtime**; GPU initialization overhead of 100–700 ms makes the CPU faster on files under 16 MB
  - Algorithm: per-voxel "Vector of Changes in the Euler Characteristic" + prefix sum, O(n)
  - [arXiv 2203.09087](https://ar5iv.labs.arxiv.org/html/2203.09087)
- **Tensor-based GPU framework (pyECT)** for weighted Euler characteristic transforms and functions on simplicial and cubical complexes of any dimension, with large speedups over prior methods — [arXiv 2511.03909](https://arxiv.org/abs/2511.03909); vectorized ECF/ECT at SoCG 2026 — [Dagstuhl LIPIcs SoCG 2026](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SoCG.2026.32)
- **Minkowski functionals** (used in cosmology among other fields) have been accelerated with **OpenACC and CUDA** on commodity GPUs — [Abertay publication](https://rke.abertay.ac.uk/en/publications/runtime-analysis-of-spatial-structure-a-cuda-implementation-of-mi/)

### Inferences
- **I/O dominates here too:** in the ECC work, disk is 75–99.7% of runtime. A fused "ML surrogate → threshold → ECC/MF/polygons" path on device avoids it.
- **Cosmology usage:** cosmology non-Gaussianity analyses (e.g. weak-lensing mass maps, CMB) compute MFs at many thresholds × many simulations; emulator suites have thousands of realisations. This is batched 2×2 counting, which the engine's classification pass already does. The value-add of polygons is low; the value-add of the 2×2 counting core is high.

### Gaps
- No 2023–2026 source with specific cosmology MF throughput numbers (e.g. maps per second for Euclid/LSST weak-lensing analyses).

---

## 9. Neuroscience: calcium imaging ROI extraction

### Takeaway
Calcium imaging ROI pipelines (Suite2p, CaImAn) run faster than real time offline. **Online, large-field-of-view** ROI extraction remains hard. ROIs are pixel masks, and polygon output is needed mainly for GUI/curation and closed-loop targeting.

### Cited Findings
- **Suite2p** processes 1 h of 30 Hz 512×512 data in **~20 min with a GPU**, i.e. faster than real time. It registers, detects cells (now via **Cellpose**), extracts traces and infers spikes — [Suite2p (CBMM)](https://cbmm.mit.edu/node/4903); [Suite2p fork](https://github.com/swkeemink/Suite2P)
- **CaImAn** supports **real-time analysis on streaming data** for 1-photon and 2-photon imaging — [search summary](https://pmc.ncbi.nlm.nih.gov/articles/PMC8940911)
- **Online analysis is "challenging, especially for large FOV"**. **CITE-On** (CNN) does online cell identification, segmentation, identity tracking and trace extraction for thousands of cells, including during mesoscopic 2-photon imaging — [Nat. Commun. 2022, s41467-022-29180-0](https://link.springer.com/10.1038/s41467-022-29180-0)

### Inferences
- **Fit:** ROI masks → polygon contours are needed for closed-loop holographic photostimulation targeting and for curation GUIs, but the frame sizes are small (512²). Speed is not the bottleneck; segmentation quality is. This is a low-priority market for a high-throughput polygonizer.

### Gaps
- No 2023–2026 numbers for online ROI extraction latency on mesoscopes (e.g. multi-thousand-cell, multi-plane rigs).

---

## Cross-cutting inferences for the user's engine

These are drawn from the cited findings above.

1. **I/O, not compute, is the bottleneck everywhere it was measured** (TempestExtremes 66% file input; GPU ECC 75–99.7% disk). The pitch is **in-line, on-device fusion** with ML inference or simulation: threshold in device memory, emit compact polygons and attributes, skip writing fields.
2. **Batched multi-threshold and multi-member is the native workload:** ensembles of 1000–7000+ members, eddy contour scans at 2 mm, MFs at many thresholds, PHI at several probability levels. A batched API (many fields × many thresholds per launch) matters more than single-image latency.
3. **Containment and nesting matter:** eddies (nested closed contours), merge trees (all thresholds) and holes in ice floes all need ring hierarchy (outer/inner, parent/child) in the output.
4. **Label-map / planar-partition mode** (shared boundaries, triple junctions) is needed for grains, and likely for ClimateNet-style multi-class masks. Binary occupancy alone does not cover it.
5. **Spans output:** LSST-style Footprints are run-length span sets, and the scanline engine could output these cheaply alongside rings.
6. **Post-processing:** operational polygons need smoothing, simplification and vertex limits (NowCastMIX, NWS). Raw stair-step rings must be simplified, so Douglas-Peucker or topology-preserving simplification on GPU is a needed companion.
7. **Tracking/linking is the other half:** overlap-based (StitchBlobs, PET AreaTracker), distance-based (Earth2Studio 300 km) and watershed+linking (tobac, WDSS-II). Polygon overlap/IoU on GPU between consecutive frames is an adjacent opportunity.
8. **Competitive reference numbers:**
   - Sparse GPU CCL: 300M hits/s (RTX 4070 Ti Super)
   - Dense GPU ECC (2×2/voxel lookup): 5–10 Gvox/s (RTX 2070)
   - HEP triggers: ~500 GPUs for 40 Tbit/s (LHCb); 350×8 GPUs for 3.5 TB/s (ALICE)

## Overall gaps

- Not researched or no sourced data: marine heatwave region extraction, oil-spill polygons, front detection, vortex/flame tracking 2023–2026, Euclid MER segmentation, micro-CT pore networks, phase-field grain tracking, TECA status, MetNet polygon products.
- No direct source found for "store contours instead of fields" as an explicit exascale reduction strategy with 2023–2026 volume numbers. The nearest found evidence is ALPINE data-reduction goals, segmented merge trees and topology-preserving compression.
