# Spatial Analysis & Spatial ML Beyond Imagery: Points / Trajectories / Fields -> Grid -> Polygons (2023–2026)

Scope note: non-imagery only. Context engine = CUDA/C++ binary-occupancy-bitmap -> closed polygon rings, with tile halos, pluggable transform, JTS-like GeometryFactory. Research done 2026-10-07. Every claim below carries a source; unsourced domain knowledge is confined to "Inferences" and labelled as such.

## Q1. Density / hotspot / field -> polygon use cases (KDE hotspots, isochrones, catchments, telecom coverage, fire perimeters)

### Takeaway
The "scalar field on a grid -> threshold -> polygon rings" pattern is the core of isochrone engines (Valhalla builds a lat/lon cost grid then contours it), national telecom coverage reporting (FCC BDC ships both propagation-model polygons and H3 res-9 hexagon versions), wildfire simulators (arrival-time rasters + perimeter polygons), and KDE crime/disease hotspot mapping. In most of these the vectorization step is a CPU marching-squares or GDAL-style step bolted onto the end of the pipeline, often followed by Douglas-Peucker simplification that can break topology.

### Cited Findings
**Isochrones / service areas**
- Valhalla (Mapzen's routing engine) builds isochrones by "first creating a 2-D grid in latitude/longitude about the location"; the grid is filled by a Dijkstra-style least-cost search that marks cells touched by road edges with time/cost, then contoured. This is the canonical "graph search -> raster -> contour polygon" pipeline. (Source is the original Mapzen launch post, pre-2023, but it describes the architecture Valhalla still uses.) — [Mapzen blog: Introducing isochrone service](https://mapzen.com/blog/introducing-isochrone-service)
- Mapbox Isochrone API limits: max 4 contours per request, max 60 minutes or 100 km, 300 requests/minute. `polygons=true` returns GeoJSON polygons ("Any contour that forms a ring is returned as a polygon"). A `denoise` parameter (0–1) drops small contours (1.0 keeps only the largest). A `generalize` parameter applies Douglas-Peucker, which the docs say "can lead to self-intersections, as well as intersections of adjacent contours." — [Mapbox Isochrone API docs](https://docs.mapbox.com/api/navigation/isochrone/)
- Marching squares is the standard algorithm for turning a 2D scalar field into contour polygons. Libraries cited in the search results include d3-contour (JS), contour-isobands (Rust) and silx (Python). — [Wikipedia: Marching squares](https://en.wikipedia.org/wiki/Marching_squares); [d3-contour](https://github.com/d3/d3-contour); [contour-isobands crate](https://docs.rs/contour-isobands/)
- CARTO published "Creating isolines with CARTO" (April 2024). Only the title and date were retrievable, not the content. — [plushcap summary of CARTO post](https://plushcap.com/analysis/carto/carto-creating-isolines-with-carto)

**Telecom coverage (FCC Broadband Data Collection)**
- FCC BDC mobile broadband and voice availability data comes as Shapefile/GeoPackage. Most datasets include **both** "a hexagonal representation of mobile coverage (using the H3 Spatial Index at resolution 9)" **and** "the raw, propagation modeled mobile coverage." — [FCC BDC Help: How to download mobile coverage data](https://help.bdc.fcc.gov/hc/en-us/articles/43909220634651-How-to-Download-Mobile-Broadband-Coverage-Data-from-the-FCC-s-National-Broadband-Map-Step-by-Step-Instructions)
- Providers must submit coverage as **polygon GIS data** derived from propagation modeling, with fixed parameters per technology/speed tier under the Broadband DATA Act. — [FCC BDC: Formatting Mobile Broadband Availability Coverage Maps](https://help.bdc.fcc.gov/hc/en-us/articles/6047425308187-Formatting-Mobile-Broadband-Availability-Coverage-Maps)
- Propagation models must use terrain and clutter data at 100 m resolution or better. The FCC uses nested H3 hexagon resolutions for verification and challenges. — [TLP law summary of FCC BDC technical webinar](https://www.tlp.law/2021/08/18/fcc-webinar-on-broadband-data-collection-technical-requirements/); [Geo-Tel: hexagonal grids for BDC](https://geo-tel.com/?p=5312)
- Utah's state GIS office posted new broadband service data in Nov 2025, which shows downstream state reuse of the BDC data. — [Utah gis.utah.gov blog 2025-11-04](https://gis.utah.gov/blog/2025-11-04-new-broadband-service-data/)

**Wildfire spread perimeters**
- FARSITE propagates the fire front as a **vector** and outputs perimeter polygons at set time steps. Raster outputs (time of arrival, fireline intensity, flame length, ROS, etc.) are interpolated from those perimeters using inverse-distance weighting. — [FARSITE tech docs: Raster Data Generation](https://owfflammaphelp62.firenet.gov/LegacyFARSITE/Tech_Raster_Data_Generation.htm); [USFS RMRS-RP-4 (Finney)](https://www.fs.usda.gov/rm/pubs/rmrs_rp004.pdf)
- Arrival time is defined as a raster of minutes from ignition to each cell. Thresholding it at time t gives the perimeter at time t, which is a raster -> polygon step. — [FlamMap/FARSITE help: Arrival Time](https://owfflammaphelp62.firenet.gov/FileTypes/Arrival_Time.htm)
- 2025 work on wildfire spread **data assimilation** (Fire-EnSF, ensemble score filter) operates on gridded fire-state fields. — [arXiv 2510.15954](https://arxiv.org/pdf/2510.15954)

**KDE crime / disease hotspots**
- Feb 2026 arXiv "Predictive Hotspot Mapping for Data-driven Crime Prediction" (IIM Ahmedabad, Delhi police-call data Oct 2019–Mar 2021) computes density on **36,263 grid boxes of ~200 m x 200 m** and labels the top 20% and 40% of boxes "red" and "yellow" zones. It never vectorizes them into smooth polygons. One computation round takes ">12 hours" on a standard PC and "<1 hour" on HPC. — [arXiv 2602.23750](https://arxiv.org/html/2602.23750v1)
- KDV-Explorer (HKBU) is a near-real-time kernel density visualization system aimed at disease outbreak, crime and traffic-accident hotspots. — [KDV-Explorer](https://www.comp.hkbu.edu.hk/~edisonchan/system/KDV-Explorer.html)
- Practitioner workflow in Esri forums: run Kernel Density, then reclassify, then Raster-to-Polygon (giving blocky rings), or contour a single threshold interval as "Regions". A Python script in the results uses `skimage.measure.find_contours` specifically to avoid "blocky boundaries from merging rectilinear grid cells." — [Esri Community: density and convert to polygon](https://community.esri.com/t5/geoprocessing-questions/density-and-convert-to-polygon/m-p/229391); [Miami U. extract_hotspot_contours.py](https://vision.csi.miamioh.edu/data/heatmap/scripts/extract_hotspot_contours.py)
- UCGIS GIS&T BoK frames raster->vector "vectorization" and hot-spot region identification as core GIS operations. — [GIS&T BoK](https://gistbok-ltb.ucgis.org/30/tag/show/26614)

### Inferences
- The engine's input contract (binary occupancy bitmap) matches the "threshold the field first" step. Isochrones, coverage, arrival-time perimeters and KDE hotspots all reduce to one or a few thresholded bitmaps per output band. Multi-band isochrones (e.g., 4 contours) would need either N binary passes or a native multi-level (isoband) mode. (Inference)
- Mapbox's own docs admit that its generalization can produce self-intersecting and mutually intersecting contours. A topology-correct, nested isoband output (rings that share edges across levels) is a concrete quality gap a GPU tracer could fill. (Inference from the Mapbox doc)
- The FCC dual product (raw propagation polygons plus H3 res-9 hexes) shows a regulatory workflow where national-scale propagation rasters must become polygons, and polygons must become hex cells, both ways. (Inference)

### Gaps
- Could not verify whether OSRM, GraphHopper or Valhalla currently use GPU anywhere in isochrone contouring. No evidence found; assume CPU.
- No public throughput numbers were found for FCC BDC submission processing (e.g., number of polygons or vertices per carrier filing).
- A search snippet mentioned "Marching Squares KDE (MSKDE)" for smooth crime-hotspot contours, but I could not trace it to a primary source. Treat as unverified.
- No sourced findings on noise maps (EU END strategic noise maps) or air-quality exceedance-zone polygons. These are likely the same raster->isoband pattern but are unverified here.
- Wastewater-surveillance catchment polygons: not researched within budget.

## Q2. Discrete global grids (H3, S2, DGGS): cells -> polygons at scale

### Takeaway
"Set of hex cells -> outline polygon" is a live, actively patched pain point. H3 core restructured it in v4.5.0 (May 2026). The Go binding was still being optimized in Sept 2026. Correctness bugs remain open: degenerate 1-vertex polygons and failures near the poles. In SQL engines the common workaround is `ST_Union`/`ST_CoverageUnion` over cell hexagons, which is orders of magnitude slower than a native boundary trace.

### Cited Findings
- **H3 v4.5.0 (May 21, 2026)** added `cellsToMultiPolygon` and refactored `cellsToLinkedMultiPolygon` to use it internally. Both now return errors for invalid cells (`E_CELL_INVALID`), mixed resolutions (`E_RES_MISMATCH`) and duplicate cells (`E_DOMAIN`) "instead of producing undefined output." This means mixed-resolution (compacted) sets are **rejected**, not outlined. — [uber/h3 releases](https://github.com/uber/h3/releases); [v4.5.0](https://github.com/uber/h3/releases/tag/v4.5.0)
- **H3 v4.2.0 (Dec 4, 2024)** added `polygonToCellsExperimental`, which supports more containment modes and uses less memory. This is the inverse direction, polygon -> cells. — [uber/h3 releases](https://github.com/uber/h3/releases)
- **h3-go PR #141 (Sept 28, 2026)** replaced Go maps with a flat linear-probing hash table keyed by DirectedEdge for `CellsToMultiPolygon`. On an Apple M3 Max: 46.25 µs -> 33.18 µs (-28%), allocations 188 -> 8 per op (-96%). It is now about 10% slower than C (33.2 vs 30.1 µs). This confirms the algorithm is a CPU hash-map edge-cancellation over directed edges. — [uber/h3-go PR #141](https://github.com/uber/h3-go/pull/141)
- h3-js benchmarks list cellsToMultiPolygon at about **760 ops/sec**, slow next to other H3 functions. A long-standing issue (#71) came from users who wanted to merge "many thousands of hexagons" into multipolygons to speed up rendering. — [uber/h3-js](https://github.com/uber/h3-js); [h3-js issue #71](https://github.com/uber/h3-js/issues/71)
- **H3 issue #1049**: 164 valid res-2 cells -> 14 polygons, 2 of which were degenerate with only 1 vertex. — [uber/h3 #1049](https://github.com/uber/h3/issues/1049)
- **H3 issue #1029**: `cells_to_multi_polygon` "will fail to produce coherent polygon output" near the poles and at high latitudes. — [uber/h3 #1029](https://github.com/uber/h3/issues/1029)
- ClickHouse PR #123911 speeds up `h3PolygonToCells` for small polygons, which shows DB vendors optimizing H3 polygon conversions directly. — [ClickHouse PR #123911](https://github.com/ClickHouse/ClickHouse/pull/123911)
- PostGIS options for dissolving cell polygons: `ST_Union` is "orders of magnitude slower than ST_Collect". `ST_CoverageUnion` skips overlay for valid coverages and is "typically one to two orders of magnitude faster" ("minutes versus seconds" at 100k polygons, according to an aggregator page). `ST_MemUnion` trades speed for memory. — [PostGIS ST_Union docs](https://postgis.net/docs/manual-3.0/ST_Union.html); [atlas.co ST_CoverageUnion](https://atlas.co/functions/postgis/st-coverage-union/) (aggregator, not primary; numbers unverified against PostGIS benchmarks); [PostGIS-devel ST_Union perf thread](https://lists.osgeo.org/pipermail/postgis-devel/2013-October/023839.html)
- DuckDB has an H3 community extension (wrapped in R by `duckh3`, 2025). No perf data was found for its cells->polygon WKT function. — [duckh3 README (CRAN)](https://cran.r-project.org/web/packages/duckh3/readme/README.html); [Mark Litwintschik: H3 + DuckDB + QGIS](https://tech.marksblogg.com/h3-duckdb-qgis.html)
- Uber built H3 to bucket marketplace events into hexagons for **surge pricing** and dispatch: "measuring supply and demand in hexagons in each city." — [Uber blog: H3 (2018)](https://www.uber.com/blog/h3/)
- FCC BDC standardizes on H3 res-9 for mobile coverage hexes (see Q1). — [FCC BDC Help](https://help.bdc.fcc.gov/hc/en-us/articles/43909220634651-How-to-Download-Mobile-Broadband-Coverage-Data-from-the-FCC-s-National-Broadband-Map-Step-by-Step-Instructions)

### Inferences
- **Direct fit for the engine:** a set of H3/S2 cells at one resolution is a sparse occupancy map on a non-Cartesian lattice. Two routes look viable. (a) Rasterize cells into a bitmap in a local projection, then trace. The output is a pixel-staircase approximation, not exact hex edges. (b) Add a hex-lattice boundary tracer that emits exact cell-vertex rings. H3 already does (b) on CPU via directed-edge cancellation (h3-go PR). A GPU edge-cancellation or "boundary edge mark + ring link" kernel is a natural port. (Inference)
- The H3 refusal to handle mixed resolutions means compacted sets must be uncompacted first, which inflates cell counts by 7^k. That inflation is exactly where GPU throughput would matter. (Inference)
- Pole and antimeridian failures (#1029) show that the "pluggable coordinate transform" must handle the sphere correctly. An equirectangular bitmap will hit the same issue. (Inference)

### Gaps
- No public benchmark found for cellsToMultiPolygon at 10^6–10^8 cells (the scale of national coverage or global AIS grids).
- S2 (`S2CellUnion` -> `S2Polygon`) and OGC DGGS API (OGC API – DGGS) conversion performance: not found within budget.
- No GPU implementation of H3 cells->polygon was found (searched cuSpatial; it has no such function, see Q5).

## Q3. Mobility, trajectories, dynamic zones (delivery, surge, UTM/drone geofences, maritime)

### Takeaway
Mobility platforms run on hex grids (Uber surge). Delivery "dynamic zones" are patented and described publicly only loosely. Drone/UTM research is moving to fine 3D airspace grids with dynamic occupancy that must be published as geofence volumes. Maritime analytics (Global Fishing Watch) relies on gridded AIS/SAR activity density. All of these produce grids or occupancy that downstream consumers want as polygons.

### Cited Findings
- **Uber surge**: supply/demand measured per H3 hexagon per city (see Q2). — [Uber H3 blog](https://www.uber.com/blog/h3/)
- **DoorDash dynamic delivery zones**: US patent US10176448B1, "Generation of dynamic delivery zones for merchants." A GTC 2020 NVIDIA talk covers ML for DoorDash logistics. I found no 2023–2026 DoorDash engineering post on zone geometry. — [Google Patents US10176448B1](https://patents.google.com/patent/US10176448); [NVIDIA GTC 2020 S22191 slides](https://developer.download.nvidia.com/video/gputechconf/gtc/2020/presentations/s22191-how-machine-learning-powers-on-demand-logistics-at-doordash.pdf)
- **UTM / U-space**: An ISPRS Congress 2026 paper (Zhao et al., Shenzhen University, July 2026) proposes "fine-grained urban low-altitude airspace gridding with dynamic event response." It uses a unified **3D grid** on a city digital twin. Static and dynamic attributes are injected into the cells: terrain, buildings, no-fly zones, wind fields, temporary restrictions, "occupancy and release information." The case study is Yuehai Subdistrict, Nanshan, Shenzhen. — [ISPRS Archives XLIX-B2-2026-611](https://isprs-archives.copernicus.org/articles/XLIX-B2-2026/611/2026/isprs-archives-XLIX-B2-2026-611-2026.pdf)
- U-space phases: U1 foundation services include **geofencing**, and U3 (advanced) adds capacity management and dynamic airspace reconfiguration. A Cyprus 2025 funding call targets "U3 U-space advanced services," including enhanced geofencing. — [Cyprus Funding Programmes Portal 2025](https://www.fundingprogrammesportal.gov.cy/en/call/fast-track-u3-u-space-advanced-services-and-cns-capabilities-agglika-2025/)
- Dynamic geofences "appear and disappear at user-specified times and are spatially and temporally deconflicted to maximize authorized airspace volume" (U. Michigan, Atkins group). — [Deep Blue: Geofencing for sUAS in complex low-altitude airspace](https://deepblue.lib.umich.edu/items/b590e1ba-1493-40fa-8b5e-9bfdafb4abdb); [SESAR GEOSAFE no-fly-zone demo](https://www.unmannedairspace.info/uncategorized/sesar-geosafe-project-demonstrates-no-fly-zones/)
- **Maritime**: Global Fishing Watch's Nature paper (Jan 2024) combined SAR and optical satellite detections, AIS and deep learning over the busiest ~15% of the ocean for 5 years. It found ~75% of industrial fishing vessels are not publicly tracked, with >90% classification accuracy. GFW publishes AIS-derived fishing effort as gridded density. Note that the detections are imagery-based; the AIS density grids are not. — [Oceana summary](https://oceana.org/?p=29432); [Maritime Executive](https://www.maritime-executive.com/article/study-75-of-the-world-s-fishing-activity-is-invisible-on-ais); [Esri ArcNews on GFW](https://www.esri.com/about/newsroom/arcnews/illuminate-fishing-activity-on-the-open-ocean)

### Inferences
- UTM grids ("occupancy and release" per cell) are literally occupancy rasters, in 3D and updated per event. Publishing them as geofence polygons or volumes to operators needs fast incremental bitmap->polygon per altitude layer. Halo/tile handling maps onto spatial partitions of a city. (Inference)
- Fishing-ground or MPA-compliance zones derived from AIS density grids are a "threshold the density, then polygonize" task at global scale. Antimeridian handling matters here. (Inference)
- Stay-point/dwell-zone polygons from trajectories are usually built by DBSCAN plus convex or concave hull. A rasterize-and-trace route would give non-convex shapes with holes. (Inference; no 2023–2026 source verified)

### Gaps
- No 2023–2026 public engineering posts found on DoorDash, Instacart or Lyft zone-geometry generation. These are likely internal.
- FAA Remote ID / LAANC grid (the UAS Facility Maps are pre-gridded at roughly 1/30th of a degree, from domain knowledge) was not verified in this session.
- Grid cell sizes in the ISPRS 2026 airspace paper: abstract read only; exact numbers not extracted.

## Q4. Spatial ML: regionalization, location/polygon embeddings, redistricting, foundation models

### Takeaway
Spatial ML has two polygon-relevant directions. (1) Regionalization and redistricting produce contiguous unit-level labelings, which are graph partitions whose dissolved outlines are polygons. (2) Location/polygon encoders (SatCLIP, GeoCLIP, Poly2Vec at ICML 2025) and dense embedding fields (AlphaEarth, 64-D at 10 m, 1.4 trillion footprints per year) produce dense per-pixel or per-location outputs that clustering turns into label rasters, which must be polygonized to become zones.

### Cited Findings
- **Max-p regionalization** is NP-hard. Approximate algorithms "do not scale for large datasets." Two scalable frameworks: SMP (SIGSPATIAL 2022: randomized initial solution, then heuristic improvement) and PAGE (ACM TSAS, Sept 2023: a parallel regionalization framework that extends spatial partitioning without violating contiguity). A 2024 SIGSPATIAL paper extends this to "scalable and expressive spatial grouping queries," covering regionalization and hotspot detection. — [SMP (NSF PAR)](https://par.nsf.gov/servlets/purl/10357505); [PAGE, ACM TSAS 2023](https://doi.org/10.1145/3609956.3609980); [SIGSPATIAL 2024 spatial grouping](https://unpaywall.org/10.1145%2F3678717.3695757)
- **Poly2Vec** (Siampou, Li, Krumm, Shahabi, Lu; ICML 2025) uses a 2D Fourier-transform encoding that unifies points, polylines and polygons into fixed-length vectors with a learned magnitude/phase fusion. It beats object-specific baselines on topology, direction and distance relations, and improves population prediction and land-use inference. — [PMLR v267](https://proceedings.mlr.press/v267/siampou25a.html); [arXiv 2408.14806](https://arxiv.org/abs/2408.14806v2); [USC InfoLab](https://infolab.usc.edu/research/poly2vec/)
- **SatCLIP** trains a location encoder Siren(SH(lat,lon)) with spherical-harmonic positional encodings, contrastively against Sentinel-2. **GeoCLIP** uses hierarchical random Fourier features. — [SatCLIP arXiv 2311.17179](https://arxiv.org/pdf/2311.17179); [Intrinsic Dimension of Earth Representations, arXiv 2511.02101](https://arxiv.org/pdf/2511.02101)
- **AlphaEarth Foundations / Satellite Embedding V1** (Google/DeepMind, July 2025): global 10 m pixels, each a unit-length **64-D** annual embedding, covering 2017–2024. It is designed for clustering, classification and change detection, has >1.4 trillion embedding footprints per year, and is also on GCS at `gs://alphaearth_foundations`. Although derived from EO data, its use is non-imagery ML: people cluster embeddings, then need zones. — [GEE catalog](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_SATELLITE_EMBEDDING_V1_ANNUAL); [DeepMind blog](https://deepmind.google/blog/alphaearth-foundations-helps-map-our-planet-in-unprecedented-detail/?hl=ca)
- **Redistricting**: ReCom Markov chains on census-block adjacency graphs (GerryChain). Ensemble evidence has been used in court in at least 11 states since the 2020 Census. A 2024 Rust ReCom/RevReCom implementation enables "block-level runs on large states in a matter of hours." Each sampled plan is a block-level labeling that must be dissolved into district polygons for display and metrics such as compactness. — [UChicago Data Science Institute: Gerry-Suite updates](https://democracy.datascience.uchicago.edu/posts/gerrysuite); [Samplers page](https://democracy.datascience.uchicago.edu/samplers); [arXiv 2210.01401](https://arxiv.org/html/2210.01401v2)

### Inferences
- **ML-to-polygon handoff:** k-means on AlphaEarth embeddings, GNN node classification over H3 cells, or regionalization over grid cells all output a **label raster or label grid**. The engine handles one binary bitmap, so multi-class output needs either per-class passes (shared edges are then duplicated and may not match) or a multi-label "coverage" mode that emits shared-edge topology. A coverage mode would be a differentiator. (Inference)
- Redistricting ensembles generate thousands to millions of plans. Compactness scores (Polsby-Popper and similar) need district perimeters. If blocks are pre-rasterized, a GPU tracer could compute perimeters and rings in batch, but blocks are irregular polygons, so a graph-dissolve approach probably dominates. Weak fit unless the engine supports label maps on a fine raster. (Inference)
- Poly2Vec-style encoders consume polygons as input. Fast raster->polygon allows ML-derived zones to be re-encoded as features in a loop. (Inference, speculative)

### Gaps
- No GPU-accelerated max-p or SKATER implementation was found (searched; only CPU/parallel frameworks such as PAGE).
- Spatial foundation models over **vector** data beyond Poly2Vec (e.g., Geo2Vec, OSM-based foundation models) were not surveyed in depth. Space2Vec appears only as a baseline lineage.
- No source found describing how practitioners polygonize AlphaEarth cluster maps (presumably GEE `reduceToVectors`, which is unverified here).

## Q5. Platforms & GPU stacks: what they do for raster/grid -> polygon

### Takeaway
GPU spatial stacks (deck.gl, Kinetica, HeavyDB, cuSpatial) GPU-accelerate **aggregation and rendering** but not exact vector polygon extraction. deck.gl's ContourLayer aggregates on GPU, while its marching-squares location is undocumented and was originally CPU. Kinetica renders contours and isochrones server-side as images. cuSpatial has no polygonize. Sedona/SedonaDB offer `RS_Polygonize` and per-pixel `RS_PixelAsPolygons` on CPU. Outside the browser, no mainstream platform exposes a GPU bitmap->polygon-ring primitive.

### Cited Findings
- **deck.gl ContourLayer**: `gpuAggregation` defaults to `true` ("aggregation is performed on GPU" when supported). The data is binned by `cellSize`, then "resulting scalar field is used to run Marching Squares." The docs do not say whether marching squares runs on GPU. Isolines feed LineLayer and isobands are filled polygons. — [deck.gl ContourLayer docs](https://deck.gl/docs/api-reference/aggregation-layers/contour-layer)
- deck.gl v6.1 ContourLayer RFC proposed marching squares that "parse a texture once for each threshold" to generate vertex data, transformed to lng/lat via cellSize/gridOrigin. This was a design proposal; current implementation status was not verified. — [deck.gl RFC v6.1 contour-layer](https://github.com/visgl/deck.gl/blob/9a77f7d30b6b5a17313a120aff47f85bd7814daa/dev-docs/RFCs/v6.1/contour-layer-rfc.md); [GPU grid layer RFC v7.1](https://github.com/visgl/deck.gl/blob/9a77f7d30b6b5a17313a120aff47f85bd7814daa/dev-docs/RFCs/v7.1/gpu-grid-layer.md)
- **Apache Sedona / SedonaDB**: `RS_Polygonize(rast, band)` returns one polygon (WKB) per connected same-value region. Connectivity rule and implementation are not documented on the page. `RS_PixelAsPolygons` (since v1.5.1) returns one polygon per pixel, which is the naive and explosive approach. — [SedonaDB RS_Polygonize](https://sedona.apache.org/sedonadb/latest/reference/sql/rs_polygonize/); [Sedona RS_PixelAsPolygons 1.9.1](https://sedona.apache.org/1.9.1/api/sql/Pixel-Functions/RS_PixelAsPolygons/)
- **CARTO**: H3 point aggregation layers bin points dynamically into hexes for rendering. The BigQuery Tiler handles "millions or even billions of rows" in-warehouse. — [CARTO H3 point aggregation](https://docs.carto.com/carto-user-manual/maps/layers/point/h3-point-aggregation); [CARTO BigQuery codelab](https://codelabs.carto.com/bigquery)
- **Kinetica**: 150+ PostGIS-compatible functions, H3 indexing, track operators, all on GPU. Server-side rendering modes include raster, **contour (isolines) and isochrones**. These are rendered outputs; I found no documented SQL function that returns contour polygons. — [Kinetica geospatial analytics](https://www.kinetica.com/features/geospatial-analytics)
- **HeavyDB**: GPU OGC types/operators. HeavyRender produces heatmaps and choropleths. No raster->polygon function was found. — [HEAVY.AI HeavyDB](https://heavy.ai/product/heavydb)
- **cuSpatial (RAPIDS)**: GPU point-in-polygon, spatial joins, trajectory reconstruction, haversine, GeoArrow. It has **no raster polygonize or contour function**. Academic GPU contour vectorization exists (boundary detection -> connectivity -> loop organization -> segment vectorization). — [cuSpatial docs](https://docs.rapids.ai/api/cuspatial); [PKU GPU vectorization paper](https://www.icst.pku.edu.cn/GraphicsLab/docs/20200416023556073458.pdf); [hgpu entry](https://hgpu.org/?p=6790)

### Inferences
- **Clear whitespace:** GPU databases and viz stacks stop at "aggregate on GPU, then render pixels or CPU contours." A CUDA library returning topologically valid rings with tile halos would slot in as (a) a Sedona/SedonaDB UDF replacing `RS_Polygonize`, (b) a cuSpatial-adjacent RAPIDS component (GeoArrow output), (c) a Kinetica/HeavyDB-style SQL function `ST_ContourPolygons(grid, threshold)`, or (d) a server-side producer of isoband vector tiles for deck.gl/kepler.gl. (Inference)
- The engine's tile+halo design matches vector-tile producers (CARTO Tiler, deck.gl MVT). Seamless per-tile isoband polygons are a known rendering need. (Inference)

### Gaps
- Not verified: Esri GeoAnalytics Engine, Foursquare (Unfolded/kepler.gl), BigQuery GIS, Wherobots-specific raster->vector functions, and whether kepler.gl's HexagonLayer/Contour use GPU marching squares.
- No benchmark data comparing GDAL Polygonize vs. Sedona `RS_Polygonize` vs. any GPU method.

## Q6. Pain points and gaps reported by practitioners

### Takeaway
Recurring complaints include blocky staircase polygons from cell-merge vectorization, slow or incorrect hex-set dissolves (H3 degenerate rings, pole failures, mixed-resolution rejection, slow JS path), topology breakage after simplification (Mapbox self-intersections), and hours-long hotspot pipelines. These map directly onto an engine that is fast, topology-aware, tile-seam-aware and transform-pluggable.

### Cited Findings
- Blocky boundaries from merging grid cells are why practitioners switch to `skimage.find_contours` for KDE hotspot polygons. — [extract_hotspot_contours.py](https://vision.csi.miamioh.edu/data/heatmap/scripts/extract_hotspot_contours.py)
- Douglas-Peucker generalization "can lead to self-intersections, as well as intersections of adjacent contours" (Mapbox, a production isochrone service). — [Mapbox Isochrone docs](https://docs.mapbox.com/api/navigation/isochrone/)
- H3 cells->polygon: degenerate 1-vertex rings (#1049), incoherent output near poles (#1029), mixed resolutions now hard-error (v4.5.0), JS path about 760 ops/sec, Go binding needing a 2026 allocation rewrite. — [#1049](https://github.com/uber/h3/issues/1049); [#1029](https://github.com/uber/h3/issues/1029); [v4.5.0](https://github.com/uber/h3/releases/tag/v4.5.0); [h3-js](https://github.com/uber/h3-js); [h3-go #141](https://github.com/uber/h3-go/pull/141)
- Users merging thousands of hexes for faster rendering (h3-js #71). — [h3-js #71](https://github.com/uber/h3-js/issues/71)
- ST_Union dissolves of large coverages are slow, and PostGIS has a long history of ST_Union performance complaints, including a funding offer on postgis-devel. — [postgis-devel 2013](https://lists.osgeo.org/pipermail/postgis-devel/2013-October/023839.html); [PostGIS trac #2354](https://trac.osgeo.org/postgis/ticket/2354)
- Crime hotspot computation takes >12 h per round on a PC for one city (Delhi) and <1 h on HPC (2026). — [arXiv 2602.23750](https://arxiv.org/html/2602.23750v1)

### Inferences
- Use-case priority for this engine, by fit and pain (inference):
  1. Isochrone/isoband polygonization (routing engines, Mapbox-like APIs): topology-correct nested bands.
  2. Hex/DGGS set -> outline at 10^6+ cells: coverage maps, H3 analytics, FCC-style filings.
  3. Real-time occupancy -> geofence polygons: UTM/U-space 3D grids, per altitude slice.
  4. ML label-raster -> zones: AlphaEarth clustering, GNN or regionalization outputs. Needs a multi-label coverage mode.
  5. Simulation fields -> perimeters: wildfire arrival time, flood depth, plume exceedance.
- Feature asks implied by these pain points: multi-threshold isobands with shared edges; hex-lattice or DGGS-aware tracing (or exact cell-vertex rings); antimeridian and pole-safe transforms; validity guarantees (no degenerate rings, correct hole nesting); seam-free tile stitching; topology-preserving simplification. (Inference)

### Gaps
- Reddit r/gis and GIS StackExchange threads specifically about tile-seam artifacts when polygonizing tiled heatmaps were not found within budget. The seam pain point is plausible but unsourced here.
- No quantitative survey of how often analysts hit raster->polygon bottlenecks.
