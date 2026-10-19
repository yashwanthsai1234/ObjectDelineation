# Geospatial and Earth Observation: Where Raster-Mask-to-Polygon Conversion Is Used and Where It Hurts (2023–2026)

Context for the report writer: the engine under evaluation is a CUDA/C++ single-pass scanline polygonizer. It uses 2x2 pixel-type classification, takes a binary occupancy bitmap with a halo for tiling, and emits closed rectilinear rings through a pluggable (georeferencing) coordinate transform. These notes cover the geospatial/EO branch only. Research ran on 2026-10-07; tool-call budget was about 25.

## Q1. Which pipelines turn segmentation masks into polygons at scale, who runs them, and how large are they?

### Takeaway
Planet-scale "segment, then polygonize" pipelines now produce billions of polygons per product. Fields of The World has about 3.2B field polygons. Overture has 2.3B+ buildings, Microsoft 1.4B, and GlobalBuildingAtlas 2.75B. Polygonization is consistently a named stage of each pipeline, and almost all of them run CPU tooling (GDAL/rasterio) after GPU inference. Annual refreshes (FTW does 2024 and 2025; Google 2.5D Temporal covers 2016–2023) multiply the volume.

### Cited Findings
**Agricultural field boundaries**
- Fields of The World (FTW) Global is described as "the first global, wall-to-wall agricultural field-boundary dataset at 10 m resolution with 3.17 billion field polygons across 241 countries and territories for 2024–2025". It is CC-BY-4.0, from Taylor Geospatial with Microsoft AI for Good, ASU, WashU, Oregon State and Clark, hosted on Source Cooperative — [Source Cooperative llms.txt](https://source.coop/ftw/global-data/llms.txt)
- That is about 1.6B polygons per year: 1.63B for 2024 and 1.58B for 2025. They were produced by the PRUE model (U-Net with an EfficientNet-B7 encoder) on cloud-free Sentinel-2 planting/harvest composites — [FTW README](https://data.source.coop/ftw/global-data/README_next.md)
- Vectors are made by "thresholding the softmax outputs for [non_field_background, field, field_boundaries] at 0.5 and polygonizing". They are published as GeoParquet in 598 country/subdivision partitions (~3.2B polygons, 219 countries). The upstream predictions are stored as Zarr — [FTW vectors llms.txt](https://data.source.coop/ftw/global-data/predictions/vectors/llms.txt)
  - Count discrepancy: 241 vs 219 countries across the two pages.
- The FTW "Field Guide" paper (arXiv 2602.08131) gives the recipe: threshold softmax at 0.5, morphological opening, then "convert rasters to polygons via connected components" (`ftw inference polygonize filtered.tif --out boundaries.gpkg`).
  - Optional filtering with ESA WorldCover masks.
  - Country-scale inference uses a 256×256 sliding window with 25% overlap and Gaussian-weighted (apodization) merging to suppress edge artifacts.
  - The paper gives no runtime or compute figures for polygonization. Median field sizes are only 0.06–0.28 ha by region, so polygons are very numerous and very small — [arXiv 2602.08131](https://arxiv.org/html/2602.08131v1)
- The FTW benchmark has 70,462 samples in 24 countries on four continents — [arXiv 2409.16252](https://arxiv.org/pdf/2409.16252); [FTW llms.txt](https://source.coop/ftw/global-data/llms.txt)
- Wherobots worked with Taylor Geospatial Engine to run the FTW model at scale on an early release of RasterFlow, its serverless inference engine (private preview announced 2025-12-10). RasterFlow "predictions are then vectorized into geometries and made available as an Iceberg table" — [Wherobots RasterFlow blog](https://wherobots.com/blog/rasterflow-earth-observation-inference-engine/); [press release](https://business.minstercommunitypost.com/minstercommunitypost/article/accwirecq-2025-12-10-wherobots-introduces-rasterflow-a-major-step-forward-in-earth-intelligence-that-makes-complex-earth-observation-data-simple-and-ai-ready)

**Building footprints**
- Microsoft Global ML Building Footprints: 1.4B footprints with tens of millions of height estimates after 25 refreshes. Sources are Airbus, Maxar, Vexcel and IGN FI imagery — [Mark Litwintschik blog, "ms-buildings-2026"](https://tech.marksblogg.com/ms-buildings-2026.html)
  - The pipeline has two stages: (1) DNN semantic segmentation, then (2) "Polygonization: converting building pixel blobs into polygons", using a custom method that approximates prediction pixels into polygons "based on the whole prediction feature space" — [same](https://tech.marksblogg.com/microsofts-global-ml-building-footprints.html)
  - Distribution: 113 GB across 30,344 gzipped CSVs. It is also on Planetary Computer as GeoParquet/Delta partitioned by Bing quadkey level 9 — [same](https://tech.marksblogg.com/microsofts-global-ml-building-footprints.html); [Planetary Computer STAC](https://planetarycomputer.microsoft.com/api/stac/v1/collections/ms-buildings)
  - Precision by region is 92–97%. Recall is only 70.9–85.9% — [same](https://tech.marksblogg.com/microsofts-global-ml-building-footprints.html)
- Google Open Buildings 2.5D Temporal gives building presence, confidence and height at an effective ~1 m resolution (raster delivered at 0.5 m). It covers 2016–2023 with annual updates and is derived from Sentinel-2. It ships as a raster, so users have to polygonize or join it to V3 polygons themselves to get yearly vector datasets — [Earth Engine catalog](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_Research_open-buildings-temporal_v1?hl=zh-cn); [Spatial Thoughts 2025](https://spatialthoughts.com/2025/03/29/building_height_gee/)
- The VIDA merged Google V3 + Microsoft + OSM dataset holds 2,705,459,584 footprints — [Source Cooperative / VIDA](https://source.coop/vida/google-microsoft-osm-open-buildings)
- Overture Maps buildings: 2.3B+ buildings as of December 2023. Sources are OSM (~600M), Esri Community Maps (~14M), Google ML roofprints (~400M) and Microsoft ML roofprints (~600M) — [Overture blog 2023](https://overturemaps.org/blog/2023/overture-buildings-theme-hits-2-3b-buildings-with-addition-of-google-open-buildings-data/)
  - The 2025-01-22 release added 228.5M ML-derived East Asian buildings — [Overture release notes 2025-01-22](https://docs.overturemaps.org/blog/2025/01/22/release-notes/)
  - The buildings theme reached GA in 2025-03-19 — [release notes 2025-03-19](https://docs.overturemaps.org/blog/2025/03/19/release-notes/)
- GlobalBuildingAtlas (TUM, arXiv 2506.04106, 2025): more than 2.75 billion building polygons, built with ML pipelines on PlanetScope imagery plus quality-based fusion with open datasets — [arXiv 2506.04106](https://arxiv.org/abs/2506.04106)
  - It includes 2.68B LoD1 3D instances, 3×3 m height resolution ("30 times finer than previous global products (90 m)"), and height RMSE of 1.5–8.9 m — [same](https://arxiv.org/abs/2506.04106)

**Floods and disasters**
- Copernicus EMS Global Flood Monitoring (GFM) is a fully automated, 24/7 global service. It processes the entire Sentinel-1 SAR influx with three flood-detection algorithms in parallel to produce flood and water extent maps — [GFM product launch](https://global-flood.emergency.copernicus.eu/news/107-global-flood-monitoring-product-launch/); [GFM poster](https://emergency.copernicus.eu/downloads/GFM_poster.pdf)
- CEMS Rapid Mapping delivers "delineation" products as vector packages (ZIP) per activation, for example EMSR-numbered events — [CEMS EMSR498](https://emergency.copernicus.eu/mapping/list-of-components/EMSR498)

**Solar and other objects**
- Wherobots' segmentation tutorial uses solar farms as its example. Inference runs as SQL (`RS_Segment`) and output is vectorized with `vectorize_mosaic` / `RS_Polygonize` — [Wherobots segmentation tutorial](https://docs.wherobots.com/latest/tutorials/wherobotsai/wherobots-inference/segmentation); [RS_Polygonize docs](https://docs.wherobots.com/reference/wherobots-db/raster-data/pixel-functions/RS_Polygonize)

### Inferences
- In this domain, the polygonization volume for one product is on the order of 10^9 polygons per annual run (FTW, buildings). Products are increasingly refreshed every year (FTW), are multi-temporal (Google 2.5D 2016–2023), or are refreshed many times (Microsoft has done 25). "Re-vectorize the planet every year, or every time the model improves" is therefore a recurring workload, not a one-time cost.
- FTW is a clean target workload. Its polygonization input is exactly a thresholded binary or ternary mask at 10 m, polygons are tiny (median <0.3 ha, i.e. tens of pixels), there are billions of them, and output goes to GeoParquet in a partitioned layout. This closely matches a binary-bitmap, halo-tiled, rectilinear-ring engine with a georeferencing transform.
- Google 2.5D Temporal ships as a raster rather than vectors, which suggests that vectorizing every year at 0.5 m globally is costly or impractical. This is a candidate demonstration: GPU-vectorize the whole time series.

### Gaps
- No public runtime or cost figures were found for the polygonization stage of FTW, Microsoft or Google pipelines. FTW's paper explicitly omits them. The share of total pipeline time spent on vectorization remains unquantified.
- No sources were found for Varuna, EU IACS/CAP parcel pipelines, mining, or deforestation-alert polygonization (GLAD/RADD alerts, which are raster) in the searches run.

---

## Q2. What are geospatial foundation models producing, and how do their outputs get vectorized? Is vectorization a known bottleneck?

### Takeaway
Foundation models output either per-pixel masks (SAM/samgeo, Prithvi/Clay fine-tunes) or dense per-pixel embeddings (AlphaEarth: 64-D at 10 m). Vectorization is handled by generic CPU tools such as `rasterio.features.shapes` followed by shapely/geopandas regularization. Nobody frames it as a headline bottleneck, but every toolkit re-implements the same mask-to-shapes-to-regularize chain on CPU.

### Cited Findings
- samgeo (segment-geospatial) converts masks with `raster_to_vector`, `raster_to_geojson`, `raster_to_shp` and `raster_to_gpkg`, using rasterio (`rasterio.features.shapes`), GeoPandas and Shapely. It then post-processes with `regularize()`, `regularize_legacy()`, `simplify()` and `smooth_vector()`. The docs give no performance figures — [samgeo docs](https://samgeo.gishub.org/common/)
- Google AlphaEarth Foundations (arXiv 2507.22291, 2025) produces an "embedding field": an annual 64-dimensional vector per 10 m pixel from 2017 onward, hosted in Earth Engine as the Satellite Embedding dataset — [arXiv 2507.22291](https://arxiv.org/pdf/2507.22291); [ArxivIQ summary](https://arxiviq.substack.com/p/alphaearth-foundations-an-embedding)
  - Typical workflows aggregate embeddings within field polygons (zonal statistics). That means AEF consumes existing polygons rather than producing them — [Wherobots AEF zonal stats](https://wherobots.com/blog/alphaearth-embeddings-zonal-statistics/)
  - Its 64-D output is much more compact than Prithvi or Clay hidden states and is designed for classical ML (RF/GBT/MLP) — [Medium, Brian Curry](https://medium.com/@brian-curry-research/geospatial-foundation-models-for-agricultural-machine-learning-a3cd0169d256)
- Wherobots RasterFlow (Dec 2025) explicitly packages the full chain "build mosaics → run inference with PyTorch CV models → vectorize results" as a product. Inputs are Zarr/COG and output is an Iceberg table — [Wherobots RasterFlow](https://wherobots.com/blog/rasterflow-earth-observation-inference-engine/)
- Research on end-to-end polygon prediction aims to bypass mask-to-polygon conversion entirely:
  - Pix2Poly (arXiv 2412.07899): sequence prediction of polygonal building footprints — [arXiv 2412.07899](https://arxiv.org/pdf/2412.07899)
  - Global Collinearity-aware Polygonizer (arXiv 2505.01385) — [arXiv 2505.01385](https://arxiv.org/pdf/2505.01385)
  - The P³ dataset (Pixels, Points, Polygons for multimodal building vectorization, arXiv 2505.15379) — [arXiv 2505.15379](https://arxiv.org/pdf/2505.15379)
  - Topology-aware agricultural parcel vectorization (arXiv 2609.07520, 2026). It states that "segmentation-based methods mainly produce parcel masks or raster boundary cues and rely on heuristic raster-to-vector conversion". Its method rebuilds parcels from a shared planar graph and reports "zero measured intrusion ratio" — [arXiv 2609.07520](https://arxiv.org/abs/2609.07520)

### Inferences
- The research community sees naive raster-to-vector conversion as a quality problem (staircase edges, no shared boundaries, redundant vertices). It does not frame it as a speed problem, and is moving toward learned polygon outputs.
- A fast exact rectilinear polygonizer could still serve as the deterministic first stage feeding regularization or topology reconstruction. It could also be the scalable baseline that learned methods are benchmarked against.
- Exact rectilinear rings preserve pixel-exact topology: adjacent regions share edges exactly. That is a property generic pipelines then lose through per-polygon simplification.
- AlphaEarth and similar embedding fields create a new pattern: unsupervised clustering or classification over embeddings produces label rasters. Turning those into objects means polygonizing a whole continent at 10 m, potentially every year. That is a fresh, high-volume demand for raster-to-vector conversion.

### Gaps
- No explicit statements were found from the Prithvi (IBM/NASA), Clay, or SatMAE teams on how their downstream outputs are vectorized, or that vectorization is a bottleneck.
- No benchmark comparing the speed of mask-to-polygon conversion against model inference in a GeoFM pipeline was found.

---

## Q3. What tools exist and what are their limits? (GDAL, rasterio, PostGIS, Sedona, cuSpatial, etc.)

### Takeaway
Every mainstream polygonizer is CPU-based and largely single-threaded per raster. This includes GDAL Polygonize, rasterio.features.shapes, Sedona/Wherobots RS_Polygonize, xarray-spatial polygonize, ArcGIS Raster to Polygon and FME. NVIDIA RAPIDS cuSpatial has no raster-to-vector conversion. Documented pains are memory on large or complex rasters, needing everything in memory (rasterio), slow output writing, and no built-in tiling or merging.

### Cited Findings
- GDAL's own polygonize documentation: the algorithm tries to minimize memory, but "if the raster has many polygons or very large/complex polygons, the memory use for holding polygon enumerations and active polygon geometries may grow to be quite large" — [gdalraster polygonize doc](https://rdrr.io/cran/gdalraster/man/polygonize.html)
- rasterio issue #630, "Polygonize for larger-than-memory datasets": polygonize requires reading the full array into memory and an in-memory dataset before calling `gdal.Polygonize`, so larger-than-memory rasters fail. The issue is still an open enhancement with no PR — [rasterio #630](https://github.com/rasterio/rasterio/issues/630)
- GDAL issue #13169 (2025-10-03), "gdal raster polygonize: writing to temporary GeoPackage is slow": writing took 2 min 4 s to GeoPackage vs ~9.6 s to Shapefile. The issue says "Writing these features is much slower than generating them". Fix PR #13179 — [GDAL #13169](https://github.com/osgeo/gdal/issues/13169)
- GDAL issue #5000 reports a problem with the `-8` (8-connectedness) option of gdal_polygonize.py — [GDAL #5000](https://github.com/OSGeo/gdal/issues/5000)
- xarray-spatial added an experimental polygonize (PR #585). It is pure NumPy CPU with no GPU, Numba or Dask path, and is "simpler and has fewer memory reallocations, at the cost of requiring more RAM" than GDAL. No benchmarks were provided — [xarray-spatial PR #585](https://github.com/xarray-contrib/xarray-spatial/pull/585)
- GlacioHack/geoutils PR #962 worked on improving the performance of rasterization, polygonize and clip, showing ongoing CPU optimization work — [geoutils PR #962](https://github.com/GlacioHack/geoutils/pull/962)
- RAPIDS cuSpatial is a "vector-based, GPU accelerated GIS library" covering spatial computation, indexing, joins and trajectories. Its API mirrors GeoPandas and leaves unsupported operations to GeoPandas; no raster-to-polygon function is documented — [cuSpatial docs](https://docs.rapids.ai/api/cuspatial/stable/); [cuSpatial GitHub](https://github.com/rapidsai/cuspatial); [Hivenet guide](https://www.hivenet.com/post/gpu-geospatial-rapids-cuspatial-joins-pip)
- Apache Sedona (and SedonaDB) and WherobotsDB provide `RS_Polygonize`, which returns a polygon plus pixel value for each connected region — [SedonaDB RS_Polygonize](https://sedona.apache.org/sedonadb/latest/reference/sql/rs_polygonize/); [Wherobots RS_Polygonize](https://docs.wherobots.com/reference/wherobots-db/raster-data/pixel-functions/RS_Polygonize)
- Safe Software (FME) users filed an idea to "improve support for raster to vector" because production workflows were dropping out of FME into another tool "because the performance is so slow" — [FME Community idea](https://knowledge.safe.com/idea/52539/improve-support-for-raster-to-vector-rastertopolyg.html)
- ArcGIS Raster to Polygon offers a "simplify" option that smooths polygons into shapes with a minimum number of segments while staying close to cell edges. Esri also ships a separate "Regularize Adjacent Building Footprint" tool to fix stair-stepped building outputs — [ArcGIS Raster to Polygon](https://pro.arcgis.com/en/pro-app/3.4/tool-reference/conversion/raster-to-polygon.htm); [Regularize Adjacent Building Footprint](https://pro.arcgis.com/en/pro-app/3.4/tool-reference/3d-analyst/regularize-adjacent-building-footprint.htm)
- In R, `rasterToPolygons` is far slower than gdal_polygonize.py. Users report GDAL "doesn't do well with ultra-large rasters", and the workaround of splitting into tiles and processing in parallel brings national high-resolution rasters "down to tens-of-minutes instead of days" — [r-sig-geo thread](https://r-mailing-lists.thecoatlessprofessor.com/lists/r-sig-geo/msg/msg-9f8c2118727f/raw)

### Inferences
- A GPU polygonizer that emits GeoArrow/GeoParquet-compatible buffers directly would fill a real hole in the RAPIDS/cuSpatial stack. The input would be a cuPy/torch mask already on the GPU after inference. Today the mask has to be copied to host and handed to GDAL.
- The GDAL #13169 finding, where writing is 10x slower than generation, warns that the end-to-end win also depends on output serialization (columnar formats, avoiding per-feature OGR writes), not just the polygon-tracing kernel.
- An exact 2x2-classification scanline approach is an alternative to the many-polygon memory growth that GDAL's docs acknowledge.

### Gaps
- No GIS StackExchange threads were retrieved directly in this search pass. The complaints above come from GitHub, mailing lists and vendor forums.
- No published benchmark of GDAL polygonize throughput (pixels/s or polygons/s) on modern hardware was found.
- PostGIS ST_DumpAsPolygons, WhiteboxTools, tippecanoe and Martin were not researched (budget). tippecanoe and Martin consume vectors and do not polygonize rasters.

---

## Q4. Tile seam / halo problems: how do people merge polygons across tiles at continent scale? What is the state of the art?

### Takeaway
The common practice is ad hoc: tile the raster, polygonize tiles in parallel, append, then dissolve/unary-union or heuristically prune duplicates in overlap zones. On the model side, overlapping windows with weighted blending suppress inference seams, but the polygon-level seam problem remains unsolved in open tooling. Research is moving toward stitching intermediate primitives first and polygonizing once.

### Cited Findings
- A Mapbox gdal-polygonize-test user polygonized a 110,000 × 260,000 raster in about 20 minutes. They split it into a 4×4 grid with gdal_translate, used multiprocessing.Pool, and appended with ogr2ogr. The open question "dissolve or union adjacent polygons after merging?" was never answered in the thread, which was archived in 2021 — [mapbox/gdal-polygonize-test #2](https://github.com/mapbox/gdal-polygonize-test/issues/2)
- In building-vectorization research, the standard practice is to patch large aerial images with 10% overlap and merge redundant polygons in overlap regions with a unary union. Traditional approaches "extract polygons from individual patches, lift them to tile coordinates, and heuristically merge or prune polygons in overlap regions". Newer ones "stitch intermediate primitives in tile coordinates and perform topology-aware polygonization once on the full tile" — [P³ / Pix2Poly-related search summaries, arXiv 2505.15379](https://arxiv.org/pdf/2505.15379); [arXiv 2412.07899](https://arxiv.org/pdf/2412.07899)
  - Caveat: these statements come from search-snippet summaries of these papers; full text was not fetched.
- A US patent, "Method and system for merging of polygons in adjacent tiles" (US 10157485), covers merging polygons that share overlapping edges at a tile boundary — [USPTO 10157485](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10157485)
- FTW handles inference seams with a 25% overlap sliding window and Gaussian apodization before polygonization. Polygon-level seam handling at country partitions is not documented — [arXiv 2602.08131](https://arxiv.org/html/2602.08131v1); [FTW vectors llms.txt](https://data.source.coop/ftw/global-data/predictions/vectors/llms.txt)
- Microsoft partitions output by Bing quadkey level 9. How cross-quadkey buildings are handled is not documented in the sources found — [marksblogg](https://tech.marksblogg.com/microsofts-global-ml-building-footprints.html)

### Inferences
- **GAP / opportunity:** no open tool offers exact, deterministic cross-tile polygon merging for polygonized rasters. The halo plus pixel-exact rectilinear rings design makes an exact approach possible. Ring segments that hit the tile boundary land on exact pixel-grid coordinates, so stitching becomes integer edge matching instead of floating-point unary_union. That is cheaper and topologically exact, with no slivers or near-coincident vertices. This looks like the strongest technical differentiator for continent-scale work.
- Gaussian-weighted inference blending (FTW) means the mask itself is seamless in the global Zarr. The remaining seam is purely in polygonization and partitioning, which is exactly the engine's layer.

### Gaps
- No published state-of-the-art benchmark for distributed polygon stitching (accuracy or runtime) was found.
- The approach Overture uses to conflate ML footprints across sources and tiles was not researched in depth.

---

## Q5. Real-time / on-board satellite processing: do satellites emit vectors to save downlink?

### Takeaway
Yes, though so far for detections (points/boxes), not delineated polygons. Planet's Pelican-4 ran Jetson Orin object detection in orbit in March 2026 and produced geo-rectified GeoTIFF and GeoJSON outputs. ESA Φsat-2 runs six onboard apps, including cloud, vessel, fire and street detection. Onboard emission of full delineation polygons (flood/burn extents) is an open frontier.

### Cited Findings
- On 2026-03-25, Planet's Pelican-4 imaged an airport over Alice Springs and ran an airplane-detection model on its onboard NVIDIA Jetson Orin, producing "geo-rectified GeoTIFF and GeoJSON outputs" at about 80% initial accuracy on raw imagery — [Stocktitan / Planet PR](https://www.stocktitan.net/news/PL/planet-successfully-runs-ai-in-msc7rawzxhoq.html); [Orbital Today 2026-04-08](https://orbitaltoday.com/2026/04/08/planets-pelican-4-runs-nvidia-powered-ai-in-space-for-the-first-time/)
  - Detection took about 0.5 s, cutting latency "from hours to minutes" and reducing downlink needs — [blockchain.news](https://blockchain.news/ainews/planet-labs-deploys-onboard-ai-for-0-5s-detections); [letsdatascience](https://letsdatascience.com/news/planet-labs-deploys-edge-ai-onboard-satellites-for-real-time-1187251a)
  - Caveat: these are secondary outlets relaying a Planet press release.
- Pelican fleet: Pelican-1 launched in 2024 and Pelicans 2–6 in 2025. Each carries a Jetson for on-orbit compute and delivers 50 cm-class imagery in 6 multispectral bands. Planet's 2026 roadmap includes "rapid downlink, and select on-orbit processing capabilities" — [Planet "Pelican 2025 in Review"](https://www.planet.com/pulse/pelican-2025-in-review/)
- ESA Φsat-1 (2020) used an Intel Movidius Myriad 2 to discard cloudy images onboard, saving about 30% of bandwidth — [Intel newsroom](https://www.intel.co.uk/content/www/us/en/newsroom/news/first-satellite-ai.html); [eeNews Europe](https://www.eenewseurope.com/en/first-satellite-with-ai-on-board-uses-intel-chip)
- ESA Φsat-2 (Open Cosmos, with Ubotica CogniSAT) runs six onboard AI apps: deep compression, cloud detection (KP Labs), street detection, vessel detection (CEiiA), fire detection and marine anomaly detection. More than 100 GB was acquired during commissioning — [ESA](https://esa.int/Applications/Observing_the_Earth/Phsat-2/New_satellite_demonstrates_the_power_of_AI_for_Earth_observation); [Open Cosmos](https://open-cosmos.com/news/phisat-2-launch); [CEOS WGISS-61 slides](https://ceos.org/document_management/Working_Groups/WGISS/Meetings/WGISS-61/7.4_Phisat-2_ESA%20V3.pdf)

### Inferences
- Onboard compute is now Jetson Orin-class, which is CUDA-capable. A compact CUDA polygonizer that turns onboard segmentation masks (fire fronts, flood extent, street networks, ice edges) into GeoJSON polygons fits naturally. Polygons are far smaller than masks, especially as rectilinear rings with run-length-style vertex counts.
- The current onboard outputs (Planet aircraft boxes, Φsat-2 vessel detection) are point or box detections. The "segment → polygon → downlink a few KB of GeoJSON" pattern for extents was not found deployed. It is a plausible gap but unconfirmed.

### Gaps
- No source confirmed any satellite downlinking delineated polygons (as opposed to detections or masks).
- No Satellogic-specific onboard AI details were found in the searches.
- Downlink bandwidth numbers (Mbps, cost per GB) were not found.

---

## Q6. Regulatory or market drivers creating new demand

### Takeaway
EUDR is the clearest regulatory driver that mandates polygons. Plots over 4 ha need perimeter polygons in WGS-84 GeoJSON, with 6 decimal places, uploaded to TRACES. Large operators must comply from 30 Dec 2026, though the date has been slipping. Other drivers (3D city/LoD1 models, flood monitoring, carbon MRV, insurance) consume polygon datasets, but this pass found less hard data on them.

### Cited Findings
- EUDR geolocation: plots over 4 ha (except cattle) must be described by polygons with multiple perimeter points. Plots under 4 ha may use a single point. Coordinates are decimal degrees with 6 decimals on WGS-84 (EPSG:4326), uploaded as GeoJSON to the EU TRACES information system — [Coolset](https://www.coolset.com/academy/eudr-geolocation-requirements-how-to-collect-and-validate-gps-polygon-data-for-your-dds); [Fairtrade geolocation guide](https://www.fairtrade.net/content/dam/fairtrade/global/fairthrive/Geolocation%20data%20requirements%20EUDR.pdf); [EFI Insight, Dec 2025](https://efi.int/sites/default/files/files/flegtredd/TF/Documents/Insight%20-%20Geolocation%20and%20the%20EUDR%20Dec25.pdf)
- EUDR timeline: per the Coolset page, obligations apply from 30 Dec 2026 for medium and large operators and 30 Jun 2027 for micro and small enterprises. In Nov 2025, Parliament and Council backed changes including a 12-month enforcement delay — [Coolset](https://www.coolset.com/academy/eudr-geolocation-requirements-how-to-collect-and-validate-gps-polygon-data-for-your-dds)
  - Note: the exact current legal status should be checked against EUR-Lex. Search summaries differed on whether the delay text was final.
- Satellite imagery is an accepted verification method for EUDR — [EUDR guide summary](https://greensutra.in/eudr-guide/)
- 3D city and LoD1 models: GlobalBuildingAtlas offers 2.68B LoD1 building instances, positioned for urbanization monitoring and SDG analysis — [arXiv 2506.04106](https://arxiv.org/abs/2506.04106)
- Emergency response: Copernicus GFM is a fully automated global SAR flood service designed to speed up flood maps and CEMS Rapid Mapping activations — [GFM launch](https://global-flood.emergency.copernicus.eu/news/107-global-flood-monitoring-product-launch/)

### Inferences
- EUDR creates a two-sided demand. Supply-chain plot polygons, often hand-drawn or GPS-traced, must be checked against raster deforestation alerts and forest masks, which in practice means polygonizing deforestation and forest-cover rasters or doing raster–polygon overlay at national scale.
- Pixel-exact rectilinear rings at 10 m may not suit compliance submission directly. Simplification is likely needed, and EUDR polygons are farmer plots, not raster-derived. The stronger fit is the verification side: polygonized forest-loss and forest masks for overlay and audit.
- Parametric insurance (flood or burn footprint triggers) and carbon MRV (forest-stand or project-boundary change) are plausible consumers of fast, repeated polygonization, but this pass did not source them.

### Gaps
- No sourced numbers were found for parametric insurance, carbon MRV, mining or solar-panel mapping volumes or pipelines.
- No primary EUR-Lex text was fetched for the EUDR amendment status.
- Commercial field-boundary vendors (Varuna, EU IACS/LPIS national agencies, Planet Field Boundaries, EOSDA) were not researched.
