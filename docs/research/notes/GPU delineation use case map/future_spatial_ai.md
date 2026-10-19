# Future Trajectory of Spatial Datasets and Spatial AI (2025–2030), and Where a GPU Raster→Vector Engine Becomes Strategic

Scope note: research done 2026-10-07. Covers trends and future directions only. Other researchers cover today's EO, robotics, biomedical, EDA and ML geometry. Each finding carries an inline source. Inferences and speculation are labeled. "The engine" means the user's CUDA scanline engine: binary bitmap with halos in, closed polygon rings out, with a pluggable coordinate transform.

---

## 1. Spatial intelligence and world models (World Labs, Niantic Spatial, Genie 3, NVIDIA Cosmos, Apple/visionOS): what do they produce, and does anything need boundary extraction at scale?

### Takeaway
Money and attention have moved hard toward "spatial intelligence" and world models: a $1B round for World Labs in Feb 2026, a $250M spin-out for Niantic Spatial, and Cosmos trained on 20M hours of video. Their native outputs are video frames, 3D Gaussian splats and coarse meshes or planes, not clean vector geometry. Turning them into usable structure (floorplans, walkable or occupied footprints, semantic object outlines) is an unsolved downstream step. A 2D occupancy-slice → polygon step is a plausible fit, but no source describes it as a recognized bottleneck yet.

### Cited Findings
- World Labs raised **$1 billion in new funding (announced Feb 18, 2026)** from AMD, Autodesk, Emerson Collective, Fidelity, NVIDIA, Sea and others. Its product **Marble** creates "spatially cohesive, high-fidelity, and persistent 3D worlds from images, video, or text". The stated targets are storytelling, creativity, robotics and scientific discovery. — [World Labs blog](https://www.worldlabs.ai/blog/funding-2026)
- World Labs launched a **World API**, a public interface for generating navigable 3D worlds with Marble from text, images, panoramas, multi-view inputs and video. — [World Labs: Announcing the World API](https://www.worldlabs.ai/blog/announcing-the-world-api)
- Marble outputs full-resolution 3D scenes with **Gaussian splat export**. Tencent's HY-World-2.0 also generates 3D Gaussian splat output. (Secondary/aggregator source.) — [Tech-Insider comparison 2026](https://tech-insider.org/genie-3-vs-marble-vs-nvidia-cosmos-world-models-2026/); [Spheron GPU infra for world models](https://www.spheron.network/blog/gpu-infrastructure-world-models-2026/)
- **Niantic Spatial** was spun out of Niantic in **May 2025** after the games were sold to Scopely. It builds "geospatial world models" (3D reconstruction, visual positioning, semantic understanding) on a **Large Geospatial Model (LGM)**. In 2026 it partnered with **Spexi Geospatial** (drone imagery) for city-scale 3D models, piloted with **ExxonMobil**, and worked with Flexion AI and **NVIDIA** on robotics simulation training. — [Wikipedia: Niantic Spatial](https://en.wikipedia.org/wiki/Niantic_Spatial)
- Niantic Spatial launched with **$250M** ($200M from Niantic plus a $50M strategic investment from Scopely). It kept **Scaniverse** (mobile Gaussian-splat capture), its VPS, and an LGM "built on **30B+ posed images**". — [Radiance Fields](https://radiancefields.com/niantic-spatial-bets-big-on-scaniverse-and-gaussian-splatting)
- In 2026 Niantic Spatial released a revamped business Scaniverse and a new **global VPS**, and plans to add **semantic understanding** so that AI can reason about objects and environments. — [GeekWire 2026](https://www.geekwire.com/2026/from-pokemon-go-to-physical-ai-niantic-spatial-unveils-its-global-3d-mapping-platform/) (search snippet only; the page returned 403)
- **Genie 3** (Google DeepMind, announced Aug 2025) generates interactive environments at **720p / 24 fps** and stays consistent for **several minutes** through "World Memory". It is closed, with a limited preview through Project Genie for Google AI Ultra subscribers. — [YourStory](https://yourstory.com/ai-story/deepmind-genie-3-ai-world-model); [Tech-Insider](https://tech-insider.org/genie-3-vs-marble-vs-nvidia-cosmos-world-models-2026/)
- **NVIDIA Cosmos** world foundation models (CES, Jan 2025) were trained on **9,000 trillion tokens from 20 million hours of video**. NVIDIA says Blackwell processes 20M hours of video in two weeks, against 40 days on Hopper and about 3 years on CPUs. Open variants exist (Cosmos-Predict 7B/14B). — [The Decoder](https://the-decoder.com/nvidias-new-cosmos-world-models-aim-to-understand-physics-through-video/); [NVIDIA press release](https://investor.nvidia.com/news/press-release-details/2025/NVIDIA-Launches-Cosmos-World-Foundation-Model-Platform-to-Accelerate-Physical-AI-Development/default.aspx)
- Extracting surfaces from Gaussians is a known weak point. Converting Gaussian opacity directly into an **occupancy field for mesh extraction** gives "unsatisfactory surfaces with visible holes" because generated Gaussians are sparse. Workarounds use block-wise local density queries. — [DreamGaussian, arXiv 2309.16653](https://arxiv.org/pdf/2309.16653)
- Apple: **RoomPlan** (iPhone/iPad LiDAR) produces parametric floor plans with walls, windows, openings, doors and furniture. RoomPlan is **not available on visionOS**, where apps only get ARKit's `PlaneDetectionProvider` and `SceneReconstructionProvider` (a "coarse" mesh). — [Apple Dev Forums](https://developer.apple.com/forums/thread/746957); [WWDC22 RoomPlan](https://developer.apple.com/videos/play/wwdc2022/10127/?time=339)

### Inferences
- (Inference) All of these systems emit dense or implicit representations (splats, meshes, video, voxels). Robotics, AEC, real estate and games need **2D plans and footprints**: navigable floor area, wall outlines, object footprints for collision and placement, and room polygons for BIM/IFC. A standard way to get there is: rasterize a height slice or top-down occupancy of the splat/mesh, then binary bitmap → polygon rings → transform to world coordinates. This maps directly onto the engine's input and output contract, including the pluggable transform from voxel/pixel to world or ECEF coordinates.
- (Inference) Niantic's plan for semantic understanding on an LGM built from 30B+ images, plus city-scale models from Spexi drones, points to **per-class semantic occupancy grids at city scale** that will need vectorizing into GIS/BIM-compatible footprints. That is a tiling problem where halo buffers matter.
- (Speculation) If world models such as Genie 3 or Cosmos are used to generate synthetic training worlds for robots, then per-frame or per-world "map extraction" (free space and obstacle polygons) becomes a high-throughput batch job on GPUs that already hold the data. A CPU round-trip for vectorization would be the odd step out.

### Gaps
- No source found that names 2D boundary or polygon extraction from world-model outputs as a recognized bottleneck. This is a predicted gap, not a documented one.
- Could not confirm Marble's exact export formats (mesh, collider mesh) from a primary source. The World Labs funding post does not list them.
- No primary data found on Android XR or Meta scene-understanding representations in 2026.
- Valuations for World Labs and Niantic Spatial in 2026 were not disclosed in the sources reviewed.

---

## 2. Embedding fields ("embedding as the new raster"): how do users get objects or regions from embedding rasters, and what are the pain points?

### Takeaway
Global annual 10 m embedding rasters now exist as free products: Google's AlphaEarth/Satellite Embedding (64-D, 2017–2025) and TESSERA (128-D). The standard workflow to get objects is: dot-product similarity, then a threshold, then a **binary mask**, then **vectorize**, then centroids or polygons. Google's own tutorial warns that the vectorization step is "memory and compute-intensive" and causes timeouts. This is the clearest documented match for the engine's exact input type: a thresholded binary bitmap, tiled, at planetary scale.

### Cited Findings
- **AlphaEarth Foundations Satellite Embedding dataset**: global, annual, **10 m**, **64 channels**, years **2017–2025**. It is distributed as Cloud-Optimized GeoTIFFs of **8192×8192 px** tiles, and each pixel vector is normalized to unit length. It fuses Sentinel-2, Landsat, radar, 3D/lidar, elevation, climate, gravity and text. — [Source Cooperative: AEF](https://source.coop/tge-labs/aef); [Google Earth Medium: AEF on GCS](https://medium.com/google-earth/alphaearth-foundations-satellite-embeddings-now-available-on-google-cloud-storage-f9ab0f7252d6)
- Google's **similarity-search tutorial** (Earth Engine) works like this: `meanDistance.gt(0.90)` gives a binary mask. Then `similarPixels.selfMask().reduceToVectors({scale: 10, eightConnected: false, maxPixels: 1e10})` gives polygons, and centroids of those polygons are the "detected objects" (grain silos in the example). The tutorial says verbatim: "**Vectorization is a memory and compute-intensive operation in Earth Engine**". For larger regions it recommends exporting as assets to avoid "**Computation Timed Out**" or "**Tile Error**". The results contain false positives that need threshold tuning. — [Earth Engine tutorial: Similarity Search](https://developers.google.com/earth-engine/tutorials/community/satellite-embedding-05-similarity-search)
- **TESSERA** is 128-D at 10 m, 2017–2025, global. The **Tessera Embeddings Explorer (TEE)** lets users draw polygon labels, whose interiors are rasterized at 10 m. That is the reverse direction, vector → raster for training labels. — [TESSERA GitHub (TEE)](https://github.com/ucam-eo/TEE); [Earth Embeddings arXiv 2608.03410](https://arxiv.org/html/2608.03410)
- On **Sept 25, 2026** Copernicus Data Space published an openEO community example for generating TESSERA pixel embeddings from Sentinel-1/2. Embedding generation is becoming a hosted service. — [Copernicus Data Space news](https://dataspace.copernicus.eu/news/2026-9-25-new-openeo-community-example-generate-tessera-pixel-embeddings-sentinel-1-and)
- "**Earth Embeddings as Products**" (Fang, Stewart, Corley, Zhu, Azizpour; **Jan 19, 2026**) catalogs Clay (5.12 km, 768-D), Major TOM (2.14–3.56 km, 2048-D), Earth Index (320 m, 384-D), Copernicus-Embed (0.25°, 768-D), Presto (10 m, Togo only), TESSERA (10 m, 128-D) and Google Satellite Embedding (10 m, 64-D). The authors call the landscape "**fragmented**": formats differ (GeoParquet for patches, GeoTIFF for pixels), there is "no unified methodology to load or evaluate them", and there is insufficient **explainability and uncertainty quantification** for downstream use. — [arXiv 2601.13134](https://arxiv.org/html/2601.13134v1)
- Downstream research uses of AEF/TESSERA are growing fast in 2025–2026: SAR flood segmentation with AEF land-cover priors ([arXiv 2606.29134](https://arxiv.org/pdf/2606.29134)), Local Climate Zone mapping in 5 Swiss cities ([arXiv 2606.20034](https://arxiv.org/pdf/2606.20034)), streetscape perception ([arXiv 2608.16310](https://arxiv.org/pdf/2608.16310)), burned-area mapping ([arXiv 2509.07852](https://arxiv.org/pdf/2509.07852)), and scalable geospatial data generation ([arXiv 2508.11739](https://arxiv.org/pdf/2508.11739)). Most of these produce segmentation **masks**, which become polygons when delivered as GIS layers.
- Google Earth AI's **Population Dynamics** embeddings are **updated monthly** (17 countries), which shows that embedding fields are moving toward higher temporal cadence. — [Google Research blog, Oct 23 2025](https://research.google/blog/google-earth-ai-unlocking-geospatial-insights-with-foundation-models-and-cross-modal-reasoning/)

### Inferences
- (Inference, back-of-envelope) Global land is about 1.49×10⁸ km², which is about **1.5×10¹² pixels per year at 10 m**. Across 9 annual layers (2017–2025) that is about 1.3×10¹³ pixel-years. Every new similarity or classification query produces a new binary bitmap of up to that size. "Query-time vectorization" is therefore a recurring, interactive workload, not a one-off ETL job.
- (Inference) The "embedding as the new raster" pattern turns object detection into **thresholding**. The expensive step is no longer the model but **mask → polygons**, followed by area/centroid/count statistics. A GPU engine that vectorizes thresholded tiles in a single pass, with halos so that objects crossing the 8192×8192 COG tile boundaries close correctly, targets the exact pain point Google documents.
- (Speculation) A strong product idea: **"interactive embedding search → instant polygons"**. A user clicks examples, the similarity raster is computed on GPU, the threshold slider re-vectorizes in milliseconds, and polygons stream out as GeoArrow/GeoParquet. Earth Engine's own guidance today is "export as asset and wait".

### Gaps
- No published benchmark found that measures vectorization time against similarity-compute time for AEF/TESSERA workflows.
- No source found on total byte volume of the AEF dataset or of TESSERA.
- No third-party tool found that specifically offers fast embedding-mask → polygon conversion. The absence is notable, but not proof that none exists.

---

## 3. Spatial LLM agents and reasoning: what geometry operations do they need?

### Takeaway
Through 2025–2026 every major GIS vendor (Google, Esri, CARTO, Wherobots) shipped agentic or MCP-connected spatial AI. The agents call tools such as Earth Engine, BigQuery, spatial joins and model inference. When an agent asks for something like "find X in imagery and intersect with Y", it needs a **fast raster-model-output → vector conversion** so it can join, count and intersect. Interactive agent loops make latency matter in a way batch ETL never did.

### Cited Findings
- **Google Earth AI** (expanded **Oct 23, 2025**): a Gemini-powered geospatial reasoning agent breaks queries into multi-step plans and calls Earth Engine, BigQuery public datasets (official boundaries), Data Commons, satellite-imagery models (open-vocabulary detector **RS-OWL-ViT-v2**) and on-the-fly model training. On a Q&A benchmark it scores **0.82 accuracy against 0.50 for Gemini 2.5 Pro and 0.39 for 2.5 Flash** (+124% on analytical tasks). Enterprise users include **CARTO**, Public Storage and Visiona. Outputs include vector intersections such as wind-zone × county polygons and infrastructure identified in imagery. — [Google Research blog](https://research.google/blog/google-earth-ai-unlocking-geospatial-insights-with-foundation-models-and-cross-modal-reasoning/)
- Google's Hurricane Helene demo: the agent identified Florida counties with given populations in the highest-wind path and computed vulnerability scores. — [R&D World](https://www.rdworldonline.com/google-earth-ai-shifts-from-map-layers-to-decision-support-with-geospatial-reasoning-agent/)
- **Esri**: the **2026 Esri UC** (San Diego, 18,000+ in-person attendees) framed ArcGIS as "**an AI system with a map inside it**" ([Cercana Systems briefing, Jul 2026](https://cercanasystems.com/2026/07/supplemental-briefing-the-2026-esri-user-conference/)). ArcGIS added **MCP support** so third-party agents can call ArcGIS tools ([Esri ArcNews](https://www.esri.com/about/newsroom/arcnews/the-next-era-of-ai-and-arcgis)). ArcGIS Online (June 2026) shipped agentic mapping apps in beta, and the Maps SDK for JS added AI components (beta). — [What's new in ArcGIS Online, June 2026](https://www.esri.com/arcgis-blog/products/arcgis-online/announcements/whats-new-in-arcgis-online-june-2026); [Esri JS AI components](https://www.esri.com/arcgis-blog/products/js-api-arcgis/developers/introducing-ai-components-beta-in-the-js-maps-sdk)
- **CARTO** published open **Agent Skills** that teach AI agents to operate CARTO through CLI and MCP Server, including an Esri→CARTO migration skill. — [CARTO blog](https://carto.com/blog/migrate-from-esri-to-carto-ai-agent/)
- **Wherobots** offers MCP and CLI "Spatial AI Coding Tools" and is listed in Anthropic's Connectors Directory. — [Wherobots AI](https://wherobots.com/products/ai/)
- Agent benchmarks: **GeoAgentBench** (2026) has 117 atomic GIS tools and 53 spatial-analysis tasks across 6 domains, and finds Plan-and-React works best ([arXiv 2604.13888](https://www.emergentmind.com/papers/2604.13888)). **GeoBenchX** tests 8 commercial LLMs with 23 geospatial functions ([arXiv 2503.18129](https://arxiv.org/html/2503.18129v2)). **GeoNatureAgent** has 93 tasks across 18 categories ([emergentmind 2606.12821](https://www.emergentmind.com/papers/2606.12821)).

### Inferences
- (Inference) Agents compose detector/segmenter output (raster masks) with authoritative vectors (counties, parcels, roads). Each such query needs a mask → polygon step before the join. In a chat loop that step should take seconds, not minutes. A GPU vectorizer exposed as an **MCP tool** ("polygonize this mask/tile set → GeoArrow") is a natural, small and composable agent primitive.
- (Inference) Agent benchmark tool inventories (117 tools in GeoAgentBench, 23 in GeoBenchX) include classic GIS ops. Raster-to-vector ("polygonize") is a standard GIS tool category, so agents will call it repeatedly, often on model outputs.
- (Speculation) As agents produce more "derived maps on demand", **provenance** (which model, threshold and tile produced this polygon) becomes a first-class need. An engine that emits per-ring metadata (source tile ID, threshold, pixel count) would differentiate.

### Gaps
- No source found that names vectorization latency as a bottleneck in agent loops. This is inferred from the GEE timeout warning in Q2 plus the interactivity requirement.
- No public tool-call statistics found from Earth AI or the ArcGIS agents (for example, how often polygonize is invoked).

---

## 4. Open spatial data growth (Overture, OSM AI mapping, national digital twins, cloud-native formats)

### Takeaway
Open vector maps are now dominated by **ML-derived geometry**. About 63% of Overture's 2.53B buildings come from Google or Microsoft ML pipelines. Formats are converging on native Parquet/Iceberg geometry types, and Overture's members position open maps explicitly as "grounding for AI". Releases are monthly, so the pipeline that turns AI rasters into conflated, ID-stable vectors runs continuously.

### Cited Findings
- **Overture 2026-09-23 release**: **2,533,842,612 buildings** (+0.17% from the prior release) and 4,486,107 building parts (+3.39%). Buildings by source: **Google Open Buildings 898,130,959**, OSM 706,818,187, **Microsoft ML Buildings 689,237,083**, Zenodo 212,222,368, Esri Community Maps 15,119,688, IGN España 12,298,053. The buildings theme is released **monthly** as GeoParquet on AWS and Azure. — [Overture 2026-09-23 release notes](https://docs.overturemaps.org/blog/2026/09/23/release-notes/)
- Overture reached **50 members (July 13, 2026)**, nearly double the 2024 count. New members include **Grab, Uber, Samsara, Maxar/Vantor, Regrid, Wherobots** and HOT. The framing is that "AI can't reason about the physical world without an open, verifiable foundation beneath it" (Will Mortenson, Executive Director). **GERS** gives stable IDs so that AI can anchor to entities. Microsoft reported "several percentage point gains in address accuracy" and "significant building coverage improvements" after adopting Overture layers. — [Overture announcement](https://overturemaps.org/announcements/2026/overture-maps-foundation-reaches-50-members-as-industry-converges-on-open-data-to-ground-ai/)
- **Google Open Buildings 2.5D Temporal**: annual 2016–2023 rasters of building **presence, fractional count and height**, at 4 m effective resolution (served at 0.5 m), covering about **58M km²** (Africa, S/SE Asia, LatAm and the Caribbean), derived from Sentinel-2. It is distributed as **rasters**, not polygons. — [Earth Engine catalog](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_Research_open-buildings-temporal_v1?hl=en); [Google Open Buildings temporal](https://sites.research.google/gr/open-buildings/temporal)
- **OSM AI-assisted mapping**: Rapid lets mappers accept AI-predicted features and integrates Esri open data. Its main developer, Bryan Housel, was laid off by Meta on **Jan 28, 2025**, and development was "mostly inactive from March to December 2025" before resuming. — [OSM Wiki: Rapid](https://wiki.OpenStreetMap.org/wiki/RapiD)
- **Cloud-native formats**: native **GEOMETRY/GEOGRAPHY types landed in Apache Parquet core in March 2025**. **Iceberg v3** includes them. "GeoParquet 2.0" means using these native types. The CNG community hopes 2026 is the year major datasets adopt them. — [Apache Parquet blog, Feb 13 2026](https://parquet.apache.org/blog/2026/02/13/native-geospatial-types-in-apache-parquet/); [CNG: GeoParquet 2.0 going native](https://cloudnativegeo.org/blog/2025/02/geoparquet-2.0-going-native/); [CNG: time of transition](https://cloudnativegeo.org/blog/2025/10/geoparquet-parquet-geospatial-types-a-time-of-transition/); [CARTO: Iceberg v3 geometry one year later](https://carto.com/blog/iceberg-v3-geometry-one-year-later/)
- **National/city digital twins**: the UK National Digital Twin Programme (since 2018) has a new **£37.6M National Digital Twin Centre in Belfast** ([Computing](https://computing.co.uk/news/4203936/belfast-spearhead-uks-digital-revolution-gbp37m-digital-twin-centre)). The UK **NUAR** (National Underground Asset Register) launched its first phase ([PBC Today](https://www.pbctoday.co.uk/news/digital-construction-news/bim-news/national-digital-twin-2/80648/)). Singapore began testing a **port digital twin with 1,000+ sensors** in H2 2025, building on Virtual Singapore ([OECD OPSI](https://oecd-opsi.org/innovations/virtual-twin-singapore/); search snippet).

### Inferences
- (Calculation) ML-derived sources (Google 898M + Microsoft 689M) are about **1.59B of 2.53B buildings, roughly 63%** of Overture buildings. The world's open building layer is mostly the output of segmentation → vectorization pipelines, and it is refreshed monthly.
- (Inference) Google's 2.5D Temporal product ships buildings as **rasters** over 58M km² × 8 years. Anyone who wants per-year footprints or change polygons must vectorize billions of pixels per year. That is a direct, existing, open workload for the engine (presence band → threshold → rings; transform to lat/lon).
- (Inference) Native Parquet/Iceberg geometry plus GeoArrow means the natural output of a GPU vectorizer is **GeoArrow buffers in GPU memory → Parquet**, with no WKT/Shapefile serialization. Ring output that maps directly onto GeoArrow's offset+coordinate layout would plug into SedonaDB, DuckDB and Iceberg lakehouses.

### Gaps
- No projection found for total open vector data volume in 2030. Overture's own announcement gives no volume stats.
- No primary 2026 data found on the Netherlands' or Singapore's national twin data volumes, or on how much of their content is AI-derived.

---

## 5. Market data and investment

### Takeaway
Forecasts put geospatial intelligence at about $63B by 2030 (11% CAGR) and digital twins at $150–196B by 2030 (40%+ CAGR, with wide divergence between vendors). The strongest real signal is capital flowing to spatial-AI and world-model companies (World Labs $1B, Niantic Spatial $250M) and NVIDIA's partnerships across the stack.

### Cited Findings
- **Geospatial Intelligence market**: **$37.13B (2025) → $62.88B (2030), CAGR 11.1%**. "GeoAI and ML platforms" are the fastest-growing software segment. North America is the largest region. Players named include Google, Esri, Hexagon, Trimble, Planet, Bentley and TomTom. — [MarketsandMarkets via GlobeNewswire, Aug 19 2026](https://www.globenewswire.com/news-release/2026/08/19/3347778/0/en/geospatial-intelligence-market-surges-to-62-88-billion-at-a-cagr-11-1-by-2030-report-by-marketsandmarkets.html)
- **Geospatial AI market**: $38.0B (2023) → $64.6B (2029). This is a low-quality press-release aggregator and should be treated with caution. — [OpenPR](https://www.openpr.com/news/4624919/global-geospatial-artificial-intelligence-market-size-to-reach)
- **Digital twin market**: **$21.14B (2025) → $149.81B (2030), CAGR 47.9%** ([MarketsandMarkets](https://www.marketsandmarkets.com/ResearchInsight/digital-twin-market-size.asp)). This conflicts with **$183B by 2030, 41.6% CAGR** ([Meticulous Research](https://www.meticulousresearch.com/pressrelease/218/digital-twin-market-2030)) and $113B–$196B in other reports ([search aggregate](https://jotup.co/node/2007279)). Treat these as order-of-magnitude figures only.
- **Funding**: World Labs **$1B** (Feb 2026; investors include NVIDIA, AMD, Autodesk) ([World Labs](https://www.worldlabs.ai/blog/funding-2026)). Niantic Spatial **$250M** at spin-out in 2025 ([Radiance Fields](https://radiancefields.com/niantic-spatial-bets-big-on-scaniverse-and-gaussian-splatting)).
- **NVIDIA as a cross-stack investor or partner**: World Labs investor; Niantic Spatial robotics-sim partner ([Wikipedia](https://en.wikipedia.org/wiki/Niantic_Spatial)); Wherobots RasterFlow on RTX PRO 6000 Blackwell (Sep 30, 2026) ([Wherobots](https://wherobots.com/blog/wherobots-nvidia-geospatial-ai-physical-world/)); Planet Pelicans carry **NVIDIA Jetson** for on-orbit compute ([Planet Pelican](https://www.planet.com/constellations/pelican)).

### Inferences
- (Inference) NVIDIA is building a GPU-native geospatial and physical-AI ecosystem, yet **discontinued cuSpatial** in June 2025 (see Q7). That leaves a hole in the open GPU geometry-primitive layer just as demand rises. A focused, high-performance raster→vector primitive could be adopted into this ecosystem through Sedona/Wherobots, NVIDIA Earth-2/Omniverse pipelines, or RAPIDS-adjacent projects.

### Gaps
- No a16z or other named VC thesis on "spatial AI" was retrieved in this pass.
- No 2024–2026 funding figures were retrieved for Wherobots, CARTO, Clay or other geospatial-infra startups.
- Market reports are paywalled vendor estimates with divergent methods. None isolates "AI-derived vector data production" as a segment.

---

## 6. 4D / real-time: live twins, continuous monitoring, crowdsourced maps, change detection as a service

### Takeaway
Sensing cadence is collapsing to daily or sub-daily: Planet Pelican at 30 cm with on-board Jetson, Mobileye REM covering 15M km of road per day, and monthly embeddings and building releases. Inference cost has dropped to dollars per thousand km² (Wherobots + SAM3). Continuous model output means continuous mask → vector → diff work. Notably, Mobileye already shows that the winning pattern is shipping **geometry, not pixels**.

### Cited Findings
- **Mobileye REM**: 8M+ vehicles logged **54.7 billion km** in the prior year. **15M+ km of roads mapped daily**. Only about **10 KB per km** is uploaded. Coverage exceeds 95% of public roads in the US and Europe, and 5 of the 10 largest automakers contribute. Stellantis integrates REM from **2027**. — [Mobileye REM](https://mobileye.com/technology/rem); [TimesLive, Jul 21 2026](https://www.timeslive.co.za/motoring/2026-07-21-mobileye-to-supply-stellantis-with-its-cloud-based-driver-assistance-tech/); [Pulse 2.0](https://pulse2.com/mobileye-to-supply-cloud-enhanced-adas-technology-for-select-stellantis-vehicles-starting-in-2027/)
- **Planet Pelican**: next-generation Pelicans target **30 cm class** imagery with 6 multispectral bands and carry **NVIDIA Jetson** for on-orbit edge compute. Pelican-12 is the **fifth Pelican launched in 2026**, and production is scaling in Berlin. Planet's fleet images Earth's whole landmass daily (Doves). — [Planet Pelican](https://www.planet.com/constellations/pelican); [Webull/Planet Pelican-12](https://www.webull.com/news/15639162848166912); [BusinessWire Pelican-6](https://www.businesswire.com/news/home/20251209493743/en/Planet-Releases-First-Light-Image-From-Pelican-6)
- **Wherobots RasterFlow + NVIDIA (Sep 30, 2026, public preview)**: object detection with **Meta SAM3 on 30 cm imagery for under $25 per 1,000 km²**. Sentinel-2 cloud-free mosaics plus inference for **under $1 per 10,000 km²**. 2–3× cheaper than alternatives. Production uses: **Fields of the World** (field boundaries; Taylor Geospatial, Microsoft AI for Good, NASA Harvest) and **FireCon** (US Forest Service wildfire containment, H2 2026). The post says **nothing about polygonization or post-processing cost**. — [Wherobots + NVIDIA blog](https://wherobots.com/blog/wherobots-nvidia-geospatial-ai-physical-world/)
- In Wherobots' Tile2Net flow, the model outputs a 4-class raster (background, road, sidewalk, crosswalk), and "a separate flow" converts classes to vector geometries based on a confidence threshold. — [Wherobots RasterFlow Tile2Net](https://wherobots.com/solution-gallery/rasterflow-tile2net/); [Wherobots raster inference GA](https://wherobots.com/blog/wherobotsai-raster-inference-is-ga/)
- Overture buildings release monthly ([Overture](https://docs.overturemaps.org/blog/2026/09/23/release-notes/)). Google Population Dynamics embeddings update monthly ([Google Research](https://research.google/blog/google-earth-ai-unlocking-geospatial-insights-with-foundation-models-and-cross-modal-reasoning/)). AEF is annual ([Source Cooperative](https://source.coop/tge-labs/aef)).
- Singapore port twin: 1,000+ sensors (H2 2025). — [OECD OPSI](https://oecd-opsi.org/innovations/virtual-twin-singapore/) (search snippet)

### Inferences
- (Inference) As inference costs fall to about $25/1,000 km² at 30 cm, the **relative** share of cost and latency in post-processing (vectorize → simplify/regularize → topology-check → conflate → diff) rises. Wherobots optimized upload/compute/download with CUDA streams for inference but does not mention vectorization, which suggests it is still CPU or Spark-side. This is unconfirmed.
- (Speculation, high strategic value) **On-orbit or edge vectorization**: Pelicans carry Jetsons, and Mobileye proves the economics of sending 10 KB/km of geometry instead of imagery. A lightweight CUDA raster→polygon engine running on Jetson-class hardware (satellite, drone, vehicle) could downlink **polygons of change** instead of image tiles, which bears directly on bandwidth-limited sensors.
- (Inference) "Change detection as a service" becomes **polygon diffing at cadence**. Each new pass yields masks, which yield rings, which are diffed against the previous version. Deterministic, topologically consistent output from a scanline engine (same input gives the same rings, and halos guarantee seamless tiles) matters for stable diffs and versioning.

### Gaps
- No source quantifies what fraction of RasterFlow or similar pipeline time or cost goes to vectorization.
- No 2026 data retrieved on drones-as-a-service volumes (Spexi is mentioned only as a Niantic partner).
- No source found confirming any operator doing on-board vectorization today.

---

## 7. Predicted bottlenecks: vectorization compute, topology QA, conflation, versioning, provenance and uncertainty

### Takeaway
There is documented evidence of (a) a quality gap between masks and polygons, (b) vectorization being memory- and compute-heavy in the dominant cloud platform (GEE), (c) a vacated GPU geometry library slot (cuSpatial discontinued 2025), and (d) a lack of uncertainty and explainability tooling for embedding-derived outputs. Conflation of multi-source ML geometry is a live problem at Overture scale. Versioning and provenance of AI polygons are mostly undocumented. That is a gap, and an opportunity.

### Cited Findings
- **Mask → polygon degradation**: high-IoU masks lose average precision after vectorization, **e.g. 79% → 67% AP on the AICrowd benchmark**. Classic Douglas–Peucker, marching squares and hand-tuned regularization "could not bridge the gap". — [Spyrosoft on HiSup](https://spyro-soft.com/blog/geospatial/closing-the-geometry-gap-in-building-footprint-extraction-our-first-look-at-hisup) (secondary source). Imperfect masks from shadows, occlusion and raster quantization yield "suboptimal polygons during the vectorization process" — [arXiv 2407.12256](https://arxiv.org/html/2407.12256v1). Research is moving toward **direct polygon prediction** (PolyWorld, BiSVP) to avoid post-hoc vectorization — [PolyWorld arXiv 2111.15491](https://arxiv.org/pdf/2111.15491); [BiSVP arXiv 2303.00300](https://arxiv.org/pdf/2303.00300)
- **Post-segmentation QA as its own research topic (2026)**: "GeoAI-based post-segmentation quality validation of building footprints via spatial feature engineering". — [arXiv 2608.09048](https://arxiv.org/pdf/2608.09048) (title only; the PDF was too large to fetch, so details are unverified)
- **GEE vectorization cost**: "Vectorization is a memory and compute-intensive operation", with timeouts and tile errors at scale. — [Earth Engine tutorial](https://developers.google.com/earth-engine/tutorials/community/satellite-embedding-05-similarity-search)
- **GPU geometry library vacuum**: RAPIDS **stopped publishing cuSpatial packages as of v25.06 (June 5, 2025)**. Development is paused, and it was removed from the RAPIDS metapackage and containers. — [RAPIDS Support Notice RSN 45](https://docs.rapids.ai/notices/rsn0045/)
- Meanwhile GPU spatial **joins** are advancing: SedonaDB runs joins on **RT cores**, with a 5.93× speedup and 59% cost reduction. — [Wherobots + NVIDIA](https://wherobots.com/blog/wherobots-nvidia-geospatial-ai-physical-world/)
- **Uncertainty and explainability**: embedding products lack "explainability and uncertainty quantification" and standardized benchmarks. — [arXiv 2601.13134](https://arxiv.org/html/2601.13134v1). GEE similarity results carry false positives that need manual threshold tuning. — [Earth Engine tutorial](https://developers.google.com/earth-engine/tutorials/community/satellite-embedding-05-similarity-search)
- **Conflation at scale**: Overture merges overlapping building sources (Google ML, Microsoft ML, OSM, Esri, national agencies) monthly and assigns stable **GERS** IDs. Microsoft reported accuracy and coverage gains from adopting the conflated layer. — [Overture release notes](https://docs.overturemaps.org/blog/2026/09/23/release-notes/); [Overture announcement](https://overturemaps.org/announcements/2026/overture-maps-foundation-reaches-50-members-as-industry-converges-on-open-data-to-ground-ai/)
- **Human-in-the-loop fragility**: Rapid, the main OSM AI-feature intake tool, went mostly dormant for about 9 months in 2025 after Meta's layoff. — [OSM Wiki: Rapid](https://wiki.OpenStreetMap.org/wiki/RapiD)

### Inferences
- (Inference) **Strategic position for the engine.** The emerging stack is: GPU inference (cheap, solved) → **[gap: GPU raster→vector, topology-correct, tile-seamless]** → GPU spatial joins (SedonaDB on RT cores, solved) → Parquet/Iceberg native geometry (standardized). The engine fills the missing middle link. The data never has to leave the GPU between model output and join.
- (Inference) The documented mask→polygon AP drop means a raw pixel-edge ring tracer alone does not solve **quality**. Its value is as a fast, exact, deterministic **first stage** that feeds regularization (HiSup or PolyWorld-style), simplification and QA. Exact pixel-boundary rings with correct holes and no tile seams are the ground truth that later stages need.
- (Inference) **Versioning and change**: deterministic scanline output with halo-consistent tiling allows content-addressable geometry (hash rings per tile), which enables cheap diffs between monthly or daily releases. This lines up with Overture's monthly cadence and GERS stable IDs.
- (Speculation) **Provenance and uncertainty polygons**: run the engine at multiple thresholds on the same similarity or probability raster to emit **nested confidence contours** (for example p>0.5, 0.7, 0.9 rings). That yields "uncertainty bands" for AI polygons cheaply, which addresses the uncertainty gap flagged by the Earth Embeddings survey. The pluggable transform could also attach model/threshold metadata per ring.

### Gaps
- No source quantifies the global compute spent on raster→vector conversion, or its share of GeoAI pipeline cost.
- No source found describing best practice for versioning AI-derived polygons, or for publishing per-polygon uncertainty, in Overture, OSM or national mapping agencies.
- No primary evidence found that topology QA (self-intersections, slivers, gaps between tiles) is a bottleneck at Overture scale. This is plausible but unverified.
- The details of arXiv 2608.09048 (post-segmentation QA) could not be retrieved.
