# Mask → Polygon Conversion in the AI / Data Infrastructure Stack (2023–2026)

Scope: where binary masks / occupancy rasters become closed polygon rings in modern AI pipelines (labeling, foundation-model segmentation serving, VLM geometry outputs, storage/analytics, video analytics, synthetic data, vectorization-as-a-service). This is for a CUDA/C++ single-pass, 2x2-classification, halo-tiled mask→ring engine with a pluggable coordinate transform. Research date: 2026-10-07. Source quality is noted where it matters.

---

## Q1. Labeling platforms: SAM-assisted polygon tools, conversion quality, COCO polygon vs RLE, auto-labeling at scale

### Takeaway
Every major labeling platform now produces SAM 2/SAM 3 masks and then converts them to polygons, usually on request. That conversion is still done the same ad-hoc way everywhere: `cv2.findContours` followed by `approxPolyDP`/Douglas-Peucker on the CPU. It causes known, recurring problems: holes get dropped or encoded badly, vertex counts blow up, rounding drift appears on mask↔polygon round trips, and COCO polygon and RLE don't interoperate. No platform documents a hole-correct, topology-preserving, GPU-batched converter.

### Cited Findings
- CVAT's SAM 2 tracker (Enterprise) has a "Convert the mask into a polygon" option the user must enable to track and edit SAM output as polygons. — [CVAT SAM 2 tracker docs](https://docs.cvat.ai/docs/enterprise/segment-anything-2-tracker); [CVAT blog: video annotation with SAM 2](https://www.cvat.ai/blog/video-annotation-sam-2)
- A community write-up of CVAT + SAM 2 auto-annotation says raw SAM 2 masks can produce ~2000-point polygons that are hard to edit. Its fix is `cv2.findContours(..., RETR_EXTERNAL, CHAIN_APPROX_SIMPLE)` and then `cv2.approxPolyDP(epsilon=2.0)`. RETR_EXTERNAL throws away holes. (Low-authority personal site; treat the "2000 points" figure as anecdotal.) — [cv-sam2-autoannotate](https://neytro-org.duckdns.org/a/cv-sam2-autoannotate.html)
- Encord offers polygons, polylines, bitmasks and keypoints with SAM 2 built in. SAM 3 was in early access in early 2026, adding text-prompt and class-name click annotation that labels all instances of a class at once. Encord's SDK documents converting polygons to bitmasks. — [Encord annotate images](https://docs.encord.com/platform-documentation/Annotate/annotate-label-editor/annotate-images); [Encord SDK: converting ontology shapes](https://docs.encord.com/sdk-documentation/projects-sdk/sdk-converting-ontology-shapes)
- V7 lets SAM create either Polygon or Mask annotations and added one-click polygon→mask conversion (changelog, Jul 2024). — [V7 SAM docs](https://docs.v7labs.com/docs/segment-anything-model-sam-for-images); [V7 changelog: convert polygon to mask](https://v7darwin.com/changelog/convert-polygon-to-mask)
- Roboflow supervision: a user reported that masks "become crooked and noisy" after YOLO→mask conversion with `force_masks=True`. The suspected cause was integer rounding in the mask→polygon→mask round trip. Roboflow staff acknowledged it and passed it to the supervision repo as a precision issue, calling it "minimally impactful on model learning, though requires investigation." — [Roboflow forum thread](https://discuss.roboflow.com/t/masks-become-crooked-and-noisy-after-passing-through-supervision-library/3269)
- Roboflow forum threads report "Smart Polygon tool not working … 'Incorrect Conversion' pop-up" and "hosted inference for YOLO instance segmentation only returns one polygon per prediction". The second is the multi-part / disconnected-mask problem. — [Smart Polygon issue](https://discuss.roboflow.com/t/issue-smart-polygon-tool-not-working-bounding-box-disappears-incorrect-conversion-pop-up/9113); [one polygon per prediction](https://discuss.roboflow.com/t/hosted-inference-for-yolo-instance-segmentation-only-returns-one-polygon-per-prediction/9648)
- supervision issue #1114 (opened Apr 12, 2024, closed as done in Q2 2024) added COCO RLE support. Before that the library could only read and write COCO polygon masks. — [supervision #1114](https://github.com/roboflow/supervision/issues/1114)
- In COCO, `segmentation` can be polygons (a list of lists) or RLE (a dict, usually when `iscrowd=1`). Converters such as Label Studio's COCO importer have to handle both. — [Label Studio community: COCO to LS](https://community.labelstud.io/t/converter-coco-to-ls/647); [Detectron #100: RLE vs polygon](https://github.com/facebookresearch/Detectron/issues/100); [cocoapi #39](https://github.com/cocodataset/cocoapi/issues/39)
- To convert an RLE mask with holes to YOLO-seg format, a third-party script makes parent (clockwise) and child (counter-clockwise) polygons with `cv2.findContours` and stitches them together. So holes have no native representation in the YOLO polygon format. — [ryouchinsa/donut (general_json2yolo.py)](https://github.com/ryouchinsa/donut)
- SAHI PR #599 converts RLE→polygon by running `mask_utils.decode()` and then `cv2.findContours()`. This is the standard idiom. — [SAHI PR #599](https://github.com/obss/sahi/pull/599)
- Roboflow publishes a dedicated guide to converting boxes, masks and polygons, which shows this is a frequent user need. — [Roboflow blog](https://blog.roboflow.com/convert-bboxes-masks-polygons/)

### Inferences
- The industry-wide "conversion kernel" is OpenCV `findContours` (Suzuki-Abe border following, CPU, one mask at a time) plus Douglas-Peucker. Its known weaknesses are exactly what the user's engine targets: hole/parent-child topology, multiple disjoint parts, deterministic closed rings, and throughput on batches of thousands of masks.
- SAM 3 "label every instance of a concept" (Encord, others) multiplies masks per image by 10–100x compared with click-per-object. If polygons are the export format, mask→polygon becomes a bulk batch job rather than a UI-latency job. A GPU batched converter that keeps masks on-device (no D2H copy of full-resolution masks) is a credible auto-labeling infrastructure piece.
- The COCO polygon vs RLE split is really an argument about lossless-with-holes (RLE) vs editable-and-compact (polygon). A converter that guarantees exact ring topology (outer plus holes, pixel-edge exact before simplification) narrows the gap.

### Gaps
- No public numbers found for Scale AI, Labelbox or SuperAnnotate on mask→polygon conversion (vertex budgets, simplification tolerances, or throughput of auto-labeling exports).
- I couldn't fetch a canonical CVAT GitHub issue that quantifies polygon-point explosion. The 2000-point claim comes from a low-authority source.
- No platform publishes how long mask→polygon export takes per million annotations.

---

## Q2. SAM 2 / SAM 3 / SAM 3.1 and video segmentation: throughput, downstream use of outlines, edge latency

### Takeaway
Foundation segmentation models now produce masks at tens of frames per second for dozens to 100+ objects (SAM 2 at about 44 FPS, SAM 3 at 30 ms per image with 100+ objects on H200, SAM 3.1 tracking 16 objects at 32 fps). Their native output is per-object full-resolution masks, never polygons. Any consumer that wants outlines (editing, AR overlays, zone logic, geometry storage, LLM context) pays for a separate conversion step. Because the conversion is serial CPU work, it scales with object count × frame rate and can fall behind the GPU model.

### Cited Findings
- SAM 2 (Aug 2024) runs at 43.8 FPS (Hiera-B+) and 30.2 FPS (Hiera-L) on one A100, using a streaming-memory transformer for video. Checkpoints range from 38.9M to 224.4M parameters. — [SAM 2 paper, arXiv 2408.00714](https://arxiv.org/pdf/2408.00714); [facebookresearch/sam2](https://github.com/facebookresearch/sam2)
- SAM 3 (released Nov 19, 2025) is Promptable Concept Segmentation: it detects, segments and tracks all instances of a concept given a text phrase and/or image exemplars. It has 473.6M parameters (3.45 GB), a "presence token" and a decoupled detector-tracker design. — [Ultralytics SAM 3 docs](https://docs.ultralytics.com/models/sam-3)
- SAM 3 speed: about 30 ms per image with 100+ detected objects on H200. Ultralytics separately measured 2,921 ms per image on an RTX PRO 6000 through its own integration, a big gap from Meta's figure that probably reflects setup and settings. Video tracking is near real-time for about 5 concurrent objects. — [Ultralytics SAM 3 docs](https://docs.ultralytics.com/models/sam-3)
- SAM 3.1 (Mar 27, 2026) adds "Object Multiplex" for multi-object video tracking, with a reported ~7x speedup at 128 objects on H100. It tracks 16 objects in one pass at 32 fps. — [Ultralytics SAM 3 docs](https://docs.ultralytics.com/models/sam-3); [Voxel51 glossary / summary](https://voxel51.com/glossary/sam-3-segment-anything-model-3)
- SAM 3 benchmarks: SA-Co/Gold cgF1 65.0; LVIS zero-shot mask AP 47.0; MOSEv2 J&F 60.1; DAVIS17 J&F 92.0. The data engine auto-annotated over 4M unique concepts. — [Ultralytics SAM 3 docs](https://docs.ultralytics.com/models/sam-3)
- SAM 3 outputs: instance masks, boxes and track IDs. No polygon output. — [Ultralytics SAM 3 docs](https://docs.ultralytics.com/models/sam-3)
- Edge: Ultralytics recommends YOLO26n-seg instead of SAM 3 on constrained devices ("231x faster"). EdgeSAM distils SAM to a CNN encoder and runs at 30+ FPS on phones and Jetson-class devices, more than 7x faster than MobileSAM/EfficientSAM on smartphones. — [Ultralytics SAM 3 docs](https://docs.ultralytics.com/models/sam-3); [EdgeSAM arXiv 2312.06660](https://arxiv.org/html/2312.06660v3)
- Voxel51 FiftyOne ships SAM 3 plugins for images, so SAM 3 masks go straight into dataset curation tools. — [FiftyOne sam3_images plugin](https://docs.voxel51.com/plugins/plugins_ecosystem/sam3_images.html)

### Inferences
- Rough arithmetic: SAM 3.1 at 128 objects × 30 fps is about 3,800 masks per second per stream. At that rate a per-mask CPU `findContours` (roughly 0.1–1 ms per mask at 1080p, typical figure, not measured here) uses 0.4–4 CPU-seconds per second per stream. A GPU converter that runs on the label/instance map already in device memory, ideally one pass over a multi-label map rather than N binary masks, removes the bottleneck.
- On Jetson Orin/Thor, the CPU cores are weak compared with the GPU. Running mask→polygon on the GPU right after EdgeSAM or YOLO-seg avoids copying large masks to the host and keeps the payload small (polygons are bytes; masks are megabytes). That matters for streaming metadata to the cloud (Kafka/MQTT).
- Downstream consumers that want outlines include: annotation editing (CVAT/Encord), AR occlusion and outline rendering, zone and intersection logic (Q5), storing tracks as geometry in spatial databases (Q4), and passing compact geometry to LLM agents.

### Gaps
- No published benchmark of the mask→polygon step's share of SAM 2/3 pipeline latency was found. This is an opportunity for a measured claim.
- No authoritative Jetson Thor numbers for SAM 2/3 or EdgeSAM were found within budget.
- The Ultralytics 2,921 ms figure and Meta's 30 ms figure are hard to reconcile from the sources. Treat the Ultralytics number as integration-specific.

---

## Q3. Vision-language models emitting geometry: Gemini, Florence-2, PaliGemma, Moondream, polygon-as-tokens

### Takeaway
VLMs are splitting into three ways of outputting geometry. (a) Quantized polygon vertex tokens: Florence-2, PolyFormer-style models, newer Gemini docs. (b) Learned mask codebook tokens: PaliGemma's 16 VQ-VAE tokens per mask. (c) Vector paths: Moondream's SVG output. Polygons win on token efficiency, interpretability and composability with text. Training such models requires converting masks to polygons (findContours + Douglas-Peucker + quantize), so mask→polygon is now a training-data preprocessing step for VLMs. Production behaviour is still inconsistent: Gemini's format has changed across versions.

### Cited Findings
- Gemini 2.5 docs describe a JSON response of `box_2d` ([y0,x0,y1,x1], normalized 0–1000), `label`, and `mask` as a base64 PNG probability map (0–255) clipped to the box. — [Simon Willison: Gemini 2.5 image segmentation (Apr 2025)](https://simonwillison.net/2025/Apr/18/gemini-image-segmentation/); [Gemini image understanding docs](https://ai.google.dev/gemini-api/docs/image-understanding)
- Google developer forum (2026) reports that Gemini 2.5 Flash returns raw codebook strings like `"<start_of_mask><seg_4><seg_20>"` instead of the documented base64 PNG. Gemini 3.5 Flash and Gemini Robotics-ER 1.6 return masks as polygons, `[[325,411],[327,471],...]`. Release notes quoted in the thread say "Image segmentation is not supported in Gemini 3.x. For segmentation workloads, continue using Gemini 2.5 Flash." Users call the masks "unusable" and ask whether segmentation has been abandoned. No Google staff reply. — [Google AI dev forum thread](https://discuss.ai.google.dev/t/image-understanding-and-segmentation-mask-support/146559)
- **Conflict:** the current Gemini image-understanding docs, as fetched 2026-10-07, describe `mask` as "a polygon of [x,y] coordinates, normalized to 0-1000" and recommend thinking level "minimal" for better results. The fetch summary named a recent 3.x Flash model in the examples; I couldn't verify the exact model string. This contradicts the "not supported in 3.x" release-note quote above. The direction is clear: Gemini has moved from PNG masks to polygons. — [Gemini image understanding docs](https://ai.google.dev/gemini-api/docs/image-understanding)
- Florence-2's `<REFERRING_EXPRESSION_SEGMENTATION>` task outputs polygons as quantized coordinate tokens `<X1><Y1>…<Xn><Yn>`. — [Roboflow Florence-2 docs](https://docs.roboflow.com/models/ja/supported-models/florence-2); [fal.ai Florence-2 RES](https://fal.ai/models/fal-ai/florence-2-large/referring-expression-segmentation/llms.txt); [TDS Florence-2 explainer](https://towardsdatascience.com/florence-2-mastering-multiple-vision-tasks-with-a-single-vlm-model-435d251976d0/)
- PaliGemma adds 128 `<seg000>…<seg127>` tokens. Each mask is 16 VQ-VAE codebook tokens decoding to a 64×64 binary mask inside the box given by 4 location tokens. — [Google Developers blog: PaliGemma architecture](https://developers.googleblog.com/gemma-explained-paligemma-architecture/)
- Moondream 3 Preview (segment skill launched Sep 2025) returns "native SVG masks (vectors, not bitmasks)" plus a box. The Mar 2026 update reports RefCOCO val 83.2 mIoU, RefCOCO+ 79.1, RefCOCOg 80.7, RefCOCO-M 88.2, and runs 40% faster. Moondream positions it as crisper than SAM 3. — [Moondream blog, 2026-03-10](https://moondream.ai/blog/segmenting-update-2026-03-10); [Moondream segment docs](https://docs.moondream.ai/skills/segment)
- A 2026 arXiv paper on VLM geometry grounding builds polygon tokens by extracting contours with OpenCV `findContours` plus Douglas-Peucker, quantizing vertices to [0,255], and wrapping them in special delimiters. It argues polygons are "substantially more token-efficient", interpretable and composable, with "a controllable token budget that scales linearly with the number of vertices rather than the image resolution." It also says Moondream segmentation "autoregressively decodes a vector path and iteratively refines the rasterized mask." — [arXiv 2604.02893](https://arxiv.org/pdf/2604.02893)
- Earlier sequential polygon predictors (Polygon-RNN++, PolyFormer) are cited as the lineage. — [arXiv 2604.02893](https://arxiv.org/pdf/2604.02893)

### Inferences
- Why VLMs prefer polygons: token cost scales with the number of vertices, not with resolution; the output is valid text/JSON; it composes with reasoning; and it can be checked by rasterizing. The cost is losing holes and multi-part shapes unless the format supports them (Gemini's single `[[x,y],...]` list appears to be one outer ring), plus quantization error (0–1000 or 0–255 grids).
- This opens two infrastructure roles for a fast, correct converter. (1) Training data factories: converting SA-1B / SA-Co-scale masks (billions) to vertex-budgeted, hole-aware polygon token sequences for VLM pretraining. That is a batch GPU ETL job where quality (topology, ring orientation, deterministic start vertex) directly affects training targets. (2) Serving: turning SAM-family masks into compact polygon context for LLM agents, which is cheaper than sending PNG masks.
- A pluggable coordinate transform maps directly onto VLM quantization grids (0–1000 normalized, 0–255 bins, box-relative coordinates), so one kernel could emit model-ready tokens.

### Gaps
- I didn't fetch VistaLLM or PolyFormer papers directly. Claims about them come via secondary citation.
- The Gemini model and format matrix is contradictory across sources (forum vs current docs). I couldn't verify the exact model name in the current docs.
- No study found that directly compares mask tokens vs polygon tokens for VLM segmentation quality at a fixed token budget.

---

## Q4. Data formats, storage and GPU spatial analytics (GeoParquet/GeoArrow, COCO, Overture, Wherobots, Sedona, cuSpatial)

### Takeaway
Polygons are what the analytics layer runs on. GeoParquet (for example Overture's 2B+ buildings) and spatial SQL engines consume them, and the vectors themselves increasingly come from ML segmentation of imagery followed by polygonization. GPU spatial joins have moved away from cuSpatial (publishing stopped June 2025) toward ray-tracing-core joins in SedonaDB 0.4 (2026). The raster→vector step (`RS_Polygonize`, GDAL-style) is still CPU and per-tile in these stacks: a gap between GPU inference and GPU joins.

### Cited Findings
- RAPIDS stopped publishing cuSpatial packages from release v25.06 (June 5, 2025). Development is paused, and it was dropped from the RAPIDS metapackage and containers (RSN 45, completed by Jul 17, 2025). — [RAPIDS RSN 45](https://docs.rapids.ai/notices/rsn0045/)
- SedonaDB (introduced Sep 2025) is a single-node, spatial-first analytical database. Version 0.4.0 (Jun 2026) adds GPU spatial joins through "RayBooster", which runs on NVIDIA RT cores. Its design: Structure-of-Arrays geometry instead of WKB, a single global BVH with geometry IDs packed into the Z axis, a DE-9IM RelateEngine on RT cores, and memory-aware scheduling. Results: up to 5.93x speedup (59% AWS cost reduction); 4.93–9.68x on heavy joins; an RTX 3090 beat an H100 on some queries (1.26 s vs 1.77 s, Q10). Requires CUDA 12+ and compute capability ≥7.5. VLDB 2026 Industry Track. The blog does not mention polygonization or raster→vector. — [SedonaDB 0.4 GPU joins blog](https://sedona.apache.org/1.9.1/blog/2026/06/26/sedonadb-04-gpu-accelerated-spatial-joins/); [SedonaDB GPU docs](https://sedona.apache.org/sedonadb/latest/gpu-acceleration/); [Introducing SedonaDB](https://sedona.apache.org/1.9.0/blog/2025/09/24/introducing-sedonadb-a-single-node-analytical-database-engine-with-geospatial-as-a-first-class-citizen/)
- Wherobots `RS_Polygonize` turns each connected region of equal-valued pixels into a polygon with its value. RasterFlow is described as a production inference engine that "detects buildings, roads, vegetation, water, and land cover from satellite imagery and returns analysis-ready vector features" ready for spatial joins, with GPU serverless inference. — [Wherobots RS_Polygonize](https://docs.wherobots.com/reference/wherobots-db/raster-data/pixel-functions/RS_Polygonize.md); [RasterFlow models](https://docs.wherobots.com/develop/rasterflow/rasterflow-models.md); [Wherobots raster inference overview](https://docs.wherobots.com/tutorials/wherobotsai/wherobots-inference/raster-inference-overview); [HPCwire: Wherobots raster inference](https://www.hpcwire.com/bigdatawire/this-just-in/wherobots-introduces-raster-inference-to-accelerate-geospatial-ai-across-industries/)
- Overture buildings are distributed as GeoParquet (e.g. 2025-01-22 release on AWS/Azure). Sources include OSM (~635M), Microsoft ML buildings (~738M), Google Open Buildings (~350M high-precision + ~650M lower-precision), Esri (~15M), IGN Spain (~13M) and East-Asia ML buildings (~228.5M added in 2025). Most of these are ML mask→polygon outputs. The theme reached 2.3B buildings in 2023 after adding Google Open Buildings. — [Overture 2025-01-22 release notes](https://docs.overturemaps.org/blog/2025/01/22/release-notes/); [Overture buildings guide](https://docs.overturemaps.org/guides/buildings); [Overture 2023: 2.3B buildings](https://overturemaps.org/blog/2023/overture-buildings-theme-hits-2-3b-buildings-with-addition-of-google-open-buildings-data/)
- Ecopia AI: AI-extracted vector features over 40M+ km² in 100+ countries; 176M+ US building footprints. — [Ecopia blog](https://ecopiatech.com/resources/blog/digitizing-the-world-with-ecopia-ai-based-mapping-system-global-partner-network); [geo-matching Building Footprints](https://geo-matching.com/products/building-footprints)
- Raster Vision ships `mask_to_building_polygons` as a vectorization utility, another example of every geo-ML framework writing its own converter. — [Raster Vision vectorization](https://docs.rastervision.io/en/0.31/_modules/rastervision/core/data/utils/vectorization.html)

### Inferences
- The pipeline is: GPU inference (masks) → CPU polygonize (GDAL/OpenCV/RS_Polygonize) → GeoParquet → GPU or RT-core spatial joins. The middle step is the one non-GPU hop. A GPU polygonizer that writes GeoArrow-style SoA coordinate buffers (which SedonaDB's RayBooster also uses instead of WKB) could make a zero-copy, all-GPU path from inference to join.
- Halo-buffered tiling matters most here. Planetary inference is tiled, and the standard failure is polygons cut at tile seams, needing a post-hoc dissolve. A converter with halos that produces seam-consistent rings is a differentiator. The pluggable transform covers pixel→CRS (affine geotransform) directly.
- With cuSpatial paused, NVIDIA's own GPU geometry stack is thin. This is an opening, not a threat.

### Gaps
- No current figures found for Lance, Zarr or vector-DB indexing of shapes, or for DuckDB spatial raster polygonization performance.
- Not covered within budget: Heavy.AI/HeavyDB, Kinetica, CARTO and Foursquare GPU geometry capabilities.
- The current (2026) Overture building total wasn't confirmed. The 2.3B (2023) figure is the latest firm number found.
- No benchmark found for RS_Polygonize throughput.

---

## Q5. Video analytics / smart cities: zone polygons, occupancy, heatmaps (NVIDIA Metropolis/DeepStream)

### Takeaway
In production video analytics, polygons are hand-drawn configuration inputs (ROIs, lines). Objects are reduced to a single point (the bottom-centre of the bbox) and tested against those zones. Turning occupancy or heatmap rasters and segmentation masks into zone polygons automatically is largely absent from mainstream SDKs.

### Cited Findings
- DeepStream `gst-nvdsanalytics` supports ROI filtering, overcrowding, direction detection and line crossing. ROIs are configured as semicolon-separated polygon coordinates ("simple polygon with more than 3 points"). All analytics use the bottom-centre of the object's bbox. Direction and line crossing need tracker IDs. Metadata includes per-ROI object counts, cumulative line-crossing counts and overcrowding flags. — [DeepStream nvdsanalytics docs](https://docs.nvidia.com/metropolis/deepstream/dev-guide/text/DS_plugin_gst-nvdsanalytics.html); [NvDsAnalyticsObjInfo (DS 9.1)](https://docs.nvidia.com/metropolis/deepstream/9.1/sdk-api/structNvDsAnalyticsObjInfo.html)

### Inferences
- Obvious gaps: (1) mask-aware zone occupancy, meaning the area fraction of an instance mask inside a zone rather than the bottom-centre point (useful for parking spots, shelf fill, crowd density); (2) automatically deriving zones from long-run occupancy or heatmap rasters (threshold, then GPU polygonize), producing "learned zones" like walkways, queue areas or dwell hotspots; (3) emitting per-frame instance outlines as compact metadata instead of masks. All three need fast on-device mask→polygon on Jetson.
- Retail heatmap vendors usually display rasters. Converting them to polygons enables SQL/GIS joins with floor plans, which the user's coordinate transform (pixel → floor-plan homography) supports directly.

### Gaps
- Not researched within budget: Metropolis microservices / VSS (video search and summarization) handling of segmentation outlines, and vendor products (Hailo, Axis, Milestone) for auto-zone generation.
- No source found on heatmap→polygon zoning products.

---

## Q6. Synthetic data and simulation (Omniverse Replicator, Isaac Sim)

### Takeaway
Synthetic data generators render perfect per-instance segmentation buffers and export them through writers (Basic, COCO, KITTI). COCO export means rasterized instance IDs must become polygons or RLE at dataset scale, which is another batch mask→polygon workload that today runs in Python writers.

### Cited Findings
- Omniverse Replicator is the core SDG extension. Writers (BasicWriter, CocoWriter, KittiWriter) take annotators such as `rgb`, `semantic_segmentation`, `instance_segmentation` and `bounding_box_3d`. CocoWriter takes `coco_categories`. — [Omniverse Replicator docs](https://docs.omniverse.nvidia.com/kit/docs/omni_replicator/latest/index.html); [Isaac Sim scene-based SDG tutorial](https://docs.isaacsim.omniverse.nvidia.com/latest/replicator_tutorials/tutorial_replicator_scene_based_sdg.html)

### Inferences
- SDG pipelines often render millions of frames with dozens to hundreds of instances each. Converting instance-ID buffers to COCO polygons, ideally in one pass over the multi-label buffer straight from the renderer's GPU memory, is a natural fit for a scanline GPU engine.

### Gaps
- I didn't confirm whether Replicator's CocoWriter emits polygons or RLE, or how fast it runs. Not checked: Isaac Sim 5.x/6.x release notes, Parallel Domain, Bifrost, Duality and other SDG vendors.

---

## Q7. Vectorization / geometry-extraction startups and services

### Takeaway
"Vectorization as a service" exists mainly in geospatial (Ecopia, Blackshark.ai, Wherobots RasterFlow) and in graphic design (Vectorizer.AI). In both cases polygonization is an internal step bundled into a full product. No company offers a general GPU mask→polygon primitive or SDK for the ML stack.

### Cited Findings
- Blackshark.ai Orca runs as a service on Azure Synapse. Its building-footprint container produces probability heatmaps "which can then be vectorized". — [Microsoft/Blackshark Orca doc](https://query.prod.cms.rt.microsoft.com/cms/api/am/binary/RWZVZB)
- Ecopia continuously ingests imagery, extracts vector features and detects change. — [Ecopia blog](https://ecopiatech.com/resources/blog/digitizing-the-world-with-ecopia-ai-based-mapping-system-global-partner-network)
- Wherobots RasterFlow returns vector features from satellite inference. — [RasterFlow models](https://docs.wherobots.com/develop/rasterflow/rasterflow-models.md)
- Planet/Sentinel Hub users have asked for raster→vector in the Processing API on community forums, which shows unmet demand. — [Planet community: raster to vector in Processing API](https://community.planet.com/analysis-apis-81/raster-to-vector-in-processing-api-3535); [Sentinel Hub forum](https://forum.sentinel-hub.com/t/raster-to-vector-in-processing-api/7664)
- Vectorizer.AI offers a bitmap-tracing API for graphics, not geospatial or ML masks. — [Vectorizer.AI about](https://tw.vectorizer.ai/about)
- Prior GPU vectorization research exists, e.g. "Real-time Image Vectorization on GPU" (boundary pixel detection → connectivity → loop organization → segment fitting), but no maintained production library surfaced. — [Real time Image Vectorization on GPU](https://www.academia.edu/144396955/Real_time_Image_Vectorization_on_GPU)

### Inferences
- The gap is a reusable, correct, GPU-native polygonization primitive: hole-aware, multi-part, seam-consistent across tiles, emitting GeoArrow/COCO/VLM-token outputs. It would plug in between SAM/VLM/segmentation inference and every downstream consumer (labeling exports, spatial databases, VLM training data, edge analytics).

### Gaps
- Picterra, Mapflow, EOSDA and similar firms' vectorization internals weren't researched. No funding or market sizing found for vectorization-specific startups.
