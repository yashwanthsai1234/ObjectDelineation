# Raster -> Vector Polygon Conversion in Fabrication, Manufacturing, Energy, Infrastructure and Other Physical-World Industries (2023-2026)

Scope note: this covers the "long tail" outside semiconductor lithography, geospatial EO, biomedical and robotics. The engine under study takes a binary occupancy bitmap with halo buffers for tiling, does one scanline pass with 2x2 pixel-type classification, and emits closed polygon rings through a pluggable coordinate transform. Budget limits: about 18 tool calls. Several sub-domains (power-line vegetation, wind blade inspection, mining/ore sorting, coral reefs, P&ID, SDF fonts, fluid level sets) got no dedicated searches and are listed under Gaps rather than filled from memory.

---

## Q1. Additive manufacturing: in-situ layer monitoring -> contour polygons vs CAD (closed-loop AM), voxel/bitmap printing, voxel<->contour formats

### Takeaway
Comparing layer contours against CAD is an active 2025-2026 research topic in laser powder bed fusion (LPBF): segment the layer image, rebuild the solidified-region contour, and compare it with the nominal CAD slice. Commercial monitoring (EOS EOSTATE Exposure OT) works on raster images and per-region thresholds rather than polygons. Multi-material voxel printing (Stratasys J750 Voxel Print) takes stacks of bitmaps as input. So both "raster in, polygon out" (monitoring) and "polygon in, raster out" (slicing to bitmaps) appear, and the 3MF Slice extension is the standard container for per-layer polygon slices.

### Cited Findings
- In LPBF in-situ inspection, powder-bed images are captured layer by layer after laser exposure. Solidified-layer contours are rebuilt by image segmentation and compared against the nominal shape to find geometric deviations. Studies change nominal CAD layer masks to inject controlled deviations between intended and reconstructed geometry. — [Graph-augmented Segmentation of Complex Shapes in LPBF (ResearchGate, 2025)](https://www.researchgate.net/publication/404249497_Graph-augmented_Segmentation_of_Complex_Shapes_in_Laser_Powder_bed_Fusion_for_Enhanced_In_Situ_Inspection)
- A 2025 Additive Manufacturing paper uses active contours driven by multiple energy terms for "high precision 3D contour detection" in layerwise LPBF monitoring. Here contour extraction is the main research problem, not an afterthought. — [ScienceDirect S2214860425001253](https://www.sciencedirect.com/science/article/abs/pii/S2214860425001253)
- 2025 work measures how well geometric defects are caught (probability of detection by defect size) across in-situ sensor architectures. — [Progress in Additive Manufacturing, 2025 (s40964-025-01235-w)](https://link.springer.com/article/10.1007/s40964-025-01235-w) (paywalled; I could not get the numbers)
- 2026: Fringe Projection Profilometry built into a commercial LPBF machine gives continuous, time-resolved in-situ "4D" monitoring of the powder bed and build surface during contour and infill scanning. — [Int. J. Adv. Manuf. Tech. 2026](https://link.springer.com/article/10.1007/s00170-026-19042-0)
- Other 2024-2026 layer-image work: real-time powder-bed defect detection with deep learning on 3D point clouds ([Virtual & Physical Prototyping 2024/25](https://www.tandfonline.com/doi/full/10.1080/17452759.2024.2449171)), multi-scale defect detection from in-situ visual monitoring ([Additive Manufacturing 2024](https://www.sciencedirect.com/science/article/abs/pii/S2214860424006225)), and wavelet-feature-fusion online powder-bed defect monitoring ([J. Intell. Manuf. 2026](https://link.springer.com/article/10.1007/s10845-026-02852-0)).
- EOS EOSTATE Exposure OT, described as the first commercial optical tomography system, records the whole build platform at high frequency in near-infrared on the EOS M 290. Areas outside a user-set "normal range" are flagged. It was developed with MTU Aero Engines to cut downstream CT inspection costs by rejecting suspect parts earlier. — [Metal AM](https://www.metal-am.com/?p=63196); [Develop3D](https://develop3d.com/develop3d-blog/eos-adds-camera-based-monitoring-software-for-3d-printing-quality-assurance/); [TCT Magazine](https://www.tctmagazine.com/eos-expansion-eostate-exposure-ot-mtu-aero-engines/)
- Stratasys GrabCAD Voxel Print (add-on for the J750 PolyJet) lets users set the material/colour of every voxel. The workflow slices the model into stacks of .BMP/.PNG images in an outside tool such as MATLAB, then packs them into a GrabCAD Voxel File (*.GCVF). — [GrabCAD Help: Guide to Voxel Printing](https://help.grabcad.com/article/230-guide-to-voxel-printing)
- The 3MF Slice Extension standardizes pre-sliced layer data. HP said it needed both the Slice and Production extensions for its Jet Fusion printers, for industrial scale and part traceability. — [TCT Magazine](https://www.tctmagazine.com/3mf-two-extensions-3d-printing-file/); [3DPrint.com](https://3dprint.com/144590/3mf-consortium-extensions)

### Inferences
- A natural closed-loop AM pipeline: OT/NIR or visible layer image -> threshold or segmentation mask -> **closed polygon rings in machine (mm) coordinates via the pluggable transform (camera homography / galvo calibration)** -> polygon Boolean XOR against the 3MF/CLI slice polygon -> deviation area and Hausdorff distance per part per layer. A polygon representation lets you diff by part ID and keep per-layer records in the same vector format the slicer already uses (3MF slice / CLI). Raster-only comparison cannot do this as cleanly.
- The engine's halo tiling fits large build plates (e.g., multi-laser 400-600 mm plates imaged at tens of microns/pixel give images of many tens of megapixels per layer). The layer recoat interval of a few seconds to tens of seconds sets the latency budget. These numbers are my estimates, not sourced.
- Reverse direction (voxel/bitmap -> contour): converting GrabCAD Voxel Print-style or volumetric-printing bitmap stacks back into 3MF slice polygons, for archiving, compression, or driving contour-scan lasers, is a format-bridging use. The bitmap stacks are huge, and contour slices compress them by orders of magnitude.

### Gaps
- No commercial product found (EOS, Materialise, Additive Assurance, Nikon SLM, Velo3D) that clearly exports in-situ melt-region **polygons** and diffs them against CAD in real time. Commercial tools seem to stay raster/statistical. This looks like a real product gap, but I could not confirm it from vendor docs within budget.
- Did not find per-layer latency budgets or image resolutions for commercial OT systems.
- Did not research volumetric (tomographic) printing or HP MJF internal bitmap pipelines in depth.

---

## Q2. Cutting & nesting: leather hides, stone slabs, textiles, sheet metal

### Takeaway
Scanning natural materials (leather hides, stone slabs) and turning defect regions into keep-out zones for nesting is a mature, commercially active workflow (Mindhive + Zünd, Lectra, SlabWise, Slabsmith, SeKon). Its core primitive is defect raster -> polygon -> nesting constraint. Throughput is about 15 s per hide at 0.17 mm precision.

### Cited Findings
- Mindhive FinishSelect: AI finished-leather inspection that detects and maps defects at **0.17 mm precision in about 15 seconds per hide**. It recognizes 20+ defect types with customer-specific "Natural Feature Guide" grading rules and processes **1,000-1,920 hides per 8-hour shift**. Its defect data feeds **MindCUT**, nesting software for Zünd cutters with multi-hide nesting. — [Zünd Dectura page / search summary](https://www.zund.com/en/cutting-systems/registration-methods/dectura); [DSI offline leather solution PDF](https://dsinternational.com/cdn/shop/files/OFFLINE_LEATHER_SOLUTION_ENG.pdf)
- Zünd markets Dectura as moving leather processing "from subjective, manual operations to a digital, data-driven workflow." — [Zünd Dectura](https://www.zund.com/en/cutting-systems/registration-methods/dectura)
- Lectra (which now owns Gerber, e.g., the Gerber Taurus leather cutter) sells leather cutting and nesting focused on cutting hide waste. — [Lectra Gerber Taurus](https://lectra.com/en/products/gerber-taurus); [Lectra/Zolano case](https://lectra.com/en/library/zolano-overcomes-furniture-production-challenges-during-the-covid-19-crisis-by-increasing)
- Stone: SlabWise offers AI slab nesting, defect mapping (defects to avoid when cutting), and DXF middleware for CNC. Slabsmith has shops photograph raw slabs at full resolution and map jobs onto the photos. SeKon TOP digitizes slab inventory and auto-nests for bridge saws and waterjets. — [SlabWise guide](https://slabwise.com/guides/what-is-slabsmith-software-used-for-in-stone-fabrication); [SlabWise digital fabrication](https://slabwise.com/guides/digital-stone-fabrication); [Stone World on SeKon](https://stoneworld.com/articles/84538-global-e-sourcing-inc-1)
- Research: a vision system detects marble surface defects semi-automatically and nests workpieces around them as an irregular packing problem. — [IP Paris: Smart marble cutting through a vision system](https://researchportal.ip-paris.fr/en/publications/smart-marble-cutting-through-a-vision-system/)
- Nesting around natural-leather defects is old art (US5258917A, "Method for nesting contours to be cut out of natural leather"). — [PatentHub](https://www.patenthub.cn/cpc/patent-120165-US5258917A-001c36866b477a32df96c1d62fd5f306.html)

### Inferences
- Every step after scanning (hide outline, defect zones by quality grade, slab edge, veining/fissure keep-outs) needs closed polygons in cutter-table coordinates. The pluggable transform would map camera pixels to table mm, including lens/perspective correction. Hides have holes (defect islands inside the hide outline), so correct ring orientation and nesting of holes matters. A 2x2-classification tracer that handles holes and diagonal (saddle) cases cleanly is a real advantage here.
- At 15 s/hide the throughput bottleneck is probably not polygonization. The pain point is more likely **precision and topology** (sub-mm boundaries, multi-grade nested zones) and **direct DXF/nesting output**, not raw speed. That makes leather/stone a weaker "GPU speed" story but a good "correctness and topology" story.

### Gaps
- Did not find how Mindhive or Lectra represent defect zones internally (polygon vs raster) or any stated vectorization latency.
- Sheet-metal "scan the remnant and re-nest" (skeleton reuse), apparel pattern digitization from photos, and glass inspection-to-cut were not researched.

---

## Q3. Industrial vision & quality inspection (battery electrodes, solar EL, steel surface)

### Takeaway
Roll-to-roll battery electrode coating inspection runs at **120-400 m/min** with line-scan cameras over the full web width and produces continuous "roll maps" of defects. That is effectively an unbounded strip image, which matches the halo-tiled scanline design well. Solar EL inspection is mostly cell-level segmentation and classification; polygon output for reporting is not well documented.

### Cited Findings
- Electrode coating inspection line speeds are typically **120-400 m/min**. High-speed line-scan cameras look for pinholes, bubbles, foreign material, agglomerates, edge defects and thickness variation over the full web width. Systems keep "roll maps and records" so bad sections stay out of good cells. — [Roboflow: Battery Electrode Coating Inspection](https://roboflow.com/ai/battery-electrode-coating-inspection); [iFactory battery cell surface defects](https://ifactoryapp.com/industries/automotive-manufacturing/battery-cell-surface-defect-ai-electrode-tab-seal) (vendor marketing pages; moderate reliability)
- ISRA Vision publishes inline battery-production inspection material (coating, calendering, slitting). — [ISRA Vision PDF](https://isravision.com/content/dam/brands/isravision/markets/battery/documents/inline-quality-inspection-battery-production-en.pdf.coredownload.pdf)
- Vieworks (line-scan/X-ray camera maker) and Robovision (Battery Show NA 2025) target battery inspection. — [Vieworks news](https://vieworks.com/en/board/news_view/17544); [Robovision Battery Show NA 2025](https://robovision.ai/resources/events/robovision-battery-show-na-2025)
- Solar EL: cell segmentation in module EL images reaches a median weighted Jaccard of 94.47% and F1 97.62% (FAU, Deitsch et al.). — [FAU open repository](https://open.fau.de/items/a3bd3a76-6a5d-4e85-81f6-5f560a75f493/full); [arXiv 1806.06530](https://arxiv.org/abs/1806.06530)
- A claim of ">3,600 cells per hour" and "85-95% detection accuracy, FPR <0.5%" for solar inline metrology comes from a patent-analytics aggregator, **low reliability**, possibly understated for modern lines. — [GreyB](https://xray.greyb.com/solar-cells/defect-metrology-passivation-ai)

### Inferences
- Rough data rate: a web 1-1.5 m wide at ~20-50 µm/pixel cross-web is ~30k-75k pixels per line. At 400 m/min (6.7 m/s) and 50 µm along-web, that is ~133k lines/s, or roughly 4-10 Gpixel/s of binary defect mask. This is my estimate, not sourced. A one-pass scanline tracer with halo tiles that streams along the web (polygons crossing tile seams get stitched) is what a continuous roll map needs: per-defect polygons in roll coordinates (meters along the coil, mm across) via the pluggable transform, for the MES/roll map and for downstream slitting/cutting decisions.
- The "streaming infinite strip" property (closing rings that span tile boundaries) is the differentiator. Most libraries assume a whole bounded image.

### Gaps
- No source found for how battery-inspection vendors encode defects (bounding boxes vs polygons vs raster roll maps), or for the latency budget between defect detection and the marking/cutting actuator.
- Steel surface (e.g., hot-strip mill) and AOI/PCB inspection polygon output not researched.

---

## Q4. Construction / AEC: floor-plan raster->vector, scan-to-BIM

### Takeaway
Floor-plan vectorization is an established research line (Raster-to-Vector, ICCV 2017; CubiCasa5k benchmark; U-Net wall segmentation followed by vectorization), and commercial tools turn scanned plans into BIM. The hard part is semantic and topological (walls as clean, axis-aligned polygons, rooms as closed regions), not tracing speed.

### Cited Findings
- "Raster-to-Vector: Revisiting Floorplan Transformation" (ICCV 2017) predicts junctions and assembles walls/doors with integer programming. — [CVF Open Access](https://openaccess.thecvf.com/content_iccv_2017/html/Liu_Raster-To-Vector_Revisiting_Floorplan_ICCV_2017_paper.html)
- Pipelines that rebuild semantic 3D models from 2D plans report state of the art on **CubiCasa5k** and generalize across plan styles. — [Papers with Code: Automatic Reconstruction of Semantic 3D Models from 2D Floor Plans](https://astro.paperswithcode.com/paper/automatic-reconstruction-of-semantic-3d)
- U-Net segmentation of walls, then vectorization to wall polygons, from complex floor plans (Universidad de Chile). — [Repositorio UChile](https://repositorio.uchile.cl/handle/2250/196842)
- Commercial: ArchLine.XP documents a "raster to BIM" conversion of floor plans to 3D models. — [ArchLine.XP help](https://help.archlinexp.com/hc/en-us/articles/20856004418449-From-raster-to-BIM-Converting-floor-plans-to-3D-models)

### Inferences
- For floor plans the tracer is a necessary but undifferentiated step. A raw pixel-boundary ring needs Manhattan regularization and simplification afterwards. A pluggable transform from pixels to real-world units, using scale bar or DPI, is useful. Throughput is not the selling point.

### Gaps
- Drone progress monitoring and scan-to-BIM from point clouds were not researched. CubiCasa's own 2025 product status was not checked.

---

## Q5. Energy & utilities: methane plumes, oil spills, power lines, solar layout, wind blades

### Takeaway
Carbon Mapper's operational L3 methane/CO2 plume pipeline (Tanager-1, EMIT) builds **72 candidate binary masks per plume**, combines them into a per-pixel probability, applies thresholds of 0.3 and 0.7 with connected-component selection near the origin, and **ships the result as a raster GeoTIFF, not a polygon**. Vector plume boundaries are therefore a missing output, and the ensemble-of-masks structure suits a GPU batch tracer well.

### Cited Findings
- Carbon Mapper L3/L4 ATBD (v3.0.1, 17 Feb 2026): the L2b concentration map is cropped to a **2,500 m radius** around the plume origin and split into **72 five-degree wedges**. Per-wedge thresholds μ_i + k·σ_all (k=2 for CH4, k=1 for CO2) give **72 candidate masks**. In each, connected pixels above threshold are grouped, clusters must have **≥3 pixels**, and the cluster closest to the origin (within 8 px) is chosen. Per-pixel plume probability p = N/72. An initial mask uses p≥0.3, keeping the connected component nearest the origin, and is refined with p≥0.7. For another sensor a "nearly identical" process dilates by 1 pixel and uses a 16-pixel origin radius. — [Carbon Mapper L3/L4 ATBD PDF](https://assets.carbonmapper.org/documents/L3_L4%20Algorithm%20Theoretical%20Basis%20Document_021726.pdf)
- Same ATBD: output 8.1 is a "**Segmented plume image (GeoTIFF)** that provides the boundaries and concentrations". The plume length used for emissions is taken from the segmented plume's **convex hull (maximum Feret diameter)**. Masking uncertainty is the standard deviation of IME estimates across the dynamic plume masks. — [Carbon Mapper ATBD](https://assets.carbonmapper.org/documents/L3_L4%20Algorithm%20Theoretical%20Basis%20Document_021726.pdf)
- QC rejects quantification for heavy retrieval artifacts, non-standard shapes (wind shear, pooling, **large gaps in plume** that make segmentation hard), plumes at image edges, and **overlapping plumes**. — [Carbon Mapper ATBD](https://assets.carbonmapper.org/documents/L3_L4%20Algorithm%20Theoretical%20Basis%20Document_021726.pdf)
- Tanager-1 has 30 m spatial resolution and detected plumes throughout 2025 (many catalogued per-plume on GEM.wiki). — [Carbon Mapper: Seeing around the clouds](https://carbonmapper.org/articles/seeing-around-the-clouds); [GEM.wiki example](https://www.gem.wiki/Alabama_State_Methane_Observation_2025-03-02,_2)

### Inferences
- Concrete, current use case: **vectorize all 72 candidate masks (and the 0.3/0.7 probability contours) per plume in one GPU batch**. That gives (a) GeoJSON plume polygons for the public data portal and regulators, (b) ring-based geometric uncertainty (polygon spread across candidates), and (c) exact convex hull and Feret length from the polygon instead of pixel sets. The pluggable coordinate transform maps sensor/pixel coordinates to lat/lon or UTM. "Plumes at image boundaries" are a QC failure today, which is where halo tiling across scene edges is directly relevant.
- Volume: at Tanager scale (thousands of plumes per year, each a small crop), raw speed is not the bottleneck for one sensor. The case gets stronger for (i) ensemble/uncertainty work (72× masks, or Monte Carlo over retrievals) and (ii) multi-constellation scaling (MethaneSAT was not researched).

### Gaps
- Oil spill (SAR) delineation, power-line vegetation encroachment (LiDAR/imagery), solar farm layout, and wind blade inspection had no dedicated searches. No findings.
- Did not confirm whether Carbon Mapper's portal ships any vector plume outlines anywhere (ATBD lists only GeoTIFF plus point attributes).

---

## Q6. Weather & scientific computing: NWS warning polygons, ProbSevere, isocontours

### Takeaway
NWS warnings are polygons, and NOAA's FACETs/Threats-in-Motion program wants warning polygons that **update every minute**, driven by gridded Probabilistic Hazard Information (PHI). ProbSevere already outputs "shapefiles contoured around radar storm cells". In the Python geoscience stack, GPU contouring is explicitly missing: xarray-spatial's contour code says stitching is "inherently sequential" and copies CuPy arrays to the CPU. The user's engine is a counterexample to that claim.

### Cited Findings
- ProbSevere (NOAA/CIMSS) combines GOES, MRMS radar, lightning and NWP to estimate the probability of severe weather up to 90 min ahead. "Model output are shapefiles contoured around radar storm cells". Gradient-boosted trees produce ProbHail, ProbWind and ProbTor. — [NOAA VLab ProbSevere briefing](https://vlab.noaa.gov/documents/10157/137122/ProbSevere_Vlab_Briefing.pdf/fb318646-32ab-6bb6-2213-1f6291fdd8b8); [CIMSS ProbSevere HWT2021](https://cimss.ssec.wisc.edu/severe_conv/psv3/ProbSevere_HWT2021_OnePager.pdf)
- Threats-in-Motion (TIM) would replace static warning polygons with continuously updating polygons that **update every minute** and move with the storm. It is the first stage of FACETs. PHI grids are continuously updating probabilistic hazard grids that can drive user-specific products (e.g., longer lead time at lower confidence for vulnerable populations), and PHI sets the spatial extent of the TIM polygon. — [NOAA GSL TIM](https://gsl.noaa.gov/tim/); [NOAA repository: NWS Severe Weather Warnings as TIM](https://repository.library.noaa.gov/view/noaa/28465); [NOAA GSL PHI news](https://gsl.noaa.gov/news-media/news/probabilistic-hazard-information); [NSSL FACETs](https://inside.nssl.noaa.gov/facets)
- xarray-spatial `contour`: CuPy input is accepted, but "CuPy backend: transfer to CPU and run numpy implementation". The rationale given is "**Contour tracing and segment stitching are inherently sequential graph operations that don't benefit from GPU parallelism**." Segments are stitched by endpoint matching into polylines and closed rings. — [xarray-spatial contour source](https://xarray-spatial.readthedocs.io/en/latest/_modules/xrspatial/contour.html)
- GPU isosurface extraction (marching-cubes family) can beat CPU versions by up to two orders of magnitude. The literature focuses on surface/triangle generation, not ring stitching. — [Efficient and Quality Contouring Algorithms on the GPU (Eurographics)](https://diglib.eg.org/items/c94e2e75-fc6a-4c55-8e77-c6ef52460288/full); [ZCU: Parallel iso-surface extraction](https://dspace5.zcu.cz/handle/11025/26435)

### Inferences
- **Strong candidate use case:** PHI-style probability grids thresholded at several levels (e.g., 10/30/50/70%), refreshed every minute over CONUS MRMS (~7000×3500 at 1 km, my estimate), polygonized into nested, closed probability rings in lat/lon via the pluggable projection transform, every minute and for ensembles (Warn-on-Forecast members). Today this is CPU contouring plus stitching, and the main open-source GPU geoscience library says it cannot be done on GPU.
- The binary-occupancy design fits multi-threshold contouring if each threshold is run as its own bitmap, batched. For smooth isolines, sub-pixel interpolation (marching squares) would be an extension. The current 2x2 classification gives pixel-edge (staircase) rings.

### Gaps
- Did not confirm the operational status of TIM/PHI in 2025-2026 (still testbed vs deployed).
- Did not research climate model output contouring tools (e.g., matplotlib/contourpy, GDAL `gdal_contour`/`gdal_polygonize`) performance numbers or GitHub complaints.
- Fluid simulation level-set interface extraction was not researched.

---

## Q7. Oceanography / polar: sea-ice floe delineation

### Takeaway
Floe size distribution (FSD) needs individual floe polygons at very large scale. One study automatically processed **4,861 MODIS images and identified over 9.4 million floes over 23 years**, and SAM-type foundation models are now used for floe segmentation.

### Cited Findings
- Automated segmentation identified and labeled individual floes in MODIS imagery: **4,861 images, >9.4 million floes over 23 years**. — [The Cryosphere 18, 5031 (2024)](https://tc.copernicus.org/articles/18/5031/2024/)
- Floe segmentation in close-range optical imagery with active contours and foundation models (SAM). Reported as suitable for real-time, scalable analysis with close agreement to benchmark FSDs. — [arXiv 2409.06641](https://arxiv.org/html/2409.06641v4)
- FSD is an important metric for characterizing sea ice cover and checking models. Segmenting individual floes is the core challenge. — [The Cryosphere 17, 5519 (2023)](https://tc.copernicus.org/articles/17/5519/2023/)

### Inferences
- Floes touch and merge, so instance masks need to become per-floe polygons with area/perimeter/shape statistics. A GPU tracer that labels and polygonizes millions of touching instances per scene, in polar-stereographic coordinates via the transform, is a scientific-scale batch use case. It overlaps with geospatial EO, which another researcher covers.

### Gaps
- Coral reef mapping, mining (core imaging, ore sorting, rock fragmentation) were not researched.

---

## Q8. Documents & cartography: historical maps, engineering drawings

### Takeaway
Historical map vectorization (Paris atlases benchmark, Sanborn fire-insurance building footprints) is an active 2024-2026 area. Closed-shape extraction there uses learned edge filters plus watershed, then polygon extraction.

### Cited Findings
- Benchmark for automatic vectorization of historical maps (Paris atlases, 19th-early 20th c.): deep edge filters (including ViTs with topology-preserving losses) plus Meyer/learnable watershed for **closed shape extraction**. Data and code are open. — [PLoS ONE 2024 / PMC10868791](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10868791/)
- Scalable open-source workflow extracting building footprints as vector polygons from **Sanborn maps** (Chicago Urban Heritage Project). Notes that earlier ML methods were computationally heavy and urban-focused. — [EPB: Urban Analytics & City Science, DOI 10.1177/23998083261466497](https://www.citedrive.com/en/discovery/computational-mapping-of-historic-built-environments-in-american-cities-from-sanborn-maps-a-case-from-the-chicago-urban-heritage-project/)
- 2025 arXiv pipeline for historical-map processing. — [arXiv 2508.03564](https://arxiv.org/pdf/2508.03564v1.pdf) (not read in detail)

### Inferences
- Watershed basins are a label image, so per-label closed-ring extraction with georeferencing (pixel -> map CRS via the pluggable transform, including polynomial/TPS georef) is exactly the last mile. Corpora like Sanborn (hundreds of thousands of sheets) make batch GPU throughput relevant, but this is a batch/archival use case, not a latency one.

### Gaps
- MapKurator (text spotting) and David Rumsey integration, and engineering drawing/P&ID digitization, were not researched.

---

## Q9. Games / graphics / VFX: sprite colliders, rotoscoping

### Takeaway
Game engines (Unity, Godot) generate polygon colliders from sprite alpha. The known pain is too many vertices and runtime cost, not generation speed. In VFX, 2026 brought spline-native AI roto (Foundry SmartRoto, $599/yr), because artists need **editable vector splines, not raster mattes**. Mask-to-spline conversion is the gap competitors fill afterwards.

### Cited Findings
- Unity 6 documents creating collision geometry (physics shape) for sprites. — [Unity 6000.3 manual](https://docs.unity3d.com/Manual/sprite/create-collision-geometry.html)
- Unity users report auto-generated polygon colliders dropping framerate (~60 to ~30 fps with ~50 sprites). Fixes are fewer points/paths or simpler primitive colliders. — [Unity Discussions](https://discussions.unity.com/t/improving-2d-collision-detection-performance-with-polygoncolliders/594670); asset-store tool [Collider2D Optimization](https://assetstore.unity.com/packages/tools/physics/collider2d-optimization-63433)
- Godot asset-store tools generate collision polygons from sprites (e.g., "Godot Sprite Collision Generator", "Breakable 2D Sprites"). — [Godot Asset Store: sprite collision generator](https://store.godotengine.org/asset/zwickarr/godot-sprite-collision-generator/); [Breakable 2D sprites](https://store.godotengine.org/asset/the-tech-umbrella/breakable-2d-sprites/)
- **Foundry SmartRoto** (July 2026): AI plugin for Nuke that moves and deforms existing roto splines across shots, "up to 4x faster". Unlike Silhouette/Mocha Pro, which produce AI mattes (convertible to splines afterwards), it works natively on splines so output stays editable. Runs locally. **$599/year** separate subscription, criticized for not being bundled into NukeX. — [CG Channel, 21 Jul 2026](https://www.cgchannel.com/2026/07/foundry-releases-new-ai-powered-smartroto-add-on-for-nuke/)
- Nukepedia "Smart Roto Simplifier" reduces dense, jittery roto data into clean editable splines. This is evidence that mask->dense-polygon->spline cleanup is a real workflow pain. — [Nukepedia](https://www.nukepedia.com/tools/python/misc/smart-roto-simplifier/)

### Inferences
- Destructible 2D terrain or breakable sprites (the Godot "Breakable 2D Sprites" type) need **runtime** re-polygonization of alpha masks every frame after damage. A GPU tracer that keeps the mask on GPU and emits collider rings without a CPU round-trip is a plausible niche, though the market is small.
- VFX: SAM 2 / matting output -> closed rings -> spline fitting, done per frame over 4K sequences with consistent topology across frames. The tracer gives the raw ring. The value-added gap is **temporally coherent point correspondence and spline fitting**, which a pixel tracer alone does not solve.

### Gaps
- SDF font generation (msdfgen etc.) and Adobe/Runway mask-to-shape features not researched.

---

### Cross-cutting summary of gaps/pain points (for the report writer)
1. **"GPU can't stitch contours"**: xarray-spatial says this explicitly and falls back to CPU ([source](https://xarray-spatial.readthedocs.io/en/latest/_modules/xrspatial/contour.html)). A one-pass GPU ring builder directly answers it, with a clear hook into the geoscience/weather Python ecosystem.
2. **Operational pipelines stop at rasters**: Carbon Mapper ships GeoTIFF masks plus convex-hull metrics, not polygons ([ATBD](https://assets.carbonmapper.org/documents/L3_L4%20Algorithm%20Theoretical%20Basis%20Document_021726.pdf)). EOS OT flags raster regions ([Metal AM](https://www.metal-am.com/?p=63196)).
3. **Ensembles and minute-scale refresh**: 72-mask plume ensembles, minute-refresh TIM/PHI polygons ([NOAA GSL](https://gsl.noaa.gov/tim/)), and per-layer AM monitoring are the latency- and batch-heavy settings where GPU polygonization pays off.
4. **Unbounded strips and tile seams**: roll-to-roll inspection at 120-400 m/min ([Roboflow](https://roboflow.com/ai/battery-electrode-coating-inspection)), and plumes "at image boundaries" rejected by QC, are where halo-tiled seam-correct ring closure matters.
5. **Downstream needs simplification/splines**: games (vertex count kills perf), VFX (editable splines), floor plans (Manhattan regularization). A raw pixel-edge ring is necessary but not sufficient, so pairing the tracer with GPU simplification is the likely product gap.
