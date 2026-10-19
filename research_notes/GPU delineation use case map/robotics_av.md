# Robotics, Autonomous Vehicles, Drones and Spatial/Embodied AI: where real-time raster/occupancy-grid to polygon conversion matters (2023–2026)

Scope note: the engine under study takes a binary occupancy bitmap (tiled with halo buffers) and does one scanline pass with 2x2 pixel-type classification. It outputs closed polygon rings (outer boundaries and holes) and supports a pluggable coordinate transform. This branch looks at robotics, AV, UAV and spatial AI. About 23 searches and fetches went into it. Where a claim could not be confirmed from a source, it is listed under Gaps.

---

## Q1. Occupancy grids and BEV/3D occupancy networks: are outputs turned into polygons or free-space contours for planners? Which planners use polygons (convex decomposition for MPC, IRIS/GCS, polygonal obstacle avoidance)? What are the latency budgets?

### Takeaway
Occupancy networks (Tesla, Occ3D-style) output voxels or grids. Planners that use convex or polygonal structure (MPC with half-space constraints, IRIS/FIRI safe corridors, Graphs of Convex Sets) need polygons or half-spaces. The step that bridges the two (grid to contours to convex pieces) is usually an unpublished, ad-hoc CPU step, often built on OpenCV findContours or straight-skeleton decomposition. Recent work keeps shifting this convex-set generation onto the GPU, but it starts from point or obstacle primitives, not from a fast raster-to-polygon front end.

### Cited Findings
- **Occ3D** set up the Occ3D-Waymo and Occ3D-nuScenes benchmarks for 3D semantic occupancy prediction. The output is a voxel occupancy and semantics grid, not polygons. — [Occ3D, NeurIPS 2023 D&B](https://neurips.cc/virtual/2023/poster/73670)
- **Tesla occupancy network** was first shipped in FSD Beta 10.69. It predicts per-voxel occupancy from 8 cameras, separates static from dynamic objects, and is reported to run at ">100 FPS (three times the camera)". This comes from an aggregator write-up of Tesla AI Day 2022, not a primary source. — [42how, Tesla FSD occupancy network explained](https://en.42how.com/2022/10/26/tesla-fsd-occupancy-network-explained/)
- **OccGen** (generative 3D occupancy) says its latency is "comparable to previous state-of-the-art methods and has not achieved a significant speed advantage". Occupancy nets are already heavy, so any post-processing has to be cheap. — [OccGen arXiv 2404.15014](https://arxiv.org/pdf/2404.15014)
- **Volvo Cars patent US12230018B2** ("Free space estimator for autonomous movement", priority June 2022, published Feb 18, 2025) has DNNs predict free-space polygon vertices and vertex distances directly from camera and ultrasonic data. A sequential Kalman filter fuses them into a "virtual polygon", which is sent to the trajectory planning module. The patent presents this as a "lightweight" alternative to "heavyweight occupancy grid determinations" for AEB and low-speed parking. It gives no latency figures. — [Google Patents US12230018B2](https://patents.google.com/patent/US12230018)
- **Deformable Radar Polygon** (2022) is a polygon-based occupancy representation in which "the interior signifies free space and the exterior represents undrivable areas". It uses reduced-vertex boundary polygons as "compact, easy-to-process" contour objects for short-range collision avoidance. — [arXiv 2203.01442](https://ar5iv.labs.arxiv.org/html/2203.01442)
- **PolySAP (2025, Engineering Applications of AI)** says the hard part of MPC planning is collision avoidance on "large and complicated grid maps". It approximates grid-map obstacles with polygons using a Straight-Skeleton-based decomposition, convexifies them, and builds an artificial potential field from them. In other words, it is an explicit grid to polygon to convex pipeline in front of MPC. — [ACM DL / EAAI 2025, 110690](https://dl.acm.org/doi/10.1016/j.engappai.2025.110690)
- **Polygonal Obstacle Avoidance combining MPC and Fuzzy Logic** (Schröder, Schöneberg, Görges, Schotten; arXiv July 2025) re-formulates "occupancy grid maps to continuously differentiable functions" for MPC. Each obstacle is a polygon written as an intersection of half-spaces, and fuzzy logic combines the edge inequalities. Only simulation results; no runtimes in the abstract. — [arXiv 2507.10310](https://arxiv.org/abs/2507.10310)
- **IRIS** (Iterative Regional Inflation by Semidefinite Programming) grows maximum-volume polytopes around convex obstacles from seed points. Its "considerable computational overhead in complex environments limits its real-time applications". — [ST-GCS arXiv 2503.00583](https://arxiv.org/abs/2503.00583)
- **FIRI** (Fast Iterative Region Inflation, Zhejiang Univ., Wang et al.) claims "orders of magnitude improvement compared to IRIS". It works directly on polytope-type obstacles ("as opposed to point-only methods that require discretization") and outputs H-representation polytopes. For voxel maps it notes that Gao et al.'s Parallel Convex Cluster Inflation is "near real-time only at low map resolution". — [FIRI arXiv 2403.02977](https://arxiv.org/html/2403.02977)
- **Space-Time Graphs of Convex Sets (ST-GCS)**, IROS 2025, covers collision-free space-time with convex sets instead of random sampling. It adds Exact Convex Decomposition to reserve the trajectories of other robots as space-time obstacles. — [arXiv 2503.00583](https://arxiv.org/abs/2503.00583); [project site](https://sites.google.com/view/stgcs)
- **Superfast configuration-space convex set computation on GPUs** (Werner, Cheng, Stewart, Tedrake, Rus; arXiv April 2025) inflates collision-free piecewise-linear paths into sequences of convex sets using GPU parallelism. It reports being 17.1x faster with 27.9% higher reliability than a nonlinear trajectory-optimization baseline. — [arXiv 2504.10783](https://arxiv.org/abs/2504.10783)
- **STAR-Filter** (Wu, Zhao, Ong, Kumar; arXiv April 2026) uses star-shaped sets as a fast filter for convex free-space (safe flight corridor) generation on "real-world noisy and large-scale data". It claims the lowest computation time among the methods compared. — [arXiv 2604.26626](https://arxiv.org/abs/2604.26626)
- Other polytope-planning work from 2026 includes "Safe Polytope-in-Polytope Motion Planning and Control with Control Barrier Functions" ([arXiv 2606.09719](https://arxiv.org/pdf/2606.09719)) and "Search-Based Spatiotemporal and Multi-Robot Motion Planning on Graphs of Space-Time Convex Sets" ([arXiv 2607.00444](https://arxiv.org/pdf/2607.00444)). The convex/polytope planning line of work is still active as of mid-2026.
- Search snippets show that other MPC and corridor planners use spatiotemporal occupancy grid maps directly, extending kinodynamic A* plus corridor-constrained optimization to "obstacles with arbitrary shapes" without polygonizing them. — [arXiv 2404.15602](https://arxiv.org/html/2404.15602)

### Inferences
- There is a clear gap in the pipeline: **GPU occupancy (nvblox, occupancy nets) → [CPU findContours / skeleton / point sampling] → convex decomposition (IRIS/FIRI/GCS) → MPC**. The middle step has no named, benchmarked GPU solution. A GPU polygonizer that outputs rings plus holes (holes being obstacles inside free space, which is exactly what the "free-space with holes" polygon input to convex decomposition needs) would let the whole pipeline stay on the GPU.
- FIRI accepts polytope obstacles directly and says point-discretization is a weakness. Exact, closed polygon rings from a grid fit FIRI-style and GCS-style inputs better than raw occupied-cell point clouds do: far fewer constraints than one per occupied cell.
- Free-space polygons work as a compact interface. Volvo predicts them directly with DNNs, and radar-polygon work uses them for short-range AEB. A deterministic raster to polygon converter could give a verifiable, non-learned version of the same interface from any occupancy head (useful for safety cases).
- Latency: at 10–100 Hz perception (Tesla reportedly >100 FPS for occupancy), any polygonization step would need to run in about 1 ms on embedded GPUs (Orin/Thor) so it doesn't take a large share of the cycle. This budget is my estimate; I found no source that publishes it.

### Gaps
- I found no public confirmation of whether Tesla, Waymo or Mobileye convert occupancy to polygons inside their planners. Their planner interfaces are proprietary.
- I found no paper that benchmarks the grid-to-polygon step itself (time per frame) inside an IRIS/GCS/MPC pipeline. Papers report convex-set times only.
- I could not confirm the "Tesla occupancy runs every 10 ms" figure from a primary source. Only the aggregator ">100 FPS" claim was found.

---

## Q2. Navigation meshes (Recast/Detour), ROS 2 Nav2 maps, costmap polygons, warehouse AMRs, Isaac ROS / nvblox: is there any GPU polygonization yet?

### Takeaway
The ROS 2/Nav2 ecosystem works mostly in **polygons → raster**: keepout zones, speed zones and collision-monitor zones start as polygons that get rasterized into masks, or are checked as polygons against points. The reverse direction (live raster → polygon) has no standard GPU tool. NVIDIA nvblox builds 3D maps and 2D costmaps on the GPU but does not output polygons. Navigation meshes are getting attention again (ETH's SE(2) NavMesh, July 2026), which brings back the need for fast online polygonal free-space extraction.

### Cited Findings
- **nvblox** (Millane, Oleynikova et al., ICRA 2024) is GPU incremental SDF mapping with "up to a 177x speed-up in surface reconstruction, and up to a 31x improvement in distance field computation" over CPU systems (voxblox). It is aimed at embedded GPUs. — [arXiv 2311.00626](https://arxiv.org/abs/2311.00626)
- **Isaac ROS Nvblox** takes depth cameras and/or 3D LiDAR and computes a 3D reconstruction and **2D costmaps** on the GPU for navigation. The docs list mesh, ESDF and occupancy layers (static and dynamic). I found no polygon or contour output. — [Isaac ROS Nvblox docs](https://nvidia-isaac-ros.github.io/repositories_and_packages/isaac_ros_nvblox/index.html); [technical details](https://nvidia-isaac-ros.github.io/v/release-3.1/concepts/scene_reconstruction/nvblox/technical_details.html)
- **isaac_mapping_ros** (Isaac ROS 4.0) builds occupancy grid maps from images with cuVSLAM and nvblox. Map quality depends on the disparity input. — [Isaac ROS map creation tutorial](https://nvidia-isaac-ros.github.io/v/release-4.0/repositories_and_packages/isaac_ros_mapping_and_localization/isaac_mapping_ros/tutorial_map_creation.html)
- **Isaac Sim Occupancy Map Generator** produces 2D occupancy maps from simulated USD scenes, for digital-twin to Nav2 map workflows. — [Isaac Sim docs](https://docs.isaacsim.omniverse.nvidia.com/5.0.0/digital_twin/ext_isaacsim_asset_generator_occupancy_map.html)
- **Nav2 Collision Monitor** is a safety node that bypasses costmaps and planners. It uses polygons relative to the robot base frame to define zones (stop, slowdown, approach), and sensor points falling in a zone trigger actions. Polygons here are hand-configured, not extracted. — [nav2_collision_monitor README (Jazzy)](https://api.nav2.org/nav2-jazzy/html/md_nav2_collision_monitor_README.html)
- **Nav2 Costmap Filters / Keepout Filter**: keepout and speed zones are "mask-based" raster layers. A pixel mask (for example, part of the map blacked out) is loaded as a filter. — [ROS Discourse, keepout zones](https://discourse.openrobotics.org/t/nav2-keep-out-zone-support-now-available-and-meeting-oct-15-cancelled/16633); [KeepoutFilter API](https://api.nav2.org/nav2-jazzy/html/classnav2__costmap__2d_1_1KeepoutFilter.html)
- **Nav2 Vector Object Server** "puts vector objects (such as polygons and circles) on OccupancyGrid raster map". It rasterizes once at startup and serves the result to costmap filters. Uses include dynamic-obstacle highlighting, sensor-noise removal and black-out areas. This is the **vector → raster** direction. — [Nav2 docs, Vector Object Server](https://docs.nav2.org/configuration/packages/map_server/configuring-vector-object-server.html)
- **SE(2) Navigation Mesh** (Shi, …, Marco Hutter; arXiv 2607.01454, July 2026) is a polygonal navmesh with yaw-dependent traversability. It is built from footprint masks over yaw-specific layers and updated incrementally online from streaming point clouds. It captures ">50% more traversable area than classical NavMeshes", uses A*-String-Pulling-A* pathfinding, and was validated on physical robots. It names the drawbacks it addresses: point clouds and volumetric maps lack surface structure, and dense triangle meshes are "computationally prohibitive" for pathfinding. — [arXiv 2607.01454](https://arxiv.org/abs/2607.01454)

### Inferences
- The SE(2) NavMesh's per-yaw "footprint mask" layers are binary rasters (traversable for heading θ) that become polygonal regions. Turning N yaw-layer masks into polygons, incrementally and tile by tile, closely matches a halo-tiled GPU raster → ring engine. This is a concrete, current, research-grade use case.
- Recast/Detour builds navmeshes by voxelizing, finding regions and tracing contours, then simplifying them into polygons. It runs on the CPU, so rebuilding tiles online on robots is the expensive part. A GPU contour stage could replace Recast's contour-tracing step for dynamic tile updates. (This is inferred from general knowledge of Recast's pipeline; I found no 2023–2026 source benchmarking it.)
- Nav2 has a full polygon → raster toolchain but no raster → polygon one. Live use cases with no tool include: exporting learned or observed keepout areas as polygons, building collision-monitor zones from nvblox dynamic-obstacle layers, and handing fleet managers polygonal occupancy summaries, which use much less bandwidth than OccupancyGrid messages.

### Gaps
- I could not find a Recast/Detour GPU port or GPU navmesh-tile rebuild benchmark from 2023–2026.
- I found no warehouse AMR vendor (Locus, OTTO/Rockwell, Geek+, etc.) documentation describing runtime grid → polygon conversion.
- I did not get exact nvblox per-frame latencies on Jetson Orin from the fetched pages.

---

## Q3. Drones/UAV: real-time field mapping, precision-ag spraying zones (prescription maps), SAR, wildfire perimeter mapping from aircraft: latency, and how perimeters are produced

### Takeaway
Wildfire perimeter mapping is the strongest, best-documented case of **raster → polygon with painful latency**. Airborne IR line-scanner imagery is still interpreted by humans into perimeter shapefiles, and those reach crews hours later (up to 12 h). Onboard ML now produces fire masks at >100 FPS, but the mask → georeferenced perimeter polygon step has not been automated or published. Precision-ag spraying uses polygonal "rate zones" made in desktop software (Pix4Dfields) after the flight, not in flight.

### Cited Findings
- **NIROPS** (USDA Forest Service National Infrared Operations) flies three fixed-wing aircraft at night. Their line scanners use rotating mirrors, and successive scan lines are assembled into high-resolution heat imagery from which "highly accurate map products" are produced. — [USGS NCAC 2023 Wildfire Mapping Support (AK Geo Summit)](https://dggs.alaska.gov/webpubs/dggs/ago/documents/2023AKGeoSummit/2023AKGeoSummit_Session14_Ganuza.pdf)
- **FireGuard** (satellite/DoD-sourced) products arrive every 10 minutes with 24/7/365 nationwide coverage. FireGuard "shows fire detections and fire spread but is not a perimeter product per se", with about ±100 m location accuracy that needs correction by interpretation. "Aircraft 3" products work like NIROPS products at 5–10 m scale. — [same USGS NCAC deck](https://dggs.alaska.gov/webpubs/dggs/ago/documents/2023AKGeoSummit/2023AKGeoSummit_Session14_Ganuza.pdf)
- Current workflow: "Trained human firespotters then label these images with the location of the fire perimeter and heat sources." That information is "often used to direct response teams, oftentimes up to 12 hours after the original flight". NIROPS delivers "geo-rectified TIFFs and shapefiles … before morning briefings". California's FIRIS provides "near-real-time perimeter and hotspot intelligence". — [Real-Time Wildfire Localization on NASA AMS, arXiv 2601.14475](https://arxiv.org/html/2601.14475v1)
- **NASA AMS onboard DL (arXiv 2601.14475, Jan 2026)** outputs **pixel-level 256×256 binary segmentation masks, not polygons or shapefiles**. Inference runs at ">100 fps on NVIDIA T1200 laptop GPUs" (classification 1.19 ms, segmentation 7.75 ms), processing 200–300 patches per 5–10 minute flight. The authors note it still "requires integration with existing operational workflows". — [arXiv 2601.14475](https://arxiv.org/html/2601.14475v1)
- Other references: the NIROPS/FIRIS local news coverage ([KRCR](https://krcrtv.com/news/local/nirops-and-firis-programs-boost-firefighting-with-advanced-data-and-mapping-tools)); NIFC/USFS Fire Imaging Technologies Users Guide ([NIFC PDF](https://www.nifc.gov/sites/default/files/document-media/Fire%20Imaging%20Technologies%20Users%20Guide.pdf)); the Esri blog on the IAA wildfire intelligence program ([Esri](https://www.esri.com/en-us/industries/blog/articles/iaa-the-leading-wildfire-intelligence-program)); NASA's firefront aerial DAQ report ([NTRS](https://ntrs.nasa.gov/api/citations/20210023767/downloads/DAQ_Report.pdf)); and EGP (wildfire.gov) release notes ([EGP](https://egp.wildfire.gov/egp/about/release-notes/)), the platform where perimeters are published.
- **Precision ag:** PIX4Dfields supports DJI Agras variable-rate and spot-spray missions. Application zones come from automatic detection or manual drawing, with up to 7 rate zones. — [Pix4D blog, PIX4Dfields 2.4](https://pix4d.com/blog/pix4dfields-2-4-revolutionize-precision-agriculture-targeted-operations). PIX4Dfields 2.8 added XAG prescriptions. — [sUAS News](https://www.suasnews.com/2024/10/pix4dfields-2-8-xag-prescriptions-ai-crop-report-and-more/)
- Users report trouble estimating herbicide amounts for DJI Agras T25 spot-spray prescriptions. This points to friction between detected weed areas and the zones the sprayer actually uses. — [Pix4D community forum](https://community.pix4d.com/t/dji-agras-t25-for-spotspraying-trouble-calculatiing-right-amount-of-herbizide-needed/34667)
- DJI Agras T50/T25 do offline field reconstruction and field planning ("Mapping mode") in the Agras app. — [DJI T50/T25 user manual](https://dl.djicdn.com/downloads/t50_t25/20240201/T50_T25_User_Manual_v1.0_EN.pdf)

### Inferences
- **Highest-signal use case in this branch:** onboard or edge conversion of fire-segmentation masks into georeferenced perimeter and hotspot polygons, live, during the flight. The pluggable coordinate transform matters here because scanner geometry and orthorectification differ by sensor. The halo-tiled design also matters because line scanners produce strips, and polygons need to be stitched across strips and tiles. The current baseline is human interpretation with latency up to about 12 h, while ML masks are already >100 FPS. That makes the missing mask → polygon step the bottleneck in an otherwise fast pipeline.
- Precision-ag spot spraying fits a "see → polygonize → spray-zone" loop. Detections need to become zones with buffers and holes (for example, avoid a waterway inside a field). Today this happens after the flight on a desktop. A real-time version on the drone or edge box could enable same-flight scout-and-spray.

### Gaps
- I found no documented latency for FIRIS's "near-real-time" perimeters, and no description of how FIRIS or Overwatch-type vendors vectorize their perimeters (manual digitization vs. automated).
- I found no 2023–2026 paper on automated GPU vectorization of airborne fire masks into perimeters.
- I did not research search-and-rescue coverage polygons, so there are no findings on it.

---

## Q4. Online HD map construction (MapTR, StreamMapNet) and raster → vector

### Takeaway
Online HD mapping has moved away from "segment BEV raster → heuristic polylines" (HDMapNet) to end-to-end vector heads (MapTR/MapTRv2), largely because the raster-to-instance post-processing was slow and heuristic. This area is evidence that naive raster → vector is a pain point. It is also a warning: the field moved to learned vector outputs rather than fixing the post-processing.

### Cited Findings
- HDMapNet "groups pixel-wise semantic segmentation results with heuristic and time-consuming post-processing to generate instances". Rasterized maps lack instance-level structure (direction, keypoints, instance IDs) and handle complex topology poorly. — [MapTR arXiv 2208.14437](https://ar5iv.labs.arxiv.org/html/2208.14437)
- MapTRv2-R18 runs at 33.7 FPS on an RTX 3090, 15x faster than VectorMapNet with +6.3 mAP. MapTRv2-R50 reaches 68.7 mAP at 14.1 FPS. VectorMapNet's autoregressive point decoding is slow. — [MapTRv2 arXiv 2308.05736](https://ar5iv.labs.arxiv.org/html/2308.05736)
- Work is still active in 2026, for example D2HDMap (non-visible driveline map priors for online vectorized HD maps). — [arXiv 2606.20725](https://arxiv.org/pdf/2606.20725)

### Inferences
- HD-map elements are mostly **open polylines** (lane dividers) plus some closed polygons (pedestrian crossings, road boundaries). A closed-ring engine fits closed classes such as drivable area, crosswalks and parking areas better than lane dividers.
- Hybrid opportunity: dense raster heads are still common, for auxiliary losses and occupancy. A very fast deterministic ring extractor could act as (a) a fast vectorization baseline that removes HDMapNet's post-processing bottleneck, (b) a label generator that turns rasterized HD-map GT into polygon targets, or (c) a consistency check between the raster and vector heads.

### Gaps
- I found no measured HDMapNet post-processing latency (ms) in the fetched sources, only the "time-consuming" description.
- StreamMapNet details were not fetched.

---

## Q5. 3DGS/NeRF/scene graphs with SAM masks → object footprints; open-vocabulary 3D mapping (ConceptGraphs, HOV-SG) → room/floor polygons; floor-plan extraction and digital twins

### Takeaway
Open-vocabulary 3D scene graphs (HOV-SG, RSS 2024) segment rooms with classic 2D raster operations: BEV wall histogram, distance transform and watershed. They produce **room masks rather than polygons** and are described as too slow for real-time. Learned floorplan reconstruction (RoomFormer CVPR 2023, PolyRoom ECCV 2024, CAGE 2025) predicts polygons straight from density rasters and reports accuracy, not latency.

### Cited Findings
- **HOV-SG** (Hierarchical Open-Vocabulary 3D Scene Graphs, RSS 2024) builds a floor → room → object hierarchy using OpenCLIP and SAM. It reports a 75% reduction in representation size versus dense open-vocab maps and uses a cross-floor Voronoi graph for multi-story navigation. — [arXiv 2403.17846](https://arxiv.org/html/2403.17846v2); [RSS 2024](https://roboticsconference.org/2024/program/papers/77/); [project](https://hovsg.github.io/)
- HOV-SG room segmentation: a per-floor BEV histogram gives a binary wall mask. The wall mask is dilated, a Euclidean distance field is computed and thresholded into isolated seed regions, and watershed then gives "2D region masks". Floors are found from height histograms (0.01 m bins) plus DBSCAN. The authors say the method "is time-consuming, rendering the method unsuitable for real-time mapping". — [arXiv 2403.17846v2](https://arxiv.org/html/2403.17846v2)
- **RoomFormer** (CVPR 2023) represents a floorplan as variable-size polygons with ordered corners, using two-level queries over a point-cloud density map. — [CVPR 2023 paper](https://openaccess.thecvf.com/content/CVPR2023/papers/Yue_Connecting_the_Dots_Floorplan_Reconstruction_Using_Two_Level_Queries_CVPR_2023_paper.pdf)
- **PolyRoom** (ECCV 2024) predicts initial room shapes from the density map and refines them into a vectorized floorplan. It addresses RoomFormer's missing or biased corners. — [arXiv 2407.10439](https://arxiv.org/html/2407.10439v1); [ECCV paper](https://www.ecva.net/papers/eccv_2024/papers_ECCV/papers/06708.pdf)
- **CAGE** (2025): continuity-aware edge network for robust floorplan reconstruction. — [arXiv 2509.15459](https://arxiv.org/html/2509.15459v2)

### Inferences
- In HOV-SG and similar scene-graph systems, room masks are produced but never become polygons. Room polygons with holes (pillars, fixed furniture) would make "is the robot in room X", room adjacency and map compression cheap, and they map directly to digital-twin and BIM formats. Running watershed-mask → ring extraction on the GPU would also address part of the "not real-time" complaint.
- Object footprints: projecting SAM or 3DGS instance masks to BEV and polygonizing them gives 2D footprints for semantic costmaps ("avoid the rug", "keep out of the kitchen") and for embodied-AI agents that reason over polygons in language ("go to the region left of the sofa"). This is plausible, but I found no paper that explicitly benchmarks this step.

### Gaps
- I did not fetch ConceptGraphs details or 3DGS/NeRF-to-footprint pipelines, so I have no sourced findings there.
- RoomFormer and PolyRoom inference latency was not found.
- I found no sourced Matterport or commercial digital-twin floorplan-extraction pipeline details.

---

## Q6. Collision checking, swept volumes and polygon-based safety fields (safety laser scanners)

### Takeaway
Industrial safety (AMR/AGV) runs on **configured polygonal protective fields** in certified scanners, with up to 128 fields and speed-based switching. These polygons are designed offline. Dynamic, map-aware or swept-path-derived field generation is the frontier, but certification limits what can be computed at runtime.

### Cited Findings
- SICK microScan3 Pro I/O supports 128 simultaneously configurable fields and 128 monitoring cases. — [SICK safety laser scanners](https://cdn.sickcn.com/ca/en/c/g187241)
- SICK S300 Expert supports 16 field sets (48 fields) and 32 monitoring cases, with incremental-encoder inputs for **speed-dependent field switching** ("the correct protective field at any speed avoids unnecessary stops"), a 270° scan angle and 0.5° resolution. — [RS Online S300 listing](https://sg.rs-online.com/web/p/light-detection-and-ranging/2856808); [plc-direct S30B-3011GB](https://plc-direct.com/en-de/products/s30b-3011gb)
- Nav2 Collision Monitor mirrors this pattern in software, with polygon stop, slowdown and approach zones around the robot base. — [nav2_collision_monitor](https://api.nav2.org/nav2-jazzy/html/md_nav2_collision_monitor_README.html)

### Inferences
- A non-safety-rated but useful use case: automatically generating the field polygon sets offline from rasterized swept-volume footprints (robot footprint × braking distance × speed × turn rate). The robot's swept footprint is rasterized for each speed/curvature bin, then polygonized and simplified into scanner field polygons. The engine would serve as a design and commissioning tool, not a runtime safety function.
- Runtime collision-checking stacks often already use rasters (costmaps, ESDF), so polygonization only pays off where downstream consumers need polygons (scanners, fleet managers, MPC half-spaces).

### Gaps
- I found no 2023–2026 source on automated generation of safety-scanner field polygons from maps or swept volumes.
- Field-switching latency figures were not found.

---

## Q7. Papers on GPU-accelerated contour extraction and polygonization for robotics

### Takeaway
GPU contour tracing is still missing from the standard toolchain: **OpenCV CUDA has no findContours**, which forces GPU → CPU → GPU round trips in robot vision pipelines. GPU connected-component labeling is mature (block-based algorithms win per a 2026 TPDS review). Block-parallel GPU contour detection is only now appearing in the literature (2026), and none of it is robotics-specific.

### Cited Findings
- OpenCV CUDA has no GPU counterpart to findContours, and no contributor is currently implementing one. Users describe having to move data from the GPU (after undistortion or Canny) to the CPU for findContours and back again, which costs significant time. — [OpenCV Q&A: Why no CUDA findContours](https://answers.opencv.org/question/181859/why-no-cuda-version-of-findcontours/); [OpenCV Q&A 84942](https://answers.opencv.org/question/84942/isnt-there-a-opencv-cuda-function-similar-to-findcontours); [NVIDIA forum, CUDA blob detection](https://forums.developer.nvidia.com/t/cuda-support-in-opencv-blob-detection/78564)
- A 2026 Journal of Supercomputing article (s11227-026-08488-4) describes block-parallel CUDA contour detection. The image is split into rectangular blocks, structures of interconnected points are extracted per block in parallel, and "triads" are connected within and between blocks to form closed contours. The data stays on the GPU. (Search-snippet level only; the full text was behind a paywall redirect, so I could not verify speedups or hole and hierarchy handling.) — [Springer s11227-026-08488-4](https://link.springer.com/10.1007/s11227-026-08488-4)
- "A State-of-the-Art Review With Code About Connected Components Labeling on GPUs" (IEEE TPDS, April 2026) finds that block-based algorithms beat all other GPU CCL solutions on 2D and 3D. — [UNIMORE IRIS PDF](https://iris.unimore.it/retrieve/a031cde1-6b22-43db-b3c7-45a5a254c791/A_State-of-the-Art_Review_With_Code_About_Connected_Components_Labeling_on_GPUs.pdf)
- "GPU-Accelerated Vision for Robots: Improving System Throughput Using OpenCV and CUDA" (Universitat Jaume I repository) is an example of robotics throughput work built around OpenCV CUDA, which lacks contours. — [UJI repository](https://repositori.uji.es/xmlui/handle/10234/189456)

### Inferences
- The engine's 2x2 classification in one scanline pass, with halo tiles, is in the same family as block-based CCL and the 2026 block-parallel contour work. Its advantages would be emitting **closed rings with an outer/hole distinction** and a **pluggable coordinate transform** (grid → map/world/geodetic), which nothing found here combines. A clear novelty claim for robotics is "the missing cv::cuda::findContours equivalent, emitting topology-correct rings directly in the map frame, for Jetson/Thor".
- **Recommended candidate use cases (ranked by fit × evidence of pain):**
  1. **Airborne wildfire perimeter vectorization** (mask → georeferenced perimeter with holes, tiled across scan strips). The latency gap is documented (up to 12 h, human digitization), and the upstream ML already exceeds 100 FPS.
  2. **GPU occupancy → convex free-space pipeline** (nvblox/occupancy-net grid → rings with holes → FIRI/IRIS/GCS seeds → MPC). Very active research area in 2025–2026, the middle step is unaddressed, and it is attractive to the Tedrake, Kumar and ZJU FAST Lab groups.
  3. **Online navmesh / SE(2) NavMesh tile rebuilds** from per-yaw footprint masks (ETH, July 2026). A direct structural match to halo-tiled raster → polygon.
  4. **Room and object polygons for open-vocab scene graphs** (HOV-SG masks → polygons). Fixes a stated "not real-time" limitation and produces BIM/digital-twin output.
  5. **Nav2 raster → polygon bridge**: the missing inverse of the Vector Object Server. A low-barrier open-source contribution with ecosystem reach.

### Gaps
- I could not access the full text of the 2026 Springer GPU contour paper, so its speedups, hole handling and benchmark sizes are unknown.
- I found no GPU contour or polygonization package in Isaac ROS, cuCIM or VPI documentation during this search. This should be confirmed before claiming "none exists". cuCIM's skimage-compatible API was not checked.
