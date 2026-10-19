# Biomedical Use Cases for GPU Mask-to-Polygon Conversion (2023–2026)

Scope: spatial transcriptomics/proteomics, digital pathology, radiotherapy, connectomics, live-cell/HCS, plus existing GPU tooling and practitioner pain points. The engine under consideration: CUDA/C++ single-pass scanline polygonizer on binary occupancy bitmaps with halo buffers for tiling, 2x2 pixel-type classification, closed rings, pluggable coordinate transform.

Note on method: about 20 search/fetch calls. Several primary pages could not be fetched (10x Xenium zarr doc returned 429, the image.sc thread returned 403, a PMC Sopa page returned 404, the Sopa bioRxiv PDF was too large). Where that happened, the facts below come from search-result snippets and are marked as such.

---

## Q1. Spatial transcriptomics: how cell boundaries are stored and polygonized, scale, and whether polygonization or point-in-polygon is a bottleneck

### Takeaway
Polygons are the standard currency in spatial omics. Xenium ships per-cell vertex arrays (float32, µm), SpatialData stores shapes as GeoPandas GeoDataFrames (GeoParquet), and Visium HD's Space Ranger 4.0 takes in and writes out GeoJSON FeatureCollections. Slides now hold about 0.4–0.85M cells each (about 1.9M across 3 Xenium 5K slides). Every pipeline that segments in tiles (Sopa, Cellpose-based flows, Proseg, Segger) has to (a) vectorize label masks, (b) resolve duplicates and conflicts at patch seams, and (c) assign transcripts to cells. Today those steps run on CPU tools: rasterio, skimage, shapely/JTS, and Dask processes.

### Cited Findings
**Storage formats**
- In SpatialData, polygons must be stored in a GeoPandas GeoDataFrame. Xenium cell boundaries are read from `cell_boundaries.parquet` as polygons by `spatialdata_io.xenium` — [spatialdata-io xenium docs](https://spatialdata.scverse.org/projects/io/en/latest/generated/spatialdata_io.xenium.html); [Bento data prep guide](https://bento-tools.readthedocs.io/en/2.1.1/tutorial_gallery/Data_Prep_Guide.html)
- A SpatialData tutorial example has 167,780 cell-boundary polygons in one GeoDataFrame — [SpatialData squidpy integration](https://spatialdata.scverse.org/en/latest/tutorials/notebooks/notebooks/examples/squidpy_integration.html)
- Xenium Onboard Analysis zarr output stores, per polygon, the vertex count (including a repeat of the first vertex). XY vertex coordinates are float32 in physical µm, and the first vertex is repeated at the end to close the ring. A polygon with 0 vertices means that cell has no polygon. (Search-snippet level, because the page fetch hit HTTP 429.) — [10x XOA output zarr](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/advanced/xoa-output-zarr)
- Visium HD, Space Ranger v4.0, introduced nucleus and cell segmentation. Accepted inputs are a labeled mask (TIFF or NPY), polygons as a GeoJSON FeatureCollection, or a CSV mapping squares to cell IDs. Outputs are `cell_segmentations.geojson` and `nucleus_segmentations.geojson`. `barcode_mappings.parquet` maps 2 µm barcodes to 8 µm and 16 µm bins and to Cell_id — [10x Space Ranger segmentation inputs](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/inputs/segmentation-inputs); [segmented outputs](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/segmented-outputs)
  - The mask-to-GeoJSON round trip is a native part of the official Visium HD workflow. Users bring Cellpose or StarDist masks, and they leave as GeoJSON.

**Scale**
- A Xenium FFPE 5K QC report covers 3 samples with 1,946,621 cells in total: cervical cancer 840,387, ovarian 407,124, breast 699,110 — [U Calgary Xenium 5K QC report](https://asoc.ucalgary.ca/TXG_Xenium/Xenium_FFPE_5k_QCreport.html)
- Sopa is described as scaling to "1TB+ images". The authors say most existing methods are not designed for "large images with millions of cells" and need HPC clusters with large memory — [Sopa PyPI](https://pypi.org/project/sopa); [Sopa, Nature Communications 2024](https://ideas.repec.org/a/nat/natcom/v15y2024i1d10.1038_s41467-024-48981-z.html)
- An image.sc user (thread "Fastest way to import qupath detections from numpy / convert numpy mask into shapely polygons") reports more than 100 images, each with more than 600,000 cells. For large images, converting the numpy array to polygons took "several hours", longer than producing the mask itself. (Search snippet only, because the thread fetch hit 403.) — [image.sc thread 121900](https://forum.image.sc/t/fastest-way-to-import-qupath-detections-from-numpy-convert-numpy-mask-into-shapely-polygons/121900)

**Tiling and seam conflicts**
- Sopa's workflow is: patchify the image for low-memory segmentation, segment (e.g., Cellpose), then "resolve segmentation conflicts at boundaries". It supports Xenium, MERSCOPE, CosMx, PhenoCycler, MACSima and Hyperion — [Sopa PyPI](https://pypi.org/project/sopa); [Vizgen Sopa page](https://vizgen.com/resources/sopa-a-technology-invariant-pipeline-for-analyses-of-image-based-spatial-omics/)
- Sopa GitHub issue #120 reports Xenium Explorer showing both Cellpose and Baysor boundaries, a sign of friction around boundary sets and exports — [sopa issue #120](https://github.com/prism-oncology/sopa/issues/120)

**Vectorization implementations in the scverse ecosystem**
- `vectorize()` (labels to shapes) in the SpatialData and Harpy/Sparrow ecosystem uses rasterio if it is installed and falls back to skimage otherwise. Rasterio is "recommended for increased performance and more precise vectorization", and the docs advise switching Dask to the `processes` scheduler for speed — [harpy.sh.vectorize](https://harpy.readthedocs.io/en/latest/generated/harpy.sh.vectorize.html); [sparrow.sh.vectorize](https://sparrow-pipeline.readthedocs.io/en/latest/generated/sparrow.sh.vectorize.html)
  - Official docs recommend working around the GIL with process-level parallelism, which shows that CPU vectorization is a known throughput limit.

**Transcript-based segmenters and their polygon outputs**
- Proseg (Jones et al., "Cell simulation as cell segmentation", Nature Methods 2025) is voxel-based. Its polygons are axis-aligned, built from horizontal and vertical segments. It writes SpatialData zarr (`sdata.shapes["cell_boundaries"]`), GeoJSON 2D consensus polygons (`--output-cell-polygons`), and non-overlapping per-z-layer polygons (`--output-cell-polygon-layers`). The README warns that earlier Xenium Ranger/Explorer versions "mangle polygons generated by Proseg" because of mis-rendering. Running out of memory is a documented failure mode, mitigated by a coarser voxel size or fewer layers — [proseg GitHub](https://github.com/dcjones/proseg); [Nature Methods 2025](https://www.nature.com/articles/s41592-025-02697-0)
- Proseg issue #165 is titled "Ultra jagged masks with empty spaces and…", a user complaint about polygon and mask quality from voxel-based output — [proseg issue #165](https://github.com/dcjones/proseg/issues/165)
- Segger (bioRxiv 2025) is a heterogeneous GNN that treats segmentation as transcript-to-cell link prediction. It reports a 30x speedup over Baysor and 3x over BIDCell, with multi-GPU support (benchmarks on 4x A100). The preprint does not say which algorithm builds the cell boundary polygons from the assigned transcripts (credited only as "G.R. implemented the cell mask boundary generation algorithm") — [Segger PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11952575/)
- Baysor, JSTA and Proseg form a "hybrid" class that combines transcripts with image priors. FastReseg (Sci Rep 2025) revises boundaries only where local transcript evidence supports a change — [FastReseg, Sci Rep 2025](https://www.nature.com/articles/s41598-025-08733-5)
- A 2026 arXiv review, "The Challenge of Cell Segmentation in Spatially Resolved Transcriptomics", and a Nature Methods 2025 commentary, "Confronting the challenge of cell segmentation in spatial transcriptomics", both treat segmentation as the main open problem — [arXiv 2606.09675](https://arxiv.org/pdf/2606.09675)

### Inferences
- Opportunity 1, a GPU `vectorize()` drop-in. SpatialData, Sopa, Harpy and Sparrow all funnel through one labels-to-shapes function backed by rasterio or skimage. A CUDA polygonizer returning a GeoArrow/GeoParquet-compatible GeoDataFrame could replace a multi-hour CPU step. The pluggable coordinate transform fits SpatialData's affine coordinate systems (pixel to µm) directly.
- Opportunity 2, seam conflicts solved by halo design. Sopa-style "resolve conflicts at boundaries" is a post-hoc fix that works on polygons with shapely and STRtree. An engine with halo buffers that emits only rings owned by a tile (for example, by an anchor-pixel rule) could remove the conflict-resolution stage.
- Opportunity 3, transcript-assignment-native segmenters (Segger, Proseg, Baysor) emit assignments or voxels rather than masks. Rasterizing those into per-cell occupancy and polygonizing on GPU gives clean, non-jagged polygons. Proseg's axis-aligned voxel output is exactly a binary occupancy bitmap. Issue #165 suggests users want smoother or cleaner rings. Corner-cutting or subpixel options in the 2x2 classifier could serve this.
- Point-in-polygon (transcript to cell) at 10^8–10^9 transcripts per slide is a natural adjacent kernel. If labels already exist, a raster lookup is cheaper than PIP. PIP matters when polygons are the source of truth (re-segmentation, manual edits, GeoJSON inputs). This is inference; I found no benchmark quantifying it.

### Gaps
- I could not confirm Xenium's internal mask-to-polygon algorithm, whether vertex counts are capped or simplified, or transcript counts per Xenium 5K slide. The 10x page returned 429.
- I could not fetch Sopa's exact conflict-resolution rule (overlap threshold) or runtime and memory numbers.
- No vendor details were found for Vizgen MERSCOPE VPT or CosMx (AtoMx) polygon storage.
- I found no published benchmark that isolates "polygonization time" or "transcript point-in-polygon time" as a share of a spatial-omics pipeline.
- Bo, ComSeg and BIDCell polygon outputs were not researched.

---

## Q2. Digital pathology WSIs: tile-seam stitching of masks into polygons, memory and time pain points

### Takeaway
Gigapixel WSI nucleus segmentation (HoVer-Net, CellViT, StarDist, HoVer-NeXt) is tile-based. Seam handling uses heuristic polygon merging (STRtree overlap merging), which the 2026 LSP-DETR paper says "frequently split[s] or duplicate[s] cells". Throughput ranges from 0.45 to 13 s/mm². QuPath users report that turning numpy masks into objects takes hours, and the BIOP Cellpose extension had to drop JTS-based label reading because of large-tile problems.

### Cited Findings
- Tile-based WSI processing needs overlapping margins. "A nucleus split across two tiles yields conflicting predictions that must be reconciled by custom heuristic stitching routines, which frequently split or duplicate cells" — [LSP-DETR, arXiv 2601.03163](https://arxiv.org/html/2601.03163v2)
- Measured inference speed on TCGA (LSP-DETR paper): LSP-DETR 0.45 s/mm², StarDist 1.42, HoVer-NeXt 2.38, LKCell 4.39, CellViT 13.24 s/mm². TCGA test slides cover 20–501 mm² — [LSP-DETR](https://arxiv.org/html/2601.03163v2)
  - At CellViT speed, a 500 mm² slide takes about 110 minutes, so post-processing is not the only cost. Still, the paper calls watershed and post-processing pipelines in HoVer-Net and CellViT "computationally expensive post-processing pipelines".
- A transformer WSI detection paper (arXiv 2502.06307) converts contours to polygons and merges tiles with an STRtree-based algorithm. It finds overlapping polygons, merges them, keeps class labels, and computes probability-weighted averages for cells spanning tiles — [arXiv 2502.06307](https://arxiv.org/html/2502.06307v1)
- LazySlide (Nature Methods 2026) is an "accessible and interoperable whole-slide image analysis" framework in the scverse/SpatialData family. It brings WSI results into the same GeoDataFrame-shapes model — [LazySlide, Nat Methods 2026](https://www.nature.com/articles/s41592-026-03044-7)
- QuPath and Python pain points on image.sc. There are recurring threads on importing label masks as QuPath annotations or detections (e.g., threads 56834, 60807, 86823, 72021 "Run cellpose alone then import masks into QuPath", 49410). Thread 121900: more than 600k cells per image and "several hours" to convert masks to polygons. A reply notes that StarDist already returns polygons/polyhedra and that reusing those vertices "would save a lot of time recomputing contours" — [image.sc 121900](https://forum.image.sc/t/fastest-way-to-import-qupath-detections-from-numpy-convert-numpy-mask-into-shapely-polygons/121900); [image.sc 72021](https://forum.image.sc/t/run-cellpose-alone-then-import-segmentation-masks-into-qupath/72021); [image.sc 56834](https://forum.image.sc/t/import-label-mask-as-annnotation-in-qupath/56834)
- BIOP's qupath-extension-cellpose release notes say problems with large tiles and large fields of view came from "the way label images were being read using JTS". The fix switched to LaRoMe's Labels2Rois approach and skipped very small images — [BIOP qupath-extension-cellpose releases](https://github.com/biop/qupath-extension-cellpose/releases)
- QuPath guidance: add objects in batch with `addObjects()` rather than one at a time, because each add triggers costly checks and events — (search snippet from QuPath docs) [QuPath docs 0.7 overview](https://qupath.readthedocs.io/en/0.7/_sources/docs/concepts/overview.md.txt)
- fast-pathology issue #36 relates to label-image and contour handling (snippet only) — [fast-pathology issue #36](https://github.com/aican-research/fast-pathology/issues/36)

### Inferences
- A WSI has about 10^5–10^6 nuclei (consistent with the more than 600k cells per image above). Mask-to-polygon plus seam merging is CPU-bound Python (OpenCV findContours per tile, shapely or JTS unions). A GPU engine that is halo-aware and produces seam-free polygons in one pass solves the specific failure the LSP-DETR authors name (split or duplicated cells).
- Foundation-model pipelines (UNI, Virchow, Prov-GigaPath) mostly output embeddings or tile-level predictions. Polygonization matters most for dense segmentation heads (CellViT, which builds on foundation-model encoders, and HoVer-Net) and for tissue and region masks exported to GeoJSON for QuPath or HistomicsTK.
- Moving GPU output straight to polygons, without a host round-trip of a full-resolution label mask, saves PCIe bandwidth and memory. That matters when a single WSI label mask at 40x is about 100k x 100k pixels.

### Gaps
- I found no paper quantifying the share of WSI pipeline wall time spent in contour extraction and merging specifically.
- I did not verify a GeoJSON file-size figure for a full WSI nucleus export (anecdotally hundreds of MB to GB; unsourced).
- HistomicsTK and Prov-GigaPath, UNI and Virchow polygon export paths were not researched in depth.

---

## Q3. Radiotherapy: RTSTRUCT contours from auto-segmentation, mask-to-contour accuracy, adaptive RT time limits

### Takeaway
DICOM RTSTRUCT stores structures as per-slice planar polygons. Every AI auto-contouring system (TotalSegmentator, MONAI, commercial tools) outputs a mask that has to be converted to contours, and every consumer converts it back. Each round trip causes software-dependent volume and feature discrepancies, worst for small organs (10–20%). Online adaptive RT (Ethos, MR-Linac) has a 15–20+ minute patient-on-couch budget that includes auto-segmentation. The conversion problem is mainly about fidelity and consistency, not speed.

### Cited Findings
- When RTSTRUCT contours are rasterized to dose-grid masks, volume discrepancies are about +1% for large structures (PTV), about +2.5% for medium (brainstem), and 10–20% for small organs (lenses, optic structures) — [BioRT issue #29](https://github.com/mesava/BioRT/issues/29)
- Whybra & Spezi (Sci Rep 2023). DICOM stores contours as coordinate points that must be converted to a binary mask. Importing the same dataset into 4 medical imaging platforms produced different masks and radiomic features, which affected patient clustering and multi-centre models — [Sci Rep 2023](https://www.nature.com/articles/s41598-023-41475-w); [PMC10475062](https://pmc.ncbi.nlm.nih.gov/articles/PMC10475062)
- PyRaDiSe (Comput Methods Programs Biomed 2023) describes the auto-segmentation pipeline as segmentation, then mask-to-contour conversion, then RTSTRUCT export. It notes that open-source RTSTRUCT converters "often rely exclusively on 2D reconstruction approaches leading to pixelated contours", and it offers slice-based and volume-based conversion — [PyRaDiSe, ScienceDirect](https://www.sciencedirect.com/science/article/pii/S016926072300041X)
- The 3D Slicer forum has an "RTSTRUCT segmentation discrepancy" thread, and Neurostars has "Convert nifti segmentation mask to DICOM RTstruct", both recurring practitioner pain points — [Slicer discourse](https://discourse.slicer.org/t/rtstruct-segmentation-discrepancy/10662); [Neurostars](https://neurostars.org/t/convert-nifti-segmentation-mask-to-dicom-rtstruct/7060)
- Clinically deployed auto-contouring research continues in 2026, with thoracic and cervical models described as "clinically deployed" — [arXiv 2609.16036](https://arxiv.org/pdf/2609.16036); [arXiv 2607.11949](https://arxiv.org/pdf/2607.11949)
- Ethos (Varian, launched 2020) is CBCT-based adaptive RT with AI auto-contouring of OARs. Treatments take about 15–20 minutes, and MR-Linac adaptive sessions are generally longer. Adaptive time covers auto-segmentation, review and editing, re-optimization, dose calculation and QA — [Varian Herlev case](https://d8-prod.varian.com/tr/node/11760); [PMC8504605](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8504605/)
- AI auto-segmentation at a 1.5 T MR-Linac is feasible and reduces delineation time — [PMC10289262](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10289262/); [SDU](https://portal.findresearcher.sdu.dk/en/publications/feasibility-and-time-gain-of-implementing-artificial-intelligence/)

### Inferences
- In RT the bottleneck is human review and editing time plus fidelity, not raw speed. A CT has about 100–300 slices at 512², which a CPU contours in under a second per structure. The value of a GPU engine here is (a) a deterministic, topology-correct ring output with explicit hole and inner-ring handling, (b) exact area preservation, where the polygon area equals the mask area under a defined convention, and (c) a coordinate-transform hook for patient (mm) coordinates. RTSTRUCT has no native hole semantics; holes are encoded with keyhole contours or vendor-specific inner contours. A converter that emits explicit outer and inner rings, with a documented area convention, targets the 10–20% small-organ discrepancy directly.
- Speed may matter for real-time MR-Linac with intrafraction tracking (cine-MR at several Hz) and for re-contouring 100+ structures (TotalSegmentator-scale) per fraction. This is inference; I found no source quantifying conversion latency inside an adaptive workflow.

### Gaps
- I found no source on commercial auto-contouring (Limbus, MIM Contour ProtégéAI, Ethos) mask-to-contour algorithm details or timings.
- I found no AAPM Task Group report specifically on mask-to-RTSTRUCT conversion standards in 2023–2026. TG-263 covers nomenclature, not geometry.
- I did not verify TotalSegmentator's RTSTRUCT export path (e.g., rt-utils) or its known issues.

---

## Q4. Connectomics and EM: meshes, skeletons and compression from petascale segmentation

### Takeaway
Connectomics works in 3D. Segmentations are turned into multi-resolution meshes and skeletons with Igneous (distributed, cloud-scale), not into 2D polygons. H01 alone is 1.4 PB. A 2D scanline polygonizer is relevant only as a slice-wise primitive, for 2D proofreading overlays, or for per-section annotation export. The natural connectomics fit is a 3D marching-cubes or surface-nets analogue, which is outside the engine's current scope.

### Cited Findings
- Igneous (Silversmith et al., Frontiers in Neural Circuits 2022) is a Python distributed framework for "economical meshing, skeletonization, image hierarchy creation" on cloud or cluster, built because EM segmentations "are now petascale and growing" — [Frontiers 2022](https://www.frontiersin.org/journals/neural-circuits/articles/10.3389/fncir.2022.977700/pdf); [Princeton record](https://collaborate.princeton.edu/en/publications/igneous-distributed-dense-3d-segmentation-meshing-neuron-skeleton/)
- H01 is a 1.4 PB human cortex dataset (about 1 mm³) with tens of thousands of reconstructed neurons, millions of fragments and 130M annotated synapses — [Google Research blog](https://research.google/blog/a-browsable-petascale-reconstruction-of-the-human-cortex/?m=1)
- MICrONS meshes are static multi-LOD or dynamic, backed by proofreading, and are downloaded through CloudVolume — [MICrONS tutorial](https://tutorial.microns-explorer.org/quickstart_notebooks/06-cloudvolume-download-mesh.html); [MICrONS data](https://microns-explorer.org/data)

### Inferences
- The engine's ideas carry over to connectomics: a halo-aware, single-pass, local-case-table extraction (2x2 becomes 2x2x2 marching cubes) with stitching-free chunk boundaries. Re-meshing after proofreading edits (FlyWire/CAVE "dynamic meshes") is latency-sensitive and touches small local chunks, a plausible GPU niche, but this is speculative.

### Gaps
- No data on zmesh or Igneous per-chunk meshing throughput, cost per teravoxel, or segmentation compression (compresso, crackle) benchmarks. Not searched within the call budget.
- FlyWire and CAVE dynamic mesh latency figures were not found.

---

## Q5. Live-cell, high-content screening, organoids: real-time shape tracking

### Takeaway
No 2023–2026 sources were found in this pass. This branch is unresearched.

### Cited Findings
- None. (The search budget went to spatial omics, pathology and RT.)

### Inferences
- Plausible fits: smart microscopy and event-driven acquisition, where shape changes trigger imaging; per-frame contour output for tracking and morphodynamics; HCS plates with millions of cells per plate. Not verified.

### Gaps
- No sources on HCS polygon export (CellProfiler, Harmony), organoid contour tracking, or real-time feedback microscopy latency budgets.

---

## Q6. Existing GPU tools for contour extraction and polygonization

### Takeaway
cuCIM (RAPIDS) provides a CuPy-based GPU port of the scikit-image API. Research and patents exist on GPU marching squares (including for optical proximity correction). None of the dominant biomedical vectorization paths found here use GPU vectorization: SpatialData, Harpy and Sparrow use rasterio or skimage, QuPath uses JTS or Labels2Rois, and the WSI tile-merge stage runs on shapely STRtree.

### Cited Findings
- cuCIM is a RAPIDS library for accelerated n-dimensional image processing and I/O. `cucim.skimage` is a CuPy-based CUDA implementation of the scikit-image API — [NVIDIA blog](https://developer.nvidia.com/blog/cucim-rapid-n-dimensional-image-processing-and-i-o-on-gpus/); [EPFL field guide](https://imaging.epfl.ch/field-guide/sections/performance_optimization/notebooks/performance_cucim.html)
- skimage `find_contours` uses marching squares with linear interpolation (iso-contours, not exact pixel-edge rings) — [skimage source](https://github.com/scikit-image/scikit-image/blob/v0.14.1/skimage/measure/_find_contours.py)
- A USPTO patent, "Determining contour edges for an image", describes GPU embarrassingly parallel marching squares for contour edge extraction and stitching, applied to optical proximity correction — [USPTO 12315131](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12315131)
- A University of Bologna thesis on GPU contour work also exists — [amslaurea thesis](https://amslaurea.unibo.it/id/eprint/29780/1/Thesis.pdf)

### Inferences
- Positioning: skimage and cuCIM-style marching squares give unordered iso-line segments, or open contours at borders, with no instance or topology semantics. Biomedical consumers need closed, oriented, per-instance outer and inner rings in physical coordinates (GeoJSON, GeoParquet, RTSTRUCT). That is the engine's output, and it is a real gap if cuCIM lacks it.
- For label images, a per-label or "occupancy = label == k" pass is needed. Handling a multi-label image in one pass (2x2 classification over label equality) would be a major differentiator for cell masks with touching cells, where binary masks merge neighbours.

### Gaps
- I did not confirm whether `cucim.skimage.measure.find_contours` exists in current cuCIM releases. Treat it as unverified.
- No benchmark was found comparing GPU and CPU polygonization on biomedical label images.
- Whether cuSpatial (RAPIDS) point-in-polygon is used in any spatial-omics pipeline was not established.

---

## Q7. Practitioner complaints: slowness, topology errors (holes, touching cells), file size

### Takeaway
The complaints are consistent. Converting masks to polygons takes hours at 600k+ cells per image. JTS label reading broke on large tiles in QuPath. Voxel-based polygons are jagged or have empty spaces (Proseg #165). Vendor viewers mis-render third-party polygons (Proseg README on Xenium Explorer). Tile seams split or duplicate cells (LSP-DETR). RT conversion is software-dependent, with 10–20% volume errors for small structures.

### Cited Findings
- "Several hours" to convert numpy masks to polygons at more than 600k cells per image, across 100+ images — [image.sc 121900](https://forum.image.sc/t/fastest-way-to-import-qupath-detections-from-numpy-convert-numpy-mask-into-shapely-polygons/121900)
- JTS label reading caused large-tile and large-FOV failures in the BIOP Cellpose QuPath extension — [BIOP releases](https://github.com/biop/qupath-extension-cellpose/releases)
- Proseg "Ultra jagged masks with empty spaces" — [proseg #165](https://github.com/dcjones/proseg/issues/165); Xenium software "mangles" Proseg polygons — [proseg README](https://github.com/dcjones/proseg)
- Heuristic tile stitching "frequently split[s] or duplicate[s] cells" — [LSP-DETR](https://arxiv.org/html/2601.03163v2)
- RT volume discrepancies of 1% (large), 2.5% (medium) and 10–20% (small) — [BioRT #29](https://github.com/mesava/BioRT/issues/29); cross-software mask discrepancies — [Whybra & Spezi 2023](https://www.nature.com/articles/s41598-023-41475-w)
- Rasterio is recommended over skimage for "more precise vectorization", which implies skimage contours are lossy for label masks — [harpy vectorize](https://harpy.readthedocs.io/en/latest/generated/harpy.sh.vectorize.html)

### Inferences (candidate "intellectually strong" use cases, ranked)
1. **Seam-free, GPU-resident polygonization for spatial omics and SpatialData.** Implement `vectorize(labels) -> GeoDataFrame` on GPU, with halo-based tile ownership that removes Sopa-style conflict resolution, and plug in the pixel-to-µm affine transform. Fuse it with transcript-to-cell assignment, either raster lookup or GPU PIP. Evidence: hours-long conversions, millions of cells per slide, an official push toward GeoJSON and GeoParquet (Xenium, Visium HD, SpatialData). This is the strongest fit.
2. **Polygonizing transcript-assignment segmenters (Proseg, Segger, Baysor).** Their native output is voxel or assignment occupancy, and they have open boundary-quality complaints (jagged rings, empty spaces, unspecified boundary algorithm in Segger). This is a clean binary-occupancy-to-ring problem, and it applies to 3D per-z-layer polygons too.
3. **WSI nucleus and tissue export at gigapixel scale.** Polygonize on GPU directly from model output tensors, avoiding the host-side label mask, OpenCV and shapely merge. This targets the "split or duplicate cells" seam failure.
4. **RTSTRUCT fidelity converter.** Use a deterministic, area-preserving, hole-aware ring extraction with a patient-coordinate transform, plus a published accuracy benchmark against the 1% / 2.5% / 10–20% discrepancy figures. The value is correctness and reproducibility, not speed, and there is regulatory friction.
5. Connectomics and live-cell work need 3D extension or remain unverified. They are lower priority for a 2D engine.

### Gaps
- No quantitative GeoJSON or GeoParquet file-size complaints were sourced.
- No sourced complaints specifically about holes or touching-cell topology in mask-to-polygon conversion (beyond Proseg's "empty spaces"). These are likely to exist on image.sc and GitHub but were not found within the budget.
- No Biostars threads were found.
