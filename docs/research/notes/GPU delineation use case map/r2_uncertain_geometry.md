# Uncertain Geometry: ensembles, probability fields and candidate masks to polygons with quantified boundary uncertainty (state of the art 2021–2026)

Research note (round 2, deep dive). Date: 2026-10-07.

Method note: the shared web-search budget ran out partway through. Remaining facts come from primary PDFs downloaded directly (arXiv/OSTI), operator documentation pages, and a small local CPU benchmark. That benchmark is my own measurement on an Apple M3 laptop using Python, numpy, contourpy, scikit-image and rasterio, run on synthetic smooth fields. Treat its numbers as order-of-magnitude only. Several operational questions stayed open and are listed under Gaps.

---

## 1. Uncertainty visualization of isocontours (contour boxplots, band depth, probabilistic marching squares/cubes): state of the art, cost, GPU, raster or polygons?

### Takeaway
The field is active and moving fast in 2024–2026. ORNL's Athawale group (Viskores, formerly VTK-m) has closed-form and GPU-parallel probabilistic marching cubes. Contour depth has gone from O(N³) to O(N²) to linear-time GPU variants (Inclusion Depth 2024, Probabilistic Inclusion Depth Dec 2025), and a distribution-free Hoeffding confidence band for isocontours came out in Aug–Sep 2026.

Almost all of these methods output **per-cell rasters**: level-crossing probability, entropy, or "the union of cells forms the band". Some draw contour boxplots as union/intersection masks rendered on screen. None of the papers I read emits **topologically valid, exportable polygons** (with holes) for the confidence bands. They are visualization research for ParaView/VisIt/Matplotlib, not GIS deliverables. GPU work covers the probability and depth computation, not band-to-polygon extraction.

### Cited Findings
**Classic foundations**
- Contour boxplots and Contour Band Depth (CBD) were introduced by Whitaker, Mirzargar & Kirby, IEEE TVCG 19(12):2713–2722, 2013. Curve boxplots followed from Mirzargar, Whitaker & Kirby, TVCG 20(12):2654–2663, 2014 (doi 10.1109/TVCG.2014.2346455). — [cited in Chaves-de-Plaza et al. Inclusion Depth](https://graphics.tudelft.nl/~klaus/papers/inclusion_depth.pdf); [cited in PID paper](https://arxiv.org/pdf/2512.15187)
- Positional uncertainty of isocontours appeared in Pöthkow & Hege, TVCG 17(10), 2011 (doi 10.1109/TVCG.2010.247). Probabilistic marching cubes appeared in Pöthkow, Weber & Hege, CGF 30(3), 2011 (doi 10.1111/j.1467-8659.2011.01942.x). — [references in PID paper, arXiv 2512.15187](https://arxiv.org/pdf/2512.15187)

**Contour depth (ensemble → median, 50% band, outliers)**
- Inclusion Depth (ID): Chaves-de-Plaza, Mody, Staring, van Egmond, Vilanova & Hildebrandt, TU Delft, IEEE TVCG 2024.
  - ID is O(N²) in the number of members. CBD is "at least O(N³)".
  - Inputs are **binary masks** (300×300 numpy arrays in the experiments), obtained "as the output of a segmentation algorithm or by thresholding scalar fields using an iso-value".
  - Timings on an M1 Pro CPU without GPU, N=100: CBD 619 s, mCBD 703 s, ID 2.3 s, eID 7.4 s. CBD was only run up to N=100 because "after this point, the CBD method took too long".
  - Radiotherapy example: on a 120-member ensemble of parotid/brain-stem segmentations (540×540 slices, thresholded at 0.8), mCBD took more than 20 minutes and eID took seconds.
  - Weather example: on 500 hPa geopotential-height contours, mCBD took 5.6 s and eID 0.2 s.
  - The authors note that the outer loop is "highly parallelizable, so it could significantly profit from GPU acceleration". They did not implement this.
  - [Inclusion Depth for Contour Ensembles (TVCG 2024)](https://graphics.tudelft.nl/~klaus/papers/inclusion_depth.pdf)
- Depth for multi-modal contour ensembles: Chaves-de-Plaza et al., CGF 2024, which also adds a linear-time ε-ID. — [Wiley CGF 15083](https://onlinelibrary.wiley.com/doi/10.1111/cgf.15083)
- Probabilistic Inclusion Depth (PID): Wu, Lin, Yu, Guo, Klötzl, Weiskopf & Zhou, arXiv 2512.15187, Dec 2025, IEEE TVCG format.
  - It supports ensembles of **soft/fuzzy masks**, "such as soft masks from modern segmentation methods", and treats contour extraction as a fuzzy decision over an isovalue distribution.
  - PID-mean is linear-time. A CUDA implementation ran on an RTX 4080 Super.
  - Ensemble-size scaling at 50³ voxels and 10k–20k members: eID 48.6–248.9 s, PID-mean CPU 4.2–11.3 s, PID-mean GPU 1.23–2.00 s (39.6–124.4× over eID).
  - Resolution scaling at 50 members and 256³–352³ voxels: eID 21.8–216 s, GPU 2.4–5.5 s.
  - eID ran out of memory beyond 20,000 members or 352³ voxels.
  - 3D boxplots are built "by taking the unions and intersections of surfaces".
  - Medical demo: 31 SegResNet models on MSD brain-tumour MRI at 224×224×144.
  - [arXiv 2512.15187](https://arxiv.org/pdf/2512.15187)
- Accelerated surface-boxplot depth with deep learning: Han, Athawale, Li & Johnson, IEEE VIS 2025. Surface boxplot cost "increases cubically with the number of ensemble members, making it impractical for integration into visualization tools". The DL surrogate gives 6× (CPU) and 15× (GPU) speedups on a 50-member 2D Red Sea ensemble. — [OSTI 2538111](https://www.osti.gov/servlets/purl/2538111)
- Coherent visualization of 2D contour ensembles with a VAE latent space (MLS depth and density plots): Wu, Lin, Yu & Zhou, arXiv 2607.24596, Jul 2026. Demonstrated on a weather ensemble. It reports depth–density coherence of R²=0.677, compared with 0.370 for EnConVis+eID. — [arXiv 2607.24596](https://arxiv.org/pdf/2607.24596)
- Alternative raster-band approaches already exist:
  - "Contour probability plots, which extract bands by thresholding a scalar field of percentages" (Kumpf et al.).
  - EnConVis, which uses pointwise kernel density to extract bands.
  - [Inclusion Depth paper, related work](https://graphics.tudelft.nl/~klaus/papers/inclusion_depth.pdf)

**Probabilistic marching squares/cubes (ORNL / Utah SCI)**
- Uncertainty visualization of the 16 marching-squares and 256 marching-cubes topology cases (Athawale, Sane & Johnson, IEEE VIS 2021) uses topology-case counts and entropy per cell. Results are shown by color mapping and entropy isosurfaces, i.e. per-cell scalar output. — [arXiv 2108.03066](https://arxiv.org/pdf/2108.03066)
- Accelerated probabilistic marching cubes by deep learning (Han, Athawale, Pugmire & Johnson, IEEE VIS 2022). It predicts the **level-crossing probability** field. Parallel PMC is up to 19× faster than serial, and the DL surrogate adds up to 10× more. The test case was the Red Sea ensemble on a 500×500×50 domain. — [arXiv 2207.07260](https://arxiv.org/pdf/2207.07260)
- Parallel filters on Viskores/VTK-m:
  - FunMC² (Wang, Athawale, Moreland, Chen, Johnson & Pugmire, EGPGV 2023).
  - FunM²C for multivariate data (Hari et al.).
  - Data-Driven Computation of PMC (Athawale, Wang, Johnson & Pugmire, EuroVis short paper 2024).
  - [references in MAGIC](https://www.osti.gov/pages/servlets/purl/3384930); [FunMC² on ResearchGate](https://www.researchgate.net/publication/371901804_FunMC2_A_Filter_for_Uncertainty_Visualization_of_Marching_Cubes_on_Multi-Core_Devices)
- MAGIC (Athawale, Moreland, Pugmire, Johnson, Rosen, Norman, Georgiadou & Entezari; IEEE TVCG, author version 2025/26).
  - It gives closed-form marching-cubes uncertainty for spatially correlated Gaussian data and avoids Monte Carlo.
  - It is integrated into Viskores, so it runs in ParaView and VisIt.
  - Speedup is "up to 15 × 39 = 585×" over single-core Monte Carlo, measured on a Frontier MI250X GPU and a 112-core EPYC.
  - It also extracts the "most probable topology" and the "expected position" of edge crossings.
  - The 2D wind-ensemble results were done in Python/Matplotlib.
  - [OSTI 3384930](https://www.osti.gov/pages/servlets/purl/3384930)
- Distribution-agnostic isocontour confidence bounds (Ouermi, Gottschling, Gorczowski & Athawale; ORNL/UIUC; arXiv 2608.18220, v2 14 Sep 2026).
  - Vertex-wise Hoeffding confidence intervals are propagated into an isocontour **confidence band**.
  - A cell is included if the isovalue falls within any vertex's interval, and "the union of these cells forms the confidence band". The output is a cell set, not polygons.
  - Runtime on a 15-member wind ensemble: 20 ms, against 21 ms for the Gaussian band and 78 ms for bootstrap.
  - It is explicitly motivated by Gaussian and bootstrap bands that "may underestimate uncertainty for ensembles with a limited number of samples".
  - [arXiv 2608.18220](https://arxiv.org/pdf/2608.18220)
- Other ORNL 2024–2025 work: critical-point uncertainty for parametric and nonparametric models (TVCG 2025), isosurface uncertainty from linear vs high-order interpolation (2024), and local divergence of 2D vector fields "parallelized with Viskores" (2025). — [arXiv 2407.18015](https://arxiv.org/pdf/2407.18015); [arXiv 2409.00043](https://arxiv.org/pdf/2409.00043); [arXiv 2510.01190](https://arxiv.org/html/2510.01190v1)

### Inferences
- **Output gap confirmed:**
  - Probabilistic marching cubes and squares produce per-cell probabilities or entropies, to be colormapped or volume-rendered.
  - Hoeffding, Gaussian and bootstrap bands produce cell unions.
  - Contour boxplots compute depth on binary masks, then draw union/intersection bands.
  - None of these ships a georeferenced polygon with holes for "the 50% band" or "the 90% inner/outer envelope". A bit-quad polygonizer that turns the band masks into rings (intersection of the deepest k members as the inner ring, union as the outer ring) is directly complementary.
- **The compute bottleneck is depth, not contouring:**
  - At N ≥ 1000 the classical depths (CBD O(N³), surface boxplot cubic) are what's slow.
  - ID/PID already fix this, and a GPU implementation exists (PID).
  - So "GPU for depth" is no longer open territory. "GPU for many-ring extraction plus band polygons plus export" still is.
- The ORNL work lives in the HPC visualization stack (Viskores, ParaView, VisIt). There is a clear gap between this and the GIS/Earth-observation stack (GeoJSON, GeoParquet, PostGIS), where none of these methods has been ported, as far as I found.

### Gaps
- I couldn't extract wall-clock numbers from FunMC² or the EuroVis 2024 data-driven PMC paper (no PDF access within budget).
- I found no paper that produces vector polygons (with holes) for contour-boxplot bands or PMC credible regions. Absence is from targeted reading, not an exhaustive search, because the search budget ran out.

---

## 2. ML segmentation uncertainty to geometry: conformal segmentation, SAM ambiguity, ensembles. Does anything produce polygon-level confidence ("conformal polygons")?

### Takeaway
The statistical machinery for **guaranteed inner/outer sets** now exists, and it is exactly the "confidence-banded polygon" concept in raster form.

- Conformal inner/outer confidence sets: Davenport 2024.
- Conformal morphological margins: CONSEMA, MICCAI 2025.
- Conformal false-positive control via nested shrunken masks: Nov 2025.
- Conformal **instance** sets built from a **threshold/parameter sweep**, applied to agricultural field delineation (MIT, Feb–Jul 2026).

All of them output pixel masks. None outputs polygons. Several generate their nested candidates by **sweeping a score threshold**, which is precisely the multi-threshold polygonization workload. So "conformal polygons" is open as a deliverable but already backed by theory.

### Cited Findings
- Davenport, "Conformal confidence sets for biomedical image segmentation" (arXiv 2410.03406, Oct 2024).
  - Inner set I(X) and outer set O(X) come with a finite-sample guarantee that they contain the true mask with the desired probability.
  - Thresholds are calibrated on the max of transformed logit scores inside and outside the ground truth.
  - Distance-transformed scores (with boundaries obtained "using the marching squares algorithm") give tight outer sets.
  - Validated on polyp segmentation.
  - [arXiv 2410.03406](https://arxiv.org/html/2410.03406v2)
- Mossina & Friedrich (IRT Saint Exupéry), "Conformal Prediction for Image Segmentation Using Morphological Prediction Sets" (CONSEMA, MICCAI 2025). It calibrates a **dilation margin** around a predicted mask so that mask plus margin contains the ground truth with user-set confidence. It needs only the predicted masks, and "the size of the margin serves as an indicator of predictive uncertainty". — [arXiv 2503.05618](https://arxiv.org/abs/2503.05618); [code](https://github.com/deel-ai-papers/consema)
- Mossina & Friedrich, "Controlling False Positives in Image Segmentation via Conformal Prediction" (arXiv 2511.15406, Nov 2025). It defines "a nested family of shrunken masks obtained either by increasing the score threshold or by applying morphological erosion". The calibrated inner mask keeps the false-positive proportion below tolerance with probability ≥1−α. It is model-agnostic and needs no retraining. — [arXiv 2511.15406](https://arxiv.org/abs/2511.15406); [code](https://github.com/deel-ai-papers/conseco)
- Lu, Kluger, Bates & Wang (MIT), "Conformal Prediction Sets for Instance Segmentation" (arXiv 2602.10045, v2 Jul 2026).
  - Given a pixel query, it returns a **set of candidate instance masks** with a guarantee that at least one has high IoU with the truth, e.g. "80% coverage at IoU > 0.7".
  - Candidates are generated by **sweeping a tunable parameter T**: "thresholds before connected-components analysis, and thresholds in watershed segmentation". In field delineation they "vary the watershed algorithm threshold T to generate multiple segmentations for each image".
  - Applied to agricultural field delineation, cell segmentation and vehicle detection. The authors describe it as "the first to construct confidence sets of diverse segmentation predictions".
  - [arXiv 2602.10045](https://arxiv.org/pdf/2602.10045)
- COMPASS (Sep 2025) does conformal prediction for medical segmentation **metrics** rather than geometry. — [arXiv 2509.22240](https://arxiv.org/abs/2509.22240)
- SAM is ambiguity-aware by design. It predicts **3 masks per prompt** ("nested masks are often at most three deep: whole, part, and subpart") and ranks them with a predicted-IoU confidence score. — [SAM, arXiv 2304.02643](https://arxiv.org/abs/2304.02643)
- Deep ensembles as sample sets:
  - PID used 31 SegResNet models for brain-tumour soft masks.
  - Inclusion Depth used a 120-member bootstrap/initialization ensemble for head-and-neck organs-at-risk (MICCAI H&N challenge data).
  - Both summarize these as contour boxplots.
  - [arXiv 2512.15187](https://arxiv.org/pdf/2512.15187); [Inclusion Depth TVCG 2024](https://graphics.tudelft.nl/~klaus/papers/inclusion_depth.pdf)
- Round-1 context:
  - Carbon Mapper builds 72 candidate wedge masks per plume, thresholds the resulting probability at 0.3/0.7, and ships a **GeoTIFF**. Masking uncertainty is reported only as the standard deviation of emission estimates. — [Carbon Mapper L3/L4 ATBD, Feb 2026](https://assets.carbonmapper.org/documents/L3_L4%20Algorithm%20Theoretical%20Basis%20Document_021726.pdf)
  - The Jan 2026 "Earth Embeddings as Products" survey says embedding products lack "explainability and uncertainty quantification". — [arXiv 2601.13134](https://arxiv.org/html/2601.13134v1)

### Inferences
- A "conformal polygon" product is well-defined from published pieces. Calibrate a threshold, erosion or dilation parameter (Davenport, Mossina), then polygonize the inner and outer masks into nested rings. The result is a pair of rings per object with a stated coverage level (e.g. 90%).
- Computationally this is cheap per object. The heavy case is **calibration**, which needs masks for every candidate parameter value over the whole calibration set (thousands of images × dozens of thresholds), and the **instance-set** variant (MIT), which needs many candidate segmentations per image.
- Field delineation for agriculture and EUDR is where the MIT instance-conformal method is already applied. That is the most direct bridge to a GIS polygon deliverable.
- No paper I read outputs vectors. All evaluate in pixel space (IoU, false-positive proportion, coverage).

### Gaps
- I didn't verify probabilistic U-Net or diffusion-segmentation papers from 2024–2026 that produce sample sets, or their typical sample counts. Search budget was exhausted.
- I found no source on the SAM 2/SAM 3 multi-mask or score behaviour beyond the original SAM paper.
- I found no published work that names or ships "conformal polygons" or vector inner/outer envelopes.

---

## 3. Operational demand: who needs uncertain polygons? (NWS FACETs/PHI, ECMWF AIFS ENS, flood, EUDR, insurance, radiotherapy)

### Takeaway
Demand is real but uneven.

**Strongest concrete demand, current products shipping uncertainty as rasters next to vector extents:**
- Copernicus GFM: flood extent is offered as GeoJSON, while likelihood ships only as COG.
- Carbon Mapper plumes.
- NWS's long-standing FACETs design, in which warnings are explicitly "derived from the exceedance of predefined probability thresholds" on ≤2 km PHI grids. This is a multi-threshold polygonization spec.
- Radiotherapy, which has the most mature research use of contour ensembles. Clinical product status is unverified.

**Weaker or unverified demand:**
- EUDR, which is a legal driver for polygons but does not require confidence.
- Insurance parametric triggers. I found no evidence.

### Cited Findings
**Weather**
- FACETs envisions "high-resolution forecasts of probabilistic hazard information conveyed on fine-resolution (≤2 kilometers) grids". Legacy watches and warnings would be "derived from the exceedance of predefined probability thresholds". Warn-on-Forecast multi-member ensembles provide the probabilistic storm-scale input. The article dates from 2014, updated Aug 2025; its 2019 "first operational iterations" target is evidently not met in that text. — [Eos FACETs article](https://eos.org/science-updates/next-generation-severe-weather-forecasting-communication)
- Threats-in-Motion (TIM): polygons "update every minute" and move with the storm. "TIM is essentially PHI, but without the probabilities (which will come later as the science evolves)." It is proposed as the first FACETs stage. — [NSSL FACETs TIM tag](https://inside.nssl.noaa.gov/facets/tag/threats-in-motion/); [NOAA repository 28465](https://repository.library.noaa.gov/view/noaa/28465); [GSL PHI news](https://gsl.noaa.gov/news-media/news/probabilistic-hazard-information)
  - The NSSL TIM explainer is dated **16 March 2021** (Kodi Berry) and is the only post under that tag that I retrieved.
  - PHI is shown as continuously updating probabilistic "plumes" in a prototype **AWIPS2 Hazard-Services-PHI (HS-PHI)**.
  - Its inputs (radar, satellite, lightning, MRMS, ProbSevere) update "at 1- to 2-minute intervals".
  - It names a forecaster-workload problem: "the workload of the forecaster is too great to keep up with the issuance of new downstream warnings".
  - [NSSL TIM post](https://inside.nssl.noaa.gov/facets/tag/threats-in-motion/)
- ProbSevere model output is "shapefiles contoured around radar storm cells", each with ProbHail, ProbWind and ProbTor. — [NOAA VLab ProbSevere briefing](https://vlab.noaa.gov/documents/10157/137122/ProbSevere_Vlab_Briefing.pdf/fb318646-32ab-6bb6-2213-1f6291fdd8b8)
- ECMWF's AIFS ENS became **operational on 1 July 2025**: 51 members at ~31 km (versus 9 km physics IFS), "over 10 times faster" than physics-based systems. The announcement mentions no polygon or warning-object products. — [ECMWF news, 2025](https://www.ecmwf.int/en/about/media-centre/news/2025/ecmwfs-ensemble-ai-forecasts-become-operational)
- HENS: 7,424 members per initialization day (SFNO), summer 2023 hindcasts. — [arXiv 2408.01581](https://arxiv.org/abs/2408.01581); [arXiv 2408.03100](https://arxiv.org/pdf/2408.03100)

**Flood**
- The Copernicus/GloFAS **Global Flood Monitoring (GFM)** output layers include Observed Flood Extent, delivered as "COG / geojson" and retrievable "as vector (geojson), raster (COG), or both". They also include **Likelihood Values (COG only)**. This is "a simplified 'appraisal of trust' in the ensemble flood mapping method" that combines **three** flood-mapping algorithms: values from 0–100, 50 separates flooded from unflooded, and values toward 100 mean higher flood confidence. — [EODC GFM Product User Manual, output layers](https://extwiki.eodc.eu/en/GFM/PUM/Products)

**Methane**
- Carbon Mapper uses 72 candidate masks to build a per-pixel probability, applies 0.3/0.7 thresholds, and ships a GeoTIFF. Plume length is taken from the convex hull of the segmented plume. — [Carbon Mapper ATBD](https://assets.carbonmapper.org/documents/L3_L4%20Algorithm%20Theoretical%20Basis%20Document_021726.pdf)

**EUDR**
- Regulation (EU) 2023/1115, Art. 2(28): for plots of land over 4 ha (commodities other than cattle), geolocation "shall be provided using polygons with sufficient latitude and longitude points to describe the perimeter of each plot of land". The Commission must review "the threshold for mandatory use of polygons". The text sets no positional-accuracy or confidence requirement. — [EUR-Lex 32023R1115](https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:32023R1115)

**Radiotherapy / medical**
- Contour ensembles of head-and-neck organs-at-risk (parotid gland, brain stem) from a 120-member model ensemble are summarized as contour boxplots for treatment-planning analysis. "In image-guided medical specialties, ensembles of segmentations are analyzed for planning the patients' treatments." — [Inclusion Depth TVCG 2024](https://graphics.tudelft.nl/~klaus/papers/inclusion_depth.pdf)
- Conformal segmentation papers (polyps, MICCAI 2025) are motivated by clinical risk, e.g. "false positives may trigger unnecessary interventions". — [arXiv 2511.15406](https://arxiv.org/abs/2511.15406); [arXiv 2503.05618](https://arxiv.org/abs/2503.05618)

### Inferences
- **Strongest real demand, ranked:**
  1. **Earth-observation hazard and emission products that already compute an ensemble or likelihood but ship only a raster**: GFM likelihood (three-algorithm ensemble), Carbon Mapper (72 masks), and embedding-similarity maps. These organizations already ship vector extents, so adding nested "likelihood ≥ 50/70/90" rings is incremental, cheap and explainable.
  2. **NWS PHI/TIM**, whose design is literally multi-threshold polygonization of probability grids every minute. Operational status is unverified (see Gaps).
  3. **Radiotherapy contour QA and treatment planning**, which has the most mature research use but is 3D and mostly served by medical vendors.
  4. **Ag field delineation, i.e. EUDR due diligence**, where MIT's instance-conformal sets apply directly. The regulation does not demand confidence, but disputes over plot boundaries plausibly would. That last point is unverified.
- For ECMWF AIFS ENS (51 members at 31 km), the raw polygonization load is trivial (see §5). The value there lies in representation (banded hazard objects), not speed.

### Gaps
- I couldn't verify the **2025–2026 operational status** of FACETs/PHI/TIM in NWS AWIPS/Hazard Services. The NSSL site failed TLS and the search budget was exhausted. Treat it as still experimental in Hazardous Weather Testbed form unless confirmed.
- No sources were found within budget on wildfire perimeter uncertainty (NIFC/FIRIS perimeters with confidence), parametric insurance triggers based on uncertain polygons, or EUDR enforcement disputes over boundary accuracy.
- I didn't check whether GFM's Product Definition Document defines how the three algorithms' likelihoods are combined. It is referenced but was not read.
- I didn't verify clinical commercial products for probabilistic CTV/PTV (e.g. in TPS vendors).

---

## 4. Data formats and standards for uncertain polygons

### Takeaway
I found no deployed standard for publishing a polygon together with its boundary uncertainty. Real products use one of two de-facto patterns: a vector extent plus a **separate likelihood raster** (GFM), or a raster only (Carbon Mapper). Classic GIScience models (egg-yolk, broad boundaries, fuzzy regions) exist in the literature but could not be verified as implemented in any format within this research.

### Cited Findings
- GFM publishes Observed Flood Extent as GeoJSON or COG and **Likelihood only as COG**, so uncertainty is decoupled from the vector. — [EODC GFM PUM](https://extwiki.eodc.eu/en/GFM/PUM/Products)
- EUDR geolocation is specified as plain polygons (lat/lon points describing the perimeter) with no uncertainty attribute. — [EUR-Lex 32023R1115](https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:32023R1115)
- Carbon Mapper's segmented plume product is a GeoTIFF that "provides the boundaries and concentrations". — [Carbon Mapper ATBD](https://assets.carbonmapper.org/documents/L3_L4%20Algorithm%20Theoretical%20Basis%20Document_021726.pdf)
- The Earth Embeddings survey (Jan 2026) calls the embedding-product landscape "fragmented" and lacking uncertainty quantification and standards. — [arXiv 2601.13134](https://arxiv.org/html/2601.13134v1)

### Inferences
- A pragmatic encoding that existing tools can read today: per object, a **MultiPolygon per confidence level** (nested rings, e.g. p≥0.9 ⊂ p≥0.5 ⊂ p≥0.1), or inner/outer conformal rings with attributes `level`, `method` (exceedance / conformal / depth-band) and `alpha`. It fits GeoJSON or GeoParquet without a new standard. Nesting (inner ⊆ outer) must be guaranteed. Polygonizing nested thresholds from one field with a consistent pixel-edge tracer gives this for free, because the threshold sets are monotone. That is a point in favour of an exact bit-quad tracer over smoothed marching squares, where independent simplification can break nesting.
- The egg-yolk / broad-boundary model (an inner "yolk" certainly inside, plus an outer "egg" possibly inside) maps one-to-one onto conformal inner/outer sets.

### Gaps
- I could not verify the OGC UncertML status page (404). I didn't check whether any OGC/ISO 19157 data-quality element covers per-feature positional uncertainty bands.
- I couldn't fetch primary sources for the egg-yolk model (Cohn & Gotts) or broad boundaries (Clementini & Di Felice) within budget. They are cited here only as known background, not verified.
- I found no GeoJSON or STAC extension for uncertain or fuzzy geometry.

---

## 5. Computational load: e.g. 1000 members × 1440×721 × 10 thresholds. How long do current tools take?

### Takeaway
At ERA5/ML-weather resolution (0.25°, about 1 Mpx), **raw contouring is not the bottleneck**. My CPU measurements give 4–16 ms per field per threshold, so 1000×10 takes about 40–160 s on one core and seconds on a multi-core machine. Ensemble statistics, including O(N²) inclusion depth, take seconds via matrix products.

The cost becomes significant only with large multipliers:
- HENS scale: 7,424 members × lead times × variables × thresholds gives hours per day of initializations.
- km-scale or 10 m grids: about 25–100× more pixels per field.
- Minute-cadence nowcasting.
- Calibration sweeps for conformal methods.

The rigorous classical methods (CBD O(N³), surface boxplots) are what's slow, at minutes to hours for N=100. Recent algorithms (ID, PID-GPU) solved that, not contouring.

### Cited Findings
**Own benchmark** (synthetic smooth fields, 721×1440 float32, 10 thresholds, Apple M3, single core):
- contourpy (Matplotlib's marching-squares engine) line contours: **3.6 ms** per member-threshold, about 56 contours per field.
- contourpy filled polygons with holes: **4.2 ms**.
- scikit-image `find_contours`: **16.0 ms**.
- rasterio/GDAL `features.shapes` polygonizing the binary mask: **10.5 ms**.
- Extrapolated to 1000 members × 10 thresholds on one core: **36 s (contourpy lines), 42 s (filled), 105 s (rasterio), 160 s (skimage)**.
- Exceedance-probability field over 20 members: about 21 ms per threshold with numpy.
- (Own measurement; not a published source.)

**Own benchmark, depth:**
- ε-inclusion depth for N=100 masks of 721×1440, written as one float32 matrix product (|cᵢ ∩ cⱼ| for all pairs): **0.075 s**. That scales to about 7.5 s for N=1000 on CPU, since it is O(N²·pixels).
- The 50% band (intersection and union of the deepest half) takes 14 ms as rasters.
- (Own measurement.)

**Published numbers:**
- CBD, N=100 on 300×300 masks (Python, M1 Pro): about 620 s. ID: 2.3 s. CBD was not attempted beyond N=100. — [Inclusion Depth TVCG 2024](https://graphics.tudelft.nl/~klaus/papers/inclusion_depth.pdf)
- PID-mean on GPU: 1.2–2.0 s for 10k–20k members on 50³ voxels. eID took 48–249 s and ran out of memory beyond 20k. — [arXiv 2512.15187](https://arxiv.org/pdf/2512.15187)
- Hoeffding, Gaussian and bootstrap isocontour bands on a 15-member wind ensemble: 20/21/78 ms. — [arXiv 2608.18220](https://arxiv.org/pdf/2608.18220)
- MAGIC closed-form plus Viskores: up to 585× over single-core Monte Carlo PMC. — [OSTI 3384930](https://www.osti.gov/pages/servlets/purl/3384930)
- Parallel PMC is 19× over serial, and the DL surrogate adds 10× more. — [arXiv 2207.07260](https://arxiv.org/pdf/2207.07260)

**Multiplier examples:**
- HENS: 7,424 members per day. — [arXiv 2408.01581](https://arxiv.org/abs/2408.01581)
- AIFS ENS: 51 members at 31 km. — [ECMWF](https://www.ecmwf.int/en/about/media-centre/news/2025/ecmwfs-ensemble-ai-forecasts-become-operational)
- FACETs PHI: ≤2 km grids, with TIM updating every minute. — [Eos](https://eos.org/science-updates/next-generation-severe-weather-forecasting-communication); [NSSL TIM](https://inside.nssl.noaa.gov/facets/tag/threats-in-motion/)

### Inferences
- Back-of-envelope for HENS, assuming 7,424 members × ~60 lead times (6-hourly to 15 days) × 10 thresholds ≈ 4.5M fields per variable per day:
  - At about 4 ms per field, that's about 5 CPU-hours per variable per day of initializations, before I/O.
  - A GPU engine doing about 1 Mpx fields in well under 1 ms each in batch would make this minutes.
  - This is my estimate. HENS's real post-processing variables and lead-time count were not verified.
- For AIFS ENS (51 members) or a 1000×10 batch at 0.25°, CPU is already fine, so **speed alone is not a compelling pitch** at that scale.
- The pitch becomes compelling for:
  - (a) many-member × many-threshold × many-lead-time sweeps (HENS-class, or conformal calibration sweeps);
  - (b) high-resolution grids, such as CONUS 1 km at about 7000×3500 (25× ERA5) at minute cadence, or 10 m EO tiles;
  - (c) keeping everything on the GPU next to the model that produced the ensemble, so there is no device-to-host round trip and no file I/O. Round 1 found this is where CPU trackers actually lose (Earth2Studio moving trackers into PyTorch).

### Gaps
- My benchmark used synthetic fields. Real precipitation and probability fields with many small, fragmented objects will produce more rings and run slower, especially for rasterio/GDAL.
- I couldn't measure GDAL Polygonize or contourpy on 7000×3500 grids, or with real HENS output.
- No published wall-clock figure was found for any operational centre's ensemble→polygon step.

---

## 6. Verdict on the hypothesis ("a GPU engine that polygonizes many thresholds/members per launch and summarizes them into confidence-banded polygons is a real, under-served, intellectually rich direction")

### Takeaway
**The hypothesis largely holds, with one important correction.**

- **Under-served: yes.** The visualization literature (contour boxplots, PMC, Hoeffding bands) and the ML literature (conformal inner/outer sets, instance-conformal sets) both stop at **rasters or cell sets**. Operational products (GFM, Carbon Mapper) ship uncertainty as a raster beside, or instead of, vector extents. No standard encoding exists.
- **Intellectually rich: yes.** There are guaranteed nested rings, conformal coverage, depth-ranked bands, topological consistency across thresholds, and Hoeffding distribution-free bounds as polygons.
- **The correction:** **raw speed at 0.25° global resolution is not the pain point.** CPU contouring of 1000×10 ERA5-size fields takes about a minute on one core, and depth via matrix product takes seconds. The value proposition should be **"correct, nested, hole-aware, exportable confidence polygons, computed in-GPU next to the model, at batch scale"**, not "faster contours". Speed matters at HENS-class multiplicity, km/10 m resolution, minute cadence, and calibration sweeps.

### Cited Findings
- Raster or cell-set outputs in the SOTA visualization methods:
  - Union of cells as the band. — [arXiv 2608.18220](https://arxiv.org/pdf/2608.18220)
  - Level-crossing probability. — [arXiv 2207.07260](https://arxiv.org/pdf/2207.07260)
  - Entropy and color mapping. — [arXiv 2108.03066](https://arxiv.org/pdf/2108.03066)
  - Unions and intersections for boxplots. — [arXiv 2512.15187](https://arxiv.org/pdf/2512.15187)
- Mask-only outputs in conformal segmentation: [arXiv 2410.03406](https://arxiv.org/html/2410.03406v2); [arXiv 2503.05618](https://arxiv.org/abs/2503.05618); [arXiv 2511.15406](https://arxiv.org/abs/2511.15406); [arXiv 2602.10045](https://arxiv.org/pdf/2602.10045)
- Operational products that ship vector extents with raster-only uncertainty: [GFM PUM](https://extwiki.eodc.eu/en/GFM/PUM/Products); [Carbon Mapper ATBD](https://assets.carbonmapper.org/documents/L3_L4%20Algorithm%20Theoretical%20Basis%20Document_021726.pdf)
- Threshold sweeps are already the mechanism for generating candidates:
  - Instance-conformal field delineation. — [arXiv 2602.10045](https://arxiv.org/pdf/2602.10045)
  - FP-control nested shrunken masks. — [arXiv 2511.15406](https://arxiv.org/abs/2511.15406)
  - FACETs warnings from probability thresholds. — [Eos](https://eos.org/science-updates/next-generation-severe-weather-forecasting-communication)
- GPU acceleration so far targets depth and probability computation, not polygons: [PID CUDA](https://arxiv.org/pdf/2512.15187); [MAGIC/Viskores](https://www.osti.gov/pages/servlets/purl/3384930)

### Inferences
- **Most defensible product shape** is a batched API:
  - Input: (fields[B], thresholds[T]) or (masks[N]).
  - Output: per-object **nested ring sets**, with guaranteed nesting, holes, and stable object IDs across levels.
  - Optional reducers:
    - (i) exceedance-probability bands from N members;
    - (ii) depth-band polygons (intersection and union of the deepest α-fraction, via eID/PID on GPU);
    - (iii) conformal inner/outer rings, given calibrated λ;
    - (iv) per-object geometric spread metrics (area or Hausdorff spread across members, centroid scatter, Feret-length spread for plumes).
- **Strongest real demand to target first:**
  - (1) Earth-observation hazard and emission products with existing ensembles or likelihoods (GFM floods, Carbon Mapper plumes, embedding-similarity maps), because they already publish vectors.
  - (2) PHI/TIM-style severe-weather polygons, if NWS operational status is confirmed.
  - (3) Ag-field delineation with conformal instance sets (EUDR-adjacent).
  - Radiotherapy is intellectually close but 3D and vendor-locked. It is a weaker target for a 2D bitmap engine.
- **Novelty risk:** the ORNL/Viskores group could add polygon export to its bands, and depth-on-GPU (PID) already exists. The distinctive contribution is the **2D GIS-grade vector output** (exact pixel-edge rings with holes, nesting-consistent across thresholds, georeferenced) plus batching. Neither visualization nor conformal researchers are building this.

### Gaps
- Willingness to pay or adopt is unverified. I found no source in which GFM, NWS or Carbon Mapper asks for vector uncertainty bands. This demand is inferred from what they ship, not from stated requirements.
- I didn't benchmark GPU contouring libraries (e.g. Viskores contour filter, cuCIM find_contours) for the batched multi-threshold case.
