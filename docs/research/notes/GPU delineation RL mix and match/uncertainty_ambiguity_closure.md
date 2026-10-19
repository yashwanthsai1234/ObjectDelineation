# Uncertain geometry and label ambiguity in segmentation: closing the open items (as of Oct 2026)

Research note, written 2026-10-12. It closes the items left open by `research_notes/GPU delineation use case map/r2_uncertain_geometry.md`, which ran out of search budget.

The note assesses two proposals:
- **(a) Confidence-banded polygons:** nested inner and outer rings built from ensembles or probability fields.
- **(b) Uncertainty-tolerant RL rewards:** no penalty inside the ambiguous boundary band, and a penalty outside it.

Method: web search plus primary-source fetches. Some primary pages could not be reached:
- The NSSL FACETs pages failed TLS verification.
- The 2024 HWT Taller-TIM summary PDF and the 2022 Tiny-TIM OPG report now return HTTP 410 Gone or redirect to HTML.
- The GSL pages returned 403.
- PMC returned a CAPTCHA.

Where a fact rests only on a search-engine summary of a page I couldn't open, it is marked **(snippet-only)**.

---

## 1. NOAA FACETs / PHI / Threats-in-Motion: operational status in 2025–2026

### Takeaway
**Verdict: still pre-operational as of Oct 2026, as far as I can find.** Hazard Services, the issuing platform, is in operational use. However, I found no NWS Service Change Notice, NWS announcement or NSSL/GSL page confirming that Tiny TIM, Taller TIM or storm-based PHI probabilities are issued operationally. The latest primary documents describe TIM as a "(potential) path toward NWS operations", contingent on NWS governance approval. PHI appears as a 2025 NOAA technical report summarising HWT experiments.

### Cited Findings
- **TIM's design:** warning polygons are "attached to the threat and move forward along with it", and warnings are "automatically cleared from locations where the threat has passed". The paper reports greater average lead times and lower departure times "with little to no impact to average false alarm time". TIM is proposed as the first stage of FACETs. — [NOAA repository 28465, "NWS Severe Weather Warnings as Threats-in-Motion"](https://repository.library.noaa.gov/view/noaa/28465)
- **Two-phase plan:** "Tiny TIM" extends warnings in time and area within today's framework (one hazard keeps one ETN, so there are fewer polygons). "Taller TIM" is the full continuously-updating moving warning. The most recent HWT experiment was "held in the spring of 2024 in person at the National Weather Center". — [2024 Hazard Services Taller-TIM HWT experiment summary (search listing; PDF now 410 Gone)](https://hwt.nssl.noaa.gov/docs/2024-Taller-TIM-HWT-Experiment-Summary-FINAL.pdf) **(snippet-only)**; [NWS OPG Tiny TIM Final Report, Sept 2022 (listing)](https://goes-r.noaa.gov/users/docs/pg-activities/Sept2022OPG_TinyTIMEvaluationReport.pdf) **(snippet-only)**
- **Status as of 2024:** an AMS 2024 abstract is titled "Advancing Warnings: Threats-in-Motion and its (Potential) Path Toward NWS Operations". Search summaries of NSSL pages state: "If TIM is approved through the NWS governance process, the transition to operations is expected to occur in two phases." — [ADS 2024AMS...10429803B](https://ui.adsabs.harvard.edu/abs/2024AMS...10429803B/abstract) **(title only; abstract fetch returned 405)**; [NSSL FACETs tag "tim"](https://inside.nssl.noaa.gov/facets/tag/tim/) **(snippet-only, TLS failure)**
- **PHI in 2025:** Calhoun, Satrio, Steeves, Sandmæl, Campbell & Berry (2025), "Creation and Use of Automated Storm-Based Probabilistic Hazard Information for NWS Forecasting and Warnings", NOAA technical report.
  - It states: "Experiments in the NOAA Hazardous Weather Testbed (HWT) have shown that storm-based PHI can support forecasters in making more effective short-term forecasts and warning decisions."
  - It describes PHI as "storm-based probabilities of individual hazards along with projected storm tracks".
  - The repository page gives no deployment date.
  - [NOAA repository 71712](https://repository.library.noaa.gov/view/noaa/71712)
- **PHI definition:** "continuously-updating probabilistic hazard grids" for "custom user-specific products", for example "longer lead times, at lower confidence, for more vulnerable populations". — [MDL VLab Hazard Services – PHI](https://vlab.noaa.gov/web/mdl/hazard-services-probabilistic-hazard-information-phi) **(snippet-only)**
- **Hazard Services in operations:** GSL describes Hazard Services as NOAA's next-generation alert/warning platform. "During the March 2025 tornado outbreaks, NWS forecasters were able to issue, update, and adjust alerts more efficiently". GSL says it "aims to make these improvements permanent at NWS forecast offices nationwide". — [GSL Hazard Services](https://gsl.noaa.gov/hazard_services/) **(snippet-only; page 403)**

### Inferences
- The operational platform (Hazard Services) exists. The probabilistic layer (PHI grids thresholded into polygons) and moving polygons (TIM) remain in testbed and governance stages. Proposal (a) should present NWS as a **design-aligned future customer**, not a current buyer.
- The FACETs/PHI concept, "warnings derived from exceedance of probability thresholds" on probability grids, is still the clearest public spec for multi-threshold, nested probabilistic polygons that update every 1–2 minutes (see the r2 note).

### Gaps
- I found no 2025–2026 NWS Service Change Notice or news release confirming Tiny TIM or PHI deployment. Absence is not proof, because the NSSL pages were unreachable (TLS failure).
- The quantitative results of the 2024 Taller-TIM HWT (forecaster counts, lead-time gains) could not be retrieved: the PDF returns 410 Gone.
- It is unclear whether Weather Ready Nation / NWS budget and staffing changes in 2025 affected the FACETs timeline. Not researched.

---

## 2. Standards for uncertain/fuzzy geometry (OGC, ISO 19157, broad boundaries / egg-yolk, GeoJSON/GeoParquet)

### Takeaway
**Verdict: no deployed interchange standard encodes a per-feature broad-boundary or banded polygon.**
- ISO 19157-1:2023 describes data quality, including positional accuracy, mostly at the dataset level.
- OGC's UncertML/QualityML work lives in testbed engineering reports.
- GeoJSON explicitly lacks any concept of uncertainty.
- I found no uncertainty field in GeoParquet.
- The theory is mature: egg-yolk (Cohn & Gotts) and broad boundaries (Clementini & Di Felice) map directly onto inner/outer rings.

Practical consequence: proposal (a) must ship as a convention, i.e. multiple geometries plus attributes, not a standard type.

### Cited Findings
- **GeoJSON:** "GeoJSON has no concept of uncertainty; imprecise or uncertain 'geo' URIs thus cannot be mapped to GeoJSON geometries" (in contrast to RFC 5870 geo URIs, which carry a `u=` uncertainty parameter). — [Wikipedia: Geo URI scheme](https://en.wikipedia.org/wiki/Geo_URI_scheme)
- **GeoParquet:** the spec covers geometry encoding, CRS, bbox and covering columns. Searching the spec and its issues surfaced no uncertainty or fuzzy-geometry extension. — [opengeospatial/geoparquet spec](https://github.com/opengeospatial/geoparquet/blob/main/format-specs/geoparquet.md)
- **ISO 19157-1:2023** ("Geographic information — Data quality — Part 1: General requirements") includes positional accuracy as a quality element (§8.3.4). — [ISO 19157-1:2023 sample PDF (iTeh)](https://cdn.standards.iteh.ai/samples/78900/96ad08103cea4bb7bf4cac7990c7dd38/ISO-19157-1-2023.pdf); [SIST EN ISO 19157-1:2023 catalog](https://standards.iteh.ai/catalog/standards/sist/ca52d0c4-ea3c-4673-9c0e-172c282966a9/sist-en-iso-19157-1-2023)
- **Per-feature uncertainty:** the critique "Geospatial data quality (ISO 19157-1): evolve or perish" argues the standard should handle per-feature or attribute uncertainty more explicitly, as distinct from dataset quality. — [ResearchGate 341111956](https://www.researchgate.net/publication/341111956_Geospatial_data_quality_ISO_19157-1_evolve_or_perish) **(snippet-only)**
- **OGC testbed work:** OGC Testbed-12 (Imagery Quality and Accuracy ER) and Testbed-13 (Data Quality Specification ER) discuss encoding quality using ISO 19157, UncertML and QualityML. These are engineering reports, not adopted standards. — [OGC 16-050](https://docs.ogc.org/per/16-050.html); [OGC 17-018](https://docs.ogc.org/per/17-018.pdf)
- **Broad boundaries:** the crisp boundary is replaced by an area. Under positional uncertainty it is "the set of all possible positions" of the unknown boundary. For inherently fuzzy objects it is "their minimum and maximum extent". Clementini & Di Felice model a broad-boundary region as "two regions, one lying inside the other". — [Objects with Broad Boundaries (ResearchGate)](https://www.researchgate.net/publication/235005951_Objects_with_Broad_Boundaries)
- **Egg-yolk:** the "egg" is the maximal extent and the "yolk" is the certain core, based on mereology/RCC. — [Cohn & Gotts, Egg-Yolk representation (ResearchGate)](https://www.researchgate.net/publication/2827200_The_Egg-Yolk'_Representation_Of_Regions_with_Indeterminate_Boundaries); [A Theory of Spatial Regions with Indeterminate Boundaries](https://www.researchgate.net/publication/2622707_A_Theory_of_Spatial_Regions_with_Indeterminate_Boundaries)
- **Vague regions:** fuzzy and vague region data types exist in the database research literature, e.g. Schneider's fuzzy spatial data types and "Definition and Implementation of Vague Objects". There, points in the broad boundary carry membership values. — [Springer: Fuzzy Spatial Data Types](https://link.springer.com/chapter/10.1007/3-540-48482-5_20); [Definition and Implementation of Vague Objects](https://www.researchgate.net/publication/2935813_Definition_and_Implementation_of_Vague_Objects); [Review: boundary of a fuzzy set in GIS (2022)](https://link.springer.com/article/10.1007/s10462-022-10331-0)
- **De-facto practice:** a vector extent plus a separate likelihood raster. GFM's ensemble likelihood is "the mean of the individual classification likelihoods" of three algorithms on a 0–100 scale. It is delivered alongside the flood-extent mask, not as vector bands (the r2 note found the likelihood is COG-only). — [Krullikowski et al., arXiv 2304.12488](https://arxiv.org/abs/2304.12488)

### Inferences
- Egg-yolk and broad boundary correspond **one-to-one** to conformal inner/outer sets (yolk = inner set, egg = outer set). A nested set of k rings generalises this to "fuzzy region by α-cuts". That gives proposal (a) a 30-year-old theoretical pedigree, which is useful for credibility.
- A workable encoding with today's tools: a FeatureCollection, or a GeoParquet table, with one row per (object_id, level). Columns would be `level` (e.g. 0.1/0.5/0.9 or α), `kind` (inner|median|outer|exceedance|conformal) and `method`, and nesting is guaranteed by construction. ISO 19157 quality reports could cite the conformal coverage as a positional-accuracy "conformance" result. That is an inference, not something the standard prescribes.

### Gaps
- I did not read the ISO 19157-1:2023 full text, which is paywalled. I can't confirm whether it has any element for feature-level positional uncertainty regions.
- OGC UncertML's current status (retired, or merged into OGC API work) is unconfirmed. The r2 note hit a 404, and no current OGC SWG page was found.
- I found no STAC extension, OGC Features extension or GeoJSON community extension for fuzzy geometry. This rests on targeted searches, not an exhaustive registry scan.

---

## 3. Domain demand: wildfire perimeters, flood extent, parametric insurance, EUDR

### Takeaway
- **Wildfire:** operational perimeters (NIFC/WFIGS, FIRIS) are published as crisp polygons with no uncertainty attribute that I could find. Research quantifies 20–166 m perimeter uncertainty, and one 2025 paper states that operational perimeter forecasts "lack uncertainty quantification".
- **Flood:** CEMS/GFM deliberately computes ensemble likelihood, but ships it as a raster.
- **Parametric insurance:** the newest products (Liberty Mutual × ICEYE, June 2026) trigger on building-footprint classification, not on burn perimeters. Radius-based wildfire covers ("did fire touch the circle") are directly exposed to perimeter error. I found no contract language that references perimeter confidence.
- **EUDR:** application now starts 30 Dec 2026 for medium and large operators (secondary sources). Validation is strictly geometric, using GeoJSON; secondary sources state that polygons with holes ("doughnut shapes") are rejected. Disputes arise from polygons leaking into 10 m forest pixels, which is a natural use for inner/outer bands, but the regulation does not demand confidence.

### Cited Findings
**Wildfire**
- **CAL FIRE/NIFC/FIRIS public perimeters** are published as a public-view perimeter layer on data.gov. — [data.gov: CA Perimeters CAL FIRE NIFC FIRIS public view](https://catalog.data.gov/dataset/ca-perimeters-nifc-firis-public-view-330dc)
- **FIRIS** is a public-private programme that delivers real-time fire perimeters from aircraft-borne IR. — [same data.gov record](https://catalog.data.gov/dataset/ca-perimeters-nifc-firis-public-view-330dc) **(snippet-only)**
- **Accuracy of NIFC perimeters vs NIROPS:** "NIFC polygons have very low accuracy at the beginning of the wildfire and there are high commission errors in some study areas". NIROPS airborne IR perimeters are "considered a very accurate representation… on the order of meters (approximately 6.3 m/pixel)". — [Chen et al., VIIRS object-based fire tracking (FEDS), PMC9151742](https://pmc.ncbi.nlm.nih.gov/articles/PMC9151742/) **(snippet-only; PMC CAPTCHA blocked the fetch)**
- **Chakravarty, "Spatial Uncertainty Quantification in Wildfire Forecasting for Climate-Resilient Emergency Planning"** (arXiv 2510.09666, Oct 2025):
  - Data: WildfireSpreadTS, 607 western-US fire events.
  - Models: UTAE with MC Dropout (20 passes), BNN, and Deep Ensembles; Deep Ensembles were best calibrated.
  - Results: "consistent uncertainty buffer zones of 28–35 meters using our centroid metric, 47–64 meters using Average Surface Distance, and 148–166 meters using Hausdorff". "Epistemic uncertainty forms coherent 20–60 m buffer zones aligned with fire perimeters".
  - Gap statement: "Despite the operational risks involved, no prior work has investigated uncertainty quantification" in this setting.
  - [arXiv 2510.09666](https://arxiv.org/pdf/2510.09666)

**Flood**
- **CEMS ensemble approach:** CEMS's request to integrate flood algorithms "resulted in the adoption of an ensemble approach" that merges three independently developed algorithms. The likelihood is intended so users can treat flagged areas "with adequate caution". — [Krullikowski et al., arXiv 2304.12488](https://arxiv.org/abs/2304.12488)
- **Communicating flood-map uncertainty:** research on "Uncertainty in three dimensions: the challenges of communicating probabilistic flood forecast maps" (HESS 2023) uses end-user surveys. It finds probabilistic/ensemble map information "could be overwhelming" for decision-makers. — [HESS 27, 3351 (2023) discussion](https://hess.copernicus.org/articles/27/3351/2023/hess-27-3351-2023-discussion.html)
- **Ensemble flood maps:** spread–skill assessment of ensemble flood maps exists (NHESS 2023), as does copula-Bayesian probabilistic inundation mapping (NHESS 2024). — [NHESS 23, 2769 (2023)](https://nhess.copernicus.org/articles/23/2769/2023/nhess-23-2769-2023.pdf); [NHESS 24, 2647 (2024)](https://nhess.copernicus.org/articles/24/2647/2024/)

**Parametric insurance**
- **Liberty Mutual × ICEYE parametric wildfire** (15 June 2026):
  - Each insured property is "classified against predefined criteria as either 'destroyed' or 'undamaged'". SAR imagery is "combined with detailed property footprint data and automated machine learning algorithms".
  - Results arrive "within hours" and payouts "within days".
  - Initial markets are the U.S. and Australia.
  - The article does not discuss basis risk or classification uncertainty.
  - [Claims Journal](https://www.claimsjournal.com/news/national/2026/06/15/338222.htm)
- **Radius-based wildfire parametrics:** these "set a specific radius around the covered property, and payout will occur if satellite imagery shows a wildfire burned any part of that area". — [CFO Brew, Jan 2025](https://www.cfobrew.com/stories/2025/01/29/why-more-businesses-are-considering-parametrics-to-protect-against-wildfire-risk) **(snippet-only)**
- **Other contract terms:** WTW structured a TNC parametric wildfire policy in 2025 that incorporates forest-management data. FlashFlood uses a minimum water-depth trigger, and Floodbase runs a municipal flood programme in California. — [WTW, Jul 2025](https://www.wtwco.com/en-us/insights/2025/07/wildfires-combine-parametric-insurance-analytics-and-forest-management-for-improved-risk-transfer); [WEF, Jan 2025](https://www.weforum.org/stories/2025/01/what-is-parametric-insurance-and-how-is-it-building-climate-resilience/) **(snippet-only)**

**EUDR**
- **Timeline:** "The obligation starts on 30 December 2026 for medium and large operators, and on 30 June 2027 for micro and small enterprises." — [Coolset EUDR geolocation guide (secondary)](https://www.coolset.com/academy/eudr-geolocation-requirements-how-to-collect-and-validate-gps-polygon-data-for-your-dds)
- **Information System rules** (secondary sources summarising Commission IS specifications):
  - GeoJSON only, in WGS84.
  - Coordinates are truncated to 6 decimals on storage.
  - Accepted types are Point, MultiPoint, Polygon (≥4 coordinate pairs) and MultiPolygon.
  - Polygons need "no intersecting sides".
  - "Doughnut shapes and figure-eight shapes are rejected".
  - Each DDS is capped at 25 MB, about 30–40k geometries.
  - A single polygon cannot cover several plots.
  - [Kopik / eudrcheck / Coolset summaries](https://eudrcheck.co.uk/guides/eudr-geolocation-requirements-polygons-points-4-hectare-rule); [Coolset](https://www.coolset.com/academy/eudr-geolocation-requirements-how-to-collect-and-validate-gps-polygon-data-for-your-dds) **(secondary; not verified against the Commission's IS technical spec)**
- **Verification and dispute mode:** authorities overlay plots on the JRC 2020 forest-cover baseline at 10 m. "A polygon that leaks beyond the real cultivation boundary can cross into a 10 meter pixel classified as forest at the cut-off date". A valid polygon "can still belong to the wrong producer, omit part of the production land, duplicate another supplier's plot". — [Coolset](https://www.coolset.com/academy/eudr-geolocation-requirements-how-to-collect-and-validate-gps-polygon-data-for-your-dds); [EUDR Navigator](https://eudr-navigator.com/blog/eudr-geolocation-verification-satellite-monitoring) (secondary)
- **Legal text (from r2):** Art. 2(28) requires polygons for plots over 4 ha and sets no positional-confidence requirement. — [EUR-Lex 32023R1115](https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:32023R1115)

### Inferences
- **Demand ranking for (a), updated:**
  1. **Flood (GFM)** and **wildfire perimeters**: likelihoods or uncertainty already exist or are quantified, and outputs are already vector.
  2. **EUDR due-diligence tooling (B2B)**: the core risk is "polygon leaks into a forest pixel", and an inner ring (certain cultivated area) versus outer ring (maximal claim) directly supports risk triage. But the submission format bans holes, so bands must be internal analytics, not the DDS payload.
  3. **Parametric insurance**: the market is moving to footprint-level classification, which shrinks the role of perimeter geometry. Radius-trigger products would benefit from "burn-perimeter outer ring ∩ radius" basis-risk analytics, but no contract I found demands it.
  4. **NWS** (see §1).
- The 20–60 m epistemic buffer on 375 m-class inputs, from Chakravarty, gives a concrete band width for wildfire inner/outer rings.

### Gaps
- I did not find NIFC/WFIGS metadata fields describing perimeter accuracy or method confidence. I did not inspect the schema directly.
- I did not verify the Commission's official EUDR IS GeoJSON technical specification or the post-2025 amendment text (delay, simplification review). The dates and "doughnut rejected" come from secondary compliance vendors.
- I found no concrete EUDR enforcement dispute case over boundary accuracy, since enforcement has not started.
- I found no parametric contract that uses probabilistic or banded perimeters.

---

## 4. Label ambiguity in segmentation: variability statistics, multi-annotator learning, boundary-tolerant metrics, and RL with disagreement bands

### Takeaway
Label ambiguity is large and well-quantified:
- Pairwise Dice is about 0.75–0.79 for skin lesions.
- LIDC-IDRI has four radiologists per nodule.
- On LVIS/COCO, inter-annotator contour distance is "seldom exceeds 1% of image diagonal".

The **boundary-tolerant metrics already set their tolerance from inter-annotator agreement**:
- Surface Dice/NSD sets τ to the **95th percentile of inter-oncologist surface distances**.
- Boundary IoU sets d to 2% of the image diagonal because two experts reach median Boundary IoU > 0.9 at that d.

Ambiguous-segmentation generators (Prob. U-Net → PHiSeg → diffusion → Schrödinger bridge) are mature, with GED on LIDC improved from 0.310 to 0.208.

**RL verdict: I found no RL/GRPO segmentation method whose reward uses annotator-disagreement bands or a "no penalty inside the band" rule.** Current RL segmentation rewards (Seg-Zero, SAM-R1, CORE-Seg) use IoU/Dice/L1 thresholds against a single ground truth. Proposal (b) is therefore **novel as an RL reward**, though it is conceptually the RL analogue of NSD-at-tolerance and of soft-label/multi-rater losses.

### Cited Findings
**Inter-annotator variability statistics**
- **LIDC-IDRI** has 875 patients and 1,575 nodules "each with complete four-expert annotations". **RIGA** has 750 fundus images, with optic disc and cup annotated by **six** ophthalmologists. **QUBIQ** brain-growth has 34 training and 5 validation subjects, with **seven** annotations per slice. QUBIQ evaluation compares averaged predictions to averaged labels with continuous Dice. — [search summary of multi-rater literature: Ji et al. CVPR 2021 and related](https://www.researchgate.net/publication/355882110_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-rater_Agreement_Modeling) **(snippet-only aggregation; counts are standard but not re-verified against dataset papers)**
- **Skin lesions (IMA++, arXiv 2508.09381, Aug 2025):**
  - 2,394 images, 5,111 masks, 15 annotators.
  - Pairwise Dice is **0.791 ± 0.215 for benign** and **0.753 ± 0.227 for malignant** lesions (p < 0.01).
  - Using agreement prediction as an auxiliary task gave "+4.2% balanced accuracy" for diagnosis.
  - [arXiv 2508.09381](https://arxiv.org/html/2508.09381v1)
- **Natural images (Boundary IoU paper, Cheng et al., CVPR 2021):**
  - "A classic study… shows that the pixel distance between two contours of the same object labeled by different annotators seldomly exceeds 1% of the image diagonal, irrespective of the object size". This is confirmed on LVIS double annotations, where "boundary ambiguity is fixed and independent of objects area".
  - On LVIS-collected double expert annotations for COCO and ADE20k images, "median Boundary IoU between the annotations of the two experts exceeds 0.9… when d equals 2% of an image diagonal (15 pixels distance on average)".
  - For Cityscapes, d is 0.5% of the diagonal.
  - Rule given: "the annotation consistency sets the lower bound on d".
  - [arXiv 2103.16562 (PDF text)](https://arxiv.org/abs/2103.16562)
- **Building footprints:** a 2026 Remote Sensing building-extraction paper notes "no independent field-survey-based accuracy assessment or formal inter-annotator agreement analysis was available for the reference labels". A 2026 open building-footprint benchmark uses a pairwise IoU matching threshold of only **0.30**. — [Remote Sens. 18(17):2971](https://doi.org/10.3390/rs18172971); [Remote Sens. 18(17):2895](https://doi.org/10.3390/rs18172895) **(snippet-only)**
- **Caution on raw IoU as agreement:** chance-corrected agreement for vision tasks is argued in KαLOS (arXiv 2603.27197), which claims raw mean IoU "reads around 0.95 even for annotators who never looked at the image" on typical detection corpora. — [arXiv 2603.27197](https://arxiv.org/html/2603.27197) **(snippet-only)**

**Methods that learn from multiple annotators / output sample sets**
- **LIDC-IDRI benchmark** (Table 1 of Diffusion Schrödinger Bridge for ambiguous segmentation, arXiv 2509.17187, Sep 2025, also listed at MICCAI 2025). GED is lower-is-better.

  | Method | GED | D_max | CI | NFE |
  |---|---|---|---|---|
  | Prob. U-Net | **0.310** | 0.907 | 0.831 | – |
  | PHiSeg | **0.270** | 0.904 | 0.736 | – |
  | CIMD (diffusion) | **0.318** | 0.841 | 0.883 | 1000 |
  | CCDM | **0.264** | 0.886 | 0.918 | 1000 |
  | SSB | **0.245** | 0.925 | 0.945 | 50 |
  | SSB++ | **0.208** | 0.943 | 0.958 | 50 |

  SSB++ reduces GED by 21.2% relative to CCDM. — [arXiv 2509.17187](https://arxiv.org/html/2509.17187); [MICCAI 2025 paper 4859](https://papers.miccai.org/miccai-2025/paper/4859_paper.pdf)
- **Lineage:** Probabilistic U-Net (Kohl et al. 2018); Stochastic Segmentation Networks; CIMD (Rahman et al., CVPR 2023); 2025–2026 diffusion variants. — [arXiv 1806.05034](https://arxiv.org/pdf/1806.05034); [arXiv 2006.06015](https://arxiv.org/pdf/2006.06015); [CVPR 2023 CIMD](https://openaccess.thecvf.com/content/CVPR2023/papers/Rahman_Ambiguous_Medical_Image_Segmentation_Using_Diffusion_Models_CVPR_2023_paper.pdf); [arXiv 2504.05977](https://arxiv.org/pdf/2504.05977); [arXiv 2603.04024 Volumetric Directional Diffusion](https://arxiv.org/html/2603.04024); [arXiv 2603.11550 PCA-Enhanced Prob. U-Net](https://arxiv.org/pdf/2603.11550)
- **Foundation-model variants:** Probabilistic SAM (arXiv 2509.05809) captures annotator variability on LIDC. Multi-rater prompting (arXiv 2404.07580) treats raters as domains. SPA (Zhu, Wu, Ouyang, Kamnitsas, Noble; arXiv 2411.15513, rev. Jul 2025) presents distinct candidates and adapts to "variability in labeling preferences on ambiguous pixels" via user feedback. — [arXiv 2509.05809](https://arxiv.org/pdf/2509.05809); [arXiv 2404.07580](https://arxiv.org/pdf/2404.07580); [arXiv 2411.15513](https://arxiv.org/abs/2411.15513)
- **Soft-label / agreement-modelling losses:** Ji et al., CVPR 2021, model multi-rater agreement for calibrated segmentation. Soft labels from averaged raters encode pixel-level difficulty. — [CVPR 2021 Ji et al.](https://openaccess.thecvf.com/content/CVPR2021/papers/Ji_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-Rater_Agreement_Modeling_CVPR_2021_paper.pdf); [Zhang et al. 2023, PMC9966563](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9966563/)

**Boundary-tolerant metrics**
- **Surface Dice / NSD** (Nikolov et al., DeepMind, arXiv 1809.04430):
  - It measures overlap of surfaces "at a specified tolerance τ".
  - "We computed these acceptable tolerances for each organ by measuring the inter-observer variation in segmentations between three different consultant oncologists". τ is set as "the 95th percentile of the distances collected across multiple segmentations".
  - It works "by strongly penalising border placement outside a tolerance determined by consultant oncologists".
  - **This is exactly the "no penalty inside the band, penalty outside" rule, used as a metric.**
  - [arXiv 1809.04430](https://arxiv.org/abs/1809.04430)
- **Boundary IoU:** IoU over the d-wide boundary strips. It is more sensitive than Mask IoU for large objects and does not over-penalise small ones. It underpins Boundary AP and Boundary PQ. The authors compared it against mF-measure (BF-score-style) with d from 0.1% to 2.1% of the diagonal. — [arXiv 2103.16562](https://arxiv.org/abs/2103.16562)

**RL and rewards in segmentation (2025–2026)**
- **Seg-Zero** uses binary tolerance rewards against a single ground truth: bbox IoU > 0.5 → 1, bbox L1 < 10 px → 1, point L1 < 100 px → 1. — [arXiv 2503.06520](https://arxiv.org/html/2503.06520)
- **SAM-R1** uses SAM inside the RL loop as the reward generator (IoU-based). — [arXiv 2505.22596](https://arxiv.org/pdf/2505.22596)
- **Other RL segmentation work:** Seg-R1, GenSeg-R1 and VideoSeg-R1 also use IoU/Dice-style rewards. — [arXiv 2506.22624](https://arxiv.org/html/2506.22624); [arXiv 2602.09701](https://arxiv.org/pdf/2602.09701); [arXiv 2511.16077](https://arxiv.org/html/2511.16077)
- **CORE-Seg (2026, GRPO, complex lesions):** the reward is format + Hungarian bbox + a "dual-granularity mask reward": "if Dice < τ (0.05)… GIoU, otherwise 1 + λ·Dice". There is **no uncertainty, multi-annotator or soft-label term**. — [arXiv 2603.05911](https://arxiv.org/html/2603.05911)
- **MedSAM-Agent:** multi-turn agentic RL for interactive medical segmentation. I found no evidence of a disagreement-band reward. — [arXiv 2602.03320](https://arxiv.org/pdf/2602.03320) **(not read in full)**
- **Older RL work:** boundary-aware rewards (global + boundary reward) and uncertainty kept in the MDP *state* exist in interactive 3D RL segmentation, but they are not annotator bands. — [arXiv 2303.10692](https://arxiv.org/pdf/2303.10692); [CVPR 2020 Liao et al.](https://openaccess.thecvf.com/content_CVPR_2020/papers/Liao_Iteratively-Refined_Interactive_3D_Medical_Image_Segmentation_With_Multi-Agent_Reinforcement_Learning_CVPR_2020_paper.pdf)

### Inferences
- **Proposal (b) has strong precedent for its *rule* but none for its *use*.**
  - NSD is the metric form of a band-tolerant reward, with τ derived from inter-observer variation.
  - Boundary IoU derives d from two-expert agreement.
  - Seg-Zero-style rewards already use fixed tolerance thresholds (10 px, 100 px).
  - Nobody I found replaces the fixed threshold with a **spatially varying band from annotators, ensembles or conformal sets**.
  - Natural formulation: reward = NSD-like score where the per-pixel tolerance is the local band width (outer minus inner ring), or a penalty computed only on predicted boundary points outside the [inner, outer] band.
- A band can be built without multiple annotators:
  - from a sample generator (PHiSeg, diffusion, SSB at 50 NFE);
  - from deep ensembles;
  - from conformal inner/outer sets (§5).
  This matters for remote sensing, where multi-annotator datasets such as building footprints are basically absent.
- Default band widths, if no annotator data exists: about 1–2% of the image diagonal (natural images), 0.5% (Cityscapes-like high-quality data), and organ-specific 95th-percentile distances (medical). Each is a citable prior.

### Gaps
- **BF score** (Csurka et al. 2013) parameters and the MATLAB default tolerance were not fetched. They are covered only indirectly via Boundary IoU's mF-measure comparison.
- I did not fetch a remote-sensing **building-footprint inter-annotator study** with numbers. I found none, and the 2026 papers explicitly say none was available for their labels. This is a real gap in the literature, and an opportunity.
- I did not retrieve COCO-specific boundary-quality numbers (e.g., LVIS paper statistics on COCO polygon coarseness) beyond the Boundary IoU passages.
- The RL-with-disagreement-band search covered arXiv/GRPO-era papers through Oct 2026 via web search only. It is not exhaustive, so a workshop paper could exist.

---

## 5. Conformal prediction for segmentation, 2024–2026: inner/outer guarantees, polygon outputs, and use as a training signal

### Takeaway
**Verdict:** the conformal segmentation literature added several variants in 2025–2026: spatially-aware CONSIGN, FDR control for instances, industrial defects, and radiomics. Conformal **risk training** (NeurIPS 2025) now differentiates through conformal risk control during training, and ConformalSAM uses conformal filtering to produce training pseudo-labels.

**Still, I found no conformal segmentation method that outputs polygons.** I also found none that uses calibrated inner/outer *bands* as a reward or tolerance region for training. Training-time uses control a scalar risk such as FNR, or filter pseudo-labels.

### Cited Findings
- Inner/outer and morphological conformal sets from the r2 note, still the core: Davenport (arXiv 2410.03406); CONSEMA, MICCAI 2025 (arXiv 2503.05618); false-positive control with nested shrunken masks (arXiv 2511.15406); instance conformal sets with a threshold sweep, applied to agricultural fields (arXiv 2602.10045). — [arXiv 2410.03406](https://arxiv.org/pdf/2410.03406); [arXiv 2503.05618](https://arxiv.org/abs/2503.05618); [arXiv 2602.10045](https://arxiv.org/abs/2602.10045)
- **CONSIGN** (arXiv 2505.14113): "applying CP directly to image segmentation ignores the spatial correlations between pixels". It incorporates spatial groupings via decomposition. — [arXiv 2505.14113](https://arxiv.org/abs/2505.14113)
- **FDR control** in medical *instance* segmentation via conformal risk control. — [arXiv 2504.04482](https://arxiv.org/pdf/2504.04482)
- **Industrial surface defects**, with statistical guarantees. — [arXiv 2504.17721](https://arxiv.org/pdf/2504.17721)
- **ConRad** (Jul 2026), conformal prediction for radiomics. — [arXiv 2607.08084](https://arxiv.org/pdf/2607.08084)
- **COMPASS** (Sep 2025), conformal intervals on segmentation *metrics*. — [arXiv 2509.22240](https://arxiv.org/abs/2509.22240)
- **Conformal Risk Training** (Yeh et al., NeurIPS 2025, arXiv 2510.08748):
  - "the first to incorporate conformal risk control differentiably during model training".
  - It controls OCE risks, including expected loss and CVaR.
  - Demonstrated on **controlling false-negative rate in tumor image segmentation**.
  - It notes that post-hoc CRC "can degrade average-case performance due to its lack of feedback to the model".
  - [arXiv 2510.08748](https://arxiv.org/abs/2510.08748); [NeurIPS 2025 PDF](https://proceedings.neurips.cc/paper_files/paper/2025/file/6559542f75b4452ebaaf82094c7defb7-Paper-Conference.pdf); [code](https://github.com/chrisyeh96/conformal-risk-training)
- **Conformal Risk Control foundation** (Angelopoulos et al.; v4 Jun 2025), including the monotone-loss framework used for FNR-controlled segmentation. — [arXiv 2208.02814](https://arxiv.org/pdf/2208.02814)
- **ConformalSAM** (arXiv 2507.15803) uses conformal prediction on a foundation segmentation model to support **semi-supervised** semantic segmentation, i.e. conformally filtered pseudo-labels as a training signal. — [arXiv 2507.15803](https://arxiv.org/pdf/2507.15803) **(title/listing only; mechanism detail from title, not read in full)**

### Inferences
- The missing pieces are exactly (a) and (b):
  - **(a)** No paper polygonizes conformal inner/outer masks. The egg-yolk encoding in §2, plus a nesting-safe tracer, gives "conformal polygons" with a stated coverage.
  - **(b)** Conformal risk training shows that gradients through calibrated sets are feasible. The step nobody has taken is to treat the calibrated [inner, outer] band as the **zero-penalty region of an RL reward** (or of a loss), so the policy is not punished for boundary placement within calibrated ambiguity.
- **Caveat for (b):** conformal bands are calibrated on the *model's* errors, not on *label* ambiguity. Using them as a reward tolerance risks a feedback loop, because a worse model gets a wider band and so a more lenient reward. Bands should be fixed from annotators or an independent ensemble, or re-calibrated with a frozen reference model. This is an inference, not a published finding.

### Gaps
- I did not read ConformalSAM or CONSIGN in full, so their exact outputs (pixel sets vs. regions) are unverified beyond abstracts and titles.
- I found no conformal method for remote-sensing *polygon* products, beyond the agricultural-field instance sets of arXiv 2602.10045.
- I found no evaluation of conformal coverage under *annotator* variability as distinct from a single reference label.

---

## 6. Overall verdict for proposals (a) and (b)

### Takeaway
- **(a) Confidence-banded polygons: novel as a deliverable, well-founded in theory, with real but mostly latent demand.**
  - The theory exists: egg-yolk, broad boundaries, conformal inner/outer sets, contour-boxplot bands.
  - No standard or product ships banded vectors. GeoJSON has no uncertainty concept, GFM ships likelihood as a raster, and NIFC/FIRIS perimeters are crisp.
  - The strongest near-term users are flood (GFM) and wildfire perimeter producers, plus EUDR due-diligence analytics. NWS PHI/TIM is still not operational.
- **(b) Uncertainty-tolerant RL rewards: no prior RL implementation found.**
  - The tolerance principle is validated by NSD, whose τ is the 95th percentile of inter-observer distance, and by Boundary IoU, whose d comes from two-expert agreement.
  - Current GRPO segmentation rewards (Seg-Zero, SAM-R1, CORE-Seg) use fixed IoU/Dice/pixel-distance thresholds against a single ground truth.

### Cited Findings
- NSD tolerance from inter-observer variation. — [arXiv 1809.04430](https://arxiv.org/abs/1809.04430)
- Boundary IoU d from expert agreement. — [arXiv 2103.16562](https://arxiv.org/abs/2103.16562)
- RL rewards are single-GT with fixed thresholds. — [Seg-Zero](https://arxiv.org/html/2503.06520); [CORE-Seg](https://arxiv.org/html/2603.05911)
- GeoJSON has no concept of uncertainty. — [Geo URI scheme](https://en.wikipedia.org/wiki/Geo_URI_scheme)
- Operational wildfire perimeter forecasts "lack uncertainty quantification". — [arXiv 2510.09666](https://arxiv.org/pdf/2510.09666)
- TIM is still on a "(potential) path toward NWS operations". — [ADS 2024AMS...10429803B](https://ui.adsabs.harvard.edu/abs/2024AMS...10429803B/abstract)

### Inferences
- The defensible novelty claim: **"band-tolerant rewards whose tolerance is spatially varying and derived from multi-annotator, ensemble or conformal bands, generalising NSD's fixed τ to a per-location band, used inside GRPO-style segmentation RL."**
- Recommended validation datasets with multi-rater bands: LIDC-IDRI (4 raters), RIGA (6), QUBIQ (7), IMA++ (2–5, 15 annotators). LVIS double annotations serve natural images. Remote sensing lacks a multi-annotator building or field dataset, so ensembles or conformal sets would have to substitute.

### Gaps
- I found no willingness-to-pay or explicit customer requirement for vector uncertainty bands from any operator (GFM, NIFC, insurers, EUDR authorities).
- I found no empirical comparison of band-tolerant vs. strict rewards, in RL or in supervised losses, on multi-rater datasets.
