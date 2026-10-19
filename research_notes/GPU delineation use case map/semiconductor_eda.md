# Semiconductor / EDA / Electronics: where raster <-> (rectilinear) polygon conversion matters (2023–2026)

Scope note: the user's engine converts a binary occupancy bitmap into closed RECTILINEAR (Manhattan) polygon rings using one scanline pass with 2x2 pixel-type classification, on CUDA. These notes map where that operation, or its inverse (rasterization), sits in chip and electronics manufacturing flows. Research date: 2026-10-07. Items before 2023 are marked [BACKGROUND].

---

## Q1. Computational lithography (ILT, curvilinear masks, OPC, cuLitho): pixel mask to polygon conversion, data volume, runtimes

### Takeaway
Leading-edge ILT is computed in the pixel domain on GPUs (cuLitho, D2S, NVIDIA research), but mask data still has to leave as polygons or curves (OASIS / the MULTIGON extension) for MDP, MRC and mask writing. That pixel-to-geometry step decides how big the data gets: a flat full reticle of curvilinear data is estimated at roughly 8 to 440 TB as piecewise-linear polygons, depending on segment length. Public cuLitho material talks about accelerating OPC/ILT simulation and says nothing about contour or polygon extraction, so that step looks under-served in public sources.

### Cited Findings
- **cuLitho production (March 18, 2024):** TSMC and Synopsys put NVIDIA cuLitho into production. Joint tests showed a **45x speedup on curvilinear flows** and **close to 60x on Manhattan-style flows**. A generative-AI step that produces a "near-perfect inverse mask" adds another **2x**, after which "traditional methods" derive the final mask. **350 NVIDIA H100 systems replace 40,000 CPU systems.** A typical mask set takes "30 million or more hours of CPU compute time." The release says nothing about mask data prep or contour handling. — [NVIDIA Newsroom](https://nvidianews.nvidia.com/news/tsmc-synopsys-nvidia-culitho); [Tom's Hardware](https://www.tomshardware.com/pc-components/cpus/nvidias-generative-ai-tool-delivers-a-radical-60x-performance-boost-for-chipmakers-tsmc-and-synopsys-are-now-using-the-culitho-software-in-production)
- **NVIDIA/IMEC paper (arXiv 2602.15036, submitted Jan 27, 2026; Mukhopadhyay, ..., Werner Gillijns, Ryan Ryoung-han Kim, Kumara Sastry, Vivek Singh):** reports **57x end-to-end acceleration**, a **35% process-window improvement** and **19% lower EPE**, with silicon validation at IMEC. The authors describe it as "the first quantified chip-scale demonstration" of AC+AI litho benefits. The primitives they rebuilt include "diffractive optics, computational geometry, multi-variant optimization, data processing," and the work covers curvilinear masks and high-NA EUV. The abstract does not mention polygon or contour extraction. — [arXiv 2602.15036](https://arxiv.org/abs/2602.15036)
- **ILT is pixel-native (D2S, JM3 Vol. 23(1), Jan–Mar 2024; Shendre, Fujimura, Niewczas, Kronmiller):** "curvilinear ILT, practically speaking, is all computed in the pixel domain." Even tools that appear to move polygon edges rasterize the mask, simulate the wafer in the Fourier domain, and then **contour** the result to compare against the target. — [D2S JM3 paper PDF](https://design2silicon.com/wp-content/uploads/2023/12/JM3-23048SS_pixel_based_computing.pdf)
- **Data-volume numbers from the same paper (flat full reticle, 104 x 132 mm, circle-array upper bound, 2 bytes/vertex):**
  - With **1 nm segments**: ~160 to **~440 TB**. With **5 nm**: ~32–88 TB. With **10 nm**: ~16–44 TB. With **20 nm**: ~8–22 TB. Values depend on the circle diameter/pitch (100/150, 100/200, 150/225, 150/300 nm).
  - An implicit-Bezier curvilinear format brings this to an estimated **10–30 TB**, about a 10x reduction.
  - A **raw 16 nm pixel format** gives about **50 x 10^12 pixels, roughly 100 TB uncompressed**, and that figure is constant whatever the design complexity.
  - A "pixel-based-computing-aware" curvilinear format could fall below 10 TB.
  — [D2S JM3 paper](https://design2silicon.com/wp-content/uploads/2023/12/JM3-23048SS_pixel_based_computing.pdf)
- **Key theoretical claim in that paper:** contours and pixel dose maps are "duals" at a given resolution limit. The maximum useful vertex count is about **contour perimeter / pixel size**. The multibeam pixel is about **16 nm on mask (≈4 nm on wafer)**, so "we do not need 1 nm segments" and a sampling interval near the writer pixel size is enough ("pixel dose equivalence"). — [D2S JM3 paper](https://design2silicon.com/wp-content/uploads/2023/12/JM3-23048SS_pixel_based_computing.pdf)
- **MULTIGON:** a native curved-geometry record that uses Bezier curves. It was developed by the Curvilinear Working Group and added as the **SEMI P49 extension of OASIS**. Its benefit depends on how many control points OPC and MPC emit, because "high control point count increase[s] runtime complexity and data volume which affects the performance of downstream processes like fracture." A 2025 SPIE paper (vol. 13687) studies MDP performance against Multigon data volume. There is also a 2023 SPIE paper on a multibeam fracture flow with Multigon input. — [SPIE 13687 136871F](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/13687/136871F/A-study-of-Multigon-mask-data-preparation-performance-and-its/10.1117/12.3071765.full); [SPIE 12751 Multibeam fracture with Multigon](https://www.spiedigitallibrary.org/conference-proceedings-of-spie/12751/127510A/Multibeam-fracture-flow-based-on-multigon-format-input/10.1117/12.2687753.short); [ADS: Multigon curvy mask manufacturability, SPIE 13655 (2025)](https://ui.adsabs.harvard.edu/abs/2025SPIE13655E..11B/abstract)
- **[BACKGROUND 2021]** SemiEngineering:
  - EUV data volume per mask is about **5x** that of optical.
  - Curvilinear data causes "a data volume explosion in the SEMI P44 (OASIS) layout format" (Peter Buck, Siemens).
  - A 2019 working group of Samsung, TSMC and Intel looked at quadratic Bezier, B-spline, polygon simplification and curvature-based fragmentation, with **2x–8x compaction vs. OASIS P44**.
  - Aki Fujimura (D2S): "an overwhelming majority of the shapes on the mask are Manhattan shapes."
  — [SemiEngineering, "The Quest For Curvilinear Photomasks"](https://semiengineering.com/the-quest-for-curvilinear-photomasks/)
- **Manhattanization is still a real step for VSB writers:** "Mask rule checkers may perform Manhattan conversion to convert a curvy ideal mask from ILT to a zigzag polygon pattern... restricting output pattern edges to be either horizontal or vertical to accommodate an e-beam mask writer." Freeform patterns turned into polygons with many jogs raise e-beam shot count and mask data volume. After Manhattanization, mask shapes change a lot and need re-optimization. — [USPTO patent 9495507 "Method for integrated circuit mask patterning"](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/9495507); [JM3 "ILT: 30 years from concept to practical, full-chip reality" (2021, BACKGROUND)](https://www.spiedigitallibrary.org/journals/Inverse-lithography-technology--30-years-from-concept-to-practical/volume-20/issue-03/030901/Inverse-lithography-technology--30-years-from-concept-to-practical/10.1117/1.JMM.20.3.030901.full)
- **Industry sentiment (12th eBeam Initiative Luminaries survey, July 2023, 47 companies):**
  - 87% say leading-edge mask shops can handle at least a limited number of curvilinear masks.
  - 75% say curvilinear ILT is needed for 2 nm / 0.33 NA EUV.
  - 95% say multi-beam writers are needed for EUV masks.
  — [SPIE 13177 131770A](https://proceedings.spiedigitallibrary.org/conference-proceedings-of-spie/13177/131770A/12th-eBeam-initiative-survey-reports-semiconductor-industry-luminaries-are-confident/10.1117/12.3031849.full); [Silicon Semiconductor](https://siliconsemiconductor.net/article/117819/Luminaries_like_high-NA_EUV_and_curvilinear_photomasks)

### Inferences
- ILT runs on GPUs in the pixel domain, mask writers write pixels, and the formats in between (OASIS/MULTIGON) are vector. That leaves two conversions in every curvilinear flow: pixel to contour after ILT, and contour to pixel in the writer. NVIDIA, Synopsys and TSMC publicize the simulation and optimization speedups but say nothing public about the contour/polygon extraction stage. That stage is a candidate bottleneck under Amdahl's law once simulation is 40–60x faster. (Inference; not confirmed by a source.)
- **Fit for the user's engine (Manhattan output from a binary bitmap):**
  - Direct fit: (a) Manhattanization of pixel ILT for **VSB-written layers and masks**, which per Fujimura are still the large majority of mask shapes; (b) fast "pixel-exact" polygon export where D2S's pixel-dose-equivalence argument applies. A staircase at the writer pixel pitch is by construction pixel-dose-equivalent to the binary raster.
  - Mismatch: production curvilinear ILT wants smooth curves or Multigon, and gray-level (sub-pixel) edges, not binary staircases. A binary staircase at 16 nm mask pixels also produces about one vertex per pixel of perimeter, which matches D2S's upper bound but is far from the most compact form.
  - An intellectually stronger angle: use the rectilinear ring extraction as the fast topology/first stage (connected components, holes, ring ordering), then run a GPU curve fit (Bezier/Multigon) on each ring. Alternatively, extend the 2x2 classifier to gray-level / marching-squares sub-pixel edges.

### Gaps
- No public runtime or share-of-total figures for the contour/polygon-extraction step inside cuLitho, Synopsys Proteus/S-Litho, Siemens Calibre or ASML Tachyon flows.
- Full text of the SPIE 2025 Multigon MDP-performance paper (13687-1F) could not be retrieved, so its specific numbers are missing.
- No source found on whether TSMC's production cuLitho flow emits Manhattan or curvilinear output for specific layers.

---

## Q2. Mask data preparation (MDP), MRC, multi-beam mask writers: rasterization vs vectorization

### Takeaway
Multi-beam mask writers (IMS Nanofabrication, NuFlare) rasterize all input into dosed pixels, and pixel-level corrections such as PLDC/MPC already run on GPUs inline. Upstream, MDP, fracture and MRC still work on vector geometry, and curvilinear MRC is an open problem.

### Cited Findings
- In multi-beam writers, "input shapes are rasterized to arrays of equal sized pixels that are each independently dosed." Pixel-level dose correction (PLDC) is done in the pixel domain with GPU acceleration, offline or inline. NuFlare's MBM-1000 has inline GPU PLDC, and NuFlare built inline pixel-based MPC into the writer "with zero additional turn-around-time." — [SPIE 10454 GPU-accelerated inline PLDC (2017, BACKGROUND)](https://ebooks.spiedigitallibrary.org/conference-proceedings-of-spie/10454/104540D/GPU-accelerated-inline-linearity-correction--pixel-level-dose-correction/10.1117/12.2281922.full); [SemiEngineering, Improving Uniformity and Linearity](https://semiengineering.com/improving-uniformity-and-linearity-for-all-masks)
- **[BACKGROUND]** The IMS multi-beam writer (2016) has **262,144 beamlets** and a write time of about **12 hours or less regardless of complexity**, against 24+ hours for single-beam VSB on complex masks. The roadmap named the IMS MBMW-301 (50 keV, 3 nm node) and the NuFlare MBM-2000. — [SemiEngineering](https://semiengineering.com/the-quest-for-curvilinear-photomasks/)
- D2S: "Most advanced masks are written by MBM writers and are therefore written using pixels." The writer pixel size fixes the smallest printable feature. Different circles rasterize to the same pixel data, so geometry finer than the pixel is wasted. Contour-based EPE checks at "infinite resolution" are "somewhat insufficient and misleading in 2D curvy region." — [D2S JM3 2024](https://design2silicon.com/wp-content/uploads/2023/12/JM3-23048SS_pixel_based_computing.pdf)
- Bezier curves bring practical problems in MDP (SPIE 13216, 2024), and mask process correction is being improved on Bezier-represented curvilinear data. — [SPIE 13216 "Practical challenges of Bezier curves in MDP"](https://neurophotonics.spiedigitallibrary.org/conference-proceedings-of-spie/13216/1321621/Practical-challenges-of-B%C3%A9zier-curves-in-mask-data-preparation/10.1117/12.3034563.short); [SPIE Photomask 13216-172](https://spie.org/photomask-technology/presentation/Enhancing-mask-process-correction-on-curvilinear-data-with-B%C3%A9zier-curve/13216-172)

### Inferences
- The machine's internal data path is raster, so the vector stage (OASIS/MULTIGON) is an interchange and verification format. A bitmap-to-polygon converter fits naturally in **reverse/verification** paths: re-deriving geometry from the writer's pixel dose map or from a pixel-domain MPC result for MRC, diffing, or archiving. It also helps **VSB legacy writers**, which still need rectangles/trapezoids. Rectilinear rings decompose straight into VSB rectangle shots, a natural follow-on kernel.
- Pixel-based MRC (width/space checks via morphology on the raster) is an alternative to vector MRC. CurvyILT (Q5) already uses differentiable morphology for MSA/MSD rules. Vector MRC would still need polygons.

### Gaps
- No public 2024–2026 specs found for current-generation writer data paths (IMS MBMW-301/401, NuFlare MBM-3000/4000), such as input format, rasterizer hardware, or throughput in GB/s.
- No source found that quantifies MRC runtime on curvilinear full-chip data.

---

## Q3. Wafer/mask inspection and metrology: SEM contour extraction, die-to-database, "SEM contour to GDS"

### Takeaway
SEM images are routinely turned into contours: Siemens Calibre SEMSuite, with ML/GPU-accelerated contouring shown with IMEC at SPIE 2024, and KLA's die-to-database (D:DB) modeling. The contours feed OPC/etch model calibration and are compared with layout. Two steps are hard: aligning noisy SEM images with clean CAD, and extracting sub-pixel curved edges. Binary rectilinear vectorization only partly fits here.

### Cited Findings
- Siemens EDA and IMEC presented "Machine learning (ML) based SEM contour extraction accelerated by GPU for etch modeling application" at SPIE Advanced Lithography 2024. The KB article requires a login, so details were not retrieved. — [Siemens Support KB000130564](https://support.sw.siemens.com/en-US/okba/KB000130564_EN_US/Machine-learning-ML-based-SEM-contour-extraction-accelerated-by-GPU-for-etch-modeling-application/index.html)
- Calibre SEMSuite extracts contours from SEM images to build OPC models for the contact layer, with CD-SEM anchor-point measurements used to calibrate (SPIE 13423, 2025). — [SPIE 13423 134230E](https://journals.spiedigitallibrary.org/conference-proceedings-of-spie/13423/134230E/SEM-contour-extraction-application-on-OPC-model-of-contact-layer/10.1117/12.3052691.full)
- A full-chip OPC model calibration flow based on contours was published at SPIE 2024 (12954). — [SPIE 12954 129540K](https://journals.spiedigitallibrary.org/conference-proceedings-of-spie/12954/129540K/A-novel-flow-of-full-chip-OPC-model-calibration-and/10.1117/12.3010034.full)
- KLA's D:DB image modeling is used for accurate contour extraction (patent "Image modeling-assisted contour extraction"). — [USPTO 12444628](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12444628)
- Aligning layout with SEM is "a crucial task for die-to-database inspection." The CAD style and the noisy grayscale SEM differ so much that plain cross-correlation often misaligns; GAN-based alignment has been proposed. — [ResearchGate: Precise Pattern Alignment for D2DB Inspection based on GAN](https://www.researchgate.net/publication/363047059_Precise_Pattern_Alignment_for_Die-to-Database_Inspection_Based_on_the_Generative_Adversarial_Network)
- In research fabs (photonics/MEMS), custom Python SEM-vs-GDS pipelines do preprocessing, contour extraction, SEM-to-GDS registration and error analysis (arXiv 2609.10741). — [arXiv 2609.10741](https://arxiv.org/pdf/2609.10741)

### Inferences
- In metrology, the contour has to be sub-pixel and curved, so the binary Manhattan engine is not the right primitive. Where it can fit:
  - Fast **polygonization of segmentation masks** from ML SEM models, after which the polygons are compared with GDS by Boolean XOR.
  - **Rasterize-GDS / vectorize-SEM** round trips at large scale (die-to-database over full dies).
  - Metal/via layers that are Manhattan by design, where snapping SEM segmentation back to rectilinear rings is a meaningful regularization ("re-Manhattanize the measured layout").

### Gaps
- No public throughput numbers found for commercial SEM contouring (contours per second, GPU vs CPU).
- No details found on how Applied Materials, ASML-HMI or KLA represent extracted contours internally.

---

## Q4. IC reverse engineering (delayering SEM images to GDS/netlist; hardware Trojan detection)

### Takeaway
This is the closest structural match to the user's engine. Academic pipelines segment SEM images of metal layers into binary track/via masks and then **vectorize them into polygons with CPU OpenCV**. The result is compared with GDSII to find Trojans. Data runs to thousands to hundreds of thousands of images per chip. Manual cleanup and stitching across images are the stated pain points.

### Cited Findings
- **Unsupervised SEM segmentation for IC layout extraction** (Journal of Cryptographic Engineering 2025, 15:21; author affiliations not verified in my read):
  - The dataset has **648 SEM images** (327 M1, 321 M2) of a **commercial 180 nm IC** at **14.65 nm/pixel** with 10% overlap. M2 is split into **15,408 patches** of 512x512.
  - Labels were made by automatic polygon labeling plus manual correction, averaging **6 min per image**.
  - The evaluation binarizes segmentation masks and then **"extract[s] track and label polygons using OpenCV"**, drops polygons under 35 pixels, and checks electrical errors with polygon intersection.
  - Stated limitations: images are processed in isolation, and tracks spanning multiple images are not stitched. "IC layout extraction from SEM images is error prone" and needs manual review.
  — [J. Cryptographic Engineering 2025 (open PDF)](https://d-nb.info/1391004175/34); [Springer page](https://link.springer.com/article/10.1007/s13389-025-00385-5)
- Prior supervised work cited there: a CNN trained on 21 SEM images of **8192x8192 pixels**, reaching 91.86% mIoU (Yu et al.), and gate annotation on 640 poly-layer images (Cheng et al.). — [same PDF](https://d-nb.info/1391004175/34)
- **"Red Team vs. Blue Team"** (IEEE S&P 2023; Puschner, Moos, Becker, Kison, Moradi, Paar; RUB/MPI-SP/UCLouvain/BKA):
  - Trojan detection on 90/65/40/28 nm ICs by **GDSII-vs-SEM-image comparison**.
  - All modifications were found at 90/65/40 nm, and **37/40 at 28 nm**. About **500 false positives** came from **more than 1.5 million standard cells**.
  - One composite image of the 65 nm chip was made from **4,225 SEM images**.
  - Pain points: superimposing and registering the images precisely ("quite a challenge"), and dust or hair hiding whole rows of cells at 28 nm.
  - The data and algorithms are public.
  — [RUB HGI news](https://hgi.rub.de/en/news/newsarchiv/hginews/detecting-manipulations-in-microchips); [IACR ePrint 2022/1720](https://eprint.iacr.org/2022/1720); [Dataset: back-side SEM images and modified GDSII](https://b2find.eudat.eu/dataset/5278c019-5a2e-5e97-afbf-cd1fde625cac)
- Other work: U-Net-based die RE for trusted microelectronics (IOP 2024), and federated learning for hardware assurance (arXiv 2604.20020, 2026). — [IOP Eng. Res. Express](https://iopscience.iop.org/article/10.1088/2631-8695/ad7c06/ampdf); [arXiv 2604.20020](https://arxiv.org/pdf/2604.20020)

### Inferences
- **Strongest niche for the user's engine:** metal layers are Manhattan by design. Segmentation produces binary masks with thousands of 8k x 8k tiles per layer, and the vectorization step today is OpenCV `findContours` on CPU followed by Python polygon operations.
  - A CUDA scanline Manhattan-ring extractor that runs straight on GPU segmentation output, with zero copy and hole-aware closed rings, would give **SEM to GDS polygons** (OASIS/GDS export) at line rate.
  - Natural extensions: cross-tile ring stitching (the gap the JCE paper names explicitly), rectilinear snapping/regularization to the design grid, and GPU XOR against the golden GDS for Trojan diffing.
- Hardware assurance (DoD trusted foundry and the EU Chips Act security programs) gives this a policy-relevant framing. That is inference and was not researched here.

### Gaps
- No published throughput or bottleneck breakdown (segmentation vs vectorization vs stitching) for full-chip RE pipelines.
- Commercial RE vendors (e.g., TechInsights) and their tooling were not researched; I found no sources on their vectorization methods.

---

## Q5. ML lithography (neural ILT, GAN-OPC, LithoBench, CurvyILT, etc.): mask to polygon conversion and differentiability

### Takeaway
Academic ML and GPU ILT works almost entirely on pixel grids and ends with **binarization plus morphology**. Polygon/contour export and MRC on vector output are mostly left out of benchmarks. Open-source tools are only starting to add contour extraction and OASIS export (OpenLithoHub, 2026), so this is a visible gap.

### Cited Findings
- **CurvyILT** (Haoyu Yang and Haoxing Ren, NVIDIA; ISPD 2025; arXiv 2411.07311):
  - Optimizes a continuous pixel mask and finishes with a **threshold binarization plus morphological opening/closing**. No polygon extraction or spline fitting is described.
  - Curvilinear mask rules (minimum shape area, minimum shape distance) are handled by **differentiable morphological operators** during optimization.
  - Runtime is **2.11 s per clip** vs MultiILT's 3.45 s, with **0.6 GB vs 7.2 GB** peak memory on an RTX A6000, using ICCAD13 10x10 µm clips and LithoBench. Code: github.com/phdyang007/curvyILT.
  — [arXiv 2411.07311](https://arxiv.org/html/2411.07311v1)
- **LithoBench** (NeurIPS 2023 Datasets & Benchmarks, CUHK): **more than 120k tiles** from real designs or synthesized from ILT testcase topologies, with ground truth from a pixel-based ILT (multi-resolution, average pooling). — [LithoBench paper](https://proceedings.neurips.cc/paper_files/paper/2023/file/604b9fa9e1c16284e6517d923cf9ff20-Paper-Datasets_and_Benchmarks.pdf)
- **OpenLithoHub** (v0.1.0a2, May 2026):
  - Open-source OPC/ILT platform with GDSII/OASIS/DEF parsing, "Manhattan & curvilinear" contour extraction with OASIS/GDSII export, and full-chip core/halo tiling.
  - Metrics include EPE, PV band, **shot count** and MRC/DRC hard-fail gates.
  - Its stated limits: CPU reference hardware, GPU timing "pending Industrial Benchmark v2," no foundry or SEM calibration, and no comparison with commercial tools.
  — [GitHub OpenLithoHub](https://github.com/OpenLithoHub/OpenLithoHub)
- A patent on ML mask synthesis notes that "a variety of post-processing techniques can be used to convert the mask output by the ML model to a polygon representation" for lithography tools. — [USPTO 11762283](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/11762283)
- Pixel-level explicit ILT tends to produce isolated and slit structures that are hard to manufacture. Level-set methods give cleaner, more Manhattan-friendly masks (arXiv 2308.12299). Other approaches parametrize curvilinear masks with **B-splines and Delaunay triangulation** (arXiv 2504.11962, 2025) to avoid pixels. Recent work also includes MaskOpt, a large-scale mask optimization dataset (arXiv 2512.20655), and generative-RL ILT (arXiv 2602.19027). I only saw titles and snippets for these. — [arXiv 2308.12299](https://arxiv.org/pdf/2308.12299); [arXiv 2504.11962](https://arxiv.org/pdf/2504.11962); [arXiv 2512.20655](https://arxiv.org/pdf/2512.20655); [arXiv 2602.19027](https://arxiv.org/html/2602.19027v1)
- Photonics inverse design has the same problem. OptiClear (arXiv 2607.03632) does differentiable curvilinear design-rule legalization. — [arXiv 2607.03632](https://arxiv.org/pdf/2607.03632)

### Inferences
- Benchmarks score pixel masks (L2, PVB, EPE), but real deployment needs polygons plus MRC plus **shot count** or vertex count. A fast GPU bitmap-to-Manhattan-polygon kernel could serve as the missing **"deployability metric" stage**: polygon count, vertex count, VSB shot count, OASIS bytes and MRC violations per clip, computed inside the training loop. Today's flows cannot do this cheaply.
- Polygonization is not differentiable, but its outputs (vertex/jog counts) could act as non-differentiable rewards for RL-based ILT (cf. arXiv 2602.19027) or as straight-through regularizers. This is a creative research angle (inference).

### Gaps
- I could not confirm how LithoBench or the NeuralILT/GAN-OPC/DiffOPC codebases export polygons. They appear to stop at pixel masks, but I did not verify each repo.

---

## Q6. PCB and advanced packaging (AOI, image-to-Gerber, LDI / maskless lithography, panel-level packaging)

### Takeaway
PCB and packaging flows run vector (Gerber/ODB++/GDS) to raster in LDI/maskless writers through raster image processors (RIP). The reverse direction, image to vector, shows up in PCB reverse engineering and AOI. Both sides are almost entirely Manhattan or 45-degree geometry, but I found few current, well-sourced specifics.

### Cited Findings
- LDI writes the circuit with a UV laser (355 or 405 nm) driven from Gerber or ODB++ data, with no film phototool. The workstation load includes **raster image processing (RIP)**, real-time geometric scaling/rotation correction, and moving and storing "enormous data." Direct writing turns CAD/GDSII/Gerber into scan paths or projected (DMD) exposure patterns. — [AllPCB LDI blog](https://www.allpcb.com/ar-SA/blog/pcb-manufacturing/laser-direct-imaging.html); [USPTO 10149390 (multi-DMD maskless)](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10149390)
- Fan-out wafer-level packaging is driving use of **LDI with adaptive patterning**, where per-die placement errors are corrected by regenerating the pattern for each panel. — [IMAPS paper](https://imapsource.org/article/90153.pdf); [CiteDrive: Mask-less LDI & Adaptive Patterning for Fan-Out](https://citedrive.com/en/discovery/mask-less-laser-direct-imaging-amp-adaptive-patterning-for-fan-out-heterogeneous-integration)
- A patented PCB RE method converts **raster images of conductor layers into vectorized CAM files**: polygons defined by boundary points, rectangles/ovals, and stroked lines. — [USPTO 12450794](https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/12450794)
- Altium's CAM editor can reverse-engineer PCB layouts and netlists from Gerber/ODB++, but the output "may require validation and manual repair." — [Altium KB](https://www.altium.com/documentation/knowledge-base/altium-designer/convert-gerber-odb-fabrication-data-back-to-pcb)

### Inferences
- Adaptive patterning in panel-level packaging means imaging the panel, measuring die shifts and regenerating routing on the fly. That is a latency-critical loop where image-derived geometry must be fast (inference). PCB AOI and RE (scan to Gerber) are clean Manhattan bitmap-to-polygon use cases, but commercial demand and pain points are poorly documented.

### Gaps
- No 2024–2026 numbers found for LDI RIP throughput, AOI-to-CAD vectorization methods (Orbotech/KLA, Camtek), or panel-level adaptive patterning compute budgets (e.g., Deca Technologies).
- I found no academic GPU work on image-to-Gerber vectorization.

---

## Q7. Academic work on GPU-parallel polygon extraction from layouts or masks (ICCAD/DAC/SPIE 2022–2026)

### Takeaway
GPU acceleration of layout geometry exists for DRC (X-Check, ICCAD 2022), litho simulation, ILT and PLDC. I found **no dedicated 2022–2026 EDA/SPIE paper on GPU raster-to-polygon (contour/ring) extraction for masks or layouts**. The step is treated as solved or folded into proprietary tools, which leaves a publishable gap.

### Cited Findings
- **X-Check (ICCAD 2022, CUHK):** GPU-accelerated DRC using **parallel sweepline algorithms** on layout polygons. Sweepline is structurally close to the user's scanline approach. — [X-Check slides](https://www.cse.cuhk.edu.hk/~byu/papers/C149-ICCAD2022-GPU-DRC-slides.pdf)
- GPU geometric pattern matching for layouts with axis-parallel polygons dates from older work. — [BGU: Applying graphics hardware to achieve extremely fast geometric pattern matching](https://cris.bgu.ac.il/en/publications/applying-graphics-hardware-to-achieve-extremely-fast-geometric-pa/)
- In 2023–2026, GPU work in this area has gone into ILT (CurvyILT, cuLitho, the arXiv 2602.15036 "computational geometry" primitives), pixel MPC/PLDC in writers, and GPU-ML SEM contouring (Siemens/IMEC 2024), as cited above. None of it publishes the polygonization kernel itself. — [arXiv 2602.15036](https://arxiv.org/abs/2602.15036); [arXiv 2411.07311](https://arxiv.org/html/2411.07311v1)

### Inferences
- A paper titled roughly "Single-pass GPU rectilinear polygon extraction for pixel-domain mask synthesis and SEM-to-layout reconstruction" could be presented at ISPD/ICCAD/SPIE Photomask, given:
  - benchmarks against OpenCV `findContours`, GDAL polygonize, and KLayout region-from-bitmap;
  - outputs measured in OASIS bytes, VSB shot count and MRC violations;
  - test data from LithoBench/CurvyILT masks and the public RUB SEM/GDS Trojan dataset.
- It does not appear to have been done publicly. That is an inference from the search results, not a proven absence.

### Gaps
- My searches did not reach SPIE Photomask/BACUS 2025 proceedings or full DAC 2025/2026 programs, so I cannot rule out an existing paper. Proprietary MDP tools (Synopsys CATS, Siemens Calibre MDP, D2S TrueMask, NCS, ASML/Brion) very likely contain GPU or multithreaded contouring that is undisclosed.
