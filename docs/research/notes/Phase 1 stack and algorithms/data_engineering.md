# Data engineering for multi-rater medical segmentation (Phase 1: SFT + GRPO of a polygon-emitting small VLM)

Scope: splits, crop and windowing rules, storage format and tooling, licences, leakage hygiene and augmentation for RIGA, LIDC-IDRI and the next-track candidates. Research date 2026-10-18. Sources are primary (papers, PDFs, repo code, official pages) unless flagged. "UNVERIFIED" means the claim comes from memory or a secondary source and was not confirmed against a primary source in this pass.

Local context checked: `data/riga_true/` holds per-expert files named `<Subset>_<sub>_image<N>-<rater>.png` (for example `Magrabia_MagrabiFemale_image4-2.png`). `ROADMAP.md` Phase 1 currently plans "RIGA: train MESSIDOR plus BinRushed; test Magrabia" and "LIDC: split by patient 70 / 10 / 20".

---

## Q1. RIGA: splits, preprocessing and metrics used by published multi-rater work

### Takeaway
The de facto RIGA multi-rater protocol, set by MRNet (CVPR 2021) and followed by TAB, MSE-Nets, the multi-source-annotation paper and the self-calibrated OD/OC paper, is: train on BinRushed (195) + MESSIDOR (460) = 655 images, test on Magrabia (95). Inputs are resized to 256×256. The metric is "soft Dice" for disc and cup: hard Dice averaged over binarization thresholds {0.1, 0.3, 0.5, 0.7, 0.9}, reported against the mean-vote soft GT and/or averaged over the 6 raters. This is the opposite of the "train BinRushed+Magrabia, test MESSIDOR" option raised in the question. The ROADMAP plan (train MESSIDOR+BinRushed, test Magrabia) already matches the literature.

### Cited Findings
- MRNet: "we select 195 samples from BinRushed and 460 samples from MESSIDOR as the training set… The Magrabia set with 95 samples is selected as" the test set. RIGA = 460 MESSIDOR + 195 BinRushed + 95 Magrabia, annotated by six glaucoma experts. — [MRNet CVPR 2021 paper (PDF)](https://openaccess.thecvf.com/content/CVPR2021/papers/Ji_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-Rater_Agreement_Modeling_CVPR_2021_paper.pdf)
- MRNet: "All training and test images are uniformly resized to the dimension of 256×256 pixels"; DeepLab-V3+ / VGG-16; batch 8; 60 epochs; LR 1e-4. — [MRNet paper](https://openaccess.thecvf.com/content/CVPR2021/papers/Ji_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-Rater_Agreement_Modeling_CVPR_2021_paper.pdf)
- MRNet metric: soft Dice (Ds) and soft IoU over thresholds (0.1, 0.3, 0.5, 0.7, 0.9). Both the predicted probability map and the soft GT are binarized at each threshold, and the scores are averaged. — [MRNet paper](https://openaccess.thecvf.com/content/CVPR2021/papers/Ji_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-Rater_Agreement_Modeling_CVPR_2021_paper.pdf)
- MRNet code: `config.SCALE_SIZE = 256`. The loader reads `Glaucoma_multirater_<split>.csv` with columns `imgName, maskName, label, center, xmin, ymin, xmax, ymax, width, height, discFlag, rater`. It keeps rows with `rater == 0` ("FinalLabel") and loads Rater1–Rater6 masks by replacing `FinalLabel` in the path. One grey-level mask per rater encodes disc (≠0) and cup (==255, after /255 → 1). Augmentations are flip and rotation; random crop is commented out. — [MRNet repo: dataset_loader_cvpr.py](https://github.com/jiwei0921/MRNet/blob/main/MRNet_Code/utils/dataset_loader_cvpr.py), [config.py](https://github.com/jiwei0921/MRNet/blob/main/MRNet_Code/config.py)
- MRNet's preprocessed RIGA release is distributed via Baidu Pan (code 1627) and Google Drive links in the README. The repo is MIT, but the README gives no dataset licence. — [MRNet GitHub](https://github.com/jiwei0921/MRNet)
- TAB (MICCAI 2023): "We use 655 samples from BinRushed and MESSIDOR for training and 95 samples from Magrabia for test [26,9,14]". Input resized to 256×256, batch 8. Metrics: Soft Dice at (0.1, 0.3, 0.5, 0.7, 0.9), reported as "Average" (mean of per-annotator soft Dice) and "Mean Voting" (vs. mean-vote GT), for disc and cup. Baselines in its tables: MR-Net, CM-Net, PADL. — [TAB arXiv 2306.01340](https://arxiv.org/abs/2306.01340)
- TAB's README points to the MRNet repo for RIGA download and details, so its data comes from MRNet's preprocessed release. — [TAB GitHub](https://github.com/Merrical/TAB)
- MSE-Nets uses the same 195 BinRushed + 460 MESSIDOR train / 95 Magrabia test split, "following the RIGA benchmark", with inputs resized to 256×256. — [MSE-Nets arXiv 2311.10380](https://arxiv.org/abs/2311.10380)
- "Learning Robust Medical Image Segmentation from Multi-source Annotations": "All images are resized to 256×256"; 655 BinRushed+MESSIDOR train, 95 Magrabia test. — [arXiv 2304.00466](https://arxiv.org/abs/2304.00466)
- The self-calibrated OD/OC multi-rater paper (Wu et al., MICCAI 2022) also resizes all images to 256×256 and reports MV and STAPLE variants. — [arXiv 2206.05092](https://arxiv.org/abs/2206.05092)
- The OD/OC noise-aware pseudo-label paper reuses the 195 BinRushed + 460 MESSIDOR / 95 Magrabia split and a 256×256 input. — [arXiv 2311.18496](https://arxiv.org/abs/2311.18496)
- D-Persona (CVPR 2024) does not use RIGA. It evaluates on NPC-170 (in-house, now released via the MMIS-2024 challenge) and LIDC-IDRI. — [D-Persona arXiv 2403.13417](https://arxiv.org/abs/2403.13417), [D-Persona GitHub](https://github.com/ycwu1997/D-Persona)
- ProSeg (arXiv 2512.00748, v2 May 2026), the successor to D-Persona, also uses only NPC + LIDC-IDRI. Metrics: GED, Dsoft, Dmax, Dmatch. — [arXiv 2512.00748](https://arxiv.org/abs/2512.00748)

### Inferences
- To be comparable on RIGA, report disc and cup soft Dice at thresholds {0.1, 0.3, 0.5, 0.7, 0.9} on the 95 Magrabia images. Use (a) mean-vote GT and (b) the average over the 6 raters, at 256×256 resolution. A polygon-emitting VLM gives a hard mask per sample, so a "probability map" must come from K sampled polygons (pixel-wise frequency) before soft Dice is computed. K must be stated, for example K=6 to mirror the rater count.
- MRNet's released crops appear to be disc-region crops (the CSV carries bbox and center columns). If our RIGA images are full fundus images, a 256×256 resize of the full image is not the same input as the papers' 256×256 disc crop. We should either use MRNet's released crops or reproduce a disc-centred crop. The exact crop size of MRNet's release is UNVERIFIED (see Gaps).
- Validation set: none of the papers documents one. Hold out a fixed slice of the 655 training images, for example stratified by subset (BinRushed1–4, MESSIDOR), for SFT/GRPO checkpoint selection. Never tune on Magrabia.
- Magrabia is split into Female/Male (47 + 48 = 95 per the literature; the brief says 47 each). Reconcile the count: MRNet, TAB and MSE-Nets all state 95.

### Gaps
- MRNet's preprocessed release (exact crop size before resize, crop method and the definition of "FinalLabel"/rater 0) was not inspected. The Google Drive/Baidu archives were not downloaded. A search result attributes OD cropping to "a pre-trained disc-aware ensemble network", but I could not tie that to a specific RIGA multi-rater paper. UNVERIFIED.
- No RIGA-specific per-rater "personalization" protocol (D-Persona style Dmatch on RIGA) was found in the sources checked. TAB reports per-annotator columns.
- Magrabia count discrepancy: the brief says 47 + 47 = 94; the papers say 95. This needs a local count.

---

## Q2. LIDC-IDRI: benchmark splits, preprocessing, small nodules, non-marking readers, pylidc vs own parser, XML quirks, licence

### Takeaway
Three incompatible "standard" LIDC protocols exist. Pick one per reported metric and say which.
- **(A) Prob-U-Net / HPU-Net (DeepMind)**: 0.5 mm in-plane resample, 180×180 lesion-centred crops, random 128×128 training crops, patient split. Public PNG release on GCS, CC BY 3.0.
- **(B) Knegt pickle (PyTorch Prob-U-Net)**: 128×128 patches, 15,096 total, 60/20/20.
- **(C) MedicalMatting / D-Persona / ProSeg**: 1,609 2D slices from 214 subjects, 128×128 centre crop, 4-fold patient-level CV. The 4 masks are re-ordered by area into "virtual annotators".

Only ≥3 mm nodules have contours. Readers who did not mark a lesion give empty masks; up to 3 of 4 can be empty. pylidc (v0.2.3, last release Dec 2021) does annotation clustering (connected components at `tol` = pixel spacing, shrinking by 0.9 until each cluster has ≤4 annotations) and consensus (`clevel`). That makes it the right reference to validate our own XML parser. Licence: CC BY 3.0 with a mandatory citation/DOI.

### Cited Findings
**Protocol A (Kohl et al. 2018, Prob-U-Net)**
- 1,018 CT scans from 1,010 patients, 4 of 12 radiologists per scan. It uses the masks from the second (unblinded) reading. Patient split 722 / 144 / 144. — [Prob-U-Net arXiv 1806.05034](https://arxiv.org/abs/1806.05034)
- Resampled to 0.5 mm × 0.5 mm in-plane (original 0.461–0.977 mm, mean 0.688 mm). 180×180 crops are centred at lesion positions "where at least one of the experts segmented a lesion". That gives 8,882 / 1,996 / 1,992 images. "Because the experts can disagree whether the lesion is abnormal tissue, up to 3 masks per image can be empty." — [Prob-U-Net paper](https://arxiv.org/abs/1806.05034)
- Masks from different graders count as the same lesion if their tightest bounding boxes overlap. Training augmentation: elastic deformation, rotation, shear, scale, and a randomly translated crop from 180×180 down to 128×128. — [Prob-U-Net paper, appendix](https://arxiv.org/abs/1806.05034)
- HPU-Net public LIDC release: 180×180 PNG crops at 0.5 mm. Splits are 8,843 images / 530 patients (train), 1,993 / 111 (val), 1,980 / 103 (test). Only polygon-annotated (≥3 mm) abnormalities are kept. Bounding-box-overlap matching across graders. Hosted at `gs://hpunet-data/lidc_crops/`. CC BY 3.0. — [deepmind-research/hierarchical_probabilistic_unet](https://github.com/deepmind/deepmind-research/tree/master/hierarchical_probabilistic_unet)
- Conflict: the patient counts in the HPU-Net release (530/111/103 = 744 patients with retained lesions) differ from the Prob-U-Net paper's 722/144/144 split (of all 1,010 patients). The image counts are also slightly different (8,843 vs 8,882 train). Treat them as two different releases. — [Prob-U-Net paper](https://arxiv.org/abs/1806.05034) vs [HPU-Net README](https://github.com/deepmind/deepmind-research/tree/master/hierarchical_probabilistic_unet)
- PHiSeg: "Similar to [4] we extracted square 2D patches of size 128 × 128 pixels such that each patch was centred on a lesion" (1,018 scans, 4 radiologists). — [PHiSeg arXiv 1906.04045](https://arxiv.org/abs/1906.04045)

**Protocol B (Knegt pickle)**
- A widely reused PyTorch Prob-U-Net ships a preprocessed LIDC of 128×128 lesion patches. A later paper reports 15,096 patches split 60:20:20 (9,058 / 3,019 / 3,019), intensities scaled to [0,1] — [Generalized Prob-U-Net arXiv 2207.12872 (via search snippet)](https://arxiv.org/abs/2207.12872); repo [stefanknegt/Probabilistic-Unet-Pytorch](https://github.com/stefanknegt/Probabilistic-Unet-Pytorch). Whether this split is patient-level is UNVERIFIED. A random patch-level split leaks adjacent slices.

**Protocol C (MedicalMatting → D-Persona → ProSeg)**
- "The LIDC-IDRI dataset contains 1609 2D thoracic CT scans belonging to 214 subjects… four binary masks… we manually rank the areas of four provided annotations… the first expert is… conservative… the fourth… largest annotation regions. Following [42], a four-fold cross-validation setting at the patient level is used." [42] = Wang et al., "Medical Matting", Computers in Biology and Medicine 158:106714, 2023. — [D-Persona arXiv 2403.13417](https://arxiv.org/abs/2403.13417)
- Inputs are "center-cropped into a fixed size of 128 × 128"; LIDC preprocessing follows [42]; flip and rotation augmentation; GED and soft Dice at thresholds {0.1, 0.3, 0.5, 0.7, 0.9}. Personalization uses Dmax and Dmatch. — [D-Persona](https://arxiv.org/abs/2403.13417)
- ProSeg repeats this: 1,609 images / 214 subjects, area-ranked virtual annotators, 4-fold CV "based on the patient ID to reduce bias… from the similarity between adjacent lung nodule slices", 128×128 centre crop. — [ProSeg arXiv 2512.00748](https://arxiv.org/abs/2512.00748)

**Small nodules, non-nodules, reading protocol, size, licence (TCIA)**
- Two-phase reading: a blinded read, then an unblinded read with access to anonymized peer marks. Each case was read by 4 of 12 radiologists, who marked "nodule ≥3 mm" (full contours), "nodule <3 mm" (centroid only) and "non-nodule ≥3 mm" (point only). — [TCIA LIDC-IDRI collection](https://www.cancerimagingarchive.net/collection/lidc-idri/)
- Size: 1,010 patients, 1,308 studies/series, 244,527 images, 133.16 GB. — [TCIA](https://www.cancerimagingarchive.net/collection/lidc-idri/)
- Known issues on TCIA: about 100 cases have inconsistent spiculation/lobulation ratings. XML files 044.xml and 191.xml have incomplete characteristic ratings. Eight patients have two timepoint scans. — [TCIA](https://www.cancerimagingarchive.net/collection/lidc-idri/)
- Licence CC BY 3.0. The required data citation is "Data From LIDC-IDRI", DOI 10.7937/K9/TCIA.2015.LO9QL9SX. An acknowledgement of the NCI and the Foundation for the NIH is also required. TCIA recommends pylidc and a standardized DICOM (SEG/SR) representation of the annotations. — [TCIA](https://www.cancerimagingarchive.net/collection/lidc-idri/)

**pylidc**
- `Scan.cluster_annotations(metric='min', tol=None, factor=0.9, min_tol=0.1, return_distance_matrix=False, verbose=True)`. Annotations are adjacent when the minimum boundary-point distance is ≤ `tol` (default `tol = scan.pixel_spacing`, in mm). Clusters are the connected components. If any cluster has >4 annotations, `tol *= factor` and it repeats. Below `min_tol` it gives up with a warning ("cannot be automatically reduced to have groups with each group having Annotations<=4"). — [pylidc source Scan.py](https://github.com/notmatthancock/pylidc/blob/master/pylidc/Scan.py)
- `pylidc.utils.consensus(anns, clevel=0.5, pad=[(20,20),(20,20),(0,0)])` returns the consensus mask, its bbox, and the per-annotation boolean masks "placed in a common frame of reference". — [pylidc consensus tutorial](https://pylidc.github.io/tuts/consensus.html)
- pylidc `Annotation` covers only "nodule >= 3mm" annotations. Contours carry an `inclusion` flag that adds or subtracts area. `contour_slice_indices` / `contour_slice_zvals` map contours to scan slices by z. — [pylidc Annotation docs](https://pylidc.github.io/annotation.html)
- Maintenance: last PyPI release 0.2.3 (2021-12-21); last repo push 2024-07-24; 11 open issues; not archived. — [PyPI pylidc](https://pypi.org/project/pylidc/), [GitHub notmatthancock/pylidc](https://github.com/notmatthancock/pylidc)

### Inferences
- Readers who did not mark: keep them as explicit empty masks, not missing. The rater-set semantics ("3 of 4 saw nothing") are the ambiguity we want the policy to reproduce. For GRPO, an empty polygon set needs a defined reward, for example IoU(∅,∅)=1, and a defined output token such as `<none>`. Store `rater_present` (read this case) separately from `rater_marked` (drew a contour) per slice.
- <3 mm nodules and non-nodules: drop them as segmentation targets, since they have no contours. Optionally keep the centroids as metadata or hard negatives. This matches HPU-Net ("only polygons ≥3 mm retained") and pylidc's scope.
- Clustering: for comparability with Prob-U-Net/HPU-Net (bbox-overlap per 2D crop) and with pylidc (3D boundary distance ≤ pixel spacing), run pylidc's `cluster_annotations` on the same scans as our parser and report agreement. pylidc's adaptive-tol rule is a sensible spec to mirror.
- Windowing: none of the benchmark papers checked documents an HU window. Prob-U-Net releases PNGs that are already windowed. Keep raw HU (int16) in storage and apply a window at render time, for example lung window W1500/L-600 (UNVERIFIED as "standard"; a common clinical convention). That lets the VLM image policy change without re-extraction.
- Crop: generate lesion-centred 180×180 crops at 0.5 mm (protocol A geometry) so one store serves A-style random 128 crops and C-style 128 centre crops. Record the crop origin in mm so polygons can be mapped back to scan coordinates for 3D work.
- Split: patient-level, with the patient ID taken from the DICOM PatientID (`LIDC-IDRI-NNNN`), not the series. The 8 patients with two timepoints must stay in one fold.
- Our ROADMAP 70/10/20 patient split is defensible, but it is not identical to any published one. For literature comparability, also run D-Persona's 4-fold patient CV on the 1,609-slice MedicalMatting subset if we can obtain its file list.

### Gaps
- The MedicalMatting preprocessed LIDC (the 1,609 slices / 214 subjects list, how slices were chosen, and HU windowing) was not inspected. The D-Persona README says it uses "its pre-processed version as MedicalMatting" but gives no link in the text retrieved.
- LIDC XML quirks not confirmed from primary sources in this pass (all UNVERIFIED, from memory):
  - (i) `inclusion=FALSE` contours mark holes/exclusions.
  - (ii) The LIDC drawing instruction placed outlines just outside the nodule, so outline pixels may not belong to the nodule. This affects mask-from-polygon rasterization by about 1 px.
  - (iii) Some XMLs have >4 `readingSession`s or duplicate reader IDs.
  - (iv) z lookup via `imageZposition` vs `imageSOP_UID` can disagree.
  - (v) The XML holds the unblinded (post-review) marks, not the blinded ones.
  These should be checked against Armato et al. 2011 (Med Phys 38:915) and the LIDC XML documentation on TCIA.
- Whether the Knegt pickle split is patient-level is unconfirmed.
- The exact HU windows used by the HPU-Net PNG release are unknown.

---

## Q3. Next-track candidate datasets: licence, access, raters, size

### Takeaway
Realistic next steps, in order:
1. **IMA++** (ISIC, 2,394 multi-annotator images, official splits, masks on Zenodo, images from ISIC with CC-0/CC-BY-NC).
2. **NuCLS** (CC0, multi-rater subsets, on HF).
3. **QUBIQ** 2D tasks (prostate 6 raters, brain-growth 7, kidney 3, brain-tumour 3; small).
4. **DRIVE** (only 2 observers on the 20 test images, so weak as multi-rater).

MSSEG-2016 (7 raters, 3D MRI, Shanoir registration) and Gleason2019 (6 pathologists, 244 TMAs, DUA, multi-class) are heavier and lower priority.

### Cited Findings
- **IMA++**: 17,684 masks over 14,967 dermoscopic images. 2,394 images have 2–5 masks: 2,130×2, 209×3, 51×4, 4×5, for 5,111 masks. There are 15–16 annotators (the paper says 15 in one place and 16 in another). Metadata covers the tool (manual polygon / semi-auto flood-fill / auto + review) and skill (expert/novice). Majority-vote and STAPLE consensus are included. Recommended multi-annotator splits are 1,675 / 240 / 479, stratified by mask count and agreement. Masks are PNG on Zenodo (DOI 10.5281/zenodo.14201692). Images must be fetched from ISIC Archive collection 482 via the ISIC API; their licences are CC-0 or CC-BY-NC. — [IMA++ arXiv 2512.21472](https://arxiv.org/abs/2512.21472), [GitHub sfu-mial/IMAplusplus](https://github.com/sfu-mial/IMAplusplus)
- The IMA++ README has a commented-out "CC BY-NC 4.0" badge. The actual Zenodo licence of the mask package was not confirmed (UNVERIFIED). — [IMAplusplus README](https://github.com/sfu-mial/IMAplusplus)
- **QUBIQ** (MICCAI 2020/2021): 2D prostate MRI (55 cases), brain-growth MRI (39), brain-tumour multimodal MRI (32), kidney CT (24); 3D pancreas CT (38), pancreatic lesion CT (21). Metric: Q-Dice, which binarizes GT and prediction at 0.1…0.9 and averages Dice. — [QUBIQ challenge paper arXiv 2405.18435](https://arxiv.org/abs/2405.18435), [qubiq21.grand-challenge.org](https://qubiq21.grand-challenge.org/)
- QUBIQ annotation counts: brain growth 7, tumour 3, kidney 3, prostate 6 — [secondary survey via search snippet, arXiv 2603.27460](https://arxiv.org/abs/2603.27460). UNVERIFIED against the challenge site. The "CC BY-NC-ND 4.0" that WebFetch extracted from the QUBIQ arXiv page may be the paper's licence, not the data's. Data licence UNVERIFIED. Download is via Grand Challenge registration (UNVERIFIED in detail).
- **MSSEG-2016**: 7 expert raters plus a consensus. Training set of 15 patients from 3 sites. Raw and preprocessed packages. Access is "open" via the Shanoir platform after account registration (select study MSSEG-2016). — [FLI-IAM MSSEG page](https://portal.fli-iam.irisa.fr/english-msseg/), [database description](https://portal.fli-iam.irisa.fr/?p=870)
- **Gleason2019**: 244 prostate TMA images (H&E), pixel-level annotations by six pathologists. Registration plus a Data Use and Confidentiality Agreement are required. — [Gleason2019 Grand Challenge](https://gleason2019.grand-challenge.org/), [arXiv 2005.04886](https://arxiv.org/abs/2005.04886)
- **NuCLS**: >220,000 labelled nuclei from TCGA breast cancer; single- and multi-rater annotations by pathologists, residents and students; CC0 1.0 with a citation request (Amgad et al., arXiv 2102.09099). A community copy is on the HF Hub (`minhanhto09/NuCLS_dataset`). — [NuCLS licence page](https://nucls.grand-challenge.org/License/), [NuCLS paper (PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9112766/), [HF copy](https://huggingface.co/datasets/minhanhto09/NuCLS_dataset)
- **DRIVE**: 40 fundus images, 565×584. 20 train with one manual mask; 20 test with two observers (first = GT, second = human baseline). Download from Grand Challenge. — [drive.grand-challenge.org](https://drive.grand-challenge.org/)

### Inferences
| Dataset | Raters | Size | Access | Licence | Fit for next track |
|---|---|---|---|---|---|
| IMA++ | 2–5 per image (15–16 annotators total) | 2,394 multi-rater images | Zenodo + ISIC API, no registration | masks: UNVERIFIED (likely CC BY-NC); images CC-0/CC-BY-NC | **High**: 2D RGB, blob-like single polygons, official splits, rater metadata |
| NuCLS | multi-rater subsets | 220k nuclei | direct/HF | CC0 | **High** for many-small-polygon outputs; different task shape (instances) |
| QUBIQ 2D | 3–7 | 24–55 cases per task | Grand Challenge registration | UNVERIFIED | **Medium**: tiny, but a classic multi-rater benchmark with Q-Dice |
| DRIVE | 2 (test only) | 40 | Grand Challenge | UNVERIFIED (research use) | **Low**: vessels are not polygon-friendly; only 2 raters |
| MSSEG-2016 | 7 | 15 train patients, 3D | Shanoir account | DUA (UNVERIFIED terms) | **Low/Medium**: 3D, many small lesions |
| Gleason2019 | 6 | 244 TMAs | registration + DUA | DUA | **Low**: multi-class, giant images |

- Non-commercial clauses (ISIC CC-BY-NC images, RIGA, likely Gleason) mean derived HF datasets must be private or gated, and model release terms should be checked.

### Gaps
- Exact data licences for QUBIQ, DRIVE, MSSEG-2016 and the IMA++ mask package, and the size in GB of each, were not confirmed.
- NuCLS multi-rater subset sizes (the number of FOVs with ≥2 raters) were not retrieved.

---

## Q4. Storage, versioning and trainer input formats (2026)

### Takeaway
- **Master store**: a Parquet-backed HF `datasets` dataset (private Hub repo or local). It has one row per (image crop × task). Images are embedded as bytes (PNG for RGB; 16-bit PNG or raw int16 + window params for CT). Each rater's geometry is a list of polygons in a JSON/struct column, with optional COCO-RLE masks. Polygons are our model's output language and the most compact form.
- **Adapter views**: thin exports to TRL (conversational `prompt`/`completion` or `messages` + `images`), verl (Parquet with `data_source`, `prompt`, `images`, `ability`, `reward_model.ground_truth`, `extra_info`) and LLaMA-Factory (ShareGPT JSON with `<image>` tokens and image paths).
- **Versioning**: content hashes (SHA-256 per file plus a manifest) committed to git, with or without DVC. DVC is maintained but now under lakeFS stewardship; that is fine for small datasets.
- Lance is the stronger choice only if we need random access over millions of rows or vector indexes. Not needed at RIGA/LIDC scale.

### Cited Findings
- TRL vision datasets: use the conversational format. The row must contain `images` (list of PIL images) or `image` (single). Message `content` is a list of `{"type":"image"}` / `{"type":"text","text":…}` dicts. Mixing text-only and VLM rows requires transformers ≥4.57.0. GRPOTrainer expects prompt-only; SFTTrainer expects language-modeling (`messages`) or prompt-completion. — [TRL dataset formats](https://huggingface.co/docs/trl/main/en/dataset_formats)
- TRL: `datasets` ≥4.7.0 has a `Json()` feature type for arbitrary JSON columns. Older versions should store JSON as strings. — [TRL dataset formats](https://huggingface.co/docs/trl/main/en/dataset_formats)
- verl Parquet schema: `data_source` (selects the reward function), `prompt` (chat-template messages), `ability`, `reward_model` (with `ground_truth`), `extra_info`. — [verl data preparation docs](https://verl.readthedocs.io/en/latest/preparation/prepare_data.html)
- verl's multimodal example (geo3k) writes `"images": images` (taken directly from the HF dataset column), `"prompt": [{"role":"user","content": prompt}]`, `"reward_model": {"style":"rule","ground_truth": answer}`. — [verl examples/data_preprocess/geo3k.py](https://github.com/volcengine/verl/blob/main/examples/data_preprocess/geo3k.py)
- LLaMA-Factory multimodal ShareGPT: `{"conversations":[{"from":"human","value":"<image>…"},{"from":"gpt","value":"…"}],"images":["path"]}`. "The number of images should be identical to the `<image>` tokens in the conversations." Registered in `dataset_info.json` with `"formatting":"sharegpt"` and a `columns` mapping. — [LLaMA-Factory data/README.md](https://github.com/hiyouga/LLaMA-Factory/blob/main/data/README.md)
- Lance: Apache-2.0, "Open Lakehouse Format for Multimodal AI… 100x faster random access, vector index, and data versioning". Actively developed (push 2026-10-19, about 7.2k stars). — [lance-format/lance](https://github.com/lance-format/lance)
- Lance is natively supported on the HF Hub (`hf://` paths, scan/filter/search without a full download). The blog dates this to around Feb 2026 (date from a secondary snippet). — [LanceDB blog: Lance × Hugging Face](https://lancedb.com/blog/lance-x-huggingface-a-new-era-of-sharing-multimodal-data)
- DVC: lakeFS acquired the DVC open-source project from Iterative.ai (announced 2025-11-18). DVC continues "as an independent open source tool focused on data scientists working with smaller datasets". — [lakeFS blog](https://lakefs.io/blog/lakefs-acquires-dvc/)

### Inferences
- Proposed row schema (one row per crop):

```
id (str, stable)          dataset (riga|lidc|...)     split (train|val|test|fold_k)
group_id (patient/subject/source image)              subset (MESSIDOR|BinRushed3|Magrabia...)
image (Image, PNG bytes)  image_sha256                  pixel_spacing_mm, crop_origin, crop_size, window (for CT)
raters: list<struct{rater_id, present: bool, polygons: list<list<[x,y]>>, holes, rle?}>
consensus: struct{mv_polygons, staple_rle?}          meta: Json (nodule_id, slice_z, cluster_id, source_xml...)
```

  Training views are generated from this deterministically, for example the SFT target = one sampled rater's polygon serialized in the model's coordinate grammar. Rater identity stays in metadata, so personalization experiments (D-Persona-style) remain possible.
- Polygons vs RLE vs PNG stacks: for RIGA (6 raters × 2 structures) and LIDC (≤4 raters, 1 structure per crop), polygons are a few hundred bytes per rater. Store them as the source of truth because they are what we parsed from XML/outlines. Rasterize on the fly for metrics, with a single shared rasterizer (the same one used in the reward) to avoid off-by-one discrepancies. Keep COCO RLE (pycocotools) only as a cached derivative for fast IoU.
- Licence-driven storage: RIGA prohibits redistribution (see Q5), so any Hub copy must be a private repo. LIDC derivatives can be public under CC BY 3.0 with attribution.
- Scale check: RIGA ≈ 750 images × 6 raters; LIDC ≈ 10–15k crops. Both fit in a few-GB Parquet. Webdataset/Zarr bring no benefit at this size. Zarr is only worth it for 3D LIDC volumes (Phase 3).

### Gaps
- HF `datasets` `Image()` Parquet embedding and the exact `datasets` version in 2026 were not fetched in this pass. The `Json()` type (≥4.7.0) is confirmed only via the TRL docs.
- Whether verl requires `<image>` placeholders in the prompt text for Qwen-VL-style models was not confirmed: the geo3k script passes `problem` text unchanged, and the source dataset may already embed `<image>`.

---

## Q5. Leakage and split hygiene

### Takeaway
Group every split by the true source unit: LIDC PatientID (including two-timepoint patients), RIGA source image across raters, IMA++ patient/lesion IDs. Run perceptual-hash near-duplicate detection across RIGA subsets and against MESSIDOR-derived public sets before freezing splits. Use group-aware stratified splitting for validation folds.

### Cited Findings
- ProSeg and D-Persona split LIDC by patient ID "to reduce bias in results arising from the similarity between adjacent lung nodule slices". — [ProSeg arXiv 2512.00748](https://arxiv.org/abs/2512.00748)
- Eight LIDC patients have two timepoint scans. — [TCIA](https://www.cancerimagingarchive.net/collection/lidc-idri/)
- RIGA usage terms: free for research and educational purposes; "copy, redistribution, and any unauthorized commercial use are prohibited". Hosted on U-Michigan Deep Blue (12.9 GB). The MESSIDOR file holds 460 originals + 460 × 6 marked images = 3,220 images. — [Deep Blue RIGA record (via search snippet; direct fetch returned 403)](https://deepblue.lib.umich.edu/data/concern/data_sets/3b591905z)
- IMA++ provides patient/lesion identifiers in its metadata and official splits stratified by segmentation count and agreement level. — [IMA++ arXiv 2512.21472](https://arxiv.org/abs/2512.21472)
- imagededup (idealo): Apache-2.0, last push 2025-08-15, not archived. It offers hashing (PHash/DHash/WHash/AHash) plus CNN-based near-duplicate finding (method list from memory; repo status verified). — [idealo/imagededup](https://github.com/idealo/imagededup)

### Inferences
- RIGA leakage risks to test:
  - (a) The same eye or patient across BinRushed1–4 sub-folders or across Magrabia Female/Male. Literature splits by subset, so cross-subset duplicates would leak only if they cross train/test, which is unlikely: Magrabia is a different hospital.
  - (b) RIGA-MESSIDOR images are a subset of public MESSIDOR. Any VLM pretraining or other fundus data (REFUGE, Drishti, ORIGA, MESSIDOR-2) may contain them. The base VLM may have seen MESSIDOR.
  - (c) Rater files of one image must never be split apart: group by `<subset>_<image>` stripping the `-<rater>` suffix.
- Use `sklearn.model_selection.StratifiedGroupKFold` (groups = patient; strata = lesion-size bin × rater-agreement bin × number of raters marking) for LIDC folds and RIGA validation. This tool is from memory: it exists in scikit-learn ≥1.0 (UNVERIFIED in this pass).
- Stratify LIDC test folds so that "reader-disagreement slices" (0 < #marking < 4, or topology disagreement) appear in each split, as the ROADMAP requires.

### Gaps
- MESSIDOR's own licence terms (the RIGA-MESSIDOR images derive from it) were not checked.
- No published near-duplicate audit of RIGA was found.

---

## Q6. Augmentation for medical crops with polygons/masks

### Takeaway
Use torchvision transforms v2 with `tv_tensors.KeyPoints`. It is BSD-licensed, and since v0.23 it natively represents polygons as `[N_polygons, N_points, 2]` and transforms them jointly with images and masks. The original albumentations (MIT) was archived in June 2025. Its successor AlbumentationsX is AGPL-3.0 or commercial, and it has keypoints but no polygon type, so vertex order and closure would be handled by hand.

### Cited Findings
- torchvision `tv_tensors.KeyPoints` (beta since torchvision 0.23) can represent polygonal chains and polygons of shape `[N_polygons, N_points, 2]`. It works with all transforms in `torchvision.transforms.v2`. — [torchvision KeyPoints docs (stable)](https://docs.pytorch.org/vision/stable/generated/torchvision.tv_tensors.KeyPoints.html)
- `albumentations-team/albumentations`: archived, MIT, last push 2025-06-25. `AlbumentationsX`: active (push 2026-10-19), "AGPL-3.0-only or commercial licensing". — [albumentations repo](https://github.com/albumentations-team/albumentations), [AlbumentationsX repo](https://github.com/albumentations-team/AlbumentationsX)
- Albumentations keypoint formats are xy, yx, xya, xys, xyas, xysa, xyz. `KeypointParams` has `remove_invisible` (default True), `label_fields` and `label_mapping`. The docs do not describe native polygon support. — [Albumentations keypoint docs](https://albumentations.ai/docs/3-basic-usage/keypoint-augmentations/)
- Literature augmentations: MRNet uses random flip and rotation (random crop disabled). D-Persona uses flip and rotation on LIDC. Prob-U-Net uses elastic, rotation, shear, scale and a random 128 crop from 180. — [MRNet loader](https://github.com/jiwei0921/MRNet/blob/main/MRNet_Code/utils/dataset_loader_cvpr.py), [D-Persona](https://arxiv.org/abs/2403.13417), [Prob-U-Net](https://arxiv.org/abs/1806.05034)

### Inferences
- Polygon-safe policy: use only geometric transforms that keep polygons valid (flip, 90° and small rotations, scale, translation/crop) via `KeyPoints`. After each transform, re-validate: clip to the canvas, drop degenerate rings, re-orient CCW (a flip reverses orientation) and re-establish a canonical start vertex, so that SFT targets stay deterministic. Elastic deformation is unsafe for polygons unless vertices are densified first. Prefer applying it to masks, then re-polygonizing with our certified polygonizer.
- Photometric transforms only on the image: CLAHE/gamma for fundus; HU-window jitter for CT (for example ±100 HU on level, ±200 on width — a heuristic, UNVERIFIED).
- To avoid the AGPL, do not adopt AlbumentationsX. A pinned archived albumentations (MIT, ≤2.0.x) is acceptable but unmaintained.

### Gaps
- Kornia's polygon/keypoint augmentation support (`kornia.augmentation` with keypoints) was not checked in this pass.
- Edge-case behaviour of torchvision `KeyPoints` with `RandomRotation(expand=False)` clipping (whether vertices outside the canvas are clamped or kept) was not verified.

---

### Candidate hypotheses to test locally
1. **RIGA count and grouping (already run, 2026-10-18; it found a problem)**: unique stems in `data/riga_true/` are MESSIDOR 460, BinRushed2 46, BinRushed3 47, BinRushed4 51, MagrabiFemale 47, MagrabiaMale 47. That is 698 images, against 750 in the literature.
   - **BinRushed1 is missing entirely.** `data/raw/riga/BinRushed/BinRushed1` has 301 files and mixes `.jpg` with some `.tif` rater files, for example `image1-3.tif`. Also, `RIGA_preprocessed/BinRushed/BinRushed1-Corrected` exists.
   - Magrabia has 94 images, not the 95 in the papers.
   - 12 images have fewer than 6 rater files: ten have 5, one has 4 and one has 2 (for example `BinRushed_BinRushed3_image15`).
   - The literature training set is 655 and the test set 95. We are short about 50 BinRushed images and 1 Magrabia image, so our numbers are not yet comparable. Fix the extractor: handle `.tif`, the BinRushed1 naming, and the `-Corrected` folder. Then re-count before freezing splits.
2. **RIGA near-duplicates**: pHash (imagededup, Hamming ≤ 6–10) across all 750 originals. Expect 0 cross-subset duplicates, especially none between Magrabia and train. Also hash against any other fundus sets in our pipeline.
3. **Our extraction vs MRNet's release**: on about 20 images, compare our per-rater disc/cup masks with MRNet's Rater1–6 PNGs (Dice > 0.97 expected). Find out whether MRNet images are disc crops and at what size. Adopt their crop if so.
4. **Literature-metric reproduction**: a trivial baseline (for example, the mean of the 6 rater masks from train as a prior, or a U-Net) evaluated with our soft-Dice implementation on Magrabia at 256×256. It should land in the published range: disc about 95–97%, cup about 80–88% (TAB Table 1 ranges). This validates the metric code before VLM runs.
5. **LIDC clustering parity**: on about 50 scans, run `pylidc.Scan.cluster_annotations()` and our XML-parser clustering. Measure (a) the same number of nodules per scan, (b) the same reader membership per cluster, (c) per-reader mask Dice > 0.95 against `pylidc.utils.consensus(..., clevel=0.5)` outputs and `Annotation.boolean_mask()`. Log scans where pylidc hits `min_tol`.
6. **LIDC rasterization convention**: check whether the contour pixels are inside the nodule mask in pylidc vs ours (the outline-outside-nodule convention). Quantify the Dice/area difference on small (3–6 mm) nodules, where 1 px matters most.
7. **Empty-reader frequency**: in our crops, measure the distribution of #readers-marking ∈ {1,2,3,4}. Compare with the Prob-U-Net statement that up to 3 masks can be empty. Make sure each split contains disagreement slices.
8. **Patient leakage check**: assert there is no PatientID overlap across LIDC splits, including the 8 two-timepoint patients. Also check for near-duplicate adjacent-slice crops across splits (should be impossible after grouping).
9. **Polygon augmentation invariants**: property-test torchvision v2 `KeyPoints` under flip/rotate/scale/crop. After canonicalization, check that rasterize(transform(poly)) ≈ transform(rasterize(poly)) (IoU > 0.98) and that polygon validity is preserved.

### Recommendation
- **RIGA split**: train BinRushed + MESSIDOR (655), test Magrabia (95), exactly as MRNet/TAB/MSE-Nets do. Carve a validation set (about 10%) out of the 655, grouped by image and stratified by subset. Report disc and cup soft Dice at {0.1, 0.3, 0.5, 0.7, 0.9}, both "Average" and "Mean Voting", at 256×256. Derive the VLM probability map from K samples. Try to obtain MRNet's preprocessed release to match its crop. The ROADMAP RIGA line is correct.
- **LIDC split**: keep the patient-level 70/10/20 split as the primary dev split, stratified with StratifiedGroupKFold on agreement and size bins. Build crops in protocol-A geometry: 0.5 mm, 180×180 lesion-centred, with raw HU stored and the window applied at render. Use only ≥3 mm contours, and keep non-marking readers as explicit empty masks. For external comparability, add (a) the HPU-Net public crops (CC BY 3.0) as a sanity benchmark and (b) D-Persona-style 4-fold patient CV with GED/Dsoft/Dmax/Dmatch if the MedicalMatting slice list can be obtained.
- **Parser**: keep our own XML parser as primary, but gate it on parity with pylidc clustering and consensus on a sample (hypothesis 5). pylidc is unmaintained since 2021 and depends on an SQLAlchemy ORM, which makes it heavy as a core dependency.
- **Storage**: use one canonical Parquet HF dataset per source, with embedded image bytes, a per-rater polygon struct column, `group_id`, split columns and a SHA-256 manifest committed in git. Generate TRL / verl / LLaMA-Factory views with small exporters. Keep RIGA private (no redistribution). Skip DVC, webdataset, Zarr and Lance at current scale; revisit Lance if we exceed about 1M rows or need embedding search for dedup.
- **Augmentation**: torchvision v2 + `tv_tensors.KeyPoints` for polygons, with a post-transform canonicalize-and-validate step. Avoid AlbumentationsX (AGPL).
- **Next track**: IMA++ first (2D, polygon-friendly, official splits, 2–5 raters), then NuCLS (CC0) for multi-instance polygons. QUBIQ 2D is a small, classic Q-Dice benchmark. Defer MSSEG-2016 and Gleason2019.
- **Attribution checklist**: for LIDC, the "Data From LIDC-IDRI" DOI 10.7937/K9/TCIA.2015.LO9QL9SX, the Armato et al. 2011 paper and the NCI/FNIH acknowledgement. For RIGA, Almazroa et al., SPIE Medical Imaging 2018 (Proc. SPIE 10579). For NuCLS, Amgad et al. For IMA++, Abhishek et al. 2025 plus ISIC image licences.
