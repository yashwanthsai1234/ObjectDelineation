# A4: evaluation, baselines, data, augmentation, licence and the precision probe

Task A4 of `experiments/phase2/RESEARCH.md`. Written 2026-10-19.

Questions:
1. Multi-rater metrics, literature numbers on LIDC/RIGA, and whether the entry gates (`experiments/phase1/RESULTS.md` §8) are sensible.
2. Which baselines to run on our splits, and how to make "IoU near MedSAM" fair.
3. Augmentation that keeps targets certified (I7), plus mixing and size of RIGA vs LIDC.
4. RIGA and LIDC licences, and whether a rented GPU is "redistribution".
5. A proper design for the precision probe.

Evidence labels:
- **[local]** measured in this repo today. Scripts and JSON are in `experiments/phase2/lab/A4/` (copied from the scratch sandbox). Every script that reads `data/` runs with `.venv/bin/python -I` from the repo root.
- **[paper] / [card] / [code] / [page]** read in the primary source this session.
- **[agent]** gathered by one of three research sub-agents from primary pages. They report the table each number came from. I spot-checked the load-bearing ones and mark them **checked**.
- **[unverified]** not confirmed.

---

## 0. Bottom line

1. **The frame gives the answer away.** `targets.frame` centres the crop on the raters' union, and the side S follows the union's size. An image-blind model passes the RIGA geometry gates: a disc at the frame centre, sized by the train median for that S. **[local] S5**
   - RIGA, median IoU_tc: disc 0.904 (gate 0.85), cup 0.783 (gate 0.75); 95–100 % of items at IoU_tc ≥ 0.5 (gate 70 %).
   - LIDC: 0.588 (gate 0.60), 65 % at ≥ 0.5 (gate 70 %), just below both.
   - With the crop window shifted at random (the union still inside), the same prior falls to LIDC 0.00, disc 0.32, cup 0.07.
   - Fix: gate on a jittered-frame val, require a paired win over this prior, and train with frame jitter.
2. **The precision probe's "≤ 2 px" is far too loose for LIDC.** **[local] S2–S4**
   - LIDC nodules have a median diameter of 10.9 px at frame resolution.
   - Readers disagree by a median ASSD of 0.93 px and HD95 of 2.0 px.
   - A uniform 2 px error on a 10.9 px disc gives IoU 0.40.
   - Moving every vertex by σ = 2 px gives IoU 0.63, and 49 % of the polygons become illegal. A 2 px shift gives 0.61.
   - Re-specify the probe (§6): on synthetic shapes with exact ground truth, median ASSD ≤ 0.6 px, IoU ≥ 0.80 and legal ≥ 95 % on 6–16 px objects, with uncentred frames.
3. **On LIDC at this resolution, shape-blind blobs are as good as a second reader.** Against each reader, an equal-area disc at that reader's centroid scores median IoU 0.806 and ASSD 0.61 px. Another reader scores 0.71 and 0.93. **[local] S3**
   - So LIDC geometry measures mostly localisation and size; topology and calibration are where the claim lives.
   - Our LIDC entry gate (0.60) is 78 % of the readers' leave-one-out IoU_tc (0.768). That is sensible only once the centring shortcut is removed.
4. **Phase 0's MedSAM on LIDC (IoU 0.33) was a pipeline artefact.**
   - Published box-prompted numbers are much higher. On LUNA16 (LIDC images, window [−1000, 400]): Dice 0.79 for MedSAM, 0.89 for SAM2-tiny, 0.88 for MedSAM2 ([paper] 2608.16832, **checked**). Elsewhere: MedSAM Dice 0.773 / IoU 0.668 on LIDC.
   - Re-run on our frames, with the box from the raters' union, MedSAM scores **median IoU_tc 0.78 on LIDC** (96.7 % ≥ 0.5; IoU vs vote 0.63), 0.94 on the RIGA disc and 0.85 on the cup. That is at or above the human leave-one-out values. Official min-max preprocessing changes nothing (0.778 vs 0.780) **[local] S7**. Phase 0's 0.33 came from its crop and scoring, not from MedSAM. With an oracle box, MedSAM is a high bar: 0.78 against the 0.60 gate.
   - The "IoU near MedSAM" gate must use MedSAM re-run on our frames, same val split, same `score_mask` path. MedSAM gets an oracle box, so add SAM 2.1 (Apache-2.0) and the image-blind prior next to it.
5. **Literature numbers are not comparable to ours without protocol work.**
   - LIDC has four incompatible protocols:
     - Kohl v1, random 60/20/20 on 128² crops;
     - Kohl v2, patient split;
     - D-Persona, 4-fold over 214 subjects;
     - CIMD's own setup.
   - CAR shows the split alone moves GED100 from 0.243 to 0.228.
   - RIGA is uniform (655 train / Magrabia 95 test, soft Dice at 0.1–0.9, 256²). Our RIGA val/test are the same Magrabia-95 hospital split, at native disc-crop resolution, not resized to 256.
   - Our human leave-one-out ceilings, measured on our frames: LIDC 0.768, disc 0.945, cup 0.815 IoU_tc.
6. **Report these metrics:**
   - GED at K = 16 (and 100 where feasible);
   - Hungarian-matched IoU (HM-IoU) at K = 16;
   - D_max;
   - soft Dice (thresholds 0.1–0.9);
   - NSD at τ = the inter-reader 95th-percentile distance (Nikolov). τ ≈ 2 px for LIDC at frame resolution; τ per structure for RIGA;
   - HD95.

   `objdel.eval` has GED, D_max and soft Dice already; HM-IoU and NSD/HD95 are missing. None of the multi-rater LIDC/RIGA papers report NSD or HD95; Metrics Reloaded makes NSD the default boundary metric.
7. **Augmentation that keeps certification is cheap.** **[local] S1**
   - Rebuilding targets from D4-transformed masks certifies 100 % (798 RIGA and 4,858 LIDC target×op cases).
   - Cost: 23.6 ms per RIGA target and 0.8 ms per LIDC target on one core. All 7 non-identity ops over the whole train set take about 11 CPU-minutes.
   - Naively transforming vertices is illegal on every reflection unless re-canonicalised: the orientation flips.
   - Even re-canonicalised, naive vertices lose the IoU ≥ 0.99 certificate on 7.3 % of RIGA cups. LIDC is pixel-exact.
   - The image must be transformed at frame resolution, before the Pillow upscale. Rotations and transposes do not commute with it (up to 16 grey levels); flips do.
8. **Which augmentation pays off.** VLM grounding reports (Qwen2/2.5/3-VL, Florence-2, PaliGemma) use little or none. PaliGemma got the same anti-overfitting effect from label smoothing, dropout and crop/zoom. Our data is small and repeated (all 6 RIGA raters are rows each epoch), so:
   - frame translation and scale jitter are the important augmentation;
   - D4 transforms for LIDC; horizontal flip only for RIGA;
   - CT window jitter (L ± 60, W ± 80 HU) and mild photometric jitter for fundus;
   - no elastic deformations or non-90° rotations in SFT v1.
9. **RIGA's terms conflict; on a Vast community host, treat RIGA as not allowed.** **[page]**
   - The Deep Blue record's licence field and its DataCite metadata say **CC BY-NC 4.0**, which permits non-commercial redistribution.
   - The bundled `RIGAdataset_readme.pdf` ("Use and Access") says "Copy, redistribution, and any unauthorized commercial use are prohibited". MESSIDOR's own terms add that use "is restricted to individuals or organizations that obtained the database directly from" ADCIS.
   - Copying to a machine we alone control is processing, not redistribution, under any ordinary reading.
   - A Vast community host is a stranger's machine with root and physical access. Vast itself only promises Docker isolation there, and steers sensitive data to Secure Cloud.
   - Recommendation: RIGA only on Vast **Secure Cloud** (ISO 27001 datacentres with DPAs), encrypted at rest, deleted on teardown. Email the depositor for written permission, and register at ADCIS for MESSIDOR.
   - LIDC is CC BY 3.0 plus the TCIA policy (no re-identification): any host is fine.
10. **Baselines to run, ranked:**
    1. the image-blind prior (free);
    2. MedSAM ViT-B, box, on our frames;
    3. SAM 2.1-tiny, box (Apache-2.0; best published zero-shot LUNA16 box score);
    4. a supervised 2D nnU-Net v2, label-sampling and majority-vote variants. This is the specialist reference, about 1–3 GPU-hours per dataset at 250 epochs [agent estimate];
    5. MedSAM2 (weights for research and education only).

    Reject:
    - SAM 3 text prompts: 11.9 % external Dice untuned;
    - BiomedParse: trained on LIDC-IDRI, so contaminated; NC licence;
    - LISA, GLaMM, MedPLIB: LIDC Dice ≤ 0.006;
    - MedSAM-Agent: its LIDC numbers contradict each other.

---

## 1. Local measurements [local]

All runs used the val split of `data/store` (manifest `37fe7665…`), Mac M3 under load, one process each. Frame resolution means the `targets.frame` crop: RIGA S = 200–500 px, LIDC S = 50–125 px, scale 1 everywhere.

### S1. Dihedral augmentation: naive vertex transform vs rebuilt targets (`aug_cost.py`)

Setup: 21 RIGA val items (114 rater targets) and 217 LIDC val items (694 targets). For each of the 7 non-identity D4 ops:
- **naive**: transform the original target's vertices, then `canonicalize`;
- **rebuilt**: `build_target` on the transformed masks, with the same `min_areas`.

| | RIGA disc | RIGA cup | LIDC nodule |
|---|---|---|---|
| target × op cases | 798 | 798 | 4,858 |
| rebuilt: certified (legal, signature, IoU ≥ 0.99) | 100 % | 100 % | 100 % |
| rebuilt vertices / original, median (p95) | 1.00 (1.00) | 1.00 (1.07) | 1.00 (1.00) |
| naive + canonicalize: legal / signature kept | 100 / 100 % | 100 / 100 % | 100 / 100 % |
| naive: pixel-identical to the transformed mask | 0 % | 0 % | **100 %** |
| naive: still meets the certificate (IoU ≥ 0.99) | 100 % | **92.7 %** | 100 % |
| naive IoU min / median | 0.990 / 0.992 | 0.987 / 0.991 | 1.0 / 1.0 |
| naive **without** canonicalize: legal | 3/7 of ops (the rotations; every reflection reverses ring orientation) | same | same |
| build time per target, 1 core | 23.6 ms | (per item, both objects) | 0.81 ms |

Image check: transforming the frame-resolution crop and then upscaling (Pillow bicubic, 512 px), against transforming the upscaled image:
- flips: identical;
- rot90, rot270, transpose, antitranspose: max 16 grey levels, mean 0.27. Pillow resamples in two separable passes with an 8-bit intermediate, so swapping axes changes the rounding.

Readings:
- I7's "99.8 % miss" is about pixel identity. In certificate terms, naive mirroring breaks 7 % of RIGA cups and nothing on LIDC: at S = 50–125 a pixel is 8–20 grid units, so τ ≤ 3 never puts an edge through a pixel centre.
- Rebuilding is simple and fast enough to precompute offline:
  - RIGA: 3,282 targets × 7 × 23.6 ms ≈ 9 min;
  - LIDC: 23,469 × 7 × 0.81 ms ≈ 2 min, single core.
- Store the frame-resolution crops (LIDC: 50–125 px) and upscale when collating. Storing 8 upscaled copies would be about 10 GB; the LIDC prepared images are 1.3 GB today.

### S2. Human ceiling on our frames (`rater_ceiling.py`, every val item)

Definitions:
- Pairs: every pair of voting raters.
- LOO: each rater against the others, using the gate's own metric. IoU_tc is the best IoU to another rater with the same size-filtered signature, else IoU to the others' vote.
- ASSD, HD95 and NSD use 4-connected boundary pixels on the frame grid.

| | LIDC nodule | RIGA disc | RIGA cup |
|---|---|---|---|
| pairs / LOO raters | 4,109 / 3,444 | 945 / 378 | 945 / 378 |
| pairwise IoU p10 / p50 / p90 | 0.49 / **0.72** / 0.86 | 0.87 / **0.92** / 0.95 | 0.47 / **0.72** / 0.86 |
| pairwise ASSD p50 (p90), px | **0.93** (1.95) | 3.73 (6.60) | 7.31 (16.5) |
| pairwise HD95 p50, px | **2.02** | 8.72 | 15.6 |
| NSD@1 px / @2 px, p50 | 0.77 / 0.97 | 0.22 / 0.36 | 0.09 / 0.16 |
| **LOO IoU_tc p10 / p50** (the gate metric's human value) | 0.55 / **0.768** | 0.913 / **0.945** | 0.685 / **0.815** |
| LOO IoU vs vote p50 | 0.731 | 0.937 | 0.767 |
| LOO share with IoU_tc ≥ 0.5 | 92.9 % | 100 % | 99.2 % |
| object equivalent diameter p10 / p50 / p90, px | 6.6 / **10.9** / 22.4 | 138 / 203 / 232 | 63 / 95 / 128 |

LIDC by diameter (reader pairs, medians):

| Diameter | 0–6 px | 6–10 | 10–15 | 15–25 | > 25 |
|---|---|---|---|---|---|
| pairs | 135 | 1,445 | 1,176 | 1,095 | 258 |
| IoU | 0.60 | 0.69 | 0.73 | 0.75 | 0.76 |
| ASSD, px | 0.69 | 0.74 | 0.89 | 1.19 | 1.71 |
| HD95, px | 1.41 | 1.41 | 2.0 | 3.0 | 5.0 |

What a uniform boundary error does to a disc of the median size (arithmetic): IoU = (r − e)² / r².

| Diameter | e = 1 px | e = 2 px |
|---|---|---|
| LIDC 10.9 px | 0.67 | **0.40** |
| RIGA cup 95 px | 0.96 | 0.92 |
| RIGA disc 203 px | 0.98 | 0.96 |

### S3. Shape-blind baselines on real LIDC (`probe_proto.py`, every 3rd val item, 1,140 reader masks)

| Prediction vs that reader | IoU p50 | ASSD p50 (p90), px | HD95 p50 | NSD@1 |
|---|---|---|---|---|
| equal-area disc at the reader's centroid | **0.806** | **0.61** (1.64) | 1.35 | 0.95 |
| ellipse inscribed in the reader's box | 0.800 | 0.62 (1.71) | 1.41 | 0.88 |
| convex hull | 0.752 | 0.78 (1.27) | 1.41 | 0.83 |
| **another reader** | 0.710 | 0.94 (1.98) | 2.20 | 0.77 |

On synthetic radial-Fourier shapes my first generator was too round: the disc baseline reached IoU 0.90. A probe built on it would not test shape, so the generator must be calibrated against real LIDC (§6).

Moving every vertex of a 64-vertex synthetic polygon by N(0, 0.5 px) made 85 % of them illegal. Dense vertices plus coordinate noise means self-intersections, which matters for the legality gate.

### S4. What coordinate error costs on real targets (`noise_model.py`)

Setup: LIDC val, every 3rd item, 1,140 certified targets with median 10 vertices (p90 23); RIGA cup, 192 targets with median 16. Perturbations: every vertex moved by N(0, σ px); a whole-polygon shift of d px; scaling about the centroid. "Legal" is the polycore legality rate; IoU is over legal outputs.

| Perturbation | LIDC legal | LIDC IoU p50 (p10) | LIDC ASSD p50 | RIGA cup legal | RIGA cup IoU p50 |
|---|---|---|---|---|---|
| vertex σ = 0.25 px | 99.3 % | 0.97 (0.93) | 0.10 | 99.0 % | 0.989 |
| σ = 0.5 px | 95.7 % | 0.90 (0.80) | 0.30 | 97.9 % | 0.985 |
| σ = 1 px | 83.7 % | 0.79 (0.64) | 0.61 | 94.8 % | 0.971 |
| σ = 2 px | **50.9 %** | **0.63** (0.41) | 1.11 | 95.3 % | 0.947 |
| shift 0.5 / 1 / 2 px | 100 % | 0.91 / 0.77 / **0.61** | 0.26 / 0.71 / 1.22 | 100 % | 0.984 / 0.971 / 0.944 |
| scale −20 / −10 / +10 / +20 % | 100 % | 0.64 / 0.81 / 0.83 / 0.69 | ≈ 0.5–1.0 | 100 % | 0.64 / 0.81 / 0.83 / 0.70 |

Readings:
- To pass the LIDC legality gate (≥ 90 %) through noise alone, per-vertex error must stay at or below about 0.6 px at frame resolution. That is 12 grid units at S = 50.
- The geometry gate (median 0.60) tolerates about a 2 px shift.
- So **legality, not IoU, sets the precision LIDC needs.** RIGA tolerates 2 px everywhere.

### S5–S6. Image-blind prior against the gates (`null_prior.py`, `null_prior_jitter.py`)

The prior: a disc at the frame centre, with diameter the train median for that frame side. LIDC: 11.2 / 27.5 / 32.7 px at S = 50 / 100 / 125. RIGA: disc 138 / 210 and cup 70 / 99 px at S = 250 / 500. Scoring is `reward.score_mask` IoU_tc on every val item, which is the gate path.

| | LIDC nodule | RIGA disc | RIGA cup |
|---|---|---|---|
| gate (median IoU_tc / share ≥ 0.5) | 0.60 / 70 % | 0.85 / 70 % | 0.75 / 70 % |
| prior, centred frames | 0.588 / 65.3 % | **0.904 / 100 %** | **0.783 / 95.2 %** |
| prior, frames shifted at random (union kept inside) | 0.000 / 1.4 % | 0.320 / 25.4 % | 0.072 / 7.9 % |

### S7. MedSAM on our frames (`medsam_check.py`)

Setup:
- 121 LIDC val items (every 9th) and 19 RIGA val items (every 3rd, with images).
- Input: `frame(item).image`, which is the LIDC 512-px upscale.
- Box: the raters' union, +1 frame px for LIDC and +5 px for RIGA, scaled to the image.
- Output: the mask resampled to the frame at pixel centres.
- Scoring: `reward.score_mask` IoU_tc on `ExpertCache.build(raters)`.
- Hardware: CPU, 2 threads, 8.6 s per image embedding (2,408 s in total). Peak RSS 3.13 GB, slightly over the 3 GB budget; MPS peaked at 4.5 GB footprint and was not used.

| Prediction | LIDC median IoU_tc (≥ 0.5) | LIDC IoU vs vote | RIGA disc IoU_tc | RIGA cup IoU_tc |
|---|---|---|---|---|
| MedSAM, HF port (rescale 1/255, as in Phase 0) | **0.780** (96.7 %) | 0.634 | 0.941 | 0.847 |
| MedSAM, official min-max | 0.778 (96.7 %) | 0.632 | 0.944 | 0.834 |
| ellipse inscribed in the same box (null with the same information) | 0.602 (81.0 %) | 0.444 | 0.867 | 0.752 |
| human LOO, for reference (S2) | 0.768 (92.9 %) | 0.731 | 0.945 | 0.815 |
| Phase 0 MedSAM (3× window ≥ 64 px, strict-majority scoring) | 0.33 (IoU vs majority) | — | 0.92 | 0.67 |

MedSAM on LIDC by union diameter, median IoU_tc:

| Union diameter | n | MedSAM | Box ellipse |
|---|---|---|---|
| < 10 px | 47 | 0.77 | 0.55 |
| 10–16 px | 34 | 0.79 | 0.63 |
| > 16 px | 40 | 0.79 | 0.67 |

Readings:
- Phase 0's 0.33 is an artefact of its crop and scoring, not of MedSAM's preprocessing.
- With an oracle box on our frames, MedSAM matches human leave-one-out agreement on LIDC and on the RIGA disc, and exceeds it on the cup.
- The union box encodes the raters' extent, so this is an upper bound for box prompting.

---

## 2. Multi-rater evaluation: metrics, literature numbers, gate sanity

### 2.1 Metrics in use (2024–2026) and what we have

| Metric | Definition / convention | Who reports it | In `objdel.eval`? |
|---|---|---|---|
| GED (D²_GED) | 2E d(S,Y) − E d(S,S′) − E d(Y,Y′), d = 1 − IoU, d(∅,∅) = 0. **Sample count must be stated** (16 / 32 / 50 / 100) | Kohl 2018 App. B; every LIDC paper | yes (V-statistic, K samples) |
| HM-IoU | Hungarian matching of samples to raters (raters duplicated to match K), mean IoU of the matched pairs. Introduced because GED "unduly rewards sample diversity" | HPU-Net 2019; CAR, CCDM, MoSE, BerDiff | **no: add it** |
| D_max | best IoU (Dice in CIMD) of any sample to any rater | CIMD, D-Persona, ProSeg | yes (per sample, best rater) |
| D_soft / soft Dice | Dice of mean(samples) ≥ t vs mean(raters) ≥ t, t ∈ {0.1, …, 0.9}; RIGA uses 5 thresholds | MRNet, TAB (RIGA); D-Persona (LIDC) | yes (5 thresholds) |
| D_match | per-annotator Dice after matching personalised outputs to annotators | D-Persona, ProSeg, ProSona | n/a (no rater identity in our policy) |
| S_NCC | NCC between the expected cross-entropy map and the rater variance | PHiSeg, ValUES | no (needs probabilities; optional, from K samples) |
| CI score | 3 S_c D_max D_a / (S_c + D_max + D_a) | CIMD | no (skip: not standard) |
| NSD@τ | (|S_P ∩ B_G(τ)| + |S_G ∩ B_P(τ)|) / (|S_P| + |S_G|); τ from inter-observer variation (95th percentile) | Nikolov 2018; Metrics Reloaded default boundary metric | **no: add it** |
| HD95 | max of the two directed 95th percentiles | Metrics Reloaded (preferred over HD for outliers) | **no: add it** |
| Boundary IoU | IoU of the d-band (2 % of the diagonal) | Cheng 2021 | not useful here: it equals mask IoU for small objects |

Sources: Kohl https://arxiv.org/abs/1806.05034; HPU-Net https://arxiv.org/abs/1905.13077 (App. B.1 criticism); CIMD https://arxiv.org/abs/2304.04745; D-Persona https://arxiv.org/html/2403.13417v1; PHiSeg https://arxiv.org/html/1906.04045v2; ValUES https://arxiv.org/abs/2401.08501; Nikolov https://arxiv.org/abs/1809.04430; Metrics Reloaded https://arxiv.org/abs/2206.01653; Metric pitfalls https://arxiv.org/abs/2302.01790; Boundary IoU https://arxiv.org/abs/2103.16562 [agent; definitions checked against M8's existing GED / soft-Dice citations].

What Metrics Reloaded says that matters for us [agent]:
- fuzzy multi-observer references are "out of scope";
- avoid Dice alone for small structures;
- pair an overlap metric with NSD, with τ set from inter-rater variability;
- prefer HD95 or MASD over HD;
- consider treating small structures as detection.

### 2.2 LIDC numbers (protocol first: they are not comparable across protocols)

| Method | Protocol | Numbers (table) | URL |
|---|---|---|---|
| Prob. U-Net (NeurIPS 2018) | v2 (patient split, 722/144/144) | GED at 1/4/8/16 samples: 0.811 / 0.388 / 0.321 / 0.287 | https://arxiv.org/abs/1806.05034 |
| HPU-Net (2019) | v2 | HM-IoU 0.53; GED 0.27 | https://arxiv.org/abs/1905.13077 |
| PHiSeg (MICCAI 2019) | v1, 100 samples | GED 0.2248, S_NCC 0.8453 (Prob. U-Net 0.2393 / 0.7749) | https://arxiv.org/html/1906.04045v2 |
| SSN (NeurIPS 2020) | v1, 100 | GED 0.225; deterministic U-Net 0.607; expert diversity 0.399 | https://arxiv.org/abs/2006.06015 |
| CAR (ICCV 2021) | v2 | GED 16/50/100: 0.264 / 0.248 / 0.243; HM-IoU16 0.592 | https://arxiv.org/abs/2006.13144 |
| MoSE (ICLR 2023) | v1 | GED 16/50/100: 0.218 / 0.195 / 0.189; M-IoU16 0.624 | https://arxiv.org/abs/2212.07328 |
| CIMD (CVPR 2023) | own split, 4 samples | GED 0.321, CI 0.759, D_max 0.915 | https://arxiv.org/abs/2304.04745 |
| CCDM (ICCV 2023) | v1 / v2 | v1 GED16 0.212, HM-IoU16 0.623; v2 GED16 0.239, HM-IoU16 0.598 | https://arxiv.org/abs/2303.08888 |
| BerDiff (2023) | standard | GED16 0.238, HM-IoU16 0.596; MedSegDiff 0.420 / 0.413 | https://arxiv.org/abs/2304.04429 |
| D-Persona (CVPR 2024) | 214 subjects, 4-fold | GED 0.1358, D_soft 90.45, D_max 91.37, D_match 91.33 | https://arxiv.org/html/2403.13417v1 |
| DiffOSeg (MICCAI 2025) | D-Persona | GED50 0.0745, D_soft 92.37 | https://arxiv.org/abs/2507.13087 |
| ProSeg (2025/26) | D-Persona | GED 0.1077, D_soft 91.62, D_max 91.46 | https://arxiv.org/html/2512.00748 |
| ATFM (AAAI 2026) | v1 | GED 16/32/100: 0.206 / 0.188 / 0.162; HM-IoU32 0.667 | https://arxiv.org/abs/2511.06857 |
| Tyche (CVPR 2024) | out of domain, Dice-distance GED | not comparable | https://arxiv.org/abs/2401.13650 |

All [agent]; the protocol split and CAR's 0.243 vs 0.228 are from CCDM §4.1 and CAR Table 1.

Human references:
- v1 expert pairwise distance 0.399, i.e. a mean pairwise IoU of about 0.60 including empty masks (SSN §5.2).
- Ours, voting raters only: pairwise IoU median 0.72 (S2).

### 2.3 RIGA numbers (train BinRushed + MESSIDOR 655, test Magrabia 95, soft Dice, 256² inputs)

| Method | Disc / cup soft Dice (%) | URL |
|---|---|---|
| MRNet (CVPR 2021) | vs majority vote 97.55 / 87.20; vs STAPLE 96.26 / 86.37 | https://openaccess.thecvf.com/content/CVPR2021/papers/Ji_Learning_Calibrated_Medical_Image_Segmentation_via_Multi-Rater_Agreement_Modeling_CVPR_2021_paper.pdf |
| TAB (MICCAI 2023) | average over raters 96.70 / 85.52; mean voting 97.82 / 88.22; per-rater cup 82.6–89.4 | https://arxiv.org/abs/2306.01340 |
| Self-Calib (MICCAI 2022) | MV 97.22 / 88.24 | https://arxiv.org/abs/2206.05092 |
| UMA-Net (2023) | 97.95 / 90.70 (compared against a self-fusion reference, so not like for like) | https://arxiv.org/abs/2304.00466 |
| Multi-rater prompting (2024) | hard Dice, 224 px input: 96.58 / 85.70 | https://arxiv.org/abs/2404.07580 |

Notes:
- All [agent].
- None of these papers states whether its 256² input is a disc crop or the whole fundus [unverified].
- None of LIDC's ambiguity papers (D-Persona, ProSeg, DiffOSeg) uses RIGA.

### 2.4 How to report ours so a reader can compare

- **Same split, our resolution.** RIGA test = Magrabia 95 (86 unique eyes), soft Dice "mean voting" at {0.1, 0.3, 0.5, 0.7, 0.9}, with the prediction map from K = 16 samples. State the difference: disc crop at native 200–500 px, not 256².
- **LIDC.** Report on our patient-level split, with the protocol stated. Give:
  - GED16 and HM-IoU16 (Kohl conventions, d(∅, ∅) = 0);
  - D_max and D_soft;
  - our human references from S2.

  Do not put our numbers in the same table column as v1/v2 numbers. Optionally evaluate on the Kohl v2 test crops (CC BY 3.0, public) for one cross-reference row [suggestion].
- **Add HM-IoU, NSD@τ and HD95 to `objdel.eval`** (Phase 2 or 9). Use τ = 2 px for LIDC, the inter-reader HD95 median at frame resolution, in Nikolov's sense. For RIGA, τ per structure: disc 9 px, cup 16 px from S2. Alternatively τ = 1 mm / 2 mm in physical units for publication.

### 2.5 Are the entry gates sensible?

| Gate (RESULTS §8) | Human LOO on our frames | Image-blind prior (centred) | Verdict |
|---|---|---|---|
| parse ≥ 98 %, legal ≥ 90 % | targets 100 % | — | fine. Legality is the binding precision constraint on LIDC (S4) |
| median IoU_tc ≥ 0.85 disc | 0.945 | **0.904** | **not informative on centred frames**; keep the number, but measure on jittered frames and require a paired win over the prior |
| median IoU_tc ≥ 0.75 cup | 0.815 | **0.783** | same |
| median IoU_tc ≥ 0.60 LIDC | 0.768 | 0.588 | 78 % of human; sensible only on jittered frames |
| ≥ 70 % with IoU_tc ≥ 0.5 | LIDC 92.9 %, RIGA ≥ 99 % | LIDC 65 %, RIGA 95–100 % | same issue |
| topology-any ≥ 90 % on unanimous | targets 100 %; LOO topology match 98.5 % on LIDC | prior: always one simple part. A1 found 0 unanimous non-simple LIDC val items, and RIGA has no topology disagreement, so the prior scores 100 % | also uninformative: a "one simple blob" policy passes |
| precision probe ≤ 2 px | reader ASSD 0.93, HD95 2.0 | the disc baseline already has ASSD 0.61 | too loose; re-specify (§6) |

Two proposals:
1. Add a **jittered val** to the gates: the same val items, window shifted with a fixed seed, union kept inside, targets and caches rebuilt.
2. Add a **null-model gate**: on centred frames the model's per-item IoU_tc beats the prior's on ≥ 60 % of items, sign test p < 0.01.

The topology gate also needs a non-trivial counterpart. A1 proposes a topology probe; from the eval side, report topology-any on **disagreement** items, where "simple blob" is wrong for some raters, plus min-mode retention and coverage, as Phase 3 will.

---

## 3. Baselines on our splits

### 3.1 Candidates ([agent] unless marked; licences checked where load-bearing)

| Model | Prompt | Licence (code / weights) | Size | Evidence on our objects | Use? |
|---|---|---|---|---|---|
| **Image-blind prior** | none (frame) | ours | 0 | S5 | **yes, mandatory null** |
| **MedSAM ViT-B** | box | Apache-2.0 / Apache-2.0 | 94 M, 375 MB | LIDC Dice 0.773 / IoU 0.668 (https://arxiv.org/html/2602.03320); LUNA16 Dice 0.79 (https://arxiv.org/html/2608.16832, **checked**); REFUGE Dice 0.846; our Phase 0 disc 0.92 / cup 0.67 / LIDC 0.33 | **yes**, re-run on our frames |
| **SAM 2.1 hiera-tiny / small** | box | Apache-2.0 / Apache-2.0 | 39–46 M, 156 MB | LUNA16 Dice 0.89 (**checked**); REFUGE 0.717 | **yes**: the licence-clean strongest box model on nodules |
| MedSAM2 | box | Apache-2.0 / "research and education purposes" only ([card], **checked**) | 39 M | LUNA16 0.88 (**checked**) | optional; research use only |
| nnU-Net v2 (2D) | none, supervised | Apache-2.0 | per dataset | REFUGE OD 0.947 / OC 0.849 (https://arxiv.org/html/2502.06220v1); multi-rater RIGA 93.1 / LIDC 3D 56.4 (https://arxiv.org/html/2606.16868) | **yes**: specialist reference, label-sampling and majority-vote variants |
| SAM 3 / 3.1 | text concept, box | SAM License (custom), gated | 848 M, 3.45 GB | text-only Dice 11.9 % on external medical data (https://arxiv.org/html/2601.10880v1) | no (box mode adds nothing over SAM 2.1; 3.45 GB, CUDA ≥ 12.6) |
| Medical SAM3 | text | base SAM License | 848 M | trained on PAPILA / other fundus; no LIDC; external 73.9 % | optional disc/cup text baseline only |
| BiomedParse v1 / v2 | text | Apache-2.0 / CC-BY-NC-SA-4.0, gated | 1.8 / 4.5 GB | **trained on LIDC-IDRI, REFUGE and G1020** (https://huggingface.co/datasets/microsoft/BiomedParseData) | no (contaminated on LIDC; NC) |
| LISA, GLaMM, MedPLIB | text | Apache-2.0 | 7–13 B | LIDC Dice 0.006 / 0.005 / 0.001 (https://arxiv.org/html/2602.03320) | no |
| UniBiomed | text | not stated | MLLM + SAM | LIDC Dice 0.50; REFUGE 0.79 | no (no clear licence) |
| MedSAM-Agent | text → SAM | — | Qwen3-VL-8B + MedSAM2 | headline LIDC 0.780 vs appendix 0.523 | no (inconsistent) |
| Florence-2 B / L | text → polygon | MIT | 0.23 / 0.77 B | none medical | A1's question (an SFT'd polygon VLM baseline) |

### 3.2 Making "IoU near MedSAM" fair

**Same crop.** The baseline sees `frame(item).image`: the RIGA crop, or the LIDC 512-px upscale of the lung-windowed crop. Its mask is resampled to the S×S frame at pixel centres (S7 does this). Phase 0 used a different crop, a 3× box window of at least 64 px.

**Same split, same scorer.** Val for the gate, Magrabia/LIDC test once for the paper. Score with `reward.score_mask` IoU_tc on `ExpertCache.build(raters)`, exactly the VLM's path.

**Same information, as far as possible.** The VLM gets a frame centred on the union; MedSAM gets a box. Report three rows:
1. MedSAM with the oracle tight union box (+1 px LIDC, +5 px RIGA): its upper bound;
2. MedSAM with a jittered box (0–20 px at 1024 scale, MedSAM's own training jitter, https://pmc.ncbi.nlm.nih.gov/articles/PMC10803759/);
3. the image-blind prior.

On jittered val frames, the VLM must localise and MedSAM still gets the box, so the comparison favours MedSAM; say so.

**Preprocessing.** Use the official MedSAM path: resize to 1024, min-max to [0, 1], no ImageNet normalisation. https://github.com/bowang-lab/MedSAM/blob/main/MedSAM_Inference.py lines 121–126 [code]. The HF port `flaviagiammarino/medsam-vit-base` rescales by 1/255 without min-max ([code] preprocessor_config.json). S7 tests both.

**Gate wording.** S7 puts MedSAM with an oracle box at 0.78 (LIDC), 0.94 (disc) and 0.85 (cup). That is at human leave-one-out level, so "near MedSAM" is stricter than the RESULTS §8 medians. Proposal:
- On **centred** val, where the frame carries about as much information as a box, the model's median IoU_tc is ≥ MedSAM-oracle-box − 0.10 per structure: LIDC ≥ 0.68, disc ≥ 0.84, cup ≥ 0.74.
- The model beats the image-blind prior by the paired test of §2.5.
- Jittered-val numbers are reported next to these and gated at the RESULTS §8 medians.
- The planner must decide whether LIDC's "near MedSAM" (0.68) or the 0.60 median binds.

Convert any published Dice to IoU before comparing: IoU = D / (2 − D).

### 3.3 Costs

| Baseline | Where | Cost |
|---|---|---|
| MedSAM / SAM 2.1, val + test (≈ 1,300 LIDC + 160 RIGA frames, two variants) | Mac CPU | measured 8.6 s per image embedding on CPU (S7, 2 threads, under load); about 3.5 h CPU for everything, or about 5 min on any Vast GPU |
| nnU-Net v2 2D, 1 fold × 2 label strategies × 2 datasets, 250 epochs | Vast 4090 | about 1–3 GPU-h per run [agent estimate from 12 h/fold at 1,000 epochs], i.e. about $2–6 in total |

The nnU-Net RIGA variant must run on Secure Cloud (§5).

---

## 4. Data and augmentation

### 4.1 What the literature does

- **VLM grounding.** No flip or crop is described for Qwen2-VL (https://arxiv.org/html/2409.12191), Qwen2.5-VL (copy-paste only; https://arxiv.org/html/2502.13923) or Qwen3-VL (https://arxiv.org/html/2511.21631, first ~100k characters read) [agent].
- **PaliGemma** (https://arxiv.org/html/2407.07726) [agent]: "we have not found any significant benefit from data augmentation". Plain resize overfits beyond 10 epochs unless crop/zoom or label smoothing plus dropout is used; they kept resize + label smoothing + dropout. The big_vision refcoco_seg transfer used label smoothing 0.3 for 100 epochs (https://github.com/google-research/big_vision/blob/main/big_vision/configs/proj/paligemma/transfers/refcoco_seg.py).
- **Florence-2** (https://arxiv.org/html/2311.06242): no augmentation in detection fine-tuning [agent].
- **Pix2Seq** (https://ar5iv.labs.arxiv.org/html/2109.10852) [agent]:
  - large-scale jitter [0.1, 3] from scratch;
  - sequence augmentation helps recall (AR 55.4 → 61.7) more than AP (43.7 → 44.7);
  - "600 bins" give zero quantisation error at 600 px. Our 1,000-unit grid on S ≤ 500 frames is exact by construction (targets.frame).
- **PolyFormer** (https://ar5iv.labs.arxiv.org/html/2302.07387): polygon vertex re-sampling augmentation, +3.72 mIoU [agent]. Not applicable as is: our targets are canonical certified rings.
- **Flips and left/right wording.** MDETR swaps the words; THflip+ skips flipping such captions (https://arxiv.org/html/2311.02536). Our prompts have no laterality words, so flips are safe [agent].
- **Medical.**
  - LIDC classification uses exact transforms only: 3-axis flips, 90° rotations, transposes (https://arxiv.org/pdf/2202.12515).
  - Prob. U-Net: elastic, rotation, shear, scale, random 128 crop from 180 (https://arxiv.org/abs/1806.05034).
  - nnU-Net default: rotation p 0.2, scale 0.7–1.4 p 0.2, no elastic, noise, blur, brightness and contrast 0.75–1.25, gamma 0.7–1.5, low-resolution simulation, mirroring (https://github.com/MIC-DKFZ/nnUNet/blob/master/nnunetv2/training/nnUNetTrainer/nnUNetTrainer.py) [agent].
  - AADG for optic disc/cup uses colour-only augmentation, because "shearing and translation … destroy the structural prior of OD/OC" (https://arxiv.org/abs/2207.13249) [agent].
  - Random window augmentation for CT: L ± 60, W ± 80 HU on raw HU (https://arxiv.org/html/2510.08116v1) [agent].

### 4.2 What this means for us

1. **Frame jitter is the augmentation that matters**, because the frame leaks the answer (S5–S6).
   - Translation: shift the window uniformly, keeping the union inside, with at least 2 px clearance.
   - Scale: draw the margin from [8 px, 60 %] so that S moves between divisors: LIDC 50 ↔ 100 ↔ 125, RIGA 200 ↔ 250 ↔ 500.
   - Rebuild targets from the re-cropped masks; that costs the same as S1.
   - Keep the frame invariants: S divides 1000, scale 1, S ≤ 500.
2. **D4 for LIDC; horizontal flip only for RIGA.**
   - Nodule appearance has no canonical orientation, and the exact LIDC transforms are standard.
   - A fundus disc crop has a physiological orientation: vessels leave vertically, and the cup sits temporally. A horizontal flip is a left/right eye swap, which occurs naturally. Vertical flips and 90° rotations make images that never occur at test time.
   - Assign a different op to each rater row of an image (all voting raters are rows in each epoch; `train/dataset.py` `sft_rows`). Then an image is never shown identically twice in one epoch.
3. **Photometric.**
   - LIDC: HU window jitter applied on raw HU, before `frame`'s windowing: L = −300 ± 60, W = 1400 ± 80, plus noise σ ≤ 10 HU. This needs `frame` to take the window as a parameter.
   - RIGA: brightness and contrast 0.75–1.25, gamma 0.7–1.5, p ≈ 0.15–0.3 (nnU-Net). No hue shift.
4. **Not in SFT v1:**
   - non-90° rotation, which resamples masks (nnU-Net uses argmax over interpolated one-hot labels) and needs recertification;
   - elastic deformation;
   - mixup, cutout or copy-paste (a pasted nodule changes the frame's ground truth);
   - random start vertex (Pix2Seq-v2). D10's canonical start is worth about 40 IoU points (A1).
5. **Precompute offline, deterministically.**
   - Per epoch e and row r, op = hash(seed, e, item, r). Store the frame-resolution PNG and target text.
   - `target_set_hash` per epoch view; K5/K6 re-checked by the builder (it raises when certification fails).
   - Cost: about 11 CPU-minutes per full D4 view (S1).

### 4.3 Size and balance, RIGA vs LIDC

- **SFT rows.** LIDC train 7,341 items / 23,469 targets; RIGA 547 items / 3,282 targets (RESULTS §3). RIGA answers are about 3× longer: Qwen median 309 vs 100 tokens.
- **Mixing options** (computed):

  | Rule | RIGA share of rows | RIGA share of answer tokens |
  |---|---|---|
  | proportional by targets | 12.3 % | 30 % |
  | proportional by items | 6.9 % | 19 % |
  | temperature T = 2 by targets (T5) | 27.2 % | 54 % |
  | temperature T = 2 by items | 21.4 % | 46 % |

  T5's temperature mixing: https://arxiv.org/abs/1910.10683 [agent].
- **Recommendation:** T = 2 by items. This is about 5 RIGA epochs per 2 LIDC epochs. Use frame jitter and the rater-row ops above for RIGA diversity, and label smoothing 0.1 as PaliGemma does.
- About 4 epochs of repetition is harmless in LM pretraining (Muennighoff, https://arxiv.org/abs/2305.16264). That is not SFT evidence, so watch RIGA val per epoch.
- **Synthetic shapes.**
  - Euclid finds low-level geometric perception poor, and curriculum training on synthetic geometry beats mixed training (https://arxiv.org/html/2412.08737).
  - VLMs are "blind" to overlapping primitives (https://arxiv.org/abs/2407.06581).
  - Formula-driven segmentation pretraining (SegRCDB, https://arxiv.org/abs/2309.17083) transfers for CNN segmenters.
  - Nobody shows synthetic-shape SFT transferring to medical polygon VLMs [agent: gap]. Use synthetic shapes as the **probe** (§6), and as a pre-SFT stage only if a model fails the probe (H-A4-6).

---

## 5. Licences

### 5.1 RIGA: what the documents say

| Document | Text | Source |
|---|---|---|
| Deep Blue record, "License" field (Wayback snapshot 2026-08-03, live page behind a Cloudflare 403) | "http://creativecommons.org/licenses/by-nc/4.0/" | https://web.archive.org/web/20260803015113/https://deepblue.lib.umich.edu/data/concern/data_sets/3b591905z [page] |
| DataCite metadata for DOI 10.7302/Z23R0R29 | `rightsList: by-nc/4.0` (twice) | https://api.datacite.org/dois/10.7302/Z23R0R29 [page] |
| `RIGAdataset_readme.pdf` (19.6 KB, in the record; Wayback 2023-01-21), section "Use and Access" | "The dataset can be used, free of charge, for research and educational purposes. **Copy, redistribution, and any unauthorized commercial use are prohibited.** Any researcher reporting results that use this dataset must acknowledge the MESSIDOR program, Magrabi Eye Center in Saudi Arabia and Bin Rushed Ophthalmic center in Saudi Aribia. Users of the RIGA dataset must cite" Almazroa et al., Proc. SPIE 10579, 105790B (2018), doi 10.1117/12.2293584 | https://web.archive.org/web/20230121082108/https://deepblue.lib.umich.edu/data/downloads/g732d957m [page; PDF text extracted locally] |
| CC BY-NC 4.0 | licensees may "copy and redistribute the material in any medium or format" for non-commercial purposes, with attribution | https://creativecommons.org/licenses/by-nc/4.0/ |
| MESSIDOR (source of 460 RIGA images) | "Copy, redistribution, and any unauthorized commercial use are prohibited"; "the use of this database is restricted to individuals or organizations that obtained the database directly from this website"; acknowledge "Kindly provided by the Messidor program partners" | https://www.adcis.net/en/third-party/messidor/ [page] |

The readme's sentence is copied from MESSIDOR's terms.

### 5.2 Is copying RIGA to a rented GPU "redistribution"?

**Plain meaning.** Redistribution means making the data available to another person or entity. Copying data onto infrastructure that we alone use, to process it ourselves, is use. "Copy … prohibited" cannot mean every copy: downloading is a copy. Read with "redistribution", it means copies made for others. On that reading, a rented machine used only by us is not redistribution. The CC BY-NC 4.0 field permits even non-commercial redistribution.

**Where the risk sits.** On a Vast community host, a third party (the host operator) has root and physical access to the disk, and the data stays on their storage unless we delete it.
- Vast promises only container isolation there. Its security FAQ says "Individual hosts may have less formal security measures" and advises "Use Secure Cloud certified providers only, Encrypt sensitive data at rest" for sensitive data (https://docs.vast.ai/guides/reference/faq/security).
- Secure Cloud partners are "ISO 27001-certified" and "sign expanded hosting agreements that include a separate Data Processing Agreement" (https://vast.ai/article/security-and-compliance-at-vast-ai).
- With such a provider the host is a contracted processor, like AWS. With a community host it is an unknown individual, and an uncontrolled copy on their disk is the case closest to "redistribution".

**The bigger, separate risks:**
1. "Free of charge, for research and educational purposes" and "unauthorized commercial use". If this project is done for a company, ask whether that counts as commercial use. This decides whether RIGA can be used at all, wherever it runs.
2. MESSIDOR's direct-download clause: our MESSIDOR images came from Deep Blue, not ADCIS.

**Recommendation (ranked):**
1. Email the RIGA depositor (contact in the readme: almazroa@med.umich.edu) and ask for written confirmation of research use on a rented GPU in an ISO 27001 datacentre. Record the reply in `docs/`. Cost: one email; it settles the readme/CC conflict.
2. Until then, RIGA trains only on Vast **Secure Cloud** instances:
   - ship the RIGA crops as an encrypted archive (age/GPG), decrypted into the container at start, key never on disk;
   - delete the instance and volume after the run;
   - log instance ID, provider and datacentre.
   - LIDC can go anywhere.
3. Register at ADCIS for MESSIDOR (free) to satisfy its direct-download clause.
4. Never publish RIGA-derived crops, masks or polygons; publish code, hashes and the split lists only. Phase 9's "certified polygon versions of the datasets" release covers LIDC only.
5. Fallback: RIGA SFT on the Mac only (547 images; the SmolVLM / LFM LoRA footprint was 6.2 GB in K12). This constrains model choice and is a last resort.

Do not use third-party RIGA mirrors (the MRNet Google Drive / Baidu copies, the 2026 Zenodo re-processing at https://zenodo.org/records/23158636). Their existence is itself redistribution, and it does not license us.

### 5.3 LIDC

- **Licence:** CC BY 3.0. Required: the data citation (DOI 10.7937/K9/TCIA.2015.LO9QL9SX), the Armato et al. paper, and the acknowledgement "The authors acknowledge the National Cancer Institute and the Foundation for the National Institutes of Health…" (https://www.cancerimagingarchive.net/collection/lidc-idri/ [page]).
- **TCIA policy:** CC BY 3.0/4.0 allows "commercial, scientific and educational purposes". Users must not "identify or contact individual participants". Mirroring public data as a service falls under the DAC rules (https://www.cancerimagingarchive.net/data-usage-policies-and-restrictions/ [page]).
- **Conclusion:** copying LIDC crops to any Vast host is allowed. Publishing derived polygons is allowed with attribution.

---

## 6. The precision probe, redesigned

### 6.1 Why the RESULTS §8 version cannot work

- **"≤ 2 px boundary error" is looser than the readers.**
  - The readers' ASSD is 0.93 px and HD95 2.0 px (S2).
  - A shape-blind disc at the right centroid scores ASSD 0.61 px (S3).
  - Errors that keep boundaries within 2 px (σ = 2 px vertex noise, 2 px shift) put LIDC IoU at 0.61–0.63 and halve legality (S4).
- **On centred frames, the probe measures the prior.** The object is at the centre, and S tells its size (S5).
- **On real data the probe would measure rater noise.** That is why it must be synthetic, with exact ground truth.
- **Synthetic shapes must be as irregular as real nodules,** or blobs pass. My first generator's disc baseline reached IoU 0.90, against 0.806 on real data (S3).

### 6.2 Design

**Purpose.** Before any LIDC SFT decision, measure whether a model plus our output format can place boundaries precisely enough, free of rater ambiguity. If it cannot, judge it on RIGA only (`phase2_entry(..., datasets=("riga",))`), or give it a synthetic pre-SFT stage first.

**Generator (Mac, numpy + polycore; deterministic seed):**
- **Ground truth.** A polygon on the 0–1000 grid, rasterised by polycore, so ground truth = target with no tie error. Its certified target is built from that raster with `build_target`, so the probe answers are exactly what SFT teaches.
- **Sizes.** Equivalent diameter drawn from the LIDC train distribution (p10 / p50 / p90 = 6.6 / 10.9 / 22.4 px at frame resolution). Frame S ∈ {50, 100, 125} by the `frame` rule, then **window jitter** (object anywhere inside, with at least 2 px clearance) and scale jitter.
- **Shape families,** calibrated so that the equal-area disc at the true centroid scores median IoU 0.80 ± 0.03, as on real LIDC (S3):
  - radial Fourier with spiculation (harmonics 2–12, amplitude scaled until calibrated);
  - lobulated unions of 2–4 discs;
  - 5 % with a hole (cavitation);
  - 5 % with two parts.

  This roughly matches LIDC's 2–5 % non-simple targets (A1).
- **Rendering.** Real LIDC train lung-window background patches, from train patients only and away from any outlined nodule; fall back to a noise texture. Interior HU: solid (−100 to +100) or ground-glass (−700 to −300). Partial-volume blur σ 0.5–1 px, CT noise σ 10–30 HU. Then `frame`'s lung window and the 512-px upscale, as for real items.

**Protocol per candidate model:**
1. **Zero-shot transfer:** the main LIDC-SFT checkpoint on 1,000 held-out probe items.
2. **Capacity:** a short probe-only SFT (20k generated items, same recipe and steps budget as one LIDC epoch, about 0.3 GPU-h on a 4090), then the same 1,000 items.

Greedy decoding plus K = 8 samples at temperature 1.

**Metrics per object (frame pixels),** stratified by diameter bins 0–6, 6–10, 10–16, 16–25, > 25:
- legal and parse rates;
- centroid error;
- relative equivalent-diameter error;
- ASSD;
- HD95;
- NSD@1 px;
- IoU;
- topology correct (holes and parts).

**Pass rule (capacity run, 6–16 px bins), with what it rests on:**

| Criterion | Value | Justification |
|---|---|---|
| legal | ≥ 95 % | S4: σ = 0.5 px vertex noise → 95.7 % legal |
| median ASSD | ≤ 0.6 px | S4: σ ≈ 1 px → 0.61. Below the reader-vs-reader 0.93 because the ground truth has no rater noise |
| median IoU | ≥ 0.80 | S4: σ = 1 px → 0.79; centroid-disc baseline on calibrated shapes 0.80 |
| median centroid error | ≤ 0.5 px | S4: a 1 px shift costs 23 IoU points |
| topology correct | ≥ 90 % on hole and two-part items | A0′ claim |

Readings:
- **Fail** → drop the model from LIDC, or try one synthetic pre-SFT stage (H-A4-6), then re-probe.
- **Pass on capacity but fail zero-shot** → the model can be precise and LIDC SFT is the bottleneck. Report it, but do not drop the model.

**Literature support:**
- NSD tolerance from inter-observer variation (Nikolov, https://arxiv.org/abs/1809.04430).
- Metrics Reloaded: NSD plus HD95/MASD for boundaries, and Boundary IoU is uninformative for small objects (https://arxiv.org/abs/2206.01653; https://arxiv.org/abs/2103.16562).
- VLM low-level geometric perception is weak and improves with targeted synthetic training (Euclid https://arxiv.org/html/2412.08737; PointArena https://arxiv.org/abs/2505.09990).
- Sensitivity to object size (https://arxiv.org/abs/2502.17422).
- Visual token stride: A1 §1.3, SmolVLM 64 px per token, i.e. 6.25 frame px at S = 50. The probe tests whether the model localises below one visual token.

---

## 7. Recommendations (ranked)

1. **Remove the centring shortcut from the gates.**
   - Add a jittered-frame val, and gate geometry on it.
   - Add the paired null-model test on centred val.
   - Train SFT with frame translation and scale jitter.
   - Evidence: S5–S6. Code: a `frame(item, offset=…, margin=…)` variant plus a target / cache rebuild for the jittered val; a Phase 2 lab item, not a Phase 1 edit.
2. **Re-specify the precision probe** as in §6: synthetic, calibrated, uncentred; legal ≥ 95 %, ASSD ≤ 0.6 px, IoU ≥ 0.80 on 6–16 px objects.
3. **Re-run MedSAM** (official min-max preprocessing) and SAM 2.1-tiny on our frames, val split, `score_mask` path. Use them plus the prior as the "near MedSAM" reference (§3.2).
4. **RIGA only on Vast Secure Cloud,** encrypted at rest, plus the depositor email, ADCIS registration and the commercial-use question to the owner (§5.2). LIDC anywhere.
5. **Augmentation:** frame jitter for both datasets; D4 for LIDC; horizontal flip for RIGA; a different op per rater row; HU window jitter; mild fundus photometric jitter. Precompute and certify offline (S1: ≈ 11 CPU-min per full view).
6. **Mixing:** T = 2 by items (RIGA ≈ 21 % of rows, ≈ 46 % of answer tokens), with label smoothing 0.1. Watch RIGA val per epoch.
7. **Add HM-IoU (K = 16), NSD@τ and HD95 to `objdel.eval`.** Report GED16 with the sample count, D_max and soft Dice. Report with the protocol named, next to our human references (S2). Never in the same column as Kohl v1/v2 numbers.
8. **nnU-Net v2 2D** (label sampling; majority vote) as the supervised specialist row, 1 fold each, about $2–6 on a 4090 (RIGA on Secure Cloud).

## 8. Rejected alternatives

| Alternative | Why not |
|---|---|
| Keep the ≤ 2 px probe | Looser than inter-reader ASSD (0.93 px); passed by a blob (0.61 px); 2 px errors give LIDC IoU 0.61–0.63 and 51 % legality (S2–S4) |
| Probe on real LIDC items | Rater noise (ASSD 0.93) dominates the error budget; cannot separate precision from ambiguity |
| Naive vertex transform for flips (with canonicalize) | Loses the IoU ≥ 0.99 certificate on 7.3 % of RIGA cups; rebuilding is cheap (S1). Acceptable only for LIDC, where it is pixel-exact, but one code path is simpler |
| Transform the upscaled LIDC image | Rotations and transposes do not commute with Pillow's two-pass bicubic (16 grey levels); transform at frame resolution |
| 90° rotations / vertical flips for RIGA | Non-physiological fundus orientation (AADG avoids geometric augmentation for disc/cup) |
| Elastic, non-90° rotation, copy-paste in SFT v1 | Masks resample and need recertification; copy-paste changes the frame's ground truth; no evidence of gain for VLM SFT |
| Rely on augmentation instead of regularisation | PaliGemma: augmentation gave no significant benefit; label smoothing and dropout fixed overfitting |
| SAM 3 / Medical SAM3 / BiomedParse / LISA-class VLM segmenters as baselines | Collapse on medical data untuned (11.9 %), LIDC contamination (BiomedParse), LIDC Dice ≤ 0.006, gated or NC licences |
| Phase 0's MedSAM numbers as the gate reference | Different crop and scoring: LIDC 0.33, against 0.78 IoU_tc on our frames (S7) and a published IoU of about 0.67 |
| Comparing our LIDC GED with published v1/v2 tables | Four incompatible protocols; even the split moves GED (CAR) |
| Boundary IoU as the probe metric | Equals mask IoU for objects under about 30 px (Boundary IoU paper) |
| RIGA on community Vast hosts | Third-party root and physical access, the readme's "copy … prohibited", Vast's own guidance to use Secure Cloud for sensitive data |
| Third-party RIGA mirrors (MRNet Drive, Zenodo 2026) | Unlicensed redistribution; our provenance must stay Deep Blue |

## 9. Hypotheses to test

Each hypothesis gives: claim; sandbox test design; pass criterion. Status: **done** = measured here; **todo** = a B-wave test.

| # | Claim | Test design | Pass criterion | Status |
|---|---|---|---|---|
| H-A4-1 | The frame leaks location and size enough to pass the RIGA geometry gates image-blind | S5: centred-disc prior sized by the train median per S, `score_mask` on every val item | prior ≥ gate on some structure | **done: PASS** (disc 0.904, cup 0.783; LIDC 0.588) |
| H-A4-2 | Frame jitter removes the shortcut | S6: same prior on frames shifted at random (union inside) | prior median IoU_tc < 0.4 on every structure | **done: PASS** (0.00 / 0.32 / 0.07) |
| H-A4-3 | A model trained without frame jitter exploits centring; with jitter it does not | Two short SFT runs (SmolVLM-256M or LFM2.5-VL-450M, LIDC, same steps, Mac MPS LoRA or one 4090 hour), with and without translation/scale jitter. Evaluate on centred and jittered val (fixed seed) | without jitter: jittered-val median IoU_tc drops ≥ 0.15 from centred. With jitter: drop ≤ 0.05, and centred-val IoU_tc ≥ the no-jitter model − 0.02 | todo |
| H-A4-4 | D4 targets rebuilt from transformed masks are 100 % certified at negligible cost | S1 on 21 RIGA + 217 LIDC val items | 100 % certified; ≤ 50 ms per target | **done: PASS** (100 %; 23.6 / 0.8 ms) |
| H-A4-5 | The full augmented train view is certified and reproducible | Precompute every train item × 8 ops (3 processes), twice | K5 = 1.0, IoU_min ≥ 0.99, identical `target_set_hash` across the two runs, wall ≤ 15 min | todo |
| H-A4-6 | A synthetic pre-SFT stage rescues a model that fails the probe | Only if some model fails §6. Probe-SFT 20k synthetic items, then LIDC SFT, vs LIDC SFT alone, same total steps | LIDC jittered-val median IoU_tc + ≥ 0.03 and legal + ≥ 3 points | todo (conditional) |
| H-A4-7 | The probe generator can be calibrated to real LIDC irregularity | Tune spiculation amplitude so the centroid-disc baseline's median IoU on synthetic shapes = real (0.806 ± 0.03). Then check S4's σ → (IoU, legal) curve on synthetic certified targets | disc baseline within ± 0.03 of 0.806; S4 curve within ± 0.05 IoU and ± 5 points legal at σ = 0.5 and 1 px | todo (S3 shows the first generator is too round: 0.90) |
| H-A4-8 | The probe discriminates precision | Feed the probe scorer synthetic "models": exact targets; σ = 0.25 / 0.5 / 1 / 2 px vertex noise; centroid disc; centred prior | exact and σ ≤ 0.5 pass; σ ≥ 1, centroid disc and prior fail | todo |
| H-A4-9 | Phase 0's MedSAM LIDC IoU 0.33 is a pipeline artefact | S7: MedSAM on our frames (official min-max vs HF port), oracle box | MedSAM median IoU_tc on LIDC ≥ 0.55 (published IoU ≈ 0.67 on consensus masks) | **done: PASS** (0.78; min-max irrelevant; Phase 0's crop / scoring was the cause) |
| H-A4-10 | SAM 2.1-tiny is at least as good as MedSAM on nodules with a box | Same as S7 with SAM 2.1-tiny (156 MB, Apache-2.0) | LIDC median IoU_tc ≥ MedSAM's − 0.02; also record RIGA | todo |
| H-A4-11 | HU window jitter does not hurt and improves robustness | LIDC SFT with and without L ± 60 / W ± 80 HU jitter; evaluate on val with the standard window and with a shifted window (L −600 / W 1500, the MedSAM2 lung window) | standard: Δ median IoU_tc ≥ −0.01; shifted window: ≥ +0.03 | todo |
| H-A4-12 | Augmentation helps RIGA more than LIDC | 2 × 2 short SFT (augmentation on/off × dataset), same steps, jittered val | RIGA Δ IoU_tc ≥ +0.02 and train–val loss gap smaller; LIDC Δ ≥ 0 | todo |
| H-A4-13 | T = 2 mixing beats proportional mixing for RIGA without hurting LIDC | Two mixed SFT runs, same steps | RIGA cup median IoU_tc (jittered val) + ≥ 0.02; LIDC within − 0.01 | todo |
| H-A4-14 | HM-IoU and NSD/HD95 implementations are correct | Unit tests (in Phase 2's eval extension): HM-IoU = brute-force assignment over all permutations (k ≤ 4, K ≤ 8, Hypothesis); NSD = 1 and HD95 = 0 on the raters' own targets; NSD matches the DeepMind surface-distance reference on 100 random pairs within 1e-6 | all equal | todo |
| H-A4-15 | nnU-Net 2D cost is as estimated | One fold, 250 epochs, LIDC label sampling, on a 4090 | ≤ 3 GPU-h; LIDC median IoU_tc on centred val ≥ the human LOO p50 − 0.10 | todo |
| H-A4-16 | The RIGA licence permits our setup | Email the depositor; owner decision on "commercial" | written yes, or Secure Cloud + encryption recorded as the policy | todo (owner) |

## 10. Sources

Primary pages read this session (WebFetch / curl):
- RIGA Deep Blue record (Wayback 2026-08-03): https://web.archive.org/web/20260803015113/https://deepblue.lib.umich.edu/data/concern/data_sets/3b591905z
- RIGA readme PDF (Wayback 2023-01-21): https://web.archive.org/web/20230121082108/https://deepblue.lib.umich.edu/data/downloads/g732d957m
- DataCite RIGA DOI metadata: https://api.datacite.org/dois/10.7302/Z23R0R29
- CC BY-NC 4.0: https://creativecommons.org/licenses/by-nc/4.0/
- MESSIDOR terms: https://www.adcis.net/en/third-party/messidor/
- TCIA LIDC-IDRI: https://www.cancerimagingarchive.net/collection/lidc-idri/
- TCIA usage policy: https://www.cancerimagingarchive.net/data-usage-policies-and-restrictions/
- Vast security FAQ: https://docs.vast.ai/guides/reference/faq/security
- Vast Secure Cloud and DPAs: https://vast.ai/article/security-and-compliance-at-vast-ai
- Vast on sensitive data: https://vast.ai/article/running-private-ai-models-without-the-risk-of-data-exposure
- Vast FAQ ("For sensitive or production workloads, we recommend checking the 'Secure Cloud' filter"): https://docs.vast.ai/guides/reference/faq/general
- LUNA16 prompt-sensitivity study (MedSAM 0.79, SAM2 0.89, MedSAM2 0.88; window [−1000, 400]): https://arxiv.org/html/2608.16832
- MedSAM2 card ("research and education purposes"): https://huggingface.co/wanglab/MedSAM2
- MedSAM official inference: https://github.com/bowang-lab/MedSAM/blob/main/MedSAM_Inference.py
- HF MedSAM preprocessor: https://huggingface.co/flaviagiammarino/medsam-vit-base/raw/main/preprocessor_config.json

Gathered by the research sub-agents (primary pages; numbers with their tables): all URLs in §2.2, §2.3, §3.1 and §4.1, plus
- MedSAM: https://pmc.ncbi.nlm.nih.gov/articles/PMC10803759/
- SAM 2: https://github.com/facebookresearch/sam2
- SAM 3: https://github.com/facebookresearch/sam3/blob/main/LICENSE, https://arxiv.org/abs/2511.16719
- BiomedParse: https://github.com/microsoft/BiomedParse, https://huggingface.co/microsoft/BiomedParse
- nnU-Net: https://github.com/MIC-DKFZ/nnUNet/blob/master/documentation/resenc_presets.md
- SimpleSeg: https://arxiv.org/abs/2601.19228
- VectorLLM: https://arxiv.org/html/2507.04664
- Molmo: https://arxiv.org/html/2409.17146
- Data Mixing Laws: https://arxiv.org/abs/2403.16952

Local evidence: `experiments/phase2/lab/A4/` (`aug_cost.py`, `rater_ceiling.py`, `probe_proto.py`, `noise_model.py`, `null_prior.py`, `null_prior_jitter.py`, `medsam_check.py` and their JSON outputs).
