# P_VIS: VLMs that output geometry, multi-rater medical segmentation, topology, test-time selection

## Skills used

- huggingface-skills:huggingface-papers: gave the fetch recipe (paper .md page, papers API). The .md page worked for 17 papers. For 12 papers the HF page was missing (404 or empty). I used arXiv HTML for those.
- nvidia-skills:nvidia-skill-finder: loaded as required. No NVIDIA product task here, so I ran no catalog lookup and installed nothing.
- engineering:debug: used as a method to isolate the failed HF fetches (reproduce, isolate, fall back to arXiv HTML).
- engineering:system-design, engineering:architecture, engineering:documentation: loaded. I used them only for the entry layout and to keep claims tied to evidence.
- engineering:code-review: loaded. No code diff to review, so no contribution.
- huggingface-skills:huggingface-best: loaded. No model-leaderboard question here, so no contribution.

## Counts

Papers read in full: 20. Abstract + tables (abstract and the numbers I could find): 9. New 2026 additions: 5 (total 34 entries).

"Read: full text" means I read the abstract, introduction, method and result prose of the converted full text. In some papers the main result tables are images. I then say so in the entry. I did not read every appendix line of the long papers.

Three IDs in the list are not on topic. 2504.09975 is OctGPT (3D shape generation). 2501.04001 is Sa2VA. 2505.22914 is cadrille (CAD code). I list them anyway, with a short note.

## Quick picture for the project

- Closest prior art stays SimpleSeg (2601.19228). Newer neighbour: Moondream Segmentation (2604.02593), which also uses RL on a rasterized path.
- Nobody found uses a VLM polygon with LIDC multi-rater scoring. The multi-rater papers use U-Net, diffusion or VAE models.
- Most useful for our main error (small off-centre nodules): 2504.05977 (random-crop LIDC), 2603.28120 (reward curriculum for sparse IoU), 2608.03471 (small-box coordinate sensitivity), 2604.27720 (zero-shot VLM boxes are near chance).

---

## A. VLMs that output geometry (polygons, points, boxes, paths)

### 2601.19228 — Towards Pixel-Level VLM Perception via Simple Points Prediction (Jan 2026)
- Read: full text
- Core idea: SimpleSeg makes a VLM write a clockwise sequence of contour points as plain text. It trains with SFT, then RL with an IoU-based reward.
- Key result: RL adds +9.7, +9.8 and +10.5 gIoU over SFT on refCOCO, refCOCO+ and refCOCOg (65.5, 60.8, 60.4 to 75.2, 70.6, 70.9). Best point density is about 221 tokens (35.6 cIoU at 78 tokens, 72.5 cIoU at 859 tokens). A hard length penalty in the reward drops refCOCO gIoU from 76.9 to 66.7. Unordered points do not form a valid polygon. Models: Qwen2.5-VL-7B and Kimi-VL.
- Use for us: copy the clockwise start-point rule and run an epsilon (point density) sweep on our nodule polygons. Do not add a hard length penalty.

### 2604.02593 — Moondream Segmentation: From Words to Masks (Apr 2026) [new addition]
- Read: full text (method and appendix; I did not recover the ablation numbers)
- Core idea: A 2B-active VLM decodes an SVG-style vector path. A rasterizer makes a coarse mask. RL optimizes the mask after rasterization. A separate refiner finishes the mask.
- Key result: cIoU 80.2 % on RefCOCO val and 62.6 % mIoU on LVIS val. The RL reward is piecewise: box IoU first, then mask IoU, with a boundary-IoU band (width 5 % of image diagonal). It penalizes false negatives more than false positives. Invalid traces are filtered from the group.
- Use for us: test a false-negative-heavy IoU term plus boundary IoU in our reward, and filter invalid polygons out of the GRPO group instead of giving them a score.

### 2507.04664 — VectorLLM: Human-like Extraction of Structured Building Contours via Multimodal LLMs (Jul 2025)
- Read: full text
- Core idea: An MLLM regresses building corner points one by one. Training is pretrain, SFT, then DPO.
- Key result: +5.6 AP on WHU, +7.1 AP on WHU-Mix, +13.6 AP on CrowdAI over the prior best. With oracle boxes, WHU AP is 90.2. DPO adds 1.8 AP and 1.4 AP75. DPO pairs come from SFT outputs with IoU below 0.8, plus labels with random point deletion or insertion.
- Use for us: a cheap DPO round on our model's own bad polygons (IoU below 0.8, or corrupted ground truth) is a fallback if GRPO stalls.

### 2608.03471 — Hi-Token: Hierarchical Coordinate Tokenization for Generative Visual Grounding (Aug 2026)
- Read: full text
- Core idea: Each token encodes axis, decimal place and digit value (12 tokens per box, 60 types). The GRPO reward adds multi-scale coordinate accuracy, strict-IoU bonuses and a validity gate for near-zero overlap.
- Key result: gains on three backbones in matched comparisons. The paper reports that Hi-Token beats an extra-tuned flat baseline on RefCOCO before RL. On small boxes Hi-R1 gets mIoU 58.0 and P@0.9 14.7 (375 small test boxes). A one-bin coordinate shift hurts small boxes more than large ones. I did not recover the exact backbone gain table.
- Use for us: small nodules are sensitive to one grid step on the 0-1000 grid. Test finer grid digits or a coarse-then-fine coordinate order, and a validity gate in the reward.

### 2603.28069 — MolmoPoint: Better Pointing for VLMs with Grounding Tokens (Mar 2026)
- Read: abstract + tables (abstract, method overview)
- Core idea: Instead of writing coordinates as text, the model emits tokens that select a coarse patch, then a subpatch, then a location. Points are produced in a fixed order with a "no more points" class.
- Key result: 70.7 % on PointBench, 61.1 % on ScreenSpot-Pro, +6.3 % on Molmo2Track. The paper also claims higher sample efficiency (no number captured).
- Use for us: a "no more vertices" token and a fixed vertex order are cheap to copy. The patch-selection head needs architecture changes, which our text-only LoRA plan avoids.

### 2510.03230 — Improving GUI Grounding with Explicit Position-to-Coordinate Mapping (Oct 2025)
- Read: abstract + tables (main tables not recovered)
- Core idea: Ruler tokens act as gridlines, so the model adjusts a coordinate instead of generating it from scratch. Interleaved MRoPE gives width and height equal weight.
- Key result: consistent gains on ScreenSpot, V2 and Pro, largest on high-resolution screens. No headline number captured. The paper says it does not reach state of the art.
- Use for us: a grid-hint idea for the 0-1000 grid (print a ruler row in the prompt). Low priority, needs a test on our crops.

### 2601.13633 — EGM: Efficient Visual Grounding Language Models (Jan 2026)
- Read: full text
- Core idea: Small VLMs fail at grounding because of language understanding, not vision. EGM adds reasoning traces (SFT) then GRPO. Reward is 0.5 x IoU + 0.5 x grounding success.
- Key result: Qwen3-VL-2B-Thinking goes from 83.6 to 85.6 (SFT) to 89.6 (RL) on the RefCOCO family average; RL adds 4.0 for 2B versus 2.6 for 8B. EGM-8B reaches 91.4 IoU at 737 ms versus 90.5 IoU at 4,320 ms for Qwen3-VL-235B. Removing reasoning costs only 0.9 for 2B (88.7 vs 89.6).
- Use for us: same 2B backbone as ours. RL gains are largest at 2B, which supports our RL path. Reasoning text helps little at 2B, so skip it to save tokens.

### 2602.09701 — RL-Driven Vision–Language Grounding for Fine-Grained Referring Segmentation (Feb 2026)
- Read: full text (method and results prose)
- Core idea: GenSeg-R1 fine-tunes Qwen3-VL (4B, 8B) with GRPO to emit a box plus two interior points. A frozen SAM 2 makes the mask. The reward is the SAM 2 mask IoU.
- Key result: 0.7127 cIoU on RefCOCOg val (+15.3 over the Qwen3-VL Instruct baseline, +3.3 over Seg-Zero-7B). On GRefCOCO, 82.40 % accuracy on no-target prompts with 0.5 % false-negative rate. Qwen3-VL uses a [0,1000] coordinate grid.
- Use for us: shows the no-target output ("empty polygon") can be trained by RL. Our nodule crops always contain a nodule, but the RIGA cup can be absent, so test an explicit empty answer.

### 2605.12953 — Seg-Agent: Test-Time Multimodal Reasoning for Training-Free Language-Guided Segmentation (May 2026)
- Read: full text (method, ablation prose)
- Core idea: No training. Candidate regions are drawn on the image (Set-of-Mark), and the MLLM selects and refines in a loop: generate, select, refine.
- Key result: +5.1 gIoU over the single-step baseline on reasoning segmentation. It matches the trained method Seg-Zero (Qwen2.5-VL-7B + SAM2-L). I did not capture the absolute numbers.
- Use for us: a render-and-reselect pass is a test-time option. Draw the model's polygon on the crop and ask a second pass to accept or fix it.

### 2606.09303 — Reason Twice: Segmentation via Candidate Discovery and Comparative Reasoning (Jun 2026)
- Read: full text (method and appendix; the main result tables were not recovered)
- Core idea: Rea2Seg first makes candidate masks from the MLLM attention maps. Then an MLLM scores the candidates and the best one is chosen (rerank).
- Key result: the paper reports gains on its new ReasonSeg-SGDR benchmark and on ReasonSeg. I did not recover the numbers.
- Use for us: a test-time selector over K sampled polygons. Our checker can score each candidate exactly, so use it as the selector on the dev split to find the upper bound of best-of-K.

### 2603.14382 — StAR: Segment Anything Reasoner (Mar 2026)
- Read: full text (tables are images; numbers not recovered)
- Core idea: A better RLVR recipe (tuning scheme, rewards, answer format), plus parallel test-time scaling: sample N answers and pick masks by IoU clustering (threshold 0.85), with the mode of the predicted instance count.
- Key result: StAR-7B with mask-level majority voting is "overall comparable" to the SAM 3 Agent (72B) on the ReasonSeg-X test set (paper text). Trained with 5k samples. ReasonSeg-X has 1,169 samples (240/156/773 train/val/test).
- Use for us: our best test-time candidate. Sample K polygons, cluster by IoU, return the medoid of the largest cluster. No new training.

### 2503.01785 — Visual-RFT: Visual Reinforcement Fine-Tuning (Mar 2025)
- Read: abstract + tables (abstract, reward formula)
- Core idea: GRPO with verifiable perception rewards. Detection reward is IoU + confidence + format.
- Key result: +24.3 % one-shot fine-grained classification (about 100 samples); +21.9 on COCO two-shot; +15.4 on LVIS.
- Use for us: early proof that an IoU reward works with few samples. Our reward is already stronger.

### 2504.07615 — VLM-R1: A Stable and Generalizable R1-style Large Vision-Language Model (Apr 2025)
- Read: full text (abstract, intro, ablation highlights)
- Core idea: Open RL framework for VLMs with rule-based rewards. It reports reward hacking in detection and a fix by reward engineering.
- Key result: Qwen2.5-VL-3B RL scores 63.16 on out-of-domain LISA-Grounding (SFT 54.82). On OVDEval it reaches 31.01 nms-AP (SFT 26.50).
- Use for us: RL transfers better than SFT to new data, which fits our off-centre shift. Watch for reward hacking (for example many tiny polygons); log failure cases.

### 2508.14153 — LENS: Learning to Segment Anything with Unified Reinforced Reasoning (Aug 2025)
- Read: abstract + tables (abstract only)
- Core idea: RL with sentence-, box- and mask-level rewards, so reasoning text and the mask are tuned together.
- Key result: 81.2 % average cIoU on RefCOCO/+/g with Qwen2.5-VL-3B, up to +5.6 over GLaMM.
- Use for us: shows a 3B VLM can match larger models with multi-level rewards. We already have box-free reward; no new action.

### 2501.04001 — Sa2VA: Marrying SAM2 with MLLM for Dense Grounded Understanding of Images and Videos (Jan 2025)
- Read: abstract + tables
- Core idea: An MLLM emits a [SEG] token that prompts SAM-2 for masks, for images and video.
- Key result: Sa2VA-8B gets cIoU 81.9 / 76.5 / 78.9 on RefCOCO / + / g, +2.4, +3.9, +4.7 over GLaMM-7B.
- Use for us: the decoder route (token to SAM) is the alternative to our text polygons. Not needed unless polygons fail.

### 2511.16719 — SAM 3: Segment Anything with Concepts (Nov 2025)
- Read: abstract + tables
- Core idea: Concept-prompted detect, segment and track model. Its data engine uses AI verifiers (fine-tuned Llama 3.2) to accept or reject masks.
- Key result: the paper claims it doubles prior accuracy on concept segmentation. On SA-Co/VEval SA-V it scores cgF1 30.3 against a human 53.1.
- Use for us: the verifier idea matches test-time selection. A generic SAM does not know nodule rater style.

### 2511.21631 — Qwen3-VL Technical Report (Nov 2025)
- Read: abstract + tables (grounding and coordinate sections)
- Core idea: Qwen3-VL family (dense 2B, 4B, 8B, 32B; MoE variants), 256K context. Coordinates use a normalized [0,1000] range.
- Key result: no single grounding number used here. The report states the [0,1000] grid improves robustness to image size and aspect ratio.
- Use for us: confirms our 0-1000 grid is the native format of Qwen3-VL-2B.

### 2504.09975 — OctGPT: Octree-based Multiscale Autoregressive Models for 3D Shape Generation (Apr 2025)
- Read: abstract + tables
- Core idea: Multiscale octree tokens for 3D shape generation.
- Key result: training 13x faster and generation 69x faster than prior autoregressive 3D models.
- Use for us: off topic. Coarse-to-fine tokens echo Hi-Token, nothing to run.

### 2505.22914 — cadrille: Multi-modal CAD Reconstruction with Online Reinforcement Learning (May 2025)
- Read: full text (method, reward, metrics)
- Core idea: A VLM writes CAD code. SFT on procedural data, then online GRPO with a programmatic reward.
- Key result: SFT-only CAD-Recode has CC3D IoU about 60 % and 10 % invalid outputs. The cadrille reward is 10 x IoU minus 10 for an invalid program. Training uses only hard examples (mean SFT reward below 7.5 over 3 samples). Online GRPO beats offline DPO.
- Use for us: a model for our reward: scale IoU, add a big negative for invalid polygons, and mine hard crops (small off-centre nodules) for RL.

---

## B. Medical VLM grounding and agents

### 2508.08177 — MedReasoner: Reinforcement Learning Drives Reasoning Grounding from Clinical Thought to Pixel-Level Precision (Aug 2025)
- Read: full text
- Core idea: An MLLM reasoner outputs box plus points with RL. A frozen MedSAM2 makes the mask. Task: implicit clinical queries (U-MRG-14K, 14K samples, 108 categories).
- Key result: overall IoU 32.42 with Lingshu-7B and the soft reward (IoU and pDice only); the strict reward scores 31.69, and a hard-threshold baseline scores 15.85. With reasoning, 32.42 vs 30.29 without. Seg-Zero-7B gets IoU 16.14. All models are weak on histology (11.66).
- Use for us: shows that soft (continuous) rewards beat strict ones in medical grounding. Our checker is already continuous.

### 2603.05911 — CORE-Seg: Reasoning-Driven Segmentation for Complex Lesions via Reinforcement Learning (Mar 2026)
- Read: full text (method and ablation prose)
- Core idea: End-to-end MLLM + prompt adapter + mask head. SFT on chain-of-thought, then GRPO. The mask reward is dual-granularity and penalizes disconnected regions.
- Key result: 37.06 % mean Dice and 18.42 % failure rate on ComLesion-14K, +14.89 points over baselines. Random rationales drop Dice to 33.85 % and raise failures to 23.64 %.
- Use for us: a connectivity penalty already sits in a medical RL reward. Our Betti reward goes further. The failure rate (no valid mask) is a metric worth logging.

### 2602.03320 — MedSAM-Agent: Empowering Interactive Medical Image Segmentation with Multi-turn Agentic Reinforcement Learning (Feb 2026)
- Read: full text (method, ablation prose)
- Core idea: A multi-turn agent calls a SAM-type tool. It trains with an SFT cold start, then RL with a process reward that punishes extra turns.
- Key result: 6 modalities, 21 datasets. In the ablation with MedSAM2, Dice is 0.794 (0.772 without the improvement reward) and about 2.11 turns. Direct RL without SFT fails on CT and MRI.
- Use for us: RL alone struggles on low-contrast CT. Our SFT-then-RL order is right. A second "fix" turn is a later option.

### 2603.28120 — MedLoc-R1: Performance-Aware Curriculum Reward Scheduling for GRPO-Based Medical Visual Grounding (Mar 2026) [new addition]
- Read: full text (method; table numbers not recovered)
- Core idea: Fixed IoU thresholds give sparse reward on small or vague medical regions. MedLoc-R1 tracks recent reward in a sliding window and tightens the IoU threshold only when the policy is ready.
- Key result: improves A@0.5 and A@0.8 and training stability over GRPO baselines on three medical grounding benchmarks, over 3 seeds. I did not capture the numbers.
- Use for us: direct fix for off-centre small nodules. Start with a loose reward and tighten as the group IoU rises.

### 2605.26621 — MedVol-R1: Reward-Driven Evidence Grounding for Volumetric Reasoning Segmentation (May 2026) [new addition]
- Read: abstract + tables (abstract and result lines)
- Core idea: A VLM picks a key slice and a 2D box. A frozen MedSAM2 propagates it to 3D. GRPO reward covers evidence selection, 2D grounding and cross-slice coherence.
- Key result: GRPO adds +4.34 DSC (AbdomenCT-1K), +2.20 (CT-ORG), +9.25 (KiTS23) over SFT. Ablation: SFT-only 85.52 to full model 89.86 DSC.
- Use for us: the largest RL gain is on the hardest set (KiTS23, tumours), a hint that RL helps small lesions most.

### 2604.27720 — Why Does Grounding Hurt Medical VQA? Benchmarking, Diagnosis, and Fine-Tuning of Vision-Language Models (Apr 2026)
- Read: full text (abstract, sections IV-C and IV-D)
- Core idea: Tests whether medical VLMs can localize. Zero-shot boxes are bad, and cropping to even a perfect box lowers VQA accuracy.
- Key result: all off-the-shelf systems have mean IoU 0.05 to 0.24 on SLAKE grounding (center-box baseline 0.10). Answer-only SFT of Qwen-2.5-VL-7B gives 0/418 parseable boxes. Adding a 15 % mix of real box supervision gives IoU 0.390 (0.360 on unseen images).
- Use for us: confirms that off-the-shelf VLMs cannot do our task and that localization supervision must stay in the mix. Cropping hides context, which links to our centred-crop training gap.

---

## C. Multi-rater and ambiguous medical segmentation (LIDC)

### 2504.05977 — Diffusion Based Ambiguous Image Segmentation (Apr 2025)
- Read: full text
- Core idea: Studies diffusion design choices for LIDC-IDRI. x- and v-prediction beat epsilon-prediction. Harder noise via input scaling helps. It adds a random 64x64 crop version of LIDC.
- Key result: GED (16 samples) 0.233 on 128x128 centred crops; on random 64x64 crops, GED 0.312 with 16 predictions (the best prior was 0.320) and dice 0.442 with 16 predictions (0.459 for 4 predictions with post-processing). The paper says centred LIDC is "too easy", because the lesion location is known. Typical failure: predictions that are too small. A relative-area post-process (threshold 0.5) helps.
- Use for us: this is the same effect as our off-centre problem. Report our metrics on this random-crop protocol too. Test a relative-area filter on our polygons.

### 2512.00748 — Disentangling Preference and Ambiguity for Controllable Medical Lesion Segmentation (Dec 2025)
- Read: full text
- Core idea: ProSeg uses two latent variables, one for rater preference and one for boundary ambiguity. It is trained by variational inference, so it can make diverse and rater-specific masks.
- Key result: LIDC-IDRI GED 0.1152, soft Dice 91.53 %, mean personalised Dice 90.25 % (per-rater 89.49 to 91.01). NPC GED 0.2272. The four LIDC "raters" are virtual: 12 real readers are assigned to 4 slots by area rank.
- Use for us: LIDC rater identity is not stable, so a rater-conditioned output is weak. Scoring against the SET of raters (our checker) is the safer choice. Use ProSeg's GED and Dice-max as extra metrics.

### 2603.04024 — Volumetric Directional Diffusion: Anchoring Uncertainty Quantification in Anatomical Consensus for Ambiguous Medical Image Segmentation (Mar 2026)
- Read: full text
- Core idea: A coarse consensus prediction anchors the diffusion. The diffusion only explores boundary changes (residuals). This keeps volumetric topology.
- Key result: 3D LIDC-IDRI MaxDice 0.8031, GED 0.2081 (3D nnUNet MaxDice 0.7622). Without the anchor, MaxDice (16 samples) drops to 0.6670 and GED rises to 0.5809. Results hold from 2 to 4 raters (MaxDice about 0.803).
- Use for us: a consensus anchor is the right frame for set-of-raters scoring. Our model could predict a consensus polygon first, then variants.

### 2605.08210 — Harmonized Feature Conditioning and Frequency-Prompt Personalization for Multi-Rater Medical Segmentation (May 2026)
- Read: full text
- Core idea: Probabilistic U-Net with a noise harmonizer and rater-style prompts in the frequency domain, plus a GED alignment loss.
- Key result: LIDC GED 0.105 (D-Persona 0.136); NPC-170 GED 0.175 (0.197). Personalised LIDC Dice-max 92.65 %, Dice-match 90.0 %, +1.61 mean Dice over D-Persona. Metrics given for K = 10 to 50 samples.
- Use for us: the SOTA table for LIDC multi-rater. Compare our Dice-max and GED against these numbers, with the caveat that they use 128x128 centred crops.

### 2605.02437 — Multi-Rater Calibrated Segmentation Models (May 2026)
- Read: full text
- Core idea: Treat voxel agreement among K raters as an ordinal level and train with a ranked probability score loss plus a binary loss.
- Key result: lowest multi-rater ECE on all four sets (CoCaHis, REFUGE disc, REFUGE cup, LongCIU). REFUGE cup 6.67 vs 12.01 for random sampling; disc 1.09 vs 3.41 (values x 10^-2). AUC is unchanged (cup 92.57 vs 92.06).
- Use for us: uses REFUGE, the same family as RIGA. If we output a confidence per polygon, the ordinal-agreement target is the calibration baseline.

### 2603.27197 — K-alpha-LOS finds Consensus: A Meta-Algorithm for Evaluating Inter-Annotator Agreement in Complex Vision Tasks (Mar 2026) [new addition]
- Read: full text (LIDC section; abstract)
- Core idea: Match annotations in space first, then compute a reliability score (Krippendorff alpha). It works from boxes to 3D segmentation.
- Key result: on LIDC-IDRI (nodules at least 3 mm, four reading sessions treated as unnamed raters) alpha is stable at about 0.37 for localization thresholds 0.2 to 0.5, then drops to about 0.33 at 0.6. Raters agree where the nodule is but not how far it extends.
- Use for us: gives the evidence that location is consistent and extent is not. This supports a miss (wrong location) being a model error, not rater noise.

---

## D. Topology-aware learning

### 2603.13054 — Topo-R1: Detecting Topological Anomalies via Vision-Language Models (Mar 2026)
- Read: full text (method and appendix; main results table not recovered)
- Core idea: Tests VLMs on finding broken or spurious connections and missing or extra branches in masks. Training: SFT, then GRPO with a reward for localization, classification and skeleton structure. Matching uses per-type Hungarian assignment on IoU, with a soft tiered IoU score.
- Key result: leading VLMs score near random. Topo-R1 (Qwen3-VL-4B) reaches F1@0.5 of 45.2 on the test set. With only 25 % of the SFT data it keeps over 91 % of the F1. The reward weights are 0.10 format, 0.85 accuracy, 0.05 topology.
- Use for us: the closest recipe for a Betti-based reward in GRPO. Verifiable Betti labels are the right target. Note the small topology weight (0.05).

### 2603.19169 — ARIADNE: A Perception-Reasoning Synergy Framework for Trustworthy Coronary Angiography Analysis (Mar 2026) [new addition]
- Read: abstract + tables (abstract, introduction)
- Core idea: DPO on Sa2VA, using Betti numbers to pick the preferred mask. A second RL stage lets the model reject ambiguous candidates.
- Key result: centerline Dice 0.838 on 1,400 angiograms, 41 % fewer false positives than geometric baselines. External checks on ARCADE and XCAD.
- Use for us: Betti numbers work as preference labels (chosen vs rejected). This is a cheaper alternative to a Betti reward in GRPO.

---

## E. SAM-style medical foundation models

### 2504.03600 — MedSAM2: Segment Anything in 3D Medical Images and Videos (Apr 2025)
- Read: abstract + tables
- Core idea: SAM 2 fine-tuned on over 455,000 3D image-mask pairs and 76,000 frames.
- Key result: the human-in-the-loop pipeline cuts manual annotation cost by more than 85 % (5,000 CT lesions, 3,984 liver MRI lesions, 251,550 echo frames).
- Use for us: MedSAM2 is the segmenter used by MedReasoner and MedVol-R1. We could use it as an optional box-to-mask baseline on the same crops.
