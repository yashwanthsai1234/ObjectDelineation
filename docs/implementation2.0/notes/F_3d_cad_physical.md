# Where the certified geometry verifier goes next: 3D voxels, CAD, Physical AI

## Skills used

| Skill | What it contributed |
|---|---|
| `nvidia-skills:nvidia-skill-finder` (run twice) | Live catalog check (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`), plus a read of the skill cards. Relevant, none installed, recommendation only: `omniverse-cad-to-simready` (CAD → SimReady orchestrator: conversion, physics, conformance, validation), `nv-segment-ct` / `nv-segment-ctmr` (VISTA3D-style label maps + "label-map evidence"), `nv-generate-ct-rflow` (synthetic CT + masks, ≈16 GB at 256³, "not safe as training data without independent quality review"), `dicom-series-to-volume` (DICOM → HU NIfTI with affine evidence), `physical-ai-physics-selection` (Newton vs PhysX vs MuJoCo), `physical-ai-neural-reconstruction` (NuRec/3DGS router; "do NOT use for SimReady"), `cdf-*` (Cosmos data factory; its simulation card is PCB/solder-defect specific, so it is a weak fit), `paidf-cosmos-predict`, `i4h-workflow-train-rl` (Isaac-for-Healthcare RL policies), `warp-eval` (evidence-gated Warp candidate check). These shaped Sections 1.2, 3.3 and Direction C. |
| `engineering:system-design` (run twice) | Requirements → components → trade-offs structure of Section 5 (verifier sketch). |
| `system-design-skills:back-of-the-envelope` | Throughput chain in Section 5.2 / 5.4: 2,048 solids × 128³ winding-number queries at 10⁹/s ≈ 4 s per GRPO step; 64³ ≈ 0.5 s; OCC execution stays the CPU bottleneck. All labelled **[est]**. |
| `engineering:architecture` | ADR in Section 5.5: Warp + fVDB + own CUDA kernels vs pure own-CUDA vs CPU OCC/trimesh. |
| `engineering:code-review` | Review of `src/objdel/checker/volume.py` as the CPU oracle for a GPU 3D verifier (Section 4.1): the Euler and 26/6 logic is sound; risks are z-gaps and slice anisotropy, not arithmetic. |
| `engineering:debug` | Failure-mode table for a certified reward that is silently wrong (Section 5.6): OCC timeouts laundered into GRPO advantages, voxel aliasing of thin walls, non-watertight inputs, frame misalignment. |
| `engineering:documentation` | Report layout: bottom line first, then evidence, then directions, then open questions; links instead of copied text. |
| `huggingface-skills:huggingface-papers` | HF paper API: cadrille (2505.22914) → GitHub col14m/cadrille, models `maksimko123/cadrille`, `maksimko123/cadrille-rl`. CAD-Recode (2412.14042) → models `filapro/cad-recode(-v1.5)` and the `filapro/cad-recode-v1.5` dataset. ConCAD (2609.34494) is not indexed on HF (no linked code or models). The papers search surfaced BenchCAD, GACO-CAD, VVR, RFW-A, SIMART and EmbodiedGen V2. |
| `huggingface-skills:huggingface-best` | Official HF benchmark list (58 datasets): **no official leaderboard covers CAD code generation, 3D medical segmentation or 3D shape**. Closest are `google/physical-engineerings-last-exam`, `VLABench`, `tiiuae/PBench`. So model choice falls back to paper results and model metadata (Section 5.7). Key metadata: cadrille 2.21 B and CAD-Recode-v1.5 1.54 B are both **CC-BY-NC-4.0 (non-commercial)**; Qwen3-VL-2B-Instruct 2.13 B is Apache-2.0; Cosmos-Reason2-2B 2.44 B uses NVIDIA's "other" license, needs ≥24 GB, 584k downloads. |
| Web (Exa search/fetch, WebFetch) | Every external number in the report; URLs inline. |

Survey date: 2026-11-08. Scope: web survey (primary sources: arXiv, CVF/MICCAI/ICLR proceedings, vendor docs, GitHub)
plus a reading of this repo (README "3D volumes", `src/objdel/checker/volume.py`, ROADMAP Phase 7,
`docs/research/reports/GPU delineation use case map.md`, `experiments/phase0/RESULTS.md`, `experiments/phase0/results/lidc/clusters.csv`).
No repo edits, no spend. Numbers marked **[est]** are my own back-of-envelope figures, not measurements.

---

## 0. Bottom line (read this first)

1. **The strongest research direction is "certified geometric rewards": one GPU verifier for 2D rings, 3D voxel
   volumes and closed triangle-mesh solids (tessellated CAD), used as an RL reward.** Every CAD-RL paper found
   in 2025–2026 scores rollouts with **CPU, approximate** geometry: cadrille uses trimesh mesh booleans per
   component pair; ConCAD estimates volumetric IoU by **Monte Carlo with 8,192 points**; CAD-Coder uses Chamfer distance;
   a Stanford CS224r project uses Open CASCADE (OCC) Betti-1 in a multiprocessing pool. None publishes a certified, batched,
   GPU-resident reward with topology terms (b0/b1/b2, watertightness, manifoldness, self-intersection). The
   repo already has the 2D half (bit-identical CUDA checker, RL that moved held-out IoU 0.511 → 0.659) and the
   3D Euler-characteristic method on CPU.
2. **For "CAD": aim at image/scan → CAD (reverse engineering), not CAD → 3D.** CAD → 3D (tessellate or voxelize
   a B-rep) is solved engineering: OCC tessellation, the Omniverse CAD converters, and NVIDIA's `omniverse-cad-to-simready`
   workflow skill. Image/scan → CAD is the open frontier, and it is RL-driven: cadrille (ICLR 2026 Oral) reaches
   92.2 % IoU and 0 % invalid on DeepCAD multi-view, and RL is what brought invalidity to 0.0–0.2 %.
   It also matches the repo's stack: a 2B VLM writes a program, and a verifier scores it.
3. **For "Physical AI": the exact-geometry entry point is assets and robot-facing geometry, not world-model video.**
   Video physics evals (Physics-IQ, VideoPhy-2, WorldScore, PDI-Bench) all reach geometry through *learned* lifting
   (SLAM, monocular depth, trackers). That lifting step loses exactness, so an "exact verifier for world
   models" is mostly hype today. Explicit-geometry gaps are real. SimReady's own requirement index lists
   "Transparent physical bodies should be watertight" and "Meshes should not overlap unnecessarily" with **no
   validator**. SimReady's FAQ says it has **no scene-level validators**. Collision-proxy quality changes robot outcomes:
   CoACD vs V-HACD raised a drawer-opening RL success rate from **49 % to 80 %**.
4. **Honest caveats.** (a) The speed pitch is weak for CAD. The CPU bottleneck is *executing the CAD program in OCC*,
   not scoring it, and one vendor claims CPU B-rep verification in <10 ms. The pitch must be **certification,
   topology, set-valued scoring and batching**, not raw speed. This is the same lesson as in the repo's use-case map.
   (b) In LIDC, 3D cavities are rare: **92 of 2,654 nodule clusters raw, 29 after the 5-px filter** have a reader with a
   cavity (`clusters.csv`). So a b2 reward on LIDC trains mostly *against false cavities and tunnels from slice
   stacking*. It does not teach real cavities. (c) "Exact" is exact for a stated discretization. For curved
   CAD solids, the honest output is a **certified interval** (Section 5.2), not a single "exact IoU".

Three candidate directions with milestones are in Section 6. My ranking: **A (certified rewards, 2D → 3D → CAD)
first, B (set-valued ambiguous geometry) as the scientific spine, C (verifier-driven data engine) only after A
produces the kernels**.

---

## 1. 3D: outputs, segmentation, topology, GPU libraries

### 1.1 How LLMs/VLMs emit 3D shape today

| Output form | Examples (date) | Notes for a verifier |
|---|---|---|
| Per-slice 2D outlines stacked into a volume | This repo; MedVol-R1 (MICCAI 2026) picks a key slice and 2D boxes, then MedSAM2 propagates to 3D. Its GRPO reward includes a **cross-slice consistency** Dice term. https://papers.miccai.org/miccai-2026/paper/2593_paper.pdf | Geometry is explicit and exact. Stacking creates z-artifacts (false tunnels/cavities) that a 3D Betti check catches. MedVol-R1 is the closest prior art to the repo's Phase 7. It scores with Dice, with no topology and no reader set. |
| Program (CadQuery, OpenSCAD, Three.js, voxel API) | CAD-Recode (ICCV 2025), cadrille (ICLR 2026), P3D-Bench (2026: 400 text, 400 image, 203 assemblies; scores executability, geometry, **topology**, MLLM judge, part structure) https://arxiv.org/html/2606.11152v2 ; VoxelCodeBench (Meta FAIR, 2026-04: 220 tasks, "executable code is far easier than spatially correct output") https://arxiv.org/html/2604.02580 ; 3DCodeBench (Google DeepMind, 2026-05: successful renders "still suffer from disconnected or floating 3D geometric components") https://arxiv.org/html/2606.01057v1 | A program executes to exact geometry. **"Disconnected/floating components" is a b0 failure**, and a certified b0 check catches it directly. |
| Voxel tokens via VQ-VAE | ShapeLLM-Omni (text → 64³ voxel tokens → rectified-flow mesh; described in the AR3D-R1 supplement, CVPR 2026) https://openaccess.thecvf.com/content/CVPR2026/supplemental/Tang_Are_We_Ready_CVPR_2026_supplemental.pdf ; CTTok (MICCAI 2026: Qwen-1.7B emits tokens from a 262,144-entry codebook of 16³ CT patches; 256×256×128 volumes) https://papers.miccai.org/miccai-2026/paper/2284_paper.pdf ; SuperVoxelGPT (2026-05: 12.8 % of uniform-voxel sequence length, Qwen2.5-0.5B) https://arxiv.org/html/2605.29655v2 | The decoder breaks exactness: the output is checkable, but the reward cannot be attributed token-exactly. |
| Mesh tokens | LLaMA-Mesh; Mesh-RFT (NeurIPS 2025: Boundary Edge Ratio + Topology Score, M-DPO, −24.6 % Hausdorff) https://proceedings.neurips.cc/paper_files/paper/2025/file/e4ca3675c3569a7c222c8351e4db8b9e-Paper-Conference.pdf ; Mesh-Pro (2026: ray-cast "bad face" integrity reward, async RL 3.75× faster) https://arxiv.org/html/2603.00526v1 | **These are mesh-quality rewards, but heuristic.** Ray-cast integrity and boundary-edge ratio are proxies for watertightness, not certificates. |
| Text-to-3D RL with learned rewards | AR3D-R1 / Hi-GRPO (CVPR 2026, rewards = HPS, UnifiedReward, Qwen2.5-VL) https://openaccess.thecvf.com/content/CVPR2026/html/Tang_Are_We_Ready_for_RL_in_Text-to-3D_Generation_A_Progressive_CVPR_2026_paper.html ; DreamCS (2025, learned 3D reward from 30k Llama-Mesh-scored meshes) https://arxiv.org/html/2506.09814v2 | 3D-generation RL runs almost entirely on **learned preference rewards**. Verifiable geometric rewards are the missing half. |

Also relevant as a precedent for "programmatic verifier → RL transfers": **VVR, Verifiable Visual Rewards** (2026-09-28)
generates scenes of geometric objects and relations with deterministic verifiers. RL on it raises SD3.5-Medium from 2.8 % to
28.3 % on VVRBench, with easy-to-hard and out-of-domain transfer. https://huggingface.co/papers/2609.35641

### 1.2 3D medical segmentation state of the art

- **nnU-Net is still the bar.** "nnU-Net Revisited" (MICCAI 2024) argues many claimed gains vanish under rigorous validation.
  https://arxiv.org/html/2404.09556v2 A 2026 thyroid study: nnU-Net v2 0.914–0.950 DSC beats every zero-shot foundation model,
  and MedSAM-2 zero-shot reaches up to 0.907. https://papers.miccai.org/miccai-2026-sat/paper/ASMUS_021.pdf
- **VISTA3D** (CVPR 2025, NVIDIA/MONAI): 127 classes of automatic + interactive 3D segmentation; supervoxel distillation from 2D backbones.
  https://openaccess.thecvf.com/content/CVPR2025/html/He_VISTA3D_A_Unified_Segmentation_Foundation_Model_For_3D_Medical_Imaging_CVPR_2025_paper.html
  NVIDIA's own review: zero-shot SAM2 on 3D CT is "far from satisfactory" and over-tracks objects into empty slices
  (e.g. liver 0.454 vs VISTA3D 0.959 on MSD03). https://arxiv.org/html/2408.11210
- **MedSAM2** (2025-04): fine-tuned SAM2.1-Tiny on 455,000+ 3D image–mask pairs and 76,000 video frames. Median DSC: CT organs 88.84 %,
  CT lesions 86.68 %. User study: 5,000 CT lesions, 3,984 liver MRI lesions, 251,550 echo frames; ">85 %" annotation cost cut.
  https://arxiv.org/html/2504.03600
- **Foundation-model audit** (FAU, 80+ datasets): MedSAM2 is "most foundational"; SAM-Med3D and VISTA3D excel on 3D CT but need more
  compute and memory. https://cris.fau.de/publications/365627564/
- **Language-to-volume**: VoxTell (2025-11, DKFZ): free text → 3D masks, 62K+ volumes, 1K+ classes. https://arxiv.org/html/2511.11450
  Medical VLMs that emit 3D masks (MedVL-SAM2, RepConTUR, PSMA PET/CT VLM) all use a **[SEG] token + dense decoder**, not explicit geometry.
- **3D medical shape data**: MedShapeNet, 100,000+ shapes (106,101 download URLs) from 23 datasets. MedShapeNetCore is watertight and manually
  checked. https://pmc.ncbi.nlm.nih.gov/articles/PMC13332746/ , https://mml.ikim.nrw/ait/projects/medshapenet/
- **NVIDIA skill catalog** (checked live, `NVIDIA/skills` `skills.sh.json`): `nv-segment-ct` (VISTA3D wrapper, emits label-map evidence),
  `nv-generate-ct-rflow` (MAISI-family synthetic CT **+ masks**; ≈16 GB VRAM at 256³; explicitly "not safe as training data … without an
  independent quality review"), `dicom-series-to-volume`, `nv-reason-ct`. None is installed. A verifier is exactly the "independent
  quality review" that `nv-generate-ct-rflow` says it needs.

### 1.3 Slice-to-volume consistency

- RT-planning software reconstructs surfaces from planar contours with undisclosed, differing algorithms. This causes DVH differences
  that are "pronounced in small structures" (SlicerRT, Sunderland 2015).
  https://labs.cs.queensu.ca/perklab/wp-content/uploads/sites/3/2024/02/Sunderland2015-manuscript.pdf
  The repo's use-case map already cites 10–20 % volume error on small organs from RTSTRUCT conversion.
- SAM2-style slice propagation creates false positives in empty slices (NVIDIA review above). MedVol-R1 adds a cross-slice Dice reward to fight
  incoherence. **Nobody found uses a 3D Betti certificate to detect z-artifacts** (a tunnel or cavity created only by stacking).
  That is the repo's H6 hypothesis.

### 1.4 3D topology in learning

- **Betti matching** (ICML 2023) and **Betti-Matching-3D** (2024): optimized C++ persistent homology, faster than Cubical Ripser. Next steps
  listed by the authors: "accelerated through GPU-based implementations". https://arxiv.org/html/2407.04683 , https://www.mdsi.tum.de/en/mdsi/research/kickstarter/segmatch3d/
- **Multiclass Betti matching** (2024): improves Betti errors on ACDC, TopCoW (13 classes), Platelet. https://arxiv.org/html/2403.11001
- Persistent-homology multi-class CMR loss "practical within high resolution 3D data for the first time" (Byrne 2022) https://pmc.ncbi.nlm.nih.gov/articles/PMC7614102/
- Topograph loss calls cost **96–603 ms** (repo map, citing https://arxiv.org/html/2411.03228v2).
- **Gap:** these are *differentiable losses on probability maps* for dense networks. A **non-differentiable, exact b0/b1/b2 reward
  for a generative policy (VLM/LLM) in 3D** has not been published. Nearest: the CS224r CAD project's B-rep Betti-1 match
  (86.7 % → 88.3 % after a short GRPO run, not significant), and CORE-Seg's reward "penalizes disconnected regions" in 2D.
  https://cs224r.stanford.edu/projects/pdfs/Gaurav%20Tyagi%20submission_416087944/CS224r_Final_Report.pdf , https://arxiv.org/html/2603.05911v1

### 1.5 GPU voxel/sparse/geometry libraries (what to build on, not rebuild)

| Library | Status (date) | Use for the verifier |
|---|---|---|
| **fVDB** (NVIDIA, SIGGRAPH/TOG 2024; now `openvdb/fvdb-core`) | Active. Release notes report sparse-grid construction 246 → 1.85 ms (~133×) at 11.5M voxels on RTX PRO 6000 Blackwell. Blog: "4× spatial scale, 3.5× faster, 10× more operators" vs prior frameworks. https://fvdb-core.readthedocs.io/stable/release_notes.html , https://developer.nvidia.com/blog/building-spatial-intelligence-from-real-world-3d-data-using-deep-learning-framework-fvdb/ | Sparse voxel container for large volumes (whole-CT, scenes). Useful once 64³ dense crops stop being enough. |
| **NanoVDB/OpenVDB** | Industry standard; NanoVDB = GPU read-only VDB | Storage and transfer format for certified volumes. |
| **NVIDIA Warp** 1.18 (2026-10-05) | New `warp.geometry` (sparse marching cubes, swept-volume meshing), `wp.Mesh(support_winding_number=True)`, BVH/mesh queries, `mesh_query_point_sign_parity`. https://github.com/NVIDIA/warp/blob/main/CHANGELOG.md | **Inside/outside and distance queries for mesh → voxel conversion**, BVH triangle tests for self-intersection. Swept volume is directly a robot-clearance primitive. |
| **Kaolin** v0.18 (2025-08-08) | Still maintained, but focus has moved to physics (Simplicits), splats and USD I/O. https://github.com/NVIDIAGameWorks/kaolin/releases | Representation conversions; not a core dependency. |
| **cuCIM / MONAI** | cuCIM has GPU morphology and labeling; MONAI hosts VISTA3D/MAISI | Data loaders and baselines for 3D medical. |
| **MinkowskiEngine / torchsparse** | Mature sparse convolution | Only needed for a learned 3D model; irrelevant to the verifier. |
| **Manifold** (elalish) | Guaranteed-manifold mesh booleans; Thrust-based (CPU or CUDA). In OpenSCAD: 5–30× faster than CGAL fast-csg, and 100–1,000× faster than the old CGAL Nef path. https://github.com/elalish/manifold/discussions/383 | **Exact-ish solid booleans for IoU on meshes.** trimesh's boolean backend in cadrille-style evals is Manifold-based (manifold3d). |
| **Generalized winding numbers on GPU** | "Antipodal method" (2026): **10⁹ queries/s** on an RTX 4090 for moderately complex meshes, full precision, handles non-manifold input. https://arxiv.org/html/2605.01536 | **Core primitive for certified mesh → voxel occupancy** (Section 5.2). |
| **CuACD** (2026-09-23, Hao Su lab) | First fully GPU-resident convex decomposition, >10× faster than CoACD-class baselines at matched quality. https://arxiv.org/abs/2609.28731 | Collision-proxy generation for the Physical-AI arm. Its warp-per-task design matches the repo's "irregular geometry on GPU" style. |
| **EMBER** (SIGGRAPH 2022) | Exact mesh booleans via homogeneous integer plane coordinates + generalized winding numbers; CPU, multi-threaded. https://doi.org/10.1145/3528223.3530181 | Reference oracle for exact booleans in the verifier's test suite. |

---

## 2. CAD: LLM/VLM CAD generation, datasets, RL with geometric verifiers

### 2.1 Methods (facts, dates, numbers)

| Method | Date / venue | Representation | Key numbers | Reward / feedback |
|---|---|---|---|---|
| **DeepCAD** | 2021 | Command tokens | 46.7 % IoU (DeepCAD test, point cloud) | – |
| **CAD-Recode** | arXiv 2024-12, ICCV 2025 | CadQuery Python, Qwen2-1.5B + 1-layer point projector | Trained on **1M procedurally generated** programs: DeepCAD 92.0 % IoU, median CD 0.16, IR 0.4 %; Fusion360 87.8 % IoU. SFT: **one H100, ~12 h, 100k iterations, batch 18**. Data generator validates with OCC `BRepCheck_Analyzer`. https://cad-recode.github.io/ , https://openaccess.thecvf.com/content/ICCV2025/supplemental/Rukhovich_CAD-Recode_Reverse_Engineering_ICCV_2025_supplemental.pdf | SFT only |
| **cadrille** | arXiv 2025-05-28, **ICLR 2026 Oral** | CadQuery; Qwen2-VL-2B; point cloud / multi-view / text | Multi-view RL: DeepCAD 92.2 % IoU, IR 0.0; Fusion360 84.6 %; real-scan CC3D 65.0 %. Point cloud RL: 90.2 / 85.0 / 67.9 %. https://m3l.am/publications/cadrille | **R = 10·IoU − 10·[invalid]**. Hard-example mining (mean SFT reward < 7.5). Dr. CPPO, G=16, N=4, 8×H100. IoU = trimesh boolean intersection of each GT component × each predicted component (`evaluate.py`). https://arxiv.org/html/2505.22914 , https://github.com/col14m/cadrille/blob/master/evaluate.py |
| **CAD-Coder** (text → CadQuery) | 2025-05 / v2 2025-10 | CadQuery; Qwen2.5-7B | 110K text–CadQuery–mesh triplets. GRPO took **146 h on 8×A800** (SFT 7 h). https://arxiv.org/html/2505.19713 | Chamfer distance + format |
| **CAD-RL / ExeCAD** | 2025-08 | CadQuery | 99.63 % executability, IoU 78.7 %. 16,540 examples. https://arxiv.org/html/2508.10118v2 | exec gate × (IoU + GPT-4o judge) |
| **ConCAD** (image → CAD) | 2026-09-28 | CadQuery | DeepCAD IoU 0.717 → 0.734 (SFT → RL), IR 0 %, G-CSR 0.913. Zero2CAD single view: IoU 0.443 → 0.507 (IR 23 %). Qwen3-VL-2B zero-shot: **IR 98 %, IoU 0.160**. https://arxiv.org/html/2609.34494 | exec gate + code-level constraint reward + **Monte Carlo IoU, 8,192 points**, max over 5 inertia-frame sign flips |
| **Ortho2CAD** | 2026 | CadQuery from orthographic drawings | 1M+ rendered drawings. RL with solid IoU reward → 100 % valid code. https://arxiv.org/html/2607.08891v2 | solid IoU via STEP |
| **CME-CAD** | CVPR 2026 | CadQuery | multi-expert RL | format + exec + **IoU (Jaccard)** https://openaccess.thecvf.com/content/CVPR2026/papers/Niu_CME-CAD_Heterogeneous_Collaborative_Multi-Expert_Reinforcement_Learning_for_CAD_Code_Generation_CVPR_2026_paper.pdf |
| **GACO-CAD** | 2025-10 | single image + depth/normal priors | RL with group length (conciseness) reward https://huggingface.co/papers/2510.17157 | geometry + length |
| **CADFusion** | ICML 2025 (MSR) | LLaMA-3-8B, command sequence | LVM score 8.96/10, IR 6.2 %; GPT-4o IR 93 %. https://arxiv.org/html/2501.19054 | **DPO from a VLM (LLaVA-OneVision-7B) scoring renders** |
| **CADCodeVerify** | ICLR 2025 | VLM-generated yes/no questions on renders | GPT-4: −7.30 % point-cloud distance, +5.5 % compile rate. CADPrompt: 200 prompts. https://proceedings.iclr.cc/paper_files/paper/2025/file/81a934cd364e18ea6fdeaf57a93c17d4-Paper-Conference.pdf | VLM feedback loop, no geometric solver |
| **FlexCAD** | ICLR 2025 (MSR) | LLaMA-3-8B + LoRA, hierarchy masking | PV 93.3 %, human realism 42.1 % vs GPT-4o 12.8 % (via https://cadarena.dev/eval-review.html) | – |
| **CAD-GPT** | 2024-12 | LLaVA-1.5-7B + spatial tokens | Image IR 1.61 %, median CD 9.77; text IR 7.43 %, CD 28.33 (same source) | – |
| **Text-to-CadQuery** | 2025-05 | CadQuery | **170K** text–CadQuery pairs on Text2CAD; exact match 58.8 → 69.3 %, CD −48.6 %. https://arxiv.org/html/2505.06507 | SFT |
| **RLCAD gym** | 2025 | B-rep → command sequence incl. revolve | 42 ms/step at 16 envs vs **1,674 ms** in the Fusion 360 Gym (39×). https://arxiv.org/html/2503.18549v1 | geometric difference in gym |
| **BenchCAD** | 2026-05 | CadQuery | 17,900 execution-verified programs, 106 industrial part families (gears, springs, drills). https://huggingface.co/papers/2605.10865 | benchmark |
| **UniCAD-MLLM** | 2026-06 | multi-task | 89.7 % IoU DeepCAD (leaderboard) https://opencodepapers-b7572d.gitlab.io/benchmarks/cad-reconstruction-on-deepcad.html | – |

**Not found or not verified** in this pass: CAD-Llama and Img2CAD numbers. They are known 2024–25 works, but I did not reach a
primary source with numbers, so they are not quoted.

### 2.2 Datasets

| Dataset | Size | Note |
|---|---|---|
| DeepCAD | 178,238 sketch-extrude models (Onshape) | everyone evaluates here |
| Text2CAD / text-annotated DeepCAD | ~170K text–sequence pairs (Text-to-CadQuery adds CadQuery code) | text conditioning |
| CAD-Recode | ~1M procedural CadQuery programs (+~1k val) https://huggingface.co/datasets/l2labs/cad-recode-v1.5 | procedurally valid by construction (OCC BRepCheck) |
| ABC | 1M+ B-rep models; "around 2 CPU years to extract triangle meshes" https://arxiv.org/html/1812.06216v2 | needs healing; large-scale verifier stress test |
| Fusion 360 Gallery | 8,625 reconstruction models | secondary test |
| CC3D | real scans | real-world gap (cadrille 65–68 % IoU) |
| BenchCAD | 17,900 programs / 106 families | industrial |
| Zero2CAD | 100K (test subset used by ConCAD) | harder; single view |

### 2.3 Where an exact GPU verifier fits CAD-RL (and where it does not)

Today's reward stack, per rollout:
`program text → Python/CadQuery → OCC builds B-rep → tessellate → (trimesh boolean | 8,192-point MC | Chamfer on samples) → scalar`.

- **Approximation noise in the reward.** **[est]** For Monte Carlo IoU with N = 8,192 points, if the union fills ~40 % of the joint box
  (~3,300 points), the standard error at IoU ≈ 0.73 is ≈ sqrt(0.73·0.27/3,300) ≈ **0.008**. ConCAD's whole RL gain on DeepCAD is
  0.717 → 0.734 (**+0.017**). GRPO ranks rollouts *within a group*, where IoU differences are often of this size. So per-rollout
  reward noise is comparable to the signal being learned. An exact or interval-certified IoU removes this noise. This is a testable,
  publishable claim: "does reward exactness change what GRPO learns?"
- **Topology is under-used.** Volume IoU cannot tell a through-hole from a blind pocket of equal volume (ConCAD makes the same
  argument for design intent). b1 (through-holes and handles) and b2 (internal voids) are exact, cheap invariants of a closed solid.
  Only one non-peer-reviewed project used B-rep Betti-1 as a CAD reward.
- **Validity is binary today** (code runs → 1). A verifier can grade *how* a solid is invalid: open boundary edges, non-manifold edges,
  self-intersecting faces, inverted shells, zero-volume slivers. This gives the repo's 13 defect classes a 3D analogue and a
  denser reward for near-miss programs.
- **Where it does not help:** OCC execution itself. Evidence that the geometric kernel is the CPU bottleneck: (a) block-decomposition RL (ECCV 2026),
  where OCC booleans take **>30 of ~36 min per model** and RL updates <2.5 min
  (https://en.papernotes.org/ECCV2026/vlm_efficiency/lightweight_online_reinforcement_learning_for_block_decomposition_of_cad_models/);
  (b) GRPO straggler analysis: one 30-s timeout in 128 completions moves GPU idle from 7.5 % to 55 % (https://taylorkolasinski.com/notes/price-the-reward-function/).
  A GPU verifier does not remove OCC. It does remove the *scoring* stage and makes scoring deterministic.
  Counter-claim from a vendor: SpacetimeML says its CPU STEP-pair verifier runs in **3.8 ms per solid** with a topology sub-score
  (vendor report, not peer-reviewed) https://spacetimeml.com/papers/SpacetimeML-CAD-CUA-Technical-Report.pdf . So "GPU is faster" alone is not a story.
- **A way around OCC for the verifier's own milestone:** restrict the program language to a subset whose solids can be built
  directly on GPU: sketch → extrude → union/cut of prisms and cylinders, i.e. the CAD-Recode generator's own op set.
  Prismatic extrusions of 2D rings are exactly what the repo's 2D ring kernels already certify. Extrude the certified ring along z
  to get a voxel prism. Unions/cuts become voxel boolean ops. This gives an OCC-free, fully GPU, exact (at grid h) reward for the
  sketch-extrude family, which covers DeepCAD's whole op set. It is a genuine bridge from the 2D checker to CAD.

---

## 3. Physical AI

### 3.1 NVIDIA stack facts (dates)

- **Cosmos Reason 1** (7B; GTC 2025; HF blog 2025-05-19). Physical AI SFT gives +10 %, RL another +5 %, average 65.7 on its benchmarks.
  https://developer.nvidia.com/blog/curating-synthetic-datasets-to-train-physical-ai-models-with-nvidia-cosmos-reason/
- **Cosmos Reason 2**: 2B/8B released 2025-12-19, 32B on 2026-04-29; 256K context; 2D/3D point localization. 2B needs 24 GB, 8B needs 32 GB.
  https://github.com/nvidia-cosmos/cosmos-reason2 , https://huggingface.co/nvidia/Cosmos-Reason2-2B
- **Cosmos Predict 2.5 / Transfer 2.5** (released 2025-10-06; 2B and 14B; up to 30-s videos; trained on 200M clips). **Cosmos 3** (2026-06-01) unifies
  language, image, video, audio and action. **GR00T N1.6** (CES 2026-01-05) integrates Cosmos Reason. **Newton** physics engine beta (2025-09-29, Warp + OpenUSD,
  with Google DeepMind and Disney; Linux Foundation). **Isaac Lab 2.3** dexterous-grasp workflow; **Isaac Lab-Arena** for policy eval.
  https://nvidianews.nvidia.com/news/nvidia-releases-new-physical-ai-models-as-global-partners-unveil-next-generation-robots ,
  https://investor.nvidia.com/news/press-release-details/2025/NVIDIA-Accelerates-Robotics-Research-and-Development-With-New-Open-Models-and-Simulation-Libraries/default.aspx
- **Cosmos Reason as reward/critic**: NVIDIA's cookbook uses it as a reward model (0–1 score) and as a best-of-N video critic; in one example, 4
  evaluations per video, score = fraction artifact-free. **The cookbook itself warns: "a proxy for physical plausibility, not a ground-truth oracle …
  can miss fine-grained physics."** https://nvidia-cosmos.github.io/cosmos-cookbook/core_concepts/evaluation/reason_as_reward.html ,
  https://nvidia-cosmos.github.io/cosmos-cookbook/recipes/end2end/gr00t-dreams/post-training.html

### 3.2 Physical-plausibility evals for world models

| Benchmark | Headline number |
|---|---|
| VideoPhy (2024) | best model (Pika) adheres to caption + physics in **19.7 %** https://arxiv.org/html/2406.03520v1 |
| VideoPhy-2 (2025) | Wan2.1-14B joint score **32.6 %**, **22 %** on the hard subset; conservation laws are hardest https://ar5iv.labs.arxiv.org/html/2503.06800 |
| T2VPhysBench (2025) | best Newtonian score 0.56, conservation 0.29 https://arxiv.org/html/2505.00337 |
| Physics-IQ (2025) | "physical understanding is severely limited" across Sora, Runway, Pika, Lumiere, SVD, VideoPoet https://physics-iq.github.io/ ; a 2026 leaderboard reports Cosmos3-Nano 40.23 and Veo 3.1 34.99 on "Physics-IQ Verified" https://elemmire1.github.io/physis-lang/ |
| WorldBench (2026) | models produce realistic *trajectories* but wrong *parameters* (e.g. not 9.8 m/s²) https://arxiv.org/html/2601.21282 |
| PhysVidBench (2026) | best model 36.2 %; caption-grounded QA correlates with humans at r = 0.45–0.69 vs 0.21–0.47 for direct VLM scoring https://arxiv.org/html/2507.15824 |
| Geometry-based: WorldScore (ICCV 2025, DROID-SLAM reprojection error), PDI-Bench (ECCV 2026 workshop: scale–depth, 3D motion, rigidity residuals via SAM 2 + MegaSaM + CoTracker3), Orbia (27 models, 720 instances, DA3 depth) | https://arxiv.org/html/2504.00983 , https://eccv2026-3d-world-models.github.io/papers/files/18-pdi-bench-quantitative-video-world-model-evaluation-for-geometric-consistency.pdf , https://orbia-bench.github.io/paper/orbia.pdf |

**Reading:** the field is moving from VLM judges toward *geometric residuals*. But every geometric eval begins with a **learned 3D lifting**
(SLAM, monocular depth, tracking), so the result is not certified. An exact verifier adds value only on **explicit** geometry, i.e. when
the world model or generator outputs meshes, voxels, USD or CAD, or in simulation where ground-truth geometry exists. Claiming
"exact physics verification of video world models" would be hype.

### 3.3 Sim-ready assets: where the exact-geometry gap is concrete

- **SimReady Foundation requirement index**: "Mesh geometry must be manifold" and "Mesh topology must be valid" have validators, but
  "**Transparent physical bodies should be watertight**" (no validator listed), "**Meshes should not overlap unnecessarily**" (none),
  "Meshes should not share the exact same space" (none), "Repeated occurrences … identical connectivity" (Manual).
  https://nvidia.github.io/simready-foundation/latest/indexes/requirements_index.html
- **SimReady FAQ: "SimReady does not yet define scene-level specifications or validators … many scene-level issues surface only at runtime."**
  https://docs.omniverse.nvidia.com/simready/latest/simready-faq.html . Validation is `simready-validate` (Python `BaseRuleChecker` rules on USD prims), and the
  result can be stamped into `customLayerData`. That gives a ready place to put a *certificate*.
- **Collision proxies matter**: CoACD vs V-HACD collision shapes on 49 SAPIEN drawers, SAC success **49 % → 80 %**, because V-HACD fills handle holes.
  **This is a b1 (handle hole) failure in the collision proxy, and a b1 check would flag it.** https://arxiv.org/html/2205.02961v1
- **Generated assets**: EmbodiedGen V2 (2026-08) runs a generate–verify–retry loop with topological repair + CoACD. EmbodiedGen notes raw
  TRELLIS output lacks "physically plausible geometry". Seed3D 1.0 (2025-10) claims watertight, manifold output imported into Isaac Sim.
  Gen2Physics (2026) is "the only approach to output watertight per-material sub-meshes". PhysX-Anything (2025-11), SIMART (2026-03),
  RFW-A (2026-10: **400K articulated assets** over 500 categories via agent-written procedural generators).
  https://arxiv.org/html/2607.07459 , https://arxiv.org/html/2510.19944v1 , https://arxiv.org/html/2608.23869 , https://huggingface.co/papers/2610.40333
  **All these pipelines verify with heuristics or repair tools. None emits a certificate (watertight + manifold + no self-intersection +
  b0/b1/b2 + proxy-covers-visual within ε).**
- **NVIDIA catalog skills relevant here** (not installed; recommendation only): `omniverse-cad-to-simready` (conversion + SimReady conformance
  orchestrator), `physical-ai-physics-selection` (Newton vs PhysX vs MuJoCo), `using-ovnewton`, `physical-ai-neural-reconstruction`,
  `warp-eval` (evidence-gated "is this hot path a Warp candidate"), `tao-finetune-cosmos-reason`. Install form, if the owner wants:
  `npx skills add nvidia/skills --skill omniverse-cad-to-simready --agent claude-code --global --yes`.

### 3.4 Synthetic data (Replicator, procedural)

- Omniverse Replicator writers and COCO export are hand-built (repo map). CAD-Recode (1M), RFW-A (400K) and VVR (10,000 tasks) show the pattern:
  **a procedural generator plus a deterministic verifier gives unlimited labelled data at any difficulty.** For geometry, the verifier is
  what makes the data trustworthy; CAD-Recode filters with OCC BRepCheck.

---

## 4. What the repo already has (verified in code)

- `src/objdel/checker/volume.py` (64 lines, CPU, NumPy/SciPy): `stack` (slice masks → volume), `euler` (V − E + F − C of the closed-cube
  complex), `betti3d` (b0 via 26-connected labels, b2 via 6-connected background components not touching the border,
  b1 = b0 + b2 − χ), `topology_any` (pred's (b0, b1, b2) ∈ set of readers' triples).
- ROADMAP Phase 7 plans the GPU version: χ from the 256 bit-cube patterns, a double certificate (bit-cube χ equals half the surface V − E + F;
  b1 from labels and χ), "mesh volume = voxel count", a whole-nodule reward, and the 2.5D input arm (H6, H6b), ~15 GPU-h.
- LIDC census: topology disagreement in **114 of 1,877 (6.1 %)** nodules after filtering; different 3D component counts in 16 (0.9 %). Clusters
  with a reader cavity: **92 / 2,654 raw, 29 / 2,654 filtered** (my count from `clusters.csv`, columns `cavity_any_0/5`).
- 2D CUDA checker is bit-identical to CPU on 195,809 cases. Full RL (Qwen3-VL-2B, 770 steps, one RTX 5090, 5 h 52 min, ~$4.56): test IoU 0.645,
  boundary F 0.816. Shared-prefix forward: 23.6 s/step.

### 4.1 Code review of `volume.py` as the oracle for a GPU 3D verifier

| # | Line | Finding | Severity |
|---|---|---|---|
| 1 | 27–30 (`stack`) | Missing slice indices between z0 and z1 become **empty planes**. A reader who skipped one interior slice gets **two 3D parts (b0 = 2)** and a false `topology_any` mismatch. Before using b0 as a reward, check how many LIDC reader volumes have interior gaps. The census's 32 raw `b0_3d_disagree` clusters may include some. | Medium (reward correctness) |
| 2 | whole module | **No voxel spacing.** LIDC slice thickness is 0.6–5 mm (ROADMAP). With 26-connectivity on anisotropic voxels, two blobs on adjacent 5-mm slices touching only at a corner count as one part. b1/b2 then mean different things at different thickness. Record spacing with each certificate and report results split by thickness (ROADMAP H6b already plans to record it). | Medium (interpretation) |
| 3 | 34–48 (`euler`) | Correct: closed-cube complex, a cell counted when any incident cube is present, padding avoids edge effects. It matches the 26/6 dual pair used for b0/b2, so b1 = b0 + b2 − χ is consistent. Good oracle for the planned 256-pattern bit-cube GPU χ. | OK |
| 4 | 57–58 | b2 = 6-connected background components of the padded volume − 1. Correct. | OK |
| 5 | 62–64 (`topology_any`) | Recomputes every reader's Betti triple per call. Cache per nodule when used as a reward (CPU cost only; irrelevant once on GPU). | Low (perf) |
| 6 | tests | `tests/test_volume.py` covers hollow cube, ring, corner contact and stacked hollow nodule (README). Add: interior z-gap, anisotropic thin-slice case, a 2×2×2 checkerboard (26-connectivity stress), and random-volume fuzz against an independent χ (e.g. bit-cube histogram) once the GPU port exists. | Suggestion |

Verdict: sound as an oracle. Fix or flag #1 before a b0 reward ships. Carry spacing (#2) in the certificate.

---

## 5. System sketch: one certified verifier, four input kinds

### 5.1 Requirements

- **Functional:** given a batch of candidates (2D rings | stacked slice masks | voxel grids | closed triangle meshes from a CAD program or a
  generator) and a *set* of references, return per candidate: a validity verdict + defect counts, (b0, b1, b2), volume, IoU / Dice / surface
  F-score against each reference (with intervals for meshes), set-valued scores (max, Hungarian-matched, topology-any), and a hash.
- **Non-functional:** device-resident (pixels and voxels never return to the host, matching the repo rule); deterministic and bit-identical to a CPU
  oracle; batched across rollouts (GRPO group × batch); latency ≪ policy step (the repo's step is 23.6 s).
- **Constraints:** one consumer GPU (RTX 4090/5090, 24–32 GB); a small team; the existing CUDA checker code style.

### 5.2 Components

```
 candidates ──► [adapter]──► canonical form ──► [certify] ──► [score vs reference set] ──► reward / eval row / data-filter verdict
                  │                               │                    │
   2D rings ──────┤ (existing ring kernels)       │ validity + defects │ IoU/Dice/F, intervals
   slice stack ───┤ stack → dense 64³–128³        │ bit-cube χ, b0,b2  │ max / Hungarian / topology-any
   voxel grid ────┤ as is (dense or fVDB sparse)  │ b1 = b0 + b2 − χ   │ GED decomposition
   mesh/STEP ─────┘ edge-sort watertight check,   │ volume = count·h³  │
                    BVH self-intersection,        │                    │
                    winding-number voxelization   └──► certificate (hash, counts, grid h, interval widths)
```

- **Mesh certification without voxels (exact combinatorics):** sort directed edge keys. Each undirected edge must appear exactly twice with opposite
  orientation (watertight + consistently oriented); edges seen ≠ 2 are open or non-manifold. This uses the same sort → pair → pointer-jump
  pattern as the repo's GPU ring linking. Self-intersection uses BVH triangle–triangle tests (Warp BVH). Then χ(solid) = χ(boundary)/2 from
  V − E + F; b0 from components plus shell nesting; b2 = number of inner shells; b1 = b0 + b2 − χ. All exact for the tessellation.
- **Certified IoU interval for curved solids (new idea):** at voxel size h, label a voxel centre *inside-certain* if its winding number is ≥ ½ and
  its distance to the surface is > h·√3/2. Label it *outside-certain* if the winding number is < ½ and the same distance holds. Otherwise it
  is *boundary*. Then |A∩B| and |A∪B| have hard lower/upper bounds, which give **IoU ∈ [lo, hi]** with width shrinking as h drops. This replaces
  "8,192-point Monte Carlo ± ~0.008" with a *guaranteed* interval. Winding-number throughput of 10⁹ queries/s (Antipodal, RTX 4090) means
  128³ = 2.1M queries per solid ≈ **2 ms [est]**. A GRPO step of 128 prompts × 16 rollouts = 2,048 solids ≈ **4 s at 128³ or ≈0.5 s at 64³ [est]**,
  small next to OCC execution and generation. Must be measured.
- **Topology stability flag (heuristic, not a proof):** compute (b0, b1, b2) on the inside-certain set and on inside-certain ∪ boundary. If they agree,
  mark the voxel topology "stable at h"; if not, report the mesh-combinatorial Betti numbers instead.

### 5.3 Trade-offs

| Choice | Gain | Cost | Revisit when |
|---|---|---|---|
| Dense 64³–128³ voxels | simple kernels, reuses the bit-cube plan | memory grows as n³; thin CAD features alias | whole-CT or scene scale → fVDB sparse |
| Mesh-combinatorial topology | exact, resolution-free | requires valid tessellation; OCC on CPU | GPU B-rep evaluation (out of scope) |
| Interval IoU vs point estimate | certified, no Monte Carlo noise | needs distance queries; boundary band costs memory | if intervals are too wide at 64³, go adaptive (octree refine only boundary cells) |
| Restricting CAD to sketch-extrude | fully GPU, exact, no OCC | loses fillets, revolves, lofts | after the first paper; add revolve via 2D ring sweep |

### 5.4 Back-of-the-envelope (all **[est]**, to be measured)

| Quantity | Chain | Result |
|---|---|---|
| Solids per GRPO step | 128 prompts × 16 rollouts (cadrille's G=16) | 2,048 |
| Winding-number queries per solid | 128³ ≈ 2×10⁶ (64³ ≈ 2.6×10⁵) | — |
| GPU time for occupancy | 2,048 × 2×10⁶ ÷ 10⁹ q/s (Antipodal, RTX 4090) | ≈ 4 s at 128³, ≈ 0.5 s at 64³ |
| Dense 64³ volume bytes | 2.6×10⁵ B (1 B/voxel) × 2,048 | ≈ 0.5 GB per step: fits a 32 GB card; bit-packing cuts it 8× |
| Bit-cube χ pass | memory-bound read of 0.5 GB at ~1 TB/s order | ~ms; labeling (b0, b2) dominates and costs more, but is still ≪ 1 s |
| OCC execution (CPU) | 2,048 programs × (unmeasured, assume 50–200 ms) ÷ 16 cores | ≈ 6–25 s: **the bottleneck**, comparable to the repo's whole 23.6 s step |

Decision these numbers force: the GPU verifier is cheap enough to run on every rollout at 64³ (and probably 128³). CAD *program execution*
must be async, with timeouts excluded from the group rather than scored 0 (see 5.6). Otherwise rewards come late or are wrong.

### 5.5 ADR-001: substrate for the 3D/CAD verifier

**Status:** Proposed. **Context:** the 2D checker is custom CUDA, bit-identical to a CPU oracle. 3D adds voxels (dense now, sparse later) and
meshes (winding numbers, BVH tests, edge sorting).

| Option | Complexity | Exactness control | Speed | Team fit |
|---|---|---|---|---|
| A. Own CUDA kernels for voxels (bit-cube χ, labeling) + **Warp** for mesh queries (winding number, BVH, distance) + **fVDB** only when sparse is needed | Medium | Full for voxels; Warp's winding number is float (needs the inside/outside band to absorb error) | High | Good (extends existing `cuda/`) |
| B. Everything in own CUDA | High | Full | Highest | Long build; reimplements BVH and winding numbers |
| C. CPU OCC + trimesh/Manifold (status quo of CAD papers) | Low | Approximate (MC, tessellation tolerance) | Low; process pool | Easy, but no contribution |

**Decision (proposed):** A. Keep the certificate-critical integer paths (χ, labels, edge pairing) in own CUDA with a CPU oracle. Use Warp for float
geometric queries and *bound their error* through the inside/outside band (Section 5.2). Keep C as the **reference oracle in tests**, not in the reward.
**Consequences:** easier: reuse of the bit-identical test harness. Harder: two GPU stacks (CUDA C++ and Warp) in one process. Use `warp-eval`'s
evidence-gated method before committing, and measure that Warp queries and torch/vLLM share memory pools without fighting
(Warp 1.18 added `WARP_ENABLE_MEMPOOLS_AT_INIT=0` for that case). Revisit if Warp winding-number error forces bands so wide that intervals are useless.

### 5.6 Failure modes: how a "certified" reward can still be wrong

| Symptom | Root cause | Detection | Guard |
|---|---|---|---|
| Reward drops to 0 on correct programs | OCC timeout/crash scored as invalid; GRPO normalization then **flips advantages** (a correct rollout scored 0 in an all-correct group gets −2.64; Kolasinski 2026) | log timeout rate per step; compare timed-out programs re-run offline | exclude timed-out rollouts from the group (mask, do not score 0); hard per-program time limit |
| IoU high but wrong part | frame/scale mismatch (cadrille normalizes to the unit box; ConCAD uses inertia frames with sign flips) | report IoU under both conventions on a fixed probe set | fix one canonical frame per task in the certificate |
| b1/b2 flicker between resolutions | thin walls alias at grid h | topology-stability flag (inside-certain vs inside ∪ boundary disagree) | fall back to mesh-combinatorial Betti numbers; never reward unstable voxel topology |
| "Valid" solid with open edges | tessellation tolerance (OCC `tessellate(0.001, 0.1)` in cadrille's eval) leaves cracks | edge-pair count ≠ 2 | tighter tolerance or B-rep-side check; report as a defect class, not a crash |
| Reward hacking via degenerate solids | zero-volume slivers or huge boxes can win under some IoU definitions | volume and bounding-box sanity bands; G-CSR-like constraint checks | graded validity term; cap reward when defects > 0 |
| GPU vs CPU disagreement | float winding number near the surface | bit-compare the integer outputs (χ, labels); interval containment test for float outputs | CI fuzz, like the 2D checker's 195,809-case gate |

### 5.7 Model choice for one 24–32 GB GPU (no official HF leaderboard covers this; choice from papers + metadata)

| Role | Candidate | Params | License | Evidence | Fit |
|---|---|---|---|---|---|
| Policy for image → CadQuery (research) | `maksimko123/cadrille-rl` | 2.21 B | **CC-BY-NC-4.0** | DeepCAD 92.2 % IoU multi-view | yes (bf16) |
| Policy (repo continuity, permissive) | `Qwen/Qwen3-VL-2B-Instruct` (or Qwen3.5-2B, Apache-2.0, per main-thread lookup) | 2.13 B | Apache-2.0 | zero-shot CAD is poor: IR 98 %, IoU 0.160 (ConCAD), so it needs SFT on CAD-Recode first | yes |
| Point cloud → CadQuery baseline | `filapro/cad-recode-v1.5` | 1.54 B | **CC-BY-NC-4.0** | DeepCAD 92.0 % IoU | yes |
| 3D CT baseline | VISTA3D via `nv-segment-ct`; MedSAM2 | — | VISTA3D NVIDIA terms; MedSAM2 CC-BY-SA-4.0 | Section 1.2 | yes |
| Physical critic (comparison only) | `nvidia/Cosmos-Reason2-2B` | 2.44 B | NVIDIA "other" | proxy, not oracle | needs ≥24 GB |

The non-commercial licenses on cadrille and CAD-Recode matter if the owner wants a commercial verifier service. Train a permissive policy (Qwen)
for any product path, and use cadrille only as a research baseline.

---

## 6. Three candidate strategic directions

### Direction A: "Certified geometric rewards": one verifier for 2D → 3D → CAD (recommended first)

**Claim to test:** *exact, topology-aware, batched geometry rewards change what RL learns, compared with approximate CPU rewards*, at equal
compute, across 2D outlines, 3D volumes and CAD solids.

**Why novel (who has / has not):**
- Has: CAD-RL with CPU geometric rewards (cadrille trimesh IoU; ConCAD MC IoU with 8,192 points; CAD-Coder Chamfer; Ortho2CAD solid IoU; CME-CAD IoU);
  mesh-quality RL with heuristic integrity rewards (Mesh-RFT, Mesh-Pro); 3D-gen RL with learned rewards (AR3D-R1, DreamCS); dense topology
  *losses* (Betti matching 3D, which wants a GPU version); MedVol-R1's cross-slice Dice reward.
- Has not (none found): a **certified** reward (validity + defects + b0/b1/b2 + interval IoU) that is **GPU-resident and batched over rollouts**, and an
  ablation of **reward exactness** (MC vs exact vs interval) and **topology terms** in GRPO for CAD or 3D medical outlines.

**First milestones on one GPU:**
1. *3D LIDC (already planned, Phase 7):* GPU bit-cube χ + labels, certificate checked against `volume.py` on all LIDC volumes; whole-nodule reward;
   measure false tunnels/cavities per stacked volume before vs after RL (H6) and with the 2.5D arm (H6b). ROADMAP estimate ~15 GPU-h.
   Metric: z-artifact rate (b1, b2 not present in any reader), per-slice IoU/NSD unchanged.
2. *Verifier vs CAD eval agreement:* run the GPU mesh certifier + interval IoU on cadrille's released DeepCAD/Fusion360 predictions (HF checkpoints exist).
   Report the share of trimesh-IoU values inside the certified interval, the interval width at 64³ and 128³, and the time per solid. Report how
   many "valid" (executable) predictions are non-watertight, self-intersecting or have spurious b0 > 1. **No training needed; ~hours on one GPU [est].**
   This milestone alone is a short paper or workshop result ("hidden invalidity in state-of-the-art CAD reconstructions").
3. *Small CAD RL:* Qwen3-VL-2B (the repo's model; ConCAD shows it zero-shot at IR 98 %), SFT on a CAD-Recode subset rendered to images, then GRPO with
   reward arms: (i) MC-8192 IoU, (ii) interval-midpoint IoU, (iii) (ii) + topology + graded validity. CAD-Recode's SFT took ~12 h on one H100
   for 100k iterations, so plan **10–20 GPU-h SFT + 10–20 GPU-h RL [est]**. Metrics: IoU, CD, IR, **b1/b2 match rate**, watertight rate, G-CSR-style
   constraint rate, and the change in within-group reward rank between arms.

**Datasets:** LIDC-IDRI (have), DeepCAD, Fusion360, CC3D (real scans), CAD-Recode 1M, BenchCAD 17.9K.
**Risks:** OCC execution stays the CPU bottleneck (mitigate: start with sketch-extrude solids built on GPU from certified rings). The exactness ablation might
show "no difference". That is still publishable as a negative result, but it weakens the pitch. CAD-RL is crowded (≥8 papers in 18 months), so the
contribution must be the verifier and the ablation, not another SOTA claim.

### Direction B: "Set-valued evaluation and rewards for ambiguous geometry" (scientific spine)

**Claim:** when the truth is a *set* (readers, occluded CAD backsides, multiple plausible futures), scoring against one reference punishes correct
answers. Set-valued scores (best-match, Hungarian-matched, topology-any, coverage of the set by K samples) are the right reward and eval, and they
can be certified.

**Why novel:** GED (Kohl 2018) and its known flaw (dominated by empty-vs-non-empty detection, https://arxiv.org/html/2407.16367v1); Dice_max / Dice_match
toolkits (MBH-Seg 2025, https://github.com/mbhseg/metrics-mbhseg2025); 3D GED 0.401 and Hungarian-matched 3D IoU 0.468 on LIDC (Viviers 2023,
https://arxiv.org/abs/2305.00950); CVPR 2026 multi-rater personalization uses GED as a regularizer. **Missing:** *topology* in set metrics; set-valued
*rewards* for generative policies; set-valued evaluation in CAD, where single-view image → CAD is ambiguous (ConCAD's Zero2CAD single-view IoU is
only 0.507, partly because the hidden geometry is unknowable).
Inter-reader variability is large: LIDC reader-vs-reader Dice 0.77 ± 0.09 on 928 nodules (RSNA 2014, https://archive.rsna.org/2014/14046034.html).
Only 34 % of LIDC nodules are outlined by all four readers (repo census).

**First milestone (one GPU or Mac):** a "LIDC-3D-Set" benchmark. Certified reader volumes, set metrics (best-match IoU, Hungarian IoU, topology-any, presence
F1, GED split into detection and shape terms). Baselines: MedSAM2 (box on the middle slice), VISTA3D point mode (`nv-segment-ct` wrapper), nnU-Net, and the
repo's VLM stacked. Cheap: mostly inference. A CAD analogue: render DeepCAD solids from *one* view, define the reference set as all solids
consistent with the visible silhouette (sample by perturbing hidden-face parameters), and score set coverage.
**Risks:** medical multi-rater evaluation is a crowded niche; reviewers may see it as incremental unless topology-any and certification change rankings.
b2 cases are rare (29 filtered clusters), so 3D topology claims rest mostly on b0 and on artifact suppression.

### Direction C: "Verifier-driven synthetic data engine" (later; product-shaped)

**Claim:** procedural generator → certified verifier → certificate stamped into the asset gives trustworthy training data for (a) 3D medical masks
(filter `nv-generate-ct-rflow`/MAISI masks for anatomical topology priors: one liver, b0 = 1, b2 = 0, vessel-tree b1 ranges), (b) CAD programs
(CAD-Recode-style generator, with b1/b2 and constraint labels), and (c) sim-ready assets (watertight, manifold, no self-overlap, collision proxy
covers the visual mesh within ε and preserves its b1, so the "handle hole" failure is gone; GPU proxies via CuACD).

**Why novel:** generators exist (CAD-Recode 1M, RFW-A 400K, EmbodiedGen V2, Seed3D, PhysX-Anything), and NVIDIA's validator framework exists
(`simready-validate`, `BaseRuleChecker`). Missing: certified checks for the rules that SimReady lists **without validators** (watertight, overlap,
coincident meshes), any scene-level validator, and a *certificate* rather than a repair attempt.
**First milestone:** run the GPU mesh certifier over a public generated-asset set (e.g. the 2,023-mesh `torusprime/cobot-atlas` from TRELLIS 2 /
Hunyuan3D, https://huggingface.co/datasets/torusprime/cobot-atlas , or EmbodiedGen outputs). Publish the defect census: share watertight, manifold,
self-intersecting, b0 > 1, proxy-fills-a-hole. Then show one downstream effect in Newton/Isaac Lab (grasp or drawer task) between certified and
uncertified assets, echoing CoACD's 49 % → 80 %.
**Risks:** crowded by well-funded teams (NVIDIA content agents, Horizon Robotics, ByteDance). A census is useful but not "history-making". The
downstream effect needs simulation engineering outside the repo's core; GPU hours are moderate, setup time is large.

### Honest ranking

| | Novelty | One-GPU feasibility | Builds on what exists | Risk of being scooped |
|---|---|---|---|---|
| A: certified rewards | high (no certified GPU reward + exactness ablation found) | high (milestone 2 needs no training) | very high (2D CUDA checker, `volume.py`, Phase 7) | medium (CAD-RL moves fast) |
| B: set-valued ambiguity | medium-high (topology-any + CAD ambiguity is new) | very high | high (LIDC census done) | low-medium |
| C: data engine | medium | medium | medium | high |

**Recommended sequence:** Phase 7 (3D LIDC) → A-milestone 2 (certify cadrille outputs, no training) → A-milestone 3 (small CAD RL with the
exactness ablation), with B's set metrics used as the evaluation throughout. Start C only once the mesh certifier exists, as its first external user.

### What is hype (explicitly)

- "Exact verifier for world-model videos": no. Geometry comes from learned lifting; the verifier can only certify what the lifter outputs.
- "GPU makes CAD-RL fast": mostly no. OCC execution dominates, and CPU scoring can be milliseconds. Speed is a side benefit.
- "b2 rewards teach cavities in LIDC": no. Real cavities are 29 of 2,654 filtered clusters. The b2 term mainly suppresses stacking artifacts.
- "Cosmos Reason as ground-truth physics reward": NVIDIA's own docs call it a proxy.
- "CAD → 3D is a research problem": no. It is conversion engineering (OCC tessellation, Omniverse CAD → SimReady). The research problem is the inverse.

---

## 7. Open questions to settle by measurement before any spend

1. OCC execution time per CadQuery program on the box CPU, p50 and p99 (the p99 sets the GRPO straggler cost).
2. Interval-IoU width at 64³ and 128³ on DeepCAD solids, and the share of trimesh IoU values inside it.
3. Share of "valid" cadrille/CAD-Recode outputs that fail watertight, manifold, self-intersection or b0 = 1.
4. z-artifact rate (b1/b2 present in no reader) of the current 2D model's stacked LIDC volumes. This is the H6 baseline, and it comes from existing checkpoints at no cost.
