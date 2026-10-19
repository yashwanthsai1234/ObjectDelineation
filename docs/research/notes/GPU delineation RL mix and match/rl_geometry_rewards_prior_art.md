# Prior art for geometric, topological and verifiable rewards in RL fine-tuning of vision and vision-language models (gap-closing sweep, as of 12 Oct 2026)

Scope and method notes:
- **The proposal under test:** "RLVR for spatial outputs." An exact GPU geometry engine turns a binary mask into closed rings with holes, using 2x2 bit-quads. It returns validity, β0 (number of components), β1 (number of holes), χ, area, perimeter, vertex count and a canonical ring hash. These would serve as reward and verifier terms for GRPO-style RL of models that output masks or polygons.
- **This file builds on round 2** (`../GPU delineation use case map/r2_geometry_in_ml_loop.md`). Papers already verified there are cited back to their sources and are not re-fetched: Seg-Zero, SAM-R1, VisionReasoner, LENS, ALTo, RLRF and LithoGRPO.
- **Evidence grades**, used throughout:
  - **[P]**: verified from the paper's full text (arXiv HTML or PDF text that I extracted and grepped, or a full-text fetch).
  - **[A]**: abstract or landing page only.
  - **[S]**: search-engine snippet only, not opened.
- **Tooling caveat.**
  - The WebFetch summarizer **invented** content for one PDF. For arXiv 2607.09492 it reported IoU and segmentation reward-hacking examples. The extracted PDF text shows the paper covers safety VQA, chart VQA and stress tests, with no such examples.
  - All load-bearing claims below (Topo-R1, FloorplanVLM, SimpleSeg, TraversRL, Mesh-RFT, the part-grounding paper, VectorLLM, SegAgent, Point-RFT, Ground-R1, DenseMLLM) were therefore re-checked against raw extracted text with grep.
  - Claims from WebFetch summaries of arXiv HTML pages that I did not grep are still marked [P], because they came from full text. They carry slightly lower confidence.
- **Exa MCP was unavailable** (server failed to connect). Searches used WebSearch only.

---

## Q1. What rewards do the 2025–2026 perception and segmentation RL papers use (output type, reward terms, how computed)?

### Takeaway
Across 30-odd papers, the reward vocabulary is almost entirely **overlap and distance terms on boxes, points or masks** (IoU, GIoU, L1, Dice, S-measure, mAP, F1 after Hungarian matching), plus **format gates**. Some papers add **count or length penalties**, and a few add **ad-hoc geometric constraints** (part-in-object box containment, canvas overflow, negative-point placement).

Mask-output pipelines almost always compute the mask with a frozen SAM/SAM2 and score it by IoU or Dice. Only two papers score a polygon the policy emits itself:
- SimpleSeg: rasterized IoU, centroid distance, format.
- FloorplanVLM: a binary "watertight closed polygon" validity gate, plus IoU.

No paper in the table rewards β0, β1/holes, χ or ring validity of the model's own output mask or polygon.

### Cited Findings

**Reward table.** The Type column gives the output type.

| Paper (date) | Type | Reward terms and computation | Grade |
|---|---|---|---|
| Visual-RFT (2503.01785) | boxes; class labels | **Detection:** R_d = R_IoU + R_conf + R_format. R_IoU is the mean IoU of matched boxes. Confidence reward: r = c if IoU≠0, else 1−c. **Classification:** R_cls = R_acc + R_format. | [P] — [arXiv html](https://arxiv.org/html/2503.01785) |
| VLM-R1 (2504.07615) | boxes | **REC:** R = IoU(b*, f(o)). **OVD:** R = s·mAP with s = min(1, L_gt/L_pred), the "odLength" term. Format reward is 0/1. | [P] — [arXiv html](https://arxiv.org/html/2504.07615) |
| UniVG-R1 (2505.14231) | boxes (multi-image) | r = IoU + format (0/1). Samples are re-weighted by difficulty with φ(mIoU); exp(1−mIoU) worked best, ahead of −log(mIoU) and (1−mIoU)². | [P] — [arXiv html](https://arxiv.org/html/2505.14231) |
| Perception-R1 (2504.07954) | boxes, points, text | Format reward is ±1. Answer reward Φ is IoU for grounding, Euclidean point distance for counting, and edit distance for OCR. **Detection** combines an object-number reward (F1), a location reward (IoU) and a binary classification reward with a missing penalty. Multiple objects are matched with Hungarian argmax ΣΦ. | [P] — [arXiv html](https://arxiv.org/html/2504.07954) |
| Vision-R1 (Zhan et al., 2503.18013) | boxes | A "criterion-driven reward function" with multi-dimensional feedback, plus a "progressive rule refinement strategy" that adjusts the criteria "to mitigate reward hacking". The terms themselves were not extracted. | [A] — [arXiv abs](https://arxiv.org/abs/2503.18013) |
| Seg-Zero (2503.06520) | box + points → SAM | Format reward plus IoU and L1 accuracy rewards. | via r2 — [alphaxiv](https://www.alphaxiv.org/abs/2503.06520.md) |
| SAM-R1 (2505.22596) | box + point → SAM mask | Tiered mask IoU: 4 / 3 / 2 / 0 at the 0.8 / 0.7 / 0.5 thresholds. Two format rewards. | via r2 — [arXiv](https://arxiv.org/html/2505.22596) |
| VisionReasoner (2505.12081) | boxes + points | Box IoU>0.5, box L1<10 px and point L1<30 px, with Hungarian matching. States explicitly: "rather than binary masks". | via r2 — [arXiv](https://arxiv.org/html/2505.12081v5) |
| LENS (2508.14153) | box → SAM2 mask | λ1·format + λ2·box IoU + λ3·mask IoU. | via r2 — [arXiv](https://arxiv.org/html/2508.14153) |
| Seg-R1 (2506.22624) | points, box, labels → SAM2 | Format (0/1) + segmentation reward = **0.7·IoU + 0.3·S-measure**. The authors report that "using S-measure alone leads to reward hacking". | [P] — [arXiv html](https://arxiv.org/html/2506.22624) |
| GenSeg-R1 (2602.09701) | box + positive/negative keypoints → SAM2 | **Base:** format/JSON, box IoU>0.5 plus L1, point distance <119 (in 0–1000 space) and inside the box, and a non-repetition penalty. **-G variant adds:** SAM2-in-loop mask IoU (weight 5.0), negative-point validity (weight 10.0) and a no-target reward (weight 10.0). | [P] — [arXiv html](https://arxiv.org/html/2602.09701v1) |
| CORE-Seg (2603.05911) | `<seg>` token → SAM decoder mask | Format; r_bbox = mIoU_matched × F1 (Hungarian); **r_mask = GIoU(box(M̂), box(M)) if Dice<0.05, else 1 + 2·Dice** ("dual-granularity", aimed at reward sparsity). | [P] — [arXiv html](https://arxiv.org/html/2603.05911) |
| SegAgent (2503.08625) | click actions → iterative mask | R(s_t) = IoU(M_t, M_target). The policy is trained by imitation plus StaR-style filtering, with a process reward model and tree search at inference. This is not GRPO. | [P] — [arXiv html](https://arxiv.org/html/2503.08625) |
| SimpleSeg (2601.19228) | **polygon as text** | GSPO with three rewards. (1) "Mask IoU reward: the direct IoU between the predicted and ground-truth mask", which is 0 if IoU<τ. (2) Negative centroid MSE. (3) Format reward, which returns zero if the polygon-coordinate format is wrong. Without clockwise ordering "the coordinate sequence fails to form a valid polygon, and no segmentation mask can be derived". | [P] — [arXiv pdf](https://arxiv.org/pdf/2601.19228) |
| ALTo (2505.16495) | adaptive mask tokens | R = 𝕀_format + IoU − α·L. | via r2 — [arXiv](https://arxiv.org/html/2505.16495v1) |
| DenseMLLM (2602.14134) | dense per-patch label tokens | DAPO with "class-label IoU (label set overlap) for semantic segmentation". | [P] — [arXiv html](https://arxiv.org/html/2602.14134v2) |
| VectorLLM (2507.04664 v1) | **building polygons as tokens** | DPO, not RLVR. Rejected answers are made by "randomly corrupt[ing] the labels … by deleting points or inserting random points". DPO adds +1.8 AP and fixes "ignoring small edge lines and duplicate nodes". | [P] — [arXiv html](https://arxiv.org/html/2507.04664v1) |
| VectorLLM++ (Research Square rs-10949232) | polygons, masks, boxes | "reinforcement learning with verifiable rewards that directly optimize geometric quality". **The reward terms are unknown.** | [S] — [Research Square](https://www.researchsquare.com/article/rs-10949232/v1) |
| DeepEyes (2505.14362) | zoom-in boxes (tool calls) | R = R_acc + R_format + 𝕀[R_acc>0]·R_tool. **No box-IoU reward**; IoU is used only for post-hoc analysis. | [P] — [arXiv html](https://arxiv.org/html/2505.14362) |
| Ground-R1 (2505.20272, latest version) | evidence boxes | Answer correctness plus format. Diagnoses a **"scale-driven bias … training rewards are dominated by large visual regions"** and introduces SRPO with scale-binned advantages. | [P] — [arXiv html](https://arxiv.org/html/2505.20272) |
| Point-RFT (2505.19702) | points inside CoT | Format (regex) + accuracy 𝕀(ŷ=y) on ChartQA. **No geometric reward** on the points. | [P] — [arXiv html](https://arxiv.org/html/2505.19702v1) |
| Topo-R1 (2603.13054) | boxes + error type (4 classes) on a given tubular mask | R = 0.10·fmt + 0.85·acc + **0.05·topo**. See Q2. | [P] — [arXiv html](https://arxiv.org/html/2603.13054v1) |
| FloorplanVLM (2602.06507) | **floorplan polygons (JSON)** | R = 0.1·R_val + 0.5·R_ext + α·0.4·R_int. See Q2. | [P] — [arXiv html](https://arxiv.org/html/2602.06507) |
| Part-level grounding OP-HRG (2607.15374) | parent-object and part boxes (→ masks) | Parent-object IoU via Hungarian matching. **Part-containment reward:** the part box must lie inside an object box, must not equal it, and must have strictly smaller area. Also reflection-consistency and improvement rewards. | [P] — [arXiv html](https://arxiv.org/html/2607.15374) |
| QwenVL-SAMTok | mask codebook tokens | "purely textual answer-matching reward for the GRPO of the mask generation task". | [S] — [emergentmind topic](https://www.emergentmind.com/topics/discrete-mask-tokenization) |
| Veason-R1, STVG-R1, VideoSeg-R1, VideoSEG-O3, R1-Track | video | See Q4. | — |
| TraversRL (2607.17479) | pedestrian graph (polylines) | Terminal **buffered IoU**, plus a stepwise r_t = ΔI_t/(ΔA_t+ε). See Q2. | [P] — [arXiv pdf](https://arxiv.org/pdf/2607.17479) |
| GeoSVG-RL (2605.25447) | SVG diagrams | GRPO (G=4). Terms: rendering-validity gate; canvas fit and overflow, R = −Area(B_all∖C)/Area(B_all); arrow-anchor accuracy (12 px threshold); text-in-box containment and 6 px padding; graph edge F1; code cleanliness. Computed with headless Chromium `getBBox()` and lxml. | [P] — [arXiv html](https://arxiv.org/html/2605.25447v1) |
| Mesh-RFT (2505.16761, NeurIPS 2025) | 3D meshes | Masked DPO with automatic **Boundary Edge Ratio** BER = \|E_∂M\|/\|E_M\| ("quantifies holes and surface discontinuities") and a Topology Score (quad-friendliness), used to build preference pairs. | [P] — [arXiv html](https://arxiv.org/html/2505.16761v1) |
| CAD-Coder (2505.19713) | CadQuery code → solid | GRPO with Chamfer distance + format. A related snippet describes executability as a "multiplicative gating factor"; which paper it comes from is uncertain. | [S] — [arXiv](https://arxiv.org/abs/2505.19713) |
| LithoGRPO (2606.00228) | ILT masks (pixels) | −Σ normalized L2, PVB, EPE and **shot count** (fast CPU ILP). | via r2 — [arXiv](https://arxiv.org/html/2606.00228) |
| RLRF (2505.20793) | SVG | Render with CairoSVG, then L2, edge L2, DreamSim, CLIP and a length penalty. | via r2 — [arXiv](https://arxiv.org/html/2505.20793v2) |

### Inferences
- **The field's reward space has two axes:** metric (IoU, Dice, L1, F1) and matching (Hungarian, thresholds, tiers).
  - Nobody computes **structural invariants of the output**: component count, hole count, χ, ring simplicity or orientation, or multi-ring legality.
  - The likely reason is architectural. Most pipelines hand geometry to SAM, so the policy never emits a mask or polygon whose structure it controls.
  - Exceptions where it could matter: SimpleSeg, FloorplanVLM, VectorLLM(++), ALTo, SAMTok and DenseMLLM, whose policies emit geometry tokens directly.
- **SimpleSeg is the most natural host for the engine.** Its policy emits polygons, invalid polygons yield no mask at all, and its reward is IoU after rasterization. The engine would add a canonical re-polygonization and β/hole/validity terms there.

### Gaps
- **Not located or not verified:**
  - **RefAM.** I recall it as training-free "attention magnets" referring segmentation, not RL, but this is unverified.
  - **"Text4Seg + RL"**: no paper found.
  - **"LISA-RL"**: no specifically named variant found.
  - **Vision-R1 (Huang et al., 2503.06749)**: math-reasoning RL, reward not checked.
- **VectorLLM++'s RLVR reward terms are the single most important unknown.** It is a polygon-emitting MLLM with RLVR "that directly optimize[s] geometric quality". If those terms include validity, self-intersection or holes, it overlaps the proposal directly. Only a Research Square landing page was reachable.
- **Not opened:** GRPO-SEG, GRASP, FineRS, AURORA, CoPRS, RSAgent and Saliency-R1. Snippets suggest IoU-family rewards; RSAgent adds a process reward based on step-wise IoU improvement [S] ([arXiv 2512.24023](https://arxiv.org/pdf/2512.24023)).

---

## Q2. Is there any topology, validity, hole, Betti, Euler or connectivity reward in RL for vision models? (falsification search)

### Takeaway
**Yes, partially, and the proposal must not claim "first topology-aware reward" or "first polygon-validity reward".** Five published works overlap.
1. **Topo-R1 (Mar 2026)** puts a topology term (clDice) into a VLM GRPO reward. It is not computed on the policy's own output geometry, and Betti numbers appear only in data curation.
2. **FloorplanVLM (Feb 2026)** uses a binary "watertight closed polygon" validity reward in GRPO for a polygon-emitting VLM. Validity rises from 90.2% to 96.1%.
3. **Mesh-RFT (NeurIPS 2025)** uses automatic hole-sensitive topology scores (BER, TS) to drive DPO for 3D mesh generators.
4. **TraversRL (Jul 2026)** improves graph connectivity through RL, but its reward is buffered IoU.
5. **Learning to Trace (NeurIPS 2026)** trains an RL tracing agent that reduces vessel connected components from 29.1 to 4.5. Its reward could not be read.

I found **no** paper using exact β0/β1/χ (or a Betti-error) of the model's own output mask or polygon, a hole-aware multi-ring validity verifier, or a canonical ring hash as an RL reward. Explicit searches for "Euler characteristic reward" and "Betti reward" returned nothing relevant.

### Cited Findings

**Topo-R1: Detecting Topological Anomalies via Vision-Language Models** (Xu, Hu, Hu, Abousamra, Yu, Lyu, Qi, Samaras, Chao Chen; arXiv 2603.13054, 13 Mar 2026) [P]
- **Task:** the VLM is given a tubular segmentation mask (vessels, roads, neurons) and outputs boxes plus error type. The types are broken connection, spurious connection, missing branch and extra branch.
- **Reward:** "R_total = w_fmt·R_fmt + w_acc·R_acc + w_topo·R_topo, with w_fmt=0.10, w_acc=0.85, and w_topo=0.05".
  - R_acc = soft-F1 detection + localization φ(IoU) + type bonus, with type-aware Hungarian matching.
  - φ is a piecewise power map with thresholds {0.3, 0.5, 0.7, 0.9} → rewards {0.25, 0.55, 0.80, 1.0}.
  - Per matched pair: "R_topo(i,j) = (1 − clDice(M_i^gt, M_j^corr))·LocPen(b_j^pred)". Here M^corr is the **corrupted (input) mask** cropped to the region. The topology term therefore scores **where** the predicted box lands, not any geometry the model produced.
  - "LocPen … discourages oversized bounding boxes that trivially capture topological changes": τ_size = 0.3, λ = 0.8. This is an explicit anti-reward-hacking guard.
- **Betti numbers:** "each type is formally verifiable via changes in the Betti numbers (β0, β1), which underpins the automated verification step in our pipeline". Betti numbers are used for **data curation**; failed injections are retried up to 20 times. They are not used in the reward.
- **Setup:** GRPO with G=4 and KL β=0.05. F1@0.5 rises from 23.0% (SFT) to 45.2% on Qwen3-VL-4B. — [arXiv 2603.13054 html](https://arxiv.org/html/2603.13054v1)

**FloorplanVLM: A Vision-Language Model for Floorplan Vectorization** (Liu, Yang, Li, Yang; arXiv 2602.06507, Feb 2026) [P]
- **Reward:** "Validity Check (R_val): A binary reward ensuring the output is valid JSON and able to form a watertight closed polygon."
  - R_ext is the IoU of the external polygon. It sets a gate α: 0.1 if R_ext<0.3, a lerp up to 1.0 at 0.7.
  - R_int covers "F1 score and IoU metric of internal walls, openings, and room loops".
  - Weights: 0.1 / 0.5 / 0.4.
- **Metric definition:** validity rate means "syntactically valid JSON and geometrically watertight (i.e., forming closed polygons without topological gaps)".
- **What it does not specify:** the paper names no validity library or algorithm and has no hole, self-intersection or Betti terms.
- **Results:** SFT→GRPO raises validity from 90.20% to 96.10%, Ext-IoU from 0.8567 to 0.9252, and Room-IoU from 0.8521 to 0.8920. — [arXiv 2602.06507 html](https://arxiv.org/html/2602.06507)

**Mesh-RFT** (arXiv 2505.16761, NeurIPS 2025) [P]
- BER = \|E_∂M\|/\|E_M\| "quantifies the integrity of the mesh". The Topology Score measures quad-friendliness.
- Both build preference pairs and face-level masks for **Masked DPO**; this is preference RL, not GRPO.
- Reported gains: HD −24.6% and TS +3.8% over the pretrained model. — [arXiv html](https://arxiv.org/html/2505.16761v1)

**TraversRL** (arXiv 2607.17479, Jul 2026) [P]
- **Rewards:** the terminal reward is "buffered IoU" (buffer r=0.01 ≈ 1 m), and the local stepwise reward is r_t = ΔI_t/(ΔA_t+ε). Both are normalized within a rollout.
- Connectivity agreement improves, but as an **evaluation** outcome, not a reward term.
- The paper cites the RLVR literature. — [arXiv pdf](https://arxiv.org/pdf/2607.17479)

**Learning to Trace: Zero-Shot Connectivity Repair for Vessel Segmentation via Procedural RL** (Geshvadi, Durupinar; NeurIPS 2026) [A]
- A memory-augmented actor-critic traces vessels to reconnect fragments.
- On DRIVE it reduces connected components from 29.1 to 4.5 and lifts centerline recall from 85.8% to 87.7%. On ISBI-12 it goes from 10.8 to 4.2 components.
- **Reward not visible** (the OpenReview page was blocked). — [NeurIPS 2026](https://neurips.cc/virtual/2026/161763); [OpenReview](https://openreview.net/forum?id=NZuAkfqSsw)

**Other geometric-constraint rewards (not topology)** [P]
- GeoSVG-RL: rendering-validity gate, canvas overflow, text containment. — [arXiv html](https://arxiv.org/html/2605.25447v1)
- OP-HRG: part-in-object containment. — [arXiv html](https://arxiv.org/html/2607.15374)
- GenSeg-R1: negative-point validity, requiring points "outside ground-truth mask, at sufficient distance from the mask boundary". — [arXiv html](https://arxiv.org/html/2602.09701v1)
- LithoGRPO: non-differentiable shot-count reward. — [arXiv](https://arxiv.org/html/2606.00228)

**Topology losses, not RL (context)**
- TopoSculpt (2509.03938) uses a "Topological Integrity Betti (TIB) constraint" with persistent-homology curriculum refinement. It is a network and loss, not RL [A]. — [arXiv abs](https://arxiv.org/abs/2509.03938)
- Fast-χ, Betti matching and Topograph are in round 1 and round 2.

**Null results.** The searches "'Euler characteristic' reward reinforcement learning", "Betti number reward reinforcement learning segmentation", "persistent homology reward reinforcement learning", "GRPO segmentation reward connected components fragmentation", "VLM polygon output RL self-intersection holes" and "rasterize polygon reward GRPO hole multi-part" turned up no topology-invariant reward on model outputs. — [WebSearch, 12 Oct 2026; no single URL]

### Inferences
- **Claims the proposal must give up:**
  - "first topology-aware reward in VLM RL" (Topo-R1 holds it);
  - "first polygon-validity reward" (FloorplanVLM holds it);
  - "first automated topology score for RL fine-tuning of geometry generators" (Mesh-RFT holds it, in 3D and with DPO).
- **Claims that still look open:**
  1. **Exact topological invariants of the policy's own output** used as reward terms, compared against GT topology:
     - per-instance β0 and β1 matching;
     - hole recall and IoU computed per ring;
     - a χ-error reward.
  2. **A hole-aware, multi-ring, multi-part validity verifier** (ring simplicity, orientation, hole-in-shell, 4/8-connectivity consistency) as a graded reward, not a binary watertightness gate.
  3. **A canonical ring hash** as an exact-match or verifier reward, or for deduplicating rollouts in a group.
  4. **One batched GPU pass** that produces all these terms for G×B rollouts on-device.
- **Relationship to Topo-R1.**
  - Topo-R1 *detects* topological errors; the proposal would *penalize producing* them.
  - Topo-R1's clDice is a soft skeleton proxy; the proposal uses exact counts.
  - Topo-R1 already computes Betti numbers to verify synthetic errors. That makes it the obvious related-work anchor and a possible baseline, for example a "Topo-R1-style clDice reward" against a "β-exact reward".
- **Relationship to FloorplanVLM.** Its validity is binary. At 96% valid, a binary term is constant across most GRPO groups and stops producing advantage. A graded validity and topology vector could keep giving signal after saturation, for example the number of invalid rings, the hole-count error and the component-count error. This is an inference, not a measurement.

### Gaps
- **Which papers could still overlap, and why each is unresolved:**
  - Learning to Trace: reward not readable.
  - VectorLLM++: RLVR terms not readable.
  - FloorplanVLM: watertightness checker not stated.
  - Topo-R1: supplementary material (implementation) not read.
- **Coverage blind spots:**
  - Chinese-language venues and journals (ISPRS, TGRS, MICCAI 2026 proceedings) were not systematically searched.
  - Papers under review without arXiv versions are invisible to this search.

---

## Q3. Uncertainty- and ambiguity-aware rewards (annotator disagreement, soft IoU, band tolerance, conformal, multi-annotator)

### Takeaway
I found **no RL or RLVR paper whose segmentation reward uses multiple annotators, LIDC-style ambiguity, or conformal sets**. What exists is ad-hoc **tolerance shaping**:
- buffered IoU (TraversRL);
- tiered or piecewise IoU maps (SAM-R1, Topo-R1);
- IoU thresholds (SimpleSeg);
- coarse-to-fine fallbacks (CORE-Seg's GIoU-when-Dice<0.05).

SimpleSeg argues that RL itself handles the many-to-one nature of contours. Multi-reference scoring exists only in **evaluation** (Groundbench). This area is open.

### Cited Findings
- **TraversRL** rewards buffered IoU with buffer radius r = 0.01 (≈ 1 m). Its gains "widen as the buffer shrinks". — [arXiv 2607.17479](https://arxiv.org/pdf/2607.17479) [P]
- **SAM-R1** uses tiered IoU (4 / 3 / 2 / 0). — [arXiv 2505.22596](https://arxiv.org/html/2505.22596) (via r2)
- **Topo-R1** uses a piecewise-power φ(IoU), "ensuring continuous reward values for policy optimization". — [arXiv 2603.13054](https://arxiv.org/html/2603.13054v1) [P]
- **CORE-Seg** uses "r_mask = GIoU(B_M̂, B_M) if Dice < τ; 1 + λ·Dice otherwise" (τ=0.05, λ=2), to address "reward sparsity". — [arXiv 2603.05911](https://arxiv.org/html/2603.05911) [P]
- **SimpleSeg:** "we are not aiming to force the model to rigidly regress fixed ground-truth coordinates … as contour sequences are inherently flexible", and "Contours are inherently many-to-one w.r.t. masks; enforcing exact token matching is suboptimal." — [arXiv 2601.19228](https://arxiv.org/pdf/2601.19228) [P]
- **Groundbench** rescoring uses alternate reference masks, eight start phases and continuous-area intersection. This is evaluation only. — [arXiv 2609.27821](https://arxiv.org/pdf/2609.27821) (via r2)
- A search snippet attributes a reward that "aggregat[es] soft IoU, soft dice, and hard IoU" to an unidentified RL segmentation paper. — [S] (results included [AURORA 2508.02149](https://arxiv.org/pdf/2508.02149) and [CoPRS 2510.11173](https://arxiv.org/pdf/2510.11173); the source was not identified)
- **Non-RL ambiguity work** that could supply reward targets:
  - Multi-rater prompting. — [arXiv 2404.07580](https://arxiv.org/pdf/2404.07580)
  - Confidence Contours, which represents uncertainty as high- and low-confidence contours evaluated on LIDC. — [arXiv 2308.07528](https://arxiv.org/pdf/2308.07528)
  - Pionono, inter- and intra-observer variability modelling. — [ICCV 2023](https://openaccess.thecvf.com/content/ICCV2023/papers/Schmidt_Probabilistic_Modeling_of_Inter-_and_Intra-observer_Variability_in_Medical_Image_ICCV_2023_paper.pdf)

### Inferences
- **A topology-plus-ambiguity reward is a plausible differentiator.** Reward a rollout if its β0/β1 fall within the set of annotator topologies, or if its rings lie inside a high/low-confidence contour band (Confidence-Contours style). The engine gives ring-level geometry for band tests and exact β for set-membership tests. I found no prior art for this.
- **Band tolerance through the engine.** Morphological dilation and erosion of the GT, then ring containment, would generalize TraversRL's buffered IoU to masks with holes.

### Gaps
- No RL paper on LIDC-IDRI, QUBIQ or other multi-rater data was found. No conformal-prediction-as-reward paper was found for vision.

---

## Q4. Temporal-consistency rewards in video segmentation, grounding and tracking RL (2025–2026)

### Takeaway
Video RL rewards are **per-frame overlap averaged over frames, or temporal-interval IoU**. Mask propagation (SAM2 or XMem) supplies the temporal consistency, not the reward. No reward penalizes temporal *structural* changes, such as identity switches inside a mask, or components or holes appearing and disappearing between frames.

### Cited Findings
- **Veason-R1** (2508.11538) [P] — [arXiv html](https://arxiv.org/html/2508.11538)
  - Keyframe reward: R_k = S_{T_k}/max_t S_t, the foreground area ratio.
  - Spatial reward: Hungarian box IoU normalized by max(N_k, N_k^gt).
  - **Unified consistency reward:** R_u = (1/T̂)·Σ_t IoU(m_t, m_t^gt), with masks from SAM2 propagation starting at the chosen keyframe.
  - All weights are 1.0.
- **STVG-R1** (2602.11730) [P] — [arXiv html](https://arxiv.org/html/2602.11730v1)
  - R = temporal IoU of the interval + spatial-consistency indicator + format.
  - The indicator is 1 only if the predicted object ID equals the GT ID *and* appears within the predicted segment. It is "not computed as per-frame box IoU".
- **VideoSeg-R1** (2511.16077, AAAI) [S] — [arXiv](https://arxiv.org/html/2511.16077v1)
  - Rewards: reasoning format, segmentation format, mask IoU, point L1, box L1, non-redundant reasoning, and a negative-point distance reward.
  - XMem propagates the keyframe mask bidirectionally "ensuring temporal consistency without requiring reasoning computation on every frame".
- **VideoSEG-O3** (2606.06819, ICML 2026) [S] — [arXiv](https://arxiv.org/pdf/2606.06819)
  - Multi-turn RL with Progressive Mask Improvement, Keyframe Mask Quality and Temporal Accuracy rewards.
- **R1-Track** (2506.21980) [S] — [arXiv](https://arxiv.org/abs/2506.21980v1)
  - GRPO on Qwen2.5-VL-3B with a GIoU "soft reward"; 0.68 AO on GOT-10k.
- **TaRoS** (2511.19356), for video *generation* GRPO [A] — [arXiv abs](https://arxiv.org/abs/2511.19356)
  - Under sustained optimization "the reward score can lose fidelity as a proxy", with saturation and shortcuts.
  - Fix: component-level assessment and adaptive down-weighting of saturated components.

### Inferences
- **A temporal-topology reward appears open.** Examples: penalize Δβ0 and Δβ1 between consecutive frames (flicker of components or holes), or use ring-hash stability for static objects. The engine's per-frame invariants are cheap enough to compute on every propagated frame. No prior art was found.

### Gaps
- I did not open the VideoSeg-R1, VideoSEG-O3 or R1-Track full texts to confirm exact formulas.
- I did not search for video instance segmentation or multi-object tracking (MOT) RL with ID-switch rewards.

---

## Q5. RLVR theory relevant to geometry rewards: why verifiable rewards work, IoU-type reward hacking, group-relative ranking, dense vs sparse, multi-term weighting

### Takeaway
1. **What RLVR does.** It mostly sharpens what the base model can already sample, and its gains can be partly spurious on Qwen-family models. Any geometry-reward paper therefore needs random-reward and format-only controls.
2. **IoU-family rewards are demonstrably hackable or biased:**
   - mAP rewards encourage over-prediction;
   - S-measure alone gets hacked;
   - oversized boxes game region-based rewards;
   - rewards are dominated by large objects.
3. **GRPO standardizes rewards within a group.** A binary validity term goes silent once most rollouts pass it.
4. **Summed multi-term rewards lose resolution under GRPO** ("reward signal collapse"). Decoupled per-term normalization (GDPO) is the 2026 fix.
5. **Dense or smoothed rewards help with near-misses.**

### Cited Findings

**Why verifiable rewards work, and limits**
- VLM-R1: visual tasks have "well-defined ground-truth annotations", making them "naturally compatible with rule-based reward mechanisms". — [arXiv 2504.07615](https://arxiv.org/html/2504.07615) [P]
- Perception-R1: "perceptual perplexity determines RL superiority over SFT". RL helps most on counting and detection, and is comparable to SFT on grounding and OCR. — [arXiv 2504.07954](https://arxiv.org/html/2504.07954) [P]
- Yue et al.: "base models consistently surpass RLVR models … as k increases"; RLVR "does not elicit fundamentally new reasoning patterns". — [arXiv 2504.13837](https://arxiv.org/abs/2504.13837v4) [S]
- Counterpoint from Wen et al.: RLVR "implicitly incentivizes correct reasoning". — [arXiv 2506.14245](https://arxiv.org/html/2506.14245v2) [S]
- Spurious rewards: random rewards gave **+21.4 points** on MATH-500 for Qwen2.5-Math-7B, against +29.1 for ground truth. They "often fail to produce gains for other model families, such as Llama3 or OLMo2". — [arXiv 2506.10947](https://arxiv.org/abs/2506.10947) [A]
- Multimodal reward hacking (2607.09492): "outcome-only rewards cause severe hacking (up to 48.1% Reward Hacking Rate)", measured on safety VQA and chart VQA. It contains **no geometric or IoU cases**; a WebFetch summary claiming otherwise was false. — [arXiv 2607.09492](https://arxiv.org/pdf/2607.09492) [P, abstract text extracted]

**Reward hacking and bias with IoU-type rewards**
- VLM-R1: under an mAP reward "the model is incentivized to predict all categories … an instance of reward hacking". The fix is s = min(1, L_gt/L_pred). — [arXiv 2504.07615](https://arxiv.org/html/2504.07615) [P]
- Seg-R1: "using S-measure alone leads to reward hacking"; the fix is 0.7·IoU + 0.3·S. — [arXiv 2506.22624](https://arxiv.org/html/2506.22624) [P]
- Topo-R1: LocPen guards against "oversized bounding boxes that trivially capture topological changes". — [arXiv 2603.13054](https://arxiv.org/html/2603.13054v1) [P]
- Ground-R1: "training rewards are dominated by large visual regions", which leaves small regions with "consistently negative advantages". — [arXiv 2505.20272](https://arxiv.org/html/2505.20272) [P]
- GenSeg-R1: a no-target reward (weight 10) penalizes hallucinated boxes on negative queries. — [arXiv 2602.09701](https://arxiv.org/html/2602.09701v1) [P]
- ALTo: R = 𝕀_fmt + IoU − α·L; SimpleSeg found that hard length constraints degrade performance. — via r2 ([ALTo](https://arxiv.org/html/2505.16495v1); [SimpleSeg](https://www.alphaxiv.org/abs/2601.19228.md))

**Group-relative advantage**
- Topo-R1 writes the standard GRPO advantage as Â_i = (R_i − μ_G)/(σ_G + ε). — [arXiv 2603.13054](https://arxiv.org/html/2603.13054v1) [P]
- LithoGRPO notes that group normalization "cancels any constant offset". — via r2 ([arXiv 2606.00228](https://arxiv.org/html/2606.00228))
- Dr. GRPO identifies "an optimization bias in GRPO, which artificially increases response length (especially for incorrect outputs)". — [arXiv 2503.20783](https://arxiv.org/abs/2503.20783) [A]

**Dense vs sparse rewards**
- Smooth Operator (SNRA + AP-GRPO): "near-miss" samples are treated like failures, causing "advantage collapse"; a sigmoid-smoothed numeric reward plus absolute-preserving GRPO fixes it, on Numerical3D-50k. — [arXiv 2601.07695](https://arxiv.org/abs/2601.07695) [A]
- Group Revision (CVPR 2026): GRPO's "sparse reward, often criterion-induced, leads to minimal learning signals when all candidate responses fail". It uses potential-based shaping from revision improvements. — [CVF open access](https://openaccess.thecvf.com/content/CVPR2026/papers/Liu_From_Failure_to_Feedback_Group_Revision_Unlocks_Hard_Cases_in_CVPR_2026_paper.pdf) [P]
- Other dense or stepwise signals:
  - TraversRL: stepwise ΔI/ΔA. — [arXiv](https://arxiv.org/pdf/2607.17479) [P]
  - RSAgent: process reward from IoU improvement. — [S]
  - Point-VLM "Geometric Reward Credit Assignment": routes field-specific rewards to their token spans, plus a reprojection-consistency verifier. — [arXiv 2604.21160](https://arxiv.org/abs/2604.21160) [A]

**Multi-term weighting**
- GDPO (NVIDIA, 2601.05242): GRPO on a *sum* of rewards suffers "reward signal collapse", where distinct reward combinations yield identical advantages. GDPO normalizes each reward per group, sums, then batch-normalizes. — [arXiv 2601.05242](https://arxiv.org/pdf/2601.05242); [GitHub NVlabs/GDPO](https://github.com/NVlabs/GDPO) [S]
- Observed weightings and gating:

  | Paper | Weights or gating |
  |---|---|
  | Topo-R1 | 0.10 / 0.85 / 0.05 |
  | FloorplanVLM | 0.1 / 0.5 / 0.4, plus an IoU-gated α on the internal term |
  | GenSeg-R1 | 5 / 10 / 10 |
  | Veason-R1 | all 1.0 |
  | GeoSVG-RL | curriculum weighting, local constraints before global canvas fit |

  — sources as cited in Q1 [P]

### Inferences
- **Design rules for the proposal**, from the evidence above:
  1. Use **validity as a multiplicative gate**, as in FloorplanVLM, CAD executability and the GeoSVG render gate, rather than an additive term. Reason: a binary additive term is constant once most rollouts pass it.
  2. Use **graded topology errors**, such as |β0−β0*| and |β1−β1*| per instance, so groups keep variance after validity saturates.
  3. Normalize each term separately (GDPO-style) so the small topology term is not swamped by IoU. Topo-R1's topology term is only 5%.
  4. Add a size or area guard so components are not merged trivially to hit β0*. An over-dilated mask can achieve β0=1; this is analogous to LocPen.
- **What GRPO advantages are invariant to.**
  - They are invariant to a positive affine transform of the *total* reward within a group.
  - They are **not** invariant to monotone nonlinear maps (tiers, φ) or to re-weighting terms, because those change both rankings and spacings.
  - So "only ordering matters" is literally true only for rank-based advantage variants. This is my derivation, not a cited result.
- **Required controls.** The spurious-reward results mean a topology-reward paper on Qwen-VL needs a random-reward control and a format-only control. Gains should also be reported on Betti-error metrics that IoU rewards do not directly optimize.

### Gaps
- I found no paper measuring IoU-reward hacking specifically for mask-emitting (not box-emitting) policies.
- I found no theory on GRPO with lexicographic or gated rewards.
- GDPO details come from snippets only.

---

## Q6. Verdict: is "an exact topology and validity verifier as the RL reward for mask- or polygon-output models" already published?

### Takeaway
**Partially anticipated; the core claim is not published (moderate-to-high confidence).**

Prior art already covers:
- (a) a topology-flavoured term in a VLM GRPO reward: Topo-R1, clDice on input masks, weight 0.05, model outputs boxes;
- (b) a binary polygon-validity ("watertight closed polygon") reward in GRPO: FloorplanVLM;
- (c) automatic hole-sensitive topology scores for preference RL of 3D geometry: Mesh-RFT;
- (d) RL that improves connectivity outcomes: Learning to Trace, TraversRL;
- (e) non-differentiable geometric rewards in general: LithoGRPO shot count, GeoSVG containment and overflow, OP-HRG containment.

Not found anywhere:
- **exact β0/β1/χ of the policy's own output** (mask or polygon) as a reward;
- **hole-aware multi-ring validity** as a graded reward;
- **per-ring hole matching**;
- a **canonical ring hash** verifier;
- **one GPU engine** batching these for G×B rollouts.

The differentiation is real but narrower than "first topology reward". It should be framed as **exact, structural, output-side, hole-aware RLVR for spatial outputs**.

### Cited Findings
- **Topo-R1:**
  - "R_topo(i,j) = (1 − clDice(M_i^gt, M_j^corr))·LocPen(b_j^pred)", with w_topo = 0.05.
  - Betti numbers serve only for "the automated verification step" of synthetic errors. — [arXiv 2603.13054](https://arxiv.org/html/2603.13054v1) [P]
- **FloorplanVLM:** "A binary reward ensuring the output is valid JSON and able to form a watertight closed polygon"; validity 90.20% → 96.10%. — [arXiv 2602.06507](https://arxiv.org/html/2602.06507) [P]
- **Mesh-RFT:** BER = \|E_∂M\|/\|E_M\| plus TS, used for M-DPO. — [arXiv 2505.16761](https://arxiv.org/html/2505.16761v1) [P]
- **SimpleSeg** (closest polygon-output VLM RL): rewards are only IoU (with threshold τ), centroid MSE and format. — [arXiv 2601.19228](https://arxiv.org/pdf/2601.19228) [P]
- **Learning to Trace:** components 29.1 → 4.5; reward unknown. — [NeurIPS 2026](https://neurips.cc/virtual/2026/161763) [A]
- **VectorLLM++:** RLVR "that directly optimize geometric quality"; terms unknown. — [Research Square](https://www.researchsquare.com/article/rs-10949232/v1) [S]

### Inferences
- **Novelty claims, from safest to riskiest:**
  1. **Safe:** "Exact topological invariants (per-instance β0, β1, χ) and hole-aware ring validity of the model's own mask or polygon output, used as verifiable GRPO reward terms, computed in a single batched GPU pass." No prior art found.
  2. **Safe:** "A ring-hash and canonicalization verifier for polygon-emitting VLMs." It addresses SimpleSeg's invalid-polygon failure and Groundbench's legality collapse (7.5% legal at N=64, via r2). No prior art found.
  3. **Defensible, with citations required:** "Graded topology rewards beat binary validity (FloorplanVLM) and soft-skeleton proxies (Topo-R1's clDice)." This should be an *experimental* claim with those two as baselines.
  4. **Do not claim:** "first topology-aware reward in VLM RL" (Topo-R1); "first validity reward for polygon outputs" (FloorplanVLM); "first RL for connectivity" (Learning to Trace, vessel-tracing DRL).
- **Natural first experiments:**
  - SimpleSeg-style polygon RL on RefCOCO subsets rich in holes and multi-part objects;
  - Topo-R1's tubular datasets with a mask-*emitting* policy;
  - FloorplanVLM-style vectorization with hole-bearing rooms (courtyards).
- **Remaining uncertainty:** roughly a 20–30% chance that an unread paper already does β-exact output rewards. This is my subjective estimate. The likeliest candidates:
  1. VectorLLM++ (RLVR geometric terms unknown);
  2. Learning to Trace (its reward may count components);
  3. FloorplanVLM's unspecified watertightness checker (it might count rings);
  4. unindexed MICCAI 2026, ISPRS or TGRS papers on RL for vessel, road or building delineation.

### Gaps
- **Next fetches that would close the remaining uncertainty:**
  - the VectorLLM++ full text, from Research Square or an arXiv version;
  - the Learning to Trace paper (OpenReview NZuAkfqSsw), for its reward definition;
  - the FloorplanVLM code or supplement, for its validity checker;
  - the Topo-R1 supplement, for its Betti computation library and any β-based reward ablation.
- **Tooling:** the Exa MCP was unavailable for the whole session, and openaccess.thecvf.com and OpenReview blocked WebFetch (CVF was reachable via curl).
