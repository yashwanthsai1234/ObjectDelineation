# RL / reward fine-tuning of generative models whose outputs are geometry: which validity and geometric rewards they use, and what those rewards cost (state as of Oct 2026)

Scope and method notes:
- This is a gap-closing pass. It builds on `../GPU delineation use case map/r2_geometry_in_ml_loop.md` (Q1 and Q3) and `semiconductor_eda.md`. LithoGRPO, RLRF, SimpleSeg, ALTo and CurvyILT are already covered there; here they are only re-used or extended.
- "The engine" means the proposed exact GPU verifier. Its input is a binary mask. It outputs closed rings with holes, validity, β0, β1, χ, area, perimeter, vertex and jog counts, a rectangle-partition shot proxy, and a canonical hash.
- How the evidence was gathered:
  - About 28 tool calls: web searches, arXiv HTML/PDF fetches, and shallow clones of three reward-code repos (LithoGRPO, cadrille, floor-plan-rlvr). The clones were inspected directly, so the code claims below come from the code itself.
  - Exa MCP failed to connect, so only WebSearch and WebFetch were used.
  - Several PDFs came back compressed and the fetch summarizer could not read the numbers. Where that happened, it is flagged.

---

## Q1. Floorplan / layout generation with RL or reward-guided diffusion; building and urban layout RL

### Takeaway
RL with verifiable geometric rewards now exists for floorplans and layouts. All of it works in **vector or box space**, with LLM policies emitting JSON polygons or boxes, and every validity check is a **CPU shapely or box-arithmetic routine**.
- The clearest example is Mila's *Generative Floor Plan Design with LLMs via RLVR* (arXiv 2605.14117). Room-polygon overlap is a hard gate: overlap means reward 0.
- Its reward code silently "repairs" invalid polygons with `buffer(0)` and drops rooms that throw exceptions. That is a validity loophole an exact verifier would close.
- No paper reports reward compute as a bottleneck. Group sizes are small (G=4), polygon counts are small (5–8 rooms), and checks take milliseconds.
- Graphic layout RL (LaySPA, TextDiffuser-RL, GeoSVG-RL) uses box-level collision, alignment and containment rewards.
- Urban RL (CityPlanner) precomputes feasibility over OSM polygons, so no geometry is checked at reward time.

### Cited Findings
**Floorplans (LLM + GRPO)**
- **Floor-plan RLVR** (Lara et al., Mila; arXiv 2605.14117, May 2026), a two-stage SFT → GRPO pipeline.
  - Output: JSON with `spaces[]`, where each room has a `floor_polygon` (ordered x,y vertices in metres).
  - Rewards: a connectivity reward r_conn ∈ [0,1] (reconstructed adjacency graph against the input bubble diagram) and a total-area reward r_TA = max(0, 1 − |A − A*|/A*), weighted equally.
  - Hard gate: "Any candidate that fails to parse as valid JSON or contains overlapping room polygons receives a reward of 0."
  - Model and training: Llama-3.3-70B-Instruct, 4-bit, LoRA r=64. GRPO with G=4, KL 0.04, about 2 h per task on 6 nodes × 4×H100 for optimization plus 1 node for rollouts.
  - Inference: best-of-10, selecting minimum overlap area. The paper reports that RLVR "reduces overlap by 65% and compatibility issues by 56% relative to SFT" across 5–8-room tasks. Overlap still grows with room count (0.13 for 8 rooms). — [arXiv 2605.14117](https://arxiv.org/html/2605.14117v1); [repo](https://github.com/ludolara/floor-plan-rlvr)
- **What the floor-plan RLVR reward code actually does** (I read `src/grpo/grpo_evaluator.py` and `reward_calculator.py` in the repo):
  - Rooms are built as `shapely.Polygon`. `if not poly.is_valid: poly = poly.buffer(0)`, so self-intersecting or bow-tie rooms are **repaired rather than penalized**.
  - Rooms whose coordinates raise an exception are skipped (`except Exception: continue`), and zero-area rooms are dropped.
  - Overlap is O(n²) pairwise `intersection().area`. `is_overlap = round(overlap_ratio, 3) != 0`, so overlaps below about 0.05% of total area count as zero.
  - Rewards are computed in a sequential Python list comprehension.
  - No hole, gap or closure check exists. Nothing checks that rooms tile the footprint, that there are no slivers, or that the plan is a single connected outline. — [floor-plan-rlvr grpo_evaluator.py](https://github.com/ludolara/floor-plan-rlvr)
- Other floorplan generators found but **not RL**, titles and snippets only:
  - GRE-Diff (Gaussian room embeddings for layout diffusion) — [arXiv 2607.08086](https://arxiv.org/pdf/2607.08086)
  - DPLAN (minimal connectivity → floorplan) — [arXiv 2606.21159](https://arxiv.org/pdf/2606.21159)
  - FloorplanVLM (vectorization VLM) — [arXiv 2602.06507](https://arxiv.org/pdf/2602.06507)
  - A multi-agent DRL residential floor plan generator (health-oriented) — [ScienceDirect S266616592600092X](https://www.sciencedirect.com/science/article/pii/S266616592600092X)
  - Round 1 found that HouseDiffusion, GSDiff and DiffPlanner are vector-native with no RL. — [arXiv 2508.13738](https://arxiv.org/html/2508.13738v1)

**Graphic and design layout RL**
- **LaySPA** (arXiv 2509.16891, Sep 2025): GRPO with KL on Qwen-7B. It uses a hybrid reward with three parts:
  - geometric validity: inverse collision rate, alignment, boundary compliance;
  - structural fidelity: spacing consistency, distribution;
  - visual quality: underlay-text pairing.

  Format and IoU terms are added, with weights λ_format 0.1, λ_quality 0.8, λ_IoU 0.1. On CGL, collision rate falls from 0.692 to 0.431 and alignment rises from 0.319 to 0.597. No reward runtime is given. — [Emergent Mind summary of LaySPA](https://www.emergentmind.com/topics/layspa)
- **TextDiffuser-RL** (arXiv 2505.19291): a PPO layout stage before diffusion, with an overlap penalty between text boxes, an out-of-bounds penalty and an alignment reward. The fetch could not extract its runtime numbers. — [arXiv 2505.19291](https://arxiv.org/pdf/2505.19291)

**Urban and building layout RL**
- **CityPlanner** (arXiv 2609.09578, 2026) uses GRPO with an LLM agent on 5,670 OpenStreetMap instances for land allocation, road construction and station placement.
  - Reward: R = S(P) if Ψ(P) = 1, else R_fail.
  - Spatial units are "urban faces enclosed by roads" from polygonized road topology, and candidates are **precomputed** in `candidates.json` with feasibility flags. No geometry is computed at reward time.
  - Atomic-task RL raises feasibility to 100% (against 89.6% for flat GRPO) and improves objectives by 20.7%. — [arXiv 2609.09578](https://arxiv.org/html/2609.09578v1)
- **U-GRPO** residential block generation (2026) optimizes 3D block layouts for building energy and outdoor thermal comfort. Gains are 1.08–2.36% lower annual energy load and 0.37–0.56% more comfortable area in Shenzhen. The rewards are simulation-based, not validity-based. — [ScienceDirect S2666165926002115](https://www.sciencedirect.com/science/article/pii/S2666165926002115)
- Older urban-community DRL (Nature Computational Science 2023) and DQN building arrangement (Sci. Rep. 2026) also exist, titles and snippets only. — [Nature CS](https://www.nature.com/articles/s43588-023-00503-5); [Sci Rep](https://www.nature.com/articles/s41598-026-40788-w)

### Inferences
- **Floorplans are the clearest "validity reward is under-specified" case.** The checks people actually want are: rooms tile the footprint with no gaps (a hole or β1 check on the union), no slivers, every room a simple polygon (no `buffer(0)` laundering), one connected outline (β0 = 1), and adjacency derived from shared edges.
  - All of these fall out of rasterizing the plan into a label map and running the engine: per-label β0 and β1, union β1 for gaps, and ring adjacency for connectivity.
  - That would turn the current lenient shapely gate into an exact one. This is my inference.
- **Cost is not the argument here.** About 8 polygons × G=4 per prompt is microseconds to milliseconds in shapely, against a 70B LLM decode. The value is correctness and anti-reward-hacking, consistent with round 2's VLM-segmentation finding.

### Gaps
- I did not find a floorplan *diffusion* model fine-tuned with RL or reward guidance (DDPO, DRaFT and similar) using validity rewards. HouseDiffusion, GSDiff and GRE-Diff appear to be supervised only. The search was moderately thorough.
- The floor-plan RLVR paper does not report reward compute time.
- No source quantifies how often `buffer(0)` repairs fired during training, which would measure the loophole.

---

## Q2. CAD generation with RL (rewards = executability, geometric validity, Chamfer, IoU)

### Takeaway
CAD-code RL became standard in 2025–2026. cadrille (ICLR 2026 oral), CAD-Coder, ConCAD, PR-CAD and RA-CAD all follow the same recipe:
1. Execute the CadQuery code in an isolated subprocess with a timeout.
2. Tessellate the result into a mesh.
3. Compute IoU or Chamfer against the ground truth on CPU.
4. Add a hard invalidity penalty.

RL drives the invalidity ratio to about 0%. No paper reports reward compute as the bottleneck. The code does show that the **reward path is a fragile CPU sandbox**: cadrille needs a per-exec process with a 3 s timeout because of "CadQuery memory leaks", and IoU is computed by trimesh mesh booleans. That makes it the heaviest per-sample geometric reward outside lithography, though it is 3D and B-rep, not a 2D mask problem.

### Cited Findings
- **cadrille** (Kolodiazhnyi et al., arXiv 2505.22914; ICLR 2026 oral) is described as "the first to explore RL fine-tuning of LLMs for CAD tasks".
  - Reward: R = 10·IoU(generated CAD, GT mesh) + r_invalid, with r_invalid = −10 for invalid code and 0 otherwise.
  - Algorithm: Dr. CPPO (Dr. GRPO + CPPO) with G=16, keeping the top-|A| samples. Hard-example mining uses only examples with mean reward < 7.5. Hardware is 8×H100.
  - Results:

    | Dataset | IoU (SFT → RL) | Invalidity ratio (SFT → RL) |
    |---|---|---|
    | DeepCAD | 87.1 → 90.2% | 2.1 → 0.0% |
    | Fusion360 | 79.8 → 85.0% | 2.8 → 0.2% |
    | CC3D (real) | 61.8 → 67.9% | 5.9 → 0.2% |

  - Online RL beat offline DPO. — [arXiv 2505.22914](https://arxiv.org/html/2505.22914v2); [ICLR 2026 oral](https://www.iclr.cc/virtual/2026/oral/10006760)
- **cadrille reward and eval code** (I read `evaluate.py` in the repo):
  - IoU is computed by trimesh boolean `intersection` over every pair of `mesh.split()` components, with volumes summed.
  - Code runs via `exec()` and is tessellated with `compound.tessellate(0.001, 0.1)`.
  - Source comment: "To overcome CadQuery memory leaks, we call each exec() in a separate Process with timeout of 3 seconds", using a 16-process non-daemon pool. Chamfer uses 8,192 points.
  - These are the evaluation scripts; the RL training reward code was not in the released repo snapshot I inspected. — [GitHub col14m/cadrille](https://github.com/col14m/cadrille)
- **CAD-Coder** (Guan et al., arXiv 2505.19713): text → CadQuery, SFT then GRPO with chain of thought.
  - Rewards: a Chamfer-distance geometric reward plus a format reward.
  - The search-snippet text says non-executable or runtime-error outputs get a **−0.1 penalty**. My PDF fetch did not confirm the exact value or the group size, hardware or runtime. — [arXiv 2505.19713](https://arxiv.org/pdf/2505.19713)
- **ConCAD** (arXiv 2609.34494, Sep 2026): image → CadQuery with Qwen3-VL-2B.
  - Reward: R = 0.1·R_exec + 0.3·R_con + 1.0·R_iou, with three terms:
    - R_exec is a binary gate: executes and produces a non-degenerate solid.
    - R_con is regex pattern-matching over 8 engineering-constraint categories with **no kernel execution**.
    - R_iou is volumetric IoU by Monte Carlo with **8,192 points** after inertia-aligned normalization.
  - Programs run "in isolated subprocesses with strict timeouts" (OpenCASCADE → triangle mesh).
  - Training: GRPO with G=4, 800 updates, on a **single A800**.
  - DeepCAD results: IoU 0.717 → 0.734, CD 0.011 → 0.009, invalidity ratio 1% → 0%.
  - New metric: G-CSR (B-rep constraint satisfaction, analytically extracted). — [arXiv 2609.34494](https://arxiv.org/html/2609.34494)
- **PR-CAD** (arXiv 2604.19773): progressive-refinement text-to-CAD. A search snippet reports invalidity ratio 0.62 and mean CD 5.87 against GPT-4o; its RL details were not fetched. **RA-CAD** (arXiv 2608.05714) learns a post-execution critique, title and snippet only. — [arXiv 2604.19773](https://arxiv.org/html/2604.19773); [arXiv 2608.05714](https://arxiv.org/html/2608.05714v1)
- **CADFusion** (MSRA; arXiv 2501.19054) alternates sequential learning with a visual-feedback stage.
  - It renders CAD and scores the renders with an LVLM, then uses **DPO instead of online RL** because "conventional RL is computationally expensive".
  - LVM score is 8.96 against 2.01 for Text2CAD. — [arXiv 2501.19054](https://arxiv.org/html/2501.19054v2)

### Inferences
- **Reward-cost profile:**
  1. Subprocess spawn plus CadQuery/OCC execution (up to 3 s timeout)
  2. Tessellation
  3. Mesh boolean IoU or point sampling

  This is plausibly 0.1–3 s per rollout on CPU, which at G=16 on 8×H100 could rival decode time. No paper measures it, so this is an inference only. CADFusion's choice of offline DPO over online RL "because RL is expensive" is weak indirect evidence that render-and-score cost matters.
- **Fit for the engine:** limited and indirect. CAD rewards are 3D B-rep or mesh. The 2D analogues where the engine fits are:
  - (a) sketch- or profile-level validity: closed loops, holes, self-intersection of each extrude profile, by rasterizing the sketch and re-polygonizing;
  - (b) slice-wise or voxel IoU and topology (β0/β1 per slice, or a 3D Euler number) as a cheaper, deterministic stand-in for trimesh booleans, which fail on non-watertight meshes and return `None` (cadrille's `except: pass`).

  The second point is an inference worth testing. Voxelized IoU on GPU is standard, but exact per-slice topology is not used in any CAD RL reward I found.

### Gaps
- No CAD-RL paper reports per-rollout reward latency or the share of step time spent in CadQuery execution.
- I found no CAD-RL reward with explicit **topological validity** (genus, β1, watertightness) terms beyond "executes and is non-degenerate".
- I did not find papers literally titled "CAD-RL", "CADCrafter RL" or "Text2CAD-RL". CADCrafter was not checked.

---

## Q3. SVG / vector-graphics RL: reward terms and what they cost (rasterization, CLIP, DINO)

### Takeaway
SVG RL is crowded:
- RLRF (NeurIPS 2025)
- Reason-SVG
- CTRL-S, the multi-task, multi-reward SVG paper
- GeoSVG-RL
- RULER
- IntroSVG
- RefineSVG
- Render-in-the-Loop
- GVR-Coder

Rewards follow one template: **render (CairoSVG or headless Chromium) → pixel, edge or embedding similarity (L2, Canny-L2, DreamSim, CLIP or Long-CLIP, DINOv2, HPSv2, VLM judge) → length penalty → renderability gate**.

The only explicitly *geometric* verifier is **GeoSVG-RL's browser-backed checker**: bbox containment, canvas overflow, anchor endpoints within 12 px, text-in-box, and edge-connectivity F1. It runs at **about 18 SVGs/s on 8 CPU workers, with no GPU**, roughly 0.44 s per SVG per worker. That is the only reward-throughput number I found in SVG RL. No SVG paper says reward compute is *the* bottleneck. RLRF's 2,048 renders per step over 3 days on 32 H100s suggests that rendering is small next to decoding (my arithmetic).

### Cited Findings
- **RLRF** (Rodriguez et al., NeurIPS 2025; arXiv 2505.20793).
  - Im2SVG rewards: normalized L2, L2-Canny, DreamSim, DreamSim-Canny, and SVG length deviation.
  - Text2SVG rewards: CLIP plus a Qwen2.5-VL-72B judge, with text elements stripped before rendering to prevent reward hacking.
  - Renderer: CairoSVG, rendering at the reference image's size and aspect ratio "to prevent the model from exploiting small viewbox dimensions".
  - Scale: 64 rollouts × 32 images = 2,048 renders per step, 500 steps on 20k samples, about 3 days on 4 nodes × 8×H100. Handling of invalid SVGs is not specified. — [alphaxiv 2505.20793](https://www.alphaxiv.org/abs/2505.20793.md); [NeurIPS 2025](https://papers.neurips.cc/paper_files/paper/2025/hash/57126328c3b40cf618a34f1c5df24d8a-Abstract-Conference.html)
- **Reason-SVG** (arXiv 2505.24499): GRPO on Qwen2.5-VL-7B with G=8, 8,000 steps on 2,000 prompts, **32×H800**.
  - Hybrid reward with four terms:

    | Term | Weight | How it is scored |
    |---|---|---|
    | Think (reasoning presence) | 0.1 | — |
    | Render (structural validity) | 0.1 | CairoSVG |
    | Semantic | 0.6 | CLIP ViT-L/14 |
    | Aesthetic | 0.2 | HPSv2 |

  - SVG validity is 99.8%, with about 12 s inference per SVG. — [arXiv 2505.24499](https://arxiv.org/html/2505.24499v2)
- **CTRL-S / "Reliable Reasoning in SVG-LLMs via Multi-Task Multi-Reward RL"** (arXiv 2603.16189).
  - Reward: R_total = R_format · (2·R_DINOv2 + 1·R_Long-CLIP + 1·R_eff).
  - R_format is 1 only if the think block is present and the SVG renders in CairoSVG; R_eff = 1 − ((1/L_gt)·max(0, L_gen − L_gt/2))².
  - Training: GRPO via verl, 16 responses per prompt, 32 GPUs, about 12 h.
  - Results: success rate 99.85%, tokens 1.1k → 346. — [arXiv 2603.16189](https://arxiv.org/html/2603.16189)
- **GeoSVG-RL** (Li, Cai, Chen, Wang; arXiv 2605.25447): layout-constrained text-to-SVG diagrams.
  - The verifier has two paths:
    - Rendering path: headless Chromium at 800×600, with text boxes via `getBBox()` and `getBoundingClientRect()`, and a 5 s render timeout.
    - Structure path: lxml parse into global-coordinate primitives.
  - Six reward families:
    - R_exec gate;
    - canvas fit and overflow area fraction;
    - anchor accuracy (endpoints within τ=12 px) and anchor error normalized by node-box diagonal;
    - text-in-box and padding (< 6 px flagged);
    - edge-connectivity F1;
    - code cleanliness.
  - **Throughput: "approximately 18 SVGs per second using 8 parallel CPU workers. No GPU acceleration."**
  - Training: GRPO with G=4 and 1,500 updates on A100s, 96 GPU-h SFT plus 128 GPU-h GRPO.
  - Results against SFT: anchor accuracy 78.6% vs 43.2%, text-in-box 83.0% vs 39.8%, edge F1 90.4 vs 47.3, render success 92.3%. — [arXiv 2605.25447](https://arxiv.org/html/2605.25447)
- Further 2026 SVG-RL papers found, titles and snippets only:
  - **RULER**: instance-aware rubric rewards — [arXiv 2609.25270](https://arxiv.org/html/2609.25270)
  - **IntroSVG**: generator–critic with rendering feedback — [arXiv 2603.09312](https://arxiv.org/html/2603.09312)
  - **Render-in-the-Loop**: visual self-feedback — [arXiv 2604.20730](https://arxiv.org/html/2604.20730v1)
  - **RefineSVG**: visual-feedback RL for image-to-SVG — [Pith 2607.27699](https://pith.science/paper/2607.27699)
  - **GVR-Coder**: document and meeting diagrams — [arXiv 2607.28073](https://arxiv.org/pdf/2607.28073)

### Inferences
- **Cost arithmetic (mine, not measured by the papers):**
  - RLRF needed about 1.02M renders over about 72 h, roughly 4 renders/s, which is trivial for CPU CairoSVG. The real reward cost is the GPU embedding models (DreamSim, CLIP, DINOv2, HPSv2) and, for Text2SVG, a **72B VLM judge**, which is likely the dominant reward cost.
  - GeoSVG-RL's 18 SVGs/s CPU verifier is the one place where a geometric verifier could plausibly limit throughput. If each update scores about 128–512 candidates, verification takes about 7–28 s per update, against unknown decode time. Batch semantics are ambiguous in the paper, so treat this as unconfirmed.
- **Fit for the engine:** SVG rewards are render-and-compare in pixel or embedding space. Nobody checks **path-level topology**: whether filled regions have the right number of holes and components, whether the rendered fill is the intended shape (even-odd versus nonzero artefacts), or a vertex and segment budget against the minimum needed.
  - The engine can provide all of these from the *rendered* raster: render once, polygonize, then compare β0/β1, area and perimeter, and ring count with the target's polygonized raster.
  - That adds a cheap, exact "structural" reward beside DINO and CLIP, computed on the same GPU raster. It does not replace rendering; it reuses the render.
- A GPU rasterizer (DiffVG, or GPU path rasterization) plus the engine would make the whole reward GPU-resident. Today CairoSVG and Chromium force a CPU round-trip.

### Gaps
- No SVG-RL paper reports CairoSVG render time per sample or the reward share of step time.
- OmniSVG and StarVector have no published RL variant that I could confirm. StarVector appears as RLRF's base model; OmniSVG appears only as a baseline.
- Handling of invalid SVGs in RLRF is unspecified.

---

## Q4. Chip design: RL for mask synthesis, OPC and ILT beyond LithoGRPO; macro placement (AlphaChip and its critics); routing RL; PCB RL — geometric terms and costs

### Takeaway
- **Mask synthesis is still the one domain where geometric reward cost is a documented bottleneck.** LithoGRPO's released shot counter is a CPU Numba-plus-PuLP/CBC ILP with a 300 s solver time limit, and its README says only "up to ~130×" faster.
- NVIDIA's generative-RL ILT puts a **GPU litho simulator plus 100-step low-resolution ILT** inside the reward. The authors call it "computationally heavier than purely learned proxies" but give no per-sample time. Its reward is EPE only, and MRC is enforced by morphological opening before evaluation, not rewarded.
- **Placement** has moved away from RL toward diffusion and flow models with *guidance or projection* for legality: ChipDiffusion (ICML 2025), FlowPlace (DAC 2026). AlphaChip's proxy cost (wirelength + density + congestion on a grid) is cheap, but critics show it correlates poorly with PPA.
- **PCB RL** is now benchmarked through KiCad DRC (PCBWorld, KDD 2026 workshop), but per-step verifier runtimes could not be extracted.

### Cited Findings
**Mask synthesis, OPC and ILT**
- **LithoGRPO released code** (I inspected the repo).
  - Only the shot counter is released so far; GRPO, the flow prior and the litho simulator are "coming soon".
  - `FastShotCounter` works in three steps:
    1. Row RLE, then enumerate maximal all-ones rectangles.
    2. Prune contained rectangles with a Numba `@njit` greedy pass.
    3. Solve minimum set cover with **PuLP CBC** (`PULP_CBC_CMD(timeLimit=300)`).
  - The reference counter is LithoBench's **adabox** decomposition (vendored, MIT).
  - The README claims "Spearman ρ ≈ 1.0" against the reference and speedups "up to ~130× in our paper's setting… the exact factor depends on hardware, mask resolution, and the ILP solver". The paper reports 134.6–491.3×, 0.16–0.30 s per 512² mask, from round 2.
  - Dependencies are torch, numba, pulp and pandas, and the code is CPU-bound. — [GitHub laiyao1/LithoGRPO](https://github.com/laiyao1/LithoGRPO); [arXiv 2606.00228](https://arxiv.org/html/2606.00228)
- **Generative RL ILT** (Yang and Ren, NVIDIA; arXiv 2602.19027).
  - The reward R_k = −EPE(M_k^ILT, Z) is computed **after** a short ILT refinement: 8× downsample, 100 ILT iterations, upsample, binarize at 0.5.
  - Training: K=16 proposals per design, 8 designs per step, on 1 DGX (8×A100).
  - The authors say embedding ILT in the reward "is computationally heavier than purely learned proxies" but "GPU-accelerated implementations and our low-resolution, few-step refinement loop keep wall-clock overhead modest". **No per-sample runtime is given.**
  - The interface "naturally supports… mask complexity/shot count, and manufacturability constraints", but the reported results use only EPE and PVB.
  - Curvilinear MRC compliance is enforced by "morphological opening and cleanup" before printability evaluation. — [arXiv 2602.19027](https://arxiv.org/html/2602.19027)
- **RL-OPC** (Bei Yu group; TCAD 2024): a DRL agent moves segments, reducing EPE by 43% on average and beating Calibre by 24.5% EPE and 2% PV band. — [RL-OPC PDF](http://www.cse.cuhk.edu.hk/~byu/papers/J97-TCAD2024-RL-OPC.pdf)

**Macro placement**
- **AlphaChip / Circuit Training critiques** (TILOS MacroPlacement).
  - The proxy cost is wirelength + density + congestion, with congestion weight 0.5 in some runs. It is computed on a grid of cells by the `Plc_client`, with no routing.
  - Main finding: when CT trains to convergence, "proxy cost improves significantly, but does not correlate with timing metrics".
  - CT runs took about 36 h (200 iterations on 8×V100).
  - In March 2025 the testcases were updated with AlphaChip and improved simulated-annealing baselines. — [TILOS MacroPlacement](https://github.com/TILOS-AI-Institute/MacroPlacement); [Updated Assessment arXiv 2302.11014](https://arxiv.org/pdf/2302.11014)
- Further context from search snippets:
  - CT from scratch takes 30–80 h on 8 GPUs, against under 2 h for CMP or AutoDMP.
  - The **2026 macro placement challenge** caps runtime at 1 h end to end per benchmark. — [partcleda/macro-place-challenge-2026](https://github.com/partcleda/macro-place-challenge-2026); [Vighnesh Iyer saga summary](https://vighneshiyer.com/research/eda-cad-vlsi/machine-learning-for-macro-placement-alphachip-the-saga/)
  - Google's rebuttal is "That Chip Has Sailed". — [arXiv 2411.10053](https://arxiv.org/pdf/2411.10053); [Wikipedia: AlphaChip controversy](https://en.wikipedia.org/wiki/AlphaChip_(controversy))
- **ChipDiffusion** (Lee et al., ICML 2025) is a diffusion placer pretrained on synthetic netlists. It uses "guided sampling in lieu of RL" to optimize **legality (no overlap, in bounds) and HPWL**. — [arXiv 2407.12282](https://arxiv.org/html/2407.12282v3); [PMLR v267](https://proceedings.mlr.press/v267/lee25y.html)
- **FlowPlace** (DAC 2026; arXiv 2604.23658) is flow matching **without RL or reward fine-tuning**.
  - Legality comes from a hard-constraint projection, greedy GPU-parallel overlap legalization at each sampling step.
  - It reports zero overlaps against ChipDiffusion's 4.63%, uses 20–50 steps against about 1000, is 10–50× faster, and runs in seconds against "hours of rollout exploration" for RL placers. — [arXiv 2604.23658](https://arxiv.org/html/2604.23658)
- Other placement RL work noted, titles only:
  - "Expertise can be helpful for RL" (ICLR 2026) criticizes RL placers' "oversimplified objectives such as macro HPWL and rectangular uniform wire density". — [ICLR 2026 PDF](https://proceedings.iclr.cc/paper_files/paper/2026/file/cbc1ad2066f0afebbcea930c5688fc1f-Paper-Conference.pdf)
  - MaskPlace (pixel-mask RL placement) — [arXiv 2211.13382](https://arxiv.org/pdf/2211.13382)
  - PlaceReasoner-Beta — [arXiv 2609.21263](https://arxiv.org/pdf/2609.21263)

**Routing and PCB**
- **PCBWorld** (Song et al., KDD 2026 Workshop on Evaluation and Trustworthiness of Agentic AI; arXiv 2607.05915) is an engine-grounded PCB routing environment on **KiCad**.
  - Its signals include DRC violations, clearance, overlap, wirelength and routability.
  - The PDF streams were compressed, so per-step verifier runtime and baseline numbers **could not be extracted**. — [arXiv 2607.05915](https://arxiv.org/pdf/2607.05915)
- Other routing and PCB RL, titles and snippets only:
  - Offline-RL detailed-routing convergence — [arXiv 2512.03594](https://arxiv.org/pdf/2512.03594)
  - Transformer-RL net ordering for detailed routing (IJCAI 2025), via search snippet
  - DRL dummy-pad placement for PCB electroplating — [ScienceDirect S0957417425022584](https://www.sciencedirect.com/science/article/abs/pii/S0957417425022584)
  - MIT thesis on DRL PCB placement — [MIT DSpace](https://dspace.mit.edu/bitstream/handle/1721.1/139247/Crocker-pcrocker-meng-eecs-2021-thesis.pdf?sequence=1&isAllowed=y)
  - DeepPCB's description of routing reward design: DRC penalties, via cost — [DeepPCB blog](https://deeppcb.ai/reinforcement-learning-pcb-routing-explained/)

### Inferences
- **Mask synthesis remains the beachhead and is getting sharper.** LithoGRPO's released counter is CPU ILP (Numba + CBC) on the exact problem the engine's rectangle-partition proxy targets. NVIDIA's ILT-RL already says shot count and MRC *should* go into the reward but uses only EPE, plausibly because of cost. That is an inference; the paper does not say why.
  - A GPU, batched, rank-faithful shot/vertex/jog/hole proxy is directly substitutable. Its success criterion is set by LithoGRPO: Spearman ρ ≈ 1 against adabox, at well under 0.16 s per 512² mask.
- **Placement is a weak fit.** Objects are rectangles; legality is box overlap, cheap on GPU (FlowPlace already does it); and cost bottlenecks lie in PPA correlation, not geometry.
- **PCB and routing are a plausible but unverified fit.** KiCad DRC is a CPU process, and an RL loop calling it per step could be slow. Rasterize copper per layer and polygonize with the engine to get clearance (via distance or erosion), shorts (β0 merges between nets) and opens (β0 splits within a net) as GPU rewards. This is unmeasured.

### Gaps
- No public per-sample runtime exists for NVIDIA's ILT-in-the-loop reward.
- PCBWorld verifier timings and baseline results were unreadable from the PDF.
- I did not find mask-synthesis RL papers other than LithoGRPO, NVIDIA GRL-ILT and RL-OPC that reward shot count, vertex count or MRC. Whether others exist from 2026 DAC/ICCAD is unknown, since ICCAD 2026 proceedings are likely not yet indexed.
- No routing-RL paper reports DRC-check cost as a bottleneck.

---

## Q5. Map, road-network and HD-map generation with RL

### Takeaway
RL here is sparse. HD maps are mostly used as **reward sources for driving policies**, such as solid-line-crossing and crosswalk legality rewards in closed-loop GRPO of a 10B driving model, not as RL-generated outputs. The one generative example is **TraversRL** (pedestrian pathway polylines, arXiv 2607.17479). Its rewards are vector-based (connectivity, topology, obstacle overlap, IoU), but the fetch could not extract specifics. No road-network generator with GRPO and topology rewards was found.

### Cited Findings
- **TraversRL** (Han, Wolfe, Howe; arXiv 2607.17479) generates traversable pedestrian-pathway **polylines** with RL.
  - The fetch summary lists connectivity, topology or junction validity, traversability, obstacle-overlap and IoU terms, computed on vectors rather than rasters.
  - The algorithm, numbers and runtimes were not extractable from the compressed PDF, so treat the reward list as low-confidence. — [arXiv 2607.17479](https://arxiv.org/pdf/2607.17479)
- **"What HD-Map Rewards Buy"** (SJSU thesis) does closed-loop GRPO of the 10B Alpamayo 1.5 driving model in Argoverse 2 3DGS reconstructions.
  - It adds a route-corridor lane-keeping reward and map-derived legality rewards: a solid-line crossing penalty and crosswalk compliance.
  - These are geometric point-in-polygon or segment-crossing tests against vector maps. — [SJSU ScholarWorks](https://scholarworks.sjsu.edu/etd_theses/5827/)
- Lane-topology HD-mapping methods are supervised, with no RL found. — [Lane topology survey arXiv 2504.01989](https://arxiv.org/pdf/2504.01989)

### Inferences
- The map-legality rewards in driving RL are point, segment and polygon predicates over a static vector map. They are cheap and do not need mask polygonization. The engine's role would be indirect: producing clean, hole-aware drivable-area and crosswalk polygons *from* segmentation maps that then serve as the reward map.
- For road-network *extraction* with RL, topology rewards (β0 of the road mask, loops = β1, APLS-like connectivity) are a natural fit for the engine. I found no paper using them as RL rewards.

### Gaps
- I found no GRPO or DDPO-style RL for HD-map or road-graph *generation* with topology or validity rewards. The search was moderate.
- TraversRL's exact reward formulas and costs were not verified.

---

## Q6. Reward libraries and "verifier suites" for spatial and geometric outputs; their runtimes; is reward compute a bottleneck anywhere besides lithography?

### Takeaway
**There is no shared geometric reward library.** Every paper ships bespoke verifiers:
- shapely pairwise overlaps (floorplans)
- subprocess CadQuery plus trimesh booleans (CAD)
- CairoSVG or Chromium plus CLIP, DINO or DreamSim (SVG)
- Numba plus PuLP ILP (lithography)
- KiCad DRC (PCB)
- precomputed feasibility tables (urban)

None is GPU-native for geometry. The only throughput figures found are LithoGRPO's (34–148 s → 0.16–0.30 s per mask) and GeoSVG-RL's (about 18 SVGs/s on 8 CPU workers).

General RL-systems work (ROLLART) shows reward or verifier stages *can* dominate when rollouts are fast: reward-GPU utilization rose from 6% to 88% and per-step rollout time halved after disaggregating reward. Apart from that, no geometry paper besides lithography *reports* reward compute as the bottleneck. CAD (sandboxed kernel exec, 3 s timeouts) and the SVG VLM judges are the likeliest hidden ones.

### Cited Findings
- **Spatial verifiers in image generation:** SpatialReward (CVPR 2026) is a verifiable spatial reward model for T2I. It uses a prompt decomposer, expert detectors for grounding and VLM chain-of-thought over relations, with SpatRelBench. It is detector- and VLM-based, not exact geometry. — [arXiv 2603.22228](https://arxiv.org/abs/2603.22228); [CVPR 2026](https://openaccess.thecvf.com/content/CVPR2026/papers/Zhou_SpatialReward_Verifiable_Spatial_Reward_Modeling_for_Fine-Grained_Spatial_Consistency_in_CVPR_2026_paper.pdf)
- **GeoSVG-RL's verifier** is the closest thing to a "geometry unit-test suite" in RL: six reward families, Chromium plus lxml, 5 s timeout, **about 18 SVGs/s on 8 CPU workers, no GPU**. — [arXiv 2605.25447](https://arxiv.org/html/2605.25447)
- **Groundbench's polygon legality audit** (distinct in-range vertices, simple positive-area polygon, no self-intersection) is an *evaluation* suite, not used in training. This is from round 2. — [arXiv 2609.27821](https://arxiv.org/pdf/2609.27821)
- **Code-level verifiers:**
  - floor-plan-rlvr: shapely with `buffer(0)` repair. — [repo](https://github.com/ludolara/floor-plan-rlvr)
  - cadrille: per-exec subprocess, 3 s timeout, 16-proc pool, trimesh booleans. — [repo](https://github.com/col14m/cadrille)
  - LithoGRPO: Numba plus CBC ILP with timeLimit=300. — [repo](https://github.com/laiyao1/LithoGRPO)
- **Reward-stage cost in RL systems generally:** ROLLART (arXiv 2512.22560) offloads stateless reward computation, which spans reward models, verifiers and sandboxes. It reports "reward-GPU utilization increased from 6% to 88%" and "per-step rollout time halved". This concerns LLM agentic RL, not geometry. — [arXiv 2512.22560](https://arxiv.org/pdf/2512.22560)
- Geometric-reasoning RLVR datasets (GeoSym127K, geometry diagram parsing) use symbolic verifiers for *math* geometry, not spatial outputs, so they are out of scope. — [arXiv 2605.16371](https://arxiv.org/pdf/2605.16371)

### Inferences
- **Reward-compute bottleneck by domain:**

  | Domain | Typical reward | Where it runs | Reported bottleneck? | Likely hidden cost |
  |---|---|---|---|---|
  | Lithography / ILT | L2, PVB, EPE, **shot count** | GPU sim + **CPU ILP** | **Yes** (LithoGRPO: 34–148 s → 0.16–0.30 s) | NVIDIA ILT-in-reward: "heavier", not quantified |
  | CAD code | exec gate + IoU/Chamfer | CPU subprocess, OCC, trimesh | No | Sandbox exec ≤3 s timeout, mesh booleans |
  | SVG | render + DINO/CLIP/DreamSim/HPSv2/VLM judge | CPU render + GPU encoders; 72B judge | No | VLM judge; Chromium verifier at 18/s per 8 workers |
  | Floorplan / layout | overlap gate, adjacency, area | CPU shapely | No | Negligible |
  | Urban | precomputed feasibility | lookup | No | n/a |
  | Placement | HPWL, density, congestion proxy | GPU/CPU grid | No (problem is proxy fidelity) | n/a |
  | PCB | KiCad DRC | CPU | Unknown | KiCad DRC per step (unmeasured) |
  | VLM segmentation (round 2) | mask/box IoU | CPU/GPU | No | n/a |

- **Positioning the engine as a reusable reward library.** The defensible pitch is **correctness, uniformity and anti-hacking first, speed second**, with lithography as the one domain where speed is proven to matter.
  - Every domain above re-implements validity checks inconsistently: shapely repair-then-check, regex constraint checks, a `try/except` that returns `None` IoU.
  - None checks holes, components or closure exactly.
  - One GPU verifier with a stable API (`validity`, `β0`, `β1`, `χ`, `area`, `perimeter`, `n_vertices`, `jogs`, `shot_proxy`, `hash`) applied to *rendered or rasterized* rollouts would give the same canonical checks across floorplans, SVG, CAD sketch profiles, PCB copper and masks.
  - The canonical hash also enables **reward caching and duplicate-rollout detection** within GRPO groups, which no paper does.

### Gaps
- I found no published "geometry verifier suite" (a shapely- or GEOS-based library) shared across RL projects.
- No paper measures per-sample shapely or GEOS validity-check time inside RL.
- The CAD sandbox and SVG VLM-judge cost shares are unmeasured, so whether they are bottlenecks is an inference.

---

## Q7. Diffusion-model reward fine-tuning (DDPO, DPOK, ReFL, DRaFT, Flow-GRPO, DanceGRPO) applied to segmentation masks, layouts or geometry with geometric rewards

### Takeaway
Diffusion and flow RL with **geometric** rewards on **geometric outputs** is still rare. **LithoGRPO** (flow matching + GRPO on binary masks, with shot-count reward) remains the main published instance.
- Flow-GRPO and DanceGRPO target T2I, video and compositional rewards from detectors or VLMs.
- Chip placers (ChipDiffusion, FlowPlace) use **guidance or projection**, not reward fine-tuning.
- Floorplan diffusion models are not RL-tuned.
- Medical work (Compositional Reward Models, GDPO) uses segmentation-verifier sub-rewards on *images*, not on masks as outputs.

Fine-tuning a mask-generating diffusion or flow model with exact topology rewards is an open slot.

### Cited Findings
- **Flow-GRPO** (arXiv 2505.05470) is the first online policy-gradient RL for flow matching, via ODE→SDE conversion. Its rewards are compositional T2I (GenEval: detector-based counting and position), OCR and human preference. — [arXiv 2505.05470](https://arxiv.org/pdf/2505.05470)
- **DanceGRPO** applies GRPO to visual generation (image and video). — [arXiv 2505.07818](https://arxiv.org/pdf/2505.07818)
- **DDPO** frames denoising as an MDP with REINFORCE; **DRaFT** back-propagates differentiable rewards; **Diffusion-DPO** handles preferences. — [DDPO](https://rl-diffusion.github.io/files/paper.pdf); [DRaFT arXiv 2309.17400](https://arxiv.org/pdf/2309.17400); [Diffusion-DPO CVPR 2024](https://openaccess.thecvf.com/content/CVPR2024/papers/Wallace_Diffusion_Model_Alignment_Using_Direct_Preference_Optimization_CVPR_2024_paper.pdf)
- **SPO** notes that early denoising steps set layout, so layout-dispreferred samples should penalize early steps. This is relevant to credit assignment for geometric rewards. — [arXiv 2406.04314](https://arxiv.org/html/2406.04314v1)
- **LithoGRPO**: flow matching + GRPO on 512² binary masks, with reward L2 + PVB + EPE + shot count, G=6, 1000 steps, 4×RTX 3090. This is from round 2. — [arXiv 2606.00228](https://arxiv.org/html/2606.00228)
- **ChipDiffusion** uses guided sampling for legality and HPWL "in lieu of RL"; **FlowPlace** uses hard-constraint projection, with no RL. — [arXiv 2407.12282](https://arxiv.org/html/2407.12282v3); [arXiv 2604.23658](https://arxiv.org/html/2604.23658)
- **Compositional Reward Models for Conditional Medical Image Generation** (arXiv 2609.05028) decomposes quality into "verifier grounded sub rewards", including segmentation, measurement and classification, and optimizes conditional diffusion with GDPO. This is a snippet only; whether any sub-reward is topological is unknown. — [arXiv 2609.05028](https://arxiv.org/html/2609.05028)
- **Discrete-diffusion policy gradients** (arXiv 2502.01384) are relevant for token-mask generators, title only. — [arXiv 2502.01384](https://arxiv.org/pdf/2502.01384)

### Inferences
- **Opportunity for the engine.** Any mask-emitting diffusion or flow model can take Flow-GRPO-style fine-tuning with rewards computed by one engine pass on the final sample. Examples are synthetic-data generators for segmentation, ILT mask generators, and floorplan or label-map diffusion.
  - Candidate rewards: β0/β1 match to a target topology, hole-free or hole-count constraints, vertex or jog budgets, a rectangle-partition shot proxy, and minimum feature or area per component.
  - Group-relative normalization only needs within-group ranking (LithoGRPO), so exact-but-cheap integer statistics are ideal reward signals.
  - Because the engine is non-differentiable, it suits DDPO, Flow-GRPO and DanceGRPO, not DRaFT or ReFL. Soft bit-quad surrogates (Fast-χ, Minkowski loss, from round 2) would be the DRaFT-compatible counterparts.
- At diffusion RL scale, reward evaluations run in the millions; for example, Flow-GRPO-style runs use groups of 16–24 over thousands of prompts and steps. A CPU polygonizer or ILP per sample would become noticeable for large masks. That is the second domain, after lithography, where engine speed could matter, but **no paper yet measures it** (inference).

### Gaps
- I found no DDPO, DPOK, ReFL, DRaFT, Flow-GRPO or DanceGRPO paper whose output is a segmentation mask or label map and whose reward is topological or polygon-based.
- I found no reward-guided floorplan or label-map diffusion paper with validity rewards.
- The specific sub-rewards in the medical CRM paper were not verified.
