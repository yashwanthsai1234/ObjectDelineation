# P_GEO: 3D / CAD / Physical AI / GPU systems annotated bibliography

## Skills used
- huggingface-skills:huggingface-papers: gave the paper-markdown URL pattern (huggingface.co/papers/<id>.md) and the arXiv HTML fallback. Five papers were not indexed on Hugging Face; arXiv HTML was used for them.
- nvidia-skills:nvidia-skill-finder: loaded. It routes NVIDIA topics to the live skills catalog. The catalog was not queried and nothing was installed. The SimReady Foundation validator pages were found by web search instead.
- engineering:code-review, engineering:debug, engineering:system-design, engineering:architecture, engineering:documentation: loaded as required. They gave only review, debug and design checklists. They were used to structure the "Use for us" lines and to check each entry for unsupported numbers.
- huggingface-skills:huggingface-best: loaded. It is a model-leaderboard finder; it did not apply to paper reading.

## Count line
Papers read in full: 20 / abstract + tables only: 3 (2503.06800, 2605.01536, 2609.28731) / new 2026 additions: 6 (read key sections: abstract, method summary, results, conclusion).

## Method notes
- Four sub-readers (Sonnet) read the 23 listed papers from local text; the librarian read the 6 new papers and checked their numbers.
- Numbers stripped by HTML or markdown conversion are marked "not stated". No number was guessed.
- "Relevance" is added for sorting.

## Top 5 for us
1. 2609.34494 ConCAD, 2. 2505.19713 CAD-Coder, 3. 2608.14498 Rollplex, 4. 2609.35641 Verifiable Visual Rewards, 5. 2607.07459 EmbodiedGen V2. Close next: 2605.10865 BenchCAD, 2610.10959 TopoGPU, 2607.08891 Ortho2CAD, 2610.37297 TRUCO.

## A. CAD (image/text to CAD code)

### 2501.19054 — Text-to-CAD Generation Through Infusing Visual Feedback in Large Language Models (January 2025)
- Read: full text
- Core idea: CADFusion trains an LLM to write CAD parametric sequences. It alternates supervised steps with DPO steps that use rendered images.
- Key result: LLaMA-3-8B with LoRA on DeepCAD text pairs (4 A6000 GPUs): LVM score 8.96, invalidity ratio 6.20%, versus GPT-4o 5.13 and 37.93%.
- Use for us: It shows that a non-differentiable renderer can drive training, and that alternating reward steps with supervised steps stops format collapse; our exact checker should give a stronger reward than their VLM judge.
- Relevance: medium

### 2505.06507 — Text-to-CadQuery: A New Paradigm for CAD Generation with Scalable Large Model Capabilities (May 2025)
- Read: full text
- Core idea: The model writes executable CadQuery code, not custom command sequences. The authors fine-tune six LLMs (124M to 7B) on 170,000 text-CadQuery pairs.
- Key result: Qwen2.5-3B reaches 69.3% Gemini-judged match (10.5% above the Text2CAD baseline) and volumetric IoU median 0.9868. The abstract states a 48.6% Chamfer Distance drop. Training is plain cross-entropy, no RL.
- Use for us: It gives a ready CadQuery dataset recipe for the SFT start of an image-to-CAD run, and shows that a small model suffices; its voxel IoU (0.02 resolution) is the metric our GPU verifier should replace exactly.
- Relevance: medium

### 2505.19713 — CAD-Coder: Text-to-CAD Generation with Chain-of-Thought and Geometric Reward (May 2025)
- Read: full text (the extracted text lost the numbers inside Tables 1 to 3; I used numbers quoted in the prose)
- Core idea: SFT on CadQuery, then GRPO with a reward that mixes Chamfer Distance and a format check. The reward scales linearly between CD 1e-5 (reward 1.0) and CD 0.5 (reward 0.01).
- Key result: Qwen2.5-7B-Instruct, 8 A800 GPUs: mean CD 6.54 (x10^3) and median CD 0.17, versus 74.55 and 0.33 for SFT only. GRPO took 146 hours; SFT took 7 hours. IR values are not stated in the text.
- Use for us: It is the closest recipe to our setup (code output, executable geometry reward, GRPO, k=8). Its CD-only reward failed after about 200 steps when code became invalid, so we must keep a validity term; its 146 GPU-hours show reward execution cost matters.
- Relevance: high

### 2510.17157 — GACO-CAD: Geometry-Augmented and Conciseness-Optimized CAD Model Generation from Single Image (October 2025)
- Read: full text (Tables 1 to 4 lost their numbers in the extracted text)
- Core idea: SFT adds depth and normal maps as extra image inputs. RL (GSPO) adds a group length reward that favours short code among accurate outputs, with weights that change linearly during training.
- Key result: Qwen2-VL-2B, one H800 for SFT and four H800 for RL: generated code is about 6.1% shorter on DeepCAD and 5.7% shorter on Fusion360. IoU, CD and IR values are not stated in the extracted text.
- Use for us: It shows a 2B VLM can do image-to-CAD with an IoU reward (I_high 0.8, I_low 0.4), and the scheduled reward weights and length reward are cheap tricks to copy into our reward.
- Relevance: medium

### 2605.10865 — BenchCAD: A Comprehensive, Industry-Standard Benchmark for Programmatic CAD (May 2026)
- Read: full text (main body and training appendix; Tables 2, 3, 6 and 10 lost most numbers in the extracted text)
- Core idea: A benchmark of 17,900 execution-verified CadQuery parts in 106 industrial families. It has four tasks: image-to-code, code editing, image QA and code QA.
- Key result: Best frontier model (Gemini 3.1 Pro thinking) scores 0.318 total on image-to-code. Code QA reaches about 0.838 versus 0.587 for image QA. A Qwen3-VL-2B trained on it scores 0.7631 on in-distribution code generation, with a clear gap on held-out families.
- Use for us: It is a ready test set and training set for a Qwen3-VL-2B CAD run, with a reward of 0.2 x essential-op recall + 0.8 x IoU, 16 rollouts and GRPO. It also supplies a rotation-invariant voxel IoU (24 rotations) for the GPU verifier.
- Relevance: high

### 2609.34494 — ConCAD: Constraint-Aware Image-to-CAD Generation with Dual-Granularity Rewards (September 2026)
- Read: full text (the ablation tables were only partly read)
- Core idea: GRPO with two rewards. One is a regex check of code constraints before execution. The other is volumetric IoU after inertia-based alignment. A B-rep constraint metric (G-CSR) checks design intent.
- Key result: Qwen3-VL-2B, one A800 for training, DeepCAD test: IoU 0.734, CD 0.009, IR 0%, G-CSR 0.913. The SFT start gave 0.717, 0.011, 1% and 0.903. On Zero2CAD, IoU rose from 0.443 to 0.507.
- Use for us: It gives a 2B-model, single-GPU RL recipe that fits our case. It warns that a pure execution reward collapses to trivial cubes. Its IoU uses 8,192 Monte Carlo points, so our exact GPU IoU is a direct upgrade.
- Relevance: high

## B. RL systems, 3D medical

### 2604.26779 — Accelerating RL Post-Training Rollouts via System-Integrated Speculative Decoding (April 2026)
- Read: full text
- Core idea: The authors put speculative decoding (EAGLE-3 draft) into NeMo-RL with a vLLM backend. It speeds up rollouts and keeps the policy's output distribution.
- Key result: 8B reasoning RL, synchronous, 8 GB200 nodes (4 GPUs each). Generation time falls from 100.0 s to 56.6 s (1.8x) on RL-Zero, and from 133.6 s to 87.0 s (1.5x) on RL-Think. Step speedup is 1.41x and 1.35x. Generation is 65-72% of step time. Draft length 3 is best; longer drafts are slower. N-gram drafting is slower than plain decoding. Async at policy lag 1: step time 75.0 s to 60.5 s (1.24x). The 2.5x at 235B is a simulator projection, not a measurement.
- Use for us: Our rollouts are short polygon outputs from a 2B VLM, so decode share and acceptance length are likely small; use its step-time formula to check the gain before any work.
- Relevance: low

### 2604.07853 — QaRL: Rollout-Aligned Quantization-Aware RL for Fast and Stable Training under Training–Inference Mismatch (April 2026)
- Read: full text
- Core idea: Run the rollout in low-bit weights, and run the learner forward in the same low-bit GEMM. A sequence-level dual-clipping loss (TBPO) drops responses with garbled "error tokens".
- Key result: Qwen3-30B-A3B-Base, math, OpenR1-Math-46K: average math score 51.2 for QaRL, 52.1 for BF16, 45.7 for quantized rollout alone. Training is 1.3x faster than BF16 (1.4x for quantized rollout alone). Qwen3-8B-Base: 48.9 vs 51.8 (BF16) vs 43.9 (quantized rollout). GPU type is not stated in the sections read.
- Use for us: It names the sampler-learner mismatch that can hurt our RL; the cheap part to copy is the mismatch weight and sequence-level clipping, not the W4 rollout. Their approach tolerates mismatch; ours needs bit-exact scoring only for the reward, which is separate.
- Relevance: low

### 2606.18967 — EfficientRollout: System-Aware Self-Speculative Decoding for RL Rollouts (June 2026)
- Read: full text
- Core idea: Make a 4-bit weight draft model from the current policy at each step. Switch speculation on only when the active batch is small, using a roofline model.
- Key result: 8x A100-80GB, veRL + vLLM, GRPO, batch 128, group 8, max 8k tokens, 100 steps. Qwen2.5-7B: rollout 82.4 s to 66.3 s (-19.6%), step 132.6 s to 115.7 s (-12.7%). Qwen2.5-14B: step -10.8%. Llama3.1-8B: step -7.9%. A learned EAGLE3 drafter made Llama3.1-8B step 25.6% slower. Reward curves match the baseline. Rollout is about 70% of step time.
- Use for us: It gives a concrete rule to decide on one RTX 5090 whether decode is long-tailed enough for speculation; with short outputs it likely is not, and the roofline toggle is a cheap test.
- Relevance: medium

### 2608.14498 — Rollplex: Cross-Phase GPU Spatial Sharing for Vision Language Model Post-Training (August 2026)
- Read: full text
- Core idea: In VLM RL, the prompt and image prefix work of reference scoring and training does not depend on the response. Run it on the same GPUs during rollout decode, using CUDA MPS.
- Key result: Qwen2.5-VL-32B, Video-R1 setup, 32 H800 GPUs: 1.23x-1.30x faster than serial colocation and 1.57x-2.24x faster than disaggregation. Decode uses under 19% of active SMs. Default MPS is within 2.7% of the best policy; 70/30 limits cost 1.30x-3.44x. Rewards track the baseline. Prompt tokens are 79%-98% of video samples (11%-19% for text RL). Gain is capped by the shorter of decode time and prefix time.
- Use for us: Our 2B VLM sees image-heavy prompts and short outputs, so prefix work is large; overlapping image encoding and prefix passes with generation on one 5090 under default MPS is the closest match to our GPU-use goal.
- Relevance: high

### 2511.11450 — VoxTell: Free-Text Promptable Universal 3D Medical Image Segmentation (November 2025)
- Read: full text
- Core idea: A 3D UNet-style model maps a free-text prompt to a volumetric mask. Text is fused into the decoder at every scale (multi-stage fusion) with deep supervision.
- Key result: Trained on 62K+ CT, MRI and PET volumes (158 public datasets, 1K+ classes), 64 A100 GPUs, about six days. Zero-shot mean Dice 70.85 vs 51.23 for SAT (next best) over 11 unseen datasets. Ablation: late fusion 51.68-55.11, multi-stage with deep supervision 62.55, plus batch 128 gives 69.43. Real report prompts: 50.2 Dice. ReXGroundingCT after fine-tuning: 28.2 Dice vs 13.1 for SAT. Limit: fails on knee MRI, not seen in training. No Betti or topology measure is reported.
- Use for us: Gives a Dice baseline and a prompt-robustness test idea for a future 3D phase; it is a segmentation network, not a polygon writer or RL reward, so it is a comparison point, not a method to copy.
- Relevance: low

## C. GPU geometry, sim-ready assets, 3D benchmarks

### 2605.01536 — The Antipodal Method: Fast, Accurate, and Robust 3D Generalized Winding Numbers (May 2026)
- Read: abstract + tables (intro, related work and results text read; method equations skimmed; most speedup numbers were stripped from the text)
- Core idea: The winding number equals a signed ray-surface intersection count plus a line integral over the surface boundary projected onto the unit sphere. This avoids surface integrals.
- Key result: Speedup factors and queries-per-second are blank in the text (stripped), so not stated here. It reports "up to three orders of magnitude" faster than the exact hierarchical method on CPU. Setup: AMD Ryzen 9 7950X3D CPU, NVIDIA RTX 4090 GPU (OpenGL). It claims 120 FPS for 4K slices on moderate meshes.
- Use for us: It gives an insidedness test for open or non-manifold meshes, so a GPU voxel verifier could fill voxels from CAD meshes that are not watertight.
- Relevance: medium

### 2609.28731 — CuACD: A Fully GPU-Resident Approximate Convex Decomposition (September 2026)
- Read: abstract + tables (intro, design and experiment text read; speedup numbers were stripped from the text)
- Core idea: Each candidate cut runs as one warp. A device-side heap allocator lets the whole decomposition run without the host sizing launches.
- Key result: "More than an order of magnitude" faster than prior search-based ACD, at matched or lower mean concavity. Exact factors are blank in the text. Setup: single RTX 4090, V-HACD, PartNet-Mobility and Objaverse subset. Zero failures on all inputs.
- Use for us: It shows how to keep irregular geometry kernels resident on one GPU (warp per task, device heap, no host loop), which we can copy for our CUDA checker.
- Relevance: medium

### 2607.07459 — EmbodiedGen V2: An Agentic, Simulation-Ready 3D World Engine for Embodied AI (July 2026)
- Read: full text
- Core idea: A pipeline turns text or image into sim-ready assets and worlds. It uses a generate-verify-retry loop with quality gates, mesh repair and convex decomposition.
- Key result: Single RTX 4090, 200 held-out assets: 96.5% human acceptance and 98.6% collision success at 2.6 min per asset. Without mesh fixing, runtime rises from 2.6 to 21.3 min. 83.3% of 150 task worlds usable. Cited RL study: simulation success 9.7% to 79.8%, real robot 21.7% to 75.0%.
- Use for us: It is a template for our Physical AI asset certification, with explicit gates and a single-attempt pass rate per gate (76.2%, 75.9%, 91.0%) that an exact checker could replace.
- Relevance: high

### 2608.23869 — Gen2Physics: Grounding Generated 3D Meshes in Physics via Multi-View Material Decomposition (August 2026)
- Read: full text (Table 1 and Table 3 values are not in the text; only the prose numbers were used)
- Core idea: It lifts 2D material masks from 11 views onto mesh faces by majority vote, then a VLM repairs the labels. HoloPart completes each part into a watertight sub-mesh.
- Key result: Material mIoU 48.27% on PartNet-Material (100 objects, 1100 masks), against 15.56% for PUGS and 14.54% for NeRF2Physics. ViT only: 22.58%. Mass estimation on ABO-500 is within 0.2% ADE of the best method. Hardware not stated in the paper.
- Use for us: Its per-face majority vote across views and VLM-critique step is a pattern for our image-to-outline work, and its watertight sub-meshes are a target for an exact watertightness check.
- Relevance: medium

### 2610.40333 — RFW-A: Agentic Procedural Generation of 400K Articulated Assets over 500 Categories (October 2026)
- Read: full text (appendix tables skimmed)
- Core idea: Agents write a reusable category-level generator program once, and procedural execution then makes many assets with no more model calls. Each generator is released only after strict validation.
- Key result: 402,410 assets, 500+ categories. Mean generator cost $10.62; $0.053 per asset over 200 assets, against $1.237 per asset for per-asset agentic synthesis (23-fold). Strict URDF validity 1.0000. Stable 1-second rollout 0.995 in four simulators (Core-200 subset). Real-robot success 6.2/10 to 8.0/10 with edited assets.
- Use for us: Its validity checks (watertight-like collision coverage, kinematic sweeps, cross-simulator rollout) and amortized-generator idea suit CAD code and asset certification with our exact checker as the gate.
- Relevance: high

### 2604.02580 — VoxelCodeBench: Benchmarking 3D World Modeling Through Code Generation (April 2026)
- Read: full text
- Core idea: Models write Python for a voxel API in Unreal Engine. Humans rate the renders on position, material, shape and visual quality across 220 tasks.
- Key result: GPT-5 reaches 87.9% shape correctness, Claude Sonnet 4.5 and GPT-5 Mini 80.4%, Gemini 3 Pro 19.5%. Geometric construction is hardest (GPT-5 66.7%). Hardware not stated in the paper.
- Use for us: It shows that human rating of renders is the weak point of code-to-3D tasks, and an exact voxel verifier (b0/b1/b2, solid IoU) could replace it as a reward.
- Relevance: medium

## D. Physics and world-model benchmarks, verifiable rewards

### 2503.06800 — VideoPhy-2: A Challenging Action-Centric Physical Commonsense Evaluation in Video Generation (March 2025)
- Read: abstract + tables
- Core idea: The paper gives 200 real-world actions with prompts for video models. Humans rate meaning match and physical sense. A 7B auto-rater scales the rating.
- Key result: Best model Wan2.1-14B scores 32.6% joint score on the full set and 21.9% on the hard subset (human rating, 590 test prompts per model). Conservation of momentum and mass are violated about 40% of the time. The auto-rater beats Gemini-2.0-Flash-Exp by 81% (meaning match) and 236% (physical sense) in relative correlation. GPU and hardware cost: not stated in the paper.
- Use for us: It shows that learned video judges are weak and human labels are costly, so an exact checker as reward is a strong choice for our CAD and sim-asset work.
- Relevance: low

### 2504.00983 — WorldScore: A Unified Evaluation Benchmark for World Generation (April 2025)
- Read: full text
- Core idea: The benchmark splits world generation into next-scene tasks with camera paths. It scores 10 metrics in three groups: controllability, quality, dynamics. The metrics use tools such as DROID-SLAM (3D consistency) and optical flow.
- Key result: 3,000 test examples, 20 models. WorldScore-Static: WonderWorld 72.69, LucidDreamer 70.40, best video model CogVideoX-I2V 62.15. Best video camera control is 40.22 (CogVideoX-T2V). Hardware cost: not stated in the paper.
- Use for us: It gives a template for a geometric 3D consistency check (reprojection error from SLAM depth), which is a soft cousin of our exact voxel and Betti verifier for sim-ready assets.
- Relevance: low

### 2505.00337 — T2VPhysBench: A First‑Principles Benchmark for Physical Consistency in Text‑to‑Video Generation (May 2025)
- Read: full text
- Core idea: 12 physical laws, 7 prompts per law, 84 prompts per model. Three students rate each video on four levels. The paper adds hint and counterfactual tests.
- Key result: Best average score is 0.42 (Wan 2.1), worst 0.19 (SD Video), across 10 models. Conservation laws score lowest (best 0.29). Hints rarely help: Newton's third law rises from 0.25 to 0.88, but Newton's first law drops from 0.12 to 0.00. The numeric value of each rating level was stripped from the text.
- Use for us: It shows that prompt hints cannot fix physics errors and that humans are the only judge here, which supports training with an exact verifier reward instead.
- Relevance: low

### 2507.15824 — PhysVidBench: Language-Grounded Evaluation of Physical Commonsense in Text-to-Video Models (July 2025)
- Read: full text
- Core idea: 383 prompts from PIQA, each with a detail-rich variant. A VLM captions the video. An LLM answers yes/no questions from the caption only.
- Key result: Best average accuracy across 12 models is 36.2%. Pearson correlation with humans on 120 videos is 0.45 to 0.69 for this method and 0.21 to 0.47 for the best direct VLM judge. Two video VLMs answer "yes" to nearly all questions. Open models ran on a single NVIDIA A40. The refinement gains (+11.1 and +11.9 points) are printed, but the starting values were stripped.
- Use for us: Direct VLM judges collapse or score poorly on physics, which backs our choice of an exact checker over a VLM judge for the reward.
- Relevance: low

### 2601.21282 — WorldBench: Disambiguating Physics for Diagnostic Evaluation of World Models (January 2026)
- Read: full text
- Core idea: Each test isolates one physics concept. Scores come from SAM2 masks against ground-truth masks (foreground mIoU), and from physical constants fitted from tracked object positions (gravity, friction, viscosity).
- Key result: Foreground mIoU averages 0.4508 for Cosmos-1 (33 frames, 7B), 0.2573 for Cosmos-1 (121 frames), and 0.3488 for Cosmos-2.5 (2B). Models give realistic paths but wrong constants. Free-fall gravity estimates run from 10.831 (Cosmos-1 AR) to 0.939 (Hunyuan) against 9.81. Dataset: 469 intuitive-physics videos, 279 parameter videos. Hardware: not stated in the paper.
- Use for us: It scores generated video with mask IoU against ground truth and fits exact physical constants; the same per-frame IoU could run as a GPU kernel inside our checker for Physical AI tests.
- Relevance: medium

### 2609.35641 — Verifiable Visual Rewards Transfer from Synthetic Scenes to Natural Prompts (September 2026)
- Read: full text
- Core idea: A task generator makes scenes of coloured shapes, and from each scene a prompt plus a deterministic pixel verifier (no VLM, no learned detector). The pass and partial-credit scores are the reward for RL (Flow-GRPO).
- Key result: SD3.5-M accuracy on VVRBench rises from 2.81% to 28.27% (trained on VVR-Easy) and to 46.60% (VVR-Matched). GenEval is 0.729 vs 0.688 for the GenEval2-trained model. Human win rate is 71.6% on 160 natural prompts. GPT-Image-2.5-Sunburst solves 21.39% of VVRBench-Challenge. Training hardware: not stated in the text I read.
- Use for us: It is the closest match to our design: an exact program verifier with partial credit and a penalty factor against reward hacking, plus a generator with set difficulty for curricula in CAD and voxel tasks.
- Relevance: high

## E. New 2026 additions (not in the given list)

### 2610.10959 — GPU-Accelerated Computation of Persistent Homology for Topological Analysis of Image Data (October 2026) [NEW]
- Read: full text (abstract, intro, contributions, experiments, training study; proofs skimmed)
- Core idea: TopoGPU streams a 2D or 3D image in chunks through GPU kernels. The CPU then reduces the boundary matrix.
- Key result: average end-to-end speedup 53.24x, maximum 198.01x, over Cubical Ripser on every benchmark. In a PARSE pulmonary artery network (128^3 volumes, 480 epochs), epoch time falls from 309.76 s to 48.56 s.
- Use for us: It gives a ready GPU path for b0, b1, b2 of voxel masks, but it is a persistence tool, not an exact integer Betti counter.
- Relevance: high

### 2610.37297 — TRUCO: A Scalable Lock-Free Algorithm for Parallel Contour Extraction (October 2026) [NEW]
- Read: full text (abstract, related work, results, conclusion; proofs skimmed)
- Core idea: Each thread traces whole contours in a horizontal stripe. Two rules (row-local writes, top-row ownership) make the output equal to Suzuki-Abe, with the same hierarchy.
- Key result: 2.3 to 17.8 times faster than the fastest sequential method with the same output, on multi-core CPUs. Output matched Suzuki-Abe in 5.2 million executions.
- Use for us: It shows how to prove a parallel contour tracer bit-identical to a reference; the ownership rule can guide a CUDA hole and outline checker. It targets CPUs, not GPUs.
- Relevance: high

### 2607.08891 — Ortho2CAD: 3D CAD generation from orthographic drawings using vision language models (July 2026) [NEW]
- Read: full text (abstract, method, RL reward, main table; appendices skimmed)
- Core idea: A small VLM writes CadQuery code from a drawing. GRPO-style RL uses the solid IoU of the executed code as reward, with reward 0 for invalid code.
- Key result: RL model reaches 100 % valid code on four test sets (100 examples each). On F360Recon the IoU is 0.5601, against 0.3697 for the SFT model.
- Use for us: It is the closest recipe to our CAD plan: execute code, compute exact solid IoU, give 0 on failure. Its IoU needs a CAD kernel on CPU.
- Relevance: high

### 2606.13368 — IterCAD: An Iterative Multimodal Agent for Visually-Grounded CAD Generation and Editing (June 2026) [NEW]
- Read: full text (abstract, RL section, main tables, ablation, conclusion)
- Core idea: A 4B agent writes CadQuery in a sandbox loop. RL uses geometry rewards and Geometry-Viable Prefix Masking, which masks tokens after an error cascade or a stalled Chamfer Distance.
- Key result: IterCAD-Draw benchmark, invalid ratio 64.60 % (Qwen3.5-4B base) to 0.30 % (full), mean Chamfer 10.48 to 5.09. Average 2.97 turns.
- Use for us: The prefix mask fixes credit assignment when an early valid geometry is followed by errors; reuse it for multi-turn CAD repair.
- Relevance: medium

### 2610.03453 — I2CD: Direct Image-to-Convex Decomposition for Simulation-Ready Collision Geometry (October 2026) [NEW]
- Read: full text (abstract, intro, simulator study; method internals skimmed)
- Core idea: A frozen Hunyuan3D-2 model plus a 38 M-parameter head emits K convex polytopes from one image. The output loads in physics engines without repair.
- Key result: On 227 held-out objects it has the highest volumetric IoU of eight reconstruct-then-decompose pipelines, 6 to 37 times faster. MuJoCo, PyBullet and PhysX silently replace a raw mesh by its convex hull, the wrong collision shape for about 60 % of objects; I2CD geometry is preserved on at least 98.4 % of objects in four engines.
- Use for us: It shows that "it loads" hides silent collision-shape swaps; a sim-ready certificate must test the shape the engine really uses.
- Relevance: high

### 2607.09207 — Bidirectional Resource Scheduling for Disaggregated and Asynchronous RL Post-Training (July 2026) [NEW]
- Read: full text (abstract, design summary, evaluation; planner math skimmed)
- Core idea: BiDiRL keeps separate rollout and training GPU pools. An idle pool temporarily borrows the other pool's role.
- Key result: Up to 1.94x throughput over veRL, AReaL and ROLL on two 32-GPU testbeds (A6000, H100), with unchanged convergence.
- Use for us: It confirms that idle time is the main loss, but it needs many GPUs; on one RTX 5090 only the idea of filling bubbles transfers.
- Relevance: medium

Extra, not an arXiv paper: Alonso-Jorda et al., Accelerated border tracking in binary images with GPUs, J. Supercomputing, April 2026 (CUDA triad linking for contours; seen as abstract only, not counted).
