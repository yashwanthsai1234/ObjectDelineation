# P_RL: RL algorithms, RL infrastructure, RL for vision grounding (annotated bibliography)

## Skills used
- huggingface-skills:huggingface-papers: gave the fetch method. I pulled `huggingface.co/papers/<id>.md` (arXiv HTML converted) for 32 papers. For 5 papers (2512.24023, 2509.21792, 2606.13156, 2512.23075, 2608.13549) the HF page returned 404 HTML, so I used `arxiv.org/html/<id>` as the skill's fallback.
- nvidia-skills:nvidia-skill-finder: loaded. No NVIDIA catalog lookup needed (this task reads papers, no NVIDIA product use). No contribution.
- engineering:code-review, engineering:debug, engineering:system-design, engineering:architecture, engineering:documentation: loaded as required. They shaped the structure: each entry names the problem it fixes (debug lens) and the system trade-off (design lens). No code was reviewed.
- huggingface-skills:huggingface-best: loaded as required. Not used (no model ranking asked).

## Counts
- Papers in the list: 37. Read "full text": 23. Read "abstract + tables": 14.
- New 2026 additions: 5 (all read: 4 with targeted full-text passages, 1 abstract only).
- Method note: I fetched full text for every paper. I did not read every page. For "full text" papers I read the abstract, the contribution list, the method passages and the results or ablation passages (by targeted extraction). For "abstract + tables" papers I read the abstract plus the key result sentences or tables I could find. Numbers below come only from those passages. If a number is not stated, I say so.
- Extraction limit: math inside HTML was stripped, so a few numbers (marked "not legible") did not survive.

## Problems in our run, and the best papers for each
- Trainer-vLLM mismatch (IS ratio min 0 in 98 % of steps): 2510.26788, 2609.20807, 2605.14220, 2601.18150, 2510.18855, 2506.13585, 2510.13786, 2602.19362, 2507.18071.
- Entropy collapse: 2505.22617, 2503.14476, 2506.01939, 2505.24864, 2508.10751.
- No signal for all-miss groups: 2605.15951, 2601.07695, 2503.20783, 2601.05242, 2510.13786, 2501.03262.
- Centred-train vs shifted-test: no paper in this set targets it. 2606.13156 warns about refinement loops. 2504.07954 and 2601.19228 show reward and data choices.
- 40 % idle GPU: 2509.19128, 2505.24298, 2506.05433 (already in use), 2509.23232, 2509.21792, 2501.12599, 2510.18855.

---

## Part A. RL algorithms

### 2503.14476 — DAPO: An Open-Source LLM Reinforcement Learning System at Scale (March 2025)
- Read: full text
- Core idea: Four fixes make GRPO work at scale: Clip-Higher, dynamic sampling, token-level loss, soft overlong punishment. The authors name entropy collapse, reward noise and instability as the failures of plain GRPO.
- Key result: Qwen2.5-32B base, AIME 2024: plain GRPO reaches 30. DAPO reaches 50, with 50 % of the training steps of DeepSeek-R1-Zero-Qwen-32B.
- Use for us: Dynamic sampling drops groups with zero gradient. It saves compute on all-miss groups but gives no signal. Clip-Higher is the cheap first test against entropy collapse.

### 2503.20783 — Understanding R1-Zero-Like Training: A Critical Perspective (March 2025)
- Read: abstract + tables
- Core idea: GRPO has two biases: length normalization and std normalization. They push up the length of wrong answers. Dr.GRPO removes both terms.
- Key result: 7B base, AIME 2024: 43.3 % with the minimalist recipe. Dr.GRPO improves token efficiency with the same reasoning accuracy.
- Use for us: Std normalization blows up tiny reward differences. In an all-miss group with near-zero IoU, Dr.GRPO keeps the advantage small instead of large. Test it as a one-line change.

### 2505.22617 — The Entropy Mechanism of Reinforcement Learning for Reasoning Language Models (May 2025)
- Read: full text
- Core idea: Entropy and score follow R = -a exp(H) + b. Entropy falls because high-probability, high-advantage tokens have positive covariance. Clip-Cov and KL-Cov limit updates on high-covariance tokens.
- Key result: Qwen2.5, DAPO-MATH training, average score (Table 2): 7B GRPO 38.6, Clip-Cov 40.4, KL-Cov 40.6. 32B GRPO 45.8, Clip-Cov 50.3, KL-Cov 52.2. Entropy stays higher.
- Use for us: Direct fix for our entropy collapse. We can log the covariance between log-probability and advantage and apply KL-Cov to the top 0.2 % tokens (paper setting for 32B: k = 2e-4).

### 2506.01939 — High-Entropy Minority Tokens Drive Effective Reinforcement Learning for LLM Reasoning (June 2025)
- Read: abstract + tables
- Core idea: Few tokens have high entropy and act as forks. RLVR changes mostly those tokens. Training only on the top 20 % entropy tokens works.
- Key result: Only the 20 % tokens: Qwen3-32B +11.04 on AIME'25 and +7.71 on AIME'24 over full-gradient updates. Qwen3-14B +4.79 and +5.21. Qwen3-8B stays comparable.
- Use for us: In polygon text, the fork tokens are probably the vertex coordinates. We can restrict the gradient to high-entropy tokens and see if entropy stays alive.

### 2505.24864 — ProRL: Prolonged Reinforcement Learning Expands Reasoning Boundaries in Large Language Models (May 2025)
- Read: full text
- Core idea: Long RL with KL control, periodic reference-policy reset and diverse tasks keeps finding new solutions.
- Key result: Nemotron-Research-Reasoning-Qwen-1.5B vs DeepSeek-R1-1.5B: average pass@1 +14.7 % (math), +13.9 % (coding). Gains continue after 2,000 steps.
- Use for us: A small 1.5B model can keep improving. The KL-reset trick gives a plan if our policy drifts from the reference after many steps.

### 2508.10751 — Pass@k Training for Adaptively Balancing Exploration and Exploitation of Large Reasoning Models (August 2025)
- Read: abstract + tables
- Core idea: Use Pass@k as the reward, with an analytic advantage. Exploration and exploitation help each other.
- Key result: Gain is shown on Enigmata with Qwen2.5-7B-Instruct. I did not extract a single headline number.
- Use for us: Our reward is a best-of-group style check on a set of doctors. A pass@k-shaped advantage can reward a group where one answer is a hit. Low priority until entropy is fixed.

### 2507.18071 — Group Sequence Policy Optimization (July 2025)
- Read: full text
- Core idea: GSPO uses a sequence-level importance ratio and clips whole sequences. It fixes the token-level noise of GRPO.
- Key result: More stable than GRPO and stabilizes MoE training. No single number extracted. The paper says sequence likelihoods tolerate train-inference precision differences better, so it can skip recomputing old-policy likelihoods.
- Use for us: Direct help for the mismatch problem. SimpleSeg (2601.19228) trains a polygon-output model with GSPO.

### 2511.20347 — Soft Adaptive Policy Optimization (November 2025)
- Read: full text
- Core idea: SAPO replaces hard clipping with a smooth, temperature-controlled gate. It down-weights only the off-policy tokens and keeps the near-on-policy ones.
- Key result: Better Pass@1 than GSPO and GRPO under the same budget on math. It also improves Qwen3-VL-30B-A3B training (AIME25, LiveCodeBench v6, ZebraLogic, MathVision). No numbers extracted.
- Use for us: A hard clip kills our whole sequence when a few tokens mismatch. SAPO keeps the good tokens. It was tested on a VLM family.

### 2507.20673 — Geometric-Mean Policy Optimization (July 2025)
- Read: abstract + tables
- Core idea: Use the geometric mean of token-level rewards instead of the arithmetic mean. It is less sensitive to extreme importance ratios.
- Key result: GMPO-7B improves average Pass@1 over GRPO by up to 4.1 % on math.
- Use for us: One-line way to damp the extreme ratios we see. Compare with GSPO and SAPO.

### 2602.05261 — Revealing and Controlling Response Length Variation in RLVR (February 2026)
- Read: full text
- Core idea: GSPO has a length bias. LUSPO removes it and stops length collapse.
- Key result: Qwen2.5-VL-7B-Instruct: LUSPO beats GSPO by 5.1 % (WeMath) and 6.0 % (LogicVista).
- Use for us: Polygon length matters for us. If we adopt GSPO, use the length-unbiased form so the model does not shrink to short polygons.

### 2501.03262 — Stabilizing Critic-Free Policy Optimization with Global Normalization (January 2025; REINFORCE++)
- Read: full text
- Core idea: Normalize advantages over the whole batch, not per prompt. The paper shows per-group normalization is biased.
- Key result: Chat-Arena-Hard: single-sample variant 46.7 vs GRPO 46.8, with shorter answers (832 vs 860 tokens). On a small dataset GRPO with local normalization overfits completely.
- Use for us: Batch-level normalization keeps a small but real signal from all-miss groups. ScaleRL (2510.13786) also uses batch-level advantage normalization.

### 2601.05242 — GDPO: Group reward-Decoupled Normalization Policy Optimization for Multi-reward RL Optimization (January 2026)
- Read: full text
- Core idea: Normalize each reward term on its own, then sum. Plain GRPO collapses different reward mixes into the same advantage.
- Key result: Qwen2.5-Instruct-1.5B tool calling: about +2.7 % average accuracy and over 4 % in format ratio over GRPO. Up to 80 % less length-exceeding on AIME.
- Use for us: Our reward has IoU, boundary F and Betti topology. Normalize each term per group before the sum so a rare topology term is not drowned.

### 2503.06639 — GRPO's Effective Loss, Dynamics, and Success Amplification (March 2025)
- Read: abstract + tables
- Core idea: Theory. With binary reward, GRPO is a contrastive loss. The success probability follows a recurrence and converges to a fixed point above the reference.
- Key result: Proof only. No benchmark number.
- Use for us: Explains why GRPO needs some successes in the group. All-miss groups do not move. Background for the reward-signal fix.

### 2504.05118 — VAPO: Efficient and Reliable Reinforcement Learning for Advanced Reasoning Tasks (April 2025)
- Read: abstract + tables
- Core idea: A value-model PPO with length-adaptive GAE, value pre-training and several DAPO tricks.
- Key result: Qwen2.5-32B, AIME 2024: 60.4. More than 10 points above DAPO and R1-Zero-Qwen-32B. 5,000 steps, no crashes across runs.
- Use for us: A value model costs memory on one RTX 5090. Does not apply now. Keep as a fallback if GRPO stays flat.

### 2501.12599 — Scaling Reinforcement Learning with LLMs (January 2025; Kimi k1.5)
- Read: full text
- Core idea: Simple RL (online mirror descent, no value model, no process reward). Partial rollouts, curriculum sampling and prioritized sampling.
- Key result: AIME 77.5, MATH-500 96.2, MathVista 74.9. Prioritized sampling: pick problem i with weight 1 - s_i (s_i = its success rate).
- Use for us: Prioritized sampling sends compute to hard crops. Partial rollouts address the long-tail idle time.

### 2512.23075 — Trust Region Masking for Long-Horizon LLM Reinforcement Learning (December 2025)
- Read: abstract + tables
- Core idea: Trust-region error bounds grow with sequence length. Masking whole sequences that break the trust region gives non-vacuous guarantees. The bounds depend on max token-level divergence.
- Key result: Theory with a stability demonstration. No number extracted.
- Use for us: Gives a principled rule for the mismatch: mask the whole sequence if any token diverges too much, and log the acceptance rate as a diagnostic.

### 2602.10693 — Variational Sequence-Level Soft Policy Optimization for Stable Off-Policy LLM Training (February 2026)
- Read: abstract + tables
- Core idea: VESPO reshapes sequence-level importance weights with a closed-form kernel, with no length normalization.
- Key result: Stable up to 64x staleness and fully asynchronous runs, on dense and MoE models.
- Use for us: Another way to tame ratios without clipping. Useful if we later go async.

---

## Part B. Train-inference mismatch and precision

### 2510.26788 — Defeating the Training-Inference Mismatch via FP16 (October 2025)
- Read: full text
- Core idea: The root cause is BF16 rounding. Switch to FP16 for both trainer and sampler. No algorithm change.
- Key result: FP16 mismatch is about 24x smaller than BF16. Plain importance-weighted policy gradient in FP16 beats all BF16 baselines (GRPO, GSPO, TIS, MIS) in the paper's sanity test. LoRA and larger models are in the tested set.
- Use for us: Cheapest test for the "IS ratio min 0 in 98 % of steps" problem. Run the trainer and vLLM in FP16 and compare the ratio logs.

### 2609.20807 — Score Centering Stabilizes Off-policy Reinforcement Learning (September 2026)
- Read: abstract + tables
- Core idea: Mismatch instability comes from a drift bias that builds each step. An additive "score centering" term cancels it.
- Key result: 0.6B to 30B models. It matches or beats importance sampling under quantization. Gap grows with mismatch (Qwen3-30B-A3B-Base, INTELLECT-2 math data).
- Use for us: Works when mismatch cannot be removed. It composes with importance sampling. Try after the FP16 test.

### 2601.18150 — FP8-RL: A Practical and Stable Low-Precision Stack for LLM Reinforcement Learning (January 2026)
- Read: full text
- Core idea: FP8 rollout (weights and KV cache) in veRL with vLLM, plus token-level TIS or MIS correction.
- Key result: Up to 44 % rollout throughput gain. FP8 without TIS loses accuracy against the BF16 baseline. FP8 with TIS matches it.
- Use for us: Confirms that a rollout correction is required when the sampler is quantized. Speed gain is for long contexts. Our sequences are short, so expect less.

### 2510.18855 — Every Step Evolves: Scaling Reinforcement Learning for Trillion-Scale Thinking Model (October 2025; Ring-1T)
- Read: full text
- Core idea: IcePop masks and clips tokens whose train-inference discrepancy is too large. C3PO++ splits long rollouts under a token budget. ASystem is the RL framework.
- Key result: Ring-1T: AIME-2025 93.4, HMMT-2025 86.72, CodeForces 2088, ARC-AGI-1 55.94.
- Use for us: IcePop is a token-level mask for the mismatch. It is the closest fix for tokens with ratio 0.

### 2506.13585 — MiniMax-M1: Scaling Test-Time Compute Efficiently with Lightning Attention (June 2025)
- Read: full text
- Core idea: CISPO clips the importance weights, not the token updates, so every token keeps its gradient. Found a training-inference precision gap at the LM head and fixed it with FP32 there.
- Key result: On Qwen2.5-32B, CISPO gives a 2x speedup over DAPO.
- Use for us: Two cheap items: CISPO as the loss, and an FP32 LM head in both trainer and vLLM.

### 2605.14220 — Diagnosing Training Inference Mismatch in LLM Reinforcement Learning (May 2026) [new]
- Read: full text (targeted)
- Core idea: VeXact is a rollout engine with zero mismatch. It shows that small token-level disagreement alone can cause collapse. Compares remedies.
- Key result: Experiments on Qwen3-1.7B GRPO and a R1D-1.5B sanity test. Plain recompute or bypass fail while the exact engine works. I did not extract a single headline number.
- Use for us: Tells us to treat our 98 % min-0 finding as a first-order bug, not noise. Gives a method: build a zero-mismatch reference and compare.

### 2602.19362 — LLMs Can Learn to Reason Via Off-Policy RL (February 2026) [new]
- Read: abstract + tables
- Core idea: OAPL accepts off-policy data from a lagged inference policy. It uses no importance ratios.
- Key result: Beats GRPO with importance sampling on competition math. Matches DeepCoder on LiveCodeBench with 3x fewer generations. Works with lags above 400 gradient steps.
- Use for us: If the ratio is unusable (min 0), this removes the ratio altogether.

---

## Part C. Scaling recipes and infrastructure

### 2510.13786 — The Art of Scaling Reinforcement Learning Compute for LLMs (October 2025)
- Read: full text
- Core idea: Fit sigmoid compute-performance curves. ScaleRL combines PipelineRL, forced length interruption, CISPO, prompt-level loss average, batch-level advantage normalization, FP32 logits, zero-variance filtering and no-positive-resampling.
- Key result: More than 400,000 GPU-hours study. ScaleRL asymptotic reward A = 0.61. One run scaled to 100,000 GPU-hours and the curve predicted it. CISPO and GSPO give a higher asymptote than DAPO. FP32 at the LM head raises the asymptote.
- Use for us: Best single recipe checklist. Most items are small code changes in TRL. Also gives a way to compare our runs by curve fit, not by one step.

### 2505.24298 — AReaL: A Large-Scale Asynchronous Reinforcement Learning System for Language Reasoning (May 2025)
- Read: full text
- Core idea: Fully decouple generation and training. Bound staleness. Use a staleness-aware PPO variant.
- Key result: Up to 2.77x training speedup over synchronous systems with the same GPU count, with equal or better final score. The paper says the gain is smaller at short context length.
- Use for us: Our idle GPU is the sync pattern. But we have one GPU, so full decoupling does not fit. The short-context caveat applies to us.

### 2509.19128 — Faster On-policy Reinforcement Learning for Long Sequence Generation (September 2025; PipelineRL)
- Read: full text
- Core idea: Generation and training run at once. In-flight weight updates change the weights during a running generation.
- Key result: About 2x faster learning than a conventional baseline. The abstract says 128 H100 GPUs. The introduction says 4 DGX-H100 nodes. Data stays highly on-policy.
- Use for us: Pattern for removing idle time while staying near on-policy. Needs more than one GPU or a clean split of the one GPU.

### 2506.05433 — Prefix Grouper: Efficient GRPO Training through Shared-Prefix Forward (June 2025)
- Read: full text
- Core idea: Encode the shared prefix once per group. Gradients are identical to standard GRPO.
- Key result: Same outputs and gradients as GRPO. Lower GPU memory at every tested group size. No speed number extracted.
- Use for us: Already in our run as the shared-prefix forward. This paper is the reference for the equivalence proof.

### 2509.23232 — SPEC-RL: Accelerating On-Policy Reinforcement Learning via Speculative Rollouts (September 2025)
- Read: full text
- Core idea: Reuse the previous epoch's rollout as a speculative prefix. Verify it with the new policy.
- Key result: Rollout time cut by 2 to 3x with no loss in score (Qwen-3-8B-Base, PPO, GRPO, DAPO).
- Use for us: We revisit the same crops across epochs. Our outputs are short, so check the gain before building.

### 2509.21792 — FastGRPO: Accelerating Policy Optimization via Concurrency-aware Speculative Decoding and Online Draft Learning (September 2025)
- Read: abstract + tables
- Core idea: Speculative decoding adapted to the high concurrency of GRPO sampling, with an online-trained draft model.
- Key result: End-to-end speedup of 2.35x to 2.72x on math datasets.
- Use for us: Needs a draft model. Poor fit for a 2B model with short outputs. Low priority.

---

## Part D. RL for vision grounding and segmentation

### 2503.06520 — Seg-Zero: Reasoning-Chain Guided Segmentation via Cognitive Reinforcement (March 2025)
- Read: full text
- Core idea: A reasoning model outputs a box and points. A frozen segmenter makes the mask. GRPO with format and accuracy rewards, no reasoning labels.
- Key result: Seg-Zero-7B reaches 57.5 on ReasonSeg zero-shot, 18 % above LISA-7B. Soft reward gives a small gain on ReasonSeg and a clear loss on RefCOCOg. A strict format reward helps out-of-domain.
- Use for us: Supports a hard-threshold reward next to our soft IoU. Test both. Strict format reward is good for polygon syntax.

### 2504.07954 — Perception-R1: Pioneering Perception Policy with Reinforcement Learning (April 2025)
- Read: full text
- Core idea: GRPO for perception tasks. Reward design matters. It matches several predictions to several ground truths by bipartite matching before reward.
- Key result: Qwen2-VL-2B-Instruct: +4.2 % RefCOCO+, +17.9 % PixMo-Count, +4.2 % PageOCR, 31.9 AP on COCO2017 val.
- Use for us: Same model size as ours. Matching is the right pattern for scoring against a set of doctors.

### 2505.12081 — VisionReasoner: Unified Reasoning-Integrated Visual Perception via Reinforcement Learning (May 2025)
- Read: abstract + tables
- Core idea: One model for detection, segmentation and counting. A unified reward and multi-object learning.
- Key result: Relative to Qwen2.5-VL: +29.1 % COCO, +22.1 % ReasonSeg, +15.3 % CountBench.
- Use for us: Confirms that rewards on box and point accuracy transfer. Little new for polygons.

### 2505.22596 — SAM-R1: Leveraging SAM for Reward Feedback in Multimodal Segmentation via Reinforcement Learning (May 2025)
- Read: abstract + tables
- Core idea: Use SAM to produce the mask that is scored as the reward.
- Key result: Strong results with only 3k training samples. Exact numbers not extracted.
- Use for us: Shows RL needs few samples. Our exact CUDA checker plays the role of SAM here and is cheaper.

### 2506.22624 — Seg-R1: Segmentation Can Be Surprisingly Simple with Reinforcement Learning (June 2025)
- Read: full text
- Core idea: GRPO trains a model to prompt SAM2. Reward is format plus mask quality.
- Key result: 0.873 S-measure on COD10K. 71.4 cIoU on RefCOCOg test and 56.7 gIoU on ReasonSeg test, zero-shot. S-measure alone led to reward hacking. The fix: 0.7 IoU + 0.3 S-measure.
- Use for us: Direct warning for our multi-term reward. Watch for hacking in the boundary F and Betti terms. Keep IoU as the main weight.

### 2509.06321 — Text4Seg++: Advancing Image Segmentation via Generative Language Modeling (September 2025)
- Read: full text
- Core idea: Write the mask as text. Run-length coding shortens it.
- Key result: Row-wise RLE cuts descriptor length by 74 % and speeds inference 3x with no loss.
- Use for us: Gives a text-output precedent. Shows output length cuts cost. Our polygons are already short.

### 2512.24023 — RSAgent: Learning to Reason and Act for Text-Guided Segmentation via Multi-Turn Tool Invocations (December 2025)
- Read: abstract + tables
- Core idea: A multi-turn agent calls a segmentation tool, sees the result and revises. SFT cold start, then agentic RL.
- Key result: 66.5 % gIoU on ReasonSeg test (+9 % over Seg-Zero-7B). 81.5 % cIoU on RefCOCOg.
- Use for us: Multi-turn refinement is costly. Read together with 2606.13156 before trying.

### 2606.13156 — Iterative Visual Thinking and the Self-Correction Mirage in VLM Grounding (June 2026)
- Read: full text
- Core idea: A box-refinement loop looks better only because the score keeps the best step per sample, which uses the answer. With fixed label-free stopping, the gain vanishes.
- Key result: Qwen3-VL-4B: stopping at step 0 matches the base and beats every shippable stopping rule. The gain size is not legible in my extract.
- Use for us: Do not add a self-correction loop to fix shifted crops. Score only with a fixed rule that does not use the label.

### 2608.13549 — Exponential Convex Calibration Dimension for the Multi-Label Jaccard Measure (August 2026)
- Read: abstract + tables
- Core idea: Theory. An exactly calibrated convex surrogate for Jaccard (IoU) needs exponentially many outputs. Approximations in polynomial size exist.
- Key result: Proofs only. No benchmark number.
- Use for us: Supports using the exact IoU as a reward instead of a smooth surrogate loss. No code change.

---

## Part E. New 2026 additions (not in the given list)

### 2601.19228 — Towards Pixel-Level VLM Perception via Simple Points Prediction (January 2026; SimpleSeg)
- Read: full text
- Core idea: A plain MLLM writes the mask outline as a sequence of point coordinates in text. SFT first, then RL with GSPO and an IoU reward.
- Key result: Reward table (RefCOCO / RefCOCO+ / RefCOCOg gIoU): IoU 76.9 / 71.9 / 72.9. IoU + distance 77.1 / 72.2 / 73.1. Adding a hard length penalty drops to 66.7 / 62.4 / 62.0. Polygon length sweet spot: 221 tokens best, 35.6 cIoU at 78 tokens, 72.5 cIoU at 859 tokens. Pre-training raises gIoU by 4.6 (SFT) and 13.0 (SFT+RL) on RefCOCO.
- Use for us: Closest match to our task (polygon text from a VLM with IoU reward). Warns against a hard length penalty. Shows a point-count sweep we can copy.

### 2605.15951 — From Failure to Feedback: Group Revision Unlocks Hard Cases in Object-Level Grounding (May 2026)
- Read: full text (targeted)
- Core idea: When all answers in a group fail, sample an initial answer, make revised candidates, and use each candidate's gain over the initial answer as a reward-shaping signal and an advantage modulation.
- Key result: Gains over GRPO-based models on referring and reasoning segmentation, REC and counting. No single number extracted.
- Use for us: Aimed at our all-miss problem. A revision gain gives signal when every raw IoU is near zero. Check cost: it adds rollouts.

### 2601.07695 — Smooth Operator: Smooth Verifiable Reward Activates Spatial Reasoning Ability of Vision-Language Models (January 2026)
- Read: abstract + tables
- Core idea: In GRPO, near-miss answers lose their advantage after relative normalization. A smooth sigmoid reward and an absolute-preserving advantage (AP-GRPO) keep the signal.
- Key result: AP-GRPO matches large supervised methods with higher data efficiency on 3D tasks (Numerical3D-50k). No number extracted.
- Use for us: Names our failure: near-zero IoU groups lose signal after normalization. A dense reward plus an absolute scale term is a small change.

### 2605.14220 and 2602.19362
- Listed in Part B (new 2026 additions on the mismatch problem).

---

## Short ranking for this project
1. 2601.19228 SimpleSeg: same task, same output type, IoU reward, GSPO.
2. 2510.26788 FP16: one config change to test the mismatch.
3. 2510.13786 ScaleRL: checklist of 8 changes with ablations.
4. 2505.22617 Entropy Mechanism: Clip-Cov or KL-Cov for entropy collapse.
5. 2605.15951 Group Revision (with 2601.07695 and 2503.20783): signal for all-miss groups.
