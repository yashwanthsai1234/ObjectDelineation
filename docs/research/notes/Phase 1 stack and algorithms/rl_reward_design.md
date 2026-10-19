# RLVR algorithms, reward design, controllers and novelty refresh for polygon-emitting VLM (Phase 1)

Scope: choose the policy-optimization variant, reward normalization/combination, Level 2 controller design and evaluation protocol for GRPO-training a ≤2B VLM that emits segmentation polygons scored by an exact geometry checker against a set of expert annotations; refresh the prior-art check for our novelty claims. Date of research: 2026-10-18.

Source-quality legend used below:
- **[fetched]**: read on the primary page (arXiv abs/html or project page) during this session.
- **[snippet]**: claim comes from a search-result summary only; treat as likely but unverified.
- **[prior knowledge, not re-fetched]**: well-known 2025 paper cited from memory with its arXiv ID. The ID and headline claim are high-confidence, but exact numbers should be checked before they go into a paper.
- PDF summaries from the fetch tool were sometimes generic. Where a summary looked templated I say so and treat it as title-level evidence only.

**Correction to the brief:** arXiv **2601.07695** is *not* an "advantage collapse" paper. It is "Smooth Operator: Smooth Verifiable Reward Activates Spatial Reasoning Ability of Vision-Language Model" (Jiao et al., Jan 2026), which proposes AP-GRPO. The advantage-collapse paper is **arXiv 2605.21125** (ICML 2026). Both are covered below.

---

## Q1. Policy-optimization variant: which default for a ~0.5–2B VLM with dense continuous rewards?

### Takeaway
No 2025–2026 paper runs a head-to-head of GRPO / Dr. GRPO / DAPO / GSPO / GMPO / CISPO on a ≤2B VLM with continuous geometric rewards, so the choice has to be argued from adjacent evidence. That evidence points to a GRPO-family estimator with:
- Dr. GRPO-style constant (not per-sequence) length aggregation;
- DAPO clip-higher;
- a very small or zero KL;
- variance-aware filtering of degenerate groups.

GSPO is the most relevant alternative because SimpleSeg, the closest polygon-as-text work, used it. Its main documented advantage (MoE stability) does not apply to a dense ≤2B model, so it should be an ablation arm, not the default. Std normalization is the biggest live risk for dense multi-term rewards (see Q2).

### Cited Findings
- **ScaleRL** (Meta, 400k+ GPU-hours of ablations): "CISPO and GSPO outperform DAPO in asymptotic performance"; CISPO was chosen as the default for robustness. KL is "often set much smaller or omitted entirely" in reasoning RL — [aweers blog summary, 2026](https://aweers.de/blog/2026/rl-for-llms/) [fetched]. Original paper: [arXiv 2510.13786](https://arxiv.org/abs/2510.13786) [prior knowledge, not re-fetched].
- **CISPO** (MiniMax-M1) clips the importance-sampling *weight* with stop-gradient instead of masking token updates. It is reported as ~2× faster than DAPO — [aweers blog](https://aweers.de/blog/2026/rl-for-llms/) [fetched]; [MiniMax-M1, arXiv 2506.13585](https://arxiv.org/abs/2506.13585) [prior knowledge].
- **Dr. GRPO** found that sample-level loss normalization "introduces a bias favoring short correct responses and long incorrect ones". Its fix is a constant denominator, and it also removes std division — [aweers blog](https://aweers.de/blog/2026/rl-for-llms/) [fetched]; [arXiv 2503.20783](https://arxiv.org/abs/2503.20783) [prior knowledge].
- **DAPO** uses token-level aggregation, asymmetric clipping (ε_low = 0.2, ε_high = 0.28), dynamic sampling (drop groups with all-equal reward), and an overlong-response penalty — [aweers blog](https://aweers.de/blog/2026/rl-for-llms/) [fetched]; [arXiv 2503.14476](https://arxiv.org/abs/2503.14476) [prior knowledge].
- **GSPO** (sequence-level, length-normalized importance ratio) is described as "best for stable MoE RL", with GRPO/CISPO best for short verifiable tasks — [TuringPost "Reasoning RL in 2026"](https://www.turingpost.com/p/reasoning-rl-in-2026) [snippet]; [arXiv 2507.18071](https://arxiv.org/abs/2507.18071) [prior knowledge].
- **SimpleSeg** (Moonshot AI + Nanjing U., arXiv 2601.19228; polygon-as-text-points segmentation) — [HF paper page](https://huggingface.co/papers/2601.19228) [fetched]:
  - RL setup: GSPO, clip range [3e-4, 4e-4], KL coefficient 0.01, 8 responses per group, batch 256, 32 GPUs, Muon optimizer.
  - Backbones: Qwen2.5-VL-7B and Kimi-VL (2.8B activated, MoE).
  - SFT → SFT+RL gains: refCOCO val 65.5 → 75.2 gIoU (+9.7); refCOCO+ +9.8; refCOCOg +10.5.
  - Rationale: "Contours are inherently many-to-one w.r.t. masks; enforcing exact token matching is suboptimal."
  - The arXiv HTML/PDF could not be fetched (404, and the PDF exceeds 10 MB). Their use of GSPO may reflect the MoE backbone.
- **GDPO's own recipe** (Qwen2.5-1.5B tool calling; DeepSeek-R1-7B math/code) — [arXiv 2601.05242](https://arxiv.org/html/2601.05242) [fetched]:
  - KL 0.0005 (MSE-form);
  - clip_low 0.2 / clip_high 0.28;
  - LR 1e-6;
  - 4 rollouts per query (tool calling), 16 (math/code).
- **VLM-R1** (3B/7B/32B) — [arXiv 2504.07615](https://arxiv.org/html/2504.07615) [fetched]:
  - KL β = 0.04 for REC, β = 0 for open-vocabulary detection;
  - group size 8; temperature 0.9.
  - The 3B model did not gain on reasoning-heavy subtasks where 7B/32B did (Negation 38.7 → 37.7 for 3B vs 39.0 → 43.1 for 7B).
- **Topo-R1** (InternVL-2.5-2B through Qwen3-VL-8B) used plain GRPO with G = 4 and a fixed KL to the SFT reference. RL-only collapsed (F1@0.5 = 3.1); SFT-only reached 23.0 and SFT+RL 45.2 — [arXiv 2603.13054](https://arxiv.org/html/2603.13054) [fetched]. Implication: SFT warm-start is mandatory at small scale.
- **AP-GRPO** ("Smooth Operator", arXiv 2601.07695) proposes "Absolute-Preserving GRPO". It maps raw spatial feedback through a dynamically parameterized sigmoid into a dense continuous reward and preserves absolute reward information that pure group ranking discards — [arXiv abs](https://arxiv.org/abs/2601.07695) [fetched, abstract only]. Relevant because our IoU-type rewards are dense and their absolute level matters.
- **Advantage collapse** — [arXiv 2605.21125, ICML 2026](https://arxiv.org/abs/2605.21125) [fetched, abstract]:
  - Defines the Advantage Collapse Rate (ACR): the share of batches with ineffective gradients from homogeneous groups. It predicts stagnation across 0.5B–14B.
  - AVSPO injects *virtual* reward samples when ACR is high, cutting collapse 58–63% and adding +4–6 pp accuracy at all scales tested.
- **Entropy-collapse remedies:**
  - Clip-Cov / KL-Cov act on high-covariance tokens, which drive entropy collapse — [arXiv 2505.22617](https://arxiv.org/abs/2505.22617v1) [snippet].
  - ProRL adds KL to a reference plus *periodic reference resets* (with optimizer re-init) on top of DAPO decoupled clipping, starting from a 1.5B model — [arXiv 2505.24864](https://arxiv.org/pdf/2505.24864) [snippet].
  - "Clip-Low Increases Entropy and Clip-High Decreases Entropy" — [HF 2509.26114](https://huggingface.co/papers/2509.26114.md) [title-level only; note that it qualifies the usual "clip-higher = more exploration" story].
- **Algorithm robustness under imperfect multimodal rewards** (Qwen-VL 2B → 32B) — [arXiv 2607.09492](https://arxiv.org/html/2607.09492) [fetched]:
  - GRPO was the most stable (reward-hacking rate 48–53%).
  - RLOO was the most vulnerable (67–68%).
  - DAPO was strongly scale-dependent (67.2% at 2B → 45.5% at 8B).
  - The authors recommend "Prefer GRPO under reward uncertainty". This bears directly on a ≤2B model: DAPO's machinery looked worse at 2B in their setting.
- **Small models and exploration:** "Smaller Models are Natural Explorers for Policy-Level Diversity in GRPO" claims smaller models (0.5B–7B tested) keep higher policy-level diversity — [arXiv 2605.30789](https://arxiv.org/pdf/2605.30789) [fetched, but the summary was generic; treat as title-level].
- **REINFORCE++** (global batch normalization instead of per-group) — [arXiv 2501.03262](https://arxiv.org/abs/2501.03262) [prior knowledge].
- **GMPO** (geometric mean of token ratios for stability) — [arXiv 2507.20673](https://arxiv.org/abs/2507.20673) [prior knowledge].
- **RLOO** (leave-one-out baseline) — [arXiv 2402.14740](https://arxiv.org/abs/2402.14740) [prior knowledge].
- None of these three was benchmarked against GRPO on VLM geometry in anything I found.

### Inferences
- For a dense ≤2B model with one sequence-level reward per polygon, GSPO's sequence ratio matches the reward granularity. Its documented win is MoE stability, which we do not need. GRPO with token-level aggregation and Dr. GRPO length handling is the safer default, given the 2B robustness result in 2607.09492 and Topo-R1's 2B success with plain GRPO.
- Polygon outputs have variable length: more vertices means more tokens. Per-sequence length normalization (original GRPO) biases toward short correct polygons and long wrong ones, which is exactly the Dr. GRPO bias. This could push the model toward under-vertexed polygons, so use a constant denominator.
- DAPO's "dynamic sampling" assumes binary rewards. With continuous IoU, exact-zero-variance groups are rare, but *near*-zero-variance groups (all IoU ≈ 0.93) are common. Filter on a std threshold, not equality. AVSPO-style virtual samples or AP-GRPO-style absolute terms are alternatives.
- KL: the evidence spans β = 0 (VLM-R1 OVD), 0.0005 (GDPO), 0.01 (SimpleSeg) and 0.04 (VLM-R1 REC). With SFT warm-start and a small model, a small fixed β (1e-3 to 1e-2) plus the option of a ProRL-style reference reset is a reasonable band. The controller can move within it.

### Gaps
- No 2025–2026 controlled comparison of GRPO vs GSPO vs CISPO vs Dr. GRPO on ≤2B **VLM** perception tasks was found. ScaleRL's ablations are on LLM reasoning at 8B+ scale.
- SimpleSeg's full RL section (why GSPO, exact IoU-reward form) could not be read. The HF page's "Mask IoU reward (range [0.0, 0.1])" looks garbled and is unverified.
- GMPO, REINFORCE++ and RLOO numbers were not re-verified this session.

---

## Q2. Multi-objective reward combination, scaling, partial credit, and reward hacking in segmentation/grounding RL

### Takeaway
Per-term group normalization (GDPO) beats summed-then-normalized GRPO when terms differ in scale and difficulty. Its own authors show, however, that weights alone cannot fix difficulty imbalance: you must **gate easy terms on hard ones**. Separately, 2026 evidence shows that z-scoring a dense, low-variance auxiliary term inside groups can catastrophically amplify it (the "dark room" collapse). Our legality, hole and topology terms are exactly that kind of term. So we need GDPO + gating + a std floor (or mean-only centring) on low-variance terms.

Every segmentation/grounding RL paper reports a geometry-specific hack, and each fix was a multiplicative penalty or gate:
- over-prediction (VLM-R1 odLength);
- oversized boxes (Topo-R1 LocPen);
- S-measure-only hacking (Seg-R1).

### Cited Findings
- **GDPO formulation and results** — [arXiv 2601.05242](https://arxiv.org/html/2601.05242) [fetched]:
  - Steps: (1) z-score each reward within its group; (2) weighted sum; (3) batch-level z-score of the sum.
  - Problem it solves: under GRPO, distinct reward combinations collapse to identical advantages, e.g. (0,1) and (0,2) — [HF page](https://huggingface.co/papers/2601.05242) [snippet].
  - Tool calling (Qwen2.5-1.5B, BFCL-v3): 32.81% vs 30.18% accuracy; format 80.66% vs 76.33%.
  - Math (R1-7B AIME): 53.1% vs 50.2%, length violations 0.2% vs 2.1%.
  - Weights alone fall short: "Reducing w_length from 1.0 to 0.25 yields only modest changes" without conditioning.
  - Their recommendation is to *condition* easier rewards on harder ones (e.g. length reward = 1 only if correct). With conditioning, GDPO gains +4.4% AIME with 16.9% fewer violations vs GRPO.
- **Dark-room collapse** from z-scored dense rewards — ["The Dark Room in the Reward Channel", arXiv 2607.21273 (Yu Wang, Aug 2026)](https://www.alphaxiv.org/overview/2607.21273) [fetched overview]:
  - Design: 74 preregistered arms; ALFWorld/WebShop/POMDP; Qwen 1.7B/4B/8B.
  - Mechanism: in all-fail groups, σ is tiny and "z-scoring cancels the shaping coefficient". A small auxiliary reward therefore becomes order-1 advantage.
  - Outcome: every standard-GRPO-normalized run collapsed; removing std normalization alone restored baseline.
  - Fixes: mean-only normalization, gate dense rewards on task success, monitor within-group variance.
  - A [ModelWire article](https://themodelwire.com/article/grpo-normalization-breaks-dense-prediction-rewards-for-llm-agents-01KY8YQ85PJEJRTXTV4G6XJR7F) [snippet] reports a 0.1 reward amplified to advantage 5–7.
- **Counter-view: std normalization can help.** It acts as an adaptive, curvature-based step size and gives strictly better convergence than unnormalized REINFORCE under conditions tied to within-prompt reward std. The benefit fades in late training — ["Why GRPO Needs Normalization", arXiv 2601.23135 (Ge et al., Jan 2026; v2 Oct 2026)](https://arxiv.org/abs/2601.23135v1) [fetched abstract].
- **Topo-R1 composite reward** — [arXiv 2603.13054](https://arxiv.org/html/2603.13054) [fetched]:
  - Weights: R = 0.10·R_fmt + 0.85·R_acc + 0.05·R_topo, with R_fmt a binary gate.
  - Partial credit: R_acc uses a smooth piecewise φ(IoU) with tiers {0.3, 0.5, 0.7, 0.9} → ceilings {0.25, 0.55, 0.80, 1.0}, exponent γ = 1.5.
  - Topology term: per matched pair (1 − clDice)·LocPen(b), as extracted by the fetch tool; this sign is unusual and should be checked in the paper.
  - LocPen = 1 if box area ≤ 0.3 of the image, else max(0, 1 − 0.8·(a − 0.3)/0.7). It suppresses "oversized boxes that trivially capture topological changes".
  - Ablation (Qwen3-VL-4B, F1@0.5): fmt only 25.1; fmt+acc 43.8; fmt+topo 31.4; all three 45.2.
- **VLM-R1 over-prediction hack** — [arXiv 2504.07615](https://arxiv.org/html/2504.07615) [fetched]:
  - Hack: an AP reward on OVD let the model list every category, because absent GT categories are excluded from COCO eval. Output length rose sharply.
  - Fix: multiplicative odLength, s = min(1, L_gt/L_pred)·mAP. COCO-filtered mAP 21.1 vs 11.8; OVDEval nms-AP 31.01 vs 21.68.
- **Seg-R1:** S-measure alone caused reward hacking, so the segmentation reward became 0.7·IoU + 0.3·S-measure — [arXiv 2506.22624](https://arxiv.org/html/2506.22624) [snippet].
- **SAM-R1 / Seg-Zero:** SAM-R1 notes that Seg-Zero's full decoupling of reasoning model and mask decoder "increases the risk of reward hacking" — [arXiv 2505.22596](https://arxiv.org/pdf/2505.22596) [snippet].
- **Multimodal reward hacking taxonomy** — [arXiv 2607.09492](https://arxiv.org/html/2607.09492) [fetched]:
  - Hack types: decision, evidence and reward-form hacking. Ambiguous boundaries raise the hacking rate by +5.1 to +9.0 pp.
  - Metrics: reward-hacking rate (RHR), "Newly Rewarded Failure Rate" (NRFR, an alarm when NRFR > RHR) and reward–oracle gap.
  - Unreliable sub-verifiers (keyword matching) *increased* hacking; "hard gates" are recommended over soft lexical bonuses.
- **Field-routed advantages:** GRCA routes each geometric field's standardized reward only to that field's tokens instead of broadcasting a summed reward. On Qwen2.5-VL-3B this raised 3D keypoint accuracy 0.64 → 0.93 while keeping IoU2D 0.894 — [arXiv 2604.21160 (Chen et al., Apr 2026)](https://arxiv.org/html/2604.21160) [fetched].
- **Partial credit vs binary:** AP-GRPO/"Smooth Operator" argues that a dense sigmoid-shaped continuous reward beats thresholded or ranked rewards for VLM spatial reasoning — [arXiv 2601.07695](https://arxiv.org/abs/2601.07695) [fetched abstract]. Topo-R1's smooth φ tiers make the same choice.

### Inferences
- **Our likely failure modes** map directly onto these findings:
  - **Legality** will be ~1.0 for most rollouts after SFT, so its within-group σ ≈ 0. GDPO's per-term z-score then turns sub-1% legality differences into order-1 advantages, the dark-room mechanism.
  - **Hole matching and topology** behave similarly on simply connected targets, where every rollout scores the same.
  - **Fix:** per-term σ floor (e.g. σ_k ← max(σ_k, σ_min,k) with σ_min set from the term's natural scale), or mean-only centring for low-variance terms, or skip a term's advantage when its group σ < ε.
- **Gating** is supported by GDPO's conditioning result, Topo-R1's binary format gate and LocPen, and the dark-room fix ("gate dense rewards on task success"). Our plan to gate topology on overlap is in line with this. Extend the same treatment to the egg–yolk band term, since being "inside the band" is trivially easy for a huge blob.
- **The guard we need** is the analogue of LocPen/odLength: a multiplicative penalty on area inflation or part-count inflation relative to the rater union. This blocks the "cover everything to satisfy union-band and any-rater topology" hack. Over-prediction of parts is the polygon analogue of VLM-R1's over-detection.
- **GRCA-style routing** is only partly applicable. IoU depends on all vertices, but per-ring or per-part terms (hole matching, per-part legality) *could* be routed to that ring's token span. This is a cheap ablation once the output format delimits rings.

### Gaps
- GDPO does not, in the fetched text, say how it handles zero-σ groups for a single term. Check the GitHub/verl implementation for an epsilon or skip.
- I found no paper that compares summed vs per-term-normalized vs multiplicative-gated rewards specifically for segmentation/polygon RL.
- SAM-R1's specific guards were not read in full.

---

## Q3. Adaptive controllers / curricula in RLVR 2025–2026, and whether anyone already uses typed verifier feedback to steer training

### Takeaway
Difficulty-based prompt selection is mature and works. It includes:
- online difficulty filtering (target p ≈ 0.5, maximize p(1−p));
- DAPO dynamic sampling;
- PODS max-variance down-sampling;
- AdaRFT and METIS curricula;
- ACR-triggered virtual samples.

Several 2026 papers *audit* verifier errors by category, and one uses a failure-mode taxonomy to *redesign rewards between runs* with an LLM. I found **no paper that reads typed, per-error-category verifier logs online during RL to adapt sampling weights, KL and stopping**. The Level 2 claim survives, but its parts all have precedents, so it should be framed as a typed extension of learnability-based sampling plus alarms, not as a new class of method.

### Cited Findings
- **Online Difficulty Filtering:**
  - Learnability ∝ p(x)(1−p(x)), maximized at p = 0.5. Balanced filtering gave up to +12% in under half the GRPO steps.
  - Venue: EACL 2026 long — [alphaXiv 2504.03380](https://www.alphaxiv.org/abs/2504.03380); [ACL Anthology](https://aclanthology.org/2026.eacl-long.30.pdf) [snippet].
- **PODS** (TMLR 04/2026) — [arXiv 2504.13818](https://arxiv.org/html/2504.13818v2) [snippet]:
  - Max-variance down-sampling of rollouts, in O(n log n).
  - For binary rewards it takes the top m/2 and bottom m/2 rollouts.
  - GRPO+PODS reaches vanilla GRPO's peak accuracy ≥1.7× faster.
  - It applies directly to continuous rewards.
- **Further difficulty predictors and schedulers** (2025–2026):
  - "Can Prompt Difficulty be Online Predicted" — [arXiv 2507.04632](https://arxiv.org/pdf/2507.04632)
  - "Efficient RLVR Scheduling via Graph-Structured Online Difficulty Estimation" — [arXiv 2608.17941](https://arxiv.org/pdf/2608.17941)
  - "Efficient Dynamic Prompt Selection" — [arXiv 2607.27610](https://arxiv.org/pdf/2607.27610)
  - "ThinkPrior: zero-rollout difficulty priors" — [arXiv 2609.09075](https://arxiv.org/pdf/2609.09075)
  - All four are title-level only [snippet].
- **METIS** ("Internalizing Curriculum Judgment") claims up to 67% faster convergence on math, code and function calling — [arXiv 2605.11235](https://arxiv.org/html/2605.11235v1) [snippet].
- **FG-EXPO** combines adaptive KL with a Gaussian curriculum — [arXiv 2605.11403](https://arxiv.org/pdf/2605.11403) [title-level snippet]. It is a direct precedent for "adaptive KL + curriculum" in one controller.
- **AdaRFT** (adaptive curriculum that targets a reward level) — [arXiv 2504.05520](https://arxiv.org/abs/2504.05520) [prior knowledge, not re-fetched].
- **ACR → AVSPO:** a closed-loop controller on one *scalar* diagnostic (advantage-collapse rate) that injects virtual samples — [arXiv 2605.21125](https://arxiv.org/abs/2605.21125) [fetched].
- **Typed-error papers found (closest to Level 2):**
  - **"When LLM Reward Design Fails: Diagnostic-Driven Refinement for Sparse Structured RL"** (Wang et al., TMLR 2026) — [arXiv 2605.28918](https://arxiv.org/abs/2605.28918) [fetched]:
    - Uses training diagnostics plus a failure-mode taxonomy (reward flooding, semantic/API misunderstanding, weak shaping) to prompt an LLM to *revise the reward function between runs*.
    - Scale: PPO on MiniGrid (DoorKey 2.3% → 97.6%).
    - Limits: no online curriculum, KL or stopping control. **Closest precedent in spirit**, but offline, LLM-in-the-loop and reward-side.
  - **"Where the Verifier Fails: A Category-Level Audit of Reward Signals in RLVR"** (Sept 2026) — [arXiv 2609.01354](https://arxiv.org/html/2609.01354) [fetched]:
    - Metamorphic testing over 307,420 verdicts; 93% of failures are whitespace/punctuation.
    - **Audit only**; no training adaptation.
  - **"Verifier Errors in RLVR: Reward Hacking, Limits of Feedback, and Selective Control"** (Moya, Thornley, Lin, Sept 2026) — [alphaXiv 2609.35677](https://www.alphaxiv.org/abs/2609.35677) [fetched]:
    - Projected Audit Correction audits accepted answers and projects corrections orthogonal to the reward gradient.
    - Qwen2-0.5B: 2.2% → 95.2% correctness with 25% audit coverage.
    - **Not typed**; does not adapt sampling, KL or stopping.
  - **"Delay, Plateau, or Collapse"** — [arXiv 2605.02909](https://arxiv.org/html/2605.02909) [snippet]:
    - The *pattern* of systematic verifier errors, not the aggregate rate, decides whether training is delayed, plateaus or collapses.
    - This supports logging per-type rather than aggregate.
  - **"Are Verifier Errors Independent Within a GRPO Group?"** — [arXiv 2609.06386](https://arxiv.org/pdf/2609.06386) [snippet]: within-group verifier error is correlated and structured by answer-form category.
  - **CURE** (curriculum-guided multi-task training for anatomy-grounded report generation) — [arXiv 2601.15408](https://arxiv.org/pdf/2601.15408) [snippet]:
    - Reweights sampling per dataset and per anatomical class by error rate.
    - It is **supervised multi-task**, not RL, but it is a precedent for "per-category error rate → sampling weight".
- **Early-stopping / overfitting signals:**
  - **Diversity collapse as overtraining.** "Understanding Diversity Collapse in RLVR via the Lens of Overtraining" (Yuan et al., June 2026) — [arXiv 2606.15455](https://arxiv.org/pdf/2606.15455) [fetched, summary generic].
    - Frames diversity collapse as overtraining.
    - Lists entropy decline, pass@k crossover and a falling count of unique answers as signals. Treat that list as title/abstract-level.
  - **Newly-created failures.** NRFR rising above RHR flags that RL is actively creating failures — [arXiv 2607.09492](https://arxiv.org/html/2607.09492) [fetched].
  - **Implicit reward overfitting.** "On the Implicit Reward Overfitting and the Low-rank Dynamics in RLVR" — [arXiv 2605.06523 via fugumt](https://fugumt.com/fugumt/paper_check/2605.06523v1) [snippet, Japanese summary site]:
    - Test performance can be good while training reward stays low.
    - It reports rank-1 concentration of the RL update.
    - Implication: train reward is not a stopping signal.
  - **Dark-room monitoring.** The dark-room paper recommends monitoring within-group variance trajectories — [arXiv 2607.21273](https://www.alphaxiv.org/overview/2607.21273) [fetched].

### Inferences
- **What Level 2 should reuse:**
  - p(1−p)-style learnability per item, generalized to continuous rewards as within-group reward variance;
  - PODS-style rollout down-sampling;
  - ACR as a health metric;
  - NRFR and a pass@k/unique-topology crossover as alarms;
  - ProRL reference reset as a KL "action".
- **What is plausibly new:** a *vector* learnability signal (per error type: legality, topology mismatch vs all raters, hole miss, band violation) used to:
  - reweight sampling toward items whose *dominant error type* is still learnable;
  - adjust per-term weights or gates and KL per phase;
  - stop on type-specific overfitting (e.g. train topology-error rate falls while held-out any-rater topology rate does not).
- **Reviewer framing:** 2605.28918 (taxonomy-driven reward redesign) and CURE (per-class error-rate reweighting) are the two citations a reviewer will raise. Distinguish by: online, inside one RL run, verifier-native types from an exact checker, no LLM in the loop.

### Gaps
- PODS, AdaRFT, METIS and FG-EXPO details beyond abstracts were not read. FG-EXPO's adaptive-KL rule in particular should be checked for overlap.
- I found no 2025–2026 paper on *adaptive KL controllers* specific to RLVR besides FG-EXPO (title) and ProRL resets. The classical PPO adaptive-KL target controller (Ziegler et al. 2019) is the obvious baseline [prior knowledge].
- I could not run an exhaustive search of OpenReview NeurIPS 2026 papers for "typed feedback curriculum". A paper may exist under different wording, e.g. "skill-aware", "error-conditioned" or "diagnostic curriculum".

---

## Q4. Ambiguity-aware learning from multiple annotators; does GRPO preserve or collapse a multimodal answer distribution?

### Takeaway
Reward-maximizing GRPO is mode-seeking: there is strong 2025–2026 evidence that it collapses to the dominant solution. So "sampled topologies match the readers' distribution" will **not** come for free from an any-rater reward. If any rater's topology earns full reward, GRPO will converge on whichever topology is easiest or most frequent and drop the others.

Fixes exist:
- distribution-matching objectives (FlowRL);
- mode-covering divergences;
- pass@k-style advantages;
- set-valued multi-answer rewards (Multi-Answer RL, MIT 2026), the closest analogue to our "any-rater" idea, though for text QA.

The segmentation side (Prob-U-Net, PHiSeg, diffusion multi-rater segmenters) learns rater distributions with supervised/likelihood objectives. I found **no RL work that rewards matching ANY rater's geometry or topology**.

### Cited Findings
- **Mode collapse under standard GRPO.** GRPO "progressively narrows its output distribution as it reinforces the dominant solution pattern"; pass@1 rises while high-k pass@k degrades — [search summary of 2601.22478 and 2509.07430](https://arxiv.org/html/2601.22478v1) [snippet].
- **Mode-covering divergences.** Forward-KL / JS-divergence regularization (DPH-JS) "preserv[es] a multi-modal distribution" and improves both pass@1 and pass@k — ["The Choice of Divergence", arXiv 2509.07430](https://arxiv.org/html/2509.07430v1) [snippet].
- **Limit of RLVR at large k.** GRPO is overtaken by the base model at large k (k = 64: 86.8 vs 87.7). Pass@k Training derives an analytical pass@k advantage — [arXiv 2508.10751](https://arxiv.org/abs/2508.10751) [snippet]. Original "limit of RLVR" paper: [arXiv 2504.13837](https://arxiv.org/abs/2504.13837) [prior knowledge].
- **FlowRL** (ICLR 2026) — [arXiv 2509.15207](https://arxiv.org/html/2509.15207v1); [ICLR proceedings](https://proceedings.iclr.cc/paper_files/paper/2026/file/f657eb3343d4bb85a1d77821d1fbe4b8-Paper-Conference.pdf) [snippet]:
  - Turns scalar rewards into a normalized target distribution with a learnable partition function and minimizes reverse KL to it.
  - "Reward-maximizing methods like GRPO concentrate on single high-reward peaks while FlowRL … matches the full reward distribution".
- **Multi-Answer RL** ("Reaching Beyond the Mode", Puri, Damani, Shenfeld, Ghassemi, Andreas, Kim, MIT) — [project page](https://multi-answer-rl.github.io/); arXiv 2603.24844 [fetched]:
  - Rewards: coverage R = Σᵢ 𝟙[aᵢ ∈ Y*] over K answers in one generation, against a *set* of valid answers Y*. The RLCR variant subtracts a Brier calibration penalty.
  - Results: standard RL mode-collapses, while Multi-Answer RL gives ~2× more unique answers (8 vs 4 on medical tasks).
  - Scope: text QA/medical/code only; no structured outputs, no inter-annotator disagreement.
  - The project page says "(2025)" but the arXiv ID implies March 2026; the date is unverified.
- **DARL** encourages diverse answers for general reasoning without verifiers — [arXiv 2601.14700](https://arxiv.org/html/2601.14700) [title-level snippet].
- **Vector Policy Optimization.** "Training for Diversity Improves Test-Time Search" — [arXiv 2605.22817](https://arxiv.org/html/2605.22817v1) [title-level snippet].
- **Measuring diversity beyond pass@k.** "What pass@k Cannot Measure: Evaluating Diversity and Capability Retention after Post-Training" — [arXiv 2610.07405](https://arxiv.org/html/2610.07405) [title-level snippet]. Relevant to how we evaluate topology diversity.
- **Supervised multi-rater segmentation (2025–2026):**
  - Multi-Rater Calibrated Segmentation treats voxel-wise annotator agreement as an ordinal target with a Ranked-Probability-Score loss — [arXiv 2605.02437](https://arxiv.org/html/2605.02437v1) [snippet].
  - DiffOSeg (multi-expert collaboration diffusion) — [arXiv 2507.13087](https://arxiv.org/pdf/2507.13087) [title-level].
  - Harmonized multi-rater personalization — [arXiv 2605.08210](https://arxiv.org/pdf/2605.08210) [title-level].
  - None uses RL.
- **Classic multi-rater / set-valued baselines** [prior knowledge, not re-fetched]:
  - Probabilistic U-Net — [arXiv 1806.05034](https://arxiv.org/abs/1806.05034)
  - PHiSeg — [arXiv 1906.04045](https://arxiv.org/abs/1906.04045)
  - Stochastic Multiple Choice Learning (winner-takes-all / min-over-hypotheses loss) — [arXiv 1606.07839](https://arxiv.org/abs/1606.07839)
  - These use min-over-raters or ELBO objectives, i.e. the supervised analogue of "any-rater".
- **Small models may stay more diverse.** "Smaller models are natural explorers for policy-level diversity in GRPO" — [arXiv 2605.30789](https://arxiv.org/pdf/2605.30789) [title-level]. This might delay collapse at ≤2B but does not prevent it.

### Inferences
- **Why collapse is expected.** Our any-rater topology term is a max-over-raters reward, the RL analogue of winner-takes-all. Under reward maximization, the policy's optimum is a point mass on the most reliably achieved rater mode. The hypothesis "sampled topologies follow the readers' distribution" is therefore *contrary* to the 2025–2026 evidence unless we add a distributional term.
- **Candidate distributional terms**, in rough order of cost:
  - (a) a group-level coverage bonus: the fraction of distinct rater topologies matched by the G rollouts, in the style of Multi-Answer RL or pass@k advantages;
  - (b) a FlowRL-style target distribution in which each rater mode gets reward mass proportional to its rater frequency;
  - (c) a forward-KL/JS regularizer to the SFT policy if SFT was trained on all raters.
- **Set-valued reward precedent.** Multi-Answer RL is the strongest precedent for a "set of valid answers Y*" reward. Cite it as neighbouring for the set-valued topology term, while noting it is text, single-sample-K-answers, and not rater-derived.
- **Evaluation implication.** We need a *distributional* metric: KL or total variation between the sampled topology histogram and the rater topology histogram per image, plus a GED-style metric for geometry. Pass@k alone (2610.07405) is insufficient.

### Gaps
- No RL paper found that rewards matching any of several annotators' segmentations, or that measures whether RL-sampled segmentations reproduce the rater distribution. This is both a novelty opportunity and a risk, since there is no prior evidence it works.
- FlowRL and pass@k-training results are on math/code LLMs; transfer to dense geometric rewards is untested.

---

## Q5. Novelty refresh: threats to (a) exact β0/β1 of the policy's own output as reward, (b) rater-derived spatially varying band, (c) topology scored against a rater set, (d) typed-error-driven controller

### Takeaway
No new paper found today directly anticipates (a), (b) or (c). New neighbours to cite:
- **CurveBench** (RLVR on topological containment trees, but the model outputs a *tree as text*, not geometry);
- **PolyTopoBench** (NeurIPS 2026 E&D, ring/hole-aware polygon metrics, no RL);
- **VectorLLM++** (RLVR "directly optimizing geometric quality" for polygons; reward terms still unreadable, so it remains the main unresolved threat to (a));
- **Multi-Answer RL** (set-valued reward, text) for (c);
- the **typed-error papers in Q3** for (d).

"Learning to Trace" (NeurIPS 2026) still could not be located.

### Cited Findings
- **VectorLLM++** ("A Unified Next-Token-Prediction MLLM for Dense Remote Sensing Perception and Vector Mapping") — [Research Square rs-10949232](https://www.researchsquare.com/article/rs-10949232/v1) [snippet; full text not readable via fetch]:
  - Seven tasks including polygon extraction.
  - Data and training: 62.9M pretraining samples from 57 datasets, 7.2M SFT instructions, then "reinforcement learning with verifiable rewards that directly optimize geometric quality".
  - 1024-token vector coordinate vocabulary.
  - **The reward terms remain unknown.** Any topology, hole or vertex-regularity term would partly threaten claim (a) for building polygons.
  - Predecessor VectorLLM: [arXiv 2507.04664](https://arxiv.org/pdf/2507.04664).
- **CurveBench** (Mohseni et al., May 13 2026) — [arXiv 2605.14068](https://arxiv.org/abs/2605.14068) [fetched abstract]; reward description from [search snippet](https://arxiv.org/pdf/2605.14068):
  - Task: recover the containment tree of nested Jordan curves from 756 images.
  - RLVR fine-tunes Qwen3-VL-8B to 33.3% on the Easy tier (Gemini 3.1 Pro: 71.1% Easy / 19.1% Hard).
  - Reward, per the snippet: "parsed into a predicted region tree … compared against the ground-truth containment tree, with the reward combining exact tree-generation correctness and node-count correctness".
  - **Distinction:** the model emits a symbolic tree, so topology is not computed from the model's own geometric output, and there is a single GT. It is still a neighbour for "topology-correctness reward" and should be cited.
- **PolyTopoBench** (Liu, Lao, Sun, Wolff, Xie, Zhao, Jiao, Mai; UT Austin / Amazon / UMD / Emory), NeurIPS 2026 E&D — [lab post](https://sites.utexas.edu/seai/2026/09/29/polytopobench-accepted-at-neurips-2026-evaluations-datasets-track/); [arXiv 2609.32856](https://arxiv.org/abs/2609.32856); [code](https://github.com/seai-lab/PolyTopoBench) [fetched lab post]:
  - About 300K polygons, 10K+ complex polygons, 42K+ interior rings.
  - "Ring-aware metrics that separately score region overlap, exterior boundaries, and interior ring structure."
  - No RL; single OSM-aligned annotation.
  - It could serve as an external eval set for hole handling. Its metrics overlap conceptually with our hole-matching term, but as evaluation, not reward.
- **Topo-R1** — [arXiv 2603.13054](https://arxiv.org/html/2603.13054) [fetched]:
  - Betti-number annotations are used for *data synthesis*.
  - The topology reward is clDice on skeletons of matched anomaly regions, with the LocPen guard.
  - Topology is not computed as β0/β1 of the policy's own output. The prior conclusion stands.
- **Geometric reward credit assignment** — [arXiv 2604.21160](https://arxiv.org/html/2604.21160) [fetched]: routed per-field geometric rewards with a reprojection-consistency verifier. Not topology, not multi-rater.
- **Set-valued reward** — Multi-Answer RL, [arXiv 2603.24844](https://multi-answer-rl.github.io/) [fetched project page]: set Y* of valid answers, text only.
- **Rater-derived band / any-annotator RL:** searches for RL rewards from inter-rater variability in segmentation found only supervised multi-rater methods ([2605.02437](https://arxiv.org/html/2605.02437v1), [2605.08210](https://arxiv.org/pdf/2605.08210)) and medical GRPO with cross-dataset reward normalization ([MedGRPO, CVPR 2026](https://arxiv.org/html/2512.06581v3) [snippet]). None uses a rater band or any-rater topology in the reward.
- **"Learning to Trace" (NeurIPS 2026):** an extended search returned no matching paper. The nearest hits were PolyTopoBench and unrelated "trace" papers ([2605.15672](https://arxiv.org/pdf/2605.15672) "VLMs Trace Without Tracking", a diagnostic). **Still unverified/unreadable.**

### Inferences
- **(a) exact β0/β1 of own output as reward:** not anticipated by anything readable. **Residual risk: VectorLLM++** (big team, polygon RLVR, "geometric quality"). Action: obtain the Research Square PDF manually, since it is too large for fetch. CurveBench should be cited as "topology-correct structured output via RLVR, but symbolic".
- **(b) rater-derived spatially varying band:** no threat found.
- **(c) topology against a rater set:**
  - No threat in geometry/segmentation.
  - Multi-Answer RL is a conceptual precedent for set-valued correctness rewards in text, so phrase our claim as "first set-valued *topological* reward from *expert rater sets* for geometric outputs".
  - Also note that GRPO collapse (Q4) may undercut the distributional story unless addressed.
- **(d) typed-error controller:** components exist (Q3). Claim novelty only for the specific online, typed, exact-verifier-driven combination, and cite 2605.28918, CURE, AVSPO, ODF and PODS.

### Gaps
- VectorLLM++ reward terms and the "Learning to Trace" reward: not readable or not found.
- PolyTopoBench's arXiv paper was not read, so it is unknown whether it includes an RL baseline. The lab post does not mention RL.
- TraversRL (2607.17479), FloorplanVLM (2602.06507) and Groundbench (2609.27821) were not re-checked today; earlier conclusions are assumed to stand.

---

## Q6. Evaluation protocol for RLVR papers in 2026

### Takeaway
Reviewers now expect:
- multiple seeds with confidence intervals;
- pass@k curves vs the base/SFT model out to large k;
- spurious-reward controls (random reward, format-only), because Qwen models gain large amounts from uninformative rewards;
- a cross-family check;
- independent held-out metrics rather than train reward.

For our setting, add a distributional (rater-histogram) metric and a hack audit (NRFR-style).

### Cited Findings
- **Spurious rewards** (arXiv 2506.10947) — [HF knowledge-base mirror](https://huggingface.co/datasets/rl-llm-wiki/knowledge-base/discussions/231); [Notion write-up](https://rethink-rlvr.notion.site/Spurious-Rewards-Rethinking-Training-Signals-in-RLVR-1f4df34dac1880948858f95aeb88872f) [snippet]:
  - Qwen2.5-Math-7B MATH-500 gains: +21.4 pts (random reward), +13.8 (format-only), +24.1 (incorrect labels), vs +29.1 (ground truth).
  - The gains do not transfer to Llama3/OLMo2.
  - Mechanism: GRPO clipping bias amplifies high-prior pretrained behaviours ("code reasoning" 65% → >90%).
  - Recommendation: validate across model families.
- **Large-k behaviour.** Base models overtake GRPO at large k (k = 64: 87.7 vs 86.8) — [arXiv 2508.10751](https://arxiv.org/abs/2508.10751) [snippet]. So report pass@k from k = 1 to 64 vs base and SFT.
- **Measuring diversity directly.** Pass@k alone misses diversity and capability-retention changes — [arXiv 2610.07405](https://arxiv.org/html/2610.07405) [title-level snippet].
- **Hack audit metrics.** RHR, NRFR, reward–oracle gap and "worse rate" separate true improvement from hacking; "higher rewards do not always mean better task performance" — [arXiv 2607.09492](https://arxiv.org/html/2607.09492) [fetched].
- **Preregistration and placebo arms are appearing in RLVR work:**
  - the dark-room paper's 74 preregistered arms, including content-free placebo auxiliary rewards — [arXiv 2607.21273](https://www.alphaxiv.org/overview/2607.21273);
  - "When the Reward Suite Is Leaky: A Preregistered Causal Contrast…" — [arXiv 2607.11022](https://arxiv.org/pdf/2607.11022) [title-level].
- **Seed variance.** "A Sober Look at Progress in LM Reasoning" (Hochlehnert et al.) showed that seed and decoding variance on small benchmarks can exceed claimed RL gains, and recommended multiple seeds with CIs — [arXiv 2504.07086](https://arxiv.org/abs/2504.07086) [prior knowledge, not re-fetched; check exact numbers].
- **RL-only vs SFT+RL at small scale.** Topo-R1 found RL-only collapses (F1 3.1 vs 45.2 with SFT warm-start), so the comparison must include SFT-only, not just base — [arXiv 2603.13054](https://arxiv.org/html/2603.13054) [fetched].

### Inferences
- **Minimum protocol:**
  - 3–5 seeds per arm, with mean ± 95% CI;
  - per-image paired bootstrap, or a Wilcoxon signed-rank test with Holm correction across metrics, comparing RL vs SFT-only;
  - pass@k (k = 1…64) on any-rater topology match and on IoU ≥ τ;
  - controls: random reward, format/legality-only reward and IoU-only reward, to isolate the topology/band terms;
  - one non-Qwen ≤2B backbone (e.g. InternVL-2B, which Topo-R1 used) to rule out the Qwen prior effect;
  - a rater-distribution metric (per-image TV/KL between the sampled-topology histogram and the rater histogram);
  - a hack audit (area inflation, part-count inflation, NRFR) on a held-out set.
- **Freeze the protocol early.** Preregistering arms and stopping rules is cheap and is becoming a reviewer expectation.

### Gaps
- I found no community-standard statistical test for RLVR papers; practice varies.

---

## Candidate hypotheses to test locally

All tests below are cheap and offline: they run on Phase 0 scored outputs or synthetic data, with no GPU training.

1. **H1 – Per-term normalization ranks errors correctly; summed GRPO does not.**
   - Setup: take Phase 0 outputs grouped per image (or synthesize groups of G = 8 by perturbing them). Compute advantages under:
     - (i) GRPO on the summed reward;
     - (ii) Dr. GRPO (mean-only, no σ);
     - (iii) GDPO (per-term z-score → sum → batch z);
     - (iv) GDPO with a per-term σ floor;
     - (v) GDPO with gating (topology/hole/band only when IoU > τ).
   - Measure: the Kendall τ between advantage order and an expert "badness" ranking (or the lexicographic order legality > topology > IoU); the share of advantage mass from near-constant terms.
   - Prediction: (iii) assigns large advantages to tiny legality differences (dark-room effect); (iv)/(v) fix it.
2. **H2 – Dark-room amplification exists in our reward vector.**
   - Measure, per term, the distribution of within-group σ and the max |advantage| contributed when σ < 0.01.
   - Prediction: legality and hole terms have σ ≈ 0 on >50% of groups after SFT. Pick σ_min from this distribution.
3. **H3 – A hack-seeking policy can beat the reward.**
   - Construct adversarial outputs: a convex hull of the rater union, a dilated union, many small parts, a polygon with all rater holes stamped in.
   - Prediction: without an area/part-inflation multiplier (an analogue of LocPen/odLength), at least one adversary scores above the median honest output. Choose the guard's threshold by sweeping on these.
4. **H4 – GRPO collapses rater-topology diversity.**
   - Toy model: a categorical policy over K topologies, with rater frequencies, e.g. [0.6, 0.3, 0.1], and success probabilities per topology.
   - Run exact GRPO updates with (a) the any-rater max reward; (b) the any-rater reward plus a group coverage bonus; (c) a FlowRL-style frequency-weighted target.
   - Prediction: (a) converges to a point mass; (b)/(c) keep entropy near the rater histogram. This tells us whether a distributional term is needed before spending GPU time.
5. **H5 – A typed controller beats scalar difficulty filtering on synthetic curves.**
   - Simulate per-item, per-error-type learning curves (logistic, with different rates; one type unlearnable or noisy, one prone to overfitting).
   - Compare item sampling by:
     - uniform;
     - scalar p(1−p) / variance (ODF/PODS-style);
     - per-type variance-weighted (Level 2).
   - Measure: steps to a target held-out error profile; false-alarm rate of the overfitting alarm (train-type-error ↓ while val-type-error flat for N evals); early-stop regret.
   - Prediction: typed sampling wins only when error types have different learnability. If not, Level 2 reduces to ODF.
6. **H6 – ACR/near-zero-variance rate after SFT.**
   - On Phase 0 outputs, compute the share of groups whose summed-reward σ is below thresholds (1e-3, 1e-2).
   - If it is high, plan for std-threshold dynamic sampling or AVSPO-style virtual samples.
7. **H7 – Length bias.**
   - Correlate vertex count with reward and with GRPO's per-sequence-normalized advantage on Phase 0 outputs.
   - Prediction: per-sequence normalization favours short polygons, which justifies the Dr. GRPO constant denominator.

## Recommendation

- **Algorithm (default):**
  - GRPO-family loss with token-level aggregation and a Dr. GRPO constant denominator;
  - DAPO clip-higher (ε_low 0.2, ε_high 0.28);
  - G = 8; small fixed KL (start ~1e-3, allowed range 0 to 1e-2) with an optional ProRL-style reference reset;
  - SFT warm-start (required: Topo-R1 RL-only collapsed);
  - std-threshold dynamic sampling (drop or down-sample groups with σ below ε, PODS-style).
- **Ablation arm:** GSPO (as in SimpleSeg). Skip CISPO and GMPO unless instability appears. Avoid RLOO, which showed the worst hacking rate in 2607.09492.
- **Reward combination:**
  - Use GDPO (per-term group z-score → weighted sum → batch z), but with:
    - a per-term σ floor, or mean-only centring, for low-variance terms (legality, holes, topology);
    - hard gates: legality as a binary gate (as Topo-R1 does for format); topology, hole and band terms only when overlap > τ (as GDPO's conditioning does);
    - a multiplicative area/part-inflation guard modelled on LocPen/odLength.
  - Keep IoU-like terms continuous, with smooth tiers like Topo-R1's φ.
  - Treat field-routed advantages (GRCA) for per-ring terms as an optional later ablation.
- **Distribution:**
  - Do not assume GRPO preserves rater topology diversity. Add a group-level coverage term over matched rater topologies (Multi-Answer-RL style) or a frequency-weighted target (FlowRL style), gated on H4.
  - Report a rater-histogram divergence metric.
- **Level 2 controller:**
  - Build on per-item, per-error-type within-group reward variance (typed p(1−p)).
  - Actions: sampling weights, per-term gate thresholds, and KL within [0, 1e-2] plus a reference reset.
  - Alarms:
    - ACR;
    - per-term σ collapse (dark room);
    - NRFR-style "newly rewarded failures" on a held-out audit set;
    - adversarial-shape score drift;
    - entropy / unique-topology decline (diversity-collapse overtraining);
    - train–val divergence per error type for early stopping.
  - Frame novelty narrowly (online, typed, exact-verifier-driven) and cite 2605.28918, CURE (2601.15408), AVSPO (2605.21125), ODF (2504.03380) and PODS (2504.13818) as neighbours.
- **Novelty:**
  - Claims (a), (b) and (c) survive today's refresh. Add CurveBench, PolyTopoBench, Multi-Answer RL and GRCA as neighbours.
  - The main outstanding risk is VectorLLM++'s unreadable RL reward: get the PDF manually before submission.
  - Keep retrying "Learning to Trace" via the NeurIPS 2026 proceedings and OpenReview.
- **Evaluation:**
  - 3–5 seeds with CIs; paired per-image tests vs SFT-only;
  - pass@k to 64 vs base/SFT;
  - random-reward, legality-only and IoU-only controls;
  - one non-Qwen backbone;
  - a hack audit;
  - preregister the arms and stopping rules.
