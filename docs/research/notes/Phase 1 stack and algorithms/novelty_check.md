# Novelty check: topology rewards, rater-distribution matching, egg-yolk band (2026-10-19)

How to read the flags:
- **[V]**: verified from the primary text (PDF or HTML read directly).
- **[P]**: verified only from an abstract, a venue page or a secondary summary.
- **[U]**: unverified.

The search was targeted, not exhaustive. Contributions (1), (3) and (4) are "not found", which is weaker than "does not exist".

---

## 1. VectorLLM++ (closest geometry RLVR prior) [V]

- **Citation:** Shunping Ji (Wuhan Univ.), *VectorLLM++: A Unified Next-Token-Prediction MLLM for Dense Remote Sensing Perception and Vector Mapping*, Research Square preprint, posted 2026-09-10.
  - DOI: https://doi.org/10.21203/rs.3.rs-10949232/v1
  - PDF: https://assets-eu.researchsquare.com/files/rs-10949232/v1_covered_cdfb9971-0996-41cf-89ad-7188625d582f.pdf (29 pp.)
- **RL setup (Sec. 4.4):** GRPO with K rollouts and group mean/std advantages. The KL term sits in the actor loss. Rollouts are pre-filtered (hard-sample mining) to drop prompts whose rewards are all 0 or all 1.
- **Reward terms (exact):**
  - **Format:** if JSON parsing fails but a fallback parser recovers the output, the reward is multiplied by 0.8. Unrecoverable output scores 0.
  - **Box (Eq. 13):** bidirectional best-match IoU per phrase, R_box(p) = ½[(1/|B_p|) Σ_b max_{b*} IoU + (1/|B*_p|) Σ_{b*} max_b IoU], then averaged over phrases.
  - **Polygon (Eq. 14):** R_poly = 0.75·R_IoU + 0.25·R_node.
    - R_IoU is the same bidirectional score, computed on rasterized polygon masks.
    - R_node is a vertex-efficiency score on Hungarian-matched pairs (matched on IoU). It penalizes having more vertices than the matched reference.
  - **Mask tokens:** prefix-consistency ratio, aggregated bidirectionally.
  - **VQA:** exact match.
- **What it does not do:**
  - It never counts components, holes or Betti numbers.
  - It uses one reference annotation y\*. Multi-annotator consensus is only mentioned as future work, for ambiguous *classes*.
  - It has no tolerance band. ("Tolerance" appears only in data-engine polygon simplification.)
  - "Topology consistency" appears only in the discussion, as a desideratum.
  - No hole or interior-ring representation is mentioned.

## 2. "Learning to Trace" [P, reward U]

- **Citation:** Mahsa Geshvadi, Funda Durupinar, *Learning to Trace: Zero-Shot Connectivity Repair for Vessel Segmentation via Procedural Reinforcement Learning*. This is a NeurIPS 2026 **workshop** paper (RL for Experimental Sciences), not main track.
  - https://neurips.cc/virtual/2026/161763
  - https://openreview.net/forum?id=NZuAkfqSsw (the PDF was blocked by an OpenReview challenge).
- **Method:** a memory-augmented actor-critic that traces tubular structures step by step, post-processing a segmenter's output. It is trained on synthetic data with an adaptive curriculum. It is not a VLM and not GRPO.
- **Results:** connected components 29.1 → 4.5 on DRIVE segmenter outputs, and centerline recall 85.8 → 87.7%. The component count is reported as an **evaluation metric**.
- **Reward:** not visible in the abstract. **Unverified** whether the reward counts components. The likely reward is per-step on-centerline/coverage against synthetic ground truth, but this is a guess.
- **Relevance:** it supports the motivation that "fragments matter even when overlap is high". It is not a per-instance topology reward on a generative model's own polygons.

## 3. Distribution-matching prior art, and where r_topo sits

### MaD-RL: the decisive prior [V]

- **Citation:** Kulkarni, Vepuri, Demir, …, Nagpal (Meta Superintelligence Labs), *MaD-RL: Matching Distributions for Calibrating LLMs with RL*, 2026-09-23.
  - Meta page: https://ai.meta.com/research/publications/mad-rl-matching-distributions-for-calibrating-llms-with-reinforcement-learning/
  - alphaXiv: https://www.alphaxiv.org/abs/2609.matching-distributions-calibrating-llms-rl
  - arXiv id not resolved; the arXiv API returned 503.
- **Setup:** GRPO with a category z = f(y), and a target P\* that "may depend on the prompt" (Sec. 2.1).
- **Reward:** derived via the score-function identity, R_D(z) = −∂D/∂P_π(z) + baseline.
  - **Reverse KL (Eq. 10):** R_RKL(z) = log P\*(z) − log P̂_π(z).
  - **Plug-in estimator (Eq. 6):** P̂_π(c) = (1/G) Σ_j 1{z_j = c}. It **includes self**, and the paper notes it is biased for KL and JS.
- **Gated combination (Eq. 17):** R_i = R_i^correct · (1 + λ_g R_i^div).
- **Results:** reverse KL fits peaked targets best. All variants cost 1.2–8.5 pass@1 points relative to correctness-only RLVR. A collapsed category gets no signal.
- **Implication:** our r_topo = 1[valid]·(1 + λ(log f(a) − log π̂_loo(a))) **is MaD-RL's gated reverse-KL reward**, with three substitutions:
  - **Target:** P\* := f, the per-image empirical rater frequency over a computed topology class.
  - **Estimator:** plug-in → leave-one-out.
  - **Gate:** correctness → set-valued validity, which here equals supp(f).

### Other distribution-matching work

- **GAPO** (Anschel et al., EMNLP 2025) [P]: within-group frequency reward 1/K − P̂(z) toward a uniform target. MaD-RL shows this is the L2 special case.
- **Distributional Alignment Games**:
  - Mohri, Schneider, Wu, arXiv [2604.27166](https://arxiv.org/abs/2604.27166), ICML 2026 [P].
  - Mohri, Schneider, Zhong, *Generalized DAG for Unbiased Answer-Level Fine-Tuning*, arXiv [2605.02435](https://arxiv.org/abs/2605.02435) [V].
  - What they do:
    - The answer-level marginal objective becomes a game. The per-sample reward is β·log q\*(z), where q\* is computed from group counts (consensus or inverse-frequency targets).
    - 2605.02435 shows that the plug-in log reward has O(1/K) Jensen bias (Eq. 7). It gives U-statistic unbiased estimators for polynomial rewards, and minimax-polynomial or Taylor corrections (O(1/K²)) for log.
  - **Implication:** the bias of finite-group log-frequency estimators, and "exclude-self"-style corrections, are already studied. LOO is still Jensen-biased: E[log π̂_loo] ≠ log π.
- **FlowRL** (Zhu et al., arXiv [2509.15207](https://arxiv.org/abs/2509.15207), ICLR 2026) [V]:
  - Target π̃(y) ∝ π_ref(y)·exp(β r̂(y)), over **whole trajectories**.
  - Trained with a trajectory-balance loss (learned log Z, 1/|y| length normalization, clipped importance weights).
  - **r_topo is not a special case.** Setting r = β⁻¹ log f(a) gives an a-marginal ∝ π_ref(a)·f(a), not f(a). FlowRL has no external target on a marginal.
  - FlowRL is the right citation only for the general idea that "reverse KL to a target beats reward maximization for diversity".
- **Outcome-based Exploration** (Song, Kempe, Munos, arXiv [2509.06941](https://arxiv.org/abs/2509.06941), NeurIPS 2025) [P]: within-batch repetition penalty on final answers, plus UCB bonuses on rare answers. This is the uniform/diversity version of group-count rewards.
- **Multi-Answer RL / "Reaching Beyond the Mode"** (Puri, Damani, Shenfeld, Ghassemi, Andreas, Kim, arXiv [2603.24844](https://arxiv.org/abs/2603.24844)) [P]:
  - The model emits a *set* of answers with confidences in one pass. The reward is set-level coverage and calibration.
  - There is no group-frequency target matching.
- **SHALA-LLM** (Wu, Wang, Ong, Liang, Picard, arXiv [2606.05376](https://arxiv.org/abs/2606.05376)) [V via HTML summary]:
  - GRPO for classification with annotator distributions (ChaosNLI, GoEmotions).
  - Reward = H̃(p_q)·[1 − D_JS(p̂, p_q)], where p̂ is the distribution the model *verbalizes* in one rollout, not group frequency.
  - This is the closest "match human rater distribution in RL" prior.
- **Fitness sharing** (Goldberg & Richardson 1987; niching) [from known literature, not re-fetched]:
  - Shared fitness f_i/m_i. Under proportional selection, the equilibrium population share is ∝ niche resource.
  - r_topo is the log-space analogue: log f − log m̂, where m̂ is the group niche share. This is a useful framing, but not a novelty.

### What is genuinely new in (2)

Only the instantiation:
- The target is the per-instance empirical distribution of *expert raters* over a *computed geometric invariant* (the topology of the parsed polygon), with no category label in the model's output.
- The validity gate equals supp(f). Because invalid outputs score 0, reverse KL is never evaluated where f = 0, so MaD-RL's ε-smoothing is unnecessary.
- The LOO estimator is a minor engineering choice. It removes the P̂ ≥ 1/G self-inclusion floor but makes singletons give log 0, so it needs add-α smoothing.

### Caveat to check: the gate can reward invalid outputs (our own derivation, not from the literature)

- With stop-gradient frequencies, the gated reward's expected update is (1 − λ·KL(π_c ‖ f))·∇V + (terms that move π_c toward f).
  - V is the valid mass.
  - π_c is the topology distribution conditional on being valid.
- So when λ·KL > 1, invalid outputs are *favored*. Equivalently, any valid rollout with λ(log f − log π̂) < −1 scores below an invalid one.
- Fixes:
  - clip the bracket at ≥ ε;
  - use λ ≤ 1/max|log f − log π̂|; or
  - give invalid outputs a reward below the valid floor.
- MaD-RL calls its gated form "a heuristic" and shares this property.

## 4. Topology, multi-rater and tolerance rewards in 2026

- **Topo-R1** (Xu et al., arXiv [2603.13054](https://arxiv.org/abs/2603.13054), Mar 2026) [V via HTML summary]:
  - A VLM trained with GRPO to **detect** topological errors in a *given* mask, using Betti-verified synthetic error injection.
  - R = 0.10·fmt + 0.85·acc + 0.05·(1 − clDice)·LocPen.
  - It does not reward the topology of the model's own segmentation, and uses no multiple raters.
- **Topology-aware RL in VLM segmentation:** Seg-Zero, SAM-R1, LENS, GenSeg-R1 and Dr. Seg use IoU, L1 and point terms, with no topology. I found no paper rewarding exact β₀/β₁ (part/hole counts) of a model's own output.
- **Multi-rater segmentation:** prior work uses distribution matching via generative models and GED/Hungarian-matched IoU (Prob. U-Net / PHiSeg line, LIDC), or calibration (*Multi-Rater Calibrated Segmentation Models*, arXiv [2605.02437](https://arxiv.org/abs/2605.02437)). None is RL, none conditions IoU on a topology class, and none is VLM.
- **"Disconnected by Design"** (Lacuna/Tiptree hypothesis page, https://lacuna.tiptreesystems.com/hypothesis/disconnected-by-design-matching-the-topology-of-expert-disagreement-in-bd0623e939bf867a) [V, not a paper]:
  - A public **idea** to cluster expert masks into topological modes and set-match samples to experts with a Betti-matching + Dice transport cost.
  - No RL, no VLM, no authors. It is not citable prior art, but the "topology of expert disagreement" idea is publicly floated.
- **Tolerance from inter-rater variability:** Nikolov et al., surface Dice, arXiv [1809.04430](https://arxiv.org/abs/1809.04430) [P].
  - The per-organ tolerance τ is the 95th percentile of inter-observer surface distances.
  - This is the direct prior for "tolerance from the raters' own variability". The egg-yolk idea traces to Cohn & Gotts (1996), regions with indeterminate boundaries.
  - I found nothing using a band between rater intersection and union, normalized by contour length, as an **RL reward**.

## 5. Recommended claim wording

| # | Verdict | Wording |
|---|---|---|
| 1 | **Safe (narrowed)** | "To our knowledge, the first RLVR reward that computes the exact per-instance topology (number of parts, holes per part) of a VLM's own polygon output, and accepts it if it matches *any* expert rater's topology. Prior geometric RLVR (VectorLLM++) rewards IoU, vertex count and format against a single reference; Topo-R1 detects topological errors in given masks; Learning to Trace repairs connectivity with a non-VLM tracer." |
| 2 | **Narrow hard; do not claim the objective** | "We instantiate the gated reverse-KL distribution-matching reward of MaD-RL (Kulkarni et al., 2026; cf. GAPO, Distributional Alignment Games) with the *per-image empirical rater distribution over a computed topological invariant* as target. The validity gate coincides with the target's support, so no smoothing is needed; we use a leave-one-out group estimate." Do **not** claim novelty for the log-ratio form, for group-count estimation, or for LOO. Cite Mohri 2026b for the estimator bias. |
| 3 | **Narrow / demote to design choice** | "IoU is computed against the subset of raters sharing the sample's topology." No prior found, but it is a small variant of matched-rater / GED-style scoring. Present it as a method component, not a headline contribution. |
| 4 | **Narrow** | "A rater-calibrated egg–yolk band reward: boundary deviation outside the band between rater intersection and union, normalized by contour length, with tolerance set from raters' leave-one-out agreement (cf. surface Dice tolerance from inter-observer variation, Nikolov et al.), used as an RL reward." Claim novelty only for *its use as an RL reward with rater-count correction*, not for inter-rater-derived tolerance. |
| 5 | **Drop as a scientific claim [U, not searched]** | Present typed-verifier logging and rollback as infrastructure or engineering ("training controller"). It is unlikely to survive review as a contribution. |

**Suggested headline:** a topology-exact, multi-rater-valid verifiable reward for VLM polygon delineation, plus an application of distribution matching (MaD-RL style) to rater topology frequencies.
