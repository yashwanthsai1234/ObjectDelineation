# H2: The RL objective and estimator, on mathematical ground

Scope: read-only study of `src/objdel/reward/{score,spec,advantages,cache}.py`, `src/objdel/train/{rewards,configs,fast}.py`,
`experiments/phase2/{train_rl,g_probe}.py`, `experiments/phase2/configs/{rl_full,rl_shifted}.yaml`, installed TRL
(`.venv/.../trl/trainer/grpo_trainer.py`), and the repo notes `docs/research/notes/{full_rl,rl_plateau_diagnosis}_2026-11-08.md`.
No edits, no git, no spend. Session transcripts, memory and `.claude/` were not read.

## Skills used

| Skill | What it contributed |
|---|---|
| nvidia-skills:nvidia-skill-finder | Checked the live NVIDIA catalog (`skills.sh.json`: Training AI, Inference AI, GPU Development groups). No strong match for RL objective/estimator math. `cuequivariance-*` handles E(3)-equivariant networks, not 2-D translation in a VLM. `nemo-rl-*` is a different RL framework. Nothing to install. |
| engineering:code-review | Correctness review of `advantages.py`/`score.py`. It found the per-term std normalization (a biased, step-varying reweighting, §2.2), the coupled leave-one-out terms that bias the baseline (§1.4, §2.2), and an exploitable margin in the topology fallback (§3.4 H1). Section 0 corrects two premises of the brief. |
| engineering:debug | Repro → isolate → diagnose → fix structure for the plateau (§2.5-2.6, §4). The "collapse" is the predicted optimum of a linear-in-policy objective. Mismatched inputs and a k-mismatched ceiling explain the gap. Each claim is paired with a measurement (§6). |
| engineering:system-design | Requirements first: the meaning of "correct" (one sample vs calibrated diversity) picks the objective (§5). Explicit trade-offs and assumptions, plus what to revisit. |
| engineering:architecture | ADR-style choice in §5: max-IoU vs GED vs Hungarian, and absolute vs anchor-relative coordinates, with the consequences of each. |
| engineering:documentation | Layout: verified facts first, then theorems with stated assumptions, then a decision and a measurement list. |
| huggingface-skills:huggingface-best | Queried HF official benchmark datasets (58). None covers multi-rater segmentation (no LIDC/GED leaderboard), so the comparison metrics (GED, HM-IoU) come from the literature, not a leaderboard. |
| huggingface-skills:huggingface-papers | Pulled the abstracts/text of 2505.22617 (entropy mechanism: ΔH ≈ −Cov(log π, Δlogit)), 2503.06639 (GRPO effective loss under std normalization) and 2503.20783 (Dr. GRPO). |
| huggingface-skills:trl-training | Checked TRL semantics in the installed source. `dr_grpo` divides by `B·max_completion_length`. `scale_rewards="none"` skips std scaling. `mask_truncated_completions` zeroes the loss mask but keeps the rewards in group statistics. Clipping uses ε_low/ε_high. β·KL is added per token. |

Also used: Exa web search (positive definiteness of the Jaccard matrix: Gower 1971, Bouchard et al. 2013, arXiv 2608.13549).

---

## 0. Facts from the code (verified) and corrections to the brief

F1. **IoU_tc excludes no pixels.** `score_numpy.score_mask` / the CUDA scorer compute
`IoU_tc(x) = max_{j ∈ S(sig(x))} IoU(x, Y_j)`, where `S(s)` is the set of voting raters whose size-filtered Betti
signature equals the prediction's. If no rater drew that signature, it falls back to `IoU(x, vote)`. Plain IoU, no
exclusion. The disagreement band is excluded in **F5**: `D = Σ_{x∖U}(d−½) + Σ_{I∖x}(d−½)`. Error pixels inside the
egg-yolk band `U∖I` (union minus core) cost nothing.

F2. **There is no scalar reward.** The advantage is (`advantages_v4`)

    A_i = Σ_k w_k (X_ik − X̄_k) / max(σ̂_k, f_k)  +  c (r_topo,i − r̄_topo),     X = terms × 1[legal]
    w = (0.3, 0.3, 0.1, 0.3) for (IoU_tc, F5, legal, g_area),  f = (0.05, 0.05, 0.25, 0.05),  c = 2

Groups whose terms all span < 0.01 and that share one topology answer get A = 0. TRL adds nothing (sum_then_normalize,
scale_rewards none, aggregator weight 1). So TRL's `dr_grpo` "no std scaling" holds at the TRL level, **but the
advantage is std-normalized per term, per group.**

F3. Loss: `dr_grpo`, i.e. `Σ_tokens / (B · L_max)`, L_max = 1536. One generation per optimizer step, one iteration, so
the PPO ratio is 1 up to the vLLM/trainer mismatch. The diagnosis note confirms "PPO clip never active". **ε_high = 0.28
(clip-higher) therefore does nothing.** β = 1e-3 (k3 KL to the SFT reference), LR 5e-5 constant, LoRA r16.

F4. `mask_truncated_completions=True` (TRL source, l.2604): truncated rows lose their loss mask. Their rewards (all terms
0) stay in the group mean.

F5. SFT trains on one row per voting rater, which is MLE of the **rater mixture**. The SFT policy is therefore a
calibrated sampler (forward-KL fit to Q). RL then moves it to a different objective (§1).

F6. The "doctor ceiling" (0.732) is the leave-one-out max over **k−1** doctors. The model's IoU_tc is a max over **k**
doctors. The two numbers are different order statistics (§1.3).

---

## 1. The objective

**Setup.** Image I. Raters Y_1..Y_k ~ Q_I i.i.d. (A1: exchangeable raters; k ∈ {2,3,4} on RL items). The policy
π(·|I) is a law on the finite mask space 𝓜 (the push-forward of text through parse+raster). Per-sample utility
u(x; Y_{1:k}). ū_I(x) := E[u(x; Y_{1:k}) | I].

### 1.1 Theorem 1 (linear objectives elicit a point, not a distribution)

For any per-sample reward, J(π) = E_{x~π} ū(x) is linear in π, so its maximizers are exactly the laws supported on
argmax ū. As a scoring rule for Q, S(π, y) = E_π u(x, y) is not strictly proper. It is not even proper unless Q itself is
supported on argmax ū.

*Proof.* E_{Y~Q} S(P, Y) = Σ_x p(x) ū(x) ≤ max ū, with equality iff supp P ⊆ argmax ū. Take P = δ_{x*}. Then Q is
optimal only if Q(argmax ū) = 1. ∎

With the KL term the optimum is the Gibbs policy (standard; e.g. Rafailov et al. 2023, Eq. 4):

    π*_β(x|I) ∝ π_ref(x|I) · exp(ū(x)/β_eff).

Late in training the effective reward is about 6·IoU (§2.2), so β = 1e-3 corresponds to an IoU resolution of about
1.7e-4. A gap of 0.01 IoU changes the odds by e^{60}. **The entropy collapse (0.27 → 0.06) is the predicted optimum of
the chosen objective, not a malfunction.** Any max-IoU, mean-IoU, Dice or F5 reward behaves this way. Diversity can only
come from a reward that is non-linear in π (§1.4) or from a large β.

### 1.2 What max-over-raters elicits

Let ū_k(x) = E max_{j≤k} IoU(x, Y_j). Then:

(a) ū_k ≥ ū_1 = E_Q IoU(x, Y), and ū_k is nondecreasing in k.

(b) Atoms. On finite 𝓜, ū_k(y) ≥ 1 − (1 − q(y))^k for y ∈ supp Q. As k → ∞, ū_k(x) → max_{y∈supp Q} IoU(x, y), and
**every** rater shape becomes a maximizer. That is exactly "any expert's shape is correct". At finite k the Bayes act is
mode-seeking: it prefers shapes that are both frequent (q(y)) and central among the other raters' shapes.

(c) It is not k-invariant. Items with more raters get larger, flatter ū. Group centering removes the per-item offset but
not the change in curvature.

(d) Topology coupling (a hidden bias). Within signature class s, IoU_tc is a max over N_s ~ Bin(k, f_s) references.
E max over N references is increasing in N. So the IoU term alone rewards picking the **majority topology** even when a
minority signature is equally well drawn. The reverse-KL topology term (§3.4 H7) is what counteracts this, so it is
load-bearing, not cosmetic.

Bottom line: E max-IoU_tc is Fisher-consistent for the metric the repo reports (one sample, best rater). It is the right
objective if the product is one outline judged by "matches some doctor". It is the wrong objective for calibrated
uncertainty.

### 1.3 The ceiling comparison has an order-statistic bias (in the model's favour)

E max_{j≤k} IoU(x, Y_j) ≥ E max over a random (k−1)-subset. A like-for-like model score is the U-statistic

    M_{k−1}(x) = (1/k) Σ_j max_{i≠j} IoU(x, Y_i)      (unbiased for the random-(k−1)-subset max)

and that is what should be compared with the doctors' leave-one-out ceiling. The true gap is **larger** than
0.732 − 0.645. The same applies to "centred test items 0.77 vs ceiling 0.73": that is not yet evidence of super-doctor
performance.

### 1.4 Alternatives

**(i) Expected IoU under a random rater, ū_1.** This is the consensus semantics. The Bayes act is a Jaccard 1-median of
Q (cf. Nowozin 2014, CVPR, for the pixel-independent case: a thresholded marginal). It penalizes a valid minority
outline. Still linear, so still a Dirac.

**(ii) GED / energy score.** Distance d = 1 − IoU (Kohl et al. 2018, Probabilistic U-Net):

    D²_GED(P, Q) = 2E d(X,Y) − E d(X,X') − E d(Y,Y').

*Theorem 2. With d = 1 − IoU, the energy score ES(P, y) = E_P d(X,y) − ½E_P d(X,X') is strictly proper on 𝓜.*

*Proof.* Let K be the IoU Gram matrix on 𝓜 (convention IoU(∅,∅)=1, IoU(∅,A)=0, as in `score_numpy.iou`). Then
D²_GED = 2(1 − pᵀKq) − (1 − pᵀKp) − (1 − qᵀKq) = (p−q)ᵀK(p−q). Bouchard, Jousselme and Doré (2013, IJAR) proved that K
restricted to the non-empty subsets is **strictly** positive definite (finite proof: arXiv 2608.13549, Lemma 4.1; PSD
since Gower 1971). The empty-set row is e_∅, so K stays strictly PD. Hence D²_GED = 0 iff p = q. Also
E_Q ES(P,Y) = ½D²_GED(P,Q) + const. ∎ (Also: 1 − IoU is a metric of negative type. 1 − Dice is not a metric.)

*Policy-gradient form (new derivation, needed for correctness).* Maximize J(π) = −E_π d̄(X) + ½E_{π⊗π} d(X,X') with
d̄(x) = (1/k)Σ_j d(x,Y_j). Because ∇E_{π⊗π} d = 2E_π[∇log π(X) · E_{X'} d(X,X')],

    ∇J = E_π[ ∇log π(X) · ( −d̄(X) + E_{X'~π} d(X, X') ) ].

So the per-sample reward is

    r_i = (1/k) Σ_j IoU(x_i, Y_j)  −  (1/(G−1)) Σ_{l≠i} IoU(x_i, x_l)        (+ const)

**with coefficient 1 on the diversity term, not ½.** With ½ the objective becomes E ū_1 − ¼E IoU(X,X'), which is
improper and under-diversifies.

*Baseline caveat.* r_l (l ≠ i) depends on x_i. The usual LOO baseline b_i = mean_{l≠i} r_l then changes the
diversity coefficient to 1 − 1/(G−1) = 6/7 at G=8, which is improper again. The correct baseline excludes every term
containing x_i:

    b_i = mean_{l≠i} ū_1(x_l) − [1/((G−1)(G−2))] Σ_{l≠m; l,m≠i} IoU(x_l, x_m).

The same O(1/G) coupling exists today in the `r_topo` term (its π_loo uses the other rollouts) under plain group-mean
centering.

**(iii) Hungarian set reward** (M hypotheses in one completion, matched to the k raters):
R = max_σ (1/M) Σ_m IoU(h_m, Y_σ(m)) = 1 − W₁ between empirical measures under cost d. Sample Wasserstein losses have
biased gradients and improper finite-sample minimizers (Bellemare et al. 2017, "The Cramér distance as a solution to
biased Wasserstein gradients"). Worked example: Q = 0.8δ_a + 0.2δ_b, d(a,b)=1, M=k=2. The set {a,a} has expected cost
0.20, the set {a,b} has 0.34, so the reward prefers the collapsed set. In general it computes a quantization of Q, not Q.
It also costs M× output tokens and a format change. Use HM-IoU as an **evaluation** metric (as in PHiSeg / stochastic
segmentation network papers), not as a reward.

**(iv) Soft-IoU / Lovász extension.** These are convex surrogates for differentiable per-pixel predictors. Berman et al.
2018 show the Lovász extension is the tight convex closure of the submodular Jaccard loss. Here the policy is a sampled
discrete text, and REINFORCE needs only a bounded measurable reward. Estimator variance depends on the law of ū under
π, not on its smoothness in x. **Not applicable, no benefit.**

**(v) Dice vs IoU.** Dice = 2J/(1+J) is strictly increasing in J, so rankings and argmax are identical (it commutes with
max). Under centering, advantages scale by the slope 2/(1+J)²: 0.64 at J=0.77 and 1.65 at J=0.1. Dice therefore
weights misses about 2.6× more relative to saturated items. Its concavity in J makes it risk-averse under uncertainty.
1 − Dice breaks the triangle inequality, so the properness of GED (Thm 2) is lost. Keep IoU. Handle misses with an
explicit localization term (§5).

**Zero-overlap flat region.** For misses IoU ≡ 0, so the IoU term cannot rank them. Only F5 orders them (spill depth
D ∝ A·dist, hence F5 ∝ 1/dist). That signal is weak and is crushed by the floor and z-scoring. 23% of offset > 0.3 test
items are misses (diagnosis note), so this is the term that matters for the held-out gap.

**Which objective is right?**
- "Any expert's shape is correct, emit one outline": E max-IoU_tc (keep), plus a localization term for misses.
- "Calibrated diversity / uncertainty": the GED energy reward of (ii), the only strictly proper option among those
  listed. Its optimum is π = Q, so its single-sample max-IoU equals the doctors' own leave-one-out level. The max-IoU
  Dirac optimum is at least that, so **one policy cannot be optimal for both**. The trade-off is intrinsic.

---

## 2. Gradient estimator

Notation: s_i = ∇log π(o_i) (sum of token scores). TRL's update is ĝ = (1/(B·L_max)) Σ_q Σ_i A_qi s_qi (+ β∇KL).
With ratio ≡ 1 (F3) this is exactly REINFORCE with baseline.

### 2.1 Mean centering is RLOO (unbiased)

r_i − r̄ = ((G−1)/G)(r_i − r̄_{−i}), and r̄_{−i} is independent of o_i. So E[(r_i − r̄)s_i] = ((G−1)/G)∇J_q: unbiased
up to a constant (Kool et al. 2019; Ahmadian et al. 2024).

### 2.2 Per-term std normalization (what the code does) is a step-varying reweighting

Write σ̃_qk = max(σ̂_qk, f_k). To first order, E ĝ ∝ Σ_q Σ_k (w_k/σ̃_qk) ∇J_qk. That is a different objective every
step, with prompts weighted by 1/σ (the Dr. GRPO argument, 2503.20783). σ̂ also depends on o_i, which adds a
finite-sample bias. For a binary term, ∇p/√(p(1−p)) = 2∇arcsin√p: std scaling ascends a variance-stabilized objective
(cf. Mroueh 2025, 2503.06639). For a continuous term it has no fixed target.

**Inference from the logs.** If every term had σ̂ ≥ floor, each z would have unit in-group sd. Σ w_k z_k would then
have sd between √Σw² = 0.53 (uncorrelated) and Σw = 1.0 (perfectly correlated), before the topology term. The observed
late advantage sd is 0.34 < 0.53. So **most late groups are floored**, and the late regime is a fixed, unnormalized
scalarization:

    A ≈ centred[ 6·IoU_tc + 6·F5 + 6·g_area + 0.4·legal ] + 2·centred[r_topo].

Early it was σ-normalized. The objective therefore changed implicitly during the run, and F5 and g_area (already near
saturation) carry the same weight as IoU. Verify with M2 (§6).

### 2.3 Optimal baselines and James–Stein: at most 12.5% to gain

The variance-optimal constant baseline is b* = E[r‖s‖²]/E[‖s‖²] (Weaver & Tao 2001; Greensmith, Bartlett & Baxter
2004). With LOO, Var(r_i − r̄_{−i}) = σ²(1 + 1/(G−1)) (assuming r is independent of ‖s‖). So even an oracle baseline
removes at most σ²/(G−1) out of σ²G/(G−1), i.e. 1/G = **12.5% at G = 8**.

James–Stein shrinkage toward a predicted per-item mean m̂_q, b = m̂ + λ(r̄_{−i} − m̂) with
λ* = τ²/(τ² + σ²/(G−1)), recovers only part of that 12.5%. Not worth the complexity.

A cheap principled tweak: if token scores are roughly orthogonal, ‖s_i‖² ∝ L_i (assumption), and
b*_i ≈ Σ_{l≠i} r_l L_l / Σ_{l≠i} L_l (a length-weighted LOO mean). The effect is likewise ≤ 12.5%.

### 2.4 Length bias

`dr_grpo` uses a constant normaliser L_max, so each sequence contributes A_i s_i/L_max: exact REINFORCE with no length
bias. GRPO's 1/|o_i| under-penalizes long negative-advantage outputs per token. Remaining effects:
- Variance grows with L (more score terms). Polygons with more vertices are noisier but unbiased.
- Truncation masking (F4) lowers the group mean by r̄/G per truncated row. That adds a constant c to every surviving
  row, and E[c·s·1{not trunc}] = c∇P(not trunc). This is a benign anti-truncation push, biased only in that direction.

### 2.5 Entropy dynamics

Lemma (Cui et al. 2025, 2505.22617, for softmax policies): ΔH ≈ −Cov_{a~π}(log π(a), Δz_a). Under PG
Δz_a ∝ η π(a) A(a). Under NPG Δz_a ∝ η A(a). When the confident tokens are also the high-advantage ones (a calibrated SFT
policy under a mode-seeking reward), the covariance stays positive and H decreases monotonically. This is exactly
Theorem 1's prediction.

Brakes available here:
- β = 1e-3 KL: negligible (§1.1).
- Clip-higher: **inactive** (F3).
- Nothing else.

Fixes and their guarantees:
- **Entropy constraint** max E ū s.t. H(π(·|I)) ≥ H₀. The problem is concave in π (linear objective, concave entropy).
  Slater's condition holds, so strong duality holds and dual ascent on the multiplier α (SAC-style auto-temperature,
  Haarnoja et al. 2018) reaches the constrained optimum π ∝ exp(ū/α). Guarantee: at the population level, in policy
  space.
- **Larger β**: π* ∝ π_ref exp(ū/β). β is a temperature in reward units. To resolve ΔIoU ≈ 0.01 without collapsing,
  β_eff ~ 0.06 in the late units (6·IoU).
- **Clip-Cov / KL-Cov** (Cui et al.): limit the step on high-covariance tokens. Empirical only; no convergence
  guarantee.
- **GED reward**: the optimum is Q, so diversity is regulated by properness. This is the only fix that is principled
  rather than a regularizer.

Two caveats. (1) On centred training items, collapse **is** convergence (train IoU 0.77 against a ceiling of about 0.73).
Adding an entropy bonus there hurts the reported metric. Entropy matters only on inputs where the policy is still wrong
(off-centre), so measure entropy per offset bucket. (2) Token entropy ≠ shape diversity: most remaining entropy is ±1
coordinate jitter. Measure E_{group} IoU(x_i, x_l) instead (M3).

### 2.6 Signal-to-noise when groups saturate

Per prompt, g_q = Cov_π(ū, s) and |g_q| ≤ σ_q‖s‖ρ_q, where ρ_q is the correlation of the reward with the score
direction. The noise is about σ_q²E‖s‖²/(G−1). So the per-prompt SNR ~ ρ_q²(G−1) **does not depend on σ_q**. The
saturation damage comes from ρ_q → 0: what reward variation remains is jitter, rasterization and floored terms,
uncorrelated with any improvable direction.

Adam normalizes gradient scale, so with low SNR the update becomes a random walk. That drifts the KL while the reward
stays flat. Observed: median KL 0.19 → 0.25–0.36, IoU slope +0.005/100 steps, gradient norm steady while the advantage
sd fell 3×. Consistent, not proof. Measure the gradient noise scale B_simple = tr Σ / ‖g‖² (McCandlish et al. 2018)
from the 16 accumulation micro-batches (M5).

### 2.7 How many samples per prompt (fixed budget)

Let V_b = tr Var_q(g_q) (between prompts) and V_w = E_q tr Var(single-sample estimate | q) (within). With cost
C = B(c_p + G c_s) (prompt/image cost c_p, per-answer cost c_s):

    Var(ĝ) = (c_p + G c_s)/C · (V_b + V_w/(G−1))   ⇒   G* = 1 + √((1 + c_p/c_s) · V_w / V_b).

(Set u = G−1, take the derivative: −(c_p+c_s)V_w/u² + c_s V_b = 0.) Shared-prefix packing lowers c_p, which lowers G*.
Across prompts, Neyman allocation (n_q ∝ per-prompt sd of the estimate) minimizes variance for a fixed total. So spend
samples on high-σ prompts. Practical version: keep a running σ̂_q per item and draw items ∝ σ̂_q. Use importance
weights 1/(Bp_q) if the uniform objective must stay unbiased; skipping saturated items biases by only O(σ_q) anyway. The
fixed ×4 disagreement oversampling is a crude Neyman rule. G also controls the O(1/G) biases of the LOO topology
frequency and of the GED U-statistics.

---

## 3. Reward composition

### 3.1 The utopia point is attainable, so scalarization cannot change the optimum

Every non-abstaining rater shape Y_j gives IoU_tc = 1 (its own signature class contains j), F5 = 1 (I ⊆ Y_j ⊆ U ⇒ D = 0),
g_area = 1 (A_min ≤ |Y_j| ≤ A_max), a valid topology with frequency f, and legality up to polygon quantization. So the
reward vector reaches its componentwise maximum at every rater shape. **All strictly monotone scalarizations (linear,
product, Chebyshev) share the maximizer set** (up to representability). Weights affect only the optimization path and
the behaviour under limited capacity. F5 and g_area have flat maximizer sets ([I, U] band, [A_min, A_max]) that
strictly contain the rater shapes, so they act as **constraints**. IoU_tc is the only term that discriminates inside the
band.

### 3.2 Multiplicative gates vs additive terms

The gate R = 1[legal]·U with U ≥ 0 gives the lexicographic order legal ≻ illegal whenever U > 0. It also zeroes the
geometric signal on illegal rollouts (∂R/∂U = g), which is correct. An additive λ·legal needs λ > sup U − inf U for the
same guarantee. Continuous products (VLM-R1's count guard) are linear scalarizations in log space (Nash-bargaining type).
Linear scalarization reaches only the convex hull of the Pareto front (Das & Dennis 1997). Chebyshev reaches all of it
(Miettinen 1999). Given §3.1 this matters only away from the optimum. Rule: hard gates for constraints, one main term,
small fixed weights for shaping.

### 3.3 Potential-based shaping is an exact no-op in sequence-level GRPO

Token MDP, deterministic transitions, terminal reward, γ = 1. The shaped return is R + Φ(s_T) − Φ(s_0) = R − Φ(s_0)
(with Φ(s_T) ≡ 0): a per-prompt constant, removed exactly by group centering (Ng, Harada & Russell 1999, plus the
baseline). Shaping helps only with per-token advantages. Φ = V^π(prefix) gives the advantage function itself (a critic,
or VinePPO-style Monte Carlo values). Dense per-vertex rewards that are not potential-based change the optimum, which is
a hacking risk.

### 3.4 Reward-hacking analysis

- **H1, topology fallback (real, bounded margin).** Emit a signature no rater drew: IoU_tc becomes IoU vs the vote, and
  r_topo = 0 (invalid). A valid rollout whose bracket is at its floor ε = 0.1 earns 2·0.1 = 0.2 advantage. With the late
  IoU weight 6, the fallback wins whenever IoU(x, vote) − max_{same sig} IoU > 0.033. Monitor `topo_any = 0` among legal
  rollouts and the share of r_topo at ε (M6).
- **H2, majority-topology pull** of the max over N_s references (§1.2d). It is opposed by r_topo. Their balance is the
  tilted law π_topo ∝ f·exp(ΔU/λ'), not f.
- **H3, sub-filter specks/holes** change IoU but not the signature: penalized, safe.
- **H4, band freedom**: anything in [I, U] has F5 = 1. That is consistent with "any expert", not a hack.
- **H5, empty or illegal** output: 0 everywhere, safe.
- **H6, rater-count dependence** (§1.2c): implicit item weighting.
- **H7, r_topo** = 1 + λ(log f − log π̂). E_π[log f − log π] = −KL(π‖f), and ∇E_π[−log π] = E[s·(−log π)]. So the
  per-sample reward log f − log π ascends −KL(π‖f) exactly. It is strictly proper for the topology law, apart from the
  clip ε, the gate (IoU ≥ 0.5) and O(1/G) LOO coupling. This is the only distribution-matching part of the current
  reward.

---

## 4. Distribution shift (centred training crops, shifted exams)

**Formalization.** Covariate shift: P_tr(I) ≠ P_te(I), with the conditional fixed up to the group action. A2
(equivariance): for a window transform t, Q_{tI} = t_#Q_I. The repo realizes this by re-rendering rater masks and caches
per window. The reward is invariant: IoU(tx, ty) = IoU(x, y) for integer translations without border clipping. F5 and
g_area are invariant because the band, T5 and the area range move with the window. Under scaling, IoU is invariant only
in the continuum, and T5 is in pixels, but the caches are recalibrated per window.

**Prop. (equivariant Bayes act).** Under A2 and invariant u, ū_{tI}(tx) = ū_I(x), so argmax ū_{tI} = t·argmax ū_I. An
optimal policy for every window exists that is equivariant.

**Importance weighting is infeasible.** The offset > 0.2 share is 0.1% in RL training rows against 68–70% in the exams.
The weight p_te/p_tr is unbounded and the Rényi divergence d₂ that controls IW variance (Cortes, Mansour & Mohri 2010)
is effectively infinite.

**Augmentation is the theorem-backed fix.** Training on the orbit minimizes the orbit-averaged risk. That equals the
test risk when the test law is the augmentation law applied to training (Chen, Dobriban & Lee 2020, JMLR). It also
reduces variance. The built `--train-views` matches the law: shifted median offset 0.247 vs 0.26 in the exams. On-policy
caveat: RL gives no gradient on inputs it never samples. Transfer to off-centre inputs through the shared LoRA has no
guarantee in either direction; here it happened to be positive (0.511 → 0.659).

**Anchor-relative coordinates.** Factorize π(o|I) = π(a|I)·π(Δ|a, I), where a is a centroid/box and Δ are offsets.
Equivariant models generalize strictly better when the symmetry holds (Elesedy & Zaidi 2021). Gains here:
1. Localization becomes a 2-number decision that can get its own reward, and misses become analysable through
   Var(R) = Var(E[R|a]) + E[Var(R|a)]. For misses the first term dominates.
2. Shape tokens get a translation-invariant target, provided visual features at the anchor are equivariant. This is an
   assumption: Qwen3-VL's 16-px patches with 2×2 merging give ~32-px token stride with 2-D RoPE, which is only
   approximately equivariant.

Honest discount: an autoregressive model already treats its first vertex as an implicit anchor, so the gain is smaller
than a naive per-vertex argument suggests. It also costs a re-SFT with a new format. **Order of operations:** run orbit
augmentation first (built, cheap, proven). Adopt anchors only if misses at offset > 0.3 stay high afterwards. A signature
to test first: regress the predicted centroid on the true centroid (M7). A slope < 1, i.e. shrinkage toward the frame
centre, is the fingerprint of a centre prior learned from centred RL.

---

## 5. Recommendation (decision record)

**Decide the semantics first.**

**S1, one outline, "matches some doctor" (current product and metric). Recommended now.**

    R(x) = 1[legal(x)] · ( IoU_tc(x) + α·L(x) ) + c·r_topo(x)
    L(x) = 1[IoU_tc(x) < 0.1] · exp(−‖c(x) − c(V)‖² / (2 s²))      # localization signal for misses only
           c(·) = centroid, V = vote mask, s = equivalent radius of U;  α ≈ 0.1  (< min useful IoU gain among hits)
    F5, g_area: constraints, not equal-weight objectives:
           either fixed small weights (≤ 0.25 of IoU's)  or  a gate  1[F5 ≥ 0.9 ∧ g_area ≥ 0.9]-style penalty
    A_i = R_i − mean_{l≠i} R_l           # RLOO / Dr. GRPO; NO per-group (per-term) std division
    (optional) divide the whole batch by one batch std: a global scalar per step, prompt weights unchanged

Rationale: Fisher-consistent for the reported metric (§1.2). It removes the step-varying implicit objective (§2.2),
gives misses a gradient (§1.4), and keeps the proper topology term (H7). Keep G = 8 unless M5 gives
V_w/V_b ≳ 50/(1 + c_p/c_s) (G* > 8). Keep β small. Train on centred plus shifted plus rescaled windows (§4). Spend
samples by σ̂_q (§2.7). Do not add an entropy bonus on saturated items. If exploration is needed on hard items, use the
entropy **constraint** with dual ascent, per offset bucket (§2.5).

Revisit if: M6 shows topology-fallback use (raise ε or make the fallback IoU = 0), or if misses at offset > 0.3 stay
above ~10% after the views run (then anchor-relative coordinates, §4).

**S2, calibrated diversity / uncertainty (e.g. present K samples, report GED).**

    r_i = 1[legal] · [ (1/k) Σ_j IoU(x_i, Y_j) − (1/(G−1)) Σ_{l≠i} IoU(x_i, x_l) ]
    b_i = mean_{l≠i} ū_1(x_l) − (1/((G−1)(G−2))) Σ_{l≠m; l,m≠i} IoU(x_l, x_m)      # x_i-free baseline
    A_i = r_i − b_i

Strictly proper (Thm 2), optimum π = Q, no entropy fix needed. Topology law: handled implicitly (IoU distinguishes
holes) or keep r_topo. Expect a lower single-sample max-IoU than S1 (§1.4). A mixture (1−γ)·S1 + γ·S2 is a deliberate
trade-off, not proper for γ < 1.

**Rejected:** Hungarian set reward (improper at finite k, M× tokens). Soft-IoU/Lovász (irrelevant for a sampled
policy). Dice (same argmax, loses the GED geometry). James–Stein baselines (≤ 12.5%). Potential-based shaping in
sequence GRPO (exact no-op).

---

## 6. What to measure (falsifiable predictions)

| # | Measurement | Prediction if this analysis is right |
|---|---|---|
| M1 | Model M_{k−1} (§1.3) against the doctors' leave-one-out ceiling, per offset bucket | Gap > 0.087 overall. Centred items fall below or at the ceiling, not above it. |
| M2 | Per step: share of groups with σ̂_k < floor_k for each term | > 80% of late groups floored for IoU/F5/g_area. The early phase is mostly unfloored. |
| M3 | Shape diversity: mean in-group IoU(x_i, x_l), against mean rater-pair IoU, per offset bucket | Model ≫ raters (over-confident) on centred items. Lower on off-centre items. |
| M4 | Views run: advantage sd, entropy and IoU slope by window kind | Shifted/rescaled rows carry most of the gradient. Their IoU slope ≫ +0.008/100 steps. Held-out offset > 0.3 misses fall. |
| M5 | Gradient noise scale from the 16 micro-batches; V_w/V_b; c_p/c_s timings | Late B_simple ≫ 256 (noise-dominated), consistent with KL drift at flat reward. G* from §2.7. |
| M6 | Share of legal rollouts with topo_any = 0; share of r_topo at ε | Non-zero and rising would confirm H1 exploitation. |
| M7 | Regression of predicted centroid on true centroid, test shifted view | Slope < 1 (centre shrinkage) before the views run, about 1 after. |
| M8 | GED (U-statistic) and HM-IoU with K = 8 samples, before/after RL | GED worsens from SFT to RL (diversity collapse) while max-IoU improves. That is the Theorem 1 trade-off, made visible. |
| M9 | Per item: ū_max − ū_1 (credit from "any expert") | Large on disagreement items. Explains why oversampling them raises σ and signal. |

## References (checked or standard)

Cui et al. 2025, *The Entropy Mechanism of RL for Reasoning LMs*, arXiv 2505.22617 · Liu et al. 2025, *Dr. GRPO*,
arXiv 2503.20783 · Mroueh 2025, *GRPO's effective loss*, arXiv 2503.06639 · Kool, van Hoof & Welling 2019 (LOO
baselines); Ahmadian et al. 2024, RLOO, arXiv 2402.14740 · Yu et al. 2025, DAPO, arXiv 2503.14476 · Weaver & Tao 2001;
Greensmith, Bartlett & Baxter 2004, JMLR (optimal baselines) · McCandlish et al. 2018, gradient noise scale, arXiv
1812.06162 · Haarnoja et al. 2018, SAC auto-temperature · Ng, Harada & Russell 1999 (PBRS) · Gneiting & Raftery 2007,
JASA (energy/kernel scores) · Kohl et al. 2018, Probabilistic U-Net (GED, d = 1 − IoU) · Gower 1971; Bouchard,
Jousselme & Doré 2013, IJAR (strict PD of the Jaccard matrix); arXiv 2608.13549 Lemma 4.1 (finite proof) · Bellemare
et al. 2017, arXiv 1705.10743 (biased Wasserstein gradients) · Berman, Triki & Blaschko 2018, Lovász-Softmax · Nowozin
2014, CVPR (optimal IoU decisions) · Das & Dennis 1997; Miettinen 1999 (scalarization) · Cortes, Mansour & Mohri 2010
(IW bounds) · Chen, Dobriban & Lee 2020, JMLR (group-theoretic augmentation) · Elesedy & Zaidi 2021, ICML (equivariance
benefit) · Rafailov et al. 2023 (KL-regularized optimum).
