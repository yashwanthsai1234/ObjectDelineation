# Math and statistics behind implementation 2.0 (research note, 2026-11-09)

Scope: the formulas used by A1, A2, A3 step 7, A4, A5, B1.9, B2.0, B3.2-B3.4, B3.7 and gates B.6. Versions: numpy 2.5.4,
scipy 1.18.1, scikit-image 0.26.0, TRL 1.17.1, PEFT 0.21.2 (all read in `.venv` site-packages).
Tags: [web] primary source, [code] file read at the pinned version, [measured] run by me on 2026-11-09 (scripts were in
the session scratchpad, so the numbers can be re-run), [derived] my algebra, checked numerically where stated,
[unverified] not checked.

## 1. Patient cluster bootstrap (A1, A2, B.6.2)

Notation: items i = 1..N, per-item value x_i (or paired difference d_i = a_i - b_i), cluster c(i) in 1..P (patient
`group_id`), N_p items in patient p, S_p = sum of x over patient p.

**Draw rule (as `signal.draws`).** Per draw b: `rng.integers(0, P, P)` picks P clusters with replacement. The draw keeps
every item of each picked cluster, so items of one patient move together [code: src/objdel/eval/signal.py:102-107].
**Vectorized equivalent** (same rng stream, so the same draws; only the summation order differs. Max |diff| vs the
concatenation means was 2.2e-16 [measured]). Bit-exact goldens must therefore keep `signal.draws`:
```python
rng = np.random.default_rng(seed)
idx = np.stack([rng.integers(0, P, P) for _ in range(B)])          # (B, P), same stream as signal.draws
M = np.zeros((B, P)); np.add.at(M, (np.arange(B)[:, None], idx), 1)  # cluster multiplicities
W = M[:, cl]                                                       # (B, N) item weights; cl = cluster index per item
mean_b = W @ x / W.sum(1)                                          # item-weighted mean per draw
pair_b = W @ (a - b) / W.sum(1)                                    # paired: same W for both arms
sub_b  = W @ (f * d) / (W @ f)                                     # subset rule, f in {0,1}; NaN if W @ f == 0
did_b  = W @ (far * d) / (W @ far) - W @ (near * d) / (W @ near)   # difference in differences, same draws
```
Memory: B x N = 2000 x 2200 float64 = 35 MB. Pitfall: `signal._ci` uses `np.nanquantile`, so a draw with an empty subset
(NaN) is silently dropped. Count the NaN draws and print the count beside every subset CI. The quantile method is numpy's
default `linear` (Hyndman-Fan type 7) [web: https://numpy.org/doc/stable/reference/generated/numpy.quantile.html].

**Percentile interval** (A1's only decision interval): [q_{alpha/2}(theta*_b), q_{1-alpha/2}(theta*_b)].
**BCa** (only if the owner asks; A1 step 5). This is scipy 1.18.1's algorithm
[web: https://raw.githubusercontent.com/scipy/scipy/v1.18.1/scipy/stats/_resampling.py, `_bca_interval`]:
- z0 = Phi^-1( (#{theta*_b < theta_hat} + #{theta*_b <= theta_hat}) / (2B) ) (ties count one half);
- a = sum_j D_j^3 / (6 (sum_j D_j^2)^{3/2}), with D_j = mean(theta_(-j)) - theta_(-j) over jackknife replicates.
  scipy writes U = (n-1) D, and the factors cancel;
- alpha_{1,2} = Phi( z0 + (z0 + z_{alpha}) / (1 - a (z0 + z_{alpha})) ). The interval reads the bootstrap
  distribution at alpha_1 and alpha_2.
- For cluster data, the jackknife must **leave out one patient**: theta_(-p) = (sum S - S_p) / (N - N_p), P
  replicates [derived]. scipy's `bootstrap` has no cluster argument. Its `paired=True` only resamples one shared
  index vector [web: same file]. So a unit test against `scipy.stats.bootstrap(method="BCa")` is valid only on
  unclustered data, as A1 says.

**CR1 cluster t (`signal.cluster_t`).** m = sum S_p / N; se^2 = P/(P-1) * sum_p (S_p - m N_p)^2 / N^2;
CI = m +- t_{P-1,(1+level)/2} se. This is Cameron-Miller eq. (11) with Stata's factor c = G/(G-1) * (N-1)/(N-K) at K = 1,
where (N-1)/(N-K) = 1 [web: https://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf eq. 12]. It
reproduces the plan's R5 val CI [0.5837, 0.6527] (se 0.01731, t_77) [measured]. With equal N_p it reduces to s_pm^2 / P,
where s_pm^2 is the sample variance of the patient means. So A1's test "equals the hand t interval of the 5 patient
means" is right [derived].
- Few-cluster caveat (P about 39 on val-select): with standard CRVE and t(G-1), simulated rejection still exceeds 0.05 at
  G = 20-30. Unbalanced clusters make it worse (G = 10, sizes 10/50: .183) [web: Cameron-Miller §VI.A, p. 24]. Our
  sizes are very unbalanced (CV of N_p = 0.99 on val [measured]). So the nominal 2.5 % one-sided error of "lo > 0" is
  optimistic.
- Optional CR2 (Bell-McCaffrey) for the intercept-only model. The closed form is
  se2^2 = sum_p (S_p - m N_p)^2 / (1 - N_p/N) / N^2. It equals the matrix CR2 to 3e-18 on a synthetic case
  [derived + measured]. It is about +9 % in variance there. [unverified: its effect on our coverage]
- **CR1 for a difference in differences** (needed when a subset has < 50 patients; the plan gives none).
  m_F = sum F_p / sum NF_p and m_E = sum E_p / sum NE_p are the far and near means of d. psi_p = (F_p - m_F NF_p) / NF -
  (E_p - m_E NE_p) / NE. se^2 = P/(P-1) sum_p psi_p^2, with df = P - 1, where P counts the patients with items in either
  subset [derived: linearisation of two ratio means]. Same pattern for any ratio-of-sums statistic.

## 2. Power and noise (B1.9, A-G3, A3 step 7)

**Design effect.** Var(mean) = sigma^2/N * Deff, Deff = 1 + (m_eff - 1) rho, m_eff = sum N_p^2 / N = m_bar (1 + CV^2)
(Eldridge, Ashby, Kerry 2006, Int J Epidemiol 35:1292) [web: https://pubmed.ncbi.nlm.nih.gov/16943232/ [unverified page
text]; the identity m_eff = m_bar(1+CV^2) is exact, derived]. Measured on R5 items with the store's `group_id`:

| split | N | P | m_bar | CV(N_p) | m_eff | sd(IoU_tc) | Deff of IoU_tc mean (CR1 var / iid var) | implied rho | Deff of miss share |
|---|---|---|---|---|---|---|---|---|---|
| val | 1084 | 78 | 13.90 | 0.99 | 27.53 | 0.2874 | 3.94 | 0.111 | 1.97 |
| test | 2200 | 158 | 13.92 | 0.94 | 26.27 | 0.2785 | 2.83 | 0.072 | 2.06 |

So 1,084 val items carry about 1084 / 3.94 = 275 independent items for a level. A paired difference has its own Deff.
Independent per-item sampling noise has Deff about 1. A patient-clustered effect (B1.9 picks patients) has Deff >> 1.

**Power of "CI lower bound > 0"** (one-sided alpha = 0.025): 1 - beta = Phi(delta/SE - z_0.975). With `cluster_t`:
`scipy.stats.nct.sf(t.ppf(0.975, P-1), P-1, delta/SE)`. Minimum effect at 80 % power: MDE = (1.96 + 0.84) SE = 2.80 SE.
SE = sqrt(Deff_d * sigma_d^2 / N).
**TOST power** (90 % CI inside +-Delta): Phi((Delta - delta)/SE - z_0.95) - Phi((-Delta - delta)/SE + z_0.95), floored at 0.
Worked numbers [estimate, from the measured inputs above and B.6.2]:
- A-G3 val TOST (Delta = 0.01). Noise only: c = 0.077 gives var(d_i) = 0.00324, Deff 1, SE = 0.00173, power = 1.00.
  If the retest shift -0.0076 seen on the 128-item pair is real, power is 0.40. Under c = 1 (var 0.0422): SE 0.0062,
  power 0. The TOST passes only if the box-to-box shift is well under 0.01.
- Gate B1 IoU rule, 25 % far-miss effect. 28 items gain about 0.68 (0.706 minus a miss's own IoU). sigma_d^2 =
  p(1-p) 0.68^2 + 0.00324 = 0.0151 (p = 28/1084), SE_iid = 0.0037. With Deff 2-4 (patient-picked), SE is 0.0053-0.0075,
  delta/SE = 3.3-2.3, power about 0.9-0.7. B1.9's simulation is the authority. This bound only checks it.

**Split-half noise algebra** (A3 step 7). Assume the K answers of an item are independent with variance sigma^2. Let
d = mean(0-3) - mean(4-7). Then Var(d) = sigma^2/2. The SD of **one** K = 8 mean is sigma/sqrt 8 = sd(d)/2. The SD of
the **difference of two independent K = 8 means** is sigma/2 = sd(d)/sqrt 2 [derived]. So `sd(d)/sqrt 2` is the retest-
difference SD, comparable to B.6.2's measured 0.0569. The per-arm SD is sd(d)/2 (0.0402 for that pair). Var bound:
a [0,1] variable with mean mu has Var <= mu(1-mu) (Bhatia-Davis), so Var(K-mean) <= mu(1-mu)/K. B1.9's c x mu(1-mu)/8
uses this bound. Its arithmetic checks: 0.00324 = 0.077 x 0.0422, and 1/0.077 = 13 [derived].
**Spearman-Brown:** r_sh = 2 r_h / (1 + r_h) is the reliability of a K = 8 mean. It is the expected correlation of two
independent K = 8 samplings of the same model on the same view [web: https://en.wikipedia.org/wiki/Spearman%E2%80%93Brown_prediction_formula].

## 3. Fair leave-one-out ceiling (A2)

Per item with doctor masks R = {m_1..m_k}, voters V ⊆ R and K answers s_j (legal flag l_j):
doctor = (1/|V|) sum_{r in V} IoU_tc(m_r ; R \ {r}); model = (1/|V|) sum_{r in V} (1/K) sum_j l_j IoU_tc(s_j ; R \ {r});
diff = model - doctor, paired per item and patient-bootstrapped. Both sides see the same k - 1 references, so the gap
from the number of references cancels. At k = 2, IoU_tc(m_1; {m_2}) = IoU(m_1, m_2) = doctor_2, so the item value is the
plain pair IoU. Caveat: IoU_tc is a max over same-signature raters, **but it falls back to the IoU vs the vote** when no
rater drew that signature [code: src/objdel/reward/score.py:7-8]. So `model_all_k >= model` (more references can only
help) is **not** guaranteed item by item. Do not assert it in a test.

## 4. GED and HM-IoU (A5)

**GED** (Kohl et al. 2018): D^2_GED = 2 E d(S,Y) - E d(S,S') - E d(Y,Y'), with d = 1 - IoU and IoU(empty, empty) = 1.
- V-statistic (the repo's `_ged`, and PHiSeg's reference code, which loops over i == j and divides by N^2 and M^2)
  [web: https://raw.githubusercontent.com/baumgach/PHiSeg-code/master/utils.py `generalised_energy_distance`]:
  `2*d[:K,K:].mean() - d[:K,:K].mean() - d[K:,K:].mean()`.
- U-statistic (unbiased): replace each within-set mean by its off-diagonal mean, `d_ss.sum()/(K*(K-1))`.
  GED_V - GED_U = U_SS/K + U_YY/k [derived]. With k = 2 doctors the rater term is half its U value. For a paired
  model-vs-model comparison on the same items the rater term cancels. The sample term (1/K) U_SS does not cancel: it
  favours the less diverse model.
- Non-negativity: the Jaccard/Tanimoto similarity is a positive-definite kernel (Bouchard, Jousselme, Dore 2013, IJAR
  54(5):615) [web: https://arxiv.org/abs/1510.02833 gives an independent proof]. The empty-set convention adds an
  orthogonal unit vector, so the kernel stays PD. So d = (1/2)||phi(x) - phi(y)||^2 and GED_V = ||mu_S - mu_Y||^2 >= 0
  (it is MMD^2 with the Jaccard kernel; HPU-Net calls GED "also referred to as MMD"). Measured: min eigenvalue of 300
  random IoU matrices -1e-15. So the repo's `max(0, .)` only absorbs rounding. A U-statistic can be negative: never clamp it.
- HPU-Net warns that GED "unduly rewards sample diversity" when samples match poorly
  [web: https://arxiv.org/abs/1905.13077 App. B.1]. That is a reason to report HM-IoU beside it.

**HM-IoU, published protocol** (HPU-Net App. B.1): "duplicate both sets so that their number of elements matches their
least common multiple", then a 1:1 Hungarian match on IoU, with the mean IoU of the matched pairs. LIDC there has always
4 graders; an empty grading is valid, and a correctly predicted absence has IoU = 1. Crops are 180 x 180 (128 x 128 fed),
with 16 samples [web: same paper, §3.2, App. C].
```python
L = K * k // math.gcd(K, k)
big = np.repeat(np.repeat(iou, L // K, 0), L // k, 1)       # (L, L); K=8, k=3 -> 24 x 24
r, c = scipy.optimize.linear_sum_assignment(big, maximize=True); hm = big[r, c].mean()
```
The repo's `hm_iou` (floor/ceil copies of each rater) equals this exactly for k in {2, 4}, but **not for k = 3**. There
it is always higher than LCM: mean +0.0146, max +0.053 on 2,000 random 8 x 3 matrices [measured]. Our k counts only the
doctors who drew, while LIDC benchmarks use 4 graders with empties. So our numbers are not comparable to the
published LIDC GED/HM-IoU values for two reasons.

## 5. Answer choice (A4)

- The medoid is minimum-Bayes-risk decoding with utility IoU over the legal pool: argmax_i (1/(n-1)) sum_{j != i}
  IoU(s_i, s_j) [web: https://arxiv.org/abs/2310.01387]. Including j = i adds the same +1/(n-1)-scaled constant to
  every candidate, so the argmax does not change [derived]. It is the sample 1-median under the Jaccard distance.
- The pixel strict-majority vote is the exact minimiser of the summed Hamming distance sum_j |m XOR m_j| over all masks
  (pixelwise median). It is **not** the Jaccard medoid or median [derived]. With even n_legal, pixels at exactly n/2
  are ties. `2*sum > n` drops them (a valid choice; state it in the docstring).

## 6. RLOO, Dr. GRPO and TRL 1.17.1 (B3.4)

- RLOO: A_i = R_i - (1/(G-1)) sum_{j != i} R_j [web: https://arxiv.org/html/2402.14740 §2.3]. It equals
  G/(G-1) (R_i - mean R). Dr. GRPO's A_i = R_i - mean R is RLOO up to that factor
  [web: https://arxiv.org/abs/2503.20783 App. A, "equivalent to RLOO up to a scaling factor"]. sum_i A_i = 0 exactly in
  real arithmetic. If all R are equal, A = 0 with no special case.
- TRL 1.17.1 `_generate_and_score_completions`, `multi_objective_aggregation="sum_then_normalize"`:
  rewards = nansum(w_k r_k), advantages = rewards - group nanmean. With `scale_rewards != "none"` it divides by
  (std + 1e-4), where std = `nanstd` with the n/(n-1) correction. NaN advantages become 0
  [code: .venv trl/trainer/grpo_trainer.py:2880-2929; web: https://github.com/huggingface/trl/blob/v1.17.1/trl/trainer/grpo_trainer.py].
  So with `scale_rewards="none"`, an RLOO A returned as the only reward (weight 1) comes back unchanged up to float32
  rounding. The tensor is float32 (relative 6e-8) and the re-centring shift is about 1e-9.
- Zero or near-zero spread: GRPO's std division turns a group with std 5e-4 into advantages of about +-0.8. This is the
  Dr. GRPO "question-level difficulty bias". RLOO/Dr. GRPO keep such a group near 0. The v5 drop rule then zeroes it
  when every raw term has ptp < 0.01. Never return `None` from the reward: an all-None row becomes NaN and is
  excluded from the baseline [code: grpo_trainer.py:2876-2883].
- `loss_type="dr_grpo"`: loss = sum(per_token_loss * mask) / (B_micro * max_completion_length) / grad_accum
  [code: grpo_trainer.py:3285-3289]. Zeroed groups still count in B_micro, so they dilute the step. That is intended,
  and it is why B3.7 rule 10 measures S.
- Adam scale argument (B3.4): with beta -> beta/r and eps -> eps/r, the loss scales by 1/r and
  m/(sqrt v + eps) is invariant. Only `max_grad_norm` breaks exact invariance. It never binds (0.0018-0.022 vs 1.0)
  [derived]. Check: 0.0048/sqrt(17,432,576) = 1.15e-6 = 115 eps. An element at 1.15 eps gets 1.15/2.15 = 0.53 of its
  step [derived].

## 7. Location reward (B3.2, B3.3, B3.7)

- loc_raw = exp(-x^2 / (2 w^2)) or 1/(1 + x^2/w^2), x = d / r, d = centroid distance in px, r = sqrt(|U|/pi).
  Both are bounded in (0, 1], equal 1 at d = 0, and do not change under a common rescale of image and masks (x is a
  ratio). Pixel-index centroids on both sides make the +0.5 offsets cancel.
- Signal range: the Gaussian slope x e^{-x^2/2} is 0.033 at x = 3 and 1.9e-5 at x = 5. The Cauchy slope
  2x/(1+x^2)^2 is 0.060 and 0.015. Width rule 5 (median x_med gets loc >= 0.1): Gaussian w >= x_med / sqrt(2 ln 10) =
  x_med / 2.146; Cauchy w >= x_med / 3 [derived].
- Order lemma (B3.3). For iou < m: R_loc = iou + alpha (1 - iou/m) L, with L = g*topo*loc_raw <= 1. This is linear in
  iou with slope 1 - alpha L/m >= 0 iff alpha <= m/L. Its sup over [0, m) is max(alpha L, m) = m when alpha <= m. Any
  hit has iou >= m. So **miss <= hit for equal f5 and g_area, with equality possible when alpha = m** (a miss at
  IoU 0 with L = 1 ties a hit at IoU m) [derived]. Use alpha < m, or test T11 with `>=`.
- Not potential-based shaping. Ng-Harada-Russell invariance needs F = gamma Phi(s') - Phi(s) with a constant terminal
  potential [web: https://www.cs.utexas.edu/~shivaram/readings/b2hd-NgHR1999.html]. One completion is one episode, so a
  term that depends on the final answer is a real reward change. Its safety is the order lemma: hits keep loc = 0 and
  their order, so argmax R is still a hit when a hit is reachable. Under the KL-regularised optimum
  pi* ∝ pi_ref exp(R/beta), the miss/hit odds move by at most exp(alpha/beta) [derived].
- Anti-gaming gates: g_area = 0 outside the doctors' area range (speck, huge outline); topo_freq > 0 (two blobs around
  the centroid); ramp (1 - iou/m) makes R continuous at m. Drop rule: B3.4 tests ptp of the **gated** loc. T14 states
  loc_raw spans 0.04. That survives the 0.01 test only if ramp*g*topo >= 0.25. Pin g_area = 1 in T14 [derived].

## 8. bf16 merge rounding (B2.0)

- bf16: 8 significand bits, ulp(W) = 2^(floor(log2|W|) - 7), eps = 2^-7. fp16: 11 bits, ulp 8x smaller. W is on the
  grid, so with t = d/ulp(W) (inside one binade): e = round(W + d) - W - d = ulp (rint(t) - t). |t| < 1/2 is erased
  (e = -d). Ties at |t| = 1/2 go to even. **At W = +-2^E with d toward zero, the spacing is ulp/2 and the erase
  threshold is ulp/4** [derived from IEEE 754 round-to-nearest-even].
- Error model: rel^2 = sum ulp_i^2 (rint(t_i) - t_i)^2 / sum d_i^2. For t ~ N(0, s^2): erased share = 2 Phi(1/(2s)) - 1
  and rel ≈ min(1, 1/(sqrt(12) s)). Measured on synthetic bf16 W ~ N(0, 0.02^2): s = 0.5: erased 0.683, rel 0.5745
  (model 0.577); s = 1: 0.383 / 0.2885; s = 2: 0.197 / 0.1443; s = 8: 0.05 / 0.0358 [measured]. The plan's step-50 row
  (0.692 / 0.605) fits s ≈ 0.48. The fp16 column (0.128) is above the Gaussian prediction (~0.074): d/ulp is
  heavy-tailed across |W| [estimate].
- Exact erased test (no ulp code): `erased = (W + d).to(bf16) == W` on non-zero d. It handles ties and binade edges.
  The seed script's ulp (`nextafter` toward +inf) takes the downward spacing for negative W and the upward spacing for
  positive W.
- **What vLLM receives.** `lora_merged_params` pushes `p.data + d.to(p.dtype)`, that is round(W + round(d)) (double
  rounding) [code: src/objdel/train/fast.py:717]. PEFT's merge does `weight.data += delta_weight` with fp32 adapters
  (PEFT autocasts adapters to fp32 by default), that is round(W + d) [code: .venv peft/tuners/lora/layer.py ~1007-1008;
  tuners_utils.py `autocast_adapter_dtype=True`]. The two differ on 0.15 % (s = 0.5) to 1.9 % (s = 8) of elements
  [measured]. B2.0 should compute the pushed values, not round(W + d).
- Merge/unmerge drift (TRL's path; U-SY2). unmerge is `weight.data -= delta_weight.to(orig_dtype)` [code: layer.py
  ~1035]. Then round(round(W + d) - round(d)) != W when |e + (d - round d)| crosses ulp/2, or when W + d crosses a
  binade. Measured per sync: 1.9e-5 (s = 0.1), 0.27 % (0.5), 0.40 % (1), 0.71 % (2), 2.6 % (8) of elements move by an
  ulp. Planted drift for U-SY2: W = 2 - 2^-7 = 1.9921875 (bf16), d = 15 * 2^-10 (fp32). Merge gives 2.0 and unmerge
  gives 1.984375 != W, a binade crossing. W = 1 with d = 3 * 2^-9 does **not** drift (merge 1.0078125, unmerge 1.0)
  [measured, torch bf16 in-place ops as PEFT runs them].
- B2.0 planted case: W = 1, d = [[2^-10, 0], [0, 2^-6]]. bf16: 1 + 2^-10 -> 1 (erased), 1 + 2^-6 exact. Erased share
  1/2. rel = 2^-10 / sqrt(2^-20 + 2^-12) = 1/sqrt(257) = 0.0623783. fp16: both exact, erased 0, rel 0 [derived].

## 9. Digital topology (checker, D2/D4b, 3D voxel phase)

- 2D, (8,4): b0 = 8-connected foreground components. b1 = 4-connected components of the background of the 1-px
  zero-padded mask, minus 1. chi = b0 - b1. 2D (4,8) is the dual. The closed cubical complex (each pixel a closed unit
  square) gives 8-connectivity: chi = V - E + F.
- Bit quads on the zero-padded mask, over all 2x2 windows: Q1 = one fg, Q3 = three fg, QD = two diagonal fg.
  chi_8 = (Q1 - Q3 - 2 QD)/4 and chi_4 = (Q1 - Q3 + 2 QD)/4 (Gray 1971, IEEE Trans. Comput. C-20:551) [unverified
  citation]. Checked: identical to b0 - b1, V - E + F and `skimage.measure.euler_number(m, 2 / 1)` on 300 random masks [measured].
- 3D, (26,6): b0 = 26-components (`ndi.label(m, np.ones((3,3,3)))`). b2 = 6-components of the padded background minus 1.
  chi = V - E + F - C of the closed voxel complex = `euler_number(m, connectivity=3)`. b1 = b0 + b2 - chi. (6,26) is the
  dual: b0 6-conn, b2 26-conn background, chi = `euler_number(m, 1)`. skimage defines only 1 and 3 in 3D
  [web: https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.euler_number]. Checked on 150
  random volumes and a voxel ring (b0 1, b2 0, chi 0, so b1 1) [measured].
- Pitfall: `ndi.label` defaults to the cross structure (4/6). Always pass the structure for 8/26.

## 10. Boundary metrics (A5)

- Boundary IoU: band(m, d) = m AND NOT erode(pad(m, 1), 3x3, d iterations), then crop. That is the official
  `mask_to_boundary` (`copyMakeBorder` 1 px zeros, `np.ones((3,3))`, `iterations=dilation`, crop, `mask - erode`)
  [web: https://github.com/bowenc0221/boundary-iou-api/blob/master/boundary_iou/utils/boundary_utils.py].
  BIoU = |Ba ∩ Bb| / |Ba ∪ Bb|. The official file has no empty-mask rule, so the plan's 1.0 / 0.0 is our convention.
- Boundary F with tolerance (DAVIS): bound_pix = ceil(0.008 * ||shape||) if bound_th < 1 else bound_th. Boundaries
  come from `seg2bmap`; dilate with `disk(bound_pix)`. P = |B_pred ∩ dil(B_gt)| / |B_pred|, R = |B_gt ∩ dil(B_pred)| /
  |B_gt|, F = 2PR/(P+R) (0 if P+R = 0). Empty cases: both empty -> P = R = 1; only pred empty -> P = 1, R = 0
  [web: https://github.com/davisvideochallenge/davis2017-evaluation/blob/master/davis2017/metrics.py]. EDT form:
  P = mean(EDT_{~B_gt}[B_pred] <= tau), because `disk(r)` holds x^2 + y^2 <= r^2 [derived].
- The repo's F5 is **not** this F: F5 = 1/(1 + (2D/P)/T5^2), a smooth depth-beyond-band score [code: score.py:9]. NSD
  is surfel-length-weighted [code: metrics.py docstring]. Never put F5 beside a published boundary F.

## 11. Plan errors (line numbers in Implementation.md)

1. **L1772 (A3 step 7) and L4766 (B1.9).** `item_sd = sd(d)/sqrt 2` is labelled "per-item K = 8 noise", but it is
   the SD of the **difference of two K = 8 means** (= the SD of one K = 4 mean). The SD of one K = 8 mean is sd(d)/2
   (§2). If B1.9 feeds `item_sd` in as per-arm noise, the variance doubles. Use item_sd^2/2 per arm. Fix the label
   (and L2251's description). `mean_sd` has the same meaning: noise of a paired difference.
2. **L5306 (B3.3) and L5294/L5392.** "alpha_loc <= miss_iou keeps a hit above a miss" holds only as `>=`. With the
   defaults alpha = m = 0.1, a miss at IoU 0 with full location ties a hit at IoU 0.1 (§7). Either require
   alpha < m (raise on `>=`), or state "not below" and test T11 with `>=`.
3. **L1867 / L1260 (A5) and the `metrics.py` HM-IoU docstring.** For k = 3, the floor/ceil assignment is not the
   published LCM protocol. It reads +0.015 on average (up to +0.053) [measured]. Use the LCM duplication (24 x 24 at
   k = 3; cheap). Also state that published LIDC values use 4 graders with empty masks.
4. **L1439 (A1 test).** "Two clusters of unequal size: nodule CI width > item CI width" is not a theorem. It needs
   positive within-cluster correlation, and a 2-cluster bootstrap has only 3 distinct means. Build the fixture with
   >= 5 clusters of equal within-cluster values.
5. **L4753/L4769 with L5657 (B1.9 rule 2).** The simulated effect (25 % of far misses converted) lands the B1 far-miss
   point at 0.75 x 27.3 % = 20.475 %, against the threshold 20.5 %. The margin is 0.1 item of 406. With any noise,
   P(rule 2) is about 0.5 by construction, so joint P(pass) <= ~0.5 and the "< 0.5 -> owner decides" trigger fires
   almost surely [estimate]. B3.7 states this effect for B3; B1.9 should state it too, and the "smallest effect with
   joint P >= 0.8" line is the useful number.
6. **L5655 (Gate B1 decomposition).** "+0.018" credits each converted miss with the full found IoU 0.706. The gain is
   0.706 minus the miss's own K-mean IoU (< 0.1). With 0.02-0.03 the gain is about +0.0175 [estimate]. Minor.
7. **L1745 / L2290 (A-G3 view identity).** It compares the point Pearson r_obs with the lower 95 % bound of r_sh, and
   ignores r_obs's own sampling error. So the false-fail rate on a true same view is not controlled (it is not 2.5 %)
   [derived]. Better: bootstrap r_obs - r_sh with the same patient draws, and fail only if its upper bound < 0.
8. **L4811 (B2.0) "erased = |d| < ulp(W)/2".** This misses ties (half erased under round-half-even) and the ulp/4
   threshold at power-of-two W with d toward zero. Its relative error uses round(W + d), while vLLM gets
   round(W + round(d)) (§8). Use the exact test `round(W + d) == W` on the pushed values. The planted test values are
   in §8.
9. **Gap, B.6.2 / A1 step 5.** If a subset (far, near) has < 50 patients, A1 requires `cluster_t`, but the
   difference-in-differences rule has no CR1 formula. Use the psi_p formula in §1. Also count NaN (empty-subset) draws.
10. **L5294 / T14 (L5393).** The drop rule tests the gated loc. T14's loc_raw span of 0.04 survives only if
   ramp x g_area x topo >= 0.25. Pin g_area = 1 and IoU = 0 in T14.
11. Checked and right: CR1 formula and the equal-N_p test (L1412, L1442, L5592); [0.5837, 0.6527]; RLOO =
   G/(G-1)(R - mean) (L5320); the B3.4 worked check +0.057; exp(-2) = 0.135 and exp(-4.5) = 0.011; Adam eps arithmetic;
   c = 0.077 and "13 times"; Spearman-Brown; DAPO soft penalty arithmetic (-0.905 at 1,181 tokens); TOST with a 90 % CI;
   k = 2 LOO equals plain IoU; boundary IoU band and the 16-px ring; B2.0 planted share 0.5.
