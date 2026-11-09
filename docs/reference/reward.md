# Reward reference: reward v4

This page defines every term of the reward as it is implemented, the gates, and the way terms combine into one advantage.
Sources: `src/objdel/reward/{spec,cache,score,advantages}.py`, `cuda/gpucheck.cpp` (`ObjTerms`), `src/objdel/train/rewards.py`.
The design reasoning is in [ADR-0002](../architecture/decisions/0002-set-valued-multi-rater-reward.md). The lab evidence is in `experiments/phase1/lab/t_reward/results/ITERATION4.md`.

## 1. Pipeline

```text
  answer text
     |  parse (strict)                    format_ok = 0, all terms 0 on failure
     v
  outline --> checker: even-odd mask, 13 defect counts       legal = no defect
     |
     v   per object, against its expert table (resident on the device)
  iou_tc, f5, g_area, signature, topo_freq, area
     |
     v   per answer (several objects: mean, min, joint topology)
  terms: format_ok legal iou_tc iou_min f5 g_area topo_any topo_freq sig_id
     |
     v   per group of 8 answers
  advantages_v4  ->  one number per answer
```

A job for the scorer is the tuple `(text, truncated, names, caches, W, H)`.
`caches` maps each object name to its `ExpertCache`.

## 2. Frozen constants (`REWARD_V4`)

Changing any value is a change of the reward definition. It needs a new lab round and it invalidates every stored expert table, because `fingerprint()` includes these values.

| Constant | Value | Meaning |
|---|---|---|
| `weights` | 0.3, 0.3, 0.1, 0.3 | For `iou_tc`, `f5`, `legal`, `g_area` |
| `floors` | 0.05, 0.05, 0.25, 0.05 | Sigma floors of the same terms |
| `c_topo` | 2.0 | Weight of the centred topology reward |
| `lam` | 0.5 | Lambda of the topology bracket |
| `eps` | 0.1 | Clip floor of the topology bracket |
| `group` | 8 | Answers per image |
| `drop_spread` | 0.01 | A group whose terms all span less than this, with one topology answer, gets advantage 0 |
| `iou_gate` | 0.5 | A topology answer counts only with `iou_tc` at least this |
| `expert_score` | 0.95 | Expected F5 of a new rater (s*) |
| `tau_floor` | 1.0 px | Floor of tau and of the F5 tolerance |
| `tau_quantile` | 95.0 | Percentile for tau |

`fingerprint()` is the SHA-256 of `REWARD_V4`, the size rule probed at areas 0, 999, 1001, 10^4 and 10^6, and the versions of numpy, OpenCV and scipy.

## 3. Expert table (`ExpertCache`)

Built once per object of one item from the rater stack `(k, H, W)`, `k >= 2`. A stack with one rater raises `ValueError`.

| Field | Definition |
|---|---|
| `union` U | Pixels drawn by at least one rater |
| `core` I | Pixels drawn by all raters |
| `vote` | Pixels drawn by more than half: `2 * count > k` |
| `perimeter` P | Boundary pixels of U (at least 1). A boundary pixel has a 4-adjacent background pixel |
| `dist_out_union` | Exact distance to U, for pixels outside U |
| `dist_in_core` | Exact distance to the complement of I, for pixels inside U |
| `tau` | `max(1, P95 of pooled rater-pair contour distances)` in pixels |
| `t5` | The F5 tolerance, see section 4.2 |
| `rater_areas` | Pixels per rater. Gives A_min and A_max |
| `min_area` | `max(5, ceil(|U| / 200))` |
| `rater_sigs` | Size-filtered signature of each rater |
| `voters` | Raters with a non-empty signature. All raters when every one abstains |
| `same_sig` | Signature to voting raters who drew it. The rater topology histogram |

The distance transform is scipy's exact Maurer transform: integer squared distances, a float64 square root, rounded to float32. OpenCV's x86 distance transform is not used because it is not bit-reproducible across CPUs.

**Abstaining raters.** A rater whose every part is below the size filter has the empty signature and abstains. Such a mask is the cross-section of a real nodule on an end slice. It is not an answer "nothing here". The rater leaves the topology histogram and the IoU_tc rater set. It stays in the band, the vote and the area range.

**Storage.** `save` writes a compressed `.npz` with the rater masks, tau, T5, the fingerprint and a SHA-256. `load` raises `ValueError` for corruption or a wrong `expected_digest`, and `StaleCacheError` (a `ValueError`) when the fingerprint changed. `load(..., rebuild=True)` recomputes tau and T5 from the verified masks.

## 4. Terms of one object

Let M be the answer's mask, A its area in pixels, and the raters R_1 to R_k.

### 4.1 `iou_tc`: IoU with the topology condition

```text
iou_tc = max over voting raters r with signature(r) == signature(M) of IoU(M, R_r)
       = IoU(M, vote)    if no voting rater drew that signature
IoU(a, b) = |a ∩ b| / |a ∪ b|      (1.0 when the union is empty)
```

The signature is `signature(mask, min_area)`: holes per part, descending, after the size filter.

### 4.2 `f5`: boundary term against the band

```text
D  = sum over error pixels of (d - 1/2)
     spill:     pixels of M outside U, d = dist_out_union
     core miss: pixels of I outside M, d = dist_in_core
f5 = 1 / (1 + (2 D / P) / T5^2)
```

Pixels between core and union are never errors. A new rater scores 0.95 on average, at any object size and rater count, because T5 is self-calibrated:

```text
q_j   = 2 D_j / P_j      for each rater j scored against the band of the other k - 1 raters
q_k   = mean_j(q_j) * S_k / S_(k-1)
T5    = max( sqrt(q_k / x*),  tau_floor * sqrt(S_k / x*) / (1.96 * sqrt(2)) ),    x* = 1 / 0.95 - 1
```

The second argument is a floor that keeps T5 above a fixed size when the raters agree almost exactly (`f5_tolerance` in `cache.py`).

`S_m` is the second moment of a new rater's excursion beyond the band of m i.i.d. N(0, 1) raters. It is frozen for m up to 8 (`S_M`) and computed by quadrature above.

### 4.3 `g_area`: area guard

```text
g_area = min(1, A_max / A) * min(1, A / A_min)       0 when A = 0
```

A_min and A_max are over all raters, with zero tolerance. It punishes inflation and deflation outside the rater range.

### 4.4 Signature and topology frequency

```text
sig        = signature of M at min_area
topo_freq  = (voting raters whose signature == sig) / (voting raters)
```

## 5. Terms of one answer

An answer has one object per structure name: `nodule` for LIDC, `disc` and `cup` for RIGA.

| Term | Definition |
|---|---|
| `format_ok` | 1 if the text parsed and is not truncated, else 0 |
| `legal` | 1 if format_ok and no defect in any object, else 0 |
| `iou_tc`, `f5`, `g_area` | Mean over the objects |
| `iou_min` | Smallest `iou_tc` over the objects. Feeds the validity gate, so a good disc cannot carry a missed cup |
| `topo_freq` | Share of voting raters whose signatures equal the answer's for every object. A rater of one object's table is rater r of the others. Tables with different rater counts raise `ValueError` |
| `topo_any` | 1 if `topo_freq > 0` |
| `sig_id` | CRC32 of the tuple of per-object signatures. Independent of `PYTHONHASHSEED`. -1 for an unparseable answer |

For an unparseable or truncated answer: `format_ok = legal = 0`, every other term 0, `sig_id = -1`.
For a parsed but illegal answer the terms are kept for logs. The legality gate zeroes them in the advantage.
An empty object list is a legal empty mask.

## 6. Gates

| Gate | Function | Rule |
|---|---|---|
| Legality | `legal_gate(T)` | True when `legal >= 0.5`. The advantage sets every geometric term of an illegal answer to 0 first |
| Topology validity | `topology_valid(T)` | `legal` and `topo_any > 0` and `iou_min >= 0.5` and `topo_freq > 0` |
| Flat group | in `advantages_v4` | A group whose terms all span less than 0.01 and whose answers share one `sig_id` gets advantage 0 (dynamic sampling) |

Evaluation uses the same gates: an illegal answer scores 0, and "passed" means `topology_valid`.

## 7. Topology reward `r_topo`

For an answer with topology answer a in a group of G = 8:

```text
valid     = topology_valid
pi_loo(a) = max(n(a) - 1, 0.5) / (G - 1)         n(a) = answers of the group with topology a
b(f, G)   = log f - E_{N ~ Binomial(G-1, f)} [ log( max(N, 0.5) / (G - 1) ) ]
r_topo    = valid * max( eps, 1 + lam * ( log f(a) - b(f(a), G) - log pi_loo(a) ) )
```

f(a) is `topo_freq`. The term compares the rater frequency of a with the model's own frequency, estimated leave-one-out inside the group.
`b` removes the bias of that estimate: at pi equal to f no rater answer is favoured. `eps` keeps every valid answer above every invalid one.
The term is mean-centred in the group and never divided by its standard deviation.

Worked values (`topology_reward`, G = 8, all answers legal with `iou_min` 0.8; raters split 0.75 and 0.25 between signature (0,) and (1,)):

| Group | r_topo of the (0,) answers | r_topo of the (1,) answers |
|---|---|---|
| 6 answers (0,), 2 answers (1,) | 1.010 | 1.207 |
| 8 answers (0,) | 0.842 | none |
| 2 answers (0,), 6 answers (1,) | 1.815 | 0.402 |

In the second row all answers are equal, so after mean-centring r_topo adds nothing. In the first row the under-sampled minority topology gets a positive centred value. In the third row the under-sampled (0,) answers win.
An answer whose topology no rater drew has `r_topo = 0` (a group with seven valid answers at 0.919 and one invalid answer gives 0.919 and 0.000).

## 8. Advantage

```text
X_k    = term_k * legal_gate                  k in (iou_tc, f5, legal, g_area)
z_k    = (X_k - group_mean) / max(group_std, floor_k)
A      = sum_k w_k z_k  +  c_topo * ( r_topo - group_mean(r_topo) )
```

- There is no batch-level z-score. One answer's geometric advantage is bounded by `sum_k w_k * sqrt(G - 1)`. `c_topo = 2` dominates it on purpose.
- The floors keep a near-constant term from amplifying numeric noise.
- A flat group gets A = 0.
- The result is a flat array `(B * G,)`, prompt by prompt. The length must be a multiple of G, else `ValueError`.

### Worked example

Raters drew discs of radius 20, 22 and 24 pixels on a 100 x 100 frame. Eight answers are 16-gons of radius r grid units (10 grid units = 1 pixel).

| Answer | r | iou_tc | f5 | g_area | advantage |
|---|---|---|---|---|---|
| 1 | 200 | 0.944 | 0.985 | 0.976 | 0.467 |
| 2 | 210 | 0.924 | 0.999 | 1.000 | 0.485 |
| 3 | 220 | 0.952 | 1.000 | 1.000 | 0.534 |
| 4 | 230 | 0.934 | 1.000 | 1.000 | 0.503 |
| 5 | 240 | 0.955 | 0.994 | 1.000 | 0.532 |
| 6 | 150 | 0.547 | 0.581 | 0.547 | -1.312 |
| 7 | 330 | 0.537 | 0.216 | 0.537 | -1.744 |
| 8 | 220 | 0.952 | 1.000 | 1.000 | 0.534 |

All answers are legal (`legal` is 1). All have the signature (0,), as every rater does. So `r_topo` is equal for all, and its centred value is 0.

Check answer 7 by hand. Group means are 0.8432 (iou_tc), 0.8469 (f5), 1.0 (legal), 0.8826 (g_area).
Group standard deviations are 0.1740, 0.2746, 0 and 0.1966. The `legal` term uses its floor 0.25.

```text
z(iou_tc) = (0.5373 - 0.8432) / 0.1740 = -1.758
z(f5)     = (0.2156 - 0.8469) / 0.2746 = -2.299
z(legal)  = 0
z(g_area) = (0.5373 - 0.8826) / 0.1966 = -1.756
A = 0.3 * (-1.758) + 0.3 * (-2.299) + 0.1 * 0 + 0.3 * (-1.756) = -1.744
```

This matches the value computed by `advantages_v4`.

## 9. Integration with TRL (`objdel.train.rewards.TRLReward`)

TRL calls every reward function once per generation batch. `TRLReward.funcs()` returns eight functions and weights `[0, 0, 0, 0, 0, 0, 0, 1]`:

| Function | Weight | Returns |
|---|---|---|
| `format_ok`, `legal`, `iou_tc`, `f5`, `g_area`, `topo_any`, `r_topo` | 0 | One term each, for logging as `rewards/<term>/mean` |
| `advantage_v4` | 1 | `advantages_v4` |

- With `multi_objective_aggregation="sum_then_normalize"` and `scale_rewards="none"`, TRL's advantage is the weighted sum minus its group mean. The aggregator's group mean is already 0, so TRL's advantage equals `advantages_v4` to float32 precision (test R2, `tests/test_train_trl.py`).
- `iou_tc`, `f5` and `g_area` are logged times `legal_gate`, so curves and the best-checkpoint metric do not credit illegal output.
- The first call scores the whole batch in one checker call. The other functions reuse it when passed the same `completions` list object.
- `prefetch` starts scoring on a background thread as soon as vLLM returns. The reward call uses the result only for exactly the same answers.
- A call must hold whole groups of 8 in order, else `ValueError`.
- Typed logs: `error/<parser code>` and `defect/<class>` shares per batch.
- A truncated answer is a format failure even if its prefix parses. Truncated means the last token id is not a stop id (`truncated(ids, stop)`).

Guard switches: `Controller.switch(reward_weights={"g_area": 0.5})` gives a weight-0 logger a weight w. It adds `w * (term - group mean)` to every advantage. Use it to raise the pressure on one term during a run.
