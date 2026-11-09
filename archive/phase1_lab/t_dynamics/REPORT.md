# t_dynamics: GRPO dynamics of set-valued topology rewards and a typed-error controller

Saved by the coordinator from the agent's final message; the agent could not write this file.
Full tables: `results/a_tables.md`, `results/b_tables.md`. Figures: `figs/`. Tests: 29 passing (`pytest -q`).
Toy: exact GRPO on per-prompt softmax policies over (β0, β1); G = 8; clip 0.2/0.28; 20 seeds; 95 % CIs.

## A. Does an "any rater counts" topology reward collapse the topology distribution? Yes.

- **Collapse (TV to the raters' histogram):** 0.11 at SFT → 0.30 after 300 visits → 0.33 after 1500. The minority mode keeps only 11 % of its rater share.
  - The share-of-raters, majority and log-frequency rewards collapse as badly or worse (TV 0.30–0.32).
  - **Mechanism:** while invalid answers exist, every group pushes all valid answers up in proportion to how often they were sampled (rich get richer). Std-normalization keeps that push at order 1.
- **IoU gap:** a systematic gap between topologies makes the collapse total (TV 0.45 = point mass). The repo's IoU is scored against the strict-majority vote mask, which creates exactly such a gap.
  - Under the repo's per-term z-scoring (GDPO), every fix fails once the gap reaches 0.1.
- **Recommended topology term:** `r_topo = 1[valid] · (1 + λ·(log f(a) − log π̂_loo(a)))`.
  - λ ≈ 0.5; π̂ is the leave-one-out share of the group giving the same answer; G ≥ 8 (16 is better).
  - **Mean-centre it** (Dr. GRPO), or add it before a joint z-score. Never z-score it per term.
  - **Score IoU against the raters who share the sample's topology.**
  - Result: invalid mass → 0 and TV ≈ 0.05–0.10, robust to a 0.1 IoU gap, a miscalibrated SFT and parameters shared across prompts.
- **Rejected:**
  - Coverage bonus (Multi-Answer RL): spreads answers evenly over the rater set (TV ≈ 0.20).
  - KL to SFT: needs a calibrated SFT, keeps invalid mass alive, and fails under an IoU gap.
- **Metric to report:** per-image TV or JS to the rater histogram, plus minimum-mode retention. Never histograms averaged over images.

## B. Level 2 controller

- **Learnability sampling p(1−p) beats uniform** (about 2× faster to the same error).
  - Band-pass difficulty filtering on stale pass rates is worse than uniform in 10 of 12 scenarios.
  - Typed samplers match learnability on final error and gain only 3–7 % in area under the error curve.
  - The "upweight stalled types" rule is never better, and is worse with an ambiguous error type.
  - Per-type balancing helps only when a high-prevalence error type plateaus because of ambiguity.
- **The typed alarm is the valuable part.** Rule: "train reward up while a held-out typed error rises above its best".
  - Fast hacks: caught in 98.6–100 % of runs.
  - False alarms: 0 % in 9 of 10 hack-free scenarios.
  - Effect: final error drops from 0.59 to about 0.10 once the guard it names is switched on.
  - Limits: it lags the held-out optimum by about 270 steps and misses slow hacks, so pair it with rollback to the best-validation checkpoint.
  - Threshold trade-off (lenient/medium/default/strict): see `results/b_tables.md`.
- **Recommended design:** learnability sampler + typed held-out monitoring with alarm + rollback, where the alarm names the guard to switch on.

## Caveats

- **The topology class is a marginal over polygon sequences.** Only group-count estimators of π̂ are deployable, and they are noisier.
- **IoU-gap sizes (0.03, 0.1) are guesses.** Measure the real gap on the multi-rater fixtures first.
- **The curriculum simulator's dynamics are invented**, so they are indicative only.
