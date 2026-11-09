# t_reward: reward core on real expert masks (iteration 1)

Saved by the coordinator from the agent's final message; the agent could not write this file.
Code: `rewardcore/`. Tests: `tests/` (52 passing). Scripts: `scripts/`. Outputs: `results/`.
Machine was heavily loaded during timing, so ratios are valid; re-check wall-clock K4 on an idle host.

## KPIs

| KPI | Target | Measured | Verdict |
|---|---|---|---|
| Topology exact | 0 mismatches | 0 / 10,000 random (1.39 M components) and 0 / 8,738 fixture masks; certificate χ8 = Σ(1 − holes) | PASS |
| K8: abs(ρ) of band vs area, leave-one-rater-out | ≤ 0.10 | F5 +0.107 [0.088, 0.126]; within subsets ≤ 0.04; current band +0.531 | borderline FAIL |
| K8: median LOO score | ≥ 0.90 | F5 0.973 | PASS |
| K9: topology ranking | ≥ 95 % | 98.7 % of 12,127 pairs; 100 % at IoU ≥ 0.5; IoU alone 53.9 % | PASS |
| K3: batch of 256 | ≤ 250 ms | 78 ms (8 threads, mixed); 106–129 ms (RIGA-only / processes) | PASS |
| K4: p99 per rollout | ≤ 5 ms | CPU p99 1.2–2.1 ms; wall-clock under load inconclusive | PASS (CPU) |
| K10: determinism | 100 % | identical SHA-256 across 8 execution settings | PASS |

## Decisions

**Topology.**
- Crop to the bounding box, pad by 1 px, label with cv2 Spaghetti (8-connected foreground, 4-connected background).
- Assign each hole to the component on the pixel directly above its first pixel.
- 5–24× faster than the current path.

**Band.** H5 is rejected: rater disagreement grows with object size (log-log slope 0.46–0.83), partly because of resolution mixing. Recommended formula:
- F5 = 1 / (1 + (2D/P)/T5²), where D is the distance-weighted error beyond the band and P is the union's contour length.
- T5 is self-calibrated per image from leave-one-rater-out scores of the raters themselves, with a rater-count correction (Gaussian model, S_m constants).

**Advantages.** The current GDPO normalization has the "dark room" problem: 1e-3 noise on a constant term moves advantages by up to 1.24 and cuts ordering accuracy to 80.8 %. Use `gdpo_floor_gated`:
- legality as a gate;
- topology gated on IoU ≥ 0.5;
- per-term σ floors (0.05 for IoU/band, 0.25 for topology/legality);
- a batch-level final z-score;
- drop groups whose spread is below about 1e-2.

**Execution.** Build per-prompt caches offline (vote, union, core, distance transforms, τ, T5, topology set). Score with a persistent pool of 8 workers, OpenCV pinned to 1 thread each.

## Risks

- F5 assumes Gaussian rater excursions; validate on another multi-rater set.
- K8 may be better defined within a protocol or adjusted for rater count.
- Single-rater items need a size prior, which is unstable.
- The IoU gate gives no topology signal on low-IoU nodules.
- IoU is scored against the vote mask, which **conflicts with t_dynamics** (that causes collapse): resolve in iteration 2.
