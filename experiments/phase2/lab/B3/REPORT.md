# B3: probes and offline augmentation (Phase 2 planning, wave B, CPU)

Pre-registered criteria: `prereg.py` (written before any run; scripts assert against it). One re-registration, marked
RE-REG: the precision-probe *generator* was redesigned after run 1 failed calibration; the criteria did not change.
(Report text returned by the B3 agent and saved by the coordinator; the harness blocked the agent's write.)

## 1. Verdicts

| Test | Verdict |
|---|---|
| H-A1-7: topology probe of ~300 held-out real LIDC items | **Rejected as stated.** No real-data design reaches ±6 points at the measured within-nodule correlation, and the ~300-item version removes every non-simple train target from SFT. **A synthetic topology probe passes:** 320 non-simple items, ±5.5 points |
| Unanimous non-simple items | val 0; train 12 items (7 nodules); test 4 items (1 nodule) |
| Precision probe (A4 design, H-A4-7, H-A4-8) | **Pass with generator 2** (generator 1 failed calibration). Exact on perfect answers, sensitive to 1-px shifts, separates every synthetic model as predicted. **But the centroid disc passes every A4 geometry criterion; only topology fails it** |
| H-A6-6 / H-A4-5: offline dihedral targets | **Pass** on certification, reproducibility and CPU time; two misses as written (wall time under load; one non-SFT target over budget) |
| H-A6-9: crop before trace | **Pass.** Target-set hash unchanged (`1af819f4…`); RIGA builds 23 % faster per target |

## 2. What the plan should adopt

1. **A synthetic topology probe decides A0 vs A0′.** 160 one-hole, 160 two-part and 160 simple controls, independent
   items rendered on real LIDC val lung patches. Metrics: non-simple hit (primary), hole recall/precision, part
   recall/precision, false non-simple rate on simple items. A0 scores exactly 0.50 [0.45, 0.55] (all two-part, no
   hole). Keep every train non-simple item in SFT (H-A1-6 needs them for oversampling). Report the val non-simple
   slice (44 items) descriptively with its ±15–23-point cluster CI. Replace the entry gate "topology-any ≥ 90 % on
   unanimous items": val has none, so it is vacuous.
2. **Precision probe rule adds a paired win over the centroid disc:** per-item IoU beats the equal-area disc at the
   true centroid on ≥ 60 % of 6–16 px items (exploratory, not pre-registered: σ ≤ 0.5 px models win 78–100 %;
   σ ≥ 0.75, the 1-px shift and the disc ≤ 49 %). Topology is scored on the topology probe (the precision probe's
   5 % hole/two-part share leaves only 62 items). The legality ≥ 95 % cut is knife-edge (σ = 0.5: 92.3 % on
   generator 1, 96.4 % on generator 2).
3. **Precompute the 8 dihedral variants offline with crop-before-trace.** Store one frame-resolution crop per item
   plus 8 target texts; transform the image at collate time at frame resolution, before the LIDC upscale (A4).
   All 8 crops would cost ~1.2 GB; one crop + 8 text sets ~220 MB. Keep a per-variant `target_set_hash`; apply the
   token-budget filter per variant. Promote the `_rings` crop into `src/` (with tests) in Phase 2. LIDC variants
   are pixel-exact and length-identical.

## 3. Results

### P1. Topology probe (H-A1-7): `census.py`, `topo_probe.py`

Census of every LIDC item under the reward's size rule and voting raters (35 s, 0.8 GB):

| Split | Non-simple items | Nodules | Patients | Unanimous | Majority non-simple |
|---|---|---|---|---|---|
| train | 301 | 100 | 89 | 12 (7 nodules) | 47 |
| val | 44 | 14 | 12 | 0 | 3 |
| test | 83 | 22 | 21 | 4 (1 nodule) | 12 |

A1's raw 366 train items drop to 301 under the size rule; one val patient holds 19 of the 44 val items. Adjacent
slices are correlated: intraclass correlation of the raters' non-simple share within nodules 0.15–0.19.

Cluster-bootstrap CI half-widths (resampling nodules) at a true rate of 0.5:

| Design | Items / nodules | ρ = 0 | ρ = 0.15 | ρ = 0.3 | ρ = 0.6 | Train non-simple targets kept | Train items removed |
|---|---|---|---|---|---|---|---|
| val | 44 / 14 | ±15.0 | ±17.0 | ±19.5 | ±23.2 | 100 % | 0 |
| val + half the train non-simple patients | 195 / 63 | ±6.8 | ±9.0 | ±10.5 | ±13.1 | 45 % | 997 |
| val + all train non-simple patients | 345 / 114 | ±5.2 | ±6.8 | ±8.0 | ±10.2 | 0 % | 2,058 |
| synthetic, independent items | 320 | ±5.5 | ±5.5 | ±5.5 | ±5.5 | 100 % | 0 |

Criterion 1 (±6) fails for every real design at the measured ICC; criterion 3 fails for every train hold-out. The
paired A0 vs A0′ difference is ±8.5 to ±21 on real designs. The synthetic probe passes at n ≥ 270.

### P2. Precision probe: `probe_lib.py`, `precision_probe.py`, `probe_rule.py`

Conventions: boundary = 4-connected inner boundary pixels (MONAI, A4 S2); ASSD pooled over both directions; HD95 =
max of the two directed 95th percentiles (implementations differ, Podobnik & Vrtovec 2024, so stated with every
number).

- P2.1 exactness: PASS. On 1,200 items the ground-truth raster scores ASSD 0, HD95 0, NSD 1; certified target texts
  are 100 % pixel-exact.
- P2.2 correctness: PASS. Max difference 0.0 vs brute force over 3,600 pairs.
- P2.3 1-px shifts: PASS. ASSD > 0 on 100 % of 4,800 shifts, centroid error exactly 1 px; median ASSD 0.70, IoU
  0.80; the pass rule fails the shift model.
- P2.4 calibration: run 1 (one amplitude parameter) FAILED (disc IoU 0.771 vs 0.806 ± 0.03; σ = 1 legality 69 % vs
  84 %; needed 16 vertices vs 11 real). RE-REG (generator only): run 2 adds an elongation term E; (E, A) calibrated on
  seed 100 to disc IoU and vertex count, legality left out as the check; E = 0.2, A = 0.15. Fresh-seed confirmation
  PASSES:

  | Check | Synthetic | Real LIDC | Tolerance |
  |---|---|---|---|
  | centroid-disc IoU | 0.822 | 0.806 (train sample 0.804) | ± 0.03 |
  | σ = 0.5: IoU / legal | 0.915 / 96.4 % | 0.90 / 95.7 % | ± 0.05 / ± 5 pts |
  | σ = 1: IoU / legal | 0.816 / 85.4 % | 0.79 / 83.7 % | ± 0.05 / ± 5 pts |

  Diameters p10/p50/p90 6.8 / 12.2 / 25.8 px vs real 6.9 / 12.2 / 26.6. Half-way threshold of the noise-free
  composite recovers the ground truth at median IoU 0.992 (p10 0.85); needed lung-only background patches.
- P2.5 discrimination: PASS. Pass: exact, σ = 0.25, σ = 0.5. Fail: σ = 1, σ = 2, centroid disc, centred prior
  (IoU 0.00), 1-px shift. Caveat: the centroid disc (legal 100 %, ASSD 0.53, IoU 0.82, centroid error 0.08 px)
  fails only on topology.
- Synthetic topology probe: exact answers 1.0 everywhere; A0 non-simple hit 0.50 [0.45, 0.55] (hole recall 0, part
  recall 1); disc 0 on every recall; σ = 0.5 non-simple hit 0.92 [0.89, 0.95], no false non-simple on simple items.

### P3. Offline dihedral targets: `dihedral.py`

Every train target (LIDC 23,469; RIGA 3,552) under the 8 dihedral ops:

- P3.1 certification: PASS. 215,208 variant targets build and re-check from text alone, 0 errors.
- P3.2 inverse transform: PASS. 100 % keep IoU ≥ 0.99 and the signature vs the original rater mask. Against the
  original reward caches: LIDC F5 and IoU_tc min 1.0; RIGA F5 min 0.975, IoU_tc min 0.99. LIDC variants are
  pixel-identical to the identity; RIGA is not (vflip 98 %, hflip 0 %, rotations ~0 % pixel-identical: top-left
  tie rule, I7) but every variant certifies.
- P3.3 identity hashes: PASS (LIDC `b4814506…`, RIGA `54bc88bf…`).
- P3.4 reproducibility: PASS. Two runs (2 processes with checks; 1 process) give identical hashes and byte-identical
  files; overall `e86e07e6…`; per-op hashes in `dihedral_run1.json`.
- P3.5 cost: CPU PASS (16.5 and 18.1 build CPU-minutes); wall time FAILS as measured (22 and 33 min at load average
  15–42 on 8 cores; the pool's parent also reads the store, as in Phase 1 `build_targets`).
- P3.6 token budget: median relative length change 0 for every op and tokenizer (SmolVLM-256M, LFM2.5-VL-450M,
  Qwen3-VL-2B, Qwen3.5-0.8B); RIGA p5–p95 −2.8 % to +5.8 %, LIDC exactly 0. Fails as written on the full set: one
  BinRushed1 target (image044, rater 4) grows from 1,052 to 2,196 SmolVLM tokens under rot270 / antitranspose (not
  in SFT). On SFT rows it holds; the only over-budget SFT target is BinRushed3/image014 (I5), over budget under every
  op.
- Disk: frame-resolution PNG crops LIDC 17.6 MB, RIGA 129 MB for the identity (~1.17 GB for all 8 ops); target
  texts ~57 MB LIDC + 15 MB RIGA for 8 ops.

### P4. Crop before trace: `crop_trace.py`

Monkeypatches `_rings` in its own process only (`src/` untouched); each target built both ways, alternating order.
Every per-file hash and the overall hash `1af819f4…` unchanged, 0 texts differ. Median paired time ratio: RIGA train
0.768 (p10–p90 0.59–1.03), RIGA val 0.772, LIDC 0.925. RIGA train total 106 s → 82 s.

## Sources

- MONAI surface distance: https://github.com/Project-MONAI/MONAI/blob/dev/monai/metrics/surface_distance.py ; edges:
  https://github.com/Project-MONAI/MONAI/blob/dev/monai/metrics/utils.py
- DeepMind surface-distance: https://github.com/google-deepmind/surface-distance/blob/master/surface_distance/metrics.py
- Metrics Reloaded: https://arxiv.org/abs/2206.01653
- Podobnik & Vrtovec, distance-metric implementation pitfalls: https://arxiv.org/abs/2410.02630
- Yeghiazaryan & Voiculescu 2018, boundary metrics: https://www.cs.ox.ac.uk/publications/publication11655-abstract.html
- Nikolov et al., NSD: https://arxiv.org/abs/1809.04430
- Field & Welsh 2007, cluster bootstrap: https://ideas.repec.org/a/bla/jorssb/v69y2007i3p369-390.html
- Kish design effect (Gabler et al. 1999): https://www150.statcan.gc.ca/n1/pub/12-001-x/1999001/article/4718-eng.pdf
- Stucki et al. 2023, Betti matching: https://proceedings.mlr.press/v202/stucki23a
- Qwen technical report (single-digit numbers): https://arxiv.org/abs/2309.16609
- Internal: `docs/research/notes/Phase 2 SFT/A1`, `A4`, `A6`; `experiments/phase1/RESULTS.md` §3, §8.

## Notes

`topo_probe_items.json` lists the rejected half-hold-out design (kept for reference). Large outputs (variant targets,
example images) are in the session scratchpad `phase2/B3/`. Nothing in `src/`, `tests/`, `experiments/phase1/`
changed.

---

## Addendum A: hole-calibrated topology probe and the C2 primary metric (judge MUST-FIX 1, 3)

Pre-registered as P5 in `prereg.py` before any of these scripts ran. Scripts: `hole_census.py`, `topo_probe_v3.py`,
`power_c2.py`. Outputs: `hole_census.json`, `topo_probe_v3{a,b,c}.json`, `topo_probe_v3d_seed{5,6}.json`,
`topo_probe_v3e_seed{7,8}.json`, `power_c2.json`. Raw per-target records: scratch `phase2/B3/hole_records.jsonl`.

### A1. Real hole and part geometry (P5.1, descriptive)
Scope: voting raters of every non-simple LIDC item, frame resolution, after the reward's size rule.

| | train | val | test |
|---|---|---|---|
| holed targets (1 / 2 / ≥ 3 holes) | 241 (204 / 32 / 5) | 44 (38 / 6 / 0) | 62 (49 / 7 / 6) |
| holes | 285 | 50 | 83 |
| hole area px, p10 / p50 / p90 | 5 / 12 / 49 | 5 / 10 / 26 | 5 / 13 / 41 |
| bins ≤ 10 / 11–30 / > 30 px | 44 / 36 / 20 % | 54 / 36 / 10 % | 39 / 45 / 17 % |
| hole as a fraction of its filled part, p50 (p90) | 3.6 % (15 %) | 3.0 % | 2.3 % |
| centroid offset / equivalent radius, p50 | 0.46 | — | — |
| wall to the outside, px, p10 / p50 / p90 | 0.41 / 1.24 / 4.0 | 0.41 / 2.0 / 4.1 | 0.41 / 1.24 / 3.9 |
| diameter of the holed part, px, p50 | 23.6 (all splits) | | |
| holed targets whose item has only 1 / 2 / 3 hole-drawing raters | 182 / 56 / 3 | 36 / 8 / 0 | 48 / 14 / 0 |
| multi-part targets; part-area ratio p50; gap px p10 / p50 / p90 | 182; 0.57; 1.2 / 6.8 / 31.5 | 5 | 44 |

Readings:
- Real LIDC holes are small: median 12 px, 44 % at or below 10 px, 3.6 % of the object.
- They sit behind thin walls (median 1.2 px; p10 0.41 px, a hole diagonally touching the outside), in large objects (about 24 px).
- They are idiosyncratic: 75 % of holed targets come from an item where only one rater drew a hole, and no hole is unanimous.
- Tertiles of hole area: 8 / 19 px. I used the judge's fixed bins (≤ 10, 11–30, > 30 px).

### A2. Generator recalibration and fidelity (P5.2)
- Outlines are generator 2 (E = 0.2, A = 0.15).
- Each hole resamples one train hole record jointly: the part's diameter, hole area, centroid offset and wall.
- Two-part items resample one train multi-part record jointly: object diameter, part-area ratio and gap.
- Every hole is at least the size rule: the ground truth must have signature (1,) both raw and at `min_area(|GT|)`.

The criterion: KS p ≥ 0.01 on six geometry variables, and bin shares within ±5 points, on a natural draw of 300 holes and 200 two-part items.

| Version | Hole placement | Seed | Lowest KS p (variable) | Bin shares | Pass |
|---|---|---|---|---|---|
| v3a | by centroid offset | 5 | 0.0000 (wall; area 0.0001) | 44 / 37 / 19 | no |
| v3b | by wall | 5 | 0.0000 (offset) | 47 / 34 / 19 | no |
| v3c | offset + outline dented to the wall, wall tolerance 0.8 px | 5 | 0.006 (wall) | 43 / 39 / 18 | no |
| v3d | wall tolerance 0.45 px | 5 (tuning seed) | 0.020 (offset) | 45 / 36 / 19 | yes |
| v3d | same | 6 (fresh) | 0.001 (area) | 42 / 39 / 20 | no |
| v3e | area tolerance max(1 px, 10 %) | 7 (fresh) | 0.004 (area; offset 0.006) | 40 / 35 / 25 | no |
| v3e | same | 8 (fresh) | 0.022 (offset) | 43 / 36 / 21 | yes |

Verdict: **P5.2 fails as pre-registered.** One of the two fresh-seed confirmations of the final version fails, and I stopped iterating so as not to tune to the test.

What is left is small:
- Probe per-bin median hole areas are 6 / 17 / 48 px against real 5 / 17 / 49. Exact 5-px holes, at the size-rule floor, are under-drawn by about 1 px.
- The low-offset tail is slightly less central (offset p10 0.25 vs 0.18).
- Wall, fraction, part ratio and gap match on every seed (p ≥ 0.11).

Because the probe is stratified by hole-area bin, the area mismatch only acts inside a bin. I judge the probe fit for C2, with this caveat stated. The plan should record P5.2 as "failed narrowly (1 px area bias in the smallest bin)", not as a pass.

### A3. Probe design and power (P5.3)
Probe composition:
- 120 one-hole items per bin, 360 in total;
- 160 two-part items, descriptive only (A0 keeps parts);
- 160 simple controls.

Pooled hole recall is the train-share-weighted mean of the bins (0.439 / 0.361 / 0.200), with a stratified bootstrap. `power_c2.py` simulates it with item difficulty SD σ_u (logit scale), a seed effect σ_v = 0.3, an item×seed term of 0.5, and 2 seeds:

| CI half-width at true recall 0.5 | greedy | mean over K = 8 | any of 8 |
|---|---|---|---|
| pooled, σ_u = 0.5 / 1.5 / 3 | ±4.1 / ±4.9 / ±5.1 | **±1.9 / ±3.1 / ±4.0** | ±3.9 / ±4.4 / ±4.8 |
| per bin (n = 120), σ_u = 1.5 | ±8.1 | ±5.1 | ±7.3 |

- **P5.3 passes:** pooled ±3.1 ≤ ±6 under the recommended metric.
- **The pooled test is the powered one.** The per-bin floor "hole recall ≥ 0.30" is a check on point estimates. It is reliable at n = 120 (mean-over-8, σ_u = 1.5):

  | True bin recall | P(point estimate < 0.30) |
  |---|---|
  | 0.50 | 0.000 |
  | 0.40 | 0.010 |
  | 0.20 | 0.998 |
  | 0.15 | 1.000 |

  So the floor rarely fails a bin whose true recall is ≥ 0.4, and almost always catches one at ≤ 0.2.

### A4. Sanity on the recalibrated probe (P5.4; v3e seed 8, seed 7 alike)
| Model | Hole recall, ≤ 10 / 11–30 / > 30 px | Pooled (95 % CI) | Part recall (descriptive) | False non-simple on controls |
|---|---|---|---|---|
| exact | 1 / 1 / 1 | 1.0 | 1.0 | 0 |
| A0 (holes emptied) | **0 / 0 / 0** | 0 | 1.0 | 0 |
| centroid disc | 0 / 0 / 0 | 0 | 0 | 0 |
| σ = 0.5 px, mean over 8 | 0.65 / 0.81 / 0.90 | 0.76 [0.73, 0.78] | 0.95 | 0 |
| σ = 0.5 px, single draw | 0.60 / 0.88 / 0.92 | 0.76 [0.72, 0.80] | — | — |
| σ = 0.5 px, any of 8 | 0.98 / 0.99 / 1.00 | 0.99 | — | — |

- All three expectations hold: exact scores 1, A0 scores 0 in every bin by construction, and the disc scores 0.
- Recall rises with hole size: a 0.5 px vertex error destroys about a third of the ≤ 10 px holes.
- Any-of-8 saturates at about 0.99 and hides that loss.

### A5. Primary metric and seeds for C2 (P5.5)
**Primary metric: mean over K = 8 samples at T = 1, per item, averaged over the 2 seeds.** Greedy and any-of-8 are reported alongside.

Why mean-over-8:
- It is the tightest: ±3.1 vs ±4.9 for greedy, so greedy would need about 2.5× the items for the same CI.
- It estimates the probability the policy draws the hole, which is what GRPO samples from.
- Any-of-8 is excluded by the pre-registered rule. It saturates (0.99 for σ = 0.5) and rewards a policy that stamps holes at random. Controls then catch that only through the false-hole rate.

Seeds: the bootstrap resamples probe items, stratified by bin, after averaging each item over seeds and samples. On the real val slice it resamples nodules instead.

**Caveat (new):** this CI is conditional on the two trained checkpoints. It does not cover seed-to-seed variation of the training recipe. Its coverage of the seed-population mean is 95 % with no seed effect, but only **54 % at σ_v = 0.3 logit** and **30 % at 0.6**. Two seeds cannot estimate that variance. So report:
- the pooled CI;
- per-seed CIs;
- a seed-discordant flag (per-seed CIs on opposite sides of the threshold);
- and take the decision with both seeds' point estimates required to meet it.

### A6. What C2 should say (replaces the hole/part wording)
A0 is A0-multi, so on hole items the paired difference (A0′ − A0) equals A0′'s hole recall; it adds no information and is reported only for completeness.

C2 passes iff all of these hold:
1. A0′ pooled hole recall (mean over 8, 2 seeds) ≥ 0.50, with 95 % CI lower bound > 0.30;
2. both seeds' pooled point estimates ≥ 0.50;
3. no bin's point estimate < 0.30;
4. false non-simple on controls ≤ 5 %;
5. legal ≥ 90 %.

C2 fails, with "holes left to GRPO", if no arm reaches pooled 0.30. Part recall is descriptive.

Also report the 44-item real val slice descriptively. Its holes are smaller (54 % ≤ 10 px), and most are single-rater; the matching real-data quantity is the K-sample hole rate against the rater hole fraction (TV), not recall.

Additional sources:
- Kish design effect and the limits of cluster sampling: https://www150.statcan.gc.ca/n1/pub/12-001-x/1999001/article/4718-eng.pdf
- Cluster bootstrap: https://ideas.repec.org/a/bla/jorssb/v69y2007i3p369-390.html
- pass@k as an optimistic, hedging-tolerant metric: Chen et al. 2021, https://arxiv.org/abs/2107.03374
