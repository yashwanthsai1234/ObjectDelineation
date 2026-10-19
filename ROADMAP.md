# Roadmap: a certified geometry checker as a verifiable reward

Written 2026-10-18, after Phase 0. Revised the same day to match the two research reports.

**Sources**
- `docs/research/reports/GPU delineation use case map.md`: the domain tree, and why the pitch is correctness and topology, not speed.
- `docs/research/reports/GPU delineation RL mix and match.md`: the reward design, prior art, safe claims, the experiment program, and the secondary ideas.
- `experiments/phase0/RESULTS.md`, `experiments/phase0/RESULTS_0b.md`: what Phase 0 measured.

## How this roadmap relates to the research

The research program had four phases:

| Research phase | Content | Status |
|---|---|---|
| 0 | Census and engine validation | done for medical; the SA-1B and natural-image census is not done |
| 1 | Polygon VLM, single ground truth, natural images | — |
| 2 | Multi-rater ambiguity, medical | — |
| 3 | Generative models and connectivity | — |

It also ran a lithography track in parallel.

Phase 0 changed four things:

| Research assumption | Phase 0 finding | Change |
|---|---|---|
| Raters may rarely disagree on topology (the main kill risk) | LIDC readers do (6–12 % of nodules; every nodule with a hole). RIGA's experts never do. | The rater-set claim survives. LIDC carries it; RIGA tests the band only. |
| Polygon VLM can start from a base model | Qwen3-VL-2B gave 1–2 valid polygons of 30 | SFT with hole-aware tokens comes first (the research's A0′ arm). |
| Band = 1 − (outside egg + missing yolk) / \|egg\| | Collapses to ≈ 0 on small nodules | Make it scale-aware. |
| Natural images first, then medical | The headline (rater sets) lives in the medical data we already have | Run medical (Track A) first, then natural images (Track B), so the claim generalises. |

Everything else in the research program is kept below:
- the arms and controls;
- the hacking guards;
- the 3D, video, connectivity and lithography tracks;
- the "do not claim" list.

## Claims

Safe claims, each with a dated "to our knowledge" footnote naming VectorLLM++ and Learning to Trace:

1. Exact per-instance β0 / β1 / χ and hole-aware ring legality **of the model's own output**, used as GRPO reward terms.
2. A boundary tolerance that varies over space and is **derived from annotators** (the egg–yolk band), used as an RL reward.
3. Topology scored against the **set** of rater topologies, not a single ground truth. Phase 0 supplies the evidence that raters disagree.

Defensible only as experimental results: "graded topology beats binary validity (FloorplanVLM) and clDice (Topo-R1)".

Do not claim:
- first topology reward in VLM RL (Topo-R1);
- first polygon-validity reward (FloorplanVLM);
- first connectivity RL (Learning to Trace);
- first geometric verifier (GeoSVG-RL, LithoGRPO);
- that GPU speed makes VLM RL feasible.

## The system

```
 image ─► policy (VLM polygon text │ mask tokens │ mask generator) ─► canonicalize ─► checker ─► reward terms ─► GRPO
                                                                          │
                     reference per item: rater set (Y, E, T) or morphological band from one mask
                                                                          │
                     Level 2 controller ◄── per-term logs, error types, hashes ──► evaluation, deployment guard
```

## The reward (research design, and what exists today)

| Term | Definition | Built? |
|---|---|---|
| format gate | unparseable → bottom of every term | yes (`polygon_io`) |
| r_legal, graded | 1 − (illegal rings + [even-odd fill ≠ nonzero fill]) / (rings + 1) | partly: ring checks exist; fill-rule XOR and graded score to add |
| r_topo | set membership in T, decaying with distance to the nearest rater answer | yes (`topo_freq`, `topo_near`); switch decay to exp(−Δ) as designed |
| r_hole | Hungarian-matched hole rings, ring IoU against consensus holes | **to build** |
| per-instance matching | match components to reference components, not just whole-mask counts | **to build** |
| r_band | pixels between yolk and egg are neutral; penalize outside the egg and inside the yolk | yes; **make scale-aware** |
| r_bnd (optional) | share of the predicted boundary inside the band | to build |
| gate | topology / legality / hole terms count only when band ≥ κ | yes (IoU gate); move to band as designed |
| GDPO | per-term group normalization, sum, batch-normalize | yes (`group_advantages`); add batch normalization |

Rules from the research:
- bands come only from raters, a frozen independent ensemble, or erode/dilate of a single mask, **never from the policy**;
- never repair invalid output (`make_valid`, `buffer(0)`);
- no hard vertex caps (SimpleSeg found they hurt). Use length-bias correction instead.

## Domains and data

The checker does not care what the image shows. Per domain, three things change: the loader, the output format, and where the reference set comes from.

| Track | Domain | Dataset | Reference | What it tests | Status |
|---|---|---|---|---|---|
| **A: multi-rater medical (headline)** | lung CT nodules | LIDC-IDRI (4 readers) | rater set | set-valued topology, 2D and 3D | have it |
| | eye fundus | RIGA (6 experts) | rater set | band (no topology ambiguity) | have it |
| | brain, prostate, kidney… | QUBIQ (2–7 raters) | rater set | multi-organ generality | to get (registration) |
| | MS lesions | MSSEG-2016 (7 raters) | rater set | many small parts (β0) | to get (registration) |
| | skin lesions | ISIC multi-annotator / IMA++ | rater set | band; pairwise Dice is only 0.75–0.79 | to get |
| | pathology | Gleason2019 (6 pathologists), NuCLS (multi-rater nuclei) | rater set | many parts, holes, touching cells | to get |
| **B: natural images (generality)** | referring segmentation | refCOCO/+/g with an enriched split of objects with holes or several parts | morphological band (d = 2 % of the diagonal) | SimpleSeg comparison at scale | to get (public) |
| | everything | SA-1B topology-stratified subset | morphological band | topology census (publishable on its own) | needs the licence link and ~11 GB |
| **C: connectivity and structure** | retinal vessels | DRIVE (2 observers on the test set) | raters / single | component counts (β0), against Learning to Trace's 29.1 → 4.5 | to get |
| | buildings | CrowdAI mapping / SpaceNet (courtyards = holes) | morphological band | holes, multi-part; no multi-rater data exists | to get |
| | floor plans | CubiCasa5k | single | gaps, holes, single outline; the floor-plan RLVR loophole | to get |
| | video | DAVIS / YouTube-VOS | single + temporal | Δβ flicker and ring-hash stability across frames | later |
| **D: GPU speed matters** | lithography | LithoGRPO setup, LithoBench | — | GPU mask-cost proxy replacing a 34–148 s / 0.16–0.30 s CPU step | parallel track |

Access details are from the research notes and must be re-checked before download.

## Phases

| # | Phase | What it does | Where | GPU-hours | Done when |
|---|---|---|---|---|---|
| 0 ✅ | Census and baselines | checker built; experts disagree; MedSAM makes structure errors IoU misses | Mac | 0 | done (`experiments/phase0/RESULTS*.md`) |
| 1 | Reward and data | **Planned in detail: [experiments/phase1/PLAN.md](experiments/phase1/PLAN.md).** C++ polygon core, exact reward v4 (topology-conditional IoU, self-calibrated band, rater-frequency topology term, area guard), SFT target builder, full RIGA/LIDC data and splits, TRL integration, eval harness | Mac | 0 | Phase 1 exit gates in the plan, §7 |
| 1c | Natural-image census | share of objects with holes or several parts in SA-1B / refCOCO | Mac | 0 | numbers for Track B's enriched split |
| 2 | **SFT** | pick the smallest VLM that learns polygons; outer-ring tokens (A0) vs hole-aware tokens (A0′) | Vast | ~6 | ≥ 90 % valid; A0′ recovers holes and parts that A0 cannot |
| 3 | **GRPO, Track A (medical)** | arms A1–A9 plus best-of-8 and DPO, on rater sets | Vast | ~50 | decision on H1–H3 |
| 4 | GRPO, Track B (natural images) | same arms with morphological bands | Vast | ~40 | generality of H1 |
| 5 | Level 2 controller | checker steers training; overfitting and hacking alarms | Vast | ~10 | H4 |
| 6 | Track C | 2 of vessels, buildings, floor plans; then video | Vast | ~20 | H5 |
| 7 | 3D (2×2×2 bit-cubes) | 3D checker; whole-nodule reward on LIDC | Mac → Vast | ~15 | certificate on all volumes; H6 |
| 8 | CUDA checker plus Track D lithography | GPU port of the checker, hash-identical to Python; lithography mask-cost proxy | cheap Vast card | ~20 (cheap) | identical hashes; Spearman ρ ≈ 1 against the reference shot counter |
| 9 | Evaluation and release | held-out sites, deployment guard, data release, paper | both | ~5 | — |
| | **Total** | | | **~150** | **Track A alone (2, 3, 5, 9): ~70** |

GPU-hours assume a 2B model on one A100/H100. A SmolVLM cuts them by roughly 2–4×. See *Compute* and *GPUs* below.

Phases 1c, 7 and 8 do not block training. They run on the Mac, or on a small GPU, while the big runs go.

### Phase 1: reward and data (Mac)

> **Superseded by [experiments/phase1/PLAN.md](experiments/phase1/PLAN.md).** That plan was chosen by online research, sandboxed tests and an independent review. Where this section and the plan differ, the plan wins. Changes from what is written below:
> - the reward uses topology-conditional IoU, the band F5 and a rater-frequency topology term;
> - polygons are rasterized and checked by a C++ core;
> - SFT targets use the potrace polygon stage, not Douglas–Peucker.

**Where the code goes** (see `README.md`):
- new reward terms go in `src/objdel/reward/`;
- dataset loaders move out of the Phase 0 scripts into a new `src/objdel/data/`, with tests;
- scripts and results go in `experiments/phase1/`.

**Reward terms**
- Build the missing terms (graded legality with the fill-rule XOR, per-instance matching, hole matching, boundary-in-band).
- Make the band scale-aware: 1 − (outside egg + missed yolk) / (|egg| + τ · perimeter(egg)), τ = 2 px.
- Add the morphological band for single-mask data.
- Brute-force tests for each term, as for the checker.

**Data**
- **Splits:**
  - RIGA: train MESSIDOR plus BinRushed; test Magrabia, an unseen hospital and MedSAM's worst case.
  - LIDC: split by patient 70 / 10 / 20, with reader-disagreement slices in each split.
- **Output format:** this must be able to express any expert answer.
  ```json
  {"nodule": [{"outer": [[x, y], ...], "holes": [[[x, y], ...]]}]}
  ```
  Coordinates are integers 0–1000 relative to the crop (Qwen-VL convention).
- **SFT targets:**
  - Simplify the exact rings with topology-preserving Douglas–Peucker.
  - Accept a ring only if the checker gives the same (β0, β1) and IoU ≥ 0.98. Otherwise lower ε.
  - No hard cap on vertices.

**Smoke test:** a few LoRA steps on the Mac, end to end.

### Phase 2: SFT

**Model: the smallest VLM that passes the gate.** Short SFT runs (about 1 GPU-hour each) on SmolVLM-256M,
SmolVLM-500M and Qwen3-VL-2B (the smallest Qwen3-VL). Keep the smallest that reaches ≥ 90 % valid polygons
with IoU near MedSAM.

SmolVLM's language model is not Qwen, so it can also serve as the non-Qwen control (Qwen bases gain even under random rewards).

| | |
|---|---|
| Arms | **A0**: outer ring only (SimpleSeg-style). **A0′**: hole- and part-aware tokens. |

- Each epoch's target is a random expert's polygon, so the model learns the spread of answers, not an average.
- Augmentation: flips, 90° rotations, scale jitter.

### Phases 3 and 4: GRPO arms (from the research program)

All arms start from the A0′ checkpoint. Screen with 1 seed, then run 3 seeds for the arms that decide.

| Arm | Reward | Role |
|---|---|---|
| A1 | random | control: Qwen gains under random reward |
| A2 | format only | control |
| A3 | IoU only | today's practice (SimpleSeg recipe) |
| A4 | IoU + binary validity | FloorplanVLM-style baseline |
| A5 | IoU + clDice on the output | Topo-R1-style baseline |
| A6 | IoU + fixed buffer | TraversRL / NSD-τ baseline |
| A7 | all terms, topology against a single ground truth or the majority | ablation |
| **A8** | **all terms, rater-set topology and band, GDPO** | **ours** |
| A9 | A8 with a band derived from the policy itself (conformal) | negative control: should drift lenient |
| + | best-of-8 with the checker on A0′; DPO on checker-ranked pairs | is RL needed at all? |

Track B runs A1–A8 with morphological bands in place of rater sets.

Settings to start from:
- G = 8, temperature 1.0, KL 0.01–0.04;
- Dr. GRPO length-bias correction;
- κ = 0.3 on the band.

Logged throughout:
- each reward term;
- |M| / |E| (the predicted mask's size relative to the egg);
- component-count histogram;
- duplicate-hash rate within each group;
- pass@k against the base model.

Hypotheses:
- **H1**: A8 cuts the β error at equal IoU (within 1 point) against A3 and A4, by more than seed variance against A1 and A2.
- **H2**: on LIDC slices where readers disagree, A8 hits a rater's topology more often than A7 without losing NSD. Its sampled topologies match the readers' distribution.
- **H3**: A8 beats best-of-8 and DPO at equal compute; A9 shows leniency drift.

### Phase 5: Level 2 controller

The controller turns the checker's typed errors into training decisions.

**Signals:**
- error-type rates (extra/missing parts, extra/missing holes, illegal rings);
- per-site validation curves.

**Actions:**
- curriculum on stuck error types;
- raise KL when validation falls while training reward rises;
- early stop.

**Overfitting checks:**
- the training–validation gap;
- per-site curves;
- memorization: an output's ring hash equals a training target's hash.

**Hacking alarms** (from the research):
- dilating until parts merge;
- specks or pinholes added to change counts;
- filling every hole;
- hugging the yolk;
- self-intersecting polygons that happen to fill well;
- vertex padding;
- IoU rising while structure errors also rise.

**H4**: A8 plus the controller reaches A8's held-out score in fewer steps, with no train/validation divergence.

### Phase 6: Track C

Pick 2 of:
- vessels: β0, Learning to Trace's benchmark;
- buildings: holes;
- floor plans: gaps and the single outline.

Then video, with temporal topology. This reward can be gamed by masks that never change, so pair it with ground-truth IoU and a flow warp.

**H5**: the checker's terms improve structure in each domain at equal overlap.

### Phase 7: 3D (2×2×2 bit-cubes)

| 2D (done) | 3D |
|---|---|
| β0, β1 | β0, tunnels β1, cavities β2 |
| χ from 16 bit-quad patterns | χ from 256 bit-cube patterns |
| rings | closed voxel-face surface |
| area = pixel count = ring area | volume = voxel count = mesh volume |

**Certificate.** Two independent routes to χ must agree:
- χ from the bit-cube histogram equals half the surface Euler number (V − E + F);
- β0 and β2 from 26/6 labeling give β1 = β0 + β2 − χ;
- the mesh volume equals the voxel count.

**3D reward.** The VLM outlines every slice of a nodule; the slices are stacked; the volume is scored against the readers' 3D answers.

**H6**: fewer inconsistencies between slices than per-slice training gives.

Further 3D use (research Phase 3): Flow-GRPO on an ambiguous LIDC mask generator, measured by GED.

### Phase 8: CUDA checker and lithography

- Port the array-only design (bit-quads, edges, successor, pointer jumping, list ranking) to CUDA, building on the original code.
- Then add the fill-rule polygon rasterizer: even-odd and nonzero in one pass.
- **Equivalence**: identical ring hashes on all 67,009 Phase 0 masks.
- **Benchmarks**: masks or volumes per second, against Python and OpenCV.
- **Lithography**: a GPU mask-cost proxy (vertices, jogs, holes, islands, rectangle-partition bound).
  - It must rank masks like LithoGRPO's shot counter (Spearman ρ ≈ 1).
  - This is the one place where speed is the product.

### Phase 9: evaluation, deployment, release

**Test sets:** held-out hospitals and patients.

**Metrics:**
- IoU, Dice, NSD, HD95;
- β error and the topology-wrong rate at IoU ≥ 0.8;
- rater-set hit rate, legality, hole recall.

**Clinical measures:**
- cup-to-disc ratio against the experts' range;
- nodule diameter and volume against the readers' range.

**Deployment guard:** certify every output; resample K times; flag for a human. Report accuracy against flag rate.

**Releases:**
- certified polygon versions of the datasets, with each item's topology set and band;
- the SA-1B topology census;
- the checker as a cross-domain verifier library.

**Paper:** H1–H6, negative results included.

## Compute (estimate, Vast AI)

GPUs are rented on Vast AI (decided 2026-10-18).
- Training: on-demand A100 / H100 80 GB, checkpointed so an interrupted instance loses little.
- CUDA development: the cheapest NVIDIA card available.
- Data: preprocessed into small crops on the Mac and uploaded; well under 1 GB, not the 25 GB raw download.

**Per-run estimate, for a 2B model** (the largest candidate) on one A100/H100:
- One GRPO run is about 300 steps of 32 prompts × 8 rollouts, each rollout about 300 output tokens.
- That is about 25 s per step, or **2–3 GPU-hours per run**.
- SFT is about 1 GPU-hour.
- A smaller model (SmolVLM) cuts these by roughly 2–4×.

The first Vast session measures the real step time and replaces these numbers.

| Phase | GPU-hours |
|---|---|
| 2 SFT: model screening (3 models) + A0 vs A0′ | ~6 |
| 3 Track A GRPO: 10 arms × 1 seed, then 3 decisive arms × 2 more seeds, plus best-of-8 and DPO | ~50 |
| 4 Track B: same pattern on a capped natural-image subset | ~40 |
| 5 Level 2 controller | ~10 |
| 6 Track C, 2 domains | ~20 |
| 7 3D | ~15 |
| 9 Evaluation | ~5 |
| **Total** | **~150 (A100/H100)** |
| 8 CUDA development | ~20 on a cheap card |

**Track A alone** (Phases 2, 3, 5, 9) is about 70 GPU-hours.

## GPUs

The GPU is chosen when we reach Phase 2, from the first session's measured memory use and speed. The tables are the starting estimate.

**Memory per model.**
- Full fine-tuning needs about 16 bytes per parameter (weights, gradients, optimizer state).
- GRPO also holds a frozen reference copy of the model and space to generate the rollouts.

| Model | SFT | GRPO | Fits on |
|---|---|---|---|
| SmolVLM-256M | ~6 GB | ~10 GB | any 16–24 GB card (RTX 3090 / 4090) |
| SmolVLM-500M | ~10 GB | ~16–20 GB | RTX 4090 / 3090 (24 GB) |
| Qwen3-VL-2B | ~35–40 GB | ~50–60 GB | A100 80 GB or H100 80 GB |
| Qwen3-VL-2B, LoRA | ~12 GB | ~20 GB | RTX 4090 (24 GB); may learn the polygon format less well |

**Vast AI cards.** Prices are rough; the marketplace changes constantly.

| Card | Memory | Rough $/hour | Use |
|---|---|---|---|
| RTX 3090 | 24 GB | 0.15–0.25 | CUDA development; SmolVLM SFT |
| RTX 4090 | 24 GB | 0.30–0.50 | Phase 2 screening; all SmolVLM training |
| A6000 / L40S | 48 GB | 0.40–0.90 | Qwen 2B with memory-saving tricks |
| A100 80 GB | 80 GB | 0.80–1.50 | Qwen 2B full GRPO |
| H100 80 GB | 80 GB | 1.60–2.50 | same as A100, about 2× faster |

**Starting plan.**
- Screen on one RTX 4090 (about $2–5).
- If a SmolVLM passes, Track A RL stays on the 4090 (about $15–40).
- If Qwen 2B is needed, Track A RL moves to one A100 80 GB (about $60–100).

**When renting**, ask for:
- at least 8 CPU cores (the checker scores 256 rollouts per step on the CPU);
- at least 50 GB of disk;
- on-demand, not interruptible, instances for training;
- frequent checkpoints.

## Kill and continue criteria

- **After Phase 2**: if A0′ cannot reach 90 % valid polygons:
  - use fewer coordinate tokens or a coarser grid;
  - try a 4B model;
  - fall back to mask tokens.
- **After Phase 3**: if A8 does not beat A3 on β error at equal IoU, the claim fails. IoU already captures the signal, and we report it.
- If A8 does not beat A1 and A2 beyond seed variance, the gain is the Qwen random-reward effect, not the checker.

## Risks

| Risk | Mitigation |
|---|---|
| Unread prior art (VectorLLM++, Learning to Trace) already uses exact β rewards | re-check before writing; claims carry dated footnotes |
| IoU captures most of the signal | A3 vs A8 decides it; a null result is reported |
| Topology errors too rare for RL to see | oversample high-IoU-wrong and disagreement items |
| Multi-rater data is small and medical | Track B (natural images, morphological bands) covers generality |
| Policy emits prompts and SAM draws the mask, so the reward grades SAM | main arms use policies that emit geometry; SAM pipeline only as Plan B |
| 16 GB Mac | Mac does data, checker, 3D and analysis; training runs on rented GPUs |
