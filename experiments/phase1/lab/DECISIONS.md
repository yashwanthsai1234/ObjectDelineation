# Phase 1 decision record (draft for the judge)

**Goal of Phase 1:** build, test and freeze everything the GPU phases need, on the Mac:
- the exact reward (C++ core plus Python orchestration);
- SFT targets;
- data, splits and storage;
- TRL integration and the Level 2 controller hooks;
- KPIs and gates.

**Evidence.** Research notes are in `docs/research/notes/Phase 1 stack and algorithms/*.md` (6 notes). Lab reports are in `phase1_lab/t_*/REPORT.md` and `t_reward/results/ITERATION2.md`.

## Module decisions

| # | Module | Decision | Rejected alternatives (why) | Evidence (lab) | KPI result |
|---|---|---|---|---|---|
| D1 | Polygon → mask | C++17 header-only `polycore` integer scanline, winding per pixel: even-odd + nonzero in one pass, top-left tie rule = oracle spec. Per-ring mixed-sign flag. nanobind, GIL released, internal thread pool. | PIL/OpenCV (do not sample pixel centres; 135–175 wrong px per polygon; PIL inflates small nodules by 1.56×). skimage/matplotlib/shapely/rasterio (even-odd only, undocumented ties, 4–45 ms). numba (5× slower, JIT warm-up). | t_polycore | K1 0 / 849.8 M px wrong; K3 8.5 ms per batch of 256; K4 p99 0.38 ms; Linux arm64 + amd64 Docker pass |
| D2 | Legality | Same C++ core: one integer pass over edge pairs (x-sweep), graded per-class counts, never repair. Reward uses presence or saturating counts. | shapely/GEOS (first error only, ignores winding, accepts pinches). Pure Python (43 ms per ring at 200 vertices, ~33 s per step). | t_polycore | K2 0 / 10,000 rings wrong |
| D3 | Mask topology | cv2 Spaghetti labeling (8-connected foreground / 4-connected background) on the bbox crop. Holes assigned per component by the pixel above. Bit-quad χ kept as the certificate. | scipy.ndimage (5–24× slower in our use). Persistent homology (same answer, far costlier). | t_reward | 0 mismatches on 10,000 random + 8,738 fixture masks; 5–24× faster |
| D4 | Band | F5 = 1 / (1 + (2D/P)/T5²): distance-weighted error beyond the egg–yolk band, normalized by the union's contour length. T5 is self-calibrated per image from the raters' own leave-one-rater-out scores, with a rater-count correction. | Current band (ρ with area +0.53). Fixed-pixel τ (H5 rejected: rater disagreement grows with size). √area tolerance (over-tolerant). | t_reward | Median 0.973; rater-count-adjusted ρ −0.07 to −0.09 (pooled +0.107) |
| D5 | Topology reward | Set-valued validity plus distribution matching: r_topo = 1[valid]·(1 + 0.5(log f − log π̂_loo)), mean-centred, c_topo = 2, G ≥ 8. | Plain any-rater (collapses: TV 0.30–0.45). Coverage bonus (spreads answers uniformly over the rater set). KL to SFT (fragile). Majority (collapses). | t_dynamics, t_reward iter. 2 | K9 98.93 %; minority answer no longer disadvantaged (0/14 items, was 10/14) |
| D6 | IoU term | IoU_tc = max IoU over raters sharing the sample's topology (all raters if none). | IoU vs the strict-majority vote (+0.13–0.15 gap between topologies, which drives collapse). Same-topology consensus (unstable). | t_reward iter. 2 | Gap +0.007 [−0.014, +0.026] |
| D7 | Advantages | GDPO with σ floors (IoU/band 0.05, legality 0.25), batch-level z-score, legality gate, topology gated on IoU ≥ 0.5. Topology term added mean-centred after normalization. Drop groups whose spread is below ~1e-2. Implemented in TRL as weight-0 per-term reward functions (logging only) plus one aggregator reward, with `scale_rewards="none"`. | Plain GDPO / TRL `normalize_then_sum` ("dark room": 1e-3 noise shifts advantages by up to 1.24; TRL amplifies a 0.001 legality difference). | t_reward, t_stack | Noise shifts advantages by ≤ 0.006 |
| D8 | RL loss | TRL 1.16 GRPOTrainer: `loss_type="dr_grpo"`, clip 0.2 / 0.28, β = 1e-3 (non-zero at init), `mask_truncated_completions`, constant learning rate. No vLLM for SmolVLM until PR #6718 merges. | verl/EasyR1 (no SmolVLM support; heavier). LLaMA-Factory (no GRPO). RLOO (most prone to hacking). GSPO kept only as an ablation. | t_stack, research | 2–4 GRPO steps run on 4 models; per-term logs present |
| D9 | Level 2 controller | Learnability sampler p(1−p) or group reward variance; typed held-out monitoring; alarm "train reward ↑ while a held-out typed error rises"; rollback to the best checkpoint; the rising type names the guard to switch on. Implemented as a TrainerCallback mutating `beta` and `reward_weights` (verified live). | Typed curriculum sampler ("upweight stalled types": never better). Band-pass difficulty filter (worse than uniform in 10/12 scenarios). | t_dynamics, t_stack | Alarm detects 98.6–100 % of fast hacks, 0 % false alarms in 9/10 scenarios |
| D10 | Output format | Compact JSON, integers 0–1000 relative to the crop, outer rings clockwise on screen, holes counter-clockwise, start at the top-most then left-most vertex, parts sorted. RIGA emits disc + cup and the rim is computed. LIDC emits `{"nodule":[{"outer","holes"}]}`. No vertex cap; `max_new_tokens` 1,536 (Qwen/SmolVLM) or 1,024 (LFM). | Spaced JSON (+22 % tokens). LFM-native 0–1 floats (2×). Special location tokens (data-hungry). Fixed-N vertices (Groundbench collapse). | output_format, t_stack, t_targets | — |
| D11 | SFT targets | `minseg_hybrid` (potrace polygon stage OR lossless digital straight segments), τ = 3 grid units, IoU ≥ 0.99, defect-localized repair, pinch bevel, exact re-certification. Crop side ≤ 500 px, preferring sides that divide 1000. Each rater is a separate target, with a random rater per epoch. | Douglas–Peucker / Visvalingam–Whyatt (more vertices, worse IoU, useless on small crops). Classic potrace ½ px (4–5× the tokens). | t_targets | K5 100 %; K6 min 0.990; K7 median 157 / 197 / 92 tokens |
| D12 | Model shortlist | Phase 2 screening, smallest first: SmolVLM-256M (longest_edge 512 → 64 image tokens), LFM2.5-VL-450M (37 % fewer tokens; licence to check), Qwen3-VL-2B (only one with a real prior: zero-shot box IoU 0.76). Pick the smallest that passes the SFT gate. | SmolVLM-500M (adds nothing over 256M). Florence-2 (no GRPO support). Qwen3.5-0.8B (stretch). | t_stack | SFT overfit valid 4/4 on 256M and 450M |
| D13 | Data | **RIGA:** literature split (train BinRushed + MESSIDOR 655, test Magrabia 95); fix the extractor (BinRushed1, missing raters). **LIDC:** patient-level 70/10/20, patients with two scans kept together, parity check against pylidc. Storage: Parquet / HF datasets with a SHA-256 manifest; exporters to the TRL chat format. | DVC/webdataset/Lance (unneeded at our scale). albumentations (archived; successor is AGPL). | t_data, data_engineering | RIGA 4,500/4,500 masks (was 4,172); 0 cross-split duplicates (37 within-split pairs grouped; Magrabia = 88 unique photos); LIDC parity with pylidc 67/67 scans, 99.57 % masks exact; deterministic group splits 11/11 tests; K11 PASS |
| D14 | C++ vs Python | C++: polycore (raster, legality; later the simplifier and the CUDA port share the header). Python + OpenCV/numpy: topology, band, caches, orchestration. numba: target-builder kernels (offline), to be ported to C++ only if needed. | All-Python (≈1,000× slower legality). All-C++ (slower iteration; OpenCV already optimal for labeling). | t_polycore, t_reward, t_targets | — |

## Novelty, as of this record

1. **Exact per-instance topology of the model's own output, scored against a rater set.** A distribution-matching term keeps the policy's topology mix at rater frequencies instead of collapsing to one answer. Coverage bonus and KL were tested and rejected.
2. **A topology-conditional geometric term** (IoU_tc) that removes the hidden majority bias. We measured that bias as a +0.13–0.15 IoU gap.
3. **A self-calibrated, rater-count-corrected band (F5)** as an RL reward.
4. **Typed-verifier alarms with rollback** as the Level 2 controller. Narrow framing: typed sampling was shown not to help.
5. **A bit-exact geometry core** (C++ / CUDA-ready) with a single specification shared by SFT targets, reward and evaluation.

## Open risks

- Topology disagreement in the fixtures is thin (14 LIDC items, all with a single-rater minority).
- No real GRPO run yet.
- K8 meets the target only once adjusted for rater count.
- Lossy targets are not idempotent on large objects.
- VectorLLM++ reward terms are still unread.
- Wall-clock benchmarks ran on a loaded Mac; Vast x86 is not yet measured.
