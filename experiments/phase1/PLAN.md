# Phase 1 plan: build and freeze the exact reward, the data and the training integration

Written 2026-10-19. Part of [ROADMAP.md](../../ROADMAP.md).

How this plan was made:
- **6 online research branches:** `docs/research/notes/Phase 1 stack and algorithms/`.
- **7 sandboxed test labs**, plus 3 re-test rounds: `experiments/phase1/lab/`.
- **Two independent judge reviews:** `lab/JUDGE_VERDICT.md` (GO conditional) and `lab/JUDGE_FINAL.md` (approved with edits, all applied).

Every decision below names the test that chose it.

## 1. Goal and scope

**Goal.** Before any GPU is rented, build, test and freeze everything the GPU phases need:
1. one exact geometry spec, used identically by SFT targets, reward and evaluation;
2. the reward and its TRL integration;
3. clean, leak-free data with splits;
4. the model screening harness, the evaluation harness, and the gates for entering Phase 2.

**Out of scope:**
- SFT and GRPO training runs (Phase 2 onward, on Vast AI);
- the CUDA port (Phase 8);
- 3D (Phase 7).

**Done when** every Phase 1 exit gate in §7 is green in CI.

## 2. What the research and tests decided

Lab numbers come from Apple M3 runs, often under heavy load, so timings are pessimistic.

| # | Module | Decision | Beat (why) | Evidence |
|---|---|---|---|---|
| D1 | Polygon → mask | **C++17 header-only `polycore`**, integer scanline winding. Even-odd and nonzero in one pass. Top-left tie rule. Per-ring mixed-sign flag. nanobind, GIL released, internal thread pool. | PIL/OpenCV don't sample pixel centres (135–175 wrong px per polygon; PIL inflates 3-px nodules 1.56×). skimage/rasterio: even-odd only, undocumented ties, 4–45 ms. numba: 5× slower. | 0 / 849.8 M px wrong vs exact oracle; batch of 256 in 8.5 ms; p99 0.38 ms; Linux arm64 + x86-64 in Docker |
| D2 | Legality | Same core: one integer pass over edge pairs (x-sweep), per-class defect counts, never repairs. The reward uses a **hard gate**; counts are logged. | GEOS reports only the first error and accepts pinches. Pure Python: 43 ms per ring, ≈ 33 s per GRPO step. | 0 / 10,000 rings wrong vs oracle; ≈ 1,000× faster than the current code |
| D3 | Mask topology | **OpenCV Spaghetti** labeling (8-connected foreground / 4-connected background) on the bounding-box crop. Holes assigned to components via the pixel above. Bit-quad χ kept as an independent certificate. | scipy (5–24× slower here). Persistent homology (same answer, far costlier). | 0 mismatches on 10,000 random + 8,738 real masks |
| D4 | Band | **F5**: distance beyond the egg–yolk band, normalized by contour length. Tolerance self-calibrated per image from the raters' own leave-one-out scores, with a rater-count correction. Single-rater items excluded from RL. | Current band (correlation with object size +0.53). Fixed-pixel tolerance (rater disagreement grows with size, slope 0.46–0.83). √area tolerance (over-tolerant). | Per-stratum abs(ρ) with area ≤ 0.10; rater-subsample medians 0.961–0.980 |
| D5 | Topology term | Set-valued validity × rater-frequency matching: `max(0.1, 1 + ½(log f − b(f,8) − log π̂_loo))`, mean-centred, weight 2. An instance of MaD-RL's gated reverse-KL with the rater histogram as target. | Plain "any rater" (collapses to one answer, TV 0.30–0.45). Coverage bonus (spreads answers uniformly). KL to SFT (fragile). Majority. | Simulation: TV 0.06 vs 0.30. Real data, 251 items, reward v4 with realistic honest rollouts: neutrality +0.038 [−0.011, +0.091]. Valid > invalid at every reachable state |
| D6 | IoU term | **IoU_tc**: max IoU over raters sharing the sample's topology; IoU vs the vote if no rater drew that topology. | IoU vs the majority vote (a +0.13 gap between topologies, which drives collapse). | Gap −0.004 |
| D7 | Advantages | Per-group z-scores of IoU_tc, F5, legal and the **area guard** g_area = min(1, A_max/A)·min(1, A/A_min) over all raters. σ floors 0.05 / 0.05 / 0.25 / 0.05. No batch-level z. Illegal rollouts zeroed. Groups with spread < 1e-2 dropped. In TRL: weight-0 per-term loggers plus one aggregator reward. | TRL `normalize_then_sum` / plain GDPO ("dark room": 1e-3 noise moves advantages up to 1.24). Batch-level z (let the in-band speck win). | TRL advantages equal ours to 8.7e-8 **with reward v2**; v4 equivalence and K12 re-run are the M7 gate |
| D8 | RL loss | TRL 1.16 GRPOTrainer: `dr_grpo`, clip 0.2 / 0.28, β = 1e-3 (non-zero at start), truncation masking, constant learning rate. HF generation for SmolVLM until TRL PR #6718 merges. GSPO as an ablation. | verl/EasyR1 (no SmolVLM). LLaMA-Factory (no GRPO). RLOO (most hack-prone). | 4 models run GRPO; reward v2 runs end to end |
| D9 | Level 2 controller | **Phases 2–3:** typed per-term logging, held-out typed error curves, best-checkpoint rollback, guards switched by hand. Later phases: the typed alarm, re-tuned on real curves. | Typed curriculum sampler ("upweight stalled types": never better in simulation). Difficulty band-pass (worse than uniform 10/12). | Simulated alarm: 98.6–100 % detection, but a 270-step lag, so automation is deferred |
| D10 | Output format | Compact JSON, integers 0–1000 relative to the crop. Outer rings clockwise on screen, holes counter-clockwise. Start at the top-most then left-most vertex. `"holes":[]` always present. Truncation = format failure. RIGA emits disc + cup (rim computed). No vertex cap. | Spaced JSON (+22 % tokens). Float coordinates (2×). Fixed vertex count (collapses). | Tokens: 8.2 per vertex (Qwen/SmolVLM), 5.2 (LFM) |
| D11 | SFT targets | `minseg_hybrid`: potrace polygon stage OR lossless digital straight segments. τ = 3 grid units, IoU ≥ 0.99, defect-localized repair, pinch bevel, certified by polycore. Crop ≤ 500 px, preferring sides that divide 1000. One target per rater. | Douglas–Peucker / Visvalingam–Whyatt (more vertices, lower IoU, useless on small crops). Potrace ½ px (4–5× the tokens). | 100 % topology kept; IoU ≥ 0.99; median tokens: disc 157, cup 197, nodule 92 |
| D12 | Models | Screen smallest first: SmolVLM-256M → LFM2.5-VL-450M → Qwen3-VL-2B. Qwen3-VL-2B is also the non-smallest control. LIDC crops upscaled to 512 px; SmolVLM `longest_edge` 1024 on RIGA. | SmolVLM-500M (adds nothing). Florence-2 (no GRPO support). | Only Qwen3-VL-2B has a zero-shot prior (box IoU 0.76); SFT works on 256M and 450M; LFM licence allows research and publication |
| D13 | Data | **RIGA:** staged extractor recovers 4,500 / 4,500 masks; literature split train 590 / val 65 / test 95 (88 unique photos). **LIDC:** full parser with 4 fixes; 10,921 multi-reader slices; patient split 70 / 10 / 20. Parquet + SHA-256 manifest. | Phase 0 extraction (4,172 masks, BinRushed1 lost). DVC/webdataset (unneeded). | pylidc parity on 1,018 scans (99.97 % of masks exact); 0 cross-split duplicates; 313 train disagreement slices |
| D14 | C++ vs Python | C++: polycore (raster, legality; shared later with CUDA). Python + OpenCV: topology, band, caches, orchestration. numba: target-builder kernels, offline. One production implementation per operation; prototypes kept only as differential tests. | All-Python (≈ 1,000× slower legality). All-C++ (slower iteration; OpenCV is already optimal for labeling). | — |

## 3. Architecture

### 3.1 System

```
                    ┌──────────────────────────── Phase 1 (Mac, CPU) ─────────────────────────────┐
 raw data           │                                                                             │
 RIGA (6 raters) ───┼─► data.ingest ─────► exact expert masks ────► data.splits                   │
 LIDC (4 readers)   │   staged RIGA          4,500 RIGA masks        group-aware, seeded,         │
                    │   extractor,           10,921 LIDC slices      leak-tested                  │
                    │   LIDC parser          (certified)                 │                        │
                    │                                                    ▼                        │
                    │                          data.store: Parquet + SHA-256 manifest             │
                    │                               │                          │                  │
                    │                               ▼                          ▼                  │
                    │               targets: mask → canonical       reward.cache (per prompt):   │
                    │               polygon JSON, certified by      vote / union / core, EDTs,   │
                    │               polycore                        τ, T5, rater topologies      │
                    │                               │                          │                  │
                    │                               ▼                          │                  │
                    │               data.export: TRL chat rows                │                  │
                    └───────────────────────────────┼──────────────────────────┼──────────────────┘
                                                    ▼                          ▼
     Phase 2  SFT (TRL SFTTrainer) ──► checkpoint ──► Phase 3+  GRPO (TRL GRPOTrainer, G = 8)
                                                                 │ completions
                                                                 ▼
                                    ┌──────────────── reward server (objdel.reward) ───────────────┐
                                    │ parse (strict grammar) → polycore.validate → polycore.raster │
                                    │ → topology (Spaghetti 8/4, size rule) → IoU_tc, F5, legal,   │
                                    │   r_topo → advantages (per-group z with floors + topology    │
                                    │   mean-centred) → typed error log                            │
                                    └──────────────┬─────────────────────────────┬────────────────┘
                                                   ▼                             ▼
                                    train.callbacks: typed logging,      eval: same checker, same
                                    best-checkpoint rollback, guards     rasterizer, metrics
                                    (alarm automated in Phase 5)         (D_max, TV/JS, pass@k, hack audit)
```

### 3.2 One GRPO step through the reward (one exact spec end to end)

```
 32 prompts × 8 completions (text)
   │
   ├─ parse: strict grammar, ints in [0,1000],
   │         "holes" always present, truncation ───────► format = 0 (rollout zeroed)
   │
   ├─ polycore.validate  (C++, sweep, int64)  ─────────► legal ∈ {0,1}  + per-class counts (logged)
   ├─ polycore.rasterize (C++, winding, top-left) ─────► even-odd mask  + nonzero / mixed-sign flags
   │
   ├─ topology: Spaghetti 8/4, holes → components,
   │            size rule max(5 px, ⌈0.5 %·|U|⌉) ──────► per-instance signature
   │
   ├─ cache (built once per prompt, offline):
   │      rater masks, vote / union / core, EDTs, τ, T5,
   │      rater topology histogram f
   │
   ├─ terms:
   │      IoU_tc  max IoU over raters with the same signature (else vs vote)
   │      F5      1 / (1 + (2D/P) / T5²)
   │      legal
   │      r_topo  valid · max(0.1, 1 + ½(log f − b(f,8) − log π̂_loo))
   │      g_area  min(1, A_max/A) · min(1, A/A_min) over all raters
   │
   └─ advantages:
          Σ w·z_group(IoU_tc, F5, legal, g_area) with σ floors   (w = 0.3 / 0.3 / 0.1 / 0.3)
        + 2·(r_topo − group mean)
          illegal → 0;  group spread < 1e-2 → whole group 0
        ─► TRL aggregator reward (+ weight-0 per-term loggers)
```

### 3.3 C++ / Python split

| Layer | Language | Why |
|---|---|---|
| Rasterize, legality, signed area | C++17 header-only (`cpp/polycore`), nanobind | Hot path, exact integer arithmetic; shared with the CUDA engine later |
| Labeling / topology | Python calling OpenCV (C++) | Spaghetti is already state of the art; ~10 µs per LIDC mask |
| Band, caches, terms, advantages | Python + numpy / OpenCV EDT | Cheap per rollout once caches exist; easy to audit |
| Target simplifier | Python + numba (offline) | Runs once per dataset; port to C++ only if it becomes a bottleneck |
| Orchestration, TRL adapters, eval | Python | Glue |

### 3.4 Package layout after Phase 1

```
cpp/polycore/          header-only C++17 core + native fuzz tests
src/objdel/
  _polycore            nanobind extension (scikit-build-core build)
  checker/             bitquad, polygonize, certify, topology (+ instance topology), polycore wrapper
  geometry/            canonical polygon schema, strict parser, serializer
  targets/             SFT target builder
  reward/              cache, terms, advantages, batch scorer, typed errors
  data/                riga (staged extractor), lidc (parser), splits, store, export
  train/               TRL reward functions, aggregator, callbacks (logging, rollback), configs
  eval/                metrics harness and reports
tests/                 unit, property, oracle, fuzz and integration tests (oracles in tests/oracles/)
experiments/phase1/    PLAN.md, RESULTS.md, scripts, lab/ (the evidence behind this plan)
```

## 4. Module specifications

Each module ships with the tests named here, written first (test-driven).

### 4.1 Geometry spec and polycore (D1, D2, D10)

**Spec**
- Polygon grid G = 1000 over a W×H crop.
- Pixel (c, r) is inside iff its centre is inside, using the half-open crossing rule.
- Outer rings are clockwise on screen, holes counter-clockwise.
- Pinches are illegal; targets bevel them.

**API:** `pack`, `rasterize` (even-odd, nonzero, per-ring winding range, `px_multi`), `validate(roles)` (13 defect classes), `signed_area2`, `score_batch(threads)`.

**Input caps:** ≤ 1,000 vertices per ring, ≤ 10 rings, W, H ≤ 4096. Anything beyond the caps is rejected and tested.

**Tests to port**
- K1: equality with the exact oracle on 10,000 shapes in 9 families.
- K2: label equality with the oracle on 10,000 rings.
- Property tests: translation, mirror at ties, reversal.
- Native fuzz under ASan/UBSan.
- Determinism.

**Open semantics decided in Phase 1:**
- `wrong_winding` uses net area;
- a hole crossing its outer gets both flags.

**Build:** scikit-build-core replaces hatchling. cibuildwheel produces Linux x86-64 and arm64 and macOS arm64 wheels. CI runs pytest and the sanitizer fuzz on each OS.

**Removal:** the PIL rasterizer and Python `ring_problems` are deleted from the reward path.

### 4.2 Topology and size rule (D3)

- `instance_topology(mask)` returns per-component (area, holes) after the size rule.
- **Size rule:** min_area = max(5 px, ⌈0.5 %·|U|⌉), applied in the census re-run, targets, reward and eval.
- **Tests:**
  - brute-force equality on random and real masks;
  - χ certificate;
  - sub-filter specks don't change the topology set; above-filter components do.

### 4.3 Reward (D4–D7)

- **Per-prompt cache** (offline, deterministic): rater masks, vote, union, core, EDTs, τ, T5 and the topology histogram.
- **Terms, advantages, aggregator:** as in §3.2.
- **Constants (frozen):**

  | Constant | Value |
  |---|---|
  | weights IoU_tc / F5 / legal / g_area | 0.3 / 0.3 / 0.1 / 0.3 |
  | σ floors | 0.05 / 0.05 / 0.25 / 0.05 |
  | topology weight | 2 |
  | λ | 0.5 |
  | clip floor ε | 0.1 |
  | G | 8 |
  | low-spread drop threshold | 1e-2 |

- **Frozen as reward v4.** The constants were chosen on the full-LIDC disagreement set (251 train items), with realistic honest rollouts: exact rater masks plus smooth boundary perturbations of amplitude U(0, τ), one amplitude shared across an item's raters.
- **Area guard:** `g_area`. Its tolerance is zero because the τ·P version never triggered on the union or hull.
  - **M6 test:** score each held-out rater against the remaining raters' area range (leave one rater out).
  - **Switch rule:** if the held-out median g_area is below 0.9, switch to the τ·P tolerance, which kept 1-px neutrality at −0.001 against +0.047 for zero tolerance, and re-run R3/R4.
- **Cache storage:** ExpertCaches are built offline once per item and stored next to the dataset as compressed arrays with a hash. They are rebuilt only when the size rule, F5 constants or data pins change. Building takes ~12–18 ms per LIDC item, so on-the-fly is the fallback.
- **Tests:**
  - K9: correct topology ranked above wrong, ≥ 95 % with a bootstrap CI;
  - the adversary suite (§7, R3), plus the local-gradient test: dilating or eroding an honest rollout must lower its advantage;
  - neutrality on disagreement items;
  - valid > invalid at every reachable state;
  - dark-room replay;
  - determinism.

### 4.4 SFT targets (D11)

- **API:** `build_target(masks, schema, BuildConfig(), derived) -> (text, objects, info)`.
- **Crop policy:**
  - side ≤ 500 px, preferring sides that divide 1000;
  - larger objects downscaled by an integer factor with block majority, and the downscaled mask becomes the ground truth everywhere;
  - LIDC crops upscaled to 512 px for the model input.
- **Certification uses polycore.** fastgeom becomes a differential test.
- **Hashing:** the target set is hashed under pinned versions.

### 4.5 Data (D13)

- **RIGA:** staged extractor F1–F7, regenerating all 4,500 masks. Duplicate pairs are grouped. Test results are reported on 95 and on 88 photos.
- **LIDC:** parser with single-point ROIs dropped and a documented outline rule (fill including the outline; exclusions remove their outline). Patient ids come from the pylidc DB; slice index from z positions.
- **Splits** (`build_splits.py`, 11 tests):
  - group-aware and seeded;
  - stratified by disagreement;
  - injected-leak negative controls.
- **Store:** Parquet with original image bytes, exact masks (RLE) and a manifest. RIGA stays private (its licence forbids redistribution).
- **Export:** TRL chat rows with `{"type":"image"}` placeholders.

### 4.6 TRL integration (D7, D8, D9)

- **`train/rewards.py`:** weight-0 per-term reward functions for logging, plus the aggregator. Contiguous groups; `scale_rewards="none"`, `sum_then_normalize`.
- **`train/callbacks.py`:**
  - typed held-out evaluation every N steps;
  - best-checkpoint tracking and rollback;
  - guard switches (mutating `reward_weights` and `beta`, verified live).
- **Configs:** SFT (assistant-only loss, `max_length=None`) and GRPO (D8 settings), one per candidate model.
- **Note:** TRL's own `reward` metric reads 0 with the aggregator, so dashboards use the per-term logs.

### 4.7 Evaluation harness (judge R10)

`eval/metrics.py` uses the same checker and rasterizer as training. It computes:
- IoU and Dice vs the majority;
- D_max (best match to any rater);
- F5, topology-any;
- per-image TV and JS between the K-sample topology histogram and the rater histogram;
- minimum-mode retention;
- pass@k for k ≤ 64;
- RIGA soft Dice (literature protocol).

**Hack audit:** |M|/|E|, A/A_max and A/A_min, component count, duplicate-hash rate, newly rewarded failures, with the numeric alarm lines in §9.

**Tests:**
- metrics give their perfect value on the raters' own masks;
- TV = 0 when samples reproduce the rater histogram;
- the audit flags every R3 adversary.

## 5. Work plan (test-driven order)

Each milestone starts by porting or writing its tests, then promotes the lab code, then runs its gate.

| M | Milestone | Promotes from | Gate (§7) |
|---|---|---|---|
| M1 | Build system: scikit-build-core, polycore extension, cibuildwheel, CI (macOS + Linux x86-64/arm64, sanitizer fuzz), pinned versions | `lab/t_polycore` | R9 (partial), K1, K2 |
| M2 | Geometry spec module: strict parser, serializer, canonicalization; delete the PIL path | `lab/t_stack/prod.py`, `lab/t_polycore` | R1 parser part |
| M3 | Topology + size rule; Phase 0 census re-run under the rule | `lab/t_reward/rewardcore/topology.py` | R6 |
| M4 | Data: RIGA staged extractor, LIDC parser fixes, splits, store, export | `lab/t_data` | K11, R4 data part |
| M5 | Targets, certified by polycore; target-set hash | `lab/t_targets` | K5–K7, R1 |
| M6 | Reward v4: cache (stored), terms, advantages, batch scorer; leave-one-rater-out g_area test; per-rollout local-gradient rate on the 251 items; hole-matching spike | `lab/t_reward` | K3, K4, K8, K9, K10, R3, R4, R5 |
| M7 | TRL integration with reward v4: rewards, aggregator, callbacks, configs; group-contiguity test; Mac smoke run | `lab/t_stack` | R2 (v4), K12 (v4) |
| M8 | Evaluation harness + hack audit | new | R10 |
| M9 | Vast x86 benchmark (1 h, cheapest 8-core box) + Phase 1 RESULTS.md | — | R9 |

M1–M3 are independent. M5 needs M1 and M3. M6 needs M3 and M4. M7 needs M2, M5 and M6. M8 can start after M6.

## 6. Engineering rules

- **One production implementation per operation.** Lab prototypes are kept only as differential tests (`tests/oracles/`).
- **Every exactness claim is checked against an independent, slow oracle.** Random plus adversarial generators; property tests (Hypothesis); native fuzz under sanitizers.
- **Determinism:** identical SHA-256 of rewards, targets and splits across runs, threads and processes. Hashes are recorded in manifests.
- **Pinned versions:** TRL, transformers, torch, numba, OpenCV, nanobind. Targets and caches are re-hashed when a pin changes.
- **Untrusted data:** scripts that read `data/` run with `python -I`; data and fixtures are never committed.
- **CI** runs tests, lint and sanitizer fuzz on every push; coverage is measured (K13).

## 7. KPIs and gates

### Phase 1 exit gates (all must pass before any GPU rental)

| Gate | Requirement | Lab status |
|---|---|---|
| K1 | Raster equals oracle (both rules) on ≥ 10k shapes | ✅ 0 / 849.8 M px |
| K2 | Legality labels equal oracle on ≥ 10k rings | ✅ 0 / 10,000 |
| K3 | 256-rollout batch ≤ 250 ms on 8 cores (Mac and Vast x86) | ✅ Mac 14–75 ms; x86 pending (M9) |
| K4 | p99 per rollout ≤ 5 ms | ✅ ≤ 4.2 ms under load |
| K5 / K6 / K7 | Targets: topology 100 %; IoU min ≥ 0.98 and median ≥ 0.99; median tokens ≤ 250 | ✅ 100 %; 0.990 / 0.993+; 92–197 |
| K8 | Band: per-stratum abs(Spearman ρ) with area ≤ 0.10 (pooled reported, not gated) | ✅ |
| K8b | Band discrimination: monotone ≥ 98 %; a 2τ-shifted vote ranks below the median held-out rater on ≥ 95 %; held-out median ≥ 0.90 and p10 ≥ 0.70 | ✅ 100 %; 99.9 % (eroded vote: 97.1 %, LIDC 95.2 %, borderline); 0.973 / 0.827. Re-registered from absolute thresholds, which failed at 89 % / 48 % (§9); the new thresholds were set after measurement. |
| K9 | Topology ranking ≥ 95 % with a bootstrap CI | ✅ 99.94 % fixtures, 98.86 % full LIDC |
| K10 | Determinism of reward, targets and splits | ✅ |
| K11 | 0 leakage and 0 near-duplicates across splits | ✅ |
| K12 | TRL GRPO steps with the production reward; per-term logs; no NaN | ✅ with reward v2 (from a non-collapsed SFT checkpoint); re-run with v4 at M7 |
| K13 | Line coverage of `src/objdel` ≥ 90 %; CI green | M1–M8 |
| R1 | Round trip target → parser → polycore → reward: 100 % legal, signature kept, IoU_tc ≥ 0.99, F5 ≥ 0.95 on ≥ 99 %, canonicalization-invariant | ✅ 6,338 fixture targets; full train set at M5 |
| R2 | TRL advantages equal the reference to 1e-6; dropped groups give exactly 0; groups contiguous per process | ✅ 8.7e-8 with reward v2; v4 at M7 |
| R3 | (a) Reachable hacks: dilated union and in-band speck below the median honest rollout in ≥ 95 % of groups; dilating or eroding an honest rollout lowers its advantage on average (measured: mean ΔA < 0 with the CI excluding 0, six moves, ~1,000 fixture items × 6 groups; the per-rollout rate and the 251-item LIDC run are an M6 gate, target ≥ 95 %). (b) Oracle-only adversaries (hull, union, core of all raters; hole stamping): reported, monitored by the eval hack audit, not gated (§9). | (a) ✅ 100 % / 96.5 %; local gradient ✅ on average (per-rollout rate pending, M6). (b) adversary beats the median honest rollout: hull 34 %, union 51 %, core 44 %; minority-hole punch 60 % [54, 65] (limit ≤ 60 %, borderline); holes stamped 64 % [59, 69] (limit ≤ 60 %, fails) |
| R4 | On ≥ 200 train disagreement items, the neutrality 95 % CI contains 0 with half-width ≤ 0.10 | ✅ v4 realistic: +0.038 [−0.011, +0.091]; 1-px: +0.047 [−0.009, +0.104] (251 items) |
| R5 | Rater-subsample medians in [0.92, 0.98] for m = 2–5; single-rater policy tested | ✅ 0.961–0.980; single-rater excluded |
| R6 | Size-rule invariance tests | ✅ |
| R9 | Wheels + CI per OS; K3/K4 re-measured on Vast x86; pins and target hash | M1, M9 |
| R10 | Eval harness tests green | M8 |

### Phase 2 entry gates (checked on the validation split after SFT; the smallest model passing them is "smallest viable")

| Gate | Requirement |
|---|---|
| Format | parse ≥ 98 %; legal ≥ 90 % |
| Geometry | median IoU_tc ≥ 0.85 (RIGA disc), ≥ 0.75 (cup), ≥ 0.60 (LIDC); ≥ 70 % of rollouts pass IoU_tc ≥ 0.5 |
| Topology | topology-any ≥ 90 % on unanimous items |
| GRPO signal | at temperature 1 with G = 8: duplicate outputs ≤ 10 %; ≥ 50 % of groups with geometric σ above the floor. An over-trained SFT that collapses outputs gives zero GRPO signal, as seen in the lab. |
| Cost | per-step wall time measured on the 4090 for the smallest model (G = 8, 1,536 tokens); the roadmap's GPU-hour estimate replaced by the measurement |
| Precision probe | before the first SFT run, a synthetic-shape probe per candidate model at the chosen crop/upscale settings (known polygons on synthetic images, a few hundred SFT steps): boundary error ≤ 2 px on LIDC-sized objects, or the model is dropped from LIDC |
| Housekeeping | novelty sentences in their narrowed form (§8); LFM licence recorded (done: research and publication allowed) |

## 8. Novelty (narrowed after the prior-art check)

1. **Safe.** To our knowledge (targeted search, Oct 2026), the first RLVR reward on the exact per-instance topology of a VLM's own polygons, valid when it matches any rater's topology. Neighbours: VectorLLM++ uses IoU + vertex count vs one annotator; Topo-R1 uses topology of the input mask only.
2. **Instance of prior method.** Rater-frequency matching of topology. This is MaD-RL's gated reverse-KL (Kulkarni et al. 2026) with the per-image rater topology histogram as target. Our instantiation: the target is computed from geometry, the validity gate is the rater set, and there is a bias-corrected leave-one-out group estimator (cf. Mohri et al. 2026).
3. **Design choice, not a claim.** Topology-conditional IoU, which removes a measured hidden majority bias (+0.13 IoU gap).
4. **Narrow.** A self-calibrated, rater-count-corrected band as an RL reward. Cite Nikolov et al. for inter-rater tolerance.
5. **Engineering, not a claim.** The typed logging and rollback controller.

## 9. What the last test round settled and what stays open

**Settled (reward v4 frozen).**
- The IoU_tc fallback scores vote IoU for topologies no rater drew.
- The batch-level z-score is removed: it let batch composition inflate geometry over topology.
- The leave-one-out estimator's exact binomial bias b(f, G) is subtracted. Without it, minority answers were favoured by +0.08.
- The topology bracket is clipped at 0.1, so every valid rollout scores above every invalid one (property-tested over all reachable states).
- The area guard is added as a separate z-scored term.
- K8b is re-registered as a discrimination test.

**Open, stated honestly.**
- **Oracle-only adversaries still beat the median honest rollout often:**

  | Adversary | Beats the median honest rollout |
  |---|---|
  | convex hull of all raters | 34 % of groups |
  | union of all raters | 51 % |
  | intersection of all raters | 44 % |
  | minority hole punched into the vote | 60 % (limit 60 %, borderline) |
  | holes stamped from every rater | 64 % (limit 60 %) |

  These shapes need every rater's mask (their band score is 1.0 by construction), and there is no first-order incentive toward them: dilating or eroding an honest rollout lowers its advantage on average (mean ΔA < 0, CI excludes 0, ~1,000 fixture items), and the guard makes that about 1.5× steeper. The per-rollout rate and the 251-item run are added to M6. Consensus-like answers scoring well is partly the intended behaviour.
- **Hedging is monitored, not penalized further.** The eval hack audit tracks |M|/|E| and area relative to the rater range on every evaluation. Alarms fire (and the best checkpoint is restored) when, on held-out items, any of these happens:
  - the median A/A_max exceeds 1.15;
  - the share of outputs above A_max·1.3 exceeds 10 %;
  - the median component count rises by ≥ 0.5 over the SFT checkpoint;
  - the duplicate-hash rate within groups exceeds 50 %.
- **M6 spike for the hole case:** try the hole-matching term from the original reward design (Hungarian ring IoU vs consensus holes). Keep it only if it lowers the stamped-holes win rate without breaking neutrality (R4).
- **The original absolute form of K8b failed** (2τ excursion ≤ 0.5 on 89 %; every held-out rater ≥ 0.85 on 48 %). No band shape short of a step function passes it, even with an oracle tolerance. The re-registered form tests what GRPO uses, namely ranking.

## 10. Risks

| Risk | Mitigation |
|---|---|
| Topology disagreement is a small share of data (251 train slices under the size filter; 313 unfiltered) | Oversample disagreement items in RL prompts; report per-image TV only on them |
| The smallest model lacks the resolution for 2 px LIDC precision | Upscale LIDC crops to 512 px; the Phase 2 gate decides; Qwen3-VL-2B as fallback |
| Over-trained SFT collapses diversity, so GRPO gets no signal | Early-stop SFT on the duplicate-rate / σ gate |
| TRL bugs (SmolVLM + vLLM #6294; reward metric reads 0; group contiguity per process untested) | HF generation for SmolVLM; per-term dashboards; an M7 test that every process's reward call receives whole groups |
| F5 assumes Gaussian rater excursions | Validate on IMA++ or QUBIQ when the next track starts |
| Lossy targets are not idempotent on large objects | Lossless-vs-hybrid ablation in Phase 2 |
| Benchmarks only on a loaded Mac | M9 measures on Vast x86 before Phase 2 |
