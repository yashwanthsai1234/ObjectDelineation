# Phase 1 results: exact reward, data and training integration

Written 2026-10-19. Plan: `PLAN.md`. Coordination: `EXECUTION.md`. Per-task notes:
`impl/`. Issue ledger: `impl/ISSUES.md`.

**Status: done.** Every Phase 1 exit gate passes on the Mac, or is re-registered or accepted with a recorded reason
(§2). The exactness and determinism gates also pass on x86-64 in Docker (§6). Two checks need hardware Phase 1 did
not rent, and are Phase 2 entry steps (I34): the first GitHub Actions run, and the native x86 timing and hash check on
the first Vast box. No GPU was rented in Phase 1.

Evidence conventions:
- `M6 §Result` means `impl/M6_reward.md`, section "Result".
- `F2 §…` means `impl/fixes/F2.md`; `B1`–`B4` are `impl/breakit/`; `S1`–`S4` are `impl/S<n>.md`.
- Test names are pytest node names under `tests/`.
- Hashes are SHA-256.

## 1. What Phase 1 delivered

**Library (`src/objdel`)**

| Module | One line | Note |
|---|---|---|
| `cpp/polycore/` + `_polycore` | Header-only C++17 core: exact integer rasterizer (even-odd + nonzero, top-left ties), 13-class legality, exact signed area, threaded `score_batch`; nanobind, abi3, GIL released | M1 |
| `checker/` | Instance topology (OpenCV Spaghetti 8/4), size rule `max(5, ⌈0.5 %·|U|⌉)`, `as_mask` normalisation, bit-quads, certify, polygonize, polycore wrapper | M3, F1 |
| `geometry/` | Output spec: strict regex parser (never raises, `MAX_CHARS` 16,384), compact serializer, canonical rings, `to_mask`, portable `pattern()` for constrained decoding | M2, F1 |
| `targets/` | Model frame (square crop, side divides 1000, LIDC upscaled to 512 px) and certified SFT targets (per-object tolerance ladder, numba kernels), `min_areas`, `target_set_hash` | M5, F1, S1 |
| `reward/` | Reward v4: stored `ExpertCache` (vote/union/core, EDTs, τ, T5, voters), per-rollout terms, `advantages_v4`, fingerprint and `StaleCacheError` | M6, F2, S2 |
| `data/` | RIGA staged extractor (4,500 / 4,500 masks), LIDC parser (pylidc parity), group-aware splits with NCC + ORB near-duplicate checks, Parquet store with SHA-256 manifest | M4, F4, S3 |
| `train/` | TRL rows (`prepare`, `sft_rows`, `grpo_rows`, cache loader), `TRLReward` aggregator, `Controller` (best checkpoint, rollback, audit alarm, live switches), configs, sampler/exporter | M7, F3, S4 |
| `eval/` | Metrics harness (IoU_tc, D_max, F5, TV/JS, pass@k, soft Dice, GED), hack audit, `alarms_for`, `phase2_entry` per dataset | M8, F2 |

**Scripts (`experiments/phase1/`, run as `.venv/bin/python -I -m experiments.phase1.<name>`)**

| Script | One line |
|---|---|
| `fetch_lidc_pixels` | Fetches CT slices for every LIDC series with contours (SOP UID, then z match) |
| `build_data` | Builds `data/store` (Parquet + manifest) and `build_report.json` (K11 lists, counts) |
| `census_sizerule` | Phase 0 topology census re-run under the size rule (M3 table) |
| `build_targets` | Builds every train/val SFT target, K5–K7 summary, target-set hash (`data/targets`) |
| `prepare_rl` | Builds `data/prepared`: PNG frames, caches, per-split JSONL records, manifest |
| `reward_gates` | K3, K4, K8, K8b, K9, K10, R3, R4, R5, AB2–AB4 tables (writes the M6 gate section) |
| `smoke_grpo` | K12: short SFT, then GRPO steps with the production reward (Mac, MPS) |
| `sample_eval` | Checkpoint → K samples per item, or TRL completions table → `eval_run` JSONL |
| `eval_run` | Completions JSONL + store → `report.json` / `report.md` / `items.jsonl` |
| `eval_demo` | Synthetic models (honest, raters, R3 adversaries) over the val store |
| `bench/x86.sh`, `bench/xarch_hashes.py` | x86-64 build, tests, wheels and cross-architecture hashes in Docker (M9) |

**Removed:** the Phase 0 PIL rasterizer `checker/polygon_io.py` (M2), `reward/annotators.py` (`AnnotatorSet`, M6),
`betti3d`/`FG8`/`BG4` (M3), and `src/objdel/models/` with the `models` extra (S4 §medsam.py, coordinator decision).
Phase 0 scripts run against commit `faaacc4` (`experiments/phase0/README.md`, F5 §I3).

## 2. Gate table (PLAN §7)

PASS = meets the plan's requirement. RE-REG = the requirement was changed after measurement, with the reason
recorded in PLAN §9. ACCEPTED = fails as written; a recorded decision accepts it with monitoring.

| Gate | Target | Measured (final) | Verdict | Evidence |
|---|---|---|---|---|
| K1 | Raster = exact oracle, both rules, ≥ 10k shapes | 0 / 849,802,293 px on 10,000 shapes (9 families); 0 on x86-64 (Docker) | PASS | M1 §Gate results; `test_polycore_oracle.py::test_k1_full_10k_shapes` (slow); M9 §Build |
| K2 | Legality labels = oracle, ≥ 10k rings | 0 / 10,000 rings; fine classes 0 / 3,000; multi-ring 0 / 1,500; 0 on x86-64 | PASS | M1 §Gate results; `test_k2_full_10k_rings` (slow); M9 §Build |
| K3 | 256 completions ≤ 250 ms (Mac; x86 is a Phase 2 entry gate) | serial: mixed 48.6 ms, RIGA-only 103.5 ms (load 9); 107.6 / 158.3 ms at load 21 | PASS (Mac) | M6 §K3, K4, K10; F2 §Gate re-run |
| K4 | p99 per rollout ≤ 5 ms | 0.53 / 0.60 ms (load 9); 2.12 / 1.06 ms (load 21) | PASS (Mac) | same |
| K5 | Target topology kept 100 % | 1.0 for disc, cup, rim and nodule, train and val (size-filtered signature) | PASS | `data/targets/summary.json` `objects.*.K5`; M5 §Gates |
| K6 | Target IoU min ≥ 0.98, median ≥ 0.99 | RIGA disc 0.9900 / 0.9922, cup 0.9900 / 0.9915 (train); LIDC 1.0 / 1.0 | PASS | `summary.json` `iou_min`, `iou_median` |
| K7 | Median tokens per object ≤ 250 | train, Qwen / SmolVLM / LFM: disc 163 / 166 / 109, cup 147 / 150 / 99, nodule 100 / 103 / 70 | PASS | `summary.json` `objects.*.tokens` |
| K8 | Per-stratum abs ρ(F5, area) ≤ 0.10 | RIGA disc +0.008, cup −0.040; LIDC k=3 +0.012, k=4 −0.023 | PASS | M6 §K8, K8b, R5; unchanged in F2 |
| K8b | Monotone ≥ 98 %; 2τ vote below median held-out ≥ 95 %; median ≥ 0.90, p10 ≥ 0.70 | 100 %; 100 % (eroded 99.2 %); 0.975 / 0.824 | PASS (RE-REG in the lab, before Phase 1: the absolute form failed at 89 % / 48 %) | M6 §K8b; PLAN §9; `test_k8b_band_discrimination` |
| K9 | Topology ranking ≥ 95 % with CI | fixtures + composites 99.9 % [99.8, 100.0]; full-LIDC disagreement 99.3 % [98.5, 99.9] (731 pairs, 228 items) | PASS | F2 §Gate re-run; `test_k9_full_lidc_disagreement_items` |
| K10 | Identical SHA-256 of rewards, targets, splits | Mac: identical across serial, 8 threads, 3 processes, OpenCV threads, fresh interpreter; store rebuilt twice byte-identical (manifest `37fe7665…`); target-set hash identical across workers. x86-64 (Docker), after M9 fixes F1–F3: run 1 = run 2 on 26 / 26 hash sections; Mac = x86 on 22 / 26, the other 4 being raw split scores whose decisions are identical (accepted, I33) | PASS (Mac and x86-64 in Docker) | M6 §K10; F4 §Store; S3 §Golden; M9 §Final run; ISSUES I30–I33; `test_k10_identical_rewards_across_runs_threads_processes_interpreters`, `test_edt_is_correctly_rounded_without_opencv`, `test_excursion_moment_is_frozen_and_equals_quad`, `test_lidc_upscale_is_pillow_bicubic` |
| K11 | 0 leakage, 0 near-duplicates across splits | All `build_report.k11` lists empty: 0 RIGA group leaks, 0 cross-split NCC or keypoint duplicates (max NCC 0.326, max inliers 16), 0 LIDC patient leaks, 0 identical CT slices | PASS | `data/store/build_report.json` `k11`; F4 §B3-8; `test_store_k11`, `test_store_k11_keypoints` (slow) |
| K12 | TRL GRPO steps with the production reward; per-term logs; no NaN | Final run on the final `data/prepared` (manifest `3747acbc…`). SFT: SmolVLM-256M, LoRA, 40 steps, 11.4 s/step, loss 0.93 → 0.46. GRPO: 3 steps, 69 s/step; every per-term log present; no NaN; 6 / 6 groups with non-zero advantages; duplicates 0 %; σ above the floor in 100 % of groups; format_ok 1.0; legal 0.44–0.63. Peak footprint 6.0 GB; 787 s wall | PASS (Mac, LoRA, micro-batch 1) | `data/runs/smoke_grpo_final/smoke.json` (`checks`, `sft`, `grpo.batches`); `/usr/bin/time -l` of the run (787.14 s real, peak footprint 6,009,933,088 B); first run: M7 §K12 smoke; I22 |
| K13 | Line coverage of `src/objdel` ≥ 90 %; CI green | `make cov` (final, Mac): 705 passed, 3 skipped, 0 failed, **93.03 %** (2,239 statements), 14 min; 91.85 % in a CI-like Linux container; the CI jobs pass in Linux containers, but the workflow has never run on GitHub (no push from Phase 1) | PASS (coverage); first GitHub run is a Phase 2 entry step (I34) | final verification (§9); F5 §B4-8; ISSUES I34 |
| R1 | Target → parser → polycore → reward: 100 % legal, signature kept, IoU_tc ≥ 0.99, F5 ≥ 0.95 on ≥ 99 %, canonicalization-invariant | Parser: 6,338 / 6,338 on every check. Prepared train sets: every voting rater's target legal, topo-any, IoU_tc ≥ 0.99, F5 ≥ 0.95; 10 / 23,469 LIDC targets come from non-voting (sub-filter) raters and are excluded from SFT | PASS | M2 §Gate; M5 §R1; S4 §Follow-up 1; `test_r1_round_trip_on_the_prepared_train_set[riga|lidc]` |
| R2 | TRL advantages = reference to 1e-6; dropped groups exactly 0; groups contiguous per process | max \|TRL − v4\| = 2.37e-7; identical and low-spread groups exactly 0; world size 1/2/4 contiguous; split-group config refused | PASS | M7 §Result; `test_r2_trl_advantages_equal_advantages_v4`, `test_every_process_receives_whole_groups` |
| R3 (a) | Dilated union and in-band speck below the median honest rollout ≥ 95 % | 100 %; 96.3 % [94.9, 97.6] (real honest rollouts) | PASS (speck CI lower bound 94.9 %) | F2 §Gate re-run; `test_r3_reachable_hacks_lose_to_the_median_honest_rollout` |
| R3 (a), local gradient | Plan: dilating or eroding an honest rollout lowers its advantage, ≥ 95 % per rollout | All moves 80.8 % (FAIL as written); moves away from the rater consensus 97.3 % (full LIDC) / 98.3 % (fixtures); every move's mean ΔA < 0 (CI −0.50 … −1.70) | **RE-REG after measurement** (I12): gate = moves away ≥ 95 % + negative mean ΔA for every move | M6 §The two gates; F2 §Gate re-run; PLAN §9; `test_local_gradient_moves_away_from_consensus_per_rollout`, `test_local_gradient_on_average_on_the_lidc_disagreement_items` |
| R3 (b) | Oracle-only adversaries: reported, not gated | hull / union / core below median 66.2 / 49.3 / 55.2 %; holes stamped win 61.6 % [56.5, 66.8] (limit 60 %); minority hole 56.6 % | reported; holes stamped ACCEPTED, monitored (I15, I16) | F2 §Gate re-run; M8 §Open issues |
| R4 | ≥ 200 disagreement items; neutrality 95 % CI contains 0, half-width ≤ 0.10 | 239 items. Realistic: +0.033 [−0.013, +0.078] in the gate run; paired re-draws +0.062 … +0.067 (CI excludes 0); 1-px +0.019 … +0.028 (accept) | **ACCEPTED, monitored** (I13): mild pro-minority lean (≈ 0.04 advantage units), v4 frozen; per-image TV/JS tracked in Phase 3 | M6 §The two gates; F2 §R4, paired; `test_r4_neutrality_on_the_lidc_disagreement_items` (half-width ≤ 0.10, \|mean\| ≤ 0.07, fixed draw) |
| R5 | Rater-subsample medians in [0.92, 0.98], m = 2–5; single-rater policy | 0.961 / 0.971 / 0.976 / 0.980 (m = 5 at the edge); k = 1 refused by `ExpertCache.build` | PASS | M6 §R5; `test_r5_rater_subsample_medians` (slow) |
| R6 | Size-rule invariance | Specks/pinholes of A−1 px leave the signature set unchanged, A px change it; filter = oracle, idempotent | PASS | M3 §Result; `test_sub_filter_specks_and_pinholes_do_not_change_the_signature_set`, `test_above_filter_components_and_holes_change_the_signature_set` |
| R9 | Wheels + CI per OS; x86-64 tests pass in Docker; pins and target hash | Wheels: manylinux_2_27/2_28 x86-64 + aarch64 (auditwheel, abi3audit strict) and macOS arm64. x86-64 Docker final run, ALL PASS: build, lint, K1/K2 at full N, ASan/UBSan fuzz, 626 fast tests, release wheel (57 tests); hashes as in K10, `target_set_hash` included. Re-run after the review fixes R1–R4: PASS | PASS (under Rosetta emulation; native x86 on Vast and GitHub CI are Phase 2 entry steps, I34) | M1 §Gate results; M9 §Final run, §Re-run after the review fixes |
| R10 | Eval harness tests green | Perfect values on the raters' own targets (85 real val items); every R3 adversary flagged; 21 / 21 mutations killed; `src/objdel/eval` 100 % covered. Since R2, an illegal answer scores 0 on every geometry metric, through the reward's own `legal_gate` | PASS | M8 §Result; fixes/R2.md §1; `test_raters_own_targets_score_perfectly_on_real_items`, `test_illegal_answers_score_zero_on_geometry` |

## 3. Data and hashes

**Store** `data/store` (git-ignored; RIGA stays local under its licence, EXECUTION §5).
- `manifest.json` SHA-256 **`37fe7665e43423333681390c597ea265ab0fda6c00e708f4aacbc7c138dc7707`**.
- Byte-identical in a second F4 build and in S3's rebuild on the swept code (S3 §Golden).
- Build: seed 0, 827 s (`build_report.json` `timing_s`).

| File | Rows | Bytes | SHA-256 |
|---|---|---|---|
| riga/train.parquet | 592 | 1,592,457,774 | 6421dd5ffc6530dc0abd1ce9848a8fc2ee4c82e6e8261d491e9c1af059738364 |
| riga/val.parquet | 63 | 164,228,791 | 01d6bb68f64c0b3cc08224ed7bf5d8e9deed1aca1aac5662a1cdcadbe29d15fb |
| riga/test.parquet | 95 | 487,025,710 | 09a59f23122b5d15ceea08819c569173b2458562f532eceed0cc056ec396794d |
| lidc/train.parquet | 7,341 | 1,959,937,188 | 23a754949f9749a8326100a758fc5214b75fb7cbc184ac6866a26e3aeb7f0467 |
| lidc/val.parquet | 1,084 | 287,925,611 | bf26148ce7e04cfc3d3bad5d3ba95fae1e6950c6a524d5c1b6adaed28cd3b2f4 |
| lidc/test.parquet | 2,200 | 588,319,426 | 5555fb2e6d9211250ff9d2ace24e5f705bd9005141be672466d66e762f243765 |

**Counts** (from `data/store/build_report.json`, except the RIGA test image count; see below).
- RIGA: 750 photographs, 4,500 certified masks.
  - Stages: A 4,456, B 41, C 3.
  - References: prime 699, maxrg6 (BinRushed1 reconstruction) 50, prime-is-a-copy 1.
- RIGA duplicates: 51 near-duplicate pairs (37 NCC + 14 keypoint-only), in 45 duplicate groups.
- RIGA split 592 / 63 / 95. The Magrabia test set has 86 unique photographs.
  - Model images: 547 / 58 / 95 (700 of 750). The 50 BinRushed1 (maxrg6) items have `image = None` (B3-1).
    `build_report.json` says 94 for test: `build_data` counted the F7 item (Magrabia Female image048, which keeps
    its image) as imageless. R1 fixed the counter (REVIEW #4); the store itself was always right, and the report is
    corrected at the next rebuild (fixes/R1.md §2).
- LIDC: 1,018 series, 1,010 patients, **10,625 items**.
  - Items per split: 7,341 / 1,084 / 2,200.
  - Patients per split: **709 / 100 / 201**.
  - Every item has CT pixels.
- LIDC rules applied at build time:
  - crowded frames dropped: 143 / 16 / 57 (B3-2);
  - copied reader layers merged: 225 / 40 / 73 (B3-3).
- LIDC disagreement items under the size rule: 299 / 44 / 85.
  - 10 / 0 / 6 of them disagree only through a sub-filter rater (B2-4).

**Prepared set** `data/prepared` (final).
- `manifest.json` SHA-256 `3747acbcda2bc0503065af47b63861dcc3474e96bda19fb9881b8b1114d14ff4`.
- Built from store `37fe7665…`. Reward fingerprint `b1eb9f51824fe9617341a34e0aaa470a4b75646e4fad035679c4ff772d96c81f`.
- Rebuilt after M9 fix F3 (LIDC upscale with Pillow bicubic, S1 §F3). Only the LIDC images and each LIDC record's
  `image_sha256` changed, so only the two LIDC jsonl hashes moved (previously `6de765f7…`, `57bf803e…`). The RIGA
  files, every `target_set_sha256`, the caches and the reward fingerprint are unchanged (`data/prepared/manifest.json`).
- Before F3, S4's sweep and follow-up reproduced the prepared set byte for byte (18,670 files, S4 §Golden comparison).

| Set | Items | Targets | Disagreement items | Skipped (no image) | jsonl SHA-256 | target_set_sha256 |
|---|---|---|---|---|---|---|
| lidc/train | 7,341 | 23,469 | 289 (3.9 %) | 0 | c19b190977cd19afc91bf6709e5478f0524f12f78c64d9a0f510fd5c0cb2420c | b481450688d3b8920562be3516ddd77401379288e4925f55c10954bd99c2c4c6 |
| lidc/val | 1,084 | 3,444 | 44 | 0 | bf73ac0db25274ab3a975400fabb8164526bccbb74628fd3ba9f40b13f7b6a13 | 138290299b5afcf72c10f45ff32171473ba8cc56c2b3b2c2853aae568d0e6c40 |
| riga/train | 547 | 3,282 | 0 | 45 (BinRushed1) | 987d73b84e689425bd33d284368bbf556c74b7f6da030fd8e9ee525aa729aa5e | 6b2f46da23b5c6da16e6f70da5513585ccff26e26c1bbbda64e980d24ec39f0f |
| riga/val | 58 | 348 | 0 | 5 (BinRushed1) | 1735e142fb607da76463b13d643549763fbeef097693f5ed951058c9018b2858 | 8a5653cc71dc40c20a52347186a894d0bdbe484a202c830deb5eabca0723ea6d |

- The 289 LIDC train disagreement items are 3.9 %, the ≈ 4 % rate the ×4 oversampling was designed for (M7 §Decisions; S4 §Cross-path reports).
- SFT rows: one per *voting* rater. 10 / 23,469 LIDC train targets come from sub-filter raters and are dropped from SFT; RIGA drops none (S4 §Follow-up 1).
- Over-budget answers in the prepared set: 1 RIGA train answer per model (S4 §Cross-path reports). It is dropped from SFT only (I5).

**Targets** `data/targets` (built from every store item, including the BinRushed1 items without images).
- `target_set_sha256` **`1af819f474536dffe5a4a33888ed343bebd3855383fc162b3ca00e0ab82e91d3`**. Its manifest records store `37fe7665…`.
- Per file (`data/targets/manifest.json`):

  | File | Targets | SHA-256 |
  |---|---|---|
  | lidc/train.jsonl | 23,469 | b481450688d3b8920562be3516ddd77401379288e4925f55c10954bd99c2c4c6 |
  | lidc/val.jsonl | 3,444 | 138290299b5afcf72c10f45ff32171473ba8cc56c2b3b2c2853aae568d0e6c40 |
  | riga/train.jsonl | 3,552 | 54bc88bfce0f7ee80d572ac88989442bb8556e3fcf6c5f0b00f2b579123a86ef |
  | riga/val.jsonl | 378 | 39b2f17131b59c8ae9a3cbadefaa770d7de8bd83a3851a104f7a1e4a51015190 |

- Quality: K5 = 1.0 for every object; legal 1.0; parse 1.0; IoU_min ≥ 0.99 (RIGA 0.9900, LIDC 1.0); 0 errors (`summary.json`).
- Config: τ = 3, step 0.25, IoU ≥ 0.99; frames with margin max(8 px, 25 %), side 50–500 px.
- Frame sides: RIGA 200 / 250 / 500 px; LIDC 50 / 100 / 125 px. No item needs a downscale (scale 1 everywhere).
- Build time: 115.7 s wall (`summary.json` `wall_s`).
- Whole-answer tokens, median / p99 / max (`summary.json` `answer_tokens`):

  | Model | RIGA train | LIDC train |
  |---|---|---|
  | SmolVLM-256M | 316 / 504 / 2,268 | 103 / 407 / 1,186 |
  | LFM2.5-VL-450M | 207 / 324 / 1,427 | 70 / 260 / 750 |
  | Qwen3-VL-2B | 309 / 497 / 2,261 | 100 / 404 / 1,181 |

  On the pre-fix store, 9 of the 10 over-budget RIGA answers were BinRushed1 reconstructions (M7 §Decisions, I5).
  Those items now have no model image and are not in the prepared set; 1 over-budget answer per model remains there
  (S4 §Cross-path reports).

**Eval demo** `data/eval/demo` (`demo.md`): regenerated on the final store (LIDC val 1,084, RIGA val 63 items) with
the R2 code. The raters' own targets score IoU_tc 1.000 (LIDC) and 0.992 (RIGA) and fire nothing. Every pure
adversary fires its check. Holes-stamped and minority-hole mixed with honest rollouts fire nothing (I16).

**Pinned stack** (installed in `.venv`, recorded in the manifests): Python 3.12.11; numpy 2.5.4; scipy 1.18.1;
opencv-python-headless 5.0.0.93; numba 0.68.0; pyarrow 26.0.0; pillow 13.0.0; torch 2.14.1; transformers 5.20.1;
trl 1.16.0; datasets 5.1.0; peft 0.21.2; accelerate 1.15.0; pylidc 0.2.3; pydicom 3.0.2; nanobind 3.1.0;
scikit-build-core 1.1.1; ruff 0.17.1 (`pyproject.toml`; EXECUTION §4).

Changing numpy, scipy or OpenCV changes the reward fingerprint. Loading the caches then raises `StaleCacheError`
until `prepare_rl` is re-run (F2 §B2-8).

## 4. Break-it and fixes

- **Break-it (wave 5):** 136 attacks, 33 breakdowns.
  - Severity: 2 critical, 2 high, 11 medium, 18 low (ISSUES.md break-it table).
  - Notes: `impl/breakit/B1_geometry.md`, `B2_reward.md`, `B3_data.md`, `B4_wiring.md`. Ledger: `impl/ISSUES.md`.
- **Fixes (wave 6):** test-first, each recorded with its regression test:

  | Fix | Scope | Findings fixed |
  |---|---|---|
  | F1 | geometry, topology, targets | B1-1…6, B2-10, I1, I4, I8, I14 |
  | F2 | reward, eval | B2-4, B2-6…9, B4-5, B4-4 (eval side), I17/I18 definitions |
  | F3 | train | B2-1…3, B2-5, B2-8 (loader side), B3-5, B4-3, B4-4 / I21, the bf16-on-CPU follow-up |
  | F4 | data | B3-1…4, B3-6…8, I10, I20 |
  | F5 | CI, docs | B4-1, B4-2, B4-6…9, I3 and I19 (documented) |

- **The two critical findings:**
  - B4-1: a CI action tag did not exist. Every action is now pinned by commit SHA (F5).
  - B4-2: CI lacked pyarrow. It is now a core dependency (F5).
- **The two high findings:**
  - B2-1: a hacked checkpoint could become the rollback target (F3).
  - B3-1: pen strokes leaked the answer in BinRushed1 images (F4).
- **Markers:** no `xfail` marker remains in `tests/` (grep, 2026-10-19). Every reproducer is now a regression test, or was replaced by one (F4 §Reproducers).
- **Sweep (wave 7, S1–S4):** behaviour-preserving, with goldens byte-identical before and after.
  - S2: reward/eval goldens and the gate tables are identical (−17 lines).
  - S3: store rebuild byte-identical (−64 lines).
  - S4: prepared set byte-identical (−29 lines, then the voting-rater follow-up).
  - S1: all six golden sections identical; 6,021 → 5,896 lines (−125; the sweep alone −169, then `min_areas` +17 and F3 +23); K1/K2 slow gates and the fuzz re-run after the C++ edit; `make cov` 92.98 % (S1 §Line counts, §Verification).
- **Ledger:** every row of `impl/ISSUES.md` is fixed, decided or accepted, or open with a Phase 2
  owner. The open rows are I34–I38 and I41–I45. The 33 break-it rows and the 12 code-review rows (RV-1 … RV-12) are
  all fixed or decided.

### Code review (wave 7)

An independent, read-only review of the swept code (`impl/REVIEW.md`) found 12 issues: 1 high,
2 medium, 9 low. Fix round R1–R4 (`impl/fixes/R1.md` … `R4.md`), each test-first:

| Finding | Fix | Regression test |
|---|---|---|
| #1 high: `fetch_lidc_pixels` used series UIDs from the untrusted XML as paths (`rmtree`) and inside the URL | `lidc.uid` enforces the DICOM UID grammar (≤ 64 chars) at the one XML entry point; the URL is quoted (R1) | `test_data_lidc.py::test_planted_series_uid_is_refused_by_every_xml_reader`, `::test_uid_refuses_everything_outside_the_dicom_grammar` |
| #2 medium: illegal answers counted toward the Phase 2 geometry gates | one legality rule, `reward.legal_gate`, shared by `advantages_v4`, `topology_valid` and eval; illegal = 0 on IoU_tc, F5, g_area, Dice-type scores (R2). Follow-up: the TRL-logged IoU_tc / F5 / g_area, and so the best-checkpoint metric, are legality-gated too (R3) | `test_eval.py::test_illegal_answers_score_zero_on_geometry`, `test_train_trl.py::test_an_illegal_high_iou_evaluation_never_becomes_best` |
| #3 medium: the controller's audit statistics never reached a log | stored in the `audit` event of `controller.jsonl` and logged through `transformers.Trainer.log` (R3) | `test_train_trl.py::test_audit_statistics_reach_the_logs_and_controller_jsonl` |
| #4 low: build report counted 94 RIGA test images (store has 95) | count `image is not None` on the built item (R1) | `test_data_build.py::test_report_counts_the_riga_images_the_store_gets` |
| #12 low: `store.load` did not check a row's split | rows must match their file's `(dataset, split)` (R1) | `test_data_store.py::test_load_refuses_a_row_of_another_file` |
| #5, #6, #7, #8, #10, #11 low | stale docstring counts; `References:` lists and two citations checked against full texts; one truncation rule; `splits.crossing`, `configs.MODELS`, `configs.lora` and `topology.bbox` reused; smoke GRPO statistics use the gate definitions; `eval_demo` uses `min_areas` (R1–R4) | per fix note |
| #9 low: polycore `score_batch` unused in production | decided: kept (§3 contract, K3 bench, Phase 8 CUDA target) (R4) | — |

Reward v4 stayed byte-identical through R1–R4: the S2/S1 goldens of caches, rollouts and harness were re-run (fixes/R2.md, R4.md §Goldens).

## 5. A/B results

Full table, with sources and D1's measurement method: `impl/AB.md`. One line each:

| AB | Winner | Result | Evidence |
|---|---|---|---|
| AB1 PIL vs polycore | polycore | polycore: 0 wrong pixels. PIL on 3,000 K1 shapes: 340 wrong px per shape (4.3 % of area); on LIDC val targets, 14.7 % of area. Batch of 256: PIL raster 16.3 ms vs polycore 8.2 ms; with validity, 502 ms vs 9.9 ms | AB.md (D1 run); M1 |
| AB2 Phase 0 / v2 / v4 | v4 | Phase 0 loses the in-band speck (52.8 %), is majority-biased (R4 −0.56) and gets K9 90.7 % on full LIDC. v2 fails realistic R4 (+0.111) and gets K9 96.2 % | AB.md; M6 §A/B decisions |
| AB3 hole-matching term | v4 (term rejected) | Raises the holes-stamped win rate (tc 66.4 %) or leaves it unchanged (all 62.7 %), and breaks R4 (1 px −0.130, realistic −0.225) | AB.md; M6 |
| AB4 g_area zero vs τ·P | zero tolerance | Held-out median 1.000 / 0.991 / 0.972 ≥ 0.9, so the switch rule does not fire; τ·P lets every hedge through | AB.md; M6 §AB4; F2 |
| AB5 lossless vs `minseg_hybrid` | hybrid | Lossless costs 3.5–3.7× the tokens on RIGA, and 13 % of RIGA val answers exceed 1,536 Qwen tokens (hybrid 0). LIDC texts are identical. Hybrid keeps K5 at 100 % and IoU ≥ 0.99 | AB.md (D1 run); M5 |
| AB6 scipy vs Spaghetti | Spaghetti | Equal on every mask. 1.2–3.6× faster than scipy with the same algorithm (the lab's 5–24× is not confirmed). M3's single-thread timings may have used 8 threads (I14) | AB.md; M3 §Sandbox |
| AB7 our vs TRL advantages | TRL aggregator path | max \|TRL − v4\| = 2.37e-7; dropped groups exactly 0; end to end ≤ 1e-5 | AB.md; M7 §Result |

## 6. x86-64

Note: `impl/M9_x86.md`, sections "Final run" and "Re-run after the review fixes". All runs are in
Docker linux/amd64 under Rosetta 2, with `--cpus 3 --memory 4g --network none`.

**Final run** (post-sweep code with F1–F3, 477 s): **ALL PASS**.
- Build, lint, slow K1/K2 at full N, ASan/UBSan fuzz, local wheel.
- 626 fast tests. The three IPP cache tests that failed in Phase A now pass.
- Release wheel `objdel-0.1.0-cp312-abi3-manylinux_2_27_x86_64.manylinux_2_28_x86_64.whl`: 57 tests, auditwheel, abi3audit.
- Hashes, Mac vs x86: 22 / 26 sections identical, plus 4 accepted (raw split scores, I33). Caches, reward terms, advantages, eval, frame images and PNGs are now identical.
- Hashes, x86 run 1 vs run 2: 26 / 26.
- The Rosetta out-of-memory test runs as an informational step only (an artefact of Rosetta's own mmap).

**Re-run after the review fixes R1–R4** (`--quick --no-wheels`): PASS on build, fast suite, hashes (same counts), fuzz
and local wheel.

**Phase A findings and their fixes:**

| Finding (Phase A) | Fix | Status |
|---|---|---|
| 1. OpenCV's IPP `distanceTransform` is not deterministic run to run on x86 | F1: scipy EDT → float32, bit-identical to the Mac output (S2) | adopted (I30) |
| 2. `excursion_moment` S_m differs by 1 ULP | F2: S_1..S_8 frozen as `float.fromhex` (S2) | adopted (I31) |
| 3. LIDC bicubic upscale differs by ±1 on 413 k px | F3: Pillow fixed-point bicubic (S1); LIDC prepared images regenerated, masks / targets / caches unchanged | adopted (I32) |
| 4. Raw descriptor / NCC / ORB scores differ; decisions identical | F4 (float64 descriptor norm) skipped: splits are computed once and pinned by the store manifest | accepted (I33) |

**Still open for Phase 2** (I34):
- the native run on the first Vast box: AVX-512 numpy kernels, and K3/K4 timing (§8);
- the CI jobs on `ubuntu-latest`.

## 7. Decisions and deviations from PLAN

| # | Decision | Reason | Evidence |
|---|---|---|---|
| 1 | Local-gradient gate re-registered after measurement: moves *away* from the rater consensus ≥ 95 % per rollout, plus every move's mean ΔA < 0 | A closeness reward must raise the advantage of moves toward the raters (64 % per rollout lower it); the original per-rollout form cannot hold | I12; M6 §The two gates; PLAN §9 |
| 2 | R4 realistic neutrality accepted, monitored | Borderline pro-minority lean, +0.03 to +0.07 depending on the pool draw. It works against collapse; v4 constants were frozen on this set | I13; F2 §R4, paired |
| 3 | Hole adversaries (R3 b) monitored, not penalized; AB3 rejected | The hole term made them win more and broke R4. Residual: invisible to the audit at realistic prevalence | I15, I16; M8 §Open issues |
| 4 | Sub-filter raters abstain from the topology vote; gate set 251 → 239 items | A ≤ 4 px speck farmed the topology bonus | B2-4; F2 §B2-4; PLAN §9 |
| 5 | Validity gate per object (`iou_min`) | The mean over objects let one bad object pass | B2-7; F2 |
| 6 | Targets: per-object tolerance ladder (step 0.25), certified on the size-filtered signature; `fastgeom` dropped | 19–21 % fewer vertices than the lab's defect-localized repair, at the same certification; pen-stroke masks otherwise fell to τ = 0.5 | M5 §Sandbox |
| 7 | Over-budget answers dropped from SFT only, budgets unchanged (1,536 / 1,024 LFM); SFT teaches only voting raters | 0.28 % of targets; raising the budget lengthens every runaway rollout. Non-voting targets would teach answers RL scores as no-topology | I5; M7 §Decisions; S4 §Follow-up 1 |
| 8 | No augmentation in Phase 1; flips must rebuild targets and caches from flipped masks | The top-left tie rule is not mirror-symmetric (99.8 % of RIGA objects) | I7; M7 §Decisions |
| 9 | Disagreement items oversampled ×4 in GRPO rows | 3.9 % prevalence; ×4 gives ≈ 14.6 % of prompts and a step without one in 0.6 % | M7 §Decisions |
| 10 | TRL chat-row export moved from M4 to M7 | SFT rows need M5 targets; GRPO rows need M6 caches | EXECUTION §2 |
| 11 | No GPU rental in Phase 1; x86 timing (K3/K4) moved to the Phase 2 entry gates | Decided 2026-10-19; Docker under Rosetta cannot judge timing | PLAN §5 M9, §7; M9 §Open issues |
| 12 | LIDC store: identical reader layers merged; frames showing another outlined nodule dropped; 10,921 → 10,625 items; all series fetched (998 with pixels) | B3-3 (420 items), B3-2 (114 frames); I9 | F4; ISSUES I9 |
| 13 | RIGA: BinRushed1 reconstructions have no model image; second near-duplicate test (ORB + RANSAC); split 592 / 63 / 95 (plan 590 / 65 / 95) | B3-1 pen-stroke leak; B3-8 same-eye leak; the val carve-out deals groups, so photo counts follow the grouping | F4 §B3-1, §B3-8, §I10 |
| 14 | Hack-audit `area_tail` relative to a baseline (+10 points); "newly rewarded failures" = fails pass or outside [A_min/1.3, 1.3·A_max] | The absolute lines fired on honest small-nodule output (15 %) and marked ≈ 40 % of honest rollouts as failures. Thresholds set after measurement on synthetic models | I17, I18; F2 §I17 / I18 |
| 15 | Parser input cap `MAX_CHARS` = 16,384 (`too_long`); only JSON whitespace stripped | Bounded rejection time; RFC 8259 | I4, B1-1; F1 |
| 16 | Constrained-decoding regex `pattern(names)` available but not adopted | It would remove format errors from RL and change the reward mix: a Phase 2 decision | M2 §Parser |
| 17 | Phase 0 scripts not ported; they run at `faaacc4` | Porting would change the code behind published numbers | I3; F5 §I3 |
| 18 | `src/objdel/models/` and the `models` extra deleted | Unused in Phase 1, the only 0 %-covered file | S4 §medsam.py |
| 19 | Gate scripts import `tests/oracles` | One harness for pytest gates and data gates; moving it would ship test code or duplicate it | I19; F5 §I19; S2 §I19 |
| 20 | bf16 only where the hardware supports it (`is_torch_bf16_gpu_available`) | TRL defaults to bf16 and fails on CPU-only Linux | F3 §Follow-up; S4 §Follow-up 2 |
| 21 | K12 run with LoRA and micro-batch 1 on MPS | Mac memory (6.2 GB footprint); full fine-tune and vLLM paths are tested only on tiny models | I22; M7 §K12 smoke |
| 22 | `pack` offsets int32; lab API trimmed (`score`, `ring_labels`, brute validate) | One dtype on the hot path; one implementation per operation | M1 §API |
| 23 | Cross-architecture determinism: M9 F1–F3 adopted (scipy EDT, frozen S_m, Pillow bicubic LIDC upscale); F4 skipped | Caches built on x86 were not reproducible run to run. F3 is a sanctioned behaviour change (LIDC model images only). Split scores never need recomputing | M9 §Coordinator decisions, §Final run; I30–I33; S1 §F3 |
| 24 | One legality rule (`reward.legal_gate`): eval scores illegal answers 0 on geometry; the TRL-logged geometry and the default best-checkpoint metric are gated the same way | REVIEW #2: illegal answers counted toward the Phase 2 geometry gates and could become `best` | fixes/R2.md §1; fixes/R3.md §Follow-up |
| 25 | RIGA on rented GPUs moved to Phase 2 planning | Licence terms; see §8 | I44; `docs/research/notes/Phase 2 SFT/A3_vast_gpus.md` §5 |

## 8. Phase 2 handoff

**Entry gates** (PLAN §7; checked on val after SFT; `objdel.eval.gates.phase2_entry` decides each dataset on its
own, B4-5):

| Gate | Requirement |
|---|---|
| Format | parse ≥ 98 %, legal ≥ 90 % |
| Geometry | median IoU_tc ≥ 0.85 (disc), 0.75 (cup), 0.60 (LIDC); ≥ 70 % of rollouts pass IoU_tc ≥ 0.5 |
| Topology | topology-any ≥ 90 % on unanimous items |
| GRPO signal | temperature 1, G = 8: duplicates ≤ 10 %; ≥ 50 % of groups with σ above the floor (needs K a multiple of 8) |
| Reward speed on x86 | K3/K4 re-measured on the first Vast box's CPU (256 rollouts ≤ 250 ms, p99 ≤ 5 ms) before any training step |
| Cost | per-step wall time on the 4090 for the smallest model (G = 8, 1,536 tokens) replaces the roadmap's estimate |
| Precision probe | synthetic-shape probe per model: boundary error ≤ 2 px on LIDC-sized objects, else drop the model from LIDC (`phase2_entry(..., datasets=("riga",))`) |

**First steps on the first Vast box:**
1. Time K3/K4 serially and with 8 threads (M6 §Open issues).
2. Run `.venv/bin/python -I -m experiments.phase1.bench.xarch_hashes --out vast.json`, then compare it with the Mac
   hashes using `--compare build/x86/hashes-mac.json vast.json`. If `reward.advantages` differs (AVX-512 SVML), set
   `NPY_DISABLE_CPU_FEATURES="X86_V4"` for the reward process, or accept bit-level differences (M9 §Residual risk).

**Open issues to watch:**

| Issue | What to watch |
|---|---|
| I13 | R4 pro-minority lean; track per-image topology TV/JS on disagreement items |
| I15 / I16 | Hole-stamping adversary wins 61.6 % and is not seen by the audit at realistic prevalence; a baseline-relative TV / retention rule at K = 64 is the next check once real samples exist |
| I17 / I18 | Audit thresholds (tail share, newly-rewarded-failure envelope, added flags) were set on synthetic models. Recalibrate on the first SFT checkpoint's samples |
| I22 / I42 | The full fine-tune (SmolVLM, LFM) and vLLM colocate (LFM, Qwen) paths are untested at scale. vLLM cannot share the pinned `.venv` (TRL notes below) |
| I34 | First GitHub Actions run; native x86 hash compare and K3/K4 timing on the first Vast box |
| I35, I37 | Audit thresholds set on synthetic models; crop margin, LIDC side, jitter, item weighting and constrained decoding chosen, not tuned |
| I43 | The frame gives the answer away. `targets.frame` centres the crop on the raters' union and sizes it from the union. So an image-blind centred disc passes the RIGA geometry gates (median IoU_tc disc 0.904 vs 0.85, cup 0.783 vs 0.75) and nearly LIDC's (0.588 vs 0.60). Under a random crop shift it collapses (LIDC 0.00, disc 0.32, cup 0.07). Owner: Phase 2 plan. Gate on shifted crops, require a paired win over the blind prior, and train with crop shift / scale jitter (`docs/research/notes/Phase 2 SFT/A4_eval_data.md` findings 1–3, S5–S6; `experiments/phase2/lab/A4/`). A4 also re-ran MedSAM on our frames (box from the raters' union): median IoU_tc LIDC 0.78, disc 0.94, cup 0.85. Phase 0's LIDC 0.33 was an artefact of its crop and scoring (A4 finding 4, S7) |
| I44 | RIGA licence (moved to Phase 2 planning). LIDC (CC BY 3.0) may go on any host. RIGA goes only on Secure Cloud (`datacenter=True`) hosts, and only with written permission from the RIGA depositor and the MESSIDOR partners. MESSIDOR-derived images are 460 of the 750 (`docs/research/notes/Phase 2 SFT/A3_vast_gpus.md` §5 Data governance). Phase 2 screening uses LIDC only |

**TRL notes** (F3 §For the coordinator; M7 §Decisions; S4 §Follow-up 2):
- Use `grpo_config(..., eval_strategy="steps", eval_steps=N, eval_on_start=True)` with an `eval_dataset`, plus
  `trainer.add_callback(Controller(trainer))`. The step-0 evaluation is the audit baseline, and rollback is on by
  default.
- `num_generations = num_generations_eval = 8`; anything else is refused (B2-5).
- `per_device_eval_batch_size` ≥ 8 per process.
- bf16 only where `is_torch_bf16_gpu_available()` is true (Ampere+ on Vast); `sample_eval` uses bf16 only on CUDA
  cards that support it.
- TRL's own `reward` metric reads 0 with the aggregator; dashboards use `rewards/<term>/mean` (PLAN §4.6).
- SmolVLM uses HF generation (TRL PR #6718 not merged). The planned vLLM colocate path for LFM and Qwen
  (M7 §Decisions) **cannot share the pinned `.venv`** (I42):
  - vLLM 0.32.0 pins `torch==2.13.0` and `transformers>=5.16.1,<5.20.0`; 0.31 pins transformers < 5.18;
  - TRL 1.16 supports vLLM ≤ 0.31;
  - the repo pins torch 2.14.1 and transformers 5.20.1.

  Phase 2 / 3 decision. Phase 2 eval sampling runs in its own environment; Phase 3 colocate needs a re-pinned stack
  or server mode (`docs/research/notes/Phase 2 SFT/A5_decoding_serving.md` recommendation 2; `A3_vast_gpus.md`
  recommendation 5).

**Rebuilding `data/` from scratch** (Mac, `-I` because data is untrusted):

| # | Step | Command | Time (from the notes) |
|---|---|---|---|
| 0 | Raw inputs, by hand | RIGA from Deep Blue in a browser (`experiments/phase0/RESULTS.md`) → `data/raw/riga`; LIDC XML → `data/raw/lidc_xml`; pylidc's bundled DB is installed with `[data]` | — |
| 1 | CT pixels | `.venv/bin/python -I -m experiments.phase1.fetch_lidc_pixels data/raw/lidc_xml data/lidc_slices` (`--refetch-incomplete` for gaps) | wall time not recorded; 784 series, 0 failures, 5.0 GB (I9) |
| 2 | Store | `.venv/bin/python -I -m experiments.phase1.build_data --out data/store --workers 2` | 827 s (`build_report.json`); 1,578 s under load (S3); peak RSS 1.74–1.85 GB |
| 3 | Prepared set | `.venv/bin/python -I -m experiments.phase1.prepare_rl --store data/store --out data/prepared --procs 2` | 302 s, peak 1.2 GB (S4 §Golden comparison) |
| 4 | Targets | `.venv/bin/python -I -m experiments.phase1.build_targets --store data/store --out data/targets --procs 3` | 115.7 s (`summary.json`) |

Check the outputs against §3:
- store manifest `37fe7665…`;
- prepared jsonl hashes (LIDC `c19b1909…`, `bf73ac0d…`);
- `target_set_sha256` `1af819f4…`.

**Vast prerequisites:**
- The Vast AI account must be funded before Phase 2 starts (coordinator, 2026-10-19).
- The first session measures step time, and the measurement replaces the roadmap's GPU-hour estimate (ROADMAP §Compute).
- Use the smallest VLM that passes the entry gates: SmolVLM-256M → LFM2.5-VL-450M → Qwen3-VL-2B (PLAN D12).

## 9. Reproduce

```bash
make setup      # .venv (Python 3.12, uv), editable install with [train,data,dev]; builds the C++ extension
make test       # fast suite (not slow)
make test-all   # every test, including the slow K1/K2 10k gates
make cov        # fast suite with coverage of src/objdel; fails below 90 % (K13)
make lint       # ruff on src, tests, experiments
make fuzz       # native polycore fuzz under ASan/UBSan
experiments/phase1/bench/x86.sh            # x86-64 in Docker: build, lint, tests, K1/K2, fuzz, wheels, hashes (~13 min)
experiments/phase1/bench/x86.sh --quick    # without the 10k gates; 2,000-shape fuzz
.venv/bin/python -I -m experiments.phase1.reward_gates            # gate tables (941–3,147 s, 3 workers)
HF_HUB_OFFLINE=1 .venv/bin/python -I -m experiments.phase1.smoke_grpo --grpo-steps 3 --out data/runs/smoke_grpo_final   # K12 (787 s, Mac)
```

**Final verification (Mac):**

| Command | Result |
|---|---|
| `make lint` | clean (re-checked 2026-10-19) |
| `make cov` | 705 passed, 3 skipped, 0 failed; coverage **93.03 %** (2,239 statements); 14 min |
| fast suite after R3's last fix | 707 passed, 3 skipped, 0 failed (fixes/R3.md §Follow-up) |
| full `pytest` with the slow gates, before the R1–R4 fixes | 676 passed, 3 skipped, 0 failed; 47 min. The slow tests are the K1/K2 polycore gates, and the C++ is unchanged since |
| x86-64 Docker | ALL PASS (§6) |

That full run also printed a 76 % coverage reading. It came from files edited during the run and is not a
measurement; `make cov` above is the coverage of record.
