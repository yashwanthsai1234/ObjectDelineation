## A1. System (Phase 1 builds everything inside the dashed box; Phases 2+ run it on Vast)

```
                   ┌───────────────────────────── Phase 1 (Mac, CPU) ─────────────────────────────┐
 raw data          │                                                                              │
 RIGA (6 raters) ──┼─► data/ingest ──► exact expert masks ──► splits (group-aware, seeded)        │
 LIDC (4 readers)  │   (staged RIGA     (certified, 4,500       │                                 │
                   │    extractor,       RIGA + LIDC slices)    ▼                                 │
                   │    LIDC parser)              store (Parquet + SHA-256 manifest)              │
                   │                                  │                    │                      │
                   │                                  ▼                    ▼                      │
                   │                      targets (mask → canonical  reward cache (per prompt:   │
                   │                      polygon JSON, certified)   vote/union/core, EDTs, τ,   │
                   │                                  │              T5, rater topologies)       │
                   │                                  ▼                    │                      │
                   │                      exporters (TRL chat rows)        │                      │
                   └──────────────────────────────────┼────────────────────┼──────────────────────┘
                                                      ▼                    ▼
   Phase 2  SFT (TRL SFTTrainer) ──► SFT checkpoint ──► Phase 3+ GRPO (TRL GRPOTrainer)
                                                            │   rollouts (G = 8 per prompt)
                                                            ▼
                                     ┌──────────── reward server (objdel.reward) ────────────┐
                                     │ parse → polycore.validate → polycore.rasterize         │
                                     │ → topology (Spaghetti) → IoU_tc, F5, legality, r_topo  │
                                     │ → advantages (GDPO floor + gates, DM-loo mean-centred) │
                                     │ → typed error log                                      │
                                     └───────────────┬───────────────────────────┬────────────┘
                                                     ▼                           ▼
                                        Level 2 controller (callback):   evaluation harness
                                        learnability sampler, typed       (same checker, same
                                        alarms, rollback, β / weights     rasterizer, metrics)
```

## A2. One GRPO step through the reward (all exact, single spec)

```
 32 prompts × 8 completions (text)
   │
   ├─ parse JSON (schema per dataset) ─────────────── unparseable → format = 0 (bottom of every term)
   │
   ├─ polycore.validate (C++, sweep, int64) ───────── defect classes → legal ∈ {0,1}, graded counts
   ├─ polycore.rasterize (C++, winding, top-left) ─── even-odd mask (+ nonzero, mixed-sign flags)
   │
   ├─ topology: cv2 Spaghetti 8/4 + hole→component ─ per-instance (β0, holes) key
   │
   ├─ cached per prompt: rater masks, vote/union/core, EDTs, τ, T5, rater topology histogram f
   │
   ├─ terms:  IoU_tc = max IoU over raters with the same topology
   │          F5     = 1/(1 + (2D/P)/T5²)            (distance beyond the egg–yolk band)
   │          legal  = no defects
   │          r_topo = 1[valid]·(1 + ½(log f(a) − log π̂_loo(a)))
   │
   └─ advantages: legality gate → per-term z with σ floors (IoU, band) → batch z
                  + c_topo · mean-centred r_topo → drop groups with spread < 1e-2
                  → TRL (aggregator reward; per-term rewards logged with weight 0)
```

## A3. Package layout after Phase 1

```
cpp/polycore/          header-only C++17 core (rasterize, validate, signed_area2, score_batch)
src/objdel/
  _polycore            nanobind extension (built by scikit-build-core)
  checker/             bitquad, polygonize, certify, topology (+ fast instance topology), polycore wrapper
  geometry/            canonical polygon schema, parse/serialize
  targets/             SFT target builder (minseg_hybrid, repair, bevel, crop policy)
  reward/              cache, terms (IoU_tc, F5, legality, r_topo), advantages, batch scorer, typed errors
  data/                riga (staged extractor), lidc (parser, patient ids), splits, store, exporters
  train/               TRL adapters: reward functions, aggregator, controller callback, configs
  eval/                metrics harness (same checker), reports
tests/                 unit + property + oracle tests (oracles vendored in tests/oracles/)
experiments/phase1/    scripts + results + RESULTS.md
cuda/                  original engine; later includes cpp/polycore
```
