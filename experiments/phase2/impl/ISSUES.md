# Phase 2 Stage 1 issues

Break-it findings (wave 2). Format: id, severity, component, repro, expected vs actual, owner.

## Seeded by the coordinator from wave-1 reports

| id | sev | component | issue | owner |
|---|---|---|---|---|
| S1-1 | medium | M1 views | Under scale jitter, 11 of 24,701 LIDC rows (side ≥ 100) are certified but not pixel-exact against their mask (plan says 100 %); root cause unknown | V |
| S1-2 | low | M1 / M2 | `build_eval_views` calls the private `dataset._prepare(job)`; a signature change in T's file silently breaks it | V / T |
| S1-3 | low | M4b | Only 71 lung-background patches on the probe-test patient half (B3 used 200); topology probe draws on 19 val patients | P |
| S1-4 | low | M4b | MedSAM / SAM 2.1-tiny baseline paths never executed (no downloads in Stage 1) | P (Mac baseline step) |
| S1-5 | info | M7 | Unverified until V0: self-destroy endpoint, `VAST_API_KEY`, `destroy -y`, invoices format, onstart length limit | O (V0) |
| S1-6 | info | M1 | RIGA brightness/contrast/gamma jitter not implemented (Phase 2 runs are LIDC-only) | V (if RIGA permitted) |
| S1-7 | medium | M2 / M8 | `view_rows` has no `epoch` column; train_sft parses it from `row_id` (contract gap, `dataset.py:278`) | T |
| S1-8 | medium | M1 / M8 | `build_views` defaults to `--epochs 1`, R1 needs 3 (2-epoch branch decays into epoch 3); default and configs must agree | V / R |
| S1-9 | low | M4 / M8 | `vision.check` leaves the model in eval mode; caller must restore train mode | E |
| S1-10 | low | M2 / M8 | `configs.sampling` uses vLLM field names; g_probe maps to GRPOConfig names — one mapping should live in configs | T |
| S1-11 | medium | M3 | macOS: a pool worker importing a module after fork crashes the pool; imports must precede `fork_pool` (undocumented, untested) | S |
| S1-12 | medium | M8 | Not done: K = 16 + HM-IoU16 on the selected branch, rotated-view attribution check, WiSE-FT report, automatic D12(c) repair loop; `--items` takes the first N not a seeded subset | R |
| S1-13 | high | M2 | Qwen3.5 full FT loads bf16 base weights: no fp32 master weights for AdamW (updates below bf16 resolution are lost). Load fp32 + `bf16=True` autocast for full FT; bf16 base only for LoRA | T |
| S1-14 | medium | tests | `tests/test_breakit_reward.py::prepared()` writes no `manifest.json`; 2 tests fail after the I45 manifest-driven `cache_loader` | T |
| S1-15 | medium | M4 | `tests/test_p2_eval.py::test_nsd_equals_the_deepmind_reference` fails in the full suite (passes alone? order / randomness) | E |
| S1-16 | low | M5 | hf_teacher logprob parity (H-A5-2) not promoted; R0 needs it | T / R |

## Break-it round (wave 2): BK1 7, BK2 12, BK3 16, BK4 20 findings

Details, repro tests and suggested fixes: `impl/BK1.md` … `BK4.md`; strict-xfail tests in
`tests/test_p2_breakit_{data,scoring,runners,ops}.py`. High: BK2-1, BK2-2, BK3-1, BK3-3, BK3-4, BK4-1 … BK4-5.

## Fix round (wave 3): assignment by file owner

| Fixer | Files | Issues |
|---|---|---|
| FX-train | `train_sft.py`, `configs/*.yaml` | BK3-1 (yaml side), BK3-2, BK3-3, BK3-4 (yaml), BK3-5, BK3-8, BK3-10, BK3-12, BK3-13, BK3-14, BK4-5, BK4-6, S1-8 |
| FX-select | `select_ckpt.py`, `g_probe.py`, new `hf_sample.py` | BK2-1, BK2-2, BK2-5, BK2-6, BK3-4 (select), BK3-6, BK3-9, BK3-11, BK3-15, BK3-16, S1-12, BK3 simplicity cuts |
| FX-T | `train/{configs,dataset,sample}.py`, `sample_vllm.py`, T's tests, `test_breakit_reward.py` fixture | BK3-1 (sft_config), BK2-3 / BK3-7, BK1-3 (check side), BK1-4, BK1-5, BK4-1, S1-7, S1-10, S1-13, S1-14, S1-16 |
| FX-V | `targets/*`, `build_views.py`, `build_eval_views.py` | BK1-1, BK1-3 (write side), BK1-6, BK4-2, S1-2, S1-8 (default), BK1 simplicity cuts |
| FX-O | `vast/**`, `docker/**` | BK4-3 … BK4-19, BK1-2, BK1-7, BK4 simplicity cuts |
| FX-EGS | `eval/{gates,metrics,vision,signal}.py`, `train/serve.py`, `prereg.py`, `c2_analyze.py`, `eval/probes.py` (dedupe only), `tests/test_data_lidc.py` | BK2-4, BK2-7 … BK2-12, S1-9, S1-11, S1-15, BK4-20, c2/probes bootstrap dedupe |

## Status after the fix round (2026-10-20)

All 55 break-it findings and the seeded S1-* items are fixed, each with a test that failed first (no strict xfail
remains in `tests/test_p2_breakit_*.py`), except:

| id | state | reason / where it goes |
|---|---|---|
| S1-2 | open (low) | `build_eval_views` still calls `dataset._prepare` on one line; needs public `dataset.prepare_record` / `set_summary` helpers |
| WiSE-FT (S1-12 part) | deferred | report-only, runs only when the D12(c) repair fires; needs a run-chain stage after g_probe |
| S1-3, S1-4 | open | probe patch supply (71 vs 200); MedSAM / SAM 2.1-tiny paths first run in the Mac baseline step before V1 |
| S1-5, BK4-18, BK4-19 | UNVERIFIED until V0 | Vast self-destroy endpoint, `VAST_API_KEY`, invoice format, bandwidth field names, /dev/shm size |
| S1-6 | deferred | RIGA brightness / contrast / gamma jitter (LIDC-only until RIGA permission) |
| policy.py size | open (low) | 524 lines vs the ~380 target; the fixes added the safety checks |

Coverage config: `source_pkgs = ["objdel"]` added to `[tool.coverage.run]` (forked reward-pool workers otherwise traced
tests/ and experiments/, reporting 77 % for a library that is at 94 %).
