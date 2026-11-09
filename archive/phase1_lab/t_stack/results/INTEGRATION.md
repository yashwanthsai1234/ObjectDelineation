# Integration: R1, R2, K12 redo, R7 (t_stack, 2026-10-19)

Production path in `prod.py`, 1 module:
- strict serializer and parser
- polycore geometry
- reward-v2 terms
- advantages with a low-spread drop
- TRL aggregator

Components are vendored unedited into `vendor/`:
- polycore wheel, built from `../t_polycore` into `.venv` (`uv pip install ./vendor/polycore_src`)
- `targets` from `../t_targets`
- `rewardcore` from `../t_reward`, the iteration-2 snapshot

**IoU_tc fallback is a parameter.** `prod.terms(..., fallback="vote" | "max_all")` defaults to `"vote"`: IoU against the strict majority when no rater drew the sample's topology. It is implemented as a wrapper, so the rewardcore source is untouched.

Added environment packages: numba 0.68.0, opencv-python-headless 5.0.0.93, polycore 0.1.0. The rest is as in `versions.txt`.

## Summary

| Item | Result | Key numbers |
|---|---|---|
| R1 round trip | **PASS** | 6,338 targets (RIGA 1,200 × 2 modes, LIDC 1,969 × 2 modes) |
| R2 TRL aggregator equivalence | **PASS** | max \|TRL − advantages_v2\| = 8.7e-8 |
| K12 redo | **PASS (pipeline)**, with a caveat on the instructed checkpoint | see below |
| R7 precision feasibility | **analysis done**: RIGA 1–2 px not plausible (and not needed), LIDC plausible only if upscaled | see below |
| LFM licence | research use and publication allowed | LFM Open License v1.0 |

## R1: single-spec round trip (`r1_roundtrip.py`, results `r1_summary.json`, `r1_rows.jsonl`)

**Pipeline.** For each item:
1. Load the fixture item.
2. Crop with `pad_div`, the chosen H6 policy. RIGA crops become 500 px; LIDC crops are 25–100 px.
3. Build the target, per rater, in two modes: hybrid (`BuildConfig()`) and lossless (`minseg_exact`, iou_min = 1).
4. Serialize with the production format, so `"holes"` is always present.
5. Parse with the strict parser.
6. Run `polycore.validate` (roles, per object) and `polycore.rasterize` (even-odd, G = 1000, crop size).
7. Score with rewardcore v2 against the item's full rater set.

Run time was 50 s on 3 processes.

| | RIGA hybrid | RIGA lossless | LIDC hybrid | LIDC lossless |
|---|---|---|---|---|
| targets / objects | 1,200 / 2,400 | 1,200 / 2,400 | 1,969 / 1,969 | 1,969 / 1,969 |
| exact text → parts round trip | 100 % | 100 % | 100 % | 100 % |
| legal (polycore, 0 defects) | 100 % | 100 % | 100 % | 100 % |
| per-instance signature = source rater's | 100 % (rim 100 %) | 100 % | 100 % | 100 % |
| IoU_tc min / median | **0.9900** / 0.9941 | **1.0** / 1.0 | 1.0 / 1.0 | 1.0 / 1.0 |
| F5 ≥ 0.95 (min) | 100 % (0.985) | 100 % (1.0) | 100 % (1.0) | 100 % (1.0) |
| polycore raster == fastgeom raster, bit for bit | 100 % | 100 % | 100 % | 100 % |
| canonicalization invariance: rotate start; reverse then re-fix winding with a new start | 100 % (all terms and rasters identical) | 100 % | 100 % | 100 % |
| pure reversal flagged illegal by polycore (wrong_winding) | 100 % | 100 % | 100 % | 100 % |

**Notes.**
- The hybrid minimum IoU_tc sits exactly at the 0.99 certification floor, by construction.
- The RIGA rim (derived, not a target) has a minimum IoU_tc of 0.966. The rim was not part of the requirement.
- Grammar checks, all failures by design:
  - whitespace
  - a missing `"holes"` key
  - leading zeros
  - floats
  - values > 1000
  - code fences
  - trailing text
  - truncation, which also counts as a failure when the prefix parses
  - caps of > 1,000 vertices per ring or > 10 rings per object

## R2: TRL aggregator equivalence (`r2_trl_equiv.py`, `r2_trl_equiv.json`)

**Exact TRL 1.16 path.**
- `GRPOTrainer` on the tiny `trl-internal-testing/tiny-Qwen2ForCausalLM-2.5`, on CPU.
- 6 weight-0 loggers: `format_ok`, `legal`, `iou_tc`, `f5`, `topo_any`, `r_topo`.
- 1 aggregator, `v2_aggregate`, with weight 1.
- `sum_then_normalize`, `scale_rewards="none"`, `num_generations=8`, B = 4 × G = 8.
- Advantages are captured from `_generate_and_score_completions` and from `trainer._logs["advantages"]`.

**Groups are contiguous.** This was asserted inside the reward function. TRL presented the groups in the order 0, 1, 3, 2, and the reference uses the same order.

| check | run A (no drop) | run B (drop_spread = 1e-2) |
|---|---|---|
| max \|TRL adv − reference\| | 8.7e-8 (reference = `rewardcore.batch.advantages_v2`) | 8.4e-8 |
| all-identical group | 0 exactly | 0 exactly |
| low-spread group (range 5e-3, 1 topology) | – | 0 exactly. Without the drop it would get \|A\| up to 0.060 |
| kept groups vs `advantages_v2` on the kept groups only | – | 0.0 |
| illegal rollout | −5.60, the group minimum | −5.24, the group minimum |
| lone-minority-topology rollout vs the mean of the others | +0.59 | +0.71 |
| per-term logs vs hand-set vectors (per rollout, and `rewards/<t>/mean`) | ≤ 8.2e-8 | ≤ 8.2e-8 |
| NaN in logs | none | none |

**Drop rule as implemented.** A group is dropped when, for every geometric term (IoU_tc, F5, legal), max − min < 1e-2 **and** every rollout gives the same topology answer. Dropped groups get advantage 0 and are excluded from the batch-level z-score.

**Two caveats.**
- TRL's own `reward` metric now equals the mean advantage, which is 0. Monitor the per-term loggers instead.
- With more than one process, each reward call sees only the local slice, so the batch z-score must be gathered across processes. This is untested; the run here was single-process.

## K12 redo: real GRPO with the production reward (`k12_grpo.py`, `k12_*.json`)

**Setup.**
- SmolVLM-256M, HF generate (no vLLM), MPS.
- LoRA r = 16, β = 1e-3, `dr_grpo`, ε 0.2 / 0.28.
- `max_completion_length` = 1,536 (D10), `mask_truncated_completions=True`, G = 8.
- 2 prompts per generation, so the batch z-score covers 16 rollouts.
- The same 4 RIGA crops (448 px) the SFT used. The rater set is all experts' disc masks.
- Truncation is detected from `completion_ids`: no EOS means a format failure.

| run | steps | per-term logs | NaN | TRL vs aggregator adv | rewards | s/step, MPS |
|---|---|---|---|---|---|---|
| (a) instructed checkpoint `runs/sft_smolvlm-256m_full` (50-step overfit), T = 1.0 | 3 | all 7 present | no | 0 | **degenerate**: format 100 %, legal 100 %, IoU_tc 0.962–0.971, entropy 0.002. Every group is identical up to < 1e-2, so all groups were dropped and advantages are 0 | 46, 8.1 GB |
| (b) 15-step SFT checkpoint (loss 0.37), T = 1.0 | 4 | all 7 present | no | ≤ 1.2e-7 | **non-degenerate** (details below) | 209, 12.9 GB (swap) |

Run (b) details:
- format_ok 50–69 %, legal 6–25 %.
- IoU_tc 0–0.54, with group σ 0.02–0.21.
- 4–7 distinct topology answers per batch.
- 0 of 2 groups dropped (except step 4: frac_zero_std 0.5).
- |A| max 2.7–3.7.

Parse errors in run (b): int > 1000 ×9, expected ']' ×8, expected ',' ×6, expected int ×1.

**Conclusion.** The production pipeline runs end to end in TRL with correct logs and no NaN. An over-trained SFT checkpoint gives GRPO no signal, and the drop rule correctly zeroes it. **Phase 2 must stop SFT before collapse.** This argues for R7's "≥ 50 % of groups with σ above the floor at T = 1" gate.

`r_topo` was 0 throughout: no rollout was legal with IoU_tc ≥ 0.5 and a rater topology at the same time.

## R7: precision feasibility (`r7_resolution.py`, `r7_resolution.json`; processors only)

"Patch" is crop pixels per 16 px encoder patch. "Token" is crop pixels per LLM image token, along one side.

| model, setting | RIGA 500 px crop | LIDC 64 px | LIDC 128 px | LIDC 64→512 upscale | LIDC 128→512 upscale |
|---|---|---|---|---|---|
| SmolVLM-256M, longest_edge 512 | 64 tok; patch 15.6, token 62.5 | 64 tok (auto-upscaled to 512); patch 2.0, token 8 | patch 4.0, token 16 | same as 64 | same as 128 |
| SmolVLM-256M, longest_edge 1024 | 320 tok; patch 7.8, token 28 | patch 1.0, token 3.6 | patch 2.0, token 7.2 | – | – |
| LFM2.5-VL-450M | 256 tok; patch 15.6, token 31 | 64 tok (auto-upscaled to 256 only); patch 4.0, token 8 | patch 8.0, token 16 | 256 tok; patch 2.0, token 4 | patch 4.0, token 8 |
| Qwen3-VL-2B | 256 tok; patch 15.6, token 31 | 64 tok (256 px); patch 4.0, token 8 | patch 8.0, token 16 | 256 tok; patch 2.0, token 4 | patch 4.0, token 8 |

**How much precision the data needs.** This is the inter-rater P95 contour distance (τ_raw) at the `pad_div` crops in the fixtures:

| structure | crop side | object diameter (median) | τ_raw median [p10–p90] |
|---|---|---|---|
| RIGA disc | 500 px | 203 px | 8.9 [5.0–14.6] px |
| RIGA cup | 500 px | 112 px | 19.5 [9.5–31.5] px |
| LIDC nodule | 25–100 px (median 40) | 13 px | 2.0 [1.0–3.6] px |

**Judgement.** This is analytic, not a trained probe.
- **RIGA.** 1–2 px is 1/8–1/16 of an encoder patch and 1/30–1/60 of a SmolVLM@512 token, so it is **not plausible**. It is also **not needed**: raters disagree by about 9 px (disc) and 20 px (cup).
  - About ¼ patch (~4 px) is plausible for LFM and Qwen at 500 px.
  - For SmolVLM, use longest_edge 1024 for RIGA. That gives 320 tokens, a 7.8 px patch and a 28 px token. The 64-token setting compresses each 62 px square into one token.
- **LIDC.** τ is about 2 px, so 1–2 px matters.
  - **Upscale every LIDC crop to 512 px** (nearest or bicubic image upscaling; the mask stays at crop resolution).
  - That puts 1 crop px at 4–8 encoder px, so 1–2 px is ¼–1 patch, which is plausible.
  - LFM and Qwen only auto-upscale to 256 (64 tokens). Upscale explicitly to 512 (256 tokens).
  - SmolVLM already upscales to 512.
- **Recommendations.**
  - Crops: RIGA 500 px, as now; LIDC upscaled ×4–8 to 512.
  - Add a Phase 2 entry probe: SFT on synthetic shapes and measure boundary error in crop pixels against these settings.

## LFM licence

LFM2.5-VL-450M ships under the **LFM Open License v1.0** (https://huggingface.co/LiquidAI/LFM2.5-VL-450M/blob/main/LICENSE).
- **Allowed:** research and non-commercial use, with no restriction on publishing results.
- **Fine-tuned weights:** may be redistributed with the licence attached, NOTICE attribution retained, and modifications marked.
- **Commercial use:** not licensed for entities with annual revenue of $10M or more. Qualified non-profits are exempt for research.
- **Verdict:** OK for this project's research and publication. Recheck before any commercial use.

## Reproduction (in `t_stack/`, `export HF_HOME=$PWD/hf HF_HUB_OFFLINE=1`)

```
uv pip install --python .venv/bin/python ./vendor/polycore_src numba opencv-python-headless
.venv/bin/python -I r1_roundtrip.py 3
.venv/bin/python -I r2_trl_equiv.py
.venv/bin/python -I k12_grpo.py --temperature 1.0 --tag t10                                   # (a) overfit checkpoint
.venv/bin/python -I t4_sft.py --model smolvlm-256m --steps 15 --tag s15
.venv/bin/python -I k12_grpo.py --init runs/sft_smolvlm-256m_s15 --temperature 1.0 --steps 4 --tag s15_t10   # (b)
.venv/bin/python -I r7_resolution.py
```
