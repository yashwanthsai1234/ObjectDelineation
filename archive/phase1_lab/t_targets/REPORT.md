# t_targets: SFT target builder (mask → compact, legal, canonical polygon text)

Saved by the coordinator from the agent's final message; the agent could not write this file.
Code: `targets/` (fastgeom, simplify, builder, crops). Tests: `tests/` (67 passing). Tables: `results/*_table.md`.

## Chosen configuration

`BuildConfig(method="minseg_hybrid", tau=3, iou_min=0.99, repair="best", max_level=3, floor="exact")`.

**Pipeline**
1. Mask → exact crack-edge rings (from `objdel.polygonize`) → a unit-step lattice path.
2. Simplify with the optimal minimum-segment cycle (potrace polygon stage, vertices on pixel corners). A segment is admissible if it is either lossless (no pixel centre changes side) or within L∞ ≤ 3 grid units. Admissibility is tested in the rounded 0–1000 frame.
3. Round to the grid (exact round-half-up).
4. Remove duplicate and collinear vertices, and bevel pinches.
5. Certify:
   - legality rule L, per object;
   - exact re-rasterization (oracle rule);
   - (β0, β1) per object and per derived object (rim);
   - IoU ≥ 0.99.
6. Repair if needed. "best" = defect-localized pinning, or per-ring tolerance halving, whichever gives fewer vertices. Fallback: exact rings.
7. Canonicalize and serialize as compact JSON.

## KPIs (4,369 objects + 1,200 derived rims)

| KPI | Target | Result | Verdict |
|---|---|---|---|
| K5: topology preserved | 100 % | 100 % (RIGA rims 98.04 % before repair) | PASS |
| K6: IoU | min ≥ 0.98, median ≥ 0.99 | min 0.990; medians disc 0.9930, cup 0.9946, nodule 1.0 | PASS |
| K7: Qwen tokens per object, median | ≤ 250 | disc 157, cup 197, nodule 92 (LFM2.5 about 37 % fewer) | PASS |
| Legality | 100 % | 100 % | PASS |
| Determinism | identical text | identical across processes | PASS |
| Idempotence | rebuild gives same text | lossless mode 100 %; chosen mode LIDC 100 %, RIGA 2 % | lossless mode only |

**Answer length.** Qwen tokens ≈ 7.86 × vertices + 11 per object.

| Dataset | p50 | p95 | p99 | max |
|---|---|---|---|---|
| RIGA item | 353 | 481 | 665 | 1,017 |
| LIDC item | 92 | 212 | 300 | 643 |

So `max_new_tokens` = 1,536 for Qwen and SmolVLM, 1,024 for LFM2.5.

**Time per target:** RIGA median 66 ms (max 0.6 s); LIDC about 1 ms.

**Method comparison.** The hybrid beats the alternatives (`results/` has the full table):
- Douglas–Peucker and Visvalingam–Whyatt use more vertices at lower IoU. With a tolerance in grid units they do nothing on small LIDC crops (30 vertices vs 10).
- Classic potrace at ½ px costs 4–5× the tokens on large objects.

## Policies

**Pinches: bevel.**
- Each pinch vertex is split into two points one grid unit along its two edges.
- The result is legal and changes no pixel for crops < 1,000 px.
- Needed by 9.9 % of RIGA items before simplification and 12 / 3,169 after.

**Crop policy (H6).**
- Crop side must be ≤ 500 px at the resolution the reward uses.
- Prefer sides that divide 1000 (50, 100, 125, 200, 250, 500), reached by growing the context window, not by resampling.
- Larger objects are downscaled by an integer factor with block majority, and the downscaled mask is the ground truth for both SFT and reward. This changed 7.8 % of rim topologies.
- Never use a crop side ≥ 1,000 px.
- Simplifying in pixel space and then rounding flips pixels on every crop whose side doesn't divide 1000 (LIDC IoU down to 0.874).

## Risks

- **Lossy targets are not idempotent on large objects.** Use the lossless mode (3.3× the tokens) if that matters.
- **Thin rims produce the longest answers.**
- **The target style changes with crop size:** the lossless branch dominates on small crops. Run an SFT ablation.
- **Two floating-point tie-breaks.** Pin library versions and hash the targets.
- **The reward must use the same legality rule L** (per object) and the same downscaled masks.
- **Admissibility is a heuristic.** Exactness comes from the post-check, so keep the post-check in any port.
- **numba becomes a dependency.** Alternatively, port to C++ next to polycore.

## Promotion

Proposed locations: `src/objdel/targets/` (and `geometry/fast.py`), `tests/targets/`, `experiments/phase1/targets/`.

API: `build_target(masks, schema, BuildConfig(), derived) -> (text, objs, info)`.
