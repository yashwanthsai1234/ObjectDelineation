# ADR-0003: Exact integer geometry on a 0-1000 grid

**Status:** Accepted
**Date:** 2026-10-19 (commit `1c54bc4`: polycore, text grammar and targets; certificates since `29868b9`, 2026-10-18)
**Deciders:** Repository owner

## Context

The model writes polygons. The reward needs to know if a polygon is legal and what pixels it covers.
Floating-point geometry has edge cases: a vertex that sits exactly on an edge, a collinear overlap, a pixel centre on a polygon edge.
An RL policy finds such cases. A checker that answers differently for them gives the policy a way to hack the reward.

## Decision

All checking uses integer arithmetic.

- Outline coordinates are integers on a 0 to 1000 grid. Grid point (x, y) is pixel point (x W / 1000, y H / 1000).
- Each shape runs in int64 or in `__int128`. The choice is made per shape from fixed bounds (`Precision` in `cuda/gpucheck.cpp`).
- A pixel is inside when its centre is inside under a half-open crossing rule. Centres on left and top edges are inside. Centres on right and bottom edges are outside.
- A ring is legal only if none of 13 named defects is present.
- Frame sides divide 1000. Then one pixel is a whole number of grid units and the round trip mask, polygon, mask has no rounding error.
- Masks are certified three ways: bit quads, labelled topology and traced rings must agree on chi8, ring counts, area and corners (`objdel.checker.certify`).
- Device float64 terms use no FMA contraction, so each term is the IEEE operation of the host formula.

## Options considered

### Option A: Exact integers (chosen)

| Dimension | Assessment |
|---|---|
| Complexity | Medium. Needs wide integers and a defect taxonomy |
| Cost | `__int128` on the GPU is slower than int64. The cost is accepted |
| Reproducibility | Identical on every CPU and GPU |

### Option B: Floating-point polygon clipping (shapely-style)

**Pros:** common libraries. **Cons:** tolerance choices. Results depend on library and architecture. Not bit-reproducible.

### Option C: Rasterise and compare only

**Pros:** simple. **Cons:** a self-crossing ring still rasterises. Legality would be undefined, and the model could write degenerate rings that pass.

## Trade-off analysis

The reward is the specification of correct output. Any ambiguity in it is a hole in the specification.
Integers close that hole at a modest speed cost.
The division rule for frame sides costs flexibility in frame size, and it removes a measured loss: crops whose side did not divide 1000 flipped pixels and dropped LIDC IoU to 0.874 (`objdel.targets.frame`, lab H6).

## Consequences

- Easier: tests use exact oracles (`tests/oracles`) and demand equality, not tolerance.
- Easier: the same polygons, masks and defects on every CPU and GPU tested: arm64 CPU, sm_86, sm_120 and RTX 4090.
- Harder: frames are limited to side <= 500 and a divisor of 1000: 50, 100, 125, 200, 250, 500.
- Harder: the reward tables hold distances that are 0 or in [1, 2^24], which the binding checks before upload.
- Revisit: pinned library versions (numpy, scipy, OpenCV, numba) are part of the reward fingerprint.

## Evidence

| Check | Result | Source |
|---|---|---|
| Certificates on real masks | 0 failures on 26,820 RIGA masks and 40,189 LIDC slices | `docs/research/records/phase0/RESULTS.md` |
| Raster and legality vs slow oracles | 0 of 849,802,293 pixels differ on 10,000 shapes; 0 of 10,000 rings differ | `experiments/phase1/RESULTS.md` (K1, K2) |
| Defect classes | All 13 flagged exactly on 8,559 fuzz polygons | `experiments/phase2/perf/p1_results.md` |

## Action items

1. [x] Exact oracles for raster, validity and topology (`tests/oracles`).
2. [x] Slow 10,000-case gate in CI (`gates` job).
3. [ ] Extend frames above 500 pixels. The target builder refuses them today.
