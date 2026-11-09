# t_polycore: exact C++ polygon core (header-only C++17 + nanobind)

Saved by the coordinator from the agent's final message; the agent could not write this file. Raw results are in `results/`.
- **Core:** `cpp/polycore/polycore.hpp` provides `rasterize`, `validate`, `signed_area2` and `score_batch`.
- **Bindings:** `src/_polycore.cpp` (nanobind, GIL released in every kernel), built with scikit-build-core.
- **Wrapper:** `python/polycore/`.
- **Tests:** 46 in `tests/`.
- **Benchmark:** `bench/bench_k3k4.py`.
- **Linux check:** `docker/`.

## KPIs (M3, load average 50–150 from other labs, so timings are pessimistic)

| KPI | Target | Measured | Verdict |
|---|---|---|---|
| K1: raster equals the oracle (both rules, per-ring winding flags) | 100 % | 0 / 849.8 M pixels wrong on 10,000 shapes (9 families incl. on-edge lattices, degenerate cases, holes in both orientations), against `oracle.py` itself; 10⁶ more natively | PASS |
| K2: labels equal `oracle.ring_problems` | 100 % | 0 / 10,000 rings; finer classes and multi-ring checks exact on 3,000 + 1,500 cases | PASS |
| K3: 256-rollout batch | ≤ 250 ms | 8.5 ms median, internal 8-thread pool (Python ThreadPoolExecutor: 14.8 ms) | PASS (~30× margin) |
| K4: p99 per rollout | ≤ 5 ms | 0.38 ms (p50 0.09 ms); worst case 3 × 200 vertices at 512²: 0.58 ms | PASS |
| K10: determinism | 100 % | identical SHA-256 across thread counts 0–16, Python pools and subprocesses | PASS |

Comparisons:

| | Time | vs polycore |
|---|---|---|
| repo pure-Python path | 46 s per batch serially | ~1,000× slower per rollout |
| numba prototype | 0.455 ms p50 per rollout | ~5× slower |

**Linux (Docker, python:3.12-slim, g++ 14):**
- arm64 native: 46/46 tests pass, native fuzz of 100k shapes clean under ASan+UBSan.
- amd64 emulated: 46/46 pass, 50k fuzz clean.

**macOS:** 1 M shapes clean under ASan+UBSan.

## Design

**Tie rule.** Identical to the oracle: pixel-centre sampling, half-open crossing, edges going down the screen count +1. This is the GPU top-left rule. It is not mirror-symmetric, so ground truth and predictions must both be rasterized by polycore.

**Scanline.** In the scaled frame (x·2W, y·2H) everything is an integer. Each edge contributes per row by exact floor division, then a row prefix sum gives the winding number, from which even-odd and nonzero come in one pass.
- Arithmetic is int64 when values are < 2³⁰, otherwise __int128 (exact for any int32 input).

**Flags.**
- `ring_wmin`/`ring_wmax` flag rings with mixed-sign winding, which catches bow-ties and figure-8s that the even-odd/nonzero XOR misses.
- `px_multi` counts pixels with |winding| ≥ 2.

**Defect classes (counts are edge pairs).**
- Per ring: fewer_than_3_points, zero_area, self_crossing, vertex_pinch, collinear_overlap, spike.
- Between rings: ring_ring_crossing, ring_ring_touching.
- Structure: hole_outside_outer, nested_holes, nested_outers, wrong_winding (outer must be clockwise on screen, hole counter-clockwise), out_of_range.

**Pair enumeration and containment.**
- Pairs are found by sort-and-sweep on x with a y-overlap test: 2.6× faster at 200 edges, 22× at 2,000.
- Containment fast path: when two rings neither touch nor cross, one witness vertex decides containment exactly.

**Build fixes found.**
- nanobind adds `-Os` after `-O3`, which roughly halves speed; `NOMINSIZE` is required.
- An O(n²) containment scan was replaced by the fast path above.

**API.**
- `pack`, `rasterize`, `validate(roles)`, `ring_labels`, `signed_area2`, `score`, `score_batch(threads=8)`.
- Bad input raises an exception and never crashes. W and H must be ≤ 4096; at most 1,024 rings.

## Proposed repo layout

| Path | Contents |
|---|---|
| `cpp/polycore/polycore.hpp` | the core, also used by `cuda/` |
| `src/objdel/_polycore.cpp` | nanobind module `objdel._polycore` |
| `src/objdel/checker/polycore.py` | wrapper; replaces PIL rasterization and the Python `ring_problems` |
| `tests/polycore/` | tests, with `oracle.py` vendored |

Build: switch to scikit-build-core, or ship a separate small wheel. CI per OS: build, sanitizer fuzz, pytest.

## Risks and open decisions

1. Re-run the benchmark on a Vast x86 host.
2. Counts are edge pairs (one pinch = 4), so the reward should use presence or a saturating transform.
3. Semantics to confirm:
   - wrong_winding uses net area;
   - a hole crossing its outer gets both crossing and outside flags.
4. 8-connected pinches in bit-quad SFT targets are flagged `vertex_pinch`. Targets must remove them, or the reward must allow them.
5. Flip augmentation of rasterized ground truth differs at on-edge centres; re-rasterize after the flip.
6. Cap input size: ≤ 1,000 vertices per ring, ≤ 10 rings.
7. CUDA port: keep int64 on the device and route inputs ≥ 2²⁸ to the host (or reject them).
8. Not done:
   - mark the module free-threaded;
   - avoid zeroing the output masks twice;
   - cibuildwheel configuration.
