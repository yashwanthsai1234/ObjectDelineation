# K2: exact scoring and evaluation kernels on the device (design, 2026-11-08)

Folder: `docs/implementation2.0/kernels/K2/` (SP = the session scratchpad). REPO was not edited.

Tags used in this document:
- [measured: file]: I ran it in this session.
- [code: file:line]: read from the source.
- [paper] or [web]: from a cited source.
- [estimate]: not yet measured.

## Skills used

| Skill | What it changed |
|---|---|
| nvidia-skills:nvidia-skill-finder | I ran a live catalog check (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`, group "GPU Development"). The only relevant entry is `compute-sanitizer-diagnose-cuda-bugs`. It is recommended for the GPU run, not installed. Its tools are in GPU_RUN.md §4: memcheck, synccheck (the warp collectives), initcheck and racecheck. No catalog skill covers EDT, CCL or geometry kernels. |
| engineering:code-review | I self-reviewed for correctness, overflow, warp convergence, edge cases and resources. Results: (1) every warp collective sits behind warp-uniform control flow; (2) an i64/i128 overflow test case was added (coordinates near ±2^31); (3) a mutation test was written (§6). |
| engineering:debug | Reproduce → isolate → diagnose, applied to a future CPU/GPU mismatch. Both builds print one FNV-1a digest of every output. GPU_RUN §3 says how to bisect by section. The CPU build doubles as a bug oracle that a debugger can step through. |
| engineering:system-design | Requirements first: exact, single-source, and pixels never leave the device. That fixed the data flow: packed words on the device, then integers, then a few float64 values made by one thread. It also showed what stays on the host: the derived formulas, and HM-IoU only when its state space is too large. |
| engineering:architecture | The ADR-style choices are in §2. EDT: separable FH/Meijster rather than PBA+ or jump flooding. Pairs: uniform grid with a reference point rather than a sweep. HM-IoU: a quota DP rather than a GPU Hungarian. |
| engineering:documentation | Most useful facts first. GPU_RUN.md is a runbook with exact commands and pass rules. Facts are tagged, and links point to files instead of copying them. |
| huggingface-skills:huggingface-best | Loaded. No HF benchmark measures geometry kernels or EDT, so it has no model or leaderboard input. The kernel speed reference is the benchmark in this folder. |
| huggingface-skills:huggingface-papers | Used for metric definitions: GED (arXiv 1806.05034), HM-IoU (1905.13077, 2006.13144), HD-p pitfalls (2410.02630) and MeshMetrics (2509.05670). The papers in §10 that are not on arXiv were checked on the web. |

## 1. Problem

The checker and scorer already run on the GPU [code: REPO/cuda/gpucheck.cpp]. Four pieces were still slow or still on the CPU:

- **Validity totals.**
  - `ShapeTotals` runs one thread per shape [code: gpucheck.cpp:626-749]. When rings touch, it does O(R²·n²) serial work, up to about 10^8 operations on one thread [code: H1 F5].
  - `EdgePairs` is O(E²) per shape [code: gpucheck.cpp:585-623].
- **Evaluation metrics.** These run on the CPU with numpy and scipy, and they need the masks on the host [code: REPO/src/objdel/eval/metrics.py].
  - Pairwise IoU uses `_pairwise_iou` (:172).
  - GED uses `_ged` (:238).
  - HM-IoU uses `hm_iou` (:189, a scipy Hungarian).
  - NSD and EDT use `_surface` (:199, scipy EDT).
  - Probe HD95 is `probes.surface` [code: REPO/src/objdel/eval/probes.py:476].
  - The cache EDT is `cache.edt` [code: REPO/src/objdel/reward/cache.py:130].
  - Owner rule: GPU everywhere, and only scores leave the device.
- **Requirement: bit-exact.** Any thread schedule, any atomic order, and the CPU or GPU build must all give the same bits.

## 2. Algorithms (the sources are in k2score.inc)

All kernels follow gpucheck's single-source pattern: `__host__ __device__` functors, launched with `for_each` (one item
per thread). This work adds `for_each_warp` (one item per warp).

On the CPU build, `for_each_warp` runs a warp of one lane. Lane loops (`for i = lane; i < n; i += size`) become plain
loops, and the collectives (`wballot`, `wany`, `wodd`, `wsum`, `wmin`, `wsync`) become the identity. So one source
serves both builds, and the CPU build is the oracle that a debugger can step through.

**2.1 ShapeTotalsWarp: validity totals with one warp per shape.** This is a drop-in for `ShapeTotals`, with the same
outputs plus two new ones.

- Phase 1, one lane per ring: the ring-ring row sums and the wrong-winding flag.
- Phase 2, the whole warp per containment test:
  - `any_witness` with many witnesses: one lane per witness (vertices and edge midpoints). Each lane runs the serial
    `locate2`. The warp leaves at the first chunk of 32 that has a hit (`__any_sync`).
  - The single-witness test splits the ring's edges over the lanes:
    - on-boundary = `any(ballot(on))`;
    - inside = `popc(ballot(crossing parity)) & 1`.
  - Lane 0 writes the per-ring flags.
- Phase 3: each defect total is a warp sum (`__shfl_xor_sync`).
  - "Rings with a nonzero count" = Σ `popc(ballot(nz))`.
  - `clean` = every total is 0.
- `__syncwarp` separates the phases.

**2.2 GridBox → GridCount → GridScan → GridFill → GridPairs: candidate edge pairs from a uniform grid.**

- The grid has g × g cells per shape over the box of its clean edges. The host picks g = ceil(√E), capped at 128.
- The cell of coordinate v is (v − lo)·g / (hi − lo + 1) in i64. This map is monotone, so an edge's cell range covers
  its whole closed box.
- An edge goes into every cell that its closed box touches. Counts, offsets and fill use integer atomics. The cell order
  of entries does not matter.
- A pair whose closed boxes meet is classified only in the cell that holds the reference point (max xmin, max ymin):
  the lower-left corner of the boxes' intersection. That point lies in both boxes, so the cell is in both edges' cell
  ranges, and the pair is found exactly once [paper: Franklin / de Magalhães uniform grid].
- The pair is then classified by the same `classify<T>(earlier, later)` call, with the same counting, as `EdgePairs`.
  "Earlier" is the smaller vertex slot. Slots grow with (ring, k), which is `EdgePairs`' order.
- Shapes below `grid_min_edges` keep `EdgePairs` (`EdgePairsRest`). So do shapes whose entry count is above
  `budget_per_edge·E + 1024`; this caps memory on long diagonal edges.

**2.3 BoundaryRows, EdtRowsH, EdtColsV: boundary and exact squared EDT.**

- **Boundary.** It equals `cache.boundary` [code: cache.py:124]: a foreground pixel with a 4-neighbour in the
  background, where the outside counts as background. It is computed with word shifts: `b = cur & ~(up & dn & left &
  right)`, plus `popc` for the per-mask count.
- **EDT pass 1, one warp per row.** Lane b takes pixel 32j + b. The nearest source at or left of it is the top set bit
  at or below b (`clz`). If there is none, it is the last source of the earlier words, carried warp-uniformly. The
  right side works the same way with `ctz`. Stores are coalesced. This is the ballot/clz 1-D EDT idea of Manduhu and
  Jones [paper: TIP 2019].
- **EDT pass 2, one thread per column.** It uses Felzenszwalb–Huttenlocher's lower envelope of parabolas [paper: ToC
  8:415, Alg. 1], in Meijster's integer form: Sep(i, u) = floor((u² − i² + g(u)² − g(i)²) / (2(u − i))). There are no
  floats. The FH order of columns then rows also works ("changing the order … yields the same result", FH §2). The
  passes were swapped so that pass 2, the heavy one, has adjacent threads on adjacent columns: loads at one row and
  all stores are coalesced.
- **Empty masks.** g = H + W means "no source in this column". An all-H + W column means the mask is empty, and every
  value is then `kNoSource`.

**2.4 Directed, Select: boundary F and HD95 for mask pairs.**

- Each mask's boundary gets its own EDT, computed once per mask.
- Directed (one item per direction and row):
  - For every boundary pixel of src, d² = EDT of dst's boundary at that pixel.
  - hits += (d² ≤ τ²).
  - d² is appended to the direction's list (`atomic_inc`; order is free).
- Select (one warp per direction): the order statistics at ranks k = ⌊p_num·(n − 1)/p_den⌋ and min(k + 1, n − 1).
  - The k-th smallest value is the largest v with #(x < v) ≤ k, built bit by bit from bit 30 down. That is 31 passes,
    each a warp sum.
  - The value at rank k + 1 is the k-th value again if #(x ≤ v_k) ≥ k + 2, else the warp minimum of the values above
    v_k.
- The host turns the integers into floats:
  - HD-p = max over the two directions of √lo + (√hi − √lo)·frac, with frac = (p_num(n − 1) − k·p_den)/p_den, which is
    numpy's linear rule.
  - BF = 2PR/(P + R).
- Csurka's strict d < θ is the same as d² ≤ ⌈θ²⌉ − 1 [paper: BMVC 2013, eq. 3].

**2.5 SetInter, SetTerms: GED and HM-IoU over sets.**

- SetInter: one warp per pair (i ≤ j) of the set's K + k masks. It is Σ `popc(a & b)` over the frame's words, then a
  warp sum and one plain store. There are no atomics.
- SetTerms: one thread per set.
  - GED is metrics._ged's formula [code: metrics.py:238]: 2·E d(S, Y) − E d(S, S′) − E d(Y, Y′), with d = 1 − IoU.
    IoU(∅, ∅) = 1, and the result is clamped at 0. Sums run in row-major order.
  - HM-IoU is an exact dynamic program over the rater quota states. Each rater is matched ⌊K/k⌋ or ⌈K/k⌉ times
    (metrics.hm_iou's protocol [code: metrics.py:189]).
  - The state is a mixed-radix vector of counts. B[state] is the best total over the samples still to place.
    Successors are numerically larger states, so one backward sweep suffices.
  - The objective is Σ ⌊2^32·IoU⌋ in u64.
  - Reconstruction takes the smallest rater that keeps the optimum, so the result is the lexicographically smallest
    optimal assignment.
  - More than 2^16 states (for example k ≥ 17 raters with K = 1) gives choice −1, and the host falls back to scipy.

**Choices not taken (ADR summary)**

- PBA+ or cuCIM `distance_transform_edt` [web: D2 note F9].
  - For: published and fast.
  - Against: it returns nearest-site indices that can differ from scipy on ties, and it would be a new dependency in
    the RL image.
  - Our separable integer EDT outputs only the unique squared distance, so ties cannot matter.
- Jump flooding: approximate.
- A sweep-line (Bentley–Ottmann): serial, and it reports pairs, while polycore needs counts.
- GPU Hungarian: an 8 × 4 matrix is too small to need it [code: A1 note R5]. The DP is exact and has no float ties.

## 3. Complexity (N = pixels, E = clean edges of a shape, C = candidate pairs)

| Kernel | Work | Depth (per item) |
|---|---|---|
| ShapeTotalsWarp | that of ShapeTotals | O(R²·n²/32), was O(R²·n²) on one thread |
| Grid* | O(E + entries + C) per shape; entries ≤ budget·E | O(longest cell list) |
| BoundaryRows | O(N/32) word operations | O(W/32) |
| EdtRowsH + EdtColsV | O(N), FH Theorem 2.2 | O(W/32) + O(H) |
| Directed | O(N/32 + boundary pixels) | O(W/32 + boundary pixels in the row) |
| Select | 32·n per direction | O(31·n/32) |
| SetInter | O(M²·N/64) per set | O(N/32/32) |
| SetTerms | O(M² + S·k), S = (⌈K/k⌉ + 1)^k quota states | serial per set |

## 4. Memory and occupancy plan (RTX 5090, sm_120)

Device limits on sm_120 [web: Blackwell Tuning Guide, fetched 2026-11-08]:
- at most 48 warps per SM;
- 64K 32-bit registers per SM;
- 128 KB of shared memory per SM, 99 KB per block;
- at most 32 blocks per SM.

No K2 kernel uses shared memory. Registers per kernel from nvcc 13.0.88 `--resource-usage -arch=sm_120`
[measured: results/ptxas.txt]: every kernel has 0 spills and a 0-byte stack.

| Kernel | Regs | Warps/SM by registers |
|---|---|---|
| EdtRowsH (warp) | 38 | 48 (full) |
| EdtColsV | 38 | 48 |
| BoundaryRows | 40 | 48 |
| Directed | 32 | 48 |
| Select (warp) | 52 | 36 |
| SetInter (warp) | 40 | 48 |
| SetTerms | 64 | 32 |
| ShapeTotalsWarp (warp) | 56 | 36 |
| GridPairs | 98 | about 19 |
| EdgePairs, the existing kernel, for comparison | 80 | about 25 |
| ShapeTotals, the existing kernel, for comparison | 66 | about 28 |

All launches use 128 threads per block (4 warps), as gpucheck does.

Scratch memory:
- EDT: 16 B per pixel (g, s, t, d as i32). The 384 × 512² benchmark needs 1.6 GB, which fits in 32 GB. Production
  works on bounding-box crops [code: cache.py `_crop`], which are far smaller.
- Lists: 4 B per boundary pixel per direction.
- Grid: 4 B per cell plus 12 B per entry.
- HM: 8 B per quota state per set.

Data stays where gpucheck leaves it: packed words in the lane arena. Only integers and a few float64 values per object
are copied back. The test driver uses pageable copies and has one host round trip per call: the list sizes, and the
grid entry counts.

## 5. Exactness argument

1. **Every K2 output is one of three kinds:**
   - an integer;
   - a boolean;
   - a float64 computed by one thread from integers, in a fixed order, with FMA contraction off. Off means
     `--fmad=false` and `-ffp-contract=off`, as the repo's CMakeLists.txt sets for gpucheck.
2. **The integers are values, not orders.** Counts and sums are integer additions, which are associative, so atomics in
   any order give the same result. Other integer outputs are a minimum (EDT), an order statistic (Select) or the
   argmax with a fixed tie-break (the HM DP). Each is defined by the input alone, so any correct schedule gives it.
3. **The Directed list order is nondeterministic, but only order-free functions read it:** counts and order
   statistics.
4. **The grid finds each pair whose closed boxes meet exactly once** (the reference-point lemma, §2.2). Each such pair
   gets the same `classify` call as in `EdgePairs`. Pairs whose boxes do not meet are rejected by `EdgePairs`' own box
   test, so they add nothing there either. The integer counts are therefore equal.
5. **Warp collectives:** "any" and the parity are symmetric functions of the lanes' flags, so their result does not
   depend on how witnesses or edges are split among lanes.

**Against the current Python paths:**
- Integer outputs, the EDT and the order statistics are equal.
- The EDT matches `cache.edt` bit for bit after `float32(sqrt(double(d²)))`. scipy computes a float64 sqrt of the same
  unique integer and then casts to float32. So a GPU cache rebuild would not change the reward fingerprint (D1 §1.5).
- The GED and HM-IoU float64 values use our summation order. numpy's mean uses pairwise summation, so its values can
  differ in the last ulp.
- The HM assignment can differ from scipy's only when two assignments tie to within K·2^−32 (§8).

## 6. Verification done (CPU, no GPU spent)

- **k2test: exact equality on random and adversarial inputs** [measured: results/run_cpu_mac.txt].
  - Scale 1: 11,901 checks, 0 failures, digest `b158d8132b18be83`.
  - Scale 4: 47,547 checks, 0 failures, digest `f92b9a3939cec10f`.
  - The references:

    | What is checked | Reference |
    |---|---|
    | EDT | brute force min over sources |
    | Boundary | per-pixel 4-neighbour test |
    | Surface | sort-based percentiles |
    | Sets | brute-force enumeration of all k^K assignments, first strict maximum |
    | Validity | gpucheck's EdgePairs + ShapeTotals, and `gpucheck::check_polygons` end to end |

  - Validity runs three modes: all grid, mixed with budget fallback, and none.
  - All 13 defect classes occur. Shapes per class at scale 1 range from 17 to 53.
  - The adversarial cases:
    - masks: W ∈ {1, 31, 32, 33, 63, 65, 100}, empty, full, single pixel, checkerboard, lines, noise, discs with
      holes;
    - rings: repeats, closing repeats, spikes, pinches, collinear overlaps, rings with fewer than 3 points, touching
      strips, edges on grid seams, out-of-range coordinates, coordinates near ±2^31 (the i128 path);
    - sets: duplicate masks (exact HM ties), empty samples, and a set that is too large for the DP.
- **ASan + UBSan build**: passes with the same digest [measured].
- **The Linux g++ CPU build gives the same digest as macOS clang** [measured: docker, python:3.12-slim-trixie].
- **nvcc 13.0.88 builds and links the GPU build for sm_120**, with no warnings [measured: results/nvcc.log.txt, empty].
  With no driver, the GPU binary stops cleanly.
- **Mutation test**: 24 mutants of k2score.inc in three rounds. Two of them hit the earlier EDT layout, before the
  passes were swapped. Scripts: `results/mut_run.sh`, `results/mut_run2.py`.
  - 18 are killed. One of them was killed only after the ±2^31 test case was added.
  - The 6 survivors were checked by hand and are equivalent. They are:
    - swapped e/f order (classify is symmetric here);
    - the warp-locate on-boundary branch, which apart rings cannot reach;
    - the diagonal of the pair matrix, which is always 0;
    - `>=` instead of `>` in the DP maximum and in the envelope pop (both give the same value);
    - `<` instead of `!=` in the reconstruction (B[state] is the maximum).

## 7. Integration points (REPO, file:function)

1. **`cuda/gpucheck.cpp`**: paste `k2score.inc` into the unnamed namespace after `ShapeTotals` (:749).
2. **In `launch` (:1110), the validity kernels:**
   - Replace `for_each(st, c.n, ShapeTotals{…})` (:1124) with `for_each_warp(st, c.n, ShapeTotalsWarp{…})`.
     `clean` and `bad_rings` are optional extras.
   - Replace `for_each(st, c.nv, EdgePairs{…})` (:1123) with the grid chain, for shapes with E ≥ 512 input vertices.
     Make it without a round trip: give each shape a capacity of budget·E_s + 1024 entries; the host knows E_s, so it
     can plan this in `plan` (:1070). Then let GridScan set `use` on the device (entries ≤ capacity).
   - Skip GridBox and GridCount when no shape in the chunk uses the grid.
3. **`cuda/gpucheck.h` and `cuda/gpucheck_py.cpp`** (`m.def`, :113-213): new calls, which keep the masks on the device.
   - `surface_words` (BF and HD95 from the masks already in `L.dev` after `check_polygons`, as `score_words` does at
     :1549).
   - `set_terms` (inter matrix, GED, HM).
   - `edt_tables` (the cache EDT at `add_tables`).
4. **`src/objdel/eval/metrics.py`**:
   - `_pairwise_iou` (:172) → the SetInter inter matrix.
   - `_ged` (:238) and `hm_iou` (:189) → SetTerms. Fall back to scipy when choice = −1.
   - `metrics()` (:276) stops building the K × H × W masks on the host for these terms.
   - NSD (`_surface`, :199) can use the same EDT on the (H + 1) × (W + 1) corner-code border. The surfel lengths are
     integer multiples of 1 and √2/2, so the sums stay exact. This is not built yet.
5. **`src/objdel/eval/probes.py:476` `surface`**: HD95 → Select. ASSD and NSD@1 come from the same Directed lists.
6. **`src/objdel/reward/cache.py:130` `edt`**, `_contour_distances` (:168), `rater_tau` (:176) and `_band` (:185):
   - Use EdtRowsH and EdtColsV with sources = the complement.
   - `kNoSource` → `_BIG`.
   - This path is offline. The bits are identical (§5), so there is no fingerprint change.
7. **`experiments/phase2/final_eval.py:23` `TERMS` and `per_item` (:26)**: add bf, hd95, ged and hm_iou. GED is lower
   is better.
8. **Tests**: port the k2test families into `tests/test_cuda_checker_gpu.py` and `tests/test_gpu_scoring.py`, with
   GPU vs CPU-build equality, as the repo does now.

## 8. Risks

1. **No GPU numbers yet.** Correctness on sm_120 is not verified until GPU_RUN.md §3–4 pass. The compile is verified.
2. **The grid helps only the tail.**
   - CPU build, 512 typical star shapes of 64–200 vertices: grid 34.1 ms against EdgePairs 12.5 ms. Long radial edges
     cover many cells [measured: results/run_cpu_mac.txt].
   - Adversarial shapes of 10 rings × 1,000 vertices: pairs 295 → 24 ms [measured: CPU build].
   - So gate the grid at E ≥ 512 (D1 also advised against the grid at E ≲ 300).
   - A crowded cell, such as many edges through a star's centre, stays O(E²) inside that cell. The budget caps
     memory, not time.
3. **ShapeTotalsWarp splits witnesses only.** One straggler shape now takes O(n²/32) on one warp. If the GPU run shows
   it still dominates, split each witness's edges too, in a 2-D split.
4. **EdtColsV gathers g(s[q])** from up to 32 rows per warp. It relies on L1/L2 hits. Fix if needed: keep g(s[q])²
   beside s[q] (+4 B/px) or run the envelope warp-cooperatively (PBA-style bands).
5. **Float formulas differ from numpy/scipy in the last ulp** (GED, HM mean, HD-p lerp), and the HM argmax can differ
   on near-ties under 2^−32·K. Report them as K2's definitions, or make the host formulas copy numpy's pairwise sum
   and `_lerp`.
6. **The HM state space grows as (⌈K/k⌉ + 1)^k.** RIGA has k = 6; K = 16 gives 4,096 states, which is fine. Above
   2^16 states the device declines.
7. **Two host round trips per call in the test driver**, and many small launches. Fix in integration: the capacity
   plan (§7) and a device scan (CUB `DeviceScan`), then CUDA Graphs if launch time shows up in nsys.
8. **Low occupancy for GridPairs** at 98 registers (about 19 of 48 warps). It is a tail-only kernel.
   `__launch_bounds__` or an i32 orient tier (H1 F7) would cut it.

## 9. Results

- CPU build, macOS M-series, 1/8-size batches [measured: results/run_cpu_mac.txt]:
  - EDT: 81 Mpx/s.
  - Boundary F + HD95: 128 pairs in 162 ms, of which the EDT is 98 %.
  - GED + HM-IoU: 4 sets in 0.23 ms.
  - Adversarial validity: 310 → 39 ms, the same defects. ShapeTotalsWarp does not help on the CPU (1 lane).
- GPU: pending GPU_RUN.md. Record the numbers here, marked [measured: bench_gpu.txt].
  - Expectation [estimate]: EDT about 2–5 ms for 100 Mpx. That is about 16–24 B/px of traffic at about 1.8 TB/s.

## 10. Citations

- P. F. Felzenszwalb, D. P. Huttenlocher, "Distance Transforms of Sampled Functions", Theory of Computing 8(19):415–428,
  2012, doi:10.4086/toc.2012.v008a019 [web: theoryofcomputing.org/articles/v008a019].
- A. Meijster, J. B. T. M. Roerdink, W. H. Hesselink, "A general algorithm for computing distance transforms in linear
  time", Mathematical Morphology and its Applications to Image and Signal Processing, 2000 [from memory, not
  re-checked].
- M. Manduhu, M. W. Jones, "A Work Efficient Parallel Algorithm for Exact Euclidean Distance Transform", IEEE TIP 2019,
  doi:10.1109/TIP.2019.2916741 (1-D EDT with ballot/ffs/clz/shfl) [web].
- T.-T. Cao et al., Parallel Banding Algorithm (PBA+), github.com/orzzzjq/Parallel-Banding-Algorithm-plus [web].
- D. Black et al., "Accelerating Euclidean Distance Transforms…", IEEE Access 2025, doi:10.1109/ACCESS.2025.3548563
  [web].
- Y. Lu et al., "GPU-Based Parallel Euclidean Distance Transform Algorithm", Mathematics 14(4):597, 2026 [web].
- M. V. A. de Magalhães, W. R. Franklin et al.:
  - uniform grid and exact predicates on GPUs, CAD 150 (2022) 103285 (wrfranklin.org/p/239);
  - the 2020 Spatial Gem (wrfranklin.org/p/238) [web].
- G. Csurka, D. Larlus, F. Perronnin, "What is a good evaluation measure for semantic segmentation?", BMVC 2013,
  doi:10.5244/C.27.32 (BF score) [web].
- L. Maier-Hein et al., Metrics Reloaded, Nat. Methods 2024 (HD95 = max of the directed 95th percentiles) [web:
  metrics-reloaded.dkfz.de/metric-library/xhd]. Implementation pitfalls: arXiv 2410.02630. MeshMetrics: arXiv
  2509.05670.
- S. Kohl et al., Probabilistic U-Net (GED), arXiv 1806.05034. HPU-Net, arXiv 1905.13077. CAR (HM-IoU16), arXiv
  2006.13144.
- NVIDIA Blackwell Tuning Guide (cc 12.0 limits) [web: docs.nvidia.com/cuda/blackwell-tuning-guide].
- CUDA Programming Guide [web: docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/compute-capabilities.html]:
  - compute capabilities: warp reduce functions need cc ≥ 8.0;
  - warp vote and shuffle semantics, the CUDA Handbook ch. 8.

## Files

| Path | Role |
|---|---|
| `k2score.inc` | The kernels: a gpucheck.cpp fragment, 722 lines |
| `k2.cpp`, `k2.h` | Host driver and API. It includes `vendor/gpucheck.cpp` and then `k2score.inc` |
| `vendor/gpucheck.{cpp,h}` | An unmodified snapshot of REPO/cuda. SHA-256 values are in `vendor/SOURCE.txt` |
| `k2ref.cpp`, `k2ref.h` | The independent brute-force references |
| `k2gen.h` | Deterministic generators (splitmix64, integer-only geometry) |
| `k2test.cpp` | The exact-equality harness and digest |
| `k2bench.cpp` | The micro-benchmark |
| `run_cpu.sh` | The CPU build, tests, sanitizers and benchmark |
| `GPU_RUN.md` | The GPU runbook |
| `results/` | Logs: `run_cpu_mac.txt`, `ptxas.txt`, `nvcc.log.txt` |
