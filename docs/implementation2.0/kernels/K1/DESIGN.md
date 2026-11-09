# K1: a work-optimal CUDA polygonizer for binary masks

Folder: `docs/implementation2.0/kernels/K1/` (SP = the session scratchpad). The REPO is not edited. Status: correct on the CPU
build. Its output is identical to `polygonize.py`, `bitquads` and `check_masks` on all 37,093 real LIDC masks (seed
rows 16 and 0), on 4,500 random masks, and on 4,453 random and adversarial masks run in 15 configurations and 3 item
orders (results/). nvcc 13.0 compiles and links it for sm_120. The GPU run (GPU_RUN.md) is
still to do.

## Skills used

| Skill | What it changed |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | Live catalog check (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`, "GPU Development" group). No skill covers CCL, contour tracing or list ranking. The only relevant one is `compute-sanitizer-diagnose-cuda-bugs`. It is recommended only, not installed, and its tools are GPU_RUN.md step 4 (memcheck, initcheck). racecheck and synccheck do not apply: K1 has no shared memory and no barriers. |
| `engineering:code-review` | Applied to K1 itself. 8 deliberate mutants (wrong connectivity, missing border union, no surround check, right turn first, corner miscount, seed-map off by one, dead-seed check off, first-write instead of min) are all caught. The last one survived forward CPU order, which led to `K1_ORDER` (items in reverse / permuted order on the CPU build). ASan + UBSan are clean. The worst-case memory review added the `max_scratch_bytes` halving guard. |
| `engineering:debug` | Defined the isolation path before any GPU spend: `seed_rows=0` separates tracing from CCL, an error flag holds the device-side invariants, and GPU_RUN.md "If something fails" lists the steps. |
| `engineering:system-design` | Requirements first (exact, batch, CPU = GPU, ≤ 512² frames, worst-case bounded), then data flow (section 2), trade-offs, and what to revisit (section 10). |
| `engineering:architecture` | ADR in section 3: dense slot grid vs edge list + list ranking vs CCL-derived rings + seeded tracing vs block-based CCL. |
| `engineering:documentation` | Reader-first layout. GPU_RUN.md is a runbook with exact commands, pass criteria and time per step. |
| `huggingface-skills:huggingface-best` | No HF benchmark measures CCL or polygonization (an earlier lookup in this job pulled all 58 official HF benchmarks [docs/implementation2.0/notes/H1_geometry_algorithms.md]). The CCL benchmark is YACCLAB [Bolelli 2024 review], cited instead. No model choice is involved. |
| `huggingface-skills:huggingface-papers` | `api/papers/2610.37297` (TRUCO) returns no record (not indexed on HF), and `api/papers/search?q=connected component labeling GPU` returns no CCL paper. So the papers come from arXiv, the journals and the authors' PDFs (section 11). |

## 1. Problem

`objdel.checker.polygonize` turns a binary mask into exact rectilinear rings: crack edges, (8,4) connectivity, outer
rings clockwise with area > 0, holes with area < 0, corners from a canonical start, rings in a fixed order. It is numpy
with an `argsort` link and pointer jumping [code: REPO/src/objdel/checker/polygonize.py:72-115]. Production calls it
on the CPU [code: REPO/src/objdel/targets/build.py:93, REPO/src/objdel/checker/certify.py:46]. The CUDA port
`check_masks` exists [code: REPO/cuda/gpucheck.cpp:751-1057], but it works on the dense slot grid S = 4(W+1)(H+1). Each
of its 2·`jump_rounds(S)` rounds touches every slot, and it runs one thread per vertex row. For a 512² frame that is
1.05 M slots × 2 × 22 rounds ≈ 46 M slot visits per mask [code: gpucheck.cpp:931-937, `jump_rounds`]. A lesion's ring
has a few hundred edges.

Goal: the same output, field for field (`quads`, `ring_start`, `pts_start`, `ring_area`, `pts`), with work
proportional to the input words plus the runs plus the boundary edges. One source for CPU and GPU, in gpucheck.cpp's
style.

## 2. Algorithm

```
bytes ──Pack──► bit words ──RunCount──► runs/row ──scan──► [sync 1: run totals]
                   │                                            │
                   └──Quads (bit slicing) ──► n1 n2 nD n3 n4    ▼
                                                     RunEmit (fg runs, bg runs, sentinel OUTSIDE per mask)
                                                            │
                                     Merge: union-find, fg 8-adjacent, bg 4-adjacent, border bg ∪ OUTSIDE
                                                            │
                                     Flatten ──► root = first run in raster order (min index)
                                                            │
       HoleKey: per hole B, min raster index of a fg pixel of S(B) whose upper neighbour is in B (atomicMin)
                                                            │
       RingMark: order bitmap (1 bit per ring at its order pixel), start bitmap (1 bit at its start pixel)
                                                            │
       Popc + scan ──► ring index = rank in the order bitmap (polygonize's ring order, no sort)
                                                            │
       RingInit: start edge in closed form (outer: first run's (x0,y,east); hole: (x0,y,south))
       Seeds: the R ring starts + (K > 0) both vertical sides of each fg run on rows y % K == 0
                                                            │
       Trace 1: each seed walks boundary edges (2x2 rule, left turn at saddles) to the next seed
                → corners, area, next seed
       RingWalk: per ring, follow the seeds; each segment gets the offset of its first corner
       scan of ring corner counts ──► output offsets
       Trace 2: each seed walks its segment again and writes its corners in place
                                                            │
                                         [sync 2: totals] [sync 3: outputs]
```

Code: `k1poly.cpp` (`Pack`/`pack_kernel` 357-400, `Quads` 402, `RunCount` 429, `RunEmit` 448, `Merge` 492, `Flatten`
523, `HoleKey` 535, `RingMark` 563, `RingInit` 604, `SeedFlag`/`SeedEmit` 631-657, `Trace` 669, `RingWalk` 727,
`chunk` 763). The union is the lock-free min-root union of Playne and Hawick (`unite`, 339). Successor and
predecessor directions come from the bit-quad (`Geo::succ_dir`, `Geo::pred_dir`, 300-321). The rule is the same as
gpucheck.cpp `InitSlots` and polygonize `link`.

## 3. Decision record (ADR-K1)

Context: an exact polygonizer for batches of masks up to 512² (any size up to 4096²). It must be identical to the
numpy reference, and the same code must run on the CPU for the Mac tests.

| Option | Work per mask | Span | Exactness | Complexity |
|---|---|---|---|---|
| A. dense slot grid + pointer jumping (now: check_masks) | O(S log S), S = 4(W+1)(H+1) | O(log S) rounds of O(W) | proven, measured | low |
| B. compacted edge list + Wyllie list ranking (H1 §4c) | O(N/32 + E log E) | O(log E) | same keys | medium |
| C. CCL-derived ring identity + seeded crack tracing (K1) | O(N/32 + F + B + E) | O(K·W + seeds per ring) | closed-form keys, proofs §4 | medium |
| D. block-based pixel CCL (BKE) + tracing | O(N) labels + E | low | same | high (2x2 block logic, bg 4-connectivity separately) |

Decision: C. It is the only work-optimal one. It also gives the component labels (fg 8, bg 4, hole → part) that
`targets/build.py` now gets from OpenCV and that the scorer's signature labeling needs. Run-based CCL fits sparse
lesion masks: few runs per row. Bolelli et al.'s review finds block-based fastest on general images and run-based
weakest for vertically elongated shapes [Bolelli 2024]. This is a measured risk (§10), not a blocker. Consequences:
three host syncs per chunk, run-sized scratch (bounded by halving), and one sequential walk per ring over its seeds.

## 4. Exactness argument

Notation: each boundary edge has a foreground pixel f on its right and a background pixel b on its left. C(f) is the
8-component of f, and B(b) is the 4-component of b. The border of the image is background (polygonize pads with 0), so
every background pixel on the image border is in one component, OUTSIDE.

Lemma 1 (one ring per adjacent pair). In the (8,4) pairing, the border between a foreground component and an
adjacent background component is one closed crack curve [Rosenfeld 1974; Kong and Rosenfeld 1989]. Along a ring, the
right pixels of consecutive edges are 4- or 8-adjacent (the left turn at saddles joins diagonal foreground). The left
pixels are 4-adjacent. So all edges of a ring share one pair (C, B), and each pair is one ring. Each C has exactly one
adjacent component that surrounds it. Rings are therefore {outer ring of each C} ∪ {hole ring of each bounded B}.
The test checks the ring counts on every mask.

Lemma 2 (outer start). Let p0 = (x0, y0) be the first raster pixel of C. Its upper neighbour is background (else it
would be in C and come earlier). That background pixel cannot be in a hole of C: a hole of C needs C pixels above it.
So the top edge of p0 is on C's outer ring. Now consider the edges leaving vertices in row y0. East edges need a C
pixel (x, y0), so x ≥ x0. South edges need a C pixel (x-1, y0), so x ≥ x0 + 1. West and north edges need C pixels in
row y0 - 1, and there are none. No vertex of the ring lies above row y0. So the minimum key (y(W+1)+x)·4+d is
(x0, y0, east), which is polygonize's canonical start [code: polygonize.py:130-133].

Lemma 3 (hole start and surround). Let q0 = (bx, by) be the first raster pixel of a bounded B. Its left neighbour and
its upper neighbour are not in B (they come earlier). They are not background either: a background pixel there would
be 4-adjacent and so in B. B is bounded, so bx ≥ 1 and by ≥ 1. Both neighbours are foreground and diagonal to each
other, so they are in one component, S(B). The edge on the left side of q0 goes south from vertex (bx, by), with q0 on
its left. Edges with B on the left that leave row by: south edges need B pixel (x, by), so x ≥ bx; west edges need
x ≥ bx + 1; east and north edges need B pixels in row by - 1. So (bx, by, south) is the minimum key of B's ring,
provided that edge is on B's ring and not on an island's ring. Take F = B plus everything B encloses. Anything B
encloses lies in rows ≥ by and, in row by, at x > bx: a vertical path upward from it would otherwise avoid B. So q0
is also F's first raster pixel, and the edge left of q0 is on F's outer border. That border is B's hole ring. At a
saddle (code 6) the north edge at the same vertex has key +3 and belongs to the upper-left pixel's ring. Then S(B) is
the component of the foreground run just left of B's first run (`RunEmit` stores it as `bleft`).

Lemma 4 (ring order). polygonize labels each ring by its minimum edge index and orders the rings by that label
[code: polygonize.py:125-127]. `boundary_edges` lists the up-side (east) edges first, in raster order [code:
polygonize.py:60-65]. Every closed ring has an east edge (its horizontal steps sum to zero). So a ring's order key is
the raster index of its first east edge, as gpucheck's `rkey` also uses [code: gpucheck.cpp:908]. For an outer ring
the key is p0 (lemma 2: every east edge of the ring is the top of a C pixel, and p0 comes first). For a hole the key
is the minimum raster index of a foreground pixel of S(B) whose upper neighbour is in B (lemma 1: an edge is on B's
ring iff its pair is (S(B), B)). `HoleKey` computes this per overlapping (fg run, bg run above) pair with
`atomicMin`. Keys are distinct, so ring index = popcount rank in the order bitmap gives polygonize's order with no
sort.

Lemma 5 (segments). The successor map is a permutation of the boundary edges. Seeds are recognised locally: the
start-bitmap bit at the edge's pixel for ring starts, and the seed map for vertical sides on seed rows. A row seed
that is also a ring start (a hole start on a seed row) is marked dead, so each edge has one seed id. Each segment runs
from a seed to the next seed, so the segments partition every ring. The walk from the ring's start seed visits them in
ring order. A segment cannot cross a seed row without a vertical edge on that row, and such an edge is a seed. So a
segment stays inside one band of K-1 pixel rows.

Lemma 6 (corners). An edge is a corner iff its direction differs from its predecessor's (polygonize `corner = d !=
d[pred]`). The ring start is always a corner: it arrives with d ∈ {2,3} and leaves with d ∈ {0,1}. A segment's first
edge gets its predecessor's direction from `pred_dir`. Corner k of a ring is written at ring base + segment offset +
k-th corner in the segment, so the output order is the ring order from the start.

Order independence (CPU = GPU). Every value is an integer. Cross-thread writes are (a) `atomicMin` / `atomicOr` /
`atomicAdd` on integers, which commute; (b) union-find links, whose final partition is the connected components
whatever the order, with root = minimum index; or (c) writes of an item's own slots. Every read of another item's data
happens in a later kernel. The CPU build checks this by running all items forward, reversed and permuted
(`K1_ORDER` 0/1/2). The GPU run checks it under the real schedule.

## 5. Complexity

N = pixels, F/B = foreground/background runs, E = boundary edges, R = rings, Σ = seeds = R + 2·F_K (F_K = fg runs on
seed rows ≈ F/K), K = seed rows (default 16).

| Step | Work | Span (per item) |
|---|---|---|
| Pack | O(N) byte reads (the input is bytes; O(N/32) if given words) | 1 warp-ballot per word |
| Quads, RunCount, bitmap popcounts | O(N/32) | O(W/32) per row |
| RunEmit | O(N/32 + F + B) | O(W/32 + runs in row) |
| Merge (+ Flatten) | O(F + B) unions × find depth | find depth (no proven bound, §10) |
| HoleKey, RingMark, RingInit, seeds | O(F + B) | O(runs in row) |
| scans (CUB) | O(rows + N/32 + F + R) | O(log) |
| Trace 1 + Trace 2 | 2E + O(Σ) | O(K·W) edges per segment (lemma 5) |
| RingWalk | O(Σ) | seeds of the ring: O(H·W/K) worst case |

Total work: O(N/32 + F + B + E), plus the unavoidable byte read. check_masks does O(S log S). For a 512² lesion frame
the ratio is about 46 M slot visits against roughly 8 K words + 100 runs + 1 K edges [estimate]. The CPU measurement
below shows ×95.

## 6. Memory and occupancy (RTX 5090, CC 12.0)

- Layout A (inputs, words, two bitmaps, their counts and scans, the seed map, rows): about 2.1 B per pixel. That is
  about 70 MB for the default chunk of 2^25 pixels (128 frames of 512²) [code: `chunk`, layout A].
- Layout B (runs, rings, seeds, corners) is sized by the run totals after sync 1. Fg run 56 B, bg run 24 B, ring 20 B,
  seed 44 B. That is tiny for lesions and about 150 B per pixel for a checkerboard (a run per pixel pair, a ring per
  pixel). Above `max_scratch_bytes` (1 GiB) the chunk is halved and redone (`run_range`), down to one mask [code].
  check_masks keeps up to about 0.5 GB of scratch per chunk [code: gpucheck.cpp `kChunkSlots`].
- Occupancy: every K1 kernel uses ≤ 40 registers and no shared memory [measured: nvcc 13.0.88 `-Xptxas -v`, sm_120,
  results/nvcc_docker_sm120.log: Trace 40, RunEmit 40, Quads 37, RunCount 36, HoleKey 36, the rest ≤ 30]. CUB's scan uses 71 registers and 33,808 B of shared memory. CC 12.0 allows 48 warps,
  32 blocks and 64 K registers per SM [doc: Blackwell Tuning Guide 1.4.1.1]. 128-thread blocks at ≤ 40 registers
  can therefore keep 48 warps resident (64 K / 1536 = 42.7). Trace is latency-bound (a dependent bit-word load per
  step). A 512² mask's words are 32 KB, so a whole 128-frame chunk (4 MB of words) fits in L2 (`device_info()` prints
  the L2 size on the box).
- Host: pinned staging buffers that only grow. One stream (non-blocking). Three syncs per chunk: run totals, counts,
  outputs.

## 7. Verification done (Mac, Apple M3, clang 21, CPU build)

- `test_k1_cpu 3000 41` (and `_o1`, `_o2`, the reversed and permuted item orders): PASS, 4,453 masks, 453,142 rings,
  0 failures each [measured: run in this job]. Each mask runs 15 K1 configurations (seed rows 0/1/2/3/16 × one chunk
  / 700-pixel chunks / 1-byte scratch limit). Each is compared field by field with `ref_polygonize` (a literal port of
  polygonize.py + bitquad) and `gpucheck::check_masks` (CPU build). The masks cover: noise at densities 0.05-0.95,
  blobs with holes, checkerboard (all saddles), nested squares, comb (one ring of about 262 K edges), spiral,
  diagonals, stripes, a frame with 1-pixel holes, 1×N, N×1, sizes around 32 and 33, and 0×W, H×0 and 0×0 masks,
  inside mixed batches.
- `cross_numpy.py` against `objdel.checker.polygonize` + `bitquads` (repo venv, `-I`): 3,000 random masks with K=16
  (311,518 rings), 1,500 with K=0 (153,724 rings), and all 37,093 LIDC crops of data/checker_eval/masks.npz (38,002
  rings, 1,476,508 corners): 0 mismatches. `check_masks` is IDENTICAL in-process [measured].
- 8 mutants: 7 caught in forward order. The 8th (first write instead of min) is caught in reversed and permuted order.
- ASan + UBSan (permuted order, 400 cases): clean.
- nvcc 13.0.88 `-arch=sm_120`: k1poly, gpucheck (as CUDA) and the three GPU binaries compile and link in a linux/amd64
  container (`nvcc_docker_check.sh`). The only warning: `Pack::operator()` is unused in the CUDA build.

## 8. Measured speed (CPU build, single thread; the GPU numbers come from GPU_RUN.md step 5)

`bench_k1_cpu 3`, median of 3, K1 against check_masks, both CPU builds, every result IDENTICAL [measured: Apple M3,
results/bench_cpu_m3.txt]:

| Workload | K1 ms | check_masks ms | speedup |
|---|---|---|---|
| lesion512 b256 (K=16) | 41.1 | 3,927 | ×95 |
| noise512 b16 (286 K rings) | 194 | 1,241 | ×6.4 |
| comb512 b4, K=16 / K=0 | 32.4 / 28.0 | 115 / 113 | ×3.6 / ×4.0 |
| spiral512 b4 | 23.6 | 168 | ×7.1 |
| concentric512 b4 | 11.7 | 130 | ×11 |
| checker512 b4 (520 K rings) | 79.3 | 116 | ×1.5 |
| real LIDC crops, all 37,093 | 218 | 4,834 | ×22 (170 K vs 7.7 K masks/s) |

On the CPU, seed rows only cost time (`pred_dir` per seed); they pay off only with parallel threads. On the GPU K1's
end-to-end time for 512² frames should be dominated by moving bytes: 256 frames are 67 MB of mask bytes, staged into
pinned memory and copied over PCIe [estimate]. The device-side stages should be far smaller. Step 5 measures both
(`h2d` against the rest).

## 9. Integration points (REPO)

1. `cuda/gpucheck.cpp:check_mask_chunk` (966-1057) and its functors `InitSlots`, `LabelRound`, `Pred`, `RankInit`,
   `RankRound`, `RingStats`, `ScatterCorners`, `GatherRings`, `SetBase` (806-964). Replace the body with K1's stages:
   same `MaskResult`, same binding (`cuda/gpucheck_py.cpp:140-159`). Archive the slot-grid functors (owner rule).
   Reuse gpucheck's `Lane` (lowest-priority stream, so the trainer's kernels go first) instead of K1's own `Engine`.
2. `src/objdel/targets/build.py:93` (`polygonize(sub, offset=(x0, y0))` per object, plus `cv2.connectedComponents`
   at :91 to give each hole its part). Make one batched `_gpucheck.check_masks` call over all crops of a build and add
   the offsets on the host. Optional extension: return `ring_part` (outer: its fg root; hole: S(B)'s root). K1 already
   has it after `Flatten`, so OpenCV can go.
3. `src/objdel/checker/certify.py:46`: the three-way Euler certification keeps numpy polygonize as one independent
   leg. Keep it as the oracle and add K1 as the device leg.
4. Scorer signature labeling (`cuda/gpucheck.cpp` `LabelInit`/`LabelMerge`/`LabelFlatten`/`Holes`, 1302-1407). This
   is pixel union-find on the crop. K1's run CCL (fg 8 + bg 4 with the OUTSIDE sentinel, hole → S(B)) gives the same
   partition from runs. It is a later candidate (separate exactness gate: the signature tests).
5. Tests: port `test_k1.cpp`'s generators and the K1_ORDER idea into `tests/test_cuda_checker.py` (mask fuzz,
   131-158). `experiments/cuda_checker/compare.py:76` and `eval.py` stay the corpus gate (195,809 cases on the
   RTX 4090 [docs/implementation2.0/notes/B_gpu_audit.md]).

## 10. Risks and what to revisit

1. Union-find depth on the GPU. There is no path compression during `Merge` (path compression there races with
   `atomicMin`), so under arbitrary scheduling chains can reach O(H) for vertical snakes. The CPU runs are depth 1
   because rows go in order. The spiral and comb lines in step 5 measure it. Fix if slow: a per-row local merge
   first, or Komura-style equivalence, or BKE blocks [Allegretti 2019].
2. Span of very long rings. RingWalk is sequential over one ring's seeds (comb512: 32 seed rows × 256 runs × 2 =
   16,384 seeds at K=16). Fix: pointer jumping or a sparse ruling set on the seed list [Helman-JáJá 2001; Rehman 2009]. That is
   O(Σ log Σ) on a list K times shorter than E.
3. Byte input. For 512² frames, staging and H2D of mask bytes likely dominate. Better integration: accept bit-packed
   rows (gpucheck's raster output is already bit-packed), or keep masks on the device. This needs a
   `polygonize_words` entry; `Pack` is the only step that reads bytes.
4. Three syncs per chunk. Sync 1 (run totals) could go with a persistent upper-bounded run budget.
5. Not yet run on a GPU. `pack_kernel` (ballot), the CUB scans and the real concurrency of atomics/union-find are
   covered only by GPU_RUN.md steps 3-4. The CPU K1_ORDER runs cover order independence, not memory-model issues;
   loads of union-find parents are `volatile`, as in gpucheck.
6. Index width: i32 inside a chunk (2^25 pixels by default; 4·NF corners ≤ 2^31 is checked). A 4096² mask runs
   alone, which is fine.
7. Ring order convention. K1 keeps polygonize's "first east edge" order. H1 F11 notes that ordering holes by their
   canonical start instead would match `shape.canonicalize`. That change is easy here (rank the start bitmap
   instead), but it changes target hashes, so it is the owner's call.
8. Divergence in Trace. One thread per seed with segments of different lengths. This is fine for lesion rings of
   hundreds of edges. If profiles show idle warps, sort the seeds by expected length or use a warp per ring.

## 11. Citations

- S. Suzuki, K. Abe, "Topological structural analysis of digitized binary images by border following", CVGIP 30(1):
  32-46, 1985. doi:10.1016/0734-189X(85)90016-7
- V. A. Kovalevsky, "Finite topology as applied to image analysis", CVGIP 46(2):141-161, 1989 (cracks).
- A. Rosenfeld, "Adjacency in digital pictures", Information and Control 26:24-33, 1974; T. Y. Kong, A. Rosenfeld,
  "Digital topology: introduction and survey", CVGIP 48(3):357-393, 1989 ((8,4) duality, one border per adjacent pair).
- D. P. Playne, K. Hawick, "A new algorithm for parallel connected-component labelling on GPUs", IEEE TPDS 29(6):
  1217-1230, 2018 (lock-free min-root union, `unite`).
- S. Allegretti, F. Bolelli, C. Grana, "Optimized block-based algorithms to label connected components on GPUs",
  IEEE TPDS 31(2):423-438, 2020 (BUF, BKE).
- F. Bolelli et al., "Connected Components Labeling on GPUs: a state-of-the-art review with code", IEEE TPDS, 2024
  (https://federicobolelli.it/media/publications/pdfs/2024tpds.pdf): block-based fastest overall, run-based weakest on
  vertically elongated shapes; HA4/HA8 run-based with `__ballot_sync`, `__ffs`, `__clz`.
- A. Hennequin, L. Lacassagne, L. Cabaret, Q. Meunier, "A new direct connected component labeling and analysis
  algorithms for GPUs", DASIP 2018 (HA4).
- F. Lemaitre, N. Maurice, L. Lacassagne, "An efficient run-based connected component labeling algorithm for
  processing holes", 2022 (https://largo.lip6.fr/~lacas/Publications/BNBW22_FLSL+holes.pdf): fg and bg run labels
  and the adjacency tree in one pass, the closest prior art to K1's fg/bg runs + surround.
- V. M. Garcia-Molla, P. Alonso-Jordá, R. García-Laguía, "Parallel border tracking in binary images using GPUs",
  J. Supercomputing 2022, doi:10.1007/s11227-021-04260-y (tiles, one thread per tile, stitching).
- P. Alonso-Jordá, R. Díaz-Cano, E. S. Quintana-Ortí, F. Folch, "Accelerated border tracking in binary images with
  GPUs", J. Supercomputing 2026, doi:10.1007/s11227-026-08488-4 (triads linked inside and across blocks).
- R. Muñoz-Salinas, F. J. Romero-Ramírez, M. J. Marín-Jiménez, "TRUCO: a scalable lock-free algorithm for parallel
  contour extraction", arXiv:2610.37297, 2026: exact Suzuki-Abe output, a contour owned by the stripe of its top row
  (K1's closed-form starts play the same role). In OpenCV `findContours` (PR 28773).
- Yin et al., "G-Contour: GPU accelerated contour tracing for VLSI layouts", ICCAD 2025
  (https://www.cse.cuhk.edu.hk/~byu/papers/C295-ICCAD2025-GContour.pdf): union-find → start points → parallel
  per-contour search, the same macro-structure as K1.
- D. R. Helman, J. JáJá, "Prefix computations on symmetric multiprocessors", JPDC 61(2), 2001 (sparse ruling set); M.
  S. Rehman, K. Kothapalli, P. J. Narayanan, "Fast and scalable list ranking on the GPU", ICS 2009; J. C. Wyllie,
  Cornell TR79-387, 1979.
- NVIDIA: Blackwell Tuning Guide §1.4.1.1 (CC 12.0 occupancy limits), CUDA Programming Guide "Compute Capabilities"
  (CC 12.x: 1,536 threads and 100 KB shared memory per SM), CUB `DeviceScan::ExclusiveSum`, warp vote functions
  (`__ballot_sync`).
- Repo sources: REPO/src/objdel/checker/polygonize.py, REPO/src/objdel/checker/bitquad.py, REPO/cuda/gpucheck.cpp,
  docs/implementation2.0/notes/H1_geometry_algorithms.md §4 (the closed-form-start idea that K1 implements and proves), docs/implementation2.0/notes/B_gpu_audit.md §3.2.

## Files

`k1poly.h`, `k1poly.cpp` (single source), `masks_gen.h` (generators), `test_k1.cpp` (exactness), `k1_dump.cpp` +
`cross_numpy.py` (numpy cross-check), `bench_k1.cpp` (speed), `build_cpu.sh`, `build_gpu.sh`, `nvcc_docker_check.sh`,
`GPU_RUN.md`, `results/` (logs).
