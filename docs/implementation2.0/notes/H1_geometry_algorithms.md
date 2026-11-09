# H1: Provably good algorithms for the geometry core, run exactly and in parallel on the GPU

Date: 2026-11-08. Scope: read-only review of `cuda/gpucheck.cpp` and `src/objdel/{checker,geometry,targets}`, plus
the literature. I ran no code on a GPU. Every speed figure below is an **estimate**, with its arithmetic shown;
none is a measurement.

## Skills used

| Skill | What it contributed |
|---|---|
| `nvidia-skills:nvidia-skill-finder` | Live catalog check (`raw.githubusercontent.com/NVIDIA/skills/main/skills.sh.json`, "GPU Development" group). The only relevant entry is **`compute-sanitizer-diagnose-cuda-bugs`**. Use it (racecheck, synccheck, memcheck) on the union-find and pointer-jumping kernels whenever they change (section 9). No catalog skill covers CCL, EDT or list ranking. I installed nothing. |
| `engineering:code-review` | Review lens applied to the code: correctness (edge cases, races, overflow), performance (O(n^2) in hot paths, wasted work), maintainability. Gives the "Findings" table (section 1). |
| `engineering:debug` | Reproduce, isolate, root-cause applied to the worst-case latency hazards: a single thread doing `ShapeTotals`, dense-slot rounds, row-serial labeling. Each finding names its root cause and the guard to add. |
| `engineering:system-design` | Requirements first (exact, bit-reproducible, CUDA for checker, scorer and polygonizer, 2D now and 3D later). Then data flow (one GPU-resident pipeline from bit words to per-object numbers), then trade-offs. Section 8 holds the pipeline diagram. |
| `engineering:architecture` | ADR format for the two big decisions: ring identity from CCL labels instead of pointer jumping, and run-based CCL on bit words (section 8). |
| `engineering:documentation` | Layout: the most useful information first (findings, then one section per primitive with fixed parts a-d), links rather than copied text. |
| `huggingface-skills:huggingface-best` | The main thread had already pulled all 58 official HF benchmarks (`MAIN_skill_lookups.md`). None measures geometry kernels (CCL, EDT, polygonization). The kernel benchmark that does exist is YACCLAB (Allegretti et al.), which is cited instead. |
| `huggingface-skills:huggingface-papers` | Paper lookup by arXiv ID: CubicalRipser 2005.12692 and the duality of image persistence 2005.04597. The HF API returned no metadata for either, so I confirmed both on arXiv via web search. |
| `mattpocock-skills:codebase-design` | The "deep module" view of the proposed design: one geometry core with a small interface (rasterize, label, trace, topology, distance). The CPU build is a second adapter at the same seam, so the seam is real (section 8). |

---

## 0. Executive summary

The repo's core is already well founded. Its predicates are exact integer arithmetic, its labels are min-index
(deterministic) and its distance sums are exact integers. Its maths has no correctness bugs that I could find.
Its algorithms, though, are **not work-optimal on the GPU** in four places, and one production path still runs on
the CPU:

1. **Mask polygonizer (`check_masks`)** works on the **dense slot grid** S = 4(W+1)(H+1) (gpucheck.cpp:795-804,
   937), not on the boundary edges E. It runs **one thread per vertex row**, with ~2·log2 S kernel rounds. For a
   256x256 mask that is 264,196 slots x 40 rounds, where ~1,000 edges x ~20 rounds would do: about 500x more work
   than needed. Worse, the production target builder does not use it at all. It calls the CPU numpy
   `polygonize` (`targets/build.py:93`), which breaks the owner's "CUDA for the polygonizer" rule.
2. **Ring identity needs no pointer jumping.** In the (8,4) topology every ring is the border between one
   8-connected foreground component and one 4-connected background component. Its canonical start has a
   **closed form**: the top-left corner of that component's first raster pixel, leaving east for an outer ring
   and south for a hole. The labels that give this are already computed by the scorer's CCL. So the
   O(E log E) min-label pass becomes O(E), with no log-depth rounds. Only the corner order still needs list
   ranking.
3. **CCL in the scorer** works row-serially: one thread per crop row, byte pixels, so its loads do not coalesce
   (gpucheck.cpp:1314-1334). It always runs **4 labelings**, where the Python path usually runs 2
   (`topology.signature`:161 skips the filter when nothing is small). Holes per component can come from a
   **per-component bit-quad Euler number** (each 2x2 window's foreground pixels are mutually 8-adjacent, so
   each window belongs to exactly one component), which removes the last background labeling.
4. **`ShapeTotals` runs one GPU thread per shape**, with O(R^2·n^2) sequential work when rings touch
   (gpucheck.cpp:626-749). That is up to ~9·10^7 dependent operations on one thread for 10 rings of 1,000
   vertices. On adversarial completions this thread straggles and stalls the whole batch.
   Pair classification is O(E^2) per shape (gpucheck.cpp:585-623). A uniform grid with reference-point
   de-duplication keeps the polycore counts exact at O(E + candidate pairs).

Everything below is exact by construction. Every output is an integer, or a correctly rounded function of
integers (sqrt), and is uniquely defined by the input. **Any** correct parallel algorithm therefore returns
bit-identical results, whatever the schedule or atomic order. This is the doctrine that lets the design pick
fast, nondeterministically scheduled algorithms freely (section 7).

---

## 1. Findings (bugs and suboptimal choices)

| # | Where | Finding | Severity |
|---|---|---|---|
| F1 | `targets/build.py:80-97`, `checker/polygonize.py:72-89` | The production polygonizer is CPU numpy: `argsort` linking costs O(E log E), plus OpenCV CCL. The CUDA `check_masks` is used only by tests and experiments (`tests/test_cuda_checker.py:152`, `experiments/cuda_checker/compare.py:76`). This breaks the owner rule. | High (policy) |
| F2 | gpucheck.cpp:795-804 (`for_row_slots`), 845-896, 931-937 | The polygonizer touches all S = 4(W+1)(H+1) slots in each of 2·`jump_rounds(max_slots)` rounds, sized by slots rather than by ring length. It runs one thread per vertex row, so a 256² mask gets only 257 threads. Work is O(S log S); the optimum is O(E). | High (perf) |
| F3 | gpucheck.cpp:1019-1056 | Three host round trips per chunk (ring count; ring info, sorted on the host; corners). Closed-form starts plus a device scan cut this to one. | Medium |
| F4 | gpucheck.cpp:1314-1334, 1602-1617 | CCL is row-per-thread on a byte copy of the mask (`CropFill` makes the data 8x larger, 1284-1299), and the loads do not coalesce. It runs 4 labelings every time, even when the size filter changes nothing. | Medium-High (perf) |
| F5 | gpucheck.cpp:626-749 (`locate2`, `any_witness`, `locate2_region`) | One thread per shape does O(R^2·n·n) work when rings are not "apart". The worst case (R=10, n=1,000) is ~10^8 serial operations on one thread, a ~0.1 s straggler (estimate: ~1 op per cycle at ~2 GHz, before memory stalls). | Medium (tail latency) |
| F6 | gpucheck.cpp:585-623 | Edge-pair classification is O(E^2) per shape, with triangular load imbalance (edge k scans E-k partners). Typical E of 100-300 makes this negligible. Adversarial E = 10^4 gives 5·10^7 pair tests per shape. | Low-Medium |
| F7 | gpucheck.cpp:537-541, 309-321 | `orient<i64>` is used whenever max|coord| < 2^28. On the 0..1000 grid every orientation fits in **int32** (\|diff\| ≤ 2,000, products ≤ 4·10^6). A 64-bit multiply is several 32-bit instructions on NVIDIA GPUs. Add an i32 tier for max|coord| < 2^15. | Low (perf) |
| F8 | gpucheck.cpp:806-842 | `InitSlots` computes the 2x2 code with up to 20 bounds-checked byte loads per vertex: `quad()` at the vertex plus up to 4 more at end vertices. Bit-sliced codes from two packed rows give 32 vertices per ~20 ALU operations. | Low (perf) |
| F9 | `targets/simplify.py:158-198` | `_best_cycle` costs O(c·n·L) (c = shortcuts over the anchor), which is O(n^3) in the worst case (long, nearly straight rings). It is offline, so it is not urgent. | Low |
| F10 | `checker/volume.py:51-59`, `eval/metrics.py:199-206`, `reward/cache.py:130-135` | 3D Betti numbers (scipy `ndimage.label`), eval NSD EDT and cache EDT are all CPU. The cache EDT is offline, once per item, which is acceptable. 3D topology and eval EDT would need CUDA under the owner rule. | Medium (policy, for 3D) |
| F11 | polygonize.py:127, gpucheck.cpp:908, 1032-1034 | Not a bug, a convention: rings are ordered by their **first east edge** in raster order. For a hole this edge is on the hole's lower boundary, not at its canonical start. `shape.canonicalize` orders by the (y, x) start sequence. Ordering by the closed-form start (section 4) would make the two agree. It would change no ring and no vertex order, only the ring order, so target hashes would need a re-check. | Info |

Checked and found **correct**: the left turn at saddles joins diagonal foreground, so (8,4) holds (verified by hand on
code 6). The corner rank `RankInit`/`RankRound` gives index(e) = #corners in [start, pred(e)]. The `Holes`
attribution argument holds. The raster window suffix sums and the half-open crossing rule are right, and so are
the i64/i128 precision bounds. The union-find retry-on-failed-`atomicMin` is the Playne-Hawick invariant. The
`rmax = S/8 + 1` and corner ≤ S/2 bounds hold. The numpy shoelace via float `bincount` is exact (|term| ≤ 4,096
per edge, so the sum stays far below 2^53). The 3D Euler formula V−E+F−C over closed cubes, with b2 = 6-connected
bounded background, is right.

---

## 2. Rasterization of polygons with holes (even-odd at pixel centres)

**a) Now.** `raster_row_warp` (gpucheck.cpp:388-415). One warp per pixel row. Each lane tests every edge of the
shape (`window_crossings`, 357-372) with the exact `RowLine::cross` (338-351), in i64/i128 scaled coordinates.
It scatters ±1 into a 1,024-column shared-memory window, takes warp suffix sums with `__shfl_down_sync`, and
writes bits with `__ballot_sync`. Work O(H·E + H·W/32) per shape. Depth O(E/32 + (W/32)·log 32). Exact:
integer crossing columns, an integer winding sum, parity read from it. Atomics are integer, so the result is
deterministic.

**b) Best known.** Scanline fill with an active edge table: O(E log E + K + HW/w) sequential, where K = Σ rows
spanned by the edges ("crossings"). On PRAM/GPU an **edge-centric** fill emits each edge's K_e crossings (counts,
then a prefix sum), scatters ±1 into per-row difference arrays, and runs a segmented scan (or XOR scan for
even-odd) per row. That is O(E + K + HW/w) work and O(log(E + K)) depth, the optimum, since every output word
and every crossing must be touched. The literature is the classic polygon scan conversion (Foley-van Dam). The
GPU form is any "sort-free coverage" rasterizer. Nothing beyond integer arithmetic is needed for exactness.

**c) Recommendation.** Keep the current kernel. With E ≤ ~300 typical and H ≤ 512, H·E ≈ 1.5·10^5 tests per shape,
which is negligible. Add a guard for adversarial inputs: when E·H > ~10^7, bucket edges by row range (sort edges
by y_min with one block-level radix sort per shape, then each warp scans only its active edges). That restores
O(E + K). Exactness is unchanged: the same `cross` predicate and integer sums. Add the i32 tier when
max|coord|·2·dim < 2^31 (it never holds for W up to 4,096 at G = 1,000, so i64 stays; this note is for the
validation side, F7).

**d) 3D.** Voxelizing a closed triangle mesh at voxel centres is the same algorithm one dimension up. For each
(x, y) column, count the signed ray crossings with triangles using exact orientation (integer 3x3 determinants;
i128 when coordinates are large), with a half-open rule for ties. For this project the 3D masks are stacks of
2D slices (`volume.stack`), so the 2D rasterizer suffices.

---

## 3. Ring validity / simple-polygon test

**a) Now.** `CleanRings` (479-520): one thread per ring, which drops repeats and computes the box and the i128
shoelace. `EdgePairs` (585-623): one thread per edge, testing it against every later edge of the shape with a
box filter and then exact `classify` (550-582; orientation, collinear overlap, spike, pinch). This gives counts
per defect kind and per ring pair. `ShapeTotals` (626-749): one thread per shape, for winding, hole-in-outer and
nesting through exact witness points (vertices and doubled-coordinate edge midpoints) with ray-crossing
`locate2`. Work O(E^2) for the pairs, O(R^2·n^2) worst for the witnesses. Depth O(E), and O(R^2 n^2) on the
`ShapeTotals` thread. Exact: integer orientation in i64 (or i128).

**b) Best known.**
- Detection (is any pair improper?): **Shamos-Hoey** 1976, O(n log n), optimal in the algebraic decision tree
  model. A simple-polygon test is in theory O(n) via Chazelle's linear triangulation (1991), which nobody
  implements.
- Reporting all k pairs: Bentley-Ottmann 1979, O((n+k) log n). Optimal is O(n log n + k): Chazelle-Edelsbrunner
  1992 and Balaban 1995 (O(n) space).
- PRAM: Goodrich, "Intersecting line segments in parallel with an output-sensitive number of processors",
  SIAM J. Comput. 20(4), 1991 (SPAA'89): O(log n) time, O(n + k/log n) processors on a CREW PRAM, builds the
  arrangement. https://epubs.siam.org/doi/10.1137/0220047
- Practice on GPU: **uniform grid binning** (W. R. Franklin's line of work). Bin edges into cells of side
  ~G/√E, test pairs inside each cell, and count a pair only in the cell holding the reference point
  (min x, min y) of box(e) ∩ box(f), so no pair is counted twice. Expected O(E + k) work for "reasonable"
  inputs, and exact with integer predicates.
- Exact predicates: integer orientation is exact with no filter. Shewchuk's adaptive predicates are needed only
  for floating input, which this project does not have.

**c) Recommendation.** The checker must reproduce polycore's **counts**, not just a yes/no, so reporting, not
detection, is the problem to solve.
1. Per shape, one thread block. Build a uniform grid over the shape's box (cell side ≈ box/√E). Each edge
   increments the counts of the cells it covers; a block scan gives the offsets; then scatter. Each candidate
   pair (e, f) in a cell is classified only when the cell contains the reference point. The integer counts are
   then identical to the O(E^2) loop. Work O(E + C), with C = candidate pairs.
2. Add an int32 `orient` tier (F7).
3. Parallelize `ShapeTotals`. One warp per (ring a, ring b) pair; lanes split the witnesses of a. Each `locate2`
   becomes a warp-cooperative crossing parity: each lane tests n/32 edges, `parity = popc(ballot(cross)) & 1`
   accumulated over chunks, and on-boundary = `any(ballot(on))`. That turns O(n^2) serial work into
   O(n^2/32) per warp across R^2 warps. Result: identical integers.
4. An optional early exit for the RL path only: if any crossing exists, the object is illegal and its terms are
   zeroed by the legality gate anyway. Keep full counts for logging and tests.

Expected speed (estimate): typical shapes are already < 10 µs. The fixes target the adversarial tail. A
1,000-vertex x 10-ring shape goes from ~0.1 s on one thread (F5) to ~1 ms across ~50 warps.

**d) 3D.** The analogue is a closed-mesh validity check: triangle-triangle intersection with exact orientation
(Devillers-Guigue), edge manifoldness and consistent orientation. The same grid-binning broad phase works. For
voxel outputs no validity test is needed: a voxel set is always a valid cubical complex.

---

## 4. Crack-edge polygonization (mask → exact rectilinear rings)

**a) Now.**
- CPU (production): `polygonize.py`. `boundary_edges` (53-65), `link` by argsort + searchsorted (72-89,
  O(E log E)), `ring_labels` min-label pointer jumping (92-99, O(E log E) work, ⌈log2 E⌉+1 rounds),
  `list_rank` Wyllie (102-115, O(E log E)), shoelace area, corner filter, then a Python loop per ring.
- GPU (tests only): `check_mask_chunk` (966-1057). The successor comes from the 2x2 code table (`InitSlots`,
  806-842; left turn at saddles) with no sort, which is good. But the arrays are dense over S slots, there is
  one thread per vertex row, 2·`jump_rounds(S)` rounds (label, then rank), and 3 host round trips.

**b) Best known.**
- Sequential: border following (Suzuki-Abe 1985; crack following, Kovalevsky 1989), O(HW) to find edges plus
  O(E) to trace. Optimal.
- List ranking (orders a ring once it is cut at its start): Wyllie 1979, O(n log n) work, O(log n) depth.
  Anderson-Miller 1990 (randomized) and Cole-Vishkin 1986 (deterministic coin tossing), O(n) work, O(log n)
  depth. Helman-JáJá 2001, a sparse ruling set, practical on SMPs. On GPU: Rehman, Kothapalli, Narayanan,
  "Fast and scalable list ranking on the GPU", ICS 2009 (recursive Helman-JáJá, ~1 s for 32 M nodes on 2009
  hardware; Wyllie loses on large lists) https://faculty.iiit.ac.in/~kkishore/ics152-rehman.pdf . Wei-JáJá,
  IPDPS 2010.
- Labeling cycles: min-label pointer jumping (Shiloach-Vishkin style), O(n log n) work. **But here no cycle
  labeling is needed at all** (next item).

**Key structural fact (enables c).** In the (8,4) pairing, the components' adjacency graph is a tree
(Rosenfeld 1974, "Adjacency in digital pictures"; Kong-Rosenfeld 1989), and the border between an adjacent
foreground component F and background component B is **one** closed crack curve. Every crack edge has a
foreground pixel f on its right and a background pixel b on its left. Hence:
- ring(e) = **outer ring of F = root(f)** when root_bg(b) = surround(F); otherwise **hole ring of B = root_bg(b)**.
  Here surround(F) = the background label of the pixel directly above F's first raster pixel (that pixel is
  background, and it is not enclosed by F: the same argument as the `Holes` kernel, gpucheck.cpp:1379-1391).
- **Canonical start in closed form** (min key (y, x, d) over the ring, as in polygonize.py:130-133):
  - outer ring of F with first raster pixel (x, y): start = vertex (x, y), direction **east**. The vertex
    code there is 8 (only BR set), so east is the only edge leaving it. No vertex of F's ring lies in row y with
    a smaller x, and none lies above row y.
  - hole ring of B with first raster pixel (x, y): start = vertex (x, y), direction **south**. Above b and left
    of b are foreground. The code is 7 or saddle 6; in the saddle case the north edge belongs to the top-left
    pixel's ring, not to B's.
  So per-ring `atomicMin` on keys is unnecessary, and so is the host sort. Ring order by start = sort of
  roots by (y, x), and roots are already in raster order.
- The hole → part grouping (what `build.py:91-96` does with OpenCV) is the foreground label of the pixel above
  the hole's root.

I proved this by hand above. It must be **certified by test**: on the existing goldens, the closed-form starts
must equal `polygonize`'s minimum keys, for every ring.

**c) Recommended CUDA polygonizer (one block per mask, everything in shared memory when E ≤ ~12 K).**

```
Input: bit rows of the mask (already resident: raster output or packed mask) and the fg(8)/bg(4) root
labels from the CCL of section 5 (also resident: the scorer computes them).
1. Codes by bit slicing (warp per vertex row y, lane = word j):
     p = TL = row[y-1] shifted right by 1 with carry, q = TR = row[y-1], r = BL = row[y] shifted, s = BR = row[y]
     outE = s & ~q ; outS = r & ~s ; outW = p & ~r ; outN = q & ~p          (has_edge, gpucheck.cpp:770-777)
     corner = ~(straightH | straightV | empty | full)                        (Q1, Q3, QD)
     deg = popc(outE)+popc(outS)+popc(outW)+popc(outN) per word -> block exclusive scan -> edge base per word
2. Emission: each set bit becomes edge e = base + rank within the word (rank by popc of the lower bits:
   the "compress-store via ballot + popc" idiom). Store (x, y, d), sorted by (y, x, d) for free.
3. Successor by the 2x2 table: the end vertex v' = v + step(d). Its code gives the out edges. Choose the left
   turn ((d+3)&3) at a saddle (codes 6/9), else the unique edge. idx(v', d') = base(word(v')) + popc(lower bits
   of the out masks of v' before d'). O(1), no search, no sort.
4. Ring id and start: ring(e) from the (root_fg(right pixel), root_bg(left pixel)) rule. start(ring) = closed
   form. When no CCL is available, fall back to min-label pointer jumping in shared memory (log2 E rounds).
5. Corner index: weighted list ranking with w = is_corner(pred), cut at start(ring):
     E ≤ 12 K: Wyllie in shared memory (⌈log2 Lmax⌉ rounds, __syncthreads between them), O(E log E) but tiny.
     larger:   sparse ruling set. Rulers = the corners with (x + y) mod s == 0 (deterministic) plus every ring
               start. Each ruler walks sequentially to the next ruler (expected s steps). Wyllie on the rulers.
               Broadcast back. O(E) work, O(s + log E) depth (Helman-JáJá / RHJ).
   Optional corner contraction first: link corner to next corner directly by a find-first-set scan along
   the row (ctz on the corner mask) or column (32x32 bit transpose). Then rank only K ≤ E corners.
6. Area: twice the signed area = shoelace with integer atomicAdd per ring (order-free, exact). For an outer
   ring it also equals 2(|F| + Σ|holes of F|), a free cross-check.
7. Output: ring offsets from a block scan over ring corner counts. Scatter corners. One copy out (or keep
   them on the device).
```

Exactness: every value is an integer that the mask determines uniquely, so the result is identical to
`polygonize` and independent of thread schedule or atomic order.

Speed (estimate): a 256² mask with a ~100 px blob has E ≈ 400-1,000 edges. Steps 1-7 touch 2 KB of bits plus
~16 B per edge, so one block finishes in ~10-30 µs. 512 masks over 170 SMs (RTX 5090) is ~3 waves, ~0.1 ms plus
one launch. The current dense design per 256² mask: 264 K slots x 40 rounds ≈ 10^7 slot visits, at only 257
threads. The work ratio is ≈ 264 K·40 / (1 K·20) ≈ 500x.

**d) 3D: surfels.** The exact 3D boundary of a voxel set is the set of unit faces between foreground and
background voxels (surfels, the 2-cells of the cubical complex; Herman's digital surfaces; Artzy-Frieder-Herman
1981). Each surfel has 4 edges. The neighbour across an edge follows from the 4 voxels around that edge: the
**same 2x2 table** applied in the plane perpendicular to the edge. Joining diagonal foreground at the ambiguous
edge gives (26,6); the vertex case needs Lachaud's tables, below. The surface id is the (foreground root,
background root) pair, exactly as in 2D. This produces exact closed surfaces with **no marching-cubes
ambiguity**.

On marching cubes: MC (Lorensen-Cline 1987) has ambiguous faces and cells. MC33 (Chernyaev 1995) fixes the
topology for trilinear interpolants, not for digital connectivity. For binary voxels the right tool is
Lachaud-Montanvert, "Continuous analogs of digital boundaries: a topological approach to iso-surfaces"
(Graphical Models 62, 2000). It gives MC-style tables consistent with (6,26), (26,6), (6,18), (18,6), and it
proves the triangulated surface keeps digital connectedness and separation.
https://jacquesolivierlachaud.github.io/publication/lachaud-00-a/lachaud-00-a.pdf

Dual contouring (Ju et al. 2002) can be non-manifold; Manifold DC (Schaefer-Ju-Warren 2007) fixes that. Neither is
needed for binary data. Use surfels for exactness and the Lachaud tables when a triangle mesh is wanted. Mesh
simplification (QEM, Garland-Heckbert 1997) comes after, certified by re-voxelization, as in 2D.

---

## 5. Connected components labeling (2D 8/4, 3D 26/6)

**a) Now.**
- Python: OpenCV Spaghetti (`topology._label`, topology.py:75-77; Bolelli et al., TIP 2020), the best
  sequential algorithm. scipy `ndimage.label` for 3D (volume.py:56-57).
- GPU scorer: union-find with `atomicMin` linking (Playne-Hawick 2018). `unite`/`find` at 1262-1282;
  `LabelInit`/`LabelMerge`/`LabelFlatten` at 1302-1347, **one thread per crop row**, iterating pixels
  serially on a byte image. Roots = min pixel index (deterministic; root = first raster pixel, which `Holes`
  relies on). The sequence fg, filter, bg, filter, fg, bg (1611-1617) is 4 labelings per call. Work O(N α)
  per labeling. Depth O(W + union chains). Memory: byte pixels plus 4-byte labels, with stride-W access across
  a warp, so loads do not coalesce.

**b) Best known.**
- Sequential: Spaghetti / BBDT block-based decision trees (Grana, Bolelli et al.), O(N).
- PRAM: Shiloach-Vishkin O(log n) with n+m processors. Optimal work O(n+m) in O(log n) (Gazit 1991,
  randomized; Halperin-Zwick).
- GPU, measured state of the art: **BKE and BUF**, Allegretti, Bolelli, Grana, "Optimized Block-Based
  Algorithms to Label Connected Components on GPUs", IEEE TPDS 31(2):423-438, 2020 (2x2-block union-find /
  Komura equivalence; 2D 8-conn and 3D; the YACCLAB GPU benchmark).
  https://iris.unimore.it/bitstream/11380/1179616/1/2018_TPDS_Optimized_Block_Based_Algorithms_to_Label_Connected_Components_on_GPUs.pdf
  Playne-Hawick, TPDS 29(6), 2018 (pixel-parallel union-find). Komura 2015 (label equivalence).
  **HA4 / FLSL** (Hennequin, Lacassagne, Cabaret, Meunier, DASIP 2018; GTC'19 S9111; ICASSP'21): one warp
  per row, segments (runs), and `__ballot_sync`/`__popc` compress-store for run detection. This is the
  closest match to this repo's warp-per-row, bit-word design.
  https://largo.lip6.fr/~lacas/Publications/DASIP18_CCL_GPU.pdf

**c) Recommendation.**
1. **Run-based union-find on the bit words** (HA4/FLSL style). Do not copy to bytes.
   - Warp per row. Run starts = `w & ~((w << 1) | carry_in)`, run ends likewise, for foreground and for
     `~w` (background). One union-find node per run, indexed by the **run's first pixel index**. Roots stay
     "min pixel index = first raster pixel", so `Holes` and determinism are unchanged.
   - Merge row y with row y-1. Foreground 8-connectivity: runs overlap after widening by 1 on each side.
     Background 4-connectivity: exact overlap. Use the existing `unite` (atomicMin retry).
   - Work O(#runs) unions instead of O(#pixels); memory ≈ bits plus runs. Coalesced word loads.
2. **Per-component Euler number replaces the final background labeling.** Every 2x2 window's foreground
   pixels are mutually 8-adjacent, so the window belongs to one component k. Then
   χ8(k) = Σ_windows of k (q1 − q3 − 2qD)/4 (Gray 1971, applied to k alone, which is valid because no window
   mixes two components), and **holes(k) = 1 − χ8(k)**. One pass: a bit-sliced code per 32 vertices,
   `atomicAdd` into the root's counter (integer, so the order does not matter).
3. **Skip the relabel when the filter changes nothing.** `Filter` sets a per-object "changed" flag. With no
   change, the signature comes from the first foreground labels plus step 2. This mirrors
   `topology.signature`:161.
   Common case: 1 foreground + 1 background labeling + 1 Euler pass, instead of 4 labelings + `Holes`.
4. Medium term: BKE (2x2 blocks) for foreground 8-connectivity on large crops. It has 4x fewer nodes and is
   the measured fastest on YACCLAB.
Expected speed (estimate): the signature stage gets 5-20x faster (4 → ~2 labelings, coalesced word access,
runs instead of pixels). Absolute numbers need a measurement on the RTX 5090. Use `compute-sanitizer
--tool racecheck` on the new kernels (section 9).

**d) 3D.** Runs along x. Merge with the runs of the 4 face-neighbour rows (6-connectivity, background) or of the 9
rows of the previous slab plus the row above (26-connectivity, foreground; widen by 1). BUF/BKE 3D (same TPDS paper)
are the measured alternatives. Roots = min voxel index stays deterministic. Cavity attribution: the voxel at z−1
above a cavity's first raster voxel is foreground of the enclosing component (same argument as 2D, using the 3D
(26,6) adjacency tree, Kong-Rosenfeld 1989).

---

## 6. Euler number and Betti numbers (2D and 3D, including b1 in 3D)

**a) Now.**
- 2D: Gray bit-quads in numpy (`bitquad.py:63-68`), χ8 = (n1 − n3 − 2nD)/4. The GPU computes quads inside
  `InitSlots` (gpucheck.cpp:814-840) with byte loads. b0 and b1 per component come from two labelings
  (`topology.instance_topology`, 105-125). `certify.py:42-58` cross-checks three ways (Euler, labels, rings).
- 3D: `volume.euler` (34-48) computes V − E + F − C of the closed-cube complex with OR-reductions (correct);
  `betti3d` (51-59) gives b0 by 26-CCL, b2 = bounded 6-connected background components, b1 = b0 + b2 − χ. CPU
  only.

**b) Best known.** χ is a sum of local terms: O(N) work, O(log N) depth (a reduction), exact. Gray 1971 for 2D;
the 2x2x2 configuration LUT for 3D (Lobregt-Toet-Groen 1980; Toriwaki-Yonekura 2002). b0 and b2 by CCL, then b1
by Euler-Poincaré. No faster asymptotic method exists, since the input must be read.

**c) Recommendation (CUDA).**
- 2D: bit-sliced counts. With p, q, r, s as in section 4, `odd = p^q^r^s`, `ge2 = (p&q)|(p&r)|(p&s)|(q&r)|(q&s)|(r&s)`:
  n1 = popc(odd & ~ge2), n3 = popc(odd & ge2), nD = popc((p&s&~q&~r)|(q&r&~p&~s)). About 20 operations per 32
  vertices. The per-component version is section 5 c2.
- 3D global χ: on x-packed words, V = popc(OR of the 8 cubes around each vertex), E and F likewise, C = popc(v).
  This is exactly `volume.euler` on 32 voxels per word.
- **Per-component 3D b1, exact:** every cell (vertex, edge, face, cube) belongs to the component of any present
  cube that contains it. All cubes containing one cell are mutually 26-adjacent, so the attribution is
  unambiguous. Then χ_k = V_k − E_k + F_k − C_k by integer atomics, and **b1_k = 1 + b2_k − χ_k** (b2_k from
  cavity attribution, section 5d).
- An independent 3D certificate (the analogue of `certify`): the total genus of the boundary surfaces of
  component k equals b1_k. This is "half lives, half dies" for compact 3-manifolds in R^3; χ(∂M) = 2χ(M) holds
  once the (26,6) boundary is made manifold by Lachaud's construction. From the surfel mesh of section 4d:
  g = (2 − χ(S))/2 per closed surface S.

**d) Persistent homology on cubical complexes (only needed if graded or soft masks enter the reward).**
- 0-dim persistence = Kruskal / union-find over cells sorted by value (the elder rule), equivalently the MST of
  the grid graph with edge weight max(f(u), f(v)). On GPU: deterministic radix sort plus parallel Borůvka with
  index tie-breaks.
- 2D, dimension 1: by **duality**, PH1 of an image equals PH0 of the complement under the dual connectivity
  (Garin, Heiss, Maggs, Bleile, Robins, arXiv:2005.04597, 2020). So all 2D persistence is union-find.
  https://arxiv.org/abs/2005.04597
- 3D: PH0 and PH2 (by duality) are union-find. Only PH1 needs matrix reduction. CubicalRipser (Kaji, Sudo,
  Ahara, arXiv:2005.12692; CPU, Z/2, 1-4D, the fastest published) https://arxiv.org/abs/2005.12692 . GUDHI's
  cubical complex. Ripser++ (GPU) is Vietoris-Rips only and does not apply.
- A cheap exact alternative to PH1 diagrams: **Betti curves**. Do one GPU pass for the Euler characteristic
  curve χ(t) (each cell is born at the max of its voxels: histogram ±1 at birth values, then a prefix sum). Then
  b1(t) = b0(t) + b2(t) − χ(t), with b0 and b2 from the union-find barcodes. All integers, all exact.

---

## 7. Boundary F-score and distance transforms

**a) Now.**
- Reward F5: the exact EDTs of the experts' union and core are **precomputed offline** (scipy Maurer,
  `cache.py:130-135`, `_band` 185-189), stored float32, resident on the device. Online `ScoreRows`
  (gpucheck.cpp:1411-1442) only looks up distances at error pixels (`ctz` loops on spill/miss words) and sums
  them as integer multiples of 2^-23 (`units`, 1407; `sum_value`, 1469-1473). That is O(HW/32 + #error pixels)
  per object, exact. **This is already the optimal online design.** No EDT per rollout.
- Eval: NSD uses `ndimage.distance_transform_edt` on the corner grid per sample (`metrics._surface`,
  199-206), on the CPU.

**b) Best known.** Exact EDT is O(N) sequentially in any dimension: Maurer-Qi-Raghavan TPAMI 2003; Meijster et al.
2000; Felzenszwalb-Huttenlocher (ToC 8, 2012; lower envelope of parabolas). On GPU, exact: **PBA**, Cao, Tang,
Mohamed, Tan, I3D 2010, "Parallel Banding Algorithm to compute exact distance transform with the GPU": optimal
linear work, 2D and 3D. https://www.comp.nus.edu.sg/~tants/pba.html (PBA+ later). Jump flooding (Rong-Tan 2006)
is **approximate**; JFA+1/1+JFA reduce errors but give no guarantee, so reject it for this project. The separable
FH/Meijster passes (column pass, then a row lower-envelope pass) parallelize trivially: O(N) work, O(W + H) depth.

**c) Recommendation.**
- GPU EDT (for eval metrics and any future online band): separable Meijster/FH with **integer squared distances**
  and integer intersection tests (cross-multiplied, no division), one thread per column then one per row
  (crops ≤ 512² make the O(W + H) depth irrelevant). Then `__double2float_rn(__dsqrt_rn((double)d2))`. CUDA's
  double sqrt is IEEE correctly rounded, and binary64 → binary32 double rounding is innocuous for sqrt, as the
  cache docstring already argues. The value then equals scipy's float32 bit for bit. The EDT value is unique, so
  nearest-feature ties do not matter. Use PBA for large 3D volumes.
- Radius tests (a boundary F at tolerance τ in the DAVIS/BSDS style) need no EDT: dilate the boundary bitmap by
  a digital disk = OR of shifted bit rows, O(τ·HW/32) word operations, exact.
- Keep the offline cache EDT on the CPU (it runs once per item), or port it to the same GPU EDT so the whole
  scorer runs on CUDA.

**d) 3D.** The separable EDT generalizes to 3 passes. **Anisotropic CT spacing** (e.g. 0.7 x 0.7 x 2.5 mm): to stay
exact, scale the spacings to integers (µm) so squared distances stay int64. scipy's `sampling` uses floats and
gives no such guarantee. Surface Dice in 3D uses surfel areas instead of surfel lengths (DeepMind
surface-distance does the same).

---

## 8. Polygon simplification

**a) Now.** `targets/simplify.py`. A potrace-style **optimal** min-# polygon: the shortcut graph from exact integer
cones, with lossless DSS or L∞-τ admissibility (`_graph`, 107-155, O(n·L)), then the lexicographic (count,
penalty) shortest cycle, one run per start in the anchor's window (`_best_cycle`, 158-198, O(c·n·L), worst
O(n^3)). CPU numba, offline (target building). Exactness comes from re-rasterizing and certifying on the checker.

**b) Best known.**
- Douglas-Peucker: O(n log n) expected, O(n^2) worst. Hershberger-Snoeyink 1992 (path hull), O(n log n) worst
  case. Neither is optimal in vertex count, and neither preserves topology (Saalfeld 1999 fixes that).
- Visvalingam-Whyatt: O(n log n) with a heap. Area-based, not optimal.
- Imai-Iri min-# (graph shortest path): O(n^2 log n), improved by **Chan-Chin** (IJCGA 6(1), 1996) to O(n^2) for
  min-# and O(n^2 log n) for min-ε.
- Digital curves: DSS recognition is linear (Debled-Rennesson-Reveillès 1995). Maximal segments in O(n)
  (Feschet-Tougne 1999; Lachaud et al.). The **minimum length polygon** in O(n) (Provençal-Lachaud, DAM
  159(18), 2011) is unique and convergent, but its vertices sit on half-pixel boundaries: it can tie at pixel
  centres, so it is not certifiably lossless here.
- Topology-preserving: de Berg-van Kreveld-Schirra 1998; Saalfeld 1999.

**c) Recommendation.** The current method is the right one: optimal min-# under exact admissibility, certified.
Two improvements:
1. If admissibility is **interval-hereditary** (lossless DSS is: every sub-segment of a DSS is a DSS), the min-#
   cyclic cover reduces to farthest reach f(i) plus **binary lifting** (pointer doubling on f). That is
   O(n log n) work and O(log n) depth, exact (greedy is optimal from a fixed start; take the min over starts in
   one greedy window). Use it to bound the count. Then run the (count, penalty) DP only among starts that reach
   that count. This is GPU-friendly: one thread per i builds f(i); one block per ring.
   The current `ok[l]` checks each endpoint separately, so the L∞-τ test may not be hereditary. Verify this
   before using the shortcut for τ > 0.
2. Bound the worst case: cap L, or use the potrace anchor trick with the fewest-crossing anchor (already done).
   Port only if targets are rebuilt often. The owner's CUDA rule names checker, scorer and polygonizer, not
   simplification.

**d) 3D.** Mesh decimation (QEM) constrained to stay within the voxel band, certified by re-voxelization. Or keep
per-slice 2D polygons (the current LIDC design).

---

## 9. Set-valued IoU against k experts

**a) Now.** `ScoreRows` (1411-1442): per row, k AND-popcount sums against the resident rater bit planes, plus vote,
spill and miss. `ObjTerms` (1476-1524): row sums, IoU_tc = max IoU over raters with an equal signature (else IoU vs
vote), F5, g_area, topo_freq. Work O(k·HW/32). Depth O(W/32 + H) per object. Integer counts, then one IEEE divide,
so the result is exact. Eval: pairwise IoU by float64 matmul on unpacked crops (`metrics._pairwise_iou`, 172-186),
Hungarian on CPU.

**b) Best known.** A lower bound of Ω(k·HW/w) word reads (every plane must be read); the current kernel attains
it. Any integer reduction order is exact.

**c) Recommendation.**
- Keep it. Optional: **bit-sliced vote counters** (carry-save adders over the k planes, Harley-Seal style,
  ⌈log2(k+1)⌉ planes). They give the "≥ t raters" level sets for all t in one pass, so IoU against every
  consensus level (the full egg-yolk family, STAPLE-like) costs O(k·HW/32).
- Pairwise IoU (GED, diversity) on GPU: popcount Gram matrix `inter_ij = Σ popc(A_i & A_j)`, or 1-bit tensor
  cores (`mma ... .b1 ... .and.popc`, sm_80+; **check the PTX ISA table for sm_120 before relying on it**, which
  I could not confirm). The result is an integer, so it is exact.
- Hungarian for small K x k stays on the host (microseconds).

**d) 3D.** The same kernel over z-stacked planes: the set IoU of volumes is a sum of per-slice integer counts.

---

## 10. Exactness doctrine (why fast, unordered GPU algorithms stay bit-reproducible)

1. Every geometric predicate is an integer sign (orientation, crossing column, parity) at the precision tier
   chosen per shape (i32/i64/i128). No floating point is used in decisions.
2. Every structural output (masks, labels as min indices, rings, starts, ranks, areas, Betti numbers) is a
   **unique function of the input**. Any correct algorithm, under any schedule or atomic order, returns the same
   integers. Nondeterministic `atomicAdd`/`atomicMin` on integers is therefore safe. Only list **orders** built
   with `atomicInc` must be canonicalized (sort or closed form), as the repo already does.
3. The only reals are correctly rounded IEEE operations on exact integers (division for IoU, sqrt for the EDT) or
   exact integer sums of dyadic values (band depth in units of 2^-23). Build without FMA contraction, as the repo
   does.
4. Testing: certify every new kernel against the CPU oracle on the goldens. Run the three-way `certify` check on
   the GPU outputs (χ from quads = b0 − b1 = #outer − #holes; Σ ring areas = popcount; corners = n1 + n3 + 2nD).
   Run `compute-sanitizer` racecheck, synccheck and memcheck (the NVIDIA catalog skill
   `compute-sanitizer-diagnose-cuda-bugs`) on the union-find and ring kernels.

---

## 11. Proposed design (system-design + ADR, short form)

```
         device-resident pipeline, one stream per lane (as now)
 polygons ──► [raster: warp/row, ballot]──► bit words ──┬──► [validity: grid broad phase, exact classify,
                                                        │     warp point-in-ring parity]──► defects
 masks ─────────────────────────────────► bit words ──┤
                                                        ├──► [CCL: run-based UF on words, fg8 + bg4]──► roots
                                                        │        │
                                                        │        ├──► [per-component χ8 by quads] ─► holes(k)
                                                        │        └──► [polygonizer: emit, 2x2 succ,
                                                        │              closed-form starts, rank, scatter]──► rings
                                                        └──► [score: k popc planes, EDT lookups]──► terms
 3D: same modules, with runs along x, surfels instead of cracks, cell-attributed χ_k, cavities by attribution.
```

**ADR-1: ring identity and start from CCL roots (status: proposed).** Options: (A) the current min-label pointer
jumping on dense slots; (B) the same on compacted edges in shared memory; (C) CCL-root identity plus closed-form
starts. C has the lowest work (O(E), no log rounds), no host sort, and reuses scorer labels. Its cost is a
dependency on CCL and a new invariant that tests must certify. Recommend C, with B as a fallback when no labels
exist.

**ADR-2: CCL on bit words with runs (status: proposed).** Options: (A) the current row-serial byte union-find;
(B) pixel-parallel Playne-Hawick; (C) run-based warp-per-row (HA4/FLSL); (D) BKE 2x2 blocks. C fits the existing
bit-word layout and the ballot idioms, cuts union work to the number of runs, and keeps min-index roots. D is
faster on large dense images. Recommend C now, D if profiling demands it.

**Module seam (codebase-design).** One geometry core with a small interface: `rasterize`, `label`, `trace`,
`topology`, `distance`. The implementation stays in `gpucheck.cpp`. The CPU build is the second adapter at the
same seam, already used by the tests, so the seam is real. Tests cross only this interface, comparing against the
CPU oracle.

**Action list (ordered by value).**
1. GPU polygonizer per section 4c; switch `targets/build.py` to it (F1, F2, F3). Certify against
   `polygonize` on all targets (hashes unchanged, or re-hashed if ADR-1 changes the ring order).
2. CCL: changed flag, per-component χ8, then run-based labeling (F4).
3. Parallelize `ShapeTotals` with warp parity; grid broad phase for pairs; i32 orient tier (F5-F7).
4. GPU exact EDT (eval NSD, then cache); 3D: runs-based 26/6 CCL, cell-attributed χ_k, surfel certificate (F10).

---

## References (verified this session unless marked "standard")

- Allegretti, Bolelli, Grana, "Optimized Block-Based Algorithms to Label Connected Components on GPUs", IEEE TPDS
  31(2):423-438, 2020. https://iris.unimore.it/handle/11380/1179616
- Hennequin, Lacassagne, Cabaret, Meunier, "A new direct connected component labeling and analysis algorithm for
  GPUs", DASIP 2018. https://largo.lip6.fr/~lacas/Publications/DASIP18_CCL_GPU.pdf ; GTC 2019 S9111
  https://developer.download.nvidia.com/video/gputechconf/gtc/2019/presentation/s9111-a-new-direct-connected-component-labeling-and-analysis-algorithm-for-gpus.pdf ;
  ICASSP 2021 (ballot/popc compress-store) https://largo.lip6.fr/~lacas/Publications/ICASSP21_CCA_GPU.pdf
- Playne, Hawick, "A new algorithm for parallel connected-component labelling on GPUs", IEEE TPDS 29(6), 2018 (standard).
- Bolelli, Allegretti, Baraldi, Grana, Spaghetti labeling, IEEE TIP 29, 2020 (cited in topology.py).
- Rehman, Kothapalli, Narayanan, "Fast and scalable list ranking on the GPU", ICS 2009, 235-243.
  https://faculty.iiit.ac.in/~kkishore/ics152-rehman.pdf
- Wyllie 1979; Anderson-Miller 1990; Cole-Vishkin 1986; Helman-JáJá 2001; Wei-JáJá IPDPS 2010 (standard).
- Shamos-Hoey FOCS 1976; Bentley-Ottmann 1979; Chazelle-Edelsbrunner JACM 1992; Balaban SoCG 1995; Chazelle 1991
  linear triangulation (standard).
- Goodrich, "Intersecting line segments in parallel with an output-sensitive number of processors", SIAM J. Comput.
  20, 1991. https://epubs.siam.org/doi/10.1137/0220047
- Cao, Tang, Mohamed, Tan, "Parallel Banding Algorithm to compute exact distance transform with the GPU", I3D 2010,
  83-90. https://www.comp.nus.edu.sg/~tants/pba.html
- Maurer, Qi, Raghavan TPAMI 2003; Felzenszwalb-Huttenlocher ToC 2012; Meijster et al. 2000; Rong-Tan JFA 2006 (standard;
  the first two are already cited in cache.py).
- Lachaud, Montanvert, "Continuous analogs of digital boundaries: a topological approach to iso-surfaces", 2000.
  https://jacquesolivierlachaud.github.io/publication/lachaud-00-a/lachaud-00-a.pdf
- Lorensen-Cline 1987; Chernyaev MC33 1995; Ju et al. DC 2002; Schaefer-Ju-Warren 2007; Garland-Heckbert 1997 (standard).
- Rosenfeld, "Adjacency in digital pictures", Information and Control 26(1):24-33, 1974; Kong-Rosenfeld CVGIP 1989.
- Gray 1971; Lobregt-Toet-Groen 1980; Toriwaki-Yonekura 2002 (3D Euler LUTs) (standard).
- Kaji, Sudo, Ahara, "Cubical Ripser", arXiv:2005.12692. https://arxiv.org/abs/2005.12692
- Garin, Heiss, Maggs, Bleile, Robins, "Duality in persistent homology of images", arXiv:2005.04597.
  https://arxiv.org/abs/2005.04597
- Chan, Chin, "Approximation of polygonal curves with minimum number of line segments or minimum error", IJCGA 6(1):59-77,
  1996. https://hub.hku.hk/handle/10722/88918
- Provençal, Lachaud, "Two linear-time algorithms for computing the minimum length polygon of a digital contour", DAM
  159(18):2229-2250, 2011. https://hal.archives-ouvertes.fr/hal-00943824
- Hershberger-Snoeyink 1992; Visvalingam-Whyatt 1993; Saalfeld 1999; de Berg-van Kreveld-Schirra 1998;
  Debled-Rennesson-Reveillès 1995; Feschet-Tougne 1999 (standard; DSS refs already in simplify.py).
- PTX `mma` b1 `.and.popc` needs sm_80+ (LLVM review https://reviews.llvm.org/D105384 ). sm_120 support is unconfirmed: check
  the PTX ISA "warp-level matrix instructions" table.
