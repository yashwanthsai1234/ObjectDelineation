# Phase D test cases, set 1 (written from D.spec.md only)

## Skills used
- nvidia-skills:nvidia-skill-finder: loaded; no NVIDIA catalog lookup run (task is test writing, not skill install). Changed nothing.
- engineering:testing-strategy: pyramid split (unit / property / integration / GPU / e2e) and "gaps" thinking.
- mattpocock-skills:tdd: tests only at spec seams (`rings.trace`, `gc.*`, `objdel.checker.*`, `metrics`, `edt`); expected values come from oracles, literals or the spec, never recomputed by the code under test (no tautologies).
- Other listed skills (code-review, debug, system-design, architecture, documentation, huggingface-best) not invoked: no design, review or HF content in this task.

## Conventions
- ORC = `tests/oracles/{polygonize_oracle,bitquad_oracle,topology_oracle}.py` (spec B1, B17). `gc` = `objdel._gpucheck`.
- Level tags: U unit, P property, I integration, G GPU (marker `gpu`), E end-to-end, S sanitizer, F perf gate.
- Tests with no marker run on the CPU build (Mac). "BOTH" = also run with `-m gpu` on the box.
- Hand-computed expectations below are literals, derived from the spec's edge-key rules (B2): key = (y*(W+1)+x)*4+d, d: 0E 1S 2W 3N; foreground lies right of every edge; outer ring area > 0.
- Fuzz seeds are fixed. Data-backed tests skip when `data/` is missing (fixtures fixture pattern in tests/conftest.py) and run Python with `-I` semantics.

---------------------------------------------------------------------------------------------------------------

## A. Polygonizer (`rings.trace`, `gc.check_masks`) - B1, B2, B5, B7

### T01 single pixel (U, BOTH)
- Setup: mask 1x1 = [[1]].
- Input: `trace([m])`.
- Expect: one ring, pts == [(0,0),(1,0),(1,1),(0,1)], area == 1, part == 0, hole False; quads: n1=1, others 0; b0=1, b1=0.
- Catches: start-corner or orientation error (clockwise on screen, east edge first), off-by-one on (W+1) key grid.

### T02 empty / degenerate shapes (U, BOTH)
- Input: masks of shape (0,0), (0,5), (5,0), (3,3) all-zero, (1,1) zero; one batch.
- Expect: 5 `MaskRings`, each `rings == []`, quads == (0,0,0,0,0), (b0,b1)==(0,0). No exception. Output equals oracle.
- Why: B5; zero-size grids give zero-length kernels / zero-size allocs (CUDA launch with 0 blocks is an error if not guarded); `start` offsets must stay consistent so later masks in the batch are not shifted.

### T03 empty mask between non-empty masks (U, BOTH)
- Input: batch [A=2x2 full, E=(0,0), B=3x1 full]. Compare each to `trace([x])` alone.
- Expect: identical per mask; `ring_start`/`pts_start` monotone, ring_start[-1]==R, pts_start[-1]==P.
- Catches: prefix-sum bug when a mask contributes zero rings.

### T04 full mask, sizes 1xN, Nx1, NxN (U, BOTH)
- Input: full masks of 1x7, 7x1, 5x5, 33x1, 1x33, 32x32, 31x31, 33x33.
- Expect: one ring of exactly 4 corners (0,0),(W,0),(W,H),(0,H); area W*H; quads n4 == (H-1)*(W-1) (interior), n2 per bit-quad oracle; equals oracle.
- Catches: word-boundary (32/64 bit packing) errors at sides 31/32/33; border treated as background.

### T05 2x2 diagonal pair: 8-connectivity (U, BOTH)
- Input: [[1,0],[0,1]].
- Expect: ONE foreground component (b0=1, b1=0). Two outer rings?  Per B2/B4 (outer rings == b0) the two pixels form ONE ring that touches itself at the saddle vertex: 1 ring, 8 corners, area 2; both ring parts == 0. Quads: nD=1.
- Counter case [[0,1],[1,0]] same with anti-diagonal. Expect also 1 ring, 8 corners.
- Catches: saddle resolved wrongly (4-conn join gives 2 rings, b0=2). The exact corner list must equal the oracle (left turn at the saddle).

### T06 diagonal background gap: 4-connected background (U, BOTH)
- Input: 3x3 ring of ones with centre 0 and one corner pixel of the ring zeroed so background is only diagonally linked to outside: [[1,1,1],[1,0,1],[1,1,0]]... (centre bg touches outside bg only via diagonal at bottom-right).
- Expect: b1 == 1 (hole NOT merged with outside under 4-conn background), one hole ring area == -1, outer ring area == 8. Matches oracle.
- Catches: background treated as 8-connected (hole disappears, b1==0).

### T07 hole in a hole (nested squares) (U, BOTH)
- Input: 11x11: ring of width 1 at border, inside 0 band, ring of ones at inset 2, centre 0 pixel at (5,5)... 4 nested levels.
- Expect: oracle equality; rings sorted by raster index of first east edge; outer rings area>0, holes<0; each hole ring `part` equals the part of the surrounding outer component; the inner island has its own part (rank by first pixel).
- Catches: part assignment for holes vs islands; ordering.

### T08 checkerboard (U, BOTH)
- Input: 16x16 and 17x17 checkerboards, both phases.
- Expect: all ones diagonal-connected -> b0 == 1; holes (4-connected bg components that are single pixels, bounded): b1 == number of interior zero cells; sum of ring areas == popcount; `b0 - b1 == chi8`; oracle equal.
- Catches: saddle logic at every vertex; ring linking when a vertex has two out-edges (successor choice must be left turn consistently).

### T09 touching borders, all four sides (U, BOTH)
- Input: 6x6 with ones along only the top row; only bottom row; only left column; only right column; four corners only (diag-free); plus rim-hole mask (frame of ones on the border with an interior zero block).
- Expect: oracle equality; frame case has 1 outer (area 36) + 1 hole (area -16); one-row cases have one 4-corner ring.
- Catches: out-of-bounds reads at x=-1/y=-1/x=W/y=H (image border must read as background; sanitizer-relevant).

### T10 offsets (U, BOTH) - B7
- Input: `polygonize(m, (ox, oy))` for m = L-shaped mask, offsets (0,0), (5,7), (-3,-2), (4096,4096), (0,0) with trace offsets list of 3 different offsets in one call.
- Expect: pts == pts0 + [ox, oy]; area, part, order unchanged; `polygonize(m,o) == trace([m],[o])[0].rings`; `trace(masks, None)` equals offsets all (0,0). Mismatched `offsets` length raises ValueError.
- Catches: offset applied to key order (it must NOT change ring order or start corner), applied twice, applied to wrong mask in batch.

### T11 non-binary values (U)
- Input: uint8 mask with values {0, 1, 2, 255}, int64 mask {0, -1, 7}, bool mask, float mask {0.0, 0.5}.
- Expect: nonzero = foreground, identical to bool(mask) result. Non-2-D input (1-D, 3-D) raises ValueError. Non-contiguous (transposed / strided slice) input equals contiguous copy.
- Catches: packing that only checks bit 0 (2 -> background bug); strided buffer passed raw to C++.

### T12 maximum size 4096x4096 (G, plus U on CPU with sparse mask)
- Input: 4096x4096 noise p=0.5; also 4096x4096 full; also 4096x1 and 1x4096; also 4097x1 -> must raise.
- Expect: noise result == oracle (CPU build may use the oracle on a 1024 crop if oracle too slow; the full 4096 case runs on the GPU box and compares to CPU build bits); full -> 1 ring of 4 corners; 4097 side raises ValueError, not crash.
- Catches: int32 overflow in key (y*(W+1)+x)*4+d: 4097*4097*4 = 67,141,636 fits, but N*5 etc. at batch level may not; scratch sizing.

### T13 noise density sweep fuzz vs oracle (P, BOTH)
- Input: seeds 0..49, sizes random in [1, 64]x[1, 64] plus fixed sizes 31, 32, 33, densities {0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95}.
- Expect: for every mask `trace` rings (pts bytes, area, part, order) and quads equal oracle exactly; hold invariants (T14).
- Why: G1 families; densities near 0.5 maximise saddles and holes.

### T14 invariants on every traced mask (P, BOTH)
Invariants (run inside T13, T15, T16 as a helper `check_invariants(mask, mr)`):
1. every ring is closed rectilinear: consecutive corners differ in exactly one coordinate; no two consecutive collinear corners with same direction (corners are true corners); first corner is the min-key edge's start.
2. sum(ring.area) == popcount(mask).
3. #outer == b0 == number of cv2 8-components; #holes == b1; b0 - b1 == chi8 from quads.
4. rings sorted strictly ascending by first-east-edge raster index.
5. sum(len(pts)) == bit-quad corners (n1 + n3 + 2*nD ... use oracle value).
6. foreground right of every edge: for each unit edge (x,y)->(x+dx,y+dy), the pixel to its right (by direction) is 1 and left pixel is 0 (outside-of-image = 0).
7. part: `part+1 == cv2.connectedComponents(m, connectivity=8)[1][py, px]` where (px,py) = pixel right of first edge (B3); parts of rings form {0..b0-1}.
- Catches: any family of ring-link or label bug independent of oracle (oracle may share a misconception with the implementation).

### T15 structured families (P, BOTH)
- Families, each at 3 sizes incl. 31/32/33: blobs with holes (random discs minus discs, seed fixed), comb (teeth width 1, gap 1, depth up to H-1), spiral (1-px arm, 1-px gap; deepest nesting, longest single ring: corners ~ O(W*H/2)), diagonal lines (both diagonals, thickness 1), stripes (rows, columns, period 2), nested squares (up to 10 levels).
- Expect: oracle equality + T14.
- Catches: very long ring (spiral) breaks pointer-jumping round count (log2 of ring length: ceil(log2 n)+1 rounds must cover n = ring edges, not n = pixels); comb gives many corners in one ring.

### T16 adversarial: single huge ring / many tiny rings (P, G)
- (a) 512x512 spiral: one ring with ~65k+ corners. (b) 512x512 isolated pixels on a 2-px lattice (every other row/col): 16,384 outer rings of area 1, each ring 4 corners, b0 = 16,384. (c) 512x512 lattice of single-pixel holes inside full mask: b1 = 255*255 or as computed by oracle.
- Expect: oracle equal; ring order = raster order of first east edge (pixels in lattice -> ring k at pixel k in raster order).
- Catches: ring sort not stable / sort by wrong key; list-ranking depth limit; part rank overflow.

### T17 ring-count and vertex limits (U, BOTH) - B16
- Input to `gc.check_polygons`/limits only apply to polygons ("at most 10 rings, 1000 vertices"): build a polygon with 11 rings, and one ring with 1001 vertices.
- Expect: ValueError (Python) naming the limit; with exactly 10 rings and exactly 1000 vertices, accepted. Mask side 4097 raises ValueError; side 4096 accepted; negative dims in `hw` raise ValueError; `start` not monotone / start[-1] != len(bits) raises ValueError; hw shape mismatch with n raises ValueError.
- Catches: boundary off-by-one on limits (>= vs >), memory-unsafe acceptance of malformed `start`.

---------------------------------------------------------------------------------------------------------------

## B. Batch invariance and determinism - B6, G6

### T20 chunk / scratch split invariance (I, CPU; BOTH)
- Setup: 40 mixed masks (empty, 1x1, noise, spiral, checkerboard, full, 1xN; sizes 1..200), fixed seed.
- Input: `gc.check_masks(..., chunk_pixels=c, max_scratch_bytes=s)` for (c, s) in {default, (4096, 65536), (1, 1<<30)?, (4096, 1<<30), (1<<25, 65536)}. If a cap is smaller than the largest single mask needs, expect ValueError (never silent truncation) - add that as T21.
- Expect: all output arrays (quads, ring_start, pts_start, ring_area, pts, ring_part, betti) byte-identical across settings and equal to the one-mask-at-a-time results concatenated.
- Catches: chunk boundary bugs (ring crossing chunk, key base offset), item index confusion within a chunk.

### T21 scratch too small for one item (U)
- Input: single 512x512 noise mask, `max_scratch_bytes=1024`.
- Expect: ValueError (or documented RuntimeError), no crash, and a following valid call still returns right results (engine state not corrupted).
- Catches: unreported partial processing, leaked arena after failure.

### T22 item order permutation (P, BOTH)
- Input: 30 masks; random permutation pi (seeds 0..9). Output for permuted batch, un-permuted per mask.
- Expect: per-mask rings, part, quads, betti identical to the original order.
- Catches: shared scratch indexed by global edge index leaking neighbours (e.g. ring ids unique per batch vs per mask).

### T23 alone vs in batch (P, BOTH)
- Input: each of 40 masks alone vs inside a batch of 40 vs inside a batch of 40 + 200 random neighbours.
- Expect: byte-identical `Ring.pts`, area, part.

### T24 repeat determinism (G, CPU too with 20)
- Input: one fixed batch (60 masks incl. spiral, 512x512 noise) run 20 times on the same process, and in a fresh process.
- Expect: byte-identical output arrays (compare `.tobytes()`) for `check_masks`, `topology`, `certify_polygons`, `score_polygons`, `eval_sets`, `add_tables_from_raters` (G6). 
- Catches: atomics ordering, racy label min-propagation, uninitialised scratch (compare after pre-filling the pool with 0xFF via a prior large call of different shape: `gc.trim(0)` and not-trim variants).

### T25 stale-scratch poisoning (I, BOTH)
- Input: call 1: large noise batch (fills arena). Call 2: small batch of sparse masks. Call 3 = call 2 after `gc.trim(0)`.
- Expect: call 2 == call 3 == oracle.
- Catches: kernels that assume zeroed scratch (initcheck-style bug detected without the sanitizer).

---------------------------------------------------------------------------------------------------------------

## C. Part labels and Betti - B3, B4

### T30 part labels vs cv2 (P, BOTH)
- As T14.7 on 300 random masks (sizes 2..80, density 0.2-0.8) and on the 4 fixed oracles of: nested squares, two islands inside a ring (parts 1 and 2 for islands, hole part = ring's part), checkerboard (all one part).
- Expect: exact B3 equality; `ring_part` int32; part is per-mask (not offset by batch position).
- Catches: part = raster rank of the *ring* rather than the *component*; hole assigned to the wrong enclosing component when two components nest differently (island inside hole of component X: island's part != X's).

### T31 component-with-island worked case (U, BOTH)
- Mask 9x9: outer frame (ring of ones, border), one hole (inside), island 1 pixel at (4,4) inside hole.
- Expect: rings in order: frame outer (part 0, area 81-... compute: 81 area frame outer ring = 9*9=81), frame hole (area -49, part 0), island (area 1, part 1). b0=2, b1=1.
- Catches: hole/island part swap.

### T32 betti = (b0, b1) and Euler (P, BOTH)
- On T13/T15 corpus: `betti[m] == (n_components_8, n_bounded_bg_components_4)` computed with cv2 on the padded mask (independent of quads); `b0-b1 == chi8`.

---------------------------------------------------------------------------------------------------------------

## D. Certify - B8, G16

### T40 certify on normal masks (U, BOTH)
- Input: certify on 200 masks (T13 corpus) with min_size in {0, 1, 2, 5, 50}.
- Expect: every Report has certified True and fields equal to the pre-Phase-D values stored as a golden (b0, b1, area, perimeter, corners, chi8) - golden generated from the archived implementation, kept in tests/goldens; `failures == []`.
- Catches: field meaning drift (perimeter, corners).

### T41 certify_many == [certify(m) ...] (U)
- Input: 30 masks incl. empty and 1x1. Expect list equality of dataclass fields, same order.

### T42 certify independence mutation (U/I) - G16
- Setup: a test-only build/flag where the ring labeller joins foreground via 4-adjacency only.
- Input: [[1,0],[0,1]].
- Expect: with the mutant `certified == False` (rings outer=2 != b0=1 from topology, with `failures` containing a "rings:" line); with normal build certified == True.
- Also: mask with a hole touching diagonally ([[1,1,1],[1,0,1],[1,1,0]]-like): the mutant must also fail.
- Catches: certificate that is tautological (b0,b1 taken from the same labeller as rings).
- If the mutant build flag is unavailable on the box, run the mutation on the CPU build (compiled twice).

### T43 certify empty mask (U)
- Input: all-zero 5x5; min_size larger than the only component (e.g. 1-pixel mask, min_size 5).
- Expect: Report b0=b1=0, area=0, corners=0, certified True, no exception (empty crop `bbox` of zero size must not crash).

### T44 certify min_size filtering (U)
- Input: mask with a 3-pixel component and a 10-pixel component and a hole of 2 pixels inside the big one; min_size=4 and 3.
- Expect: min_size=4 drops the 3-pixel component; hole smaller than min_size is filled (per B13 "fills holes"): b1 == 0 after fill. min_size=3 keeps all (size >= min_size kept; boundary value tested exactly: components of size exactly min_size stay if the pre-Phase-D semantic is `>=`; assert against ORC).
- Catches: off-by-one in size filter, hole filling semantics.

---------------------------------------------------------------------------------------------------------------

## E. Topology - B17

### T50 topology vs cv2 oracle (P, BOTH)
- Input: T13 corpus + T15 families; `min_area` in {0, 1, 2, 3, 7, 100, 10**6}.
- Expect: `instance_topology`, `betti`, `size_filter`, `signature` equal to `topology_oracle` in values AND types (python int / tuple / np.ndarray dtype and shape of size_filter output: same dtype bool/uint8 as oracle; compare `type()` of each element).
- Catches: dtype drift (np.int64 vs int), component order (raster order of first pixel), signature sort (descending), hole counting after filter.

### T51 size_filter semantics (U)
- Input: 5x5 mask: component A size 2, B size 3, hole of size 1 in B; min_area 0,1,2,3,4.
- Expect (derive from ORC, assert literals as well): min_area=0 -> mask unchanged (no hole filling); min_area=2: A stays; =3: A removed; =4: B removed too -> empty.  Hole of size 1 is filled whenever min_area >= 2. Literal output arrays hard-coded.
- Catches: hole filling applied when min_area==0 or 1; strictness of `<` vs `<=`.

### T52 `signatures`/`topologies` batch = loop (U, BOTH)
- Input: 40 masks (incl. empty and 0x0). Expect `signatures(masks, m) == [signature(x, m) for x in masks]`, `topologies(masks) == [instance_topology(x)...]`; exactly one engine call each (spy on `gc.topology` call count with monkeypatch counting wrapper).

### T53 topology want_masks (U)
- `gc.topology(..., want_masks=True)['masks']` equals ORC `size_filter` result for each mask; default call has no 'masks' key.

### T54 signature with holes and zero-area-hole components (U)
- Mask: component with 3 holes; component without hole; expect signature (3, 0) descending; mask with no components: ().

---------------------------------------------------------------------------------------------------------------

## F. certify_polygons - B11

### T60 legal / inter / union / equal (I, BOTH)
- Setup: 300 fuzz polygon shapes (SQ, HOLE from existing tests; 2,853-case fuzz corpus) with W,H in {16, 64, 512}, G in {1000, 1<<20, 1<<30}; references = random masks, the shape's own raster, and empty reference.
- Expect: `legal[i] == (check_polygons defects all zero)`; `inter == popcount(cand & ref)`, `union == popcount(cand | ref)` computed with numpy from the oracle raster; `equal[i] == np.array_equal(cand, ref)`.
- Edge: candidate identical to ref -> inter == union == area, equal True; candidate empty & ref empty -> inter=union=0, equal True (not False due to 0/0); illegal shape: legal False, raster per spec (empty) -> inter 0.
- Catches: popcount of padding bits past W, wrong ref_of mapping when two shapes share a reference, equal computed only on area.

### T61 reference shared by many shapes and out-of-order (I, BOTH)
- ref_of = [2, 0, 0, 1, 2, ...]; expect result independent of permutation of references (T22-like).

### T62 sig / ref_sig vs `signature` (I, BOTH)
- Expect `sig[sig_start[i]:sig_start[i+1]] == signature(raster_i, min_area[i])` exactly, including when min_area removes all components (empty slice) and when raster is empty.

### T63 pair signatures (RIGA rim) (I, BOTH)
- Input: pairs [(0,1)] disc and cup; cup partially outside disc; cup identical to disc (a & ~b empty); cup empty. pair_min_area in {0, 5, 50}.
- Expect: `pair_sig == signature(a & ~b, pm)`, `pair_ref_sig == signature(ref_a & ~ref_b, pm)`; P=0 with `pairs=None` gives empty pair arrays with start == [0] (length 1). Pairs referencing out-of-range shape index raise ValueError.
- Catches: &~ vs ^ ; wrong reference pairing; pair index into wrong batch slot.

### T64 certify_polygons argument validation (U)
- Expect ValueError: W/H 0 or 4097 (polygon sides are 1..4096 here), G 0 or 2**30+1; ref_of out of range; ref_hw not equal to candidate (W,H)? (spec: refs have own hw: assert documented behaviour: mismatch raises); min_area negative; k raters <2 or >32 for tables. Boundary accepted: G=1, G=2**30, W=H=4096 (G).

---------------------------------------------------------------------------------------------------------------

## G. Scorer, validity, eval_sets, metrics - B12, B13, B14

### T70 validity totals regression (I, BOTH)
- Input: 144,133 LIDC polygons fixture subset (e.g. 2,000 sampled with fixed seed) + 2,853 fuzz polygons; `check_polygons` and `score_polygons` 13 defect totals, per-ring counts, ring-pair matrix.
- Expect: equal to arrays recorded before Phase D (golden npz generated by the pre-Phase-D commit, committed in tests/goldens). 0 mismatches. 
- Catches: unintended behavioural change from the refactor.

### T71 scorer bit-identity incl. size filter (I, BOTH)
- Input: score_polygons / score_masks on 200 items with min_area = 0 (filter removes nothing), a value that removes small components, and one that fills holes.
- Expect: float outputs compared by `.view(int64)` / `tobytes` equal to golden; NOT allclose.
- Catches: summation order change; fmad contraction (C1: `-ffp-contract=off`/`--fmad=false`); NaN vs 0.

### T72 metrics() equality (I)
- Input: saved val answers (G9) or 50 synthetic items with K samples (including duplicates, illegal samples, format_ok False).
- Expect: every field equal to golden exactly except `nsd` (|a-b| <= 1e-12); `metrics` never calls `objdel.geometry.to_masks` (monkeypatch it to raise AssertionError); also patch `cv2.connectedComponents*`, `cv2.connectedComponentsWithStats`, `scipy.ndimage.label`, `distance_transform_edt` to raise (G11).
- Catches: hidden CPU fallback; field drift.

### T73 duplicates rule (U)
- Cases (each with K=3, k irrelevant): (a) two samples with equal raw masks, equal flags -> duplicate; (b) equal masks but one `legal=False` -> NOT duplicate; (c) equal masks, different `format_ok` -> not; (d) masks differ by one pixel -> not; (e) both empty raw masks and flags equal -> duplicate; (f) object-wise: two objects (disc, cup) where disc equal and cup differs -> not duplicate; (g) illegal sample with raw mask equal to legal's raw mask but flags differ -> not.
- Expect: use `inter_raw[i,j] == area_raw[i] == area_raw[j]` per object, results as listed; and equals the old "raw sample masks (before the legal gate)" definition from golden.
- Catches: comparing gated (drawn) masks rather than raw (illegal -> empty would make two illegals appear equal AND an illegal vs empty-legal appear equal).

### T74 eval_sets output layout and exactness (I, BOTH)
- Input: n_sets=3, K=4, k=3 (M=8) small masks; compute expected arrays with numpy from rasters: `inter[s,a,b] == popcount(a&b)`, diagonal == area, symmetry `inter == inter.T`, `inter_raw` with illegal sample, `soft`, `surf`, `nsd_out`.
- Hard-coded literal case: W=H=4 with sample = 2x2 square at (1,1), rater identical: surf counts hand computed: 4 corners... [literal from independent numpy surfel implementation in tests/oracles; do not copy the kernel].
- Edge: K=1; k=2 (min), k=32 (max, G); k=1 and k=33 raise ValueError (B16); all samples illegal -> drawn masks empty: area 0, soft P_t = 0, nsd all surfels of rater count (empty other surface rule), surf[sample]=0.
- Catches: M ordering (samples, raters, vote), diagonal, nsd direction swap, empty-surface rule.

### T75 soft cut and nsd_cut host helpers (U, exhaustive)
- soft cut: for n in 1..40 and t in {0, 1e-9, 0.1, 1/3, 0.5, 2/3, 0.7, 0.9, 1.0, 1+1e-12 excluded}: result == min c in [0,n] with float64(c)/n >= t, verified by brute-force independent loop; t=0 -> 0; t=1 -> n; t = c/n exactly (e.g. 1/3 with n=3 -> 1, 0.7 with n=10: float 7/10 >= 0.7 True -> 7; `0.1*3/3` float traps).
- nsd_cut(tau): q largest with sqrt(float64(q)) <= tau: tau=0 -> 0; tau=1 -> 1; tau=1.4142135623730951 -> 2; tau=np.nextafter(sqrt(2), 0) -> 1; tau=2.9999999999999996 -> 8; tau=3 -> 9; tau=1e6 -> 10**12 (or 10**12 - check int overflow), negative tau -> raise or documented (assert ValueError per spec "q >= 0"), tau=nan/inf -> ValueError.
- Catches: float-boundary off by one that would silently shift thresholds (reward changes).

### T76 eval_sets scratch cap invariance (I, BOTH) - B6
- Input: 40-set batch with max_scratch_bytes in {256<<20, 1<<20, enough for exactly one set}.
- Expect: all output arrays byte-identical; cap smaller than one set raises ValueError.

### T77 eval_sets vs metrics from old path (I)
- On 20 real val items: fields computed through the new device path equal the golden metrics (T72); sampled item with all K samples identical -> duplicate share 1.0 (as computed before).

---------------------------------------------------------------------------------------------------------------

## H. EDT and expert tables - B15, B18, B19

### T80 edt bit-equality vs scipy (P, BOTH)
- Input: 200 random masks (sizes 1..120, densities 0.01..0.99), single-zero-pixel masks at corners/centre, single-nonzero-pixel mask, 1xN, Nx1, 512x512 crop, 4096x4096 sparse.
- Expect: `gc.edt(m).tobytes() == scipy.ndimage.distance_transform_edt(m).astype(np.float32).tobytes()`.
- Catches: ties and float32 rounding (sqrt on int64 squared distance then cast: must round like float64 sqrt -> float32, not float32 sqrt); distances > 2^24 squares.

### T81 edt with no zero pixel (U, BOTH)
- Input: all-ones 7x7, 1x1 [[1]], 0xN.
- Expect: every element == `objdel.reward.cache._BIG` (same dtype float32, imported from the module, plus a hard-coded value check taken from the pre-Phase-D source). 0-size input returns an empty (H,W) array, no crash.
- Also all-zero mask -> all zeros.

### T82 `reward.spec.fingerprint()` unchanged (U)
- Expect fingerprint equals the pre-Phase-D literal string recorded in the test.

### T83 add_tables_from_raters == add_tables from host (I, BOTH) - B18
- Setup: 30 caches with k in {2, 3, 4, 32}, H,W in {1x1, 17x33, 128x128, 512x512}, some raters empty, some identical (voter flags), min_area in {0, 1, 5, 50}, t5 values.
- Input: build host inputs with the oracle (tests/oracles host table build), `gc.add_tables`; and `gc.add_tables_from_raters`.
- Expect: returned arrays (perimeter, rater_area, sig, voter, n_vote) equal host-derived ones; then read the stored tables back (test-only dump hook or via scoring) and compare byte-equal: distance map, union, core, vote, etc.; reward terms scored against both handles bit-identical.
- Edge: all raters empty -> n_vote 0, vote plane empty; k=2 where raters disagree (tie votes: vote threshold rule exactly as host); raters with different H,W in one list.
- Catches: vote tie-break, EDT on-device vs scipy, voter flag threshold.

### T84 add_tables_many (U, BOTH)
- Input: 5 argument tuples; expect handles identical in effect to 5 `add_tables` calls; `stats()['waits']` increases by exactly 1 per batch; empty list returns empty int64 array with no wait; `resident()` grows by 5; `drop_tables` each, resident back to start.

### T85 table lifecycle (U)
- drop_tables(invalid handle) raises ValueError, double drop raises; a stale handle used in eval_sets raises (not use-after-free).

### T86 ExpertCache lazy build (U) - B18, C11
- Monkeypatch `gc` entry points with a counter that raises: `ExpertCache(...)` construction and `.load(path)` must make 0 calls; first `.device()` makes tables (counter > 0); a second `.device()` makes none. Fork test: construct cache in parent, `os.fork()` child calls nothing on CUDA (on CPU build simulate with the counter).

### T87 k limits (U) - B16
- k=1 and k=33 raise ValueError from `add_tables_from_raters` (nothing leaked: `resident()` unchanged).

---------------------------------------------------------------------------------------------------------------

## I. Targets and hashes - B9, B10

### T90 target hash regression (E, data-backed)
- Input: for each of the 7 manifests in B9, rebuild every set (`build_target`) with `-I`, recompute `target_set_sha256`.
- Expect: equal to recorded hash (G4). 0 diffs. Run both on CPU build and (once) GPU.
- Catches: ring order / corner order / start corner drift, tie handling in the Phase D polygonizer.

### T91 build_target text equality on synthetic corpus (I)
- 300 fixtures (fuzz masks, with holes, multi-component, min_area bigger than components, RIGA disc/cup/rim): `.text` equals text from archived pre-Phase-D implementation (golden file of sha256 per case).

### T92 build_target min_area required (U) - B10
- Input: min_area dict missing "lidc"/name -> ValueError; missing "rim" for RIGA -> ValueError; extra keys ignored or raise per old behaviour (assert old); `min_area=0` accepted; negative raises.

---------------------------------------------------------------------------------------------------------------

## J. Engine counters, errors, no-CPU geometry - B16, B19, G11, G15

### T100 waits == 3 and zero allocation after warm-up (I, BOTH CPU+GPU)
- Setup: `gc.reset_stats()`, one `score_polygons` warm-up; reset; 10 calls same shape.
- Expect: per call waits == 3 (30 total); device_allocs == device_releases == 0 after warm-up; host_allocs, host_releases per spec 0 (arena stable); `drop_tables` makes 0 waits. A call of a larger shape allocates (counter increases) - proves counters live.
- Also `timing=True` call: dict has `timing`; default call has no `timing` key, and waits unchanged by default.
- Catches: hidden cudaMemcpy sync, per-call cudaMalloc.

### T101 stats after errors (U)
- A failing call (ValueError) leaves waits/allocs counters sane (non-negative), and the next valid call works and `arena_bytes` does not grow unboundedly across 100 failing calls.

### T102 overflow guard (U)
- Input: `start` claims > INT_MAX items (use numpy views with `start` values only, no allocation: bits small but n*... use lazy array lengths via a zero-stride broadcast or crafted `hw`: H=W=4096 and n=129 (129*16.7M = 2.16e9 > INT_MAX) with `bits` = np.broadcast_to / np.lib.stride_tricks (no memory) when binding accepts; else mark G).
- Expect: ValueError/std::invalid_argument, never wrap negative.

### T103 CUDA failure message (G)
- Input: provoke cudaErrorMemoryAllocation by requesting scratch beyond the free memory via huge valid batch (e.g. 4096x4096 x 600 masks, cap raised) or fault injection test hook.
- Expect: `RuntimeError` whose text contains the CUDA error string ("out of memory"); engine usable afterwards.

### T104 no device (U, CUDA build, no GPU) 
- Run the nvcc build with `CUDA_VISIBLE_DEVICES=""`: importing is fine; first engine call raises clean RuntimeError (no segfault, process exit code 0 after catching). CPU build: `device_count() == 0` and everything works.

### T105 no CPU geometry (E) - G11
- Patch `cv2.connectedComponents`, `cv2.connectedComponentsWithStats`, `scipy.ndimage.label`, `scipy.ndimage.distance_transform_edt` to raise; run `build_target` on LIDC and RIGA fixtures, `ExpertCache.build/load/device`, `certify`, `metrics()`.
- Expect: no error. Static assertions: `grep -R "tests.oracles\|tests/oracles" src/objdel` empty; files `src/objdel/checker/polygonize.py`, `bitquad.py` do not exist; `from objdel.checker import bitquads` raises ImportError; `trace, polygonize, Ring, MaskRings, certify, certify_many, as_mask, betti, instance_topology, min_area, signature, size_filter, signatures, topologies` all importable.

### T106 Ring / MaskRings value semantics (U)
- `Ring` frozen: assigning raises FrozenInstanceError; `hole` property == (area < 0); pts dtype int64 shape (c,2); `MaskRings.quads` is a 5-tuple of python ints; b0/b1 python ints.

---------------------------------------------------------------------------------------------------------------

## K. GPU-only equality, sanitizers, performance

### T110 module is the CUDA build (G)
- `gc.build() == "cuda"` on a GPU host (fail, not skip). `gc.device_count() == nvidia-smi -L` count.

### T111 CPU vs GPU bit-equality (G) - G1, G3
- Input: corpus of 6,024 fuzz masks, 2,853 fuzz polygons, 37,093 crops of `data/checker_eval/masks.npz`, 144,133 LIDC polygons.
- Expect: GPU arrays byte-equal to CPU-build arrays saved from the Mac (npz golden), for `check_masks`, `topology`, `certify_polygons`, `score_polygons`, `eval_sets`, `add_tables_from_raters`; 0 mismatches; 37,093 crops give 38,002 rings and 1,476,508 corners (hard-coded regression numbers).
- Catches: any float contraction, atomics order, uninitialised read.

### T112 sanitizer corpus (S) - G5
- `tests/test_gpu_sanitize.py` small fixed corpus: T01, T02, T05, T09 (borders), 4096x1, 1x4096, spiral 64x64, 20 fuzz masks, certify_polygons 8 shapes with RIGA pair, eval_sets 2 sets, add_tables_from_raters 2 caches (k=2 and k=32).
- Expect: `ERROR SUMMARY: 0 errors` under memcheck, initcheck, synccheck, racecheck, each < 300 s.
- Catches: OOB reads at borders (T09), uninit scratch (T25), missing __syncthreads in pointer-jumping.

### T113 make gpu-check DRY=1 (U, Mac)
- Expect: exit 0, prints the chain including memcheck/initcheck/synccheck/racecheck with `timeout 300`; runs `pytest --collect-only -q -m gpu` and collected count > 0; if any named test file is renamed/missing, exit non-zero (rename `tests/test_gpu_sanitize.py` temporarily in a tmp copy of the repo). Marker `gpu` and `sanitizer` are registered (no PytestUnknownMarkWarning with `-W error`).

### T114 polygonizer speed gate (F, G) - G7
- Workloads: 256 frames 512x512, 37,093 crops, 1 spiral 4096x4096?, 40 small masks. Alternate old/new 7 runs each after warm-up; median.
- Expect: outputs byte-identical; new median <= 1.05 x old on all; <= 0.9 x old on the 256-frame and crops workloads. Use CUDA events, not wall-clock only; report spread.

### T115 reward-path speed gate (F, G) - G8
- 256 completions replay: new `score_polygons` median <= 1.05 x old, waits == 3.

### T116 audit scratch bound (G) - G14
- 6 RL steps with audit callback: `stats()['arena_bytes'] <= eval_sets cap + table bytes` sampled each step; no OOM.

### T117 post-warm-up allocation in training (G) - G15
- steps 2-6: device_allocs delta == 0 and device_releases delta == 0.

---------------------------------------------------------------------------------------------------------------

## Pytest sketches for the 5 most important tests

Most important: (1) oracle + invariants fuzz, (2) batch/chunk invariance, (3) certify independence mutation, (4) edt bit-equality, (5) tables device == host. (Hash regression T90 and CPU-GPU equality T111 are covered by existing gates and are data-driven; they are next in priority.)

```python
# tests/test_phase_d_core.py
import numpy as np, pytest, cv2
from objdel import _gpucheck as gc
from objdel.checker import trace, polygonize
from tests.oracles import polygonize_oracle as PO, bitquad_oracle as BQ

def invariants(m, mr):
    m = (np.asarray(m) != 0).astype(np.uint8)
    n, lab = cv2.connectedComponents(m, connectivity=8)
    assert sum(r.area for r in mr.rings) == int(m.sum())
    assert sum(not r.hole for r in mr.rings) == mr.b0 == n - 1
    assert sum(r.hole for r in mr.rings) == mr.b1
    for r in mr.rings:
        d = np.diff(np.vstack([r.pts, r.pts[:1]]), axis=0)
        assert ((d[:, 0] == 0) ^ (d[:, 1] == 0)).all()          # axis-aligned, no zero edge
        x, y = r.pts[0]                                          # first corner = start of first east edge
        assert lab[y, x] == r.part + 1                           # pixel right of the first (east) edge is below it
    assert sorted(set(r.part for r in mr.rings)) == list(range(mr.b0))

def fuzz_masks(seed=0):
    rng = np.random.default_rng(seed)
    out = [np.zeros((0, 0), np.uint8), np.zeros((3, 3), np.uint8), np.ones((1, 1), np.uint8),
           np.eye(2, dtype=np.uint8), np.ones((1, 33)), np.ones((33, 1))]
    out += [np.indices((n, n)).sum(0) % 2 for n in (16, 17, 31, 32, 33)]           # checkerboards
    for p in (0.05, 0.3, 0.5, 0.7, 0.95):
        for n in (31, 32, 33, 40):
            out.append((rng.random((n, rng.integers(1, 70))) < p).astype(np.uint8))
    return out

def test_trace_matches_oracle_and_invariants():
    masks = fuzz_masks()
    got = trace(masks)
    assert len(got) == len(masks)
    for m, mr in zip(masks, got):
        want = PO.polygonize(np.asarray(m) != 0)
        assert [r.area for r in mr.rings] == [r.area for r in want]
        for a, b in zip(mr.rings, want):
            assert np.array_equal(a.pts, b.pts)                                    # same corners, same start
        assert mr.quads == tuple(BQ.bitquads(np.asarray(m) != 0)[:5]) or True      # replace with real oracle accessor
        if m.size:
            invariants(m, mr)

@pytest.mark.parametrize("cp,sb", [(4096, 65536), (1 << 25, 65536), (4096, 1 << 30)])
def test_split_and_order_invariance(cp, sb):
    masks = fuzz_masks(1)
    def pack(ms):
        from experiments.cuda_checker import compare as C
        return C.pack_masks(ms)
    p = pack(masks)
    ref = gc.check_masks(*p, )                                            # defaults
    out = gc.check_masks(*p, chunk_pixels=cp, max_scratch_bytes=sb)
    for k in ref:
        assert ref[k].tobytes() == out[k].tobytes(), k
    # alone vs batch
    for i in (0, 5, 17, len(masks) - 1):
        alone = trace([masks[i]])[0]
        inbatch = trace(masks)[i]
        assert [r.pts.tobytes() for r in alone.rings] == [r.pts.tobytes() for r in inbatch.rings]
        assert [(r.area, r.part) for r in alone.rings] == [(r.area, r.part) for r in inbatch.rings]

def test_certify_independence_mutation(mutant_build):
    """mutant_build: fixture that loads the test-only build whose labeller joins foreground by 4-adjacency."""
    from objdel.checker import certify
    diag = np.eye(2, dtype=np.uint8)
    assert certify(diag).certified is True
    with mutant_build():
        r = certify(diag)
        assert r.certified is False and any(f.startswith("rings") for f in r.failures)
        assert (r.b0, r.b1) == (1, 0)          # b0 from the pixel label chain, not the mutant rings

@pytest.mark.parametrize("seed", range(5))
def test_edt_bit_equal_scipy(seed):
    from scipy.ndimage import distance_transform_edt
    rng = np.random.default_rng(seed)
    for _ in range(40):
        h, w = rng.integers(1, 90, 2)
        m = (rng.random((h, w)) < rng.uniform(0.02, 0.98)).astype(np.uint8)
        if m.all():
            continue
        assert gc.edt(m).tobytes() == distance_transform_edt(m).astype(np.float32).tobytes()
    from objdel.reward import cache
    allone = np.ones((7, 5), np.uint8)
    assert (gc.edt(allone) == cache._BIG).all() and gc.edt(allone).dtype == np.float32

def test_device_tables_equal_host_tables(host_tables):
    """host_tables: helper in tests/oracles building the add_tables inputs (dist, union, core, vote, areas, flags,
    sigs, n_vote, perimeter) from raters on the host; gc.dump_table(handle) is a test-only readback."""
    rng = np.random.default_rng(3)
    for k in (2, 3, 32):
        raters = (rng.random((k, 40, 51)) < 0.4).astype(np.uint8)
        t5, ma = np.array([0.5]), np.array([3], np.int64)
        dev = gc.add_tables_from_raters([raters], t5, ma)
        host = gc.add_tables(*host_tables(raters, 0.5, 3))
        assert dev["n_vote"][0] == host_tables.n_vote
        assert gc.dump_table(dev["handles"][0]) == gc.dump_table(host)            # byte-equal planes
        gc.drop_tables(dev["handles"][0]); gc.drop_tables(host)

def test_waits_and_allocations(score_args):
    gc.score_polygons(*score_args)           # warm-up
    gc.reset_stats()
    for _ in range(10):
        s0 = gc.stats()["waits"]
        gc.score_polygons(*score_args)
        assert gc.stats()["waits"] - s0 == 3
    st = gc.stats()
    assert st["device_allocs"] == 0 and st["device_releases"] == 0
```

Notes on the sketches: `gc.dump_table`, `host_tables`, and `mutant_build` are test-support hooks the spec does not define; D.spec.md must either provide them or the tests must compare via reward terms scored against both handles (B18 second clause). The quad assertion line has an `or True` placeholder to be replaced by the real oracle accessor once `bitquad_oracle` is read; do not ship with it.

## Gaps and ambiguities found in the spec (for the spec owner)
1. B3 says `part+1` equals the label at the pixel "right of the ring's first edge": the first edge is east, so right of an east edge (foreground right-hand side, y down) is the pixel BELOW it, (x, y). Tests use pixel (x, y) of the first corner; confirm.
2. B16 "A scan over more than INT_MAX items raises": no input size is given that a test can build cheaply; T102 needs a zero-copy trick or a test hook.
3. `check_masks` behaviour when one mask alone exceeds `max_scratch_bytes` is unspecified (T21 assumes ValueError).
4. Whether size filter keeps `area >= min_area` or `> min_area` is not stated; tests defer to the oracle (T44, T51).
5. B18 needs a way to read stored tables (dump hook) or a reward-term proxy; unspecified.
6. T05: "outer rings count == b0" with a self-touching ring is assumed; spec B4 supports it, but explicit sample outputs would help.
7. No limit given on batch count n; none tested.
