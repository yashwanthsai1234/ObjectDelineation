# Phase D test cases, set 2 (written from D.spec.md only; adversarial, property, performance, failure modes)

## Skills used
- nvidia-skills:nvidia-skill-finder: CUDA/GPU in scope, but no catalog install (needs owner approval). Used only to keep sanitizer, nvcc and RTX 5090 gates explicit.
- engineering:testing-strategy: pyramid split. Many unit/property tests on the CPU build (run on the Mac in `make test`), few `gpu`-marked integration tests, gates G1-G15 mapped to cases.
- mattpocock-skills:tdd: tests only at public seams (`gc.*`, `objdel.checker.*`, `metrics`, `build_target`, `ExpertCache`). Expected values come from independent oracles (`tests/oracles/`, cv2, scipy) or literals, never from recomputing the engine's own algorithm.
- engineering:code-review / debug / system-design / architecture / documentation / huggingface-skills:huggingface-best: not needed for test authoring; consulted only for failure-mode thinking (stale state, ordering, limits). No changes to REPO.

## Conventions
- "oracle" = `tests/oracles/polygonize_oracle.polygonize`, `bitquad_oracle.bitquads`, `topology_oracle.*` (cv2), `scipy.ndimage`. Tests never import oracles from `src/objdel` (C5).
- Marks: no mark = runs on the CPU build, Mac, in `make test`. `gpu` = RTX 5090 + nvcc build only. `slow` = needs data/`-m slow`. `sanitizer` = fixed small corpus.
- `gc` = `objdel._gpucheck`. `rings` = `objdel.checker.rings`. Seeds are fixed (numpy `default_rng(seed)`), so failures reproduce.
- Helper `canon(rings)` = `[(r.pts.tolist(), r.area, r.part) for r in rings]`. Helper `pack(masks)` = (bits uint8[N], hw int32[n,2], start int64[n+1]) with row-major flattening.
- Bit-equality means `np.array_equal` / `tobytes()` equality, never `allclose`, except `nsd` (1e-12 abs, B14).

---------------------------------------------------------------------------------------------------------------

## A. Polygonizer / rings (B1-B7)

### 1. trace_equals_oracle_on_fuzz_families
- Level: property. Marks: none (CPU), also `gpu`.
- Setup: families from G1: noise at density 0.05,0.1,...,0.95 (sides 1..70 random), blobs with holes, checkerboard (odd/even phase), nested squares (depth 1..8), comb, spiral, diagonals (both slopes), stripes (horizontal/vertical), 1xN, Nx1, sides 31/32/33 (word boundaries), empty. 25 seeds per family.
- Input: `trace(masks)` on the whole list in one call.
- Expected: for every mask: `canon(trace_i.rings) == canon(oracle.polygonize(mask))`; `trace_i.quads == tuple(bitquads(mask))`. Exact, including ring order, start corner, signed area.
- Why: core exactness (G1). Catches wrong start corner, wrong ring order, wrong hole sign, off-by-one in word-boundary handling (31/32/33).

### 2. single_pixel_every_position
- Level: unit. Marks: none.
- Setup: 5x5 mask with one set pixel at each of 25 positions; also 1x1 set.
- Expected: exactly one ring, 4 corners, area +1, part 0; quads == (n1=1,0,0,0,0)-style counts equal to oracle; `betti == (1,0)`. Corner list starts at the minimum edge key: for pixel (x,y) the first corner is (x,y) with direction east.
- Why: smallest legal shape; touches every border combination; catches border-as-background bugs.

### 3. empty_and_degenerate_shapes
- Level: unit. Marks: none.
- Input: masks of shape (0,0), (0,7), (7,0), (1,0), (1,1) all zero, (3,3) all zero, (1,5) zeros, (5,1) zeros, plus one non-empty in the same batch.
- Expected: no rings, `quads == (0,0,0,0,0)`, `betti == (0,0)` for each empty one, in the same call as the non-empty one that still equals oracle. `ring_start` diff == 0 for empty. Nothing raises.
- Why: B5. Zero-size masks produce zero-length slices; common source of div-by-zero, `-1` index, or null pointer passed to the kernel.

### 4. full_mask_and_stripe_sizes
- Level: unit. Marks: none.
- Input: all-ones masks of shape (1,1),(1,64),(64,1),(31,31),(32,32),(33,33),(64,64),(65,3). Odd-width stripes.
- Expected: one outer ring with 4 corners, area H*W, part 0, b=(1,0), equal oracle.
- Why: collinear corner collapse at all borders; word-length padding must not be read as foreground (padding bits).

### 5. holes_and_nested_components
- Level: unit. Marks: none.
- Input: (a) 7x7 ring (frame) with a 1-pixel hole; (b) frame inside frame inside frame (3 outer, 3 holes... check count) on 21x21; (c) hole touching the image border only diagonally (must NOT be a hole, border counts as background); (d) background region connected only diagonally (4-conn background splits it so two holes); (e) foreground diagonal pair enclosing a diagonal-only gap.
- Expected: equals oracle; for (a) rings = [outer +48, hole -1]; `parts`: hole.part == outer.part; (c) b1 == 0; (d) b1 == 2; (e) b0 == 1 under 8-connectivity.
- Why: B2 connectivity pair (8 fg / 4 bg). A swapped pair is the classic polygonizer bug.

### 6. part_labels_vs_cv2
- Level: property. Marks: none.
- Input: 300 random masks (density 0.2-0.7, sides 1..64) plus nested squares and spiral.
- Expected: for each ring, with first edge key giving (x,y): `ring.part + 1 == cv2.connectedComponents(mask.astype(uint8), connectivity=8)[1][y, x]` (foreground pixel right of the first east edge is (x, y) = first corner). Outer ring: components count of distinct parts among outer rings == b0; each part value appears exactly once among outer rings; each hole's part is the part of the outer ring that contains it (verify by flood: part of the nearest enclosing foreground neighbour of the hole's first pixel).
- Why: B3/G2. Independent oracle (cv2) so the test is not tautological.

### 7. betti_euler_identities
- Level: property. Marks: none.
- Input: same 300 masks as test 6.
- Expected: `b0 - b1 == chi8` computed from quads as `(n1 - n3 - 2*nD) / 4`; number of rings with area>0 == b0; with area<0 == b1; `sum(area) == mask.sum()`; `betti == (cv2 8-comp count, scipy.label bg 4-comp count with a padded border minus 1)`.
- Why: B4. Independent checks (chi from quad counts, components from cv2/scipy). Catches a labeller that merges across diagonals only in one direction.

### 8. offsets_add_to_every_corner
- Level: unit. Marks: none.
- Input: 20 masks, offsets from {(0,0),(5,-3),(-7,100),(2**20,2**20),(-1000,-1000)}; `trace(masks, offsets)` and `polygonize(mask, offset)`.
- Expected: corners == oracle corners + (ox, oy) exactly; `area`, `part`, ring order unchanged; `polygonize(m,o) == trace([m],[o])[0].rings`; `pts.dtype == int64`. A negative offset yields negative coordinates (no wrap to uint).
- Why: B7. Offsets are applied host side; int32 engine output plus large offset overflows if added in int32.

### 9. batch_composition_independence
- Level: property. Marks: none, plus `gpu` repeat with 20 repeats.
- Input: 40 mixed masks (empty, 1x1, 1xN, noise 512x512, spiral 1024x1024, nested). Run (i) each alone, (ii) all together, (iii) reversed, (iv) 5 random permutations, (v) with `chunk_pixels` in {1, 64, 4096, 65536, 1<<25} and `max_scratch_bytes` in {4096, 65536, 1<<30} (grid product).
- Expected: output arrays of each mask byte-identical in every configuration (compare the dict slices, after un-permuting). `chunk_pixels=1` and tiny scratch must still work (a chunk holds at least one mask) or raise ValueError with a clear message, never return wrong data. The same for `eval_sets` under `max_scratch_bytes` in {1<<12, 1<<16, 1<<28}.
- Why: B6. Catches cross-mask contamination via shared scratch (stale labels from a previous mask), chunk-boundary off-by-one, a mask larger than the scratch budget.

### 10. stale_scratch_poisoning
- Level: adversarial. Marks: none, `gpu`.
- Setup: run a large dense noise batch (many rings) to dirty scratch; then immediately run a small batch of sparse masks, then an empty-only batch; repeat in the opposite order.
- Expected: each result equals a fresh-process result (compare against the oracle and against the same batch run first). With a pool allocator that does not zero memory, any kernel that reads uninitialised scratch fails here.
- Why: stream-ordered pool reuse (C4) means uninitialised reads give non-deterministic results; initcheck also covers this.

### 11. maximum_size_and_limits
- Level: integration. Marks: `slow`, `gpu` for 4096x4096 (CPU build allowed but slow).
- Input: 4096x4096 noise 0.5 and a 4096x4096 checkerboard (worst case about 8.4 million rings of 4 corners... checkerboard diagonal touch is one 8-component each; count by oracle on a 512 sub-sample), 4096x1, 1x4096, then 4097x1, 1x4097, 4097x4097.
- Expected: sizes <= 4096 equal the oracle (use the sub-tiled oracle where the oracle is too slow, comparing ring count, corner total and area sum, plus the full oracle for 4096x1); sizes 4097 raise `ValueError` (not crash, not truncate); the error text names the limit 4096.
- Why: B16 and the largest possible scratch. Catches int32 overflow of keys `(y*(W+1)+x)*4+d` (4097*4097*4 = 67 million fits; guard against future change) and chunk logic when one mask exceeds the chunk.

### 12. dtype_and_layout_robustness
- Level: unit. Marks: none.
- Input: masks as bool, uint8 {0,1}, uint8 {0,255}, int64 {0,-3,7}, float32 {0.0,0.5}; non-contiguous (transposed, `[::2, ::3]` slices, Fortran order); read-only array; 3-D input (1,5,5); list-of-lists.
- Expected: any nonzero = foreground, so output identical to `mask != 0` contiguous bool; a 3-D or 0-D input raises `ValueError`; input arrays are not modified (`array_equal` before/after, `flags.writeable` preserved).
- Why: B1 says any nonzero. Catches stride bugs (`.tobytes` vs `ascontiguousarray`) and mutation of caller data.

### 13. bad_input_never_crashes
- Level: failure mode. Marks: none.
- Input to `gc.check_masks`: `hw` with negative entries, `hw` product != `start` diff, `start` not starting at 0, not monotone, `start[-1] != len(bits)`, `bits` values 2..255 (valid, nonzero), wrong dtype (int64 bits, float hw), `n=0` (empty arrays), mismatched lengths, `hw` columns swapped.
- Expected: `ValueError` (never segfault; run each case in a subprocess to catch an abort and assert return code 0 and the exception type). `n=0` returns empty arrays with correct shapes (`ring_start == [0]`, `pts_start == [0]`).
- Why: B16. Raw pointer/array bindings are the likely site for memory faults.

---------------------------------------------------------------------------------------------------------------

## B. check / certify (B8)

### 14. certify_independent_of_ring_labeller (mutation test)
- Level: adversarial. Marks: none (needs the test-only build flag; skip if the flag is absent, fail if present but certify returns True).
- Setup: build with the test-only "4-adjacency ring labeller" option. Mask 2x2 with pixels (0,0) and (1,1).
- Expected: normal build: `certify(mask).certified is True`. Mutated build: `certified is False`. Also with a mutated build, `check_masks` betti shows b0 == 2 while `gc.topology` shows b0 == 1; `certify` must expose this disagreement.
- Why: G16. Without independence the certificate checks the labeller against itself.

### 15. certify_fields_and_failure_cases
- Level: unit. Marks: none.
- Input: 100 fuzz masks, `min_size` in {0,1,2,5,50}. Also masks with: diagonal-only touching; hole fully removed by `min_size` (hole smaller than min_size); component removed by size filter; empty mask.
- Expected: `certified` equals `(quad_chi8 == b0-b1) and (outer,holes)==(b0,b1) and area_sum==quad_area==mask.sum() and corners==quad_corners` computed independently from the oracles (cv2 for b0, b1). Report field names and types equal the pre-Phase-D `Report` (snapshot of `dataclasses.fields`). `certify_many(masks) == [certify(m) for m in masks]`. Empty mask certifies True.
- Why: B8. The truth table is computed outside the engine; a swapped condition is caught.

---------------------------------------------------------------------------------------------------------------

## C. Topology and signatures (B17)

### 16. topology_matches_cv2_oracle_all_min_area
- Level: property. Marks: none.
- Input: 500 masks (fuzz + nested + holes + diagonal), `min_area` in {0,1,2,3,10,100,10**9} per mask (vectorised by one list).
- Expected: `instance_topology`, `betti`, `size_filter`, `signature` equal the oracle in value AND type (`type(x)`, element types: Python int vs numpy int; `size_filter` returns a bool/uint8 mask of the same dtype/shape as oracle). Components in raster order of their first pixel; signature sorted descending. `min_area=0` keeps everything; `min_area` larger than every component gives empty topology and an all-zero mask.
- Why: B17 including type equality. Callers use these as dict keys / hashes (target_set_sha256).

### 17. size_filter_hole_fill_boundary
- Level: edge. Marks: none.
- Input: a 9x9 square with a hole of area exactly `m`, `m-1`, `m+1` pixels for `min_area = m`; a component of area exactly `min_area`; a 1-pixel component touching the border; two components whose hole totals differ.
- Expected: "holes smaller than min_area are filled; area == min_area is kept" exactly as the cv2 oracle says (comparison direction). `comp_holes` counts only surviving holes. Boundary values `min_area` in {area-1, area, area+1} produce all three oracle outcomes.
- Why: off-by-one `<` vs `<=` is the dominant bug in size filters, and it changes target texts (B9).

### 18. signatures_and_topologies_single_engine_call
- Level: integration. Marks: none.
- Setup: `gc.reset_stats()`; call `signatures(masks, 3)` with 200 masks; `topologies(masks)`.
- Expected: results equal per-mask `signature`/`instance_topology`; each batch function triggers exactly one engine call (count by wrapping `gc.topology` with a counter spy or by `gc.stats()` waits == 1 for the CPU-count of waits). `signatures([]) == []`.
- Why: spec 3.1 states one engine call each. A Python loop fallback would silently pass value tests.

---------------------------------------------------------------------------------------------------------------

## D. Targets and hashes (B9, B10)

### 19. recorded_target_hashes_reproduce
- Level: regression. Marks: `data` (skipped where files are missing).
- Setup: for each manifest listed in B9 (7 files), rebuild each set's targets with the public pipeline and hash.
- Expected: `target_set_sha256` byte-identical to the recorded value. Zero tolerance. Also a smaller always-run copy: `tests/goldens/` text of 50 fixed LIDC and RIGA fixture items equals the stored pre-Phase-D text.
- Why: G4. Any ring-order or start-corner drift changes the text the model is trained on.

### 20. build_target_min_area_required
- Level: failure mode. Marks: none.
- Input: `build_target` with `min_area` missing a name; for RIGA without "rim"; with an extra unknown name; with `min_area={}`; with negative value; value of wrong type (str).
- Expected: missing name or missing "rim" -> `ValueError` naming the missing key; extra key: same behaviour as before Phase D (record in the golden); negative/str -> `ValueError`/`TypeError` as before; no partial output.
- Why: B10.

---------------------------------------------------------------------------------------------------------------

## E. certify_polygons / validity / scorer (B11-B13)

### 21. certify_polygons_against_set_algebra
- Level: property. Marks: none, `gpu`.
- Input: 2,853-style fuzz polygons (use 300), reference masks random and nested; G random in {1, 7, 1000, 2**30}; `min_area` in {0, 1, 4, 20}; pairs for RIGA rim `(a,b)` including a==b, and shapes with no pair.
- Expected: `legal == (check_polygons(...) has no defect)`; `inter == popcount(cand & ref)`, `union == popcount(cand | ref)`, `equal == array_equal`, `inter + union` identity: `union == area_c + area_r - inter` (independent check); `sig == signature(cand, min_area)`; `pair_sig == signature(raster_a & ~raster_b, pair_min_area)`. The candidate raster is produced by an independent rasteriser (`score_numpy.py`/raster oracle).
- Why: B11. Catches wrong sentinel for illegal shapes (illegal shape must still return legal False, inter/union defined).

### 22. certify_polygons_identical_and_disjoint
- Level: edge. Marks: none.
- Input: candidate == reference; candidate empty; reference empty; both empty; candidate disjoint; candidate superset.
- Expected: equal=True and inter==union==area only in the first; both empty -> `equal True`, inter 0, union 0 (no divide); superset -> inter == area_ref.
- Why: empty sets are where `equal` implemented as `inter==union` can be wrong (spec says equal exactly when masks equal; both empty is equal).

### 23. validity_regression_all_13_defect_classes
- Level: regression. Marks: none.
- Input: the existing 13-defect corpus (`test_every_defect_class_is_exercised`) and `tests/goldens/polycore.npz` (if present).
- Expected: per-class totals, per-ring counts and ring-pair matrix are identical to the golden recorded before Phase D, for `check_polygons` and `score_polygons`; every class has count >= 1 in the corpus (guard that the corpus still exercises each class).
- Why: B12.

### 24. scorer_bit_identity_with_and_without_size_filter
- Level: regression. Marks: none, `gpu`.
- Input: 200 fixed items where (a) the filter removes nothing, (b) removes small components, (c) fills a hole, (d) removes everything; scored with `score_polygons` and `score_masks`.
- Expected: output arrays `tobytes()` equal the stored golden floats from the pre-Phase-D build; CPU build equals GPU build bit for bit.
- Why: B13/C2. Float results from integer inputs by one thread: a reordered sum shows up as a 1-ulp change.

### 25. illegal_and_hostile_polygon_inputs
- Level: failure mode. Marks: none.
- Input: >10 rings, >1000 vertices in a ring, W or H in {0, 4097}, G in {0, 2**30+1}, NaN-free but out-of-range coordinates (negative, >G), `ring_start` not monotone, `roles` of wrong length or value 99, `k` rater count 1 and 33 for tables.
- Expected: `ValueError` for limit violations (B16); legal-but-bad geometry (out-of-range coordinates) returns a defect result as before Phase D, not an exception (compare to the pre-Phase-D golden behaviour: "every input that the parser and binding accept"). Exactly at the limit (10 rings, 1000 vertices, W=H=4096, G=2**30, k=2 and k=32) succeeds.
- Why: B16/B12. Boundaries on both sides of every limit.

---------------------------------------------------------------------------------------------------------------

## F. eval_sets and metrics (B14)

### 26. eval_sets_inter_matrix_vs_numpy
- Level: property. Marks: none, `gpu`.
- Input: 6 sets, K in {1, 8, 16}, k in {2, 4, 32}, random drawn masks, 64x64 .. 512x512, one set with all empty, one with all identical.
- Expected: `inter[s,a,b] == (M_a & M_b).sum()` computed with numpy; matrix symmetric; diagonal == area; `inter_raw` / `area_raw` equal numpy on RAW sample masks; illegal samples (legal=False) have empty drawn rows/cols in `inter` but their raw values keep their pixels. `surf` equals the independent 2x2-code surfel counts (numpy, classes sqrt(2)/2,1,sqrt(2)); empty mask -> zeros.
- Why: B14 inputs. The raw-versus-gated split is the subtle point; a bug there changes duplicate share.

### 27. soft_cut_and_nsd_cut_exactness
- Level: unit. Marks: none (host helpers).
- Input: `soft_cut(n, t)` for n in 1..40 and t in {0, 0.1, 0.3, 0.5, 1/3, 2/3, 0.7, 0.9, 1.0, 1e-12}; `nsd_cut(tau)` for tau in {0, 0.5, 1, 1.4142135623730951, 1.5, 2, 3.0000000000000004, 2.9999999999999996, 100.7}.
- Expected: soft_cut = smallest c in [0,n] with `float64(c)/n >= t` verified by brute force loop in the test (independent, not by formula); never above n; `nsd_cut` = largest integer q>=0 with `sqrt(float64(q)) <= tau`, brute force up to q = 1e5; tau = 0 -> 0.
- Why: spec 3.2 host helper contracts; float rounding at 0.3*10 and sqrt boundaries is where rounding errors bite.

### 28. duplicate_definition
- Level: unit. Marks: none.
- Input: items of K=4 samples: (a) identical raw masks and flags; (b) identical raw but different `legal`; (c) different `format_ok`; (d) raw masks subset of each other (inter == area_i but area_j larger); (e) both empty raw masks (equal); (f) multi-object: one object equal, another differs.
- Expected: duplicates exactly when flags equal and for every object `inter_raw[i,j] == area_raw[i] == area_raw[j]`; (d) is NOT a duplicate; (e) is a duplicate; (f) not a duplicate.
- Why: B14 duplicate share is a headline metric.

### 29. metrics_regression_and_no_to_masks
- Level: regression. Marks: none (fixture) and `data` (saved val answers, G9).
- Setup: monkeypatch `objdel.geometry.to_masks` to raise; cv2 connected-component functions and scipy label/edt to raise (G11).
- Expected: `metrics()` returns every field of the pre-Phase-D golden; all fields exactly equal; `nsd` within 1e-12; the set of field names equal. No call to `to_masks`.
- Why: B14/G9/G11.

---------------------------------------------------------------------------------------------------------------

## G. EDT and expert tables (B15, B18)

### 30. edt_bit_identity_property
- Level: property. Marks: none, `gpu`.
- Input: 300 random binary maps (sides 1..300, density 0.001..0.999), single zero pixel, single nonzero pixel, a zero pixel at each corner, thin lines, a full-zero map, an all-nonzero map, 1xN, Nx1, a very elongated 4096x2 map.
- Expected: `edt(nz).tobytes() == scipy.ndimage.distance_transform_edt(nz).astype(float32).tobytes()` (bytes, so 0.0 vs -0.0 and NaN would be caught); all-nonzero returns `_BIG` everywhere (import the existing constant from `objdel.reward.cache`, assert its value is unchanged by recording a literal); `reward.spec.fingerprint()` equals the stored fingerprint literal.
- Why: B15. EDT is the classic place for 1-ulp differences (sqrt of integer squared distance; the distance in float64 then cast to float32).

### 31. edt_ties_and_diagonal_distances
- Level: adversarial. Marks: none.
- Input: maps where many nearest zeros are equidistant (a circle of zeros around a centre, zeros at (±3,±4) -> distance 5, (±5,0),(0,±5)); large distances: a 4096x4096 map with one zero at a corner (squared distance up to 2*4095^2 = 33.5 million, exceeding float32 24-bit mantissa for exact integers).
- Expected: bytes equal to scipy; the largest value equals `float32(sqrt(2)*4095)` as scipy computes it.
- Why: float32 exactness above 2^24; an int32/float32 squared-distance accumulator would differ.

### 32. device_tables_equal_host_tables
- Level: integration. Marks: none, `gpu`; `data` variant for 1,000 real caches (G10).
- Input: m=20 caches with k in {2,3,4,32} raters each, mixed sizes, empty raters (an all-zero rater plane), a rater plane all ones, `t5` values {0, .25, .5, .75, 1}, `min_area` in {0, 3, 100}. Build via `gc.add_tables_from_raters` and via `gc.add_tables` with host-built inputs from the oracle host build in `tests/oracles/`.
- Expected: for each handle the stored tables (distance map, union, core, vote, rater areas, voter flags, signatures, n_vote, perimeter) are byte-equal (read back via a test-only dump function or via the reward terms); reward terms for 50 candidate masks scored against both equal bit for bit. `add_tables_many` equals one-by-one `add_tables` and makes exactly 1 wait per batch (B19 counters).
- Why: B18. Verifies the move of table construction to the device by an independent host build.

### 33. expert_cache_makes_no_engine_call
- Level: failure mode. Marks: none.
- Setup: spy on every `gc` attribute (wrap with a counter); also `os.fork()` after `ExpertCache.build`/`load`.
- Expected: construct and `load` -> zero calls to `gc` functions that touch the device (counts 0, also `gc.stats()` unchanged); first `device()` creates tables (>=1 call); a second `device()` makes none; child process after fork can use `device()` without hanging or CUDA error (GPU variant: child process created before first `device()`).
- Why: B18/C11. CUDA contexts are not fork-safe.

### 34. table_handle_lifecycle
- Level: failure mode. Marks: none.
- Input: `drop_tables` on a handle twice, on an unknown handle (-1, 10**12), `resident()` count before/after; `eval_sets` with a dropped handle; `add_tables_from_raters` with raters of k=1, k=33, an H or W of 0, mismatched H/W between raters.
- Expected: `ValueError` (never crash, never use freed memory); `resident()` returns to the baseline after drop; double drop raises or is a documented no-op, same as before Phase D (record the golden); a failed `add_tables_from_raters` leaves `resident()` unchanged (no leaked partial tables).
- Why: B16/B19. Cleanup on error paths is a common leak.

---------------------------------------------------------------------------------------------------------------

## H. Counters, determinism, failure modes (B6, B16, B19, G6, G15)

### 35. score_polygons_host_waits_and_allocations
- Level: integration. Marks: none (CPU counts the same events), `gpu`.
- Setup: batch of 64 polygons of one shape; `gc.reset_stats()`; call 1 (warm-up) then `reset_stats()`; calls 2..10; then `drop_tables` on the cache.
- Expected: call 1: `waits == 3`; calls 2..10, each: `waits == 3`, `device_allocs == 0`, `device_releases == 0`; `drop_tables`: `waits == 0`. `timing=True` call may add waits but the default call adds none: assert `waits` after a default call equals 3 even when interleaved with `timing=True` calls. Different batch size (larger) after warm-up may allocate; assert it does not exceed `arena_bytes` ceiling and the next same-shape call allocates 0.
- Why: B19/G8/G15. Hidden synchronisation is the performance bug the whole design avoids.

### 36. determinism_20_repeats
- Level: GPU. Marks: `gpu`.
- Input: one fixed batch per function: `check_masks`, `score_polygons`, `eval_sets`, `certify_polygons`, `topology`, `add_tables_from_raters`; 20 runs each, with a cuda stream idle/busy mix (run a second thread submitting `eval_sets` concurrently in the "busy" variant if the binding permits; else skip that variant).
- Expected: all 20 outputs byte-identical (`pickle`/`tobytes` hash set size 1). Also identical to the CPU build.
- Why: G6/C2. Races in atomics-based label propagation or float atomics give run-to-run drift.

### 37. cuda_error_and_no_device
- Level: failure mode. Marks: none; `gpu` variant.
- Input: (a) CUDA build on a no-device host (simulate with `CUDA_VISIBLE_DEVICES=-1` in a subprocess); (b) request that exceeds memory (mask batch with `chunk_pixels` huge and a 4096x4096 x many in the unsplit setting, or `max_scratch_bytes=1<<45`); (c) a mask list totalling > INT_MAX items (use a virtual test hook or call with `start` int64 array claiming N > 2**31 via a check on the argument before any allocation).
- Expected: (a) clean `RuntimeError`/documented error, subprocess exit code 0 after catching, never SIGSEGV; (b) `RuntimeError` including the CUDA error text (`out of memory` or similar), and the process still runs a following normal call correctly (no sticky error); (c) `ValueError`/`OverflowError` raised before any work, never a wrapped negative count.
- Why: B16.

### 38. engine_threading_and_reentrancy
- Level: adversarial. Marks: none.
- Setup: 4 Python threads each calling `check_masks` on different batches 50 times; one thread calls `trim(0)` and `reset_stats()` between.
- Expected: each result equals the single-thread result; no deadlock (timeout 60 s); `stats()` counters non-negative; `trim` never frees memory in use.
- Why: the arena/pool and global stats are shared state; the GIL may be released in the binding.

### 39. sanitizer_corpus (sanitizer)
- Level: GPU. Marks: `gpu`, `sanitizer` (file `tests/test_gpu_sanitize.py`).
- Input: fixed small corpus: 40 masks (empty, 1x1, 31/32/33 sides, nested, noise 64x64 and one 512x512), 60 polygons across the 13 defect classes, 3 eval_sets, 4 tables caches.
- Expected: `ERROR SUMMARY: 0 errors` under memcheck, initcheck, synccheck, racecheck, each inside 300 s; the corpus covers each public entry point once. Includes a deliberately out-of-range-sized case that raises ValueError before the kernel (sanitizer log has no error).
- Why: G5. Padding bits, last-word reads, and warp-synchronous shortcuts show up here.

---------------------------------------------------------------------------------------------------------------

## I. No-CPU-geometry and repo rules (G11, C5, C6)

### 40. no_cpu_geometry_patched
- Level: end-to-end. Marks: none.
- Setup: monkeypatch `cv2.connectedComponents`, `cv2.connectedComponentsWithStats`, `scipy.ndimage.label`, `scipy.ndimage.distance_transform_edt` to raise; LIDC and RIGA fixtures.
- Expected: `build_target` (LIDC, RIGA), `ExpertCache.build`/`load`/`device`, `certify`, `metrics()` run without error, with the same results as unpatched runs (computed beforehand). Static checks: `grep -R "tests.oracles\|tests/oracles" src/objdel` is empty; files `src/objdel/checker/polygonize.py` and `bitquad.py` do not exist; `objdel.checker` has no attribute `bitquads`; `from objdel.checker import trace, polygonize, Ring, MaskRings, certify, certify_many, as_mask, betti, instance_topology, min_area, signature, size_filter, signatures, topologies` works.
- Why: G11. The patching is the only way to prove the CPU paths are gone.

### 41. gpu_marker_and_dry_run_gate
- Level: integration. Marks: none.
- Expected: `pytest --collect-only -q -m gpu` collects >= 1 test on the Mac; `make gpu-check DRY=1` prints the chain (make cuda, pytest -m gpu, memcheck/initcheck/synccheck/racecheck with `timeout 300`, bench_checker, old-vs-new replay) and exits 0; exits non-zero when one named test file is renamed (rename test in a tmp copy). On a GPU host with `gc.build() != "cuda"`, the gpu tests FAIL (not skip).
- Why: spec 3.3.

---------------------------------------------------------------------------------------------------------------

## J. Performance gates (G7, G8, G14), GPU box only

### 42. polygonizer_speed_vs_old
- Level: performance. Marks: `gpu`, `slow`.
- Setup: one RTX 5090 rental, alternating old/new runs, 1 warm-up, 15 timed runs per side; workloads: 256 frames 512x512 noise, 256 frames 512x512 LIDC-like blobs, 37,093 LIDC crops, 4096x4096 sparse, 4096x4096 dense, 100,000 tiny (8x8) masks, one 1x4096 mask.
- Expected: outputs identical to the old build (hash equal); median new <= 1.05 x old on every workload; <= 0.9 x old on the 256x512x512 and on the 37,093 crops. Use `cudaEvent` or wall time with `torch.cuda.synchronize`-like fence; record P95 and min.
- Why: G7. Alternating runs remove thermal and clock drift bias; the gate measures what the rewrite is for.

### 43. work_proportional_scaling
- Level: performance (property). Marks: `gpu`, `slow`; a CPU variant with op counters if exposed.
- Setup: mask 4096x4096 with (a) all zero, (b) one pixel, (c) one solid rectangle (few runs, few edges), (d) 1000 small dots, (e) noise 0.5.
- Expected: time(a) <= ~ time of words-scan only: assert time(a) < 0.3 x time(e) and time(b) < 0.3 x time(e); time grows roughly with runs + edges, not with area: time(c) < 0.3 x time(e) (a solid rectangle has the same number of words but few edges); report ratios. Fail when any "cheap" case exceeds the ratio.
- Why: spec section 1 states work proportional to words + runs + boundary edges. A kernel that does per-pixel work shows no difference between (c) and (e).

### 44. reward_path_replay_speed
- Level: performance. Marks: `gpu`, `slow`.
- Setup: replay set of 256 completions on their caches; old vs new alternating.
- Expected: new `score_polygons` median <= 1.05 x old, `waits == 3` (B19), results bit-identical to the old.
- Why: G8.

### 45. audit_memory_bound_in_rl
- Level: end-to-end performance. Marks: `gpu`, `slow`.
- Setup: 6 RL steps with the audit callback on.
- Expected: `gc.stats()['arena_bytes'] <= eval_sets scratch cap (256 MiB default) + table bytes`; no OOM; steps 2-6: `device_allocs` delta == 0 and `device_releases` delta == 0 (G15); audit time and peak reserved memory recorded in the log.
- Why: G14/G15.

---------------------------------------------------------------------------------------------------------------

## K. Corpus-scale regression (G1, G2, G3), data-backed

### 46. checker_eval_crops_exact
- Level: regression. Marks: `data`, `slow`, `gpu` variant.
- Input: `data/checker_eval/masks.npz` read with `python -I` (untrusted data). 37,093 crops.
- Expected: total rings 38,002, total corners 1,476,508 (literals from spec G1); every crop equals the oracle (compare ring lists); part labels equal cv2 B3 on every crop; CPU build == GPU build.
- Why: G1/G2/G3; the literals pin the counts so an oracle bug cannot pass quietly.

---------------------------------------------------------------------------------------------------------------

# Pytest sketches (the 5 most important)

Importance order: (1) polygonizer exactness + part labels, (2) batch/chunk independence, (3) certificate independence mutation, (4) host waits/allocation counters, (5) EDT bit identity. All five are CPU-runnable except where `gpu` is stated.

```python
# tests/test_phase_d_core.py
import numpy as np
import cv2
import pytest
import scipy.ndimage as ndi

from objdel import _gpucheck as gc
from objdel.checker import rings as R
from tests.oracles.polygonize_oracle import polygonize as oracle_poly
from tests.oracles.bitquad_oracle import bitquads as oracle_quads

def canon(rs):
    return [(r.pts.tolist(), int(r.area), int(r.part)) for r in rs]

def families(seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for d in np.arange(0.05, 0.96, 0.1):
        for side in (1, 2, 31, 32, 33, 64, 65):
            out.append(rng.random((side, int(rng.integers(1, 70)))) < d)
    out += [np.zeros((0, 0), bool), np.zeros((3, 3), bool), np.ones((33, 33), bool),
            (np.indices((16, 16)).sum(0) % 2).astype(bool)]            # checkerboard
    sq = np.zeros((21, 21), bool)
    for i in range(0, 10, 2):
        sq[i:21 - i, i:21 - i] = True
        sq[i + 1:20 - i, i + 1:20 - i] = False                        # nested frames with holes
    out.append(sq)
    return out

# 1. exactness vs oracle + part labels vs cv2 + Euler identities
def test_trace_equals_oracle_parts_and_euler():
    masks = families()
    got = R.trace(masks)                                              # one call
    for m, g in zip(masks, got):
        exp = oracle_poly(m)
        assert [(p.tolist(), a) for p, a, _ in canon(g.rings)] == [
            (e.pts.tolist(), int(e.area)) for e in exp]               # corners, start, order, sign
        assert g.quads == tuple(int(v) for v in oracle_quads(m))
        if m.size == 0 or not m.any():
            assert g.rings == [] and (g.b0, g.b1) == (0, 0)
            continue
        n, lab = cv2.connectedComponents(m.astype(np.uint8), connectivity=8)
        for r in g.rings:
            x, y = r.pts[0]                                           # first corner = first east edge start
            assert r.part + 1 == lab[y, x]
        outer = [r for r in g.rings if r.area > 0]
        assert len(outer) == g.b0 == n - 1
        assert len({r.part for r in outer}) == len(outer)
        assert g.b1 == len(g.rings) - len(outer)
        bg, nb = ndi.label(~np.pad(m, 1, constant_values=False))      # 4-connected background
        assert g.b1 == nb - 1
        assert sum(r.area for r in g.rings) == int(m.sum())
        assert g.b0 - g.b1 == (g.quads[0] - g.quads[3] - 2 * g.quads[2]) // 4

# 2. batch / chunk / scratch / order independence
@pytest.mark.parametrize("chunk", [1, 64, 4096, 65536, 1 << 25])
@pytest.mark.parametrize("scratch", [4096, 65536, 1 << 30])
def test_results_independent_of_batching(chunk, scratch):
    masks = families(3)[:40]
    bits, hw, st = pack(masks)
    ref = [gc.check_masks(*pack([m])) for m in masks]                 # each alone
    out = gc.check_masks(bits, hw, st, chunk_pixels=chunk, max_scratch_bytes=scratch)
    for i, r in enumerate(ref):
        a, b = out["ring_start"][i:i + 2]
        assert out["ring_area"][a:b].tobytes() == r["ring_area"].tobytes()
        pa, pb = out["pts_start"][a], out["pts_start"][b]
        assert out["pts"][pa:pb].tobytes() == r["pts"].tobytes()
        assert out["ring_part"][a:b].tobytes() == r["ring_part"].tobytes()
        assert out["quads"][i].tobytes() == r["quads"][0].tobytes()
        assert out["betti"][i].tobytes() == r["betti"][0].tobytes()
    perm = np.random.default_rng(0).permutation(len(masks))           # item order independence
    out2 = gc.check_masks(*pack([masks[j] for j in perm]))
    assert out2["quads"][np.argsort(perm)].tobytes() == out["quads"].tobytes()

# 3. certificate independence (needs the test-only mutated build; skipped when absent)
@pytest.mark.skipif(not hasattr(gc, "_test_four_adjacent_labeller"), reason="mutation build absent")
def test_certify_catches_a_wrong_ring_labeller():
    from objdel.checker import certify
    m = np.zeros((2, 2), bool); m[0, 0] = m[1, 1] = True
    assert certify(m, 0).certified is True                            # normal build
    gc._test_four_adjacent_labeller(True)
    try:
        assert certify(m, 0).certified is False                       # topology chain (pixels) still says 1 component
    finally:
        gc._test_four_adjacent_labeller(False)

# 4. counters: 3 waits, then 0 allocations, drop_tables 0 waits (CPU build counts the same events)
def test_score_polygons_waits_and_allocations(score_batch_and_cache):  # fixture: args tuple + table handle
    args, handles = score_batch_and_cache
    gc.score_polygons(*args)                                          # warm-up
    for _ in range(9):
        gc.reset_stats()
        gc.score_polygons(*args)
        s = gc.stats()
        assert (s["waits"], s["device_allocs"], s["device_releases"]) == (3, 0, 0)
    gc.reset_stats()
    for h in handles:
        gc.drop_tables(h)
    assert gc.stats()["waits"] == 0

# 5. EDT bit identity and _BIG
from objdel.reward import cache as rc

def test_edt_bits_equal_scipy():
    rng = np.random.default_rng(5)
    maps = [rng.random((int(rng.integers(1, 300)), int(rng.integers(1, 300)))) > d
            for d in np.linspace(0.001, 0.999, 60)]
    z = np.ones((9, 9), np.uint8); z[0, 8] = 0
    maps += [z, np.ones((1, 7), np.uint8), np.zeros((5, 1), np.uint8)]
    for nz in maps:
        nz = nz.astype(np.uint8)
        got = gc.edt(nz)
        assert got.dtype == np.float32
        if nz.all():
            assert (got == rc._BIG).all()
        else:
            assert got.tobytes() == ndi.distance_transform_edt(nz).astype(np.float32).tobytes()
```

`pack(masks)` returns `(bits, hw, start)` per section Conventions. Fixture `score_batch_and_cache` builds 64 polygons of one shape and their cache tables.

## Coverage map (gate to cases)
G1: 1,2,3,4,5,11,46. G2: 6,46. G3: 21,26,32,36. G4: 19. G5: 39. G6: 36. G7: 42,43. G8: 44. G9: 29. G10: 30,32. G11: 33,40. G12: 20-25 plus `make cov` (not a test). G13: no test (ledger). G14/G15: 35,45. G16: 14. B16 limits: 11,13,25,34,37. B17: 16,17,18. B19: 35.

## Known gaps (cannot be tested from the spec alone)
- Spec gives no name for the mutated-build switch (test 14) or a table dump function (test 32); sketch assumes `gc._test_four_adjacent_labeller` and read-back through reward terms.
- Spec leaves open whether double `drop_tables` raises (test 34): the golden is the pre-Phase-D behaviour.
- `chunk_pixels=1` behaviour when one mask exceeds the chunk is unspecified: test accepts correct output or a clear ValueError.
