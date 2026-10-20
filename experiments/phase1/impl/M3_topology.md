# M3: per-instance topology and the size rule

Plan: D3, §4.2, gate R6. Promoted from `lab/t_reward/rewardcore/topology.py`.

**Code:**
- `src/objdel/checker/topology.py`
- `src/objdel/checker/__init__.py`
- `tests/test_topology.py`
- `tests/oracles/topology_oracle.py`
- `experiments/phase1/census_sizerule.py`

## Result

| Gate / check | Result |
|---|---|
| Exactness vs slow oracle (scipy, each component labelled on its own) | 0 mismatches on 10,000 random masks. Shapes 1×1 to 200×200, 10 families. |
| Exactness on real fixture masks, each mask and its complement | 0 mismatches on 2,400 RIGA, 1,969 LIDC and 35,124 full-LIDC masks |
| Exactness vs pure-Python BFS oracle | 0 mismatches on 300 random + 300 Hypothesis small masks |
| χ certificate: bit-quad χ8 = Σ(1 − holes_k) | Holds on every mask above |
| Production vs lab code (topology and filtered signature) | 0 / 78,986 fixture masks and complements (sandbox) |
| R6 size-rule invariance | Pass: specks/pinholes of A−1 px leave the signature set unchanged; A px change it. Filter = oracle; idempotent. |
| Determinism | Identical arrays across repeated calls and OpenCV 1 vs 4 threads (parallel labeling path) |
| `test_topology.py` | 45 tests pass, about 35–60 s on the loaded Mac. The full-LIDC case is marked `slow`. |
| `make lint` on M3 files | Clean |

## Research

**Spaghetti labeling and OpenCV.**
- The algorithm is Bolelli, Allegretti, Baraldi, Grana, IEEE TIP 29:1999–2012, 2020 ([record](https://iris.unimore.it/handle/11380/1181437), [PDF](https://iris.unimore.it/retrieve/e31e124d-c7be-987f-e053-3705fe0a095a/2018_TIP_Spaghetti%20Labeling%20Directed%20Acyclic%20Graphs%20for%20Block-Based%20Connected%20Components%20Labeling.pdf)). It is a block-based decision forest compressed into a DAG. On YACCLAB it beats SAUF, BBDT, DRAG and the others.
- Flag mapping, read from the OpenCV 5.x [imgproc.hpp](https://github.com/opencv/opencv/blob/5.x/modules/imgproc/include/opencv2/imgproc.hpp) and the dispatcher `connectedComponents_sub1` in [connectedcomponents.cpp](https://github.com/opencv/opencv/blob/5.x/modules/imgproc/src/connectedcomponents.cpp):

  | Flag | 8-connectivity | 4-connectivity |
  |---|---|---|
  | `CCL_SPAGHETTI` / `CCL_BOLELLI` / `CCL_DEFAULT` | Spaghetti | Spaghetti4C (Bolelli 2021) |
  | `CCL_BBDT` / `CCL_GRANA` | BBDT | SAUF (falls back) |
  | `CCL_SAUF` / `CCL_WU` | SAUF | SAUF |

- **Stats.** Every path supports stats, because one templated `StatsOp` serves all of them.
- **Parallel labeling.** The parallel variants (Bolelli 2017) run when `getNumThreads() > 1` and `rows / threads >= 2`.
- **Types.** Labels must be `CV_16U` or `CV_32S`.
- **Input types.** The 5.x branch also accepts `CV_Bool` input. The installed 5.0.0 wheel rejects numpy bool ("iDepth == CV_8U || iDepth == CV_8S", observed), so the code passes zero-copy `uint8` views. The 5.0.0 Python API also rejects the `ccltype=` keyword (lab note), so arguments are positional.
- **OpenCV 5 migration** ([wiki mirror](https://lab.trianta.dev/mirrors/opencv.opencv/wiki/OpenCV-4-to-5-migration)): 1-D numpy arrays now map to 1-D `Mat`. We always pass 2-D arrays, and 1×N and N×1 masks are tested.

**Alternatives.**
- [cc3d](https://github.com/seung-lab/connected-components-3d): SAUF-based. It is not installed; the earlier local benchmark in `docs/research/notes/Phase 1 stack and algorithms/mask_algorithms.md` found it about 3× slower than cv2 at 64².
- scikit-image `label`: costs about the same as scipy.
- `skimage.measure.euler_number`: about 2× slower than the numpy bit-quads.
- `skimage.morphology.remove_small_holes` ([docs](https://scikit-image.org/docs/stable/api/skimage.morphology.html), source read locally). Two pitfalls, so not used:
  - It does not pad, so background touching the border counts as a "hole".
  - In 0.26 the new `max_size` removes holes **≤** the size, while the old `area_threshold` removed holes **<** it.
- `findContours(RETR_CCOMP)` (Suzuki–Abe 1985, [OpenIMAJ notes](https://openimaj.org/apidocs/org/openimaj/image/contour/SuzukiContourProcessor.html)): gives exact hole counts per outer contour (0 mismatches in the sandbox). It does not give hole pixel counts, which the size rule needs: `contourArea` is a polygon area through pixel centres.

**Digital topology.** The (8, 4) pair (8-connected foreground, 4-connected background) is the pair for which a digital Jordan curve theorem holds, so holes are well defined.
- A. Rosenfeld, "Connectivity in digital pictures", J. ACM 17(1):146–160, 1970 ([JACM index](https://projects.csail.mit.edu/jacm/Authors/rosenfeldazriel.html)).
- T. Y. Kong, A. Rosenfeld, "Digital topology: introduction and survey", CVGIP 48(3):357–393, 1989 ([doi](https://doi.org/10.1016/0734-189X(89)90147-3)).

**Size rule.**
- The filter is a binary area opening followed by its dual, the area closing. Source: L. Vincent, "Morphological area openings and closings for grey-scale images", Shape in Picture, NATO ASI F 126, 1994 ([doi](https://www.doi.org/10.1007/978-3-662-03039-4_13)).
- In a binary area closing, a background component's area is its own pixel count.

## Sandbox (scratch `impl/M3/`, Apple M3, single thread, load average about 2.6)

**Exactness.** 3,000 random masks plus 1,600 RIGA and 1,600 LIDC fixture masks/complements gave 0 mismatches for every candidate.

**Speed.** Median µs per mask:

| set | A: lab (2× Spaghetti + pixel-above) | D: `WithStats` + hole top row | B: findContours CCOMP (no hole areas) | S: A with scipy labels |
|---|---|---|---|---|
| LIDC fixture | 9.1 | 11.0 | 10.8 | 21.7 |
| RIGA fixture | 110 | 186 | 99 | 205 |
| RIGA complement | 828 | 2,142 | 680 | 1,591 |
| 256² blob | 182 | 302 | 147 | 327 |
| 256² noisy | 383 | 752 | 538 | 704 |
| 1024² blob | 2,083 | 4,180 | 1,766 | 3,842 |

**Where the time goes, and the variants it rules out:**
- **`np.bincount` dominates.** Over the int32 label image it costs more than the labeling itself. On a RIGA crop: label8 12 µs, label4 19 µs, bincount 53 µs each.
- **Single component.** The production code skips the foreground bincount when there is one component (`count_nonzero`). It bincounts the background only when a hole exists.
- **`connectedComponentsWithStats` is slower** (variant D): OpenCV's stats op also accumulates bounding boxes and centroids.
- **Cropping.** `cv2.boundingRect` on a `uint8` view plus `cv2.copyMakeBorder` crops faster than numpy `any()` for small masks: LIDC 6.4 vs 9.2 µs. It was checked against numpy on 20,000 random shapes and layouts, including empty, 1×N, Fortran-ordered and strided masks.
- **Threads.** OpenCV's default 8 threads vs 1 thread:
  - LIDC: 11.6 vs 9.2 µs.
  - RIGA: 136 vs 110 µs.
  - 1024² blob: 2.0 vs 2.5 ms.

  Small crops are faster single-threaded. Results are identical either way, and a test checks it.

**Production vs lab vs scipy.** Median µs; separate run:

| set | production | lab | scipy, same algorithm | scipy, per-component oracle | production `signature(m, 7)` | lab `filtered_signature` |
|---|---|---|---|---|---|---|
| LIDC fixture | 6.7 | 8.2 | 22.6 | 31.1 | 8.8 | 11.9 |
| RIGA fixture | 61.5 | 57.9 | 209 | 856 | 63.3 | 61.8 |
| RIGA complement | 833 | 849 | 1,631 | 1,628 | 858 | 892 |
| 256² blob | 126 | 131 | 332 | 450 | 137 | 140 |
| 256² noisy | 602 | 623 | 721 | 2,190 | 1,820 | 1,892 |
| 1024² blob | 1,162 | 1,122 | 4,160 | 6,285 | 1,191 | 1,168 |

**The lab's "5–24× vs scipy" is not confirmed for an equal algorithm:**
- Against scipy running the same two-labeling algorithm, Spaghetti is 1.2–3.6× faster.
- Against the per-component scipy approach it is 2–14× faster.
- Production matches the lab within noise; it is about 20 % faster on small LIDC masks because of the cv2 crop.

## Decisions

1. **Algorithm: the lab's.**
   - Crop to the bounding box, pad by 1 pixel.
   - Label the foreground with Spaghetti (8-connected) and the background with Spaghetti4C (4-connected).
   - Assign each hole to the component on the pixel directly above its first raster pixel (`np.unique` over the boundary-sized candidate set).
   - It is the fastest exact option that also returns hole pixel counts.
2. **Canonical output.**
   - `InstanceTopology(areas, holes, hole_areas)` is sorted: components by (area, holes) descending, hole areas descending.
   - The result is a function of the mask alone, not of label numbering, which can differ on OpenCV's parallel path. This serves K10.
   - The lab's `labels`/`origin`/`full_labels`/`as_list`/method `signature(min_area)` are dropped: unused, and the method only filtered components, not holes.
3. **How a hole's size is measured.**
   - A hole's size is its own background pixel count. It is measured after sub-filter components are removed, so a speck inside a hole counts toward the hole.
   - Islands inside a hole do not count toward its size.
   - Filling a sub-filter hole merges the islands inside it into the enclosing component.
   - This is the lab and Phase 0 rule. The reasons to keep it:
     - It is the complement dual of the component rule, which also counts only a component's own pixels.
     - It keeps the census numbers comparable.
     - The alternative, enclosed area, would keep a hairline moat (32 px) around a large island (49 px) as a hole. That is exactly the annotation crack the filter should suppress. A test covers this case.
   - The filter is idempotent: every surviving component and hole has ≥ min_area px. Tested with Hypothesis.
4. **`min_area`** = `max(5, -(-|U| // 200))`, exact integer arithmetic. It equals the float formula for every |U| ≤ 2·10⁶ (checked).
5. **`signature(mask, min_area)`**:
   - It computes the topology once.
   - It re-runs on `size_filter(mask)` only if some component or hole is below the filter (the lab shortcut).
   - `signature(mask)` is the raw signature.
6. **Reconciliation with the old module:**
   - `betti(mask)` is kept, now derived from `instance_topology`, because `certify.py` and `reward/annotators.py` use it (`annotators.py` was deleted in M6; `certify.py` remains the library user).
   - `size_filter` is replaced by the OpenCV version with the same semantics; the scipy version moved to the oracle.
   - `betti3d` is removed: no user outside `experiments/phase0`; Phase 1 is 2D.
   - `FG8`/`BG4` are removed: only `experiments/phase0` used them.
   - The scipy import is gone from `objdel.checker`.
   - `bitquad.py`, `certify.py` and `polygonize.py` are unchanged.

## Tests (`tests/test_topology.py`)

- **Oracle equality.**
  - 10,000 random masks from 10 families:
    - i.i.d. pixels at 5–95 % density;
    - smooth noise and its complement;
    - blobs with 1-px holes;
    - checkerboards;
    - nested rings (island in hole in island);
    - diagonal-only rings;
    - dots inside a frame;
    - noisy crop edges;
    - 1-px lines.
  - Fixture masks and their complements (`@pytest.mark.data`, `fixtures_dir`).
  - Each check compares (area, holes) in canonical order, `hole_areas`, Σ areas = pixel count, `betti`, and χ8.
- **Size rule (R6):**
  - The `min_area` table.
  - Sub-filter specks and pinholes leave the rater signature set unchanged, while raw signatures see them.
  - Filter-size specks and holes change the set.
  - A speck inside a hole counts toward the hole.
  - The moat case.
  - `size_filter` and `signature` equal the oracle on 2,000 random masks plus 300 Hypothesis masks, and the filter is idempotent.
- **Edge cases:**
  - empty masks: 0×0, 0×4, 1×1, 1×7, 7×1;
  - full masks;
  - masks touching every border;
  - 1×N and N×1;
  - dtypes `uint8`/`int64`, Fortran-ordered and strided inputs;
  - 3-D input raises `ValueError`.
- **Determinism:** repeated calls, and OpenCV 1 vs 4 threads on 512×384 masks.
- **Hypothesis properties:** pure-Python BFS equality, and invariance under the 8 grid symmetries plus translation/padding (topology, signature and filter).

## Census re-run under the size rule

Command: `.venv/bin/python -I -m experiments.phase1.census_sizerule` (26 s, 3 processes).

An item is a slice (LIDC) or an image structure (RIGA) with k ≥ 2 raters. An item counts as a disagreement when the raters' keys differ. The keys are:
- `p0_raw`: global (b0, b1).
- `p0_fixed`: global (b0, b1) after Phase 0's fixed filter (5 px LIDC, 100 px RIGA).
- `sig_raw`: per-instance signature.
- `sig_rule`: signature with `min_area(|U|)`.

| items | n | p0_raw | p0_fixed | sig_raw | sig_rule | masks changed by the rule | median min_area |
|---|---|---|---|---|---|---|---|
| LIDC full, all slices | 10,921 | 515 (4.7 %) | 456 (4.2 %) | 515 (4.7 %) | 446 (4.1 %) | 159 | 5 |
| LIDC full, train (no duplicated reader) | 7,562 | 293 (3.9 %) | 256 (3.4 %) | 293 (3.9 %) | **251 (3.3 %)** | 96 | 5 |
| LIDC full, val (no duplicated reader) | 1,117 | 51 (4.6 %) | 42 (3.8 %) | 51 (4.6 %) | **37 (3.3 %)** | 34 | 5 |
| LIDC full, test (no duplicated reader) | 2,004 | 96 (4.8 %) | 87 (4.3 %) | 96 (4.8 %) | **87 (4.3 %)** | 18 | 5 |
| LIDC full, nodules (any slice) | 1,879 | 174 (9.3 %) | 154 (8.2 %) | 174 (9.3 %) | 154 (8.2 %) | 159 | 5 |
| fixture lidc_experts | 602 | 14 (2.3 %) | 5 (0.8 %) | 14 (2.3 %) | 5 (0.8 %) | 9 | 5 |
| fixture riga_experts | 400 | 0 | 0 | 0 | 0 | 0 | 92 |
| RIGA true (Phase 0b), disc | 698 | 0 | 0 | 0 | 0 | 0 | 186 |
| RIGA true (Phase 0b), cup | 698 | 0 | 0 | 0 | 0 | 0 | 56 |
| RIGA true (Phase 0b), rim | 698 | 106 (15.2 %) | 0 | 106 (15.2 %) | 0 | 163 | 158 |

*Superseded for LIDC by the final store (`fixes/F4.md`: copied reader outlines merged, B3-3; crowded frames dropped,
B3-2): 10,625 items; lab-split train disagreement 296 (5 px) / 291 (rule) / 242 (rule, no duplicated reader), was
313 / 308 / 251 (`data/store/build_report.json`). B2-4 (F2) further drops items whose only disagreement is an
abstaining sub-filter rater (the 251-item gate set becomes 239). The rows above are the lab fixture `lidc_full.npz`.*

**Comparison with Phase 0** (`experiments/phase0/RESULTS.md`, `results/lidc/summary.json`):

- **The lab numbers are reproduced exactly:**
  - full-LIDC disagreement items under the rule: train 251, val 37, test 87 (ITERATION3 R6);
  - fixtures: 14 → 5 disagreement items, 9 rater masks changed.
- **LIDC slices, raw.** Phase 0 reported 632 / 10,916 (5.8 %); the re-run gives 515 / 10,921 (4.7 %).
- **LIDC slices, filtered.** Phase 0 reported 317 / 10,898 (2.9 %) at a fixed 5 px; the re-run gives 446 (4.1 %) under the rule and 456 with Phase 0's fixed 5 px. The re-run is not a pure re-run of the Phase 0 script, for three reasons:
  - **Different rasterization.** Phase 0 rasterized the XML with PIL (outline included, exclusion interiors only). The fixtures come from the Phase 1 parser, which reaches pylidc parity. The raw differences come from that.
  - **Presence counts here.** Phase 0 dropped a reader whose slice became empty after filtering. Here an empty filtered slice has signature `()` and counts as a disagreement, as in the lab's reward and census.
  - **The Phase 0 script cannot run any more.** It imports `betti3d`, removed here.
- **Rule vs fixed 5 px on LIDC.** The rule gives 10 fewer disagreements (446 vs 456). The cause is the 314 slices with |U| > 1,000 px, where the filter rises above 5 px.
- **Per-instance vs global.** The per-instance signature and the global (b0, b1) give the same counts on LIDC: 515 raw, and 456 at a fixed 5 px (checked separately). No disagreement on LIDC depends on how holes are spread across parts.
- **LIDC nodules.**
  - Re-run: 154 / 1,879 (8.2 %) of nodules have at least one disagreeing slice under the rule.
  - Phase 0: 114 / 1,877 (6.1 %) "any topology disagreement" at 5 px.
  - The definitions differ: Phase 0 used per-reader "any hole slice" / "any multi-part slice" flags, plus 3-D components and reader splits, not a slice-by-slice comparison.
- **RIGA.**
  - Disc and cup: 0 disagreement under every rule, as in Phase 0/0b.
  - Rim, raw: 106 images (15.2 %) disagree, from slivers where the cup touches the disc edge.
  - Rim, under the rule or the fixed 100 px: 0. This matches Phase 0b's 0 at 100 px.

## Contract changes (EXECUTION.md §3 updated)

- **`InstanceTopology`** has the fields `areas`, `holes`, `hole_areas` (int64 arrays) in canonical order. It has no methods.
- **`betti(mask)`** is kept. `betti3d`, `FG8` and `BG4` are removed.
- **`size_filter`** is documented as an area opening (8-connected) followed by an area closing (4-connected).
- **New exports** from `objdel.checker`: `InstanceTopology`, `instance_topology`, `min_area`, `signature`.
- *Later (F1, I1 / B1-5):* `as_mask(mask, ndim=2)`, the one normalisation every function calls (any nonzero byte is
  foreground).

## Open issues

- **OpenCV dependency.** `objdel.checker` now imports `cv2`, but `opencv-python-headless` is only in the `models` extra. It must move to core dependencies, pinned at 5.0.0 (coordinator). *(Resolved: core dependency `opencv-python-headless==5.0.0.93`; the `models` extra is gone.)*
- **Signature reuse in M6.** `signature(m, a)` and `instance_topology(m)` compute the topology separately. If M6 needs the raw and filtered signature of the same rollout and K3 gets tight, add an optional precomputed `InstanceTopology` argument. Not added now (no evidence of need).
- **Filter order.** The filter removes components first, then fills holes. The order matters only at tiny scale; for example, a 4-px diamond ring around a 1-px hole is removed at min_area 5, but the opposite order would keep it at 5 px. Kept for consistency with Phase 0 and the lab.
- **Timings.** All timings come from a loaded Mac. The x86 re-measurement belongs to M9.
- **OpenCV threads.** Callers that score in worker pools should call `cv2.setNumThreads(0)` (corrected in fix F1, issue I14: `setNumThreads(1)` is a no-op on the GCD backend of the macOS wheels, where `getNumThreads()` stays 8; 0 is sequential on every backend). Sequential is faster on small crops, and the results are identical. The "1 vs 4 threads" determinism runs above therefore compared 8 vs 8 threads on macOS; the test now uses 0 vs 4.
