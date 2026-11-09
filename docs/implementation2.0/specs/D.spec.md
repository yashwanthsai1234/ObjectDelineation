# Phase D spec: CUDA core (contract for test writers, revision 2)

REPO = $REPO. All paths below are REPO-relative.
`gc` = `objdel._gpucheck` (one C++ source, `cuda/gpucheck.cpp`. The CPU build is on the Mac and the nvcc build is on
a CUDA box. Both builds must return the same bits). Revision 2 follows the review of D.md (log in D.md §9).

## 1. Goal

The checker (validity and topology), the scorer and the polygonizer run in the CUDA engine in production. The mask
polygonizer does work proportional to the input words plus runs plus boundary edges. Expert tables are built on the
device. Evaluation keeps pixels on the device and returns only integers. The CPU build of the same source gives
identical results.

## 2. Required behaviours

B1. `rings.trace` returns, for each mask, the rings that the numpy oracle `tests/oracles/polygonize_oracle.polygonize`
    returns. Same ring order, same corner lists from the same start corner, same signed areas (> 0 outer, < 0 hole).
    It also returns bit-quad counts equal to `tests/oracles/bitquad_oracle.bitquads` (n1, n2, nD, n3, n4).
B2. Ring order: rings sort by the raster index of their first east edge. Each ring starts at its minimum edge key
    (y·(W+1) + x)·4 + d (0 east, 1 south, 2 west, 3 north). Foreground lies right of every edge. Connectivity: foreground
    8-connected, background 4-connected. The image border counts as background.
B3. `Ring.part` = the 0-based raster rank (by first pixel) of the foreground 8-component that the ring bounds. An outer
    ring gives its own component. A hole gives the component that surrounds it. `part + 1` equals the label of
    `cv2.connectedComponents(mask_uint8, connectivity=8)` at the foreground pixel right of the ring's first edge.
B4. `check_masks` `betti[m]` = (number of foreground components, number of bounded background components).
    `b0 - b1 == chi8` from the bit-quads. The number of outer rings is b0, and the number of hole rings is b1.
B5. Empty masks (no foreground pixel), 0-width or 0-height masks, 1x1, 1xN, Nx1, and full masks are valid inputs.
    They give the oracle's output (no rings for empty masks).
B6. Results do not depend on batch composition, chunking, scratch splitting, item order or the thread schedule. The
    same mask gives the same output alone and in any batch, with any `chunk_pixels` / `max_scratch_bytes` test value
    (for example 4096 and 65536 on 40 mixed masks), and in every run (byte-identical over 20 repeats on the GPU).
    The same holds for `eval_sets` under any scratch cap.
B7. `polygonize(mask, offset)` = `trace([mask], [offset])[0].rings`. Offsets add to every corner (x + ox, y + oy).
B8. `certify(mask, min_size)` keeps its `Report` fields and meanings. `certified` is True exactly when bit-quad
    chi8 == b0 - b1, (outer, holes) == (b0, b1), the ring area sum == the bit-quad area == the pixel count, and the
    corner count == the bit-quad corners. Here b0 and b1 come from the pixel label chain (`gc.topology`, B17), never
    from `check_masks`' `betti` (which shares the rings' labeller). `certify_many` returns the list of `certify`
    results. With a test-only build in which the ring labeller joins foreground only by 4-adjacency, `certify` of a
    mask with two diagonally touching pixels returns `certified == False`; the normal build returns True.
B9. Target texts do not change. `build_target(...)` returns the same `.text` as before Phase D for every input. Every
    recorded `target_set_sha256` in these files reproduces byte-identically:
    `data/prepared/manifest.json` (sets lidc/train, lidc/val, riga/train, riga/val), `data/targets/manifest.json`,
    `data/eval_views/lidc-val-shifted/manifest.json`, `data/eval_views/lidc-val-rescaled/manifest.json`,
    `data/eval_views/lidc-val-rotated/manifest.json`, `data/views/lidc-train-shifted/manifest.json`,
    `data/views/lidc-train-rescaled/manifest.json`.
B10. `build_target` requires `min_area` for every name (and "rim" for RIGA). A missing entry raises `ValueError`.
B11. `gc.certify_polygons` per shape: `legal` equals "no defect" from `check_polygons`. `inter` and `union` equal the
    popcounts of (candidate & reference) and (candidate | reference). `equal` is True exactly when the masks are equal.
    `sig` and `ref_sig` equal `objdel.checker.topology.signature(mask, min_area)` for the candidate raster and the
    reference. Per pair (a, b) (RIGA rim): `pair_sig` equals `signature(raster_a & ~raster_b, pair_min_area)` and
    `pair_ref_sig` equals `signature(ref_a & ~ref_b, pair_min_area)`.
B12. Validity totals (all 13 defect classes, per-ring counts, the ring-pair matrix) from `check_polygons` and
    `score_polygons` equal the values before Phase D on every input that the parser and binding accept.
B13. Scorer outputs (`score_polygons`, `score_masks`) are bit-identical to the values before Phase D, including when
    the size filter removes nothing and when it removes components or fills holes.
B14. `metrics()` (`objdel.eval.metrics`) returns the same fields as before. All fields are equal (iou_vote,
    dice_vote, d_max, precision_vote, recall_vote, area, a_amax, a_amin, m_e, components, extra_components, ged,
    hm_iou, soft_dice, duplicate share), except `nsd`, which is equal within 1e-12 absolute. `metrics()` never calls
    `objdel.geometry.to_masks`. Duplicates are defined as before: equal format and legal flags and equal RAW sample
    masks (before the legal gate) for every object.
B15. `objdel.reward.cache.edt(nonzero)` is bit-identical to `scipy.ndimage.distance_transform_edt(nonzero).astype(float32)`.
    It returns the existing `_BIG` constant everywhere when `nonzero` has no zero pixel. `reward.spec.fingerprint()`
    does not change.
B16. Errors: inputs outside the limits raise `ValueError` (Python) / `std::invalid_argument` (C++), never crash.
    Limits: mask sides 0..4096, polygon W and H 1..4096, G 1..2^30, at most 10 rings and 1,000 vertices per ring,
    2 <= k <= 32 raters for tables. A scan over more than INT_MAX items raises (never wraps). A CUDA failure raises
    `RuntimeError` with the CUDA error text. On a box with no device, the CUDA build raises a clean error (it does not
    crash).
B17. `objdel.checker.topology` `instance_topology`, `betti`, `size_filter` and `signature` return the same values and
    types as the cv2 versions in `tests/oracles/topology_oracle.py`, for every mask and every `min_area` >= 0
    (components in raster order of their first pixel; signature sorted descending).
B18. Device expert tables: for every cache, the tables built by `gc.add_tables_from_raters` are byte-equal to the
    tables `gc.add_tables` stores from the host-built inputs (distance map, union, core, vote, rater areas, voter flags,
    signatures, n_vote, perimeter), and the reward terms scored against them are bit-identical.
    Constructing or loading an `ExpertCache` makes no engine call; the first `device()` call makes the tables.
B19. Engine counters (`gc.stats()`): a `score_polygons` call makes exactly 3 host waits. After one warm-up call,
    repeated calls of the same shape make 0 device allocations and 0 device releases, and `drop_tables` makes 0 waits.
    The CPU build counts the same events.

## 3. Interfaces

### 3.1 `src/objdel/checker/rings.py`

```python
@dataclass(frozen=True)
class Ring:
    pts: np.ndarray      # (c, 2) int64 corners (x, y), canonical start first, offset applied
    area: int            # signed pixel area: > 0 outer ring, < 0 hole
    part: int            # B3
    @property
    def hole(self) -> bool: ...   # area < 0

@dataclass(frozen=True)
class MaskRings:
    rings: list[Ring]                         # B1/B2 order
    quads: tuple[int, int, int, int, int]     # n1, n2, nD, n3, n4
    b0: int                                   # from check_masks betti (B4)
    b1: int

def trace(masks: Sequence[np.ndarray], offsets: Sequence[tuple[int, int]] | None = None) -> list[MaskRings]: ...
    # masks: 2-D arrays, any nonzero value = foreground, sides 0..4096. One engine call per trace call.
def polygonize(mask: np.ndarray, offset: tuple[int, int] = (0, 0)) -> list[Ring]: ...
```

`objdel.checker` exports `trace`, `polygonize`, `Ring`, `MaskRings`, `certify`, `certify_many`, and the topology
functions (`as_mask`, `betti`, `instance_topology`, `min_area`, `signature`, `size_filter`, `signatures`,
`topologies`). It no longer exports `bitquads`.

`objdel.checker.topology` adds `signatures(masks, min_area) -> list[tuple[int, ...]]` and
`topologies(masks) -> list[InstanceTopology]` (one engine call each).

### 3.2 Engine bindings (`gc`)

```text
check_masks(bits: uint8[N], hw: int32[n, 2] (H, W), start: int64[n + 1],
            *, chunk_pixels: int = 1 << 25, max_scratch_bytes: int = 1 << 30) -> dict   # keywords: tests only
    quads int64[n, 5]; ring_start int64[n + 1]; pts_start int64[R + 1]; ring_area int64[R]; pts int32[P, 2];
    ring_part int32[R]  (B3, per-mask rank, no offset);  betti int64[n, 2]  (b0, b1)

topology(bits: uint8[N], hw: int32[n, 2], start: int64[n + 1], min_area: int64[n], want_masks: bool = False) -> dict
    comp_start int64[n + 1]; comp_area int64[C]; comp_holes int32[C]      # after the size filter (min_area 0 = none)
    hole_start int64[C + 1]; hole_area int64[Hn]
    sig_start int64[n + 1]; sig int32[S]                                  # B17
    masks uint8[N]  (only with want_masks; filtered masks, offline callers only)

certify_polygons(xy: int32[V, 2], ring_start: int64, shape_start: int64, roles: int32 | None,
                 W: int32[n], H: int32[n], G: int,
                 ref_bits: uint8[N], ref_hw: int32[n_ref, 2], ref_start: int64[n_ref + 1], ref_of: int32[n],
                 min_area: int64[n],
                 pairs: int32[P, 2] | None = None, ref_pairs: int32[P, 2] | None = None,
                 pair_min_area: int64[P] | None = None) -> dict
    legal bool[n]; inter int64[n]; union int64[n]; equal bool[n];        # shape i vs reference ref_of[i]
    sig_start int64[n + 1], sig int32[S]  (candidate signatures, holes per component, descending);
    ref_sig_start int64[n + 1], ref_sig int32[S']
    pair_sig_start int64[P + 1], pair_sig int32[Q]; pair_ref_sig_start int64[P + 1], pair_ref_sig int32[Q']   # B11

eval_sets(xy, ring_start, shape_start, roles, W, H, G,        # polygons of n_sets * K sample objects, set-major
          tables: int64[n_sets]  (table handles: k rater planes + vote plane),
          K: int, legal: bool[n_sets * K]  (illegal sample -> empty drawn mask, as metrics' `drawn`),
          soft_cuts: int32[n_sets, T, 2]  (sample cut, rater cut per threshold),
          nsd_cut: int64[n_sets], *, max_scratch_bytes: int = 256 << 20) -> dict
    # M = K + k + 1 masks per set, in order: K drawn samples, k raters, the vote mask
    inter int64[n_sets, M, M]          # popcount(a & b); diagonal = areas
    inter_raw int64[n_sets, K, K]      # popcount between RAW sample masks (before the legal gate)
    area_raw int64[n_sets, K]          # raw sample areas
    soft int64[n_sets, T, 3]           # |P_t|, |R_t|, |P_t & R_t|; P_t = (#drawn samples at pixel >= sample cut),
                                       # R_t = (#raters at pixel >= rater cut)
    surf int64[n_sets, K + k, 3]       # surfel counts per mask by length class (sqrt(2)/2, 1, sqrt(2)) on the
                                       # (H+1) x (W+1) corner grid (2x2 code not in {0, 15})
    nsd_out int64[n_sets, K, k, 2, 3]  # per (sample, rater), direction (sample->rater, rater->sample), class:
                                       # surfels whose squared corner distance to the other surface > nsd_cut
                                       # (empty other surface: all surfels count)

add_tables_from_raters(raters: list[uint8[k_i, H_i, W_i]], t5: float64[m], min_area: int64[m]) -> dict   # B18
    handles int64[m]; perimeter int64[m]; rater_area_start int64[m + 1], rater_area int64[sum k];
    sig_start int64[sum k + 1], sig int32[S]; voter uint8[sum k]; n_vote int64[m]
add_tables_many(list of add_tables argument tuples) -> int64[m]        # one H2D copy and one wait per batch
add_tables(...), drop_tables(handle), resident()                      # unchanged signatures

edt(nonzero: uint8[H, W]) -> float32[H, W]   # B15

stats() -> dict    # waits, device_allocs, device_releases, host_allocs, host_releases, arena_bytes (B19)
reset_stats() -> None
trim(bytes: int) -> None                     # pool trim (tests, end of run); no-op on the CPU build
```

`score_polygons`, `score_masks`, `check_polygons` and `check_masks` accept `timing: bool = False`; when True the result
dict also has `timing` (per-stage milliseconds). Default calls add no synchronisation.

Host helper contracts (`src/objdel/eval/device.py`):
- the soft cut for n masks and threshold t = the smallest integer c in [0, n] with `float64(c) / n >= t`;
- `nsd_cut(tau)` = the largest integer q >= 0 with `sqrt(float64(q)) <= tau`;
- samples i and j of one item are duplicates exactly when their `format_ok` and `legal` flags are equal and, for every
  object name, `inter_raw[i, j] == area_raw[i] == area_raw[j]`.

### 3.3 Test markers and targets

- pytest markers: `gpu` (needs a CUDA device and the nvcc build; skipped where `nvidia-smi -L` lists no GPU, failing
  if `gc.build() != "cuda"` on a GPU host), `sanitizer` (run under compute-sanitizer; the tests live in
  `tests/test_gpu_sanitize.py` and use a small fixed corpus).
- `make gpu-check` (box only) = `make cuda` + `pytest -m gpu` + the sanitizer chain (memcheck, initcheck, synccheck,
  racecheck; each under `timeout 300`) + `experiments/phase2/perf/bench_checker.py` and the old-vs-new replay.
  No `-n` flag (pytest-xdist is not a dependency).
- `make gpu-check DRY=1` (Mac) prints the chain and runs `pytest --collect-only -q -m gpu`; it exits non-zero if any
  named test file is missing.
- No GitHub workflow runs the GPU chain in Phase D.

## 4. Constraints

C1. One source file for the CPU build and the GPU build (`cuda/gpucheck.cpp`). Keep `-ffp-contract=off` / `--fmad=false`.
C2. All cross-thread results are integers. Floats are made from integers by one thread in a fixed order, or on the host.
C3. Pixels never return to the host in `score_*`, `eval_sets`, `certify_polygons`, `add_tables_from_raters` and
    `topology` (unless `want_masks=True`, offline only). Only per-object, per-component or per-set numbers return.
    `check_masks` returns rings (not pixels). `edt` returns a map for the offline cache build only.
C4. The engine uses a lowest-priority non-blocking stream, and it adds no host synchronisation to `score_polygons`.
    Device memory on the CUDA build comes from a stream-ordered pool (no `cudaMalloc`/`cudaFree` after warm-up).
C5. Library code lives in `src/objdel` with tests in `tests/`. The numpy polygonizer, bit-quads, the cv2 topology
    functions, the `_kept` interval loop and the host table build live only in `tests/oracles/`. No `src/objdel` module
    imports `tests.oracles`.
C6. Old code is archived under `archive/` in the commit that replaces it. No dead copy stays beside new code.
C7. `make test` and `make lint` pass. Coverage stays >= 90 %.
C8. Toolchain: `cuda-toolkit[nvcc,cccl,nvtx]==13.0.3`, `nvidia-cuda-sanitizer-api==13.4.92`. Target sm_120 (RTX 5090).
    Compile check for sm_80/86/89/90/120 (`cuda/docker_check.sh`, NVTX ranges included).
C9. GPU spend for the phase <= the ceiling in D.md §8 ($1.40 including destroyed launches; $1.25 if the owner does not
    approve Q4) on Vast AI, North America hosts only, with a preflight before each rental.
C10. Python that reads `data/` runs with `-I`. RIGA data never leaves the Mac.
C11. `ExpertCache` construction and `load` make no engine call (fork safety).

## 5. Acceptance gates

| # | Gate | Pass rule |
|---|---|---|
| G1 | Polygonizer exactness | `trace` == oracle on the fuzz families (noise 0.05-0.95, blobs with holes, checkerboard, nested squares, comb, spiral, diagonals, stripes, 1xN, Nx1, sides 31/32/33, empty) and on the 37,093 crops of `data/checker_eval/masks.npz` (38,002 rings, 1,476,508 corners), on the CPU build and on the GPU |
| G2 | Part labels | B3 holds on the fuzz families and on the 37,093 crops |
| G3 | Engine equality | `pytest -m gpu` passes on an RTX 5090 with 0 mismatches on the corpus of 2,853 fuzz polygons, 6,024 fuzz masks, 144,133 LIDC polygons and 37,093 LIDC masks; B6 split test passes on the CPU build |
| G4 | Target hashes | B9: every listed `target_set_sha256` is byte-identical |
| G5 | Sanitizers | memcheck, initcheck, synccheck and racecheck report `ERROR SUMMARY: 0 errors` on `pytest -m "gpu and sanitizer"`, each tool inside 300 s |
| G6 | Determinism | 20 runs of `check_masks`, `score_polygons`, `eval_sets`, `certify_polygons`, `topology` and `add_tables_from_raters` on one fixed batch are byte-identical |
| G7 | Polygonizer speed | on one RTX 5090 rental, new vs pre-Phase-D `check_masks`, same inputs, alternating runs: outputs identical; new median <= 1.05 x old on every workload, and <= 0.9 x old on 256 frames of 512x512 and on the 37,093 crops |
| G8 | Reward path | same rental, replay set of 256 completions on their caches: new `score_polygons` median <= 1.05 x old; B19 waits == 3 |
| G9 | Eval equality | B14 on the saved val answers |
| G10 | Cache EDT and tables | B15 and B18 on 1,000 real caches |
| G11 | No CPU geometry | with cv2 connected-component functions, `scipy.ndimage.label` and `scipy.ndimage.distance_transform_edt` patched to raise, `build_target` (LIDC and RIGA fixture), `ExpertCache.build`/`load`/`device`, `certify` and `metrics()` run without error; no `src/objdel` import of `tests.oracles`; no `src/objdel/checker/polygonize.py` or `bitquad.py` |
| G12 | Repo rules | `make test`, `make lint` and `make cov` (>= 90 %) pass |
| G13 | Budget | C9, counted over every Phase D ledger row |
| G14 | Audit inside RL | 6 RL steps with the audit callback on: engine `arena_bytes` <= the `eval_sets` scratch cap + table bytes; no OOM; audit time and peak reserved memory recorded |
| G15 | No sync allocation | B19 on the CPU build; on the box, 0 sync allocations or frees over RL steps 2-6 |
| G16 | Certificate independence | B8 mutation case |
