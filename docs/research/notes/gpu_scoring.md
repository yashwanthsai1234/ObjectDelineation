# Scoring on the GPU checker

Owner's rule: the checker stays on the GPU. Every check and every score runs in the CUDA checker: polygon -> mask ->
defects -> IoU_tc, F5, g_area, size-filtered topology signature. Masks are never copied to the host for scoring;
only per-object numbers come back. The CPU checker (polycore) leaves every runtime path.

## Device API (`cuda/gpucheck.h`, module `objdel._gpucheck`)

| Call | In | Out |
|---|---|---|
| `add_tables(raters, union, core, vote, dist, rater_area, voter, sig_start, sig, n_vote, perimeter, t5, min_area)` | one `ExpertCache` | handle |
| `drop_tables(handle)`, `resident()` | | |
| `score_polygons(xy, ring_start, shape_start, roles, W, H, G, tables)` | packed polygons, one table handle per object | per object: `defects` (13), `area`, `iou_tc`, `f5`, `g_area`, `topo_freq`, signature (`sig_start`, `sig`) |
| `score_masks(bits, hw, start, tables)` | masks that do not come from a completion (baselines, reward gates) | the same, defects 0 |

`check_polygons` (bit-packed masks) and `check_masks` (polygonizer) stay for callers that need pixels or rings
(below). Python: `objdel.reward.score_many / score_batch / score_completion` make one `score_polygons` call per
batch; `score_mask` wraps `score_masks`. The per-completion combination (mean over objects, joint topology over
RIGA's disc and cup, `sig_id`) is arithmetic on the returned per-object numbers.

## Expert tables resident on the device

`ExpertCache` is still built offline on the host (scipy EDT, `rater_tau`, `f5_tolerance`): data preparation, not
checking. `ExpertCache.device()` uploads its tables at the first score in a process and keeps them until the cache is
garbage-collected (`weakref.finalize` -> `drop_tables`; a forked child uploads its own copy and never frees the
parent's). Uploaded: k rater bit planes and union / core / vote planes (rows of 32-pixel words, LSB = leftmost, as
the raster), one float32 distance map (`dist_in_core` inside the union, `dist_out_union` outside it: each is 0 on
the other side, so one map holds both), rater areas, voter flags, rater signatures, and the scalars `n_vote`,
`perimeter`, `t5`, `min_area`. A 50 x 50 LIDC object is about 13 KB, a 100 x 100 one about 51 KB: all 7,341 LIDC
train items take about 140 MB. `add_tables` refuses a distance value that is neither 0 nor in [1, 2^24] (the
precondition of the exactness argument).

## One scoring call

1. Raster and validate (the existing kernels); the masks stay in the lane's arena, only the defects are copied.
2. Per row: first and last foreground column and count; per object: bounding box and area. Round trip 1 (5 numbers
   per object) gives the host the crop layout.
3. Signature, per object on its mask's bounding box padded by one background pixel, one item per crop row:
   label the foreground (8-connected), drop components below `min_area`; label the background (4-connected), fill
   holes below `min_area`; relabel both; count every hole for the component above its first raster pixel; gather each
   component's hole count and sort it descending. This is `objdel.checker.topology.signature` step by step
   (`size_filter`, then `instance_topology`).
4. Per frame row: intersections with each rater and the vote, spill (mask outside the union) and core-miss pixel
   counts and distance sums. Per object: reduce the rows, match the signature against the voting raters' signatures,
   and compute the terms. Round trip 2: the per-object numbers and signature lengths; round trip 3: the signatures.

Labeling is union-find with lock-free unions (Playne and Hawick, "A new algorithm for parallel connected-component
labelling on GPUs", IEEE TPDS 29(6), 2018): a root is always linked under a smaller index, so the final root of a
component is its smallest pixel index, its first raster pixel. Consequences used above: the outside is the
background component of pixel 0 (the padding); a hole's root is its first raster pixel, whose upper neighbour is
foreground of the enclosing component (an island inside the hole lies strictly below the hole's top row), exactly
the rule of `instance_topology`. The CPU build runs the same functors as plain loops.

Why not the ring polygonizer of `check_masks`: rings give b0, b1 and signed areas, but the size filter needs every
component's pixel area and every hole's own background area, and assigning a hole ring to its component needs the
nesting of rings. Labeling gives both directly, with one array per pixel of the crop.

## Exactness

Results are bit-identical to the host scorer this replaced (`tests/oracles/score_numpy.py`, a frozen copy of the old
`score_mask`):

- Areas and intersections are integer counts. IoU = inter / union is one float64 division of exactly representable
  integers, as Python's int / int (correctly rounded).
- Band depth D = (sum of spill distances - n_spill / 2) + (sum of miss distances - n_miss / 2). The distances are
  float32 values that are 0 or >= 1 (EDT of an integer grid, or the 1e6 sentinel), so each is a multiple of 2^-23.
  The device sums each as an exact integer count of 2^-23 (int64 per row, int128 over rows; order-free). While a
  sum is below 2^30, every partial sum of numpy's float64 pairwise summation is a multiple of 2^-23 below 2^30, hence
  exact (53 bits), so numpy's float64 result equals the exact sum, which is the device's. The rest is the same IEEE
  operations in the same order as the Python formulas: `1 / (1 + (2 D / P) / (t5 t5))`, the `min(1, .)` products of
  g_area, `len(same) / len(voters)`, Python's first-maximum rule for IoU_tc.
- No fused multiply-add: `-ffp-contract=off` (C++) and `--fmad=false` (nvcc) in `CMakeLists.txt`, so the CPU build
  and the GPU build round every operation once, as Python does. CUDA's float64 division is IEEE round-to-nearest.
- The signature is a topological invariant; the labeling, filter and hole rule follow `signature`'s definitions.

Where it cannot be exact: a band-depth sum of 2^30 px or more. With frames of at most 500 px the largest real sum is
500^2 x 707 < 2^28; only the 1e6 sentinel gets there (an empty rater union with more than 1,073 spill pixels, or a
core covering the whole frame). The device result is then the deterministic rounding of the exact sum, and numpy's
pairwise sum may differ in the last bits. Neither occurs in a prepared cache.

Verified: `tests/test_gpu_scoring.py` (hypothesis rater stacks and masks on frames of 1-70 px, edge cases at the size
filter, the 1e6 sentinel, polygon completions end to end; the GPU chain runs it on the GPU build before training);
`tests/test_reward_one_checker_call.py` (reward values recorded from the polycore + numpy code, reproduced exactly);
a stress run of 26,268 masks (frames up to 300 px) and 8,082 fuzz completions, 0 differences.

## What still reads pixels

- Exam metrics (`objdel.eval.metrics`): the per-sample terms come from the device; the item-level comparisons between
  samples and raters (pairwise IoU for GED and HM-IoU, soft Dice, NSD, duplicate hashes) still read the samples'
  masks from one `to_masks` call per item. NSD sums irrational surfel lengths in float64, so a device port must
  either replicate numpy's pairwise order or change the metric's last bits; not done yet.
- Building SFT targets and probe shapes (`objdel.targets`, `objdel.eval.probes`) and stacking slices into volumes
  (`objdel.checker.volume`): data generation, not scoring; they take the GPU checker's masks (`to_masks`).

## The CPU process pool

Gone from the reward: `score_batch` has no executor, `TRLReward` no pool. `serve.fork_pool` remains only for the
exam metrics in `select_ckpt` and `g_probe` (host work per item); each worker uploads the caches it scores to its own
device context.

## polycore, archived

polycore (`cpp/polycore`, `src/objdel/_polycore.cpp`, `objdel.checker.polycore`, its tests and
`experiments/cuda_checker/prepare_data.py`) is in `archive/cpu_checker/`, with a README. It is gone from the build
(`CMakeLists.txt`, `Makefile` `EXT_SRC` and `fuzz`, the cibuildwheel test command, CI) and from every runtime
import. Lab prototypes under `experiments/*/lab/` keep their own vendored copies and are not run.

- `gpucheck.h` now holds the caps, the defect enum and names, and the input checks that came from `polycore.hpp`.
- `objdel.geometry`:
  - `to_mask` is `to_masks` of one object, so targets and probes are built from the GPU checker's masks.
  - `pack_rings` packs rings, with the same int32 checks as before.
  - `signed_area2` is exact Python integer arithmetic.
- Tests:
  - Before the archive, polycore's results were recorded for every case the tests and `eval.py` compare with:
    11,739 cases (defects and mask SHA-256, keyed per case) in `tests/goldens/polycore.npz`, written by
    `archive/cpu_checker/record_goldens.py`. No test imports polycore.
  - The K1/K2 oracle tests now check the CUDA checker directly against the slow exact Python oracles
    (`tests/test_checker_oracle.py`).
  - `eval.py` can no longer shrink a polygon mismatch, because the sub-cases have no reference. Mask mismatches
    still shrink.
