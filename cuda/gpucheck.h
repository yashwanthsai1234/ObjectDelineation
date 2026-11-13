// gpucheck.h -- batched exact checker: polygon legality (13 defect classes) and the even-odd mask, the mask polygonizer
// of objdel.checker (bit-quads, rings by run union-find, closed-form ring starts and seeded crack tracing: the K1-draft
// design), and the reward's scorer (score_polygons / score_masks: IoU_tc, F5, g_area and the size-filtered signature
// against expert tables resident on the device; only per-object numbers leave it). One source, gpucheck.cpp: the C++
// compiler builds it for the CPU, nvcc (OBJDEL_CUDA=ON, `make cuda`) for the GPU. Results are bit-identical to the CPU
// checker it replaced (polycore, archive/cpu_checker), checker.polygonize and the host scorer
// (docs/research/notes/gpu_scoring.md).
//
// Inputs are packed like data/checker_eval: many shapes or masks per call, one
// kernel launch per stage. Inputs must already be checked (caps, offsets): the binding does that with the check_*
// functions below. Calls are serialized by one mutex (the device memory, staging buffers and streams are shared and
// grow-only), so any thread may call them.
#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace gpucheck {

using i32 = std::int32_t;
using i64 = std::int64_t;
using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;

// Input caps. A shape is a list of rings of integer vertices on a G x G grid over a W x H crop; grid point (x, y) is
// pixel point (x W / G, y H / G). Pixel (c, r) is inside iff its centre is inside under the half-open crossing rule
// (centres on left / top edges inside, on right / bottom edges outside). Outer rings are clockwise on screen (positive
// shoelace sum x_i y_{i+1} - x_{i+1} y_i, y down), holes counter-clockwise.
constexpr i32 kMaxDim = 4096;
constexpr i32 kMaxRings = 10;
constexpr i32 kMaxRingVertices = 1000;
constexpr i64 kMaxGrid = i64(1) << 30;

// Defect classes of `defects`; crossing / touching / overlap classes count edge pairs, ring-ring classes once per pair.
enum Defect : int {
  kFewerThan3Points = 0,  // < 3 distinct vertices after dropping repeats (no other self label then)
  kZeroArea,              // shoelace sum == 0
  kSelfCrossing,          // # non-adjacent edge pairs of the ring that cross properly
  kVertexPinch,           // # non-adjacent pairs touching at a point (vertex=vertex or vertex on edge)
  kCollinearOverlap,      // # non-adjacent pairs overlapping along a segment of positive length
  kSpike,                 // # adjacent pairs folding back on each other (collinear backtrack)
  kRingRingCrossing,      // # edge pairs (this ring x other rings) crossing properly
  kRingRingTouching,      // # edge pairs (this ring x other rings) touching / overlapping
  kHoleOutsideOuter,      // hole with a witness point strictly outside its outer ring (a hole that crosses its outer
                          // gets this flag and ring_ring_crossing)
  kNestedHoles,           // # sibling holes this hole overlaps / is nested with
  kNestedOuters,          // # other outer parts this outer overlaps / is nested with
  kWrongWinding,          // outer with net area < 0 or hole with net area > 0 (net shoelace area of the cleaned ring,
                          // so a bow-tie is judged by its larger lobe)
  kOutOfRange,            // # vertices with x or y outside [0, G]
  kNumDefects
};

inline constexpr const char* kDefectNames[kNumDefects] = {
    "fewer_than_3_points", "zero_area",          "self_crossing",
    "vertex_pinch",        "collinear_overlap",  "spike",
    "ring_ring_crossing",  "ring_ring_touching", "hole_outside_outer",
    "nested_holes",        "nested_outers",      "wrong_winding",
    "out_of_range"};

// roles[r] == -1: ring r is an outer ring; roles[r] == p: ring r is a hole of outer ring p.
inline void check_roles(const i32* roles, i32 nrings) {
  for (i32 r = 0; r < nrings; ++r) {
    const i32 p = roles[r];
    if (p != -1 && (p < 0 || p >= nrings || p == r || roles[p] != -1))
      throw std::invalid_argument("roles[r] must be -1 (outer) or the index of an outer ring (hole)");
  }
}

inline void check_grid(i64 G) {
  if (G < 1 || G > kMaxGrid) throw std::invalid_argument("G must be in [1, 2^30]");
}

inline void check_dims(i64 W, i64 H) {
  if (W < 1 || H < 1 || W > kMaxDim || H > kMaxDim) throw std::invalid_argument("W and H must be in [1, 4096]");
}

struct Polygons {
  i64 n = 0;                          // shapes
  const i32* xy = nullptr;            // vertices x0, y0, x1, y1, ... of all rings
  const i64* ring_start = nullptr;    // ring r = vertices [ring_start[r], ring_start[r + 1])
  const i64* shape_start = nullptr;   // shape s = rings [shape_start[s], shape_start[s + 1]), at most 10
  const i32* roles = nullptr;         // per ring: -1 outer, p = hole of the shape's ring p (check_roles)
  const i32* W = nullptr;             // per shape crop, 1 <= W, H <= 4096
  const i32* H = nullptr;
  i64 G = 1000;                       // grid
};

struct PolygonResult {
  std::vector<i64> word_start;   // shape s = words [word_start[s], word_start[s + 1]) = H rows of ceil(W / 32) words
  std::vector<u32> evenodd;      // even-odd mask, bit-packed per row: pixel x is bit x % 32 (LSB first) of word x / 32;
                                 // the bits past W in a row's last word are 0
  std::vector<i64> defects;      // (n, kNumDefects) defect counts per shape, in Defect order
};

struct Masks {
  i64 n = 0;
  const u8* bits = nullptr;    // all masks row-major; nonzero = foreground
  const i64* start = nullptr;  // mask m = bits [start[m], start[m + 1]) = H * W
  const i32* H = nullptr;      // 0 <= W, H <= 4096
  const i32* W = nullptr;
};

struct MaskResult {
  std::vector<i64> quads;       // (n, 5): n1, n2, nD, n3, n4 (objdel.checker.bitquad)
  std::vector<i64> ring_start;  // mask m = rings [ring_start[m], ring_start[m + 1]), in polygonize's order
  std::vector<i64> pts_start;   // ring r = corners [pts_start[r], pts_start[r + 1])
  std::vector<i64> ring_area;   // signed area per ring: > 0 outer, < 0 hole
  std::vector<i32> pts;         // corners x0, y0, x1, y1, ...
  std::vector<i32> ring_part;   // per ring: 0-based raster rank (by first pixel) in its mask of the foreground
                                // 8-component it bounds (an outer ring's own, a hole's surrounding one)
  std::vector<i64> betti;       // (n, 2): foreground components, bounded background components (from the same run
                                // labels as the rings, so never an independent check of them)
};

// check_masks' work split. The defaults are the production values; tests set others (the binding's test keywords).
constexpr i64 kMaxChunkPixels = i64(1) << 30;
struct MaskOptions {
  i32 seed_rows = 16;                     // >= 0: every K-th row's run sides also start a trace (0 = one per ring)
  i64 chunk_pixels = i64(1) << 25;        // 1..2^30: masks per chunk by (W + 1)(H + 1); a larger mask goes alone
  i64 max_scratch_bytes = i64(1) << 30;   // >= 1: a chunk whose run scratch is larger is halved; one such mask raises
};

// The expert tables of one object of one item (objdel.reward.ExpertCache), as host arrays. add_tables copies them to
// the device once; they stay resident until drop_tables. The binding checks shapes and the dist precondition.
struct Tables {
  i32 H = 0, W = 0, k = 0;          // frame (1 <= W, H <= 4096) and raters (2 <= k <= 32)
  const u8* raters = nullptr;       // (k, H, W), nonzero = inside
  const u8* uni = nullptr;          // (H, W) union, core (intersection) and strict-majority vote of the raters
  const u8* core = nullptr;
  const u8* vote = nullptr;
  const float* dist = nullptr;      // (H, W) distance to the complement of the core inside the union, to the union
                                    // outside it (ExpertCache.dist_in_core, dist_out_union); every value 0 or in [1, 2^24]
  const i64* rater_area = nullptr;  // (k)
  const u8* voter = nullptr;        // (k) 1 = the rater votes on topology (ExpertCache.voters)
  const i64* sig_start = nullptr;   // (k + 1): rater r's signature = sig[sig_start[r], sig_start[r + 1])
  const i32* sig = nullptr;         // holes per component, descending (ExpertCache.rater_sigs)
  i64 n_vote = 0, perimeter = 1;    // vote area, contour length of the union
  double t5 = 1;                    // F5 tolerance
  i32 min_area = 0;                 // the signature's size filter
};

// The rater planes of one cache (objdel.reward.ExpertCache) for add_tables_from_raters, which builds its Tables on
// the device. The binding checks the caps (1 <= W, H <= 4096, 2 <= k <= 32, t5 > 0, min_area >= 0).
struct Raters {
  i32 H = 0, W = 0, k = 0;
  const u8* px = nullptr;  // (k, H, W), nonzero = inside
  double t5 = 1;
  i32 min_area = 0;
};

struct BuiltTables {                     // per cache of an add_tables_from_raters call
  std::vector<i64> handles, perimeter, n_vote;
  std::vector<i64> rater_area_start;     // (m + 1): cache c's raters = [rater_area_start[c], rater_area_start[c + 1])
  std::vector<i64> rater_area;           // per rater of all caches
  std::vector<u8> voter;                 // per rater: 1 = votes on topology
  std::vector<i64> sig_start;            // (raters + 1): rater r's signature = sig[sig_start[r], sig_start[r + 1])
  std::vector<i32> sig;                  // holes per component, descending
};

struct Scores {                // per object (objdel.reward.score: ObjectTerms), nothing per pixel
  std::vector<i64> defects;    // (n, 13) validate totals (score_polygons; zero for score_masks)
  std::vector<i64> area;       // pixels of the mask
  std::vector<i64> sum_x, sum_y;  // sums of the column and row indices of its pixels (the centroid times the area)
  std::vector<double> iou_tc, f5, g_area, topo_freq;
  std::vector<i64> sig_start;  // (n + 1): object s's signature = sig[sig_start[s], sig_start[s + 1])
  std::vector<i32> sig;        // holes per component of the size-filtered mask, descending
};

// The ladder certificate of objdel.targets (certify_polygons): shape s of a polygon call against reference plane
// ref_of[s] of `refs` (same frame), both size-filtered at min_area[s]. Pair p (the RIGA rim) is shapes
// (shapes[2p], shapes[2p + 1]) and reference planes (refs[2p], refs[2p + 1]): the signatures of a & ~b at min_area[p].
struct Pairs {
  i64 n = 0;
  const i32* shapes = nullptr;
  const i32* refs = nullptr;
  const i64* min_area = nullptr;
};

struct Certified {                 // an illegal shape draws nothing (empty raster)
  std::vector<u8> legal, equal;    // per shape: no defect; raster == reference
  std::vector<i64> inter, uni;     // per shape: popcounts of raster & reference, raster | reference
  std::vector<i64> sig_start, ref_sig_start, pair_sig_start, pair_ref_sig_start;  // (n + 1) or (P + 1)
  std::vector<i32> sig, ref_sig, pair_sig, pair_ref_sig;  // holes per component, descending
};

// objdel.checker.topology of many masks (topology): each mask size-filtered at min_area[m] (0 = none), then its
// 8-components in raster order of their first pixel and the bounded 4-background components (holes).
struct Topology {
  std::vector<i64> comp_start, comp_area;  // (n + 1): mask m = components [comp_start[m], comp_start[m + 1]); areas
  std::vector<i32> comp_holes;             // holes per component
  std::vector<i64> hole_start, hole_area;  // (C + 1): component c's hole areas, in raster order of their first pixels
  std::vector<i64> sig_start;              // (n + 1)
  std::vector<i32> sig;                    // holes per component, descending (topology.signature)
  std::vector<u8> masks;                   // want_masks only: the filtered masks, laid out as the input (0 / 1)
};

// eval_sets (objdel.eval.device): set s = shapes s K .. s K + K - 1 of the polygons (the K samples of one object of
// one item; an unparsed sample is a shape with no ring), against the tables tables[s] (frame = the shapes' frame;
// every set's tables have the same k raters, else std::invalid_argument). M = K + k + 1 masks per set, in order: the K
// drawn samples (raster AND legal), the k raters, the vote. Inputs checked by the binding (1 <= K <= 256,
// 1 <= T <= kMaxThresholds, sample cuts in [0, K], rater cuts >= 0, nsd_cut >= 0, 1 <= band_d <= 31,
// max_scratch_bytes >= 1); here: rater cuts <= k, and a set whose scratch is above the cap raises
// std::invalid_argument. Results do not depend on the chunk split.
constexpr i32 kMaxThresholds = 16;
struct EvalSets {
  i64 n = 0;                    // sets
  i32 K = 0, T = 0;             // samples, soft-Dice thresholds
  const i64* tables = nullptr;  // (n) table handles
  const u8* legal = nullptr;    // (n K) 0: the sample draws nothing (metrics' `drawn`)
  const i32* cuts = nullptr;    // (n, T, 2) per threshold: the sample cut, the rater cut
  const i64* nsd_cut = nullptr; // (n) the largest squared corner distance within the NSD tolerance
  const i32* band_d = nullptr;  // (n) boundary-IoU band width, px
  i64 max_scratch_bytes = i64(256) << 20;
};
struct EvalCounts {                  // per set, row-major
  i32 k = 0;                         // raters (of every set's tables)
  std::vector<i64> inter;            // (n, M, M) popcount(a & b) of the drawn samples, raters and vote
  std::vector<i64> inter_raw;        // (n, K, K) the same between the raw sample rasters (before the legal gate)
  std::vector<i64> band_inter;       // (n, M, M) the same between their boundary-IoU bands (metrics.band)
  std::vector<i64> soft;             // (n, T, 3) |P_t|, |R_t|, |P_t & R_t|: P_t = drawn samples >= cut, R_t = raters >= cut
  std::vector<i64> surf;             // (n, K + k, 3) surfels of the drawn samples and raters by length class
                                     // (sqrt(2)/2, 1, sqrt(2)) on the (H + 1) x (W + 1) corner grid
  std::vector<i64> nsd;              // (n, K, k, 2, 3) per (sample, rater), direction (sample -> rater, rater ->
                                     // sample), class: surfels whose squared distance to the other surface > nsd_cut
                                     // (all of them when the other surface is empty)
  std::vector<double> vote_iou_tc;   // (n) IoU_tc of the pixel vote of the legal drawn samples (2 x count > legal)
};
EvalCounts eval_sets(const Polygons& in, const EvalSets& sets);
// Per-stage device time of one call in milliseconds, stages in the order they first ran (timing=True in the binding).
// Read after the call's last wait: asking for it adds no host wait.
using Times = std::vector<std::pair<std::string, double>>;

// Engine counters, process-wide (gpucheck.cpp, "counters"). Events since the last reset_stats():
struct Stats {
  i64 waits = 0;        // host waits for a lane (stream synchronizations inside a call)
  i64 drain_syncs = 0;  // stream synchronizations of both lanes on every exit of a call (Drain)
  i64 label_pairs = 0;  // crop rows the scorer's second label pair relabels
  i64 sync_alloc = 0, sync_free = 0;  // allocations that synchronize the device: cudaMalloc(Host) / cudaFree(Host);
                                      // on the CPU build, the host staging buffers
  i64 pool_allocs = 0, pool_frees = 0;  // stream-ordered pool allocations; on the CPU build, device-kind memory
  i64 pinned_grows = 0;  // growths of a staging buffer after its one sizing in the first call (logged to stderr)
  i64 eval_scratch_peak = 0;  // the largest scratch block of one eval_sets chunk (bytes)
  // State, kept by reset_stats():
  i64 arena_bytes = 0;        // capacity of the grow-only device arenas of both lanes
  i64 table_bytes = 0;        // bytes of the resident expert tables
  i64 live_device_bytes = 0;  // device bytes allocated and not released (== arena_bytes + table_bytes)
  i64 pool_reserved = 0;          // bytes the engine's pool holds from the driver (CUDA build, read live)
  i64 default_pool_reserved = 0;  // bytes the device's default pool holds (torch, vLLM; the engine never uses it)
};

bool cuda_build();   // compiled by nvcc
int device_count();  // CUDA devices (0 on a CPU build or without a driver)
Stats stats();
void reset_stats();  // zeroes the events, keeps the state
// times != nullptr: also the per-stage times of the call.
PolygonResult check_polygons(const Polygons& in, Times* times = nullptr);
MaskResult check_masks(const Masks& in, Times* times = nullptr, const MaskOptions& opt = {});
std::vector<std::string> test_hooks();  // the test-only compile definitions of this build (none in a shipped one)
i64 add_tables(const Tables& t);  // -> handle
std::vector<i64> add_tables_many(const std::vector<Tables>& ts);  // one host wait per batch (none when empty)
// Tables built on the device from rater planes (one upload, one host wait per batch; none when empty): the same
// tables add_tables stores from the host-built arrays, and the scalars the host needs.
BuiltTables add_tables_from_raters(const std::vector<Raters>& in);
// The exact Euclidean distance (float32(sqrt(float64(d2))), d2 the integer squared distance) from each pixel of an
// H x W map (row-major, 0 <= H, W <= kMaxEdtDim; checked by the binding) to its nearest zero pixel, 1e6 everywhere
// when there is none (objdel.reward.cache.edt, the offline cache build; the only entry that returns pixels).
// kMaxEdtDim: the largest squared distance 2 x 32766^2 stays below the int32 no-source mark; the cache build pads a
// cropped 4096 frame by 1 px on each side, so the frame bound kMaxDim is too small here.
constexpr i32 kMaxEdtDim = 32767;
std::vector<float> edt(const u8* nonzero, i32 H, i32 W);
void drop_tables(i64 handle);     // no host wait
void trim(i64 bytes);             // the engine's pool keeps <= `bytes` unused (cudaMemPoolTrimTo); CPU build: no-op
i64 resident();                    // tables on the device
// Object s is scored against the tables tables[s]; its frame must be theirs (std::invalid_argument otherwise).
Scores score_polygons(const Polygons& in, const i64* tables, Times* times = nullptr);
Scores score_masks(const Masks& in, const i64* tables, Times* times = nullptr);
// Inputs checked by the binding: ref_of[s] < refs.n with the shape's frame, min_area >= 0, pair indices in range
// and each pair's two planes of one frame.
Certified certify_polygons(const Polygons& in, const Masks& refs, const i32* ref_of, const i64* min_area,
                           const Pairs& pairs = {});
// min_area >= 0 per mask (checked by the binding); std::invalid_argument beyond 2^31 - 1 crop pixels in one call.
Topology topology(const Masks& in, const i64* min_area, bool want_masks = false);

}  // namespace gpucheck
