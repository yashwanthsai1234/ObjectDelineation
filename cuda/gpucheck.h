// gpucheck.h -- batched exact checker: polygon legality (13 defect classes) and the even-odd mask, the mask polygonizer of objdel.checker (bit-quads, rings by successor linking, pointer jumping and list ranking), and
// the reward's scorer (score_polygons / score_masks: IoU_tc, F5, g_area and the size-filtered signature against
// expert tables resident on the device; only per-object numbers leave it). One source, gpucheck.cpp: the C++ compiler
// builds it for the CPU, nvcc (OBJDEL_CUDA=ON, `make cuda`) for the GPU. Results are bit-identical to the CPU checker
// it replaced (polycore, archive/cpu_checker), checker.polygonize and the host scorer
// (docs/research/notes/gpu_scoring.md).
//
// Inputs are packed like data/checker_eval (experiments/cuda_checker): many shapes or masks per call, one
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

struct Scores {                // per object (objdel.reward.score: ObjectTerms), nothing per pixel
  std::vector<i64> defects;    // (n, 13) validate totals (score_polygons; zero for score_masks)
  std::vector<i64> area;       // pixels of the mask
  std::vector<double> iou_tc, f5, g_area, topo_freq;
  std::vector<i64> sig_start;  // (n + 1): object s's signature = sig[sig_start[s], sig_start[s + 1])
  std::vector<i32> sig;        // holes per component of the size-filtered mask, descending
};

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
  // State, kept by reset_stats():
  i64 arena_bytes = 0;        // capacity of the grow-only device arenas of both lanes
  i64 table_bytes = 0;        // bytes of the resident expert tables
  i64 live_device_bytes = 0;  // device bytes allocated and not released (== arena_bytes + table_bytes)
  i64 contexts = 0;           // 1 once a call has opened a CUDA context in this process (0 on the CPU build)
};

bool cuda_build();   // compiled by nvcc
int device_count();  // CUDA devices (0 on a CPU build or without a driver)
Stats stats();
void reset_stats();  // zeroes the events, keeps the state
// times != nullptr: also the per-stage times of the call.
PolygonResult check_polygons(const Polygons& in, Times* times = nullptr);
MaskResult check_masks(const Masks& in, Times* times = nullptr);
i64 add_tables(const Tables& t);  // -> handle
void drop_tables(i64 handle);
i64 resident();                    // tables on the device
// Object s is scored against the tables tables[s]; its frame must be theirs (std::invalid_argument otherwise).
Scores score_polygons(const Polygons& in, const i64* tables, Times* times = nullptr);
Scores score_masks(const Masks& in, const i64* tables, Times* times = nullptr);

}  // namespace gpucheck
