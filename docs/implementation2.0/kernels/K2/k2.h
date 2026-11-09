// k2.h -- K2 exact scoring / evaluation kernels (k2score.inc) behind a small host API, for the test harness and the
// micro-benchmark. One source: k2.cpp is compiled by the C++ compiler (CPU build) or by nvcc -x cu (GPU build); both
// builds return the same bits. Masks use gpucheck's packed layout (gpucheck.h PolygonResult::evenodd): mask m = H[m]
// rows of ceil(W[m] / 32) uint32 words at word_start[m], pixel x = bit x % 32 of word x / 32, bits past W are 0.
#pragma once

#include <cstdint>
#include <vector>

#include "vendor/gpucheck.h"

namespace k2 {

using gpucheck::i32;
using gpucheck::i64;
using gpucheck::u32;
using gpucheck::u8;
using u64 = std::uint64_t;

constexpr i32 kNoSource = 0x7FFFFFFF;        // squared distance when the mask has no source pixel
constexpr i64 kMaxHmStates = i64(1) << 16;  // HM-IoU quota states per set on the device

struct Bits {
  i64 n = 0;
  const u32* words = nullptr;
  const i64* word_start = nullptr;  // n + 1
  const i32* H = nullptr;
  const i32* W = nullptr;
};

struct Times {  // device milliseconds per stage (GPU: CUDA events; CPU build: steady clock)
  int n = 0;
  const char* name[16] = {};
  double ms[16] = {};
  void add(const char* s, double v) {
    if (n < 16) name[n] = s, ms[n] = v, ++n;
  }
  double total() const {
    double t = 0;
    for (int i = 0; i < n; ++i) t += ms[i];
    return t;
  }
};

bool cuda_build();

// ---- validity: gpucheck's path (EdgePairs + ShapeTotals, one thread per shape) against K2's (uniform-grid pairs +
// one warp per shape). Shapes with fewer than grid_min_edges input vertices, or whose grid would hold more than
// budget_per_edge entries per input vertex, keep EdgePairs.
struct Validity {
  std::vector<i64> defects_old, defects_new;  // (n, 13)
  std::vector<i32> counts_old, counts_new;    // (rings, 13), after the totals step
  std::vector<i32> pairs_old, pairs_new;      // (n, 200): ring-ring cross / touch counts
  std::vector<u8> clean;                      // new path: every total 0
  std::vector<i32> bad_rings;                 // new path: rings with a nonzero count
  std::vector<u8> used_grid;                  // per shape
  i64 grid_entries = 0;
  Times t_old, t_new;
};
Validity validity(const gpucheck::Polygons& in, i64 grid_min_edges, i64 budget_per_edge);

// ---- boundary and EDT
std::vector<u32> boundary(const Bits& m, std::vector<i64>* count = nullptr);  // same layout as the input
std::vector<i32> edt_sq(const Bits& src, Times* t = nullptr);  // per pixel (masks in order, row-major): squared
                                                               // distance to the nearest set pixel, kNoSource if none

// ---- boundary F and HD95: pair p compares masks a[p] and b[p] (same frame). Direction A = boundary pixels of a
// against the boundary of b; direction B the reverse. All outputs are integers.
struct PairQuery {
  i64 n = 0;
  const i64* a = nullptr;
  const i64* b = nullptr;
  i64 tau2 = 4;               // boundary F tolerance, squared px (d^2 <= tau2)
  i64 p_num = 95, p_den = 100;  // percentile p_num / p_den
};
struct Surface {
  std::vector<i64> n_a, n_b;          // boundary pixels
  std::vector<i64> hit_a, hit_b;      // boundary pixels within tau of the other boundary
  std::vector<i64> rank_a, rank_b;    // k = floor(p (n - 1)), -1 when n = 0
  std::vector<i64> lo_a, hi_a, lo_b, hi_b;  // squared distances at ranks k and min(k + 1, n - 1)
  Times t;
};
Surface surface(const Bits& m, const PairQuery& q);
double boundary_f(i64 n_a, i64 n_b, i64 hit_a, i64 hit_b);  // host formulas from the integers
double hd_percentile(i64 n_a, i64 rank_a, i64 lo_a, i64 hi_a, i64 n_b, i64 rank_b, i64 lo_b, i64 hi_b,
                     i64 p_num, i64 p_den);

// ---- sets: set s = masks [mask_start[s], mask_start[s + 1]) of m, the first K[s] samples then k >= 1 raters, one frame
struct SetQuery {
  i64 n = 0;
  const i64* mask_start = nullptr;
  const i32* K = nullptr;
};
struct SetStats {
  std::vector<i64> inter_start, inter;   // M x M intersection counts per set (diagonal = areas)
  std::vector<i64> choice_start;
  std::vector<i32> choice;               // HM matching: rater of each sample (-1 if the state space is too large)
  std::vector<u64> hm_q;                 // optimal sum of floor(2^32 IoU)
  std::vector<double> hm_iou, ged;
  Times t;
};
SetStats sets(const Bits& m, const SetQuery& q);

}  // namespace k2
