// k1poly.h -- K1: a work-optimal mask polygonizer (binary masks -> exact rectilinear rings), single source for the
// CPU (C++ compiler) and the GPU (nvcc -x cu). Output is identical, field by field, to gpucheck::check_masks
// (REPO/cuda/gpucheck.cpp) and to objdel.checker.polygonize (REPO/src/objdel/checker/polygonize.py):
//   rings per mask in polygonize's order (by the ring's first east edge in raster order), each ring's corners from its
//   canonical start (minimum (y (W + 1) + x) 4 + d edge key), signed area (> 0 outer ring, < 0 hole), plus the
//   bit-quad counts n1, n2, nD, n3, n4 of objdel.checker.bitquad.
// Algorithm (DESIGN.md): bit-pack -> runs per row -> union-find on runs (foreground 8-connected, background
// 4-connected) -> one ring per foreground component (outer) and per bounded background component (hole), its start
// and order key in closed form from the component labels -> rank rings by a popcount scan over a bitmap (no sort) ->
// crack-edge tracing from seeds (ring starts plus, optionally, the vertical edges on every K-th row) that visits only
// boundary edges -> corner offsets by a scan -> corners written in place.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace k1poly {

using i32 = std::int32_t;
using i64 = std::int64_t;
using u8 = std::uint8_t;
using u32 = std::uint32_t;

constexpr i32 kMaxDim = 4096;

struct Masks {                 // the layout of gpucheck::Masks
  i64 n = 0;
  const u8* bits = nullptr;    // all masks row-major; nonzero = foreground
  const i64* start = nullptr;  // mask m = bits [start[m], start[m + 1]) = H * W
  const i32* H = nullptr;      // 0 <= W, H <= 4096
  const i32* W = nullptr;
};

struct Result {                 // the layout of gpucheck::MaskResult
  std::vector<i64> quads;       // (n, 5): n1, n2, nD, n3, n4
  std::vector<i64> ring_start;  // mask m = rings [ring_start[m], ring_start[m + 1]), in polygonize's order
  std::vector<i64> pts_start;   // ring r = corners [pts_start[r], pts_start[r + 1])
  std::vector<i64> ring_area;   // signed area per ring: > 0 outer, < 0 hole
  std::vector<i32> pts;         // corners x0, y0, x1, y1, ...
};

struct Options {
  // Seed rows: every boundary edge that is a vertical side of a foreground run on a row y with y % seed_rows == 0
  // also starts a trace segment, so one long ring is traced by many threads (span ~ seed_rows * W steps instead of
  // the ring length). 0 = one trace per ring.
  i32 seed_rows = 16;
  // Pixels per device chunk (masks are processed in chunks of at most this many pixels; a larger mask goes alone).
  i64 chunk_pixels = i64(1) << 25;
  // Device scratch sized by the run counts (layout B) above this splits the chunk in halves (one mask always runs).
  // Typical masks need little; the worst case (checkerboard: a run per pixel pair) needs ~150 bytes per pixel.
  i64 max_scratch_bytes = i64(1) << 30;
};

struct Timing {  // optional per-call stage times in milliseconds (GPU build with events; CPU build wall clock)
  double h2d = 0, ccl = 0, rings = 0, trace = 0, d2h = 0, total = 0;
};

bool cuda_build();
std::string device_info();  // GPU name, SMs, L2, CC (CUDA build); "cpu" otherwise
Result polygonize_masks(const Masks& in, const Options& opt = {}, Timing* timing = nullptr);

}  // namespace k1poly
