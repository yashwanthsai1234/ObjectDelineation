// gpucheck.h -- batched exact checker: the parts of cpp/polycore the reward uses (even-odd and nonzero masks,
// validate totals) and the mask polygonizer of objdel.checker (bit-quads, rings by successor linking, pointer
// jumping and list ranking). One source, gpucheck.cpp: the C++ compiler builds it for the CPU, nvcc
// (OBJDEL_CUDA=ON, `make cuda`) for the GPU. Results are bit-identical to polycore / checker.polygonize.
//
// Inputs are packed like data/checker_eval (experiments/cuda_checker): many shapes or masks per call, one
// kernel launch per stage. Inputs must already be checked (caps, offsets): the binding does that with polycore's
// check_* functions.
#pragma once

#include <cstdint>
#include <vector>

namespace gpucheck {

using i32 = std::int32_t;
using i64 = std::int64_t;
using u8 = std::uint8_t;

struct Polygons {
  i64 n = 0;                          // shapes
  const i32* xy = nullptr;            // vertices x0, y0, x1, y1, ... of all rings
  const i64* ring_start = nullptr;    // ring r = vertices [ring_start[r], ring_start[r + 1])
  const i64* shape_start = nullptr;   // shape s = rings [shape_start[s], shape_start[s + 1]), at most 10
  const i32* roles = nullptr;         // per ring: -1 outer, p = hole of the shape's ring p (polycore roles)
  const i32* W = nullptr;             // per shape crop, 1 <= W, H <= 4096
  const i32* H = nullptr;
  i64 G = 1000;                       // grid
};

struct PolygonResult {
  std::vector<i64> mask_start;   // shape s = pixels [mask_start[s], mask_start[s + 1]), row-major (H, W)
  std::vector<u8> evenodd, nonzero;
  std::vector<i64> defects;      // (n, 13) = polycore validate(...)["total"], DEFECTS order
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

bool cuda_build();   // compiled by nvcc
int device_count();  // CUDA devices (0 on a CPU build or without a driver)
PolygonResult check_polygons(const Polygons& in);
MaskResult check_masks(const Masks& in);

}  // namespace gpucheck
