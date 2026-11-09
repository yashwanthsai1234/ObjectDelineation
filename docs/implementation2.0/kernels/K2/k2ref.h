// k2ref.h -- independent CPU references for the K2 kernels: brute force on unpacked bool masks, written without
// k2score.inc, so a shared bug cannot hide. The validity reference is the existing CPU path itself (gpucheck's
// EdgePairs + ShapeTotals, and gpucheck::check_polygons), compared inside k2test.cpp.
#pragma once

#include <vector>

#include "k2.h"

namespace k2ref {

using k2::i32;
using k2::i64;
using k2::u32;
using k2::u64;

struct Mask {  // unpacked, row-major
  i32 H = 0, W = 0;
  std::vector<unsigned char> px;
  bool at(i32 x, i32 y) const { return x >= 0 && y >= 0 && x < W && y < H && px[size_t(y) * W + x] != 0; }
};

std::vector<Mask> unpack(const k2::Bits& b);
void pack(const std::vector<Mask>& ms, std::vector<u32>& words, std::vector<i64>& word_start, std::vector<i32>& H,
          std::vector<i32>& W);

Mask boundary(const Mask& m);                 // foreground with a 4-neighbour in the background (outside = background)
std::vector<i32> edt_brute(const Mask& src);  // min over set pixels of dx^2 + dy^2; k2::kNoSource if none

struct Directed {
  i64 n = 0, hit = 0, rank = -1, lo = -1, hi = -1;
};
Directed directed(const Mask& from, const Mask& to, i64 tau2, i64 p_num, i64 p_den);  // sort-based

struct Set {
  std::vector<i64> inter;  // M x M
  std::vector<i32> choice;
  u64 hm_q = 0;
  double hm_iou = 0, ged = 0;
};
Set set_stats(const std::vector<Mask>& masks, i32 K);  // brute-force enumeration of all assignments (small K only)

}  // namespace k2ref
