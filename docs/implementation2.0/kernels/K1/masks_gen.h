// masks_gen.h -- random and adversarial binary masks for the K1 tests and benchmark, plus a batch packer.
#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "k1poly.h"

namespace k1gen {

using k1poly::i32;
using k1poly::i64;
using k1poly::u8;

struct Mask {
  i32 H = 0, W = 0;
  std::vector<u8> px;  // row-major, 0/1
  std::string name;
  u8& at(i32 x, i32 y) { return px[size_t(y) * W + x]; }
};

inline Mask blank(i32 H, i32 W, const std::string& name) {
  Mask m;
  m.H = H, m.W = W, m.name = name;
  m.px.assign(size_t(H) * W, 0);
  return m;
}

inline Mask noise(std::mt19937_64& rng, i32 H, i32 W, double p) {
  Mask m = blank(H, W, "noise p=" + std::to_string(p));
  std::bernoulli_distribution b(p);
  for (auto& v : m.px) v = b(rng);
  return m;
}

// Union of random discs and rectangles, minus some smaller discs: blob-like shapes with holes and islands (a rough
// stand-in for lesion masks; the real benchmark input is the LIDC targets, GPU_RUN.md step 5).
inline Mask blobs(std::mt19937_64& rng, i32 H, i32 W, int k, double rmax) {
  Mask m = blank(H, W, "blobs k=" + std::to_string(k));
  std::uniform_real_distribution<double> u(0, 1);
  for (int t = 0; t < 2 * k; ++t) {
    const double cx = u(rng) * W, cy = u(rng) * H, r = (0.3 + 0.7 * u(rng)) * rmax * (t < k ? 1.0 : 0.45);
    const u8 val = t < k ? 1 : 0;
    const i32 x0 = std::max(0, i32(cx - r - 1)), x1 = std::min(W - 1, i32(cx + r + 1));
    const i32 y0 = std::max(0, i32(cy - r - 1)), y1 = std::min(H - 1, i32(cy + r + 1));
    for (i32 y = y0; y <= y1; ++y)
      for (i32 x = x0; x <= x1; ++x)
        if ((x + 0.5 - cx) * (x + 0.5 - cx) + (y + 0.5 - cy) * (y + 0.5 - cy) <= r * r) m.at(x, y) = val;
  }
  return m;
}

inline Mask checker(i32 H, i32 W) {  // every vertex inside is a saddle
  Mask m = blank(H, W, "checker");
  for (i32 y = 0; y < H; ++y)
    for (i32 x = 0; x < W; ++x) m.at(x, y) = (x + y) & 1;
  return m;
}

inline Mask concentric(i32 H, i32 W) {  // nested square rings: holes inside parts inside holes ...
  Mask m = blank(H, W, "concentric");
  for (i32 y = 0; y < H; ++y)
    for (i32 x = 0; x < W; ++x) {
      const i32 d = std::min(std::min(x, y), std::min(W - 1 - x, H - 1 - y));
      m.at(x, y) = (d / 2) % 2 == 0;
    }
  return m;
}

inline Mask comb(i32 H, i32 W) {  // one ring of length ~ H * W: vertical teeth joined at the bottom
  Mask m = blank(H, W, "comb");
  for (i32 y = 0; y < H; ++y)
    for (i32 x = 0; x < W; ++x) m.at(x, y) = (y == H - 1) || (x % 2 == 0);
  return m;
}

inline Mask spiral(i32 H, i32 W) {  // one long thin spiral (long ring, deep union-find chains)
  Mask m = blank(H, W, "spiral");
  i32 x0 = 0, y0 = 0, x1 = W - 1, y1 = H - 1;
  while (x0 <= x1 && y0 <= y1) {
    for (i32 x = x0; x <= x1; ++x) m.at(x, y0) = 1;
    for (i32 y = y0; y <= y1; ++y) m.at(x1, y) = 1;
    if (y1 > y0) for (i32 x = x0; x <= x1; ++x) m.at(x, y1) = 1;
    if (x0 + 2 <= x1) for (i32 y = y0 + 2; y <= y1; ++y) m.at(x0, y) = 1;
    if (x0 + 2 <= x1 - 2 && y0 + 2 <= y1) m.at(x0 + 1, y0 + 2) = 1;
    x0 += 2, y0 += 2, x1 -= 2, y1 -= 2;
  }
  return m;
}

inline Mask diagonals(i32 H, i32 W) {  // diagonal staircases: 8-connected parts, many saddles, 4-split background
  Mask m = blank(H, W, "diagonals");
  for (i32 y = 0; y < H; ++y)
    for (i32 x = 0; x < W; ++x) m.at(x, y) = ((x + y) % 5 == 0) || ((x - y + 1000) % 7 == 0);
  return m;
}

inline Mask stripes(i32 H, i32 W, bool vertical) {
  Mask m = blank(H, W, vertical ? "vstripes" : "hstripes");
  for (i32 y = 0; y < H; ++y)
    for (i32 x = 0; x < W; ++x) m.at(x, y) = vertical ? (x % 2 == 0) : (y % 2 == 0);
  return m;
}

inline Mask full(i32 H, i32 W, u8 v) {
  Mask m = blank(H, W, v ? "full" : "empty");
  std::fill(m.px.begin(), m.px.end(), v);
  return m;
}

inline Mask frame_with_holes(i32 H, i32 W) {  // border-touching foreground with many 1-pixel holes
  Mask m = full(H, W, 1);
  m.name = "frame_holes";
  for (i32 y = 1; y + 1 < H; y += 2)
    for (i32 x = 1; x + 1 < W; x += 2) m.at(x, y) = 0;
  return m;
}

struct Batch {  // gpucheck / k1poly packed layout
  std::vector<u8> bits;
  std::vector<i64> start{0};
  std::vector<i32> H, W;
  void add(const Mask& m) {
    bits.insert(bits.end(), m.px.begin(), m.px.end());
    start.push_back(i64(bits.size()));
    H.push_back(m.H), W.push_back(m.W);
  }
  k1poly::Masks view() const { return {i64(H.size()), bits.data(), start.data(), H.data(), W.data()}; }
};

}  // namespace k1gen
